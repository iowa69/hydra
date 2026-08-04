"""Translated-search gene calling and point-mutation detection.

Assembly contigs are searched with ``blastx`` against the AMRFinderPlus protein
reference. The same alignments serve two purposes: they yield acquired-gene
calls (with the AMR / STRESS / VIRULENCE element types that ``--plus`` exposes),
and they carry the gapped alignment strings the protein point-mutation caller
needs. Organism-specific DNA mutations (23S rRNA, promoters, *pbp4*) are called
from a separate ``blastn`` against the organism's reference loci.
"""

from __future__ import annotations

from pathlib import Path

from ..config import Config
from ..db.manager import DatabaseStore
from ..records import Hit
from ..seqio import read_fasta
from ..utils import LOG, natural_key
from .blast import (Hsp, blast, check_db_exists, deduplicate, interval_length,
                    merge_hsps, merge_intervals)
from .mutations import MutationCatalog, evaluate, multi_residue, walk_alignment
from .nucl import QueryBatch


class ProteinScreener:
    """AMRFinderPlus-equivalent translated search plus point-mutation calling."""

    def __init__(self, store: DatabaseStore, config: Config):
        self.store = store
        self.config = config
        self._catalogs: dict[str | None, MutationCatalog] = {}

    def catalog(self, organism: str | None) -> MutationCatalog:
        """Mutation catalogue for one organism, built once and reused."""
        if organism not in self._catalogs:
            self._catalogs[organism] = MutationCatalog(self.store.root, organism)
        return self._catalogs[organism]

    # ------------------------------------------------------------------ genes
    def screen(self, batch: QueryBatch, workdir: Path,
               organism_by_sample: dict[str, str | None] | None = None,
               threads: int | None = None) -> dict[str, list[Hit]]:
        """Call acquired genes and protein point mutations for every sample.

        Each sample is evaluated against the mutation catalogue of *its own*
        organism, so a mixed-species batch is still typed correctly from one
        translated search.
        """
        handle = self.store.handle("amrfinderplus")
        check_db_exists(handle.fasta, "prot")
        threads = threads or self.config.threads
        thresholds = self.config.thresholds
        organism_by_sample = organism_by_sample or {}
        organisms = sorted({o for o in organism_by_sample.values() if o})
        LOG.info("translated search of %d contigs against AMRProt (%d proteins)%s",
                 batch.n_contigs, handle.n_sequences,
                 f"; organisms: {', '.join(organisms)}" if organisms else "")

        out_tab = workdir / "amrfinderplus.blastx.tsv"
        hsps = blast(
            "blastx", batch.path, handle.fasta, out_tab,
            # A contig carrying several resistance genes matches hundreds of
            # alleles of each; a low cap makes BLAST discard whole genes at the
            # preliminary stage, so this stays as generous as the nucleotide path.
            threads=threads, evalue=1e-10, max_target_seqs=10000,
            with_alignment=True,
            # blastx-fast uses a longer seed word: ~4x quicker with no measurable
            # loss on references this similar. Composition-based statistics are
            # off because every alignment here is to a curated full-length protein.
            task="blastx-fast",
            extra=["-query_gencode", "11", "-comp_based_stats", "0"],
        )
        results: dict[str, list[Hit]] = {}
        gene_hits = self._call_genes(batch, hsps, handle, organism_by_sample, thresholds)
        for sample, hits in gene_hits.items():
            results.setdefault(sample, []).extend(hits)
        mutation_hits = self._call_protein_mutations(batch, hsps, handle, organism_by_sample,
                                                     thresholds)
        for sample, hits in mutation_hits.items():
            results.setdefault(sample, []).extend(hits)
        return results

    def _call_genes(self, batch: QueryBatch, hsps: list[Hsp], handle,
                    organism_by_sample: dict[str, str | None],
                    thresholds) -> dict[str, list[Hit]]:
        merged = merge_hsps(hsps, query_scale=3.0)
        per_sample: dict[str, list[tuple]] = {}
        for hit in merged:
            meta = handle.meta_for(hit.sseqid)
            if not meta:
                continue
            piece = batch.id_map.get(hit.qseqid)
            if piece is None:
                continue
            catalog = self.catalog(organism_by_sample.get(piece.sample))
            accession = meta.get("accession", "")
            if accession in catalog.suppress:
                continue
            # Reference proteins that exist only to anchor point mutations (gyrA,
            # rpoB, the PBPs) are core housekeeping genes present in every genome;
            # they are not acquired resistance and are left to the mutation caller.
            if meta.get("element_subtype") == "POINT" or meta.get("part_of_gene") == "mutation":
                continue
            element_type = meta.get("element_type", "AMR")
            if element_type != "AMR" and not self.config.plus:
                continue
            if hit.coverage_pct < thresholds.protein_partial_coverage:
                continue
            if hit.identity_pct < thresholds.protein_min_identity:
                continue
            per_sample.setdefault(piece.sample, []).append((hit, piece, meta))

        out: dict[str, list[Hit]] = {}
        for sample, entries in per_sample.items():
            entries = deduplicate(
                entries,
                key_span=lambda e: (e[0].query_start + e[1].offset, e[0].query_end + e[1].offset),
                key_seq=lambda e: e[1].contig,
                key_score=lambda e: (round(e[0].match_score, 4), round(e[0].identity_pct, 4),
                                     e[0].bitscore),
                key_tiebreak=lambda e: natural_key(e[2].get("gene", "")),
                overlap_fraction=1.0 if self.config.report_overlaps
                else self.config.overlap_fraction,
            )
            hits: list[Hit] = []
            for merged_hit, piece, meta in entries:
                exact = merged_hit.identity_pct >= 99.999 and merged_hit.coverage_pct >= 99.999
                partial = merged_hit.coverage_pct < thresholds.protein_min_coverage
                if partial and merged_hit.coverage_pct < thresholds.protein_partial_coverage:
                    continue
                method = "EXACTX" if exact else ("PARTIALX" if partial else "BLASTX")
                hits.append(Hit(
                    sample=sample, database="amrfinderplus",
                    gene=meta.get("gene", merged_hit.sseqid),
                    accession=meta.get("accession", ""),
                    product=meta.get("product", ""),
                    element_type=meta.get("element_type", "AMR"),
                    element_subtype=meta.get("element_subtype", "AMR"),
                    drug_class=meta.get("class", ""), subclass=meta.get("subclass", ""),
                    sequence=piece.contig, start=merged_hit.query_start + piece.offset,
                    end=merged_hit.query_end + piece.offset,
                    strand=merged_hit.strand, coverage=merged_hit.coverage_string,
                    coverage_pct=merged_hit.coverage_pct, identity_pct=merged_hit.identity_pct,
                    gaps=merged_hit.gaps, bitscore=merged_hit.bitscore,
                    method=method, resolution="PARTIAL" if partial else "COMPLETE",
                ))
            # Mosaic PBP-style calls: identity *below* the catalogued cutoff is the signal.
            hits.extend(self._call_susceptible(
                sample, entries, self.catalog(organism_by_sample.get(sample))))
            hits.sort(key=lambda h: (h.sequence, h.start))
            out[sample] = hits
        return out

    @staticmethod
    def _call_susceptible(sample: str, entries, catalog: MutationCatalog) -> list[Hit]:
        """Report proteins whose divergence from a susceptible reference confers resistance."""
        best: dict[str, tuple] = {}
        for merged_hit, piece, meta in entries:
            accession = meta.get("accession", "")
            if accession not in catalog.susceptible:
                continue
            current = best.get(accession)
            if current is None or merged_hit.bitscore > current[0].bitscore:
                best[accession] = (merged_hit, piece, meta)
        hits: list[Hit] = []
        for accession, (merged_hit, piece, meta) in best.items():
            rule = catalog.susceptible[accession]
            if merged_hit.identity_pct >= rule["cutoff"]:
                continue
            hits.append(Hit(
                sample=sample, database="amrfinderplus", gene=rule["gene"], accession=accession,
                product=rule["name"].replace("_", " "), element_type="AMR",
                element_subtype="POINT", drug_class=rule["class"], subclass=rule["subclass"],
                sequence=piece.contig, start=merged_hit.query_start + piece.offset,
                end=merged_hit.query_end + piece.offset,
                strand=merged_hit.strand, coverage=merged_hit.coverage_string,
                coverage_pct=merged_hit.coverage_pct, identity_pct=merged_hit.identity_pct,
                gaps=merged_hit.gaps, bitscore=merged_hit.bitscore,
                method="SUSCEPTIBLEX", resolution="POINT",
                note=f"divergent from susceptible reference (<{rule['cutoff']:.0f}% identity)",
            ))
        return hits

    # -------------------------------------------------------------- mutations
    def _call_protein_mutations(self, batch: QueryBatch, hsps: list[Hsp], handle,
                                organism_by_sample: dict[str, str | None],
                                thresholds) -> dict[str, list[Hit]]:
        # Index the best HSP set per (sample, contig, reference protein).
        by_target: dict[tuple[str, str, str], list[Hsp]] = {}
        for hsp in hsps:
            meta = handle.meta_for(hsp.sseqid)
            accession = meta.get("accession", "")
            piece = batch.id_map.get(hsp.qseqid)
            if piece is None:
                continue
            catalog = self.catalog(organism_by_sample.get(piece.sample))
            if accession not in catalog.protein:
                continue
            by_target.setdefault((piece.sample, piece.contig, piece.offset,
                                  hsp.sseqid), []).append(hsp)
        if not by_target:
            return {}

        out: dict[str, list[Hit]] = {}
        # Keep only the best contig per (sample, protein) so a duplicated gene is
        # not reported twice with conflicting alleles.
        best_locus: dict[tuple[str, str], tuple[float, tuple]] = {}
        for key, group in by_target.items():
            sample, _contig, _offset, sseqid = key
            score = sum(h.bitscore for h in group)
            marker = (sample, sseqid)
            if marker not in best_locus or score > best_locus[marker][0]:
                best_locus[marker] = (score, key)

        for _marker, (_score, key) in best_locus.items():
            sample, contig, offset, sseqid = key
            group = by_target[key]
            meta = handle.meta_for(sseqid)
            accession = meta["accession"]
            entries = self.catalog(organism_by_sample.get(sample)).protein.get(accession, [])
            if not entries:
                continue
            reference_length = int(meta.get("length", 0) or 0)
            covered = interval_length(merge_intervals(
                (min(h.sstart, h.send), max(h.sstart, h.send)) for h in group))
            if reference_length and 100.0 * covered / reference_length < thresholds.mutation_min_coverage:
                continue
            calls = self._match_entries(group, entries, reference_length, nucleotide=False)
            for entry, obs, hsp in calls:
                out.setdefault(sample, []).append(Hit(
                    sample=sample, database="amrfinderplus", gene=entry.gene,
                    accession=accession, product=entry.name.replace("_", " "),
                    element_type="AMR", element_subtype="POINT",
                    drug_class=entry.drug_class, subclass=entry.subclass,
                    sequence=contig, start=min(hsp.qstart, hsp.qend) + offset,
                    end=max(hsp.qstart, hsp.qend) + offset,
                    strand=hsp.strand, coverage=f"{entry.position}/{reference_length or '?'}",
                    coverage_pct=100.0 * covered / reference_length if reference_length else 0.0,
                    identity_pct=hsp.pident, gaps=hsp.gaps, bitscore=hsp.bitscore,
                    method="POINTX", resolution="POINT",
                    note=f"{entry.symbol}" + (f" ({obs.observed})" if obs.observed else ""),
                ))
        for sample in out:
            out[sample].sort(key=lambda h: (h.gene, h.note))
        return out

    @staticmethod
    def _match_entries(group: list[Hsp], entries, reference_length: int, *, nucleotide: bool):
        """Test every catalogue entry against the best-covering HSP for its position."""
        calls = []
        for entry in entries:
            position = entry.position
            if position < 0 and reference_length:
                position = reference_length + position + 1
            span = max(1, len(entry.ref)) if entry.is_deletion else 1
            best = None
            for hsp in group:
                lo, hi = min(hsp.sstart, hsp.send), max(hsp.sstart, hsp.send)
                if lo <= position <= hi and (best is None or hsp.bitscore > best.bitscore):
                    best = hsp
            if best is None or not best.qseq:
                continue
            if span > 1:
                obs = multi_residue(best.qseq, best.sseq, best.sstart, best.send,
                                    position, span, nucleotide=nucleotide)
            else:
                seen = walk_alignment(best.qseq, best.sseq, best.sstart, best.send,
                                      {position: None}, nucleotide=nucleotide)
                obs = seen[position]
            is_mutant, note = evaluate(entry, obs)
            if is_mutant:
                calls.append((entry, obs, best))
            elif note.startswith("reference mismatch"):
                LOG.debug("skipped %s: %s", entry.symbol, note)
        return calls

    def call_dna_mutations(self, batch: QueryBatch, workdir: Path, organism: str,
                           threads: int | None = None) -> dict[str, list[Hit]]:
        """Call organism-specific DNA point mutations (23S rRNA, promoters, ...)."""
        catalog = MutationCatalog(self.store.root, organism)
        if not catalog.dna or catalog.dna_reference is None:
            return {}
        threads = threads or self.config.threads
        thresholds = self.config.thresholds
        reference_lengths = {header.split()[0]: len(seq)
                             for header, seq in read_fasta(catalog.dna_reference)}
        LOG.info("screening %s DNA mutation loci (%d catalogued positions)",
                 organism, catalog.n_dna)
        out_tab = workdir / f"mutations.{organism}.blastn.tsv"
        hsps = blast(
            "blastn", batch.path, catalog.dna_reference, out_tab,
            threads=threads, evalue=1e-20, task="blastn",
            perc_identity=thresholds.mutation_min_identity - 10.0,
            max_target_seqs=1000, with_alignment=True,
        )
        by_target: dict[tuple[str, str, str], list[Hsp]] = {}
        for hsp in hsps:
            if hsp.sseqid not in catalog.dna:
                continue
            piece = batch.id_map.get(hsp.qseqid)
            if piece is None:
                continue
            by_target.setdefault((piece.sample, piece.contig, piece.offset,
                                  hsp.sseqid), []).append(hsp)

        # Multi-copy loci (rRNA operons) appear on several contigs; keep the best
        # copy per sample so the assembly consensus is reported once.
        best_locus: dict[tuple[str, str], tuple[float, tuple]] = {}
        for key, group in by_target.items():
            sample, _contig, _offset, sseqid = key
            score = sum(h.bitscore for h in group)
            marker = (sample, sseqid)
            if marker not in best_locus or score > best_locus[marker][0]:
                best_locus[marker] = (score, key)

        out: dict[str, list[Hit]] = {}
        for _marker, (_score, key) in best_locus.items():
            sample, contig, offset, sseqid = key
            group = by_target[key]
            # Merge, never sum: a multi-copy locus such as an rRNA operon puts
            # several alignments on the same reference, and summing their
            # lengths would report coverage above 100% and pass the gate with
            # fragments that never touch the catalogued position.
            reference_length = reference_lengths.get(sseqid, 0)
            covered = interval_length(merge_intervals(
                (min(h.sstart, h.send), max(h.sstart, h.send)) for h in group))
            if reference_length and 100.0 * covered / reference_length < thresholds.mutation_min_coverage:
                continue
            entries = catalog.dna[sseqid]
            calls = self._match_entries(group, entries, reference_length, nucleotide=True)
            for entry, obs, hsp in calls:
                out.setdefault(sample, []).append(Hit(
                    sample=sample, database="amrfinderplus", gene=entry.gene,
                    accession=sseqid.split("@")[0], product=entry.name.replace("_", " "),
                    element_type="AMR", element_subtype="POINT",
                    drug_class=entry.drug_class, subclass=entry.subclass,
                    sequence=contig, start=min(hsp.qstart, hsp.qend) + offset,
                    end=max(hsp.qstart, hsp.qend) + offset,
                    strand=hsp.strand,
                    coverage=f"{entry.position}/{reference_length or '?'}",
                    coverage_pct=100.0 * covered / reference_length if reference_length else 0.0,
                    identity_pct=hsp.pident, gaps=hsp.gaps, bitscore=hsp.bitscore,
                    method="POINTN", resolution="POINT",
                    note=f"{entry.symbol}" + (f" ({obs.observed})" if obs.observed else ""),
                ))
        for sample in out:
            out[sample].sort(key=lambda h: (h.gene, h.note))
        return out
