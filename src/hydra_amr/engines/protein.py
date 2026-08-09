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
from ..db.manager import DatabaseStore, makeblastdb
from ..records import Hit
from ..seqio import read_fasta
from ..utils import LOG, natural_key
from .blast import (Hsp, blast, check_db_exists, deduplicate, interval_length,
                    merge_hsps, merge_intervals, overlap_length)
from .mutations import MutationCatalog, evaluate, multi_residue, walk_alignment
from .nucl import QueryBatch


#: The K. pneumoniae species complex, where every rule below is established.
_KPSC = {"Klebsiella_pneumoniae", "Klebsiella_oxytoca",
         "Klebsiella_variicola", "Klebsiella_quasipneumoniae"}

#: Genes whose *loss* confers resistance, scoped to the organisms where that is
#: established, with the drug class each carries in the AMRFinderPlus catalogue.
#: Reporting one of these for an organism that does not have the gene, or whose
#: resistance to that drug works another way, would be inventing a result.
LOSS_OF_FUNCTION: dict[str, dict] = {
    # mgrB represses PhoPQ; losing it is the commonest route to colistin resistance
    # in K. pneumoniae. 47 residues, so a modest identity floor is enough -- it has
    # no close paralogue to be confused with.
    "mgrB": {
        "organisms": _KPSC, "class": "COLISTIN", "subclass": "COLISTIN",
        "intact_coverage": 90.0, "min_partial_coverage": 20.0, "min_identity": 80.0,
    },
    # The porins. Losing either raises the beta-lactam and carbapenem MIC, and with a
    # carbapenemase present that is the difference between susceptible and not.
    #
    # These need a high identity floor, unlike mgrB: ompK35 and ompK36 are homologous
    # (~58% identical) and each aligns to the other's reference over its full length.
    # At an 80% floor a paralogous full-length hit masks a real truncation, and both
    # genes read as intact in every genome.
    "ompK35": {
        "organisms": _KPSC, "class": "CARBAPENEM", "subclass": "CARBAPENEM",
        "intact_coverage": 90.0, "min_partial_coverage": 20.0, "min_identity": 90.0,
    },
    "ompK36": {
        "organisms": _KPSC, "class": "CARBAPENEM", "subclass": "CARBAPENEM",
        "intact_coverage": 90.0, "min_partial_coverage": 20.0, "min_identity": 90.0,
    },
    # cirA imports cefiderocol on its iron-siderophore route; without it the drug
    # does not get in.
    "cirA": {
        "organisms": _KPSC, "class": "CEFIDEROCOL", "subclass": "CEFIDEROCOL",
        "intact_coverage": 90.0, "min_partial_coverage": 20.0, "min_identity": 90.0,
    },
    # nfsB is the nitroreductase that activates nitrofurantoin.
    "nfsB": {
        "organisms": _KPSC, "class": "NITROFURAN", "subclass": "NITROFURANTOIN",
        "intact_coverage": 90.0, "min_partial_coverage": 20.0, "min_identity": 90.0,
    },
    # ramR represses ramA; losing it derepresses the AcrAB efflux pump.
    "ramR": {
        "organisms": _KPSC, "class": "TIGECYCLINE", "subclass": "TIGECYCLINE",
        "intact_coverage": 90.0, "min_partial_coverage": 20.0, "min_identity": 90.0,
    },
}



def loss_of_function_call(gene: str, organism: str | None,
                          coverage: float, identity: float) -> bool:
    """Whether a truncated reference should be reported as resistance.

    Kept separate from the search so the rule can be read and tested on its own.
    The organism gate is the important half: mgrB disruption confers colistin
    resistance in Klebsiella, and reporting it for an organism that has no mgrB --
    or whose colistin resistance works another way -- would be inventing a result
    from an alignment that happens to be short.
    """
    rule = LOSS_OF_FUNCTION.get(gene)
    if rule is None:
        return False
    if not organism or organism not in rule["organisms"]:
        return False
    if not rule["min_partial_coverage"] <= coverage < rule["intact_coverage"]:
        return False
    return identity >= rule["min_identity"]


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
    def _lof_reference_fasta(self, handle, workdir: Path) -> Path | None:
        """Write the loss-of-function reference proteins to their own small FASTA.

        AMRProt holds a plain reference for each of these genes and one entry per
        catalogued mutation in it, distinguished by a position in the accession
        (``WP_002911375.1`` against ``WP_002911375.1:47:mgrB_W47R``). Only the plain
        reference is wanted: the mutation entries are the same protein and would
        each report the same truncation again.
        """
        wanted = {seqid: meta.get("gene") for seqid, meta in handle.meta.items()
                  if meta.get("gene") in LOSS_OF_FUNCTION
                  and ":" not in (meta.get("accession") or "")}
        if not wanted:
            return None
        dest = workdir / "lossoffunction.faa"
        written = 0
        with open(handle.fasta) as src, open(dest, "w") as out:
            keep = False
            for line in src:
                if line.startswith(">"):
                    keep = line[1:].split()[0] in wanted
                    written += keep
                if keep:
                    # A trailing stop is part of the record, not of the protein, and
                    # counting it in the length would understate every coverage.
                    out.write(line.replace("*", "") if not line.startswith(">") else line)
        return dest if written else None

    def _loss_of_function(self, batch: QueryBatch, handle, workdir: Path,
                          organism_by_sample: dict[str, str | None],
                          threads: int) -> dict[str, list[Hit]]:
        """Report a regulator that is present but broken, where losing it confers
        resistance.

        Most resistance is a gene arriving. A few are a gene breaking, and those are
        invisible to a screen that reports only what it finds intact: in
        *K. pneumoniae* the commonest route to colistin resistance is *mgrB*, a
        47-residue repressor of PhoPQ, disrupted by an insertion sequence. There is no
        acquired gene to find and no catalogued substitution.

        This deliberately reads the *raw* HSPs rather than the merged hits. Merging is
        done on the subject axis, so a gene split in two by an inserted element adds
        back up to near-full coverage and looks intact -- which is right for a gene
        split by an assembly gap and exactly wrong here. The signal is the longest
        single ungapped stretch of the reference that still aligns.

        Only truncation counts, never absence: a gene missing from a draft assembly is
        far more often a contig boundary than a deletion. On the phenotyped isolates
        here, 22 of 100 colistin-resistant genomes carry a truncated *mgrB* and none
        of 200 susceptible ones do.
        """
        if not any(o in rule["organisms"]
                   for o in organism_by_sample.values() if o
                   for rule in LOSS_OF_FUNCTION.values()):
            return {}
        query = self._lof_reference_fasta(handle, workdir)
        if query is None:
            return {}
        # Its own search, and a permissive one. An e-value is scaled by database
        # size: a 30-residue stretch of a 47-residue protein scores 4e-12 against
        # this one reference and misses the main search's 1e-10 cutoff entirely
        # once it is competing with 9998 proteins. Searching the handful of
        # loss-of-function references on their own is what makes a partial hit
        # visible, and it costs one BLAST over a few hundred residues.
        makeblastdb(query, "prot")
        hsps = blast(
            "blastx", batch.path, query, workdir / "lossoffunction.blastx.tsv",
            threads=threads, evalue=1e-3, max_target_seqs=10000,
            extra=["-query_gencode", "11", "-comp_based_stats", "0"],
        )
        best: dict[tuple, tuple] = {}
        for hsp in hsps:
            meta = handle.meta_for(hsp.sseqid)
            if not meta:
                continue
            gene = meta.get("gene", hsp.sseqid)
            if gene not in LOSS_OF_FUNCTION or not hsp.slen:
                continue
            piece = batch.id_map.get(hsp.qseqid)
            if piece is None:
                continue
            # Reject paralogues before choosing the best alignment, not after.
            # ompK35 and ompK36 each align to the other's reference over its full
            # length at ~58% identity, so picking the longest hit first and checking
            # identity afterwards lets the paralogue win: the real gene's truncated
            # 75% hit is discarded in favour of a 105% one that is a different gene,
            # and the locus reads as intact in every genome.
            if hsp.pident < LOSS_OF_FUNCTION[gene]["min_identity"]:
                continue
            coverage = 100.0 * (abs(hsp.send - hsp.sstart) + 1) / hsp.slen
            key = (piece.sample, gene)
            if key not in best or coverage > best[key][0]:
                best[key] = (coverage, hsp, piece, meta)

        out: dict[str, list[Hit]] = {}
        for (sample, gene), (coverage, hsp, piece, meta) in best.items():
            rule = LOSS_OF_FUNCTION[gene]
            # Organism-scoped, like the mutation catalogues. mgrB disruption means
            # colistin resistance in Klebsiella; in an organism where that has not
            # been established -- or where the species call names none -- it means
            # nothing and is not reported.
            if not loss_of_function_call(gene, organism_by_sample.get(sample),
                                         coverage, hsp.pident):
                continue
            out.setdefault(sample, []).append(Hit(
                sample=sample, database="protein", gene=gene,
                accession=meta.get("accession", ""), product=meta.get("product", ""),
                element_type="AMR", element_subtype="DISRUPTION",
                drug_class=rule["class"], subclass=rule["subclass"],
                sequence=piece.contig, start=min(hsp.qstart, hsp.qend) + piece.offset,
                end=max(hsp.qstart, hsp.qend) + piece.offset, strand=hsp.strand,
                coverage=f"{min(hsp.sstart, hsp.send)}-{max(hsp.sstart, hsp.send)}/{hsp.slen}",
                coverage_pct=coverage, identity_pct=hsp.pident,
                gaps=hsp.gaps, bitscore=hsp.bitscore,
                method="DISRUPTX", resolution="PARTIAL",
                note=(f"{gene} truncated: longest aligned stretch covers "
                      f"{coverage:.0f}% of the reference; loss of function confers "
                      f"{rule['class'].lower()} resistance"),
            ))
        return out

    def screen(self, batch: QueryBatch, workdir: Path,
               organism_by_sample: dict[str, str | None] | None = None,
               threads: int | None = None) -> dict[str, list[Hit]]:
        """Call acquired genes and protein point mutations for every sample.

        Each sample is evaluated against the mutation catalogue of *its own*
        organism, so a mixed-species batch is still typed correctly from one
        translated search.
        """
        handle = self.store.handle("protein")
        check_db_exists(handle.fasta, "prot")
        threads = threads or self.config.threads
        thresholds = self.config.thresholds
        organism_by_sample = organism_by_sample or {}
        organisms = sorted({o for o in organism_by_sample.values() if o})
        LOG.info("translated search of %d contigs against AMRProt (%d proteins)%s",
                 batch.n_contigs, handle.n_sequences,
                 f"; organisms: {', '.join(organisms)}" if organisms else "")

        out_tab = workdir / "protein.blastx.tsv"
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
        for sample, hits in self._loss_of_function(batch, handle, workdir,
                                                   organism_by_sample, threads).items():
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
                    sample=sample, database="protein",
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
                sample=sample, database="protein", gene=rule["gene"], accession=accession,
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
        # Index HSPs per (sample, contig, reference protein). AMRProt carries a
        # plain copy of each catalogued protein and one copy per mutation in it,
        # the latter marked by a position in the accession
        # (``WP_002911375.1`` against ``WP_002911375.1:47:mgrB_W47R``). Both are
        # indexed under the plain accession, because only the mutant copy aligns
        # across a variant at the very end of a short protein: mgrB is 47
        # residues and the wild-type copy stops at 43, leaving the catalogued
        # position uncovered and the mutation uncalled.
        by_target: dict[tuple[str, str, str], list[Hsp]] = {}
        for hsp in hsps:
            meta = handle.meta_for(hsp.sseqid)
            accession = (meta.get("accession", "") or "").split(":")[0]
            piece = batch.id_map.get(hsp.qseqid)
            if piece is None:
                continue
            catalog = self.catalog(organism_by_sample.get(piece.sample))
            if accession not in catalog.protein:
                continue
            by_target.setdefault((piece.sample, piece.contig, piece.offset,
                                  accession), []).append(hsp)
        if not by_target:
            return {}

        # Every copy of the gene is examined. Keeping only the best-aligning
        # contig reports whichever copy scores highest and misses a mutation
        # carried by another copy; one row per mutation is enforced on the call
        # below instead, so a duplicated gene still yields a single finding.
        out: dict[str, list[Hit]] = {}
        best_call: dict[tuple[str, str], float] = {}
        for key in sorted(by_target, key=lambda k: -sum(h.bitscore for h in by_target[k])):
            sample, contig, offset, accession = key
            group = by_target[key]
            entries = self.catalog(organism_by_sample.get(sample)).protein.get(accession, [])
            if not entries:
                continue
            reference_length = max(
                (int(handle.meta_for(h.sseqid).get("length", 0) or 0) for h in group),
                default=0)
            covered = interval_length(merge_intervals(
                (min(h.sstart, h.send), max(h.sstart, h.send)) for h in group))
            if reference_length and 100.0 * covered / reference_length < thresholds.mutation_min_coverage:
                continue
            calls = self._match_entries(group, entries, reference_length, nucleotide=False)
            for entry, obs, hsp in calls:
                # One row per mutation per sample, from the strongest copy that
                # carries it; groups are visited best-scoring first.
                marker = (sample, entry.symbol)
                if marker in best_call:
                    continue
                best_call[marker] = hsp.bitscore
                out.setdefault(sample, []).append(Hit(
                    sample=sample, database="protein", gene=entry.gene,
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
    def _copies_spanning(group: list[Hsp], position: int) -> list[Hsp]:
        """One best alignment per distinct copy of the locus covering *position*.

        A genome often carries the same reference locus more than once, and the
        copies need not agree: a chromosome can hold one wild-type blaSHV
        promoter and two mutant ones. Testing only the highest-scoring alignment
        reads whichever copy aligns best and reports nothing when that copy
        happens to be the wild type, so every copy is tested and the mutation is
        called if any of them carries it.

        Copies are separated on the query, and only the best alignment of each is
        kept, so a fragmentary or spurious HSP cannot vote on the residue.
        """
        spanning = [hsp for hsp in group
                    if min(hsp.sstart, hsp.send) <= position <= max(hsp.sstart, hsp.send)
                    and hsp.qseq]
        copies: list[Hsp] = []
        for hsp in sorted(spanning, key=lambda h: h.bitscore, reverse=True):
            lo, hi = min(hsp.qstart, hsp.qend), max(hsp.qstart, hsp.qend)
            length = hi - lo + 1
            if any(overlap_length((lo, hi), (min(k.qstart, k.qend), max(k.qstart, k.qend)))
                   > 0.5 * length for k in copies):
                continue
            copies.append(hsp)
        return copies

    @staticmethod
    def _match_entries(group: list[Hsp], entries, reference_length: int, *, nucleotide: bool):
        """Test every catalogue entry against each copy of the locus it sits in."""
        calls = []
        for entry in entries:
            position = entry.position
            if position < 0 and reference_length:
                position = reference_length + position + 1
            span = max(1, len(entry.ref)) if entry.is_deletion else 1
            for hsp in ProteinScreener._copies_spanning(group, position):
                if span > 1:
                    obs = multi_residue(hsp.qseq, hsp.sseq, hsp.sstart, hsp.send,
                                        position, span, nucleotide=nucleotide)
                else:
                    seen = walk_alignment(hsp.qseq, hsp.sseq, hsp.sstart, hsp.send,
                                          {position: None}, nucleotide=nucleotide)
                    obs = seen[position]
                is_mutant, note = evaluate(entry, obs)
                if is_mutant:
                    # One row per mutation, from the strongest copy carrying it.
                    # A 23S rRNA operon puts the same change on seven copies and
                    # it is one finding, not seven.
                    calls.append((entry, obs, hsp))
                    break
                if note.startswith("reference mismatch"):
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

        # Every copy of the locus is examined, not just the best-aligning one.
        # Copies of a gene need not agree: a chromosome can carry one wild-type
        # blaSHV promoter and two mutant ones, and keeping only the strongest
        # alignment reports the wild type and misses the mutation entirely.
        # Multi-copy loci such as rRNA operons still yield one row per mutation,
        # because the de-duplication now happens on the call below rather than
        # on the locus here.
        out: dict[str, list[Hit]] = {}
        best_call: dict[tuple[str, str], float] = {}
        for key in sorted(by_target, key=lambda k: -sum(h.bitscore for h in by_target[k])):
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
                # One row per mutation per sample, from whichever copy supports
                # it most strongly; groups are visited best-scoring first.
                marker = (sample, entry.symbol)
                if marker in best_call:
                    continue
                best_call[marker] = hsp.bitscore
                out.setdefault(sample, []).append(Hit(
                    sample=sample, database="protein", gene=entry.gene,
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
