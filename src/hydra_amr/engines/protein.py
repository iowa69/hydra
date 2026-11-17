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
from ..seqio import read_fasta, read_fasta_dict, read_fasta_subset
from ..utils import LOG, natural_key, revcomp, translate
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
    # oqxR represses the OqxAB efflux pump, so losing it raises the MIC of
    # everything the pump exports at once.
    "oqxR": {
        "organisms": _KPSC,
        "class": "NITROFURAN/PHENICOL/QUINOLONE/TETRACYCLINE",
        "subclass": "NITROFURANTOIN/PHENICOL/QUINOLONE/TIGECYCLINE",
        "intact_coverage": 90.0, "min_partial_coverage": 20.0, "min_identity": 90.0,
    },
}



def query_frame(hsp: Hsp) -> int:
    """The reading frame a translated alignment sits in, 1..3 or -1..-3."""
    if hsp.qend >= hsp.qstart:
        return ((hsp.qstart - 1) % 3) + 1
    return -(((hsp.qlen - hsp.qstart) % 3) + 1)


def frameshifted(hsps: list[Hsp], reference_length: int) -> bool:
    """True when the alignments at one locus change reading frame.

    A frameshift is invisible to a coverage test. An indel of one or two bases
    near the start of a gene destroys everything downstream, but the downstream
    part still aligns to the reference in a different frame, so the pieces tile
    the reference and the gene reads as complete: ompK35 broken at residue 29
    still covers 92% of its reference in a single stretch and passes any
    "intact" threshold.

    Two alignments of one reference, adjacent at one locus but in different
    reading frames, are the frameshift itself. Requiring them to be adjacent is
    what separates this from a second, intact copy of the gene elsewhere, which
    may sit in any frame for reasons that have nothing to do with disruption.
    """
    ordered = sorted(hsps, key=lambda h: min(h.qstart, h.qend))
    for left, right in zip(ordered, ordered[1:]):
        gap = min(right.qstart, right.qend) - max(left.qstart, left.qend)
        # Within one gene's span: the pieces of a broken gene stay together.
        if not -60 <= gap <= reference_length * 3:
            continue
        if query_frame(left) != query_frame(right):
            return True
    return False


def premature_stop(hsp: Hsp, reference_length: int, cutoff: float = 0.95) -> int:
    """Reference position of the first in-frame stop, or 0 if there is none.

    ``blastx`` aligns straight through a stop codon and writes it as ``*`` in the
    translated query, so a nonsense mutation leaves full-length coverage at high
    identity and is invisible to any coverage or identity test: ompK35 stopped at
    residue 138 of 359 still aligns over the whole reference at 99.7%.

    Only a stop before *cutoff* of the reference counts. The last codon of a
    protein is a stop in every genome and means nothing.
    """
    if not hsp.qseq or not hsp.sseq or not reference_length:
        return 0
    limit = reference_length * cutoff
    step = 1 if hsp.send >= hsp.sstart else -1
    position = hsp.sstart
    for query_residue, subject_residue in zip(hsp.qseq, hsp.sseq):
        if subject_residue == "-":
            continue
        if query_residue == "*" and position <= limit:
            return position
        position += step
    return 0


#: Shortest unaligned end worth reading the DNA for.
EDGE_MIN_MISSING = 5


def edge_frameshift(reference: str, query: str, hsp: Hsp, probe: int = 6,
                    min_missing: int = EDGE_MIN_MISSING, max_missing: int = 40) -> str:
    """Which end of a gene falls out of frame -- ``"N"``, ``"C"`` or ``""``.

    :func:`frameshifted` needs two alignments of one reference in different
    frames, and blastx supplies them only when both pieces score. A gene
    frameshifted within its first few codons has no second piece: the fragment
    is 10 or 15 residues long and never reaches the reporting threshold, so the
    lesion comes back as a *single* clean alignment that simply starts late.
    ``ompK35`` broken at codon 12 aligns from residue 15 to 359 at 100%
    identity, which no coverage or identity test can tell from an intact gene
    whose N-terminus went unaligned.

    What separates them is the DNA that blastx could not score. The two ends
    fail differently, so they are tested differently:

    ``N``  the reference's first residues are still there, immediately upstream,
           but in one of the other two reading frames. The start codon exists
           and the body of the gene is not in frame with it.
    ``C``  translating in frame past the end of the alignment runs into a stop
           before the reference ends, reading residues that are not the
           reference's. A gene whose last residues are merely absent, with no
           novel sequence before the stop, is a short allele and is not
           reported.

    Both are bounded by *max_missing*, and the bound is the whole point of the
    test rather than a safety margin. What justifies reading the DNA at all is
    that the missing piece is too short for blastx to score; a piece of 40
    residues or more would have been reported as its own alignment, and
    :func:`frameshifted` would already have two frames to compare. Past that
    bound the unaligned part is not an unscored fragment but a gene in pieces,
    which is the truncation rule's business: without the bound a locus aligning
    over 46 of 359 residues was called frameshifted on the strength of a
    six-residue probe found somewhere in 900 bases of flank.
    """
    if not reference or not query or not hsp.slen:
        return ""
    # A translated search walks the protein forwards; anything else is not an
    # alignment these tests can read.
    if hsp.send < hsp.sstart:
        return ""
    plus = hsp.qend >= hsp.qstart

    missing = reference[:hsp.sstart - 1]
    if min_missing <= len(missing) <= max_missing:
        width = 3 * (len(missing) + 6)
        flank = (query[max(0, hsp.qstart - 1 - width):hsp.qstart - 1] if plus
                 else revcomp(query[hsp.qstart:hsp.qstart + width]))
        # Trim from the far end so the last codon of the flank abuts the
        # alignment: frame 0 is then the alignment's own frame.
        flank = flank[len(flank) % 3:]
        if len(flank) >= 3 * min_missing:
            head = missing[:probe]
            if head not in translate(flank, False) and any(
                    head in translate(flank[shift:], False) for shift in (1, 2)):
                return "N"

    missing = reference[hsp.send:]
    if min_missing <= len(missing) <= max_missing:
        width = 3 * (len(missing) + 30)
        tail = (query[hsp.qend:hsp.qend + width] if plus
                else revcomp(query[max(0, hsp.qend - 1 - width):hsp.qend - 1]))
        peptide = translate(tail, False)
        stop = peptide.find("*")
        if 0 < stop < len(missing) and peptide[:stop] != missing[:stop]:
            return "C"
    return ""


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
            threads=threads, evalue=1e-3, max_target_seqs=10000, with_alignment=True,
            extra=["-query_gencode", "11", "-comp_based_stats", "0"],
        )
        best: dict[tuple, tuple] = {}
        #: Every qualifying alignment per locus, kept so a change of reading
        #: frame can be seen. The best-coverage hit alone cannot show one.
        by_contig: dict[tuple, list[Hsp]] = {}
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
            by_contig.setdefault((piece.sample, gene, piece.contig), []).append(hsp)

        # The flank tests below read sequence, not alignments. The batch FASTA
        # holds every contig of every sample in the run, so only the contigs
        # carrying a locus with an unaligned end are loaded, and a run where no
        # locus has one never opens the file at all.
        references = read_fasta_dict(query)
        wanted_contigs = {hsp.qseqid for _cov, hsp, _piece, _meta in best.values()
                          if max(hsp.sstart - 1, hsp.slen - hsp.send) >= EDGE_MIN_MISSING}
        contigs = read_fasta_subset(batch.path, wanted_contigs) if wanted_contigs else {}

        out: dict[str, list[Hit]] = {}
        for (sample, gene), (coverage, hsp, piece, meta) in best.items():
            rule = LOSS_OF_FUNCTION[gene]
            broken_frame = any(
                frameshifted(group, hsp.slen)
                for (s, g, _contig), group in by_contig.items()
                if s == sample and g == gene)
            stop_at = premature_stop(hsp, hsp.slen)
            edge = edge_frameshift(references.get(hsp.sseqid, ""),
                                   contigs.get(hsp.qseqid, ""), hsp)
            # Organism-scoped, like the mutation catalogues. mgrB disruption means
            # colistin resistance in Klebsiella; in an organism where that has not
            # been established -- or where the species call names none -- it means
            # nothing and is not reported.
            organism = organism_by_sample.get(sample)
            truncated = loss_of_function_call(gene, organism, coverage, hsp.pident)
            # A frameshift is disruption whatever the coverage says, but the
            # organism gate and the identity floor still apply: a paralogue
            # aligning in two frames is not this gene breaking.
            credible = (organism in rule["organisms"]
                        and hsp.pident >= rule["min_identity"])
            broken = truncated or ((broken_frame or stop_at or edge) and credible)
            if not broken:
                continue
            if truncated:
                reason, detail = "truncated", (
                    f"longest aligned stretch covers {coverage:.0f}% of the reference")
            elif stop_at:
                reason, detail = "nonsense", (
                    f"stop codon at reference residue {stop_at} of {hsp.slen}")
            elif broken_frame:
                reason, detail = "frameshifted", (
                    "alignments at this locus change reading frame")
            elif edge == "N":
                reason, detail = "frameshifted", (
                    f"the reference's first {hsp.sstart - 1} residues lie upstream "
                    f"in a different reading frame")
            else:
                reason, detail = "frameshifted", (
                    f"translation past reference residue {hsp.send} of {hsp.slen} "
                    f"reads novel sequence and stops early")
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
                note=(f"{gene} {reason}: {detail}; loss of function confers "
                      f"{rule['class'].lower()} resistance"),
                # For engines/lesions.py, which names the lesion from the sequence.
                evidence={"kind": "disruption", "hsp": hsp, "offset": piece.offset,
                          "contig": piece.contig, "reason": reason,
                          "reference": references.get(hsp.sseqid, ""),
                          "pieces": [(h, batch.id_map[h.qseqid])
                                     for (s, g, _contig), group in by_contig.items()
                                     if s == sample and g == gene for h in group]},
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
            # One lesion, one row. A catalogued nonsense mutation is found twice
            # over -- once by the mutation caller as ompK36_W125Ter and once here
            # as a stop at residue 125 - and reporting both says the gene is
            # broken twice.
            catalogued = {(hit.gene, hit.note.split()[0]) for hit in mutation_hits.get(sample, [])}
            stops = {gene for gene, symbol in catalogued if symbol.endswith(("Ter", "*"))}
            results.setdefault(sample, []).extend(
                hit for hit in hits
                if not ("nonsense" in hit.note and hit.gene in stops))
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
            calls = self._match_entries(group, entries, reference_length, nucleotide=False,
                                        min_identity=thresholds.mutation_min_identity)
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
                    evidence={"kind": "protein", "hsp": hsp, "offset": offset,
                              "entry": entry, "position": entry.position},
                ))
        for sample in out:
            out[sample].sort(key=lambda h: (h.gene, h.note))
        return out

    @staticmethod
    def _copies_spanning(group: list[Hsp], position: int,
                         min_identity: float = 0.0) -> list[Hsp]:
        """One best alignment per distinct copy of the locus covering *position*.

        A genome often carries the same reference locus more than once, and the
        copies need not agree: a chromosome can hold one wild-type blaSHV
        promoter and two mutant ones. Testing only the highest-scoring alignment
        reads whichever copy aligns best and reports nothing when that copy
        happens to be the wild type, so every copy is tested and the mutation is
        called if any of them carries it.

        Since every copy gets a vote, each has to be a credible copy. An
        alignment below *min_identity* is a paralogue or a chance match rather
        than the gene, and a residue read off one is not evidence: a topoisomerase
        aligning to parC at 36% identity will happily show a serine at position
        80 and turn a mutation call into a coin toss.

        Copies are separated on the query, and only the best alignment of each is
        kept, so a fragmentary HSP cannot outvote the real one.
        """
        spanning = [hsp for hsp in group
                    if min(hsp.sstart, hsp.send) <= position <= max(hsp.sstart, hsp.send)
                    and hsp.qseq and hsp.pident >= min_identity]
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
    def _match_entries(group: list[Hsp], entries, reference_length: int, *,
                       nucleotide: bool, min_identity: float = 0.0):
        """Test every catalogue entry against each copy of the locus it sits in."""
        calls = []
        for entry in entries:
            position = entry.position
            if position < 0 and reference_length:
                position = reference_length + position + 1
            span = max(1, len(entry.ref)) if entry.is_deletion else 1
            for hsp in ProteinScreener._copies_spanning(group, position, min_identity):
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
        reference_seqs = {header.split()[0]: seq
                          for header, seq in read_fasta(catalog.dna_reference)}
        reference_lengths = {name: len(seq) for name, seq in reference_seqs.items()}
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
            calls = self._match_entries(group, entries, reference_length, nucleotide=True,
                                        min_identity=thresholds.mutation_min_identity)
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
                    evidence={"kind": "dna", "hsp": hsp, "offset": offset, "entry": entry,
                              "position": (entry.position if entry.position > 0
                                           else reference_length + entry.position + 1),
                              "reference": reference_seqs.get(sseqid, "")},
                ))
        for sample in out:
            out[sample].sort(key=lambda h: (h.gene, h.note))
        return out
