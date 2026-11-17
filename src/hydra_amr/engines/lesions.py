"""What exactly changed: every mutation call described from the isolate's own sequence.

The callers in :mod:`.protein` decide *whether* a catalogued mutation or a gene
disruption is present. This module never revisits that decision. For each call
already made it reads the assembly at the place the alignment points to and
states the change the way AMRFinderPlus and Kleborate name it:

* the codon (or base) the isolate carries, read from the contig at coordinates
  derived by walking the gapped alignment -- never by arithmetic that ignores
  indels -- with the strand handled and the codon re-translated as a check;
* the reference codon or base, taken from a reference nucleotide sequence when
  the database has one (the organism-specific DNA loci: rRNA, promoters), and
  otherwise *inferred* by parsimony and labelled as inferred: AMRFinderPlus
  ships proteins, not coding sequences, for the point-mutation genes;
* for a disrupted gene, the causal lesion: premature stop, frameshift (with the
  net indel, where the new frame stops, and the HGVS-style ``fsTer`` symbol),
  interruption by inserted sequence, loss of the start codon -- or, when the
  gene simply runs off the end of a contig, that this is an assembly boundary
  and not evidence of loss of function;
* the resistance the call is associated with, in words.

Anything that cannot be resolved to a codon is said so in ``unresolved`` rather
than guessed. Nomenclature follows AMRFinderPlus (``gyrA_S83I``,
``mgrB_Q30Ter``, ``cirA_Y253CfsTer5``, ``ompK36_D135DGD``) with HGVS
(``p.Ser83Ile``, ``c.248C>T``) alongside, and Kleborate v3's calls
(``GyrA:p.Ser83Ile``, ``OmpK35:p.Trp133fs``) for the K. pneumoniae loci
Kleborate types.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

from ..records import CONTIG_END_RESOLUTION, CONTIG_END_SUBTYPE, Hit, MutationDetail
from ..seqio import read_fasta_subset
from ..utils import LOG, revcomp, translate

_BASES = "TCAG"
#: The bacterial code (NCBI table 11 differs from table 1 only in its starts).
CODE: dict[str, str] = {a + b + c: translate(a + b + c, False)
                        for a in _BASES for b in _BASES for c in _BASES}
CODONS_FOR: dict[str, list[str]] = {}
for _codon, _aa in CODE.items():
    CODONS_FOR.setdefault(_aa, []).append(_codon)
#: Start codons of table 11 (the set Kleborate also accepts for mgrB).
START_CODONS = frozenset({"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"})

AA3 = {"A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe", "G": "Gly", "H": "His",
       "I": "Ile", "K": "Lys", "L": "Leu", "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln",
       "R": "Arg", "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr", "*": "Ter",
       "X": "Xaa", "U": "Sec"}

#: Kleborate v3 locus names, exactly as its KpSC AMR module prints them.
KLEBORATE_LOCI = {"gyrA": "GyrA", "parC": "ParC", "ompK35": "OmpK35", "ompK36": "OmpK36",
                  "mgrB": "mgrb", "pmrB": "pmrb"}
#: QRDR positions Kleborate reports (klebsiella_pneumo_complex__amr/qrdr_mutations.py).
KLEBORATE_QRDR = {"gyrA": (83, 87), "parC": (80, 84)}
#: The AMRFinderPlus taxgroup of the K. pneumoniae species complex, the only one
#: Kleborate's AMR module types.
KLEBORATE_ORGANISM = "Klebsiella_pneumoniae"

#: An inserted or deleted stretch longer than this is reported as an
#: interruption or a large deletion rather than a small frameshifting indel.
LARGE_INDEL = 30

TYPE_LABELS = {
    "substitution": "Amino-acid substitution",
    "nonsense": "Nonsense (premature stop codon)",
    "insertion": "In-frame insertion",
    "deletion": "In-frame deletion",
    "rRNA": "rRNA base change",
    "promoter": "Promoter / upstream base change",
    "nucleotide": "Nucleotide change",
    "dna_insertion": "Nucleotide insertion",
    "dna_deletion": "Nucleotide deletion",
    "disruption_stop": "Gene disruption: premature stop codon",
    "disruption_frameshift": "Gene disruption: frameshift",
    "disruption_interrupted": "Gene disruption: interrupted by inserted sequence",
    "disruption_deletion": "Gene disruption: internal deletion",
    "disruption_start": "Gene disruption: start codon lost",
    "disruption_5prime": "Gene disruption: 5' part of the gene absent",
    "disruption_unresolved": "Gene disruption: lesion not resolved",
    "contig_end": ("Gene cut by a contig end: NOT evidence of disruption, "
                   "loss of function not assessed"),
    "divergent": "Divergent (mosaic) allele",
}


# ------------------------------------------------------------------ helpers
def resistance_text(drug_class: str, subclass: str, *, loss_of_function: bool = False) -> str:
    """'quinolone resistance', 'carbapenem resistance (beta-lactam class)' ..."""
    cls = (drug_class or "").strip()
    sub = (subclass or "").strip()
    if not cls and not sub:
        return ""

    def words(text: str) -> str:
        return text.replace("/", ", ").lower()

    if sub and cls and sub.upper() != cls.upper():
        text = f"{words(sub)} resistance ({words(cls)} class)"
    else:
        text = f"{words(sub or cls)} resistance"
    return text + ("; loss of function" if loss_of_function else "")


#: What a contig-end row is marked as, in its own columns. Defined in
#: :mod:`hydra_amr.records`, next to the two fields they are written into, and
#: named here so a reader of this module can see what disclose_contig_end()
#: sets without going to look.
_CONTIG_END_COLUMNS = (CONTIG_END_SUBTYPE, CONTIG_END_RESOLUTION)


def disclose_contig_end(hit: Hit) -> bool:
    """Mark a row whose gene turned out to be cut by a contig end.

    Called after the lesion has been read from the assembly, which is the first
    moment anything knows that the "disruption" is an assembly boundary. Three
    things have to be visible: that the gene runs off the end of the contig,
    that this is not evidence of disruption, and that loss of function could not
    be assessed. The resistance field stops asserting a loss of function that
    was never established, and says what the finding would mean if it were.

    Returns whether the row was marked, so a caller can count them.
    """
    detail = hit.mutation
    if detail is None or detail.call_type != "contig_end":
        return False
    hit.element_subtype = CONTIG_END_SUBTYPE
    hit.resolution = CONTIG_END_RESOLUTION
    disclosure = (f"{hit.gene} runs off the end of contig {hit.sequence}, so the gene is cut "
                  f"by the assembly: this is NOT evidence of disruption and loss of function "
                  f"could not be assessed. Reported rather than suppressed, but it is not a "
                  f"confirmed loss-of-function call")
    hit.note = f"{disclosure}. The call it was reported from: {hit.note}" if hit.note \
        else disclosure
    would_be = resistance_text(hit.drug_class, hit.subclass)
    detail.resistance = ("loss of function not established: the gene is cut by a contig end"
                         + (f", so the {would_be} a confirmed loss of function would mean "
                            f"cannot be reported" if would_be else ""))
    return True


def infer_reference_codon(observed: str, ref_aa: str) -> list[str]:
    """Codons for *ref_aa* closest to *observed* (fewest base differences).

    A parsimony inference, not an observation: the caller labels it as such.
    """
    candidates = CODONS_FOR.get(ref_aa.upper(), [])
    if not candidates or len(observed) != 3:
        return []
    best = min(_differences(c, observed) for c in candidates)
    return sorted(c for c in candidates if _differences(c, observed) == best)


def _differences(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def base_changes(reference: str, observed: str, first_cds_position: int) -> tuple[str, str]:
    """('c.248C>T', 'C>T at codon position 2') for one codon, HGVS-style."""
    diffs = [i for i, (x, y) in enumerate(zip(reference, observed)) if x != y]
    if not diffs:
        return "", ""
    if len(diffs) == 1:
        i = diffs[0]
        return (f"c.{first_cds_position + i}{reference[i]}>{observed[i]}",
                f"{reference[i]}>{observed[i]} at codon position {i + 1}")
    lo, hi = diffs[0], diffs[-1]
    return (f"c.{first_cds_position + lo}_{first_cds_position + hi}delins{observed[lo:hi + 1]}",
            "; ".join(f"{reference[i]}>{observed[i]} at codon position {i + 1}" for i in diffs))


class ContigReader:
    """A contig read in a gene's orientation, with 1-based contig coordinates.

    *top* is the contig coordinate of a base as the gene reads it: on the minus
    strand the gene runs towards lower coordinates and every base is complemented.
    """

    def __init__(self, sequence: str, strand: str):
        self.seq = sequence.upper()
        self.minus = strand == "-"
        self.length = len(sequence)

    def step(self, top: int, k: int) -> int:
        """The coordinate *k* bases further along the gene from *top*."""
        return top - k if self.minus else top + k

    def read(self, top: int, n: int) -> str:
        """Up to *n* bases from *top* onwards in gene orientation (shorter at a contig end)."""
        if n <= 0:
            return ""
        if self.minus:
            lo, hi = max(1, top - n + 1), min(top, self.length)
            return revcomp(self.seq[lo - 1:hi]) if hi >= lo else ""
        lo, hi = max(1, top), min(top + n - 1, self.length)
        return self.seq[lo - 1:hi] if hi >= lo else ""

    def codon(self, top: int) -> str | None:
        """The codon whose first base is at *top*, or None when it leaves the contig."""
        lo, hi = (top - 2, top) if self.minus else (top, top + 2)
        if lo < 1 or hi > self.length:
            return None
        return self.read(top, 3)

    def span(self, top: int, n: int) -> tuple[int, int]:
        """(low, high) contig coordinates of *n* bases from *top* in gene orientation."""
        end = self.step(top, n - 1)
        return (min(top, end), max(top, end))

    def remaining(self, top: int) -> int:
        """Bases left on the contig from *top* onwards in gene orientation, *top* included."""
        return top if self.minus else self.length - top + 1


def _location(contig: str, lo: int, hi: int, strand: str) -> str:
    return f"{contig}:{lo}-{hi} ({strand})" if lo != hi else f"{contig}:{lo} ({strand})"


@dataclass
class _Column:
    subject: int | None   # reference residue/base number, None in a subject gap
    query: int | None     # 0-based index of the query residue/base, None in a query gap
    q: str
    s: str


def alignment_columns(hsp, *, nucleotide: bool) -> list[_Column]:
    """Every column of a BLAST alignment with its subject number and query index.

    Subject numbers descend for a nucleotide alignment to the reverse strand of
    the subject; query indexes count query residues (codons for blastx, bases
    for blastn) from the alignment's first query position.
    """
    step = 1 if hsp.send >= hsp.sstart else -1
    if not nucleotide:
        step = 1
    columns: list[_Column] = []
    subject = hsp.sstart
    query = 0
    for q, s in zip(hsp.qseq, hsp.sseq):
        columns.append(_Column(subject if s != "-" else None, query if q != "-" else None, q, s))
        if s != "-":
            subject += step
        if q != "-":
            query += 1
    return columns


def _protein_frame(hsp, offset: int) -> tuple[str, int]:
    """(strand, contig coordinate of the first base of the first aligned codon)."""
    strand = "+" if hsp.qend >= hsp.qstart else "-"
    return strand, hsp.qstart + offset


# --------------------------------------------------------- protein mutations
def describe_protein(hit: Hit, evidence: dict, contig: str | None,
                     organism: str | None) -> MutationDetail:
    """A catalogued protein mutation (POINTX): codon, base change and location."""
    entry = evidence["entry"]
    hsp = evidence["hsp"]
    position = evidence["position"]
    strand, first_top = _protein_frame(hsp, evidence["offset"])
    detail = MutationDetail(symbol=entry.symbol,
                            resistance=resistance_text(hit.drug_class, hit.subclass))
    detail.gene_location = _location(hit.sequence, min(hit.start, hit.end),
                                     max(hit.start, hit.end), strand)
    change = entry.change_type
    ref, alt = entry.ref.upper(), entry.alt.upper()
    gene = entry.gene or hit.gene
    if contig is None:
        detail.call_type = {"nonsense": "nonsense", "deletion": "deletion",
                            "insertion": "insertion"}.get(change, "substitution")
        detail.type_label = TYPE_LABELS[detail.call_type]
        detail.unresolved = f"contig {hit.sequence} was not found in the assembly file"
        return detail
    reader = ContigReader(contig, strand)
    columns = alignment_columns(hsp, nucleotide=False)
    by_subject = {c.subject: c for c in columns if c.subject is not None}

    def codon_top(query_index: int) -> int:
        return reader.step(first_top, 3 * query_index)

    def read_codon(query_index: int, expected: str) -> tuple[str | None, str]:
        top = codon_top(query_index)
        codon = reader.codon(top)
        if codon is None:
            return None, "the codon runs past the end of the contig"
        if CODE.get(codon, "X") != expected.upper() and expected.upper() != "X":
            return codon, (f"the codon read from the assembly ({codon}) does not encode the "
                           f"aligned residue {expected}")
        return codon, ""

    if change in ("substitution", "nonsense"):
        detail.call_type = change
        detail.type_label = TYPE_LABELS[change]
        detail.protein_change = f"{ref}{position}{'Ter' if alt == '*' else alt}"
        detail.protein_hgvs = f"p.{AA3.get(ref, ref)}{position}{AA3.get(alt, alt)}"
        column = by_subject.get(position)
        if column is None or column.query is None:
            detail.unresolved = (f"reference residue {position} falls in a gap of the alignment; "
                                 "there is no codon to read")
            return detail
        codon, problem = read_codon(column.query, column.q)
        if codon is None:
            detail.unresolved = problem
            return detail
        top = codon_top(column.query)
        lo, hi = reader.span(top, 3)
        detail.contig, detail.start, detail.end, detail.strand = hit.sequence, lo, hi, strand
        detail.location = _location(hit.sequence, lo, hi, strand)
        detail.observed = codon
        if problem:
            detail.unresolved = problem
        observed_aa = CODE.get(codon, "X")
        if observed_aa != alt:
            detail.unresolved = "; ".join(filter(None, [
                detail.unresolved,
                f"the codon {codon} encodes {observed_aa}, not the catalogued {alt}"]))
        candidates = infer_reference_codon(codon, ref)
        cds = 3 * (position - 1) + 1
        if len(candidates) == 1:
            detail.reference_inferred = candidates[0]
            change_text, detail.base_change = base_changes(candidates[0], codon, cds)
            detail.nucleotide_change = f"{change_text} (inferred)"
            detail.nucleotide_basis = ("inferred: the database has no reference coding sequence "
                                       f"for this protein; {candidates[0]} is the only {ref} "
                                       f"codon one step from {codon}")
        elif candidates:
            detail.reference_inferred = " or ".join(candidates)
            detail.nucleotide_basis = (f"not determinable: {len(candidates)} {ref} codons are "
                                       f"equally close to {codon}, and the database has no "
                                       "reference coding sequence")
        detail.description = (
            f"{gene} codon {position} reads {codon} ({AA3.get(observed_aa, observed_aa)}) in this "
            f"isolate where the reference has {AA3.get(ref, ref)}"
            + (f" (reference codon {detail.reference_inferred}, inferred)"
               if len(candidates) == 1 else "")
            + f": {detail.protein_hgvs}, {detail.resistance or 'catalogued resistance mutation'}.")
    elif change == "deletion":
        detail.call_type = "deletion"
        detail.type_label = TYPE_LABELS["deletion"]
        n = max(1, len(ref))
        end_position = position + n - 1
        if n == 1:
            detail.protein_change = f"{ref}{position}del"
            detail.protein_hgvs = f"p.{AA3.get(ref, ref)}{position}del"
        else:
            detail.protein_change = f"{ref}{position}del"
            detail.protein_hgvs = (f"p.{AA3.get(ref[0], ref[0])}{position}_"
                                   f"{AA3.get(ref[-1], ref[-1])}{end_position}del")
        detail.nucleotide_change = f"c.{3 * (position - 1) + 1}_{3 * end_position}del"
        detail.nucleotide_basis = ("positions from the catalogue numbering; the deleted bases "
                                   "themselves are not in the assembly and the database has no "
                                   "reference coding sequence")
        detail.observed = "deleted"
        before = by_subject.get(position - 1)
        after = by_subject.get(end_position + 1)
        if before and after and before.query is not None and after.query is not None:
            last = reader.step(codon_top(before.query), 2)
            nxt = codon_top(after.query)
            lo, hi = min(last, nxt), max(last, nxt)
            detail.contig, detail.start, detail.end, detail.strand = hit.sequence, lo, hi, strand
            detail.location = (f"{hit.sequence}: between {last} and {nxt} ({strand}), "
                               "the codons either side of the deletion")
            detail.observed = (f"deleted ({n} codon{'s' if n > 1 else ''}, {3 * n} bp); "
                               f"flanking codons {reader.codon(codon_top(before.query))}"
                               f"|{reader.codon(nxt)}")
        else:
            detail.unresolved = "the residues flanking the deletion are not both aligned"
        detail.description = (f"{n} residue{'s' if n > 1 else ''} ({ref}) from position "
                              f"{position} of {gene} {'are' if n > 1 else 'is'} missing in this "
                              f"isolate (in-frame deletion): {detail.resistance}.")
    elif change == "insertion":
        detail.call_type = "insertion"
        detail.type_label = TYPE_LABELS["insertion"]
        inserted = entry.inserted_bases
        detail.protein_change = entry.symbol.rpartition("_")[2]
        located = _locate_insertion(columns, position, inserted)
        nxt = by_subject.get(position + 1)
        next_aa = (nxt.s.upper() if nxt else "")
        detail.protein_hgvs = (f"p.{AA3.get(ref, ref)}{position}_{AA3.get(next_aa, next_aa)}"
                               f"{position + 1}ins{''.join(AA3.get(a, a) for a in inserted)}"
                               if next_aa else f"p.{position}_{position + 1}ins"
                               f"{''.join(AA3.get(a, a) for a in inserted)}")
        if located is None:
            detail.unresolved = "the inserted residues could not be located in the alignment"
        else:
            anchor_index, first_inserted = located
            codons = [reader.codon(codon_top(first_inserted + i)) for i in range(len(inserted))]
            if any(c is None for c in codons):
                detail.unresolved = "the inserted codons run past the end of the contig"
            else:
                bases = "".join(codons)
                translated = "".join(CODE.get(c, "X") for c in codons)
                if translated != inserted:
                    detail.unresolved = (f"the inserted codons read {bases} ({translated}), "
                                         f"not {inserted}")
                top = codon_top(first_inserted)
                lo, hi = reader.span(top, 3 * len(inserted))
                detail.contig, detail.start, detail.end, detail.strand = hit.sequence, lo, hi, strand
                detail.location = _location(hit.sequence, lo, hi, strand)
                detail.observed = " ".join(codons)
                cds = 3 * position
                detail.nucleotide_change = f"c.{cds}_{cds + 1}ins{bases}"
                detail.nucleotide_basis = ("inserted bases read from the assembly; placed after "
                                           f"codon {position} as the catalogue names it (in a "
                                           "repeat the exact point is ambiguous)")
                if organism == KLEBORATE_ORGANISM and gene == "ompK36":
                    detail.kleborate = _kleborate_ompk36_insertion(reader, first_top, columns,
                                                                   codon_top)
        detail.description = (f"{len(inserted)} residue{'s' if len(inserted) > 1 else ''} "
                              f"({inserted}) inserted after {AA3.get(ref, ref)}{position} of "
                              f"{gene}: {detail.resistance}.")
    else:
        detail.call_type = "substitution"
        detail.type_label = "Complex change"
        detail.unresolved = f"catalogue change {entry.ref}->{entry.alt} is not a simple substitution"
    if not detail.kleborate and organism == KLEBORATE_ORGANISM:
        detail.kleborate = _kleborate_point(gene, change, ref, alt, position)
    return detail


def _locate_insertion(columns: list[_Column], position: int, inserted: str):
    """(query index of the anchor residue, query index of the first inserted residue).

    The catalogue names an insertion after the 3'-most residue that describes it
    (HGVS); the aligner may have placed the gap further left inside a repeat, so
    each inserted block is shifted right exactly as :mod:`.mutations` does before
    comparing (``_shift_insertions_right``). Shifting a block one place moves the
    inserted window one query residue to the right.
    """
    subject_residue: dict[int, str] = {}
    blocks: list[tuple[int, list[int], str]] = []
    pending_q: list[int] = []
    pending_r: list[str] = []
    last = None
    for column in columns:
        if column.subject is None:
            if column.query is not None:
                pending_q.append(column.query)
                pending_r.append(column.q.upper())
            continue
        if pending_q and last is not None:
            blocks.append((last, pending_q, "".join(pending_r)))
        pending_q, pending_r = [], []
        subject_residue[column.subject] = column.s.upper()
        last = column.subject
    for anchor, queries, block in blocks:
        shift = 0
        while block and subject_residue.get(anchor + 1) == block[0]:
            anchor += 1
            block = block[1:] + block[0]
            shift += 1
        if anchor != position or block != inserted.upper():
            continue
        first = queries[0] + shift
        return first - 1, first
    return None


def _kleborate_point(gene: str, change: str, ref: str, alt: str, position: int) -> str:
    if gene in KLEBORATE_QRDR and change == "substitution" and position in KLEBORATE_QRDR[gene]:
        return f"{KLEBORATE_LOCI[gene]}:p.{AA3.get(ref, ref)}{position}{AA3.get(alt, alt)}"
    if gene in ("ompK35", "ompK36", "mgrB", "pmrB") and change == "nonsense" and position > 1:
        # Kleborate names any truncation after the last residue before the break.
        return ""
    return ""


def _kleborate_ompk36_insertion(reader: ContigReader, first_top: int, columns, codon_top) -> str:
    """Kleborate's OmpK36 loop-3 call, computed the way Kleborate computes it.

    Kleborate translates the gene and looks for GDGDTY (GD) or GDTDTY (TD); the
    insertion is numbered from that motif in the isolate's own translation.
    """
    n_query = max((c.query for c in columns if c.query is not None), default=-1) + 1
    protein = "".join(CODE.get(reader.codon(codon_top(i)) or "", "X") for i in range(n_query))
    for motif, kind, name in (("GDGDTY", "GD", "GlyAsp"), ("GDTDTY", "TD", "ThrAsp")):
        index = protein.find(motif)
        if index >= 0:
            # Kleborate numbers from its own translation, which starts at the
            # gene's first codon; an alignment starting later is offset by it.
            first_subject = next((c.subject for c in columns if c.subject is not None), 1)
            start = index + motif.index(kind) + first_subject
            return f"OmpK36:p.{start}_{start + 1}ins{name}"
    return ""


# ------------------------------------------------------------- DNA mutations
def describe_dna(hit: Hit, evidence: dict, contig: str | None) -> MutationDetail:
    """A catalogued nucleotide mutation (POINTN): rRNA, promoter and similar loci."""
    entry = evidence["entry"]
    hsp = evidence["hsp"]
    position = evidence["position"]          # 1-based in the reference record
    reference_seq = (evidence.get("reference") or "").upper()
    offset = evidence["offset"]
    detail = MutationDetail(symbol=entry.symbol,
                            resistance=resistance_text(hit.drug_class, hit.subclass))
    minus = hsp.send < hsp.sstart
    strand = "-" if minus else "+"
    detail.gene_location = _location(hit.sequence, min(hit.start, hit.end),
                                     max(hit.start, hit.end), strand)
    locus = (hsp.sseqid.split("@")[1] if "@" in hsp.sseqid else entry.gene).replace("_", " ")
    rrna = any(tag in (entry.gene + hsp.sseqid) for tag in ("16S", "23S", "rRNA", "rrs", "rrl"))
    promoter = entry.position < 0 or "promoter" in hsp.sseqid.lower()
    change = entry.change_type
    if change == "insertion":
        detail.call_type = "dna_insertion"
    elif change == "deletion":
        detail.call_type = "dna_deletion"
    else:
        detail.call_type = "rRNA" if rrna else ("promoter" if promoter else "nucleotide")
    detail.type_label = TYPE_LABELS[detail.call_type]
    number = f"{entry.position}"
    prefix = "n." if rrna else "c."
    ref, alt = entry.ref.upper(), entry.alt.upper()
    if change == "substitution":
        detail.nucleotide_change = f"{prefix}{number}{ref}>{alt}"
    elif change == "insertion":
        detail.nucleotide_change = f"{prefix}{number}_{entry.position + 1}ins{entry.inserted_bases}"
    elif change == "deletion":
        detail.nucleotide_change = f"{prefix}{number}del{ref}"
    detail.protein_change = ""
    if contig is None:
        detail.unresolved = f"contig {hit.sequence} was not found in the assembly file"
        return detail
    columns = alignment_columns(hsp, nucleotide=True)
    column = next((c for c in columns if c.subject == position), None)
    if column is None:
        detail.unresolved = f"reference base {position} is not covered by the alignment"
        return detail
    if 1 <= position <= len(reference_seq):
        detail.reference = reference_seq[position - 1]
        detail.nucleotide_basis = (f"reference base from the catalogue's reference sequence "
                                   f"({hsp.sseqid.split('@')[0]}, {locus})")
    else:
        aligned_ref = column.s.upper()
        detail.reference = revcomp(aligned_ref) if minus else aligned_ref
        detail.nucleotide_basis = "reference base from the aligned reference sequence"
    if column.query is None:
        detail.observed = "deleted"
        detail.unresolved = ("" if change == "deletion"
                             else f"reference base {position} falls in a gap of the alignment")
        return detail
    coordinate = hsp.qstart + column.query + offset       # blastn queries always ascend
    base = contig.upper()[coordinate - 1] if 1 <= coordinate <= len(contig) else ""
    observed = revcomp(base) if minus else base
    aligned = (revcomp(column.q.upper()) if minus else column.q.upper())
    detail.observed = observed
    detail.contig, detail.start, detail.end, detail.strand = (hit.sequence, coordinate,
                                                              coordinate, strand)
    detail.location = _location(hit.sequence, coordinate, coordinate, strand)
    if observed != aligned:
        detail.unresolved = (f"the base read from the assembly ({observed}) differs from the "
                             f"aligned base ({aligned})")
    if change == "substitution":
        detail.base_change = f"{detail.reference}>{observed}"
        if observed != alt:
            detail.unresolved = "; ".join(filter(None, [
                detail.unresolved, f"the isolate carries {observed}, not the catalogued {alt}"]))
    elif change == "insertion":
        detail.base_change = f"{entry.inserted_bases} inserted after {detail.reference}"
    where = (f"{abs(entry.position)} bases upstream of the {entry.gene} start"
             if entry.position < 0 else f"position {entry.position} of {locus}")
    detail.description = (f"{locus}: the isolate carries {observed} at {where} where the "
                          f"reference has {detail.reference} ({detail.nucleotide_change}): "
                          f"{detail.resistance or 'catalogued resistance mutation'}.")
    return detail


# ---------------------------------------------------------------- disruptions
@dataclass
class _Piece:
    """One alignment of the reference protein, in contig coordinates."""

    contig: str
    strand: str
    first_top: int          # contig coordinate of the first base of the first codon
    sub_lo: int
    sub_hi: int
    columns: list
    identity: float

    def top_of_query(self, index: int) -> int:
        return self.first_top - 3 * index if self.strand == "-" else self.first_top + 3 * index

    def subject_tops(self) -> dict[int, int]:
        """Reference residue -> contig coordinate of its codon's first base."""
        return {c.subject: self.top_of_query(c.query) for c in self.columns
                if c.subject is not None and c.query is not None}

    @property
    def gene_lo(self) -> int:
        last = self.top_of_query(max((c.query for c in self.columns if c.query is not None),
                                     default=0))
        ends = (self.first_top, last - 2 if self.strand == "-" else last + 2)
        return min(ends)

    @property
    def gene_hi(self) -> int:
        last = self.top_of_query(max((c.query for c in self.columns if c.query is not None),
                                     default=0))
        ends = (self.first_top, last - 2 if self.strand == "-" else last + 2)
        return max(ends)


def _pieces(evidence: dict) -> list[_Piece]:
    out = []
    for hsp, piece in evidence.get("pieces") or [(evidence["hsp"], None)]:
        offset = piece.offset if piece is not None else evidence["offset"]
        contig = piece.contig if piece is not None else evidence.get("contig", "")
        strand, first_top = _protein_frame(hsp, offset)
        out.append(_Piece(contig=contig, strand=strand, first_top=first_top,
                          sub_lo=min(hsp.sstart, hsp.send), sub_hi=max(hsp.sstart, hsp.send),
                          columns=alignment_columns(hsp, nucleotide=False), identity=hsp.pident))
    return out


def _translate_from(reader: ContigReader, top: int, max_codons: int) -> tuple[str, bool]:
    """Translate codon by codon until a stop; (protein incl. '*', reached a stop)."""
    protein = []
    for i in range(max_codons):
        codon = reader.codon(reader.step(top, 3 * i))
        if codon is None:
            return "".join(protein), False
        aa = CODE.get(codon, "X")
        protein.append(aa)
        if aa == "*":
            return "".join(protein), True
    return "".join(protein), False


def describe_disruption(hit: Hit, evidence: dict, contigs: Mapping[str, str],
                        organism: str | None) -> MutationDetail:
    """The lesion behind a loss-of-function call, read from the assembly."""
    gene = hit.gene
    reference = (evidence.get("reference") or "").replace("*", "").upper()
    length = len(reference) or int(evidence["hsp"].slen or 0)
    detail = MutationDetail(resistance=resistance_text(hit.drug_class, hit.subclass,
                                                       loss_of_function=True))
    pieces = [p for p in _pieces(evidence) if p.contig == hit.sequence]
    best_strand, _ = _protein_frame(evidence["hsp"], evidence["offset"])
    detail.gene_location = _location(hit.sequence, min(hit.start, hit.end),
                                     max(hit.start, hit.end), best_strand)
    contig = contigs.get(hit.sequence)
    if contig is None or not reference:
        detail.call_type = "disruption_unresolved"
        detail.type_label = TYPE_LABELS[detail.call_type]
        detail.symbol = gene
        detail.unresolved = ("contig not found in the assembly file" if contig is None
                             else "reference protein not available")
        return detail
    reader = ContigReader(contig, best_strand)
    # The pieces of this copy of the gene: same contig and strand, near the call.
    window = 3 * length * 2 + 5000
    centre = (min(hit.start, hit.end) + max(hit.start, hit.end)) // 2
    local = sorted((p for p in pieces if p.strand == best_strand
                    and abs((p.gene_lo + p.gene_hi) // 2 - centre) <= window),
                   key=lambda p: (p.sub_lo, -p.identity))
    elsewhere = [p for p in _pieces(evidence) if p not in local]
    if not local:
        local = [p for p in _pieces(evidence) if p.contig == hit.sequence][:1]
    first = local[0]
    tops = first.subject_tops()
    start_top = None
    if first.sub_lo == 1 and 1 in tops:
        start_top = tops[1]
    elif 1 < first.sub_lo <= 40:
        start_top = _find_start_upstream(reader, first, reference)

    # The gene as the aligned pieces place it, which for a broken gene is wider
    # than the single alignment the call was reported from.
    lo = min([p.gene_lo for p in local] + ([start_top] if start_top else []))
    hi = max([p.gene_hi for p in local] + ([start_top] if start_top else []))
    detail.gene_location = _location(hit.sequence, lo, hi, best_strand) + (
        f", {len(local)} aligned parts" if len(local) > 1 else "")
    if start_top is None:
        return _describe_missing_5prime(hit, detail, reader, first, local, elsewhere,
                                        reference, length, organism)

    start_codon = reader.codon(start_top) or ""
    protein, stopped = _translate_from(reader, start_top, 3 * length + 50)
    body = protein[:-1] if stopped else protein

    # Register: isolate residue i (0-based, from the start codon) <-> reference
    # residue, following the first alignment's columns (they carry any in-frame
    # indel) and continuing one-for-one past its end.
    register = _register(first, start_top, reader)

    if start_codon not in START_CODONS:
        detail.call_type = "disruption_start"
        detail.type_label = TYPE_LABELS["disruption_start"]
        detail.symbol = f"{gene}_M1?"
        detail.protein_change = "M1?"
        detail.protein_hgvs = "p.Met1?"
        detail.observed = start_codon
        lo, hi = reader.span(start_top, 3)
        _set_location(detail, hit.sequence, lo, hi, best_strand)
        detail.description = (f"The {gene} start codon reads {start_codon} in this isolate, "
                              "which is not a start codon: no protein is expected to be made "
                              f"from it ({detail.resistance}).")
        if organism == KLEBORATE_ORGANISM and gene == "mgrB":
            detail.kleborate = "mgrB$"
        return detail

    if not stopped:
        if reader.remaining(start_top) < 3 * length:
            return _contig_end(hit, detail, reader, start_top, length, len(body), "3'")
        detail.call_type = "disruption_unresolved"
        detail.type_label = TYPE_LABELS[detail.call_type]
        detail.symbol = gene
        detail.unresolved = ("translating from the start codon reads past the reference length "
                             "without a stop; no disruptive lesion identified")
        return detail

    # First reference residue the isolate does not reproduce, scanning in
    # register from the end of the first alignment (where BLAST stopped
    # extending) -- the HGVS 'first amino acid changed'.
    stop_index = len(body)                      # isolate residue index of the stop
    stop_ref = register(stop_index)             # reference residue the stop replaces
    if stop_ref is not None and stop_ref > 0.95 * length:
        detail.call_type = "disruption_unresolved"
        detail.type_label = TYPE_LABELS[detail.call_type]
        detail.symbol = gene
        detail.unresolved = (f"the isolate's reading frame from the start codon runs to residue "
                             f"{stop_ref} of {length} and stops normally; the call rests on the "
                             "alignment alone and no disruptive lesion is visible in the sequence")
        return detail

    # Where the first alignment ends, as an isolate residue index; a first
    # alignment in another frame than the start codon (an N-terminal
    # frameshift) maps to nothing and the scan starts at the start codon.
    scan_from = 0
    last_query = max((c.query for c in first.columns if c.query is not None), default=0)
    distance = abs(first.top_of_query(last_query) - start_top)
    if distance % 3 == 0:
        scan_from = distance // 3
    divergence = None
    if stop_index < scan_from:
        # A stop inside the first alignment: blastx aligned straight through it.
        divergence = stop_index
    else:
        for i in range(scan_from, stop_index + 1):
            r = register(i)
            if r is None or r > length:
                divergence = i
                break
            aa = protein[i] if i < len(protein) else ""
            if aa != reference[r - 1]:
                divergence = i
                break
    if divergence is None:
        divergence = stop_index
    d_ref = register(divergence) or (divergence + 1)
    ref_aa = reference[d_ref - 1] if 0 < d_ref <= length else "X"
    new_aa = protein[divergence] if divergence < len(protein) else "*"
    intact = d_ref - 1
    detail.intact = f"residues 1-{intact} of {length} ({100.0 * intact / length:.0f}%)"

    # What lies beyond: the next piece of the gene in a different frame (frameshift),
    # far away (interrupted) or nowhere (C-terminal part absent / novel sequence).
    downstream = next((p for p in local[1:] if p.sub_lo > first.sub_hi - 5), None)
    net = None
    junction = None
    if downstream is not None:
        last_ref = first.sub_hi
        last_top = first.subject_tops().get(last_ref)
        next_ref = downstream.sub_lo
        next_top = downstream.subject_tops().get(next_ref)
        if last_top is not None and next_top is not None:
            last_end = reader.step(last_top, 2)
            between = (next_top - last_end - 1) if best_strand == "+" else (last_end - next_top - 1)
            net = between - 3 * (next_ref - last_ref - 1)
            junction = (last_ref, last_end, next_ref, next_top)

    nonsense = new_aa == "*" and _continues_in_frame(reader, start_top, divergence, register,
                                                     reference)
    if nonsense:
        detail.call_type = "disruption_stop"
        detail.protein_change = f"{ref_aa}{d_ref}Ter"
        detail.protein_hgvs = f"p.{AA3.get(ref_aa, ref_aa)}{d_ref}Ter"
        top = reader.step(start_top, 3 * divergence)
        codon = reader.codon(top) or ""
        lo, hi = reader.span(top, 3)
        _set_location(detail, hit.sequence, lo, hi, best_strand)
        detail.observed = codon
        candidates = infer_reference_codon(codon, ref_aa)
        if len(candidates) == 1:
            detail.reference_inferred = candidates[0]
            change_text, detail.base_change = base_changes(
                candidates[0], codon, 3 * (d_ref - 1) + 1)
            detail.nucleotide_change = f"{change_text} (inferred)"
            detail.nucleotide_basis = ("inferred: no reference coding sequence in the database; "
                                       f"{candidates[0]} is the only {ref_aa} codon one step "
                                       f"from {codon}")
        elif candidates:
            detail.reference_inferred = " or ".join(candidates)
            detail.nucleotide_basis = "not determinable: several reference codons equally close"
        detail.description = (f"{gene} codon {d_ref} ({AA3.get(ref_aa, ref_aa)} in the reference) "
                              f"is a stop codon ({codon}) in this isolate: the protein ends after "
                              f"{intact} of {length} residues ({100.0 * intact / length:.0f}%), "
                              f"{detail.resistance}.")
    else:
        # AMRFinderPlus counts the residues of the new frame before the stop
        # (cirA_Y253CfsTer5, ompK35_Y108TerfsTer0); HGVS counts the stop itself
        # as the last position, so the long form is one higher.
        novel = stop_index - divergence
        if new_aa == "*":
            detail.protein_change = f"{ref_aa}{d_ref}TerfsTer0"
            detail.protein_hgvs = f"p.{AA3.get(ref_aa, ref_aa)}{d_ref}Ter"
        else:
            detail.protein_change = f"{ref_aa}{d_ref}{new_aa}fsTer{novel}"
            detail.protein_hgvs = (f"p.{AA3.get(ref_aa, ref_aa)}{d_ref}"
                                   f"{AA3.get(new_aa, new_aa)}fsTer{novel + 1}")
        ter = novel + 1
        stop_top = reader.step(start_top, 3 * stop_index)
        if net is not None and abs(net) > LARGE_INDEL:
            detail.call_type = ("disruption_interrupted" if net > 0 else "disruption_deletion")
        else:
            detail.call_type = "disruption_frameshift"
        if junction is not None:
            last_ref, last_end, next_ref, next_top = junction
            lo, hi = min(last_end, next_top), max(last_end, next_top)
            _set_location(detail, hit.sequence, lo, hi, best_strand)
        else:
            top = reader.step(start_top, 3 * divergence)
            lo, hi = reader.span(top, 3)
            _set_location(detail, hit.sequence, lo, hi, best_strand)
        stop_codon = reader.codon(stop_top) or ""
        stop_lo, stop_hi = reader.span(stop_top, 3)
        detail.observed = (f"new frame from codon {d_ref} stops at {stop_codon} "
                           f"({hit.sequence}:{stop_lo}-{stop_hi})")
        if detail.call_type == "disruption_frameshift":
            if net is not None and net % 3:
                kind = "deletion" if net < 0 else "insertion"
                detail.base_change = f"net {abs(net)}-bp {kind}"
                inferred = _infer_small_indel(reader, first, downstream, net, reference)
                if len(inferred) == 1:
                    notation, where, span = inferred[0]
                    detail.nucleotide_change = f"{notation} (inferred)"
                    detail.nucleotide_basis = ("inferred: the database has no reference coding "
                                               "sequence; this is the only base change that "
                                               "restores the reference protein across the "
                                               f"junction ({where})")
                    _set_location(detail, hit.sequence, span[0], span[1], best_strand)
                elif inferred:
                    detail.nucleotide_change = " or ".join(n for n, _w, _s in inferred) + \
                        " (inferred, ambiguous)"
                    detail.nucleotide_basis = (
                        "inferred, ambiguous: the database has no reference coding sequence and "
                        f"{len(inferred)} different base changes restore the reference protein "
                        "across the junction ("
                        + "; ".join(f"{n}: {w}" for n, w, _s in inferred) + ")")
                else:
                    detail.nucleotide_basis = ("the exact bases are not determinable without a "
                                               "reference coding sequence")
                event = (f"a net {abs(net)}-bp {kind} between reference codons "
                         f"{junction[0]} and {junction[2]}")
            elif net is not None:
                event = (f"a net {abs(net)}-bp change between codons {junction[0]} and "
                         f"{junction[2]} that keeps the frame of the downstream part")
                detail.unresolved = ("the downstream part is in frame with the upstream part; "
                                     "the early stop lies in sequence between them")
            else:
                event = "a frameshift (the downstream part of the gene does not align, so the " \
                        "size of the indel is not measurable)"
                detail.base_change = "indel size not measurable"
            detail.description = (f"{gene} is frameshifted by {event}: the reading frame from the "
                                  f"start codon follows the reference for {intact} of {length} "
                                  f"residues ({100.0 * intact / length:.0f}%), then reads "
                                  f"{ter - 1} novel residue{'s' if ter - 1 != 1 else ''} and "
                                  f"stops ({detail.protein_hgvs}); {detail.resistance}.")
        elif detail.call_type == "disruption_interrupted":
            detail.base_change = f"{net} bp inserted"
            detail.nucleotide_basis = ("inserted sequence measured from the two aligned parts of "
                                       "the gene; the inserted element is not identified (no "
                                       "insertion-sequence database)")
            detail.description = (f"{gene} is interrupted by about {net} bp of inserted sequence "
                                  f"between reference codons {junction[0]} and {junction[2]} "
                                  f"({detail.location}); the protein made from the start codon "
                                  f"keeps {intact} of {length} residues "
                                  f"({100.0 * intact / length:.0f}%) and runs into the insert "
                                  f"({detail.protein_hgvs}); {detail.resistance}.")
        else:
            detail.base_change = f"{-net} bp deleted"
            detail.description = (f"{gene} carries an internal deletion of about {-net} bp "
                                  f"between reference codons {junction[0]} and {junction[2]}; "
                                  f"the protein keeps {intact} of {length} residues "
                                  f"({detail.protein_hgvs}); {detail.resistance}.")
    detail.type_label = TYPE_LABELS[detail.call_type]
    detail.symbol = f"{gene}_{detail.protein_change}"
    if organism == KLEBORATE_ORGANISM and gene in ("ompK35", "ompK36", "mgrB", "pmrB") \
            and intact >= 1:
        last = reference[intact - 1]
        detail.kleborate = f"{KLEBORATE_LOCI[gene]}:p.{AA3.get(last, last)}{intact}fs"
    if elsewhere:
        detail.unresolved = "; ".join(filter(None, [detail.unresolved, (
            "other alignments to this gene elsewhere in the assembly: "
            + ", ".join(f"{p.contig}:{p.gene_lo}-{p.gene_hi} (residues {p.sub_lo}-{p.sub_hi})"
                        for p in elsewhere[:3]))]))
    return detail


def _set_location(detail: MutationDetail, contig: str, lo: int, hi: int, strand: str) -> None:
    detail.contig, detail.start, detail.end, detail.strand = contig, lo, hi, strand
    detail.location = _location(contig, lo, hi, strand)


def _register(first: _Piece, start_top: int, reader: ContigReader):
    """isolate residue index (from the start codon) -> reference residue number."""
    mapping: dict[int, int] = {}
    tops = first.subject_tops()
    for residue, top in tops.items():
        distance = (start_top - top) if reader.minus else (top - start_top)
        if distance % 3 == 0:
            mapping[distance // 3] = residue
    if not mapping:
        return lambda i: i + 1
    last_index = max(mapping)
    last_residue = mapping[last_index]
    first_index = min(mapping)
    first_residue = mapping[first_index]

    def register(i: int) -> int | None:
        if i in mapping:
            return mapping[i]
        if i > last_index:
            return last_residue + (i - last_index)
        if i < first_index:
            return first_residue - (first_index - i)
        return None
    return register


def _continues_in_frame(reader: ContigReader, start_top: int, stop_index: int, register,
                        reference: str, window: int = 12) -> bool:
    """After a stop, does the same frame go on encoding the reference? (nonsense, not fs)"""
    r = register(stop_index)
    if r is None:
        return False
    matches = total = 0
    for k in range(1, window + 1):
        codon = reader.codon(reader.step(start_top, 3 * (stop_index + k)))
        if codon is None or r + k > len(reference):
            break
        total += 1
        matches += CODE.get(codon, "X") == reference[r + k - 1]
    return total >= 4 and matches >= 0.6 * total


def _infer_small_indel(reader: ContigReader, first: _Piece, downstream: _Piece | None,
                       net: int, reference: str, flank: int = 6,
                       max_candidates: int = 3) -> list[tuple[str, str, tuple[int, int]]]:
    """The single-base (or two-base) change that restores the reference protein.

    Without a reference coding sequence the deleted or inserted base is not
    observable. It can be *inferred*: try every way of undoing a net ``net``-bp
    indel in the junction between the two aligned parts and keep those whose
    translation reproduces the reference residues on both sides. Equivalent
    placements inside a homopolymer collapse to one event, reported at the 3'-most
    position as HGVS requires. Returns one (``c.`` notation, contig position text,
    contig span) per distinct change that fits -- several when synonymous codons
    make the change ambiguous -- and an empty list when none (or too many) fit.
    """
    if downstream is None or abs(net) not in (1, 2):
        return []
    tops = first.subject_tops()
    last_ref = first.sub_hi
    anchor_ref = max(first.sub_lo, last_ref - flank)
    if anchor_ref not in tops:
        return []
    anchor_top = tops[anchor_ref]
    down_tops = downstream.subject_tops()
    end_ref = min(downstream.sub_hi, downstream.sub_lo + flank)
    if end_ref not in down_tops:
        return []
    end_top = reader.step(down_tops[end_ref], 2)
    span = abs(end_top - anchor_top) + 1
    isolate = reader.read(anchor_top, span)
    target = reference[anchor_ref - 1:end_ref]
    solutions: dict[str, list[tuple[int, str]]] = {}
    if net < 0:
        for i in range(1, len(isolate)):
            for bases in _all_bases(-net):
                candidate = isolate[:i] + bases + isolate[i:]
                if translate(candidate[:3 * len(target)], False) == target:
                    solutions.setdefault(candidate, []).append((i, bases))
    else:
        for i in range(1, len(isolate) - net):
            candidate = isolate[:i] + isolate[i + net:]
            if translate(candidate[:3 * len(target)], False) == target:
                solutions.setdefault(candidate, []).append((i, isolate[i:i + net]))
    if not solutions or len(solutions) > max_candidates:
        return []
    cds_anchor = 3 * (anchor_ref - 1) + 1
    out = []
    for placements in solutions.values():
        i, bases = max(placements)                  # 3'-most equivalent placement
        first_i = min(p for p, _ in placements)
        if net < 0:
            position = cds_anchor + i               # reference number of the deleted base
            lo_pos = cds_anchor + first_i
            notation = (f"c.{position}del{bases}" if -net == 1
                        else f"c.{position}_{position + 1}del{bases}")
            left = reader.step(anchor_top, first_i - 1)
            right = reader.step(anchor_top, i)
            run = (f", one {bases} of the run c.{lo_pos}-{position}" if lo_pos != position
                   else "")
            where = f"between contig positions {min(left, right)} and {max(left, right)}{run}"
        else:
            position = cds_anchor + i - 1
            notation = f"c.{position}_{position + 1}ins{bases}"
            left = reader.step(anchor_top, first_i)
            right = reader.step(anchor_top, i + net - 1)
            where = f"inserted {bases} at contig {min(left, right)}-{max(left, right)}"
        out.append((notation, where, (min(left, right), max(left, right))))
    return sorted(out)


def _all_bases(n: int):
    if n == 1:
        yield from "ACGT"
    else:
        for a in "ACGT":
            for b in "ACGT":
                yield a + b


def _find_start_upstream(reader: ContigReader, first: _Piece, reference: str) -> int | None:
    """Find the gene's start codon upstream of an alignment that begins late.

    A frameshift in the first few codons leaves the true N-terminus upstream in a
    different frame (the 'N' edge in :func:`.protein.edge_frameshift`). The
    upstream flank is searched, in all three frames, for the reference's first
    residues; the occurrence closest to the alignment wins.
    """
    missing = first.sub_lo - 1
    head = reference[:min(6, missing)]
    if len(head) < 3:
        return None
    anchor = first.subject_tops().get(first.sub_lo)
    if anchor is None:
        return None
    upstream = (reader.length - anchor) if reader.minus else (anchor - 1)
    width = min(3 * missing + 30, upstream)
    if width < 3 * len(head):
        return None
    upstream_top = reader.step(anchor, -width)
    flank = reader.read(upstream_top, width)
    best = None
    for shift in range(3):
        protein = translate(flank[shift:], False)
        index = protein.rfind(head)
        if index >= 0:
            offset = shift + 3 * index
            if best is None or offset > best:
                best = offset
    return None if best is None else reader.step(upstream_top, best)


def _contig_end(hit: Hit, detail: MutationDetail, reader: ContigReader, top: int, length: int,
                kept: int, side: str) -> MutationDetail:
    detail.call_type = "contig_end"
    detail.type_label = TYPE_LABELS["contig_end"]
    detail.symbol = f"{hit.gene} (contig end)"
    end = 1 if reader.minus else reader.length
    edge = "start" if end == 1 else "end"
    detail.location = f"{hit.sequence}:{end} (contig {edge})"
    detail.contig, detail.start, detail.end, detail.strand = (hit.sequence, end, end,
                                                              "-" if reader.minus else "+")
    detail.intact = f"{kept} of {length} residues on this contig ({100.0 * kept / length:.0f}%)"
    detail.description = (f"The {hit.gene} reading frame runs off the {edge} of contig "
                          f"{hit.sequence} ({side} part missing): this is an assembly boundary, "
                          "not evidence that the gene is disrupted. The loss-of-function call "
                          "should not be relied on (AMRFinderPlus reports such genes as "
                          "PARTIAL_CONTIG_END).")
    detail.unresolved = "gene cut by a contig end; the rest of the gene is not in this contig"
    return detail


def _describe_missing_5prime(hit: Hit, detail: MutationDetail, reader: ContigReader,
                             first: _Piece, local: list[_Piece], elsewhere: list[_Piece],
                             reference: str, length: int, organism: str | None) -> MutationDetail:
    """The alignment starts late and no start codon can be found upstream."""
    anchor = first.subject_tops().get(first.sub_lo)
    missing = first.sub_lo - 1
    kept = first.sub_hi - first.sub_lo + 1
    if anchor is not None:
        upstream_bases = (reader.length - anchor) if reader.minus else (anchor - 1)
        if upstream_bases < 3 * missing:
            return _contig_end(hit, detail, reader, anchor, length, kept, "5'")
    last_missing = reference[missing - 1] if missing >= 1 else "X"
    detail.call_type = "disruption_5prime"
    detail.type_label = TYPE_LABELS["disruption_5prime"]
    detail.protein_change = f"M1_{last_missing}{missing}del"
    detail.protein_hgvs = f"p.Met1_{AA3.get(last_missing, 'Xaa')}{missing}del"
    detail.symbol = f"{hit.gene}_{detail.protein_change}"
    if anchor is not None:
        lo, hi = reader.span(anchor, 3)
        _set_location(detail, hit.sequence, lo, hi, first.strand)
        detail.location = (f"{hit.sequence}:{lo}-{hi} ({first.strand}), first aligned codon "
                           f"(reference residue {first.sub_lo})")
    detail.intact = (f"residues {first.sub_lo}-{first.sub_hi} of {length} aligned "
                     f"({100.0 * kept / length:.0f}%); residues 1-{missing} absent at this locus")
    where = ""
    if elsewhere:
        where = (" The missing part aligns elsewhere: "
                 + ", ".join(f"{p.contig}:{p.gene_lo}-{p.gene_hi} (residues {p.sub_lo}-"
                             f"{p.sub_hi})" for p in elsewhere[:3])
                 + ", so the gene is in pieces (an insertion or rearrangement).")
    detail.description = (f"The first {missing} residues of {hit.gene}, start codon included, "
                          "are not upstream of the aligned part: the sequence there is "
                          f"unrelated to {hit.gene} (a 5' deletion, or an insertion that "
                          f"separated the start of the gene).{where} Associated with "
                          f"{detail.resistance}.")
    detail.unresolved = ("the breakpoint is located but the inserted or deleted sequence is "
                         "not identified")
    return detail


# ------------------------------------------------------------ divergent alleles
def describe_divergent(hit: Hit) -> MutationDetail:
    detail = MutationDetail(symbol=hit.gene, call_type="divergent",
                            type_label=TYPE_LABELS["divergent"],
                            resistance=resistance_text(hit.drug_class, hit.subclass))
    detail.gene_location = _location(hit.sequence, min(hit.start, hit.end),
                                     max(hit.start, hit.end), hit.strand)
    detail.description = (f"{hit.gene} differs from the susceptible reference at many positions "
                          f"({hit.identity_pct:.1f}% identity; {hit.note}): a mosaic allele "
                          f"associated with {detail.resistance}, not a single mutation.")
    detail.unresolved = "not a single mutation; no codon-level change applies"
    return detail


# ------------------------------------------------------------------- reads
def _read_support(alt: str, alt_count: int, depth: int, fraction: float) -> str:
    return f"{alt} in {alt_count}/{depth} reads (allele fraction {fraction:.3f})"


def describe_read_site(call, organism: str | None = None) -> MutationDetail:
    """A catalogued DNA site measured from reads mapped to the catalogue's reference locus.

    There is no contig: the position is on the reference the reads were mapped
    to, and the reference base is that reference's own. What the reads carry is
    a count per base, so it is reported as such rather than as a codon.
    """
    locus = call.reference.split("@")[1].replace("_", " ") if "@" in call.reference \
        else call.gene
    rrna = any(tag in (call.gene + call.reference) for tag in ("16S", "23S", "rRNA"))
    promoter = call.position < 0 or "promoter" in call.reference.lower()
    kind = "rRNA" if rrna else ("promoter" if promoter else "nucleotide")
    detail = MutationDetail(symbol=call.symbol, call_type=kind, type_label=TYPE_LABELS[kind],
                            resistance=resistance_text(call.drug_class, call.subclass))
    prefix = "n." if rrna else "c."
    if len(call.alt_base) == 1 and len(call.ref_base) == 1:
        detail.nucleotide_change = f"{prefix}{call.position}{call.ref_base}>{call.alt_base}"
        detail.base_change = f"{call.ref_base}>{call.alt_base}"
    detail.reference = call.ref_base
    detail.nucleotide_basis = f"reference base from the catalogue's reference sequence ({locus})"
    detail.observed = _read_support(call.alt_base, call.alt_count, call.depth,
                                    call.allele_fraction)
    detail.location = (f"{call.reference.split('@')[0]}:{call.position} on the reference "
                       f"{locus} (reads mapped to the reference; no contig)")
    status = {"fixed": "fixed (all copies)", "heteroresistant": "heteroresistant (a minority of "
              "copies or cells)", "absent": "absent", "low-depth": "not assessable (low depth)"}
    detail.description = (f"{locus}: {call.alt_count} of {call.depth} reads carry {call.alt_base} "
                          f"at position {call.position} where the reference has {call.ref_base} "
                          f"({detail.nucleotide_change or call.symbol}), "
                          f"{status.get(call.status, call.status)}: "
                          f"{detail.resistance or 'catalogued resistance mutation'}.")
    return detail


def describe_read_variant(variant, entry, gene: str, reference_name: str,
                          drug_class: str = "", subclass: str = "") -> MutationDetail:
    """A difference between the reads and the closest reference gene they were mapped to.

    The reference here is a real nucleotide sequence (the database allele), so
    the base change is observed, not inferred. The isolate codon is not given:
    reads support each base separately, and a minority allele has no consensus
    codon of its own.
    """
    detail = MutationDetail(resistance=resistance_text(drug_class, subclass))
    protein = variant.protein_change
    if protein:
        ref_aa, alt_aa = variant.ref_aa, variant.alt_aa
        detail.call_type = "nonsense" if alt_aa == "*" else "substitution"
        detail.protein_change = (f"{ref_aa}{variant.codon_position}"
                                 f"{'Ter' if alt_aa == '*' else alt_aa}")
        detail.protein_hgvs = (f"p.{AA3.get(ref_aa, ref_aa)}{variant.codon_position}"
                               f"{AA3.get(alt_aa, alt_aa)}")
    else:
        detail.call_type = "nucleotide"
    detail.type_label = TYPE_LABELS[detail.call_type]
    detail.symbol = entry.symbol if entry is not None else \
        f"{gene}_{detail.protein_change or variant.nucleotide_change}"
    detail.nucleotide_change = f"c.{variant.position}{variant.ref_base}>{variant.alt_base}"
    detail.reference = variant.ref_base
    detail.base_change = f"{variant.ref_base}>{variant.alt_base}" + (
        f" at codon position {(variant.position - 1) % 3 + 1}" if protein else "")
    detail.nucleotide_basis = f"reference base from the reference gene {reference_name}"
    detail.observed = _read_support(variant.alt_base, variant.alt_count, variant.depth,
                                    variant.allele_fraction)
    detail.location = (f"{reference_name}:{variant.position} on the closest reference gene "
                       "(reads mapped to the reference; no contig)")
    detail.description = (f"{gene}: {variant.alt_count} of {variant.depth} reads carry "
                          f"{variant.alt_base} at position {variant.position} of {reference_name} "
                          f"where it has {variant.ref_base}"
                          + (f" ({detail.protein_hgvs})" if protein else "")
                          + (f": catalogued as {entry.symbol}, {detail.resistance}."
                             if entry is not None else ": not in the mutation catalogue."))
    if entry is None:
        detail.unresolved = "not a catalogued resistance mutation"
    return detail


# --------------------------------------------------------------------- driver
def describe(hit: Hit, contigs: Mapping[str, str], organism: str | None) -> MutationDetail | None:
    evidence = hit.evidence or {}
    kind = evidence.get("kind")
    if kind == "protein":
        return describe_protein(hit, evidence, contigs.get(hit.sequence), organism)
    if kind == "dna":
        return describe_dna(hit, evidence, contigs.get(hit.sequence))
    if kind == "disruption":
        return describe_disruption(hit, evidence, contigs, organism)
    if hit.method == "SUSCEPTIBLEX":
        return describe_divergent(hit)
    return None


def describe_mutations(hits: Sequence[Hit], assembly: Path | str | None,
                       organism: str | None) -> int:
    """Attach a :class:`MutationDetail` to every mutation call among *hits*.

    Only the contigs carrying a call are read from the assembly. A failure to
    describe one call is logged and recorded on that call; it never removes or
    alters the call itself. Returns the number of calls described.
    """
    wanted = [h for h in hits if h.evidence or h.method == "SUSCEPTIBLEX"]
    if not wanted:
        return 0
    contigs: dict[str, str] = {}
    if assembly is not None:
        names = {h.sequence for h in wanted}
        for h in wanted:
            for _hsp, piece in (h.evidence or {}).get("pieces") or []:
                if piece is not None:
                    names.add(piece.contig)
        try:
            contigs = read_fasta_subset(assembly, names)
        except OSError as exc:
            LOG.warning("could not read %s to describe mutations: %s", assembly, exc)
    done = 0
    for hit in wanted:
        try:
            hit.mutation = describe(hit, contigs, organism)
        except Exception as exc:  # noqa: BLE001 - a description must never sink a call
            LOG.warning("could not describe %s on %s: %s", hit.gene, hit.sequence, exc)
            LOG.debug("traceback", exc_info=True)
            hit.mutation = MutationDetail(symbol=hit.gene, call_type="unresolved",
                                          type_label="Not described",
                                          unresolved=f"description failed: {exc}")
        hit.evidence = None
        # A gene cut by a contig end is kept, and the row is made to say so.
        disclose_contig_end(hit)
        done += hit.mutation is not None
    return done
