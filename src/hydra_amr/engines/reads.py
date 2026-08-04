"""Detection directly from raw reads, including heteroresistance.

Assembly consensus hides minority alleles. When a resistance mutation sits in
only some copies of a multi-copy locus — the classic case being the 23S rRNA
mutations that confer linezolid resistance across four to six rRNA operons — the
assembled contig carries the wild-type base and every assembly-only caller
reports a susceptible genome. Hydra therefore aligns the reads themselves to the
reference loci and reports the *allele fraction* at each catalogued position, so
a mutation present in one operon out of five is visible as a ~20% minority
allele rather than being lost.

The same alignment pass yields depth-and-breadth based acquired-gene calls, so
a FASTQ-only sample still gets a full resistance profile.
"""

from __future__ import annotations

import math
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from ..config import Config
from ..db.manager import DatabaseStore
from ..records import Hit
from ..utils import LOG, HydraError, have, natural_key, require, run, translate
from .mutations import MutationCatalog, reference_index

#: Per-base sequencing error rate assumed by the minority-allele significance test.
BACKGROUND_ERROR_RATE = 0.002

#: pysam's "all" stepper drops unmapped, secondary, QC-failed and duplicate
#: reads while keeping the ones whose mate maps elsewhere. The "samtools"
#: stepper would additionally discard those orphans, which on a gene-sized
#: reference is most of the data; "nofilter" would keep everything, and
#: silently ignores any flag_filter passed alongside it.
PILEUP_STEPPER = "all"
#: Multi-mapping reads carry a low mapping quality; counting them at full
#: weight invents minority alleles at every position a paralogue overlaps.
MIN_MAPPING_QUALITY = 20
#: Deep enough for an rRNA operon collapsed onto one reference, without letting
#: one column allocate unbounded memory.
MAX_PILEUP_DEPTH = 200_000
#: Allele fractions are only claimed once a site is this deep. The gene-detection
#: floor is much lower, because "is this gene here" needs far less evidence than
#: "what fraction of the population carries this base".
MIN_VARIANT_DEPTH = 20
#: Share of reads carrying an indel before a position is treated as indel-bearing.
INDEL_EVIDENCE_FRACTION = 0.15
#: A consensus base supported by only this share of the reads is not a clean
#: majority; several such positions in one locus mean a mixed sample.
MIXED_ALLELE_RANGE = (0.20, 0.80)

#: rRNA operon copy numbers, used to translate an allele fraction into an
#: estimated number of mutated operons.
RRNA_OPERON_COUNTS: dict[str, int] = {
    "Staphylococcus_aureus": 5,
    "Enterococcus_faecalis": 4,
    "Enterococcus_faecium": 6,
    "Streptococcus_pneumoniae": 4,
    "Escherichia": 7,
    "Klebsiella_pneumoniae": 8,
    "Klebsiella_oxytoca": 8,
    "Salmonella": 7,
    "Campylobacter": 3,
    "Clostridioides_difficile": 11,
    "Acinetobacter_baumannii": 6,
    "Neisseria_gonorrhoeae": 4,
    "Streptococcus_pyogenes": 6,
    "Streptococcus_agalactiae": 7,
    "Pseudomonas_aeruginosa": 4,
}


@dataclass
class ReadSet:
    """One sample's reads."""

    sample: str
    r1: Path
    r2: Path | None = None
    interleaved: bool = False
    single: bool = False

    @property
    def files(self) -> list[Path]:
        return [self.r1] + ([self.r2] if self.r2 else [])


@dataclass
class Variant:
    """One difference between the reads and the reference they were mapped to."""

    reference: str
    position: int          # 1-based on the reference
    ref_base: str
    alt_base: str
    depth: int
    alt_count: int
    allele_fraction: float
    #: Amino-acid consequence when the reference is an in-frame coding sequence.
    codon_position: int | None = None
    ref_aa: str = ""
    alt_aa: str = ""
    #: The run's --fixed-allele-fraction, so "fixed" means the same here as it
    #: does for the catalogued-site calls.
    fixed_threshold: float = 0.9

    @property
    def is_fixed(self) -> bool:
        return self.allele_fraction >= self.fixed_threshold

    @property
    def nucleotide_change(self) -> str:
        return f"{self.ref_base}{self.position}{self.alt_base}"

    @property
    def protein_change(self) -> str:
        if not self.ref_aa or self.codon_position is None:
            return ""
        return f"{self.ref_aa}{self.codon_position}{self.alt_aa}"

    @property
    def is_synonymous(self) -> bool:
        return bool(self.ref_aa) and self.ref_aa == self.alt_aa


@dataclass
class LocusConsensus:
    """The sequence the reads support at one reference locus."""

    reference: str
    sequence: str          # '-' where depth was insufficient
    depth: float
    breadth: float
    variants: list[Variant] = field(default_factory=list)
    #: Positions where a substantial share of the reads carry an insertion or
    #: deletion. A consensus is always the reference's length, so an indel
    #: cannot be represented in it - only flagged.
    indel_sites: list[int] = field(default_factory=list)
    #: Positions carrying a well-supported intermediate allele: the signature of
    #: a mixed sample when they pile up across a housekeeping locus.
    mixed_sites: list[int] = field(default_factory=list)

    @property
    def called_bases(self) -> int:
        return sum(1 for base in self.sequence if base != "-")

    @property
    def has_indel_evidence(self) -> bool:
        return bool(self.indel_sites)


@dataclass
class SiteCall:
    """Allele fractions observed at one catalogued position."""

    reference: str
    position: int
    symbol: str
    gene: str
    ref_base: str
    alt_base: str
    depth: int
    ref_count: int
    alt_count: int
    allele_fraction: float
    p_value: float
    drug_class: str = ""
    subclass: str = ""
    name: str = ""
    status: str = "absent"   # fixed | heteroresistant | absent | low-depth
    estimated_operons: float | None = None
    counts: dict[str, int] = field(default_factory=dict)


#: A pileup token is a base, optionally followed by an indel: ``A``, ``A+2GT``,
#: ``A-1N``, or ``*`` where a deletion covers the position itself.
_PILEUP_TOKEN = re.compile(r"^(?P<base>[ACGTN*])(?:(?P<op>[+-])(?P<length>\d+)(?P<seq>[ACGTN]*))?$")


def _base_counts(tokens: list[str]) -> dict[str, int]:
    """Count the aligned base of each pileup token, ignoring any indel suffix."""
    counts: dict[str, int] = {}
    for token in tokens:
        match = _PILEUP_TOKEN.match(token)
        if match is None:
            continue
        base = match.group("base")
        if base in "ACGT":
            counts[base] = counts.get(base, 0) + 1
    return counts


def _count_alt(tokens: list[str], entry) -> int | None:
    """Reads supporting *entry*'s alternate allele, or None if it is not callable."""
    change = entry.change_type
    if change in ("substitution", "nonsense"):
        return _base_counts(tokens).get(entry.alt.upper()[:1], 0)
    if change == "deletion":
        # samtools marks a base removed by a deletion with '*'.
        return sum(1 for token in tokens if token.startswith("*"))
    if change == "insertion":
        inserted = entry.inserted_bases
        count = 0
        for token in tokens:
            match = _PILEUP_TOKEN.match(token)
            if match is None or match.group("op") != "+":
                continue
            if not inserted or match.group("seq") == inserted:
                count += 1
        return count
    return None


def _looks_like_cds(sequence: str) -> bool:
    """True when *sequence* can be read as a single coding frame.

    Several databases hold promoter regions, rRNA and intergenic replicons whose
    length happens to be a multiple of three. Translating those would invent
    amino-acid changes and, worse, discard real nucleotide variants as
    "synonymous" in a frame that does not exist.
    """
    if len(sequence) < 60 or len(sequence) % 3:
        return False
    protein = translate(sequence, table_start_is_met=False)
    if not protein:
        return False
    body = protein[:-1] if protein.endswith("*") else protein
    return "*" not in body


def _indel_count(tokens: list[str]) -> int:
    """Reads whose alignment inserts or deletes sequence at this position."""
    total = 0
    for token in tokens:
        if token.startswith("*"):
            total += 1
            continue
        match = _PILEUP_TOKEN.match(token)
        if match is not None and match.group("op"):
            total += 1
    return total


def _annotate_codons(consensus: LocusConsensus, reference: str) -> None:
    """Attach the amino-acid consequence to each variant of a coding sequence."""
    for variant in consensus.variants:
        index = variant.position - 1
        codon_start = (index // 3) * 3
        codon = reference[codon_start:codon_start + 3]
        if len(codon) < 3:
            continue
        mutated = list(codon)
        mutated[index - codon_start] = variant.alt_base
        variant.codon_position = codon_start // 3 + 1
        variant.ref_aa = translate(codon, table_start_is_met=False)
        variant.alt_aa = translate("".join(mutated), table_start_is_met=False)


def binomial_upper_tail(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p).

    Used to ask whether the minority allele exceeds what sequencing error alone
    would produce. The upper tail is summed directly rather than complementing
    the lower tail: for a convincing call the tail is far below 1e-16, and
    ``1 - cumulative`` would collapse every such p-value onto the same
    floating-point noise floor.
    """
    if n <= 0 or k <= 0:
        return 1.0
    if k > n:
        return 0.0
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    log_p = math.log(p)
    log_q = math.log1p(-p)
    log_n_fact = math.lgamma(n + 1)
    total = 0.0
    for i in range(k, n + 1):
        log_term = (log_n_fact - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                    + i * log_p + (n - i) * log_q)
        term = math.exp(log_term)
        total += term
        # Terms decay geometrically past the mode; stop once they cannot matter.
        if i > n * p and term < total * 1e-12:
            break
    return max(0.0, min(1.0, total))


class ReadMapper:
    """Aligns reads to reference FASTAs and summarises coverage and alleles."""

    def __init__(self, store: DatabaseStore, config: Config):
        self.store = store
        self.config = config
        self._pysam = None

    @property
    def pysam(self):
        if self._pysam is None:
            try:
                import pysam  # noqa: PLC0415 - optional dependency, imported on demand
            except ImportError as exc:
                raise HydraError(
                    "read mode needs pysam.\n"
                    "  conda install -c bioconda -c conda-forge pysam"
                ) from exc
            self._pysam = pysam
        return self._pysam

    @staticmethod
    def check_tools() -> None:
        require("minimap2", "read alignment")
        require("samtools", "BAM sorting and indexing")

    def align(self, reads: ReadSet, reference: Path, out_bam: Path, threads: int = 1,
              preset: str = "sr") -> Path:
        """Align *reads* to *reference*, returning the sorted, indexed BAM path.

        The alignment is streamed straight into ``samtools sort``: writing a SAM
        first would cost several gigabytes of temporary space per sample, most of
        it the unmapped reads, which ``--sam-hit-only`` discards anyway. Only a
        few genes' worth of the run is ever of interest here.
        """
        self.check_tools()
        out_bam = Path(out_bam)
        out_bam.parent.mkdir(parents=True, exist_ok=True)
        map_threads = max(1, threads - 1) if threads > 1 else 1
        sort_threads = max(1, threads // 4)
        cmd = ["minimap2", "-ax", preset, "-t", str(map_threads), "--secondary=no",
               "--sam-hit-only", "-N", "5", str(reference), str(reads.r1)]
        if reads.r2:
            cmd.append(str(reads.r2))
        sort_cmd = ["samtools", "sort", "-@", str(sort_threads), "-o", str(out_bam), "-"]
        LOG.debug("exec: %s | %s", " ".join(cmd), " ".join(sort_cmd))
        # Both stderr streams go to files rather than pipes. minimap2 emits a
        # warning per malformed read; on a pipe those would fill the 64 KiB
        # buffer, block minimap2, starve samtools of input and hang the run.
        map_log = out_bam.with_suffix(".minimap2.log")
        sort_log = out_bam.with_suffix(".samtools.log")
        with open(map_log, "wb") as map_err, open(sort_log, "wb") as sort_err:
            mapper = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=map_err)
            try:
                sorter = subprocess.Popen(sort_cmd, stdin=mapper.stdout, stderr=sort_err)
            except FileNotFoundError as exc:
                mapper.kill()
                mapper.wait()
                raise HydraError("samtools not found on PATH") from exc
            # Close our copy so minimap2 sees EPIPE if samtools exits early.
            mapper.stdout.close()
            sorter.wait()
            mapper.wait()

        def _tail(path: Path) -> str:
            try:
                return path.read_text(errors="replace")[-2000:].strip()
            except OSError:
                return ""

        # Report the downstream failure first: when samtools dies, minimap2's
        # non-zero exit is only the SIGPIPE that followed.
        if sorter.returncode != 0:
            raise HydraError(f"samtools sort failed (exit {sorter.returncode}):\n"
                             f"{_tail(sort_log)}")
        if mapper.returncode != 0:
            raise HydraError(f"minimap2 failed (exit {mapper.returncode}):\n{_tail(map_log)}")
        run(["samtools", "index", str(out_bam)])
        for log in (map_log, sort_log):
            log.unlink(missing_ok=True)
        return out_bam

    # ------------------------------------------------------------------ genes
    def call_genes(self, sample: str, bam: Path, db_name: str, handle,
                   ) -> list[Hit]:
        """Report reference genes covered deeply and broadly enough by the reads."""
        pysam = self.pysam
        thresholds = self.config.thresholds
        hits: list[Hit] = []
        with pysam.AlignmentFile(str(bam), "rb") as align:
            references = {}
            for stats in align.get_index_statistics():
                if stats.mapped > 0:
                    references[stats.contig] = stats.mapped
            for reference, n_reads in references.items():
                length = align.get_reference_length(reference)
                if not length:
                    continue
                coverage = align.count_coverage(reference, quality_threshold=0,
                                                read_callback="nofilter")
                depths = [coverage[0][i] + coverage[1][i] + coverage[2][i] + coverage[3][i]
                          for i in range(length)]
                covered = sum(1 for d in depths if d >= 1)
                breadth = 100.0 * covered / length
                if breadth < thresholds.min_gene_breadth:
                    continue
                deep = [d for d in depths if d >= 1]
                mean_depth = sum(deep) / len(deep) if deep else 0.0
                # The depth cut-off is applied after family collapsing: reads are
                # split between near-identical alleles, so any single allele can
                # sit below it while the locus as a whole is well covered.
                identity = self._reference_identity(align, reference, length)
                if identity < thresholds.min_identity:
                    continue
                meta = handle.meta_for(reference)
                hit = Hit(
                    sample=sample, database=db_name,
                    gene=meta.get("gene", reference), accession=meta.get("accession", ""),
                    product=meta.get("product", ""),
                    element_type=meta.get("element_type", handle.element_type),
                    element_subtype=meta.get("element_subtype", ""),
                    drug_class=meta.get("class", ""), subclass=meta.get("subclass", ""),
                    sequence=reference, start=1, end=length, strand="+",
                    coverage=f"1-{length}/{length}", coverage_pct=breadth,
                    identity_pct=identity, gaps=0, depth=mean_depth,
                    method="READS", resolution="COMPLETE" if breadth >= 99.0 else "PARTIAL",
                    note=f"{n_reads} reads",
                )
                hits.append((hit, meta.get("fam_id", "")))
        collapsed = self._collapse_families(hits)
        return [h for h in collapsed if (h.depth or 0.0) >= thresholds.min_depth]

    @staticmethod
    def _reference_identity(align, reference: str, length: int) -> float:
        """Weighted identity of the reads aligned to one reference."""
        matches = 0
        total = 0
        for read in align.fetch(reference):
            if read.is_unmapped:
                continue
            nm = read.get_tag("NM") if read.has_tag("NM") else 0
            aligned = read.query_alignment_length or 0
            if aligned <= 0:
                continue
            matches += max(0, aligned - nm)
            total += aligned
        return (100.0 * matches / total) if total else 0.0

    @staticmethod
    def _collapse_families(hits: list[tuple[Hit, str]]) -> list[Hit]:
        """Reduce each gene family to one call, recovering its split read depth.

        A gene database holds hundreds of near-identical alleles - 14 variants of
        ``aac(6')-Ib``, dozens of ``blaTEM`` - and a read matching the locus maps
        to exactly one of them. Reporting every allele that picked up reads turns
        one gene into fourteen, and dividing the depth between them can push all
        of them below the depth threshold so the locus disappears entirely.
        Grouping by the curated gene family and summing the depth reports one
        gene at its true coverage.
        """
        groups: dict[str, list[tuple[Hit, str]]] = {}
        for hit, fam_id in hits:
            groups.setdefault(fam_id or hit.gene, []).append((hit, fam_id))
        out: list[Hit] = []
        for key, members in groups.items():
            best = max(members, key=lambda m: (m[0].coverage_pct, m[0].identity_pct,
                                               m[0].depth or 0.0))[0]
            if len(members) > 1:
                total_depth = sum(m[0].depth or 0.0 for m in members)
                names = sorted({m[0].gene for m in members}, key=natural_key)
                best.depth = total_depth
                best.note = (f"{best.note}; {len(members)} alleles of family {key} shared "
                             f"the reads ({', '.join(names[:4])}"
                             f"{', ...' if len(names) > 4 else ''})")
            out.append(best)
        return sorted(out, key=lambda h: h.gene)

    # -------------------------------------------------------- consensus/variants
    def consensus(self, bam: Path, references: dict[str, str], min_depth: int = 5,
                  min_base_quality: int = 13, min_allele_fraction: float = 0.1,
                  wanted: set[str] | None = None, min_alt_reads: int = 3,
                  max_p_value: float = 1e-3, min_variant_depth: int = MIN_VARIANT_DEPTH,
                  fixed_allele_fraction: float = 0.9) -> dict[str, LocusConsensus]:
        """Read the sequence the reads support at each reference, and its variants.

        This is what lets a FASTQ-only sample be typed and have its resistance
        genes inspected: reads are mapped to the closest reference Hydra holds,
        the majority base is taken at every position, and every difference from
        that reference is reported with the fraction of reads supporting it.
        """
        pysam = self.pysam
        out: dict[str, LocusConsensus] = {}
        with pysam.AlignmentFile(str(bam), "rb") as align:
            available = set(align.references)
            targets = [r for r in (wanted or available) if r in available and r in references]
            for reference in targets:
                sequence = references[reference].upper()
                length = len(sequence)
                bases = ["-"] * length
                depths = [0] * length
                variants: list[Variant] = []
                indel_sites: list[int] = []
                mixed_sites: list[int] = []
                for column in align.pileup(reference, min_base_quality=min_base_quality,
                                           min_mapping_quality=MIN_MAPPING_QUALITY,
                                           max_depth=MAX_PILEUP_DEPTH, stepper=PILEUP_STEPPER,
                                           truncate=True):
                    index = column.reference_pos
                    if index >= length:
                        continue
                    tokens = [t.upper() for t in column.get_query_sequences(
                        mark_matches=False, mark_ends=False, add_indels=True) if t]
                    counts = _base_counts(tokens)
                    indels = _indel_count(tokens)
                    # Deleted positions count toward the depth: without them a
                    # deletion looks like missing coverage rather than evidence.
                    depth = sum(counts.values()) + indels
                    depths[index] = depth
                    if depth < min_depth or not counts:
                        continue
                    if indels >= max(min_alt_reads, INDEL_EVIDENCE_FRACTION * depth):
                        # A consensus is fixed to the reference's length, so an
                        # indel cannot be written into it. Record the position so
                        # callers know this locus is not simply a substitution
                        # variant of the reference.
                        indel_sites.append(index + 1)
                    best, best_count = max(counts.items(), key=lambda kv: (kv[1], kv[0]))
                    bases[index] = best
                    reference_base = sequence[index]
                    if (depth >= MIN_VARIANT_DEPTH
                            and MIXED_ALLELE_RANGE[0] <= best_count / depth <= MIXED_ALLELE_RANGE[1]):
                        mixed_sites.append(index + 1)
                    if depth < min_variant_depth:
                        continue
                    for base, count in counts.items():
                        fraction = count / depth
                        if base == reference_base or fraction < min_allele_fraction:
                            continue
                        # A single mismatching read at 25x is 4% of the depth and
                        # is almost always a sequencing error, so a minority
                        # allele has to clear both a read count and the error
                        # model before it is called.
                        if count < min_alt_reads and fraction < fixed_allele_fraction:
                            continue
                        p_value = binomial_upper_tail(count, depth, BACKGROUND_ERROR_RATE)
                        if p_value > max_p_value and fraction < fixed_allele_fraction:
                            continue
                        variants.append(Variant(
                            reference=reference, position=index + 1,
                            ref_base=reference_base, alt_base=base, depth=depth,
                            alt_count=count, allele_fraction=fraction,
                            fixed_threshold=fixed_allele_fraction))
                covered = sum(1 for d in depths if d >= min_depth)
                deep = [d for d in depths if d > 0]
                consensus = LocusConsensus(
                    reference=reference, sequence="".join(bases),
                    indel_sites=indel_sites, mixed_sites=mixed_sites,
                    depth=(sum(deep) / len(deep)) if deep else 0.0,
                    breadth=100.0 * covered / length if length else 0.0,
                    variants=variants,
                )
                if _looks_like_cds(sequence):
                    _annotate_codons(consensus, sequence)
                out[reference] = consensus
        return out

    def variant_hits(self, sample: str, database: str, consensus: LocusConsensus, meta: dict,
                     catalog: MutationCatalog | None = None,
                     synonymous: bool = False) -> list[Hit]:
        """Turn read-derived variants in one gene into reportable hits.

        Every hit records the reference it was measured against, because that is
        what the call means: this is variation relative to the closest sequence
        in the database, not relative to the isolate's own assembled gene.
        """
        gene = meta.get("gene", consensus.reference)
        reference_name = meta.get("accession") or consensus.reference
        catalogued: dict[str, object] = {}
        if catalog is not None:
            # Keyed on the gene symbol, not the accession: the catalogue records
            # protein accessions while every nucleotide database carries a
            # nucleotide one, so an accession lookup can never match.
            wanted = gene.lower()
            for entries in catalog.protein.values():
                for entry in entries:
                    if entry.gene.lower() == wanted:
                        catalogued[entry.symbol.rsplit("_", 1)[-1].upper()] = entry

        def make(variant: Variant, entry, note: str, method: str) -> Hit:
            return Hit(
                sample=sample, database=database, gene=gene,
                accession=meta.get("accession", ""), product=meta.get("product", ""),
                element_type=meta.get("element_type", "AMR"), element_subtype="POINT",
                drug_class=getattr(entry, "drug_class", "") or meta.get("class", ""),
                subclass=getattr(entry, "subclass", "") or meta.get("subclass", ""),
                sequence=consensus.reference, start=variant.position, end=variant.position,
                strand="+", coverage=f"{variant.position}/{len(consensus.sequence)}",
                coverage_pct=consensus.breadth,
                identity_pct=reference_identity,
                depth=float(variant.depth), allele_fraction=variant.allele_fraction,
                method=method, resolution="POINT", note=note,
            )

        called = consensus.called_bases or 1
        reference_identity = 100.0 * (1 - len(consensus.variants) / called)
        hits: list[Hit] = []
        fixed: list[Variant] = []
        for variant in consensus.variants:
            if variant.is_synonymous and not synonymous:
                continue
            change = variant.protein_change or variant.nucleotide_change
            entry = catalogued.get(change.upper())
            if entry is None and variant.is_fixed:
                # A fixed difference just means the database has no exact copy of
                # this allele. Interesting in aggregate, not one row per base.
                fixed.append(variant)
                continue
            status = "FIXED" if variant.is_fixed else "HETERORESISTANT"
            note = (f"{gene} {change} vs closest reference {reference_name}; "
                    f"AF={variant.allele_fraction:.3f}; "
                    f"{variant.alt_count}/{variant.depth} reads; {status}")
            note += (f"; catalogued as {entry.symbol}" if entry is not None
                     else "; not in the mutation catalogue")
            hits.append(make(variant, entry, note,
                             "POINTR" if entry is not None else "VARIANTR"))

        if fixed:
            changes = ", ".join(v.protein_change or v.nucleotide_change for v in fixed[:6])
            first = fixed[0]
            note = (f"{gene} differs from its closest reference {reference_name} at "
                    f"{len(fixed)} protein-changing position(s): {changes}"
                    f"{', ...' if len(fixed) > 6 else ''}; this is a variant allele, "
                    f"not a minority population")
            summary = make(first, None, note, "ALLELER")
            summary.resolution = "ALLELE"
            summary.start, summary.end = 1, len(consensus.sequence)
            summary.allele_fraction = None
            hits.append(summary)
        return hits

    # -------------------------------------------------------------- mutations
    def call_sites(self, bam: Path, catalog: MutationCatalog, reference_lengths: dict[str, int],
                   organism: str | None = None, min_base_quality: int = 13) -> list[SiteCall]:
        """Measure the allele fraction at every catalogued DNA mutation position."""
        pysam = self.pysam
        thresholds = self.config.thresholds
        operons = RRNA_OPERON_COUNTS.get(organism or "", 0)
        calls: list[SiteCall] = []
        with pysam.AlignmentFile(str(bam), "rb") as align:
            available = set(align.references)
            for reference, entries in catalog.dna.items():
                if reference not in available:
                    continue
                length = reference_lengths.get(reference) or align.get_reference_length(reference)
                wanted: dict[int, list] = {}
                for entry in entries:
                    index = reference_index(entry.position, length)
                    if 0 <= index < length:
                        wanted.setdefault(index, []).append(entry)
                if not wanted:
                    continue
                lo, hi = min(wanted), max(wanted)
                observed: dict[int, list[str]] = {}
                for column in align.pileup(reference, start=lo, stop=hi + 1, truncate=True,
                                           min_base_quality=min_base_quality,
                                           min_mapping_quality=MIN_MAPPING_QUALITY,
                                           max_depth=MAX_PILEUP_DEPTH, stepper=PILEUP_STEPPER):
                    if column.reference_pos not in wanted:
                        continue
                    # Indel annotations are kept: an entry such as ampC_G-15GG is
                    # an inserted base, and reading only the first character would
                    # make its alternate allele identical to the reference.
                    observed[column.reference_pos] = [
                        token.upper() for token in column.get_query_sequences(
                            mark_matches=False, mark_ends=False, add_indels=True) if token
                    ]
                for index, entry_list in wanted.items():
                    tokens = observed.get(index, [])
                    counts = _base_counts(tokens)
                    depth = len(tokens)
                    for entry in entry_list:
                        alt = entry.alt.upper()
                        ref = entry.ref.upper()
                        ref_count = counts.get(ref[:1], 0)
                        alt_count = _count_alt(tokens, entry)
                        if alt_count is None:
                            LOG.debug("skipping %s: %s change is not callable from reads",
                                      entry.symbol, entry.change_type)
                            continue
                        fraction = (alt_count / depth) if depth else 0.0
                        if depth < thresholds.min_depth:
                            status = "low-depth"
                        elif fraction >= thresholds.fixed_allele_fraction:
                            status = "fixed"
                        elif (fraction >= thresholds.min_allele_fraction
                              and alt_count >= thresholds.min_allele_reads):
                            status = "heteroresistant"
                        else:
                            status = "absent"
                        p_value = binomial_upper_tail(alt_count, depth, BACKGROUND_ERROR_RATE) \
                            if depth else 1.0
                        if operons and status in ("fixed", "heteroresistant"):
                            estimated = fraction * operons
                        else:
                            estimated = None
                        calls.append(SiteCall(
                            reference=reference, position=entry.position, symbol=entry.symbol,
                            gene=entry.gene, ref_base=ref[:1], alt_base=alt, depth=depth,
                            ref_count=ref_count, alt_count=alt_count, allele_fraction=fraction,
                            p_value=p_value, drug_class=entry.drug_class,
                            subclass=entry.subclass, name=entry.name, status=status,
                            estimated_operons=estimated, counts=dict(counts),
                        ))
        return calls

    def site_hits(self, sample: str, calls: list[SiteCall], report_absent: bool = False,
                  organism: str | None = None) -> list[Hit]:
        """Turn site calls into reportable hits."""
        operons = RRNA_OPERON_COUNTS.get(organism or "", 0)
        hits: list[Hit] = []
        for call in calls:
            if call.status in ("absent", "low-depth") and not report_absent:
                continue
            note = f"{call.symbol}; AF={call.allele_fraction:.3f}; {call.alt_count}/{call.depth} reads"
            if call.status == "heteroresistant":
                note += "; HETERORESISTANT"
                if call.estimated_operons and operons:
                    note += f"; ~{call.estimated_operons:.1f}/{operons} operons"
            elif call.status == "fixed":
                note += "; FIXED"
            if call.p_value < 0.05:
                note += f"; p={call.p_value:.2g}"
            hits.append(Hit(
                sample=sample, database="amrfinderplus", gene=call.gene,
                accession=call.reference.split("@")[0], product=call.name.replace("_", " "),
                element_type="AMR", element_subtype="POINT",
                drug_class=call.drug_class, subclass=call.subclass,
                sequence=call.reference, start=call.position, end=call.position, strand="+",
                coverage=f"{call.position}", coverage_pct=100.0,
                identity_pct=100.0 * (1.0 - call.allele_fraction) if call.status == "absent"
                else 100.0,
                depth=float(call.depth), allele_fraction=call.allele_fraction,
                method="POINTR", resolution="POINT", note=note,
            ))
        return hits


def pair_reads(paths: list[Path]) -> list[ReadSet]:
    """Group FASTQ paths into paired samples using common R1/R2 conventions."""
    from ..utils import sample_name_from_path

    buckets: dict[str, dict[str, Path]] = {}
    singles: list[Path] = []
    markers = (("_R1_001", "_R2_001"), ("_R1", "_R2"), ("_1", "_2"), (".R1", ".R2"), (".1", ".2"))
    for path in paths:
        stem = sample_name_from_path(path)
        matched = False
        for first, second in markers:
            if stem.endswith(first):
                buckets.setdefault(stem[: -len(first)], {})["r1"] = path
                matched = True
                break
            if stem.endswith(second):
                buckets.setdefault(stem[: -len(second)], {})["r2"] = path
                matched = True
                break
        if not matched:
            singles.append(path)
    out: list[ReadSet] = []
    for sample, pair in sorted(buckets.items()):
        r1 = pair.get("r1")
        r2 = pair.get("r2")
        if r1 and r2:
            out.append(ReadSet(sample=sample, r1=r1, r2=r2))
        elif r1 or r2:
            only = r1 or r2
            out.append(ReadSet(sample=sample, r1=only, single=True))
    for path in singles:
        out.append(ReadSet(sample=sample_name_from_path(path), r1=path, single=True))
    return out


def assemble(reads: ReadSet, outdir: Path, threads: int = 1) -> Path:
    """Assemble reads when an assembler is available, returning the contigs path."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if have("skesa"):
        contigs = outdir / f"{reads.sample}.contigs.fasta"
        cmd = ["skesa", "--cores", str(threads), "--contigs_out", str(contigs)]
        if reads.r2:
            cmd += ["--reads", f"{reads.r1},{reads.r2}"]
        else:
            cmd += ["--reads", str(reads.r1)]
        run(cmd)
        return contigs
    if have("spades.py"):
        spades_dir = outdir / f"{reads.sample}.spades"
        cmd = ["spades.py", "--isolate", "-t", str(threads), "-o", str(spades_dir)]
        if reads.r2:
            cmd += ["-1", str(reads.r1), "-2", str(reads.r2)]
        else:
            cmd += ["-s", str(reads.r1)]
        run(cmd)
        contigs = spades_dir / "contigs.fasta"
        if not contigs.exists():
            raise HydraError(f"SPAdes finished but produced no contigs for {reads.sample}")
        final = outdir / f"{reads.sample}.contigs.fasta"
        shutil.copy2(contigs, final)
        return final
    raise HydraError(
        "--assemble needs an assembler on PATH (skesa or spades.py).\n"
        "  conda install -c bioconda skesa"
    )
