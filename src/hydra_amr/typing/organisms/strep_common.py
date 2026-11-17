"""Machinery shared by the two streptococcal typing modules.

Both *Streptococcus pyogenes* and *Streptococcus agalactiae* are typed the same
way underneath: a short reference-element BLAST for the antigen locus (*emm*,
*cps*), and a protein-level screen of a curated gene panel for the virulence,
regulatory and surface-protein repertoire. The two organism modules differ in
which panel they screen and how the answers are named, so everything that is
not organism-specific lives here.

Three decisions in this file are worth stating up front, because the calls rest
on them.

**Reference proteins come from RefSeq annotation of named reference strains,
not from a vendored FASTA.** :func:`fetch_panel` downloads the feature table and
the protein FASTA of each pinned assembly and picks the protein whose RefSeq
gene symbol (or, failing that, product name) matches the panel entry. The
accession it picked is carried into every call, so a result can always be traced
to a sequence someone else can download. This also means a gene the current
RefSeq annotation does not name is *absent from the panel*, which is reported as
"not determined", never as a negative.

**Presence is measured at the protein level.** A translated search (blastx:
assembly as query, reference proteins as subject) tolerates the synonymous
divergence that separates streptococcal lineages, and it gives the residue-level
alignment that a lesion call needs. Nucleotide identity would make a real gene
look absent at the edges of a species.

**A gene that is not intact is not simply "present".** :func:`describe_lesion`
reads the alignment and says what is wrong with the gene - premature stop,
truncation, internal deletion, split across contigs - because for *covR*/*covS*
and *rocA* the lesion *is* the result, and for everything else a fragment at a
contig end must not be reported as an intact gene.
"""

from __future__ import annotations

import gzip
import json
import re
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

from ...engines.blast import OUTFMT_FIELDS, blast, merge_hsps
from ...seqio import read_fasta, write_fasta
from ...utils import LOG, HydraError, require, run
from .base import OrganismDb, Reference

__all__ = [
    "GeneSpec", "GeneHit", "OrganismDb", "Reference", "PanelSource",
    "download", "fetch_panel", "make_blast_db", "screen_panel",
    "describe_lesion", "span_string", "pct", "contig_lengths",
    "NCBI_GENOMES", "USER_AGENT",
]

#: Upstream is told who is calling, as :mod:`hydra_amr.db.fetch` already does.
USER_AGENT = "hydra-amr (+https://github.com/iowa69/hydra)"
TIMEOUT = 180
NCBI_GENOMES = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

#: Identity floor for a translated hit to count as the same gene. Streptococcal
#: surface proteins diverge, so this is deliberately permissive; what separates
#: "present" from "a different gene" here is the *reference* coverage plus the
#: panel's own ``min_identity`` where a gene needs a stricter floor.
DEFAULT_MIN_IDENTITY = 80.0
#: Reference coverage below which a hit is a fragment rather than a gene.
DEFAULT_MIN_COVERAGE = 60.0
#: Reference coverage at or above which a gene is called intact.
INTACT_COVERAGE = 95.0
#: Internal gap, in residues, still treated as allelic in-frame variation rather
#: than a lesion. The repeat regions of streptococcal surface proteins differ in
#: length between lineages by a few residues at a time.
INFRAME_INDEL = 10
#: How close to the end of a contig a fragment must sit before the assembly,
#: rather than the gene, is blamed for it. One read length is the natural scale.
CONTIG_EDGE_BP = 300

# --------------------------------------------------------------------- panels
@dataclass(frozen=True)
class PanelSource:
    """One RefSeq assembly the reference proteins are taken from."""

    accession: str      # GCF_000006785.2
    asm_name: str       # ASM678v2
    strain: str         # SF370
    note: str = ""

    @property
    def stem(self) -> str:
        return f"{self.accession}_{self.asm_name}"

    @property
    def ftp_dir(self) -> str:
        digits = self.accession.split("_", 1)[1].split(".")[0]
        prefix = self.accession.split("_", 1)[0]
        return (f"{NCBI_GENOMES}/{prefix}/{digits[0:3]}/{digits[3:6]}/"
                f"{digits[6:9]}/{self.stem}")


@dataclass(frozen=True)
class GeneSpec:
    """One gene the panel screens for, and how to name what is found."""

    name: str
    #: RefSeq gene symbols that identify this gene, lower-cased. Tried first.
    symbols: tuple[str, ...]
    #: Assembly the reference protein is taken from.
    source: str
    #: What the gene is, in words the report can print.
    product: str
    #: ``prophage`` | ``chromosome`` | ``variable`` | ``regulator`` | ``surface``
    location: str = "chromosome"
    #: What the call means clinically; printed beside the gene.
    meaning: str = ""
    #: Regular expression on the RefSeq product name, used when no symbol matches.
    product_re: str = ""
    #: Expected reference protein length in residues, as a sanity check on the
    #: pick. A mismatch is logged and the pick is still used, never silently.
    expect_aa: int = 0
    min_identity: float = DEFAULT_MIN_IDENTITY
    min_coverage: float = DEFAULT_MIN_COVERAGE
    #: Demand that the alignment reach both ends of the reference protein. Set
    #: for the regulators, where losing either terminus is loss of function: the
    #: MGAS5005 *covS* allele aligns from reference residue 23 to 500 - 95.6% of
    #: the protein, which a coverage rule alone calls intact - and the missing
    #: N-terminus is the published single-base deletion that makes it a
    #: pseudogene in RefSeq and defines the strain's invasive phenotype.
    strict_ends: bool = False
    #: Reference key in the module's REFERENCES that this entry rests on.
    reference: str = ""

    @property
    def key(self) -> str:
        return self.name


@dataclass
class GeneHit:
    """What the screen found for one panel gene in one assembly."""

    spec: GeneSpec
    accession: str = ""          # the reference protein the call is against
    contig: str = ""
    start: int = 0
    end: int = 0
    strand: str = "+"
    identity_pct: float = 0.0
    coverage_pct: float = 0.0
    ref_len: int = 0
    subject_spans: list[tuple[int, int]] = field(default_factory=list)
    n_hsps: int = 0
    internal_stops: int = 0
    frames: int = 1
    #: ``intact`` | ``truncated`` | ``internal_stop`` | ``fragment`` |
    #: ``split`` | ``absent``
    state: str = "absent"
    lesion: str = ""
    at_contig_edge: bool = False

    @property
    def present(self) -> bool:
        return self.state != "absent"

    @property
    def intact(self) -> bool:
        return self.state == "intact"

    def evidence(self) -> str:
        """One cell of evidence: where it is, how well it matches, what state."""
        if not self.present:
            return "not found"
        where = f"{self.contig}:{self.start}-{self.end}({self.strand})"
        cover = f"{self.coverage_pct:.0f}% of {self.ref_len} aa"
        out = f"{where} {self.identity_pct:.1f}% aa id, {cover} of {self.accession}"
        if self.lesion:
            out += f"; {self.lesion}"
        return out


# ------------------------------------------------------------------ transport
def download(url: str, dest: Path, *, attempts: int = 3) -> Path:
    """Fetch *url* to *dest*, retrying transport failures but not a 404."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    last: Exception | None = None
    for attempt in range(attempts):
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
                payload = response.read()
            if not payload:
                raise HydraError(f"{url} returned an empty file")
            dest.write_bytes(payload)
            return dest
        except urllib.error.HTTPError as exc:
            if 400 <= exc.code < 500 and exc.code != 429:
                raise HydraError(f"could not download {url}: {exc}") from exc
            last = exc
        except (urllib.error.URLError, OSError) as exc:
            last = exc
        if attempt + 1 < attempts:
            import time  # noqa: PLC0415 - only the retry path needs it

            time.sleep(2.0 * (attempt + 1))
    raise HydraError(f"could not download {url}: {last}")


def _gunzip(src: Path, dest: Path) -> Path:
    dest.write_bytes(gzip.decompress(src.read_bytes()))
    return dest


# ------------------------------------------------------------- panel assembly
def _feature_table(source: PanelSource, cache: Path) -> list[dict]:
    """RefSeq's own feature table for one assembly, as a list of CDS rows."""
    gz = cache / f"{source.accession}.feature_table.txt.gz"
    plain = cache / f"{source.accession}.feature_table.txt"
    if not plain.exists():
        if not gz.exists():
            download(f"{source.ftp_dir}/{source.stem}_feature_table.txt.gz", gz)
        _gunzip(gz, plain)
    rows: list[dict] = []
    with open(plain, encoding="utf-8") as handle:
        header = handle.readline().lstrip("#").strip().split("\t")
        for line in handle:
            values = line.rstrip("\n").split("\t")
            if len(values) < len(header):
                continue
            row = dict(zip(header, values))
            if row.get("feature") == "CDS":
                rows.append(row)
    return rows


def _proteins(source: PanelSource, cache: Path) -> dict[str, str]:
    """The assembly's protein FASTA, keyed by accession."""
    gz = cache / f"{source.accession}.protein.faa.gz"
    plain = cache / f"{source.accession}.protein.faa"
    if not plain.exists():
        if not gz.exists():
            download(f"{source.ftp_dir}/{source.stem}_protein.faa.gz", gz)
        _gunzip(gz, plain)
    return {header.split()[0]: seq for header, seq in read_fasta(plain)}


def fetch_panel(panel: tuple[GeneSpec, ...], sources: dict[str, PanelSource],
                out_dir: Path, *, name: str = "panel") -> Path:
    """Build the panel's protein BLAST database from RefSeq annotation.

    Idempotent: an existing, non-empty database with a manifest covering every
    panel gene is left alone. The manifest records which protein accession each
    gene resolved to, which is what the typing results cite.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = out_dir / f"{name}.faa"
    manifest = out_dir / f"{name}.manifest.json"
    if fasta.exists() and manifest.exists():
        try:
            known = json.loads(manifest.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            known = {}
        if all(spec.name in known for spec in panel):
            return fasta

    cache = out_dir / "refseq"
    cache.mkdir(parents=True, exist_ok=True)
    records: list[tuple[str, str]] = []
    resolved: dict[str, dict] = {}
    missing: list[str] = []
    for spec in panel:
        source = sources.get(spec.source)
        if source is None:
            raise HydraError(f"panel gene {spec.name} names an unknown source "
                             f"{spec.source!r}")
        picked = _pick(spec, _feature_table(source, cache))
        if picked is None:
            missing.append(spec.name)
            continue
        seqs = _proteins(source, cache)
        protein = seqs.get(picked["product_accession"])
        if not protein:
            missing.append(spec.name)
            continue
        if spec.expect_aa and abs(len(protein) - spec.expect_aa) > 2:
            LOG.warning("%s: RefSeq %s is %d aa, the panel expected %d",
                        spec.name, picked["product_accession"], len(protein),
                        spec.expect_aa)
        records.append((f"{spec.name}|{picked['product_accession']}", protein))
        resolved[spec.name] = {
            "accession": picked["product_accession"],
            "strain": source.strain,
            "assembly": source.accession,
            "locus_tag": picked.get("locus_tag", ""),
            "product": picked.get("name", ""),
            "length_aa": len(protein),
        }
    if missing:
        # Loud, and recorded: a gene the panel could not resolve is reported as
        # "not determined" rather than quietly becoming a negative call.
        LOG.warning("%s: %d panel gene(s) unresolved in RefSeq annotation: %s",
                    name, len(missing), ", ".join(missing))
    if not records:
        raise HydraError(f"no panel protein could be resolved for {name}")
    write_fasta(fasta, records)
    manifest.write_text(json.dumps({"genes": resolved, "unresolved": missing},
                                   indent=1, sort_keys=True), encoding="utf-8")
    make_blast_db(fasta, "prot")
    return fasta


def panel_manifest(out_dir: Path, name: str = "panel") -> dict:
    path = out_dir / f"{name}.manifest.json"
    if not path.exists():
        return {"genes": {}, "unresolved": []}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {"genes": {}, "unresolved": []}


def _pick(spec: GeneSpec, rows: list[dict]) -> dict | None:
    """The CDS row that is *spec* in one assembly's RefSeq annotation.

    Symbols first, product name second. A pseudogene has no protein accession
    in the feature table, so it cannot supply a reference and is skipped - which
    is correct: a truncated allele is a poor reference for detecting truncation.
    """
    wanted = {s.lower() for s in spec.symbols}
    by_symbol = [r for r in rows
                 if (r.get("symbol") or "").lower() in wanted
                 and r.get("product_accession")]
    if by_symbol:
        # A symbol can be shared by paralogues - SF370 annotates three hylP
        # copies - so the expected length decides between them when it is known.
        if spec.expect_aa:
            return min(by_symbol,
                       key=lambda r: abs(int(r.get("product_length") or 0) - spec.expect_aa))
        return max(by_symbol, key=lambda r: int(r.get("product_length") or 0))
    if spec.product_re:
        pattern = re.compile(spec.product_re, re.I)
        by_product = [r for r in rows
                      if pattern.search(r.get("name") or "")
                      and r.get("product_accession")]
        if spec.expect_aa and by_product:
            return min(by_product,
                       key=lambda r: abs(int(r.get("product_length") or 0) - spec.expect_aa))
        if by_product:
            return max(by_product, key=lambda r: int(r.get("product_length") or 0))
    return None


def make_blast_db(fasta: Path, dbtype: str) -> None:
    """Index *fasta* for BLAST, skipping the work when the index is current."""
    suffix = ".pin" if dbtype == "prot" else ".nin"
    index = Path(f"{fasta}{suffix}")
    if index.exists() and index.stat().st_mtime >= fasta.stat().st_mtime:
        return
    require("makeblastdb", "building a typing database")
    run(["makeblastdb", "-in", str(fasta), "-dbtype", dbtype,
         "-out", str(fasta), "-title", fasta.stem])


# -------------------------------------------------------------- the gene scan
def contig_lengths(assembly: Path) -> dict[str, int]:
    return {header.split()[0]: len(seq) for header, seq in read_fasta(assembly)}


def screen_panel(assembly: Path, fasta: Path, workdir: Path,
                 panel: tuple[GeneSpec, ...], *, threads: int = 1,
                 tag: str = "panel") -> tuple[dict[str, GeneHit],
                                              dict[str, list[GeneHit]]]:
    """blastx the assembly against the panel and describe every gene found.

    Returns the best hit per panel gene *and* every qualifying copy of it. Both
    are needed: the repertoire question is "does this isolate carry *speC*", for
    which one answer per gene is right, while describing prophage regions needs
    every copy of every phage module gene, since a group A streptococcus carries
    several prophages and each one has its own terminase and portal. Keeping only
    the best copy would let the module describe at most one prophage.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    by_name = {spec.name: spec for spec in panel}
    out_tab = workdir / f"{tag}.blastx.tsv"
    hsps = blast(
        "blastx", assembly, fasta, out_tab, threads=threads,
        evalue=1e-10, max_target_seqs=10000, with_alignment=True,
    )
    merged = merge_hsps(hsps, query_scale=3.0)
    # The alignment strings live on the HSPs, not on the merged hit, so keep
    # them to count premature stops afterwards. They have to be matched to the
    # merged hit by position: on a closed chromosome every HSP against one
    # reference shares the same query name, and pooling them would count a
    # spurious out-of-frame match elsewhere on the chromosome as a stop codon
    # inside this gene.
    aligned: dict[tuple[str, str], list] = {}
    for hsp in hsps:
        aligned.setdefault((hsp.qseqid, hsp.sseqid), []).append(hsp)

    lengths = contig_lengths(assembly)
    every: dict[str, list[GeneHit]] = {}
    for hit in merged:
        gene, _, accession = hit.sseqid.partition("|")
        spec = by_name.get(gene)
        if spec is None:
            continue
        if hit.identity_pct < spec.min_identity or hit.coverage_pct < spec.min_coverage:
            continue
        candidate = GeneHit(
            spec=spec, accession=accession, contig=hit.qseqid,
            start=hit.query_start, end=hit.query_end, strand=hit.strand,
            identity_pct=hit.identity_pct, coverage_pct=hit.coverage_pct,
            ref_len=hit.slen, subject_spans=list(hit.subject_intervals),
            n_hsps=hit.n_hsps,
        )
        own = _own_hsps(aligned.get((hit.qseqid, hit.sseqid), []), hit)
        candidate.internal_stops = _count_internal_stops(own, hit.slen)
        candidate.frames = _count_frames(own)
        describe_lesion(candidate, lengths.get(hit.qseqid, 0))
        every.setdefault(gene, []).append(candidate)

    best: dict[str, GeneHit] = {}
    for gene, found in every.items():
        found.sort(key=_rank, reverse=True)
        best[gene] = found[0]
        if len(found) > 1:
            note = (f"{len(found) - 1} further copy/copies of {gene} in this "
                    f"assembly: " + ", ".join(
                        f"{h.contig}:{h.start}-{h.end}" for h in found[1:4]))
            best[gene].lesion = (f"{best[gene].lesion}; {note}"
                                 if best[gene].lesion else note)
    for spec in panel:
        best.setdefault(spec.name, GeneHit(spec=spec))
        every.setdefault(spec.name, [])
    return best, every


def _own_hsps(hsps: list, hit) -> list:
    """The HSPs that belong to one merged hit: same strand, inside its span."""
    lo, hi = hit.query_start, hit.query_end
    return [h for h in hsps
            if h.strand == hit.strand and h.q_lo >= lo - 3 and h.q_hi <= hi + 3]


def _count_frames(hsps: list) -> int:
    """How many reading frames one gene's alignment used.

    This is how a frameshift shows itself in a translated search. The MGAS5005
    *covS* allele carries a single-base deletion (Sumby 2006): blastx aligns the
    gene in two pieces, consecutive on the reference but in different reading
    frames on the same strand, with no stop codon between them. Coverage alone
    would call that gene intact, and it is not - it is the lesion that defines
    the strain's invasive phenotype.
    """
    return len({h.q_lo % 3 for h in hsps}) if hsps else 1


def _rank(hit: GeneHit) -> tuple:
    return (hit.state == "intact", hit.coverage_pct, hit.identity_pct)


def _count_internal_stops(hsps: list, ref_len: int) -> int:
    """Stop codons inside the aligned part of the gene, from the translation.

    blastx translates the query, so a premature stop appears as ``*`` in the
    query string of the alignment. A stop aligned against the last few residues
    of the reference is the real stop codon and does not count.
    """
    total = 0
    for hsp in hsps:
        qseq = getattr(hsp, "qseq", "") or ""
        if "*" not in qseq:
            continue
        pos = min(hsp.sstart, hsp.send)
        for char, ref in zip(qseq, getattr(hsp, "sseq", "") or ""):
            if ref != "-":
                pos += 1
            if char == "*" and pos < ref_len - 2:
                total += 1
    return total


def describe_lesion(hit: GeneHit, contig_len: int) -> GeneHit:
    """Classify a hit as intact or lesioned, and say what the lesion is.

    The distinction that matters clinically is between a gene this isolate has
    broken and a gene this *assembly* has broken. A fragment whose alignment
    stops within a few hundred bases of a contig end is the second; the state is
    still reported, but the note names the assembly as the reason and the caller
    must not read it as a loss-of-function allele.
    """
    near_edge = (hit.start <= CONTIG_EDGE_BP
                 or (contig_len and contig_len - hit.end <= CONTIG_EDGE_BP))
    hit.at_contig_edge = bool(near_edge)
    missing_head = hit.subject_spans[0][0] - 1 if hit.subject_spans else 0
    missing_tail = (hit.ref_len - hit.subject_spans[-1][1]) if hit.subject_spans else 0
    internal_gap = max(
        (b[0] - a[1] - 1 for a, b in zip(hit.subject_spans, hit.subject_spans[1:])),
        default=0)

    if hit.internal_stops:
        hit.state = "internal_stop"
        hit.lesion = (f"{hit.internal_stops} premature stop codon(s) in the "
                      f"translated gene; product truncated")
        return hit
    if hit.frames > 1 and not near_edge:
        hit.state = "frameshift"
        hit.lesion = (f"the gene aligns in {hit.frames} reading frames on one "
                      f"strand with no stop between them, which is a frameshift; "
                      f"{hit.coverage_pct:.0f}% of the reference is covered")
        return hit
    ends_ok = (not hit.spec.strict_ends
               or (missing_head <= 5 and missing_tail <= 5))
    if (hit.coverage_pct >= INTACT_COVERAGE and internal_gap <= INFRAME_INDEL
            and ends_ok):
        # An in-frame indel of a few residues is allelic variation in these
        # species - the repeat regions of the collagen-like proteins and the M
        # protein differ in length between lineages - and calling it a lesion
        # would flag a normal allele as broken. It is still reported.
        hit.state = "intact"
        hit.lesion = (f"in-frame indel of {internal_gap} residue(s) relative to "
                      f"the reference allele" if internal_gap > 3 else "")
        return hit
    detail: list[str] = []
    if missing_tail > 3:
        detail.append(f"last {missing_tail} reference residues not aligned")
    if missing_head > 3:
        detail.append(f"first {missing_head} reference residues not aligned, "
                      f"consistent with a frameshift or a lost start codon")
    if internal_gap > INFRAME_INDEL:
        detail.append(f"internal deletion of {internal_gap} residues")
    if hit.n_hsps > 1:
        detail.append(f"{hit.n_hsps} alignment blocks")
    if hit.spec.strict_ends and not ends_ok and hit.coverage_pct >= INTACT_COVERAGE \
            and not near_edge:
        hit.state = "truncated"
        hit.lesion = "; ".join(detail)
        return hit
    if near_edge:
        hit.state = "fragment"
        detail.append("alignment reaches a contig end, so the assembly and not "
                      "the gene may be truncated")
    else:
        hit.state = "truncated"
    hit.lesion = "; ".join(detail)
    return hit


# ----------------------------------------------------------------- formatting
def span_string(spans: list[tuple[int, int]], length: int) -> str:
    joined = ",".join(f"{lo}-{hi}" for lo, hi in spans)
    return f"{joined}/{length}" if joined else f"-/{length}"


def pct(value: float) -> str:
    return f"{value:.1f}"


def blast_fields() -> tuple[str, ...]:
    """The tabular fields Hydra always asks BLAST for, for the docs to quote."""
    return OUTFMT_FIELDS
