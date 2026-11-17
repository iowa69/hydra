"""Helpers shared by the *Haemophilus influenzae* and *Neisseria meningitidis*
organism typing modules.

Both organisms are typed the same way: a small reference set of protein and
nucleotide sequences is screened against the assembly, and every call is made
from the resulting evidence rather than from gene names alone. The helpers here
are therefore deliberately generic -- locate a gene, pull its coding sequence
out of the assembly, translate it, say how it is broken, and match a peptide
against a downloaded allele definition -- and the organism modules hold the
decision rules.

Only BLAST+ is used, so the same code runs inside the frozen Windows build.
"""

from __future__ import annotations

import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

from ...engines.blast import MergedHit, blast, merge_hsps
from ...seqio import read_fasta, write_fasta
from ...utils import LOG, HydraError, require, revcomp, run, translate

#: Upstream is asked to identify the client, as the rest of Hydra does.
USER_AGENT = "hydra-amr (+https://github.com/iowa69/hydra)"
TIMEOUT = 180

#: Bacterial start codons. An alternative start is translated as methionine,
#: which is how reference protein records store it.
START_CODONS = ("ATG", "GTG", "TTG", "ATT", "CTG")
STOP_CODONS = ("TAA", "TAG", "TGA")


# --------------------------------------------------------------------- fetch
def download(url: str, dest: Path, *, attempts: int = 4) -> Path:
    """Download *url* to *dest*, retrying transport failures and throttling.

    A PubMLST install is a few dozen requests and a dropped one would leave a
    locus silently empty, so a transient failure is retried; a 404 is answered
    once, because retrying a wrong URL only wastes the user's time.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    last: Exception | None = None
    for attempt in range(attempts):
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
                payload = response.read()
            if not payload:
                raise HydraError(f"{url} returned an empty response")
            dest.write_bytes(payload)
            return dest
        except urllib.error.HTTPError as exc:
            last = HydraError(f"could not download {url}: {exc}")
            if exc.code not in (429, 500, 502, 503, 504):
                raise last from exc
            time.sleep(5.0 * (attempt + 1))
        except (urllib.error.URLError, OSError) as exc:
            last = HydraError(f"could not download {url}: {exc}")
            time.sleep(1.5 * (attempt + 1))
    raise last if last else HydraError(f"could not download {url}")


#: Read-only REST interface to the PubMLST sequence-definition databases.
PUBMLST_REST = "https://rest.pubmlst.org"


def pubmlst_locus_fasta(database: str, locus: str, dest: Path) -> Path:
    """Download every allele of one PubMLST locus as FASTA."""
    locus_path = urllib.parse.quote(locus, safe="")
    url = f"{PUBMLST_REST}/db/{database}/loci/{locus_path}/alleles_fasta"
    return download(url, dest)


def pubmlst_profiles(database: str, scheme: int, dest: Path) -> Path:
    """Download one PubMLST scheme's profile table as tab-separated text."""
    url = f"{PUBMLST_REST}/db/{database}/schemes/{scheme}/profiles_csv"
    return download(url, dest)


def ncbi_cds_fasta(accession: str, dest: Path) -> Path:
    """Download the annotated coding sequences of one NCBI nucleotide record.

    ``fasta_cds_na`` carries the ``[gene=...]``/``[protein=...]`` qualifiers in
    its description lines, which is what lets a reference locus be picked out by
    the name the literature uses for it without parsing GenBank flat files.
    """
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={urllib.parse.quote(accession)}"
           "&rettype=fasta_cds_na&retmode=text")
    return download(url, dest)


def ncbi_protein_fasta(accessions: Sequence[str], dest: Path) -> Path:
    """Download protein records from NCBI by accession, in one request."""
    return _ncbi_fasta("protein", accessions, dest)


def ncbi_nucleotide_fasta(accessions: Sequence[str], dest: Path) -> Path:
    """Download nucleotide records from NCBI by accession, in one request."""
    return _ncbi_fasta("nuccore", accessions, dest)


def _ncbi_fasta(database: str, accessions: Sequence[str], dest: Path) -> Path:
    ids = urllib.parse.quote(",".join(accessions))
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db={database}&id={ids}&rettype=fasta&retmode=text")
    download(url, dest)
    text = dest.read_text(errors="replace")
    if not text.lstrip().startswith(">"):
        # E-utilities answers a throttled or malformed request with a JSON or
        # XML error body and HTTP 200, which would otherwise be written out as
        # if it were sequence.
        raise HydraError(f"NCBI returned no FASTA for {', '.join(accessions)}: "
                         f"{text.strip()[:200]}")
    return dest


_QUALIFIER = re.compile(r"\[(\w+)=([^\]]*)\]")


def cds_qualifiers(header: str) -> dict[str, str]:
    """Parse the ``[key=value]`` qualifiers of an NCBI CDS FASTA header."""
    return {key: value for key, value in _QUALIFIER.findall(header)}


# ----------------------------------------------------------------- sequences
def load_contigs(assembly: Path) -> dict[str, str]:
    """Load an assembly, keyed by the first whitespace token of each header."""
    contigs: dict[str, str] = {}
    for header, sequence in read_fasta(assembly):
        contigs[header.split()[0]] = sequence.upper()
    if not contigs:
        raise HydraError(f"no sequences found in {assembly}")
    return contigs


def makeblastdb(fasta: Path, dbtype: str = "nucl") -> Path:
    """Build a BLAST index beside *fasta* and return the database path."""
    require("makeblastdb", "organism typing")
    run(["makeblastdb", "-in", str(fasta), "-dbtype", dbtype, "-out", str(fasta)])
    return fasta


def best_reading_frame(nucleotide: str) -> str:
    """Translate a reference CDS in whichever orientation has no internal stop.

    Some published reference sets store a gene as the plus-strand slice of the
    genome it was taken from rather than in the gene's own orientation. Picking
    the frame by trial means the protein reference is right either way instead of
    being 30 stop codons long.
    """
    candidates = []
    for rank, sequence in ((0, nucleotide), (1, revcomp(nucleotide))):
        trimmed = sequence[:len(sequence) - len(sequence) % 3]
        protein = translate(trimmed).rstrip("*")
        # A coding sequence runs from a start codon to a stop codon. Preferring
        # the orientation that does, and the given orientation on a tie, stops a
        # reverse complement whose junk frame happens to be a residue longer
        # from winning over the real reading frame.
        framed = 0 if (trimmed[:3] in START_CODONS
                       and trimmed[-3:] in STOP_CODONS) else 1
        candidates.append((protein.count("*"), framed, rank, protein))
    candidates.sort()
    stops, _framed, _rank, protein = candidates[0]
    return protein if stops == 0 else ""


def write_subset(records: Iterable[tuple[str, str]], dest: Path) -> int:
    """Write (header, sequence) pairs to *dest*, returning how many were written."""
    items = list(records)
    dest.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(dest, items, wrap=60)
    return len(items)


# ------------------------------------------------------------------ screening
@dataclass
class GeneHit:
    """One reference locus located in the assembly, with its evidence."""

    gene: str
    identity_pct: float
    coverage_pct: float
    contig: str
    start: int                  # 1-based, ascending, on the contig
    end: int
    strand: str
    #: Reference-span notation, e.g. ``1-861/861``.
    coverage: str = ""
    bitscore: float = 0.0
    n_hsps: int = 1
    reference: str = ""         # the reference record the hit was made against
    #: Translated subject as BLAST aligned it (tblastn only); ``*`` marks a stop
    #: codon inside the alignment.
    aligned_subject: str = ""
    aligned_query: str = ""
    query_start: int = 0        # first reference residue/base the hit covers
    query_end: int = 0
    note: str = ""

    @property
    def span(self) -> tuple[int, int]:
        return (self.start, self.end)

    def evidence(self) -> str:
        return (f"{self.identity_pct:.1f}% id, {self.coverage_pct:.0f}% cov, "
                f"{self.contig}:{self.start}-{self.end}{self.strand}")


def _hit_from_merged(hit: MergedHit, gene: str, contig_of: dict[str, str]) -> GeneHit:
    return GeneHit(
        gene=gene,
        identity_pct=hit.identity_pct,
        coverage_pct=hit.coverage_pct,
        contig=contig_of.get(hit.qseqid, hit.qseqid),
        start=hit.query_start,
        end=hit.query_end,
        strand=hit.strand,
        coverage=hit.coverage_string,
        bitscore=hit.bitscore,
        n_hsps=hit.n_hsps,
        reference=hit.sseqid,
    )


def screen_proteins(assembly_db: Path, proteins: Path, workdir: Path, *,
                    threads: int = 1, min_identity: float = 70.0,
                    min_coverage: float = 60.0, tag: str = "prot") -> list[GeneHit]:
    """tblastn a reference protein set against an assembly's BLAST database.

    Translated search is the right tool for a presence/absence call on a
    divergent locus: the capsule biosynthesis genes of two serogroups can sit
    below the useful limit of blastn while their products still align over their
    whole length.
    """
    out = workdir / f"{tag}.tblastn.tsv"
    # Single-threaded on purpose. BLAST+ 2.16.0 segfaults deterministically on a
    # translated search when -max_target_seqs is combined with -num_threads > 1
    # and the search returns no hits, which is exactly what happens when a
    # sample of the wrong species is screened against a small protein reference
    # set. Hydra's blast() always passes -max_target_seqs, so the only fix
    # available without editing shared code is to drop the thread count; these
    # reference sets are a few hundred kilobases, so the cost is small.
    del threads
    hsps = blast("tblastn", proteins, assembly_db, out, threads=1,
                 evalue=1e-10, max_target_seqs=10000, with_alignment=True)
    if not hsps:
        return []
    # tblastn reports the protein as the query, so merging has to be done on the
    # *protein* axis: that is the reference whose coverage decides presence.
    merged = merge_hsps([_swap(h) for h in hsps], query_scale=3.0)
    hits: list[GeneHit] = []
    best_hsp = _best_hsp_by_pair(hsps)
    for hit in merged:
        if hit.identity_pct < min_identity or hit.coverage_pct < min_coverage:
            continue
        gene = hit.sseqid.split()[0]
        found = GeneHit(
            gene=gene, identity_pct=hit.identity_pct, coverage_pct=hit.coverage_pct,
            contig=hit.qseqid, start=hit.query_start, end=hit.query_end,
            strand=hit.strand, coverage=hit.coverage_string, bitscore=hit.bitscore,
            n_hsps=hit.n_hsps, reference=gene,
        )
        anchor = best_hsp.get((hit.qseqid, hit.sseqid))
        if anchor is not None:
            found.aligned_query = anchor.qseq
            found.aligned_subject = anchor.sseq
            found.query_start = min(anchor.qstart, anchor.qend)
            found.query_end = max(anchor.qstart, anchor.qend)
        hits.append(found)
    return hits


def _swap(hsp):
    """Exchange query and subject of a translated HSP.

    ``merge_hsps`` merges on the subject axis and Hydra's convention is that the
    reference is the subject. For tblastn the reference protein is BLAST's
    query, so the two are exchanged before merging and the protein becomes the
    subject again, which makes ``coverage_pct`` mean coverage of the reference
    protein in both directions.
    """
    from dataclasses import replace

    return replace(
        hsp,
        qseqid=hsp.sseqid, sseqid=hsp.qseqid,
        qstart=hsp.sstart, qend=hsp.send,
        sstart=hsp.qstart, send=hsp.qend,
        slen=hsp.qlen, qlen=hsp.slen,
        qseq=hsp.sseq, sseq=hsp.qseq,
        strand="+" if hsp.send >= hsp.sstart else "-",
    )


def _best_hsp_by_pair(hsps) -> dict[tuple[str, str], object]:
    """Highest-scoring HSP per (contig, protein) pair, keeping its alignment."""
    best: dict[tuple[str, str], object] = {}
    for hsp in hsps:
        key = (hsp.sseqid, hsp.qseqid)   # contig, protein
        current = best.get(key)
        if current is None or hsp.bitscore > current.bitscore:
            best[key] = hsp
    return best


def screen_nucleotide(assembly: Path, reference_db: Path, workdir: Path, *,
                      threads: int = 1, min_identity: float = 80.0,
                      min_coverage: float = 60.0, tag: str = "nucl",
                      task: str = "blastn") -> list[GeneHit]:
    """blastn an assembly against a reference nucleotide database.

    Every copy that passes the thresholds is returned, because copy number is
    itself a result here: a duplicated *H. influenzae* cap locus and a single
    one are different calls.
    """
    out = workdir / f"{tag}.blastn.tsv"
    hsps = blast("blastn", assembly, reference_db, out, threads=threads,
                 evalue=1e-20, task=task,
                 perc_identity=max(0.0, min_identity - 10.0), max_target_seqs=10000)
    merged = merge_hsps(hsps)
    hits: list[GeneHit] = []
    for hit in merged:
        if hit.identity_pct < min_identity or hit.coverage_pct < min_coverage:
            continue
        hits.append(_hit_from_merged(hit, hit.sseqid.split()[0], {}))
    return hits


def best_per_gene(hits: Sequence[GeneHit]) -> dict[str, GeneHit]:
    """Keep the strongest hit per gene name, ranked on matched reference."""
    best: dict[str, GeneHit] = {}
    for hit in hits:
        score = hit.coverage_pct * hit.identity_pct
        current = best.get(hit.gene)
        if current is None or score > current.coverage_pct * current.identity_pct:
            best[hit.gene] = hit
    return best


# ------------------------------------------------------------ coding sequence
@dataclass
class Cds:
    """A coding sequence pulled out of an assembly around a BLAST hit."""

    nucleotide: str
    protein: str
    contig: str
    start: int
    end: int
    strand: str
    complete: bool              # a start codon and an in-frame stop were found
    truncated_at_contig_end: bool = False
    internal_stops: int = 0
    note: str = ""


def extract_cds(contigs: dict[str, str], contig: str, start: int, end: int,
                strand: str, *, ref_offset: int = 0,
                max_extend: int = 3000) -> Cds | None:
    """Recover the coding sequence overlapping a hit, in the hit's frame.

    The hit's own span is not the gene: a translated search stops at the first
    stop codon and starts wherever the alignment did. The frame comes from the
    hit; the start codon is then the in-frame candidate nearest to where the
    reference says the gene should begin (*ref_offset* is the number of
    reference residues preceding the aligned region), bounded by the nearest
    upstream in-frame stop. Taking the furthest upstream start instead would
    lengthen the product whenever a random GTG or TTG happens to sit in frame
    just before the real initiation codon, and every residue number downstream
    would then be wrong.

    A gene running off the end of a contig is reported as truncated rather than
    as absent or as a loss of function, which an assembly cannot distinguish.
    """
    sequence = contigs.get(contig)
    if not sequence:
        return None
    length = len(sequence)
    start = max(1, start)
    end = min(length, end)
    if strand == "-":
        local = revcomp(sequence)
        lo = length - end + 1
        hi = length - start + 1
    else:
        local = sequence
        lo, hi = start, end

    frame_start = lo - 1                          # 0-based, in frame with the hit
    expected = frame_start - 3 * max(0, ref_offset)
    # Collect the in-frame start codons upstream of the hit, stopping at the
    # first in-frame stop codon, then take the one closest to where the
    # reference puts the start of the gene.
    candidates: list[int] = []
    limit = max(0, expected - max_extend)
    position = frame_start
    ran_off_contig = False
    while position >= limit:
        codon = local[position:position + 3]
        if len(codon) < 3:
            break
        if codon in STOP_CODONS and position != frame_start:
            break
        if codon in START_CODONS:
            candidates.append(position)
        position -= 3
    else:
        ran_off_contig = limit == 0
    if position < 0:
        ran_off_contig = True
    begin = frame_start
    found_start = bool(candidates)
    if candidates:
        begin = min(candidates, key=lambda p: (abs(p - expected), p))
    hit_contig_start = ran_off_contig and not candidates

    # Walk downstream from the chosen start to the first in-frame stop.
    internal_stops = 0
    stop_at: int | None = None
    position = begin
    while position + 3 <= len(local):
        codon = local[position:position + 3]
        if codon in STOP_CODONS:
            if position + 3 < hi:
                # A stop before the end of the alignment is a real internal stop.
                internal_stops += 1
            stop_at = position + 3
            break
        position += 3
    truncated = stop_at is None
    nucleotide = local[begin:stop_at if stop_at is not None else len(local) - (len(local) - begin) % 3]
    if not nucleotide:
        return None
    protein = translate(nucleotide).rstrip("*")
    if strand == "-":
        abs_start = length - (stop_at if stop_at is not None else len(local)) + 1
        abs_end = length - begin
    else:
        abs_start = begin + 1
        abs_end = (stop_at if stop_at is not None else len(local))
    return Cds(
        nucleotide=nucleotide, protein=protein, contig=contig,
        start=abs_start, end=abs_end, strand=strand,
        complete=found_start and stop_at is not None,
        truncated_at_contig_end=truncated or hit_contig_start,
        internal_stops=internal_stops,
    )


def substitutions(reference: str, query: str, aligned_query: str = "",
                  aligned_subject: str = "", offset: int = 0) -> list[str]:
    """List ``N526K``-style differences of *query* against *reference*.

    When the two proteins are the same length the comparison is positional,
    which is the usual case for a housekeeping gene such as *ftsI* or *gyrA*.
    Otherwise the gapped alignment strings from a translated BLAST search are
    walked instead, so an indel does not shift every downstream residue and turn
    one substitution into fifty.
    """
    out: list[str] = []
    if reference and query and len(reference) == len(query):
        for index, (ref, obs) in enumerate(zip(reference, query), start=1):
            if ref != obs:
                out.append(f"{ref}{index}{obs}")
        return out
    if aligned_query and aligned_subject:
        position = offset
        for ref, obs in zip(aligned_query, aligned_subject):
            if ref != "-":
                position += 1
            if ref != "-" and obs != "-" and ref != obs:
                out.append(f"{ref}{position}{obs}")
        return out
    return out


def peptide_index(fasta: Path, locus: str = "") -> dict[str, str]:
    """Map allele/peptide sequence -> identifier, from a PubMLST FASTA.

    PubMLST names a record ``<locus>_<id>``; the identifier alone is what the
    nomenclature uses (``P1.7-2`` is PorA VR1 allele ``7-2``).
    """
    index: dict[str, str] = {}
    for header, sequence in read_fasta(fasta):
        name = header.split()[0]
        identifier = name
        if locus and name.lower().startswith(locus.lower()):
            identifier = name[len(locus):].lstrip("_-")
        elif "_" in name:
            identifier = name.rsplit("_", 1)[1]
        clean = sequence.strip().upper().rstrip("*")
        if clean:
            index.setdefault(clean, identifier)
    return index


def match_peptide(protein: str, index: dict[str, str]) -> tuple[str, str]:
    """Identify a peptide allele, exactly or as the nearest definition.

    Returns ``(identifier, note)``. An exact match carries no note; otherwise
    the closest definition of the same length is named with the number of
    differences, because a peptide one residue away from a known one is a new
    peptide and must not be reported under the old one's number.
    """
    if not protein:
        return "", "not determined: peptide not recovered from the assembly"
    clean = protein.strip().upper().rstrip("*")
    exact = index.get(clean)
    if exact:
        return exact, ""
    best_id = ""
    best_diff = 10 ** 6
    for sequence, identifier in index.items():
        if len(sequence) != len(clean):
            continue
        diff = sum(1 for a, b in zip(sequence, clean) if a != b)
        if diff < best_diff:
            best_diff, best_id = diff, identifier
    if best_id:
        return "new", (f"peptide not in the downloaded definitions; closest is {best_id} "
                       f"({best_diff} residue difference"
                       f"{'s' if best_diff != 1 else ''})")
    return "new", "peptide not in the downloaded definitions and no same-length neighbour"


def find_variable_region(protein: str, index: dict[str, str], *,
                         window: tuple[float, float] = (0.0, 1.0)) -> tuple[str, str]:
    """Identify a variable region by finding its peptide inside *protein*.

    PubMLST defines PorA VR1/VR2 and FetA VR as short peptide sequences, so the
    variant is read off the assembly by locating the defined peptide in the
    translated gene rather than by guessing at the region's coordinates, which
    move with every indel in the loop. The longest definition that occurs in the
    expected part of the protein wins; a shorter one can be a substring of it.
    """
    if not protein:
        return "", "not determined: gene product not recovered from the assembly"
    lo = int(len(protein) * window[0])
    hi = int(len(protein) * window[1])
    region = protein[lo:hi]
    best_id = ""
    best_len = 0
    for sequence, identifier in index.items():
        if len(sequence) > best_len and sequence in region:
            best_id, best_len = identifier, len(sequence)
    if best_id:
        return best_id, ""
    return "new", ("no defined variable-region peptide found in this product; "
                   "the region is new or the gene is incomplete")


def profile_table(path: Path) -> tuple[list[str], dict[tuple[str, ...], str], dict[str, str]]:
    """Load a PubMLST ``profiles_csv`` table.

    Returns the locus column names, a profile -> ST map, and the per-ST
    clonal-complex column when the scheme publishes one. Which columns are loci
    is decided by position: PubMLST puts the ST first, the loci next and the
    metadata columns last, and the metadata names are known.
    """
    meta_names = {"clonal_complex", "species", "mlst_clade", "lineage", "cc",
                  "clade", "mlst_cluster", "subspecies", "st_complex"}
    loci: list[str] = []
    table: dict[tuple[str, ...], str] = {}
    complexes: dict[str, str] = {}
    if not path.exists():
        return loci, table, complexes
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        indexes = [i for i, name in enumerate(header)
                   if i > 0 and name.strip().lower() not in meta_names]
        loci = [header[i].strip() for i in indexes]
        cc_index = next((i for i, name in enumerate(header)
                         if name.strip().lower() in ("clonal_complex", "cc")), None)
        for line in handle:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) <= max(indexes, default=0):
                continue
            st = row[0].strip()
            table[tuple(row[i].strip() for i in indexes)] = st
            if cc_index is not None and len(row) > cc_index:
                value = row[cc_index].strip()
                if value:
                    complexes[st] = value
    return loci, table, complexes


@dataclass
class MlstProfile:
    """A seven-locus sequence type read from the assembly."""

    scheme: str
    sequence_type: str = "-"
    clonal_complex: str = ""
    alleles: dict[str, str] = field(default_factory=dict)
    loci_found: int = 0
    loci_total: int = 0
    note: str = ""


def call_mlst(assembly: Path, allele_db: Path, profiles: Path, scheme: str, *,
              workdir: Path, threads: int = 1,
              min_identity: float = 90.0) -> MlstProfile:
    """Assign an ST from a single-scheme allele database.

    The organism modules need the clonal complex as part of their own output
    (a meningococcal fine type is normally quoted with it), and the scheme is
    known from the species, so the lookup is done here rather than by inspecting
    Hydra's scheme-agnostic MLST result. Allele naming follows PubMLST's own
    ``<locus>_<id>`` convention.
    """
    loci, table, complexes = profile_table(profiles)
    out = workdir / f"{scheme}.mlst.tsv"
    hsps = blast("blastn", assembly, allele_db, out, threads=threads,
                 evalue=1e-20, task="megablast", word_size=20,
                 perc_identity=min_identity - 5.0, max_target_seqs=100000)
    best: dict[str, tuple[float, float, str]] = {}
    for hit in merge_hsps(hsps):
        if hit.identity_pct < min_identity:
            continue
        name = hit.sseqid.split()[0]
        if "_" not in name:
            continue
        locus, allele = name.rsplit("_", 1)
        key = (hit.coverage_pct >= 100.0 and hit.identity_pct >= 100.0,
               hit.coverage_pct, hit.identity_pct)
        current = best.get(locus)
        if current is None or key > (current[0] >= 100.0 and current[1] >= 100.0,
                                     current[0], current[1]):
            best[locus] = (hit.coverage_pct, hit.identity_pct, allele)
    if not loci:
        loci = sorted(best)
    alleles: dict[str, str] = {}
    exact = 0
    for locus in loci:
        found = best.get(locus)
        if found is None:
            alleles[locus] = "-"
        elif found[0] >= 100.0 and found[1] >= 100.0:
            alleles[locus] = found[2]
            exact += 1
        elif found[0] >= 100.0:
            alleles[locus] = f"~{found[2]}"
        else:
            alleles[locus] = f"{found[2]}?"
    st = table.get(tuple(alleles[locus] for locus in loci), "-")
    note = ""
    if st == "-":
        missing = [locus for locus in loci if alleles[locus] == "-"]
        inexact = [locus for locus in loci
                   if alleles[locus] != "-" and not alleles[locus].isdigit()]
        parts = []
        if missing:
            parts.append(f"{', '.join(missing)} not found")
        if inexact:
            parts.append(f"novel allele at {', '.join(inexact)}")
        note = ("; ".join(parts) + "; no ST assigned") if parts else \
            "allele combination absent from the profile table (novel ST)"
    return MlstProfile(scheme=scheme, sequence_type=st,
                       clonal_complex=complexes.get(st, ""), alleles=alleles,
                       loci_found=exact, loci_total=len(loci), note=note)


def describe_gene(hit: GeneHit | None, cds: Cds | None = None) -> str:
    """One human-readable phrase for a locus and the state it was found in."""
    if hit is None:
        return "not detected"
    parts = [f"{hit.identity_pct:.1f}% identity over {hit.coverage_pct:.0f}% of the reference"]
    parts.append(f"{hit.contig}:{hit.start}-{hit.end}({hit.strand})")
    if cds is not None:
        if cds.internal_stops:
            parts.append(f"{cds.internal_stops} internal stop codon(s)")
        if cds.truncated_at_contig_end:
            parts.append("runs into a contig end")
        if not cds.complete:
            parts.append("no complete open reading frame recovered")
    if hit.n_hsps > 1:
        parts.append(f"{hit.n_hsps} alignment blocks")
    return "; ".join(parts)


def log_missing(what: str, why: str) -> None:
    LOG.debug("organism typing: %s not determined (%s)", what, why)
