"""Thin, typed wrapper around the BLAST+ programs plus HSP-merging helpers.

Hydra always asks BLAST for the same tabular fields so every caller can rely on
the same keys, and merges the HSPs of a query/subject pair on the *subject*
axis: reference coverage is what matters when deciding whether a gene is
present, and a gene split by an assembly gap still adds up to full coverage.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Sequence

from ..utils import LOG, HydraError, require, run

#: Tabular fields requested from BLAST, in order.
OUTFMT_FIELDS = (
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "slen", "qlen", "gaps", "nident",
)
#: Appended when the caller needs the gapped alignment strings (mutation calling).
ALIGNMENT_FIELDS = ("qseq", "sseq")
OUTFMT = "6 " + " ".join(OUTFMT_FIELDS)

_INT_FIELDS = {"length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
               "slen", "qlen", "gaps", "nident"}
_FLOAT_FIELDS = {"pident", "evalue", "bitscore"}


@dataclass
class Hsp:
    """One high-scoring segment pair, with subject coordinates normalised."""

    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    slen: int
    qlen: int
    gaps: int
    nident: int
    strand: str
    qseq: str = ""
    sseq: str = ""

    @property
    def sub_lo(self) -> int:
        return min(self.sstart, self.send)

    @property
    def sub_hi(self) -> int:
        return max(self.sstart, self.send)

    @property
    def q_lo(self) -> int:
        return min(self.qstart, self.qend)

    @property
    def q_hi(self) -> int:
        return max(self.qstart, self.qend)


def parse_tabular(path: str | Path, with_alignment: bool = False) -> list[Hsp]:
    """Parse a BLAST tabular file written with :data:`OUTFMT`."""
    fields = OUTFMT_FIELDS + (ALIGNMENT_FIELDS if with_alignment else ())
    out: list[Hsp] = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < len(fields):
                continue
            values: dict = {}
            for name, raw in zip(fields, parts):
                if name in _INT_FIELDS:
                    try:
                        values[name] = int(raw)
                    except ValueError:
                        values[name] = 0
                elif name in _FLOAT_FIELDS:
                    try:
                        values[name] = float(raw)
                    except ValueError:
                        values[name] = 0.0
                else:
                    values[name] = raw
            strand = "+" if values["send"] >= values["sstart"] else "-"
            out.append(Hsp(strand=strand, **values))
    return out


def blast(
    program: str,
    query: str | Path,
    db: str | Path,
    out: str | Path,
    *,
    threads: int = 1,
    evalue: float = 1e-20,
    perc_identity: float | None = None,
    max_target_seqs: int = 10000,
    task: str | None = None,
    culling_limit: int | None = None,
    extra: Sequence[str] = (),
    word_size: int | None = None,
    with_alignment: bool = False,
) -> list[Hsp]:
    """Run a BLAST program and return the parsed HSPs.

    ``perc_identity`` is only accepted by the nucleotide programs; it is ignored
    for blastx/tblastn where BLAST+ rejects the flag. Set ``with_alignment`` to
    also retrieve the gapped alignment strings needed for mutation calling.
    """
    require(program, "sequence search")
    fields = OUTFMT_FIELDS + (ALIGNMENT_FIELDS if with_alignment else ())
    cmd: list[str] = [
        program, "-query", str(query), "-db", str(db), "-out", str(out),
        "-outfmt", "6 " + " ".join(fields), "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs), "-num_threads", str(max(1, threads)),
    ]
    if task:
        cmd += ["-task", task]
    if program in ("blastn", "megablast", "dc-megablast"):
        cmd += ["-dust", "no"]
        if perc_identity is not None:
            cmd += ["-perc_identity", str(perc_identity)]
    else:
        cmd += ["-seg", "no"]
    if word_size:
        cmd += ["-word_size", str(word_size)]
    if culling_limit:
        cmd += ["-culling_limit", str(culling_limit)]
    cmd += list(extra)
    run(cmd)
    hsps = parse_tabular(out, with_alignment=with_alignment)
    if program in ("blastx", "tblastn"):
        # The subject of a translated search is a protein, whose coordinates
        # always ascend; the reading frame lives in the query coordinates.
        hsps = [replace(h, strand="+" if h.qend >= h.qstart else "-") for h in hsps]
    LOG.debug("%s: %d HSPs against %s", program, len(hsps), Path(db).name)
    return hsps


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge 1-based inclusive intervals, returning them sorted and disjoint."""
    ordered = sorted((min(a, b), max(a, b)) for a, b in intervals)
    if not ordered:
        return []
    merged = [list(ordered[0])]
    for lo, hi in ordered[1:]:
        if lo <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [(lo, hi) for lo, hi in merged]


def interval_length(intervals: Iterable[tuple[int, int]]) -> int:
    return sum(hi - lo + 1 for lo, hi in intervals)


def overlap_length(a: tuple[int, int], b: tuple[int, int]) -> int:
    lo = max(min(a), min(b))
    hi = min(max(a), max(b))
    return max(0, hi - lo + 1)


@dataclass
class MergedHit:
    """All HSPs between one query sequence and one subject, combined."""

    qseqid: str
    sseqid: str
    slen: int
    qlen: int
    identity_pct: float
    coverage_pct: float
    covered_bases: int
    subject_intervals: list[tuple[int, int]]
    query_start: int
    query_end: int
    strand: str
    gaps: int
    bitscore: float
    evalue: float
    n_hsps: int

    @property
    def coverage_string(self) -> str:
        """Reference-span notation, e.g. ``1-861/861`` or ``1-400,500-861/861``."""
        spans = ",".join(f"{lo}-{hi}" for lo, hi in self.subject_intervals)
        return f"{spans}/{self.slen}"

    @property
    def match_score(self) -> float:
        """Percentage of the reference matched identically.

        Ranking alleles by bit score favours whichever reference happens to be
        longest, which picks the wrong variant name when several near-identical
        alleles align to the same locus. This combines coverage and identity
        into one length-independent number instead.
        """
        return self.coverage_pct * self.identity_pct / 100.0


def _query_gap(a: Hsp, b: Hsp) -> int:
    """Distance between two HSPs on the query, 0 when they touch or overlap."""
    if a.q_hi < b.q_lo:
        return b.q_lo - a.q_hi
    if b.q_hi < a.q_lo:
        return a.q_lo - b.q_hi
    return 0


def merge_hsps(hsps: Sequence[Hsp], *, max_hsp_gap: int | None = None,
               query_scale: float = 1.0, span_slack: int = 5000) -> list[MergedHit]:
    """Group HSPs by (query, subject) and merge them on the subject axis.

    HSPs on opposing strands are kept apart, since a real gene copy does not
    change orientation mid-alignment; the stronger orientation wins.

    Merging is also constrained on the query axis. Without that, two truncated
    copies of the same gene hundreds of kilobases apart on one contig - each
    covering a different half of the reference - would add up to a single
    "complete" gene spanning the whole region. ``query_scale`` is the query
    length expected per unit of reference (3 for a translated search).
    """
    grouped: dict[tuple[str, str, str], list[Hsp]] = {}
    for hsp in hsps:
        grouped.setdefault((hsp.qseqid, hsp.sseqid, hsp.strand), []).append(hsp)

    best_by_pair: dict[tuple[str, str], MergedHit] = {}
    for (qseqid, sseqid, strand), group in grouped.items():
        group.sort(key=lambda h: h.bitscore, reverse=True)
        slen = group[0].slen
        # A real locus cannot occupy much more query sequence than the reference
        # it matches; allow generous slack for assembly gaps and repeats.
        allowed_span = int(slen * query_scale * 2) + span_slack
        kept: list[Hsp] = []
        used: list[tuple[int, int]] = []
        for hsp in group:
            span = (hsp.sub_lo, hsp.sub_hi)
            # Skip an HSP that mostly repeats subject territory already claimed.
            if any(overlap_length(span, prev) > 0.5 * (span[1] - span[0] + 1) for prev in used):
                continue
            if kept:
                if max_hsp_gap is not None and min(_query_gap(hsp, k) for k in kept) > max_hsp_gap:
                    continue
                q_lo = min([hsp.q_lo] + [k.q_lo for k in kept])
                q_hi = max([hsp.q_hi] + [k.q_hi for k in kept])
                if q_hi - q_lo + 1 > allowed_span:
                    continue
            kept.append(hsp)
            used.append(span)
        if not kept:
            continue
        intervals = merge_intervals((h.sub_lo, h.sub_hi) for h in kept)
        covered = interval_length(intervals)
        # Weight each HSP's identity by the reference it newly covers, so a pair
        # of partly overlapping HSPs does not count their shared bases twice.
        claimed: list[tuple[int, int]] = []
        weighted_ident = 0.0
        weighted_len = 0.0
        for hsp in kept:
            span = (hsp.sub_lo, hsp.sub_hi)
            hsp_len = span[1] - span[0] + 1
            already = sum(overlap_length(span, prev) for prev in claimed)
            fresh = max(0, hsp_len - already)
            if fresh and hsp_len:
                share = fresh / hsp_len
                weighted_ident += hsp.nident * share
                weighted_len += hsp.length * share
            claimed.append(span)
        identity = (100.0 * weighted_ident / weighted_len) if weighted_len else 0.0
        merged = MergedHit(
            qseqid=qseqid, sseqid=sseqid, slen=slen, qlen=kept[0].qlen,
            identity_pct=identity,
            coverage_pct=(100.0 * covered / slen) if slen else 0.0,
            covered_bases=covered, subject_intervals=intervals,
            query_start=min(h.q_lo for h in kept), query_end=max(h.q_hi for h in kept),
            strand=strand, gaps=sum(h.gaps for h in kept),
            bitscore=sum(h.bitscore for h in kept), evalue=min(h.evalue for h in kept),
            n_hsps=len(kept),
        )
        key = (qseqid, sseqid)
        current = best_by_pair.get(key)
        if current is None or merged.bitscore > current.bitscore:
            best_by_pair[key] = merged
    return list(best_by_pair.values())


def deduplicate(hits: Sequence, *, key_span, key_seq, key_score,
                overlap_fraction: float = 0.5, key_tiebreak=None) -> list:
    """Keep the best-scoring hit among those overlapping on the same sequence.

    Used to stop a single locus being reported once per near-identical database
    allele (the classic ``blaTEM-1A``/``blaTEM-1B`` duplication). *key_tiebreak*
    makes the winner deterministic when scores are identical, which they often
    are for allele families: without it the reported name would depend on the
    order BLAST emitted its hits in.
    """
    ordered = list(hits)
    if key_tiebreak is not None:
        ordered.sort(key=key_tiebreak)
    ordered = sorted(ordered, key=key_score, reverse=True)
    kept: list = []
    claimed: dict[str, list[tuple[int, int]]] = {}
    for hit in ordered:
        seq = key_seq(hit)
        span = key_span(hit)
        length = max(1, span[1] - span[0] + 1)
        conflict = False
        for prev in claimed.get(seq, []):
            if overlap_length(span, prev) >= overlap_fraction * length:
                conflict = True
                break
        if conflict:
            continue
        kept.append(hit)
        claimed.setdefault(seq, []).append(span)
    return kept


def check_db_exists(fasta: str | Path, dbtype: str = "nucl") -> None:
    """Raise unless *fasta* has a usable BLAST index beside it.

    A large database is split into numbered volumes (``.00.nin``) with an alias
    file (``.nal``/``.pal``), so all three shapes count as installed.
    """
    suffix = ".nin" if dbtype == "nucl" else ".pin"
    alias = ".nal" if dbtype == "nucl" else ".pal"
    fasta = Path(fasta)
    candidates = [f"{fasta}{suffix}", f"{fasta}.00{suffix}", f"{fasta}{alias}"]
    if not any(Path(candidate).exists() for candidate in candidates):
        raise HydraError(f"BLAST index missing for {fasta}; run 'hydra db import --force'")
