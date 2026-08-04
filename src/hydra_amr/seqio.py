"""FASTA / FASTQ reading and lightweight sequence statistics."""

from __future__ import annotations

import gzip
import lzma
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Sequence

from .utils import HydraError, is_gzip, smart_open

FASTA_EXTS = (".fasta", ".fa", ".fna", ".ffn", ".fsa", ".seq", ".fas", ".contigs")
FASTQ_EXTS = (".fastq", ".fq")


def _base_ext(path: str | Path) -> str:
    name = str(path).lower()
    for comp in (".gz", ".bz2", ".xz"):
        if name.endswith(comp):
            name = name[: -len(comp)]
            break
    return os.path.splitext(name)[1]


def looks_like_fastq(path: str | Path) -> bool:
    if _base_ext(path) in FASTQ_EXTS:
        return True
    return _sniff(path) == "fastq"


def looks_like_fasta(path: str | Path) -> bool:
    if _base_ext(path) in FASTA_EXTS:
        return True
    return _sniff(path) == "fasta"


def _sniff(path: str | Path) -> str | None:
    """Peek at the first non-blank line to decide FASTA vs FASTQ."""
    try:
        with smart_open(path, "rt") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    return "fasta"
                if line.startswith("@"):
                    return "fastq"
                return None
    except (OSError, UnicodeDecodeError):
        return None
    return None


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield (header, sequence) pairs. Header keeps the full description line."""
    header: str | None = None
    chunks: list[str] = []
    with smart_open(path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def read_fasta_dict(path: str | Path, key: str = "id") -> dict[str, str]:
    """Load a FASTA into a dict keyed by first whitespace token (``id``) or full header."""
    out: dict[str, str] = {}
    for header, seq in read_fasta(path):
        name = header.split()[0] if key == "id" else header
        out[name] = seq
    return out


def write_fasta(path: str | Path, records: Sequence[tuple[str, str]], wrap: int = 60) -> None:
    with open(path, "w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            if wrap and wrap > 0:
                for i in range(0, len(seq), wrap):
                    handle.write(seq[i:i + wrap] + "\n")
            else:
                handle.write(seq + "\n")


def read_fastq_head(path: str | Path, max_records: int = 10000) -> Iterator[tuple[str, str, str]]:
    """Yield up to *max_records* (name, seq, qual) tuples from a FASTQ."""
    with smart_open(path, "rt") as handle:
        count = 0
        while count < max_records:
            name = handle.readline()
            if not name:
                return
            seq = handle.readline().strip()
            handle.readline()
            qual = handle.readline().strip()
            if not qual:
                return
            yield name[1:].strip(), seq, qual
            count += 1


@dataclass
class AssemblyStats:
    """Basic contiguity/composition metrics for an assembly."""

    contigs: int = 0
    total_length: int = 0
    n50: int = 0
    n90: int = 0
    largest_contig: int = 0
    gc_percent: float = 0.0
    ambiguous_bases: int = 0
    lengths: list[int] = field(default_factory=list, repr=False)

    def as_dict(self) -> dict:
        return {
            "contigs": self.contigs,
            "total_length": self.total_length,
            "n50": self.n50,
            "n90": self.n90,
            "largest_contig": self.largest_contig,
            "gc_percent": round(self.gc_percent, 2),
            "ambiguous_bases": self.ambiguous_bases,
        }


def assembly_stats(path: str | Path) -> AssemblyStats:
    stats = AssemblyStats()
    gc = 0
    at = 0
    for _, seq in read_fasta(path):
        length = len(seq)
        if length == 0:
            continue
        stats.contigs += 1
        stats.total_length += length
        stats.lengths.append(length)
        upper = seq.upper()
        gc += upper.count("G") + upper.count("C")
        at += upper.count("A") + upper.count("T")
        stats.ambiguous_bases += length - (
            upper.count("A") + upper.count("C") + upper.count("G") + upper.count("T")
        )
    if not stats.lengths:
        return stats
    stats.lengths.sort(reverse=True)
    stats.largest_contig = stats.lengths[0]
    known = gc + at
    stats.gc_percent = (100.0 * gc / known) if known else 0.0
    running = 0
    for length in stats.lengths:
        running += length
        if not stats.n50 and running >= stats.total_length * 0.5:
            stats.n50 = length
        if not stats.n90 and running >= stats.total_length * 0.9:
            stats.n90 = length
            break
    return stats


def fastq_stats(path: str | Path, sample_records: int = 20000) -> dict:
    """Estimate read count, length and quality from the head of a FASTQ.

    Read count is extrapolated from compressed/uncompressed byte ratios, so it is
    an estimate; ``reads_estimated`` marks whether extrapolation was used.
    """
    lengths: list[int] = []
    qual_sum = 0
    qual_bases = 0
    bytes_consumed = 0
    n_records = 0
    for name, seq, qual in read_fastq_head(path, sample_records):
        n_records += 1
        lengths.append(len(seq))
        qual_sum += sum(ord(c) - 33 for c in qual)
        qual_bases += len(qual)
        # '@', '+', four newlines: name + seq + qual + 6 bytes per record.
        bytes_consumed += len(name) + len(seq) + len(qual) + 6
    if not lengths:
        return {"reads": 0, "mean_read_length": 0, "mean_quality": 0.0, "reads_estimated": False}
    mean_len = sum(lengths) / len(lengths)
    file_size = os.path.getsize(path)
    estimated = False
    if n_records >= sample_records:
        estimated = True
        # gzip on FASTQ typically lands near 3.5x; uncompressed uses the byte count directly
        effective = file_size * (3.5 if is_gzip(path) else 1.0)
        per_record = bytes_consumed / n_records
        reads = int(effective / per_record) if per_record else n_records
    else:
        reads = n_records
    return {
        "reads": reads,
        "mean_read_length": round(mean_len, 1),
        "mean_quality": round(qual_sum / qual_bases, 1) if qual_bases else 0.0,
        "reads_estimated": estimated,
        "bases": int(reads * mean_len),
    }


class _Unreadable(HydraError):
    """A file that cannot be decompressed or decoded."""


def _read_error(path: str | Path, exc: Exception) -> HydraError:
    name = str(path)
    if isinstance(exc, (gzip.BadGzipFile, EOFError, lzma.LZMAError)) or "bz2" in type(exc).__module__:
        return _Unreadable(
            f"{name} could not be decompressed ({exc}).\n"
            f"The file is either truncated or not compressed despite its name.")
    if isinstance(exc, UnicodeDecodeError):
        return _Unreadable(
            f"{name} is not readable as text ({exc}).\n"
            f"If it is compressed, give it the matching .gz, .bz2 or .xz suffix.")
    return _Unreadable(f"{name} could not be read: {exc}")


def validate_fasta(path: str | Path) -> None:
    """Raise HydraError when *path* is not a readable, non-empty FASTA."""
    p = Path(path)
    if not p.exists():
        raise HydraError(f"input not found: {path}")
    if p.stat().st_size == 0:
        raise HydraError(f"input file is empty: {path}")
    try:
        # Scan on: a leading placeholder record with no sequence does not make
        # the rest of the assembly unusable.
        for _, seq in read_fasta(path):
            if seq:
                return
        if _sniff(path) == "fastq":
            raise HydraError(
                f"{path} is named like an assembly but its contents are FASTQ.\n"
                f"Pass reads with --reads, or with -1/-2 for a pair.")
    except (OSError, EOFError, UnicodeDecodeError, lzma.LZMAError) as exc:
        raise _read_error(path, exc) from exc
    raise HydraError(f"no sequence records found in FASTA: {path}")


def validate_fastq(path: str | Path) -> None:
    """Raise HydraError when *path* is not a readable, non-empty FASTQ."""
    p = Path(path)
    if not p.exists():
        raise HydraError(f"input not found: {path}")
    if p.stat().st_size == 0:
        raise HydraError(f"input file is empty: {path}")
    try:
        for _name, seq, qual in read_fastq_head(path, max_records=1):
            if seq and len(seq) == len(qual):
                return
            raise HydraError(f"first FASTQ record in {path} is malformed "
                             f"(sequence and quality lengths differ)")
    except (OSError, EOFError, UnicodeDecodeError, lzma.LZMAError) as exc:
        raise _read_error(path, exc) from exc
    if _sniff(path) == "fasta":
        raise HydraError(f"{path} is named like reads but its contents are FASTA.\n"
                         f"Pass assemblies with -a/--assembly.")
    raise HydraError(f"no reads found in FASTQ: {path}")
