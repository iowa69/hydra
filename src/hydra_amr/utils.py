"""Shared low-level helpers: logging, subprocess, temp files, sequence math."""

from __future__ import annotations

import gzip
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Sequence

LOG = logging.getLogger("hydra")

_LEVEL_COLOURS = {
    "DEBUG": "\033[36m",
    "INFO": "\033[32m",
    "WARNING": "\033[33m",
    "ERROR": "\033[31m",
    "CRITICAL": "\033[1;31m",
}
_RESET = "\033[0m"


class _Formatter(logging.Formatter):
    def __init__(self, colour: bool):
        super().__init__("[%(asctime)s] %(levelname)s %(message)s", datefmt="%H:%M:%S")
        self.colour = colour

    def format(self, record: logging.LogRecord) -> str:
        text = super().format(record)
        if self.colour:
            prefix = _LEVEL_COLOURS.get(record.levelname, "")
            if prefix:
                text = text.replace(record.levelname, f"{prefix}{record.levelname}{_RESET}", 1)
        return text


def setup_logging(verbosity: int = 0, quiet: bool = False) -> None:
    """Configure the package logger. verbosity: 0=INFO, >=1=DEBUG."""
    level = logging.DEBUG if verbosity >= 1 else logging.INFO
    if quiet:
        level = logging.WARNING
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(_Formatter(colour=sys.stderr.isatty()))
    LOG.handlers[:] = [handler]
    LOG.setLevel(level)
    LOG.propagate = False


class HydraError(RuntimeError):
    """Fatal, user-facing error. The CLI prints these without a traceback."""


def die(message: str, code: int = 1) -> None:
    raise HydraError(message)


def run(
    cmd: Sequence[str],
    *,
    stdin_data: str | None = None,
    stdout_path: str | Path | None = None,
    check: bool = True,
    capture: bool = True,
    env: dict | None = None,
) -> subprocess.CompletedProcess:
    """Run an external command, logging it at DEBUG level.

    Raises HydraError with the captured stderr when the command fails.
    """
    cmd = [str(c) for c in cmd]
    LOG.debug("exec: %s", " ".join(cmd))
    started = time.time()
    stdout_handle = None
    try:
        if stdout_path is not None:
            stdout_handle = open(stdout_path, "w")
            stdout_target = stdout_handle
        else:
            stdout_target = subprocess.PIPE if capture else None
        proc = subprocess.run(
            cmd,
            input=stdin_data,
            stdout=stdout_target,
            stderr=subprocess.PIPE if capture else None,
            text=True,
            env=env,
        )
    except FileNotFoundError as exc:
        raise HydraError(
            f"required executable not found: {cmd[0]}\n"
            f"Install it, or re-create the environment with:\n"
            f"  conda install -c bioconda -c conda-forge {conda_package(cmd[0])}"
        ) from exc
    finally:
        if stdout_handle is not None:
            stdout_handle.close()
    LOG.debug("done: %s (%.2fs, rc=%d)", cmd[0], time.time() - started, proc.returncode)
    if check and proc.returncode != 0:
        err = (proc.stderr or "").strip()
        raise HydraError(f"{cmd[0]} failed (exit {proc.returncode}):\n{err[:4000]}")
    return proc


#: Executables whose conda package is named differently.
CONDA_PACKAGES = {
    "blastn": "blast", "blastx": "blast", "tblastn": "blast",
    "makeblastdb": "blast", "blastdbcmd": "blast", "blast_formatter": "blast",
    "spades.py": "spades", "bcftools": "bcftools",
}


def conda_package(tool: str) -> str:
    """The conda package that provides *tool*."""
    return CONDA_PACKAGES.get(tool, tool)


def have(tool: str) -> bool:
    """True when *tool* is on PATH."""
    return shutil.which(tool) is not None


def require(tool: str, why: str = "") -> str:
    path = shutil.which(tool)
    if path is None:
        extra = f" ({why})" if why else ""
        raise HydraError(
            f"required executable '{tool}' not found on PATH{extra}.\n"
            f"  conda install -c bioconda -c conda-forge {conda_package(tool)}"
        )
    return path


@contextmanager
def tempdir(prefix: str = "hydra.", keep: bool = False, parent: str | Path | None = None):
    path = Path(tempfile.mkdtemp(prefix=prefix, dir=str(parent) if parent else None))
    try:
        yield path
    finally:
        if keep:
            LOG.debug("keeping temp dir %s", path)
        else:
            shutil.rmtree(path, ignore_errors=True)


def smart_open(path: str | Path, mode: str = "rt"):
    """Open a plain, gzip, bzip2 or xz file transparently.

    Text reads replace undecodable bytes rather than raising. Sequence data is
    ASCII, but description lines are free text and some published databases carry
    stray bytes in them (VFDB has a non-breaking space); losing a character in a
    product description beats aborting the run with a UnicodeDecodeError.
    """
    name = str(path)
    text = {"errors": "replace"} if "b" not in mode else {}
    if name.endswith(".gz"):
        return gzip.open(name, mode, **text)
    if name.endswith(".bz2"):
        import bz2  # noqa: PLC0415 - only needed for this suffix

        return bz2.open(name, mode, **text)
    if name.endswith(".xz"):
        import lzma  # noqa: PLC0415 - only needed for this suffix

        return lzma.open(name, mode, **text)
    return open(name, mode, **text)


def is_gzip(path: str | Path) -> bool:
    try:
        with open(path, "rb") as handle:
            return handle.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def human_bytes(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024.0:
            return f"{n:.0f}{unit}" if unit == "B" else f"{n:.1f}{unit}"
        n /= 1024.0
    return f"{n:.1f}PB"


def human_time(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, secs = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)}m{secs:04.1f}s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)}h{int(minutes)}m"


_COMPLEMENT = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


_CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(seq: str, table_start_is_met: bool = True) -> str:
    """Translate a nucleotide sequence in frame 1 using the standard code.

    Alternative bacterial start codons (GTG/TTG/ATT/CTG) are rendered as M when
    *table_start_is_met* is set, matching how reference protein records are stored.
    """
    seq = seq.upper().replace("U", "T")
    out = []
    for i in range(0, len(seq) - 2, 3):
        out.append(_CODONS.get(seq[i:i + 3], "X"))
    if table_start_is_met and out:
        out[0] = "M"
    return "".join(out)


def chunked(items: Iterable, size: int):
    """Yield lists of at most *size* items."""
    batch = []
    for item in items:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


def cpu_count() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def free_bytes(path: str | Path) -> int:
    """Bytes available on the filesystem holding *path* (0 when unknown)."""
    try:
        stats = os.statvfs(str(path))
    except (OSError, AttributeError):
        return 0
    return stats.f_bavail * stats.f_frsize


def check_workspace(path: str | Path, input_bytes: int, factor: float = 1.5,
                    floor: int = 512 * 1024 * 1024) -> None:
    """Warn when the working directory may not have room for the intermediates.

    Tabular BLAST output for a large batch runs to a multiple of the input, and
    on many systems the temporary directory is a memory-backed filesystem that
    is much smaller than the disk. Running out of space part-way through wastes
    the whole run, so say something before it starts.
    """
    available = free_bytes(path)
    if not available:
        return
    needed = max(floor, int(input_bytes * factor))
    if available < needed:
        LOG.warning(
            "only %s free on the filesystem holding %s, and this run may write about %s "
            "of intermediates. Use --tmpdir to point at somewhere roomier.",
            human_bytes(available), path, human_bytes(needed))


def natural_key(text: str) -> tuple:
    """Sort key that orders embedded numbers numerically.

    ``blaSHV-2`` sorts before ``blaSHV-11``, which plain string ordering gets
    backwards. Used to break ties between equally-scoring database alleles so
    the reported name is the canonical low-numbered one and does not depend on
    the order BLAST happened to emit its hits in.
    """
    parts: list = []
    number = ""
    for ch in str(text):
        if ch.isdigit():
            number += ch
        else:
            if number:
                parts.append((1, int(number), ""))
                number = ""
            parts.append((0, 0, ch))
    if number:
        parts.append((1, int(number), ""))
    return tuple(parts)


def safe_name(text: str) -> str:
    """Filesystem/column-safe token."""
    keep = []
    for ch in text:
        keep.append(ch if (ch.isalnum() or ch in "._-") else "_")
    return "".join(keep).strip("_") or "sample"


def sample_name_from_path(path: str | Path, strip_read_suffix: bool = False) -> str:
    """Derive a sample name from a file path by stripping known extensions."""
    name = Path(path).name
    for ext in (".gz", ".bz2", ".xz"):
        if name.endswith(ext):
            name = name[: -len(ext)]
    for ext in (".fasta", ".fa", ".fna", ".ffn", ".fsa", ".seq", ".fastq", ".fq", ".contigs", ".scaffolds"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    if strip_read_suffix:
        for suffix in ("_R1_001", "_R2_001", "_R1", "_R2", "_1", "_2", ".R1", ".R2", ".1", ".2"):
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
    return name or "sample"
