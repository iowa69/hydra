"""Runtime configuration: database locations and screening thresholds."""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field
from pathlib import Path

from .utils import HydraError, cpu_count

PACKAGE_DIR = Path(__file__).resolve().parent
BUNDLED_DATA = PACKAGE_DIR / "data"

#: Environment variable that overrides the database root.
DB_ENV_VAR = "HYDRA_DB"

#: Cell values that ``--cell`` accepts for pivoted matrices.
CELL_MODES = ("binary", "identity", "coverage", "count", "genes", "depth", "fraction", "symbol")

#: Output formats ``--format`` accepts. ``genes`` and ``elements`` are the two
#: flat layouts kept for interoperability with existing downstream scripts.
OUTPUT_FORMATS = ("tsv", "csv", "json", "html", "xlsx", "genes", "elements")

#: Older spellings, still accepted so existing command lines keep working.
FORMAT_ALIASES = {"abricate": "genes", "amrfinder": "elements"}

#: Fields a pivot matrix can be keyed on, for both ``--rows`` and ``--columns``.
MATRIX_FIELDS = ("sample", "gene", "class", "subclass", "database", "element_type", "product")

#: Element types a hit can carry.
ELEMENT_TYPES = ("AMR", "VIRULENCE", "STRESS", "PLASMID", "SEROTYPE")


def default_db_dir() -> Path:
    """Resolve the database root.

    Order: ``$HYDRA_DB`` -> ``~/.hydra/db`` -> bundled ``data/db``.
    """
    env = os.environ.get(DB_ENV_VAR)
    if env:
        return Path(env).expanduser().resolve()
    user_dir = Path.home() / ".hydra" / "db"
    if user_dir.exists():
        return user_dir
    bundled = BUNDLED_DATA / "db"
    if bundled.exists():
        return bundled
    return user_dir


@dataclass
class Thresholds:
    """Identity/coverage cut-offs, expressed as percentages."""

    min_identity: float = 80.0
    min_coverage: float = 60.0
    #: A hit shorter than this fraction of the reference is flagged PARTIAL.
    partial_coverage: float = 90.0
    #: Translated (protein) search defaults, matching AMRFinderPlus behaviour.
    protein_min_identity: float = 90.0
    protein_min_coverage: float = 90.0
    #: Coverage floor below which a protein hit is dropped rather than called partial.
    protein_partial_coverage: float = 50.0
    #: Point mutations demand near-full-length, high-identity reference alignment.
    mutation_min_identity: float = 90.0
    mutation_min_coverage: float = 90.0
    #: Read mode: minimum depth at a locus before a call is trusted.
    min_depth: int = 5
    #: Read mode: minimum breadth of coverage (%) for a gene to be reported.
    min_gene_breadth: float = 80.0
    #: Heteroresistance: minimum minor-allele fraction to report.
    min_allele_fraction: float = 0.02
    #: Heteroresistance: minimum reads supporting the minor allele.
    min_allele_reads: int = 3
    #: Fraction at or above which a mutation is considered fixed rather than hetero.
    fixed_allele_fraction: float = 0.90
    #: Names the user set on the command line. Per-database defaults never
    #: override these, so ``--min-identity 60`` is honoured even where a database
    #: would normally impose a stricter cut-off.
    explicit: set = field(default_factory=set)

    def set_explicit(self, name: str, value) -> None:
        setattr(self, name, value)
        self.explicit.add(name)

    def was_set(self, name: str) -> bool:
        return name in self.explicit

    def validate(self) -> None:
        for name in ("min_identity", "min_coverage", "partial_coverage",
                     "protein_min_identity", "protein_min_coverage", "protein_partial_coverage",
                     "mutation_min_identity", "mutation_min_coverage", "min_gene_breadth"):
            value = getattr(self, name)
            if not 0.0 <= value <= 100.0:
                raise HydraError(f"--{name.replace('_', '-')} must be between 0 and 100 (got {value})")
        for name in ("min_allele_fraction", "fixed_allele_fraction"):
            value = getattr(self, name)
            if not 0.0 <= value <= 1.0:
                raise HydraError(f"--{name.replace('_', '-')} must be between 0 and 1 (got {value})")
        if self.min_depth < 1:
            raise HydraError("--min-depth must be >= 1")
        if self.min_allele_reads < 1:
            raise HydraError("--min-allele-reads must be >= 1")
        if self.min_allele_fraction >= self.fixed_allele_fraction:
            raise HydraError(
                "--min-allele-fraction must be below --fixed-allele-fraction "
                f"({self.min_allele_fraction} >= {self.fixed_allele_fraction})"
            )


@dataclass
class Config:
    """Everything a run needs that is not a per-sample input."""

    db_dir: Path = field(default_factory=default_db_dir)
    databases: list[str] = field(default_factory=lambda: ["ncbi"])
    threads: int = field(default_factory=cpu_count)
    jobs: int = 0  # 0 -> auto: derived from threads
    thresholds: Thresholds = field(default_factory=Thresholds)
    organism: str | None = None
    auto_organism: bool = True
    plus: bool = False
    keep_temp: bool = False
    tmp_dir: Path | None = None
    #: Report every overlapping hit rather than only the best one per locus.
    report_overlaps: bool = False
    #: Maximum reference-length fraction two hits may share before deduplication.
    overlap_fraction: float = 0.5

    def __post_init__(self) -> None:
        self.db_dir = Path(self.db_dir)
        if self.threads < 1:
            self.threads = 1
        if self.jobs < 0:
            self.jobs = 0

    def worker_layout(self, n_samples: int) -> tuple[int, int]:
        """Split the thread budget into (parallel samples, threads per sample)."""
        if n_samples <= 0:
            return 1, self.threads
        if self.jobs > 0:
            jobs = min(self.jobs, n_samples)
        else:
            # BLAST scales poorly past ~4 threads; prefer sample-level parallelism.
            jobs = min(n_samples, max(1, self.threads // 2), 16)
        jobs = max(1, jobs)
        per = max(1, self.threads // jobs)
        return jobs, per

    def as_dict(self) -> dict:
        data = asdict(self)
        data["db_dir"] = str(self.db_dir)
        data["tmp_dir"] = str(self.tmp_dir) if self.tmp_dir else None
        return data
