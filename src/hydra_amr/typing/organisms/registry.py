"""Discovery and dispatch for the organism-specific typing modules.

Each module in this package types one organism and exposes the same handful of
names (``KEY``, ``TITLE``, ``SPECIES``, ``MLST_SCHEMES``, ``DATABASES``,
``REFERENCES``, ``fetch`` and ``run``). This module finds them, decides which
ones a sample should be put through, and runs them.

The selection rule, in order:

1. the species call, when Hydra made one and a module claims that species;
2. failing that, the MLST scheme, which is evidence of its own: a sample that
   typed against ``saureus`` is a *S. aureus* whatever the sketch said.

A module whose database is not installed is skipped with a line saying which
database is missing and how to install it, never with a silent absence and
never with a fabricated negative.
"""

from __future__ import annotations

import importlib
import pkgutil
import re
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Iterable, Sequence

from ...records import MlstCall, SpeciesCall, TypingResult
from ...utils import LOG, HydraError
from .base import OrganismDb, Reference

#: The organism modules, named explicitly. A frozen Windows build cannot walk
#: the package directory, so the list is the source of truth and the directory
#: scan below only adds modules that were dropped in later.
MODULE_NAMES: tuple[str, ...] = (
    "ecoli",
    "hinfluenzae",
    "nmeningitidis",
    "sagalactiae",
    "salmonella",
    "saureus",
    "spneumoniae",
    "spyogenes",
)

#: Files in this package that are shared helpers, not organism modules.
_NOT_A_MODULE = frozenset({"base", "registry", "strep_common"})

#: Where a module's data lives, relative to the Hydra database root. Keeping it
#: in its own subdirectory means ``db_dir/<KEY>/`` -- what every module's
#: ``fetch()`` and ``run()`` expect -- cannot collide with Hydra's own
#: ``nucl/``, ``prot/``, ``mlst/`` and ``typing/`` trees.
DB_SUBDIR = "organisms"

#: Readable headings for the schemes the modules report, used by the reports in
#: place of the bare scheme identifier. A scheme missing from here falls back to
#: its identifier with underscores turned into spaces, which reads acceptably
#: ("sag pilus"), so this table is presentation only and never load-bearing.
SCHEME_TITLES: dict[str, str] = {
    # Staphylococcus aureus
    "spa": "spa type (Ridom)",
    "agr": "agr group",
    "clonal_complex": "MLST clonal complex",
    "sccmec": "SCCmec type",
    "pvl": "Panton-Valentine leucocidin",
    "tst": "Toxic shock toxin (tst)",
    "enterotoxins": "Enterotoxin gene repertoire",
    "haemolysins": "Haemolysin gene repertoire",
    "hlb": "beta-haemolysin (hlb) status",
    "capsule": "Capsule serotype locus",
    "iec": "Immune evasion cluster",
    "acme": "Arginine catabolic mobile element",
    "bla": "beta-lactamase locus",
    # Streptococcus pyogenes
    "emm": "emm type (M protein)",
    "emm_virulome": "Virulence gene repertoire",
    "spy_prophage": "Prophage content",
    "spy_regulators": "Regulatory loci",
    "spy_adhesins": "Adhesins",
    # Streptococcus agalactiae
    "cps": "Capsular serotype",
    "alp": "Alpha-like protein (Alp) type",
    "sag_surface": "Surface protein repertoire",
    "sag_pilus": "Pilus island type",
    "sag_markers": "Marker gene repertoire",
    "sag_pbp2x": "PBP2X substitutions",
    # Streptococcus pneumoniae
    "serotype": "Capsular serotype",
    "pbp_type": "PBP type (1A:2B:2X)",
    "beta_lactam_mic": "Predicted beta-lactam MIC",
    "spn_resistance": "Non-beta-lactam resistance determinants",
    "spn_pili": "Pilus islets",
    "spn_mlst": "MLST",
    # Haemophilus influenzae / Neisseria meningitidis
    "capsule_type": "Capsule type",
    "serogroup": "Capsule serogroup",
    "pbp3": "PBP3 (ftsI) substitution pattern",
    "beta_lactam": "beta-lactam interpretation",
    "finetype": "Fine type",
    "bast": "Bexsero antigen sequence type",
    "mlst": "MLST",
    "amr": "Resistance determinants",
    # Salmonella
    "serovar": "Serovar (antigenic formula)",
    "salmonella_mlst": "MLST / eBurst group",
    "antigens": "Antigenic formula",
    "spi": "Pathogenicity islands",
    "markers": "Marker gene repertoire",
    # Escherichia coli / Shigella
    "phylogroup": "Phylogroup (Clermont)",
    "fimh": "fimH allele",
    "st131_subclone": "ST131 subclone",
    "stx": "Shiga toxin subtype",
    "pathotype": "Pathotype",
    "intimin": "Intimin (eae) subtype",
    "shigella_eiec": "Shigella / EIEC cluster",
    "serotype_oh": "O:H serotype",
    "virulence": "Virulence gene repertoire",
}

#: CDC's emm database also holds the group C and G types of
#: *Streptococcus dysgalactiae* subsp. *equisimilis*, named ``stG``, ``stC``,
#: ``stL`` and ``emmG``. The *S. pyogenes* module recognises them deliberately
#: and reports them under their own names, which is right -- but a reader who
#: sees one of them under a heading reading "emm type" would take it for a group
#: A emm type. The heading is corrected here instead.
_NON_PYOGENES_EMM = re.compile(r"^(st[GCL]|emmG)", re.IGNORECASE)

_NON_PYOGENES_EMM_NOTE = (
    "{call} is a group C/G type: the matching record comes from the group C and "
    "G section of the CDC emm database, which holds the M-protein gene types of "
    "Streptococcus dysgalactiae subsp. equisimilis. It is NOT a Streptococcus "
    "pyogenes (group A) emm type; confirm the species before reporting it"
)


@dataclass(frozen=True)
class OrganismModule:
    """One organism typing module, with the metadata needed to dispatch it."""

    key: str
    title: str
    species: tuple[str, ...]
    mlst_schemes: tuple[str, ...]
    module: ModuleType

    # -------------------------------------------------------------- metadata
    @property
    def databases(self) -> tuple[OrganismDb, ...]:
        return tuple(getattr(self.module, "DATABASES", ()))

    @property
    def references(self) -> tuple[Reference, ...]:
        return tuple(getattr(self.module, "REFERENCES", ()))

    def root(self, db_dir: Path | str) -> Path:
        """Where this module's data lives under *db_dir*."""
        return Path(db_dir) / self.key

    def missing_databases(self, db_dir: Path | str) -> list[str]:
        """The declared databases whose directory is absent under *db_dir*."""
        root = self.root(db_dir)
        return [db.name for db in self.databases if not (root / db.name).exists()]

    def installed(self, db_dir: Path | str) -> bool:
        """True when the module has something to work with.

        A module is run when its own directory exists and holds at least one of
        its databases. Modules degrade on purpose -- *H. influenzae* still types
        PBP3 without the *cap* database -- so requiring every database would
        withhold calls that can be made, while requiring none would run a
        module against an empty directory and report a page of
        "not determined".
        """
        root = self.root(db_dir)
        if not root.is_dir():
            return False
        if not self.databases:
            return True
        return any((root / db.name).exists() for db in self.databases)

    # ------------------------------------------------------------- execution
    def fetch(self, db_dir: Path | str) -> None:
        Path(db_dir).mkdir(parents=True, exist_ok=True)
        self.module.fetch(Path(db_dir))

    def run(self, assembly: Path | str, db_dir: Path | str, *, threads: int = 1,
            workdir: Path | str | None = None,
            mlst: MlstCall | None = None) -> list[TypingResult]:
        """Type one assembly with this module.

        ``mlst`` is Hydra's own MLST call. Modules that can use it accept an
        ``mlst`` keyword and are given it, which saves typing the housekeeping
        loci a second time; the rest are called without it.
        """
        kwargs: dict = {"threads": threads, "workdir": workdir}
        if mlst is not None and _accepts_mlst(self.module.run):
            kwargs["mlst"] = mlst
        results = self.module.run(Path(assembly), Path(db_dir), **kwargs)
        return relabel(list(results))


def _accepts_mlst(func) -> bool:
    try:
        import inspect

        return "mlst" in inspect.signature(func).parameters
    except (TypeError, ValueError):  # pragma: no cover - builtins only
        return False


def _load(name: str) -> OrganismModule | None:
    try:
        module = importlib.import_module(f"{__package__}.{name}")
    except Exception as exc:  # pragma: no cover - a broken module must not stop a run
        LOG.warning("organism typing module %s could not be imported: %s", name, exc)
        return None
    if not hasattr(module, "KEY") or not hasattr(module, "run"):
        return None
    return OrganismModule(
        key=str(module.KEY),
        title=str(getattr(module, "TITLE", module.KEY)),
        species=tuple(getattr(module, "SPECIES", ())),
        mlst_schemes=tuple(getattr(module, "MLST_SCHEMES", ())),
        module=module,
    )


_CACHE: dict[str, OrganismModule] | None = None


def modules(refresh: bool = False) -> dict[str, OrganismModule]:
    """Every organism typing module, keyed by ``KEY`` and sorted by key."""
    global _CACHE
    if _CACHE is not None and not refresh:
        return _CACHE
    names = list(MODULE_NAMES)
    # Pick up a module added to the package but not yet listed above. This is a
    # no-op in a frozen build, where iter_modules sees nothing.
    try:
        for info in pkgutil.iter_modules([str(Path(__file__).parent)]):
            if (info.name not in names and not info.name.startswith("_")
                    and info.name not in _NOT_A_MODULE):
                names.append(info.name)
    except OSError:  # pragma: no cover
        pass
    found: dict[str, OrganismModule] = {}
    for name in names:
        entry = _load(name)
        if entry is not None:
            found[entry.key] = entry
    _CACHE = dict(sorted(found.items()))
    return _CACHE


def get(key: str) -> OrganismModule:
    try:
        return modules()[key]
    except KeyError:
        known = ", ".join(sorted(modules())) or "none"
        raise HydraError(f"unknown organism typing module '{key}'. "
                         f"Known modules: {known}") from None


def database_keys() -> tuple[str, ...]:
    """The names ``hydra db`` uses for the organism typing databases."""
    return tuple(f"typing-{key}" for key in modules())


def key_from_database(name: str) -> str | None:
    """``typing-saureus`` -> ``saureus``; anything else -> ``None``."""
    if name.startswith("typing-") and name[len("typing-"):] in modules():
        return name[len("typing-"):]
    return None


# --------------------------------------------------------------------- select
def _norm(text: str) -> str:
    return re.sub(r"\s+", " ", str(text or "").strip().lower())


def select(species: SpeciesCall | None = None, mlst: MlstCall | None = None,
           *, only: Sequence[str] | None = None) -> list[OrganismModule]:
    """The modules that should be run for one sample, in key order.

    *only* restricts the answer to those module keys, for a user who asked for
    a particular module.
    """
    available = modules()
    wanted = set(only) if only else None
    chosen: dict[str, OrganismModule] = {}

    name = _norm(getattr(species, "name", "")) if species else ""
    binomial = ""
    if species is not None:
        binomial = _norm(f"{species.genus} {species.species}")
    for key, entry in available.items():
        claims = {_norm(s) for s in entry.species}
        if not claims:
            continue
        if name in claims or (binomial and binomial in claims):
            chosen[key] = entry
            continue
        # "Escherichia coli O157:H7 str. Sakai" is still an E. coli, and a
        # subspecies call ("Salmonella enterica subsp. enterica") still names
        # the species in front. Match a claim that the call starts with, at a
        # word boundary, so "Streptococcus pyogenes" never matches
        # "Streptococcus pneumoniae".
        for claim in claims:
            if name.startswith(claim + " "):
                chosen[key] = entry
                break

    if not chosen and mlst is not None and mlst.scheme and mlst.scheme != "-":
        scheme = _norm(mlst.scheme)
        for key, entry in available.items():
            if scheme in {_norm(s) for s in entry.mlst_schemes}:
                chosen[key] = entry

    if wanted is not None:
        chosen = {k: v for k, v in chosen.items() if k in wanted}
    return [chosen[k] for k in sorted(chosen)]


# ------------------------------------------------------------------- relabel
def relabel(results: Iterable[TypingResult]) -> list[TypingResult]:
    """Give every result a readable heading, and correct the misleading ones.

    Presentation only: no call, allele or note is replaced, only added to.
    """
    out: list[TypingResult] = []
    for result in results:
        title = SCHEME_TITLES.get(result.scheme)
        if title is None:
            title = result.scheme.replace("_", " ")
        note = result.note
        if result.scheme == "emm" and _NON_PYOGENES_EMM.match(str(result.call or "")):
            title = "M-protein gene type (group C/G, S. dysgalactiae)"
            extra = _NON_PYOGENES_EMM_NOTE.format(call=result.call)
            note = f"{note}; {extra}" if note else extra
        if title == result.title and note == result.note:
            out.append(result)
            continue
        out.append(TypingResult(
            scheme=result.scheme, call=result.call, lineage=result.lineage,
            alleles=dict(result.alleles), loci_found=result.loci_found,
            loci_total=result.loci_total, note=note, title=title))
    return out


# ----------------------------------------------------------------- top level
def type_assembly(assembly: Path | str, db_dir: Path | str, *,
                  species: SpeciesCall | None = None,
                  mlst: MlstCall | None = None,
                  threads: int = 1, workdir: Path | str | None = None,
                  only: Sequence[str] | None = None,
                  warn: list[str] | None = None) -> list[TypingResult]:
    """Run every module that applies to one assembly and pool the results.

    Anything that stops a module -- a missing database, a module that raises --
    is recorded in *warn* and in the returned results, so the report says what
    was not done and why rather than leaving a gap.
    """
    db_dir = Path(db_dir)
    out: list[TypingResult] = []
    for entry in select(species, mlst, only=only):
        if not entry.installed(db_dir):
            message = (f"{entry.title} typing skipped: its database is not "
                       f"installed. Install it with:  "
                       f"hydra db download typing-{entry.key}")
            LOG.info("%s", message)
            if warn is not None:
                warn.append(message)
            continue
        missing = entry.missing_databases(db_dir)
        try:
            results = entry.run(assembly, db_dir, threads=threads, workdir=workdir,
                                mlst=mlst)
        except HydraError as exc:
            message = f"{entry.title} typing failed: {exc}"
            LOG.warning("%s", message)
            if warn is not None:
                warn.append(message)
            continue
        except Exception as exc:  # pragma: no cover - a module bug, not a run-ender
            message = f"{entry.title} typing failed unexpectedly: {exc}"
            LOG.warning("%s", message, exc_info=True)
            if warn is not None:
                warn.append(message)
            continue
        if missing and warn is not None:
            warn.append(f"{entry.title} typing ran without its "
                        f"{', '.join(missing)} database(s); the calls that need "
                        f"them read 'not determined'")
        out.extend(results)
    return out
