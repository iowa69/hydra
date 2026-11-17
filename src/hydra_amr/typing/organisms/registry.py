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

import gzip
import importlib
import pkgutil
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Iterable, Sequence

from ...records import MlstCall, SpeciesCall, TypingResult
from ...utils import LOG, HydraError, tempdir
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
    "pvl": "Panton-Valentine leucocidin (lukS/lukF-PV)",
    "tst": "Toxic shock toxin (tst)",
    "enterotoxins": "Enterotoxin gene repertoire",
    "haemolysins": "Haemolysin gene repertoire",
    "hlb": "beta-haemolysin (hlb) status",
    "iec": "Immune evasion cluster (IEC type)",
    "acme": "Arginine catabolic mobile element",
    "blaZ": "blaZ penicillinase",
    "bla_operon": "blaZ operon (blaZ/blaR1/blaI)",
    "capsule": "Capsule serotype (cap locus)",
    "saureus": "S. aureus typing",
    # Streptococcus pyogenes
    "emm": "emm type (M protein)",
    "emm_virulome": "Virulence and superantigen repertoire",
    "spy_prophage": "Prophage content",
    "spy_regulators": "Regulatory loci",
    "spy_adhesins": "Adhesins",
    # Streptococcus agalactiae
    "sag_alp": "Alpha-like protein (Alp) type",
    "sag_surface": "Surface protein repertoire",
    "sag_pilus": "Pilus island type",
    "sag_markers": "Marker gene repertoire",
    "sag_pbp2x": "PBP2X substitutions",
    # Streptococcus pneumoniae
    "cps_serotype": "Capsular serotype",
    "pbp_type": "PBP type (1A:2B:2X)",
    "beta_lactam_mic": "Predicted beta-lactam MIC",
    "target_substitutions": "Resistance target substitutions",
    "macrolide_mechanism": "Macrolide resistance mechanism",
    "other_determinants": "Other resistance determinants",
    "pilus_islets": "Pilus islets",
    # Haemophilus influenzae
    "cap": "Capsule type (cap locus)",
    "pbp3": "PBP3 (ftsI) substitution pattern",
    "beta_lactam": "beta-lactam interpretation",
    # Neisseria meningitidis
    "finetype": "Fine type",
    "bast": "Bexsero antigen sequence type (BAST)",
    "amr_alleles": "Resistance target alleles",
    # shared by more than one module, so resolved per module below
    "cps": "Capsule locus",
    "mlst": "MLST",
    # Salmonella
    "salmonella_serovar": "Serovar (antigenic formula)",
    "salmonella_subspecies": "Subspecies",
    "salmonella_h2_locus": "H2 (fljB) locus",
    "salmonella_spi": "Pathogenicity islands",
    "salmonella_mlst": "MLST / eBurst group",
    "salmonella_virulence_plasmid": "Virulence plasmid",
    "salmonella_typhoidal": "Typhoidal markers",
    # Escherichia coli / Shigella
    "phylogroup": "Phylogroup (Clermont)",
    "fimH": "fimH allele (FimTyper)",
    "ST131_subclone": "ST131 subclone",
    "stx": "Shiga toxin subtype",
    "pathotype": "Pathotype",
    "intimin": "Intimin (eae) subtype",
    "ExPEC_VF": "ExPEC / UPEC virulence factors",
    "shigella_eiec": "Shigella / EIEC cluster",
}

#: Where two modules use the same scheme identifier for different things.
#: ``cps`` is the capsular serotype of a group B streptococcus and the capsule
#: serogroup of a meningococcus, and the two must not share a heading.
SCHEME_TITLES_BY_MODULE: dict[tuple[str, str], str] = {
    ("sagalactiae", "cps"): "Capsular serotype",
    ("nmeningitidis", "cps"): "Capsule serogroup",
    ("hinfluenzae", "mlst"): "MLST (H. influenzae scheme)",
    ("nmeningitidis", "mlst"): "MLST and clonal complex",
}

#: The calls that answer "what is this isolate?" for each organism, with the
#: short prefix a one-line summary uses. The order is the order a report should
#: read them in, and the whole point is brevity: this goes in a table cell and
#: in the title of a sheet, where there is room for one line and not for
#: twenty, so the full set of calls stays in the Typing sheet beneath it.
HEADLINE_SCHEMES: dict[str, tuple[tuple[str, str], ...]] = {
    "saureus": (("spa", "spa"), ("clonal_complex", ""), ("agr", "agr")),
    "spyogenes": (("emm", ""),),
    "sagalactiae": (("cps", "serotype"), ("sag_alp", "Alp")),
    "spneumoniae": (("cps_serotype", "serotype"), ("pbp_type", "PBP")),
    "hinfluenzae": (("cap", "capsule"), ("pbp3", "PBP3")),
    "nmeningitidis": (("cps", "serogroup"), ("finetype", "")),
    "salmonella": (("salmonella_serovar", ""),),
    "ecoli": (("pathotype", ""), ("phylogroup", "phylogroup")),
}

#: Which files on disk belong to each declared database, relative to
#: ``db_dir/<KEY>/``. The modules do not share a layout -- some give each
#: database a directory, some write files side by side, and *S. pneumoniae*
#: puts two declared databases in one directory -- so a name is not a path and
#: cannot be guessed. A trailing ``/`` is a directory, and a trailing ``*``
#: picks up the BLAST index files written beside a FASTA.
#:
#: This is what lets a packager ship one of a module's databases and download
#: another: the Windows build bundles the *S. aureus* marker panel and fetches
#: the Ridom *spa* database on first launch, and it needs to know which files
#: are which.
DB_LAYOUT: dict[str, dict[str, tuple[str, ...]]] = {
    "ecoli": {
        "fimtyper": ("fimtyper/",),
        "stxtyper": ("stxtyper/",),
        "virulencefinder": ("virulencefinder/",),
        "shigeifinder": ("shigeifinder/",),
    },
    "hinfluenzae": {
        "cap": ("cap/",),
        "ftsI": ("ftsI/",),
        "bla": ("bla/",),
        "mlst": ("mlst/",),
    },
    "nmeningitidis": {
        "capsule": ("capsule/",),
        "peptides": ("peptides/",),
        "mlst": ("mlst/",),
        "amr": ("amr/",),
        "antigens": ("antigens/",),
    },
    "sagalactiae": {
        "cps": ("cps.fna*", "cps.LICENSE", "cps.version"),
        "panel": ("panel.faa*", "panel.manifest.json", "pbp2x.faa*"),
        "surface_targets": ("surface_targets.fna*", "surface_targets.version"),
        "alp_alleles": ("alp_alleles.fna*",),
    },
    "salmonella": {
        "antigens": ("antigens/antigens.fna*", "antigens/antigens.sources.tsv"),
        "serovars": ("antigens/serovars.csv",),
        "spi": ("spi/",),
        "markers": ("markers/",),
        "mlst": ("mlst/",),
    },
    "saureus": {
        "spa": ("spa/",),
        "markers": ("markers/",),
        "mlst": ("mlst/",),
    },
    "spneumoniae": {
        "cps": ("cps/cps_loci.fna*",),
        "cdc-serotype-probes": ("cps/sero_probes.fna*",),
        "cdc-pbp": ("pbp/",),
        "cdc-resistance-probes": ("res/",),
    },
    "spyogenes": {
        "emm": ("emm.tfa*", "emm.version", "emm_clusters.tsv"),
        "panel": ("panel.faa*", "panel.manifest.json"),
    },
}

#: Provenance files a module writes beside its databases. They are small, they
#: belong with whatever was installed, and they are what a user reads to find
#: out where their data came from, so they travel with any subset.
DB_PROVENANCE: dict[str, tuple[str, ...]] = {
    "ecoli": ("provenance.tsv",),
    "hinfluenzae": ("version.json",),
    "nmeningitidis": ("version.json",),
    "sagalactiae": ("cps.LICENSE",),
    "salmonella": ("provenance.json",),
    "saureus": ("manifest.json",),
    "spneumoniae": ("SOURCES.tsv",),
    "spyogenes": (),
}

#: Files whose presence makes a module's ``fetch()`` return immediately,
#: treating the install as complete. They are completion stamps, not data, and
#: a packager that ships only part of a module's databases must leave them out:
#: with the stamp present the first-launch download sees "already installed"
#: and never fetches the rest. Only these two modules have one -- the other six
#: guard per file, so a missing database is still fetched.
DB_STAMPS: dict[str, tuple[str, ...]] = {
    "hinfluenzae": ("version.json",),
    "nmeningitidis": ("version.json",),
}

#: Download caches a module keeps so a re-fetch does not re-download. They are
#: the largest thing in two of these trees (12 MB of RefSeq protein FASTA for
#: *S. pyogenes*) and nothing reads them at run time, so a packaged copy leaves
#: them out.
DB_CACHES: dict[str, tuple[str, ...]] = {
    "sagalactiae": ("refseq/",),
    "spyogenes": ("refseq/",),
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

    def layout(self) -> dict[str, tuple[str, ...]]:
        """This module's :data:`DB_LAYOUT` entry."""
        return DB_LAYOUT.get(self.key, {})

    def files_for(self, name: str, db_dir: Path | str) -> list[Path]:
        """Every file on disk that belongs to the declared database *name*."""
        root = self.root(db_dir)
        out: list[Path] = []
        for pattern in self.layout().get(name, (f"{name}/",)):
            if pattern.endswith("/"):
                base = root / pattern.rstrip("/")
                if base.is_dir():
                    out.extend(sorted(p for p in base.rglob("*") if p.is_file()))
            elif pattern.endswith("*"):
                parent = (root / pattern[:-1]).parent
                stem = Path(pattern[:-1]).name
                if parent.is_dir():
                    out.extend(sorted(p for p in parent.iterdir()
                                      if p.is_file() and p.name.startswith(stem)))
            else:
                path = root / pattern
                if path.is_file():
                    out.append(path)
        return out

    def provenance_files(self, db_dir: Path | str, complete: bool = True) -> list[Path]:
        """The small files recording where this module's data came from.

        With *complete* false the completion stamps are left out, because a
        partial copy that carries one would never be completed: see
        :data:`DB_STAMPS`.
        """
        root = self.root(db_dir)
        names = list(DB_PROVENANCE.get(self.key, ()))
        if not complete:
            names = [n for n in names if n not in DB_STAMPS.get(self.key, ())]
        return [root / name for name in names if (root / name).is_file()]

    def stamp_files(self, db_dir: Path | str) -> list[Path]:
        root = self.root(db_dir)
        return [root / name for name in DB_STAMPS.get(self.key, ())
                if (root / name).is_file()]

    def cache_dirs(self, db_dir: Path | str) -> list[Path]:
        root = self.root(db_dir)
        return [root / name.rstrip("/") for name in DB_CACHES.get(self.key, ())
                if (root / name.rstrip("/")).exists()]

    def missing_databases(self, db_dir: Path | str) -> list[str]:
        """The declared databases with no files under *db_dir*."""
        return [db.name for db in self.databases
                if not self.files_for(db.name, db_dir)]

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
        return any(self.files_for(db.name, db_dir) for db in self.databases)

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
        return relabel(list(results), organism=self.title, key=self.key)


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
def relabel(results: Iterable[TypingResult], organism: str = "",
            key: str = "") -> list[TypingResult]:
    """Give every result a readable heading, and correct the misleading ones.

    Presentation only: no call, allele or note is replaced, only added to.
    *organism* names the module the calls came from, so a report can group a
    genome's calls by it, and *key* is that module's key, which decides the
    heading where two modules share a scheme identifier.
    """
    out: list[TypingResult] = []
    for result in results:
        title = SCHEME_TITLES_BY_MODULE.get((key, result.scheme))
        if title is None:
            title = SCHEME_TITLES.get(result.scheme)
        if title is None:
            title = result.scheme.replace("__", " ").replace("_", " ")
        note = result.note
        if result.scheme == "emm" and _NON_PYOGENES_EMM.match(str(result.call or "")):
            title = "M-protein gene type (group C/G, S. dysgalactiae)"
            extra = _NON_PYOGENES_EMM_NOTE.format(call=result.call)
            note = f"{note}; {extra}" if note else extra
        organism_name = organism or result.organism
        if (title == result.title and note == result.note
                and organism_name == result.organism):
            out.append(result)
            continue
        out.append(TypingResult(
            scheme=result.scheme, call=result.call, lineage=result.lineage,
            alleles=dict(result.alleles), loci_found=result.loci_found,
            loci_total=result.loci_total, note=note, title=title,
            organism=organism_name))
    return out


def plain_fasta(assembly: Path, workdir: Path) -> Path:
    """A path BLAST+ can read: a gzipped assembly is expanded into *workdir*.

    BLAST+ cannot read a compressed FASTA, and the modules hand the assembly
    straight to ``makeblastdb`` and ``blastn``. Hydra accepts ``.gz`` input --
    which is how most assemblies arrive -- so without this every organism
    module fails on a gzipped genome with "makeblastdb failed". Expanded once
    here rather than in eight modules.

    The magic bytes decide, not the suffix: a gzipped file named ``.fna`` is
    common, and a plain file named ``.gz`` is not unknown either.
    """
    assembly = Path(assembly)
    try:
        with open(assembly, "rb") as handle:
            compressed = handle.read(2) == b"\x1f\x8b"
    except OSError as exc:
        raise HydraError(f"cannot read {assembly}: {exc}") from exc
    if not compressed:
        return assembly
    workdir.mkdir(parents=True, exist_ok=True)
    plain = workdir / (assembly.name[:-3] if assembly.suffix == ".gz"
                       else assembly.name + ".fna")
    if plain.exists() and plain.stat().st_size:
        return plain
    LOG.debug("organism typing: expanding %s to %s", assembly.name, plain)
    with gzip.open(assembly, "rb") as src, open(plain, "wb") as out:
        shutil.copyfileobj(src, out)
    if not plain.stat().st_size:
        raise HydraError(f"{assembly} expanded to an empty file")
    return plain


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
    chosen = select(species, mlst, only=only)
    if not chosen:
        return []
    temporary: str | None = None
    if workdir is None:
        manager = tempdir(prefix="hydra.organisms.")
        temporary = manager.__enter__()
        work = Path(temporary)
    else:
        work = Path(workdir)
        work.mkdir(parents=True, exist_ok=True)
    try:
        return _type_with(chosen, assembly, db_dir, work, mlst, threads, warn)
    finally:
        if temporary is not None:
            manager.__exit__(None, None, None)


def _type_with(chosen: list[OrganismModule], assembly: Path | str, db_dir: Path,
               work: Path, mlst: MlstCall | None, threads: int,
               warn: list[str] | None) -> list[TypingResult]:
    out: list[TypingResult] = []
    # Expanded on first use, not up front: a module whose database is missing
    # is skipped before the assembly is touched at all.
    plain: list[Path] = []

    def readable() -> Path:
        if not plain:
            plain.append(plain_fasta(Path(assembly), work))
        return plain[0]

    for entry in chosen:
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
            results = entry.run(readable(), db_dir, threads=threads,
                                workdir=work / entry.key, mlst=mlst)
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


def headline(results: Iterable[TypingResult]) -> str:
    """One short line naming what the organism-specific typing concluded.

    The headline calls of every module that reported, in the order
    :data:`HEADLINE_SCHEMES` gives them. A call that could not be made is kept
    and shown as such -- "serovar not determined" is a result, and dropping it
    would leave the line reading as though the isolate had been fully typed.
    """
    by_scheme: dict[str, TypingResult] = {}
    keys: list[str] = []
    # Which module made the call, not which modules claim the scheme name:
    # ``cps`` belongs to both the group B streptococcus and the meningococcus,
    # and matching by name alone put a meningococcus under both headings.
    by_title = {entry.title: key for key, entry in modules().items()}
    for result in results:
        if not result.organism:
            continue
        by_scheme.setdefault(result.scheme, result)
        key = by_title.get(result.organism)
        if key in HEADLINE_SCHEMES and key not in keys:
            keys.append(key)
    parts: list[str] = []
    for key in keys:
        for scheme, prefix in HEADLINE_SCHEMES[key]:
            found = by_scheme.get(scheme)
            if found is None:
                continue
            call = str(found.call or "-")
            if prefix:
                parts.append(f"{prefix} {call}")
            elif call.startswith("not determined"):
                # Without a prefix a bare "not determined" says nothing about
                # which call is missing, and this line is often all a reader
                # sees, so the scheme names itself.
                label = (SCHEME_TITLES_BY_MODULE.get((key, scheme))
                         or SCHEME_TITLES.get(scheme, scheme))
                parts.append(f"{label}: {call}")
            else:
                parts.append(call)
    return " \u00b7 ".join(parts)
