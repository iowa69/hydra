"""Catalogue of the reference databases Hydra knows how to install and use."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..utils import HydraError

#: Element types Hydra assigns to hits.
AMR = "AMR"
VIRULENCE = "VIRULENCE"
STRESS = "STRESS"
PLASMID = "PLASMID"
SEROTYPE = "SEROTYPE"


@dataclass(frozen=True)
class DbSpec:
    """Static description of one reference database."""

    name: str
    kind: str  # nucl | prot | mutation | mlst | typing | sketch
    title: str
    element_type: str = AMR
    #: Where an `import` can find it inside a conda environment.
    conda_env_hint: tuple[str, ...] = ()
    conda_rel_path: str = ""
    #: Canonical upstream source, shown by `hydra db info` and used by `db download`.
    url: str = ""
    citation: str = ""
    licence: str = ""
    #: Default identity/coverage overrides for this database.
    default_identity: float | None = None
    default_coverage: float | None = None
    #: Databases in the "standard" preset.
    standard: bool = False
    aliases: tuple[str, ...] = field(default_factory=tuple)
    notes: str = ""


DATABASES: dict[str, DbSpec] = {}


def _register(spec: DbSpec) -> None:
    DATABASES[spec.name] = spec


# --------------------------------------------------------------------------
# Acquired-gene nucleotide databases (abricate-compatible sources)
# --------------------------------------------------------------------------
_register(DbSpec(
    name="ncbi", kind="nucl", title="NCBI AMRFinderPlus reference gene catalog (nucleotide)",
    element_type=AMR, conda_env_hint=("abricate",), conda_rel_path="db/ncbi",
    url="https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest",
    citation="Feldgarden et al. 2021, Sci Rep 11:12728", licence="Public domain",
    standard=True, notes="Curated, low false-positive rate. The default database.",
))
_register(DbSpec(
    name="card", kind="nucl", title="CARD - Comprehensive Antibiotic Resistance Database",
    element_type=AMR, conda_env_hint=("abricate",), conda_rel_path="db/card",
    url="https://card.mcmaster.ca/latest/data",
    citation="Alcock et al. 2023, Nucleic Acids Res 51:D690", licence="CARD academic licence",
    standard=True, notes="Broad ontology coverage including efflux and regulator genes.",
))
_register(DbSpec(
    name="resfinder", kind="nucl", title="ResFinder acquired resistance genes",
    element_type=AMR, conda_env_hint=("abricate",), conda_rel_path="db/resfinder",
    url="https://bitbucket.org/genomicepidemiology/resfinder_db",
    citation="Bortolaia et al. 2020, J Antimicrob Chemother 75:3491", licence="Apache-2.0",
    standard=True,
))
_register(DbSpec(
    name="argannot", kind="nucl", title="ARG-ANNOT antibiotic resistance genes",
    element_type=AMR, conda_env_hint=("abricate",), conda_rel_path="db/argannot",
    url="https://www.mediterranee-infection.com/acces-ressources/base-de-donnees/arg-annot-2/",
    citation="Gupta et al. 2014, Antimicrob Agents Chemother 58:212",
))
_register(DbSpec(
    name="megares", kind="nucl", title="MEGARes resistance, biocide and metal genes",
    element_type=AMR, conda_env_hint=("abricate",), conda_rel_path="db/megares",
    url="https://www.meglab.org/megares/",
    citation="Bonin et al. 2023, Nucleic Acids Res 51:D744",
    notes="Includes biocide and metal resistance determinants.",
))
_register(DbSpec(
    name="vfdb", kind="nucl", title="VFDB core virulence factors",
    element_type=VIRULENCE, conda_env_hint=("abricate",), conda_rel_path="db/vfdb",
    url="http://www.mgc.ac.cn/VFs/download.htm",
    citation="Liu et al. 2022, Nucleic Acids Res 50:D912", standard=True,
))
_register(DbSpec(
    name="vfdb_full", kind="nucl", title="VFDB full dataset (setB, all predicted VFs)",
    element_type=VIRULENCE, conda_env_hint=("abricate",), conda_rel_path="db/vfdb_full",
    url="http://www.mgc.ac.cn/VFs/download.htm",
    citation="Liu et al. 2022, Nucleic Acids Res 50:D912",
    notes="Large (~70 MB); use when broad virulence screening is required.",
))
_register(DbSpec(
    name="ecoli_vf", kind="nucl", title="Escherichia coli virulence factors",
    element_type=VIRULENCE, conda_env_hint=("abricate",), conda_rel_path="db/ecoli_vf",
    url="https://github.com/phac-nml/ecoli_vf",
))
_register(DbSpec(
    name="plasmidfinder", kind="nucl", title="PlasmidFinder replicon typing",
    element_type=PLASMID, conda_env_hint=("abricate",), conda_rel_path="db/plasmidfinder",
    url="https://bitbucket.org/genomicepidemiology/plasmidfinder_db",
    citation="Carattoli et al. 2014, Antimicrob Agents Chemother 58:3895",
    default_identity=95.0, default_coverage=60.0,
))
_register(DbSpec(
    name="ecoh", kind="nucl", title="E. coli O- and H-antigen serotyping loci",
    element_type=SEROTYPE, conda_env_hint=("abricate",), conda_rel_path="db/ecoh",
    url="https://github.com/katholt/srst2",
    default_identity=90.0, default_coverage=80.0,
))

# --------------------------------------------------------------------------
# Protein + point-mutation reference (AMRFinderPlus data directory)
# --------------------------------------------------------------------------
_register(DbSpec(
    name="protein", kind="prot", title="NCBI protein & point-mutation reference (AMRFinderPlus data)",
    element_type=AMR, conda_env_hint=("amrfinder", "amrfinderplus", "ncbi-amrfinderplus", "linezolid-amr", "klebo"),
    conda_rel_path="share/amrfinderplus/data/latest",
    url="https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest",
    citation="Feldgarden et al. 2021, Sci Rep 11:12728", licence="Public domain",
    standard=True, aliases=("amrfinderplus", "amrfinder", "afp", "point"),
    notes="Supplies translated-search AMR calls, --plus stress/virulence classes, "
          "and the organism-specific protein and DNA point-mutation catalogues.",
))

# --------------------------------------------------------------------------
# Typing resources
# --------------------------------------------------------------------------
_register(DbSpec(
    name="pubmlst", kind="mlst", title="PubMLST 7-locus MLST schemes",
    element_type="TYPING", conda_env_hint=("mlst",), conda_rel_path="db",
    url="https://rest.pubmlst.org/db",
    citation="Jolley, Bray & Maiden 2018, Wellcome Open Res 3:124",
    standard=True, aliases=("mlst",),
))
_register(DbSpec(
    name="lineage", kind="typing", title="Lineage / sublineage typing loci (Kleborate-derived and extensions)",
    element_type="TYPING", conda_env_hint=("klebo", "kleborate"),
    conda_rel_path="lib/*/site-packages/kleborate/modules",
    url="https://github.com/klebgenomics/Kleborate",
    citation="Lam et al. 2021, Nat Commun 12:4188",
    standard=True, aliases=("kleborate", "typing"),
    notes="Yersiniabactin/colibactin/aerobactin/salmochelin/rmp sublineage schemes plus "
          "Hydra's generic per-species virulence and resistance scoring rules.",
))
_register(DbSpec(
    name="species", kind="sketch", title="Species identification sketches",
    element_type="TYPING", conda_env_hint=("klebo", "kleborate"),
    conda_rel_path="lib/*/site-packages/kleborate/modules/enterobacterales__species/data",
    url="https://github.com/klebgenomics/Kleborate",
    citation="Ondov et al. 2016, Genome Biol 17:132",
    aliases=("mash",),
    notes="Optional Mash sketches; Hydra also identifies species from MLST and marker loci.",
))


#: Convenience groups usable anywhere a database name is accepted.
DB_GROUPS: dict[str, tuple[str, ...]] = {
    "all": tuple(n for n, s in DATABASES.items() if s.kind in ("nucl", "prot")),
    "standard": tuple(n for n, s in DATABASES.items() if s.standard and s.kind in ("nucl", "prot")),
    "amr": ("ncbi", "card", "resfinder", "argannot", "megares"),
    "virulence": ("vfdb", "ecoli_vf"),
    "nucl": tuple(n for n, s in DATABASES.items() if s.kind == "nucl"),
    "core": ("ncbi", "vfdb"),
}


def resolve_names(names: list[str]) -> list[str]:
    """Expand group names and aliases into a de-duplicated database name list."""
    alias_map = {}
    for spec in DATABASES.values():
        for alias in spec.aliases:
            alias_map[alias] = spec.name
    out: list[str] = []
    for raw in names:
        for token in str(raw).split(","):
            token = token.strip()
            if not token:
                continue
            key = token.lower()
            if key in DB_GROUPS:
                candidates = list(DB_GROUPS[key])
            elif key in alias_map:
                candidates = [alias_map[key]]
            elif key in DATABASES:
                candidates = [key]
            else:
                known = ", ".join(sorted(set(list(DATABASES) + list(DB_GROUPS))))
                raise HydraError(f"unknown database '{token}'. Known names and groups: {known}")
            for candidate in candidates:
                if candidate not in out:
                    out.append(candidate)
    return out


#: What "hydra db download" installs when no names are given: enough for
#: --preset standard to run, MLST included.
DEFAULT_DOWNLOADS = ("protein", "ncbi", "card", "vfdb", "pubmlst", "lineage", "species")

#: Roughly how long each download takes, for the message printed before starting.
#: PubMLST is a thousand small files and is the only one worth warning about.
SLOW_DOWNLOADS = {"pubmlst": "takes several minutes: about 1200 files",
                  "lineage": "pulls a 70 MB release archive, shared with species"}


def protein_dir(db_root: Path | str) -> Path:
    """Where the protein reference is installed under ``db_root``.

    Stores built before the database was renamed hold it under the old
    directory name, and renaming the database must not orphan them.
    """
    root = Path(db_root)
    legacy = root / "prot" / "amrfinderplus"
    current = root / "prot" / "protein"
    return legacy if legacy.exists() and not current.exists() else current


def spec_for(name: str) -> DbSpec:
    if name not in DATABASES:
        raise HydraError(f"unknown database '{name}'")
    return DATABASES[name]
