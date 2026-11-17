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
    #: A nucleotide database that is not an acquired-gene database. SCCmec
    #: references are whole cassettes: they belong to `hydra run`'s typing step,
    #: not to `-d all`, where they would report a 30 kb element as a gene hit.
    typing_only: bool = False
    aliases: tuple[str, ...] = field(default_factory=tuple)
    notes: str = ""


#: The PubMLST terms, quoted where they decide what may be redistributed.
#:
#: PubMLST publishes no open licence. Its Terms and Conditions
#: (https://pubmlst.org/terms-conditions, page dated 2026-05-13, retrieved
#: 2026-11-17) split the data by submission date, and the deciding sentence is
#: quoted here in full because a paraphrase would not do: everything Hydra
#: ships or redistributes rests on it.
#:
#: The right depends entirely on the download being anonymous. The
#: unauthenticated REST API answers every request with
#: ``"Please note that you are currently restricted to accessing data that was
#: submitted on or prior to 2024-12-31. Please authenticate to access the full
#: dataset."``, and an anonymous ``?added_after=2025-01-01`` query returns zero
#: records, so an anonymous download can only return pre-2025 data. Sending an
#: API key or any other credential would pull post-2024 data and with it
#: Section 4's clause forbidding redistribution -- so no fetcher in Hydra sends
#: one, and none ever should.
PUBMLST_TERMS = (
    "PubMLST Terms and Conditions (https://pubmlst.org/terms-conditions, page dated "
    "2026-05-13, retrieved 2026-11-17); not a named licence but a University of Oxford "
    "terms-of-use document. The deciding sentence: \"All data submitted on or before "
    "31 December 2024 may continue to be downloaded, used, and redistributed without "
    "restriction, subject only to proper scientific citation and acknowledgement.\" The "
    "terms add that \"The restrictions described in Sections 4 and 5 do not apply to "
    "pre-2025 data\", and \"Anonymous (non-logged in) users may view only pre 2025 "
    "data\". This snapshot is therefore fetched anonymously and contains only pre-2025 "
    "data; it must never be fetched with an API key or any other credential. Required "
    "acknowledgement: \"This publication made use of the PubMLST website "
    "(https://pubmlst.org/) sited at the University of Oxford. The development of that "
    "website was funded by the Wellcome Trust.\" Cite Jolley, Bray & Maiden 2018, "
    "Wellcome Open Res 3:124."
)

#: The acknowledgement PubMLST asks to appear wherever its data is used.
PUBMLST_ACKNOWLEDGEMENT = (
    "This publication made use of the PubMLST website (https://pubmlst.org/) sited at "
    "the University of Oxford. The development of that website was funded by the "
    "Wellcome Trust."
)

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
    licence=PUBMLST_TERMS,
    standard=True, aliases=("mlst",),
    notes="Downloaded anonymously, which is what makes the snapshot "
          "redistributable: see PUBMLST_TERMS and db/fetch.py's "
          "verify_pubmlst_is_pre_2025().",
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
    name="sccmec", kind="nucl", title="SCCmec cassette reference elements",
    element_type="TYPING", typing_only=True,
    url="https://bitbucket.org/genomicepidemiology/sccmecfinder_db/get/HEAD.tar.gz",
    citation="Kaya et al. 2018, mSphere 3:e00612-17", licence="Apache 2.0",
    aliases=("sccmec_finder",),
    notes="Whole-cassette references for SCCmec types I-XIII. Used by the typing "
          "step for staphylococci, not by gene screening.",
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


# --------------------------------------------------------------------------
# Organism-specific typing databases, one per module under
# ``hydra_amr.typing.organisms``. Registered from the modules themselves so the
# source, version, licence and redistribution status shown by ``hydra db info``
# are the ones the module author recorded, never a second copy of them.
# --------------------------------------------------------------------------
def _register_organism_databases() -> None:
    try:
        from ..typing.organisms import registry as organisms
    except Exception:  # pragma: no cover - never let this break the db CLI
        return
    for key, entry in organisms.modules().items():
        parts = entry.databases
        sources = "; ".join(f"{db.name}: {db.source}" for db in parts)
        licences = " | ".join(f"{db.name}: {db.licence}" for db in parts)
        redistribution = ", ".join(f"{db.name}={db.redistribution}" for db in parts)
        citations = "; ".join(sorted({ref.key for ref in entry.references}))
        _register(DbSpec(
            name=f"typing-{key}",
            kind="organism",
            title=f"{entry.title} organism typing data",
            element_type="TYPING",
            url=sources,
            citation=citations,
            licence=licences,
            notes=(f"Installed by the module's own fetch(). {len(parts)} "
                   f"database(s): {', '.join(db.name for db in parts)}. "
                   f"Redistribution: {redistribution}. Triggered by the species "
                   f"call ({', '.join(entry.species) or '-'}) or the MLST scheme "
                   f"({', '.join(entry.mlst_schemes) or '-'})."),
            aliases=(key,),
        ))


_register_organism_databases()


#: The organism typing database names, for the CLI and the groups below.
ORGANISM_DATABASES: tuple[str, ...] = tuple(
    name for name, spec in DATABASES.items() if spec.kind == "organism")


#: Convenience groups usable anywhere a database name is accepted.
DB_GROUPS: dict[str, tuple[str, ...]] = {
    "all": tuple(n for n, s in DATABASES.items() if s.kind in ("nucl", "prot") and not s.typing_only),
    "standard": tuple(n for n, s in DATABASES.items() if s.standard and s.kind in ("nucl", "prot") and not s.typing_only),
    "amr": ("ncbi", "card", "resfinder", "argannot", "megares"),
    "virulence": ("vfdb", "ecoli_vf"),
    "nucl": tuple(n for n, s in DATABASES.items() if s.kind == "nucl" and not s.typing_only),
    "core": ("ncbi", "vfdb"),
    "organisms": ORGANISM_DATABASES,
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
