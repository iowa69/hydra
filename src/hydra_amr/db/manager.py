"""Installation, normalisation and lookup of Hydra's reference databases.

Every acquired-gene database is normalised into the same on-disk shape so the
screening engines never have to care where the sequences came from:

    <db_dir>/nucl/<name>/sequences.fna   short unique ids (``<name>_000123``)
    <db_dir>/nucl/<name>/meta.tsv        id -> gene/accession/product/class/...
    <db_dir>/nucl/<name>/sequences.fna.n*   BLAST index

Drug-class annotation is transferred onto every nucleotide database from the
AMRFinderPlus family table when it is available, which gives databases such as
ARG-ANNOT or CARD a consistent ``class``/``subclass`` column they do not ship.
"""

from __future__ import annotations

import glob
import json
import os
import re
import shutil
import time
from dataclasses import dataclass
from pathlib import Path

from ..seqio import read_fasta, write_fasta
from ..utils import LOG, HydraError, human_bytes, require, run
from .registry import (AMR, DATABASES, PLASMID, SEROTYPE, STRESS, VIRULENCE, DbSpec,
                       protein_dir, spec_for)

MANIFEST_NAME = "manifest.json"

META_COLUMNS = ("seqid", "gene", "accession", "product", "class", "subclass",
                "element_type", "element_subtype", "fam_id")

# Keyword -> drug class, used when a database ships no structured annotation.
_CLASS_KEYWORDS: tuple[tuple[str, str], ...] = (
    ("carbapenem", "BETA-LACTAM"), ("beta-lactam", "BETA-LACTAM"), ("metallo-beta", "BETA-LACTAM"),
    ("cephalospor", "BETA-LACTAM"), ("penicillin", "BETA-LACTAM"), ("oxacillin", "BETA-LACTAM"),
    ("ampicillin", "BETA-LACTAM"), ("cephamycin", "BETA-LACTAM"), ("monobactam", "BETA-LACTAM"),
    ("aminoglycoside", "AMINOGLYCOSIDE"), ("streptomycin", "AMINOGLYCOSIDE"),
    ("gentamicin", "AMINOGLYCOSIDE"), ("kanamycin", "AMINOGLYCOSIDE"), ("tobramycin", "AMINOGLYCOSIDE"),
    ("amikacin", "AMINOGLYCOSIDE"), ("apramycin", "AMINOGLYCOSIDE"), ("neomycin", "AMINOGLYCOSIDE"),
    ("macrolide", "MACROLIDE"), ("erythromycin", "MACROLIDE"), ("azithromycin", "MACROLIDE"),
    ("lincosamide", "LINCOSAMIDE"), ("clindamycin", "LINCOSAMIDE"),
    ("streptogramin", "STREPTOGRAMIN"), ("virginiamycin", "STREPTOGRAMIN"),
    ("oxazolidinone", "OXAZOLIDINONE"), ("linezolid", "OXAZOLIDINONE"),
    ("tetracycline", "TETRACYCLINE"), ("tigecycline", "TETRACYCLINE"), ("doxycycline", "TETRACYCLINE"),
    ("quinolone", "QUINOLONE"), ("fluoroquinolone", "QUINOLONE"), ("ciprofloxacin", "QUINOLONE"),
    ("nalidixic", "QUINOLONE"),
    ("phenicol", "PHENICOL"), ("chloramphenicol", "PHENICOL"), ("florfenicol", "PHENICOL"),
    ("sulphonamide", "SULFONAMIDE"), ("sulfonamide", "SULFONAMIDE"), ("sulfamethoxazole", "SULFONAMIDE"),
    ("trimethoprim", "TRIMETHOPRIM"),
    ("glycopeptide", "GLYCOPEPTIDE"), ("vancomycin", "GLYCOPEPTIDE"), ("teicoplanin", "GLYCOPEPTIDE"),
    ("colistin", "COLISTIN"), ("polymyxin", "COLISTIN"),
    ("rifampin", "RIFAMYCIN"), ("rifamycin", "RIFAMYCIN"), ("rifampicin", "RIFAMYCIN"),
    ("fosfomycin", "FOSFOMYCIN"), ("bleomycin", "BLEOMYCIN"), ("nitroimidazole", "NITROIMIDAZOLE"),
    ("mupirocin", "MUPIROCIN"), ("fusidic", "FUSIDIC-ACID"), ("bacitracin", "BACITRACIN"),
    ("pleuromutilin", "PLEUROMUTILIN"), ("lipopeptide", "LIPOPEPTIDE"), ("daptomycin", "LIPOPEPTIDE"),
    ("quaternary ammonium", "QUATERNARY-AMMONIUM"), ("biocide", "BIOCIDE"),
    ("mercury", "METAL"), ("copper", "METAL"), ("arsenic", "METAL"), ("silver", "METAL"),
    ("tellurite", "METAL"), ("nickel", "METAL"), ("zinc", "METAL"),
    ("efflux", "EFFLUX"), ("multidrug", "MULTIDRUG"),
)


@dataclass
class DbHandle:
    """A ready-to-use installed database."""

    name: str
    kind: str
    path: Path
    element_type: str
    fasta: Path | None = None
    meta: dict[str, dict] | None = None
    version: str = ""
    n_sequences: int = 0

    def meta_for(self, seqid: str) -> dict:
        if self.meta is None:
            return {}
        return self.meta.get(seqid, {})


def _normalise_gene(name: str) -> str:
    """Reduce a gene label to a comparable key across databases.

    Each database decorates its gene names differently:
    ``blaKPC-2_1_AY034847`` (ResFinder, with variant and accession),
    ``ARR-2_1`` (ResFinder, variant only) and ``AAC(2')-IIa`` (CARD, cased)
    all have to collapse onto the AMRFinderPlus symbol so curated drug-class
    and gene-family annotation can be transferred onto every database.
    """
    name = re.sub(r"\s+", "", name.strip())
    # ARG-ANNOT prefixes every gene with its drug-class code: (AGly)AAC(6')-Isa
    name = re.sub(rf"^\((?:{'|'.join(_ARGANNOT_CLASSES)})\)", "", name, flags=re.IGNORECASE)
    parts = name.split("_")
    # ResFinder: trailing accession, preceded by a variant number.
    if len(parts) >= 3 and re.fullmatch(r"[A-Z]{1,2}\d{5,}(\.\d+)?", parts[-1]):
        parts = parts[:-2] or parts[:1]
    # ResFinder: bare trailing variant number (ARR-2_1). Only strip it when what
    # comes before is not itself purely numeric, so 1567214_ble is left alone.
    if len(parts) >= 2 and parts[-1].isdigit() and not parts[-2].isdigit():
        parts = parts[:-1]
    name = "_".join(parts).lower()
    # Note: doubled primes are NOT collapsed. aph(3'')-Ia (streptomycin) and
    # aph(3')-Ia (kanamycin) are different genes, and merging their keys would
    # give one of them the other's drug class and hide it in read mode.
    name = re.sub(r"^(bla|aac|aph|ant|aad)_", r"\1", name)
    return name


#: ARG-ANNOT encodes the drug class as a prefix on every gene name.
_ARGANNOT_CLASSES = {
    "agly": "AMINOGLYCOSIDE", "bla": "BETA-LACTAM", "col": "COLISTIN",
    "colistin": "COLISTIN", "fcyn": "FOSFOMYCIN", "fcd": "FOSFOMYCIN",
    "fos": "FOSFOMYCIN", "flq": "QUINOLONE", "gly": "GLYCOPEPTIDE",
    "mls": "MACROLIDE", "phe": "PHENICOL", "rif": "RIFAMYCIN",
    "sul": "SULFONAMIDE", "tet": "TETRACYCLINE", "tmt": "TRIMETHOPRIM",
    "ntmdz": "NITROIMIDAZOLE", "bcl": "BACITRACIN", "oxzln": "OXAZOLIDINONE",
    "mupirocin": "MUPIROCIN", "fus": "FUSIDIC-ACID",
    "tetracenomycinc": "TETRACYCLINE", "trim": "TRIMETHOPRIM",
}

#: Beta-lactamase families that CARD and MEGARes name without the ``bla`` prefix
#: AMRFinderPlus uses (``OXA-48`` vs ``blaOXA-48``).
_BLA_FAMILIES = (
    "oxa", "tem", "shv", "ctx-m", "cmy", "imp", "vim", "ndm", "kpc", "ges", "per",
    "veb", "act", "mir", "dha", "fox", "adc", "carb", "pse", "lap", "sfo", "tla",
    "bel", "acc", "cfx", "mox", "oxy", "sme", "nmc", "ime", "spm", "gim", "sim",
    "aim", "dim", "bic", "tmb", "fri", "lcr", "slb", "rob", "hera", "ec", "cmh",
)


def _lookup_keys(gene: str) -> list[str]:
    """Keys to try when transferring annotation onto a gene name.

    CARD and MEGARes call beta-lactamases ``OXA-48``/``NDM-1`` where
    AMRFinderPlus calls them ``blaOXA-48``/``blaNDM-1``; without the second key
    the curated subclass (which is what separates a carbapenemase from a plain
    penicillinase) is lost for well over half of both databases.
    """
    key = _normalise_gene(gene)
    keys = [key]
    if not key.startswith("bla"):
        family = key.split("-")[0]
        if family in _BLA_FAMILIES or any(
                key.startswith(f"{fam}-") for fam in _BLA_FAMILIES):
            keys.append("bla" + key)
    return keys


def _infer_class(text: str) -> str:
    prefix = re.match(r"^\(([A-Za-z]+)\)", text.strip())
    if prefix:
        klass = _ARGANNOT_CLASSES.get(prefix.group(1).lower())
        if klass:
            return klass
    low = text.lower()
    for keyword, klass in _CLASS_KEYWORDS:
        if keyword in low:
            return klass
    return ""


def _parse_abricate_header(header: str) -> dict:
    """Parse ``db~~~gene~~~accession~~~product`` (with graceful degradation).

    Several databases - VFDB, ARG-ANNOT, EcOH - only carry three ``~~~`` fields
    and append the description to the accession after a space, so the accession
    field has to be split again or the product is lost and the accession becomes
    a sentence.
    """
    fields = header.split("~~~")
    if len(fields) >= 4:
        return {
            "gene": fields[1].strip(),
            "accession": fields[2].strip(),
            "product": "~~~".join(fields[3:]).strip(),
        }
    if len(fields) == 3:
        accession, _, product = fields[2].strip().partition(" ")
        return {"gene": fields[1].strip(), "accession": accession.strip(),
                "product": product.strip()}
    token, _, rest = header.partition(" ")
    return {"gene": token.strip(), "accession": token.strip(), "product": rest.strip()}


def parse_amrprot_header(header: str) -> dict:
    """Parse an AMRFinderPlus ``AMRProt.fa`` header.

    Layout: ``accession|fusion_part|fusion_total|gene_symbol|fam_id|part_of_gene|
    reportable|subclass|class|product_name``. Databases up to the 2025 releases
    prefixed the line with a numeric GI, so that field is dropped when present.
    """
    fields = header.split("|")
    if fields and fields[0].strip().isdigit():
        fields = fields[1:]
    if len(fields) < 10:
        fields = fields + [""] * (10 - len(fields))
    return {
        "accession": fields[0].strip(),
        "fusion_part": fields[1].strip(),
        "fusion_total": fields[2].strip(),
        "gene": fields[3].strip(),
        "fam_id": fields[4].strip(),
        # 'mutation' marks a reference protein that exists only to anchor point
        # mutations; 'hydrolase' and friends describe the gene product.
        "part_of_gene": fields[5].strip(),
        "reportable": fields[6].strip(),
        "subclass": fields[7].strip().replace("_", " ").upper() if fields[7].strip() else "",
        "class": fields[8].strip().replace("_", " ").upper() if fields[8].strip() else "",
        "product": fields[9].replace("_", " ").strip(),
    }


def read_fam_table(path: Path) -> dict[str, dict]:
    """Load ``fam.tsv`` keyed by node id."""
    out: dict[str, dict] = {}
    if not path.exists():
        return out
    with open(path) as handle:
        header = handle.readline().lstrip("#").rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in handle:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < len(header):
                row += [""] * (len(header) - len(row))
            node = row[idx.get("node_id", 0)]
            out[node] = {
                "parent": row[idx["parent_node_id"]] if "parent_node_id" in idx else "",
                "gene": row[idx["gene_symbol"]] if "gene_symbol" in idx else "",
                "type": row[idx["type"]] if "type" in idx else "",
                "subtype": row[idx["subtype"]] if "subtype" in idx else "",
                "class": row[idx["class"]] if "class" in idx else "",
                "subclass": row[idx["subclass"]] if "subclass" in idx else "",
                "family_name": row[idx["family_name"]] if "family_name" in idx else "",
                "reportable": row[idx["reportable"]] if "reportable" in idx else "",
            }
    return out


#: Databases renamed after release, mapped old name -> current name.
LEGACY_DB_NAMES = {"amrfinderplus": "protein"}


class DatabaseStore:
    """Owns ``$HYDRA_DB`` — installation, manifest bookkeeping and lookup."""

    def __init__(self, db_dir: Path | str):
        self.root = Path(db_dir)
        self._manifest: dict | None = None
        self._cache: dict[str, DbHandle] = {}

    # ------------------------------------------------------------------ manifest
    @property
    def manifest_path(self) -> Path:
        return self.root / MANIFEST_NAME

    def manifest(self, refresh: bool = False) -> dict:
        if self._manifest is None or refresh:
            data: dict | None = None
            if self.manifest_path.exists():
                # An import running in another process rewrites this file; retry
                # briefly rather than failing a concurrent analysis run.
                last: Exception | None = None
                for attempt in range(5):
                    try:
                        with open(self.manifest_path) as handle:
                            data = json.load(handle)
                        break
                    except (OSError, json.JSONDecodeError) as exc:
                        last = exc
                        time.sleep(0.05 * (attempt + 1))
                if data is None:
                    raise HydraError(f"database manifest is unreadable ({last}); "
                                     f"re-run 'hydra db import' to rebuild it")
            self._manifest = data if isinstance(data, dict) else {}
            self._manifest.setdefault("hydra_db_version", 1)
            # A manifest written by a future version may lack the key entirely.
            if not isinstance(self._manifest.get("databases"), dict):
                self._manifest["databases"] = {}
            # Stores built before the protein reference was renamed. The recorded
            # path is left alone, so the files stay where they were installed.
            for old, new in LEGACY_DB_NAMES.items():
                databases = self._manifest["databases"]
                if old in databases and new not in databases:
                    databases[new] = databases.pop(old)
        return self._manifest

    def _save_manifest(self) -> None:
        """Write the manifest atomically so a reader never sees a partial file."""
        self.root.mkdir(parents=True, exist_ok=True)
        temporary = self.manifest_path.with_name(f".{self.manifest_path.name}.{os.getpid()}")
        with open(temporary, "w") as handle:
            json.dump(self.manifest(), handle, indent=2, sort_keys=True)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, self.manifest_path)

    def record(self, name: str, **fields) -> None:
        # Re-read first: another process may have registered a database since
        # this one cached the manifest, and a blind write would drop it.
        if self.manifest_path.exists():
            try:
                self.manifest(refresh=True)
            except HydraError:
                pass
        entry = self.manifest()["databases"].setdefault(name, {})
        entry.update(fields)
        entry["installed"] = time.strftime("%Y-%m-%d %H:%M:%S")
        self._save_manifest()

    def installed(self) -> dict[str, dict]:
        return dict(self.manifest()["databases"])

    def is_installed(self, name: str) -> bool:
        entry = self.manifest()["databases"].get(name)
        if not entry:
            return False
        path = self.root / entry.get("path", "")
        return path.exists()

    def require_installed(self, name: str) -> dict:
        entry = self.manifest()["databases"].get(name)
        if not entry:
            raise HydraError(
                f"database '{name}' is not installed in {self.root}.\n"
                f"Install it with:  hydra db import         (from local conda environments)\n"
                f"                  hydra db download {name}  (from the upstream source)"
            )
        return entry

    # ------------------------------------------------------------------ discovery
    @staticmethod
    def conda_envs() -> list[Path]:
        """All conda environment roots visible on this machine."""
        roots: list[Path] = []
        prefix = os.environ.get("CONDA_PREFIX")
        if prefix:
            roots.append(Path(prefix))
            envs = Path(prefix).parent
            if envs.name == "envs":
                roots.append(envs.parent)
        for base in ("~/miniconda3", "~/anaconda3", "~/miniforge3", "~/mambaforge",
                     "/opt/conda", "/usr/local/conda"):
            path = Path(base).expanduser()
            if path.exists():
                roots.append(path)
        out: list[Path] = []
        seen: set[Path] = set()
        for root in roots:
            candidates = [root] + sorted((root / "envs").glob("*")) if (root / "envs").exists() else [root]
            for cand in candidates:
                if cand.is_dir() and cand not in seen:
                    seen.add(cand)
                    out.append(cand)
        return out

    def find_source(self, spec: DbSpec, search_paths: list[Path] | None = None) -> Path | None:
        """Locate the on-disk source for *spec* in conda environments or explicit paths.

        A directory the user named with ``--source`` always wins. Otherwise the
        *newest* database version wins: several conda environments commonly hold
        different AMRFinderPlus releases, and importing a year-old catalogue
        because its path happened to be shorter would be silently wrong.
        """
        explicit: list[Path] = []
        for base in (search_paths or []):
            explicit.append(Path(base))
            if spec.conda_rel_path:
                explicit.append(Path(base) / spec.conda_rel_path)
        for cand in explicit:
            if cand.exists() and self._source_valid(spec, cand):
                return cand

        candidates: list[Path] = []
        for env in self.conda_envs():
            rel = spec.conda_rel_path
            if not rel:
                continue
            if "*" in rel:
                candidates.extend(Path(p) for p in glob.glob(str(env / rel)))
            else:
                candidates.append(env / rel)

        def rank(path: Path) -> tuple:
            version = self._detect_version(path)
            hinted = any(hint in str(path) for hint in spec.conda_env_hint) \
                if spec.conda_env_hint else False
            # Sorted descending, so a newer version and then a hinted environment
            # come first; path length only breaks a genuine tie.
            return (version if version != "unknown" else "", 1 if hinted else 0,
                    -len(str(path)))

        usable = [c for c in sorted({c for c in candidates})
                  if c.exists() and self._source_valid(spec, c)]
        if not usable:
            return None
        usable.sort(key=rank, reverse=True)
        if len(usable) > 1:
            LOG.debug("%s: chose %s (version %s) from %d candidates", spec.name, usable[0],
                      self._detect_version(usable[0]), len(usable))
        return usable[0]

    @staticmethod
    def _source_valid(spec: DbSpec, path: Path) -> bool:
        if spec.kind == "nucl":
            return (path / "sequences").exists() or (path / "sequences.fa").exists() or \
                   any(path.glob("*.fa")) or any(path.glob("*.fna"))
        if spec.kind == "prot":
            return (path / "AMRProt.fa").exists() or (path / "AMRProt").exists()
        if spec.kind == "mlst":
            return (path / "pubmlst").is_dir() and (path / "scheme_species_map.tab").exists()
        if spec.kind == "typing":
            return path.is_dir() and any(path.glob("*mlst*")) or any(path.glob("klebsiella*"))
        if spec.kind == "sketch":
            return any(path.glob("*.msh"))
        return path.exists()

    # ------------------------------------------------------------------ install
    def import_all(self, only: list[str] | None = None, search_paths: list[Path] | None = None,
                   force: bool = False) -> dict[str, str]:
        """Import every database that can be found locally. Returns name -> status."""
        results: dict[str, str] = {}
        names = only or list(DATABASES)
        for name in names:
            spec = spec_for(name)
            if self.is_installed(name) and not force:
                results[name] = "already installed"
                continue
            source = self.find_source(spec, search_paths)
            if source is None:
                results[name] = "source not found"
                continue
            try:
                self.import_one(spec, source)
                results[name] = f"installed from {source}"
            except HydraError as exc:
                results[name] = f"FAILED: {exc}"
                LOG.warning("import of %s failed: %s", name, exc)
        return results

    def download(self, name: str, force: bool = False, cache: Path | None = None) -> None:
        """Fetch *name* from its upstream source and import it.

        Staging happens in a temporary directory that is removed afterwards, so
        a failed download never leaves a half-built database behind: the import
        only runs once every file has arrived.
        """
        from .fetch import stage  # imported lazily: only a download needs urllib

        spec = spec_for(name)
        if self.is_installed(name) and not force:
            LOG.info("%s is already installed", name)
            return
        work = Path(cache) if cache else self.root / ".download" / name
        if work.exists():
            shutil.rmtree(work)
        work.mkdir(parents=True, exist_ok=True)
        try:
            source = stage(name, work, progress=lambda msg: LOG.debug("%s", msg))
            self.import_one(spec, source)
        finally:
            shutil.rmtree(work, ignore_errors=True)

    def import_one(self, spec: DbSpec, source: Path) -> None:
        LOG.info("importing %s from %s", spec.name, source)
        if spec.kind == "nucl":
            self._import_nucl(spec, source)
        elif spec.kind == "prot":
            self._import_prot(spec, source)
        elif spec.kind == "mlst":
            self._import_mlst(spec, source)
        elif spec.kind == "typing":
            self._import_typing(spec, source)
        elif spec.kind == "sketch":
            self._import_sketch(spec, source)
        else:
            raise HydraError(f"do not know how to import a '{spec.kind}' database")

    def _class_lookup(self) -> dict[str, tuple[str, str, str, str, str]]:
        """gene-key -> (class, subclass, element_type, element_subtype, fam_id).

        Built from the AMRFinderPlus protein reference so that every other
        database inherits its curated drug-class annotation and, importantly for
        read mode, its gene-family grouping: ``blaTEM-1`` and ``blaTEM-217``
        both belong to family ``blaTEM`` and describe one locus, not two.
        """
        cached = getattr(self, "_class_lookup_cache", None)
        if cached is not None:
            return cached
        lookup: dict[str, tuple[str, str, str, str, str]] = {}
        prot_dir = protein_dir(self.root)
        fam = read_fam_table(prot_dir / "fam.tsv")
        # Read the normalised metadata table, not AMRProt.fa: the installed FASTA
        # carries short synthetic ids, so its headers no longer hold annotation.
        for row in self.load_meta(prot_dir / "meta.tsv").values():
            famrow = fam.get(row.get("fam_id", ""), {})
            etype = row.get("element_type") or famrow.get("type") or (
                AMR if row.get("class") else "")
            esub = row.get("element_subtype") or famrow.get("subtype") or ""
            entry = (row.get("class", ""), row.get("subclass", ""), etype, esub,
                     row.get("fam_id", ""))
            for key in {_normalise_gene(row.get("gene", "")),
                        _normalise_gene(row.get("fam_id", ""))}:
                if key and key not in lookup:
                    lookup[key] = entry
        for node, row in fam.items():
            key = _normalise_gene(row.get("gene") or node)
            if key and key not in lookup and (row.get("class") or row.get("type")):
                lookup[key] = (row.get("class", ""), row.get("subclass", ""),
                               row.get("type", ""), row.get("subtype", ""), node)
        self._class_lookup_cache = lookup
        return lookup

    def _import_nucl(self, spec: DbSpec, source: Path) -> None:
        src_fasta = None
        for candidate in ("sequences", "sequences.fa", "sequences.fna"):
            if (source / candidate).exists():
                src_fasta = source / candidate
                break
        if src_fasta is None:
            hits = sorted(list(source.glob("*.fa")) + list(source.glob("*.fna")))
            if hits:
                src_fasta = hits[0]
        if src_fasta is None:
            raise HydraError(f"no sequence file found in {source}")

        dest = self.root / "nucl" / spec.name
        dest.mkdir(parents=True, exist_ok=True)
        out_fasta = dest / "sequences.fna"
        class_lookup = self._class_lookup()

        records: list[tuple[str, str]] = []
        rows: list[dict] = []
        index = 0
        seen_seqids: set[str] = set()
        for header, seq in read_fasta(src_fasta):
            if not seq:
                continue
            info = _parse_abricate_header(header)
            seqid = f"{spec.name}_{index:06d}"
            index += 1
            seen_seqids.add(seqid)
            annotation = ("", "", "", "", "")
            for gene_key in _lookup_keys(info["gene"]):
                if gene_key in class_lookup:
                    annotation = class_lookup[gene_key]
                    break
            klass, subclass, etype, esub, fam_id = annotation
            if not klass:
                klass = _infer_class(f"{info['gene']} {info['product']}")
            if not etype:
                etype = spec.element_type
            default_subtype = {AMR: "AMR", VIRULENCE: "VIRULENCE", STRESS: "STRESS",
                               PLASMID: "REPLICON", SEROTYPE: "SEROTYPE"}.get(
                spec.element_type, spec.element_type)
            if spec.element_type in (VIRULENCE, PLASMID, SEROTYPE):
                # Never let an AMR lookup relabel an explicitly non-AMR database:
                # a virulence gene must not come out as (VIRULENCE, POINT), which
                # the mutation caller reads as a mutation reference.
                etype = spec.element_type
                esub = default_subtype
                if klass and klass not in ("", default_subtype):
                    subclass = ""
                    klass = ""
            if not esub:
                esub = default_subtype
            # 'POINT' marks a mutation-reference protein; an acquired gene in a
            # nucleotide database is never one.
            if esub == "POINT":
                esub = default_subtype
            records.append((seqid, seq))
            rows.append({
                "seqid": seqid, "gene": info["gene"], "accession": info["accession"],
                "product": info["product"], "class": klass, "subclass": subclass,
                "element_type": etype, "element_subtype": esub, "fam_id": fam_id,
            })
        if not records:
            raise HydraError(f"{src_fasta} contained no usable sequences")

        write_fasta(out_fasta, records, wrap=0)
        self._write_meta(dest / "meta.tsv", rows)
        makeblastdb(out_fasta, "nucl")
        version = self._detect_version(source)
        self.record(spec.name, kind="nucl", path=str(Path("nucl") / spec.name),
                    source=str(source), version=version, sequences=len(records),
                    element_type=spec.element_type, title=spec.title)

    def _import_prot(self, spec: DbSpec, source: Path) -> None:
        dest = self.root / "prot" / spec.name
        dest.mkdir(parents=True, exist_ok=True)
        src_prot = source / "AMRProt.fa"
        if not src_prot.exists():
            src_prot = source / "AMRProt"
        if not src_prot.exists():
            raise HydraError(f"AMRProt.fa not found under {source}")

        records: list[tuple[str, str]] = []
        rows: list[dict] = []
        fam = read_fam_table(source / "fam.tsv")
        for i, (header, seq) in enumerate(read_fasta(src_prot)):
            if not seq:
                continue
            info = parse_amrprot_header(header)
            seqid = f"afp_{i:06d}"
            famrow = fam.get(info["fam_id"], {})
            is_mutation_reference = info["part_of_gene"].lower() == "mutation"
            records.append((seqid, seq))
            rows.append({
                "seqid": seqid, "gene": info["gene"], "accession": info["accession"],
                "product": info["product"], "class": info["class"], "subclass": info["subclass"],
                "element_type": famrow.get("type", AMR) or AMR,
                "element_subtype": ("POINT" if is_mutation_reference
                                    else (famrow.get("subtype", "AMR") or "AMR")),
                "fam_id": info["fam_id"], "reportable": info["reportable"] or famrow.get("reportable", ""),
                "fusion_part": info["fusion_part"], "fusion_total": info["fusion_total"],
                "length": str(len(seq)), "part_of_gene": info["part_of_gene"],
            })
        write_fasta(dest / "AMRProt.fa", records, wrap=0)
        self._write_meta(dest / "meta.tsv", rows,
                         columns=META_COLUMNS + ("fam_id", "reportable", "fusion_part",
                                                 "fusion_total", "length", "part_of_gene"))
        makeblastdb(dest / "AMRProt.fa", "prot")

        for aux in ("fam.tsv", "AMRProt-mutation.tsv", "AMRProt-suppress.tsv",
                    "AMRProt-susceptible.tsv", "taxgroup.tsv", "version.txt",
                    "database_format_version.txt"):
            if (source / aux).exists():
                copy_data(source / aux, dest / aux)

        mut_dir = self.root / "mutation" / "dna"
        mut_dir.mkdir(parents=True, exist_ok=True)
        organisms = []
        for fa in sorted(source.glob("AMR_DNA-*.fa")):
            organism = fa.name[len("AMR_DNA-"):-len(".fa")]
            tsv = source / f"AMR_DNA-{organism}.tsv"
            copy_data(fa, mut_dir / f"{organism}.fna")
            if tsv.exists():
                copy_data(tsv, mut_dir / f"{organism}.tsv")
            makeblastdb(mut_dir / f"{organism}.fna", "nucl")
            organisms.append(organism)

        version = self._detect_version(source)
        self.record(spec.name, kind="prot", path=str(Path("prot") / spec.name),
                    source=str(source), version=version, sequences=len(records),
                    organisms=organisms, element_type=spec.element_type, title=spec.title)
        self._class_lookup_cache = None  # rebuild now that AMRProt is available

    def _import_mlst(self, spec: DbSpec, source: Path) -> None:
        dest = self.root / "mlst"
        dest.mkdir(parents=True, exist_ok=True)
        pubmlst_src = source / "pubmlst"
        if not pubmlst_src.is_dir():
            raise HydraError(f"pubmlst directory not found under {source}")
        pubmlst_dest = dest / "pubmlst"
        if pubmlst_dest.exists():
            shutil.rmtree(pubmlst_dest)
        shutil.copytree(pubmlst_src, pubmlst_dest,
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        for aux in ("scheme_species_map.tab", "VERSION.txt"):
            if (source / aux).exists():
                copy_data(source / aux, dest / aux)

        # Build one combined BLAST database over all schemes: >scheme.locus_allele
        combined = dest / "blast" / "mlst.fna"
        combined.parent.mkdir(parents=True, exist_ok=True)
        records: list[tuple[str, str]] = []
        schemes: list[str] = []
        for scheme_dir in sorted(p for p in pubmlst_dest.iterdir() if p.is_dir()):
            scheme = scheme_dir.name
            n_before = len(records)
            for tfa in sorted(scheme_dir.glob("*.tfa")):
                locus = tfa.stem
                for header, seq in read_fasta(tfa):
                    if not seq:
                        continue
                    allele = header.split()[0]
                    allele_id = allele[len(locus):].lstrip("_-") if allele.startswith(locus) else allele
                    records.append((f"{scheme}.{locus}.{allele_id or allele}", seq))
            if len(records) > n_before:
                schemes.append(scheme)
        if not records:
            raise HydraError(f"no MLST alleles found under {pubmlst_src}")
        write_fasta(combined, records, wrap=0)
        makeblastdb(combined, "nucl")
        version = self._detect_version(source)
        self.record(spec.name, kind="mlst", path="mlst", source=str(source), version=version,
                    sequences=len(records), schemes=len(schemes), title=spec.title)

    def _import_typing(self, spec: DbSpec, source: Path) -> None:
        dest = self.root / "typing" / spec.name
        dest.mkdir(parents=True, exist_ok=True)
        copied: list[str] = []
        combined: list[tuple[str, str]] = []
        for module_dir in sorted(p for p in Path(source).iterdir() if p.is_dir()):
            data_dir = module_dir / "data"
            if not data_dir.is_dir():
                continue
            fastas = sorted(list(data_dir.glob("*.fasta")) + list(data_dir.glob("*.fa")))
            profiles = data_dir / "profiles.tsv"
            if not fastas:
                continue
            module = module_dir.name
            target = dest / module
            target.mkdir(parents=True, exist_ok=True)
            records: list[tuple[str, str]] = []
            for fa in fastas:
                locus = fa.stem
                for header, seq in read_fasta(fa):
                    if not seq:
                        continue
                    allele = header.split()[0]
                    allele_id = allele[len(locus):].lstrip("_-") if allele.startswith(locus) else allele
                    records.append((f"{locus}.{allele_id or allele}", seq))
            if not records:
                continue
            duplicates = len(records) - len({seqid for seqid, _seq in records})
            if duplicates:
                LOG.warning("%s: %d duplicate allele id(s); the upstream scheme has two "
                            "records with the same name, so their hits are ambiguous",
                            module, duplicates)
            write_fasta(target / "alleles.fna", records, wrap=0)
            if profiles.exists():
                copy_data(profiles, target / "profiles.tsv")
            # '#' separates the module from locus.allele so a single BLAST pass
            # can serve every scheme at once.
            combined.extend((f"{module}#{seqid}", seq) for seqid, seq in records)
            copied.append(module)
        if not copied:
            raise HydraError(f"no typing schemes with allele data found under {source}")
        blast_dir = dest / "blast"
        blast_dir.mkdir(parents=True, exist_ok=True)
        write_fasta(blast_dir / "alleles.fna", combined, wrap=0)
        makeblastdb(blast_dir / "alleles.fna", "nucl")
        self.record(spec.name, kind="typing", path=str(Path("typing") / spec.name),
                    source=str(source), modules=copied, sequences=len(combined),
                    title=spec.title)

    def _import_sketch(self, spec: DbSpec, source: Path) -> None:
        dest = self.root / "sketch" / spec.name
        dest.mkdir(parents=True, exist_ok=True)
        found = []
        for msh in sorted(Path(source).glob("*.msh")):
            copy_data(msh, dest / msh.name)
            found.append(msh.name)
        if not found:
            raise HydraError(f"no .msh sketches found under {source}")
        self.record(spec.name, kind="sketch", path=str(Path("sketch") / spec.name),
                    source=str(source), sketches=found, title=spec.title)

    @staticmethod
    def _detect_version(source: Path) -> str:
        for candidate in ("version.txt", "VERSION.txt", "database_version.txt"):
            path = source / candidate
            if path.exists():
                try:
                    return path.read_text().strip().splitlines()[0]
                except (OSError, IndexError):
                    pass
        # AMRFinderPlus data dirs are named by release date
        if re.fullmatch(r"\d{4}-\d{2}-\d{2}(\.\d+)?", source.name):
            return source.name
        resolved = source.resolve()
        if re.fullmatch(r"\d{4}-\d{2}-\d{2}(\.\d+)?", resolved.name):
            return resolved.name
        return "unknown"

    @staticmethod
    def _write_meta(path: Path, rows: list[dict], columns: tuple[str, ...] = META_COLUMNS) -> None:
        with open(path, "w") as handle:
            handle.write("\t".join(columns) + "\n")
            for row in rows:
                handle.write("\t".join(str(row.get(col, "")).replace("\t", " ") for col in columns) + "\n")

    # ------------------------------------------------------------------ lookup
    def load_meta(self, path: Path) -> dict[str, dict]:
        meta: dict[str, dict] = {}
        if not path.exists():
            return meta
        with open(path) as handle:
            columns = handle.readline().rstrip("\n").split("\t")
            for line in handle:
                if not line.strip():
                    continue
                values = line.rstrip("\n").split("\t")
                if len(values) < len(columns):
                    values += [""] * (len(columns) - len(values))
                row = dict(zip(columns, values))
                meta[row["seqid"]] = row
        return meta

    def handle(self, name: str) -> DbHandle:
        """Return an installed database, loading its metadata once and caching it."""
        if name in self._cache:
            return self._cache[name]
        entry = self.require_installed(name)
        path = self.root / entry["path"]
        kind = entry.get("kind", "nucl")
        if kind == "nucl":
            fasta = path / "sequences.fna"
        elif kind == "prot":
            fasta = path / "AMRProt.fa"
        else:
            fasta = None
        if fasta is not None and not fasta.exists():
            raise HydraError(f"database '{name}' is registered but its sequences are missing "
                             f"({fasta}). Re-run: hydra db import --force {name}")
        handle = DbHandle(
            name=name, kind=kind, path=path,
            element_type=entry.get("element_type", AMR), fasta=fasta,
            meta=self.load_meta(path / "meta.tsv") if kind in ("nucl", "prot") else None,
            version=entry.get("version", ""), n_sequences=int(entry.get("sequences", 0) or 0),
        )
        self._cache[name] = handle
        return handle

    def mutation_organisms(self) -> list[str]:
        mut_dir = self.root / "mutation" / "dna"
        if not mut_dir.is_dir():
            return []
        return sorted(p.stem for p in mut_dir.glob("*.fna"))

    def summary_rows(self) -> list[dict]:
        rows = []
        for name, entry in sorted(self.installed().items()):
            path = self.root / entry.get("path", "")
            size = 0
            if path.exists():
                for root, _dirs, files in os.walk(path):
                    for fname in files:
                        try:
                            size += os.path.getsize(os.path.join(root, fname))
                        except OSError:
                            pass
            spec = DATABASES.get(name)
            rows.append({
                "name": name,
                "kind": entry.get("kind", ""),
                "sequences": entry.get("sequences", entry.get("schemes", "")),
                "version": entry.get("version", ""),
                "size": human_bytes(size),
                "installed": entry.get("installed", ""),
                "title": entry.get("title", spec.title if spec else ""),
            })
        return rows


BUNDLE_MANIFEST = "hydra-bundle.json"


def _directory_digest(path: Path) -> str:
    """Content hash of a directory tree: sorted relative paths plus file bytes."""
    import hashlib

    digest = hashlib.sha256()
    files: list[tuple[str, str]] = []
    for root, dirs, names in os.walk(path):
        dirs.sort()
        for fname in sorted(names):
            full = os.path.join(root, fname)
            files.append((os.path.relpath(full, path).replace(os.sep, "/"), full))
    for relative, full in sorted(files):
        digest.update(relative.encode())
        digest.update(b"\0")
        with open(full, "rb") as handle:
            for block in iter(lambda h=handle: h.read(1 << 20), b""):
                digest.update(block)
    return digest.hexdigest()


def create_bundle(store: "DatabaseStore", output: Path, names: list[str] | None = None,
                  compress: str = "gz") -> Path:
    """Pack installed databases into a single portable archive.

    A lab builds the archive once on a machine that has the source tools, then
    every other machine installs the identical databases with
    ``hydra db download --from-file``. Each database's content hash is recorded
    and checked on install, so a truncated or tampered download is caught.
    """
    import io
    import tarfile

    output = Path(output)
    installed = store.installed()
    selected = [n for n in (names or list(installed)) if n in installed]
    if not selected:
        raise HydraError("no installed databases to bundle; run 'hydra db import' first")
    mode = {"gz": "w:gz", "bz2": "w:bz2", "xz": "w:xz", "none": "w"}.get(compress)
    if mode is None:
        raise HydraError(f"unknown compression '{compress}'; use gz, bz2, xz or none")

    bundle = {"hydra_bundle_version": 1, "created": time.strftime("%Y-%m-%d %H:%M:%S"),
              "databases": {}, "checksums": {}}
    LOG.info("bundling %d databases into %s", len(selected), output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(output, mode) as tar:
        for name in selected:
            rel = installed[name].get("path", "")
            path = store.root / rel
            if not rel or not path.exists():
                LOG.warning("skipping %s: %s is missing", name, path)
                continue
            # Record only what is actually in the archive, so the manifest can
            # never advertise a database whose files were never added.
            bundle["databases"][name] = installed[name]
            bundle["checksums"][name] = _directory_digest(path)
            tar.add(path, arcname=rel)
        if not bundle["databases"]:
            raise HydraError("none of the selected databases have files on disk")
        manifest_bytes = json.dumps(bundle, indent=2, sort_keys=True).encode()
        info = tarfile.TarInfo(BUNDLE_MANIFEST)
        info.size = len(manifest_bytes)
        info.mtime = int(time.time())
        tar.addfile(info, io.BytesIO(manifest_bytes))
    LOG.info("wrote %s (%s)", output, human_bytes(output.stat().st_size))
    return output


def _safe_members(tar, root: Path, wanted_prefixes: list[str]):
    """Archive members that belong to *wanted_prefixes* and cannot escape *root*.

    A bundle can come from anywhere - a URL, a shared drive - so every member is
    checked: no absolute paths, no ``..`` traversal, and no links, which would
    otherwise let an archive write through a symlink to any file the user owns.
    """
    resolved_root = root.resolve()
    for member in tar.getmembers():
        if member.name == BUNDLE_MANIFEST:
            continue
        name = member.name.replace("\\", "/")
        if member.issym() or member.islnk():
            raise HydraError(f"refusing to extract link member {member.name!r} from the bundle")
        if not (member.isfile() or member.isdir()):
            raise HydraError(f"refusing to extract special member {member.name!r}")
        if name.startswith("/") or os.path.isabs(name) or ".." in Path(name).parts:
            raise HydraError(f"refusing to extract {member.name!r}: path escapes {root}")
        target = (root / name).resolve()
        if target != resolved_root and resolved_root not in target.parents:
            raise HydraError(f"refusing to extract {member.name!r}: path escapes {root}")
        if any(name == prefix or name.startswith(prefix.rstrip("/") + "/")
               for prefix in wanted_prefixes):
            yield member


def install_bundle(store: "DatabaseStore", archive: Path, force: bool = False) -> list[str]:
    """Unpack a bundle created by :func:`create_bundle` into the database root."""
    import tarfile

    archive = Path(archive)
    if not archive.exists():
        raise HydraError(f"bundle not found: {archive}")
    store.root.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:*") as tar:
        try:
            manifest_member = tar.getmember(BUNDLE_MANIFEST)
        except KeyError as exc:
            raise HydraError(f"{archive} is not a Hydra database bundle "
                             f"(no {BUNDLE_MANIFEST})") from exc
        with tar.extractfile(manifest_member) as handle:
            bundle = json.load(handle)
        databases = bundle.get("databases", {})
        selected = []
        for name in databases:
            if store.is_installed(name) and not force:
                LOG.info("%s is already installed; use --force to replace it", name)
                continue
            selected.append(name)
        if not selected:
            return []
        # Extract only the databases being installed: the others must keep the
        # files that match the manifest entry already recorded for them.
        prefixes = [databases[name].get("path", "") for name in selected]
        prefixes = [p for p in prefixes if p]
        members = list(_safe_members(tar, store.root, prefixes))
        LOG.info("extracting %d database(s) from %s", len(selected), archive.name)
        tar.extractall(store.root, members=members)

    checksums = bundle.get("checksums", {})
    for name in selected:
        expected = checksums.get(name)
        path = store.root / databases[name].get("path", "")
        if expected:
            actual = _directory_digest(path)
            if actual != expected:
                raise HydraError(
                    f"checksum mismatch for '{name}' after extracting {archive.name}: "
                    f"the bundle is truncated or has been modified. Nothing was recorded "
                    f"for it; delete {path} and fetch the bundle again.")
        else:
            LOG.warning("bundle carries no checksum for '%s'; installing unverified", name)
        store.record(name, **databases[name])
    return selected


def download_bundle(url: str, destination: Path) -> Path:
    """Fetch a bundle over HTTP(S) to *destination*."""
    import urllib.error
    import urllib.request

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    LOG.info("downloading %s", url)
    try:
        with urllib.request.urlopen(url, timeout=60) as response:  # noqa: S310 - user-supplied URL
            total = int(response.headers.get("Content-Length") or 0)
            written = 0
            with open(destination, "wb") as out:
                while True:
                    block = response.read(1 << 20)
                    if not block:
                        break
                    out.write(block)
                    written += len(block)
                    if total:
                        LOG.debug("  %s / %s", human_bytes(written), human_bytes(total))
    except urllib.error.URLError as exc:
        raise HydraError(
            f"could not download {url}: {exc}\n"
            f"If this machine has no network access, build a bundle on one that does:\n"
            f"  hydra db bundle -o hydra-db.tar.gz\n"
            f"then install it here with:\n"
            f"  hydra db download --from-file hydra-db.tar.gz"
        ) from exc
    LOG.info("downloaded %s (%s)", destination, human_bytes(destination.stat().st_size))
    return destination


def copy_data(source: Path, destination: Path) -> None:
    """Copy a reference file and stamp it with the current time.

    ``copy2`` would preserve the source mtime, and upstream releases delivered by
    conda, tar or rsync keep theirs too. An index built from the *previous*
    release would then look newer than the freshly copied FASTA and be reused -
    silently applying the new mutation coordinates to the old sequences.
    """
    shutil.copy2(source, destination)
    os.utime(destination, None)


def makeblastdb(fasta: Path, dbtype: str, force: bool = False) -> None:
    """Build a BLAST index next to *fasta* (skipped when demonstrably current)."""
    fasta = Path(fasta)
    suffix = ".nin" if dbtype == "nucl" else ".pin"
    candidates = [Path(f"{fasta}{suffix}"), Path(f"{fasta}.00{suffix}")]
    if not force:
        for candidate in candidates:
            if candidate.exists() and candidate.stat().st_mtime >= fasta.stat().st_mtime:
                return
    require("makeblastdb", "building BLAST indexes")
    run(["makeblastdb", "-in", str(fasta), "-dbtype", dbtype, "-out", str(fasta), "-hash_index"])
