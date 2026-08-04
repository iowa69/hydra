"""Lineage / sublineage typing and genome-level risk scores.

The allele-based schemes (yersiniabactin, colibactin, aerobactin, salmochelin,
*rmpA*) come from Kleborate and are typed exactly as MLST is: best allele per
locus, then a profile lookup that also yields the lineage label.

The scores are Hydra's generalisation of Kleborate's. The virulence score keeps
Kleborate's published 0-5 definition for the *Klebsiella pneumoniae* species
complex, and applies analogous species-specific rules elsewhere. The resistance
score is computed from AMRFinderPlus drug-class annotation rather than a
Klebsiella-specific gene list, so it is meaningful for any Gram-negative.
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from ..config import Config
from ..db.manager import DatabaseStore
from ..records import Hit, SpeciesCall, TypingResult
from ..utils import LOG
from ..engines.blast import blast, check_db_exists, merge_hsps
from ..engines.nucl import QueryBatch

MIN_ALLELE_IDENTITY = 90.0
MIN_ALLELE_COVERAGE = 80.0

#: Human-readable names and the genera each scheme applies to. A module absent
#: from this table still runs, but only when the species is unknown.
MODULE_INFO: dict[str, dict] = {
    "klebsiella__ybst": {"label": "ybST", "locus": "Yersiniabactin",
                         "genera": ("Klebsiella", "Raoultella", "Enterobacter", "Citrobacter",
                                    "Escherichia", "Serratia")},
    "klebsiella__cbst": {"label": "cbST", "locus": "Colibactin",
                         "genera": ("Klebsiella", "Escherichia", "Citrobacter", "Enterobacter")},
    "klebsiella__abst": {"label": "AbST", "locus": "Aerobactin",
                         "genera": ("Klebsiella", "Escherichia", "Enterobacter", "Salmonella")},
    "klebsiella__smst": {"label": "SmST", "locus": "Salmochelin",
                         "genera": ("Klebsiella", "Escherichia", "Enterobacter")},
    "klebsiella__rmst": {"label": "RmST", "locus": "RmpADC", "genera": ("Klebsiella",)},
    "klebsiella__rmpa2": {"label": "rmpA2", "locus": "rmpA2", "genera": ("Klebsiella",)},
    "klebsiella_pneumo_complex__wzi": {"label": "wzi", "locus": "wzi",
                                       "genera": ("Klebsiella",)},
    "escherichia__ectyper": {"label": "serotype", "locus": "O/H antigen",
                             "genera": ("Escherichia", "Shigella")},
    "escherichia__pathovar": {"label": "pathovar_markers", "locus": "pathovar",
                              "genera": ("Escherichia", "Shigella")},
    "escherichia__mlst_achtman": {"label": "ST_achtman", "locus": "MLST",
                                  "genera": ("Escherichia", "Shigella")},
    "escherichia__mlst_pasteur": {"label": "ST_pasteur", "locus": "MLST",
                                  "genera": ("Escherichia", "Shigella")},
    "escherichia__mlst_lee": {"label": "ST_lee", "locus": "MLST",
                              "genera": ("Escherichia", "Shigella")},
}

#: Modules that ship allele FASTAs but are not lineage schemes. The ``__amr``
#: modules are resistance gene catalogues, which Hydra screens itself, and the
#: complex-specific ``__mlst`` modules duplicate the PubMLST schemes.
EXCLUDED_MODULES: frozenset[str] = frozenset({
    "klebsiella_pneumo_complex__amr",
    "escherichia__amr",
    "klebsiella_pneumo_complex__mlst",
    "klebsiella_oxytoca_complex__mlst",
})

#: Loci that describe virulence rather than lineage. Everything else in
#: MODULE_INFO is a typing scheme and must not raise a virulence score.
VIRULENCE_LOCI = frozenset({
    "Yersiniabactin", "Colibactin", "Aerobactin", "Salmochelin", "RmpADC", "rmpA2",
})

#: Virulence-score rules keyed by the genus (or genus+species) they apply to.
#: Each rule maps a score onto the set of loci that must be present.
VIRULENCE_RULES: dict[str, list[tuple[int, tuple[str, ...], tuple[str, ...]]]] = {
    # (score, all-of, none-of) evaluated highest score first
    "Klebsiella": [
        (5, ("Colibactin", "Aerobactin"), ()),
        (4, ("Aerobactin", "Yersiniabactin"), ()),
        (3, ("Aerobactin",), ()),
        (2, ("Colibactin",), ()),
        (1, ("Yersiniabactin",), ()),
    ],
    "Escherichia": [
        (3, ("Colibactin", "Aerobactin"), ()),
        (2, ("Aerobactin",), ()),
        (1, ("Yersiniabactin",), ()),
    ],
}


@dataclass
class SchemeProfile:
    loci: list[str]
    table: dict[tuple[str, ...], str]
    lineage: dict[str, str]
    lineage_column: str = ""


def _load_profiles(path: Path) -> SchemeProfile:
    loci: list[str] = []
    table: dict[tuple[str, ...], str] = {}
    lineage: dict[str, str] = {}
    lineage_column = ""
    if not path.exists():
        return SchemeProfile(loci, table, lineage)
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for name in header[1:]:
            if "lineage" in name.lower() or name.lower() in ("clonal_complex", "cc"):
                lineage_column = name
                continue
            loci.append(name)
        for line in handle:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < len(loci) + 1:
                continue
            st = row[0]
            table[tuple(row[1:len(loci) + 1])] = st
            if lineage_column and len(row) > len(loci) + 1:
                lineage[st] = row[len(loci) + 1].strip()
    return SchemeProfile(loci, table, lineage, lineage_column)


class LineageTyper:
    """Types every installed lineage scheme in a single BLAST pass."""

    def __init__(self, store: DatabaseStore, config: Config):
        self.store = store
        self.config = config
        entry = store.require_installed("lineage")
        self.root = store.root / entry["path"]
        self.modules: list[str] = [m for m in entry.get("modules", [])
                                   if m not in EXCLUDED_MODULES]
        self.fasta = self.root / "blast" / "alleles.fna"
        self._profiles: dict[str, SchemeProfile] = {}

    def profiles(self, module: str) -> SchemeProfile:
        if module not in self._profiles:
            self._profiles[module] = _load_profiles(self.root / module / "profiles.tsv")
        return self._profiles[module]

    def applicable_modules(self, species: SpeciesCall | None) -> list[str]:
        """Modules worth reporting for a species (all of them when species is unknown)."""
        genus = species.genus if species else ""
        if not genus:
            return list(self.modules)
        out = []
        for module in self.modules:
            info = MODULE_INFO.get(module)
            if info is None or genus in info["genera"]:
                out.append(module)
        return out

    def type_batch(self, batch: QueryBatch, workdir: Path, threads: int | None = None,
                   species_by_sample: dict[str, SpeciesCall] | None = None,
                   ) -> dict[str, list[TypingResult]]:
        check_db_exists(self.fasta, "nucl")
        threads = threads or self.config.threads
        out_tab = workdir / "lineage.blastn.tsv"
        LOG.info("lineage typing: %d schemes", len(self.modules))
        hsps = blast(
            "blastn", batch.path, self.fasta, out_tab,
            threads=threads, evalue=1e-20, task="megablast", word_size=20,
            perc_identity=MIN_ALLELE_IDENTITY - 5.0, max_target_seqs=10000,
        )
        merged = merge_hsps(hsps)

        # sample -> module -> locus -> (allele, identity, coverage, bitscore)
        found: dict[str, dict[str, dict[str, tuple[str, float, float, float]]]] = \
            defaultdict(lambda: defaultdict(dict))
        for hit in merged:
            if hit.identity_pct < MIN_ALLELE_IDENTITY or hit.coverage_pct < MIN_ALLELE_COVERAGE:
                continue
            piece = batch.id_map.get(hit.qseqid)
            if piece is None:
                continue
            module, _, rest = hit.sseqid.partition("#")
            parts = rest.split(".")
            if len(parts) < 2:
                continue
            locus, allele = parts[0], ".".join(parts[1:])
            if not self.profiles(module).table:
                # Schemes with no profile table (serotyping, pathovar markers)
                # keep every allele in one FASTA, so the file stem is the same
                # for all of them. Grouping on it would keep a single hit and
                # make, for instance, the O antigen impossible to ever report.
                locus = _marker_group(allele)
            sample = piece.sample
            candidate = (allele, hit.identity_pct, hit.coverage_pct, hit.bitscore)
            current = found[sample][module].get(locus)
            if current is None or _better(candidate, current):
                found[sample][module][locus] = candidate

        results: dict[str, list[TypingResult]] = {}
        species_by_sample = species_by_sample or {}
        for sample in sorted(batch.samples()):
            modules = self.applicable_modules(species_by_sample.get(sample))
            results[sample] = self._call_sample(found.get(sample, {}), modules)
        return results

    def _call_sample(self, by_module: dict[str, dict[str, tuple]],
                     modules: list[str]) -> list[TypingResult]:
        out: list[TypingResult] = []
        for module in modules:
            hits = by_module.get(module, {})
            info = MODULE_INFO.get(module, {})
            label = info.get("label", module.split("__")[-1])
            profile = self.profiles(module)
            scheme_loci = profile.loci or sorted(hits)
            if not hits:
                out.append(TypingResult(scheme=label, call="-", lineage="-",
                                        loci_found=0, loci_total=len(scheme_loci)))
                continue
            alleles: dict[str, str] = {}
            exact = 0
            for locus in scheme_loci:
                hit = hits.get(locus)
                if hit is None:
                    alleles[locus] = "-"
                elif hit[1] >= 100.0 and hit[2] >= 100.0:
                    alleles[locus] = hit[0]
                    exact += 1
                else:
                    alleles[locus] = f"{hit[0]}*"
            call = "-"
            lineage = "-"
            note = ""
            if profile.table:
                key = tuple(alleles[locus] for locus in scheme_loci)
                st = profile.table.get(key)
                if st is not None:
                    call = st
                    lineage = profile.lineage.get(st, "-")
                elif exact or any(v != "-" for v in alleles.values()):
                    call = "NA"
                    note = "locus present, allele combination not in profile table"
                    lineage = _majority_lineage(profile, alleles, scheme_loci)
            elif hits:
                # No profile table: report every marker or antigen found, best
                # allele per group. A single-locus scheme such as wzi has one
                # group and so still reports just its allele.
                groups = sorted(hits.items(), key=lambda kv: _marker_sort_key(kv[0]))
                antigens = [key for key, _v in groups if _ANTIGEN_NAME.fullmatch(key)]
                if antigens and len(antigens) == len(groups):
                    # Conventional serotype notation, e.g. O121:H7.
                    call = ":".join(antigens)
                else:
                    call = "/".join(key if key != value[0] else str(value[0])
                                    for key, value in groups)
                lineage = "-"
            present = sum(1 for value in alleles.values() if value != "-")
            out.append(TypingResult(scheme=label, call=call, lineage=lineage or "-",
                                    alleles=alleles, loci_found=present,
                                    loci_total=len(scheme_loci), note=note))
        return out


_ANTIGEN = re.compile(r"^([OH]\d+)[-_]")
_ANTIGEN_NAME = re.compile(r"[OH]\d+")


def _marker_group(allele: str) -> str:
    """The thing an allele of a profile-less scheme identifies.

    ``O88-4-wzx`` and ``O88-2-wzy`` are two genes of the same O antigen, while
    ``H7-6-fliC`` is a different antigen and has to stay separate; a pathovar
    marker such as ``ipaH_c`` is its own group.
    """
    match = _ANTIGEN.match(allele)
    if match:
        return match.group(1)
    return allele.split("-")[0]


def _marker_sort_key(group: str) -> tuple:
    """O antigens before H antigens, each in numeric order, then everything else."""
    match = re.fullmatch(r"([OH])(\d+)", group)
    if match:
        return (0 if match.group(1) == "O" else 1, int(match.group(2)), group)
    return (2, 0, group)


def _majority_lineage(profile: SchemeProfile, alleles: dict[str, str],
                      scheme_loci: list[str]) -> str:
    """Lineage for an unlisted allele combination, or ``-`` when it is ambiguous.

    Only a combination one allele away from a known profile is annotated, and
    only when every profile that close agrees on the lineage. A looser
    nearest-neighbour rule is wrong often enough to mislead: on the *E. coli*
    Achtman scheme, accepting a half-match names the wrong clonal complex for a
    third of novel profiles.
    """
    if not profile.lineage:
        return "-"
    observed = [alleles[locus].rstrip("*") for locus in scheme_loci]
    if any(value == "-" for value in observed):
        return "-"
    best_shared = -1
    nearest: set[str] = set()
    for key, st in profile.table.items():
        shared = sum(1 for a, b in zip(key, observed) if a == b)
        if shared > best_shared:
            best_shared = shared
            nearest = {profile.lineage.get(st, "")}
        elif shared == best_shared:
            nearest.add(profile.lineage.get(st, ""))
    # One differing allele at most, and no disagreement among equally close profiles.
    if best_shared < len(scheme_loci) - 1 or len(nearest) != 1:
        return "-"
    lineage = next(iter(nearest)).strip()
    return f"{lineage} (single-locus variant)" if lineage else "-"


def _better(candidate: tuple, current: tuple) -> bool:
    def rank(value: tuple) -> tuple:
        exact = value[1] >= 100.0 and value[2] >= 100.0
        return (exact, value[2], value[1], value[3])
    return rank(candidate) > rank(current)


# ---------------------------------------------------------------------- scores
def _has_class(hits: list[Hit], *, drug_class: str = "", subclass: str = "",
               gene_prefix: str = "") -> bool:
    for hit in hits:
        if hit.element_type != "AMR":
            continue
        if drug_class and drug_class.upper() in (hit.drug_class or "").upper():
            return True
        if subclass and subclass.upper() in (hit.subclass or "").upper():
            return True
        if gene_prefix and hit.gene.lower().startswith(gene_prefix.lower()):
            return True
    return False


#: Beta-lactamase families that are chromosomal, intrinsic or narrow-spectrum.
#: They are catalogued under a cephalosporin subclass but are not ESBLs, and
#: counting them would make every E. coli (which all carry blaEC) score 1.
_NOT_ESBL_FAMILIES = (
    "blaec", "ampc", "blaampc", "blaact", "blamir", "blacmy", "bladha", "blafox",
    "blamox", "blaadc", "blacfe", "blalat", "blabil", "blaoxa-1", "blaoxy",
    "blaokp", "blalen", "blasrt", "blacepa", "blapdc", "blaz",
)
#: Gene-name fragments that identify a carbapenemase when curated subclass
#: annotation is missing, as it is for databases with no AMRFinderPlus family.
_CARBAPENEMASE_HINTS = ("kpc", "ndm", "vim", "imp", "oxa-23", "oxa-24", "oxa-40",
                        "oxa-48", "oxa-51", "oxa-58", "oxa-181", "oxa-232", "ges-",
                        "spm", "gim", "sim", "bic", "dim", "aim", "smb", "tmb", "frl")
#: Same, for extended-spectrum beta-lactamases.
_ESBL_HINTS = ("ctx-m", "veb", "per", "tla", "bel", "sfo", "shv-2", "shv-5", "shv-12",
               "tem-3", "tem-10", "tem-52", "ges-1", "ges-9")


def _family_of(hit: Hit) -> str:
    return (hit.gene or "").lower()


def _strip_bla(name: str) -> str:
    """Remove a literal ``bla`` prefix. ``str.lstrip`` would remove characters."""
    return name[3:] if name.startswith("bla") else name


def _in_family(name: str, families: tuple[str, ...]) -> bool:
    """True when *name* is a member of one of *families*.

    Matching stops at a family boundary: ``blaCTX-M-15`` must not be read as a
    member of ``blaACT`` just because ``ct`` happens to start its name.
    """
    key = _strip_bla(name)
    for family in families:
        stem = _strip_bla(family)
        if key == stem or key.startswith(stem + "-"):
            return True
    return False


def resistance_score(hits: list[Hit]) -> tuple[int, dict[str, bool]]:
    """Kleborate's 0-3 resistance score, generalised via drug-class annotation.

    0 = neither ESBL nor carbapenemase; 1 = ESBL only; 2 = carbapenemase;
    3 = carbapenemase plus colistin resistance.

    Intrinsic and AmpC-type beta-lactamases are excluded from the ESBL test, and
    gene names are consulted when a database supplies no curated subclass, so the
    score does not depend on which databases were screened.
    """
    beta_lactam = [h for h in hits if h.element_type == "AMR"
                   and ("BETA-LACTAM" in (h.drug_class or "").upper()
                        or _family_of(h).startswith("bla"))]
    has_carb = any("CARBAPENEM" in (h.subclass or "").upper() for h in beta_lactam) or \
        any(hint in _family_of(h) for h in beta_lactam for hint in _CARBAPENEMASE_HINTS)

    def is_esbl(hit: Hit) -> bool:
        name = _family_of(hit)
        if _in_family(name, _NOT_ESBL_FAMILIES):
            return False
        if "CEPHALOSPORIN" in (hit.subclass or "").upper():
            return True
        return any(hint in name for hint in _ESBL_HINTS)

    has_esbl = any(is_esbl(h) for h in beta_lactam)
    has_col = _has_class(hits, drug_class="COLISTIN", subclass="COLISTIN", gene_prefix="mcr")
    if has_carb and has_col:
        score = 3
    elif has_carb:
        score = 2
    elif has_esbl:
        score = 1
    else:
        score = 0
    return score, {"esbl": has_esbl, "carbapenemase": has_carb, "colistin": has_col}


def virulence_score(species: SpeciesCall, typing: list[TypingResult]) -> tuple[int, str]:
    """Species-aware virulence score from the lineage typing results."""
    present = {t.scheme: t for t in typing}
    locus_present: set[str] = set()
    for module, info in MODULE_INFO.items():
        label = info.get("label")
        result = present.get(label)
        if result is not None and result.call not in ("-", ""):
            locus_present.add(info.get("locus", label))
    # Only genuine virulence loci count; MLST, serotype and wzi are typing
    # schemes, and including them would score an ordinary isolate for having a
    # sequence type.
    locus_present &= VIRULENCE_LOCI
    rules = VIRULENCE_RULES.get(species.genus)
    if rules is None:
        # No published rule for this species: fall back to a count of loci found.
        return min(len(locus_present), 5), "virulence loci detected"
    for score, required, forbidden in rules:
        if all(locus in locus_present for locus in required) and \
           not any(locus in locus_present for locus in forbidden):
            return score, f"{species.genus} rule"
    return 0, f"{species.genus} rule"
