"""Named option bundles so common jobs need one flag instead of ten.

A preset only sets defaults: anything given explicitly on the command line wins.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .utils import HydraError


@dataclass(frozen=True)
class Preset:
    name: str
    summary: str
    #: Values applied to the parsed namespace when the user did not set them.
    options: dict = field(default_factory=dict)
    detail: str = ""


PRESETS: dict[str, Preset] = {}


def _add(preset: Preset) -> None:
    PRESETS[preset.name] = preset


_add(Preset(
    name="fast",
    summary="quickest useful screen: NCBI genes only, no translated search or typing",
    options={"db": ["ncbi"], "mlst": False, "typing": False, "protein": False,
             "point_mutations": False, "format": ["tsv"]},
    detail="Use for large batches when only acquired AMR genes are needed.",
))
_add(Preset(
    name="standard",
    summary="the default: NCBI + VFDB genes, translated AMR search, point mutations, MLST, typing",
    options={"db": ["ncbi", "vfdb", "amrfinderplus"], "mlst": True, "typing": True,
             "protein": True, "point_mutations": True, "format": ["tsv", "html"]},
    detail="Balanced profile suitable for routine isolate characterisation.",
))
_add(Preset(
    name="deep",
    summary="everything: all nucleotide databases, --plus elements, typing and mutations",
    options={"db": ["all"], "plus": True, "mlst": True, "typing": True, "protein": True,
             "point_mutations": True, "format": ["tsv", "json", "html"]},
    detail="Slowest and most complete; expect redundant calls across databases.",
))
_add(Preset(
    name="surveillance",
    summary="multi-sample recap: presence/absence matrices, class summary and HTML report",
    options={"db": ["ncbi", "vfdb", "amrfinderplus"], "mlst": True, "typing": True,
             "protein": True, "point_mutations": True,
             "format": ["tsv", "csv", "html", "json"], "cell": "binary"},
    detail="Adds every tabular artefact plus the clustered HTML heatmaps.",
))
_add(Preset(
    name="amr",
    summary="resistance only: NCBI, CARD and ResFinder plus point mutations",
    options={"db": ["ncbi", "card", "resfinder", "amrfinderplus"], "typing": False,
             "protein": True, "point_mutations": True, "format": ["tsv", "html"]},
))
_add(Preset(
    name="virulence",
    summary="virulence only: VFDB and the lineage typing schemes",
    options={"db": ["vfdb", "ecoli_vf"], "typing": True, "protein": False,
             "point_mutations": False, "format": ["tsv", "html"]},
))
_add(Preset(
    name="linezolid",
    summary="linezolid heteroresistance from reads: 23S allele fractions plus cfr/optrA/poxtA",
    options={"db": ["ncbi", "amrfinderplus"], "mlst": True, "typing": False, "protein": True,
             "point_mutations": True, "heteroresistance": True,
             "min_allele_fraction": 0.02, "min_allele_reads": 3, "min_depth": 20,
             "format": ["tsv", "html", "json"]},
    detail="Needs paired reads. Reports the fraction of rRNA operons carrying each "
           "catalogued 23S mutation, so minority resistant alleles invisible to an "
           "assembly consensus are still called.",
))
_add(Preset(
    name="gram-positive",
    summary="staphylococci, enterococci and streptococci: AMR, mutations and heteroresistance",
    options={"db": ["ncbi", "vfdb", "amrfinderplus"], "mlst": True, "typing": True,
             "protein": True, "point_mutations": True, "heteroresistance": True,
             "format": ["tsv", "html"]},
))
_add(Preset(
    name="enterobacterales",
    summary="Klebsiella/E. coli focus: AMR, virulence loci, lineage typing and scores",
    options={"db": ["ncbi", "vfdb", "plasmidfinder", "amrfinderplus"], "mlst": True,
             "typing": True, "protein": True, "point_mutations": True, "plus": True,
             "format": ["tsv", "html", "json"]},
))
_add(Preset(
    name="abricate",
    summary="abricate-compatible behaviour and output columns",
    options={"db": ["ncbi"], "mlst": False, "typing": False, "protein": False,
             "point_mutations": False, "min_identity": 80.0, "min_coverage": 0.0,
             "format": ["abricate"]},
    detail="Matches abricate's thresholds and emits its exact column layout.",
))
_add(Preset(
    name="amrfinder",
    summary="AMRFinderPlus-compatible behaviour and output columns",
    options={"db": ["amrfinderplus"], "mlst": True, "typing": False, "protein": True,
             "point_mutations": True, "plus": True, "format": ["amrfinder"]},
    detail="Translated search with organism-aware point mutations, reported in "
           "AMRFinderPlus's column layout.",
))
_add(Preset(
    name="plasmid",
    summary="plasmid replicon typing and the AMR genes carried alongside",
    options={"db": ["plasmidfinder", "ncbi"], "mlst": False, "typing": False,
             "protein": False, "point_mutations": False, "format": ["tsv", "html"]},
))


def get(name: str) -> Preset:
    key = name.strip().lower()
    if key not in PRESETS:
        known = ", ".join(sorted(PRESETS))
        raise HydraError(f"unknown preset '{name}'. Available presets: {known}")
    return PRESETS[key]


def describe() -> str:
    width = max(len(name) for name in PRESETS)
    lines = []
    for name, preset in PRESETS.items():
        lines.append(f"  {name:<{width}}  {preset.summary}")
        if preset.detail:
            lines.append(f"  {'':<{width}}  {preset.detail}")
    return "\n".join(lines)
