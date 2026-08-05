"""The record types every engine produces and every reporter consumes."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

#: Column order for the canonical long-format table.
HIT_COLUMNS = (
    "sample", "database", "element_type", "element_subtype", "gene", "accession",
    "product", "class", "subclass", "sequence", "start", "end", "strand",
    "coverage", "coverage_pct", "identity_pct", "gaps", "depth", "allele_fraction",
    "method", "resolution", "primary", "note",
)


@dataclass
class Hit:
    """One detected element in one sample.

    ``method`` follows the AMRFinderPlus convention (EXACTP, BLASTP, PARTIALP,
    POINTP, POINTN, BLASTN, PARTIALN, READS ...) so downstream users can filter
    on call confidence.
    """

    sample: str
    database: str
    gene: str
    accession: str = ""
    product: str = ""
    element_type: str = "AMR"
    element_subtype: str = "AMR"
    drug_class: str = ""
    subclass: str = ""
    sequence: str = ""          # contig / reference name the hit sits on
    start: int = 0
    end: int = 0
    strand: str = "+"
    coverage: str = ""          # human-readable span on the reference, e.g. 1-861/861
    coverage_pct: float = 0.0
    identity_pct: float = 0.0
    gaps: int = 0
    bitscore: float = 0.0
    depth: float | None = None            # read mode
    allele_fraction: float | None = None  # read mode / heteroresistance
    method: str = "BLASTN"
    resolution: str = "COMPLETE"          # COMPLETE | PARTIAL | POINT | INTERNAL_STOP
    #: False when another database already reported this same locus. Redundant
    #: hits stay in the long table but are excluded from counts and matrices.
    primary: bool = True
    note: str = ""

    def as_row(self) -> dict[str, Any]:
        return {
            "sample": self.sample,
            "database": self.database,
            "element_type": self.element_type,
            "element_subtype": self.element_subtype,
            "gene": self.gene,
            "accession": self.accession,
            "product": self.product,
            "class": self.drug_class,
            "subclass": self.subclass,
            "sequence": self.sequence,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "coverage": self.coverage,
            "coverage_pct": round(self.coverage_pct, 2),
            "identity_pct": round(self.identity_pct, 2),
            "gaps": self.gaps,
            "depth": None if self.depth is None else round(self.depth, 2),
            "allele_fraction": None if self.allele_fraction is None else round(self.allele_fraction, 4),
            "method": self.method,
            "resolution": self.resolution,
            "primary": self.primary,
            "note": self.note,
        }


@dataclass
class MlstCall:
    """MLST result for one sample."""

    scheme: str = "-"
    sequence_type: str = "-"
    alleles: dict[str, str] = field(default_factory=dict)
    similarity: float = 0.0     # mean identity across matched loci
    loci_found: int = 0
    loci_total: int = 0
    novel_alleles: list[str] = field(default_factory=list)
    note: str = ""
    candidates: list[dict] = field(default_factory=list)
    #: What the call was made from: an assembly, or reads mapped to the loci.
    source: str = "assembly"

    @property
    def profile(self) -> str:
        if not self.alleles:
            return "-"
        return " ".join(f"{locus}({allele})" for locus, allele in sorted(self.alleles.items()))

    def as_row(self) -> dict[str, Any]:
        return {
            "mlst_scheme": self.scheme,
            "ST": self.sequence_type,
            "mlst_source": self.source if self.scheme != "-" else "-",
            "mlst_profile": self.profile,
            "mlst_loci_found": f"{self.loci_found}/{self.loci_total}" if self.loci_total else "-",
            "mlst_note": self.note,
        }


@dataclass
class SpeciesCall:
    """Species identification for one sample."""

    name: str = "unknown"
    genus: str = ""
    species: str = ""
    confidence: str = "none"    # strong | good | weak | none
    evidence: str = ""
    distance: float | None = None
    #: AMRFinderPlus taxgroup this maps to, when one exists.
    organism: str | None = None

    def as_row(self) -> dict[str, Any]:
        return {
            "species": self.name,
            "species_confidence": self.confidence,
            "species_evidence": self.evidence,
            "organism_db": self.organism or "-",
        }


@dataclass
class TypingResult:
    """One lineage or sublineage typing call."""

    scheme: str
    call: str = "-"
    lineage: str = "-"
    alleles: dict[str, str] = field(default_factory=dict)
    loci_found: int = 0
    loci_total: int = 0
    note: str = ""

    def as_row(self) -> dict[str, Any]:
        return {
            f"{self.scheme}": self.call,
            f"{self.scheme}_lineage": self.lineage,
        }


@dataclass
class SampleResult:
    """Everything Hydra determined about one sample."""

    sample: str
    input_type: str = "assembly"     # assembly | reads | assembly+reads
    inputs: list[str] = field(default_factory=list)
    hits: list[Hit] = field(default_factory=list)
    mlst: MlstCall = field(default_factory=MlstCall)
    species: SpeciesCall = field(default_factory=SpeciesCall)
    typing: list[TypingResult] = field(default_factory=list)
    scores: dict[str, Any] = field(default_factory=dict)
    qc: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    #: Wall-clock of the whole run this sample belonged to, not this sample's own
    #: cost. Assemblies are screened as one batch -- a single BLAST pass per
    #: database over every contig -- so there is no per-sample time to report, and
    #: dividing the total by the sample count would invent one.
    runtime_seconds: float = 0.0

    def hits_of(self, element_type: str, primary_only: bool = True) -> list[Hit]:
        return [h for h in self.hits if h.element_type == element_type
                and (h.primary or not primary_only)]

    def gene_set(self, element_type: str | None = None, primary_only: bool = True) -> set[str]:
        return {h.gene for h in self.hits
                if (element_type is None or h.element_type == element_type)
                and (h.primary or not primary_only)}

    def class_set(self, primary_only: bool = True) -> set[str]:
        return {h.drug_class for h in self.hits
                if h.drug_class and h.element_type == "AMR" and (h.primary or not primary_only)}

    def summary_row(self) -> dict[str, Any]:
        row: dict[str, Any] = {"sample": self.sample, "input_type": self.input_type}
        row.update(self.species.as_row())
        row.update(self.mlst.as_row())
        for typing in self.typing:
            row.update(typing.as_row())
        primary = [h for h in self.hits if h.primary]
        row.update({
            "amr_genes": len({h.gene for h in primary if h.element_type == "AMR"}),
            "amr_classes": len(self.class_set()),
            "virulence_genes": len({h.gene for h in primary if h.element_type == "VIRULENCE"}),
            "stress_genes": len({h.gene for h in primary if h.element_type == "STRESS"}),
            "plasmid_replicons": len({h.gene for h in primary if h.element_type == "PLASMID"}),
            "point_mutations": len([h for h in primary if h.resolution == "POINT"
                                    and h.element_type == "AMR" and h.method != "VARIANTR"]),
            # Only catalogued resistance mutations count: a read-derived variant
            # in some other gene is an observation, not a resistance call.
            "heteroresistant_sites": len([h for h in self.hits if h.method == "POINTR"
                                          and h.element_type == "AMR"
                                          and "HETERORESISTANT" in h.note.upper()]),
            "gene_variants": len([h for h in self.hits if h.method in ("VARIANTR", "ALLELER")]),
        })
        row.update(self.scores)
        row.update({f"qc_{k}": v for k, v in self.qc.items()})
        row["warnings"] = "; ".join(self.warnings)
        return row

    def as_dict(self) -> dict[str, Any]:
        return {
            "sample": self.sample,
            "input_type": self.input_type,
            "inputs": self.inputs,
            "species": asdict(self.species),
            "mlst": asdict(self.mlst),
            "typing": [asdict(t) for t in self.typing],
            "scores": self.scores,
            "qc": self.qc,
            "hits": [h.as_row() for h in self.hits],
            "warnings": self.warnings,
            "runtime_seconds": round(self.runtime_seconds, 2),
        }
