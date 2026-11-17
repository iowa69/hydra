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
) + (
    # What exactly mutated, for mutation and disruption calls (engines/lesions.py);
    # empty for every other hit.
    "mutation", "mutation_type", "protein_change", "protein_change_hgvs",
    "nucleotide_change", "nucleotide_change_basis", "observed_codon", "reference_codon",
    "reference_codon_inferred", "base_change", "mutation_location", "gene_location",
    "intact_reference", "kleborate_call", "resistance", "mutation_description",
    "unresolved",
)


@dataclass
class MutationDetail:
    """What exactly changed at a mutation or disruption call, read from the isolate.

    Filled by :mod:`hydra_amr.engines.lesions` after the call is made; it never
    affects which calls are made. Field by field:

    ``symbol``             AMRFinderPlus-style name (``gyrA_S83I``, ``ompK35_G134AfsTer11``)
    ``call_type``          machine-readable kind (substitution, nonsense, insertion, deletion,
                           rRNA, promoter, disruption_frameshift, contig_end, divergent ...)
    ``type_label``         the same in words
    ``protein_change``     short form (``S83I``); ``protein_hgvs`` long form (``p.Ser83Ile``)
    ``nucleotide_change``  HGVS-style (``c.248C>T``, ``n.2059A>G``)
    ``nucleotide_basis``   where the reference base comes from: the reference sequence,
                           or "inferred" by parsimony when the database has none
    ``observed``           codon(s) or base read from the assembly
    ``reference``          reference codon/base taken from a reference sequence
    ``reference_inferred`` reference codon inferred by parsimony (never presented as fact)
    ``base_change``        which base changed (``C>T at codon position 2``)
    ``contig/start/end/strand``, ``location``  the mutated codon/base on the contig
    ``gene_location``      the gene's span on the contig
    ``intact``             disruptions: how much of the reference is intact
    ``kleborate``          Kleborate v3's name for the same change, for the loci it types
    ``resistance``         the resistance the call is associated with, in words
    ``description``        one plain sentence
    ``unresolved``         anything that could not be resolved, said plainly
    """

    symbol: str = ""
    call_type: str = ""
    type_label: str = ""
    protein_change: str = ""
    protein_hgvs: str = ""
    nucleotide_change: str = ""
    nucleotide_basis: str = ""
    observed: str = ""
    reference: str = ""
    reference_inferred: str = ""
    base_change: str = ""
    contig: str = ""
    start: int = 0
    end: int = 0
    strand: str = ""
    location: str = ""
    gene_location: str = ""
    intact: str = ""
    kleborate: str = ""
    resistance: str = ""
    description: str = ""
    unresolved: str = ""

    def as_row(self) -> dict[str, Any]:
        return {
            "mutation": self.symbol, "mutation_type": self.type_label,
            "protein_change": self.protein_change, "protein_change_hgvs": self.protein_hgvs,
            "nucleotide_change": self.nucleotide_change,
            "nucleotide_change_basis": self.nucleotide_basis,
            "observed_codon": self.observed, "reference_codon": self.reference,
            "reference_codon_inferred": self.reference_inferred,
            "base_change": self.base_change, "mutation_location": self.location,
            "gene_location": self.gene_location, "intact_reference": self.intact,
            "kleborate_call": self.kleborate, "resistance": self.resistance,
            "mutation_description": self.description, "unresolved": self.unresolved,
        }


_NO_MUTATION = {key: "" for key in MutationDetail().as_row()}


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
    #: What exactly mutated (mutation and disruption calls only); see MutationDetail.
    mutation: MutationDetail | None = None
    #: The alignment a mutation call was made from, kept only until
    #: engines/lesions.py has described it. Never written out.
    evidence: dict | None = field(default=None, repr=False, compare=False)

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
            **(self.mutation.as_row() if self.mutation is not None else _NO_MUTATION),
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
    #: Readable heading for the scheme, set by the organism typing registry.
    #: Presentation only; the scheme identifier stays the stable key.
    title: str = ""

    def as_row(self) -> dict[str, Any]:
        return {
            f"{self.scheme}": self.call,
            f"{self.scheme}_lineage": self.lineage,
        }

    @property
    def heading(self) -> str:
        """What a reader should see as the name of this scheme."""
        return self.title or self.scheme

    @property
    def loci(self) -> str:
        """``6/7``, or ``-`` when the module counts no loci for this scheme."""
        return f"{self.loci_found}/{self.loci_total}" if self.loci_total else "-"

    def evidence(self) -> str:
        """The ``alleles`` mapping as one readable line.

        This is where the module put the alleles, residues, identities,
        coordinates and database versions the call rests on. Nothing is
        abbreviated: a call a clinical microbiologist cannot check is a call
        they cannot use.
        """
        return "; ".join(f"{key}={value}" for key, value in self.alleles.items())

    def as_detail_row(self, sample: str = "", species: str = "") -> dict[str, Any]:
        """One row of the long-form typing table: the call and all its evidence."""
        row: dict[str, Any] = {}
        if sample:
            row["sample"] = sample
        if species:
            row["species"] = species
        row.update({
            "scheme": self.scheme,
            "scheme_name": self.heading,
            "call": self.call,
            "lineage": self.lineage,
            "loci": self.loci,
            "evidence": self.evidence(),
            "caveats": self.note,
        })
        return row


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
