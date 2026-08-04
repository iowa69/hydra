"""Tabular views over a set of :class:`~hydra.records.SampleResult`."""

from __future__ import annotations

from typing import Sequence

import pandas as pd

from ..config import CELL_MODES
from ..records import HIT_COLUMNS, SampleResult
from ..utils import HydraError

#: Columns of the abricate-compatible output, in abricate's own order.
ABRICATE_COLUMNS = ("#FILE", "SEQUENCE", "START", "END", "STRAND", "GENE", "COVERAGE",
                    "COVERAGE_MAP", "GAPS", "%COVERAGE", "%IDENTITY", "DATABASE",
                    "ACCESSION", "PRODUCT", "RESISTANCE")

#: Columns of the AMRFinderPlus-compatible output.
AMRFINDER_COLUMNS = ("Name", "Contig id", "Start", "Stop", "Strand", "Element symbol",
                     "Element name", "Scope", "Element type", "Element subtype",
                     "Class", "Subclass", "Method", "Target length",
                     "Reference sequence length", "% Coverage of reference",
                     "% Identity to reference", "Accession of closest sequence")


def long_table(results: Sequence[SampleResult]) -> pd.DataFrame:
    """One row per detected element."""
    rows = [hit.as_row() for result in results for hit in result.hits]
    if not rows:
        return pd.DataFrame(columns=list(HIT_COLUMNS))
    frame = pd.DataFrame(rows)
    for column in HIT_COLUMNS:
        if column not in frame.columns:
            frame[column] = None
    return frame[list(HIT_COLUMNS)].sort_values(
        ["sample", "element_type", "class", "gene"], kind="stable").reset_index(drop=True)


def summary_table(results: Sequence[SampleResult]) -> pd.DataFrame:
    """One row per sample: species, ST, typing calls, counts and QC."""
    rows = [result.summary_row() for result in results]
    if not rows:
        return pd.DataFrame(columns=["sample"])
    frame = pd.DataFrame(rows)
    lead = [c for c in ("sample", "input_type", "species", "species_confidence",
                        "mlst_scheme", "ST") if c in frame.columns]
    rest = [c for c in frame.columns if c not in lead]
    return frame[lead + rest]


def mlst_table(results: Sequence[SampleResult]) -> pd.DataFrame:
    """Per-sample MLST detail including the individual allele calls."""
    rows = []
    for result in results:
        row = {
            "sample": result.sample,
            "scheme": result.mlst.scheme,
            "ST": result.mlst.sequence_type,
            "loci_exact": result.mlst.loci_found,
            "loci_total": result.mlst.loci_total,
            "note": result.mlst.note,
        }
        for locus, allele in result.mlst.alleles.items():
            row[locus] = allele
        rows.append(row)
    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=["sample", "scheme", "ST"])


def typing_table(results: Sequence[SampleResult]) -> pd.DataFrame:
    """Per-sample lineage typing calls."""
    rows = []
    for result in results:
        row = {"sample": result.sample, "species": result.species.name}
        for typing in result.typing:
            row[typing.scheme] = typing.call
            row[f"{typing.scheme}_lineage"] = typing.lineage
        row.update(result.scores)
        rows.append(row)
    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=["sample"])


def _cell_value(hits: list, mode: str):
    """Reduce the hits sharing one matrix cell to a single value."""
    if not hits:
        return None
    if mode == "binary":
        return 1
    if mode == "count":
        return len(hits)
    if mode == "genes":
        return len({h.gene for h in hits})
    if mode == "identity":
        return round(max(h.identity_pct for h in hits), 2)
    if mode == "coverage":
        return round(max(h.coverage_pct for h in hits), 2)
    if mode == "depth":
        depths = [h.depth for h in hits if h.depth is not None]
        return round(max(depths), 2) if depths else None
    if mode == "fraction":
        fractions = [h.allele_fraction for h in hits if h.allele_fraction is not None]
        return round(max(fractions), 4) if fractions else None
    if mode == "symbol":
        best = max(hits, key=lambda h: (h.identity_pct, h.coverage_pct))
        return f"{best.identity_pct:.0f}/{best.coverage_pct:.0f}"
    raise HydraError(f"unknown cell mode '{mode}'; choose from {', '.join(CELL_MODES)}")


def matrix(results: Sequence[SampleResult], *, rows: str = "sample", columns: str = "gene",
           cell: str = "binary", element_types: Sequence[str] | None = None,
           databases: Sequence[str] | None = None, absent=0,
           drop_empty: bool = True, primary_only: bool = True) -> pd.DataFrame:
    """Pivot the hits into a rows x columns matrix.

    ``rows`` and ``columns`` accept ``sample``, ``gene``, ``class``, ``subclass``,
    ``database``, ``element_type`` or ``product``. Redundant cross-database calls
    for the same locus are excluded unless *primary_only* is turned off, so a
    gene present once is counted once however many databases reported it.
    """
    if cell not in CELL_MODES:
        raise HydraError(f"unknown --cell '{cell}'; choose from {', '.join(CELL_MODES)}")
    field_map = {
        "sample": lambda h: h.sample,
        "gene": lambda h: h.gene,
        "class": lambda h: h.drug_class or "unclassified",
        "subclass": lambda h: h.subclass or "unclassified",
        "database": lambda h: h.database,
        "element_type": lambda h: h.element_type,
        "product": lambda h: h.product or h.gene,
    }
    for axis, name in (("rows", rows), ("columns", columns)):
        if name not in field_map:
            raise HydraError(f"--{axis} must be one of {', '.join(sorted(field_map))} "
                             f"(got '{name}')")
    row_key = field_map[rows]
    col_key = field_map[columns]
    # A pivot keyed on database exists to compare databases, so it must keep
    # every database's own view of a locus.
    drop_redundant = primary_only and "database" not in (rows, columns)

    wanted_types = set(element_types) if element_types else None
    wanted_dbs = set(databases) if databases else None
    buckets: dict[tuple, list] = {}
    all_rows: list[str] = []
    all_cols: list[str] = []
    for result in results:
        label = result.sample
        if rows == "sample" and label not in all_rows:
            all_rows.append(label)
        for hit in result.hits:
            if drop_redundant and not hit.primary:
                continue
            if wanted_types and hit.element_type not in wanted_types:
                continue
            if wanted_dbs and hit.database not in wanted_dbs:
                continue
            r = row_key(hit)
            c = col_key(hit)
            if r not in all_rows:
                all_rows.append(r)
            if c not in all_cols:
                all_cols.append(c)
            buckets.setdefault((r, c), []).append(hit)

    if not all_cols and drop_empty:
        return pd.DataFrame(index=pd.Index(all_rows, name=rows))
    data = {}
    for c in sorted(all_cols):
        column_values = []
        for r in all_rows:
            value = _cell_value(buckets.get((r, c), []), cell)
            column_values.append(absent if value is None else value)
        data[c] = column_values
    frame = pd.DataFrame(data, index=pd.Index(all_rows, name=rows))
    frame.columns.name = columns
    return frame


def _coverage_glyph(hit, width: int = 15) -> str:
    """abricate's COVERAGE_MAP: ``=`` where the reference is covered, ``.`` where not."""
    spans, _, total = hit.coverage.rpartition("/")
    try:
        length = int(total)
    except ValueError:
        return ""
    if length <= 0:
        return ""
    covered = [False] * width
    for span in spans.split(","):
        start, _, end = span.partition("-")
        try:
            lo, hi = int(start), int(end or start)
        except ValueError:
            continue
        for cell in range(width):
            position = (cell * length) // width + 1
            if lo <= position <= hi:
                covered[cell] = True
    return "".join("=" if flag else "." for flag in covered)


def abricate_table(results: Sequence[SampleResult]) -> pd.DataFrame:
    """Drop-in replacement for abricate's tabular output."""
    rows = []
    for result in results:
        source = result.inputs[0] if result.inputs else result.sample
        for hit in result.hits:
            resistance = hit.subclass or hit.drug_class or ""
            # abricate puts the full "start-end/reflen" span in COVERAGE and an
            # alignment glyph in COVERAGE_MAP; scripts parse the former.
            rows.append({
                "#FILE": source, "SEQUENCE": hit.sequence, "START": hit.start, "END": hit.end,
                "STRAND": hit.strand, "GENE": hit.gene,
                "COVERAGE": hit.coverage,
                "COVERAGE_MAP": _coverage_glyph(hit), "GAPS": hit.gaps,
                "%COVERAGE": round(hit.coverage_pct, 2), "%IDENTITY": round(hit.identity_pct, 2),
                "DATABASE": hit.database, "ACCESSION": hit.accession, "PRODUCT": hit.product,
                "RESISTANCE": resistance,
            })
    return pd.DataFrame(rows, columns=list(ABRICATE_COLUMNS))


def amrfinder_table(results: Sequence[SampleResult]) -> pd.DataFrame:
    """Output shaped like AMRFinderPlus's report."""
    rows = []
    for result in results:
        for hit in result.hits:
            reference_length = ""
            if hit.coverage and "/" in hit.coverage:
                reference_length = hit.coverage.rsplit("/", 1)[1]
            rows.append({
                "Name": hit.sample, "Contig id": hit.sequence, "Start": hit.start,
                "Stop": hit.end, "Strand": hit.strand, "Element symbol": hit.gene,
                "Element name": hit.product, "Scope": "core" if hit.database == "ncbi" else "plus",
                "Element type": hit.element_type, "Element subtype": hit.element_subtype,
                "Class": hit.drug_class, "Subclass": hit.subclass, "Method": hit.method,
                "Target length": abs(hit.end - hit.start) + 1,
                "Reference sequence length": reference_length,
                "% Coverage of reference": round(hit.coverage_pct, 2),
                "% Identity to reference": round(hit.identity_pct, 2),
                "Accession of closest sequence": hit.accession,
            })
    return pd.DataFrame(rows, columns=list(AMRFINDER_COLUMNS))


def heteroresistance_table(results: Sequence[SampleResult]) -> pd.DataFrame:
    """Point mutations called from reads, with their allele fractions."""
    rows = []
    for result in results:
        for hit in result.hits:
            if hit.allele_fraction is None:
                continue
            status = "heteroresistant" if "HETERORESISTANT" in hit.note.upper() else (
                "fixed" if "FIXED" in hit.note.upper() else "detected")
            rows.append({
                "sample": result.sample, "gene": hit.gene, "mutation": hit.note.split(";")[0],
                "class": hit.drug_class, "subclass": hit.subclass,
                "allele_fraction": round(hit.allele_fraction, 4),
                "depth": hit.depth, "status": status, "reference": hit.sequence,
                "position": hit.start, "detail": hit.note,
            })
    return pd.DataFrame(rows, columns=["sample", "gene", "mutation", "class", "subclass",
                                       "allele_fraction", "depth", "status", "reference",
                                       "position", "detail"])


def class_summary(results: Sequence[SampleResult]) -> pd.DataFrame:
    """Sample x drug-class matrix of distinct resistance genes.

    Distinct genes rather than hits: a gene present on three replicons is one
    resistance determinant, not three.
    """
    return matrix(results, rows="sample", columns="class", cell="genes",
                  element_types=["AMR"], absent=0)
