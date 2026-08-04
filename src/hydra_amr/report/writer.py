"""Dispatch results to the requested output formats."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Sequence

import pandas as pd

from ..records import SampleResult
from ..utils import LOG, HydraError
from . import html as html_report
from .tables import (abricate_table, amrfinder_table, class_summary, heteroresistance_table,
                     long_table, matrix, mlst_table, summary_table, typing_table)

#: Suffixes used for each artefact, appended to the output prefix.
ARTEFACTS = {
    "long": "",
    "summary": ".summary",
    "matrix": ".matrix",
    "classes": ".classes",
    "mlst": ".mlst",
    "typing": ".typing",
    "hetero": ".heteroresistance",
    "abricate": ".abricate",
    "amrfinder": ".amrfinder",
}


def _write_frame(frame: pd.DataFrame, path: Path, fmt: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fmt == "csv":
        frame.to_csv(path, index=frame.index.name is not None)
    elif fmt == "tsv":
        frame.to_csv(path, sep="\t", index=frame.index.name is not None)
    elif fmt == "xlsx":
        try:
            frame.to_excel(path, index=frame.index.name is not None)
        except ImportError as exc:
            raise HydraError("xlsx output needs openpyxl.\n"
                             "  conda install -c conda-forge openpyxl") from exc
    else:
        raise HydraError(f"cannot write a table as '{fmt}'")
    LOG.info("wrote %s", path)


def write_outputs(results: Sequence[SampleResult], *, outdir: Path | None, prefix: str = "hydra",
                  formats: Sequence[str] = ("tsv",), cell: str = "binary",
                  matrix_rows: str = "sample", matrix_columns: str = "gene",
                  element_types: Sequence[str] | None = None,
                  databases: Sequence[str] = (), command: str = "",
                  meta: dict | None = None, to_stdout: bool = False,
                  title: str = "Hydra report") -> list[Path]:
    """Write every requested artefact; returns the paths created."""
    written: list[Path] = []
    long = long_table(results)

    if to_stdout:
        long.to_csv(sys.stdout, sep="\t", index=False)
        ignored = [f for f in formats if f not in ("tsv",)]
        if ignored:
            LOG.warning("--stdout writes only the long table; ignoring --format %s. "
                        "Drop --stdout and pass --outdir to get those files.",
                        ", ".join(ignored))
        return written

    if outdir is None:
        raise HydraError("no output directory set; pass --outdir or --stdout")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    base = outdir / prefix

    table_formats = [f for f in formats if f in ("tsv", "csv", "xlsx")]
    compat_only = {"abricate", "amrfinder"} & set(formats)
    if not table_formats and not (set(formats) & {"json", "html"}) and not compat_only:
        table_formats = ["tsv"]

    frames = {
        "long": long,
        "summary": summary_table(results),
        "matrix": matrix(results, rows=matrix_rows, columns=matrix_columns, cell=cell,
                         element_types=element_types),
        "classes": class_summary(results),
        "mlst": mlst_table(results),
        "typing": typing_table(results),
        "hetero": heteroresistance_table(results),
    }

    for fmt in table_formats:
        for name, frame in frames.items():
            # Always emit the long and summary tables, even empty, so downstream
            # tooling still sees a header row.
            if (frame is None or frame.empty) and name not in ("long", "summary"):
                continue
            path = Path(f"{base}{ARTEFACTS[name]}.{fmt}")
            _write_frame(frame, path, fmt)
            written.append(path)

    # The compatibility layouts are defined as tab-separated, whatever the other
    # table formats are, and are written whenever they are requested.
    for name, builder in (("abricate", abricate_table), ("amrfinder", amrfinder_table)):
        if name in formats:
            path = Path(f"{base}{ARTEFACTS[name]}.tsv")
            _write_frame(builder(results), path, "tsv")
            written.append(path)

    if "json" in formats:
        path = Path(f"{base}.json")
        payload = {
            "hydra_version": meta.get("hydra_version") if meta else None,
            "command": command,
            "databases": list(databases),
            "parameters": meta or {},
            "samples": [result.as_dict() for result in results],
        }
        with open(path, "w") as handle:
            json.dump(payload, handle, indent=2, default=str)
        LOG.info("wrote %s", path)
        written.append(path)

    if "html" in formats:
        path = Path(f"{base}.html")
        html_report.write_html(path, results, cell=cell, title=title, databases=databases,
                               command=command, extra=meta, element_types=element_types,
                               matrix_rows=matrix_rows, matrix_columns=matrix_columns)
        written.append(path)

    return written
