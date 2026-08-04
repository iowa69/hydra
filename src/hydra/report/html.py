"""Self-contained HTML report with heatmaps and recaps.

Everything is inlined - no CDN, no network access at view time - so a report can
be attached to an email or archived with the run. Rows and columns of the
heatmaps are ordered by hierarchical clustering when SciPy is available and by
prevalence otherwise.
"""

from __future__ import annotations

import html
import json
from datetime import datetime
from pathlib import Path
from typing import Sequence

import pandas as pd

from .. import __version__
from ..records import SampleResult
from ..utils import LOG
from .tables import (class_summary, heteroresistance_table, long_table, matrix,
                     mlst_table, summary_table, typing_table)

_CSS = """
:root{
  --bg:#ffffff; --fg:#16181d; --muted:#5f6672; --line:#e3e6ec; --panel:#f7f8fa;
  --accent:#2563eb; --warn:#b45309; --danger:#b91c1c; --ok:#15803d;
  --c0:#f2f5fa; --c1:#dbe6f6; --c2:#b9d0ee; --c3:#8fb4e2; --c4:#6394d4; --c5:#3a72c0; --c6:#1e4f96;
  --hetero:#c2410c;
}
@media (prefers-color-scheme: dark){
  :root{ --bg:#111318; --fg:#e8eaee; --muted:#9aa3b2; --line:#272b34; --panel:#181b22;
    --c0:#1a1d24; --c1:#1e2b3f; --c2:#223d5e; --c3:#2a5081; --c4:#3566a4; --c5:#4a80c4; --c6:#6b9de0; }
}
:root[data-theme="dark"]{
  --bg:#111318; --fg:#e8eaee; --muted:#9aa3b2; --line:#272b34; --panel:#181b22;
  --c0:#1a1d24; --c1:#1e2b3f; --c2:#223d5e; --c3:#2a5081; --c4:#3566a4; --c5:#4a80c4; --c6:#6b9de0;
}
:root[data-theme="light"]{
  --bg:#ffffff; --fg:#16181d; --muted:#5f6672; --line:#e3e6ec; --panel:#f7f8fa;
  --c0:#f2f5fa; --c1:#dbe6f6; --c2:#b9d0ee; --c3:#8fb4e2; --c4:#6394d4; --c5:#3a72c0; --c6:#1e4f96;
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
  font:14px/1.55 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
.wrap{max-width:1400px;margin:0 auto;padding:28px 20px 80px}
header{display:flex;align-items:baseline;gap:14px;flex-wrap:wrap;margin-bottom:6px}
h1{font-size:26px;margin:0;letter-spacing:-.02em}
h2{font-size:17px;margin:34px 0 10px;letter-spacing:-.01em}
h3{font-size:14px;margin:20px 0 8px;color:var(--muted);font-weight:600;
  text-transform:uppercase;letter-spacing:.06em}
.sub{color:var(--muted);font-size:13px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));gap:10px;margin:18px 0 6px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:12px 14px}
.card .n{font-size:24px;font-weight:650;letter-spacing:-.02em}
.card .l{color:var(--muted);font-size:12px;margin-top:2px}
.scroll{overflow-x:auto;border:1px solid var(--line);border-radius:10px;background:var(--panel)}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:6px 9px;text-align:left;border-bottom:1px solid var(--line);white-space:nowrap}
thead th{position:sticky;top:0;background:var(--panel);z-index:2;font-weight:600;
  border-bottom:2px solid var(--line);cursor:pointer;user-select:none}
thead th:hover{color:var(--accent)}
tbody tr:hover{background:rgba(37,99,235,.06)}
.hm{border-collapse:separate;border-spacing:0;font-size:11px}
.hm th.rowlab{position:sticky;left:0;background:var(--panel);z-index:1;text-align:left;
  max-width:220px;overflow:hidden;text-overflow:ellipsis;border-right:1px solid var(--line)}
.hm thead th{height:auto;padding:4px 2px;font-weight:500;font-size:10px;
  writing-mode:vertical-rl;transform:rotate(180deg);max-height:180px;cursor:default}
.hm td{padding:0;border:0;width:16px;height:16px;min-width:16px}
.hm td div{width:100%;height:100%;border-right:1px solid var(--bg);border-bottom:1px solid var(--bg)}
.v0{background:var(--c0)}.v1{background:var(--c1)}.v2{background:var(--c2)}
.v3{background:var(--c3)}.v4{background:var(--c4)}.v5{background:var(--c5)}.v6{background:var(--c6)}
.vh{background:var(--hetero)}
.legend{display:flex;gap:6px;align-items:center;color:var(--muted);font-size:12px;margin:8px 0 0}
.legend i{width:16px;height:12px;display:inline-block;border-radius:2px}
.pill{display:inline-block;padding:1px 7px;border-radius:999px;font-size:11px;font-weight:600}
.pill.ok{background:rgba(21,128,61,.14);color:var(--ok)}
.pill.warn{background:rgba(180,83,9,.16);color:var(--warn)}
.pill.danger{background:rgba(185,28,28,.15);color:var(--danger)}
.pill.mute{background:rgba(120,125,135,.16);color:var(--muted)}
.controls{display:flex;gap:8px;margin:10px 0;flex-wrap:wrap}
input[type=search]{padding:6px 10px;border:1px solid var(--line);border-radius:8px;
  background:var(--bg);color:var(--fg);min-width:240px;font-size:13px}
button.tgl{padding:6px 12px;border:1px solid var(--line);border-radius:8px;background:var(--bg);
  color:var(--fg);cursor:pointer;font-size:13px}
button.tgl:hover{border-color:var(--accent);color:var(--accent)}
footer{margin-top:50px;color:var(--muted);font-size:12px;border-top:1px solid var(--line);padding-top:14px}
details{margin:8px 0}
summary{cursor:pointer;color:var(--accent);font-size:13px}
.empty{color:var(--muted);padding:14px;font-style:italic}
"""

_JS = """
document.querySelectorAll('table.sortable thead th').forEach(function(th){
  th.addEventListener('click',function(){
    // cellIndex is the column within this table; a document-wide counter would
    // index past the end of every table after the first.
    var idx=th.cellIndex;
    var table=th.closest('table');
    if(!table||!table.tBodies.length) return;
    var tb=table.tBodies[0];
    var rows=Array.from(tb.rows);
    var asc=!(th.dataset.asc==='1');
    th.closest('tr').querySelectorAll('th').forEach(function(o){o.dataset.asc='';});
    th.dataset.asc=asc?'1':'0';
    rows.sort(function(a,b){
      var x=a.cells[idx]?a.cells[idx].innerText.trim():'';
      var y=b.cells[idx]?b.cells[idx].innerText.trim():'';
      if(x===''&&y!=='') return 1;
      if(y===''&&x!=='') return -1;
      var nx=parseFloat(x),ny=parseFloat(y);
      var both=!isNaN(nx)&&!isNaN(ny)&&x!==''&&y!=='';
      var c=both?(nx-ny):x.localeCompare(y,undefined,{numeric:true});
      return asc?c:-c;
    });
    rows.forEach(function(r){tb.appendChild(r);});
  });
});
document.querySelectorAll('input[data-filter]').forEach(function(inp){
  var table=document.getElementById(inp.dataset.filter);
  if(!table||!table.tBodies.length){ inp.disabled=true; return; }
  inp.addEventListener('input',function(){
    var q=inp.value.toLowerCase();
    Array.from(table.tBodies[0].rows).forEach(function(r){
      r.style.display=r.innerText.toLowerCase().indexOf(q)>-1?'':'none';
    });
  });
});
"""


def _esc(value) -> str:
    return html.escape("" if value is None else str(value))


def _bucket(value: float, vmax: float) -> int:
    """Map a value onto one of seven colour steps."""
    if value is None or value <= 0 or vmax <= 0:
        return 0
    ratio = min(1.0, value / vmax)
    return max(1, min(6, int(round(ratio * 6))))


def _format_cell(value) -> str:
    """Render one table value without losing precision or leaking ``nan``.

    Allele fractions are the reason this is not just ``f"{value:.2f}"``: a 0.4%
    minority allele is a real signal, and rounding it to two decimals would
    print it as 0.
    """
    if value is None:
        return ""
    if isinstance(value, float):
        if value != value:  # NaN
            return ""
        if value == int(value) and abs(value) < 1e15:
            return str(int(value))
        text = f"{value:.6f}".rstrip("0").rstrip(".")
        return text or "0"
    if isinstance(value, bool):
        return "yes" if value else "no"
    text = str(value)
    return "" if text in ("nan", "None", "NaT") else text


def _frame_to_html(frame: pd.DataFrame, table_id: str, sortable: bool = True,
                   max_rows: int | None = None) -> str:
    if frame is None or frame.empty:
        return '<div class="empty">nothing to show</div>'
    shown = frame if max_rows is None else frame.head(max_rows)
    classes = "sortable" if sortable else ""
    out = [f'<div class="scroll"><table id="{table_id}" class="{classes}"><thead><tr>']
    for column in shown.columns:
        out.append(f"<th>{_esc(column)}</th>")
    out.append("</tr></thead><tbody>")
    for _, row in shown.iterrows():
        out.append("<tr>")
        for value in row:
            out.append(f"<td>{_esc(_format_cell(value))}</td>")
        out.append("</tr>")
    out.append("</tbody></table></div>")
    if max_rows is not None and len(frame) > max_rows:
        dropped_samples = ""
        if "sample" in frame.columns:
            missing = set(frame["sample"]) - set(shown["sample"])
            if missing:
                dropped_samples = (f"; {len(missing)} sample(s) fall entirely beyond the "
                                   f"cut-off and are not shown here")
        out.append(f'<div class="sub">showing the first {max_rows} of {len(frame)} rows'
                   f'{dropped_samples}. The complete table is in the TSV/CSV output.</div>')
    return "".join(out)


def _order_axes(frame: pd.DataFrame) -> pd.DataFrame:
    """Cluster rows/columns when SciPy is present, else order by prevalence."""
    if frame.empty or frame.shape[0] < 3 or frame.shape[1] < 3:
        return frame
    numeric = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage  # noqa: PLC0415
        from scipy.spatial.distance import pdist  # noqa: PLC0415

        def order(values):
            distance = pdist(values, metric="jaccard" if set(values.flatten()) <= {0.0, 1.0}
                             else "euclidean")
            if not len(distance) or not distance.any():
                return list(range(values.shape[0]))
            return list(leaves_list(linkage(distance, method="average")))

        row_order = order(numeric.values)
        col_order = order(numeric.values.T)
        return frame.iloc[row_order, col_order]
    except Exception as exc:  # noqa: BLE001 - clustering is a nicety, never fatal
        LOG.debug("clustering unavailable (%s); ordering by prevalence", exc)
        col_rank = numeric.astype(bool).sum(axis=0).sort_values(ascending=False)
        row_rank = numeric.astype(bool).sum(axis=1).sort_values(ascending=False)
        return frame.loc[row_rank.index, col_rank.index]


def _heatmap(frame: pd.DataFrame, cell_mode: str, title: str, max_cols: int = 260) -> str:
    if frame is None or frame.empty or frame.shape[1] == 0:
        return '<div class="empty">no data for this heatmap</div>'
    frame = _order_axes(frame)
    truncated = False
    if frame.shape[1] > max_cols:
        keep = frame.apply(pd.to_numeric, errors="coerce").fillna(0).astype(bool).sum(axis=0)
        frame = frame[keep.sort_values(ascending=False).head(max_cols).index]
        frame = _order_axes(frame)
        truncated = True
    numeric = frame.apply(pd.to_numeric, errors="coerce")
    vmax = float(numeric.max().max()) if numeric.size and numeric.notna().any().any() else 0.0
    if cell_mode == "binary":
        vmax = 1.0
    # Some cell modes are textual (``symbol`` renders "identity/coverage"), so
    # presence has to be judged from the text rather than from a magnitude.
    textual = cell_mode == "symbol" or vmax <= 0.0
    out = ['<div class="scroll"><table class="hm"><thead><tr><th class="rowlab"></th>']
    for column in frame.columns:
        out.append(f"<th>{_esc(column)}</th>")
    out.append("</tr></thead><tbody>")
    for index, row in frame.iterrows():
        out.append(f'<tr><th class="rowlab" title="{_esc(index)}">{_esc(index)}</th>')
        for column, value in row.items():
            label = _format_cell(value)
            if textual:
                present = bool(label) and label not in ("0", "-")
                step = 4 if present else 0
            else:
                try:
                    number = float(value)
                except (TypeError, ValueError):
                    number = 0.0
                if number != number:
                    number = 0.0
                step = _bucket(number, vmax)
                if number == 0:
                    label = ""
            tip = f"{index} / {column}" + (f" = {label}" if label else " = absent")
            out.append(f'<td title="{_esc(tip)}"><div class="v{step}"></div></td>')
        out.append("</tr>")
    out.append("</tbody></table></div>")
    legend = ['<div class="legend"><span>absent</span>']
    for step in range(7):
        legend.append(f'<i class="v{step}"></i>')
    unit = {"binary": "present", "identity": "% identity", "coverage": "% coverage",
            "count": "copies", "genes": "distinct genes", "depth": "read depth",
            "fraction": "allele fraction", "symbol": "identity/coverage"}.get(
        cell_mode, cell_mode)
    scale = "" if (textual or vmax <= 1) else f" (max {vmax:g})"
    legend.append(f"<span>{_esc(unit)}{scale}</span></div>")
    note = ('<div class="sub">columns limited to the %d most prevalent</div>' % max_cols
            if truncated else "")
    return f'<h3>{_esc(title)}</h3>' + "".join(out) + "".join(legend) + note


def render(results: Sequence[SampleResult], *, cell: str = "binary",
           title: str = "Hydra report", databases: Sequence[str] = (),
           command: str = "", extra: dict | None = None,
           element_types: Sequence[str] | None = None,
           matrix_rows: str = "sample", matrix_columns: str = "gene") -> str:
    """Build the complete HTML document."""
    summary = summary_table(results)
    long = long_table(results)
    hetero = heteroresistance_table(results)
    n_samples = len(results)
    # Count what the tables count: one locus reported by three databases is one
    # element, and the heatmap below has one column for it.
    primary = [h for r in results for h in r.hits if h.primary]
    n_hits = len(primary)
    amr_genes = len({h.gene for h in primary if h.element_type == "AMR"})
    vir_genes = len({h.gene for h in primary if h.element_type == "VIRULENCE"})
    n_point = len([h for h in primary if h.resolution == "POINT"])
    n_hetero = len(hetero[hetero["status"] == "heteroresistant"]) if not hetero.empty else 0
    n_species = len({r.species.name for r in results if r.species.name != "unknown"})
    n_st = len({f"{r.mlst.scheme}:{r.mlst.sequence_type}" for r in results
                if r.mlst.sequence_type not in ("-", "")})

    wanted = {t.upper() for t in element_types} if element_types else None
    show_amr = wanted is None or "AMR" in wanted
    show_virulence = wanted is None or "VIRULENCE" in wanted
    amr_matrix = matrix(results, rows=matrix_rows, columns=matrix_columns, cell=cell,
                        element_types=["AMR"]) if show_amr else pd.DataFrame()
    vir_matrix = matrix(results, rows=matrix_rows, columns=matrix_columns, cell=cell,
                        element_types=["VIRULENCE"]) if show_virulence else pd.DataFrame()
    other_types = sorted(wanted - {"AMR", "VIRULENCE"}) if wanted else []
    classes = class_summary(results) if show_amr else pd.DataFrame()

    cards = [
        (n_samples, "samples"), (n_hits, "elements detected"), (amr_genes, "distinct AMR genes"),
        (vir_genes, "virulence genes"), (n_point, "point mutations"),
        (n_hetero, "heteroresistant sites"), (n_species, "species"), (n_st, "sequence types"),
    ]
    card_html = "".join(
        f'<div class="card"><div class="n">{value}</div><div class="l">{_esc(label)}</div></div>'
        for value, label in cards
    )

    parts: list[str] = []
    parts.append("<h2>Sample overview</h2>")
    parts.append('<div class="controls"><input type="search" data-filter="tbl-summary" '
                 'placeholder="filter samples..."></div>')
    parts.append(_frame_to_html(summary, "tbl-summary"))

    axis = f"{matrix_columns} by {matrix_rows}"
    if show_amr:
        parts.append("<h2>Resistance heatmap</h2>")
        parts.append(_heatmap(amr_matrix, cell, f"AMR {axis} ({cell})"))
        parts.append(_heatmap(classes, "genes", "Drug classes by sample (distinct genes)"))
    if show_virulence and not vir_matrix.empty and vir_matrix.shape[1]:
        parts.append("<h2>Virulence heatmap</h2>")
        parts.append(_heatmap(vir_matrix, cell, f"Virulence {axis} ({cell})"))
    for element_type in other_types:
        extra_matrix = matrix(results, rows=matrix_rows, columns=matrix_columns, cell=cell,
                              element_types=[element_type])
        if not extra_matrix.empty and extra_matrix.shape[1]:
            parts.append(f"<h2>{_esc(element_type.title())} heatmap</h2>")
            parts.append(_heatmap(extra_matrix, cell, f"{element_type} {axis} ({cell})"))

    if not hetero.empty:
        parts.append("<h2>Heteroresistance and point mutations from reads</h2>")
        parts.append('<div class="sub">Allele fractions measured directly from the reads. '
                     'A fraction well below 1.0 indicates the mutation is present in only some '
                     'copies of a multi-copy locus, which an assembly consensus would miss.</div>')
        styled = hetero.copy()
        styled["status"] = styled["status"].map(
            lambda s: {"heteroresistant": "HETERORESISTANT", "fixed": "fixed"}.get(s, s))
        parts.append(_frame_to_html(styled, "tbl-hetero"))

    mlst = mlst_table(results)
    if not mlst.empty:
        parts.append("<h2>MLST</h2>")
        parts.append(_frame_to_html(mlst, "tbl-mlst"))

    typing = typing_table(results)
    if not typing.empty and typing.shape[1] > 2:
        parts.append("<h2>Lineage typing and scores</h2>")
        parts.append(_frame_to_html(typing, "tbl-typing"))

    parts.append("<h2>All detected elements</h2>")
    parts.append('<div class="controls"><input type="search" data-filter="tbl-long" '
                 'placeholder="filter genes, classes, samples..."></div>')
    parts.append(_frame_to_html(long, "tbl-long", max_rows=3000))

    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    db_line = ", ".join(databases) if databases else "-"
    meta_json = json.dumps(extra or {}, indent=2, default=str)
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{_esc(title)}</title><style>{_CSS}</style></head>
<body><div class="wrap">
<header><h1>{_esc(title)}</h1>
<span class="sub">Hydra v{__version__} &middot; {generated}</span></header>
<div class="sub">Databases: {_esc(db_line)}</div>
<div class="cards">{card_html}</div>
{''.join(parts)}
<details><summary>Run parameters</summary><pre class="sub">{_esc(command)}

{_esc(meta_json)}</pre></details>
<footer>Generated by Hydra v{__version__}. Cell values show
<strong>{_esc(cell)}</strong>. Absent elements are blank.</footer>
</div><script>{_JS}</script></body></html>"""


def write_html(path: Path, results: Sequence[SampleResult], **kwargs) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render(results, **kwargs), encoding="utf-8")
    LOG.info("wrote HTML report: %s", path)
