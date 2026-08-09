#!/usr/bin/env python3
"""Write the README performance block from metrics.json, or check it still matches.

The figures and the table said different things about the same measurement -- one
caption read 120,380 gene calls where the table read 116,381 -- because both were
typed by hand from the same run and only one was updated when the run was redone.
Generating both from one file removes the possibility rather than the symptom.

  sync_readme.py METRICS_JSON README.md          rewrite the block
  sync_readme.py METRICS_JSON README.md --check  exit 1 if it would change

--check is what CI runs, so a hand-edited number in the README fails the build
instead of reaching a reader.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

BEGIN = "<!-- METRICS:BEGIN -->"
END = "<!-- METRICS:END -->"


def render(metrics: dict) -> str:
    cohort = metrics["cohort"]
    lines = [BEGIN, "", "## Performance", ""]
    lines.append(cohort["blurb"])
    lines.append("")
    lines.append("| | Hydra | Reference tool |")
    lines.append("|---|---|---|")
    for row in metrics["readme_performance"]:
        lines.append(f'| {row["metric"]} | **{row["hydra"]}** | {row["rival"]} |')
    lines.append("")
    for figure in metrics["readme_figures"]:
        lines.append(f'<p align="center"><img src="{figure["src"]}" '
                     f'alt="{figure["alt"]}" width="820"></p>')
        lines.append("")
    lines.append(END)
    return "\n".join(lines)


def splice(readme: str, block: str) -> str:
    start = readme.index(BEGIN)
    stop = readme.index(END) + len(END)
    return readme[:start] + block + readme[stop:]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("metrics", type=Path)
    parser.add_argument("readme", type=Path)
    parser.add_argument("--check", action="store_true",
                        help="report whether the block is current, change nothing")
    args = parser.parse_args()

    metrics = json.loads(args.metrics.read_text())
    readme = args.readme.read_text()
    if BEGIN not in readme or END not in readme:
        print(f"error: {args.readme} has no {BEGIN} / {END} markers", file=sys.stderr)
        return 2

    updated = splice(readme, render(metrics))
    if args.check:
        if updated == readme:
            print("README performance block matches metrics.json")
            return 0
        print("README performance block is stale; run sync_readme.py without --check",
              file=sys.stderr)
        return 1
    if updated != readme:
        args.readme.write_text(updated)
        print(f"updated {args.readme}")
    else:
        print(f"{args.readme} already current")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
