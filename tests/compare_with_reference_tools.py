#!/usr/bin/env python3
"""Concordance report: Hydra versus abricate and mlst on the same genomes.

Both reference tools are run separately; this script only compares their output
files with Hydra's, so it can be re-run without repeating any analysis.

    compare_with_reference_tools.py --hydra-dir RESULTS --mlst ref/mlst.tsv \\
        --abricate ref/abricate_ncbi.tsv
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def sample_of(path: str) -> str:
    name = Path(path).name
    for ext in (".fasta", ".fna", ".fa"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def load_hydra_mlst(path: Path) -> dict[str, tuple[str, str]]:
    out = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            out[row["sample"]] = (row["scheme"], row["ST"])
    return out


def load_reference_mlst(path: Path) -> dict[str, tuple[str, str]]:
    out = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            out[sample_of(fields[0])] = (fields[1], fields[2])
    return out


def load_hydra_genes(path: Path, database: str) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["database"] != database:
                continue
            out[row["sample"]].add(row["gene"])
    return out


def load_abricate(path: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[0].startswith("#"):
                continue
            out[sample_of(fields[0])].add(fields[5])
    return out


def load_hydra_loci(path: Path, database: str) -> dict[str, list[tuple]]:
    """sample -> [(contig, start, end, gene, identity)] for one database."""
    out: dict[str, list[tuple]] = defaultdict(list)
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["database"] != database:
                continue
            out[row["sample"]].append((row["sequence"], int(row["start"]), int(row["end"]),
                                       row["gene"], float(row["identity_pct"])))
    return out


def load_abricate_loci(path: Path) -> dict[str, list[tuple]]:
    out: dict[str, list[tuple]] = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11 or fields[0].startswith("#"):
                continue
            try:
                start, end = int(fields[2]), int(fields[3])
                identity = float(fields[10])
            except ValueError:
                continue
            out[sample_of(fields[0])].append((fields[1], start, end, fields[5], identity))
    return out


def _overlaps(a: tuple, b: tuple) -> bool:
    if a[0] != b[0]:
        return False
    lo = max(min(a[1], a[2]), min(b[1], b[2]))
    hi = min(max(a[1], a[2]), max(b[1], b[2]))
    shared = hi - lo + 1
    shortest = min(abs(a[2] - a[1]), abs(b[2] - b[1])) + 1
    return shared >= 0.5 * shortest


def compare_loci(hydra: dict[str, list[tuple]], reference: dict[str, list[tuple]],
                 identity_floor: float = 0.0) -> dict:
    """Match hits by genomic position rather than by gene name."""
    both = hydra_only = ref_only = 0
    name_same = name_diff = 0
    diffs: list[str] = []
    for sample in sorted(set(hydra) | set(reference)):
        h = [x for x in hydra.get(sample, []) if x[4] >= identity_floor]
        r = [x for x in reference.get(sample, []) if x[4] >= identity_floor]
        matched_r: set[int] = set()
        for hit in h:
            partner = None
            for index, other in enumerate(r):
                if index in matched_r:
                    continue
                if _overlaps(hit, other):
                    partner = index
                    break
            if partner is None:
                hydra_only += 1
                continue
            matched_r.add(partner)
            both += 1
            if hit[3] == r[partner][3]:
                name_same += 1
            else:
                name_diff += 1
                if len(diffs) < 10:
                    diffs.append(f"    {sample} {hit[0]}:{hit[1]}  "
                                 f"hydra {hit[3]} ({hit[4]:.1f}%) vs "
                                 f"abricate {r[partner][3]} ({r[partner][4]:.1f}%)")
        ref_only += len(r) - len(matched_r)
    return {"both": both, "hydra_only": hydra_only, "ref_only": ref_only,
            "name_same": name_same, "name_diff": name_diff, "diffs": diffs}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hydra-dir", type=Path, required=True)
    parser.add_argument("--prefix", default="hydra")
    parser.add_argument("--mlst", type=Path)
    parser.add_argument("--abricate", type=Path)
    parser.add_argument("--database", default="ncbi",
                        help="Hydra database to compare with abricate (default: ncbi)")
    parser.add_argument("--show", type=int, default=8, help="example discordances to print")
    args = parser.parse_args()

    if args.mlst and args.mlst.exists():
        hydra = load_hydra_mlst(args.hydra_dir / f"{args.prefix}.mlst.tsv")
        reference = load_reference_mlst(args.mlst)
        shared = sorted(set(hydra) & set(reference))
        scheme_same = st_same = st_both_called = st_hydra_only = st_ref_only = 0
        examples = []
        for sample in shared:
            h_scheme, h_st = hydra[sample]
            r_scheme, r_st = reference[sample]
            if h_scheme == r_scheme:
                scheme_same += 1
            if h_st == r_st:
                st_same += 1
            elif h_st not in ("-", "") and r_st in ("-", ""):
                st_hydra_only += 1
                if len(examples) < args.show:
                    examples.append(f"    {sample}: hydra {h_scheme}/{h_st}, mlst {r_scheme}/{r_st}")
            elif h_st in ("-", "") and r_st not in ("-", ""):
                st_ref_only += 1
                if len(examples) < args.show:
                    examples.append(f"    {sample}: hydra {h_scheme}/{h_st}, mlst {r_scheme}/{r_st}")
            else:
                st_both_called += 1
                if len(examples) < args.show:
                    examples.append(f"    {sample}: hydra {h_scheme}/{h_st}, mlst {r_scheme}/{r_st}")
        print(f"MLST ({len(shared)} genomes compared with the reference `mlst` tool)")
        print(f"  same scheme        {scheme_same}/{len(shared)} "
              f"({100*scheme_same/max(1,len(shared)):.1f}%)")
        print(f"  same ST            {st_same}/{len(shared)} "
              f"({100*st_same/max(1,len(shared)):.1f}%)")
        print(f"  ST only by hydra   {st_hydra_only}")
        print(f"  ST only by mlst    {st_ref_only}")
        print(f"  different ST       {st_both_called}")
        for line in examples:
            print(line)
        print()

    if args.abricate and args.abricate.exists():
        hydra_loci = load_hydra_loci(args.hydra_dir / f"{args.prefix}.tsv", args.database)
        reference_loci = load_abricate_loci(args.abricate)
        for floor, label in ((0.0, "all hits"), (98.0, "hits >=98% identity")):
            stats = compare_loci(hydra_loci, reference_loci, floor)
            union = stats["both"] + stats["hydra_only"] + stats["ref_only"]
            print(f"Loci detected, database '{args.database}', {label}")
            print(f"  found by both      {stats['both']}/{union} "
                  f"({100*stats['both']/max(1,union):.1f}%)")
            print(f"  hydra only         {stats['hydra_only']}")
            print(f"  abricate only      {stats['ref_only']}")
            named = stats["name_same"] + stats["name_diff"]
            print(f"  same allele name   {stats['name_same']}/{named} "
                  f"({100*stats['name_same']/max(1,named):.1f}% of shared loci)")
            for line in stats["diffs"][:5]:
                print(line)
            print()

        hydra = load_hydra_genes(args.hydra_dir / f"{args.prefix}.tsv", args.database)
        reference = load_abricate(args.abricate)
        samples = sorted(set(hydra) | set(reference))
        shared_total = hydra_only_total = ref_only_total = 0
        per_sample = []
        for sample in samples:
            h = hydra.get(sample, set())
            r = reference.get(sample, set())
            shared_total += len(h & r)
            hydra_only_total += len(h - r)
            ref_only_total += len(r - h)
            if h - r or r - h:
                per_sample.append((sample, sorted(h - r), sorted(r - h)))
        union = shared_total + hydra_only_total + ref_only_total
        print(f"Acquired genes, database '{args.database}' "
              f"({len(samples)} genomes compared with abricate)")
        print(f"  called by both     {shared_total} "
              f"({100*shared_total/max(1,union):.1f}% of the union)")
        print(f"  hydra only         {hydra_only_total}")
        print(f"  abricate only      {ref_only_total}")
        for sample, extra, missing in per_sample[:args.show]:
            print(f"    {sample}")
            if extra:
                print(f"      hydra only:    {', '.join(extra[:12])}")
            if missing:
                print(f"      abricate only: {', '.join(missing[:12])}")
        if len(per_sample) > args.show:
            print(f"    ... and {len(per_sample) - args.show} more genomes with differences")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
