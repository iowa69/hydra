#!/usr/bin/env python3
"""Build a synthetic heteroresistance positive control.

Reads are simulated from a species 23S rRNA reference, with a defined fraction
carrying a catalogued linezolid-resistance mutation. Because rRNA operons are
multi-copy, a strain in which only some operons are mutated produces exactly
this signal, and the assembly consensus would show the wild-type base.

Usage:
    make_heteroresistance_control.py REFERENCE.fna OUTDIR --position 2577 \\
        --ref-base G --alt-base T --fraction 0.2 --depth 400
"""

from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path


def read_first_record(path: Path) -> tuple[str, str]:
    name = None
    chunks: list[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:]
                continue
            if name is not None:
                chunks.append(line)
    if name is None:
        raise SystemExit(f"no FASTA record in {path}")
    return name, "".join(chunks)


COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference", type=Path)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--position", type=int, required=True,
                        help="1-based position of the mutation in the reference")
    parser.add_argument("--ref-base", required=True)
    parser.add_argument("--alt-base", required=True)
    parser.add_argument("--fraction", type=float, default=0.2,
                        help="fraction of reads carrying the mutation")
    parser.add_argument("--depth", type=int, default=400)
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--insert", type=int, default=350)
    parser.add_argument("--error-rate", type=float, default=0.001)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--sample", default="hetero_control")
    args = parser.parse_args()

    random.seed(args.seed)
    _name, sequence = read_first_record(args.reference)
    sequence = sequence.upper()
    index = args.position - 1
    if not 0 <= index < len(sequence):
        raise SystemExit(f"position {args.position} is outside the {len(sequence)} bp reference")
    observed = sequence[index]
    if observed != args.ref_base.upper():
        raise SystemExit(f"reference base at {args.position} is {observed}, "
                         f"not {args.ref_base.upper()}")
    mutant = sequence[:index] + args.alt_base.upper() + sequence[index + 1:]

    args.outdir.mkdir(parents=True, exist_ok=True)
    r1_path = args.outdir / f"{args.sample}_R1.fastq.gz"
    r2_path = args.outdir / f"{args.sample}_R2.fastq.gz"

    n_pairs = max(1, int(len(sequence) * args.depth / (2 * args.read_length)))
    quality = "I" * args.read_length
    n_mutant = 0
    with gzip.open(r1_path, "wt") as r1, gzip.open(r2_path, "wt") as r2:
        for i in range(n_pairs):
            is_mutant = random.random() < args.fraction
            template = mutant if is_mutant else sequence
            start = random.randint(0, max(0, len(template) - args.insert))
            fragment = template[start:start + args.insert]
            if len(fragment) < args.read_length:
                continue
            forward = fragment[:args.read_length]
            reverse = revcomp(fragment[-args.read_length:])
            forward = add_errors(forward, args.error_rate)
            reverse = add_errors(reverse, args.error_rate)
            if is_mutant and start <= index < start + len(fragment):
                n_mutant += 1
            r1.write(f"@{args.sample}:{i}/1\n{forward}\n+\n{quality}\n")
            r2.write(f"@{args.sample}:{i}/2\n{reverse}\n+\n{quality}\n")

    print(f"wrote {r1_path} and {r2_path}")
    print(f"{n_pairs} pairs, {n_mutant} fragments spanning the site carry "
          f"{args.ref_base}{args.position}{args.alt_base}")
    print(f"expected allele fraction at the site: ~{args.fraction:.2f}")
    return 0


def add_errors(seq: str, rate: float) -> str:
    if rate <= 0:
        return seq
    out = []
    for base in seq:
        if random.random() < rate:
            out.append(random.choice([b for b in "ACGT" if b != base]))
        else:
            out.append(base)
    return "".join(out)


if __name__ == "__main__":
    raise SystemExit(main())
