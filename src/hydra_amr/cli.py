"""Hydra's command line."""

from __future__ import annotations

import argparse
import importlib.util
import os
import shlex
import sys
from dataclasses import replace
from pathlib import Path

from . import __version__, presets
from .config import (CELL_MODES, ELEMENT_TYPES, FORMAT_ALIASES, MATRIX_FIELDS,
                     OUTPUT_FORMATS, Config, Thresholds, default_db_dir)
from .db.fetch import can_fetch
from .db.manager import DatabaseStore, create_bundle, download_bundle, install_bundle
from .db.registry import (DATABASES, DB_GROUPS, DEFAULT_DOWNLOADS, SLOW_DOWNLOADS,
                          protein_dir, resolve_names, spec_for)
from .engines.reads import ReadSet, pair_reads
from .pipeline import Pipeline, RunOptions
from .report.writer import write_outputs
from .seqio import looks_like_fasta, looks_like_fastq
from .utils import (LOG, HydraError, cpu_count, sample_name_from_path, setup_logging,
                    tempdir)

EPILOG = f"""
examples:
  # one isolate, everything on
  hydra run -a isolate.fasta -o results/

  # a whole directory of assemblies, presence/absence matrix + HTML heatmaps
  hydra run assemblies/ -o results/ --preset surveillance

  # paired reads: linezolid heteroresistance from 23S allele fractions
  hydra run -1 s_R1.fq.gz -2 s_R2.fq.gz --organism Staphylococcus_aureus \\
      --preset linezolid -o results/

  # single-database screen straight to stdout as a long table
  hydra screen -d card assemblies/*.fasta --stdout

  # a flat one-row-per-gene table, written to a file
  hydra screen -d card assemblies/*.fasta -o results/ -f genes

  # install the reference databases from local conda environments
  hydra db import

presets:
{presets.describe()}

Hydra v{__version__}
"""


# --------------------------------------------------------------------- helpers
def _prepare_tmpdir(path: Path | None) -> Path | None:
    """Create --tmpdir if it is missing, rather than failing deep inside tempfile."""
    if path is None:
        return None
    try:
        path.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise HydraError(f"cannot use --tmpdir {path}: {exc}") from exc
    return path


def _bool_pair(group, name: str, dest: str, help_on: str, help_off: str) -> None:
    """Add matching ``--x`` / ``--no-x`` flags that default to unset."""
    group.add_argument(f"--{name}", dest=dest, action="store_true", default=None, help=help_on)
    group.add_argument(f"--no-{name}", dest=dest, action="store_false", default=None,
                       help=help_off)


def _add_common(parser: argparse.ArgumentParser) -> None:
    runtime = parser.add_argument_group("runtime")
    runtime.add_argument("-t", "--threads", type=int, default=None,
                         help=f"CPU threads to use (default: all {cpu_count()})")
    runtime.add_argument("-j", "--jobs", type=int, default=None,
                         help="read-mapping samples to process concurrently "
                              "(default: auto from --threads; assembly screening is "
                              "batched and always uses all threads)")
    runtime.add_argument("--db-dir", type=Path, default=None,
                         help=f"database directory (default: $HYDRA_DB or {default_db_dir()})")
    runtime.add_argument("--tmpdir", type=Path, default=None,
                         help="directory for intermediate files (default: system temp)")
    runtime.add_argument("--keep-temp", action="store_true",
                         help="keep intermediate BLAST/BAM files for debugging")
    runtime.add_argument("-v", "--verbose", action="count", default=0,
                         help="verbose logging (repeat for more)")
    runtime.add_argument("-q", "--quiet", action="store_true", help="only warnings and errors")


def _add_inputs(parser: argparse.ArgumentParser) -> None:
    inputs = parser.add_argument_group("inputs")
    inputs.add_argument("inputs", nargs="*", type=Path, metavar="INPUT",
                        help="assemblies, FASTQ files or directories to scan")
    inputs.add_argument("-a", "--assembly", dest="assemblies", action="append", type=Path,
                        default=[], metavar="FASTA", help="assembly FASTA (repeatable)")
    inputs.add_argument("-1", "--r1", dest="r1", action="append", type=Path, default=[],
                        metavar="FASTQ", help="forward reads (repeatable, pairs with --r2)")
    inputs.add_argument("-2", "--r2", dest="r2", action="append", type=Path, default=[],
                        metavar="FASTQ", help="reverse reads (repeatable)")
    inputs.add_argument("--reads", dest="reads", action="append", type=Path, default=[],
                        metavar="FASTQ", help="reads to auto-pair by filename (repeatable)")
    inputs.add_argument("--input-list", type=Path, default=None, metavar="TSV",
                        help="tab-separated sample sheet: sample, assembly, R1, R2")
    inputs.add_argument("--name", dest="names", action="append", default=[], metavar="NAME",
                        help="override sample names, in input order (repeatable)")


def _add_databases(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("databases")
    group.add_argument("-d", "--db", dest="db", action="append", default=None, metavar="NAME",
                       help="database or group to screen against; repeatable or comma-separated. "
                            f"Groups: {', '.join(sorted(DB_GROUPS))}")
    group.add_argument("--list-databases", action="store_true",
                       help="print the installed databases and exit")


def _add_analysis(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("analysis")
    group.add_argument("--preset", default=None, metavar="NAME",
                       help=f"option bundle: {', '.join(sorted(presets.PRESETS))}")
    group.add_argument("-O", "--organism", default=None, metavar="NAME",
                       help="organism for point mutations, e.g. Staphylococcus_aureus "
                            "(default: detected from MLST/Mash)")
    group.add_argument("--list-organisms", action="store_true",
                       help="print the organisms with point-mutation support and exit")
    _bool_pair(group, "auto-organism", "auto_organism",
               "detect the organism automatically (default)",
               "never infer the organism; only use --organism")
    _bool_pair(group, "plus", "plus",
               "also report stress-response and virulence elements from the "
               "protein reference",
               "report acquired resistance only, even if the preset asked for --plus")
    _bool_pair(group, "mlst", "mlst", "run MLST (default)", "skip MLST")
    group.add_argument("--scheme", default=None, metavar="NAME",
                       help="force a PubMLST scheme instead of choosing automatically")
    _bool_pair(group, "typing", "typing", "run lineage typing (default)", "skip lineage typing")
    _bool_pair(group, "protein", "protein",
               "run the translated AMR search (default)", "skip the translated search")
    _bool_pair(group, "point-mutations", "point_mutations",
               "call resistance point mutations (default)", "skip point mutations")
    _bool_pair(group, "heteroresistance", "heteroresistance",
               "measure allele fractions from reads (default when reads are given)",
               "skip heteroresistance analysis")
    _bool_pair(group, "reads-mlst", "reads_mlst",
               "type from reads by mapping them to the scheme's loci when no assembly "
               "produced an ST (default)",
               "never type from reads")
    _bool_pair(group, "reads-variants", "reads_variants",
               "report differences between the reads and the closest reference of each "
               "resistance gene (default when reads are given)",
               "skip read-derived gene variants")
    group.add_argument("--report-synonymous", action="store_true", default=None,
                       help="also report read-derived variants that do not change the protein")
    group.add_argument("--assemble", action="store_true", default=None,
                       help="assemble read-only samples first (needs skesa or spades.py)")


def _add_thresholds(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("thresholds")
    group.add_argument("--min-identity", type=float, default=None, metavar="PCT",
                       help="minimum %% identity for nucleotide hits (default: 80)")
    group.add_argument("--min-coverage", type=float, default=None, metavar="PCT",
                       help="minimum %% of the reference covered (default: 60)")
    group.add_argument("--protein-min-identity", type=float, default=None, metavar="PCT",
                       help="minimum %% identity for translated hits (default: 90)")
    group.add_argument("--protein-min-coverage", type=float, default=None, metavar="PCT",
                       help="coverage above which a translated hit is complete (default: 90)")
    group.add_argument("--min-contig-length", type=int, default=None, metavar="BP",
                       help="ignore contigs shorter than this (default: 0)")
    group.add_argument("--report-overlaps", action="store_true", default=None,
                       help="report every overlapping hit instead of the best one per locus")
    reads = parser.add_argument_group("read and heteroresistance thresholds")
    reads.add_argument("--min-depth", type=int, default=None, metavar="N",
                       help="minimum read depth for a call (default: 5)")
    reads.add_argument("--min-gene-breadth", type=float, default=None, metavar="PCT",
                       help="minimum %% of a gene covered by reads (default: 80)")
    reads.add_argument("--min-allele-fraction", type=float, default=None, metavar="FRAC",
                       help="lowest minority allele fraction to report (default: 0.02)")
    reads.add_argument("--fixed-allele-fraction", type=float, default=None, metavar="FRAC",
                       help="fraction at or above which a mutation is called fixed "
                            "rather than heteroresistant (default: 0.90)")
    reads.add_argument("--min-allele-reads", type=int, default=None, metavar="N",
                       help="minimum reads supporting a minority allele (default: 3)")
    reads.add_argument("--min-base-quality", type=int, default=None, metavar="Q",
                       help="minimum base quality in the pileup (default: 13)")
    reads.add_argument("--report-absent-sites", action="store_true", default=None,
                       help="also report catalogued sites where no mutation was found")


def _add_outputs(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("output")
    group.add_argument("-o", "--outdir", type=Path, default=None, metavar="DIR",
                       help="write results here (created if missing)")
    group.add_argument("--prefix", default="hydra", metavar="NAME",
                       help="basename for output files (default: hydra)")
    group.add_argument("-f", "--format", dest="format", action="append", default=None,
                       metavar="FMT",
                       help=f"output formats, repeatable or comma-separated: "
                            f"{', '.join(OUTPUT_FORMATS)} (default: tsv,html)")
    group.add_argument("--stdout", action="store_true",
                       help="write the long-format table to stdout instead of files")
    group.add_argument("--cell", default=None, choices=list(CELL_MODES), metavar="MODE",
                       help=f"matrix cell values: {', '.join(CELL_MODES)} (default: binary)")
    group.add_argument("--rows", default="sample", metavar="FIELD", choices=MATRIX_FIELDS,
                       help=f"matrix row field: {', '.join(MATRIX_FIELDS)} (default: sample)")
    group.add_argument("--columns", default="gene", metavar="FIELD", choices=MATRIX_FIELDS,
                       help=f"matrix column field: {', '.join(MATRIX_FIELDS)} (default: gene)")
    group.add_argument("--element-types", default=None, metavar="LIST",
                       help=f"restrict the matrix and heatmaps to these element types: "
                            f"{', '.join(ELEMENT_TYPES)} (comma-separated)")
    group.add_argument("--title", default=None, metavar="TEXT",
                       help="title for the HTML report")


# ------------------------------------------------------------------- collection
def _scan_directory(path: Path) -> tuple[list[Path], list[Path]]:
    fasta: list[Path] = []
    fastq: list[Path] = []
    for entry in sorted(path.rglob("*")):
        if not entry.is_file():
            continue
        if looks_like_fastq(entry):
            fastq.append(entry)
        elif looks_like_fasta(entry):
            fasta.append(entry)
    return fasta, fastq


def collect_inputs(args) -> tuple[dict[str, Path], dict[str, ReadSet]]:
    """Turn the input flags into {sample: assembly} and {sample: ReadSet}."""
    assemblies: dict[str, Path] = {}
    read_paths: list[Path] = []
    readsets: dict[str, ReadSet] = {}

    def add_assembly(path: Path, name: str | None = None) -> None:
        if not path.exists():
            raise HydraError(f"input not found: {path}")
        sample = name or sample_name_from_path(path)
        if sample in assemblies and assemblies[sample] != path:
            LOG.warning("two assemblies map to sample name '%s'; keeping %s",
                        sample, assemblies[sample])
            return
        assemblies[sample] = path

    for path in args.inputs:
        if path.is_dir():
            found_fasta, found_fastq = _scan_directory(path)
            if not found_fasta and not found_fastq:
                LOG.warning("no FASTA or FASTQ files found under %s", path)
            for entry in found_fasta:
                add_assembly(entry)
            read_paths.extend(found_fastq)
        elif path.exists():
            if looks_like_fastq(path):
                read_paths.append(path)
            elif looks_like_fasta(path):
                add_assembly(path)
            else:
                raise HydraError(f"cannot tell whether {path} is FASTA or FASTQ; "
                                 f"pass it with --assembly or --reads")
        else:
            raise HydraError(f"input not found: {path}")

    for path in args.assemblies:
        add_assembly(path)
    read_paths.extend(args.reads)

    if args.r1 or args.r2:
        if args.r2 and len(args.r1) != len(args.r2):
            raise HydraError(f"--r1 was given {len(args.r1)} time(s) but --r2 "
                             f"{len(args.r2)} time(s); they must pair up")
        for index, first in enumerate(args.r1):
            second = args.r2[index] if index < len(args.r2) else None
            for candidate in (first, second):
                if candidate is not None and not candidate.exists():
                    raise HydraError(f"input not found: {candidate}")
            sample = sample_name_from_path(first, strip_read_suffix=True)
            readsets[sample] = ReadSet(sample=sample, r1=first, r2=second,
                                       single=second is None)

    for readset in pair_reads(read_paths):
        for path in readset.files:
            if not path.exists():
                raise HydraError(f"input not found: {path}")
        readsets.setdefault(readset.sample, readset)

    if args.input_list:
        assemblies_extra, readsets_extra = _read_sample_sheet(args.input_list)
        assemblies.update(assemblies_extra)
        readsets.update(readsets_extra)

    if args.names:
        # Rename assemblies and read sets together so a sample with both keeps
        # them paired under one name.
        ordered = list(dict.fromkeys(list(assemblies) + list(readsets)))
        if len(args.names) > len(ordered):
            raise HydraError(f"--name given {len(args.names)} times but only "
                             f"{len(ordered)} sample(s) were collected")
        mapping = dict(zip(ordered, args.names))
        renamed_assemblies = {mapping.get(k, k): v for k, v in assemblies.items()}
        renamed_reads = {}
        for key, readset in readsets.items():
            new = mapping.get(key, key)
            renamed_reads[new] = replace(readset, sample=new)
        # An assembly and a read set sharing a name is the intended way to pair
        # them; two assemblies sharing one is silent data loss.
        for label, before, after in (("assemblies", assemblies, renamed_assemblies),
                                     ("read sets", readsets, renamed_reads)):
            if len(after) < len(before):
                lost = sorted(set(before) - set(after))
                raise HydraError(
                    f"--name would give more than one of the {label} the same name, "
                    f"discarding {', '.join(lost)}. Give every sample a distinct name.")
        assemblies, readsets = renamed_assemblies, renamed_reads

    return assemblies, readsets


def _read_sample_sheet(path: Path) -> tuple[dict[str, Path], dict[str, ReadSet]]:
    if not path.exists():
        raise HydraError(f"--input-list not found: {path}")
    assemblies: dict[str, Path] = {}
    readsets: dict[str, ReadSet] = {}
    with open(path) as handle:
        for lineno, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = [f.strip() for f in line.split("\t")]
            if len(fields) == 1:
                fields = [f.strip() for f in line.split(",")]
            sample = fields[0]
            if not sample:
                raise HydraError(f"{path}:{lineno}: empty sample name")
            assembly = Path(fields[1]) if len(fields) > 1 and fields[1] else None
            r1 = Path(fields[2]) if len(fields) > 2 and fields[2] else None
            r2 = Path(fields[3]) if len(fields) > 3 and fields[3] else None
            if assembly:
                if not assembly.exists():
                    raise HydraError(f"{path}:{lineno}: assembly not found: {assembly}")
                assemblies[sample] = assembly
            if r1:
                if not r1.exists():
                    raise HydraError(f"{path}:{lineno}: reads not found: {r1}")
                if r2 and not r2.exists():
                    raise HydraError(f"{path}:{lineno}: reads not found: {r2}")
                readsets[sample] = ReadSet(sample=sample, r1=r1, r2=r2, single=r2 is None)
            if not assembly and not r1:
                raise HydraError(f"{path}:{lineno}: sample '{sample}' has no assembly or reads")
    return assemblies, readsets


# ---------------------------------------------------------------- option merge
def _apply_preset(args) -> dict:
    """Fill unset options from the chosen preset. Returns the applied values."""
    name = args.preset or ("standard" if getattr(args, "_default_preset", True) else None)
    if not name:
        return {}
    preset = presets.get(name)
    applied = {}
    for key, value in preset.options.items():
        if getattr(args, key, None) is None:
            setattr(args, key, value)
            applied[key] = value
    return applied


def _build_config(args, store_dir: Path) -> Config:
    thresholds = Thresholds()
    for attribute in ("min_identity", "min_coverage", "protein_min_identity",
                      "protein_min_coverage", "min_depth", "min_gene_breadth",
                      "min_allele_fraction", "fixed_allele_fraction", "min_allele_reads"):
        value = getattr(args, attribute, None)
        if value is not None:
            thresholds.set_explicit(attribute, value)
    thresholds.validate()

    if args.threads is not None and args.threads < 1:
        raise HydraError(f"--threads must be at least 1 (got {args.threads})")
    if args.jobs is not None and args.jobs < 1:
        raise HydraError(f"--jobs must be at least 1 (got {args.jobs})")
    config = Config(
        db_dir=store_dir,
        threads=args.threads if args.threads else cpu_count(),
        jobs=args.jobs or 0,
        thresholds=thresholds,
        organism=args.organism,
        auto_organism=True if args.auto_organism is None else args.auto_organism,
        plus=bool(args.plus),
        keep_temp=args.keep_temp,
        tmp_dir=_prepare_tmpdir(args.tmpdir),
        report_overlaps=bool(getattr(args, "report_overlaps", False)),
    )
    return config


def _resolve_formats(args) -> list[str]:
    raw = args.format or ["tsv", "html"]
    out: list[str] = []
    for item in raw:
        for token in str(item).split(","):
            token = token.strip().lower()
            if not token:
                continue
            token = FORMAT_ALIASES.get(token, token)
            if token not in OUTPUT_FORMATS:
                raise HydraError(f"unknown --format '{token}'; choose from "
                                 f"{', '.join(OUTPUT_FORMATS)}")
            if token not in out:
                out.append(token)
    return out or ["tsv"]


def _resolve_databases(args, store: DatabaseStore) -> list[str]:
    names = resolve_names(args.db or ["standard"])
    missing = [n for n in names if not store.is_installed(n)]
    available = [n for n in names if store.is_installed(n)]
    if missing:
        for name in missing:
            LOG.warning("database '%s' is not installed; skipping it", name)
    if not available:
        raise HydraError(
            f"none of the requested databases are installed in {store.root}.\n"
            f"Install them with:  hydra db import      (from local conda environments)\n"
            f"                    hydra db download    (from upstream sources)\n"
            f"Requested: {', '.join(names)}")
    return available


def _mutation_counts(store: DatabaseStore) -> dict[str, tuple[int, int]]:
    """How many protein and DNA mutations each organism has catalogued.

    A bare list of names does not say whether asking for an organism will do
    anything; these counts do.
    """
    counts: dict[str, list[int]] = {}
    root = protein_dir(store.root)
    table = root / "AMRProt-mutation.tsv"
    if table.exists():
        with open(table, errors="replace") as handle:
            header = handle.readline().split("\t")
            try:
                column = header.index("taxgroup")
            except ValueError:
                column = 0
            for line in handle:
                fields = line.split("\t")
                if len(fields) > column:
                    counts.setdefault(fields[column].strip(), [0, 0])[0] += 1
    for path in sorted((store.root / "mutation" / "dna").glob("*.tsv")):
        with open(path, errors="replace") as handle:
            rows = sum(1 for line in handle if line.strip()) - 1
        counts.setdefault(path.stem, [0, 0])[1] += max(rows, 0)
    return {name: (value[0], value[1]) for name, value in counts.items()}


def _organism_choices(store: DatabaseStore) -> list[str]:
    organisms = set(store.mutation_organisms())
    taxgroup = protein_dir(store.root) / "taxgroup.tsv"
    if taxgroup.exists():
        with open(taxgroup) as handle:
            handle.readline()
            for line in handle:
                if line.strip():
                    organisms.add(line.split("\t")[0].strip())
    return sorted(organisms)


# ------------------------------------------------------------------- commands
def cmd_run(args) -> int:
    _apply_preset(args)
    store_dir = args.db_dir or default_db_dir()
    store = DatabaseStore(store_dir)

    if args.list_databases:
        return _print_databases(store)
    if args.list_organisms:
        choices = _organism_choices(store)
        if not choices:
            print("no organisms available; install the protein database first "
                  "(hydra db import)")
            return 1
        counts = _mutation_counts(store)
        print("Catalogued point mutations per organism. "
              "Pass one to -O/--organism.\n")
        print(f"{'ORGANISM':34s} PROTEIN   DNA")
        for name in choices:
            protein, dna = counts.get(name, (0, 0))
            print(f"{name:34s} {protein:7d}  {dna:4d}")
        print(f"\n{len(choices)} organisms. Without -O, Hydra picks the catalogue "
              f"from its own species call;\nan organism with no catalogued mutations "
              f"is still screened for acquired genes.")
        return 0

    assemblies, readsets = collect_inputs(args)
    if not assemblies and not readsets:
        raise HydraError("no inputs given. Pass assemblies, FASTQ files or a directory; "
                         "see 'hydra run --help'")

    config = _build_config(args, store_dir)
    databases = _resolve_databases(args, store)
    formats = _resolve_formats(args)

    # Argument-only checks come first: they need nothing but the command line, so
    # failing them here reports the mistake before any startup noise.
    element_types = None
    if args.element_types:
        element_types = [t.strip().upper() for t in args.element_types.split(",") if t.strip()]
        unknown = [t for t in element_types if t not in ELEMENT_TYPES]
        if unknown:
            raise HydraError(f"unknown --element-types {', '.join(unknown)}; choose from "
                             f"{', '.join(ELEMENT_TYPES)}")
    if not args.stdout and args.outdir is None:
        raise HydraError("no output directory set; pass --outdir, or --stdout to write the "
                         "long table to standard output")
    if any(sep in args.prefix for sep in ("/", os.sep)):
        raise HydraError(f"--prefix is a file basename, not a path (got '{args.prefix}')")
    if "xlsx" in formats and importlib.util.find_spec("openpyxl") is None:
        # Checked before the analysis runs: discovering it afterwards throws
        # away the whole run.
        raise HydraError("xlsx output needs openpyxl.\n"
                         "  conda install -c conda-forge openpyxl")

    if args.organism:
        choices = _organism_choices(store)
        if choices and args.organism not in choices:
            close = [c for c in choices if c.lower().startswith(args.organism.split("_")[0].lower())]
            hint = f" Did you mean: {', '.join(close[:5])}?" if close else ""
            raise HydraError(f"unknown --organism '{args.organism}'.{hint}\n"
                             f"List them with: hydra run --list-organisms")

    options = RunOptions(
        databases=databases,
        mlst=True if args.mlst is None else args.mlst,
        typing=True if args.typing is None else args.typing,
        protein=True if args.protein is None else args.protein,
        point_mutations=True if args.point_mutations is None else args.point_mutations,
        heteroresistance=True if args.heteroresistance is None else args.heteroresistance,
        reads_mlst=True if args.reads_mlst is None else args.reads_mlst,
        reads_variants=True if args.reads_variants is None else args.reads_variants,
        report_synonymous=bool(args.report_synonymous),
        assemble_reads=bool(args.assemble),
        force_scheme=args.scheme,
        min_contig_length=args.min_contig_length or 0,
        report_absent_sites=bool(args.report_absent_sites),
        min_base_quality=13 if args.min_base_quality is None else args.min_base_quality,
    )
    if args.scheme and options.mlst:
        typer_root = store.root / "mlst" / "pubmlst"
        if typer_root.exists() and not (typer_root / args.scheme).is_dir():
            available = sorted(p.name for p in typer_root.iterdir() if p.is_dir())
            close = [s for s in available if s.startswith(args.scheme[:3])]
            hint = f"\nDid you mean: {', '.join(close[:6])}?" if close else ""
            raise HydraError(
                f"unknown --scheme '{args.scheme}'.{hint}\n"
                f"{len(available)} schemes are installed, for example: "
                f"{', '.join(available[:8])}...\n"
                f"List them all with:  hydra db info pubmlst")
    elif args.scheme and not options.mlst:
        LOG.warning("--scheme has no effect because MLST is switched off")

    LOG.info("Hydra v%s | %d assemblies, %d read sets | databases: %s | %d threads",
             __version__, len(assemblies), len(readsets), ", ".join(databases), config.threads)

    if args.cell in ("depth", "fraction") and not readsets:
        LOG.warning("--cell %s only has values for read-derived hits, and this run has no "
                    "read input; the matrix will be empty", args.cell)

    command = "hydra " + " ".join(shlex.quote(a) for a in sys.argv[1:])
    with tempdir(prefix="hydra.", keep=config.keep_temp, parent=config.tmp_dir) as workdir:
        pipeline = Pipeline(config, store, options)
        results = pipeline.run(assemblies, readsets, workdir)
        if config.keep_temp:
            LOG.info("intermediate files kept in %s", workdir)

        meta = {
            "hydra_version": __version__,
            "databases": {name: store.installed().get(name, {}).get("version", "")
                          for name in databases},
            "thresholds": vars(config.thresholds),
            "preset": args.preset or "standard",
            "organism": args.organism or "auto",
        }
        written = write_outputs(
            results, outdir=args.outdir, prefix=args.prefix, formats=formats,
            cell=args.cell or "binary", matrix_rows=args.rows, matrix_columns=args.columns,
            element_types=element_types, databases=databases, command=command, meta=meta,
            to_stdout=args.stdout, title=args.title or "Hydra report",
        )
    if not args.stdout:
        _print_recap(results, written)
    return 0


def _print_recap(results, written) -> None:
    total_hits = sum(len(r.hits) for r in results)
    hetero = sum(1 for r in results for h in r.hits if "HETERORESISTANT" in h.note.upper())
    print(f"\n{len(results)} sample(s), {total_hits} element(s) detected", file=sys.stderr)
    if hetero:
        print(f"{hetero} heteroresistant site(s) - see the heteroresistance table",
              file=sys.stderr)
    warned = [r for r in results if r.warnings]
    if warned:
        print(f"{len(warned)} sample(s) raised warnings:", file=sys.stderr)
        for result in warned[:10]:
            print(f"  {result.sample}: {result.warnings[0]}", file=sys.stderr)
        if len(warned) > 10:
            print(f"  ... and {len(warned) - 10} more (see the summary table)", file=sys.stderr)
    for path in written:
        print(f"  {path}", file=sys.stderr)


def cmd_screen(args) -> int:  # noqa: D401 - see below
    """Acquired-gene screening only.

    The analysis stages default to off here, but an explicit flag still wins:
    ``hydra screen --mlst`` is a reasonable thing to ask for.
    """
    args.preset = args.preset or "genes"
    args._default_preset = False
    _apply_preset(args)
    if args.db is None:
        args.db = ["ncbi"]
    for name, default in (("mlst", False), ("typing", False),
                          ("point_mutations", False), ("heteroresistance", False),
                          ("reads_mlst", False), ("reads_variants", False)):
        if getattr(args, name) is None:
            setattr(args, name, default)
    if args.protein is None:
        args.protein = "protein" in resolve_names(args.db)
    return cmd_run(args)


def cmd_db(args) -> int:
    store_dir = args.db_dir or default_db_dir()
    store = DatabaseStore(store_dir)
    action = args.db_action

    if action == "list":
        return _print_databases(store)

    if action == "info":
        for name in args.names or sorted(DATABASES):
            spec = spec_for(name)
            entry = store.installed().get(name, {})
            print(f"{name}")
            print(f"  title       {spec.title}")
            print(f"  kind        {spec.kind}")
            print(f"  installed   {'yes, ' + entry.get('installed', '') if entry else 'no'}")
            if entry:
                print(f"  path        {store.root / entry.get('path', '')}")
                print(f"  version     {entry.get('version', 'unknown')}")
                print(f"  sequences   {entry.get('sequences', '-')}")
                print(f"  source      {entry.get('source', '-')}")
            print(f"  upstream    {spec.url or '-'}")
            if spec.citation:
                print(f"  citation    {spec.citation}")
            if spec.licence:
                print(f"  licence     {spec.licence}")
            if spec.notes:
                print(f"  notes       {spec.notes}")
            if entry:
                for label, values in _contents(store, name, entry).items():
                    print(f"  {label:<11} {len(values)} available:")
                    for line in _wrap_columns(values):
                        print(f"              {line}")
            print()
        return 0

    if action == "import":
        only = resolve_names(args.names) if args.names else None
        search = [Path(p) for p in (args.source or [])]
        LOG.info("importing databases into %s", store.root)
        # The protein reference first: its family table annotates every other database.
        ordered = None
        if only is None:
            ordered = ["protein"] + [n for n in DATABASES if n != "protein"]
        elif "protein" in only:
            ordered = ["protein"] + [n for n in only if n != "protein"]
        else:
            ordered = only
        results = store.import_all(only=ordered, search_paths=search, force=args.force)
        width = max(len(n) for n in results) if results else 10
        installed = 0
        for name, status in results.items():
            marker = "ok " if "installed from" in status or "already" in status else "-- "
            if "installed from" in status or "already" in status:
                installed += 1
            print(f"{marker}{name:<{width}}  {status}")
        print(f"\n{installed}/{len(results)} databases available in {store.root}")
        if installed == 0:
            print("\nNothing was found locally. Either install the source tools "
                  "(abricate, ncbi-amrfinderplus, mlst, kleborate) into conda "
                  "environments, or point at a directory with --source.")
            return 1
        return 0

    if action == "download":
        return _download(store, args)

    if action == "update":
        return _update(store, args)

    if action == "bundle":
        return _bundle(store, args)

    if action == "remove":
        if not args.names:
            raise HydraError("give at least one database name to remove")
        import shutil
        for name in resolve_names(args.names):
            entry = store.installed().get(name)
            if not entry:
                print(f"-- {name}: not installed")
                continue
            path = store.root / entry["path"]
            if path.exists():
                shutil.rmtree(path, ignore_errors=True)
            store.manifest()["databases"].pop(name, None)
            store._save_manifest()
            print(f"ok {name}: removed")
        return 0

    if action == "check":
        problems = 0
        stale: list[str] = []
        if not store.installed():
            print(f"no databases installed in {store.root}\n\n"
                  f"Install them with:  hydra db import")
            return 1
        have_amrfinder = store.is_installed("protein")
        for name, entry in sorted(store.installed().items()):
            path = store.root / entry.get("path", "")
            if not path.exists():
                print(f"MISSING {name}: {path}")
                problems += 1
                continue
            try:
                if entry.get("kind") in ("nucl", "prot"):
                    handle = store.handle(name)
                    if not handle.meta:
                        print(f"WARN    {name}: metadata table is empty")
                        problems += 1
                        continue
                    # A nucleotide database imported before the protein reference was
                    # installed carries none of its curated annotation, and
                    # nothing re-annotates it later.
                    if (have_amrfinder and entry.get("kind") == "nucl"
                            and DATABASES.get(name) and DATABASES[name].element_type == "AMR"
                            and not any(row.get("fam_id") for row in handle.meta.values())):
                        stale.append(name)
                print(f"ok      {name}: {entry.get('sequences', '-')} sequences, "
                      f"version {entry.get('version', 'unknown')}")
            except HydraError as exc:
                print(f"BROKEN  {name}: {exc}")
                problems += 1
        if stale:
            print(f"\nSTALE   {', '.join(stale)}: imported before the protein reference, "
                  f"so they carry no curated drug-class or gene-family annotation.\n"
                  f"        Refresh them with:  hydra db import --force {' '.join(stale)}")
        if problems:
            print(f"\n{problems} problem(s) found; re-run 'hydra db import --force'")
            return 1
        if stale:
            return 1
        print(f"\nall databases in {store.root} look healthy")
        return 0

    raise HydraError(f"unknown db action '{action}'")


def _download(store: DatabaseStore, args) -> int:
    """Install a prebuilt database bundle, or explain where the sources live."""
    source = getattr(args, "from_file", None)
    url = getattr(args, "url", None)
    if source or url:
        if source:
            archive = Path(source)
        else:
            cache = store.root.parent / "cache"
            archive = download_bundle(url, cache / Path(url).name)
        installed = install_bundle(store, archive, force=getattr(args, "force", False))
        if not installed:
            print("nothing installed; every database in the bundle is already present "
                  "(use --force to replace them)")
            return 0
        print(f"installed {len(installed)} database(s): {', '.join(sorted(installed))}")
        return _print_databases(store)

    names = resolve_names(args.names) if args.names else list(DEFAULT_DOWNLOADS)
    fetchable = [n for n in names if can_fetch(n)]
    if fetchable and not getattr(args, "list_only", False):
        return _fetch_databases(store, fetchable, names, args)

    print("Hydra does not redistribute third-party sequence data, so there are three\n"
          "ways to get the databases.\n")
    print("1. From tools already installed on this machine (quickest):")
    print("     hydra db import\n")
    print("2. From a prebuilt bundle, for machines with no source tools or no network:")
    print("     hydra db bundle -o hydra-db.tar.gz     # on a machine that has them")
    print("     hydra db download --from-file hydra-db.tar.gz")
    print("     hydra db download --url https://example.org/hydra-db.tar.gz\n")
    print("3. From the upstream sources, unpacked into a directory:")
    print("     hydra db import --source DIR\n")
    for name in names:
        spec = spec_for(name)
        print(f"   {name}\n     {spec.title}\n     {spec.url or 'no public URL recorded'}")
        if spec.licence:
            print(f"     licence: {spec.licence}")
    print("\nThe source tools themselves:")
    print("  conda create -n hydra-db -c bioconda -c conda-forge \\\n"
          "      abricate ncbi-amrfinderplus mlst kleborate")
    print("  amrfinder -u && hydra db import")
    return 0


def _fetch_databases(store: DatabaseStore, fetchable: list[str], asked: list[str],
                     args) -> int:
    """Download and import each database that has an automatic source."""
    force = getattr(args, "force", False)
    manual = [n for n in asked if not can_fetch(n)]
    todo = [n for n in fetchable if force or not store.is_installed(n)]
    if todo:
        print(f"downloading {len(todo)} database(s) into {store.root}: {', '.join(todo)}")
        for name in todo:
            if name in SLOW_DOWNLOADS:
                print(f"  {name} {SLOW_DOWNLOADS[name]}")
        print()
    done, failed, skipped = [], [], []
    for name in fetchable:
        if store.is_installed(name) and not force:
            skipped.append(name)
            continue
        try:
            store.download(name, force=force)
        except HydraError as exc:
            LOG.error("%s: %s", name, exc)
            failed.append(name)
            continue
        done.append(name)

    if skipped:
        print(f"already installed, left alone: {', '.join(skipped)} (--force to replace)")
    if done:
        print(f"installed {len(done)} database(s): {', '.join(done)}")
    if manual:
        print(f"\nno automatic source for: {', '.join(manual)}")
        print("  these are published as landing pages rather than versioned files;")
        print("  'hydra db download --list' says where each one lives, and")
        print("  'hydra db import' converts a copy already on this machine.")
    if failed:
        print(f"\nfailed: {', '.join(failed)}")
        return 1
    if done or skipped:
        print()
        _print_databases(store)
    return 0



def _update(store: DatabaseStore, args) -> int:
    """Refresh installed databases from their upstream sources.

    Never runs on its own. Hydra will not reach the network unless asked, because a
    database that changes underneath a running study changes its results: two
    isolates screened a week apart should be comparable, and silently pulling a new
    CARD release between them means they are not. Updating is a decision, so it is a
    command.

    Only databases that are already installed and have an automatic source are
    touched. Anything imported from a local conda environment is left alone -- its
    source is that environment, and `hydra db import --force` is what refreshes it.
    """
    installed = store.installed()
    if not installed:
        print(f"no databases installed in {store.root}\n\n"
              f"Install them first with:  hydra db download")
        return 1

    wanted = resolve_names(args.names) if args.names else sorted(installed)
    refreshable, skipped = [], []
    for name in wanted:
        if name not in installed:
            skipped.append((name, "not installed"))
        elif not can_fetch(name):
            skipped.append((name, "no automatic source; use 'hydra db import --force'"))
        else:
            refreshable.append(name)

    for name, why in skipped:
        print(f"-- {name}: {why}")
    if not refreshable:
        print("\nnothing to update")
        return 0

    print(f"{'database':16s} {'installed':22s} {'action'}")
    for name in refreshable:
        entry = installed.get(name, {})
        stamp = entry.get("installed") or entry.get("version") or "unknown"
        print(f"{name:16s} {str(stamp)[:22]:22s} "
              f"{'would be refreshed' if args.dry_run else 'refreshing'}")
    if args.dry_run:
        print(f"\n{len(refreshable)} database(s) would be refreshed. "
              f"Run without --dry-run to do it.")
        return 0

    # Replace one at a time. A failure part-way leaves every other database as it
    # was rather than a store half-way between two releases.
    failed = []
    for name in refreshable:
        try:
            store.download(name, force=True)
            print(f"ok {name}: refreshed")
        except HydraError as exc:
            failed.append(name)
            LOG.error("%s: %s", name, exc)
    print(f"\n{len(refreshable) - len(failed)} refreshed, {len(failed)} failed")
    if failed:
        print("the databases that failed are unchanged, not half-written")
    return 1 if failed else 0

def _bundle(store: DatabaseStore, args) -> int:
    names = resolve_names(args.names) if args.names else None
    output = Path(args.output)
    create_bundle(store, output, names=names, compress=args.compress)
    print(f"wrote {output}\n\nInstall it elsewhere with:\n"
          f"  hydra db download --from-file {output.name}")
    return 0


def _contents(store: DatabaseStore, name: str, entry: dict) -> dict[str, list[str]]:
    """The named things an installed database provides, for ``hydra db info``."""
    kind = entry.get("kind", "")
    path = store.root / entry.get("path", "")
    if kind == "mlst":
        pubmlst = path / "pubmlst"
        if pubmlst.is_dir():
            return {"schemes": sorted(p.name for p in pubmlst.iterdir() if p.is_dir())}
    if kind == "typing":
        from .typing.lineage import EXCLUDED_MODULES  # noqa: PLC0415 - avoids a cycle at import

        modules = sorted(entry.get("modules", []))
        used = [m for m in modules if m not in EXCLUDED_MODULES]
        listing = {"schemes": used}
        skipped = [m for m in modules if m in EXCLUDED_MODULES]
        if skipped:
            listing["not typed"] = skipped
        return listing
    if kind == "prot":
        organisms = store.mutation_organisms()
        if organisms:
            return {"organisms": organisms}
    return {}


def _wrap_columns(values: list[str], width: int = 88) -> list[str]:
    lines: list[str] = []
    current = ""
    for value in values:
        candidate = f"{current}, {value}" if current else value
        if len(candidate) > width:
            lines.append(current)
            current = value
        else:
            current = candidate
    if current:
        lines.append(current)
    return lines


def _print_databases(store: DatabaseStore) -> int:
    rows = store.summary_rows()
    if not rows:
        print(f"No databases installed in {store.root}\n\nInstall them with:  hydra db import")
        return 1
    headers = ("NAME", "KIND", "SEQUENCES", "VERSION", "SIZE", "TITLE")
    widths = [max(len(h), max((len(str(r[k])) for r in rows), default=0))
              for h, k in zip(headers, ("name", "kind", "sequences", "version", "size", "title"))]
    print(f"database directory: {store.root}\n")
    print("  ".join(h.ljust(w) for h, w in zip(headers, widths)))
    for row in rows:
        print("  ".join(str(row[k]).ljust(w) for k, w in
                        zip(("name", "kind", "sequences", "version", "size", "title"), widths)))
    not_installed = [n for n in DATABASES if n not in {r["name"] for r in rows}]
    if not_installed:
        print(f"\nnot installed: {', '.join(sorted(not_installed))}")
    return 0


def cmd_presets(args) -> int:
    if args.name:
        preset = presets.get(args.name)
        print(f"{preset.name}\n  {preset.summary}")
        if preset.detail:
            print(f"  {preset.detail}")
        print("\n  options set:")
        for key, value in sorted(preset.options.items()):
            print(f"    {key:<22} {value}")
        return 0
    print("presets:\n")
    print(presets.describe())
    print("\nUse one with:  hydra run --preset NAME ...")
    return 0


# ---------------------------------------------------------------------- parser
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hydra",
        description="Hydra - AMR, virulence, point mutations, heteroresistance, "
                    "MLST and lineage typing from assemblies or raw reads.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"hydra {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    run_parser = subparsers.add_parser(
        "run", help="full analysis of assemblies and/or reads",
        description="Screen assemblies and reads for AMR genes, virulence factors, "
                    "point mutations and heteroresistance, and type them.",
        epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
    _add_inputs(run_parser)
    _add_databases(run_parser)
    _add_analysis(run_parser)
    _add_thresholds(run_parser)
    _add_outputs(run_parser)
    _add_common(run_parser)
    run_parser.set_defaults(func=cmd_run, _default_preset=True)

    screen_parser = subparsers.add_parser(
        "screen", help="acquired-gene screening only, as a flat gene table",
        description="Screen assemblies against nucleotide gene databases. "
                    "Typing and point mutations are off by default here, but the "
                    "flags still work: --mlst or --point-mutations turns one back "
                    "on for a single run.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    _add_inputs(screen_parser)
    _add_databases(screen_parser)
    _add_analysis(screen_parser)
    _add_thresholds(screen_parser)
    _add_outputs(screen_parser)
    _add_common(screen_parser)
    screen_parser.set_defaults(func=cmd_screen, _default_preset=False)

    db_parser = subparsers.add_parser(
        "db", help="install and inspect reference databases",
        description="Hydra keeps its databases in $HYDRA_DB (default ~/.hydra/db).")
    db_sub = db_parser.add_subparsers(dest="db_action", metavar="ACTION")
    for action, help_text in (
        ("list", "show installed databases"),
        ("import", "build databases from local conda environments or a directory"),
        ("download", "install a prebuilt bundle, or show where the sources live"),
        ("update", "refresh installed databases from upstream (never automatic)"),
        ("bundle", "pack the installed databases into a portable archive"),
        ("info", "show details for one or more databases"),
        ("check", "verify installed databases are usable"),
        ("remove", "delete an installed database"),
    ):
        sub = db_sub.add_parser(action, help=help_text)
        sub.add_argument("names", nargs="*", metavar="NAME",
                         help="database names (default: all)")
        if action == "import":
            sub.add_argument("--source", action="append", default=[], metavar="DIR",
                             help="extra directory to search for database sources")
            sub.add_argument("--force", action="store_true",
                             help="re-import databases that are already installed")
        if action == "download":
            sub.add_argument("--from-file", dest="from_file", type=Path, metavar="ARCHIVE",
                             help="install a bundle already on this machine")
            sub.add_argument("--url", metavar="URL",
                             help="download a bundle over HTTP(S) and install it")
            sub.add_argument("--list", dest="list_only", action="store_true",
                             help="only print where each database comes from, "
                                  "downloading nothing")
            sub.add_argument("--force", action="store_true",
                             help="replace databases that are already installed")
        if action == "update":
            sub.add_argument("--dry-run", action="store_true",
                             help="list what would be refreshed and change nothing")
        if action == "bundle":
            sub.add_argument("-o", "--output", type=Path, default=Path("hydra-db.tar.gz"),
                             metavar="ARCHIVE", help="archive to write (default: hydra-db.tar.gz)")
            sub.add_argument("--compress", default="gz", choices=("gz", "bz2", "xz", "none"),
                             help="compression to use (default: gz)")
        sub.add_argument("--db-dir", type=Path, default=None,
                         help=f"database directory (default: {default_db_dir()})")
        sub.add_argument("-v", "--verbose", action="count", default=0)
        sub.add_argument("-q", "--quiet", action="store_true")
        sub.set_defaults(func=cmd_db, db_action=action)
    db_parser.set_defaults(func=cmd_db, db_action=None)

    preset_parser = subparsers.add_parser("presets", help="list the available presets")
    preset_parser.add_argument("name", nargs="?", help="show one preset in detail")
    preset_parser.add_argument("-v", "--verbose", action="count", default=0)
    preset_parser.add_argument("-q", "--quiet", action="store_true")
    preset_parser.set_defaults(func=cmd_presets)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging(getattr(args, "verbose", 0), getattr(args, "quiet", False))

    if not getattr(args, "command", None):
        parser.print_help()
        return 1
    if args.command == "db" and not getattr(args, "db_action", None):
        print("usage: hydra db {list,import,download,bundle,info,check,remove}\n\n"
              "Start with:  hydra db import", file=sys.stderr)
        return 1

    try:
        return args.func(args)
    except HydraError as exc:
        LOG.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        LOG.error("interrupted")
        return 130
    except BrokenPipeError:
        os._exit(0)


if __name__ == "__main__":
    sys.exit(main())
