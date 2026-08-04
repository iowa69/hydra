"""End-to-end orchestration: from inputs to :class:`SampleResult` objects."""

from __future__ import annotations

import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

from .config import Config
from .db.manager import DatabaseStore
from .db.registry import DATABASES
from .engines.mutations import MutationCatalog
from .engines.nucl import NucleotideScreener, QueryBatch, build_query_batch
from .engines.protein import ProteinScreener
from .engines.reads import ReadMapper, ReadSet, assemble
from .records import Hit, SampleResult
from .seqio import (assembly_stats, fastq_stats, read_fasta, validate_fasta,
                    validate_fastq)
from .typing.lineage import LineageTyper, resistance_score, virulence_score
from .typing.mlst import MlstTyper
from .typing.species import SpeciesIdentifier
from .utils import LOG, HydraError, check_workspace, human_time


@dataclass
class RunOptions:
    """Which analyses to run, and how."""

    databases: list[str] = field(default_factory=lambda: ["ncbi"])
    mlst: bool = True
    typing: bool = True
    protein: bool = True
    point_mutations: bool = True
    reads_genes: bool = True
    #: Type from reads by mapping them to the scheme's loci when no assembly did.
    reads_mlst: bool = True
    #: Report differences between the reads and the closest reference of each gene.
    reads_variants: bool = True
    report_synonymous: bool = False
    heteroresistance: bool = True
    assemble_reads: bool = False
    force_scheme: str | None = None
    min_contig_length: int = 0
    report_absent_sites: bool = False
    min_base_quality: int = 13


class Pipeline:
    """Runs the requested analyses over a batch of samples."""

    def __init__(self, config: Config, store: DatabaseStore, options: RunOptions):
        self.config = config
        self.store = store
        self.options = options
        self._species_identifier: SpeciesIdentifier | None = None

    def _identifier(self) -> SpeciesIdentifier:
        """One identifier for the run: it loads the scheme table and sketch list once."""
        if self._species_identifier is None:
            self._species_identifier = SpeciesIdentifier(self.store)
        return self._species_identifier

    # ------------------------------------------------------------------ entry
    def run(self, assemblies: dict[str, Path], readsets: dict[str, ReadSet],
            workdir: Path) -> list[SampleResult]:
        started = time.time()
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)

        assemblies = dict(assemblies)
        if self.options.assemble_reads:
            assemblies.update(self._assemble_missing(assemblies, readsets, workdir))

        samples = sorted(set(assemblies) | set(readsets))
        if not samples:
            raise HydraError("no input samples")
        results: dict[str, SampleResult] = {}
        for sample in samples:
            has_assembly = sample in assemblies
            has_reads = sample in readsets
            kind = ("assembly+reads" if has_assembly and has_reads
                    else "assembly" if has_assembly else "reads")
            inputs = []
            if has_assembly:
                inputs.append(str(assemblies[sample]))
            if has_reads:
                inputs.extend(str(p) for p in readsets[sample].files)
            results[sample] = SampleResult(sample=sample, input_type=kind, inputs=inputs)

        for readset in readsets.values():
            for path in readset.files:
                validate_fastq(path)

        organism_by_sample: dict[str, str | None] = {}
        batch: QueryBatch | None = None
        if assemblies:
            LOG.info("preparing %d assemblies", len(assemblies))
            for sample, path in assemblies.items():
                validate_fasta(path)
            batch = build_query_batch(assemblies, workdir / "query.fna",
                                      min_contig_length=self.options.min_contig_length)
            if batch.n_contigs == 0:
                raise HydraError("no contigs survived filtering; lower --min-contig-length")
            LOG.info("%d contigs, %.1f Mb total", batch.n_contigs, batch.total_bases / 1e6)
            check_workspace(workdir, batch.total_bases)
            self._assembly_qc(assemblies, results)
            organism_by_sample = self._identify(batch, assemblies, results, workdir)
            self._screen_assemblies(batch, results, organism_by_sample, workdir)

        if readsets:
            self._process_reads(readsets, results, organism_by_sample, workdir)

        for result in results.values():
            self._finalise(result)
            result.runtime_seconds = time.time() - started

        LOG.info("analysed %d samples in %s", len(results), human_time(time.time() - started))
        return [results[sample] for sample in samples]

    # ------------------------------------------------------------------ steps
    def _assemble_missing(self, assemblies: dict[str, Path], readsets: dict[str, ReadSet],
                          workdir: Path) -> dict[str, Path]:
        out: dict[str, Path] = {}
        pending = [s for s in readsets if s not in assemblies]
        if not pending:
            return out
        LOG.info("assembling %d read-only samples", len(pending))
        for sample in pending:
            try:
                out[sample] = assemble(readsets[sample], workdir / "assembly",
                                       threads=self.config.threads)
                LOG.info("assembled %s", sample)
            except HydraError as exc:
                LOG.error("assembly of %s failed: %s", sample, exc)
        return out

    def _assembly_qc(self, assemblies: dict[str, Path], results: dict[str, SampleResult]) -> None:
        for sample, path in assemblies.items():
            stats = assembly_stats(path)
            results[sample].qc.update(stats.as_dict())
            if stats.contigs > 1000:
                results[sample].warnings.append(
                    f"fragmented assembly ({stats.contigs} contigs); "
                    f"genes spanning contig breaks may be reported as partial")
            if stats.total_length < 500_000:
                results[sample].warnings.append(
                    f"assembly is only {stats.total_length/1e6:.2f} Mb; results may be incomplete")

    def _identify(self, batch: QueryBatch, assemblies: dict[str, Path],
                  results: dict[str, SampleResult], workdir: Path) -> dict[str, str | None]:
        """MLST, species call and the AMRFinderPlus organism each sample maps to."""
        organism_by_sample: dict[str, str | None] = {}
        identifier = self._identifier()

        # Sketch-based species first: it is independent of MLST, and its genus
        # keeps scheme selection on the right organism when neighbouring genera
        # share housekeeping alleles.
        sketches: dict[str, object] = {}
        if identifier.sketches:
            for sample, path in assemblies.items():
                call = identifier.sketch_species(path, threads=self.config.threads)
                if call is not None:
                    sketches[sample] = call

        if self.options.mlst and self.store.is_installed("pubmlst"):
            typer = MlstTyper(self.store, self.config)
            scheme_genus = {name: entry.genus for name, entry in identifier.scheme_table.items()
                            if entry.genus}
            genus_by_sample = {sample: call.genus for sample, call in sketches.items()
                               if getattr(call, "genus", "")}
            calls = typer.type_batch(batch, workdir, threads=self.config.threads,
                                     force_scheme=self.options.force_scheme,
                                     genus_by_sample=genus_by_sample,
                                     scheme_genus=scheme_genus)
            for sample, call in calls.items():
                if sample in results:
                    results[sample].mlst = call
        elif self.options.mlst:
            LOG.warning("MLST requested but the pubmlst database is not installed; skipping")

        for sample, result in results.items():
            if sample not in assemblies:
                continue
            if not (self.config.organism or self.config.auto_organism):
                organism_by_sample[sample] = None
                continue
            result.species = identifier.identify(sample, result.mlst, assemblies.get(sample),
                                                 threads=self.config.threads,
                                                 sketch=sketches.get(sample))
            if self.config.organism:
                result.species.organism = self.config.organism
                result.species.evidence = ((result.species.evidence + "; ")
                                           if result.species.evidence else "") + "organism set by user"
                organism_by_sample[sample] = self.config.organism
            else:
                organism_by_sample[sample] = result.species.organism
            if organism_by_sample.get(sample):
                LOG.debug("%s -> organism %s", sample, organism_by_sample[sample])
        return organism_by_sample

    def _screen_assemblies(self, batch: QueryBatch, results: dict[str, SampleResult],
                           organism_by_sample: dict[str, str | None], workdir: Path) -> None:
        nucl_dbs = [name for name in self.options.databases
                    if DATABASES.get(name) and DATABASES[name].kind == "nucl"]
        want_protein = self.options.protein and "protein" in self.options.databases

        if nucl_dbs:
            screener = NucleotideScreener(self.store, self.config)
            for sample, hits in screener.screen(batch, nucl_dbs, workdir,
                                                threads=self.config.threads).items():
                if sample in results:
                    results[sample].hits.extend(hits)

        if want_protein:
            if not self.store.is_installed("protein"):
                LOG.warning("protein database not installed; skipping translated search")
            else:
                protein = ProteinScreener(self.store, self.config)
                for sample, hits in protein.screen(batch, workdir, organism_by_sample,
                                                   threads=self.config.threads).items():
                    if sample in results:
                        results[sample].hits.extend(hits)
                if self.options.point_mutations:
                    self._dna_mutations(batch, results, organism_by_sample, protein, workdir)

        if self.options.typing and self.store.is_installed("lineage"):
            typer = LineageTyper(self.store, self.config)
            species_by_sample = {name: result.species for name, result in results.items()}
            for sample, typing in typer.type_batch(
                    batch, workdir, threads=self.config.threads,
                    species_by_sample=species_by_sample).items():
                if sample in results:
                    results[sample].typing = typing
        elif self.options.typing:
            LOG.debug("lineage database not installed; skipping lineage typing")

    def _dna_mutations(self, batch: QueryBatch, results: dict[str, SampleResult],
                       organism_by_sample: dict[str, str | None], protein: ProteinScreener,
                       workdir: Path) -> None:
        """Run one DNA-mutation search per organism group present in the batch."""
        groups: dict[str, list[str]] = defaultdict(list)
        for sample, organism in organism_by_sample.items():
            if organism:
                groups[organism].append(sample)
        installed = set(self.store.mutation_organisms())
        for organism, samples in sorted(groups.items()):
            if organism not in installed:
                LOG.debug("no DNA mutation reference for %s", organism)
                continue
            sub_batch = self._sub_batch(batch, set(samples), workdir, organism)
            if sub_batch is None:
                continue
            hits = protein.call_dna_mutations(sub_batch, workdir, organism,
                                              threads=self.config.threads)
            for sample, sample_hits in hits.items():
                if sample in results:
                    results[sample].hits.extend(sample_hits)

    @staticmethod
    def _sub_batch(batch: QueryBatch, samples: set[str], workdir: Path,
                   tag: str) -> QueryBatch | None:
        """Slice a batch down to a subset of samples, reusing the combined query FASTA."""
        keep = {qid for qid, piece in batch.id_map.items() if piece.sample in samples}
        if not keep:
            return None
        path = workdir / f"query.{tag}.fna"
        wanted = False
        n_contigs = 0
        total = 0
        with open(batch.path) as source, open(path, "w") as target:
            for line in source:
                if line.startswith(">"):
                    qid = line[1:].strip().split()[0]
                    wanted = qid in keep
                    if wanted:
                        n_contigs += 1
                if wanted:
                    target.write(line)
                    if not line.startswith(">"):
                        total += len(line.strip())
        return QueryBatch(path=path,
                          id_map={q: v for q, v in batch.id_map.items() if q in keep},
                          lengths={q: v for q, v in batch.lengths.items() if q in keep},
                          n_samples=len(samples), n_contigs=n_contigs, total_bases=total)

    def _process_reads(self, readsets: dict[str, ReadSet], results: dict[str, SampleResult],
                       organism_by_sample: dict[str, str | None], workdir: Path) -> None:
        mapper = ReadMapper(self.store, self.config)
        try:
            mapper.check_tools()
        except HydraError as exc:
            LOG.error("skipping read analysis: %s", exc)
            for sample in readsets:
                results[sample].warnings.append(f"read analysis skipped: {exc}")
            return

        bam_dir = workdir / "bam"
        bam_dir.mkdir(parents=True, exist_ok=True)
        nucl_dbs = [name for name in self.options.databases
                    if DATABASES.get(name) and DATABASES[name].kind == "nucl"]
        jobs, per_job_threads = self.config.worker_layout(len(readsets))
        if jobs > 1 and len(readsets) > 1:
            LOG.info("mapping reads for %d samples, %d at a time (%d threads each)",
                     len(readsets), jobs, per_job_threads)
            with ThreadPoolExecutor(max_workers=jobs) as pool:
                futures = {
                    pool.submit(self._process_one_readset, sample, reads, results[sample],
                                organism_by_sample, mapper, bam_dir, nucl_dbs,
                                per_job_threads): sample
                    for sample, reads in readsets.items()
                }
                for future in as_completed(futures):
                    sample = futures[future]
                    try:
                        future.result()
                    except Exception as exc:  # noqa: BLE001 - one sample must not sink the run
                        LOG.error("read analysis of %s failed: %s", sample, exc)
                        results[sample].warnings.append(f"read analysis failed: {exc}")
            return
        for sample, reads in readsets.items():
            try:
                self._process_one_readset(sample, reads, results[sample], organism_by_sample,
                                          mapper, bam_dir, nucl_dbs, self.config.threads)
            except Exception as exc:  # noqa: BLE001 - one sample must not sink the run
                LOG.error("read analysis of %s failed: %s", sample, exc)
                results[sample].warnings.append(f"read analysis failed: {exc}")

    def _process_one_readset(self, sample: str, reads: ReadSet, result: SampleResult,
                             organism_by_sample: dict[str, str | None], mapper: ReadMapper,
                             bam_dir: Path, nucl_dbs: list[str], threads: int) -> None:
        """Everything the reads of one sample contribute."""
        result.qc.update({f"reads_{k}": v for k, v in fastq_stats(reads.r1).items()})
        organism = organism_by_sample.get(sample) or self.config.organism
        if not organism and self.config.auto_organism and result.input_type == "reads":
            # A FASTQ-only sample has no assembly to type, so sketch the reads
            # directly rather than making the user name the organism.
            call = self._identifier().identify_reads(reads.r1, threads=threads)
            if call is not None:
                result.species = call
                organism = call.organism
                LOG.info("%s: reads identified as %s (mash d=%.4f)",
                         sample, call.name, call.distance or 0.0)

        # A sample that also has an assembly already has its acquired genes called
        # from contigs; its reads only add allele-fraction resolution.
        # A sample that also has an assembly already has its genes called from
        # contigs, so its reads are mapped only to add allele fractions.
        reads_only = result.input_type == "reads"
        if nucl_dbs and ((self.options.reads_genes and reads_only)
                         or self.options.reads_variants):
            for db_name in nucl_dbs:
                handle = self.store.handle(db_name)
                bam = mapper.align(reads, handle.fasta, bam_dir / f"{sample}.{db_name}.bam",
                                   threads=threads)
                if self.options.reads_genes and reads_only:
                    result.hits.extend(mapper.call_genes(sample, bam, db_name, handle))
                if self.options.reads_variants:
                    self._variants_from_reads(sample, result, mapper, bam, db_name, handle,
                                              organism, reads_only=reads_only)

        if self.options.reads_mlst and not (result.mlst.scheme != "-"
                                            and result.mlst.source == "assembly"):
            self._mlst_from_reads(sample, reads, result, mapper, bam_dir, threads)

        if not self.options.heteroresistance:
            return
        if not organism:
            result.warnings.append(
                "heteroresistance skipped: no organism determined; "
                "pass --organism or provide an assembly for species identification")
            return
        catalog = MutationCatalog(self.store.root, organism)
        if catalog.dna_reference is None or not catalog.dna:
            LOG.debug("no DNA mutation reference for %s", organism)
            return
        lengths = {header.split()[0]: len(seq)
                   for header, seq in read_fasta(catalog.dna_reference)}
        bam = mapper.align(reads, catalog.dna_reference,
                           bam_dir / f"{sample}.mutations.bam", threads=threads)
        calls = mapper.call_sites(bam, catalog, lengths, organism=organism,
                                  min_base_quality=self.options.min_base_quality)
        result.hits.extend(mapper.site_hits(
            sample, calls, report_absent=self.options.report_absent_sites, organism=organism))
        hetero = [c for c in calls if c.status == "heteroresistant"]
        if hetero:
            LOG.info("%s: %d heteroresistant site(s): %s", sample, len(hetero),
                     ", ".join(f"{c.symbol}@{c.allele_fraction:.0%}" for c in hetero))
        low = [c for c in calls if c.status == "low-depth"]
        if low and len(low) == len(calls):
            result.warnings.append(
                f"all {len(calls)} mutation sites below --min-depth ({self.config.thresholds.min_depth}x); "
                f"heteroresistance cannot be assessed")

    def _mlst_from_reads(self, sample: str, reads: ReadSet, result: SampleResult,
                         mapper: ReadMapper, bam_dir: Path, threads: int) -> None:
        """Type a sample from its reads by mapping them to the scheme's loci."""
        if not self.store.is_installed("pubmlst"):
            return
        typer = MlstTyper(self.store, self.config)
        identifier = self._identifier()
        genus = result.species.genus
        table = identifier.scheme_table
        installed = self.store.root / "mlst" / "pubmlst"
        if self.options.force_scheme:
            schemes = [self.options.force_scheme]
        else:
            genus_match = [name for name, entry in table.items()
                           if genus and entry.genus == genus and (installed / name).is_dir()]
            # Prefer the schemes for this exact species: mapping against a sister
            # species' loci as well only splits the reads.
            species = result.species.species
            exact = [name for name in genus_match
                     if species and table[name].species == species]
            schemes = exact or genus_match
        share_groups = {name: f"{table[name].genus} {table[name].species}".strip()
                        for name in schemes if name in table}
        if not schemes:
            result.warnings.append(
                "read MLST skipped: no scheme matches the species call; "
                "pass --scheme to choose one")
            return
        LOG.info("%s: typing from reads against %s", sample, ", ".join(schemes))
        reference = bam_dir / f"{sample}.mlst.fna"
        index = typer.representative_reference(schemes, reference, share_groups=share_groups)
        sequences = {name: seq for name, seq in read_fasta(reference)}
        bam = mapper.align(reads, reference, bam_dir / f"{sample}.mlst.bam", threads=threads)
        consensus = mapper.consensus(
            bam, sequences, min_depth=self.config.thresholds.min_depth,
            min_base_quality=self.options.min_base_quality,
            min_allele_fraction=self.config.thresholds.min_allele_fraction,
            min_alt_reads=self.config.thresholds.min_allele_reads,
        )
        call = typer.type_from_reads(consensus, index, genus=genus,
                                     force_scheme=self.options.force_scheme)
        if call.scheme == "-" and result.mlst.scheme != "-":
            return
        result.mlst = call
        if "mixed or contaminated" in call.note:
            result.warnings.append(
                "MLST loci carry intermediate allele fractions: this sample looks mixed or "
                "contaminated, and every read-derived call for it should be treated as such")
        LOG.info("%s: read MLST %s/%s (%d/%d loci)", sample, call.scheme, call.sequence_type,
                 call.loci_found, call.loci_total)

    def _variants_from_reads(self, sample: str, result: SampleResult, mapper: ReadMapper,
                             bam: Path, db_name: str, handle, organism: str | None,
                             reads_only: bool = True) -> None:
        """Report mutations in the resistance genes, measured against the closest reference."""
        if reads_only:
            wanted = {hit.sequence for hit in result.hits
                      if hit.database == db_name and hit.method == "READS"}
        else:
            # Genes were called from the assembly, so the references worth
            # inspecting are the ones those calls named.
            called = {hit.gene for hit in result.hits if hit.database == db_name}
            wanted = {seqid for seqid, meta in (handle.meta or {}).items()
                      if meta.get("gene") in called}
        if not wanted:
            return
        sequences = {name: seq for name, seq in read_fasta(handle.fasta) if name in wanted}
        consensus = mapper.consensus(
            bam, sequences, min_depth=self.config.thresholds.min_depth,
            min_base_quality=self.options.min_base_quality,
            min_allele_fraction=self.config.thresholds.min_allele_fraction,
            wanted=wanted, min_alt_reads=self.config.thresholds.min_allele_reads,
            fixed_allele_fraction=self.config.thresholds.fixed_allele_fraction,
        )
        catalog = MutationCatalog(self.store.root, organism) if organism else None
        for reference, locus in consensus.items():
            meta = handle.meta_for(reference)
            result.hits.extend(mapper.variant_hits(
                sample, db_name, locus, meta, catalog=catalog,
                synonymous=self.options.report_synonymous))

    #: Databases ranked by how much they can say about a hit. AMRFinderPlus and
    #: NCBI carry curated class/subclass annotation, so their call wins when
    #: several databases describe the same locus.
    DATABASE_RANK = {"protein": 0, "ncbi": 1, "resfinder": 2, "card": 3,
                     "argannot": 4, "megares": 5, "vfdb": 1, "vfdb_full": 2,
                     "ecoli_vf": 3, "plasmidfinder": 1, "ecoh": 1}

    def _finalise(self, result: SampleResult) -> None:
        """Deduplicate hits, resolve cross-database redundancy, and score."""
        seen: set[tuple] = set()
        unique: list[Hit] = []
        for hit in result.hits:
            key = (hit.database, hit.gene, hit.sequence, hit.start, hit.end, hit.method, hit.note)
            if key in seen:
                continue
            seen.add(key)
            unique.append(hit)
        result.hits = unique
        self._resolve_loci(result)

        score, flags = resistance_score([h for h in result.hits if h.primary])
        result.scores["resistance_score"] = score
        result.scores["has_esbl"] = flags["esbl"]
        result.scores["has_carbapenemase"] = flags["carbapenemase"]
        result.scores["has_colistin_resistance"] = flags["colistin"]
        if result.typing:
            virulence, basis = virulence_score(result.species, result.typing)
            result.scores["virulence_score"] = virulence
            result.scores["virulence_score_basis"] = basis

    def _resolve_loci(self, result: SampleResult) -> None:
        """Mark one hit per genomic locus as primary when databases overlap.

        Screening several databases finds the same gene several times: abricate
        users are used to seeing every row, but counting them all would inflate
        every summary. The rows are kept and the duplicates flagged instead.
        Point mutations and read-derived calls are never demoted - they describe
        something an acquired-gene hit at the same coordinates does not.
        """
        if self.config.report_overlaps:
            return
        clusters: dict[tuple[str, str], list[Hit]] = {}
        for hit in result.hits:
            if hit.resolution == "POINT" or hit.method in ("READS", "POINTR"):
                continue
            clusters.setdefault((hit.sequence, hit.element_type), []).append(hit)
        for hits in clusters.values():
            hits.sort(key=lambda h: (h.start, h.end))
            claimed: list[tuple[int, int, Hit]] = []
            for hit in hits:
                span = (min(hit.start, hit.end), max(hit.start, hit.end))
                length = max(1, span[1] - span[0] + 1)
                overlapping = None
                for index, (lo, hi, other) in enumerate(claimed):
                    shared = min(hi, span[1]) - max(lo, span[0]) + 1
                    other_length = max(1, hi - lo + 1)
                    if shared >= 0.6 * min(length, other_length):
                        overlapping = index
                        break
                if overlapping is None:
                    claimed.append((span[0], span[1], hit))
                    continue
                lo, hi, other = claimed[overlapping]
                if self._prefer(hit, other):
                    other.primary = False
                    claimed[overlapping] = (span[0], span[1], hit)
                else:
                    hit.primary = False

    def _prefer(self, candidate: Hit, current: Hit) -> bool:
        """True when *candidate* is the better representative of a shared locus."""
        def rank(hit: Hit) -> tuple:
            return (
                -self.DATABASE_RANK.get(hit.database, 9),
                1 if hit.resolution == "COMPLETE" else 0,
                1 if hit.drug_class else 0,
                hit.coverage_pct,
                hit.identity_pct,
            )
        return rank(candidate) > rank(current)
