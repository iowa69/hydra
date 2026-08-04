"""Acquired-gene screening against nucleotide databases.

This is Hydra's abricate-equivalent, with two differences that matter in
practice: every sample in a run is screened in a *single* BLAST invocation per
database (one process start, one database load, one index scan instead of N),
and HSPs are merged on the reference axis so a gene broken across a contig
boundary is still reported at its true coverage.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..config import Config
from ..db.manager import DatabaseStore
from ..db.registry import spec_for
from ..records import Hit
from ..seqio import read_fasta
from ..utils import LOG, HydraError, natural_key
from .blast import MergedHit, blast, check_db_exists, deduplicate, merge_hsps


#: Contigs longer than this are split so BLAST can spread them across threads.
DEFAULT_CHUNK_SIZE = 300_000
#: Neighbouring chunks share this much sequence, so no gene falls in a seam.
DEFAULT_CHUNK_OVERLAP = 20_000


@dataclass(frozen=True)
class QueryPiece:
    """Where one BLAST query sequence came from."""

    sample: str
    contig: str
    offset: int = 0   # 0-based start of this piece within the original contig


@dataclass
class QueryBatch:
    """A combined FASTA over many samples, with a map back to the originals."""

    path: Path
    id_map: dict[str, QueryPiece]        # blast id -> origin
    lengths: dict[str, int]              # blast id -> piece length
    n_samples: int
    n_contigs: int
    total_bases: int

    def samples(self) -> set[str]:
        return {piece.sample for piece in self.id_map.values()}


def build_query_batch(samples: dict[str, Path], out_path: Path,
                      min_contig_length: int = 0,
                      chunk_size: int = DEFAULT_CHUNK_SIZE,
                      chunk_overlap: int = DEFAULT_CHUNK_OVERLAP) -> QueryBatch:
    """Concatenate the assemblies of many samples into one BLAST query.

    Contig names are replaced with short synthetic ids so BLAST parsing stays
    cheap and no original header can collide across samples. Long contigs are
    split into overlapping chunks: BLAST parallelises over query *sequences*, so
    a closed 5 Mb chromosome would otherwise pin the whole search to one thread.
    The overlap is far longer than any resistance gene, so every gene is intact
    in at least one chunk and the duplicate copies are removed downstream.
    """
    if chunk_size and chunk_overlap >= chunk_size:
        # Chunking is an internal performance detail, so a chunk size smaller
        # than the default overlap adapts rather than failing the run.
        chunk_overlap = max(0, chunk_size // 4)
        LOG.debug("chunk overlap reduced to %d to fit a chunk size of %d",
                  chunk_overlap, chunk_size)
    id_map: dict[str, QueryPiece] = {}
    lengths: dict[str, int] = {}
    n_contigs = 0
    total = 0
    with open(out_path, "w") as handle:
        for s_index, (sample, fasta) in enumerate(samples.items()):
            wrote_any = False
            for c_index, (header, seq) in enumerate(read_fasta(fasta)):
                if not seq or len(seq) < min_contig_length:
                    continue
                original = header.split()[0] if header.split() else header
                n_contigs += 1
                total += len(seq)
                wrote_any = True
                if not chunk_size or len(seq) <= chunk_size:
                    qid = f"q{s_index}_{c_index}"
                    id_map[qid] = QueryPiece(sample, original, 0)
                    lengths[qid] = len(seq)
                    handle.write(f">{qid}\n{seq}\n")
                    continue
                step = chunk_size - chunk_overlap
                for p_index, start in enumerate(range(0, len(seq), step)):
                    piece = seq[start:start + chunk_size]
                    if not piece:
                        break
                    qid = f"q{s_index}_{c_index}_{p_index}"
                    id_map[qid] = QueryPiece(sample, original, start)
                    lengths[qid] = len(piece)
                    handle.write(f">{qid}\n{piece}\n")
                    if start + chunk_size >= len(seq):
                        break
            if not wrote_any:
                LOG.warning("sample %s contributed no contigs (empty or all below "
                            "--min-contig-length)", sample)
    return QueryBatch(path=out_path, id_map=id_map, lengths=lengths,
                      n_samples=len(samples), n_contigs=n_contigs, total_bases=total)


class NucleotideScreener:
    """Screens a batch of assemblies against one or more nucleotide databases."""

    def __init__(self, store: DatabaseStore, config: Config):
        self.store = store
        self.config = config

    def screen(self, batch: QueryBatch, db_names: list[str], workdir: Path,
               threads: int | None = None) -> dict[str, list[Hit]]:
        """Return sample -> hits for every requested nucleotide database."""
        results: dict[str, list[Hit]] = {sample: [] for sample in batch.samples()}
        threads = threads or self.config.threads
        for name in db_names:
            try:
                hits = self._screen_one(batch, name, workdir, threads)
            except HydraError as exc:
                LOG.error("database %s: %s", name, exc)
                raise
            for sample, sample_hits in hits.items():
                results.setdefault(sample, []).extend(sample_hits)
        return results

    def _screen_one(self, batch: QueryBatch, name: str, workdir: Path,
                    threads: int) -> dict[str, list[Hit]]:
        handle = self.store.handle(name)
        if handle.kind != "nucl":
            raise HydraError(f"database '{name}' is not a nucleotide database")
        check_db_exists(handle.fasta, "nucl")
        spec = spec_for(name) if name in ("ncbi", "card", "vfdb", "vfdb_full", "resfinder",
                                          "argannot", "megares", "plasmidfinder", "ecoh",
                                          "ecoli_vf") else None
        thresholds = self.config.thresholds
        min_identity = thresholds.min_identity
        min_coverage = thresholds.min_coverage
        # Some databases need stricter cut-offs than the global default (replicon
        # typing, for one), but only when the user has not asked for a value.
        if spec is not None:
            if spec.default_identity is not None and not thresholds.was_set("min_identity"):
                min_identity = spec.default_identity
            if spec.default_coverage is not None and not thresholds.was_set("min_coverage"):
                min_coverage = spec.default_coverage

        out_tab = workdir / f"{name}.blastn.tsv"
        LOG.info("screening %d contigs from %d samples against %s (%s sequences)",
                 batch.n_contigs, batch.n_samples, name, handle.n_sequences or "?")
        hsps = blast(
            "blastn", batch.path, handle.fasta, out_tab,
            threads=threads,
            evalue=1e-20,
            perc_identity=max(0.0, min_identity - 5.0),  # merge can raise identity; filter later
            task="blastn",
            max_target_seqs=10000,
        )
        merged = merge_hsps(hsps)
        LOG.debug("%s: %d HSPs -> %d merged query/subject pairs", name, len(hsps), len(merged))

        per_sample: dict[str, list[tuple[MergedHit, dict]]] = {}
        for hit in merged:
            if hit.identity_pct < min_identity or hit.coverage_pct < min_coverage:
                continue
            piece = batch.id_map.get(hit.qseqid)
            if piece is None:
                continue
            meta = handle.meta_for(hit.sseqid)
            per_sample.setdefault(piece.sample, []).append(
                (hit, {"contig": piece.contig, "offset": piece.offset, **meta}))

        out: dict[str, list[Hit]] = {}
        for sample, entries in per_sample.items():
            # Deduplicate on true contig coordinates: chunking can surface the same
            # locus twice when it falls inside the overlap between two chunks.
            entries = deduplicate(
                entries,
                key_span=lambda e: (e[0].query_start + e[1].get("offset", 0),
                                    e[0].query_end + e[1].get("offset", 0)),
                key_seq=lambda e: e[1]["contig"],
                key_score=lambda e: (round(e[0].match_score, 4), round(e[0].identity_pct, 4),
                                     e[0].bitscore),
                key_tiebreak=lambda e: natural_key(e[1].get("gene", "")),
                overlap_fraction=1.0 if self.config.report_overlaps
                else self.config.overlap_fraction,
            )
            hits: list[Hit] = []
            for merged_hit, meta in entries:
                partial = merged_hit.coverage_pct < thresholds.partial_coverage
                offset = meta.get("offset", 0)
                hits.append(Hit(
                    sample=sample,
                    database=name,
                    gene=meta.get("gene") or merged_hit.sseqid,
                    accession=meta.get("accession", ""),
                    product=meta.get("product", ""),
                    element_type=meta.get("element_type") or handle.element_type,
                    element_subtype=meta.get("element_subtype", ""),
                    drug_class=meta.get("class", ""),
                    subclass=meta.get("subclass", ""),
                    sequence=meta["contig"],
                    start=merged_hit.query_start + offset,
                    end=merged_hit.query_end + offset,
                    strand=merged_hit.strand,
                    coverage=merged_hit.coverage_string,
                    coverage_pct=merged_hit.coverage_pct,
                    identity_pct=merged_hit.identity_pct,
                    gaps=merged_hit.gaps,
                    bitscore=merged_hit.bitscore,
                    method="PARTIALN" if partial else "BLASTN",
                    resolution="PARTIAL" if partial else "COMPLETE",
                    note="" if merged_hit.n_hsps == 1 else f"{merged_hit.n_hsps} alignment blocks",
                ))
            hits.sort(key=lambda h: (h.sequence, h.start))
            out[sample] = hits
        total = sum(len(v) for v in out.values())
        LOG.info("%s: %d hits across %d samples", name, total, len(out))
        return out
