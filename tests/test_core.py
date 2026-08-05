"""Unit tests for the parts of Hydra that do not need a database or BLAST."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from hydra_amr.cli import collect_inputs
from hydra_amr.db.manager import _infer_class, _normalise_gene, parse_amrprot_header
from hydra_amr.db.fetch import (can_fetch, card_header, cge_header, ncbi_header,
                                pubmlst_organism, pubmlst_scheme_name, vfdb_header)
from hydra_amr.db.registry import DB_GROUPS, resolve_names
from hydra_amr.engines.blast import (Hsp, deduplicate, interval_length, merge_hsps,
                                     merge_intervals)
from hydra_amr.engines.mutations import (AlignedObservation, MutationEntry, _parse_symbol,
                                         evaluate, reference_index, walk_alignment)
from hydra_amr.engines.nucl import build_query_batch
from hydra_amr.engines.reads import (LocusConsensus, ReadMapper, Variant, _annotate_codons,
                                     _base_counts, _count_alt, _indel_count, _looks_like_cds,
                                     binomial_upper_tail, pair_reads)
from hydra_amr.pipeline import Pipeline
from hydra_amr.records import Hit, SampleResult, SpeciesCall, TypingResult
from hydra_amr.typing.lineage import _marker_group, resistance_score, virulence_score
from hydra_amr.typing.species import SchemeOrganism, _inherit_organism_by_species
from hydra_amr.typing.mlst import MlstTyper, SchemeProfiles
from hydra_amr.report.html import _format_cell
from hydra_amr.report.tables import (GENE_COLUMNS, _coverage_glyph, class_summary,
                                     gene_table, long_table, matrix, summary_table)
from hydra_amr.report.writer import write_outputs
from hydra_amr.utils import HydraError, natural_key, revcomp, sample_name_from_path, translate


# --------------------------------------------------------------------- utils
def test_revcomp_handles_ambiguity():
    assert revcomp("ACGT") == "ACGT"
    assert revcomp("AACCGGTT") == "AACCGGTT"
    assert revcomp("ATGCN") == "NGCAT"


def test_translate_uses_bacterial_start():
    # GTG is a valid bacterial start codon and is rendered as methionine.
    assert translate("GTGAAATAA").startswith("M")
    assert translate("ATGAAATAA") == "MK*"


@pytest.mark.parametrize("path,expected", [
    ("sample.fasta", "sample"),
    ("sample.fna.gz", "sample"),
    ("/tmp/dir/isolate.contigs.fasta", "isolate.contigs"),
    ("reads_R1_001.fastq.gz", "reads_R1_001"),
])
def test_sample_name_from_path(path, expected):
    assert sample_name_from_path(path) == expected


def test_sample_name_strips_read_suffix():
    assert sample_name_from_path("x_R1_001.fastq.gz", strip_read_suffix=True) == "x"
    assert sample_name_from_path("x_1.fq.gz", strip_read_suffix=True) == "x"


# ------------------------------------------------------------------ database
def test_parse_amrprot_header_current_format():
    header = ("WP_004199234.1|1|1|blaKPC-2|blaKPC|hydrolase|2|CARBAPENEM|BETA-LACTAM|"
              "carbapenem-hydrolyzing_class_A_beta-lactamase_KPC-2")
    info = parse_amrprot_header(header)
    assert info["accession"] == "WP_004199234.1"
    assert info["gene"] == "blaKPC-2"
    assert info["class"] == "BETA-LACTAM"
    assert info["subclass"] == "CARBAPENEM"
    assert info["part_of_gene"] == "hydrolase"


def test_parse_amrprot_header_legacy_gi_prefix():
    """Databases before 2026 prefixed the header with a numeric GI."""
    header = ("0|WP_004199234.1|1|1|blaKPC-2|blaKPC|hydrolase|2|CARBAPENEM|BETA-LACTAM|"
              "carbapenem-hydrolyzing_class_A_beta-lactamase_KPC-2")
    info = parse_amrprot_header(header)
    assert info["accession"] == "WP_004199234.1"
    assert info["class"] == "BETA-LACTAM"
    assert info["subclass"] == "CARBAPENEM"


def test_parse_amrprot_header_marks_mutation_references():
    header = "WP_000116450.1|1|1|gyrA|gyrA|mutation|2|||DNA_gyrase_subunit_A_GyrA"
    assert parse_amrprot_header(header)["part_of_gene"] == "mutation"


def test_normalise_gene_strips_resfinder_suffixes():
    assert _normalise_gene("blaKPC-2_1_AY034847") == _normalise_gene("blaKPC-2")


def test_infer_class_from_product_text():
    assert _infer_class("carbapenem-hydrolyzing beta-lactamase") == "BETA-LACTAM"
    assert _infer_class("linezolid resistance protein") == "OXAZOLIDINONE"
    assert _infer_class("hypothetical protein") == ""


def test_resolve_names_expands_groups_and_aliases():
    assert resolve_names(["ncbi"]) == ["ncbi"]
    assert set(resolve_names(["amr"])) == set(DB_GROUPS["amr"])
    assert resolve_names(["kleborate"]) == ["lineage"]
    assert resolve_names(["ncbi,vfdb"]) == ["ncbi", "vfdb"]


def test_resolve_names_rejects_unknown():
    with pytest.raises(HydraError, match="unknown database"):
        resolve_names(["not-a-database"])


# --------------------------------------------------------------------- blast
def _hsp(**kwargs) -> Hsp:
    defaults = dict(qseqid="q", sseqid="s", pident=100.0, length=100, mismatch=0, gapopen=0,
                    qstart=1, qend=100, sstart=1, send=100, evalue=0.0, bitscore=200.0,
                    slen=100, qlen=1000, gaps=0, nident=100, strand="+")
    defaults.update(kwargs)
    return Hsp(**defaults)


def test_merge_intervals():
    assert merge_intervals([(1, 10), (5, 20), (30, 40)]) == [(1, 20), (30, 40)]
    assert merge_intervals([(10, 20), (21, 30)]) == [(10, 30)]
    assert interval_length([(1, 10), (20, 29)]) == 20


def test_merge_hsps_combines_split_gene():
    """A gene split across an assembly gap still reaches full reference coverage."""
    a = _hsp(qstart=1, qend=50, sstart=1, send=50, length=50, nident=50, slen=100, bitscore=100)
    b = _hsp(qstart=60, qend=109, sstart=51, send=100, length=50, nident=50, slen=100, bitscore=100)
    merged = merge_hsps([a, b])
    assert len(merged) == 1
    assert merged[0].coverage_pct == 100.0
    assert merged[0].identity_pct == 100.0
    assert merged[0].n_hsps == 2


def test_merge_hsps_keeps_strands_apart():
    plus = _hsp(sstart=1, send=100, strand="+", bitscore=200)
    minus = _hsp(sstart=100, send=1, strand="-", bitscore=50)
    merged = merge_hsps([plus, minus])
    assert len(merged) == 1
    assert merged[0].strand == "+"


def test_merge_hsps_refuses_to_join_distant_query_fragments():
    """Two truncated copies far apart on a contig are not one complete gene."""
    left = _hsp(qstart=10_000, qend=10_450, sstart=1, send=450, length=450, nident=450,
                slen=900, bitscore=800)
    right = _hsp(qstart=250_000, qend=250_450, sstart=451, send=900, length=450, nident=450,
                 slen=900, bitscore=790)
    merged = merge_hsps([left, right])
    # Two separate partial hits, each half the reference - not one full-length gene.
    assert len(merged) == 2
    assert all(m.n_hsps == 1 for m in merged)
    assert all(m.coverage_pct == pytest.approx(50.0) for m in merged)


def test_merge_hsps_still_joins_a_gene_split_by_an_assembly_gap():
    left = _hsp(qstart=1000, qend=1450, sstart=1, send=450, length=450, nident=450,
                slen=900, bitscore=800)
    right = _hsp(qstart=1500, qend=1950, sstart=451, send=900, length=450, nident=450,
                 slen=900, bitscore=790)
    merged = merge_hsps([left, right])
    assert merged[0].n_hsps == 2
    assert merged[0].coverage_pct == pytest.approx(100.0)


def test_merge_hsps_does_not_double_count_overlapping_identity():
    a = _hsp(sstart=1, send=100, length=100, nident=100, slen=200, bitscore=200)
    b = _hsp(sstart=80, send=200, length=121, nident=60, slen=200, bitscore=100)
    merged = merge_hsps([a, b])[0]
    assert merged.coverage_pct == pytest.approx(100.0)
    # The shared 21 bases must not be counted once as identical and once as not.
    assert 45.0 < merged.identity_pct < 100.0


def test_merge_hsps_identity_is_length_weighted():
    a = _hsp(sstart=1, send=50, length=50, nident=50, slen=100, bitscore=100)
    b = _hsp(sstart=51, send=100, length=50, nident=40, slen=100, bitscore=90)
    merged = merge_hsps([a, b])[0]
    assert merged.identity_pct == pytest.approx(90.0)


def test_match_score_is_not_length_biased():
    """A shorter, better-matching allele must beat a longer, worse one."""
    short = merge_hsps([_hsp(sseqid="dfrA1", sstart=1, send=474, length=474, nident=473,
                             slen=474, bitscore=800)])[0]
    long = merge_hsps([_hsp(sseqid="dfr1_rpt", sstart=1, send=600, length=600, nident=593,
                            slen=600, bitscore=1000)])[0]
    assert long.bitscore > short.bitscore          # bit score prefers the longer one
    assert short.match_score > long.match_score    # match score prefers the better one


def test_natural_key_orders_allele_numbers():
    names = ["blaSHV-11", "blaSHV-2", "blaSHV-67"]
    assert sorted(names, key=natural_key) == ["blaSHV-2", "blaSHV-11", "blaSHV-67"]


def test_deduplicate_tiebreak_is_deterministic():
    """Equal-scoring alleles must resolve the same way regardless of input order."""
    forward = [("blaSHV-67", 1, 100, 10.0), ("blaSHV-11", 1, 100, 10.0)]
    reverse = list(reversed(forward))
    picks = []
    for items in (forward, reverse):
        kept = deduplicate(items, key_span=lambda x: (x[1], x[2]),
                           key_seq=lambda _x: "c", key_score=lambda x: x[3],
                           key_tiebreak=lambda x: natural_key(x[0]),
                           overlap_fraction=0.5)
        picks.append(kept[0][0])
    assert picks == ["blaSHV-11", "blaSHV-11"]


def test_deduplicate_keeps_best_of_overlapping():
    items = [("a", 1, 100, 10.0), ("b", 5, 95, 20.0), ("c", 500, 600, 5.0)]
    kept = deduplicate(items, key_span=lambda x: (x[1], x[2]), key_seq=lambda _x: "contig",
                       key_score=lambda x: x[3], overlap_fraction=0.5)
    assert {k[0] for k in kept} == {"b", "c"}


# ----------------------------------------------------------------- mutations
def test_parse_symbol_variants():
    assert _parse_symbol("gyrA_S83L", 83) == ("gyrA", "S", "L")
    assert _parse_symbol("23S_G2576T", 2576) == ("23S", "G", "T")
    assert _parse_symbol("pbp4_T-266A", -266) == ("pbp4", "T", "A")
    gene, ref, alt = _parse_symbol("mgrB_Q30STOP", 30)
    assert (gene, ref, alt) == ("mgrB", "Q", "*")


def test_reference_index_negative_positions_count_from_the_end():
    """Promoter coordinates are negative and measured back from the record end."""
    assert reference_index(1, 301) == 0
    assert reference_index(301, 301) == 300
    assert reference_index(-266, 301) == 35


def test_walk_alignment_reads_query_residue_at_subject_position():
    #            subject 1..5
    qseq = "ACGTA"
    sseq = "ACCTA"
    seen = walk_alignment(qseq, sseq, 1, 5, {3: None}, nucleotide=True)
    assert seen[3].reference == "C"
    assert seen[3].observed == "G"


def test_walk_alignment_handles_subject_gaps():
    qseq = "AC-GT"
    sseq = "ACGGT"
    seen = walk_alignment(qseq, sseq, 1, 5, {3: None, 4: None}, nucleotide=True)
    assert seen[3].observed == "-"
    assert seen[4].observed == "G"


def test_walk_alignment_complements_minus_strand():
    """BLAST prints minus-strand subjects complemented and descending."""
    # Subject 5..1 on the minus strand: displayed letters are complements.
    qseq = "TGCAT"
    sseq = "TGCAT"
    seen = walk_alignment(qseq, sseq, 5, 1, {5: None}, nucleotide=True)
    assert seen[5].reference == "A"   # complement of displayed T
    assert seen[5].observed == "A"


def _entry(**kwargs) -> MutationEntry:
    defaults = dict(kind="dna", taxgroup="X", key="k", position=10, symbol="g_A10T",
                    reported_symbol="g_A10T", ref="A", alt="T", drug_class="C",
                    subclass="S", name="n", gene="g")
    defaults.update(kwargs)
    return MutationEntry(**defaults)


def test_evaluate_calls_matching_alternate_allele():
    obs = AlignedObservation(observed="T", reference="A", aligned=True)
    assert evaluate(_entry(), obs)[0] is True


def test_evaluate_rejects_wild_type():
    obs = AlignedObservation(observed="A", reference="A", aligned=True)
    assert evaluate(_entry(), obs)[0] is False


def test_evaluate_refuses_when_reference_disagrees():
    """A misplaced alignment must not be allowed to call a mutation."""
    obs = AlignedObservation(observed="T", reference="G", aligned=True)
    called, note = evaluate(_entry(), obs)
    assert called is False
    assert "reference mismatch" in note


def test_evaluate_requires_coverage():
    obs = AlignedObservation(observed="", reference="", aligned=False)
    assert evaluate(_entry(), obs)[0] is False


def test_evaluate_deletion():
    obs = AlignedObservation(observed="--", reference="WR", aligned=True)
    assert evaluate(_entry(ref="WR", alt="del"), obs)[0] is True


def test_nonsense_mutations_use_the_catalogue_spelling():
    """The shipped catalogue writes stop codons as 'Ter', older ones as 'STOP'."""
    assert _parse_symbol("nfsB_W94Ter", 94) == ("nfsB", "W", "*")
    assert _parse_symbol("mgrB_Q30STOP", 30) == ("mgrB", "Q", "*")
    assert _entry(ref="W", alt="*").is_stop is True
    assert _entry(ref="W", alt="*").change_type == "nonsense"


@pytest.mark.parametrize("ref,alt,expected", [
    ("G", "T", "substitution"),
    ("G", "GG", "insertion"),      # ampC_G-15GG inserts a base
    ("T", "TGT", "insertion"),     # ampC_T-14TGT
    ("A", "del", "deletion"),
    ("W", "*", "nonsense"),
    ("G", "G", "complex"),         # would otherwise call every wild type mutant
])
def test_change_type_classification(ref, alt, expected):
    assert _entry(ref=ref, alt=alt).change_type == expected


def test_insertion_entry_is_never_called_from_a_substitution_view():
    """ampC_G-15GG must not match a wild-type G at that position."""
    entry = _entry(ref="G", alt="GG")
    obs = AlignedObservation(observed="G", reference="G", aligned=True)
    called, note = evaluate(entry, obs)
    assert called is False
    assert "insertion" in note


def test_inserted_bases():
    assert _entry(ref="T", alt="TGT").inserted_bases == "GT"
    assert _entry(ref="G", alt="T").inserted_bases == ""


def test_count_alt_substitution():
    tokens = ["A", "A", "T", "T", "T"]
    assert _count_alt(tokens, _entry(ref="A", alt="T")) == 3


def test_count_alt_insertion_requires_the_inserted_bases():
    """A wild-type pileup must not support an insertion entry."""
    wild = ["G", "G", "G", "G"]
    mutant = ["G+1G", "G+1G", "G", "G"]
    entry = _entry(ref="G", alt="GG")
    assert _count_alt(wild, entry) == 0
    assert _count_alt(mutant, entry) == 2


def test_count_alt_deletion_counts_the_placeholder():
    assert _count_alt(["A", "*", "*", "A"], _entry(ref="A", alt="del")) == 2


def test_count_alt_returns_none_for_uncallable_changes():
    assert _count_alt(["A"], _entry(ref="A", alt="A")) is None


def test_base_counts_ignore_indel_suffixes():
    assert _base_counts(["A+2GT", "A", "T-1N", "*"]) == {"A": 2, "T": 1}


# --------------------------------------------------------------------- reads
def test_binomial_upper_tail_is_monotonic_and_finite():
    assert binomial_upper_tail(0, 100, 0.002) == 1.0
    high = binomial_upper_tail(1, 500, 0.002)
    low = binomial_upper_tail(50, 500, 0.002)
    assert 0.0 < low < high <= 1.0


def test_binomial_upper_tail_does_not_saturate_in_the_far_tail():
    """The far tail must stay distinguishable, not collapse to floating-point noise."""
    a = binomial_upper_tail(21, 464, 0.002)
    b = binomial_upper_tail(93, 464, 0.002)
    assert b < a
    assert b > 0.0
    assert math.log10(a) - math.log10(b) > 50


def test_pair_reads_groups_by_convention(tmp_path):
    names = ["s1_R1_001.fastq.gz", "s1_R2_001.fastq.gz", "s2_1.fq.gz", "s2_2.fq.gz"]
    paths = []
    for name in names:
        path = tmp_path / name
        path.write_bytes(b"")
        paths.append(path)
    sets = {rs.sample: rs for rs in pair_reads(paths)}
    assert set(sets) == {"s1", "s2"}
    assert sets["s1"].r2 is not None
    assert sets["s2"].r2 is not None


# --------------------------------------------------------------------- batch
def _write_fasta(path, records):
    with open(path, "w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n{seq}\n")


def test_build_query_batch_chunks_long_contigs(tmp_path):
    fasta = tmp_path / "a.fasta"
    _write_fasta(fasta, [("chr", "A" * 1000)])
    batch = build_query_batch({"a": fasta}, tmp_path / "q.fna",
                              chunk_size=400, chunk_overlap=100)
    assert batch.n_contigs == 1
    assert len(batch.id_map) > 1
    offsets = sorted(p.offset for p in batch.id_map.values())
    assert offsets[0] == 0
    # Every chunk must overlap its neighbour by the requested amount.
    assert offsets[1] == 300
    assert all(p.contig == "chr" for p in batch.id_map.values())


def test_build_query_batch_leaves_short_contigs_whole(tmp_path):
    fasta = tmp_path / "a.fasta"
    _write_fasta(fasta, [("c1", "ACGT" * 10), ("c2", "TTTT" * 10)])
    batch = build_query_batch({"a": fasta}, tmp_path / "q.fna", chunk_size=1000)
    assert len(batch.id_map) == 2
    assert all(p.offset == 0 for p in batch.id_map.values())


def test_build_query_batch_honours_min_contig_length(tmp_path):
    fasta = tmp_path / "a.fasta"
    _write_fasta(fasta, [("short", "ACGT"), ("long", "A" * 500)])
    batch = build_query_batch({"a": fasta}, tmp_path / "q.fna", min_contig_length=100)
    assert batch.n_contigs == 1


def test_build_query_batch_clamps_oversized_overlap(tmp_path):
    """A chunk size below the default overlap adapts instead of failing."""
    fasta = tmp_path / "a.fasta"
    _write_fasta(fasta, [("c", "A" * 100)])
    batch = build_query_batch({"a": fasta}, tmp_path / "q.fna",
                              chunk_size=40, chunk_overlap=100)
    offsets = sorted(p.offset for p in batch.id_map.values())
    assert offsets[0] == 0
    assert offsets[1] == 30  # 40 - (40 // 4)
    assert batch.n_contigs == 1


# -------------------------------------------------------------------- report
def _result(sample: str, genes, element_type="AMR", database="ncbi") -> SampleResult:
    result = SampleResult(sample=sample)
    for index, gene in enumerate(genes):
        result.hits.append(Hit(sample=sample, database=database, gene=gene,
                               element_type=element_type, drug_class="BETA-LACTAM",
                               identity_pct=99.0, coverage_pct=100.0,
                               sequence="c1", start=index * 1000 + 1, end=index * 1000 + 500))
    return result


def test_matrix_binary_and_identity():
    results = [_result("s1", ["blaKPC-2", "tet(A)"]), _result("s2", ["tet(A)"])]
    binary = matrix(results, cell="binary")
    assert binary.loc["s1", "blaKPC-2"] == 1
    assert binary.loc["s2", "blaKPC-2"] == 0
    identity = matrix(results, cell="identity")
    assert identity.loc["s1", "tet(A)"] == pytest.approx(99.0)


def test_matrix_rejects_unknown_axis_and_cell():
    results = [_result("s1", ["a"])]
    with pytest.raises(HydraError):
        matrix(results, columns="nonsense")
    with pytest.raises(HydraError):
        matrix(results, cell="nonsense")


def test_matrix_excludes_redundant_hits():
    """A locus reported by two databases counts once."""
    result = _result("s1", ["blaKPC-2"])
    duplicate = Hit(sample="s1", database="card", gene="blaKPC-2", element_type="AMR",
                    sequence="c1", start=1, end=500, identity_pct=98.0, coverage_pct=100.0)
    duplicate.primary = False
    result.hits.append(duplicate)
    counts = matrix([result], cell="count")
    assert counts.loc["s1", "blaKPC-2"] == 1
    with_redundant = matrix([result], cell="count", primary_only=False)
    assert with_redundant.loc["s1", "blaKPC-2"] == 2


def test_long_and_summary_tables_have_stable_columns():
    results = [_result("s1", ["a", "b"]), _result("s2", [])]
    long = long_table(results)
    assert list(long.columns)[:4] == ["sample", "database", "element_type", "element_subtype"]
    assert len(long) == 2
    summary = summary_table(results)
    assert set(summary["sample"]) == {"s1", "s2"}
    assert summary.set_index("sample").loc["s1", "amr_genes"] == 2
    assert summary.set_index("sample").loc["s2", "amr_genes"] == 0


# ----------------------------------------------------------------- input names
class _Args:
    def __init__(self, **kwargs):
        defaults = dict(inputs=[], assemblies=[], r1=[], r2=[], reads=[],
                        input_list=None, names=[])
        defaults.update(kwargs)
        self.__dict__.update(defaults)


def _touch(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def test_name_renames_samples_in_order(tmp_path):
    """--name applies positionally to the assemblies then the read sets."""
    fasta = _touch(tmp_path / "asm_v2.fasta", ">c\nACGT\n")
    r1 = _touch(tmp_path / "sample_R1.fastq", "@r\nACGT\n+\nIIII\n")
    r2 = _touch(tmp_path / "sample_R2.fastq", "@r\nACGT\n+\nIIII\n")
    assemblies, readsets = collect_inputs(
        _Args(assemblies=[fasta], r1=[r1], r2=[r2], names=["isolate1"]))
    assert list(assemblies) == ["isolate1"]
    assert list(readsets) == ["sample"]


def test_name_can_pair_an_assembly_with_its_reads(tmp_path):
    """Giving both the same name is how a sample gets assembly+reads together."""
    fasta = _touch(tmp_path / "asm_v2.fasta", ">c\nACGT\n")
    r1 = _touch(tmp_path / "sample_R1.fastq", "@r\nACGT\n+\nIIII\n")
    r2 = _touch(tmp_path / "sample_R2.fastq", "@r\nACGT\n+\nIIII\n")
    assemblies, readsets = collect_inputs(
        _Args(assemblies=[fasta], r1=[r1], r2=[r2], names=["isolate1", "isolate1"]))
    assert list(assemblies) == ["isolate1"]
    assert list(readsets) == ["isolate1"]


def test_name_still_refuses_to_collapse_two_assemblies(tmp_path):
    first = _touch(tmp_path / "a.fasta", ">c\nACGT\n")
    second = _touch(tmp_path / "b.fasta", ">c\nACGT\n")
    with pytest.raises(HydraError, match="same name"):
        collect_inputs(_Args(assemblies=[first, second], names=["b"]))


# ------------------------------------------------------- read-derived variants
def _variant(**kwargs) -> Variant:
    defaults = dict(reference="ref", position=10, ref_base="G", alt_base="T",
                    depth=100, alt_count=20, allele_fraction=0.2)
    defaults.update(kwargs)
    return Variant(**defaults)


def test_variant_protein_change_from_codon_annotation():
    consensus = LocusConsensus(reference="gene", sequence="A" * 9, depth=50.0, breadth=100.0,
                               variants=[_variant(position=5, ref_base="A", alt_base="T")])
    _annotate_codons(consensus, "ATGGCATTT")
    variant = consensus.variants[0]
    assert variant.codon_position == 2
    assert variant.ref_aa == "A"          # GCA
    assert variant.alt_aa == "V"          # GTA
    assert variant.protein_change == "A2V"
    assert variant.is_synonymous is False


def test_synonymous_change_is_flagged():
    consensus = LocusConsensus(reference="g", sequence="A" * 6, depth=50.0, breadth=100.0,
                               variants=[_variant(position=6, ref_base="A", alt_base="G")])
    _annotate_codons(consensus, "ATGGCA")
    assert consensus.variants[0].is_synonymous is True


def test_variant_hits_summarise_fixed_differences(tmp_path):
    """A gene that simply differs from the closest reference is one row, not twenty."""
    mapper = ReadMapper.__new__(ReadMapper)
    consensus = LocusConsensus(
        reference="ncbi_000001", sequence="A" * 300, depth=40.0, breadth=100.0,
        variants=[_variant(position=p, allele_fraction=1.0, alt_count=40,
                           ref_aa="A", alt_aa="V", codon_position=p // 3 + 1)
                  for p in (10, 40, 70)],
    )
    hits = mapper.variant_hits("s1", "ncbi", consensus, {"gene": "tet(A)",
                                                         "accession": "NG_1"})
    assert len(hits) == 1
    assert hits[0].method == "ALLELER"
    assert "3 protein-changing position" in hits[0].note
    assert "closest reference NG_1" in hits[0].note


def test_looks_like_cds_rejects_a_non_coding_reference():
    """Promoter and rRNA references can be length-divisible by three by chance."""
    coding = "ATG" + "GCA" * 25 + "TAA"
    assert _looks_like_cds(coding) is True
    # Internal stop codons mean this frame is not a coding sequence.
    assert _looks_like_cds("ATG" + "TAA" * 25 + "TAA") is False
    assert _looks_like_cds("ACGT") is False


def test_fixed_threshold_follows_the_run_setting():
    strict = _variant(allele_fraction=0.8, fixed_threshold=0.9)
    lenient = _variant(allele_fraction=0.8, fixed_threshold=0.5)
    assert strict.is_fixed is False
    assert lenient.is_fixed is True


def test_indel_count_sees_deletions_and_insertions():
    assert _indel_count(["A", "A+2GT", "*", "T-1N", "C"]) == 3


def test_variant_hits_report_minority_alleles_individually():
    mapper = ReadMapper.__new__(ReadMapper)
    consensus = LocusConsensus(
        reference="r", sequence="A" * 300, depth=100.0, breadth=100.0,
        variants=[_variant(position=10, allele_fraction=0.25, ref_aa="S", alt_aa="L",
                           codon_position=4)],
    )
    hits = mapper.variant_hits("s1", "ncbi", consensus, {"gene": "gyrA", "accession": "X"})
    assert len(hits) == 1
    assert hits[0].method == "VARIANTR"
    assert hits[0].allele_fraction == pytest.approx(0.25)
    assert "HETERORESISTANT" in hits[0].note
    assert "S4L" in hits[0].note


# -------------------------------------------------------------------- typing
def _profile_table(rows):
    return {tuple(r[1:]): r[0] for r in rows}


def test_resolve_profile_fills_a_single_missing_locus():
    """E. faecium CC17 deletes pstS, which the scheme records as allele 0."""
    loci = ["atpA", "ddl", "gdh", "purK", "gyd", "pstS", "adk"]
    table = _profile_table([("1478", "9", "1", "1", "1", "1", "0", "1")])
    alleles = {"atpA": "9", "ddl": "1", "gdh": "1", "purK": "1", "gyd": "1",
               "pstS": "-", "adk": "1"}
    st, resolved, note = MlstTyper._resolve_profile(table, loci, alleles)
    assert st == "1478"
    assert resolved["pstS"] == "0"
    assert "pstS" in note
    # Nothing else in the table fits, so no rival STs are listed.
    assert "also fit" not in note


def test_resolve_profile_flags_competing_profiles():
    """17 H. influenzae STs differ only at fucK; say so rather than imply certainty."""
    loci = ["adk", "atpG", "frdB", "fucK", "mdh", "pgi", "recA"]
    table = _profile_table([
        ("18", "1", "1", "1", "1", "1", "1", "1"),
        ("66", "1", "1", "1", "7", "1", "1", "1"),
        ("2449", "1", "1", "1", "0", "1", "1", "1"),
    ])
    alleles = {locus: "1" for locus in loci}
    alleles["fucK"] = "-"
    st, _resolved, note = MlstTyper._resolve_profile(table, loci, alleles)
    # PubMLST's 0 convention is what the field uses, so the ST is still reported,
    # but the note has to name the STs it could equally have been.
    assert st == "2449"
    assert "also fit ST 18, 66" in note
    assert "fucK" in note


def test_scheme_profile_takes_loci_from_the_allele_files(tmp_path):
    """PubMLST can put metadata columns between the loci and clonal_complex."""
    scheme_dir = tmp_path / "pubmlst" / "aphagocytophilum"
    scheme_dir.mkdir(parents=True)
    loci = ["pheS", "glyA", "fumC", "mdh", "sucA", "dnaN", "atpA"]
    for locus in loci:
        (scheme_dir / f"{locus}.tfa").write_text(f">{locus}_1\nACGT\n")
    (scheme_dir / "aphagocytophilum.txt").write_text(
        "ST\t" + "\t".join(loci) + "\tclonal_complex\tMLST_cluster\n"
        "1\t" + "\t".join(["1"] * 7) + "\tST-1 complex\tcluster A\n")
    profiles = SchemeProfiles(tmp_path)
    found_loci, table, extra = profiles.load("aphagocytophilum")
    assert found_loci == loci
    assert table[tuple(["1"] * 7)] == "1"
    assert extra["1"] == "ST-1 complex"


def test_marker_group_separates_o_and_h_antigens():
    assert _marker_group("O88-4-wzx") == "O88"
    assert _marker_group("O88-2-wzy") == "O88"
    assert _marker_group("H7-6-fliC-origin") == "H7"
    assert _marker_group("ipaH_c") == "ipaH_c"


def test_marker_group_leaves_numbered_alleles_to_their_locus():
    """wzi alleles are numbered: they are one locus, not 500 separate markers."""
    assert _marker_group("ipaH_c") == "ipaH_c"
    assert _marker_group("O88-4-wzx") == "O88"
    # A bare number is an allele id; the caller keeps the locus as the group.
    assert _marker_group("64") == "64"
    assert "64".split(".")[0].isdigit() is True
    assert "ipaH_c".split(".")[0].isdigit() is False


def test_merge_hsps_keeps_two_copies_of_one_gene():
    """A duplicated gene is two hits; keying only on query/subject would hide one."""
    first = _hsp(qstart=20_001, qend=20_861, sstart=1, send=861, length=861, nident=861,
                 slen=861, bitscore=1500)
    second = _hsp(qstart=60_862, qend=61_722, sstart=1, send=861, length=861, nident=861,
                  slen=861, bitscore=1500)
    merged = merge_hsps([first, second])
    assert len(merged) == 2
    assert {m.coverage_pct for m in merged} == {100.0}


def _amr_hit(gene, drug_class="BETA-LACTAM", subclass="") -> Hit:
    return Hit(sample="s", database="ncbi", gene=gene, element_type="AMR",
               drug_class=drug_class, subclass=subclass)


def test_resistance_score_ignores_intrinsic_ampc():
    """Every E. coli carries blaEC; it must not make every E. coli an ESBL producer."""
    score, flags = resistance_score([_amr_hit("blaEC-15", subclass="CEPHALOSPORIN")])
    assert flags["esbl"] is False
    assert score == 0


@pytest.mark.parametrize("gene", ["blaCTX-M-15", "blaSHV-12", "blaVEB-1", "blaPER-1"])
def test_resistance_score_detects_a_real_esbl(gene):
    """An ESBL family must not be swallowed by an AmpC family with a shared prefix."""
    score, flags = resistance_score([_amr_hit(gene, subclass="CEPHALOSPORIN")])
    assert flags["esbl"] is True
    assert score == 1


@pytest.mark.parametrize("gene", ["blaEC-15", "blaACT-9", "blaCMY-2", "blaOXA-1", "blaADC-30"])
def test_resistance_score_excludes_ampc_and_narrow_spectrum(gene):
    _score, flags = resistance_score([_amr_hit(gene, subclass="CEPHALOSPORIN")])
    assert flags["esbl"] is False


def test_resistance_score_finds_carbapenemases_without_curated_subclass():
    """A database with no AMRFinderPlus family still has to score correctly."""
    score, flags = resistance_score([_amr_hit("NDM-1", subclass="")])
    assert flags["carbapenemase"] is True
    assert score == 2
    score, flags = resistance_score([_amr_hit("NDM-1"), _amr_hit("mcr-1", "COLISTIN")])
    assert score == 3


def test_virulence_score_ignores_typing_loci():
    """An ST, a serotype and a wzi allele are not virulence."""
    species = SpeciesCall(name="Escherichia coli", genus="Escherichia", species="coli")
    typing = [TypingResult(scheme="ST_achtman", call="5082"),
              TypingResult(scheme="serotype", call="O121:H7"),
              TypingResult(scheme="wzi", call="12")]
    assert virulence_score(species, typing)[0] == 0


def test_virulence_score_counts_real_loci():
    species = SpeciesCall(name="Klebsiella pneumoniae", genus="Klebsiella")
    typing = [TypingResult(scheme="AbST", call="1"), TypingResult(scheme="cbST", call="2")]
    assert virulence_score(species, typing)[0] == 5


def test_empty_results_produce_empty_tables():
    assert long_table([]).empty
    assert matrix([]).empty


def test_format_cell_preserves_small_allele_fractions():
    """Rounding an allele fraction to two decimals would erase the signal."""
    assert _format_cell(0.0043) == "0.0043"
    assert _format_cell(0.9987) == "0.9987"
    assert _format_cell(99.0) == "99"
    assert _format_cell(1.5) == "1.5"


def test_format_cell_hides_missing_values():
    assert _format_cell(float("nan")) == ""
    assert _format_cell(None) == ""
    assert _format_cell("nan") == ""


def test_class_summary_counts_distinct_genes_not_hits():
    """A gene on three replicons is one resistance determinant."""
    result = SampleResult(sample="s1")
    for start in (1, 5000, 9000):
        result.hits.append(Hit(sample="s1", database="ncbi", gene="aadA1", element_type="AMR",
                               drug_class="AMINOGLYCOSIDE", sequence=f"c{start}",
                               start=start, end=start + 500))
    summary = class_summary([result])
    assert summary.loc["s1", "AMINOGLYCOSIDE"] == 1
    counts = matrix([result], columns="class", cell="count", element_types=["AMR"])
    assert counts.loc["s1", "AMINOGLYCOSIDE"] == 3


def test_coverage_glyph_marks_covered_regions():
    full = Hit(sample="s", database="d", gene="g", coverage="1-100/100")
    half = Hit(sample="s", database="d", gene="g", coverage="1-50/100")
    assert set(_coverage_glyph(full)) == {"="}
    glyph = _coverage_glyph(half)
    assert "=" in glyph and "." in glyph
    assert _coverage_glyph(Hit(sample="s", database="d", gene="g", coverage="")) == ""


def test_gene_table_has_a_stable_column_layout():
    result = _result("s1", ["blaKPC-2"])
    result.hits[0].coverage = "1-861/861"
    frame = gene_table([result])
    assert list(frame.columns) == list(GENE_COLUMNS)
    # COVERAGE holds the full span; COVERAGE_MAP is the alignment glyph.
    assert frame.iloc[0]["COVERAGE"] == "1-861/861"
    assert set(frame.iloc[0]["COVERAGE_MAP"]) <= {"=", "."}


def test_writer_always_writes_flat_layouts_as_tsv(tmp_path):
    """--format genes must produce a file, even alongside html/json."""
    results = [_result("s1", ["blaKPC-2"])]
    written = write_outputs(results, outdir=tmp_path, formats=["genes", "json"])
    names = {p.name for p in written}
    assert "hydra.genes.tsv" in names
    assert "hydra.json" in names
    content = (tmp_path / "hydra.genes.tsv").read_text()
    assert "\t" in content.splitlines()[0]


def test_writer_honours_every_requested_table_format(tmp_path):
    results = [_result("s1", ["blaKPC-2"])]
    written = write_outputs(results, outdir=tmp_path, formats=["tsv", "csv", "elements"])
    names = {p.name for p in written}
    assert {"hydra.tsv", "hydra.csv", "hydra.elements.tsv"} <= names


# ------------------------------------------------------- upstream header formats
def test_ncbi_cds_header_conversion():
    header = ("AAA16360.1|L11078.1|1|1|stxA2b|stxA2b|Shiga_toxin_Stx2b_subunit_A "
              "L11078.1:177-1136")
    assert ncbi_header(header) == "ncbi~~~stxA2b~~~L11078.1~~~L11078.1:177-1136"


def test_ncbi_cds_header_rejects_short_records():
    assert ncbi_header("something|else") == ""


def test_card_header_conversion():
    header = "gb|GQ343019.1|+|132-1023|ARO:3002999|CblA-1 [mixed culture bacterium AX_gF3SD01_15]"
    assert card_header(header) == ("card~~~CblA-1~~~GQ343019.1~~~"
                                   "mixed culture bacterium AX_gF3SD01_15")


def test_cge_header_splits_from_the_right():
    """Gene names contain underscores, so only the last two fields are structural."""
    assert cge_header("blaNDM-19_1_MF370080", "resfinder", "beta-lactam") == \
        "resfinder~~~blaNDM-19~~~MF370080~~~beta-lactam"
    # PlasmidFinder separates the accession with a double underscore, so the copy
    # number stays on the replicon name - which is how PlasmidFinder names them.
    assert cge_header("pKPC-CAV1321_1__CP011611", "plasmidfinder", "enterobacteriales") == \
        "plasmidfinder~~~pKPC-CAV1321_1~~~CP011611~~~enterobacteriales"


def test_cge_header_keeps_a_descriptive_middle_field_out_of_the_gene_name():
    assert cge_header("repUS46_1_SAP099B017(SAP099B)_GQ900449", "plasmidfinder", "gram") == \
        "plasmidfinder~~~repUS46_1~~~GQ900449~~~gram"


def test_resfinder_copy_suffixes_normalise_away():
    """The copy number is cosmetic: both spellings must match the same locus."""
    assert _normalise_gene("blaNDM-19") == _normalise_gene("blaNDM-19_1")


def test_vfdb_header_conversion():
    header = ("VFG037176(gb|WP_001081735) (plc1) phospholipase C "
              "[Phospholipase C (VF0470) - Exotoxin (VFC0235)] [Acinetobacter baumannii ACICU]")
    assert vfdb_header(header) == "vfdb~~~plc1~~~WP_001081735~~~phospholipase C"


def test_vfdb_keeps_upstreams_accession_as_a_gene_symbol():
    """VFDB names a record after its accession when no gene symbol is known."""
    header = "VFG000710(gb|AAC38364) (AAC38364) Orf1 [Ler (VF0189) - Regulation (VFC0301)]"
    assert vfdb_header(header) == "vfdb~~~AAC38364~~~AAC38364~~~Orf1"


def test_every_database_a_full_install_needs_can_be_fetched():
    for name in ("ncbi", "card", "vfdb", "protein", "resfinder", "plasmidfinder",
                 "pubmlst", "lineage", "species"):
        assert can_fetch(name), name
    # Published only as a landing page, so it still has to be fetched by hand.
    assert not can_fetch("argannot")


# --------------------------------------------------------------------- PubMLST
def _profiles(slug: str, scheme: int = 1) -> str:
    return f"https://rest.pubmlst.org/db/pubmlst_{slug}_seqdef/schemes/{scheme}/profiles_csv"


def test_pubmlst_scheme_name_uses_the_database_slug():
    assert pubmlst_scheme_name("Klebsiella pneumoniae species complex",
                               _profiles("klebsiella")) == "klebsiella"
    assert pubmlst_scheme_name("Enterococcus faecium", _profiles("efaecium")) == "efaecium"


def test_pubmlst_scheme_name_numbers_a_databases_second_scheme():
    assert pubmlst_scheme_name("Acinetobacter baumannii#1",
                               _profiles("abaumannii")) == "abaumannii"
    assert pubmlst_scheme_name("Acinetobacter baumannii#2",
                               _profiles("abaumannii", 2)) == "abaumannii_2"


def test_pubmlst_scheme_name_keeps_the_spelling_already_in_use():
    """The derived name for these differs from what everything else calls them."""
    assert pubmlst_scheme_name("Escherichia coli#1",
                               _profiles("escherichia")) == "ecoli_achtman_4"
    assert pubmlst_scheme_name("Escherichia coli#2", _profiles("ecoli")) == "ecoli"


def test_pubmlst_scheme_name_rejects_an_unparseable_url():
    assert pubmlst_scheme_name("Whatever", "https://example.org/nothing") == ""


@pytest.mark.parametrize("label,expected", [
    ("Escherichia coli#1", ("Escherichia", "coli")),
    ("Klebsiella pneumoniae species complex", ("Klebsiella", "pneumoniae")),
    ("Achromobacter spp.", ("Achromobacter", "")),
    ("Campylobacter jejuni", ("Campylobacter", "jejuni")),
])
def test_pubmlst_organism_parsing(label, expected):
    assert pubmlst_organism(label) == expected


def test_scheme_inherits_its_taxgroup_from_its_species():
    """A scheme name the curated table has never seen must still map to an organism."""
    table = _inherit_organism_by_species({
        "ecoli_achtman_4": SchemeOrganism("Escherichia", "coli", "Escherichia"),
        "escherichia_9": SchemeOrganism("Escherichia", "coli", ""),
    })
    assert table["escherichia_9"].organism == "Escherichia"


def test_species_inheritance_does_not_guess_across_an_ambiguous_genus():
    table = _inherit_organism_by_species({
        "kpneumoniae": SchemeOrganism("Klebsiella", "pneumoniae", "Klebsiella_pneumoniae"),
        "koxytoca": SchemeOrganism("Klebsiella", "oxytoca", "Klebsiella_oxytoca"),
        "kother": SchemeOrganism("Klebsiella", "variicola", ""),
    })
    assert table["kother"].organism == ""


def test_one_unreadable_assembly_does_not_discard_the_rest_of_the_batch(tmp_path):
    """A truncated download in a thousand-genome run must not cost the other 999."""
    good = tmp_path / "good.fasta"
    good.write_text(">contig\nACGTACGTACGTACGT\n")
    bad = tmp_path / "bad.fasta"
    bad.write_text("this is not a FASTA file\n")

    results = {"good": SampleResult(sample="good"), "bad": SampleResult(sample="bad")}
    kept = Pipeline._drop_unreadable(
        Pipeline.__new__(Pipeline), {"good": good, "bad": bad}, results)

    assert set(kept) == {"good"}
    # The skipped sample is still named and carries the reason, so it cannot be
    # mistaken for a genome that was screened and found to carry nothing.
    assert results["bad"].warnings
    assert "not analysed" in results["bad"].warnings[0]
    assert not results["good"].warnings


def test_a_batch_of_entirely_unreadable_assemblies_is_still_an_error(tmp_path):
    bad = tmp_path / "bad.fasta"
    bad.write_text("not a FASTA either\n")
    results = {"bad": SampleResult(sample="bad")}

    with pytest.raises(HydraError):
        Pipeline._drop_unreadable(Pipeline.__new__(Pipeline), {"bad": bad}, results)
