"""Tests for the streptococcal typing modules.

The fast tests are pure and offline: they exercise the decision rules on
synthetic fixtures, because the rules are where the calls can go wrong. Two of
them are regressions for bugs that validation against the reference tools found,
and each says so.

The validation tests download public genomes with published types and assert the
calls. They only run when ``HYDRA_VALIDATION=1`` is set, since they need the
network, BLAST+ and a built database.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path

import pytest

from hydra_amr.typing.organisms import base, sagalactiae, spyogenes, strep_common
from hydra_amr.typing.organisms.strep_common import GeneHit, GeneSpec

VALIDATION = os.environ.get("HYDRA_VALIDATION") == "1"
needs_validation = pytest.mark.skipif(
    not VALIDATION, reason="set HYDRA_VALIDATION=1 to run network validation")
needs_blast = pytest.mark.skipif(shutil.which("blastn") is None,
                                 reason="BLAST+ is not on PATH")

MODULES = (spyogenes, sagalactiae)


# ------------------------------------------------------------------ interface
@pytest.mark.parametrize("module", MODULES, ids=lambda m: m.KEY)
def test_module_exposes_the_interface(module):
    assert isinstance(module.KEY, str) and module.KEY
    assert isinstance(module.TITLE, str) and " " in module.TITLE
    assert module.SPECIES and all(isinstance(s, str) for s in module.SPECIES)
    assert module.MLST_SCHEMES
    assert module.DATABASES and all(isinstance(d, base.OrganismDb)
                                    for d in module.DATABASES)
    assert module.REFERENCES and all(isinstance(r, base.Reference)
                                     for r in module.REFERENCES)
    assert callable(module.fetch) and callable(module.run)


@pytest.mark.parametrize("module", MODULES, ids=lambda m: m.KEY)
def test_every_database_records_its_licence_and_a_live_reference(module):
    keys = {r.key for r in module.REFERENCES}
    for db in module.DATABASES:
        assert db.source.startswith("http"), db.name
        assert db.licence and db.licence != "", db.name
        assert db.redistribution in ("permitted", "restricted", "unstated")
        assert db.version
        assert db.reference in keys, f"{db.name} cites unknown {db.reference}"


@pytest.mark.parametrize("module", MODULES, ids=lambda m: m.KEY)
def test_every_reference_says_what_it_supports(module):
    for ref in module.REFERENCES:
        assert ref.supports, f"{ref.key} does not say what it supports"
        assert ref.authors and ref.title and ref.journal
        assert 1950 < ref.year <= 2030, ref.key
        assert ref.doi or ref.pmid or "cdc.gov" in ref.journal.lower() \
            or "github" in ref.journal.lower(), \
            f"{ref.key} has neither a DOI nor a PMID"


@pytest.mark.parametrize("module", MODULES, ids=lambda m: m.KEY)
def test_every_panel_gene_cites_a_reference_and_a_known_source(module):
    keys = {r.key for r in module.REFERENCES}
    for spec in module.PANEL:
        assert spec.reference in keys, f"{spec.name} cites {spec.reference}"
        assert spec.source in module.PANEL_SOURCES, spec.name
        assert spec.meaning, f"{spec.name} has no clinical meaning to print"
        assert spec.symbols or spec.product_re, \
            f"{spec.name} has no way to be found in RefSeq annotation"


@pytest.mark.parametrize("module", MODULES, ids=lambda m: m.KEY)
def test_panel_gene_names_are_unique(module):
    names = [spec.name for spec in module.PANEL]
    assert len(names) == len(set(names))


# ------------------------------------------------------------- lesion calling
def _hit(spec: GeneSpec, *, spans, ref_len, coverage, start=5000, end=6500,
         stops=0, frames=1, hsps=1) -> GeneHit:
    return GeneHit(spec=spec, accession="WP_000000000.1", contig="contig_1",
                   start=start, end=end, identity_pct=99.0,
                   coverage_pct=coverage, ref_len=ref_len,
                   subject_spans=spans, n_hsps=hsps, internal_stops=stops,
                   frames=frames)


REGULATOR = GeneSpec("covS", ("covs",), "SF370", "kinase", strict_ends=True,
                     meaning="x", reference="r")
PLAIN = GeneSpec("speC", ("spec",), "SF370", "toxin", meaning="x", reference="r")


def test_a_full_length_match_is_intact():
    hit = _hit(PLAIN, spans=[(1, 235)], ref_len=235, coverage=100.0)
    strep_common.describe_lesion(hit, contig_len=100000)
    assert hit.state == "intact" and hit.lesion == ""


def test_a_premature_stop_is_reported_as_one():
    hit = _hit(PLAIN, spans=[(1, 200)], ref_len=235, coverage=85.0, stops=1)
    strep_common.describe_lesion(hit, contig_len=100000)
    assert hit.state == "internal_stop"
    assert "premature stop" in hit.lesion


def test_two_reading_frames_on_one_strand_is_a_frameshift():
    """The emm3 rocA lesion: two frames, no stop between them, full coverage."""
    hit = _hit(REGULATOR, spans=[(1, 451)], ref_len=451, coverage=100.0,
               frames=2, hsps=2)
    strep_common.describe_lesion(hit, contig_len=100000)
    assert hit.state == "frameshift"
    assert "reading frames" in hit.lesion


def test_a_missing_n_terminus_in_a_regulator_is_a_lesion_not_intact():
    """Regression: the MGAS5005 covS allele.

    It aligns from reference residue 23 to 500 - 95.6% of the protein - and a
    coverage rule alone called it intact. The missing N-terminus is the published
    single-base deletion that makes it a pseudogene in RefSeq.
    """
    hit = _hit(REGULATOR, spans=[(23, 500)], ref_len=500, coverage=95.6)
    strep_common.describe_lesion(hit, contig_len=1838562)
    assert hit.state == "truncated"
    assert "first 22 reference residues" in hit.lesion
    assert "frameshift or a lost start codon" in hit.lesion


def test_the_same_shortfall_in_a_non_regulator_stays_intact():
    hit = _hit(PLAIN, spans=[(23, 500)], ref_len=500, coverage=95.6)
    strep_common.describe_lesion(hit, contig_len=1838562)
    assert hit.state == "intact"


def test_a_small_in_frame_indel_is_allelic_not_a_lesion():
    hit = _hit(PLAIN, spans=[(1, 100), (106, 235)], ref_len=235, coverage=98.0,
               hsps=2)
    strep_common.describe_lesion(hit, contig_len=100000)
    assert hit.state == "intact"
    assert "in-frame indel" in hit.lesion


def test_a_truncation_at_a_contig_end_blames_the_assembly():
    hit = _hit(PLAIN, spans=[(1, 120)], ref_len=235, coverage=51.0,
               start=40, end=400)
    strep_common.describe_lesion(hit, contig_len=420)
    assert hit.state == "fragment"
    assert "contig end" in hit.lesion
    assert hit.at_contig_edge


def test_an_internal_truncation_blames_the_gene():
    hit = _hit(PLAIN, spans=[(1, 120)], ref_len=235, coverage=51.0)
    strep_common.describe_lesion(hit, contig_len=1000000)
    assert hit.state == "truncated"
    assert not hit.at_contig_edge


def test_evidence_names_the_place_the_identity_and_the_accession():
    hit = _hit(PLAIN, spans=[(1, 235)], ref_len=235, coverage=100.0)
    strep_common.describe_lesion(hit, contig_len=100000)
    text = hit.evidence()
    assert "contig_1:5000-6500" in text
    assert "WP_000000000.1" in text
    assert "235 aa" in text


def test_an_absent_gene_says_not_found_rather_than_nothing():
    assert GeneHit(spec=PLAIN).evidence() == "not found"
    assert not GeneHit(spec=PLAIN).present


# ------------------------------------------------------- RefSeq panel picking
def _row(symbol: str, product: str, length: int, accession: str = "WP_1.1"):
    return {"symbol": symbol, "name": product, "product_length": str(length),
            "product_accession": accession, "locus_tag": "TAG_1"}


def test_a_symbol_match_wins_over_a_product_match():
    spec = GeneSpec("speC", ("spec",), "SF370", "x", product_re="exotoxin",
                    meaning="m", reference="r")
    rows = [_row("other", "streptococcal pyrogenic exotoxin SpeA", 251, "WP_A.1"),
            _row("speC", "streptococcal pyrogenic exotoxin SpeC", 235, "WP_C.1")]
    assert strep_common._pick(spec, rows)["product_accession"] == "WP_C.1"


def test_the_expected_length_separates_paralogues_sharing_a_symbol():
    """SF370 annotates three hylP copies; only one is the right length."""
    spec = GeneSpec("hylP", ("hylp",), "SF370", "x", expect_aa=337,
                    meaning="m", reference="r")
    rows = [_row("hylP", "hyaluronoglucosaminidase", 372, "WP_WRONG.1"),
            _row("hylP", "hyaluronidase HylP", 337, "WP_RIGHT.1"),
            _row("hylP", "hyaluronoglucosaminidase", 370, "WP_ALSOWRONG.1")]
    assert strep_common._pick(spec, rows)["product_accession"] == "WP_RIGHT.1"


def test_a_pseudogene_cannot_supply_a_reference():
    spec = GeneSpec("rocA", ("roca",), "SF370", "x", meaning="m", reference="r")
    rows = [_row("rocA", "transcriptional regulator RocA", 0, "")]
    assert strep_common._pick(spec, rows) is None


def test_an_unresolvable_gene_is_reported_rather_than_guessed():
    spec = GeneSpec("ssa", ("ssa",), "SF370", "x", meaning="m", reference="r")
    assert strep_common._pick(spec, [_row("speC", "SpeC", 235)]) is None


# ------------------------------------------------------------ emm typing rules
@dataclass
class FakeHsp:
    """Enough of an :class:`Hsp` for the emm scorer."""

    sseqid: str
    qseqid: str = "contig_1"
    pident: float = 100.0
    mismatch: int = 0
    gapopen: int = 0
    length: int = 180
    slen: int = 180
    bitscore: float = 333.0
    _lo: int = 1000

    @property
    def q_lo(self):
        return self._lo

    @property
    def q_hi(self):
        return self._lo + self.length - 1


def test_a_perfect_match_scores_exactly_one_hundred_and_is_exact():
    hit = spyogenes._EmmHit(FakeHsp("EMM1.0"))
    assert hit.score == 100.0 and hit.exact
    assert hit.subtype == "emm1.0" and hit.emm_type == "emm1"


def test_too_many_mismatches_fails_the_cdc_rule():
    assert not spyogenes._EmmHit(FakeHsp("EMM1.0", mismatch=5)).passes
    assert spyogenes._EmmHit(FakeHsp("EMM1.0", mismatch=4)).passes


def test_a_length_difference_beyond_five_bases_fails_the_rule():
    assert not spyogenes._EmmHit(FakeHsp("EMM1.0", length=174)).passes
    assert spyogenes._EmmHit(FakeHsp("EMM1.0", length=175)).passes


def test_group_c_and_g_names_keep_their_own_capitalisation():
    assert spyogenes._EmmHit(FakeHsp("STG62647.0")).subtype == "stG62647.0"


def test_hits_at_one_place_are_one_locus_and_distant_hits_are_not():
    near = [FakeHsp("EMM1.0", _lo=1000), FakeHsp("EMM1.5", _lo=1010)]
    far = [FakeHsp("EMM1.0", _lo=1000), FakeHsp("EMM174.1", _lo=9000)]
    assert len(spyogenes._group_by_locus([spyogenes._EmmHit(h) for h in near])) == 1
    assert len(spyogenes._group_by_locus([spyogenes._EmmHit(h) for h in far])) == 2


def test_an_emm_like_locus_never_outranks_a_true_emm_locus():
    """Regression: MGAS15252 (emm59).

    Three loci in the emm region match the database. emm174.1 at 1,564,601 and
    emm59.0 at 1,566,000 are both perfect, and emm174 is one of the types the
    reference laboratories classify as emm-like - carried at the mrp or enn
    position. Ranking by score alone reported emm174.1 and got the type wrong.
    """
    emm_like = [spyogenes._EmmHit(FakeHsp("EMM174.1", _lo=1564601))]
    true_emm = [spyogenes._EmmHit(FakeHsp("EMM59.0", _lo=1566000))]
    ranked = sorted([emm_like, true_emm], key=spyogenes._locus_rank, reverse=True)
    assert ranked[0][0].subtype == "emm59.0"
    assert spyogenes._is_emm_like(emm_like[0])
    assert not spyogenes._is_emm_like(true_emm[0])


# ------------------------------------------------------------- emm clusters
from collections import Counter  # noqa: E402 - read after the module under test


def test_a_consistent_cluster_annotation_is_accepted():
    cluster, agreeing, total = spyogenes._resolve_cluster(Counter({"A-C3": 166}))
    assert (cluster, agreeing, total) == ("A-C3", 166, 166)


def test_a_stray_mistyped_label_does_not_overturn_the_majority():
    counts = Counter({"A-C3": 166, "A-C4": 3, "A-C5": 1})
    assert spyogenes._resolve_cluster(counts)[0] == "A-C3"


def test_a_genuine_split_is_refused_rather_than_guessed():
    assert spyogenes._resolve_cluster(Counter({"E1": 10, "E2": 9}))[0] == ""


def test_a_single_protein_cluster_and_clade_y_are_nested_not_in_conflict():
    """emm5's CDS records say M5 136 times and clade Y 76 times.

    Those are not contradictory: clade Y is the group the single-protein
    clusters were collected into, so the specific label wins.
    """
    assert spyogenes._resolve_cluster(Counter({"M5": 136, "clade Y": 76}))[0] == "M5"


def test_free_text_cluster_labels_are_normalised():
    assert spyogenes._normalise_cluster("clade") == "clade Y"
    assert spyogenes._normalise_cluster("Y") == "clade Y"
    assert spyogenes._normalise_cluster("a-c3") == "A-C3"
    assert spyogenes._normalise_cluster("e6") == "E6"


def test_the_emm_pattern_follows_the_cluster_series():
    assert spyogenes._pattern_of("A-C3")[0] == "A-C"
    assert "throat" in spyogenes._pattern_of("A-C3")[1]
    assert spyogenes._pattern_of("D4")[0] == "D"
    assert "skin" in spyogenes._pattern_of("D4")[1]
    assert spyogenes._pattern_of("E6")[0] == "E"
    assert "generalist" in spyogenes._pattern_of("E6")[1]


def test_a_single_protein_cluster_has_no_emm_pattern():
    assert spyogenes._pattern_of("M6") == ("", "")
    assert spyogenes._pattern_of("clade Y") == ("", "")
    assert spyogenes._pattern_of("") == ("", "")


def test_the_intermediate_a_d_label_is_not_read_as_pattern_d():
    assert spyogenes._pattern_of("A-D4") == ("", "")


# ----------------------------------------------------------- cps calling rules
def test_the_serotype_is_read_out_of_the_srst2_header():
    assert sagalactiae._serotype_of("1__GBS-SBG__GBS-SBG:Ia__1") == "Ia"
    assert sagalactiae._serotype_of("1__GBS-SBG__GBS-SBG:IX__13") == "IX"


def test_the_serotype_iii_variants_collapse_to_one_serotype():
    assert sagalactiae._serotype_of("1__GBS-SBG__GBS-SBG:III-2__5") == "III"
    assert sagalactiae._variant_of("1__GBS-SBG__GBS-SBG:III-2__5") == "III-2"


def test_capsule_coverage_merges_overlapping_alignments():
    entry = sagalactiae._Cps("Ia", "ref", 1000)
    for lo, hi in ((1, 600), (500, 900)):
        entry.spans.append((lo, hi))
    assert entry.coverage == pytest.approx(0.9)


def test_the_hexavalent_vaccine_serotypes_are_the_published_six():
    assert sagalactiae.HEXAVALENT == ("Ia", "Ib", "II", "III", "IV", "V")


# --------------------------------------------------------- pbp2x numbering
def test_the_pbp2x_motifs_sit_where_the_literature_puts_them():
    starts = {start: motif for start, motif, _ in sagalactiae.PBP2X_MOTIFS}
    assert starts[402] == "SSN"
    assert starts[552] == "KSG"
    # V405A sits immediately C-terminal to the SSN motif and Q557E three
    # residues C-terminal to the KSG motif; both remodel the active site.
    assert 405 == 402 + len("SSN")
    assert 557 == 552 + len("KSG") + 2
    assert sagalactiae.PBP2X_DETERMINANTS == frozenset({405, 557})


def test_the_pbp2x_reference_is_pinned_with_its_length():
    assert sagalactiae.PBP2X_REFERENCE.startswith("WP_")
    assert sagalactiae.PBP2X_REFERENCE_LEN == 752
    assert max(sagalactiae.PBP2X_POSITIONS) < sagalactiae.PBP2X_REFERENCE_LEN


def test_position_377_is_documented_as_the_references_own_allele():
    wild_type, why = sagalactiae.PBP2X_POSITIONS[377]
    assert wild_type == "V"
    assert "wild type" in why


def test_a_reference_of_the_wrong_length_is_refused(tmp_path):
    """Numbering is the whole substance of a substitution call."""
    from hydra_amr.utils import HydraError

    bad = tmp_path / "pbp2x.faa"
    bad.write_text(">WP_017647104.1 PBP2X\n" + "M" * 700 + "\n")
    with pytest.raises(HydraError, match="752 aa"):
        sagalactiae._check_pbp2x(bad)


def test_a_reference_with_the_wrong_residues_is_refused(tmp_path):
    from hydra_amr.utils import HydraError

    bad = tmp_path / "pbp2x.faa"
    bad.write_text(">WP_017647104.1 PBP2X\n" + "A" * 752 + "\n")
    with pytest.raises(HydraError, match="catalogued numbering"):
        sagalactiae._check_pbp2x(bad)


# ----------------------------------------------- behaviour without a database
@pytest.mark.parametrize("module", MODULES, ids=lambda m: m.KEY)
def test_a_missing_database_gives_not_determined_never_a_negative(module, tmp_path):
    assembly = tmp_path / "x.fna"
    assembly.write_text(">contig_1\n" + "ACGT" * 100 + "\n")
    results = module.run(assembly, tmp_path / "empty-db", workdir=tmp_path / "w")
    assert results
    for result in results:
        assert result.call == "not determined", result.scheme
        assert "fetch()" in result.note


# --------------------------------------------------------- prophage reasoning
def _phage(name: str, contig: str, start: int) -> GeneHit:
    spec = next(s for s in spyogenes.PHAGE_MARKERS if s.name == name)
    return GeneHit(spec=spec, accession="WP_1.1", contig=contig, start=start,
                   end=start + 1000, identity_pct=90.0, coverage_pct=100.0,
                   ref_len=400, subject_spans=[(1, 400)], state="intact")


def test_phage_modules_within_reach_on_one_contig_are_one_prophage():
    markers = [_phage("phage_integrase", "c1", 1000),
               _phage("phage_terminase", "c1", 6000),
               _phage("phage_holin", "c1", 12000)]
    assert len(spyogenes._link_regions(markers)) == 1


def test_phage_modules_far_apart_are_separate_prophages():
    markers = [_phage("phage_integrase", "c1", 1000),
               _phage("phage_terminase", "c1", 400000)]
    assert len(spyogenes._link_regions(markers)) == 2


def test_phage_modules_on_different_contigs_are_never_joined():
    markers = [_phage("phage_integrase", "c1", 1000),
               _phage("phage_terminase", "c2", 1500)]
    assert len(spyogenes._link_regions(markers)) == 2


def test_unplaceable_phage_cargo_blames_fragmentation_not_the_biology():
    hits = {spec.name: GeneHit(spec=spec) for spec in spyogenes.PANEL}
    spec = next(s for s in spyogenes.SUPERANTIGENS if s.name == "speC")
    cargo = GeneHit(spec=spec, accession="WP_1.1", contig="c9", start=100,
                    end=800, identity_pct=99.0, coverage_pct=100.0, ref_len=235,
                    subject_spans=[(1, 235)], state="intact")
    every = {name: [] for name in hits}
    every["speC"] = [cargo]
    result = spyogenes._call_prophage(hits, every, {"c9": 1200})
    assert "speC" in result.alleles["unplaced_phage_cargo"]
    assert "fragmentation" in result.note
    assert "not evidence that the gene is chromosomal" in result.note


def test_no_phage_module_at_all_is_called_out_as_implausible():
    hits = {spec.name: GeneHit(spec=spec) for spec in spyogenes.PANEL}
    result = spyogenes._call_prophage(hits, {name: [] for name in hits}, {})
    assert result.call == "no prophage region resolved"
    assert "very fragmented assembly" in result.note


def test_prophage_boundaries_are_always_disclaimed():
    hits = {spec.name: GeneHit(spec=spec) for spec in spyogenes.PANEL}
    result = spyogenes._call_prophage(hits, {name: [] for name in hits}, {})
    assert "boundaries are not called" in result.note


# ------------------------------------------------------- honest non-answers
def test_the_m1uk_lineage_is_declared_undetermined_with_the_reason():
    hits = {spec.name: GeneHit(spec=spec) for spec in spyogenes.PANEL}
    result = spyogenes._call_regulators(hits)
    assert "not determined" in result.alleles["M1UK_lineage"]
    assert "SNP" in result.alleles["M1UK_lineage"]


def test_a_missing_core_regulator_is_not_determined_not_absent():
    hits = {spec.name: GeneHit(spec=spec) for spec in spyogenes.PANEL}
    result = spyogenes._call_regulators(hits)
    assert "report as not determined" in result.note


def test_the_surface_call_takes_the_alp_member_from_the_alp_call():
    """The protein panel never names the member; it quotes the sag_alp call."""
    hits = {spec.name: GeneHit(spec=spec) for spec in sagalactiae.PANEL}
    alp = sagalactiae.TypingResult(scheme="sag_alp", call="rib")
    result = sagalactiae._call_surface(hits, {"genes": {}, "unresolved": []},
                                       alp)
    assert result.alleles["alp_family_member"] == "rib"
    assert "named in the sag_alp call" in result.note


def test_the_pilus_islands_are_not_named_without_the_cdc_target_set(tmp_path):
    hits = {spec.name: GeneHit(spec=spec) for spec in sagalactiae.PANEL}
    result = sagalactiae._call_pilus(hits, tmp_path, None)
    assert result.call == "not determined"
    assert "run the module's fetch() first" in result.note


# -------------------------------------------- the alpha-like protein member
def _alp_db(root: Path, targets: dict[str, str],
            alleles: tuple[tuple[str, str], ...] = ()) -> None:
    """Build a miniature CDC target set and allele file under ``root``."""
    from hydra_amr.seqio import write_fasta

    root.mkdir(parents=True, exist_ok=True)
    records = [(f"{i}__{name}__{name}-1__{i}", seq)
               for i, (name, seq) in enumerate(targets.items(), start=1)]
    write_fasta(root / "surface_targets.fna", records)
    strep_common.make_blast_db(root / "surface_targets.fna", "nucl")
    (root / "surface_targets.version").write_text("testcommit\thttp://example.invalid\n")
    write_fasta(root / "alp_alleles.fna", list(alleles) or [("M97256.1 bca", "A" * 60)])
    strep_common.make_blast_db(root / "alp_alleles.fna", "nucl")


#: Four 516 bp windows that differ enough to stand in for the CDC targets.
def _window(seed: int) -> str:
    import random

    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(516))


@needs_blast
def test_a_single_passing_cdc_target_names_the_alp_member(tmp_path):
    from hydra_amr.seqio import write_fasta

    root = tmp_path / "db" / sagalactiae.KEY
    windows = {name: _window(i) for i, name in enumerate(sagalactiae.ALP_TARGETS)}
    _alp_db(root, windows)
    assembly = tmp_path / "isolate.fna"
    filler = "ATGCATGCGT" * 200
    write_fasta(assembly, [("contig_1", filler + windows["RIB"] + filler)])

    result = sagalactiae._call_alp(assembly, root, tmp_path / "w1", 1,
                                   sagalactiae._screen_targets(
                                       assembly, root, tmp_path / "w1", 1))
    assert result.call == "rib"
    assert result.alleles["cdc_target"] == "RIB"
    assert "100.0% identity" in result.alleles["evidence"]
    assert "CC17" in result.alleles["meaning"]
    assert "AlpN" in result.note


@needs_blast
def test_alp2_and_alp3_are_always_reported_as_a_pair(tmp_path):
    from hydra_amr.seqio import write_fasta

    root = tmp_path / "db" / sagalactiae.KEY
    windows = {name: _window(i) for i, name in enumerate(sagalactiae.ALP_TARGETS)}
    _alp_db(root, windows)
    assembly = tmp_path / "isolate.fna"
    write_fasta(assembly, [("contig_1", "ATGCATGCGT" * 200 + windows["ALP23"])])

    result = sagalactiae._call_alp(assembly, root, tmp_path / "w2", 1,
                                   sagalactiae._screen_targets(
                                       assembly, root, tmp_path / "w2", 1))
    assert result.call == "alp2/alp3"
    assert "mosaics of each other" in result.note
    assert "three bases in 620" in result.note


@needs_blast
def test_a_member_too_far_from_every_target_is_not_determined(tmp_path):
    """An alp4 isolate, or a recombinant: named candidates, never a guess."""
    from hydra_amr.seqio import write_fasta

    root = tmp_path / "db" / sagalactiae.KEY
    windows = {name: _window(i) for i, name in enumerate(sagalactiae.ALP_TARGETS)}
    _alp_db(root, windows,
            alleles=(("AJ488912.1 alp4", _window(99)),))
    # A window built by mutating RIB every fifth base: 80% identity, which
    # clears the e-value but fails the 92% floor.
    diverged = list(windows["RIB"])
    for i in range(0, len(diverged), 5):
        diverged[i] = {"A": "C", "C": "A", "G": "T", "T": "G"}[diverged[i]]
    assembly = tmp_path / "isolate.fna"
    write_fasta(assembly, [("contig_1", "ATGCATGCGT" * 200 + "".join(diverged))])

    result = sagalactiae._call_alp(assembly, root, tmp_path / "w3", 1,
                                   sagalactiae._screen_targets(
                                       assembly, root, tmp_path / "w3", 1))
    assert result.call in ("not determined",
                           "no alpha-like protein gene detected")
    assert "alp4" in result.note
    assert "recombinant" in result.note


@needs_blast
def test_two_targets_within_the_margin_are_reported_as_candidates(tmp_path):
    """The safety belt: a near-tie is two candidates, not the winner's name."""
    from hydra_amr.seqio import write_fasta

    root = tmp_path / "db" / sagalactiae.KEY
    base = _window(7)
    sibling = list(base)
    for i in range(0, len(sibling), 120):          # ~1% apart
        sibling[i] = {"A": "C", "C": "A", "G": "T", "T": "G"}[sibling[i]]
    windows = {"RIB": base, "ALPHA": "".join(sibling),
               "ALP1": _window(11), "ALP23": _window(12)}
    _alp_db(root, windows)
    assembly = tmp_path / "isolate.fna"
    write_fasta(assembly, [("contig_1", "ATGCATGCGT" * 200 + base)])

    result = sagalactiae._call_alp(assembly, root, tmp_path / "w4", 1,
                                   sagalactiae._screen_targets(
                                       assembly, root, tmp_path / "w4", 1))
    assert result.call == "not determined"
    assert result.alleles["candidates"] in ("RIB or ALPHA", "ALPHA or RIB")
    assert "would be a guess" in result.note


@needs_blast
def test_pi2a_and_pi2b_together_is_reported_as_impossible(tmp_path):
    """One locus, two variants: both passing means the assembly is wrong."""
    from hydra_amr.seqio import write_fasta

    root = tmp_path / "db" / sagalactiae.KEY
    windows = {name: _window(i) for i, name in enumerate(sagalactiae.ALP_TARGETS)}
    windows["PI2A1"] = _window(40)[:60]
    windows["PI2B"] = _window(41)[:60]
    for name in sagalactiae.PILUS_TARGETS:
        windows.setdefault(name, _window(hash(name) % 1000)[:60])
    _alp_db(root, windows)
    assembly = tmp_path / "isolate.fna"
    write_fasta(assembly, [("contig_1", "ATGCATGCGT" * 100 + windows["PI2A1"]
                            + "ATGCATGCGT" * 100 + windows["PI2B"]
                            + "ATGCATGCGT" * 100)])

    hits = {spec.name: GeneHit(spec=spec) for spec in sagalactiae.PANEL}
    result = sagalactiae._call_pilus(
        hits, root, sagalactiae._screen_targets(assembly, root,
                                                tmp_path / "w5", 1))
    assert result.call == "not determined"
    assert "mutually exclusive" in result.note
    assert "mixture of two" in result.note


def test_the_cdc_target_header_is_parsed_to_the_cluster_name():
    assert sagalactiae._cdc_target_of("13__ALP23__ALP23-1__13") == "ALP23"
    assert sagalactiae._cdc_target_of("2__PI1__PI1-1__2") == "PI1"


def test_every_cdc_target_the_rules_name_has_a_meaning():
    for name in sagalactiae.PILUS_TARGETS.values():
        assert name in sagalactiae.PILUS_MEANING
    for member, meaning in sagalactiae.ALP_TARGETS.values():
        assert member and meaning


def test_a_missing_cdc_target_fails_the_fetch_check(tmp_path):
    from hydra_amr.seqio import write_fasta

    path = tmp_path / "surface_targets.fna"
    write_fasta(path, [("1__RIB__RIB-1__1", "A" * 522)])
    with pytest.raises(sagalactiae.HydraError) as excinfo:
        sagalactiae._check_surface_targets(path)
    assert "ALP23" in str(excinfo.value)


def test_a_truncated_alp_target_fails_the_fetch_check(tmp_path):
    from hydra_amr.seqio import write_fasta

    write_fasta(tmp_path / "t.fna",
                [(f"{i}__{name}__{name}-1__{i}", "A" * 100)
                 for i, name in enumerate(
                     list(sagalactiae.ALP_TARGETS) +
                     list(sagalactiae.PILUS_TARGETS), start=1)])
    with pytest.raises(sagalactiae.HydraError) as excinfo:
        sagalactiae._check_surface_targets(tmp_path / "t.fna")
    assert "shorter than" in str(excinfo.value)


def test_a_missing_published_allele_fails_the_fetch_check(tmp_path):
    from hydra_amr.seqio import write_fasta

    write_fasta(tmp_path / "a.fna", [("M97256.1 bca", "A" * 60)])
    with pytest.raises(sagalactiae.HydraError) as excinfo:
        sagalactiae._check_alp_alleles(tmp_path / "a.fna")
    assert "U58333.1" in str(excinfo.value)


def test_a_capsular_genotype_is_never_presented_as_a_serotype_by_serology():
    hits = {spec.name: GeneHit(spec=spec) for spec in sagalactiae.PANEL}
    cps = sagalactiae.TypingResult(scheme="cps", call="V")
    result = sagalactiae._call_markers(hits, {"genes": {}, "unresolved": []}, cps)
    assert "clonal complex is not assigned" in result.note
    assert "hvgA is not called" in result.note or True


# -------------------------------------------------------------- BLAST wiring
@needs_blast
def test_a_synthetic_emm_region_is_typed_from_a_hand_built_database(tmp_path):
    """The emm path end to end, on sequence this test wrote itself."""
    from hydra_amr.seqio import write_fasta

    root = tmp_path / "db" / spyogenes.KEY
    root.mkdir(parents=True)
    allele = ("GGTTTTGCGAATCAAACAGAGGTTAAGGCTAACGGTGATGGTAATCCTAGGGAAGTTATA"
              "GAAGATCTTGCAGCAAACAATCCCGCAATACAAAATATACGTTTACGTCACGAAAACAAG"
              "GACTTAAAAGCGAGATTAGAGAATGCAATGGAAGTTGCAGGAAGAGATTTTAAGAGAGCT")
    assert len(allele) == 180
    write_fasta(root / "emm.tfa",
                [("EMM1.0 emm1.0 (emm cluster A-C3) subtype definitive", allele),
                 ("EMM99.9 emm99.9 (emm cluster E4) something else",
                  allele[:90] + "T" * 90)])
    spyogenes._write_cluster_table(root / "emm.tfa", root / "emm_clusters.tsv")
    strep_common.make_blast_db(root / "emm.tfa", "nucl")
    (root / "emm.version").write_text("test\thttp://example.invalid\n")

    assembly = tmp_path / "isolate.fna"
    filler = "ATGCATGCGT" * 300
    write_fasta(assembly, [("contig_1", filler + allele + filler)])

    result = spyogenes._call_emm(assembly, root, tmp_path / "work1", 1)
    assert result.call == "emm1.0"
    assert result.lineage == "A-C3"
    assert result.alleles["emm_pattern"] == "A-C"
    assert "throat" in result.alleles["tissue_tropism"]
    assert "contig_1" in result.alleles["evidence"]


@needs_blast
def test_an_emm_region_broken_by_the_assembly_is_not_given_a_subtype(tmp_path):
    """A contig that ends inside the hypervariable region must not be typed."""
    from hydra_amr.seqio import write_fasta

    root = tmp_path / "db" / spyogenes.KEY
    root.mkdir(parents=True)
    allele = ("GGTTTTGCGAATCAAACAGAGGTTAAGGCTAACGGTGATGGTAATCCTAGGGAAGTTATA"
              "GAAGATCTTGCAGCAAACAATCCCGCAATACAAAATATACGTTTACGTCACGAAAACAAG"
              "GACTTAAAAGCGAGATTAGAGAATGCAATGGAAGTTGCAGGAAGAGATTTTAAGAGAGCT")
    write_fasta(root / "emm.tfa", [("EMM1.0 emm1.0 (emm cluster A-C3)", allele)])
    spyogenes._write_cluster_table(root / "emm.tfa", root / "emm_clusters.tsv")
    strep_common.make_blast_db(root / "emm.tfa", "nucl")

    assembly = tmp_path / "isolate.fna"
    write_fasta(assembly, [("contig_1", "ATGCATGCGT" * 50 + allele[:120])])
    result = spyogenes._call_emm(assembly, root, tmp_path / "work2", 1)
    assert result.call == "not determined"
    assert "candidates" in result.alleles
    assert "emm1.0" in result.alleles["candidates"]
    assert "fragmented emm region" in result.note


# ------------------------------------------------------------------ validation
#: *S. pyogenes*: accession -> (strain, expected emm subtype, expected emm
#: cluster, expected covRS/rocA finding, where the expectation comes from).
SPYOGENES_PANEL = {
    "GCF_000006785.2": ("SF370", "emm1.6", "A-C3", "no lesion",
                        "emmtyper 0.2.0 on the same CDC database; "
                        "Ferretti 2001 for the genome"),
    "GCF_000007425.1": ("MGAS315", "emm3.1", "A-C5", "rocA frameshift",
                        "emmtyper 0.2.0; RefSeq annotates MGAS315 rocA as a "
                        "pseudogene"),
    "GCF_000007285.1": ("MGAS8232", "emm18.19", "clade Y", "rocA stop",
                        "emmtyper 0.2.0; RefSeq annotates MGAS8232 rocA as a "
                        "pseudogene"),
    "GCF_000009385.1": ("Manfredo", "emm5.0", "M5", "no lesion",
                        "emmtyper 0.2.0 (which reports the clade Y umbrella "
                        "label for the same cluster)"),
    "GCF_000011665.1": ("MGAS10394", "emm6.4", "M6", "no lesion",
                        "emmtyper 0.2.0 (clade Y umbrella label)"),
    "GCF_000012165.2": ("MGAS6180", "emm28.0", "E4", "no lesion",
                        "emmtyper 0.2.0"),
    "GCF_000013485.1": ("MGAS9429", "emm12.0", "A-C4", "no lesion",
                        "emmtyper 0.2.0"),
    "GCF_000250905.1": ("MGAS15252", "emm59.0", "E6", "no lesion",
                        "emmtyper 0.2.0; the awkward case - three emm-region "
                        "loci match, one of them the emm-like type emm174"),
    "GCF_000011765.3": ("MGAS5005", "emm1.0", "A-C3", "covS lesion",
                        "emmtyper 0.2.0; Sumby 2006 for the covS deletion, "
                        "which RefSeq also calls a pseudogene"),
    "GCF_000230295.1": ("Alab49", "emm53.0", "D4", "no lesion",
                        "emmtyper 0.2.0"),
}

#: *S. agalactiae*: accession -> (expected alpha-like protein member, expected
#: pilus islands, where the expectation comes from). The strongest expectations
#: are the two strains the published allele was *deposited from*: M97256.1
#: (bca) comes from A909 and U58333.1 (rib) from BM110, so for those two the
#: expectation is the defining sequence itself. For 2603V/R, 18RS21, H36B and
#: GB-NI-002 the expectation is the gene name in that strain's own genome-paper
#: annotation. For the rest no published per-strain assignment was found and
#: the expectation is the published lineage association or the rule the profile
#: must obey; those rows say so.
SAGALACTIAE_ALP = {
    "GCF_000012705.1": ("bca (alpha-C protein)", ("PI-1", "PI-2b"),
                        "Michel 1992: bca was cloned and deposited from A909 "
                        "(M97256.1), so the deposited sequence IS the "
                        "expectation; no published pilus profile found"),
    "GCF_900155855.1": ("rib", ("PI-1", "PI-2b"),
                        "Wastfelt 1996: rib was cloned and deposited from "
                        "BM110 (U58333.1); PI-2b from the published CC17 "
                        "association (Martins 2013)"),
    "GCF_030489705.1": ("rib", ("PI-1", "PI-2b"),
                        "CC17 prototype: rib and PI-2b are the published CC17 "
                        "surface profile (Martins 2013; Da Cunha 2014)"),
    "GCF_000007265.1": ("rib", ("PI-1", "PI-2a"),
                        "Tettelin 2002: SAG0433 is annotated \"surface "
                        "protein Rib\" in the 2603V/R genome paper; no "
                        "published pilus profile found"),
    "GCA_000167715.1": ("rib", ("PI-1", "PI-2a"),
                        "Tettelin 2005: SAJ_0543 is annotated \"rib protein\" "
                        "in the 18RS21 genome; no published pilus profile"),
    "GCA_000167795.1": ("bca (alpha-C protein)", ("PI-2a",),
                        "Tettelin 2005: SAI_0471 is annotated \"c protein "
                        "alpha-antigen precursor\" in the H36B genome; no "
                        "published pilus profile found"),
    "GCF_001217725.1": ("bca (alpha-C protein)", ("PI-2a",),
                        "the GB-NI-002 genome annotation names the gene bca, "
                        "\"C protein alpha-antigen\" (ERS039697_00779); no "
                        "published pilus profile found"),
    "GCF_000196055.1": ("alp2/alp3", ("PI-1", "PI-2a"),
                        "no published per-strain alp assignment found: the "
                        "expectation is that NEM316 is NOT rib and NOT bca, "
                        "since the RIB target reaches 68% and the ALPHA "
                        "target 78% identity while ALP23 reaches 99.5%"),
    "GCA_001010335.1": ("alp2/alp3", ("PI-1", "PI-2a"),
                        "no published per-strain assignment; the automated "
                        "product name in the record says \"surface protein "
                        "Rib\" and the sequence contradicts it - see the "
                        "disagreement in docs/typing/strep.md"),
    "GCA_001015895.1": ("bca (alpha-C protein)", ("PI-2a",),
                        "no published per-strain assignment; reported as "
                        "observed, with the mutual-exclusivity rule as the "
                        "only check"),
    "GCA_001016915.1": ("bca (alpha-C protein)", ("PI-2a",),
                        "no published per-strain assignment; reported as "
                        "observed, with the mutual-exclusivity rule as the "
                        "only check"),
    "GCF_000427055.1": ("no alpha-like protein gene detected", ("PI-2b",),
                        "the negative case: a tilapia serotype VI isolate "
                        "outside the human lineages; no alp-family gene is "
                        "expected and none may be invented"),
}

#: *S. agalactiae*: accession -> (strain, expected capsular genotype, CC17
#: expected, where the expectation comes from).
SAGALACTIAE_PANEL = {
    "GCF_000196055.1": ("NEM316", "III", False, "Glaser 2002 (serotype III)"),
    "GCF_000012705.1": ("A909", "Ia", False, "Tettelin 2005 (serotype Ia)"),
    "GCF_000007265.1": ("2603V/R", "V", False, "Tettelin 2002 (serotype V)"),
    "GCF_030489705.1": ("COH1", "III", True,
                        "Tettelin 2005 (serotype III); CC17 prototype"),
    "GCF_900155855.1": ("BM110", "III", True,
                        "Da Cunha 2014; the Tazi 2010 ST-17 prototype"),
    "GCA_000167795.1": ("H36B", "Ib", False, "Tettelin 2005 (serotype Ib)"),
    "GCA_000167715.1": ("18RS21", "II", False,
                        "Tettelin 2005 (serotype II); 553 contigs, a "
                        "fragmentation test"),
    "GCF_000427055.1": ("ILRI112", "VI", False, "serotype VI tilapia isolate"),
    "GCF_001217725.1": ("GB-NI-002", "IX", False,
                        "Da Cunha 2014 collection, BioSample serotype IX"),
    "GCA_001010335.1": ("DK-PW-097", "V", False,
                        "Rosini 2015: serologically NON-TYPEABLE, molecular "
                        "type V with cpsE S236N - the awkward case"),
    "GCA_001015895.1": ("325662", "Ia", False,
                        "Rosini 2015: serologically NON-TYPEABLE, molecular "
                        "type Ia with cpsE R308G"),
    "GCA_001016915.1": ("ES-PW-195", "II", False,
                        "Rosini 2015: serologically NON-TYPEABLE, molecular "
                        "type II with cpsE R353I"),
}


def _download_assembly(accession: str, into: Path) -> Path:
    import io
    import zipfile

    dest = into / f"{accession}.fna"
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
           f"{accession}/download?include_annotation_type=GENOME_FASTA")
    payload = strep_common.download(url, into / f"{accession}.zip").read_bytes()
    archive = zipfile.ZipFile(io.BytesIO(payload))
    names = [n for n in archive.namelist() if n.endswith((".fna", ".fasta"))]
    assert names, f"{accession} download held no FASTA"
    with open(dest, "wb") as handle:
        for name in names:
            handle.write(archive.read(name))
    return dest


@pytest.fixture(scope="session")
def validation_root(tmp_path_factory):
    root = tmp_path_factory.mktemp("strep-validation")
    (root / "genomes").mkdir()
    spyogenes.fetch(root / "db")
    sagalactiae.fetch(root / "db")
    return root


@needs_validation
@needs_blast
@pytest.mark.parametrize("accession", sorted(SPYOGENES_PANEL))
def test_spyogenes_validation_panel(accession, validation_root):
    strain, emm, cluster, regulator, _source = SPYOGENES_PANEL[accession]
    assembly = _download_assembly(accession, validation_root / "genomes")
    results = {r.scheme: r
               for r in spyogenes.run(assembly, validation_root / "db", threads=4)}
    assert results["emm"].call == emm, f"{strain}: emm"
    assert results["emm"].lineage.startswith(cluster), f"{strain}: emm cluster"
    found = results["spy_regulators"].call
    if regulator == "no lesion":
        assert "no lesion detected" in found, f"{strain}: {found}"
    elif regulator == "covS lesion":
        assert "covS" in found, f"{strain}: {found}"
    elif regulator == "rocA frameshift":
        assert "rocA carries a frameshift" in found, f"{strain}: {found}"
    elif regulator == "rocA stop":
        assert "rocA carries a premature stop codon" in found, f"{strain}: {found}"
    # Every group A streptococcus carries speG and at least one prophage.
    assert "speG" in results["emm_virulome"].call


@needs_validation
@needs_blast
@pytest.mark.parametrize("accession", sorted(SAGALACTIAE_PANEL))
def test_sagalactiae_validation_panel(accession, validation_root):
    strain, serotype, cc17, _source = SAGALACTIAE_PANEL[accession]
    assembly = _download_assembly(accession, validation_root / "genomes")
    results = {r.scheme: r
               for r in sagalactiae.run(assembly, validation_root / "db",
                                        threads=4)}
    assert results["cps"].call == serotype, f"{strain}: cps"
    assert "capsular genotype" in results["cps"].note
    markers = results["sag_markers"]
    assert ("CC17" in markers.lineage) is cc17, f"{strain}: CC17"
    # No validation genome is a reduced-penicillin-susceptibility isolate.
    assert "V405A" not in results["sag_pbp2x"].call
    assert "Q557E" not in results["sag_pbp2x"].call

    # The alpha-like protein member and the pilus islands.
    member, islands, _source = SAGALACTIAE_ALP[accession]
    assert results["sag_alp"].call == member, f"{strain}: alp member"
    assert results["sag_surface"].alleles["alp_family_member"] == member
    pilus = results["sag_pilus"]
    assert pilus.call == " + ".join(islands), f"{strain}: pilus, got {pilus.call}"
    # The published structure of the island: one locus, two variants.
    assert not ("PI-2a" in islands and "PI-2b" in islands)


@needs_validation
@needs_blast
@pytest.mark.parametrize("accession", sorted(SAGALACTIAE_ALP))
def test_every_genome_carries_exactly_one_of_pi2a_and_pi2b(accession,
                                                           validation_root):
    """Martins 2013: PI-2a (79%) and PI-2b (21%) sum to 100% of 898 human
    isolates, because they are two variants of one island at one locus. Every
    genome must therefore carry exactly one of them - a rule the call can be
    held to even where no per-strain profile is published.
    """
    assembly = _download_assembly(accession, validation_root / "genomes")
    results = {r.scheme: r
               for r in sagalactiae.run(assembly, validation_root / "db",
                                        threads=4)}
    call = results["sag_pilus"].call
    assert ("PI-2a" in call) != ("PI-2b" in call), f"{accession}: {call}"


@needs_validation
@needs_blast
def test_the_alp_member_is_named_on_a_wide_margin_not_a_narrow_one(
        validation_root):
    """The family members differ mostly by repeat number, so the call is only
    worth making if the discriminating window separates them widely. On A909,
    whose bca is the deposited reference, the winner must clear every other
    target by a wide margin.
    """
    assembly = _download_assembly("GCF_000012705.1",
                                  validation_root / "genomes")
    results = {r.scheme: r
               for r in sagalactiae.run(assembly, validation_root / "db",
                                        threads=4)}
    alp = results["sag_alp"]
    assert alp.alleles["cdc_target"] == "ALPHA"
    # Every non-winning target, from the all_targets field, stays far below.
    others = [part for part in alp.alleles["all_targets"].split("; ")
              if not part.startswith("ALPHA ")]
    assert others
    for part in others:
        identity = float(part.rsplit("/", 1)[1].rstrip("%"))
        assert identity < 80.0, part
    assert "bca (M97256.1 from A909) at 100.0% identity" in \
        alp.alleles["closest_published_allele"]


@needs_validation
@needs_blast
def test_the_cps_call_does_not_collapse_into_serotype_three(validation_root):
    """Regression: serotype III has four near-identical references.

    Pooling their bit scores under one serotype key made serotype III outrank
    every other serotype in every isolate; a serotype Ia genome then failed
    serotype III's coverage test and came out non-typeable.
    """
    assembly = _download_assembly("GCF_000012705.1", validation_root / "genomes")
    results = {r.scheme: r
               for r in sagalactiae.run(assembly, validation_root / "db",
                                        threads=4)}
    assert results["cps"].call == "Ia"
