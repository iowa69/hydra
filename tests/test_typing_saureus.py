"""Tests for the Staphylococcus aureus typing module.

The fast tests use small synthetic fixtures and need neither network nor
BLAST+. They cover the decision rules, which is where this module can be wrong
in a way that matters: the spa repeat tiling, the Ridom ``NT,`` sentinel, the
lesion descriptions for the bla regulators, the immune-evasion-cluster table
and the clonal-complex rule.

The validation tests download public genomes of reference strains and assert
the published types. They are skipped unless ``HYDRA_VALIDATION=1``, and they
need BLAST+ on PATH and a built database directory (``HYDRA_SAUREUS_DB``, or
they build one, which takes a few minutes and about 6 MB).
"""

from __future__ import annotations

import os
import urllib.request
import zipfile
from pathlib import Path

import pytest

from hydra_amr.records import TypingResult
from hydra_amr.typing.organisms import saureus as sa
from hydra_amr.typing.organisms.base import OrganismDb, Reference

VALIDATION = os.environ.get("HYDRA_VALIDATION") == "1"
needs_validation = pytest.mark.skipif(
    not VALIDATION, reason="set HYDRA_VALIDATION=1 to run the download validation")


# --------------------------------------------------------------- module shape
def test_module_exposes_the_agreed_interface():
    assert sa.KEY == "saureus"
    assert sa.TITLE == "Staphylococcus aureus"
    assert "Staphylococcus aureus" in sa.SPECIES
    assert "saureus" in sa.MLST_SCHEMES
    assert all(isinstance(db, OrganismDb) for db in sa.DATABASES)
    assert all(isinstance(ref, Reference) for ref in sa.REFERENCES)
    assert callable(sa.fetch) and callable(sa.run)


def test_every_database_names_a_reference_that_exists():
    keys = {ref.key for ref in sa.REFERENCES}
    for db in sa.DATABASES:
        assert db.reference in keys, f"{db.name} cites unknown {db.reference}"
        assert db.redistribution in ("permitted", "restricted", "unstated")
        assert db.licence and db.version and db.source


def test_every_reference_is_identifiable_and_says_what_it_supports():
    for ref in sa.REFERENCES:
        assert ref.doi or ref.pmid, f"{ref.key} has neither DOI nor PMID"
        assert ref.supports, f"{ref.key} does not say what it supports"
        assert ref.authors and ref.title and ref.journal and ref.year


def test_marker_panel_is_coherent():
    names = [m.name for m in sa.MARKERS]
    assert len(names) == len(set(names)), "duplicate marker names"
    for marker in sa.MARKERS:
        assert marker.stop > marker.start, marker.name
        assert marker.strand in (1, 2), marker.name
        # Mostly RefSeq; the agr group IV window is a GenBank locus record,
        # because no complete genome has a published group IV assignment.
        assert marker.accession[:2] in ("NC", "AF"), marker.name
        assert marker.role, f"{marker.name} has no role"
        assert marker.product, f"{marker.name} has no product description"
    # Every gene set the callers read must be in the panel.
    for group in (sa.EGC_GENES, sa.HAEMOLYSINS, sa.IEC_GENES,
                  sa.ENTEROTOXINS_TYPED, ("blaZ", "blaR1", "blaI")):
        for gene in group:
            assert gene in names, f"{gene} is read by a caller but not in MARKERS"


def test_the_spa_database_is_flagged_as_not_redistributable():
    """The Ridom policy forbids shipping a copy inside a product."""
    spa = next(db for db in sa.DATABASES if db.name == "spa")
    assert spa.redistribution == "restricted"
    assert "non-commercial" in spa.licence
    assert "not be vendored" in spa.notes or "NOT be vendored" in spa.notes


# ------------------------------------------------------------- spa repeat data
def test_load_repeats_strips_the_r_prefix_and_survives_crlf(tmp_path):
    """The Ridom FASTA is CRLF with blank lines between records."""
    path = tmp_path / "sparepeats.fasta"
    path.write_bytes(b">r01\r\nGAGGAAGACAACAACAAGCCTAGC\r\n\r\n"
                     b">r02\r\nAAAGAAGACAACAAAAAACCTGGC\r\n\r\n"
                     b">r100\r\nAAAGAAGACGGCAACAAACCTGGC\r\n")
    repeats = sa._load_repeats(path)
    assert repeats["GAGGAAGACAACAACAAGCCTAGC"] == "01"
    assert repeats["AAAGAAGACAACAAAAAACCTGGC"] == "02"
    assert repeats["AAAGAAGACGGCAACAAACCTGGC"] == "100"
    assert all(not code.startswith("r") for code in repeats.values())


def test_load_spa_types_drops_the_NT_sentinel(tmp_path):
    """``NT,`` has an empty succession; keeping it would type every failure NT."""
    path = tmp_path / "spatypes.txt"
    path.write_bytes(b"t008,11-19-12-21-17-34-24-34-22-25\r\n"
                     b"t002,26-23-17-34-17-20-17-12-17-16\r\n"
                     b"NT,\r\n")
    types = sa._load_spa_types(path)
    assert types["11-19-12-21-17-34-24-34-22-25"] == "t008"
    assert "" not in types, "the empty succession must not map to a type"
    assert "NT" not in types.values()


def test_tile_repeats_reads_a_clean_array():
    repeats = {"AAAAAAAAAAAAAAAAAAAAAAAA": "11",
               "CCCCCCCCCCCCCCCCCCCCCCCC": "19",
               "GGGGGGGGGGGGGGGGGGGGGGGG": "12"}
    lengths = (24,)
    region = "AAAAAAAAAAAAAAAAAAAAAAAA" "CCCCCCCCCCCCCCCCCCCCCCCC" \
             "GGGGGGGGGGGGGGGGGGGGGGGG"
    codes, gaps, lead, trail = sa._tile_repeats(region, repeats, lengths)
    assert codes == ["11", "19", "12"]
    assert (gaps, lead, trail) == (0, 0, 0)


def test_tile_repeats_prefers_the_longest_match():
    """Short repeats are prefixes of long ones, so order of trial decides."""
    short = "ACGT" * 6              # 24 bp
    long = short + "TTT"            # 27 bp, with the 24-mer as its prefix
    repeats = {short: "44", long: "30"}
    lengths = tuple(sorted({len(short), len(long)}, reverse=True))
    codes, gaps, lead, trail = sa._tile_repeats(long, repeats, lengths)
    assert codes == ["30"], "a 24-bp prefix must not win over the real 27-bp repeat"
    assert (gaps, lead, trail) == (0, 0, 0)


def test_tile_repeats_counts_gaps_and_flanks():
    repeats = {"AAAAAAAAAAAAAAAAAAAAAAAA": "11",
               "CCCCCCCCCCCCCCCCCCCCCCCC": "19"}
    region = ("TTTTT"                              # 5 bp unrecognised lead
              "AAAAAAAAAAAAAAAAAAAAAAAA"
              "GTGTGTGT"                           # 8 bp unrecognised gap
              "CCCCCCCCCCCCCCCCCCCCCCCC"
              "TTTTTTT")                           # 7 bp unrecognised trail
    codes, gaps, lead, trail = sa._tile_repeats(region, repeats, (24,))
    assert codes == ["11", "19"]
    assert gaps == 1
    assert lead == 5
    assert trail == 7


def test_tile_repeats_on_an_array_with_no_known_repeat():
    codes, gaps, lead, trail = sa._tile_repeats("ACGT" * 20, {"AAAA": "01"}, (4,))
    assert codes == []
    assert lead == 80


# --------------------------------------------------------------- hit fixtures
def _hit(name, *, coverage=100.0, identity=100.0, ref_len=900, contig="contig1",
         contig_length=100000, strand="+", ref_blocks=None, contig_blocks=None,
         insert_gap=0, insert_at=0, n_alignments=1):
    return sa.MarkerHit(
        name=name, reference=f"NC_000000.1:1-{ref_len}", strain="synthetic",
        locus_tag="SYN_0001", role="test", product=f"{name} product",
        reference_length=ref_len, contig=contig, contig_length=contig_length,
        strand=strand, identity_pct=identity, coverage_pct=coverage,
        reference_blocks=ref_blocks or [(1, int(ref_len * coverage / 100))],
        contig_blocks=contig_blocks or [(1000, 1000 + int(ref_len * coverage / 100))],
        insert_gap=insert_gap, insert_at=insert_at, n_alignments=n_alignments)


def test_marker_hit_present_partial_and_absent_bands():
    assert _hit("x").present
    assert not _hit("x", identity=85.0).present, "81% lukF-PV cross-hits must fail"
    assert _hit("x", coverage=60.0).partial
    assert not _hit("x", coverage=60.0).present
    assert not _hit("x", coverage=10.0).partial


def test_marker_hit_detects_a_contig_edge():
    edge = _hit("x", contig_blocks=[(5, 500)], contig_length=600)
    middle = _hit("x", contig_blocks=[(1000, 1500)], contig_length=100000)
    assert edge.at_contig_edge()
    assert not middle.at_contig_edge()


# ------------------------------------------------------------- PVL and toxins
def test_pvl_needs_both_subunits():
    both = sa._pvl({"lukS-PV": _hit("lukS-PV"), "lukF-PV": _hit("lukF-PV")})
    assert both.call == "positive"
    assert both.loci_found == 2


def test_pvl_with_one_subunit_is_not_determined_not_positive():
    one = sa._pvl({"lukS-PV": _hit("lukS-PV")})
    assert one.call == "not determined"
    assert "contig break" in one.note


def test_pvl_rejects_the_lukF_family_cross_hit():
    """lukF-PV aligns to lukF/hlgB at about 81% in every PVL-negative genome."""
    result = sa._pvl({"lukF-PV": _hit("lukF-PV", identity=81.0)})
    assert result.call == "negative"


def test_pvl_negative_states_why_rather_than_a_bare_no():
    result = sa._pvl({})
    assert result.call == "negative"
    assert "absent" in result.note and "lukS-PV" in result.note


def test_tst_reports_a_partial_hit_as_undetermined():
    assert sa._tst({"tst": _hit("tst")}).call == "positive"
    assert sa._tst({}).call == "negative"
    assert sa._tst({"tst": _hit("tst", coverage=55.0)}).call == "not determined"


# ----------------------------------------------------------------- hlb lesion
def test_hlb_with_a_prophage_inside_it_is_not_reported_as_present():
    """Coverage alone says 100%; the two alignment pieces say otherwise."""
    hit = _hit("hlb", ref_len=993, insert_gap=43036, insert_at=183, n_alignments=2,
               contig_blocks=[(2084489, 2084671), (2127707, 2128530)])
    result = sa._hlb({"hlb": hit}, {})
    assert result.call == "disrupted by prophage integration"
    assert result.lineage == "Sa3int prophage"
    assert "43036" in result.note
    assert "183" in result.note


def test_hlb_contiguous_is_intact():
    result = sa._hlb({"hlb": _hit("hlb", ref_len=993)}, {})
    assert result.call == "intact"
    assert "no insertion" in result.note


def test_hlb_absent_is_flagged_as_unusual():
    result = sa._hlb({}, {})
    assert result.call == "absent"
    assert "core-genome" in result.note


# -------------------------------------------------- immune evasion cluster
@pytest.mark.parametrize("genes,expected", [
    (("sea", "sak", "chp", "scn"), "type A"),
    (("sak", "chp", "scn"), "type B"),
    (("chp", "scn"), "type C"),
    (("sea", "sak", "scn"), "type D"),
    (("sak", "scn"), "type E"),
    (("sep", "sak", "chp", "scn"), "type F"),
    (("sep", "sak", "scn"), "type G"),
])
def test_iec_types_follow_the_published_table(genes, expected):
    best = {gene: _hit(gene) for gene in genes}
    result = sa._iec(best, TypingResult(scheme="hlb", call="disrupted"))
    assert result.call == expected


def test_iec_without_scn_is_absent():
    result = sa._iec({"sea": _hit("sea")},
                     TypingResult(scheme="hlb", call="intact"))
    assert result.call == "absent"
    assert "agrees with the hlb call" in result.note


def test_iec_with_an_undescribed_combination_is_not_forced_into_a_type():
    best = {gene: _hit(gene) for gene in ("scn", "sea", "sep", "chp")}
    result = sa._iec(best, TypingResult(scheme="hlb", call="disrupted"))
    assert result.call == "not determined"
    assert "does not match any of the types" in result.note


def test_iec_notices_it_contradicts_an_intact_hlb():
    best = {gene: _hit(gene) for gene in ("sak", "chp", "scn")}
    result = sa._iec(best, TypingResult(scheme="hlb", call="intact"))
    assert result.call == "type B"
    assert "disagrees with the hlb call" in result.note


# ----------------------------------------------------------------- capsule
def test_capsule_serotypes_are_mutually_exclusive():
    five = {g: _hit(g) for g in ("cap5H", "cap5I", "cap5J")}
    assert sa._capsule(five).call == "5"
    eight = {g: _hit(g) for g in ("cap8H", "cap8I", "cap8J")}
    assert sa._capsule(eight).call == "8"


def test_capsule_warns_that_an_intact_locus_is_not_an_expressed_capsule():
    five = {g: _hit(g) for g in ("cap5H", "cap5I", "cap5J")}
    assert "not the same as" in sa._capsule(five).note


def test_capsule_with_a_split_locus_is_not_determined():
    partial = {"cap5H": _hit("cap5H"), "cap5I": _hit("cap5I", coverage=50.0)}
    assert sa._capsule(partial).call == "not determined"
    assert sa._capsule({}).call == "not determined"


# -------------------------------------------------------------------- ACME
def test_acme_separates_the_element_from_the_core_genome_operon():
    """The two arcA alleles are only 76% identical; that is the whole test."""
    present = sa._acme({"arcA-ACME": _hit("arcA-ACME"),
                        "arcA-native": _hit("arcA-native")})
    assert present.call == "arc cluster present"
    absent = sa._acme({"arcA-ACME": _hit("arcA-ACME", identity=76.0),
                       "arcA-native": _hit("arcA-native")})
    assert absent.call == "absent"
    assert "must not be" in absent.note


# ------------------------------------------------------------ bla lesions
def _orf(residues, *, stop_at=None, trailing_stop=True):
    """A synthetic coding sequence of *residues* codons."""
    codons = ["AAA"] * residues
    if stop_at is not None:
        codons[stop_at] = "TAA"
    return "ATG" + "".join(codons[1:]) + ("TAA" if trailing_stop else "")


def test_lesion_absent_when_there_is_no_hit():
    lesion = sa._describe_lesion(None, "blaR1", {}, 585)
    assert lesion.state == "absent"
    assert "not detected" in lesion.detail


def test_lesion_intact_reads_the_reading_frame():
    gene = _orf(126)
    contigs = {"c1": "N" * 100 + gene + "N" * 100}
    hit = _hit("blaI", ref_len=len(gene), contig="c1",
               contig_length=len(contigs["c1"]),
               contig_blocks=[(101, 100 + len(gene))])
    lesion = sa._describe_lesion(hit, "blaI", contigs, 126)
    assert lesion.state == "intact"
    assert "126 residues" in lesion.detail


def test_lesion_finds_a_premature_stop():
    gene = _orf(585, stop_at=100)
    contigs = {"c1": "N" * 100 + gene + "N" * 100}
    hit = _hit("blaR1", ref_len=len(gene), contig="c1",
               contig_length=len(contigs["c1"]),
               contig_blocks=[(101, 100 + len(gene))])
    lesion = sa._describe_lesion(hit, "blaR1", contigs, 585)
    assert lesion.state == "internal stop"
    assert "residue 101" in lesion.detail
    assert "no functional product" in lesion.detail


def test_lesion_finds_a_frameshift():
    gene = _orf(585)[:-1]           # one base short: frame no longer a multiple of 3
    contigs = {"c1": "N" * 100 + gene + "N" * 100}
    hit = _hit("blaR1", ref_len=len(gene), contig="c1",
               contig_length=len(contigs["c1"]),
               contig_blocks=[(101, 100 + len(gene))])
    lesion = sa._describe_lesion(hit, "blaR1", contigs, 585)
    assert lesion.state == "frameshifted"
    assert "multiple of three" in lesion.detail


def test_lesion_separates_a_contig_break_from_a_real_truncation():
    """Both lose 3' sequence; only one is evidence the gene is truncated."""
    gene = _orf(585)
    partial = len(gene) // 2
    # Cut by the end of the contig.
    broken = {"c1": gene[:partial]}
    at_edge = _hit("blaR1", ref_len=len(gene), coverage=50.0, contig="c1",
                   contig_length=partial, contig_blocks=[(1, partial)],
                   ref_blocks=[(1, partial)])
    lesion = sa._describe_lesion(at_edge, "blaR1", broken, 585)
    assert lesion.state == "truncated by contig break"
    assert "assembly break" in lesion.detail

    # Truncated with the contig carrying on past it.
    inside = {"c1": "N" * 100 + gene[:partial] + "N" * 5000}
    mid = _hit("blaR1", ref_len=len(gene), coverage=50.0, contig="c1",
               contig_length=len(inside["c1"]),
               contig_blocks=[(101, 100 + partial)], ref_blocks=[(1, partial)])
    lesion = sa._describe_lesion(mid, "blaR1", inside, 585)
    assert lesion.state == "truncated"
    assert "really is truncated" in lesion.detail


def test_lesion_reports_an_insertion_as_interrupted():
    hit = _hit("blaR1", ref_len=1758, insert_gap=5000, insert_at=400,
               n_alignments=2)
    lesion = sa._describe_lesion(hit, "blaR1", {}, 585)
    assert lesion.state == "interrupted"
    assert "inserted into the gene" in lesion.detail


def _blaz_protein(res128="T", res216="S"):
    """A synthetic 281-residue BlaZ with the class A motifs in place."""
    aa = ["A"] * 281
    aa[62] = "S"                                   # Ambler Ser70
    aa[120], aa[121], aa[122] = "S", "D", "N"      # Ambler SDN loop 130-132
    aa[224], aa[225], aa[226] = "K", "S", "G"      # Ambler K[ST]G 234-236
    aa[sa.BLAZ_AMBLER_128] = res128
    aa[sa.BLAZ_AMBLER_216] = res216
    return "".join(aa)


def _blaz_hit(protein):
    """A MarkerHit over a contig carrying the coding sequence for *protein*."""
    codon = {"A": "GCA", "S": "TCA", "D": "GAT", "N": "AAT", "K": "AAA",
             "G": "GGA", "T": "ACA", "L": "CTA", "W": "TGG"}
    gene = "".join(codon[a] for a in protein) + "TAA"
    contigs = {"c1": "N" * 100 + gene + "N" * 100}
    hit = _hit("blaZ", ref_len=len(gene), contig="c1",
               contig_length=len(contigs["c1"]),
               contig_blocks=[(101, 100 + len(gene))])
    return hit, contigs


@pytest.mark.parametrize("res128,res216,expected", [
    ("T", "S", "A"), ("K", "N", "B"), ("T", "N", "C"),
    ("A", "S", "D"), ("L", "S", "E"), ("T", "T", "F"),
])
def test_blaz_allotype_reads_ambler_128_and_216(res128, res216, expected):
    hit, contigs = _blaz_hit(_blaz_protein(res128, res216))
    allotype, detail = sa._blaz_allotype(hit, contigs)
    assert allotype == expected
    assert "Ambler 128" in detail and "Ambler 216" in detail


def test_blaz_allotype_refuses_an_undescribed_residue_pair():
    hit, contigs = _blaz_hit(_blaz_protein("W", "W"))
    allotype, detail = sa._blaz_allotype(hit, contigs)
    assert allotype == ""
    assert "not one of the described allotype combinations" in detail


def test_blaz_allotype_checks_the_class_a_motifs_before_reading_residues():
    """Ambler numbering is not a fixed offset, so the frame must be verified."""
    protein = list(_blaz_protein("T", "S"))
    protein[62] = "A"                      # break the catalytic serine
    hit, contigs = _blaz_hit("".join(protein))
    allotype, detail = sa._blaz_allotype(hit, contigs)
    assert allotype == ""
    assert "numbering cannot be trusted" in detail


def test_blaz_allotype_not_assigned_without_a_complete_gene():
    assert sa._blaz_allotype(None, {})[0] == ""
    interrupted = _hit("blaZ", insert_gap=5000, n_alignments=2)
    assert sa._blaz_allotype(interrupted, {})[0] == ""


def test_blaz_allotype_table_covers_only_described_combinations():
    assert set(sa.BLAZ_ALLOTYPE.values()) == {"A", "B", "C", "D", "E", "F"}
    assert sa.BLAZ_ALLOTYPE[("T", "S")] == "A", "type A is Thr128/Ser216"
    assert sa.BLAZ_ALLOTYPE[("K", "N")] == "B", "Lys128 is type B, not A"


def test_bla_operon_predicts_an_inducible_penicillinase():
    gene = {name: _orf(sa.BLA_EXPECTED[name]) for name in sa.BLA_EXPECTED}
    contigs = {}
    best = {}
    for name, seq in gene.items():
        contigs[name] = "N" * 100 + seq + "N" * 100
        best[name] = _hit(name, ref_len=len(seq), contig=name,
                          contig_length=len(contigs[name]),
                          contig_blocks=[(101, 100 + len(seq))])
    results = {r.scheme: r for r in sa._bla(best, contigs)}
    assert results["blaZ"].call == "present"
    assert results["blaZ"].lineage == "intact"
    assert results["bla_operon"].call == "blaR1 intact, blaI intact"
    assert "inducible" in results["bla_operon"].note


def test_bla_operon_flags_the_susceptible_blaZ_positive_configuration():
    """blaZ with a dead sensor and a live repressor cannot be induced.

    This is the case the owner asked for: a blaZ gene call alone would say
    resistant, and the isolate tests penicillin-susceptible.
    """
    blaz = _orf(sa.BLA_EXPECTED["blaZ"])
    blai = _orf(sa.BLA_EXPECTED["blaI"])
    blar = _orf(sa.BLA_EXPECTED["blaR1"], stop_at=50)
    contigs, best = {}, {}
    for name, seq in (("blaZ", blaz), ("blaI", blai), ("blaR1", blar)):
        contigs[name] = "N" * 100 + seq + "N" * 100
        best[name] = _hit(name, ref_len=len(seq), contig=name,
                          contig_length=len(contigs[name]),
                          contig_blocks=[(101, 100 + len(seq))])
    results = {r.scheme: r for r in sa._bla(best, contigs)}
    assert results["blaZ"].call == "present"
    operon = results["bla_operon"]
    assert operon.call == "blaR1 internal stop, blaI intact"
    assert "cannot be induced" in operon.note
    assert "susceptible" in operon.note
    assert "zone-edge" in operon.note


def test_bla_operon_without_blaZ_is_not_applicable():
    results = {r.scheme: r for r in sa._bla({}, {})}
    assert results["blaZ"].call == "absent"
    assert results["bla_operon"].call == "not applicable"
    assert "zone-edge" in results["blaZ"].note


def test_bla_operon_without_the_repressor_is_constitutive():
    blaz = _orf(sa.BLA_EXPECTED["blaZ"])
    blar = _orf(sa.BLA_EXPECTED["blaR1"])
    contigs, best = {}, {}
    for name, seq in (("blaZ", blaz), ("blaR1", blar)):
        contigs[name] = "N" * 100 + seq + "N" * 100
        best[name] = _hit(name, ref_len=len(seq), contig=name,
                          contig_length=len(contigs[name]),
                          contig_blocks=[(101, 100 + len(seq))])
    operon = {r.scheme: r for r in sa._bla(best, contigs)}["bla_operon"]
    assert operon.call == "blaI absent, blaR1 intact"
    assert "constitutive" in operon.note


# --------------------------------------------------------------- agr grouping
def test_agr_group_is_decided_on_coverage_as_well_as_identity():
    """A 99%-identity hit over 8% of the window must not beat a full-length one.

    This is the real failure seen on MRSA252, whose agr group III window
    matches over its whole length at 98.6% while a 107-bp stretch matches the
    group II reference at 99%.
    """
    best = {"agr-II": _hit("agr-II", identity=99.0, coverage=8.0, ref_len=1396),
            "agr-III": _hit("agr-III", identity=98.6, coverage=100.0, ref_len=1382)}
    result = sa._agr(best)
    assert result.call == "III"
    assert "next closest" in result.note


def test_agr_below_the_floor_mentions_the_missing_group_IV():
    best = {"agr-I": _hit("agr-I", identity=81.0, coverage=64.0)}
    result = sa._agr(best)
    assert result.call == "not determined"
    assert "group IV" in result.note


def test_agr_with_no_window_says_so():
    result = sa._agr({})
    assert result.call == "not determined"
    assert "absent" in result.note


def test_agr_notes_a_truncated_agrA():
    best = {"agr-I": _hit("agr-I", ref_len=1382),
            "agrA": _hit("agrA", coverage=45.0)}
    result = sa._agr(best)
    assert result.call == "I"
    assert "non-functional" in result.note


# ------------------------------------------------------------- clonal complex
@pytest.fixture
def profile_db(tmp_path):
    """A tiny PubMLST-shaped profile table."""
    root = tmp_path / "saureus" / "mlst"
    root.mkdir(parents=True)
    rows = [
        "ST\tarcC\taroE\tglpF\tgmk\tpta\ttpi\tyqiL",
        "1\t1\t1\t1\t1\t1\t1\t1",                  # CC1 founder
        "5\t1\t4\t1\t4\t12\t1\t10",                # CC5 founder
        "8\t3\t3\t1\t1\t4\t4\t3",                  # CC8 founder
        "250\t3\t3\t1\t1\t4\t4\t16",               # single-locus variant of ST8
        "30\t2\t2\t2\t2\t6\t3\t2",                 # CC30 founder
        "36\t2\t2\t2\t2\t3\t3\t2",                 # single-locus variant of ST30
        "9999\t91\t92\t93\t94\t95\t96\t97",        # unrelated to every founder
        "7777\t1\t1\t1\t1\t1\t1\t55",              # SLV of ST1 only
    ]
    (root / "profiles.tsv").write_text("\n".join(rows) + "\n")
    return tmp_path


def test_clonal_complex_names_the_founder(profile_db):
    result = sa.clonal_complex("8", profile_db)
    assert result.call == "CC8"
    assert "is the founder of" in result.note
    assert result.alleles["loci_shared_with_founder"] == "7/7"


def test_clonal_complex_accepts_a_single_locus_variant(profile_db):
    for st, want in (("250", "CC8"), ("36", "CC30"), ("7777", "CC1")):
        result = sa.clonal_complex(st, profile_db)
        assert result.call == want, st
        assert "single-locus variant" in result.note
        assert result.alleles["loci_shared_with_founder"] == "6/7"


def test_clonal_complex_refuses_to_force_a_distant_st(profile_db):
    result = sa.clonal_complex("9999", profile_db)
    assert result.call == "no established complex"
    assert result.lineage == "ST9999"
    assert "rather than forced into the nearest one" in result.note


def test_clonal_complex_needs_a_numeric_st(profile_db):
    for st in ("-", "", "novel", None):
        result = sa.clonal_complex(st, profile_db)
        assert result.call == "not determined"
        assert "numeric sequence type" in result.note


def test_clonal_complex_reports_an_unknown_st_as_such(profile_db):
    result = sa.clonal_complex("424242", profile_db)
    assert result.call == "not determined"
    assert "not in the installed PubMLST profile table" in result.note


def test_clonal_complex_without_a_database_says_so(tmp_path):
    result = sa.clonal_complex("8", tmp_path)
    assert result.call == "not determined"
    assert "not installed" in result.note


def test_cc_founders_are_unique_and_sorted():
    assert len(sa.CC_FOUNDERS) == len(set(sa.CC_FOUNDERS))
    assert list(sa.CC_FOUNDERS) == sorted(sa.CC_FOUNDERS)


# ------------------------------------------------------------- run() guards
def test_run_without_a_database_explains_itself(tmp_path):
    assembly = tmp_path / "a.fna"
    assembly.write_text(">c1\nACGT\n")
    results = sa.run(assembly, tmp_path / "missing", workdir=tmp_path / "w")
    assert len(results) == 1
    assert results[0].call == "not determined"
    assert "not installed" in results[0].note


def test_never_reports_a_bare_negative_without_a_reason():
    """Every 'negative' or 'absent' call must say what the evidence was."""
    for result in (sa._pvl({}), sa._tst({}), sa._acme({}), sa._hlb({}, {}),
                   sa._capsule({}), sa._iec({}, TypingResult(scheme="hlb",
                                                             call="intact"))):
        assert result.note, f"{result.scheme} gave a bare {result.call}"
        assert len(result.note) > 40, f"{result.scheme} note is too thin"


# ------------------------------------------------- validation on real genomes
#: accession -> (strain, expected calls). Sources for every expectation are in
#: docs/typing/saureus.md.
PANEL: dict[str, tuple[str, dict[str, str]]] = {
    "GCF_000013465.1": ("USA300_FPR3757", {
        "spa": "t008", "agr": "I", "capsule": "5", "pvl": "positive",
        "tst": "negative", "clonal_complex": "CC8",
        "hlb": "disrupted by prophage integration", "iec": "type B",
        "acme": "arc cluster present", "blaZ": "absent"}),
    "GCF_000009645.1": ("N315", {
        "spa": "t002", "agr": "II", "capsule": "5", "pvl": "negative",
        "tst": "positive", "clonal_complex": "CC5",
        "hlb": "disrupted by prophage integration", "iec": "type F",
        "acme": "absent", "blaZ": "present"}),
    # MW2's expected spa type is the one its own deposited closed genome
    # encodes. Rouha 2018 Table 1 lists strain MW2 as t125; the X region of
    # NC_003923.1 tiles with no gaps to 07-23-23-21-16-34-33-13 = t128, and
    # t125 (07-23-23-23-23-21-13) is not a repeat-count variant of it. See
    # docs/typing/saureus.md; this asserts the module reads the record right.
    "GCF_000011265.1": ("MW2", {
        "spa": "t128", "agr": "III", "capsule": "8", "pvl": "positive",
        "tst": "negative", "clonal_complex": "CC1",
        "hlb": "disrupted by prophage integration", "acme": "absent"}),
    "GCF_000010465.1": ("Newman", {
        "spa": "t008", "agr": "I", "capsule": "5", "pvl": "negative",
        "clonal_complex": "CC8", "iec": "type A", "acme": "absent"}),
    "GCF_000011505.1": ("MRSA252", {
        "spa": "t018", "agr": "III", "capsule": "8", "pvl": "negative",
        "clonal_complex": "CC30", "acme": "absent", "blaZ": "present"}),
    "GCF_000011525.1": ("MSSA476", {
        "spa": "t607", "agr": "III", "capsule": "8", "pvl": "negative",
        "clonal_complex": "CC1", "acme": "absent", "blaZ": "present"}),
    "GCF_000012045.1": ("COL", {
        "spa": "t008", "agr": "I", "capsule": "5", "pvl": "negative",
        "clonal_complex": "CC8", "hlb": "intact", "iec": "absent",
        "acme": "absent", "blaZ": "absent"}),
    "GCF_000027045.1": ("TW20", {
        "spa": "t037", "agr": "I", "capsule": "8", "pvl": "negative",
        "clonal_complex": "CC239", "acme": "absent", "blaZ": "present"}),
    # The awkward case: RF122 is bovine ST151 and its spa is a published
    # pseudogene, so the t-number must arrive carrying caveats, and lukM /
    # lukF'-PV must not be mistaken for PVL.
    "GCF_000009005.1": ("RF122", {
        "spa": "t529", "capsule": "8", "pvl": "negative", "tst": "positive",
        "clonal_complex": "CC151", "hlb": "intact", "iec": "absent",
        "acme": "absent", "blaZ": "absent"}),
    "GCF_000013425.1": ("NCTC8325", {
        "spa": "t211", "agr": "I", "capsule": "5", "pvl": "negative",
        "clonal_complex": "CC8", "hlb": "disrupted by prophage integration",
        "iec": "type B", "acme": "absent", "blaZ": "absent"}),
}

DATASETS = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
            "{acc}/download?include_annotation_type=GENOME_FASTA&filename={acc}.zip")


def _download_genome(accession: str, into: Path) -> Path:
    fasta = into / f"{accession}.fna"
    if fasta.exists():
        return fasta
    into.mkdir(parents=True, exist_ok=True)
    archive = into / f"{accession}.zip"
    request = urllib.request.Request(
        DATASETS.format(acc=accession), headers={"User-Agent": sa.USER_AGENT})
    with urllib.request.urlopen(request, timeout=300) as response:
        archive.write_bytes(response.read())
    with zipfile.ZipFile(archive) as zf:
        name = next(n for n in zf.namelist()
                    if n.endswith("_genomic.fna") and "cds_from" not in n)
        fasta.write_bytes(zf.read(name))
    archive.unlink()
    return fasta


@pytest.fixture(scope="session")
def validation_db(tmp_path_factory):
    existing = os.environ.get("HYDRA_SAUREUS_DB")
    if existing:
        return Path(existing)
    root = tmp_path_factory.mktemp("saureus-db")
    sa.fetch(root)
    return root


@pytest.fixture(scope="session")
def genome_dir(tmp_path_factory):
    existing = os.environ.get("HYDRA_SAUREUS_GENOMES")
    return Path(existing) if existing else tmp_path_factory.mktemp("saureus-genomes")


@needs_validation
@pytest.mark.parametrize("accession", list(PANEL))
def test_validation_panel(accession, validation_db, genome_dir, tmp_path):
    strain, expected = PANEL[accession]
    fasta = _download_genome(accession, genome_dir)
    results = {r.scheme: r for r in sa.run(fasta, validation_db, threads=2,
                                           workdir=tmp_path)}
    wrong = []
    for scheme, want in expected.items():
        got = results[scheme].call
        if got != want:
            wrong.append(f"{scheme}: expected {want!r}, got {got!r} "
                         f"({results[scheme].note[:160]})")
    assert not wrong, f"{strain} ({accession}):\n  " + "\n  ".join(wrong)


@needs_validation
@pytest.mark.parametrize("accession,expected", [
    ("GCF_000009645.1", "C"), ("GCF_000011505.1", "A"),
    ("GCF_000011525.1", "C"), ("GCF_000027045.1", "A"),
])
def test_validation_blaz_allotypes_are_canonical(accession, expected,
                                                 validation_db, genome_dir,
                                                 tmp_path):
    """Each blaZ-positive panel genome must pass the motif checks and type."""
    fasta = _download_genome(accession, genome_dir)
    results = {r.scheme: r for r in sa.run(fasta, validation_db, threads=2,
                                           workdir=tmp_path)}
    blaz = results["blaZ"]
    assert blaz.call == "present"
    assert blaz.alleles["blaZ_allotype"] == expected
    assert blaz.lineage == f"allotype {expected}"


@needs_validation
def test_validation_rf122_spa_call_carries_its_caveats(validation_db, genome_dir,
                                                       tmp_path):
    """A real t-number off a degraded array must not arrive unqualified."""
    fasta = _download_genome("GCF_000009005.1", genome_dir)
    results = {r.scheme: r for r in sa.run(fasta, validation_db, threads=2,
                                           workdir=tmp_path)}
    spa = results["spa"]
    assert spa.call == "t529"
    assert "only 2 repeat(s)" in spa.note
    assert "spa gene itself is not intact" in spa.note
    # lukM / lukF'-PV of the bovine lineage must not read as PVL.
    assert results["pvl"].call == "negative"


@needs_validation
def test_validation_spa_successions_are_in_the_ridom_table(validation_db,
                                                           genome_dir, tmp_path):
    """A t-number must come from the table, never from a partial succession."""
    types = sa._load_spa_types(validation_db / "saureus" / "spa" / "spatypes.txt")
    for accession, (strain, expected) in PANEL.items():
        if "spa" not in expected:
            continue
        fasta = _download_genome(accession, genome_dir)
        results = {r.scheme: r for r in sa.run(fasta, validation_db, threads=2,
                                               workdir=tmp_path / accession)}
        spa = results["spa"]
        assert types.get(spa.lineage) == spa.call, strain
        assert spa.alleles["repeat_succession"] == spa.lineage
