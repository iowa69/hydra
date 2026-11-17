"""Tests for the Escherichia coli / Shigella typing module.

The unit tests build small synthetic fixtures and never touch the network. The
validation tests at the bottom download public reference genomes and assert the
calls against reference-tool output; they run only with ``HYDRA_VALIDATION=1``
set, because they need about 50 MB of downloads and BLAST+.
"""

from __future__ import annotations

import os
import re
import shutil
import urllib.request
import zipfile
from pathlib import Path

import pytest

from hydra_amr.engines.blast import Hsp
from hydra_amr.typing.organisms import ecoli
from hydra_amr.typing.organisms.base import OrganismDb, Reference


# ------------------------------------------------------------ module contract
def test_module_exposes_the_agreed_api():
    assert ecoli.KEY == "ecoli"
    assert ecoli.TITLE
    assert isinstance(ecoli.SPECIES, tuple) and "Escherichia coli" in ecoli.SPECIES
    assert isinstance(ecoli.MLST_SCHEMES, tuple) and ecoli.MLST_SCHEMES
    assert all(isinstance(db, OrganismDb) for db in ecoli.DATABASES)
    assert all(isinstance(ref, Reference) for ref in ecoli.REFERENCES)
    assert callable(ecoli.fetch) and callable(ecoli.run)


def test_shigella_species_trigger_the_module():
    # Shigella is an E. coli clade, every scheme here applies to it, and the
    # shigella_eiec call exists to correct the species label.
    for name in ("Shigella flexneri", "Shigella sonnei", "Shigella dysenteriae",
                 "Shigella boydii"):
        assert name in ecoli.SPECIES


def test_every_reference_is_citable():
    for ref in ecoli.REFERENCES:
        assert ref.doi or ref.pmid, f"{ref.key} has neither a DOI nor a PMID"
        assert ref.supports, f"{ref.key} does not say what it supports"
        assert ref.authors.endswith("."), f"{ref.key} author list looks truncated"
        assert 1990 < ref.year < 2030


def test_every_database_names_a_real_reference_and_states_its_licence():
    keys = {ref.key for ref in ecoli.REFERENCES}
    for db in ecoli.DATABASES:
        assert db.reference in keys, f"{db.name} cites unknown {db.reference}"
        assert db.licence, f"{db.name} states no licence"
        assert db.redistribution in ("permitted", "restricted", "unstated")
        assert db.source.startswith("https://")
        assert db.version


# -------------------------------------------------------------- phylogroup
def _concrete(primer: str) -> str:
    """A concrete sequence matching a degenerate primer (first base of each set)."""
    out = []
    for base in primer:
        expanded = ecoli.IUPAC.get(base)
        out.append(expanded[1] if expanded else base)
    return "".join(out)


def _amplicon(name: str, filler: str = "A") -> str:
    """A synthetic PCR product of the published length for one reaction."""
    forward, reverse, size = ecoli.PRIMERS[name]
    head = _concrete(forward)
    tail = ecoli._revcomp(_concrete(reverse))
    middle = filler * (size - len(head) - len(tail))
    assert len(head) + len(middle) + len(tail) == size
    return head + middle + tail


def test_iupac_expansion_matches_a_concrete_sequence():
    pattern = ecoli._primer_regex("AAYRCYATTCGCCAGCTTGC")
    assert pattern.fullmatch(_concrete("AAYRCYATTCGCCAGCTTGC"))
    assert pattern.fullmatch("AACGCCATTCGCCAGCTTGC")
    assert not pattern.fullmatch("AAGGCCATTCGCCAGCTTGC")  # Y cannot be G


def test_revcomp_handles_ambiguity_codes():
    assert ecoli._revcomp("ACGT") == "ACGT"
    assert ecoli._revcomp("AAAA") == "TTTT"
    assert ecoli._revcomp("ARYG") == "CRYT"


def test_in_silico_pcr_finds_a_product_of_the_right_length():
    contig = ("G" * 500) + _amplicon("arpA") + ("G" * 500)
    result = ecoli._in_silico_pcr([(">c1", contig)], "arpA", allow_partial=False)
    assert result.hit
    assert "400 bp" in result.evidence


def test_in_silico_pcr_finds_the_product_on_the_reverse_strand():
    contig = ecoli._revcomp(("G" * 300) + _amplicon("chuA") + ("G" * 300))
    result = ecoli._in_silico_pcr([(">c1", contig)], "chuA", allow_partial=False)
    assert result.hit


def test_in_silico_pcr_rejects_a_product_of_the_wrong_length():
    forward, reverse, _size = ecoli.PRIMERS["yjaA"]
    contig = (("G" * 100) + _concrete(forward) + ("A" * 900)
              + ecoli._revcomp(_concrete(reverse)) + ("G" * 100))
    result = ecoli._in_silico_pcr([(">c1", contig)], "yjaA", allow_partial=False)
    assert not result.hit
    assert "expected 211 bp" in result.evidence


def test_in_silico_pcr_accepts_a_partial_product_only_when_asked():
    # One primer, and the contig ends before a full product could fit: a
    # draft assembly broken through the amplicon.
    forward, _reverse, _size = ecoli.PRIMERS["TspE4.C2"]
    contig = ("G" * 200) + _concrete(forward) + ("A" * 20)
    assert not ecoli._in_silico_pcr([(">c", contig)], "TspE4.C2",
                                    allow_partial=False).hit
    partial = ecoli._in_silico_pcr([(">c", contig)], "TspE4.C2",
                                   allow_partial=True)
    assert partial.hit and "partial" in partial.evidence


def test_quadruplex_table_covers_every_pattern():
    # All sixteen +/- combinations of the four markers must have an entry, or a
    # real isolate could fall through the table.
    for arpa in (0, 1):
        for chua in (0, 1):
            for yjaa in (0, 1):
                for tspe4 in (0, 1):
                    assert (arpa, chua, yjaa, tspe4) in ecoli.QUADRUPLEX


def test_phylogroup_calls_group_a_from_the_published_pattern():
    contigs = [(">c1", ("G" * 200) + _amplicon("arpA") + ("G" * 200)
                + _amplicon("trpBA_control") + ("G" * 200))]
    result = ecoli._phylogroup(contigs)
    # arpA+ chuA- yjaA- TspE4.C2- is group A, and no follow-up reaction is run.
    assert result.call == "A"
    assert result.alleles["arpA"] == "+"
    assert result.alleles["chuA"] == "-"
    assert result.loci_found == 1
    assert result.loci_total == 4
    assert "commensal" in result.lineage


def test_phylogroup_calls_b1_and_reports_the_marker_pattern():
    contigs = [(">c1", _amplicon("arpA") + ("G" * 300) + _amplicon("TspE4.C2"))]
    result = ecoli._phylogroup(contigs)
    assert result.call == "B1"
    assert result.alleles["TspE4.C2"] == "+"
    assert "arpA+" in result.alleles["quadruplex_pattern"]


def test_phylogroup_resolves_c_with_the_trpa_reaction():
    contigs = [(">c1", _amplicon("arpA") + ("G" * 200) + _amplicon("yjaA")
                + ("G" * 200) + _amplicon("trpBA_control") + ("G" * 200)
                + _amplicon("trpA_gpC"))]
    result = ecoli._phylogroup(contigs)
    assert result.call == "C"
    assert result.alleles["trpA_gpC"] == "+"


def test_phylogroup_falls_back_to_a_when_the_c_reaction_is_negative():
    contigs = [(">c1", _amplicon("arpA") + ("G" * 200) + _amplicon("yjaA")
                + ("G" * 200) + _amplicon("trpBA_control"))]
    assert ecoli._phylogroup(contigs).call == "A"


def test_phylogroup_resolves_e_with_the_arpa_group_e_reaction():
    contigs = [(">c1", _amplicon("arpA") + ("G" * 200) + _amplicon("chuA")
                + ("G" * 200) + _amplicon("trpBA_control") + ("G" * 200)
                + _amplicon("arpA_gpE"))]
    assert ecoli._phylogroup(contigs).call == "E"


def test_phylogroup_resolves_g_with_the_ybgd_reaction():
    contigs = [(">c1", _amplicon("chuA") + ("G" * 200) + _amplicon("ybgD_gpG")
                + ("G" * 200) + _amplicon("trpBA_control"))]
    result = ecoli._phylogroup(contigs)
    assert result.call == "G"
    assert "high virulence" in result.lineage


def test_phylogroup_refuses_to_refine_without_the_internal_control():
    # arpA+ chuA+ is D/E, and the group E reaction is read against trpBA. With
    # no control product the reaction cannot be read, and guessing D would be
    # wrong half the time.
    contigs = [(">c1", _amplicon("arpA") + ("G" * 200) + _amplicon("chuA"))]
    result = ecoli._phylogroup(contigs)
    assert result.call == "D/E (control failed)"
    assert "internal control" in result.note


def test_phylogroup_says_not_determined_when_every_marker_is_negative():
    result = ecoli._phylogroup([(">c1", "ACGT" * 2000)])
    assert result.call == "not determined"
    assert "cryptic clade" in result.note
    assert result.loci_found == 0


# ---------------------------------------------------------------------- stx
def test_stx_header_parsing():
    assert ecoli._parse_stx_header("AAA98151.1|stxA1a|stxA1a") == ("A", "1a", "1a")
    # stx2a, stx2c and stx2d share one generalised class.
    assert ecoli._parse_stx_header("ABR09924.1|stxA2a|stxA2") == ("A", "2a", "2")
    assert ecoli._parse_stx_header("BAB82998.1|stxB2b|stxB2b") == ("B", "2b", "2b")
    assert ecoli._parse_stx_header("rubbish") is None


def _hsp(sseq: str, qseq: str, sstart: int) -> Hsp:
    return Hsp(qseqid="c1", sseqid="x|stxA2a|stxA2", pident=99.0,
               length=len(sseq), mismatch=0, gapopen=0, qstart=1,
               qend=3 * len(sseq), sstart=sstart,
               send=sstart + len(sseq.replace("-", "")) - 1, evalue=0.0,
               bitscore=600.0, slen=320, qlen=5000, gaps=0,
               nident=len(sseq), strand="+", qseq=qseq, sseq=sseq)


def test_residue_mapping_reads_the_reference_position():
    # Reference positions 10..14, so position 12 is the third column.
    hsp = _hsp(sseq="ABCDE", qseq="ABXDE", sstart=10)
    assert ecoli._residues_from_hsp(hsp, (10, 12, 14)) == {10: "A", 12: "X", 14: "E"}


def test_residue_mapping_skips_gaps_in_the_reference():
    hsp = _hsp(sseq="AB-CD", qseq="ABXYD", sstart=1)
    # The gap column consumes no reference position, so 3 is the 4th column.
    assert ecoli._residues_from_hsp(hsp, (3,)) == {3: "Y"}


def _stx_hit(subunit: str, subtype: str, *, coverage: float = 100.0,
             identity: float = 100.0, residues: dict[int, str] | None = None,
             lo: int = 1, hi: int = 960, contig: str = "c1",
             strand: str = "+") -> ecoli._StxHit:
    reflen = 320 if subunit == "A" else 90
    stx_class = "2" if subtype in ("2a", "2c", "2d") else subtype
    return ecoli._StxHit(
        subunit=subunit, subtype=subtype, stx_class=stx_class,
        reference="ACC1", contig=contig, lo=lo, hi=hi, strand=strand,
        nident=int(reflen * coverage / 100.0 * identity / 100.0), reflen=reflen,
        coverage_pct=coverage, identity_pct=identity, residues=residues or {})


def test_stx1a_operon_is_called_from_a_complete_pair():
    call, status, note = ecoli._call_stx_operon(
        _stx_hit("A", "1a"), _stx_hit("B", "1a"), gap=12)
    assert (call, status) == ("stx1a", "COMPLETE")
    assert note == ""


@pytest.mark.parametrize("res313,res319,res35,expected", [
    ("F", "K", "D", "stx2a"),
    ("S", "E", "D", "stx2a"),
    ("F", "K", "N", "stx2c"),
    ("F", "E", "N", "stx2c"),
    ("S", "E", "N", "stx2d"),
])
def test_stx2acd_is_separated_on_three_residues(res313, res319, res35, expected):
    a = _stx_hit("A", "2a", residues={313: res313, 319: res319})
    b = _stx_hit("B", "2a", residues={35: res35})
    call, status, _note = ecoli._call_stx_operon(a, b, gap=12)
    assert (call, status) == (expected, "COMPLETE")


def test_stx2acd_with_an_unlisted_residue_pattern_is_novel_not_guessed():
    a = _stx_hit("A", "2a", residues={313: "W", 319: "W"})
    b = _stx_hit("B", "2a", residues={35: "W"})
    call, status, note = ecoli._call_stx_operon(a, b, gap=12)
    assert call == "stx2 novel"
    assert status == "COMPLETE_NOVEL"
    assert "313=W" in note


def test_a_partial_operon_is_never_given_a_subtype():
    a = _stx_hit("A", "2a", coverage=61.0, residues={313: "F", 319: "K"})
    b = _stx_hit("B", "2a", residues={35: "D"})
    call, status, note = ecoli._call_stx_operon(a, b, gap=12)
    assert status == "PARTIAL"
    assert call == "stx2 partial"
    assert "not subtyped" in note
    assert "61.0 %" in note


def test_an_operon_of_two_different_types_is_reported_as_novel():
    call, status, note = ecoli._call_stx_operon(
        _stx_hit("A", "2b"), _stx_hit("B", "2e"), gap=10)
    assert status == "COMPLETE_NOVEL"
    assert "different stx types" in note


def test_too_long_an_intergenic_region_is_not_one_operon():
    call, status, note = ecoli._call_stx_operon(
        _stx_hit("A", "1a"), _stx_hit("B", "1a"),
        gap=ecoli.STX_INTERGENIC_MAX + 1)
    assert status == "COMPLETE_NOVEL"
    assert "intergenic" in note


def test_identity_below_the_subtype_cutoff_is_novel():
    # stx1 needs 98.3 % combined identity; 95 % is a novel subtype, not stx1a.
    call, status, note = ecoli._call_stx_operon(
        _stx_hit("A", "1a", identity=95.0), _stx_hit("B", "1a", identity=95.0),
        gap=10)
    assert status == "COMPLETE_NOVEL"
    assert "98.3 %" in note


def test_the_published_identity_cutoffs_are_the_ones_implemented():
    assert ecoli.STX_IDENTITY["1a"] == 0.983
    assert ecoli.STX_IDENTITY["2k"] == 0.985
    assert ecoli.STX_IDENTITY["2l"] == 0.985
    assert ecoli.STX_IDENTITY["2"] == 0.98
    # Every subtype in the current nomenclature has a cut-off.
    for subtype in ("1a", "1c", "1d", "1e", "2b", "2e", "2f", "2g", "2h",
                    "2i", "2j", "2k", "2l", "2m", "2n", "2o"):
        assert subtype in ecoli.STX_IDENTITY


def test_pairing_keeps_subunits_on_the_same_contig_and_strand():
    a = _stx_hit("A", "2a", lo=100, hi=1059)
    b_right = _stx_hit("B", "2a", lo=1080, hi=1349)
    b_other_contig = _stx_hit("B", "2a", lo=1080, hi=1349, contig="c2")
    b_other_strand = _stx_hit("B", "2a", lo=1080, hi=1349, strand="-")
    operons, orphans = ecoli._pair_stx([a, b_right])
    assert len(operons) == 1 and operons[0][2] == 20
    assert not orphans
    for wrong in (b_other_contig, b_other_strand):
        operons, orphans = ecoli._pair_stx([a, wrong])
        assert not operons
        assert len(orphans) == 2


def test_a_lone_subunit_is_an_orphan_not_an_operon():
    operons, orphans = ecoli._pair_stx([_stx_hit("A", "2a")])
    assert not operons
    assert len(orphans) == 1 and orphans[0].subunit == "A"


# -------------------------------------------------- pathotype and virulence
def _genes(*names: str) -> dict[str, ecoli._Gene]:
    return {name: ecoli._Gene(name=name, identity_pct=100.0, coverage_pct=100.0,
                              contig="c1", span="1-900(+)", reference=f"{name}:1:X")
            for name in names}


class _FakeStx:
    def __init__(self, call: str):
        self.call = call


def _no_shigella():
    from hydra_amr.records import TypingResult
    return TypingResult(scheme="shigella_eiec", call="not Shigella/EIEC",
                        alleles={"organism": "Escherichia coli"})


def _shigella(organism: str):
    from hydra_amr.records import TypingResult
    return TypingResult(scheme="shigella_eiec", call="C3",
                        alleles={"organism": organism})


def test_typical_epec_needs_eae_and_bfpa():
    result = ecoli._pathotype(_genes("eae-a01-alpha_1", "bfpA", "tir"),
                              _FakeStx("none detected"), _no_shigella())
    assert result.call == "typical EPEC"
    assert "bfpA" in result.lineage


def test_atypical_epec_is_eae_without_bfpa():
    result = ecoli._pathotype(_genes("eae-g01-gamma_1", "tir"),
                              _FakeStx("none detected"), _no_shigella())
    assert result.call == "atypical EPEC"
    assert "bfpA absent" in result.lineage


def test_ehec_is_stx_plus_eae():
    result = ecoli._pathotype(_genes("eae-g01-gamma_1", "tir", "ehxA"),
                              _FakeStx("stx1a + stx2a"), _no_shigella())
    assert result.call == "EHEC"
    assert "enterohaemorrhagic" in result.lineage


def test_stec_without_eae_says_so():
    result = ecoli._pathotype(_genes("ehxA"), _FakeStx("stx2a"), _no_shigella())
    assert result.call == "STEC (eae-negative)"


def test_the_stec_eaec_hybrid_is_recognised_and_not_double_reported():
    result = ecoli._pathotype(_genes("aggR", "aggA", "aatA"),
                              _FakeStx("stx2a"), _no_shigella())
    assert result.call == "STEC/EAEC hybrid"
    assert "O104:H4" in result.lineage


def test_etec_from_either_toxin():
    for toxin in ("eltIAB-1", "estah-STa2"):
        result = ecoli._pathotype(_genes(toxin), _FakeStx("none detected"),
                                  _no_shigella())
        assert result.call == "ETEC"


def test_typical_eaec_is_aggr_positive():
    result = ecoli._pathotype(_genes("aggR", "aggA", "aap"),
                              _FakeStx("none detected"), _no_shigella())
    assert result.call == "EAEC (typical)"


def test_aggr_negative_eaec_is_flagged_as_contested():
    result = ecoli._pathotype(_genes("aafA", "aatA"),
                              _FakeStx("none detected"), _no_shigella())
    assert "aggR-negative" in result.call
    assert "contested" in result.lineage


def test_shigella_and_eiec_are_separated_by_the_cluster_call():
    shig = ecoli._pathotype(_genes("ipaH", "ipaD", "virF"),
                            _FakeStx("none detected"), _shigella("Shigella"))
    assert shig.call == "Shigella"
    eiec = ecoli._pathotype(_genes("ipaH", "ipaD", "virF"),
                            _FakeStx("none detected"), _shigella("EIEC"))
    assert eiec.call == "EIEC"
    both = ecoli._pathotype(_genes("ipaH"), _FakeStx("none detected"),
                            _shigella("ambiguous"))
    assert both.call == "EIEC/Shigella"
    assert "not separated" in both.lineage


def test_expec_needs_two_of_the_five_marker_groups():
    one = ecoli._pathotype(_genes("iutA"), _FakeStx("none detected"),
                           _no_shigella())
    assert "ExPEC" not in one.call
    two = ecoli._pathotype(_genes("iutA", "kpsMII_K1"),
                           _FakeStx("none detected"), _no_shigella())
    assert two.call == "ExPEC"
    assert "2 of the 5" in two.lineage


def test_afa_is_not_called_daec_in_an_expec_genome():
    # afa/draBC is itself an ExPEC marker; EC958 and its like must not be
    # reported as possible DAEC.
    result = ecoli._pathotype(_genes("afaB", "afaC", "iutA", "papC"),
                              _FakeStx("none detected"), _no_shigella())
    assert result.call == "ExPEC"
    assert "not evidence of DAEC" in result.lineage


def test_afa_alone_is_a_daec_candidate_not_a_call():
    result = ecoli._pathotype(_genes("afaB", "afaC"),
                              _FakeStx("none detected"), _no_shigella())
    assert result.call == "possible DAEC"
    assert "candidate, not a call" in result.lineage


def test_a_genome_with_no_markers_gets_an_explanation_not_a_bare_negative():
    result = ecoli._pathotype({}, _FakeStx("none detected"), _no_shigella())
    assert result.call == "no pathotype markers"
    assert "does not mean the isolate is harmless" in result.note


# -------------------------------------------------------------- intimin, VFs
def test_intimin_subtype_names_are_read_from_the_allele():
    assert ecoli._eae_subtype("eae-g01-gamma_1") == "gamma1"
    assert ecoli._eae_subtype("eae-b01a-beta_1a") == "beta1a"
    assert ecoli._eae_subtype("eae-a03-zeta") == "zeta"
    assert ecoli._eae_subtype("eae-e07-xi") == "xi"


def test_intimin_reports_the_subtype_with_its_evidence():
    result = ecoli._intimin(_genes("eae-g01-gamma_1", "tir"))
    assert result.call == "gamma1"
    assert result.alleles["identity_pct"] == "100.00"
    assert result.alleles["tir"] == "tir"


def test_intimin_marks_a_nearest_neighbour_rather_than_naming_a_subtype():
    genes = {"eae-g01-gamma_1": ecoli._Gene("eae-g01-gamma_1", 93.0, 100.0,
                                            "c1", "1-900(+)", "x")}
    result = ecoli._intimin(genes)
    assert result.call == "gamma1*"
    assert "nearest neighbour" in result.note


def test_intimin_absence_is_explained():
    result = ecoli._intimin(_genes("tir"))
    assert result.call == "eae absent"
    assert "attaching-and-effacing" in result.note


def test_missing_tir_beside_a_present_eae_is_flagged():
    result = ecoli._intimin(_genes("eae-g01-gamma_1"))
    assert "no tir allele" in result.note


def test_the_expec_panel_applies_both_published_rules():
    result = ecoli._expec_panel(_genes(
        "papC", "papA_F11", "sfaD", "kpsMII_K1", "iutA",
        "vat", "fyuA", "chuA", "yfcV", "hlyA", "usp"))
    assert "ExPEC (4/5)" in result.call
    assert "UPEC by the four-gene rule" in result.call
    assert result.alleles["UPEC_four_gene_rule"].startswith("4/4")
    assert "adhesins" in result.alleles


def test_the_expec_panel_states_what_was_not_met():
    result = ecoli._expec_panel(_genes("chuA", "fyuA"))
    assert result.call.startswith("not ExPEC (0/5 marker groups)")
    assert "2/4" in result.alleles["UPEC_four_gene_rule"]


# ----------------------------------------------------------------- fimH, ST131
def test_fimh_allele_names_are_spelled_as_the_literature_does():
    assert ecoli._fimh_name("FimH30") == "fimH30"
    assert ecoli._fimh_name("fimH1") == "fimH1"
    assert ecoli._fimh_name("FimH3363:1:X") == "fimH3363"


def test_st131_subclone_is_only_claimed_for_a_published_clade_allele():
    from hydra_amr.records import TypingResult
    for allele, expected in (("fimH30", "H30 (clade C) if ST131"),
                             ("fimH22", "H22 (clade B) if ST131"),
                             ("fimH41", "H41 (clade A) if ST131")):
        result = ecoli._st131_subclone(TypingResult(scheme="fimH", call=allele))
        assert result.call == expected
        assert "conditional on the ST" in result.note


def test_st131_subclone_refuses_a_clade_for_an_unlisted_allele():
    from hydra_amr.records import TypingResult
    result = ecoli._st131_subclone(TypingResult(scheme="fimH", call="fimH27"))
    assert result.call == "not applicable"
    assert "read the MLST call" in result.note


def test_st131_h30_names_the_unresolvable_c1_c2_split():
    from hydra_amr.records import TypingResult
    result = ecoli._st131_subclone(TypingResult(scheme="fimH", call="fimH30"))
    assert "H30Rx" in result.note
    assert "core-genome SNPs" in result.note


def test_st131_subclone_handles_a_missing_fimh_call():
    from hydra_amr.records import TypingResult
    for call in ("fimH not found", "not determined"):
        result = ecoli._st131_subclone(TypingResult(scheme="fimH", call=call))
        assert result.call == "not determined"


# ------------------------------------------------------------ Shigella / EIEC
def test_cluster_table_covers_every_cluster_in_the_published_definitions():
    # Every cluster this module can call must have a Shigella/EIEC verdict and
    # a description, or the call would be a bare cluster code.
    for cluster, (organism, description) in ecoli.SHIGEI_CLUSTERS.items():
        assert organism in ("Shigella", "EIEC"), cluster
        assert description


def test_published_cluster_overlaps_are_resolved_and_others_are_not():
    assert ecoli._resolve_cluster_overlap(["C1", "CSB12"]) == "CSB12"
    assert ecoli._resolve_cluster_overlap(["C5", "C3"]) == "C3"
    assert ecoli._resolve_cluster_overlap(["C2", "CSS"]) == "C2"
    # No published rule for this pair, so it stays ambiguous.
    assert ecoli._resolve_cluster_overlap(["C4", "C9"]) == ""


# -------------------------------------------------------------------- errors
def test_run_refuses_an_empty_assembly(tmp_path):
    from hydra_amr.utils import HydraError
    empty = tmp_path / "empty.fna"
    empty.write_text("")
    with pytest.raises(HydraError):
        ecoli.run(empty, tmp_path, workdir=tmp_path)


def test_missing_database_is_named_in_the_error(tmp_path):
    from hydra_amr.utils import HydraError
    with pytest.raises(HydraError) as excinfo:
        ecoli._require_db(tmp_path, "fimtyper/fimH.fsa", "the FimTyper database")
    assert "fetch()" in str(excinfo.value)
    assert "FimTyper" in str(excinfo.value)


# ------------------------------------------------------------------ validation
#: accession -> (name, expected calls, where the expectation comes from).
#: Phylogroups are EzClermont 1.0.0 output; stx subtypes are NCBI StxTyper
#: v1.0.45 output; the rest is the published characterisation of the strain.
VALIDATION: tuple[tuple[str, str, dict[str, str], str], ...] = (
    ("GCA_000005845.2", "MG1655",
     {"phylogroup": "A", "stx": "none detected",
      "pathotype": "no pathotype markers", "intimin": "eae absent",
      "shigella_eiec": "not Shigella/EIEC"},
     "E. coli K-12 MG1655; phylogroup by EzClermont 1.0.0, stx by StxTyper 1.0.45"),
    ("GCA_000008865.2", "Sakai",
     {"phylogroup": "E", "stx": "stx1a + stx2a", "pathotype": "EHEC",
      "intimin": "gamma1", "shigella_eiec": "not Shigella/EIEC"},
     "O157:H7 Sakai RIMD 0509952 (Hayashi 2001, PMID 11258796); phylogroup by "
     "EzClermont, stx by StxTyper, intimin gamma1 published for O157:H7 "
     "(Oswald 2000, PMID 10603369)"),
    ("GCA_000732965.1", "EDL933",
     {"phylogroup": "E", "stx": "stx1a + stx2a", "pathotype": "EHEC",
      "intimin": "gamma1"},
     "O157:H7 EDL933; phylogroup by EzClermont, stx by StxTyper"),
    ("GCA_000285655.3", "EC958",
     {"phylogroup": "B2", "fimH": "fimH30",
      "ST131_subclone": "H30 (clade C) if ST131", "stx": "none detected",
      "pathotype": "ExPEC"},
     "O25b:H4-ST131 H30Rx EC958; phylogroup by EzClermont, fimH30 published "
     "(Price 2013, PMID 24345742)"),
    ("GCA_000210475.1", "H10407",
     {"phylogroup": "A", "stx": "none detected", "pathotype": "ETEC"},
     "ETEC H10407 (Crossman 2010, PMID 20802035); phylogroup by EzClermont"),
    ("GCA_000027125.1", "EAEC042",
     {"phylogroup": "D", "stx": "none detected",
      "pathotype": "EAEC (typical)"},
     "EAEC 042 (O44:H18), the prototypical aggR-positive EAEC; phylogroup by "
     "EzClermont"),
    ("GCA_000299455.1", "O104H4",
     {"phylogroup": "B1", "stx": "stx2a",
      "pathotype": "STEC/EAEC hybrid + ExPEC"},
     "O104:H4 2011 German outbreak clone; phylogroup by EzClermont, stx2a by "
     "StxTyper"),
    ("GCA_000006925.2", "Sflex301",
     {"phylogroup": "A", "stx": "none detected", "pathotype": "Shigella",
      "shigella_eiec": "C3 (Shigella)"},
     "Shigella flexneri 2a str. 301 (Jin 2002, PMID 12384590); cluster C3 is "
     "the S. flexneri 1-5/X/Y cluster of Zhang 2021"),
    ("GCA_000010485.1", "SE15",
     {"phylogroup": "B2", "fimH": "fimH41",
      "ST131_subclone": "H41 (clade A) if ST131", "stx": "none detected"},
     "E. coli SE15, an ST131 clade A commensal; phylogroup by EzClermont"),
)

_NCBI_DOWNLOAD = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/"
                  "accession/{acc}/download?include_annotation_type=GENOME_FASTA")

validation = pytest.mark.skipif(
    os.environ.get("HYDRA_VALIDATION") != "1",
    reason="set HYDRA_VALIDATION=1 to download reference genomes and validate")


def _genome(accession: str, name: str, into: Path) -> Path:
    target = into / f"{name}.fna"
    if target.exists():
        return target
    archive = into / f"{name}.zip"
    request = urllib.request.Request(
        _NCBI_DOWNLOAD.format(acc=accession),
        headers={"User-Agent": ecoli.USER_AGENT})
    with urllib.request.urlopen(request, timeout=600) as response, \
            open(archive, "wb") as handle:
        shutil.copyfileobj(response, handle)
    with zipfile.ZipFile(archive) as zf:
        members = [m for m in zf.namelist() if m.endswith(".fna")]
        assert members, f"{accession} download held no FASTA"
        with zf.open(members[0]) as src, open(target, "wb") as dst:
            shutil.copyfileobj(src, dst)
    archive.unlink()
    return target


@pytest.fixture(scope="session")
def validation_root(tmp_path_factory) -> Path:
    root = Path(os.environ.get("HYDRA_VALIDATION_DIR")
                or tmp_path_factory.mktemp("ecoli-validation"))
    (root / "genomes").mkdir(parents=True, exist_ok=True)
    ecoli.fetch(root / "db")
    return root


@validation
@pytest.mark.parametrize("accession,name,expected,source", VALIDATION,
                         ids=[v[1] for v in VALIDATION])
def test_validation_genomes(validation_root, accession, name, expected, source):
    assembly = _genome(accession, name, validation_root / "genomes")
    work = validation_root / "work" / name
    work.mkdir(parents=True, exist_ok=True)
    results = {r.scheme: r for r in ecoli.run(assembly, validation_root / "db",
                                              threads=2, workdir=work)}
    for scheme, want in expected.items():
        assert scheme in results, f"{name}: {scheme} not reported"
        assert results[scheme].call == want, (
            f"{name} {scheme}: expected {want!r} (source: {source}), "
            f"got {results[scheme].call!r}")


@validation
def test_validation_every_call_states_its_evidence(validation_root):
    """No call may be a bare word: the owner rejected output that names a gene
    without saying what was found. Every scheme must carry either the evidence
    in ``alleles`` or an explanation in ``note``."""
    assembly = _genome("GCA_000008865.2", "Sakai", validation_root / "genomes")
    work = validation_root / "work" / "Sakai"
    work.mkdir(parents=True, exist_ok=True)
    for result in ecoli.run(assembly, validation_root / "db", threads=2,
                            workdir=work):
        assert result.call and result.call != "-", result.scheme
        assert result.alleles or result.note, (
            f"{result.scheme} reports {result.call!r} with no evidence "
            f"and no explanation")


@validation
def test_validation_fetch_is_idempotent(validation_root):
    before = {p: p.stat().st_mtime for p in
              (validation_root / "db").rglob("*") if p.is_file()}
    ecoli.fetch(validation_root / "db")
    after = {p: p.stat().st_mtime for p in
             (validation_root / "db").rglob("*") if p.is_file()}
    assert set(before) <= set(after)
    # provenance.tsv is rewritten every time; nothing else should be.
    changed = {p.name for p in before if before[p] != after.get(p)}
    assert changed <= {"provenance.tsv"}


@validation
def test_validation_provenance_records_every_database(validation_root):
    provenance = (validation_root / "db" / ecoli.KEY / "provenance.tsv")
    text = provenance.read_text(encoding="utf-8")
    for db in ecoli.DATABASES:
        assert db.name in text
        assert db.source in text
    assert len(text.strip().splitlines()) == len(ecoli.DATABASES) + 1
