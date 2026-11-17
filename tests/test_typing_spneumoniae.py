"""Tests for the *Streptococcus pneumoniae* typing module.

The fast tests exercise the decision rules on synthetic probe and locus
evidence: no network, no BLAST, no reference data. They are where the rules are
pinned, because the rules are the part a future edit can quietly break.

The validation tests (``HYDRA_VALIDATION=1``) build the real databases and type
public assemblies of characterised strains. See ``docs/typing/spneumoniae.md``
for the expectations and where each one comes from.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from hydra_amr.typing.organisms import spneumoniae as spn
from hydra_amr.typing.organisms.base import OrganismDb, Reference

# --------------------------------------------------------------------- helpers


def probe(gene: str, *, exact: bool = True, cluster: str = "", allele: str = "",
          coverage: float = 100.0) -> spn.ProbeHit:
    """A screening probe hit, present and either exact or one mismatch off."""
    return spn.ProbeHit(
        cluster=cluster or gene, gene=gene, allele=allele or f"{gene}-1",
        identity=100.0 if exact else 99.0, coverage=coverage,
        contig="contig1", start=1, end=100, strand="+")


def probes(*specs) -> dict[str, spn.ProbeHit]:
    """``probes("WZY4", ("WCJE9V", False))`` -> the probe dictionary."""
    out: dict[str, spn.ProbeHit] = {}
    for spec in specs:
        if isinstance(spec, tuple):
            out[spec[0]] = probe(spec[0], exact=spec[1])
        else:
            out[spec] = probe(spec)
    return out


def locus(serotype: str, coverage: float = 100.0, identity: float = 98.0,
          length: int = 20000, intervals=None) -> spn.LocusMatch:
    return spn.LocusMatch(
        serotype=serotype, length=length, coverage=coverage, identity=identity,
        intervals=intervals if intervals is not None else [(1, length)])


# ----------------------------------------------------------- module contract


def test_module_exposes_the_interface_the_integrator_calls():
    assert spn.KEY == "spneumoniae"
    assert spn.TITLE == "Streptococcus pneumoniae"
    assert spn.SPECIES == ("Streptococcus pneumoniae",)
    assert spn.MLST_SCHEMES == ("spneumoniae",)
    assert all(isinstance(db, OrganismDb) for db in spn.DATABASES)
    assert all(isinstance(ref, Reference) for ref in spn.REFERENCES)
    assert callable(spn.fetch) and callable(spn.run)


def test_every_database_records_its_licence_and_redistribution_status():
    for db in spn.DATABASES:
        assert db.source.startswith("https://")
        assert db.version and db.licence
        assert db.redistribution in {"permitted", "restricted", "unstated"}
        assert db.reference in {ref.key for ref in spn.REFERENCES}


def test_every_reference_says_what_it_supports_and_is_identifiable():
    keys = [ref.key for ref in spn.REFERENCES]
    assert len(keys) == len(set(keys)), "duplicate reference keys"
    for ref in spn.REFERENCES:
        assert ref.supports, f"{ref.key} does not say what it supports"
        assert ref.authors and ref.title and ref.journal
        # Everything except the CLSI standard is a journal article with an id.
        if ref.key != "CLSI":
            assert ref.doi or ref.pmid, f"{ref.key} has neither DOI nor PMID"


def test_every_source_file_is_attributed_to_a_declared_database():
    names = {db.name for db in spn.DATABASES}
    for _url, _path, _dbtype, name in spn._FILES:
        assert name in names


# ------------------------------------------------------------------- vaccines


def test_vaccine_formulations_have_the_valency_their_name_claims():
    expected = {"PCV7": 7, "PCV10-GSK": 10, "PCV13": 13, "PCV15": 15,
                "PCV20": 20, "PPSV23": 23}
    for name, count in expected.items():
        assert len(spn.VACCINES[name]) == count, name
    # PCV21 names 22 serotypes for 21 valencies: 15B and 15C share one
    # component in the licensed formulation.
    assert len(spn.VACCINES["PCV21"]) == 22
    assert {"15B", "15C"} <= set(spn.VACCINES["PCV21"])


def test_each_formulation_contains_the_one_before_it():
    assert set(spn.VACCINES["PCV7"]) < set(spn.VACCINES["PCV13"])
    assert set(spn.VACCINES["PCV13"]) < set(spn.VACCINES["PCV15"])
    assert set(spn.VACCINES["PCV15"]) < set(spn.VACCINES["PCV20"])


def test_vaccine_membership_of_a_single_serotype():
    included, partial, cross = spn.vaccine_membership(("19A",))
    assert (partial, cross) == ([], [])
    assert included == ["PCV13", "PCV15", "PCV20", "PCV21", "PPSV23"]


def test_a_non_vaccine_serotype_is_reported_as_in_nothing():
    assert spn.vaccine_membership(("23B",))[0] == ["PCV21"]
    assert spn.vaccine_membership(("16F",))[0] == ["PCV21"]


def test_6c_is_in_no_formulation_but_cross_reacts_with_the_6a_conjugate():
    included, partial, cross = spn.vaccine_membership(("6C",))
    assert included == []
    assert partial == []
    assert cross == ["PCV13 (cross-reactive with 6C)",
                     "PCV15 (cross-reactive with 6C)",
                     "PCV20 (cross-reactive with 6C)"]


def test_an_unresolved_pair_reports_partial_cover_not_full_cover():
    """The clinical point of not resolving 9A/9V: 9V is in six formulations
    and 9A is in none, so neither "covered" nor "not covered" is true."""
    included, partial, _cross = spn.vaccine_membership(("9A", "9V"))
    assert included == []
    assert partial == ["PCV7 (covers 9V only)", "PCV10-GSK (covers 9V only)",
                       "PCV13 (covers 9V only)", "PCV15 (covers 9V only)",
                       "PCV20 (covers 9V only)", "PPSV23 (covers 9V only)"]


def test_serogroup_is_read_off_the_serotype_name():
    assert spn.serogroup_of("19A") == "19"
    assert spn.serogroup_of("6C") == "6"
    assert spn.serogroup_of("1") == "1"
    assert spn.serogroup_of("35B") == "35"
    assert "tts" in spn.serogroup_of("37")


# ------------------------------------------------------- serotype probe rules


def test_a_single_target_probe_names_its_serotype():
    assert ("4", "WZY4 exact, 100.0% of WZY4-1") in \
        spn._serotype_from_probes(probes("WZY4"))


def test_a_single_target_probe_is_accepted_with_mismatches():
    """The CDC rule for these probes allows differences within the probe."""
    calls = dict(spn._serotype_from_probes(probes(("WZY23F", False))))
    assert "23F" in calls


def test_serogroup_19_probes_must_match_exactly():
    """19A and 19F wzy alleles differ by a few bases, so a near match is not
    evidence and the rule requires no differences at all."""
    assert spn._serotype_from_probes(probes(("WZY19F", False))) == []
    assert dict(spn._serotype_from_probes(probes("WZY19F"))).get("19F")


@pytest.mark.parametrize("wcip,wcin,expected", [
    ("WCIP6AC", "WCIN6AB", "6A"),
    ("WCIP6AC", "WCIN6CD", "6C"),
    ("WCIP6BD", "WCIN6AB", "6B"),
    ("WCIP6BD", "WCIN6CD", "6D"),
])
def test_serogroup_6_is_resolved_by_wcip_crossed_with_wcin(wcip, wcin, expected):
    """wciP residue 195 separates 6A/6C from 6B/6D, and the wciN allele
    separates 6A/6B from 6C/6D; together they give all four serotypes."""
    calls = dict(spn._serotype_from_probes(probes(wcip, wcin)))
    assert expected in calls
    assert wcip in calls[expected] and wcin in calls[expected]


def test_a_non_exact_wcip_does_not_call_a_serogroup_6_serotype():
    assert spn._serotype_from_probes(probes(("WCIP6AC", False), "WCIN6AB")) == []


def test_15b_and_15c_are_always_reported_together():
    """The discriminator is a TA repeat whose length an assembler collapses,
    so neither member may be claimed from an assembly - with the wciZ intact
    or not."""
    for exact in (True, False):
        calls = dict(spn._serotype_from_probes(
            probes("WZY15B", ("WCIZ15B", exact))))
        assert "15B/15C" in calls
        assert "15B" not in calls and "15C" not in calls


def test_15b_and_15c_are_reported_together_when_wciz_is_absent_entirely():
    calls = dict(spn._serotype_from_probes(probes("WZY15B")))
    assert "15B/15C" in calls


@pytest.mark.parametrize("wzy,gene,single,pair", [
    ("WZY9V", "WCJE9V", "9V", "9A/9V"),
    ("WZY7F", "WCWD7F", "7F", "7A/7F"),
    ("WZY18C", "WCIX18C", "18C", "18B/18C"),
])
def test_an_intact_acetyltransferase_gives_one_serotype_and_a_damaged_one_a_pair(
        wzy, gene, single, pair):
    assert single in dict(spn._serotype_from_probes(probes(wzy, gene)))
    assert pair in dict(spn._serotype_from_probes(probes(wzy, (gene, False))))


def test_35b_becomes_a_pair_only_when_wcig_differs():
    assert "35B" in dict(spn._serotype_from_probes(probes("WZY35B", "WCIG35B")))
    assert "35B/35D" in dict(
        spn._serotype_from_probes(probes("WZY35B", ("WCIG35B", False))))


def test_33f_becomes_33a_only_on_an_exact_wcje():
    assert "33A" in dict(spn._serotype_from_probes(probes("WZY33F", "WCJE33A")))
    assert "33F" in dict(spn._serotype_from_probes(probes("WZY33F")))


def test_9l_and_9n_are_separated_by_wcja():
    assert "9L" in dict(spn._serotype_from_probes(probes("WZY9N", "WCJA9L")))
    assert "9N" in dict(spn._serotype_from_probes(probes("WZY9N", "WCJA9N")))


def test_22a_and_22f_are_separated_by_wcwa():
    assert "22A" in dict(spn._serotype_from_probes(probes("WZY22F", "WCWA22A")))
    assert "22F" in dict(spn._serotype_from_probes(probes("WZY22F", "WCWA22F")))


def test_7c_needs_wchf_and_7b_is_its_absence():
    assert "7C" in dict(spn._serotype_from_probes(probes("WZY7C", "WCHF7C")))
    assert "7B" in dict(spn._serotype_from_probes(probes("WZY7C")))


def test_groups_the_probe_set_cannot_split_are_named_as_groups():
    for gene, expected in (("WZY11B:11C", "11B/11C"),
                           ("WZY24F:24A:24B", "24A/24B/24F"),
                           ("WCII12A:12B:46", "12A/12B/46"),
                           ("WCYE25A:25F", "25A/25F"),
                           ("WZY35C:42", "35C/42")):
        assert expected in dict(spn._serotype_from_probes(probes(gene)))


def test_serotype_37_is_called_from_tts_outside_the_capsule_locus():
    assert "37" in dict(spn._serotype_from_probes(probes("TTS")))


# -------------------------------------------------------- the serotype result


def test_a_clean_serotype_call_carries_its_evidence_and_vaccine_membership():
    result = spn.call_serotype([locus("19A"), locus("19F", coverage=78.0)],
                               probes("WZY19A"))
    assert result.scheme == "cps_serotype"
    assert result.call == "19A"
    assert result.lineage == "serogroup 19"
    assert result.alleles["serotype_source"] == "CDC serotype probe rules"
    assert "PCV13" in result.alleles["vaccines_included_in"]
    assert "19A locus" in result.alleles["cps_locus_best_match"]
    assert "19F locus" in result.alleles["cps_locus_runner_up"]


def test_a_pair_call_explains_in_the_note_why_it_is_a_pair():
    result = spn.call_serotype([locus("15B")],
                               probes("WZY15B", ("WCIZ15B", False)))
    assert result.call == "15B/15C"
    assert "wciZ" in result.note
    assert "TA" in result.note
    assert "reported as a group" in result.note


def test_a_serotype_outside_the_probe_set_says_so_and_does_not_pretend():
    result = spn.call_serotype([locus("17F"), locus("17A", coverage=60.0)], {})
    assert result.call == "17F"
    assert result.alleles["serotype_source"] == "whole cps locus coverage only"
    assert "No CDC probe rule covers this serotype" in result.note
    assert "not a serologically verified serotype" in result.note


def test_serogroups_with_no_published_split_are_reported_as_the_group():
    for serotype, expected in (("32F", "32A/32F"), ("24F", "24A/24B/24F")):
        result = spn.call_serotype([locus(serotype)], {})
        assert result.call == expected, serotype
        assert result.alleles["group_reason"]


def test_two_near_tied_loci_in_one_serogroup_are_reported_together():
    result = spn.call_serotype(
        [locus("16F", coverage=99.0), locus("16A", coverage=98.5)], {})
    assert result.call == "16F/16A"
    assert "within" in result.alleles["group_reason"]


def test_no_capsule_locus_at_all_is_reported_as_non_encapsulated():
    result = spn.call_serotype([locus("14", coverage=22.0)], {})
    assert result.call == "non-encapsulated"
    assert result.lineage == "null capsule"
    assert "not covered by any conjugate or polysaccharide vaccine" in result.note
    # The clades are named so the reader knows what was not determined.
    assert "NCC1" in result.note and "pspK" in result.note


def test_a_locus_present_but_gutted_is_also_non_encapsulated():
    """R6's shape: flanks on both sides of a hole where the serotype-specific
    genes were, and no serotype probe matching."""
    result = spn.call_serotype(
        [locus("2", coverage=63.6, length=20602,
               intervals=[(1, 2557), (10063, 20602)])], {})
    assert result.call == "non-encapsulated"
    assert result.lineage == "disrupted serogroup 2 locus"
    assert "2558-10062" in result.alleles["cps_locus_gaps"]


def test_a_probe_rule_overrides_a_short_locus_match_rather_than_being_lost():
    result = spn.call_serotype([locus("14", coverage=70.0)], probes("WZY14"))
    assert result.call == "14"
    assert "coverage_warning" in " ".join(result.alleles)


def test_disagreeing_rules_are_both_reported_not_silently_resolved():
    result = spn.call_serotype([locus("19A")], probes("WZY19A", "WZY23F"))
    assert set(result.call.split("/")) == {"19A", "23F"}
    assert "More than one CDC rule fired" in result.note


def test_11a_warns_that_11e_would_be_reported_as_11a():
    result = spn.call_serotype([locus("11A")], probes("WZY11A"))
    assert result.call == "11A"
    assert "11E" in result.note and "wcjE" in result.note


# -------------------------------------------------------------- PBP + the MIC


def test_the_pbp_regions_are_the_published_coordinates():
    by_gene = {region.gene: region for region in spn.PBP_REGIONS}
    assert (by_gene["pbp1a"].first_residue, by_gene["pbp1a"].length) == (371, 276)
    assert (by_gene["pbp2b"].first_residue, by_gene["pbp2b"].length) == (379, 277)
    assert (by_gene["pbp2x"].first_residue, by_gene["pbp2x"].length) == (229, 358)


def test_the_active_site_motifs_sit_where_the_literature_puts_them():
    motifs = {region.gene: dict(region.motifs) for region in spn.PBP_REGIONS}
    assert motifs["pbp1a"] == {"STMK": 370, "SSN": 466, "KTG": 557}
    # pbp2b's SXXK motif is SVVK, not STMK, and its KTG begins at 615.
    assert motifs["pbp2b"] == {"SVVK": 386, "SSN": 443, "KTG": 615}
    assert motifs["pbp2x"] == {"STMK": 337, "SSN": 395, "KSG": 547}


def test_pbp1a_catalytic_serine_lies_outside_the_typed_fragment():
    """S370 is one residue before the region, so it can never be read from an
    allele; the module must not claim otherwise."""
    region = next(r for r in spn.PBP_REGIONS if r.gene == "pbp1a")
    assert region.motifs[0] == ("STMK", 370)
    assert region.first_residue == 371


def test_every_key_site_lies_inside_its_region():
    for region in spn.PBP_REGIONS:
        last = region.first_residue + region.length - 1
        for position, _wild, label in region.key_sites:
            assert region.first_residue <= position <= last, \
                f"{region.gene} {position}"
            assert label


def test_a_complete_triplet_is_rendered_pbp1a_colon_2b_colon_2x():
    calls = {"1A": spn.PbpCall(locus="1A", allele="13"),
             "2B": spn.PbpCall(locus="2B", allele="16"),
             "2X": spn.PbpCall(locus="2X", allele="47")}
    result = spn.call_pbp(calls)
    assert result.call == "13:16:47"
    assert result.lineage == "in reference set"
    assert result.alleles["pbp1a"] == "13"
    assert result.loci_found == 3


def test_a_novel_region_is_reported_new_and_never_rounded_to_a_number():
    calls = {
        "1A": spn.PbpCall(locus="1A", allele="NEW", identity=99.3,
                          coverage=100.0, nearest="281",
                          substitutions=["T371S", "T574N"],
                          motif_changes=["T371S, in the STMK motif (370-373)"]),
        "2B": spn.PbpCall(locus="2B", allele="348"),
        "2X": spn.PbpCall(locus="2X", allele="620"),
    }
    result = spn.call_pbp(calls)
    assert result.call == "NEW:348:620"
    assert result.lineage == "novel allele present"
    assert "281" in result.alleles["pbp1a_nearest"]
    assert "T371S" in result.alleles["pbp1a_substitutions_vs_allele_0"]
    assert "STMK" in result.alleles["pbp1a_active_site"]
    assert "not rounded to the nearest allele" in result.note


def test_a_novel_region_with_no_active_site_change_says_that_explicitly():
    calls = {"1A": spn.PbpCall(locus="1A", allele="NEW", nearest="4",
                               substitutions=["A500V"]),
             "2B": spn.PbpCall(locus="2B", allele="0"),
             "2X": spn.PbpCall(locus="2X", allele="0")}
    result = spn.call_pbp(calls)
    assert "no substitution in the SXXK" in result.alleles["pbp1a_active_site"]


def test_an_unrecovered_region_is_a_dash_and_the_triplet_is_incomplete():
    calls = {"1A": spn.PbpCall(locus="1A", allele="-", note="not recovered"),
             "2B": spn.PbpCall(locus="2B", allele="57"),
             "2X": spn.PbpCall(locus="2X", allele="56")}
    result = spn.call_pbp(calls)
    assert result.call == "-:57:56"
    assert result.lineage == "incomplete"
    assert "Only 2 of the 3 regions were typed" in result.note


def test_a_known_triplet_yields_the_reference_mic_and_both_penicillin_verdicts(
        tmp_path):
    table = tmp_path / "mic.csv"
    table.write_text(
        '"APT","PEN","AMO","MER","TAX","CFT","CFX"\n'
        '"13-16-47","2","1","0.5","1","1","4"\n', encoding="utf-8")
    calls = {"1A": spn.PbpCall(locus="1A", allele="13"),
             "2B": spn.PbpCall(locus="2B", allele="16"),
             "2X": spn.PbpCall(locus="2X", allele="47")}
    result = spn.call_mic("13:16:47", calls, table)
    assert result.alleles["penicillin_predicted_mic_mg_L"] == "2"
    verdicts = result.alleles["penicillin_predicted_category"]
    # Penicillin at 2 mg/L is resistant for meningitis, susceptible otherwise:
    # the two criteria must both be stated, because the module cannot know
    # which applies.
    assert "meningitis (S <=0.06, R >=0.12 mg/L): R" in verdicts
    assert "non-meningitis parenteral (S <=2, I 4, R >=8 mg/L): S" in verdicts
    assert result.alleles["cefuroxime_predicted_category"].endswith("R")
    assert result.alleles["basis"].startswith("predicted")
    assert "not measured" in result.note


def test_the_predicted_mic_never_claims_to_be_a_measurement(tmp_path):
    table = tmp_path / "mic.csv"
    table.write_text('"APT","PEN"\n"0-0-0","0.03"\n', encoding="utf-8")
    calls = {k: spn.PbpCall(locus=k, allele="0") for k in ("1A", "2B", "2X")}
    result = spn.call_mic("0:0:0", calls, table)
    assert "*predicted*" in result.note
    assert "must not be reported as a susceptibility test result" in result.note
    assert "Li et al. 2016" in result.note


def test_a_novel_region_gets_no_mic_at_all(tmp_path):
    table = tmp_path / "mic.csv"
    table.write_text('"APT","PEN"\n"0-0-0","0.03"\n', encoding="utf-8")
    calls = {"1A": spn.PbpCall(locus="1A", allele="NEW"),
             "2B": spn.PbpCall(locus="2B", allele="0"),
             "2X": spn.PbpCall(locus="2X", allele="0")}
    result = spn.call_mic("NEW:0:0", calls, table)
    assert result.call == "not determined"
    assert "must be measured" in result.note
    assert not any(k.endswith("_mic_mg_L") for k in result.alleles)


def test_a_triplet_absent_from_the_table_is_declined_not_extrapolated(tmp_path):
    table = tmp_path / "mic.csv"
    table.write_text('"APT","PEN"\n"0-0-0","0.03"\n', encoding="utf-8")
    calls = {"1A": spn.PbpCall(locus="1A", allele="196"),
             "2B": spn.PbpCall(locus="2B", allele="0"),
             "2X": spn.PbpCall(locus="2X", allele="2")}
    result = spn.call_mic("196:0:2", calls, table)
    assert result.call == "not determined"
    assert result.lineage == "PBP type not in the reference table"
    assert "Extrapolating from a neighbouring triplet is not valid" in result.note


@pytest.mark.parametrize("code,mic,expected", [
    ("PEN", 0.03, ["meningitis", "S"]),
    ("PEN", 0.12, ["meningitis", "R"]),
    ("PEN", 8.0, ["non-meningitis", "R"]),
    ("MER", 0.25, ["S"]),
    ("MER", 0.5, ["I"]),
    ("MER", 1.0, ["R"]),
    ("CFX", 0.5, ["S"]),
    ("CFX", 2.0, ["R"]),
])
def test_clsi_breakpoints_land_where_m100_puts_them(code, mic, expected):
    verdicts = " ".join(spn._interpret(code, mic))
    for token in expected:
        assert token in verdicts


# ---------------------------------------------------------------- other calls


def test_both_pilus_islets_are_reported_from_their_marker_probes():
    result = spn.call_pilus(probes("RRGA", "PITB"))
    assert result.call == "PI-1+PI-2"
    assert result.loci_found == 2
    result = spn.call_pilus(probes("RRGA"))
    assert result.call == "PI-1"
    assert "not detected" in result.alleles["PI-2"]
    result = spn.call_pilus({})
    assert result.call == "none"
    assert "commonest state" in result.note


def test_ermb_predicts_clindamycin_resistance_and_mef_does_not():
    ermb = spn.call_macrolide(probes("ERMB"))
    assert ermb.call == "erm(B)"
    assert "clindamycin resistant" in ermb.lineage
    assert "MLS-B" in ermb.lineage

    mef = spn.call_macrolide(probes("MEF"))
    assert mef.call == "mef"
    assert "M phenotype" in mef.lineage
    assert "clindamycin active" in mef.lineage
    # Honest about what the probe cannot do.
    assert "does not separate mef(A)" in mef.alleles["mef"] or \
        "does not separate mef(A)" in mef.note


def test_ermb_plus_mef_is_reported_as_both_mechanisms():
    result = spn.call_macrolide(probes("ERMB", "MEF"))
    assert "erm(B)" in result.call and "mef" in result.call
    assert "resistant" in result.lineage


def test_erm_tr_warns_that_a_susceptible_clindamycin_result_is_unreliable():
    result = spn.call_macrolide(probes("ERMTR"))
    assert "inducible" in result.alleles["erm(TR) (erm(A) subclass TR)"] or \
        "inducible" in result.note
    assert "should not be used" in result.note


def test_lsac_gives_clindamycin_resistance_without_macrolide_resistance():
    result = spn.call_macrolide(probes("LSAC"))
    assert "lsa(C)" in result.call
    assert "LS-A" in result.lineage


def test_no_macrolide_determinant_is_not_reported_as_plain_susceptible():
    result = spn.call_macrolide({})
    assert result.call == "no acquired determinant detected"
    assert "not a susceptibility result" in result.note
    assert "23S rRNA" in result.note


def test_tetracycline_and_phenicol_determinants_are_interpreted():
    result = spn.call_other_determinants(probes("TETM", "CAT"))
    assert "tet(M)" in result.call and "cat" in result.call
    assert "doxycycline" in result.alleles["tet(M)"]
    assert spn.call_other_determinants({}).call == "none detected"


# -------------------------------------------- target-gene substitution calling


def test_substitutions_are_numbered_on_the_full_protein():
    region = next(r for r in spn.TARGET_REGIONS if r.gene == "folA")
    observed = list(region.peptide)
    observed[19] = "L"            # the 20th residue of an 81-100 window
    subs = spn._substitutions("".join(observed), region.peptide, region.offset,
                              region.neutral)
    assert subs == ["I100L"]


def test_the_polymorphisms_the_cdc_rules_ignore_are_not_called_resistance():
    region = next(r for r in spn.TARGET_REGIONS if r.gene == "parC")
    observed = list(region.peptide)
    observed[50] = "D"            # G51D locally, a known neutral variant
    assert spn._substitutions("".join(observed), region.peptide, region.offset,
                              region.neutral) == []


@pytest.mark.parametrize("gene,index,substitution", [
    ("gyrA", 10, "S81F"),
    ("parC", 1, "S79F"),
    ("parC", 5, "D83N"),
])
def test_the_canonical_qrdr_substitutions_come_out_with_their_known_names(
        gene, index, substitution):
    region = next(r for r in spn.TARGET_REGIONS if r.gene == gene)
    observed = list(region.peptide)
    observed[index] = substitution[-1]
    subs = spn._substitutions("".join(observed), region.peptide, region.offset,
                              region.neutral)
    assert subs == [substitution]


def test_a_wild_type_genome_is_reported_as_wild_type_over_named_residues():
    extracted = {
        region.probe: spn.Extracted(
            probe=region.probe, contig="c", start=1, end=1, strand="+",
            nucleotides=_back_translate(region.peptide), coverage=100.0,
            identity=100.0, complete=True)
        for region in spn.TARGET_REGIONS
    }
    result = spn.call_target_substitutions(extracted, {})
    assert result.call == "wild type at every target screened"
    assert result.alleles["gyrA"] == "wild type over residues 71-96"
    assert result.alleles["folA"] == "wild type over residues 81-100"
    assert "folP insertion" in result.note


def test_a_region_that_was_not_recovered_is_not_determined_rather_than_negative():
    result = spn.call_target_substitutions({}, {})
    assert "not determined" in result.alleles["gyrA"]
    assert "Not determined" in result.note


def test_an_implausibly_divergent_region_is_rejected_as_a_bad_extraction():
    region = next(r for r in spn.TARGET_REGIONS if r.gene == "gyrA")
    nonsense = "A" * len(region.peptide)
    extracted = {region.probe: spn.Extracted(
        probe=region.probe, contig="c", start=1, end=1, strand="+",
        nucleotides=_back_translate(nonsense), coverage=100.0, identity=100.0,
        complete=True)}
    result = spn.call_target_substitutions(extracted, {})
    assert "the extraction, not the allele, is wrong" in result.alleles["gyrA"]


_ONE_CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT",
    "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTT", "M": "ATG", "N": "AAT",
    "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT",
    "W": "TGG", "Y": "TAT",
}


def _back_translate(peptide: str) -> str:
    return "".join(_ONE_CODON[aa] for aa in peptide)


def test_the_fixture_back_translation_round_trips():
    """Guards the tests above: if this breaks they are testing nothing."""
    from hydra_amr.utils import translate
    for region in spn.TARGET_REGIONS:
        assert translate(_back_translate(region.peptide),
                         table_start_is_met=False) == region.peptide


# ---------------------------------------------------------------- end to end


def test_run_refuses_to_guess_when_the_data_is_not_installed(tmp_path):
    from hydra_amr.utils import HydraError
    with pytest.raises(HydraError, match="not installed"):
        spn.run(tmp_path / "x.fna", tmp_path / "db")


def test_gaps_reports_the_holes_in_a_locus():
    assert spn._gaps([(1, 100), (201, 300)], 300) == [(101, 200)]
    assert spn._gaps([(50, 100)], 200) == [(101, 200), (1, 49)]


# --------------------------------------------------------- validation (network)
#: Each entry: accession -> (strain, expected serotype, expected PBP triplet or
#: None, expected penicillin category or None, where the expectation comes
#: from). The reasoning behind each is in docs/typing/spneumoniae.md.
VALIDATION: dict[str, tuple[str, str, str | None, str | None]] = {
    "GCF_000006885.1": ("TIGR4", "4", "0:0:3", "S"),
    "GCF_000014365.1": ("D39", "2", None, None),
    "GCF_000007045.1": ("R6", "non-encapsulated", None, None),
    "GCF_000026665.1": ("ATCC 700669 (Spain23F-1)", "23F", None, "R"),
    "GCF_000147095.1": ("670-6B", "6B", None, "R"),
    "GCF_000019265.1": ("Hungary19A-6", "19A", None, None),
    "GCF_000019025.1": ("Taiwan19F-14", "19F", None, "R"),
    "GCF_000019985.1": ("CGSP14", "14", None, None),
    "GCF_000210955.1": ("OXC141", "3", None, "S"),
    "GCF_000019005.1": ("P1031", "1", None, "S"),
    "GCF_000018965.1": ("70585", "5", None, "S"),
    "GCF_033096205.1": ("PZ900700054", "6A", None, None),
    "GCF_033096125.1": ("PZ900700012", "6C", None, None),
    "GCF_022558605.1": ("PZ900700608", "15B/15C", None, None),
    "GCF_033096045.1": ("PZ900700005", "9A/9V", None, None),
    "GCF_033096325.1": ("PZ900701114", "24A/24B/24F", None, None),
    # Submitted to NCBI as the unresolved pair "35B/35D"; the CDC rule resolves
    # it to 35B on an exact wciG match. The narrowing is deliberate and is
    # discussed in docs/typing/spneumoniae.md; the test pins what the published
    # rule gives, not what the submitter could not decide.
    "GCF_036322935.1": ("PF2/TKY092604", "35B", None, None),
    "GCF_022069005.1": ("TVO_1901948", "8", None, None),
    "GCF_002813955.1": ("Korean 11A isolate", "11A", None, None),
    "GCF_000180515.1": ("SPNA45", "3", None, None),
}

validation = pytest.mark.skipif(
    os.environ.get("HYDRA_VALIDATION") != "1",
    reason="set HYDRA_VALIDATION=1 to download public genomes and validate")


@validation
@pytest.mark.parametrize("accession", sorted(VALIDATION))
def test_public_genomes_type_as_published(accession, tmp_path_factory):
    strain, serotype, triplet, penicillin = VALIDATION[accession]
    shared = tmp_path_factory.getbasetemp() / "spn-validation"
    db_dir = shared / "db"
    spn.fetch(db_dir)
    assembly = _download_assembly(accession, shared / "genomes")
    results = {r.scheme: r for r in spn.run(assembly, db_dir, threads=2,
                                            workdir=shared / accession)}

    assert results["cps_serotype"].call == serotype, \
        f"{strain}: expected {serotype}, got {results['cps_serotype'].call}"
    if triplet is not None:
        assert results["pbp_type"].call == triplet, strain
    if penicillin is not None:
        mic = results["beta_lactam_mic"]
        assert "penicillin" in mic.lineage, f"{strain}: no penicillin verdict"
        # The meningitis criterion is the strict one and is reported first.
        assert mic.lineage.split(":")[1].strip().startswith(penicillin), \
            f"{strain}: expected penicillin {penicillin}, got {mic.lineage}"
    # Whatever the call, it must carry evidence and an interpretation.
    for result in results.values():
        assert result.note, f"{strain}: {result.scheme} has no note"


@validation
def test_the_databases_build_and_record_their_provenance(tmp_path_factory):
    db_dir = tmp_path_factory.getbasetemp() / "spn-validation" / "db"
    spn.fetch(db_dir)
    sources = (db_dir / spn.KEY / "SOURCES.tsv").read_text(encoding="utf-8")
    for _url, path, _dbtype, _name in spn._FILES:
        assert path in sources
    # Idempotent: a second call must not fail or re-download.
    spn.fetch(db_dir)


_AGENT = {"User-Agent": "hydra-amr validation"}
_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _urlopen(url: str, timeout: int = 180) -> bytes:
    """GET *url*, backing off properly when NCBI throttles us.

    Without this the suite fails on HTTP 429 rather than on a typing result,
    which would look like a validation failure and is not one.
    """
    import time  # noqa: PLC0415
    import urllib.error  # noqa: PLC0415
    import urllib.request  # noqa: PLC0415

    last: Exception | None = None
    for attempt in range(6):
        if attempt:
            time.sleep(3.0 * attempt)
        try:
            request = urllib.request.Request(url, headers=_AGENT)
            with urllib.request.urlopen(request, timeout=timeout) as handle:
                return handle.read()
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code not in (429, 500, 502, 503, 504):
                break
        except (urllib.error.URLError, OSError) as exc:
            last = exc
    pytest.skip(f"could not fetch {url}: {last}")
    raise AssertionError("unreachable")


def _download_assembly(accession: str, into: Path) -> Path:
    """Fetch one assembly's FASTA from the NCBI genomes tree over HTTPS."""
    import gzip  # noqa: PLC0415
    import json  # noqa: PLC0415
    import time  # noqa: PLC0415
    import urllib.parse  # noqa: PLC0415

    into.mkdir(parents=True, exist_ok=True)
    out = into / f"{accession}.fna"
    if out.exists() and out.stat().st_size:
        return out
    query = urllib.parse.urlencode({"db": "assembly", "term": accession,
                                    "retmax": 1, "retmode": "json"})
    payload = json.loads(_urlopen(f"{_EUTILS}/esearch.fcgi?{query}"))
    uids = payload["esearchresult"]["idlist"]
    if not uids:
        pytest.skip(f"NCBI returned no assembly for {accession}")
    # Entrez allows three requests a second without a key, and the suite is
    # parametrised, so pace the second call rather than risk a 429.
    time.sleep(1.0)
    query = urllib.parse.urlencode({"db": "assembly", "id": uids[0],
                                    "retmode": "json"})
    record = json.loads(_urlopen(f"{_EUTILS}/esummary.fcgi?{query}"))
    record = record["result"][uids[0]]
    path = record.get("ftppath_refseq") or record.get("ftppath_genbank")
    if not path:
        pytest.skip(f"no FTP path published for {accession}")
    path = path.replace("ftp://", "https://")
    url = f"{path}/{path.rsplit('/', 1)[1]}_genomic.fna.gz"
    out.write_bytes(gzip.decompress(_urlopen(url, timeout=300)))
    return out
