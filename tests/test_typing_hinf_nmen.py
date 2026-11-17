"""Tests for the *Haemophilus influenzae* and *Neisseria meningitidis* typing modules.

The unit tests are fast, offline and run on synthetic fixtures. They exercise
the decision rules, which is where a typing module goes wrong: the thresholds
and the BLAST plumbing are Hydra's and are tested elsewhere.

The validation tests download public genomes with published types and assert the
calls. They only run with ``HYDRA_VALIDATION=1`` and need the module databases
installed, because they need the network and a few minutes.
"""

from __future__ import annotations

import os
import shutil
import urllib.request
import zipfile
from pathlib import Path

import pytest

from hydra_amr.typing.organisms import hinfluenzae as hi
from hydra_amr.typing.organisms import nmeningitidis as nm
from hydra_amr.typing.organisms import _hinf_nmen_common as common
from hydra_amr.typing.organisms.base import OrganismDb, Reference
from hydra_amr.utils import revcomp, translate


# ------------------------------------------------------------------ interface
@pytest.mark.parametrize("module", [hi, nm])
def test_module_exposes_the_interface(module):
    assert isinstance(module.KEY, str) and module.KEY
    assert isinstance(module.TITLE, str) and " " in module.TITLE
    assert isinstance(module.SPECIES, tuple) and module.SPECIES
    assert isinstance(module.MLST_SCHEMES, tuple)
    assert callable(module.fetch) and callable(module.run)
    assert all(isinstance(db, OrganismDb) for db in module.DATABASES)
    assert all(isinstance(ref, Reference) for ref in module.REFERENCES)


@pytest.mark.parametrize("module", [hi, nm])
def test_every_database_declares_its_licence_and_a_reference(module):
    keys = {ref.key for ref in module.REFERENCES}
    for db in module.DATABASES:
        assert db.licence and db.licence != ""
        assert db.redistribution in ("permitted", "restricted", "unstated")
        assert db.reference in keys, f"{db.name} cites unknown reference {db.reference}"
        assert db.version


@pytest.mark.parametrize("module", [hi, nm])
def test_every_reference_says_what_it_supports(module):
    for ref in module.REFERENCES:
        assert ref.supports, f"{ref.key} does not say what it supports"
        assert ref.doi or ref.pmid, f"{ref.key} has neither a DOI nor a PMID"
        assert 1950 < ref.year <= 2027


# ------------------------------------------------------- shared helpers
def test_best_reading_frame_prefers_the_real_coding_frame():
    # Both orientations of a simple gene can be free of internal stops, so the
    # helper must prefer the one that actually runs start codon to stop codon
    # rather than whichever translation comes out a residue longer.
    gene = "ATG" + "AAA" * 20 + "TAA"
    assert common.best_reading_frame(gene) == translate(gene).rstrip("*")


def test_best_reading_frame_finds_a_reference_stored_on_the_wrong_strand():
    # This is the real case: the hicap bex genes are stored as the plus-strand
    # slice of the genome they came from, so frame 1 is full of stop codons and
    # only the reverse complement translates cleanly.
    # TTA, CTA and TCA are codons whose reverse complements are stop codons, so
    # a gene carrying them in frame is unreadable on the opposite strand.
    gene = "ATG" + "TTAGCTCTAAAATCAGGT" * 2 + "TAA"
    stored = revcomp(gene)
    assert translate(stored).count("*") > 0, "fixture must be broken in frame 1"
    assert common.best_reading_frame(stored) == translate(gene).rstrip("*")


def test_best_reading_frame_rejects_a_sequence_with_no_clean_frame():
    # Stop codons in frame 1 of both orientations: TAA/TTA and TAG/CTA.
    both_broken = "ATGTAACTATTAGGGTAGCTATTATAA"
    assert translate(both_broken).count("*") > 0
    assert translate(revcomp(both_broken)).count("*") > 0
    assert common.best_reading_frame(both_broken) == ""


def test_substitutions_are_named_from_the_reference():
    reference = "MKDSAN"
    assert common.substitutions(reference, "MKDSAN") == []
    assert common.substitutions(reference, "MKNSAN") == ["D3N"]
    assert common.substitutions(reference, "MKNSAK") == ["D3N", "N6K"]


def test_substitutions_fall_back_to_the_alignment_when_lengths_differ():
    # One residue inserted in the query: positional comparison would report
    # every downstream residue as substituted, which is the bug this guards.
    found = common.substitutions(
        "MKDSAN", "MKXDSAN",
        aligned_query="MK-DSAN", aligned_subject="MKXDSAN", offset=0)
    assert found == []


def test_match_peptide_names_an_exact_definition_and_flags_a_new_one():
    index = {"ACDEFG": "1", "ACDEFH": "2"}
    assert common.match_peptide("ACDEFG", index) == ("1", "")
    identifier, note = common.match_peptide("ACDEFW", index)
    assert identifier == "new"
    assert "closest is" in note and "1 residue difference" in note
    assert common.match_peptide("", index)[0] == ""


def test_find_variable_region_takes_the_longest_contained_definition():
    index = {"KLSSTNAK": "7", "KLSSTNAKTG": "7-2"}
    protein = "MRAAKLSSTNAKTGNKVEV"
    assert common.find_variable_region(protein, index) == ("7-2", "")
    identifier, note = common.find_variable_region("MRAAQQQQ", index)
    assert identifier == "new"
    assert "new or the gene is incomplete" in note


def test_extract_cds_recovers_the_gene_on_both_strands(tmp_path):
    gene = "ATG" + "GCTAAAGGTTTCCATACC" * 6 + "TAA"
    protein = translate(gene).rstrip("*")
    contigs = {"c1": "TTTTTT" + gene + "GGGGGG",
               "c2": "GGGGGG" + revcomp(gene) + "TTTTTT"}
    forward = common.extract_cds(contigs, "c1", 7, 6 + len(gene) - 3, "+", ref_offset=0)
    assert forward is not None and forward.protein == protein and forward.complete
    reverse = common.extract_cds(contigs, "c2", 10, 6 + len(gene), "-", ref_offset=0)
    assert reverse is not None and reverse.protein == protein


def test_extract_cds_reports_an_internal_stop_rather_than_hiding_it():
    gene = "ATG" + "GCTAAAGGT" * 3 + "TAA" + "GCTAAAGGT" * 3 + "TAA"
    contigs = {"c1": "TTTTTT" + gene}
    cds = common.extract_cds(contigs, "c1", 7, 6 + len(gene), "+", ref_offset=0)
    assert cds is not None
    assert cds.internal_stops == 1


def test_profile_table_separates_loci_from_the_clonal_complex(tmp_path):
    path = tmp_path / "profiles.txt"
    path.write_text("ST\tabcZ\tadk\tclonal_complex\n"
                    "33\t8\t10\tST-32 complex\n"
                    "99\t1\t2\t\n")
    loci, table, complexes = common.profile_table(path)
    assert loci == ["abcZ", "adk"]
    assert table[("8", "10")] == "33"
    assert complexes == {"33": "ST-32 complex"}


# --------------------------------------------------- H. influenzae rules
def test_cap_regions_are_assigned_correctly():
    assert hi._region_of("bexA") == ("I", "")
    assert hi._region_of("hcsB") == ("III", "")
    assert hi._region_of("bcs3") == ("II", "b")
    assert hi._region_of("ecs8") == ("II", "e")
    assert hi._region_of("IS1016V6") == ("?", "")


def test_region_two_gene_sets_match_the_published_cassettes():
    # Counts from the papers that named each cassette: four acs, four bcs,
    # four ccs, five dcs, eight ecs, three fcs.
    assert [len(hi.REGION_TWO[s]) for s in "abcdef"] == [4, 4, 4, 5, 8, 3]
    all_genes = [g for genes in hi.REGION_TWO.values() for g in genes]
    assert len(all_genes) == len(set(all_genes))


def _cap_gene(gene, *, serotype="", contig="c1", start=1, end=1000, broken=False):
    region, inferred = hi._region_of(gene)
    return hi.CapGene(gene=gene, serotype=serotype or inferred, region=region,
                      contig=contig, start=start, end=end, strand="+",
                      identity=99.0, coverage=60.0 if broken else 100.0,
                      covered_bases=600 if broken else 1000, score=99.0, broken=broken)


def test_serotype_is_scored_on_the_most_complete_region_two_set():
    # A serotype b locus also matches two serotype a genes, because acs1/bcs1
    # and acs2/bcs2 are 85-95% identical. The complete set must win.
    cluster = [_cap_gene(g) for g in hi.REGION_TWO["b"]]
    cluster += [_cap_gene("acs1", start=2000, end=3000)]
    serotype, _score, per_type = hi._serotype_score(cluster)
    assert serotype == "b"
    assert len(per_type["b"]) == 4


def test_duplication_needs_more_than_one_duplicated_gene():
    single = [_cap_gene("bcs1"), _cap_gene("bcs1", start=20000, end=21000)]
    duplicated, description, note = hi._copy_number(single)
    assert duplicated is False
    assert "too little to call" in note

    many = single + [_cap_gene("bcs2"), _cap_gene("bcs2", start=22000, end=23000)]
    duplicated, description, note = hi._copy_number(many)
    assert duplicated is True
    assert "at least 2 copies" in description


def test_copy_number_ignores_bexa_because_it_is_the_deleted_gene():
    genes = [_cap_gene("bexA"), _cap_gene("bexA", start=20000, end=21000)]
    duplicated, _description, _note = hi._copy_number(genes)
    assert duplicated is False


@pytest.mark.parametrize("found,stage,group", [
    (set(), 0, ""),
    ({"D350N", "V547I"}, 0, ""),
    ({"N526K"}, 1, "II"),
    ({"R517H"}, 1, "I"),
    ({"Y528H"}, 1, "Misc1"),
    ({"N526H"}, 1, "Misc2"),
    ({"N526K", "S385T"}, 2, "III(-)"),
    ({"R517H", "S385T"}, 2, "III-like(-)"),
    ({"N526K", "S385T", "L389F"}, 3, "III(+)"),
    ({"R517H", "S385T", "L389F"}, 3, "III-like(+)"),
    ({"N526K", "S385T", "L389F", "M377I", "D350N"}, 3, "III(+)"),
])
def test_pbp3_stage_follows_the_published_decision_tree(found, stage, group):
    assert hi._pbp3_group(found) == (stage, group)


def test_pbp3_group_two_subgroups_follow_the_residue_at_502_and_449():
    assert hi._pbp3_subgroup("II", {"N526K"}) == "IIa"
    assert hi._pbp3_subgroup("II", {"N526K", "A502V"}) == "IIb"
    assert hi._pbp3_subgroup("II", {"N526K", "A502T"}) == "IIc"
    assert hi._pbp3_subgroup("II", {"N526K", "I449V"}) == "IId"
    assert hi._pbp3_subgroup("I", {"R517H"}) == ""


def test_named_substitutions_keeps_only_scheme_positions_with_the_right_reference():
    subs = ["D350N", "N526K", "Q100R", "A350N"]
    named = hi._named_substitutions(subs)
    assert named == ["D350N", "N526K"]


def test_pbp3_positions_are_unique_and_cover_the_essential_residues():
    positions = [position for position, _residue in hi.PBP3_POSITIONS]
    assert len(positions) == len(set(positions))
    for essential in ("R517H", "N526K", "N526H", "Y528H"):
        position = int("".join(c for c in essential if c.isdigit()))
        assert position in positions


# ------------------------------------------------- N. meningitidis rules
def test_capsule_gene_table_is_consistent():
    names = set(nm.CAPSULE_LOCI.values())
    for group, essential in nm.SEROGROUP_ESSENTIAL.items():
        for gene in essential:
            assert gene in names, f"{group} needs {gene}, which has no PubMLST locus"
    for gene in nm.TRANSPORT_ESSENTIAL + ("tex", "cnl"):
        assert gene in names
    for gene, group in nm.SEROGROUP_UNIQUE.items():
        assert gene in names
        assert group in nm.SEROGROUP_ESSENTIAL or group in ("E", "L", "X", "Z", "H")


def test_b_c_w_and_y_share_the_sialic_acid_genes_and_differ_only_by_polymerase():
    shared = {"cssA", "cssB", "cssC"}
    for group in ("B", "C", "W", "Y"):
        assert shared <= set(nm.SEROGROUP_ESSENTIAL[group])
    polymerases = {group: set(nm.SEROGROUP_ESSENTIAL[group]) - shared
                   for group in ("B", "C", "W", "Y")}
    assert polymerases == {"B": {"csb"}, "C": {"csc"}, "W": {"csw"}, "Y": {"csy"}}
    for gene, group in (("csb", "B"), ("csc", "C"), ("csw", "W"), ("csy", "Y")):
        assert nm.SEROGROUP_UNIQUE[gene] == group


def test_bast_has_five_components_in_the_scheme_order():
    assert nm.BAST_ORDER == ("fHbp_peptide", "NHBA_peptide", "NadA_peptide",
                             "PorA_VR1", "PorA_VR2")
    assert set(nm.BAST_PEPTIDE_LOCI) <= set(nm.BAST_ORDER)
    assert "PorA_VR1" in nm.BAST_ORDER, "PorA VR1 is part of BAST despite not being in MATS"


def test_pubmlst_locus_names_are_case_sensitive_as_spelled():
    # The lower-case spellings answer HTTP 404, so a typo here would silently
    # leave the fine type undetermined.
    assert nm.VR_LOCI == ("PorA_VR1", "PorA_VR2", "FetA_VR")
    assert nm.BAST_PEPTIDE_LOCI == ("fHbp_peptide", "NHBA_peptide", "NadA_peptide")


def test_clonal_complex_is_rendered_as_the_designation_spells_it():
    assert nm._short_complex("ST-32 complex") == "cc32"
    assert nm._short_complex("ST-41/44 complex") == "cc41/44"
    assert nm._short_complex("") == ""
    assert nm._short_complex("-") == ""


def test_mendevar_categories_are_explained_not_just_repeated():
    for value in ("exact match", "cross-reactive", "none", "insufficient data"):
        rendered = nm._mendevar(value, "Bexsero")
        assert value in rendered and "Bexsero" in rendered
        assert len(rendered) > len(value) + 20, "the category must be explained"


def test_finetype_string_follows_the_published_nomenclature():
    from hydra_amr.records import TypingResult

    capsule = TypingResult(scheme="cps", call="B", lineage="serogroup B")
    mlst = TypingResult(scheme="mlst", call="ST-33", lineage="ST-32 complex")
    peptides = {"PorA_VR1": ("7-2", ""), "PorA_VR2": ("4", ""), "FetA_VR": ("F1-5", "")}
    result = nm._call_finetype(capsule, peptides, mlst)
    assert result.call == "B:P1.7-2,4:F1-5"
    assert result.alleles["full_designation"] == "B:P1.7-2,4:F1-5:ST-33(cc32)"
    assert result.loci_found == 3


def test_finetype_marks_a_missing_region_rather_than_inventing_one():
    from hydra_amr.records import TypingResult

    capsule = TypingResult(scheme="cps", call="NG", lineage="genogroup B")
    mlst = TypingResult(scheme="mlst", call="not determined")
    peptides = {"PorA_VR1": ("", "not determined: PorA not recovered"),
                "PorA_VR2": ("4", ""), "FetA_VR": ("", "not determined")}
    result = nm._call_finetype(capsule, peptides, mlst)
    assert result.call == "NG:P1.?,4:?"
    assert "not determined" in result.note
    assert result.loci_found == 1


# ------------------------------------------------------------- validation
#: Accession, strain, expected capsule call, and where the expectation comes
#: from. Every row was checked against hicap 1.0.4 run locally, and the
#: published serotype where one exists.
HINF_VALIDATION = (
    ("GCF_000210875.1", "10810", "b", "Watts2019 Table 2 (serologically confirmed b); "
                                      "hicap predicts type_b"),
    ("GCF_900478275.1", "NCTC 13377", "b", "Watts2019 Table 2; hicap predicts type_b"),
    ("GCF_001856725.1", "NML-Hia-1", "a", "Watts2019 Table 2; hicap predicts type_a"),
    ("GCF_003351605.1", "M12125", "c", "Watts2019 Table 2; hicap predicts type_c"),
    ("GCF_900407385.1", "PTHi-10983", "d", "Watts2019 Table S1, the only cap-d genome in "
                                           "the survey; hicap predicts type_d"),
    ("GCF_001975845.1", "hi467", "e", "Watts2019 Table 2; hicap predicts type_e"),
    ("GCF_000465255.1", "KR494", "f", "Watts2019 Table 2; hicap predicts type_f"),
    ("GCF_900475755.1", "NCTC 11426", "f", "Watts2019 Table 2; hicap predicts type_f"),
    ("GCF_000027305.1", "Rd KW20", "NT", "Watts2019 Table 2 'No cap locus'; hicap writes "
                                         "no report. Retains a full-length IS1016, so it "
                                         "is also the test that IS1016 is not treated as a "
                                         "capsule marker"),
    ("GCF_000012185.1", "86-028NP", "NT", "Watts2019 Table 2 'No cap locus'"),
    ("GCF_000016465.1", "PittEE", "NT", "Watts2019 Table 2 'No cap locus'"),
    ("GCF_000165525.1", "R2866", "NT", "Watts2019 Table S1 (predicted serotype '-')"),
)

#: Accession, strain, expected serogroup call, expected genogroup, source.
NMEN_VALIDATION = (
    ("GCF_000009465.1", "FAM18", "C", "C", "Bentley 2007 genome paper and the NCBI "
                                           "BioSample serotype field: serogroup C"),
    ("GCF_000026965.1", "8013", "C", "C", "Rusniok 2009 genome paper: serogroup C, "
                                          "ST-177, cc18"),
    ("GCF_000191445.1", "H44/76", "B", "B", "the reference serogroup B strain, published "
                                            "fine type B:P1.7,16:ST-32(cc32)"),
    ("GCF_000191525.1", "NZ-05/33", "B", "B", "the New Zealand epidemic strain, published "
                                              "fine type B:P1.7-2,4, cc41/44"),
    ("GCF_000191465.1", "M01-240149", "B", "B", "serogroup B disease isolate from the "
                                                "Broad MenB project"),
    ("GCF_000191485.1", "M01-240355", "B", "B", "serogroup B disease isolate, ST-213"),
    ("GCF_000152165.1", "alpha710", "B", "B", "Schoen 2008; PubMLST and the BioSample "
                                              "record both give serogroup B"),
    ("GCF_000083565.1", "alpha14", "cnl", "cnl", "Schoen 2008 PNAS; the ENA sample record "
                                                 "gives serotype cnl. The capsule-null "
                                                 "case"),
    ("GCF_000009105.1", "Z2491", "A", "A", "Parkhill 2000 Nature: serogroup A"),
    ("GCF_000008805.1", "MC58", "B", "B", "Tettelin 2000 Science; BioSample "
                                          "SAMN02603997 serogroup=B"),
)

_DOWNLOAD = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
             "{accession}/download?include_annotation_type=GENOME_FASTA")

validation = pytest.mark.skipif(
    os.environ.get("HYDRA_VALIDATION") != "1",
    reason="set HYDRA_VALIDATION=1 to run the validation tests (network, several minutes)")


def _db_dir() -> Path:
    root = os.environ.get("HYDRA_DB_DIR")
    if not root:
        pytest.skip("set HYDRA_DB_DIR to the directory holding the installed typing data")
    return Path(root)


def _genome(accession: str, cache: Path) -> Path:
    target = cache / f"{accession}.fna"
    if target.exists():
        return target
    archive = cache / f"{accession}.zip"
    request = urllib.request.Request(_DOWNLOAD.format(accession=accession),
                                     headers={"User-Agent": common.USER_AGENT})
    with urllib.request.urlopen(request, timeout=300) as response:
        archive.write_bytes(response.read())
    with zipfile.ZipFile(archive) as handle:
        names = [name for name in handle.namelist() if name.endswith(".fna")]
        assert len(names) == 1, f"{accession}: unexpected archive contents {names}"
        target.write_bytes(handle.read(names[0]))
    archive.unlink()
    return target


@validation
@pytest.mark.parametrize("accession,strain,expected,source", HINF_VALIDATION,
                         ids=[row[1] for row in HINF_VALIDATION])
def test_validate_hinfluenzae_capsule(accession, strain, expected, source, tmp_path_factory):
    cache = tmp_path_factory.mktemp("hinf", numbered=False)
    db_dir = _db_dir()
    if not (db_dir / hi.KEY / "cap" / "cap_genes.fna").exists():
        pytest.skip(f"{hi.KEY} data not installed in {db_dir}")
    genome = _genome(accession, cache)
    results = {result.scheme: result for result in hi.run(genome, db_dir, threads=2)}
    assert results["cap"].call == expected, (
        f"{strain} ({accession}): expected {expected} from {source}, "
        f"got {results['cap'].call} - {results['cap'].note}")
    # Every call must carry its evidence; a bare serotype letter is what the
    # owner rejected.
    assert results["cap"].alleles
    assert results["cap"].note
    shutil.rmtree(cache, ignore_errors=True)


@validation
def test_validate_hinfluenzae_hib_duplication_is_not_called_capsule_deficient(
        tmp_path_factory):
    """Strain 10810 carries a duplicated cap-b locus with one bexA truncated.

    That is the normal, virulent Hib genotype, so it must be reported as
    encapsulated. Calling it capsule-deficient would be a clinical error.
    """
    cache = tmp_path_factory.mktemp("hinf_dup", numbered=False)
    db_dir = _db_dir()
    if not (db_dir / hi.KEY / "cap" / "cap_genes.fna").exists():
        pytest.skip(f"{hi.KEY} data not installed in {db_dir}")
    result = {r.scheme: r for r in hi.run(_genome("GCF_000210875.1", cache), db_dir,
                                          threads=2)}["cap"]
    assert result.call == "b"
    assert "encapsulated" in result.lineage
    assert "IS1016-bexA" in result.lineage or "IS1016-bexA" in result.note
    shutil.rmtree(cache, ignore_errors=True)


@validation
@pytest.mark.parametrize("accession,strain,expected,genogroup,source", NMEN_VALIDATION,
                         ids=[row[1] for row in NMEN_VALIDATION])
def test_validate_nmeningitidis_serogroup(accession, strain, expected, genogroup, source,
                                          tmp_path_factory):
    cache = tmp_path_factory.mktemp("nmen", numbered=False)
    db_dir = _db_dir()
    if not (db_dir / nm.KEY / "capsule" / "cps_proteins.faa").exists():
        pytest.skip(f"{nm.KEY} data not installed in {db_dir}")
    genome = _genome(accession, cache)
    results = {result.scheme: result for result in nm.run(genome, db_dir, threads=2)}
    capsule = results["cps"]
    # A non-groupable call is acceptable only if the genogroup is still right:
    # the published serogroup is a serological result and an assembly can show a
    # lesion the serology of the sequenced stock did not have.
    assert capsule.call == expected or (
        capsule.call == "NG" and capsule.alleles.get("genogroup") == genogroup), (
        f"{strain} ({accession}): expected {expected} from {source}, got "
        f"{capsule.call} (genogroup {capsule.alleles.get('genogroup')}) - {capsule.note}")
    assert capsule.alleles and capsule.note
    shutil.rmtree(cache, ignore_errors=True)


@validation
def test_validate_nmeningitidis_published_fine_types(tmp_path_factory):
    """Two strains whose full designation is published, checked field by field."""
    cache = tmp_path_factory.mktemp("nmen_fine", numbered=False)
    db_dir = _db_dir()
    if not (db_dir / nm.KEY / "peptides" / "PorA_VR1.fas").exists():
        pytest.skip(f"{nm.KEY} data not installed in {db_dir}")
    expectations = {
        # H44/76 is the reference serogroup B strain, published as
        # B:P1.7,16:ST-32(cc32).
        "GCF_000191445.1": {"PorA_VR1": "7", "PorA_VR2": "16", "ST": "ST-32",
                            "cc": "ST-32 complex"},
        # NZ-05/33 is the New Zealand epidemic strain, B:P1.7-2,4, cc41/44.
        "GCF_000191525.1": {"PorA_VR1": "7-2", "PorA_VR2": "4", "cc": "ST-41/44 complex"},
    }
    for accession, expected in expectations.items():
        results = {r.scheme: r for r in nm.run(_genome(accession, cache), db_dir,
                                               threads=2)}
        fine, mlst = results["finetype"], results["mlst"]
        assert fine.alleles["PorA_VR1"] == expected["PorA_VR1"], (
            f"{accession}: PorA VR1 {fine.alleles['PorA_VR1']}, "
            f"expected {expected['PorA_VR1']}")
        assert fine.alleles["PorA_VR2"] == expected["PorA_VR2"]
        if "ST" in expected:
            assert mlst.call == expected["ST"]
        assert mlst.lineage == expected["cc"]
    shutil.rmtree(cache, ignore_errors=True)
