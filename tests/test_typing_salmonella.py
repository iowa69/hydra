"""Tests for the *Salmonella enterica* typing module.

The fast tests build a miniature database and a synthetic assembly from it, so
they need BLAST+ on PATH but no network and no real reference data. The
validation tests download public genomes and assert the published serovar of
each; they only run with ``HYDRA_VALIDATION=1`` set.
"""

from __future__ import annotations

import gzip
import json
import os
import random
import shutil
import urllib.request
from pathlib import Path

import pytest

from hydra_amr.seqio import write_fasta
from hydra_amr.typing.organisms import base
from hydra_amr.typing.organisms import salmonella as sal
from hydra_amr.utils import HydraError

HAVE_BLAST = shutil.which("blastn") is not None and shutil.which("makeblastdb") is not None
blast_only = pytest.mark.skipif(not HAVE_BLAST, reason="BLAST+ is not on PATH")
validation_only = pytest.mark.skipif(
    os.environ.get("HYDRA_VALIDATION") != "1",
    reason="set HYDRA_VALIDATION=1 to run the download-backed validation tests")


# --------------------------------------------------------------------- shape
def test_module_exposes_the_agreed_interface():
    assert sal.KEY == "salmonella"
    assert sal.TITLE == "Salmonella enterica"
    assert "Salmonella enterica" in sal.SPECIES
    assert "senterica_achtman_2" in sal.MLST_SCHEMES
    assert all(isinstance(d, base.OrganismDb) for d in sal.DATABASES)
    assert all(isinstance(r, base.Reference) for r in sal.REFERENCES)
    assert callable(sal.fetch) and callable(sal.run)


def test_every_database_records_its_licence_and_a_reference():
    keys = {r.key for r in sal.REFERENCES}
    for db in sal.DATABASES:
        assert db.source and db.version and db.licence
        assert db.redistribution in ("permitted", "restricted", "unstated")
        assert db.reference in keys, f"{db.name} cites an unknown reference"


def test_every_reference_says_what_it_supports_and_is_identifiable():
    for ref in sal.REFERENCES:
        assert ref.supports, f"{ref.key} does not say what it supports"
        assert ref.authors and ref.title and ref.journal and ref.year
        # A monograph has neither; everything else must be resolvable.
        if ref.key != "Grimont2007":
            assert ref.doi or ref.pmid, f"{ref.key} has no DOI and no PMID"


# ------------------------------------------------------------ pure functions
def test_clean_id_keeps_antigen_labels_blast_safe():
    assert sal._clean_id("g,m,[s],t") == "g,m,[s],t"
    assert sal._clean_id("e,n,x, z15") == "e,n,x,z15"
    assert sal._clean_id("a|b") == "a/b"


def test_similar_o_groups_uses_the_tables_own_spelling():
    # SISTR writes this group as ("S", "62"), but the WHO 2007 table spells the
    # second one "O:62", so a literal copy would never fire.
    assert sal._similar_o_groups("S") == ("O:62",)
    assert sal._similar_o_groups("O:62") == ("S",)
    assert set(sal._similar_o_groups("D1")) == {"A", "D2"}
    assert sal._similar_o_groups("B") == ()
    assert sal._similar_o_groups(None) == ()


def _rows():
    return [
        sal.SerovarRow("Typhimurium", "1,4,[5],12", "i", "1,2", "", "B", False,
                       "enterica"),
        sal.SerovarRow("I 1,4,[5],12:i:-", "1,4,[5],12", "i", "-", "", "B", True,
                       "enterica"),
        sal.SerovarRow("Heidelberg", "1,4,[5],12", "r", "1,2", "", "B", False,
                       "enterica"),
        sal.SerovarRow("Enteritidis", "1,9,12", "g,m", "-", "", "D1", True,
                       "enterica"),
        sal.SerovarRow("Gallinarum", "1,9,12", "-", "-", "", "D1", True,
                       "enterica"),
        sal.SerovarRow("II 1,9,12:g,m:e,n,x", "1,9,12", "g,m", "e,n,x", "", "D1",
                       False, "salamae"),
    ]


def test_match_serovars_on_a_complete_formula():
    assert [r.serovar for r in sal._match_serovars(_rows(), "B", "i", "1,2")] \
        == ["Typhimurium"]


def test_match_serovars_treats_absent_h2_differently_from_undetermined():
    rows = _rows()
    # A genuine negative only fits a monophasic row or one the scheme allows to
    # lose H2 -- it must not pull in biphasic Typhimurium.
    absent = [r.serovar for r in sal._match_serovars(rows, "B", "i", "-")]
    assert absent == ["I 1,4,[5],12:i:-"]
    # Undetermined fits anything, which is why the module says so out loud.
    undetermined = [r.serovar for r in sal._match_serovars(rows, "B", "i", None)]
    assert undetermined == ["Typhimurium", "I 1,4,[5],12:i:-"]


def test_match_serovars_separates_subspecies_by_formula_alone():
    rows = sal._match_serovars(_rows(), "D1", "g,m", "e,n,x")
    assert [r.subspecies for r in rows] == ["salamae"]


def test_non_motile_siblings_are_findable_so_they_can_be_caveated():
    # Gallinarum and Pullorum are D1:-:- and carry Enteritidis' flagellin
    # genes; the module must be able to name them beside an Enteritidis call.
    assert [r.serovar for r in sal._match_serovars(_rows(), "D1", "-", "-")] \
        == ["Gallinarum"]


def test_formula_prefixes_the_subspecies_numeral():
    assert sal._formula("1,4,[5],12", "i", "1,2", "enterica") == "I 1,4,[5],12:i:1,2"
    assert sal._formula("62", "z4,z23", "-", "arizonae") == "IIIa 62:z4,z23:-"
    assert sal._formula("1,9,12", "g,m", "-") == "1,9,12:g,m:-"


def test_o_antigen_for_picks_the_groups_usual_string():
    assert sal._o_antigen_for(_rows(), "B") == "1,4,[5],12"


def test_island_key_sorts_spi_numbers_numerically():
    names = ["SPI-10", "SPI-2", "CS54_island", "SPI-1"]
    assert sorted(names, key=sal._island_key) == \
        ["SPI-1", "SPI-2", "SPI-10", "CS54_island"]


def test_runners_up_reports_the_margin_not_just_the_name():
    def ev(label, identity):
        return sal.Evidence(label, label, "ref", identity, 100.0, 1500, "c1", 1,
                            1500, "+", "1-1500/1500", 5000, 1)
    best = ev("1,2", 100.0)
    text = sal._runners_up([best, ev("1,5", 99.3), ev("1,5", 98.0)], best)
    assert len(text) == 1 and "1,5" in text[0] and "99.3" in text[0]


def test_built_rejects_a_fasta_without_a_blast_index(tmp_path):
    fasta = tmp_path / "x.fna"
    fasta.write_text(">a\nACGT\n")
    assert not sal._built(fasta)
    (tmp_path / "x.fna.nin").write_text("")
    assert sal._built(fasta)


def test_evidence_spots_a_hit_sitting_on_a_contig_boundary():
    mid = sal.Evidence("i", "i", "r", 100.0, 100.0, 1500, "c", 2000, 3500, "+",
                       "1-1500/1500", 9000, 1)
    edge = sal.Evidence("i", "i", "r", 100.0, 60.0, 900, "c", 40, 940, "+",
                        "1-900/1500", 1000, 1)
    assert not mid.at_contig_edge
    assert edge.at_contig_edge
    assert "at contig edge" in edge.describe()


def test_load_serovars_strips_stray_whitespace(tmp_path):
    path = tmp_path / "serovars.csv"
    path.write_text(
        "Serovar,O_antigen,H1,H2,H_other,Serogroup,can_h2_be_missing,subspecies\n"
        "Typhimurium,\"1,4,[5],12\",i,\"1,2\",,B,FALSE,enterica\n"
        "Blegdam,\"9,12\",\"g,m,q\",-,,\"R \",TRUE,enterica\n")
    rows = sal._load_serovars(path)
    assert rows[1].serogroup == "R"
    assert rows[1].h2_may_be_absent is True
    assert rows[0].h2_may_be_absent is False


def test_load_serovars_refuses_an_empty_table(tmp_path):
    path = tmp_path / "empty.csv"
    path.write_text("Serovar,O_antigen,H1,H2,H_other,Serogroup,"
                    "can_h2_be_missing,subspecies\n")
    with pytest.raises(HydraError):
        sal._load_serovars(path)


def test_run_explains_itself_when_the_database_is_not_installed(tmp_path):
    with pytest.raises(HydraError) as excinfo:
        sal.run(tmp_path / "nothing.fna", tmp_path / "db")
    assert "fetch()" in str(excinfo.value)


# ------------------------------------------------- synthetic end-to-end tests
_RANDOM = random.Random(20261117)


def _filler(length: int) -> str:
    return "".join(_RANDOM.choice("ACGT") for _ in range(length))


def _gene(length: int) -> str:
    """A pseudo-gene long enough to clear the module's length thresholds."""
    return _filler(length)


@pytest.fixture(scope="module")
def mini_db(tmp_path_factory):
    """A miniature Salmonella database: four antigen alleles and the markers.

    The sequences are random, so nothing here tests biology. What it tests is
    the decision logic end to end: that a formula becomes a serovar, that an
    absent fljAB region becomes a monophasic call and that a disrupted one
    becomes "not determined" rather than a negative.
    """
    if not HAVE_BLAST:
        pytest.skip("BLAST+ is not on PATH")
    root = tmp_path_factory.mktemp("salmdb") / sal.KEY
    sequences = {
        "wzx_B": _gene(1293), "wzy_B": _gene(1221),
        "fliC_i": _gene(1488), "fljB_1_2": _gene(1521),
        "fliC_r": _gene(1488),
        "fljA": _gene(540), "hin": _gene(573), "iroB": _gene(1116),
        "fljB_marker": None, "invA": _gene(2058), "oafA": _gene(1830),
        "tviB": _gene(1278), "vexC": _gene(696), "cdtB": _gene(810),
        "pltA": _gene(729), "pltB": _gene(414), "sdf-I": _gene(333),
        "spvR": _gene(894), "spvA": _gene(768), "spvB": _gene(1776),
        "spvC": _gene(726), "spvD": _gene(651),
    }
    # The fljB marker gene and the fljB antigen allele are the same locus.
    sequences["fljB_marker"] = sequences["fljB_1_2"]

    (root / "antigens").mkdir(parents=True)
    write_fasta(root / "antigens" / "antigens.fna", [
        ("wzx|1|B", sequences["wzx_B"]),
        ("wzy|1|B", sequences["wzy_B"]),
        ("fliC|1|i", sequences["fliC_i"]),
        ("fliC|2|r", sequences["fliC_r"]),
        ("fljB|1|1,2", sequences["fljB_1_2"]),
    ], wrap=0)
    (root / "antigens" / "serovars.csv").write_text(
        "Serovar,O_antigen,H1,H2,H_other,Serogroup,can_h2_be_missing,subspecies\n"
        "Typhimurium,\"1,4,[5],12\",i,\"1,2\",,B,FALSE,enterica\n"
        "\"I 1,4,[5],12:i:-\",\"1,4,[5],12\",i,-,,B,TRUE,enterica\n"
        "Heidelberg,\"1,4,[5],12\",r,\"1,2\",,B,FALSE,enterica\n")
    (root / "markers").mkdir(parents=True)
    write_fasta(root / "markers" / "markers.fna", [
        ("fliC|h1_locus|synthetic", sequences["fliC_i"]),
        ("fljB|h2_locus|synthetic", sequences["fljB_marker"]),
        ("fljA|h2_locus|synthetic", sequences["fljA"]),
        ("hin|h2_locus|synthetic", sequences["hin"]),
        ("iroB|h2_locus_flank|synthetic", sequences["iroB"]),
        ("invA|genus|synthetic", sequences["invA"]),
        ("oafA|o_antigen|synthetic", sequences["oafA"]),
        ("sdf-I|enteritidis|synthetic", sequences["sdf-I"]),
        ("tviB|viaB|synthetic", sequences["tviB"]),
        ("vexC|viaB|synthetic", sequences["vexC"]),
        ("cdtB|typhoid_toxin|synthetic", sequences["cdtB"]),
        ("pltA|typhoid_toxin|synthetic", sequences["pltA"]),
        ("pltB|typhoid_toxin|synthetic", sequences["pltB"]),
        ("spvR|spv|synthetic", sequences["spvR"]),
        ("spvA|spv|synthetic", sequences["spvA"]),
        ("spvB|spv|synthetic", sequences["spvB"]),
        ("spvC|spv|synthetic", sequences["spvC"]),
        ("spvD|spv|synthetic", sequences["spvD"]),
    ], wrap=0)
    for fasta in (root / "antigens" / "antigens.fna", root / "markers" / "markers.fna"):
        sal._makeblastdb(fasta)
    return root.parent, sequences


def _assembly(path: Path, parts, pad: int = 800) -> Path:
    """Lay the given sequences out on one contig with random spacers."""
    contig = _filler(pad)
    for piece in parts:
        contig += piece + _filler(pad)
    write_fasta(path, [("contig_1", contig)], wrap=70)
    return path


def _by_scheme(results):
    return {r.scheme: r for r in results}


@blast_only
def test_biphasic_assembly_is_named_from_its_formula(mini_db, tmp_path):
    db_dir, seq = mini_db
    path = _assembly(tmp_path / "biphasic.fna", [
        seq["wzx_B"], seq["wzy_B"], seq["fliC_i"], seq["invA"], seq["oafA"],
        seq["fljA"], seq["fljB_1_2"], seq["hin"], seq["iroB"],
    ])
    results = _by_scheme(sal.run(path, db_dir, workdir=tmp_path / "work"))
    serovar = results["salmonella_serovar"]
    assert serovar.call == "Typhimurium"
    assert serovar.alleles["antigenic_formula"] == "I 1,4,[5],12:i:1,2"
    assert serovar.alleles["O_serogroup"] == "serogroup B"
    assert serovar.loci_found == 4
    assert results["salmonella_subspecies"].call == "I (enterica)"
    assert results["salmonella_h2_locus"].call == "intact"
    # Evidence, not just a verdict: the allele, its identity and where it sits.
    assert "100.0% id" in serovar.alleles["fliC"]
    assert "contig_1" in serovar.alleles["fliC"]


@blast_only
def test_deleted_fljab_region_is_called_a_monophasic_variant(mini_db, tmp_path):
    db_dir, seq = mini_db
    path = _assembly(tmp_path / "monophasic.fna", [
        seq["wzx_B"], seq["wzy_B"], seq["fliC_i"], seq["invA"], seq["oafA"],
        seq["iroB"],
    ])
    results = _by_scheme(sal.run(path, db_dir, workdir=tmp_path / "work"))
    serovar = results["salmonella_serovar"]
    assert results["salmonella_h2_locus"].call == "fljBA-hin absent"
    assert serovar.alleles["H2_fljB"] == "-"
    assert serovar.alleles["antigenic_formula"] == "I 1,4,[5],12:i:-"
    assert serovar.call == "I 1,4,[5],12:i:-"
    assert "Typhimurium (H2 1,2)" in serovar.alleles["biphasic_counterparts"]
    assert "monophasic variant" in serovar.note


@blast_only
def test_fljb_lost_while_flja_and_hin_remain_is_not_reported_as_negative(
        mini_db, tmp_path):
    """The distinction the reference multiplex exists to make.

    fljA and hin are present, so the operon is there and only fljB is missing.
    That is a disrupted gene or an assembly gap, not a monophasic serovar, and
    H2 must come back undetermined rather than "-".
    """
    db_dir, seq = mini_db
    path = _assembly(tmp_path / "disrupted.fna", [
        seq["wzx_B"], seq["wzy_B"], seq["fliC_i"], seq["invA"],
        seq["fljA"], seq["hin"], seq["iroB"],
    ])
    results = _by_scheme(sal.run(path, db_dir, workdir=tmp_path / "work"))
    serovar = results["salmonella_serovar"]
    assert results["salmonella_h2_locus"].call == "partial"
    assert serovar.alleles["H2_fljB"] == "not determined"
    assert "not determined rather than negative" in serovar.note
    assert serovar.alleles["antigenic_formula"].endswith(":i:?")
    # Both candidates consistent with any second phase are offered.
    assert "Typhimurium" in serovar.call and "1,4,[5],12:i:-" in serovar.call


@blast_only
def test_fljb_truncated_at_a_contig_break_is_flagged(mini_db, tmp_path):
    """A gene cut by the end of a contig must not silently become a negative."""
    db_dir, seq = mini_db
    write_fasta(tmp_path / "broken.fna", [
        ("contig_1", _filler(400) + seq["wzx_B"] + _filler(400) + seq["fliC_i"]
         + _filler(400) + seq["invA"] + _filler(400) + seq["fljA"]
         + _filler(200) + seq["fljB_1_2"][:700]),
    ], wrap=70)
    results = _by_scheme(sal.run(tmp_path / "broken.fna", db_dir,
                                 workdir=tmp_path / "work"))
    h2 = results["salmonella_h2_locus"]
    assert h2.call == "partial"
    assert "contig boundary" in h2.note
    # The antigen is still reported - throwing it away would lose information -
    # but it is labelled provisional, because the conserved end of a flagellin
    # is shared between antigen types and half an allele can belong to another.
    reported = results["salmonella_serovar"].alleles["H2_fljB"]
    assert reported.startswith("1,2 (provisional")
    assert "contig boundary" in reported


@blast_only
def test_the_whole_spv_operon_is_reported_with_its_meaning(mini_db, tmp_path):
    db_dir, seq = mini_db
    path = _assembly(tmp_path / "pslt.fna", [
        seq["wzx_B"], seq["fliC_i"], seq["fljB_1_2"], seq["fljA"], seq["hin"],
        seq["spvR"], seq["spvA"], seq["spvB"], seq["spvC"], seq["spvD"],
    ])
    results = _by_scheme(sal.run(path, db_dir, workdir=tmp_path / "work"))
    plasmid = results["salmonella_virulence_plasmid"]
    assert plasmid.call == "spv locus complete"
    assert plasmid.loci_found == 5 and plasmid.loci_total == 5
    assert "ADP-ribosyltransferase" in plasmid.note


@blast_only
def test_absent_spv_is_reported_as_an_expected_negative(mini_db, tmp_path):
    db_dir, seq = mini_db
    path = _assembly(tmp_path / "nospv.fna",
                     [seq["wzx_B"], seq["fliC_i"], seq["invA"]])
    plasmid = _by_scheme(sal.run(path, db_dir, workdir=tmp_path / "work"))[
        "salmonella_virulence_plasmid"]
    assert plasmid.call == "spv absent"
    assert "expected negative" in plasmid.note


@blast_only
def test_vi_capsule_and_typhoid_toxin_are_reported_together(mini_db, tmp_path):
    db_dir, seq = mini_db
    path = _assembly(tmp_path / "typhoidal.fna", [
        seq["wzx_B"], seq["fliC_i"], seq["invA"],
        seq["tviB"], seq["vexC"], seq["cdtB"], seq["pltA"], seq["pltB"],
    ])
    results = _by_scheme(sal.run(path, db_dir, workdir=tmp_path / "work"))
    typhoidal = results["salmonella_typhoidal"]
    assert "typhoid toxin" in typhoidal.call
    assert "cdtB, pltA and pltB are all present" in typhoidal.note
    # Only two of the ten viaB genes are in the mini database, so Vi must come
    # back as not determined rather than present.
    assert "not determined rather than absent" in typhoidal.note


@blast_only
def test_a_genome_with_no_antigen_locus_says_so(mini_db, tmp_path):
    db_dir, _seq = mini_db
    write_fasta(tmp_path / "empty.fna", [("contig_1", _filler(30000))], wrap=70)
    results = _by_scheme(sal.run(tmp_path / "empty.fna", db_dir,
                                 workdir=tmp_path / "work"))
    serovar = results["salmonella_serovar"]
    assert serovar.call == "not determined"
    assert serovar.loci_found == 0
    assert "may not be Salmonella" in serovar.note
    assert results["salmonella_subspecies"].call == "not determined"


@blast_only
def test_a_gzipped_assembly_is_accepted(mini_db, tmp_path):
    """BLAST+ cannot read a gzipped FASTA, so the module expands it first."""
    db_dir, seq = mini_db
    plain = _assembly(tmp_path / "gz.fna",
                      [seq["wzx_B"], seq["wzy_B"], seq["fliC_i"], seq["fljA"],
                       seq["fljB_1_2"], seq["hin"]])
    packed = tmp_path / "gz.fna.gz"
    with open(plain, "rb") as src, gzip.open(packed, "wb") as out:
        shutil.copyfileobj(src, out)
    plain.unlink()
    results = _by_scheme(sal.run(packed, db_dir, workdir=tmp_path / "work"))
    assert results["salmonella_serovar"].call == "Typhimurium"


@blast_only
def test_no_call_is_ever_a_bare_negative(mini_db, tmp_path):
    """Every result must carry a note; a blank explanation is the bug the
    owner rejected earlier output for."""
    db_dir, seq = mini_db
    path = _assembly(tmp_path / "plain.fna",
                     [seq["wzx_B"], seq["fliC_i"], seq["invA"]])
    for result in sal.run(path, db_dir, workdir=tmp_path / "work"):
        assert result.note, f"{result.scheme} returned a call with no note"
        assert result.alleles, f"{result.scheme} returned no evidence"


# --------------------------------------------------------- validation (network)
#: Public assemblies whose serovar is published, with the citation for it.
#: "expected_call" is what this module should report, which for a shared
#: formula is the whole candidate list, not one name.
VALIDATION = {
    "GCF_000006945.2": {
        "strain": "LT2",
        "published_serovar": "Typhimurium (1,4,[5],12:i:1,2)",
        "source": "McClelland et al. 2001, Nature 413:852 (PMID 11677609)",
        "serovar": "Typhimurium",
        "formula": "I 1,4,[5],12:i:1,2",
        "st": "ST19",
    },
    "GCF_000009505.1": {
        "strain": "P125109",
        "published_serovar": "Enteritidis PT4 (1,9,12:g,m:-)",
        "source": "Thomson et al. 2008, Genome Res 18:1624 (PMID 18583645)",
        "serovar": "Enteritidis",
        "formula": "I 1,9,12:g,m:-",
        "st": "ST11",
    },
    "GCF_000195995.1": {
        "strain": "CT18",
        "published_serovar": "Typhi (9,12,Vi:d:-)",
        "source": "Parkhill et al. 2001, Nature 413:848 (PMID 11677608)",
        "serovar": "Typhi",
        "formula": "I 9,12[Vi]:d:-",
        "st": "ST2",
        "typhoidal": "Vi capsule",
    },
    "GCF_000020705.1": {
        "strain": "SL476",
        "published_serovar": "Heidelberg (1,4,[5],12:r:1,2)",
        "source": "Fricke et al. 2011, J Bacteriol 193:3556 (PMID 21602358)",
        "serovar": "Heidelberg",
        "formula": "I 1,4,[5],12:r:1,2",
        "st": "ST15",
    },
    "GCF_000016045.1": {
        "strain": "SL254",
        "published_serovar": "Newport (6,8,20:e,h:1,2)",
        "source": "Fricke et al. 2011, J Bacteriol 193:3556 (PMID 21602358)",
        "serovar": "Newport|Bardo",
        "st": "ST45",
    },
    "GCF_010919335.1": {
        "strain": "119944",
        "published_serovar": "Infantis (6,7,14:r:1,5)",
        "source": "Aviv et al. 2014, Environ Microbiol 16:977 (PMID 24320043)",
        "serovar": "Infantis",
        "formula": "I 6,7,14:r:1,5",
        "st": "ST32",
    },
    "GCF_002210025.1": {
        "strain": "TW-Stm6",
        "published_serovar": "Typhimurium var. monophasic 4,5,12:i:-",
        "source": "NCBI assembly record GCF_002210025.1; the monophasic "
                  "phenotype and ST34 clone are Petrovska et al. 2016, "
                  "Emerg Infect Dis 22:617 (PMID 26982594)",
        "serovar": "I 1,4,[5],12:i:-",
        "formula": "I 1,4,[5],12:i:-",
        "st": "ST34",
        "h2_locus": "fljBA-hin absent",
    },
    "GCF_000009525.1": {
        "strain": "287/91",
        "published_serovar": "Gallinarum (1,9,12:-:-), non-motile",
        "source": "Thomson et al. 2008, Genome Res 18:1624 (PMID 18583645)",
        # The antigen genes are Enteritidis'. This is the documented limit, not
        # a failure: motility is a phenotype.
        "serovar": "Enteritidis",
        "formula": "I 1,9,12:g,m:-",
        "st": "ST331",
        "note_contains": "Gallinarum",
    },
    "GCF_000018705.1": {
        "strain": "SPB7",
        "published_serovar": "Paratyphi B (1,4,[5],12:b:1,2)",
        "source": "NCBI assembly record GCF_000018705.1; the Paratyphi B / "
                  "var. Java pair share one formula per Grimont & Weill 2007",
        "serovar": "Paratyphi B|Paratyphi B var. Java",
        "formula": "I 1,4,[5],12:b:1,2",
        "st": "ST307",
    },
    "GCF_000011885.1": {
        "strain": "ATCC 9150",
        "published_serovar": "Paratyphi A (1,2,12:a:[1,5])",
        "source": "McClelland et al. 2004, Nat Genet 36:1268 (PMID 15531882)",
        "serovar": "Paratyphi A",
        "formula": "I 1,2,12:a:1,5",
        "st": "ST85",
    },
    "GCF_000018625.1": {
        "strain": "RSK2980",
        "published_serovar": "subsp. arizonae (IIIa) 62:z4,z23:-",
        "source": "NCBI assembly record GCF_000018625.1 (RefSeq organism name)",
        # A known discordance, kept in the suite on purpose: wzx/wzy place this
        # genome in O group S (O:41), not O:62. See docs/typing/salmonella.md.
        "serovar": "IIIa 41:z4,z23:-|IV 41:z4,z23:-|Waycross",
        "expect_o_group_ambiguity": "IIIa 62:z4,z23:-",
        "st": "ST2131",
    },
}


def _download_assembly(accession: str, into: Path) -> Path:
    dest = into / f"{accession}.fna"
    if dest.exists():
        return dest
    url = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
           f"{accession}/dataset_report")
    request = urllib.request.Request(url, headers={"User-Agent": sal._USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as handle:
        report = json.load(handle)["reports"][0]
    name = report["assembly_info"]["assembly_name"]
    digits = accession.split("_")[1].split(".")[0]
    parts = "/".join(digits[i:i + 3] for i in range(0, 9, 3))
    fna = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/"
           f"{accession.split('_')[0]}/{parts}/{accession}_{name}/"
           f"{accession}_{name}_genomic.fna.gz")
    archive = into / f"{accession}.fna.gz"
    sal._download(fna, archive)
    with gzip.open(archive, "rb") as src, open(dest, "wb") as out:
        shutil.copyfileobj(src, out)
    archive.unlink()
    return dest


@pytest.fixture(scope="module")
def real_db(tmp_path_factory):
    """The real database, built once. Honours HYDRA_SALMONELLA_DB if set."""
    existing = os.environ.get("HYDRA_SALMONELLA_DB")
    if existing:
        return Path(existing)
    root = tmp_path_factory.mktemp("salmonella-db")
    sal.fetch(root)
    return root


@validation_only
@blast_only
@pytest.mark.parametrize("accession", sorted(VALIDATION))
def test_published_serovars_are_reproduced(accession, real_db, tmp_path):
    expected = VALIDATION[accession]
    genomes = Path(os.environ.get("HYDRA_SALMONELLA_GENOMES", tmp_path))
    genomes.mkdir(parents=True, exist_ok=True)
    assembly = _download_assembly(accession, genomes)
    results = _by_scheme(sal.run(assembly, real_db, threads=4,
                                 workdir=tmp_path / "work"))
    serovar = results["salmonella_serovar"]
    assert serovar.call == expected["serovar"], (
        f"{accession} ({expected['strain']}, published "
        f"{expected['published_serovar']}, source {expected['source']}): "
        f"expected {expected['serovar']!r}, got {serovar.call!r}")
    if "formula" in expected:
        assert serovar.alleles["antigenic_formula"] == expected["formula"]
    if "st" in expected:
        assert results["salmonella_mlst"].call == expected["st"]
    if "h2_locus" in expected:
        assert results["salmonella_h2_locus"].call == expected["h2_locus"]
    if "typhoidal" in expected:
        assert expected["typhoidal"] in results["salmonella_typhoidal"].call
    if "note_contains" in expected:
        assert expected["note_contains"] in serovar.note
    if "expect_o_group_ambiguity" in expected:
        assert expected["expect_o_group_ambiguity"] in \
            serovar.alleles["o_group_ambiguity"]


@validation_only
@blast_only
def test_fetch_is_idempotent(real_db):
    """Running fetch twice must not re-download or break the install."""
    before = {p.name: p.stat().st_size
              for p in (real_db / sal.KEY).rglob("*") if p.is_file()}
    sal.fetch(real_db)
    after = {p.name: p.stat().st_size
             for p in (real_db / sal.KEY).rglob("*") if p.is_file()}
    assert set(before) == set(after)


@validation_only
@blast_only
def test_the_installed_database_records_its_own_provenance(real_db):
    provenance = json.loads((real_db / sal.KEY / "provenance.json").read_text())
    assert provenance["key"] == sal.KEY
    names = {entry["name"] for entry in provenance["databases"]}
    assert names == {db.name for db in sal.DATABASES}
    for entry in provenance["databases"]:
        assert entry["licence"] and entry["redistribution"]
