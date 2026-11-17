"""The integration layer: discovery, dispatch, database wiring and reporting.

These tests do not type anything -- each module has its own tests for that.
They cover the parts that only exist once the modules are wired together: which
module a sample is sent to, what happens when its database is missing, and
whether the evidence and the caveats reach the output.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from hydra_amr.db.registry import DATABASES, ORGANISM_DATABASES, resolve_names
from hydra_amr.engines.blast import TRANSLATED_PROGRAMS, safe_threads
from hydra_amr.records import MlstCall, SampleResult, SpeciesCall, TypingResult
from hydra_amr.report.tables import (TYPING_COLUMNS, summary_table, typing_overview,
                                     typing_table)
from hydra_amr.typing.organisms import registry

EXPECTED_MODULES = {"ecoli", "hinfluenzae", "nmeningitidis", "sagalactiae",
                    "salmonella", "saureus", "spneumoniae", "spyogenes"}


# ------------------------------------------------------------------ discovery
def test_every_module_is_discovered():
    assert set(registry.modules()) == EXPECTED_MODULES


def test_shared_helpers_are_not_mistaken_for_modules():
    # strep_common and _hinf_nmen_common live in the same package and must not
    # be dispatched: they have no KEY and no run().
    assert "strep_common" not in registry.modules()
    assert "_hinf_nmen_common" not in registry.modules()


def test_every_module_exposes_the_agreed_interface():
    for key, entry in registry.modules().items():
        assert entry.key == key
        assert entry.title
        assert entry.species or entry.mlst_schemes
        assert entry.databases, f"{key} declares no databases"
        assert entry.references, f"{key} declares no references"
        assert callable(entry.module.fetch)
        assert callable(entry.module.run)
        for db in entry.databases:
            assert db.redistribution in ("permitted", "restricted", "unstated")
            assert db.licence
            assert db.source


def test_module_names_list_matches_what_is_on_disk():
    # The explicit list is what a frozen Windows build dispatches from; a module
    # added without being listed would silently never run there.
    assert set(registry.MODULE_NAMES) == EXPECTED_MODULES


# ------------------------------------------------------------------ selection
@pytest.mark.parametrize("species,expected", [
    ("Staphylococcus aureus", {"saureus"}),
    ("Streptococcus pyogenes", {"spyogenes"}),
    ("Streptococcus agalactiae", {"sagalactiae"}),
    ("Streptococcus pneumoniae", {"spneumoniae"}),
    ("Haemophilus influenzae", {"hinfluenzae"}),
    ("Neisseria meningitidis", {"nmeningitidis"}),
    ("Salmonella enterica", {"salmonella"}),
    ("Escherichia coli", {"ecoli"}),
    ("Shigella sonnei", {"ecoli"}),
])
def test_species_selects_its_module(species, expected):
    call = SpeciesCall(name=species, genus=species.split()[0],
                       species=species.split()[-1])
    assert {m.key for m in registry.select(call, None)} == expected


def test_a_strain_level_species_name_still_selects():
    call = SpeciesCall(name="Escherichia coli O157:H7 str. Sakai")
    assert [m.key for m in registry.select(call, None)] == ["ecoli"]


def test_a_sister_species_is_not_confused_with_its_neighbour():
    # "Streptococcus pneumoniae" must never match a "Streptococcus pyogenes"
    # claim by prefix, and vice versa.
    for name, key in (("Streptococcus pneumoniae", "spneumoniae"),
                      ("Streptococcus pyogenes", "spyogenes")):
        assert [m.key for m in registry.select(SpeciesCall(name=name), None)] == [key]


def test_s_dysgalactiae_does_not_run_the_pyogenes_module():
    """A group C/G streptococcus is not a GAS and must not be typed as one.

    The CDC emm database holds the stG/stC/stL types of *S. dysgalactiae*, so
    running the *S. pyogenes* module on one would produce an emm-shaped call for
    an organism none of its schemes were defined in.
    """
    for name in ("Streptococcus dysgalactiae",
                 "Streptococcus dysgalactiae subsp. equisimilis"):
        call = SpeciesCall(name=name, genus="Streptococcus",
                           species="dysgalactiae")
        assert registry.select(call, None) == []


def test_mlst_scheme_selects_when_the_species_call_is_empty():
    assert [m.key for m in registry.select(SpeciesCall(),
                                           MlstCall(scheme="saureus",
                                                    sequence_type="8"))] == ["saureus"]
    assert [m.key for m in registry.select(None,
                                           MlstCall(scheme="senterica_achtman_2",
                                                    sequence_type="19"))] == ["salmonella"]


def test_the_species_call_wins_over_the_mlst_scheme():
    # Neighbouring organisms share housekeeping alleles; the species call is the
    # stronger evidence when both are present.
    chosen = registry.select(SpeciesCall(name="Staphylococcus aureus"),
                             MlstCall(scheme="ecoli", sequence_type="131"))
    assert [m.key for m in chosen] == ["saureus"]


def test_unknown_species_and_no_mlst_selects_nothing():
    assert registry.select(SpeciesCall(name="unknown"), MlstCall()) == []
    assert registry.select(None, None) == []


def test_only_restricts_the_selection():
    call = SpeciesCall(name="Escherichia coli")
    assert registry.select(call, None, only=["saureus"]) == []
    assert [m.key for m in registry.select(call, None, only=["ecoli"])] == ["ecoli"]


# ------------------------------------------------- missing databases are said
def test_a_module_with_no_database_is_reported_not_silently_skipped(tmp_path):
    warn: list[str] = []
    out = registry.type_assembly(tmp_path / "nothing.fna", tmp_path / "db",
                                 species=SpeciesCall(name="Staphylococcus aureus"),
                                 warn=warn)
    assert out == []
    assert len(warn) == 1
    assert "not installed" in warn[0]
    assert "hydra db download typing-saureus" in warn[0]


def test_installed_needs_at_least_one_declared_database(tmp_path):
    entry = registry.get("saureus")
    root = tmp_path / "saureus"
    assert not entry.installed(tmp_path)
    root.mkdir(parents=True)
    assert not entry.installed(tmp_path), "an empty directory is not an install"
    (root / "markers").mkdir()
    assert entry.installed(tmp_path)
    assert entry.missing_databases(tmp_path) == ["spa", "mlst"]


def test_a_module_that_raises_does_not_end_the_run(tmp_path, monkeypatch):
    entry = registry.get("saureus")
    (tmp_path / "saureus" / "markers").mkdir(parents=True)

    def boom(*_args, **_kwargs):
        from hydra_amr.utils import HydraError

        raise HydraError("the marker panel is corrupt")

    monkeypatch.setattr(entry.module, "run", boom)
    warn: list[str] = []
    out = registry.type_assembly(tmp_path / "x.fna", tmp_path,
                                 species=SpeciesCall(name="Staphylococcus aureus"),
                                 warn=warn)
    assert out == []
    assert "the marker panel is corrupt" in warn[0]


# ------------------------------------------------------- MLST pass-through
def test_modules_that_can_reuse_hydras_mlst_call_accept_it():
    # These four re-derived the seven housekeeping loci before integration.
    import inspect

    for key in ("saureus", "hinfluenzae", "nmeningitidis", "salmonella"):
        params = inspect.signature(registry.get(key).module.run).parameters
        assert "mlst" in params, f"{key}.run() cannot be given Hydra's MLST call"


def test_saureus_names_the_clonal_complex_from_hydras_call():
    from hydra_amr.typing.organisms import saureus

    st, alleles, problem = saureus._st_from_call(
        MlstCall(scheme="saureus", sequence_type="8",
                 alleles={"arcC": "3", "aroE": "3"}))
    assert st == "8"
    assert alleles["arcC"] == "3"
    assert problem == ""


def test_saureus_refuses_an_mlst_call_from_another_scheme():
    from hydra_amr.typing.organisms import saureus

    st, _alleles, problem = saureus._st_from_call(
        MlstCall(scheme="ecoli", sequence_type="131"))
    assert st == ""
    assert "not the" in problem and "S. aureus" in problem


def test_saureus_no_longer_downloads_a_second_copy_of_the_mlst_alleles():
    import inspect

    from hydra_amr.typing.organisms import saureus

    source = inspect.getsource(saureus.fetch)
    assert "alleles_fasta" not in source
    assert "profiles.tsv" in source or "PUBMLST_PROFILES" in source


def test_salmonella_builds_the_mlst_result_from_hydras_call(tmp_path):
    from hydra_amr.typing.organisms import salmonella

    profiles = tmp_path / "profiles.tsv"
    profiles.write_text(
        "ST\taroC\tdnaN\themD\thisD\tpurE\tsucA\tthrA\teBG\n"
        "19\t10\t19\t12\t9\t5\t9\t2\t1\n")
    serovar = TypingResult(scheme="serovar", call="Typhimurium",
                           alleles={"antigenic_formula": "1,4,[5],12:i:1,2"})
    out = salmonella._mlst_from_call(
        MlstCall(scheme="senterica_achtman_2", sequence_type="19",
                 alleles={"aroC": "10"}, loci_found=7, loci_total=7),
        profiles, serovar)
    assert out is not None
    assert out.call == "ST19"
    assert out.lineage == "eBG1"
    assert "not typed twice" in out.note
    assert salmonella._mlst_from_call(MlstCall(), profiles, serovar) is None
    assert salmonella._mlst_from_call(
        MlstCall(scheme="saureus", sequence_type="8"), profiles, serovar) is None


def test_hinf_nmen_profile_from_call_checks_the_scheme(tmp_path):
    from hydra_amr.typing.organisms import _hinf_nmen_common as common

    profiles = tmp_path / "profiles.txt"
    profiles.write_text("ST\tabcZ\tadk\taroE\tfumC\tgdh\tpdhC\tpgm\tclonal_complex\n"
                        "33\t1\t3\t2\t1\t1\t1\t3\tST-32 complex\n")
    good = common.profile_from_call(
        MlstCall(scheme="neisseria", sequence_type="33", loci_found=7, loci_total=7),
        "neisseria", profiles, accepted=("neisseria",))
    assert good is not None
    assert good.sequence_type == "33"
    assert good.clonal_complex == "ST-32 complex"
    assert common.profile_from_call(MlstCall(scheme="saureus", sequence_type="8"),
                                    "neisseria", profiles,
                                    accepted=("neisseria",)) is None
    assert common.profile_from_call(None, "neisseria", profiles) is None


# ------------------------------------------------------------ the BLAST+ fix
def test_translated_searches_are_pinned_to_one_thread():
    for program in TRANSLATED_PROGRAMS:
        assert safe_threads(program, 16) == 1
        assert safe_threads(program, 1) == 1


def test_nucleotide_searches_keep_their_threads():
    for program in ("blastn", "megablast", "dc-megablast"):
        assert safe_threads(program, 16) == 16
    assert safe_threads("blastn", 0) == 1


def test_no_module_pins_its_own_thread_count_any_more():
    """The workaround lived in three modules; the fix lives in blast() now."""
    root = Path(__file__).resolve().parents[1] / "src" / "hydra_amr" / "typing" / "organisms"
    for path in sorted(root.glob("*.py")):
        text = path.read_text()
        assert "BLASTX_THREADS" not in text, path.name
        for program in ("blastx", "tblastn"):
            assert f'blast("{program}"' not in text or "threads=1," not in text, path.name


# ------------------------------------------------------- evidence in reports
def _sample() -> SampleResult:
    result = SampleResult(sample="ISOLATE-1")
    result.species = SpeciesCall(name="Staphylococcus aureus", confidence="strong")
    result.typing = registry.relabel([
        TypingResult(scheme="spa", call="t008", lineage="-",
                     alleles={"repeats": "11-19-12-21-17-34-24-34-22-25",
                              "evidence": "contig_1:120-410 100.0% identity"},
                     loci_found=1, loci_total=1,
                     note="the Ridom spa database is a non-commercial licence"),
        TypingResult(scheme="agr", call="not determined", lineage="-",
                     note="no agrBDCA locus reached the identity floor"),
    ])
    return result


def test_the_long_typing_table_carries_the_evidence_and_the_caveats():
    frame = typing_table([_sample()])
    assert list(frame.columns) == list(TYPING_COLUMNS)
    spa = frame[frame["scheme"] == "spa"].iloc[0]
    assert spa["call"] == "t008"
    assert "11-19-12" in spa["evidence"]
    assert "100.0% identity" in spa["evidence"]
    assert "non-commercial" in spa["caveats"]
    assert spa["scheme_name"] == "spa type (Ridom)"
    assert spa["loci"] == "1/1"


def test_a_not_determined_call_says_why():
    frame = typing_table([_sample()])
    agr = frame[frame["scheme"] == "agr"].iloc[0]
    assert agr["call"] == "not determined"
    assert "identity floor" in agr["caveats"]


def test_the_overview_puts_the_caveat_beside_the_call():
    frame = typing_overview([_sample()])
    row = frame.iloc[0]
    assert row["spa"] == "t008"
    assert "non-commercial" in row["spa_note"]
    # A call with nothing to qualify gets no note column of its own.
    assert "agr" in frame.columns


def test_the_summary_still_names_every_call():
    frame = summary_table([_sample()])
    assert frame.iloc[0]["spa"] == "t008"
    assert frame.iloc[0]["agr"] == "not determined"


def test_an_empty_typing_table_still_has_its_header():
    frame = typing_table([SampleResult(sample="x")])
    assert frame.empty
    assert list(frame.columns) == list(TYPING_COLUMNS)


def test_the_html_report_shows_the_evidence_and_the_caveats():
    from hydra_amr.report.html import render

    page = render([_sample()])
    assert "Organism typing: the calls and their evidence" in page
    assert "11-19-12-21-17-34-24-34-22-25" in page
    assert "non-commercial" in page
    assert "identity floor" in page


def test_the_writer_emits_both_typing_artefacts(tmp_path):
    from hydra_amr.report.writer import write_outputs

    written = write_outputs([_sample()], outdir=tmp_path, prefix="t",
                            formats=("tsv",))
    names = {p.name for p in written}
    assert "t.typing.tsv" in names
    assert "t.typing.overview.tsv" in names
    body = (tmp_path / "t.typing.tsv").read_text()
    assert "11-19-12-21-17-34-24-34-22-25" in body
    assert "non-commercial" in body


# ---------------------------------------------------- the group C/G relabel
@pytest.mark.parametrize("call", ["stG62647.0", "stC74a.23", "stL1929.3",
                                  "emmG1.0"])
def test_a_group_cg_emm_record_is_not_presented_as_a_gas_emm_type(call):
    [out] = registry.relabel([TypingResult(scheme="emm", call=call)])
    assert "emm type" not in out.heading
    assert "dysgalactiae" in out.heading
    assert "NOT a Streptococcus" in out.note
    assert call in out.note
    # The call itself is untouched: the module's nomenclature is correct.
    assert out.call == call


def test_a_real_emm_type_keeps_the_emm_heading():
    [out] = registry.relabel([TypingResult(scheme="emm", call="emm1.0")])
    assert out.heading == "emm type (M protein)"
    assert out.note == ""


def test_relabel_never_changes_a_call_or_its_alleles():
    before = TypingResult(scheme="cps", call="III", alleles={"locus": "cps"},
                          note="keep me")
    [after] = registry.relabel([before])
    assert (after.call, after.lineage, after.alleles) == (
        before.call, before.lineage, before.alleles)
    assert after.note == "keep me"
    assert after.heading == "Capsular serotype"


def test_an_unknown_scheme_still_gets_a_readable_heading():
    [out] = registry.relabel([TypingResult(scheme="some_new_scheme", call="x")])
    assert out.heading == "some new scheme"


# --------------------------------------------------------- database wiring
def test_each_module_is_installable_through_hydra_db():
    assert set(ORGANISM_DATABASES) == {f"typing-{k}" for k in EXPECTED_MODULES}
    for name in ORGANISM_DATABASES:
        spec = DATABASES[name]
        assert spec.kind == "organism"
        assert spec.url, f"{name} records no source"
        assert spec.licence, f"{name} records no licence"
        assert "Redistribution:" in spec.notes


def test_the_organism_databases_are_resolvable_by_name_and_as_a_group():
    assert resolve_names(["typing-saureus"]) == ["typing-saureus"]
    assert resolve_names(["saureus"]) == ["typing-saureus"]
    assert set(resolve_names(["organisms"])) == set(ORGANISM_DATABASES)


def test_the_organism_databases_stay_out_of_the_gene_screening_groups():
    for group in ("all", "standard", "nucl", "amr", "virulence", "core"):
        assert not set(resolve_names([group])) & set(ORGANISM_DATABASES)


def test_an_organism_database_can_be_downloaded_but_not_imported(tmp_path):
    from hydra_amr.db.manager import DatabaseStore
    from hydra_amr.db.registry import spec_for
    from hydra_amr.utils import HydraError

    store = DatabaseStore(tmp_path)
    with pytest.raises(HydraError, match="hydra db download"):
        store.import_one(spec_for("typing-saureus"), tmp_path)
    assert store.organism_dir == tmp_path / "organisms"


def test_db_download_knows_how_to_fetch_an_organism_database():
    from hydra_amr.db.fetch import can_fetch

    assert can_fetch("typing-saureus")
    assert not can_fetch("argannot")   # a landing page, still manual


def test_install_organism_records_the_licence_it_was_given(tmp_path, monkeypatch):
    from hydra_amr.db.manager import DatabaseStore

    store = DatabaseStore(tmp_path)
    entry = registry.get("spyogenes")

    def fake_fetch(db_dir: Path) -> None:
        for db in entry.databases:
            (Path(db_dir) / entry.key / db.name).mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(entry.module, "fetch", fake_fetch)
    store.install_organism("typing-spyogenes")
    recorded = store.installed()["typing-spyogenes"]
    assert recorded["kind"] == "organism"
    assert recorded["path"] == "organisms/spyogenes"
    assert "emm" in recorded["databases"]
    assert recorded["missing"] == []
    assert "public domain" in recorded["licence"]
    assert store.is_installed("typing-spyogenes")


# ----------------------------------------------------------- pipeline gating
def _pipeline(tmp_path, **option_overrides):
    from hydra_amr.config import Config
    from hydra_amr.db.manager import DatabaseStore
    from hydra_amr.pipeline import Pipeline, RunOptions

    return Pipeline(Config(), DatabaseStore(tmp_path), RunOptions(**option_overrides))


def test_no_typing_turns_organism_typing_off_too(tmp_path, monkeypatch):
    called: list[str] = []
    monkeypatch.setattr(registry, "select",
                        lambda *a, **k: called.append("selected") or [])
    pipe = _pipeline(tmp_path, typing=False)
    pipe._type_organisms({"s": tmp_path / "s.fna"}, {"s": SampleResult(sample="s")},
                         tmp_path)
    assert called == []


def test_organism_typing_can_be_turned_off_on_its_own(tmp_path, monkeypatch):
    called: list[str] = []
    monkeypatch.setattr(registry, "select",
                        lambda *a, **k: called.append("selected") or [])
    pipe = _pipeline(tmp_path, typing=True, organism_typing=False)
    pipe._type_organisms({"s": tmp_path / "s.fna"}, {"s": SampleResult(sample="s")},
                         tmp_path)
    assert called == []


def test_the_pipeline_records_why_a_module_was_skipped(tmp_path):
    result = SampleResult(sample="s")
    result.species = SpeciesCall(name="Staphylococcus aureus")
    pipe = _pipeline(tmp_path)
    pipe._type_organisms({"s": tmp_path / "s.fna"}, {"s": result}, tmp_path)
    assert result.typing == []
    assert any("not installed" in w for w in result.warnings)


def test_the_pipeline_types_the_species_it_was_given(tmp_path, monkeypatch):
    seen: dict = {}

    def fake_type_assembly(assembly, db_dir, **kwargs):
        seen.update(kwargs)
        seen["db_dir"] = db_dir
        return [TypingResult(scheme="spa", call="t008")]

    monkeypatch.setattr(registry, "type_assembly", fake_type_assembly)
    result = SampleResult(sample="s")
    result.species = SpeciesCall(name="Staphylococcus aureus")
    result.mlst = MlstCall(scheme="saureus", sequence_type="8")
    pipe = _pipeline(tmp_path)
    pipe._type_organisms({"s": tmp_path / "s.fna"}, {"s": result}, tmp_path)
    assert [t.call for t in result.typing] == ["t008"]
    assert seen["db_dir"] == tmp_path / "organisms"
    assert seen["mlst"].sequence_type == "8"
    assert seen["species"].name == "Staphylococcus aureus"
