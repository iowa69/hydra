"""*Streptococcus pneumoniae*: capsular serotype, PBP type, and the rest of what
a reference laboratory reports from an assembly.

Six schemes are called, in the order a clinical report needs them.

**Capsular serotype** is decided in the two stages the reference laboratories
use. The first asks which of the 106 sequenced *cps* loci the genome carries,
measured as coverage of the reference locus - that settles the serogroup, tells
a non-encapsulated genome from an encapsulated one, and is the evidence trail.
The second applies the CDC Streptococcus Laboratory's serotype-discriminating
probe set and its published decision rules, because within a serogroup the
serotypes differ by a single residue of *wciP*, or by whether one
acetyltransferase is intact, and whole-locus similarity cannot see either.

Where the discriminator is missing or ambiguous the pair is reported as a pair -
``9A/9V``, ``18B/18C``, ``35B/35D``, ``7A/7F`` - never resolved by guessing.
``15B`` and ``15C`` differ only by the length of a TA repeat in *wciZ* that
shifts reversibly in culture, so that call carries its own warning.

**PBP typing** reproduces the CDC pbp1a/pbp2b/pbp2x transpeptidase-region
allele triplet: the three regions are extracted from the assembly, translated,
and matched at 100% identity over their full length against the CDC allele sets.
A region with no exact match is reported ``NEW``, with the substitutions it
carries in the transpeptidase motifs, rather than rounded to the nearest number.

**Beta-lactam MIC** is then read from the CDC reference table of PBP type to
MIC, and interpreted against the CLSI breakpoints the CDC pipeline applies.
This is a *prediction* from a genotype, and every field that carries it says so.

The three remaining schemes are the pilus islets, the macrolide-lincosamide
mechanism (which Hydra's AMR screen detects as genes but does not interpret as a
phenotype), and the fluoroquinolone/co-trimoxazole target substitutions.
"""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

from ...engines.blast import Hsp, blast, check_db_exists, merge_hsps
from ...records import TypingResult
from ...seqio import read_fasta, write_fasta
from ...utils import LOG, HydraError, revcomp, translate
from .base import OrganismDb, Reference

KEY = "spneumoniae"
TITLE = "Streptococcus pneumoniae"
#: Species names, as Hydra's SpeciesCall reports them, that trigger the module.
SPECIES: tuple[str, ...] = ("Streptococcus pneumoniae",)
#: PubMLST scheme names that also trigger it.
MLST_SCHEMES: tuple[str, ...] = ("spneumoniae",)


# --------------------------------------------------------------- data sources
#: Pinned commits, so a database build is reproducible and a silent upstream
#: change cannot alter a call. Both repositories are plain Git over HTTPS.
_CPS_COMMIT = "master"
_CDC_COMMIT = "master"
_GITHUB = "https://raw.githubusercontent.com"
_CPS_RAW = f"{_GITHUB}/GlobalPneumoSeq/cps_reference_database/{_CPS_COMMIT}"
_CDC_RAW = f"{_GITHUB}/BenJamesMetcalf/Spn_Scripts_Reference/{_CDC_COMMIT}"

#: (url, path under db_dir/<KEY>/, BLAST dbtype or "" for a plain data file,
#: the DATABASES entry recording its provenance)
_FILES: tuple[tuple[str, str, str, str], ...] = (
    (f"{_CPS_RAW}/all.fasta", "cps/cps_loci.fna", "nucl", "cps"),
    (f"{_CDC_RAW}/SPN_Reference_DB/SPN_Sero_Gene-DB_Final.fasta",
     "cps/sero_probes.fna", "nucl", "cdc-serotype-probes"),
    (f"{_CDC_RAW}/SPN_Reference_DB/SPN_Res_Gene-DB_Final.fasta",
     "res/res_probes.fna", "nucl", "cdc-resistance-probes"),
    (f"{_CDC_RAW}/SPN_Reference_DB/MOD_bLactam_resistance.fasta",
     "pbp/pbp_regions.fna", "nucl", "cdc-pbp"),
    (f"{_CDC_RAW}/SPN_Reference_DB/SPN_bLactam_1A-DB.faa",
     "pbp/pbp_1A.faa", "prot", "cdc-pbp"),
    (f"{_CDC_RAW}/SPN_Reference_DB/SPN_bLactam_2B-DB.faa",
     "pbp/pbp_2B.faa", "prot", "cdc-pbp"),
    (f"{_CDC_RAW}/SPN_Reference_DB/SPN_bLactam_2X-DB.faa",
     "pbp/pbp_2X.faa", "prot", "cdc-pbp"),
    (f"{_CDC_RAW}/bLactam_MIC_Rscripts/newDB/Ref_PBPtype_MIC.csv",
     "pbp/pbp_type_mic.csv", "", "cdc-pbp"),
)

DATABASES: tuple[OrganismDb, ...] = (
    OrganismDb(
        name="cps",
        source="https://github.com/GlobalPneumoSeq/cps_reference_database",
        version="master, retrieved 2026-11-17 (repository last pushed 2025-08-13)",
        licence="none stated: the repository carries no LICENSE file and no "
                "licence statement in its README",
        redistribution="unstated",
        reference="Bentley2006",
        notes="all.fasta, 106 capsular (cps) locus reference sequences, 2.1 MB. "
              "The capsular loci themselves were sequenced by Bentley et al. "
              "2006 and deposited in GenBank; this repository is the curated, "
              "per-serotype set the Global Pneumococcal Sequencing project uses. "
              "Because no licence is stated, the integrator should download it "
              "rather than vendor it into the Windows zip unless the GPS "
              "project is asked and agrees.",
    ),
    OrganismDb(
        name="cdc-serotype-probes",
        source="https://github.com/BenJamesMetcalf/Spn_Scripts_Reference",
        version="master, retrieved 2026-11-17 (repository last pushed 2026-08-31)",
        licence="none stated: the repository carries no LICENSE file. It is the "
                "reference data of the CDC Streptococcus Laboratory; works of "
                "the United States federal government are generally not subject "
                "to domestic copyright (17 U.S.C. 105), but the repository "
                "itself makes no such statement.",
        redistribution="unstated",
        reference="Metcalf2016Sero",
        notes="SPN_Sero_Gene-DB_Final.fasta, 72 serotype-discriminating probes "
              "in 61 clusters, 14 KB. Includes the RRGA and PITB probes used "
              "for the pilus islets.",
    ),
    OrganismDb(
        name="cdc-pbp",
        source="https://github.com/BenJamesMetcalf/Spn_Scripts_Reference",
        version="master, retrieved 2026-11-17; PBP allele sets as of the "
                "2026-08-31 database update",
        licence="none stated (see cdc-serotype-probes)",
        redistribution="unstated",
        reference="Li2016",
        notes="SPN_bLactam_1A/2B/2X-DB.faa (533 / 744 / 1186 transpeptidase "
              "alleles, 0.81 MB total), MOD_bLactam_resistance.fasta (the three "
              "nucleotide reference fragments that locate the regions) and "
              "Ref_PBPtype_MIC.csv (307 PBP types with their reference MICs).",
    ),
    OrganismDb(
        name="cdc-resistance-probes",
        source="https://github.com/BenJamesMetcalf/Spn_Scripts_Reference",
        version="master, retrieved 2026-11-17",
        licence="none stated (see cdc-serotype-probes)",
        redistribution="unstated",
        reference="Metcalf2016Res",
        notes="SPN_Res_Gene-DB_Final.fasta, 19 probes, 3.2 KB: the macrolide, "
              "lincosamide, tetracycline, phenicol, fluoroquinolone and "
              "folate-pathway targets the CDC pipeline screens.",
    ),
)


# ---------------------------------------------------------------- thresholds
#: A probe counts as present at these levels. They are the BLAST equivalents of
#: the ``srst2 --min_coverage 99.8 --max_divergence 5`` the CDC pipeline runs
#: its serotype and resistance screens with.
PROBE_MIN_COVERAGE = 99.8
PROBE_MIN_IDENTITY = 95.0
#: Several CDC rules turn on whether a probe matched with *no* differences at
#: all ("diffs" in their output), which on an assembly means exactly this.
EXACT_IDENTITY = 100.0
EXACT_COVERAGE = 100.0

#: Fraction of a *cps* reference locus that must be covered before the genome is
#: accepted as carrying an intact capsule type. On the validation panel every
#: encapsulated genome covered its reference by 91.8% or more (91.8, 93.2, 96.1,
#: 97.5, 99.9, 99.9, 100, 100, 100, 100), while R6 - whose capsule locus carries
#: a 7.5 kb internal deletion - covered the serotype 2 reference by 63.6%. The
#: floor sits in that gap.
CPS_MIN_COVERAGE = 80.0
CPS_MIN_IDENTITY = 85.0
#: Below this, nothing recognisable as a capsule locus is present at all: the
#: conserved 5' *wzg-wzh-wzd-wze* block and the *dexB*/*aliA* flanks are shared
#: by every locus, so even a capsule-free genome covers a fifth to a third of
#: any reference.
CPS_ABSENT_COVERAGE = 35.0
#: Two whole-locus matches this close cannot be separated by coverage alone.
CPS_TIE_MARGIN = 2.0

#: The transpeptidase region must be recovered this completely before it is
#: typed. Both figures are the CDC pipeline's own: ``ExtractGene.pl`` is called
#: with ``-L 0.95`` and no ``-I``, and its documented default identity
#: threshold is 0.5. A permissive identity floor is not a licence to invent an
#: allele - the allele call itself still needs a 100% full-length protein match
#: - and it is necessary, because a mosaic *pbp1a* recombined from *S. mitis*
#: is only about 77% identical to the susceptible reference fragment. At an 80%
#: floor the two most penicillin-resistant genomes on the validation panel
#: failed to type at pbp1a at all.
PBP_MIN_REGION_COVERAGE = 95.0
PBP_MIN_REGION_IDENTITY = 50.0
#: Floors below which a translated region is not even called a novel allele,
#: again the CDC typer's own (identity >= 50% over >= 50% of the allele).
PBP_NEW_MIN_IDENTITY = 50.0
PBP_NEW_MIN_COVERAGE = 50.0

#: Thresholds for recovering the short target-gene regions (gyrA, parC, folA,
#: rplD). They are lower than the PBP ones because the probes are 36-210 bases
#: and a mosaic *folA* in a co-trimoxazole-resistant isolate aligns over only
#: 53 of the 60 reference bases at 83% identity; the unaligned flank is then
#: extended back from the reference length, as the CDC pipeline does.
TARGET_MIN_COVERAGE = 80.0
TARGET_MIN_IDENTITY = 75.0
#: A translated target region differing from the reference at more than this
#: fraction of its residues is not a divergent allele, it is the wrong region
#: or a frameshifted extraction, and is reported as not determined.
TARGET_MAX_DIVERGENT_FRACTION = 0.4


# ------------------------------------------------------------ vaccine formulae
#: Serotypes in each licensed formulation. Every list below was read from the
#: product's own FDA package insert or review memorandum on 2026-11-17, except
#: PCV7, PCV10-GSK and PCV13, which were taken from the Global Pneumococcal
#: Sequencing project's curated valency table and cross-checked against the
#: Prevnar 20 and Vaxneuvance inserts (the Prevnar 20 insert states the seven
#: PCV7 serotypes verbatim as its otitis media indication; the Vaxneuvance
#: insert's fifteen serotypes are PCV13 plus 22F and 33F).
VACCINES: dict[str, tuple[str, ...]] = {
    "PCV7": ("4", "6B", "9V", "14", "18C", "19F", "23F"),
    "PCV10-GSK": ("1", "4", "5", "6B", "7F", "9V", "14", "18C", "19F", "23F"),
    "PCV13": ("1", "3", "4", "5", "6A", "6B", "7F", "9V", "14", "18C", "19A",
              "19F", "23F"),
    "PCV15": ("1", "3", "4", "5", "6A", "6B", "7F", "9V", "14", "18C", "19A",
              "19F", "22F", "23F", "33F"),
    "PCV20": ("1", "3", "4", "5", "6A", "6B", "7F", "8", "9V", "10A", "11A",
              "12F", "14", "15B", "18C", "19A", "19F", "22F", "23F", "33F"),
    "PCV21": ("3", "6A", "7F", "8", "9N", "10A", "11A", "12F", "15A", "15B",
              "15C", "16F", "17F", "19A", "20A", "22F", "23A", "23B", "24F",
              "31", "33F", "35B"),
    "PPSV23": ("1", "2", "3", "4", "5", "6B", "7F", "8", "9N", "9V", "10A",
               "11A", "12F", "14", "15B", "17F", "18C", "19A", "19F", "20",
               "22F", "23F", "33F"),
}

#: Serotypes a formulation does not contain but is expected to cover, because
#: the capsule differs from an included one only in O-acetylation or in one
#: sugar and the antibodies cross-react. Reported separately from membership,
#: never folded into it.
VACCINE_CROSS_REACTIVE: dict[str, tuple[str, ...]] = {
    # The 6A conjugate raises antibody that opsonises 6C, whose capsule differs
    # from 6A only in the wciN-encoded sugar (Cooper 2011).
    "6C": ("PCV13", "PCV15", "PCV20"),
    # 15C is de-O-acetylated 15B and the two interconvert by a wciZ repeat shift
    # (van Selm 2003); PCV20 lists 15B only, PCV21 lists both.
    "15C": ("PCV20", "PPSV23"),
    # Danish type 20 was split into 20A and 20B after PPSV23 was formulated
    # (Calix 2012), so a "20" in that insert covers both.
    "20A": ("PPSV23",),
    "20B": ("PPSV23",),
}

#: Serotype to serogroup: read off the name, with the exceptions named.
_SEROGROUP_EXCEPTIONS = {"37": "37 (tts, outside the cps locus)"}

#: Serogroups whose members no current method separates from sequence, so a
#: whole-locus match inside one of them is reported as the group. Each entry
#: gives the members and why they cannot be split.
_UNSPLITTABLE_GROUPS: dict[str, tuple[tuple[str, ...], str]] = {
    "32": (("32A", "32F"),
           "no nucleotide difference has been found within the coding regions "
           "of the 32A and 32F loci (Kapatai et al. 2016, on five isolates); "
           "neither PneumoCaT nor PneumoKITy attempts the split"),
    "24": (("24A", "24B", "24F"),
           "the loci do differ, but no validated rule for separating them "
           "exists: PneumoCaT ships no serogroup 24 variant data and SeroBA "
           "reports serogroup 24 rather than a serotype"),
}


def serogroup_of(serotype: str) -> str:
    """The serogroup a serotype belongs to, as the field writes it."""
    if serotype in _SEROGROUP_EXCEPTIONS:
        return _SEROGROUP_EXCEPTIONS[serotype]
    match = re.match(r"^(\d+)", serotype)
    return match.group(1) if match else "-"


def vaccine_membership(
        serotypes: tuple[str, ...]
) -> tuple[list[str], list[str], list[str]]:
    """Which formulations include *serotypes*, and which only cross-react.

    A pair such as ``9A/9V`` is reported as *included* only where every member
    of the pair is in the formulation, and as *partial* where some are: that
    difference is the whole clinical point of not resolving the pair.
    """
    included: list[str] = []
    partial: list[str] = []
    for name, contents in VACCINES.items():
        hits = [s for s in serotypes if s in contents]
        if len(hits) == len(serotypes):
            included.append(name)
        elif hits:
            partial.append(f"{name} (covers {'/'.join(hits)} only)")
    already = set(included) | {entry.split(" ")[0] for entry in partial}
    cross = [f"{name} (cross-reactive with {serotype})"
             for serotype in serotypes
             for name in VACCINE_CROSS_REACTIVE.get(serotype, ())
             if name not in already]
    return included, partial, sorted(set(cross))


# -------------------------------------------------- serotype discrimination
#: Probes whose presence alone names a serotype, and what it names. Transcribed
#: from the ``%singlTarg95`` table of the CDC ``SPN_Serotyper.pl``: a match at
#: the screening thresholds is enough, differences within the probe allowed.
_SINGLE_TARGET: dict[str, str] = {
    "WZY1": "1", "WZY2": "2", "WCHE3": "3", "WZY4": "4", "WZY5": "5",
    "WZY8": "8", "WCRG10A": "10A", "GTF10F": "10F", "WZY11A": "11A",
    "WZY11B:11C": "11B/11C", "WCII12A:12B:46": "12A/12B/46", "WCII12F": "12F",
    "WZY13": "13", "WZY14": "14", "WZY16F": "16F", "GTF17F": "17F",
    "WZY20": "20", "WZY21": "21", "WZY23A": "23A", "WZY23B": "23B",
    "WZY23F": "23F", "WZY24F:24A:24B": "24A/24B/24F", "WCYE25A:25F": "25A/25F",
    "WZY28A": "28A", "WZY31": "31", "TTS": "37", "WZY34": "34",
    "WZY35A": "35A", "WZY35C:42": "35C/42", "WCRO35F": "35F",
    "WHAI47F": "47F", "WCYV38": "38",
}
#: Probes that must match with no differences at all before they name a
#: serotype (``%singlTarg100``), because the serogroup's *wzy* alleles differ by
#: only a few bases and a near match is not evidence.
_SINGLE_TARGET_EXACT: dict[str, str] = {
    "WZY18A": "18A", "WZY18F": "18F", "WZY19A": "19A", "WZY19AVAR": "19A",
    "WZY19F": "19F", "WZY19FVAR": "19F",
}

#: Probes the pilus islets are read from.
_PILUS_PROBES = {"RRGA": "PI-1", "PITB": "PI-2"}

#: What each serotype call means for the pairs that an assembly cannot always
#: separate, keyed by the call this module emits.
_SEROTYPE_CAVEATS: dict[str, str] = {
    "15B/15C": (
        "15B and 15C differ only by the number of TA dinucleotide repeats in "
        "the O-acetyltransferase wciZ - eight in 15B, nine in 15C, which "
        "frameshifts the gene - and the repeat count drifts by slipped-strand "
        "mispairing. Reads show the mixture directly (PneumoCaT routinely "
        "finds 20-80% of reads carrying the frameshift and reports 15B/C); an "
        "assembly collapses the repeat to a single consensus, which is the "
        "variant class short-read assemblers resolve worst. The two are also "
        "very hard to separate with SSI sera on the bench. 15B/C is an "
        "accepted reported serotype; neither 15B nor 15C alone should be "
        "claimed from an assembly"
    ),
    "9A/9V": (
        "9A is 9V with the O-acetyltransferase wcjE inactivated, in 9A by a "
        "single-base deletion at wcjE position 721. The wcjE sequence here "
        "differs from the reference, but this module does not determine "
        "whether the change abolishes function, so both are reported"
    ),
    "18B/18C": (
        "18B is 18C with wciX inactivated by a TGA stop at codon 168-171. The "
        "wciX sequence here differs from the reference without a determination "
        "of whether function is lost. Of this family of pairs this is the most "
        "assembly-tractable, because the lesion is a substitution creating a "
        "stop rather than an indel; whole-locus coverage can never separate "
        "them, since the 18B and 18C reference loci are the same length"
    ),
    "35B/35D": (
        "35D is 35B with the O-acetyltransferase wciG inactivated; the wciG "
        "sequence here differs from the reference without a determination of "
        "whether function is lost"
    ),
    "7A/7F": (
        "7A is 7F with wcwD inactivated by a single-base insertion at position "
        "587. The wcwD sequence here differs from the reference without a "
        "determination of whether function is lost. No pure 7A has yet been "
        "resolved by PneumoCaT either, on the few isolates tried, so a clean "
        "7A call is not currently achievable from sequence"
    ),
    "11B/11C": "the CDC probe set does not separate 11B from 11C",
    "12A/12B/46": "the CDC probe set does not separate 12A, 12B and 46",
    "24A/24B/24F": "no validated rule separates 24A, 24B and 24F from "
                   "sequence; PneumoCaT ships no serogroup 24 variant data and "
                   "SeroBA reports serogroup 24 rather than a serotype",
    "32A/32F": "no nucleotide difference has been found within the coding "
               "regions of the 32A and 32F loci",
    "25A/25F": "the CDC probe set does not separate 25A from 25F (38 has its "
               "own wcyV probe). PneumoCaT does publish six wcyC SNPs for the "
               "split, but PneumoKITy, the current UKHSA screening tool, "
               "dropped them and reports 25F/25A",
    "35C/42": "the CDC probe set does not separate 35C from 42",
}
#: Serotypes whose detection this module cannot distinguish from a near
#: relative at all, because the distinguishing lesion is heterogeneous rather
#: than at a fixed position. Reported as a note on the relative's call.
_SILENT_RELATIVES: dict[str, str] = {
    "11A": "A serotype 11E isolate would be reported here as 11A. 11E is 11A "
           "with wcjE disrupted, but every 11E described carries a different, "
           "independently arisen lesion - IS1515 insertions in either "
           "orientation, deletions of 332 bp and 43 bp, a two-guanine "
           "insertion, a point mutation - so there is no position to query, "
           "and no serotyping tool implements it. O-acetylation is also lost "
           "only partially in some isolates, so the 11A/11E boundary is not "
           "a clean one",
    "33F": "A serotype 37 isolate also carries a 33F-like cps locus; 37 is "
           "called from tts, which lies outside the locus, and is screened "
           "separately here. An isolate lacking wcjE entirely is not 33A",
}


# ---------------------------------------------------------------- PBP typing
@dataclass(frozen=True)
class PbpRegion:
    """One transpeptidase region: where it sits, and what is conserved in it."""

    locus: str          # 1A | 2B | 2X
    gene: str           # pbp1a | pbp2b | pbp2x
    probe: str          # sequence name in pbp_regions.fna
    first_residue: int  # position in the full protein of the region's residue 1
    length: int         # residues in the typed fragment
    #: Active-site motifs, as (motif, absolute position of its first residue).
    motifs: tuple[tuple[str, int], ...]
    #: Substitutions of established beta-lactam relevance, as
    #: (absolute position, wild-type residue, the substitutions reported).
    key_sites: tuple[tuple[int, str, str], ...]


#: The three regions, with the coordinates the CDC allele numbering is defined
#: on. Numbering is on the full-length TIGR4 proteins (GenBank AE005672.3):
#: PBP1a = SP_0369 / UniProt Q04707, 719 residues; PBP2b = SP_1673 / P0A3M5,
#: 680 residues; PBP2x = SP_0336 / P14677, 750 residues. The CDC reference
#: fragments are exact, gap-free substrings of those proteins at 371-647,
#: 379-656 and 229-587, which is the 914-site window the published alignment
#: table uses; the *allele* sets are one residue shorter at each locus
#: (276 / 277 / 358, i.e. 371-646, 379-655 and 229-586), which is what the
#: nucleotide fragments in ``MOD_bLactam_resistance.fasta`` encode (830, 833
#: and 1076 bases = 276, 277 and 358 whole codons) and what a triplet must be
#: read from. Both were confirmed against the modal allele length per locus.
#:
#: The motifs are those the CDC paper defines, after Hakenbeck et al. 2012,
#: each located by substring search in the TIGR4 protein: S370TMK / S466SN /
#: K557TG in PBP1a, S386-V-V-K389 / S443SN / K615TG in PBP2b, and S337TMK /
#: S395SN / K547SG in PBP2x. Note two places where the published text and the
#: protein disagree and the protein wins: PBP2b's SXXK motif is SVVK, not
#: "S386TMK", and its KTG begins at 615, not 614. Note also that PBP1a's
#: catalytic serine S370 sits one residue *before* the typed fragment, so it
#: cannot be read from an allele at all.
#:
#: ``key_sites`` is the CDC pipeline's own empirical list: the 27 positions at
#: which every isolate with a penicillin MIC of 4 mg/L or more differed from
#: PBP type 2-0-2 (Li et al. 2016, Table S1). Positions that additionally have
#: a structural or mutagenesis study behind them carry that citation in the
#: label. Positions with no primary source are deliberately not labelled as
#: resistance mutations, however polymorphic they are in the allele set.
_1A_HOTSPOT = "one of the 27 positions altered in every highly penicillin-" \
              "resistant isolate (Li 2016 Table S1)"
PBP_REGIONS: tuple[PbpRegion, ...] = (
    PbpRegion(
        locus="1A", gene="pbp1a", probe="1A-S2", first_residue=371, length=276,
        motifs=(("STMK", 370), ("SSN", 466), ("KTG", 557)),
        key_sites=(
            (371, "T", ("T371A/T371S, the first variable residue of the STMK "
                       "motif; reverting it lowers the penicillin MIC "
                       "threefold (Smith & Klugman 1998)")),
            (397, "E", _1A_HOTSPOT), (405, "N", _1A_HOTSPOT),
            (540, "T", _1A_HOTSPOT), (546, "N", _1A_HOTSPOT),
            (550, "A", _1A_HOTSPOT),
            (574, "T", ("the TSQF574-577 to NTGY block substitution "
                       "(Smith & Klugman 1998; Contreras-Martel 2006)")),
            (575, "S", "the TSQF574-577 to NTGY block substitution"),
            (576, "Q", "the TSQF574-577 to NTGY block substitution"),
            (577, "F", "the TSQF574-577 to NTGY block substitution"),
            (583, "L", _1A_HOTSPOT), (585, "A", _1A_HOTSPOT),
        ),
    ),
    PbpRegion(
        locus="2B", gene="pbp2b", probe="2B-S2", first_residue=379, length=277,
        motifs=(("SVVK", 386), ("SSN", 443), ("KTG", 615)),
        key_sites=(
            (427, "Q", ("one of the 27 positions altered in every highly "
                       "penicillin-resistant isolate (Li 2016 Table S1)")),
            (446, "T", ("T446A, present in every cefotaxime-resistant strain "
                       "examined and adjacent to the S443SN motif; it cuts "
                       "penicillin G affinity by 60% in vitro "
                       "(Pagliero et al. 2004)")),
            (476, "E", ("E476G, at the active-site cleft "
                       "(Contreras-Martel et al. 2009)")),
            (489, "T", ("one of the 27 positions altered in every highly "
                       "penicillin-resistant isolate (Li 2016 Table S1)")),
        ),
    ),
    PbpRegion(
        locus="2X", gene="pbp2x", probe="2X-S2", first_residue=229, length=358,
        motifs=(("STMK", 337), ("SSN", 395), ("KSG", 547)),
        key_sites=(
            (254, "R", ("one of the 27 positions altered in every highly "
                       "penicillin-resistant isolate (Li 2016 Table S1)")),
            (256, "M", "one of the 27 positions (Li 2016 Table S1)"),
            (338, "T", ("T338A/T338P, in the buried cavity beside the catalytic "
                       "serine; directly implicated in resistance "
                       "(Gordon et al. 2000; Mouz et al. 1998) and the "
                       "single highest-impact reversion in a resistant "
                       "clinical isolate (Smith & Klugman 2005)")),
            (339, "M", ("M339F, which further lowers beta-lactam "
                       "susceptibility in resistant strains "
                       "(Chesnel et al. 2003)")),
            (371, "I", ("I371T, the other highest-impact reversion in a "
                       "resistant clinical isolate (Smith & Klugman 2005)")),
            (382, "G", "one of the 27 positions (Li 2016 Table S1)"),
            (384, "R", ("R384G, one of six substitutions whose collective "
                       "reversion restored susceptibility "
                       "(Smith & Klugman 2005)")),
            (401, "T", "one of the 27 positions (Li 2016 Table S1)"),
            (444, "N", "one of the 27 positions (Li 2016 Table S1)"),
            (531, "S", "one of the 27 positions (Li 2016 Table S1)"),
            (550, "T", ("T550A, directly related to resistance in the crystal "
                       "structure (Gordon et al. 2000)")),
            (552, "Q", ("Q552E, directly related to resistance in the crystal "
                       "structure (Gordon et al. 2000)")),
            (565, "V", "one of the 27 positions (Li 2016 Table S1)"),
            (576, "S", "one of the 27 positions (Li 2016 Table S1)"),
        ),
    ),
)
_REGION_BY_LOCUS = {region.locus: region for region in PBP_REGIONS}

#: The MIC table's column order, with the drug each column is.
PBP_MIC_DRUGS: tuple[tuple[str, str], ...] = (
    ("PEN", "penicillin"), ("AMO", "amoxicillin"), ("MER", "meropenem"),
    ("TAX", "cefotaxime"), ("CFT", "ceftriaxone"), ("CFX", "cefuroxime"),
)

#: CLSI interpretive criteria, transcribed from the log2 thresholds in
#: ``MIC_format_with_SIR.R`` of the CDC pipeline, which is the only place the
#: CDC publishes both the meningitis and the non-meningitis sets - the two
#: papers each print only one of them, and they print different ones. Each
#: entry gives the breakpoints as log2(MIC in mg/L): susceptible at or below
#: the first, resistant at or above the second, intermediate between. Where the
#: two are adjacent there is no intermediate category, which is correct for
#: penicillin under the meningitis criterion.
_CLSI: dict[str, tuple[float, float, str]] = {
    # Penicillin has two sets of criteria; both are reported, because which one
    # applies is a clinical question about the site of infection, not a
    # laboratory one, and the module cannot know the answer.
    "PEN/meningitis": (-4.0, -3.0, "penicillin, meningitis (S <=0.06, R >=0.12 mg/L)"),
    "PEN/non-meningitis": (1.0, 3.0, "penicillin, non-meningitis parenteral "
                                     "(S <=2, I 4, R >=8 mg/L)"),
    "AMO": (1.0, 3.0, "amoxicillin (S <=2, I 4, R >=8 mg/L)"),
    "MER": (-2.0, 0.0, "meropenem (S <=0.25, I 0.5, R >=1 mg/L)"),
    "TAX/meningitis": (-1.0, 1.0, "cefotaxime, meningitis (S <=0.5, I 1, R >=2 mg/L)"),
    "TAX/non-meningitis": (0.0, 2.0, "cefotaxime, non-meningitis "
                                     "(S <=1, I 2, R >=4 mg/L)"),
    "CFT/meningitis": (-1.0, 1.0, "ceftriaxone, meningitis (S <=0.5, I 1, R >=2 mg/L)"),
    "CFT/non-meningitis": (0.0, 2.0, "ceftriaxone, non-meningitis "
                                     "(S <=1, I 2, R >=4 mg/L)"),
    "CFX": (-1.0, 1.0, "cefuroxime (S <=0.5, I 1, R >=2 mg/L)"),
}

#: Published agreement of the CDC PBP-type MIC prediction with broth
#: microdilution, quoted in the note that carries every prediction.
#:
#: Neither paper publishes per-antibiotic agreement as numbers - in both they
#: are bar charts with no source-data table - so what is quoted here is what
#: the papers state in words, for all six beta-lactams together, and the two
#: per-antibiotic extremes from the one validation paragraph that names them.
#: Nothing here is interpolated from a figure.
PBP_MIC_ACCURACY = (
    "As published: for isolates whose PBP type was in the training set, the "
    "mode-MIC model gave essential agreement (within one doubling dilution of "
    "the measured MIC) above 98%, category agreement above 94%, a major "
    "discrepancy rate below 3% and a very major discrepancy rate below 2%, for "
    "all six beta-lactams across 2,528 isolates (Li et al. 2016). On an "
    "independent set of 1,781 isolates, essential agreement ranged from 97.6% "
    "(penicillin, mode-MIC) to 99.8% (ceftriaxone) and category agreement from "
    "95.7% (penicillin, mode-MIC) to 99.1% (ceftriaxone) (Li et al. 2017). PBP "
    "type explained 94-99% of the variance in log2 MIC, against 73-91% for the "
    "MLST sequence type."
)


# ------------------------------------------------------ macrolide and friends
#: Probe -> (determinant, mechanism, what it predicts). The phenotype each
#: mechanism implies is the clinical point: Hydra's AMR screen already reports
#: these as genes, but reports them without saying what they mean for therapy.
_MACROLIDE: dict[str, tuple[str, str, str]] = {
    "ERMB": (
        "erm(B)", "23S rRNA A2058 dimethylation (target modification)",
        "MLS-B phenotype: erythromycin, clarithromycin and azithromycin "
        "resistant, and clindamycin resistant as well. Macrolides and "
        "clindamycin are both unusable. erm(B) is usually constitutive in "
        "pneumococcus but can be inducible; the distinction is not visible in "
        "the sequence of the gene alone.",
    ),
    "ERMTR": (
        "erm(TR) (erm(A) subclass TR)", "23S rRNA A2058 dimethylation",
        "inducible MLS-B phenotype: erythromycin resistant. Clindamycin may "
        "test susceptible while being clinically unreliable, because the "
        "methylase is induced by macrolide exposure; clindamycin should not be "
        "used on the basis of a susceptible result alone.",
    ),
    "MEF": (
        "mef", "ATP-independent efflux of 14- and 15-membered macrolides",
        "M phenotype: erythromycin, clarithromycin and azithromycin resistant, "
        "typically at lower MICs than erm(B) gives. Clindamycin and the "
        "16-membered macrolides remain active. This probe does not separate "
        "mef(A) (Tn1207.1) from mef(E) (the Mega element) - the two give the "
        "same phenotype and differ epidemiologically; Hydra's AMR gene screen "
        "names the specific allele, and msr(D)/mel, where its reference "
        "database distinguishes them.",
    ),
    "LSAC": (
        "lsa(C)", "ABC-F ribosomal protection",
        "LS-A phenotype: clindamycin (and lincomycin, streptogramin A, "
        "pleuromutilin) resistant with erythromycin unaffected. Clindamycin "
        "resistance without macrolide resistance.",
    ),
}
#: Probes that qualify an erm(B) call rather than standing on their own.
_ERMB_CONTEXT = {
    "ERMBUPST": "erm(B) upstream/leader region present and divergent from the "
                "reference",
    "ERMBS": "erm(B) leader region identical to the reference",
}
#: Other determinants the CDC probe set screens, with what they predict.
_OTHER_DETERMINANTS: dict[str, tuple[str, str]] = {
    "TETM": ("tet(M)", "tetracycline and doxycycline resistance (ribosomal "
                       "protection); usually carried on Tn916-family elements "
                       "that also carry erm(B)"),
    "TETO": ("tet(O)", "tetracycline resistance (ribosomal protection)"),
    "CAT": ("cat", "chloramphenicol resistance (acetyltransferase)"),
}

#: Target-gene regions whose substitutions are called, with the reference
#: peptide, the offset that turns a position in that peptide into a position in
#: the protein, the polymorphisms the CDC pipeline treats as neutral (given in
#: peptide-local numbering), and what the region means.
@dataclass(frozen=True)
class TargetRegion:
    probe: str
    gene: str
    peptide: str
    offset: int
    neutral: frozenset[str]
    drug: str
    meaning: str


TARGET_REGIONS: tuple[TargetRegion, ...] = (
    TargetRegion(
        probe="GYRA", gene="gyrA",
        peptide="VMGKYHPHGDSSIYEAMVRMAQWWSY", offset=70, neutral=frozenset(),
        drug="fluoroquinolone",
        meaning="gyrA quinolone-resistance-determining region (residues 71-96). "
                "A substitution at S81 (S81F or S81Y) is the first-step "
                "mutation that raises levofloxacin and moxifloxacin MICs; with "
                "a parC substitution as well the isolate is usually fully "
                "fluoroquinolone resistant.",
    ),
    TargetRegion(
        probe="PARC", gene="parC",
        peptide=(("DSSIYDAMVRMSQNWKNREILVEMHGNNGSMDGDPPAAMRYTEARLSEIA"
                 "GYLLQDIEKKTVPFAWNFDD")),
        offset=77, neutral=frozenset({"G51D", "K60N", "R18C"}),
        drug="fluoroquinolone",
        meaning="parC quinolone-resistance-determining region (residues 78-147). "
                "S79F/S79Y and D83N/D83Y are the established first-step "
                "substitutions; alone they raise ciprofloxacin and "
                "levofloxacin MICs but often leave moxifloxacin active.",
    ),
    TargetRegion(
        probe="FOLA", gene="folA",
        peptide="QDVQSVLDWYQDQEKNLYII", offset=80,
        neutral=frozenset({"D12A", "D12V", "D8G"}),
        drug="trimethoprim-sulfamethoxazole",
        meaning="folA (dihydrofolate reductase) residues 81-100. I100L is the "
                "substitution that confers trimethoprim resistance. "
                "Sulfamethoxazole resistance additionally requires a 1-2 codon "
                "insertion in folP, which this module does not call.",
    ),
    TargetRegion(
        probe="RPLD1", gene="rplD",
        peptide="KPWRQKGTGRAR", offset=62, neutral=frozenset(),
        drug="macrolide",
        meaning="ribosomal protein L4 residues 63-74. Substitutions in this "
                "motif raise macrolide MICs without erm or mef, and are the "
                "commonest cause of macrolide resistance in isolates negative "
                "for both.",
    ),
)
_TARGET_BY_PROBE = {region.probe: region for region in TARGET_REGIONS}


# ------------------------------------------------------------- bibliography
REFERENCES: tuple[Reference, ...] = (
    Reference(
        key="Bentley2006",
        authors="Bentley SD, Aanensen DM, Mavroidi A, et al.",
        year=2006,
        title="Genetic analysis of the capsular biosynthetic locus from all 90 "
              "pneumococcal serotypes",
        journal="PLoS Genetics 2:e31",
        doi="10.1371/journal.pgen.0020031", pmid="16532061",
        supports="the cps locus reference set the serotype call's first stage "
                 "measures coverage against",
    ),
    Reference(
        key="Geno2015",
        authors="Geno KA, Gilbert GL, Song JY, et al.",
        year=2015,
        title="Pneumococcal capsules and their types: past, present, and future",
        journal="Clinical Microbiology Reviews 28:871-899",
        doi="10.1128/CMR.00024-15", pmid="26085553",
        supports="the serotype/serogroup nomenclature, and which within-serogroup "
                 "pairs differ only by O-acetylation or by one glycosyltransferase",
    ),
    Reference(
        key="Metcalf2016Sero",
        authors="Metcalf BJ, Gertz RE Jr, Gladstone RA, et al.",
        year=2016,
        title="Strain features and distributions in pneumococci from children "
              "with invasive disease before and after 13-valent conjugate "
              "vaccine implementation in the USA",
        journal="Clinical Microbiology and Infection 22:60.e9-60.e29",
        doi="10.1016/j.cmi.2015.08.027", pmid="26363404",
        supports="the CDC whole-genome serotyping probe set, the pipeline it is "
                 "used in, and the extraction of the PBP transpeptidase regions",
    ),
    Reference(
        key="Metcalf2016Res",
        authors="Metcalf BJ, Chochua S, Gertz RE Jr, et al.",
        year=2016,
        title="Using whole genome sequencing to identify resistance "
              "determinants and predict antimicrobial resistance phenotypes "
              "for year 2015 invasive pneumococcal disease isolates recovered "
              "in the United States",
        journal="Clinical Microbiology and Infection 22:1002.e1-1002.e8",
        doi="10.1016/j.cmi.2016.08.001", pmid="27542334",
        supports="the macrolide, lincosamide, tetracycline, fluoroquinolone and "
                 "folate-pathway determinant set, the reference peptides and "
                 "offsets the substitutions are numbered on, and the "
                 "determinant-to-phenotype interpretation",
    ),
    Reference(
        key="Kapatai2016",
        authors="Kapatai G, Sheppard CL, Al-Shahib A, et al.",
        year=2016,
        title="Whole genome sequencing of Streptococcus pneumoniae: "
              "development, evaluation and verification of targets for "
              "serogroup and serotype prediction using an automated pipeline",
        journal="PeerJ 4:e2477",
        doi="10.7717/peerj.2477", pmid="27672516",
        supports="which serotype distinctions sequence can and cannot make, "
                 "including 32A/32F, serogroup 24, 15B/15C and 7A/7F",
    ),
    Reference(
        key="Epping2018",
        authors="Epping L, van Tonder AJ, Gladstone RA, et al.",
        year=2018,
        title="SeroBA: rapid high-throughput serotyping of Streptococcus "
              "pneumoniae from whole genome sequence data",
        journal="Microbial Genomics 4:e000186",
        doi="10.1099/mgen.0.000186", pmid="29870330",
        supports="the serogroup 24 limitation, and the reads-based serotyping "
                 "tool this module is the assembly-only counterpart of",
    ),
    Reference(
        key="Li2016",
        authors="Li Y, Metcalf BJ, Chochua S, et al.",
        year=2016,
        title="Penicillin-binding protein transpeptidase signatures for "
              "tracking and predicting beta-lactam resistance levels in "
              "Streptococcus pneumoniae",
        journal="mBio 7:e00756-16",
        doi="10.1128/mBio.00756-16", pmid="27302760",
        supports="the pbp1a/pbp2b/pbp2x transpeptidase allele triplet, its "
                 "numbering, and the PBP-type-to-MIC reference table",
    ),
    Reference(
        key="Li2017",
        authors="Li Y, Metcalf BJ, Chochua S, et al.",
        year=2017,
        title="Validation of beta-lactam minimum inhibitory concentration "
              "predictions for pneumococcal isolates with newly encountered "
              "penicillin binding protein (PBP) sequences",
        journal="BMC Genomics 18:621",
        doi="10.1186/s12864-017-4017-7", pmid="28810827",
        supports="the behaviour of the prediction for PBP types absent from the "
                 "reference table, which is why this module reports such a "
                 "triplet without an MIC rather than extrapolating",
    ),
    Reference(
        key="Mavroidi2004",
        authors="Mavroidi A, Godoy D, Aanensen DM, Robinson DA, "
                "Hollingshead SK, Spratt BG",
        year=2004,
        title="Evolutionary genetics of the capsular locus of serogroup 6 "
              "pneumococci",
        journal="Journal of Bacteriology 186:8181-8192",
        doi="10.1128/JB.186.24.8181-8192.2004", pmid="15576766",
        supports="wciP as the 6A/6B determinant (Ser195 in 6A/6C, Asn195 in "
                 "6B/6D) and the serogroup 6 locus structure",
    ),
    Reference(
        key="Park2007Discovery",
        authors="Park IH, Pritchard DG, Cartee R, Brandao A, Brandileone MC, "
                "Nahm MH",
        year=2007,
        title="Discovery of a new capsular serotype (6C) within serogroup 6 of "
              "Streptococcus pneumoniae",
        journal="Journal of Clinical Microbiology 45:1225-1233",
        doi="10.1128/JCM.02199-06", pmid="17267625",
        supports="serotype 6C as a serotype",
    ),
    Reference(
        key="Park2007",
        authors="Park IH, Park S, Hollingshead SK, Nahm MH",
        year=2007,
        title="Genetic basis for the new pneumococcal serotype, 6C",
        journal="Infection and Immunity 75:4482-4489",
        doi="10.1128/IAI.00510-07", pmid="17576753",
        supports="wciN as the 6A/6C determinant, and therefore the four-way "
                 "6A/6B/6C/6D rule",
    ),
    Reference(
        key="Bratcher2010",
        authors="Bratcher PE, Kim KH, Kang JH, Hong JY, Nahm MH",
        year=2010,
        title="Identification of natural pneumococcal isolates expressing "
              "serotype 6D by genetic, biochemical and serological "
              "characterization",
        journal="Microbiology 156:555-560",
        doi="10.1099/mic.0.034116-0", pmid="19942663",
        supports="serotype 6D, completing the wciP x wciN rule",
    ),
    Reference(
        key="vanSelm2003",
        authors="van Selm S, van Cann LM, Kolkman MA, van der Zeijst BA, "
                "van Putten JP",
        year=2003,
        title="Genetic basis for the structural difference between "
              "Streptococcus pneumoniae serotype 15B and 15C capsular "
              "polysaccharides",
        journal="Infection and Immunity 71:6192-6198",
        doi="10.1128/IAI.71.11.6192-6198.2003", pmid="14573636",
        supports="the wciZ TA-repeat count (eight units in 15B, nine and a "
                 "frameshift in 15C) and therefore why 15B and 15C are "
                 "reported as a pair. The paper calls the gene cps15bM/cps15cM; "
                 "the name wciZ comes from Bentley et al. 2006",
    ),
    Reference(
        key="Calix2010",
        authors="Calix JJ, Nahm MH",
        year=2010,
        title="A new pneumococcal serotype, 11E, has a variably inactivated "
              "wcjE gene",
        journal="Journal of Infectious Diseases 202:29-38",
        doi="10.1086/653123", pmid="20507232",
        supports="wcjE inactivation as the mechanism behind the O-acetylation "
                 "pairs (11A/11E, 9V/9A, 33F/33A)",
    ),
    Reference(
        key="Calix2012",
        authors="Calix JJ, Porambo RJ, Brady AM, et al.",
        year=2012,
        title="Biochemical, genetic, and serological characterization of two "
              "capsule subtypes among Streptococcus pneumoniae serotype 20 "
              "strains: discovery of a new pneumococcal serotype",
        journal="Journal of Biological Chemistry 287:27885-27894",
        doi="10.1074/jbc.M112.380451", pmid="22736767",
        supports="the split of Danish type 20 into 20A and 20B, and therefore "
                 "how a PPSV23 'type 20' membership is reported",
    ),
    Reference(
        key="Barocchi2006",
        authors="Barocchi MA, Ries J, Zogaj X, et al.",
        year=2006,
        title="A pneumococcal pilus influences virulence and host inflammatory "
              "responses",
        journal="Proceedings of the National Academy of Sciences USA "
                "103:2857-2862",
        doi="10.1073/pnas.0511017103", pmid="16481624",
        supports="pilus islet 1 (the rlrA islet, rrgA-rrgC) and its clinical "
                 "relevance",
    ),
    Reference(
        key="Bagnoli2008",
        authors="Bagnoli F, Moschioni M, Donati C, et al.",
        year=2008,
        title="A second pilus type in Streptococcus pneumoniae is prevalent in "
              "emerging serotypes and mediates adhesion to host cells",
        journal="Journal of Bacteriology 190:5480-5492",
        doi="10.1128/JB.00384-08", pmid="18515415",
        supports="pilus islet 2 (pitA-pitB-sipA-srtG1) and the pitB marker",
    ),
    Reference(
        key="Cooper2011",
        authors="Cooper D, Yu X, Sidhu M, Nahm MH, Fernsten P, Jansen KU",
        year=2011,
        title="The 13-valent pneumococcal conjugate vaccine (PCV13) elicits "
              "cross-functional opsonophagocytic killing responses in humans "
              "to Streptococcus pneumoniae serotypes 6C and 7A",
        journal="Vaccine 29:7207-7211",
        doi="10.1016/j.vaccine.2011.06.056", pmid="21689707",
        supports="the 6A-to-6C and 7F-to-7A cross-reactivity reported alongside "
                 "vaccine membership",
    ),
    Reference(
        key="Hakenbeck2012",
        authors="Hakenbeck R, Bruckner R, Denapaite D, Maurer P",
        year=2012,
        title="Molecular mechanisms of beta-lactam resistance in Streptococcus "
              "pneumoniae",
        journal="Future Microbiology 7:395-410",
        doi="10.2217/fmb.12.2", pmid="22393892",
        supports="the definition of the SXXK, SXN and KTG/KSG active-site "
                 "motifs and their positions, which the CDC scheme adopts",
    ),
    Reference(
        key="Smith1998",
        authors="Smith AM, Klugman KP",
        year=1998,
        title="Alterations in PBP 1A essential for high-level penicillin "
              "resistance in Streptococcus pneumoniae",
        journal="Antimicrobial Agents and Chemotherapy 42:1329-1333",
        doi="10.1128/AAC.42.6.1329", pmid="9624469",
        supports="pbp1a residue 371 and the TSQF574-577 block as positions "
                 "that determine the penicillin MIC",
    ),
    Reference(
        key="ContrerasMartel2006",
        authors="Contreras-Martel C, Job V, Di Guilmi AM, Vernet T, "
                "Dideberg O, Dessen A",
        year=2006,
        title="Crystal structure of penicillin-binding protein 1a (PBP1a) "
              "reveals a mutational hotspot implicated in beta-lactam "
              "resistance in Streptococcus pneumoniae",
        journal="Journal of Molecular Biology 355:684-696",
        doi="10.1016/j.jmb.2005.10.030", pmid="16316661",
        supports="the structural basis of the pbp1a T371 and TSQF574-577 "
                 "hotspot",
    ),
    Reference(
        key="Pagliero2004",
        authors="Pagliero E, Chesnel L, Hopkins J, et al.",
        year=2004,
        title="Biochemical characterization of Streptococcus pneumoniae "
              "penicillin-binding protein 2b and its implication in "
              "beta-lactam resistance",
        journal="Antimicrobial Agents and Chemotherapy 48:1848-1855",
        doi="10.1128/AAC.48.5.1848-1855.2004", pmid="15105143",
        supports="pbp2b T446A, its proximity to the S443SN motif and its 60% "
                 "reduction in penicillin G affinity",
    ),
    Reference(
        key="ContrerasMartel2009",
        authors="Contreras-Martel C, Dahout-Gonzalez C, Martins A dos S, "
                "Kotnik M, Dessen A",
        year=2009,
        title="PBP active site flexibility as the key mechanism for "
              "beta-lactam resistance in pneumococci",
        journal="Journal of Molecular Biology 387:899-909",
        doi="10.1016/j.jmb.2009.02.024", pmid="19233207",
        supports="the structural role of pbp2b T446 and E476",
    ),
    Reference(
        key="Gordon2000",
        authors="Gordon E, Mouz N, Duee E, Dideberg O",
        year=2000,
        title="The crystal structure of the penicillin-binding protein 2x from "
              "Streptococcus pneumoniae and its acyl-enzyme form: implication "
              "in drug resistance",
        journal="Journal of Molecular Biology 299:477-485",
        doi="10.1006/jmbi.2000.3740", pmid="10860753",
        supports="pbp2x T338, T550 and Q552 as the residues directly related "
                 "to resistance",
    ),
    Reference(
        key="Mouz1998",
        authors="Mouz N, Gordon E, Di Guilmi AM, et al.",
        year=1998,
        title="Identification of a structural determinant for resistance to "
              "beta-lactam antibiotics in Gram-positive bacteria",
        journal="Proceedings of the National Academy of Sciences USA "
                "95:13403-13406",
        doi="10.1073/pnas.95.23.13403", pmid="9811812",
        supports="the buried cavity beside the pbp2x catalytic serine that "
                 "T338 lines",
    ),
    Reference(
        key="Chesnel2003",
        authors="Chesnel L, Pernot L, Lemaire D, et al.",
        year=2003,
        title="The structural modifications induced by the M339F substitution "
              "in PBP2x from Streptococcus pneumoniae further decreases the "
              "susceptibility to beta-lactams of resistant strains",
        journal="Journal of Biological Chemistry 278:44448-44456",
        doi="10.1074/jbc.M305948200", pmid="12923202",
        supports="pbp2x M339F",
    ),
    Reference(
        key="Smith2005",
        authors="Smith AM, Klugman KP",
        year=2005,
        title="Amino acid mutations essential to production of an altered PBP "
              "2X conferring high-level beta-lactam resistance in a clinical "
              "isolate of Streptococcus pneumoniae",
        journal="Antimicrobial Agents and Chemotherapy 49:4622-4627",
        doi="10.1128/AAC.49.11.4622-4627.2005", pmid="16251304",
        supports="pbp2x T338P, I371T and R384G as the substitutions whose "
                 "reversion restores susceptibility",
    ),
    Reference(
        key="Park2012",
        authors="Park IH, Kim KH, Andrade AL, Briles DE, McDaniel LS, Nahm MH",
        year=2012,
        title="Nontypeable pneumococci can be divided into multiple cps types, "
              "including one type expressing the novel gene pspK",
        journal="mBio 3:e00035-12",
        doi="10.1128/mBio.00035-12", pmid="22532557",
        supports="the null-capsule clade definitions quoted when a "
                 "non-encapsulated genotype is reported: NCC1 is pspK-positive "
                 "and aliC/aliD-negative, NCC2 carries aliC and aliD, NCC3 "
                 "carries aliD only",
    ),
    Reference(
        key="Keller2016",
        authors="Keller LE, Robinson DA, McDaniel LS",
        year=2016,
        title="Nonencapsulated Streptococcus pneumoniae: emergence and "
              "pathogenesis",
        journal="mBio 7:e01792-15",
        doi="10.1128/mBio.01792-15", pmid="27006456",
        supports="what a non-encapsulated genotype means clinically",
    ),
    Reference(
        key="CLSI",
        authors="Clinical and Laboratory Standards Institute",
        year=2013,
        title="Performance standards for antimicrobial susceptibility testing; "
              "twenty-third informational supplement (M100-S23)",
        journal="CLSI, Wayne PA",
        supports="the penicillin, amoxicillin, meropenem, cefotaxime, "
                 "ceftriaxone and cefuroxime breakpoints, including the "
                 "separate meningitis criteria, that the predicted MIC is "
                 "interpreted against; this is the edition both CDC papers "
                 "cite, and the thresholds are unchanged in current M100",
    ),
    Reference(
        key="Gladstone2019",
        authors="Gladstone RA, Lo SW, Lees JA, et al.",
        year=2019,
        title="International genomic definition of pneumococcal lineages, to "
              "contextualise disease, antibiotic resistance and vaccine impact",
        journal="EBioMedicine 43:338-346",
        doi="10.1016/j.ebiom.2019.04.021", pmid="31003929",
        supports="the GPSC lineage scheme, which this module explains it does "
                 "not assign and why",
    ),
)


# ---------------------------------------------------------------------- fetch
def _download(url: str, dest: Path) -> None:
    """Fetch *url* to *dest*, retrying transport failures a few times."""
    import time  # noqa: PLC0415 - only a download needs it
    import urllib.error  # noqa: PLC0415
    import urllib.request  # noqa: PLC0415

    dest.parent.mkdir(parents=True, exist_ok=True)
    agent = "hydra-amr (+https://github.com/iowa69/hydra)"
    last: Exception | None = None
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": agent})
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = response.read()
            if not payload:
                raise HydraError(f"{url} returned an empty file")
            dest.write_bytes(payload)
            return
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code not in (429, 500, 502, 503, 504):
                break
        except (urllib.error.URLError, OSError) as exc:
            last = exc
        time.sleep(1.5 * (attempt + 1))
    raise HydraError(f"could not download {url}: {last}")


def fetch(db_dir: Path) -> None:
    """Download and index this module's data into ``db_dir/spneumoniae/``.

    Idempotent: a file already present with a BLAST index newer than itself is
    left alone, so re-running costs nothing and a half-finished download is
    completed rather than repeated from the start.
    """
    from ...db.manager import makeblastdb  # noqa: PLC0415 - avoids a cycle

    root = Path(db_dir) / KEY
    root.mkdir(parents=True, exist_ok=True)
    for url, relative, dbtype, _db in _FILES:
        dest = root / relative
        if not dest.exists() or dest.stat().st_size == 0:
            LOG.info("%s: downloading %s", KEY, relative)
            _download(url, dest)
        if dbtype:
            makeblastdb(dest, dbtype)
    _write_provenance(root)


def _write_provenance(root: Path) -> None:
    """Record where every file came from, beside the files themselves.

    The integrator decides from this whether a database may be shipped inside
    the Windows zip, so it lives with the data rather than only in the docs.
    """
    by_name = {db.name: db for db in DATABASES}
    lines = ["file\tsource\tversion\tlicence\tredistribution"]
    for url, relative, _dbtype, name in _FILES:
        db = by_name[name]
        flat = " ".join(db.licence.split())
        lines.append(f"{relative}\t{url}\t{db.version}\t{flat}"
                     f"\t{db.redistribution}")
    (root / "SOURCES.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")


# -------------------------------------------------------------------- helpers
@dataclass
class ProbeHit:
    """One screening probe, as found in the assembly."""

    cluster: str
    gene: str
    allele: str
    identity: float
    coverage: float
    contig: str
    start: int
    end: int
    strand: str

    @property
    def present(self) -> bool:
        return (self.coverage >= PROBE_MIN_COVERAGE
                and self.identity >= PROBE_MIN_IDENTITY)

    @property
    def exact(self) -> bool:
        """No differences at all: the CDC rules' "diffs" column being empty."""
        return self.identity >= EXACT_IDENTITY and self.coverage >= EXACT_COVERAGE

    @property
    def evidence(self) -> str:
        mark = "exact" if self.exact else f"{self.identity:.1f}% identity"
        return f"{mark}, {self.coverage:.1f}% of {self.allele}"


def _parse_srst2_header(header: str) -> tuple[str, str, str]:
    """``>30__WCIP6AC__WCIP6AC-1__37`` -> (cluster, gene, allele)."""
    parts = header.split()[0].split("__")
    if len(parts) >= 3:
        return parts[0], parts[1], parts[2]
    return "", header.split()[0], header.split()[0]


def _probe_hits(assembly: Path, fasta: Path, out: Path,
                threads: int) -> dict[str, ProbeHit]:
    """Screen the assembly for a probe set, keeping the best allele per cluster.

    The CDC rules read one gene per probe cluster, because that is what srst2
    reports: serogroup 19's four *wzy* alleles sit in one cluster and only the
    closest is considered. Keeping them all would let a genome satisfy the 19A
    and the 19F rule at once.
    """
    check_db_exists(fasta, "nucl")
    # Probes run from 35 to 1479 bases. The default e-value would discard a
    # perfect 35-mer, and megablast's word size cannot seed one at all.
    hsps = blast("blastn", assembly, fasta, out, threads=threads, evalue=1e-5,
                 task="blastn", word_size=11, perc_identity=PROBE_MIN_IDENTITY - 5.0,
                 max_target_seqs=10000)
    best_per_cluster: dict[str, ProbeHit] = {}
    for hit in merge_hsps(hsps):
        cluster, gene, allele = _parse_srst2_header(hit.sseqid)
        candidate = ProbeHit(
            cluster=cluster, gene=gene, allele=allele,
            identity=hit.identity_pct, coverage=hit.coverage_pct,
            contig=hit.qseqid, start=hit.query_start, end=hit.query_end,
            strand=hit.strand,
        )
        if not candidate.present:
            continue
        current = best_per_cluster.get(cluster or gene)
        if current is None or _better_probe(candidate, current):
            best_per_cluster[cluster or gene] = candidate
    return {hit.gene: hit for hit in best_per_cluster.values()}


def _better_probe(candidate: ProbeHit, current: ProbeHit) -> bool:
    return ((candidate.identity, candidate.coverage, candidate.allele)
            > (current.identity, current.coverage, current.allele))


@dataclass
class LocusMatch:
    """How much of one reference *cps* locus the assembly covers."""

    serotype: str
    length: int
    coverage: float
    identity: float
    intervals: list[tuple[int, int]] = field(default_factory=list)

    @property
    def score(self) -> float:
        return self.coverage * self.identity / 100.0

    @property
    def evidence(self) -> str:
        spans = ",".join(f"{lo}-{hi}" for lo, hi in self.intervals[:4])
        if len(self.intervals) > 4:
            spans += ",..."
        return (f"{self.coverage:.1f}% of the {self.serotype} locus "
                f"({spans}/{self.length}) at {self.identity:.1f}% identity")


def _cps_loci(assembly: Path, fasta: Path, out: Path, threads: int) -> list[LocusMatch]:
    """Rank the reference *cps* loci by how completely the assembly covers them.

    Coverage is summed over every contig, not taken from the best contig: a
    draft assembly routinely breaks a 20 kb capsule locus across two or three
    of them, and the locus is no less present for that.
    """
    check_db_exists(fasta, "nucl")
    hsps = blast("blastn", assembly, fasta, out, threads=threads, evalue=1e-20,
                 task="megablast", perc_identity=CPS_MIN_IDENTITY,
                 max_target_seqs=10000)
    spans: dict[str, list[tuple[int, int, int, int]]] = {}
    lengths: dict[str, int] = {}
    for hsp in hsps:
        serotype = hsp.sseqid.split()[0]
        lengths[serotype] = hsp.slen
        spans.setdefault(serotype, []).append(
            (hsp.sub_lo, hsp.sub_hi, hsp.nident, hsp.length))
    out_matches: list[LocusMatch] = []
    for serotype, hits in spans.items():
        length = lengths[serotype]
        intervals = _merge((lo, hi) for lo, hi, _n, _l in hits)
        covered = sum(hi - lo + 1 for lo, hi in intervals)
        # Weight identity by alignment length, so one short spurious HSP cannot
        # drag the figure for a locus that is otherwise a clean match.
        total_len = sum(length for _lo, _hi, _n, length in hits) or 1
        nident = sum(n for _lo, _hi, n, _l in hits)
        out_matches.append(LocusMatch(
            serotype=serotype, length=length,
            coverage=100.0 * covered / length if length else 0.0,
            identity=100.0 * nident / total_len, intervals=intervals))
    out_matches.sort(key=lambda m: (m.score, m.serotype), reverse=True)
    return out_matches


def _gaps(intervals: list[tuple[int, int]], length: int) -> list[tuple[int, int]]:
    """The stretches of a reference that *intervals* leave uncovered."""
    out: list[tuple[int, int]] = []
    position = 1
    for lo, hi in intervals:
        if lo > position:
            out.append((position, lo - 1))
        position = max(position, hi + 1)
    if position <= length:
        out.append((position, length))
    return sorted(out, key=lambda span: span[0] - span[1])


def _merge(intervals) -> list[tuple[int, int]]:
    ordered = sorted((min(a, b), max(a, b)) for a, b in intervals)
    if not ordered:
        return []
    merged = [list(ordered[0])]
    for lo, hi in ordered[1:]:
        if lo <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [(lo, hi) for lo, hi in merged]


@dataclass
class Extracted:
    """A reference-defined region recovered from the assembly."""

    probe: str
    contig: str
    start: int
    end: int
    strand: str
    nucleotides: str
    coverage: float
    identity: float
    complete: bool
    note: str = ""


def _extract_regions(assembly: Path, fasta: Path, out: Path, threads: int,
                     wanted: set[str], min_coverage: float,
                     min_identity: float) -> dict[str, Extracted]:
    """Recover the full span of each named reference region from the assembly.

    BLAST finds where the region sits; the ends are then extended by however
    much of the reference the alignment did not reach, so a region whose first
    or last few bases diverge is still extracted at full length and in frame.
    This is what the CDC pipeline's ``ExtractGene.pl`` does, and it is why an
    allele carrying a substitution in its first codon is still typed.
    """
    check_db_exists(fasta, "nucl")
    # The e-value has to be loose. A 54-base alignment at 83% identity - which
    # is what a mosaic folA gives against this 60-base reference - scores about
    # 4e-10 against a 2 Mb genome, so a 1e-10 cut-off silently discards the one
    # hit that matters. Identity, coverage and the divergence check downstream
    # are what reject a spurious match; the e-value only has to let it through.
    hsps = blast("blastn", assembly, fasta, out, threads=threads, evalue=1e-3,
                 task="blastn", word_size=11, perc_identity=min_identity - 10.0,
                 max_target_seqs=10000)
    best: dict[str, tuple[float, object]] = {}
    for hit in merge_hsps(hsps):
        name = hit.sseqid.split()[0]
        probe = name.split("__")[1] if "__" in name else name
        if probe not in wanted:
            continue
        score = hit.coverage_pct * hit.identity_pct
        if probe not in best or score > best[probe][0]:
            best[probe] = (score, hit)

    contigs: dict[str, str] = {}
    for header, seq in read_fasta(assembly):
        contigs[header.split()[0]] = seq

    found: dict[str, Extracted] = {}
    for probe, (_score, hit) in best.items():
        if hit.coverage_pct < min_coverage or hit.identity_pct < min_identity:
            continue
        contig = contigs.get(hit.qseqid)
        if contig is None:
            continue
        ref_lo = min(lo for lo, _hi in hit.subject_intervals)
        ref_hi = max(hi for _lo, hi in hit.subject_intervals)
        missing_left, missing_right = ref_lo - 1, hit.slen - ref_hi
        if hit.strand == "-":
            missing_left, missing_right = missing_right, missing_left
        lo = hit.query_start - missing_left
        hi = hit.query_end + missing_right
        complete = lo >= 1 and hi <= len(contig)
        lo, hi = max(1, lo), min(len(contig), hi)
        piece = contig[lo - 1:hi]
        if hit.strand == "-":
            piece = revcomp(piece)
        note = "" if complete else (
            f"the region runs off the end of contig {hit.qseqid}; "
            f"{hi - lo + 1} of {hit.slen} reference bases recovered")
        found[probe] = Extracted(
            probe=probe, contig=hit.qseqid, start=lo, end=hi, strand=hit.strand,
            nucleotides=piece.upper(), coverage=hit.coverage_pct,
            identity=hit.identity_pct, complete=complete, note=note)
    return found


def _substitutions(observed: str, reference: str, offset: int,
                   neutral: frozenset[str]) -> list[str]:
    """Differences between two equal-length peptides, in absolute numbering.

    The polymorphisms the CDC pipeline lists as neutral are dropped, as it does:
    they are population variation, and reporting them as resistance markers is
    exactly the kind of bare gene name the owner rejected.
    """
    out = []
    for i in range(min(len(observed), len(reference))):
        if observed[i] == reference[i]:
            continue
        local = f"{reference[i]}{i + 1}{observed[i]}"
        if local in neutral:
            continue
        out.append(f"{reference[i]}{i + 1 + offset}{observed[i]}")
    return out


# ---------------------------------------------------------- serotype calling
def _serotype_from_probes(probes: dict[str, ProbeHit]) -> list[tuple[str, str]]:
    """Apply the CDC decision rules, returning (serotype, evidence) candidates.

    Transcribed rule for rule from ``SPN_Serotyper.pl``. The rules are not
    mutually exclusive by construction, so every one that fires is returned and
    the caller decides; disagreement between two of them is itself a finding.
    """
    out: list[tuple[str, str]] = []

    def hit(name: str) -> ProbeHit | None:
        return probes.get(name)

    for gene, serotype in _SINGLE_TARGET.items():
        probe = hit(gene)
        if probe is not None:
            out.append((serotype, f"{gene} {probe.evidence}"))
    for gene, serotype in _SINGLE_TARGET_EXACT.items():
        probe = hit(gene)
        if probe is not None and probe.exact:
            out.append((serotype, f"{gene} {probe.evidence}"))

    # 35B, and 35D when its wciG O-acetyltransferase is not intact.
    wzy35b = hit("WZY35B")
    if wzy35b is not None:
        wcig = hit("WCIG35B")
        if wcig is not None and not wcig.exact:
            out.append(("35B/35D", f"WZY35B {wzy35b.evidence}; "
                                   f"WCIG35B {wcig.evidence}"))
        else:
            out.append(("35B", f"WZY35B {wzy35b.evidence}"
                               + (f"; WCIG35B {wcig.evidence}" if wcig else "")))

    # Serogroup 6: wciP residue 195 (6A/6C vs 6B/6D) x wciN (6A/6B vs 6C/6D).
    wcip_ac, wcip_bd = hit("WCIP6AC"), hit("WCIP6BD")
    wcin_ab, wcin_cd = hit("WCIN6AB"), hit("WCIN6CD")
    if wcip_ac is not None and wcip_ac.exact:
        if wcin_ab is not None:
            out.append(("6A", f"WCIP6AC {wcip_ac.evidence}; "
                              f"WCIN6AB {wcin_ab.evidence}"))
        elif wcin_cd is not None:
            out.append(("6C", f"WCIP6AC {wcip_ac.evidence}; "
                              f"WCIN6CD {wcin_cd.evidence}"))
    if wcip_bd is not None and wcip_bd.exact:
        if wcin_ab is not None:
            out.append(("6B", f"WCIP6BD {wcip_bd.evidence}; "
                              f"WCIN6AB {wcin_ab.evidence}"))
        elif wcin_cd is not None:
            out.append(("6D", f"WCIP6BD {wcip_bd.evidence}; "
                              f"WCIN6CD {wcin_cd.evidence}"))

    wzy7c = hit("WZY7C")
    if wzy7c is not None:
        wchf = hit("WCHF7C")
        serotype = "7C" if wchf is not None else "7B"
        out.append((serotype, f"WZY7C {wzy7c.evidence}"
                              + (f"; WCHF7C {wchf.evidence}" if wchf
                                 else "; WCHF7C absent")))

    wzy7f = hit("WZY7F")
    if wzy7f is not None:
        wcwd = hit("WCWD7F")
        if wcwd is not None and wcwd.exact:
            out.append(("7F", f"WZY7F {wzy7f.evidence}; WCWD7F {wcwd.evidence}"))
        elif wcwd is not None:
            out.append(("7A/7F", f"WZY7F {wzy7f.evidence}; "
                                 f"WCWD7F {wcwd.evidence}"))

    wzy9n = hit("WZY9N")
    if wzy9n is not None:
        wcja_l, wcja_n = hit("WCJA9L"), hit("WCJA9N")
        if wcja_l is not None and wcja_l.exact:
            out.append(("9L", f"WZY9N {wzy9n.evidence}; WCJA9L {wcja_l.evidence}"))
        elif wcja_n is not None and wcja_n.exact:
            out.append(("9N", f"WZY9N {wzy9n.evidence}; WCJA9N {wcja_n.evidence}"))

    wzy9v = hit("WZY9V")
    if wzy9v is not None:
        wcje = hit("WCJE9V")
        if wcje is not None and wcje.exact:
            out.append(("9V", f"WZY9V {wzy9v.evidence}; WCJE9V {wcje.evidence}"))
        elif wcje is not None:
            out.append(("9A/9V", f"WZY9V {wzy9v.evidence}; WCJE9V {wcje.evidence}"))

    wzy15a = hit("WZY15A")
    if wzy15a is not None:
        wciz15f = hit("WCIZ15F")
        if wciz15f is not None and wciz15f.exact:
            out.append(("15F", f"WZY15A {wzy15a.evidence}; "
                               f"WCIZ15F {wciz15f.evidence}"))
        elif wciz15f is None:
            out.append(("15A", f"WZY15A {wzy15a.evidence}; WCIZ15F absent"))

    wzy15b = hit("WZY15B")
    if wzy15b is not None:
        # The CDC rule reads 15B from an intact wciZ and 15C from a disrupted
        # one. That works on reads, where the pileup shows the TA-repeat length
        # and its mixture, but not on an assembly: the discriminator is a short
        # tandem repeat whose length the assembler collapses to one consensus,
        # and it is exactly the feature short-read assemblers resolve worst.
        # So the pair is reported, with the wciZ evidence, and the choice is
        # left to whoever has the reads or the serology.
        wciz = hit("WCIZ15B")
        where = f"WCIZ15B {wciz.evidence}" if wciz else "WCIZ15B absent"
        out.append(("15B/15C", f"WZY15B {wzy15b.evidence}; {where}"))

    wzy18c = hit("WZY18C")
    if wzy18c is not None and wzy18c.exact:
        wcix = hit("WCIX18C")
        if wcix is not None and wcix.exact:
            out.append(("18C", f"WZY18C {wzy18c.evidence}; WCIX18C {wcix.evidence}"))
        elif wcix is not None:
            out.append(("18B/18C", f"WZY18C {wzy18c.evidence}; "
                                   f"WCIX18C {wcix.evidence}"))

    wzy22f = hit("WZY22F")
    if wzy22f is not None:
        wcwa_a, wcwa_f = hit("WCWA22A"), hit("WCWA22F")
        if wcwa_a is not None:
            out.append(("22A", f"WZY22F {wzy22f.evidence}; WCWA22A {wcwa_a.evidence}"))
        elif wcwa_f is not None:
            out.append(("22F", f"WZY22F {wzy22f.evidence}; WCWA22F {wcwa_f.evidence}"))

    wzy33f = hit("WZY33F")
    if wzy33f is not None:
        wcje33a = hit("WCJE33A")
        if wcje33a is not None and wcje33a.exact:
            out.append(("33A", f"WZY33F {wzy33f.evidence}; "
                               f"WCJE33A {wcje33a.evidence}"))
        else:
            out.append(("33F", f"WZY33F {wzy33f.evidence}"))
    return out


def call_serotype(loci: list[LocusMatch], probes: dict[str, ProbeHit]) -> TypingResult:
    """Combine the whole-locus match with the CDC probe rules into one call."""
    alleles: dict[str, str] = {}
    best = loci[0] if loci else None
    runner_up = loci[1] if len(loci) > 1 else None
    if best is not None:
        alleles["cps_locus_best_match"] = best.evidence
    if runner_up is not None:
        alleles["cps_locus_runner_up"] = runner_up.evidence

    # No capsule locus at all: a non-encapsulated (null-capsule) genotype.
    if best is None or best.coverage < CPS_ABSENT_COVERAGE:
        where = best.evidence if best else "no cps alignment at all"
        return TypingResult(
            scheme="cps_serotype", call="non-encapsulated", lineage="null capsule",
            alleles=alleles, loci_found=0, loci_total=1,
            note="No capsular biosynthetic locus is present: the closest "
                 f"reference reaches only {where}, which is the conserved "
                 "wzg-wzh-wzd-wze block and the dexB/aliA flanks that every "
                 "genome carries. This is a non-encapsulated (null-capsule) "
                 "genotype, which is serologically non-typeable and is not "
                 "covered by any conjugate or polysaccharide vaccine. Such "
                 "isolates cause conjunctivitis and otitis media and colonise "
                 "readily, but rarely cause invasive disease. Which "
                 "null-capsule clade it is is not determined: the clades are "
                 "NCC1 (pspK present, aliC and aliD absent), NCC2 (aliC and "
                 "aliD present, pspK absent) and NCC3 (aliD only), and this "
                 "module does not carry those three gene references. Note also "
                 "that a laboratory capsule-deletion mutant such as R6 reads "
                 "the same way as a natural non-encapsulated isolate here.",
        )

    probe_calls = _serotype_from_probes(probes)

    # A locus that is recognisably present but far from complete is a
    # disrupted capsule locus, not a serotype. R6's reads exactly this way: it
    # covers the serotype 2 reference in two blocks either side of a 7.5 kb
    # hole where the serotype-specific genes should be, and no serotype probe
    # fires because wzy has gone with them. Calling it "serotype 2" on the
    # strength of the flanks would be reporting a capsule the genome cannot
    # make. If a serotype probe does fire, the probe is believed and the
    # coverage shortfall is recorded as a warning instead.
    if best.coverage < CPS_MIN_COVERAGE and not probe_calls:
        gaps = _gaps(best.intervals, best.length)
        alleles["cps_locus_gaps"] = (
            ", ".join(f"{lo}-{hi} ({hi - lo + 1} bp)" for lo, hi in gaps[:5])
            or "none contiguous")
        return TypingResult(
            scheme="cps_serotype", call="non-encapsulated",
            lineage=f"disrupted serogroup {serogroup_of(best.serotype)} locus",
            alleles=alleles, loci_found=0, loci_total=1,
            note="The capsule locus is recognisable but incomplete: "
                 f"{best.evidence}, with the region(s) listed above missing, "
                 "and no serotype-specific probe matched. That is a disrupted "
                 "capsule locus, so the isolate is non-encapsulated and "
                 "serologically non-typeable, whatever the flanking sequence "
                 "resembles. A laboratory capsule-deletion derivative such as "
                 "R6 reads this way, and so does a natural isolate that has "
                 "lost its capsule genes. No conjugate or polysaccharide "
                 "vaccine covers a non-encapsulated isolate. The null-capsule "
                 "clade (NCC1, pspK; NCC2, aliC and aliD; NCC3, aliD only) is "
                 "not determined: this module does not carry those references.",
        )
    if best.coverage < CPS_MIN_COVERAGE:
        alleles["cps_locus_coverage_warning"] = (
            f"only {best.coverage:.1f}% of the closest locus is covered, below "
            f"the {CPS_MIN_COVERAGE:.0f}% an intact locus reaches; the serotype "
            "below rests on the probe rules, not on the locus match")
    if best.serotype == "6E":
        alleles["cps_locus_note"] = (
            "the locus matched is the 6E genotype, which is a serogroup 6 cps "
            "variant that is serologically typed as 6B; the serotype reported "
            "is the serological one")

    # Prefer the probe rules: within a serogroup they see what coverage cannot.
    if probe_calls:
        serotypes = sorted({call for call, _ev in probe_calls})
        for i, (call, evidence) in enumerate(sorted(probe_calls), start=1):
            alleles[f"cdc_rule_{i}" if len(probe_calls) > 1 else "cdc_rule"] = \
                f"{call}: {evidence}"
        if len(serotypes) == 1:
            call = serotypes[0]
            conflict = ""
        else:
            # Two rules fired. Where one is a superset of the other (a pair and
            # one of its members) the pair is the honest answer; otherwise the
            # genome is reported as ambiguous between them rather than resolved.
            call = "/".join(sorted({part for s in serotypes
                                    for part in s.split("/")}))
            conflict = (f" More than one CDC rule fired ({', '.join(serotypes)}); "
                        "they are reported together because no evidence here "
                        "chooses between them.")
        members = tuple(call.split("/"))
        included, partial, cross = vaccine_membership(members)
        alleles["serotype_source"] = "CDC serotype probe rules"
        alleles["vaccines_included_in"] = ", ".join(included) or "none"
        if partial:
            alleles["vaccines_partial"] = ", ".join(partial)
        if cross:
            alleles["vaccines_cross_reactive"] = ", ".join(cross)
        note = _serotype_note(call, members, included, partial, cross, best,
                              runner_up, from_probes=True) + conflict
        return TypingResult(
            scheme="cps_serotype", call=call,
            lineage=f"serogroup {serogroup_of(members[0])}",
            alleles=alleles, loci_found=1, loci_total=1, note=note)

    # No rule fired: the serotype is outside the CDC probe set, so the call
    # rests on whole-locus similarity and must say what it cannot separate.
    group = _UNSPLITTABLE_GROUPS.get(serogroup_of(best.serotype))
    if group is not None and best.serotype in group[0]:
        members = group[0]
        alleles["group_reason"] = group[1]
    else:
        members = (best.serotype,)
        siblings = [m.serotype for m in loci[1:4]
                    if abs(best.score - m.score) <= CPS_TIE_MARGIN
                    and serogroup_of(m.serotype) == serogroup_of(best.serotype)]
        if siblings:
            members = tuple([best.serotype] + siblings)
            alleles["group_reason"] = (
                "the runner-up loci are within "
                f"{CPS_TIE_MARGIN:.0f} points of the best match on "
                "coverage x identity, and no discriminating rule covers them")
    call = "/".join(members)
    included, partial, cross = vaccine_membership(members)
    alleles["serotype_source"] = "whole cps locus coverage only"
    alleles["vaccines_included_in"] = ", ".join(included) or "none"
    if partial:
        alleles["vaccines_partial"] = ", ".join(partial)
    if cross:
        alleles["vaccines_cross_reactive"] = ", ".join(cross)
    return TypingResult(
        scheme="cps_serotype", call=call,
        lineage=f"serogroup {serogroup_of(best.serotype)}",
        alleles=alleles, loci_found=1, loci_total=1,
        note=_serotype_note(call, members, included, partial, cross, best,
                            runner_up, from_probes=False))


def _serotype_note(call: str, members: tuple[str, ...], included: list[str],
                   partial: list[str], cross: list[str], best: LocusMatch,
                   runner_up: LocusMatch | None, *, from_probes: bool) -> str:
    parts: list[str] = []
    if len(members) > 1:
        parts.append(f"Serotype {call}, reported as a group because the "
                     "discriminating evidence is missing or does not separate "
                     "its members.")
    else:
        parts.append(f"Serotype {call}, serogroup {serogroup_of(members[0])}.")
    if from_probes:
        parts.append(("Called by the CDC Streptococcus Laboratory's serotype "
                     f"probe rules; the whole-locus match ({best.evidence}) "
                     "corroborates it."))
    else:
        parts.append(("No CDC probe rule covers this serotype, so the call rests "
                     f"on whole-locus similarity alone ({best.evidence})")
                     + (f", with {runner_up.serotype} next at "
                        f"{runner_up.coverage:.1f}% coverage" if runner_up else "")
                     + ". Treat it as a serogroup-level call confirmed to the "
                       "locus, not a serologically verified serotype.")
    if included:
        parts.append("Included in " + ", ".join(included) + ".")
    else:
        parts.append(("Not included in any licensed conjugate or polysaccharide "
                     "vaccine formulation."))
    if partial:
        parts.append("Partly covered: " + "; ".join(partial) + ".")
    if cross:
        parts.append("Expected cross-protection: " + "; ".join(cross) + ".")
    caveat = _SEROTYPE_CAVEATS.get(call)
    if caveat:
        parts.append(f"Caveat: {caveat}.")
    silent = _SILENT_RELATIVES.get(call)
    if silent:
        parts.append(f"Caveat: {silent}.")
    return " ".join(parts)


# -------------------------------------------------------------- PBP + the MIC
@dataclass
class PbpCall:
    """One transpeptidase region's allele, or the reason there is not one."""

    locus: str
    allele: str                     # a number, "NEW", or "-"
    identity: float = 0.0
    coverage: float = 0.0
    nearest: str = ""
    substitutions: list[str] = field(default_factory=list)
    motif_changes: list[str] = field(default_factory=list)
    note: str = ""
    protein: str = ""


def _type_one_pbp(region: PbpRegion, extracted: Extracted | None, db: Path,
                  work: Path) -> PbpCall:
    """Translate one extracted region and match it against the CDC allele set.

    An allele is only named on a 100% identical, full-length protein match, as
    the CDC typer requires; anything else is ``NEW``. Forcing a near match onto
    the nearest number would put a penicillin-susceptible MIC against a mosaic
    PBP, which is the one error this call must not make.
    """
    if extracted is None:
        return PbpCall(locus=region.locus, allele="-",
                       note=f"the {region.gene} transpeptidase region was not "
                            "recovered from the assembly at the required 95% of "
                            "the reference fragment; not determined")
    protein = translate(extracted.nucleotides, table_start_is_met=False)
    protein = protein[:region.length]
    stop = protein.find("*")
    if stop >= 0:
        return PbpCall(locus=region.locus, allele="-", protein=protein,
                       note=f"the translated {region.gene} region carries a stop "
                            f"codon at residue {region.first_residue + stop}; "
                            "the extraction or the assembly is wrong here, so no "
                            "allele is called")
    query = work / f"pbp_{region.locus}_query.faa"
    write_fasta(query, [(f"{region.gene}_query", protein)], wrap=0)
    check_db_exists(db, "prot")
    # Deliberately single-threaded, whatever the caller asked for. BLAST+ 2.16.0
    # has been seen to crash on a protein search when -max_target_seqs is given
    # together with -num_threads > 1 and the search returns nothing, which is
    # exactly the shape of this call on a sample that is not a pneumococcus.
    # Hydra's blast() always passes -max_target_seqs, so the only lever here is
    # the thread count, and it costs nothing: one 300-residue query against at
    # most 1,186 short proteins is instantaneous on one core.
    hsps = blast("blastp", query, db, work / f"pbp_{region.locus}.blastp.tsv",
                 threads=1, evalue=1e-5, max_target_seqs=5000,
                 with_alignment=True)
    if not hsps:
        return PbpCall(locus=region.locus, allele="-", protein=protein,
                       note=f"the translated {region.gene} region matched no "
                            "sequence in the CDC allele set; not determined")
    hsps.sort(key=lambda h: (h.bitscore, h.pident, h.length), reverse=True)
    best = hsps[0]
    allele = best.sseqid.split("||")[0]
    # Substitutions are numbered against allele 0, the sequence of the
    # penicillin-susceptible reference the whole scheme is anchored on, not
    # against whichever allele happens to be nearest: "T338A relative to the
    # susceptible reference" means something, "A12T relative to allele 417"
    # does not.
    zero = next((h for h in hsps if h.sseqid.split("||")[0] == "0"), None)
    exact = (best.pident >= 100.0 and best.length == best.slen
             and best.length == len(protein) and best.gaps == 0)
    if exact:
        return PbpCall(locus=region.locus, allele=allele, identity=best.pident,
                       coverage=100.0, nearest=allele, protein=protein,
                       note=extracted.note)
    coverage = 100.0 * best.length / best.slen if best.slen else 0.0
    if best.pident < PBP_NEW_MIN_IDENTITY or coverage < PBP_NEW_MIN_COVERAGE:
        return PbpCall(locus=region.locus, allele="-", protein=protein,
                       identity=best.pident, coverage=coverage, nearest=allele,
                       note=f"the translated {region.gene} region is only "
                            f"{best.pident:.1f}% identical over {coverage:.0f}% "
                            f"of its closest allele ({allele}); too divergent to "
                            "call, and probably not this gene")
    subs, motifs = _pbp_substitutions(region, zero or best)
    note = extracted.note
    if zero is None:
        note = (note + " " if note else "") + (
            "the substitutions are numbered against the nearest allele rather "
            "than against allele 0, which this region did not align to")
    return PbpCall(locus=region.locus, allele="NEW", identity=best.pident,
                   coverage=coverage, nearest=allele, substitutions=subs,
                   motif_changes=motifs, protein=protein, note=note)


def _pbp_substitutions(region: PbpRegion, hsp: Hsp) -> tuple[list[str], list[str]]:
    """Substitutions a novel allele carries, and which touch the active site.

    Numbering is on the full-length TIGR4 protein, against allele 0 of the CDC
    set, which is the penicillin-susceptible sequence the scheme is anchored on.
    """
    subs: list[str] = []
    notable: list[tuple[int, str]] = []
    motif_spans = [(name, start, start + len(name) - 1)
                   for name, start in region.motifs]
    key = {position: label for position, _wild, label in region.key_sites}
    # Walk the gapped alignment, tracking the position in the reference allele.
    ref_position = hsp.sstart
    for query_aa, ref_aa in zip(hsp.qseq, hsp.sseq):
        if ref_aa == "-":
            continue
        absolute = region.first_residue + ref_position - 1
        if query_aa != ref_aa and query_aa != "-":
            change = f"{ref_aa}{absolute}{query_aa}"
            subs.append(change)
            # One substitution can be both inside an active-site motif and at a
            # named resistance position; say so once, not twice.
            reasons = [f"in the {name} motif ({lo}-{hi})"
                       for name, lo, hi in motif_spans if lo <= absolute <= hi]
            if absolute in key:
                reasons.append(key[absolute])
            if reasons:
                notable.append((absolute, f"{change}, {'; '.join(reasons)}"))
        ref_position += 1
    return subs, [text for _position, text in sorted(set(notable))]


def call_pbp(calls: dict[str, PbpCall]) -> TypingResult:
    """The pbp1a:pbp2b:pbp2x triplet, with the evidence behind each number."""
    ordered = [calls.get(region.locus) or PbpCall(locus=region.locus, allele="-")
               for region in PBP_REGIONS]
    triplet = ":".join(call.allele for call in ordered)
    alleles: dict[str, str] = {}
    notes: list[str] = []
    novel: list[str] = []
    for region, call in zip(PBP_REGIONS, ordered):
        if call.allele not in ("-", "NEW"):
            alleles[region.gene] = call.allele
            alleles[f"{region.gene}_region"] = (
                f"residues {region.first_residue}-"
                f"{region.first_residue + region.length - 1}, exact match to "
                f"CDC allele {call.allele}")
        elif call.allele == "NEW":
            novel.append(region.gene)
            alleles[region.gene] = "NEW"
            alleles[f"{region.gene}_nearest"] = (
                f"closest CDC allele {call.nearest} at {call.identity:.1f}% "
                f"identity over {call.coverage:.0f}% of its length")
            if call.substitutions:
                alleles[f"{region.gene}_substitutions_vs_allele_0"] = (
                    f"{len(call.substitutions)} substitutions against allele 0, "
                    "the penicillin-susceptible reference: "
                    + ", ".join(call.substitutions[:40])
                    + (f" (+{len(call.substitutions) - 40} more)"
                       if len(call.substitutions) > 40 else ""))
            alleles[f"{region.gene}_active_site"] = (
                "; ".join(call.motif_changes) if call.motif_changes
                else "no substitution in the SXXK, SXN or KTG/KSG motifs or at "
                     "the established resistance sites")
        else:
            alleles[region.gene] = "-"
        if call.note:
            notes.append(f"{region.gene}: {call.note}")

    found = sum(1 for call in ordered if call.allele != "-")
    parts = [f"PBP type {triplet}, the pbp1a:pbp2b:pbp2x transpeptidase-region "
             "allele triplet in the CDC Streptococcus Laboratory numbering "
             "(Li et al. 2016)."]
    if novel:
        parts.append(
            f"{', '.join(novel)} carries a transpeptidase region with no exact "
            "match in the reference set and is reported NEW, not rounded to the "
            "nearest allele. The substitutions it carries, and any that fall in "
            "an active-site motif, are listed in the evidence. A novel region "
            "means the MIC cannot be read from the reference table.")
    if found < 3:
        parts.append(f"Only {found} of the 3 regions were typed; the triplet is "
                     "incomplete and no MIC is predicted from it.")
    parts.extend(notes)
    return TypingResult(scheme="pbp_type", call=triplet,
                        lineage="novel allele present" if novel
                        else ("incomplete" if found < 3 else "in reference set"),
                        alleles=alleles, loci_found=found, loci_total=3,
                        note=" ".join(parts))


def _mic_table(path: Path) -> dict[str, dict[str, str]]:
    """The CDC PBP-type to MIC reference table, keyed by its own ``APT`` code."""
    table: dict[str, dict[str, str]] = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (row.get("APT") or "").strip()
            if key:
                table[key] = row
    return table


def _interpret(code: str, mic: float) -> list[str]:
    """CLSI categories for one predicted MIC, in every context that has one."""
    out = []
    for name, (s_max, r_min, label) in _CLSI.items():
        if name.split("/")[0] != code:
            continue
        log2 = round(math.log2(mic)) if mic > 0 else -10
        if log2 <= s_max:
            category = "S"
        elif log2 >= r_min:
            category = "R"
        else:
            category = "I"
        out.append(f"{label}: {category}")
    return out


def call_mic(triplet: str, calls: dict[str, PbpCall], table_path: Path) -> TypingResult:
    """Read the beta-lactam MICs the PBP type predicts, and interpret them.

    Every string this produces says "predicted"; none of it is a measurement.
    A triplet the reference table does not hold gets no MIC at all: Li et al.
    2017 showed newly encountered PBP sequences need the published
    random-forest model to be predicted, and that model is an R object this
    module cannot evaluate.
    """
    complete = all(call.allele not in ("-", "NEW") for call in calls.values()) \
        and len(calls) == 3
    alleles: dict[str, str] = {"pbp_type": triplet, "basis": "predicted from the "
                               "PBP type, not measured"}
    if not complete:
        reason = ("one or more transpeptidase regions is novel or was not typed"
                  if calls else "no transpeptidase region was typed")
        return TypingResult(
            scheme="beta_lactam_mic", call="not determined", lineage="-",
            alleles=alleles, loci_found=0, loci_total=len(PBP_MIC_DRUGS),
            note="No beta-lactam MIC is predicted, because " + reason + ". "
                 "The CDC reference table is keyed on the complete "
                 "pbp1a:pbp2b:pbp2x triplet; Li et al. 2017 validated a "
                 "random-forest predictor for triplets outside that table, but "
                 "it is distributed as an R model and is not reimplemented here. "
                 "Beta-lactam susceptibility must be measured for this isolate.")
    if not table_path.exists():
        raise HydraError(f"PBP MIC table missing at {table_path}; "
                         "run the module's fetch() first")
    table = _mic_table(table_path)
    code = triplet.replace(":", "-")
    row = table.get(code)
    alleles["reference_table_key"] = code
    if row is None:
        return TypingResult(
            scheme="beta_lactam_mic", call="not determined",
            lineage="PBP type not in the reference table", alleles=alleles,
            loci_found=0, loci_total=len(PBP_MIC_DRUGS),
            note=f"PBP type {triplet} is a combination of known alleles that the "
                 f"CDC reference table of {len(table)} PBP types does not hold, "
                 "so no MIC is read from it. Extrapolating from a neighbouring "
                 "triplet is not valid - the three regions interact - and the "
                 "published alternative (Li et al. 2017) is a random-forest "
                 "model distributed as an R object. Beta-lactam susceptibility "
                 "must be measured for this isolate.")
    categories: list[str] = []
    for code_short, drug in PBP_MIC_DRUGS:
        raw = (row.get(code_short) or "").strip()
        if not raw:
            continue
        try:
            mic = float(raw)
        except ValueError:
            continue
        verdicts = _interpret(code_short, mic)
        alleles[f"{drug}_predicted_mic_mg_L"] = raw
        if verdicts:
            alleles[f"{drug}_predicted_category"] = "; ".join(verdicts)
            categories.extend(verdicts)
    pen = alleles.get("penicillin_predicted_mic_mg_L", "")
    headline = f"penicillin MIC {pen} mg/L (predicted)" if pen else "predicted"
    return TypingResult(
        scheme="beta_lactam_mic", call=headline,
        lineage="; ".join(v for v in categories if v.startswith("penicillin"))
        or "-",
        alleles=alleles, loci_found=sum(1 for k in alleles if k.endswith("_mic_mg_L")),
        loci_total=len(PBP_MIC_DRUGS),
        note="These MICs are *predicted* from the PBP transpeptidase type, not "
             "measured, and must not be reported as a susceptibility test "
             f"result. They are read from the CDC reference table entry for PBP "
             f"type {triplet} and interpreted against the CLSI M100 breakpoints "
             "the CDC pipeline applies, including the separate meningitis "
             "criteria for penicillin, cefotaxime and ceftriaxone. "
             + PBP_MIC_ACCURACY +
             " A prediction is not a substitute for testing when the isolate is "
             "from a sterile site or the patient is on a beta-lactam.")


# ----------------------------------------------------------- the rest of them
def call_pilus(probes: dict[str, ProbeHit]) -> TypingResult:
    """Which pilus islets the genome carries."""
    present, alleles = [], {}
    for gene, islet in _PILUS_PROBES.items():
        probe = probes.get(gene)
        if probe is not None:
            present.append(islet)
            alleles[islet] = f"{gene.lower()} present ({probe.evidence})"
        else:
            alleles[islet] = f"{gene.lower()} not detected"
    call = "+".join(sorted(present)) if present else "none"
    note = (("Pilus islet 1 (the rlrA islet: rrgA-rrgB-rrgC with srtB-D) and "
            "islet 2 (pitA-pitB-sipA-srtG1) are adhesive surface structures "
            "carried by a minority of lineages. They are detected here from the "
            "CDC rrgA and pitB marker probes, so this is islet presence, not a "
            "full islet reconstruction and not rrgB clade typing. PI-1 "
            "contributes to adherence and to the inflammatory response and is "
            "enriched in antibiotic-resistant lineages; both islets are protein "
            "vaccine candidates. Neither is used to guide treatment."))
    if not present:
        note += ((" Both marker probes are absent: this genome carries neither "
                 "islet, which is the commonest state."))
    return TypingResult(scheme="pilus_islets", call=call,
                        lineage="-", alleles=alleles,
                        loci_found=len(present), loci_total=len(_PILUS_PROBES),
                        note=note)


def call_macrolide(probes: dict[str, ProbeHit]) -> TypingResult:
    """The macrolide-lincosamide mechanism, and the phenotype it implies."""
    alleles: dict[str, str] = {}
    found: list[tuple[str, str, str]] = []
    for gene, (name, mechanism, phenotype) in _MACROLIDE.items():
        probe = probes.get(gene)
        if probe is None:
            continue
        found.append((name, mechanism, phenotype))
        alleles[name] = f"detected ({probe.evidence})"
    for gene, meaning in _ERMB_CONTEXT.items():
        probe = probes.get(gene)
        if probe is not None:
            alleles[gene.lower()] = f"{meaning} ({probe.evidence})"

    if not found:
        return TypingResult(
            scheme="macrolide_mechanism", call="no acquired determinant detected",
            lineage="-", alleles=alleles, loci_found=0, loci_total=len(_MACROLIDE),
            note="None of erm(B), erm(TR), mef or lsa(C) was detected. That "
                 "predicts macrolide and clindamycin susceptibility, but it is "
                 "not a susceptibility result: macrolide resistance also arises "
                 "from substitutions in 23S rRNA domain V and in ribosomal "
                 "proteins L4 (rplD) and L22 (rplV), which this scheme does not "
                 "screen; the ribosomal-target scheme reports the L4 motif "
                 "separately.")

    names = [name for name, _m, _p in found]
    if any(n.startswith("erm") for n in names) and "mef" in names:
        phenotype = (("MLS-B plus efflux: macrolides and clindamycin both "
                     "resistant, usually at high MIC"))
    elif any(n.startswith("erm") for n in names):
        phenotype = "MLS-B: macrolides and clindamycin resistant"
    elif "mef" in names:
        phenotype = (("M phenotype: 14- and 15-membered macrolides resistant, "
                     "clindamycin active"))
    else:
        phenotype = "LS-A: clindamycin resistant, macrolides unaffected"
    note = ("Mechanism, with what it means for therapy. "
            + " ".join(f"{name} - {mechanism}. {implies}"
                        for name, mechanism, implies in found)
            + " Hydra's AMR gene screen lists these determinants as genes; this "
              "call is the pneumococcal interpretation of them, which is what "
              "decides whether clindamycin can be used.")
    return TypingResult(scheme="macrolide_mechanism", call=", ".join(names),
                        lineage=phenotype, alleles=alleles,
                        loci_found=len(found), loci_total=len(_MACROLIDE),
                        note=note)


def call_other_determinants(probes: dict[str, ProbeHit]) -> TypingResult:
    """Tetracycline and chloramphenicol determinants, interpreted."""
    alleles: dict[str, str] = {}
    found: list[tuple[str, str]] = []
    for gene, (name, meaning) in _OTHER_DETERMINANTS.items():
        probe = probes.get(gene)
        if probe is not None:
            found.append((name, meaning))
            alleles[name] = f"detected ({probe.evidence}): {meaning}"
    if not found:
        return TypingResult(
            scheme="other_determinants", call="none detected", lineage="-",
            alleles=alleles, loci_found=0, loci_total=len(_OTHER_DETERMINANTS),
            note="Neither tet(M), tet(O) nor cat was detected, which predicts "
                 "tetracycline and chloramphenicol susceptibility.")
    return TypingResult(
        scheme="other_determinants",
        call=", ".join(name for name, _meaning in found), lineage="-",
        alleles=alleles, loci_found=len(found),
        loci_total=len(_OTHER_DETERMINANTS),
        note="Acquired determinants outside the beta-lactam and macrolide "
             "classes, with the resistance each predicts. "
             + " ".join(f"{name}: {meaning}" for name, meaning in found))


def call_target_substitutions(extracted: dict[str, Extracted],
                              probes: dict[str, ProbeHit]) -> TypingResult:
    """gyrA, parC, folA and rplD substitutions, in the published numbering."""
    alleles: dict[str, str] = {}
    hits: list[str] = []
    missing: list[str] = []
    for region in TARGET_REGIONS:
        piece = extracted.get(region.probe)
        if piece is None:
            if region.probe not in probes:
                missing.append(f"{region.gene} (region not recovered)")
                alleles[region.gene] = "not determined: region not recovered"
            continue
        protein = translate(piece.nucleotides, table_start_is_met=False)
        protein = protein[:len(region.peptide)]
        if len(protein) < len(region.peptide):
            missing.append(f"{region.gene} (region truncated)")
            alleles[region.gene] = (
                f"not determined: only {len(protein)} of "
                f"{len(region.peptide)} residues recovered")
            continue
        subs = _substitutions(protein, region.peptide, region.offset,
                              region.neutral)
        first = region.offset + 1
        last = region.offset + len(region.peptide)
        # Too many differences is not a divergent allele, it is the wrong
        # region or a frameshifted extraction. Printing twenty substitutions
        # from a mis-framed translation would be worse than printing nothing.
        if len(subs) > TARGET_MAX_DIVERGENT_FRACTION * len(region.peptide):
            missing.append(f"{region.gene} (recovered region implausible)")
            alleles[region.gene] = (
                f"not determined: the recovered region differs from the "
                f"reference at {len(subs)} of {len(region.peptide)} residues, "
                "which means the extraction, not the allele, is wrong")
            continue
        if subs:
            hits.append(f"{region.gene} {','.join(subs)}")
            alleles[region.gene] = (
                f"{', '.join(subs)} (residues {first}-{last}; {region.meaning})")
        else:
            alleles[region.gene] = (
                f"wild type over residues {first}-{last}")
    call = "; ".join(hits) if hits else "wild type at every target screened"
    parts = ["Substitutions in the chromosomal targets that reference "
             "laboratories read from an assembly, numbered on the full protein "
             "exactly as the CDC pipeline numbers them."]
    for region in TARGET_REGIONS:
        value = alleles.get(region.gene, "")
        if value and not value.startswith("wild type") \
                and not value.startswith("not determined"):
            parts.append(f"{region.gene}: {region.meaning}")
    if not hits:
        parts.append(("No substitution was found at any screened target, which "
                     "predicts fluoroquinolone and trimethoprim susceptibility. "
                     "Only the regions listed are screened: a substitution "
                     "outside them, and the folP insertion that carries "
                     "sulfonamide resistance, would not be seen."))
    if missing:
        parts.append("Not determined: " + ", ".join(missing) + ".")
    return TypingResult(
        scheme="target_substitutions", call=call, lineage="-", alleles=alleles,
        loci_found=len(TARGET_REGIONS) - len(missing),
        loci_total=len(TARGET_REGIONS), note=" ".join(parts))


# ------------------------------------------------------------------------ run
def run(assembly: Path, db_dir: Path, *, threads: int = 1,
        workdir: Path | None = None) -> list[TypingResult]:
    """Type one assembly, returning one :class:`TypingResult` per scheme."""
    assembly = Path(assembly)
    root = Path(db_dir) / KEY
    if not root.is_dir():
        raise HydraError(f"{TITLE} typing data not installed at {root}; "
                         "run the module's fetch() first")
    work = Path(workdir) if workdir is not None else assembly.parent
    work.mkdir(parents=True, exist_ok=True)

    results: list[TypingResult] = []

    sero_probes = _probe_hits(assembly, root / "cps" / "sero_probes.fna",
                              work / "spn.sero.blastn.tsv", threads)
    loci = _cps_loci(assembly, root / "cps" / "cps_loci.fna",
                     work / "spn.cps.blastn.tsv", threads)
    results.append(call_serotype(loci, sero_probes))

    regions = _extract_regions(
        assembly, root / "pbp" / "pbp_regions.fna",
        work / "spn.pbp_regions.blastn.tsv", threads,
        wanted={region.probe for region in PBP_REGIONS},
        min_coverage=PBP_MIN_REGION_COVERAGE,
        min_identity=PBP_MIN_REGION_IDENTITY)
    pbp_calls: dict[str, PbpCall] = {}
    for region in PBP_REGIONS:
        pbp_calls[region.locus] = _type_one_pbp(
            region, regions.get(region.probe),
            root / "pbp" / f"pbp_{region.locus}.faa", work)
    pbp_result = call_pbp(pbp_calls)
    results.append(pbp_result)
    results.append(call_mic(pbp_result.call, pbp_calls,
                            root / "pbp" / "pbp_type_mic.csv"))

    results.append(call_pilus(sero_probes))

    res_probes = _probe_hits(assembly, root / "res" / "res_probes.fna",
                             work / "spn.res.blastn.tsv", threads)
    results.append(call_macrolide(res_probes))
    results.append(call_other_determinants(res_probes))
    target_regions = _extract_regions(
        assembly, root / "res" / "res_probes.fna",
        work / "spn.targets.blastn.tsv", threads,
        wanted={region.probe for region in TARGET_REGIONS},
        min_coverage=TARGET_MIN_COVERAGE, min_identity=TARGET_MIN_IDENTITY)
    results.append(call_target_substitutions(target_regions, res_probes))
    return results
