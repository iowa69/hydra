"""Organism-specific typing for *Haemophilus influenzae*.

What is called, in the order a clinical laboratory needs it:

``cap``
    Capsule type a-f from the *cap* locus, separating a genuinely non-typeable
    isolate (no capsule locus at all) from a capsule-deficient encapsulated one,
    and reporting *bexA* status and the duplication/deletion state of the locus
    as far as an assembly can show it. This is the call that decides whether a
    case is vaccine-preventable disease.
``pbp3``
    The *ftsI* (PBP3) substitution pattern, named residue by residue and grouped
    by the published classification, which is what separates a beta-lactamase
    negative ampicillin-resistant isolate from a susceptible one.
``beta_lactam``
    The BLNAS / BLNAR / BLPAR / BLPACR category, from the beta-lactamase gene
    and the PBP3 pattern together.
``mlst``
    Seven-locus sequence type, with the clonal grouping only where PubMLST
    publishes one.

Everything is read from an assembly with BLAST+ and pure Python, and every call
carries the evidence it rests on.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import date
from pathlib import Path

from ...engines.blast import blast, deduplicate, merge_hsps, overlap_length
from ...records import TypingResult
from ...seqio import read_fasta
from ...utils import LOG, HydraError, tempdir
from . import _hinf_nmen_common as common
from .base import OrganismDb, Reference

KEY = "hinfluenzae"
TITLE = "Haemophilus influenzae"
SPECIES: tuple[str, ...] = (
    "Haemophilus influenzae",
    "Haemophilus influenzae biogroup aegyptius",
    "Haemophilus haemolyticus",      # typed so an NTHi/H. haemolyticus mix-up is visible
)
MLST_SCHEMES: tuple[str, ...] = ("hinfluenzae",)

# --------------------------------------------------------------------- data
#: hicap's curated *cap* reference set, pinned to the commit it was read from.
HICAP_COMMIT = "8b7f1385ec5ff167b6498e136ceacd0f04f34af1"
HICAP_RAW = (f"https://raw.githubusercontent.com/scwatts/hicap/{HICAP_COMMIT}"
             "/hicap/database")
HICAP_FILES = ("bexA", "bexB", "bexC", "bexD", "hcsA", "hcsB",
               "type_a", "type_b", "type_c", "type_d", "type_e", "type_f",
               "IS1016V6")

#: PBP3 numbering reference. Every residue number in the published *ftsI*
#: classifications is counted on this protein.
FTSI_REFERENCE = "P45059.1"

#: The two acquired beta-lactamases of *H. influenzae*, taken from the NCBI
#: Reference Gene Catalogue (the nomenclature authority AMRFinderPlus is built
#: on) as its single-gene RefSeq reference records.
BLA_REFERENCES = {
    "blaTEM-1": "NG_050145.1",
    "blaROB-1": "NG_049977.1",
}

PUBMLST_DB = "pubmlst_hinfluenzae_seqdef"
PUBMLST_MLST_SCHEME = 1
MLST_LOCI = ("adk", "atpG", "frdB", "fucK", "mdh", "pgi", "recA")

DATA_VERSION = "2026-11-17"

# ------------------------------------------------------------- cap structure
#: Region II gene sets, serotype by serotype, as hicap's database names them.
#: Region I is *bexDCBA* (capsule export) and region III is *hcsAB*; region II
#: is the serotype-specific biosynthesis cassette between them.
REGION_TWO: dict[str, tuple[str, ...]] = {
    "a": ("acs1", "acs2", "acs3", "acs4"),
    "b": ("bcs1", "bcs2", "bcs3", "bcs4"),
    "c": ("ccs1", "ccs2", "ccs3", "ccs4"),
    "d": ("dcs1", "dcs2", "dcs3", "dcs4", "dcs5"),
    "e": ("ecs1", "ecs2", "ecs3", "ecs4", "ecs5", "ecs6", "ecs7", "ecs8"),
    "f": ("fcs1", "fcs2", "fcs3"),
}
REGION_ONE: tuple[str, ...] = ("bexA", "bexB", "bexC", "bexD")
REGION_THREE: tuple[str, ...] = ("hcsA", "hcsB")

#: What each capsule type means at the bedside.
CAPSULE_MEANING = {
    "a": "serotype a; invasive disease described, mainly in indigenous North American "
         "and Arctic populations; not covered by the Hib conjugate vaccine",
    "b": "serotype b (Hib); the vaccine-preventable serotype, historically the cause of "
         "most paediatric meningitis and epiglottitis; a case in a vaccinated child is "
         "notifiable in most countries and prompts a vaccine-failure investigation",
    "c": "serotype c; rare, usually non-invasive",
    "d": "serotype d; rare, usually non-invasive",
    "e": "serotype e; occasional invasive disease, including in adults",
    "f": "serotype f; the commonest encapsulated type in many post-vaccine settings, "
         "with invasive disease mainly in adults and the immunocompromised",
}

#: hicap's published thresholds, used here so the two tools agree on what counts
#: as a gene (Watts & Holt 2019; hicap 1.0.4 defaults).
GENE_IDENTITY = 70.0
GENE_COVERAGE = 80.0
BROKEN_IDENTITY = 80.0
BROKEN_MIN_BASES = 60
#: Two cap genes further apart than this are not one locus.
LOCUS_GAP = 10000
#: IS1016 filter, from hicap's own sanity filter on its IS search.
IS_MIN_IDENTITY = 60.0
IS_MIN_BASES = 200

# ------------------------------------------------------------------- PBP3
#: Residues of PBP3 whose substitution the literature names, with the reference
#: residue expected at each. Numbering is counted on the Rd KW20 PBP3; the
#: fetch step checks every one of these against the downloaded record, so a
#: changed reference cannot silently shift the numbering.
PBP3_POSITIONS: tuple[tuple[int, str], ...] = (
    (350, "D"), (357, "S"), (368, "A"), (377, "M"), (385, "S"), (389, "L"),
    (392, "P"), (437, "A"), (449, "I"), (490, "G"), (502, "A"), (511, "V"),
    (517, "R"), (526, "N"), (528, "Y"), (530, "A"), (532, "T"), (547, "V"),
    (551, "D"), (554, "A"), (569, "N"),
)
#: The substitutions near the 512-KTG motif that define a resistant PBP3. R517H
#: and N526K are mutually exclusive and neither is subordinate to the other;
#: N526H and Y528H were added to the scheme later.
PBP3_ESSENTIAL = {
    "R517H": "near the KTG motif; found in about 9% of PBP3-mediated resistance",
    "N526K": "near the KTG motif; found in about 84% of PBP3-mediated resistance",
    "N526H": "near the KTG motif; a later addition to the scheme",
    "Y528H": "near the KTG motif; a later addition to the scheme",
}

#: The rPBP3 stage system, as published. Each entry is
#: (stage, group, required substitutions, forbidden substitutions). The first
#: row whose requirements are met and whose exclusions are absent wins, so the
#: rows are ordered from the highest stage down.
PBP3_GROUPS: tuple[tuple[int, str, frozenset[str], frozenset[str]], ...] = (
    (3, "III(+)", frozenset({"S385T", "L389F", "N526K"}), frozenset({"R517H"})),
    (3, "III-like(+)", frozenset({"S385T", "L389F", "R517H"}), frozenset({"N526K"})),
    (2, "III(-)", frozenset({"S385T", "N526K"}), frozenset({"L389F", "R517H"})),
    (2, "III-like(-)", frozenset({"S385T", "R517H"}), frozenset({"L389F", "N526K"})),
    (1, "II", frozenset({"N526K"}), frozenset({"S385T"})),
    (1, "I", frozenset({"R517H"}), frozenset({"S385T"})),
    (1, "Misc1", frozenset({"Y528H"}), frozenset({"S385T", "N526K", "R517H"})),
    (1, "Misc2", frozenset({"N526H"}), frozenset({"S385T", "N526K", "R517H"})),
)

#: What each stage means for treatment.
PBP3_STAGE_MEANING = {
    3: "stage 3, high-level rPBP3: L389F on top of S385T. In a national collection of "
       "cefotaxime-resistant isolates 78.5% were stage 3, and 86% of ceftriaxone-resistant "
       "and 26% of meropenem-resistant isolates (meningitis breakpoints) sat here. "
       "Ampicillin, amoxicillin-clavulanate, cefuroxime and the oral cephalosporins should "
       "be reported resistant, and a third-generation cephalosporin or a carbapenem for "
       "meningitis needs an MIC rather than a susceptible-by-default assumption",
    2: "stage 2, high-level rPBP3: S385T near the SSN motif in addition to the essential "
       "KTG-motif substitution. Cefuroxime and oral cephalosporin activity is reduced and "
       "cefotaxime MICs are raised relative to stage 1",
    1: "stage 1, low-level rPBP3: the essential KTG-motif substitution alone. Ampicillin "
       "MICs sit at or just above the breakpoint, so ampicillin, amoxicillin, "
       "ampicillin-sulbactam and amoxicillin-clavulanate should be reported resistant "
       "whatever the measured MIC; even here 6% of isolates were cefotaxime-non-susceptible "
       "and 20% meropenem-non-susceptible at meningitis breakpoints",
    0: "stage 0: no substitution at any position the scheme uses, so PBP3-mediated "
       "beta-lactam resistance is not predicted (susceptible PBP3, sPBP3)",
}

#: The finer group II subgroups, defined on the residue at 502 and 449. They
#: are not a severity ordering; they are reported because the literature and
#: the reference laboratories use them.
PBP3_SUBGROUPS: tuple[tuple[str, frozenset[str]], ...] = (
    ("IIb", frozenset({"A502V"})),
    ("IIc", frozenset({"A502T"})),
    ("IId", frozenset({"I449V"})),
)

# ------------------------------------------------------------- bibliography
REFERENCES: tuple[Reference, ...] = (
    Reference(
        key="Watts2019", authors="Watts SC, Holt KE", year=2019,
        title="hicap: in silico serotyping of the Haemophilus influenzae capsule locus",
        journal="J Clin Microbiol 57:e00190-19", doi="10.1128/JCM.00190-19",
        pmid="30944197",
        supports="the cap reference gene set, the identity/coverage thresholds, the "
                 "most-complete-region-II serotyping rule and the reference tool the "
                 "capsule call was validated against"),
    Reference(
        key="Kroll1990", authors="Kroll JS, Loynds B, Brophy LN, Moxon ER", year=1990,
        title="The bex locus in encapsulated Haemophilus influenzae: a chromosomal region "
              "involved in capsule polysaccharide export",
        journal="Mol Microbiol 4:1853-1862", doi="10.1111/j.1365-2958.1990.tb02034.x",
        pmid="2082145",
        supports="region I of the cap locus as the four-gene operon bexD, bexC, bexB, bexA"),
    Reference(
        key="Kroll1988", authors="Kroll JS, Hopkins I, Moxon ER", year=1988,
        title="Capsule loss in H. influenzae type b occurs by recombination-mediated "
              "disruption of a gene essential for polysaccharide export",
        journal="Cell 53:347-356", doi="10.1016/0092-8674(88)90155-9", pmid="3284653",
        supports="bexA as the gene whose disruption abolishes capsule export, and the "
                 "duplicated cap locus with bexA in the bridge region"),
    Reference(
        key="Kroll1993", authors="Kroll JS, Moxon ER, Loynds BM", year=1993,
        title="An ancestral mutation enhancing the fitness and increasing the virulence of "
              "Haemophilus influenzae type b",
        journal="J Infect Dis 168:172-176", doi="10.1093/infdis/168.1.172", pmid="8515106",
        supports="the IS1016-bexA partial deletion as the founder mutation of the virulent "
                 "serotype b lineages, so a truncated bexA in a duplicated locus is the "
                 "normal Hib genotype"),
    Reference(
        key="Corn1993", authors="Corn PG, Anders J, Takala AK, Kayhty H, Hoiseth SK",
        year=1993,
        title="Genes involved in Haemophilus influenzae type b capsule expression are "
              "frequently amplified",
        journal="J Infect Dis 167:356-364", doi="10.1093/infdis/167.2.356", pmid="8421169",
        supports="cap locus copy number: two copies is the normal serotype b state, three "
                 "to five copies occur in about a third of invasive isolates and raise "
                 "capsule production"),
    Reference(
        key="Satola2003", authors="Satola SW, Schirmer PL, Farley MM", year=2003,
        title="Complete sequence of the cap locus of Haemophilus influenzae serotype b and "
              "nonencapsulated b capsule-negative variants",
        journal="Infect Immun 71:3639-3644", doi="10.1128/IAI.71.6.3639-3644.2003",
        pmid="12761153",
        supports="region III as hcsA and hcsB, and the reference sequences of the "
                 "single-copy bexA-negative capsule-deficient b-minus locus"),
    Reference(
        key="Satola2007", authors="Satola SW, Collins JT, Napier R, Farley MM", year=2007,
        title="Capsule gene analysis of invasive Haemophilus influenzae: accuracy of "
              "serotyping and prevalence of IS1016 among nontypeable isolates",
        journal="J Clin Microbiol 45:3230-3238", doi="10.1128/JCM.00794-07", pmid="17699642",
        supports="the definitions this module calls on: an isolate with a cap gene but no "
                 "bexA is capsule-deficient, one lacking both is non-typeable, and IS1016 "
                 "is present in 9.5% of non-typeable isolates so it is not a cap marker"),
    Reference(
        key="Kapogiannis2005",
        authors="Kapogiannis BG, Satola S, Keyserling HL, Farley MM", year=2005,
        title="Invasive infections with Haemophilus influenzae serotype a containing an "
              "IS1016-bexA partial deletion: possible association with virulence",
        journal="Clin Infect Dis 41:e97-e103", doi="10.1086/498028", pmid="16267724",
        supports="the IS1016-bexA partial deletion and locus duplication occurring in "
                 "serotype a as well as b, so the rule is not applied to b alone"),
    Reference(
        key="Ubukata2001",
        authors="Ubukata K, Shibasaki Y, Yamamoto K, Chiba N, Hasegawa K, Takeuchi Y, "
                "Sunakawa K, Inoue M, Konno M", year=2001,
        title="Association of amino acid substitutions in penicillin-binding protein 3 "
              "with beta-lactam resistance in beta-lactamase-negative "
              "ampicillin-resistant Haemophilus influenzae",
        journal="Antimicrob Agents Chemother 45:1693-1699",
        doi="10.1128/AAC.45.6.1693-1699.2001", pmid="11353613",
        supports="the original PBP3 groups: R517H defines group I, N526K group II, and "
                 "M377I + S385T + L389F added to N526K defines group III"),
    Reference(
        key="Dabernat2002",
        authors="Dabernat H, Delmas C, Seguy M, Pelissier R, Faucon G, Bennamani S, "
                "Pasquier C", year=2002,
        title="Diversity of beta-lactam resistance-conferring amino acid substitutions in "
              "penicillin-binding protein 3 of Haemophilus influenzae",
        journal="Antimicrob Agents Chemother 46:2208-2218",
        doi="10.1128/AAC.46.7.2208-2218.2002", pmid="12069976",
        supports="the group II subgroups: IIa is N526K alone, IIb adds A502V, IIc A502T "
                 "and IId I449V, and R517H and N526K never occur together"),
    Reference(
        key="GarciaCobos2007",
        authors="Garcia-Cobos S, Campos J, Lazaro E, Roman F, Cercenado E, Garcia-Rey C, "
                "Perez-Vazquez M, Oteo J, de Abajo F", year=2007,
        title="Ampicillin-resistant non-beta-lactamase-producing Haemophilus influenzae in "
              "Spain: recent emergence of clonal isolates with increased resistance to "
              "cefotaxime and cefixime",
        journal="Antimicrob Agents Chemother 51:2564-2573", doi="10.1128/AAC.00354-07",
        pmid="17470649",
        supports="the group III-like category (SSN-motif substitutions with R517H rather "
                 "than N526K) and the mutual exclusivity of R517H and N526K"),
    Reference(
        key="Hasegawa2004",
        authors="Hasegawa K, Chiba N, Kobayashi R, Murayama SY, Iwata S, Sunakawa K, "
                "Ubukata K", year=2004,
        title="Rapidly increasing prevalence of beta-lactamase-nonproducing, "
              "ampicillin-resistant Haemophilus influenzae type b in patients with "
              "meningitis",
        journal="Antimicrob Agents Chemother 48:1509-1514",
        doi="10.1128/AAC.48.5.1509-1514.2004", pmid="15105098",
        supports="the low-BLNAR / BLNAR / BLPACR-I / BLPACR-II categories and the rule "
                 "that S385T on top of the essential substitution separates BLNAR from "
                 "low-BLNAR"),
    Reference(
        key="Skaare2014",
        authors="Skaare D, Anthonisen IL, Caugant DA, Jenkins A, Steinbakk M, Strand L, "
                "Sundsfjord A, Tveten Y, Kristiansen BE", year=2014,
        title="Multilocus sequence typing and ftsI sequencing: a powerful tool for "
              "surveillance of penicillin-binding protein 3-mediated beta-lactam "
              "resistance in nontypeable Haemophilus influenzae",
        journal="BMC Microbiol 14:131", doi="10.1186/1471-2180-14-131", pmid="24884375",
        supports="the rPBP3 / sPBP3 distinction, the high-level and low-level grouping of "
                 "III, III-like, II and I, and the transpeptidase-domain numbering on Rd "
                 "KW20"),
    Reference(
        key="Skaare2025",
        authors="Skaare D, Anthonisen IL, Zecic N, Jenkins A, Caugant DA, Ranheim TE, "
                "Sundsfjord A, Hegstad K", year=2025,
        title="Genomic surveillance of cefotaxime-resistant Haemophilus influenzae in "
              "Norway",
        journal="Front Microbiol 16:1601390", doi="10.3389/fmicb.2025.1601390",
        pmid="40800114",
        supports="the rPBP3 stage system this module implements: stages 0 to 3 on S385T, "
                 "L389F, R517H, N526K, N526H and Y528H, and the cephalosporin and "
                 "carbapenem consequences of each stage"),
    Reference(
        key="Osaki2005",
        authors="Osaki Y, Sanbongi Y, Ishikawa M, Kataoka H, Suzuki T, Maeda K, Ida T",
        year=2005,
        title="Genetic approach to study the relationship between penicillin-binding "
              "protein 3 mutations and Haemophilus influenzae beta-lactam resistance by "
              "using site-directed mutagenesis and gene recombinants",
        journal="Antimicrob Agents Chemother 49:2834-2839",
        doi="10.1128/AAC.49.7.2834-2839.2005", pmid="15980357",
        supports="the causal, not merely associative, role of the individual PBP3 "
                 "substitutions in raising beta-lactam MICs"),
    Reference(
        key="Meats2003",
        authors="Meats E, Feil EJ, Stringer S, Cody AJ, Goldstein R, Kroll JS, Popovic T, "
                "Spratt BG", year=2003,
        title="Characterization of encapsulated and noncapsulated Haemophilus influenzae "
              "and determination of phylogenetic relationships by multilocus sequence "
              "typing",
        journal="J Clin Microbiol 41:1623-1636", doi="10.1128/JCM.41.4.1623-1636.2003",
        pmid="12682154",
        supports="the seven-locus MLST scheme (adk, atpG, frdB, fucK, mdh, pgi, recA) and "
                 "its use for both encapsulated and non-typeable isolates"),
    Reference(
        key="Jolley2018", authors="Jolley KA, Bray JE, Maiden MCJ", year=2018,
        title="Open-access bacterial population genomics: BIGSdb software, the "
              "PubMLST.org website and their applications",
        journal="Wellcome Open Res 3:124", doi="10.12688/wellcomeopenres.14826.1",
        pmid="30345391",
        supports="the PubMLST database the MLST alleles, profiles and clonal-complex "
                 "designations are downloaded from; PubMLST asks that this be cited"),
    Reference(
        key="Feldgarden2021",
        authors="Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, "
                "Haft DH, et al.", year=2021,
        title="AMRFinderPlus and the Reference Gene Catalog facilitate examination of the "
              "genomic links among antimicrobial resistance, stress response, and "
              "virulence",
        journal="Sci Rep 11:12728", doi="10.1038/s41598-021-91456-0", pmid="34135355",
        supports="the NCBI Reference Gene Catalogue the blaTEM-1 and blaROB-1 reference "
                 "records and their nomenclature are taken from"),
)

DATABASES: tuple[OrganismDb, ...] = (
    OrganismDb(
        name="cap",
        source=f"https://github.com/scwatts/hicap (hicap/database/*.fasta) at commit "
               f"{HICAP_COMMIT}",
        version=f"hicap 1.0.4, commit {HICAP_COMMIT}; retrieved {DATA_VERSION}",
        licence="GNU General Public License v3.0. The repository's LICENSE is the verbatim "
                "GPLv3 text and the paper states the tool is 'freely available at "
                "https://github.com/scwatts/hicap under the GNU General Public License v3 "
                "(GPL3)'. The FASTA files are package data of that GPLv3 work.",
        redistribution="permitted",
        reference="Watts2019",
        notes="Redistribution is permitted under GPLv3, which is copyleft: bundling these "
              "files inside an MIT-licensed Windows zip raises a licence-compatibility "
              "question the integrator must settle. Two clean options: (1) download on "
              "first launch, as this module's fetch() does, and ship nothing; or (2) "
              "rebuild an equivalent database from the public-domain GenBank records the "
              "hicap paper's Table 1 cites, which carry no licence at all: CP005967.1 "
              "(KR494, bexA-D and hcsA-B, with the exact coordinates in hicap's own FASTA "
              "headers), Z37516.2 (cap-a, acs1-4), AF549213.1 (cap-b, bcs1-4, the "
              "duplicated locus), HQ651151.1 (cap-c, ccs1-4), HQ424464.1 (cap-d, dcs1-5), "
              "FM882247.1 (cap-e, ecs1-8), AF549211.1 (cap-f, fcs1-3) and X59756.1 "
              "(IS1016-V6). Note that those records are annotated inconsistently - ccs4 is "
              "a misc_feature, fcs3 has no /gene qualifier, and acs2-4 appear as orf2-4 - "
              "so a rebuild needs an alias table and must be re-validated.",
    ),
    OrganismDb(
        name="ftsI",
        source="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&"
               f"id={FTSI_REFERENCE}",
        version=f"{FTSI_REFERENCE}; retrieved {DATA_VERSION}",
        licence="none stated. NCBI states that it 'does not restrict the use or "
                "distribution of the data it contains'; individual submitters may claim "
                "patent or copyright on a sequence, but no licence is attached to this "
                "record.",
        redistribution="permitted",
        reference="Skaare2014",
        notes="One protein record of 610 aa, the PBP3 of H. influenzae Rd KW20. Every "
              "residue number in the rPBP3 scheme is counted on it, and fetch() verifies "
              "all 21 scheme positions against the downloaded sequence rather than "
              "trusting it.",
    ),
    OrganismDb(
        name="bla",
        source="NCBI Reference Gene Catalogue reference records NG_050145.1 (blaTEM-1) "
               "and NG_049977.1 (blaROB-1), fetched from NCBI nuccore",
        version=f"NG_050145.1, NG_049977.1; retrieved {DATA_VERSION}",
        licence="none stated; NCBI RefSeq reference records for antimicrobial resistance.",
        redistribution="permitted",
        reference="Feldgarden2021",
        notes="Two nucleotide records, about 2 kb in total. These are the same reference "
              "records AMRFinderPlus uses, so Hydra's own AMR screen and this module agree "
              "on the gene names.",
    ),
    OrganismDb(
        name="mlst",
        source=f"https://rest.pubmlst.org/db/{PUBMLST_DB} (loci "
               f"{', '.join(MLST_LOCI)} and scheme {PUBMLST_MLST_SCHEME} profiles)",
        version=f"anonymous snapshot retrieved {DATA_VERSION}; PubMLST serves "
                f"unauthenticated clients only data submitted on or before 2024-12-31",
        licence="PubMLST Terms and Conditions (https://pubmlst.org/terms-conditions, "
                "version dated 2026-05-13). The deciding sentence: \"All data submitted on "
                "or before 31 December 2024 may continue to be downloaded, used, and "
                "redistributed without restriction, subject only to proper scientific "
                "citation and acknowledgement.\"",
        redistribution="permitted",
        reference="Jolley2018",
        notes="About 3 MB. Redistribution is permitted for the pre-2025 snapshot an "
              "anonymous download returns, and only for that: the sentence quoted above "
              "places exactly that data outside every restriction, subject to the citation "
              "and the acknowledgement below. Hydra's own pubmlst database holds the same "
              "H. influenzae scheme, and this module uses Hydra's MLST call when it has "
              "one, so a packager that ships Hydra's copy needs only the profile table "
              "from here and not a second copy of the allele sets. "
              "The redistribution right depends on the download being "
              "anonymous: the terms restrict post-2024 data, and the unauthenticated API "
              "serves only pre-2025 data, so the two boundaries coincide. fetch() "
              "therefore sends no credentials and must not be given an API key. PubMLST "
              "requires this acknowledgement wherever the data is used: \"This publication "
              "made use of the PubMLST website (https://pubmlst.org/) sited at the "
              "University of Oxford. The development of that website was funded by the "
              "Wellcome Trust.\"",
    ),
)


# ----------------------------------------------------------------- fetching
def fetch(db_dir: Path) -> None:
    """Download and build this module's reference data into ``db_dir/hinfluenzae``.

    Idempotent: a complete install with the current data version is left alone.
    """
    root = Path(db_dir) / KEY
    stamp = root / "version.json"
    if stamp.exists():
        try:
            if json.loads(stamp.read_text()).get("version") == DATA_VERSION:
                LOG.debug("%s data already installed at %s", KEY, root)
                return
        except (OSError, ValueError):
            pass
    root.mkdir(parents=True, exist_ok=True)
    _fetch_cap(root)
    _fetch_ftsi(root)
    _fetch_bla(root)
    _fetch_mlst(root)
    stamp.write_text(json.dumps({
        "version": DATA_VERSION,
        "retrieved": date.today().isoformat(),
        "hicap_commit": HICAP_COMMIT,
        "ftsI_reference": FTSI_REFERENCE,
    }, indent=2))


def _fetch_cap(root: Path) -> None:
    """Assemble the *cap* reference set from hicap's pinned database."""
    work = root / "cap"
    work.mkdir(parents=True, exist_ok=True)
    records: list[tuple[str, str]] = []
    is_records: list[tuple[str, str]] = []
    for name in HICAP_FILES:
        dest = work / f"{name}.fasta"
        common.download(f"{HICAP_RAW}/{name}.fasta", dest)
        for header, sequence in read_fasta(dest):
            gene = header.split()[0]
            target = is_records if name == "IS1016V6" else records
            target.append((gene, sequence))
        dest.unlink()
    if len(records) < 30:
        raise HydraError(f"{KEY}: the cap reference set downloaded incomplete "
                         f"({len(records)} genes)")
    common.write_subset(records, work / "cap_genes.fna")
    common.makeblastdb(work / "cap_genes.fna", "nucl")
    common.write_subset(is_records, work / "IS1016.fna")
    common.makeblastdb(work / "IS1016.fna", "nucl")

    # A translated copy of the same references. The export genes (bexABCD) and
    # region III (hcsAB) are conserved across serotypes at the protein level but
    # not at the nucleotide level: the reference bexA comes from a serotype f
    # strain and blastn does not align it to a serotype b bexA at all, which
    # would report an intact gene as deleted. The region I genes are stored on
    # the plus strand of the source record rather than in their own orientation,
    # so the reading frame is chosen as the one that translates without an
    # internal stop.
    proteins: list[tuple[str, str]] = []
    for gene, sequence in records:
        protein = common.best_reading_frame(sequence)
        if protein:
            proteins.append((gene, protein))
    if len(proteins) != len(records):
        missing = sorted({g for g, _ in records} - {g for g, _ in proteins})
        raise HydraError(f"{KEY}: no clean reading frame for cap gene(s) "
                         f"{', '.join(missing)}")
    common.write_subset(proteins, work / "cap_proteins.faa")


def _fetch_ftsi(root: Path) -> None:
    """Download the PBP3 numbering reference."""
    work = root / "ftsI"
    work.mkdir(parents=True, exist_ok=True)
    raw = work / "reference.raw.faa"
    common.ncbi_protein_fasta([FTSI_REFERENCE], raw)
    records = [("ftsI", sequence) for _header, sequence in read_fasta(raw)]
    raw.unlink()
    if not records or len(records[0][1]) < 500:
        raise HydraError(f"{KEY}: the PBP3 reference {FTSI_REFERENCE} did not download")
    common.write_subset(records, work / "ftsI.faa")
    for position, residue in PBP3_POSITIONS:
        observed = records[0][1][position - 1:position]
        if observed != residue:
            raise HydraError(
                f"{KEY}: {FTSI_REFERENCE} carries {observed or '?'} at position "
                f"{position}, not the expected {residue}; the residue numbering of the "
                f"published PBP3 groups would not apply to it")


def _fetch_bla(root: Path) -> None:
    """Download the two beta-lactamase references *H. influenzae* carries."""
    work = root / "bla"
    work.mkdir(parents=True, exist_ok=True)
    records: list[tuple[str, str]] = []
    for gene, accession in sorted(BLA_REFERENCES.items()):
        raw = work / f"{gene}.raw.fna"
        common.ncbi_nucleotide_fasta([accession], raw)
        sequences = [sequence for _header, sequence in read_fasta(raw)]
        raw.unlink()
        if not sequences or len(sequences[0]) < 500:
            raise HydraError(f"{KEY}: the beta-lactamase reference {accession} "
                             f"({gene}) did not download")
        records.append((gene, sequences[0]))
    common.write_subset(records, work / "bla.fna")
    common.makeblastdb(work / "bla.fna", "nucl")


def _fetch_mlst(root: Path) -> None:
    """Download the PubMLST *H. influenzae* seven-locus scheme."""
    work = root / "mlst"
    work.mkdir(parents=True, exist_ok=True)
    records: list[tuple[str, str]] = []
    for locus in MLST_LOCI:
        raw = work / f"{locus}.raw.fas"
        common.pubmlst_locus_fasta(PUBMLST_DB, locus, raw)
        for header, sequence in read_fasta(raw):
            name = header.split()[0]
            if not name.startswith(f"{locus}_"):
                name = f"{locus}_{name}"
            records.append((name, sequence))
        raw.unlink()
    if not records:
        raise HydraError(f"{KEY}: no MLST alleles downloaded from PubMLST")
    common.write_subset(records, work / "mlst.fna")
    common.makeblastdb(work / "mlst.fna", "nucl")
    common.pubmlst_profiles(PUBMLST_DB, PUBMLST_MLST_SCHEME, work / "profiles.txt")


# ------------------------------------------------------------------ typing
@dataclass
class CapGene:
    """One cap gene copy located on the assembly."""

    gene: str
    serotype: str               # "" for region I and III
    region: str                 # "I" | "II" | "III"
    contig: str
    start: int
    end: int
    strand: str
    identity: float
    coverage: float
    covered_bases: int
    score: float
    broken: bool
    #: Which search found it: ``blastn`` on the reference CDS, or ``tblastn`` on
    #: its product. The export genes are only conserved at the protein level
    #: across serotypes, so the second search is not a formality.
    method: str = "blastn"


def _region_of(gene: str) -> tuple[str, str]:
    if gene in REGION_ONE:
        return "I", ""
    if gene in REGION_THREE:
        return "III", ""
    for serotype, genes in REGION_TWO.items():
        if gene in genes:
            return "II", serotype
    return "?", ""


def run(assembly: Path, db_dir: Path, *, threads: int = 1,
        workdir: Path | None = None, mlst=None) -> list[TypingResult]:
    """Type one *H. influenzae* assembly, returning one result per scheme."""
    assembly = Path(assembly)
    root = Path(db_dir) / KEY
    if workdir is not None:
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        return _run_in(assembly, root, workdir, threads, mlst)
    with tempdir(prefix="hydra.hinfluenzae.") as work:
        return _run_in(assembly, root, work, threads, mlst)


def _run_in(assembly: Path, root: Path, work: Path, threads: int,
            mlst=None) -> list[TypingResult]:
    results: list[TypingResult] = []
    contigs = common.load_contigs(assembly)
    local = work / "assembly.fna"
    common.write_subset(contigs.items(), local)
    common.makeblastdb(local, "nucl")

    results.append(_call_capsule(local, contigs, root, work, threads))
    pbp3 = _call_pbp3(local, contigs, root, work, threads)
    results.append(pbp3)
    results.append(_call_beta_lactam(local, root, work, threads, pbp3))
    results.append(_call_mlst(local, root, work, threads, mlst))
    return results


# ------------------------------------------------------------------ capsule
def _call_capsule(assembly: Path, contigs: dict[str, str], root: Path,
                  work: Path, threads: int) -> TypingResult:
    """Call the capsule type and the state of the *cap* locus."""
    fasta = root / "cap" / "cap_genes.fna"
    if not fasta.exists():
        return TypingResult(scheme="cap", call="not determined",
                            note=f"cap reference data missing; run 'hydra db download {KEY}'")
    genes = _find_cap_genes(assembly, root, work, threads)
    if not genes:
        is_only = _find_is1016(assembly, root, work, threads)
        return TypingResult(
            scheme="cap", call="NT", lineage="non-typeable (no capsule locus)",
            loci_found=0, loci_total=len(REGION_ONE) + len(REGION_THREE),
            alleles={
                "cap_locus": f"no bex, hcs or region II gene reached "
                             f"{GENE_IDENTITY:.0f}% identity over {GENE_COVERAGE:.0f}% of "
                             f"its reference, by either blastn or tblastn",
                "IS1016_copies": f"{len(is_only)}; IS1016 is not a capsule marker and its "
                                 f"presence does not make an isolate encapsulated",
            },
            note="no part of the cap locus was detected. By the definition of Satola et al. "
                 "2007, an isolate 'lacking both bexA and any of the cap genes' is "
                 "non-typeable, so this is a genuinely acapsular isolate rather than a "
                 "capsule-deficient one. Non-typeable H. influenzae causes otitis media, "
                 "sinusitis, exacerbations of COPD and invasive disease in neonates, the "
                 "elderly and the immunocompromised, and is not vaccine-preventable. "
                 "Confirm the species, because H. haemolyticus is also acapsular and is "
                 "commonly misidentified as non-typeable H. influenzae")
    clusters = _cluster(genes)
    best = max(clusters, key=lambda c: (_serotype_score(c)[1], len(c)))
    serotype, score, per_type = _serotype_score(best)
    region_two = [g for g in best if g.region == "II"]
    region_one = {g.gene: g for g in genes if g.region == "I"}
    region_three = {g.gene: g for g in genes if g.region == "III"}

    alleles: dict[str, str] = {}
    for gene in REGION_ONE:
        hit = region_one.get(gene)
        alleles[gene] = _gene_evidence(hit)
    for gene in REGION_THREE:
        hit = region_three.get(gene)
        alleles[gene] = _gene_evidence(hit)
    expected = REGION_TWO.get(serotype, ())
    found_two = {g.gene: g for g in region_two}
    for gene in expected:
        alleles[f"regionII_{gene}"] = _gene_evidence(found_two.get(gene))

    is_hits = _find_is1016(assembly, root, work, threads)
    alleles["IS1016_copies"] = (
        f"{len(is_hits)} (IS1016 is not a capsule marker: it flanks the division I cap "
        f"loci but is also carried by about 9.5% of genuinely non-typeable isolates, and "
        f"a full-length copy remains in the cap-deleted laboratory strain Rd KW20)")
    duplicated, copies, copy_note = _copy_number(genes)
    alleles["cap_locus_copies"] = copies
    bexa = region_one.get("bexA")
    alleles["bexA_status"] = _bexa_status(bexa)

    missing_two = [g for g in expected if g not in found_two]
    broken = sorted({g.gene for g in best if g.broken} |
                    {g.gene for g in genes if g.region in ("I", "III") and g.broken})
    contigs_used = sorted({g.contig for g in best} | {g.contig for g in region_one.values()}
                          | {g.contig for g in region_three.values()})
    alleles["locus_contigs"] = ", ".join(contigs_used) or "-"
    alleles["serotype_support"] = "; ".join(
        f"{name}: {len(hits)}/{len(REGION_TWO[name])} region II genes"
        for name, hits in sorted(per_type.items(), key=lambda kv: -len(kv[1]))
        if hits) or "none"

    notes: list[str] = []
    call = serotype or "not determined"
    bexa_intact = bexa is not None and not bexa.broken and bexa.coverage >= 95.0

    if not region_two:
        # Region I and/or III without any region II gene: the biosynthesis
        # cassette has gone. That is a capsule-deficient isolate, not an
        # acapsular one, and the two are different organisms clinically.
        call = "capsule-deficient"
        lineage = "cap locus present, region II absent, serotype not nameable"
        notes.append(
            "region I and/or region III genes are present but no serotype-specific region "
            "II gene was found, so the isolate carries a remnant capsule locus and cannot "
            "make a capsule. By the definition of Satola et al. 2007 an isolate 'containing "
            "a cap gene but not bexA' is a capsule-deficient variant, and only one 'lacking "
            "both bexA and any of the cap genes' is non-typeable")
    elif bexa_intact:
        lineage = f"encapsulated, serotype {serotype}"
        notes.append(f"a complete region II cassette for serotype {serotype} with an intact "
                     f"bexA: this isolate is expected to express capsule")
    elif len(contigs_used) > 1 and not duplicated:
        # Without the copy number there is no way to tell the normal Hib
        # arrangement from the capsule-minus state, and guessing either way
        # would be a clinical error in opposite directions.
        lineage = (f"serotype {serotype}; capsule expression not determined "
                   f"(bexA truncated, copy number unresolved)")
        notes.append(
            "bexA is absent or truncated and the locus is split across contigs, so the "
            "copy number cannot be read. A truncated bexA in one copy of a duplicated "
            "locus is the normal arrangement of virulent serotype a and b lineages, while "
            "the same finding in a single-copy locus is the capsule-minus state; these "
            "cannot be separated from this assembly. Resolve with long reads, or with the "
            "PCR assays for the IS1016-bexA junction, before reporting capsule status")
    elif duplicated:
        # The textbook Hib locus is a tandem duplication in which one copy has
        # lost bexA and the IS1016 junction with it. Calling that
        # capsule-deficient would mislabel the normal, fully virulent genotype.
        lineage = f"encapsulated, serotype {serotype} (IS1016-bexA deletion variant)"
        notes.append(
            "bexA is absent or truncated in one copy of a duplicated locus, with at least "
            "one copy of the rest of the locus intact. This is the IS1016-bexA partial "
            "deletion: a roughly 1.2 kb deletion removing the 3' end of IS1016 and about "
            "83% of one bexA copy, described as the founder mutation of the virulent "
            "serotype b lineages and also found in virulent serotype a. It is the normal, "
            "capsule-expressing arrangement and not a defect; a complete cap-b duplication "
            "without bexA truncation has not been observed")
    else:
        # The same shape as the meningococcal module: the serotype first and the
        # lesion stated with it, so the row cannot be read as an encapsulated
        # isolate and the genotype is not hidden in another column.
        call = common.qualified_call(serotype, common.DISRUPTED)
        lineage = f"capsule-deficient, genotypically serotype {serotype}"
        notes.append(
            f"a single copy of the locus with bexA absent or truncated: by the Satola et al. "
            f"2007 definition this is a capsule-deficient ({serotype}-minus) variant. It "
            f"carries the serotype {serotype} genotype but cannot export capsular "
            f"polysaccharide, so it is serologically non-typeable while being, by descent, "
            f"an encapsulated lineage. A slide-agglutination report of 'non-typeable' does "
            f"not contradict this call")
    if serotype:
        lineage = lineage + "; " + CAPSULE_MEANING.get(serotype, "")
    if missing_two or broken:
        notes.append("incomplete elsewhere in the locus: " + "; ".join(filter(None, [
            f"region II genes not found: {', '.join(missing_two)}" if missing_two else "",
            f"region I only {len(region_one)}/{len(REGION_ONE)} genes"
            if len(region_one) < len(REGION_ONE) else "",
            f"truncated or fragmented: {', '.join(broken)}" if broken else "",
        ])))
    if len(contigs_used) > 1:
        notes.append(f"the locus is split across {len(contigs_used)} contigs, so gene order, "
                     f"copy number and the extent of any deletion cannot be read from this "
                     f"assembly")
    if copy_note:
        notes.append(copy_note)
    if len(per_type) > 1 and sum(1 for hits in per_type.values() if hits) > 1:
        rivals = ", ".join(f"{name} ({len(hits)})" for name, hits in sorted(per_type.items())
                           if hits and name != serotype)
        if rivals:
            notes.append(f"region II genes of other serotypes also matched ({rivals}); the "
                         f"serotype is the one with the most complete region II set, not the "
                         f"best single hit, because the a/b, c/f and d/e gene pairs are "
                         f"84-94% identical and all clear the 70% identity threshold")

    found = len(region_one) + len(region_three) + len(found_two)
    total = len(REGION_ONE) + len(REGION_THREE) + (len(expected) or 0)
    return TypingResult(scheme="cap", call=call, lineage=lineage or "-",
                        alleles=alleles, loci_found=found, loci_total=total,
                        note=". ".join(notes))


def _gene_evidence(hit: CapGene | None) -> str:
    if hit is None:
        return "not detected"
    state = "truncated/fragment" if hit.broken else "complete"
    unit = "amino-acid" if hit.method == "tblastn" else "nucleotide"
    return (f"{state}, {hit.identity:.1f}% {unit} identity over {hit.coverage:.0f}% of "
            f"the reference, {hit.contig}:{hit.start}-{hit.end}({hit.strand}), "
            f"by {hit.method}")


def _bexa_status(hit: CapGene | None) -> str:
    """Report bexA, whose loss is what separates capsule-positive from minus.

    bexA encodes the ATP-binding subunit of the capsule export machinery, so
    losing it abolishes capsule expression whatever the rest of the locus says.
    It is a poor serotype discriminator (0.17-0.85% divergence between serotypes
    a to d) and is not used here to name the serotype.
    """
    if hit is None:
        return ("not detected: bexA-negative. Read together with the copy number: in one "
                "copy of a duplicated locus this is the normal IS1016-bexA deletion of "
                "virulent serotype a and b lineages, while in a single-copy locus it is "
                "the capsule-minus state")
    if hit.broken or hit.coverage < 95.0:
        return (f"partial: {hit.coverage:.0f}% of the reference at {hit.identity:.1f}% "
                f"identity, by {hit.method}; consistent with the IS1016-bexA partial "
                f"deletion, which removes about 83% of one bexA copy")
    return (f"intact: {hit.coverage:.0f}% coverage at {hit.identity:.1f}% identity, "
            f"by {hit.method}")


def _find_cap_genes(assembly: Path, root: Path, work: Path,
                    threads: int) -> list[CapGene]:
    """Locate every cap gene copy on the assembly, by two searches.

    A nucleotide search is done first, with hicap's thresholds, so the two tools
    agree on what counts as a gene. It is then backed by a translated search,
    because the reference export genes come from a single serotype f strain and
    blastn will not align that *bexA* to a serotype b *bexA* at all: strain
    10810, whose annotated *bexA* is intact, comes out bexA-negative on the
    nucleotide search alone. A gene the nucleotide search missed is taken from
    the translated search and labelled as such.

    Each physical locus keeps only its best-matching reference. Serotypes a/b,
    c/f and d/e share region II genes at 85-95% nucleotide identity, so without
    that step one gene is reported once per serotype and the serotype count
    becomes meaningless.
    """
    out = work / "cap.blastn.tsv"
    hsps = blast("blastn", assembly, root / "cap" / "cap_genes.fna", out,
                 threads=threads, evalue=1e-10, task="blastn",
                 perc_identity=60.0, max_target_seqs=10000)
    candidates: list[CapGene] = []
    for hit in merge_hsps(hsps):
        complete = hit.identity_pct >= GENE_IDENTITY and hit.coverage_pct >= GENE_COVERAGE
        fragment = hit.identity_pct >= BROKEN_IDENTITY and hit.covered_bases >= BROKEN_MIN_BASES
        if not (complete or fragment):
            continue
        region, serotype = _region_of(hit.sseqid)
        if region == "?":
            continue
        candidates.append(CapGene(
            gene=hit.sseqid, serotype=serotype, region=region, contig=hit.qseqid,
            start=hit.query_start, end=hit.query_end, strand=hit.strand,
            identity=hit.identity_pct, coverage=hit.coverage_pct,
            covered_bases=hit.covered_bases,
            score=hit.coverage_pct * hit.identity_pct / 100.0,
            broken=not complete, method="blastn",
        ))
    # One name per physical locus: the best-matching reference wins.
    nucleotide = deduplicate(
        candidates,
        key_span=lambda g: (g.start, g.end),
        key_seq=lambda g: g.contig,
        key_score=lambda g: g.score,
        key_tiebreak=lambda g: g.gene,
        overlap_fraction=0.5,
    )
    protein = _find_cap_proteins(assembly, root, work, threads)
    claimed: dict[str, list[tuple[int, int]]] = {}
    for gene in nucleotide:
        claimed.setdefault(gene.contig, []).append((gene.start, gene.end))
    extra: list[CapGene] = []
    for gene in protein:
        length = max(1, gene.end - gene.start + 1)
        if any(overlap_length((gene.start, gene.end), prev) > 0.5 * length
               for prev in claimed.get(gene.contig, [])):
            continue
        extra.append(gene)
        claimed.setdefault(gene.contig, []).append((gene.start, gene.end))
    return nucleotide + extra


def _find_cap_proteins(assembly: Path, root: Path, work: Path,
                       threads: int) -> list[CapGene]:
    """Translated search for the cap genes, as a backstop for divergent copies."""
    proteins = root / "cap" / "cap_proteins.faa"
    if not proteins.exists():
        return []
    lengths = {header.split()[0]: len(sequence)
               for header, sequence in read_fasta(proteins)}
    hits = common.screen_proteins(assembly, proteins, work, threads=threads,
                                  min_identity=GENE_IDENTITY, min_coverage=0.0,
                                  tag="cap_prot")
    candidates: list[CapGene] = []
    for hit in hits:
        region, serotype = _region_of(hit.gene)
        if region == "?":
            continue
        residues = hit.coverage_pct / 100.0 * lengths.get(hit.gene, 0)
        complete = hit.coverage_pct >= GENE_COVERAGE
        fragment = (hit.identity_pct >= BROKEN_IDENTITY
                    and residues >= BROKEN_MIN_BASES / 3.0)
        if not (complete or fragment):
            continue
        candidates.append(CapGene(
            gene=hit.gene, serotype=serotype, region=region, contig=hit.contig,
            start=hit.start, end=hit.end, strand=hit.strand,
            identity=hit.identity_pct, coverage=hit.coverage_pct,
            covered_bases=int(residues * 3),
            score=hit.coverage_pct * hit.identity_pct / 100.0,
            broken=not complete, method="tblastn",
        ))
    return deduplicate(
        candidates,
        key_span=lambda g: (g.start, g.end),
        key_seq=lambda g: g.contig,
        key_score=lambda g: g.score,
        key_tiebreak=lambda g: g.gene,
        overlap_fraction=0.5,
    )


def _cluster(genes: list[CapGene]) -> list[list[CapGene]]:
    """Group cap genes into candidate loci by position on a contig."""
    by_contig: dict[str, list[CapGene]] = {}
    for gene in genes:
        by_contig.setdefault(gene.contig, []).append(gene)
    clusters: list[list[CapGene]] = []
    for contig_genes in by_contig.values():
        contig_genes.sort(key=lambda g: g.start)
        current = [contig_genes[0]]
        for gene in contig_genes[1:]:
            if gene.start - max(g.end for g in current) <= LOCUS_GAP:
                current.append(gene)
            else:
                clusters.append(current)
                current = [gene]
        clusters.append(current)
    if not clusters:
        return [[]]
    # A fragmented assembly leaves pieces of one locus on several contigs; they
    # are still one locus for the purpose of naming the serotype.
    return clusters + [genes] if len(clusters) > 1 else clusters


def _serotype_score(cluster: list[CapGene]) -> tuple[str, float, dict[str, list[CapGene]]]:
    """Name the serotype of one candidate locus and say how well it is supported."""
    per_type: dict[str, list[CapGene]] = {name: [] for name in REGION_TWO}
    for gene in cluster:
        if gene.region == "II" and gene.serotype:
            per_type[gene.serotype].append(gene)
    scored = sorted(
        ((len(hits), sum(h.score for h in hits), name) for name, hits in per_type.items()),
        reverse=True)
    count, score, name = scored[0]
    return (name if count else "", score, per_type)


def _find_is1016(assembly: Path, root: Path, work: Path, threads: int) -> list:
    """Count IS1016 copies, using hicap's own filter on its IS search."""
    fasta = root / "cap" / "IS1016.fna"
    if not fasta.exists():
        return []
    hits = common.screen_nucleotide(assembly, fasta, work, threads=threads,
                                    min_identity=IS_MIN_IDENTITY, min_coverage=0.0,
                                    tag="is1016")
    return [h for h in hits
            if h.identity_pct >= IS_MIN_IDENTITY
            and (h.coverage_pct / 100.0) * 712 >= IS_MIN_BASES]


def _copy_number(genes: list[CapGene]) -> tuple[bool, str, str]:
    """How many copies of the locus the assembly shows, and what that is worth.

    Returns (duplicated, description, caveat). A tandem duplication of
    near-identical sequence collapses in a short-read assembly, so two copies
    seen in the assembly are evidence of a duplication while one copy is not
    evidence against one. Serotype b normally carries two copies and 23/66
    invasive isolates in the published survey carried three to five, each extra
    copy raising capsule production.
    """
    counts: dict[str, int] = {}
    for gene in genes:
        if gene.gene.startswith("bexA"):
            # bexA is the gene deleted from one copy, so counting it would
            # report the duplication as a single copy.
            continue
        counts[gene.gene] = counts.get(gene.gene, 0) + 1
    duplicated = sorted(name for name, count in counts.items() if count > 1)
    if len(duplicated) >= 2:
        copies = max(counts.values())
        return (True,
                f"at least {copies} copies of the locus "
                f"({len(duplicated)} genes found more than once)",
                f"{', '.join(duplicated)} are present in more than one copy, which is "
                f"direct evidence of a duplicated or amplified locus; the count is a "
                f"lower bound because an assembler collapses identical repeats")
    if duplicated:
        # One gene twice is as likely to be a misassembled repeat end as a real
        # duplication, so it is reported without being called one.
        return (False, f"one copy resolved, with {duplicated[0]} found twice",
                f"{duplicated[0]} was found in two copies but no other cap gene was, which "
                f"is too little to call a duplicated locus")
    return (False, "one copy resolved",
            "only one copy of each cap gene is assembled. A short-read assembly collapses "
            "tandem repeats, so this neither shows nor excludes the duplication typical of "
            "serotypes a and b; confirm by PCR or long reads if the copy number matters")


# --------------------------------------------------------------------- PBP3
def _call_pbp3(assembly: Path, contigs: dict[str, str], root: Path,
               work: Path, threads: int) -> TypingResult:
    """Call the PBP3 substitution pattern from *ftsI*."""
    proteins = root / "ftsI" / "ftsI.faa"
    if not proteins.exists():
        return TypingResult(scheme="pbp3", call="not determined",
                            note=f"PBP3 reference missing; run 'hydra db download {KEY}'")
    reference = next((seq for _h, seq in read_fasta(proteins)), "")
    hits = common.screen_proteins(assembly, proteins, work, threads=threads,
                                  min_identity=80.0, min_coverage=70.0, tag="ftsI")
    if not hits:
        return TypingResult(
            scheme="pbp3", call="not determined", loci_total=1,
            note="ftsI was not found at >=80% identity over >=70% of PBP3; the gene is "
                 "essential, so this means the assembly is incomplete at this locus or "
                 "the species call is wrong, not that the gene is absent")
    hit = max(hits, key=lambda h: h.coverage_pct * h.identity_pct)
    cds = common.extract_cds(contigs, hit.contig, hit.start, hit.end, hit.strand,
                             ref_offset=max(0, hit.query_start - 1))
    protein = cds.protein if cds else ""
    alleles: dict[str, str] = {
        "ftsI_locus": hit.evidence(),
        "ftsI_reference": f"{FTSI_REFERENCE} (PBP3 of H. influenzae Rd KW20, "
                          f"{len(reference)} aa; all residue numbers are counted on it)",
    }
    caveats: list[str] = []
    if cds is None or not protein:
        return TypingResult(scheme="pbp3", call="not determined", alleles=alleles,
                            loci_found=1, loci_total=1,
                            note="ftsI was located but no open reading frame could be "
                                 "recovered from the assembly, so the substitutions "
                                 "cannot be read")
    if cds.internal_stops:
        caveats.append(f"{cds.internal_stops} internal stop codon(s) in the assembled ftsI; "
                       f"PBP3 is essential, so this is almost certainly an assembly or "
                       f"sequencing error rather than a real truncation")
    if cds.truncated_at_contig_end:
        caveats.append("the ftsI reading frame runs into a contig end, so residues beyond "
                       "the break were not examined")
    if len(protein) != len(reference):
        caveats.append(f"the assembled PBP3 is {len(protein)} aa against the reference's "
                       f"{len(reference)} aa, so residues were numbered through the "
                       f"gapped alignment rather than by position")

    subs = common.substitutions(reference, protein,
                                aligned_query=hit.aligned_query,
                                aligned_subject=hit.aligned_subject,
                                offset=max(0, hit.query_start - 1))
    named = _named_substitutions(subs)
    alleles["ftsI_substitutions"] = ", ".join(subs) if subs else "none"
    alleles["pbp3_scheme_positions"] = ", ".join(named) if named else "none"
    found = set(named)
    essential = sorted(found & set(PBP3_ESSENTIAL))
    alleles["essential_substitution"] = (
        "; ".join(f"{sub}: {PBP3_ESSENTIAL[sub]}" for sub in essential)
        if essential else
        "none; without a substitution near the 512-KTG motif the pattern is not a "
        "resistant PBP3 whatever else is substituted")
    stage, group = _pbp3_group(found)
    subgroup = _pbp3_subgroup(group, found)
    if group:
        call = f"stage {stage} group {group}" + (f" ({subgroup})" if subgroup else "")
        lineage = "resistant PBP3 (rPBP3)"
        defining = sorted(next(req for st, gp, req, _f in PBP3_GROUPS
                               if gp == group and st == stage) & found)
        parts = [PBP3_STAGE_MEANING[stage],
                 f"group {group} is defined by {', '.join(defining)}, at the residues "
                 f"flanking the 379-SSN and 512-KTG motifs of the transpeptidase domain"]
        if subgroup:
            parts.append(f"within group II the residue at 502/449 places it in subgroup "
                         f"{subgroup}, which is a descriptive subdivision and not a "
                         f"severity grade")
        parts.append(f"all substitutions found at scheme positions: {', '.join(named)}")
        note = ". ".join(parts)
    elif named:
        call = "stage 0, substituted outside the scheme"
        lineage = "susceptible PBP3 (sPBP3) with incidental substitutions"
        note = (f"substitutions were found at positions the scheme lists "
                f"({', '.join(named)}) but none of them is an essential KTG-motif "
                f"substitution, so this is a susceptible PBP3. Isolates of this kind are "
                f"reported as beta-lactamase-negative ampicillin-susceptible unless the "
                f"MIC says otherwise. {PBP3_STAGE_MEANING[0]}")
    else:
        call = "stage 0, wild type"
        lineage = "susceptible PBP3 (sPBP3)"
        note = PBP3_STAGE_MEANING[0] + ("; this is the Rd KW20 reference sequence over "
                                        "every position the scheme uses")
    if caveats:
        note = note + ". " + "; ".join(caveats)
    return TypingResult(scheme="pbp3", call=call, lineage=lineage, alleles=alleles,
                        loci_found=1, loci_total=1, note=note)


def _named_substitutions(subs: list[str]) -> list[str]:
    """Keep the substitutions at positions the published groups are built on."""
    wanted = {position: residue for position, residue in PBP3_POSITIONS}
    out: list[str] = []
    for sub in subs:
        digits = "".join(ch for ch in sub if ch.isdigit())
        if not digits:
            continue
        position = int(digits)
        if position in wanted and sub.startswith(wanted[position]):
            out.append(sub)
    return sorted(out, key=lambda s: int("".join(ch for ch in s if ch.isdigit())))


def _pbp3_group(found: set[str]) -> tuple[int, str]:
    """Match a substitution set against the rPBP3 stage definitions."""
    for stage, group, required, forbidden in PBP3_GROUPS:
        if required <= found and not (forbidden & found):
            return stage, group
    return 0, ""


def _pbp3_subgroup(group: str, found: set[str]) -> str:
    """The group II subdivision, where one applies."""
    if group != "II":
        return ""
    for label, required in PBP3_SUBGROUPS:
        if required <= found:
            return label
    return "IIa"


# ------------------------------------------------------------- beta-lactams
def _call_beta_lactam(assembly: Path, root: Path, work: Path, threads: int,
                      pbp3: TypingResult) -> TypingResult:
    """Combine the beta-lactamase gene and the PBP3 pattern into one category."""
    fasta = root / "bla" / "bla.fna"
    alleles: dict[str, str] = {}
    if not fasta.exists():
        bla_found: list = []
        alleles["beta_lactamase"] = "not determined: reference data missing"
    else:
        bla_found = common.screen_nucleotide(assembly, fasta, work, threads=threads,
                                             min_identity=90.0, min_coverage=80.0,
                                             tag="bla")
        if bla_found:
            alleles["beta_lactamase"] = "; ".join(
                f"{h.gene} ({h.identity_pct:.1f}% identity, {h.coverage_pct:.0f}% coverage, "
                f"{h.contig}:{h.start}-{h.end})" for h in sorted(bla_found, key=lambda h: h.gene))
        else:
            alleles["beta_lactamase"] = "no blaTEM-1 or blaROB-1 found at >=90% identity"
    alleles["pbp3"] = f"{pbp3.call} ({pbp3.lineage})"

    resistant_pbp3 = pbp3.lineage == "resistant PBP3 (rPBP3)"
    undetermined = pbp3.call == "not determined" or "not determined" in alleles["beta_lactamase"]
    if undetermined:
        return TypingResult(
            scheme="beta_lactam", call="not determined", alleles=alleles,
            note="the category needs both the beta-lactamase status and the PBP3 pattern, "
                 "and one of them could not be read from this assembly")
    if bla_found and resistant_pbp3:
        roman = "II" if pbp3.call.startswith(("stage 2", "stage 3")) else "I"
        call, meaning = f"BLPACR (BLPACR-{roman})", (
            "beta-lactamase positive, amoxicillin-clavulanate resistant: the enzyme is "
            "inhibited by clavulanate but the altered PBP3 is not, so neither ampicillin "
            "nor amoxicillin-clavulanate is reliable; a third-generation cephalosporin or "
            "a carbapenem is the usual choice for invasive disease")
    elif bla_found:
        call, meaning = "BLPAR", (
            "beta-lactamase positive, ampicillin resistant: ampicillin and amoxicillin "
            "fail, amoxicillin-clavulanate is expected to work, and this is the commonest "
            "resistance mechanism in the species")
    elif resistant_pbp3:
        high = pbp3.call.startswith(("stage 2", "stage 3"))
        call, meaning = ("BLNAR" if high else "low-BLNAR"), (
            "beta-lactamase negative, ampicillin resistant: resistance rests on the "
            "altered PBP3, so ampicillin, amoxicillin and amoxicillin-clavulanate should "
            "all be reported as resistant whatever the measured MIC, and oral "
            "cephalosporin activity is reduced. " +
            ("The S385T substitution puts this in the BLNAR class proper rather than "
             "low-BLNAR" if high else
             "Without S385T this is the low-BLNAR class: the ampicillin MIC often sits "
             "at the breakpoint and a disc test can read it as susceptible, which is "
             "exactly why the genotype is reported"))
    else:
        call, meaning = "BLNAS", (
            "beta-lactamase negative, ampicillin susceptible: no acquired beta-lactamase "
            "and no resistant PBP3 pattern")
    return TypingResult(scheme="beta_lactam", call=call, lineage="-", alleles=alleles,
                        note=meaning)


# --------------------------------------------------------------------- MLST
def _call_mlst(assembly: Path, root: Path, work: Path, threads: int,
               mlst=None) -> TypingResult:
    """Seven-locus sequence type from the PubMLST *H. influenzae* scheme."""
    allele_db = root / "mlst" / "mlst.fna"
    profiles = root / "mlst" / "profiles.txt"
    # Hydra's own MLST call first: it types these loci already.
    profile = common.profile_from_call(mlst, "hinfluenzae", profiles,
                                       accepted=MLST_SCHEMES)
    if profile is None:
        if not allele_db.exists():
            return TypingResult(
                scheme="mlst", call="not determined",
                note=f"MLST data missing; run 'hydra db download typing-{KEY}', "
                     f"or install Hydra's pubmlst database so the loci are typed "
                     f"once in the MLST step")
        profile = common.call_mlst(assembly, allele_db, profiles, "hinfluenzae",
                                   workdir=work, threads=threads)
    alleles = dict(profile.alleles)
    note = profile.note
    lineage = profile.clonal_complex or "-"
    if profile.sequence_type != "-" and not profile.clonal_complex:
        note = "; ".join(part for part in [note, (
            f"ST{profile.sequence_type} is not assigned to a clonal complex in the PubMLST "
            f"profile table, so none is reported; the complexes it does publish are named "
            f"after their central type (for example 'ST-6 complex', the serotype b "
            f"division I lineage)")] if part)
    call = f"ST{profile.sequence_type}" if profile.sequence_type != "-" else "not determined"
    return TypingResult(scheme="mlst", call=call, lineage=lineage, alleles=alleles,
                        loci_found=profile.loci_found, loci_total=profile.loci_total,
                        note=note)
