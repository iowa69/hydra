"""*Streptococcus pyogenes* (group A streptococcus) typing.

What this module reports, in the order a clinical report wants it:

1. **emm type and subtype** from the 5' hypervariable region of *emm*, against
   the CDC *emm* sequence database, with the **emm cluster** and, where the
   cluster allows it, the **emm pattern** and its tissue-tropism grouping.
2. **The phage-encoded virulence repertoire**: the superantigens and DNases
   carried on prophages, separated from the chromosomal superantigens, each with
   the evidence it rests on; plus a description of the prophage regions the
   assembly resolves, and an explicit statement where it does not.
3. **The regulatory lesions that matter clinically**: *covR*/*covS* (*csrR*/
   *csrS*) and *rocA*, described as lesions - premature stop, truncation,
   internal deletion - rather than as presence or absence, because a
   loss-of-function *covS* is the finding, and *covS* is present in every group A
   streptococcus.
4. **The adhesin repertoire**: fibronectin-binding and collagen-binding
   adhesins, the FCT-region pilus genes and the Mga-regulon surface proteins.

Everything is called from an assembly with blastn and blastx only.

Three things this module deliberately does not do, because an assembly and the
rules available cannot support them honestly:

* It does not assign the **M1UK** lineage. M1UK is defined by a set of SNPs
  relative to M1global and is assigned phylogenetically; the module reports the
  *nga*/*slo* operon genes and their state, and says that the lineage call needs
  the SNP set, which is not implemented here.
* It does not report a **clonal complex**. The PubMLST *S. pyogenes* scheme
  carries no clonal-complex field, so there is no table to look one up in.
* It does not name **speA or smeZ allele numbers**. The allele sets are defined
  in the literature but are not published as a downloadable, versioned
  reference; the module reports the residue differences from the named reference
  protein instead, which is the evidence an allele number would have rested on.

See ``docs/typing/spyogenes.md`` for the rules, thresholds and bibliography.
"""

from __future__ import annotations

import re
from collections import Counter
from datetime import date
from pathlib import Path

from ...engines.blast import blast, check_db_exists
from ...records import TypingResult
from ...seqio import read_fasta
from ...utils import LOG, tempdir
from .base import OrganismDb, Reference
from .strep_common import (
    CONTIG_EDGE_BP,
    GeneHit,
    GeneSpec,
    PanelSource,
    contig_lengths,
    download,
    fetch_panel,
    make_blast_db,
    panel_manifest,
    screen_panel,
)

KEY = "spyogenes"
TITLE = "Streptococcus pyogenes"
SPECIES: tuple[str, ...] = ("Streptococcus pyogenes",)
MLST_SCHEMES: tuple[str, ...] = ("spyogenes",)

# --------------------------------------------------------------------- sources
#: The CDC trimmed *emm* subtype file. One record per named subtype, each the
#: 180 bases (5' hypervariable region) that the CDC protocol compares.
CDC_EMM_URL = ("https://ftp.cdc.gov/pub/infectious_diseases/biotech/tsemm/"
               "alltrimmed.tfa")

#: RefSeq assemblies the reference proteins are taken from. Each is a strain
#: whose genome paper is in ``REFERENCES``, so a reference protein can be traced
#: to a published, closed genome rather than to a draft of unknown provenance.
PANEL_SOURCES: dict[str, PanelSource] = {
    "SF370": PanelSource("GCF_000006785.2", "ASM678v2", "SF370",
                         "serotype M1; Ferretti 2001"),
    "MGAS315": PanelSource("GCF_000007425.1", "ASM742v1", "MGAS315",
                           "serotype M3; Beres 2002"),
    "MGAS8232": PanelSource("GCF_000007285.1", "ASM728v1", "MGAS8232",
                            "serotype M18; Smoot 2002"),
    "MGAS10394": PanelSource("GCF_000011665.1", "ASM1166v1", "MGAS10394",
                             "serotype M6; Banks 2004"),
    "MGAS5005": PanelSource("GCF_000011765.3", "ASM1176v2", "MGAS5005",
                            "serotype M1; Sumby 2005"),
    "MGAS6180": PanelSource("GCF_000012165.2", "ASM1216v2", "MGAS6180",
                            "serotype M28; Green 2005"),
    "Alab49": PanelSource("GCF_000230295.1", "ASM23029v1", "Alab49",
                          "serotype M53, emm pattern D; Bessen 2015"),
}

# ----------------------------------------------------------------- emm scheme
#: emmtyper's filters, which implement the CDC subtype-calling rule: an HSP is a
#: candidate subtype only if it differs from the reference by at most four
#: mismatches, at most two gap openings, and at most five bases of length.
EMM_MIN_IDENTITY = 95.0
EMM_MAX_MISMATCH = 4
EMM_MAX_GAPOPEN = 2
EMM_MAX_ALIGN_DIFF = 5
EMM_CULLING_LIMIT = 5
#: Hits further apart than this are separate loci, so they get separate answers.
#: *emm* sits in a tandem array with *mrp* and *enn*, which cross-match.
EMM_CLUSTER_DISTANCE = 500

#: Types that PHE/UKHSA classify as *emm*-like rather than true *emm*: the
#: sequence lies in the emm database but the locus is *mrp* or *enn*, or the
#: organism is a group C/G streptococcus. A call on one of these is flagged.
#: Taken from emmtyper, which is the tool this scheme is validated against.
EMM_LIKE: frozenset[str] = frozenset({
    "EMM51", "EMM134", "EMM138", "EMM149", "EMM156", "EMM159", "EMM164",
    "EMM167", "EMM170", "EMM174", "EMM202", "EMM205", "EMM236", "EMM240",
})

#: Cluster prefix -> emm pattern and the tissue-tropism grouping it carries.
#: Only the A-C, D and E cluster series map onto an emm pattern; a
#: single-protein cluster does not, and is reported as undetermined.
PATTERN_BY_CLUSTER: tuple[tuple[str, str, str], ...] = (
    ("A-C", "A-C", "throat-associated (emm pattern A-C)"),
    ("D", "D", "skin-associated (emm pattern D)"),
    ("E", "E", "generalist, throat and skin (emm pattern E)"),
)

_CLUSTER_IN_HEADER = re.compile(r"emm[- ]cluster\s+([A-Za-z0-9-]+)", re.I)
#: Minimum share of a type's subtype records that must agree before the cluster
#: annotation is used. The CDC descriptions are free text and a handful carry a
#: mistyped cluster; one stray label must not overturn a hundred consistent ones.
CLUSTER_AGREEMENT = 2 / 3


# ------------------------------------------------------------------ the panel
def _sag(name: str, symbols: tuple[str, ...], source: str, aa: int,
         location: str, meaning: str, product_re: str = "") -> GeneSpec:
    return GeneSpec(name=name, symbols=symbols, source=source, expect_aa=aa,
                    location=location, meaning=meaning, product_re=product_re,
                    product="superantigen" if name.startswith(("spe", "sme", "ssa"))
                    else "secreted virulence factor",
                    min_identity=85.0, min_coverage=70.0,
                    reference="Commons2014")


#: The superantigens and the phage-associated DNases. ``location`` is the
#: published genomic context of the gene, not something measured here: it is
#: what tells the reader whether a positive is expected to travel on a prophage.
SUPERANTIGENS: tuple[GeneSpec, ...] = (
    _sag("speA", ("spea",), "MGAS315", 251, "prophage",
         "streptococcal pyrogenic exotoxin A; the classic scarlet-fever toxin, "
         "carried on a prophage and strongly associated with invasive disease"),
    _sag("speC", ("spec",), "SF370", 235, "prophage",
         "streptococcal pyrogenic exotoxin C; prophage-encoded superantigen"),
    _sag("speH", ("speh",), "SF370", 236, "prophage",
         "prophage-encoded superantigen, usually paired with speI"),
    _sag("speI", ("spei",), "SF370", 225, "prophage",
         "prophage-encoded superantigen, usually paired with speH"),
    _sag("speK", ("spek",), "MGAS315", 259, "prophage",
         "prophage-encoded superantigen (also called speH-like/SpeK/SpeL group)"),
    _sag("speL", ("spel",), "MGAS8232", 262, "prophage",
         "prophage-encoded superantigen, described in serotype M18"),
    _sag("speM", ("spem",), "MGAS8232", 237, "prophage",
         "prophage-encoded superantigen, described in serotype M18"),
    _sag("speG", ("speg",), "SF370", 234, "chromosome",
         "chromosomal superantigen; present in almost every group A "
         "streptococcus, so its absence is the notable finding"),
    _sag("speJ", ("spej",), "SF370", 232, "chromosome",
         "chromosomal superantigen, characteristic of the emm1 lineage"),
    _sag("smeZ", ("smez",), "SF370", 233, "chromosome",
         "streptococcal mitogenic exotoxin Z; chromosomal, the most potent "
         "group A streptococcal superantigen, allelically variable"),
)

#: SpeB is a cysteine protease, not a superantigen, and is reported separately so
#: that a reader does not count it among the superantigens.
PROTEASES_AND_DNASES: tuple[GeneSpec, ...] = (
    _sag("speB", ("speb",), "SF370", 398, "chromosome",
         "SpeB / streptopain, a chromosomal cysteine protease (not a "
         "superantigen); loss of SpeB activity accompanies covRS lesions"),
    _sag("sdaB", ("sdab",), "SF370", 271, "chromosome",
         "DNase B / streptodornase B (mitogenic factor SpeF); chromosomal"),
    _sag("mf2", ("mf2",), "SF370", 252, "prophage",
         "prophage-encoded DNase Mf2 (Spd3)"),
    _sag("sda3", ("sda3",), "SF370", 266, "prophage",
         "prophage-encoded streptodornase Sda3"),
    _sag("sda1", ("sda1",), "MGAS5005", 390, "prophage",
         "prophage-encoded DNase Sda1; degrades neutrophil extracellular traps "
         "and is a determinant of the invasive M1T1 phenotype"),
    _sag("sdaD", ("sda",), "MGAS8232", 332, "prophage",
         "prophage-encoded streptodornase A (SdaD/SdaA group)"),
)

#: The two-component system and the regulators whose lesions change phenotype.
REGULATORS: tuple[GeneSpec, ...] = (
    GeneSpec("covR", ("covr", "csrr"), "SF370", "response regulator CovR/CsrR",
             location="regulator", expect_aa=228, min_identity=90.0,
             min_coverage=50.0, strict_ends=True, reference="Federle1999",
             meaning="represses capsule, streptokinase, SpyCEP and the "
                     "nga/slo operon; a lesion de-represses them"),
    GeneSpec("covS", ("covs", "csrs"), "SF370",
             "sensor histidine kinase CovS/CsrS",
             location="regulator", expect_aa=500, min_identity=90.0,
             min_coverage=50.0, strict_ends=True, reference="Sumby2006",
             meaning="loss of function switches the isolate to the invasive, "
                     "SpeB-negative, hypercapsular phenotype selected during "
                     "invasive infection"),
    GeneSpec("rocA", ("roca",), "SF370", "regulator of cov, RocA",
             location="regulator", expect_aa=451, min_identity=90.0,
             min_coverage=50.0, strict_ends=True, reference="Biswas2003",
             meaning="a partner of CovRS; naturally truncated in some "
                     "serotypes, which raises CovR-repressed virulence factors"),
    GeneSpec("mga", ("mga",), "SF370", "multigene regulator Mga",
             location="regulator", expect_aa=529, min_identity=85.0,
             min_coverage=60.0, strict_ends=True, reference="Hondorp2007",
             meaning="activates the emm-region surface protein genes"),
    GeneSpec("ropB", ("ropb",), "SF370", "RopB / Rgg transcriptional regulator",
             location="regulator", expect_aa=280, min_identity=85.0,
             min_coverage=60.0, strict_ends=True, reference="Hondorp2007",
             meaning="required for SpeB production; lesions give a "
                     "SpeB-negative phenotype independent of covRS"),
)

#: The nga/ifs/slo operon. The M1UK lineage carries promoter variants that raise
#: NADase and SLO expression; this module reports the genes and their state and
#: leaves the lineage call to a SNP-based method.
NGA_OPERON: tuple[GeneSpec, ...] = (
    GeneSpec("nga", ("nga",), "SF370", "NAD+ glycohydrolase (NADase)",
             location="chromosome", expect_aa=447, min_identity=85.0,
             min_coverage=70.0, reference="Zhu2015",
             meaning="NADase; high expression, together with SLO, marks the "
                     "invasive emm1 and emm89 lineages"),
    GeneSpec("slo", ("slo",), "SF370", "streptolysin O",
             location="chromosome", expect_aa=571, min_identity=85.0,
             min_coverage=70.0, reference="Zhu2015",
             meaning="cholesterol-dependent cytolysin, co-transcribed with nga"),
    GeneSpec("ifs", (), "MGAS315", "NADase inhibitor Ifs",
             location="chromosome", expect_aa=161, min_identity=80.0,
             min_coverage=60.0, product_re=r"NAD glycohydrolase inhibitor",
             reference="Zhu2015",
             meaning="the intracellular inhibitor that protects the cell from "
                     "its own NADase"),
)

ADHESINS: tuple[GeneSpec, ...] = (
    GeneSpec("emm", ("emm",), "SF370", "M protein", location="surface",
             expect_aa=484, min_identity=70.0, min_coverage=40.0,
             reference="Beall1996",
             meaning="the M protein itself; the locus the emm type is read from"),
    GeneSpec("scpA", ("scpa/b", "scpa"), "SF370", "C5a peptidase ScpA",
             location="surface", expect_aa=1181, min_identity=85.0,
             min_coverage=70.0, reference="Walker2014",
             meaning="cleaves C5a and so blunts neutrophil recruitment"),
    GeneSpec("enn", ("enn",), "Alab49", "M-related protein Enn",
             location="surface", expect_aa=356, min_identity=70.0,
             min_coverage=40.0, reference="Walker2014",
             meaning="an emm-region paralogue; its presence marks the "
                     "three-gene emm region of emm pattern D and E strains"),
    GeneSpec("fbp54", ("fbp54",), "SF370", "fibronectin-binding protein Fbp54",
             location="surface", expect_aa=550, min_identity=85.0,
             min_coverage=70.0, reference="Walker2014",
             meaning="fibronectin-binding adhesin, present in nearly all "
                     "isolates"),
    GeneSpec("prtF1", ("prtf1", "sfbi"), "MGAS10394",
             "fibronectin-binding protein PrtF1/SfbI", location="surface",
             expect_aa=615, min_identity=80.0, min_coverage=60.0,
             reference="Kreikemeyer2004",
             meaning="fibronectin-binding adhesin of the FCT region; mediates "
                     "epithelial invasion"),
    GeneSpec("fbaB", ("fbab", "prtf2"), "MGAS8232",
             "fibronectin-binding protein FbaB/PrtF2/SfbII",
             location="surface", expect_aa=696, min_identity=80.0,
             min_coverage=60.0, reference="Kreikemeyer2004",
             meaning="the second fibronectin-binding adhesin; associated with "
                     "invasive isolates"),
    GeneSpec("scl1", ("scl1", "scla"), "MGAS6180",
             "collagen-like surface protein Scl1/SclA", location="surface",
             expect_aa=455, min_identity=70.0, min_coverage=40.0,
             reference="Lukomski2000",
             meaning="collagen-like adhesin; the repeat region varies in length "
                     "between lineages, so coverage below 100% is expected"),
    GeneSpec("scl2", ("slc2", "scl2", "sclb"), "MGAS315",
             "collagen-like surface protein Scl2/SclB", location="surface",
             expect_aa=573, min_identity=70.0, min_coverage=40.0,
             reference="Lukomski2000",
             meaning="the second collagen-like adhesin; also repeat-variable"),
    GeneSpec("cpa", ("cpa",), "SF370", "FCT-region pilus tip adhesin Cpa",
             location="surface", expect_aa=757, min_identity=75.0,
             min_coverage=50.0, reference="Falugi2008",
             meaning="collagen-binding pilus tip adhesin of the FCT region"),
    GeneSpec("fctA", ("fcta", "tee1", "tee3", "tee6", "tee18"), "SF370",
             "FCT-region pilus backbone (T antigen)", location="surface",
             expect_aa=340, min_identity=70.0, min_coverage=40.0,
             reference="Falugi2008",
             meaning="the pilus backbone protein, which is the T antigen; its "
                     "sequence type is the basis of T typing, which this "
                     "module does not call"),
    GeneSpec("fctB", ("fctb",), "SF370", "FCT-region minor pilin FctB",
             location="surface", expect_aa=215, min_identity=75.0,
             min_coverage=50.0, reference="Falugi2008",
             meaning="minor pilin of the FCT region"),
    GeneSpec("srtB", ("srtb",), "MGAS6180", "class B pilus sortase SrtB",
             location="surface", expect_aa=241, min_identity=80.0,
             min_coverage=60.0, reference="Falugi2008",
             meaning="polymerises the pilus backbone"),
    GeneSpec("grab", ("grab",), "SF370",
             "alpha-2-macroglobulin-binding protein GRAB", location="surface",
             expect_aa=217, min_identity=80.0, min_coverage=60.0,
             reference="Walker2014",
             meaning="binds alpha-2-macroglobulin and protects SpeB from "
                     "host protease inhibitors"),
    GeneSpec("ska", ("ska",), "SF370", "streptokinase", location="chromosome",
             expect_aa=440, min_identity=85.0, min_coverage=70.0,
             reference="Walker2014",
             meaning="plasminogen activator; the allele determines host "
                     "specificity"),
    GeneSpec("scpC", ("scpc",), "SF370", "SpyCEP / ScpC IL-8 protease",
             location="chromosome", expect_aa=1647, min_identity=85.0,
             min_coverage=70.0, reference="Walker2014",
             meaning="cleaves IL-8 and so prevents neutrophil recruitment; a "
                     "vaccine candidate"),
    GeneSpec("ideS", ("ides",), "SF370", "IdeS / Mac IgG protease",
             location="chromosome", expect_aa=339, min_identity=85.0,
             min_coverage=70.0, reference="Walker2014",
             meaning="cleaves IgG and blocks opsonophagocytosis"),
    GeneSpec("hasA", ("hasa",), "SF370", "hyaluronan synthase HasA",
             location="chromosome", expect_aa=419, min_identity=85.0,
             min_coverage=70.0, reference="Walker2014",
             meaning="capsule synthase; an intact has operon is needed for the "
                     "hyaluronic acid capsule, and emm4 and emm22 lack it"),
    GeneSpec("sic", ("sic",), "SF370",
             "streptococcal inhibitor of complement Sic", location="surface",
             expect_aa=308, min_identity=70.0, min_coverage=40.0,
             reference="Walker2014",
             meaning="complement inhibitor restricted to emm1, emm57 and "
                     "emm12; hypervariable, so identity below 100% is normal"),
)

#: Prophage structural and lysogeny markers, used to describe prophage regions.
#: These are not a published prophage caller; they are the gene classes a
#: prophage must carry, so that a region can be described as a prophage with the
#: evidence shown. Each class is listed with the function it stands for.
PHAGE_MARKERS: tuple[GeneSpec, ...] = (
    GeneSpec("phage_integrase", (), "MGAS315", "phage integrase",
             location="phage", expect_aa=191, min_identity=55.0,
             min_coverage=40.0, product_re=r"^site-specific integrase$",
             reference="Beres2002",
             meaning="lysogeny module; marks the attachment end of a prophage"),
    GeneSpec("phage_terminase", (), "SF370", "phage terminase large subunit",
             location="phage", expect_aa=403, min_identity=55.0,
             min_coverage=40.0, product_re=r"terminase large subunit",
             reference="Beres2002", meaning="DNA packaging module"),
    GeneSpec("phage_portal", (), "SF370", "phage portal protein",
             location="phage", expect_aa=500, min_identity=55.0,
             min_coverage=40.0, product_re=r"^phage portal protein$",
             reference="Beres2002", meaning="head-tail connector"),
    GeneSpec("phage_capsid", (), "SF370", "phage major capsid protein",
             location="phage", expect_aa=303, min_identity=55.0,
             min_coverage=40.0, product_re=r"^phage major capsid protein$",
             reference="Beres2002", meaning="head morphogenesis module"),
    GeneSpec("phage_tape_measure", (), "SF370", "phage tape measure protein",
             location="phage", expect_aa=1086, min_identity=55.0,
             min_coverage=40.0, product_re=r"^tape measure protein$",
             reference="Beres2002", meaning="tail module"),
    GeneSpec("phage_holin", (), "SF370", "phage holin",
             location="phage", expect_aa=75, min_identity=55.0,
             min_coverage=40.0, product_re=r"^phage holin$",
             reference="Beres2002", meaning="lysis module"),
    GeneSpec("hylP", ("hylp",), "SF370", "phage hyaluronidase HylP",
             location="phage", expect_aa=337, min_identity=55.0,
             min_coverage=40.0, reference="Beres2002",
             meaning="phage-associated hyaluronidase; the tail-spike lysin that "
                     "lets the phage through the capsule"),
    GeneSpec("phage_minor_capsid", (), "SF370", "phage minor capsid protein",
             location="phage", expect_aa=497, min_identity=55.0,
             min_coverage=40.0, product_re=r"^phage minor capsid protein$",
             reference="Beres2002", meaning="head morphogenesis module"),
    GeneSpec("phage_tail", (), "SF370", "phage tail protein",
             location="phage", expect_aa=785, min_identity=55.0,
             min_coverage=35.0, product_re=r"^phage tail protein$",
             reference="Beres2002", meaning="tail module"),
    GeneSpec("phage_tail_spike", (), "SF370", "phage tail spike protein",
             location="phage", expect_aa=714, min_identity=55.0,
             min_coverage=35.0, product_re=r"^phage tail spike protein$",
             reference="Beres2002", meaning="host recognition"),
    GeneSpec("phage_head_morphogenesis", (), "SF370",
             "phage head morphogenesis protein", location="phage",
             expect_aa=541, min_identity=55.0, min_coverage=40.0,
             product_re=r"^phage head morphogenesis protein$",
             reference="Beres2002", meaning="head morphogenesis module"),
    GeneSpec("phage_repressor", (), "SF370", "phage repressor protein",
             location="phage", expect_aa=207, min_identity=55.0,
             min_coverage=40.0, product_re=r"^phage repressor protein$",
             reference="Beres2002",
             meaning="lysogeny module; maintains the prophage state"),
    GeneSpec("phage_recombination", (), "SF370",
             "phage recombination protein Bet", location="phage",
             expect_aa=263, min_identity=55.0, min_coverage=40.0,
             product_re=r"^phage recombination protein Bet$",
             reference="Beres2002", meaning="recombination module"),
    GeneSpec("phage_major_tail", (), "SF370",
             "phage major tail protein", location="phage", expect_aa=188,
             min_identity=55.0, min_coverage=40.0,
             product_re=r"^phage major tail protein",
             reference="Beres2002", meaning="tail module"),
)

PANEL: tuple[GeneSpec, ...] = (
    SUPERANTIGENS + PROTEASES_AND_DNASES + REGULATORS + NGA_OPERON
    + ADHESINS + PHAGE_MARKERS
)

#: Classes a region must carry before it is described as a prophage region.
PHAGE_CLASSES_FOR_REGION = 3
#: How far apart two phage-module genes can sit and still be one prophage.
#: Group A streptococcal prophages are 33-48 kb (Beres 2002, Banks 2002), so a
#: 20 kb reach joins the modules of one element without bridging two.
PHAGE_LINK_BP = 20000

# -------------------------------------------------------------------- metadata
DATABASES: tuple[OrganismDb, ...] = (
    OrganismDb(
        name="emm",
        source=CDC_EMM_URL,
        version="retrieved at fetch time; the file carries no internal version, "
                "so fetch() records the HTTP Last-Modified date in emm.version",
        licence="No licence statement accompanies the file. CDC's own terms "
                "state: \"Unless otherwise noted, all materials on this site "
                "are in the public domain and may be used and reprinted without "
                "special permission.\"",
        redistribution="permitted",
        reference="CDCstreplab",
        notes="One record per named emm subtype, each the 180-base 5' "
              "hypervariable region the CDC protocol compares. About 0.9 MB. "
              "The emm cluster is parsed out of the CDC record descriptions.",
    ),
    OrganismDb(
        name="panel",
        source=f"{'https://ftp.ncbi.nlm.nih.gov/genomes/all'} (RefSeq feature "
               f"tables and protein FASTA for "
               f"{', '.join(s.accession for s in PANEL_SOURCES.values())})",
        version="pinned assembly accession.version per strain, see PANEL_SOURCES",
        licence="NCBI: \"Please be aware that data in the Database may contain "
                "errors ... NLM does not claim copyright on the data in the "
                "Database\"; RefSeq records are free to use and redistribute.",
        redistribution="permitted",
        reference="OLeary2016",
        notes="Reference proteins for the virulence, regulatory and adhesin "
              "panel, selected from RefSeq annotation by gene symbol. The built "
              "database is about 60 kB; the download is about 10 MB and is "
              "cached under db_dir/spyogenes/refseq/.",
    ),
)

REFERENCES: tuple[Reference, ...] = (
    Reference(
        key="Beall1996",
        authors="Beall B, Facklam R, Thompson T",
        year=1996,
        title="Sequencing emm-specific PCR products for routine and accurate "
              "typing of group A streptococci",
        journal="Journal of Clinical Microbiology 34:953-958",
        doi="10.1128/jcm.34.4.953-958.1996",
        pmid="8815115",
        supports="the emm typing scheme: typing from the 5' variable region of "
                 "emm rather than by M serotyping",
    ),
    Reference(
        key="Beall1998",
        authors="Beall B, Facklam RR, Elliott JA, et al.",
        year=1998,
        title="Streptococcal emm types associated with T-agglutination types and "
              "the use of conserved emm gene restriction fragment patterns for "
              "subtyping group A streptococci",
        journal="Journal of Medical Microbiology 47:893-898",
        doi="10.1099/00222615-47-10-893",
        pmid="9788813",
        supports="emm subtyping and the relationship between emm type and "
                 "T-agglutination pattern",
    ),
    Reference(
        key="CDCstreplab",
        authors="Centers for Disease Control and Prevention, Streptococcus "
                "Laboratory",
        year=2026,
        title="Streptococcus laboratory: emm typing protocol and emm sequence "
              "database (alltrimmed.tfa)",
        journal="CDC, Atlanta; https://www.cdc.gov/strep-lab/php/protocols/",
        supports="the emm sequence database, the trimmed 180-base hypervariable "
                 "region, and the type/subtype naming the module reports",
    ),
    Reference(
        key="SandersonSmith2014",
        authors="Sanderson-Smith M, De Oliveira DM, Guglielmini J, et al.",
        year=2014,
        title="A systematic and functional classification of Streptococcus "
              "pyogenes that serves as a new tool for molecular typing and "
              "vaccine development",
        journal="The Journal of Infectious Diseases 210:1325-1338",
        doi="10.1093/infdis/jiu260",
        pmid="24799598",
        supports="the emm cluster system (clusters A-C1 to A-C5, D1 to D5, E1 to "
                 "E6, single-protein clusters and clade Y) that the module reports",
    ),
    Reference(
        key="Bessen2010",
        authors="Bessen DE, Lizano S",
        year=2010,
        title="Tissue tropisms in group A streptococcal infections",
        journal="Future Microbiology 5:623-638",
        doi="10.2217/fmb.10.28",
        pmid="20353302",
        supports="the emm pattern (A-C, D, E) tissue-tropism grouping the module "
                 "derives from the emm cluster",
    ),
    Reference(
        key="Frost2020",
        authors="Frost HR, Davies MR, Velusamy S, et al.",
        year=2020,
        title="Updated emm-typing protocol for Streptococcus pyogenes",
        journal="Clinical Microbiology and Infection 26:946.e5-946.e8",
        doi="10.1016/j.cmi.2020.02.026",
        pmid="32120034",
        supports="the current reference-laboratory emm typing protocol, and the "
                 "definition of emm-like loci that stops an emm-like type being "
                 "reported as the emm type",
    ),
    Reference(
        key="emmtyper",
        authors="Mu A, Kwong JC, Isles NS, et al. (Microbiological Diagnostic "
                "Unit Public Health Laboratory)",
        year=2019,
        title="emmtyper: emm automatic isolate labeller, version 0.2.0",
        journal="Software; https://github.com/MDUPHL/emmtyper (GPL-3.0)",
        supports="the reference tool this module's emm calls were validated "
                 "against, and the filter thresholds it applies - at most 4 "
                 "mismatches, at most 2 gap openings, at most 5 bases of "
                 "alignment-to-subject length difference, 95% identity and 500 "
                 "bp position clustering - which this module reimplements. No "
                 "journal article describes the tool; it is cited as software",
    ),
    Reference(
        key="Commons2014",
        authors="Commons RJ, Smeesters PR, Proft T, et al.",
        year=2014,
        title="Streptococcal superantigens: categorization and clinical "
              "associations",
        journal="Trends in Molecular Medicine 20:48-62",
        doi="10.1016/j.molmed.2013.10.004",
        pmid="24210845",
        supports="which superantigens are chromosomal and which prophage-encoded, "
                 "and the clinical associations the module prints",
    ),
    Reference(
        key="Beres2002",
        authors="Beres SB, Sylva GL, Barbian KD, et al.",
        year=2002,
        title="Genome sequence of a serotype M3 strain of group A Streptococcus: "
              "phage-encoded toxins, the high-virulence phenotype, and clone "
              "emergence",
        journal="Proceedings of the National Academy of Sciences USA "
                "99:10078-10083",
        doi="10.1073/pnas.152298499",
        pmid="12122206",
        supports="the prophage-encoded toxin content of MGAS315 and the prophage "
                 "module structure (integrase, packaging, head, tail, lysis) the "
                 "region description rests on",
    ),
    Reference(
        key="Ferretti2001",
        authors="Ferretti JJ, McShan WM, Ajdic D, et al.",
        year=2001,
        title="Complete genome sequence of an M1 strain of Streptococcus pyogenes",
        journal="Proceedings of the National Academy of Sciences USA 98:4658-4663",
        doi="10.1073/pnas.071559398",
        pmid="11296296",
        supports="SF370, the source of most reference proteins, and its four "
                 "prophages with their toxin cargo",
    ),
    Reference(
        key="Smoot2002",
        authors="Smoot JC, Barbian KD, Van Gompel JJ, et al.",
        year=2002,
        title="Genome sequence and comparative microarray analysis of serotype "
              "M18 group A Streptococcus strains associated with acute rheumatic "
              "fever outbreaks",
        journal="Proceedings of the National Academy of Sciences USA "
                "99:4668-4673",
        doi="10.1073/pnas.062526099",
        pmid="11917108",
        supports="MGAS8232 and the speL/speM reference proteins taken from it",
    ),
    Reference(
        key="Banks2004",
        authors="Banks DJ, Porcella SF, Barbian KD, et al.",
        year=2004,
        title="Progress toward characterization of the group A Streptococcus "
              "metagenome: complete genome sequence of a macrolide-resistant "
              "serotype M6 strain",
        journal="The Journal of Infectious Diseases 190:727-738",
        doi="10.1086/422697",
        pmid="15272401",
        supports="MGAS10394, the source of the prtF1/SfbI reference protein",
    ),
    Reference(
        key="Sumby2005",
        authors="Sumby P, Porcella SF, Madrigal AG, et al.",
        year=2005,
        title="Evolutionary origin and emergence of a highly successful clone of "
              "serotype M1 group A Streptococcus involved multiple horizontal "
              "gene transfer events",
        journal="The Journal of Infectious Diseases 192:771-782",
        doi="10.1086/432514",
        pmid="16088826",
        supports="MGAS5005 and the sda1 reference protein taken from it",
    ),
    Reference(
        key="Sumby2006",
        authors="Sumby P, Whitney AR, Graviss EA, et al.",
        year=2006,
        title="Genome-wide analysis of group A streptococci reveals a mutation "
              "that modulates global phenotype and disease specificity",
        journal="PLoS Pathogens 2:e5",
        doi="10.1371/journal.ppat.0020005",
        pmid="16446783",
        supports="the covS lesion in MGAS5005 and the clinical meaning of a "
                 "covRS lesion, which is the rule the regulator call implements",
    ),
    Reference(
        key="Green2005",
        authors="Green NM, Zhang S, Porcella SF, et al.",
        year=2005,
        title="Genome sequence of a serotype M28 strain of group A Streptococcus: "
              "potential new insights into puerperal sepsis and bacterial "
              "disease specificity",
        journal="The Journal of Infectious Diseases 192:760-770",
        doi="10.1086/430618",
        pmid="16088825",
        supports="MGAS6180, the source of the scl1 and srtB reference proteins",
    ),
    Reference(
        key="Bessen2015",
        authors="Bessen DE, McShan WM, Nguyen SV, et al.",
        year=2015,
        title="Molecular epidemiology and genomics of group A Streptococcus",
        journal="Infection, Genetics and Evolution 33:393-418",
        doi="10.1016/j.meegid.2014.10.011",
        pmid="25460818",
        supports="the Alab49 reference genome, the emm pattern framework and the "
                 "emm-region gene arrangements the enn call rests on",
    ),
    Reference(
        key="Federle1999",
        authors="Federle MJ, McIver KS, Scott JR",
        year=1999,
        title="A response regulator that represses transcription of several "
              "virulence operons in the group A streptococcus",
        journal="Journal of Bacteriology 181:3649-3657",
        doi="10.1128/JB.181.12.3649-3657.1999",
        pmid="10368137",
        supports="CovR/CsrR as a repressor of capsule and other virulence "
                 "operons, which is what a covR lesion de-represses",
    ),
    Reference(
        key="Biswas2003",
        authors="Biswas I, Scott JR",
        year=2003,
        title="Identification of rocA, a positive regulator of covR expression in "
              "the group A streptococcus",
        journal="Journal of Bacteriology 185:3081-3090",
        doi="10.1128/JB.185.10.3081-3090.2003",
        pmid="12730168",
        supports="rocA as a positive regulator of covR, and therefore the meaning "
                 "of a rocA truncation",
    ),
    Reference(
        key="Hondorp2007",
        authors="Hondorp ER, McIver KS",
        year=2007,
        title="The Mga virulence regulon: infection where the grass is greener",
        journal="Molecular Microbiology 66:1056-1065",
        doi="10.1111/j.1365-2958.2007.06006.x",
        pmid="18001346",
        supports="Mga and RopB as regulators of the emm-region surface proteins "
                 "and of SpeB",
    ),
    Reference(
        key="Zhu2015",
        authors="Zhu L, Olsen RJ, Nasser W, et al.",
        year=2015,
        title="A molecular trigger for intercontinental epidemics of group A "
              "Streptococcus",
        journal="The Journal of Clinical Investigation 125:3545-3559",
        doi="10.1172/JCI82478",
        pmid="26258415",
        supports="the nga/ifs/slo operon, the promoter variants that raise NADase "
                 "and SLO production, and their role in epidemic emm89 and emm1 "
                 "lineages",
    ),
    Reference(
        key="Lynskey2019",
        authors="Lynskey NN, Jauneikaite E, Li HK, et al.",
        year=2019,
        title="Emergence of dominant toxigenic M1T1 Streptococcus pyogenes clone "
              "during increased scarlet fever activity in England: a "
              "population-based molecular epidemiological study",
        journal="The Lancet Infectious Diseases 19:1209-1218",
        doi="10.1016/S1473-3099(19)30446-3",
        pmid="31519541",
        supports="the definition of the M1UK lineage, which is why this module "
                 "reports the nga/slo operon but does not assign the lineage",
    ),
    Reference(
        key="Walker2014",
        authors="Walker MJ, Barnett TC, McArthur JD, et al.",
        year=2014,
        title="Disease manifestations and pathogenic mechanisms of group A "
              "Streptococcus",
        journal="Clinical Microbiology Reviews 27:264-301",
        doi="10.1128/CMR.00101-13",
        pmid="24696436",
        supports="the clinical meaning printed for the secreted virulence "
                 "factors and surface adhesins in the panel",
    ),
    Reference(
        key="Kreikemeyer2004",
        authors="Kreikemeyer B, Klenk M, Podbielski A",
        year=2004,
        title="The intracellular status of Streptococcus pyogenes: role of "
              "extracellular matrix-binding proteins and their regulation",
        journal="International Journal of Medical Microbiology 294:177-188",
        doi="10.1016/j.ijmm.2004.06.017",
        pmid="15493828",
        supports="the fibronectin-binding adhesins PrtF1/SfbI and FbaB/PrtF2 and "
                 "their role in epithelial invasion",
    ),
    Reference(
        key="Lukomski2000",
        authors="Lukomski S, Nakashima K, Abdi I, et al.",
        year=2000,
        title="Identification and characterization of the scl gene encoding a "
              "group A Streptococcus extracellular protein virulence factor with "
              "similarity to human collagen",
        journal="Infection and Immunity 68:6542-6553",
        doi="10.1128/IAI.68.12.6542-6553.2000",
        pmid="11083763",
        supports="the collagen-like surface proteins Scl1/SclA and Scl2/SclB",
    ),
    Reference(
        key="Falugi2008",
        authors="Falugi F, Zingaretti C, Pinto V, et al.",
        year=2008,
        title="Sequence variation in group A Streptococcus pili and association of "
              "pilus backbone types with Lancefield T serotypes",
        journal="The Journal of Infectious Diseases 198:1834-1841",
        doi="10.1086/593176",
        pmid="18928376",
        supports="the FCT-region pilus genes, and the fact that the pilus "
                 "backbone is the T antigen",
    ),
    Reference(
        key="OLeary2016",
        authors="O'Leary NA, Wright MW, Brister JR, et al.",
        year=2016,
        title="Reference sequence (RefSeq) database at NCBI: current status, "
              "taxonomic expansion, and functional annotation",
        journal="Nucleic Acids Research 44:D733-D745",
        doi="10.1093/nar/gkv1189",
        pmid="26553804",
        supports="RefSeq, the source of every reference protein in the panel",
    ),
)


# ------------------------------------------------------------------- the fetch
def fetch(db_dir: Path) -> None:
    """Download and build this module's data into ``db_dir/spyogenes/``.

    Idempotent: an existing, indexed emm database with a cluster table and a
    panel manifest covering every gene is left alone.
    """
    root = Path(db_dir) / KEY
    root.mkdir(parents=True, exist_ok=True)

    emm_fasta = root / "emm.tfa"
    clusters = root / "emm_clusters.tsv"
    version = root / "emm.version"
    if not (emm_fasta.exists() and clusters.exists()):
        LOG.info("%s: downloading the CDC emm database", KEY)
        download(CDC_EMM_URL, emm_fasta)
        stamp = _last_modified(CDC_EMM_URL) or date.today().isoformat()
        version.write_text(f"{stamp}\t{CDC_EMM_URL}\n", encoding="utf-8")
        _write_cluster_table(emm_fasta, clusters)
    make_blast_db(emm_fasta, "nucl")

    LOG.info("%s: building the virulence panel from RefSeq annotation", KEY)
    fetch_panel(PANEL, PANEL_SOURCES, root, name="panel")


def _last_modified(url: str) -> str:
    """The HTTP Last-Modified date of *url*, which is the only version it has."""
    import urllib.request  # noqa: PLC0415 - only fetch() needs it

    request = urllib.request.Request(url, method="HEAD",
                                     headers={"User-Agent": "hydra-amr"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return (response.headers.get("Last-Modified") or "").strip()
    except OSError:
        return ""


def _write_cluster_table(fasta: Path, dest: Path) -> None:
    """Derive an emm type -> emm cluster table from the CDC record descriptions.

    The CDC curators annotate the emm cluster inside the free-text description of
    many subtype records ("emm1.10 (emm cluster A-C3) subtype definitive
    sequence"). Collecting those per type gives a table from the primary source
    rather than from a third-party copy, but it also inherits the occasional
    mistyped label, so a type's cluster is only accepted when at least two thirds
    of its annotated subtypes agree. The counts are written out so the decision
    can be audited.
    """
    votes: dict[str, Counter] = {}
    for header, _seq in read_fasta(fasta):
        token = header.split()[0]
        stem = token.rsplit(".", 1)[0].upper()
        match = _CLUSTER_IN_HEADER.search(header)
        if match:
            votes.setdefault(stem, Counter())[_normalise_cluster(match.group(1))] += 1
    lines = ["#emm_type\tcluster\tagreeing\ttotal_annotated"]
    for stem in sorted(votes):
        cluster, best, total = _resolve_cluster(votes[stem])
        lines.append(f"{stem}\t{cluster}\t{best}\t{total}")
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    LOG.info("%s: emm cluster table covers %d emm types", KEY, len(votes))


_SINGLE_PROTEIN = re.compile(r"^M\d+$")


def _resolve_cluster(counts: Counter) -> tuple[str, int, int]:
    """Decide one emm type's cluster from the votes its subtype records cast.

    Two labels in the CDC descriptions are not in conflict even when they
    disagree: a single-protein cluster (``M5``) and ``clade Y`` are nested, since
    clade Y is the group the single-protein clusters were collected into
    (Sanderson-Smith 2014). Where those are the only two labels, the specific one
    wins and the nesting is noted downstream. Everything else needs a two-thirds
    majority, because a free-text description can carry a mistyped label and one
    stray vote must not overturn a hundred consistent ones.
    """
    total = sum(counts.values())
    if not total:
        return "", 0, 0
    labels = set(counts)
    nested = {label for label in labels
              if _SINGLE_PROTEIN.match(label) or label == "clade Y"}
    if labels == nested and len(labels) > 1:
        specific = sorted((label for label in labels if _SINGLE_PROTEIN.match(label)),
                          key=lambda label: -counts[label])
        if specific:
            return specific[0], counts[specific[0]], total
        return "clade Y", counts["clade Y"], total
    cluster, best = counts.most_common(1)[0]
    if not cluster or best / total < CLUSTER_AGREEMENT:
        return "", best, total
    return cluster, best, total


def _normalise_cluster(label: str) -> str:
    """Tidy a free-text cluster label into the published nomenclature."""
    text = label.strip()
    if text.lower() in ("clade", "y"):
        return "clade Y"
    if text.lower() == "outlier":
        return "outlier"
    upper = text.upper()
    if upper.startswith("A-C") or upper.startswith("A-D"):
        return upper
    if re.fullmatch(r"[DE]\d+", upper):
        return upper
    if re.fullmatch(r"M\d+", upper):
        return upper
    return text


def _load_cluster_table(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path.exists():
        return out
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1]:
                out[parts[0].upper()] = parts[1]
    return out


# --------------------------------------------------------------------- the run
def run(assembly: Path, db_dir: Path, *, threads: int = 1,
        workdir: Path | None = None) -> list[TypingResult]:
    """Type one assembly. One :class:`TypingResult` per scheme."""
    assembly = Path(assembly)
    root = Path(db_dir) / KEY
    if workdir is not None:
        work = Path(workdir)
        work.mkdir(parents=True, exist_ok=True)
        return _run_in(assembly, root, work, threads)
    with tempdir(prefix="hydra.spyogenes.") as temp:
        return _run_in(assembly, root, Path(temp), threads)


def _run_in(assembly: Path, root: Path, work: Path,
            threads: int) -> list[TypingResult]:
    results = [_call_emm(assembly, root, work, threads)]
    panel_fasta = root / "panel.faa"
    if not panel_fasta.exists():
        note = (f"the virulence panel is not installed under {root}; "
                f"run the module's fetch() first")
        for scheme in ("emm_virulome", "spy_prophage", "spy_regulators",
                       "spy_adhesins"):
            results.append(TypingResult(scheme=scheme, call="not determined",
                                        note=note))
        return results
    check_db_exists(panel_fasta, "prot")
    hits, every = screen_panel(assembly, panel_fasta, work, PANEL,
                               threads=threads, tag=KEY)
    manifest = panel_manifest(root, "panel")
    lengths = contig_lengths(assembly)
    results.append(_call_superantigens(hits, manifest))
    results.append(_call_prophage(hits, every, lengths))
    results.append(_call_regulators(hits))
    results.append(_call_adhesins(hits))
    return results


# ------------------------------------------------------------------ emm typing
class _EmmHit:
    """One candidate emm subtype, scored the way the reference tool scores it."""

    __slots__ = ("allele", "contig", "start", "end", "identity", "mismatch",
                 "gapopen", "align_len", "ref_len", "bitscore", "score")

    def __init__(self, hsp) -> None:
        self.allele = hsp.sseqid
        self.contig = hsp.qseqid
        self.start = hsp.q_lo
        self.end = hsp.q_hi
        self.identity = hsp.pident
        self.mismatch = hsp.mismatch
        self.gapopen = hsp.gapopen
        self.align_len = hsp.length
        self.ref_len = hsp.slen
        self.bitscore = hsp.bitscore
        # emmtyper's score: a perfect, full-length, gap-free match scores 100.
        self.score = self.identity - (self.gapopen
                                      + abs(self.align_len - self.ref_len))

    @property
    def passes(self) -> bool:
        return (self.mismatch <= EMM_MAX_MISMATCH
                and self.gapopen <= EMM_MAX_GAPOPEN
                and abs(self.ref_len - self.align_len) <= EMM_MAX_ALIGN_DIFF)

    @property
    def exact(self) -> bool:
        return self.score >= 100.0

    @property
    def subtype(self) -> str:
        """``EMM1.6`` -> ``emm1.6``; the lower-case form is the one in use."""
        name = self.allele
        if name[:3].upper() == "EMM":
            return "emm" + name[3:]
        # The CDC file stores record ids in upper case but writes the names in
        # their conventional form in the description: stG62647.0, stL1929.3,
        # stC74a.23. These are the group C and G streptococcal st types.
        if name[:2].upper() == "ST":
            return "st" + name[2:]
        return name[0].lower() + name[1:]

    @property
    def emm_type(self) -> str:
        return self.subtype.rsplit(".", 1)[0]


def _call_emm(assembly: Path, root: Path, work: Path,
              threads: int) -> TypingResult:
    fasta = root / "emm.tfa"
    if not fasta.exists():
        return TypingResult(
            scheme="emm", call="not determined",
            note=f"the CDC emm database is not installed under {root}; "
                 f"run the module's fetch() first")
    check_db_exists(fasta, "nucl")
    work.mkdir(parents=True, exist_ok=True)
    clusters = _load_cluster_table(root / "emm_clusters.tsv")
    version = ""
    stamp = root / "emm.version"
    if stamp.exists():
        version = stamp.read_text(encoding="utf-8").split("\t")[0].strip()

    hsps = blast("blastn", assembly, fasta, work / "emm.blastn.tsv",
                 threads=threads, evalue=1e-20,
                 perc_identity=EMM_MIN_IDENTITY,
                 culling_limit=EMM_CULLING_LIMIT, max_target_seqs=10000)
    candidates = [_EmmHit(h) for h in hsps]
    passing = [c for c in candidates if c.passes]
    lengths = contig_lengths(assembly)

    if not passing:
        return _emm_undetermined(candidates, lengths, version)

    # One answer per locus: emm sits in a tandem array with mrp and enn, which
    # cross-match the emm database, so hits are grouped by position and only the
    # best in each group is an answer. Choosing between the groups is the step
    # that matters. In MGAS15252 (emm59) three loci match the database perfectly
    # or nearly so - emm174.1 at 1,564,601, emm59.0 at 1,566,000 and emm141.3 at
    # 1,567,472 - and emm174 is one of the types the reference laboratories
    # classify as *emm-like*: its sequence is in the database but it is carried at
    # the mrp or enn position. Ranking by score alone would report emm174.1 and
    # get the type wrong, so an emm-like locus never outranks a true emm locus.
    loci = _group_by_locus(passing)
    loci.sort(key=_locus_rank, reverse=True)
    best_group = loci[0]
    best = max(best_group, key=lambda c: (not _is_emm_like(c), c.score, c.bitscore))
    ties = sorted({c.subtype for c in best_group
                   if c.score >= best.score - 1e-9 and c.subtype != best.subtype})

    emm_type = best.emm_type
    cluster = clusters.get(best.allele.rsplit(".", 1)[0].upper(), "")
    pattern, tropism = _pattern_of(cluster)
    contig_len = lengths.get(best.contig, 0)
    at_edge = (best.start <= CONTIG_EDGE_BP
               or (contig_len and contig_len - best.end <= CONTIG_EDGE_BP))

    alleles = {
        "emm_type": emm_type,
        "emm_subtype": best.subtype if best.exact else f"{best.subtype}*",
        "emm_cluster": cluster or "not determined",
        "emm_pattern": pattern or "not determined",
        "tissue_tropism": tropism or "not determined",
        "evidence": (f"{best.contig}:{best.start}-{best.end} "
                     f"{best.identity:.1f}% identity over {best.align_len}/"
                     f"{best.ref_len} bases, {best.mismatch} mismatch(es), "
                     f"{best.gapopen} gap opening(s)"),
        "database": f"CDC alltrimmed.tfa ({version or 'version not recorded'})",
        "other_emm_region_loci": "; ".join(_describe_locus(g) for g in loci[1:])
                                 or "none",
    }
    if _SINGLE_PROTEIN.match(cluster or ""):
        alleles["emm_cluster"] = (
            f"{cluster} (a single-protein cluster; these are collected as "
            f"clade Y)")

    notes: list[str] = []
    if not best.exact:
        notes.append(
            f"no exact match: the closest database subtype {best.subtype} "
            f"differs by {best.mismatch} mismatch(es) and "
            f"{abs(best.ref_len - best.align_len)} base(s) of length, so the "
            f"emm type is {emm_type} but the subtype is provisional and may be "
            f"a new subtype - submit to the CDC for assignment")
    if ties:
        notes.append("equally good database subtype(s): " + ", ".join(ties))
    if best.allele.rsplit(".", 1)[0].upper() in EMM_LIKE:
        notes.append(
            f"{emm_type} is an emm-like type: the sequence is in the emm "
            f"database but the locus is usually mrp or enn, or the organism is a "
            f"group C/G streptococcus - confirm the species before reporting")
    if at_edge:
        notes.append(
            f"the emm region reaches within {CONTIG_EDGE_BP} bp of a contig end "
            f"({best.contig} is {contig_len} bp), so the hypervariable region "
            f"may be incomplete and the subtype unresolved; re-assemble or type "
            f"by the CDC PCR protocol to confirm")
    if not cluster:
        notes.append(
            f"no emm cluster could be read for {emm_type}: the CDC database "
            f"either carries no cluster annotation for it or its records "
            f"disagree")
    elif not pattern:
        notes.append(
            f"emm cluster {cluster} is a single-protein cluster, which does not "
            f"map onto an emm pattern, so the tissue-tropism grouping is not "
            f"determined")

    return TypingResult(
        scheme="emm",
        call=alleles["emm_subtype"],
        lineage=cluster or "-",
        alleles=alleles,
        loci_found=1,
        loci_total=1,
        note="; ".join(notes),
    )


def _emm_undetermined(candidates: list[_EmmHit], lengths: dict[str, int],
                      version: str) -> TypingResult:
    """No candidate passed the filters: say so, with the nearest sequences."""
    if not candidates:
        return TypingResult(
            scheme="emm", call="not determined", loci_total=1,
            alleles={"database": f"CDC alltrimmed.tfa ({version or 'unknown'})"},
            note="no part of the assembly matched the emm database at 95% "
                 "identity; the isolate may not be a pyogenic streptococcus, or "
                 "the emm region may be missing from the assembly")
    near = sorted(candidates, key=lambda c: c.score, reverse=True)[:4]
    best = near[0]
    contig_len = lengths.get(best.contig, 0)
    return TypingResult(
        scheme="emm", call="not determined", lineage="-", loci_found=0,
        loci_total=1,
        alleles={
            "candidates": "; ".join(
                f"{c.subtype} ({c.identity:.1f}% id, {c.align_len}/{c.ref_len} "
                f"bases, {c.mismatch} mismatch(es))" for c in near),
            "best_locus": f"{best.contig}:{best.start}-{best.end} "
                          f"(contig {contig_len} bp)",
            "database": f"CDC alltrimmed.tfa ({version or 'unknown'})",
        },
        note="the emm hypervariable region is present but no database subtype "
             "fits within the CDC calling rule (at most 4 mismatches, 2 gap "
             "openings and 5 bases of length difference). The candidates above "
             "are the nearest sequences. This is what a fragmented emm region "
             "looks like when the assembly breaks inside the hypervariable "
             "region, and it is also what a genuinely new emm subtype looks "
             "like: re-assemble, or type by the CDC PCR protocol, before "
             "reporting a type",
    )


def _describe_locus(group: list[_EmmHit]) -> str:
    """One secondary emm-region locus, named and placed."""
    best = max(group, key=lambda c: (c.score, c.bitscore))
    tag = " [emm-like: usually the mrp or enn position]" if _is_emm_like(best) else ""
    mark = "" if best.exact else "*"
    return f"{best.subtype}{mark}@{best.contig}:{best.start}-{best.end}{tag}"


def _is_emm_like(hit: _EmmHit) -> bool:
    return hit.allele.rsplit(".", 1)[0].upper() in EMM_LIKE


def _locus_rank(group: list[_EmmHit]) -> tuple:
    """Rank one locus: a true emm locus first, then exactness, then score."""
    true_emm = [c for c in group if not _is_emm_like(c)]
    pool = true_emm or group
    return (bool(true_emm), max(c.exact for c in pool),
            max(c.score for c in pool), max(c.bitscore for c in pool))


def _group_by_locus(hits: list[_EmmHit]) -> list[list[_EmmHit]]:
    """Group candidate hits that sit at the same place in the assembly.

    Single-linkage on position with a 500 bp reach, which is what the reference
    tool's hierarchical clustering achieves for hits of a 180 bp reference.
    """
    groups: list[list[_EmmHit]] = []
    for hit in sorted(hits, key=lambda h: (h.contig, h.start)):
        placed = False
        for group in groups:
            if group[0].contig != hit.contig:
                continue
            if any(abs(hit.start - other.start) <= EMM_CLUSTER_DISTANCE
                   or abs(hit.end - other.end) <= EMM_CLUSTER_DISTANCE
                   for other in group):
                group.append(hit)
                placed = True
                break
        if not placed:
            groups.append([hit])
    return groups


def _pattern_of(cluster: str) -> tuple[str, str]:
    """emm pattern and tissue-tropism grouping for an emm cluster, or blanks."""
    if not cluster:
        return "", ""
    for prefix, pattern, tropism in PATTERN_BY_CLUSTER:
        if cluster.upper().startswith(prefix):
            # A-D4 is an intermediate label and must not be read as pattern D.
            if cluster.upper().startswith("A-D"):
                return "", ""
            return pattern, tropism
    return "", ""


# ------------------------------------------------- superantigens and prophages
def _call_superantigens(hits: dict[str, GeneHit], manifest: dict) -> TypingResult:
    """The superantigen and DNase repertoire, split by genomic context."""
    unresolved = set(manifest.get("unresolved") or ())
    phage_sag = [s for s in SUPERANTIGENS if s.location == "prophage"]
    chrom_sag = [s for s in SUPERANTIGENS if s.location == "chromosome"]
    dnases = list(PROTEASES_AND_DNASES)

    alleles: dict[str, str] = {}
    present_phage = _present(hits, phage_sag, alleles, unresolved)
    present_chrom = _present(hits, chrom_sag, alleles, unresolved)
    present_dnase = _present(hits, dnases, alleles, unresolved)

    sags = present_phage + present_chrom
    call = "+".join(sags) if sags else "none detected"
    alleles["superantigen_profile"] = call
    alleles["phage_encoded_superantigens"] = (
        "+".join(present_phage) if present_phage else "none detected")
    alleles["chromosomal_superantigens"] = (
        "+".join(present_chrom) if present_chrom else "none detected")
    alleles["dnases_and_proteases"] = (
        "+".join(present_dnase) if present_dnase else "none detected")

    notes: list[str] = []
    for spec in phage_sag + chrom_sag + dnases:
        hit = hits.get(spec.name)
        if hit is None or not hit.present:
            continue
        if hit.state != "intact":
            notes.append(f"{spec.name} is not intact: {hit.lesion or hit.state}")
        elif hit.identity_pct < 100.0:
            # The allele numbers published for speA and smeZ are not available as
            # a versioned reference set, so the difference is reported instead of
            # being given a number it might not deserve.
            diffs = round((100.0 - hit.identity_pct) * hit.ref_len / 100.0)
            notes.append(
                f"{spec.name} differs from the reference allele "
                f"{hit.accession} at about {diffs} residue(s) "
                f"({hit.identity_pct:.1f}% aa identity): a variant allele, not "
                f"named here because no versioned allele set is published")
    missing = sorted(unresolved & {s.name for s in phage_sag + chrom_sag + dnases})
    if missing:
        notes.append(
            "not determined (no reference protein in the current RefSeq "
            "annotation of the panel strains): " + ", ".join(missing))
    if "speG" in alleles and alleles["speG"] == "not detected":
        notes.append(
            "speG was not detected; it is present in almost every group A "
            "streptococcus, so check the assembly before reporting its absence")

    total = len(phage_sag) + len(chrom_sag) + len(dnases) - len(missing)
    return TypingResult(
        scheme="emm_virulome",
        call=call,
        lineage=f"{len(present_phage)} phage-encoded / {len(present_chrom)} "
                f"chromosomal superantigen(s)",
        alleles=alleles,
        loci_found=len(sags) + len(present_dnase),
        loci_total=max(total, 0),
        note="; ".join(notes),
    )


def _present(hits: dict[str, GeneHit], specs: list[GeneSpec],
             alleles: dict[str, str], unresolved: set[str]) -> list[str]:
    """Record each spec's result in *alleles* and return the ones found."""
    found: list[str] = []
    for spec in specs:
        hit = hits.get(spec.name)
        if spec.name in unresolved:
            alleles[spec.name] = "not determined (no reference sequence)"
            continue
        if hit is None or not hit.present:
            alleles[spec.name] = "not detected"
            continue
        alleles[spec.name] = f"{hit.state} - {hit.evidence()}"
        found.append(spec.name)
    return found


def _call_prophage(hits: dict[str, GeneHit], every: dict[str, list[GeneHit]],
                   lengths: dict[str, int]) -> TypingResult:
    """Describe the prophage regions, and say where the assembly cannot.

    A prophage is called where at least three different phage module classes
    (lysogeny, packaging, head, tail, lysis) sit within reach of one another on
    one contig. Boundaries are not called: an assembly rarely resolves a
    prophage's attachment sites, so the region is reported as the span the
    marker genes cover, with the cargo genes that fall inside it.
    """
    # Every copy, not the best copy: each prophage has its own terminase and
    # portal, so one hit per marker gene would describe at most one prophage.
    markers = [h for spec in PHAGE_MARKERS
               for h in every.get(spec.name, []) if h.present]
    cargo = [h for spec in SUPERANTIGENS + PROTEASES_AND_DNASES
             if spec.location == "prophage"
             for h in every.get(spec.name, []) if h.present]

    regions = _link_regions(markers)
    alleles: dict[str, str] = {}
    described = 0
    orphan_cargo: list[GeneHit] = []
    claimed: set[str] = set()
    for index, region in enumerate(sorted(regions, key=lambda r: -len(r)), 1):
        classes = sorted({h.spec.name for h in region})
        if len(classes) < PHAGE_CLASSES_FOR_REGION:
            continue
        described += 1
        contig = region[0].contig
        lo = min(h.start for h in region)
        hi = max(h.end for h in region)
        inside = [h for h in cargo
                  if h.contig == contig
                  and lo - PHAGE_LINK_BP <= h.start <= hi + PHAGE_LINK_BP]
        claimed.update(h.spec.name for h in inside)
        contig_len = lengths.get(contig, 0)
        edge = (lo <= CONTIG_EDGE_BP
                or (contig_len and contig_len - hi <= CONTIG_EDGE_BP))
        alleles[f"prophage_{index}"] = (
            f"{contig}:{lo}-{hi} ({hi - lo + 1} bp spanned by "
            f"{len(classes)} phage module genes: {', '.join(classes)}); "
            f"cargo: {'+'.join(sorted(h.spec.name for h in inside)) or 'none in range'}"
            + ("; the region reaches a contig end, so its true extent and its "
               "attachment site are not resolved" if edge else ""))

    for hit in cargo:
        if hit.spec.name not in claimed:
            orphan_cargo.append(hit)
    if orphan_cargo:
        alleles["unplaced_phage_cargo"] = "; ".join(
            f"{h.spec.name} on {h.contig}:{h.start}-{h.end} "
            f"(contig {lengths.get(h.contig, 0)} bp)" for h in orphan_cargo)

    found_classes = sorted({h.spec.name for h in markers})
    alleles["phage_module_genes"] = ", ".join(found_classes) or "none detected"

    notes: list[str] = []
    if orphan_cargo:
        notes.append(
            f"{len(orphan_cargo)} prophage-encoded gene(s) could not be placed "
            f"in a prophage region: the contig they sit on carries too few phage "
            f"module genes. This is contig fragmentation, not evidence that the "
            f"gene is chromosomal - a draft assembly breaks inside prophages, "
            f"which are repetitive, far more often than elsewhere")
    if not described and markers:
        notes.append(
            f"phage module genes are present ({', '.join(found_classes)}) but no "
            f"contig carries {PHAGE_CLASSES_FOR_REGION} different module "
            f"classes within {PHAGE_LINK_BP // 1000} kb, so no prophage region "
            f"can be described from this assembly; the prophage content is not "
            f"resolved")
    if not markers:
        notes.append(
            "no phage module gene was detected. A group A streptococcus "
            "normally carries one to eight prophages, so this is more likely a "
            "very fragmented assembly, or a panel whose phage references are too "
            "distant, than a genuinely prophage-free isolate")
    notes.append(
        "prophage boundaries are not called: attachment sites are rarely "
        "resolved in a draft assembly, so each region is the span its marker "
        "genes cover")

    return TypingResult(
        scheme="spy_prophage",
        call=(f"{described} prophage region(s) described" if described
              else "no prophage region resolved"),
        lineage=("+".join(sorted(h.spec.name for h in cargo))
                 if cargo else "no phage-encoded cargo detected"),
        alleles=alleles,
        loci_found=len(found_classes),
        loci_total=len(PHAGE_MARKERS),
        note="; ".join(notes),
    )


def _link_regions(markers: list[GeneHit]) -> list[list[GeneHit]]:
    """Single-linkage grouping of phage marker hits along each contig."""
    groups: list[list[GeneHit]] = []
    for hit in sorted(markers, key=lambda h: (h.contig, h.start)):
        if groups and groups[-1][0].contig == hit.contig and \
                hit.start - max(h.end for h in groups[-1]) <= PHAGE_LINK_BP:
            groups[-1].append(hit)
        else:
            groups.append([hit])
    return groups


# ----------------------------------------------------------------- regulators
def _call_regulators(hits: dict[str, GeneHit]) -> TypingResult:
    """covR/covS/rocA and friends, reported as lesions rather than presence."""
    alleles: dict[str, str] = {}
    lesions: list[str] = []
    notes: list[str] = []
    for spec in REGULATORS + NGA_OPERON:
        hit = hits.get(spec.name)
        if hit is None or not hit.present:
            alleles[spec.name] = "not detected"
            notes.append(
                f"{spec.name} was not detected at all. It is a core gene, so "
                f"this means the locus is missing from the assembly, not that "
                f"the isolate lacks it - report as not determined")
            continue
        if hit.intact:
            alleles[spec.name] = f"intact ({hit.evidence()})"
            continue
        alleles[spec.name] = f"{hit.state} - {hit.evidence()}"
        if hit.state == "fragment":
            notes.append(
                f"{spec.name} is incomplete but its alignment reaches a contig "
                f"end, so the assembly and not the gene may be broken; the "
                f"lesion is not determined")
            continue
        wording = {
            "internal_stop": "a premature stop codon",
            "frameshift": "a frameshift",
            "truncated": "a truncation",
        }.get(hit.state, hit.state)
        # Only the covRS axis drives the summary. A lesion in ifs is common and
        # means something different (more free NADase), so it is reported in its
        # own field rather than as a covRS finding.
        if spec.name in ("covR", "covS", "rocA"):
            lesions.append(f"{spec.name} carries {wording}")
        else:
            notes.append(f"{spec.name} carries {wording}: {spec.meaning}")
        notes.append(f"{spec.name}: {hit.lesion}. {spec.meaning}")

    cov = [name for name in ("covR", "covS")
           if any(lesion.startswith(name) for lesion in lesions)]
    if cov:
        call = f"covRS lesion ({', '.join(cov)})"
        notes.insert(0,
                     "A loss-of-function lesion in covR or covS de-represses "
                     "capsule, SpyCEP, streptokinase and the nga/slo operon and "
                     "abolishes SpeB activity. This is the switch selected during "
                     "invasive infection, and such isolates are typically "
                     "recovered from sterile sites; it is not a transmissible "
                     "trait and it does not change antimicrobial choice")
    elif lesions:
        call = "; ".join(lesions)
    else:
        call = "no lesion detected in covR, covS or rocA"

    alleles["summary"] = call
    # The M1UK question, answered honestly rather than guessed at.
    alleles["M1UK_lineage"] = (
        "not determined: M1UK is defined by a set of SNPs relative to M1global "
        "and assigned phylogenetically (Lynskey 2019); this module reports the "
        "nga/ifs/slo genes and their state but does not implement the SNP set")
    return TypingResult(
        scheme="spy_regulators",
        call=call,
        lineage="-",
        alleles=alleles,
        loci_found=sum(1 for spec in REGULATORS + NGA_OPERON
                       if (h := hits.get(spec.name)) is not None and h.present),
        loci_total=len(REGULATORS) + len(NGA_OPERON),
        note="; ".join(notes),
    )


# ------------------------------------------------------------------- adhesins
def _call_adhesins(hits: dict[str, GeneHit]) -> TypingResult:
    alleles: dict[str, str] = {}
    found: list[str] = []
    notes: list[str] = []
    for spec in ADHESINS:
        hit = hits.get(spec.name)
        if hit is None or not hit.present:
            alleles[spec.name] = "not detected"
            continue
        alleles[spec.name] = f"{hit.state} - {hit.evidence()}"
        found.append(spec.name)
        if hit.state == "fragment":
            notes.append(f"{spec.name} reaches a contig end; presence is "
                         f"established but the allele is not")
    fibronectin = [n for n in ("prtF1", "fbaB", "fbp54") if n in found]
    collagen = [n for n in ("scl1", "scl2", "cpa") if n in found]
    alleles["fibronectin_binding"] = "+".join(fibronectin) or "none detected"
    alleles["collagen_binding"] = "+".join(collagen) or "none detected"
    alleles["fct_region_pilus"] = "+".join(
        n for n in ("cpa", "fctA", "fctB", "srtB") if n in found) or "none detected"
    notes.append(
        "the FCT region type (FCT-1 to FCT-9) and the T antigen type are not "
        "called: both need the pilus backbone allele to be matched against a "
        "published allele set, and no versioned, downloadable set exists")
    notes.append(
        "scl1, scl2, sic and emm carry variable-length repeats, so reference "
        "coverage below 100% is expected and is not a lesion")
    return TypingResult(
        scheme="spy_adhesins",
        call="+".join(found) if found else "none detected",
        lineage="-",
        alleles=alleles,
        loci_found=len(found),
        loci_total=len(ADHESINS),
        note="; ".join(notes),
    )
