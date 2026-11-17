"""Organism-specific typing for *Staphylococcus aureus*.

What this module reports, and why each call is here:

``spa``
    The repeat succession of the *spa* polymorphic X region, translated into the
    Ridom repeat codes and the spa type (t-number). spa typing is the single
    most widely used short-turnaround MRSA genotype in European reference
    laboratories, and the t-number is the unit outbreak reports are written in.
    The repeat region is a tandem array, which is exactly what a short-read
    assembler collapses, so this module refuses to call a type unless both
    flanking anchors are present on one contig and the succession tiles
    completely (see :func:`_spa`).

``agr``
    The *agr* group (I-IV). The group determines which autoinducing peptide the
    isolate makes and which groups it cross-inhibits; group is lineage-linked
    and is reported alongside the clone in most *S. aureus* virulence work.

``blaZ`` / ``bla_operon``
    The staphylococcal penicillinase and its regulators. ``blaZ`` alone does not
    settle penicillin susceptibility: expression needs the BlaR1 sensor to
    cleave the BlaI repressor, so a *blaZ*-positive isolate whose *blaR1* is
    truncated can test penicillin-susceptible. This module therefore states, per
    gene, whether it is intact, truncated, internally stopped, frameshifted, cut
    by a contig break, or absent - not merely present or absent.

``pvl``, ``tst``, ``enterotoxins``, ``haemolysins``, ``hlb``
    The toxin repertoire. ``hlb`` is reported as intact or as disrupted by
    prophage integration: the beta-haemolysin-converting phages of the Sa3int
    group integrate inside *hlb*, inactivating it while delivering the
    immune-evasion cluster, and "hlb present" on coverage alone is wrong for
    roughly nine in ten human isolates.

``capsule``
    Capsular serotype 5 versus 8, from the type-specific *capHIJ* genes.

``iec``
    The immune-evasion cluster type (A-G) carried by the same Sa3int prophage.

``acme``
    The arginine catabolic mobile element *arc* cluster, discriminated from the
    core-genome arginine deiminase operon that every *S. aureus* carries.

``clonal_complex``
    The clonal complex, from the 7-locus MLST profile, by the single-locus
    variant rule against the named founder sequence types.

Everything rests on BLAST+ and pure Python, so the same code runs inside the
frozen Windows build. SCCmec is deliberately absent: Hydra's existing
:mod:`hydra_amr.typing.cassette` types the cassette already.
"""

from __future__ import annotations

import json
import re
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

from ...engines.blast import (Hsp, blast, interval_length, merge_intervals,
                              overlap_length)
from ...records import TypingResult
from ...seqio import read_fasta, write_fasta
# Imported under another name: this module's own public entry point is run().
from ...utils import LOG, HydraError, revcomp, translate
from ...utils import run as run_command
from .base import OrganismDb, Reference

KEY = "saureus"
TITLE = "Staphylococcus aureus"

#: Species names, as :class:`hydra_amr.records.SpeciesCall` reports them, that
#: switch this module on. *S. argenteus* and *S. schweitzeri* are deliberately
#: excluded: they split out of the *S. aureus* complex, their MLST profiles are
#: in the same scheme, but none of the schemes below was defined in them.
SPECIES: tuple[str, ...] = ("Staphylococcus aureus",)

#: PubMLST scheme names that also switch it on.
MLST_SCHEMES: tuple[str, ...] = ("saureus",)


# --------------------------------------------------------------------- limits
#: Identity floor for a marker gene to count as present. Below this a hit is a
#: related family member, not the gene: *lukF-PV* aligns to the *hlgB*/*lukF*
#: family at 81% across every PVL-negative genome on the validation panel, and
#: *seb*, *sea* and *sep* cross-hit their relatives at 75-85%.
MIN_IDENTITY = 90.0
#: Coverage floor for "present". Between this and :data:`PARTIAL_COVERAGE` a
#: gene is reported as partial rather than silently present or absent.
MIN_COVERAGE = 90.0
PARTIAL_COVERAGE = 40.0

#: Separation on the contig, in bases, above which two alignment blocks of one
#: reference gene mean an insertion rather than an assembly artefact. The
#: Sa3int prophages insert 40-45 kb into *hlb*; a short indel is a few bases.
INSERT_GAP = 1000

#: Identity floor for assigning an *agr* group. Between-group identity over the
#: *agrBDC* window is 73-81% on the panel, within-group 98.5-100%.
AGR_MIN_IDENTITY = 92.0
AGR_MIN_COVERAGE = 80.0

#: Ridom spa primers. The repeat region is only typed when both are found on one
#: contig in the right relative orientation.
SPA_FWD = "TAAAGACGATCCTTCGGTGAGC"   # spa-1113f, J01786:1092-1113
SPA_REV = "CAGCAGTAGTGCCGTTTGCTT"    # spa-1514r, J01786:1534-1514
#: Width the extracted region must fall in, as the Ridom protocol's amplicon
#: plus room for the longest known succession (27 repeats).
SPA_MIN_REGION = 60
SPA_MAX_REGION = 1500
#: A repeat array ending this close to a contig end cannot be trusted to be
#: complete, whatever it tiles to.
SPA_EDGE_MARGIN = 50

USER_AGENT = "hydra-amr (+https://github.com/iowa69/hydra)"
TIMEOUT = 120

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
RIDOM_REPEATS = "https://spa.ridom.de/dynamic/sparepeats.fasta"
RIDOM_TYPES = "https://spa.ridom.de/dynamic/spatypes.txt"
PUBMLST_PROFILES = ("https://rest.pubmlst.org/db/pubmlst_saureus_seqdef/"
                    "schemes/1/profiles_csv")


# ---------------------------------------------------------------- marker panel
@dataclass(frozen=True)
class Marker:
    """One reference gene, pinned to an accession and a coordinate range.

    Pinning to coordinates rather than to a gene name matters: RefSeq calls
    *hla* ``hyl``, leaves *chp* unnamed, and names the capsule *H*/*I*/*J* genes
    ``cap8H``/``cap8I``/``cap8J`` in serotype-5 genomes as well as serotype-8
    ones, so a name-driven build would mistype the capsule.
    """

    name: str
    accession: str
    start: int
    stop: int
    strand: int          # 1 plus, 2 minus
    strain: str          # the strain the reference allele is taken from
    locus_tag: str
    role: str = ""       # which call it feeds
    product: str = ""


#: The reference panel. Every coordinate was read from the RefSeq annotation of
#: the named genome and round-tripped against NCBI efetch.
MARKERS: tuple[Marker, ...] = (
    # ---- Panton-Valentine leukocidin
    Marker("lukS-PV", "NC_007793.1", 1547149, 1548087, 2, "USA300", "SAUSA300_RS07545",
           "pvl", "Panton-Valentine leukocidin subunit S"),
    Marker("lukF-PV", "NC_007793.1", 1546170, 1547153, 2, "USA300", "SAUSA300_RS07540",
           "pvl", "Panton-Valentine leukocidin subunit F"),
    # ---- toxic shock syndrome toxin
    Marker("tst", "NC_002745.2", 2061076, 2061780, 1, "N315", "SA_RS10445",
           "tst", "toxic shock syndrome toxin TSST-1"),
    # ---- haemolysins and leukocidins
    Marker("hla", "NC_007793.1", 1156367, 1157326, 2, "USA300", "SAUSA300_RS05720",
           "haemolysin", "alpha-haemolysin (alpha-toxin)"),
    Marker("hlb", "NC_002951.2", 2063726, 2064718, 1, "COL", "SACOL_RS10470",
           "hlb", "beta-haemolysin (sphingomyelinase C)"),
    Marker("hld", "NC_007793.1", 2146572, 2146652, 2, "USA300", "SAUSA300_RS16125",
           "haemolysin", "delta-haemolysin"),
    Marker("hlgA", "NC_007793.1", 2542127, 2543056, 1, "USA300", "SAUSA300_RS13070",
           "haemolysin", "gamma-haemolysin subunit A"),
    Marker("hlgB", "NC_007793.1", 2544572, 2545549, 1, "USA300", "SAUSA300_RS13080",
           "haemolysin", "gamma-haemolysin subunit B"),
    Marker("hlgC", "NC_007793.1", 2543623, 2544570, 1, "USA300", "SAUSA300_RS13075",
           "haemolysin", "gamma-haemolysin subunit C"),
    Marker("lukE", "NC_007793.1", 1956304, 1957239, 2, "USA300", "SAUSA300_RS09685",
           "haemolysin", "leukocidin LukE"),
    Marker("lukD", "NC_007793.1", 1955319, 1956302, 2, "USA300", "SAUSA300_RS09680",
           "haemolysin", "leukocidin LukD"),
    Marker("lukA", "NC_007793.1", 2128768, 2129784, 2, "USA300", "SAUSA300_RS10845",
           "haemolysin", "leukocidin LukAB/LukGH subunit A"),
    Marker("lukB", "NC_007793.1", 2129806, 2130861, 2, "USA300", "SAUSA300_RS10850",
           "haemolysin", "leukocidin LukAB/LukGH subunit B"),
    # ---- enterotoxins and enterotoxin-like superantigens
    Marker("sea", "NC_009641.1", 2094197, 2094970, 2, "Newman", "NWMN_RS10835",
           "enterotoxin", "staphylococcal enterotoxin A"),
    Marker("seb", "NC_002951.2", 916478, 917278, 1, "COL", "SACOL_RS04655",
           "enterotoxin", "staphylococcal enterotoxin B"),
    Marker("sec", "NC_002745.2", 2058050, 2058850, 2, "N315", "SA_RS10435",
           "enterotoxin", "staphylococcal enterotoxin C"),
    Marker("seg", "NC_002745.2", 1876310, 1877086, 2, "N315", "SA_RS09240",
           "enterotoxin", "staphylococcal enterotoxin G (egc)"),
    Marker("seh", "NC_003923.1", 61524, 62249, 1, "MW2", "MW_RS00265",
           "enterotoxin", "staphylococcal enterotoxin H"),
    Marker("sei", "NC_002745.2", 1879088, 1879816, 2, "N315", "SA_RS09260",
           "enterotoxin", "staphylococcal enterotoxin I (egc)"),
    Marker("sek", "NC_007793.1", 883304, 884032, 2, "USA300", "SAUSA300_RS04320",
           "enterotoxin", "staphylococcal enterotoxin K"),
    Marker("sel", "NC_002745.2", 2057159, 2057881, 1, "N315", "SA_RS10430",
           "enterotoxin", "staphylococcal enterotoxin L"),
    Marker("sem", "NC_002745.2", 1879851, 1880570, 2, "N315", "SA_RS09265",
           "enterotoxin", "staphylococcal enterotoxin M (egc)"),
    Marker("sen", "NC_002745.2", 1877369, 1878124, 2, "N315", "SA_RS09245",
           "enterotoxin", "staphylococcal enterotoxin N (egc)"),
    Marker("seo", "NC_002745.2", 1880851, 1881633, 2, "N315", "SA_RS09270",
           "enterotoxin", "staphylococcal enterotoxin O (egc)"),
    Marker("sep", "NC_002745.2", 2011380, 2012162, 2, "N315", "SA_RS10125",
           "enterotoxin", "staphylococcal enterotoxin P"),
    Marker("seq", "NC_007793.1", 884056, 884784, 2, "USA300", "SAUSA300_RS04325",
           "enterotoxin", "staphylococcal enterotoxin Q"),
    Marker("seu", "NC_002745.2", 1878163, 1878934, 2, "N315", "SA_RS15170",
           "enterotoxin", "staphylococcal enterotoxin U (egc)"),
    Marker("selx", "NC_007793.1", 422037, 422648, 1, "USA300", "SAUSA300_RS01970",
           "enterotoxin", "staphylococcal superantigen-like toxin X"),
    # ---- immune evasion cluster
    Marker("scn", "NC_007793.1", 2085861, 2086211, 2, "USA300", "SAUSA300_RS10525",
           "iec", "staphylococcal complement inhibitor SCIN"),
    Marker("chp", "NC_007793.1", 2086896, 2087345, 1, "USA300", "SAUSA300_RS10530",
           "iec", "chemotaxis inhibitory protein CHIPS"),
    Marker("sak", "NC_007793.1", 2088425, 2088916, 2, "USA300", "SAUSA300_RS10540",
           "iec", "staphylokinase"),
    # ---- beta-lactamase operon
    Marker("blaZ", "NC_003140.1", 7945, 8790, 2, "N315", "SA_RS00040",
           "bla", "staphylococcal penicillinase BlaZ"),
    Marker("blaR1", "NC_003140.1", 8897, 10654, 1, "N315", "SA_RS00045",
           "bla", "beta-lactam sensor/transducer BlaR1"),
    Marker("blaI", "NC_003140.1", 10644, 11024, 1, "N315", "SA_RS00050",
           "bla", "penicillinase repressor BlaI"),
    # ---- capsule: the type-specific genes only
    Marker("cap5H", "NC_002745.2", 174660, 175286, 1, "N315", "SA_RS00925",
           "capsule", "capsular polysaccharide type 5 Cap5H"),
    Marker("cap5I", "NC_002745.2", 175291, 176400, 1, "N315", "SA_RS00930",
           "capsule", "capsular polysaccharide type 5 Cap5I"),
    Marker("cap5J", "NC_002745.2", 176414, 177580, 1, "N315", "SA_RS00935",
           "capsule", "capsular polysaccharide type 5 Cap5J"),
    Marker("cap8H", "NC_003923.1", 152869, 153948, 1, "MW2", "MW_RS00695",
           "capsule", "capsular polysaccharide type 8 Cap8H"),
    Marker("cap8I", "NC_003923.1", 153941, 155335, 1, "MW2", "MW_RS00700",
           "capsule", "capsular polysaccharide type 8 Cap8I"),
    Marker("cap8J", "NC_003923.1", 155332, 155889, 1, "MW2", "MW_RS00705",
           "capsule", "capsular polysaccharide type 8 Cap8J"),
    Marker("capA", "NC_002745.2", 167409, 168077, 1, "N315", "SA_RS00890",
           "capsule", "capsular polysaccharide type 5/8 CapA"),
    # ---- ACME, against the core-genome operon it is confused with
    Marker("arcA-ACME", "NC_007793.1", 73113, 74348, 2, "USA300", "SAUSA300_RS00330",
           "acme", "ACME arginine deiminase ArcA"),
    Marker("arcA-native", "NC_007793.1", 2781714, 2782949, 2, "USA300", "SAUSA300_RS14295",
           "acme", "core-genome arginine deiminase ArcA"),
    # ---- adhesins, biofilm, protein A, agr regulator
    Marker("icaA", "NC_007793.1", 2826611, 2827849, 1, "USA300", "SAUSA300_RS14455",
           "adhesin", "poly-N-acetylglucosamine synthase IcaA"),
    Marker("fnbA", "NC_007793.1", 2632745, 2635801, 2, "USA300", "SAUSA300_RS13530",
           "adhesin", "fibronectin-binding protein A"),
    Marker("fnbB", "NC_007793.1", 2629241, 2632063, 2, "USA300", "SAUSA300_RS13525",
           "adhesin", "fibronectin-binding protein B"),
    Marker("cna", "NC_003923.1", 2798376, 2801927, 2, "MW2", "MW_RS14170",
           "adhesin", "collagen adhesin Cna"),
    Marker("clfA", "NC_007793.1", 859148, 861949, 1, "USA300", "SAUSA300_RS04165",
           "adhesin", "clumping factor A"),
    Marker("clfB", "NC_007793.1", 2774340, 2777039, 2, "USA300", "SAUSA300_RS14270",
           "adhesin", "clumping factor B"),
    Marker("sdrC", "NC_007793.1", 611147, 613990, 1, "USA300", "SAUSA300_RS02915",
           "adhesin", "serine-aspartate repeat protein SdrC"),
    Marker("sdrD", "NC_007793.1", 614357, 618502, 1, "USA300", "SAUSA300_RS02920",
           "adhesin", "serine-aspartate repeat protein SdrD"),
    Marker("sdrE", "NC_007793.1", 618896, 622360, 1, "USA300", "SAUSA300_RS02925",
           "adhesin", "serine-aspartate repeat protein SdrE"),
    Marker("sbi", "NC_007793.1", 2540281, 2541591, 1, "USA300", "SAUSA300_RS13060",
           "adhesin", "immunoglobulin-binding protein Sbi"),
    Marker("spa", "NC_007793.1", 128167, 129693, 2, "USA300", "SAUSA300_RS00585",
           "spa", "staphylococcal protein A"),
    Marker("agrA", "NC_007793.1", 2148984, 2149700, 1, "USA300", "SAUSA300_RS10950",
           "agr", "quorum-sensing response regulator AgrA"),
    # ---- agr group windows: agrB start to agrC start + 650, the interval that
    #      carries agrB, agrD and the hypervariable 5' half of agrC.
    Marker("agr-I", "NC_007793.1", 2146942, 2148323, 1, "USA300", "agrBDC",
           "agr", "agr group I agrB-agrD-agrC' window"),
    Marker("agr-II", "NC_002745.2", 2079440, 2080835, 1, "N315", "agrBDC",
           "agr", "agr group II agrB-agrD-agrC' window"),
    Marker("agr-III", "NC_003923.1", 2108050, 2109431, 1, "MW2", "agrBDC",
           "agr", "agr group III agrB-agrD-agrC' window"),
    # Group IV has no complete reference genome with a published assignment, so
    # the window comes from the agr locus deposited with the paper that defined
    # the group: AF288215, the exfoliatin-producing strain RN5881, in which
    # agrB is 718..1281 and agrC starts at 1505.
    Marker("agr-IV", "AF288215.1", 718, 2155, 1, "RN5881", "agrBDC",
           "agr", "agr group IV agrB-agrD-agrC' window"),
)

#: Enterotoxin genes of the enterotoxin gene cluster, which travels as a unit on
#: *vSa*beta and is therefore reported as a cluster rather than gene by gene.
EGC_GENES = ("seg", "sei", "sem", "sen", "seo", "seu")

#: Founder sequence types of the clonal complexes in general use. A query ST is
#: assigned to a complex when it is that founder or a single-locus variant of
#: it. Connected-component clustering at the same 6/7 threshold is *not* used:
#: on the November 2026 PubMLST table it merges 4456 of 9503 STs - CC1, CC5,
#: CC8, CC15, CC80, CC97 and CC239 among them - into one component, which would
#: report ST8 as "CC5".
CC_FOUNDERS: tuple[int, ...] = (
    1, 5, 7, 8, 9, 12, 15, 20, 22, 25, 30, 45, 50, 59, 72, 78, 80, 88, 93, 97,
    101, 121, 126, 130, 133, 151, 152, 182, 188, 239, 361, 395, 398, 425, 479,
    522, 700, 705, 772, 1021,
)

MLST_LOCI = ("arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL")


# ------------------------------------------------------------------ provenance
DATABASES: tuple[OrganismDb, ...] = (
    OrganismDb(
        name="spa",
        source="https://spa.ridom.de/dynamic/sparepeats.fasta and "
               "https://spa.ridom.de/dynamic/spatypes.txt",
        version="retrieved 2026-11-17; upstream Last-Modified 2026-11-12 "
                "(865 repeat codes, 22886 spa types)",
        licence='Ridom spaServer Server Policy (https://spa.ridom.de/policy.shtml): '
                '"The nomenclature data (spa-repeat and -type assignments) of the '
                'SpaServer is the property of Ridom GmbH." '
                '"Individuals may download the database and use the content for '
                'non-commercial purposes only. Copies of the SpaServer database may '
                'not be used without permission to produce database driven web server '
                'or any other product or service." '
                'Attribution is required: "this website must be acknowledged in all '
                'publications".',
        redistribution="restricted",
        reference="Harmsen2003",
        notes="Must NOT be vendored into the Windows zip: the policy forbids using a "
              "copy of the database to produce a product or service. fetch() downloads "
              "it to the user's own database directory at install time, which is the "
              "permitted 'individuals may download' case. HTTPS only - the port-80 "
              "virtual host returns 404 for /dynamic/.",
    ),
    OrganismDb(
        name="markers",
        source="NCBI RefSeq, via efetch coordinate slices of NC_007793.1 "
               "(USA300_FPR3757), NC_002745.2 and NC_003140.1 (N315 and pN315), "
               "NC_003923.1 (MW2), NC_002951.2 (COL), NC_009641.1 (Newman)",
        version="accession-and-coordinate pinned; see MARKERS in this module",
        licence='NCBI: "Information that is created by or for the US government on '
                'this site is within the public domain." '
                '(https://www.ncbi.nlm.nih.gov/home/about/policies/)',
        redistribution="permitted",
        reference="Diep2006",
        notes="About 70 kB of sequence in total. Pinned to coordinates rather than "
              "gene names because RefSeq names hla 'hyl', leaves chp unnamed, and "
              "labels the capsule H/I/J genes cap8* in serotype 5 genomes too.",
    ),
    OrganismDb(
        name="mlst",
        source=PUBMLST_PROFILES + " and "
               "https://rest.pubmlst.org/db/pubmlst_saureus_seqdef/loci/<locus>/"
               "alleles_fasta",
        version="retrieved 2026-11-17; scheme 1 last_updated 2026-11-13, 12583 "
                "profile records (9503 reachable without authentication)",
        licence='PubMLST terms and conditions (https://pubmlst.org/terms-conditions, '
                'page updated 2026-05-13; the older /about/terms URL now 404s). The '
                'page names no open licence. It splits the data by submission date: '
                '"All data submitted on or before 31 December 2024 may continue to be '
                'downloaded, used, and redistributed without restriction, subject '
                'only to proper scientific citation and acknowledgement." For later '
                'data: "Downloading nomenclature data, including allele sequences, '
                'allele definitions, profiles, or sequence types, is permitted only '
                'for non-commercial academic research or public health surveillance '
                'use", and "Restricted access also applies to commercial services and '
                'software tools that enable the use or the downloading of data ... '
                'The publishers of such commercial tools will require a specific '
                'licence contract."',
        redistribution="permitted",
        reference="Enright2000",
        notes="Redistribution is permitted ONLY because fetch() downloads "
              "anonymously, and anonymous access returns just the pre-2025 data that "
              "the terms place under no restriction - 9503 of the scheme's 12583 "
              "profiles. Adding an API key to reach the rest would pull this database "
              "under the post-2024 restrictions, and a commercial Hydra would then "
              "need a licence contract from PubMLST. The visible cost of staying "
              "anonymous is that an ST assigned after 2024-12-31 reads as "
              "'not determined', with that reason given. The scheme carries no "
              "clonal_complex field, so the complex is derived here; see CC_FOUNDERS. "
              "If the integrator passes Hydra's existing MlstCall into "
              "clonal_complex() this database is not needed at all.",
    ),
)


REFERENCES: tuple[Reference, ...] = (
    Reference(
        key="Harmsen2003",
        authors="Harmsen D, Claus H, Witte W, Rothganger J, Claus H, Turnwald D, Vogel U.",
        year=2003,
        title="Typing of methicillin-resistant Staphylococcus aureus in a university "
              "hospital setting by using novel software for spa repeat determination "
              "and database management",
        journal="J Clin Microbiol 41:5442-5448",
        doi="10.1128/JCM.41.12.5442-5448.2003",
        pmid="14662923",
        supports="the spa repeat-code nomenclature and the t-number type assignments; "
                 "the citation the Ridom spaServer policy requires",
    ),
    Reference(
        key="Shopsin1999",
        authors="Shopsin B, Gomez M, Montgomery SO, et al.",
        year=1999,
        title="Evaluation of protein A gene polymorphic region DNA sequencing for "
              "typing of Staphylococcus aureus strains",
        journal="J Clin Microbiol 37:3556-3563",
        doi="10.1128/JCM.37.11.3556-3563.1999",
        pmid="10523551",
        supports="spa repeat-region sequencing as a typing method, and the repeat "
                 "structure the tiling in _spa() relies on",
    ),
    Reference(
        key="Frenay1996",
        authors="Frenay HM, Bunschoten AE, Schouls LM, et al.",
        year=1996,
        title="Molecular typing of methicillin-resistant Staphylococcus aureus on the "
              "basis of protein A gene polymorphism",
        journal="Eur J Clin Microbiol Infect Dis 15:60-64",
        doi="10.1007/BF01586186",
        pmid="8641305",
        supports="the original description of spa polymorphic-region typing, before "
                 "the Ridom repeat codes",
    ),
    Reference(
        key="Mollerup2022",
        authors="Mollerup S, Worning P, Petersen A, Bartels MD.",
        year=2022,
        title="spa Typing of Methicillin-Resistant Staphylococcus aureus Based on "
              "Whole-Genome Sequencing: the Impact of the Assembler",
        journal="Microbiol Spectr 10:e0218922",
        doi="10.1128/spectrum.02189-22",
        pmid="36350148",
        supports="the rule that a collapsed repeat array is typically two repeats "
                 "short and that the error concentrates in successions with repeated "
                 "identical repeats; the reason _spa() reports the raw succession and "
                 "refuses a call at a contig edge",
    ),
    Reference(
        key="Bartels2014",
        authors="Bartels MD, Petersen A, Worning P, et al.",
        year=2014,
        title="Comparing whole-genome sequencing with Sanger sequencing for spa typing "
              "of methicillin-resistant Staphylococcus aureus",
        journal="J Clin Microbiol 52:4305-4308",
        doi="10.1128/JCM.01979-14",
        pmid="25297335",
        supports="the measured rate of lost 24-bp repeats in short-read assemblies, "
                 "quantifying the caveat attached to every spa call",
    ),
    Reference(
        key="Bletz2015",
        authors="Bletz S, Mellmann A, Rothganger J, Harmsen D.",
        year=2015,
        title="Ensuring backwards compatibility: traditional genotyping efforts in the "
              "era of whole genome sequencing",
        journal="Clin Microbiol Infect 21:347.e1-347.e4",
        doi="10.1016/j.cmi.2014.11.005",
        pmid="25658529",
        supports="that spa typing accuracy from WGS is dominated by the assembly "
                 "pipeline, not the typing step",
    ),
    Reference(
        key="Enright2000",
        authors="Enright MC, Day NP, Davies CE, Peacock SJ, Spratt BG.",
        year=2000,
        title="Multilocus sequence typing for characterization of methicillin-resistant "
              "and methicillin-susceptible clones of Staphylococcus aureus",
        journal="J Clin Microbiol 38:1008-1015",
        doi="10.1128/JCM.38.3.1008-1015.2000",
        pmid="10698988",
        supports="the 7-locus S. aureus MLST scheme the clonal-complex call reads",
    ),
    Reference(
        key="Feil2004",
        authors="Feil EJ, Li BC, Aanensen DM, Hanage WP, Spratt BG.",
        year=2004,
        title="eBURST: inferring patterns of evolutionary descent among clusters of "
              "related bacterial genotypes from multilocus sequence typing data",
        journal="J Bacteriol 186:1518-1530",
        doi="10.1128/JB.186.5.1518-1530.2004",
        pmid="14973027",
        supports="the single-locus-variant definition of a clonal complex and the "
                 "founder-naming convention used by clonal_complex()",
    ),
    Reference(
        key="Diep2006",
        authors="Diep BA, Gill SR, Chang RF, et al.",
        year=2006,
        title="Complete genome sequence of USA300, an epidemic clone of "
              "community-acquired meticillin-resistant Staphylococcus aureus",
        journal="Lancet 367:731-739",
        doi="10.1016/S0140-6736(06)68231-7",
        pmid="16517273",
        supports="the USA300_FPR3757 reference genome most markers are taken from, "
                 "and the arginine catabolic mobile element the acme call detects",
    ),
    Reference(
        key="Kuroda2001",
        authors="Kuroda M, Ohta T, Uchiyama I, et al.",
        year=2001,
        title="Whole genome sequencing of meticillin-resistant Staphylococcus aureus",
        journal="Lancet 357:1225-1240",
        doi="10.1016/S0140-6736(00)04403-2",
        pmid="11418146",
        supports="the N315 reference genome and plasmid pN315, source of the tst, "
                 "egc, capsule type 5 and bla operon references",
    ),
    Reference(
        key="Baba2002",
        authors="Baba T, Takeuchi F, Kuroda M, et al.",
        year=2002,
        title="Genome and virulence determinants of high virulence community-acquired "
              "MRSA",
        journal="Lancet 359:1819-1827",
        doi="10.1016/S0140-6736(02)08713-5",
        pmid="12044378",
        supports="the MW2 reference genome, source of the PVL-positive, agr group III "
                 "and capsule type 8 references",
    ),
    Reference(
        key="Holden2004",
        authors="Holden MT, Feil EJ, Lindsay JA, et al.",
        year=2004,
        title="Complete genomes of two clinical Staphylococcus aureus strains: evidence "
              "for the rapid evolution of virulence and drug resistance",
        journal="Proc Natl Acad Sci U S A 101:9786-9791",
        doi="10.1073/pnas.0402521101",
        pmid="15213324",
        supports="the MRSA252 (EMRSA-16) and MSSA476 reference genomes used in "
                 "validation",
    ),
    Reference(
        key="Gill2005",
        authors="Gill SR, Fouts DE, Archer GL, et al.",
        year=2005,
        title="Insights on evolution of virulence and resistance from the complete "
              "genome analysis of an early methicillin-resistant Staphylococcus aureus "
              "strain and a biofilm-producing methicillin-resistant Staphylococcus "
              "epidermidis strain",
        journal="J Bacteriol 187:2426-2438",
        doi="10.1128/JB.187.7.2426-2438.2005",
        pmid="15774886",
        supports="the COL reference genome, the only panel member with an intact hlb, "
                 "from which the hlb reference allele is taken",
    ),
    Reference(
        key="vanWamel2006",
        authors="van Wamel WJ, Rooijakkers SH, Ruyken M, van Kessel KP, van Strijp JA.",
        year=2006,
        title="The innate immune modulators staphylococcal complement inhibitor and "
              "chemotaxis inhibitory protein of Staphylococcus aureus are located on "
              "beta-hemolysin-converting bacteriophages",
        journal="J Bacteriol 188:1310-1315",
        doi="10.1128/JB.188.4.1310-1315.2006",
        pmid="16452413",
        supports="the immune-evasion-cluster type A-G definitions, and that the same "
                 "prophage that carries them inactivates hlb",
    ),
    Reference(
        key="Ambler1991",
        authors="Ambler RP, Coulson AF, Frere JM, et al.",
        year=1991,
        title="A standard numbering scheme for the class A beta-lactamases",
        journal="Biochem J 276(Pt 1):269-270",
        doi="10.1042/bj2760269",
        pmid="2039479",
        supports="the numbering the blaZ allotype positions 128 and 216 are "
                 "expressed in",
    ),
    Reference(
        key="Voladri1996",
        authors="Voladri RK, Tummuru MK, Kernodle DS.",
        year=1996,
        title="Structure-function relationships among wild-type variants of "
              "Staphylococcus aureus beta-lactamase: importance of amino acids 128 "
              "and 216",
        journal="J Bacteriol 178:7248-7253",
        doi="10.1128/jb.178.24.7248-7253.1996",
        pmid="8955409",
        supports="the origin of the blaZ allotype scheme and of positions 128 and "
                 "216; only the abstract was obtainable, so the residue table comes "
                 "from Aarris2025",
    ),
    Reference(
        key="Carvajal2020",
        authors="Carvajal LP, Rincon S, Echeverri AM, et al.",
        year=2020,
        title="Novel insights into the classification of staphylococcal "
              "beta-lactamases in relation to the cefazolin inoculum effect",
        journal="Antimicrob Agents Chemother 64:e02511-19",
        doi="10.1128/AAC.02511-19",
        pmid="32071048",
        supports="that Ambler 128 and 216 are sequential residues 119 and 207 of the "
                 "281-residue preproenzyme, and the cefazolin inoculum effect "
                 "association reported with allotypes A and C",
    ),
    Reference(
        key="Aarris2025",
        authors="Aarris M, et al.",
        year=2025,
        title="Genetic variation in the blaZ gene leading to the BORSA phenotype in "
              "Staphylococcus aureus",
        journal="Antibiotics (Basel) 14:449",
        doi="10.3390/antibiotics14050449",
        pmid="40426516",
        supports="Table 1, the complete blaZ allotype A-F residue table implemented "
                 "in BLAZ_ALLOTYPE",
    ),
    Reference(
        key="Olsen2006",
        authors="Olsen JE, Christensen H, Aarestrup FM.",
        year=2006,
        title="Diversity and evolution of blaZ from Staphylococcus aureus and "
              "coagulase-negative staphylococci",
        journal="J Antimicrob Chemother 57:450-460",
        doi="10.1093/jac/dki492",
        pmid="16449305",
        supports="the diversity of blaZ (69 sequence types, 11 BlaZ protein types) "
                 "and its separate plasmid and chromosomal lines; it does NOT define "
                 "allotypes A-D, contrary to common citation",
    ),
    Reference(
        key="Jarraud2000",
        authors="Jarraud S, Lyon GJ, Figueiredo AM, et al.",
        year=2000,
        title="Exfoliatin-producing strains define a fourth agr specificity group in "
              "Staphylococcus aureus",
        journal="J Bacteriol 182:6517-6522",
        doi="10.1128/JB.182.22.6517-6522.2000",
        pmid="11053400",
        supports="the definition of agr group IV and GenBank AF288215 (strain "
                 "RN5881), the source of the group IV reference window",
    ),
    Reference(
        key="Jarraud2002",
        authors="Jarraud S, Mougel C, Thioulouse J, et al.",
        year=2002,
        title="Relationships between Staphylococcus aureus genetic background, "
              "virulence factors, agr groups (alleles), and human disease",
        journal="Infect Immun 70:631-641",
        doi="10.1128/IAI.70.2.631-641.2002",
        pmid="11796592",
        supports="the conventional agr typing region - the 3' end of agrB, all of "
                 "agrD and the 5' end of agrC - which is the window used here",
    ),
    Reference(
        key="Dufour2002",
        authors="Dufour P, Jarraud S, Vandenesch F, et al.",
        year=2002,
        title="High genetic variability of the agr locus in Staphylococcus species",
        journal="J Bacteriol 184:1180-1186",
        doi="10.1128/jb.184.4.1180-1186.2002",
        pmid="11807079",
        supports="that AgrB, AgrC and AgrD vary strikingly between the groups while "
                 "AgrA is conserved, which is why agrA is excluded from the group call",
    ),
    Reference(
        key="Gilot2002",
        authors="Gilot P, Lina G, Cochard T, Poutrel B.",
        year=2002,
        title="Analysis of the genetic variability of genes encoding the RNA "
              "III-activating components Agr and TRAP in a population of "
              "Staphylococcus aureus strains isolated from cows with mastitis",
        journal="J Clin Microbiol 40:4060-4067",
        doi="10.1128/JCM.40.11.4060-4067.2002",
        pmid="12409375",
        supports="the published agr group assignments for COL (I), N315 (II), "
                 "MRSA252 (III) and MSSA476 (III) used as validation expectations",
    ),
    Reference(
        key="Raghuram2022",
        authors="Raghuram V, Alexander AM, Loo HQ, Petit RA 3rd, Goldberg JB, Read TD.",
        year=2022,
        title="Species-wide phylogenomics of the Staphylococcus aureus Agr operon "
              "revealed convergent evolution of frameshift mutations",
        journal="Microbiol Spectr 10:e0133421",
        doi="10.1128/spectrum.01334-21",
        pmid="35315691",
        supports="AgrVATE, the reference tool whose motif database was used to "
                 "re-type USA300, Newman and RF122 as an independent cross-check",
    ),
    Reference(
        key="Tan2022",
        authors="Tan L, Li SR, Jiang B, Hu XM, Li S.",
        year=2022,
        title="Therapeutic targeting of the Staphylococcus aureus accessory gene "
              "regulator (agr) system",
        journal="Front Microbiol 13:700894",
        doi="10.3389/fmicb.2022.700894",
        pmid="35295312",
        supports="the published agr group assignments of the prototype strains "
                 "Newman (I), N315 (II) and MW2 (III), used as the group references",
    ),
    Reference(
        key="Baba2008",
        authors="Baba T, Bae T, Schneewind O, Takeuchi F, Hiramatsu K.",
        year=2008,
        title="Genome sequence of Staphylococcus aureus strain Newman and comparative "
              "analysis of staphylococcal genomes",
        journal="J Bacteriol 190:300-310",
        doi="10.1128/JB.01000-07",
        pmid="17951380",
        supports="Newman as the fully susceptible MSSA case: no SCCmec, "
                 "PVL-negative, no drug resistance genes",
    ),
    Reference(
        key="Holden2010",
        authors="Holden MT, Lindsay JA, Corton C, et al.",
        year=2010,
        title="Genome sequence of a recently emerged, highly transmissible, "
              "multi-antibiotic- and antiseptic-resistant variant of "
              "methicillin-resistant Staphylococcus aureus, sequence type 239 (TW)",
        journal="J Bacteriol 192:888-892",
        doi="10.1128/JB.01255-09",
        pmid="19948800",
        supports="TW20 (ST239), the CC8/CC30 hybrid that exercises the "
                 "clonal-complex rule's ambiguity handling",
    ),
    Reference(
        key="HerronOlson2007",
        authors="Herron-Olson L, Fitzgerald JR, Musser JM, Kapur V.",
        year=2007,
        title="Molecular correlates of host specialization in Staphylococcus aureus",
        journal="PLoS One 2:e1120",
        doi="10.1371/journal.pone.0001120",
        pmid="17971880",
        supports="RF122 as the awkward validation case: spa is a pseudogene, agrB and "
                 "agrD match group II while the agrC 5' end differs from all four "
                 "groups, tst-positive via SaPIbov1, and lukM/lukF' rather than PVL",
    ),
    Reference(
        key="Rouha2018",
        authors="Rouha H, Weber S, Janesch P, et al.",
        year=2018,
        title="Disarming Staphylococcus aureus from destroying human cells by "
              "simultaneously neutralizing six cytotoxins with two human monoclonal "
              "antibodies",
        journal="Virulence 9:231-247",
        doi="10.1080/21505594.2017.1391447",
        pmid="29099326",
        supports="Table 1: the published spa types and lukSF-PV status for MRSA252, "
                 "MW2 and Newman. Its MW2 row (t125) is the disagreement investigated "
                 "in docs/typing/saureus.md",
    ),
    Reference(
        key="Beckley2026",
        authors="Beckley B, Amarh V.",
        year=2026,
        title="StaphScope: a genomic typing tool for Staphylococcus aureus",
        journal="BMC Genomics 27",
        doi="10.1186/s12864-026-12609-x",
        pmid="41645058",
        supports="Table 1: the published spa and ST expectations for USA300 (t008), "
                 "N315 (t002), MRSA252 (t018), TW20 (t037) and NCTC 8325 (t211)",
    ),
)


# ------------------------------------------------------------------- downloads
def _get(url: str, dest: Path, attempts: int = 4) -> Path:
    """Download *url* to *dest*, retrying transport failures."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    last: Exception | None = None
    for attempt in range(attempts):
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
                payload = response.read()
            if not payload:
                raise HydraError(f"{url} returned an empty file")
            dest.write_bytes(payload)
            return dest
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code not in (429, 500, 502, 503, 504):
                raise HydraError(f"could not download {url}: {exc}") from exc
        except (urllib.error.URLError, OSError) as exc:
            last = exc
        if attempt + 1 < attempts:
            time.sleep(1.5 * (attempt + 1))
    raise HydraError(f"could not download {url}: {last}")


def _efetch_slice(marker: Marker, attempts: int = 5) -> str:
    """One reference gene, as the accession's own bases over its own range.

    NCBI answers a burst of unauthenticated requests with HTTP 429, so this
    backs off and retries rather than failing the whole build for one gene.
    """
    url = (f"{NCBI_EFETCH}?db=nuccore&id={marker.accession}&rettype=fasta"
           f"&retmode=text&seq_start={marker.start}&seq_stop={marker.stop}"
           f"&strand={marker.strand}")
    text = ""
    last: Exception | None = None
    for attempt in range(attempts):
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
                text = response.read().decode("utf-8", "replace")
            break
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code not in (429, 500, 502, 503, 504):
                raise HydraError(f"{marker.name}: could not fetch "
                                 f"{marker.accession}: {exc}") from exc
        except (urllib.error.URLError, OSError) as exc:
            last = exc
        time.sleep(2.0 * (attempt + 1))
    else:
        raise HydraError(f"{marker.name}: could not fetch {marker.accession} "
                         f"after {attempts} attempts: {last}")
    seq = "".join(line.strip() for line in text.splitlines()
                  if line and not line.startswith(">"))
    expected = marker.stop - marker.start + 1
    if len(seq) != expected:
        raise HydraError(f"{marker.name}: {marker.accession}:{marker.start}-"
                         f"{marker.stop} returned {len(seq)} bases, expected {expected}")
    return seq.upper()


def fetch(db_dir: Path) -> None:
    """Download and build this module's data into ``db_dir/saureus/``.

    Idempotent: anything already present and complete is left alone, so the
    call is cheap to repeat and safe to run on a partially built directory.
    """
    root = Path(db_dir) / KEY
    root.mkdir(parents=True, exist_ok=True)

    # ---- Ridom spa data. Downloaded, never vendored; see DATABASES.
    spa_dir = root / "spa"
    repeats = spa_dir / "sparepeats.fasta"
    types = spa_dir / "spatypes.txt"
    if not (repeats.exists() and repeats.stat().st_size > 10_000):
        LOG.info("saureus: downloading the Ridom spa repeat codes")
        _get(RIDOM_REPEATS, repeats)
    if not (types.exists() and types.stat().st_size > 100_000):
        LOG.info("saureus: downloading the Ridom spa type table")
        _get(RIDOM_TYPES, types)

    # ---- marker panel, one coordinate slice per gene
    marker_dir = root / "markers"
    panel = marker_dir / "markers.fna"
    have: dict[str, str] = {}
    if panel.exists():
        have = {header.split("|")[0]: seq for header, seq in read_fasta(panel)}
    if set(have) != {m.name for m in MARKERS}:
        LOG.info("saureus: building the %d-marker reference panel from NCBI",
                 len(MARKERS))
        records: list[tuple[str, str]] = []
        for marker in MARKERS:
            seq = have.get(marker.name)
            if not seq:
                seq = _efetch_slice(marker)
                time.sleep(0.4)      # NCBI asks for at most three requests a second
            records.append((_marker_header(marker), seq))
        marker_dir.mkdir(parents=True, exist_ok=True)
        write_fasta(panel, records)
    _makeblastdb(panel)

    # ---- PubMLST, for the clonal-complex call only
    mlst_dir = root / "mlst"
    profiles = mlst_dir / "profiles.tsv"
    if not (profiles.exists() and profiles.stat().st_size > 10_000):
        LOG.info("saureus: downloading the PubMLST saureus profile table")
        _get(PUBMLST_PROFILES, profiles)
    alleles = mlst_dir / "alleles.fna"
    if not (alleles.exists() and alleles.stat().st_size > 1_000_000):
        LOG.info("saureus: downloading the 7 PubMLST saureus allele sets")
        records = []
        for locus in MLST_LOCI:
            raw = mlst_dir / f"{locus}.tfa"
            _get("https://rest.pubmlst.org/db/pubmlst_saureus_seqdef/loci/"
                 f"{locus}/alleles_fasta", raw)
            for header, seq in read_fasta(raw):
                name = header.split()[0]
                allele = name.rsplit("_", 1)[-1]
                records.append((f"{locus}_{allele}", seq))
            raw.unlink()
        write_fasta(alleles, records)
    _makeblastdb(alleles)

    (root / "manifest.json").write_text(json.dumps({
        "key": KEY,
        "built": time.strftime("%Y-%m-%d"),
        "databases": [db.name for db in DATABASES],
        "markers": len(MARKERS),
        "sources": {db.name: db.source for db in DATABASES},
        "versions": {db.name: db.version for db in DATABASES},
        "licences": {db.name: db.licence for db in DATABASES},
        "redistribution": {db.name: db.redistribution for db in DATABASES},
    }, indent=2) + "\n")


def _marker_header(marker: Marker) -> str:
    return (f"{marker.name}|{marker.accession}:{marker.start}-{marker.stop}"
            f"|{marker.strain}|{marker.locus_tag}|{marker.role}|{marker.product}")


def _makeblastdb(fasta: Path) -> None:
    if not fasta.exists():
        raise HydraError(f"saureus: {fasta} was not built")
    if Path(f"{fasta}.nin").exists():
        return
    run_command(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl",
                 "-out", str(fasta)])


# ------------------------------------------------------------------ hit model
@dataclass
class MarkerHit:
    """One reference marker, as found in one assembly.

    Coverage and identity are measured on the *reference* axis, which is what
    decides whether a gene is there. The contig blocks are kept separately
    because a reference that aligns in two widely separated blocks is the
    signature of something inserted into the gene, not of a poor match.
    """

    name: str
    reference: str                      # accession:start-stop
    strain: str = ""
    locus_tag: str = ""
    role: str = ""
    product: str = ""
    reference_length: int = 0
    contig: str = ""
    contig_length: int = 0
    strand: str = "+"
    identity_pct: float = 0.0
    coverage_pct: float = 0.0
    reference_blocks: list[tuple[int, int]] = field(default_factory=list)
    contig_blocks: list[tuple[int, int]] = field(default_factory=list)
    #: Largest jump on the contig between consecutive reference blocks.
    insert_gap: int = 0
    #: Where that jump sits on the reference, as a 1-based base position.
    insert_at: int = 0
    #: How many alignments made up the hit. Not the same as the number of
    #: merged reference blocks: a prophage insertion leaves two alignments that
    #: overlap by the duplicated attachment site, which merge into one block.
    n_alignments: int = 1

    @property
    def present(self) -> bool:
        return (self.coverage_pct >= MIN_COVERAGE
                and self.identity_pct >= MIN_IDENTITY)

    @property
    def partial(self) -> bool:
        return (not self.present and self.identity_pct >= MIN_IDENTITY
                and self.coverage_pct >= PARTIAL_COVERAGE)

    @property
    def evidence(self) -> str:
        spans = ",".join(f"{lo}-{hi}" for lo, hi in self.reference_blocks)
        return (f"{self.identity_pct:.1f}% identity over {spans}/"
                f"{self.reference_length} of {self.reference} "
                f"on {self.contig}")

    def at_contig_edge(self, margin: int = SPA_EDGE_MARGIN) -> bool:
        if not self.contig_blocks or not self.contig_length:
            return False
        lo = min(lo for lo, _ in self.contig_blocks)
        hi = max(hi for _, hi in self.contig_blocks)
        return lo <= margin or hi >= self.contig_length - margin


def _contig_lengths(assembly: Path) -> dict[str, int]:
    return {header.split()[0]: len(seq) for header, seq in read_fasta(assembly)}


def _merge_on_reference(hsps: list[Hsp], lengths: dict[str, int]) -> list[MarkerHit]:
    """Group HSPs by (marker, contig, strand) and merge on the reference axis.

    Hydra's :func:`~hydra_amr.engines.blast.merge_hsps` merges on the subject
    axis, which is the right axis when the reference is the subject. Here the
    markers are the *query* - that is what puts the hit coordinates on the
    assembly's contigs, which the lesion calls need - so the merge is done here
    over the same interval helpers.
    """
    grouped: dict[tuple[str, str, str], list[Hsp]] = {}
    for hsp in hsps:
        grouped.setdefault((hsp.qseqid, hsp.sseqid, hsp.strand), []).append(hsp)

    hits: list[MarkerHit] = []
    for (qseqid, sseqid, strand), group in grouped.items():
        fields = qseqid.split("|")
        name = fields[0]
        # Keep the strongest HSPs, skipping any that re-covers reference ground
        # already claimed: two copies of one reference must not add up to more
        # than 100% coverage.
        group.sort(key=lambda h: h.bitscore, reverse=True)
        kept: list[Hsp] = []
        claimed: list[tuple[int, int]] = []
        for hsp in group:
            span = (hsp.q_lo, hsp.q_hi)
            length = span[1] - span[0] + 1
            if any(overlap_length(span, prev) > 0.5 * length for prev in claimed):
                continue
            kept.append(hsp)
            claimed.append(span)
        if not kept:
            continue
        qlen = kept[0].qlen
        ref_blocks = merge_intervals((h.q_lo, h.q_hi) for h in kept)
        covered = interval_length(ref_blocks)
        # Weight identity by the reference each HSP newly covers, so partly
        # overlapping HSPs do not count shared bases twice.
        seen: list[tuple[int, int]] = []
        ident = span_len = 0.0
        for hsp in kept:
            span = (hsp.q_lo, hsp.q_hi)
            width = span[1] - span[0] + 1
            fresh = max(0, width - sum(overlap_length(span, p) for p in seen))
            if fresh and width:
                share = fresh / width
                ident += hsp.nident * share
                span_len += hsp.length * share
            seen.append(span)

        # The insertion test: walk the HSPs in reference order and measure how
        # far the contig jumps between them.
        ordered = sorted(kept, key=lambda h: h.q_lo)
        gap = at = 0
        for left, right in zip(ordered, ordered[1:]):
            jump = abs(right.sub_lo - left.sub_hi)
            if jump > gap:
                gap, at = jump, left.q_hi

        hits.append(MarkerHit(
            name=name,
            reference=fields[1] if len(fields) > 1 else "",
            strain=fields[2] if len(fields) > 2 else "",
            locus_tag=fields[3] if len(fields) > 3 else "",
            role=fields[4] if len(fields) > 4 else "",
            product=fields[5] if len(fields) > 5 else "",
            reference_length=qlen,
            contig=sseqid,
            contig_length=lengths.get(sseqid, 0),
            strand=strand,
            identity_pct=(100.0 * ident / span_len) if span_len else 0.0,
            coverage_pct=(100.0 * covered / qlen) if qlen else 0.0,
            reference_blocks=ref_blocks,
            contig_blocks=merge_intervals((h.sub_lo, h.sub_hi) for h in kept),
            insert_gap=gap,
            insert_at=at,
            n_alignments=len(kept),
        ))
    return hits


def _best_per_marker(hits: list[MarkerHit]) -> dict[str, MarkerHit]:
    """One hit per marker: the one covering most of the reference, most closely."""
    best: dict[str, MarkerHit] = {}
    for hit in hits:
        score = hit.coverage_pct * hit.identity_pct
        current = best.get(hit.name)
        if current is None or score > current.coverage_pct * current.identity_pct:
            best[hit.name] = hit
    return best


# ---------------------------------------------------------------- spa typing
def _load_repeats(path: Path) -> dict[str, str]:
    """Repeat sequence -> Ridom code, with the ``r`` prefix the type table omits."""
    out: dict[str, str] = {}
    for header, seq in read_fasta(path):
        code = header.strip().lstrip("r")
        if code and seq:
            out[seq.strip().upper()] = code
    return out


def _load_spa_types(path: Path) -> dict[str, str]:
    """Repeat succession -> spa type.

    The ``NT,`` sentinel at the end of the Ridom table has an empty succession;
    keeping it would make every failed extraction look up as type ``NT``.
    """
    out: dict[str, str] = {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            name, _, profile = line.strip().partition(",")
            profile = profile.strip()
            if name and profile:
                out[profile] = name.strip()
    return out


def _tile_repeats(region: str, repeats: dict[str, str],
                  lengths: tuple[int, ...]) -> tuple[list[str], int, int, int]:
    """Walk *region* left to right, taking the longest repeat that matches.

    Longest-match-first matters because the Ridom repeats are 18, 21, 24, 27 and
    30 bases long and the shorter ones are prefixes of longer ones; taking
    whichever length happened to be tried first would mis-tile the array.

    Returns the codes, the number of unrecognised gaps between recognised
    repeats, and the unrecognised lead and trail around the array.
    """
    codes: list[str] = []
    index = 0
    first: int | None = None
    end = 0
    gaps = 0
    while index < len(region):
        taken = 0
        for width in lengths:
            code = repeats.get(region[index:index + width])
            if code is not None:
                taken = width
                if first is None:
                    first = index
                elif index != end:
                    gaps += 1
                codes.append(code)
                break
        index += taken if taken else 1
        if taken:
            end = index
    if first is None:
        return [], 0, len(region), 0
    return codes, gaps, first, len(region) - end


def _spa_region(assembly: Path, workdir: Path, db: Path,
                threads: int) -> tuple[str, str, int, int, int, str]:
    """Excise the repeat region between the two Ridom primers.

    Returns (region, contig, lo, hi, contig_length, problem). A problem string
    means no region could be excised, and says why.
    """
    query = workdir / "spa_primers.fna"
    write_fasta(query, [("spa-1113f", SPA_FWD), ("spa-1514r", SPA_REV)])
    hsps = blast("blastn", query, db, workdir / "spa_primers.tsv",
                 threads=threads, evalue=1.0, task="blastn-short", word_size=11,
                 max_target_seqs=1000)
    # A primer is an anchor only if nearly all of it matched.
    anchors: dict[str, list[Hsp]] = {"spa-1113f": [], "spa-1514r": []}
    for hsp in hsps:
        if hsp.qseqid in anchors and hsp.length >= hsp.qlen - 2 and hsp.pident >= 90.0:
            anchors[hsp.qseqid].append(hsp)
    if not anchors["spa-1113f"] and not anchors["spa-1514r"]:
        return "", "", 0, 0, 0, ("neither spa flanking primer site was found; the spa "
                                 "gene is absent or too diverged to anchor")
    if not anchors["spa-1113f"] or not anchors["spa-1514r"]:
        found = "spa-1113f" if anchors["spa-1113f"] else "spa-1514r"
        return "", "", 0, 0, 0, (f"only the {found} anchor was found, so the repeat "
                                 f"region is incomplete - most likely split by a "
                                 f"contig break")

    contigs = {header.split()[0]: seq for header, seq in read_fasta(assembly)}
    candidates: list[tuple[int, str, str, int, int, int]] = []
    for fwd in anchors["spa-1113f"]:
        for rev in anchors["spa-1514r"]:
            if fwd.sseqid != rev.sseqid:
                continue
            sequence = contigs.get(fwd.sseqid, "")
            if not sequence:
                continue
            if fwd.strand == "+" and rev.strand == "-":
                lo, hi = fwd.sub_lo, rev.sub_hi
                region = sequence[lo - 1:hi].upper()
            elif fwd.strand == "-" and rev.strand == "+":
                lo, hi = rev.sub_lo, fwd.sub_hi
                region = revcomp(sequence[lo - 1:hi].upper())
            else:
                continue
            if not (SPA_MIN_REGION <= len(region) <= SPA_MAX_REGION):
                continue
            candidates.append((len(region), region, fwd.sseqid, lo, hi, len(sequence)))
    if not candidates:
        return "", "", 0, 0, 0, ("both spa anchors were found but never together on "
                                 "one contig in the right orientation, so the repeat "
                                 "region is split by a contig break")
    if len({c[2] for c in candidates}) > 1 or len(candidates) > 1:
        LOG.debug("saureus: %d candidate spa regions; taking the longest",
                  len(candidates))
    candidates.sort(reverse=True)
    _, region, contig, lo, hi, contig_length = candidates[0]
    return region, contig, lo, hi, contig_length, ""


#: Successions this short are rare - the Ridom table averages 9.5 repeats - and
#: are the signature of a deleted or degraded array rather than a short type.
SPA_SHORT_ARRAY = 3


def _spa_gene_state(best: dict[str, MarkerHit], contigs: dict[str, str]) -> str:
    """What is wrong with the *spa* gene itself, if anything.

    A repeat succession can tile cleanly to a real t-number while the gene it
    came from encodes nothing: *spa* is a pseudogene in the bovine ST151
    lineage. Reporting the t-number without that is the kind of bare call this
    module exists to avoid.
    """
    lesion = _describe_lesion(best.get("spa"), "spa", contigs, 508)
    if lesion.state in ("intact", "absent", "not determined"):
        return ""
    return lesion.detail


def _spa(assembly: Path, root: Path, workdir: Path, db: Path,
         threads: int, best: dict[str, MarkerHit] | None = None,
         contigs: dict[str, str] | None = None) -> TypingResult:
    repeats_file = root / "spa" / "sparepeats.fasta"
    types_file = root / "spa" / "spatypes.txt"
    if not repeats_file.exists() or not types_file.exists():
        return TypingResult(scheme="spa", call="not determined",
                            note="the Ridom spa repeat and type files are not "
                                 "installed; run the saureus fetch step")
    repeats = _load_repeats(repeats_file)
    types = _load_spa_types(types_file)
    lengths = tuple(sorted({len(s) for s in repeats}, reverse=True))

    region, contig, lo, hi, contig_length, problem = _spa_region(
        assembly, workdir, db, threads)
    if problem:
        return TypingResult(scheme="spa", call="not determined",
                            loci_total=1,
                            note=f"{problem}. spa typing needs the whole tandem "
                                 f"array; reads or a long-read assembly would "
                                 f"resolve it")

    codes, gaps, lead, trail = _tile_repeats(region, repeats, lengths)
    succession = "-".join(codes)
    alleles = {
        "repeat_succession": succession or "-",
        "repeats": str(len(codes)),
        "region": f"{contig}:{lo}-{hi} ({len(region)} bp)",
        "unrecognised_lead_bp": str(lead),
        "unrecognised_trail_bp": str(trail),
    }
    if not codes:
        return TypingResult(scheme="spa", call="not determined", loci_total=1,
                            alleles=alleles,
                            note=f"the {len(region)} bp between the spa anchors "
                                 f"contained no recognised Ridom repeat; the array "
                                 f"may be deleted or the region misassembled")

    # A succession running into a contig end cannot be known to be complete.
    at_edge = lo <= SPA_EDGE_MARGIN or hi >= contig_length - SPA_EDGE_MARGIN
    caveats: list[str] = []
    if gaps:
        caveats.append(f"{gaps} unrecognised stretch(es) inside the array")
    if lead > 80 or trail > 100:
        caveats.append(f"{lead} bp before and {trail} bp after the array were not "
                       f"recognised repeats")
    if at_edge:
        caveats.append(f"the region reaches within {SPA_EDGE_MARGIN} bp of the end of "
                       f"{contig} ({contig_length} bp), so repeats may be missing")
    if len(codes) < SPA_SHORT_ARRAY:
        caveats.append(f"only {len(codes)} repeat(s): successions this short are rare "
                       f"(the Ridom table averages 9.5) and usually mean the array is "
                       f"deleted or degraded rather than genuinely short")
    gene_problem = _spa_gene_state(best or {}, contigs or {})
    if gene_problem:
        caveats.append(f"the spa gene itself is not intact - {gene_problem} - so "
                       f"protein A is not produced and this t-number describes a "
                       f"non-functional gene")

    call = types.get(succession, "")
    if not call:
        return TypingResult(
            scheme="spa", call="new succession", lineage=succession,
            loci_found=1, loci_total=1, alleles=alleles,
            note="; ".join([
                f"{len(codes)} repeats, succession {succession}, which is not in the "
                f"Ridom spaServer table - a new type, or a misassembled array",
                "submit it to spaServer for a t-number",
                *caveats]),
        )
    if gaps or at_edge:
        return TypingResult(
            scheme="spa", call="not determined", lineage=succession,
            loci_found=1, loci_total=1, alleles=alleles,
            note="; ".join([
                f"the succession {succession} matches {call}, but it is not reported "
                f"as a type because the evidence is incomplete",
                *caveats,
                "a collapsed array is typically two repeats short (Mollerup 2022), so "
                "a shortened succession can still match a real type",
            ]),
        )
    return TypingResult(
        scheme="spa", call=call, lineage=succession,
        loci_found=1, loci_total=1, alleles=alleles,
        note="; ".join([
            f"{len(codes)} Ridom repeats, succession {succession}, read from "
            f"{contig}:{lo}-{hi}",
            *caveats,
            "spa types are tandem arrays and short-read assemblies lose repeats in "
            "about 1-15% of isolates depending on the assembler (Mollerup 2022), so "
            "confirm before acting on a single-repeat difference",
        ]),
    )


# ---------------------------------------------------------------- agr grouping
def _agr(best: dict[str, MarkerHit]) -> TypingResult:
    windows = {name: hit for name, hit in best.items() if name.startswith("agr-")}
    ranked = sorted(
        ((hit.coverage_pct * hit.identity_pct / 100.0, name, hit)
         for name, hit in windows.items()), reverse=True)
    agr_a = best.get("agrA")
    alleles = {name: f"{hit.identity_pct:.1f}% identity over "
                     f"{hit.coverage_pct:.0f}% of the window"
               for _, name, hit in ranked}
    if agr_a is not None:
        alleles["agrA"] = (f"{agr_a.identity_pct:.1f}% identity over "
                           f"{agr_a.coverage_pct:.0f}%")

    if not ranked:
        return TypingResult(scheme="agr", call="not determined", loci_total=1,
                            alleles=alleles,
                            note="no agrB-agrD-agrC window aligned; the agr locus is "
                                 "absent, deleted or on a contig edge")
    score, name, hit = ranked[0]
    group = name.split("-", 1)[1]
    if hit.identity_pct < AGR_MIN_IDENTITY or hit.coverage_pct < AGR_MIN_COVERAGE:
        runner = f", next {ranked[1][1]} at {ranked[1][2].identity_pct:.1f}%" \
            if len(ranked) > 1 else ""
        return TypingResult(
            scheme="agr", call="not determined", loci_total=1, alleles=alleles,
            note=f"the closest reference window is {name} at "
                 f"{hit.identity_pct:.1f}% identity over {hit.coverage_pct:.0f}% of "
                 f"it{runner}, below the {AGR_MIN_IDENTITY:.0f}% / "
                 f"{AGR_MIN_COVERAGE:.0f}% floor. Group IV is not in this release's "
                 f"reference set, so a group IV isolate lands here")
    margin = ""
    if len(ranked) > 1:
        margin = (f"; next closest {ranked[1][1]} at "
                  f"{ranked[1][2].identity_pct:.1f}% identity over "
                  f"{ranked[1][2].coverage_pct:.0f}%")
    caveat = ""
    if agr_a is not None and agr_a.coverage_pct < MIN_COVERAGE:
        caveat = (f". agrA covers only {agr_a.coverage_pct:.0f}% of its reference, so "
                  f"the locus may be truncated and the quorum-sensing system "
                  f"non-functional")
    return TypingResult(
        scheme="agr", call=group, lineage=f"agr group {group}",
        loci_found=1, loci_total=1, alleles=alleles,
        note=f"agr group {group} from the agrB-agrD-agrC' window: "
             f"{hit.identity_pct:.1f}% identity over {hit.coverage_pct:.0f}% of "
             f"{hit.reference} ({hit.strain}){margin}. The group sets which "
             f"autoinducing peptide the isolate makes and which groups it "
             f"cross-inhibits{caveat}",
    )


# ------------------------------------------------------- bla operon and lesions
@dataclass
class Lesion:
    """What is wrong with one gene, in words, with the evidence for it."""

    state: str          # intact | truncated | internal stop | frameshift | ...
    detail: str
    protein_length: int = 0
    expected_length: int = 0


def _describe_lesion(hit: MarkerHit | None, name: str, assembly_seq: dict[str, str],
                     expected_aa: int) -> Lesion:
    """Name what is wrong with a gene rather than only whether it is there."""
    if hit is None or hit.coverage_pct < PARTIAL_COVERAGE:
        return Lesion("absent", f"{name} not detected", 0, expected_aa)
    if hit.identity_pct < MIN_IDENTITY:
        return Lesion("not determined",
                      f"{name}: closest match only {hit.identity_pct:.1f}% identical "
                      f"over {hit.coverage_pct:.0f}% of {hit.reference}, which is a "
                      f"related gene rather than {name}", 0, expected_aa)

    sequence = assembly_seq.get(hit.contig, "")
    lo = min(lo for lo, _ in hit.contig_blocks)
    hi = max(hi for _, hi in hit.contig_blocks)
    edge = hit.at_contig_edge(margin=10)

    if hit.coverage_pct < MIN_COVERAGE:
        missing_5 = hit.reference_blocks[0][0] - 1
        missing_3 = hit.reference_length - hit.reference_blocks[-1][1]
        where = []
        if missing_5 > 0:
            where.append(f"{missing_5} bases missing from the 5' end")
        if missing_3 > 0:
            where.append(f"{missing_3} bases missing from the 3' end")
        reason = ("the alignment runs into the end of the contig, so this is an "
                  "assembly break rather than a demonstrated truncation"
                  if edge else
                  "the contig continues past the alignment, so the gene really is "
                  "truncated")
        return Lesion(
            "truncated by contig break" if edge else "truncated",
            f"{name} covers {hit.coverage_pct:.0f}% of {hit.reference} "
            f"({', '.join(where) or 'interior gap'}); {reason}",
            0, expected_aa)

    if hit.insert_gap > INSERT_GAP:
        return Lesion(
            "interrupted",
            f"{name} is complete but its reference aligns in "
            f"{hit.n_alignments} pieces that sit {hit.insert_gap} bp apart "
            f"on {hit.contig}, with the break at base {hit.insert_at} of "
            f"{hit.reference_length}: something is inserted into the gene",
            0, expected_aa)

    # Full-length and contiguous: read the open reading frame itself.
    gene = sequence[lo - 1:hi].upper()
    if hit.strand == "-":
        gene = revcomp(gene)
    protein = translate(gene)
    stop = protein.find("*")
    if len(gene) % 3:
        return Lesion(
            "frameshifted",
            f"{name} aligns over {hit.coverage_pct:.0f}% of {hit.reference} but spans "
            f"{len(gene)} bases on the contig, not a multiple of three: an indel has "
            f"shifted the frame",
            0, expected_aa)
    if stop != -1 and stop < expected_aa - 1:
        return Lesion(
            "internal stop",
            f"{name} carries a premature stop codon at residue {stop + 1} of an "
            f"expected {expected_aa}, so no functional product is made",
            stop, expected_aa)
    length = len(protein.rstrip("*"))
    return Lesion(
        "intact",
        f"{name} intact: {length} residues, {hit.identity_pct:.1f}% identity over "
        f"{hit.coverage_pct:.0f}% of {hit.reference}",
        length, expected_aa)


#: Expected protein lengths of the operon's three genes, from the pN315 alleles.
BLA_EXPECTED = {"blaZ": 281, "blaR1": 585, "blaI": 126}

#: BlaZ allotype from the residues at Ambler positions 128 and 216. The
#: allotypes were first defined kinetically (Voladri 1996); the residue pair
#: that distinguishes them is tabulated against Ambler's class A standard
#: numbering (Ambler 1991) by Aarris 2025, Table 1.
BLAZ_ALLOTYPE = {("T", "S"): "A", ("K", "N"): "B", ("T", "N"): "C",
                 ("A", "S"): "D", ("L", "S"): "E", ("T", "T"): "F"}

#: Ambler numbering is NOT a fixed offset from the BlaZ sequence: the offset is
#: +7 at the catalytic serine and +9 at the SDN loop, the KSG triad and both
#: typing positions. So the two positions are hard-coded, as 0-based indices
#: into the 281-residue preproenzyme (signal peptide 1-24 included), per
#: Carvajal 2020, which states Ambler 128/216 are "119 and 207 based on our
#: numbering that includes the signal peptide".
BLAZ_AMBLER_128 = 118
BLAZ_AMBLER_216 = 206
#: Class A beta-lactamase motifs, as a check that the frame and the numbering
#: are right before any residue is read off. Ambler Ser70, Ser130-Asp-Asn132
#: and the K[ST]G triad at 234-236.
BLAZ_ANCHORS = ((62, "S"), (120, "SDN"), (224, ("KSG", "KTG")))


def _blaz_allotype(hit: MarkerHit | None,
                   contigs: dict[str, str]) -> tuple[str, str]:
    """The BlaZ allotype, and the evidence for it.

    Returns (allotype, detail). The allotype is "" when it cannot be read,
    which is reported as such rather than guessed: the type decides nothing
    about penicillin on its own, but it does predict the cefazolin inoculum
    effect, so a wrong letter is worse than none.
    """
    if hit is None or not hit.present or hit.insert_gap > INSERT_GAP:
        return "", "no complete, uninterrupted blaZ to type"
    sequence = contigs.get(hit.contig, "")
    lo = min(lo for lo, _ in hit.contig_blocks)
    hi = max(hi for _, hi in hit.contig_blocks)
    gene = sequence[lo - 1:hi].upper()
    if hit.strand == "-":
        gene = revcomp(gene)
    protein = translate(gene).rstrip("*")
    if len(protein) != 281:
        return "", (f"the translated blaZ is {len(protein)} residues, not the 281 of "
                    f"the preproenzyme the allotype positions are numbered against")
    for index, expected in BLAZ_ANCHORS:
        found = protein[index:index + len(expected if isinstance(expected, str)
                                          else expected[0])]
        ok = (found == expected if isinstance(expected, str) else found in expected)
        if not ok:
            return "", (f"the class A beta-lactamase motif expected at residue "
                        f"{index + 1} reads {found!r}, not "
                        f"{expected!r}: the numbering cannot be trusted, so no "
                        f"allotype is assigned")
    one = protein[BLAZ_AMBLER_128]
    two = protein[BLAZ_AMBLER_216]
    allotype = BLAZ_ALLOTYPE.get((one, two), "")
    detail = (f"residue {BLAZ_AMBLER_128 + 1} (Ambler 128) is {one} and residue "
              f"{BLAZ_AMBLER_216 + 1} (Ambler 216) is {two}")
    if not allotype:
        return "", (f"{detail}, which is not one of the described allotype "
                    f"combinations")
    return allotype, detail


def _bla(best: dict[str, MarkerHit], contigs: dict[str, str]) -> list[TypingResult]:
    lesions = {name: _describe_lesion(best.get(name), name, contigs,
                                      BLA_EXPECTED[name])
               for name in ("blaZ", "blaR1", "blaI")}
    alleles = {name: f"{lesion.state} - {lesion.detail}"
               for name, lesion in lesions.items()}
    blaz = lesions["blaZ"]
    hit = best.get("blaZ")

    results: list[TypingResult] = []
    if blaz.state == "absent":
        results.append(TypingResult(
            scheme="blaZ", call="absent", loci_total=1, alleles=alleles,
            note="no blaZ penicillinase gene detected. With no other "
                 "beta-lactamase and no mecA/mecC, benzylpenicillin would be "
                 "expected to test susceptible; a penicillinase-negative result "
                 "should still be confirmed by a zone-edge or nitrocefin test, "
                 "which is what EUCAST and CLSI ask for"))
    elif blaz.state == "intact":
        allotype, detail = _blaz_allotype(hit, contigs)
        alleles["blaZ_allotype_evidence"] = detail
        if allotype:
            alleles["blaZ_allotype"] = allotype
            meaning = (f"blaZ present and intact, allotype {allotype} "
                       f"({detail}). The allotype is the kinetic class of the "
                       f"penicillinase; it does not change the benzylpenicillin "
                       f"call, but types A and C are the ones associated with a "
                       f"strong cefazolin inoculum effect, which matters if "
                       f"cefazolin is being considered for a deep-seated "
                       f"infection. Penicillin resistance is predicted, subject "
                       f"to the regulator status under bla_operon")
        else:
            meaning = (f"blaZ present and intact - {blaz.detail}. The allotype "
                       f"could not be assigned: {detail}. Penicillin resistance "
                       f"is predicted, subject to the regulator status under "
                       f"bla_operon")
        results.append(TypingResult(
            scheme="blaZ", call="present",
            lineage=f"allotype {allotype}" if allotype else "intact",
            loci_found=1, loci_total=1, alleles=alleles, note=meaning))
    else:
        results.append(TypingResult(
            scheme="blaZ", call="present", lineage=blaz.state,
            loci_found=1, loci_total=1, alleles=alleles,
            note=f"blaZ is present but {blaz.state}: {blaz.detail}. A disrupted "
                 f"blaZ does not make penicillinase, so this is not on its own a "
                 f"penicillin-resistance call"))

    # The regulators, and what they mean for the penicillin phenotype.
    r1, bi = lesions["blaR1"], lesions["blaI"]
    if blaz.state == "absent":
        call = "not applicable"
        meaning = ("there is no blaZ for these regulators to control, so their state "
                   "does not affect the penicillin phenotype")
    elif r1.state == "intact" and bi.state == "intact":
        call = "blaR1 intact, blaI intact"
        meaning = ("the operon is complete and inducible: BlaR1 senses the "
                   "beta-lactam and cleaves the BlaI repressor, so penicillinase is "
                   "produced on exposure and benzylpenicillin should test resistant")
    elif r1.state == "absent" and bi.state == "absent":
        call = "blaR1 absent, blaI absent"
        meaning = ("blaZ sits without either regulator. Loss of the BlaI repressor "
                   "leaves blaZ constitutively expressed, so penicillin resistance is "
                   "still expected, but the zone edge may be sharp rather than fuzzy")
    elif bi.state in ("absent", "truncated", "internal stop", "frameshifted",
                      "interrupted") and r1.state == "intact":
        call = f"blaI {bi.state}, blaR1 intact"
        meaning = ("the repressor is gone while the sensor remains, which "
                   "de-represses blaZ: expect constitutive penicillinase and "
                   "penicillin resistance")
    elif r1.state in ("absent", "truncated", "internal stop", "frameshifted",
                      "interrupted") and bi.state == "intact":
        call = f"blaR1 {r1.state}, blaI intact"
        meaning = ("the sensor is lost but the repressor is intact, so blaZ cannot be "
                   "induced. This is the configuration behind blaZ-positive isolates "
                   "that test benzylpenicillin-susceptible, and it is the reason a "
                   "blaZ gene call alone must not be read as resistance. Confirm "
                   "with a penicillin zone-edge test")
    else:
        call = f"blaR1 {r1.state}, blaI {bi.state}"
        meaning = ("the regulator configuration does not match a described pattern; "
                   "the per-gene evidence is in the alleles field")
    truncation_note = ""
    if "contig break" in r1.state or "contig break" in bi.state:
        truncation_note = (". One of the regulators is cut by a contig break rather "
                           "than demonstrably truncated, so this needs a better "
                           "assembly before it is acted on")
    results.append(TypingResult(
        scheme="bla_operon", call=call,
        lineage=f"blaZ {blaz.state}",
        loci_found=sum(1 for lesion in lesions.values() if lesion.state != "absent"),
        loci_total=3, alleles=alleles,
        note=meaning + truncation_note))
    return results


# --------------------------------------------------------------------- toxins
def _presence(best: dict[str, MarkerHit], names: tuple[str, ...]
              ) -> tuple[list[str], list[str], dict[str, str]]:
    """Split a gene set into present, partial, and the evidence for each."""
    present: list[str] = []
    partial: list[str] = []
    alleles: dict[str, str] = {}
    for name in names:
        hit = best.get(name)
        if hit is None:
            alleles[name] = "not detected"
            continue
        alleles[name] = hit.evidence
        if hit.present:
            present.append(name)
        elif hit.partial:
            partial.append(name)
            alleles[name] = f"PARTIAL - {hit.evidence}"
        else:
            alleles[name] = f"below threshold - {hit.evidence}"
    return present, partial, alleles


def _pvl(best: dict[str, MarkerHit]) -> TypingResult:
    present, partial, alleles = _presence(best, ("lukS-PV", "lukF-PV"))
    if len(present) == 2:
        return TypingResult(
            scheme="pvl", call="positive", lineage="lukS-PV + lukF-PV",
            loci_found=2, loci_total=2, alleles=alleles,
            note="both Panton-Valentine leukocidin subunits are present and complete. "
                 "PVL is the marker of the community-associated MRSA lineages and is "
                 "associated with skin and soft-tissue abscesses and with necrotising "
                 "pneumonia. Note that lukF-PV aligns to the lukF/hlgB family at about "
                 f"81% identity in PVL-negative genomes, so a {MIN_IDENTITY:.0f}% "
                 f"identity floor is applied")
    if not present and not partial:
        return TypingResult(
            scheme="pvl", call="negative", loci_total=2, alleles=alleles,
            note="neither lukS-PV nor lukF-PV reached the identity and coverage "
                 "floors. The phage-borne PVL locus is absent")
    found = ", ".join(present + [f"{p} (partial)" for p in partial])
    return TypingResult(
        scheme="pvl", call="not determined", loci_found=len(present), loci_total=2,
        alleles=alleles,
        note=f"only {found} was found. The two subunits are co-transcribed from one "
             f"phage-borne operon, so a single subunit means the locus is split by a "
             f"contig break or partly deleted, not that the isolate is PVL-positive")


def _tst(best: dict[str, MarkerHit]) -> TypingResult:
    hit = best.get("tst")
    alleles = {"tst": hit.evidence if hit else "not detected"}
    if hit is not None and hit.present:
        return TypingResult(
            scheme="tst", call="positive", loci_found=1, loci_total=1, alleles=alleles,
            note="tst present and complete, so the isolate can make toxic shock "
                 "syndrome toxin 1. tst sits on a SaPI pathogenicity island; in a "
                 "patient with shock and a rash this is the determinant to report")
    if hit is not None and hit.partial:
        return TypingResult(
            scheme="tst", call="not determined", loci_total=1, alleles=alleles,
            note=f"tst covers only {hit.coverage_pct:.0f}% of its reference; the SaPI "
                 f"is probably split by a contig break")
    return TypingResult(scheme="tst", call="negative", loci_total=1, alleles=alleles,
                        note="tst not detected; no TSST-1 production expected")


ENTEROTOXINS = ("sea", "seb", "sec", "sed", "see", "seg", "seh", "sei", "sek",
                "sel", "sem", "sen", "seo", "sep", "seq", "seu", "selx")
#: The enterotoxin genes this release actually carries a reference for.
ENTEROTOXINS_TYPED = tuple(n for n in ENTEROTOXINS if n not in ("sed", "see"))


def _enterotoxins(best: dict[str, MarkerHit]) -> TypingResult:
    present, partial, alleles = _presence(best, ENTEROTOXINS_TYPED)
    egc = [g for g in EGC_GENES if g in present]
    classical = [g for g in present if g in ("sea", "seb", "sec", "sed", "see")]
    notes = []
    if egc:
        notes.append(f"the enterotoxin gene cluster (egc) is represented by "
                     f"{', '.join(egc)} of {', '.join(EGC_GENES)}; egc genes travel "
                     f"together on the vSa-beta island and are weakly expressed, so "
                     f"they matter less for food poisoning than the classical genes")
    if classical:
        notes.append(f"the classical staphylococcal food-poisoning enterotoxin "
                     f"gene(s) {', '.join(classical)} are present, which is the "
                     f"finding that counts in a food-poisoning investigation")
    if partial:
        notes.append(f"partial: {', '.join(partial)} - split by a contig break or "
                     f"diverged")
    notes.append("sed, see, sej and ser are not covered by this release's reference "
                 "panel, so a negative here is a negative only for the genes listed "
                 "in the alleles field")
    return TypingResult(
        scheme="enterotoxins",
        call=", ".join(present) if present else "none of the typed genes detected",
        lineage=f"egc {len(egc)}/{len(EGC_GENES)}" if egc else "-",
        loci_found=len(present), loci_total=len(ENTEROTOXINS_TYPED), alleles=alleles,
        note=". ".join(notes))


HAEMOLYSINS = ("hla", "hld", "hlgA", "hlgB", "hlgC", "lukE", "lukD", "lukA", "lukB")


def _haemolysins(best: dict[str, MarkerHit]) -> TypingResult:
    present, partial, alleles = _presence(best, HAEMOLYSINS)
    missing = [n for n in HAEMOLYSINS if n not in present and n not in partial]
    notes = [f"{len(present)} of {len(HAEMOLYSINS)} cytolysin genes complete"]
    if missing:
        notes.append(f"not detected: {', '.join(missing)}")
    if partial:
        notes.append(f"partial: {', '.join(partial)}")
    notes.append("hla (alpha-toxin) and the LukAB/LukGH and LukED leukocidins are "
                 "core-genome genes in almost every isolate, so absence here usually "
                 "means a fragmented assembly rather than a real deletion; hlb is "
                 "reported separately because a prophage usually sits inside it")
    return TypingResult(
        scheme="haemolysins", call=", ".join(present) if present else "none detected",
        loci_found=len(present), loci_total=len(HAEMOLYSINS), alleles=alleles,
        note=". ".join(notes))


def _hlb(best: dict[str, MarkerHit], contigs: dict[str, str]) -> TypingResult:
    hit = best.get("hlb")
    if hit is None or hit.coverage_pct < PARTIAL_COVERAGE:
        return TypingResult(scheme="hlb", call="absent", loci_total=1,
                            alleles={"hlb": hit.evidence if hit else "not detected"},
                            note="no beta-haemolysin gene detected at all, which is "
                                 "unusual: hlb is a core-genome gene. Suspect a "
                                 "fragmented assembly")
    alleles = {"hlb": hit.evidence,
               "alignments": str(hit.n_alignments),
               "contig_gap_bp": str(hit.insert_gap)}
    if hit.insert_gap > INSERT_GAP and hit.coverage_pct >= MIN_COVERAGE:
        return TypingResult(
            scheme="hlb", call="disrupted by prophage integration",
            lineage="Sa3int prophage", loci_found=1, loci_total=1, alleles=alleles,
            note=f"the hlb reference aligns in {hit.n_alignments} pieces "
                 f"whose contig positions are {hit.insert_gap} bp apart, with the "
                 f"break at base {hit.insert_at} of {hit.reference_length}. That is "
                 f"the signature of a beta-haemolysin-converting (Sa3int) prophage "
                 f"integrated inside hlb: the gene is inactivated and the phage "
                 f"usually delivers the immune-evasion cluster, which the iec call "
                 f"reports. Beta-haemolysin is therefore not produced, despite hlb "
                 f"being {hit.coverage_pct:.0f}% covered - which is why coverage "
                 f"alone must not be read as 'hlb present'")
    if hit.coverage_pct >= MIN_COVERAGE and hit.identity_pct >= MIN_IDENTITY:
        return TypingResult(
            scheme="hlb", call="intact", loci_found=1, loci_total=1, alleles=alleles,
            note=f"hlb aligns as one contiguous block over {hit.coverage_pct:.0f}% of "
                 f"{hit.reference} at {hit.identity_pct:.1f}% identity, with no "
                 f"insertion. The gene is uninterrupted, so beta-haemolysin "
                 f"(sphingomyelinase C) is produced and no Sa3int prophage occupies "
                 f"the locus - the minority state among human isolates")
    return TypingResult(
        scheme="hlb", call="not determined", loci_found=1, loci_total=1,
        alleles=alleles,
        note=f"hlb covers only {hit.coverage_pct:.0f}% of its reference at "
             f"{hit.identity_pct:.1f}% identity. A prophage insertion cannot be told "
             f"from a contig break on this evidence")


# -------------------------------------------------------------------- capsule
def _capsule(best: dict[str, MarkerHit]) -> TypingResult:
    five = ("cap5H", "cap5I", "cap5J")
    eight = ("cap8H", "cap8I", "cap8J")
    present5, partial5, alleles5 = _presence(best, five)
    present8, partial8, alleles8 = _presence(best, eight)
    alleles = {**alleles5, **alleles8}
    cap_a = best.get("capA")
    if cap_a is not None:
        alleles["capA"] = cap_a.evidence

    def score(names):
        return sum(best[n].coverage_pct * best[n].identity_pct / 100.0
                   for n in names if n in best)

    if len(present5) == 3 and len(present8) < 3:
        return TypingResult(
            scheme="capsule", call="5", lineage="serotype 5",
            loci_found=3, loci_total=3, alleles=alleles,
            note="all three type-specific genes cap5H, cap5I and cap5J are complete "
                 "and the serotype 8 alleles are not, so the cap locus is serotype 5. "
                 "Serotypes 5 and 8 cover most clinical isolates and the type-specific "
                 "genes HIJ are what distinguishes them; the rest of the locus is "
                 "shared. Note that carrying an intact locus is not the same as "
                 "expressing capsule - several lineages, USA300 among them, carry cap5 "
                 "but are phenotypically non-encapsulated, which this module does not "
                 "assess")
    if len(present8) == 3 and len(present5) < 3:
        return TypingResult(
            scheme="capsule", call="8", lineage="serotype 8",
            loci_found=3, loci_total=3, alleles=alleles,
            note="all three type-specific genes cap8H, cap8I and cap8J are complete "
                 "and the serotype 5 alleles are not, so the cap locus is serotype 8. "
                 "Carrying an intact locus is not the same as expressing capsule, "
                 "which this module does not assess")
    if present5 or present8 or partial5 or partial8:
        return TypingResult(
            scheme="capsule", call="not determined",
            loci_found=len(present5) + len(present8), loci_total=3, alleles=alleles,
            note=f"the type-specific capsule genes are incomplete: serotype 5 score "
                 f"{score(five):.0f}, serotype 8 score {score(eight):.0f} out of 300. "
                 f"Either the locus is split across contigs or this is one of the "
                 f"capsule-null lineages whose cap region differs from both "
                 f"references")
    return TypingResult(
        scheme="capsule", call="not determined", loci_total=3, alleles=alleles,
        note="neither the serotype 5 nor the serotype 8 type-specific genes were "
             "found. Some lineages carry a replaced or deleted cap locus and are "
             "non-encapsulated; this needs the whole locus to resolve")


# ------------------------------------------------------- immune evasion cluster
#: Immune-evasion-cluster types, as the pattern of genes each carries
#: (van Wamel 2006). Order matters only for reporting.
IEC_TYPES: tuple[tuple[str, frozenset], ...] = (
    ("A", frozenset({"sea", "sak", "chp", "scn"})),
    ("B", frozenset({"sak", "chp", "scn"})),
    ("C", frozenset({"chp", "scn"})),
    ("D", frozenset({"sea", "sak", "scn"})),
    ("E", frozenset({"sak", "scn"})),
    ("F", frozenset({"sep", "sak", "chp", "scn"})),
    ("G", frozenset({"sep", "sak", "scn"})),
)
IEC_GENES = ("sea", "sep", "sak", "chp", "scn")


def _iec(best: dict[str, MarkerHit], hlb: TypingResult) -> TypingResult:
    present, partial, alleles = _presence(best, IEC_GENES)
    found = frozenset(present)
    if "scn" not in found:
        note = ("scn, the one gene every immune-evasion cluster carries, was not "
                "detected, so no IEC type is assigned")
        if hlb.call == "intact":
            note += (". This agrees with the hlb call: with no Sa3int prophage in hlb "
                     "there is no cluster to carry")
        return TypingResult(scheme="iec", call="absent", loci_total=len(IEC_GENES),
                            alleles=alleles, note=note)
    exact = [name for name, genes in IEC_TYPES if genes == found]
    if len(exact) == 1:
        extra = ""
        if hlb.call == "intact":
            extra = (". This disagrees with the hlb call of intact: the cluster "
                     "normally arrives on a prophage inside hlb, so check the "
                     "assembly")
        return TypingResult(
            scheme="iec", call=f"type {exact[0]}",
            lineage="+".join(sorted(found)),
            loci_found=len(present), loci_total=len(IEC_GENES), alleles=alleles,
            note=f"immune-evasion cluster type {exact[0]}: "
                 f"{', '.join(sorted(found))} present and "
                 f"{', '.join(g for g in IEC_GENES if g not in found) or 'nothing'} "
                 f"absent. The cluster is carried by the beta-haemolysin-converting "
                 f"prophage and blocks complement (scn), neutrophil chemotaxis (chp) "
                 f"and opsonisation (sak); its type is a stable lineage and "
                 f"host-adaptation marker{extra}")
    return TypingResult(
        scheme="iec", call="not determined", lineage="+".join(sorted(found)),
        loci_found=len(present), loci_total=len(IEC_GENES), alleles=alleles,
        note=f"scn is present, so a cluster is there, but the gene combination "
             f"{'+'.join(sorted(found)) or 'scn only'} does not match any of the "
             f"types A-G defined by van Wamel 2006"
             + (f"; {', '.join(partial)} are only partial, which may be the reason"
                if partial else ""))


# ----------------------------------------------------------------------- ACME
def _acme(best: dict[str, MarkerHit]) -> TypingResult:
    acme = best.get("arcA-ACME")
    native = best.get("arcA-native")
    alleles = {
        "arcA-ACME": acme.evidence if acme else "not detected",
        "arcA-native": native.evidence if native else "not detected",
    }
    if acme is not None and acme.present:
        return TypingResult(
            scheme="acme", call="arc cluster present", lineage="ACME",
            loci_found=1, loci_total=1, alleles=alleles,
            note="the arginine catabolic mobile element arc cluster is present: the "
                 "ACME arcA allele is matched at "
                 f"{acme.identity_pct:.1f}% identity over "
                 f"{acme.coverage_pct:.0f}% of it, well above the 76% it shares with "
                 "the core-genome arginine deiminase every S. aureus carries. ACME "
                 "sits next to SCCmec in USA300 and is thought to aid skin "
                 "colonisation. The ACME type (I, II or III) is not assigned here: "
                 "that needs the opp3 cluster, which this release does not type")
    if acme is not None and acme.partial:
        return TypingResult(
            scheme="acme", call="not determined", loci_total=1, alleles=alleles,
            note=f"the ACME arcA allele covers only {acme.coverage_pct:.0f}% of its "
                 f"reference, which cannot be separated from a contig break")
    return TypingResult(
        scheme="acme", call="absent", loci_total=1, alleles=alleles,
        note="no ACME arc cluster. The arcA hit reported under arcA-native is the "
             "core-genome arginine deiminase, which every S. aureus carries and which "
             "is only 76% identical to the ACME allele - the two must not be "
             "confused")


# ------------------------------------------------------------- clonal complex
def _load_profiles(path: Path) -> dict[int, tuple[int, ...]]:
    out: dict[int, tuple[int, ...]] = {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        handle.readline()
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) < 8 or not parts[0].isdigit():
                continue
            if any(not value.isdigit() for value in parts[1:8]):
                continue
            out[int(parts[0])] = tuple(int(value) for value in parts[1:8])
    return out


def clonal_complex(sequence_type: str, db_dir: Path) -> TypingResult:
    """Assign a clonal complex to a sequence type.

    The complex is the founder ST the query ST is identical to, or differs from
    at a single locus - eBURST's single-locus-variant relation, with the
    complex named after its founder. Connected-component clustering at the same
    threshold is not used because the November 2026 PubMLST table merges nearly
    half of all STs into one component (see :data:`CC_FOUNDERS`).

    The integrator can call this directly with Hydra's existing
    :class:`~hydra_amr.records.MlstCall` sequence type, which avoids typing the
    seven loci twice.
    """
    profiles_file = Path(db_dir) / KEY / "mlst" / "profiles.tsv"
    if not profiles_file.exists():
        return TypingResult(scheme="clonal_complex", call="not determined",
                            note="the PubMLST saureus profile table is not installed")
    if not sequence_type or not str(sequence_type).strip().isdigit():
        return TypingResult(
            scheme="clonal_complex", call="not determined",
            alleles={"ST": str(sequence_type or "-")},
            note=f"no numeric sequence type to work from (ST reported as "
                 f"{sequence_type or '-'}); a novel or incomplete MLST profile cannot "
                 f"be placed in a clonal complex without first being assigned an ST "
                 f"by PubMLST")
    want = int(sequence_type)
    profiles = _load_profiles(profiles_file)
    me = profiles.get(want)
    if me is None:
        return TypingResult(
            scheme="clonal_complex", call="not determined",
            alleles={"ST": str(want)},
            note=f"ST{want} is not in the installed PubMLST profile table "
                 f"({len(profiles)} profiles). Unauthenticated PubMLST access is "
                 f"capped at profiles submitted up to 2024-12-31, so a recently "
                 f"assigned ST will be missing")
    scored = []
    for founder in CC_FOUNDERS:
        theirs = profiles.get(founder)
        if theirs is None:
            continue
        shared = sum(1 for a, b in zip(me, theirs) if a == b)
        scored.append((shared, -founder, founder))
    scored.sort(reverse=True)
    if not scored:
        return TypingResult(scheme="clonal_complex", call="not determined",
                            alleles={"ST": str(want)},
                            note="no founder sequence type was found in the profile "
                                 "table")
    best_shared = scored[0][0]
    winners = [founder for shared, _, founder in scored if shared == best_shared]
    alleles = {
        "ST": str(want),
        "profile": "-".join(str(a) for a in me),
        "loci_shared_with_founder": f"{best_shared}/7",
        "nearest_founders": ", ".join(
            f"CC{founder} ({shared}/7)" for shared, _, founder in scored[:3]),
    }
    if best_shared < 6:
        return TypingResult(
            scheme="clonal_complex", call="no established complex",
            lineage=f"ST{want}", loci_found=7, loci_total=7, alleles=alleles,
            note=f"ST{want} shares at most {best_shared} of 7 loci with any founder "
                 f"of an established clonal complex (nearest: {alleles['nearest_founders']}), "
                 f"so it is neither a founder nor a single-locus variant of one. It is "
                 f"reported as ST{want} without a complex rather than forced into the "
                 f"nearest one")
    if len(winners) > 1:
        names = " or ".join(f"CC{founder}" for founder in sorted(winners))
        return TypingResult(
            scheme="clonal_complex", call=names, lineage=f"ST{want}",
            loci_found=7, loci_total=7, alleles=alleles,
            note=f"ST{want} is a single-locus variant of more than one founder "
                 f"({names}, {best_shared}/7 loci each), so the complex is ambiguous "
                 f"on the single-locus-variant rule alone")
    founder = winners[0]
    relation = ("is the founder of" if founder == want and best_shared == 7
                else "is a single-locus variant of ST%d, the founder of" % founder)
    return TypingResult(
        scheme="clonal_complex", call=f"CC{founder}", lineage=f"ST{want}",
        loci_found=7, loci_total=7, alleles=alleles,
        note=f"ST{want} {relation} CC{founder} ({best_shared}/7 MLST loci shared). "
             f"The clonal complex is how S. aureus clones are named in surveillance - "
             f"CC5, CC8, CC22, CC30 and CC45 account for most hospital MRSA - and is "
             f"a better unit than the ST for comparing isolates across studies")


def _sequence_type(assembly: Path, root: Path, workdir: Path,
                   threads: int) -> tuple[str, dict[str, str], str]:
    """Type the seven MLST loci, for the clonal-complex call.

    Only exact, full-length allele matches count: a 99% match is a new allele,
    and calling it the nearest known one would invent a sequence type.
    """
    alleles_file = root / "mlst" / "alleles.fna"
    profiles_file = root / "mlst" / "profiles.tsv"
    if not alleles_file.exists() or not profiles_file.exists():
        return "", {}, "the PubMLST saureus allele sets are not installed"
    db = workdir / "assembly"
    hsps = blast("blastn", alleles_file, db, workdir / "mlst.tsv",
                 threads=threads, evalue=1e-40, task="megablast",
                 perc_identity=100.0, max_target_seqs=10000)
    exact: dict[str, set[str]] = {}
    for hsp in hsps:
        if hsp.length != hsp.qlen or hsp.nident != hsp.qlen:
            continue
        locus, _, allele = hsp.qseqid.rpartition("_")
        if locus in MLST_LOCI:
            exact.setdefault(locus, set()).add(allele)
    called: dict[str, str] = {}
    for locus in MLST_LOCI:
        hits = exact.get(locus)
        if hits and len(hits) == 1:
            called[locus] = next(iter(hits))
        elif hits:
            called[locus] = "/".join(sorted(hits, key=int))
    if len(called) < len(MLST_LOCI):
        missing = [locus for locus in MLST_LOCI if locus not in called]
        return "", called, (f"no exact full-length allele at {', '.join(missing)}, so "
                            f"the profile is incomplete and no ST can be assigned")
    if any("/" in value for value in called.values()):
        return "", called, "more than one exact allele matched at a locus"
    profiles = _load_profiles(profiles_file)
    want = tuple(int(called[locus]) for locus in MLST_LOCI)
    for st, profile in profiles.items():
        if profile == want:
            return str(st), called, ""
    return "", called, ("the allele profile is complete but matches no ST in the "
                        "installed profile table, so it is a new profile")


# ------------------------------------------------------------------------ run
def run(assembly: Path, db_dir: Path, *, threads: int = 1,
        workdir: Path | None = None) -> list[TypingResult]:
    """Type one *S. aureus* assembly, one :class:`TypingResult` per scheme."""
    assembly = Path(assembly)
    root = Path(db_dir) / KEY
    workdir = Path(workdir) if workdir is not None else assembly.parent
    workdir.mkdir(parents=True, exist_ok=True)

    panel = root / "markers" / "markers.fna"
    if not panel.exists():
        return [TypingResult(
            scheme="saureus", call="not determined",
            note=f"the S. aureus marker panel is not installed at {panel}; run the "
                 f"saureus fetch step")]

    # One BLAST database of the assembly, then every query against it. Putting
    # the markers on the query side is what puts hit coordinates on the
    # assembly's own contigs, which the lesion calls need.
    db = workdir / "assembly"
    run_command(["makeblastdb", "-in", str(assembly), "-dbtype", "nucl",
                 "-out", str(db)])

    lengths = _contig_lengths(assembly)
    hsps = blast("blastn", panel, db, workdir / "saureus.markers.tsv",
                 threads=threads, evalue=1e-10, task="megablast",
                 max_target_seqs=10000)
    best = _best_per_marker(_merge_on_reference(hsps, lengths))
    contigs = {header.split()[0]: seq for header, seq in read_fasta(assembly)}

    results: list[TypingResult] = [
        _spa(assembly, root, workdir, db, threads, best, contigs),
        _agr(best),
    ]
    results.extend(_bla(best, contigs))
    results.append(_pvl(best))
    results.append(_tst(best))
    results.append(_enterotoxins(best))
    results.append(_haemolysins(best))
    hlb = _hlb(best, contigs)
    results.append(hlb)
    results.append(_capsule(best))
    results.append(_iec(best, hlb))
    results.append(_acme(best))

    sequence_type, mlst_alleles, problem = _sequence_type(
        assembly, root, workdir, threads)
    if sequence_type:
        results.append(clonal_complex(sequence_type, db_dir))
    else:
        results.append(TypingResult(
            scheme="clonal_complex", call="not determined",
            alleles={locus: allele for locus, allele in mlst_alleles.items()},
            loci_found=len(mlst_alleles), loci_total=len(MLST_LOCI),
            note=f"{problem}. The clonal complex is derived from the MLST profile, so "
                 f"it cannot be given without one"))
    return results
