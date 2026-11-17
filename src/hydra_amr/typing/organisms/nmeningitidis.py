"""Organism-specific typing for *Neisseria meningitidis*.

What is called, in the order a public-health laboratory needs it:

``cps``
    Serogroup from the capsule locus: A, B, C, W, X, Y, E, Z, the capsule-null
    locus (cnl), and the inactivating events that make a genogroupable isolate
    serologically non-groupable. This is the call that decides the vaccine
    question and the public-health response.
``finetype``
    PorA VR1, PorA VR2 and FetA VR in the standard nomenclature, rendered with
    the serogroup and the sequence type as the full designation recommended for
    meningococcal typing, for example ``B: P1.19,15: F5-1: ST-33 (cc32)``.
``bast``
    Bexsero antigen sequence type: the fHbp, NHBA, NadA, PorA VR1 and PorA VR2
    peptides, with the fHbp variant group and subfamily, the NadA variant and
    the MenDeVAR vaccine-reactivity index.
``mlst``
    Seven-locus sequence type and clonal complex.
``amr_alleles``
    *penA*, *gyrA* and *rpoB* allele identifiers, the unit meningococcal
    resistance is reported in.

Everything is read from an assembly with BLAST+ and pure Python, and every call
carries the evidence it rests on.
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

from ...records import TypingResult
from ...seqio import read_fasta
from ...utils import LOG, HydraError, tempdir
from . import _hinf_nmen_common as common
from .base import OrganismDb, Reference

KEY = "nmeningitidis"
TITLE = "Neisseria meningitidis"
SPECIES: tuple[str, ...] = ("Neisseria meningitidis",)
MLST_SCHEMES: tuple[str, ...] = ("neisseria",)

PUBMLST_DB = "pubmlst_neisseria_seqdef"
#: Scheme ids read from the live API, not guessed: 1 = MLST, 2 = Finetyping
#: antigens (no primary key, so no profile table), 53 = BAST.
MLST_SCHEME = 1
BAST_SCHEME = 53
MLST_LOCI = ("abcZ", "adk", "aroE", "fumC", "gdh", "pdhC", "pgm")

#: Peptide loci, spelled exactly as PubMLST spells them. The lower-case forms
#: answer HTTP 404.
VR_LOCI = ("PorA_VR1", "PorA_VR2", "FetA_VR")
BAST_PEPTIDE_LOCI = ("fHbp_peptide", "NHBA_peptide", "NadA_peptide")
#: The five BAST components, in the order the scheme's profile table uses.
BAST_ORDER = ("fHbp_peptide", "NHBA_peptide", "NadA_peptide", "PorA_VR1", "PorA_VR2")
AMR_LOCI = ("penA", "gyrA", "rpoB")

#: Full-length reference products, to find each antigen gene in the assembly.
#: The PorA and FetA PubMLST loci hold only the variable-region peptide, so the
#: gene itself has to be located against a complete product.
PROTEIN_REFERENCES = {
    "PorA": "P0DH58.1",      # major outer membrane protein P.IA (class 1 protein, PorA)
    "FetA": "WP_010981016.1",  # TonB-dependent siderophore receptor FetA/FrpB
}

#: Peptides contained in the two licensed protein-based vaccines, from the
#: MenDeVAR paper. Used only to explain the index, never to compute it.
VACCINE_PEPTIDES = {
    "Bexsero (4CMenB)": "fHbp peptide 1, NHBA peptide 2, NadA peptide 8, "
                        "PorA VR2 peptide 4 (from the MeNZB outer-membrane vesicle)",
    "Trumenba (rLP2086)": "fHbp peptides 45 and 55",
}

DATA_VERSION = "2026-11-17"

# ----------------------------------------------------------- capsule locus
#: The capsule locus genes, keyed by the PubMLST locus that defines them. The
#: current gene names and the PubMLST loci that carry them are those used by the
#: published CDC meningococcal capsule characterisation tool; the regions and
#: the renaming of the former names (siaA-C to cssA-C, siaD to csb/csc/csw/csy,
#: lipA/lipB to ctrE/ctrF, 29eA-H to cseA-G, xcbA-C to csxA-C, lcbA-C to
#: cslA-C) are those of the capsule locus nomenclature paper.
CAPSULE_LOCI: dict[str, str] = {
    # Region A, serogroup-specific biosynthesis
    "NEIS2157": "csaA", "NEIS2158": "csaB", "NEIS2159": "csaC", "NEIS2160": "csaD",
    "NEIS2161": "csb", "NEIS0051": "csc", "NEIS2162": "csw", "NEIS2163": "csy",
    "NEIS2165": "cseA", "NEIS2166": "cseB", "NEIS2167": "cseC", "NEIS2168": "cseD",
    "NEIS2169": "cseE", "NEIS2170": "cseF", "NEIS2171": "cseG",
    "NEIS2187": "csxA", "NEIS2188": "csxB", "NEIS2189": "csxC",
    "NEIS2184": "cslA", "NEIS2185": "cslB", "NEIS2186": "cslC",
    "NEIS2173": "cszA/cshA", "NEIS2174": "cszB/cshB",
    "NEIS2175": "cszC", "NEIS2176": "cszD",
    "NEIS2177": "cshC", "NEIS2178": "cshD",
    # Region A, the sialic acid pathway shared by B, C, W and Y
    "NEIS0054": "cssA", "NEIS0053": "cssB", "NEIS0052": "cssC",
    "NEIS0050": "cssE", "NEIS2164": "cssF", "NEIS0049": "ctrG",
    # Region C, the ABC capsule transporter; region B, translocation
    "NEIS0055": "ctrA", "NEIS0056": "ctrB", "NEIS0057": "ctrC", "NEIS0058": "ctrD",
    "NEIS0066": "ctrE", "NEIS0067": "ctrF",
    # Flanking region D/D' and the capsule null locus
    "NEIS0059": "tex", "NEIS2743": "cnl",
}

#: Genes that must all be present and intact for a serogroup to be expressed.
#: Region C (ctrA-D) and region B (ctrE/ctrF) are needed by every encapsulated
#: serogroup; the region A set is what names it.
SEROGROUP_ESSENTIAL: dict[str, tuple[str, ...]] = {
    "A": ("csaA", "csaB", "csaC", "csaD"),
    "B": ("csb", "cssA", "cssB", "cssC"),
    "C": ("csc", "cssA", "cssB", "cssC"),
    "W": ("csw", "cssA", "cssB", "cssC"),
    "Y": ("csy", "cssA", "cssB", "cssC"),
    "E": ("cseA", "cseB", "cseC", "cseD", "cseE", "cseG"),
    "X": ("csxA", "csxB", "csxC"),
    "L": ("cslA", "cslB", "cslC"),
    "Z": ("cszA/cshA", "cszB/cshB", "cszC", "cszD"),
    "H": ("cszA/cshA", "cszB/cshB", "cshC", "cshD"),
}
#: Transport genes every encapsulated serogroup needs, plus the flanking marker.
TRANSPORT_ESSENTIAL: tuple[str, ...] = ("ctrA", "ctrB", "ctrC", "ctrD", "ctrE", "ctrF")
#: Genes whose absence does not stop capsule expression but is worth reporting:
#: cssE O-acetylates the serogroup C capsule and cssF the W and Y capsules, so
#: losing one changes the serology without changing the serogroup.
SEROGROUP_ACCESSORY: dict[str, tuple[str, ...]] = {
    "C": ("cssE", "ctrG"), "W": ("cssF", "ctrG"), "Y": ("cssF", "ctrG"),
    "B": ("ctrG",),
}
#: A gene found only in one serogroup. The call needs at least one of these,
#: because cssA-C alone fit B, C, W and Y equally.
SEROGROUP_UNIQUE: dict[str, str] = {
    "csaA": "A", "csaB": "A", "csaC": "A", "csaD": "A",
    "csb": "B", "csc": "C", "csw": "W", "csy": "Y",
    "cseA": "E", "cseB": "E", "cseC": "E", "cseD": "E", "cseE": "E",
    "cseF": "E", "cseG": "E",
    "csxA": "X", "csxB": "X", "csxC": "X",
    "cslA": "L", "cslB": "L", "cslC": "L",
    "cszC": "Z", "cszD": "Z", "cshC": "H", "cshD": "H",
}
#: Minimum identity and coverage for a capsule gene to count as present.
CPS_IDENTITY = 80.0
CPS_COVERAGE = 80.0
#: Below this coverage a gene found at all is reported as fragmented.
CPS_INTACT_COVERAGE = 95.0
#: A gap this large inside an otherwise good alignment is an insertion, which in
#: the capsule locus is usually an insertion sequence.
CPS_INSERTION_BASES = 200

#: Serogroup meaning at the bedside and for the vaccine question.
SEROGROUP_MEANING = {
    "A": "serogroup A; the historical cause of African meningitis-belt epidemics, now "
         "almost eliminated there by MenAfriVac; covered by the quadrivalent ACWY "
         "conjugate vaccines",
    "B": "serogroup B; the commonest cause of invasive disease in Europe, the Americas and "
         "Australasia; not covered by ACWY conjugates, and the target of the "
         "protein-based vaccines Bexsero and Trumenba, for which the BAST and MenDeVAR "
         "results below are the relevant report",
    "C": "serogroup C; epidemic and outbreak-associated, covered by MenC and ACWY "
         "conjugate vaccines; a case in a vaccinated person prompts a vaccine-failure "
         "investigation",
    "W": "serogroup W; the hypervirulent ST-11 complex W lineage has caused a global "
         "resurgence with atypical presentations (pneumonia, septic arthritis, "
         "gastrointestinal illness); covered by ACWY conjugates",
    "X": "serogroup X; sporadic disease and outbreaks in the African meningitis belt; "
         "no licensed monovalent vaccine, though it is included in the pentavalent "
         "ACWXY conjugate",
    "Y": "serogroup Y; disproportionately causes pneumonia and disease in older adults; "
         "covered by ACWY conjugates",
    "E": "serogroup E (formerly 29E); a rare cause of invasive disease, mainly in the "
         "immunocompromised; not covered by any licensed vaccine",
    "Z": "serogroup Z; rare, mainly carriage and occasional disease in the "
         "immunocompromised; not covered by any licensed vaccine",
}

# ------------------------------------------------------------- bibliography
REFERENCES: tuple[Reference, ...] = (
    Reference(
        key="Harrison2013",
        authors="Harrison OB, Claus H, Jiang Y, Bennett JS, Bratcher HB, Jolley KA, et al.",
        year=2013,
        title="Description and nomenclature of Neisseria meningitidis capsule locus",
        journal="Emerg Infect Dis 19:566-573", doi="10.3201/eid1904.111799",
        pmid="23628376",
        supports="the capsule locus region structure (order D-A-C-E-D'-B), the current "
                 "gene nomenclature used here (cssA-C formerly siaA-C, csb/csc/csw/csy "
                 "formerly siaD, ctrE/ctrF formerly lipA/lipB, cseA-G formerly 29eA-H, "
                 "csxA-C formerly xcbA-C, cslA-C formerly lcbA-C) and the rule that the "
                 "polymerase alone separates B, C, W and Y, which all share cssA-C"),
    Reference(
        key="Marjuki2019",
        authors="Marjuki H, Topaz N, Rodriguez-Rivera LD, Ramos E, Potts CC, Chen A, "
                "Retchless AC, Doho GH, Chang HY, Wang X", year=2019,
        title="Whole-genome sequencing for characterization of capsule locus and "
              "prediction of serogroup of invasive meningococcal isolates",
        journal="J Clin Microbiol 57:e01609-18", doi="10.1128/JCM.01609-18",
        pmid="30567750",
        supports="the serogroup prediction rules this module implements: the essential "
                 "gene set per serogroup, the requirement for a serogroup-specific gene, "
                 "and the treatment of partial genes, missing genes, internal stops and "
                 "insertion-sequence disruptions as reasons to report a genogroup with a "
                 "non-groupable phenotype"),
    Reference(
        key="Hammerschmidt1996a",
        authors="Hammerschmidt S, Hilse R, van Putten JP, Gerardy-Schahn R, Unkmeir A, "
                "Frosch M", year=1996,
        title="Modulation of cell surface sialic acid expression in Neisseria meningitidis "
              "via a transposable genetic element",
        journal="EMBO J 15:192-198", doi="10.1002/j.1460-2075.1996.tb00428.x",
        pmid="8598202",
        supports="reversible inactivation of cssA (formerly siaA) by insertion and precise "
                 "excision of IS1301, which is why an insertion-disrupted capsule gene is "
                 "reported as a reversible lesion rather than a permanent loss"),
    Reference(
        key="Hammerschmidt1996b",
        authors="Hammerschmidt S, Muller A, Sillmann H, Muhlenhoff M, Borrow R, Fox A, "
                "van Putten J, Zollinger WD, Gerardy-Schahn R, Frosch M", year=1996,
        title="Capsule phase variation in Neisseria meningitidis serogroup B by "
              "slipped-strand mispairing in the polysialyltransferase gene (siaD): "
              "correlation with bacterial invasion and the outbreak of meningococcal "
              "disease",
        journal="Mol Microbiol 20:1211-1220", doi="10.1111/j.1365-2958.1996.tb02641.x",
        pmid="8809773",
        supports="the poly-cytidine tract in csb (siaD) at nucleotide 89, where seven dC "
                 "residues give an expressing capsule and an insertion or deletion of one "
                 "frameshifts the gene; the basis for reading a frameshifted polymerase as "
                 "phase-variable off"),
    Reference(
        key="Claus2002", authors="Claus H, Maiden MCJ, Maag R, Frosch M, Vogel U",
        year=2002,
        title="Many carried meningococci lack the genes required for capsule synthesis and "
              "transport",
        journal="Microbiology 148:1813-1819", doi="10.1099/00221287-148-6-1813",
        pmid="12055301",
        supports="the capsule null locus (cnl), which replaces region A entirely, and its "
                 "prevalence among carriage isolates"),
    Reference(
        key="Jolley2007", authors="Jolley KA, Brehony C, Maiden MCJ", year=2007,
        title="Molecular typing of meningococci: recommendations for target choice and "
              "nomenclature",
        journal="FEMS Microbiol Rev 31:89-96", doi="10.1111/j.1574-6976.2006.00057.x",
        pmid="17168996",
        supports="the designation this module assembles: 'serogroup: PorA type: FetA type: "
                 "sequence type (clonal complex)', thus 'B: P1.19,15: F5-1: ST-33 (cc32)'"),
    Reference(
        key="Russell2004",
        authors="Russell JE, Jolley KA, Feavers IM, Maiden MCJ, Suker J", year=2004,
        title="PorA variable regions of Neisseria meningitidis",
        journal="Emerg Infect Dis 10:674-678", doi="10.3201/eid1004.030247",
        pmid="15200858",
        supports="the PorA VR nomenclature: VR families grouped at about 80% amino-acid "
                 "identity around a prototype, named family-variant, and rendered with the "
                 "P1. prefix"),
    Reference(
        key="Thompson2003", authors="Thompson EAL, Feavers IM, Maiden MCJ", year=2003,
        title="Antigenic diversity of meningococcal enterobactin receptor FetA, a vaccine "
              "component",
        journal="Microbiology 149:1849-1858", doi="10.1099/mic.0.26131-0", pmid="12855736",
        supports="the FetA variable region and its F<family>-<variant> nomenclature"),
    Reference(
        key="Brehony2016",
        authors="Brehony C, Rodrigues CMC, Borrow R, Smith A, Cunney R, Moxon ER, "
                "Maiden MCJ", year=2016,
        title="Distribution of Bexsero antigen sequence types (BASTs) in invasive "
              "meningococcal disease isolates: implications for immunisation",
        journal="Vaccine 34:4690-4697", doi="10.1016/j.vaccine.2016.08.015", pmid="27521232",
        supports="the Bexsero antigen sequence type and its five components in order - "
                 "fHbp, NHBA, NadA, PorA VR1 and PorA VR2 - with PorA VR1 included for "
                 "extra discrimination although it is not part of the MATS assay"),
    Reference(
        key="Rodrigues2020",
        authors="Rodrigues CMC, Jolley KA, Smith A, Cameron JC, Feavers IM, Maiden MCJ",
        year=2020,
        title="Meningococcal deduced vaccine antigen reactivity (MenDeVAR) index: a rapid "
              "and accessible tool that exploits genomic data in public health and "
              "clinical microbiology applications",
        journal="J Clin Microbiol 59:e02161-20", doi="10.1128/JCM.02161-20", pmid="33055180",
        supports="the four MenDeVAR categories and their meaning, and the vaccine antigen "
                 "peptides (Bexsero fHbp 1, NHBA 2, NadA 8, PorA VR2 4; Trumenba fHbp 45 "
                 "and 55) this module quotes"),
    Reference(
        key="Masignani2003",
        authors="Masignani V, Comanducci M, Giuliani MM, Bambini S, Adu-Bobie J, "
                "Arico B, et al.", year=2003,
        title="Vaccination against Neisseria meningitidis using three variants of the "
              "lipoprotein GNA1870",
        journal="J Exp Med 197:789-799", doi="10.1084/jem.20021911", pmid="12642606",
        supports="the fHbp variant groups 1, 2 and 3 reported alongside the peptide "
                 "identifier"),
    Reference(
        key="Fletcher2004",
        authors="Fletcher LD, Bernfield L, Barniak V, Farley JE, Howell A, Knauf M, et al.",
        year=2004,
        title="Vaccine potential of the Neisseria meningitidis 2086 lipoprotein",
        journal="Infect Immun 72:2088-2100", doi="10.1128/IAI.72.4.2088-2100.2004",
        pmid="15039331",
        supports="the fHbp subfamilies A and B, the second nomenclature reported for the "
                 "same peptide"),
    Reference(
        key="Serruto2010",
        authors="Serruto D, Spadafina T, Ciucchi L, Lewis LA, Ram S, Tontini M, et al.",
        year=2010,
        title="Neisseria meningitidis GNA2132, a heparin-binding protein that induces "
              "protective immunity in humans",
        journal="Proc Natl Acad Sci USA 107:3770-3775", doi="10.1073/pnas.0915162107",
        pmid="20133713",
        supports="NHBA as a vaccine antigen and the peptide nomenclature used for it"),
    Reference(
        key="Bambini2014",
        authors="Bambini S, De Chiara M, Muzzi A, Mora M, Lucidarme J, Brehony C, et al.",
        year=2014,
        title="Neisseria adhesin A variation and revised nomenclature scheme",
        journal="Clin Vaccine Immunol 21:966-971", doi="10.1128/CVI.00825-13",
        pmid="24807056",
        supports="the current NadA variant names (NadA-1, NadA-2/3, NadA-4/5, NadA-6 and "
                 "hybrid variant) that this module reports with the NadA peptide"),
    Reference(
        key="Maiden1998",
        authors="Maiden MCJ, Bygraves JA, Feil E, Morelli G, Russell JE, Urwin R, et al.",
        year=1998,
        title="Multilocus sequence typing: a portable approach to the identification of "
              "clones within populations of pathogenic microorganisms",
        journal="Proc Natl Acad Sci USA 95:3140-3145", doi="10.1073/pnas.95.6.3140",
        pmid="9501229",
        supports="the seven-locus Neisseria MLST scheme (abcZ, adk, aroE, fumC, gdh, pdhC, "
                 "pgm), which was defined on the meningococcus"),
    Reference(
        key="Jolley2018", authors="Jolley KA, Bray JE, Maiden MCJ", year=2018,
        title="Open-access bacterial population genomics: BIGSdb software, the "
              "PubMLST.org website and their applications",
        journal="Wellcome Open Res 3:124", doi="10.12688/wellcomeopenres.14826.1",
        pmid="30345391",
        supports="the PubMLST database every allele definition, peptide definition, "
                 "profile and clonal complex used here is downloaded from; PubMLST asks "
                 "that this be cited"),
)

#: Quoted once, used by every PubMLST database entry below.
_PUBMLST_LICENCE = (
    "PubMLST Terms and Conditions (https://pubmlst.org/terms-conditions, page version "
    "dated 2026-05-13, retrieved 2026-11-17). It is not a named licence but a University "
    "of Oxford terms-of-use document. The deciding sentence: \"All data submitted on or "
    "before 31 December 2024 may continue to be downloaded, used, and redistributed "
    "without restriction, subject only to proper scientific citation and "
    "acknowledgement.\" The terms add that \"The restrictions described in Sections 4 and "
    "5 do not apply to pre-2025 data\", and those are the sections that would otherwise "
    "forbid redistribution: \"not to redistribute, publicly post, publish, transfer, or "
    "aggregate significant portions of the data into any public, semi-public, or "
    "commercial database, knowledgebase, software system, or service\"."
)
_PUBMLST_NOTES = (
    "Redistribution is permitted for the pre-2025 snapshot an anonymous download returns, "
    "and only for that: the sentence quoted above places exactly that data outside every "
    "restriction, subject to citation and acknowledgement, which are given below. "
    "The redistribution right rests entirely on the download being anonymous. The "
    "unauthenticated API answers every request with the message \"Please note that you "
    "are currently restricted to accessing data that was submitted on or prior to "
    "2024-12-31. Please authenticate to access the full dataset.\", so an anonymous "
    "download can only return pre-2025 data, which the terms place outside the "
    "restrictions. fetch() therefore sends no credentials, and an API key must never be "
    "added: authenticating pulls post-2024 data and with it Section 4's no-redistribution "
    "clause and Section 5's statement that \"commercial services and software tools that "
    "enable the use or the downloading of data ... will require a specific licence "
    "contract\". PubMLST requires this acknowledgement wherever the data is used: \"This "
    "publication made use of the PubMLST website (https://pubmlst.org/) sited at the "
    "University of Oxford. The development of that website was funded by the Wellcome "
    "Trust.\" The terms also reserve the right to \"modify nomenclature, allele "
    "definitions, or database structure\" without notice, so the snapshot date matters."
)

DATABASES: tuple[OrganismDb, ...] = (
    OrganismDb(
        name="capsule",
        source=f"https://rest.pubmlst.org/db/{PUBMLST_DB}/loci/<locus>/alleles_fasta for "
               f"the {len(CAPSULE_LOCI)} capsule loci NEIS0045-NEIS2743 listed in "
               f"CAPSULE_LOCI",
        version=f"anonymous snapshot retrieved {DATA_VERSION}",
        licence=_PUBMLST_LICENCE, redistribution="permitted",
        reference="Harrison2013",
        notes="About 2 MB. The choice of loci and the gene-name mapping follow the CDC "
              "capsule characterisation tool and the nomenclature paper. " + _PUBMLST_NOTES,
    ),
    OrganismDb(
        name="peptides",
        source=f"https://rest.pubmlst.org/db/{PUBMLST_DB}/loci/<locus>/alleles_fasta for "
               f"{', '.join(VR_LOCI + BAST_PEPTIDE_LOCI)}, the extended attributes of "
               f"fHbp_peptide and NadA_peptide, and scheme {BAST_SCHEME} (BAST) profiles",
        version=f"anonymous snapshot retrieved {DATA_VERSION}",
        licence=_PUBMLST_LICENCE, redistribution="permitted",
        reference="Brehony2016",
        notes="About 1.6 MB. The locus names are case-sensitive: the lower-case forms "
              "answer HTTP 404. The BAST profile table carries the two MenDeVAR columns, "
              "so the vaccine-reactivity index comes free with it; there is no per-allele "
              "MenDeVAR endpoint. " + _PUBMLST_NOTES,
    ),
    OrganismDb(
        name="mlst",
        source=f"https://rest.pubmlst.org/db/{PUBMLST_DB} (loci "
               f"{', '.join(MLST_LOCI)} and scheme {MLST_SCHEME} profiles)",
        version=f"anonymous snapshot retrieved {DATA_VERSION}",
        licence=_PUBMLST_LICENCE, redistribution="permitted",
        reference="Jolley2018",
        notes="About 7 MB, the largest part of this module's data. The profile table "
              "carries a clonal_complex column populated for about 61% of sequence types, "
              "written as 'ST-32 complex', which this module renders as cc32. " +
              _PUBMLST_NOTES,
    ),
    OrganismDb(
        name="amr",
        source=f"https://rest.pubmlst.org/db/{PUBMLST_DB}/loci/<locus>/alleles_fasta for "
               f"{', '.join(AMR_LOCI)}",
        version=f"anonymous snapshot retrieved {DATA_VERSION}",
        licence=_PUBMLST_LICENCE, redistribution="permitted",
        reference="Jolley2018",
        notes="About 1.7 MB. Allele identifiers only; no allele-to-phenotype mapping is "
              "shipped or implied. " + _PUBMLST_NOTES,
    ),
    OrganismDb(
        name="antigens",
        source="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="
               + ",".join(sorted(PROTEIN_REFERENCES.values())),
        version=f"{', '.join(sorted(PROTEIN_REFERENCES.values()))}; "
                f"retrieved {DATA_VERSION}",
        licence="none stated. NCBI states that it 'does not restrict the use or "
                "distribution of the data it contains'; no licence is attached to these "
                "records.",
        redistribution="permitted",
        reference="Russell2004",
        notes="Two protein records, about 2 kB. They exist only to locate porA and fetA in "
              "the assembly, because the PubMLST PorA and FetA loci hold the variable "
              "region peptide and not the whole product.",
    ),
)


# ----------------------------------------------------------------- fetching
def fetch(db_dir: Path) -> None:
    """Download and build this module's reference data into ``db_dir/nmeningitidis``.

    PubMLST is queried anonymously and without an API key, deliberately: the
    unauthenticated interface serves only data submitted on or before
    2024-12-31, and the terms place exactly that data outside the
    redistribution restrictions. Authenticating would pull post-2024 data and
    with it the clause forbidding redistribution.
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
    _fetch_antigen_references(root)
    _fetch_peptides(root)
    _fetch_mlst(root)
    _fetch_amr(root)
    _fetch_capsule(root)
    stamp.write_text(json.dumps({
        "version": DATA_VERSION,
        "retrieved": date.today().isoformat(),
        "pubmlst_database": PUBMLST_DB,
        "pubmlst_access": "anonymous; serves data submitted on or before 2024-12-31",
    }, indent=2))


def _fetch_antigen_references(root: Path) -> None:
    """Download the full-length PorA and FetA products."""
    work = root / "antigens"
    work.mkdir(parents=True, exist_ok=True)
    raw = work / "references.raw.faa"
    common.ncbi_protein_fasta(list(PROTEIN_REFERENCES.values()), raw)
    sequences = [sequence for _header, sequence in read_fasta(raw)]
    raw.unlink()
    names = list(PROTEIN_REFERENCES)
    if len(sequences) != len(names):
        raise HydraError(f"{KEY}: expected {len(names)} antigen reference proteins, "
                         f"got {len(sequences)}")
    # efetch answers a comma-separated id list in the order it was asked.
    common.write_subset(list(zip(names, sequences)), work / "antigens.faa")


def _fetch_peptides(root: Path) -> None:
    """Download the variable-region and vaccine-antigen peptide definitions."""
    work = root / "peptides"
    work.mkdir(parents=True, exist_ok=True)
    for locus in VR_LOCI + BAST_PEPTIDE_LOCI:
        common.pubmlst_locus_fasta(PUBMLST_DB, locus, work / f"{locus}.fas")
    common.pubmlst_profiles(PUBMLST_DB, BAST_SCHEME, work / "bast_profiles.txt")
    # The extended attributes carry the nomenclature a clinician reads: the
    # fHbp variant group and subfamily, and the NadA variant.
    for locus, fields in (("fHbp_peptide", ("Novartis_variant_family", "Novartis_variant",
                                            "Pfizer_subfamily", "Pfizer_variant")),
                          ("NadA_peptide", ("NadA_variant",))):
        _fetch_extended(work, locus, fields)


def _fetch_extended(work: Path, locus: str, fields: tuple[str, ...]) -> None:
    """Save the extended attributes of one peptide locus as a small table."""
    raw = work / f"{locus}.extended.json"
    url = (f"{common.PUBMLST_REST}/db/{PUBMLST_DB}/loci/{locus}/alleles"
           f"?return_all=1&include_records=1&extended=1")
    try:
        common.download(url, raw)
        payload = json.loads(raw.read_text())
    except (HydraError, ValueError) as exc:
        LOG.warning("%s: extended attributes for %s unavailable (%s); the variant group "
                    "will be reported as not determined", KEY, locus, exc)
        return
    finally:
        if raw.exists():
            raw.unlink()
    rows = payload.get("alleles") if isinstance(payload, dict) else payload
    out = work / f"{locus}.extended.tsv"
    with open(out, "w") as handle:
        handle.write("allele_id\t" + "\t".join(fields) + "\n")
        for record in rows or []:
            if not isinstance(record, dict):
                continue
            allele = str(record.get("allele_id", "")).strip()
            if not allele:
                continue
            values = [str(record.get(field, "") or "").strip() for field in fields]
            if any(values):
                handle.write(allele + "\t" + "\t".join(values) + "\n")


def _fetch_mlst(root: Path) -> None:
    """Download the PubMLST *Neisseria* seven-locus scheme."""
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
    common.pubmlst_profiles(PUBMLST_DB, MLST_SCHEME, work / "profiles.txt")


def _fetch_capsule(root: Path) -> None:
    """Build the capsule reference set from the PubMLST capsule loci.

    Every allele of every capsule locus is kept, and both the gene and its
    product are stored. One allele per locus is not enough: the
    serogroup-specific polymerases are the most polymorphic genes in the locus,
    and the serogroup B strain MC58 aligns to only 43% of the first *csb*
    allele, which would be read as a disrupted gene in a strain that expresses
    capsule perfectly well. Presence is decided on the product, because a
    translated search tolerates the synonymous divergence; the nucleotide copy
    is kept so a hit can be confirmed and an insertion measured in bases.
    """
    work = root / "capsule"
    work.mkdir(parents=True, exist_ok=True)
    genes: list[tuple[str, str]] = []
    proteins: list[tuple[str, str]] = []
    loci_found = 0
    for locus, gene in sorted(CAPSULE_LOCI.items()):
        raw = work / f"{locus}.raw.fas"
        try:
            common.pubmlst_locus_fasta(PUBMLST_DB, locus, raw)
        except HydraError as exc:
            LOG.warning("%s: capsule locus %s (%s) unavailable (%s)", KEY, locus, gene, exc)
            continue
        alleles = [(header.split()[0], sequence) for header, sequence in read_fasta(raw)]
        raw.unlink()
        if not alleles:
            continue
        loci_found += 1
        # One reference per distinct product. The conserved transport genes have
        # thousands of nucleotide alleles that translate to a handful of
        # proteins, and keeping them all would make this database forty times
        # larger for no gain in a translated search. The lowest-numbered allele
        # of each distinct product is kept, and the nucleotide copy kept with it.
        # The conserved transport and flanking genes need only a handful of
        # references for a presence call; the serogroup-specific region A genes
        # are the discriminators and keep all of theirs.
        cap = 250 if (gene in SEROGROUP_UNIQUE or gene.startswith("css")
                      or gene == "cnl") else 25
        seen: set[str] = set()
        kept = 0
        for name, sequence in sorted(alleles, key=lambda item: _allele_rank(item[0])):
            if kept >= cap:
                break
            protein = common.best_reading_frame(sequence)
            if not protein or protein in seen:
                continue
            seen.add(protein)
            allele = name.rsplit("_", 1)[-1]
            label = f"{gene}|{locus}_{allele}"
            genes.append((label, sequence))
            proteins.append((label, protein))
            kept += 1
        LOG.debug("%s: %s (%s) %d alleles -> %d distinct products",
                  KEY, gene, locus, len(alleles), kept)
    if loci_found < len(CAPSULE_LOCI) // 2:
        raise HydraError(f"{KEY}: only {loci_found} of {len(CAPSULE_LOCI)} capsule loci "
                         f"downloaded; the serogroup could not be called from that")
    common.write_subset(genes, work / "cps_genes.fna")
    common.makeblastdb(work / "cps_genes.fna", "nucl")
    common.write_subset(proteins, work / "cps_proteins.faa")


def _allele_rank(name: str) -> tuple:
    """Sort allele record names so the lowest-numbered allele comes first."""
    number = name.rsplit("_", 1)[-1]
    return (0, int(number), "") if number.isdigit() else (1, 0, name)


def _fetch_amr(root: Path) -> None:
    """Download the *penA*, *gyrA* and *rpoB* allele definitions."""
    work = root / "amr"
    work.mkdir(parents=True, exist_ok=True)
    records: list[tuple[str, str]] = []
    for locus in AMR_LOCI:
        raw = work / f"{locus}.raw.fas"
        common.pubmlst_locus_fasta(PUBMLST_DB, locus, raw)
        for header, sequence in read_fasta(raw):
            name = header.split()[0]
            if not name.startswith(f"{locus}_"):
                name = f"{locus}_{name}"
            records.append((name, sequence))
        raw.unlink()
    common.write_subset(records, work / "amr.fna")
    common.makeblastdb(work / "amr.fna", "nucl")


# ------------------------------------------------------------------- typing
def run(assembly: Path, db_dir: Path, *, threads: int = 1,
        workdir: Path | None = None, mlst=None) -> list[TypingResult]:
    """Type one *N. meningitidis* assembly, returning one result per scheme."""
    assembly = Path(assembly)
    root = Path(db_dir) / KEY
    if workdir is not None:
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        return _run_in(assembly, root, workdir, threads, mlst)
    with tempdir(prefix="hydra.nmeningitidis.") as work:
        return _run_in(assembly, root, work, threads, mlst)


def _run_in(assembly: Path, root: Path, work: Path, threads: int,
            mlst_call=None) -> list[TypingResult]:
    contigs = common.load_contigs(assembly)
    local = work / "assembly.fna"
    common.write_subset(contigs.items(), local)
    common.makeblastdb(local, "nucl")

    capsule = _call_capsule(local, contigs, root, work, threads)
    peptides = _call_peptides(local, contigs, root, work, threads)
    mlst = _call_mlst(local, root, work, threads, mlst_call)
    finetype = _call_finetype(capsule, peptides, mlst)
    bast = _call_bast(peptides, root)
    amr = _call_amr(local, root, work, threads)
    return [capsule, finetype, bast, mlst, amr]


# ---------------------------------------------------------------- capsule
def _call_capsule(assembly: Path, contigs: dict[str, str], root: Path, work: Path,
                  threads: int) -> TypingResult:
    """Call the serogroup, and say why an isolate with a locus is non-groupable."""
    proteins = root / "capsule" / "cps_proteins.faa"
    if not proteins.exists():
        return TypingResult(scheme="cps", call="not determined",
                            note=f"capsule reference data missing; "
                                 f"run 'hydra db download {KEY}'")
    lengths = {header.split()[0]: len(sequence)
               for header, sequence in read_fasta(proteins)}
    hits = common.screen_proteins(assembly, proteins, work, threads=threads,
                                  min_identity=CPS_IDENTITY, min_coverage=0.0,
                                  tag="cps")
    found: dict[str, dict] = {}
    for hit in hits:
        gene = hit.gene.split("|")[0]
        if hit.coverage_pct < 40.0:
            continue
        cds = common.extract_cds(contigs, hit.contig, hit.start, hit.end, hit.strand,
                                 ref_offset=max(0, hit.query_start - 1))
        residues = hit.coverage_pct / 100.0 * lengths.get(hit.gene, 0)
        span = abs(hit.end - hit.start) + 1
        expected = 3 * lengths.get(hit.gene, 0)
        record = {
            "hit": hit, "cds": cds,
            "intact": hit.coverage_pct >= CPS_INTACT_COVERAGE,
            "insertion": span - expected if span - expected >= CPS_INSERTION_BASES else 0,
            "residues": residues,
        }
        current = found.get(gene)
        score = hit.coverage_pct * hit.identity_pct
        if current is None or score > (current["hit"].coverage_pct
                                       * current["hit"].identity_pct):
            found[gene] = record
    present = set(found)

    # Score each serogroup on its own region A set plus the shared transporter.
    scored: list[tuple[int, int, str]] = []
    for group, essential in SEROGROUP_ESSENTIAL.items():
        region_a = sum(1 for gene in essential if gene in present)
        unique = sum(1 for gene in essential
                     if gene in present and SEROGROUP_UNIQUE.get(gene) == group)
        if region_a:
            scored.append((unique, region_a, group))
    scored.sort(reverse=True)

    alleles: dict[str, str] = {}
    transport = {gene: found.get(gene) for gene in TRANSPORT_ESSENTIAL}
    for gene in TRANSPORT_ESSENTIAL + ("tex", "cnl"):
        alleles[gene] = _cps_evidence(found.get(gene))

    if not scored or scored[0][0] == 0:
        # No serogroup-specific gene. Either the capsule null locus, or a locus
        # whose region A is gone, or not a meningococcus.
        if "cnl" in present:
            alleles["region_A"] = "absent; the capsule null locus (cnl) is present instead"
            return TypingResult(
                scheme="cps", call="cnl", lineage="capsule null locus (cnl)",
                alleles=alleles, loci_found=len(present),
                loci_total=len(TRANSPORT_ESSENTIAL) + 2,
                note="the capsule null locus replaces region A, so the isolate has no "
                     "capsule biosynthesis genes at all and can never be grouped. cnl "
                     "isolates are common in asymptomatic carriage, are usually avirulent "
                     "because the capsule is the main defence against complement, and are "
                     "not covered by any polysaccharide vaccine. Invasive disease caused by "
                     "a cnl isolate is reported but rare and is associated with "
                     "complement deficiency in the patient, so a cnl isolate from a sterile "
                     "site should prompt a complement screen")
        if any(transport.values()):
            alleles["region_A"] = ("no serogroup-specific biosynthesis gene found, although "
                                   "capsule transport genes are present")
            return TypingResult(
                scheme="cps", call="NG", lineage="non-groupable: region A absent",
                alleles=alleles, loci_found=len(present),
                loci_total=len(TRANSPORT_ESSENTIAL) + 4,
                note="the capsule transport genes are present but no serogroup-specific "
                     "biosynthesis gene was found, so the serogroup cannot be named and no "
                     "capsule can be made. A whole-region-A deletion is one of the ways a "
                     "carriage isolate becomes non-groupable")
        return TypingResult(
            scheme="cps", call="not determined", lineage="-", alleles=alleles,
            loci_found=0, loci_total=len(TRANSPORT_ESSENTIAL) + 4,
            note="no capsule locus gene was found at all. Confirm the species: this is "
                 "what a commensal Neisseria looks like against a meningococcal capsule "
                 "reference")

    unique, region_a, group = scored[0]
    essential = SEROGROUP_ESSENTIAL[group]
    for gene in essential + SEROGROUP_ACCESSORY.get(group, ()):
        alleles[gene] = _cps_evidence(found.get(gene))
    rivals = [f"{other} ({count}/{len(SEROGROUP_ESSENTIAL[other])} region A genes)"
              for _u, count, other in scored[1:4] if count]
    alleles["genogroup"] = group
    alleles["region_A_support"] = (f"{region_a}/{len(essential)} genes of serogroup {group}, "
                                   f"{unique} of them found only in {group}")
    if rivals:
        alleles["other_serogroups_matched"] = ", ".join(rivals)

    faults: list[str] = []
    for gene in essential + TRANSPORT_ESSENTIAL:
        record = found.get(gene)
        if record is None:
            faults.append(f"{gene} not found")
            continue
        cds = record["cds"]
        if record["insertion"]:
            faults.append(
                f"{gene} interrupted by an insertion of about {record['insertion']} bp "
                f"(in serogroups W and Y this is typically IS1301, which inserts in the "
                f"capsule locus and switches expression off)")
        elif not record["intact"]:
            faults.append(f"{gene} fragmented ({record['hit'].coverage_pct:.0f}% of the "
                          f"reference product aligned)")
        if cds is not None and cds.internal_stops:
            faults.append(f"{gene} carries {cds.internal_stops} internal stop codon(s), so "
                          f"its product is truncated")
        if cds is not None and cds.truncated_at_contig_end:
            faults.append(f"{gene} runs into a contig end, so it could not be read in full")

    meaning = SEROGROUP_MEANING.get(group, f"serogroup {group}")
    if not faults:
        return TypingResult(
            scheme="cps", call=group, lineage=f"serogroup {group}; {meaning}",
            alleles=alleles, loci_found=len(present),
            loci_total=len(essential) + len(TRANSPORT_ESSENTIAL),
            note=f"every gene serogroup {group} needs is present and intact in regions A, "
                 f"B and C, so the isolate is expected to express a serogroup {group} "
                 f"capsule and to agglutinate with {group} antiserum")
    contig_only = all("runs into a contig end" in fault or "could not be read" in fault
                      for fault in faults)
    if contig_only:
        return TypingResult(
            scheme="cps", call=group, lineage=f"serogroup {group}; {meaning}",
            alleles=alleles, loci_found=len(present),
            loci_total=len(essential) + len(TRANSPORT_ESSENTIAL),
            note=f"genogroup {group}. " + "; ".join(faults) +
                 ". The only faults found are at contig boundaries, which is an assembly "
                 "limitation rather than evidence of a disrupted gene, so the serogroup is "
                 "reported as determined; expression cannot be confirmed from this assembly")
    return TypingResult(
        scheme="cps", call="NG",
        lineage=f"genogroup {group}, serologically non-groupable",
        alleles=alleles, loci_found=len(present),
        loci_total=len(essential) + len(TRANSPORT_ESSENTIAL),
        note=f"the isolate carries the serogroup {group} capsule locus but cannot express "
             f"it: " + "; ".join(faults) +
             f". This is the common pattern in carriage isolates, and it matters in two "
             f"directions: the laboratory will report it as non-groupable while the "
             f"genogroup is {group}, and the lesion can revert, since an insertion-sequence "
             f"or slipped-strand switch is reversible. Report the genogroup and the lesion, "
             f"not the serogroup alone. {meaning}")


def _cps_evidence(record: dict | None) -> str:
    """One phrase for a capsule gene and the state it was found in."""
    if record is None:
        return "not detected"
    hit = record["hit"]
    cds = record["cds"]
    parts = [f"{hit.identity_pct:.1f}% amino-acid identity over {hit.coverage_pct:.0f}% of "
             f"the reference product", f"{hit.contig}:{hit.start}-{hit.end}({hit.strand})"]
    if record["insertion"]:
        parts.append(f"about {record['insertion']} bp of inserted sequence")
    if cds is not None and cds.internal_stops:
        parts.append(f"{cds.internal_stops} internal stop codon(s)")
    if cds is not None and cds.truncated_at_contig_end:
        parts.append("runs into a contig end")
    return "; ".join(parts)


# --------------------------------------------------------------- peptides
def _call_peptides(assembly: Path, contigs: dict[str, str], root: Path, work: Path,
                   threads: int) -> dict[str, tuple[str, str]]:
    """Identify every peptide antigen, returning {locus: (identifier, note)}.

    PorA and FetA are typed on a short variable region, so the gene is located
    against a full-length product, its coding sequence pulled out of the
    assembly and translated, and the variable region then read off by finding
    the defined peptide inside that translation. fHbp, NHBA and NadA are typed
    on the whole product, so the recovered protein is matched against the
    peptide definitions directly.
    """
    out: dict[str, tuple[str, str]] = {}
    peptides = root / "peptides"
    products = _recover_products(assembly, contigs, root, work, threads)

    for locus, window in (("PorA_VR1", (0.0, 0.6)), ("PorA_VR2", (0.3, 1.0))):
        path = peptides / f"{locus}.fas"
        protein, why = products.get("PorA", ("", "PorA not recovered"))
        if not path.exists():
            out[locus] = ("", f"not determined: {locus} definitions not installed")
        elif not protein:
            out[locus] = ("", f"not determined: {why}")
        else:
            out[locus] = common.find_variable_region(
                protein, common.peptide_index(path, locus), window=window)
    path = peptides / "FetA_VR.fas"
    protein, why = products.get("FetA", ("", "FetA not recovered"))
    if not path.exists():
        out["FetA_VR"] = ("", "not determined: FetA_VR definitions not installed")
    elif not protein:
        out["FetA_VR"] = ("", f"not determined: {why}")
    else:
        out["FetA_VR"] = common.find_variable_region(
            protein, common.peptide_index(path, "FetA_VR"))

    for locus in BAST_PEPTIDE_LOCI:
        path = peptides / f"{locus}.fas"
        if not path.exists():
            out[locus] = ("", f"not determined: {locus} definitions not installed")
            continue
        index = common.peptide_index(path, locus)
        protein, why = products.get(locus, ("", ""))
        if protein:
            # The fHbp definitions start at the cysteine of the mature
            # lipoprotein, after the signal peptide is cleaved, so the defined
            # peptide is a substring of the translated gene rather than equal to
            # it. Matching the longest contained definition handles that and
            # still gives an exact answer for NHBA and NadA, whose definitions
            # include the signal peptide.
            identifier, note = common.find_variable_region(protein, index)
            if identifier != "new":
                out[locus] = (identifier, note)
                continue
            out[locus] = common.match_peptide(protein, index)
            continue
        if not protein:
            # fHbp, NHBA and NadA are not universal: nadA in particular is
            # absent from most carriage isolates and from whole lineages, and
            # absence is a result, not a failure.
            reason = why or "no hit to any defined peptide"
            out[locus] = ("0", f"gene not detected ({reason}); for NadA this is common "
                               f"and expected, and PubMLST records an absent locus as 0")
            continue
        out[locus] = common.match_peptide(protein, index)
    return out


def _recover_products(assembly: Path, contigs: dict[str, str], root: Path, work: Path,
                      threads: int) -> dict[str, tuple[str, str]]:
    """Pull each antigen's product out of the assembly, with its caveat."""
    out: dict[str, tuple[str, str]] = {}
    references = work / "antigen_references.faa"
    records: list[tuple[str, str]] = []
    antigens = root / "antigens" / "antigens.faa"
    if antigens.exists():
        records.extend((header.split()[0], sequence)
                       for header, sequence in read_fasta(antigens))
    # The peptide definitions themselves are full-length products for these
    # three antigens, so their first allele is the natural search reference.
    for locus in BAST_PEPTIDE_LOCI:
        path = root / "peptides" / f"{locus}.fas"
        if not path.exists():
            continue
        first = next(((h, s) for h, s in read_fasta(path)), None)
        if first:
            records.append((locus, first[1]))
    if not records:
        return out
    common.write_subset(records, references)
    hits = common.screen_proteins(assembly, references, work, threads=threads,
                                  min_identity=65.0, min_coverage=50.0, tag="antigens")
    best = common.best_per_gene(hits)
    for name, hit in best.items():
        cds = common.extract_cds(contigs, hit.contig, hit.start, hit.end, hit.strand,
                                 ref_offset=max(0, hit.query_start - 1))
        if cds is None or not cds.protein:
            out[name] = ("", f"{name} located at {hit.evidence()} but no reading frame "
                             f"could be recovered")
            continue
        caveats = []
        if cds.internal_stops:
            caveats.append(f"{cds.internal_stops} internal stop codon(s)")
        if cds.truncated_at_contig_end:
            caveats.append("runs into a contig end")
        out[name] = (cds.protein, "; ".join(caveats))
    for name in ("PorA", "FetA") + BAST_PEPTIDE_LOCI:
        out.setdefault(name, ("", f"{name} not found at >=65% identity over >=50% of its "
                                  f"reference product"))
    return out


# --------------------------------------------------------------- fine type
def _call_finetype(capsule: TypingResult, peptides: dict[str, tuple[str, str]],
                   mlst: TypingResult) -> TypingResult:
    """Assemble the standard meningococcal designation.

    The nomenclature is ``serogroup: PorA type: FetA type: sequence type
    (clonal complex)``. PubMLST publishes no profile table for the fine-typing
    scheme, so the string is assembled here from the three independent
    variable-region calls.
    """
    vr1, vr1_note = peptides.get("PorA_VR1", ("", "not determined"))
    vr2, vr2_note = peptides.get("PorA_VR2", ("", "not determined"))
    feta, feta_note = peptides.get("FetA_VR", ("", "not determined"))
    serogroup = capsule.call if capsule.call in SEROGROUP_MEANING else (
        "cnl" if capsule.call == "cnl" else "NG")
    pora = f"P1.{vr1 or '?'},{vr2 or '?'}"
    fet = feta or "?"
    call = f"{serogroup}:{pora}:{fet}"
    full = call
    if mlst.call != "not determined":
        cc = _short_complex(mlst.lineage)
        full = f"{call}:{mlst.call}" + (f"({cc})" if cc else "")
    alleles = {
        "PorA_VR1": vr1 or "not determined",
        "PorA_VR2": vr2 or "not determined",
        "FetA_VR": feta or "not determined",
        "full_designation": full,
        "nomenclature": "serogroup : PorA type : FetA type : sequence type (clonal "
                        "complex), as recommended for molecular typing of meningococci; "
                        "PorA variable regions are surface loops 1 and 4 of the class 1 "
                        "porin and are named family-variant, so P1.7-2,4 is VR1 family 7 "
                        "variant 2 with VR2 family 4 prototype",
    }
    notes = [note for note in (vr1_note, vr2_note, feta_note) if note]
    if not notes:
        notes.append("the fine type is what an outbreak investigation compares between "
                     "isolates; it resolves far below the serogroup and is the unit used "
                     "for case clustering in meningococcal surveillance")
    return TypingResult(
        scheme="finetype", call=call, lineage=full, alleles=alleles,
        loci_found=sum(1 for value in (vr1, vr2, feta) if value and value != "new"),
        loci_total=3, note=". ".join(notes))


def _short_complex(clonal_complex: str) -> str:
    """Turn PubMLST's ``ST-41/44 complex`` into the ``cc41/44`` of the designation."""
    text = (clonal_complex or "").strip()
    if not text or text == "-":
        return ""
    if text.lower().startswith("st-") and text.lower().endswith("complex"):
        return "cc" + text[3:-len("complex")].strip()
    return text


# -------------------------------------------------------------------- BAST
def _call_bast(peptides: dict[str, tuple[str, str]], root: Path) -> TypingResult:
    """Bexsero antigen sequence type, with the MenDeVAR index."""
    values: dict[str, str] = {}
    notes: list[str] = []
    for locus in BAST_ORDER:
        identifier, note = peptides.get(locus, ("", "not determined"))
        values[locus] = identifier or "not determined"
        if note:
            notes.append(f"{locus}: {note}")
    alleles = {
        "fHbp_peptide": values["fHbp_peptide"],
        "NHBA_peptide": values["NHBA_peptide"],
        "NadA_peptide": values["NadA_peptide"],
        "PorA_VR1": values["PorA_VR1"],
        "PorA_VR2": values["PorA_VR2"],
    }
    extended = _extended_annotation(root, "fHbp_peptide", values["fHbp_peptide"])
    if extended:
        family = extended.get("Novartis_variant_family", "")
        subfamily = extended.get("Pfizer_subfamily", "")
        alleles["fHbp_variant_group"] = (
            f"variant group {family}" if family else "not assigned") + (
            f", subfamily {subfamily}" if subfamily else "")
        if extended.get("Novartis_variant"):
            alleles["fHbp_variant"] = extended["Novartis_variant"]
        if extended.get("Pfizer_variant"):
            alleles["fHbp_Pfizer_variant"] = extended["Pfizer_variant"]
    nada = _extended_annotation(root, "NadA_peptide", values["NadA_peptide"])
    if nada.get("NadA_variant"):
        alleles["NadA_variant"] = nada["NadA_variant"]

    bast, reactivity = _bast_profile(root, values)
    alleles.update(reactivity)
    alleles["vaccine_peptides"] = "; ".join(f"{name}: {contents}"
                                            for name, contents in VACCINE_PEPTIDES.items())
    call = f"BAST-{bast}" if bast else "not in the profile table"
    if not bast:
        notes.append("this combination of the five peptides is not in the downloaded BAST "
                     "profile table, so it has no BAST number; the component peptides "
                     "above are the report, and a new combination is assigned a number "
                     "only by PubMLST")
    notes.append("the BAST components are the five peptides contained in Bexsero: fHbp, "
                 "NHBA, NadA and PorA VR1 and VR2. PorA VR1 is in the scheme for extra "
                 "discrimination although it is not part of the MATS assay")
    found = sum(1 for locus in BAST_ORDER
                if values[locus] not in ("", "not determined", "new"))
    return TypingResult(scheme="bast", call=call,
                        lineage=reactivity.get("MenDeVAR_Bexsero", "-"),
                        alleles=alleles, loci_found=found, loci_total=len(BAST_ORDER),
                        note=". ".join(notes))


def _extended_annotation(root: Path, locus: str, allele: str) -> dict[str, str]:
    """Look one peptide up in a saved extended-attribute table."""
    path = root / "peptides" / f"{locus}.extended.tsv"
    if not path.exists() or not allele or not allele.isdigit():
        return {}
    with open(path) as handle:
        fields = handle.readline().rstrip("\n").split("\t")[1:]
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if row and row[0] == allele:
                return {name: value for name, value in zip(fields, row[1:]) if value}
    return {}


def _bast_profile(root: Path, values: dict[str, str]) -> tuple[str, dict[str, str]]:
    """Look the five peptides up in the BAST profile table."""
    path = root / "peptides" / "bast_profiles.txt"
    if not path.exists():
        return "", {}
    wanted = tuple(values[locus] for locus in BAST_ORDER)
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            columns = [header.index(locus) for locus in BAST_ORDER]
        except ValueError:
            return "", {}
        bexsero = (header.index("MenDeVAR_Bexsero_reactivity")
                   if "MenDeVAR_Bexsero_reactivity" in header else None)
        trumenba = (header.index("MenDeVAR_Trumenba_reactivity")
                    if "MenDeVAR_Trumenba_reactivity" in header else None)
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(row) <= max(columns):
                continue
            if tuple(row[index].strip() for index in columns) != wanted:
                continue
            out: dict[str, str] = {}
            if bexsero is not None and len(row) > bexsero:
                out["MenDeVAR_Bexsero"] = _mendevar(row[bexsero].strip(), "Bexsero")
            if trumenba is not None and len(row) > trumenba:
                out["MenDeVAR_Trumenba"] = _mendevar(row[trumenba].strip(), "Trumenba")
            return row[0].strip(), out
    return "", {}


#: What each MenDeVAR category means. The index is deduced from the peptides,
#: not measured on the isolate, and says nothing about an individual's response.
MENDEVAR_MEANING = {
    "exact match": "at least one antigen is an exact peptide match to a vaccine "
                   "component, so the isolate is expected to be covered",
    "cross-reactive": "no exact match, but at least one antigen has been shown "
                      "experimentally to cross-react; coverage is likely but not certain",
    "none": "no exact match and every antigen tested has been shown experimentally not to "
            "cross-react, so the isolate is not expected to be covered",
    "insufficient data": "the antigens have not been tested, or were tested in too few "
                         "isolates, so no prediction can be made; this is the commonest "
                         "category and is not evidence of either coverage or escape",
}


def _mendevar(value: str, vaccine: str) -> str:
    meaning = MENDEVAR_MEANING.get(value.lower(), "")
    text = value or "not determined"
    return f"{text} ({vaccine}): {meaning}" if meaning else f"{text} ({vaccine})"


# --------------------------------------------------------------------- MLST
def _call_mlst(assembly: Path, root: Path, work: Path, threads: int,
               mlst=None) -> TypingResult:
    """Seven-locus sequence type and clonal complex."""
    allele_db = root / "mlst" / "mlst.fna"
    profiles = root / "mlst" / "profiles.txt"
    # Hydra's own MLST call first: it types these loci already.
    profile = common.profile_from_call(mlst, "neisseria", profiles,
                                       accepted=MLST_SCHEMES)
    if profile is None:
        if not allele_db.exists():
            return TypingResult(
                scheme="mlst", call="not determined",
                note=f"MLST data missing; run 'hydra db download typing-{KEY}', "
                     f"or install Hydra's pubmlst database so the loci are typed "
                     f"once in the MLST step")
        profile = common.call_mlst(assembly, allele_db, profiles, "neisseria",
                                   workdir=work, threads=threads)
    note = profile.note
    if profile.sequence_type != "-" and not profile.clonal_complex:
        note = "; ".join(part for part in [note, (
            f"ST{profile.sequence_type} is not assigned to a clonal complex in the PubMLST "
            f"profile table, which leaves about 39% of sequence types unassigned; the "
            f"designation is then written without the complex")] if part)
    elif profile.clonal_complex:
        note = "; ".join(part for part in [note, (
            f"the clonal complex is the hyperinvasive lineage the isolate belongs to and "
            f"is what links apparently unrelated cases; reported as "
            f"{_short_complex(profile.clonal_complex)} in the standard designation")] if part)
    call = f"ST-{profile.sequence_type}" if profile.sequence_type != "-" else "not determined"
    return TypingResult(scheme="mlst", call=call,
                        lineage=profile.clonal_complex or "-",
                        alleles=dict(profile.alleles),
                        loci_found=profile.loci_found, loci_total=profile.loci_total,
                        note=note)


# ---------------------------------------------------------------------- AMR
def _call_amr(assembly: Path, root: Path, work: Path, threads: int) -> TypingResult:
    """Report the *penA*, *gyrA* and *rpoB* allele identifiers.

    The allele number is the unit meningococcal resistance surveillance is
    reported in: penA alleles with an altered PBP2 give reduced penicillin
    susceptibility, gyrA alleles carrying a quinolone-resistance substitution
    give ciprofloxacin resistance, and rpoB alleles give rifampicin resistance.
    This module reports which allele the isolate carries, and does not convert
    the allele into a phenotype: that mapping is not published as a table this
    module could implement, and inventing one would be worse than naming the
    allele and saying so.
    """
    allele_db = root / "amr" / "amr.fna"
    if not allele_db.exists():
        return TypingResult(scheme="amr_alleles", call="not determined",
                            note=f"allele data missing; run 'hydra db download {KEY}'")
    hits = common.screen_nucleotide(assembly, allele_db, work, threads=threads,
                                    min_identity=95.0, min_coverage=90.0,
                                    tag="amr", task="megablast")
    best: dict[str, tuple[float, float, str, str]] = {}
    for hit in hits:
        name = hit.gene
        if "_" not in name:
            continue
        locus, allele = name.rsplit("_", 1)
        if locus not in AMR_LOCI:
            continue
        key = (hit.coverage_pct, hit.identity_pct)
        current = best.get(locus)
        if current is None or key > (current[0], current[1]):
            best[locus] = (hit.coverage_pct, hit.identity_pct, allele, hit.coverage)
    alleles: dict[str, str] = {}
    called: list[str] = []
    for locus in AMR_LOCI:
        found = best.get(locus)
        if found is None:
            alleles[locus] = ("not determined: no allele reached 95% identity over 90% of "
                              "its length")
            continue
        coverage, identity, allele, span = found
        if coverage >= 100.0 and identity >= 100.0:
            alleles[locus] = f"allele {allele} (exact match)"
            called.append(f"{locus} {allele}")
        else:
            alleles[locus] = (f"closest allele {allele} at {identity:.2f}% identity over "
                              f"{coverage:.0f}% of its length ({span}); not an exact match, "
                              f"so this is a new allele and PubMLST would assign it a new "
                              f"number")
            called.append(f"{locus} ~{allele}")
    return TypingResult(
        scheme="amr_alleles", call=", ".join(called) if called else "not determined",
        lineage="-", alleles=alleles, loci_found=len(best), loci_total=len(AMR_LOCI),
        note="allele identifiers only. An exact match names the PubMLST allele; anything "
             "less is a new allele and is reported as the closest one with its distance. "
             "Converting an allele into a susceptibility prediction needs the published "
             "allele-to-phenotype tables and is not attempted here, so this result "
             "supplements phenotypic testing rather than replacing it")
