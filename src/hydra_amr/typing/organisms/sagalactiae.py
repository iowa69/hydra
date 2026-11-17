"""*Streptococcus agalactiae* (group B streptococcus) typing.

What this module reports, in the order a clinical report wants it:

1. **Capsular serotype** (Ia, Ib, II-IX, or non-typeable) read from the *cps*
   locus, by the published GBS-SBG rules: the reference serotype whose
   serotype-determining *cps* region the assembly covers best, with the
   uncertainty flags the reference tool defines.
2. **Surface proteins and pilus islands**: the alpha-like protein family, the
   serine-rich repeat adhesins *srr1* and *srr2*, and the pilus-island genes.
3. **The hypervirulent-clone markers** and the rest of the repertoire used in
   neonatal-disease and vaccine reporting.

Two honesty notes that the calls depend on, stated here rather than buried:

**A genotypic serotype is not a serological serotype.** GBS-SBG assigns the
capsular genotype encoded in *cps*. An isolate that is non-typeable by
serology usually carries an intact, typeable *cps* locus and will be given a
serotype here; the discrepancy is in capsule expression, not in the genotype.
The module says which it is reporting.

**The alpha-like protein family is reported as a family, not by member name.**
*bca* (alpha-C) and *bac* (beta-C) are separately annotated and are called by
name. *rib*, *alp1*/epsilon, *alp2*, *alp3* and *alp4* share their N-terminal
domain and their repeat units, RefSeq annotates all of them as "alpha-like
surface protein" without distinguishing them, and no versioned, downloadable
allele set for them exists. Naming one of them from a best BLAST hit would be
guessing, so the module reports that an Alp-family gene is present, gives the
closest reference and the identity, and says the member is not determined.

See ``docs/typing/strep.md`` for the rules, thresholds and bibliography.
"""

from __future__ import annotations

from pathlib import Path

from ...engines.blast import blast, check_db_exists
from ...records import TypingResult
from ...seqio import read_fasta
from ...utils import LOG, HydraError, tempdir
from .base import OrganismDb, Reference
from .strep_common import (
    BLASTX_THREADS,
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

KEY = "sagalactiae"
TITLE = "Streptococcus agalactiae"
SPECIES: tuple[str, ...] = ("Streptococcus agalactiae",)
MLST_SCHEMES: tuple[str, ...] = ("sagalactiae",)

#: The GBS-SBG curated capsular reference file: one sequence per serotype
#: covering the serotype-determining region of the *cps* locus, plus four
#: variants of serotype III.
GBS_SBG_URL = ("https://raw.githubusercontent.com/swainechen/GBS-SBG/main/"
               "GBS-SBG.fasta")
GBS_SBG_LICENCE_URL = ("https://raw.githubusercontent.com/swainechen/GBS-SBG/"
                       "main/LICENSE")

# ------------------------------------------------------- GBS-SBG calling rules
#: "The default cutoffs correspond to those used in SRST2 - at least 90%
#: identity over at least 90% of the total length of the reference sequence is
#: required, and the highest identity and coverage determine the best call."
SBG_MIN_IDENTITY = 90.0
SBG_MIN_COVERAGE = 0.90
#: Per-HSP filters GBS-SBG applies before anything is summed. Both are strict
#: inequalities in the reference implementation, and both matter: without the
#: length filter a long contig throws dozens of 15-30 bp matches into the bit
#: score sum and the ranking stops meaning anything.
SBG_HSP_MIN_IDENTITY = 90.0
SBG_HSP_MIN_LENGTH = 100
#: The uncertainty flags GBS-SBG raises, with the tool's own thresholds.
SBG_FLAG_COVERAGE = 0.99
SBG_FLAG_IDENTITY = 99.0
SBG_FLAG_MAX_HSPS = 1
SBG_FLAG_MAX_CONTIGS = 1
#: A runner-up whose summed bit score reaches this share of the winner's makes
#: the call uncertain, and at 0.90 the runner-up is reported alongside it.
SBG_FLAG_NEXT_RATIO = 0.95
SBG_REPORT_NEXT_RATIO = 0.90

PANEL_SOURCES: dict[str, PanelSource] = {
    "NEM316": PanelSource("GCF_000196055.1", "ASM19605v1", "NEM316",
                          "serotype III, ST23; Glaser 2002"),
    "A909": PanelSource("GCF_000012705.1", "ASM1270v1", "A909",
                        "serotype Ia, ST7; Tettelin 2005"),
    "2603VR": PanelSource("GCF_000007265.1", "ASM726v1", "2603V/R",
                          "serotype V, ST110; Tettelin 2002"),
    "COH1": PanelSource("GCF_030489705.1", "ASM3048970v1", "COH1",
                        "serotype III, CC17; the hypervirulent neonatal clone"),
}


def _g(name: str, symbols: tuple[str, ...], source: str, product: str, aa: int,
       meaning: str, *, location: str = "surface", product_re: str = "",
       ident: float = 80.0, cover: float = 60.0,
       reference: str = "Rosa-Fraile2014") -> GeneSpec:
    return GeneSpec(name=name, symbols=symbols, source=source, product=product,
                    expect_aa=aa, meaning=meaning, location=location,
                    product_re=product_re, min_identity=ident,
                    min_coverage=cover, reference=reference)


#: The alpha-like protein family. Three reference sequences of different lengths
#: are carried so that the closest can be reported, but none of them is used to
#: name a family member - see the module docstring.
ALP_FAMILY: tuple[GeneSpec, ...] = (
    _g("bca", ("bca",), "A909", "alpha-C protein Bca", 1020,
       "the alpha-C protein; one of the two components of the GBS-NN/NN2 "
       "protein vaccine candidate, and a protective antigen",
       ident=70.0, cover=35.0, reference="Michel1992"),
    _g("bac", ("bac",), "A909", "beta-C protein Bac (IgA-binding)", 1164,
       "the beta-C protein, which binds IgA Fc and factor H; carried with bca "
       "in the classic c-protein-positive serotype Ia/Ib strains",
       ident=70.0, cover=35.0, reference="Jerlstrom1991"),
    _g("alp_family_long", (), "2603VR",
       "alpha-like surface protein (long reference)", 1389,
       "an alpha-like protein family member; the family includes rib, alp1 "
       "(epsilon), alp2, alp3 and alp4, which share their N-terminal domain and "
       "repeat units and are not distinguished here",
       product_re=r"^alpha-like surface protein$", ident=70.0, cover=25.0,
       reference="Creti2004"),
    _g("alp_family_mid", (), "NEM316",
       "alpha-like surface protein (mid-length reference)", 1126,
       "an alpha-like protein family member (the NEM316 allele, which the "
       "literature calls rib)",
       product_re=r"^alpha-like surface protein$", ident=70.0, cover=25.0,
       reference="Creti2004"),
    _g("alp_family_short", (), "COH1",
       "alpha-like surface protein (short reference)", 757,
       "an alpha-like protein family member (the COH1/CC17 allele)",
       product_re=r"^alpha-like surface protein$", ident=70.0, cover=25.0,
       reference="Creti2004"),
)

#: The pilus islands. PI-2a is named because RefSeq annotates its genes as such.
#: The second island locus is detected but not named - see ``_call_pilus``.
PILUS: tuple[GeneSpec, ...] = (
    _g("pilA", ("pila",), "NEM316", "PI-2a pilus tip adhesin PilA", 901,
       "the PI-2a ancillary protein 1 / tip adhesin; a pilus-based vaccine "
       "candidate", ident=75.0, cover=50.0, reference="Margarit2009"),
    _g("pilB", ("pilb",), "NEM316", "PI-2a pilus backbone PilB", 674,
       "the PI-2a backbone protein", ident=75.0, cover=50.0,
       reference="Margarit2009"),
    _g("pilC", ("pilc",), "NEM316", "PI-2a pilus ancillary protein PilC", 308,
       "the PI-2a ancillary protein 2", ident=75.0, cover=50.0,
       reference="Margarit2009"),
    _g("srtC4", ("srtc4",), "NEM316", "PI-2a pilus assembly sortase SrtC4", 294,
       "the PI-2a class C sortase that polymerises the backbone",
       ident=75.0, cover=50.0, reference="Margarit2009"),
    _g("pilus_backbone_2", (), "NEM316",
       "SpaH/EbpB family pilus backbone of the second pilus island locus", 554,
       "the backbone of the second pilus island; present in PI-1 and in PI-2b, "
       "which this module does not distinguish",
       product_re=r"^SpaH/EbpB family LPXTG-anchored major pilin$",
       ident=70.0, cover=45.0, reference="Rosini2006"),
    _g("pilus_backbone_3", (), "A909",
       "SpaH/EbpB family pilus backbone, second allele", 502,
       "the second allele of the SpaH/EbpB pilus backbone",
       product_re=r"^SpaH/EbpB family LPXTG-anchored major pilin$",
       ident=70.0, cover=45.0, reference="Rosini2006"),
    _g("srtC1", ("srtc1",), "NEM316", "class C pilus sortase SrtC1", 305,
       "class C sortase of the second pilus island locus", ident=75.0,
       cover=50.0, reference="Rosini2006"),
    _g("srtC2", ("srtc2",), "NEM316", "class C pilus sortase SrtC2", 283,
       "class C sortase of the second pilus island locus", ident=75.0,
       cover=50.0, reference="Rosini2006"),
)

#: Virulence and clone markers reported alongside the serotype.
MARKERS: tuple[GeneSpec, ...] = (
    _g("srr1", ("srr1",), "NEM316",
       "serine-rich repeat glycoprotein adhesin Srr1", 1310,
       "the serine-rich repeat adhesin of most lineages; binds fibrinogen and "
       "keratin 4", ident=75.0, cover=30.0, reference="Six2015"),
    _g("srr2", (), "COH1", "serine-rich repeat glycoprotein adhesin Srr2", 1085,
       "Srr2 replaces Srr1 in the hypervirulent CC17 lineage, so its presence "
       "is a CC17 marker and is associated with neonatal meningitis",
       product_re=r"serine-rich repeat glycoprotein adhesin Srr2",
       ident=75.0, cover=30.0, reference="Six2015"),
    _g("hvgA", (), "COH1", "hypervirulent GBS adhesin HvgA (Gbs2018C)", 522,
       "the CC17 hypervirulence and meningeal-tropism marker. HvgA is the C "
       "variant of the gbs2018 locus; the A variant at the same locus is BibA, "
       "which every group B streptococcus carries. The two share their N- and "
       "C-terminal domains and are only about 38% identical overall, so the "
       "call rests on covering the divergent central core, and the coverage "
       "floor is set high for that reason",
       # The threshold sits in the middle of a wide, measured gap. On this
       # module's validation panel the two hvgA-positive CC17 genomes match the
       # reference at 85.4% and 100.0% amino-acid identity over the full 522
       # residues, while the bibA-carrying genomes reach only 57.9% and 60.8%
       # over 82-85% of it - the shared N- and C-terminal domains and nothing
       # of the divergent core.
       product_re=r"^MucBP domain-containing protein$", ident=75.0, cover=80.0,
       reference="Tazi2010"),
    _g("lmb", ("lmb",), "NEM316", "laminin-binding protein Lmb", 306,
       "laminin-binding adhesin; part of the lmb-scpB region used in "
       "surveillance", ident=85.0, cover=70.0, reference="Rosa-Fraile2014"),
    _g("scpB", ("scpa/b", "scpb"), "NEM316", "C5a peptidase ScpB", 1150,
       "cleaves C5a and so blunts neutrophil recruitment; a near-universal GBS "
       "gene, so its absence is the notable finding", ident=85.0, cover=70.0,
       reference="Rosa-Fraile2014"),
    _g("cfb", ("cfb",), "NEM316", "CAMP factor Cfb", 255,
       "CAMP factor, the basis of the CAMP test used to identify GBS",
       ident=80.0, cover=60.0, reference="Rosa-Fraile2014"),
    _g("cylE", ("cyle",), "NEM316", "beta-haemolysin/cytolysin CylE", 0,
       "the beta-haemolysin/cytolysin, which gives GBS its haemolysis and its "
       "orange pigment", product_re=r"CylE|beta-h[ae]molysin", ident=80.0,
       cover=60.0, reference="Rosa-Fraile2014"),
    _g("sip", ("sip",), "NEM316", "surface immunogenic protein Sip", 434,
       "a conserved protective antigen present in all serotypes; a "
       "serotype-independent vaccine candidate", ident=85.0, cover=70.0,
       reference="Brodeur2000"),
    _g("fbsA", ("fbsa",), "NEM316", "fibrinogen-binding protein FbsA", 410,
       "fibrinogen-binding adhesin; repeat-variable, so coverage below 100% is "
       "expected", ident=70.0, cover=30.0, reference="Rosa-Fraile2014"),
    _g("fbsB", ("fbsb",), "NEM316", "fibrinogen-binding protein FbsB", 635,
       "the second fibrinogen-binding adhesin", ident=75.0, cover=50.0,
       reference="Rosa-Fraile2014"),
    _g("pbp2x", ("pbp2x", "pbp2X"), "NEM316", "penicillin-binding protein 2X",
       752, "the primary target of penicillin in GBS; substitutions in the "
            "transpeptidase domain are the mechanism of reduced penicillin "
            "susceptibility. This module reports whether the gene differs from "
            "the wild-type reference and by how much; it does not call named "
            "substitutions, because that needs a verified residue catalogue",
       location="chromosome", ident=90.0, cover=90.0, reference="Kimura2008"),
    _g("neuA", ("neua",), "NEM316", "N-acylneuraminate cytidylyltransferase",
       413, "sialic acid activation; part of the cps locus and required for "
            "the sialylated capsule", location="chromosome", ident=85.0,
       cover=70.0, reference="Cieslewicz2005"),
    _g("cpsA", ("cpsa",), "NEM316", "capsule regulator CpsA/LytR", 485,
       "the first gene of the cps operon; its presence shows the locus is "
       "there even when the serotype-determining region is not resolved",
       location="chromosome", ident=85.0, cover=70.0,
       reference="Cieslewicz2005"),
)

PANEL: tuple[GeneSpec, ...] = ALP_FAMILY + PILUS + MARKERS

#: Serotypes the hexavalent capsular conjugate vaccine candidate covers.
HEXAVALENT = ("Ia", "Ib", "II", "III", "IV", "V")

#: The PBP2X reference the substitution positions are numbered against. Pinning
#: it is not optional: the positions in the Japanese PRGBS literature are given
#: against this protein, and it already carries V at 377, so an isolate with V377
#: is wild type here while those papers call V377 a substitution. NCBI's curated
#: AMRFinderPlus mutation catalogue uses this same accession for its single
#: S. agalactiae entry, pbp2x_Q557E.
PBP2X_REFERENCE = "WP_017647104.1"
PBP2X_REFERENCE_LEN = 752

#: Transpeptidase-domain positions reported in reduced-penicillin-susceptible
#: group B streptococcus, with the wild-type residue at each position in
#: PBP2X_REFERENCE. Only substitutions at these catalogued positions are named;
#: anything else is reported as an uncatalogued difference.
PBP2X_POSITIONS: dict[int, tuple[str, str]] = {
    329: ("G", "G329V reported with the V405A background in serotype Ia/ST1 "
               "isolates at penicillin MIC 0.25 mg/L"),
    377: ("V", "the reference already carries V here; the literature's I377V "
               "is therefore wild type against this accession"),
    394: ("T", "T394A, a transpeptidase-domain substitution reported in PRGBS"),
    395: ("F", "F395L, reported with V405A in serotype III/ST1 isolates"),
    398: ("G", "G398A, reported with V405A and with Q557E"),
    405: ("V", "V405A, immediately C-terminal to the SSN motif (402-404) and "
               "the commonest PRGBS substitution; on its own it accompanies a "
               "penicillin MIC of about 0.25 mg/L"),
    412: ("Q", "Q412L, reported in a serotype III/ST464 isolate"),
    429: ("G", "G429D, reported with V405A"),
    433: ("R", "R433H, reported with V405A"),
    438: ("H", "H438Y, reported in two PRGBS lineages"),
    510: ("V", "V510I, reported with V405A"),
    557: ("Q", "Q557E, three residues C-terminal to the KSG motif (552-554); the "
               "only substitution NCBI curates for this species, and the one "
               "that escalates the penicillin MIC to about 1 mg/L when it "
               "accompanies V405A"),
    648: ("G", "G648A, reported with V405A"),
}
#: The two substitutions with direct MIC evidence behind them. V405A alone
#: accompanies a penicillin MIC of about 0.25 mg/L and V405A with Q557E about
#: 1 mg/L; Q557E is the only one NCBI curates. Everything else in
#: PBP2X_POSITIONS is reported as an accompanying substitution.
PBP2X_DETERMINANTS: frozenset[int] = frozenset({405, 557})

#: The catalytic and substrate-binding motifs, so a lesion in one can be named.
PBP2X_MOTIFS: tuple[tuple[int, str, str], ...] = (
    (344, "STMK", "the SXXK catalytic motif; Ser344 is the nucleophile"),
    (402, "SSN", "the SSN substrate-binding motif"),
    (552, "KSG", "the KSG substrate-binding motif, three residues N-terminal "
                 "to position 557"),
)

DATABASES: tuple[OrganismDb, ...] = (
    OrganismDb(
        name="cps",
        source=GBS_SBG_URL,
        version="the repository has no release tags; fetch() records the commit "
                "SHA of main at retrieval time in cps.version",
        licence="MIT License, Copyright (c) 2021 Swaine Chen: \"Permission is "
                "hereby granted, free of charge, to any person obtaining a copy "
                "of this software and associated documentation files ... to "
                "deal in the Software without restriction\"",
        redistribution="permitted",
        reference="Tiruvayipati2021",
        notes="13 sequences, about 90 kB: one per serotype Ia, Ib, II, IV-IX "
              "and four variants of serotype III, each covering the "
              "serotype-determining region of the cps locus.",
    ),
    OrganismDb(
        name="panel",
        source="https://ftp.ncbi.nlm.nih.gov/genomes/all (RefSeq feature tables "
               "and protein FASTA for "
               + ", ".join(s.accession for s in PANEL_SOURCES.values()) + ")",
        version="pinned assembly accession.version per strain, see PANEL_SOURCES",
        licence="NCBI: \"NLM does not claim copyright on the data in the "
                "Database\"; RefSeq records are free to use and redistribute.",
        redistribution="permitted",
        reference="OLeary2016",
        notes="Reference proteins for the surface-protein, pilus and marker "
              "panel, selected from RefSeq annotation by gene symbol or product "
              "name. The built database is about 30 kB; the download is about "
              "6 MB and is cached under db_dir/sagalactiae/refseq/.",
    ),
)

REFERENCES: tuple[Reference, ...] = (
    Reference(
        key="Tiruvayipati2021",
        authors="Tiruvayipati S, Tang WY, Barkham TMS, Chen SL",
        year=2021,
        title="GBS-SBG - GBS Serotyping by Genome Sequencing",
        journal="Microbial Genomics 7:000688",
        doi="10.1099/mgen.0.000688",
        pmid="34895403",
        supports="the capsular serotyping database, the serotype-determining "
                 "cps region it contains, and the calling and uncertainty rules "
                 "this module implements",
    ),
    Reference(
        key="Inouye2014",
        authors="Inouye M, Dashnow H, Raven LA, et al.",
        year=2014,
        title="SRST2: Rapid genomic surveillance for public health and hospital "
              "microbiology labs",
        journal="Genome Medicine 6:90",
        doi="10.1186/s13073-014-0090-6",
        pmid="25422674",
        supports="the 90% identity over 90% of the reference length cut-off that "
                 "GBS-SBG inherits and this module applies",
    ),
    Reference(
        key="Cieslewicz2005",
        authors="Cieslewicz MJ, Chaffin D, Glusman G, et al.",
        year=2005,
        title="Structural and genetic diversity of group B Streptococcus "
              "capsular polysaccharides",
        journal="Infection and Immunity 73:3096-3103",
        doi="10.1128/IAI.73.5.3096-3103.2005",
        pmid="15845517",
        supports="the structure of the cps locus, which genes determine the "
                 "serotype, and the molecular basis of the serotypes Ia to IX",
    ),
    Reference(
        key="Slotved2007",
        authors="Slotved HC, Kong F, Lambertsen L, et al.",
        year=2007,
        title="Serotype IX, a proposed new Streptococcus agalactiae serotype",
        journal="Journal of Clinical Microbiology 45:2929-2936",
        doi="10.1128/JCM.00117-07",
        pmid="17634306",
        supports="serotype IX, the most recently described serotype in the "
                 "database",
    ),
    Reference(
        key="Glaser2002",
        authors="Glaser P, Rusniok C, Buchrieser C, et al.",
        year=2002,
        title="Genome sequence of Streptococcus agalactiae, a pathogen causing "
              "invasive neonatal disease",
        journal="Molecular Microbiology 45:1499-1513",
        doi="10.1046/j.1365-2958.2002.03126.x",
        pmid="12354221",
        supports="NEM316 (serotype III, ST23), the source of most reference "
                 "proteins in the panel",
    ),
    Reference(
        key="Tettelin2002",
        authors="Tettelin H, Masignani V, Cieslewicz MJ, et al.",
        year=2002,
        title="Complete genome sequence and comparative genomic analysis of an "
              "emerging human pathogen, serotype V Streptococcus agalactiae",
        journal="Proceedings of the National Academy of Sciences USA "
                "99:12391-12396",
        doi="10.1073/pnas.182380799",
        pmid="12200547",
        supports="2603V/R (serotype V), a panel source and a validation genome",
    ),
    Reference(
        key="Tettelin2005",
        authors="Tettelin H, Masignani V, Cieslewicz MJ, et al.",
        year=2005,
        title="Genome analysis of multiple pathogenic isolates of Streptococcus "
              "agalactiae: implications for the microbial pan-genome",
        journal="Proceedings of the National Academy of Sciences USA "
                "102:13950-13955",
        doi="10.1073/pnas.0506758102",
        pmid="16172379",
        supports="A909 (serotype Ia) and the other reference strains, and the "
                 "serotypes published for them",
    ),
    Reference(
        key="DaCunha2014",
        authors="Da Cunha V, Davies MR, Douarre PE, et al.",
        year=2014,
        title="Streptococcus agalactiae clones infecting humans were selected "
              "and fixed through the extensive use of tetracycline",
        journal="Nature Communications 5:4544",
        doi="10.1038/ncomms5544",
        pmid="25088811",
        supports="the global GBS collection the validation genomes with "
                 "published serotypes are drawn from, and the clonal-complex "
                 "structure of human GBS",
    ),
    Reference(
        key="Rosini2006",
        authors="Rosini R, Rinaudo CD, Soriani M, et al.",
        year=2006,
        title="Identification of novel genomic islands coding for antigenic "
              "pilus-like structures in Streptococcus agalactiae",
        journal="Molecular Microbiology 61:126-141",
        doi="10.1111/j.1365-2958.2006.05225.x",
        pmid="16824100",
        supports="the pilus islands and their gene content, which is why the "
                 "second island locus is detected but not named PI-1 or PI-2b",
    ),
    Reference(
        key="Margarit2009",
        authors="Margarit I, Rinaudo CD, Galeotti CL, et al.",
        year=2009,
        title="Preventing bacterial infections with pilus-based vaccines: the "
              "group B Streptococcus paradigm",
        journal="The Journal of Infectious Diseases 199:108-115",
        doi="10.1086/595564",
        pmid="19086816",
        supports="the PI-2a gene set (pilA, pilB, pilC), the distribution of the "
                 "islands across lineages, and their use as vaccine candidates",
    ),
    Reference(
        key="Tazi2010",
        authors="Tazi A, Disson O, Bellais S, et al.",
        year=2010,
        title="The surface protein HvgA mediates group B streptococcus "
              "hypervirulence and meningeal tropism in neonates",
        journal="The Journal of Experimental Medicine 207:2313-2322",
        doi="10.1084/jem.20092594",
        pmid="20956545",
        supports="hvgA as the CC17 hypervirulence and meningeal-tropism marker, "
                 "and therefore why the module states that it cannot call it",
    ),
    Reference(
        key="Six2015",
        authors="Six A, Bellais S, Bouaboud A, et al.",
        year=2015,
        title="Srr2, a multifaceted adhesin expressed by ST-17 hypervirulent "
              "group B Streptococcus involved in binding to both fibrinogen and "
              "plasminogen",
        journal="Molecular Microbiology 97:1209-1222",
        doi="10.1111/mmi.13097",
        pmid="26094503",
        supports="srr2 as the ST-17/CC17-restricted adhesin, which is the marker "
                 "this module reports in place of hvgA",
    ),
    Reference(
        key="Michel1992",
        authors="Michel JL, Madoff LC, Olson K, et al.",
        year=1992,
        title="Large, identical, tandem repeating units in the C protein alpha "
              "antigen gene, bca, of group B streptococci",
        journal="Proceedings of the National Academy of Sciences USA "
                "89:10060-10064",
        doi="10.1073/pnas.89.21.10060",
        pmid="1438195",
        supports="bca, the alpha-C protein, and its tandem repeat structure, "
                 "which is why reference coverage below 100% is expected",
    ),
    Reference(
        key="Jerlstrom1991",
        authors="Jerlstrom PG, Chhatwal GS, Timmis KN",
        year=1991,
        title="The IgA-binding beta antigen of the c protein complex of Group B "
              "streptococci: sequence determination of its gene and detection of "
              "two binding regions",
        journal="Molecular Microbiology 5:843-849",
        doi="10.1111/j.1365-2958.1991.tb00757.x",
        pmid="1857207",
        supports="bac, the beta-C protein, and its IgA-binding function",
    ),
    Reference(
        key="Creti2004",
        authors="Creti R, Fabretti F, Orefici G, von Hunolstein C",
        year=2004,
        title="Multiplex PCR assay for direct identification of group B "
              "streptococcal alpha-protein-like protein genes",
        journal="Journal of Clinical Microbiology 42:1326-1329",
        doi="10.1128/JCM.42.3.1326-1329.2004",
        pmid="15004110",
        supports="the alpha-like protein family (bca, rib, alp1/epsilon, alp2, "
                 "alp3, alp4) and the sequence similarity between its members "
                 "that stops this module naming one of them",
    ),
    Reference(
        key="Brodeur2000",
        authors="Brodeur BR, Boyer M, Charlebois I, et al.",
        year=2000,
        title="Identification of group B streptococcal Sip protein, which elicits "
              "cross-protective immunity",
        journal="Infection and Immunity 68:5610-5618",
        doi="10.1128/IAI.68.10.5610-5618.2000",
        pmid="10992461",
        supports="Sip as a conserved, serotype-independent protective antigen",
    ),
    Reference(
        key="Kimura2008",
        authors="Kimura K, Suzuki S, Wachino J, et al.",
        year=2008,
        title="First molecular characterization of group B streptococci with "
              "reduced penicillin susceptibility",
        journal="Antimicrobial Agents and Chemotherapy 52:2890-2897",
        doi="10.1128/AAC.00185-08",
        pmid="18490507",
        supports="PBP2X substitutions as the mechanism of reduced penicillin "
                 "susceptibility in GBS, and therefore why pbp2x divergence is "
                 "reported even though named substitutions are not called",
    ),
    Reference(
        key="Rosini2015",
        authors="Rosini R, Campisi E, De Chiara M, et al.",
        year=2015,
        title="Genomic analysis reveals the molecular basis for capsule loss in "
              "the group B Streptococcus population",
        journal="PLoS One 10:e0125985",
        doi="10.1371/journal.pone.0125985",
        pmid="25946017",
        supports="why a capsular genotype is not a serological serotype: 89 of "
                 "96 serologically non-typeable strains carry an inactivating "
                 "mutation in cps, most often a premature stop in cpsE, which a "
                 "coverage-and-identity rule cannot see. This is the basis of "
                 "the caveat printed with every cps call",
    ),
    Reference(
        key="Lamy2006",
        authors="Lamy MC, Dramsi S, Billoet A, et al.",
        year=2006,
        title="Rapid detection of the \"highly virulent\" group B Streptococcus "
              "ST-17 clone",
        journal="Microbes and Infection 8:1714-1722",
        doi="10.1016/j.micinf.2006.02.008",
        pmid="16822689",
        supports="the three variants (A, B, C) of the gbs2018 locus, of which "
                 "HvgA is variant C and BibA variant A; this is why the hvgA "
                 "call needs a high coverage floor over the divergent core",
    ),
    Reference(
        key="Martins2013",
        authors="Martins ER, Andreu A, Melo-Cristino J, Ramirez M",
        year=2013,
        title="Distribution of pilus islands in Streptococcus agalactiae that "
              "cause human infections: insights into evolution and implication "
              "for vaccine development",
        journal="Clinical and Vaccine Immunology 20:313-316",
        doi="10.1128/CVI.00529-12",
        pmid="23269415",
        supports="the pilus-island profiles seen in human disease and the "
                 "association of PI-2b with CC17, which the pilus call is "
                 "interpreted against",
    ),
    Reference(
        key="Koide2022",
        authors="Koide S, Nagano Y, Takizawa S, Sakaguchi K",
        year=2022,
        title="Genomic traits associated with virulence and antimicrobial "
              "resistance of invasive group B Streptococcus isolates with "
              "reduced penicillin susceptibility from elderly adults",
        journal="Microbiology Spectrum 10:e0056822",
        doi="10.1128/spectrum.00568-22",
        pmid="35639001",
        supports="the PBP2X substitution table the module names calls from, and "
                 "the observation that V405A alone accompanies a penicillin MIC "
                 "near 0.25 mg/L while V405A with Q557E reaches about 1 mg/L",
    ),
    Reference(
        key="Dahesh2008",
        authors="Dahesh S, Hensler ME, Van Sorge NM, et al.",
        year=2008,
        title="Point mutation in the group B streptococcal pbp2x gene conferring "
              "decreased susceptibility to beta-lactam antibiotics",
        journal="Antimicrobial Agents and Chemotherapy 52:2915-2918",
        doi="10.1128/AAC.00461-08",
        pmid="18541727",
        supports="that a single pbp2x point mutation is sufficient to lower "
                 "beta-lactam susceptibility in this species",
    ),
    Reference(
        key="Rosa-Fraile2014",
        authors="Rosa-Fraile M, Dramsi S, Spellerberg B",
        year=2014,
        title="Group B streptococcal haemolysin and pigment, a tale of twins",
        journal="FEMS Microbiology Reviews 38:932-946",
        doi="10.1111/1574-6976.12071",
        pmid="24617549",
        supports="the beta-haemolysin/cytolysin and the clinical meaning printed "
                 "for the virulence genes in the marker panel",
    ),
    Reference(
        key="Carreras-Abad2020",
        authors="Carreras-Abad C, Ramkhelawon L, Heath PT, Le Doare K",
        year=2020,
        title="A vaccine against group B Streptococcus: recent advances",
        journal="Infection and Drug Resistance 13:1263-1272",
        doi="10.2147/IDR.S203454",
        pmid="32425562",
        supports="the hexavalent capsular conjugate serotype set and the "
                 "protein-based vaccine candidates the module reports coverage "
                 "against",
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
    """Download and build this module's data into ``db_dir/sagalactiae/``."""
    root = Path(db_dir) / KEY
    root.mkdir(parents=True, exist_ok=True)
    cps = root / "cps.fna"
    if not cps.exists():
        LOG.info("%s: downloading the GBS-SBG capsular reference", KEY)
        download(GBS_SBG_URL, cps)
        download(GBS_SBG_LICENCE_URL, root / "cps.LICENSE")
        (root / "cps.version").write_text(
            f"{_sbg_commit()}\t{GBS_SBG_URL}\n", encoding="utf-8")
    make_blast_db(cps, "nucl")

    pbp = root / "pbp2x.faa"
    if not pbp.exists():
        LOG.info("%s: downloading the PBP2X reference %s", KEY, PBP2X_REFERENCE)
        download("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                 f"?db=protein&id={PBP2X_REFERENCE}&rettype=fasta&retmode=text",
                 pbp)
    _check_pbp2x(pbp)
    make_blast_db(pbp, "prot")

    LOG.info("%s: building the surface-protein panel from RefSeq annotation", KEY)
    fetch_panel(PANEL, PANEL_SOURCES, root, name="panel")


def _check_pbp2x(path: Path) -> None:
    """Verify the downloaded reference is the protein the positions assume.

    Numbering is the whole substance of a substitution call, so the wild-type
    residue at every catalogued position is checked against the sequence that
    was actually downloaded. If upstream ever serves a different protein under
    this accession the module fails here rather than reporting substitutions
    against the wrong numbering.
    """
    records = list(read_fasta(path))
    if not records:
        raise HydraError(f"{path} holds no sequence")
    seq = records[0][1]
    if len(seq) != PBP2X_REFERENCE_LEN:
        raise HydraError(
            f"{PBP2X_REFERENCE} is {len(seq)} aa, not the "
            f"{PBP2X_REFERENCE_LEN} aa the substitution positions assume")
    wrong = [f"{pos}: expected {wt}, found {seq[pos - 1]}"
             for pos, (wt, _why) in sorted(PBP2X_POSITIONS.items())
             if seq[pos - 1] != wt]
    for start, motif, _why in PBP2X_MOTIFS:
        if seq[start - 1:start - 1 + len(motif)] != motif:
            wrong.append(f"{start}-{start + len(motif) - 1}: expected motif "
                         f"{motif}, found {seq[start - 1:start - 1 + len(motif)]}")
    if wrong:
        raise HydraError(f"{PBP2X_REFERENCE} does not match the catalogued "
                         f"numbering: " + "; ".join(wrong))


def _sbg_commit() -> str:
    """The commit the capsular reference was taken at, since there is no tag."""
    import json  # noqa: PLC0415 - only fetch() needs it
    import urllib.request  # noqa: PLC0415

    url = "https://api.github.com/repos/swainechen/GBS-SBG/commits/main"
    try:
        request = urllib.request.Request(
            url, headers={"User-Agent": "hydra-amr",
                          "Accept": "application/vnd.github+json"})
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = json.load(response)
        return str(payload.get("sha", ""))[:12] or "commit not recorded"
    except (OSError, ValueError):
        from datetime import date  # noqa: PLC0415

        return f"retrieved {date.today().isoformat()}"


# --------------------------------------------------------------------- the run
def run(assembly: Path, db_dir: Path, *, threads: int = 1,
        workdir: Path | None = None) -> list[TypingResult]:
    assembly = Path(assembly)
    root = Path(db_dir) / KEY
    if workdir is not None:
        work = Path(workdir)
        work.mkdir(parents=True, exist_ok=True)
        return _run_in(assembly, root, work, threads)
    with tempdir(prefix="hydra.sagalactiae.") as temp:
        return _run_in(assembly, root, Path(temp), threads)


def _run_in(assembly: Path, root: Path, work: Path,
            threads: int) -> list[TypingResult]:
    results = [_call_cps(assembly, root, work, threads)]
    panel_fasta = root / "panel.faa"
    if not panel_fasta.exists():
        note = (f"the surface-protein panel is not installed under {root}; "
                f"run the module's fetch() first")
        for scheme in ("sag_surface", "sag_pilus", "sag_markers"):
            results.append(TypingResult(scheme=scheme, call="not determined",
                                        note=note))
        return results
    check_db_exists(panel_fasta, "prot")
    hits, _every = screen_panel(assembly, panel_fasta, work, PANEL,
                                threads=threads, tag=KEY)
    manifest = panel_manifest(root, "panel")
    results.append(_call_surface(hits, manifest))
    results.append(_call_pilus(hits))
    results.append(_call_markers(hits, manifest, results[0]))
    results.append(_call_pbp2x(assembly, root, work, threads))
    return results


# --------------------------------------------------------- PBP2X substitutions
def _call_pbp2x(assembly: Path, root: Path, work: Path,
                threads: int) -> TypingResult:
    """Name the PBP2X substitutions associated with reduced penicillin MICs.

    This is a flag for confirmatory phenotypic testing, never a resistance
    call. Group B streptococcus remains penicillin-susceptible in close to all
    isolates - the US surveillance figure is about one per cent with elevated
    MICs or pbp2x variants - and penicillin stays the treatment of choice unless
    an MIC says otherwise.
    """
    fasta = root / "pbp2x.faa"
    if not fasta.exists():
        return TypingResult(
            scheme="sag_pbp2x", call="not determined", loci_total=1,
            note=f"the PBP2X reference is not installed under {root}; run the "
                 f"module's fetch() first")
    check_db_exists(fasta, "prot")
    work.mkdir(parents=True, exist_ok=True)
    # Single-threaded on purpose: see strep_common.BLASTX_THREADS. pbp2x is a
    # one-protein database, so a negative sample is exactly the case that
    # crashes blastx when it is given several threads.
    hsps = blast("blastx", assembly, fasta, work / "pbp2x.blastx.tsv",
                 threads=BLASTX_THREADS, evalue=1e-40, with_alignment=True)
    if not hsps:
        return TypingResult(
            scheme="sag_pbp2x", call="not determined", loci_total=1,
            alleles={"reference": PBP2X_REFERENCE},
            note="pbp2x did not align to the reference at all; the gene is "
                 "essential, so this means it is missing from the assembly")
    best = max(hsps, key=lambda h: h.bitscore)
    observed: dict[int, str] = {}
    for hsp in sorted(hsps, key=lambda h: h.bitscore, reverse=True):
        pos = min(hsp.sstart, hsp.send) - 1
        for query_aa, ref_aa in zip(hsp.qseq, hsp.sseq):
            if ref_aa == "-":
                continue
            pos += 1
            observed.setdefault(pos, query_aa)

    covered = sorted(p for p in PBP2X_POSITIONS if p in observed)
    uncovered = sorted(p for p in PBP2X_POSITIONS if p not in observed)
    found: list[str] = []
    alleles: dict[str, str] = {
        "reference": f"{PBP2X_REFERENCE} ({PBP2X_REFERENCE_LEN} aa wild type)",
        "evidence": (f"{best.qseqid}:{best.q_lo}-{best.q_hi}, "
                     f"{best.pident:.1f}% aa identity over reference residues "
                     f"{min(best.sstart, best.send)}-"
                     f"{max(best.sstart, best.send)} of {PBP2X_REFERENCE_LEN}"),
    }
    accompanying: list[str] = []
    for pos in covered:
        wt, why = PBP2X_POSITIONS[pos]
        seen = observed[pos]
        if seen in (wt, "-", "*", "X"):
            continue
        name = f"{wt}{pos}{seen}"
        alleles[name] = why
        if pos == 377:
            # The reference carries V377 and most isolates carry I377, so this
            # difference is the reference's own allele showing, not a finding.
            alleles["V377I"] = (
                "the reference accession carries V at 377 where most isolates "
                "carry I, so this difference is wild type in the literature's "
                "numbering and is not a reduced-susceptibility substitution")
            continue
        if pos in PBP2X_DETERMINANTS:
            found.append(name)
        else:
            accompanying.append(name)
    alleles["accompanying_substitutions"] = (
        "+".join(accompanying) + " - reported alongside V405A in the PRGBS "
        "literature but not on their own associated with a raised MIC"
        if accompanying else "none")
    for start, motif, why in PBP2X_MOTIFS:
        seen = "".join(observed.get(start + i, "?") for i in range(len(motif)))
        if seen != motif:
            alleles[f"motif_{start}"] = (
                f"the {motif} motif at {start} reads {seen} in this isolate; "
                f"{why}. A change here is in the active site itself")

    notes: list[str] = []
    if uncovered:
        notes.append(
            f"{len(uncovered)} catalogued position(s) not covered by the "
            f"alignment ({', '.join(str(p) for p in uncovered)}), so "
            f"substitutions there are not determined")
    if found:
        escalating = [n for n in found if n.startswith("Q557")]
        call = "+".join(found)
        notes.insert(0,
                     "these are substitutions reported in group B streptococci "
                     "with reduced penicillin susceptibility. They are a flag "
                     "for confirmatory phenotypic MIC testing, not a resistance "
                     "call: penicillin remains the treatment of choice unless an "
                     "MIC says otherwise, and about 1% of isolates in "
                     "surveillance carry such variants")
        if escalating and any(n.startswith("V405") for n in found):
            notes.insert(1,
                         "V405A together with Q557E is the combination reported "
                         "at the highest penicillin MIC (about 1 mg/L); request "
                         "an MIC")
        elif escalating:
            notes.insert(1,
                         "Q557E is the only substitution NCBI curates for this "
                         "species; request an MIC")
    elif accompanying:
        call = "no reduced-susceptibility determinant"
        notes.insert(0,
                     f"neither V405A nor Q557E is present. The isolate does "
                     f"carry {'+'.join(accompanying)} at other catalogued "
                     f"positions; those are reported alongside V405A in the "
                     f"literature but are not on their own associated with a "
                     f"raised penicillin MIC")
    else:
        call = "no catalogued substitution"
        notes.insert(0,
                     "the transpeptidase domain matches the wild-type reference "
                     "at every catalogued position, which is the expected result "
                     "for a penicillin-susceptible isolate")
    return TypingResult(
        scheme="sag_pbp2x", call=call, lineage="-", alleles=alleles,
        loci_found=len(covered), loci_total=len(PBP2X_POSITIONS),
        note="; ".join(notes),
    )


# ----------------------------------------------------------- capsular serotype
class _Cps:
    """One candidate serotype, scored the way GBS-SBG scores it."""

    __slots__ = ("serotype", "ref", "ref_len", "bitscore", "spans", "contigs",
                 "n_hsps", "ident_num", "ident_den", "variants")

    def __init__(self, serotype: str, ref: str, ref_len: int) -> None:
        self.serotype = serotype
        self.ref = ref
        self.ref_len = ref_len
        self.bitscore = 0.0
        self.spans: list[tuple[int, int]] = []
        self.contigs: set[str] = set()
        self.n_hsps = 0
        self.ident_num = 0.0
        self.ident_den = 0.0
        self.variants: dict[str, float] = {}

    def add(self, hsp, variant: str = "") -> None:
        if variant:
            self.variants[variant] = self.variants.get(variant, 0.0) + hsp.bitscore
        self.bitscore += hsp.bitscore
        self.spans.append((min(hsp.sstart, hsp.send), max(hsp.sstart, hsp.send)))
        self.contigs.add(hsp.qseqid)
        self.n_hsps += 1
        self.ident_num += hsp.pident * hsp.length
        self.ident_den += hsp.length

    @property
    def coverage(self) -> float:
        """Fraction of the reference covered, merging overlapping alignments."""
        if not self.ref_len:
            return 0.0
        covered = 0
        end = 0
        for lo, hi in sorted(self.spans):
            if lo > end:
                covered += hi - lo + 1
                end = hi
            elif hi > end:
                covered += hi - end
                end = hi
        return covered / self.ref_len

    @property
    def identity(self) -> float:
        return (self.ident_num / self.ident_den) if self.ident_den else 0.0


def _call_cps(assembly: Path, root: Path, work: Path,
              threads: int) -> TypingResult:
    fasta = root / "cps.fna"
    if not fasta.exists():
        return TypingResult(
            scheme="cps", call="not determined", loci_total=1,
            note=f"the GBS-SBG capsular reference is not installed under {root}; "
                 f"run the module's fetch() first")
    check_db_exists(fasta, "nucl")
    work.mkdir(parents=True, exist_ok=True)
    version = ""
    stamp = root / "cps.version"
    if stamp.exists():
        version = stamp.read_text(encoding="utf-8").split("\t")[0].strip()

    hsps = blast("blastn", assembly, fasta, work / "cps.blastn.tsv",
                 threads=threads, evalue=1e-20, max_target_seqs=10000)
    # Accumulate per database record, not per serotype. The four serotype III
    # references are 99.3-99.6% identical to one another, so pooling them under
    # one "III" key sums four bit scores into one and serotype III then outranks
    # every other serotype in every isolate - which made a serotype Ia genome
    # rank as III, fail III's coverage test, and come out non-typeable.
    by_reference: dict[str, _Cps] = {}
    for hsp in hsps:
        if hsp.pident <= SBG_HSP_MIN_IDENTITY or hsp.length <= SBG_HSP_MIN_LENGTH:
            continue
        entry = by_reference.setdefault(
            hsp.sseqid,
            _Cps(_serotype_of(hsp.sseqid), hsp.sseqid, hsp.slen))
        entry.add(hsp, _variant_of(hsp.sseqid))

    ranked = sorted(by_reference.values(), key=lambda c: c.bitscore, reverse=True)
    database = f"GBS-SBG ({version or 'version not recorded'})"
    if not ranked:
        return TypingResult(
            scheme="cps", call="not determined", loci_total=1,
            alleles={"database": database},
            note="no part of the assembly aligned to the GBS-SBG capsular "
                 "reference; the isolate may not be Streptococcus agalactiae")

    best = ranked[0]
    # The runner-up that matters is the best reference of a *different* serotype.
    # Comparing against another serotype III variant would raise the "close
    # runner-up" flag on every serotype III isolate, which is the commonest
    # clinical serotype, and the flag would stop meaning anything.
    runner = next((c for c in ranked[1:] if c.serotype != best.serotype), None)
    sibling = next((c for c in ranked[1:] if c.serotype == best.serotype), None)
    passing = best.coverage >= SBG_MIN_COVERAGE and best.identity >= SBG_MIN_IDENTITY

    flags: list[str] = []
    if best.coverage < SBG_FLAG_COVERAGE:
        flags.append(f"coverage {best.coverage * 100:.1f}% of the reference, "
                     f"below {SBG_FLAG_COVERAGE * 100:.0f}%")
    if best.identity < SBG_FLAG_IDENTITY:
        flags.append(f"identity {best.identity:.1f}%, below "
                     f"{SBG_FLAG_IDENTITY:.0f}%")
    if best.n_hsps > SBG_FLAG_MAX_HSPS:
        flags.append(f"{best.n_hsps} alignments needed to cover the reference")
    if len(best.contigs) > SBG_FLAG_MAX_CONTIGS:
        flags.append(f"the locus is split across {len(best.contigs)} contigs")
    next_ratio = (runner.bitscore / best.bitscore) if runner and best.bitscore else 0.0
    if runner and next_ratio > SBG_FLAG_NEXT_RATIO:
        flags.append(f"the next best serotype {runner.serotype} reaches "
                     f"{next_ratio * 100:.0f}% of the winning bit score")

    subtype = (max(best.variants, key=lambda v: best.variants[v])
               if best.variants else "")
    if sibling is not None:
        alleles_sibling = (f"{_variant_of(sibling.ref)} "
                           f"({sibling.coverage * 100:.1f}% coverage, "
                           f"{sibling.identity:.1f}% identity)")
    else:
        alleles_sibling = "none"
    alleles = {
        "serotype": best.serotype if passing else "non-typeable",
        "cps_variant": (f"{subtype} (the serotype III reference variants differ "
                        f"only by SNPs in the conserved part of cps and are "
                        f"99.3-99.6% identical, so this is indicative, not a "
                        f"subtype call)"
                        if passing and subtype != best.serotype else "not applicable"),
        "reference": best.ref,
        "evidence": (f"{best.coverage * 100:.1f}% of the {best.ref_len} bp "
                     f"serotype-determining cps region covered at "
                     f"{best.identity:.1f}% identity, from {best.n_hsps} "
                     f"alignment(s) on {len(best.contigs)} contig(s)"),
        "next_best": (f"{runner.serotype} ({runner.coverage * 100:.1f}% "
                      f"coverage, {runner.identity:.1f}% identity, "
                      f"{next_ratio * 100:.0f}% of the winning bit score)"
                      if runner and next_ratio > SBG_REPORT_NEXT_RATIO
                      else "no close runner-up"),
        "same_serotype_runner_up": alleles_sibling,
        "uncertainty": "; ".join(flags) if flags else "none",
        "database": database,
        "vaccine_coverage": (
            f"covered by the hexavalent capsular conjugate "
            f"({', '.join(HEXAVALENT)})"
            if passing and _base_serotype(best.serotype) in HEXAVALENT
            else "not covered by the hexavalent capsular conjugate"
            if passing else "not determined"),
    }

    notes: list[str] = [
        "this is the capsular genotype encoded in cps, not a serological "
        "serotype: an isolate reported as non-typeable by serology usually "
        "carries an intact, typeable cps locus and is given a serotype here"]
    if not passing:
        notes.insert(0,
                     f"genotypically non-typeable: the best reference "
                     f"({best.serotype}) reaches only "
                     f"{best.coverage * 100:.1f}% coverage at "
                     f"{best.identity:.1f}% identity, below the GBS-SBG "
                     f"threshold of {SBG_MIN_COVERAGE * 100:.0f}% coverage and "
                     f"{SBG_MIN_IDENTITY:.0f}% identity. Either the cps locus "
                     f"is disrupted - deletions and insertion sequences in cps "
                     f"are the usual cause - or the assembly has broken inside "
                     f"it. Check the cpsA and neuA calls: if they are present "
                     f"and the serotype is not, the locus is there but its "
                     f"serotype-determining region is not resolved")
    if flags:
        notes.append("uncertainty flags raised: " + "; ".join(flags))
    if best.serotype.startswith("III"):
        notes.append("serotype III is represented by four reference variants in "
                     "the database (III-1 to III-4); all four are reported as "
                     "serotype III")

    return TypingResult(
        scheme="cps",
        call=alleles["serotype"],
        lineage=alleles["serotype"],
        alleles=alleles,
        loci_found=1 if passing else 0,
        loci_total=1,
        note="; ".join(notes),
    )


def _serotype_of(sseqid: str) -> str:
    """``1__GBS-SBG__GBS-SBG:III-2__5`` -> ``III``.

    The SRST2 header convention GBS-SBG uses puts the serotype in the third
    field. The four serotype III references are variants of one serotype and are
    collapsed, which is what the reference tool reports.
    """
    parts = sseqid.split("__")
    name = parts[2] if len(parts) > 2 else sseqid
    if ":" in name:
        name = name.split(":", 1)[1]
    return _base_serotype(name)


def _variant_of(sseqid: str) -> str:
    """The database record's own serotype label, keeping the III-n variant."""
    parts = sseqid.split("__")
    name = parts[2] if len(parts) > 2 else sseqid
    return name.split(":", 1)[1] if ":" in name else name


def _base_serotype(name: str) -> str:
    return name.split("-", 1)[0] if name.startswith("III-") else name


# --------------------------------------------------------- surface and markers
def _call_surface(hits: dict[str, GeneHit], manifest: dict) -> TypingResult:
    """The alpha-like protein family, reported as a family."""
    unresolved = set(manifest.get("unresolved") or ())
    alleles: dict[str, str] = {}
    named: list[str] = []
    for spec in ALP_FAMILY[:2]:          # bca and bac are separately annotated
        hit = hits.get(spec.name)
        if spec.name in unresolved:
            alleles[spec.name] = "not determined (no reference sequence)"
        elif hit is None or not hit.present:
            alleles[spec.name] = "not detected"
        else:
            alleles[spec.name] = f"{hit.state} - {hit.evidence()}"
            named.append(spec.name)

    family = [(spec, hits.get(spec.name)) for spec in ALP_FAMILY[2:]]
    present = [(spec, hit) for spec, hit in family
               if hit is not None and hit.present]
    if present:
        spec, hit = max(present, key=lambda pair: (pair[1].coverage_pct,
                                                   pair[1].identity_pct))
        alleles["alp_family"] = (
            f"an alpha-like protein family gene is present; closest reference "
            f"{hit.accession} ({spec.product}) at {hit.identity_pct:.1f}% aa "
            f"identity over {hit.coverage_pct:.0f}% of {hit.ref_len} aa, on "
            f"{hit.contig}:{hit.start}-{hit.end}")
        alleles["alp_family_member"] = (
            "not determined: rib, alp1/epsilon, alp2, alp3 and alp4 share their "
            "N-terminal domain and repeat units, and no versioned allele set is "
            "published, so the family member is not named")
        named.append("alp-family")
    else:
        alleles["alp_family"] = "not detected"
        alleles["alp_family_member"] = (
            "not determined: no alpha-like protein family gene was detected. "
            "Even when one is, the member is not named, because rib, "
            "alp1/epsilon, alp2, alp3 and alp4 share their N-terminal domain "
            "and repeat units and no versioned allele set is published")

    call = "+".join(named) if named else "none detected"
    notes = [
        "bca and bac carry large tandem repeats whose number varies between "
        "strains, so reference coverage well below 100% is normal and is not a "
        "lesion",
        "the alpha-like protein family member is not named; see the "
        "alp_family_member field for why",
    ]
    return TypingResult(
        scheme="sag_surface", call=call, lineage="-", alleles=alleles,
        loci_found=len(named), loci_total=3, note="; ".join(notes),
    )


def _call_pilus(hits: dict[str, GeneHit]) -> TypingResult:
    """The pilus island profile: PI-2a by name, the other locus unnamed."""
    alleles: dict[str, str] = {}
    pi2a_genes = ("pilA", "pilB", "pilC", "srtC4")
    other_genes = ("pilus_backbone_2", "pilus_backbone_3", "srtC1", "srtC2")
    for spec in PILUS:
        hit = hits.get(spec.name)
        alleles[spec.name] = (f"{hit.state} - {hit.evidence()}"
                              if hit is not None and hit.present
                              else "not detected")
    pi2a = [g for g in pi2a_genes
            if (h := hits.get(g)) is not None and h.present]
    other = [g for g in other_genes
             if (h := hits.get(g)) is not None and h.present]

    islands: list[str] = []
    if len(pi2a) >= 2:
        islands.append("PI-2a")
    if len(other) >= 2:
        islands.append("a second pilus island (not named)")
    alleles["PI-2a"] = (f"present ({'+'.join(pi2a)})" if len(pi2a) >= 2
                        else f"not detected (genes found: "
                             f"{'+'.join(pi2a) or 'none'})")
    alleles["second_island"] = (
        f"present ({'+'.join(other)}); PI-1 and PI-2b share the SpaH/EbpB "
        f"backbone and the class C sortases in the RefSeq annotation used here, "
        f"so this locus is detected but not assigned to PI-1 or PI-2b"
        if len(other) >= 2
        else f"not detected (genes found: {'+'.join(other) or 'none'})")
    notes = [
        "PI-1 and PI-2b are not distinguished: RefSeq annotates the backbone of "
        "both as a SpaH/EbpB family major pilin, and no versioned reference set "
        "that separates them is downloadable, so naming one from a best BLAST "
        "hit would be a guess",
        "PI-2a and PI-2b occupy the same chromosomal locus and are mutually "
        "exclusive, so an isolate carrying PI-2a does not carry PI-2b",
    ]
    return TypingResult(
        scheme="sag_pilus",
        call=" + ".join(islands) if islands else "no pilus island detected",
        lineage="-", alleles=alleles,
        loci_found=len(pi2a) + len(other), loci_total=len(PILUS),
        note="; ".join(notes),
    )


def _call_markers(hits: dict[str, GeneHit], manifest: dict,
                  cps: TypingResult) -> TypingResult:
    unresolved = set(manifest.get("unresolved") or ())
    alleles: dict[str, str] = {}
    found: list[str] = []
    notes: list[str] = []
    for spec in MARKERS:
        hit = hits.get(spec.name)
        if spec.name in unresolved:
            alleles[spec.name] = (
                f"not determined: {spec.meaning}"
                if spec.name == "hvgA" else
                "not determined (no reference protein in the current RefSeq "
                "annotation of the panel strains)")
            continue
        if hit is None or not hit.present:
            alleles[spec.name] = "not detected"
            continue
        alleles[spec.name] = f"{hit.state} - {hit.evidence()}"
        found.append(spec.name)

    srr2 = hits.get("srr2")
    cc17 = srr2 is not None and srr2.present
    alleles["CC17_marker"] = (
        "srr2 present: this is the ST-17 / CC17-restricted adhesin, so the "
        "isolate is probably the hypervirulent neonatal-meningitis clone. "
        "Confirm with the MLST sequence type; hvgA, the usual CC17 marker, is "
        "not called by this module"
        if cc17 else
        "srr2 not detected, which argues against CC17; srr1 is the adhesin of "
        "the other lineages")
    pbp = hits.get("pbp2x")
    if pbp is not None and pbp.present:
        if pbp.identity_pct >= 100.0:
            alleles["pbp2x_divergence"] = (
                "identical to the NEM316 wild-type PBP2X, so no substitution "
                "associated with reduced penicillin susceptibility is present")
        else:
            diffs = round((100.0 - pbp.identity_pct) * pbp.ref_len / 100.0)
            alleles["pbp2x_divergence"] = (
                f"about {diffs} residue difference(s) from the NEM316 "
                f"wild-type PBP2X ({pbp.identity_pct:.1f}% aa identity). "
                f"Named substitutions (for example in the transpeptidase "
                f"domain) are not called here, so this is a flag to inspect "
                f"the allele, not a susceptibility prediction")
            notes.append(
                "pbp2x differs from the wild-type reference; group B "
                "streptococcus remains penicillin-susceptible in almost all "
                "settings, and reduced susceptibility must be confirmed "
                "phenotypically")
    alleles["vaccine_relevance"] = (
        f"capsular: {cps.alleles.get('vaccine_coverage', 'not determined')}; "
        f"protein candidates detected: "
        + (", ".join(g for g in ("bca", "sip", "pilA") if g in found) or "none"))
    notes.append(
        "hvgA is not called: RefSeq does not name the gene in the reference "
        "assemblies this panel is built from and no versioned reference set is "
        "downloadable, so srr2 is reported as the CC17 marker instead")
    notes.append(
        "the clonal complex is not assigned here: it follows from the MLST "
        "sequence type, which Hydra reports separately")
    return TypingResult(
        scheme="sag_markers",
        call="+".join(found) if found else "none detected",
        lineage="CC17 (srr2-positive)" if cc17 else "-",
        alleles=alleles, loci_found=len(found),
        loci_total=len(MARKERS) - len(unresolved & {s.name for s in MARKERS}),
        note="; ".join(notes),
    )
