"""*Salmonella enterica*: serovar prediction and the typing a reference laboratory adds to it.

The serovar is the unit every *Salmonella* report is written in, but it is a
phenotype: an agglutination result against the somatic (O) and the two phases of
flagellar (H1, H2) antigens, read off the White-Kauffmann-Le Minor (WKL) table of
antigenic formulae. From an assembly we can only call the *genes* that encode
those antigens - the ``rfb`` cluster's ``wzx``/``wzy`` for O, ``fliC`` for H1 and
``fljB`` for H2 - and then look the resulting formula up in the same table. That
indirection is the source of every failure mode in the literature, so this module
reports the formula it found, the evidence behind each antigen, every serovar the
formula fits, and what could not be decided, rather than a bare name.

What it calls, in the order it is reported:

``salmonella_serovar``
    The antigenic formula and the serovar name(s) the WKL scheme gives it.
``salmonella_subspecies``
    Subspecies I-VI, inferred from which rows of the WKL table the formula fits.
``salmonella_h2_locus``
    The state of the ``fljBA``/``hin`` region, which is what separates a genuine
    monophasic variant (4,[5],12:i:- and friends) from an assembly that merely
    lost ``fljB`` to a contig break.
``salmonella_spi``
    The *Salmonella* pathogenicity island repertoire.
``salmonella_virulence_plasmid``
    The ``spv`` locus and the rest of the pSLT virulence plasmid backbone.
``salmonella_typhoidal``
    Vi capsule (``viaB``: ``tviA``-``tviE``, ``vexA``-``vexE``), typhoid toxin
    (``cdtB``/``pltA``/``pltB``), SPI-7 and the H:z66 flagellar variant.
``salmonella_mlst``
    The 7-locus Achtman sequence type and its eBurst group, which is the lineage
    call to read *beside* the serovar: most serovars are polyphyletic.

Every threshold and rule is set out in ``docs/typing/salmonella.md`` with the
paper it comes from; :data:`REFERENCES` carries the same list in code.
"""

from __future__ import annotations

import csv
import json
import re
import shutil
import tarfile
import time
import urllib.error
import urllib.request
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from ...engines.blast import MergedHit, blast, deduplicate, merge_hsps
from ...records import TypingResult
from ...seqio import read_fasta, write_fasta
from ...utils import LOG, HydraError, smart_open, tempdir
from ...utils import run as _run_tool
from .base import OrganismDb, Reference

KEY = "salmonella"
TITLE = "Salmonella enterica"

#: Species names (as Hydra's SpeciesCall reports them) that trigger the module.
#: *S. bongori* is included because the WKL table covers its serovars too, and
#: because a mis-sketched *bongori* genome should still be told it is not
#: *enterica* rather than silently skipped.
SPECIES: tuple[str, ...] = (
    "Salmonella enterica",
    "Salmonella bongori",
    "Salmonella sp.",
)

#: PubMLST scheme names that also trigger it. ``senterica_achtman_2`` is the
#: spelling Hydra's PubMLST importer produces for the 7-locus Achtman scheme;
#: ``salmonella`` is scheme 1 of the same database.
MLST_SCHEMES: tuple[str, ...] = ("senterica_achtman_2", "salmonella")

# ---------------------------------------------------------------------------
# Databases
# ---------------------------------------------------------------------------
_SISTR_DB_URL = ("https://zenodo.org/records/15078586/files/"
                 "SISTR_V_1.1.3_db.tar.gz?download=1")
_SISTR_MEMBERS = (
    "data/Salmonella-serotype_serogroup_antigen_table-WHO_2007.csv",
    "data/antigens/wzx.fasta",
    "data/antigens/wzy.fasta",
    "data/antigens/fliC.fasta",
    "data/antigens/fljB.fasta",
)
_SPIFINDER_URL = "https://bitbucket.org/genomicepidemiology/spifinder_db/get/HEAD.tar.gz"
_PUBMLST_BASE = "https://rest.pubmlst.org/db/pubmlst_salmonella_seqdef"
_MLST_LOCI = ("aroC", "dnaN", "hemD", "hisD", "purE", "sucA", "thrA")
_NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

DATABASES: tuple[OrganismDb, ...] = (
    OrganismDb(
        name="antigens",
        source=_SISTR_DB_URL,
        version="SISTR cmd database 1.1.3, Zenodo doi:10.5281/zenodo.15078586, "
                "published 2025-03-24",
        licence="Apache-2.0 (Zenodo record licence: apache2.0; sistr_cmd itself "
                "states 'license=Apache 2.0')",
        redistribution="permitted",
        reference="Yoshida2016",
        notes="wzx.fasta, wzy.fasta, fliC.fasta and fljB.fasta: the O serogroup "
              "and H1/H2 antigen allele sets, each sequence labelled with the "
              "antigen it encodes. 2.3 MB of a 50 MB archive; the rest "
              "(cgMLST330 and a Mash sketch) is not downloaded. The archive is "
              "fetched whole because Zenodo serves no per-member access.",
    ),
    OrganismDb(
        name="serovars",
        source=_SISTR_DB_URL,
        version="Salmonella-serotype_serogroup_antigen_table-WHO_2007.csv from "
                "SISTR cmd database 1.1.3 (file dated 2024-09-04); content "
                "derives from the WKL scheme 9th edition (2007)",
        licence="Apache-2.0 for the CSV as distributed by SISTR. The underlying "
                "table is Grimont & Weill, 'Antigenic Formulae of the Salmonella "
                "Serovars', 9th ed., WHO Collaborating Centre, Institut Pasteur, "
                "which states no licence; it is published for free download at "
                "pasteur.fr.",
        redistribution="permitted",
        reference="Grimont2007",
        notes="2,660 serovars with O antigen, serogroup, H1, H2, whether H2 may "
              "be absent, and subspecies. This is the table the serovar name is "
              "read out of; nothing else in the module invents a name.",
    ),
    OrganismDb(
        name="spi",
        source="https://bitbucket.org/genomicepidemiology/spifinder_db",
        version="commit db102668b704, files dated 2021-01-20, retrieved 2026-11-17",
        licence="Apache-2.0 ('Licensed under the Apache License, Version 2.0', "
                "README.md)",
        redistribution="permitted",
        reference="Roer2016",
        notes="84 pathogenicity-island reference sequences (SPI-1..SPI-14, CS54, "
              "C63PI, SGI1, SGI2, HPI, SESS-LEE), 1.3 MB. Derived in turn from "
              "PAIDB v2.0 (Yoon2015).",
    ),
    OrganismDb(
        name="markers",
        source="NCBI Nucleotide, records NC_003197.2 (S. Typhimurium LT2 "
               "chromosome), NC_003277.2 (pSLT), NC_003198.1 (S. Typhi CT18 "
               "chromosome), AM419040.1 (pBSSB1), AF370707.1 (Sdf-I)",
        version="retrieved from NCBI on first fetch; the record versions above "
                "are pinned in the request",
        licence="none stated. NCBI's disclaimer reads: 'Some NCBI data ... may "
                "be protected by copyright ... NCBI itself places no additional "
                "restrictions on the use or distribution of the data contained "
                "therein.' GenBank sequence records are treated as freely "
                "redistributable.",
        redistribution="permitted",
        reference="McClelland2001",
        notes="34 single genes selected by gene name / locus tag: the fljBA-hin "
              "region, oafA, invA, the spv locus and pSLT backbone, the viaB "
              "locus, typhoid toxin, SPI-7 pilL, sefA, sopE, the pBSSB1 z66 "
              "fljB allele and the Enteritidis-specific Sdf-I fragment. About "
              "40 kB once extracted.",
    ),
    OrganismDb(
        name="mlst",
        source=_PUBMLST_BASE + " (scheme 2, the 7-locus Achtman scheme)",
        version="PubMLST allele and profile download, dated at fetch time; "
                "unauthenticated downloads are served data submitted on or "
                "before 2024-12-31",
        licence="PubMLST Terms & Conditions section 3: data submitted on or "
                "before 31 December 2024 'may continue to be downloaded, used, "
                "and redistributed without restriction, subject only to proper "
                "scientific citation'. Section 4 restricts data submitted on or "
                "after 1 January 2025 to 'non-commercial academic research or "
                "public health surveillance'; the unauthenticated endpoint this "
                "module uses does not serve that data.",
        redistribution="permitted",
        reference="Jolley2018",
        notes="7 locus allele files (~6 MB) and the profile table, whose "
              "clonal_complex column carries the eBurst group (eBG). Kept "
              "separate from Hydra's own PubMLST install so the module can "
              "report ST and eBG from its own evidence.",
    ),
)

# ---------------------------------------------------------------------------
# Thresholds. Each one is justified in docs/typing/salmonella.md.
# ---------------------------------------------------------------------------
#: wzx/wzy: SISTR discards a serogroup hit below 88% identity, below 300 aligned
#: bases, or between 300 and 500 bases unless it is >=99% identical.
O_MIN_IDENTITY = 88.0
O_MIN_BASES = 300
O_SHORT_BASES = 500
O_SHORT_IDENTITY = 99.0
#: fliC/fljB. A confident H call needs a near-complete allele; a provisional one
#: needs most of it. The allele sets differ by a handful of bases between some
#: antigen types (the g-complex), so identity is held high.
H_CONFIDENT_BASES = 1000
H_CONFIDENT_IDENTITY = 99.0
H_PROVISIONAL_BASES = 700
H_PROVISIONAL_IDENTITY = 97.0
#: A hit that stops this close to the end of its contig may be truncated by the
#: assembly rather than by biology.
CONTIG_EDGE_BP = 300
#: Single marker genes: presence, and the band below it that is reported as
#: partial rather than absent.
MARKER_IDENTITY = 90.0
MARKER_COVERAGE = 80.0
MARKER_PARTIAL_COVERAGE = 40.0
#: Pathogenicity islands are mosaic, so coverage is reported rather than
#: thresholded into presence alone.
SPI_IDENTITY = 85.0
SPI_PRESENT_COVERAGE = 90.0
SPI_PARTIAL_COVERAGE = 50.0
SPI_PARTIAL_CEILING = SPI_PRESENT_COVERAGE
#: Below this a hit is the islands' shared backbone genes, not the island.
SPI_TRACE_COVERAGE = 15.0
#: MLST: an allele must match exactly to be named.
MLST_IDENTITY = 100.0
MLST_COVERAGE = 100.0

#: O groups whose wzx/wzy alleles are too close to separate reliably. Taken
#: from SISTR's serogroup similarity groups, with the spelling the WHO 2007
#: table actually uses ("O:62", not "62"). These are reported as an explicit
#: ambiguity beside the call rather than merged into it, because widening every
#: D1 call to include A and D2 would cost more than it buys.
O_GROUP_SIMILARITY: tuple[tuple[str, ...], ...] = (
    ("E1", "E4"),
    ("A", "D1", "D2"),
    ("C1", "F"),
    ("S", "O:62"),
)

_SUBSPECIES_ROMAN = {
    "enterica": "I",
    "salamae": "II",
    "arizonae": "IIIa",
    "diarizonae": "IIIb",
    "houtenae": "IV",
    "bongori": "V",
    "indica": "VI",
}

# ---------------------------------------------------------------------------
# Marker genes, selected from pinned reference records by gene name or locus tag
# ---------------------------------------------------------------------------
#: (accession, ((gene, locus_tag, group), ...)). A CDS matches on either the
#: ``[gene=]`` qualifier or the ``[locus_tag=]`` qualifier, because RefSeq
#: renumbered the Typhi CT18 locus tags while the gene names stayed put.
_MARKER_CDS: tuple[tuple[str, tuple[tuple[str, str, str], ...]], ...] = (
    ("NC_003197.2", (
        ("fliC", "STM1959", "h1_locus"),
        ("fljA", "STM2770", "h2_locus"),
        ("fljB", "STM2771", "h2_locus"),
        ("hin", "STM2772", "h2_locus"),
        ("iroB", "STM2773", "h2_locus_flank"),
        ("oafA", "STM2232", "o_antigen"),
        ("invA", "STM2896", "genus"),
    )),
    ("NC_003277.2", (
        ("spvR", "PSLT041", "spv"),
        ("spvA", "PSLT040", "spv"),
        ("spvB", "PSLT039", "spv"),
        ("spvC", "PSLT038", "spv"),
        ("spvD", "PSLT037", "spv"),
        ("pefA", "PSLT018", "pslt"),
        ("rcK", "PSLT009", "pslt"),
        ("traT", "PSLT103", "pslt"),
        ("samA", "PSLT055", "pslt"),
    )),
    ("NC_003198.1", (
        ("tviA", "STY_RS22270", "viaB"),
        ("tviB", "STY_RS22265", "viaB"),
        ("tviC", "STY_RS22260", "viaB"),
        ("tviD", "STY_RS22255", "viaB"),
        ("tviE", "STY_RS22245", "viaB"),
        ("vexA", "STY_RS22240", "viaB"),
        ("vexB", "STY_RS22235", "viaB"),
        ("vexC", "STY_RS22230", "viaB"),
        ("vexD", "STY_RS22225", "viaB"),
        ("vexE", "STY_RS22220", "viaB"),
        ("cdtB", "STY_RS08895", "typhoid_toxin"),
        ("pltA", "STY_RS08910", "typhoid_toxin"),
        ("pltB", "STY_RS08915", "typhoid_toxin"),
        ("sopE", "STY_RS22015", "effector"),
        ("pilL", "STY_RS21690", "spi7"),
        ("sefA", "STY_RS23155", "fimbriae"),
    )),
    ("AM419040.1", (
        ("fljB", "BSSB1_031", "z66"),
    )),
)
#: Whole records taken as one marker sequence.
_MARKER_WHOLE: tuple[tuple[str, str, str], ...] = (
    ("AF370707.1", "sdf-I", "enteritidis"),
)
#: The marker name used in reports, when it must differ from the gene name.
_MARKER_RENAME = {("AM419040.1", "fljB"): "fljB_z66"}

#: What each marker group means, printed in the note so a reader does not have
#: to know the gene names.
_GROUP_MEANING = {
    "spv": "Salmonella plasmid virulence operon (spvRABCD): systemic virulence "
           "in non-typhoidal serovars",
    "pslt": "pSLT virulence-plasmid backbone",
    "viaB": "Vi capsular polysaccharide biosynthesis and export (viaB region of "
            "SPI-7)",
    "typhoid_toxin": "typhoid toxin (A2B5: CdtB plus PltA/PltB)",
}

REFERENCES: tuple[Reference, ...] = (
    Reference(
        key="Grimont2007",
        authors="Grimont PAD, Weill FX.",
        year=2007,
        title="Antigenic Formulae of the Salmonella Serovars, 9th edition",
        journal="WHO Collaborating Centre for Reference and Research on "
                "Salmonella, Institut Pasteur, Paris, 166 pp.",
        doi="",
        pmid="",
        supports="The White-Kauffmann-Le Minor table itself: every serovar name "
                "this module reports is read out of it, and a formula with no "
                "name in it is reported as a formula.",
    ),
    Reference(
        key="Guibourdenche2010",
        authors="Guibourdenche M, Roggentin P, Mikoleit M, et al.",
        year=2010,
        title="Supplement 2003-2007 (No. 47) to the White-Kauffmann-Le Minor "
              "scheme",
        journal="Research in Microbiology 161:26-29",
        doi="10.1016/j.resmic.2009.10.002",
        pmid="19840847",
        supports="Serovars added to the scheme after the 9th edition; the "
                 "reason the shipped table can lag the current nomenclature.",
    ),
    Reference(
        key="IssenhuthJeanjean2014",
        authors="Issenhuth-Jeanjean S, Roggentin P, Mikoleit M, et al.",
        year=2014,
        title="Supplement 2008-2010 (no. 48) to the White-Kauffmann-Le Minor "
              "scheme",
        journal="Research in Microbiology 165:526-530",
        doi="10.1016/j.resmic.2014.07.004",
        pmid="25049166",
        supports="The most recent supplement to the scheme; cited with the "
                 "caveat that the shipped table is the 2007 edition.",
    ),
    Reference(
        key="Liu2014",
        authors="Liu B, Knirel YA, Feng L, et al.",
        year=2014,
        title="Structural diversity in Salmonella O antigens and its genetic "
              "basis",
        journal="FEMS Microbiology Reviews 38:56-89",
        doi="10.1111/1574-6976.12034",
        pmid="23848592",
        supports="That wzx and wzy inside the rfb cluster determine the O "
                 "serogroup, which is why those two genes are what the O "
                 "antigen is called from.",
    ),
    Reference(
        key="Zhang2015",
        authors="Zhang S, Yin Y, Jones MB, et al.",
        year=2015,
        title="Salmonella serotype determination utilizing high-throughput "
              "genome sequencing data",
        journal="Journal of Clinical Microbiology 53:1685-1692",
        doi="10.1128/JCM.00323-15",
        pmid="25762776",
        supports="SeqSero: the determinant set (rfb cluster, fliC and fljB "
                 "alleles) that in-silico serovar prediction rests on.",
    ),
    Reference(
        key="Zhang2019",
        authors="Zhang S, den Bakker HC, Li S, et al.",
        year=2019,
        title="SeqSero2: rapid and improved Salmonella serotype determination "
              "using whole-genome sequencing data",
        journal="Applied and Environmental Microbiology 85:e01746-19",
        doi="10.1128/AEM.01746-19",
        pmid="31540993",
        supports="wzx/wzy for O and fliC/fljB for H; the use of sdf to separate "
                 "Enteritidis from Gallinarum and of oafA for O:5-negative "
                 "Typhimurium; the 3% of genomes that hold more than one "
                 "serotype.",
    ),
    Reference(
        key="Yoshida2016",
        authors="Yoshida CE, Kruczkiewicz P, Laing CR, et al.",
        year=2016,
        title="The Salmonella In Silico Typing Resource (SISTR): an open "
              "web-accessible tool for rapidly typing and subtyping draft "
              "Salmonella genome assemblies",
        journal="PLoS ONE 11:e0147101",
        doi="10.1371/journal.pone.0147101",
        pmid="26800248",
        supports="The antigen allele database and the WHO 2007 serovar table "
                 "this module installs, and the BLAST thresholds for wzx/wzy "
                 "and fliC/fljB that it follows.",
    ),
    Reference(
        key="Agron2001",
        authors="Agron PG, Walker RL, Kinde H, et al.",
        year=2001,
        title="Identification by subtractive hybridization of sequences "
              "specific for Salmonella enterica serovar Enteritidis",
        journal="Applied and Environmental Microbiology 67:4984-4991",
        doi="10.1128/AEM.67.11.4984-4991.2001",
        pmid="11679316",
        supports="The Sdf-I fragment (GenBank AF370707) used to tell commonly "
                 "circulating Enteritidis from Gallinarum and Pullorum, which "
                 "share its antigen genes.",
    ),
    Reference(
        key="Uelze2020",
        authors="Uelze L, Borowiak M, Deneke C, et al.",
        year=2020,
        title="Performance and accuracy of four open-source tools for in silico "
              "serotyping of Salmonella spp. based on whole-genome short-read "
              "sequencing data",
        journal="Applied and Environmental Microbiology 86:e02265-19",
        doi="10.1128/AEM.02265-19",
        pmid="31862714",
        supports="The accuracy ceiling quoted in the report (SISTR 94%, "
                 "SeqSero2 87% over 1,624 reference-laboratory isolates), that "
                 "monophasic variants are where tools disagree, and that the "
                 "low GC content of wzx/wzy makes O antigen recovery "
                 "library-dependent.",
    ),
    Reference(
        key="Banerji2020",
        authors="Banerji S, Simon S, Tille A, et al.",
        year=2020,
        title="Genome-based Salmonella serotyping as the new gold standard",
        journal="Scientific Reports 10:4333",
        doi="10.1038/s41598-020-61254-1",
        pmid="32152449",
        supports="That antigen-based and MLST-based calls should be read "
                 "together (98% and 95% alone, >99% combined), and that both "
                 "misidentify monophasic Typhimurium.",
    ),
    Reference(
        key="Diep2019",
        authors="Diep B, Barretto C, Portmann AC, et al.",
        year=2019,
        title="Salmonella serotyping; comparison of the traditional method to a "
              "microarray-based method and an in silico platform using whole "
              "genome sequencing data",
        journal="Frontiers in Microbiology 10:2554",
        doi="10.3389/fmicb.2019.02554",
        pmid="31781065",
        supports="That a third of in-silico calls fit more than one serovar, "
                 "which is why this module reports a candidate list rather than "
                 "one name.",
    ),
    Reference(
        key="Robertson2018",
        authors="Robertson J, Yoshida C, Kruczkiewicz P, et al.",
        year=2018,
        title="Comprehensive assessment of the quality of Salmonella whole "
              "genome sequence data available in public sequence databases "
              "using the Salmonella in silico Typing Resource (SISTR)",
        journal="Microbial Genomics 4:e000151",
        doi="10.1099/mgen.0.000151",
        pmid="29338812",
        supports="That public-database serovar labels are themselves often "
                 "wrong, which is why the validation set uses strains with a "
                 "published serovar rather than submitter metadata alone.",
    ),
    Reference(
        key="Soyer2009",
        authors="Soyer Y, Moreno Switt A, Davis MA, et al.",
        year=2009,
        title="Salmonella enterica serotype 4,5,12:i:-, an emerging Salmonella "
              "serotype that represents multiple distinct clones",
        journal="Journal of Clinical Microbiology 47:3546-3556",
        doi="10.1128/JCM.00546-09",
        pmid="19741087",
        supports="That 4,[5],12:i:- is not one lineage, so a monophasic antigen "
                 "call is not by itself an identification of the European "
                 "clone.",
    ),
    Reference(
        key="Garcia2016",
        authors="García P, Malorny B, Rodicio MR, et al.",
        year=2016,
        title="Horizontal acquisition of a multidrug-resistance module (R-type "
              "ASSuT) is responsible for the monophasic phenotype in a "
              "widespread clone of Salmonella serovar 4,[5],12:i:-",
        journal="Frontiers in Microbiology 7:680",
        doi="10.3389/fmicb.2016.00680",
        pmid="27242707",
        supports="The mechanism this module looks for: in the widespread clone "
                 "the fljAB region is physically replaced by a resistance "
                 "module, so fljA and hin go with fljB.",
    ),
    Reference(
        key="Prendergast2013",
        authors="Prendergast DM, Hand D, Nic Ghallchoir E, et al.",
        year=2013,
        title="A multiplex real-time PCR assay for the identification and "
              "differentiation of Salmonella enterica serovar Typhimurium and "
              "monophasic serovar 4,[5],12:i:-",
        journal="International Journal of Food Microbiology 166:48-53",
        doi="10.1016/j.ijfoodmicro.2013.05.031",
        pmid="23827807",
        supports="That the fljAB region, not fljB alone, is what a reference "
                 "method interrogates to confirm a monophasic variant; the "
                 "basis of the salmonella_h2_locus call.",
    ),
    Reference(
        key="Petrovska2016",
        authors="Petrovska L, Mather AE, AbuOun M, et al.",
        year=2016,
        title="Microevolution of monophasic Salmonella Typhimurium during "
              "epidemic, United Kingdom, 2005-2010",
        journal="Emerging Infectious Diseases 22:617-624",
        doi="10.3201/eid2204.150531",
        pmid="26982594",
        supports="The ST34 epidemic monophasic clone, reported when the "
                 "monophasic formula and ST34 coincide.",
    ),
    Reference(
        key="Achtman2012",
        authors="Achtman M, Wain J, Weill FX, et al.",
        year=2012,
        title="Multilocus sequence typing as a replacement for serotyping in "
              "Salmonella enterica",
        journal="PLoS Pathogens 8:e1002776",
        doi="10.1371/journal.ppat.1002776",
        pmid="22737074",
        supports="The eBurst group (eBG) reported beside the ST, and the finding "
                 "that most serovars are polyphyletic - the reason the ST is "
                 "reported as a separate call and not folded into the serovar.",
    ),
    Reference(
        key="Alikhan2018",
        authors="Alikhan NF, Zhou Z, Sergeant MJ, Achtman M.",
        year=2018,
        title="A genomic overview of the population structure of Salmonella",
        journal="PLoS Genetics 14:e1007261",
        doi="10.1371/journal.pgen.1007261",
        pmid="29621240",
        supports="Serovar-to-lineage concordance at genome scale; the authority "
                 "for checking a serovar/ST disagreement rather than resolving "
                 "it here.",
    ),
    Reference(
        key="Chattaway2021",
        authors="Chattaway MA, Langridge GC, Wain J.",
        year=2021,
        title="Salmonella nomenclature in the genomic era: a time for change",
        journal="Scientific Reports 11:7494",
        doi="10.1038/s41598-021-86243-w",
        pmid="33820940",
        supports="Why a serovar name and a genomic lineage are reported side by "
                 "side instead of one being collapsed into the other.",
    ),
    Reference(
        key="Jolley2018",
        authors="Jolley KA, Bray JE, Maiden MCJ.",
        year=2018,
        title="Open-access bacterial population genomics: BIGSdb software, the "
              "PubMLST.org website and their applications",
        journal="Wellcome Open Research 3:124",
        doi="10.12688/wellcomeopenres.14826.1",
        pmid="30345391",
        supports="The PubMLST allele and profile data the MLST call uses; the "
                 "citation PubMLST asks for.",
    ),
    Reference(
        key="Roer2016",
        authors="Roer L, Hendriksen RS, Leekitcharoenphon P, et al.",
        year=2016,
        title="Is the evolution of Salmonella enterica subsp. enterica linked to "
              "restriction-modification systems?",
        journal="mSystems 1:e00009-16",
        doi="10.1128/mSystems.00009-16",
        pmid="27822532",
        supports="SPIFinder and its pathogenicity-island reference database, "
                 "which the salmonella_spi call is made against.",
    ),
    Reference(
        key="Yoon2015",
        authors="Yoon SH, Park YK, Kim JF.",
        year=2015,
        title="PAIDB v2.0: exploration and analysis of pathogenicity and "
              "resistance islands",
        journal="Nucleic Acids Research 43:D624-D630",
        doi="10.1093/nar/gku985",
        pmid="25336619",
        supports="The island boundaries and names the SPIFinder database "
                 "inherits.",
    ),
    Reference(
        key="McClelland2001",
        authors="McClelland M, Sanderson KE, Spieth J, et al.",
        year=2001,
        title="Complete genome sequence of Salmonella enterica serovar "
              "Typhimurium LT2",
        journal="Nature 413:852-856",
        doi="10.1038/35101614",
        pmid="11677609",
        supports="The source of the fljBA-hin region, oafA, invA and the pSLT "
                 "spv marker sequences, and the published serovar of the "
                 "validation strain LT2.",
    ),
    Reference(
        key="Parkhill2001",
        authors="Parkhill J, Dougan G, James KD, et al.",
        year=2001,
        title="Complete genome sequence of a multiple drug resistant Salmonella "
              "enterica serovar Typhi CT18",
        journal="Nature 413:848-852",
        doi="10.1038/35101607",
        pmid="11677608",
        supports="The source of the viaB, typhoid toxin, SPI-7 and sefA marker "
                 "sequences, and the published serovar of validation strain "
                 "CT18.",
    ),
    Reference(
        key="Baker2007",
        authors="Baker S, Hardy J, Sanderson KE, et al.",
        year=2007,
        title="A novel linear plasmid mediates flagellar variation in "
              "Salmonella Typhi",
        journal="PLoS Pathogens 3:e59",
        doi="10.1371/journal.ppat.0030059",
        pmid="17500588",
        supports="The H:z66 flagellar antigen and its linear plasmid pBSSB1 "
                 "(GenBank AM419040), the source of the fljB z66 allele.",
    ),
    Reference(
        key="Thomson2008",
        authors="Thomson NR, Clayton DJ, Windhorst D, et al.",
        year=2008,
        title="Comparative genome analysis of Salmonella Enteritidis PT4 and "
              "Salmonella Gallinarum 287/91 provides insights into evolutionary "
              "and host adaptation pathways",
        journal="Genome Research 18:1624-1637",
        doi="10.1101/gr.077404.108",
        pmid="18583645",
        supports="The published serovars of validation strains P125109 "
                 "(Enteritidis) and 287/91 (Gallinarum), and that Gallinarum "
                 "carries the Enteritidis antigen genes while being "
                 "non-motile.",
    ),
    Reference(
        key="Fricke2011",
        authors="Fricke WF, Mammel MK, McDermott PF, et al.",
        year=2011,
        title="Comparative genomics of 28 Salmonella enterica isolates: evidence "
              "for CRISPR-mediated adaptive sublineage evolution",
        journal="Journal of Bacteriology 193:3556-3568",
        doi="10.1128/JB.00297-11",
        pmid="21602358",
        supports="The published serovars of validation strains SL476 "
                 "(Heidelberg) and SL254 (Newport).",
    ),
    Reference(
        key="Aviv2014",
        authors="Aviv G, Tsyba K, Steck N, et al.",
        year=2014,
        title="A unique megaplasmid contributes to stress tolerance and "
              "pathogenicity of an emergent Salmonella enterica serovar "
              "Infantis strain",
        journal="Environmental Microbiology 16:977-994",
        doi="10.1111/1462-2920.12351",
        pmid="24320043",
        supports="The published serovar of validation strain 119944 (Infantis) "
                 "and the pESI megaplasmid that lineage carries.",
    ),
    Reference(
        key="McClelland2004",
        authors="McClelland M, Sanderson KE, Clifton SW, et al.",
        year=2004,
        title="Comparison of genome degradation in Paratyphi A and Typhi, human-"
              "restricted serovars of Salmonella enterica that cause typhoid",
        journal="Nature Genetics 36:1268-1274",
        doi="10.1038/ng1470",
        pmid="15531882",
        supports="The published serovar of validation strain ATCC 9150 "
                 "(Paratyphi A).",
    ),
)


# ---------------------------------------------------------------------------
# fetch
# ---------------------------------------------------------------------------
_USER_AGENT = "hydra-amr (+https://github.com/iowa69/hydra)"
_TIMEOUT = 180


def _download(url: str, dest: Path, attempts: int = 4) -> Path:
    """Fetch *url* to *dest*, retrying transport failures."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
            with urllib.request.urlopen(request, timeout=_TIMEOUT) as response, \
                    open(dest, "wb") as handle:
                shutil.copyfileobj(response, handle)
            if dest.stat().st_size == 0:
                raise HydraError(f"{url} returned an empty file")
            return dest
        except urllib.error.HTTPError as exc:
            last = exc
            if 400 <= exc.code < 500 and exc.code not in (429,):
                raise HydraError(f"could not download {url}: {exc}") from exc
            time.sleep(5.0 * (attempt + 1) if exc.code == 429 else 2.0 * (attempt + 1))
        except (urllib.error.URLError, OSError) as exc:
            last = exc
            time.sleep(2.0 * (attempt + 1))
    raise HydraError(f"could not download {url}: {last}")


def _makeblastdb(fasta: Path) -> None:
    _run_tool(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl",
               "-title", fasta.stem])


def _built(*paths: Path) -> bool:
    """True when every output exists and the FASTAs among them are indexed.

    A half-finished fetch that wrote the FASTA but died before makeblastdb must
    not be treated as complete, or the next run reports a missing index instead
    of finishing the job.
    """
    for path in paths:
        if not path.exists():
            return False
        if path.suffix in (".fna", ".fasta") and not any(
                Path(f"{path}{suffix}").exists() for suffix in (".nin", ".00.nin", ".nal")):
            return False
    return True


def _clean_id(text: str) -> str:
    """A BLAST-safe sequence id component: no whitespace, no pipe."""
    return re.sub(r"\s+", "", str(text)).replace("|", "/")


def _stage_antigens(root: Path, work: Path) -> None:
    """Install the wzx/wzy/fliC/fljB allele sets and the WKL serovar table."""
    out_fasta = root / "antigens" / "antigens.fna"
    out_table = root / "antigens" / "serovars.csv"
    if _built(out_fasta, out_table):
        return
    archive = work / "sistr_db.tar.gz"
    _download(_SISTR_DB_URL, archive)
    extracted = work / "sistr"
    extracted.mkdir(parents=True, exist_ok=True)
    wanted = set(_SISTR_MEMBERS)
    with tarfile.open(archive) as tar:
        members = [m for m in tar.getmembers() if m.name in wanted]
        missing = wanted - {m.name for m in members}
        if missing:
            raise HydraError("SISTR database archive is missing "
                             + ", ".join(sorted(missing)))
        for member in members:
            target = (extracted / member.name).resolve()
            if not str(target).startswith(str(extracted.resolve())):
                raise HydraError(f"unsafe path in SISTR archive: {member.name}")
        try:
            tar.extractall(extracted, members=members, filter="data")
        except TypeError:  # Python < 3.12 has no filter argument
            tar.extractall(extracted, members=members)
    records: list[tuple[str, str]] = []
    provenance: list[str] = []
    for locus in ("wzx", "wzy", "fliC", "fljB"):
        path = extracted / "data" / "antigens" / f"{locus}.fasta"
        for index, (header, sequence) in enumerate(read_fasta(path), start=1):
            antigen = _clean_id(header.rsplit("|", 1)[-1])
            if not antigen:
                continue
            new_id = f"{locus}|{index}|{antigen}"
            records.append((new_id, sequence.upper()))
            provenance.append(f"{new_id}\t{header}")
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(out_fasta, records, wrap=0)
    (out_fasta.parent / "antigens.sources.tsv").write_text(
        "\n".join(provenance) + "\n", encoding="utf-8")
    shutil.copyfile(
        extracted / "data" / "Salmonella-serotype_serogroup_antigen_table-WHO_2007.csv",
        out_table)
    _makeblastdb(out_fasta)
    archive.unlink(missing_ok=True)


def _stage_spi(root: Path, work: Path) -> None:
    """Install the SPIFinder pathogenicity-island references."""
    out_fasta = root / "spi" / "spi.fna"
    if _built(out_fasta):
        return
    archive = work / "spifinder.tar.gz"
    _download(_SPIFINDER_URL, archive)
    extracted = work / "spifinder"
    extracted.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            target = (extracted / member.name).resolve()
            if not str(target).startswith(str(extracted.resolve())):
                raise HydraError(f"unsafe path in SPIFinder archive: {member.name}")
        try:
            tar.extractall(extracted, filter="data")
        except TypeError:
            tar.extractall(extracted)
    found = sorted(extracted.rglob("SPI.fsa"))
    if not found:
        raise HydraError("SPIFinder archive did not contain SPI.fsa")
    records = []
    for header, sequence in read_fasta(found[0]):
        # SPIFinder ids read "SPI-1;4;Salmonella-Typhimurium-LT2;NC_003197_P3":
        # island, variant number, source strain, source accession.
        parts = header.split(";")
        island = _clean_id(parts[0]) or "unknown"
        variant = _clean_id(parts[1]) if len(parts) > 1 else "1"
        accession = _clean_id(parts[-1]) if len(parts) > 3 else "unknown"
        records.append((f"{island}|{variant}|{accession}", sequence.upper()))
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(out_fasta, records, wrap=0)
    notes = sorted(extracted.rglob("notes.txt"))
    if notes:
        shutil.copyfile(notes[0], out_fasta.parent / "spi.notes.txt")
    _makeblastdb(out_fasta)
    archive.unlink(missing_ok=True)


def _parse_cds_fasta(path: Path):
    """Yield (qualifiers, sequence) for an NCBI ``fasta_cds_na`` download."""
    for header, sequence in read_fasta(path):
        quals = dict(re.findall(r"\[([A-Za-z_]+)=([^\]]*)\]", header))
        yield quals, sequence.upper()


def _stage_markers(root: Path, work: Path) -> None:
    """Build the single-gene marker set from pinned GenBank records."""
    out_fasta = root / "markers" / "markers.fna"
    if _built(out_fasta):
        return
    records: list[tuple[str, str]] = []
    provenance: list[str] = ["marker\tgroup\tsource_accession\tlocus_tag\tlength"]
    for accession, wanted in _MARKER_CDS:
        local = work / f"{accession}.cds.fna"
        _download(f"{_NCBI_EFETCH}?db=nuccore&id={accession}"
                  f"&rettype=fasta_cds_na&retmode=text", local)
        by_gene: dict[str, tuple[str, str]] = {}
        by_tag: dict[str, tuple[str, str]] = {}
        for quals, sequence in _parse_cds_fasta(local):
            if quals.get("gene"):
                by_gene.setdefault(quals["gene"], (quals.get("locus_tag", ""), sequence))
            if quals.get("locus_tag"):
                by_tag[quals["locus_tag"]] = (quals["locus_tag"], sequence)
        for gene, locus_tag, group in wanted:
            found = by_tag.get(locus_tag) or by_gene.get(gene)
            if found is None:
                raise HydraError(f"{accession} has neither {gene} nor {locus_tag}; "
                                 f"the record may have been re-annotated")
            name = _MARKER_RENAME.get((accession, gene), gene)
            records.append((f"{name}|{group}|{accession}", found[1]))
            provenance.append(f"{name}\t{group}\t{accession}\t{found[0]}\t{len(found[1])}")
        local.unlink(missing_ok=True)
        time.sleep(0.4)  # NCBI asks for no more than three requests a second
    for accession, name, group in _MARKER_WHOLE:
        local = work / f"{accession}.fna"
        _download(f"{_NCBI_EFETCH}?db=nuccore&id={accession}"
                  f"&rettype=fasta&retmode=text", local)
        sequences = [s for _h, s in read_fasta(local)]
        if not sequences:
            raise HydraError(f"{accession} returned no sequence")
        records.append((f"{name}|{group}|{accession}", sequences[0].upper()))
        provenance.append(f"{name}\t{group}\t{accession}\t-\t{len(sequences[0])}")
        local.unlink(missing_ok=True)
        time.sleep(0.4)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(out_fasta, records, wrap=0)
    (out_fasta.parent / "markers.sources.tsv").write_text(
        "\n".join(provenance) + "\n", encoding="utf-8")
    _makeblastdb(out_fasta)


def _stage_mlst(root: Path, work: Path) -> None:
    """Install the 7-locus Achtman scheme: alleles plus the profile table."""
    out_fasta = root / "mlst" / "mlst.fna"
    out_profiles = root / "mlst" / "profiles.tsv"
    if _built(out_fasta, out_profiles):
        return
    records: list[tuple[str, str]] = []
    for locus in _MLST_LOCI:
        local = work / f"{locus}.tfa"
        _download(f"{_PUBMLST_BASE}/loci/{locus}/alleles_fasta", local)
        for header, sequence in read_fasta(local):
            name = header.split()[0]
            allele = name[len(locus):].lstrip("_-") if name.startswith(locus) else name
            records.append((f"{locus}|{_clean_id(allele)}", sequence.upper()))
        local.unlink(missing_ok=True)
        time.sleep(0.5)  # PubMLST throttles
    if not records:
        raise HydraError("PubMLST returned no Salmonella MLST alleles")
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(out_fasta, records, wrap=0)
    _download(f"{_PUBMLST_BASE}/schemes/2/profiles_csv", out_profiles)
    _makeblastdb(out_fasta)


def fetch(db_dir: Path) -> None:
    """Download and build this module's data into ``db_dir/salmonella/``.

    Idempotent: each part is skipped when its output is already in place, so a
    failed download can be retried without starting over.
    """
    root = Path(db_dir) / KEY
    root.mkdir(parents=True, exist_ok=True)
    with tempdir(prefix="hydra.salmonella.", parent=root) as work:
        _stage_antigens(root, work)
        _stage_spi(root, work)
        _stage_markers(root, work)
        _stage_mlst(root, work)
    (root / "provenance.json").write_text(
        json.dumps({
            "key": KEY,
            "built": time.strftime("%Y-%m-%d"),
            "databases": [
                {"name": d.name, "source": d.source, "version": d.version,
                 "licence": d.licence, "redistribution": d.redistribution,
                 "reference": d.reference, "notes": d.notes}
                for d in DATABASES
            ],
        }, indent=1) + "\n",
        encoding="utf-8")
    LOG.info("Salmonella typing database ready in %s", root)


# ---------------------------------------------------------------------------
# Evidence helpers
# ---------------------------------------------------------------------------
@dataclass
class Evidence:
    """One locus recovered from the assembly, with everything it rests on."""

    name: str                 # antigen allele, marker gene or island
    label: str                # what it encodes: the antigen or island name
    reference: str            # the database sequence matched
    identity: float
    coverage: float
    covered_bases: int
    contig: str
    start: int
    end: int
    strand: str
    reference_span: str
    contig_length: int
    n_hsps: int

    @property
    def at_contig_edge(self) -> bool:
        return (self.start <= CONTIG_EDGE_BP
                or self.contig_length - self.end <= CONTIG_EDGE_BP)

    def describe(self) -> str:
        text = (f"{self.reference} {self.identity:.1f}% id, "
                f"{self.coverage:.1f}% cov ({self.reference_span}), "
                f"{self.contig}:{self.start}-{self.end}{self.strand}")
        if self.n_hsps > 1:
            text += f", {self.n_hsps} segments"
        if self.at_contig_edge:
            text += ", at contig edge"
        return text


def _evidence(hit: MergedHit, label: str, name: str) -> Evidence:
    return Evidence(
        name=name, label=label, reference=hit.sseqid,
        identity=hit.identity_pct, coverage=hit.coverage_pct,
        covered_bases=hit.covered_bases, contig=hit.qseqid,
        start=hit.query_start, end=hit.query_end, strand=hit.strand,
        reference_span=hit.coverage_string, contig_length=hit.qlen,
        n_hsps=hit.n_hsps,
    )


def _best(hits: list[MergedHit], label_of) -> dict[str, Evidence]:
    """Best hit per label, ranked by matched reference fraction then identity."""
    out: dict[str, Evidence] = {}
    for hit in hits:
        label, name = label_of(hit)
        if label is None:
            continue
        candidate = _evidence(hit, label, name)
        current = out.get(name)
        if current is None or (candidate.coverage * candidate.identity,
                               candidate.identity) > (current.coverage * current.identity,
                                                      current.identity):
            out[name] = candidate
    return out


def _split_id(sseqid: str, parts: int) -> list[str]:
    pieces = sseqid.split("|")
    while len(pieces) < parts:
        pieces.append("")
    return pieces


# ---------------------------------------------------------------------------
# The WKL serovar table
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SerovarRow:
    serovar: str
    o_antigen: str
    h1: str
    h2: str
    h_other: str
    serogroup: str
    h2_may_be_absent: bool
    subspecies: str


def _load_serovars(path: Path) -> list[SerovarRow]:
    rows: list[SerovarRow] = []
    with open(path, newline="", encoding="utf-8", errors="replace") as handle:
        for raw in csv.DictReader(handle):
            rows.append(SerovarRow(
                serovar=(raw.get("Serovar") or "").strip(),
                o_antigen=(raw.get("O_antigen") or "").strip(),
                h1=(raw.get("H1") or "").strip(),
                h2=(raw.get("H2") or "").strip(),
                h_other=(raw.get("H_other") or "").strip(),
                serogroup=(raw.get("Serogroup") or "").strip(),
                h2_may_be_absent=(raw.get("can_h2_be_missing") or "").strip().upper()
                == "TRUE",
                subspecies=(raw.get("subspecies") or "").strip(),
            ))
    if not rows:
        raise HydraError("the Salmonella serovar table is empty")
    return rows


def _match_serovars(rows: list[SerovarRow], serogroup: str | None,
                    h1: str | None, h2: str | None) -> list[SerovarRow]:
    """Rows of the WKL table whose formula the observed antigens fit.

    ``h2 is None`` means *not determined* and matches any row; ``h2 == "-"``
    means the second phase is genuinely absent and matches only rows that are
    monophasic or that the scheme marks as able to lose H2.
    """
    out = []
    for row in rows:
        if serogroup is not None and row.serogroup != serogroup:
            continue
        if h1 is not None and row.h1 != h1:
            continue
        if h2 is not None:
            if h2 == "-":
                if not (row.h2 == "-" or row.h2_may_be_absent):
                    continue
            elif row.h2 != h2:
                continue
        out.append(row)
    return out


def _formula(o_antigen: str, h1: str, h2: str, subspecies: str = "") -> str:
    prefix = f"{_SUBSPECIES_ROMAN.get(subspecies, '')} " if subspecies else ""
    return f"{prefix}{o_antigen}:{h1}:{h2}".strip()


def _similar_o_groups(serogroup: str | None) -> tuple[str, ...]:
    """O groups this one cannot be cleanly separated from by wzx/wzy."""
    if not serogroup:
        return ()
    for group in O_GROUP_SIMILARITY:
        if serogroup in group:
            return tuple(g for g in group if g != serogroup)
    return ()


def _o_antigen_for(rows: list[SerovarRow], serogroup: str) -> str:
    """The commonest O antigen string written for a serogroup in the table."""
    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        if row.serogroup == serogroup and row.o_antigen:
            counts[row.o_antigen] += 1
    if not counts:
        return serogroup
    return max(counts.items(), key=lambda kv: (kv[1], -len(kv[0])))[0]


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------
def _require(path: Path, what: str) -> Path:
    if not path.exists():
        raise HydraError(f"{what} is missing from the Salmonella typing database "
                         f"({path}); run the module's fetch() first")
    return path


def run(assembly: Path, db_dir: Path, *, threads: int = 1,
        workdir: Path | None = None) -> list[TypingResult]:
    """Type one *Salmonella* assembly. One TypingResult per scheme."""
    assembly = Path(assembly)
    root = Path(db_dir) / KEY
    antigen_db = _require(root / "antigens" / "antigens.fna", "the antigen allele set")
    serovar_table = _require(root / "antigens" / "serovars.csv", "the serovar table")
    marker_db = _require(root / "markers" / "markers.fna", "the marker gene set")
    spi_db = root / "spi" / "spi.fna"
    mlst_db = root / "mlst" / "mlst.fna"
    rows = _load_serovars(serovar_table)

    if workdir is None:
        with tempdir(prefix="hydra.salmonella.run.") as scratch:
            return _run_in(assembly, root, antigen_db, marker_db, spi_db, mlst_db,
                           rows, scratch, threads)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    return _run_in(assembly, root, antigen_db, marker_db, spi_db, mlst_db, rows,
                   workdir, threads)


def _plain_fasta(assembly: Path, work: Path) -> Path:
    """A path blastn can read: compressed input is expanded into *work* first.

    BLAST+ cannot read a gzipped FASTA, and an assembly handed straight from a
    directory of ``.fna.gz`` files would otherwise fail with a BLAST parse error
    that says nothing about compression.
    """
    if assembly.suffix.lower() not in (".gz", ".bz2", ".xz"):
        return assembly
    plain = work / (assembly.stem if assembly.stem.endswith((".fna", ".fa", ".fasta"))
                    else assembly.stem + ".fna")
    if not plain.exists():
        with smart_open(assembly, "rt") as src, open(plain, "w", encoding="utf-8") as out:
            shutil.copyfileobj(src, out)
    return plain


def _run_in(assembly: Path, root: Path, antigen_db: Path, marker_db: Path,
            spi_db: Path, mlst_db: Path, rows: list[SerovarRow], work: Path,
            threads: int) -> list[TypingResult]:
    assembly = _plain_fasta(assembly, work)
    antigens = _blast_antigens(assembly, antigen_db, work, threads)
    markers = _blast_markers(assembly, marker_db, work, threads)
    serovar, subspecies = _call_serovar(rows, antigens, markers)
    results = [serovar, subspecies, _call_h2_locus(antigens, markers)]
    if spi_db.exists():
        results.append(_call_spi(assembly, spi_db, work, threads))
    if mlst_db.exists():
        results.append(_call_mlst(assembly, mlst_db, root / "mlst" / "profiles.tsv",
                                  work, threads, serovar))
    results.append(_call_virulence_plasmid(markers))
    results.append(_call_typhoidal(markers, antigens))
    return results


# ------------------------------------------------------------------ antigens
@dataclass
class AntigenCalls:
    serogroup: str | None = None
    serogroup_note: str = ""
    h1: str | None = None
    h2: str | None = None
    h1_note: str = ""
    h2_note: str = ""
    #: locus -> best Evidence, locus -> confident|provisional, locus -> runners-up
    evidence: dict = field(default_factory=dict)
    confidence: dict = field(default_factory=dict)
    alternatives: dict = field(default_factory=dict)


def _blast_antigens(assembly: Path, db: Path, work: Path,
                    threads: int) -> AntigenCalls:
    hsps = blast("blastn", assembly, db, work / "salmonella.antigens.tsv",
                 threads=threads, evalue=1e-20, task="megablast", word_size=20,
                 perc_identity=80.0, max_target_seqs=20000)
    merged = merge_hsps(hsps)
    per_locus: dict[str, list[MergedHit]] = defaultdict(list)
    for hit in merged:
        locus = _split_id(hit.sseqid, 3)[0]
        per_locus[locus].append(hit)

    calls = AntigenCalls(evidence={}, confidence={}, alternatives={})
    for locus in ("wzx", "wzy"):
        best, alts = _pick_o(per_locus.get(locus, []))
        if best is not None:
            calls.evidence[locus] = best
            calls.alternatives[locus] = alts
    for locus in ("fliC", "fljB"):
        best, level, alts = _pick_h(per_locus.get(locus, []))
        if best is not None:
            calls.evidence[locus] = best
            calls.confidence[locus] = level
            calls.alternatives[locus] = alts

    _combine_serogroup(calls)
    _combine_h(calls, "fliC", "h1")
    _combine_h(calls, "fljB", "h2")
    return calls


def _runners_up(pool: list[Evidence], best: Evidence, limit: int = 3) -> list[str]:
    """Other antigens the same locus also matched, with the margin shown.

    The margin is the point: an fljB locus that matches 1,2 at 100% and 1,5 at
    99.3% is not the same finding as one that matches both at 100%, and a reader
    deciding whether to confirm by serology needs to see which it was.
    """
    seen: dict[str, Evidence] = {}
    for ev in pool:
        if ev.label == best.label:
            continue
        current = seen.get(ev.label)
        if current is None or ev.identity > current.identity:
            seen[ev.label] = ev
    ordered = sorted(seen.values(), key=lambda e: e.identity, reverse=True)
    return [f"{e.label} at {e.identity:.1f}% identity over {e.covered_bases} bases"
            for e in ordered[:limit]]


def _pick_o(hits: list[MergedHit]) -> tuple[Evidence | None, list[str]]:
    """Best serogroup hit for wzx or wzy, with the rules SISTR publishes."""
    kept: list[Evidence] = []
    for hit in hits:
        antigen = _split_id(hit.sseqid, 3)[2]
        if not antigen:
            continue
        if hit.identity_pct < O_MIN_IDENTITY or hit.covered_bases < O_MIN_BASES:
            continue
        if hit.covered_bases < O_SHORT_BASES and hit.identity_pct < O_SHORT_IDENTITY:
            continue
        kept.append(_evidence(hit, antigen, antigen))
    if not kept:
        return None, []
    kept.sort(key=lambda e: (e.coverage >= 99.9 and e.identity >= 99.9,
                             e.coverage * e.identity, e.identity), reverse=True)
    best = kept[0]
    return best, _runners_up(kept, best)


def _pick_h(hits: list[MergedHit]) -> tuple[Evidence | None, str, list[str]]:
    """Best fliC or fljB allele, and how much weight the call can carry."""
    confident: list[Evidence] = []
    provisional: list[Evidence] = []
    for hit in hits:
        antigen = _split_id(hit.sseqid, 3)[2]
        if not antigen:
            continue
        ev = _evidence(hit, antigen, antigen)
        if (hit.covered_bases >= H_CONFIDENT_BASES
                and hit.identity_pct >= H_CONFIDENT_IDENTITY):
            confident.append(ev)
        elif (hit.covered_bases >= H_PROVISIONAL_BASES
              and hit.identity_pct >= H_PROVISIONAL_IDENTITY):
            provisional.append(ev)
    pool, level = (confident, "confident") if confident else (provisional, "provisional")
    if not pool:
        return None, "none", []
    pool.sort(key=lambda e: (e.identity, e.coverage, e.covered_bases), reverse=True)
    best = pool[0]
    return best, level, _runners_up(pool, best)


def _combine_serogroup(calls: AntigenCalls) -> None:
    wzx = calls.evidence.get("wzx")
    wzy = calls.evidence.get("wzy")
    if wzx is None and wzy is None:
        calls.serogroup_note = ("neither wzx nor wzy reached 88% identity over "
                                "300 bases, so the O serogroup is not determined")
        return
    if wzx is not None and wzy is not None:
        if wzx.label == wzy.label:
            calls.serogroup = wzx.label
            return
        stronger, weaker = ((wzx, wzy) if wzx.coverage * wzx.identity
                            >= wzy.coverage * wzy.identity else (wzy, wzx))
        calls.serogroup = stronger.label
        calls.serogroup_note = (
            f"wzx and wzy disagree (wzx O{wzx.label}, wzy O{wzy.label}); the "
            f"better-matched gene ({'wzx' if stronger is wzx else 'wzy'}) was "
            f"taken, but O{weaker.label} cannot be excluded")
        return
    only = wzx if wzx is not None else wzy
    which = "wzx" if wzx is not None else "wzy"
    other = "wzy" if wzx is not None else "wzx"
    calls.serogroup = only.label
    calls.serogroup_note = (f"serogroup called from {which} alone; {other} was "
                            f"not recovered")


def _combine_h(calls: AntigenCalls, locus: str, field: str) -> None:
    ev = calls.evidence.get(locus)
    if ev is None:
        setattr(calls, field, None)
        return
    setattr(calls, field, ev.label)
    notes = []
    if calls.confidence.get(locus) == "provisional":
        notes.append(f"{locus} matched only {ev.covered_bases} bases at "
                     f"{ev.identity:.1f}% identity, below the threshold for a "
                     f"confident allele call")
    if ev.coverage < 99.0 and ev.at_contig_edge:
        notes.append(f"{locus} is incomplete and runs to a contig boundary, so "
                     f"the allele may be truncated by the assembly")
    alts = calls.alternatives.get(locus) or []
    if alts:
        notes.append(f"{locus} also matched {', '.join(alts)}")
    setattr(calls, f"{field}_note", "; ".join(notes))


# ------------------------------------------------------------------- markers
def _blast_markers(assembly: Path, db: Path, work: Path,
                   threads: int) -> dict[str, Evidence]:
    hsps = blast("blastn", assembly, db, work / "salmonella.markers.tsv",
                 threads=threads, evalue=1e-10, task="megablast", word_size=16,
                 perc_identity=80.0, max_target_seqs=5000)
    merged = merge_hsps(hsps)
    # fliC, fljB and the pBSSB1 z66 allele are homologous flagellins that share
    # their conserved termini, so the fljB reference matches the fliC locus over
    # 50-60% of its length at ~92% identity. Left alone that makes fljB look
    # "partial" in a genome that has no fljB at all, which is precisely the
    # wrong answer for a monophasic isolate. One place on the assembly can only
    # belong to one gene: keep the better-matching marker there.
    merged = deduplicate(
        merged,
        key_span=lambda h: (h.query_start, h.query_end),
        key_seq=lambda h: h.qseqid,
        key_score=lambda h: h.match_score,
        key_tiebreak=lambda h: h.sseqid,
    )

    def label_of(hit: MergedHit):
        name, group, _acc = _split_id(hit.sseqid, 3)
        if not name or hit.identity_pct < MARKER_IDENTITY:
            return None, None
        if hit.coverage_pct < MARKER_PARTIAL_COVERAGE:
            return None, None
        return group, name

    return _best(merged, label_of)


def _marker_state(markers: dict[str, Evidence], name: str) -> str:
    ev = markers.get(name)
    if ev is None:
        return "absent"
    if ev.coverage >= MARKER_COVERAGE:
        return "present"
    return "partial"


def _marker_rows(markers: dict[str, Evidence], names) -> dict[str, str]:
    out: dict[str, str] = {}
    for name in names:
        ev = markers.get(name)
        out[name] = ev.describe() if ev else "not detected"
    return out


# ------------------------------------------------------------------- serovar
def _call_serovar(rows: list[SerovarRow], antigens: AntigenCalls,
                  markers: dict[str, Evidence]) -> tuple[TypingResult, TypingResult]:
    serogroup = antigens.serogroup
    h1 = antigens.h1
    h2_present = antigens.h2
    h2_locus_absent = all(_marker_state(markers, g) == "absent"
                          for g in ("fljB", "fljA", "hin"))
    fljB_edge = (antigens.evidence.get("fljB") is not None
                 and antigens.evidence["fljB"].at_contig_edge
                 and antigens.evidence["fljB"].coverage < 99.0)

    notes: list[str] = []
    # Decide what the second phase actually is: present, genuinely absent, or
    # undetermined. A bare negative is never reported.
    if h2_present:
        h2: str | None = h2_present
    elif h2_locus_absent:
        h2 = "-"
        notes.append("fljB, fljA and hin are all absent, so the second-phase "
                     "flagellar operon is deleted and H2 is genuinely negative")
    else:
        h2 = None
        states = ", ".join(f"{g} {_marker_state(markers, g)}"
                           for g in ("fljB", "fljA", "hin"))
        notes.append(f"no fljB allele was called but the fljBA-hin region is not "
                     f"cleanly absent ({states}); H2 is not determined rather "
                     f"than negative")
    if fljB_edge:
        notes.append("the fljB match stops at a contig boundary")

    # Which rows of the WKL table the formula fits. With neither an O group nor
    # an H1 allele there is nothing to look up: an absent fljAB region on its
    # own would otherwise match every monophasic serovar in the scheme and
    # produce a list of names out of no evidence at all.
    matches: list[SerovarRow] = []
    widened = ""
    if serogroup is not None or h1 is not None:
        matches = _match_serovars(rows, serogroup, h1, h2)
        if not matches and h2 == "-":
            # A monophasic variant of a biphasic serovar: the scheme may have no
            # row for it, so look the biphasic parents up and name the variant
            # explicitly rather than reporting nothing.
            parents = _match_serovars(rows, serogroup, h1, None)
            if parents:
                widened, matches = "monophasic", parents
        if not matches and h2 is not None and h2 != "-":
            relaxed = _match_serovars(rows, serogroup, h1, None)
            if relaxed:
                widened, matches = "h2_ignored", relaxed
        if not matches and h1 is not None:
            relaxed = _match_serovars(rows, serogroup, None, None)
            if relaxed:
                widened, matches = "h1_ignored", relaxed

    o_antigen = (_o_antigen_for(rows, serogroup) if serogroup else "-")
    if matches and widened in ("", "monophasic", "h2_ignored"):
        # Prefer the exact O antigen string the matched row(s) are written with:
        # "1,4,[5],12" is what a laboratory writes, not the commonest string in
        # the whole O group. Only when the rows disagree does the group-level
        # string stand in.
        strings = sorted({row.o_antigen for row in matches if row.o_antigen})
        if len(strings) == 1:
            o_antigen = strings[0]
        elif strings:
            # The candidate serovars are written with different O factors -
            # Newport is 6,8,20 and Bardo, which shares its H antigens, is 8.
            # wzx/wzy call the group, not the individual factors, so both
            # strings are reported rather than one being chosen.
            o_antigen = " or ".join(strings[:3])
    vi = _marker_state(markers, "tviB") == "present" and \
        _marker_state(markers, "vexC") == "present"
    if vi and "Vi" not in o_antigen:
        o_antigen = f"{o_antigen},Vi"

    subspecies_set = sorted({row.subspecies for row in matches if row.subspecies})
    subspecies = subspecies_set[0] if len(subspecies_set) == 1 else ""

    formula = _formula(o_antigen, h1 or "-", h2 if h2 is not None else "?",
                       subspecies)
    names = [row.serovar for row in matches]

    call, lineage, confidence = _name_serovar(names, widened, h2, h1, serogroup,
                                              formula, notes)
    _serovar_caveats(notes, rows, antigens, matches, names, markers, h2, widened)

    alleles: dict[str, str] = {"antigenic_formula": formula}
    for locus in ("wzx", "wzy", "fliC", "fljB"):
        ev = antigens.evidence.get(locus)
        alleles[locus] = ev.describe() if ev else "not detected"
    alleles["O_serogroup"] = (f"serogroup {serogroup}" if serogroup
                              else "not determined")
    alleles["H1_fliC"] = _h_with_confidence(antigens, "fliC", h1)
    alleles["H2_fljB"] = (_h_with_confidence(antigens, "fljB", h2)
                          if h2 is not None else "not determined")
    alleles["Vi_capsule"] = ("viaB present" if vi else "viaB not detected")
    alleles["sdf-I"] = _marker_rows(markers, ("sdf-I",))["sdf-I"]
    alleles["invA"] = _marker_rows(markers, ("invA",))["invA"]
    alleles["oafA"] = _marker_rows(markers, ("oafA",))["oafA"]
    oaf = markers.get("oafA")
    if oaf is not None and oaf.coverage < 99.0:
        notes.append(f"oafA is present but only {oaf.coverage:.1f}% covered. A "
                     f"short internal deletion in oafA is what makes a "
                     f"Typhimurium O:5-negative (var. Copenhagen), so this is "
                     f"worth looking at; this module reports the evidence and "
                     f"does not make an O:5 call, because that branch was not "
                     f"validated here")
    if len(names) > 1:
        alleles["serovars_sharing_this_formula"] = ", ".join(names[:12]) + (
            ", ..." if len(names) > 12 else "")
    if h2 == "-" and serogroup and h1 and h1 != "-":
        # Only a serovar the scheme has no *name* for is a "variant" of
        # something: Enteritidis and Typhi are monophasic in their own right and
        # must not be described as variants of a biphasic relative. The scheme
        # writes a nameless one as a bare formula, "I 1,4,[5],12:i:-".
        nameless = any(":" in row.serovar and row.serovar.endswith(":-")
                       for row in matches)
        subspecies_here = {row.subspecies for row in matches}
        parents = sorted({f"{row.serovar} (H2 {row.h2})" for row
                          in _match_serovars(rows, serogroup, h1, None)
                          if row.h2 not in ("-", "")
                          and row.subspecies in subspecies_here
                          and (not o_antigen or row.o_antigen == o_antigen)})
        if parents and (nameless or widened == "monophasic"):
            alleles["biphasic_counterparts"] = ", ".join(parents[:8])
            text = (f"the second phase is absent while {', '.join(parents[:4])} "
                    f"share this O group and H1 and are biphasic, so this is a "
                    f"monophasic variant of that set")
            if serogroup == "B" and h1 == "i":
                text += (". Monophasic 4,[5],12:i:- has arisen independently "
                         "several times, so the antigens do not identify the "
                         "lineage: read the ST - ST34 is the widespread "
                         "epidemic clone, ST19 the ordinary Typhimurium "
                         "background")
            else:
                text += (". Read the ST beside it: a monophasic formula says "
                         "nothing about which lineage lost its second phase")
            notes.append(text)
    siblings = _similar_o_groups(serogroup)
    if siblings:
        extra: list[str] = []
        for group in siblings:
            for row in _match_serovars(rows, group, h1, h2):
                extra.append(f"{row.serovar} (O group {group})")
        alleles["o_group_ambiguity"] = (
            f"wzx/wzy cannot reliably separate O group(s) "
            f"{', '.join(siblings)} from {serogroup}"
            + (f"; were the O group one of those, the same H antigens would fit "
               f"{', '.join(extra[:8])}" if extra else
               "; no serovar with these H antigens sits in those groups"))
        if extra:
            notes.append(
                f"O groups {', '.join((serogroup,) + tuple(siblings))} share "
                f"closely related wzx/wzy alleles, so the O group rests on a "
                f"gene call that cannot fully separate them; "
                f"{', '.join(extra[:6])} would fit the same H antigens and must "
                f"be excluded by agglutination")

    loci_found = sum(1 for locus in ("wzx", "wzy", "fliC", "fljB")
                     if antigens.evidence.get(locus) is not None)
    serovar_result = TypingResult(
        scheme="salmonella_serovar", call=call, lineage=lineage,
        alleles=alleles, loci_found=loci_found, loci_total=4,
        note="; ".join(n for n in notes if n),
    )
    subspecies_result = _call_subspecies(matches, subspecies_set, formula,
                                         serogroup, h1, h2)
    serovar_result.alleles["confidence"] = confidence
    return serovar_result, subspecies_result


def _h_with_confidence(antigens: AntigenCalls, locus: str,
                       value: str | None) -> str:
    """An H antigen with how much the allele behind it can be trusted.

    A flagellin's termini are conserved across antigen types, so an allele
    recovered only in part - especially one cut off by a contig end - can be
    the conserved half of a different antigen. The value is still reported,
    because losing it would be worse, but it is labelled.
    """
    if value is None:
        return "not determined"
    if value == "-":
        return "-"
    ev = antigens.evidence.get(locus)
    if ev is None:
        return value
    if antigens.confidence.get(locus) == "provisional" or ev.coverage < 99.0:
        tail = " (provisional: " + (
            f"only {ev.coverage:.0f}% of the allele matched"
            + (", and the match runs to a contig boundary" if ev.at_contig_edge
               else "")) + ")"
        return value + tail
    return value


def _name_serovar(names: list[str], widened: str, h2, h1, serogroup,
                  formula: str, notes: list[str]) -> tuple[str, str, str]:
    """Turn the matched rows into a reported name, lineage field and confidence."""
    if not names:
        if serogroup is None and h1 is None:
            notes.append("no O or H antigen gene was recovered; this assembly "
                         "may not be Salmonella, or the antigen loci may be "
                         "missing from it")
            return "not determined", "-", "none"
        notes.append("the formula found is not in the White-Kauffmann-Le Minor "
                     "table; it is reported as a formula, which is what a "
                     "reference laboratory would do with a new antigen "
                     "combination")
        return formula, "-", "formula only"
    if widened == "monophasic":
        parents = ", ".join(names[:4])
        notes.append(f"the second phase is absent, which the scheme has no row "
                     f"for: this is a monophasic variant of {parents}")
        if len(names) == 1:
            return f"{names[0]} monophasic variant", names[0], "monophasic variant"
        return f"monophasic variant of {parents}", parents, "monophasic variant"
    if widened == "h2_ignored":
        notes.append(f"H2 {h2} does not occur with this O group and H1 in the "
                     f"scheme; the serovar below is the match ignoring H2, which "
                     f"needs confirming by serology")
        return "|".join(names[:4]), (names[0] if names else "-"), "H2 unmatched"
    if widened == "h1_ignored":
        notes.append("the H1 allele found does not occur in this O group in the "
                     "scheme; only the O group could be matched")
        return formula, "-", "O group only"
    if len(names) == 1:
        return names[0], names[0], "antigens agree with one serovar"
    return "|".join(names), names[0], f"{len(names)} serovars share this formula"


def _serovar_caveats(notes: list[str], rows: list[SerovarRow],
                     antigens: AntigenCalls, matches: list[SerovarRow],
                     names: list[str], markers: dict[str, Evidence],
                     h2, widened: str) -> None:
    """Add the caveats the literature says every call needs."""
    for text in (antigens.serogroup_note, antigens.h1_note, antigens.h2_note):
        if text:
            notes.append(text)
    if len(names) > 1 and widened == "":
        notes.append(f"the antigens cannot separate {', '.join(names[:6])}: they "
                     f"share one formula in the scheme and are told apart "
                     f"biochemically or by lineage, not by agglutination")
    # A non-motile serovar still carries the flagellin genes of its motile
    # relatives, so gene content cannot exclude it. Gallinarum and Pullorum
    # sitting behind an Enteritidis call is the textbook case, and SeqSero2 adds
    # the Sdf-I fragment precisely to break that tie.
    if antigens.h1 and antigens.h1 != "-" and antigens.serogroup:
        siblings = sorted({row.serovar for row in
                           _match_serovars(rows, antigens.serogroup, "-", "-")})
        if siblings:
            text = (f"motility is a phenotype: the non-motile serovar(s) "
                    f"{', '.join(siblings[:4])} of O group {antigens.serogroup} "
                    f"carry the same flagellin genes as the call above and "
                    f"cannot be excluded from an assembly, only from a motility "
                    f"and agglutination result")
            if {"Gallinarum", "Pullorum"} & set(siblings):
                text += (f". Sdf-I is {_marker_state(markers, 'sdf-I')}: present "
                         f"supports commonly circulating Enteritidis, absent "
                         f"leaves Gallinarum and Pullorum open")
            notes.append(text)
    if h2 is None:
        notes.append("because H2 is undetermined, the serovar list above is the "
                     "set consistent with any second phase")


def _call_subspecies(matches: list[SerovarRow], subspecies_set: list[str],
                     formula: str, serogroup, h1, h2) -> TypingResult:
    alleles = {
        "antigenic_formula": formula,
        "method": ("the subspecies of the White-Kauffmann-Le Minor rows the "
                   "observed formula fits; not an independent genomic marker"),
    }
    if not matches:
        return TypingResult(
            scheme="salmonella_subspecies", call="not determined", lineage="-",
            alleles=alleles, loci_found=0, loci_total=1,
            note="the formula matched no row of the scheme, so no subspecies "
                 "could be read from it. Subspecies here is inferred from the "
                 "antigenic formula; a genome whose formula is unlisted needs "
                 "an independent method (ANI, cgMLST or the biochemical key).",
        )
    if len(subspecies_set) == 1:
        name = subspecies_set[0]
        roman = _SUBSPECIES_ROMAN.get(name, "?")
        label = ("Salmonella bongori" if name == "bongori"
                 else f"Salmonella enterica subsp. {name}")
        alleles["subspecies_rows"] = f"{len(matches)} row(s), all subsp. {name}"
        note = (f"{label} (group {roman}). Inferred from the scheme: every "
                f"serovar with this formula is {name}. This is not a "
                f"subspecies-specific marker, so a formula shared across "
                f"subspecies would be reported as undetermined.")
        if name == "bongori":
            note += (" S. bongori is now a separate species; V is its historic "
                     "subspecies designation.")
        return TypingResult(
            scheme="salmonella_subspecies", call=f"{roman} ({name})",
            lineage=label, alleles=alleles, loci_found=1, loci_total=1, note=note)
    romans = ", ".join(f"{_SUBSPECIES_ROMAN.get(s, '?')} ({s})"
                       for s in subspecies_set)
    alleles["subspecies_rows"] = f"{len(matches)} row(s) across {len(subspecies_set)} subspecies"
    return TypingResult(
        scheme="salmonella_subspecies", call="not determined",
        lineage=romans, alleles=alleles, loci_found=1, loci_total=1,
        note=f"the formula occurs in more than one subspecies in the scheme "
             f"({romans}), so it does not decide the subspecies. Separate them "
             f"by ANI, cgMLST or the biochemical key.",
    )


# ------------------------------------------------------------------ H2 locus
def _call_h2_locus(antigens: AntigenCalls,
                   markers: dict[str, Evidence]) -> TypingResult:
    """Report the state of the fljBA-hin region, not just whether fljB was seen.

    An absent fljB *call* and an absent fljB *gene* are different findings, and
    the reference multiplex interrogates the region rather than the gene for
    exactly that reason (Prendergast2013). The widespread monophasic clone has
    the region physically replaced (Garcia2016), so fljA and hin go with fljB.
    """
    states = {gene: _marker_state(markers, gene)
              for gene in ("fljB", "fljA", "hin", "iroB")}
    alleles = _marker_rows(markers, ("fljB", "fljA", "hin", "iroB", "fliC"))
    alleles["fljB_antigen_allele"] = (
        antigens.evidence["fljB"].describe() if antigens.evidence.get("fljB")
        else "no fljB antigen allele called")
    present = [g for g in ("fljB", "fljA", "hin") if states[g] == "present"]
    found = len(present)
    all_absent = all(states[g] == "absent" for g in ("fljB", "fljA", "hin"))

    if found == 3:
        call = "intact"
        note = ("fljB, fljA and hin are all present: the second-phase flagellar "
                "operon is intact and a biphasic serovar is expected")
    elif all_absent:
        call = "fljBA-hin absent"
        note = ("fljB, fljA and hin are all absent: there is no second-phase "
                "flagellar operon, so H2 is genuinely negative. That is normal "
                "for the many serovars the scheme lists as monophasic (Typhi, "
                "Enteritidis, Dublin, Gallinarum, Senftenberg and others); it is "
                "only a monophasic *variant* when the scheme entry for the "
                "serovar is biphasic. In the widespread 4,[5],12:i:- clone the "
                "region is replaced by an acquired resistance module, so check "
                "the AMR table for the ASSuT genes (blaTEM-1, strA/strB, sul2, "
                "tet(B)) and the MLST call for ST34")
    else:
        call = "partial"
        missing = [g for g in ("fljB", "fljA", "hin") if states[g] != "present"]
        note = (f"the region is incomplete: {', '.join(missing)} not fully "
                f"recovered while {', '.join(present)} is present. An absent "
                f"fljB with fljA or hin still there is not a clean deletion - it "
                f"is a disrupted gene, a short contig, or an assembly gap, and "
                f"H2 should be confirmed by serology or by inspecting the region")
    edges = [g for g in ("fljB", "fljA", "hin")
             if markers.get(g) is not None and markers[g].at_contig_edge]
    if edges:
        note += (f". {', '.join(edges)} sits within {CONTIG_EDGE_BP} bp of a "
                 f"contig boundary, so absence of its neighbours may be an "
                 f"assembly artefact")
    if _marker_state(markers, "fljB_z66") == "present":
        note += (". The pBSSB1 z66 fljB allele is present: this is the H:z66 "
                 "flagellar variant carried on a linear plasmid, seen in "
                 "Indonesian S. Typhi")
    return TypingResult(
        scheme="salmonella_h2_locus", call=call,
        lineage=antigens.h2 or ("-" if found == 0 else "not determined"),
        alleles=alleles, loci_found=found, loci_total=3, note=note,
    )


# ----------------------------------------------------------------------- SPI
def _call_spi(assembly: Path, db: Path, work: Path, threads: int) -> TypingResult:
    hsps = blast("blastn", assembly, db, work / "salmonella.spi.tsv",
                 threads=threads, evalue=1e-20, task="megablast", word_size=28,
                 perc_identity=80.0, max_target_seqs=5000)
    merged = merge_hsps(hsps)
    best: dict[str, Evidence] = {}
    for hit in merged:
        island, _variant, accession = _split_id(hit.sseqid, 3)
        if not island or hit.identity_pct < SPI_IDENTITY:
            continue
        ev = _evidence(hit, island, island)
        ev.reference = f"{hit.sseqid} ({accession})"
        current = best.get(island)
        if current is None or ev.coverage > current.coverage:
            best[island] = ev
    catalogue = _spi_catalogue(db)
    named = {k: e for k, e in best.items() if k not in _SPI_UNNAMED}
    present = sorted(k for k, e in named.items() if e.coverage >= SPI_PRESENT_COVERAGE)
    partial = sorted(k for k, e in named.items()
                     if SPI_PARTIAL_COVERAGE <= e.coverage < SPI_PARTIAL_CEILING)
    alleles: dict[str, str] = {}
    for island, ev in sorted(best.items(), key=lambda kv: _island_key(kv[0])):
        if ev.coverage >= SPI_PRESENT_COVERAGE:
            status = "present"
        elif ev.coverage >= SPI_PARTIAL_COVERAGE:
            status = "partial"
        elif ev.coverage >= SPI_TRACE_COVERAGE:
            status = "trace only - treat as absent"
        else:
            continue
        alleles[island] = (f"{status}: {ev.coverage:.0f}% of {ev.reference} at "
                           f"{ev.identity:.1f}% identity, {ev.contig}:"
                           f"{ev.start}-{ev.end}")
    for island in sorted(catalogue - set(alleles), key=_island_key):
        alleles[island] = "absent: no hit above 85% identity"
    alleles["interpretation"] = (
        "SPI-1 (invasion T3SS-1) and SPI-2 (intracellular survival T3SS-2) are "
        "core to S. enterica and their absence is a quality-control finding, not "
        "a serovar finding. SPI-7 carries the viaB Vi capsule locus and is "
        "restricted to Typhi, Paratyphi C and Dublin. SGI1 is the "
        "multidrug-resistance island of DT104 and of some Newport and "
        "Kentucky lineages.")
    call = ", ".join(present) if present else "none at full coverage"
    note = (f"{len(present)} island(s) at >={SPI_PRESENT_COVERAGE:.0f}% coverage"
            + (f"; {len(partial)} partial ({', '.join(partial)}) - islands are "
               f"mosaic and a partial hit usually means a divergent or "
               f"rearranged copy, not a fragment" if partial else "")
            + ". Coverage is against the SPIFinder reference for each island, "
              "some of which are themselves partial sequences, so a coverage "
              "figure is comparable between genomes but not an absolute "
              "measure of island length.")
    return TypingResult(scheme="salmonella_spi", call=call,
                        lineage=", ".join(partial) if partial else "-",
                        alleles=alleles, loci_found=len(present),
                        loci_total=len(catalogue), note=note)


#: Entries in the SPIFinder database with no island name of their own; they are
#: kept as evidence but never named in a headline call.
_SPI_UNNAMED = frozenset({"Not_named", "Not"})
_SPI_CATALOGUE: dict[str, frozenset] = {}


def _spi_catalogue(db: Path) -> frozenset:
    """Every island name in the installed SPIFinder database.

    The denominator has to be the database's island count, not the number this
    genome happened to hit, or ``loci_found/loci_total`` would mean something
    different for every sample.
    """
    key = str(db)
    cached = _SPI_CATALOGUE.get(key)
    if cached is None:
        names = set()
        with open(db, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    island = line[1:].split("|")[0].strip()
                    if island and island not in _SPI_UNNAMED:
                        names.add(island)
        cached = frozenset(names)
        _SPI_CATALOGUE[key] = cached
    return cached


def _island_key(name: str) -> tuple:
    match = re.match(r"SPI-(\d+)", name)
    return (0, int(match.group(1)), "") if match else (1, 0, name)


# ----------------------------------------------------------- virulence plasmid
def _call_virulence_plasmid(markers: dict[str, Evidence]) -> TypingResult:
    spv = ("spvR", "spvA", "spvB", "spvC", "spvD")
    backbone = ("pefA", "rcK", "traT", "samA")
    alleles = _marker_rows(markers, spv + backbone)
    found = [g for g in spv if _marker_state(markers, g) == "present"]
    extra = [g for g in backbone if _marker_state(markers, g) == "present"]
    alleles["interpretation"] = _GROUP_MEANING["spv"]
    if len(found) == len(spv):
        call = "spv locus complete"
        note = ("the whole spvRABCD operon is present: this genome carries the "
                "Salmonella virulence plasmid (pSLT-like, IncFII). SpvB is an "
                "actin ADP-ribosyltransferase and SpvC a phosphothreonine lyase; "
                "the locus is associated with systemic infection in "
                "non-typhoidal serovars and is characteristic of Typhimurium, "
                "Enteritidis, Dublin, Choleraesuis, Abortusovis and Gallinarum/"
                "Pullorum")
    elif found:
        call = f"spv partial ({len(found)}/{len(spv)})"
        where = markers[found[0]]
        note = (f"only {', '.join(found)} of the spv operon was recovered, at "
                f"{where.contig}:{where.start} on a {where.contig_length:,} bp "
                f"sequence. Whether that is a virulence plasmid that did not "
                f"assemble cleanly or a chromosomal copy of the locus cannot be "
                f"decided from gene presence alone - read the coordinates above "
                f"against Hydra's plasmid replicon table before calling it a "
                f"plasmid")
    else:
        call = "spv absent"
        note = ("no spv gene was detected, so no Salmonella virulence plasmid of "
                "the pSLT type. Most serovars, including Typhi, Paratyphi A and "
                "Infantis, do not carry one, so this is an expected negative "
                "rather than a quality problem")
    if extra:
        note += f". pSLT backbone genes also present: {', '.join(extra)}"
    elif found:
        note += (". No pSLT backbone gene (pefA, rck, traT, samA) was found "
                 "beside the spv genes, which is a reason not to report a pSLT "
                 "type virulence plasmid on this evidence alone")
    return TypingResult(scheme="salmonella_virulence_plasmid", call=call,
                        lineage=", ".join(found) if found else "-",
                        alleles=alleles, loci_found=len(found),
                        loci_total=len(spv), note=note)


# ----------------------------------------------------------------- typhoidal
def _call_typhoidal(markers: dict[str, Evidence],
                    antigens: AntigenCalls) -> TypingResult:
    via_b = ("tviA", "tviB", "tviC", "tviD", "tviE",
             "vexA", "vexB", "vexC", "vexD", "vexE")
    toxin = ("cdtB", "pltA", "pltB")
    others = ("pilL", "sefA", "sopE", "fljB_z66")
    alleles = _marker_rows(markers, via_b + toxin + others)
    vi = [g for g in via_b if _marker_state(markers, g) == "present"]
    tox = [g for g in toxin if _marker_state(markers, g) == "present"]
    alleles["interpretation"] = (
        f"{_GROUP_MEANING['viaB']}; {_GROUP_MEANING['typhoid_toxin']}; pilL is "
        f"the SPI-7 type IVB pilus; sefA the SEF14 fimbrial subunit.")

    parts: list[str] = []
    notes: list[str] = []
    if len(vi) >= 8:
        parts.append("Vi capsule")
        notes.append(f"the viaB region is present ({len(vi)}/10 genes: "
                     f"{', '.join(vi)}): the isolate expresses or can express the "
                     f"Vi capsular polysaccharide, which puts O:Vi in the "
                     f"formula and matters for Vi-based vaccines and for the "
                     f"Widal/Vi serology. In S. enterica viaB sits on SPI-7 and "
                     f"is found in Typhi, Paratyphi C and Dublin")
    elif vi:
        notes.append(f"only {len(vi)}/10 viaB genes were found ({', '.join(vi)}); "
                     f"that is not a Vi locus. Partial recovery of a 14 kb "
                     f"region usually means a contig break, so Vi is reported as "
                     f"not determined rather than absent")
    else:
        notes.append("no viaB gene was detected: no Vi capsule locus")
    if len(tox) == 3:
        parts.append("typhoid toxin")
        notes.append("cdtB, pltA and pltB are all present: the genome encodes "
                     "typhoid toxin, the A2B5 genotoxin of the typhoidal "
                     "serovars and of some non-typhoidal ones such as Javiana "
                     "and Montevideo")
    elif tox:
        notes.append(f"typhoid toxin incomplete ({', '.join(tox)} of cdtB, pltA, "
                     f"pltB); all three subunits are needed for the holotoxin")
    if _marker_state(markers, "fljB_z66") == "present":
        parts.append("H:z66")
        notes.append("the pBSSB1 fljB allele is present: this is the H:z66 "
                     "antigen, a second flagellar phase carried on a linear "
                     "plasmid rather than on the chromosome, and the reason a "
                     "z66-positive S. Typhi agglutinates as biphasic")
    if _marker_state(markers, "pilL") == "present":
        notes.append("SPI-7 pilL is present (type IVB pilus)")
    call = ", ".join(parts) if parts else "no typhoidal marker"
    return TypingResult(
        scheme="salmonella_typhoidal", call=call,
        lineage="typhoidal marker set" if parts else "-",
        alleles=alleles, loci_found=len(vi) + len(tox),
        loci_total=len(via_b) + len(toxin),
        note=". ".join(notes),
    )


# ---------------------------------------------------------------------- MLST
def _load_profiles(path: Path) -> tuple[list[str], dict[tuple[str, ...], str],
                                        dict[str, str]]:
    loci: list[str] = []
    table: dict[tuple[str, ...], str] = {}
    ebg: dict[str, str] = {}
    if not path.exists():
        return loci, table, ebg
    with open(path, encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        indexes = [i for i, name in enumerate(header) if name in _MLST_LOCI]
        loci = [header[i] for i in indexes]
        cc = next((i for i, name in enumerate(header)
                   if name.lower() in ("clonal_complex", "ebg", "cc")), None)
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(row) <= max(indexes, default=0):
                continue
            st = row[0]
            table[tuple(row[i] for i in indexes)] = st
            if cc is not None and len(row) > cc and row[cc].strip():
                ebg[st] = row[cc].strip()
    return loci, table, ebg


def _call_mlst(assembly: Path, db: Path, profiles: Path, work: Path,
               threads: int, serovar: TypingResult) -> TypingResult:
    hsps = blast("blastn", assembly, db, work / "salmonella.mlst.tsv",
                 threads=threads, evalue=1e-20, task="megablast", word_size=20,
                 perc_identity=90.0, max_target_seqs=20000)
    merged = merge_hsps(hsps)
    best: dict[str, tuple[str, float, float]] = {}
    for hit in merged:
        locus, allele = _split_id(hit.sseqid, 2)[:2]
        if locus not in _MLST_LOCI or not allele:
            continue
        exact = (hit.identity_pct >= MLST_IDENTITY and hit.coverage_pct >= MLST_COVERAGE)
        rank = (exact, hit.coverage_pct, hit.identity_pct,
                -(int(allele) if allele.isdigit() else 1 << 30))
        current = best.get(locus)
        if current is None or rank > current[3]:
            best[locus] = (allele, hit.identity_pct, hit.coverage_pct, rank)

    loci, table, ebg = _load_profiles(profiles)
    order = loci or list(_MLST_LOCI)
    alleles: dict[str, str] = {}
    exact_loci = 0
    for locus in order:
        found = best.get(locus)
        if found is None:
            alleles[locus] = "-"
            continue
        allele, identity, coverage, _rank = found
        if identity >= MLST_IDENTITY and coverage >= MLST_COVERAGE:
            alleles[locus] = allele
            exact_loci += 1
        elif coverage >= MLST_COVERAGE:
            alleles[locus] = f"~{allele}"
        else:
            alleles[locus] = f"{allele}?"
    st = table.get(tuple(alleles[locus] for locus in order), "-")
    group = ebg.get(st, "")

    notes: list[str] = []
    if st == "-":
        inexact = [locus for locus in order if not alleles[locus].isdigit()]
        notes.append(f"no ST assigned: {', '.join(inexact)} "
                     f"{'is' if len(inexact) == 1 else 'are'} not an exact "
                     f"allele match" if inexact
                     else "the allele combination is not in the profile table "
                          "(novel ST)")
    call = f"ST{st}" if st != "-" else "not determined"
    lineage = f"eBG{group}" if group else "-"
    if group:
        notes.append(f"eBurst group {group}. An eBG is a group of STs that share "
                     f"at least six of seven alleles with another member; it is "
                     f"the unit Achtman et al. proposed in place of the serovar, "
                     f"and it is usually a tighter lineage than the serovar name")
    notes.append(f"the serovar predicted from the antigen genes is "
                 f"'{serovar.call}'. Serovar and lineage are reported separately "
                 f"on purpose: most serovars are polyphyletic and some eBGs hold "
                 f"several serovars, so a disagreement between the two is a "
                 f"finding to look up in EnteroBase or PubMLST, not an error to "
                 f"average away")
    formula = serovar.alleles.get("antigenic_formula", "")
    monophasic_i = ":i:-" in formula and "4" in formula.split(":")[0]
    if st == "34" and monophasic_i:
        notes.append("ST34 with a 4,[5],12:i:- formula is the widespread "
                     "epidemic monophasic Typhimurium clone (ASSuT/DT193-DT120 "
                     "in the older phage typing), the single commonest "
                     "multidrug-resistant Salmonella in European pigs and "
                     "people")
    elif st == "19" and monophasic_i:
        notes.append("a 4,[5],12:i:- formula on the ST19 Typhimurium background "
                     "is a monophasic isolate that is not the ST34 epidemic "
                     "clone; monophasic Typhimurium has arisen independently "
                     "more than once")
    elif st == "313":
        notes.append("ST313 is the invasive non-typhoidal Salmonella lineage of "
                     "sub-Saharan Africa: a Typhimurium lineage that behaves "
                     "like a bloodstream pathogen rather than a "
                     "gastroenteritis one, and is commonly multidrug resistant")
    alleles["profile"] = " ".join(f"{locus}({alleles[locus]})" for locus in order)
    alleles["scheme"] = ("PubMLST Salmonella scheme 2, the 7-locus Achtman "
                         "scheme (aroC dnaN hemD hisD purE sucA thrA)")
    return TypingResult(scheme="salmonella_mlst", call=call, lineage=lineage,
                        alleles=alleles, loci_found=exact_loci,
                        loci_total=len(order), note="; ".join(notes))
