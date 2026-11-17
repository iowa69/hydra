"""*Escherichia coli* / *Shigella* typing: phylogroup, fimH, Stx subtype, pathotype.

What this module adds to what Hydra already does
------------------------------------------------
Hydra already serotypes *E. coli* (EcOH, ``serotype``), screens the Kleborate
pathovar marker set (``pathovar_markers``) and runs the Achtman and Pasteur
7-locus MLST schemes. None of that is repeated here. This module answers the
questions those calls leave open:

``phylogroup``
    Which of the eight *E. coli* phylogroups (A, B1, B2, C, D, E, F, G) or
    *Escherichia* cryptic clades the isolate belongs to, by in-silico equivalent
    of the Clermont quadruplex PCR. Reported with the marker pattern.
``fimH``
    The FimTyper fimH allele, which with the ST names the clinically used
    subclones (ST131-*fimH*30 = "H30").
``ST131_subclone``
    The ST131 clade implied by that allele, stated as conditional on the ST.
``stx``
    stx1 and stx2 *subtypes* by the NCBI StxTyper algorithm, because the subtype
    and not the mere presence of stx drives the HUS risk assessment.
``pathotype``
    EPEC (typical/atypical) / EHEC / STEC / ETEC / EAEC / EIEC / DAEC / ExPEC,
    from the defining marker combinations, with the rule applied spelled out.
``intimin``
    The eae (intimin) subtype in the Greek-letter nomenclature, and tir.
``ExPEC_VF``
    The extraintestinal and urinary virulence-factor repertoire, scored by the
    ExPEC definition of Johnson & Russo and the four-gene UPEC rule of Spurbeck.
``shigella_eiec``
    *Shigella* cluster vs enteroinvasive *E. coli*, by the ShigEiFinder
    cluster-specific gene sets. Hydra's species call cannot make this
    distinction - *Shigella* is *E. coli* by any whole-genome measure - so a
    genome called "*Escherichia coli*" upstream may be reported here as
    *S. flexneri*, and that is the call to act on.

Everything runs from a draft assembly with BLAST+ and pure Python. Where a draft
assembly cannot settle a question (the ST131 C1/C2 split, a stx operon broken
across a contig boundary) the module says so and names the candidates rather
than guessing.
"""

from __future__ import annotations

import json
import re
import shutil
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from ...engines.blast import Hsp, blast, deduplicate, merge_hsps
from ...records import TypingResult
from ...seqio import read_fasta
from ...utils import LOG, HydraError, run as run_command
from .base import OrganismDb, Reference

KEY = "ecoli"
TITLE = "Escherichia coli"

#: Species names (as Hydra's SpeciesCall reports them) that trigger the module.
#: *Shigella* is included deliberately: every *Shigella* is an *E. coli* by
#: genome content, every scheme here applies to it, and the ``shigella_eiec``
#: call exists precisely to correct the species label.
SPECIES: tuple[str, ...] = (
    "Escherichia coli",
    "Escherichia fergusonii",
    "Escherichia albertii",
    "Escherichia marmotae",
    "Escherichia ruysiae",
    "Shigella flexneri",
    "Shigella sonnei",
    "Shigella dysenteriae",
    "Shigella boydii",
)

#: PubMLST scheme names that also trigger it (Achtman and Pasteur, as Hydra
#: names them in :data:`hydra_amr.db.fetch.PUBMLST_ALIASES`).
MLST_SCHEMES: tuple[str, ...] = ("ecoli_achtman_4", "ecoli")

USER_AGENT = "hydra-amr (+https://github.com/iowa69/hydra)"
TIMEOUT = 180

# --------------------------------------------------------------------- sources
FIMTYPER_URL = ("https://bitbucket.org/genomicepidemiology/fimtyper_db/raw/"
                "1dd6cdd40bdb967cdbcb049e6e835c1df840f4af/fimH.fsa")
STXTYPER_URL = "https://raw.githubusercontent.com/ncbi/stxtyper/v1.0.45/stx.prot"
VFDB_URL = ("https://bitbucket.org/genomicepidemiology/virulencefinder_db/raw/"
            "master/virulence_ecoli.fsa")
VFDB_NOTES_URL = ("https://bitbucket.org/genomicepidemiology/virulencefinder_db/raw/"
                  "master/notes.txt")
SHIGEI_GENES_URL = ("https://raw.githubusercontent.com/LanLab/ShigEiFinder/"
                    "v1.3.5/shigeifinder/resources/genes.fasta")
SHIGEI_JSON_URL = ("https://raw.githubusercontent.com/LanLab/ShigEiFinder/"
                   "v1.3.5/shigeifinder/genes.json")

DATABASES: tuple[OrganismDb, ...] = (
    OrganismDb(
        name="fimtyper",
        source=FIMTYPER_URL,
        version="commit 1dd6cdd (2023-11-15), 3312 alleles",
        licence=('No licence file in the repository. The FimTyper paper states '
                 '"FimTyper ... is publicly available at '
                 'https://cge.cbs.dtu.dk/services/FimTyper/"; the database '
                 'repository itself states no terms.'),
        redistribution="unstated",
        reference="Roer2017",
        notes="489 bp fimH typing fragment; allele names FimH1..FimH3363.",
    ),
    OrganismDb(
        name="stxtyper",
        source=STXTYPER_URL,
        version="StxTyper v1.0.45 (stx.prot, 160 reference proteins)",
        licence=('PUBLIC DOMAIN NOTICE, NCBI: "This software/database is a '
                 '\\"United States Government Work\\" under the terms of the '
                 'United States Copyright Act ... thus cannot be copyrighted. '
                 'This software/database is freely available to the public for '
                 'use."'),
        redistribution="permitted",
        reference="Prasad2026",
        notes="StxA and StxB reference proteins for stx1a/c/d/e and stx2a-stx2o.",
    ),
    OrganismDb(
        name="virulencefinder",
        source=VFDB_URL,
        version="VirulenceFinder database 2.1.0 (virulence_ecoli.fsa)",
        licence=('No licence file in the repository. The CGE services are '
                 'published for public use; the database repository states no '
                 'terms.'),
        redistribution="unstated",
        reference="Joensen2014",
        notes=("Supplies the eae/intimin subtype alleles (Greek-letter "
               "nomenclature), the DEC pathotype markers and the ExPEC/UPEC "
               "virulence-factor panel. 7.7 MB."),
    ),
    OrganismDb(
        name="shigeifinder",
        source=SHIGEI_GENES_URL,
        version="ShigEiFinder v1.3.5 (genes.fasta + genes.json)",
        licence=("GNU General Public License v3.0 (LICENSE in the ShigEiFinder "
                 "repository). Redistribution is permitted under the GPL-3.0 "
                 "terms, which are copyleft: the integrator must decide whether "
                 "shipping it inside the Windows zip is acceptable, or whether "
                 "it is downloaded on first use."),
        redistribution="permitted",
        reference="Zhang2021",
        notes=("Cluster-specific genes and O/H antigen genes; only the "
               "cluster-specific gene sets and ipaH are used here. 539 kB."),
    ),
)

# ------------------------------------------------------------------ phylogroup
#: IUPAC ambiguity codes, expanded to a character class for the primer regexes.
IUPAC = {
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

#: Clermont primer pairs and expected amplicon sizes, as published in Clermont
#: *et al.* 2013 (quadruplex, group C and group E reactions) and Clermont *et
#: al.* 2019 (the *ybgD* reaction for phylogroup G), transcribed from the
#: EzClermont implementation so the in-silico result matches the reference tool
#: primer for primer. ``(forward, reverse, amplicon_bp)``.
PRIMERS: dict[str, tuple[str, str, int]] = {
    # quadruplex
    "arpA": ("AAYRCYATTCGCCAGCTTGC", "TCTCCMCATACYGYACGCTA", 400),
    "chuA": ("ATGGTACYGGRCGAACCAAC", "TRCCRCCAGTRCCAAAGACA", 288),
    "yjaA": ("YAAACKTGAAGTGTCAGGAG", "ARTRCGTTCCTCAACCTGTG", 211),
    "TspE4.C2": ("CACKATTYGTAAGRYCATCC", "AGTTTATCGCTGCGGGTCGC", 152),
    # group-specific follow-ups
    "arpA_gpE": ("RATKCMATYTTGTCRAAATATGCC", "GAAARKRAAAADAMYYYYCAAGAG", 301),
    "trpA_gpC": ("AGTTYTAYGCCSVRWGCGAG", "TCWGYDCYVGTYACGCCC", 219),
    "ybgD_gpG": ("GTTGACTAARCGYAGGTCGA", "KATGYDGCYGATKAAGGATC", 177),
    # internal control for the group C / group E reactions
    "trpBA_control": ("CGGSGATAAAGAYATYTTCAC", "GCAACGYGSCBWKRCGGAAG", 489),
}

#: The quadruplex decision table of Clermont *et al.* 2013, extended with the
#: phylogroup G and cryptic-clade rows EzClermont adds. Key is
#: ``(arpA, chuA, yjaA, TspE4.C2)``.
QUADRUPLEX: dict[tuple[int, int, int, int], str] = {
    (1, 0, 0, 0): "A",
    (1, 0, 0, 1): "B1",
    (0, 1, 1, 0): "B2",
    (0, 1, 1, 1): "B2",
    (1, 0, 1, 0): "A/C",
    (1, 1, 0, 0): "D/E",
    (1, 1, 0, 1): "D/E",
    (0, 1, 0, 0): "F/G",
    (0, 1, 0, 1): "B2/G",
    (1, 1, 1, 0): "E/cryptic",
    (0, 0, 1, 0): "cryptic",
    (0, 0, 0, 0): "U/cryptic",
    (0, 0, 0, 1): "U",
    (0, 0, 1, 1): "U",
    (1, 0, 1, 1): "U",
    (1, 1, 1, 1): "U",
}

#: What a phylogroup means, in one clause, for the report.
PHYLOGROUP_MEANING = {
    "A": "commensal/generalist; most ETEC and laboratory K-12 lineages",
    "B1": "commensal and animal-associated; most non-O157 STEC",
    "B2": "extraintestinal pathogenic lineages (ST131, ST73, ST95); "
          "the richest virulence-gene content",
    "C": "sister group of A; often ETEC and avian lineages",
    "D": "extraintestinal pathogenic lineages (ST69, ST405) and EAEC",
    "E": "the O157:H7 / O55:H7 EHEC lineages and their ancestors",
    "F": "sister group of B2; extraintestinal, includes ST648 and ST354",
    "G": "recently defined; high virulence and resistance potential, "
         "includes ST117 avian/extraintestinal lineages",
    "cryptic clade": "an Escherichia cryptic clade (clade I-V), not E. coli "
                     "sensu stricto; environmental, rarely clinical",
}

# -------------------------------------------------------------------------- stx
#: Minimum combined (A+B) operon identity for a subtype call, per StxTyper
#: v1.0.45 (``stxClass2identity`` in ``stxtyper.cpp``). Keyed by the stx class,
#: where "2" is the generalised stx2a/stx2c/stx2d class.
STX_IDENTITY: dict[str, float] = {
    "1a": 0.983, "1c": 0.983, "1d": 0.983, "1e": 0.983,
    "2": 0.980, "2b": 0.980, "2e": 0.980, "2f": 0.980, "2g": 0.980,
    "2h": 0.980, "2i": 0.980, "2j": 0.980, "2k": 0.985, "2l": 0.985,
    "2m": 0.980, "2n": 0.980, "2o": 0.980,
}
#: Maximum intergenic distance between stxA and stxB of one operon, in bp.
#: StxTyper uses 36 ("max. intergenic region in the reference set + 2").
STX_INTERGENIC_MAX = 36
#: Only HSPs at or above this identity are considered at all (StxTyper PD-5192).
STX_MIN_IDENTITY = 80.0

#: The three residues that separate stx2a, stx2c and stx2d, which are otherwise
#: one generalised class. Positions are 1-based in the *reference* protein:
#: StxA2 residues 313 and 319, StxB2 residue 35. Rules are tried in this order.
STX2ACD_RULES: tuple[tuple[str, tuple[str, ...], tuple[str, ...], tuple[str, ...]], ...] = (
    ("2a", ("F", "S"), ("K", "E"), ("D",)),
    ("2c", ("F",), ("K", "E"), ("N",)),
    ("2d", ("S",), ("E",), ("N",)),
)
STX_A_POS = (313, 319)
STX_B_POS = 35

#: What each subtype means for the clinical risk assessment. The subtype, not
#: the presence of stx, is what the HUS literature stratifies on.
STX_MEANING = {
    "stx1a": "classical Stx1; HUS risk much lower than stx2a",
    "stx1c": "mainly sheep-associated, mild or asymptomatic infection",
    "stx1d": "rare, low virulence",
    "stx1e": "rare, not associated with human disease",
    "stx2a": "the subtype most strongly associated with HUS and with "
             "bloody diarrhoea; highest-risk call in this scheme",
    "stx2b": "usually mild diarrhoea, HUS uncommon",
    "stx2c": "associated with HUS, classically with stx2a in O157:H7; "
             "on its own a lower risk than stx2a",
    "stx2d": "activatable Stx2d; HUS risk comparable to stx2a in eae-negative "
             "strains",
    "stx2e": "porcine oedema disease; human cases rare and mild",
    "stx2f": "pigeon-associated; mild human disease",
    "stx2g": "bovine; human disease rare",
}

# ------------------------------------------------------------------- virulence
#: Identity and reference-coverage thresholds for calling a virulence gene
#: present. These are the VirulenceFinder web defaults (90 % identity, 60 %
#: coverage); the eae subtype call uses a stricter rule, below.
VF_MIN_IDENTITY = 90.0
VF_MIN_COVERAGE = 60.0
#: An eae allele is only named as a subtype at this identity or better: the 65
#: intimin alleles are 85-98 % identical to one another, so a loose threshold
#: names the wrong Greek letter.
EAE_SUBTYPE_IDENTITY = 97.0

#: Marker groups the pathotype rules are written in terms of. A group is present
#: when any of its genes is. Names are the VirulenceFinder gene names, matched
#: as a prefix, so ``eltIAB`` covers the 33 numbered LT-I alleles.
MARKER_GROUPS: dict[str, tuple[str, ...]] = {
    "eae": ("eae",),
    "tir": ("tir",),
    "bfpA": ("bfpA",),
    "perA": ("perA",),
    "LT": ("eltIAB", "eltI-", "eltIIAB"),
    "ST": ("estah-STa", "estap-STa", "estb-STb"),
    # ETEC colonisation factors. Named individually: a bare "cs" prefix also
    # matched csgA, the curli subunit, which every E. coli carries and which is
    # not a colonisation factor antigen.
    "CF": ("cfaA", "cfaB", "cfaC", "cfaD", "cfaE", "lngA", "cmaH", "cooA",
           "cstH", "cotA", "csbA", "cseA", "cssA", "cswA"),
    "aggR": ("aggR",),
    "AAF": ("aggA", "aafA", "agg3A", "agg4A", "agg5A"),
    "aat": ("aatA",),
    "aai": ("aaiC",),
    "astA": ("astA",),
    "aap": ("aap",),
    "ipaH": ("ipaH",),
    "invasion_plasmid": ("ipaD", "virF", "sepA", "senB"),
    "afa_dra": ("afaA", "afaB", "afaC", "afaD", "afaE"),
    "pap": ("papA", "papC"),
    "sfa_foc": ("sfaD", "sfaE", "sfaS", "focC", "focG", "focI", "focCsfaE"),
    "kpsMII": ("kpsMII",),
    "iutA": ("iutA", "iucC"),
    "ehec_accessory": ("ehxA", "katP", "espP", "toxB", "nleA", "nleB", "nleC"),
}

#: The ExPEC molecular definition: an isolate is ExPEC when at least two of
#: these five marker groups are present (Johnson & Russo).
EXPEC_GROUPS = ("pap", "sfa_foc", "afa_dra", "kpsMII", "iutA")
#: The four-gene rule of Spurbeck *et al.* 2012: an isolate carrying all four
#: efficiently colonises the urinary tract.
UPEC_FOUR = ("vat", "fyuA", "chuA", "yfcV")

#: The extraintestinal/urinary virulence-factor panel, grouped as the
#: urinary-sepsis literature groups it. Reported as a repertoire, not a score,
#: because no single published score is standard.
VF_PANEL: dict[str, tuple[str, ...]] = {
    "adhesins": ("papA", "papC", "sfaD", "sfaE", "sfaS", "focG", "focCsfaE",
                 "afaB", "afaC", "afaE", "iha", "yfcV"),
    "toxins": ("hlyA", "hlyE", "hlyF", "cnf1", "cnf2", "sat", "usp", "vat",
               "astA", "cdt-IB"),
    "iron_uptake": ("chuA", "fyuA", "iroN", "iutA", "iucC", "ireA", "sitA"),
    "protectins_capsule": ("kpsMII", "kpsMIII", "kpsE", "traT", "ompT", "iss"),
    "other": ("ibeA", "clbB", "espC", "espI"),
}

# -------------------------------------------------------------- ST131 subclone
#: fimH allele -> ST131 clade. Only the three alleles the ST131 clade papers
#: actually name are listed: fimH41 for clade A, fimH22 for clade B and fimH30
#: for clade C (Price *et al.* 2013, Johnson *et al.* 2014, Matsumura *et al.*
#: 2017). Other alleles occur inside ST131 but no published table maps them to a
#: clade, and guessing would put a clade letter on an isolate that has none.
#:
#: These alleles are common outside ST131 too - fimH41 in particular - so the
#: label is strictly conditional on the ST, and the call says so in words rather
#: than relying on the reader to check.
ST131_FIMH_CLADE: dict[str, tuple[str, str]] = {
    "fimH41": ("A", "clade A; serotype O16:H5, usually fluoroquinolone-susceptible"),
    "fimH22": ("B", "clade B; serotype O25b:H4, usually fluoroquinolone-susceptible"),
    "fimH30": ("C", "clade C (the H30 subclone); serotype O25b:H4, "
                    "fluoroquinolone-resistant, the pandemic ST131 lineage"),
}

# ----------------------------------------------------------- Shigella / EIEC
#: What each ShigEiFinder cluster is, and whether it is Shigella or EIEC.
#: Taken from Zhang *et al.* 2021 table 1 and the serotype sets in the
#: distributed ``genes.json``.
SHIGEI_CLUSTERS: dict[str, tuple[str, str]] = {
    "C1": ("Shigella", "cluster 1: S. boydii and S. dysenteriae serotypes "
                       "(and S. flexneri 6 relatives)"),
    "C2": ("Shigella", "cluster 2: S. boydii and S. dysenteriae serotypes"),
    "C3": ("Shigella", "cluster 3: S. flexneri serotypes 1-5, X, Y and "
                       "S. boydii provisional serotypes"),
    "C4": ("EIEC", "EIEC cluster 4 (O28ac, O136, O121 and relatives)"),
    "C5": ("EIEC", "EIEC cluster 5 (O121:H30, O124:H30, O164:H30 and relatives)"),
    "C6": ("EIEC", "EIEC cluster 6 (O143:H26, O167:H26, O112ac:H26)"),
    "C7": ("EIEC", "EIEC cluster 7 (O144:H25)"),
    "C8": ("EIEC", "EIEC cluster 8 (O96:H19)"),
    "C9": ("EIEC", "EIEC cluster 9 (O8:H19)"),
    "C10": ("EIEC", "EIEC cluster 10 (O135:H30)"),
    "CSS": ("Shigella", "S. sonnei"),
    "CSD1": ("Shigella", "S. dysenteriae serotype 1 (Shiga's bacillus)"),
    "CSD8": ("Shigella", "S. dysenteriae serotype 8"),
    "CSD10": ("Shigella", "S. dysenteriae serotype 10"),
    "CSB12": ("Shigella", "S. boydii serotype 12"),
    "CSB13": ("Shigella", "S. boydii serotype 13"),
    "CSB13-atypical": ("Shigella", "S. boydii serotype 13, atypical lineage"),
}
#: Cluster-specific gene presence threshold, and the three genes ShigEiFinder
#: additionally requires a minimum aligned length for.
SHIGEI_MIN_COVERAGE = 50.0
SHIGEI_MIN_ALIGNED = {"C1_gene_2": 251, "C1_gene_4": 254, "C2_gene_1": 393}
#: ipaH is called present on far less: it is multi-copy and often partial.
SHIGEI_IPAH_COVERAGE = 10.0

# ------------------------------------------------------------------ references
REFERENCES: tuple[Reference, ...] = (
    Reference(
        key="Clermont2013",
        authors="Clermont O, Christenson JK, Denamur E, Gordon DM.",
        year=2013,
        title=("The Clermont Escherichia coli phylo-typing method revisited: "
               "improvement of specificity and detection of new phylo-groups"),
        journal="Environmental Microbiology Reports 5:58-65",
        doi="10.1111/1758-2229.12019", pmid="23757131",
        supports=("the quadruplex primers, the group C and group E reactions "
                  "and the decision table behind the `phylogroup` call"),
    ),
    Reference(
        key="Clermont2000",
        authors="Clermont O, Bonacorsi S, Bingen E.",
        year=2000,
        title="Rapid and simple determination of the Escherichia coli phylogenetic group",
        journal="Applied and Environmental Microbiology 66:4555-4558",
        doi="10.1128/AEM.66.10.4555-4558.2000", pmid="11010916",
        supports="the original triplex scheme the 2013 method revises",
    ),
    Reference(
        key="Clermont2019",
        authors=("Clermont O, Dixit OVA, Vangchhia B, Condamine B, Dion S, "
                 "Bridier-Nahmias A, Denamur E, Gordon D."),
        year=2019,
        title=("Characterization and rapid identification of phylogroup G in "
               "Escherichia coli, a lineage with high virulence and antibiotic "
               "resistance potential"),
        journal="Environmental Microbiology 21:3107-3117",
        doi="10.1111/1462-2920.14713", pmid="31188527",
        supports="phylogroup G and the ybgD reaction that identifies it",
    ),
    Reference(
        key="Waters2020",
        authors="Waters NR, Abram F, Brennan F, Holmes A, Pritchard L.",
        year=2020,
        title=("Easy phylotyping of Escherichia coli via the EzClermont web app "
               "and command-line tool"),
        journal="Access Microbiology 2:acmi000143",
        doi="10.1099/acmi.0.000143", pmid="33195978",
        supports=("the in-silico PCR implementation this module reproduces, and "
                  "the reference tool the phylogroup call was validated against"),
    ),
    Reference(
        key="Roer2017",
        authors=("Roer L, Tchesnokova V, Allesoe R, Muradova M, Chattopadhyay S, "
                 "Ahrenfeldt J, Thomsen MCF, Lund O, Hansen F, Hammerum AM, "
                 "Sokurenko E, Hasman H."),
        year=2017,
        title=("Development of a web tool for Escherichia coli subtyping based "
               "on fimH alleles"),
        journal="Journal of Clinical Microbiology 55:2538-2543",
        doi="10.1128/JCM.00737-17", pmid="28592545",
        supports="the fimH allele nomenclature and database used by `fimH`",
    ),
    Reference(
        key="Scheutz2012",
        authors=("Scheutz F, Teel LD, Beutin L, Pierard D, Buvens G, Karch H, "
                 "Mellmann A, Caprioli A, Tozzoli R, Morabito S, Strockbine NA, "
                 "Melton-Celsa AR, Sanchez M, Persson S, O'Brien AD."),
        year=2012,
        title=("Multicenter evaluation of a sequence-based protocol for "
               "subtyping Shiga toxins and standardizing Stx nomenclature"),
        journal="Journal of Clinical Microbiology 50:2951-2963",
        doi="10.1128/JCM.00860-12", pmid="22760050",
        supports="the stx subtype nomenclature (stx1a/c/d, stx2a-g) `stx` reports",
    ),
    Reference(
        key="Prasad2026",
        authors=("Prasad AB, Abromaitis S, Brover V, Feldgarden M, Joensen KG, "
                 "Kapsak CJ, Lindsey RL, Li L, Michelacci V, Schjorring S, "
                 "Scheutz F, Klimke W."),
        year=2026,
        title=("Algorithm and software to type Stx operons accurately from "
               "assembled genomic sequence"),
        journal="Microorganisms 14:1607",
        doi="10.3390/microorganisms14081607", pmid="42654954",
        supports=("the StxTyper algorithm reimplemented in `stx`: the operon "
                  "identity cut-offs, the 36 bp intergenic limit and the "
                  "stx2a/2c/2d discriminating residues; and the reference "
                  "database and tool it was validated against"),
    ),
    Reference(
        key="Lindsey2023",
        authors=("Lindsey RL, Prasad A, Feldgarden M, Gonzalez-Escalona N, "
                 "Kapsak C, Klimke W, Melton-Celsa A, Smith P, Souvorov A, "
                 "Truong J, Scheutz F."),
        year=2023,
        title=("Identification and characterization of ten Escherichia coli "
               "strains encoding novel Shiga toxin 2 subtypes, Stx2n as well "
               "as Stx2j, Stx2m, and Stx2o, in the United States"),
        journal="Microorganisms 11:2561",
        doi="10.3390/microorganisms11102561", pmid="37894219",
        supports="the stx2j-stx2o subtypes in the current nomenclature",
    ),
    Reference(
        key="Nataro1998",
        authors="Nataro JP, Kaper JB.",
        year=1998,
        title="Diarrheagenic Escherichia coli",
        journal="Clinical Microbiology Reviews 11:142-201",
        doi="10.1128/CMR.11.1.142", pmid="9457432",
        supports="the pathotype definitions (EPEC, EHEC, ETEC, EAEC, EIEC, DAEC)",
    ),
    Reference(
        key="Croxen2013",
        authors=("Croxen MA, Law RJ, Scholz R, Keeney KM, Wlodarska M, "
                 "Finlay BB."),
        year=2013,
        title="Recent advances in understanding enteric pathogenic Escherichia coli",
        journal="Clinical Microbiology Reviews 26:822-880",
        doi="10.1128/CMR.00022-13", pmid="24092857",
        supports=("the current marker combinations for each pathotype, "
                  "including the typical/atypical EPEC split"),
    ),
    Reference(
        key="Boisen2020",
        authors=("Boisen N, Osterlund MT, Joensen KG, Santiago AE, Mandomando I, "
                 "Cravioto A, Chattopadhyay S, Gonyar LA, Jensen SG, Nataro JP, "
                 "Scheutz F, Hazen TH, Rasko DA."),
        year=2020,
        title=("Redefining enteroaggregative Escherichia coli (EAEC): genomic "
               "characterization of epidemiological EAEC strains"),
        journal="PLoS Neglected Tropical Diseases 14:e0008613",
        doi="10.1371/journal.pntd.0008613", pmid="32898134",
        supports="the aggR-based EAEC definition used by `pathotype`",
    ),
    Reference(
        key="Joensen2014",
        authors=("Joensen KG, Scheutz F, Lund O, Hasman H, Kaas RS, Nielsen EM, "
                 "Aarestrup FM."),
        year=2014,
        title=("Real-time whole-genome sequencing for routine typing, "
               "surveillance, and outbreak detection of verotoxigenic "
               "Escherichia coli"),
        journal="Journal of Clinical Microbiology 52:1501-1510",
        doi="10.1128/JCM.03617-13", pmid="24574290",
        supports=("the VirulenceFinder database and its 90 % identity / 60 % "
                  "coverage presence rule, used for the pathotype markers, the "
                  "eae subtypes and the ExPEC/UPEC panel"),
    ),
    Reference(
        key="Oswald2000",
        authors=("Oswald E, Schmidt H, Morabito S, Karch H, Marches O, "
                 "Caprioli A."),
        year=2000,
        title=("Typing of intimin genes in human and animal enterohemorrhagic "
               "and enteropathogenic Escherichia coli: characterization of a "
               "new intimin variant"),
        journal="Infection and Immunity 68:64-71",
        doi="10.1128/IAI.68.1.64-71.2000", pmid="10603369",
        supports="the Greek-letter intimin (eae) subtype nomenclature",
    ),
    Reference(
        key="Zhang2002",
        authors=("Zhang WL, Kohler B, Oswald E, Beutin L, Karch H, Morabito S, "
                 "Caprioli A, Suerbaum S, Schmidt H."),
        year=2002,
        title="Genetic diversity of intimin genes of attaching and effacing Escherichia coli strains",
        journal="Journal of Clinical Microbiology 40:4486-4492",
        doi="10.1128/JCM.40.12.4486-4492.2002", pmid="12454140",
        supports="the extension of the intimin nomenclature beyond the first variants",
    ),
    Reference(
        key="Russo2000",
        authors="Russo TA, Johnson JR.",
        year=2000,
        title=("Proposal for a new inclusive designation for extraintestinal "
               "pathogenic isolates of Escherichia coli: ExPEC"),
        journal="Journal of Infectious Diseases 181:1753-1754",
        doi="10.1086/315418", pmid="10823778",
        supports="the ExPEC designation itself",
    ),
    Reference(
        key="Johnson2003",
        authors=("Johnson JR, Kuskowski MA, Owens K, Gajewski A, Winokur PL."),
        year=2003,
        title=("Phylogenetic origin and virulence genotype in relation to "
               "resistance to fluoroquinolones and/or extended-spectrum "
               "cephalosporins and cephamycins among Escherichia coli isolates "
               "from animals and humans"),
        journal="Journal of Infectious Diseases 188:759-768",
        doi="10.1086/377455", pmid="12934193",
        supports=("the five-marker ExPEC rule (two or more of papAH/papC, "
                  "sfa/focDE, afa/draBC, kpsM II, iutA) applied by `ExPEC_VF`"),
    ),
    Reference(
        key="Spurbeck2012",
        authors=("Spurbeck RR, Dinh PC Jr, Walk ST, Stapleton AE, Hooton TM, "
                 "Nolan LK, Kim KS, Johnson JR, Mobley HLT."),
        year=2012,
        title=("Escherichia coli isolates that carry vat, fyuA, chuA, and yfcV "
               "efficiently colonize the urinary tract"),
        journal="Infection and Immunity 80:4115-4122",
        doi="10.1128/IAI.00752-12", pmid="22966046",
        supports="the four-gene UPEC rule reported by `ExPEC_VF`",
    ),
    Reference(
        key="Price2013",
        authors=("Price LB, Johnson JR, Aziz M, Clabots C, Johnston B, "
                 "Tchesnokova V, Nordstrom L, Billig M, Chattopadhyay S, "
                 "Stegger M, Andersen PS, Pearson T, Riddell K, Rogers P, "
                 "Scholes D, Kahl B, Keim P, Sokurenko EV."),
        year=2013,
        title=("The epidemic of extended-spectrum-beta-lactamase-producing "
               "Escherichia coli ST131 is driven by a single highly pathogenic "
               "subclone, H30-Rx"),
        journal="mBio 4:e00377-13",
        doi="10.1128/mBio.00377-13", pmid="24345742",
        supports=("the H30 and H30Rx subclone definitions, and the fact that "
                  "the H30Rx split rests on core-genome SNPs rather than on "
                  "fimH, which is why `ST131_subclone` stops at clade C"),
    ),
    Reference(
        key="Matsumura2017",
        authors=("Matsumura Y, Pitout JDD, Peirano G, DeVinney R, Noguchi T, "
                 "Yamamoto M, Gomi R, Matsuda T, Nakano S, Nagao M, Tanaka M, "
                 "Ichiyama S."),
        year=2017,
        title="Rapid identification of different Escherichia coli sequence type 131 clades",
        journal="Antimicrobial Agents and Chemotherapy 61:e00179-17",
        doi="10.1128/AAC.00179-17", pmid="28584160",
        supports=("the ST131 clade A/B/C scheme and its association with the "
                  "fimH alleles fimH41, fimH22 and fimH30"),
    ),
    Reference(
        key="Johnson2014",
        authors=("Johnson JR, Johnston B, Clabots C, Kuskowski MA, Pendyala S, "
                 "DebRoy C, Nowicki B, Rice J, Sokurenko EV."),
        year=2014,
        title=("Rapid and specific detection, molecular epidemiology, and "
               "experimental virulence of the O16 subgroup within Escherichia "
               "coli sequence type 131"),
        journal="Journal of Clinical Microbiology 52:1358-1365",
        doi="10.1128/JCM.03502-13", pmid="24501035",
        supports="the O16:H5 / fimH41 clade A subgroup of ST131",
    ),
    Reference(
        key="Zhang2021",
        authors="Zhang X, Payne M, Nguyen T, Kaur S, Lan R.",
        year=2021,
        title=("Cluster-specific gene markers enhance Shigella and "
               "enteroinvasive Escherichia coli in silico serotyping"),
        journal="Microbial Genomics 7:000704",
        doi="10.1099/mgen.0.000704", pmid="34889728",
        supports=("the ShigEiFinder cluster-specific gene sets, the all-genes "
                  "presence rule and the Shigella-vs-EIEC assignment in "
                  "`shigella_eiec`"),
    ),
    Reference(
        key="Yassine2023",
        authors=("Yassine I, Lefevre S, Hansen EE, Ruckly C, Carle I, "
                 "Lejay-Collin M, Fabre L, Rapoud D, Coignard B, "
                 "Weill FX, Bonacorsi S, Mariani-Kurkdjian P."),
        year=2023,
        title=("ShigaPass: an in silico tool predicting Shigella serotypes "
               "from whole-genome sequencing assemblies"),
        journal="Microbial Genomics 9:000961",
        doi="10.1099/mgen.0.000961", pmid="36951906",
        supports=("the alternative published in-silico Shigella solution, "
                  "named in the docs as the method not chosen and why"),
    ),
)

# --------------------------------------------------------------------- fetching
def _download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT) as response, \
                open(dest, "wb") as handle:
            shutil.copyfileobj(response, handle)
    except (urllib.error.URLError, OSError) as exc:
        raise HydraError(f"could not download {url}: {exc}") from exc
    if not dest.exists() or dest.stat().st_size == 0:
        raise HydraError(f"{url} returned an empty file")


def _makeblastdb(fasta: Path, dbtype: str) -> None:
    suffix = ".nin" if dbtype == "nucl" else ".pin"
    if (fasta.parent / (fasta.name + suffix)).exists():
        return
    run_command(["makeblastdb", "-in", str(fasta), "-dbtype", dbtype,
                 "-title", fasta.stem, "-out", str(fasta)])


def fetch(db_dir: Path) -> None:
    """Download and index this module's data into ``db_dir/ecoli/``.

    Idempotent: a file already present with a BLAST index beside it is left
    alone, so re-running costs nothing.
    """
    root = Path(db_dir) / KEY
    jobs = [
        ("fimtyper/fimH.fsa", FIMTYPER_URL, "nucl"),
        ("stxtyper/stx.prot", STXTYPER_URL, "prot"),
        ("virulencefinder/virulence_ecoli.fsa", VFDB_URL, "nucl"),
        ("virulencefinder/notes.txt", VFDB_NOTES_URL, ""),
        ("shigeifinder/genes.fasta", SHIGEI_GENES_URL, "nucl"),
        ("shigeifinder/genes.json", SHIGEI_JSON_URL, ""),
    ]
    for relative, url, dbtype in jobs:
        dest = root / relative
        if not dest.exists() or dest.stat().st_size == 0:
            LOG.info("ecoli typing: fetching %s", relative)
            _download(url, dest)
        if dbtype:
            _makeblastdb(dest, dbtype)
    provenance = root / "provenance.tsv"
    lines = ["name\tsource\tversion\tlicence\tredistribution\treference"]
    for entry in DATABASES:
        licence = entry.licence.replace("\t", " ").replace("\n", " ")
        lines.append("\t".join([entry.name, entry.source, entry.version,
                                licence, entry.redistribution, entry.reference]))
    provenance.write_text("\n".join(lines) + "\n", encoding="utf-8")


class _MissingDb(HydraError):
    """One of this module's databases is absent.

    Its own subclass so :func:`run` can turn it into a "not determined" call
    for the one scheme that needed it, instead of losing every other call the
    assembly could still have been given. A packaged build may deliberately
    carry only part of this module's data -- ShigEiFinder's gene sets are
    GPL-3.0 and are downloaded rather than shipped -- and the phylogroup,
    fimH and pathotype calls must survive that.
    """


def _require_db(root: Path, relative: str, what: str) -> Path:
    path = root / relative
    if not path.exists():
        raise _MissingDb(
            f"{KEY} typing: {what} is not installed ({path}); "
            f"run the module's fetch() or 'hydra db download typing-{KEY}'")
    return path


# ------------------------------------------------------------- 1. phylogroup
def _primer_regex(primer: str) -> re.Pattern[str]:
    return re.compile("".join(IUPAC.get(base, base) for base in primer),
                      re.IGNORECASE)


_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                            "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class _PcrResult:
    hit: bool
    evidence: str


def _in_silico_pcr(contigs: list[tuple[str, str]], name: str,
                   allow_partial: bool) -> _PcrResult:
    """One Clermont PCR reaction, run against the assembly as EzClermont does.

    Both primers are searched on both orientations of every contig; a product is
    called when they are found the expected distance apart (within 20 bp). When
    the assembly is fragmented, a single primer whose contig ends before the
    expected product length could have is accepted as a partial product, which
    is what the reference tool does and what keeps a draft assembly from losing
    a marker to a contig break.
    """
    forward, reverse, size = PRIMERS[name]
    fwd = _primer_regex(forward)
    rev = _primer_regex(_revcomp(reverse))
    best = ""
    for header, sequence in contigs:
        contig = header.split()[0]
        for strand, seq in (("+", sequence), ("-", _revcomp(sequence))):
            f_match = fwd.search(seq)
            r_match = rev.search(seq)
            if f_match and r_match:
                product = r_match.end() - f_match.start()
                if abs(product - size) < 20:
                    return _PcrResult(True, f"{product} bp on {contig}({strand})")
                best = best or f"product {product} bp, expected {size} bp"
            elif allow_partial:
                match = f_match or r_match
                if match is not None and len(seq) - match.start() < size:
                    return _PcrResult(
                        True, f"partial, one primer {len(seq) - match.start()} bp "
                              f"from the end of {contig}({strand})")
    return _PcrResult(False, best or "no product")


def _phylogroup(contigs: list[tuple[str, str]]) -> TypingResult:
    allow_partial = len(contigs) > 4
    reactions = {name: _in_silico_pcr(contigs, name, allow_partial)
                 for name in ("trpBA_control", "arpA", "chuA", "yjaA", "TspE4.C2")}
    quad = ("arpA", "chuA", "yjaA", "TspE4.C2")
    pattern = tuple(1 if reactions[name].hit else 0 for name in quad)
    call = QUADRUPLEX.get(pattern, "U")
    notes: list[str] = []

    control_ok = reactions["trpBA_control"].hit
    if call in ("A/C", "D/E", "E/cryptic") and not control_ok:
        # The group C and group E reactions are read against this internal
        # control; without it a negative cannot be distinguished from a failed
        # reaction, so the reference tool refuses to refine and so does this.
        notes.append("the trpBA internal control gave no product, so the "
                     "group C / group E reaction cannot be read; "
                     f"pattern is consistent with {call}")
        final = f"{call} (control failed)"
    else:
        final, extra = _refine_phylogroup(call, contigs, allow_partial, reactions)
        notes.extend(extra)

    alleles = {name: ("+" if reactions[name].hit else "-") for name in quad}
    alleles.update({f"{name}_evidence": reactions[name].evidence for name in quad})
    for name in ("arpA_gpE", "trpA_gpC", "ybgD_gpG"):
        if name in reactions:
            alleles[name] = "+" if reactions[name].hit else "-"
    alleles["trpBA_control"] = "+" if control_ok else "-"
    alleles["quadruplex_pattern"] = "".join(
        f"{n.split('.')[0]}{'+' if reactions[n].hit else '-'} " for n in quad).strip()

    meaning = PHYLOGROUP_MEANING.get(final.split(" ")[0], "")
    if final.startswith("cryptic"):
        meaning = PHYLOGROUP_MEANING["cryptic clade"]
    if final.startswith("U"):
        notes.append("the marker pattern is not in the Clermont table; this "
                     "happens with cryptic clades, with Escherichia species "
                     "other than E. coli, and with mixed or contaminated "
                     "assemblies - check the species call")
    return TypingResult(
        scheme="phylogroup",
        call=final,
        lineage=meaning or "-",
        alleles=alleles,
        loci_found=sum(pattern),
        loci_total=len(quad),
        note="; ".join(notes),
    )


def _refine_phylogroup(call: str, contigs: list[tuple[str, str]],
                       allow_partial: bool,
                       reactions: dict[str, _PcrResult]) -> tuple[str, list[str]]:
    """Resolve the ambiguous quadruplex patterns with the follow-up reactions."""
    notes: list[str] = []
    if call in ("D/E", "E/cryptic"):
        result = _in_silico_pcr(contigs, "arpA_gpE", allow_partial)
        reactions["arpA_gpE"] = result
        if result.hit:
            return "E", notes
        return ("D" if call == "D/E" else "cryptic clade"), notes
    if call == "A/C":
        result = _in_silico_pcr(contigs, "trpA_gpC", allow_partial)
        reactions["trpA_gpC"] = result
        return ("C" if result.hit else "A"), notes
    if call in ("F/G", "B2/G"):
        result = _in_silico_pcr(contigs, "ybgD_gpG", allow_partial)
        reactions["ybgD_gpG"] = result
        if result.hit:
            return "G", notes
        return ("F" if call == "F/G" else "B2"), notes
    if call == "cryptic":
        return "cryptic clade", notes
    if call == "U/cryptic":
        notes.append("all four quadruplex markers negative: either an "
                     "Escherichia cryptic clade or an assembly the primers do "
                     "not match; not resolved")
        return "not determined", notes
    return call, notes


# ----------------------------------------------------------------- 2. fimH
def _fimh(assembly: Path, root: Path, workdir: Path, threads: int) -> TypingResult:
    fasta = _require_db(root, "fimtyper/fimH.fsa", "the FimTyper fimH database")
    hsps = blast("blastn", assembly, fasta, workdir / "ecoli.fimh.tsv",
                 threads=threads, evalue=1e-20, task="megablast",
                 perc_identity=90.0, max_target_seqs=10000)
    all_merged = merge_hsps(hsps)
    merged = [hit for hit in all_merged if hit.coverage_pct >= 60.0]
    if not merged:
        if not all_merged:
            # No fimH-like sequence anywhere, which is a different finding
            # from a fragment too short to type, and is reported as such.
            return TypingResult(
                scheme="fimH", call="fimH not found", lineage="-",
                loci_found=0, loci_total=1,
                note=("no sequence anywhere in the assembly aligns to any of "
                      "the 3312 fimH typing fragments at 90 % identity, so "
                      "this is an absence of the gene in the assembly rather "
                      "than an untypable allele. fimH is present in most but "
                      "not all E. coli, so confirm against the assembly's "
                      "completeness before reporting it as a true deletion"))
        best_partial = max(all_merged, key=lambda h: h.match_score)
        return TypingResult(
            scheme="fimH", call="not determined", lineage="-",
            alleles={"nearest_allele": _fimh_name(best_partial.sseqid),
                     "identity_pct": f"{best_partial.identity_pct:.2f}",
                     "coverage_pct": f"{best_partial.coverage_pct:.1f}",
                     "contig": best_partial.qseqid},
            loci_found=0, loci_total=1,
            note=(f"fimH-like sequence is present but only "
                  f"{best_partial.coverage_pct:.1f} % of the nearest allele is "
                  f"covered, below the 60 % needed to type it. A contig break "
                  f"through fimH or an insertion in the gene does this; no "
                  f"allele is called from a fragment"))
    # Several alleles differ by one base, so rank by the fraction of the
    # reference matched identically and break ties on the longer reference:
    # a short allele that is a substring of the right one otherwise wins.
    best = max(merged, key=lambda h: (h.match_score, h.slen, h.bitscore))
    allele = _fimh_name(best.sseqid)
    exact = best.identity_pct >= 100.0 and best.coverage_pct >= 100.0
    ties = sorted({_fimh_name(h.sseqid) for h in merged
                   if h.identity_pct >= 100.0 and h.coverage_pct >= 100.0
                   and h.slen == best.slen})
    alleles = {
        "allele": allele,
        "identity_pct": f"{best.identity_pct:.2f}",
        "coverage": best.coverage_string,
        "contig": best.qseqid,
        "contig_span": f"{best.query_start}-{best.query_end}({best.strand})",
    }
    note = ""
    if not exact:
        allele_call = f"{allele}*"
        note = (f"no exact match in the FimTyper database; nearest allele "
                f"{allele} at {best.identity_pct:.2f} % identity over "
                f"{best.coverage_pct:.1f} % of the reference. Treat as a novel "
                f"or untyped allele, not as {allele}")
    elif len(ties) > 1:
        allele_call = "/".join(ties)
        note = ("several alleles match the typing fragment exactly and cannot "
                "be separated from an assembly")
    else:
        allele_call = allele
    clade = ST131_FIMH_CLADE.get(allele)
    lineage = ""
    if exact and clade:
        lineage = f"ST131 clade {clade[0]} allele" if clade else ""
    return TypingResult(
        scheme="fimH", call=allele_call, lineage=lineage or "-",
        alleles=alleles, loci_found=1, loci_total=1, note=note)


def _fimh_name(sseqid: str) -> str:
    """``FimH30`` from the database header, spelled as the literature does."""
    name = sseqid.split(":")[0].strip()
    match = re.fullmatch(r"[Ff]im[Hh](\d+)", name)
    return f"fimH{match.group(1)}" if match else name


def _st131_subclone(fimh: TypingResult) -> TypingResult:
    """The ST131 subclone label the fimH allele implies, if the ST is 131.

    The subclone names in clinical use - H30, H30Rx, H22, H41 - are fimH
    allele labels applied inside ST131. This module cannot confirm the ST
    (Hydra's MLST module does that), so the call is reported as conditional,
    and the C1/C2 (H30R/H30Rx) split is reported as unresolvable from an
    assembly by this method, which it is: it rests on core-genome SNPs.
    """
    allele = fimh.call.rstrip("*").split("/")[0]
    if not re.fullmatch(r"fimH\d+", allele):
        return TypingResult(
            scheme="ST131_subclone", call="not determined", lineage="-",
            alleles={"fimH": fimh.call}, loci_found=0, loci_total=1,
            note=(f"no fimH allele was called ({fimh.call}), so the subclone "
                  f"label cannot be derived. An ST131 isolate with no "
                  f"fimH allele would be unusual: read the MLST call"))
    entry = ST131_FIMH_CLADE.get(allele)
    if entry is None:
        return TypingResult(
            scheme="ST131_subclone", call="not applicable", lineage="-",
            alleles={"fimH": fimh.call}, loci_found=0, loci_total=1,
            note=(f"fimH allele {fimh.call} is not one of the ST131 clade "
                  f"alleles (fimH41 = clade A, fimH22 = clade B, "
                  f"fimH30 = clade C), so no ST131 subclone is implied. This "
                  f"says nothing about the ST: read the MLST call."))
    clade, meaning = entry
    subclone = f"H{allele[4:]}"
    note = (f"conditional on the ST: {subclone} is an ST131 subclone label and "
            f"this module cannot determine the ST, so read Hydra's Achtman ST. "
            f"If it is not 131, {allele} carries no subclone meaning and this "
            f"row should be ignored - {allele} is a common allele outside "
            f"ST131.")
    if clade == "C":
        note += (" The C1 (H30R) versus C2 (H30Rx) split is not resolvable "
                 "from fimH: it is defined on core-genome SNPs (Price et al. "
                 "2013), so both remain candidates. C2/H30Rx is the subclade "
                 "that carries blaCTX-M-15 - check the AMR calls, which "
                 "support but do not define it.")
    return TypingResult(
        scheme="ST131_subclone", call=f"{subclone} (clade {clade}) if ST131",
        lineage=meaning, alleles={"fimH": fimh.call, "clade": clade},
        loci_found=1, loci_total=1, note=note)


# ------------------------------------------------------------------- 3. stx
@dataclass
class _StxHit:
    subunit: str        # "A" or "B"
    subtype: str        # "1a", "2c", ...
    stx_class: str      # "1a".."2o", with 2a/2c/2d collapsed to "2"
    reference: str      # the reference accession
    contig: str
    lo: int
    hi: int
    strand: str
    nident: int
    reflen: int
    coverage_pct: float
    identity_pct: float
    residues: dict[int, str]


def _parse_stx_header(sseqid: str) -> tuple[str, str, str] | None:
    """``AAA98151.1|stxA1a|stxA1a`` -> ("A", "1a", "1")."""
    parts = sseqid.split("|")
    if len(parts) < 2:
        return None
    family = parts[1]
    match = re.fullmatch(r"stx([AB])(\d[a-z]?)", family)
    if match is None:
        return None
    subunit, subtype = match.group(1), match.group(2)
    stx_class = "2" if subtype in ("2a", "2c", "2d") else subtype
    return subunit, subtype, stx_class


def _residues_from_hsp(hsp: Hsp, wanted: tuple[int, ...]) -> dict[int, str]:
    """Query residues aligned to the given 1-based reference positions."""
    out: dict[int, str] = {}
    if not hsp.sseq or not hsp.qseq:
        return out
    ref_pos = min(hsp.sstart, hsp.send)
    target = set(wanted)
    for sub_aa, query_aa in zip(hsp.sseq, hsp.qseq):
        if sub_aa == "-":
            continue
        if ref_pos in target:
            out[ref_pos] = query_aa
        ref_pos += 1
    return out


def _stx(assembly: Path, root: Path, workdir: Path, threads: int) -> TypingResult:
    fasta = _require_db(root, "stxtyper/stx.prot", "the StxTyper reference proteins")
    # A stx-negative isolate - which is most isolates - is the no-hit case
    # that used to segfault blastx with several threads; blast() now pins a
    # translated search to one thread itself (engines.blast.safe_threads), so
    # the requested thread count is passed straight through.
    hsps = blast("blastx", assembly, fasta, workdir / "ecoli.stx.tsv",
                 threads=threads, evalue=1e-10, max_target_seqs=10000,
                 with_alignment=True)
    hsps = [h for h in hsps if h.pident >= STX_MIN_IDENTITY]
    if not hsps:
        return TypingResult(
            scheme="stx", call="none detected", lineage="-", loci_found=0,
            loci_total=0,
            note=("no stxA or stxB sequence found at 80 % amino-acid identity "
                  "or better against the 160 StxTyper reference proteins; this "
                  "isolate is not a Shiga-toxin producer as far as the "
                  "assembly shows"))

    merged = merge_hsps(hsps, query_scale=3)
    # One reference per place on the contig: the 160 references are many
    # near-identical alleles of the same handful of genes.
    kept = deduplicate(
        merged,
        key_span=lambda h: (h.query_start, h.query_end),
        key_seq=lambda h: h.qseqid,
        key_score=lambda h: (h.match_score, h.bitscore),
        key_tiebreak=lambda h: h.sseqid,
        overlap_fraction=0.5,
    )
    by_span: dict[tuple[str, int, int], Hsp] = {}
    for hsp in hsps:
        by_span.setdefault((hsp.qseqid, hsp.sseqid), hsp)
        if hsp.bitscore > by_span[(hsp.qseqid, hsp.sseqid)].bitscore:
            by_span[(hsp.qseqid, hsp.sseqid)] = hsp

    hits: list[_StxHit] = []
    for hit in kept:
        parsed = _parse_stx_header(hit.sseqid)
        if parsed is None:
            continue
        subunit, subtype, stx_class = parsed
        wanted = STX_A_POS if subunit == "A" else (STX_B_POS,)
        source = by_span.get((hit.qseqid, hit.sseqid))
        residues = _residues_from_hsp(source, wanted) if source else {}
        hits.append(_StxHit(
            subunit=subunit, subtype=subtype, stx_class=stx_class,
            reference=hit.sseqid.split("|")[0], contig=hit.qseqid,
            lo=hit.query_start, hi=hit.query_end, strand=hit.strand,
            nident=int(round(hit.identity_pct * hit.covered_bases / 100.0)),
            reflen=hit.slen, coverage_pct=hit.coverage_pct,
            identity_pct=hit.identity_pct, residues=residues,
        ))

    operons, orphans = _pair_stx(hits)
    calls: list[str] = []
    alleles: dict[str, str] = {}
    notes: list[str] = []
    for index, (a_hit, b_hit, gap) in enumerate(operons, start=1):
        label, status, note = _call_stx_operon(a_hit, b_hit, gap)
        tag = f"operon{index}"
        calls.append(label if status == "COMPLETE" else f"{label} [{status}]")
        combined = _combined_identity(a_hit, b_hit)
        alleles[f"{tag}_call"] = label
        alleles[f"{tag}_operon"] = status
        alleles[f"{tag}_identity_pct"] = f"{combined * 100:.2f}"
        alleles[f"{tag}_stxA_nearest_ref"] = (
            f"stxA{a_hit.subtype} ({a_hit.reference}) "
            f"{a_hit.identity_pct:.1f}% id / {a_hit.coverage_pct:.1f}% cov")
        alleles[f"{tag}_stxB_nearest_ref"] = (
            f"stxB{b_hit.subtype} ({b_hit.reference}) "
            f"{b_hit.identity_pct:.1f}% id / {b_hit.coverage_pct:.1f}% cov")
        alleles[f"{tag}_location"] = (f"{a_hit.contig}:{min(a_hit.lo, b_hit.lo)}-"
                                      f"{max(a_hit.hi, b_hit.hi)}({a_hit.strand})")
        alleles[f"{tag}_intergenic_bp"] = str(gap)
        if a_hit.stx_class == "2" and b_hit.stx_class == "2":
            # These three residues only mean anything inside the generalised
            # stx2a/2c/2d class; StxA1 is 316 residues long and has no 319.
            parts = [f"StxA{p}={a_hit.residues.get(p, '?')}" for p in STX_A_POS]
            parts.append(f"StxB{STX_B_POS}={b_hit.residues.get(STX_B_POS, '?')}")
            alleles[f"{tag}_subtyping_residues"] = " ".join(parts)
            if a_hit.subtype != b_hit.subtype:
                notes.append(
                    f"{tag}: the nearest stxA reference is {a_hit.subtype} and "
                    f"the nearest stxB reference is {b_hit.subtype}. That is "
                    f"expected and not a conflict: stx2a, stx2c and stx2d are "
                    f"one class whose references cross-match, and the subtype "
                    f"is decided by the three residues above, not by which "
                    f"reference scored highest")
        if note:
            notes.append(f"{tag}: {note}")
    for hit in orphans:
        tag = f"{'stxA' if hit.subunit == 'A' else 'stxB'}_orphan"
        calls.append(f"stx{hit.stx_class[0]} partial")
        alleles[f"{tag}_nearest"] = (f"stx{hit.subunit}{hit.subtype} "
                                     f"{hit.identity_pct:.1f}% id / "
                                     f"{hit.coverage_pct:.1f}% cov")
        alleles[f"{tag}_location"] = f"{hit.contig}:{hit.lo}-{hit.hi}({hit.strand})"
        notes.append(
            f"a lone stx{hit.subunit} subunit was found with no partner "
            f"subunit within {STX_INTERGENIC_MAX} bp. This is what a "
            f"prophage-interrupted locus, a contig break through the operon or "
            f"a cryptic remnant looks like; no subtype is called from half an "
            f"operon")

    if not calls:
        return TypingResult(
            scheme="stx", call="not determined", lineage="-", alleles=alleles,
            loci_found=0, loci_total=0,
            note=("stx-like sequence was found but no stxA/stxB pair could be "
                  "assembled into an operon; see the alleles for what matched"))
    label = " + ".join(sorted(calls))
    meanings = [STX_MEANING[c] for c in sorted({c for c in calls if c in STX_MEANING})]
    return TypingResult(
        scheme="stx", call=label, lineage="; ".join(meanings) or "-",
        alleles=alleles, loci_found=2 * len(operons) + len(orphans),
        loci_total=2 * len(operons) + len(orphans),
        note="; ".join(notes))


def _pair_stx(hits: list[_StxHit]) -> tuple[list[tuple[_StxHit, _StxHit, int]],
                                            list[_StxHit]]:
    """Pair stxA with the stxB immediately downstream of it on the same strand."""
    a_hits = sorted((h for h in hits if h.subunit == "A"), key=lambda h: (h.contig, h.lo))
    b_hits = sorted((h for h in hits if h.subunit == "B"), key=lambda h: (h.contig, h.lo))
    used: set[int] = set()
    operons: list[tuple[_StxHit, _StxHit, int]] = []
    for a_hit in a_hits:
        best: tuple[int, int] | None = None
        for index, b_hit in enumerate(b_hits):
            if index in used or b_hit.contig != a_hit.contig:
                continue
            if b_hit.strand != a_hit.strand:
                continue
            gap = (b_hit.lo - a_hit.hi - 1) if a_hit.strand == "+" \
                else (a_hit.lo - b_hit.hi - 1)
            # A partial operon may have lost sequence, so allow twice the
            # reference maximum before giving up, as StxTyper does (PD-4897).
            if 0 <= gap <= STX_INTERGENIC_MAX * 2 and (best is None or gap < best[1]):
                best = (index, gap)
        if best is not None:
            used.add(best[0])
            operons.append((a_hit, b_hits[best[0]], best[1]))
    orphans = [h for h in a_hits if h not in [o[0] for o in operons]]
    orphans += [b for index, b in enumerate(b_hits) if index not in used]
    return operons, orphans


def _combined_identity(a_hit: _StxHit, b_hit: _StxHit) -> float:
    total = a_hit.reflen + b_hit.reflen
    if not total:
        return 0.0
    return (a_hit.nident + b_hit.nident) / total


def _call_stx_operon(a_hit: _StxHit, b_hit: _StxHit,
                     gap: int) -> tuple[str, str, str]:
    """Subtype one stxA/stxB pair by the StxTyper rules. Returns (call, status, note)."""
    complete = a_hit.coverage_pct >= 100.0 and b_hit.coverage_pct >= 100.0
    superclass = a_hit.stx_class[0]
    if not complete:
        missing = []
        if a_hit.coverage_pct < 100.0:
            missing.append(f"stxA covers {a_hit.coverage_pct:.1f} % of its reference")
        if b_hit.coverage_pct < 100.0:
            missing.append(f"stxB covers {b_hit.coverage_pct:.1f} % of its reference")
        return (
            f"stx{superclass} partial", "PARTIAL",
            ("not subtyped: " + " and ".join(missing) + ". A subtype is only "
             "called when both subunits align over the whole reference coding "
             "sequence, because the subtype-defining residues sit near the ends "
             "of the proteins. The nearest references are "
             f"stx{superclass}{a_hit.subtype[1:]} (A) and "
             f"stx{superclass}{b_hit.subtype[1:]} (B); treat this as "
             "stx-positive, subtype undetermined - a truncated, "
             "phage-interrupted or contig-broken locus"))
    if a_hit.stx_class != b_hit.stx_class:
        return (f"stx{superclass} novel", "COMPLETE_NOVEL",
                (f"the two subunits match different stx types "
                 f"(stxA{a_hit.subtype}, stxB{b_hit.subtype}); this is either "
                 f"a mosaic operon or a novel subtype, and no established "
                 f"subtype name applies"))
    if gap > STX_INTERGENIC_MAX:
        return (f"stx{superclass} novel", "COMPLETE_NOVEL",
                (f"the stxA-stxB intergenic region is {gap} bp, longer than the "
                 f"{STX_INTERGENIC_MAX} bp seen in any reference operon, so "
                 f"these are probably not one operon"))
    combined = _combined_identity(a_hit, b_hit)
    threshold = STX_IDENTITY.get(a_hit.stx_class, 0.98)
    if combined < threshold:
        return (f"stx{superclass} novel", "COMPLETE_NOVEL",
                (f"combined operon identity {combined * 100:.2f} % is below the "
                 f"{threshold * 100:.1f} % cut-off for "
                 f"stx{a_hit.stx_class}, so this is a novel subtype rather "
                 f"than an established one"))
    if a_hit.stx_class != "2":
        return f"stx{a_hit.stx_class}", "COMPLETE", ""
    # stx2a, stx2c and stx2d are one class and are separated on three residues.
    res_a1 = a_hit.residues.get(STX_A_POS[0], "")
    res_a2 = a_hit.residues.get(STX_A_POS[1], "")
    res_b = b_hit.residues.get(STX_B_POS, "")
    for subtype, first, second, third in STX2ACD_RULES:
        if res_a1 in first and res_a2 in second and res_b in third:
            return f"stx{subtype}", "COMPLETE", ""
    return ("stx2 novel", "COMPLETE_NOVEL",
            (f"the operon is a complete stx2 of the stx2a/2c/2d class but the "
             f"discriminating residues (StxA {STX_A_POS[0]}={res_a1 or '?'}, "
             f"{STX_A_POS[1]}={res_a2 or '?'}; StxB {STX_B_POS}={res_b or '?'}) "
             f"match none of the three published patterns, so it is stx2 of an "
             f"undefined subtype"))


# ------------------------------------- 4. pathotype, intimin, ExPEC/UPEC panel
@dataclass
class _Gene:
    name: str
    identity_pct: float
    coverage_pct: float
    contig: str
    span: str
    reference: str


def _virulence_genes(assembly: Path, root: Path, workdir: Path,
                     threads: int) -> dict[str, _Gene]:
    fasta = _require_db(root, "virulencefinder/virulence_ecoli.fsa",
                        "the VirulenceFinder E. coli database")
    hsps = blast("blastn", assembly, fasta, workdir / "ecoli.vf.tsv",
                 threads=threads, evalue=1e-20, task="megablast", word_size=20,
                 perc_identity=VF_MIN_IDENTITY - 5.0, max_target_seqs=20000)
    best: dict[str, _Gene] = {}
    for hit in merge_hsps(hsps):
        if hit.identity_pct < VF_MIN_IDENTITY or hit.coverage_pct < VF_MIN_COVERAGE:
            continue
        name = hit.sseqid.split(":")[0]
        gene = _Gene(name=name, identity_pct=hit.identity_pct,
                     coverage_pct=hit.coverage_pct, contig=hit.qseqid,
                     span=f"{hit.query_start}-{hit.query_end}({hit.strand})",
                     reference=hit.sseqid)
        current = best.get(name)
        if current is None or (gene.identity_pct, gene.coverage_pct) > \
                (current.identity_pct, current.coverage_pct):
            best[name] = gene
    return best


def _group_present(genes: dict[str, _Gene], group: str) -> list[str]:
    """Names of the genes of *group* that were found."""
    prefixes = MARKER_GROUPS[group]
    return sorted(name for name in genes
                  if any(name.startswith(prefix) for prefix in prefixes))


#: How many gene names one evidence cell may list before it is summarised.
#: eae has 65 reference alleles that all cross-match, and listing them produced
#: a 1200-character cell that named nothing - the same failure the rmpA2 marker
#: collapse in ``typing.lineage`` exists to prevent.
_EVIDENCE_LIMIT = 6


def _summarise(names: list[str], genes: dict[str, _Gene]) -> str:
    """One evidence cell for a marker group, best matches first."""
    if not names:
        return "-"
    ranked = sorted(
        names,
        key=lambda n: (genes[n].identity_pct * genes[n].coverage_pct,
                       genes[n].identity_pct),
        reverse=True)
    shown = ranked[:_EVIDENCE_LIMIT]
    cells = [f"{n} ({genes[n].identity_pct:.1f}%id/{genes[n].coverage_pct:.0f}%cov)"
             for n in shown]
    hidden = len(ranked) - len(shown)
    if hidden:
        cells.append(f"+{hidden} further allele matches")
    return ", ".join(cells)


def _pathotype(genes: dict[str, _Gene], stx: TypingResult,
               shigella: TypingResult) -> TypingResult:
    """The diarrhoeagenic/extraintestinal pathotype, with the rule applied."""
    found = {group: _group_present(genes, group) for group in MARKER_GROUPS}
    stx_positive = stx.call not in ("none detected", "not determined", "-", "")
    labels: list[str] = []
    rules: list[str] = []

    eae = bool(found["eae"])
    if stx_positive:
        if eae:
            labels.append("EHEC")
            rules.append("stx present and eae (LEE) present -> "
                         "enterohaemorrhagic E. coli")
        elif found["aggR"] or found["AAF"]:
            labels.append("STEC/EAEC hybrid")
            rules.append("stx present with aggR or an AAF adhesin and no eae -> "
                         "the enteroaggregative-haemorrhagic hybrid pathotype "
                         "of the 2011 O104:H4 outbreak")
        else:
            labels.append("STEC (eae-negative)")
            rules.append("stx present, eae absent -> Shiga toxin-producing "
                         "E. coli without the LEE")
    elif eae:
        if found["bfpA"]:
            labels.append("typical EPEC")
            rules.append("eae present, bfpA (EAF plasmid bundle-forming pilus) "
                         "present, stx absent -> typical EPEC")
        else:
            labels.append("atypical EPEC")
            rules.append("eae present, bfpA absent, stx absent -> atypical EPEC")

    if found["LT"] or found["ST"]:
        labels.append("ETEC")
        toxins = []
        if found["LT"]:
            toxins.append("LT (eltAB)")
        if found["ST"]:
            toxins.append("ST (estA)")
        rules.append(f"{' and '.join(toxins)} present -> enterotoxigenic E. coli")

    hybrid = any(label.startswith("STEC/EAEC") for label in labels)
    if found["aggR"] and not hybrid:
        # The hybrid label already says "EAEC"; adding a second EAEC row would
        # report one finding twice.
        labels.append("EAEC (typical)")
        rules.append("aggR present -> typical enteroaggregative E. coli")
    elif found["AAF"] and (found["aat"] or found["aai"]) and not hybrid:
        labels.append("EAEC (probable, aggR-negative)")
        rules.append("an AAF adhesin with aatA or aaiC but no aggR -> "
                     "aggR-negative EAEC; the aggR-negative form is of "
                     "contested pathogenicity")

    if found["ipaH"]:
        verdict = shigella.alleles.get("organism", "")
        if verdict == "Shigella":
            labels.append("Shigella")
            rules.append("ipaH present and a Shigella cluster-specific gene set "
                         "complete -> Shigella, not EIEC (see shigella_eiec)")
        elif verdict == "EIEC":
            labels.append("EIEC")
            rules.append("ipaH present and an EIEC cluster-specific gene set "
                         "complete -> enteroinvasive E. coli (see shigella_eiec)")
        else:
            labels.append("EIEC/Shigella")
            rules.append("ipaH present but no cluster-specific gene set is "
                         "complete, so EIEC and Shigella are not separated")

    expec_groups = [g for g in EXPEC_GROUPS if found[g]]
    is_expec = len(expec_groups) >= 2
    # DAEC is only worth raising when nothing else explains the Afa/Dr operon.
    # afa/draBC is itself one of the five ExPEC markers, so in a genome that
    # already meets the ExPEC rule the operon is an extraintestinal adhesin,
    # and calling that isolate "possible DAEC" would be wrong.
    if not labels and not is_expec and found["afa_dra"]:
        labels.append("possible DAEC")
        rules.append("an Afa/Dr adhesin operon with no other diarrhoeagenic "
                     "marker and no second ExPEC marker. DAEC is defined on "
                     "the afaE-III/daaE adhesin variant plus a "
                     "diffuse-adherence phenotype, neither of which an "
                     "assembly settles, so this is a candidate, not a call")
    if is_expec:
        labels.append("ExPEC")
        extra = ""
        if found["afa_dra"]:
            extra = ("; the afa/dra operon counted here is an extraintestinal "
                     "adhesin in this context, not evidence of DAEC")
        rules.append(f"{len(expec_groups)} of the 5 ExPEC marker groups present "
                     f"({', '.join(expec_groups)}); the rule is two or more of "
                     f"papAH/papC, sfa/focDE, afa/draBC, kpsM II, iutA{extra}")

    alleles = {group: _summarise(names, genes)
               for group, names in found.items() if names}
    alleles["ExPEC_marker_groups"] = f"{len(expec_groups)}/5"
    if not labels:
        return TypingResult(
            scheme="pathotype", call="no pathotype markers", lineage="-",
            alleles=alleles, loci_found=len(alleles), loci_total=len(MARKER_GROUPS),
            note=("none of the defining marker combinations for EPEC, EHEC, "
                  "STEC, ETEC, EAEC, EIEC, DAEC or ExPEC was found. This is "
                  "what a commensal or an opportunist without the classical "
                  "pathotype determinants looks like; it does not mean the "
                  "isolate is harmless, only that it carries no pathotype "
                  "signature"))
    return TypingResult(
        scheme="pathotype", call=" + ".join(labels),
        lineage="; ".join(rules), alleles=alleles,
        loci_found=len(alleles), loci_total=len(MARKER_GROUPS),
        note=("pathotype is assigned from marker combinations only; adherence "
              "phenotype and the clinical picture decide the final label"))


def _intimin(genes: dict[str, _Gene]) -> TypingResult:
    """The eae (intimin) subtype in the Greek-letter nomenclature, plus tir."""
    eae_hits = {name: gene for name, gene in genes.items()
                if name.startswith("eae")}
    tir_hits = {name: gene for name, gene in genes.items()
                if name.startswith("tir")}
    if not eae_hits:
        return TypingResult(
            scheme="intimin", call="eae absent", lineage="-",
            alleles={"tir": ", ".join(sorted(tir_hits)) or "-"},
            loci_found=0, loci_total=1,
            note=("no eae allele reached 90 % identity over 60 % of a "
                  "reference, so the locus of enterocyte effacement intimin "
                  "gene is absent and no attaching-and-effacing lesion is "
                  "predicted"))
    name, gene = max(eae_hits.items(),
                     key=lambda kv: (kv[1].identity_pct * kv[1].coverage_pct,
                                     kv[1].identity_pct))
    subtype = _eae_subtype(name)
    alleles = {
        "eae_allele": name,
        "identity_pct": f"{gene.identity_pct:.2f}",
        "coverage_pct": f"{gene.coverage_pct:.1f}",
        "contig": gene.contig,
        "contig_span": gene.span,
        "tir": ", ".join(sorted(tir_hits)) or "-",
        "runner_up_eae_alleles": _summarise(
            sorted(k for k in eae_hits if k != name), genes),
    }
    runners = sorted((kv for kv in eae_hits.items() if kv[0] != name),
                     key=lambda kv: kv[1].identity_pct, reverse=True)
    if runners and runners[0][1].identity_pct >= gene.identity_pct - 1.0:
        alleles["margin_to_next_allele"] = (
            f"{gene.identity_pct - runners[0][1].identity_pct:.2f} "
            f"percentage points over {_eae_subtype(runners[0][0])}")
    note = ""
    call = subtype
    if gene.identity_pct < EAE_SUBTYPE_IDENTITY:
        call = f"{subtype}*"
        note = (f"the best eae match is only {gene.identity_pct:.2f} % "
                f"identical to the {subtype} reference; the 65 intimin alleles "
                f"differ from one another by as little as 2 %, so below "
                f"{EAE_SUBTYPE_IDENTITY:.0f} % the Greek letter is a nearest "
                f"neighbour rather than a subtype call")
    if gene.coverage_pct < 95.0:
        note = (note + "; " if note else "") + (
            f"only {gene.coverage_pct:.1f} % of the reference allele is "
            f"covered, so the subtype rests on part of the gene")
    if not tir_hits:
        note = (note + "; " if note else "") + (
            "eae is present but no tir allele was found, which is unusual - "
            "tir sits in the same LEE operon, so check for a contig break")
    return TypingResult(
        scheme="intimin", call=call,
        lineage=("intimin subtype; determines the host-cell Tir binding "
                 "specificity and tracks with the EPEC/EHEC serotype"),
        alleles=alleles, loci_found=1, loci_total=1, note=note)


def _eae_subtype(name: str) -> str:
    """``eae-g01-gamma_1`` -> ``gamma1``."""
    parts = name.split("-", 2)
    if len(parts) < 3:
        return name
    return parts[2].replace("_", "")


def _expec_panel(genes: dict[str, _Gene]) -> TypingResult:
    """The extraintestinal/urinary virulence-factor repertoire."""
    found: dict[str, list[str]] = {}
    for category, members in VF_PANEL.items():
        hits = sorted(name for name in genes
                      if any(name.startswith(member) for member in members))
        if hits:
            found[category] = hits
    expec_groups = [g for g in EXPEC_GROUPS if _group_present(genes, g)]
    upec_four = [g for g in UPEC_FOUR
                 if any(name.startswith(g) for name in genes)]
    total = sum(len(v) for v in found.values())
    alleles = {category: ", ".join(hits) for category, hits in found.items()}
    alleles["ExPEC_groups_present"] = ", ".join(expec_groups) or "-"
    alleles["ExPEC_rule"] = f"{len(expec_groups)}/5 (>=2 required)"
    alleles["UPEC_four_gene_rule"] = f"{len(upec_four)}/4: " + \
        (", ".join(upec_four) or "none")
    calls = []
    if len(expec_groups) >= 2:
        calls.append(f"ExPEC ({len(expec_groups)}/5)")
    if len(upec_four) == 4:
        calls.append("UPEC by the four-gene rule")
    if calls:
        call = " + ".join(calls)
    else:
        call = (f"not ExPEC ({len(expec_groups)}/5 marker groups), "
                f"not UPEC by the four-gene rule ({len(upec_four)}/4); "
                f"{total} virulence factors found")
    note = ("the ExPEC rule (two or more of papAH/papC, sfa/focDE, afa/draBC, "
            "kpsM II, iutA) is a molecular definition of extraintestinal "
            "potential, not a prediction of infection; the four-gene rule "
            "(vat, fyuA, chuA, yfcV) identified efficient urinary-tract "
            "colonisers experimentally. Both are screens: a urinary isolate is "
            "UPEC by where it was grown, whatever the genotype")
    if not calls and total:
        note += ("; the factors listed are reported individually because a "
                 "repertoire below the thresholds still matters clinically")
    return TypingResult(
        scheme="ExPEC_VF", call=call,
        lineage=("extraintestinal/urinary virulence repertoire, grouped as the "
                 "urosepsis literature groups it"),
        alleles=alleles, loci_found=total,
        loci_total=sum(len(v) for v in VF_PANEL.values()), note=note)


# ---------------------------------------------------- 5. Shigella versus EIEC
def _shigella_eiec(assembly: Path, root: Path, workdir: Path,
                   threads: int) -> TypingResult:
    fasta = _require_db(root, "shigeifinder/genes.fasta",
                        "the ShigEiFinder cluster-specific genes")
    table = _require_db(root, "shigeifinder/genes.json",
                        "the ShigEiFinder cluster definitions")
    try:
        definitions = json.loads(table.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise HydraError(f"could not read {table}: {exc}") from exc

    hsps = blast("blastn", assembly, fasta, workdir / "ecoli.shigei.tsv",
                 threads=threads, evalue=1e-20, task="megablast",
                 perc_identity=80.0, max_target_seqs=20000)
    coverage: dict[str, float] = {}
    aligned: dict[str, int] = {}
    for hit in merge_hsps(hsps):
        name = hit.sseqid.split(":")[-1] if ":" in hit.sseqid else hit.sseqid
        name = name.split()[0]
        if hit.coverage_pct > coverage.get(name, -1.0):
            coverage[name] = hit.coverage_pct
            aligned[name] = hit.covered_bases

    def present(gene: str) -> bool:
        cov = coverage.get(gene)
        if cov is None:
            return False
        if gene == "ipaH":
            return cov > SHIGEI_IPAH_COVERAGE
        minimum = SHIGEI_MIN_ALIGNED.get(gene)
        if minimum is not None:
            return aligned.get(gene, 0) >= minimum
        return cov >= SHIGEI_MIN_COVERAGE

    complete = [name for name, entry in definitions.items()
                if name != "sporadic"
                and all(present(gene) for gene in entry["cluster-genes"])]
    ipah = present("ipaH")
    alleles: dict[str, str] = {
        "ipaH": "+" if ipah else "-",
        "clusters_complete": ", ".join(sorted(complete)) or "-",
    }
    for name in sorted(complete):
        genes = definitions[name]["cluster-genes"]
        alleles[f"{name}_genes"] = ", ".join(
            f"{g} {coverage.get(g, 0.0):.0f}% cov" for g in genes)

    if len(complete) == 1:
        cluster = complete[0]
        organism, description = SHIGEI_CLUSTERS.get(
            cluster, ("unknown", "cluster not in this module's table"))
        alleles["organism"] = organism
        note = ""
        if not ipah:
            note = ("the cluster-specific genes are complete but ipaH was not "
                    "found: either the invasion plasmid has been lost in "
                    "culture (common in stored Shigella) or this is a "
                    "non-invasive relative of the cluster")
        if organism == "Shigella":
            note = (note + "; " if note else "") + (
                "Hydra's species call will usually read this genome as "
                "Escherichia coli, because Shigella is a clade within E. coli "
                "by any whole-genome measure. This call, not the species call, "
                "is the one to report and to act on clinically")
        return TypingResult(
            scheme="shigella_eiec", call=f"{cluster} ({organism})",
            lineage=description, alleles=alleles,
            loci_found=len(definitions[cluster]["cluster-genes"]),
            loci_total=len(definitions[cluster]["cluster-genes"]), note=note)

    if len(complete) > 1:
        # ShigEiFinder resolves a handful of known overlaps; anything else is
        # reported as ambiguous rather than picked arbitrarily.
        resolved = _resolve_cluster_overlap(complete)
        if resolved:
            organism, description = SHIGEI_CLUSTERS.get(
                resolved, ("unknown", "cluster not in this module's table"))
            alleles["organism"] = organism
            return TypingResult(
                scheme="shigella_eiec", call=f"{resolved} ({organism})",
                lineage=description, alleles=alleles,
                loci_found=len(definitions[resolved]["cluster-genes"]),
                loci_total=len(definitions[resolved]["cluster-genes"]),
                note=(f"the gene sets of {', '.join(sorted(complete))} are all "
                      f"complete; ShigEiFinder's published precedence rule "
                      f"resolves this to {resolved}"))
        alleles["organism"] = "ambiguous"
        return TypingResult(
            scheme="shigella_eiec", call=",".join(sorted(complete)),
            lineage="more than one cluster gene set is complete",
            alleles=alleles, loci_found=len(complete), loci_total=len(complete),
            note=("several cluster-specific gene sets are complete and no "
                  "published rule separates them; a mixed culture or a "
                  "contaminated assembly does this. Not resolved"))

    alleles["organism"] = "Escherichia coli" if not ipah else "unclustered"
    if ipah:
        return TypingResult(
            scheme="shigella_eiec", call="Shigella/EIEC, unclustered",
            lineage="invasive, but outside the defined clusters",
            alleles=alleles, loci_found=0, loci_total=0,
            note=("ipaH is present, so this carries the Shigella/EIEC "
                  "invasion determinant, but no cluster-specific gene set is "
                  "complete. Zhang et al. describe such isolates as "
                  "unclustered Shigella/EIEC; it may also be an E. coli that "
                  "has acquired ipaH alone"))
    return TypingResult(
        scheme="shigella_eiec", call="not Shigella/EIEC", lineage="-",
        alleles=alleles, loci_found=0, loci_total=0,
        note=("ipaH absent and no cluster-specific gene set complete, so this "
              "is neither Shigella nor enteroinvasive E. coli"))


#: The overlaps ShigEiFinder resolves by a fixed precedence, from
#: ``determine_cluster`` in the published implementation.
_CLUSTER_PRECEDENCE: tuple[tuple[frozenset[str], str], ...] = (
    (frozenset({"C1", "CSB12"}), "CSB12"),
    (frozenset({"C5", "C3"}), "C3"),
    (frozenset({"C5", "C8"}), "C8"),
    (frozenset({"C1", "C2"}), "C2"),
    (frozenset({"C1", "CSD1"}), "CSD1"),
    (frozenset({"C2", "CSS"}), "C2"),
    (frozenset({"CSB12", "CSD1"}), "CSB12"),
)


def _resolve_cluster_overlap(complete: list[str]) -> str:
    found = frozenset(complete)
    for pair, winner in _CLUSTER_PRECEDENCE:
        if pair == found:
            return winner
    return ""


# ----------------------------------------------------------------------- run
def run(assembly: Path, db_dir: Path, *, threads: int = 1,
        workdir: Path | None = None) -> list[TypingResult]:
    """Type one assembly, returning one :class:`TypingResult` per scheme."""
    assembly = Path(assembly)
    root = Path(db_dir) / KEY
    work = Path(workdir) if workdir is not None else assembly.parent
    work.mkdir(parents=True, exist_ok=True)

    contigs = list(read_fasta(assembly))
    if not contigs:
        raise HydraError(f"{assembly} holds no sequences")

    def attempt(scheme: str, call, *args):
        """Run one scheme; a database it needs but has not got is reported.

        Returns the :class:`TypingResult`, or a "not determined" one naming the
        missing database. Any other failure is left to propagate: that is a
        bug, and a bug must not be dressed up as a typing result.
        """
        try:
            return call(*args)
        except _MissingDb as exc:
            LOG.info("%s: %s", KEY, exc)
            return TypingResult(scheme=scheme, call="not determined", lineage="-",
                                note=str(exc))

    results: list[TypingResult] = []
    results.append(_phylogroup(contigs))

    fimh = attempt("fimH", _fimh, assembly, root, work, threads)
    results.append(fimh)
    results.append(_st131_subclone(fimh))

    stx = attempt("stx", _stx, assembly, root, work, threads)
    results.append(stx)

    shigella = attempt("shigella_eiec", _shigella_eiec, assembly, root, work, threads)
    try:
        genes = _virulence_genes(assembly, root, work, threads)
    except _MissingDb as exc:
        LOG.info("%s: %s", KEY, exc)
        for scheme in ("pathotype", "intimin", "ExPEC_VF"):
            results.append(TypingResult(scheme=scheme, call="not determined",
                                        lineage="-", note=str(exc)))
    else:
        results.append(_pathotype(genes, stx, shigella))
        results.append(_intimin(genes))
        results.append(_expec_panel(genes))
    results.append(shigella)
    return results
