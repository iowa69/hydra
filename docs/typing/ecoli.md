# *Escherichia coli* / *Shigella* typing

Module: `hydra_amr.typing.organisms.ecoli` (`KEY = "ecoli"`).
Triggered by the species names in `SPECIES` (the four clinically relevant
*Shigella* species included) and by the PubMLST schemes `ecoli_achtman_4`
(Achtman) and `ecoli` (Pasteur).

Everything below runs from a draft assembly using BLAST+ and pure Python only.

---

## 1. What is called, and what each call means

### Relation to what Hydra already does

Hydra already reports, from other modules, the **O:H serotype** (EcOH /
`serotype`), the **Kleborate pathovar marker set** (`pathovar_markers`), the
**Achtman and Pasteur STs**, and the yersiniabactin / colibactin / aerobactin /
salmochelin siderophore schemes. None of that is repeated here.

| This module's call | Relation to the existing call |
|---|---|
| `phylogroup` | New. Nothing in Hydra assigned a Clermont phylogroup. |
| `fimH` | New. Completes the clinically used subclone naming, which needs ST **and** fimH. |
| `ST131_subclone` | New, and explicitly conditional on Hydra's Achtman ST. |
| `stx` | **Supersedes** the stx entries in `pathovar_markers`, which report presence of a Kleborate marker without a subtype. A subtype is the clinically actionable result. |
| `pathotype` | **Supersedes** `pathovar_markers` as a *conclusion*: `pathovar_markers` lists markers, this applies the published rules and names the pathotype and the rule. The marker list remains useful as raw evidence. |
| `intimin` | New. `pathovar_markers` reports eae presence; this names the variant. |
| `ExPEC_VF` | New. Partly overlaps the siderophore schemes (`iucC`/`iutA` is aerobactin); the siderophore schemes remain the allele-level answer, this is the clinical repertoire. |
| `shigella_eiec` | New, and it **contradicts the species call by design**. *Shigella* is a clade inside *E. coli*, so Hydra's species call reports *Shigella* genomes as *Escherichia coli*. This call is the one to report. |

### `phylogroup` — A, B1, B2, C, D, E, F, G, cryptic clade
The phylogroup is the coarsest and most load-bearing division of the species.
It predicts the kind of disease an isolate causes before any virulence gene is
read: B2 and D hold the extraintestinal lineages (ST131, ST73, ST95, ST69), E
holds the O157:H7 EHEC lineage, A and B1 hold most commensals, most ETEC and
most non-O157 STEC, and G is a recently defined group with high virulence and
resistance potential. The *Escherichia* cryptic clades are not *E. coli* at all
and are reported as such.

The call is reported **with the marker pattern**, not as a bare letter:
`arpA+ chuA+ yjaA- TspE4-` and the amplicon size and coordinate of each
product, so the call can be checked against the published table.

### `fimH` — the FimTyper allele
fimH subtyping separates lineages that MLST cannot. With the ST it produces the
subclone names in clinical and epidemiological use: ST131-*fimH*30 is the H30
subclone responsible for the global rise in ESBL-producing *E. coli*. The allele
is reported with its identity, the reference span covered and the contig
coordinates; an inexact match is marked `*` and named as a novel or untyped
allele, never passed off as the nearest allele.

### `ST131_subclone` — H30 / H22 / H41 and the clade
Reported as **`H30 (clade C) if ST131`**: the label is an ST131 label and this
module cannot determine the ST. The three alleles the clade papers name are
implemented (fimH41 = clade A, fimH22 = clade B, fimH30 = clade C); any other
allele returns `not applicable` rather than a guessed clade.

The **C1 (H30R) / C2 (H30Rx)** split is deliberately **not** attempted. H30Rx
was defined on core-genome SNPs (Price *et al.* 2013), not on any single marker
determinable by alignment to a gene database, so both remain candidates and the
note says so. `bla`CTX-M-15, which Hydra reports separately, supports C2 but
does not define it.

### `stx` — stx1 and stx2 subtypes
The subtype, not the presence of stx, drives the HUS risk assessment: **stx2a**
is the subtype most strongly associated with HUS and bloody diarrhoea, stx2d is
comparable in *eae*-negative strains, stx1a is a considerably lower risk, and
stx2b/2e/2f/2g are associated with mild or non-human disease. The module reports
the subtype, the combined operon identity, the nearest reference allele for each
subunit with identity and coverage, the operon's coordinates, the intergenic
distance and — for the stx2a/2c/2d class — the three discriminating residues it
read.

A **partial or interrupted locus is never given a subtype.** A pair covering
less than the whole reference coding sequence is reported as `stx1 partial` /
`stx2 partial` with the reason, because the subtype-defining residues sit near
the ends of the proteins. A lone stxA or stxB with no partner within the
intergenic limit is reported as an orphan subunit — what a prophage-interrupted
locus, a contig break through the operon or a cryptic remnant looks like.

### `pathotype` — EPEC / EHEC / STEC / ETEC / EAEC / EIEC / DAEC / ExPEC
The pathotype, and **the rule that produced it**, in the `lineage` field. Each
rule is stated in gene terms so it can be checked (see §2).

### `intimin` — the eae subtype
The Greek-letter intimin variant (gamma1, beta1a, epsilon1, theta, ...), which
determines Tir binding specificity and tracks with the EPEC/EHEC serotype:
gamma1 with O157:H7, beta1 with O26:H11, epsilon with O103:H2. Reported with the
identity, the coverage, the coordinates, the runner-up alleles and the margin
over the runner-up where it is thin, plus `tir`.

### `ExPEC_VF` — the extraintestinal / urinary repertoire
The virulence factors, grouped as the urosepsis literature groups them
(adhesins, toxins, iron uptake, protectins and capsule), plus the two published
rules: the five-marker **ExPEC** definition and the four-gene **UPEC** rule.
Both are reported as screens, with the explicit caveat that a urinary isolate
is UPEC by where it was grown, whatever the genotype.

### `shigella_eiec` — *Shigella* cluster vs EIEC
The ShigEiFinder cluster, the *Shigella*-or-EIEC verdict and what the cluster
contains. This is the call Hydra's species identification cannot make.

---

## 2. Method, exact rules and thresholds

### Phylogroup — in-silico Clermont PCR
Implemented as in-silico PCR, reproducing EzClermont primer for primer, so the
result matches the reference tool rather than approximating it.

Primers (Clermont *et al.* 2013 for the quadruplex, group C and group E
reactions; Clermont *et al.* 2019 for *ybgD*), with expected amplicon sizes:

| reaction | forward | reverse | product |
|---|---|---|---|
| arpA | `AAYRCYATTCGCCAGCTTGC` | `TCTCCMCATACYGYACGCTA` | 400 bp |
| chuA | `ATGGTACYGGRCGAACCAAC` | `TRCCRCCAGTRCCAAAGACA` | 288 bp |
| yjaA | `YAAACKTGAAGTGTCAGGAG` | `ARTRCGTTCCTCAACCTGTG` | 211 bp |
| TspE4.C2 | `CACKATTYGTAAGRYCATCC` | `AGTTTATCGCTGCGGGTCGC` | 152 bp |
| arpA group E | `RATKCMATYTTGTCRAAATATGCC` | `GAAARKRAAAADAMYYYYCAAGAG` | 301 bp |
| trpA group C | `AGTTYTAYGCCSVRWGCGAG` | `TCWGYDCYVGTYACGCCC` | 219 bp |
| ybgD group G | `GTTGACTAARCGYAGGTCGA` | `KATGYDGCYGATKAAGGATC` | 177 bp |
| trpBA control | `CGGSGATAAAGAYATYTTCAC` | `GCAACGYGSCBWKRCGGAAG` | 489 bp |

IUPAC codes are expanded to character classes; both primers are searched on
both orientations of every contig; a product is positive when they are found
the expected distance apart **within 20 bp**.

Quadruplex decision table (`arpA, chuA, yjaA, TspE4.C2`):

| pattern | call | follow-up |
|---|---|---|
| `+ - - -` | **A** | — |
| `+ - - +` | **B1** | — |
| `- + + -` / `- + + +` | **B2** | — |
| `+ - + -` | A/C | trpA group C: + → **C**, − → **A** |
| `+ + - -` / `+ + - +` | D/E | arpA group E: + → **E**, − → **D** |
| `- + - -` | F/G | ybgD: + → **G**, − → **F** |
| `- + - +` | B2/G | ybgD: + → **G**, − → **B2** |
| `+ + + -` | E/cryptic | arpA group E: + → **E**, − → **cryptic clade** |
| `- - + -` | **cryptic clade** | — |
| `- - - -` | **not determined** | cryptic clade or non-matching assembly |
| `- - - +`, `- - + +`, `+ - + +`, `+ + + +` | **U (unknown)** | not in the published table |

The group C and group E reactions are read against the **trpBA internal
control**. If the control gives no product, the module refuses to refine and
reports e.g. `D/E (control failed)` — a negative that cannot be distinguished
from a failed reaction is not a negative.

**Partial products.** When the assembly has more than four contigs, a single
primer sitting closer to the end of its contig than the expected product length
counts as a partial product, as EzClermont does. Without this, a contig break
through an amplicon silently loses a marker and changes the phylogroup.

### fimH
`blastn -task megablast -perc_identity 90`, HSPs merged on the reference axis,
candidate alleles filtered to ≥ 60 % reference coverage. The best allele is the
one with the highest *fraction of the reference matched identically*
(coverage × identity), ties broken on the longer reference — ranking by bit
score picks a short allele that is a substring of the right one. An allele is
named only at **100 % identity over 100 % of the reference**; anything less is
marked `*` and described as novel or untyped. Several alleles matching exactly
and inseparably are all reported, joined by `/`.

Three distinct negative results are distinguished: `fimH absent`/`fimH not
found` (nothing in the assembly aligns at all), `not determined` (a fimH-like
fragment below 60 % coverage, i.e. a contig break or an insertion), and a
definite allele.

### stx — the StxTyper algorithm
`blastx` of the assembly against the 160 StxTyper reference proteins, HSPs
≥ 80 % amino-acid identity, merged on the reference axis with `query_scale=3`,
then one reference per place on the contig.

An stxA hit is paired with the nearest stxB **on the same contig, on the same
strand, downstream**, with an intergenic gap of 0 to 72 bp (twice the reference
maximum, as StxTyper allows for partial operons).

A subtype is declared only when **all** of:

1. both subunits align over **100 %** of their reference coding sequences
   (`COMPLETE`);
2. both subunits are of the same stx class, where **stx2a, stx2c and stx2d
   count as one class** (`"2"`);
3. the intergenic region is **≤ 36 bp** (the reference-set maximum + 2);
4. the combined identity (identical residues of A + B divided by the summed
   reference lengths) is at least:
   **0.983** for stx1a/1c/1d/1e, **0.985** for stx2k and stx2l, **0.980** for
   every other type.

Inside the generalised stx2 class the subtype is then read from three residues —
StxA **313** and **319** and StxB **35**, 1-based in the reference protein:

| StxA 313 | StxA 319 | StxB 35 | subtype |
|---|---|---|---|
| F or S | K or E | D | **stx2a** |
| F | K or E | N | **stx2c** |
| S | E | N | **stx2d** |

Rules are tried in that order. Any other pattern is `stx2 novel`
(`COMPLETE_NOVEL`) with the observed residues printed, never forced into one of
the three.

Mixed-class, over-long-intergenic and below-threshold operons are all
`COMPLETE_NOVEL` with the reason. Partial operons and orphan subunits are
reported as such.

### Virulence genes, pathotype, intimin and the ExPEC panel
One `blastn -task megablast` pass against the VirulenceFinder *E. coli*
database. A gene is **present** at **≥ 90 % identity over ≥ 60 % of the
reference**, which is the VirulenceFinder web default. An eae allele is named as
a subtype only at **≥ 97 % identity** — the 65 intimin alleles differ from one
another by as little as 2 %, so a looser threshold names the wrong Greek
letter; below it the letter is marked `*` and described as a nearest neighbour.

Evidence cells list at most six gene names, best first, with
`+N further allele matches` for the rest. Listing all 57 cross-matching eae
alleles produced a 1200-character cell that named nothing.

Pathotype rules, each reported in the `lineage` field of the result:

| pathotype | rule |
|---|---|
| **EHEC** | stx present **and** eae present |
| **STEC/EAEC hybrid** | stx present **and** (aggR or an AAF adhesin) **and** eae absent |
| **STEC (eae-negative)** | stx present, eae absent, no EAEC marker |
| **typical EPEC** | eae present, **bfpA** present, stx absent |
| **atypical EPEC** | eae present, bfpA absent, stx absent |
| **ETEC** | LT (`eltIAB`/`eltIIAB`) and/or ST (`estah-STa`/`estap-STa`/`estb-STb`) |
| **EAEC (typical)** | **aggR** present |
| **EAEC (probable, aggR-negative)** | an AAF adhesin + (aatA or aaiC), aggR absent; flagged as a contested entity |
| **Shigella** / **EIEC** | ipaH present, decided by the `shigella_eiec` cluster call |
| **EIEC/Shigella** | ipaH present but no cluster gene set complete — not separated |
| **possible DAEC** | an Afa/Dr operon, **no** other DEC marker **and** the ExPEC rule not met; reported as a candidate because DAEC needs the afaE-III/daaE variant and a diffuse-adherence phenotype, neither of which an assembly settles |
| **ExPEC** | ≥ 2 of: papAH/papC, sfa/focDE, afa/draBC, kpsM II, iutA |

afa/draBC is itself one of the five ExPEC markers, so in a genome that already
meets the ExPEC rule the operon is counted as an extraintestinal adhesin and
DAEC is **not** raised. Without that rule EC958, a urosepsis reference strain,
was reported as "possible DAEC".

**UPEC four-gene rule** (Spurbeck *et al.* 2012): all four of *vat*, *fyuA*,
*chuA*, *yfcV* present.

### *Shigella* / EIEC
`blastn -perc_identity 80` against the ShigEiFinder gene set. A
cluster-specific gene counts as present at **> 50 %** reference coverage, with
the published minimum aligned lengths for `C1_gene_2` (251 bp), `C1_gene_4`
(254 bp) and `C2_gene_1` (393 bp); **ipaH** counts at **> 10 %**, because it is
multi-copy and often partial. A cluster is called when **all** of its
cluster-specific genes are present. Where more than one cluster is complete,
only the seven precedence pairs the published implementation resolves are
resolved (C1+CSB12 → CSB12, C5+C3 → C3, C5+C8 → C8, C1+C2 → C2, C1+CSD1 →
CSD1, C2+CSS → C2, CSB12+CSD1 → CSB12); anything else is reported as ambiguous,
since that is what a mixed culture looks like.

ipaH-positive with no complete cluster is reported as
`Shigella/EIEC, unclustered`.

Serotyping is **not** re-derived here: Hydra's EcOH serotyper already does it,
and ShigEiFinder's O/H antigen gene sets would duplicate it with a second,
possibly disagreeing answer.

---

## 3. Limits when working from assemblies

* **The phylogroup rests on amplicons, so it is sensitive to assembly
  fragmentation.** The partial-product rule mitigates a break through an
  amplicon but cannot rescue an amplicon lost entirely. A `U (unknown)` pattern
  most often means a cryptic clade, a non-*coli* *Escherichia*, or a
  contaminated assembly — check the species call before anything else.
* **The ST131 C1/C2 (H30R/H30Rx) split is not determinable by this method.** It
  is reported as unresolved with both candidates named.
* **The ST131 subclone label is conditional on the ST**, which this module does
  not compute. fimH41 in particular is common outside ST131: ETEC H10407
  (ST48) carries it.
* **A stx operon split across contigs cannot be subtyped**, and is reported as
  partial rather than guessed from the half that assembled. Likewise an operon
  interrupted by a phage insertion.
* **Gene presence is not gene function.** An intact ORF by alignment may carry a
  frameshift or a premature stop this module does not check; stx operon
  completeness is checked, the rest of the virulence panel is not.
* **The pathotype is a genotype.** Adherence phenotype (HEp-2 aggregative
  adherence for EAEC, diffuse adherence for DAEC), toxin production and the
  clinical picture decide the final label. DAEC in particular is reported only
  as a candidate.
* **Copy number is not reported.** Two stx2a prophages are reported as two
  operons; a single operon present in two copies in a collapsed repeat is
  reported once.
* **Reads are not used.** Nothing here needs them, but a low-level mixture that
  reads would reveal is invisible to an assembly; an ambiguous
  `shigella_eiec` result is the usual sign.

---

## 4. Database provenance and licence

Built into `<db_dir>/ecoli/` by `fetch()`; idempotent; total **14 MB** installed
(4.1 MB of source FASTA plus BLAST indices). `fetch()` also writes
`provenance.tsv` recording every row of this table.

| database | source | version | licence | redistribution |
|---|---|---|---|---|
| **fimtyper** `fimH.fsa`, 3312 alleles, 1.6 MB | `bitbucket.org/genomicepidemiology/fimtyper_db` (pinned to commit `1dd6cdd`) | commit 1dd6cdd, 2023-11-15 | **No licence file in the repository.** The repository states no terms; the FimTyper service is published for public use. | **unstated** |
| **stxtyper** `stx.prot`, 160 proteins, 45 kB | `github.com/ncbi/stxtyper` tag `v1.0.45` | StxTyper v1.0.45 | NCBI PUBLIC DOMAIN NOTICE: *"This software/database is a 'United States Government Work' under the terms of the United States Copyright Act ... thus cannot be copyrighted. This software/database is freely available to the public for use."* | **permitted** |
| **virulencefinder** `virulence_ecoli.fsa`, 7.7 MB | `bitbucket.org/genomicepidemiology/virulencefinder_db` | database 2.1.0 | **No licence file in the repository.** The repository states no terms. | **unstated** |
| **shigeifinder** `genes.fasta` + `genes.json`, 548 kB | `github.com/LanLab/ShigEiFinder` tag `v1.3.5` | ShigEiFinder v1.3.5 | **GNU General Public License v3.0** (`LICENSE` in the repository). | **permitted under GPL-3.0** |

**For the integrator's decision.** Two of the four databases state no licence at
all. They are freely downloadable and are what the published reference tools
use, but with no grant of redistribution they should be **downloaded on first
use rather than shipped inside the Windows zip** unless permission is obtained.
The ShigEiFinder data is GPL-3.0, which permits redistribution but is copyleft;
whether bundling it in a zip alongside an MIT-licensed application is acceptable
is a decision above this module. Only the NCBI StxTyper database is
unambiguously safe to ship. `fetch()` downloads all four, so a
download-on-first-use policy needs no code change.

The VirulenceFinder URL is the repository's `master`, which upstream does not
tag; the version recorded is the database's own `VERSION` file (2.1.0). Pinning
to a commit is possible and is the right change if reproducibility matters more
than currency.

---

## 5. Validation

Nine public NCBI assemblies, covering phylogroups A, B1, B2, D and E, six
pathotypes and one *Shigella*. Phylogroups are compared against **EzClermont
1.0.0** output and stx subtypes against **NCBI StxTyper v1.0.45** output, both
run locally on the same assemblies; everything else is compared against the
published characterisation of the strain, cited. Run with
`HYDRA_VALIDATION=1 pytest tests/test_typing_ecoli.py`.

| accession | strain | scheme | expected | source of the expectation | obtained |
|---|---|---|---|---|---|
| GCA_000005845.2 | K-12 MG1655 | phylogroup | A | EzClermont 1.0.0 | **A** ✓ |
| | | stx | none | StxTyper 1.0.45 | **none detected** ✓ |
| | | pathotype | none | commensal laboratory strain | **no pathotype markers** ✓ |
| | | intimin | eae absent | no LEE in K-12 | **eae absent** ✓ |
| | | fimH | — | — | fimH27 (exact) |
| GCA_000008865.2 | O157:H7 Sakai | phylogroup | E | EzClermont 1.0.0 | **E** ✓ |
| | | stx | stx1a + stx2a | StxTyper 1.0.45 | **stx1a + stx2a** ✓ |
| | | pathotype | EHEC | Hayashi 2001 | **EHEC** ✓ |
| | | intimin | gamma1 | Oswald 2000 (O157:H7 = intimin γ) | **gamma1** (100 % id, 5.3 points over γ6) ✓ |
| GCA_000732965.1 | O157:H7 EDL933 | phylogroup | E | EzClermont 1.0.0 | **E** ✓ |
| | | stx | stx1a + stx2a | StxTyper 1.0.45 | **stx1a + stx2a** ✓ |
| | | pathotype / intimin | EHEC / gamma1 | as Sakai | **EHEC / gamma1** ✓ |
| GCA_000285655.3 | EC958 (O25b:H4-ST131 H30Rx) | phylogroup | B2 | EzClermont 1.0.0 | **B2** ✓ |
| | | fimH | fimH30 | Price 2013 (H30Rx subclone), Forde 2014 | **fimH30** (100 %/100 %) ✓ |
| | | ST131_subclone | H30, clade C | Price 2013 | **H30 (clade C) if ST131** ✓ |
| | | pathotype | ExPEC | urosepsis reference strain | **ExPEC** (4/5 groups) ✓ |
| | | stx | none | StxTyper 1.0.45 | **none detected** ✓ |
| GCA_000210475.1 | ETEC H10407 | phylogroup | A | EzClermont 1.0.0 | **A** ✓ |
| | | pathotype | ETEC | Crossman 2010 | **ETEC** (LT + ST) ✓ |
| | | stx | none | StxTyper 1.0.45 | **none detected** ✓ |
| GCA_000027125.1 | EAEC 042 (O44:H18) | phylogroup | D | EzClermont 1.0.0 | **D** ✓ |
| | | pathotype | EAEC | Chaudhuri 2010 (prototypical aggR⁺ EAEC) | **EAEC (typical)**, aggR⁺ ✓ |
| | | fimH | — | — | **fimH not found** (see note) |
| GCA_000299455.1 | O104:H4 2011 outbreak | phylogroup | B1 | EzClermont 1.0.0 | **B1** ✓ |
| | | stx | stx2a | StxTyper 1.0.45 | **stx2a** ✓ |
| | | pathotype | STEC/EAEC hybrid | the defining feature of this clone | **STEC/EAEC hybrid + ExPEC** ✓ |
| GCA_000006925.2 | *S. flexneri* 2a str. 301 | shigella_eiec | *S. flexneri* cluster | Zhang 2021 (C3 = *S. flexneri* 1–5/X/Y) | **C3 (Shigella)** ✓ |
| | | pathotype | Shigella | Jin 2002 | **Shigella** ✓ |
| | | phylogroup | A | EzClermont 1.0.0 | **A** ✓ |
| GCA_000010485.1 | SE15 (ST131 clade A commensal) | phylogroup | B2 | EzClermont 1.0.0 | **B2** ✓ |
| | | fimH / ST131_subclone | fimH41 / clade A | Johnson 2014 (O16/clade A = fimH41) | **fimH41 / H41 (clade A) if ST131** ✓ |
| | | pathotype | none | commensal isolate | **no pathotype markers** ✓ |

**No disagreement with either reference tool on any genome.** 9/9 phylogroups
match EzClermont 1.0.0; 9/9 stx results match StxTyper v1.0.45, including both
double-operon genomes and the six negatives.

### Findings from validation, reported rather than tuned away

* **EAEC 042 and O104:H4 return `fimH not found`.** Not a threshold artefact: a
  separate `blastn -task blastn` search of the whole assembly with three
  divergent fimH references found **no alignment at any identity**, while the
  same probe hits MG1655 and *S. flexneri* 301. The gene is genuinely absent
  from both assemblies at the FimTyper typing fragment. The note says so and
  stops short of asserting a published cause, which could not be verified in
  the time available — worth checking before reporting it as a deletion.
* **Sakai's stx2a operon matches a stx2c reference for subunit A and a stx2a
  reference for subunit B.** That is the intended behaviour of the algorithm,
  not a conflict: stx2a/2c/2d references cross-match and the subtype is decided
  by the three residues (observed `StxA313=F StxA319=K StxB35=D` → stx2a,
  agreeing with StxTyper). The module now explains this in the note so the
  evidence does not read as self-contradictory.
* **fimH could not be validated against the FimTyper reference tool.**
  `fimtyper.pl` requires BioPerl (`Bio::SeqIO`), which is not installed and
  which the brief forbids adding to the read-only environment. fimH is
  therefore validated against the published allele for EC958 (fimH30) and SE15
  (fimH41), and against the exactness of the match (7 of 9 genomes gave a
  100 % identity / 100 % coverage hit to a single allele). The remaining
  alleles — MG1655 fimH27, Sakai/EDL933 fimH36, *S. flexneri* fimH52 — are
  exact database matches but have no published allele number this module's
  author could verify, so they are reported as obtained, not as validated.
* **A BLAST+ 2.16.0 bug forced the stx search to run single-threaded.**
  `blastx` segfaults (SIGSEGV) whenever `-max_target_seqs` is given together
  with `-num_threads > 1` and the search finds nothing. Hydra's
  `engines.blast.blast()` always passes `-max_target_seqs`, so **every
  stx-negative genome crashed**. Reproduced deterministically at every
  `max_target_seqs` from 160 to 10000, and never at `-num_threads 1`. The stx
  search therefore runs with `threads=1` (about a second against a 5 Mb
  genome). This affects any Hydra module doing a translated search against a
  small protein database on a negative sample — see the integrator notes.

### Schemes considered and left out

| left out | why |
|---|---|
| ST131 C1 / C2 (H30R / H30Rx) split | Defined on core-genome SNPs, not on a marker determinable by gene alignment. Reported as unresolved with both candidates. |
| ST131 C1-M27 subclade | Defined by the M27PP1 prophage-like region (Matsumura 2016); the region's sequence was not obtainable from a primary source in the time available, so it is not implemented rather than approximated. |
| Matsumura 2017 clade PCR (clade-specific amplicons and the *prfC* / *mgtA* SNPs) | The primer and SNP table could not be extracted from the paper within the time budget. This is the natural next step and would resolve clades A/B/C independently of fimH. |
| *Shigella*/EIEC **serotyping** | Hydra's EcOH serotyper already serotypes these genomes; adding ShigEiFinder's antigen gene sets would produce a second, possibly disagreeing serotype. Only the cluster and the *Shigella*/EIEC verdict are added. |
| ShigaPass (Yassine 2023) | An equally valid published alternative, chosen against because ShigEiFinder's cluster-specific gene sets answer the *Shigella*-vs-EIEC question directly — which is the question Hydra's species call gets wrong — while ShigaPass is primarily a serotype predictor. |
| LEE island completeness / nle effector profiling | Would need an island-boundary definition and a curated effector set; the EHEC accessory markers (*ehxA*, *katP*, *espP*, *toxB*, *nleA/B/C*) are reported instead. |
| Frameshift and premature-stop detection across the virulence panel | Only done for the stx operon, where the algorithm requires it. Extending it to the whole panel is worthwhile but was not validatable in the time. |
| Capsule (K antigen) typing beyond *kpsM* group II/III | The `kpsMII_K1`/`K5`/`K52` alleles in the VirulenceFinder set are reported as found, but no full K-antigen scheme is implemented and none was validated. |

---

## 6. Bibliography

Every reference below was verified at the time of writing against PubMed
(NCBI E-utilities `esummary`) or Crossref: the PMID or DOI resolves and the
first author, title, journal and year match. References that could not be
verified are not cited.

**Phylogroup**

1. Clermont O, Bonacorsi S, Bingen E. (2000) Rapid and simple determination of
   the *Escherichia coli* phylogenetic group. *Applied and Environmental
   Microbiology* 66:4555–4558. doi:10.1128/AEM.66.10.4555-4558.2000,
   PMID 11010916. — *the original triplex scheme the 2013 method revises.*
2. Clermont O, Christenson JK, Denamur E, Gordon DM. (2013) The Clermont
   *Escherichia coli* phylo-typing method revisited: improvement of specificity
   and detection of new phylo-groups. *Environmental Microbiology Reports*
   5:58–65. doi:10.1111/1758-2229.12019, PMID 23757131. — *the quadruplex
   primers, the group C and group E reactions and the decision table.*
3. Clermont O, Dixit OVA, Vangchhia B, Condamine B, Dion S, Bridier-Nahmias A,
   Denamur E, Gordon D. (2019) Characterization and rapid identification of
   phylogroup G in *Escherichia coli*, a lineage with high virulence and
   antibiotic resistance potential. *Environmental Microbiology* 21:3107–3117.
   doi:10.1111/1462-2920.14713, PMID 31188527. — *phylogroup G and the ybgD
   reaction.*
4. Waters NR, Abram F, Brennan F, Holmes A, Pritchard L. (2020) Easy
   phylotyping of *Escherichia coli* via the EzClermont web app and
   command-line tool. *Access Microbiology* 2:acmi000143.
   doi:10.1099/acmi.0.000143, PMID 33195978. — *the in-silico PCR
   implementation reproduced here, and the reference tool validated against.*

**fimH and ST131**

5. Roer L, Tchesnokova V, Allesoe R, Muradova M, Chattopadhyay S, Ahrenfeldt J,
   Thomsen MCF, Lund O, Hansen F, Hammerum AM, Sokurenko E, Hasman H. (2017)
   Development of a web tool for *Escherichia coli* subtyping based on fimH
   alleles. *Journal of Clinical Microbiology* 55:2538–2543.
   doi:10.1128/JCM.00737-17, PMID 28592545. — *the fimH allele nomenclature and
   database.*
6. Price LB, Johnson JR, Aziz M, Clabots C, Johnston B, Tchesnokova V,
   Nordstrom L, Billig M, Chattopadhyay S, Stegger M, Andersen PS, Pearson T,
   Riddell K, Rogers P, Scholes D, Kahl B, Keim P, Sokurenko EV. (2013) The
   epidemic of extended-spectrum-β-lactamase-producing *Escherichia coli* ST131
   is driven by a single highly pathogenic subclone, H30-Rx. *mBio* 4:e00377-13.
   doi:10.1128/mBio.00377-13, PMID 24345742. — *the H30 and H30Rx definitions,
   and the reason the C1/C2 split is not attempted here.*
7. Matsumura Y, Pitout JDD, Peirano G, DeVinney R, Noguchi T, Yamamoto M,
   Gomi R, Matsuda T, Nakano S, Nagao M, Tanaka M, Ichiyama S. (2017) Rapid
   identification of different *Escherichia coli* sequence type 131 clades.
   *Antimicrobial Agents and Chemotherapy* 61:e00179-17.
   doi:10.1128/AAC.00179-17, PMID 28584160. — *the clade A/B/C scheme and its
   association with fimH41, fimH22 and fimH30.*
8. Johnson JR, Johnston B, Clabots C, Kuskowski MA, Pendyala S, DebRoy C,
   Nowicki B, Rice J, Sokurenko EV. (2014) Rapid and specific detection,
   molecular epidemiology, and experimental virulence of the O16 subgroup
   within *Escherichia coli* sequence type 131. *Journal of Clinical
   Microbiology* 52:1358–1365. doi:10.1128/JCM.03502-13, PMID 24501035. —
   *the O16:H5 / fimH41 clade A subgroup.*

**Stx subtyping**

9. Scheutz F, Teel LD, Beutin L, Piérard D, Buvens G, Karch H, Mellmann A,
   Caprioli A, Tozzoli R, Morabito S, Strockbine NA, Melton-Celsa AR,
   Sanchez M, Persson S, O'Brien AD. (2012) Multicenter evaluation of a
   sequence-based protocol for subtyping Shiga toxins and standardizing Stx
   nomenclature. *Journal of Clinical Microbiology* 50:2951–2963.
   doi:10.1128/JCM.00860-12, PMID 22760050. — *the stx subtype nomenclature.*
10. Prasad AB, Abromaitis S, Brover V, Feldgarden M, Joensen KG, Kapsak CJ,
    Lindsey RL, Li L, Michelacci V, Schjørring S, Scheutz F, Klimke W. (2026)
    Algorithm and software to type Stx operons accurately from assembled
    genomic sequence. *Microorganisms* 14:1607. doi:10.3390/microorganisms14081607,
    PMID 42654954. — *the StxTyper algorithm reimplemented here (identity
    cut-offs, the 36 bp intergenic limit, the stx2a/2c/2d residues), the
    reference database, and the reference tool validated against.*
11. Lindsey RL, Prasad A, Feldgarden M, Gonzalez-Escalona N, Kapsak C,
    Klimke W, Melton-Celsa A, Smith P, Souvorov A, Truong J, Scheutz F. (2023)
    Identification and characterization of ten *Escherichia coli* strains
    encoding novel Shiga toxin 2 subtypes, Stx2n as well as Stx2j, Stx2m, and
    Stx2o, in the United States. *Microorganisms* 11:2561.
    doi:10.3390/microorganisms11102561, PMID 37894219. — *the stx2j–stx2o
    subtypes in the current nomenclature.*

**Pathotype, intimin and the virulence panel**

12. Nataro JP, Kaper JB. (1998) Diarrheagenic *Escherichia coli*. *Clinical
    Microbiology Reviews* 11:142–201. doi:10.1128/CMR.11.1.142, PMID 9457432. —
    *the pathotype definitions.*
13. Croxen MA, Law RJ, Scholz R, Keeney KM, Wlodarska M, Finlay BB. (2013)
    Recent advances in understanding enteric pathogenic *Escherichia coli*.
    *Clinical Microbiology Reviews* 26:822–880. doi:10.1128/CMR.00022-13,
    PMID 24092857. — *the current marker combinations, including the
    typical/atypical EPEC split.*
14. Boisen N, Østerlund MT, Joensen KG, Santiago AE, Mandomando I, Cravioto A,
    Chattopadhyay S, Gonyar LA, Jensen SG, Nataro JP, Scheutz F, Hazen TH,
    Rasko DA. (2020) Redefining enteroaggregative *Escherichia coli* (EAEC):
    genomic characterization of epidemiological EAEC strains. *PLoS Neglected
    Tropical Diseases* 14:e0008613. doi:10.1371/journal.pntd.0008613,
    PMID 32898134. — *the aggR-based EAEC definition.*
15. Joensen KG, Scheutz F, Lund O, Hasman H, Kaas RS, Nielsen EM,
    Aarestrup FM. (2014) Real-time whole-genome sequencing for routine typing,
    surveillance, and outbreak detection of verotoxigenic *Escherichia coli*.
    *Journal of Clinical Microbiology* 52:1501–1510. doi:10.1128/JCM.03617-13,
    PMID 24574290. — *the VirulenceFinder database and its 90 % / 60 %
    presence rule.*
16. Oswald E, Schmidt H, Morabito S, Karch H, Marchès O, Caprioli A. (2000)
    Typing of intimin genes in human and animal enterohemorrhagic and
    enteropathogenic *Escherichia coli*: characterization of a new intimin
    variant. *Infection and Immunity* 68:64–71. doi:10.1128/IAI.68.1.64-71.2000,
    PMID 10603369. — *the Greek-letter intimin nomenclature; and the γ intimin
    of O157:H7 used as a validation expectation.*
17. Zhang WL, Köhler B, Oswald E, Beutin L, Karch H, Morabito S, Caprioli A,
    Suerbaum S, Schmidt H. (2002) Genetic diversity of intimin genes of
    attaching and effacing *Escherichia coli* strains. *Journal of Clinical
    Microbiology* 40:4486–4492. doi:10.1128/JCM.40.12.4486-4492.2002,
    PMID 12454140. — *the extension of the intimin nomenclature.*
18. Russo TA, Johnson JR. (2000) Proposal for a new inclusive designation for
    extraintestinal pathogenic isolates of *Escherichia coli*: ExPEC. *Journal
    of Infectious Diseases* 181:1753–1754. doi:10.1086/315418, PMID 10823778. —
    *the ExPEC designation.*
19. Johnson JR, Kuskowski MA, Owens K, Gajewski A, Winokur PL. (2003)
    Phylogenetic origin and virulence genotype in relation to resistance to
    fluoroquinolones and/or extended-spectrum cephalosporins and cephamycins
    among *Escherichia coli* isolates from animals and humans. *Journal of
    Infectious Diseases* 188:759–768. doi:10.1086/377455, PMID 12934193. —
    *the five-marker ExPEC rule (≥ 2 of papAH/papC, sfa/focDE, afa/draBC,
    kpsM II, iutA).*
20. Spurbeck RR, Dinh PC Jr, Walk ST, Stapleton AE, Hooton TM, Nolan LK,
    Kim KS, Johnson JR, Mobley HLT. (2012) *Escherichia coli* isolates that
    carry vat, fyuA, chuA, and yfcV efficiently colonize the urinary tract.
    *Infection and Immunity* 80:4115–4122. doi:10.1128/IAI.00752-12,
    PMID 22966046. — *the four-gene UPEC rule.*

***Shigella* / EIEC**

21. Zhang X, Payne M, Nguyen T, Kaur S, Lan R. (2021) Cluster-specific gene
    markers enhance *Shigella* and enteroinvasive *Escherichia coli* in silico
    serotyping. *Microbial Genomics* 7:000704. doi:10.1099/mgen.0.000704,
    PMID 34889728. — *the cluster-specific gene sets, the all-genes presence
    rule and the Shigella-vs-EIEC assignment.*
22. Yassine I, Lefèvre S, Hansen EE, Ruckly C, Carle I, Lejay-Collin M,
    Fabre L, Rapoud D, Coignard B, Weill FX, Bonacorsi S,
    Mariani-Kurkdjian P. (2023) ShigaPass: an in silico tool predicting
    *Shigella* serotypes from whole-genome sequencing assemblies. *Microbial
    Genomics* 9:000961. doi:10.1099/mgen.0.000961, PMID 36951906. — *the
    alternative published method, named above as the one not chosen and why.*

**Validation strains**

23. Hayashi T, Makino K, Ohnishi M, *et al.* (2001) Complete genome sequence of
    enterohemorrhagic *Escherichia coli* O157:H7 and genomic comparison with a
    laboratory strain K-12. *DNA Research* 8:11–22. doi:10.1093/dnares/8.1.11,
    PMID 11258796. — *O157:H7 Sakai.*
24. Crossman LC, Chaudhuri RR, Beatson SA, *et al.* (2010) A commensal gone
    bad: complete genome sequence of the prototypical enterotoxigenic
    *Escherichia coli* strain H10407. *Journal of Bacteriology* 192:5822–5831.
    doi:10.1128/JB.00710-10, PMID 20802035. — *ETEC H10407.*
25. Chaudhuri RR, Sebaihia M, Hobman JL, *et al.* (2010) Complete genome
    sequence and comparative metabolic profiling of the prototypical
    enteroaggregative *Escherichia coli* strain 042. *PLoS ONE* 5:e8801.
    doi:10.1371/journal.pone.0008801, PMID 20098708. — *EAEC 042.*
26. Jin Q, Yuan Z, Xu J, *et al.* (2002) Genome sequence of *Shigella flexneri*
    2a: insights into pathogenicity through comparison with genomes of
    *Escherichia coli* K12 and O157. *Nucleic Acids Research* 30:4432–4441.
    doi:10.1093/nar/gkf566, PMID 12384590. — *S. flexneri 2a str. 301.*
27. Forde BM, Ben Zakour NL, Stanton-Cook M, Phan MD, Totsika M, Peters KM,
    Chan KG, Schembri MA, Upton M, Beatson SA. (2014) The complete genome
    sequence of *Escherichia coli* EC958: a high quality reference sequence for
    the globally disseminated multidrug resistant *E. coli* O25b:H4-ST131
    clone. *PLoS ONE* 9:e104400. doi:10.1371/journal.pone.0104400,
    PMID 25126841. — *EC958.*
28. Toh H, Oshima K, Toyoda A, *et al.* (2010) Complete genome sequence of the
    wild-type commensal *Escherichia coli* strain SE15, belonging to
    phylogenetic group B2. *Journal of Bacteriology* 192:1165–1166.
    doi:10.1128/JB.01543-09, PMID 20008064. — *SE15, and independent
    confirmation of the B2 phylogroup call.*
