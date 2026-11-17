# Typing *Haemophilus influenzae* and *Neisseria meningitidis*

Two modules, `hydra_amr.typing.organisms.hinfluenzae` and
`hydra_amr.typing.organisms.nmeningitidis`, with shared helpers in
`_hinf_nmen_common.py`. Both work from a draft assembly using only BLAST+ and
pure Python, so they run unchanged inside the frozen Windows build.

Every call states the result in the nomenclature the field uses, the evidence it
rests on, what it means, and where it is uncertain. Where an assembly cannot
settle a question, the answer is "not determined" with the reason, never a bare
negative.

---

## 1. *Haemophilus influenzae*

### 1.1 `cap` — capsule type a to f

**What is called.** The serotype letter, or `NT` for a genuinely non-typeable
isolate, or `capsule-deficient` when a remnant locus carries no
serotype-specific genes. The `lineage` field carries the state of the locus and
the `alleles` field the per-gene evidence, the *bexA* status, the IS1016 count
and the resolved copy number.

**Why it matters.** Serotype b is the vaccine-preventable serotype; a Hib case
in a vaccinated child is notifiable in most countries and triggers a
vaccine-failure investigation. Serotype a is the serotype of greatest current
concern in Arctic and indigenous North American populations and is not covered
by the Hib conjugate. Serotype f is now the commonest encapsulated type in many
post-vaccine settings. Non-typeable isolates are the commonest cause of
*H. influenzae* disease overall and are not vaccine-preventable.

**The locus.** Three regions. Region I is the export operon *bexDCBA*; region II
is the serotype-specific biosynthesis cassette (*acs1-4*, *bcs1-4*, *ccs1-4*,
*dcs1-5*, *ecs1-8*, *fcs1-3*); region III is *hcsAB*. IS1016 flanks the division
I loci.

**Method.**

1. `blastn -task blastn` of the assembly against the reference cassette, with
   hicap's own published thresholds: a gene is **complete** at ≥70% nucleotide
   identity over ≥80% of the reference, and a **fragment** at ≥80% identity over
   ≥60 reference bases. `-task blastn`, not megablast: the cross-serotype and
   truncated-gene alignments sit at 84-96% identity over short spans, and a
   107 bp *bexA* remnant is invisible to megablast's default word size.
2. A second, translated pass (`tblastn`) over the same references. This is not a
   formality: the reference *bexA* comes from the serotype f strain KR494 and
   blastn does not align it to a serotype b *bexA* at all, so strain 10810 —
   whose annotated *bexA* is intact — is reported bexA-negative on the
   nucleotide search alone. A gene the nucleotide search missed is taken from
   the translated search and labelled `by tblastn`.
3. Overlapping hits on the assembly are reduced to one per physical locus, best
   match wins. This matters because the *acs1*/*bcs1*, *acs2*/*bcs2*,
   *ccs1*/*fcs1*, *dcs1*/*ecs1*, *dcs2*/*ecs2* and *dcs5*/*ecs8* pairs are
   84-94% identical and all clear the 70% threshold; without it one gene is
   counted once per serotype.
4. Genes are clustered into candidate loci at ≤10 kb spacing, and the serotype
   is **the one with the most complete region II set**, not the best single hit.
5. IS1016 is searched separately and excluded from every completeness and
   serotype calculation.

**Deciding the state of the locus.** Taken from the definitions of Satola et al.
2007: an isolate "containing a *cap* gene but not *bexA*" is capsule-deficient,
and one "lacking both *bexA* and any of the *cap* genes" is non-typeable.

| Finding | Call | Meaning |
| --- | --- | --- |
| No region I, II or III gene above threshold | `NT` | Genuinely acapsular. IS1016 may still be present — it is in about 9.5% of non-typeable isolates, and a full-length copy remains in the *cap*-deleted laboratory strain Rd KW20 — so IS1016 is never treated as a capsule marker. |
| Region II present, *bexA* intact | serotype letter, `encapsulated` | Expected to express capsule. |
| Region II present, *bexA* absent or truncated, **locus duplicated** | serotype letter, `encapsulated (IS1016-bexA deletion variant)` | The normal arrangement of the virulent serotype a and b lineages, not a defect. The IS1016-*bexA* partial deletion removes about 1.2 kb — the 3' end of IS1016 and roughly 83% of one *bexA* copy — and is the founder mutation of the virulent Hib lineages. A complete *cap*-b duplication without *bexA* truncation has not been observed. |
| Region II present, *bexA* absent or truncated, **single copy, single contig** | serotype letter, `capsule-deficient` | The b-minus (or a-minus) state: the serotype genotype without capsule export. Serologically non-typeable, so a slide-agglutination report of "non-typeable" does not contradict the call. |
| Region II present, *bexA* truncated, **locus split across contigs** | serotype letter, `capsule expression not determined` | The copy number cannot be read, and the two states above are opposite clinical conclusions. Resolve with long reads or the IS1016-*bexA* junction PCR. |
| Region I and/or III present, no region II gene | `capsule-deficient` | Remnant locus, serotype not nameable. |

**Copy number.** Reported as a lower bound, because an assembler collapses
identical tandem repeats. A duplication is only called when **two or more**
different *cap* genes are found in more than one copy; one gene twice is as
likely to be a misassembled repeat end. *bexA* is excluded from the count, since
it is the gene deleted from one copy. Two copies is the normal serotype b state;
three to five copies occur in about a third of invasive isolates and raise
capsule production.

**Limits from an assembly.** Copy number and the extent of a deletion cannot be
resolved when the locus spans contigs. Slide-agglutination serotyping
disagreements are not necessarily errors in the sequence: serotyping has been
estimated to be wrong about 17.5% of the time.

### 1.2 `pbp3` — the *ftsI* substitution pattern

**What is called.** The rPBP3 stage and group, with the group II subgroup where
one applies, for example `stage 1 group II (IIb)` or `stage 3 group III(+)`. The
`alleles` field lists every substitution found, the subset at scheme positions,
and which essential substitution the call rests on.

**Numbering.** Counted on the PBP3 of *H. influenzae* Rd KW20 (610 aa, NCBI
protein `P45059.1`, the same sequence as GenBank U32793 HI1132 / AAC22787.1).
Typing uses the transpeptidase domain, residues 338-573, around the 379-SSN and
512-KTG motifs. `fetch()` verifies all 21 scheme positions against the
downloaded record, so a changed reference cannot silently shift the numbering.

**The scheme** (Skaare et al. 2025, Table 1; it subsumes the Ubukata,
García-Cobos and Hotomi schemes):

| Stage | Group | S385T | L389F | R517H | N526K | N526H | Y528H |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 3 | III(+) | + | + | − | + | − | − |
| 3 | III-like(+) | + | + | + | − | − | − |
| 2 | III(−) | + | − | − | + | − | − |
| 2 | III-like(−) | + | − | + | − | − | − |
| 1 | II | − | − | − | + | − | − |
| 1 | I | − | − | + | − | − | − |
| 1 | Misc1 | − | − | − | − | − | + |
| 1 | Misc2 | − | − | − | + (as N526H) | + | − |
| 0 | none | − | − | − | − | − | − |

R517H and N526K are mutually exclusive and neither is subordinate to the other;
N526K accounts for about 84% and R517H about 9% of KTG-motif resistance. S385T
is the escalator to stage 2 and L389F to stage 3. Group II is subdivided as
IIa (N526K alone), IIb (+A502V), IIc (+A502T) and IId (+I449V) — a descriptive
subdivision, not a severity ordering.

**What each stage means.** Stage 1 (low-rPBP3): ampicillin MICs at or just above
the breakpoint, so ampicillin, amoxicillin, ampicillin-sulbactam and
amoxicillin-clavulanate are reported resistant whatever the measured MIC; even
here 6% of isolates were cefotaxime-non-susceptible and 20%
meropenem-non-susceptible at meningitis breakpoints. Stage 2: cefuroxime and
oral cephalosporin activity reduced. Stage 3: in a national collection of
cefotaxime-resistant isolates 78.5% were stage 3, and 86% of
ceftriaxone-resistant and 26% of meropenem-resistant isolates sat here.

**Method.** `tblastn` of the reference PBP3 against the assembly locates *ftsI*;
the coding sequence is recovered in the hit's frame, walked back to the start
codon the reference predicts and forward to the first in-frame stop, and
translated. Substitutions are read positionally when the product is the
reference length and from the gapped alignment otherwise, so an indel does not
turn one substitution into fifty. Internal stops and contig-end truncations are
reported as caveats, not as loss of function: PBP3 is essential, so an internal
stop is almost always a sequencing or assembly artefact.

### 1.3 `beta_lactam` — the reporting category

Combines the beta-lactamase gene and the PBP3 stage:

| Beta-lactamase | PBP3 | Call |
| --- | --- | --- |
| absent | stage 0 | `BLNAS` |
| absent | stage 1 | `low-BLNAR` |
| absent | stage 2 or 3 | `BLNAR` |
| present | stage 0 | `BLPAR` |
| present | stage 1 | `BLPACR (BLPACR-I)` |
| present | stage 2 or 3 | `BLPACR (BLPACR-II)` |

*blaTEM-1* and *blaROB-1* are detected at ≥90% identity over ≥80% coverage
against the NCBI Reference Gene Catalogue reference records, so the gene names
agree with Hydra's own AMR screen. If either input is undetermined the category
is `not determined` rather than guessed.

### 1.4 `mlst` — sequence type and clonal complex

Seven loci (*adk*, *atpG*, *frdB*, *fucK*, *mdh*, *pgi*, *recA*) against the
PubMLST *H. influenzae* scheme. The profile table publishes a `clonal_complex`
column for about two thirds of sequence types, written as "ST-6 complex" — the
serotype b division I lineage. Where a sequence type has no complex, that is
said rather than left blank.

---

## 2. *Neisseria meningitidis*

### 2.1 `cps` — serogroup

**What is called.** The serogroup letter (A, B, C, E, W, X, Y, Z, H), `cnl` for
the capsule null locus, `NG` when the locus is present but cannot be expressed,
or `not determined` when no capsule gene is found at all. When the call is `NG`
the genogroup is in `alleles["genogroup"]` and the specific lesion is in the
note — that distinction is the whole point of the call.

**The locus.** Six regions in the order D-A-C-E-D'-B. Region A holds the
serogroup-specific biosynthesis genes; region C the ABC transporter *ctrA-D*;
region B the translocation genes *ctrE*/*ctrF*. The module uses the current
nomenclature, with the former names kept here for readers of older literature:

| Current | Former | Serogroup |
| --- | --- | --- |
| `csaA`-`csaD` | *mynA-D*, *sacA-D* | A |
| `cssA`, `cssB`, `cssC` | *siaA-C*, *synA-C*, *neuA-C* | shared by B, C, W, Y |
| `csb` | *siaD*B, *synD* | B |
| `csc` | *siaD*C, *synE* | C |
| `csw` | *siaD*W, *synG* | W |
| `csy` | *siaD*Y, *synF* | Y |
| `cssE` | *oatC* | C O-acetyltransferase |
| `cssF` | *oatWY* | W and Y O-acetyltransferase |
| `cseA`-`cseG` | *29eA-H* | E (formerly 29E) |
| `csxA`-`csxC` | *xcbA-C* | X |
| `cslA`-`cslC` | *lcbA-C* | L |
| `ctrE`, `ctrF` | *lipA*, *lipB* | all encapsulated |
| `ctrA`-`ctrD`, `ctrG`, `tex` | unchanged | all encapsulated |

B, C, W and Y all carry *cssA-C* and are separated **only** by the polymerase,
so the call requires at least one serogroup-specific gene; *cssA-C* alone fit all
four equally and give `NG: region A absent` rather than a guess.

**Method.** `tblastn` of the assembly against one reference per distinct product
of every capsule locus, at ≥80% amino-acid identity. Each serogroup is scored on
how much of its own region A set is present, the top scorer must include at
least one serogroup-unique gene, and the coding sequence of every essential gene
is then recovered and checked. A gene is **intact** at ≥95% coverage of its
reference product.

**Reasons an isolate with a locus is reported `NG`.** Each is named
individually, as a lesion rather than a verdict:

- a gene interrupted by an insertion of ≥200 bp — in serogroups W and Y this is
  typically IS1301, which inserts into the capsule locus and switches expression
  off. Such insertions are reversible: precise excision of IS1301 from *cssA*
  restoring capsule expression has been measured at about 4 × 10⁻⁴.
- a fragmented gene (<95% of the reference product aligned);
- an internal stop codon, the signature of the frameshift that slipped-strand
  mispairing produces — in serogroup B this is the poly-cytidine tract in *csb*
  at nucleotide 89, where seven dC residues give an expressing capsule and the
  insertion or deletion of one abolishes it;
- an essential gene not found, including a region C or B transport gene.

Faults found **only** at contig boundaries do not produce an `NG` call: that is
an assembly limitation, so the serogroup is reported with the caveat that
expression cannot be confirmed.

**`cnl`.** The capsule null locus replaces region A entirely, so there are no
biosynthesis genes to lose and the isolate can never be grouped. cnl isolates
are common in carriage and usually avirulent, because the capsule is the main
defence against complement; invasive disease caused by one is rare and
associated with complement deficiency in the patient, so a cnl isolate from a
sterile site is a reason to screen complement.

**Not implemented.** Serogroups I and K: PubMLST has no loci for *csiA-E* and
*cskA-E* under those names, so they cannot be typed and are not claimed.

### 2.2 `finetype` — PorA and FetA variable regions

**What is called.** `serogroup:P1.VR1,VR2:FetA`, with the full designation
including the sequence type in `lineage`, following the recommended nomenclature
*serogroup: PorA type: FetA type: sequence type (clonal complex)*, thus
`B:P1.7-2,4:F1-5:ST-42(cc41/44)`.

PorA VR1 and VR2 are surface loops 1 and 4 of the class 1 porin. VR families are
grouped at about 80% amino-acid identity around a prototype and named
*family*-*variant*, so `P1.7-2,4` is VR1 family 7 variant 2 with VR2 family 4
prototype (a prototype carries no hyphen). FetA VR identifiers always carry both
parts with an `F` prefix, as in `F1-5`.

**Method.** PubMLST defines these loci as short **peptides** — VR1 about 13-20
residues, VR2 about 12-20, FetA VR about 30-35 — not as whole genes. So *porA*
and *fetA* are located in the assembly with `tblastn` against a full-length
reference product, the coding sequence is recovered and translated, and the
variant is read off by finding the longest defined peptide **inside** that
translation. That avoids guessing at the region's coordinates, which move with
every indel in the loop. VR1 is sought in the N-terminal 60% of PorA and VR2 in
the C-terminal 70%.

PubMLST publishes no profile table for the fine-typing scheme (it has no primary
key), so there is no "fine type number": the string is assembled from the three
independent variable-region calls. A region that cannot be matched is rendered
`?` and the reason is in the note.

### 2.3 `bast` — Bexsero antigen sequence type

**What is called.** `BAST-<n>` where the five-peptide combination is in the
downloaded profile table, otherwise the components themselves. The components,
in the scheme's order, are **fHbp, NHBA, NadA, PorA VR1 and PorA VR2** — five,
not four: PorA VR1 is in the scheme for extra discrimination although it is not
part of the MATS assay.

Also reported: the fHbp **variant group** (1, 2 or 3) and **subfamily** (A or B)
— two nomenclatures for the same peptide, from two vaccine developers — and the
**NadA variant** under the current revised scheme (NadA-1, NadA-2/3, NadA-4/5,
NadA-6, hybrid variant). An absent gene is reported as peptide `0`, which is how
PubMLST records it; *nadA* is absent from most carriage isolates and from whole
lineages, so its absence is a result and not a failure.

**MenDeVAR.** The BAST profile table carries the Meningococcal Deduced Vaccine
Antigen Reactivity index for Bexsero and for Trumenba, and this module reports it
with its meaning spelled out rather than as a bare word:

| Category | Meaning |
| --- | --- |
| exact match | at least one antigen is an exact peptide match to a vaccine component; coverage expected |
| cross-reactive | no exact match, but at least one antigen has been shown experimentally to cross-react; coverage likely, not certain |
| none | no exact match and every antigen tested shown not to cross-react; coverage not expected |
| insufficient data | untested, or tested in too few isolates; **no prediction possible** — the commonest category, and not evidence either way |

The index is deduced from the peptide sequences, not measured on the isolate, and
says nothing about an individual person's response. Bexsero contains fHbp
peptide 1, NHBA peptide 2, NadA peptide 8 and PorA VR2 peptide 4; Trumenba
contains fHbp peptides 45 and 55.

**Method note.** The fHbp peptide definitions begin at the cysteine of the mature
lipoprotein, after signal-peptide cleavage, so the defined peptide is a substring
of the translated gene rather than equal to it. All five components are therefore
matched as the longest contained definition, which is exact for NHBA and NadA
too.

### 2.4 `mlst` — sequence type and clonal complex

Seven loci (*abcZ*, *adk*, *aroE*, *fumC*, *gdh*, *pdhC*, *pgm*). The PubMLST
profile table gives a clonal complex for about 61% of sequence types, written
`ST-32 complex`, which the module renders `cc32` for the designation (and
`cc41/44` for `ST-41/44 complex`). The clonal complex is the hyperinvasive
lineage, and it is what links apparently unrelated cases.

### 2.5 `amr_alleles` — *penA*, *gyrA*, *rpoB*

The PubMLST allele identifier for each, which is the unit meningococcal
resistance surveillance is reported in: *penA* alleles with an altered PBP2 give
reduced penicillin susceptibility, *gyrA* alleles carrying a
quinolone-resistance substitution give ciprofloxacin resistance, *rpoB* alleles
give rifampicin resistance. An exact match names the allele; anything less is a
new allele, reported as the closest one with its distance.

**This module does not convert an allele into a phenotype.** That mapping is not
published as a table this module could implement, and inventing one would be
worse than naming the allele and saying so. The result supplements phenotypic
testing; it does not replace it.

---

## 3. Database provenance and licence

| Module | Database | Source | Size | Licence | Redistribution |
| --- | --- | --- | --- | --- | --- |
| hinfluenzae | `cap` | hicap `hicap/database/*.fasta` at commit `8b7f138` | 70 kB | GPL-3.0 | permitted, but copyleft - see below |
| hinfluenzae | `ftsI` | NCBI protein `P45059.1` | 1 kB | none stated | permitted |
| hinfluenzae | `bla` | NCBI Reference Gene Catalogue `NG_050145.1`, `NG_049977.1` | 2 kB | none stated | permitted |
| hinfluenzae | `mlst` | PubMLST `pubmlst_hinfluenzae_seqdef` | ~3 MB | PubMLST terms | **restricted** - first-launch download |
| nmeningitidis | `capsule` | PubMLST capsule loci | ~10 MB | PubMLST terms | **restricted** - first-launch download |
| nmeningitidis | `peptides` | PubMLST VR and vaccine-antigen peptides, BAST profiles | ~1.6 MB | PubMLST terms | **restricted** - first-launch download |
| nmeningitidis | `mlst` | PubMLST `pubmlst_neisseria_seqdef` scheme 1 | ~7 MB | PubMLST terms | **restricted** - first-launch download |
| nmeningitidis | `amr` | PubMLST *penA*, *gyrA*, *rpoB* | ~1.7 MB | PubMLST terms | **restricted** - first-launch download |
| nmeningitidis | `antigens` | NCBI proteins `P0DH58.1`, `WP_010981016.1` | 2 kB | none stated | permitted |

### PubMLST

The terms (<https://pubmlst.org/terms-conditions>, page version dated
2026-05-13, retrieved 2026-11-17) are not a named licence but a University of
Oxford terms-of-use document. The deciding sentence:

> All data submitted on or before 31 December 2024 may continue to be
> downloaded, used, and redistributed without restriction, subject only to
> proper scientific citation and acknowledgement.

and

> The restrictions described in Sections 4 and 5 do not apply to pre-2025 data.

Sections 4 and 5 are the ones that would otherwise forbid redistribution
("not to redistribute, publicly post, publish, transfer, or aggregate
significant portions of the data into any public, semi-public, or commercial
database, knowledgebase, software system, or service") and require a commercial
licence contract for software tools that download the data.

Per the owner's licence policy a source that expressly restricts redistribution stays a first-launch download, so every PubMLST database above is marked
`restricted` and is **not** bundled in the Windows zip. Read strictly, the sentence quoted above permits redistribution of exactly what an anonymous
download can obtain, so bundling the pre-2025 snapshot would also be defensible; the conservative reading is the one applied.

**The right to redistribute depends on the download being anonymous.** The
unauthenticated API answers every request with:

> Please note that you are currently restricted to accessing data that was
> submitted on or prior to 2024-12-31. Please authenticate to access the full
> dataset.

so an anonymous download can only return pre-2025 data, which the terms place
outside the restrictions. `fetch()` therefore sends no credentials, and **an API
key must never be added**: authenticating pulls post-2024 data and with it the
no-redistribution clause.

PubMLST requires this acknowledgement wherever the data is used:

> This publication made use of the PubMLST website (<https://pubmlst.org/>)
> sited at the University of Oxford. The development of that website was funded
> by the Wellcome Trust.

together with a citation of Jolley, Bray & Maiden 2018.

The terms also reserve the right to "modify nomenclature, allele definitions, or
database structure" without notice, so the retrieval date in
`db_dir/<key>/version.json` is part of the result.

**Alternative freely redistributable sources for the same definitions.** There
is no Zenodo or figshare deposit of the PubMLST *Neisseria* allele definitions
(there is one for *Streptococcus pyogenes* cgMLST under CC-BY-4.0, with
PubMLST's own developer as a named author, which is a precedent for asking for
one). The `meningotype` tool bundles pre-2025 snapshots of the PorA, FetA, fHbp,
NHBA, NadA and BAST files, but its GPL-3.0 is the software licence and cannot
re-license PubMLST's data, so it is a convenient mirror and not an independent
right. Two derived tables **are** independently CC-BY-4.0 and may be
redistributed with attribution: the BAST definition in Brehony et al. 2016 and
the per-peptide vaccine-reactivity table (Table 2) of Rodrigues et al. 2020. The
PorA VR tables of Russell et al. 2004 are a US Government publication in the
public domain.

### hicap

GPL-3.0; the FASTA files are package data of that GPLv3 work, so redistribution
is permitted but copyleft. Bundling them inside an MIT-licensed Windows zip is a
licence-compatibility question, and `fetch()` downloads them rather than shipping
them. The clean alternative is to rebuild an equivalent database from the
public-domain GenBank records the hicap paper cites: `CP005967.1` (KR494,
*bexA-D* and *hcsA-B*, with the exact coordinates in hicap's own FASTA headers),
`Z37516.2` (*cap*-a), `AF549213.1` (*cap*-b, the duplicated locus),
`HQ651151.1` (*cap*-c), `HQ424464.1` (*cap*-d), `FM882247.1` (*cap*-e),
`AF549211.1` (*cap*-f) and `X59756.1` (IS1016-V6). Those records are annotated
inconsistently — *ccs4* is a `misc_feature`, *fcs3* has no `/gene` qualifier,
*acs2-4* appear as *orf2-4* — so a rebuild needs an alias table and must be
re-validated.

---

## 4. Validation

Run with `HYDRA_VALIDATION=1 HYDRA_DB_DIR=<dir> pytest
tests/test_typing_hinf_nmen.py`.

### 4.1 *H. influenzae* capsule — against hicap 1.0.4 run locally

hicap 1.0.4 was installed from its pinned commit with Prodigal 2.6.3 and
BLAST+ 2.16.0 and run on each assembly; its `predicted_serotype` column is the
expectation. A genome for which hicap writes no report is its way of saying "no
*cap* locus". The published serotype is given where the hicap paper reports one.

| Accession | Strain | Expected | Source of expectation | Hydra returned |
| --- | --- | --- | --- | --- |
| GCF_000210875.1 | 10810 | b, duplicated | hicap `type_b`, attributes `duplicated,truncated_genes`; published serotype b | **b**, encapsulated (IS1016-*bexA* deletion variant), ≥2 copies, IS1016 ×3 |
| GCF_900478275.1 | NCTC 13377 | b | hicap `type_b`, `truncated_genes`, not duplicated; published serotype b | **b**, capsule-deficient (single copy, *bexA* partial) |
| GCA_003497005.1 | M21328 | b | hicap `type_b`, `fragmented_locus` | **b** |
| GCA_003492045.1 | M26329 | b | hicap `type_b`, `duplicated` | **b** |
| GCA_900407865.1 | PTHi-1539 | b | hicap `type_b`, `duplicated` | **b** |
| Hi84 | Hi84 (hicap example) | b | hicap `type_b`, full gene complement | **b**, encapsulated |
| GCF_001856725.1 | NML-Hia-1 | a | hicap `type_a`, `duplicated`; published serotype a | **a**, ≥2 copies |
| Hi75 | Hi75 (hicap example) | a | hicap `type_a` | **a** |
| GCF_003351605.1 | M12125 | c | hicap `type_c`; published serotype c | **c** |
| GCF_900407385.1 | PTHi-10983 | d | hicap `type_d`; the only *cap*-d genome in the 691-genome survey | **d** |
| GCF_001975845.1 | hi467 | e | hicap `type_e`; published serotype e | **e** |
| GCF_000465255.1 | KR494 | f | hicap `type_f`; published serotype f | **f** |
| GCF_900475755.1 | NCTC 11426 | f | hicap `type_f`; published serotype f | **f** |
| GCF_000027305.1 | Rd KW20 | no *cap* locus | hicap writes no report; published as a *cap*-deleted derivative of a serotype d strain | **NT**, IS1016 present (the trap: IS1016 is not a capsule marker) |
| GCF_000012185.1 | 86-028NP | no *cap* locus | hicap writes no report | **NT** |
| GCF_000016465.1 | PittEE | no *cap* locus | hicap writes no report | **NT** |
| GCF_000165525.1 | R2866 | no *cap* locus | hicap writes no report | **NT** |
| GCF_000165575.1 | R2846 | no *cap* locus | hicap writes no report | **NT** |
| GCF_000197875.1 | F3031 | no *cap* locus | hicap writes no report | **NT** |
| GCF_008831525.1 | PittGG | no *cap* locus | hicap writes no report | **NT** |

20 of 20 serotype calls agree with hicap. Two deliberate differences in
*presentation*, both supported by the literature rather than by tuning:

- hicap counts a truncated *bexA* as "found" and reports `region_I_genes 4/4`;
  this module reports it as partial and uses it, together with the copy number,
  to separate the normal IS1016-*bexA* deletion of a duplicated locus (10810)
  from the capsule-minus state of a single-copy locus (NCTC 13377).
- hicap reports a duplication from any single gene seen twice; this module
  requires two or more genes, which is why M21328 is not called duplicated here.
  hicap also does not call it duplicated.

### 4.2 *H. influenzae* other schemes

| Accession | Strain | `pbp3` | `beta_lactam` | `mlst` |
| --- | --- | --- | --- | --- |
| GCF_000210875.1 | 10810 | stage 0, wild type | BLPAR | ST6 (ST-6 complex) |
| GCF_000012185.1 | 86-028NP | stage 0, substituted outside the scheme (A502T) | BLPAR | ST33 (ST-3 complex) |
| GCF_000016465.1 | PittEE | stage 1 group II (IIa) | low-BLNAR | ST107 (ST-107 complex) |
| GCF_001975845.1 | hi467 | stage 3 group III(+) | BLNAR | not determined |
| GCF_000027305.1 | Rd KW20 | stage 0, wild type | BLNAS | ST1621 (ST-10 complex) |
| GCF_000465255.1 | KR494 | stage 0, substituted outside the scheme | BLNAS | ST124 (ST-124 complex) |
| GCF_001856725.1 | NML-Hia-1 | stage 0, substituted outside the scheme | BLNAS | ST23 (ST-23 complex) |

Independent checks on these: ST6/ST-6 complex is the serotype b division I
lineage and is returned for both serotype b reference strains; ST-23 complex is
the serotype a lineage and is returned for NML-Hia-1; ST-124 complex is the
serotype f lineage and is returned for KR494. The *blaTEM-1* calls on 86-028NP
and R2866 were confirmed by checking the RefSeq annotation of the hit region,
which annotates "broad-spectrum class A beta-lactamase TEM-1" at exactly the
coordinates reported.

### 4.3 *N. meningitidis*

Expectations are the published serogroup of each reference strain, and for two
strains the full published fine type. No small pip-installable reference tool
could be used: `meningotype` and the CDC capsule characterisation tool are both
GitHub-only installs, and the latter has no licence statement.

| Accession | Strain | Expected | Source of expectation | Hydra returned |
| --- | --- | --- | --- | --- |
| GCF_000009105.1 | Z2491 | A | Parkhill et al. 2000 *Nature*, serogroup A in the title | **A**, `A:P1.7,13-1:F1-5:ST-4(cc4)`, BAST-639 |
| GCF_000008805.1 | MC58 | B | Tettelin et al. 2000 *Science*; BioSample SAMN02603997 `serogroup=B` | **B (capsule gene disrupted)** - *csb* truncated; see the note below |
| GCF_000191445.1 | H44/76 | B, P1.7,16, ST-32, cc32 | the reference serogroup B strain; published fine type B:P1.7,16:ST-32(cc32) | **B**, `B:P1.7,16:F3-3:ST-32(cc32)`, BAST-642 (fHbp 1, NHBA 3) |
| GCF_000191525.1 | NZ-05/33 | B, P1.7-2,4, cc41/44 | the New Zealand epidemic strain | **B**, `B:P1.7-2,4:F1-5:ST-42(cc41/44)`, BAST-223 |
| GCF_000191465.1 | M01-240149 | B | serogroup B disease isolate | **B**, `B:P1.7-2,4:F1-5`, ST-41, cc41/44 |
| GCF_000191485.1 | M01-240355 | B | serogroup B disease isolate, ST-213 | **B**, ST-213, ST-213 complex, **BAST-2255** |
| GCF_000191505.1 | M04-240196 | B | serogroup B disease isolate | **B**, ST-269, ST-269 complex |
| GCF_000152165.1 | alpha710 | B | Schoen et al. 2008; PubMLST and BioSample both give serogroup B | **B**, cc41/44 |
| GCF_000009465.1 | FAM18 | C | Bentley et al. 2007; NCBI BioSource `serotype=C` | **C**, `C:P1.5,2:F1-30:ST-11(cc11)`, BAST-3 |
| GCF_000026965.1 | 8013 | C | Rusniok et al. 2009 | **C**, ST-177, ST-18 complex |
| GCF_000083565.1 | alpha14 | cnl | Schoen et al. 2008; ENA sample `serotype=cnl` | **cnl** |

All 11 calls name the published serogroup, covering A, B, C and the capsule-null
case. MC58 names it with its lesion attached — `B (capsule gene disrupted)` —
because the gene really is truncated in the deposited sequence; see the note
below.

**How a disrupted capsule is worded.** The serogroup comes first and the lesion
is stated with it, so the call reads `B (capsule gene disrupted)` rather than a
bare `NG` with the genogroup left in another column for the reader to find. Which
lesion is detected is unchanged, and the evidence and the caveats are unchanged:
the gene-by-gene faults stay in the note and in the evidence, and the note still
says that the laboratory will report the isolate as non-groupable, that the
lesion can revert, and that the genogroup and the lesion belong together. A
capsule gene cut by a contig end is a third, separate wording —
`B (capsule gene at a contig end, expression not assessed)` — because an
assembly boundary is neither an intact capsule nor a lesion. The fine type keeps
the standard designation, whose first field is the serological result and so is
`NG` for an isolate that cannot express its locus; the reason travels with it in
the `serogroup_field` evidence, pointing at the `cps` call.

The *H. influenzae* module uses the same shape: a single-copy locus with *bexA*
absent or truncated is `b (capsule gene disrupted)`, not a bare `b`, so a
capsule-deficient isolate cannot be read as an encapsulated one. The normal
virulent Hib arrangement — a duplicated locus with one *bexA* truncated — is
still reported as plain `b`, and the test that insists on that is unchanged.

**The MC58 lesion, investigated in the sequence.** MC58 is published as
serogroup B, and this module reports a disrupted capsule because *csb*
(the polysialyltransferase, formerly *siaD*) aligns over only 45% of the
reference product. That is not a threshold artefact: RefSeq's own annotation of
the same assembly annotates the polysialyltransferase CDS at
`complement(76613..77290)`, 678 bp, against a reference *csb* of about 1488 bp.
The gene really is truncated in the deposited sequence. *csb* is also the gene
carrying the poly-cytidine tract at nucleotide 89 whose slipped-strand mispairing
switches serogroup B capsule on and off, and MC58 is the laboratory background in
which that phase variation was characterised, so a frameshifted *csb* in a stock
published as serogroup B is the expected finding rather than a surprising one.
The call is therefore left as it stands, and is now worded as the owner asked:
`B (capsule gene disrupted)` — the serogroup, which is what the published record
says, and the lesion, which is what the sequence says. The validation test is
stricter than it was: it no longer accepts a bare `NG` with a matching genogroup
hidden in another column, but requires the call itself to name the published
serogroup. The published serogroup is a serological result on a particular stock,
and an assembly can carry a reversible lesion that stock did not express, which
is exactly why both halves have to be in the call.
Two independent cross-checks landed exactly: H44/76 returns the published
`P1.7,16` with ST-32/cc32, and NZ-05/33 returns the published `P1.7-2,4` with
cc41/44 — `7-2,4` is also the PorA of NZ98/254, the Bexsero outer-membrane
vesicle strain.

**Investigated disagreements.** Two were found and fixed in the sequence rather
than tuned away:

1. With one reference allele per capsule locus, MC58 aligned to only 43% of the
   first *csb* allele and Z2491's *ctrC* fell to 77% identity, so both were
   reported `NG`. The cause was the reference, not the genomes: the capsule
   polymerases are the most polymorphic genes in the locus. Keeping one
   reference per *distinct product* of every allele fixed both without changing
   any threshold.
2. All fHbp peptides came back as new. The cause was that the PubMLST fHbp
   definitions start at the cysteine of the mature lipoprotein, so the defined
   peptide is a substring of the translated gene. Matching the longest contained
   definition fixed it.

### 4.4 *N. meningitidis* serogroups W, X, Y, E and Z

These five were implemented and unit-tested but had never been run against a
genome whose serogroup is published. They have now been. Every expectation below
comes from a paper, a reference-collection record, or both; where the record is
an NCBI BioSample the field is named, and where a paper states the method the
method is quoted.

| Accession | Strain | Expected | Source of expectation | Hydra returned |
| --- | --- | --- | --- | --- |
| GCF_000800275.1 | M7124 | W | BioSample SAMN03020254 `serogroup=W`; BioProject PRJNA260327 "serogroup W ST-11 Hajj reference strain"; Mustapha et al. 2016 *Genome Biol Evol* 8:2065 | **W**, `W:P1.5,2:F1-1:ST-11(cc11)` |
| GCF_001697425.1 | M07149 | W | BioSample SAMN04883361 `serogroup=W`; Retchless et al. 2016 *mSphere* 1:e00201-16, Methods: "the serogroup phenotype was confirmed using slide agglutination" | **W**, `W:P1.5,2:F1-1:ST-11(cc11)` |
| GCF_009664695.1 | AUSMDU00010537 | W | BioSample SAMN13191640, titled "*Neisseria meningitidis* ST11 Serogroup W", `serotype=W`; Baines et al. 2020 *Microb Genom* 6:e000471 Table 1 | **W**, `W:P1.5,2:F1-1:ST-11(cc11)` |
| GCF_015828725.1 | N.971/92 | X | Cassiolato et al. 2023 *Int Microbiol* 26:611, "six *N. meningitidis* serogroup X (MenX) isolates"; BioSample SAMN13841843 `Serogroup=X`; BioProject PRJNA601116 | **X (capsule gene disrupted)**, genogroup X, ST-15458 — see the note below |
| GCF_008330885.1 | 95-134 | Y | Stenmark et al. 2020 *Sci Rep* 10:3644 Table 1, "the single molecule real-time sequenced *N. meningitidis* serogroup Y isolates (n = 8)" | **Y**, `Y:P1.5-2,10-1:F4-1:ST-23(cc23)` |
| GCF_009664675.1 | AUSMDU00005726 | Y | BioSample SAMN13191641, titled "*Neisseria meningitidis* ST1655 Serogroup Y", `serotype=Y`; Baines et al. 2020 Table 1 | **Y**, `Y:P1.5-1,10-1:F4-1:ST-1655(cc23)` |
| GCF_016007825.1 | M47305 | E | Topaz et al. 2020 *Microb Genom* 6:mgen000477 supplementary Table S1, `Serogroup=E` (Burkina Faso 2017, ST-9945) | **E**, `E:P1.21,16-14:F3-6:ST-9945` |
| GCA_037365905.1 | M38141 | Z | Whaley et al. 2022 *Front Microbiol* 13:815044 supplementary Table 1, "slide agglutination serogrouping" = Z, genogroup Z, all capsule genes intact — the only Z among 1514 isolates | **Z**, `Z:P1.5-1,10-4:F5-1:ST-4983(cc103)` |

8 of 8, and the ten rows of section 4.3 are unchanged, so the whole set is
18 of 18. None of the five is provisional any longer.

**The defect this found: every serogroup W isolate was reported as serogroup Y.**
All three W genomes came back `Y` on the first run. The cause was not a
threshold. B, C, W and Y share the whole sialic-acid pathway and are told apart
by one gene, and *csw* (W) and *csy* (Y) are 99% identical at the amino-acid
level, so both clear any identity cut-off a translated search can set: in either
genome both genes are "present", both serogroups score four of four region A
genes, and both score one serogroup-specific gene. The ranking then fell through
to its last sort key, which was the serogroup's *name*, reverse-sorted — so "Y"
beat "W", every time, for every isolate.

The evidence that separates them had been collected and then thrown away. The
isolate's own sequence is unambiguous:

| | *csw* (W) | *csy* (Y) |
| --- | --- | --- |
| M7124, published W | **100.00%** | 99.04% |
| AUSMDU00010537, published W | **100.00%** | 99.04% |
| 95-134, published Y | 98.94% | **100.00%** |
| AUSMDU00005726, published Y | 98.94% | **100.00%** |

`rank_serogroups()` now ranks on (serogroup-specific genes present, region A
genes present, how well the isolate matches the genes found only in that
serogroup), and the serogroup's name decides nothing. The two scores go into the
call's evidence as `serogroup_discrimination`, so a reader sees why W rather than
Y. When the two are within `CPS_DISCRIMINATION_MARGIN` (0.25 points of
identity × coverage) the row says the assembly cannot separate them and asks for
serology or a serogroup-specific PCR — the margin has to be well under one point
because one point is all that separates a real W from a real Y. No threshold was
changed and no gene detection was changed; A, B, C and the capsule-null case all
return exactly what they returned before.

**The *csxA* stop in N.971/92, reported not tuned away.** The serogroup is right,
and the module additionally reports that *csxA* carries one internal stop codon,
so the call is `X (capsule gene disrupted)` rather than a plain `X`. The gene is
100% covered at 97.9% identity on a 2.19 Mb, contig-level draft, and the isolate
is a cerebrospinal-fluid isolate published as serogroup X, so the likeliest
explanation is a sequencing or assembly error in the draft rather than a genome
that cannot make capsule. It is left exactly as found: the row names the gene,
the lesion and the caveat that such lesions revert, which is the information a
reader needs in order to disbelieve it. Suppressing it on a guess about assembly
quality would be the worse error.

**What could not be validated, and what that costs.** Nothing is provisional,
but two rows are weaker than the rest and the docs should not pretend otherwise:

* **Serogroup E and serogroup Z have no closed genome.** There is no
  complete/closed assembly of either anywhere in GenBank or RefSeq, so both rows
  are contig-level. No structured serogroup field in NCBI, ENA or BV-BRC ever
  equals `E`, `29E` or `Z`; the expectations had to come from the supplementary
  tables of papers.
* **The serogroup Z row is GenBank-only (GCA).** The only RefSeq candidates for
  Z are labelled from PubMLST's `genogroup`/`capsule_group` field, which is
  derived from the capsule locus — using one to validate a capsule-locus caller
  would be circular. `GCA_037365905.1` is used instead because Whaley et al.
  state its serogroup from *slide agglutination*, which is independent of the
  sequence.
* **Serogroup L and serogroup H** are implemented and remain unvalidated: no
  public assembly with a published L or H serogroup was found. They keep the
  provisional status this section used to give all five, and their
  `SEROGROUP_MEANING` entries are not claimed to be validated.

---

## 5. Bibliography

Every reference below was verified against PubMed and Crossref at the time of
writing: the PMID and DOI resolve, and the first author, title, journal and year
match. Each entry also appears in the module's `REFERENCES` tuple with the
specific call or rule it supports.

### *H. influenzae*

1. **Watts SC, Holt KE.** hicap: in silico serotyping of the *Haemophilus
   influenzae* capsule locus. *J Clin Microbiol* 2019;57:e00190-19.
   PMID 30944197, doi:10.1128/JCM.00190-19. — *The cap reference gene set, the
   identity and coverage thresholds, the most-complete-region-II serotyping
   rule, and the reference tool the capsule call was validated against.*
2. **Kroll JS, Loynds B, Brophy LN, Moxon ER.** The *bex* locus in encapsulated
   *Haemophilus influenzae*: a chromosomal region involved in capsule
   polysaccharide export. *Mol Microbiol* 1990;4:1853-62. PMID 2082145,
   doi:10.1111/j.1365-2958.1990.tb02034.x. — *Region I as bexD, bexC, bexB,
   bexA.*
3. **Kroll JS, Hopkins I, Moxon ER.** Capsule loss in *H. influenzae* type b
   occurs by recombination-mediated disruption of a gene essential for
   polysaccharide export. *Cell* 1988;53:347-56. PMID 3284653,
   doi:10.1016/0092-8674(88)90155-9. — *bexA as the export gene whose disruption
   abolishes capsule, and the duplicated locus with bexA in the bridge region.*
4. **Kroll JS, Moxon ER, Loynds BM.** An ancestral mutation enhancing the
   fitness and increasing the virulence of *Haemophilus influenzae* type b.
   *J Infect Dis* 1993;168:172-6. PMID 8515106, doi:10.1093/infdis/168.1.172. —
   *The IS1016-bexA deletion as the founder mutation of the virulent b lineages;
   why a truncated bexA in a duplicated locus is the normal Hib genotype.*
5. **Corn PG, Anders J, Takala AK, Käyhty H, Hoiseth SK.** Genes involved in
   *Haemophilus influenzae* type b capsule expression are frequently amplified.
   *J Infect Dis* 1993;167:356-64. PMID 8421169, doi:10.1093/infdis/167.2.356. —
   *Copy number: two copies normal, three to five in about a third of invasive
   isolates, each raising capsule production.*
6. **Satola SW, Schirmer PL, Farley MM.** Complete sequence of the *cap* locus of
   *Haemophilus influenzae* serotype b and nonencapsulated b capsule-negative
   variants. *Infect Immun* 2003;71:3639-44. PMID 12761153,
   doi:10.1128/IAI.71.6.3639-3644.2003. — *Region III as hcsA and hcsB; the
   reference single-copy bexA-negative b-minus loci.*
7. **Satola SW, Collins JT, Napier R, Farley MM.** Capsule gene analysis of
   invasive *Haemophilus influenzae*: accuracy of serotyping and prevalence of
   IS1016 among nontypeable isolates. *J Clin Microbiol* 2007;45:3230-8.
   PMID 17699642, doi:10.1128/JCM.00794-07. — *The definitions the state table
   uses, and IS1016 in 9.5% of non-typeable isolates.*
8. **Kapogiannis BG, Satola S, Keyserling HL, Farley MM.** Invasive infections
   with *Haemophilus influenzae* serotype a containing an IS1016-*bexA* partial
   deletion: possible association with virulence. *Clin Infect Dis*
   2005;41:e97-e103. PMID 16267724, doi:10.1086/498028. — *The deletion and
   duplication occur in serotype a too, so the rule is not applied to b alone.*
9. **Ubukata K, Shibasaki Y, Yamamoto K, et al.** Association of amino acid
   substitutions in penicillin-binding protein 3 with beta-lactam resistance in
   beta-lactamase-negative ampicillin-resistant *Haemophilus influenzae*.
   *Antimicrob Agents Chemother* 2001;45:1693-9. PMID 11353613,
   doi:10.1128/AAC.45.6.1693-1699.2001. — *The original groups I, II and III.*
10. **Dabernat H, Delmas C, Seguy M, et al.** Diversity of beta-lactam
    resistance-conferring amino acid substitutions in penicillin-binding protein
    3 of *Haemophilus influenzae*. *Antimicrob Agents Chemother*
    2002;46:2208-18. PMID 12069976, doi:10.1128/AAC.46.7.2208-2218.2002. —
    *Subgroups IIa to IId, and the mutual exclusivity of R517H and N526K.*
11. **García-Cobos S, Campos J, Lázaro E, et al.** Ampicillin-resistant
    non-beta-lactamase-producing *Haemophilus influenzae* in Spain.
    *Antimicrob Agents Chemother* 2007;51:2564-73. PMID 17470649,
    doi:10.1128/AAC.00354-07. — *The group III-like category.*
12. **Hasegawa K, Chiba N, Kobayashi R, et al.** Rapidly increasing prevalence of
    beta-lactamase-nonproducing, ampicillin-resistant *Haemophilus influenzae*
    type b in patients with meningitis. *Antimicrob Agents Chemother*
    2004;48:1509-14. PMID 15105098, doi:10.1128/AAC.48.5.1509-1514.2004. —
    *low-BLNAR, BLNAR, BLPACR-I and BLPACR-II, and the role of S385T.*
13. **Skaare D, Anthonisen IL, Caugant DA, et al.** Multilocus sequence typing
    and *ftsI* sequencing: a powerful tool for surveillance of
    penicillin-binding protein 3-mediated beta-lactam resistance in nontypeable
    *Haemophilus influenzae*. *BMC Microbiol* 2014;14:131. PMID 24884375,
    doi:10.1186/1471-2180-14-131. — *rPBP3 versus sPBP3, the high/low grouping,
    and the transpeptidase-domain numbering on Rd KW20.*
14. **Skaare D, Anthonisen IL, Zecic N, et al.** Genomic surveillance of
    cefotaxime-resistant *Haemophilus influenzae* in Norway. *Front Microbiol*
    2025;16:1601390. PMID 40800114, doi:10.3389/fmicb.2025.1601390. — *The rPBP3
    stage system implemented here, and the cephalosporin and carbapenem
    consequences of each stage.*
15. **Osaki Y, Sanbongi Y, Ishikawa M, et al.** Genetic approach to study the
    relationship between penicillin-binding protein 3 mutations and *Haemophilus
    influenzae* beta-lactam resistance. *Antimicrob Agents Chemother*
    2005;49:2834-9. PMID 15980357, doi:10.1128/AAC.49.7.2834-2839.2005. —
    *Causal, not merely associative, role of the individual substitutions.*
16. **Meats E, Feil EJ, Stringer S, et al.** Characterization of encapsulated and
    noncapsulated *Haemophilus influenzae* and determination of phylogenetic
    relationships by multilocus sequence typing. *J Clin Microbiol*
    2003;41:1623-36. PMID 12682154, doi:10.1128/JCM.41.4.1623-1636.2003. — *The
    seven-locus MLST scheme.*
17. **Feldgarden M, Brover V, Gonzalez-Escalona N, et al.** AMRFinderPlus and the
    Reference Gene Catalog facilitate examination of the genomic links among
    antimicrobial resistance, stress response, and virulence. *Sci Rep*
    2021;11:12728. PMID 34135355, doi:10.1038/s41598-021-91456-0. — *The
    reference records and nomenclature for blaTEM-1 and blaROB-1.*

### *N. meningitidis*

18. **Harrison OB, Claus H, Jiang Y, et al.** Description and nomenclature of
    *Neisseria meningitidis* capsule locus. *Emerg Infect Dis* 2013;19:566-73.
    PMID 23628376, doi:10.3201/eid1904.111799. — *Region structure D-A-C-E-D'-B,
    the current gene names with their former names, and the rule that the
    polymerase alone separates B, C, W and Y.*
19. **Marjuki H, Topaz N, Rodriguez-Rivera LD, et al.** Whole-genome sequencing
    for characterization of capsule locus and prediction of serogroup of
    invasive meningococcal isolates. *J Clin Microbiol* 2019;57:e01609-18.
    PMID 30567750, doi:10.1128/JCM.01609-18. — *The serogroup prediction rules
    implemented here: essential gene sets, the requirement for a
    serogroup-specific gene, and partial genes, internal stops and
    insertion-sequence disruptions as reasons to report a genogroup with a
    non-groupable phenotype.*
20. **Hammerschmidt S, Hilse R, van Putten JP, et al.** Modulation of cell
    surface sialic acid expression in *Neisseria meningitidis* via a
    transposable genetic element. *EMBO J* 1996;15:192-8. PMID 8598202,
    doi:10.1002/j.1460-2075.1996.tb00428.x. — *Reversible inactivation of cssA
    by IS1301 insertion and precise excision.*
21. **Hammerschmidt S, Müller A, Sillmann H, et al.** Capsule phase variation in
    *Neisseria meningitidis* serogroup B by slipped-strand mispairing in the
    polysialyltransferase gene (*siaD*). *Mol Microbiol* 1996;20:1211-20.
    PMID 8809773, doi:10.1111/j.1365-2958.1996.tb02641.x. — *The poly-cytidine
    tract in csb at nucleotide 89, seven dC residues meaning an expressing
    capsule.*
22. **Claus H, Maiden MCJ, Maag R, Frosch M, Vogel U.** Many carried
    meningococci lack the genes required for capsule synthesis and transport.
    *Microbiology* 2002;148:1813-19. PMID 12055301,
    doi:10.1099/00221287-148-6-1813. — *The capsule null locus and its
    prevalence in carriage.*
23. **Jolley KA, Brehony C, Maiden MCJ.** Molecular typing of meningococci:
    recommendations for target choice and nomenclature. *FEMS Microbiol Rev*
    2007;31:89-96. PMID 17168996, doi:10.1111/j.1574-6976.2006.00057.x. — *The
    designation serogroup: PorA type: FetA type: sequence type (clonal
    complex).*
24. **Russell JE, Jolley KA, Feavers IM, Maiden MCJ, Suker J.** PorA variable
    regions of *Neisseria meningitidis*. *Emerg Infect Dis* 2004;10:674-8.
    PMID 15200858, doi:10.3201/eid1004.030247. — *The VR nomenclature: families
    grouped at about 80% identity around a prototype, named family-variant,
    rendered with the P1. prefix.*
25. **Thompson EAL, Feavers IM, Maiden MCJ.** Antigenic diversity of
    meningococcal enterobactin receptor FetA, a vaccine component.
    *Microbiology* 2003;149:1849-58. PMID 12855736, doi:10.1099/mic.0.26131-0. —
    *The FetA variable region and its F<family>-<variant> nomenclature.*
26. **Brehony C, Rodrigues CMC, Borrow R, et al.** Distribution of Bexsero
    antigen sequence types (BASTs) in invasive meningococcal disease isolates:
    implications for immunisation. *Vaccine* 2016;34:4690-7. PMID 27521232,
    doi:10.1016/j.vaccine.2016.08.015. — *BAST and its five components in order,
    with PorA VR1 included for discrimination although not part of MATS.*
27. **Rodrigues CMC, Jolley KA, Smith A, et al.** Meningococcal Deduced Vaccine
    Antigen Reactivity (MenDeVAR) Index. *J Clin Microbiol* 2020;59:e02161-20.
    PMID 33055180, doi:10.1128/JCM.02161-20. — *The four MenDeVAR categories and
    the vaccine antigen peptides quoted.*
28. **Masignani V, Comanducci M, Giuliani MM, et al.** Vaccination against
    *Neisseria meningitidis* using three variants of the lipoprotein GNA1870.
    *J Exp Med* 2003;197:789-99. PMID 12642606, doi:10.1084/jem.20021911. —
    *The fHbp variant groups 1, 2 and 3.*
29. **Fletcher LD, Bernfield L, Barniak V, et al.** Vaccine potential of the
    *Neisseria meningitidis* 2086 lipoprotein. *Infect Immun* 2004;72:2088-100.
    PMID 15039331, doi:10.1128/IAI.72.4.2088-2100.2004. — *The fHbp subfamilies
    A and B.*
30. **Serruto D, Spadafina T, Ciucchi L, et al.** *Neisseria meningitidis*
    GNA2132, a heparin-binding protein that induces protective immunity in
    humans. *Proc Natl Acad Sci USA* 2010;107:3770-5. PMID 20133713,
    doi:10.1073/pnas.0915162107. — *NHBA as a vaccine antigen and its peptide
    nomenclature.*
31. **Bambini S, De Chiara M, Muzzi A, et al.** *Neisseria* adhesin A variation
    and revised nomenclature scheme. *Clin Vaccine Immunol* 2014;21:966-71.
    PMID 24807056, doi:10.1128/CVI.00825-13. — *The current NadA variant names
    NadA-1, NadA-2/3, NadA-4/5, NadA-6 and hybrid variant.*
32. **Maiden MCJ, Bygraves JA, Feil E, et al.** Multilocus sequence typing: a
    portable approach to the identification of clones within populations of
    pathogenic microorganisms. *Proc Natl Acad Sci USA* 1998;95:3140-5.
    PMID 9501229, doi:10.1073/pnas.95.6.3140. — *The seven-locus Neisseria MLST
    scheme, defined on the meningococcus.*

### Both

33. **Jolley KA, Bray JE, Maiden MCJ.** Open-access bacterial population
    genomics: BIGSdb software, the PubMLST.org website and their applications.
    *Wellcome Open Res* 2018;3:124. PMID 30345391,
    doi:10.12688/wellcomeopenres.14826.1. — *The database every allele, peptide,
    profile and clonal complex is downloaded from. PubMLST asks that this be
    cited.*

### Verified but not implemented

The following were read and are cited here because they bear on calls this
module deliberately does **not** make, so that the omission is on the record:

- **Davis GS, Sandstedt SA, Patel M, Marrs CF, Gilsdorf JR.** Use of *bexB* to
  detect the capsule locus in *Haemophilus influenzae*. *J Clin Microbiol*
  2011;49:2594-601. PMID 21525217, doi:10.1128/JCM.00089-11. — *bexB is the best
  single marker gene for the capsule locus and bexA the worst, which is why
  bexA is used here only for the expression state and never to name a serotype.*
- **Falla TJ, Crook DW, Brophy LN, et al.** PCR for capsular typing of
  *Haemophilus influenzae*. *J Clin Microbiol* 1994;32:2382-6. PMID 7814470,
  doi:10.1128/jcm.32.10.2382-2386.1994. — *The PCR capsular typing this module's
  results will be compared against in laboratories that still run it.*
