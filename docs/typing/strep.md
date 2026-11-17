# Streptococcal typing: *S. pyogenes* (group A) and *S. agalactiae* (group B)

Two modules, `hydra_amr.typing.organisms.spyogenes` and
`hydra_amr.typing.organisms.sagalactiae`, with their shared machinery in
`strep_common`. Assemblies in; blastn, blastx and makeblastdb out; nothing else.

Every call states four things, because a gene name on its own is not a result:
**what was found** in the nomenclature the field uses, **the evidence** (contig,
coordinates, identity, coverage, the reference accession), **what it means**
clinically, and **what is not known**. Where the evidence is missing the answer
is "not determined" and the reason, never a bare negative.

---

## 1. *Streptococcus pyogenes*

### 1.1 emm type and subtype — scheme `emm`

**What is reported.** The emm type and subtype (`emm1.0`), the emm cluster
(`A-C3`), the emm pattern (`A-C`) and its tissue-tropism grouping, the other
emm-region loci found, and the evidence.

**What it means.** The emm type is the primary epidemiological label for group A
streptococcus and the unit in which outbreaks, rheumatic-fever risk and the
30-valent M-protein vaccine's coverage are all counted. The subtype resolves
within a type and is what a reference laboratory reports. The **emm cluster**
groups types whose M proteins share functional properties - the same host-ligand
binding, and cross-reactive immunity - so it says more about biology than the
type number does. The **emm pattern** separates throat specialists (pattern A-C)
from skin specialists (pattern D) from generalists (pattern E), which is the
tropism framework behind the epidemiology of pharyngitis against impetigo.

**Method.** blastn of the assembly against the CDC trimmed emm database, with
`-dust no`, `-perc_identity 95` and `-culling_limit 5`. Each HSP is then held to
the CDC subtype-calling rule as the reference tool implements it:

| Filter | Threshold |
| --- | --- |
| mismatches | ≤ 4 |
| gap openings | ≤ 2 |
| \|subject length − alignment length\| | ≤ 5 bases |
| identity | ≥ 95 % |

Surviving hits are scored `identity − (gap openings + |alignment length −
subject length|)`, so a perfect, full-length, gap-free match scores exactly
100 and is reported without a qualifier; anything else gets a `*` and a note
saying by how much it differs and that it may be a new subtype.

**Choosing between loci.** *emm* sits in a tandem array with *mrp* and *enn*,
and both of those match the emm database. Hits are therefore grouped by position
(single linkage, 500 bp reach) and each group is one locus with one answer. The
group that wins is not simply the highest-scoring one: a set of types that the
reference laboratories classify as **emm-like** is carried at the *mrp* or *enn*
position, so an emm-like locus never outranks a true emm locus. This matters.
MGAS15252 (emm59) has three emm-region loci matching the database, two of them
perfectly - `emm174.1` at 1,564,601 and `emm59.0` at 1,566,000 - and ranking by
score alone reports `emm174.1` and gets the type wrong.

**The emm cluster table.** The CDC curators annotate the cluster inside the
free-text description of many subtype records (`emm1.10 (emm cluster A-C3)
subtype definitive sequence`). `fetch()` collects those per type, so the table
comes from the primary source rather than from a third-party copy. It also
inherits the occasional mistyped label, so a type's cluster is accepted only
when at least two thirds of its annotated subtypes agree; the counts are written
to `emm_clusters.tsv` so the decision can be audited. One case is not a
disagreement at all: a single-protein cluster (`M5`) and `clade Y` are nested,
because clade Y is the group the single-protein clusters were collected into, so
the specific label wins. The current CDC file annotates 171 of 407 emm types.

**The emm pattern.** Derived from the cluster series: `A-C*` → pattern A-C
(throat), `D*` → pattern D (skin), `E*` → pattern E (generalist). A
single-protein cluster or clade Y does not map onto a pattern, and the
intermediate `A-D4` label is not read as pattern D; both are reported as not
determined with the reason.

**Limits from an assembly.** The hypervariable region is 180 bases at the 5' end
of a gene whose 3' end is repetitive, which is exactly where assemblers break. A
call whose alignment reaches within 300 bp of a contig end carries that warning.
A region broken *inside* the hypervariable region produces no passing hit at
all, and the module then reports `not determined` with the nearest candidate
sequences, because a fragmented emm region and a genuinely new emm subtype look
identical from an assembly and must not be told apart by guessing.

### 1.2 Phage-encoded virulence repertoire — schemes `emm_virulome`, `spy_prophage`

**What is reported.** Each superantigen and DNase as present or not, split into
prophage-encoded and chromosomal, each with its contig, coordinates, amino-acid
identity, reference coverage and the reference accession; then the prophage
regions the assembly resolves, with the cargo that falls inside each, and an
explicit statement of what it could not resolve.

| Gene | Location | Note |
| --- | --- | --- |
| speA | prophage | the scarlet-fever toxin; strongly associated with invasive disease |
| speC, speH, speI, speK, speL, speM | prophage | prophage-encoded superantigens; speH and speI usually travel together |
| speG | chromosome | present in almost every isolate, so its *absence* is the finding |
| speJ | chromosome | characteristic of the emm1 lineage |
| smeZ | chromosome | the most potent group A streptococcal superantigen |
| speB | chromosome | a cysteine protease, **not** a superantigen; reported separately so it is not counted among them |
| sdaB (DNase B / SpeF) | chromosome | |
| mf2 (Spd3), sda3, sda1, sdaD | prophage | sda1 degrades neutrophil extracellular traps and is a determinant of the invasive M1T1 phenotype |

**Method.** One blastx pass of the assembly against a protein panel built from
RefSeq annotation (§4). Presence is measured at the protein level because the
translated search tolerates the synonymous divergence between lineages and
yields the residue alignment a lesion call needs. `location` is the *published*
genomic context of each gene, not something measured here; it tells the reader
whether a positive is expected to travel on a prophage.

**Allele variants.** Where a gene differs from its reference the module reports
*how many residues* differ and the identity, and says the variant is not named.
The published allele series for *speA* and *smeZ* are real but are not
distributed as a versioned, downloadable reference set, so an allele number here
would be an assertion the evidence does not carry. The residue difference is the
evidence the number would have rested on.

**Prophage regions.** A prophage is described where at least three different
phage module classes sit within 20 kb of one another on one contig. The classes
are lysogeny (integrase, repressor), packaging (terminase), head (portal, major
and minor capsid, head morphogenesis), tail (tape measure, major tail, tail
spike, tail protein) and lysis (holin, HylP); group A streptococcal prophages
are 33–48 kb, so 20 kb joins the modules of one element without bridging two.
Every copy of every marker is used, not the best copy - each prophage has its
own terminase and portal, so keeping one hit per gene would describe at most one
prophage however many the isolate carries.

**Boundaries are never called**, and the module says so on every call: an
assembly rarely resolves a prophage's attachment sites, so a region is reported
as the span its marker genes cover and nothing more. This is not a published
prophage-boundary caller and is not presented as one.

**Where fragmentation wins.** A prophage-encoded gene on a contig that carries
too few phage module genes is reported under `unplaced_phage_cargo`, with its
contig and that contig's length, and the note says plainly that this is contig
fragmentation and **not** evidence that the gene is chromosomal — draft
assemblies break inside prophages, which are repetitive, far more often than
elsewhere. An isolate with no phage module gene at all is flagged as implausible
rather than reported as prophage-free.

### 1.3 Regulatory lesions — scheme `spy_regulators`

**What is reported.** *covR*, *covS*, *rocA*, *mga*, *ropB* and the *nga/ifs/slo*
operon, each as **intact** or as a named lesion: premature stop, frameshift,
truncation, or a fragment the assembly broke. The summary names only the covRS
axis (*covR*, *covS*, *rocA*), because that is the axis with a clinical reading.

**What it means.** *covS* is present in every group A streptococcus, so
"present" is not a result. A loss-of-function lesion in *covR* or *covS*
de-represses capsule, SpyCEP, streptokinase and the *nga/slo* operon and
abolishes SpeB activity: it is the switch selected *during* invasive infection,
and such isolates are typically recovered from sterile sites. It is not a
transmissible trait and it does not change antimicrobial choice. *rocA* is a
positive regulator of *covR*, so a *rocA* truncation raises CovR-repressed
virulence factors; it is naturally truncated in some serotypes rather than
acquired.

**How a lesion is recognised.** From the blastx alignment:

| Signal | Call |
| --- | --- |
| `*` in the translated query, inside the reference | premature stop |
| the gene aligns in more than one reading frame on one strand with no stop between the pieces | frameshift |
| reference coverage < 95 %, or (for a regulator) either terminus unaligned by > 5 residues | truncation |
| any of the above with the alignment within 300 bp of a contig end | fragment — the assembly is blamed, and the lesion is not determined |
| internal gap ≤ 10 residues | intact, with the in-frame indel reported |

Two of those rules were added because validation found them missing. The
MGAS5005 *covS* allele aligns from reference residue 23 to 500 — 95.6 % of the
protein — and a coverage rule alone called it intact; the missing N-terminus is
the published single-base deletion that RefSeq also annotates as a pseudogene.
The MGAS315 (emm3) *rocA* allele covers 100 % of the reference with no stop
codon, in two reading frames: a frameshift that coverage cannot see.

**M1UK is not called.** M1UK is defined by a set of SNPs relative to M1global
and is assigned phylogenetically. The module reports the *nga/ifs/slo* genes and
their state and states, in the result itself, that the lineage call needs the
SNP set and is not implemented. Reporting the operon's presence as an M1UK call
would be wrong: every group A streptococcus has the operon.

### 1.4 Adhesins — scheme `spy_adhesins`

Fibronectin-binding (*prtF1*/SfbI, *fbaB*/PrtF2/SfbII, *fbp54*), collagen-binding
(*scl1*, *scl2*, *cpa*), the FCT-region pilus (*cpa*, *fctA*, *fctB*, *srtB*) and
the Mga-regulon surface proteins (*emm*, *enn*, *scpA*, *sic*, *grab*), plus
*ska*, *scpC*/SpyCEP, *ideS*/Mac and *hasA*.

*scl1*, *scl2*, *sic* and *emm* carry variable-length repeats, so reference
coverage below 100 % is normal and is not reported as a lesion. The **FCT region
type (FCT-1 to FCT-9) and the T antigen type are not called**: both need the
pilus backbone allele matched against a published allele set, and no versioned,
downloadable set exists.

---

## 2. *Streptococcus agalactiae*

### 2.1 Capsular genotype — scheme `cps`

**What is reported.** The capsular genotype (Ia, Ib, II–IX, or non-typeable),
the serotype III reference variant where one fits, the coverage and identity it
rests on, the next-best serotype, the GBS-SBG uncertainty flags, and whether the
serotype is covered by the hexavalent capsular conjugate.

**What it means.** Capsular serotype is the basis of GBS epidemiology and of
every capsular conjugate vaccine in development; serotype III dominates neonatal
invasive disease and serotype IV is rising. The hexavalent candidate covers Ia,
Ib, II, III, IV and V.

**Method.** blastn against the GBS-SBG reference set, applying the published
rules of that tool:

| Step | Rule |
| --- | --- |
| per-HSP filter | identity > 90 % **and** alignment length > 100 bases (both strict) |
| accumulate | sum bit scores, merge HSP intervals in *reference* coordinates across contigs, take the length-weighted mean identity over the merged footprint |
| positive call | merged coverage ≥ 90 % **and** identity ≥ 90 % |
| rank | highest summed bit score |
| non-typeable | no serotype satisfies both thresholds |

Uncertainty flags, with the reference tool's own thresholds: coverage < 99 %,
identity < 99 %, more than one HSP needed, more than one contig needed, or a
next-best serotype reaching > 95 % of the winning bit score.

**Two deliberate departures from the reference tool, both validated.**

1. Accumulation is **per database record, then collapsed to a serotype**, not
   per serotype directly. The four serotype III references are 99.3–99.6 %
   identical to one another; pooling their bit scores under one `III` key sums
   four scores into one, serotype III then outranks every other serotype in
   every isolate, and a serotype Ia genome fails serotype III's coverage test
   and comes out non-typeable. This was a real bug here, caught by validation.
2. The close-runner-up flag compares against the best reference of a
   **different** serotype. Comparing against another serotype III variant raises
   the flag on every serotype III isolate — the commonest clinical serotype — and
   the flag stops carrying information. The same-serotype runner-up is still
   reported, in its own field.

Serotypes VII and IX are 91.2 % identical over their full length and both clear
the 90/90 threshold against each other, so the bit-score ranking, not the
threshold, is what separates them. "Passes the threshold" is never the call.

**This is a genotype, not a serotype, and the module says so on every call.**
An isolate that is non-typeable by serology usually carries an intact, typeable
*cps* locus: 89 of 96 serologically non-typeable strains in the population study
carry an inactivating mutation in *cps*, most often a premature stop in *cpsE*,
and a point mutation is invisible to a coverage-and-identity rule. Three
published serologically non-typeable strains in the validation panel are called
V, Ia and II here at full coverage and > 99.7 % identity — correctly, as
genotypes. The discrepancy is in capsule *expression*.

**Genotypic non-typeability** is reported with the best reference's coverage and
identity and a note pointing at the *cpsA* and *neuA* calls: if those are present
and the serotype is not, the locus is there but its serotype-determining region
is not resolved. Assembly-based calling loses about 5 % of isolates to
non-typeable against read mapping, all from incomplete *cps* assembly, so
non-typeable is always read as "non-typeable **or** incomplete *cps* assembly".

### 2.2 Surface proteins — scheme `sag_surface`

*bca* (alpha-C protein) and *bac* (beta-C protein, IgA-binding) are annotated
separately in RefSeq and are called by name. Both carry large tandem repeats
whose number varies between strains, so reference coverage well below 100 % is
normal and is not a lesion.

The alpha-like protein family member itself is named in scheme `sag_alp`
(§2.2a); `sag_surface` quotes that call in its `alp_family_member` field rather
than making a second, weaker attempt at it from the protein panel.

### 2.2a Alpha-like protein family member — scheme `sag_alp`

**What is called.** Which member of the alpha-like protein family the isolate
carries: `bca (alpha-C protein)`, `alp1 (epsilon)`, `alp2/alp3`, `rib`, or
`not determined` with the candidates and the reason.

**Why it matters clinically.** The AlpN protein-vaccine candidate is built from
the N-terminal domains of the alpha-C protein and Rib (Pawlowski 2022), so the
member an isolate carries is a vaccine-coverage statement: *bca* and *rib* are
homotypic targets, while Alp1, Alp2 and Alp3 share those domains and give
heterotypic cross-reactivity without being a homotypic match. The four
categories reported here — ALPHA, ALP1, ALP23, RIB — are also the four used in
US national GBS surveillance (McGee 2021), so a call is directly comparable to
published distributions.

**Why not from the whole gene.** The members differ mostly by how many copies
of their tandem repeat they carry. Aligning the published *bca*, *rib*, *alp1*,
*alp2*, *alp3* and *alp4* genes against a finished GBS genome returns between
100 % and 2 159 % "coverage" of each, because every repeat matches every other
repeat, so a whole-gene ranking sorts the references by repeat count, not by
member. This was measured on all twelve validation genomes before the method
was chosen.

**The method.** The CDC *Streptococcus* Laboratory's GBS surface-gene target
set (§4.5) carries four *alp* targets — ALPHA (516 bp), ALP1 (516 bp), ALP23
(620 bp) and RIB (522 bp) — each a short window in the member-specific
N-terminal region. The assembly is searched with `blastn` (task `blastn`, dust
off) and, per target, the union of its alignments gives its coverage and the
best alignment its identity. A target passes when **at least 99 % of it is
covered at no more than 8 % divergence** — the thresholds the CDC pipeline runs
these targets at (`srst2 --min_coverage 99.0 --max_divergence 8`, Metcalf 2017;
read depth, the pipeline's third criterion, has no meaning for an assembly).
The passing target with the highest identity names the member, **provided it
beats the next passing target by 5 percentage points of identity**. On the
validation panel exactly one target ever passed and the observed margin to the
next target was 21–32 points, so that belt is wide, not tuned; it exists so
that a divergent or recombinant *alp* gene comes out as `not determined` with
its candidates. A wrong member is worse than none.

Each call also reports `closest_published_allele`: the deposited allele
(§4.6) the region most resembles, from the single best alignment per allele.
This is a second opinion and the only route by which *alp4* can be named; it is
never the call.

**What is still not named, and why.**

* **alp2 against alp3.** They are mosaics of each other (Lachenauer 2000) and
  inside the CDC window they differ by **three bases in 620** — measured here
  against AF208158.1 (*alp2*) and AF245663.1 (*alp3*). Three bases is ordinary
  allelic variation, so the window matching one of them exactly does not
  establish which gene it is. The CDC target set calls the pair ALP23 for this
  reason and so does this module; the closer of the two appears in
  `closest_published_allele` as evidence only. Separating them would need the
  full gene, which the repeats make unreliable.
* **alp4.** It is not in the CDC target set. Its deposited sequence
  (AJ488912.1) reaches at most 71.6 % identity to any CDC target, so an *alp4*
  isolate fails the threshold and is reported as an unnamed family member with
  *alp4* among the candidates, rather than rounded to the nearest named target.

### 2.3 Pilus islands — scheme `sag_pilus`

**What is called.** Which pilus islands the isolate carries, by name: `PI-1`,
`PI-2a`, `PI-2b`, or combinations of PI-1 with one of the other two.

**The method.** The same CDC target set carries nine pilus targets, each a
60 bp window unique to one island: PI1 and PI1B (two alleles of PI-1), PI2A1 to
PI2A4 (four alleles of PI-2a) and PI2B and PI2B2 (two of PI-2b). They are
screened and thresholded exactly as the *alp* targets are. The RefSeq protein
panel's operon genes — *pilA*, *pilB*, *pilC*, *srtC4* for PI-2a, and the
SpaH/EbpB backbone with *srtC1* and *srtC2* for the other locus — are reported
beside the call as corroborating evidence, because a 60 bp window is a presence
call and says nothing about whether the whole operon is intact.

**The rule the call is held to.** PI-2a and PI-2b are two variants of one
island at one chromosomal locus (Rosini 2006), so an isolate carries one of
them and never both. If targets for both pass, the module reports **`not
determined`** and says the assembly is either a mixture of two isolates or
carries a mis-joined contig — it does not pick one. Among 898 human GBS
isolates PI-1 was found in 70 %, PI-2a in 79 % and PI-2b in 21 %, every isolate
carried at least one island, and PI-2a and PI-2b summed to 100 % (Martins
2013); those are the frequencies a profile is read against, and a genome with
no island at all raises an uncertainty flag.

**Clinical meaning.** PI-2b is strongly associated with the hypervirulent
neonatal clone CC17 and is read together with the `sag_markers` CC17 call;
PI-2a's backbone PilB and tip adhesin PilA are pilus-vaccine candidates
(Margarit 2009); PI-1 sits on a mobile element and is gained and lost within
lineages, so its absence is not a lineage statement.

### 2.4 Hypervirulence and neonatal-disease markers — scheme `sag_markers`

| Marker | What it means |
| --- | --- |
| **hvgA** | the CC17 hypervirulence and meningeal-tropism marker |
| **srr2** | replaces *srr1* in CC17; an independent CC17 marker |
| srr1 | the serine-rich repeat adhesin of the other lineages |
| lmb, scpB, cfb, sip, fbsA, fbsB | the core surface and virulence repertoire; *sip* is a serotype-independent vaccine candidate |
| cpsA, neuA | the *cps* locus itself, so a non-typeable call can be interpreted |
| pbp2x | see §2.5 |

**hvgA needs care, and gets it.** HvgA is the C variant of the *gbs2018* locus;
the A variant at the same locus is BibA, which every group B streptococcus
carries. The two share their N- and C-terminal domains and are only about 38 %
identical overall, so a full-length query cross-hits *bibA* in every genome
unless coverage is enforced. The module requires ≥ 75 % amino-acid identity over
≥ 80 % of the 522-residue reference. That threshold sits in the middle of a wide,
measured gap: on the validation panel the two hvgA-positive CC17 genomes match at
85.4 % and 100.0 % identity over the full reference, while the bibA-carrying
genomes reach only 57.9 % and 60.8 % over 82–85 % of it — the shared termini and
none of the divergent core. Across 20 genomes, *hvgA* and *srr2* were positive in
exactly the two CC17 strains and negative in the other eighteen, with no overlap.

**The clonal complex is not assigned here.** It follows from the MLST sequence
type, which Hydra reports separately. PubMLST does publish a clonal-complex
field for this species, which would let the integrator join ST to CC; see §6.

### 2.5 PBP2X substitutions — scheme `sag_pbp2x`

**What is reported.** Named substitutions in the PBP2X transpeptidase domain,
numbered against a pinned reference, split into the two with direct MIC evidence
(**V405A**, **Q557E**) and the others reported alongside them in the literature;
plus any change in the three active-site motifs, and any catalogued position the
alignment did not cover.

**What it means.** PBP2X is the primary penicillin target in this species.
V405A sits immediately C-terminal to the SSN substrate-binding motif (402–404)
and Q557E three residues C-terminal to the KSG motif (552–554); both remodel the
active site. V405A alone accompanies a penicillin MIC near 0.25 mg/L, V405A with
Q557E about 1 mg/L, and Q557E is the only substitution NCBI curates for this
species. **This is a flag for confirmatory phenotypic MIC testing, never a
resistance call.** Group B streptococcus remains penicillin-susceptible in close
to all isolates — about 1 % in US surveillance carry elevated MICs or *pbp2x*
variants — and penicillin stays the treatment of choice unless an MIC says
otherwise. The module's note says this on every positive call.

**Numbering is pinned and checked.** The reference is `WP_017647104.1`, 752
residues, which is the accession NCBI's own curated mutation catalogue uses for
its single *S. agalactiae* entry (`pbp2x_Q557E`). `fetch()` verifies the
downloaded sequence: the length, the wild-type residue at all 13 catalogued
positions, and the three motifs. If upstream ever serves a different protein
under that accession the module **fails** rather than report substitutions
against the wrong numbering. One consequence is reported explicitly: the
reference already carries **V** at position 377, so an isolate with I377 differs
from it, while the literature calls I377 the wild type and V377 the
substitution. That difference is reported in its own field as the reference's own
allele showing, not as a finding.

**Named substitutions are limited to the catalogued positions.** An
uncatalogued difference is not named, and the module does not report *gyrA* or
*parC* quinolone-resistance positions at all: no primary source for the specific
residues was verified in time.

---

## 3. Thresholds in one table

| Scheme | Programme | Identity | Coverage | Other |
| --- | --- | --- | --- | --- |
| emm | blastn | ≥ 95 % | alignment within 5 bases of the 180 bp reference | ≤ 4 mismatches, ≤ 2 gap openings, culling limit 5, 500 bp locus clustering |
| cps | blastn | per HSP > 90 %; merged ≥ 90 % | merged ≥ 90 % of the reference | HSP length > 100 bases; ranked by summed bit score |
| gene panels | blastx | 55–90 % per gene (see `GeneSpec.min_identity`) | 25–90 % per gene | intact at ≥ 95 % coverage; in-frame indel tolerance 10 residues; contig-edge window 300 bp |
| hvgA | blastx | ≥ 75 % | ≥ 80 % of 522 aa | separates HvgA from BibA |
| pbp2x | blastx | — | — | substitutions read from the alignment at 13 catalogued positions |
| sag_alp | blastn | ≥ 92 % (≤ 8 % divergence) | ≥ 99 % of the CDC target | winner must beat the next passing target by ≥ 5 points of identity, else `not determined` |
| sag_pilus | blastn | ≥ 92 % (≤ 8 % divergence) | ≥ 99 % of the CDC target | PI-2a and PI-2b together ⇒ `not determined` |

Translated searches run **single-threaded on purpose**
(`strep_common.BLASTX_THREADS`). A sibling agent reported that BLAST+ 2.16.0
`blastx` segfaults when `-max_target_seqs` is combined with `-num_threads > 1`
and the search returns no hits, which is exactly the shape of these searches on a
negative sample; `hydra_amr.engines.blast.blast()` always passes
`-max_target_seqs`. The crash could **not** be reproduced in this environment at
1, 2, 4 or 8 threads against either panel database, so the single thread is a
cheap precaution rather than a confirmed fix — the panels are tens of kilobytes
and the cost is a second or two per assembly. The nucleotide searches are
unaffected and use every thread.

---

## 4. Databases, provenance and licences

### 4.1 `spyogenes/emm` — the CDC emm sequence database

* **Source:** `https://ftp.cdc.gov/pub/infectious_diseases/biotech/tsemm/alltrimmed.tfa`
* **Contents:** about 0.9 MB, one record per named emm subtype (3,443 records at
  retrieval), each the 180-base 5' hypervariable region the CDC protocol
  compares. 407 distinct emm/st types.
* **Version:** the file carries none, so `fetch()` records the HTTP
  `Last-Modified` date in `emm.version` and every call prints it.
* **Licence:** no licence statement accompanies the file. CDC's own terms state:
  *"Unless otherwise noted, all materials on this site are in the public domain
  and may be used and reprinted without special permission."*
* **Redistribution:** permitted.

### 4.2 `sagalactiae/cps` — GBS-SBG

* **Source:** `https://raw.githubusercontent.com/swainechen/GBS-SBG/main/GBS-SBG.fasta`
* **Contents:** 13 sequences, about 90 kB: one per serotype Ia, Ib, II and IV–IX,
  plus four variants of serotype III, each about 6.5–8.3 kb spanning the
  serotype-determining region of *cps*.
* **Version:** the repository has no release tags, so `fetch()` records the
  commit SHA of `main` at retrieval time in `cps.version`.
* **Licence:** **MIT**, Copyright (c) 2021 Swaine Chen: *"Permission is hereby
  granted, free of charge, to any person obtaining a copy of this software and
  associated documentation files ... to deal in the Software without
  restriction"*. The `LICENSE` file is downloaded alongside the data as
  `cps.LICENSE`.
* **Redistribution:** permitted.

### 4.3 `sagalactiae/pbp2x` — the pinned PBP2X reference

* **Source:** NCBI EFetch for protein `WP_017647104.1`.
* **Contents:** one 752-residue protein, about 1 kB.
* **Licence:** NCBI — *"NLM does not claim copyright on the data in the
  Database"*. **Redistribution:** permitted.

### 4.4 `*/panel` — reference proteins from RefSeq annotation

Rather than vendoring a FASTA of unknown provenance, `fetch()` downloads the
RefSeq **feature table** and **protein FASTA** of a pinned set of reference
assemblies and selects, for each panel gene, the protein whose RefSeq gene symbol
(or, failing that, product name) matches — using the expected length to separate
paralogues that share a symbol, as SF370's three *hylP* copies do. The accession
it picked is written to `panel.manifest.json` and carried into every call, so any
result can be traced to a sequence someone else can download.

| Module | Assemblies | Strains |
| --- | --- | --- |
| spyogenes | GCF_000006785.2, GCF_000007425.1, GCF_000007285.1, GCF_000011665.1, GCF_000011765.3, GCF_000012165.2, GCF_000230295.1 | SF370 (M1), MGAS315 (M3), MGAS8232 (M18), MGAS10394 (M6), MGAS5005 (M1), MGAS6180 (M28), Alab49 (M53) |
| sagalactiae | GCF_000196055.1, GCF_000012705.1, GCF_000007265.1, GCF_030489705.1 | NEM316 (III/ST23), A909 (Ia/ST7), 2603V/R (V), COH1 (III/CC17) |

* **Licence:** NCBI RefSeq — *"NLM does not claim copyright on the data in the
  Database"*. **Redistribution:** permitted.
* **Sizes:** the built protein databases are about 60 kB (spyogenes, 49 genes)
  and 30 kB (sagalactiae, 24 genes). The downloads are about 10 MB and 6 MB and
  are cached under `db_dir/<KEY>/refseq/`.

A gene the current RefSeq annotation does not name is **absent from the panel**,
logged at `fetch()` time, recorded in the manifest's `unresolved` list, and
reported as "not determined" — never as a negative. At the time of writing this
affects *cylE* in *S. agalactiae*, and in *S. pyogenes* the genes *ssa*, *sof*,
*spd1*, *spd3* and *sla*, which are not covered (see §6).

### 4.5 `sagalactiae/surface_targets` — the CDC GBS surface-gene target set

* **Source:**
  `https://raw.githubusercontent.com/BenJamesMetcalf/GBS_Scripts_Reference/master/GBS_Reference_DB/GBS_Surface_Gene-DB_Final.fasta`
* **Contents:** 15 targets, 3.4 kB. Four *alp* windows (ALPHA 516 bp, ALP1
  516 bp, ALP23 620 bp, RIB 522 bp), nine 60 bp pilus-island windows (PI1,
  PI1B, PI2A1–PI2A4, PI2B, PI2B2) and SRR1, SRR2 and HVGA. Only the *alp* and
  pilus targets are used: *srr1*, *srr2* and *hvgA* are already called from the
  RefSeq protein panel with their lesions described, which a 60 bp nucleotide
  window cannot do.
* **What it is:** the reference file the CDC *Streptococcus* Laboratory's own
  GBS typing pipeline screens (Metcalf 2017), and which the Sanger Institute's
  `GBS-Typer-sanger-nf` inherits — that pipeline's README states it "has been
  adapted from Ben Metcalf's GBS Typer pipeline" at this repository.
* **Version:** the repository carries no release tags. `fetch()` records the
  commit SHA of `master` in `surface_targets.version`; it was **`ac8fe4a5cb59`,
  2026-08-31**, when this module was written.
* **Licence:** **none stated.** The repository has no `LICENSE` file and no
  licence statement, and the GitHub API reports no licence for it. It is a US
  CDC work product published openly for public-health use, and the Sanger
  Institute redistributes the same file, but **no permission is stated in
  writing**. **Redistribution: unstated** — the integrator should weigh this
  when deciding whether to ship the file inside the Windows zip or fetch it at
  install time. It is 3.4 kB, so fetching it is cheap.
* **Integrity check at fetch time:** `_check_surface_targets()` fails if any
  target the rules name is missing, or if an *alp* target is shorter than
  400 bp. The calls are nothing but these targets and the thresholds, so a
  silent upstream rename must not become a silent negative.

### 4.6 `sagalactiae/alp_alleles` — the published alpha-like protein genes

* **Source:** NCBI EFetch for nucleotide `M97256.1`, `U58333.1`, `AH013348.2`,
  `AF208158.1`, `AF245663.1`, `AJ488912.1`.
* **Contents:** about 16 kB. The six genes as deposited:

  | Accession | Gene | Source strain | Deposited with |
  | --- | --- | --- | --- |
  | M97256.1 | *bca* | **A909** | Michel 1992 |
  | U58333.1 | *rib* | **BM110** | Wästfelt 1996 |
  | AH013348.2 | *alp1* | 515 | unpublished direct submission, Lachenauer & Madoff |
  | AF208158.1 | *alp2* | CJB-110 | Lachenauer 2000 |
  | AF245663.1 | *alp3* | — | Lachenauer 2000 |
  | AJ488912.1 | *alp4* | — | Creti 2004 |

  Two of these come from strains in the validation panel, which is why those
  two rows of the validation table are the strongest expectations available:
  the deposited sequence *is* the published assignment for that strain.
* **Used for:** the `closest_published_allele` second opinion, and naming
  *alp4*, which the CDC target set does not cover. Never for the call itself.
* **Licence:** NCBI — *"NLM does not claim copyright on the data in the
  Database"*. **Redistribution:** permitted.

---

## 5. Validation

### 5.1 *S. pyogenes* — 10 genomes

Expectation source: **emmtyper 0.2.0** run against a BLAST database built from
the *same* current CDC file, so the comparison tests the algorithm rather than
the database vintage; plus the published genome for each strain, and RefSeq's own
pseudogene annotation as an independent check on the regulator calls.

| Accession | Strain | Expected emm (emmtyper) | Obtained | Expected cluster | Obtained | Regulator expected | Obtained |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GCF_000006785.2 | SF370 | emm1.6 | **emm1.6** | A-C3 | **A-C3** | no lesion | **no lesion** |
| GCF_000007425.1 | MGAS315 | emm3.1 | **emm3.1** | A-C5 | **A-C5** | *rocA* lesion (RefSeq pseudogene) | **rocA frameshift** |
| GCF_000007285.1 | MGAS8232 | emm18.19 | **emm18.19** | clade Y | **clade Y** | *rocA* lesion (RefSeq pseudogene) | **rocA premature stop** |
| GCF_000009385.1 | Manfredo | emm5.0 | **emm5.0** | clade Y (umbrella) | **M5** | no lesion | **no lesion** |
| GCF_000011665.1 | MGAS10394 | emm6.4 | **emm6.4** | clade Y (umbrella) | **M6** | no lesion | **no lesion** |
| GCF_000012165.2 | MGAS6180 | emm28.0 | **emm28.0** | E4 | **E4** | no lesion | **no lesion** |
| GCF_000013485.1 | MGAS9429 | emm12.0 | **emm12.0** | A-C4 | **A-C4** | no lesion | **no lesion** |
| GCF_000250905.1 | MGAS15252 | emm59.0 | **emm59.0** | E6 | **E6** | no lesion | **no lesion** |
| GCF_000011765.3 | MGAS5005 | emm1.0 | **emm1.0** | A-C3 | **A-C3** | *covS* deletion (Sumby 2006; RefSeq pseudogene) | **covRS lesion (covS)** |
| GCF_000230295.1 | Alab49 | emm53.0 | **emm53.0** | D4 | **D4** | no lesion | **no lesion** |

10/10 agreement on emm type and subtype. Coverage: emm types 1, 3, 5, 6, 12, 18,
28, 53, 59 and 121/91 (two further draft genomes typed during development), and
clusters A-C3, A-C4, A-C5, D4, E4, E6, M5, M6 and clade Y.

**The two differences, both investigated and neither tuned away.**

* *Clusters M5 and M6 against emmtyper's "Single protein cluster clade Y".*
  Not a disagreement: emmtyper's bundled table reports the clade Y umbrella and
  the CDC records name the specific single-protein cluster. The module reports
  the CDC label, which is the primary source, and prints "a single-protein
  cluster; these are collected as clade Y" alongside it.
* *MGAS15252 (the awkward case).* emmtyper reports `EMM59.0` as the answer and
  `EMM174.1*` as a suspect imposter. The first version of this module reported
  `emm174.1`, because three loci in the emm region match the database and
  `emm174.1` at 1,564,601 is a perfect match, as is `emm59.0` at 1,566,000. The
  sequence explains it: *emm174* is one of the emm-like types carried at the
  *mrp* or *enn* position. The ranking rule was corrected to prefer a true emm
  locus, and the regression is pinned by a test.

Six strains also validate the regulator scheme against RefSeq's independent
pseudogene annotation: SF370, MGAS6180 and Alab49 have intact *covR*, *covS* and
*rocA*; MGAS5005 has a *covS* pseudogene; MGAS8232 and MGAS315 have *rocA*
pseudogenes. 6/6 agreement, with the module naming the lesion type in each case.

The prophage scheme reproduces the published prophage counts of the closed
genomes it was run on: 5 regions in MGAS8232 (published: Φ8232.1–8232.5), 4 in
MGAS315, 3 in SF370 (published: 4; the fourth carries no superantigen cargo and
falls below the three-module floor). It is reported as a *description*, so this is
corroboration rather than a pass/fail criterion.

### 5.2 *S. agalactiae* — 12 genomes in the panel, 20 run

Expectation source: the published serotype for each strain with its citation.
The three serologically non-typeable strains are the awkward case and come from
the population study of capsule loss.

| Accession | Strain | Expected serotype | Source of the expectation | Obtained | CC17 expected | Obtained |
| --- | --- | --- | --- | --- | --- | --- |
| GCF_000196055.1 | NEM316 | III | Glaser 2002 | **III** (variant III-3) | no | **no** |
| GCF_000012705.1 | A909 | Ia | Tettelin 2005 | **Ia** | no | **no** |
| GCF_000007265.1 | 2603V/R | V | Tettelin 2002 | **V** | no | **no** |
| GCF_030489705.1 | COH1 | III | Tettelin 2005; CC17 prototype | **III** (variant III-2) | yes | **yes** (*hvgA*+, *srr2*+) |
| GCF_900155855.1 | BM110 | III | Da Cunha 2014; Tazi 2010 ST-17 prototype | **III** (variant III-2) | yes | **yes** (*hvgA*+, *srr2*+) |
| GCA_000167795.1 | H36B | Ib | Tettelin 2005 | **Ib** | no | **no** |
| GCA_000167715.1 | 18RS21 | II | Tettelin 2005 (553 contigs) | **II** | no | **no** |
| GCF_000427055.1 | ILRI112 | VI | tilapia isolate, serotype VI | **VI** | no | **no** |
| GCF_001217725.1 | GB-NI-002 | IX | Da Cunha 2014 collection, BioSample serotype IX | **IX** | no | **no** |
| GCA_001010335.1 | DK-PW-097 | **serologically non-typeable**; molecular type V, *cpsE* S236N + *cpsM* L94W | Rosini 2015 | **V** (full coverage, 99.9 % identity) | no | **no** |
| GCA_001015895.1 | 325662 | **serologically non-typeable**; molecular type Ia, *cpsE* R308G | Rosini 2015 | **Ia** | no | **no** |
| GCA_001016915.1 | ES-PW-195 | **serologically non-typeable**; molecular type II, *cpsE* R353I | Rosini 2015 | **II** | no | **no** |

12/12 agreement with the published molecular type. Eight further genomes were
typed during development (CJB111 → V, 515 → Ia, CNCTC 10/84 → V, BG-NI-009 → II,
BE-NI-005 → Ib, CUHK_GBS62_918 → VI, ILRI005 → V, strain 290183 → Ia), all
agreeing with their BioSample serotype where one is recorded; 20/20 overall.
Serotypes covered: Ia, Ib, II, III (variants III-2 and III-3), V, VI and IX.

#### 5.2a Alpha-like protein member and pilus islands — the same 12 genomes

Expectation source, strongest first. **(a)** For A909 and BM110 the published
*bca* and *rib* sequences were cloned from, and deposited from, that very
strain (M97256.1 from A909, Michel 1992; U58333.1 from BM110, Wästfelt 1996),
so the defining sequence *is* the expectation. **(b)** For 2603V/R, 18RS21,
H36B and GB-NI-002 the expectation is the gene name in that strain's own genome
annotation: SAG0433 "surface protein Rib" (Tettelin 2002), SAJ_0543 "rib
protein" and SAI_0471 "c protein alpha-antigen precursor" (Tettelin 2005), and
`bca` / "C protein alpha-antigen" (ERS039697_00779). **(c)** For COH1 the
expectation is the published CC17 surface profile, *rib* with PI-2b. **(d)**
For the remaining four no published per-strain assignment was found; those rows
state that and are held only to the rules every isolate must obey.

| Accession | Strain | Expected Alp member | Source class | Obtained | Evidence | Expected PI | Obtained |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GCF_000012705.1 | A909 | *bca* | (a) deposited **from A909** | **bca (alpha-C protein)** | ALPHA 100 % of 516 bp at 100.0 % | — | **PI-1 + PI-2b** |
| GCF_900155855.1 | BM110 | *rib* | (a) deposited **from BM110** | **rib** | RIB 100 % of 522 bp at 100.0 % | PI-2b (CC17) | **PI-1 + PI-2b** |
| GCF_030489705.1 | COH1 | *rib* | (c) CC17 profile | **rib** | RIB 100 % at 100.0 % | PI-2b (CC17) | **PI-1 + PI-2b** |
| GCF_000007265.1 | 2603V/R | *rib* | (b) SAG0433 "surface protein Rib" | **rib** | RIB 100 % at 99.8 % | — | **PI-1 + PI-2a** |
| GCA_000167715.1 | 18RS21 | *rib* | (b) SAJ_0543 "rib protein" | **rib** | RIB 100 % at 100.0 % | — | **PI-1 + PI-2a** |
| GCA_000167795.1 | H36B | *bca* | (b) SAI_0471 "c protein alpha-antigen" | **bca (alpha-C protein)** | ALPHA 100 % at 100.0 % | — | **PI-2a** |
| GCF_001217725.1 | GB-NI-002 | *bca* | (b) `bca` in the record | **bca (alpha-C protein)** | ALPHA 100 % at 100.0 % | — | **PI-2a** |
| GCF_000196055.1 | NEM316 | not *rib*, not *bca* | (d) sequence only | **alp2/alp3** | ALP23 100 % of 620 bp at 99.5 %; RIB 68.3 %, ALPHA 77.6 % | — | **PI-1 + PI-2a** |
| GCA_001010335.1 | DK-PW-097 | — | (d) sequence only; see disagreement below | **alp2/alp3** | ALP23 100 % of 620 bp at 100.0 % | — | **PI-1 + PI-2a** |
| GCA_001015895.1 | 325662 | — | (d) sequence only | **bca (alpha-C protein)** | ALPHA 100 % at 100.0 % | — | **PI-2a** |
| GCA_001016915.1 | ES-PW-195 | — | (d) sequence only | **bca (alpha-C protein)** | ALPHA 100 % at 100.0 % | — | **PI-2a** |
| GCF_000427055.1 | ILRI112 | none (tilapia isolate, outside the human lineages) | (d) the negative case | **no alpha-like protein gene detected** | no CDC *alp* target aligned at all | — | **PI-2b** |

**11/11 agreement where an expectation exists, and the negative case stays
negative.** In every genome exactly one *alp* target passed; the margin to the
next target was 21.4–32.2 percentage points of identity, so nothing was decided
on a narrow call. The independent cross-check — ranking the six deposited
alleles by best alignment — agreed with the CDC call in all twelve.

**Pilus islands: 12/12 on the published rule.** Every genome carried exactly
one of PI-2a and PI-2b and none carried both, which is what the one-locus,
two-variants structure (Rosini 2006) and the 79 % / 21 % split summing to 100 %
(Martins 2013) require. Observed: PI-2a in 8/12, PI-2b in 4/12, PI-1 in 7/12 —
consistent with the published 79 % / 21 % / 70 % at this sample size. Both CC17
strains carried PI-2b, as published. For most individual strains **no published
pilus profile was found**, so those are reported as observed; the rule and the
CC17 association are the validation, and the honest limit is that this scheme
is validated as a *rule-conforming* call rather than against twelve published
per-strain profiles.

**One disagreement, investigated and not tuned away.** DK-PW-097's
alpha-like-protein CDS is annotated "surface protein Rib" in its assembly
record, which contradicts the `alp2/alp3` call. The annotation is automated
(`"Derived by automated computational analysis"`) and the sequence contradicts
it decisively: the gene matches the CDC RIB target at **68.3 %** identity and
the ALP23 target at **100 %** over the full 620 bp, and its closest deposited
allele is *alp3* at 99.9 %. Automated pipelines name every alpha-like-protein
family gene either "surface protein Rib" or "alpha-like surface protein" from
an HMM match, without resolving the member — the same reason this module does
not read the member off the annotation. Where the sequence and the automated
product name agree (2603V/R, 18RS21), the annotation is corroboration; where
they disagree, the sequence decides. This is recorded rather than hidden
because it means **an Alp call read from a RefSeq or PGAP product name is not
reliable**, which is worth knowing beyond this module.

**No disagreement to investigate** in the capsular and CC17 calls, but two
findings worth stating. The three
serologically non-typeable strains are called V, Ia and II at full coverage and
> 99.7 % identity — correct as *genotypes*, and precisely why the output is
labelled capsular genotype. And no genome in the panel is a reduced-penicillin-
susceptibility isolate, so the *pbp2x* scheme was validated only in the negative
direction (no V405A or Q557E called anywhere) plus the `fetch()`-time check that
the reference's wild-type residues match all 13 catalogued positions and all
three motifs. A positive PRGBS control is the missing validation and is listed
in §6.

### 5.3 Running the validation

```
HYDRA_VALIDATION=1 pytest tests/test_typing_strep.py
```

Without that variable the 62 fast tests run offline in under a second; the
validation tests download the genomes and build the databases.

---

## 6. Limits, and what is deliberately not shipped

| Not called | Why |
| --- | --- |
| **M1UK lineage** | defined by a SNP set relative to M1global and assigned phylogenetically; the SNP list was not obtained from a primary source in time. The *nga/ifs/slo* genes are reported instead. |
| **emm pattern for single-protein clusters** | the cluster does not map onto a pattern, and no type→pattern table was verified. Reported as not determined. |
| **FCT region type, T antigen type** | need the pilus backbone allele matched against an allele set that is not published in versioned, downloadable form. |
| **speA and smeZ allele numbers** | the allele series are published but not distributed as a versioned reference set. Residue differences from the named reference are reported instead. |
| **ssa, sof, spd1, spd3, sla** | not named by gene symbol in the current RefSeq annotation of the panel strains, so no traceable reference protein could be selected. Listed as unresolved in the manifest. |
| ***alp2* against *alp3*** | they are mosaics of each other and differ by three bases in 620 inside the only reliable discriminating window. Reported together as `alp2/alp3`, which is also the CDC and national-surveillance convention; the closer deposited allele is given as evidence. Separating them needs the full gene, which the tandem repeats make unreliable. |
| ***alp4*** | not in the CDC target set, and at most 71.6 % identical to any of its targets. An *alp4* isolate is reported as an unnamed family member with *alp4* named among the candidates — never rounded to the nearest target. |
| **Whether a named pilus operon is intact** | the islands are named from 60 bp windows, which is a presence call only. The RefSeq operon genes are reported beside the call; a truncated *pilB* in a PI-2a-positive isolate would not currently be called as such. |
| **Per-strain published pilus profiles** | the island call conforms to the published rules (exactly one of PI-2a/PI-2b; CC17 ⇒ PI-2b) in 12/12 genomes, but published per-strain profiles were found for only the two CC17 strains. A panel with published PI profiles for each strain is the missing validation. |
| **Clonal complex** | follows from the MLST sequence type, which Hydra reports separately. The PubMLST *S. agalactiae* scheme does publish a clonal-complex field; *S. pyogenes* does not. |
| **gyrA / parC quinolone positions** | no primary source for the specific residues was verified. |
| **Prophage boundaries and attachment sites** | an assembly does not resolve them. Regions are reported as the span their marker genes cover, with that stated. |

Known missing validation: a reduced-penicillin-susceptibility *S. agalactiae*
control (a V405A or Q557E isolate), and a serotype VII genome — VII and IX are
91.2 % identical over their full length, so a VII control is the single most
valuable addition to the capsular panel.

---

## 7. Bibliography

Every reference below was checked against PubMed at the time of writing: the
PMID resolves and the title, first author, journal and year match. 48 of the 50
carry a PMID; the two that do not are a reference-laboratory protocol page and a
software repository, cited as such.

### Streptococcus pyogenes (group A)

- **Banks DJ, Porcella SF, Barbian KD, et al.** (2004). Progress toward characterization of the group A Streptococcus metagenome: complete genome sequence of a macrolide-resistant serotype M6 strain. *The Journal of Infectious Diseases 190:727-738*. PMID [15272401](https://pubmed.ncbi.nlm.nih.gov/15272401/); DOI [10.1086/422697](https://doi.org/10.1086/422697).
  <br>*Supports:* MGAS10394, the source of the prtF1/SfbI reference protein.
- **Beall B, Facklam R, Thompson T** (1996). Sequencing emm-specific PCR products for routine and accurate typing of group A streptococci. *Journal of Clinical Microbiology 34:953-958*. PMID [8815115](https://pubmed.ncbi.nlm.nih.gov/8815115/); DOI [10.1128/jcm.34.4.953-958.1996](https://doi.org/10.1128/jcm.34.4.953-958.1996).
  <br>*Supports:* the emm typing scheme: typing from the 5' variable region of emm rather than by M serotyping.
- **Beall B, Facklam RR, Elliott JA, et al.** (1998). Streptococcal emm types associated with T-agglutination types and the use of conserved emm gene restriction fragment patterns for subtyping group A streptococci. *Journal of Medical Microbiology 47:893-898*. PMID [9788813](https://pubmed.ncbi.nlm.nih.gov/9788813/); DOI [10.1099/00222615-47-10-893](https://doi.org/10.1099/00222615-47-10-893).
  <br>*Supports:* emm subtyping and the relationship between emm type and T-agglutination pattern.
- **Beres SB, Sylva GL, Barbian KD, et al.** (2002). Genome sequence of a serotype M3 strain of group A Streptococcus: phage-encoded toxins, the high-virulence phenotype, and clone emergence. *Proceedings of the National Academy of Sciences USA 99:10078-10083*. PMID [12122206](https://pubmed.ncbi.nlm.nih.gov/12122206/); DOI [10.1073/pnas.152298499](https://doi.org/10.1073/pnas.152298499).
  <br>*Supports:* the prophage-encoded toxin content of MGAS315 and the prophage module structure (integrase, packaging, head, tail, lysis) the region description rests on.
- **Bessen DE, Lizano S** (2010). Tissue tropisms in group A streptococcal infections. *Future Microbiology 5:623-638*. PMID [20353302](https://pubmed.ncbi.nlm.nih.gov/20353302/); DOI [10.2217/fmb.10.28](https://doi.org/10.2217/fmb.10.28).
  <br>*Supports:* the emm pattern (A-C, D, E) tissue-tropism grouping the module derives from the emm cluster.
- **Bessen DE, McShan WM, Nguyen SV, et al.** (2015). Molecular epidemiology and genomics of group A Streptococcus. *Infection, Genetics and Evolution 33:393-418*. PMID [25460818](https://pubmed.ncbi.nlm.nih.gov/25460818/); DOI [10.1016/j.meegid.2014.10.011](https://doi.org/10.1016/j.meegid.2014.10.011).
  <br>*Supports:* the Alab49 reference genome, the emm pattern framework and the emm-region gene arrangements the enn call rests on.
- **Biswas I, Scott JR** (2003). Identification of rocA, a positive regulator of covR expression in the group A streptococcus. *Journal of Bacteriology 185:3081-3090*. PMID [12730168](https://pubmed.ncbi.nlm.nih.gov/12730168/); DOI [10.1128/JB.185.10.3081-3090.2003](https://doi.org/10.1128/JB.185.10.3081-3090.2003).
  <br>*Supports:* rocA as a positive regulator of covR, and therefore the meaning of a rocA truncation.
- **Centers for Disease Control and Prevention, Streptococcus Laboratory** (2026). Streptococcus laboratory: emm typing protocol and emm sequence database (alltrimmed.tfa). *CDC, Atlanta; https://www.cdc.gov/strep-lab/php/protocols/*.
  <br>*Supports:* the emm sequence database, the trimmed 180-base hypervariable region, and the type/subtype naming the module reports.
- **Commons RJ, Smeesters PR, Proft T, et al.** (2014). Streptococcal superantigens: categorization and clinical associations. *Trends in Molecular Medicine 20:48-62*. PMID [24210845](https://pubmed.ncbi.nlm.nih.gov/24210845/); DOI [10.1016/j.molmed.2013.10.004](https://doi.org/10.1016/j.molmed.2013.10.004).
  <br>*Supports:* which superantigens are chromosomal and which prophage-encoded, and the clinical associations the module prints.
- **Falugi F, Zingaretti C, Pinto V, et al.** (2008). Sequence variation in group A Streptococcus pili and association of pilus backbone types with Lancefield T serotypes. *The Journal of Infectious Diseases 198:1834-1841*. PMID [18928376](https://pubmed.ncbi.nlm.nih.gov/18928376/); DOI [10.1086/593176](https://doi.org/10.1086/593176).
  <br>*Supports:* the FCT-region pilus genes, and the fact that the pilus backbone is the T antigen.
- **Federle MJ, McIver KS, Scott JR** (1999). A response regulator that represses transcription of several virulence operons in the group A streptococcus. *Journal of Bacteriology 181:3649-3657*. PMID [10368137](https://pubmed.ncbi.nlm.nih.gov/10368137/); DOI [10.1128/JB.181.12.3649-3657.1999](https://doi.org/10.1128/JB.181.12.3649-3657.1999).
  <br>*Supports:* CovR/CsrR as a repressor of capsule and other virulence operons, which is what a covR lesion de-represses.
- **Ferretti JJ, McShan WM, Ajdic D, et al.** (2001). Complete genome sequence of an M1 strain of Streptococcus pyogenes. *Proceedings of the National Academy of Sciences USA 98:4658-4663*. PMID [11296296](https://pubmed.ncbi.nlm.nih.gov/11296296/); DOI [10.1073/pnas.071559398](https://doi.org/10.1073/pnas.071559398).
  <br>*Supports:* SF370, the source of most reference proteins, and its four prophages with their toxin cargo.
- **Frost HR, Davies MR, Velusamy S, et al.** (2020). Updated emm-typing protocol for Streptococcus pyogenes. *Clinical Microbiology and Infection 26:946.e5-946.e8*. PMID [32120034](https://pubmed.ncbi.nlm.nih.gov/32120034/); DOI [10.1016/j.cmi.2020.02.026](https://doi.org/10.1016/j.cmi.2020.02.026).
  <br>*Supports:* the current reference-laboratory emm typing protocol, and the definition of emm-like loci that stops an emm-like type being reported as the emm type.
- **Green NM, Zhang S, Porcella SF, et al.** (2005). Genome sequence of a serotype M28 strain of group A Streptococcus: potential new insights into puerperal sepsis and bacterial disease specificity. *The Journal of Infectious Diseases 192:760-770*. PMID [16088825](https://pubmed.ncbi.nlm.nih.gov/16088825/); DOI [10.1086/430618](https://doi.org/10.1086/430618).
  <br>*Supports:* MGAS6180, the source of the scl1 and srtB reference proteins.
- **Hondorp ER, McIver KS** (2007). The Mga virulence regulon: infection where the grass is greener. *Molecular Microbiology 66:1056-1065*. PMID [18001346](https://pubmed.ncbi.nlm.nih.gov/18001346/); DOI [10.1111/j.1365-2958.2007.06006.x](https://doi.org/10.1111/j.1365-2958.2007.06006.x).
  <br>*Supports:* Mga and RopB as regulators of the emm-region surface proteins and of SpeB.
- **Kreikemeyer B, Klenk M, Podbielski A** (2004). The intracellular status of Streptococcus pyogenes: role of extracellular matrix-binding proteins and their regulation. *International Journal of Medical Microbiology 294:177-188*. PMID [15493828](https://pubmed.ncbi.nlm.nih.gov/15493828/); DOI [10.1016/j.ijmm.2004.06.017](https://doi.org/10.1016/j.ijmm.2004.06.017).
  <br>*Supports:* the fibronectin-binding adhesins PrtF1/SfbI and FbaB/PrtF2 and their role in epithelial invasion.
- **Lukomski S, Nakashima K, Abdi I, et al.** (2000). Identification and characterization of the scl gene encoding a group A Streptococcus extracellular protein virulence factor with similarity to human collagen. *Infection and Immunity 68:6542-6553*. PMID [11083763](https://pubmed.ncbi.nlm.nih.gov/11083763/); DOI [10.1128/IAI.68.12.6542-6553.2000](https://doi.org/10.1128/IAI.68.12.6542-6553.2000).
  <br>*Supports:* the collagen-like surface proteins Scl1/SclA and Scl2/SclB.
- **Lynskey NN, Jauneikaite E, Li HK, et al.** (2019). Emergence of dominant toxigenic M1T1 Streptococcus pyogenes clone during increased scarlet fever activity in England: a population-based molecular epidemiological study. *The Lancet Infectious Diseases 19:1209-1218*. PMID [31519541](https://pubmed.ncbi.nlm.nih.gov/31519541/); DOI [10.1016/S1473-3099(19)30446-3](https://doi.org/10.1016/S1473-3099(19)30446-3).
  <br>*Supports:* the definition of the M1UK lineage, which is why this module reports the nga/slo operon but does not assign the lineage.
- **Mu A, Kwong JC, Isles NS, et al. (Microbiological Diagnostic Unit Public Health Laboratory)** (2019). emmtyper: emm automatic isolate labeller, version 0.2.0. *Software; https://github.com/MDUPHL/emmtyper (GPL-3.0)*.
  <br>*Supports:* the reference tool this module's emm calls were validated against, and the filter thresholds it applies - at most 4 mismatches, at most 2 gap openings, at most 5 bases of alignment-to-subject length difference, 95% identity and 500 bp position clustering - which this module reimplements. No journal article describes the tool; it is cited as software.
- **O'Leary NA, Wright MW, Brister JR, et al.** (2016). Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. *Nucleic Acids Research 44:D733-D745*. PMID [26553804](https://pubmed.ncbi.nlm.nih.gov/26553804/); DOI [10.1093/nar/gkv1189](https://doi.org/10.1093/nar/gkv1189).
  <br>*Supports:* RefSeq, the source of every reference protein in the panel.
- **Sanderson-Smith M, De Oliveira DM, Guglielmini J, et al.** (2014). A systematic and functional classification of Streptococcus pyogenes that serves as a new tool for molecular typing and vaccine development. *The Journal of Infectious Diseases 210:1325-1338*. PMID [24799598](https://pubmed.ncbi.nlm.nih.gov/24799598/); DOI [10.1093/infdis/jiu260](https://doi.org/10.1093/infdis/jiu260).
  <br>*Supports:* the emm cluster system (clusters A-C1 to A-C5, D1 to D5, E1 to E6, single-protein clusters and clade Y) that the module reports.
- **Smoot JC, Barbian KD, Van Gompel JJ, et al.** (2002). Genome sequence and comparative microarray analysis of serotype M18 group A Streptococcus strains associated with acute rheumatic fever outbreaks. *Proceedings of the National Academy of Sciences USA 99:4668-4673*. PMID [11917108](https://pubmed.ncbi.nlm.nih.gov/11917108/); DOI [10.1073/pnas.062526099](https://doi.org/10.1073/pnas.062526099).
  <br>*Supports:* MGAS8232 and the speL/speM reference proteins taken from it.
- **Sumby P, Porcella SF, Madrigal AG, et al.** (2005). Evolutionary origin and emergence of a highly successful clone of serotype M1 group A Streptococcus involved multiple horizontal gene transfer events. *The Journal of Infectious Diseases 192:771-782*. PMID [16088826](https://pubmed.ncbi.nlm.nih.gov/16088826/); DOI [10.1086/432514](https://doi.org/10.1086/432514).
  <br>*Supports:* MGAS5005 and the sda1 reference protein taken from it.
- **Sumby P, Whitney AR, Graviss EA, et al.** (2006). Genome-wide analysis of group A streptococci reveals a mutation that modulates global phenotype and disease specificity. *PLoS Pathogens 2:e5*. PMID [16446783](https://pubmed.ncbi.nlm.nih.gov/16446783/); DOI [10.1371/journal.ppat.0020005](https://doi.org/10.1371/journal.ppat.0020005).
  <br>*Supports:* the covS lesion in MGAS5005 and the clinical meaning of a covRS lesion, which is the rule the regulator call implements.
- **Walker MJ, Barnett TC, McArthur JD, et al.** (2014). Disease manifestations and pathogenic mechanisms of group A Streptococcus. *Clinical Microbiology Reviews 27:264-301*. PMID [24696436](https://pubmed.ncbi.nlm.nih.gov/24696436/); DOI [10.1128/CMR.00101-13](https://doi.org/10.1128/CMR.00101-13).
  <br>*Supports:* the clinical meaning printed for the secreted virulence factors and surface adhesins in the panel.
- **Zhu L, Olsen RJ, Nasser W, et al.** (2015). A molecular trigger for intercontinental epidemics of group A Streptococcus. *The Journal of Clinical Investigation 125:3545-3559*. PMID [26258415](https://pubmed.ncbi.nlm.nih.gov/26258415/); DOI [10.1172/JCI82478](https://doi.org/10.1172/JCI82478).
  <br>*Supports:* the nga/ifs/slo operon, the promoter variants that raise NADase and SLO production, and their role in epidemic emm89 and emm1 lineages.

### Streptococcus agalactiae (group B)

- **Brodeur BR, Boyer M, Charlebois I, et al.** (2000). Identification of group B streptococcal Sip protein, which elicits cross-protective immunity. *Infection and Immunity 68:5610-5618*. PMID [10992461](https://pubmed.ncbi.nlm.nih.gov/10992461/); DOI [10.1128/IAI.68.10.5610-5618.2000](https://doi.org/10.1128/IAI.68.10.5610-5618.2000).
  <br>*Supports:* Sip as a conserved, serotype-independent protective antigen.
- **Carreras-Abad C, Ramkhelawon L, Heath PT, Le Doare K** (2020). A vaccine against group B Streptococcus: recent advances. *Infection and Drug Resistance 13:1263-1272*. PMID [32425562](https://pubmed.ncbi.nlm.nih.gov/32425562/); DOI [10.2147/IDR.S203454](https://doi.org/10.2147/IDR.S203454).
  <br>*Supports:* the hexavalent capsular conjugate serotype set and the protein-based vaccine candidates the module reports coverage against.
- **Cieslewicz MJ, Chaffin D, Glusman G, et al.** (2005). Structural and genetic diversity of group B Streptococcus capsular polysaccharides. *Infection and Immunity 73:3096-3103*. PMID [15845517](https://pubmed.ncbi.nlm.nih.gov/15845517/); DOI [10.1128/IAI.73.5.3096-3103.2005](https://doi.org/10.1128/IAI.73.5.3096-3103.2005).
  <br>*Supports:* the structure of the cps locus, which genes determine the serotype, and the molecular basis of the serotypes Ia to IX.
- **Creti R, Fabretti F, Orefici G, von Hunolstein C** (2004). Multiplex PCR assay for direct identification of group B streptococcal alpha-protein-like protein genes. *Journal of Clinical Microbiology 42:1326-1329*. PMID [15004110](https://pubmed.ncbi.nlm.nih.gov/15004110/); DOI [10.1128/JCM.42.3.1326-1329.2004](https://doi.org/10.1128/JCM.42.3.1326-1329.2004).
  <br>*Supports:* the alpha-like protein family (bca, rib, alp1/epsilon, alp2, alp3, alp4), the equivalence of alp1 and epsilon, and the deposited alp4 sequence AJ488912.1, which is the allele the CDC target set does not cover.
- **Lachenauer CS, Creti R, Michel JL, Madoff LC** (2000). Mosaicism in the alpha-like protein genes of group B streptococci. *Proceedings of the National Academy of Sciences USA 97:9630-9635*. PMID [10944228](https://pubmed.ncbi.nlm.nih.gov/10944228/); DOI [10.1073/pnas.97.17.9630](https://doi.org/10.1073/pnas.97.17.9630).
  <br>*Supports:* alp2 and alp3 as mosaics of one another and of bca and rib, which is why alp2 and alp3 are reported as a pair and never separated; and the deposited alp2 (AF208158.1, from CJB-110) and alp3 (AF245663.1) sequences the module checks the call against.
- **Wästfelt M, Stålhammar-Carlemalm M, Delisse AM, et al.** (1996). Identification of a family of streptococcal surface proteins with extremely repetitive structure. *The Journal of Biological Chemistry 271:18892-18897*. PMID [8702550](https://pubmed.ncbi.nlm.nih.gov/8702550/); DOI [10.1074/jbc.271.31.18892](https://doi.org/10.1074/jbc.271.31.18892).
  <br>*Supports:* rib and its repeat structure, and the deposited rib sequence U58333.1 — which comes from validation strain BM110 and is therefore that strain's published expectation.
- **Metcalf BJ, Chochua S, Gertz RE Jr, et al.** (2017). Short-read whole genome sequencing for determination of antimicrobial resistance mechanisms and capsular serotypes of current invasive Streptococcus agalactiae recovered in the USA. *Clinical Microbiology and Infection 23:574.e7-574.e14*. PMID [28257899](https://pubmed.ncbi.nlm.nih.gov/28257899/); DOI [10.1016/j.cmi.2017.02.021](https://doi.org/10.1016/j.cmi.2017.02.021).
  <br>*Supports:* the CDC GBS whole-genome typing pipeline whose surface-gene target set and thresholds (99 % target coverage, at most 8 % divergence) this module applies to assemblies.
- **McGee L, Chochua S, Li Z, et al.** (2021). Multistate, population-based distributions of candidate vaccine targets, clonal complexes, and resistance features of invasive group B streptococci within the United States, 2015-2017. *Clinical Infectious Diseases 72:1004-1013*. PMID [32060499](https://pubmed.ncbi.nlm.nih.gov/32060499/); DOI [10.1093/cid/ciaa151](https://doi.org/10.1093/cid/ciaa151).
  <br>*Supports:* the reporting convention this module follows — ALPHA, ALP1, ALP23 and RIB as the four alpha-like protein categories in national surveillance, with alp2 and alp3 reported together.
- **Pawlowski A, Lannergård J, Gonzalez-Miro M, et al.** (2022). A group B Streptococcus alpha-like protein subunit vaccine induces functionally active antibodies in humans targeting homotypic and heterotypic strains. *Cell Reports Medicine 3:100511*. PMID [35243418](https://pubmed.ncbi.nlm.nih.gov/35243418/); DOI [10.1016/j.xcrm.2022.100511](https://doi.org/10.1016/j.xcrm.2022.100511).
  <br>*Supports:* why naming the family member matters clinically — the AlpN candidate is built from the N-terminal domains of the alpha-C protein and Rib, so which member an isolate carries is a vaccine-coverage statement.
- **Da Cunha V, Davies MR, Douarre PE, et al.** (2014). Streptococcus agalactiae clones infecting humans were selected and fixed through the extensive use of tetracycline. *Nature Communications 5:4544*. PMID [25088811](https://pubmed.ncbi.nlm.nih.gov/25088811/); DOI [10.1038/ncomms5544](https://doi.org/10.1038/ncomms5544).
  <br>*Supports:* the global GBS collection the validation genomes with published serotypes are drawn from, and the clonal-complex structure of human GBS.
- **Dahesh S, Hensler ME, Van Sorge NM, et al.** (2008). Point mutation in the group B streptococcal pbp2x gene conferring decreased susceptibility to beta-lactam antibiotics. *Antimicrobial Agents and Chemotherapy 52:2915-2918*. PMID [18541727](https://pubmed.ncbi.nlm.nih.gov/18541727/); DOI [10.1128/AAC.00461-08](https://doi.org/10.1128/AAC.00461-08).
  <br>*Supports:* that a single pbp2x point mutation is sufficient to lower beta-lactam susceptibility in this species.
- **Glaser P, Rusniok C, Buchrieser C, et al.** (2002). Genome sequence of Streptococcus agalactiae, a pathogen causing invasive neonatal disease. *Molecular Microbiology 45:1499-1513*. PMID [12354221](https://pubmed.ncbi.nlm.nih.gov/12354221/); DOI [10.1046/j.1365-2958.2002.03126.x](https://doi.org/10.1046/j.1365-2958.2002.03126.x).
  <br>*Supports:* NEM316 (serotype III, ST23), the source of most reference proteins in the panel.
- **Inouye M, Dashnow H, Raven LA, et al.** (2014). SRST2: Rapid genomic surveillance for public health and hospital microbiology labs. *Genome Medicine 6:90*. PMID [25422674](https://pubmed.ncbi.nlm.nih.gov/25422674/); DOI [10.1186/s13073-014-0090-6](https://doi.org/10.1186/s13073-014-0090-6).
  <br>*Supports:* the 90% identity over 90% of the reference length cut-off that GBS-SBG inherits and this module applies.
- **Jerlstrom PG, Chhatwal GS, Timmis KN** (1991). The IgA-binding beta antigen of the c protein complex of Group B streptococci: sequence determination of its gene and detection of two binding regions. *Molecular Microbiology 5:843-849*. PMID [1857207](https://pubmed.ncbi.nlm.nih.gov/1857207/); DOI [10.1111/j.1365-2958.1991.tb00757.x](https://doi.org/10.1111/j.1365-2958.1991.tb00757.x).
  <br>*Supports:* bac, the beta-C protein, and its IgA-binding function.
- **Kimura K, Suzuki S, Wachino J, et al.** (2008). First molecular characterization of group B streptococci with reduced penicillin susceptibility. *Antimicrobial Agents and Chemotherapy 52:2890-2897*. PMID [18490507](https://pubmed.ncbi.nlm.nih.gov/18490507/); DOI [10.1128/AAC.00185-08](https://doi.org/10.1128/AAC.00185-08).
  <br>*Supports:* PBP2X substitutions as the mechanism of reduced penicillin susceptibility in GBS, and therefore why pbp2x divergence is reported even though named substitutions are not called.
- **Koide S, Nagano Y, Takizawa S, Sakaguchi K** (2022). Genomic traits associated with virulence and antimicrobial resistance of invasive group B Streptococcus isolates with reduced penicillin susceptibility from elderly adults. *Microbiology Spectrum 10:e0056822*. PMID [35639001](https://pubmed.ncbi.nlm.nih.gov/35639001/); DOI [10.1128/spectrum.00568-22](https://doi.org/10.1128/spectrum.00568-22).
  <br>*Supports:* the PBP2X substitution table the module names calls from, and the observation that V405A alone accompanies a penicillin MIC near 0.25 mg/L while V405A with Q557E reaches about 1 mg/L.
- **Lamy MC, Dramsi S, Billoet A, et al.** (2006). Rapid detection of the "highly virulent" group B Streptococcus ST-17 clone. *Microbes and Infection 8:1714-1722*. PMID [16822689](https://pubmed.ncbi.nlm.nih.gov/16822689/); DOI [10.1016/j.micinf.2006.02.008](https://doi.org/10.1016/j.micinf.2006.02.008).
  <br>*Supports:* the three variants (A, B, C) of the gbs2018 locus, of which HvgA is variant C and BibA variant A; this is why the hvgA call needs a high coverage floor over the divergent core.
- **Margarit I, Rinaudo CD, Galeotti CL, et al.** (2009). Preventing bacterial infections with pilus-based vaccines: the group B Streptococcus paradigm. *The Journal of Infectious Diseases 199:108-115*. PMID [19086816](https://pubmed.ncbi.nlm.nih.gov/19086816/); DOI [10.1086/595564](https://doi.org/10.1086/595564).
  <br>*Supports:* the PI-2a gene set (pilA, pilB, pilC), the distribution of the islands across lineages, and their use as vaccine candidates.
- **Martins ER, Andreu A, Melo-Cristino J, Ramirez M** (2013). Distribution of pilus islands in Streptococcus agalactiae that cause human infections: insights into evolution and implication for vaccine development. *Clinical and Vaccine Immunology 20:313-316*. PMID [23269415](https://pubmed.ncbi.nlm.nih.gov/23269415/); DOI [10.1128/CVI.00529-12](https://doi.org/10.1128/CVI.00529-12).
  <br>*Supports:* the pilus-island frequencies this module quotes (PI-1 70 %, PI-2a 79 %, PI-2b 21 % of 898 human isolates, every isolate carrying at least one island and PI-2a and PI-2b summing to 100 %, so they are mutually exclusive), and the association of PI-2b with CC17, against which the pilus call is validated.
- **Michel JL, Madoff LC, Olson K, et al.** (1992). Large, identical, tandem repeating units in the C protein alpha antigen gene, bca, of group B streptococci. *Proceedings of the National Academy of Sciences USA 89:10060-10064*. PMID [1438195](https://pubmed.ncbi.nlm.nih.gov/1438195/); DOI [10.1073/pnas.89.21.10060](https://doi.org/10.1073/pnas.89.21.10060).
  <br>*Supports:* bca, the alpha-C protein, and its tandem repeat structure, which is why reference coverage below 100% is expected.
- **O'Leary NA, Wright MW, Brister JR, et al.** (2016). Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. *Nucleic Acids Research 44:D733-D745*. PMID [26553804](https://pubmed.ncbi.nlm.nih.gov/26553804/); DOI [10.1093/nar/gkv1189](https://doi.org/10.1093/nar/gkv1189).
  <br>*Supports:* RefSeq, the source of every reference protein in the panel.
- **Rosa-Fraile M, Dramsi S, Spellerberg B** (2014). Group B streptococcal haemolysin and pigment, a tale of twins. *FEMS Microbiology Reviews 38:932-946*. PMID [24617549](https://pubmed.ncbi.nlm.nih.gov/24617549/); DOI [10.1111/1574-6976.12071](https://doi.org/10.1111/1574-6976.12071).
  <br>*Supports:* the beta-haemolysin/cytolysin and the clinical meaning printed for the virulence genes in the marker panel.
- **Rosini R, Rinaudo CD, Soriani M, et al.** (2006). Identification of novel genomic islands coding for antigenic pilus-like structures in Streptococcus agalactiae. *Molecular Microbiology 61:126-141*. PMID [16824100](https://pubmed.ncbi.nlm.nih.gov/16824100/); DOI [10.1111/j.1365-2958.2006.05225.x](https://doi.org/10.1111/j.1365-2958.2006.05225.x).
  <br>*Supports:* the pilus islands PI-1, PI-2a and PI-2b, their gene content, and that PI-2a and PI-2b are variants of one island at one locus, which is why an isolate carries one of the two and never both.
- **Rosini R, Campisi E, De Chiara M, et al.** (2015). Genomic analysis reveals the molecular basis for capsule loss in the group B Streptococcus population. *PLoS One 10:e0125985*. PMID [25946017](https://pubmed.ncbi.nlm.nih.gov/25946017/); DOI [10.1371/journal.pone.0125985](https://doi.org/10.1371/journal.pone.0125985).
  <br>*Supports:* why a capsular genotype is not a serological serotype: 89 of 96 serologically non-typeable strains carry an inactivating mutation in cps, most often a premature stop in cpsE, which a coverage-and-identity rule cannot see. This is the basis of the caveat printed with every cps call.
- **Six A, Bellais S, Bouaboud A, et al.** (2015). Srr2, a multifaceted adhesin expressed by ST-17 hypervirulent group B Streptococcus involved in binding to both fibrinogen and plasminogen. *Molecular Microbiology 97:1209-1222*. PMID [26094503](https://pubmed.ncbi.nlm.nih.gov/26094503/); DOI [10.1111/mmi.13097](https://doi.org/10.1111/mmi.13097).
  <br>*Supports:* srr2 as the ST-17/CC17-restricted adhesin, which is the marker this module reports in place of hvgA.
- **Slotved HC, Kong F, Lambertsen L, et al.** (2007). Serotype IX, a proposed new Streptococcus agalactiae serotype. *Journal of Clinical Microbiology 45:2929-2936*. PMID [17634306](https://pubmed.ncbi.nlm.nih.gov/17634306/); DOI [10.1128/JCM.00117-07](https://doi.org/10.1128/JCM.00117-07).
  <br>*Supports:* serotype IX, the most recently described serotype in the database.
- **Tazi A, Disson O, Bellais S, et al.** (2010). The surface protein HvgA mediates group B streptococcus hypervirulence and meningeal tropism in neonates. *The Journal of Experimental Medicine 207:2313-2322*. PMID [20956545](https://pubmed.ncbi.nlm.nih.gov/20956545/); DOI [10.1084/jem.20092594](https://doi.org/10.1084/jem.20092594).
  <br>*Supports:* hvgA as the CC17 hypervirulence and meningeal-tropism marker, and therefore why the module states that it cannot call it.
- **Tettelin H, Masignani V, Cieslewicz MJ, et al.** (2002). Complete genome sequence and comparative genomic analysis of an emerging human pathogen, serotype V Streptococcus agalactiae. *Proceedings of the National Academy of Sciences USA 99:12391-12396*. PMID [12200547](https://pubmed.ncbi.nlm.nih.gov/12200547/); DOI [10.1073/pnas.182380799](https://doi.org/10.1073/pnas.182380799).
  <br>*Supports:* 2603V/R (serotype V), a panel source and a validation genome.
- **Tettelin H, Masignani V, Cieslewicz MJ, et al.** (2005). Genome analysis of multiple pathogenic isolates of Streptococcus agalactiae: implications for the microbial pan-genome. *Proceedings of the National Academy of Sciences USA 102:13950-13955*. PMID [16172379](https://pubmed.ncbi.nlm.nih.gov/16172379/); DOI [10.1073/pnas.0506758102](https://doi.org/10.1073/pnas.0506758102).
  <br>*Supports:* A909 (serotype Ia) and the other reference strains, and the serotypes published for them.
- **Tiruvayipati S, Tang WY, Barkham TMS, Chen SL** (2021). GBS-SBG - GBS Serotyping by Genome Sequencing. *Microbial Genomics 7:000688*. PMID [34895403](https://pubmed.ncbi.nlm.nih.gov/34895403/); DOI [10.1099/mgen.0.000688](https://doi.org/10.1099/mgen.0.000688).
  <br>*Supports:* the capsular serotyping database, the serotype-determining cps region it contains, and the calling and uncertainty rules this module implements.
