# *Salmonella enterica* typing

`src/hydra_amr/typing/organisms/salmonella.py`

The serovar is the unit every *Salmonella* report is written in, and it is a
**phenotype**: an agglutination result against the somatic (O) antigen and the
two phases of the flagellar (H) antigen, read against the
White-Kauffmann-Le Minor (WKL) table of antigenic formulae. An assembly cannot
give that. What it can give is the *genes* that encode those antigens, and the
formula they imply, looked up in the same table.

Every known failure mode of in-silico serotyping lives in that gap. This module
is built so the gap is visible: it reports the antigenic formula, the allele and
coordinates behind each antigen, **every** serovar the formula fits, and what it
could not decide. It never reports a bare negative.

---

## 1. What it calls

| scheme | call | what it means clinically |
|---|---|---|
| `salmonella_serovar` | the antigenic formula and the serovar name(s) it fits | the identification the clinical and public-health report is written in: outbreak linkage, source attribution, notification category |
| `salmonella_subspecies` | I, II, IIIa, IIIb, IV, VI (or *S. bongori*) | the WKL scheme is organised by subspecies; subspecies other than I (*enterica*) are usually reptile- or environment-associated and rare in human disease |
| `salmonella_h2_locus` | `intact` / `fljBA-hin absent` / `partial` | separates a real monophasic variant (4,[5],12:i:- and friends) from an assembly that merely lost *fljB* to a contig break. The distinction is the whole reason the reference multiplex targets the *region*, not the gene |
| `salmonella_spi` | the pathogenicity-island repertoire | SPI-1 and SPI-2 are core: their absence is a QC finding. SPI-7 carries the Vi capsule and is restricted to Typhi, Paratyphi C and Dublin. SGI1 is the DT104 multidrug-resistance island |
| `salmonella_virulence_plasmid` | `spv locus complete` / `partial` / `absent` | the *spv* operon on a pSLT-type virulence plasmid is associated with systemic, bloodstream infection in non-typhoidal serovars |
| `salmonella_typhoidal` | Vi capsule, typhoid toxin, H:z66 | Vi decides Vi-based vaccine relevance and Vi serology; typhoid toxin is the A2B5 genotoxin of the typhoidal serovars; H:z66 is a plasmid-borne second flagellar phase that makes a *S.* Typhi agglutinate as biphasic |
| `salmonella_mlst` | ST and eBurst group (eBG) | the lineage call to read **beside** the serovar. Most serovars are polyphyletic and some eBGs hold several serovars, so the ST is often the more reliable statement about what the isolate *is* |

### Why the serovar and the ST are both reported, and neither is folded into the other

Achtman *et al.* (2012) showed that most serovars are distributed across several
unrelated eBGs and that many eBGs contain more than one serovar. Chattaway
*et al.* (2021) argue from that for keeping the two apart. So this module
reports the serovar from the antigen genes, the ST and eBG from the MLST
profile, and says in the note that a disagreement between them is a finding to
look up in EnteroBase or PubMLST, not an error to average away.

A few ST-specific statements are made where the association is unambiguous and
published: **ST34** with a 4,[5],12:i:- formula is the epidemic monophasic
Typhimurium clone (Petrovska 2016); **ST19** with the same formula is a
monophasic isolate on the ordinary Typhimurium background, which is a different
thing (Soyer 2009); **ST313** is the invasive non-typhoidal lineage of
sub-Saharan Africa. No exhaustive ST-to-serovar table is shipped, because no
machine-readable primary source for one could be verified in the time available
(see *Not implemented*).

---

## 2. Method, with the exact rules

One pass of `blastn` (`-task megablast`) per database, with the assembly's
contigs as the query and the reference set as the BLAST database, so coverage is
measured against the reference. HSPs are merged on the subject axis with
`hydra_amr.engines.blast.merge_hsps`, which is how Hydra already handles a
reference split across a contig gap.

### 2.1 O antigen (serogroup)

Called from `wzx` and `wzy`, the two O-unit genes of the *rfb* cluster that
determine the serogroup (Liu 2014). Thresholds follow SISTR (Yoshida 2016):

* identity &ge; **88 %**, and
* aligned reference bases &ge; **300**, and
* if the aligned length is **300-500 bp**, identity must be &ge; **99 %**.

Both genes are searched. If they agree, that is the group. If only one is
recovered, the call is made from it and the note says which. **If they disagree,
the better-matched gene wins and the note names the other group as not
excluded** - it is never silently dropped.

### 2.2 H antigens (phase 1 and phase 2)

Phase 1 from `fliC`, phase 2 from `fljB`, against the SISTR allele sets whose
sequences are labelled with the WKL antigen they encode.

* **confident**: &ge; 1000 aligned bases at &ge; 99.0 % identity
* **provisional**: &ge; 700 aligned bases at &ge; 97.0 % identity
* otherwise: not determined

A provisional call is reported with the label `(provisional: only N% of the
allele matched)` and, where it applies, `and the match runs to a contig
boundary`. The reason is biological: a flagellin's N- and C-termini are
conserved *across* antigen types, so half an allele can belong to a different
antigen. The value is still reported, because discarding it would lose more than
it protects, but it is marked.

Runners-up are reported **with their margin** - `1,5 at 99.3% identity over
1500 bases` - because an *fljB* locus that matches 1,2 at 100 % and 1,5 at
99.3 % is not the same finding as one that matches both at 100 %.

### 2.3 Is *fljB* absent, or just missing?

This is the question the whole monophasic literature turns on, and an absent
*fljB* **call** is not an absent *fljB* **gene**. The module therefore types the
region, not the gene: `fljB`, `fljA` and `hin` from *S.* Typhimurium LT2, with
`iroB` as the downstream flank.

| observation | H2 reported as | `salmonella_h2_locus` |
|---|---|---|
| an *fljB* allele is called | that antigen | `intact` if `fljA` and `hin` are present too |
| `fljB`, `fljA` and `hin` all absent | `-` (genuinely negative) | `fljBA-hin absent` |
| `fljB` not called but `fljA` or `hin` present | **not determined** | `partial` |

In the widespread monophasic clone the *fljAB* region is physically replaced by
an acquired ASSuT resistance module (García 2016), so *fljA* and *hin* go with
*fljB* - which is why all three are required for a clean negative. When the
region is absent the note points the reader at Hydra's AMR table for the ASSuT
genes (*bla*TEM-1, *strA/strB*, *sul2*, *tet*(B)) and at the ST.

`fliC`, `fljB` and the pBSSB1 z66 allele are homologous, so the marker hits are
deduplicated on the assembly first (`engines.blast.deduplicate`): without that
the LT2 *fljB* reference matches the *fliC* locus over ~55 % of its length at
~92 % identity, and a genome with no *fljB* at all reports *fljB* as "partial".
That was a real bug, caught on *S.* Typhi CT18 during validation.

### 2.4 Serovar lookup

The formula is matched against the WHO 2007 table on
(serogroup, H1, H2, subspecies-free). The key rule is that **H2 absent and H2
undetermined are different queries**:

* `H2 = "-"` matches only rows that are monophasic or that the scheme marks
  `can_h2_be_missing`. It will not pull in biphasic Typhimurium.
* `H2 = not determined` matches any row, and the note says the candidate list is
  "the set consistent with any second phase".

**All** matching serovars are reported, joined by `|`. Diep *et al.* (2019)
found a third of in-silico calls fit more than one serovar; collapsing them to
one name is where the information is lost.

If no row matches an H2-negative formula, the biphasic counterparts with the
same O group and H1 are looked up and the call becomes
`<serovar> monophasic variant`. If no row matches at all, the **formula itself**
is reported - which is what a reference laboratory does with a new antigen
combination - not "unknown".

With neither an O group nor an H1 allele, nothing is looked up and the call is
`not determined`, with the note that the assembly may not be *Salmonella*.

### 2.5 Subspecies

Read from the subspecies column of the WKL rows the formula fits. If they all
agree, that is the subspecies; if they do not, it is `not determined` with the
candidates listed.

**This is an inference from the scheme, not a genomic marker**, and the result
says so in its `method` field. A genome whose formula is unlisted, or whose
formula occurs in more than one subspecies, needs ANI, cgMLST or the biochemical
key. See the RSK2980 row in the validation table for what this looks like in
practice.

### 2.6 Marker genes

Single genes: identity &ge; **90 %**; coverage &ge; **80 %** = present,
40-80 % = partial, below that not reported. 34 markers, each taken from a pinned
GenBank record by gene name or locus tag (`markers.sources.tsv` in the installed
database records which).

* *fljBA-hin* region: `fliC`, `fljA`, `fljB`, `hin`, `iroB`
* O antigen: `oafA` (evidence only - see *Limits*)
* genus confirmation: `invA`
* *spv* locus: `spvR`, `spvA`, `spvB`, `spvC`, `spvD`; pSLT backbone `pefA`,
  `rck`, `traT`, `samA`
* Vi capsule (*viaB*): `tviA`-`tviE`, `vexA`-`vexE`. Vi is called present at
  &ge; 8 of 10 genes; 1-7 genes is reported as **not determined**, because
  partial recovery of a 14 kb region usually means a contig break
* typhoid toxin: `cdtB`, `pltA`, `pltB` - all three required
* SPI-7 `pilL`, `sefA`, `sopE`
* H:z66: the `fljB` allele of the linear plasmid pBSSB1 (Baker 2007)
* `sdf-I`: the Enteritidis-specific difference fragment (Agron 2001), used
  exactly as SeqSero2 uses it (Zhang 2019) - see §3.1

### 2.7 Pathogenicity islands

Whole-island references from the CGE SPIFinder database. Identity &ge; 85 %;
coverage &ge; 90 % = present, 50-90 % = partial, 15-50 % = `trace only - treat
as absent`, below 15 % not reported. `loci_total` is the island count in the
installed database, not the number this genome happened to hit, so
`loci_found/loci_total` means the same thing for every sample.

Coverage is against the SPIFinder reference for each island, and several of those
references are themselves partial sequences, so a coverage figure is comparable
between genomes but is **not** an absolute measure of island length. The note
says so.

### 2.8 MLST

The 7-locus Achtman scheme (`aroC dnaN hemD hisD purE sucA thrA`) from PubMLST,
with its own copy of the alleles and profile table so the module does not depend
on Hydra's PubMLST install. An allele must match at **100 % identity over 100 %
of its length** to be named; otherwise it is reported `~allele` (full length,
novel) or `allele?` (partial) and no ST is assigned. The eBG comes from the
profile table's `clonal_complex` column, which is what PubMLST stores the
eBurst group in.

---

## 3. Limits when working from an assembly

### 3.1 Motility is a phenotype, and non-motile serovars carry the genes anyway

*S.* Gallinarum and *S.* Pullorum are **1,9,12:-:-** - they do not agglutinate
with H antisera - but their genomes carry the same *fliC* `g,m` allele as
*S.* Enteritidis. No gene-based method can distinguish them from Enteritidis,
and this module calls *S.* Gallinarum 287/91 "Enteritidis". It says so: for any
O group containing a non-motile serovar, the note names it and states that it
cannot be excluded from an assembly, only from a motility and agglutination
result. For serogroup D1 it additionally reports the `sdf-I` state, which is the
one genetic handle there is (present supports commonly circulating Enteritidis;
absent leaves Gallinarum and Pullorum open).

### 3.2 Serovars that share one formula cannot be separated by antigens at all

116 groups of formulae in the WHO 2007 table are shared by two or more serovars
(233 serovars in total). The clinically important ones:

| formula | serovars | separated by |
|---|---|---|
| C1:c:1,5 | Choleraesuis, Paratyphi C, Typhisuis | Vi (Paratyphi C), biochemistry, lineage |
| B:b:1,2 | Paratyphi B, Paratyphi B var. Java | d-tartrate fermentation |
| C2-C3:e,h:1,2 | Newport, Bardo | O factors / lineage |
| D1:-:- | Gallinarum, Pullorum | biochemistry |
| C2-C3:d:1,2 | Muenchen, Virginia | biochemistry |
| C2-C3:z10:e,n,x | Hadar, Istanbul | biochemistry |
| D1:l,v:1,5 | Panama, Houston | biochemistry |

These are reported as a candidate list, with the note saying they are told apart
biochemically or by lineage, not by agglutination. **SeqSero2 resolves the
Paratyphi B / var. Java pair genetically** with an STM3356 SNP (Zhang 2019);
that is not implemented here because the exact SNP definition could not be
verified from a primary source in the time available.

### 3.3 O groups that `wzx`/`wzy` cannot separate

Four sets of O groups carry closely related *wzx*/*wzy* alleles (taken from
SISTR's serogroup similarity groups): **E1/E4**, **A/D1/D2**, **C1/F** and
**S/O:62**. When the called group is in one of these, the result carries an
`o_group_ambiguity` field naming the sibling groups and the serovars that would
fit the same H antigens if the O group were one of them. The headline call is
**not** widened, because widening every D1 call to include A and D2 would cost
more precision than it buys.

This is a real limit, and it bites: see RSK2980 in the validation table.

### 3.4 `wzx`/`wzy` are GC-poor

Uelze *et al.* (2020) found that O antigen recovery depends on the sequencing
library preparation, because *wzx* and *wzy* have unusually low GC content. In
8 of 11 validation genomes *wzy* was not recovered at all and the group was
called from *wzx* alone - which the note states every time. An assembly with
patchy coverage over the *rfb* cluster will return an undetermined O group, and
that is the correct answer, not a bug.

### 3.5 O:5 / var. Copenhagen is reported as evidence, not as a call

SeqSero2 flags O:5-negative Typhimurium (var. Copenhagen) from a 7 bp deletion
in `oafA`. This module reports `oafA` coverage, identity and coordinates, and
notes when coverage is below 99 % that a short internal deletion is what makes a
Typhimurium O:5-negative - but it does **not** make the call, because no
O:5-negative genome was available to validate that branch.

### 3.6 Everything else

* A mixed or contaminated assembly can produce two antigen alleles at one locus.
  Zhang *et al.* (2019) found reads from more than one serotype in 3 % of 5,794
  genomes. The runners-up list with its margins is the signal to look for; the
  module does not attempt to deconvolute a mixture.
* The shipped serovar table is the **9th edition (2007)**. Supplements 47 (2010)
  and 48 (2014) added serovars since. A genome of a post-2007 serovar will be
  reported as a formula, correctly, but without its modern name.
* Published serovar metadata in NCBI is itself often wrong (Robertson 2018),
  which is why the validation set below uses strains with a *published* serovar
  and cites the paper.
* Accuracy ceiling: against 1,624 reference-laboratory isolates, SISTR typed
  94 % and SeqSero2 87 % correctly (Uelze 2020). Banerji *et al.* (2020) found
  SeqSero 98 % and 7-gene MLST 95 % alone, >99 % combined - which is the
  argument for reading this module's serovar and ST calls together.

---

## 4. Databases: provenance and licence

| name | source | version | licence | redistribution |
|---|---|---|---|---|
| `antigens` | Zenodo `doi:10.5281/zenodo.15078586`, SISTR cmd database 1.1.3 | published 2025-03-24 | Apache-2.0 (Zenodo record licence `apache2.0`; `sistr_cmd` states `license='Apache 2.0'`) | permitted |
| `serovars` | the same archive: `Salmonella-serotype_serogroup_antigen_table-WHO_2007.csv` | file dated 2024-09-04; content from the WKL scheme 9th ed. (2007) | Apache-2.0 as distributed by SISTR. The underlying table is Grimont & Weill 2007, which **states no licence**; it is published for free download at pasteur.fr | permitted |
| `spi` | `https://bitbucket.org/genomicepidemiology/spifinder_db` | commit `db102668b704`, files dated 2021-01-20, retrieved 2026-11-17 | Apache-2.0 - quoting `README.md`: *"Licensed under the Apache License, Version 2.0"* | permitted |
| `markers` | NCBI Nucleotide: `NC_003197.2`, `NC_003277.2`, `NC_003198.1`, `AM419040.1`, `AF370707.1` | retrieved at first fetch; record versions pinned in the request | **none stated**. NCBI's disclaimer: *"Some NCBI data ... may be protected by copyright ... NCBI itself places no additional restrictions on the use or distribution of the data contained therein."* | permitted |
| `mlst` | PubMLST `pubmlst_salmonella_seqdef` scheme 2 | downloaded at fetch time; unauthenticated access is served data submitted on or before 2024-12-31 | PubMLST Terms & Conditions §3: data submitted on or before 31 December 2024 *"may continue to be downloaded, used, and redistributed without restriction, subject only to proper scientific citation"*. §4 restricts data submitted on or after 1 January 2025 to *"non-commercial academic research or public health surveillance"*; the endpoint this module uses does not serve that data | permitted |

**Sizes.** 13 MB installed. The SISTR archive is 50 MB on the wire because
Zenodo serves no per-member access; only 2.4 MB of it is kept (four antigen
FASTAs and the serovar table) and the archive is deleted. cgMLST330, the Mash
sketch and the full cgMLST allele set - 350 MB of the SISTR database - are not
downloaded, and will not be: the project owner has decided against *Salmonella*
cgMLST ("no cgmlst that is another program"). Core-genome MLST belongs to a
dedicated tool, and this module reports the antigenic formula with its evidence
instead. See section 6.

**Note for whoever decides what ships inside the Windows zip.** Every component
is Apache-2.0 or unrestricted public sequence data, so all five can be
redistributed; Hydra itself is MIT, and Apache-2.0 is compatible with that in
one direction (an MIT project may ship Apache-2.0 data, with the notice). The
PubMLST component is only redistributable because the unauthenticated endpoint
caps the data at 2024-12-31; if the fetch is ever authenticated, that stops
being true and the MLST database must be downloaded on the user's machine
instead of shipped. **SeqSero2 was deliberately not used**: it is GPL-2.0, and
its database and decision tables could not be vendored into an MIT project
without taking the copyleft with them.

---

## 5. Validation

Eleven public assemblies, each with a published serovar and its citation.
Built database 13 MB; about 2 s per genome on four threads.

Run with:

```
HYDRA_VALIDATION=1 pytest tests/test_typing_salmonella.py
```

| accession | strain | expected serovar (source) | obtained formula | obtained serovar | subsp. | ST/eBG | verdict |
|---|---|---|---|---|---|---|---|
| GCF_000006945.2 | LT2 | Typhimurium 1,4,[5],12:i:1,2 (McClelland 2001) | `I 1,4,[5],12:i:1,2` | Typhimurium | I | ST19 / eBG1 | **agree** |
| GCF_000009505.1 | P125109 | Enteritidis PT4 1,9,12:g,m:- (Thomson 2008) | `I 1,9,12:g,m:-` | Enteritidis | I | ST11 / eBG4 | **agree** |
| GCF_000195995.1 | CT18 | Typhi 9,12,Vi:d:- (Parkhill 2001) | `I 9,12[Vi]:d:-` | Typhi | I | ST2 / eBG13 | **agree** |
| GCF_010919335.1 | 119944 | Infantis 6,7,14:r:1,5 (Aviv 2014) | `I 6,7,14:r:1,5` | Infantis | I | ST32 / eBG31 | **agree** |
| GCF_000020705.1 | SL476 | Heidelberg 1,4,[5],12:r:1,2 (Fricke 2011) | `I 1,4,[5],12:r:1,2` | Heidelberg | I | ST15 / eBG26 | **agree** |
| GCF_000016045.1 | SL254 | Newport 6,8,20:e,h:1,2 (Fricke 2011) | `I 6,8,20 or 8:e,h:1,2` | Newport\|Bardo | I | ST45 / eBG3 | **agree** - Bardo shares the formula and is reported with it |
| GCF_000011885.1 | ATCC 9150 | Paratyphi A 1,2,12:a:[1,5] (McClelland 2004) | `I 1,2,12:a:1,5` | Paratyphi A | I | ST85 / eBG11 | **agree** |
| GCF_002210025.1 | TW-Stm6 | Typhimurium var. monophasic 4,5,12:i:- (NCBI record; clone per Petrovska 2016) | `I 1,4,[5],12:i:-` | `I 1,4,[5],12:i:-` | I | **ST34** / eBG1 | **agree** - *fljBA-hin* absent; note names Typhimurium (H2 1,2) as the biphasic counterpart and ST34 as the epidemic clone |
| GCF_000009525.1 | 287/91 | **Gallinarum 1,9,12:-:-, non-motile** (Thomson 2008) | `I 1,9,12:g,m:-` | Enteritidis | I | ST331 / eBG4 | **expected ambiguity** - the antigen *genes* are Enteritidis'. Note names Gallinarum and Pullorum as not excludable and reports `sdf-I` absent, which is the correct genetic signal. §3.1 |
| GCF_000018705.1 | SPB7 | **Paratyphi B 1,4,[5],12:b:1,2** (NCBI record; pair per Grimont & Weill 2007) | `I 1,4,[5],12:b:1,2` | Paratyphi B\|Paratyphi B var. Java | I | ST307 / eBG5 | **expected ambiguity** - identical formula, separated only by d-tartrate. §3.2 |
| GCF_000018625.1 | RSK2980 | **subsp. *arizonae* (IIIa) 62:z4,z23:-** (NCBI RefSeq organism name) | `41:z4,z23:-` | IIIa 41:z4,z23:-\|IV 41:z4,z23:-\|Waycross | not determined | ST2131 / eBG266.2 | **DISAGREE on the O group** - investigated below |

### The RSK2980 discordance, investigated

The published O group is **62**; the module called **S** (O:41). This is not a
threshold artefact: `wzx` matched the group-S allele at 99.8 % identity over
100 % of its length (`NC_010067.1:813868-815154`) and `wzy` the group-S allele at
99.5 % over 100 % (`NC_010067.1:815141-816361`) - the two genes are adjacent in
the *rfb* cluster and agree with each other. No group-O:62 allele in the SISTR
set matched better. So the genome's *rfb* cluster really is closer to the
group-S references in that database than to the O:62 references, and this is the
documented **S / O:62** similarity pair (§3.3), which SISTR also carries.

Two consequences, both reported rather than tuned away:

1. The `o_group_ambiguity` field on this genome reads *"wzx/wzy cannot reliably
   separate O group(s) O:62 from S; were the O group one of those, the same H
   antigens would fit IIIa 62:z4,z23:- (O group O:62)"* - i.e. the correct answer
   is named in the output, as a candidate, with the reason.
2. The subspecies comes back `not determined`, because the three group-S
   candidates span *arizonae* (IIIa) and *houtenae* (IV) and *enterica*. That is
   the honest result of inferring subspecies from a formula (§2.5), and it is why
   an independent subspecies method is on the list below.

A note on SISTR's own copy of this rule: it spells the pair `['S', '62']`, but
the WHO 2007 table spells the group `O:62`, so SISTR's membership test can never
fire for it. This module uses the table's spelling. There is a unit test for
that (`test_similar_o_groups_uses_the_tables_own_spelling`).

### Supporting observations from the same run

These were not the target of the validation but they are a strong check that the
marker and island calls are not noise:

* **SPI-7 was detected in CT18 and in no other genome** - correct, SPI-7 is
  restricted to Typhi, Paratyphi C and Dublin.
* **Vi capsule (*viaB*) only in CT18**; **typhoid toxin in CT18 and Paratyphi A
  ATCC 9150 only** - correct for the typhoidal serovars.
* **The complete *spv* operon only in LT2**, which is the genome that carries
  pSLT; absent from the other ten, which do not.
* In RSK2980 four of five *spv* genes were found at 97-99 % identity, contiguous
  at `NC_010067.1:1438410-1443656` on a 4.6 Mb sequence with **no** pSLT
  backbone gene. The module reports the coordinates and declines to call it a
  plasmid, pointing the reader at Hydra's replicon table.

### Reference tool

**Not validated against a reference tool's own output.** SeqSero2 needs bwa,
samtools, SPAdes and bcftools (conda, not pip) and is GPL-2.0; `sistr_cmd` needs
numpy, pandas, PyTables, scipy and pycurl plus Mash, and the disk budget for
this work was 600 MB. Validation is therefore against the **published serovar of
each reference strain with its citation**, which the brief allows and which is
in some ways a stronger test - it does not inherit the reference tool's own
errors. The antigen allele database and the WHO 2007 serovar table *are* SISTR's,
and the published thresholds are SISTR's, so the method is the same method; what
is not tested is agreement with SISTR's cgMLST refinement step.

---

## 6. Not implemented, and why

| left out | why |
|---|---|
| cgMLST330 and the Mash-based serovar refinement | **decided against by the project owner** ("no cgmlst that is another program"): core-genome MLST is a different tool's job, not Hydra's, so it will not be added however well the zip could carry it. The 350 MB of database was the practical objection; this is the decision. The antigen-only ceiling stands, and section 5 says where it lies |
| An exhaustive ST/eBG -> serovar table | no machine-readable primary source could be verified. The eBG -> serovar table of Achtman 2012 is in a figure/table image, and PubMLST's isolate database would need ~21,000 requests to derive one. The eBG itself is reported; only the named cross-check is missing |
| GenoTyphi *S.* Typhi genotyping (4.3.1/H58 etc.) | **not requested.** The standard for typhoid genomic surveillance, and feasible from an assembly with BLAST over the SNP-flanking regions, so it stays on this list as the single highest-value addition for a typhoid-facing laboratory -- but it has not been asked for and is not being built on spec |
| Paratyphi B dT+ / var. Java discrimination (STM3356 SNP) | SeqSero2 does this; the exact SNP position and alleles could not be verified from a primary source in the time |
| O:5 / var. Copenhagen call from the `oafA` 7 bp deletion | evidence is reported, the call is not: no O:5-negative genome was available to validate it |
| pESI megaplasmid typing for *S.* Infantis | no complete pESI reference whose accession could be verified as *the* published sequence; using the validation genome's own plasmid as the reference would be circular |
| *S.* Enteritidis global/African clade and Typhimurium ST313 sub-lineage SNPs | SNP schemes, same constraint as GenoTyphi |
| Prophage repertoire (Gifsy-1/2/3, Fels-1/2, ST64B, BTP1, mTmV) | detectable, but no curated reference set with a verifiable licence was found in the time; SPIFinder covers the islands, not the prophages |
| Fluoroquinolone and azithromycin resistance determinants | **deliberately not duplicated.** `gyrA`/`parC` QRDR substitutions and the *acrB* R717Q/L azithromycin mutation are in the AMRFinderPlus mutation catalogue under the `Salmonella` taxgroup, which Hydra's existing protein/mutation engines already screen against when the species call maps to `Salmonella`. The integrator should confirm that `--organism Salmonella` is in fact passed; if it is not, that is a pipeline fix, not a typing module |

---

## 7. Bibliography

Every reference below was checked against PubMed (E-utilities `esummary`) and
Crossref at the time of writing: the PMID and DOI resolve, and the first author,
year, journal, volume and pages match. The same list is in `REFERENCES` in the
module, each entry carrying the rule it supports.

**The scheme and its nomenclature**

1. **Grimont PAD, Weill FX.** (2007) *Antigenic Formulae of the Salmonella
   Serovars*, 9th edition. WHO Collaborating Centre for Reference and Research on
   Salmonella, Institut Pasteur, Paris, 166 pp.
   &mdash; the table every serovar name in the output is read out of. A monograph:
   no DOI or PMID; freely downloadable from pasteur.fr.
2. **Guibourdenche M, Roggentin P, Mikoleit M, et al.** (2010) Supplement
   2003-2007 (No. 47) to the White-Kauffmann-Le Minor scheme. *Research in
   Microbiology* 161:26-29. doi:10.1016/j.resmic.2009.10.002, PMID 19840847.
   &mdash; serovars added after the 9th edition.
3. **Issenhuth-Jeanjean S, Roggentin P, Mikoleit M, et al.** (2014) Supplement
   2008-2010 (no. 48) to the White-Kauffmann-Le Minor scheme. *Research in
   Microbiology* 165:526-530. doi:10.1016/j.resmic.2014.07.004, PMID 25049166.
   &mdash; the most recent supplement; the caveat in §3.6 that the shipped table
   is the 2007 edition.
4. **Chattaway MA, Langridge GC, Wain J.** (2021) Salmonella nomenclature in the
   genomic era: a time for change. *Scientific Reports* 11:7494.
   doi:10.1038/s41598-021-86243-w, PMID 33820940.
   &mdash; why a serovar name and a genomic lineage are reported side by side.

**The genetic basis of the antigens**

5. **Liu B, Knirel YA, Feng L, et al.** (2014) Structural diversity in Salmonella
   O antigens and its genetic basis. *FEMS Microbiology Reviews* 38:56-89.
   doi:10.1111/1574-6976.12034, PMID 23848592.
   &mdash; that *wzx* and *wzy* in the *rfb* cluster determine the O serogroup.
6. **Agron PG, Walker RL, Kinde H, et al.** (2001) Identification by subtractive
   hybridization of sequences specific for Salmonella enterica serovar
   Enteritidis. *Applied and Environmental Microbiology* 67:4984-4991.
   doi:10.1128/AEM.67.11.4984-4991.2001, PMID 11679316.
   &mdash; the Sdf-I fragment (GenBank AF370707) used in §3.1.
7. **Baker S, Hardy J, Sanderson KE, et al.** (2007) A novel linear plasmid
   mediates flagellar variation in Salmonella Typhi. *PLoS Pathogens* 3:e59.
   doi:10.1371/journal.ppat.0030059, PMID 17500588.
   &mdash; the H:z66 antigen and its plasmid pBSSB1 (GenBank AM419040).

**The methods this module reimplements, and their limits**

8. **Zhang S, Yin Y, Jones MB, et al.** (2015) Salmonella serotype determination
   utilizing high-throughput genome sequencing data. *Journal of Clinical
   Microbiology* 53:1685-1692. doi:10.1128/JCM.00323-15, PMID 25762776.
   &mdash; SeqSero; the determinant set (*rfb* cluster, *fliC*, *fljB*).
9. **Zhang S, den Bakker HC, Li S, et al.** (2019) SeqSero2: rapid and improved
   Salmonella serotype determination using whole-genome sequencing data.
   *Applied and Environmental Microbiology* 85:e01746-19.
   doi:10.1128/AEM.01746-19, PMID 31540993.
   &mdash; *wzx*/*wzy* and *fliC*/*fljB*; `sdf` for Enteritidis vs Gallinarum;
   `oafA` for O:5-negative Typhimurium; 3 % of genomes hold more than one
   serotype.
10. **Yoshida CE, Kruczkiewicz P, Laing CR, et al.** (2016) The Salmonella In
    Silico Typing Resource (SISTR). *PLoS ONE* 11:e0147101.
    doi:10.1371/journal.pone.0147101, PMID 26800248.
    &mdash; the antigen allele database and the WHO 2007 table installed here,
    and the *wzx*/*wzy* and *fliC*/*fljB* thresholds followed in §2.1-2.2.
11. **Uelze L, Borowiak M, Deneke C, et al.** (2020) Performance and accuracy of
    four open-source tools for in silico serotyping of Salmonella spp. based on
    whole-genome short-read sequencing data. *Applied and Environmental
    Microbiology* 86:e02265-19. doi:10.1128/AEM.02265-19, PMID 31862714.
    &mdash; the accuracy ceiling (SISTR 94 %, SeqSero2 87 % over 1,624 isolates);
    that monophasic variants are where tools disagree; the low-GC *wzx*/*wzy*
    library effect in §3.4.
12. **Banerji S, Simon S, Tille A, et al.** (2020) Genome-based Salmonella
    serotyping as the new gold standard. *Scientific Reports* 10:4333.
    doi:10.1038/s41598-020-61254-1, PMID 32152449.
    &mdash; that antigen and MLST calls should be read together, and that both
    misidentify monophasic Typhimurium.
13. **Diep B, Barretto C, Portmann AC, et al.** (2019) Salmonella serotyping;
    comparison of the traditional method to a microarray-based method and an in
    silico platform using whole genome sequencing data. *Frontiers in
    Microbiology* 10:2554. doi:10.3389/fmicb.2019.02554, PMID 31781065.
    &mdash; a third of in-silico calls fit more than one serovar; the reason for
    the candidate list in §2.4.
14. **Robertson J, Yoshida C, Kruczkiewicz P, et al.** (2018) Comprehensive
    assessment of the quality of Salmonella whole genome sequence data available
    in public sequence databases using SISTR. *Microbial Genomics* 4:e000151.
    doi:10.1099/mgen.0.000151, PMID 29338812.
    &mdash; public serovar metadata is itself often wrong; the validation design
    in §5.

**Monophasic variants and the *fljAB* region**

15. **Soyer Y, Moreno Switt A, Davis MA, et al.** (2009) Salmonella enterica
    serotype 4,5,12:i:-, an emerging Salmonella serotype that represents multiple
    distinct clones. *Journal of Clinical Microbiology* 47:3546-3556.
    doi:10.1128/JCM.00546-09, PMID 19741087.
    &mdash; 4,[5],12:i:- is not one lineage, so the antigens do not identify the
    clone.
16. **García P, Malorny B, Rodicio MR, et al.** (2016) Horizontal acquisition of
    a multidrug-resistance module (R-type ASSuT) is responsible for the
    monophasic phenotype in a widespread clone of Salmonella serovar
    4,[5],12:i:-. *Frontiers in Microbiology* 7:680.
    doi:10.3389/fmicb.2016.00680, PMID 27242707.
    &mdash; the *fljAB* region is physically replaced, so *fljA* and *hin* go with
    *fljB*: the basis of the three-gene rule in §2.3.
17. **Prendergast DM, Hand D, Nic Ghallchoir E, et al.** (2013) A multiplex
    real-time PCR assay for the identification and differentiation of Salmonella
    enterica serovar Typhimurium and monophasic serovar 4,[5],12:i:-.
    *International Journal of Food Microbiology* 166:48-53.
    doi:10.1016/j.ijfoodmicro.2013.05.031, PMID 23827807.
    &mdash; that a reference method interrogates the region, not the gene; the
    reason `salmonella_h2_locus` exists as its own call.
18. **Petrovska L, Mather AE, AbuOun M, et al.** (2016) Microevolution of
    monophasic Salmonella Typhimurium during epidemic, United Kingdom, 2005-2010.
    *Emerging Infectious Diseases* 22:617-624. doi:10.3201/eid2204.150531,
    PMID 26982594.
    &mdash; the ST34 epidemic clone named when the monophasic formula and ST34
    coincide.

**Lineage, MLST and the eBurst group**

19. **Achtman M, Wain J, Weill FX, et al.** (2012) Multilocus sequence typing as
    a replacement for serotyping in Salmonella enterica. *PLoS Pathogens*
    8:e1002776. doi:10.1371/journal.ppat.1002776, PMID 22737074.
    &mdash; the eBG concept, and that most serovars are polyphyletic.
20. **Alikhan NF, Zhou Z, Sergeant MJ, Achtman M.** (2018) A genomic overview of
    the population structure of Salmonella. *PLoS Genetics* 14:e1007261.
    doi:10.1371/journal.pgen.1007261, PMID 29621240.
    &mdash; serovar-to-lineage concordance at genome scale; the authority a
    serovar/ST disagreement should be checked against.
21. **Jolley KA, Bray JE, Maiden MCJ.** (2018) Open-access bacterial population
    genomics: BIGSdb software, the PubMLST.org website and their applications.
    *Wellcome Open Research* 3:124. doi:10.12688/wellcomeopenres.14826.1,
    PMID 30345391.
    &mdash; the allele and profile data the MLST call uses; the citation PubMLST
    asks for.

**Pathogenicity islands**

22. **Roer L, Hendriksen RS, Leekitcharoenphon P, et al.** (2016) Is the evolution
    of Salmonella enterica subsp. enterica linked to restriction-modification
    systems? *mSystems* 1:e00009-16. doi:10.1128/mSystems.00009-16,
    PMID 27822532.
    &mdash; SPIFinder and its island reference database, which `salmonella_spi`
    is called against.
23. **Yoon SH, Park YK, Kim JF.** (2015) PAIDB v2.0: exploration and analysis of
    pathogenicity and resistance islands. *Nucleic Acids Research* 43:D624-D630.
    doi:10.1093/nar/gku985, PMID 25336619.
    &mdash; the island boundaries and names SPIFinder inherits.

**Reference sequences and the validation strains**

24. **McClelland M, Sanderson KE, Spieth J, et al.** (2001) Complete genome
    sequence of Salmonella enterica serovar Typhimurium LT2. *Nature*
    413:852-856. doi:10.1038/35101614, PMID 11677609.
    &mdash; source of the *fljBA-hin*, `oafA`, `invA` and pSLT *spv* marker
    sequences, and the published serovar of LT2.
25. **Parkhill J, Dougan G, James KD, et al.** (2001) Complete genome sequence of
    a multiple drug resistant Salmonella enterica serovar Typhi CT18. *Nature*
    413:848-852. doi:10.1038/35101607, PMID 11677608.
    &mdash; source of the *viaB*, typhoid toxin, SPI-7 and `sefA` markers, and
    the published serovar of CT18.
26. **Thomson NR, Clayton DJ, Windhorst D, et al.** (2008) Comparative genome
    analysis of Salmonella Enteritidis PT4 and Salmonella Gallinarum 287/91.
    *Genome Research* 18:1624-1637. doi:10.1101/gr.077404.108, PMID 18583645.
    &mdash; the published serovars of P125109 and 287/91, and that Gallinarum
    carries the Enteritidis antigen genes while being non-motile (§3.1).
27. **Fricke WF, Mammel MK, McDermott PF, et al.** (2011) Comparative genomics of
    28 Salmonella enterica isolates: evidence for CRISPR-mediated adaptive
    sublineage evolution. *Journal of Bacteriology* 193:3556-3568.
    doi:10.1128/JB.00297-11, PMID 21602358.
    &mdash; the published serovars of SL476 (Heidelberg) and SL254 (Newport).
28. **Aviv G, Tsyba K, Steck N, et al.** (2014) A unique megaplasmid contributes
    to stress tolerance and pathogenicity of an emergent Salmonella enterica
    serovar Infantis strain. *Environmental Microbiology* 16:977-994.
    doi:10.1111/1462-2920.12351, PMID 24320043.
    &mdash; the published serovar of strain 119944, and the pESI megaplasmid of
    that lineage.
29. **McClelland M, Sanderson KE, Clifton SW, et al.** (2004) Comparison of
    genome degradation in Paratyphi A and Typhi, human-restricted serovars of
    Salmonella enterica that cause typhoid. *Nature Genetics* 36:1268-1274.
    doi:10.1038/ng1470, PMID 15531882.
    &mdash; the published serovar of ATCC 9150 (Paratyphi A).
