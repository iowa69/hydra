# Validation

The full record behind the summary in the [README](../README.md): every arm, every comparator, and the reasoning behind each measurement.

## Overview

Measured on 69 genomes (*K. pneumoniae*, *E. coli*, *S. aureus*, *E. faecium*,
*Capnocytophaga*), against independent established implementations run on the
same genomes with the same databases. `tests/compare_with_reference_tools.py`
names the comparators and reproduces the whole table:

| | Result |
|---|---|
| MLST sequence type vs an independent PubMLST typer | **68/69 (98.6%)** identical |
| Loci detected vs an independent screen of the same database | **713/723 (98.6%)** found by both |
| Allele name on shared loci, hits ≥98% identity | **98.5%** identical |
| *E. faecium* ST1478 panel (18 genomes, ST known independently) | **18/18** correct |

The single MLST difference is a *Klebsiella* genome that the comparator types
with the *E. coli* scheme (`ecoli_achtman_4`/ST14464) and Hydra types as
`klebsiella`/ST340 — the EnteroBase *E. coli* scheme has grown alleles that also
match *Klebsiella*, so Hydra confirms the scheme against the species call before
trusting it. Where allele names differ, Hydra generally reports the
higher-identity allele: it ranks candidates by the fraction of the reference
matched identically rather than by bit score, which is biased toward whichever
reference is longest.

Heteroresistance is validated against a synthetic positive control
(`tests/make_heteroresistance_control.py`) — reads simulated from a 23S
reference with a known fraction carrying `23S_G2577T`:

| Simulated | Measured | Called |
|---|---|---|
| 20% (1 of 5 operons) | 0.2004, 93/464 reads, p=2e-152 | heteroresistant, ~1.0/5 operons |
| 5% | 0.0453, 21/464 reads, p=1.1e-21 | heteroresistant |

Reproduce any of it with `tests/compare_with_reference_tools.py`.

### At scale: 5375 *Klebsiella* samples

The table above is 69 mixed genomes. What follows is a single species in depth —
every *Klebsiella* genome and read set on one machine, screened in full detection
mode and compared against the tools people already run, on the same files.

| Arm | Samples | What they are | Completed |
|---|---|---|---|
| A | 667 | closed reference genomes, chromosome **and** plasmids | 667/667 |
| B | 1279 | clinical draft assemblies, 871 with measured EUCAST phenotype | 1279/1279 |
| C | 2933 | published genomes with the sequence type known in advance | 2933/2933 |
| D | 496 | clinical isolates with assembly **and** paired reads | 496/496 |
| | **5375** | | **no failures** |

Every sample, with what was called on it, is listed in
[`docs/validation_klebsiella_samples.tsv`](docs/validation_klebsiella_samples.tsv).

Headline, across the arms that carry independent truth: **2933/2933 sequence types
correct against a recorded answer key**, **666/667, 1278/1279 and 2931/2933 against
`mlst`**, and the whole pipeline at **11.7 s a genome against 18.1 s for the four
reference tools combined**.

Those four sequence-type differences are worth naming individually rather than
rounding away, because they do not all point the same way. Three are genomes `mlst`
declined to type at all because a locus appeared twice — `gapA(51,51)`, `rpoB(4,4)` —
where both copies are the same allele and the profile is unambiguous; Hydra typed each
correctly, confirmed by Kleborate and, on arm C, by the recorded ST. **The fourth is
Hydra getting it wrong**: `KP79_ST512-1LV` matched the EnteroBase *E. coli* scheme at
8/8 exact loci, and Hydra reported *Escherichia coli* for a *Klebsiella* isolate. See
[below](#known-gap-a-scheme-can-outvote-the-species-sketch).

Where Hydra does not lead — genotype-to-phenotype accuracy — the table below says so.

#### 667 closed reference genomes

**667 closed *K. pneumoniae* complex genomes, each a finished chromosome with its
plasmids**, screened in full detection mode (`--preset deep` — nine nucleotide
databases, translated search, point mutations, MLST across every installed scheme,
species and lineage typing) and compared against the same tools run on the same
files. **667/667 completed, no failures.**

| | Result |
|---|---|
| Sequence type vs `mlst` 2.35.0 | **666/667 (99.85%)** identical |
| Sequence type vs Kleborate 3.2.4 | **666/667 (99.85%)** identical |
| Virulence score vs Kleborate | **650/667 (97.45%)** identical |
| Resistance score vs Kleborate | 548/667 (82.16%) identical |
| Acquired AMR genes vs AMRFinderPlus 4.2.7 (same catalogue), gene-family level | **0.805** mean Jaccard |
| Resistance point mutations vs AMRFinderPlus | **0.858** mean Jaccard — 556 shared, **1** only Hydra, 147 only AMRFinderPlus |
| Species | all four members of the complex separated, **667/667** at "strong" confidence |
| Sequence type assigned | 666/667, **363 distinct STs** |
| Genomes completed | **667/667**, no failures |

Every genome screened, with what was called on it, is listed in
[`docs/validation_klebsiella_samples.tsv`](docs/validation_klebsiella_samples.tsv).

The one sequence-type difference is not an error. On `ERR11578643` the comparator
finds *rpoB* twice — `rpoB(4,4)` — and declines to call an ST at all. Both copies are
allele 4, so the profile is unambiguous; Hydra calls ST20 from 7/7 exact loci, and
Kleborate independently calls ST20 as well.

The resistance score agrees less often than the virulence score because it is the one
place a database difference shows through: Kleborate scores from its own curated set,
Hydra from drug classes it carries across thirteen databases, and CARD ships no drug
class for most of its efflux entries (see below).

Where the acquired-gene calls differ, it is usually the name rather than the gene.
Hydra reports the specific allele where AMRFinderPlus reports the generic one —
`fosA5_fam` and `fosA9` (598 calls) against `fosA` (616), `oqxB20`/`oqxB25` against
`oqxB`, `aac(6')-Ib'` against `aac(6')-Ib`. Gene-family agreement of 0.805 therefore
understates it, because a suffix like `_fam` does not reduce to the bare gene.

**Two things this comparison had to avoid**, both of which produce a flattering number
that means nothing. Hydra keeps every overlapping HSP and flags the best per locus, so
counting all rows inflates its gene set roughly sevenfold. And in `--preset deep` the
primary flag is resolved across all thirteen databases at once, so per-database counts
are meaningless — `resfinder` produced 136 rows and zero primaries here, because card
and megares reached those loci first. Genes are therefore compared against
AMRFinderPlus, which shares Hydra's catalogue and nomenclature, and per database only
via `hydra screen -d NAME`, which runs one database the way abricate does.

### One database at a time, against abricate

`hydra screen -d NAME` runs a single database the way abricate does, and writes
abricate's own column layout, so the two are directly comparable on the same 667
genomes. Gene-family level, because where no reference matches exactly the two tools
break the tie to different allele numbers from identical alignments.

| Database | Mean Jaccard | Gene instances abricate found and Hydra did not |
|---|---|---|
| ecoh | **1.000** | 0 |
| *E. coli* VF | **1.000** | 0 |
| MEGARes | **0.988** | 37 |
| ResFinder | **0.985** | 29 |
| ARG-ANNOT | **0.984** | 52 |
| PlasmidFinder | **0.939** | 0 |
| VFDB | **0.903** | 0 |
| CARD | 0.730 | 20 † |
| NCBI | 0.418 | **1** ‡ |

Across nine databases and 667 genomes, the total abricate found and Hydra did not is
139 gene instances — under one per five genomes.

† CARD as-is scores 0.612 with 4011 missed, and 3992 of those are six entries —
`Klebsiella_pneumoniae_KpnH`, `acrA`, `KpnG`, `OmpK37`, `KpnF`, `KpnE` — that exist in
abricate's packaged CARD and **not in the copy `hydra db download` fetches from CARD**.
They are intrinsic efflux and porin genes present in ~100% of *K. pneumoniae*, so they
separate no two isolates. Excluding them leaves 20. The two CARD snapshots differ in
both directions; neither tool failed to detect anything.

‡ NCBI scores low while missing exactly one gene in 667 genomes, which is the shape of
a database difference rather than a detection difference: Hydra's `ncbi` is the
AMRFinderPlus catalogue at 9712 sequences and carries arsenic, mercury and efflux
stress entries; abricate's is 8232 and carries acquired resistance only. Every extra
call is a gene abricate's database cannot contain.

### Tuning the thresholds

Thresholds only ever remove hits, and Hydra writes `%IDENTITY` and `%COVERAGE` on every
row, so the grid is applied to one permissive run rather than re-running the pipeline
64 times per database. Both tools are filtered at the same point — filtering one side
only would measure the filter.

| Database | Best `--min-identity` / `--min-coverage` | Jaccard | at defaults (80/60) |
|---|---|---|---|
| PlasmidFinder | 99 / 80 | 0.939 | 0.746 (**+0.193**) |
| *E. coli* VF | 92.5 / 80 | 1.000 | 0.879 (**+0.121**) |
| NCBI | 99 / 80 | 0.418 | 0.309 (+0.109) |
| VFDB | 80 / 80 | 0.903 | 0.887 (+0.015) |
| ResFinder | 95 / 90 | 0.985 | 0.974 (+0.010) |
| MEGARes, ARG-ANNOT, CARD, ecoh | — | — | within 0.002 of default |

There is no single best setting: the databases that gain are the ones whose entries are
short or near-identical between alleles, where a 60% coverage floor lets a partial
match through. `--min-coverage 80` helps almost everywhere; the identity floor is
worth raising only for PlasmidFinder replicons and NCBI. The defaults are already at or
near the optimum for the five databases that carry most acquired resistance.

### What it costs

Same machine, 24 cores, the same 667 genomes:

| | Total | Per genome |
|---|---|---|
| abricate 1.4.0, nine databases | 2966 s | 4.45 s |
| mlst 2.35.0 | 61 s | 0.09 s |
| Kleborate 3.2.4 | 3446 s | 5.17 s |
| AMRFinderPlus 4.2.7 | 5571 s | 8.35 s |
| **all four together** | **12,044 s** | **18.06 s** |
| **Hydra `--preset deep`** | **7816 s** | **11.72 s** |

Hydra does the whole job in about two thirds of what the four tools cost between them,
and adds point mutations, lineage typing and a species call to what they produce. The
figure is if anything pessimistic: unrelated jobs held about a third of the cores
during part of Hydra's run, and none during the comparators'.

### 1279 clinical isolates

A second arm, drawn from a hospital carbapenem-resistant *K. pneumoniae* collection —
draft assemblies rather than finished genomes, which is what most people actually
have. **1279/1279 completed, 9.7 s a genome.**

| | Result |
|---|---|
| Sequence type vs `mlst` 2.35.0 | **1279/1279 (100%)** with this release's scheme fix; 1278/1279 before it |
| Virulence score vs Kleborate 3.2.4 | **1277/1278 (99.92%)** |
| Sequence type vs Kleborate | 1238/1278 (96.87%) |
| Resistance score vs Kleborate | 1122/1278 (87.79%) |
| Acquired AMR genes vs AMRFinderPlus | 0.652 mean Jaccard |
| Point mutations vs AMRFinderPlus | 2641 shared, **6** only Hydra, 1577 only AMRFinderPlus |

Point-mutation agreement is worse here than on finished genomes (0.68 against 0.86),
in one direction: AMRFinderPlus reports 1577 that Hydra does not, Hydra 6 that
AMRFinderPlus does not. These are draft assemblies, and a mutation whose locus falls
across a contig break is one Hydra's alignment declines to call. It is the same
conservatism that makes it right on closed genomes, and on drafts it costs recall.

The collection is clonal — 41 distinct STs across 1279 isolates, 1265 of them
carbapenemase-positive — so this arm tests depth on one lineage rather than breadth.
Hydra also called one isolate *Escherichia coli* and one "unknown" in a collection
labelled entirely *Klebsiella*.

### 2933 public genomes, with the sequence type known in advance

The third arm is the one with an answer key. Each genome is a published *Klebsiella*
chromosome whose sequence type was recorded independently, and each was rebuilt into a
whole genome first: the chromosome re-joined with its own plasmids, matched on the
strain name in the NCBI defline — 31,189 plasmids returned to 2705 of the 2933
chromosomes. Without that step the arm would be chromosome-only, and most acquired
resistance in *Klebsiella* travels on a plasmid.

| | Result |
|---|---|
| Sequence type vs the recorded ST | **2933/2933 (100.00%)** |
| Sequence type vs Kleborate 3.2.4 | **2933/2933 (100.00%)** |
| Sequence type vs `mlst` 2.35.0 | 2931/2933 (99.93%); of those it typed, **2931/2931 (100%)** |
| Virulence score vs Kleborate | **2917/2933 (99.45%)** |
| Resistance score vs Kleborate | 2453/2933 (83.63%) |
| Acquired AMR genes vs AMRFinderPlus | 0.719 mean Jaccard |
| Point mutations vs AMRFinderPlus | 4269 shared, **17** only Hydra, 2156 only AMRFinderPlus |
| Sequence type assigned | **2933/2933**, 539 distinct STs |
| Species | 2915 *K. pneumoniae*, 15 *variicola*, 3 *quasipneumoniae* |
| Genomes completed | **2933/2933**, no failures, 10.8 s a genome |

The two genomes `mlst` and Hydra disagree on are the same case as on the closed
references, and here the answer key settles it. On `CP110566` the comparator finds
`gapA(51,51)` and on `CP152714` `rpoB(4,4)` — a duplicated locus — and calls no ST at
all. Both copies are the same allele in both genomes, so the profile is unambiguous.
Hydra calls ST491 and ST20, which are the sequence types recorded for those genomes.
Three arms, three instances, same cause, and independent truth agreeing with Hydra
each time.

539 sequence types is the breadth this arm was for, against 363 on the closed
references and 41 on the clinical collection. The commonest are the global
carbapenem-resistant lineages — ST11 (588), ST258 (187), ST15 (144), ST307 (128) — and
ST23 (107), the hypervirulent one. 128 genomes score 5 for virulence, carrying
yersiniabactin, colibactin and aerobactin together.

Median plasmid replicons is 4 here against 2 on the closed references, which is the
plasmid re-joining showing up in the calls rather than in the input.

### Genotype against measured phenotype

Agreement between tools says who resembles whom, not who is right. 871 of these
isolates carry a laboratory EUCAST S/I/R result, so every tool can be scored against
the same measurements. Each is asked the same question using **its own drug annotation**
— a hand-written gene list would encode our judgement and flatter whichever tool
shared it. Intermediates are excluded.

Three drugs are excluded, and the reason is that they cannot separate any two tools.
755 of 773 isolates are meropenem-resistant and 832 of 838 ciprofloxacin-resistant, so
"call everything resistant" earns a perfect sensitivity and a zero specificity and
every tool scores 0.5. Fosfomycin fails the same way from the other end: *fosA* is
intrinsic to *K. pneumoniae*, every tool reports it in every isolate, so every tool
predicts resistance for 100% of them. One is a constant phenotype and the other a
constant prediction; neither is decided by looking at the answer.

| Tool | Mean balanced accuracy | Very major errors | Major errors |
|---|---|---|---|
| AMRFinderPlus 4.2.7 | **0.633** | 284 | 549 |
| **Hydra** | **0.629** | 291 | 547 |
| abricate 1.4.0 (ncbi) | 0.606 | 528 | 469 |
| Kleborate 3.2.4 | 0.587 | 586 | 905 |
| abricate (resfinder) | 0.577 | 744 | 465 |

Hydra and AMRFinderPlus are level, and on four of the five drugs they are identical to
three decimals. That is the expected result rather than a coincidence: both search the
same AMRProt catalogue, so where Hydra adds nothing beyond it they should agree
exactly, and they do. Hydra's remaining margin over abricate and Kleborate comes from
the databases it screens alongside that catalogue.

**An earlier version of this table put Hydra third at 0.580, with 1236 major errors.
That was a measurement error, not a result.** Hydra keeps every overlapping allele at
a locus and flags the best one; the comparison was reading the rejected alternatives
too. On the gentamicin-susceptible isolates Hydra's primary call is `aac(6')-Ib-AKT`,
annotated amikacin/kanamycin/tobramycin and identical to AMRFinderPlus — while the
demoted alternative at the same locus, `aac(6')-Ib'`, is annotated gentamicin. 386 of
those 391 rows were non-primary. Counting them turned correct calls into gentamicin
false positives, and only Hydra was affected, because abricate and AMRFinderPlus emit
one row per locus and have no rejected alternatives to misread. Reading the primary
call removed 689 major errors.

Read the absolute values with care. Every tool sits between 0.58 and 0.63, because
predicting a phenotype from gene presence is genuinely hard: `aac(6')-Ib` is carried by
plenty of amikacin-susceptible isolates, and tigecycline resistance is efflux-regulatory
rather than an acquired gene. Colistin is where Hydra loses most ground — 77 very major
errors against Kleborate's 64 — because resistance there is usually *mgrB* disrupted by
an insertion, and a gene that is absent or interrupted is not something a
presence-and-mutation screen reports. None of these tools was built to be an MIC
predictor, and none of them is one.

### 496 isolates with assembly and reads together

The last arm is the one the other tools cannot be compared on, because they read
assemblies only. 496 clinical isolates were screened with their assembly *and* their
paired reads, which is what makes allele fractions available at all.

| | Result |
|---|---|
| Samples completed | **496/496** (495 assembly+reads, 1 assembly-only — see below) |
| Sequence type vs the same isolates screened from the assembly alone | **496/496 identical** |
| Read-derived variants in resistance genes | 675,202 across 495 samples |
| Point mutations | 970 across 435 samples |
| Catalogued heteroresistance calls | **4** |

Adding reads never changed a sequence type. That is the result to want from this arm:
the read path is additional evidence, not a second opinion that quietly disagrees with
the first.

**Four heteroresistance calls, not 675,202.** The read pass annotates every
intermediate-frequency variant it sees with its allele fraction, and most of those are
a gene differing from its closest reference at 3-5% — sequencing noise or a mixed
culture, not resistance. A heteroresistance call means a *catalogued* resistance
mutation at intermediate frequency, and across 496 *Klebsiella* isolates there were
four, all `blaSHV_C-112A`. That is a property of the catalogue rather than of the
method: *Klebsiella* has exactly one catalogued DNA mutation position, where
*S. aureus* has the 23S rRNA sites this feature was built for. The synthetic control
above shows the method recovering a 5% minority allele in the organism it applies to.

**One isolate lost its reads and kept its assembly.** `KP88_R2.fastq.gz` in this
collection is zero bytes — truncated on download in 2020, one file out of 508. It is
reported as an `assembly` sample carrying the reason in its warnings, rather than
silently absent or silently half-analysed. Before the fix in this release it ended the
whole 25-sample batch.

### Known gap: a scheme can outvote the species sketch

One isolate in 5375 was reported as the wrong genus. `KP79_ST512-1LV` matches the
EnteroBase *E. coli* MLST scheme at 8/8 exact loci — that scheme has accumulated
alleles which also match *Klebsiella*, which is the failure this pipeline was already
designed to guard against — while its Mash sketch says *Klebsiella pneumoniae* at
d=0.0210. The sketch overrides the scheme only below d=0.02, so it missed by a
thousandth, and the isolate was reported as *Escherichia coli* **at "strong"
confidence**, with the *Escherichia* mutation catalogue applied to it. `mlst` also
declined to type this genome, and its allele calls — `gapA(54,144)`, `infB(3,123)` —
suggest a mixed or contaminated assembly, which is why neither method is on firm
ground.

This release downgrades that call: when the scheme and the sketch name different
genera, the confidence becomes `weak` and the evidence string says which two genera
disagree, instead of reporting a contested call as settled.

**It does not change which answer wins.** A sketch outside the override threshold has
not earned the call, and moving that threshold would re-decide the species of every
genome in this validation, so it is left for a release that can be validated on its
own terms. Until then the call is contested and now says so — `species_confidence`
is the field to filter on.

### Known gap: CARD carries no drug class for most efflux entries

1693 of 6052 CARD entries reach Hydra without a drug class, and because the
unannotated ones are largely efflux pumps and their regulators — what *Klebsiella*
mostly carries — 81% of CARD hits here are unclassified. `hydra db download` fetches
CARD's own FASTA, whose defline gives the organism and no drug class; CARD publishes
class separately in its ARO ontology, which abricate's packaged copy has already
joined. `ncbi` and `protein` come from the AMRFinderPlus catalogue and are ~100%
annotated, so acquired resistance still reaches the class summary; CARD's efflux
entries do not.

