<h1 align="center">Hydra</h1>

<p align="center">
  <strong>One pass over one set of databases: acquired resistance and virulence genes,
  resistance point mutations, heteroresistance from reads, MLST, and lineage typing —
  for assemblies or raw reads, one sample or a thousand.</strong>
</p>

<p align="center">
  <a href="#install">Install</a> ·
  <a href="#quick-start">Quick start</a> ·
  <a href="#what-hydra-does">What it does</a> ·
  <a href="#heteroresistance">Heteroresistance</a> ·
  <a href="#outputs">Outputs</a> ·
  <a href="#presets">Presets</a> ·
  <a href="#command-reference">Commands</a> ·
  <a href="#performance">Performance</a> ·
  <a href="#how-it-works">Methods</a>
</p>

---

## Why

Characterising a bacterial isolate normally means running four or five tools, each
with its own database, its own thresholds, and its own output layout, then writing
a script to glue the answers together. Hydra does the whole job in one command and
one table, and adds the one thing none of those tools can do from an assembly:
**it finds resistance mutations that are present in only some copies of a gene.**

## What Hydra does

| | |
|---|---|
| **Acquired genes** | NCBI, CARD, ResFinder, ARG-ANNOT, MEGARes, VFDB, PlasmidFinder, EcOH, *E. coli* VF — screened together, with drug-class annotation transferred onto every database |
| **Translated search** | AMRFinderPlus protein reference, with `--plus` stress-response and virulence elements |
| **Point mutations** | Organism-specific protein and DNA catalogues: *gyrA*, *parC*, *rpoB*, *pmrB*, 16S/23S rRNA, promoters such as `pbp4_T-266A`, and mosaic-PBP calls |
| **Heteroresistance** | Allele fractions measured from reads at every catalogued position, with an estimate of how many rRNA operons carry the mutation |
| **MLST** | All 167 installed PubMLST schemes searched at once; the scheme is chosen automatically, no `--species` needed |
| **Species** | Mash sketches and the MLST result combined — the sketch also runs on raw reads, so a FASTQ-only sample still gets an organism |
| **Lineage typing** | Yersiniabactin, colibactin, aerobactin, salmochelin and *rmp* sublineages, *wzi*, *E. coli* serotyping (`O121:H7`), pathovar markers and the Achtman/Pasteur/Lee MLST schemes |
| **Scores** | Kleborate's 0–5 virulence score, and a 0–3 resistance score generalised to any Gram-negative via drug-class annotation |
| **Reads** | FASTQ in directly — gene detection by depth and breadth, MLST called from reads mapped to the scheme's loci, mutations in each resistance gene measured against its closest reference, optional assembly |

## Install

```bash
conda create -n hydra -c conda-forge -c bioconda hydra-amr
conda activate hydra
```

Then install the reference databases. If you already have `abricate`,
`ncbi-amrfinderplus`, `mlst` or `kleborate` in conda environments on the machine,
Hydra finds and converts them:

```bash
hydra db import
hydra db list
```

Otherwise `hydra db download` prints the upstream source and licence of every
database so you can fetch them, then `hydra db import --source DIR` converts them.
Databases live in `$HYDRA_DB` (default `~/.hydra/db`).

<details>
<summary>Install from source</summary>

```bash
git clone https://github.com/iowa69/hydra && cd hydra
conda create -n hydra -c conda-forge -c bioconda \
    python=3.11 blast minimap2 samtools mash pysam pandas numpy scipy
conda activate hydra
pip install -e .
```
</details>

## Quick start

```bash
# one isolate, everything on
hydra run -a isolate.fasta -o results/

# a directory of assemblies, with heatmaps and matrices
hydra run assemblies/ -o results/ --preset surveillance

# paired reads: linezolid heteroresistance from 23S allele fractions
hydra run -1 s_R1.fq.gz -2 s_R2.fq.gz --organism Staphylococcus_aureus \
    --preset linezolid -o results/

# assembly and reads together: consensus calls plus allele-fraction resolution
hydra run -a s.fasta -1 s_R1.fq.gz -2 s_R2.fq.gz -o results/

# abricate-style single-database screen straight to stdout
hydra screen -d card assemblies/*.fasta --stdout
```

Inputs can be given as files, directories (scanned recursively, FASTA and FASTQ
sorted out automatically), `--r1/--r2` pairs, or a sample sheet:

```bash
hydra run --input-list samples.tsv -o results/
```
```
# sample    assembly            R1                  R2
KP001      KP001.fasta         KP001_R1.fq.gz      KP001_R2.fq.gz
KP002      KP002.fasta
KP003                          KP003_R1.fq.gz      KP003_R2.fq.gz
```

## Heteroresistance

Resistance to linezolid is usually **heteroresistant**: the 23S rRNA mutation that
confers it sits in only some of the four to six rRNA operons a genome carries.
The assembly consensus therefore shows the wild-type base, and every
assembly-only caller reports a susceptible isolate.

Hydra aligns the reads to the organism's reference loci and reports the observed
allele fraction at each catalogued position:

```bash
hydra run -1 s_R1.fq.gz -2 s_R2.fq.gz -O Staphylococcus_aureus --preset linezolid -o out/
```

```
sample  gene  mutation    class          allele_fraction  depth  status
s       23S   23S_G2577T  OXAZOLIDINONE  0.2004           464    heteroresistant
```

with the detail column recording `AF=0.200; 93/464 reads; HETERORESISTANT;
~1.0/5 operons; p=2.0e-152`. The p-value is the probability of seeing that many
minority-allele reads from sequencing error alone.

A mutation at or above `--fixed-allele-fraction` (default 0.90) is reported as
fixed; between `--min-allele-fraction` (default 0.02) and that, as
heteroresistant. Use `--report-absent-sites` to see the depth and fraction at
every catalogued position, including the negative ones — useful for confirming a
site was actually interrogated rather than missed.

`-O/--organism` is optional: a FASTQ-only sample is identified by sketching its
reads, so it reaches the right mutation catalogue on its own. Pass `--organism`
when the species is outside the installed sketches, or to force a choice.

The same machinery works for any catalogued DNA mutation in any supported
organism, not only linezolid: `hydra run --list-organisms` lists them.

## Reads without an assembly

A FASTQ-only sample is not a second-class input. Hydra sketches its reads to
identify the species, then:

* **types it** — reads are mapped to one representative allele per locus of the
  schemes for that species, the consensus is read off the pileup, and each locus
  is matched back against every allele of the scheme. On seven isolates with
  closed references this reproduced the assembly-based sequence type exactly,
  7/7, including *E. coli* ST5082 and *K. pneumoniae* ST37;
* **screens its resistance genes** — genes are called by depth and breadth, then
  each gene's reads are compared with the closest reference Hydra holds, and
  every difference is reported with the fraction of reads supporting it.

Both are labelled as read-derived. `mlst_source` says `reads` rather than
`assembly`, and every variant row names the reference it was measured against,
because that is what the call means: variation relative to the nearest sequence
in the database, not to the isolate's own assembled gene. A variant carried by
all the reads is summarised as one row — the gene is simply a different allele —
while a *minority* variant is reported on its own, since that is a mixed
population and the assembly consensus would hide it.

```
gene      change   AF     depth  note
gyrA      S83L     0.24   112    vs closest reference NG_050497.1; HETERORESISTANT; catalogued as gyrA_S83L
tet(A)    -        -      40     differs from its closest reference NG_048164.1 at 3 protein-changing positions
```

Variant calling is gated by `--min-allele-reads` and a binomial test against a
sequencing-error background, so a single mismatching read at 25× is not reported
as a mutation. Add `--report-synonymous` for silent changes, or
`--no-reads-variants` / `--no-reads-mlst` to switch either off.

## Outputs

`-o DIR` writes a set of tables named after `--prefix` (default `hydra`):

| File | Contents |
|---|---|
| `hydra.tsv` | long format, one row per detected element |
| `hydra.summary.tsv` | one row per sample: species, ST, typing, counts, scores, QC |
| `hydra.matrix.tsv` | pivoted sample × gene matrix |
| `hydra.classes.tsv` | sample × drug-class recap |
| `hydra.mlst.tsv` | scheme, ST and every allele call |
| `hydra.typing.tsv` | lineage schemes and scores |
| `hydra.heteroresistance.tsv` | allele fractions at mutation sites |
| `hydra.html` | self-contained report with clustered heatmaps |
| `hydra.json` | everything, machine-readable |

Choose formats with `-f/--format` (`tsv`, `csv`, `json`, `html`, `xlsx`, plus
`abricate` and `amrfinder` for drop-in compatible column layouts).

Control what the matrix contains:

```bash
--cell binary      # 1/0 presence-absence (default)
--cell identity    # % identity of the best hit
--cell coverage    # % of the reference covered
--cell count       # number of copies
--cell depth       # mean read depth
--cell fraction    # allele fraction
--cell symbol      # "identity/coverage"

--rows sample --columns gene      # default
--rows sample --columns class     # drug-class recap
--rows gene   --columns sample    # transposed
--element-types AMR,VIRULENCE     # restrict to certain element types
```

When several databases report the same locus, every row is kept in the long
table but only one is marked `primary`; counts, matrices and scores use the
primary rows, so a gene present once is counted once no matter how many
databases found it. `--report-overlaps` turns this off.

## Presets

```bash
hydra presets            # list them
hydra presets linezolid  # show one in detail
```

| Preset | For |
|---|---|
| `fast` | quickest useful screen: NCBI genes only |
| `standard` | the default: genes, translated search, mutations, MLST, typing |
| `deep` | every database, `--plus` elements, everything on |
| `surveillance` | multi-sample: all matrices, class recap, HTML heatmaps |
| `amr` | resistance only |
| `virulence` | VFDB and the lineage schemes |
| `linezolid` | 23S heteroresistance from reads, plus *cfr*/*optrA*/*poxtA* |
| `gram-positive` | staphylococci, enterococci, streptococci |
| `enterobacterales` | *Klebsiella*/*E. coli*: AMR, virulence loci, lineage scores |
| `plasmid` | replicon typing |
| `abricate` | abricate's thresholds and column layout |
| `amrfinder` | AMRFinderPlus's behaviour and column layout |

A preset only sets defaults; anything you pass explicitly wins.

## Command reference

```
hydra run       full analysis of assemblies and/or reads
hydra screen    acquired-gene screening only (abricate-style)
hydra db        list | import | download | info | check | remove
hydra presets   list the available presets
```

<details>
<summary><code>hydra run</code> options</summary>

**Inputs** — `INPUT...` (files or directories), `-a/--assembly`, `-1/--r1`,
`-2/--r2`, `--reads`, `--input-list`, `--name`

**Databases** — `-d/--db` (names, or the groups `all`, `standard`, `amr`,
`virulence`, `nucl`, `core`), `--list-databases`

**Analysis** — `--preset`, `-O/--organism`, `--list-organisms`,
`--auto-organism/--no-auto-organism`, `--plus`, `--mlst/--no-mlst`, `--scheme`,
`--typing/--no-typing`, `--protein/--no-protein`,
`--point-mutations/--no-point-mutations`,
`--heteroresistance/--no-heteroresistance`,
`--reads-mlst/--no-reads-mlst`, `--reads-variants/--no-reads-variants`,
`--report-synonymous`, `--assemble`

**Thresholds** — `--min-identity`, `--min-coverage`, `--protein-min-identity`,
`--protein-min-coverage`, `--min-contig-length`, `--report-overlaps`

**Reads** — `--min-depth`, `--min-gene-breadth`, `--min-allele-fraction`,
`--fixed-allele-fraction`, `--min-allele-reads`, `--min-base-quality`,
`--report-absent-sites`

**Output** — `-o/--outdir`, `--prefix`, `-f/--format`, `--stdout`, `--cell`,
`--rows`, `--columns`, `--element-types`, `--title`

**Runtime** — `-t/--threads`, `-j/--jobs`, `--db-dir`, `--tmpdir`, `--keep-temp`,
`-v/--verbose`, `-q/--quiet`
</details>

## Performance

Hydra screens every sample in a run with a *single* search per database rather
than one process per sample, and splits long contigs into overlapping chunks so
BLAST can spread a closed chromosome across all cores. On 24 cores:

| | Hydra | Reference tool |
|---|---|---|
| One 5.8 Mb *K. pneumoniae* genome, full pipeline | **13 s** | 17 s (AMRFinderPlus, AMR only) |
| 69 mixed genomes, full pipeline | **6 min** (5.3 s/genome) | — |
| 1.5 Gbp of paired reads, gene calls + 23S allele fractions | **9 s** | — |

The full pipeline means two nucleotide databases, translated search, protein and
DNA point mutations, MLST across all 167 schemes, species identification and
lineage typing.

## Validation

Measured on 69 genomes (*K. pneumoniae*, *E. coli*, *S. aureus*, *E. faecium*,
*Capnocytophaga*), against the tools Hydra replaces:

| | Result |
|---|---|
| MLST sequence type vs `mlst` | **68/69 (98.6%)** identical |
| Loci detected vs `abricate` (same database) | **713/723 (98.6%)** found by both |
| Allele name on shared loci, hits ≥98% identity | **98.5%** identical |
| *E. faecium* ST1478 panel (18 genomes, ST known) | **18/18** correct |

The single MLST difference is a *Klebsiella* genome that `mlst` types with the
*E. coli* scheme (`ecoli_achtman_4`/ST14464) and Hydra types as
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

## How it works

**Reference coverage, not query coverage.** HSPs between one contig and one
reference gene are merged on the *reference* axis, so a gene interrupted by an
assembly gap still reports its true coverage instead of appearing twice at half
length.

**One search, many samples.** Contigs from every sample go into one query file
with synthetic ids; results are mapped back afterwards. This removes N process
starts and N database loads, and gives BLAST enough query sequences to use every
core.

**Organism-aware mutations.** Species is inferred per sample, so a mixed-species
batch is still evaluated against each sample's own mutation catalogue from a
single translated search. Reference proteins that exist only to anchor mutations
(*gyrA*, *rpoB*, the PBPs) are excluded from acquired-gene calls — they are
housekeeping genes present in every genome.

**Alignment-verified calls.** A mutation is only called when the alignment
actually lands on the catalogued residue: if the reference base at that position
disagrees with the catalogue, the call is refused rather than guessed. Negative
DNA positions are resolved from the end of the reference record, which is how
promoter mutations are catalogued.

**Profile-aware MLST.** Which columns of a profile table are loci is decided by
the scheme's own allele files, not by guessing which trailing column names look
like metadata — PubMLST inserts columns such as `MLST_cluster` between the loci
and `clonal_complex`, and a guess shifts every allele by one. Alleles tying on
identity and coverage are resolved by lowest allele id, matching PubMLST
convention. A locus that is not found is retried as allele `0`, which is how
lineages such as *E. faecium* CC17 — where *pstS* is deleted — get their ST; when
other profiles fit the alleles that *were* found and differ only at the absent
locus, the note names them, so an ST resting on an absent locus is never quietly
presented as certain.

**Scores that mean the same thing everywhere.** The resistance score ignores
intrinsic and AmpC-type β-lactamases when testing for ESBL — otherwise every
*E. coli*, which all carry *blaEC*, would score 1 — and falls back to gene names
when a database ships no curated subclass, so the score does not change with the
`--db` list. The virulence score counts virulence loci only; a sequence type, a
serotype or a *wzi* allele is not virulence.

## Testing

```bash
python -m pytest tests/ -q
```

`tests/make_heteroresistance_control.py` builds a synthetic positive control:
reads simulated from a 23S reference with a defined fraction carrying a
resistance mutation. Hydra recovers 0.20 as 0.2004 and 0.05 as 0.0453.

`tests/compare_with_reference_tools.py` reports concordance with `abricate` and
`mlst` on the same genomes.

## Citing

Hydra orchestrates and reimplements; the databases and schemes are other
people's work, and they should be cited. `hydra db info NAME` prints the
citation and licence for any installed database, and `hydra db list` shows the
exact version in use. The full list is in [LICENSE](LICENSE).

## Licence

MIT for the code. Every bundled database keeps its own licence — see
[LICENSE](LICENSE).
