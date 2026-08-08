<p align="center">
  <img src="docs/logo.svg" alt="Hydra" width="120" height="120">
</p>

<h1 align="center">Hydra</h1>

<p align="center">
  <strong>Acquired resistance and virulence genes, resistance point mutations,
  heteroresistance from reads, MLST and lineage typing —<br>
  one command, one table, one pass over one set of databases.</strong>
</p>

<p align="center">
  <a href="#install">Install</a> ·
  <a href="#how-you-would-run-it">Examples</a> ·
  <a href="#what-it-reports">Output</a> ·
  <a href="#validation">Validation</a> ·
  <a href="#commands">Commands</a>
</p>

---

Characterising a bacterial isolate normally means four or five tools, each with its
own database, thresholds and output layout, and a script to glue the answers
together. Hydra does the whole job at once — and finds resistance mutations present
in only some copies of a gene, which an assembly consensus cannot show at all.

## Install

```bash
conda create -n hydra -c conda-forge -c bioconda \
    python=3.11 pip blast minimap2 samtools mash pysam pandas numpy scipy openpyxl
conda activate hydra
pip install https://github.com/iowa69/hydra/archive/refs/tags/v1.3.0.tar.gz

hydra db download     # genes, mutations, all PubMLST schemes, lineage, sketches
```

`hydra db check` reports anything missing. Databases live in `$HYDRA_DB`
(default `~/.hydra/db`).

## How you would run it

**One isolate, everything on.**

```bash
hydra run -a isolate.fasta -o results/
```

**A directory of assemblies, with matrices and clustered heatmaps.** Directories are
scanned recursively and FASTA/FASTQ sorted out automatically.

```bash
hydra run assemblies/ -o results/ --preset surveillance
```

**Raw reads only.** No assembly needed — genes are called by depth and breadth, and
MLST from reads mapped to the scheme's loci.

```bash
hydra run -1 sample_R1.fq.gz -2 sample_R2.fq.gz -o results/
```

**Assembly and reads together.** The consensus calls, plus allele fractions at every
catalogued position — this is what finds a mutation carried by some copies of a
multi-copy locus.

```bash
hydra run -a sample.fasta -1 sample_R1.fq.gz -2 sample_R2.fq.gz -o results/
```

**A staphylococcus.** Species, ST, SCCmec type and resistance in one pass — nothing
needs to be told it is a staphylococcus.

```bash
hydra run -a mrsa.fasta -o results/          # SCCmec type II(2A), ST5, ...
```

**Linezolid heteroresistance from 23S allele fractions.**

```bash
hydra run -1 s_R1.fq.gz -2 s_R2.fq.gz --organism Staphylococcus_aureus \
    --preset linezolid -o results/
```

**A thousand samples from a sample sheet.** Columns are `sample`, `assembly`, `R1`,
`R2`; leave a field empty when a sample has no assembly or no reads.

```bash
hydra run --input-list samples.tsv -o results/ --preset deep -t 24
```

**One database, as a flat gene table.** Writes abricate's column layout, so it drops
into scripts that already read one.

```bash
hydra screen -d card assemblies/*.fasta --stdout
```

**Ask it what it can do.** Nothing needs to be told the organism, the databases or
the catalogued mutations.

```bash
hydra presets                  # every preset and what it turns on
hydra run --list-databases     # what is installed, how big, which version
hydra run --list-organisms     # every -O value, with its mutation counts
```

### Presets

| Preset | For |
|---|---|
| `fast` | quickest useful screen: NCBI genes only |
| `standard` | the default: genes, translated search, mutations, MLST, typing |
| `deep` | every database, `--plus` elements, everything on |
| `surveillance` | multi-sample: all matrices, class recap, HTML heatmaps |
| `amr` · `virulence` · `plasmid` | one job each |
| `linezolid` | 23S heteroresistance from reads, plus *cfr*/*optrA*/*poxtA* |
| `enterobacterales` · `gram-positive` | organism-group defaults |
| `genes` · `elements` | flat gene table · typed element table |

A preset only sets defaults; anything passed explicitly wins.

## What it reports

One long table plus the artefacts a batch needs — presence/absence matrix, drug-class
recap, MLST, lineage typing, a self-contained HTML report, and JSON carrying the
version, command and databases used.

| | |
|---|---|
| **Acquired genes** | NCBI, CARD, ResFinder, ARG-ANNOT, MEGARes, VFDB, PlasmidFinder, EcOH, *E. coli* VF — screened together, drug-class annotation carried across |
| **Point mutations** | *gyrA*, *parC*, *rpoB*, *pmrB*, 16S/23S rRNA, promoters such as `pbp4_T-266A`, mosaic-PBP calls |
| **Heteroresistance** | allele fractions from reads at every catalogued position, with the number of rRNA operons carrying the mutation |
| **MLST** | every installed PubMLST scheme searched at once, scheme chosen automatically |
| **Species** | Mash sketches and the MLST result combined; separates *K. variicola* and *K. quasipneumoniae* from *K. pneumoniae* |
| **Lineage** | yersiniabactin, colibactin, aerobactin, salmochelin, *rmp*, *wzi*, *E. coli* serotyping, pathovar markers |
| **SCCmec** | cassette type for staphylococci, from whole-element references (types I–XIII) |
| **Broken genes** | *mgrB*, *ompK35/36*, *cirA*, *ramR* — where losing a gene is the resistance mechanism, not gaining one |
| **Scores** | 0–5 virulence and 0–3 resistance, generalised to any Gram-negative |

## Validation

Measured on **5375 *Klebsiella* samples** — 667 closed reference genomes, 1279 clinical
isolates, 2933 published genomes whose sequence type was recorded independently, and
496 isolates screened from assembly and reads together. Every one completed.

Every comparison below gave both tools the same genomes, the same database, the same
thresholds and the same overlap reporting.

<p align="center"><img src="docs/figures/fig-headtohead.svg" alt="Head to head against the reference tools" width="820"></p>

<p align="center"><img src="docs/figures/fig-recall.svg" alt="Gene recall against abricate, per database" width="820"></p>

Across nine databases and 667 genomes, **123 gene instances were found by abricate and
not by Hydra — 0.18 per genome.** Five databases match exactly.

<p align="center"><img src="docs/figures/fig-speed.svg" alt="Runtime against four reference tools" width="820"></p>

### Sequence typing

**2933 of 2933** published genomes match the sequence type recorded for them, and
**2933 of 2933** match Kleborate 3.2.4. Against `mlst` 2.35.0 the counts are 666/667,
1279/1279 and 2931/2933 — and all three differences are genomes `mlst` declined to
type because a locus appeared twice (`gapA(51,51)`, `rpoB(4,4)`) where both copies are
the same allele and the profile is unambiguous. Hydra typed each correctly, and the
recorded sequence type agrees.

### Colistin, where losing a gene is the mechanism

Colistin resistance in *K. pneumoniae* is usually *mgrB* — a 47-residue repressor —
broken by an insertion sequence. There is no acquired gene to find and no catalogued
substitution, so Hydra detects the truncation itself. On the 766 isolates with a
measured EUCAST result:

| | sensitivity | specificity | balanced accuracy | very major errors |
|---|---|---|---|---|
| **Hydra** | **0.440** | 0.973 | **0.706** | **56** |
| Kleborate 3.2.4 | 0.360 | 0.973 | 0.666 | 64 |

### Typing that depends on the organism

Some answers are not allele profiles. SCCmec is a cassette that is either there or
not, so the measurement is coverage of a whole reference element. On strains with
published types — and none of them in the reference database, which would make the
test circular:

| Genome | Published | Hydra |
|---|---|---|
| N315 | SCCmec II | **II(2A)** |
| USA300_FPR3757 | SCCmec IV | **IV(2B)** |
| NCTC 8325 | none, MSSA | **none** |
| Newman | none, MSSA | **none** |

A methicillin-susceptible genome still covers 21% of the nearest reference — the
*orfX* flank every *S. aureus* carries — against 98–100% for a real cassette, so the
floor sits in the middle of that gap. The scheme runs only for staphylococci.

### Genes whose loss is the mechanism

Most resistance is a gene arriving; some is a gene breaking, and a screen that
reports what it finds intact reports nothing. Hydra detects the disruption itself for
*mgrB* (colistin), *ompK35* and *ompK36* (carbapenem, porin loss), *cirA*
(cefiderocol) and *ramR* (tigecycline).

This finds structural disruption — a frameshift or an inserted element breaks the
reading frame and truncates the alignment. On 300 closed genomes it recovers 28 of
the 44 events AMRFinderPlus reports, where the previous release found none of them.
The 16 it does not are clean nonsense substitutions: a single premature stop leaves
the DNA aligning full length, so coverage cannot see it and the catalogued-mutation
path is what finds those.

### Heteroresistance

Validated against a synthetic control — reads simulated from a 23S reference with a
known fraction carrying `23S_G2577T`:

| Simulated | Measured | Called |
|---|---|---|
| 20% (1 of 5 operons) | 0.2004, 93/464 reads, p=2e-152 | heteroresistant, ~1.0/5 operons |
| 5% | 0.0453, 21/464 reads, p=1.1e-21 | heteroresistant |

Full detail, every comparator and every sample screened:
**[docs/VALIDATION.md](docs/VALIDATION.md)** ·
**[docs/validation_klebsiella_samples.tsv](docs/validation_klebsiella_samples.tsv)**

## Commands

```
hydra run       full analysis of assemblies and/or reads
hydra screen    acquired-gene screening only, as a flat gene table
hydra db        list | import | download | bundle | info | check | remove
hydra presets   list the available presets
```

<details>
<summary><code>hydra run</code> options</summary>

**Inputs** — `INPUT...`, `-a/--assembly`, `-1/--r1`, `-2/--r2`, `--reads`,
`--input-list`, `--name`

**Databases** — `-d/--db` (names, or the groups `all`, `standard`, `amr`,
`virulence`, `nucl`, `core`), `--list-databases`

**Analysis** — `--preset`, `-O/--organism`, `--list-organisms`,
`--auto-organism/--no-auto-organism`, `--plus/--no-plus`, `--mlst/--no-mlst`,
`--scheme`, `--typing/--no-typing`, `--protein/--no-protein`,
`--point-mutations/--no-point-mutations`, `--heteroresistance/--no-heteroresistance`,
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

## Licence

MIT. The reference databases carry their own licences and citations — `hydra db
download --list` prints each one, and `hydra db info NAME` prints them individually.
