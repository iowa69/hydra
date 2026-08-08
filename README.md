<p align="center">
  <img src="docs/logo.svg" alt="Hydra" width="120" height="120">
</p>

<h1 align="center">Hydra</h1>

<p align="center">
  <strong>Acquired resistance and virulence genes, resistance point mutations,
  heteroresistance from reads, MLST, SCCmec and lineage typing —<br>
  one command, one table, one pass over one set of databases.</strong>
</p>

<p align="center">
  <a href="#install">Install</a> ·
  <a href="#usage">Usage</a> ·
  <a href="#flags">Flags</a> ·
  <a href="#outputs">Outputs</a> ·
  <a href="#how-it-was-measured">Validation</a>
</p>

---

## Performance

Measured on **5375 *Klebsiella* samples**: 667 closed reference genomes, 1279 clinical
isolates, 2933 published genomes whose sequence type was recorded independently, and
496 screened from assembly and reads together. Every one completed.

| | Hydra | Reference tool |
|---|---|---|
| Sequence type vs the recorded answer, 2933 genomes | **2933 / 2933 — 100%** | — |
| Sequence type vs Kleborate 3.2.4 | **2933 / 2933 — 100%** | 2933 / 2933 |
| Sequence type vs `mlst` 2.35.0 | **99.85 – 100%** | 3 genomes it will not type, Hydra types |
| Gene recall vs abricate 1.4.0, nine databases | **116,380 / 116,381 — 100%** | abricate 100% |
| Distinct sequence types called | **539** | — |
| Virulence score vs Kleborate | **99.45%** | reference |
| Colistin, very major errors | **56** | Kleborate 64 |
| Colistin, balanced accuracy | **0.706** | Kleborate 0.666 |
| Point mutations vs AMRFinderPlus 4.2.7 | **556 shared** | 1 unique to Hydra |
| SCCmec vs strains with published types | **4 / 4** | — |
| Heteroresistance, 5% minority allele | **0.0453, p=1.1e-21** | not attempted by any |
| Whole job, per genome | **11.7 s** | 18.1 s for four tools |
| Documented flags behaving correctly | **69 / 69** | — |

<p align="center"><img src="docs/figures/fig-headtohead.svg" alt="Head to head" width="820"></p>

<p align="center"><img src="docs/figures/fig-speed.svg" alt="Runtime" width="820"></p>

<p align="center"><img src="docs/figures/fig-recall.svg" alt="Gene recall" width="820"></p>

## Install

```bash
conda create -n hydra -c conda-forge -c bioconda \
    python=3.11 pip blast minimap2 samtools mash pysam pandas numpy scipy openpyxl
conda activate hydra
pip install https://github.com/iowa69/hydra/archive/refs/tags/v1.4.0.tar.gz

hydra db download     # genes, mutations, every PubMLST scheme, SCCmec, lineage, sketches
hydra db check        # verify everything loads
```

Databases live in `$HYDRA_DB`, default `~/.hydra/db`. Hydra never reaches the network
on its own — `hydra db update` is a command you run, so a database cannot change
underneath a study.

<p align="center"><img src="docs/figures/fig-example.svg" alt="One command, one table" width="820"></p>

## Usage

```bash
hydra run -a isolate.fasta -o results/                       # one isolate, everything on
hydra run assemblies/ -o results/ --preset surveillance      # a directory, with heatmaps
hydra run -1 s_R1.fq.gz -2 s_R2.fq.gz -o results/            # reads only, no assembly
hydra run -a s.fasta -1 s_R1.fq.gz -2 s_R2.fq.gz -o results/ # both: adds allele fractions
hydra run --input-list samples.tsv -o results/ --preset deep -t 24   # a thousand samples
hydra run -a mrsa.fasta -o results/                          # a staphylococcus: adds SCCmec
hydra screen -d card assemblies/*.fasta --stdout             # one database, flat gene table
```

`--input-list` is `sample`, `assembly`, `R1`, `R2`, tab-separated. Leave a field empty
when a sample has no assembly or no reads:

```
KP001   KP001.fasta   KP001_R1.fq.gz   KP001_R2.fq.gz
KP002   KP002.fasta
KP003                 KP003_R1.fq.gz   KP003_R2.fq.gz
```

Ask it what it can do — nothing needs to be told the organism or the databases:

```bash
hydra presets                  # every preset and what it turns on
hydra run --list-databases     # what is installed, how big, which version
hydra run --list-organisms     # every -O value, with its mutation counts
hydra db info card             # source, licence and citation for one database
```

### Presets

A preset only sets defaults; anything passed explicitly wins.

| Preset | Turns on |
|---|---|
| `fast` | NCBI genes only, no translated search or typing |
| `standard` | default: NCBI + VFDB, translated search, mutations, MLST, typing |
| `deep` | every database, `--plus` elements, everything on |
| `surveillance` | multi-sample: all matrices, class recap, HTML heatmaps |
| `amr` | NCBI, CARD, ResFinder + point mutations |
| `virulence` | VFDB and the lineage schemes |
| `plasmid` | replicon typing and the AMR genes alongside |
| `linezolid` | 23S heteroresistance from reads, plus *cfr*/*optrA*/*poxtA* |
| `enterobacterales` | *Klebsiella*/*E. coli*: AMR, virulence loci, lineage scores |
| `gram-positive` | staphylococci, enterococci, streptococci |
| `genes` | permissive thresholds, flat one-row-per-gene layout |
| `elements` | translated search with mutations, as a typed element table |

## Flags

### Inputs

| Flag | Does |
|---|---|
| `INPUT...` | assemblies, FASTQ files or directories to scan (recursive; FASTA/FASTQ sorted automatically) |
| `-a`, `--assembly FASTA` | assembly FASTA (repeatable) |
| `-1`, `--r1 FASTQ` | forward reads (repeatable, pairs with `--r2`) |
| `-2`, `--r2 FASTQ` | reverse reads (repeatable) |
| `--reads FASTQ` | reads to auto-pair by filename (repeatable) |
| `--input-list TSV` | sample sheet: sample, assembly, R1, R2 |
| `--name NAME` | override sample names, in input order (repeatable) |

### Databases

| Flag | Does |
|---|---|
| `-d`, `--db NAME` | database or group; repeatable or comma-separated. Groups: `all`, `standard`, `amr`, `virulence`, `nucl`, `core` |
| `--list-databases` | print the installed databases and exit |

### Analysis

| Flag | Does |
|---|---|
| `--preset NAME` | option bundle (see above) |
| `-O`, `--organism NAME` | organism for point mutations, e.g. `Staphylococcus_aureus` (default: detected from MLST/Mash) |
| `--list-organisms` | print organisms with point-mutation support and exit |
| `--auto-organism` / `--no-auto-organism` | detect the organism automatically (default) / only use `-O` |
| `--plus` / `--no-plus` | also report stress-response and virulence elements / acquired resistance only |
| `--mlst` / `--no-mlst` | run MLST (default) / skip it |
| `--scheme NAME` | force a PubMLST scheme instead of choosing automatically |
| `--typing` / `--no-typing` | run lineage and SCCmec typing (default) / skip |
| `--protein` / `--no-protein` | run the translated AMR search (default) / skip |
| `--point-mutations` / `--no-point-mutations` | call resistance point mutations (default) / skip |
| `--heteroresistance` / `--no-heteroresistance` | measure allele fractions from reads (default with reads) / skip |
| `--reads-mlst` / `--no-reads-mlst` | type from reads when no assembly produced an ST (default) / never |
| `--reads-variants` / `--no-reads-variants` | report read-vs-reference differences per gene (default with reads) / skip |
| `--report-synonymous` | also report read-derived variants that do not change the protein |
| `--assemble` | assemble read-only samples first (needs `skesa` or `spades.py`) |

### Thresholds

| Flag | Default | Does |
|---|---|---|
| `--min-identity PCT` | 80 | minimum % identity for nucleotide hits |
| `--min-coverage PCT` | 60 | minimum % of the reference covered |
| `--protein-min-identity PCT` | 90 | minimum % identity for translated hits |
| `--protein-min-coverage PCT` | 90 | coverage above which a translated hit is complete |
| `--min-contig-length BP` | 0 | ignore contigs shorter than this |
| `--report-overlaps` | off | report every overlapping hit, including one inside another, instead of the best per locus |

### Reads and heteroresistance

| Flag | Default | Does |
|---|---|---|
| `--min-depth N` | 5 | minimum read depth for a call |
| `--min-gene-breadth PCT` | 80 | minimum % of a gene covered by reads |
| `--min-allele-fraction FRAC` | 0.02 | lowest minority allele fraction to report |
| `--fixed-allele-fraction FRAC` | 0.90 | at or above this a mutation is fixed, not heteroresistant |
| `--min-allele-reads N` | 3 | minimum reads supporting a minority allele |
| `--min-base-quality Q` | 13 | minimum base quality in the pileup |
| `--report-absent-sites` | off | also report catalogued sites where no mutation was found |

### Output

| Flag | Default | Does |
|---|---|---|
| `-o`, `--outdir DIR` | — | write results here (created if missing) |
| `--prefix NAME` | `hydra` | basename for output files |
| `-f`, `--format FMT` | `tsv,html` | `tsv`, `csv`, `json`, `html`, `xlsx`, `genes`, `elements`; repeatable or comma-separated |
| `--stdout` | off | write the long table to stdout instead of files |
| `--cell MODE` | `binary` | matrix values: `binary`, `identity`, `coverage`, `count`, `genes`, `depth`, `fraction`, `symbol` |
| `--rows FIELD` | `sample` | matrix rows: `sample`, `gene`, `class`, `subclass`, `database`, `element_type`, `product` |
| `--columns FIELD` | `gene` | matrix columns: same choices as `--rows` |
| `--element-types LIST` | all | restrict matrix and heatmaps to `AMR`, `VIRULENCE`, `STRESS`, `PLASMID`, `SEROTYPE` |
| `--title TEXT` | — | title for the HTML report |

### Runtime

| Flag | Default | Does |
|---|---|---|
| `-t`, `--threads N` | all cores | CPU threads to use |
| `-j`, `--jobs N` | auto | read-mapping samples processed concurrently (assembly screening is batched and always uses all threads) |
| `--db-dir DIR` | `$HYDRA_DB` | database directory |
| `--tmpdir DIR` | system temp | directory for intermediate files |
| `--keep-temp` | off | keep intermediate BLAST/BAM files for debugging |
| `-v`, `--verbose` | off | verbose logging (repeat for more) |
| `-q`, `--quiet` | off | only warnings and errors |

### Commands

```
hydra run       full analysis of assemblies and/or reads
hydra screen    acquired-gene screening only, as a flat gene table
hydra db        list | import | download | update | bundle | info | check | remove
hydra presets   list the available presets
```

```bash
hydra db download            # fetch everything with a stable source
hydra db import              # convert copies already in conda environments
hydra db update --dry-run    # what would change, and when each was installed
hydra db update              # refresh from upstream; never automatic
hydra db bundle -o db.tar.gz # pack for an offline machine
hydra db check               # verify every installed database still loads
```

## Outputs

One run writes a table per question. Every block below is real output.

### `hydra.summary.tsv` — one row per isolate

Species and the evidence for it, sequence type with its allele profile, capsule
marker, counts by category, the scores, and assembly QC. 48 columns.

```
sample                ERR11578019
species               Klebsiella pneumoniae        species_confidence   strong
species_evidence      MLST scheme klebsiella (7/7 exact loci); Mash Klebsiella pneumoniae
mlst_scheme           klebsiella                   ST                   3419
mlst_profile          gapA(3) infB(1) mdh(1) pgi(26) phoE(4) rpoB(4) tonB(39)
wzi                   178                          organism_db          Klebsiella_pneumoniae
amr_genes             42                           amr_classes          8
virulence_genes       82                           stress_genes         32
plasmid_replicons     3                            point_mutations      0
resistance_score      1                            virulence_score      1
has_esbl              True                         has_carbapenemase    False
qc_contigs            3    qc_n50 5228768    qc_total_length 5457457    qc_gc 57.29
```

### `hydra.tsv` — one row per element

Acquired genes, point mutations, disrupted genes, plasmid replicons, virulence,
stress and serotype markers all share one schema, so nothing needs a second parser.

```
sample        database       element_type  subtype     gene         class                %cov    %id
DRR199175     card           AMR           AMR         ArnT                              100.0   99.52
ERR1015312    protein        AMR           POINT       ramR         TETRACYCLINE         100.0   98.97
ERR10447218   protein        AMR           DISRUPTION  ramR         TIGECYCLINE           68.4   100.0
DRR199175     plasmidfinder  PLASMID       REPLICON    Col440I_1                         100.0   96.49
DRR199175     ecoli_vf       VIRULENCE     VIRULENCE   ECS88_3547                        100.0   90.06
ERR1015312    protein        STRESS        BIOCIDE     qacEdelta1   QUATERNARY AMMONIUM  100.0   100.0
DRR199175     protein        STRESS        METAL       fieF                              100.0   100.0
ERR10447212   ecoh           SEROTYPE      SEROTYPE    wzm-O9                            100.0   96.31
```

Columns: `sample database element_type element_subtype gene accession product class
subclass sequence start end strand coverage coverage_pct identity_pct gaps depth
allele_fraction method resolution primary note`.

`primary` marks the best hit per genomic locus, so a gene found by several databases
counts once. `--report-overlaps` reports the rest as well.

### `hydra.matrix.tsv` — presence/absence, samples × genes

101 samples × 407 genes in the run below. `--cell` changes what fills it, `--rows`
and `--columns` change what it pivots on.

```
sample        (Bla)ampH  APECO1_3698  ArnT  BAER      # --cell binary (default)
DRR199175     1          0            1     1
ERR1015312    1          1            1     1

sample        MerP_Gneg  aac(3)-IId   aadA2           # --cell identity
DRR199175     87.73      100.0        0.0
```

### `hydra.classes.tsv` — elements per drug class

```
sample        AMINOGLYCOSIDE  BETA-LACTAM  BIOCIDE  COLISTIN  EFFLUX  FOSFOMYCIN
DRR199175     3               5            6        0         3       1
ERR1015312    3               4            6        0         3       1
```

### `hydra.mlst.tsv` — scheme, ST and the allele profile

```
sample        scheme      ST   source     loci_exact  loci_total  gapA  infB  mdh
DRR199175     klebsiella  101  assembly   7           7           2     6     1
ERR1015312    klebsiella  15   assembly   7           7           1     1     1
```

### `hydra.typing.tsv` — lineage, capsule, SCCmec and scores

22 columns: `AbST` aerobactin, `cbST` colibactin, `RmST`/`rmpA2` hypermucoidy, `SmST`
salmochelin, `ybST` yersiniabactin, `wzi` capsule, `SCCmec` for staphylococci, each
with its lineage, plus the scores.

```
sample        ybST  ybST_lineage      AbST  AbST_lineage  RmST_lineage    wzi   virulence_score
ERR11578027   41    ybt 12; ICEKp10   3     iuc 2         rmp 2; KpVP-2   257   5
```

### `hydra.html` — the same content as a report

Self-contained, no external assets, 7 sections: sample overview, resistance heatmap,
virulence heatmap, MLST, lineage typing and scores, all detected elements. Heatmaps
are clustered. `--title` names it.

### `hydra.json` — everything, with its provenance

`hydra_version`, `command`, `databases`, `parameters` and then one record per sample
carrying `species`, `mlst`, `typing`, `scores`, `qc`, `hits` and `warnings` — so a
result can be traced back to the version and command that produced it.

### `-f genes` and `-f elements` — drop-in layouts

`genes` writes abricate's exact columns (`#FILE SEQUENCE START END STRAND GENE
COVERAGE COVERAGE_MAP GAPS %COVERAGE %IDENTITY DATABASE ACCESSION PRODUCT
RESISTANCE`). `elements` writes AMRFinderPlus's. Both drop into scripts that already
read one.

## How it was measured

Every comparison gave both tools the same genomes, the same database, the same
thresholds and the same overlap reporting, and gene calls were matched on locus
rather than on label — the same gene is held under different synonyms across
redundant databases.

Full detail, every comparator and every sample screened:
**[docs/VALIDATION.md](docs/VALIDATION.md)** ·
**[docs/validation_klebsiella_samples.tsv](docs/validation_klebsiella_samples.tsv)**

## Licence

MIT. The reference databases carry their own licences and citations —
`hydra db download --list` prints each one, `hydra db info NAME` prints them
individually.
