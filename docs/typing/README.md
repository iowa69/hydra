# Organism-specific typing

Hydra types a bacterial assembly twice. The first pass is the same for every
organism: acquired resistance genes, point mutations, virulence, plasmids, MLST,
species. The second pass is organism-specific, and it is what this directory
documents: the schemes a clinical or public-health laboratory actually reports
for that pathogen, applied automatically to the species that was found.

## What is typed, for which organism

| Module | Organism | Schemes reported |
|---|---|---|
| [`saureus`](saureus.md) | *Staphylococcus aureus* | spa type, *agr* group, MLST clonal complex, SCCmec, PVL, *tst*, enterotoxin and haemolysin repertoires, *hlb* status, capsule locus, immune evasion cluster, ACME, β-lactamase locus |
| [`spyogenes`](strep.md) | *Streptococcus pyogenes* | emm type and cluster, virulence/superantigen repertoire, prophage content, regulatory loci, adhesins |
| [`sagalactiae`](strep.md) | *Streptococcus agalactiae* | capsular serotype, alpha-like protein (Alp) type, surface proteins, pilus islands, marker genes, PBP2X substitutions |
| [`spneumoniae`](spneumoniae.md) | *Streptococcus pneumoniae* | capsular serotype, PBP type (1A:2B:2X) with the CDC reference MIC where the triplet is in the table, non-β-lactam resistance determinants, pilus islets |
| [`hinfluenzae`](hinf_nmen.md) | *Haemophilus influenzae* | capsule type (*cap* locus), PBP3 (*ftsI*) substitution pattern, β-lactam interpretation, MLST |
| [`nmeningitidis`](hinf_nmen.md) | *Neisseria meningitidis* | capsule serogroup, fine type, Bexsero antigen sequence type, MLST and clonal complex, resistance alleles |
| [`salmonella`](salmonella.md) | *Salmonella enterica* | serovar from the antigenic formula, subspecies, H2 locus, pathogenicity islands, virulence plasmid, typhoidal markers, MLST/eBurst group |
| [`ecoli`](ecoli.md) | *Escherichia coli* and *Shigella* | phylogroup, *fimH* allele, ST131 subclone, Shiga toxin subtype, pathotype, intimin subtype, *Shigella*/EIEC cluster, virulence repertoire |

Each page states what every call means clinically, the exact rules and
thresholds, the limits of working from a draft assembly, the database provenance
and licence, the validation table, and the bibliography.

## How a module is chosen

`hydra_amr.typing.organisms.registry` decides, per sample:

1. **the species call** — a module declares the species names it applies to, and
   a strain-level name (`Escherichia coli O157:H7 str. Sakai`) still matches the
   species in front of it;
2. **failing that, the MLST scheme** — a sample that typed against `saureus` is
   a *S. aureus* whatever the sketch said.

A module whose database is not installed is **skipped with the reason recorded
against the sample**, never silently and never as a negative result. Nothing is
typed from reads: every scheme here is read off contigs.

A sample selects no module when neither its species nor its MLST scheme is
claimed by one. That is deliberate: *Streptococcus dysgalactiae* is not put
through the *S. pyogenes* module, because none of its schemes was defined in a
group C/G streptococcus — even though the CDC emm database holds the group C/G
`stG`/`stC`/`stL` types, and the *S. pyogenes* module names them correctly when
it meets one in a genome called as *S. pyogenes*. In that case the report heading
says *S. dysgalactiae*, not "emm type", so the call cannot be read as a group A
type.

## Where the evidence goes

Every call carries three things beyond the call itself, and all three reach the
output:

* **evidence** — the alleles, residues, identity and coverage, contig
  coordinates and database version the call was read from (`TypingResult.alleles`);
* **caveats** — what is provisional, what disagreed, what could not be
  determined (`TypingResult.note`). "Not determined" is never a negative result:
  the note says why;
* **loci found / total**, where the scheme counts loci.

They appear in the `*.typing.tsv` table (one row per call), in the `Typing`
sheet of the Excel workbook, and in the *Organism typing* section of the HTML
report. The `*.typing.overview.tsv` table and the `Summary` sheet carry the
headline calls, each with its caveat beside it.

## Installing the data

```
hydra db download organisms          # every module
hydra db download typing-saureus     # one of them
hydra db info typing-saureus         # source, version, licence, redistribution
```

Each module's `fetch()` downloads and builds its own data under
`$HYDRA_DB/organisms/<key>/`, and records the source, version, licence and
redistribution status of every database in the manifest. Hydra does not
redistribute third-party sequence data; the Windows application bundles only the
databases whose terms permit it, and downloads the rest on first launch.

PubMLST data is fetched **anonymously and must never be given an API key**:
anonymous access returns only data submitted on or before 31 December 2024,
which PubMLST's terms permit redistributing without restriction. The visible
cost is that a sequence type assigned after that date reads as "not determined",
with that reason given.

## Running it

Organism typing runs automatically as part of `hydra run`. It is gated by the
existing `--typing` flag, and has its own:

```
hydra run --no-organism-typing  ...     # skip it
hydra run --organism-module saureus ... # only this module
```
