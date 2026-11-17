# *Staphylococcus aureus* typing

Module: `src/hydra_amr/typing/organisms/saureus.py` — `KEY = "saureus"`.
Triggered by species `Staphylococcus aureus` or PubMLST scheme `saureus`.
Assemblies in; BLAST+ and pure Python only.

SCCmec is **not** typed here. Hydra's `hydra_amr.typing.cassette.CassetteTyper`
already types the whole cassette, and duplicating it would produce two calls
that can disagree. See "Interaction with the SCCmec typer" below.

---

## 1. What is called, and what each call means

| Scheme | Call | Clinical meaning |
|---|---|---|
| `spa` | a Ridom t-number, `new succession`, or `not determined` | The protein A repeat type: the short-turnaround MRSA genotype most European reference laboratories report, and the unit outbreak reports are written in. |
| `agr` | `I`–`IV` | Which autoinducing peptide the isolate makes, and which groups it cross-inhibits. Lineage-linked; reported alongside the clone in most virulence work. |
| `blaZ` | `present` / `absent`, with the allotype (`A`-`F`) and the lesion state (`intact`, `truncated`, `internal stop`, `frameshifted`, `interrupted`) | The staphylococcal penicillinase and its kinetic class. Allotypes A and C are those associated with a strong cefazolin inoculum effect, which matters if cefazolin is being considered for a deep-seated infection. |
| `bla_operon` | the state of *blaR1* and *blaI*, with the predicted phenotype | **Penicillin susceptibility does not follow from *blaZ* alone.** BlaR1 must sense the β-lactam and cleave the BlaI repressor for *blaZ* to be induced, so a *blaZ*-positive isolate with a dead *blaR1* and an intact *blaI* can test penicillin-susceptible. This is the call that stops a bare gene hit being read as resistance. |
| `pvl` | `positive` / `negative` / `not determined` | Panton-Valentine leukocidin. Marker of the community-associated lineages; associated with recurrent skin and soft-tissue abscesses and with necrotising pneumonia. |
| `tst` | `positive` / `negative` | TSST-1, on a SaPI island. The determinant to report in a patient with shock and a rash. |
| `enterotoxins` | the genes found, plus the egc count | Food-poisoning risk (the classical *sea*–*see*) and superantigen repertoire. The *egc* cluster is reported as a cluster because it travels as one unit on νSaβ and is weakly expressed. |
| `haemolysins` | the cytolysin genes found | *hla*, LukAB/LukGH and LukED are core-genome; absence usually means a fragmented assembly, which the note says. |
| `hlb` | `intact` / `disrupted by prophage integration` | β-haemolysin. Reported separately from the other haemolysins because a prophage normally sits inside it (below). |
| `capsule` | `5` / `8` / `not determined` | Capsular serotype from the type-specific *capHIJ* genes. |
| `iec` | `type A`–`type G` / `absent` | Immune-evasion cluster, on the same prophage as the *hlb* insertion. A stable lineage and host-adaptation marker. |
| `acme` | `arc cluster present` / `absent` | The arginine catabolic mobile element, next to SCCmec in USA300. |
| `clonal_complex` | `CCn`, `no established complex`, `not determined` | How *S. aureus* clones are named in surveillance. CC5, CC8, CC22, CC30 and CC45 account for most hospital MRSA. |

Every call carries its evidence in `TypingResult.alleles` (per-gene identity,
reference coverage as a span over the reference length, accession and
coordinates) and its caveats in `note`. No scheme returns a bare negative: a
negative always states what was looked for and why it was not found.

---

## 2. Method, rules and thresholds

### 2.1 One BLAST pass, measured on the reference axis

`fetch()` builds a 58-sequence nucleotide panel, each entry pinned to an NCBI
accession and a coordinate range (see §4). `run()` makes a BLAST database of
the assembly and queries the panel against it, so hit coordinates land on the
assembly's own contigs — which the lesion calls need.

Hydra's `merge_hsps` merges on the *subject* axis, which is right when the
reference is the subject. Here the reference is the query, so `_merge_on_reference`
does the same job over the same interval helpers (`merge_intervals`,
`interval_length`, `overlap_length`): HSPs are grouped by (marker, contig,
strand), an HSP that re-covers reference ground already claimed is dropped so
two gene copies cannot sum past 100% coverage, and identity is weighted by the
reference each HSP *newly* covers so overlapping HSPs do not double-count.

Thresholds, and why:

| Constant | Value | Reason |
|---|---|---|
| `MIN_IDENTITY` | 90% | *lukF-PV* aligns to the *lukF*/*hlgB* family at **81%** in every PVL-negative genome on the panel; *seb* cross-hits at 75–76%, *sea* at 84%, *sep* at 84–85%. A coverage-only rule would call every genome PVL-positive. |
| `MIN_COVERAGE` | 90% | Present means nearly the whole reference. |
| `PARTIAL_COVERAGE` | 40% | Between 40% and 90% a gene is reported **partial**, not silently present or absent. |
| `INSERT_GAP` | 1000 bp | Above this, two alignment pieces of one reference gene mean an insertion, not an indel. |
| `AGR_MIN_IDENTITY` / `AGR_MIN_COVERAGE` | 92% / 80% | Between-group identity over the *agrBDC* window is 73–81%; within-group 98.5–100%. |

Marker choice is pinned to **coordinates, not gene names**, because RefSeq
names *hla* `hyl`, leaves *chp* unnamed, and labels the capsule *H*/*I*/*J*
genes `cap8H`/`cap8I`/`cap8J` **in serotype 5 genomes as well as serotype 8
ones**. A name-driven build would mistype the capsule.

### 2.2 spa

1. The two Ridom protocol primers are located with `blastn -task blastn-short`:
   `spa-1113f` = `TAAAGACGATCCTTCGGTGAGC` and `spa-1514r` = `CAGCAGTAGTGCCGTTTGCTT`
   (Ridom *DNA Sequencing of the spa Gene* v1.1; numbered against GenBank
   J01786 at 1092–1113 and 1534–1514). A primer counts as an anchor only if
   ≥ its length − 2 bases matched at ≥90% identity.
2. **Both anchors must be on one contig in the right relative orientation.**
   One anchor alone, or two that never pair, is reported as `not determined`
   with the reason — most often a contig break. Minus-strand products are
   reverse-complemented before tiling.
3. The excised region must be 60–1500 bp.
4. The region is tiled **longest-match-first**. The Ridom repeats are 18, 21,
   24, 27 and 30 bp and the shorter ones are prefixes of the longer ones, so
   trying lengths in arbitrary order mis-tiles the array.
5. The hyphen-joined succession is looked up in `spatypes.txt`.

Refusals and caveats, all of which fire in practice:

- **Gaps inside the array**, or a region reaching within 50 bp of a contig end
  → `not determined`, with the matching type named in the note but not
  reported as the call.
- **Succession not in the table** → `new succession`, with the succession, and
  advice to submit it to spaServer.
- **Fewer than 3 repeats** → the type is reported with a warning: the Ridom
  table averages 9.5 repeats, so an array this short usually means a deletion.
- **The *spa* gene itself not intact** → the type is reported with a warning
  that protein A is not produced. This fires on RF122, whose *spa* is a
  published pseudogene (Herron-Olson 2007).

Two parsing traps in the Ridom files, both handled and both unit-tested: the
files are **CRLF** with blank lines between FASTA records, and the type table
ends with a `NT,` sentinel whose succession is **empty** — keeping it would make
`types[""] == "NT"`, so every failed extraction would silently type as `NT`.
Repeat codes in the type table are the FASTA header minus its leading `r`.

### 2.3 agr

A single window per group — *agrB* start to *agrC* start + 650 bp, covering
*agrB*, *agrD* and the hypervariable 5′ half of *agrC* — is taken from a genome
with a published group, and the group is the window with the best
**coverage × identity**, not the best identity.

Ranking on identity alone is wrong, and MRSA252 shows why: a 107 bp stretch of
it matches the group II window at 99.1%, while the group III window matches
over its whole 1382 bp at 98.6%. Identity alone calls group II; coverage ×
identity calls group III, which is what CC30 is. This is the same failure
Hydra's own `MergedHit.match_score` docstring describes, and the rule is the
same one.

The window for groups I–III comes from a complete genome with a published
group (Tan 2022). Group IV has **no** complete genome with a published
assignment, so its window is taken from the *agr* locus deposited with the
paper that defined the group: `AF288215`, the exfoliatin-producing strain
RN5881 (Jarraud 2000), in which *agrB* is 718..1281 and *agrC* starts at 1505.
Adding it did not disturb any of the ten panel calls.

The region used — the 3′ end of *agrB*, all of *agrD* and the 5′ end of *agrC*
— is the conventional typing region: Jarraud 2002's Pan-1/Pan-2 primers amplify
"a 1,234-bp *agr* fragment encompassing the 3′ end of *agrB*, all of *agrD*, and
the 5′ end of *agrC*". **AgrA is not used**, because it is conserved across the
groups; the variability is in AgrB, AgrD and the N-terminal half of AgrC.

### 2.4 The bla operon: allotype, lesions and the phenotype

#### The blaZ allotype

The allotypes were first defined **kinetically**, by β-lactamase
substrate-hydrolysis profile (Voladri 1996); the residue pair that
distinguishes them is tabulated against **Ambler's class A standard numbering**
(Ambler 1991):

| allotype | Ambler 128 | Ambler 216 |
|---|---|---|
| A | Thr | Ser |
| B | Lys | Asn |
| C | Thr | Asn |
| D | Ala | Ser |
| E | Leu | Ser |
| F | Thr | Thr |

**Ambler numbering is not a fixed offset from the BlaZ sequence.** Measured
against the 281-residue preproenzyme (signal peptide 1–24 included, UniProt
P00807), the offset is +7 at the catalytic serine (Ambler 70 → 63) but **+9**
at the SDN loop (130 → 121), at the K[ST]G triad (234 → 236 → 225–227) and at
both typing positions. So Ambler 128 and 216 are hard-coded as sequential
residues **119 and 207**, the mapping Carvajal 2020 states in its Methods
("119 and 207 based on our numbering that includes the signal peptide").

Before any residue is read, the translated protein must be 281 residues and
must carry the three class A motifs at their expected places — catalytic Ser at
119−56 = 63, `SDN` at 121–123, `K[ST]G` at 225–227. If a motif is missing the
frame or the numbering cannot be trusted and **no allotype is assigned**. A
residue pair outside the table is reported as not one of the described
combinations, not forced to the nearest letter.

Note that **Olsen 2006 does not define allotypes A–D** — it reports 69 sequence
types translating to 11 BlaZ protein types — so it is not cited for this rule.

#### Lesion-level description

For *blaZ*, *blaR1* and *blaI* the module does not report presence. It reports,
per gene, one of:

| State | How it is decided |
|---|---|
| `absent` | reference coverage < 40% |
| `not determined` | identity < 90% — a related gene, not this one |
| `truncated by contig break` | coverage < 90% **and** the alignment runs into a contig end (within 10 bp) |
| `truncated` | coverage < 90% and the contig continues past the alignment, so the truncation is real |
| `interrupted` | full coverage, but the reference aligns in pieces > 1 kb apart on the contig |
| `frameshifted` | full coverage, but the contig span is not a multiple of three |
| `internal stop` | full coverage, contiguous, but the reading frame stops before the expected residue |
| `intact` | full coverage, contiguous, frame runs to length |

The distinction between `truncated` and `truncated by contig break` is the one
that matters: both lose sequence, only one is evidence about the gene.

The predicted phenotype follows from the combination:

| *blaZ* | *blaR1* | *blaI* | Reported meaning |
|---|---|---|---|
| absent | — | — | penicillin expected susceptible; confirm with a zone-edge or nitrocefin test |
| intact | intact | intact | inducible operon; penicillin expected resistant |
| intact | intact | lost/disrupted | repressor gone → constitutive penicillinase; resistant, zone edge may be sharp |
| intact | lost/disrupted | intact | **cannot be induced — this is the configuration behind *blaZ*-positive, penicillin-susceptible isolates.** Confirm with a zone-edge test |
| intact | absent | absent | constitutive; resistant |

A regulator cut by a contig break adds an explicit "needs a better assembly
before it is acted on".

### 2.5 hlb and the Sa3int prophage

`hlb` covers 96–100% of its reference in nine of the ten validation genomes,
*including* the eight where the gene is inactivated. Coverage cannot see the
problem, because the prophage inserts **inside** the gene and both flanks still
align. The insertion itself is the signal:

```
USA300     hlb q    1- 183  ->  contig 2084489-2084671
           hlb q  170- 993  ->  contig 2127707-2128530     43036 bp apart
COL        hlb q    1- 993  ->  contig 2063726-2064718     one piece
```

The split is at the same place in every disrupted genome — reference base
170–183, the 14 bp attachment-site core `TGTATCCAAACTGG` (*hlb* codons 57–61)
duplicated at both prophage ends on integration — and the two pieces sit
42.6–44.6 kb apart. That insertion point was derived here by alignment and then
confirmed independently by a separate derivation from the same genomes. So: **two alignment pieces, ≥90%
reference coverage, contig separation > 1 kb → disrupted by prophage
integration**, with the separation and the break position in the note. One
contiguous piece at full coverage → intact.

The `iec` call cross-checks this and says so when the two disagree, since the
cluster normally arrives on the prophage that does the disrupting.

### 2.6 Capsule

Serotypes 5 and 8 share most of the *cap* locus; the discriminator is the
type-specific *capH*, *capI*, *capJ* block. The module carries both allele sets
and requires **all three** genes of one set complete with the other set
incomplete. On the panel the two sets are mutually exclusive at 99–100% versus
no hit.

Reported explicitly: an intact locus is **not** an expressed capsule. Several
lineages, USA300 among them, carry *cap5* and are phenotypically
non-encapsulated. This module does not assess expression and says so.

### 2.7 Immune evasion cluster

`scn` must be present for any type to be assigned. The gene combination is then
matched against van Wamel 2006:

| Type | *sea* | *sep* | *sak* | *chp* | *scn* |
|---|---|---|---|---|---|
| A | ● | | ● | ● | ● |
| B | | | ● | ● | ● |
| C | | | | ● | ● |
| D | ● | | ● | | ● |
| E | | | ● | | ● |
| F | | ● | ● | ● | ● |
| G | | ● | ● | | ● |

A combination that matches none is `not determined`, not forced into the
nearest type.

### 2.8 ACME

*S. aureus* carries an arginine deiminase operon in its **core genome**, so
detecting `arcA` is not detecting ACME. The ACME allele
(`NC_007793.1:73113-74348`, in USA300's ACME next to SCCmec) and the core-genome
allele (`NC_007793.1:2781714-2782949`) are **75.9% identical over their full
1236 bp**, measured here. A 90% identity floor separates them cleanly: on the
panel only USA300 matches the ACME allele (100%), and the other nine match it
at 75% — i.e. not at all. The ACME **type** (I/II/III) is not assigned, because
that needs the *opp3* cluster, which this release does not type.

### 2.9 Clonal complex

PubMLST's *S. aureus* scheme 1 carries **no `clonal_complex` field** (checked:
`fields` lists only `ST`), and the isolate database has no such field either.
The complex is therefore derived.

**Naive eBURST-style clustering was tried and rejected.** Taking connected
components of the single-locus-variant graph at the conventional 6/7 threshold
over the November 2026 profile table (9503 complete profiles) gives:

```
components: 1046   largest sizes: 4456, 775, 568, 487, 433, 327, 203, 137 ...
```

The largest component holds **4456 of 9503 STs** and contains ST1, ST5, ST8,
ST15, ST80, ST97, ST239, ST250 and ST254 together, with ST5 as its predicted
founder. It would report ST8 as "CC5". Connected components are unusable at
this database size.

What is used instead is the **named-founder single-locus-variant rule**: a query
ST is assigned to a complex when it *is* a founder ST of an established complex,
or shares 6 of 7 loci with one. On the same table this assigns 6132/9503 STs
(64.5%), leaves the rest as `no established complex` rather than forcing them,
and resolves every validation strain correctly — including ST36 → CC30,
ST250 → CC8 and ST254 → CC8, which are single-locus variants of their founders,
not the founders themselves. An ST that is an SLV of two founders is reported as
`CC8 or CC239`-style, not silently resolved; 208 STs are in that position.
The 40 founder STs are in `CC_FOUNDERS`.

`clonal_complex(sequence_type, db_dir)` is **public**, so the integrator can
pass Hydra's existing `MlstCall.sequence_type` straight in and skip the module's
own MLST pass (§7).

---

## 3. Limits when working from assemblies

- **spa is a tandem array, and short-read assemblers collapse it.** Mollerup
  2022 measured 84.8% / 97.6% (SPAdes, two settings) and 98.6% (SKESA) correct
  on 1754 MRSA, and found misassembled types were "generally two *spa* repeats
  shorter", concentrated in types with repetitions of the same repeat. Bartels
  2014 found 97% agreement with Sanger on 699 isolates, discordance mostly from
  lost 24 bp repeats. Bletz 2015 showed the pipeline dominates: 99.1% versus
  58.9% on the same 423 isolates with different software. **Every spa call
  therefore carries the raw succession** so a collapse is auditable, and a call
  is refused outright at a contig edge. A single-repeat difference between two
  isolates should not be acted on without confirmation.
- The whole module is blind to **expression and regulation**. An intact *cap5*
  locus is not a capsule; an intact *hla* is not α-toxin production; *agr*
  group is not *agr* function. N315 is a published *agr*-defective strain and
  this module correctly types it as group II while saying nothing about its
  function — RefSeq flags its *agrC* as a pseudogene. Only `agrA` coverage is
  checked, so *agr* dysfunction is **not** systematically detected.
- **Plasmid and phage genes depend on the assembly containing them.** The *bla*
  operon is usually plasmid-borne; a plasmid dropped by the assembler reads as
  `absent`.
- For *spa*, the "interrupted" lesion mechanism is reported loosely: a gene with
  a tandem repeat region fragments its own alignment, so the piece count is
  evidence that the gene is not reference-like rather than proof of an
  insertion. The caveat wording reflects that.
- `not determined` is used wherever the assembly cannot resolve the question. It
  is never replaced by a negative.

---

## 4. Database provenance and licence

| Database | Source | Version | Licence | Redistribution |
|---|---|---|---|---|
| `spa` | `https://spa.ridom.de/dynamic/sparepeats.fasta`, `.../spatypes.txt` | retrieved 2026-11-17; upstream Last-Modified 2026-11-12; 865 repeat codes, 22 886 types | Ridom spaServer Server Policy | **restricted** |
| `markers` | NCBI RefSeq, efetch coordinate slices of NC_007793.1, NC_002745.2, NC_003140.1, NC_003923.1, NC_002951.2, NC_009641.1 | accession-and-coordinate pinned in `MARKERS` | NCBI public domain | permitted |
| `mlst` | PubMLST `pubmlst_saureus_seqdef` scheme 1 profiles + 7 allele sets | retrieved 2026-11-17; scheme `last_updated` 2026-11-13; 9503 of 12 583 profiles (anonymous) | PubMLST terms, pre-2025 data | permitted — but see below |

Total built size: **6.2 MB** (58 markers ≈ 70 kB, 7492 MLST alleles ≈ 3.5 MB,
Ridom files ≈ 0.86 MB, plus BLAST indices).

### The Ridom licence decides how spa ships

From `https://spa.ridom.de/policy.shtml`, verbatim:

> "The nomenclature data (spa-repeat and -type assignments) of the SpaServer is
> the property of Ridom GmbH."

> "Individuals may download the database and use the content for non-commercial
> purposes only. Copies of the SpaServer database may not be used without
> permission to produce database driven web server or any other product or
> service."

> "However, this website must be acknowledged in all publications."

> "The preferred format of the acknowledgement is: 'This publication made use of
> the spa typing website (https://spaserver.ridom.de/) that is developed by
> Ridom GmbH and curated by SeqNet.org.'"

**Consequence: the two Ridom files must not be vendored into the HYDRA-Win
zip.** "Copies … may not be used without permission to produce … any other
product or service" is the binding sentence. `fetch()` downloads them into the
user's own database directory, which is the permitted "individuals may download"
case. If Hydra is ever distributed commercially, written permission from Ridom
GmbH is needed. Reports that use a spa call should carry the acknowledgement
above and cite Harmsen 2003.

Note **HTTPS only**: the port-80 virtual host at `spa.ridom.de` returns 404 for
`/dynamic/`.

### The PubMLST terms changed in May 2026, and the date split is the point

`https://pubmlst.org/about/terms` now 404s; the terms are at
`https://pubmlst.org/terms-conditions` (updated 2026-05-13) and **name no open
licence**. They split the data by submission date:

> "All data submitted on or before 31 December 2024 may continue to be
> downloaded, used, and redistributed without restriction, subject only to
> proper scientific citation and acknowledgement."

> "Data submitted on or after 1 January 2025 are governed by the licensing and
> usage restrictions described in Sections 4 and 5."

For that later data:

> "Downloading nomenclature data, including allele sequences, allele
> definitions, profiles, or sequence types, is permitted only for
> non-commercial academic research or public health surveillance use."

> "Restricted access also applies to commercial services and software tools
> that enable the use or the downloading of data … The publishers of such
> commercial tools will require a specific licence contract."

**`fetch()` downloads anonymously, and that is deliberate.** Anonymous access
returns only the pre-2025 data, which the terms place under no restriction, so
this sub-database is redistributable. Adding an API key to reach the remaining
~3080 profiles would pull it under the post-2024 restrictions, and a commercial
Hydra would then need a licence contract from PubMLST. The visible cost of
staying anonymous is that an ST assigned after 2024-12-31 reads as
`not determined`, with that reason given in the note. The earlier
`redistribution: permitted` is correct **only** on that basis.

---

## 5. Validation

Ten complete public genomes, run through `run()` as shipped. Reproduce with:

```
HYDRA_VALIDATION=1 pytest tests/test_typing_saureus.py
```

### 5.1 Obtained

| strain | spa | agr | CC | capsule | pvl | tst | hlb | iec | acme | blaZ |
|---|---|---|---|---|---|---|---|---|---|---|
| USA300_FPR3757 `GCF_000013465.1` | t008 | I | CC8 | 5 | positive | negative | disrupted | type B | present | absent |
| N315 `GCF_000009645.1` | t002 | II | CC5 | 5 | negative | positive | disrupted | type F | absent | present |
| MW2 `GCF_000011265.1` | **t128** | III | CC1 | 8 | positive | negative | disrupted | type D | absent | absent |
| Newman `GCF_000010465.1` | t008 | I | CC8 | 5 | negative | negative | disrupted | type A | absent | absent |
| MRSA252 `GCF_000011505.1` | t018 | III | CC30 | 8 | negative | negative | disrupted | type A | absent | present |
| MSSA476 `GCF_000011525.1` | t607 | III | CC1 | 8 | negative | negative | disrupted | type D | absent | present |
| COL `GCF_000012045.1` | t008 | I | CC8 | 5 | negative | negative | **intact** | **absent** | absent | absent |
| TW20 `GCF_000027045.1` | t037 | I | CC239 | 8 | negative | negative | disrupted | type D | absent | present |
| NCTC 8325 `GCF_000013425.1` | t211 | I | CC8 | 5 | negative | negative | disrupted | type B | absent | absent |
| RF122 `GCF_000009005.1` | t529 ⚠ | II | CC151 | 8 | negative | positive | **intact** | **absent** | absent | absent |

### 5.2 Expected, with the source of each expectation

| Call | Strain | Expected | Source of the expectation | Obtained | Agrees |
|---|---|---|---|---|---|
| blaZ allotype | N315, MRSA252, MSSA476, TW20 | one of the four canonical residue pairs, with the class A motifs intact | Aarris 2025 Table 1 for the table; an independent survey of 60 distinct 281-aa NCBI BlaZ proteins found **only** the four canonical pairs (T/S n=28, T/N n=13, A/S n=7, K/N n=7) and the motifs in 60/60 | C, A, C, A — all four passed the motif checks | yes |
| spa | USA300 | t008 | Beckley & Amarh 2026 Table 1 (lineage-level) | t008 | yes |
| spa | N315 | t002 | Beckley & Amarh 2026 Table 1, citing GCA_000009645.1 | t002 | yes |
| spa | MRSA252 | t018 | Rouha 2018 Table 1 and Beckley & Amarh 2026 Table 1 | t018 | yes |
| spa | Newman | t008 | Rouha 2018 Table 1 ("Newman ST254-t008 MSSA") | t008 | yes |
| spa | TW20 | t037 | Beckley & Amarh 2026 Table 1, citing GCA_000027045.1 | t037 | yes |
| spa | NCTC 8325 | t211 | Beckley & Amarh 2026 Table 1, citing GCA_000013425.1 | t211 | yes |
| spa | MW2 | t125 | Rouha 2018 Table 1 | **t128** | **no — see §5.3** |
| spa | RF122 | pseudogene, expected to fail or mis-call | Herron-Olson 2007 (*spa* carries a premature stop codon in ET3-1) | t529, flagged twice | behaves as intended |
| agr | Newman | I | Tan 2022 | I | yes |
| agr | N315 | II | Tan 2022 | II | yes |
| agr | MW2 | III | Tan 2022 | III | yes |
| agr | RF122 | *agrB*/*agrD* group II, *agrC* 5′ unlike all four groups | Herron-Olson 2007 | II (98.0%) | consistent with the *agrB*/*agrD* evidence |
| agr | COL, N315, MRSA252, MSSA476 | I, II, III, III | Gilot 2002, by *agr* locus PCR | I, II, III, III | yes |
| agr | USA300, Newman, RF122 | I, I, II | independent re-typing of these accessions with AgrVATE's own 31-mer motif database (Raghuram 2022), zero cross-group hits | I, I, II | yes |
| pvl | USA300 | positive | Diep 2006 | positive | yes |
| pvl | MW2 | positive | Baba 2002 | positive | yes |
| pvl | Newman | negative | Baba 2008; Rouha 2018 Table 1 | negative | yes |
| pvl | MRSA252 | negative | Rouha 2018 Table 1 | negative | yes |
| pvl | MSSA476 | negative | Holden 2004 | negative | yes |
| pvl | RF122 | negative (carries *lukM*/*lukF*′, not PVL) | Herron-Olson 2007 | negative | yes — the bovine leukocidins did not cross-call |
| tst | N315 | positive | Kuroda 2001 | positive | yes |
| tst | RF122 | positive (SaPIbov1) | Herron-Olson 2007 | positive | yes |
| blaZ | COL | absent ("penicillinase-negative") | Gill 2005 | absent | yes |
| blaZ | MSSA476 | present (pSAS1) | Holden 2004 | present | yes |
| blaZ | N315 | present (pN315) | Kuroda 2001 | present | yes |
| blaZ | MRSA252 | present | Holden 2004 | present | yes |
| blaZ | Newman | absent | Baba 2008 ("the absence of drug resistance genes reflects the general antibiotic-susceptible phenotype") | absent | yes |
| ST | all | ST8 / 5 / 1 / 254 / 36 / 1 / 250 / 239 / 8 / 151 | PubMLST isolate records 37463, 35874, 35870, 35868, 104, 39207, 39197, 39192, 37464, 39194 | same, all 7 loci exact | yes |
| hlb | COL | intact (no φSa3; see `iec` absent) | RefSeq annotation of NC_002951.2 is the only panel member **not** flagged pseudo | intact | yes |
| hlb | 8 human strains | disrupted | RefSeq flags *sphingomyelin phosphodiesterase* as pseudo with a truncated 5′ end in all eight | disrupted | yes |

**Honest labelling of the weaker rows.** Capsule, clonal complex and IEC type
are *not* published per-strain for these strains — searched and not found. Those
expectations are derived:

- **capsule** from the independent RefSeq naming rule (strains whose *capHIJ*
  slot is named `cap8H/I/J` are serotype 8; strains where it holds unnamed
  `CatB-related O-acetyltransferase` / `glycosyltransferase` /
  `O-antigen ligase` products carry *cap5HIJ*), which agrees with the module on
  all ten, and from the lineage conventions (CC5 and CC8 → serotype 5; CC1 and
  CC30 → serotype 8).
- **clonal complex** from the documented founder rule plus the standard
  lineage names (ST8→CC8, ST5→CC5, ST1→CC1, ST36→CC30, ST250/ST254→CC8,
  ST151→CC151). PubMLST has no field to check against.
- **IEC type** from the gene patterns, cross-checked against the `hlb` call:
  COL and RF122 are the two genomes with an intact *hlb*, and they are exactly
  the two with no cluster. That internal consistency across two independent
  signals is the evidence.

These are regression assertions in the test suite, and are labelled in this
table as derived rather than published. They should not be read as validated
against the literature.

### 5.3 The MW2 disagreement, investigated

Rouha 2018 Table 1 lists strain MW2 (ATCC BAA-1707) as **t125**. The module
reports **t128** from `GCF_000011265.1`. Investigated in the sequence:

```
MW2 region NC_003923.1:99753-100099, 347 bp, plus-strand orientation
  pos   73 +24 r07  GAGGAAGACAACAACAAACCTGGT
  pos   97 +24 r23  AAAGAAGACGGCAACAAACCTGGC
  pos  121 +24 r23  AAAGAAGACGGCAACAAACCTGGC
  pos  145 +24 r21  AAAGAAGACAACAACAAGCCTGGC
  pos  169 +24 r16  AAAGAAGACGGCAACAAACCTGGT
  pos  193 +24 r34  AAAGAAGACAACAAAAAACCTGGT
  pos  217 +24 r33  AAAGAAGATGGCAACAAGCCTGGC
  pos  241 +24 r13  AAAGAAGACAACAACAAACCTGGT
  73 bp lead + 8x24 bp contiguous repeats + 82 bp trail = 347 of 347 bp
```

Findings:

1. Every base of the region is accounted for. There are no gaps, no
   unrecognised bases inside the array, and both flanking anchors are present
   at a mid-chromosome position.
2. `NC_003923.1` is a **closed, finished genome**, so assembler collapse cannot
   be the explanation.
3. `t128` = `07-23-23-21-16-34-33-13`; `t125` = `07-23-23-23-23-21-13`. These
   are **not** repeat-count variants of one another — t125 lacks `16-34-33`
   entirely — so no gain or loss of repeats converts one into the other.
4. MSSA476, the other ST1 reference genome, gives `07-16-23-21-16-34-33-13`
   (t607) by the same method: the same `…-21-16-34-33-13` family as MW2,
   differing at one repeat. The two closed ST1 genomes corroborate each other.

Conclusion: the deposited MW2 genome encodes t128. The t125 attribution refers
to a different working stock, or is a lineage-level attribution rather than a
measurement on this sequence — Rouha 2018's companion rows are lineage
associations, and the same paper lists the other USA300 reference strain
TCH1516 as t622, so the ST1/USA400 and ST8/USA300 lineages are not
spa-monotypic. **Not tuned away**; the test asserts t128 and this section
records why.

### 5.4 Noted but unresolved

- **USA300_FPR3757 *blaZ* absent.** The assembly includes all three plasmids
  (pUSA01, pUSA02, pUSA03) and carries no *bla* operon at all, so the call is
  right for this record. No published statement either way was found, and
  USA300 is commonly penicillin-resistant. Flagged for the owner rather than
  claimed as a validated negative.
- **COL *tst*.** Gill 2005 Table 1 lists νSa1 content as "seb, tsst, ear" with
  COL coordinates, but the current RefSeq annotation of NC_002951.2 contains no
  toxic-shock toxin gene, and the module reports *tst* negative. Unresolved
  conflict between the 2005 paper and the current annotation.
- **NCTC 8325** has only 219 gene symbols for 2767 CDS, so annotation-derived
  negatives from it are weak evidence.

---

## 6. Not shipped, and why

Left out deliberately. The brief's rule is not to ship an unvalidated call.

| Not shipped | Reason |
|---|---|
| **BlaZ allotypes beyond the residue pair** | Olsen 2006's 11 BlaZ protein signatures are not implemented; only the kinetic allotype from Ambler 128/216 is. No strain with a *published* allotype was available as a validation genome (ATCC 29213 is the published type A reference), so the four panel calls are validated against the rule and the 60-sequence residue survey, not against a published per-strain allotype. |
| **Direct cefazolin-inoculum-effect prediction** | Allotype A/C is reported with that association stated, but the effect is not predicted as a phenotype: it needs more than the allotype. |
| **Exfoliative toxins *eta*, *etb*, *etd*, *ete*** | None of the ten validation genomes carries them, so the calls could not be validated against anything but the annotation the references came from — circular. Resolved accessions ready to add: *eta* `WP_457786830`/`WP_001065781`, *etd* `WP_000164635` (IPG reports give nucleotide accessions and coordinates). Needs an SSSS isolate with a published toxin status. |
| **Enterotoxins *sed*, *see*, *sej*, *ser*** | Plasmid-borne; absent from the panel genomes; same validation problem. The `enterotoxins` note states that a negative covers only the genes listed. |
| ***edinA/B/C*, *lukM*** | Same reason. RF122 carries *lukM* and is available as a validation genome if it is added. |
| **ACME type I/II/III** | Needs the *opp3* cluster; the ACME *arc* cluster presence call is shipped and validated, the type is not. |
| **agr functional status** | Only `agrA` coverage is checked. Systematic detection of *agr* dysfunction (e.g. N315's pseudogenised *agrC*) needs per-gene lesion calls across *agrABCD*, which is a straightforward extension of `_describe_lesion`. |
| **BURP clustering / spa-CC** | Mellmann 2007; needs the whole spaServer type set and a clustering pass. Out of scope for a per-isolate call. |
| **SCCmec** | Already typed by `hydra_amr.typing.cassette`. |

---

## 7. For the integrator

Nothing outside the four added files was changed. Needed in shared files:

1. **Register the module.** `src/hydra_amr/typing/organisms/` is picked up by
   `[tool.setuptools.packages.find] where = ["src"]` with no pyproject change,
   but a dispatcher is needed that selects modules by
   `SpeciesCall.name ∈ SPECIES` or `MlstCall.scheme ∈ MLST_SCHEMES` and calls
   `fetch()` / `run()`.
2. **Pass the existing MLST call in.** `run()`'s signature is fixed by the brief
   and receives no `MlstCall`, so the module types the seven MLST loci itself
   purely to derive the clonal complex — duplicating work Hydra has already
   done. `clonal_complex(sequence_type, db_dir)` is public for exactly this:
   call it with `result.mlst.sequence_type` and the `mlst` sub-database (3.5 MB
   of alleles) becomes unnecessary. A future `run(..., mlst=None)` keyword would
   be the clean fix.
3. **Do not vendor the Ridom spa files** into the Windows zip — §4. The
   `markers` and `mlst` databases are freely redistributable; `spa` is not.
   `DATABASES[0].redistribution == "restricted"` is machine-readable for this.
4. **Surface the Ridom acknowledgement** in any report containing a spa call.
5. **Interaction with the SCCmec typer.** No overlap: this module emits no
   SCCmec call and `cassette.py` emits none of these schemes, so both can run.
   Two couplings are worth wiring: (a) *mecA* status from the cassette typer
   plus `bla_operon` from here together give the whole β-lactam picture, and the
   `blaZ` note already defers to "no mecA/mecC" for its penicillin statement;
   (b) `acme` is only interpretable next to the cassette call, since ACME sits
   adjacent to SCCmec in USA300 — reporting them side by side is the useful
   layout. Also note `mecA` is an **ambiguous RefSeq gene symbol**: it denotes
   both PBP2a and an unrelated ClpC adaptor present in every strain including
   MSSA, so any gene-symbol matching must use the product string.
6. **Running the tests inside a nested clone** needs
   `-o pythonpath=<clone>/src --rootdir=<clone>`, because an outer
   `pyproject.toml` with `[tool.pytest.ini_options] pythonpath = ["src"]` is
   found first and shadows the clone. Not an issue in the real repository root.

---

## 8. Bibliography

Every reference below was verified against PubMed and Crossref on 2026-11-17:
the PMID and DOI resolve, and first author, title, journal and year match.
Three PMIDs that would have been cited from memory were **wrong** and were
corrected — recorded here so they are not reintroduced: Herron-Olson 2007 is
**17971880** (not 18030345), Schijffelen 2010 is **20546576** (not 20170493),
Holden 2010 is **19948800** (not 20139188).

**spa typing**

1. **Harmsen D, Claus H, Witte W, Rothgänger J, Claus H, Turnwald D, Vogel U.**
   Typing of methicillin-resistant *Staphylococcus aureus* in a university
   hospital setting by using novel software for *spa* repeat determination and
   database management. *J Clin Microbiol* 2003;41:5442–5448.
   PMID 14662923 · doi:10.1128/JCM.41.12.5442-5448.2003
   — *the repeat-code nomenclature and the t-number assignments; the citation the
   Ridom policy requires.* Note the journal wording is "**by using** novel
   software"; the Ridom policy page renders it "using a novel software". Note
   also that BURP was introduced later (Mellmann 2007, *BMC Microbiol* 7:98),
   **not** in this paper.
2. **Shopsin B, Gomez M, Montgomery SO, et al.** Evaluation of protein A gene
   polymorphic region DNA sequencing for typing of *Staphylococcus aureus*
   strains. *J Clin Microbiol* 1999;37:3556–3563.
   PMID 10523551 · doi:10.1128/JCM.37.11.3556-3563.1999
   — *spa repeat-region sequencing as a method, and the repeat structure the
   tiling relies on.*
3. **Frénay HME, Bunschoten AE, Schouls LM, et al.** Molecular typing of
   methicillin-resistant *Staphylococcus aureus* on the basis of protein A gene
   polymorphism. *Eur J Clin Microbiol Infect Dis* 1996;15:60–64.
   PMID 8641305 · doi:10.1007/BF01586186
   — *the original protein A polymorphic-region typing, before the Ridom codes.*

**spa from draft assemblies — the caveat on every call**

4. **Mollerup S, Worning P, Petersen A, Bartels MD.** *spa* typing of
   methicillin-resistant *Staphylococcus aureus* based on whole-genome
   sequencing: the impact of the assembler. *Microbiol Spectr*
   2022;10:e0218922. PMID 36350148 · doi:10.1128/spectrum.02189-22
   — *that a collapsed array is typically two repeats short and that the error
   concentrates in successions with repeated identical repeats; the basis for
   reporting the raw succession and refusing a call at a contig edge.*
5. **Bartels MD, Petersen A, Worning P, et al.** Comparing whole-genome
   sequencing with Sanger sequencing for *spa* typing of methicillin-resistant
   *Staphylococcus aureus*. *J Clin Microbiol* 2014;52:4305–4308.
   PMID 25297335 · doi:10.1128/JCM.01979-14
   — *the measured rate of lost 24 bp repeats in short-read assemblies.*
6. **Bletz S, Mellmann A, Rothgänger J, Harmsen D.** Ensuring backwards
   compatibility: traditional genotyping efforts in the era of whole genome
   sequencing. *Clin Microbiol Infect* 2015;21:347.e1–347.e4.
   PMID 25658529 · doi:10.1016/j.cmi.2014.11.005
   — *that spa typing accuracy from WGS is dominated by the assembly pipeline.*

**MLST and clonal complexes**

7. **Enright MC, Day NP, Davies CE, Peacock SJ, Spratt BG.** Multilocus sequence
   typing for characterization of methicillin-resistant and
   methicillin-susceptible clones of *Staphylococcus aureus*. *J Clin Microbiol*
   2000;38:1008–1015. PMID 10698988 · doi:10.1128/JCM.38.3.1008-1015.2000
   — *the 7-locus scheme the clonal-complex call reads.*
8. **Feil EJ, Li BC, Aanensen DM, Hanage WP, Spratt BG.** eBURST: inferring
   patterns of evolutionary descent among clusters of related bacterial
   genotypes from multilocus sequence typing data. *J Bacteriol*
   2004;186:1518–1530. PMID 14973027 · doi:10.1128/JB.186.5.1518-1530.2004
   — *the single-locus-variant definition of a complex and the founder-naming
   convention; also the method whose naive connected-component form is rejected
   in §2.9.*

**blaZ allotypes and the bla operon**

9. **Ambler RP, Coulson AF, Frère JM, et al.** A standard numbering scheme for
   the class A β-lactamases. *Biochem J* 1991;276(Pt 1):269–270.
   PMID 2039479 · doi:10.1042/bj2760269
   — *the numbering the allotype positions 128 and 216 are expressed in.*
10. **Voladri RK, Tummuru MK, Kernodle DS.** Structure-function relationships
    among wild-type variants of *Staphylococcus aureus* β-lactamase: importance
    of amino acids 128 and 216. *J Bacteriol* 1996;178:7248–7253.
    PMID 8955409 · doi:10.1128/jb.178.24.7248-7253.1996
    — *the origin of the allotype scheme and of positions 128 and 216.* Only the
    abstract was obtainable (PMC holds it as a scanned PDF behind a CAPTCHA), so
    the residue table implemented here is taken from ref 12, not from this paper.
11. **Carvajal LP, Rincon S, Echeverri AM, et al.** Novel insights into the
    classification of staphylococcal β-lactamases in relation to the cefazolin
    inoculum effect. *Antimicrob Agents Chemother* 2020;64:e02511-19.
    PMID 32071048 · doi:10.1128/AAC.02511-19
    — *the mapping of Ambler 128/216 onto sequential residues 119/207 of the
    281-residue preproenzyme, and the cefazolin-inoculum-effect association.*
    Its Discussion contains a parenthetical that contradicts its own Methods;
    the Methods mapping is the one implemented, and it is confirmed by the
    motif positions.
12. **Aarris M, et al.** Genetic variation in the *blaZ* gene leading to the
    BORSA phenotype in *Staphylococcus aureus*. *Antibiotics (Basel)*
    2025;14:449. PMID 40426516 · doi:10.3390/antibiotics14050449
    — *Table 1, "blaZ variants in S. aureus according to Ambler's standard
    numbering scheme": the complete allotype A–F residue table implemented here.*
13. **Olsen JE, Christensen H, Aarestrup FM.** Diversity and evolution of *blaZ*
    from *Staphylococcus aureus* and coagulase-negative staphylococci.
    *J Antimicrob Chemother* 2006;57:450–460. PMID 16449305 ·
    doi:10.1093/jac/dki492
    — *cited for the diversity of blaZ (69 sequence types, 11 BlaZ protein
    types) and for the separate plasmid and chromosomal evolutionary lines.*
    **Not** the source of allotypes A–D, contrary to common citation.

**Immune evasion cluster and the hlb-converting prophage**

14. **van Wamel WJ, Rooijakkers SH, Ruyken M, van Kessel KP, van Strijp JA.** The
   innate immune modulators staphylococcal complement inhibitor and chemotaxis
   inhibitory protein of *Staphylococcus aureus* are located on
   β-hemolysin-converting bacteriophages. *J Bacteriol* 2006;188:1310–1315.
   PMID 16452413 · doi:10.1128/JB.188.4.1310-1315.2006
   — *the IEC type A–G definitions, and that the prophage carrying them
   inactivates hlb.*

**agr**

15. **Jarraud S, Lyon GJ, Figueiredo AM, et al.** Exfoliatin-producing strains
    define a fourth *agr* specificity group in *Staphylococcus aureus*.
    *J Bacteriol* 2000;182:6517–6522. PMID 11053400 ·
    doi:10.1128/JB.182.22.6517-6522.2000
    — *the definition of agr group IV, and GenBank AF288215 (strain RN5881),
    from which the group IV reference window is taken.* This, not Dufour 2002,
    is the group IV paper.
16. **Jarraud S, Mougel C, Thioulouse J, et al.** Relationships between
    *Staphylococcus aureus* genetic background, virulence factors, *agr* groups
    (alleles), and human disease. *Infect Immun* 2002;70:631–641.
    PMID 11796592 · doi:10.1128/IAI.70.2.631-641.2002
    — *the conventional typing region: the Pan-1/Pan-2 primers amplify "a
    1,234-bp agr fragment encompassing the 3′ end of agrB, all of agrD, and the
    5′ end of agrC", which is the window used here.*
17. **Dufour P, Jarraud S, Vandenesch F, et al.** High genetic variability of the
    *agr* locus in *Staphylococcus* species. *J Bacteriol* 2002;184:1180–1186.
    PMID 11807079 · doi:10.1128/jb.184.4.1180-1186.2002
    — *that AgrB, AgrC and AgrD are "strikingly variable among the four groups"
    while AgrA is "much more conserved", which is why agrA is excluded from the
    group call.* Note the title is "Staphylococcus **Species**", not
    "S. aureus" — a common miscitation.
18. **Gilot P, Lina G, Cochard T, Poutrel B.** Analysis of the genetic
    variability of genes encoding the RNA III-activating components Agr and
    TRAP in a population of *Staphylococcus aureus* strains isolated from cows
    with mastitis. *J Clin Microbiol* 2002;40:4060–4067. PMID 12409375 ·
    doi:10.1128/JCM.40.11.4060-4067.2002
    — *the published agr group assignments for COL (I), N315 (II), MRSA252 (III)
    and MSSA476 (III) used as validation expectations.*
19. **Raghuram V, Alexander AM, Loo HQ, Petit RA 3rd, Goldberg JB, Read TD.**
    Species-wide phylogenomics of the *Staphylococcus aureus* Agr operon
    revealed convergent evolution of frameshift mutations. *Microbiol Spectr*
    2022;10:e0133421. PMID 35315691 · doi:10.1128/spectrum.01334-21
    — *AgrVATE, the reference tool whose 31-mer motif database was used to
    re-type USA300, Newman and RF122 independently as a cross-check.* AgrVATE
    is MIT-licensed but needs non-Python external programs, so it is not a
    dependency here.
20. **Tan L, Li SR, Jiang B, Hu XM, Li S.** Therapeutic targeting of the
    *Staphylococcus aureus* accessory gene regulator (*agr*) system.
    *Front Microbiol* 2022;13:700894. PMID 35295312 ·
    doi:10.3389/fmicb.2022.700894
    — *the published agr group assignments for the prototype strains Newman (I),
    N315 (II) and MW2 (III) used as the group references and expectations.*

**Reference genomes the markers and the validation panel come from**

21. **Diep BA, Gill SR, Chang RF, et al.** Complete genome sequence of USA300, an
    epidemic clone of community-acquired meticillin-resistant *Staphylococcus
    aureus*. *Lancet* 2006;367:731–739. PMID 16517273 ·
    doi:10.1016/S0140-6736(06)68231-7
    — *USA300_FPR3757, source of most markers; PVL-positive; and the ACME the
    acme call detects.*
22. **Kuroda M, Ohta T, Uchiyama I, et al.** Whole genome sequencing of
    meticillin-resistant *Staphylococcus aureus*. *Lancet* 2001;357:1225–1240.
    PMID 11418146 · doi:10.1016/S0140-6736(00)04403-2
    — *N315 and plasmid pN315: the tst, egc, capsule 5 and bla operon
    references.*
23. **Baba T, Takeuchi F, Kuroda M, et al.** Genome and virulence determinants of
    high virulence community-acquired MRSA. *Lancet* 2002;359:1819–1827.
    PMID 12044378 · doi:10.1016/S0140-6736(02)08713-5
    — *MW2: the PVL-positive, agr group III and capsule 8 references.*
24. **Baba T, Bae T, Schneewind O, Takeuchi F, Hiramatsu K.** Genome sequence of
    *Staphylococcus aureus* strain Newman and comparative analysis of
    staphylococcal genomes. *J Bacteriol* 2008;190:300–310. PMID 17951380 ·
    doi:10.1128/JB.01000-07
    — *Newman: MSSA, PVL-negative, drug-susceptible; the fully susceptible MSSA
    case on the panel.*
25. **Holden MTG, Feil EJ, Lindsay JA, et al.** Complete genomes of two clinical
    *Staphylococcus aureus* strains: evidence for the rapid evolution of
    virulence and drug resistance. *Proc Natl Acad Sci U S A*
    2004;101:9786–9791. PMID 15213324 · doi:10.1073/pnas.0402521101
    — *MRSA252 (ST36, EMRSA-16) and MSSA476 (ST1): ST, blaZ and PVL
    expectations.*
26. **Gill SR, Fouts DE, Archer GL, et al.** Insights on evolution of virulence
    and resistance from the complete genome analysis of an early
    methicillin-resistant *Staphylococcus aureus* strain and a
    biofilm-producing methicillin-resistant *Staphylococcus epidermidis*
    strain. *J Bacteriol* 2005;187:2426–2438. PMID 15774886 ·
    doi:10.1128/JB.187.7.2426-2438.2005
    — *COL: the only panel member with an intact hlb, from which the hlb
    reference allele is taken; and the published penicillinase-negative status.*
27. **Holden MTG, Lindsay JA, Corton C, et al.** Genome sequence of a recently
    emerged, highly transmissible, multi-antibiotic- and antiseptic-resistant
    variant of methicillin-resistant *Staphylococcus aureus*, sequence type 239
    (TW). *J Bacteriol* 2010;192:888–892. PMID 19948800 ·
    doi:10.1128/JB.01255-09
    — *TW20 (ST239), the CC8/CC30 hybrid that tests the clonal-complex rule's
    ambiguity handling.*
28. **Herron-Olson L, Fitzgerald JR, Musser JM, Kapur V.** Molecular correlates of
    host specialization in *Staphylococcus aureus*. *PLoS One* 2007;2:e1120.
    PMID 17971880 · doi:10.1371/journal.pone.0001120
    — *RF122 (ET3-1, ST151): the awkward case. spa is a pseudogene carrying a
    premature stop codon; agrB/agrD match group II while the agrC 5′ end differs
    from all four groups; tst-positive via SaPIbov1; carries lukM/lukF′ rather
    than PVL.*

**Published type tables used as expectations**

29. **Rouha H, Weber S, Janesch P, et al.** Disarming *Staphylococcus aureus* from
    destroying human cells by simultaneously neutralizing six cytotoxins with
    two human monoclonal antibodies. *Virulence* 2018;9:231–247.
    PMID 29099326 · doi:10.1080/21505594.2017.1391447
    — *Table 1: the MLST–SCCmec–spa types and lukSF-PV status for MRSA252
    (ST36-II-t018), MW2 (ST1-IV-t125) and Newman (ST254-t008). The MW2 row is
    the disagreement investigated in §5.3.* Crossref dates it 2017, PubMed 2018.
30. **Beckley B, Amarh V.** StaphScope: a genomic typing tool for *Staphylococcus
    aureus*. *BMC Genomics* 2026;27. PMID 41645058 ·
    doi:10.1186/s12864-026-12609-x
    — *Table 1, "Validation of staphscope using well-characterized reference
    genomes": the spa and ST expectations for USA300 (t008), N315 (t002),
    MRSA252 (t018), TW20 (t037) and NCTC 8325 (t211).* Its own "expected type"
    citations are lineage-level papers, so these are inherited associations
    rather than primary strain measurements; two of its rows use non-canonical
    re-sequenced assemblies.

**Database sources**

21. Ridom spaServer, `https://spaserver.ridom.de/` — repeat codes and type
    assignments; 865 repeats and 22 886 types at retrieval. Server Policy at
    `https://spa.ridom.de/policy.shtml`. Acknowledgement required (§4).
22. PubMLST *Staphylococcus aureus* database,
    `https://rest.pubmlst.org/db/pubmlst_saureus_seqdef` — scheme 1 profile
    table and allele sets; terms at `https://pubmlst.org/about/terms`.
23. NCBI RefSeq — the marker coordinate slices; accessions and ranges are listed
    per marker in `MARKERS` in the module.

**Reference tool considered for validation**

31. **Sanchez-Herrero JF, Sullivan M.** spaTyper: staphylococcal protein A (*spa*)
    characterization pipeline. Zenodo, 2020. doi:10.5281/zenodo.4063625
    (concept DOI 10.5281/zenodo.4063624) — no peer-reviewed paper exists.
    *Not used as the validation oracle.* It is unmaintained since 2021, its
    pip package ships a January 2021 snapshot (813 repeats, 19 696 types versus
    865 and 22 886 live) and never contacts Ridom, its hard-coded `http://`
    URLs now return 404, and its default mode does no primer anchoring — on a
    deliberately truncated repeat region it returned a real-looking type
    (t6058) with no warning. Validation here is against the published types in
    refs 19–20 and the primary strain papers instead. Its licensing chain is
    also inconsistent (LGPL-3.0 LICENSE, "GPLv3" in `setup.py`) and its upstream
    `mjsull/spa_typing` has no licence at all, which is a further reason the
    algorithm here is an independent implementation from the Ridom protocol
    rather than a port.
