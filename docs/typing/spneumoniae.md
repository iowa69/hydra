# *Streptococcus pneumoniae* typing

Module: `src/hydra_amr/typing/organisms/spneumoniae.py` · key `spneumoniae` ·
triggered by the species call *Streptococcus pneumoniae* or the PubMLST scheme
`spneumoniae`.

Assembly in. BLAST+ (`blastn`, `blastp`, `makeblastdb`) and pure Python out; no
other program, no compiled dependency, nothing that needs reads.

---

## 1. What is called, and what each call means

| Scheme | Example call | What it is for |
|---|---|---|
| `cps_serotype` | `19A`, `15B/15C`, `non-encapsulated` | the capsular serotype with its serogroup and its membership of every licensed vaccine formulation |
| `pbp_type` | `13:16:47`, `NEW:348:620` | the pbp1a:pbp2b:pbp2x transpeptidase allele triplet in CDC numbering |
| `beta_lactam_mic` | `penicillin MIC 2 mg/L (predicted)` | the beta-lactam MICs that PBP type predicts, with CLSI categories |
| `pilus_islets` | `PI-1+PI-2`, `none` | the two adhesive pilus islets |
| `macrolide_mechanism` | `erm(B)`, `mef` | which macrolide/lincosamide mechanism is present and the phenotype it implies |
| `other_determinants` | `tet(M), cat` | tetracycline and phenicol determinants, interpreted |
| `target_substitutions` | `folA I100L` | gyrA, parC, folA and rplD substitutions in published numbering |

### 1.1 Capsular serotype

The serotype decides whether the isolate is a vaccine type, and that is the
single most consequential thing a pneumococcal reference laboratory reports.
The call therefore carries, for every serotype or group:

* the **serogroup**, in the `lineage` field;
* which of **PCV7, PCV10 (GSK/Synflorix), PCV13, PCV15, PCV20, PCV21 and
  PPSV23** contain it;
* where a formulation does not contain it but is expected to cross-protect
  (the 6A conjugate against 6C; PCV20's 15B against 15C; PPSV23's "type 20"
  against both 20A and 20B), that is reported **separately** from membership
  and never folded into it;
* for a pair or group, which formulations cover **only part** of it. A `9A/9V`
  call reports "PCV13 (covers 9V only)", because 9V is in six formulations and
  9A is in none, and neither "covered" nor "not covered" would be true.

**Non-encapsulated genotypes** are reported as `non-encapsulated`, in two
distinguishable states: no recognisable locus at all (`lineage: null capsule`),
and a locus present but gutted (`lineage: disrupted serogroup N locus`), which
names the missing coordinates. Both mean serologically non-typeable and covered
by no vaccine. Which null-capsule clade it is — NCC1 (*pspK*), NCC2 (*aliC* and
*aliD*), NCC3 (*aliD* only) — is **not** determined; the module does not carry
those references, and says so in the note. A laboratory capsule-deletion
derivative such as R6 reads the same way as a natural non-encapsulated isolate.

### 1.2 PBP typing and the predicted MIC

`pbp_type` is the triplet the CDC Streptococcus Laboratory's scheme defines, and
is reported as `pbp1a:pbp2b:pbp2x`. A region with no exact match is `NEW`, with
the substitutions it carries listed and those falling in an active-site motif or
at a named resistance position called out.

`beta_lactam_mic` is a **prediction from a genotype**. Every string it produces
says so, and both the penicillin meningitis and non-meningitis criteria are
reported, because which applies is a clinical question about the site of
infection that the module cannot answer. It is read from a published lookup
table and is declined, not extrapolated, whenever the triplet is not in it.

### 1.3 The rest

`macrolide_mechanism` exists because Hydra's AMR screen already finds *erm(B)*
and *mef* as genes but does not say what they mean. The mechanism decides
whether clindamycin can be used: *erm(B)* and *erm(TR)* methylate the ribosome
and take clindamycin with them (MLS-B); *mef* is an efflux pump that leaves
clindamycin and the 16-membered macrolides active (M phenotype); *lsa(C)* does
the opposite and leaves the macrolides active (LS-A). With *erm(TR)* the
methylase is inducible, so a clindamycin-susceptible result is clinically
unreliable and the note says so.

---

## 2. Method, rules and thresholds

### 2.1 Serotype, stage one: the whole capsule locus

The assembly is searched (`blastn -task megablast`, identity floor 85%) against
the 106 reference *cps* loci of the Global Pneumococcal Sequencing project's
`cps_reference_database`. For each reference, the aligned subject intervals are
merged **across all contigs** — a draft assembly routinely breaks a 20 kb
capsule locus into two or three pieces, and the locus is no less present for
that — and coverage, length-weighted identity and the gaps are recorded.
References are ranked by coverage × identity.

| Threshold | Value | Why |
|---|---|---|
| `CPS_MIN_IDENTITY` | 85% | identity floor for alignments that count |
| `CPS_MIN_COVERAGE` | 80% | an intact locus. Every encapsulated genome on the validation panel covered its reference by 91.8% or more (91.8, 93.2, 96.1, 97.5, 98.1, 99.9, 99.9, 100 ×6); R6, whose locus carries a 7.5 kb internal deletion, covered the serotype 2 reference by 63.6%. The floor sits in that gap. |
| `CPS_ABSENT_COVERAGE` | 35% | below this nothing recognisable as a capsule locus is present. The conserved 5′ *wzg-wzh-wzd-wze* block and the *dexB*/*aliA* flanks are shared by every locus, so even a capsule-free genome covers a fifth to a third of any reference. |
| `CPS_TIE_MARGIN` | 2 points | two whole-locus matches this close on coverage × identity cannot be separated by coverage alone |

Stage one alone is **not** a serotype call. Validation showed exactly why: the
serotype 6B strain 670-6B covers the **6E** reference 100% at 100% identity, a
Korean 11A isolate covers the **11E** reference, a serotype 24F isolate covers
the **24C** reference, and a 35B isolate covers the **35D** reference. In each
case the probe rules of stage two gave the serologically correct answer and the
locus match is reported as corroborating evidence, not as the call.

### 2.2 Serotype, stage two: the CDC discriminating probe rules

The assembly is searched against the CDC Streptococcus Laboratory's 72
serotype-discriminating probes (`blastn -task blastn -word_size 11`, e-value
1e-5). The CDC pipeline runs this screen with `srst2 --min_coverage 99.8
--max_divergence 5`, so the BLAST equivalents are:

* **present**: ≥ 99.8% of the probe covered at ≥ 95% identity;
* **exact** (an empty "diffs" column in the CDC output): 100% identity over
  100% of the probe.

The probes are grouped into clusters and **one gene is taken per cluster**, as
srst2 reports: serogroup 19's four *wzy* alleles sit in one cluster and only the
closest counts. Keeping them all would let one genome satisfy the 19A and the
19F rule at once.

The rules are then those of `SPN_Serotyper.pl`, transcribed rule for rule:

| Rule | Criterion | Result |
|---|---|---|
| single target, differences allowed | probe present | `WZY1`→1, `WZY2`→2, `WCHE3`→3, `WZY4`→4, `WZY5`→5, `WZY8`→8, `WCRG10A`→10A, `GTF10F`→10F, `WZY11A`→11A, `WZY13`→13, `WZY14`→14, `WZY16F`→16F, `GTF17F`→17F, `WZY20`→20, `WZY21`→21, `WZY23A`→23A, `WZY23B`→23B, `WZY23F`→23F, `WZY28A`→28A, `WZY31`→31, `WZY34`→34, `WZY35A`→35A, `WCRO35F`→35F, `WCYV38`→38, `WHAI47F`→47F, `TTS`→37 |
| single target, no differences | probe exact | `WZY18A`→18A, `WZY18F`→18F, `WZY19A`/`WZY19AVAR`→19A, `WZY19F`/`WZY19FVAR`→19F |
| groups the probes cannot split | probe present | `WZY11B:11C`→11B/11C, `WCII12A:12B:46`→12A/12B/46, `WZY24F:24A:24B`→24A/24B/24F, `WCYE25A:25F`→25A/25F, `WZY35C:42`→35C/42 |
| serogroup 6 | `WCIP6AC` exact × `WCIN6AB` → **6A**; `WCIP6AC` exact × `WCIN6CD` → **6C**; `WCIP6BD` exact × `WCIN6AB` → **6B**; `WCIP6BD` exact × `WCIN6CD` → **6D** | *wciP* residue 195 (Ser in 6A/6C, Asn in 6B/6D) crossed with the *wciN* allele |
| 35B / 35D | `WZY35B` present; `WCIG35B` exact → **35B**, otherwise → **35B/35D** | *wciG* O-acetyltransferase intact or not |
| 7B / 7C | `WZY7C` present; `WCHF7C` present → **7C**, absent → **7B** | |
| 7A / 7F | `WZY7F` present; `WCWD7F` exact → **7F**, otherwise → **7A/7F** | *wcwD* insT at 587 in 7A |
| 9L / 9N | `WZY9N` present; `WCJA9L` exact → **9L**, `WCJA9N` exact → **9N** | |
| 9A / 9V | `WZY9V` present; `WCJE9V` exact → **9V**, otherwise → **9A/9V** | *wcjE* delG at 721 in 9A |
| 15A / 15F | `WZY15A` present; `WCIZ15F` exact → **15F**, `WCIZ15F` absent → **15A** | |
| 15B / 15C | `WZY15B` present → **15B/15C**, always | see §3; this is the one place the module deliberately **departs** from the CDC rule |
| 18B / 18C | `WZY18C` exact; `WCIX18C` exact → **18C**, otherwise → **18B/18C** | *wciX* TGA stop at codon 168–171 in 18B |
| 22A / 22F | `WZY22F` present; `WCWA22A` → **22A**, `WCWA22F` → **22F** | |
| 33A / 33F | `WZY33F` present; `WCJE33A` exact → **33A**, otherwise → **33F** | *wcjE* insTA at 433 in 33A |

Every rule that fires is reported. Where more than one fires, the results are
reported **together** with the note "More than one CDC rule fired", because
disagreement between two rules is a finding, not something to resolve by
precedence. If no rule fires, the stage-one locus match becomes the call, the
note says so in those words, and the serogroups with no published split are
reported as groups (see §3).

### 2.3 PBP typing

Numbering throughout is on the full-length **TIGR4** proteins
(GenBank `AE005672.3`): PBP1a = `SP_0369` / UniProt `Q04707`, 719 aa;
PBP2b = `SP_1673` / `P0A3M5`, 680 aa; PBP2x = `SP_0336` / `P14677`, 750 aa.

| Locus | Typed region | Length | Active-site motifs |
|---|---|---|---|
| pbp1a | residues 371–646 | 276 aa | S370TMK (S370 is **outside** the region), S466SN, K557TG |
| pbp2b | residues 379–655 | 277 aa | S386-V-V-K389, S443SN, K615TG |
| pbp2x | residues 229–586 | 358 aa | S337TMK, S395SN, K547SG |

Two notes on those coordinates. First, the CDC reference fragments
(`Ref_PBP_3.faa`) span 371–647, 379–656 and 229–587, which is the 914-residue
window the published alignment table uses, but the **allele** sets are one
residue shorter at each locus; the triplet must be read from the shorter
window, and it is. Second, the published motif list gives pbp2b's SXXK motif as
"S386TMK" and its KTG as beginning at 614. Neither matches the protein: the
motif is **SVVK** at 386–389 and the KTG begins at **615**. The protein wins.

The procedure, which is `ExtractGene.pl` + `PBP-Gene_Typer.pl`:

1. `blastn` the assembly against the three nucleotide reference fragments
   (`MOD_bLactam_resistance.fasta`). Accept the best hit at ≥ 95% coverage of
   the fragment (`-L 0.95`) and ≥ 50% identity (the documented `-I` default).
   The 50% floor matters: a mosaic *pbp1a* recombined from *S. mitis* is only
   about **77%** identical to the susceptible reference, and at an 80% floor the
   two most penicillin-resistant genomes on the validation panel failed to type
   at pbp1a at all.
2. Extend the contig coordinates by whatever length of the reference the
   alignment did not reach, so a region whose first or last codons diverge is
   still extracted at full length and in frame. Reverse-complement if needed.
3. Translate in frame 1 and trim to the region length.
4. `blastp` against the CDC allele set for that locus (533 pbp1a, 744 pbp2b,
   1186 pbp2x alleles, numbered from 0).
5. **An allele is named only on 100% identity over the full length of the
   allele, gap-free** — the CDC criterion. Anything above 50% identity over 50%
   of the allele is `NEW`. Below that the region is `-` (not determined, and
   probably not this gene). A stop codon in the translation is also `-`.

For a `NEW` region the substitutions are reported **against allele 0**, the
penicillin-susceptible sequence the scheme is anchored on — not against the
nearest allele, because "T338A relative to the susceptible reference" means
something and "A12T relative to allele 417" does not. The nearest allele and
its identity are reported separately.

The positions singled out as resistance-relevant are the CDC pipeline's own
empirical list: the 27 transpeptidase positions at which every isolate with a
penicillin MIC ≥ 4 mg/L differed from PBP type 2-0-2. Those with a structural
or mutagenesis study behind them carry that citation:

* **pbp1a** T371, E397, N405, T540, N546, A550, T574, S575, Q576, F577, L583,
  A585 — with T371A/S and the TSQF574-577→NTGY block from Smith & Klugman 1998
  and Contreras-Martel 2006;
* **pbp2b** Q427, T446, E476, T489 — T446A from Pagliero 2004 (present in every
  cefotaxime-resistant strain examined, 60% loss of penicillin G affinity),
  T446/E476 structurally from Contreras-Martel 2009;
* **pbp2x** R254, M256, T338, M339, I371, G382, R384, T401, N444, S531, T550,
  Q552, L565, S576 — T338/T550/Q552 from Gordon 2000, T338 mechanistically from
  Mouz 1998, M339F from Chesnel 2003, T338P/I371T/R384G from Smith & Klugman
  2005.

A position that is polymorphic in the allele set but has no primary source is
deliberately **not** labelled a resistance mutation. `L546V` and `T451A` were
considered and dropped for that reason (T451A appears in 1 of 744 pbp2b
alleles).

### 2.4 The predicted MIC

The triplet keys `Ref_PBPtype_MIC.csv`, the CDC reference table of **307** PBP
types with their MICs for penicillin, amoxicillin, meropenem, cefotaxime,
ceftriaxone and cefuroxime. Each MIC is then interpreted against the CLSI
breakpoints, transcribed from the log2 thresholds in the CDC pipeline's
`MIC_format_with_SIR.R` — the only place the CDC publishes **both** the
meningitis and the non-meningitis criteria (the two papers print one each, and
they print different ones):

| Drug | S | I | R |
|---|---|---|---|
| penicillin, meningitis | ≤ 0.06 | — | ≥ 0.12 |
| penicillin, non-meningitis parenteral | ≤ 2 | 4 | ≥ 8 |
| amoxicillin | ≤ 2 | 4 | ≥ 8 |
| meropenem | ≤ 0.25 | 0.5 | ≥ 1 |
| cefotaxime, meningitis | ≤ 0.5 | 1 | ≥ 2 |
| cefotaxime, non-meningitis | ≤ 1 | 2 | ≥ 4 |
| ceftriaxone, meningitis | ≤ 0.5 | 1 | ≥ 2 |
| ceftriaxone, non-meningitis | ≤ 1 | 2 | ≥ 4 |
| cefuroxime | ≤ 0.5 | 1 | ≥ 2 |

**No MIC is produced** when any region is `NEW` or `-`, or when the complete
triplet is absent from the table. Li et al. 2017 validated a random-forest
predictor for exactly that case, but it is distributed as an R model object and
is not reimplemented here, and extrapolating from a neighbouring triplet is not
valid because the three regions interact. On the validation panel 10 of 20
genomes got a MIC; the other 10 got "not determined" with the reason. That is
the intended behaviour.

**Accuracy, as published.** Neither CDC paper publishes per-antibiotic
agreement as numbers — in both, they are bar charts with no source-data table —
so what the module quotes is what the papers state in words:

* for isolates whose PBP type was in the training set, the mode-MIC model gave
  essential agreement (within one doubling dilution) **> 98%**, category
  agreement **> 94%**, major discrepancy **< 3%** and very major discrepancy
  **< 2%**, for all six beta-lactams across 2,528 isolates (Li et al. 2016);
* on an independent set of 1,781 isolates, essential agreement ranged from
  **97.6%** (penicillin, mode-MIC) to **99.8%** (ceftriaxone) and category
  agreement from **95.7%** (penicillin) to **99.1%** (ceftriaxone)
  (Li et al. 2017);
* PBP type explained **94–99%** of the variance in log2 MIC, against 73–91% for
  the MLST sequence type (Li et al. 2016, Table 2).

Definitions are the FDA ones the papers use: essential agreement = within ±1
twofold dilution; category agreement = same S/I/R; major = reference S and
prediction R; very major = reference R and prediction S.

### 2.4a Why a triplet absent from the table gets no MIC: the evidence

About half the validation panel gets no MIC, and the obvious way to close that
— take the MIC of the nearest-neighbour PBP type, which is roughly what a
"mode-MIC model for unseen types" means — **was searched for in the literature
before being considered, and the literature says not to do it.** This section
records that search and its outcome, so the gap is documented rather than
quietly filled.

**What was searched.** PubMed and Crossref, for a published and validated
method predicting pneumococcal beta-lactam MICs from novel, previously
uncharacterised or mosaic PBP transpeptidase types: the CDC work on PBP
transpeptidase signatures (Li 2016), the CDC paper written specifically about
new PBP types (Li 2017), the CDC surveillance applications (Metcalf 2016), the
current production pipeline paper (Hung 2025), independent non-US validations
(Eriksen 2023), and machine-learning MIC prediction from pneumococcal genomes
up to 2026. The full text of Li 2017 was read, not only its abstract.

**What exists.** Two models are published, both defined in Li 2016 and both
evaluated for *new* PBP types in Li 2017:

* **Mode MIC (MM).** Per PBP type, the modal measured MIC (the highest of the
  modal values on a tie). A type absent from the training set is handled by a
  published, explicit nearest-neighbour rule — Li 2017's Methods, verbatim: *"A
  new PBP type that was not seen in the training dataset was approximated by a
  training PBP type showing highest amino acid identity."* So the
  nearest-neighbour idea is **not** an unpublished heuristic; it is part of the
  published model, and it has been validated.
* **Random forest (RF).** Six regression models, one per antibiotic, over the
  amino acid at each transpeptidase-domain position, with log2 MIC as a
  continuous response.

**What the validation of the implementable one says.** Li 2017 evaluated both
with a leave-one-**type**-out design built specifically to measure performance
on new PBP types, over 2,528 isolates (59 serotypes, 403 MLSTs). Quoting the
paper's own Results:

> "Results from the RF models showed higher than 97% EA for all 6 antibiotics,
> while results from the MM models showed EA lower than 90% for AMO and TAX
> … only the CA for TAX based on the MM model predictions was lower than 90%
> … Both RF and MM methods generated similar maj rates that were below 3%. The
> largest difference between the two models was the vmj rates. **Results from
> the MM model showed unacceptably high vmj rates for 4 of the 6 antibiotics,
> with 3 of them higher than 50%.**"

A very major discrepancy is a **resistant isolate predicted as susceptible** —
the error that gets a patient the wrong drug. The paper's FDA acceptance
criteria are essential and category agreement > 89.9 %, major ≤ 3 %, and for
the very major rate an upper 95 % confidence limit ≤ 7.5 % and a lower limit
≤ 1.5 %. On new PBP types the mode-MIC model fails that criterion for four of
the six beta-lactams, three of them by more than an order of magnitude. The
paper traces the failure to a single bad neighbour: isolates of PBP type
13-11-16 were approximated by the closely related type whose MICs are much
lower, and that one nearest-neighbour error destroyed the resistant calls.

Li 2017's own conclusion names the other model: *"MIC prediction using the RF
method could be an accurate alternative of phenotypic susceptibility testing
even in the presence of previously uncharacterized PBP types."*

**Why the validated method is not implemented here.** RF, the model the authors
endorse for new types, cannot be built under this project's rules (assemblies
in, pure Python plus BLAST+, must run inside a frozen Windows app). The CDC and
Pathogenwatch production implementations ship the six trained models as
serialized R objects totalling roughly 210 MB, require R with the
`randomForest` package **pinned to one release** because the authors state the
models were serialized with it and a newer release would change the
predictions, and the pipeline also uses an external multiple aligner. None of
that is portable into this tool, and a reimplementation of someone's random
forest without their trained model is a different model with no validation.

**So the position is:** the only MIC model implementable here is the one whose
published validation, on exactly the isolates it would be used for, is a
failure for four of six antibiotics in the most dangerous error direction. The
module therefore reports the triplet, the per-region substitutions and
`not determined` for the MIC, with the reason — and does not extrapolate. The
figures above are the evidence for that decision, not a limitation of the
implementation.

**What would change this.** Any of: the CDC publishing the random forest as
portable data rather than a serialized R object; a published MIC table
extending beyond the 307 PBP types of the 1998–2015 US invasive collection; or
a published, validated nearest-neighbour method with acceptable very-major
rates for the drugs that matter. A restricted extension — nearest-neighbour
mode MIC for **only** the antibiotics whose published very-major rate on new
types is acceptable — is arguable from these same figures, but it is a clinical
decision with a measured downside, it needs validation against published
per-strain MICs for isolates carrying novel types, and it is not taken here on
a developer's judgement.

**Independent context, for whoever revisits this.** The reference table covers
307 PBP types from US invasive isolates, 1998–2015, so novel types are commoner
elsewhere: Eriksen 2023 found no recognised profile in 37.5 % of 88 Danish
pneumococcal isolates, and in the largest independent evaluation of the
production random forest (Hung 2025, 1,273 Global Pneumococcal Sequencing
genomes with broth-dilution MICs) penicillin was the weakest of the six drugs.
Declining to predict will therefore fire often outside US invasive
surveillance, and that is the correct behaviour, not a defect.

### 2.5 Pilus islets, macrolide mechanism, other determinants

`pilus_islets` reads PI-1 from the CDC `RRGA` probe and PI-2 from `PITB`. This
is islet **presence from one marker gene each** — not a reconstruction of the
islet and not *rrgB* clade typing.

`macrolide_mechanism` and `other_determinants` screen the CDC resistance probe
set at the same thresholds as the serotype probes. The *mef* probe does **not**
separate *mef(A)* (Tn*1207.1*) from *mef(E)* (the Mega element); the two give the
same phenotype and differ epidemiologically, and the note says so and points at
Hydra's own AMR gene table for the specific allele and for *msr(D)*/*mel*.

### 2.6 Target-gene substitutions

For *gyrA*, *parC*, *folA* and *rplD* the probe region is located by `blastn`
(≥ 80% coverage, ≥ 75% identity, e-value 1e-3), extended to the full reference
length, translated in frame 1 and compared residue by residue with the CDC
reference peptide. Positions are reported on the full protein using the CDC
pipeline's own offsets:

| Gene | Reference peptide | Residues | Neutral polymorphisms dropped |
|---|---|---|---|
| gyrA | `VMGKYHPHGDSSIYEAMVRMAQWWSY` | 71–96 | none |
| parC | `DSSIYDAMVRMSQNWKNREILVEMHGNNGSMDGDPPAAMRYTEARLSEIAGYLLQDIEKKTVPFAWNFDD` | 78–147 | G128D, K137N, R95C |
| folA | `QDVQSVLDWYQDQEKNLYII` | 81–100 | D92A, D92V, D88G |
| rplD | `KPWRQKGTGRAR` | 63–74 | none |

So `gyrA S81F`, `parC S79F`, `parC D83N` and `folA I100L` come out with the
names the literature uses. The e-value has to be loose: a 54-base alignment at
83% identity — which is what a mosaic *folA* gives against this 60-base
reference — scores about 4e-10 against a 2 Mb genome, so a 1e-10 cut-off
silently discards the one hit that matters. A region differing from the
reference at more than 40% of its residues is reported as not determined, on
the grounds that the extraction, not the allele, is wrong.

---

## 3. Limits of working from an assembly

These are stated in the calls themselves, not only here.

**Reported as a pair or group, always.** The discriminating lesion exists but
the module does not determine whether it abolishes function, so both members
are reported: `9A/9V` (*wcjE* delG721), `18B/18C` (*wciX* TGA at codon
168–171), `35B/35D` (*wciG*), `7A/7F` (*wcwD* insT587). Of these, 18B/18C is
the most assembly-tractable because the lesion is a substitution creating a
stop rather than an indel — and whole-locus coverage can *never* separate them,
since the 18B and 18C reference loci are the same length. No pure 7A has yet
been resolved from sequence by PneumoCaT either, on the few isolates tried.

**15B / 15C: reported as `15B/15C`, always.** This is the one place the module
departs from the CDC rule, and it departs towards honesty. The two differ only
in the number of TA dinucleotide repeats in the O-acetyltransferase *wciZ* —
eight in 15B, nine in 15C, which frameshifts the gene — and the count drifts by
slipped-strand mispairing. Reads show the mixture directly (PneumoCaT routinely
finds 20–80% of reads carrying the frameshift and reports `15B/C`); an assembly
collapses a short tandem repeat to one consensus, which is the variant class
short-read assemblers resolve worst. The two are also very hard to separate
with SSI sera on the bench, and `15B/C` is an accepted reported serotype.
Validation: a strain whose BioSample serotype attribute is `15B` is reported
`15B/15C`, with an exact *wciZ* match recorded as the evidence.

**No published split exists.** `32A/32F` — no nucleotide difference has been
found within the coding regions (on five isolates), and neither PneumoCaT nor
PneumoKITy attempts it. `24A/24B/24F` — the loci do differ, but no validated
rule exists; PneumoCaT ships no serogroup 24 variant data and SeroBA reports
serogroup 24. Both are reported as the group with the reason in
`alleles["group_reason"]`.

**Not separated by the probe set.** `11B/11C`, `12A/12B/46`, `25A/25F`,
`35C/42`. For 25F vs 25A, PneumoCaT does publish six *wcyC* SNPs, but
PneumoKITy — UKHSA's current screening tool — dropped them, so the split is
treated as unavailable rather than implemented on one tool's authority.

**Silent relatives.** A serotype **11E** isolate is reported as `11A`, and the
note says so. 11E is 11A with *wcjE* disrupted, but every 11E described carries
a different, independently arisen lesion — IS*1515* insertions in either
orientation, deletions of 332 bp and 43 bp, a two-guanine insertion, a point
mutation — so there is no position to query, no tool implements it, and
O-acetylation is lost only partially in some isolates, so the boundary is not
clean. Likewise a **37** isolate carries a 33F-like locus and is caught only by
the separate *tts* probe, which lies outside the capsule locus.

**Stage-one-only serotypes** (16A/16F, 17A/17F, 33C, 35B/35F, 47A/47F and
anything else no probe rule covers) rest entirely on whole-locus identity and
degrade on fragmented assemblies. The note says "treat it as a serogroup-level
call confirmed to the locus, not a serologically verified serotype".

**PBP and MIC.** A mosaic region is `NEW` rather than a number, and then no MIC
is predicted. Half the validation panel got no MIC. This is a real limitation of
shipping only the published lookup table; see §7.

**Not attempted at all.** GPSC lineage assignment, null-capsule clade typing,
*pspA*/*pspC* typing, *rrgB* clade typing, *folP* sulfonamide insertions, 23S
rRNA and *rplV* substitutions. §7 says why for each.

---

## 4. Databases, provenance and licence

`fetch()` writes everything under `<db_dir>/spneumoniae/`, builds the BLAST
indexes, and records the provenance in `SOURCES.tsv` beside the data. Total
2.9 MB of reference data, 4.9 MB with indexes.

| Name | Files | Source | Version | Licence | Redistribution |
|---|---|---|---|---|---|
| `cps` | `cps/cps_loci.fna` (106 loci, 2.1 MB) | <https://github.com/GlobalPneumoSeq/cps_reference_database> (`all.fasta`) | `master`, retrieved 2026-11-17; repository last pushed 2025-08-13 | **none stated** — no LICENSE file and no licence statement in the README | **unstated** |
| `cdc-serotype-probes` | `cps/sero_probes.fna` (72 probes in 61 clusters, 14 KB) | <https://github.com/BenJamesMetcalf/Spn_Scripts_Reference> (`SPN_Reference_DB/SPN_Sero_Gene-DB_Final.fasta`) | `master`, retrieved 2026-11-17; repository last pushed 2026-08-31 | **none stated** — no LICENSE file. It is the reference data of the CDC Streptococcus Laboratory; works of the United States federal government are generally not subject to domestic copyright (17 U.S.C. §105), but the repository itself makes no such statement. | **unstated** |
| `cdc-pbp` | `pbp/pbp_regions.fna`, `pbp/pbp_1A.faa`, `pbp/pbp_2B.faa`, `pbp/pbp_2X.faa`, `pbp/pbp_type_mic.csv` (0.82 MB) | same repository | `master`, retrieved 2026-11-17; PBP allele sets as of the 2026-08-31 database update | none stated (as above) | **unstated** |
| `cdc-resistance-probes` | `res/res_probes.fna` (19 probes, 3.2 KB) | same repository (`SPN_Res_Gene-DB_Final.fasta`) | `master`, retrieved 2026-11-17 | none stated (as above) | **unstated** |

**For the integrator.** Neither upstream states a licence, so neither can be
assumed redistributable and both are downloaded by `fetch()` rather than
vendored. Two things would change that: a licence statement from the GPS
project for `cps_reference_database`, and a written confirmation from the CDC
Streptococcus Laboratory (or reliance on 17 U.S.C. §105) for
`Spn_Scripts_Reference`. The `cps` database is the only one over a megabyte,
so if only one can be shipped, ship the CDC files (0.84 MB total) and download
the *cps* loci.

The *cps* locus sequences themselves were produced by Bentley et al. 2006 and
are in GenBank; the GPS repository is the curated per-serotype rearrangement of
them. A licence-free fallback would be to fetch the GenBank records directly,
at the cost of redoing the curation.

Both repositories are pinned to `master` rather than a commit, which is a
weakness: a silent upstream change would alter calls. The commit identifiers are
in the module's `_CPS_COMMIT` and `_CDC_COMMIT` constants and should be set to
fixed SHAs once the integrator decides the release cadence.

---

## 5. Validation

20 public assemblies. `HYDRA_VALIDATION=1 pytest
tests/test_typing_spneumoniae.py` downloads them and asserts the calls.

Two kinds of expectation source are used, and they are not equally strong:

* **(P)** the published characterisation of the strain, cited below;
* **(B)** the serotype as submitted to NCBI, either as a BioSample `serotype`
  attribute or as the Assembly record's strain field. This is the submitting
  laboratory's own declaration, not a peer-reviewed result.

| Accession | Strain | Expected serotype | Src | Obtained | PBP type | Predicted penicillin MIC (meningitis / non-meningitis) |
|---|---|---|---|---|---|---|
| GCF_000006885.1 | TIGR4 | 4 | P | **4** ✓ | 0:0:3 | 0.03 — S / S |
| GCF_000014365.1 | D39 | 2 | P | **2** ✓ | 2:0:NEW | not determined (novel pbp2x) |
| GCF_000007045.1 | R6 | non-encapsulated (capsule-locus deletion derivative of D39) | P | **non-encapsulated**, "disrupted serogroup 2 locus", gap 2558–10062 ✓ | 196:0:2 | not determined (triplet absent from the table) |
| GCF_000026665.1 | ATCC 700669 / Spain23F-1 (PMEN1) | 23F | P | **23F** ✓ | 15:12:18 | 4 — **R** / I |
| GCF_000147095.1 | 670-6B | 6B | P | **6B** ✓ (locus matched the 6E reference; the *wciP*/*wciN* rule gave 6B) | 34:57:56 | 1 — **R** / S |
| GCF_000019265.1 | Hungary19A-6 | 19A | P | **19A** ✓ | NEW:348:620 | not determined (mosaic pbp1a) |
| GCF_000019025.1 | Taiwan19F-14 | 19F | P | **19F** ✓ | 13:16:47 | 2 — **R** / S |
| GCF_000019985.1 | CGSP14 | 14 | P | **14** ✓ | 15:31:18 | not determined (triplet absent from the table) |
| GCF_000210955.1 | OXC141 | 3 | P | **3** ✓ | 2:3:2 | 0.03 — S / S |
| GCF_000019005.1 | P1031 | 1 | P | **1** ✓ | 2:44:0 | 0.03 — S / S |
| GCF_000018965.1 | 70585 | 5 | P | **5** ✓ | 23:0:0 | 0.03 — S / S |
| GCF_033096205.1 | PZ900700054 | 6A | B | **6A** ✓ (*wciP6AC* exact × *wciN6AB*) | 13:31:73 | 2 — **R** / S |
| GCF_033096125.1 | PZ900700012 | 6C | B | **6C** ✓ (*wciP6AC* exact × *wciN6CD*; locus matched the 6D reference) | 2:673:1066 | not determined |
| GCF_022558605.1 | PZ900700608 | 15B | B | **15B/15C** — deliberately broader, see §3 | 2:26:498 | not determined |
| GCF_033096045.1 | PZ900700005 | 9V | B | **9A/9V** — broader; *wcjE* 99.6% identical to the reference, so function is undetermined | 1:4:493 | not determined |
| GCF_036322935.1 | PF2 / TKY092604 | 35B/35D (the submitter's own unresolved pair) | B | **35B** — *narrower* than the expectation; *wciG* matched at 100%, which is the CDC criterion for 35B. See the note below. | 4:7:7 | 2 — **R** / S |
| GCF_033096325.1 | PZ900701114 | 24F | B | **24A/24B/24F** — broader by design, no published split (§3). Locus matched the 24C reference. | 17:15:22 | 0.5 — **R** / S |
| GCF_022069005.1 | TVO_1901948 | 8 | B | **8** ✓ | 3:6:NEW | not determined (novel pbp2x) |
| GCF_002813955.1 | Korean 11A isolate | 11A | B | **11A** ✓ (locus matched the 11E reference; the *wzy11A* probe gave 11A) | NEW:NEW:NEW | not determined (all three regions mosaic) |
| GCF_000180515.1 | SPNA45 | 3 | B | **3** ✓ | 62:0:2 | not determined (triplet absent from the table) |

**Coverage of the required cases.** Vaccine type: 4, 6B, 14, 19F, 23F (PCV7
onwards). Non-vaccine type: 24A/24B/24F and 9A/9V are in nothing as called;
8 and 11A are PCV20/PCV21/PPSV23 but not PCV13; 6C is in nothing and is
reported as cross-protected by PCV13/15/20. Hard-to-separate serogroup: all
four serogroup 6 branches are exercised (6A, 6B, 6C) plus 15B/15C, 9A/9V,
24A/24B/24F, 35B. Non-encapsulated: R6. Penicillin-susceptible: TIGR4, OXC141,
P1031, 70585 (predicted 0.03 mg/L, S). Penicillin-non-susceptible: ATCC 700669
(4 mg/L), Taiwan19F-14 (2 mg/L), 670-6B (1 mg/L), PZ900700054 (2 mg/L) — all
resistant by the meningitis criterion.

**Independent corroboration of the PBP chain.** The CDC coordinate reference
`Ref_PBP_3.faa` is an exact substring of the TIGR4 proteins, and its three
fragments correspond to alleles 1A-0, 2B-0 and 2X-3. The module typed TIGR4 as
**0:0:3** from the assembly, which is the same answer arrived at independently.

**Disagreements, investigated not tuned away.**

1. *GCF_036322935.1, expected 35B/35D, obtained 35B.* The submitter reported the
   pair; the module resolved it. The CDC rule is explicit — `WZY35B` present
   with `WCIG35B` matching exactly is 35B, and 35D is called only when *wciG*
   differs — and here *wciG* matched at 100% identity over 100% of the probe.
   The module is therefore applying the published rule correctly and is more
   specific than the submitted label. It is *not* claiming that the submitter
   was wrong: if the isolate is in fact 35D by an *wciG* lesion the probe does
   not cover, this would be a false 35B. The threshold was not changed, because
   doing so would also turn every true 35B into a pair.
2. *GCF_033096045.1, expected 9V, obtained 9A/9V.* `WCJE9V` matched at 99.6%
   identity over the full 1,035-base gene — four mismatches, not an obvious
   inactivation. The CDC rule treats any difference as grounds for the pair, and
   the module follows it. A refinement that checked ORF integrity (internal
   stop, frameshift) rather than exact identity would very likely resolve this
   to 9V; see §7. Reporting 9V on a near-exact *wciE* match instead would be
   choosing convenience over the published rule.
3. *GCF_022558605.1, expected 15B, obtained 15B/15C.* Deliberate, see §3.

**Validation gaps, stated plainly.**

* No natural non-encapsulated isolate was validated — R6 is a laboratory
  capsule-deletion derivative. The two read the same way here, but an NCC1/2/3
  isolate with a *pspK* or *aliC*/*aliD* cassette in place of *cps* has not been
  tested.
* The 6D branch of the serogroup 6 rule, and the 18B/18C, 7A/7F, 33A, 9L/9N,
  15A/15F, 22A and 7B/7C rules, are covered by unit tests on synthetic probe
  evidence but by no real genome.
* No measured MIC was available for any panel strain; the penicillin column
  above is the module's own prediction, checked against each strain's published
  **category** (susceptible vs non-susceptible), not against a measured value.
* `pilus_islets`, `macrolide_mechanism`, `other_determinants` and
  `target_substitutions` were not validated against an independent tool. Their
  results are internally consistent and biologically coherent — PI-1 in TIGR4
  (the strain the islet was described in), *erm(B)* in the three MLS-B clones,
  *mef* in Taiwan19F-14, *tet(M)* in every Tn*916*-carrying clone, *folA* I100L
  in exactly the six co-trimoxazole-resistant multidrug clones and in none of
  the susceptible strains — but "coherent" is not "validated", and they should
  be checked against the CDC pipeline's own output before being relied on.

---

## 6. Bibliography

Every reference below was verified against PubMed at the time of writing: the
PMID resolves and the first author, title, journal and year match. They are also
in the module's `REFERENCES`, each with the call it supports.

**Capsular serotype — scheme, nomenclature and data**

1. Bentley SD, Aanensen DM, Mavroidi A, et al. Genetic analysis of the capsular
   biosynthetic locus from all 90 pneumococcal serotypes. *PLoS Genetics* 2006;
   2:e31. doi:10.1371/journal.pgen.0020031. PMID 16532061.
   — *the cps locus reference set stage one measures coverage against.*
2. Geno KA, Gilbert GL, Song JY, et al. Pneumococcal capsules and their types:
   past, present, and future. *Clinical Microbiology Reviews* 2015;
   28:871–899. doi:10.1128/CMR.00024-15. PMID 26085553.
   — *the serotype/serogroup nomenclature, and which within-serogroup pairs
   differ only by O-acetylation or by one glycosyltransferase.*
3. Metcalf BJ, Gertz RE Jr, Gladstone RA, et al. Strain features and
   distributions in pneumococci from children with invasive disease before and
   after 13-valent conjugate vaccine implementation in the USA. *Clinical
   Microbiology and Infection* 2016; 22:60.e9–60.e29.
   doi:10.1016/j.cmi.2015.08.027. PMID 26363404.
   — *the CDC whole-genome serotyping probe set, the pipeline it is used in, and
   the extraction of the PBP transpeptidase regions.*
4. Kapatai G, Sheppard CL, Al-Shahib A, et al. Whole genome sequencing of
   *Streptococcus pneumoniae*: development, evaluation and verification of
   targets for serogroup and serotype prediction using an automated pipeline.
   *PeerJ* 2016; 4:e2477. doi:10.7717/peerj.2477. PMID 27672516.
   — *which serotype distinctions sequence can and cannot make, including
   32A/32F, serogroup 24, 15B/15C and 7A/7F.*
5. Epping L, van Tonder AJ, Gladstone RA, et al. SeroBA: rapid high-throughput
   serotyping of *Streptococcus pneumoniae* from whole genome sequence data.
   *Microbial Genomics* 2018; 4:e000186. doi:10.1099/mgen.0.000186.
   PMID 29870330.
   — *the serogroup 24 limitation, and the reads-based tool this module is the
   assembly-only counterpart of.*

**Capsular serotype — the within-serogroup rules**

6. Mavroidi A, Godoy D, Aanensen DM, Robinson DA, Hollingshead SK, Spratt BG.
   Evolutionary genetics of the capsular locus of serogroup 6 pneumococci.
   *Journal of Bacteriology* 2004; 186:8181–8192.
   doi:10.1128/JB.186.24.8181-8192.2004. PMID 15576766.
   — *wciP as the 6A/6B determinant (Ser195 in 6A/6C, Asn195 in 6B/6D).*
7. Park IH, Pritchard DG, Cartee R, Brandao A, Brandileone MC, Nahm MH.
   Discovery of a new capsular serotype (6C) within serogroup 6 of
   *Streptococcus pneumoniae*. *Journal of Clinical Microbiology* 2007;
   45:1225–1233. doi:10.1128/JCM.02199-06. PMID 17267625.
   — *serotype 6C as a serotype.*
8. Park IH, Park S, Hollingshead SK, Nahm MH. Genetic basis for the new
   pneumococcal serotype, 6C. *Infection and Immunity* 2007; 75:4482–4489.
   doi:10.1128/IAI.00510-07. PMID 17576753.
   — *wciN as the 6A/6C determinant, and therefore the four-way rule.*
9. Bratcher PE, Kim KH, Kang JH, Hong JY, Nahm MH. Identification of natural
   pneumococcal isolates expressing serotype 6D by genetic, biochemical and
   serological characterization. *Microbiology* 2010; 156:555–560.
   doi:10.1099/mic.0.034116-0. PMID 19942663.
   — *serotype 6D, completing the wciP × wciN rule.*
10. van Selm S, van Cann LM, Kolkman MA, van der Zeijst BA, van Putten JP.
    Genetic basis for the structural difference between *Streptococcus
    pneumoniae* serotype 15B and 15C capsular polysaccharides. *Infection and
    Immunity* 2003; 71:6192–6198. doi:10.1128/IAI.71.11.6192-6198.2003.
    PMID 14573636.
    — *the wciZ TA-repeat count (eight units in 15B, nine and a frameshift in
    15C), and therefore why 15B and 15C are reported as a pair. The paper calls
    the gene cps15bM/cps15cM; the name wciZ comes from Bentley 2006.*
11. Calix JJ, Nahm MH. A new pneumococcal serotype, 11E, has a variably
    inactivated *wcjE* gene. *Journal of Infectious Diseases* 2010; 202:29–38.
    doi:10.1086/653123. PMID 20507232.
    — *wcjE inactivation as the mechanism behind the O-acetylation pairs, and
    the heterogeneity that makes 11E unpredictable from sequence.*
12. Calix JJ, Porambo RJ, Brady AM, et al. Biochemical, genetic, and
    serological characterization of two capsule subtypes among *Streptococcus
    pneumoniae* serotype 20 strains: discovery of a new pneumococcal serotype.
    *Journal of Biological Chemistry* 2012; 287:27885–27894.
    doi:10.1074/jbc.M112.380451. PMID 22736767.
    — *the split of Danish type 20 into 20A and 20B, and therefore how a PPSV23
    "type 20" membership is reported.*

**Non-encapsulated pneumococci**

13. Park IH, Kim KH, Andrade AL, Briles DE, McDaniel LS, Nahm MH. Nontypeable
    pneumococci can be divided into multiple *cps* types, including one type
    expressing the novel gene *pspK*. *mBio* 2012; 3:e00035-12.
    doi:10.1128/mBio.00035-12. PMID 22532557.
    — *the null-capsule clade definitions quoted when a non-encapsulated
    genotype is reported: NCC1 is pspK-positive and aliC/aliD-negative, NCC2
    carries aliC and aliD, NCC3 carries aliD only.*
14. Keller LE, Robinson DA, McDaniel LS. Nonencapsulated *Streptococcus
    pneumoniae*: emergence and pathogenesis. *mBio* 2016; 7:e01792-15.
    doi:10.1128/mBio.01792-15. PMID 27006456.
    — *what a non-encapsulated genotype means clinically.*

**Vaccines**

15. Cooper D, Yu X, Sidhu M, Nahm MH, Fernsten P, Jansen KU. The 13-valent
    pneumococcal conjugate vaccine (PCV13) elicits cross-functional
    opsonophagocytic killing responses in humans to *Streptococcus pneumoniae*
    serotypes 6C and 7A. *Vaccine* 2011; 29:7207–7211.
    doi:10.1016/j.vaccine.2011.06.056. PMID 21689707.
    — *the 6A-to-6C and 7F-to-7A cross-reactivity reported alongside
    membership.*

The formulations themselves were read from the products' own regulatory
documents on 2026-11-17, not from a review:

* **PCV15 (VAXNEUVANCE)** — FDA package insert, <https://www.fda.gov/media/150819/download>:
  serotypes 1, 3, 4, 5, 6A, 6B, 7F, 9V, 14, 18C, 19A, 19F, 22F, 23F, 33F.
* **PCV20 (PREVNAR 20)** — FDA package insert, <https://www.fda.gov/media/149987/download>:
  serotypes 1, 3, 4, 5, 6A, 6B, 7F, 8, 9V, 10A, 11A, 12F, 14, 15B, 18C, 19A,
  19F, 22F, 23F, 33F. The same insert names the seven original PCV7 serotypes
  verbatim as its otitis media indication: 4, 6B, 9V, 14, 18C, 19F, 23F.
* **PCV21 (CAPVAXIVE / V116)** — FDA Summary Basis for Regulatory Action,
  17 June 2024, <https://www.fda.gov/media/180070/download>: serotypes 3, 6A,
  7F, 8, 9N, 10A, 11A, 12F, 15A, 15B, 15C, 16F, 17F, 19A, 20A, 22F, 23A, 23B,
  24F, 31, 33F, 35B — 22 serotype names for 21 valencies, because 15B and 15C
  share one component.
* **PPSV23 (PNEUMOVAX 23)** — FDA package insert,
  <https://www.fda.gov/media/80547/download>: "the 23 serotypes contained in
  the vaccine (1, 2, 3, 4, 5, 6B, 7F, 8, 9N, 9V, 10A, 11A, 12F, 14, 15B, 17F,
  18C, 19F, 19A, 20, 22F, 23F, and 33F)".
* **PCV7, PCV10 (GSK/Synflorix) and PCV13** — taken from the Global
  Pneumococcal Sequencing project's curated valency table,
  <https://github.com/GlobalPneumoSeq/gps-database-processor/blob/HEAD/data/pcv_valency.csv>,
  cross-checked against the two inserts above (its PCV15 and PCV20 rows match
  them exactly, and PCV13 is the Vaxneuvance list minus 22F and 33F).
  PCV10 (Pneumosil, Serum Institute of India) is **deliberately omitted**: that
  table's row for it lists nine serotypes for a ten-valent product, and the
  discrepancy could not be resolved against a regulatory source in time.

**PBP typing and the MIC prediction**

16. Li Y, Metcalf BJ, Chochua S, et al. Penicillin-binding protein
    transpeptidase signatures for tracking and predicting beta-lactam
    resistance levels in *Streptococcus pneumoniae*. *mBio* 2016; 7:e00756-16.
    doi:10.1128/mBio.00756-16. PMID 27302760.
    — *the pbp1a/pbp2b/pbp2x transpeptidase allele triplet, its numbering, the
    PBP-type-to-MIC reference table, the 27-position resistance list (Table S1),
    and the published accuracy figures.*
17. Li Y, Metcalf BJ, Chochua S, et al. Validation of beta-lactam minimum
    inhibitory concentration predictions for pneumococcal isolates with newly
    encountered penicillin binding protein (PBP) sequences. *BMC Genomics*
    2017; 18:621. doi:10.1186/s12864-017-4017-7. PMID 28810827.
    — *the behaviour of the prediction for PBP types absent from the reference
    table, which is why such a triplet is reported without an MIC rather than
    extrapolated; and the independent validation figures. This is also the
    paper whose leave-one-type-out results are the evidence, quoted in §2.4a,
    for not building a nearest-neighbour mode-MIC model: on new PBP types that
    model showed essential agreement below 90 % for amoxicillin and cefotaxime
    and unacceptably high very-major-discrepancy rates for 4 of the 6
    beta-lactams, and the paper's conclusion endorses the random forest, which
    is not implementable under this project's rules.*
18. Hung HCH, Kumar N, Dyster V, et al. GPS Pipeline: portable, scalable
    genomic pipeline for *Streptococcus pneumoniae* surveillance from Global
    Pneumococcal Sequencing Project. *Nature Communications* 2025; 16:8345.
    doi:10.1038/s41467-025-64018-5. PMID 40993121.
    — *the current production implementation of the CDC PBP MIC prediction, and
    the largest independent evaluation of it (1,273 genomes with broth-dilution
    MICs), in which penicillin was the weakest of the six drugs; context in
    §2.4a for how far the prediction travels outside US surveillance.*
19. Eriksen HB, Fuursted K, Jensen A, et al. Predicting beta-lactam
    susceptibility from the genome of *Streptococcus pneumoniae* and other
    mitis group streptococci. *Frontiers in Microbiology* 2023; 14:1120023.
    doi:10.3389/fmicb.2023.1120023. PMID 36937294.
    — *an independent, non-US evaluation reporting that 37.5 % of 88 Danish
    pneumococcal isolates carried no recognised PBP profile, which is why
    declining to predict fires often outside US invasive surveillance (§2.4a).*
20. Hakenbeck R, Brückner R, Denapaite D, Maurer P. Molecular mechanisms of
    beta-lactam resistance in *Streptococcus pneumoniae*. *Future Microbiology*
    2012; 7:395–410. doi:10.2217/fmb.12.2. PMID 22393892.
    — *the SXXK, SXN and KTG/KSG active-site motif definitions and positions
    the CDC scheme adopts.*
21. Smith AM, Klugman KP. Alterations in PBP 1A essential for high-level
    penicillin resistance in *Streptococcus pneumoniae*. *Antimicrobial Agents
    and Chemotherapy* 1998; 42:1329–1333. doi:10.1128/AAC.42.6.1329.
    PMID 9624469.
    — *pbp1a residue 371 and the TSQF574-577 block.*
22. Contreras-Martel C, Job V, Di Guilmi AM, Vernet T, Dideberg O, Dessen A.
    Crystal structure of penicillin-binding protein 1a (PBP1a) reveals a
    mutational hotspot implicated in beta-lactam resistance in *Streptococcus
    pneumoniae*. *Journal of Molecular Biology* 2006; 355:684–696.
    doi:10.1016/j.jmb.2005.10.030. PMID 16316661.
    — *the structural basis of the pbp1a T371 / TSQF574-577 hotspot.*
23. Pagliero E, Chesnel L, Hopkins J, et al. Biochemical characterization of
    *Streptococcus pneumoniae* penicillin-binding protein 2b and its
    implication in beta-lactam resistance. *Antimicrobial Agents and
    Chemotherapy* 2004; 48:1848–1855. doi:10.1128/AAC.48.5.1848-1855.2004.
    PMID 15105143.
    — *pbp2b T446A, its proximity to the S443SN motif, and its 60% reduction in
    penicillin G affinity.*
24. Contreras-Martel C, Dahout-Gonzalez C, Martins A dos S, Kotnik M, Dessen A.
    PBP active site flexibility as the key mechanism for beta-lactam resistance
    in pneumococci. *Journal of Molecular Biology* 2009; 387:899–909.
    doi:10.1016/j.jmb.2009.02.024. PMID 19233207.
    — *the structural role of pbp2b T446 and E476.*
25. Gordon E, Mouz N, Duée E, Dideberg O. The crystal structure of the
    penicillin-binding protein 2x from *Streptococcus pneumoniae* and its
    acyl-enzyme form: implication in drug resistance. *Journal of Molecular
    Biology* 2000; 299:477–485. doi:10.1006/jmbi.2000.3740. PMID 10860753.
    — *pbp2x T338, T550 and Q552 as the residues directly related to
    resistance.*
26. Mouz N, Gordon E, Di Guilmi AM, et al. Identification of a structural
    determinant for resistance to beta-lactam antibiotics in Gram-positive
    bacteria. *Proceedings of the National Academy of Sciences USA* 1998;
    95:13403–13406. doi:10.1073/pnas.95.23.13403. PMID 9811812.
    — *the buried cavity beside the pbp2x catalytic serine that T338 lines.*
27. Chesnel L, Pernot L, Lemaire D, et al. The structural modifications induced
    by the M339F substitution in PBP2x from *Streptococcus pneumoniae* further
    decreases the susceptibility to beta-lactams of resistant strains. *Journal
    of Biological Chemistry* 2003; 278:44448–44456. doi:10.1074/jbc.M305948200.
    PMID 12923202.
    — *pbp2x M339F.*
28. Smith AM, Klugman KP. Amino acid mutations essential to production of an
    altered PBP 2X conferring high-level beta-lactam resistance in a clinical
    isolate of *Streptococcus pneumoniae*. *Antimicrobial Agents and
    Chemotherapy* 2005; 49:4622–4627. doi:10.1128/AAC.49.11.4622-4627.2005.
    PMID 16251304.
    — *pbp2x T338P, I371T and R384G as the substitutions whose reversion
    restores susceptibility.*
29. Clinical and Laboratory Standards Institute. Performance standards for
    antimicrobial susceptibility testing; twenty-third informational supplement
    (M100-S23). CLSI, Wayne PA, 2013.
    — *the breakpoints the predicted MIC is interpreted against, including the
    separate meningitis criteria. This is the edition both CDC papers cite; the
    thresholds are unchanged in current M100.*

**Resistance determinants and target substitutions**

30. Metcalf BJ, Chochua S, Gertz RE Jr, et al. Using whole genome sequencing to
    identify resistance determinants and predict antimicrobial resistance
    phenotypes for year 2015 invasive pneumococcal disease isolates recovered in
    the United States. *Clinical Microbiology and Infection* 2016;
    22:1002.e1–1002.e8. doi:10.1016/j.cmi.2016.08.001. PMID 27542334.
    — *the macrolide, lincosamide, tetracycline, fluoroquinolone and
    folate-pathway determinant set, the reference peptides and offsets the
    substitutions are numbered on, and the determinant-to-phenotype
    interpretation.*

**Pilus islets**

31. Barocchi MA, Ries J, Zogaj X, et al. A pneumococcal pilus influences
    virulence and host inflammatory responses. *Proceedings of the National
    Academy of Sciences USA* 2006; 103:2857–2862. doi:10.1073/pnas.0511017103.
    PMID 16481624.
    — *pilus islet 1 (the rlrA islet, rrgA–rrgC) and its clinical relevance.*
32. Bagnoli F, Moschioni M, Donati C, et al. A second pilus type in
    *Streptococcus pneumoniae* is prevalent in emerging serotypes and mediates
    adhesion to host cells. *Journal of Bacteriology* 2008; 190:5480–5492.
    doi:10.1128/JB.00384-08. PMID 18515415.
    — *pilus islet 2 (pitA–pitB–sipA–srtG1) and the pitB marker.*

**Lineage (cited for what is *not* done)**

33. Gladstone RA, Lo SW, Lees JA, et al. International genomic definition of
    pneumococcal lineages, to contextualise disease, antibiotic resistance and
    vaccine impact. *EBioMedicine* 2019; 43:338–346.
    doi:10.1016/j.ebiom.2019.04.021. PMID 31003929.
    — *the GPSC scheme, which this module explains it does not assign and why.*

---

## 7. Considered and left out

**GPSC (Global Pneumococcal Sequence Cluster) — not implemented.** GPSCs are
defined by PopPUNK clustering against the GPS reference database, and PopPUNK
needs numpy, scipy and pp-sketchlib, all compiled, which the frozen Windows app
rules out. The downloadable
`GPS_v9_external_clusters.csv` (716 KB,
<https://gps-project.cog.sanger.ac.uk/GPS_v9_external_clusters.csv>) maps the
GPS project's internal sample identifiers to GPSC numbers and is useless without
the clustering step. No authoritative ST→GPSC table is published; the
relationship is many-to-many and a lookup would misassign recombinant isolates
in exactly the lineages that matter. **Owner's decision needed**: a documented,
explicitly approximate ST→GPSC hint could be built from the GPS metadata, or
GPSC could be left to the GPS project's own pipeline.

**Null-capsule clade typing (NCC1/2/3) — not implemented, and the single most
valuable addition.** It needs three gene references; verified GenBank
accessions exist and are small: *pspK* MW205001 (1,923 bp) and MW205002;
*aliC* MW205003/4/5 (1,968 bp); *aliD* MW205006 (truncated allele) and
MW205007 (1,959 bp). NCC1 is *pspK*-positive, NCC2 carries *aliC* and *aliD*,
NCC3 carries *aliD* only (Park 2012). Adding the flanking *dexB* and *aliA*
anchors would also distinguish the cassette architecture from a plain deletion.
Left out because no validated non-encapsulated clinical genome was to hand and
the brief forbids shipping an unvalidated call.

**ORF-integrity checking for the O-acetylation pairs — the second most valuable
addition.** 9A/9V, 18B/18C, 35B/35D, 7A/7F and 33A/33F all turn on whether a
named gene is intact, and the CDC rule approximates that with "is the probe an
exact match". The full-length genes are already in the probe set
(`WCJE9V` 1,035 bp, `WCIX18C` 1,002 bp, `WCIG35B` 996 bp, `WCWD7F` 963 bp,
`WCIZ15F` 978 bp), so the region could be extracted, translated and tested for
a premature stop or a frameshift — the machinery for that is already in this
module. That would turn four of the pair calls into definite calls with real
evidence, and would very likely have resolved the 9V genome in §5 correctly.
Not done because none of 9A, 18B, 35D or 7A was available to validate against.

**The mode-MIC nearest-neighbour model — looked up in the literature, and NOT
built.** This was listed here as the next addition, on the reasoning that it
would roughly double MIC coverage. The evidence was then checked, and it does
not support building it: the published leave-one-type-out validation of exactly
this model on exactly these isolates shows essential agreement below 90 % for
amoxicillin and cefotaxime and *"unacceptably high vmj rates for 4 of the 6
antibiotics, with 3 of them higher than 50%"* — resistant isolates predicted as
susceptible (Li 2017). The model the authors endorse for new PBP types is the
random forest, which cannot be built under this project's rules: it ships as
~210 MB of serialized R objects with a pinned `randomForest` release because
the predictions depend on it. **Nothing was implemented.** The full evidence,
with the quotations and the FDA criteria, is in §2.4a, so the decision can be
re-examined rather than rediscovered.

**Deliberately out of scope**

* *folP* sulfonamide resistance: the lesion is a 1–2 codon insertion at
  residues 57–70 whose position varies, and the CDC pipeline handles it with a
  separate extraction path. *folA* I100L alone predicts trimethoprim
  resistance, which is reported; co-trimoxazole resistance is not.
* 23S rRNA domain V and *rplV* (L22) substitutions: multi-copy (23S) or needing
  an extraction path not yet written. *rplD* (L4) is reported.
* *pspA* family/clade and *pspC* group typing: real schemes, but highly
  repetitive loci that assemble badly, and reference laboratories do not report
  them for clinical purposes.
* *rrgB* clade typing (I/II/III): would need the clade reference alleles, and
  the islet-presence call already carries the clinically relevant information.
* Prophage content: pneumococcal prophages are abundant but, unlike in
  *S. pyogenes*, do not carry the main virulence or resistance determinants, so
  there is no established clinical report built on them.
* Tn*916*-family composite element typing (Tn*2009*, Tn*2010*, Tn*2017*, Mega
  class I/II): would explain *why erm(B)* and *tet(M)* travel together, which is
  epidemiologically useful and clinically not actionable. The constituent genes
  are already reported.
