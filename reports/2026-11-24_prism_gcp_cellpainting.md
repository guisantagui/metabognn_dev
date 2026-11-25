# Report 10 — PRISM, Panorama GCP and the Way cell-painting set: what each can actually test

*2026-11-24. Assessment of `data/validation/{prism,panorama_gcp,way}/`, downloaded 2026-11-23/24 with the
rationales in each `notes.md`. No association has been run yet; this is a feasibility and design review.
Continues report 09, whose section 6 asked for a second (phenotype, cell line) pair and whose section 8
disqualified naia for that purpose.*

**Summary.** **PRISM is the valuable one**, and not for the design its rationale proposes: its phenotype is
viability, which is the confounder that has dominated reports 07-09, so a direct "does the metabolic score
predict growth inhibition" test is close to uninterpretable. What PRISM makes possible instead is a
**within-drug, between-cell-line differential** design that cancels the confounder by construction, over
**821 drugs in 5 LINCS core lines**. **Panorama GCP has the most elegant rationale of anything we have and a
compound library that largely defeats it** - 117 compounds, dominated by direct inhibitors of the very enzymes
whose marks are being read. **The Way data is not actually downloaded**: all 20 profile files are 133-byte
Git LFS pointers.

## 1. PRISM — 1,449 compounds x 482 cell lines, and HEPG2 is missing

`secondary-screen-dose-response-curve-parameters.csv`: 701,004 rows, **1,449 distinct BRD compounds**, **482
cell lines**, 8-step 4-fold dilution from 10 uM, with `auc`, `ic50`, `ec50`, `slope`, `r2` per drug x line,
plus `moa`, `target`, `smiles`, `phase`.

Mapping is **exact and needs no structure work**: `broad_id` is a BRD id (strip to the first two segments) and
`ccle_name` names the line. **1,103 of 1,449 compounds are in `lincs_unique_perts.tsv`.**

| LINCS core | CCLE name | in PRISM | drugs with an AUC | also in LINCS |
|---|---|---|---|---|
| A375 | A375_SKIN | yes | 1,302 | 988 |
| A549 | A549_LUNG | yes | 1,378 | **1,053** |
| HT29 | HT29_LARGE_INTESTINE | yes | 1,364 | 1,045 |
| MCF7 | MCF7_BREAST | yes | 1,345 | 1,043 |
| PC3 | PC3_PROSTATE | yes | 1,331 | 1,008 |
| **HEPG2** | HEPG2_LIVER | **ABSENT** | - | - |
| VCAP | VCAP_PROSTATE | **ABSENT** | - | - |
| HCC515 | HCC515_LUNG | **ABSENT** | - | - |
| HA1E | - | absent by construction (immortalised normal, not a CCLE cancer line) | - | - |

**HEPG2's absence is the headline limitation.** It is checked, not assumed: PRISM carries 14 liver lines
(SKHEP1, HUH7, HEP3B217, JHH4/5/6, ...) and HEPG2 is not among them. Every positive result in reports 07-09
lives in HepG2, so **PRISM cannot extend them directly** - it can only test whether the same machinery works
in other lines.

### 1.1 Why the obvious design does not work

The rationale in `notes.md` is "predicted collapse of nucleotide/amino-acid production capacity predicts
antiproliferative effect per cell line". The problem is that **the phenotype is viability, and viability is
the confounder**. Reports 07-09 spent most of their effort showing that associations between MetChange scores
and a phenotype survive controlling for responsiveness and cell death - report 09, 9.6 turns on exactly this,
with expression's association dying under a cytotoxicity match while MetChange's survived. **You cannot match
on cytotoxicity when cytotoxicity is the phenotype.**

Concretely: any compound that is generally bioactive lowers AUC, and generally bioactive compounds produce
large MetChange signatures in every set. A naive test would likely return a significant association for
*every* metabolite set, including the glycan negative controls, and the size-matched random-set control would
be the only thing standing - which is exactly the situation report 08, 3.1 showed we handle badly.

### 1.2 The design that does work: within-drug, between-cell differences

PRISM's value is the **cell-line resolution**, which permits a contrast no previous dataset allowed. For each
drug *d* and each pair of cell lines *X*, *Y*:

```
phenotype difference:  AUC(d, X) - AUC(d, Y)
predictor difference:  setscore(d | X labels) - setscore(d | Y labels)
```

and ask whether the two correlate across drugs. **Differencing removes the drug's general potency** -
a compound that is simply toxic everywhere contributes zero to both sides. What remains is the *selectivity*
question: when a drug hits one line harder than another, do that line's MetChange labels say so?

> **CORRECTED after running it (section 8.3).** The original wording here said the confounder is "cancelled by
> construction rather than matched away". That is only half right. Differencing cancels confounding that is
> **shared** between the two cells; it does **not** cancel a **cell-specific** coupling between score and
> sensitivity. Those couplings turned out to range from +0.036 to +0.126 across the five lines, and the
> difference between them became the largest signal in the result. Read section 8 before relying on anything
> in 1.2.

This is the symmetric test report 09, section 6 asked for, and it is the first version of it with real power:

| | report 09, 9 (naia) | PRISM differential |
|---|---|---|
| unit | 25 case drugs | **821 drugs x 10 cell pairs** |
| confounder handling | matched on measured cytotoxicity, 15 pairs survive | **cancelled by differencing** |
| cell-matched? | no (SH-SY5Y labels impossible) | **yes, 5 lines** |

Coverage, computed on the 5 covered cores:

- **1,062 drugs have an AUC in all five lines; 821 of those are in LINCS.** (The common-drug restriction of
  report 09, section 2 is essential here and already learned: coverage differs per line.)
- AUC over those 821: median 0.91, IQR 0.81-1.07 (1.0 = no effect). 317 drugs reach AUC < 0.8 in at least one
  line, only 98 in all five - so most activity *is* selective rather than uniform.
- Between-cell spread per drug (SD of AUC across the five lines): median 0.090, IQR 0.046-0.142, max 0.46.
  **361 of 821 drugs have SD > 0.1**, and those are the informative ones: a drug with no between-cell spread
  carries no signal for this design whichever way the prediction goes.

**What it needs: MetChange labels for at least two PRISM-covered cores.** Locally only A549 has a
`metchange/rank/log2FC.csv` (HEPG2 has one too but HEPG2 is absent from PRISM). The HPC has all eight. Two
lines give one contrast; all five give ten.

**Pre-register the direction before running it.** The prediction is that reduced predicted production capacity
for biomass precursors - nucleotides, amino acids - in line *X* relative to *Y* goes with *lower* AUC (more
growth inhibition) in *X* relative to *Y*. Writing that down first is what report 09, 9.4 criticises itself
for not doing, and with ten cell pairs and ~160 set x variant x statistic combinations the post-hoc risk here
is larger, not smaller.

## 2. Panorama GCP — the best rationale we have, and a library that defeats it

28 plates, 49 histone-mark probes x 92 wells each, 24 h, doses 0.1-10 uM. `pert_id`, `pert_iname`,
`pert_dose`, `cell_id` and `pr_gcp_histone_mark` are all present, so **mapping is exact**.

**7 cell lines, and the same compounds in all of them:** MCF7, A549, A375, PC3 (LINCS cores) plus NPC, YAPC and
Astrocyte. Per-cell distinct compounds: 114-117. **That crossed structure is exactly what we have been
missing** - one compound set, several cell lines, a phenotype measured in each.

The rationale in `notes.md` is the strongest design logic in the project so far, because the **signs differ by
metabolite and are predicted in advance**: acetyl-CoA positive on acetylation, NAD+ negative on the same
marks, SAM positive on methylation, 2-OG negative on methylation while succinate/fumarate push it up. Opposite
signs on the same readout is a far harder test to pass by accident than any single correlation, and it is
immune to the "everything moves together" failure that viability phenotypes suffer.

### 2.1 The problem: 117 compounds, mostly direct epigenetic-enzyme inhibitors

The library is built around the enzymes whose marks the assay reads:

| class | compounds |
|---|---|
| HDAC inhibitors | vorinostat, trichostatin A, belinostat, MS-275, CI-994, RGFP966 |
| HMT inhibitors | EPZ-5676, EPZ004777, EPZ-5687, GSK126, CPI-169, UNC-0321, UNC-0646, BIX-01294, BIX-01338, UNC1215 |
| KDM inhibitors | GSK-J4, methylstat |
| BET inhibitors | JQ1-S, GSK525762A, GSK1210151A |
| HAT inhibitor | C646 |
| sirtuin inhibitors | EX527, salermide |
| DNMT inhibitors | decitabine, zebularine |

**The rationale requires acetylation to report acetyl-CoA availability. In this compound set it will
overwhelmingly report direct enzyme inhibition.** Vorinostat raises H3K9ac because it blocks HDACs, not
because acetyl-CoA rose. Run naively, the test would "succeed" for reasons that have nothing to do with
metabolism.

### 2.2 How to salvage it, and the two things that get better

1. **Exclude the direct epigenetic drugs and test the remainder.** That leaves roughly 70-80 compounds -
   mostly kinase inhibitors, but with genuinely metabolic members: **niclosamide** (protonophore uncoupler,
   should drop citrate-derived acetyl-CoA and so acetylation), **pravastatin** (mevalonate pathway),
   **bafilomycin A1** (V-ATPase), **resveratrol**, **curcumin**, **2-DG-adjacent** nothing, unfortunately.
   Small, but the crossed-sign structure means each compound yields several independent sign predictions
   rather than one, which buys back some of the lost power.
2. **The excluded drugs become a positive control for the assay**, not waste: the HDAC inhibitors *must*
   raise acetylation and the HMT inhibitors *must* lower the corresponding methylation. If they do not, the
   phenotype extraction is wrong and nothing else in the analysis is interpretable. No previous validation
   dataset gave us an internal assay control of that quality.
3. **The sirtuin inhibitors check the sign convention.** EX527 and salermide block NAD+-consuming
   deacetylases, so they should raise acetylation - the same direction predicted for *low* NAD+ production
   capacity. That is a built-in test of the NAD+ arm's sign, independent of the metabolic prediction.

Binding constraint is n: ~115 compounds per cell line before exclusions, ~70-80 after. Against that, four of
the seven lines are LINCS cores with labels, and the design yields multiple pre-specified signs per compound.
**Worth running after PRISM, and only with the epigenetic drugs split out as a separate positive-control
arm.**

## 3. Way — not downloaded

`broadinstitute-lincs-cell-painting-74f7a0d/` is the repository skeleton only. **All 20 profile files are
133-byte Git LFS pointer stubs**, for example:

```
version https://git-lfs.github.com/spec/v1
oid sha256:aa7c3f36874afceba2ef2d684dbd3ce653574bae9971e32639f74320e008bd90
size 17889995
```

`lincs-cell-painting-v1.zip` (110 MB) is the same commit and carries the same 133-byte stubs, so unzipping it
again will not help. **What is real and useful**: `metadata/moa/repurposing_info*.tsv`, `metadata/platemaps/`
and `metadata/clue_manifest/` - the MoA annotation tables, which are worth keeping regardless.

To get the data: `git lfs pull` inside that directory, or fetch the profile files directly from the Zenodo
record (5008187) rather than the GitHub tarball.

**Why it is worth getting.** The batch is `2016_04_01_a549_48hr_batch1` - **A549, 48 h, ~1,500 compounds** -
and A549 is the one non-HepG2 core line with local MetChange labels. That is a properly powered cell-matched
test with Mito-channel features, and unlike nyffeler (report 09, 8.3, 87 drugs) the overlap would be large.
The caveat from nyffeler carries over: cell painting is morphology, not metabolism, so only the `Mito_*`
features bear on the hypothesis, and a mito-selectivity contrast (mito features against the rest) is the
version that controls for general bioactivity internally.

## 4. Recommended order

1. **PRISM differential**, A549 vs one other core to start (needs that core's labels from scratch), direction
   pre-registered, restricted to the 821 common LINCS drugs and reporting the 361 cell-selective ones
   separately. This is the best-powered cell-matched test available and the confounder is cancelled rather
   than controlled.
2. **`git lfs pull` the Way profiles** - cheap, and it unlocks a 1,500-compound A549 cell-matched set.
3. **Panorama GCP**, with the epigenetic inhibitors split out as a positive-control arm and the crossed signs
   written down in advance.
4. ~~**xu** (report 09, 8.2) remains unrun and is still the cheapest falsification test of the MMP result.~~
   **WITHDRAWN - run, and underpowered by construction: see section 5.** Only 11 vs 12 compounds reach the
   HEPG2 panel, where even a HepG2-sized effect would give p ~ 0.21.

## 5. xu, run: the falsification test is underpowered by construction and cannot decide anything

*Built with `scripts/mc_phenotype_assoc/build_xu_phenotype.py` (new), scored against HEPG2 MetChange labels
with the unchanged `mc_phenotype_assoc.py`. Phenotype
`results/mc_phenotype_assoc/phenotypes_ALL/xu_mpp_ac16.tsv`, results
`results/mc_phenotype_assoc/HEPG2/xu_mpp/`.*

The phenotype extraction is clean. Of 247 compounds: 61 `active`, 123 `inactive`, 63 `inconclusive` (dropped).
**All 61 actives have a negative curve rank (median -3.33) and negative efficacy (median -53.6%)** - the
builder verifies this rather than assuming it - so the case group is directionally pure membrane-potential
*decrease*, exactly parallel to Attene-Ramos `MMP_down`.

**The mapping is what kills it.** Only **11 MMP_down and 12 inactive** compounds reach the HEPG2 LINCS panel,
and the ceiling is barely better:

| panel used | MMP_down | inactive |
|---|---|---|
| HEPG2 (5,494 compounds) | **11** | **12** |
| A375 (13,314 compounds) | 14 | 20 |
| full LINCS (34,419 perturbations) - a ceiling, not a runnable option | 17 | 25 |

What that n can detect, for a true effect the size of the HepG2 Attene result (AUROC 0.350):

| | n | p an AUROC of 0.350 would produce |
|---|---|---|
| HEPG2 panel (what we have) | 11/12 | **0.21** |
| A375 panel | 14/20 | 0.14 |
| full-LINCS ceiling | 17/25 | 0.098 |
| Attene in HEPG2, for scale | 144/475 | 4.7e-08 |

**So the test cannot reach significance even if the effect is exactly as large as the one it is meant to
check.** It is uninformative by construction, not by outcome.

The result, for the record:

| set | variant | AUROC | MW p | perm | random-set |
|---|---|---|---|---|---|
| tca | cyto | 0.295 | 0.10 | 0.10 | 0.16 |
| tca | mito_nohub | 0.311 | 0.13 | 0.13 | 0.18 |
| glycolysis | cyto | 0.379 | 0.34 | 0.34 | 0.44 |
| tca_endpoints | cyto_nohub | 0.402 | 0.44 | 0.45 | 0.46 |
| oxphos | mito | 0.470 | 0.83 | 0.84 | 0.87 |
| energy_currency | mito | 0.614 | 0.37 | 0.38 | 0.38 |

Nothing survives any control. **Read the direction of the non-result carefully, because it is not the
reassuring one.** A *positive* xu result would mean HepG2 labels predict membrane potential in a cell line
they know nothing about, i.e. the signal travels with the compounds. The point estimates lean weakly that way
(tca 0.295-0.311, the same pathway and a similar magnitude to HepG2's own 0.350), and `oxphos mito` is null
here (0.470) just as in HepG2 (0.471). At 11 against 12 this is a hint and nothing more - but it hints in the
direction that would further qualify report 09, section 3, not support it. It is also consistent with report
09, 9.3, which already refuted the strong "HepG2 labels are inert elsewhere" reading; the new wrinkle would be
that for the *MMP* phenotype the transfer runs through **TCA**, the same pathway as in HepG2, whereas the
SH-SY5Y ATP transfer ran through glycolysis.

### 5.1 Why the overlap is so bad, and a mistake to record

Xu screened **Tox21 environmental chemicals** - the actives include 4-nonylphenol (branched),
2-ethylhexylparaben and kadethrin: parabens, pesticides and industrial chemicals. LINCS is a drug-focused
library. This is the **same chemical-space mismatch that produced the Hallinger null** (report 07), and it was
predictable.

**I recommended xu as "the cheapest falsification test" in report 09, 8.2 and again in report 10, 4 on the
strength of a whole-LINCS overlap of 53 compounds, without checking the HEPG2-panel-specific count.** That is
the same error as the naia one recorded in report 09, 8.1 - 683 whole-LINCS matches collapsing to 10 in the
cell line that mattered - made again a day later. The rule that follows: **a validation dataset's usable n is
its overlap with the specific label panel, computed before the dataset is ranked, never its overlap with all
of LINCS.**

**Consequence for the plan of section 4:** item 4 is withdrawn. xu is done and cannot be improved by better
structure matching - `map_smiles.R` with salt stripping and stereo handling might add a handful of compounds
to 11/12, which does not change the power calculation. The falsification role it was meant to play falls to
the **PRISM differential** (section 1.2), which has 821 drugs and cancels the viability confounder rather than
leaving it uncontrolled as xu does - and note xu could not have controlled it anyway, since Xu reports **no
viability readout**, so unlike naia there is no measured cytotoxicity to match on.

## 6. Label generations verified, and how conserved a signature is between cell lines

*Prompted by a provenance scare on 2026-11-24. MetChange labels exist in more than one generation
distinguished only by date: **the 2026-11-12 files are the corrected ones**, the superseded buggy generation
was 2026-11-02. Nothing in a file's name or content identifies which it is, so the matrix mtime is stated with
every number below. A withdrawal issued earlier that day against the conservation result was mistaken - it had
used the correct Nov-12 files - and is reversed here.*

### 6.1 The section-3 A549 null reproduces exactly from the local corrected labels

Re-running the Attene-Ramos MMP association against
`results/labels/A549/metchange/rank` (`log2FC.csv`, 4,427,631,022 bytes, **mtime 2026-11-12 13:17**) and
comparing with the HPC cross-cell run used in report 09, section 3:

| set | variant | HPC run (report 09, 3) | local Nov-12 labels |
|---|---|---|---|
| tca | mito_nohub | 0.477 (p = 0.34) | **0.477 (p = 0.34)** |
| tca_endpoints | cyto_nohub | 0.471 (p = 0.22) | **0.471 (p = 0.22)** |
| tca | cyto | 0.461 (p = 0.10) | **0.461 (p = 0.10)** |
| oxphos | mito | 0.518 (p = 0.46) | **0.518 (p = 0.46)** |
| energy_currency | mito | 0.497 (p = 0.92) | **0.497 (p = 0.92)** |
| glycolysis | cyto | 0.463 (p = 0.12) | **0.463 (p = 0.12)** |

**All 40 shared set x variant combinations agree: median AUROC difference +0.0000, max |difference| 0.0000,
40 of 40 identical to three decimals.** Contrast sizes are identical too (179 cases / 803 controls).

Two consequences. **The A549 null of section 3 is confirmed** on labels known to be the corrected generation.
And because the numbers are bit-identical, **the HPC cross-cell runs were themselves computed on this same
generation**, so the set-level conservation table in 6.3 - which is derived from those runs' per-drug scores -
is also on corrected labels.

One discrepancy worth recording rather than hiding: the local matrix carries **9,830 perturbations / 1,085
drugs** against the HPC run's **10,150 / 1,118**, about 3% fewer. The missing columns belong to groups outside
this contrast (`MMP_up`, `cytotox`, `incon`), which is why the `MMP_down vs inactive` counts match exactly. The
two matrices are therefore not literally the same file, and anything that uses the other contrasts should be
run on one or the other, not compared across them.

### 6.2 The question: how much of a signature is the compound, and how much is the cell?

This arose from a correct objection: a compound's mechanism does not change between cell lines, and the TCA
cycle, glycolysis and OXPHOS are common to all of them, so **some transfer of a prediction between cell lines
is the expected baseline, not a finding.** That reframes reports 09 and 10:

- section 9.3 called the transfer of HEPG2 labels to an SH-SY5Y ATP phenotype a refutation of the "strong
  reading". It is better described as the expected baseline being observed.
- Report 10, section 5 read the weak xu lean as "the damaging direction". Same correction: a lean toward
  transfer is what shared mechanism predicts.
- **What then needs explaining is section 3** - why A549 labels fail to predict a HepG2 MMP phenotype
  (0.477, p = 0.34) if the compounds have the same mechanism in both cells.

The question is therefore quantitative, and measurable from the label matrices directly.
`scripts/mc_phenotype_assoc/signature_conservation.py` (new) computes it.

### 6.3 Conservation is real, weak, and set-dependent

**Full signatures**, HEPG2 vs A549, 250 shared drugs x 5,707 metabolites, perturbations averaged per drug:

| | median r | IQR |
|---|---|---|
| **same drug** | **+0.297** | -0.640 .. +0.808 |
| different drug (the control) | +0.031 | -0.707 .. +0.719 |

Same-vs-different **AUROC 0.551, p = 0.0055**; only **10%** of same-drug correlations exceed the different-drug
95th percentile. The median gap of +0.266 looks large but both distributions span nearly the whole range, so
**a drug's signature in one cell line is only marginally more predictive of its signature in another than a
random drug's is.** Read the gap and the AUROC, not the median: if some metabolites move under any
perturbation, every pair of signatures correlates.

**Set-level scores** (the quantity the association tests actually use), HEPG2 vs A549, from the cross-cell
runs' per-drug scores:

| set | members | r, all drugs | r, responsive only |
|---|---|---|---|
| energy_currency all | 69 | +0.179 | **+0.399** |
| glycolysis cyto_nohub | 27 | +0.037 | +0.474 |
| oxphos mito | 13 | +0.171 | +0.337 |
| tca mito_nohub | 19 | +0.020 | **+0.113** |
| tca all | 73 | +0.077 | +0.145 |

**Set size does not explain the differences.** Across all 40 set x variant combinations,
Spearman(n_members, conservation) = **-0.058 (p = 0.72)** over all drugs and **-0.191 (p = 0.24)**
responsive-only. So the attenuation explanation - small sets give noisy per-drug scores, which depresses
correlation - is ruled out, and TCA's low conservation is not a size artefact.

### 6.4 What this does and does not support

**Supported.** Conservation exists, which is why transfer happens; and it is **weak**, which makes section 3's
cell-specificity unsurprising rather than anomalous. A549's TCA score for a given drug correlates with HepG2's
at r = +0.11 to +0.15, so A549 labels carrying little information about a HepG2 TCA phenotype needs no appeal
to HepG2 being special. The depth explanation was already dead (report 09, 9.3: A375 has 2.4x HepG2's panel and
is the most null line); this supplies a positive account in its place.

**An earlier overstatement, withdrawn.** It was claimed in conversation that TCA is "not conserved at all" and
that there is a clean complementarity - what transfers does not predict, what predicts does not transfer. On
the best-powered pair, HEPG2 vs A549, **TCA conservation is positive (+0.11 to +0.32), merely the lowest of
the sets.** The negative values that prompted the stronger claim (-0.371 for HCC515, -0.154 for A375) come
from pairs with only 45-61 responsive drugs and cannot carry it. The ordering holds directionally; the
complementarity is a hypothesis, not a result.

**Not addressed.** Only HEPG2-vs-X pairs were computed, never X-vs-Y. Responsive-only n is 31-68 per pair. And
conservation of the *z* statistic is partly conservation of relative pathway positioning, since `z` is already
standardised within a perturbation - a different normalisation could give a different answer.

## 7. How to run the PRISM differential on the cluster

*Written so that **no large file has to move**. The PRISM download is ~800 MB and the label matrices are
GB-scale; the only things that travel are an 86 kB phenotype table going up and a ~3 MB scores file per cell
line coming back.*

### 7.1 The four pieces

| script | runs where | why |
|---|---|---|
| `build_prism_phenotype.py` | **laptop** (where PRISM is) | writes one 86 kB table, 821 drugs x 5 cells |
| `launch_prism_slurm.sh` | **cluster** | array 0-4, one task per PRISM-covered core line |
| `launch_prism_cell.sh` | cluster (called per task, or standalone with `cell=`) | scores one cell against its MetChange labels |
| `prism_differential.py` | **either** | the actual test; needs only the small files |

### 7.2 The three commands

```bash
# 1. laptop - build the table and send that one file
python3 -I scripts/mc_phenotype_assoc/build_prism_phenotype.py     --out results/mc_phenotype_assoc/phenotypes_ALL/prism_viability.tsv
scp results/mc_phenotype_assoc/phenotypes_ALL/prism_viability.tsv     <host>:/scratch/users/gsantamaria/test_large_files/metabognn_dev/results/mc_phenotype_assoc/phenotypes_ALL/

# 2. cluster - score every covered cell line (submit FROM scripts/mc_phenotype_assoc, see below)
cd scripts/mc_phenotype_assoc && sbatch launch_prism_slurm.sh
#    one line only:   cells="A549" sbatch --array=0-0 launch_prism_slurm.sh
#    no slurm:        cell=A549 bash scripts/mc_phenotype_assoc/launch_prism_cell.sh

# 3. either machine - the test
python3 -I scripts/mc_phenotype_assoc/prism_differential.py     --out results/mc_phenotype_assoc/cross_cell/prism_differential.tsv
```

Submit step 2 from `scripts/mc_phenotype_assoc`: Slurm expands `--output` before the job starts so it cannot
hold a variable, which is why `logs/` is relative to the submit directory and kept in git (report 08, 8.3).
A missing cell line fails only its own array task, so partial runs are fine - `prism_differential.py` uses
whichever cells it finds and names the ones it did not.

To bring the results home, download only
`results/mc_phenotype_assoc/<CELL>/prism_viability/per_drug_scores.tsv` - about 3 MB each, five of them.

### 7.3 What to read and what to ignore

Each run also writes `set_results.tsv`, the sensitive-vs-resistant group contrast. **That is the confounded
naive test and is not a result** (1.1): the phenotype is viability, viability is the confounder of reports
07-09, and nothing can be matched on cytotoxicity when cytotoxicity *is* the phenotype. Both launchers say so
in their headers and in their output, and the per-cell resampling is turned down to the minimum
(`n_perm=200`, `n_random_sets=10`) precisely because that contrast is not what is being computed - the
differential builds its own null.

`prism_differential.py` reports, per set x variant x cell pair:

- **`r_diff`** - the differential correlation, the number that matters;
- **`perm_p`** - drug identity of the predictor differences permuted, which breaks the drug-to-drug
  correspondence while preserving both marginals;
- **`bh_q`** - BH over the whole family, which is large: 10 unordered pairs x ~160 set x variant x statistic;
- **`r_selective`** - the same restricted to the **355** drugs with between-cell AUC spread > 0.1, pre-specified
  because a drug with no between-cell spread carries no signal for this design either way;
- **`naive_r_in_x`** - the undifferenced correlation, printed only as a foil: it is the confounded quantity and
  is expected to be larger and meaningless.

It also prints the fraction of combinations carrying the predicted sign, where 50% is no signal.

### 7.4 The prediction, recorded before the test is run

Report 09, 9.4 criticises itself for choosing a pathway post hoc over 160 combinations. With ten cell pairs on
top of that the risk is worse here, so the direction is fixed now and is in both scripts' docstrings:

> Lower predicted production capacity for biomass precursors - nucleotides, amino acids - in cell X than in
> cell Y should mean **more** growth inhibition in X, hence **lower** AUC in X. Both differences are negative
> together, so **`r_diff` is expected to be POSITIVE**, and most strongly for the biomass-precursor sets rather
> than uniformly across all of them.

A negative `r_diff` contradicts the prediction; it is not the same finding with the sign flipped. And a
positive result spread evenly over every set, negative controls included, is the signature of a confounder
that survived differencing rather than of metabolic specificity.

### 7.5 Coverage, restated

821 drugs have an AUC in all five covered lines and a LINCS `pert_id`; **355** of them have between-cell AUC
spread > 0.1. The groups in the table are a cell-independent median split on mean AUC and exist **only** to
satisfy `mc_phenotype_assoc.py`'s case/control interface - they are not a hypothesis, and keeping them
cell-independent is what guarantees every drug is scored in every cell line, which the differential requires.

### 7.6 The chain is validated, and a bug the validation caught

The per-cell launcher was run end to end on A549 with the local corrected labels
(`results/labels/A549/metchange/rank`, mtime 2026-11-12): 683 of the 821 PRISM drugs have A549 labels,
`per_drug_scores.tsv` is **2.1 MB**, confirming that bringing results home is cheap.

`prism_differential.py` was then exercised on two pseudo-cells in a throwaway directory, so nothing fabricated
was written into `results/`:

- **A null case** - A549's real scores against the same scores with Gaussian jitter added, which gives a
  non-degenerate difference vector carrying no information about the real AUC differences. The script returned
  what it should: **0 of 40 combinations at BH q < 0.05**, all q = 0.924, and **55% carrying the predicted sign
  against the 50% that means no signal**.
- **A degenerate case** - identical scores in both pseudo-cells, so the difference vector is exactly zero.

**The degenerate case exposed a real bug.** With `dS` constant, the correlation is `nan`, and
`|null| >= |r|` is then False for every permutation, so the permutation p-value came out as
`1/(1 + n_perm) = 0.0196` - **a spuriously significant value for an input containing no information at all.**
This is not hypothetical: a set whose members are all absent from both label matrices produces exactly that
vector, and `mc_phenotype_assoc.py` already reports that 23 of 41 sets have members missing from the A549
matrix.

Fixed: combinations whose difference vector is constant or all-missing are now **skipped and listed by name**,
non-finite permutation values are dropped before the p-value is formed, and if nothing is scorable the script
exits with an explanation rather than a table. Re-tested: the degenerate case now skips all 40 combinations
and says why, and the null case is unchanged.

Worth recording as a near-miss of the same family as the others in these reports - the HEPG2-labels-for-A549
fallback, the whole-LINCS-versus-panel overlap, the stale label matrix. **Each would have produced a
confident number rather than an error**, and each was caught by testing a case whose answer was known in
advance rather than by reading the code.

### 7.7 The micromamba environment, and a silent-success bug it exposed

The first cluster run failed with **"No module named numpy"**: the dependencies live in the micromamba
environment `metabolic_drug_repurposing`, and **a Slurm job does not inherit an activated environment**, so
`python3` was the system interpreter. The laptop's system python3 does have numpy and scipy, which is why
nothing showed up in local testing.

Fixed in one place rather than per launcher. `scripts/mc_phenotype_assoc/python_env.sh` (new) resolves an
interpreter and is sourced by `mc_phenotype_assoc_common.sh`, `submit_assoc.sh` and `run_naia_assoc.sh` -
the three places that each had their own `python3` default, which is how this could bite only one of them.

**The root prefix is the part that matters.** The environment lives under a shared project directory, not the
default `~/micromamba`, so the working procedure here is

```bash
export MAMBA_ROOT_PREFIX=/work/projects/car_t_methylome/software/micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate metabolic_drug_repurposing
```

and without that export micromamba looks in the wrong place and reports no such environment. The resolver now
reproduces exactly that, with the root recorded once as `mamba_root_prefix` in `project_paths.conf`
(repeatable, and only directories that exist are tried, so a path for another machine is harmless - the same
convention `scratch_proj_dir` already uses).

Resolution order, first that works winning: an exported `$py`; the current `python3`; **micromamba asked about
itself** (`micromamba run -n`, then `info --base` to supply a candidate root); then each candidate root with
`MAMBA_ROOT_PREFIX` exported first, trying hook-then-activate before falling back to
`<root>/envs/<env>/bin/python3` called directly, which needs neither the hook nor micromamba on `PATH`;
finally an explicit `mamba_env_prefix`. `provenance` now prints the interpreter and how it was found, so a log
says which Python produced the numbers.

Asking micromamba about itself was added after the first version **failed its own test**: with micromamba on
`PATH` and configured but its root absent from the candidate list, nothing was tried at all. Verified against
a stub micromamba in three configurations - on `PATH` with an unlisted root, `MAMBA_ROOT_PREFIX` exported with
no micromamba on `PATH`, and neither available - which resolve, resolve, and refuse respectively.

`require_python` is the loud check, and the launchers call it **before** doing work: it prints what it tried
and the three ways to fix it, then refuses. Nothing is attempted that is going to fail partway through.

**The fix exposed a worse bug than the original.** With no usable interpreter,
`launch_assoc_slurm_generic.sh` **exited 0** - a "successful" Slurm array task with no output. The arm
scripts run several associations in sequence and never aggregate their exit status, so every `run_assoc` was
failing in turn and the script still finished cleanly. A crash is recoverable; a green array task that
produced nothing is the kind of thing that gets discovered weeks later. `require_python || exit 1` now sits at
the top of `launch_assoc_metchange.sh` and `launch_assoc_expr.sh`, and both tasks exit 1 as they should.

Verified: with an unsatisfiable requirement every entry point exits 1 -
`launch_prism_cell.sh`, `prepare_assoc_inputs.sh`, `submit_assoc.sh`, `run_naia_assoc.sh` and both array
tasks of `launch_assoc_slurm_generic.sh` - and with the real dependencies present the resolution is unchanged
on the laptop (`python3 on PATH`).

This is the fourth member of the family catalogued in 7.6: **the failure mode that matters is not the error,
it is the confident-looking success.**

## 8. PRISM differential, run: the prediction FAILS, and the design is weaker than section 1.2 claimed

*All five covered cell lines scored on the cluster against their own MetChange labels
(`run_info.json` confirms `/scratch/.../labels/<CELL>/metchange/rank` for each; 679-760 drugs per cell),
then `prism_differential.py` over 40 set x variant combinations x 10 cell pairs = 400 tests.
Outputs: `results/mc_phenotype_assoc/cross_cell/prism_differential{,_standardised}.tsv`.*

### 8.1 The result

**The pre-registered prediction is not supported.** Of 400 combinations, **171 (43%) carry the predicted
positive sign** against the 50% that means no signal, and of the **13 that pass BH q < 0.05, none has the
predicted sign** - every one is negative:

| pair | set | variant | n | r_diff | perm p | BH q | r_selective |
|---|---|---|---|---|---|---|---|
| A375/MCF7 | glycolysis | cyto_nohub | 684 | **-0.187** | 0.0005 | 3.3e-04 | -0.244 |
| A375/MCF7 | glycolysis | cyto | 684 | -0.156 | 0.0005 | 0.008 | -0.221 |
| A375/PC3 | energy_currency | all | 735 | -0.146 | 0.0005 | 0.009 | -0.218 |
| A375/PC3 | adenylates | all | 735 | -0.141 | 0.0005 | 0.009 | -0.228 |
| A375/PC3 | oxphos | mito | 735 | -0.133 | 0.0005 | 0.015 | -0.180 |

Per the pre-registration in 7.4, a negative `r_diff` **contradicts** the prediction rather than supporting it
with a flipped sign. Effect sizes are small throughout: |r| <= 0.19.

### 8.2 It is one cell line, and standardising does not remove it

**All 13 significant results are pairs containing A375**, and the mean `r_diff` per pair sorts almost
perfectly by whether A375 is in it:

| cell pair | mean r_diff | fraction positive |
|---|---|---|
| A375/PC3 | -0.061 | 18% |
| A375/MCF7 | -0.054 | 5% |
| A375/A549 | -0.043 | 8% |
| A375/HT29 | -0.026 | 35% |
| MCF7/PC3 | -0.009 | 30% |
| A549/MCF7 | -0.006 | 40% |
| HT29/PC3 | +0.008 | 65% |
| HT29/MCF7 | +0.017 | 70% |
| A549/HT29 | +0.021 | 80% |
| A549/PC3 | **+0.035** | 78% |

The first suspicion was arithmetic: differencing only cancels the confounder if neither cell dominates the
difference, and A375 has the largest score dispersion of the five (median within-set SD 0.328 against
0.263-0.303), so `dS` would lean towards it and `r_diff` would inherit A375's own within-cell correlation.
Consistent with that, `corr(r_diff, naive_r_in_x) = +0.40` overall and **+0.45** restricted to A375.

So `--standardise-within-cell` was added, z-scoring each cell's scores across drugs per set before
differencing. **It barely changes anything**: 9 significant instead of 13, still all negative, still all
A375 pairs, 42% positive sign. Variance dominance is not the explanation.

### 8.3 What is actually going on, and the correction to section 1.2

The within-cell correlation `Spearman(set score, AUC)` - the quantity the design was meant to avoid - differs
systematically between cell lines:

| | A375 | A549 | HT29 | MCF7 | PC3 |
|---|---|---|---|---|---|
| glycolysis cyto_nohub | **-0.155** | -0.003 | +0.090 | +0.059 | +0.006 |
| tca_endpoints mito_nohub | +0.160 | -0.018 | +0.230 | +0.313 | +0.123 |
| mean over the 40 combinations | **+0.071** | +0.036 | +0.126 | +0.126 | +0.120 |
| combinations at p < 0.01 | 22 | 11 | **27** | **26** | **28** |

A375's coupling is the weakest and most often negative (32% of sets), where HT29, MCF7 and PC3 sit near +0.12
with 26-28 of 40 sets significant. The A375/MCF7 glycolysis result follows directly: -0.155 in A375 against
+0.059 in MCF7, and `r_diff = -0.180`.

**This is the correction.** Section 1.2 claimed the differential "removes the drug's general potency by
construction, so the confounder is cancelled rather than matched away". That is **half true, and the half it
gets wrong matters**:

- differencing does cancel what is **shared** between the two cells - a compound that is simply toxic
  everywhere has `dP ~ 0` and contributes nothing, exactly as argued;
- it does **not** cancel a **cell-specific** coupling between score and sensitivity. When cell X couples at
  -0.155 and cell Y at +0.059, `corr(dS, dP)` measures *the difference in their confounding*, not the biology.

So `r_diff` is not a clean estimate. The design removes drug-level confounding and leaves cell-level
confounding untouched, and with five cell lines whose couplings range from +0.036 to +0.126 that residue is
the largest signal in the table.

### 8.4 The one thing that does point the predicted way, and why it cannot be claimed

The **within-cell** correlations are mostly **positive**, which *is* the predicted direction: a drug that
lowers predicted capacity (low score) should inhibit growth more (low AUC), giving a positive correlation.
Four of five lines show it, at mean +0.036 to +0.126.

It cannot be claimed as support, for the reason the differential existed in the first place: **26 to 28 of 40
set x variant combinations are significant at p < 0.01 in HT29, MCF7 and PC3.** A hypothesis about
biomass precursors predicts specific sets, not nearly all of them. Near-uniform significance across
glycolysis, TCA, OXPHOS, PPP, adenylates and energy currency alike is the signature of general bioactivity -
cytotoxic drugs lower AUC and move every metabolite set - which is precisely what reports 07-09 spent their
effort excluding and what 1.1 predicted would happen here. The negative-control sets were not included in
this run, and should be next time: if `negctrl_o_glycan` also correlates at +0.12, the matter is settled.

### 8.5 FIRST: the stated hypothesis was not the one tested

`data/validation/prism/notes.md` states the hypothesis as **"Biomass precursors (nucleotides, amino acids) ->
proliferation."** The 11 sets in `metabolite_sets.tsv` are `adenylates, energy_currency, etc_endpoints,
glycolysis, glycolysis_endpoints, oxphos, ppp, pyruvate, redox_couples, tca, tca_endpoints` - **all central
carbon and energy metabolism. There is no nucleotide set and no amino-acid set**; `grep -i
'nucleotide|amino|purine|pyrimidine'` returns nothing. `ppp` feeds nucleotide synthesis and `adenylates` is
AMP/ADP/ATP, and that is the whole of the overlap.

**So section 8 tests an energy-metabolism proxy, not the biomass-precursor hypothesis.** Everything in 8.1-8.4
stands as stated - the design limitation of 8.3 in particular is independent of which sets are used - but the
verdict below must be read as applying to energy metabolism only.

There is a second mismatch. The `z` statistic is standardised **within each perturbation**, so it measures the
set's displacement relative to *that drug's own average feature*, i.e. pathway selectivity within the response.
A drug that uniformly collapses everything scores ~0 on every set. For "collapse of biomass precursors stops
proliferation" the relevant quantity is plausibly *total* disruption, which `z` removes by construction.
`--mode value` also produces `auroc`, and an unstandardised summary would need adding.

### 8.6 Verdict and what would be needed

**PRISM does not support an energy-metabolism version of the hypothesis, and does not test the
biomass-precursor version at all.** The pre-registered differential is null-to-contradicted (0 of 400
significant combinations with the predicted sign), the significant negatives reduce to one cell line's
anomalous within-cell coupling, and the within-cell signal that does run the right way is non-specific enough
to be bioactivity.

This was the best-powered cell-matched test available - 821 drugs, five lines - so the negative carries more
weight than the small-n nulls of report 09, 8.1 and report 10, 5. It does not bear on the HepG2 MMP and
SH-SY5Y ATP results of report 09, which used different phenotypes and, crucially, phenotypes that are not
themselves the confounder.

What a usable version would need:

0. ~~**Nucleotide and amino-acid sets**, which do not exist yet~~ - **DONE 2026-11-24, see 9.5.** This was
   described as needing curation and a scope decision. It needed neither: `SUBSYSTEMS` in
   `build_metabolite_sets.py` is a map from set name to Human-GEM subsystem names, and **Human-GEM has 19
   nucleotide / amino-acid subsystems out of 147**, so the sets are a lookup.
1. **Negative-control sets in the run**, to quantify how much of the within-cell correlation is non-specific.
   One extra `--sets` argument; should have been in the first run.
2. **Cell-level confounding handled explicitly** - for instance regressing each cell's score on its own AUC
   and differencing the residuals, or matching drugs on overall potency within each cell before differencing.
   Differencing alone is not enough, which 8.3 now establishes.
3. **A phenotype that is not viability.** This remains the real obstacle: reports 07-09 established that
   viability and responsiveness are the confounders, and PRISM's phenotype *is* viability. The problem
   identified in 1.1 was not solved by the differential, only moved.

## 9. The comparison is not fair: growth inhibition has many mechanisms

*Raised as an objection to section 8 and it is correct, quantitatively decisive, and changes the conclusion
from "PRISM says no" to "PRISM cannot say". Added `--moa-regex` to `prism_differential.py`; results in
`results/mc_phenotype_assoc/cross_cell/prism_differential_moa.tsv`.*

### 9.1 Only 5% of the library can test the hypothesis

A hypothesis about metabolite production capacity applies to drugs that act **through** metabolism. Growth
inhibition does not require that: a tubulin poison, a topoisomerase inhibitor or a kinase inhibitor lowers AUC
with no primary metabolic lesion, and for those drugs there is no reason the metabolic difference should track
the viability difference.

The PRISM repurposing set is overwhelmingly of that kind. Its commonest MoAs are EGFR inhibitor (22), tubulin
polymerization inhibitor (17), HDAC inhibitor (17), MEK inhibitor (14), topoisomerase inhibitor (13), mTOR
inhibitor (13), CDK inhibitor (13), glucocorticoid receptor agonist (13). **Drugs with a metabolic MoA: 45 of
821, or 5%** - ribonucleotide reductase (4), thymidylate synthase (4), DHFR (3), dihydroorotate dehydrogenase
(3), IMPDH (2), purine antagonists (2), and singletons.

### 9.2 What pooling across mechanisms does to the test

Simulated: 45 drugs with a genuine correlation, embedded in 776 with none, Spearman over all 821.

| true r in the 45 drugs | pooled r reported | P(pooled test reaches p < 0.05) |
|---|---|---|
| 0.0 | -0.000 | 6% |
| 0.2 | +0.011 | 6% |
| 0.4 | +0.022 | 8% |
| 0.6 | +0.030 | 14% |
| **0.8** | **+0.045** | **24%** |

**The actual run reported median |r_diff| = 0.030.** That is indistinguishable from a true r of 0.6-0.8 among
the applicable drugs. So the pooled test of section 8 **cannot separate "no effect" from "a very large
mechanism-specific effect"**, and its null is close to uninformative for the stated hypothesis. This is a
worse problem than the set mismatch of 8.5 and the cell-level confounding of 8.3, because it is about whether
the question was askable at all with this design.

### 9.3 The stratified test: also underpowered, but the aggregate comparison is not

`--moa-regex` now reports each combination for the matching drugs and for the complement, the complement being
the built-in control: **the hypothesis predicts an effect in the matching group and not in the complement.**

Per cell pair only **30** drugs match (fewer than 45, since a drug must be scored in both cells). Individually
that is hopeless after correction: BH over 400 tests needs p < 1.25e-04, which at n = 30 requires |r| ~ 0.65.
Unsurprisingly **0 of 400 reach BH q < 0.05 in the matching group** - the smallest q is 0.995 - and the largest
values are incoherent in sign across pairs for the same set (`HT29/PC3 etc_endpoints` +0.415 against
`HT29/PC3 energy_currency` -0.362), which is what n = 30 noise looks like.

The informative comparison is the aggregate, which pools over all 400 and is therefore far better powered:

| | matching drugs (hypothesis applies) | complement (control) |
|---|---|---|
| predicted positive sign | 184 of 400 (**46%**) | 168 of 400 (**42%**) |
| median r | **-0.015** | **-0.010** |

**The hypothesis predicts the matching group to exceed the complement. It does not** - 46% against 42% is
nothing, and the medians are both slightly negative and essentially equal. Caveat: the 400 tests are not
independent (10 pairs sharing 5 cells, nested set variants), so the effective n behind that comparison is much
smaller than 400 and it should not be read as a precise estimate - only as the absence of the predicted
separation.

Also worth stating: **the amino-acid arm of the hypothesis is untested for lack of drugs.** Almost all 45
matching compounds are nucleotide- or one-carbon-directed; the library has essentially no amino-acid-metabolism
agents.

### 9.4 Revised verdict on PRISM

Sections 8.5, 8.3 and 9.2 compound rather than overlap, and together they change the conclusion:

1. the sets tested were central carbon and energy metabolism, **not** the nucleotides and amino acids the
   hypothesis names (8.5);
2. differencing cancels drug-level but **not cell-level** confounding, and the residue was the largest signal
   in the table (8.3);
3. **95% of the library cannot test the hypothesis**, and pooling them makes a true r of 0.8 look like
   r = 0.045 (9.2);
4. the stratified test has no power after correction at 30 drugs per pair, and the one well-powered
   comparison - matching versus complement - is flat (9.3).

**So PRISM does not refute the hypothesis; it is structurally a weak instrument for it.** Report 10, section
1.2 called it "a far better-powered version of the cross-cell-line test" on the strength of 821 drugs and 5
cell lines. That counted compounds without asking how many were *relevant* compounds, which is the same error
as counting whole-LINCS overlap without asking about the cell-line panel (5.1, and the rule recorded there).
**The lesson generalises: usable n is the number of drugs for which the hypothesis makes a prediction, not the
number of drugs in the file.**

What PRISM remains genuinely good for: a well-powered **cell-line-selectivity** readout for any question that
applies to bioactive compounds generally rather than to metabolically-acting ones, and the 8.3 diagnostic of
how much cell-specific score-to-sensitivity coupling each label set carries - which is worth knowing
independently of any hypothesis.

## 9.5 The biomass-precursor sets existed all along; they were missing from our set file, not from Human-GEM

*Correction. Section 8.5 said "there is no nucleotide set and no amino-acid set" - true of
`metabolite_sets.tsv`, which had 11 hand-written central-carbon and energy sets - but section 9.4 then listed
building them as needing curation and a scope decision. **That was wrong.** `SUBSYSTEMS` in
`build_metabolite_sets.py` maps a set name to Human-GEM subsystem names and the builder takes every metabolite
in a reaction of those subsystems. **Human-GEM has 19 nucleotide / amino-acid subsystems of 147**, so this was
a lookup, not a curation problem.*

The subsystems, with their reaction counts:

| | reactions | | reactions |
|---|---|---|---|
| Phenylalanine, tyrosine and tryptophan biosynthesis | 157 | Glycine, serine and threonine metabolism | 61 |
| Nucleotide metabolism | 119 | Valine, leucine, and isoleucine metabolism | 53 |
| Purine metabolism | 108 | Cysteine and methionine metabolism | 44 |
| Arginine and proline metabolism | 99 | Alanine, aspartate and glutamate metabolism | 41 |
| Tyrosine metabolism | 78 | Lysine metabolism | 32 |
| Pyrimidine metabolism | 75 | Histidine metabolism, Tryptophan, Beta-alanine, Phenylalanine, Urea cycle, Metabolism of other amino acids | 6-22 each |
| Folate metabolism | 64 | | |

### Added to `build_metabolite_sets.py`

**Subsystem sets:** `nucleotide` (192 metabolites, `all`), `purine` (150), `pyrimidine` (140), `folate` (97),
and one combined `amino_acid` (823) rather than fourteen separate ones, to stop the BH family doubling.
`Amino sugar and nucleotide sugar metabolism` is deliberately excluded - it is glycan-adjacent and belongs
nearer the `NEG_CONTROL` sets than to biosynthesis of either precursor class.

**Curated endpoint sets**, in the style of the existing `*_endpoints`:

| set | members (`all`) | note |
|---|---|---|
| **`dntp_endpoints`** | 14 | **the tightest test of the nucleotide arm** - dATP/dGTP/dCTP/dTTP are the direct products of ribonucleotide reductase and thymidylate synthase, the best-represented metabolic targets in the PRISM library (RNR 4 drugs, TS 4, DHFR 3), and unlike the NTPs they are **not** currency hubs, so `all_nohub` is also 14 |
| `ntp_endpoints` | 20 | ATP/GTP/CTP/UTP; no `nohub` variant exists because all four are hubs, which is itself worth knowing |
| `nmp_endpoints` | 39 | AMP/GMP/CMP/UMP/IMP/dTMP/dUMP |
| `amino_acid_endpoints` | 89 | the 20 proteinogenic amino acids |

### Two naming traps, both hit while doing this

1. **The amino acids carry no `L-` prefix in Human-GEM** - `alanine`, not `L-alanine`. Probing with the
   prefixed form returned 1 of 20 and looked like the model was missing them.
2. **Do not check metabolite names in `data/networks/combined/nodes.tsv`.** That is the *built* network, and
   currency metabolites and hubs are filtered out during construction, so ATP, GTP, CTP, UTP and every amino
   acid are absent from it **by design**. Checking there suggested the GEM lacked ATP. The authoritative
   source is `data/GEMs/Human-GEM/model/Human-GEM.yml`, parsed the way the builder parses it; note
   `model/metabolites.tsv` carries only cross-reference IDs and no `metName` column.

The builder warns on an unmatched subsystem and on an unmatched curated metabolite, so a typo is loud rather
than silent - but only if the name is checked against the right file.

### Next step

The sets are built locally. To use them the label matrices must be rescored, which means the cluster: copy the
regenerated `metabolite_sets.tsv` across (it is small) or rebuild it there with
`rebuild_sets=true cell=ALL bash scripts/mc_phenotype_assoc/prepare_assoc_inputs.sh`, then re-run
`launch_prism_slurm.sh`. The PRISM phenotype table does not change.

Note for interpretation when that lands: adding 9 sets enlarges the BH family, so the correction gets
stricter, and `amino_acid` at 823 metabolites is a seventh of the whole matrix - a set that large is close to
a global measure and should be read with that in mind. `dntp_endpoints` is the set to look at first, because it
is small, mechanistically specific, and matched to the drugs that are actually present.

## 10. With the biomass sets and target-based drug selection, the NUCLEOTIDE arm behaves as predicted and the AMINO-ACID arm contradicts

*The cluster rescore with the 9 new sets landed the same day (`run_info.json`: 82 sets scored per cell,
679-760 drugs). Drug selection by annotated target rather than MoA free text was suggested before any of these
results existed. Outputs: `prism_differential_targets.tsv` (energy + biomass targets) and
`prism_differential_nucleotide_targets.tsv` (nucleotide-pathway targets only).*

### 10.1 The effect appears only when the drugs are mechanism-matched, and grows as the match tightens

Mean `r_diff` over cell pairs, with the pre-registered prediction being **positive**:

| set | all 821 drugs | energy+biomass targets (~80/pair) | **nucleotide-pathway targets only (~41/pair)** | complement |
|---|---|---|---|---|
| **`dntp_endpoints`** | -0.036 | +0.099 (60% of pairs +) | **+0.136 (80% +)** | -0.046 |
| `ntp_endpoints` | -0.040 | +0.110 (60% +) | **+0.154 (80% +)** | -0.051 |
| `nmp_endpoints` | -0.039 | +0.083 | +0.098 (68% +) | -0.047 |
| `nucleotide` | -0.044 | +0.086 | +0.105 (65% +) | -0.053 |
| `pyrimidine` | -0.041 | +0.080 | +0.091 (68% +) | -0.050 |
| `purine` | -0.051 | +0.055 | +0.081 (55% +) | -0.060 |
| `folate` | -0.012 | +0.016 | +0.018 (53% +) | -0.016 |
| **`amino_acid`** | -0.018 | **-0.078 (27% +)** | **-0.100 (20% +)** | -0.012 |
| `amino_acid_endpoints` | +0.001 | -0.056 (30% +) | -0.066 (30% +) | +0.008 |

Selection sizes: 90 drugs match the energy+biomass subsystems with `effect=inh`, 44 match the
nucleotide-pathway subsystems only.

Three features of that table matter:

1. **The sign flips with mechanism matching.** Pooled over all drugs every nucleotide set is slightly
   negative; restricted to drugs that inhibit a nucleotide-pathway enzyme they are positive. That is the
   direction predicted in 7.4 before any of this was run.
2. **The effect grows as the match tightens** (-0.036 -> +0.099 -> +0.136 for `dntp_endpoints`) while n falls
   from 684 to 41, which is the opposite of what a chance finding does.
3. **The complement behaves as a control should** - about -0.05 throughout, i.e. no effect in the drugs where
   the hypothesis does not apply.

### 10.2 The aggregate statistics, and which ones are legitimate

**Only `dntp_endpoints` and `ntp_endpoints` have a single variant**, so their 10 rows are genuinely one per
cell pair. The other sets have 6 compartment/hub variants each, so their "60 pairs" are nested and their sign
counts are inflated - those rows must not be read as 60 independent tests.

For the two clean sets, nucleotide-targeting drugs, over the 10 cell pairs:

| | `dntp_endpoints` | `ntp_endpoints` |
|---|---|---|
| pairs with the predicted positive sign | **8 of 10** | **8 of 10** |
| sign test (one-sided) | p = 0.055 | p = 0.055 |
| mean Fisher z | +0.140 | +0.159 |
| p from the across-pair spread | 0.005 | 0.001 |
| **matched vs complement, paired Wilcoxon over the same 10 pairs** | **p = 0.0020** | **p = 0.0010** |
| strongest single pair | HT29/MCF7 +0.374 (p = 0.018) | HT29/MCF7 +0.396 (p = 0.012) |

The **paired Wilcoxon is the right statistic**: it tests the crossed prediction directly - matched drugs
exceeding the complement - and being paired on cell pair it removes pair-level effects, which is what defeated
the pooled analysis of section 8. **It is also one test per set rather than per pair**, so its multiplicity
burden is 9 sets rather than 820 combinations; see 10.3 for the corrected q-values and for why the same
statistic must be ignored on the six-variant sets.

The amino-acid arm fails the same test decisively in the opposite direction: 2 of 10 pairs positive, mean
Fisher z -0.102, paired Wilcoxon p = 1.000 against the prediction.

### 10.3 What this is and is not

**Is:** a mechanism-matched, directionally pre-registered association in the nucleotide arm, with a
correctly-behaving control group, which strengthens monotonically as the mechanistic match tightens. On a
phenotype (viability) that is itself the confounder of reports 07-09, and in a design where the pooled version
is provably unable to see it (9.2).

**Correction, and what is and is not corrected.** Every per-pair `p_moa` quoted above is **uncorrected**, and
an earlier draft described "p = 0.018 against a BH threshold of ~6e-05", which conflated two different things.
The 6.1e-05 is α/m = 0.05/820, i.e. the BH threshold a p-value must beat *at rank 1*; it is not a q-value.
The actual corrected numbers are:

| family | smallest uncorrected p | **actual BH q** | rank-1 threshold |
|---|---|---|---|
| per-pair tests, 820 = sets x variants x pairs | 0.0019 | **0.62** | 6.1e-05 |
| pre-registered `dntp_endpoints` alone, 10 tests | 0.0176 | **0.176** | 0.005 |

So **no individual cell pair survives**, and not narrowly: q = 0.62 over the full family, and still 0.176 over
the pre-registered single-set family of 10. For p = 0.018 to pass BH over 820 it would need rank >= 296, and
only 12 tests sit at p <= 0.018.

**The aggregate statistic belongs to a much smaller family, and it does survive.** The paired Wilcoxon is one
test *per set*, not per pair, so m = 9:

| set | rows | p | **BH q over the 9 sets** |
|---|---|---|---|
| `ntp_endpoints` | 10 | 9.8e-04 | **0.0018** |
| `dntp_endpoints` | 10 | 0.0020 | **0.0029** |
| `folate` | 60 | 0.033 | 0.043 |
| `amino_acid_endpoints` | 40 | 0.9985 | 1.000 |
| `amino_acid` | 60 | 1.000 | 1.000 |
| `pyrimidine`, `nucleotide`, `nmp_endpoints`, `purine` | 60 | 7e-10 to 1e-07 | **not usable - see below** |

**The four sets with q near 1e-09 must be discarded, not celebrated.** Their Wilcoxon runs over 60 rows that
are 10 cell pairs x 6 nested compartment/hub variants of the same set, so the test treats pseudo-replicates as
independent and the p-value is inflated by roughly the variant multiplicity. That is precisely why 10.2
restricts the aggregate claim to `dntp_endpoints` and `ntp_endpoints`, the only sets with one row per cell
pair.

**And the surviving aggregate still rests on 10 cell pairs that share only 5 cell lines**, so its effective n
is below 10; halving the degrees of freedom puts it near p = 0.03, q near 0.05. It is a real pre-specified
result at the edge of significance, not a comfortable one.

**Pre-registration, honestly apportioned - and the phrase "pre-specified aggregate" withdrawn.**

Fixed before any result existed, with a written record: the **positive direction** (7.4, in
`prism_differential.py`'s docstring before the differential was first run), **`dntp_endpoints` as the set to
read first** (9.5 and STATUS, written when the sets were built and before the rescore landed), and the
**complement-as-control** design (in the `--moa-regex` help text before any stratified run).

**Not** fixed in advance: the narrowing of the drug set from 90 to 44, chosen after seeing the 90-drug result;
and - the one that matters most - **the aggregate statistic itself**, chosen after it became clear the
per-pair tests would not survive correction. Calling the result a "pre-specified aggregate" was therefore
wrong. It is an aggregate test on a pre-specified set, in a pre-specified direction, **with a statistic
selected after the fact.**

That is not a small distinction, because three defensible aggregates on the same data span a factor of 25:

| statistic | what it tests | p |
|---|---|---|
| sign test, 8 of 10 pairs positive | how many pairs point the predicted way; ignores magnitude and ignores the control | **0.055** |
| mean Fisher z with the across-pair SE | whether the average r differs from zero; ignores the control | 0.005 |
| paired Wilcoxon, matched vs complement | whether matched exceeds complement - the crossed prediction | **0.0020** |

**10.2 led with 0.0020, the most favourable of the three.** The Wilcoxon remains the best choice *by design*,
being the only one that uses the control group the crossed prediction depends on - but it is also the one
least likely to have been chosen had it come out worse, and the most assumption-free statistic, the sign test,
gives p = 0.055. **Read 0.055 as the conservative figure and 0.0020 as the design-appropriate one, and treat
the gap between them as the cost of having picked the statistic afterwards.**

The sequence all drugs -> MoA -> any GEM target -> nucleotide targets is four looks at one dataset, so the
monotone trend of 10.1 is reassuring but is not a substitute for a fresh test (10.4, item 1).

**And half the hypothesis is contradicted.** `notes.md` names nucleotides *and* amino acids. The amino-acid
sets run significantly the wrong way at every level of matching, and get worse as the match tightens. Either
amino-acid capacity is not limiting for proliferation in these lines, or the amino-acid sets - 823 metabolites,
a seventh of the matrix - are too broad to mean anything. `amino_acid_endpoints` (89, the 20 proteinogenic
acids) is also negative, which argues against the second explanation.

### 10.4 What would settle it

1. **A held-out test.** The forking path of 10.3 is the main threat. The cleanest resolution is the other
   cell lines: HEPG2, VCAP and HCC515 are absent from PRISM, but the PRISM primary screen covers far more
   compounds at a single dose, and the CTRP / GDSC viability sets cover other lines - a pre-specified
   `dntp_endpoints` test there, with the direction and the drug-selection rule fixed in advance, needs no new
   decisions.
2. **Drug-specific target-proximal sets**, the refinement that was identified but not yet built: score each
   drug against the metabolites within one or two reactions of *its own* target rather than against a shared
   nucleotide set. That collapses 820 tests to 10 and is the single biggest available power gain.
3. **Mechanism grouping.** `--group-jaccard` is implemented and unused here: 237 of 747 drugs share a target
   vector with another drug, so the per-pair n of 41 contains pseudo-replicates. Grouping will lower n and
   should be run before any number from 10.2 is quoted as final.
4. **A375.** Its pairs are the two that go negative for `dntp_endpoints`, consistent with the cell-level
   coupling anomaly of 8.3. Worth checking whether excluding A375 strengthens the result - and worth being
   suspicious of a result that improves by dropping a cell line.

## 10.5 The PRISM runs have no expression arm, which is the gap that matters most

The project's thesis is not that MetChange labels associate with metabolic phenotypes - it is that they do so
**in a way the transcriptome does not**, which is what justifies the stoichiometric layer before any
investment in multi-cell-line FNN training. Two instances support that: report 09, section 4 (membrane
potential: MetChange 0.365 matched, expression 0.500) and report 09, 9.6 (ATP: MetChange 0.173 matched,
expression 0.329, not significant).

**Both of those use HEPG2 labels.** There are no SH-SY5Y labels and none are possible (report 09, 8.1), so
9.6's SH-SY5Y is the phenotype's cell line only. The two instances therefore vary phenotype, phenotype cell
line, compound library and the nature of the confounder control - but **not the label cell line**. Anything
peculiar to HepG2's label generation is common to both.

**The PRISM arm could close that gap and currently does not.** `results/mc_phenotype_assoc/<cell>/` shows an
`_expr` counterpart for every Attene and Hallinger run; `prism_viability` has none. So section 10's nucleotide
result, whatever its eventual strength, says nothing about whether expression would track the same thing
equally well - and it is the one analysis in the project that would test the thesis on **five non-HepG2 cell
lines at once**.

The work needed is small: `launch_prism_cell.sh` runs only the MetChange arm, and an expression arm follows
the existing `launch_assoc_expr.sh` pattern. `gene_sets.tsv` already carries the nucleotide and amino-acid
**gene** sets from the 9.5 rebuild (`--emit genes`), so nothing else is missing. The prediction to fix in
advance, before running it: **if the nucleotide association is metabolic rather than transcriptional, the
expression arm should be weaker than the MetChange arm on the same drugs, sets and cell pairs** - and if it is
equally strong, the nucleotide result ceases to support the thesis regardless of its own p-value.

### 10.6 Implemented: the expression arm, and the conditional statistic it enables

**Framing first, because it changes the statistic.** The claim is not that MetChange beats the transcriptome
in a horse race - expression being informative too is expected and unremarkable, since MetChange is computed
*from* expression. The claim is that **integrating stoichiometry and the metabolic network yields information
expression alone does not**. That is a question about *increment*, and the right statistic is conditional, not
comparative.

> **Which means the increment question has not yet been asked anywhere in this project.** Every result so far
> used an unconditional statistic: sections 8, 9 and 10 are plain Spearman of the differences, and the two
> existing MetChange-versus-expression comparisons (report 09 sections 4 and 9.6) are **side-by-side** - greedy
> covariate matching followed by Mann-Whitney, run separately on each label type. MetChange surviving a
> matching control that expression fails is suggestive of increment but does not measure it: two statistics can
> both be non-null while sharing all of their information, and one can fail a control for reasons unrelated to
> the other. `r_partial_mc_given_expr` would be the first direct measurement. That makes it more valuable than
> a confirmation of existing results - and it means the present evidence for the central thesis is weaker *in
> kind*, not merely in quantity, than "MetChange 0.365 against expression 0.500" suggests.

`launch_prism_cell.sh` now runs **two arms**: MetChange labels with the metabolite sets, and the expression
labels with the GPR gene sets, on the same drugs and cell pairs. `arm=metchange` or `arm=expr` runs one alone;
the default is both. `launch_prism_slurm.sh` passes it through. The expression arm skips itself with a message
if the cell has no `expr` label directory or if `gene_sets.tsv` is missing, rather than failing the task.

`prism_differential.py --expr-run prism_viability_expr` then reports, per cell pair and set:

| column | meaning |
|---|---|
| `r_diff` | the MetChange differential, as before |
| `r_expr` | the same differential computed on the matching **gene** set |
| **`r_partial_mc_given_expr`** | **MetChange differential with the gene set partialled out** - the statistic that matters |
| `p_partial` | its p-value, partial Spearman on ranks, df = n-3 |

Only the 10 subsystem sets have same-named gene counterparts; the 10 curated endpoint sets do not, so
`--expr-set-map` pairs each to its parent subsystem. The defaults are `dntp_endpoints`, `ntp_endpoints` and
`nmp_endpoints` -> `nucleotide`; `amino_acid_endpoints` -> `amino_acid`; `tca_endpoints` -> `tca`;
`glycolysis_endpoints` -> `glycolysis`; `etc_endpoints`, `energy_currency`, `adenylates`, `redox_couples` ->
`oxphos`. The summary line reports what fraction of the MetChange magnitude survives conditioning, where 100%
means expression explains none of it and 0% means all of it.

Verified on a synthetic expression arm in a throwaway directory - the mechanics run, the set pairing resolves
and the summary computes. **That was a plumbing test on fabricated data and is not a result.** With no
expression arm present the script says so and continues with the unconditional analysis.

### 10.7 Pre-registered prediction for the expression arm, fixed before it is run

Written here before the arm exists, so that §10.3's complaint about a post-hoc statistic does not recur:

1. **Primary statistic:** `r_partial_mc_given_expr` on `dntp_endpoints`, nucleotide-pathway-targeting drugs,
   aggregated across the 10 cell pairs by paired Wilcoxon against the complement - the same specification as
   10.2, with the partial correlation substituted for the raw one.
2. **Prediction:** the MetChange differential **retains most of its magnitude** once the matching gene set is
   partialled out, i.e. the retention fraction stays well above 50% and the paired test remains in the
   predicted positive direction.
3. **What refutes it:** retention near or below 50%, or a partial correlation that loses significance while
   the raw one keeps it. Either would mean the nucleotide association is transcriptional information wearing a
   metabolic label, and §10 would stop supporting the thesis **whatever its own p-value**.
4. **What is NOT claimed either way:** `r_expr` being large. Expression carrying the same signal is
   compatible with the thesis; expression carrying *all* of it is not.
5. **Secondary, reported but not primary:** the same three columns for the subsystem sets `nucleotide`,
   `purine`, `pyrimidine`, `folate`, whose gene counterparts are exact rather than mapped. These are the
   like-for-like comparisons (same subsystem, metabolites against its own GPR genes) and are the cleanest
   test of increment, though they lack `dntp_endpoints`' mechanistic sharpness.

## 11. Jaccard grouping both ways, and why the section-10 result must now be frozen rather than refined further

*All three analyses requested: the aggregation form of Jaccard grouping, the stratified form, and the
target-proximal sets. Two were runnable locally from the existing per-drug scores; the third is not.*

### 11.1 Jaccard as aggregation: it substantially weakens section 10

Collapsing each mechanism group to its median before correlating, so 8 near-identical drugs count once:

| `dntp_endpoints`, nucleotide-targeting drugs | ungrouped | Jaccard >= 0.5 | Jaccard >= 0.7 |
|---|---|---|---|
| n (drugs -> groups) | 40 | 23 | 30 |
| mean r | **+0.136** | +0.055 | +0.056 |
| cell pairs positive | **80%** | 60% | 60% |
| paired Wilcoxon p | **0.0020** | **0.080** | **0.053** |

`nucleotide`, `purine` and `pyrimidine` go from positive to **negative** under grouping. The effect size
halving is not a power artefact - that is the estimate itself moving - so specific groups of near-identical
drugs were carrying it.

**Which ones: the fluoropyrimidines.** Among the 44 matched drugs, Jaccard >= 0.5 gives 26 groups with sizes
8, 4, 4, 3, 3, 2 and 20 singletons. The group of 8 is `{DPYD,TYMS}` - capecitabine, carmofur, doxifluridine,
floxuridine, ftorafur, trifluridine and relatives, all 5-FU prodrugs or close analogues. Collapsing one
chemotype family from 8 observations to 1 is most of the drop.

### 11.2 Jaccard as stratification: this library cannot support it for the metabolic arm

The stratified form - correlate *within* each mechanism group, to see which mechanisms carry an association -
is the better-motivated design, and it is not feasible here. **No matched group reaches 10 drugs**: the
largest is the 8 fluoropyrimidines, then 4 IMPDH inhibitors, 4 polymerase/nucleoside analogues, 3, 3, 2 and
20 singletons. A within-group correlation on 8 drugs is not interpretable.

What it does deliver is a **negative-control panel**, because the groups that *are* large enough are all
non-metabolic:

| group | n drugs | median r over the 10 pairs | pairs positive |
|---|---|---|---|
| MET | 35 | +0.016 | 50% |
| ATR, MTOR | 30 | -0.012 | 50% |
| FGFR1/2/3... | 19 | -0.039 | 50% |
| HDAC1/2/3... | 16 | +0.103 | 70% |
| AR, CYP2C8, ESR1... | 16 | +0.057 | 50% |
| MAP2K1, MAP2K2 | 13 | -0.217 | 20% |
| HSPA1A | 13 | -0.157 | 30% |
| TUBA1A/1B/1C... | 10 | +0.236 | 60% |
| AURKA/B/C... | 10 | +0.194 | 70% |

Median of the group medians: **+0.016** over 9 non-metabolic mechanisms. Individual values scatter from
-0.217 to +0.236 with sign fractions of 20-70%, i.e. noise at this n. That is a genuine, if weak,
reassurance: **the design does not manufacture associations in mechanisms where the hypothesis does not
apply.** It is weak because the per-group estimates are plainly unstable at 10-35 drugs.

### 11.3 A contamination found while doing this, and the result it produces

Listing the 44 matched drugs individually showed that **16 are not antimetabolites at all**. Nine are
**phosphodiesterase inhibitors** - papaverine, anagrelide, enoximone, cilostamide, zardaverine, YM-976,
trequinsin, avanafil, arctigenin - which qualify because cAMP/cGMP hydrolysis is filed under purine
metabolism in Human-GEM. A PDE inhibitor does not deplete nucleotide pools. The others: tofogliflozin (GCK),
adaprev (HK1), PJ-34 (EEF2), spermidine and RITA (TXNRD1/2), curcumin and diphenyleneiodonium (XDH), aspirin
(ADA/PDE). This is the same failure as the PI3K-in-inositol-phosphate problem of 10.1, one level subtler:
**GEM subsystem membership is bookkeeping, not mechanism.**

`--exclude-target-genes '^PDE|^TXNRD|^XDH|^GCK$|^HK1$|^EEF2$'` leaves **28 drugs**, of which the genuine
antimetabolite classes are TYMS (6), DHFR (3), RRM1/2/2B (5), IMPDH1/2 (6), DHODH (3), plus ADA, ADK, CDA,
TYMP, PNP, HPRT1, CMPK1, FPGS singletons:

| `dntp_endpoints`, clean 28 | ungrouped | Jaccard >= 0.7 | Jaccard >= 0.5 |
|---|---|---|---|
| n | 28 | 18 groups | 14 groups |
| mean r | **+0.291** | +0.240 | +0.187 |
| cell pairs positive | **10 of 10** | 90% | 70% |
| paired Wilcoxon p | 0.00098 | **0.0020** | 0.019 |

The two thresholds bracket it: 0.7 collapses only near-identical target vectors and retains 0.240 at 90% of
pairs, 0.5 merges more aggressively and gives 0.187 at 70% (p = 0.019). Both keep the direction and both are
far better behaved than the contaminated 44-drug version, which fell to p = 0.080 at Jaccard 0.5. That the
clean set survives pseudo-replication control where the contaminated one did not is the one piece of internal
evidence here that is not simply a consequence of narrowing.

**Exact values, and a correction to an earlier version of this paragraph.** Read from
`results/mc_phenotype_assoc/cross_cell/prism_clean28_*.tsv`:

| run | n | mean r | pairs with r_moa > 0 | pairs with r_moa > r_other | Wilcoxon p | floor at n=10 |
|---|---|---|---|---|---|---|
| ungrouped | 28 | +0.291 | 10/10 | **10/10** | **0.000977** | 0.000977 |
| Jaccard >= 0.7 | 18 | +0.240 | 9/10 | **9/10** | **0.001953** | 0.000977 |
| Jaccard >= 0.5 | 14 | +0.187 | 7/10 | **8/10** | **0.018555** | 0.000977 |

*The ungrouped p-value is saturated.* **0.000977 is exactly the floor of a one-sided Wilcoxon at n = 10 pairs**
(1/2^10), so it conveys only "all ten cell pairs favoured the matched group" and carries no effect-size
information. The same pattern would read 1e-06 with 20 pairs and could not go below 0.0039 with 8. The honest
form is "10 of 10 pairs, the strongest this test can report at this n".

> **CORRECTION.** An earlier version of this paragraph claimed the Jaccard-0.7 row *also* sat on that floor and
> was therefore "indistinguishable from the ungrouped row by p-value", and used it to illustrate that
> `frac r>0` and the Wilcoxon can disagree. **Both claims were wrong, and both came from a number I typed into
> a demonstration instead of reading from the output.** The Jaccard-0.7 p is **0.001953 - exactly twice the
> floor**, meaning 9 of 10 differences favoured the matched group, not 10. The two rows are distinguishable.
> And in all three rows here the two statistics in fact track each other closely (10/10 and 10/10, 9/10 and
> 9/10, 7/10 and 8/10), so the disagreement I illustrated did not occur. They remain *different* statistics -
> `frac r>0` counts `r_moa > 0` while the Wilcoxon tests `r_moa > r_other`, and with the complement near -0.05
> a slightly negative `r_moa` could still win - but that is a possibility, not something observed here.

`amino_acid` remains negative ungrouped (-0.084, 35% of pairs). Under grouping it moves to +0.032 with a
nominal p of 0.023, which is incoherent with its ungrouped value and is a direct warning that **the grouped
estimates at n = 14 are too noisy to read precisely.**

### 11.4 This sequence can no longer establish anything, and the result must be frozen

The effect has grown at every refinement:

| selection | matched n | mean r on `dntp_endpoints` |
|---|---|---|
| all drugs | 684 | **-0.036** |
| MoA regex | 30 | (sign-incoherent, 0 of 400 at BH q<0.05) |
| any energy+biomass GEM target | 80 | **+0.099** |
| nucleotide-pathway target | 40 | **+0.136** |
| nucleotide target, contaminants removed | 28 | **+0.291** |

**That is the textbook signature of a forking path, and the fact that each individual step was defensible
does not change it.** Every narrowing was chosen after seeing the previous result, and every one moved in the
direction of a larger effect. A PDE inhibitor genuinely is not an antimetabolite - but I removed them *after*
the number existed, not before, and would not have gone looking had the number been satisfactory.

**So none of 10.1-11.3 is evidence for the hypothesis in a confirmatory sense.** What the sequence has
produced is something else and genuinely useful: a **fully specified hypothesis with no free parameters
left.**

```
drug set   : the 28 PRISM compounds with an inhibited target in a Human-GEM
             nucleotide / purine / pyrimidine / folate subsystem, after excluding
             targets matching ^PDE|^TXNRD|^XDH|^GCK$|^HK1$|^EEF2$
metabolite : dntp_endpoints (dATP, dGTP, dCTP, dTTP), variant `all`
statistic  : per cell pair, Spearman of (setscore difference, AUC difference);
             aggregated across pairs by paired Wilcoxon against the complement
grouping   : Jaccard >= 0.5 on target vectors, group medians (the conservative form)
direction  : positive
secondary  : the same on amino_acid / amino_acid_endpoints, predicted NULL or
             negative, as the internal contrast
```

Frozen here. **The next number quoted for this hypothesis must come from data none of the above touched** -
PRISM primary, CTRP once its usable n is checked, or GDSC - run once, with this specification and nothing
added.

### 11.5 Target-proximal sets: not runnable locally, and what it would take

This is the third request and the only one that cannot be done from the per-drug scores. Those files contain
the 87 predefined set x variant combinations; a target-proximal set is **per drug** - the metabolites within
one or two reactions of *that drug's own* target - so it is a different predictor and needs a scoring pass
against the label matrices, i.e. the cluster.

It is also now the most valuable remaining analysis, for a reason 11.1-11.3 makes concrete. With a shared
`dntp_endpoints` the 28 drugs are scored against four metabolites regardless of mechanism, so a DHFR
inhibitor, a DHODH inhibitor and a TYMS inhibitor are all judged on dNTP pools. Target-proximal sets would
score each on its own neighbourhood - folate/one-carbon for DHFR, orotate for DHODH, dTMP/dUMP for TYMS -
which is both mechanistically correct and **collapses the test family from 820 combinations to roughly 10**,
the single largest power gain available.

What it needs: a new scoring mode that takes a drug -> metabolite-set mapping built from the network
(`nodes.tsv` gives reaction GPRs and reaction-metabolite membership, so one or two hops from a target gene is
a short graph walk), and then one score per drug against its own set. The sensible shape is a
`--drug-sets FILE` argument to `mc_phenotype_assoc.py` with columns `drug, node_id`, scored like the existing
sets but selected per column. Until that exists, 11.4's frozen specification stands as the testable claim.

### 11.6 What the test actually is, level by level, and the dependency it rests on

Spelled out because three different n's are in play and they are easy to conflate.

**Where the 28 drugs come from.** 821 PRISM compounds in the phenotype table (an AUC in all five covered
cells and a LINCS `pert_id`) -> 747 with at least one annotated target -> 44 with an **inhibited** target that
is a Human-GEM gene whose reactions fall in a `nucleotide|purine|pyrimidine|folate` subsystem -> **28** after
excluding `^PDE|^TXNRD|^XDH|^GCK$|^HK1$|^EEF2$`.

**Level 1, within a cell pair: Spearman across drugs.** Two correlations per pair, one on the 28 matched drugs
and one on the ~650 complement drugs:

| cell pair | rho matched (28) | rho complement (~650) | difference |
|---|---|---|---|
| A375/A549 | 0.259 | -0.048 | +0.307 |
| A375/HT29 | 0.238 | -0.097 | +0.335 |
| A375/MCF7 | 0.502 | -0.053 | +0.555 |
| A375/PC3 | 0.094 | -0.092 | +0.186 |
| A549/HT29 | 0.349 | +0.030 | +0.319 |
| A549/MCF7 | 0.190 | -0.091 | +0.281 |
| A549/PC3 | 0.288 | -0.070 | +0.359 |
| HT29/MCF7 | 0.385 | +0.008 | +0.377 |
| HT29/PC3 | 0.169 | -0.056 | +0.225 |
| MCF7/PC3 | 0.434 | +0.009 | +0.425 |

**Level 2, across cell pairs: Wilcoxon signed-rank on those 10 differences**, one-sided, `greater`. So the
test's **n is 10 cell pairs** - not 28 drugs and not the 820 set x variant x pair combinations. W = 55, the
maximum possible at n = 10, hence p = 1/2^10.

**The dependency this rests on, which earlier sections did not state.** The 10 pairs are drawn from **5 cell
lines**, so each cell appears in 4 of them; and **the same 28 drugs are used in every pair**. The 10
differences are therefore not independent in either direction, while the Wilcoxon assumes they are - so its
p-value is optimistic on top of being saturated.

Leave-one-cell-out, keeping only the 6 pairs that exclude a given cell:

| dropped | pairs | mean rho matched | mean rho complement | mean difference | p |
|---|---|---|---|---|---|
| A375 | 6 | 0.302 | -0.028 | +0.331 | 0.0156 |
| A549 | 6 | 0.304 | -0.047 | +0.350 | 0.0156 |
| HT29 | 6 | 0.295 | -0.057 | +0.352 | 0.0156 |
| MCF7 | 6 | 0.233 | -0.056 | +0.288 | 0.0156 |
| PC3 | 6 | 0.320 | -0.042 | +0.362 | 0.0156 |

Per-cell mean rho matched: A375 +0.273, A549 +0.271, HT29 +0.285, MCF7 +0.378, PC3 +0.246. **No single cell
drives the result**, which is the useful part of this check.

And note every one of those p-values is *also* at its own floor (1/2^6 = 0.0156, all six differences
positive). **Every p-value in this analysis is saturated**: they report sign consistency and nothing else.
**The defensible summary of the frozen hypothesis is therefore the effect size and its stability - matched
rho about 0.25 to 0.38 per cell line against a complement pinned near -0.05 - and not any p-value quoted from
it.** That is also why the confirmatory cohort needs more *cell lines*, not merely more drugs: with 5 cells
the aggregate test cannot resolve below 0.000977 however strong the effect.

## 12. The expression arm, run: the pre-registered conditional prediction is confirmed

*All five expression arms completed and downloaded the same day (`run_info.json`: `labels = <CELL>/expr`,
20,209-20,214 genes, 679-760 drugs). Output:
`results/mc_phenotype_assoc/cross_cell/prism_clean28_conditional.tsv`.*

**An implementation gap found and fixed first.** The conditional columns as originally written computed the
partial correlation over **all ~680 drugs**, not within the matched subgroup - so they answered a different
question (does expression explain the pooled, confounded, near-null signal?) than the one registered in 10.7.
`r_partial_moa` / `r_expr_moa` / `r_partial_other` now compute it within the matched group and within the
complement separately. The numbers below are the registered statistic.

### 12.1 The result

`dntp_endpoints`, the 28 matched drugs, per cell pair:

| cell pair | r MetChange | r expression | **r MetChange \| expression** | complement (MC\|expr) |
|---|---|---|---|---|
| A375/A549 | 0.259 | +0.009 | **0.260** | -0.047 |
| A375/HT29 | 0.238 | -0.074 | **0.226** | -0.095 |
| A375/MCF7 | 0.502 | -0.031 | **0.520** | -0.036 |
| A375/PC3 | 0.094 | -0.093 | **0.081** | -0.064 |
| A549/HT29 | 0.349 | -0.078 | **0.380** | +0.018 |
| A549/MCF7 | 0.190 | -0.130 | **0.170** | -0.085 |
| A549/PC3 | 0.288 | +0.094 | **0.348** | -0.070 |
| HT29/MCF7 | 0.385 | +0.138 | **0.458** | +0.025 |
| HT29/PC3 | 0.169 | -0.103 | **0.136** | -0.023 |
| MCF7/PC3 | 0.434 | +0.185 | **0.441** | +0.014 |

| | mean r | pairs > 0 |
|---|---|---|
| MetChange | +0.291 | 10/10 |
| expression, `nucleotide` gene set | **-0.008** | 4/10 |
| **MetChange given expression** | **+0.302** | **10/10** |

**Retention: 104%.** The registered prediction (10.7) was that retention would stay well above 50%; it is at
100%. Paired Wilcoxon of the partial, matched against complement over the 10 pairs: 10/10 differences
favouring matched, p at the n=10 floor.

### 12.2 What this does and does not license

**It is the strongest result in the project for the central thesis**, and the first measurement of *increment*
rather than a side-by-side comparison (10.6). On five **non-HepG2** cell lines, with a mechanism-matched drug
set, MetChange's cell-specific prediction tracks cell-specific sensitivity and **the expression of the same
pathway's genes carries none of that information**.

**But read why retention is 104% rather than, say, 70%.** It is not that expression carries a strong signal
which MetChange survives - it is that **the gene-set differential is flat (-0.008, 4 of 10 pairs positive)**,
so conditioning on it removes nothing because there is nothing to remove. That is the cleanest possible form
of "the transcriptome alone does not see this", and it is a weaker statement than "MetChange adds to a strong
transcriptional signal". The comparator is fair - the `nucleotide` gene set is the GPR genes of the same
subsystems whose metabolites `dntp_endpoints` summarises - but the finding is an *absence* on the expression
side, not a defeat.

**The forking path is untouched by this.** The conditional statistic was registered in advance; **the 28-drug
set was not** - it is the product of five successive narrowings (11.4). So the right description is: *within a
post-hoc-selected drug set, a pre-registered conditional test behaves exactly as the thesis requires.* Only a
fresh cohort removes the remaining threat, and 12.3 now makes that much more feasible than it looked.

All p-values remain saturated at n = 10 pairs; quote the effect sizes.

### 12.3 93 more cell lines are available, and that is the real headline for scaling

> **WRONG, corrected in section 15.** The 93 figure comes from `cellinfo_beta.txt` membership, which says
> nothing about how many of the relevant drugs each cell line was actually profiled with. Measured against
> `lincs_unique_perts_bycell.tsv`, **only two or three cell lines are worth computing, not 93**: the non-core
> LINCS panels are modally 120 or 723 perturbations from a different chemical space and carry 5 or fewer of
> the 31 nucleotide-targeting drugs. Read section 15 instead of this one.


PRISM covers 568 cell lines; 99 of them appear in `cellinfo_beta.txt`; **6 have MetChange labels (A375, A549,
HEPG2, HT29, MCF7, PC3) and 93 do not.** Among the 93: HCT116, HUH7, JHH5, JHH7, DU145, BT474, 22RV1, HCC827,
HCC44, and a deep panel of endometrial, ovarian, bladder and colon lines.

That changes the limiting factor identified in 11.6. With 5 cells the aggregate test cannot resolve below
p = 0.000977 however strong the effect, because it has only 10 pairs. **With 20 cell lines it would have 190
pairs; with 50, 1,225.** The confirmatory cohort does not need a different dataset - it needs more label runs
on the dataset already in hand.

**The gating check is LINCS depth, not `cellinfo` membership** - the lesson from SHSY5Y, which appears in
`cellinfo` with 224 perturbations and only 10 Naia compounds (8.1). Before committing to label runs, count
distinct `trt_cp` `pert_id`s per candidate cell in `siginfo_beta.txt` and keep those with enough overlap with
the 821 PRISM drugs. HUH7 is worth particular attention: it is a liver line, and HEPG2 - where every other
result in this project lives - is absent from PRISM entirely (1.1).

## 13. The specificity question, tested properly: PRISM cannot answer it

*Raised as the objection that other subsystem-targeting drug groups might show the same association. It is the
right control. Answer: the library contains exactly one testable metabolic drug class, so specificity can
neither be established nor refuted here. Ten runs, one per subsystem, drugs selected by inhibited GEM target
and scored on the metabolite set of the same subsystem - the diagonal of a drug-group x set design.*

### 13.1 The table

| subsystem | matched drugs | overlap with the 28 | mean r | pairs > 0 | complement | Wilcoxon p | MC \| expr |
|---|---|---|---|---|---|---|---|
| nucleotide | 28 | 28 | **+0.254** | 10/10 | -0.058 | 0.00098 | +0.253 |
| pyrimidine | 20 | 20 | **+0.306** | 10/10 | -0.058 | 0.00098 | +0.281 |
| purine | 13 | 13 | **+0.152** | 8/10 | -0.058 | 0.0029 | +0.151 |
| amino_acid | 48 | 3 | -0.061 | 4/10 | -0.014 | 0.78 | -0.087 |
| folate | 4 | 4 | not scorable | | | | |
| glycolysis | 4 | 0 | not scorable | | | | |
| oxphos | 4 | 3 | not scorable | | | | |
| pyruvate | 4 | 0 | not scorable | | | | |
| tca | 3 | 0 | not scorable | | | | |
| **ppp** | **1** | 0 | not scorable | | | | |

### 13.2 The three positives are one finding, not three

`purine` (13 drugs) and `pyrimidine` (20) are **entirely inside** the 28-drug nucleotide group, and 27 of the
28 fall in one or the other. They are the same drug set re-sliced, so they do not corroborate each other. The
apparent breadth across three subsystems is an artefact of the subsystems overlapping in their drugs.

### 13.3 Six of the nine comparators have 1 to 4 drugs, and are contaminated with the same drugs

`ppp` has **one** drug. `tca` has three - diphenyleneiodonium (a NOX inhibitor), disulfiram (ALDH) and
hexachlorophene (an uncoupler) - none of which inhibits a TCA enzyme. `oxphos` has four, three of them
*already in the nucleotide 28*: brequinar and teriflunomide are DHODH inhibitors, and DHODH sits in the OXPHOS
subsystem because it feeds ubiquinone. `folate`'s four are the antifolates, also in the 28. `glycolysis`'s four
are a VEGFR inhibitor, an ALDH inhibitor, a NOX inhibitor and adaprev.

A repurposing library has dozens of clinical antimetabolites and essentially no approved drugs against GAPDH,
citrate synthase or complex I. That asymmetry, not the biology, is what this table mostly measures.

### 13.4 `amino_acid` looked like the comparator and is not one

It is the one group that clears the size bar - **48 drugs, larger than the nucleotide 28, and nearly
independent of it (3 shared)** - and it is null (-0.061, 4 of 10 pairs, p = 0.78). That is not a power
artefact. It is also **not interpretable**, because the group is almost entirely contamination:

| route into the group | n | examples |
|---|---|---|
| **CYP2C19** | 13 | imatinib, gefitinib, vemurafenib, axitinib, rifampicin, amoxicillin, propranolol, ketoconazole, flutamide... - drugs **metabolised by** it, not inhibitors of amino-acid metabolism |
| **MAO A/B** | 6 | iproniazid, piperine, bifemelane, vanoxerine, salidroside - monoamine catabolism, which sits in the tyrosine and tryptophan subsystems |
| **SAM-dependent methyltransferases** | 2 | BIX-01294 (EHMT1/2), E-7438 (EZH2) - epigenetic drugs, in cysteine/methionine metabolism because their targets consume SAM |
| sirtuins | 2 | tenovin-6, bisindolylmaleimide-ix (NAD-consuming) |
| AKR1C1, AOX1, TYR, NOS, MIF, ALDH, COMT, HPGDS | ~20 | aspirin, curcumin, disulfiram, raloxifene, monobenzone, nocodazole... |
| **genuinely amino-acid metabolism** | **~1** | **nitisinone** (HPD, the tyrosinaemia drug) |

So the null means "this group contains no amino-acid-directed drugs", not "amino-acid capacity is not
limiting".

**And the comparison was asymmetric in the favourable direction.** The nucleotide arm was cleaned from 44 to
28 by excluding non-mechanistic targets (11.3); this arm was left uncleaned at 48. Comparing a cleaned arm
against an uncleaned one is biased toward the nucleotide result, and the asymmetry was mine.

### 13.5 The methodological conclusion, which generalises beyond PRISM

**Target-based selection on GEM subsystems needs a curated enzyme whitelist, not a subsystem regex.**
Human-GEM's subsystems absorb, by correct bookkeeping:

- **drug metabolism** - every CYP sits in an amino-acid or xenobiotic subsystem, so any drug with an annotated
  CYP interaction qualifies, and "target" in these tables includes metabolising enzymes;
- **monoamine catabolism** - MAO in tyrosine/tryptophan, so every MAO inhibitor qualifies;
- **SAM-dependent methylation** - every methyltransferase, hence every epigenetic methyltransferase inhibitor,
  qualifies via cysteine/methionine metabolism;
- **NAD and nucleotide second messengers** - sirtuins and phosphodiesterases via purine metabolism (11.3);
- **phosphoinositide signalling** - PI3K via inositol phosphate metabolism (10.1).

Three separate times now the same failure has produced a plausible-looking drug group: PI3K under inositol
phosphate, PDE under purine, and CYP/MAO/methyltransferase under amino acid. **The subsystem a gene's
reactions belong to does not tell you whether a drug acting on it perturbs that pathway's flux.**

### 13.6 Where the specificity question stands

Unanswered, and not answerable with this library. PRISM offers 28 cleanable nucleotide-directed drugs, roughly
one genuine amino-acid-directed drug, and one to four per remaining subsystem. The result of 12.1 therefore
stands as *"for the one metabolic drug class this library contains, MetChange's cell-specific prediction works
and the matching gene expression does not"* - and cannot be sharpened to *"specifically nucleotide"* on this
evidence.

Two routes out, in order of cost:

1. **More cell lines, same drugs** (12.3): 93 PRISM lines are in LINCS without labels. This does not address
   specificity, but it fixes the saturated p-values and would make the nucleotide finding solid in itself.
2. **A library with broader metabolic coverage.** Tox21 and ToxCast have the mitochondrial and metabolic
   actives that a repurposing library lacks - and bring back the chemical-space mismatch with LINCS that
   produced the Hallinger and Xu nulls (5.1). Checking the usable n for a *cleaned, whitelisted* target
   selection against those libraries is the next concrete step, and must be done before any further
   analysis is designed around them.

## 14. A principled target filter, and why it is for the NEXT library rather than this one

*Implemented in response to 13.5. `scripts/mc_phenotype_assoc/build_target_relevance.py` (new) plus three
flags on `prism_differential.py`. **Deliberately not run-and-quoted on PRISM** - that would be an eighth look
at the same data (11.4).*

### 14.1 Producer direction: derived from stoichiometry, and it rederives the hand-made list

`build_target_relevance.py` emits `gene, set, hops`: the shortest number of reaction steps from a gene's
reactions to a member of a metabolite set, **following production only**, with currency metabolites excluded
from the traversal for the same reason the network builder removes them.

For `dntp_endpoints` it rules out **1,509 of 2,522 GEM genes (60%) as never reached** - and the never-reached
set includes **PDE, PI3K and EZH2**, i.e. exactly the contaminations that 10.1, 11.3 and 13.4 found by
inspection. The hand-made exclusion list is recovered from the stoichiometry.

**Two things it cannot do, measured rather than assumed.**

*Hop distance does not separate biosynthetic routes from incidental ones.* Distances to a dNTP:

| hop | genes include |
|---|---|
| 1 | NME1/2, RRM1/2 |
| 4 | TK1, TYMS, **EHMT2** (histone methyltransferase) |
| 5 | DHFR, IMPDH1, TYMP, ADA, **CYP2C19** |
| 6 | DHODH, **MAOA** |
| never | PDE*, PIK3*, EZH2 |

Any cutoff including DHODH includes MAO-A; any including DHFR includes CYP2C19; any including TYMS includes
G9a. Gene counts: 19 at hop 1, 43 at 2, 117 at 3, 402 at 4, 837 at 5, 1,415 at 6.

*It only discriminates for small, specific sets.* Drugs selected by inhibited-target reachability:

| set | H<=2 | H<=3 | H<=4 | H<=6 | H<=6 and n_sources>=3 |
|---|---|---|---|---|---|
| **dntp_endpoints** | 8 | 14 | 42 | 116 | 45 |
| nucleotide (192 metabolites) | 235 | 235 | 235 | 235 | 105 |
| glycolysis / tca / oxphos / ppp | ~233 | ~233 | ~233 | ~233 | ~104 |

The subsystem sets hold 192-823 metabolites, many near-currency, so nearly every GEM gene is within a few
hops of producing *something* in them. **The filter is therefore only meaningful with endpoint sets** - which
is a design constraint worth accepting, since `dntp_endpoints` (4 metabolites) is the right object for this
question and `nucleotide` (192) is not.

### 14.2 Two different failure modes, and a correction to 11.3

| failure | example | fix |
|---|---|---|
| **enzyme is in the GEM for non-biosynthetic reasons** | PDE, CYP2C19, MAO-A, EHMT2 | class exclusion (14.3) |
| **the drug-target annotation is weak or wrong** | adaprev->HK1, tofogliflozin->GCK, PJ-34->EEF2 | `--min-target-sources` |

The three implausible annotations all carry **n_sources = 2** from ChemPert-CLUE and RepurposingHub, while
every genuine antimetabolite annotation carries **3-5** including ChEMBL and DrugCentral (methotrexate->DHFR 5,
gemcitabine->RRM1 5, floxuridine->TYMS 4, brequinar->DHODH 4, merimepodib->IMPDH1 3). So `n_sources >= 3`
removes all three.

**Correction: GCK, HK1 and EEF2 should not have been in the 11.3 exclusion list.** They are misplaced
*annotations*, not misplaced enzymes - a genuine hexokinase inhibitor would be a legitimate glycolysis drug,
and excluding HK1 as a class would discard it. 11.3's result is unaffected (the three drugs fail the source
threshold anyway) but the reasoning was wrong.

### 14.3 The class exclusion list, with its principle

Enzymes present in the GEM for **catabolic, cofactor-consuming, signalling or xenobiotic** reasons rather than
as biosynthetic routes to a precursor. Marked where the producer filter already handles it:

| class | pattern | handled by 14.1? |
|---|---|---|
| cyclic-nucleotide hydrolases | `^PDE` | **yes**, never reached |
| lipid / inositol kinases | `^PIK3`, `^PIP[45]K`, `^SPHK`, `^DGK` | **yes**, never reached |
| xenobiotic metabolism | `^CYP`, `^UGT`, `^SULT`, `^NAT[12]$`, `^FMO` | no - CYP2C19 at hop 5 |
| SAM-dependent methyltransferases | `^EHMT`, `^EZH`, `^DNMT`, `^SETD`, `^KMT`, `^PRMT`, `^COMT$`, `^NNMT$` | no - EHMT2 at hop 4 |
| monoamine / aldehyde catabolism | `^MAO`, `^ALDH`, `^ADH[0-9]`, `^AOX1$`, `^DBH$` | no - MAOA at hop 6 |
| NAD-consuming signalling | `^SIRT`, `^PARP`, `^CD38$` | no |
| redox / detoxification | `^TXNRD`, `^GSR$`, `^GPX`, `^GST`, `^NQO1$`, `^SOD[0-9]`, `^CAT$` | no |
| eicosanoids | `^PTGS`, `^ALOX`, `^PTGES`, `^HPGDS$`, `^TBXAS1$` | no |
| nitric oxide | `^NOS[0-9]` | no |
| carbonic anhydrases | `^CA[0-9]` | no |
| steroidogenesis | `^HSD`, `^SRD5A`, `^STS$` | no |

**Deliberately absent, and why:**

- **NAMPT** - NAD *synthesis* is biosynthetic, so a NAMPT inhibitor genuinely limits a cofactor pool.
- **Transporters** (`SLC*`, `ABC*`, `ATP1A/B*`, `AQP*`) - context-dependent rather than wrong. SLC7A11
  (cystine) and SLC1A5 (glutamine) genuinely limit precursor supply, so a blanket exclusion would discard real
  signal. Flag and decide per question.
- **GCK, HK1, EEF2** - see 14.2; these are annotation failures.

### 14.4 The cost, and why this is for the next library

Applying the whole stack to PRISM - subsystem match, class exclusion, `n_sources >= 3`, producer reachability
at hop <= 4 - leaves **25 drugs overall and 12 for `dntp_endpoints`**, against 28 for the hand-cleaned group of
11.3. The hop-4 cutoff drops the DHFR, IMPDH and DHODH drugs, which sit at 5 and 6.

So the clean filter **halves the power** on the one set where the effect was visible. And it cannot help with
the specificity question at all, because **a better filter cannot create drugs**: PRISM still has one PPP
drug, three TCA drugs and four each for glycolysis and OXPHOS (13.1). Cleaning only shrinks groups further.

**Hence the filter is not run on PRISM here.** It is specified now so that it can be applied *before* looking
at results on a library where the drugs exist - Tox21 or ToxCast, which carry the mitochondrial and metabolic
actives a repurposing collection lacks, and which bring back the chemical-space mismatch with LINCS that
produced the Hallinger and Xu nulls (5.1). The order that follows from this: **check usable n under this
filter stack on those libraries first, before designing any analysis around them.**

Usage, for the record:

```bash
python3 -I scripts/mc_phenotype_assoc/build_target_relevance.py \
    --out results/mc_phenotype_assoc/metabolite_sets/target_relevance.tsv

python3 -I scripts/mc_phenotype_assoc/prism_differential.py \
    --targets data/drug_targets/combined/drug_targets.tsv \
    --nodes data/networks/combined/nodes.tsv \
    --target-subsystems 'nucleotide|purine|pyrimidine|folate' \
    --exclude-target-classes '<the patterns in 14.3>' \
    --min-target-sources 3 \
    --relevance results/mc_phenotype_assoc/metabolite_sets/target_relevance.tsv --max-hops 4 \
    --expr-run prism_viability_expr
```

### 14.5 CORRECTION: use the class exclusion alone; the producer filter is a diagnostic, not a selector

*14.1-14.4 recommended the producer-reachability filter as the primary mechanism. That was wrong, and the
objection that identified it is the right one: if the filter halves the power on drugs that genuinely do
affect dNTP production, it is doing harm.*

**Why the hop cutoff fails, stated properly.** DHFR, IMPDH and DHODH *are* on producing routes to dNTPs -
DHFR makes the THF that TYMS needs for dTMP, IMPDH runs IMP to XMP to GMP to dGTP, DHODH feeds orotate into
UMP and onward to dCTP and dTTP. They sit at hops 5-6 because **real biosynthetic routes are long**, not
because they are irrelevant. Hop distance measures path length, and path length is not relevance.

**The decisive comparison.** With PDE and the lipid kinases moved *into* the class list instead of being
delegated to the producer filter:

| rule | drugs | retains the 28 of 11.3 | extras |
|---|---|---|---|
| **subsystem + class exclusion** | **33** | **all 28** | PJ-34, adaprev, tofogliflozin (bad annotations); curcumin, diphenyleneiodonium (both via XDH) |
| subsystem + class + producer-reachable, any depth | 31 | **27** - loses ABT-702 | 4 |
| subsystem + producer-reachable, no class | 33 | 27 | |
| subsystem + class + `n_sources >= 3` | 25 | **19** - loses carmofur, cladribine, AVN-944 | |

**The class exclusion alone reproduces the hand-cleaned group exactly.** The producer filter adds nothing and
costs ABT-702, a genuine adenosine kinase inhibitor. And `n_sources >= 3` is too aggressive - it discards
carmofur, cladribine and AVN-944, which are real antimetabolites - so it belongs as a sensitivity check, not
as part of the rule.

**The revised rule, and it is simpler than 14.3-14.4 proposed:**

1. target in a matching GEM subsystem, `effect == inh`;
2. **exclude the classes of 14.3, with `^PDE`, `^PIK3|^PIP[45]K|^SPHK|^DGK` and `^XDH` added to the list**
   rather than delegated elsewhere - XDH is purine catabolism and is how curcumin and diphenyleneiodonium
   enter;
3. nothing else. No hop cutoff, no source threshold in the primary rule.

That leaves 31 drugs, retaining all 28, with three known bad annotations as residue - about 10% noise, which
`--min-target-sources 3` can be used to quantify as a sensitivity check without being baked in.

**What the producer filter is still for.** It is a **diagnostic that derived the class list**, not a selector.
Running it is what established that PDE, PI3K and EZH2 are *never* on a producing route to a dNTP - which is
the justification for putting them in the class list in the first place, replacing three rounds of inspection
(10.1, 11.3, 13.4) with a stoichiometric criterion. `build_target_relevance.py` and `--relevance` stay for
that purpose, and for checking a new class before adding it. They should not be used to choose drugs.

## 15. Which cell lines are actually worth computing: two, not ninety-three

*Measured from `data/LINCS/lincs_unique_perts_bycell.tsv` (227 cell lines, 201,295 rows) against the PRISM
dose-response table, before committing to MetChange label runs. This corrects 12.3.*

### 15.1 The measurement

89 cell lines appear in both the per-cell LINCS table and PRISM. Coverage of the **31** nucleotide-targeting
drugs (the class-excluded rule of 14.5) is what decides whether the section-12 finding can be tested in a
given cell:

| nucleotide drugs covered (of 31) | cell lines |
|---|---|
| 29-30 | **5** - PC3, A375, MCF7, A549, HT29: the cores, labels already computed |
| **24** | **1 - YAPC** |
| **18** | **1 - MDAMB231** |
| 11 | 1 - U2OS |
| <= 6 | **81** |

| cell | LINCS panel | PRISM overlap | nucleotide drugs | labels |
|---|---|---|---|---|
| PC3 / A375 / MCF7 / A549 / HT29 | 12,156-13,892 | 861-927 | 29-30 | yes |
| **YAPC** | **1,716** | **605** | **24** | no |
| MDAMB231 | 1,357 | 480 | 18 | no |
| U2OS | 20,836 | 438 | 11 | no |
| HEC108, CAL29, BEN, OVTOKO, HEC1A, ... | ~723 | ~220 | 5 | no |

**The cause is panel size and composition, not cell-line availability.** Non-core LINCS panels are modally
**120 perturbations (13 cell lines) or 723-724 (17 cell lines)** - a small fixed plate - against 12,000+ for
the five cores. U2OS is the instructive exception: 20,836 perturbations and still only 11 of the 31, because
its panel is a different chemical space, not a smaller version of the same one.

### 15.2 What each addition buys

The binding constraint identified in 11.6 was the saturated Wilcoxon, which cannot resolve below 1/2^(n pairs):

| adding | cells | pairs | Wilcoxon floor |
|---|---|---|---|
| nothing | 5 | 10 | 9.8e-04 |
| **+ YAPC** | **6** | **15** | **3.1e-05** |
| + MDAMB231 | 7 | 21 | 4.8e-07 |
| + U2OS | 8 | 28 | 3.7e-09 |

**YAPC alone breaks the saturation**, and at 24 of 31 drugs its per-pair correlations will be comparable in
quality to the existing five. MDAMB231 at 18 is a reasonable second. U2OS at 11 would add pairs but with
noticeably noisier per-pair estimates, and a cell carrying 5 of 31 drugs yields a correlation on five points,
which is not worth a label run.

**Recommendation: compute YAPC, then MDAMB231 if the first is informative. Do not compute the other 86.**

Note also that YAPC is one of the seven Panorama GCP cell lines (section 2), so a label run there serves two
analyses rather than one.

### 15.3 The error, which has now happened five times

12.3's "93 cell lines available" was derived from `cellinfo_beta.txt` membership and never checked against
per-cell drug coverage. That is the same mistake as:

- Hallinger - identified tier-1 compounds, almost none in LINCS (report 07);
- naia - 683 whole-LINCS structure matches, **10** in the SHSY5Y panel (report 09, 8.1);
- xu - 53 whole-LINCS matches, **11 vs 12** in the HEPG2 panel (report 10, 5.1);
- the PRISM subsystem comparators - 1 to 4 drugs each, and contaminated (13.1, 13.3);
- and now this.

The rule recorded after the second instance was *"a validation dataset's usable n is its overlap with the
specific label panel, computed before the dataset is ranked"*. It applies to **cell lines** exactly as it does
to datasets, and the general form is: **count the units that can actually carry the prediction, for the
specific comparison being proposed, before quoting an opportunity.** `lincs_unique_perts_bycell.tsv` makes
that a one-command check for any future cell-line question and should be the first thing consulted.
