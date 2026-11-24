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

and ask whether the two correlate across drugs. **Differencing removes the drug's general potency entirely** -
a compound that is simply toxic everywhere contributes zero to both sides - so the confounder is cancelled by
construction rather than matched away. What remains is strictly the *selectivity* question: when a drug hits
one line harder than another, do that line's MetChange labels say so?

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
Resolution order, first that works winning: an exported `$py`; the current `python3`; the micromamba
environment (its interpreter located by asking `micromamba run`, so no layout is guessed, with the shell hook
as a fallback when the environment needs its own `LD_LIBRARY_PATH`); then `$mamba_env_prefix/bin/python3` or
the usual micromamba roots. `project_paths.conf` gained an optional `mamba_env_prefix` for the case where
micromamba is not on `PATH` inside batch jobs. `provenance` now prints the interpreter and how it was found,
so a log says which Python produced the numbers.

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
