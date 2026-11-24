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
