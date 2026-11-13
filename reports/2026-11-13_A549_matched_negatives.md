# Report 01 — Matched negatives on A549 MetChange, and a dose/time leak

*2026-11-13. Runs: `results/predict_targets/A549/metchange_rank{,_shuff,_matched}` (5-fold CV by drug, active
perturbations, MetChange `rank` signatures, FNN ensemble without VAE).*

## 1. The question

The FNN ranks a drug's targets well, but almost as well when every perturbation's signature is replaced by another
one's (the shuffled control). So most of its performance is a **per-protein prior**: it has learned which proteins
tend to be drug targets (hubs, highly expressed, kinases...), and that knowledge does not depend on the signature.

We tried one fix, **matched negatives** (version (a) of the plan): for each training row "protein g is a target in
perturbation p", add rows "protein g is *not* a target in perturbation q", with q from another drug. Within such a
group every per-protein feature is identical, so the model can only separate the rows through the signature.

Two questions:

1. Does the matched model use the signature (real better than its shuffled control)?
2. Is it better than the unmatched model at that, and what does it cost?

## 2. Runs compared

| run | folder | matched negatives | epoch selection | signatures |
|---|---|---|---|---|
| unmatched | `metchange_rank` | 0 | validation AUPRC | real |
| unmatched, shuffled | `metchange_rank_shuff` | 0 | validation AUPRC | shuffled |
| matched | `metchange_rank_matched` | K = 4, weight 0.5 | validation `spec` | real |
| matched, shuffled | `metchange_rank_matched_shuff` (on the cluster) | K = 4, weight 0.5 | validation `spec` | shuffled |

Same seed, folds, drugs, network and features otherwise. Note that the matched runs change **two** things at once
(sampling and epoch selection); section 7 comes back to this.

## 3. The metrics, and why there are several

All are computed on **held-out drugs** (drugs never seen in training, in the test fold), as a *percentile*: the rank
of the correct answer divided by the number of candidates. **Lower is better, 0.5 is chance.** `topX_pct` is the share
of cases where the correct answer is in the top X%.

| metric | question it asks | what a pure per-protein prior gets |
|---|---|---|
| `rank_prob` | Given perturbation p, where does its own drug rank among all held-out drugs? (the repurposing use case) | **better than 0.5**: it ranks drugs with "popular" targets high for every perturbation |
| `spec_prob` | Given drug d, does its own perturbation p score d's targets higher than other drugs' perturbations score them? | 0.5: the same protein gets the same score everywhere, so p cannot beat q |
| `rank_probc` | `rank_prob` after standardising each candidate drug's scores over the test perturbations | 0.5: the per-drug offset is removed |
| `rank_baseline` | ranking by how often each target appears among training drugs (no model) | ~0.5 here |
| `*_logreg` | the same for the logistic regression | — |

**Why two families.** `rank_prob` is what we want in the end, but it mixes prior and signature. `spec` and
`rank_probc` were added to isolate the signature: they cancel anything that is *constant per protein or per drug*.
As section 5 shows, they do **not** cancel things that vary *per perturbation* without being the signature (dose,
time). "Prior-free" is not the same as "leak-free".

## 4. How to read a comparison (the checks, in order)

1. **Against the run's own shuffled control, not against 0.5.** The shuffled control keeps everything except the
   link between a perturbation and its signature: same features, training, prior, dose, time. Whatever the shuffled
   run achieves, the model achieves *without* the signature. The signal attributable to the signature is
   *real − own control*.
2. **Paired per drug.** Both runs score the same held-out drugs, so compare drug by drug (median difference,
   fraction of drugs better/worse, Wilcoxon signed-rank test). Drugs, not perturbations, are the independent units:
   the doses of one drug share targets and much of their signature, so counting them separately inflates n and
   makes p-values far too small (compare the per-perturbation p = 8e-34 with the per-drug p = 0.02 below).
3. **Per fold.** A real effect should point the same way in most of the 5 folds; a pooled difference driven by one
   fold is fragile.
4. **Look at the control itself.** A control far from 0.5 on a metric meant to be at 0.5 is not a nuisance, it is a
   finding: something other than the signature is carrying information (section 5).
5. **Training diagnostics.** The epoch each ensemble member kept (`per_fold_best_epoch`) tells whether training is
   stable; erratic choices (epoch 1 in one fold, 26 in another) make every comparison noisier.

Tools: `scripts/predict_targets/compare_runs.py` (checks 1–3; output in
`results/predict_targets/A549/mc_matched_vs_unmatched/`) and `check_dose_time_leak.py` (section 5).

## 5. Results

### 5.1 Each run against its own shuffled control

Median percentile, active perturbations / per drug (1,358 held-out drugs, pooled over folds):

| metric | unmatched | unmatched shuffled | matched | matched shuffled |
|---|---|---|---|---|
| `rank_prob` | **0.305** / **0.409** | 0.317 / 0.418 | 0.394 / 0.432 | 0.380 / 0.440 |
| `spec_prob` | 0.447 / 0.471 | **0.464 / 0.464** | 0.433 / 0.443 | **0.495 / 0.496** |
| `rank_probc` | 0.404 / 0.441 | 0.426 / 0.429 | 0.419 / 0.411 | 0.484 / 0.476 |
| `rank_logreg` | 0.403 / 0.468 | 0.401 / 0.470 | 0.454 / 0.483 | 0.453 / 0.482 |

Paired per drug, real vs own control (Wilcoxon p):

| metric | unmatched vs its shuffle | matched vs its shuffle |
|---|---|---|
| `rank_prob` | −0.011, **p = 0.007** | +0.009, p = 0.59 |
| `spec_prob` | −0.011, p = 0.49 | −0.045, **p = 0.001** |
| `rank_probc` | −0.022, p = 0.32 | −0.019, **p = 0.0005** |

Per fold, `spec_prob` per perturbation:

| fold | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| unmatched | 0.472 | 0.418 | 0.441 | 0.484 | 0.414 |
| unmatched shuffled | 0.487 | 0.462 | 0.477 | 0.446 | 0.455 |
| matched | 0.402 | 0.406 | 0.466 | 0.459 | 0.429 |
| matched shuffled | 0.500 | 0.473 | 0.498 | 0.497 | 0.509 |

### 5.2 Matched against unmatched directly

Paired per drug: `rank_prob` matched worse (+0.016, p = 0.023; per perturbation 0.305 → 0.394); `spec_prob`
−0.011, p = 0.21; `rank_probc` −0.015, p = 0.29. Per-drug `spec` correlates only 0.23 (Spearman) between the two
runs: *which* drugs look specific changes a lot from run to run.

### 5.3 The unmatched shuffled control is below 0.5 on `spec`: dose and time

The unmatched shuffled control scores `spec` 0.446–0.487 in all five folds. A shuffled model cannot distinguish
perturbations through their signature, so it must do it through the other per-perturbation inputs: relative dose,
time, and signature size. Signature size moves with the shuffled signature, which leaves **dose and time**.

`check_dose_time_leak.py` tests this on the saved held-out scores (`test_scores_heldout.npz`). For each perturbation
it replaces the model's scores by the average scores of *other drugs'* perturbations at the same time and relative
dose, a predictor that knows nothing but dose and time, and computes `spec` of that. (Its `spec(model)` column
reproduces the cluster's per-fold numbers exactly, a check that the scores are read correctly.)

| run | `spec(model)` per fold | `spec(dose/time only)` per fold | variance explained by dose/time |
|---|---|---|---|
| unmatched shuffled | 0.487 0.462 0.477 0.446 0.455 | **0.443 0.495 0.441 0.432 0.433** | 0.18–0.26 |
| unmatched | 0.472 0.418 0.441 0.484 0.414 | 0.437 0.422 0.449 0.515 0.424 | 0.16–0.27 |
| matched | 0.402 0.406 0.466 0.459 0.429 | 0.444 0.414 0.473 0.434 0.477 | 0.16–0.28 |
| matched shuffled | 0.500 0.473 0.498 0.497 0.509 | 0.466 0.459 0.488 0.471 0.498 | 0.20–0.27 |

In the shuffled run, knowing only dose and time gives `spec` ≈ 0.44, as good as the model itself, and dose/time
groups explain about a fifth of the score variance. **The model's scores carry drug identity through dose and
time.**

How this can happen: in LINCS the experimental design is not random with respect to the drug. Different screens
used different designs (dose-response series at 6 h and 24 h, single 10 µM doses, 48 h time points), and the
screens differ in which compound classes they contain. The design can therefore say something about the drug
*class*, and so about its likely targets. (Which screens and classes drive it here has not been checked; a
breakdown of the leak by dose/time design would show it.) The FNN can use it because its inputs include products of features (the degree-2
polynomial expansion): `rel_log_dose × protein feature` lets the dose design shift the scores of specific proteins.
The logistic regression, which only adds features, cannot: a dose/time term shifts all proteins equally and cancels
in `spec`. That is consistent with its `spec_logreg` being at 0.49–0.51 in every run.

The matched shuffled run shows the difference between *depending on* dose/time and *exploiting* it. Its scores
depend on dose/time as much as the unmatched control's (R² ~0.2), but that dependence carries much less drug
identity: the dose/time-only predictor reaches 0.46–0.50 instead of ~0.44, and the model's own `spec` is 0.50. Most
of the dependence is presumably generic (e.g. 6 h perturbations scoring lower for every protein), which cancels in
`spec`. The faint remainder (~0.47) only shows after averaging perturbations within a dose/time group, which removes
noise and so reveals effects too weak to move the per-perturbation metric.

In the real runs signature and design are correlated too (24 h signatures differ from 6 h ones, high doses give
stronger signatures), so the dose/time-only predictor cannot separate leak from signal there. That is what the
shuffled control is for.

### 5.4 `spec` within dose / time strata

To take dose and time out of the metric, `compare_runs.py` now recomputes `spec` from the saved scores with each
perturbation compared only to other drugs' perturbations **at the same time point** (`spec_prob_time`) or **at the
same time and relative dose** (`spec_prob_timedose`). Within a stratum, dose and time are identical, so they cannot
help. Two checks came before trusting it:

- With no stratum, the recomputation reproduces the cluster's `spec_prob_heldout` exactly (3,970 / 3,970).
- **Chance level of a rank with few competitors.** The first version used the same formula as the training script,
  (1 + #higher + ½·#ties) / (n + 1). Under chance its average is (1 + n/2) / (n + 1): 0.5006 for the ~800
  competitors of the unstratified `spec`, but 0.58 for n = 5, which happens within time × dose strata. The shuffled
  controls showed it (0.517 and 0.512 instead of 0.5). The stratified version uses the mid-rank
  (#higher + ½·#ties + ½) / (n + 1), which is 0.5 on average for any n.

Perturbations with fewer than 5 competitors in their stratum are left out (3,511 of 3,970 kept for time × dose).
Results, each run against its own shuffled control:

| run | stratum | per perturbation, real / shuffled | folds real better | per drug, real / shuffled | per-drug p |
|---|---|---|---|---|---|
| unmatched | none | 0.447 / 0.464 | 4/5 | 0.471 / 0.464 | 0.49 |
| unmatched | time | 0.445 / 0.472 | 4/5 | 0.495 / 0.503 | 0.019 |
| unmatched | time × dose | 0.467 / 0.497 | 3/5 | 0.482 / 0.495 | 0.60 |
| matched | none | 0.433 / 0.495 | 5/5 | 0.443 / 0.496 | 0.001 |
| matched | time | 0.437 / 0.492 | 5/5 | 0.479 / 0.495 | 0.082 |
| matched | time × dose | **0.457 / 0.504** | **5/5** | 0.485 / 0.497 | 0.16 |
| noStatic | time × dose | 0.500 / 0.500 | 3/5 | 0.495 / 0.494 | 0.50 |

- **All shuffled controls are at 0.50 within time × dose**: the stratified metric is blind to the design, as
  intended, and the unmatched control's 0.464 in 5.1 was indeed dose / time.
- **The matched model keeps most of its advantage**: 4.7 percentile points per perturbation (6.2 unstratified), in
  all five folds.
- **Per drug the effect is much weaker** (p = 0.16). The per-drug value is the mean of the drug's per-perturbation
  percentiles, which stays near 0.5 unless the drug's perturbations are consistently good. A gap that is clear per
  perturbation but weak per drug means the signal comes from some perturbations of some drugs, not from most drugs.
  Per-perturbation p-values would look highly significant, but they would count the doses of one drug as independent
  evidence; per drug is the honest test.
- **The unmatched model also has a small signature effect** (3 points per perturbation within time × dose), which
  was hidden in 5.1 because its control was inflated by dose / time. It is not significant per drug.
- **noStatic shows no signature specificity at all**, although it is the run where the signature improves plain
  ranking most (`rank_prob` real vs shuffled 0.427 vs 0.472 per drug, p = 0.001). Without the protein features the
  signature helps rank drugs in a way that is the same for every perturbation (e.g. strong signatures pointing to
  drug classes with many targets), not in a way specific to the drug.

## 6. Interpretation

1. **The unmatched model's apparent specificity is mostly dose/time, not signature.** Its `spec` (0.447) is no
   better than its own shuffled control (0.464), and a dose/time-only predictor matches it. Within dose / time
   strata a small signature effect appears (3 points per perturbation, 5.4), not significant per drug.
2. **The matched model does use the signature, but the evidence is per perturbation more than per drug.** It is
   better than its own control on `spec` (−4.5 points per drug, p = 0.001) and `rank_probc` (p = 0.0005), in all
   five folds, and its shuffled control sits at 0.50. Within time × dose strata it keeps 4.7 points per perturbation
   in all five folds, but per drug the difference is small (p = 0.16): the signal is real but concentrated in part of
   the drugs. The effect is modest: 14% vs 11% of perturbations put their own drug in the top 10%.
3. **It costs plain ranking.** Moving half of the negative class's weight to the matched rows weakens the prior:
   `rank_prob` 0.305 → 0.394 per perturbation. This is the expected trade-off; the plan's answer is a
   prior + residual model, or keeping the prior in the cross-entropy and adding specificity as a separate loss term
   (version (b)).
4. **Comparing matched and unmatched directly (5.2) was misleading.** It suggested "no difference in specificity",
   because the unmatched number was inflated by dose/time. The comparison that answers "does it use the signature"
   is each run against its own control. (An earlier reading in this project, that the unmatched model "already reads
   the signature", came from exactly this mistake.)

## 7. Open points and caveats

- **Why the matched shuffled control has no leak is not settled.** The matched runs changed both the sampling and
  the epoch selection (`spec` instead of AUPRC). One explanation is that selecting on validation `spec` keeps early
  epochs, before the model has learned the dose/time shortcut. Consistent with it: the matched shuffled ensemble
  kept epochs 2–5 in most folds and members (once 10 and 13), against 1–24 in the unmatched shuffled one. Not
  proof: the 2×2 runs in section 8 separate the two changes.
- **Validation scores used for selection are optimistic.** In the matched shuffled run, validation `spec` of the
  kept epochs is 0.47–0.50 although the signatures carry nothing (test: 0.50). Picking the best of many noisy epochs
  makes the winner look better than chance (the *winner's curse*). This is why judgements here are made on the
  held-out test folds only, never on the validation metric that chose the epoch.
- **Per-drug results are unstable** (correlation 0.23 between runs). The signature effect is real on average but
  should not yet be read drug by drug.
- **Training is still erratic.** Kept epochs range from 1 to 26 across folds and ensemble members in both runs;
  validation `spec` per member varies from 0.40 to 0.49.
- **Leak-robust evaluation** is now in `compare_runs.py` (5.4). The training script's own `spec` (and the validation
  `spec` used for epoch selection) is still unstratified, so it still rewards dose / time.
- **Relative dose is approximate.** It is computed from the perturbation names, against the drug's highest *scored*
  dose; the model uses the highest dose in the perturbation table. Differences only matter for drugs with inactive
  doses left out by `--score-active-only`.

## 8. Next steps

1. ~~Run `check_dose_time_leak.py` on `metchange_rank_matched_shuff`~~ (done, section 5.3: no exploited leak).
2. ~~Add the dose/time-stratified `spec` to the evaluation~~ (done, 5.4; `results/.../A549/stratified_spec/`).
3. On the cluster: the unmatched pair with `--no-globals` (real + shuffled). Note: `--no-static` (the noStatic runs)
   removes the per-protein features, not dose / time; no `--no-globals` run exists yet. If its shuffled `spec`
   returns to 0.5 unstratified, the model was using the dose / time *inputs*; if not, the design also reaches it
   through the signature itself.
4. Separate the two changes of the matched runs: matched with `--select-metric auprc`, and/or unmatched with
   `--select-metric spec`, each with its shuffled control.
5. Then version (b): keep the cross-entropy on the original sampling (the prior stays intact) and add the
   specificity loss on matched groups, weighted by λ.

## Glossary (DL concepts that came up)

- **Negative sampling.** With ~19,000 proteins per perturbation and a handful of targets, training on every row
  would be dominated by negatives; a subset is sampled. *Which* negatives are sampled decides what the model learns:
  random negatives teach "which proteins are targets at all" (the prior); matched negatives teach "in which
  perturbations is this protein a target" (the signature).
- **Shortcut learning.** A model uses the easiest feature that predicts the label in the training data, whether or
  not it is the one we care about: here the protein prior, then dose/time. High test performance does not show that
  the model learned the intended thing; controls do.
- **Leakage.** Information about the answer reaching the model through a route that will not exist, or is not
  wanted, at application time. Dose/time are available for new compounds, so this is a *confound* rather than a
  classic leak, but it inflates the metric meant to measure signature use.
- **Negative control (shuffled signatures).** Retrain with the causal link broken (signatures permuted between
  perturbations) and keep everything else. What survives did not come from the signature.
- **Model selection / early stopping.** Training runs for many epochs and the epoch with the best validation score
  is kept. The validation metric therefore shapes the final model as much as the loss does: selecting on AUPRC
  rewards the prior; selecting on `spec` rewards signature use.
- **Ensemble.** Three FNNs with different batch sizes / regularisation, probabilities averaged. Reduces variance
  between members, which matters here because individual members are unstable.
- **Paired tests and units of analysis.** When two models score the same items, compare item by item; choose as
  items the independent units (drugs), not repeated measurements of them (doses).
