# Report 02 — Training stability was the bottleneck, not the sampling

*2026-11-14. Runs: `results/predict_targets/A549/metchange_rank_{matched_auprc,noGlobals_lr1e-3,matched_noGlobals_lr1e-3}`
and their `_shuff` controls (5-fold CV by drug, active perturbations, MetChange `rank` signatures, FNN ensemble
without VAE). Follows report 01 (`2026-11-13_A549_matched_negatives.md`), whose metric definitions and reading
checklist apply here unchanged.*

> **Correction (report 03, `2026-11-14_seed_tf_expr.md`).** The "points" quoted below (8.6 per perturbation, 5.6 per
> drug) are differences between the medians of the real and shuffled runs. The per-drug Wilcoxon p-values test the
> paired per-drug differences, whose median (`median_diff` in `paired.tsv`) is **4.2 points** for `spec_prob_timedose`
> (p = 4.5e-5), 6.8 for `rank_probc` and 1.1 for `rank_prob`. Conclusions unchanged. Report 03 also shows that a second
> seed gives a smaller effect (2.9 points, p = 0.047) and that the 60-epoch cap did not bind.

## 1. The question

Report 01 left three things open:

1. The matched runs changed sampling **and** epoch selection at once → `mc_matched_auprc` separates them.
2. `mc_noGlobals` showed a signature effect of only ~2 points, and **training stopped after 1–3 epochs in three of
   five folds**. Was the small effect the data, or the training? → `mc_noGlobals_lr1e-3`.
3. Does the matched sampling still help once dose/time inputs are gone and training runs properly? →
   `mc_matched_noGlobals_lr1e-3`.

## 2. Runs

| run | matched | dose/time inputs | epoch selection | lr | patience | epoch cap |
|---|---|---|---|---|---|---|
| `metchange_rank_matched_auprc` | K = 4, w 0.5 | yes | AUPRC | 0.005 | 5 | 30 |
| `metchange_rank_noGlobals_lr1e-3` | no | **no** | AUPRC | **0.001** | **10** | **60** |
| `metchange_rank_matched_noGlobals_lr1e-3` | K = 4, w 0.5 | **no** | `spec` | **0.001** | **10** | **60** |

Each with its own `--shuffle-signatures` control. The `lr1e-3` runs are array jobs, one fold per task; same splits
and seeds as the sequential runs. For reference, from report 01: `metchange_rank` (unmatched, dose/time inputs,
lr 0.005), `metchange_rank_noGlobals` (same at lr 0.005), `metchange_rank_matched` (matched, dose/time inputs,
selection on `spec`, lr 0.005).

## 3. Headline result

**Training was the limiting factor.** The same unmatched model that showed 2 points of signature specificity at
lr 0.005 shows 8.6 points at lr 0.001, in all five folds, with per-drug p = 4.5e-5.

`spec_prob_timedose` (stratified by time × dose, the leak-free metric), real / shuffled:

| run | per perturbation | folds real better | per drug | per-drug p |
|---|---|---|---|---|
| `mc` (report 01) | 0.467 / 0.497 | 3/5 | 0.482 / 0.495 | 0.60 |
| `mc_matched` (report 01) | 0.457 / 0.504 | 5/5 | 0.485 / 0.497 | 0.16 |
| `mc_noGlobals` (report 01) | 0.480 / 0.500 | 3/5 | 0.486 / 0.498 | 0.057 |
| `mc_matched_auprc` | 0.450 / 0.497 | 4/5 | 0.474 / 0.495 | 0.45 |
| **`mc_noGlobals_lr1e-3`** | **0.414 / 0.500** | **5/5** | **0.444 / 0.500** | **4.5e-5** |
| `mc_matched_noGlobals_lr1e-3` | 0.441 / 0.504 | 5/5 | 0.462 / 0.503 | 1.3e-3 |

Per fold, `mc_noGlobals_lr1e-3`: real 0.410 0.417 0.440 0.397 0.408; shuffled 0.496 0.498 0.509 0.511 0.500. The
shuffled controls sit at 0.50 as they should (no dose/time inputs), so for these runs the plain `spec` and the
stratified `spec` agree.

**Plain ranking improved at the same time**, which the earlier trade-off did not suggest. `rank_prob` per
perturbation, real / shuffled:

| run | per perturbation | per drug | per-drug p |
|---|---|---|---|
| `mc` (report 01) | 0.305 / 0.317 | 0.409 / 0.418 | 0.007 |
| `mc_noGlobals` (report 01) | 0.389 / 0.293 | 0.450 / 0.426 | 0.51 |
| `mc_matched_auprc` | 0.321 / 0.330 | 0.400 / 0.414 | 0.042 |
| **`mc_noGlobals_lr1e-3`** | **0.262 / 0.274** | 0.401 / 0.410 | 9.4e-4 |
| `mc_matched_noGlobals_lr1e-3` | 0.342 / 0.373 | 0.412 / 0.436 | 0.038 |

`mc_noGlobals_lr1e-3` has the best `rank_prob` of any run so far (0.262 per perturbation, against 0.305 for the
original `mc`), and 23% of perturbations put their own drug in the top 10% (`mc`: 14%). So removing the dose/time
inputs costs nothing once training works — report 01's apparent cost (0.389) was the failed folds.

**Training now runs.** Epochs kept per ensemble member:

| run | fold 0 | fold 1 | fold 2 | fold 3 | fold 4 |
|---|---|---|---|---|---|
| `mc_noGlobals` (lr 0.005) | 15, 7, 22 | 3, 5, 2 | 2, 1, 16 | 1, 1, 1 | 1, 8, 3 |
| `mc_noGlobals_lr1e-3` | 53, 21, 11 | 47, 1, 41 | 21, 22, 26 | 21, 2, 21 | 21, 40, 43 |

No fold stops at epoch 1–3 any more (two single members still do). One member kept epoch 53 of 60, so the cap may
still bind; `--ffn-epochs 100` would show whether anything is left.

## 4. Matched negatives are no longer needed

With dose/time inputs gone and training stable, the matched sampling **costs** performance instead of adding any:

| metric (per perturbation, real / shuffled) | `noGlobals_lr1e-3` | `matched_noGlobals_lr1e-3` |
|---|---|---|
| `spec_prob_timedose` | **0.414** / 0.500 | 0.441 / 0.504 |
| `rank_prob` | **0.262** / 0.274 | 0.342 / 0.373 |

Both are better than their controls in 5/5 folds, so this is not a failure of the matched run — it is simply worse
on both axes. Interpretation: matched negatives were a way of forcing the model off the per-protein prior, and
removing the dose/time shortcut plus training to convergence achieves the same thing without sacrificing the prior
(which `rank_prob` needs). Note the two runs also differ in epoch selection (`spec` vs AUPRC), so the comparison is
not clean; but since the matched run is worse on *both* metrics, no choice of selection metric rescues it.

`mc_matched_auprc` answers report 01's other open question: with dose/time inputs still present and selection on
AUPRC, the matched sampling gives unstratified `spec` 0.419 / 0.478 (per drug p = 8.8e-4) but stratified only
0.450 / 0.497 (per drug p = 0.45). So most of what looked like matched-negative specificity in report 01 was the
dose/time shortcut, and the rest came from selecting epochs on `spec`.

## 5. What this means

1. **The model does read the signature.** On the leak-free metric, with drugs as the unit, `mc_noGlobals_lr1e-3` is
   5.6 points better than its own shuffled control (difference of medians; paired median difference 4.2 points, p = 4.5e-5), in every fold. This is the first result in the
   project where per-drug evidence is strong on a metric that cannot be explained by the per-protein prior or by the
   dose/time design.
2. **The recipe is: drop the dose/time inputs, train to convergence.** Not matched negatives, not `spec`-based
   epoch selection. Both of those were fixes for symptoms.
3. **The effect is still modest in absolute terms.** 0.444 vs 0.500 per drug means the signature moves a drug's
   own perturbations from chance to slightly better than chance. Useful, not yet a discriminating tool.
4. **Per-drug plain ranking is still mostly prior.** `rank_prob` per drug 0.401 vs 0.410 control: the large
   per-perturbation gain (0.262 vs 0.274) does not survive averaging within drug. A query will therefore still
   return broadly similar drugs regardless of what is asked; the signature re-orders within that.

## 6. Caveats

- **Single seed per run.** Fold consistency (5/5) is reassuring but is not a replicate; the same splits are reused.
- **Two changes at once in the lr runs** (lr 0.001 **and** patience 10 **and** cap 60). They are all "train longer",
  so separating them is low value, but it is not a single-variable comparison.
- **The epoch cap may still bind** (one member at 53/60).
- **`matched_noGlobals_lr1e-3` differs from `noGlobals_lr1e-3` in two ways** (matched sampling and `spec` selection).
- **Report 01's conclusions about matched negatives are superseded** for the no-dose/time setting; they remain
  correct as statements about runs that kept the dose/time inputs.
- The training script's own validation `spec` is still unstratified, which no longer matters for `--no-globals` runs
  (no dose/time inputs) but still does for any run that keeps them.

## 7. Next steps

1. `--ffn-epochs 100` on `mc_noGlobals_lr1e-3` (real + shuffled) to check the cap is not limiting.
2. A second seed of `mc_noGlobals_lr1e-3` + control, to size run-to-run variability before building on it.
3. Freeze this as the reference recipe (`--no-globals --lr 0.001 --patience 10 --ffn-epochs 60`) and port it to
   `tf` and `expr`, each with its shuffled control — the first interpretable comparison across signature types.
4. Prototype path (see `STATUS.md` section 6): query → ranked perturbations, the signature-similarity baseline, and
   the Schuhknecht A549 benchmark. These do not need the GPU and judge the model on the actual task.
5. Version (b) (grouped specificity loss) is **demoted**: the problem it was designed to solve appears to be
   addressable by training properly. Revisit only if the per-drug effect plateaus.
