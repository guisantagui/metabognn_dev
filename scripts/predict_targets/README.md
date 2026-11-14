# Instructions for launching training scripts

All commands are run from this folder, on the cluster login node. Each run writes to
`/home/users/gsantamaria/projects/metabognn_dev/results/predict_targets/A549/<out_dir>`. Pull the repository on the
cluster before submitting, so that the jobs (and the comparison jobs) use the current scripts.

## Training scripts with shuffle controls

`submit_with_shuffle_control.sh` submits, for each variant, a training run on the real signatures
(`launch_predict_targets_A549_<variant>.sh`) and its control with the signatures shuffled among perturbations
(`launch_predict_targets_A549_<variant>_shuff.sh`, `--shuffle-signatures`) in parallel. It then submits a
`compare_runs.py` job (`launch_compare_runs.sh`) that starts when both have ended. The comparison starts even if one of
the runs failed, and then compares the folds finished in both. It is written to
`<real out_dir>/compare_real_vs_shuffled/` and `compare_runs_<jobid>.out`.

The model learns a target prior per protein, so a run must be compared with its shuffled control, not with chance.

```bash
bash submit_with_shuffle_control.sh tf                                # real + shuffled + comparison (TF data)
bash submit_with_shuffle_control.sh mc_noGlobals mc_matched_auprc     # several variants in one call
bash submit_with_shuffle_control.sh all                               # tf, expr and mc
bash submit_with_shuffle_control.sh expr --real-job <id>              # shuffled + comparison only: the real run is
                                                                      # already queued/running as <id> (0: finished)
```

Every variant name is checked before anything is submitted, so a typo submits nothing. Use `--real-job` when the
real run is already queued, running or finished: resubmitting its launcher would overwrite its out_dir. The script needs
bash ≥ 4 (`declare -A`), which the cluster has and macOS does not.

### Variants

| variant | what changes (vs `mc` unless stated) | out_dir (real; control adds `_shuff`) | wrapper |
|---|---|---|---|
| `tf` | TF-activity signatures | `tf_active200` | |
| `expr` | expression signatures | `expr` | |
| `mc` | MetChange signatures, random negatives, epochs selected on validation AUPRC | `metchange_rank` | |
| `mc_noStatic` | `--no-static`: no per-protein features (baseline activity and its propagation, degree, strength) | `metchange_rank_noStatic` | not in the submit script (launch both by hand) |
| `mc_noGlobals` | `--no-globals`: no per-perturbation features (relative dose, time, signature size) | `metchange_rank_noGlobals` | `submit_train_w_shuff_mc_noGlobals_matched_auprc.sh` |
| `mc_matched` | `--matched-negatives 4 --matched-weight 0.5 --select-metric spec` | `metchange_rank_matched` | `submit_train_w_shuff_mc_matched.sh` |
| `mc_matched_auprc` | as `mc_matched` but `--select-metric auprc` | `metchange_rank_matched_auprc` | `submit_train_w_shuff_mc_noGlobals_matched_auprc.sh` |
| `mc_matched_noGlobals` | as `mc_matched` plus `--no-globals` | `metchange_rank_matched_noGlobals` | `submit_train_w_shuff_mc_matched_noGlobals.sh` |
| `mc_noGlobals_lr1e-3` | as `mc_noGlobals` plus `--lr 0.001 --patience 10 --ffn-epochs 60`; array job | `metchange_rank_noGlobals_lr1e-3` | `submit_train_w_shuff_mc_noGlobals_lr1e-3.sh` |
| `mc_matched_noGlobals_lr1e-3` | as `mc_matched_noGlobals` plus the same training settings; array job | `metchange_rank_matched_noGlobals_lr1e-3` | `submit_train_w_shuff_mc_noGlobals_lr1e-3.sh` |
| `mc_noGlobals_lr1e-3_ep100` | as `mc_noGlobals_lr1e-3` but `--ffn-epochs 100`; array job, 6 h per task | `metchange_rank_noGlobals_lr1e-3_ep100` | `submit_train_w_shuff_mc_noGlobals_lr1e-3_ep100_seed1.sh` |
| `mc_noGlobals_lr1e-3_seed1` | as `mc_noGlobals_lr1e-3` plus `--seed 1` (new folds, init and permutation: a full replicate); array job | `metchange_rank_noGlobals_lr1e-3_seed1` | `submit_train_w_shuff_mc_noGlobals_lr1e-3_ep100_seed1.sh` |
| `tf_noGlobals_lr1e-3` | as `tf` plus the reference recipe (`--no-globals --lr 0.001 --patience 10 --ffn-epochs 60`); array job, 6 h per task | `tf_active200_noGlobals_lr1e-3` | `submit_train_w_shuff_tf_expr_noGlobals_lr1e-3.sh` |
| `expr_noGlobals_lr1e-3` | as `expr` plus the reference recipe; array job, 6 h per task | `expr_noGlobals_lr1e-3` | `submit_train_w_shuff_tf_expr_noGlobals_lr1e-3.sh` |
| `tf_noGlobals_lr1e-3_min15` | as `tf_noGlobals_lr1e-3` but `--min-epochs 15 --patience 15 --ffn-epochs 80`; array job, 8 h per task | `tf_active200_noGlobals_lr1e-3_min15` | `submit_train_w_shuff_noGlobals_lr1e-3_min15.sh` |
| `expr_noGlobals_lr1e-3_min15` | as `expr_noGlobals_lr1e-3` with the same training settings | `expr_noGlobals_lr1e-3_min15` | `submit_train_w_shuff_noGlobals_lr1e-3_min15.sh` |
| `mc_noGlobals_lr1e-3_min15` | as `mc_noGlobals_lr1e-3` with the same training settings | `metchange_rank_noGlobals_lr1e-3_min15` | `submit_train_w_shuff_noGlobals_lr1e-3_min15.sh` |

**`--min-epochs N`** (default 1 = plain early stopping, so older launchers behave exactly as before): epochs before N
are trained and logged but can neither be kept nor stop training; `--patience` counts from epoch N. Why: validation
AUPRC often peaks at epoch 1, dips, and recovers after the first learning-rate halvings (`StepLR` halves the lr every
10 epochs), so with patience 10 about half of the `tf` / `expr` FNNs kept epoch 1–2. Raising `--ffn-epochs` alone does
nothing for those (and little in general: after 60 epochs the lr is 1/64 of its start value).

Note: `--no-static` and `--no-globals` are different ablations. `noStatic` keeps dose, time and signature size.

**Array jobs.** The `_lr1e-3` launchers run one fold per task (`#SBATCH --array=0-4`, `--fold "$SLURM_ARRAY_TASK_ID"`, 4 h each instead of 12 h for all folds; 6 h for `_ep100` and the `tf` / `expr` ones, 8 h for `_min15`). Each task uses the same split and seed as the sequential run. Short GPU requests usually start sooner through backfill; check with `sbatch --test-only <launcher>`, which prints the estimated start time without submitting. Logs are `<job name>_<task>.out/.err`. The task that finishes last pools the folds (`cv_*`); if one task failed and was resubmitted alone (`sbatch --array=<i> <launcher>`), pool with `--aggregate-only --n-folds 5 --out-dir <out_dir>` (the comparison does not need it). A comparison submitted by the submit script waits for all tasks of both arrays.

### Why the `mc_noGlobals`, `mc_matched_auprc` and `mc_matched_noGlobals` pairs

Results so far are in `reports/2026-11-13_A549_matched_negatives.md`.

- **`mc_noGlobals`**: the unmatched model's shuffled control scores better than chance on `spec` (≈0.45 instead of 0.5),
  and a predictor that only knows dose and time does about as well. The unmatched model's apparent specificity
  therefore comes mostly from dose and time. Without these features the shuffled control should sit at 0.5 even
  without stratifying. Dose and time can still reach the model through the signature itself (its size and strength),
  so the stratified `spec` (below) is still the number to read.
- **`mc_matched_auprc`**: `mc_matched` changes two things at once compared with `mc`: the sampling (matched negatives)
  and the epoch selection (`spec` instead of AUPRC). This run has the same epoch selection as `mc`, so:
  - `mc_matched_auprc` vs `mc`: effect of the matched sampling alone;
  - `mc_matched` vs `mc_matched_auprc`: effect of selecting epochs on `spec`.

  Compare each run with its own shuffled control first. A direct comparison between two real runs can mislead,
  because their controls differ.
- **`mc_matched_noGlobals`**: `mc_noGlobals` showed that the dose/time shortcut goes through the dose/time inputs
  (without them the shuffled control's `spec` is 0.50 unstratified; report section 5.5). This pair measures the
  matched model's specificity without that shortcut, so its plain `spec` is leak-free. It also makes the
  validation `spec` used for epoch selection unable to reward dose / time through the inputs.
- **`mc_noGlobals_lr1e-3`, `mc_matched_noGlobals_lr1e-3`**: in `mc_noGlobals` training stopped after 1–3 epochs in several folds, and the signature effect appeared only in the folds that trained longer (report section 5.5). A lower learning rate with more patience and a higher epoch cap tests whether stable training gives a larger, more consistent effect.

## Matched negatives (`--matched-negatives`, `--matched-weight`, `--select-metric`)

(a) Plain (implemented): same loss, different sampling. Fold the matched rows in as ordinary negatives. Rows stay
independent; only the training distribution changes. What it buys you is decorrelation: without it a protein like
EGFR appears almost only as a positive, so protein identity predicts the label. Per chunk you have roughly

┌─────────────────────────────────────────────┬───────────┐
│                                             │   rows    │
├─────────────────────────────────────────────┼───────────┤
│ positives                                   │ ~100–300  │
├─────────────────────────────────────────────┼───────────┤
│ random negatives (--neg-per-pert 200, P=32) │ 6400      │
├─────────────────────────────────────────────┼───────────┤
│ matched negatives at K=4                    │ ~400–1200 │
└─────────────────────────────────────────────┴───────────┘

After adding them, P(target | protein = EGFR) in the training set drops from ~1 to ~1/5. That's a substantial change
for a small number of rows — but it also means they're easily swamped, so they need their own sample weight
(--matched-weight), not just to be dumped into class 0.

- `--matched-negatives K` (default 0 = off): for each target row (p, g), up to K rows (q, g) of the same protein g in
  perturbations q of other drugs in the chunk, with label 0. Excluded: other doses of the same drug, drugs that also
  target g, and unknown-direction labels. Only the signature can separate these rows from the positive. With K > 0,
  the train/validation perturbations are shuffled before chunking, so a chunk mixes drugs.
- `--matched-weight` (default 0.5, in [0, 1]): the share of the class-0 weight given to the matched rows.
- `--select-metric {auprc, macro_f1, spec, auprc_matched}`: the validation metric that picks the epoch kept.
  `spec` = mean over validation drugs of the spec percentile (lower is better; it ranks each drug among the other
  drugs for its own targets, so it is insensitive to the per-protein prior). `auprc_matched` = AUPRC on the matched
  validation rows only (needs `--matched-negatives`). `auprc` rewards the per-protein prior. Note: this `spec` is not
  stratified by dose/time, so it can still reward them.

`metrics.json` records `epoch_selection` (metric, best epoch, score, history per epoch) and the validation metrics.
The pooled CV summary adds `per_fold_spec_val`, `per_fold_auprc_matched` and `per_fold_best_epoch`.

(b) Grouped (not implemented): keep each group {(p,g), (q₁,g), …, (q_K,g)} together as a unit and put a loss on the
group: a softmax over the group (which of these K+1 perturbations targets g?) or a pairwise ranking loss on the positive
against each matched negative. Here the batch becomes (n_groups, K+1, F) instead of (n_rows, F). The group loss is
blind to per-protein offsets, so it would be added to the cross-entropy (not replace it):
L = CE + λ·L_spec, with CE on the original sampling. Check it against the shuffled control and `--no-globals`, because
dose and time do not cancel within a group.

## Checks after a run

- **`compare_runs.py`** (run automatically by the submit script; by hand:
  `python compare_runs.py --run real=<out_dir> --run shuffled=<out_dir_shuff> --control shuffled`): percentiles per
  drug and per perturbation, and paired tests of each run against `--control` per drug (drugs are the independent
  units; a drug's doses are not independent evidence). Besides the unstratified `spec`, it reports `spec` stratified by
  dose/time: `spec_<prob|logreg>_time_heldout` ranks each perturbation only among other drugs' perturbations at the
  same time, and `spec_<prob|logreg>_timedose_heldout` at the same time and relative dose. This removes the dose/time
  shortcut. These use mid-rank percentiles, which average exactly 0.5 under chance even with the ~5 competitors of a
  small stratum (the plain formula averages (1 + n/2)/(n + 1), i.e. 0.58 at n = 5).
- **`check_dose_time_leak.py`**: `python check_dose_time_leak.py <results/predict_targets/A549> <run> [<run> ...]`.
  Per fold: `spec` of the model, `spec` of a predictor that knows only (time, relative dose), and the share of score
  variance explained by (time, relative dose). If the dose/time-only predictor matches the model, the model's
  specificity is dose/time. Reads `fold*/test_scores_heldout.npz` and parses `<drug>_<dose>_<time>` perturbation names.
