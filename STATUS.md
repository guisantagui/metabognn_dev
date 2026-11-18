# Project state — 2026-11-18

Handoff note for continuing in a new session. Everything described here is committed up to
HEAD `d11744b`, except the files listed in section 8.

> **2026-11-18 — recipe state (report 06, not frozen).** Best complete run: **expr `_qsigtop` 10.5 points**
> (q < 0.05 then top 100 by |log2FC|, 1,358 drugs, one seed). `_qsigtop` matches `_min15` for all three signature
> types (±2 points, p ≥ 0.12) and beats `_qsig` by 5–8, so **requiring significance is free**; figure:
> `results/predict_targets/A549/performance_figure/` (`scripts/predict_targets/plot_performance_panels.py`); **SD scaling is
> dropped** (`_sdscale` 6.7, `_qsigtop_sdscale` ≤ `_qsigtop`). MetChange `_qsigtop` 10.3 is provisional (fold 0
> complete since 2026-11-18, 10.2 points) and covers only the 758 drugs with significant metabolites, so **`_min15`
> stays the primary MetChange run**. One seed for every `_qsigtop` run. See section 6.0.

> **2026-11-17 (evening) — selection settled, Schuhknecht null explained** (report 05, sections 7–8).
> `_qsigtop` (q < 0.05, then top 100 by |log2FC|) equals the unfiltered top 100 for expr (10.5 points),
> TF (7.8) and MetChange (10.3, folds 1–4; fold 0 missing): significance can be required at no cost.
> SD scaling does not help (`_qsigtop_sdscale` 8.1, `_sdscale` 7.1 provisional). The Schuhknecht
> benchmark is null in both seeds and both controls. Cause: on the limma labels MetChange does not agree
> with measured Schuhknecht metabolomics for the same drug (drug-specific ρ −0.013, below the
> expression-producers baseline +0.033, p = 7e-6; specificity permutation p = 1.0).
>
> Earlier today: corrected MetChange matches expr and TF (report 05, section 5): per-drug gap 8.3 / 7.7
> (seeds 0 / 1), vs expr 7.8 / 10.4 and TF 9.0 / 9.3, not significantly different.

**Read first:** `reports/2026-11-18_recipe_state.md` (report 06 — where the recipe stands, what is settled and
what is pending), then `reports/2026-11-17_signature_selection_scale.md` (report 05 — q-values unused; MetChange's
top-100 is chosen by metabolite scale; `_qsig` / `_sdscale`; section 5: corrected MetChange results), then `reports/2026-11-16_labels_network_coverage.md` (report 04 — the label bug in the
MetChange runs, expr on `combined`, full coverage check), then `reports/2026-11-14_seed_tf_expr.md` (report 03 — second seed, MetChange vs expr
vs TF, early stopping and lr schedule; sections 9–10 have the current recipe, the two-seed
result and the metabolic-coverage check), then `reports/2026-11-14_lr_noglobals_matched.md` (report 02 —
the current recipe), then `reports/2026-11-13_A549_matched_negatives.md` (report 01 — metric
definitions, the dose/time leak, how to read a comparison). `scripts/predict_targets/README.md`
documents every launcher and how to submit runs with their controls.

## Goal and standing constraints

Adapt the Menacho et al. 2026 GNN drug-target framework to LINCS signatures. A549 first,
then the 9 core cell lines (A375, A549, HA1E, HCC515, HEPG2, HT29, MCF7, PC3, VCAP).

- Do not commit; the user commits.
- Keep gene symbols in the prediction pipeline.
- KEGG-derived edges must not be redistributed.
- Uncharacterized drugs are not negatives (no annotation is not evidence of no target).
- **Priority is MetChange**; expr is a welcome byproduct. TF reproduces the published Menacho
  setting and is the **positive control** (the level to reach), not a result.
- Compare every run against **its own shuffled-signature control**, never against chance, and
  test **per drug**, not per perturbation (doses of one drug are not independent).
- Long write-ups go in `reports/`, not into chat (see section 9).

## 1. Data

LINCS L1000, A549, from the Ma'ayan lab **CycleGAN** release: transcriptome-wide, but most
genes are **inferred** from the 978 landmarks. Consequence to keep in mind: the inferred genes
give the coverage MetChange needs over the GEM's GPR genes, but add no independent
measurement — the information ceiling is set by the 978 landmarks, and imputation smooths
toward typical profiles. Two cheap, unrun checks that would size that ceiling (section 6.3).

## 2. Labels

**limma reference-level sign bug — fixed.** `model.matrix(~ perturbation)` took whichever of
the perturbation name and `DMSO_nan_<time>` sorted first as the reference, so for ~55% of
A549 perturbations limma returned `control - treated`. Correction: `log2FC -> -log2FC`,
`basal -> basal + log2FC` (`fix_limma_reference_sign.py`), with R's own level order
(`limma_reference_order.R`) deciding, not the data heuristic. Corrected expression labels are
in `expr/`; PC3 was already correct; MetChange has `rank_signfixed`. `limma_DE.R` now sets the
control reference explicitly and writes `reference_level.txt`, so new runs need no correction.

**But the home MetChange copy was never corrected** (report 04, section 1; checked 2026-11-16 with
`compute_labels/compare_label_dirs.py`). `<home>/results/labels/A549/metchange/rank/` vs the scratch
rerun: 19,884 perturbations identical, 24,512 exactly negated (r = −1, all features), and the
negated set is exactly `flipped_R = TRUE`; q / p-values and significant sets identical. The
correction was only ever copied into scratch (`move_corr_files_mc.sh`). Scratch is correct by
construction (fixed `limma_DE.R`); **confirmed 2026-11-16**: its `reference_level.txt` reads
"reference level: ctrl (explicit)" (written only by `limma_DE.R` ≥ `85526fe`, which includes the fix). The expr pair of that
comparison wrote no output (skipped or failed; `logs/cmp_labels.out` not yet looked at).

**MetChange limma rerun over the 9 cell lines — in progress, verify.** The jobs segfaulted in
`fread` (data.table crashes on the ~150k-column MetChange matrix). Fixed in
`launch_limma_DE_generic_mc.sh` (commit `de76198`): `mc` mode now passes `--reader lines`, and
output goes to `labels/<cell>/metchange/<penalty>/` instead of `labels/<cell>/metchange/`.
Resubmit with `launch_limma_mc_celllist.sh` and check the first `.out` for the
"N features read" progress lines.

**MetChange concatenation.** `mc_res_concat_list_launcher.sh` is a Slurm array, one task per
cell line (`--array=1-9`, `-c 8`, `POLARS_MAX_THREADS`, skips cell lines that already have a
`*_concat.csv` unless `overwrite=true`). Note: `concat_mc_res.py`'s docstring claims it
computes deltas against matched controls and averages per perturbation; the code only
concatenates batches and drops all-missing metabolites. Only the description is wrong.

## 3. Where the model stands

Two problems were found and are now understood:

1. **Per-protein prior.** The model learned which proteins tend to be drug targets, which does
   not depend on the signature. `spec` / `rank_probc` cancel it (report 01, section 3).
2. **Dose/time shortcut.** In LINCS the experimental design tracks the screen, and so the drug
   class. The FNN exploited it through the `rel_log_dose` / `log_time` inputs and their
   products with protein features, which inflated `spec` even in shuffled runs (report 01,
   sections 5.3–5.5). Fixes: `--no-globals`, and the dose/time-stratified `spec` in
   `compare_runs.py`.

**Current recipe** (report 03, section 9): `--no-globals --lr 0.001 --min-epochs 15
--patience 15 --ffn-epochs 80` (StepLR), unmatched, epochs selected on AUPRC. The `_min15` runs.

**Gap convention** (report 03, section 1): per-drug gap = `median_diff` of `paired.tsv`, the
median of the per-drug real − control differences, in percentile points; it is what the
Wilcoxon p tests. Report 02's "5.6 points" was a difference of medians (paired: 4.2).

Per-drug gap (points), p — real vs own control, current recipe, seeds 0 / 1. MetChange on the
corrected labels (report 05, section 5); expr / TF from report 03, 9–10:

| metric | MetChange | expr | TF (positive control) |
|---|---|---|---|
| `spec_prob_timedose` (leak-free) | 8.3 / 7.7 (6e-12 / 4e-11) | 7.8 / 10.4 (2e-11 / 1e-14) | 9.0 / 9.3 (8e-15 / 7e-15) |
| `rank_probc` | 9.3 / 7.7 | 8.0 / 11.5 | 10.4 / 10.3 |
| `rank_prob` | 2.5 / 2.5 (1e-11 / 5e-9) | 3.3 / 3.6 | 3.3 / 3.3 |

Paired by drug, MetChange vs expr: p = 0.79 / 0.085; vs TF: 0.21 / 0.12; expr vs TF: 0.36 / 0.68.
**All three signatures are indistinguishable.** (Sign-buggy MetChange was 4.2 / 2.7, 6.4 / 4.4,
1.3 / 0.7; corrected vs buggy, seed 1: p = 2.5e-5.)

- **The warm-up (`--min-epochs 15`) fixed the collapse**: no FNN keeps epoch 1–2 any more
  (before: 7–8/15 for TF and expr). Expr was training-limited and now matches the positive
  control, all five folds at 0.30–0.32 (`spec_prob_timedose` per perturbation); TF also gained
  (6.8 → 9.0).
- **Schedule variants** (patience 25, `ReduceLROnPlateau`) made no difference on sign-buggy
  MetChange, nor on expr / TF; keep StepLR. Not retested on corrected MetChange.
- **Resolved: "MetChange reads the signature half as well as expr"** (reports 03–04) was the
  label sign bug. Network ruled out independently (below).
- **Expr on the `combined` network performs as well as on `gene`** (report 04, section 2; seed 0,
  same folds, `expr_combined_noGlobals_lr1e-3_min15`). Per-drug gap over own control:
  `spec_prob_timedose` 7.7 (p = 3e-12) vs 7.8 on `gene`, `rank_probc` 8.2 vs 8.0, `rank_prob` 3.0
  vs 3.3; real better in 5/5 folds. The MetChange network is not why MetChange lags expr.
- **Robust across two seeds: MetChange ≈ expr ≈ TF** (corrected labels). Per-drug gaps agree
  only moderately between seeds (MetChange Spearman 0.32): read effects in aggregate, not drug by drug.
- **Metabolic coverage** (full data, report 04 section 3): in DE terms it is not the issue
  (every active perturbation has ~14–15% of its DE genes in Human-GEM, the background share).
  In MetChange's own output it is all-or-nothing: 57% of active perturbations have no
  significant feature among the 100 the model reads, 35% have all 100. Corrected MetChange
  (report 05, 5): 14.6 / 14.7 points on fully covered perturbations (expr on the same: 12.1 /
  10.5), 5.4 / 5.6 on uncovered ones (expr 7.3 / 11.8). MetChange beats expr where it has a
  response; uncovered perturbations are what pull its average down to the expr level.
- The epoch cap does not bind for MetChange; for TF / expr 8–12/15 members reach epoch 80, but
  at an lr ≤ 1.6e-5, where more epochs change little.

**Matched negatives are superseded.** With dose/time inputs gone and training to convergence,
matched sampling is worse on both metrics (report 02, section 4). What looked like
matched-negative specificity in report 01 was mostly the dose/time shortcut plus `spec`-based
epoch selection. Version (b) (grouped specificity loss) is demoted accordingly.

**Plain ranking:** with corrected labels the MetChange signature moves per-drug `rank_prob` by
2.5 points (0.381 vs 0.415 control), close to expr / TF (3.3–3.6). The prior still dominates the
absolute ranking for all three (real 0.38–0.40 per drug vs 0.5 chance).

## 4. Code in `scripts/predict_targets/`

`metabolic_target_prediction.py`:

- Prior-free metrics in `metrics.json` / `cv_metrics.json` / rank tables: `spec_prob_heldout`,
  `spec_logreg_heldout`, `rank_probc_heldout`, `rank_logregc_heldout`; `test_scores_heldout.npz`
  per fold (every test perturbation's score for every held-out drug).
- `--no-static` (drop per-protein features) / `--no-globals` (drop dose, time, signature size).
  **Different ablations** — `noStatic` keeps dose and time.
- `--score-active-only`, `--min-de-genes` with `--de-counts`.
- `--matched-negatives K`, `--matched-weight`, `--select-metric {auprc,macro_f1,spec,auprc_matched}`;
  `epoch_selection` and `per_fold_{spec_val,auprc_matched,best_epoch}` in the metrics.
- `--min-epochs N` (default 1 = plain early stopping, previous behaviour; verified by replaying
  the 150 recorded member histories): epochs < N are trained and logged but can neither be kept
  nor stop training; `--patience` counts from N. Recorded as `min_epochs` in `epoch_selection`.
- `--lr-schedule {step,plateau}` (default `step` = `StepLR`, `--lr-step-size 10 --lr-gamma 0.5`,
  as before); `plateau` = `ReduceLROnPlateau` on the validation selection score from epoch 1,
  `--plateau-patience` (3, must be < `--patience`), `--min-lr` (1e-6). Each epoch's lr is in
  `history`. Checked with a stub-torch test of `train_ffn` only; first real run is on the
  cluster (report 03, 8.2).
- `--fold i` with `--n-folds 5` for array jobs; the last task pools (`--aggregate-only` to pool
  by hand).

Tooling:

- `compare_runs.py` — summary / per fold / paired per-drug Wilcoxon / ECDFs, plus `spec`
  **stratified** by time and by time × dose, with mid-rank percentiles (unbiased for small
  strata; the plain formula averages 0.58 at n = 5 — report 01, section 5.4).
- `check_metabolic_coverage.py` — per-perturbation metabolic coverage (expression DE in
  Human-GEM; MetChange's own significant features and significant features among the top-n
  read) and the run's spec split by it, vs its control (report 03, 10.2–10.3).
- `check_dose_time_leak.py` — `spec` of the model vs `spec` of a dose/time-only predictor, and
  the score variance explained by dose/time groups.
- **Labels (2026-11-15):** every training launcher reads its labels from scratch,
  `${project_dir_scratch}/results/labels/<cell>/{expr,metchange/rank}` (MetChange used the home
  copy until then). The scratch copy is the rerun of the fixed `limma_DE.R`; the home copy was
  ~~corrected with `fix_limma_reference_sign.py`~~ — **wrong for MetChange**: the home
  `metchange/rank` is the uncorrected output (report 04; the earlier "spot check confirms it" was
  not true). The expr launchers always read scratch.
  Label pipeline, all writing to scratch: `parse_lincs_data/launch_get_lincs_celllist.sh`
  (expression limma) → `metchange/metchange_list_launcher.sh` → `metchange/mc_res_concat_list_launcher.sh`
  → `compute_labels/launch_limma_mc_celllist.sh`.
- **Output location (2026-11-15):** every launcher writes to
  `/scratch/users/gsantamaria/test_large_files/metabognn_dev/results/predict_targets/A549/<out_dir>`
  (`${project_dir_scratch}`), and `submit_with_shuffle_control.sh` /
  `compare_runs_mc_vs_shuffled_A549.sh` read from there. Earlier runs are in the same folders
  under `/home/users/gsantamaria/projects/metabognn_dev/results/predict_targets/A549/`.
- `submit_with_shuffle_control.sh` — submits a variant's real run, its shuffled control and a
  dependent comparison; takes several variants in one call, validates names first, `--real-job`
  for an already-queued real run. Needs bash ≥ 4 (cluster yes, macOS no).
- Wrappers: `submit_train_w_shuff_mc_matched.sh`,
  `submit_train_w_shuff_mc_noGlobals_matched_auprc.sh`,
  `submit_train_w_shuff_mc_matched_noGlobals.sh`,
  `submit_train_w_shuff_mc_noGlobals_lr1e-3.sh`,
  `submit_train_w_shuff_mc_noGlobals_lr1e-3_ep100_seed1.sh` (done, report 03),
  `submit_train_w_shuff_tf_expr_noGlobals_lr1e-3.sh` (done, report 03),
  `submit_train_w_shuff_noGlobals_lr1e-3_min15.sh`,
  `submit_train_w_shuff_mc_noGlobals_lr1e-3_pat25_plateau.sh`,
  `submit_train_w_shuff_tf_expr_noGlobals_lr1e-3_plateau.sh` (done, report 03),
  `submit_train_w_shuff_noGlobals_lr1e-3_min15_seed1.sh` (done, report 03 section 10),
  `submit_train_w_shuff_expr_combined_noGlobals_lr1e-3_min15.sh` (6.1).
- Variants registered in the submit script: `tf`, `expr`, `mc`, `mc_matched`,
  `mc_matched_auprc`, `mc_matched_noGlobals`, `mc_noGlobals`, `mc_noGlobals_lr1e-3`,
  `mc_matched_noGlobals_lr1e-3`, `mc_noGlobals_lr1e-3_ep100`, `mc_noGlobals_lr1e-3_seed1`,
  `tf_noGlobals_lr1e-3`, `expr_noGlobals_lr1e-3`, `mc_noGlobals_lr1e-3_min15`,
  `tf_noGlobals_lr1e-3_min15`, `expr_noGlobals_lr1e-3_min15`, `mc_noGlobals_lr1e-3_min15_pat25`,
  `{mc,tf,expr}_noGlobals_lr1e-3_plateau`, `{mc,tf,expr}_noGlobals_lr1e-3_min15_seed1`,
  `expr_combined_noGlobals_lr1e-3_min15`, `{mc,tf,expr}_noGlobals_lr1e-3_min15_qsig`.
  `mc_noStatic` is not registered
  (launch both by hand).
- The `_lr1e-3` launchers are **array jobs** (`--array=0-4`; 4 h per task, 6 h for `_ep100` and
  the `tf` / `expr` ones, 8 h for `_min15`, `_min15_pat25`, `_plateau`). `sbatch --test-only` prints the estimated start time without
  submitting.

## 5. Runs in `results/predict_targets/A549/`

Analysed (in reports): `metchange_rank{,_shuff}`, `metchange_rank_matched{,_shuff}`,
`metchange_rank_noStatic{,_shuff}`, `metchange_rank_noGlobals{,_shuff}`,
`metchange_rank_matched_auprc`, `metchange_rank_noGlobals_lr1e-3{,_shuff}`,
`metchange_rank_matched_noGlobals_lr1e-3{,_shuff}`, and (report 03)
`metchange_rank_noGlobals_lr1e-3_{ep100,seed1}{,_shuff}`, `tf_active200_noGlobals_lr1e-3{,_shuff}`,
`expr_noGlobals_lr1e-3{,_shuff}`. Comparisons live in each real run's
`compare_real_vs_shuffled/`; cross-run comparisons in `mc_matched_vs_unmatched/`,
`mc_real_vs_suff/`, `stratified_spec/`.

Not on the laptop: `metchange_rank_matched_auprc_shuff` (its comparison ran on the cluster),
and `metchange_rank_matched_noGlobals{,_shuff}` (plain, lr 0.005) — launched 2026-11-13, not
yet synced or analysed. It is now of low value: report 02 supersedes it.

Analysed in report 03, section 9 (seed 0, same splits as `*_noGlobals_lr1e-3`):
`metchange_rank_noGlobals_lr1e-3_{min15,min15_pat25,plateau}{,_shuff}`,
`expr_noGlobals_lr1e-3_{min15,plateau}{,_shuff}`,
`tf_active200_noGlobals_lr1e-3_{min15,plateau}{,_shuff}`. Report 03, section 10 (seed 1):
`{metchange_rank,expr,tf_active200}_noGlobals_lr1e-3_min15_seed1{,_shuff}`. Report 04:
`expr_combined_noGlobals_lr1e-3_min15{,_shuff}` (seed 0), `metabolic_coverage_full{,_topnbins}/`.
All MetChange runs listed here used the sign-buggy home labels (report 04).
Corrected (report 05, section 5): `metchange_rank_noGlobals_lr1e-3_min15{,_seed1}{,_shuff}` (the laptop
copies are now the corrected runs; the seed-0 buggy pair exists only on scratch as `_signbug`),
`metabolic_coverage_fixedlabels/`. Report 05, 5.2: `metchange_rank_noGlobals_lr1e-3_min15_qsig{,_shuff}`. On scratch, `metchange_rank_noGlobals_lr1e-3_min15{,_seed1}{,_shuff}` were renamed to
`..._signbug{,_shuff}` (2026-11-16), and `..._min15_pat25{,_shuff}` were deleted 2026-11-17 (laptop
copies remain).

Pre-dating the current metrics, cannot be re-checked (no `test_scores_heldout.npz`, no FNN
shuffle control): `expr`, `tf`, `tf_active200`, the `*_linOnly_shuff` runs.

## 6. Next steps, in priority order

**6.1 Corrected MetChange: done** (report 05, section 5: MetChange ≈ expr ≈ TF, both seeds).
Leftovers:
1. Read `compute_labels/logs/cmp_labels.out` for why the expr pair of the label comparison wrote
   nothing (fine if skipped as same provenance; the expr runs read scratch).
2. Check other consumers of the home MetChange copy (similarity checks, `metpath`,
   `mc_schuhk_assoc`) and the other cell lines' MetChange labels (`compare_label_dirs.py`).
3. The recipe parts established only on buggy MetChange (lr 0.001 vs 0.005, matched negatives
   "worse") are low priority now that MetChange matches the positive control.

**6.1b Signature selection** (no MetChange-vs-expr gap left to explain; these may still improve it).
- Network: **done, ruled out** (report 04, section 2: expr on `combined` 7.7 vs 7.8 points on `gene`).
- Coverage: **done, real** (reports 04–05): 57% of active perturbations have no significant
  MetChange feature among the 100 read. Corrected MetChange: 14.6–14.7 points on fully covered
  perturbations (beats expr there), 5.4–5.6 on uncovered. `_qsig` is in effect the MetChange
  activity filter (keeps 1,791 of 3,970 perturbations, 758 drugs).
- **MetChange `_qsig` done: worse** (report 05, 5.2): 6.1 points vs 11.9–14.2 for `_min15` on the same
  perturbations (p ≤ 3e-3). Half the training data dropped + random 100 of thousands of significant
  features. Magnitude-based top-100 stays; `_sdscale` is the test of scale.
- **Expr `_qsig` done: signal gone** (report 05, 5.3): 1.6 points (p = 0.22) vs 7.8 / 10.4 for `_min15`,
  with no perturbation dropped (all have ≥ 200 significant genes). Random 100 of the significant
  features loses the signal; the largest changes carry it. TF `_qsig` still running.
- **Significance (launchers ready, 2026-11-17):** every run so far, MetChange, expr and TF, used
  `--ignore-qvalue --top-n 100` (top 100 by |log2FC|, q-values unused; added 9 Nov in `5a28a0f`
  with the activity filter, no recorded reason). The published setting is significant features,
  capped at 100. `bash submit_train_w_shuff_noGlobals_lr1e-3_min15_qsig.sh` runs `_min15_qsig`
  (`--qvalue 0.05 --max-sources 100`) for all three, with controls. Perturbations with < 3
  significant features drop out, so compare signature types on shared drugs / perturbations.
  `_qsig` samples the significant features at random and lost the signal (report 05, 5.2–5.3).
  **`_qsigtop` (ready, 2026-11-17):** q < 0.05, then the top 100 of those by |log2FC|
  (`--qvalue-top-n`, new). Reason: of the top 100 by |log2FC|, a median of 1 is significant for MetChange
  and 84 for expr (report 05, 5.4). `bash submit_train_w_shuff_noGlobals_lr1e-3_min15_qsigtop.sh` (also
  MetChange `_qsigtop_sdscale`: plus `--signature-scale sd`).
  **`_qsigtop` done (report 05, section 7): no cost.** expr 10.5, TF 7.8, MetChange 10.3 points (folds 1–4;
  **fold 0 of the MetChange real run missing — rerun it**;
  `..._min15_qsigtop/` has fold1–4 only and no `cv_metrics.json`, while its `_shuff` control and the expr / TF
  `_qsigtop` runs have all 5, then re-run `compare_runs.py` for it); paired vs `_min15` within ±2 points (p ≥ 0.12), +5–8 over `_qsig`.
  `_qsigtop_sdscale` 8.1, not better than `_qsigtop` (−2.0, p = 0.12). MetChange `_qsigtop` covers only
  the 758 drugs with significant metabolites.
- **Feature scale (report 05):** MetChange log2FC is not standardised per metabolite (the MetChange
  paper z-scores against controls), so the current top 100 is mostly a fixed set of high-variance
  metabolites (selection frequency vs metabolite SD: Spearman 0.89; 33 metabolites in > 50% of
  signatures; different drugs share 22% of features vs 2–3% ranked by p). Expr is affected much
  less (5%). Fixes: `_qsig` (scale-free via significance) and **`--signature-scale sd`**
  (implemented: per-feature SD across perturbations, not centred, floored; z-scoring against
  controls rejected, many controls are constant in many features).
  `bash submit_train_w_shuff_mc_noGlobals_lr1e-3_min15_sdscale.sh` (MetChange, seed 0, with control).
  **`_sdscale` done** (fold 3 rerun 2026-11-17, compare re-run on all 5 folds, 1,358 drugs): 6.7 points, vs
  `_min15` −1.6 / +0.4 (p = 0.46 / 0.66), vs `_qsigtop_sdscale` −1.3 (p = 0.62). **SD scaling is neutral: drop
  it from the recipe.**
- Representation (code + GPU): continuous MetChange values or a larger `--top-n` instead of the
  top-100 signs (`--signature-weights`, section 7).

**6.2 Prototype path** (no GPU; this is what "working prototype" needs):
1. **Query path** — nothing yet turns a user's metabolic query into a ranked drug list. With
   `--no-globals` the model can score a query directly (a query has no dose or time), and its
   predicted targets map to drugs that hit them.
2. **Signature-similarity baseline** — rank perturbations by similarity of their MetChange
   signature to the query, no network. If this matches the model, the network is not yet
   earning its place. Needed for the paper either way.
3. **Schuhknecht A549 benchmark** — query with the measured metabolic change of a drug and
   check whether that drug (or its mechanism class) ranks high. This is the only validation on
   the actual task; `spec` only measures whether the signature is read at all. Partly in
   `scripts/mc_schuhk_assoc`.
   **Compartments (to handle when running it):** training signatures are compartment-specific
   (`MAM00003c`; 3,747 of 4,086 MetChange features are exact network nodes, weight ±1 each, the
   ±1/n split in `select_sources` never triggers). Intracellular metabolomics gives one value per
   metabolite, not per compartment, so a Schuhknecht feature must be spread ±1/n over its
   **intracellular** copies only (all compartments except `e`). Today `feature_node_map` maps a
   compartment-free id to *every* copy, including `e`: needs an option to exclude
   extracellular. Train/query mismatch: training sources are single nodes, query sources are
   spread; if the benchmark is weak, consider also training on compartment-collapsed signatures.
   Inputs now in `data/validation/schuhknecht/` (maps from `scripts/schuhknecht_map/`); counts and
   design in report 05, section 6: 473 non-ambiguous intracellular ions → 959 nodes (412 ions /
   801 nodes are MetChange training sources, cluster labels); 646 of 750 matched drugs have targets, 293 were
   CV-evaluated in A549, ~353 never entered training. Use `cpd_info_wSchuhknecht.csv` (all
   compartments), not `metabolite_match_table_non_ambig_match.csv` (`c` only).
   **Script ready (2026-11-17):** `benchmark_schuhknecht.py`, `sbatch launch_benchmark_schuhknecht_A549.sh`
   (report 05, 6.1). Must run on the cluster: locally the rebuilt features do not reproduce the
   training scores (local network / baseline / targets differ from the cluster copies), and the
   script refuses to score in that case. Copy `data/validation/schuhknecht/` to the cluster first.
   Cluster reproduction check passed (r = 1.000000, 8/8). Local MetChange labels are an older
   version (4,086 features vs 5,707 on the cluster): local counts based on them are approximate.
   **Result: chance in both seeds, real and control** (report 05, 6.2 and 8.1).
   Next: does MetChange agree with Schuhknecht at all, on the limma labels:
   `scripts/schuhknecht_mc_assoc/` (`sbatch launch_schuhknecht_mc_assoc_A549.sh`: prep, overlap
   enrichment with MetChange q-values, the R benchmark vs the expression baseline; R ran on the laptop,
   no R on the cluster). **Done (report 05, 8.2): no agreement.** MetChange drug-specific ρ −0.013 vs
   expression producers +0.033 (p = 7e-6); MetChange significant for only 85 of 750 matched drugs, and
   then for most metabolites; overlap enrichment equal to the size-matched expression baseline; sign
   agreement 57% vs 54% mismatched on Schuhknecht-significant metabolites. The benchmark cannot work
   with MetChange-trained models; options: train / calibrate on measured metabolomics, or drop the
   Schuhknecht query as validation of the MetChange path.

**6.0 Recipe as it stands (report 06, 2026-11-18). Not frozen** — kept here to take over from.

| signature | `_min15` (top 100 by \|log2FC\|) | `_qsigtop` (q < 0.05, then top 100) | `_qsig` | `_sdscale` |
|---|---|---|---|---|
| MetChange | 8.3 / 7.7 (seeds 0 / 1; 1,358 drugs) | **10.2** (758 drugs, complete) | 6.1 (758) | 6.7 (1,358) |
| expr | 7.8 / 10.4 | **10.5** (1,358) | 1.6 | — |
| TF | 9.0 / 9.3 | 7.8 (1,346) | 4.2 | — |

(per-drug gap over own shuffled control, `spec_prob_timedose`, points; MetChange `_qsigtop_sdscale` 8.1 / 758.)

Flags, as in the `_qsigtop` launchers (`launch_predict_targets_A549_{mc,tf,expr}_noGlobals_lr1e-3_min15_qsigtop{,_shuff}.sh`,
wrapper `submit_train_w_shuff_noGlobals_lr1e-3_min15_qsigtop.sh`):
`--signature-sign -1 --mode metchange --qvalue 0.05 --qvalue-top-n --top-n 100 --de-counts <expr>/fit_info_nsig.csv
--min-de-genes 200 --score-active-only --lr 0.001 --patience 15 --ffn-epochs 80 --min-epochs 15 --no-vae
--no-globals --n-folds 5`, **no `--signature-scale`**; `networks/gene` for expr / TF, `combined` for MetChange,
corrected scratch labels.

**Settled:** (1) requiring significance is free — `_qsigtop` ≈ `_min15` (±2 points, p ≥ 0.12) and ≫ `_qsig`
(+5–8, p ≤ 5e-5), and every feature read is then significant; expr and TF switch at no cost. (2) **SD scaling
out** (6.7 on full coverage, p ≥ 0.46 vs `_min15`; no gain on top of significance). (3) **Threshold-then-sample
out** (`_qsig`: expr 1.6 points, p = 0.22) — significance must rank, not sample. (4) **MetChange trade-off:**
`_qsigtop` ties `_min15` on its 758 covered drugs but cannot score the other 600, so **`_min15` remains the
primary MetChange run** and `_qsigtop` the significance-backed variant.

**Performance figure (report 06, section 5):** `scripts/predict_targets/plot_performance_panels.py` →
ECDF real vs own control per signature, forest of per-drug gaps with bootstrap CIs, top-k recovery — all on
`_qsigtop` (no other recipe to explain; `--extra sig=run:label` adds comparison rows to the forest when
wanted). Presentation format: 11 x 5.95 in (2:1), runs labelled only MetChange / expression / TF, no panel letters,
transparent background (`--background surface` for opaque). Light ink by default; `--theme dark` writes a
`..._dark` variant with the dark steps of the same hues for a dark slide. Output
`results/predict_targets/A549/performance_figure/performance_panels_<metric>.{png,pdf,tsv}`. The coverage-split
panel was dropped: it only exists for `_min15` (in `_qsigtop` every feature read is significant), so it is the
motivation for the recipe, not a result of it.
**Metric:** headline `rank_probc_heldout` (the deployment task, prior-controlled; MetChange / expr / TF
`_qsigtop` = 10.0 / 8.5 / 8.8 points), `spec_prob_timedose_heldout` beside it as the dose/time-leakage-controlled
check (10.2 / 10.5 / 7.8; the metric of reports 01–05), `rank_prob_heldout` in the supplement (5.8 / 3.4 / 4.0).
**Watch out:** `rank_*_heldout` in `cv_drug_ranks.tsv` are **raw ranks** — divide by `n_candidates_heldout` for a
percentile (only the stratified `spec_*` are percentiles already). Report 05 section 7's `rank_probc` numbers are
in rank units ×100, so their magnitudes are not comparable with the `spec` points; signs and p-values hold.
Re-run when a second `_qsigtop` seed exists.

**Pending / caveats:**
- **MetChange `_qsigtop` fold 0: done** (landed 2026-11-18 00:40, `compare_runs.py` re-run). 10.2 points,
  758 drugs, p < 1e-4; paired vs `_min15` s0 / s1 +0.3 / +1.1 (p = 0.76 / 0.82), vs `_qsig` +6.0 (p = 0.0023),
  vs `_qsigtop_sdscale` +1.8 (p = 0.12). All runs in the table are now complete (5 folds).
- **One seed only** for all `_qsigtop` runs; needed before calling it the recipe (the expr seed spread,
  7.8 vs 10.4, is larger than most differences being compared).
- `_sdscale` is complete (6.7) and can be dropped from future comparisons.
- Open, not selection: representation (continuous MetChange values or a larger `--top-n` instead of 100 signs,
  `--signature-weights`), and validation (6.2b).

**6.0b How to improve performance (2026-11-18).** The current recipe beats its shuffled control by 8-10 points
(median percentile ~0.36-0.40 vs ~0.47-0.50; true drug in the top 10% for ~12% of drugs vs ~5% control): real,
reproducible, not strong. **Selection is exhausted** — `_qsigtop` ≈ `_min15` ≈ `_sdscale` and only `_qsig` is
worse, so *which* 100 features are read barely matters; do not spend more there. In order:

1. **Diagnostic — a predictive similarity baseline scored on the model's metric (no GPU).**
   *Not* the model-free check, which is **already done**: `check_signature_target_similarity.py` →
   `results/signature_similarity/A549/` (expr) and `.../A549/metchange/` ask whether drugs sharing a target
   have similar signatures, and the answer is clearly yes:

   | | expr (3,970 perts, 20,210 genes) | MetChange (2,022 perts, 5,707 features) |
   |---|---|---|
   | same drug, other dose/time (reproducibility ceiling) | AUROC 0.848 | 0.736 |
   | different drugs sharing a target (same time point) | 0.672 (0.693) | 0.561 (0.568) |
   | on drug consensus signatures | 0.590 | 0.529 |
   | 10-NN: neighbours of other drugs sharing a target | 26.6% vs 2.7% expected (**9.9x**) | 15.1% vs 2.8% (**5.5x**) |

   (permuted-label nulls all ~0.50 / 1.0x.) So the signatures do carry target information, and expr carries more
   of it than MetChange — consistent with the model results.
   **The predictive version is now written (2026-11-18), not yet run:** `scripts/predict_targets/
   baseline_signature_knn.py` + `sbatch launch_baseline_signature_knn_A549.sh` (cluster only: each fold is checked
   against the run's held-out drugs and test perturbations, and the laptop inputs differ — older labels, no `gene`
   network). For each held-out perturbation, the k = 25 most similar **training** drugs vote their targets, and the
   protein scores go through the model's own `DrugRanker` / `signature_specificity` / `centered_ranks` on the
   run's own folds, queries and candidates, against the run's `_shuff` (shuffle replayed). Two inputs:
   - `sources` = exactly the model's input (selected features, signs, on nodes; cosine) → **network vs no
     network on the same input**;
   - `profile` = the full continuous signature (Spearman) → **is there more in the data than the model reads?**
     (directly informs item 3).
   One call per run scores the run and its `_shuff`, both representations and every k from a single read of the
   labels; the launcher is a 3-task array (one per signature type), so wall time is one expr read, ~30-50 min.
   Output is run-shaped, so `compare_runs.py` and the performance figure read it unchanged; the launcher draws
   `performance_figure/vs_baseline/` with the baseline rows under each model row. Smoke-tested end to end on the
   laptop with `--allow-mismatch` (code path only; numbers meaningless there).
   **Reading:** baseline ≈ model on `sources` → the network adds nothing over similarity of the same input;
   `profile` > `sources` → the 100-sign input is the bottleneck (go to item 3); model ≫ both → the network earns
   its place and the effort goes into data (item 4). Not a baseline: ranking a held-out drug by similarity to its
   *own* other doses (that is the same-drug ceiling, AUROC 0.85). Also not this: `rank_baseline_heldout` in the
   run outputs is the protein-prior baseline (`prior_scores`), signature-independent, identical in real and
   shuffled runs (0.508 / 0.508).
2. **Diagnostic — split the 10 points into prior vs signature.** The model learns a protein prior (validation
   AUPRC 0.348 real vs 0.280 shuffled, base rate 0.027; full proteome 0.0116 vs 0.0074). Retrain with
   `--no-static` (per-protein features removed, so signature only) and, the other way round, with the signature
   removed, and see how the gap splits. If signature-only is near chance, the weak link is signature → target
   and **more data will not fix it**; if it holds up, the prior is just adding a constant and the signature path
   is worth investing in.
3. **Representation: stop throwing the magnitudes away.** `metabolic_target_prediction.py:372` is
   `sigs[p] = (f, np.sign(v[keep].to_numpy()))` — the signature the model sees is **100 features, ±1 each**. It
   never sees how much anything changed, nor feature 101. The `--top-n 100` default arrived with the activity
   filter in `5a28a0f` with no recorded reason. Two small changes, the best effort-to-payoff ratio on this list
   (and the natural reading of why selection makes no difference: the bottleneck is that there are only 100
   signs):
   - **`--signature-weights`** (open in section 7, not implemented): continuous weights instead of signs —
     scaled log2FC, or signed −log10 p. Keep `--signature-scale sd` available here: scaling was neutral for
     *selection* but matters once the value itself is the input.
   - **`--top-n` sweep**: 100 → 300 → 1000 → all significant, each with its `_shuff` control.
4. **Biggest data lever: train across cell lines.** Everything so far is A549 alone (~1,358 drugs); LINCS has
   hundreds of lines and a GNN this size is probably data-starved. Pool 3-5 well-sampled lines with cell
   identity as a covariate: a drug's targets are the same everywhere, so the label is shared while the
   signature varies, which is exactly the contrast the model needs. Costs a label run per line
   (`scripts/compute_labels`), so it is the expensive item — but the one most likely to move the number
   substantially. Ties into 6.3 (data ceiling).
5. **Objective, if 2 points at the signature path:** the FFN trains with per-protein cross-entropy while we
   report a ranking metric; a pairwise / listwise ranking loss optimises what is actually measured.

**Already tried, do not repeat:** matched negatives (`--matched-negatives`, report 01) — once the dose/time
inputs were removed and training ran to convergence, they no longer help (report 02); signature scaling
(`_sdscale`, neutral, 6.0); threshold-then-sample selection (`_qsig`, worse).

**6.2b Validation strategy after the Schuhknecht null (report 05, sections 9–9.3).** Inferred labels do not
agree with measured metabolomics of the same drug (8.2). **Key distinction (9.3): pool metabolomics does not
validate the MetChange labels at all.** A MetChange label claims production / consumption *capacity*; the ground
truth for that is flux. An enzyme inhibitor acts post-translationally, so a correct MetChange need not mirror the
substrate-up / product-down pattern in pools — disagreement is not evidence against it, agreement is a bonus.
Schuhknecht / Sarvin tests are **characterisation**, not validation. Schuhknecht training is deprioritised (one
cancer cell line, ~1,500 drugs vs thousands of LINCS perturbations in hundreds of cell lines; pools cannot
separate more production from less consumption). The expr model is not restricted to metabolic genes: the expr
launchers use all DE genes on `networks/gene` (no GEM layer), so it is a generic transcriptome → target model and
the stoichiometric reading is MetChange's specific contribution.

The chain and what tests each link (report 05, 9.3): (a) signature → target = **done** (CV target recovery,
8–10 points over control); (b) transcriptome → MetChange = the LINCS genetic-perturbation test below; (c) enzyme
transcript → enzyme protein = drug proteomics atlas, overlap not checked; (d) capacity → flux = the claim itself,
prospective tracing only.

**To do, in order:**
1. **LINCS genetic perturbations of metabolic enzymes (link b, the one real label test).** `trt_sh` / `trt_xpr`
   (knockdown / knockout) and `trt_oe` (overexpression) signatures of GEM enzymes, in the cell lines with the
   most coverage. Two readouts: (i) **label-level** — does knockdown of an enzyme read as less production of its
   products (and overexpression as more), including 1–2 reactions away from the perturbed gene? (ii)
   **model-level** — query the trained model with those MetChange signatures: does the perturbed enzyme rank near
   the top? Exact ground truth, no polypharmacology, no new data modality. Partly circular (MetChange reads the
   perturbed transcript itself), so report the 1–2-reaction-away version separately.
   *First step (no GPU):* count how many GEM enzymes have `trt_sh` / `trt_xpr` / `trt_oe` signatures per cell
   line, and whether A549 or another line has the best coverage. Then compute MetChange labels for them
   (`scripts/compute_labels`, same recipe as the drug labels) on the cluster.
2. **MetChange around the target reaction for known metabolic inhibitors** (characterisation, cluster: needs the
   **corrected** labels, the laptop copy is the old uncorrected one). LINCS drugs whose target is a GEM enzyme
   (MTX / DHFR, MPA / IMPDH, leflunomide-teriflunomide / DHODH, statins / HMGCR, pemetrexed / TYMS, …): restricted
   to the target reaction's substrates and products, is MetChange specific to the drug vs other drugs, and in
   which direction (block vs compensation)? Expression read the same way as comparison. Tests the compensation
   hypothesis for the negative slope in 8.2.
3. **Local Schuhknecht test** (characterisation, laptop, runs on the corrected matched matrices in
   `results/mc_schuhk_assoc/limma_rank_*`): drugs whose annotated target **is a GEM enzyme** (see the warning
   below: *not* the `validated_metabolic` column), only the target reaction's substrates / products. Does
   Schuhknecht show substrate up / product down at 48 h, and what does MetChange say at the same metabolites?
   Caveat: few drugs will have an annotated, measured substrate or product.
4. **Knowledge-based enzyme benchmark** (link a, interpretable): expected substrate-up / product-down query
   around the enzymes of drugs with simple documented mechanisms (DHFR, DHODH, IMPDH, HMGCR, TYMS, PAICS, UMPS,
   complex I) — do the known inhibitors rank at the top? Self-built from the GEM, so partly circular with the
   network; cheap and easy to read.
5. **Inborn errors of metabolism** (link a, external ground truth): deficiency → biomarker metabolites
   (Shlomi et al., Mol Syst Biol 2009; IEMbase, OMIM / HMDB). Query = biomarker profile, answer = the deficient
   enzyme. Caveat: plasma / urine pools.
6. **Prospective**: isotope tracing / exometabolomics for a few top predictions (link d, the only direct test of
   the claim).

**0. (Do first, cheap) Protein-level target rank in the Schuhknecht benchmark.** `scores.npz` only holds the
drug x drug matrices, but `benchmark_schuhknecht.py` already computes the ~924 candidate-protein scores
internally. **To do: save them and report, per drug, the rank / percentile of its true target protein(s)** among
the candidates (and the mean over the 5 folds for never-trained drugs), with the same permutation null and the
same group split. Small code change + one rerun (~2 min per run on the GPU, runs are on scratch). Why it matters:
the model's output *is* protein scores; ranking drugs is a downstream aggregation, and "lands near the right
enzyme" is both the real task and far better powered than drug retrieval. Report it for all 646 drugs and
separately for the GEM-enzyme-target subset. If the model does not get near the right enzyme here, no new dataset
will rescue it.

**Warning: `validated_metabolic` is LLM-generated and is not "targets a metabolic enzyme".** It flags drugs
*known to have metabolic effects*, which may be indirect (signalling, transport, toxicity), not drugs that
directly inhibit an enzyme. Do not use it as the metabolic-target subset in any of the tests above. **Derive the
subset properly instead:** intersect each drug's annotated targets (the target table the training / benchmark
scripts already use: `--targets`, gene symbols) with the **GEM enzyme genes** (`data/networks/combined`
`gpr_ensembl` / the GEM gene list used by `check_metabolic_coverage.py --gem-genes`), and keep drugs whose target
is a GEM enzyme *and* whose mechanism is inhibition (ChEMBL `mechanism_of_action` / DrugBank actions where
available). Report `validated_metabolic` only as a secondary, clearly-labelled annotation. The same correction
applies to report 05 sections 6, 6.2, 9.1 and 9.2, where the 70 `validated_metabolic` drugs were used as the
metabolic subset, and to the seed-1 `validated_metabolic` spec 0.418 (p = 0.085), which is therefore not even
about enzyme-targeting drugs.

**Sarvin et al. (Nat Commun 2026, report 05, 9.2).** 76 known drugs, 30 metabolic enzymes, HCT116 (a LINCS cell
line), 6 doses × 6 times (10 min–24 h) — the best-matched public dataset (Schuhknecht is A549 at 48 h only).
**But figshare 29268887 holds raw data only:** `Data_mzxml.zip` (8.9 GB mzXML) + `Sample description.xlsx`
(downloaded to `data/validation/sarvin/`: 76 drugs by catalogue number / name, control coding). No processed
intensity matrix, no annotation, no drug → target table (that is Supplementary Data 1 of the paper). Full
untargeted processing (XCMS / MZmine + HMDB annotation + control normalisation) is days of work and would not
reproduce their values exactly; a **targeted extraction of the ~30 target enzymes' substrates / products** from
the mzXML is the lighter option, and enough for test 2–3 style comparisons. Supplementary Data 1 is needed either
way.
**Usable as a Schuhknecht-style benchmark? Only weakly, and not first.** Crude name match of the 76 drugs:
33 in LINCS, **22 evaluated in our A549 CV** (methotrexate, brequinar, leflunomide, mycophenolate mofetil,
pemetrexed, merimepodib, 3 statins, raltitrexed, pralatrexate, TEPP-46, AGI-6780, cytarabine, gemcitabine,
capecitabine, clofarabine, metformin, lonidamine, fomepizole, disulfiram, indoximod). SMILES matching should add
5-FU, mycophenolic acid, 2-deoxy-D-glucose, telaglenastat and the other statins: ~45 in LINCS, ~30 scorable.
Drug-ranking with n ~ 30 has far too little power (Schuhknecht had 646 and was flat), and it is HCT116, so it
needs either cross-cell-line transfer or HCT116 labels + a second training run. **Its value is the enzyme ground
truth** (30 known target enzymes), i.e. the protein-level metric of item 0, which works at n ~ 30 because each
drug is one trial against ~924 candidates. Processing only the useful conditions (8 h / 24 h at ~10 uM, ~300
mzXML files instead of 2,736) keeps it bounded. **Order: do item 0 on Schuhknecht first; process Sarvin only if
that shows signal.** Either way it is a benchmark of the tool, not validation of the labels (9.3): a positive
result counts, a null is uninformative.
Their model predicts the target *pathway* (373 candidates: top 1 29%, top 3 47%, AUC 0.87) from a
time × dose metabolomics matrix, and the enzyme was pinned down by rescue experiments — complementary to ours
(proteins from transcriptomes, free for thousands of LINCS compounds).

R on the cluster: conda env `r-4.3.1` (`launch_schuhknecht_mc_assoc_A549.sh` activates it for the R step).

**6.3 Size the data ceiling** (no GPU, ~a day):
1. PCA of the A549 expression matrix: how many components explain 95% of variance across
   perturbations. That is the budget of distinguishable perturbation states.
2. Fraction of the GPR genes MetChange uses that are landmarks vs CycleGAN-inferred. If mostly
   inferred, the metabolic labels are a smooth function of the 978 landmarks.

**6.4 Then** freeze the recipe (done 2026-11-15: `_min15`), and only after that revisit version (b), the prior + residual
model, or per-protein standardisation of the dynamic features. (Done on 2026-11-14: cap check —
does not bind; second seed — replicates, smaller; port to `tf` / `expr` — report 03.)

## 7. Open issues

- **Per-drug plain ranking is still prior-dominated** (section 3). This is the main scientific
  gap for the tool, not a bug.
- ~~MetChange reads the signature about half as well as expr~~ — **resolved**: the label sign bug
  (report 05, section 5; corrected MetChange ≈ expr ≈ TF).
- **Laptop expression / MetChange labels are an old partial set** (`results/labels/A549/expr/`
  and `metchange/rank/`, 2 Nov: 699 of the 3,970 scored perturbations). Analyses needing labels
  must run on the cluster or after syncing the scratch copy. (The MetChange `_min15` control's
  low folds 0–1 did not recur in seed 1.)
- **MetChange similarity-check anomaly — probably explained** (report 04, 1.1). TOP1 / JAK1 / XIAP
  shared-target pairs were anticorrelated and the split matched flip status exactly: consistent
  with the uncorrected home labels (camptothecin flipped, topotecan not). ATP1A1 (mixed-flip,
  r = +0.92) is not explained. Confirm by rerunning the similarity check on the scratch labels.
- **Model does not use edge signs.** `propagation_operator` gates on weight × activity only;
  the two channels propagate independently.
- **Signatures are discretized** (`np.sign`, with `--ignore-qvalue --top-n 100`), while the
  similarity check uses continuous Spearman. A `--signature-weights {sign,value,t}` option was
  proposed, not implemented.
- **`--no-globals` also drops `log_n_sig`** (signature size), which is *not* design information
  — it is computed from the signature and travels with it in the shuffled control. Splitting
  the flag (e.g. `--no-dose-time`) would keep the most direct effect-strength input. Small
  change, not done.
- **Validation `spec` used for epoch selection is unstratified**, so it can reward dose/time in
  any run that keeps those inputs. Harmless for `--no-globals` runs.

## 8. Uncommitted

- Report 04 (`reports/2026-11-16_labels_network_coverage.md`) and the correction banner on report 03.
- `scripts/predict_targets/check_metabolic_coverage.py`: `--coverage` (reuse a `coverage.tsv`),
  fixed 0 / 1–99 / 100 bins for `n_sig_topn`; README lines.
- `scripts/predict_targets/submit_train_w_shuff_mc_min15_fixedlabels.sh` (new).
- `launch_predict_targets_A549_{mc,tf,expr}_noGlobals_lr1e-3_min15_qsig{,_shuff}.sh`, wrapper
  `submit_train_w_shuff_noGlobals_lr1e-3_min15_qsig.sh`, variants registered in
  `submit_with_shuffle_control.sh`, README row.
- `--qvalue-top-n` in `metabolic_target_prediction.py`; `_qsigtop` launchers (mc, tf, expr, with `_shuff`),
  MetChange `_qsigtop_sdscale`,
  wrapper `submit_train_w_shuff_noGlobals_lr1e-3_min15_qsigtop.sh`, registered, README row.
- Report 05 (`reports/2026-11-17_signature_selection_scale.md`).
- Report 05, section 5 (corrected MetChange results); `results/.../metabolic_coverage_fixedlabels/`.
- `metabolic_target_prediction.py`: `--signature-scale {none,sd,mad}`, `--signature-scale-floor`;
  launchers `launch_predict_targets_A549_mc_noGlobals_lr1e-3_min15_sdscale{,_shuff}.sh`, wrapper
  `submit_train_w_shuff_mc_noGlobals_lr1e-3_min15_sdscale.sh`, variant registered, README row.
- This `STATUS.md` update.

## 9. Working notes

- A safety classifier has repeatedly truncated long prose answers in this project (three times
  on 2026-11-13). Code, tables, numbers and tool output were never affected; the pattern was
  broad written assessments. Keeping write-ups in `reports/` and chat answers short avoids it.
  Starting a new session helps because it drops accumulated context — hence this handoff.
- Laptop has `python3` (no `python`), no torch, bash 3.2 (no `declare -A`), and no network
  data, so only unit-level tests and result analysis are possible there.
