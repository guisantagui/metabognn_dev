# Project state — 2026-11-17

Handoff note for continuing in a new session. Everything described here is committed up to
HEAD `d11744b`, except the files listed in section 8.

> **2026-11-17 — corrected MetChange matches expr and TF** (report 05, section 5). Retrained on the
> corrected scratch labels (`_min15`, seeds 0 / 1), MetChange's per-drug gap is 8.3 / 7.7 points
> (buggy: 4.2 / 2.7), not significantly different from expr (7.8 / 10.4) or TF (9.0 / 9.3) paired
> by drug (p = 0.09–0.79). On perturbations with a full MetChange response it beats expr on the
> same perturbations (14.6–14.7 vs 10.5–12.1). The deficit of reports 03–04 was the sign bug
> (report 04: the home MetChange labels were the uncorrected limma output). Open: whether a
> scale-free signature selection adds to this (`_sdscale`, `_qsig`, running).

**Read first:** `reports/2026-11-17_signature_selection_scale.md` (report 05 — q-values unused; MetChange's
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
- **Feature scale (report 05):** MetChange log2FC is not standardised per metabolite (the MetChange
  paper z-scores against controls), so the current top 100 is mostly a fixed set of high-variance
  metabolites (selection frequency vs metabolite SD: Spearman 0.89; 33 metabolites in > 50% of
  signatures; different drugs share 22% of features vs 2–3% ranked by p). Expr is affected much
  less (5%). Fixes: `_qsig` (scale-free via significance) and **`--signature-scale sd`**
  (implemented: per-feature SD across perturbations, not centred, floored; z-scoring against
  controls rejected, many controls are constant in many features).
  `bash submit_train_w_shuff_mc_noGlobals_lr1e-3_min15_sdscale.sh` (MetChange, seed 0, with control).
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
   Schuhknecht ids also need mapping to `base_met` (MAM) ids.

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
