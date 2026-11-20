# Project state — 2026-11-19

Handoff note for continuing in a new session. Everything described here is committed up to
HEAD `d11744b`, except the files listed in section 8.

> **2026-11-19 — validation datasets compiled and mapped (report 07).** `map_smiles.R` rewritten (salt
> stripping, all pairs kept, stereochemistry labelled not assumed). **Attene-Ramos (mito membrane potential,
> HepG2, ~10k compounds) is the best-powered external set: 904 LINCS HepG2 compounds mapped**, HepG2 MetChange
> labels already on scratch. **Hallinger (Seahorse) is limited to its 270 named tier-2 actives** — the 1,042-
> compound tier-1 screen was never registered in invitrodb (checked: analytical_qc, assay_annotations, pubchem
> per-aeid, sc1_sc2 from the 7.5 GB summary; 0/1,051 spids), so only the authors have the plate maps. Currency
> metabolites (ATP, NAD+/NADH, O2) **are** MetChange features but **not** network nodes: label-level tests can
> use them, model queries cannot. See section 6.2c.

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

**Read first:** `reports/2026-11-20_phenotype_confounders.md` (report 08 — are the report 07 phenotype
associations chemistry or metabolism: the ToxCast-subset finding that explains Hallinger's null, mechanism-class
stratification, the absent potency gradient, and what may/may not be claimed), then
`reports/2026-11-19_validation_datasets.md` (report 07 — external validation sets, structure
mapping, what each can test), then `reports/2026-11-18_recipe_state.md` (report 06 — where the recipe stands, what is settled and
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

**6.0 Running the phenotype association (refactored 2026-11-20, report 08 section 8).**
`sbatch scripts/mc_phenotype_assoc/launch_assoc_HEPG2_slurm.sh` (or `..._A549_slurm.sh`) -> exports `cell` ->
`launch_assoc_slurm_generic.sh` (2-task array: 0 MetChange, 1 expression) -> `launch_assoc_metchange.sh` /
`launch_assoc_expr.sh` (no cell line in either) -> `mc_phenotype_assoc_common.sh`. One-off per cell:
`cell=<CELL> bash scripts/mc_phenotype_assoc/prepare_assoc_inputs.sh`, which needs
`lincs_unique_perts_<CELL>.tsv` plus `cell=<CELL> bash scripts/benchmarking/map_{attene,hallinger}_lincs_perts.sh`.
Paths (all under scratch on the cluster, repo on the laptop): phenotypes **per cell** at
`results/mc_phenotype_assoc/<cell>/phenotypes/`, sets **per project** at
`results/mc_phenotype_assoc/metabolite_sets/`, outputs at `results/mc_phenotype_assoc/<cell>/`. Sets are now
built WITHOUT `--features` so they are genuinely cell-independent; absent members are dropped at scoring time
and reported. Verified equivalent on HEPG2 (same 40 sets, same numbers). **An A549 run is the cross-cell-line
CONTROL, not a replication** - the phenotypes are measured in HepG2 regardless of which cell line is mapped.
The slurm wrappers are **submittable from anywhere in the repo** (they walk up from `SLURM_SUBMIT_DIR` to find
the root and export `repo_dir`; `$0` is unusable there because Slurm spools the script) and log to
`./mc_pheno_<cell>_%A_%a.out` in the submit dir. **Bug found and fixed 2026-11-20 while testing the A549
path:** `build_phenotype_tables.py` fell back to the legacy unsuffixed mapping for ANY cell, so `cell=A549`
silently built A549 tables from the HEPG2 mapping (identical 849/147/480 counts). The fallback is now allowed
only for HEPG2; other cells require `mapping_<cell>/` and the run stops with the commands to build it.

**Running ALL cell lines (report 08, 8.4-8.5).** `bash scripts/mc_phenotype_assoc/submit_assoc.sh` discovers
which cell lines actually have limma labels (`check_label_availability.py` over `core_cell_lines.tsv`; checks
`log2FC.csv` required, `qvals.csv` / `fit_info_nsig.csv` optional, and warns when neither exists because then
`--mode sig` AND the responsiveness stratification are unavailable), writes a `(cell, arm)` manifest, sizes
the Slurm array to it and submits `launch_assoc_manifest_slurm.sh`. `cells="HEPG2 A549"` restricts,
`dry_run=true` previews. Cell lines missing labels or phenotypes are skipped, not failed. **Logs land in
`scripts/mc_phenotype_assoc/logs/`** (kept in git: Slurm opens `--output` before the job runs, so a mkdir
inside the script is too late); `submit_assoc.sh` passes that path absolutely, the single-cell wrappers use a
relative `logs/` so submit those from `scripts/mc_phenotype_assoc`.
**Project root is configured in ONE place:** `scripts/mc_phenotype_assoc/project_paths.conf` (candidate
`scratch_proj_dir` lines, first existing wins, repo as fallback), sourced by common and submit_assoc and
parsed by check_label_availability.py. It exists because the path had been hardcoded in all three and the
checker's copy defaulted to the repo when `$labels_base` was unset - so on the cluster it read `$HOME` and
reported "no directory" for every cell line whose labels are on scratch. Fixed 2026-11-20.
**One mapping instead of nine:** `cell=ALL bash prepare_assoc_inputs.sh` maps against
`lincs_unique_perts.tsv` (the union; all 5,284 HEPG2 pertnames verified inside it) and writes
`<assoc_base>/phenotypes_ALL/`, which is cell-independent because the association script keeps only the drugs
with labels in the cell line being run. Common prefers a per-cell dir and falls back to the union, printing
which. Union mapping: Hallinger 176 compounds vs 115, Attene 1,865 vs 918; union tables hold 1,763 Attene and 177
Hallinger drugs. **NOT a drop-in replacement** - restricted to drugs with HEPG2 labels the union table gives
144 vs 147 MMP-down, 475 vs 480 inactive, 27 vs 29 Hallinger-active, because mapping against 34,419 perts
instead of 5,494 lets compounds match further pert_ids and a few aggregated labels resolve differently.
**So for ANY cross-cell comparison use `use_union=true bash submit_assoc.sh`**, which forces one identical
table onto every cell line; otherwise cells mix per-cell and union tables and a 2-3% assignment difference
would look like a cell-line effect. The per-cell HEPG2 table stays the one reproducing reports 07-08.
**Data note:** `lincs_unique_perts.tsv` originally had 4 header fields for 5 data columns (no `V1`), which R
tolerates as row names but silently misaligns every Python `DictReader`. Regenerated with `V1` 2026-11-20 -
check other files written by that code path.

**6.0a DO THIS FIRST — finish cleaning the lumped / placeholder species out of the network.**
The existing filter (a hand-made string-match list of aggregated species, plus the currency metabolites)
catches the obvious cases but has gaps. Found on 2026-11-20 while chasing the O-glycan result
(report 08, 3.3): **`Protein-Linked Serine Or Threonine Residue (O-Glycosylation Site)` (MAM03395g) is still a
metabolite node**, as are the **lumped glycan classes** `core 2`-`core 8`, `T-antigen`, `Tn-antigen`,
`sialyl-T`/`sialyl-Tn`/`disialyl-T antigen`, `glcnac-alpha-1,4-core 1/2`, `F1alpha`, `G00031`, `G00032`.
The three true `[protein]` entities ARE correctly absent, so the filter works, it is just incomplete.

Why it matters: these are not single chemical species, so a MetChange "production capacity" for them has no
measurable referent; they are low-degree, so their value is close to a readout of a few local enzymes; and
they can be selected as sources when significant and large. In report 08, 3.3 a 21-member set of exactly
these carried an association as strong as the TCA one (per-member median |AUROC-0.5| 0.097 vs 0.100).

A second, name-based candidate list is in **`data/networks/combined/lumped_candidates_for_review.tsv`**
(560 distinct names over 772 metabolite nodes, 9.8%, by category) — useful only for catching aggregates that
*do* have a clean formula, e.g. the lipoprotein particles. **Its patterns over-match on purpose — it is a
review list, not a deletion list** (it flags things like a tetrahydrobenzo-phenanthro-oxirene-diol purely for
containing a bracket). Categories in it where the aggregate reading is clear:

| category | nodes | examples |
|---|---|---|
| glycan core / antigen classes | 66 | `core 2`-`core 8`, `T-antigen`, `Tn-antigen`, `I-Antigen`, `F1alpha` |
| KEGG glycan ids | 37 | `G00005`, `G00020`, `G00031`, `G00032` |
| lipoprotein particles | 15 | `HDL`, `LDL`, `VLDL`, `IDL`, `chylomicron`, `* remnant` |
| protein-linked residues / sites | a few of the 49 | the O-glycosylation-site residue above |
| explicit pool, polynucleotide | 5 | `Gm4-Pool`, `polynucleotide` |

Needs judgement, probably KEEP (chemically defined despite the pattern match): `[ACP]` thioesters
(`(2E)-decenoyl-[ACP]`...), the GAG and glycolipid series with specified linkages
(`III3Fuc-nLc4Cer`, `chondroitin sulfate A (GalNAc4S-GlcA)`...), chain-specified lipids
(`cholesterol-ester-10,13,16-docosa`).

**Criterion to use — verified, not just suggested: the Human-GEM molecular formula.** Keep a node if it
denotes a single chemical species with a fully specified formula; drop it if the formula is empty, contains an
`R` (unspecified side chain) or an `X` (attachment point), or is a `(...)n` repeat unit. The formula is in
`Human-GEM.yml` per metabolite (NOT in `metabolites.tsv`, which has no formula column).

**This catches every culprit found on 2026-11-20**, all via the `X` atom: `T-antigen` C14H24NO10**X**,
`Tn-antigen` C8H14NO5**X**, `core 2` / `G00032` / `F1alpha` / `glcnac-alpha-1,4-core 1` C22H37N2O15**X**,
`core 6` / `G00031` C16H27N2O10**X**, `sialyl-T` C25H40N2O18**X**, `disialyl-T` C36H56N3O26**X**, and the
O-glycosylation-site residue whose formula is simply **`XH`**. 11 of 11 flagged.

**`data/networks/combined/lumped_by_formula.tsv`** is that list applied to the network: **659 of 7,841
metabolite nodes (8.4%), 398 distinct names** — 410 with an `R` atom, 246 with an `X`, 3 with no formula.
Prefer this file over the regex list below; it is principled and complete where the string match was neither.

Decide separately what to do with the `R`-atom group, which is a different kind of case: generic lipids such
as `1-alkyl-2-acylglycerol`, `1-phosphatidyl-1D-myo-inositol-3-phosphate`, `1-(1-alkenyl)-sn-glycerol`. They
are unspecified in composition for the same reason, so the same argument applies, but there are 410 of them
and dropping them changes the lipid coverage of the network substantially.

**Knock-on effects to plan for before running it:** the label matrices do NOT need recomputing (MetChange is
computed per metabolite and limma ran over the whole matrix, which is why these appear in `log2FC.csv` at
all), but **`select_sources` output changes, so model runs and any analysis keyed on node identity would need
redoing**. Decide the criterion first, apply once, and note the network version in the run metadata.


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
   **RESULT (2026-11-18, report 06 section 6.1): the kNN baseline ranks the true drug BETTER than the model on
   `rank_probc`** (median percentile 0.29-0.31 vs 0.38-0.41, controls all ~0.50; the kNN wins for 55-61% of drugs,
   p ≤ 0.012 MetChange, ≤ 7e-14 expr / TF). Full profile ≈ same input (only MetChange +1.7). Model and kNN succeed
   on different drugs (per-drug Spearman 0.18-0.31) and a rank-average of the two gives 14-16 points. **Caveat:**
   with k = 25, 68-78% of candidates score exactly 0, so `rank_prob` and `spec` are not comparable and the gap
   statistics are noisy; **next is a dense baseline** (every training drug votes, soft weight) before this goes in
   any presentation. **Dense version written and launched 2026-11-18** (`--beta 2 5 10`, launcher default; report 06
   6.2): no zeros except token-level orphans (~7-12% of candidates, flagged per query); check the conclusion
   holds across B and on `rank_prob` / `spec` too. If it holds, the network does not earn its place with the current data, and the priorities
   below change (stacking the kNN vote into the model, or rethinking the architecture, before more data).
   **The predictive version was written (2026-11-18) and run the same day:** `scripts/predict_targets/
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
7. **Seahorse (OCR / ECAR) as a query** (design in report 06 section 7, nothing run). Real flux readout, so it
   tests the capacity claim that pools cannot; coarse (whole-cell respiration / glycolysis only). Query: not
   "every TCA + ETC metabolite up" (at higher flux intermediates often stay flat or fall) but the **endpoints**
   (O2 / NADH / acetyl-CoA consumed, CO2 / ATP produced; glucose consumed, lactate produced), and the two rates
   **together** — OXPHOS up (OCR up, ECAR down), glycolytic shift (the reverse), both down (toxicity, a negative
   control). Two benchmarks: **(b) enrichment of each drug's significant MetChange metabolites in glycolysis /
   TCA-ETC vs its measured OCR / ECAR direction** (one value per drug, well powered, reuses
   `scripts/schuhknecht_mc_assoc/`, tests the labels) **first**, then **(a) query -> drug rank** (the tool's task,
   but only 2-3 queries, so it is a rank-enrichment test of the known positives against a permutation null).
   **Binding constraint: data.** Public Seahorse at screen scale is scarce and usually only in summary tables; for
   (a) curated known ETC / OXPHOS / glycolysis inhibitors and mitochondrial-toxicity screens (Tox21,
   glucose-galactose) supply the positives without raw rates. Check what exists for a LINCS cell line first.

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

**6.2c External validation datasets (report 07, 2026-11-19).**

| dataset | readout | cell | mapped to LINCS | status |
|---|---|---|---|---|
| Schuhknecht | metabolomics (pools) | A549 | 750 | done, null (report 05, 8.2) |
| **Attene-Ramos** | **mito membrane potential** | **HepG2** | **890** | **next: section 6 of report 07** |
| Hallinger | Seahorse OCR (flux) | HepG2 | **115 (29 active / 86 inactive)** | **all 1,042 tier-1 chemicals recovered, report 07 sections 8-9** |
| Sarvin | metabolomics, dose x time | HCT116 | ~22-30 | raw mzXML only; enzyme ground truth |
| Lackner | isotope labelling | breast | 1 drug class | case study |
| LINCS shRNA / CRISPR / OE | genetic perturbation | many | — | still the one test of the labels (6.2b item 1) |

> **2026-11-19 (later) — Hallinger solved, two mapping bugs fixed (report 07, sections 8-9).** The paper's
> Supplemental Table 2 lists **all 1,042 tier-1 chemicals with activity calls** (243 active / 799 inactive) and
> the 243 actives' **mechanism and inhibited ETC complex**, so the EPAPLT plate maps are not needed; parsed by
> `scripts/benchmarking/parse_hallinger_supp.py`. **Hallinger -> LINCS HepG2 is 115 compounds, 29 active / 86
> inactive** (19 ETC inhibitors, 3 uncouplers, 4 ATP synthase; rotenone = complex I), up from 42 actives-only.
> **`map_smiles.R` had two bugs** (9.2): `.I` under a filter gave subset positions used as table rows, so right
> structures carried wrong ids; and explicit `[H]` in DSSTox SMILES broke canonicalisation, losing true matches
> (rotenone, estradiol). Both fixed. **Attene-Ramos was re-mapped: 890 compounds, and 0 of the old 1,509 pairs
> survive** - anything built on the old `mapped_drugs.csv` must be redone (old file kept as
> `mapped_drugs_BUGGY_2026-11-19.csv.bak`). **Schuhknecht checked and is clean** (report 07, section 10): that
> mapping is separate code which matches on values, so the report 05 null stands - 998 of its 1,040 pairs
> reproduce and the only drug lost is Gramicidin, a `NA`-matched-`NA` artefact. Re-mapping it would *gain*
> 212 drugs (1,251 total), since explicit `[H]` (228 of 1,520 Schuhknecht SMILES) and salts blocked them.

> **2026-11-19 (later) — first non-null external association (report 07, section 11).** Label-level test of the
> HepG2 MetChange labels against measured mitochondrial phenotypes, `scripts/mc_phenotype_assoc/` (40 Human-GEM
> metabolite sets, per-drug `z` / rank statistic, Mann-Whitney + bootstrap + permutation + **size-matched
> random-set control**, no imposed sign, model activity filter deliberately NOT applied).
> **Attene-Ramos: MMP-down drugs have lower MetChange TCA capacity** (`tca mito_nohub` AUROC 0.364, q = 7e-6,
> random-set p = 0.010; non-hub variants strongest, glycolysis weak at 0.43, **oxphos / ETC flat**), and
> **MMP-up drugs have higher energy-currency capacity** (0.663, q = 2e-4). The 24 cytotoxic drugs are flat, so
> it is not a cytotoxicity readout. **Hallinger is null** (29/86; 19 ETCi vs 86 — same direction, CIs cross 0.5,
> uninformative at that n). **Confound, carry it with the result:** MMP-active drugs are 3.4x more likely to have
> a responsive MetChange signature (31% vs 9%), and within the *controls alone* responsive drugs already show the
> low-TCA pattern (0.355, p = 0.0015). Stratified, TCA holds at ~0.38 in both strata (non-responsive p = 2e-4) and
> MMP-up/currency is clean (0.698 non-responsive, confounder absent). **Decisive next control: the same test on the
> HepG2 expr labels with GPR gene sets** - needs `log2FC.csv` for HEPG2/expr (the laptop copy has only
> `qvals.csv` / `basal.csv`), then `launch_assoc_expr_HEPG2.sh` runs it; sets via
> `build_metabolite_sets.py --emit genes` (Ensembl both sides, no mapping).

> **2026-11-19 (later still) — three modalities (report 07, section 12).** `--mode value|sign|sig`;
> refactor verified (all 80 shared statistics reproduce exactly). **The directional reading confirms TCA and is
> the cleanest evidence**: `sign_dev tca cyto` AUROC 0.415 (q = 0.041, rand = 0.020), non-responsive stratum
> 0.432 (p = 0.03) and **the responsiveness confounder is absent** (0.453, p = 0.30) where it was present for
> the magnitude statistic. `sign_mean` null, so TCA is down *relative to each drug's own mean sign*.
> **The `sig` modality is a responsiveness artefact - do not use it pooled**: `tca mito_nohub` q = 4e-10 pooled
> but AUROC exactly 0.500 (p = 1, all ties) among non-responsive drugs and 0.091 (p = 2e-71) for responsiveness
> within controls alone; restricted to responsive drugs nothing survives (best p = 0.079). Always pair
> `--mode sig` with `--min-nsig 1`. **The apparent cytotoxic signal is the same artefact** - the clean cytotox
> control of section 11 holds only for value/sign. Expr baseline started but only in `sig` (no `log2FC.csv`):
> no gene set beats its random-set control (rand >= 0.25), suggestive but not like-for-like.
> **Launchers:** `scripts/mc_phenotype_assoc/launch_assoc_{metchange,expr}_HEPG2.sh` +
> `launch_assoc_HEPG2_slurm.sh` (2-task array); config in `mc_phenotype_assoc_common.sh`, which auto-detects
> the scratch base on the cluster.

> **2026-11-20 — confounder analysis of the phenotype association (report 08).** **The association is carried
> entirely by the drug-like part of the library**: restricted to the 115 compounds shared with Hallinger
> (= ToxCast) the Attene effect is 0.447 (p = 0.39, n = 38/55), outside them 0.352 (p = 1.8e-06, n = 109/425).
> **That, not n, is why Hallinger is null** - its point estimates for the Attene-signal sets are 0.47-0.55, not
> 0.37, and power at 29/86 for 0.371 is 54%. **Mechanism class does NOT explain it:** within 9 coarse MoA
> families `tca mito_nohub` is 0.380 (Holm p = 0.0054) and `energy_currency mito` 0.657 (Holm p = 0.0010),
> essentially unattenuated; the negative control `oxphos mito` stays null everywhere. Narrow (fine) MoA strata
> give p = 0.13-0.57 but are 27 tiny strata in a receptor-drug population - untestable, not negative.
> **No potency gradient survives correction** (best Holm 0.076, and the negative-control set shows an
> equal-sized gradient), so there is no dose-response support. Similarity matching is uninformative (median
> Tanimoto 0.27). **An earlier within-MoA reading of mine ("fits a compound-class effect") was an artefact of
> double-counted strata, an unpaired test and no correction - withdrawn, see report 08 section 3.**
> **Dose/potency line exhausted (report 08, 5.1-5.3):** `Ratio Potency` is the qHTS AC50 (half of the compound's
> own fitted max, ~15-point titration to ~92 uM); across actives it spans 3.9 decades but the interquartile
> range is only 0.75, and **the max LINCS dose has an interquartile spread of 0.00 decades (~all 10 uM)**.
> **48% of MMP-down cases were profiled BELOW their own AC50** - a real design fact - but splitting on it gives
> identical effects (`tca mito_nohub` 0.361 exposed vs 0.367 sub-threshold, p = 0.68, rho = +0.002), and where a
> gradient appears (`tca cyto sign_dev` rho = -0.193) the negative control `oxphos mito` shows a stronger one
> (-0.217). The test is also NOT independent of the potency test (constant dose => driven by AC50 alone) and the
> MMP AC50 is not a threshold for a 6-24 h metabolic change. So this dataset cannot provide a dose-response test
> at all; no further potency analysis of it is worthwhile.
> **DECLARED NEGATIVE CONTROLS change the interpretation (report 08, 3.1).** `oxphos mito` was NEVER a valid
> negative control - it is the most on-target set for this phenotype, was designated post hoc because it came
> out null, and is 11/13 hubs; it is a POWER DIAGNOSTIC and is relabelled as such in reports 07 and 08. Five
> proper controls (glycan / glycosaminoglycan subsystems, criteria fixed before looking):
> **4 of 5 are null, but `negctrl_o_glycan` gives AUROC 0.596-0.631 (p = 4e-4 to 1e-6), beats its size-matched
> random null (rand = 0.0099, same as TCA's), and survives coarse-MoA (0.613, Holm 0.0097) and similarity
> matching.** Same magnitude as TCA (|0.13|), opposite direction. **So beating a size-matched null is NOT
> sufficient to implicate a pathway, and the defensible claim is a MULTI-PATHWAY difference between
> phenotype-active and inactive compounds (TCA down, O-glycan up, 34/40 primary sets below 0.5), not a
> TCA-specific metabolic response.** **Tested the obvious rescue and it FAILED (report 08, 3.2):**
> if o_glycan were coupled to TCA via the hexosamine pathway its nucleotide sugars (UDP-GlcNAc/GalNAc/Gal,
> CMP-NeuAc, serine) would carry the signal - they are **null** (0.488, p = 0.66) and the signal is entirely
> in the **glycan dead ends** (core 2-8, T/Tn/sialyl antigens, degree 2-12: 0.649, p = 4e-08).
> **But the same split on TCA goes the other way:** `tca_core_cycle` (citrate, AKG, succinate, fumarate,
> malate, OAA, isocitrate, cis-aconitate, succinyl-CoA; degree 8-39) gives 0.367 (p = 1.1e-06, rand = 0.005),
> stronger than `tca_peripheral` 0.393. So the two associations are **different in kind** - TCA in high-degree
> central metabolites (a capacity interpretation means something), o_glycan only in low-degree terminal
> structures (MetChange there is ~a readout of a few glycosyltransferase genes). **General diagnostic from now
> on: for any significant set, check whether the signal sits in its high-degree or its low-degree members.**
> **Refined (report 08, 3.3): per member the two effects are the SAME SIZE** (median |AUROC-0.5| 0.100 TCA vs
> 0.097 o_glycan; 8/9 vs 21/21 members at p < 0.05), so o_glycan is not the weaker finding - the difference is
> interpretability. **Important caveat on `rand_p` throughout reports 07-08:** the o_glycan cores are
> near-collinear (mean pairwise rho 0.799 vs 0.511 for TCA), i.e. ~1 effective d.o.f. measured 21 times, so a
> near-collinear set beats a size-matched random set partly BECAUSE its average is less noisy. The random-set
> control assumes comparable within-set correlation; always report per-member statistics alongside set ones.
> The o_glycan set is mucin-type O-glycan cores / Thomsen-Friedenreich antigens, mostly **lumped glycan classes**
> plus 4 **GEM placeholders** (`[protein]`, `[protein]-L-serine`, the O-glycosylation-site residue) which carry
> the same AUROC as the real glycans - so it reads as **glycosyltransferase expression**, not metabolite capacity
> (its nucleotide-sugar precursors are null, so substrate supply is unchanged).
> **Placeholders cannot reach the model** (user, 2026-11-20): the 291/5,707 pool and placeholder features are in
> the limma matrix only because MetChange is computed for every model metabolite and limma ran on the whole
> matrix; they are filtered during network construction along with the currency metabolites (hubs, no real
> biology), so they cannot be sources or queries. Verified for this set: the three `[protein]` entities are
> absent from `networks/combined/nodes.tsv` - but 18/21 remain nodes, including the O-glycosylation-site residue
> and the lumped glycan classes, so lumped pseudo-species are NOT filtered. Cheap check outstanding: how often
> pool/lumped names appear among the 100 selected sources.
> **Decisive next experiment: the LINCS genetic perturbations - a knockdown has no chemotype, no MoA and no
> dose-response requirement, so every confounder here is unavailable by construction.**

> **2026-11-19 (FINAL, supersedes the banner below) — responsiveness-matched design: MetChange survives, the
> expression baseline does not (report 07, section 16).** `matched_check.py` (1:1 greedy matching on
> log1p(n_sig), caliper 0.25, plus a residual-on-n_sig variant; balance: 140/147 matched, imbalance p = 0.64
> MetChange and 139/147, p = 0.97 expr). **MetChange `tca mito_nohub` z: pooled 0.364 -> matched 0.371
> (p = 2e-04, q = 0.0039) -> residual 0.388; `tca cyto` sign_dev 0.415 -> 0.422 (p = 0.025). Expression
> `tca gene` z: 0.413 (p = 0.0014) -> matched 0.486 (p = 0.70, q = 0.87); `oxphos gene` 0.410 -> 0.450
> (p = 0.15).** So **the expression association is essentially all responsiveness and the MetChange TCA one is
> not** - both designs agree, in both statistics, and the `oxphos mito` negative control stays null throughout
> (so it is not lost power). **This corrects section 13.5's "the baseline matches MetChange".** The pooled
> effects were similar in size but different in kind. Also: `rand_p` does NOT depend on set size (rho = +0.10,
> p = 0.52), and at matched absolute size MetChange still beats its null (42 members, rand = 0.0099) where expr
> does not (45, rand = 0.209) - so the universe-fraction worry of 13.3 is weak, and the matched design does not
> use random sets at all. **Claim as it stands (narrowed by 16.5): for MMP-DOWN (147 drugs) the MetChange TCA
> signal survives matching (0.371, p = 2e-04, q = 0.004) and the expression one does NOT (0.486, p = 0.70);
> for MMP-UP (75 drugs) BOTH survive** - MetChange via energy currency (0.651, p = 0.0014, **magnitude only**,
> sign_dev null) and expr via pyruvate genes (0.637-0.640, p = 0.003-0.004, both statistics). So the
> stoichiometric reading is distinguishable from the transcriptome on ONE phenotype, not both. MetChange
> `oxphos mito` is nominally significant for MMP-up only (0.594, p = 0.046), the hub-dominated set, consistent
> with the currency result rather than independent of it. Still not done: a second cell line for the confounder
> structure. **Null-shape comparison DONE (report 07, section 17):** at matched universe fraction the expression
> null is **1.3-1.6x wider** (median |AUROC-0.5| 0.047 vs 0.029 at ~0.22%; 0.042 vs 0.033 at ~0.74%), flat in
> size within each space. So random GENE sets separate the groups better than random metabolite sets - the
> transcriptomic group difference is **diffuse**, MetChange's is **concentrated**: MetChange `tca mito_nohub`
> |effect| 0.136 vs its null p95 0.090 (clears), expr `tca gene` 0.087 vs 0.121 and `oxphos gene` 0.090 vs 0.112
> (both inside). **Claim to use:** the stoichiometric projection CONCENTRATES phenotype-relevant signal into an
> interpretable pathway coordinate where gene space has it diffuse and responsiveness-confounded. **Stoichiometry DOES add information** (report 07, 17.4 corrects an
> earlier objection of mine in 17.3): MetChange = LP(expression, GPR, stoichiometric matrix, directionality) via
> COBRA/Gurobi, and the GEM is a prior NOT derivable from LINCS - a gene-level pathway mean cannot encode that in
> A->B->C a fall in B interrupts flux to C. The comparison is well posed for this, since the expression arm
> already uses the GEM's GPR genes and subsystems, so what differs is stoichiometry/topology/LP. Defensible form:
> *on MMP-down, expression + stoichiometry separates the phenotype (0.371 matched, q = 0.004, clears its null)
> where aggregating the SAME genes over the SAME subsystems does not (0.486, p = 0.70, inside its null)*.
> Surviving caveats: MMP-down only (MMP-up both survive); the comparator is subsystem-MEAN aggregation, so the
> claim is not "beats any expression-based method" (a topology-aware gene method such as SPIA could narrow it);
> the MoA/chemotype confounder; one cell line. **Do NOT claim MetChange "properly
> models" metabolic change** - this is a group-level shift (|AUROC-0.5| ~ 0.14), and the per-drug test (report 05,
> 8.2, rho = -0.013) points the other way. **Untested and the likeliest alternative:** MMP-active compounds may be
> enriched for a MoA/chemotype whose transcriptional response projects onto TCA; matching on n_sig equalises how
> much a drug responds, not what kind of compound it is. A MoA- or structure-matched analysis is the next step
> before presenting this.

> **2026-11-19 — the expression baseline is run; POOLED it matches MetChange (report 07, section 13; see the
> matched-design banner above, which supersedes this reading).**
> HEPG2 expr `log2FC.csv` arrived; `value`+`sign` on the GPR gene sets, same contrasts. **MMP-down: expr
> `oxphos gene` AUROC 0.410 (z, q = 0.0035) / 0.421 (sign_dev, q = 0.018) vs MetChange TCA 0.359 / 0.415;
> MMP-up: expr `pyruvate gene` 0.615 vs MetChange currency 0.660.** So **on raw separation the stoichiometric
> reading has NO advantage over reading the genes** - the kNN verdict (report 06, 6.1), now for the labels.
> MetChange looks more specific on two criteria only: it beats size-matched random sets (rand 0.005-0.02 vs
> 0.13-0.16) and its directional result is not responsiveness-confounded (0.453 ns vs expr 0.335-0.363,
> p <= 2e-7). **Both criteria are partly artefacts of comparing two feature spaces** (20,186 genes / sets of
> 23-106 vs 5,707 metabolites / sets of 5-73; and the confounder test has 265 vs 215 controls for expr but only
> 44 vs 436 for MetChange, so MetChange's "no confound" has CI ~[0.36, 0.55]). **Most interesting: the two label
> types point at different pathways** - expr at the OXPHOS *genes* (on-target, 106 genes no hubs), MetChange at
> TCA (its own oxphos set is 24/26 hubs and null everywhere). Defensible claim for a talk: the MetChange labels
> carry a TCA-capacity signal separating MMP-active from inactive drugs, at a magnitude the transcriptome also
> reaches via OXPHOS genes. Hallinger null in both arms (n = 29). Next: responsiveness-MATCHED design, random
> sets matched on universe fraction, a second cell line, and the orthogonal `APR_HepG2_MitoMembPot` readout.
> **Set construction and a multiple-testing correction (report 07, section 14):** 11 base sets x 3 compartments
> x 2 hub variants = 66 candidates, 26 below the 5-member floor (structurally: currency sets have no non-hub
> member, `oxphos` is 24/26 hubs, `ppp mito` is empty), leaving 40 - of which **39 are distinct**; the gene file
> has 10 rows but **5 distinct sets** (every `gene_nohub` equals its `gene` set). Duplicates are now collapsed,
> with an `aliases` column. **These duplicates did NOT inflate any q-value**: BH is invariant to duplicating the
> whole family (`q = p*m/rank`, both scale), so the section 13 expression q-values are correct as published; only
> the MetChange family (40 -> 39, non-uniform) shifts, in the 4th significant figure. Statistics are per set
> (features aggregated within a set within a perturbation, then per drug), never pooled across subsystems - which
> is what makes the tca/pyruvate/glycolysis/oxphos gradient readable.
> **Run inventory complete (report 07, section 15):** both phenotypes x both label types x all three modalities.
> The last gap, `hallinger_expr_sig`, is null (0.582, p = 0.25), so **every Hallinger arm is null in both label
> types and all modalities** (n = 20-29). Clarification: in the POOLED expression `sig` run `glycolysis gene`
> does survive BH (q = 0.016) - section 12.3 said otherwise, which was true only of the responsive run. It is
> still not evidence (rand = 0.144; responsiveness confounder 0.230 at p = 1.5e-26; q = 0.106 once restricted),
> so the conclusion holds but the earlier statement was incomplete. Also note `attene_expr_sig*` (702 drugs) and
> `attene_expr_valuesign` (726) differ because only the latter carried the cytotox contrast.

**Mapping tool:** `scripts/benchmarking/map_smiles.R` (rewritten 2026-11-19): largest-fragment stripping (salts),
all pairs kept, NA never matched, unique SMILES canonicalised once, and `--stereo both` = exact pass then a
stereo-insensitive retry of the unmatched, every pair labelled with `match_type` / `ambiguous_stereo` so the cost
of merging enantiomers is measured, not assumed. `data/validation/hallinger/DSSTox_CCD_dump_12092025` (CSV)
gives DTXSID / CASRN -> SMILES / InChIKey offline from a dated release. Clowder downloads:
`https://clowder.edap-cluster.com/api/files/<file-id>/blob` (**HEAD 404s, GET works**).

**Watch out: 5,785 of 34,419 LINCS compounds (564 of 5,494 in HepG2) have neither SMILES nor InChIKey** — mostly
unnamed `BRD-…` codes. They can never be mapped to an external dataset by structure, so the mappable HepG2
denominator is 4,930, not 5,494.

**Currency metabolites** (ATP, ADP, NAD+, NADH, NADPH, O2, H+, CoA) **are MetChange features but not network
nodes** (they are among the 339 dropped before `select_sources`). So an energy-metabolism phenotype can be tested
at the **label** level but cannot be built as a model **query**; queries must use the glycolysis / TCA proxies
(lactate, pyruvate, succinate, citrate, fumarate, acetyl-CoA, glucose). This is why the enrichment test comes
before the ranking test for MMP / respiration phenotypes.

**To do (Attene-Ramos, report 07 section 6):**
1. **Viability columns of `AID_720637_datatable.csv`** (laptop): is the cytotoxicity split possible? Most
   compounds that collapse MMP just kill cells; the informative contrast is **MMP-down without viability loss**
   vs inactive. Without this the test predicts cytotoxicity, not metabolism.
2. **Real n** (HPC): mapped compounds x *active* HepG2 perturbations — not 904.
3. **Enrichment test** (HPC, HepG2 MetChange labels at `<scratch>/results/labels/HEPG2/metchange/rank/`, ~1.5 GB
   per matrix so do not download): per drug, summarise MetChange over glycolytic vs TCA/ETC metabolites
   (currency metabolites included), compare MMP-active vs inactive (Mann-Whitney / AUROC), with the expression
   baseline read the same way, and **no assumed sign** (the compensation hypothesis of report 05, 9.2).
4. Aggregation rule for the 1,379 Attene-Ramos CIDs collapsing onto 904 compounds; run first with `exact` +
   unambiguous pairs, then test whether the 225 ambiguous stereo pairs change anything.
5. **Hallinger:** email Simmons (EPA/CCTE) for the `EPAPLT0033` / `EPAPLT0034` plate maps — the only remaining
   source. The 17 GB MySQL `sample` table is a long shot (tier 1 was never registered) and not worth it on spec.

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
