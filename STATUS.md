# Project state — 2026-11-14

Handoff note for continuing in a new session. Everything described here is committed up to
HEAD `85ce09b` (reports 01–02 included), except the files listed in section 8.

**Read first:** `reports/2026-11-14_lr_noglobals_matched.md` (report 02 — current best recipe
and the main result), then `reports/2026-11-13_A549_matched_negatives.md` (report 01 — metric
definitions, the dose/time leak, how to read a comparison). `scripts/predict_targets/README.md`
documents every launcher and how to submit runs with their controls.

## Goal and standing constraints

Adapt the Menacho et al. 2026 GNN drug-target framework to LINCS signatures. A549 first,
then the 9 core cell lines (A375, A549, HA1E, HCC515, HEPG2, HT29, MCF7, PC3, VCAP).

- Do not commit; the user commits.
- Keep gene symbols in the prediction pipeline.
- KEGG-derived edges must not be redistributed.
- Uncharacterized drugs are not negatives (no annotation is not evidence of no target).
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

**Current best recipe** (report 02): `--no-globals --lr 0.001 --patience 10 --ffn-epochs 60`,
unmatched, epochs selected on AUPRC. Run `metchange_rank_noGlobals_lr1e-3`:

| metric (per perturbation / per drug) | real | shuffled | per-drug p |
|---|---|---|---|
| `spec_prob_timedose` (leak-free) | 0.414 / 0.444 | 0.500 / 0.500 | 4.5e-5 |
| `rank_prob` | 0.262 / 0.401 | 0.274 / 0.410 | 9.4e-4 |

Real beats its control in 5/5 folds; best `rank_prob` of any run so far; 23% of perturbations
rank their own drug in the top 10%. **Training stability was the bottleneck, not the sampling.**

**Matched negatives are superseded.** With dose/time inputs gone and training to convergence,
matched sampling is worse on both metrics (report 02, section 4). What looked like
matched-negative specificity in report 01 was mostly the dose/time shortcut plus `spec`-based
epoch selection. Version (b) (grouped specificity loss) is demoted accordingly.

**Still true:** per-drug plain ranking remains mostly prior (0.401 vs 0.410), so a query will
return broadly similar drugs whatever is asked, with the signature re-ordering within that.

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
- `--fold i` with `--n-folds 5` for array jobs; the last task pools (`--aggregate-only` to pool
  by hand).

Tooling:

- `compare_runs.py` — summary / per fold / paired per-drug Wilcoxon / ECDFs, plus `spec`
  **stratified** by time and by time × dose, with mid-rank percentiles (unbiased for small
  strata; the plain formula averages 0.58 at n = 5 — report 01, section 5.4).
- `check_dose_time_leak.py` — `spec` of the model vs `spec` of a dose/time-only predictor, and
  the score variance explained by dose/time groups.
- `submit_with_shuffle_control.sh` — submits a variant's real run, its shuffled control and a
  dependent comparison; takes several variants in one call, validates names first, `--real-job`
  for an already-queued real run. Needs bash ≥ 4 (cluster yes, macOS no).
- Wrappers: `submit_train_w_shuff_mc_matched.sh`,
  `submit_train_w_shuff_mc_noGlobals_matched_auprc.sh`,
  `submit_train_w_shuff_mc_matched_noGlobals.sh`,
  `submit_train_w_shuff_mc_noGlobals_lr1e-3.sh`,
  `submit_train_w_shuff_mc_noGlobals_lr1e-3_ep100_seed1.sh` (6.1),
  `submit_train_w_shuff_tf_expr_noGlobals_lr1e-3.sh` (6.4).
- Variants registered in the submit script: `tf`, `expr`, `mc`, `mc_matched`,
  `mc_matched_auprc`, `mc_matched_noGlobals`, `mc_noGlobals`, `mc_noGlobals_lr1e-3`,
  `mc_matched_noGlobals_lr1e-3`, `mc_noGlobals_lr1e-3_ep100`, `mc_noGlobals_lr1e-3_seed1`,
  `tf_noGlobals_lr1e-3`, `expr_noGlobals_lr1e-3`. `mc_noStatic` is not registered (launch both
  by hand).
- The `_lr1e-3` launchers are **array jobs** (`--array=0-4`; 4 h per task, 6 h for `_ep100` and
  the `tf` / `expr` ones). `sbatch --test-only` prints the estimated start time without
  submitting.

## 5. Runs in `results/predict_targets/A549/`

Analysed (in reports): `metchange_rank{,_shuff}`, `metchange_rank_matched{,_shuff}`,
`metchange_rank_noStatic{,_shuff}`, `metchange_rank_noGlobals{,_shuff}`,
`metchange_rank_matched_auprc`, `metchange_rank_noGlobals_lr1e-3{,_shuff}`,
`metchange_rank_matched_noGlobals_lr1e-3{,_shuff}`. Comparisons live in each real run's
`compare_real_vs_shuffled/`; cross-run comparisons in `mc_matched_vs_unmatched/`,
`mc_real_vs_suff/`, `stratified_spec/`.

Not on the laptop: `metchange_rank_matched_auprc_shuff` (its comparison ran on the cluster),
and `metchange_rank_matched_noGlobals{,_shuff}` (plain, lr 0.005) — launched 2026-11-13, not
yet synced or analysed. It is now of low value: report 02 supersedes it.

Launchers written 2026-11-14, **not yet submitted**: `metchange_rank_noGlobals_lr1e-3_ep100`,
`metchange_rank_noGlobals_lr1e-3_seed1`, `tf_active200_noGlobals_lr1e-3`,
`expr_noGlobals_lr1e-3`, each with `_shuff` (sections 6.1, 6.4).

Pre-dating the current metrics, cannot be re-checked (no `test_scores_heldout.npz`, no FNN
shuffle control): `expr`, `tf`, `tf_active200`, the `*_linOnly_shuff` runs.

## 6. Next steps, in priority order

**6.1 Confirm the new result is solid** (cheap, GPU) — launchers ready, submit with
`bash submit_train_w_shuff_mc_noGlobals_lr1e-3_ep100_seed1.sh`:
1. `--ffn-epochs 100` on `mc_noGlobals_lr1e-3` + control — one member kept epoch 53 of 60.
   Check the kept epochs: if none approaches 100 and the metrics match, the 60 cap is fine.
2. A second seed (`--seed 1`) of the same pair, to size run-to-run variability before building
   on it. The seed also redraws the drug folds, so compare per-perturbation / per-drug metrics
   and the real − shuffled gap with seed 0, not fold by fold.

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

**6.3 Size the data ceiling** (no GPU, ~a day):
1. PCA of the A549 expression matrix: how many components explain 95% of variance across
   perturbations. That is the budget of distinguishable perturbation states.
2. Fraction of the GPR genes MetChange uses that are landmarks vs CycleGAN-inferred. If mostly
   inferred, the metabolic labels are a smooth function of the 978 landmarks.

**6.4 Then** freeze the recipe and port it to `tf` and `expr` with controls (first
interpretable cross-signature comparison) — launchers ready (`tf_noGlobals_lr1e-3`,
`expr_noGlobals_lr1e-3`), submit with `bash submit_train_w_shuff_tf_expr_noGlobals_lr1e-3.sh`.
Can go in parallel with 6.1, but read them against 6.1's seed-to-seed spread. Only after that
revisit version (b), the prior + residual model, or per-protein standardisation of the dynamic
features.

## 7. Open issues

- **Per-drug plain ranking is still prior-dominated** (section 3). This is the main scientific
  gap for the tool, not a bug.
- **Expression and TF FNN shuffle controls are missing** for the old recipe (`expr_shuff`,
  `tf_shuff` launchers exist). Superseded by the `*_noGlobals_lr1e-3` pairs (6.4), which come
  with their controls; without those, cross-signature comparison says nothing.
- **MetChange similarity-check anomaly, unresolved.** TOP1 / JAK1 / XIAP shared-target pairs are
  anticorrelated (mean r -0.24 / -0.23 / -0.35) and the split matches flip status exactly,
  although ATP1A1 drugs are mixed-flip yet r = +0.92. Checks: confirm `labels_dir`; grep
  `overlap_drug_pairs.tsv` for camptothecin-topotecan in both runs; rerun the fixed
  `limma_DE.R --perts` on those perturbations and run `check_limma_sign_fix.py compare`.
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

- New launchers in `scripts/predict_targets/`: `launch_predict_targets_A549_{mc_noGlobals_lr1e-3_ep100,
  mc_noGlobals_lr1e-3_seed1,tf_noGlobals_lr1e-3,expr_noGlobals_lr1e-3}{,_shuff}.sh` and the two
  wrappers above.
- `submit_with_shuffle_control.sh` (four variants registered) and `README.md` (variants table).
- This `STATUS.md` update.

## 9. Working notes

- A safety classifier has repeatedly truncated long prose answers in this project (three times
  on 2026-11-13). Code, tables, numbers and tool output were never affected; the pattern was
  broad written assessments. Keeping write-ups in `reports/` and chat answers short avoids it.
  Starting a new session helps because it drops accumulated context — hence this handoff.
- Laptop has `python3` (no `python`), no torch, bash 3.2 (no `declare -A`), and no network
  data, so only unit-level tests and result analysis are possible there.
