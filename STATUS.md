# Project state — 2026-11-14

Handoff note for continuing in a new session. Everything described here is committed up to
HEAD `a53a99c`, except the files listed in section 8.

**Read first:** `reports/2026-11-14_seed_tf_expr.md` (report 03 — second seed, MetChange vs expr
vs TF, why FNNs keep epoch 1), then `reports/2026-11-14_lr_noglobals_matched.md` (report 02 —
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
unmatched, epochs selected on AUPRC. **Training stability was the bottleneck, not the sampling.**

**Gap convention** (report 03, section 1): per-drug gap = `median_diff` of `paired.tsv`, the
median of the per-drug real − control differences, in percentile points; it is what the
Wilcoxon p tests. Report 02's "5.6 points" was a difference of medians (paired: 4.2).

Per-drug gap (points), p — real vs own control (report 03):

| metric | MetChange seed 0 | MetChange seed 1 | expr | TF (positive control) |
|---|---|---|---|---|
| `spec_prob_timedose` (leak-free) | 4.2, 4.5e-5 | 2.9, 0.047 | 4.1, 1.7e-3 | 6.8, 1.2e-8 |
| `rank_probc` | 6.8, 5.7e-10 | 3.5, 1.5e-3 | 6.4, 4.0e-6 | 9.5, 4.2e-15 |
| `rank_prob` | 1.1, 9.4e-4 | 0.5, 0.09 | 1.1, 0.19 | 3.3, 5.0e-8 |
| FNN members keeping epoch ≤ 2 | 2/15 | 1/15 | 8/15 | 7/15 |

- MetChange reads the signature in both seeds (5/5 folds each), at roughly half to two-thirds
  of the TF level; the seed-to-seed spread (2.9–4.2) is as large as the differences between
  signature types.
- The epoch cap does not bind: cap 100 gives byte-identical scores (training is deterministic;
  no member reached 60). `StepLR` halves the lr every 10 epochs, so by epoch 60 it is 1/64.
- **New bottleneck: early stopping.** Validation AUPRC peaks at epoch 1, dips, and recovers
  after the first lr halvings; patience 10 stops before the recovery. About half of the TF /
  expr members keep epoch 1–2; expr folds 0–1 collapsed (fold 1 `rank_prob` 0.689) while its
  folds 2 and 4 are as strong as TF. Fix: `--min-epochs` (section 4), runs in 6.1.

**Matched negatives are superseded.** With dose/time inputs gone and training to convergence,
matched sampling is worse on both metrics (report 02, section 4). What looked like
matched-negative specificity in report 01 was mostly the dose/time shortcut plus `spec`-based
epoch selection. Version (b) (grouped specificity loss) is demoted accordingly.

**Still true:** per-drug plain ranking remains mostly prior for MetChange (0.401 vs 0.410), so a
query will return broadly similar drugs whatever is asked, with the signature re-ordering within
that. TF is the only signature where the signature clearly moves plain ranking (control 0.372 vs
real 0.247 per perturbation; different network, so compare gaps, not absolute values).

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
  `submit_train_w_shuff_mc_noGlobals_lr1e-3_ep100_seed1.sh` (done, report 03),
  `submit_train_w_shuff_tf_expr_noGlobals_lr1e-3.sh` (done, report 03),
  `submit_train_w_shuff_noGlobals_lr1e-3_min15.sh` (6.1).
- Variants registered in the submit script: `tf`, `expr`, `mc`, `mc_matched`,
  `mc_matched_auprc`, `mc_matched_noGlobals`, `mc_noGlobals`, `mc_noGlobals_lr1e-3`,
  `mc_matched_noGlobals_lr1e-3`, `mc_noGlobals_lr1e-3_ep100`, `mc_noGlobals_lr1e-3_seed1`,
  `tf_noGlobals_lr1e-3`, `expr_noGlobals_lr1e-3`, `mc_noGlobals_lr1e-3_min15`,
  `tf_noGlobals_lr1e-3_min15`, `expr_noGlobals_lr1e-3_min15`. `mc_noStatic` is not registered
  (launch both by hand).
- The `_lr1e-3` launchers are **array jobs** (`--array=0-4`; 4 h per task, 6 h for `_ep100` and
  the `tf` / `expr` ones, 8 h for `_min15`). `sbatch --test-only` prints the estimated start time without
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

Launchers written 2026-11-14, **not yet submitted**: `metchange_rank_noGlobals_lr1e-3_min15`,
`expr_noGlobals_lr1e-3_min15`, `tf_active200_noGlobals_lr1e-3_min15`, each with `_shuff` (6.1).

Pre-dating the current metrics, cannot be re-checked (no `test_scores_heldout.npz`, no FNN
shuffle control): `expr`, `tf`, `tf_active200`, the `*_linOnly_shuff` runs.

## 6. Next steps, in priority order

**6.1 Fix early stopping** (GPU) — launchers ready, submit with
`bash submit_train_w_shuff_noGlobals_lr1e-3_min15.sh` (MetChange, expr, TF + controls, with
`--min-epochs 15 --patience 15 --ffn-epochs 80`). Check: no member keeps epoch 1–2; expr folds
0–1 recover; whether MetChange's gap to TF narrows or TF's grows. If expr recovers, MetChange vs
expr (the expression it is computed from) is the comparison that matters for the paper. Then a
second seed of the winning recipe before building on it (report 03, section 3: the seed spread
is about as large as the differences between signature types).

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

**6.4 Then** freeze the recipe, and only after that revisit version (b), the prior + residual
model, or per-protein standardisation of the dynamic features. (Done on 2026-11-14: cap check —
does not bind; second seed — replicates, smaller; port to `tf` / `expr` — report 03.)

## 7. Open issues

- **Per-drug plain ranking is still prior-dominated** (section 3). This is the main scientific
  gap for the tool, not a bug.
- **Early stopping keeps epoch 1 in many FNN members** (report 03, section 5); `--min-epochs`
  added, runs pending (6.1). The old-recipe `expr` / `tf` runs still lack FNN controls, but are
  superseded by the `*_noGlobals_lr1e-3` pairs.
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

- `reports/2026-11-14_seed_tf_expr.md` (report 03, new); correction note in report 02 (its
  "points" were differences of medians).
- `metabolic_target_prediction.py`: `--min-epochs`.
- New launchers `launch_predict_targets_A549_{mc,tf,expr}_noGlobals_lr1e-3_min15{,_shuff}.sh` and
  wrapper `submit_train_w_shuff_noGlobals_lr1e-3_min15.sh`; `submit_with_shuffle_control.sh`
  (three variants registered) and `README.md` (variants table, `--min-epochs`).
- This `STATUS.md` update.

## 9. Working notes

- A safety classifier has repeatedly truncated long prose answers in this project (three times
  on 2026-11-13). Code, tables, numbers and tool output were never affected; the pattern was
  broad written assessments. Keeping write-ups in `reports/` and chat answers short avoids it.
  Starting a new session helps because it drops accumulated context — hence this handoff.
- Laptop has `python3` (no `python`), no torch, bash 3.2 (no `declare -A`), and no network
  data, so only unit-level tests and result analysis are possible there.
