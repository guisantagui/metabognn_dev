# Project state — 2026-11-12

Handoff note for continuing in a new session. Everything described here is committed
(HEAD `e03685b`); the working tree is clean.

## Goal and standing constraints

Adapt the Menacho et al. 2026 GNN drug-target framework to LINCS signatures. A549 first,
then the 9 core cell lines (A375, A549, HA1E, HCC515, HEPG2, HT29, MCF7, PC3, VCAP).

- Do not commit; the user commits.
- Keep gene symbols in the prediction pipeline.
- KEGG-derived edges must not be redistributed.
- Uncharacterized drugs are not negatives (no annotation is not evidence of no target).

## 1. Labels

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

**MetChange concatenation.** `mc_res_concat_list_launcher.sh` is now a Slurm array, one task
per cell line (`--array=1-9`, `-c 8`, `POLARS_MAX_THREADS`, skips cell lines that already have
a `*_concat.csv` unless `overwrite=true`). Note: `concat_mc_res.py`'s docstring claims it
computes deltas against matched controls and averages per perturbation; the code only
concatenates batches and drops all-missing metabolites. Only the description is wrong.

## 2. Main finding: the model learns a per-protein prior, not the signature

A549, 5-fold CV by drug, active perturbations (>= 200 DE genes), held-out ranking,
per drug, median percentile (chance 0.5, lower is better):

| run | real | shuffled signatures |
|---|---|---|
| MetChange FNN | **0.409** | 0.418 |
| MetChange logreg | 0.468 | 0.470 |
| expression FNN | 0.467 | not run |
| expression logreg | 0.455 | 0.456 |
| TF FNN | 0.484 | not run |
| TF logreg | 0.442 | 0.441 |
| target-frequency baseline | 0.51 | 0.51 |

- **The linear model learns nothing from the signatures** — real and shuffled are identical for
  all three signature types. Its edge over the frequency baseline comes from the static and
  global features.
- **The MetChange FNN keeps almost all of its performance under shuffling.** Paired per drug:
  55% better, 44% worse, median gain 1.1 percentile points, Wilcoxon p = 0.007. Real wins in
  3 of 5 folds. Small, real, but nearly everything else is prior.
- Supporting evidence: within-drug SD of the percentile is 0.05 (real) / 0.07 (shuffled), so a
  drug's score barely changes across its own perturbations; Spearman real vs shuffled per drug
  is 0.75.
- **Target classes are prior, not signature:** EGFR 0.03/0.04, HSP90AA1 0.04/0.02, TOP2A
  0.14/0.15, ATP1A1 0.21/0.19, NR3C1 0.54/0.39 (real/shuffled). An earlier reading that
  MetChange "adds" HSP90/TOP2A/ATP1A1/glucocorticoids was wrong.
- Where the real signatures do help (small classes, noisy): ROCK1, PRKACA, MAPK11, CALM1, ACHE,
  calcium channels.
- **Consequence:** the bar is the shuffled FNN (~0.42 per drug), not 0.5 and not the linear
  model. Cross-signature comparison (MetChange vs expression vs TF) is not yet interpretable,
  because the three use different networks (combined / gene / tf) and only MetChange has an
  FNN shuffle control.

Results live in `results/predict_targets/A549/{expr,tf_active200,metchange_rank}` and the
`*_shuff` / `*_linOnly_shuff` folders; the MetChange real-vs-shuffled comparison is in
`mc_real_vs_suff/`.

## 3. Code added this session (commit `e03685b`, plus earlier ones)

In `scripts/predict_targets/metabolic_target_prediction.py`:

- **Prior-free metrics**, written next to the existing ranks in `metrics.json`,
  `cv_metrics.json` and the rank tables. Both give 0.5 for a model that ignores the signature;
  verified on synthetic data (pure prior 0.49/0.50, planted signal ~0.00).
  - `spec_prob_heldout` / `spec_logreg_heldout` (+ `n_spec_heldout`): the rank of the score a
    perturbation gives its own drug among the scores other drugs' perturbations give that same
    drug, within the same activity class. Cancels the per-drug prior.
  - `rank_probc_heldout` / `rank_logregc_heldout`: the held-out ranking after standardising each
    candidate drug's scores over the test perturbations.
  - `test_scores_heldout.npz` per fold: every test perturbation's score for every held-out drug,
    so further prior-free metrics need no rescoring.
- `--no-static` / `--no-globals`: ablations dropping the per-protein static features (activity
  and its propagation, degree, strength) and the per-perturbation globals (dose, time,
  signature size).
- `--score-active-only`: scores only active test perturbations (~6x faster; `test_per_*` then
  equals `test_per_*_active`). Added to the six A549 launchers using `--min-de-genes 200`.

New tooling in the same folder:

- `compare_runs.py`: compares runs side by side (summary, per fold, paired per-drug Wilcoxon,
  ECDF plots). Reads unfinished runs fold by fold; defaults to active perturbations, so runs
  with and without `--score-active-only` are comparable.
- `launch_compare_runs.sh`, `submit_with_shuffle_control.sh` (submits real + shuffle in
  parallel plus a dependent comparison job; `--real-job <id>` when the real run is already
  queued).
- Launchers: `*_shuff` for tf/expr/mc, `tf_lr1e-3`, `mc_noStatic`, `mc_noStatic_shuff`.
  `launch_predict_targets_A549.sh` was renamed to `launch_predict_targets_A549_mc.sh`.

**Not yet exercised end to end:** `--score-active-only`, `--no-static`, and the new metrics.
There is no torch or network data on the laptop, so only unit-level tests were possible.

## 4. Running now

Four MetChange runs with the updated script: `mc`, `mc_shuff`, `mc_noStatic`,
`mc_noStatic_shuff`. They overwrite `metchange_rank` and `metchange_rank_shuff`.

When they finish, compare all four in one call and judge on the prior-free metrics:

```bash
sbatch launch_compare_runs.sh \
    --run real=<A549>/metchange_rank \
    --run real_noStatic=<A549>/metchange_rank_noStatic \
    --run shuffled_noStatic=<A549>/metchange_rank_noStatic_shuff \
    --run shuffled=<A549>/metchange_rank_shuff \
    --control shuffled --out-dir <A549>/mc_static_ablation
```

- Signature signal: `spec_prob_heldout` clearly below 0.5 in the real runs, ~0.5 in the shuffled.
- Share of performance that is prior: the gap between `real` and `real_noStatic` on
  `rank_prob_heldout`.
- If `spec` is below 0.5 in the *shuffled* run, dose/time/signature size leak drug identity;
  the ablation for that is `--no-globals`.

## 5. Next step: make the model use the signature

Agreed plan, not yet implemented.

**(1) Matched negatives.** For each positive row (perturbation p, protein g), add negatives
(q, g): the same protein in other perturbations of the same chunk whose drug does not annotate
g. Within such a group every per-protein feature is identical, so protein identity carries no
information and only the signature can separate the rows. Implementation notes:

- The rows already exist. `iter_chunks` builds `X` of shape `(P * n_prot, F)`, perturbation-major,
  so the matched row is `j' * n_prot + g`. The change is confined to `sample_rows`.
- Filtering on `y[j' * n_prot + g] == 0` automatically excludes other doses of the same drug
  (same label vector) and unknown-direction targets (class 3).
- Chunks currently follow the perturbation table order, which groups doses of one drug together;
  the perturbation order must be shuffled before chunking or a chunk may have no valid partners.
- Keep random negatives too (they teach the useful prior) and give the matched ones their own
  sample weight, or they are swamped: per chunk roughly 100-300 positives, 6400 random
  negatives, ~400-1200 matched at K=4.
- Stronger variant, if needed: keep each group together and put a softmax or pairwise ranking
  loss on it, so protein identity cancels exactly rather than on average. Batches then become
  `(n_groups, K+1, F)`.

**(2) Select epochs on validation specificity** instead of validation AUPRC, which rewards the
prior. Should also calm the training instability.

**(3) Standardise each dynamic feature per protein** over training perturbations ("how unusual
is the signal arriving at g here"), removing the hub effect that `--no-static` cannot reach.

**(4) Prior + residual:** a static-only model learns the prior as a fixed offset, the signature
network is trained on top. Lets you report prior+residual and residual alone.

Across cell types these keep baseline information useful — matching within a cell type stops
baseline from being an identity shortcut while still letting the model learn that a signature
around g only means something where g is expressed.

## 6. Open issues

- **FNN training is unstable.** Validation AUPRC oscillates (0.03-0.22); some ensemble members
  keep their epoch-1 or epoch-2 weights; expression and TF collapse in fold 3 (median 0.55-0.57,
  worse than chance). Probable cause: lr 0.005 with Adam on ~2,100 polynomial features.
  `launch_predict_targets_A549_tf_lr1e-3.sh` (lr 0.001, patience 10) is written but not run.
  If it helps, the real runs need rerunning, and their shuffle controls with them.
- **Expression and TF FNN shuffle controls are missing** (`expr_shuff`, `tf_shuff` launchers
  exist). Without them the comparison across signature types says nothing.
- **MetChange similarity-check anomaly, unresolved.** TOP1 / JAK1 / XIAP shared-target pairs are
  anticorrelated (mean r -0.24 / -0.23 / -0.35) and the split matches flip status exactly,
  although ATP1A1 drugs are mixed-flip yet r = +0.92. Suggested checks: confirm `labels_dir`;
  grep `overlap_drug_pairs.tsv` for camptothecin-topotecan in both runs; rerun the fixed
  `limma_DE.R --perts` on those perturbations and run `check_limma_sign_fix.py compare`.
- **Model does not use edge signs.** `propagation_operator` gates on weight x activity only; the
  two channels propagate independently. Signed propagation is a possible future option.
- **Signatures are discretized** (`np.sign`), with `--ignore-qvalue --top-n 100` in the
  launchers. The similarity check, by contrast, uses continuous Spearman. A
  `--signature-weights {sign,value,t}` option was proposed, not implemented.
