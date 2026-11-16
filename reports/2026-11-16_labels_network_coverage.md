# Report 04 — The MetChange training labels still had the limma sign bug; network and coverage checks

*2026-11-16. Inputs: `results/labels/A549/compare_corrected_vs_rerun/` (`scripts/compute_labels/compare_label_dirs.py`),
`results/predict_targets/A549/expr_combined_noGlobals_lr1e-3_min15{,_shuff}`,
`results/predict_targets/A549/metabolic_coverage_full{,_topnbins}/` (`check_metabolic_coverage.py`). Metric
definitions and the gap convention: report 01 section 3, report 03 section 1.*

**Summary.**
1. **The MetChange labels every MetChange model read until `d11744b` are the uncorrected limma output.** In 24,512 of
   44,396 A549 perturbations (55%) the home copy has log2FC exactly negated relative to the scratch rerun, and these
   are exactly the perturbations R's level order flags. Every MetChange result in reports 01–03 was trained with those
   signs flipped. The expr and TF runs are not affected.
2. **The network is not why MetChange lags expr:** expr signatures on `combined` give the same gap as on `gene`
   (7.7 vs 7.8 points).
3. **MetChange's own coverage is bimodal:** 57% of active perturbations have no significant feature among the 100 the
   model reads, 35% have all 100. On fully covered perturbations the (sign-buggy) MetChange runs reach a gap of
   7.1–7.3 points in both seeds, against 0.4–2.9 on uncovered ones.

## 1. Corrected vs rerun MetChange labels

`compare_label_dirs.py --a <home>/results/labels/A549/metchange/rank --b <scratch>/results/labels/A549/metchange/rank`.
A = the copy the MetChange launchers read until `d11744b` (described in `STATUS.md` as corrected with
`fix_limma_reference_sign.py`); B = the copy they read now (rerun of the fixed `limma_DE.R`).

| | value |
|---|---|
| perturbations compared (all in both copies) | 44,396 |
| `identical` (max \|ΔlogFC\| = 0) | 19,884 |
| `negated` (r = −1, every feature opposite sign) | 24,512 |
| anything else | 0 |
| q-values, p-values: max \|A − B\| | 0 (all perturbations) |
| significant sets (q < 0.05), top-100 feature sets | identical in all |
| NaN pattern | identical |
| `basal`: max \|A − B\| | = max \|log2FC\| on negated perturbations, 0 elsewhere |

**The negated set is exactly R's flip set:** all 24,512 are `flipped_R = TRUE` in `metchange/rank_signfixed/flipped_R.tsv`
(and `expr_signfixed/flipped_R.tsv`); none of the 19,884 identical ones is. `basal` differing by exactly log2FC is the
signature of `fix_limma_reference_sign.py` (`log2FC → −log2FC`, `basal → basal + log2FC`) having been applied to one
copy and not the other.

**Which copy is wrong.** B was produced by the fixed `limma_DE.R`, which sets the control as reference explicitly, so it
is correct by construction. Even if B were instead the `rank_signfixed` copy that `move_corr_files_mc.sh` placed there,
it would be corrected. The only way B is the buggy one is if it were the raw pre-fix output, which its
`reference_level.txt` rules out: it reads "reference level: ctrl (explicit)" (checked 2026-11-16), a line
only `limma_DE.R` from `85526fe` on writes, after the fix of `dfc8bda`. The rerun finished on 2026-11-12, after both. No script in the repository ever wrote corrected MetChange files
to the home copy (`move_corr_files_mc.sh` works on scratch only). So **A, the home copy, is the uncorrected output**,
and `STATUS.md`'s statement that it was corrected (and that a spot check confirmed it) was wrong.

Two checks that did not settle it, for the record: `basal` is not constant across unflipped perturbations of one time
point (each perturbation has its own control subset), so "basal = control mean" cannot pick the correct copy; and
MetChange signatures of TOP1 inhibitors are too noisy for a biological sign test (topotecan vs irinotecan, neither
flipped, range from −0.8 to +0.8, median 0.2).

**The expression comparison produced no output** (`compare_corrected_vs_rerun/expr/` is empty). Either the launcher
skipped it (same `reference_level.txt` on both sides) or it failed; the log `logs/cmp_labels.out` was not synced. The
expr and TF runs read the scratch copy throughout, so this does not affect them as long as that copy is correct.

### 1.1 What it means for the MetChange results

- **In every MetChange model so far, 55% of perturbations had sign-flipped signatures.** Which ones is decided by
  whether the perturbation name sorts before `DMSO_nan_<time>`, i.e. in practice by the drug name: all doses of
  afatinib, camptothecin, belinostat, … flipped; all doses of topotecan, temozolomide, … not. Drugs sharing a target
  therefore gave the model opposite signatures whenever their names fall on opposite sides of "DMSO". The model has to
  fit target labels to signatures that partly contradict each other.
- **The shuffled controls are unaffected in expectation** (a permuted signature is uninformative either way), so every
  MetChange gap in reports 01–03 is a gap of sign-corrupted signatures over a valid control: most likely an
  **underestimate**. That MetChange still showed 2.7–4.2 points means it read something despite the corruption.
- **The MetChange vs expr comparison of report 03 (sections 9–10) is not valid** until MetChange is retrained on the
  corrected labels. Its central conclusion ("MetChange reads the signature about half as well as expr") may be mostly
  this bug.
- **Explains the similarity-check anomaly** (`STATUS.md` section 7): TOP1 / JAK1 / XIAP shared-target pairs were
  anticorrelated and "the split matches flip status exactly", e.g. camptothecin (flipped) vs topotecan (not).
- **Unaffected:** expr and TF runs, the q-value based activity filter and all coverage counts (q-values identical),
  the training-recipe conclusions (warm-up, schedule, `--no-globals`; these were also seen on expr / TF), and the
  dose/time leak analysis.

## 2. Expression signatures on the `combined` network

`expr_combined_noGlobals_lr1e-3_min15` = `expr_noGlobals_lr1e-3_min15` (seed 0, same folds) with
`networks/combined` (the MetChange network) instead of `networks/gene`. Per-drug gap over own control (points),
Wilcoxon p:

| metric | expr on `gene` | expr on `combined` |
|---|---|---|
| `spec_prob_timedose` | 7.8, 2e-11 | 7.7, 3e-12 |
| `rank_probc` | 8.0, 4e-16 | 8.2, 1e-13 |
| `rank_prob` | 3.3, 5e-15 | 3.0, 2e-12 |
| `rank_prob` per drug, real / control | 0.381 / 0.437 | 0.390 / 0.435 |
| `spec_prob_timedose` per perturbation, real / control | — | 0.319 / 0.500 |

Per fold (`spec_prob_timedose`, per perturbation) real 0.311–0.328, controls 0.485–0.504; real better per drug in 5/5
folds. **The network is not the cause of the MetChange deficit**: the same signatures read through `combined` give the
same signal as through `gene`.

## 3. Metabolic coverage on the full data

`check_metabolic_coverage.py` on the cluster, all 3,565 active perturbations with labels (3,511 with a time × dose
stratum), runs `_min15` seeds 0 and 1. Output `metabolic_coverage_full/`. The fixed `n_sig_topn` split (section 3.2) is
in `metabolic_coverage_full_topnbins/` (rerun locally with the new `--coverage` option, which reuses `coverage.tsv`).

### 3.1 Coverage

| | 5% | 25% | median | 75% | 95% |
|---|---|---|---|---|---|
| DE genes | 265 | 900 | 2,929 | 7,336 | 13,711 |
| DE genes in Human-GEM | 38 | 135 | 436 | 1,086 | 2,005 |
| share in Human-GEM | 0.111 | 0.140 | 0.147 | 0.158 | 0.191 |
| `n_sig_mc` (MetChange features with q < 0.05) | 0 | 0 | 1 | 4,391 | 5,079 |
| `n_sig_topn` (of the 100 read) | 0 | 0 | 0 | 100 | 100 |

- **The laptop subset's picture holds on the full data, more sharply.** In DE terms every active perturbation changes
  GEM genes in proportion (share 14–15%, Spearman `n_de` vs `n_de_gem` 0.995).
- **MetChange's response is bimodal, essentially all or nothing:** 48% of active perturbations have no significant
  MetChange feature at all, 57% none among the 100 the model reads, 35% all 100; only 7% are in between.
- Expression and MetChange activity agree only moderately (Spearman `n_de` vs `n_sig_mc` 0.55, vs `n_sig_topn` 0.47):
  many strongly DE perturbations produce no significant metabolite change.

This answers the question that motivated the check: **most perturbations called active from expression have no
propagated metabolic effect that MetChange calls significant**, and for them the model reads the 100 largest values of
noise. Note that the coverage counts come from q-values, which the sign bug does not touch.

### 3.2 Gap by `n_sig_topn` (fixed bins)

Per-drug gap over own control, `spec_prob_timedose` (points), Wilcoxon p. MetChange rows are the sign-buggy runs (1.1).

| run | 0 of 100 (≈2,010 perts, 1,045 drugs) | 1–99 (≈250, 174) | 100 of 100 (≈1,245, 645) |
|---|---|---|---|
| MetChange s0 | 2.9, 0.07 | 7.4, 0.12 | **7.3, 3e-6** |
| MetChange s1 | 0.4, 0.71 | 0.0, 0.75 | **7.1, 4e-6** |
| expr s0 | 7.3, 7e-9 | 13.3, 2e-5 | 12.1, 5e-8 |
| expr s1 | 11.8, 5e-14 | 16.9, 1e-6 | 10.5, 2e-8 |
| expr on `combined` s0 | 7.1, 3e-9 | 17.9, 2e-5 | 9.9, 1e-8 |

(The quantile bins of `metabolic_coverage_full/` put all of `n_sig_topn` into one bin; `check_metabolic_coverage.py`
now uses fixed 0 / 1–99 / 100 bins for it.)

- **MetChange works where it has something to read.** On fully covered perturbations its gap is 7.1–7.3 points in both
  seeds (p ≈ 4e-6), about where expr was over all perturbations; on uncovered ones it is at or near zero.
- **Expr gains on the uncovered perturbations too** (7–12 points): those perturbations carry drug information in their
  expression that does not reach MetChange's significant metabolite changes. Coverage explains part of the deficit,
  not all of it: on fully covered perturbations expr still leads (10–12 vs 7 points).
- Both readings are of sign-buggy MetChange runs; with correct signs the MetChange numbers can only be re-read after
  retraining. The coverage split itself stays valid.
- The Spearman between per-perturbation gain and coverage is small for all runs (0.01–0.09), as on the subset.

## 4. Conclusions

1. **The MetChange model results of reports 01–03 used labels with the limma sign bug in 55% of perturbations.** They
   must be redone on the corrected (scratch) labels before any MetChange vs expr / TF comparison is drawn. The current
   launchers already read scratch (since `d11744b`).
2. **The network is ruled out** as the cause of the MetChange vs expr gap (expr on `combined` ≈ expr on `gene`).
3. **Coverage is a real, separate limitation:** more than half of the active perturbations have no significant
   MetChange response. An activity filter from MetChange's own output (e.g. `n_sig_topn` ≥ 1, or ≥ 100) is the natural
   next restriction, evaluated on the same perturbations for expr.

## 5. Next steps

1. ~~Confirm provenance~~ (done: scratch is the fixed rerun). Check `logs/cmp_labels.out` for why the expr pair produced
   no output. If the expr pair was skipped because both copies say "sign-fixed", that is fine (the expr runs read
   scratch); otherwise rerun it.
2. **Retrain MetChange `_min15`, seeds 0 and 1, with controls**, on the scratch labels:
   `submit_train_w_shuff_mc_min15_fixedlabels.sh` (keeps the buggy runs as `<run>_signbug{,_shuff}`). Then redo report 03's MetChange vs expr / TF comparison and the coverage split
   (`check_metabolic_coverage.py --coverage metabolic_coverage_full/coverage.tsv`, no labels needed).
3. **Check the other cell lines' MetChange labels** the same way before training on them (`compare_label_dirs.py`), and
   whether any other consumer (similarity checks, `metpath`, `mc_schuhk_assoc`) read the home copy.
4. Then the MetChange activity filter (3), and the representation (top-100 signs vs continuous values).
