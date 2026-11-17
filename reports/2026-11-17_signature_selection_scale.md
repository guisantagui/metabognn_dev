# Report 05 — Signature selection: q-values unused, and MetChange's top 100 is chosen by feature scale

*2026-11-17. No new model runs (scratch maintenance). Label analysis on the laptop copies
(`results/labels/A549/{metchange/rank,expr}/`, the old partial set of 2 Nov: 7,636 perturbations, all, not only
active ones; signs of the MetChange copy are the buggy ones, magnitudes and p-values are unaffected). Launchers
for the q-thresholded signatures written (section 3).*

## 1. How a signature is built now

`load_signatures` (`metabolic_target_prediction.py`) with `--ignore-qvalue --top-n 100`, for MetChange, expr and TF
alike: per perturbation, the 100 features with the largest |log2FC| (after mapping to network nodes), then `np.sign`.
`select_sources` turns each into a ±1/n weight on its node(s). So:

- **q-values are never used** for the signature (only the expression DE counts, for the activity filter).
- **log2FC magnitudes are used at exactly one point: choosing the 100 features.** Nothing standardises them, per
  feature or across perturbations.
- `--ignore-qvalue --top-n 100` came in on 9 Nov (`5a28a0f`) together with the activity filter, with no recorded
  reason. The first launcher (2 Nov) thresholded at q < 0.05. The published setting is significant DE TFs, max 100,
  so the TF runs are not quite the published positive control either.

For MetChange, `log2FC` is limma's mean difference of MetChange `rank` scores (not a log ratio): values in the tens
to thousands, with very different spread per metabolite. The MetChange paper z-scores perturbation scores against the
controls (`(pert − mean(ctrl)) / sd(ctrl)`) precisely so that metabolites can be compared. Our pipeline skips this,
which did not matter for signs but does matter for the top-100 selection.

## 2. What the selection picks

Top-100 sets per perturbation under three rankings: the current |log2FC|; |log2FC| divided by the feature's SD across
perturbations (a z-like proxy; the controls' SD is not in the label files); and −log10 p (≈ |t|, scale-free). Mean
Jaccard of the top-100 sets of two perturbations of different drugs (3,000 random pairs).

| | MetChange: current | MetChange: /SD | MetChange: p | expr: current | expr: /SD | expr: p |
|---|---|---|---|---|---|---|
| features ever in a top 100 | 2,249 of 4,086 | 3,869 | 3,818 | 5,790 of 19,293 | 16,877 | 18,191 |
| features in the top 100 of > 50% of perturbations | **33** | 0 | 0 | 0 | 0 | 0 |
| mean Jaccard, different drugs | **0.218** | 0.028 | 0.024 | 0.048 | 0.007 | 0.004 |
| Spearman (feature SD, how often selected) | **+0.89** | −0.57 | −0.25 | +0.76 | −0.27 | −0.07 |

(expr: a random 1,000 of the 7,636 perturbations.)

- **MetChange's top 100 is largely a fixed set of high-variance metabolites.** How often a metabolite is selected is
  almost entirely its scale (Spearman 0.89); 33 metabolites are in more than half of all signatures, and two unrelated
  drugs share on average 22% of their selected features. Ranked scale-free, they share 2–3%.
- What then distinguishes MetChange signatures is mostly the **signs** of that common set, not *which* metabolites
  changed. Much of the information about which metabolites a drug moves is discarded at selection.
- **Expr has the same bias, much weaker**: selection still tracks gene SD (0.76), but no gene is in most signatures and
  different drugs share 5%. log2FC is a genuine log ratio there, so scales are comparable across genes.
- Current vs p-based top 100 for MetChange: median overlap 0 (95th percentile 0.22). The two rankings pick almost
  disjoint features, partly because most perturbations have no significant feature (report 04, section 3), so their
  p-ranking is noise too.

Caveat: laptop subset, all perturbations (active and inactive), buggy MetChange signs (irrelevant here: only
magnitudes and p-values are used).

## 3. What follows

1. **The `_qsig` runs (written today) are scale-free by construction:** significant features (q < 0.05), random cap
   of 100, so the selection no longer depends on log2FC magnitude. They also drop perturbations with < 3 significant
   features, which is the MetChange coverage filter. Launchers
   `launch_predict_targets_A549_{mc,tf,expr}_noGlobals_lr1e-3_min15_qsig{,_shuff}.sh`, wrapper
   `submit_train_w_shuff_noGlobals_lr1e-3_min15_qsig.sh` (MetChange on the corrected scratch labels).
2. If `_qsig` loses too many MetChange perturbations: keep top-100, but rank by p (or by the limma t-statistic, which
   `limma_DE.R` does not write yet) instead of |log2FC|: scale-free without dropping perturbations. Small code change
   (`--rank-by {value,pvalue}`), not done.
3. Z-scoring against the controls (MetChange paper) is not used: many controls have no variability in many
   features, which then drop out. **Instead: standardise each feature across perturbations** (implemented,
   section 4).
4. Read every MetChange result so far, including the corrected retrain running now, with this in mind: its
   signatures are mostly signs on a common set of high-variance metabolites.

## 4. `--signature-scale` (implemented)

`metabolic_target_prediction.py --signature-scale {none,sd,mad}` (default `none`, earlier runs unchanged): each
feature of `--signatures` divided by its spread across all perturbations of the label file (active and inactive; no
target labels involved, so no leakage), **not centred** (centring would turn the sign into "above the average
perturbation" instead of "vs control"). The spread is floored at the `--signature-scale-floor` quantile (0.05) of all
features' spreads, so near-constant features are not inflated. It only changes which features `--top-n` selects; signs
are unchanged.

Check on the laptop MetChange labels (function extracted from the script; no torch there):

| | none | sd | mad |
|---|---|---|---|
| mean Jaccard of top-100 sets, different drugs | 0.218 | 0.029 | 0.027 |
| features in > 50% of top-100 sets | 33 | 0 | 0 |
| features floored | — | 201 of 4,086 | 203 |
| signs identical to unscaled | — | yes | yes |

Run: `bash submit_train_w_shuff_mc_noGlobals_lr1e-3_min15_sdscale.sh` (MetChange `_min15` + `--signature-scale sd`,
seed 0, corrected labels, with control). Compare with the corrected `_min15` retrain (same seed and folds) and with
`_qsig`. Expr could get the same option, but its scale bias is much smaller (section 2).
