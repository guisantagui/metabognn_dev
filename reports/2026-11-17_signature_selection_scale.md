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

## 5. MetChange retrained on the corrected labels (`_min15`, seeds 0 and 1)

Runs `metchange_rank_noGlobals_lr1e-3_min15{,_seed1}{,_shuff}`, current recipe, scratch labels (rerun of the fixed
`limma_DE.R`; `config_fold0.json` checked), same seeds and folds as the buggy runs of report 03 (now
`..._signbug` on scratch; the laptop copy of the seed-0 buggy pair was overwritten by this download, its numbers are
in reports 03–04). Unscaled `--ignore-qvalue --top-n 100` selection, i.e. still with the scale problem of section 2.

Per-drug gap over own control (points), Wilcoxon p:

| metric | MetChange s0: buggy → corrected | MetChange s1: buggy → corrected | expr s0 / s1 | TF s0 / s1 |
|---|---|---|---|---|
| `spec_prob_timedose` | 4.2 → **8.3** (6e-12) | 2.7 → **7.7** (4e-11) | 7.8 / 10.4 | 9.0 / 9.3 |
| `rank_probc` | 6.4 → **9.3** | 4.4 → **7.7** | 8.0 / 11.5 | 10.4 / 10.3 |
| `rank_prob` | 1.3 → **2.5** (1e-11) | 0.7 → **2.5** (5e-9) | 3.3 / 3.6 | 3.3 / 3.3 |
| `rank_prob` per drug, real / control | 0.381 / 0.415 | 0.376 / 0.414 | 0.381 / 0.437, 0.390 / 0.441 | 0.397 / 0.445, 0.383 / 0.452 |

Real beats control per drug in 5/5 folds in both seeds (fold medians: real 0.39–0.43, control 0.47–0.53).

**Paired by drug between signatures** (`spec_prob_timedose`, as in report 03, 10.1):

| comparison | seed | median gap A / B | paired median difference | drugs A > B | p |
|---|---|---|---|---|---|
| expr vs MetChange | 0 | 7.8 / 8.3 | +1.0 | 51% | 0.79 |
| expr vs MetChange | 1 | 10.4 / 7.7 | +2.2 | 53% | 0.085 |
| TF vs MetChange | 0 | 9.0 / 8.3 | +2.4 | 53% | 0.21 |
| TF vs MetChange | 1 | 9.3 / 7.7 | +2.5 | 53% | 0.12 |
| MetChange corrected vs buggy | 1 | 7.7 / 2.7 | +3.0 | 55% | 2.5e-5 |

(Report 03's "expr and TF beat MetChange, p = 0.023 / 7e-5 and 9e-4 / 7e-5" is withdrawn: it was the sign bug.)

**Coverage split, corrected runs** (`metabolic_coverage_fixedlabels/`, `--coverage` reusing report 04's counts):

| `n_sig_topn` | perts / drugs | MetChange s0 | MetChange s1 | expr s0 | expr s1 |
|---|---|---|---|---|---|
| 0 of 100 | 2,010 / 1,045 | 5.4 (5e-4) | 5.6 (3e-4) | 7.3 | 11.8 |
| 1–99 | 250 / 174 | 8.4 (1e-3) | 13.5 (1e-4) | 13.3 | 16.9 |
| 100 of 100 | 1,245 / 645 | **14.6** (<1e-6) | **14.7** (<1e-6) | 12.1 | 10.5 |

Per-perturbation Spearman of MetChange's gain with `n_sig_topn`: 0.10 / 0.11 (expr: 0.01 / 0.02).

### 5.1 Reading

1. **The MetChange deficit was the sign bug.** On corrected labels MetChange reads the signature as well as expr and
   TF: gap 7.7–8.3 points, not significantly different from either, paired by drug, in both seeds. Plain per-drug
   ranking now moves too (2.5 points, p < 1e-8; buggy: 0.7–1.3), close to expr / TF (3.3–3.6).
2. **Where MetChange has a response it outperforms expr**: 14.6–14.7 points on fully covered perturbations vs 10.5–12.1
   for expr on the same ones. Its weaker spot is the 57% of active perturbations with no significant MetChange
   response (5.4–5.6 vs 7.3–11.8), which pulls its overall gap down to the expr level. Unlike before, its gain now
   tracks coverage (Spearman 0.10 vs ~0.01 for expr).
3. Even uncovered perturbations carry some signal (5.4–5.6 points, p < 1e-3): the top-100 values are not pure noise
   when no feature passes q < 0.05.
4. These runs still use the scale-biased selection (section 2). `_sdscale` and `_qsig` test whether a scale-free
   selection adds to this.

### 5.2 MetChange `_qsig`: worse than top-100

`metchange_rank_noGlobals_lr1e-3_min15_qsig{,_shuff}` (q < 0.05, random cap 100, corrected labels, seed 0). The
significance filter keeps 1,791 of 3,970 perturbations (758 of 1,358 drugs). Expr `_qsig`: section 5.3. TF `_qsig`:
still running.

Against its own control: `spec_prob_timedose` 6.1 points (p = 1e-3, 56% of drugs better, 5/5 folds), `rank_probc`
7.7 (4e-6), `rank_prob` 2.4 (9e-7; per drug 0.409 vs 0.454).

Against the corrected `_min15` runs on the **same drugs / the same perturbations** (each drug held out once per run, so
per-drug gaps pair across runs even though the folds differ):

| reference | restricted to | drugs | reference gap | `_qsig` gap | paired difference | p |
|---|---|---|---|---|---|---|
| MetChange `_min15` s0 | same drugs | 758 | 11.3 | 6.1 | −4.1 | 0.01 |
| MetChange `_min15` s0 | same perturbations | 735 | 14.2 | 5.9 | −6.7 | 6e-4 |
| MetChange `_min15` s1 | same drugs | 758 | 10.4 | 6.1 | −3.4 | 0.012 |
| MetChange `_min15` s1 | same perturbations | 736 | 11.9 | 5.9 | −5.3 | 3e-3 |
| expr `_min15` s0 / s1 | same perturbations | 735 / 736 | 11.8 / 10.5 | 5.9 | −5.6 / −4.6 | 6e-3 / 0.02 |

**Reading.** On the perturbations it keeps, `_qsig` reads the signature about half as well as the current top-100
selection. Two changes are confounded: (a) half the training perturbations are gone, and the dropped (uncovered) ones
do carry signal (section 5, 5.4–5.6 points); (b) covered perturbations have thousands of significant features (median
`n_sig_mc` ≈ 4,400 among covered ones), and `_qsig` keeps a random 100 of them instead of the largest changes. So
magnitude carries information despite the scale bias; the significance threshold is not the way to remove it.
`_sdscale` (top 100 by scaled magnitude, all perturbations) is the cleaner test. The published TF setting (`_qsig` for
TF) is still worth having as the faithful positive control.

### 5.3 Expr `_qsig`: the signal mostly disappears

`expr_noGlobals_lr1e-3_min15_qsig{,_shuff}` (seed 0). Expr activity is defined from expr DE, so **no perturbation is
dropped** (1,358 drugs, all 3,970 perturbations). Every perturbation has ≥ 200 significant genes (median 2,929,
10th percentile 369), so the cap always applies: a random 100 of the significant genes.

| metric (vs own control) | expr `_min15` s0 / s1 | expr `_qsig` |
|---|---|---|
| `spec_prob_timedose` | 7.8 / 10.4 | 1.6 (p = 0.22) |
| `rank_probc` | — | 1.2 (p = 0.024) |
| `rank_prob` | — | 0.2 (p = 0.42) |

Paired by drug (all 1,358), `_min15` − `_qsig` on `spec_prob_timedose`: +6.0 (s0, p = 3e-5) and +8.7 (s1, p = 2e-6).

**Reading.** This separates the two confounds in 5.2. Expr `_qsig` loses no training data and still loses almost
all of its signal, so the loss comes from **which** 100 features are read: a random 100 of the significant genes
carries almost nothing, and the largest changes carry the signal. For MetChange the same mechanism most likely
explains the drop in 5.2 (MetChange `_qsig` keeps more signal, 6.1 points, probably because 100 of ~4,400 MetChange
features still covers a broad share of a coherent metabolic response). Conclusion: **feature selection by magnitude
stays**. If significance is used, it should rank (top 100 by p or t), not threshold-then-sample. Whether
`--signature-scale sd` helps is the open test (`_sdscale`).

### 5.4 Significant, then strongest: `_qsigtop`

`_qsig` sampled 100 significant features at random, so 5.2–5.3 do not test the concern behind it: the current top-100
by |log2FC| can include large but non-significant changes. How often it does:

| among the top 100 by \|log2FC\| | MetChange (3,565 active perturbations) | expr (sample of 607 active perturbations) |
|---|---|---|
| significant (q < 0.05), median | 1 | 84 |
| perturbations with all 100 significant | 35.5% | 23.7% |
| perturbations with none significant | 57.3% | 0% (min 4) |
| perturbations with < 90 significant | — | 55% |

(Expr: local `results/labels/A549/expr`, which matched only 607 of the coverage table's perturbation names; q and
|log2FC| are unaffected by the sign correction.) So the concern is real for both, and for MetChange it is the
coverage split of report 04: for most perturbations the signature read is entirely non-significant. Those
signatures nevertheless carry signal (5.4–5.6 points over control, section 5), so non-significant is not the same
as noise; but a signature backed by significance is the cleaner claim.

New option `--qvalue-top-n` (`metabolic_target_prediction.py`): features with q ≤ `--qvalue`, then the `--top-n` of
them by |value| (no random cap). Launchers `launch_predict_targets_A549_{mc,tf,expr}_noGlobals_lr1e-3_min15_qsigtop{,_shuff}.sh`,
wrapper `submit_train_w_shuff_noGlobals_lr1e-3_min15_qsigtop.sh`.

Expectations: MetChange keeps the same 1,791 perturbations as `_qsig`; on the 35% fully covered ones the signature
equals `_min15`'s, so the comparison on shared perturbations isolates the cost of dropping the uncovered half from
training. Expr keeps all perturbations and replaces on average ~26 non-significant genes per signature with the next
strongest significant ones: if it matches `_min15`, significance can be required at no cost (preferred for
interpretation); if it is better, it should become the recipe.

MetChange also runs as `_qsigtop_sdscale` (`--qvalue-top-n` with `--signature-scale sd`): significant metabolites,
ranked by SD-scaled |log2FC|. Scaling does not change q-values, so it keeps the same perturbations as `_qsigtop`; it
only changes which significant metabolites are read where more than 100 are significant (the covered perturbations,
median ~4,400). Reading: `_qsigtop_sdscale` vs `_qsigtop` = scale effect on the same perturbations; vs `_min15_sdscale`
= the effect of requiring significance. Same wrapper.


### 5.5 BH over all 4,086 features vs over the 3,747 mapped ones

MetChange q-values are BH over all 4,086 compartment-specific features; 339 (currency etc.) are not network nodes
and are dropped before selection. Recomputing BH on the 3,747 mapped features only (local MetChange p-values, 607
active perturbations that matched by name; the local BH reproduces `qvals.csv` to 5e-15): median significant mapped
features per perturbation 315 → 319; 2,923 extra calls on 606,928 (+0.5%); 6 of 607 perturbations cross the
`--min-sources 3` threshold. The unmapped rows are 7.7% of all significant calls. Negligible: keep `qvals.csv`.

## 6. Schuhknecht benchmark: inputs (`data/validation/schuhknecht/`)

Generated on the cluster by `scripts/schuhknecht_map/` (`map_gem_cpds_with_table.ipynb`: GEM ↔ ion map;
`create_drugmap.py`: LINCS ↔ Schuhknecht drug/condition match for the earlier MetChange-vs-Schuhknecht correlation).

**Metabolites.** The GEM's KEGG / HMDB ids are matched to the ion annotations (`Metainformation_Ions`, "Annotation IDs").
Multiple matches were resolved by hand (30 cases); `dup_after_shuhknecht_map` = several GEM metabolites map to the
same ion (isomers / ion annotated with several candidates). `metabolite_match_table_non_ambig_match.csv` is the
non-ambiguous map **restricted to compartment `c`** (435 rows); the full table `cpd_info_wSchuhknecht.csv` has all
compartments. Using the full table, non-ambiguous, intracellular (not `e`), network nodes:

| | ions | network nodes |
|---|---|---|
| non-ambiguous, intracellular | 473 (445 with a `c` copy) | 959 (1–6 per ion: 177 / 160 / 98 / 24 / 12 / 2) |
| … of which MetChange signature features (nodes that were sources in training) | 306 | 571 |
| ambiguous (dropped so far), intracellular | 148 | 803 |

**Drugs.** `drug_match_table.csv`: 750 Schuhknecht drugs matched to LINCS A549 by SMILES (68 exact time/dose,
682 nearest time; 72 flagged `validated_metabolic`). The `matched_column` refers to the old `metabo_pert`
averaged-delta columns; it is not needed for the benchmark (the query is the Schuhknecht profile itself). Of the
750: 646 have known targets (candidates for ranking); 293 are among the A549 drugs evaluated in CV (active
perturbations), so they were training positives in 4 of 5 folds; the other ~353 never entered training.

**Design notes for the benchmark script (not written yet):**
1. Query sources: Schuhknecht log2FC (48 h, 10 µM; sign already "positive = more", as `--signature-sign -1` makes
   the MetChange signatures), select as in training (top 100 by |log2FC| among q < 0.05 if `_qsigtop` wins), each ion
   spread ±1/n over its intracellular copies. Restrict to the 571 nodes that were training sources (306 ions) by
   default; all 959 as a sensitivity check (sources the model never saw). Ambiguous ions: optionally spread over all
   their candidate metabolites.
2. Model per drug: the fold model in which the drug was held out (the 293); the 5-fold ensemble average for drugs that
   never entered training. Report the two groups separately.
3. Ranking: score all 646 candidate drugs from the protein scores (as `rank_prob` / `rank_probc`), rank of the true
   drug, against a shuffled-query control (Schuhknecht profiles permuted across drugs) and the signature-similarity
   baseline.
4. Train/query mismatch: training sources are single compartment nodes with weight ±1, query sources are spread. If
   the benchmark is weak, train a variant on compartment-collapsed MetChange signatures.

### 6.1 Benchmark script (written, not yet validly run)

`scripts/predict_targets/benchmark_schuhknecht.py` + `launch_benchmark_schuhknecht_A549.sh` implement the design
above (README, "Schuhknecht benchmark"). Candidates and queries: the 646 matched drugs with targets (all have ≥ 3
ions; median 100 of the 306 trained-node ions). Metrics `spec`, `rankc`, `rank` per drug, by group (`heldout_cv`
293, `never_trained` 353, `validated_metabolic` 70), with a permutation null on the query → drug assignment.

**Local runs are invalid.** A first local run gave chance everywhere (spec median 0.50–0.51 for the real runs and
their controls, both seeds), but the script's reproduction check fails locally: rebuilding held-out LINCS
perturbations of fold 0 (label columns identical between the home and scratch copies) gives drug scores with
r = 0.32–0.68 against the saved `test_scores_heldout.npz`, and the protein-level probabilities differ too. The static
per-protein features already disagree with the training ranges (e.g. `activity_down_hop1` max 0.2519 locally vs
0.2523 among the training rows, impossible with identical inputs), so the local network / baseline / target files
are not the ones the cluster runs used. Those outputs were deleted; the benchmark has to run on the cluster, where it
checks itself before scoring. (Also: torch on MPS crashes on this macOS; the script uses CUDA or CPU, ~2 min per
run on CPU.)

