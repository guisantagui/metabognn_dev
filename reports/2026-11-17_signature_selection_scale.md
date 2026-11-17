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
and are dropped before selection. Recomputing BH on the 3,747 mapped features only (local MetChange p-values — an older label version with 4,086 features; the training labels on the cluster have
5,707 — 607 active perturbations that matched by name; the local BH reproduces `qvals.csv` to 5e-15): median significant mapped
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
| … of which MetChange signature features (nodes that were sources in training) | 412 | 801 |
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
perturbations of fold 0 gives drug scores with r = 0.32–0.68 against the saved `test_scores_heldout.npz`, and the
protein-level probabilities differ too. Two local inputs differ from training: (a) the local MetChange labels are an
uncorrected copy (the corrected labels exist only on the cluster), so the signature features cannot match; (b) the
signature-independent per-protein features also disagree with the training ranges (e.g. `activity_down_hop1` max
0.2519 over all proteins locally vs 0.2523 among the training rows, impossible with identical inputs), so the local
network or baseline differ from the cluster copies as well. Those outputs were deleted; the benchmark has to run on the cluster, where it
checks itself before scoring. (Also: torch on MPS crashes on this macOS; the script uses CUDA or CPU, ~2 min per
run on CPU.)

**Cluster verification (2026-11-17).** On the cluster, with the run's own inputs, the reproduction check passes:
8/8 held-out LINCS perturbations of fold 0, r = 1.000000, max |diff| 8.8e-5 (GPU training vs CPU here). The
training MetChange labels (scratch, corrected) have 5,707 features, 5,218 of them network nodes; the local copy is an
older version (4,086 / 3,747). With the training labels, 412 non-ambiguous Schuhknecht ions map to 801 trained
intracellular nodes (section 6 table corrected; the local count was 306 / 571). Full benchmark:
`sbatch launch_benchmark_schuhknecht_A549.sh`.


### 6.2 Schuhknecht benchmark, cluster results (seed 0) and the association pipeline

Full benchmark, `metchange_rank_noGlobals_lr1e-3_min15` seed 0 (reproduction check OK): chance on every metric
(all 646 drugs: `spec` 0.500, p = 0.49; `rankc` 0.519; `rank` 0.499), its `_shuff` control likewise (0.532 / 0.503 /
0.501). Never-trained drugs: `spec` 0.469 (p = 0.11; control 0.505). Seed 1: also chance (all drugs `spec` 0.485,
p = 0.18; `rankc` 0.491; `rank` 0.507; `validated_metabolic` (70) `spec` 0.418, p = 0.085, `rankc` 0.435, p = 0.13, not
significant and not seen in seed 0: 0.492 / 0.478). Seed 1 control: cut by the time limit after scoring, rerun. The model reads LINCS MetChange
signatures but not measured metabolomics; first suspect: little shared signal between MetChange (inferred
transcriptome, ≤ 24 h) and Schuhknecht (measured, 48 h) for the same drug.

To test that directly on the labels the model uses, the earlier association analysis (averaged deltas,
`scripts/mc_schuhk_assoc`) is rebuilt on the limma labels in `scripts/schuhknecht_mc_assoc/`:
`prepare_assoc_inputs.py` (matched MetChange log2FC / q, Schuhknecht log2FC / p / q, expression producers baseline
from the expression log2FC via the network's GPRs and reaction → metabolite edges, drug match re-derived on the label
columns as in `create_drugmap.py`; compartment `c` as before, or `intracellular` with copies averaged),
`overlap_enrichment.py` (new: MetChange-significant vs Schuhknecht-significant metabolites per drug, Fisher / CMH,
mismatched-drug permutation for overlap, direction-concordance and a size-matched expression baseline), the two R
scripts made non-interactive (inputs as arguments, one label set allowed, hard-coded paths and exploratory lines
removed), and `launch_schuhknecht_mc_assoc_A549.sh`. All steps run end to end locally on the old local labels
(numbers not meaningful). Local counts: 306 metabolites (`c`) / 327 (`intracellular`), 750 drugs; with the local
labels ~610 of the 750 matched perturbations have no significant MetChange metabolite, so the overlap test rests on
the ~115–140 drugs with calls in both datasets.

## 7. `_qsigtop` results (q < 0.05, then the top 100 by |log2FC|)

Seed 0, corrected labels, each against its own control. Per-drug gap (`spec_prob_timedose`, points, control −
real), then paired by drug against the earlier selections (same drugs; "same perts" = per-drug means over shared
perturbations only). `rank_prob` / `rank_probc` are raw ranks, so they only compare runs with the same candidate set
(the `_qsig*` MetChange runs have 758 drugs, the others 1,358).

| run | drugs | gap vs own control (p) | vs `_min15` s0 | vs `_min15` s1 | vs `_qsig` |
|---|---|---|---|---|---|
| expr `_qsigtop` | 1,358 | **10.5** (4e-16) | +2.2 (p 0.12) | −0.3 (0.74) | +8.0 (3e-8) |
| TF `_qsigtop` | 1,346 | **7.8** (2e-11) | −0.3 (0.38) | −0.2 (0.44) | +5.1 (5e-5) |
| MetChange `_qsigtop` (folds 1–4 only) | 606 | **10.3** (3e-10) | +0.9 (0.39); same perts −0.1 (0.85) | +0.6 (0.72); same perts −1.9 (0.99) | +6.0 (0.01) |
| MetChange `_qsigtop_sdscale` | 758 | **8.1** (4e-7) | −0.8 (0.20); same perts −3.7 (8e-3) | −0.8 (0.22); same perts −3.2 (0.06) | +3.8 (0.10) |

- MetChange `_qsigtop` **fold 0 is missing** (no `fold0/` in the download, cluster `compare` summary also on folds
  1–4), most likely a failed array task as for `_sdscale` fold 3. Numbers above are on 4 of 5 folds.
- `_qsigtop_sdscale` vs `_qsigtop` (same 606 drugs): −2.0 points (p = 0.12); `rank_prob` −4.4 ranks (p = 0.015).
- Cross-signature, same drugs: expr `_qsigtop` vs TF `_qsigtop` +1.5 (p = 0.29); MetChange `_qsigtop_sdscale` vs
  expr `_qsigtop` −1.9 (p = 0.07; same perts −3.1, p = 0.05).
- MetChange `_sdscale` (4 of 5 folds, fold 3 still to rerun): 7.1 points; vs `_min15` −1.9 (s0, p = 0.36) and +0.7
  (s1, p = 0.82).

**Reading.**
1. **Requiring significance costs nothing.** Significant-then-strongest equals the unfiltered top 100 for all three
   signatures (paired differences within ±2 points, p ≥ 0.12), and beats the random-significant `_qsig` by 5–8
   points. So the `_qsig` losses were the random sampling, not the significance filter.
2. MetChange `_qsigtop` trains on half the perturbations (the covered ones) and does as well on them as `_min15`,
   which trained on all: the uncovered perturbations add nothing to the covered ones. It cannot score the
   uncovered half, so `_min15` stays the run with full drug coverage.
3. **SD scaling does not help**, neither alone (`_sdscale`, provisional) nor on top of significance
   (`_qsigtop_sdscale` ≤ `_qsigtop`). The scale bias of section 2 changes which metabolites are read but not how well
   the drug is recovered; the high-variance metabolites carry as much drug information as the rest.
4. Recipe: `_qsigtop` is the interpretable default (every signature feature is significant) at no cost; for
   MetChange it is restricted to the 758 covered drugs. Expr and TF can switch at no cost. Pending: MetChange
   `_qsigtop` fold 0, `_sdscale` fold 3, and a second seed if `_qsigtop` becomes the recipe.

## 8. Schuhknecht: benchmark complete, and the association on the limma labels

### 8.1 Benchmark, seed 1 control

`metchange_rank_noGlobals_lr1e-3_min15_seed1_shuff` (rerun, reproduction check OK): `spec` 0.503 (all), never_trained
0.472 (p = 0.12), validated_metabolic 0.483. The never_trained dip (0.469 / 0.480 in the real seeds) appears in the
controls too (0.505 s0, 0.472 s1), so it is not signal. **The benchmark is null in both seeds, real and control.**

### 8.2 MetChange vs Schuhknecht on the limma labels (`results/mc_schuhk_assoc/limma_rank_{c,intracellular}`)

750 matched drugs (68 exact dose/time, 682 nearest time), 416 metabolites (`c`), MetChange sign flipped to
production. The two compartment settings agree; numbers for `c`.

**Coverage.** MetChange has significant metabolites (q ≤ 0.05) for **85 of 750** matched drugs, and when it has
any it has most (median 365 of 416). Schuhknecht: median 7 significant metabolites per drug (57 drugs none). So the
overlap test rests on 78 drugs, and on those MetChange's call set is nearly the whole metabolite list.

**Overlap enrichment** (`overlap_enrichment/summary.tsv`):

| | MetChange | expression producers, size-matched |
|---|---|---|
| overlap with Schuhknecht-significant, observed / mismatched-drug null | 1,002 / 882 (z 1.9, p = 0.018) | 783 / 700 (z 1.8, p = 0.020) |
| sign-concordant overlaps, observed / null | 485 / 442 (p = 0.29) | 375 / 352 (p = 0.25) |
| directional score (concordant − discordant) | −32 vs +2 (p = 0.59) | −33 vs +5 (p = 0.80) |
| pooled sign concordance among overlaps | 48% | — |

Per-drug Fisher: 1 of 78 at q ≤ 0.1. Overlap is slightly above the mismatched-drug null, but equally so for the
expression baseline of the same size, and the directions do not agree. Paired MetChange vs expression directional
score: p = 0.45.

**Whole-profile checks (all 750 drugs, laptop, on the same matrices):**

| | MetChange | expression producers |
|---|---|---|
| per-drug Spearman with Schuhknecht log2FC, median (mismatched null mean) | −0.007 (−0.001) | +0.008 |
| same, metabolites centred across drugs | −0.001 (Wilcoxon p = 0.38) | — |
| sign agreement on Schuhknecht-significant metabolites (459 drugs with ≥ 5), observed / mismatched null | **0.569 / 0.541 ± 0.014 (p = 0.02)** | 0.394 / 0.400 (p = 0.80) |

(`intracellular`: 0.569 / 0.542, p = 0.03; expression 0.425 / 0.432.)

**R benchmark** (`metchange_vs_expression_benchmark.R`, run on the laptop from the cluster matrices, the cluster
had no R; `limma_rank_c/benchmark/`, log `benchmark.log`). Default selection: the 188 drugs in the top quartile of
Schuhknecht response, 355 metabolites. `rhoc` = per-drug Spearman with the Schuhknecht log2FC after centring each
metabolite across drugs (the drug-specific part).

| | MetChange | expression producers |
|---|---|---|
| mean `rhoc`, 188 drugs (p vs 0) | **−0.013** (0.08) | **+0.033** (5e-6) |
| MetChange − expression per drug | −0.046, MetChange better in 35% of drugs (p = 7e-6) | |
| specificity permutation (mismatched drugs), MetChange | observed −0.013 vs null 0.000, p = 1.0 | |
| pooled mixed model, slope (R² marginal) | −0.054 (0.0012) | +0.032 (0.0008) |
| MetChange adds to expression (LRT) | p = 0.0025, **with a negative slope** | |
| same, all 750 / 408 / 78 drugs | MetChange − expression −0.012 / −0.026 / −0.058, all p < 0.02 | |
| gap vs expression, by number of Schuhknecht-significant metabolites (quartiles) | +0.013, −0.011, −0.008, −0.044 | |
| direction AUROC on Schuhknecht-significant metabolites (90 drugs) | 0.504 | 0.521 |
| sign accuracy, top 5% / 10% / 20% of MetChange values | 0.552 (p = 0.008) / 0.511 / 0.494 | 0.492 / 0.507 / 0.528 |

`intracellular` (394 metabolites, 188 drugs) is the same: MetChange `rhoc` −0.015 vs expression +0.037,
difference −0.052 (p = 2e-6), specificity p = 1.0, pooled slope −0.060 (LRT p = 7e-4), sign accuracy at the top
5% / 10% / 20%: 0.548 (p = 0.008) / 0.502 / 0.500.

**Reading.**
1. On the limma labels, **MetChange does not agree with Schuhknecht for the same drug.** Across the profile its
   drug-specific correlation is zero to slightly negative (−0.013; specificity permutation p = 1.0) and below the
   naive expression-producers baseline (+0.033), increasingly so for the drugs with the strongest measured
   response. The earlier "MetChange adds to expression" (delta-based labels) reappears only as an LRT with a
   negative slope, i.e. not as agreement. The only positive hints are sign agreement on the most extreme values:
   on Schuhknecht-significant metabolites 57% vs 54% for a mismatched drug (p = 0.02), and for MetChange's top 5%
   values 55% (p = 0.008), vanishing at the top 10–20%. Small, inconsistent across thresholds, not a usable signal.
   (The sign convention is the training one: corrected labels, `−log2FC` = production; with the opposite sign the
   correlation would be +0.013 and the sign agreement 43% vs 46%, no better.)
2. The expression baseline itself is weak (ρ ≈ 0.03), so measured 48 h metabolomics and LINCS ≤ 24 h
   transcriptomics share little drug-specific signal in either representation.
3. MetChange's significance is not informative here: it is on for 11% of drugs and then covers most metabolites.
4. This explains the null benchmark (8.1). The model does read MetChange signatures (sections 5, 7), but the
   measured profile of the same drug does not resemble them, so a MetChange-trained model cannot recover the drug
   from measured metabolomics. The bottleneck is the label (MetChange vs measured), not the model. Contributing:
   timepoints (≤ 24 h vs 48 h), nearest-time match for 91% of drugs, MetChange inferred from partly CycleGAN-inferred
   transcriptomes. For the Schuhknecht query to be meaningful, the model would need to be trained or calibrated on
   measured metabolomics (e.g. Schuhknecht itself, cross-validated over drugs), or the query compared at the level
   where the two agree best.

## 9. What can validate the tool, given 8.2

The intended claim, that the tool finds compounds that change metabolism in a given way, cannot be validated
through MetChange or expression labels: neither agrees with measured metabolomics of the same drug at the whole-profile
level (MetChange ρ ≈ −0.01, expression producers ≈ 0.03). Part of that is data mismatch (48 h vs ≤ 24 h,
nearest-time match for 91% of drugs, partly inferred transcriptomes), but it rules out validation through the
inferred labels. An expr model supports only the weaker claim of compounds that cause a given transcriptional
change in metabolic genes. MetChange is a network reading of the same transcriptomes, not a predictor of metabolite
levels.

Ways forward, in order of strength:

1. **Train on measured metabolomics (Schuhknecht as training data, not as a query).** About 1,500 drugs in A549,
   more than the LINCS A549 set. Use the same network model with Schuhknecht profiles as signatures, cross-validated
   over drugs, against its own shuffled control. That directly tests "a measured metabolic change predicts the
   targets of the drugs that cause it", with no inference step in the labels. LINCS / MetChange / expr then become
   pretraining or added features, judged by whether they help on held-out Schuhknecht drugs. The query building
   from `benchmark_schuhknecht.py` and the feature × perturbation signature input of the training script can be
   reused.
   **First check (laptop, ~1 h):** are the Schuhknecht profiles of drugs that share targets more similar than those
   of random drug pairs? If not, this path fails at the data level.
2. **Knowledge-based queries at the enzyme level.** Build the expected substrate-up / product-down change around
   the target of drugs with simple, documented metabolic mechanisms (methotrexate / DHFR, brequinar or leflunomide /
   DHODH, mycophenolate / IMPDH, statins / HMGCR, 5-FU or pemetrexed / TYMS, complex I inhibitors), and check that
   their inhibitors rank near the top. This tests the network reasoning, not the labels, and needs no
   high-throughput data. The Schuhknecht `validated_metabolic` drugs (70 with targets) are the natural set; the
   measured profile can also be restricted to the metabolites around the target.
3. **Prospective experiments:** targeted metabolomics or isotope tracing (or Seahorse) for a few predictions chosen
   by 1 or 2. This is the only direct proof of the claim.

Plan: run the shared-target similarity check (1) first; if it holds, set up Schuhknecht training, with (2) as a
second, interpretable benchmark.

Also: `launch_schuhknecht_mc_assoc_A549.sh` now runs the R step in the conda env `r-4.3.1` on the cluster (Python
steps first, then `conda activate r-4.3.1`).

### 9.1 Revised after discussion

- **Schuhknecht as training data: deprioritised.** One cancer cell line and ~1,500 drugs, against hundreds of cell
  lines and thousands of perturbations in LINCS. And intracellular metabolomics measures pools: a larger pool can come
  from more production or less consumption. MetChange's intended claim is about production / consumption capacity
  (e.g. drugs whose transcriptional response supports more NAD+ production), which pool measurements cannot test.
  The Schuhknecht null is therefore not evidence against that claim, only against reading MetChange as pool change.
  Pathway-level metabolomics can point to the pathway a drug hits, but not to the target (Sarvin et al., Nat Commun
  2026, as described by the user).
- **The expr model is not restricted to metabolic genes:** `launch_predict_targets_A549_expr_*` uses all DE genes on
  `networks/gene` (protein network with transcriptional edges, no GEM layer). It is a generic transcriptome → target
  model, the same class as the TF positive control. MetChange's specific value is the stoichiometric reading
  (production vs consumption of a metabolite), and that is what needs validating.
- **Validation split into two links:** (a) signature → target, already validated by target recovery in CV
  (sections 5, 7); (b) transcriptome → MetChange → production / consumption, which belongs to MetChange's own
  validation. Enzyme-level ground truth for the full chain:
  1. **LINCS genetic perturbations of metabolic enzymes** (shRNA `trt_sh`, CRISPR `trt_xpr`, overexpression `trt_oe`)
     in the same cell lines. Compute MetChange labels for them, query the trained model, and check whether the
     perturbed enzyme ranks near the top. The ground truth is exact, with no polypharmacology and no new data
     modality. **Most direct next step.**
  2. **Inborn errors of metabolism:** enzyme deficiency → known biomarker metabolites, the standard GEM
     benchmark (Shlomi et al., Mol Syst Biol 2009; IEMbase, OMIM / HMDB). The query is the biomarker profile, the
     answer the deficient enzyme. Caveat: biomarkers are plasma / urine pools.
  3. **Derived queries from the GEM:** for drugs annotated as inhibitors of a GEM enzyme (ChEMBL mechanism /
     DrugBank actions), the expected production / consumption change around the enzyme's reaction. Self-built, so
     partly circular with the network.
  4. Flux-level measurements (isotope tracing, exometabolomics uptake / secretion, Seahorse) for selected
     predictions, prospectively.
