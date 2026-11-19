# Report 06 — Where the recipe stands (not frozen): signature selection settled, scaling dropped

*2026-11-18. No new runs today. Summary of the state after the selection / scale experiments of report 05
(`reports/2026-11-17_signature_selection_scale.md`, sections 5, 7) so that the work can be taken over later.
Metric and gap convention: report 01 section 3, report 03 section 1. Gap = per-drug median percentile of the own
shuffled-signature control minus that of the run, in points, on `spec_prob_timedose_heldout`; paired by drug
between runs.*

## 1. Where we are

| signature | `_min15` (top 100 by \|log2FC\|) | `_qsigtop` (q < 0.05, then top 100) | `_qsig` (random 100 significant) | `_sdscale` |
|---|---|---|---|---|
| MetChange | 8.3 / 7.7 (seeds 0 / 1; 1,358 drugs) | **10.2** (758 drugs, complete) | 6.1 (758) | 6.7 (1,358) |
| expr | 7.8 / 10.4 | **10.5** (1,358) | 1.6 | — |
| TF | 9.0 / 9.3 | 7.8 (1,346) | 4.2 | — |

MetChange `_qsigtop_sdscale`: 8.1 (758 drugs).

**Best complete result: expr `_qsigtop`, 10.5 points** (p = 4e-16, 1,358 drugs, 5 folds, one seed), with
MetChange `_qsigtop` at 10.2 on its 758 covered drugs.

**MetChange `_qsigtop` fold 0 landed 2026-11-18 00:40**, `compare_runs.py` re-run: **10.2 points** over its own
control (`rank_probc` 10.9, `rank_prob` 5.8, all p < 1e-4; 61% of drugs better). Paired by drug: vs `_min15` s0
+0.3 (p = 0.76) and s1 +1.1 (p = 0.82) on its 758 drugs — on shared perturbations −1.0 / −2.7 (p = 0.32 / 0.47);
vs `_qsig` +6.0 (p = 0.0023); vs `_qsigtop_sdscale` +1.8 (p = 0.12). The 4-fold provisional value was 10.3, so
nothing changed with fold 0: **`_qsigtop` ties `_min15` on the drugs it covers**, as reported below.

## 2. Recipe as it stands (the `_qsigtop` launchers)

```
--signature-sign -1 --mode metchange
--qvalue 0.05 --qvalue-top-n --top-n 100
--de-counts <expr labels>/fit_info_nsig.csv --min-de-genes 200 --score-active-only
--lr 0.001 --patience 15 --ffn-epochs 80 --min-epochs 15
--no-vae --no-globals --n-folds 5
```

no `--signature-scale`; network `data/networks/gene` for expr / TF, `combined` for MetChange; labels from the
corrected scratch copies. Launchers `launch_predict_targets_A549_{mc,tf,expr}_noGlobals_lr1e-3_min15_qsigtop{,_shuff}.sh`,
wrapper `submit_train_w_shuff_noGlobals_lr1e-3_min15_qsigtop.sh`. **Not frozen** (section 4).

## 3. What is settled

1. **Requiring significance is free, so keep it.** `_qsigtop` matches `_min15` for all three signature types
   (paired differences within ±2 points, p ≥ 0.12) and beats `_qsig` by 5–8 points (p ≤ 5e-5). Every feature read
   is then significant, which is the easier claim to defend. Expr and TF can switch with no trade-off.
2. **SD scaling is out.** `_sdscale` on full coverage is 6.7 points, vs `_min15` −1.6 / +0.4 (p = 0.46 / 0.66),
   and `_qsigtop_sdscale` ≤ `_qsigtop` (−1.3, p = 0.62). The metabolite-scale bias of report 05 section 2 changes
   *which* metabolites are read but not how well the drug is recovered.
3. **Threshold-then-sample is out.** `_qsig` loses most of the signal (expr 1.6 points, p = 0.22); the magnitude
   ranking does the work. If significance is used it must rank, not sample.
4. **The MetChange trade-off.** `_qsigtop` scores only the 758 drugs that have significant metabolites, and on
   those it ties `_min15`, which trained on all 1,358. So the uncovered half adds nothing to the covered one but
   also cannot be scored. **Primary MetChange run stays `_min15`** (full coverage); `_qsigtop` is the
   significance-backed variant on its subset. Expr and TF have no such trade-off.

## 4. Caveats and pending

- ~~MetChange `_qsigtop` fold 0 missing~~ **done 2026-11-18** (section 1): all runs in the table are complete,
  5 folds each.
- **One seed only** for every `_qsigtop` run. A second seed is needed before calling `_qsigtop` the recipe,
  as the expr seed spread shows (`_min15` 7.8 vs 10.4 across seeds is larger than most of the differences
  being compared).
- **Not frozen on purpose.** The selection question is answered; the open ones are representation (continuous
  MetChange values or a larger `--top-n` instead of 100 signs, `--signature-weights`) and everything in
  STATUS 6.2b (validation).
- `_sdscale` is complete and can be dropped from future comparisons.

## 5. Performance figure

`scripts/predict_targets/plot_performance_panels.py` (new), run on the laptop. Three panels, all on the
`_qsigtop` runs, each against its own `_shuff` control. Presentation format (11 x 5.95 in, 2:1): the runs are
labelled only **MetChange / expression / TF** (the selection is the same everywhere, so `_qsigtop` is not
written on the figure), no A / B / C panel letters, and a **transparent background** (`--background surface`
for an opaque one). The ink is near-black, i.e. for a light slide; `--theme dark` re-renders with the dark
steps of the same hues and light ink for a dark slide, writing `..._dark.{png,pdf}` (both themes are
transparent, and the dark one is a selected palette, not an inverted light one). Output
`results/predict_targets/A549/performance_figure/performance_panels_<metric>.{png,pdf}` + `.tsv` with the
numbers behind every panel.

- **A** ECDF of the percentile per drug, real vs own shuffled control, one panel per signature — the full
  distribution, not a summary. The real curve is above the control everywhere for all three.
- **B** Forest of the per-drug gap with bootstrap CIs (5,000 resamples of the median): the three runs only,
  so no other recipe has to be explained. On `rank_probc` the CIs overlap for all three (MetChange 10.0
  [8.2-13.0], TF 8.8 [6.7-10.7], expression 8.5 [7.4-10.9]) — none is distinguishable from the others.
  `--extra sig=run:label` puts `_min15` or ablation rows back when a selection comparison is wanted.
- **C** Top-k recovery at 1 / 5 / 10 / 25 %, real bars with each run's own control marked on the bar.

**The coverage-split panel was dropped.** It can only be the `_min15` runs (in `_qsigtop` every feature read is
significant by construction), so it is the *motivation* for `_qsigtop`, not a result of the current recipe, and
it is not how the work is framed. It stays available in `metabolic_coverage_fixedlabels/spec_by_coverage.tsv`
and in report 05 section 5.

### 5.1 Which metric to show

All three are rendered (`--metric`); they answer different questions and the choice changes the ordering:

| metric | question | MetChange / expr / TF `_qsigtop`, gap in points |
|---|---|---|
| `rank_probc_heldout` | **the deployment task:** given a query, where does the true drug rank among candidates, after standardising each candidate across queries (which removes the protein prior) | **10.0 / 8.5 / 8.8** |
| `spec_prob_timedose_heldout` | is the signature read at all: do the drug's targets score higher under its own query than under other drugs' queries of the same time × relative dose stratum (so the design cannot carry drug identity) | 10.2 / 10.5 / 7.8 |
| `rank_prob_heldout` | the same ranking without the per-candidate standardisation, so the protein prior still contributes | 5.8 / 3.4 / 4.0 |

**Use `rank_probc` as the headline** — it is the task the tool is for and it is prior-controlled — and
`spec_prob_timedose` beside it as the leakage-controlled check; `rank_prob` belongs in the supplement.
`spec_prob_timedose` was the metric of reports 01–05 because it is the one that is not affected by the
dose / time design, which is why the earlier tables quote it; the two agree on every conclusion so far
(`_qsigtop` ≈ `_min15`, scaling neutral, `_qsig` worse), they differ only in which signature looks best.

Caveat for all three: the per-drug percentile normalises the candidate-set size, but MetChange `_qsigtop`
ranks within a different (758-drug) candidate pool than expr / TF (1,358), so cross-signature differences are
read with that in mind — the paired comparisons in section 1 are on shared drugs and are the ones to trust.

**A bug found while adding this** (worth knowing for any reuse of `cv_drug_ranks.tsv`): the `rank_*_heldout`
columns are **raw ranks**, not percentiles — `compare_runs.py` divides them by `n_candidates_heldout`. Only the
stratified `spec_*` columns are already percentiles. The earlier scratch comparisons of `rank_probc` in report 05
section 7 are in rank units ×100, not points, which is why they read in the thousands; the conclusions there
(signs and p-values) are unaffected, but the magnitudes are not comparable with the `spec` ones. The plotting
script now normalises correctly.

Re-run when a second `_qsigtop` seed exists; `--extra sig=run:label` adds forest rows.

## 6. Signature-similarity baseline (written, not yet run)

The model-free check (`check_signature_target_similarity.py`, `results/signature_similarity/A549/`) already shows
that drugs sharing a target have similar signatures (expr: AUROC 0.69 at the same time point, 10-NN enrichment
9.9x; MetChange 0.57, 5.5x). What was missing is a *predictive* baseline on the model's own metric, to say
whether the network adds anything over that similarity.

`scripts/predict_targets/baseline_signature_knn.py`: for each held-out perturbation of a run, the k (25) most
similar **training** drugs vote their targets, similarity-weighted and encoded like the model's prior; the protein
scores then go through the model's own `DrugRanker.prob_scores`, `signature_specificity` and `centered_ranks`, on
the run's folds, queries and candidates, against the run's `_shuff` control with the shuffle replayed. Same
contrast, same candidates, same code: the gaps are directly comparable with the figure in section 5.

- `--representation sources`: the model's exact input (selected features, signs, on network nodes; cosine) —
  network vs no network on the same information.
- `--representation profile`: the full continuous signature (Spearman) — whether the data hold more than the
  100-sign input keeps (bears directly on STATUS 6.0b item 3).

Design choices: held-out drugs are never references (ranking a drug by similarity to its *own* other doses would
be the same-drug ceiling, AUROC 0.85 for expr, not a baseline); references are training drugs only, active
perturbations, as the model trains on; a drug's similarity is that of its most similar perturbation (`--unit
drug`), so drugs with many doses get one vote. The output is laid out like a run, so `compare_runs.py` and
`plot_performance_panels.py --extra` read it unchanged.

Each fold is checked first: the held-out drugs and test perturbations rebuilt from the run's own inputs must
match the run's (`fold<i>/test_scores_heldout.npz`, `test_drug_ranks.tsv`), otherwise it stops. On the laptop
they do not (older labels: 835 perturbations with significant MetChange features vs thousands on the cluster; no
`gene` network for expr / TF), so it runs on the cluster: `sbatch launch_baseline_signature_knn_A549.sh`, a 3-task array
(one per `_qsigtop` run). Each task makes **one** call that scores the run and its `_shuff`, both representations
and every k from a single read of the labels (the sides differ only by the replayed permutation; the similarity
matrix of a fold is reused for every k), then `compare_runs.py` per pair and `performance_figure/vs_baseline/`.
Wall time is one expr label read, ~30-50 min (the first version re-read the labels for every side and
representation, ~2.5-3 h; it was cancelled and replaced the same day). The single-read version reproduces the
first one's real-side numbers exactly on the smoke test. Smoke-tested end to end here with `--allow-mismatch`: both representations,
`compare_runs.py` (including the dose/time-stratified spec from the `.npz`) and the figure row — code path only.

### 6.1 Results (cluster run, downloaded 2026-11-18)

All 12 baseline directories (3 runs x sources / profile x real / `_shuff`, k = 25). **Input check passed for every
fold of every side**: held-out drugs identical to the run's candidates, every test perturbation found. The
baselines are on exactly the model's queries and candidates.

**Ties.** With k = 25 a query votes only for its neighbours' targets, so 68-78% of candidates get a score of exactly
0 (model: 0%, largest tie 1.5%). A true drug in that block takes its average rank (~0.63), which is why the
baseline's `rank_prob` control sits at 0.61-0.64 instead of 0.5, and `spec` is pinned at 0.5. **`rank_prob` and `spec`
are therefore not comparable for this baseline.** `rank_probc` standardises each candidate over the queries, which
breaks the ties without bias (controls 0.49-0.51), so the comparison below uses it only.

`rank_probc_heldout`, per drug (mean over its active perturbations), gap over own control with bootstrap CI;
real / control median percentile:

| | model | kNN, same input | kNN, full profile |
|---|---|---|---|
| MetChange (758) | 10.0 [8.2, 13.0]; 0.375 / 0.500 | 8.0 [6.4, 9.4]; 0.350 / 0.512 | 10.1 [8.7, 11.8]; **0.287** / 0.509 |
| expression (1,358) | 8.5 [7.5, 10.9]; 0.388 / 0.491 | 7.6 [6.6, 8.9]; 0.309 / 0.498 | 7.8 [6.3, 10.1]; **0.305** / 0.512 |
| TF (1,346) | 8.8 [6.7, 10.6]; 0.406 / 0.495 | 8.6 [7.8, 10.2]; 0.289 / 0.498 | 8.7 [7.3, 10.6]; **0.294** / 0.497 |

Paired by drug, model gap minus kNN gap: MetChange −3.8 (same input, p = 0.012) / −4.8 (profile, p = 3e-6); expr
−8.3 / −9.5 (p ≤ 7e-14); TF −8.2 / −8.5 (p ≤ 6e-13). The model wins for only 39-45% of drugs. Full profile vs same
input: +1.7 for MetChange (p = 0.03), +0.4 / +0.3 for expr / TF (n.s.).

**Complementarity.** Per-drug gaps of model and kNN correlate weakly (Spearman 0.18-0.31): they succeed on different
drugs. A plain rank-average of the two scores (within-query ranks, then the same `rank_probc`) has a median gap of
14-16 points (MetChange 15.9 [13.0, 19.3], expr 14.1, TF 13.7) and beats the model per drug by +5.2-5.8 (p ≤ 2e-7),
but not the kNN alone per drug (−2.0 to −4.6). Gap statistics involving the sparse kNN are noisy because of the
ties (the median gaps and the paired medians do not order the methods consistently), so these ensemble numbers
are indicative.

**Reading.**
1. **On the deployment metric, a no-network nearest-neighbour baseline ranks the true drug better than the
   trained model**: median percentile 0.29-0.31 vs 0.38-0.41 with both controls at ~0.50, and it wins for
   55-61% of drugs. By the reading rule of STATUS 6.0b, the network is not earning its place on this metric with the
   current data.
2. The 100-sign input is not the bottleneck for the kNN (full profile ≈ same input, except +1.7 for MetChange),
   which weakens the premise of STATUS 6.0b item 3 at least for similarity-based scoring.
3. The model is not redundant with the kNN: their per-drug successes barely overlap, and combining them helps
   the model. The model carries something the kNN does not, but by itself ranks worse.
4. **Not final until the ties are removed.** Next: a dense version of the baseline (every training drug votes,
   with a soft weight, e.g. softmax of the similarity), so no candidate scores exactly 0 and `rank_prob` /
   `spec` become comparable as well. Until then, do not present the model's performance without this caveat:
   the first question a reviewer will ask is exactly this baseline.

### 6.2 Dense version (written, launched 2026-11-18)

`baseline_signature_knn.py --beta B`: every training drug with active perturbations votes, weight softmax(B x z),
z = the query's similarities z-scored over the voting drugs (scale-free, so the same B means the same concentration
for cosine and Spearman), computed in float64 so a sharp B cannot underflow back to zeros. The only remaining exact
zeros are **orphans**: held-out drugs none of whose directional target tokens (gene x activation / inhibition) is
carried by any reference drug, so no target-transfer method can ever vote for them, while the network still scores
them. Measured on the MetChange folds: 7-12% of the candidates (vs 68-78% zero scores with top-25). Each test
perturbation carries `target_shared_with_references`, so results can be given with and without orphans (and the
orphans are one place where the model can legitimately beat any similarity baseline).

Smoke test: the zero-score candidates per query equal the token-level orphan count exactly, also at B = 10; the
top-25 numbers reproduce the earlier run. Launcher default is now `BETA="2 5 10"` (top-k off, already done;
`K=25` re-adds it); the figure's baseline rows use the first B. Read B = 2 (diffuse) vs 10 (close to top-k but dense)
as a robustness check: the conclusion of 6.1 holds only if it holds across B, and on `rank_prob` / `spec` as well
as `rank_probc`.

### 5.2 Statistics on the forest panel (added)

Each row now shows its p-value under the label. **Test:** Wilcoxon signed-rank, two-sided, on the per-drug gaps
(control percentile − real percentile, one value per drug) against 0. Paired, because every drug is scored in the
real run and in its own control with the same fold and candidates; drugs with a gap of exactly 0 are left out
(scipy default). **Interval:** percentile bootstrap of the median gap, 5,000 resamples of the drugs (with replacement,
original size, seed 0); the bar spans the 2.5th-97.5th percentile of the resampled medians. Drugs are the resampled
unit, so the bar covers which drugs happen to be in the set, not training randomness or fold assignment (the expr
seed spread, 7.8 vs 10.4, is wider than these bars).

| metric | MetChange | expression | TF |
|---|---|---|---|
| `rank_probc` | 10.0 [8.2, 13.0], p = 5e-13 | 8.5 [7.5, 10.9], p = 1e-16 | 8.8 [6.7, 10.6], p = 2e-17 |
| `spec_prob_timedose` | 10.2 [6.8, 13.5], p = 4e-11 | 10.5 [8.1, 12.6], p = 4e-16 | 7.8 [5.4, 9.8], p = 2e-11 |
| `rank_prob` | 5.8 [4.6, 7.3], p = 6e-20 | 3.4 [2.5, 4.3], p = 4e-15 | 4.0 [3.2, 4.9], p = 3e-16 |

The p-values differ slightly from `compare_runs.py`'s `paired.tsv`, which averages each drug's scores before
ranking; the figure averages per-perturbation percentiles. Both are paired by drug.

## 7. Seahorse as a query (design, nothing run)

Seahorse (Agilent XF) measures two **rates**, not pools, which is what MetChange's claim is about (report 05,
9.3): **OCR** (pmol O2/min, mitochondrial respiration) and **ECAR** (mpH/min, acidification: lactate + H+ export,
with a CO2 contribution; as PER / glycoPER once the buffer capacity is applied). The assays derive basal and
maximal respiration, ATP-linked respiration, proton leak, spare capacity, coupling efficiency (Mito Stress),
glycolysis / capacity / reserve (Glyco Stress), mitochondrial vs glycolytic ATP (Real-Time ATP Rate) and fuel
dependency on glucose / glutamine / fatty acids (Mito Fuel Flex). Coarse (whole-cell respiration and glycolysis,
no pathway resolution), but a genuine flux readout, medium throughput and cheap per compound.

**Query construction.** The first idea — OCR up -> "increased production" of every TCA and ETC metabolite, ECAR up
-> every glycolytic metabolite — is the right shape but likely too broad: at higher flux the intermediates of a
pathway often stay flat or fall while throughput rises, so a uniformly positive query may not resemble any real
MetChange signature. Two tighter versions:
- **endpoints only:** OCR up = O2 / NADH / FADH2 / acetyl-CoA consumed, CO2 / ATP produced; ECAR up = glucose
  consumed, lactate + H+ produced;
- **or only what the instrument actually senses** (lactate, CO2, ATP), without reconstructing the pathway.

**Use the two rates together.** Their ratio is the informative part: *OXPHOS up* (OCR up, ECAR down), *glycolytic
shift* (ECAR up, OCR down) and *both down* (general toxicity — effectively a negative control, since most
cytotoxic compounds land there). The opposing pattern is far more specific than either rate alone.

**Two possible benchmarks, testing different things:**

| | (a) query -> drug rank | (b) enrichment in the MetChange signatures |
|---|---|---|
| what | the Schuhknecht design: 2-3 phenotype queries, where do the drugs known to cause that phenotype rank among candidates | per drug, are its significant MetChange metabolites enriched for glycolysis / TCA-ETC with the right sign, and does that agree with the measured OCR / ECAR direction |
| tests | the tool's actual task | the labels (link b of 9.3) |
| statistic | rank enrichment (AUROC / GSEA) of the known positives, permutation null over which drugs are positive; no per-drug value, so power comes from the number of labelled drugs | one value per drug: properly powered, reuses `scripts/schuhknecht_mc_assoc/` |

**Order: (b) first** — cheap, better powered, and it says whether the labels carry the energy-metabolism signal at
all; if they do not, (a) cannot work. Then (a) as the headline, since it is the tool's task.

**Data is the binding constraint, check it before building anything.** Public Seahorse at drug-screen scale is
scarce and usually published as summary tables, not deposited. For (a) the raw rates are not even needed: curated
**known ETC / OXPHOS / glycolysis inhibitors** and **mitochondrial-toxicity screens** (Tox21,
glucose-galactose assays) give labelled positives directly, which is all the rank-enrichment test requires.
Searching for a LINCS-cell-line Seahorse dataset is the first step.
