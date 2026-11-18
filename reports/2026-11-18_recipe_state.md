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
`gene` network for expr / TF), so it runs on the cluster: `sbatch launch_baseline_signature_knn_A549.sh` (the three
`_qsigtop` runs and their controls, both representations, then `compare_runs.py` per pair and
`performance_figure/vs_baseline/`). Smoke-tested end to end here with `--allow-mismatch`: both representations,
`compare_runs.py` (including the dose/time-stratified spec from the `.npz`) and the figure row — code path only.
