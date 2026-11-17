# Report 06 — Where the recipe stands (not frozen): signature selection settled, scaling dropped

*2026-11-18. No new runs today. Summary of the state after the selection / scale experiments of report 05
(`reports/2026-11-17_signature_selection_scale.md`, sections 5, 7) so that the work can be taken over later.
Metric and gap convention: report 01 section 3, report 03 section 1. Gap = per-drug median percentile of the own
shuffled-signature control minus that of the run, in points, on `spec_prob_timedose_heldout`; paired by drug
between runs.*

## 1. Where we are

| signature | `_min15` (top 100 by \|log2FC\|) | `_qsigtop` (q < 0.05, then top 100) | `_qsig` (random 100 significant) | `_sdscale` |
|---|---|---|---|---|
| MetChange | 8.3 / 7.7 (seeds 0 / 1; 1,358 drugs) | **10.3** (606 drugs, folds 1–4, provisional) | 6.1 (758) | 6.7 (1,358) |
| expr | 7.8 / 10.4 | **10.5** (1,358) | 1.6 | — |
| TF | 9.0 / 9.3 | 7.8 (1,346) | 4.2 | — |

MetChange `_qsigtop_sdscale`: 8.1 (758 drugs).

**Best complete result: expr `_qsigtop`, 10.5 points** (p = 4e-16, 1,358 drugs, 5 folds, one seed).

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

- **MetChange `_qsigtop` is provisional:** fold 0 of the real run is missing (array task failure; the `_shuff`
  control and the expr / TF runs have all 5 folds), so 10.3 is on folds 1–4 / 606 drugs and there is no
  `cv_metrics.json`. **Fold 0 relaunched 2026-11-18, running.** When it lands: re-run `compare_runs.py` for that
  run and redo the paired comparisons (`_qsigtop` vs `_min15` s0 / s1 and vs `_qsigtop_sdscale`).
- **One seed only** for every `_qsigtop` run. A second seed is needed before calling `_qsigtop` the recipe,
  as the expr seed spread shows (`_min15` 7.8 vs 10.4 across seeds is larger than most of the differences
  being compared).
- **Not frozen on purpose.** The selection question is answered; the open ones are representation (continuous
  MetChange values or a larger `--top-n` instead of 100 signs, `--signature-weights`) and everything in
  STATUS 6.2b (validation).
- `_sdscale` is complete and can be dropped from future comparisons.
