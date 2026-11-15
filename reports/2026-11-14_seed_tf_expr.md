# Report 03 — Second seed, epoch cap, and the recipe on TF and expression signatures

*2026-11-14. Runs: `results/predict_targets/A549/{metchange_rank_noGlobals_lr1e-3_ep100, metchange_rank_noGlobals_lr1e-3_seed1,
tf_active200_noGlobals_lr1e-3, expr_noGlobals_lr1e-3}` and their `_shuff` controls. All use the reference recipe of
report 02 (`--no-globals --lr 0.001 --patience 10 --ffn-epochs 60`, unmatched, AUPRC selection, 5-fold CV by drug,
active perturbations, `--ignore-qvalue --top-n 100`, no VAE), array jobs. Metric definitions: report 01, section 3.*

**Framing.** MetChange is the signature this project is about; expression is a welcome byproduct. TF signatures
reproduce the setting of Menacho et al. 2026, which is published: they are the **positive control** that the
pipeline works where it should, and the level MetChange is measured against, not a result in themselves.

## 1. How the numbers are read

Every percentile is the rank of the correct answer among its candidates divided by the number of candidates (0 =
first, 0.5 = chance, lower is better). Per drug, a drug's per-perturbation percentiles are averaged within its fold
(`compare_runs.py`); the 1,358 held-out drugs are the independent units. Real and control score the same drugs, so
each drug gives one difference real − control.

**Gap** in this report = `median_diff` of `paired.tsv`: the median of those per-drug differences, in percentile points
(1 point = 0.01), sign flipped so that positive = real better. It is the quantity the Wilcoxon p-value tests. Report 02
quoted the difference of the two medians instead (e.g. "5.6 points"), which is not what its p-value tests; see the
correction note there. Medians below are the medians of each run on its own (`summary.tsv`).

## 2. The epoch cap did not bind (`_ep100`)

`metchange_rank_noGlobals_lr1e-3_ep100` (cap 100) gives **byte-identical** held-out scores to the cap-60 run, real and
control. Training is deterministic for a given seed, and no member reached epoch 60: the member that kept epoch 53
stopped by patience at 63, the rest earlier.

Raising the cap alone would not change much anyway: `StepLR(step_size=10, gamma=0.5)` halves the learning rate every
10 epochs, so at epoch 60 it is 1/64 of the starting value (1.6e-5) and at epoch 100 about 1e-6. What limits training
is early stopping, not the cap (section 5).

## 3. Second seed: the effect replicates, smaller

`--seed 1` redraws the drug folds, validation drugs, negative sampling, initialisation and the control's permutation:
a full replicate. Compare pooled and per-drug values, not folds.

| `spec_prob_timedose` | seed 0 | seed 1 |
|---|---|---|
| per perturbation, real / control | 0.414 / 0.500 | 0.425 / 0.488 |
| per drug, real / control | 0.444 / 0.500 | 0.455 / 0.493 |
| per-drug gap (points), p | 4.2, 4.5e-5 | 2.9, 0.047 |
| folds real better (per drug) | 5/5 | 5/5 |
| drugs better / worse | 55% / 45% | 53% / 46% |

Same direction in every fold of both seeds, but the per-drug gap is 2.9–4.2 points and the p-value spans three orders
of magnitude between seeds. Report 02's p = 4.5e-5 was at the favourable end. The seed-1 control also sits below 0.5
(0.488 per perturbation), which narrows its gap. `rank_probc` replicates (gap 6.8 → 3.5, p 5.7e-10 → 1.5e-3); plain
per-drug `rank_prob` does not (1.1 points p = 9.4e-4 → 0.5 points p = 0.09).

## 4. MetChange, expression and TF side by side

Per-drug gap (points) and Wilcoxon p, real vs own shuffled control:

| metric | MetChange s0 | MetChange s1 | expr | TF (control) |
|---|---|---|---|---|
| `spec_prob_timedose` | 4.2, 4.5e-5 | 2.9, 0.047 | 4.1, 1.7e-3 | 6.8, 1.2e-8 |
| `rank_probc` | 6.8, 5.7e-10 | 3.5, 1.5e-3 | 6.4, 4.0e-6 | 9.5, 4.2e-15 |
| `rank_prob` | 1.1, 9.4e-4 | 0.5, 0.09 | 1.1, 0.19 | 3.3, 5.0e-8 |

Medians, real / control:

| | MetChange s0 | MetChange s1 | expr | TF |
|---|---|---|---|---|
| `spec_prob_timedose` per perturbation | 0.414 / 0.500 | 0.425 / 0.488 | 0.417 / 0.500 | 0.341 / 0.488 |
| `spec_prob_timedose` per drug | 0.444 / 0.500 | 0.455 / 0.493 | 0.461 / 0.496 | 0.417 / 0.491 |
| `rank_prob` per perturbation | 0.262 / 0.274 | 0.271 / 0.284 | 0.391 / 0.397 | 0.247 / 0.372 |
| `rank_prob` per drug | 0.401 / 0.410 | 0.404 / 0.410 | 0.454 / 0.459 | 0.405 / 0.443 |
| own drug in top 10% (`rank_prob`, per perturbation) | 23% / 18% | 23% / 20% | 15% / 13% | 24% / 14% |
| folds real better, `spec_prob_timedose` per drug | 5/5 | 5/5 | 3/5 | 5/5 |
| FNN members keeping epoch ≤ 2 (of 15) | 2 | 1 | 8 | 7 |

Reading:

1. **The positive control works.** TF has the largest gap on every metric, in 5/5 folds, and is the only signature
   where the signature clearly moves *plain* per-drug ranking (`rank_prob` 3.3 points, p = 5e-8). The pipeline does
   what the published setting says it should.
2. **MetChange reads the signature at roughly half to two-thirds of the TF level** on `spec_prob_timedose` (2.9–4.2 vs
   6.8 points) and `rank_probc`, with the cleanest training of the three (1–2 collapsed members of 15).
3. **The gap to TF is largest on plain ranking.** MetChange's real and shuffled runs both rank well (0.262 / 0.274 per
   perturbation): mostly prior. In TF the control ranks much worse (0.372) and the signature does the work. Caveat:
   TF uses a different network (`networks/tf`), so its prior differs; compare gaps, not absolute `rank_prob`.
4. **Expression is training-limited, not signal-limited.** Per fold, `spec_prob_timedose` per perturbation is
   0.484 / 0.500 / 0.328 / 0.452 / 0.309 (control ~0.50). Folds 2 and 4 are as strong as TF; in fold 1 all three
   members kept epoch 1–2 and `rank_prob` is 0.689, worse than chance.
5. Even for TF the real run beats its control on only 57% of drugs (MetChange 53–55%): the effect is a consistent shift
   across drugs, not most drugs being found.

## 5. Why so many members keep epoch 1

Validation AUPRC is very noisy early on. A typical collapsed member (expr, fold 1): 0.131 at epoch 1, then
0.071 0.061 0.091 0.041 … 0.051 by epoch 11, stopped. Members that survive show the same early dip and then climb
past 0.2 around epochs 10–14, after the first learning-rate halving (e.g. MetChange fold 1, member 1: 0.074 … 0.218 at
epoch 11, kept epoch 47). With patience 10 the run stops just before that recovery, so the noisy epoch-1 peak wins.

The ensemble then averages a 1-epoch member (close to the prior) with converged ones, which dilutes the signature
effect. This affects TF and expr most, so TF's numbers are probably an underestimate of the positive control too.

**Fix:** `--min-epochs N` (new, default 1 = previous behaviour, checked by replaying the 150 recorded member histories
of all AUPRC-selected runs: kept epochs and stopping points reproduced exactly). Epochs before N are trained and
logged but can neither be kept nor stop training; patience counts from N.

## 6. Caveats

- One seed for TF and expr; MetChange's two seeds already differ by 1.3 points of per-drug gap.
- TF's own `rank_prob` is not comparable in absolute terms (different network and prior).
- The per-drug value is a mean over a drug's perturbations, which stays near 0.5 unless they agree; per-perturbation
  gaps are larger but count doses of one drug as independent.

## 7. Next steps

1. The early-stopping and schedule runs of section 8 (`_min15`, `_min15_pat25`, `_plateau`). Check that members no
   longer keep epoch 1–2, that expr folds 0–1 recover, and whether MetChange's gap to TF narrows or TF's grows.
2. If expr recovers, it becomes the comparison that matters for the paper: MetChange vs the expression it is computed
   from, both against their own controls.
3. A second seed for whichever recipe wins, before building on it (section 3 shows the seed-to-seed spread is about as
   large as the differences between signature types).
4. The prototype path of `STATUS.md` section 6.2 (query path, signature-similarity baseline, Schuhknecht benchmark) is
   unchanged and does not need the GPU.

## 8. Follow-up runs (launched 2026-11-14/15; results in section 9)

### 8.1 What the kept epochs say about the schedule

`train_ffn` used `StepLR(step_size=10, gamma=0.5)`: the lr halves every 10 epochs. Over the four MetChange runs of this
report (seeds 0 and 1, real and control; 60 FNN members):

- **29/60 kept the first epoch after a halving** (11, 21, 31, 41 or 51; chance ~5/60). Validation AUPRC improves in a
  jump when the lr drops and not in the following nine epochs.
- **52/60 stopped by patience**, not at the cap. With patience 10 = the step size, a member whose best epoch is 21
  stops at 31, the first epoch of the next lr level, if that epoch does not beat it.
- Only 3/15 seed-0 real members (1/15 seed 1) kept an epoch before 15, so `--min-epochs 15` barely affects MetChange;
  it is aimed at TF / expr (7/15 and 8/15 members at epoch 1–2).

So the schedule, not convergence, decides when MetChange training stops. A larger step size would keep the lr high
(and validation noisy) for longer; a larger gamma makes each jump smaller. Neither addresses that coupling, hence two
other tests.

### 8.2 Code additions (`metabolic_target_prediction.py`)

All defaults reproduce the earlier runs.

| option | default | effect |
|---|---|---|
| `--min-epochs N` | 1 | epochs < N are trained and logged but can neither be kept nor stop training; `--patience` counts from N |
| `--lr-schedule {step,plateau}` | `step` | `step`: `StepLR`; `plateau`: `ReduceLROnPlateau` on the validation selection score, stepped from epoch 1 (independent of `--min-epochs`) |
| `--lr-step-size`, `--lr-gamma` | 10, 0.5 | `StepLR` step and factor; `--lr-gamma` is also the plateau factor |
| `--plateau-patience` | 3 | epochs without improvement before a reduction (the lr drops after the 4th); must be < `--patience` |
| `--min-lr` | 1e-6 | floor of the plateau schedule |

`metrics.json` now records each epoch's lr in `history`, and `min_epochs` / `lr_schedule` in `epoch_selection`.

Checks (no torch on the laptop): with the defaults, replaying the recorded validation histories of all 150
AUPRC-selected members on the laptop reproduces every kept epoch and stopping point; a stub-torch test of `train_ffn`
confirms the scheduler calls (`StepLR` halving after epoch 10 by default, plateau reduction after 4 non-improving
epochs from epoch 1, `--min-epochs` ignoring the epoch-1 peak, custom step size / gamma passed through). Not yet run
end to end with real torch: the first cluster task of each `_plateau` run is the first real test.

### 8.3 Runs

All seed 0 (same folds and splits as `*_noGlobals_lr1e-3`), `--no-globals --lr 0.001`, array jobs, each with its
shuffled control and comparison.

| run (`results/predict_targets/A549/`) | schedule | early stopping | wrapper |
|---|---|---|---|
| `metchange_rank_noGlobals_lr1e-3_min15` | StepLR | min 15, patience 15, cap 80 | `submit_train_w_shuff_noGlobals_lr1e-3_min15.sh` |
| `expr_noGlobals_lr1e-3_min15` | StepLR | same | same |
| `tf_active200_noGlobals_lr1e-3_min15` | StepLR | same | same |
| `metchange_rank_noGlobals_lr1e-3_min15_pat25` | StepLR | min 15, **patience 25**, cap 80 | `submit_train_w_shuff_mc_noGlobals_lr1e-3_pat25_plateau.sh` |
| `metchange_rank_noGlobals_lr1e-3_plateau` | **plateau** (patience 3) | min 15, patience 20, cap 80 | same |
| `expr_noGlobals_lr1e-3_plateau` | **plateau** | same | `submit_train_w_shuff_tf_expr_noGlobals_lr1e-3_plateau.sh` |
| `tf_active200_noGlobals_lr1e-3_plateau` | **plateau** | same | same |

What each answers:

- **`_min15` (TF, expr):** does a warm-up stop members keeping epoch 1–2, and do expr folds 0–1 recover?
  For MetChange it mainly tests patience 15 vs 10 (expected: little change, see 8.1).
- **`_min15_pat25` (MetChange):** does seeing two more halvings after the best epoch help? If members keep later
  epochs and the per-drug gap grows, the schedule was limiting; if not, MetChange training is done and its gap to TF
  is the signature.
- **`_plateau` (all three):** does tying lr reductions to validation do better than a fixed schedule? Compare each with
  its `_min15` run (same warm-up, StepLR).

### 8.4 What to read from the results

For each run, against its own control: per-drug `median_diff` and p for `spec_prob_timedose`, `rank_probc`,
`rank_prob`; folds better; and from `metrics.json`: kept epochs, history length (patience vs cap), and the lr at the
kept epoch. Compare with `*_noGlobals_lr1e-3` (same seed and splits, so fold-by-fold comparisons are valid within
seed 0). Any winner still needs a second seed (section 3: the seed spread is ~1.3 points of per-drug gap).

## 9. Results of the follow-up runs (2026-11-15)

All seed 0, same folds as `*_noGlobals_lr1e-3` ("base"). Gap = per-drug `median_diff` in points, positive = real
better than its own control; p = Wilcoxon per drug; folds = folds where real beats control per drug.

| run | `spec_prob_timedose` gap, p, folds | `rank_probc` gap, p | `rank_prob` gap, p, folds | `rank_prob` per drug, real / control |
|---|---|---|---|---|
| MetChange base | 4.2, 4e-5, 5/5 | 6.8, 6e-10 | 1.1, 9e-4, 4/5 | 0.401 / 0.410 |
| MetChange `_min15` | 4.2, 1e-4, 5/5 | 6.4, 7e-9 | 1.3, 1e-3, 4/5 | 0.396 / 0.407 |
| MetChange `_min15_pat25` | 3.6, 3e-4, 5/5 | 5.4, 4e-8 | 1.3, 1e-3, 4/5 | 0.392 / 0.411 |
| MetChange `_plateau` | 3.9, 2e-4, 5/5 | 6.4, 5e-8 | 1.5, 5e-4, 4/5 | 0.396 / 0.410 |
| expr base | 4.1, 2e-3, 3/5 | 6.4, 4e-6 | 1.1, 0.19, 3/5 | 0.454 / 0.459 |
| **expr `_min15`** | **7.8, 2e-11, 5/5** | **8.0, 4e-16** | **3.3, 5e-15, 5/5** | **0.381 / 0.437** |
| expr `_plateau` | 7.4, 7e-12, 5/5 | 8.0, 8e-16 | 3.1, 1e-12, 5/5 | 0.386 / 0.439 |
| TF base | 6.8, 1e-8, 5/5 | 9.5, 4e-15 | 3.3, 5e-8, 4/5 | 0.405 / 0.443 |
| TF `_min15` | 9.0, 8e-15, 5/5 | 10.4, 9e-20 | 3.3, 2e-14, 5/5 | 0.397 / 0.445 |
| TF `_plateau` | 8.7, 5e-12, 5/5 | 9.7, 8e-17 | 3.5, 2e-17, 5/5 | 0.388 / 0.445 |

Training (real runs; 15 FNN members each):

| run | kept epoch ≤ 2 | median kept epoch (range) | members at the cap | lr at kept epoch (median) |
|---|---|---|---|---|
| MetChange base | 2 | 21 (1–53) | 1 | — |
| MetChange `_min15` | 0 | 40 (21–61) | 0 | 1.3e-4 |
| MetChange `_min15_pat25` | 0 | 47 (21–67) | 5 | 6.3e-5 |
| MetChange `_plateau` | 0 | 36 (23–55) | 0 | 3.1e-5 |
| expr base | 8 | 2 (1–58) | 2 | — |
| expr `_min15` | 0 | 66 (43–80) | 8 | 1.6e-5 |
| expr `_plateau` | 0 | 62 (43–79) | 10 | 3.9e-6 |
| TF base | 7 | 12 (1–58) | 4 | — |
| TF `_min15` | 0 | 76 (21–80) | 12 | 7.8e-6 |
| TF `_plateau` | 0 | 66 (45–80) | 9 | 1.0e-6 |

### 9.1 Reading

1. **The warm-up fixes the collapse.** No member keeps epoch 1–2 in any `_min15` or `_plateau` run (base: 7–8/15 for
   TF and expr).
2. **Expression was training-limited, and now matches the positive control.** Per fold, `spec_prob_timedose` per
   perturbation goes from 0.484 / 0.500 / 0.328 / 0.452 / 0.309 to 0.295 / 0.319 / 0.316 / 0.308 / 0.301 (controls
   0.47–0.50). Its per-drug gap (7.8) is close to TF's (9.0), and its per-drug `rank_prob` (0.381) is the best of any
   run: the signature now moves plain ranking (3.3 points, p = 5e-15, 5/5 folds), as in TF.
3. **TF gains too** (`spec_prob_timedose` gap 6.8 → 9.0), so the positive control was also underestimated.
4. **MetChange does not move.** Warm-up, patience 25 and plateau all leave the per-drug gap at 3.6–4.2 points
   (base 4.2) and `rank_prob` at 1.1–1.5 points. Kept epochs move later (median 21 → 36–47), as intended, but the
   later epochs only refine what the base run had. **MetChange training is done**; the remaining gap to TF and expr
   is the signature, not the optimisation.
5. **Plateau ≈ StepLR with warm-up**, for all three signatures. On the cluster it reduced the lr every ~5–6
   epochs (e.g. MetChange fold 0: epochs 10, 15, 21, 27, 32, …), i.e. a faster decay than StepLR, down to the 1e-6
   floor. No reason to prefer it; `_min15` (StepLR) stays the simpler recipe.
6. **The cap now binds for TF and expr** (8–12 of 15 members at epoch 80), but at an lr of ~1e-5 or less, where
   further epochs change little. Not worth chasing before other questions.

### 9.2 What this means

- **New reference recipe:** `--no-globals --lr 0.001 --min-epochs 15 --patience 15 --ffn-epochs 80` (StepLR).
- **The central comparison is now MetChange vs expr:** MetChange is computed from the same expression, yet it reads
  the signature about half as well (4.2 vs 7.8 points per drug on `spec_prob_timedose`) and barely moves plain ranking
  (1.3 vs 3.3). Under this setup, the MetChange transformation loses drug-specific information that the expression
  carries. Candidate causes, none tested: the representation (`--ignore-qvalue --top-n 100`, then `np.sign`: 100
  discretised reactions), the network (MetChange uses `networks/combined`, expr `networks/gene`), and what MetChange
  computes (inconsistency over the GEM's GPR genes, many of them CycleGAN-inferred, `STATUS.md` 6.3).
- **Caveats:** one seed for every run in this section (seed spread ~1.3 points, section 3). The MetChange `_min15`
  control is below 0.5 in folds 0–1 (0.453, 0.468 per perturbation), not seen in the base control; per drug the
  control is 0.504. Worth a look if it recurs.

### 9.3 Next steps

1. **Second seed of the new recipe** for MetChange, expr and TF, each with its control, to confirm the
   MetChange-vs-expr gap is larger than the seed spread.
2. **Locate where MetChange loses information** (no GPU first): the representation (top-n, signs vs continuous
   values: the `--signature-weights` option proposed in `STATUS.md` section 7), expr on the `combined` network, and the
   landmark vs inferred share of the GPR genes (`STATUS.md` 6.3).
3. The prototype path (`STATUS.md` 6.2) with expr as a second reference alongside MetChange.
