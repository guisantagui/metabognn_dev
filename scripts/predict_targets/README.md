# Intructions for launching training scripts

## Training scripts with shuffle controls
This script launches simultaneosly a training script working with the real labels, and a job with shuffled labels. Then, it calls compare_runs.py to obtain a comparison of the results. It has different working modes.


```bash
bash submit_with_shuffle_control.sh tf                   # Submits real + shuffled + comparison (on tf data)
bash submit_with_shuffle_control.sh expr --real-job <id> # Submits real shuffle + comparison (on expression data). In case the real job is already running.
bash submit_with_shuffle_control.sh mc   --real-job <id> # Submits real shuffle + comparison (on MetChange data). In case the real job is already running.
```

(a) Plain — same loss, different sampling. Fold the matched rows in as ordinary negatives. Rows stay independent; only the training distribution changes. What it buys you is decorrelation: today a protein like EGFR appears almost only as a positive, so protein identity predicts the label. Per chunk you have roughly

┌─────────────────────────────────────────────┬───────────┐
│                                             │   rows    │
├─────────────────────────────────────────────┼───────────┤
│ positives                                   │ ~100–300  │
├─────────────────────────────────────────────┼───────────┤
│ random negatives (--neg-per-pert 200, P=32) │ 6400      │
├─────────────────────────────────────────────┼───────────┤
│ matched negatives at K=4                    │ ~400–1200 │
└─────────────────────────────────────────────┴───────────┘

After adding them, P(target | protein = EGFR) in the training set drops from ~1 to ~1/5. That's a substantial change for a small number of rows — but it also means they're easily swamped, so they need their own sample weight (--matched-weight), not just to be dumped into class 0.

(b) Grouped — the matrix version you're describing. Keep each group {(p,g), (q₁,g), …, (q_K,g)} together as a unit and put a loss on the group: a softmax over the group (which of these K+1 perturbations targets g?) or a pairwise ranking loss on the positive against each matched negative. Here the batch genuinely becomes (n_groups, K+1, F) instead of (n_rows, F).
