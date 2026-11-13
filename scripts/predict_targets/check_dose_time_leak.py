#!/usr/bin/env python
"""Usage: python check_dose_time_leak.py <results/predict_targets/A549> <run> [<run> ...]

Does dose / time carry drug identity into the held-out scores? For each run and fold (test_scores_heldout.npz):
  spec(model)      median per-perturbation spec percentile of the model's scores (chance 0.5)
  spec(dose/time)  the same after replacing each perturbation's score row by the mean row of OTHER drugs'
                   perturbations with the same (time, relative dose): a predictor that knows only dose and time
  R2               share of each candidate drug's score variance (over test perturbations) explained by
                   (time, relative dose) groups, averaged over candidates"""
import sys
import numpy as np, pandas as pd

UNIT = {"pm": 1e-6, "nm": 1e-3, "um": 1.0, "µm": 1.0, "mm": 1e3, "m": 1e6}

def parse(perts):
    d = pd.Series(perts).str.rsplit("_", n=2, expand=True)
    d.columns = ["drug", "dose", "time"]
    m = d.dose.str.extract(r"([\d.]+)\s*([a-zA-Zµ]+)")
    d["drug"] = d.drug.str.lower()
    d["dose_um"] = pd.to_numeric(m[0], errors="coerce") * m[1].str.lower().map(UNIT)
    return d

def spec_pct(S, drugs, col):
    r = np.empty(len(S))
    for i in range(len(S)):
        c = col[drugs[i]]; ref = S[drugs != drugs[i], c]
        r[i] = (1 + (ref > S[i, c]).sum() + 0.5 * (ref == S[i, c]).sum()) / (len(ref) + 1)
    return r

root = sys.argv[1]
for run in sys.argv[2:]:
    rows = []
    for f in range(5):
        z = np.load(f"{root}/{run}/fold{f}/test_scores_heldout.npz", allow_pickle=True)
        S, col = z["prob"], {d: j for j, d in enumerate(z["drug"])}
        g = parse(z["pert"])
        # relative dose as in pert_globals, but against the drug's max dose within the scored perturbations
        g["rel"] = np.round(np.log10(g.dose_um / g.groupby("drug").dose_um.transform("max")), 1)
        drugs, key = g.drug.to_numpy(), (g.time + "|" + g.rel.astype(str)).to_numpy()
        Sg = np.empty_like(S)
        for i in range(len(S)):
            m = (key == key[i]) & (drugs != drugs[i])
            Sg[i] = S[m].mean(0) if m.any() else S[drugs != drugs[i]].mean(0)
        r2 = np.mean([1 - pd.Series(S[:, c]).groupby(key).transform(lambda x: x - x.mean()).var() / S[:, c].var()
                      for c in range(S.shape[1]) if S[:, c].var() > 0])
        rows.append((np.median(spec_pct(S, drugs, col)), np.median(spec_pct(Sg, drugs, col)), r2, len(set(key))))
    o = np.array(rows)
    print(f"{run:24s} spec(model) {o[:, 0].round(3)} | spec(dose/time) {o[:, 1].round(3)} | "
          f"R2 dose/time {o[:, 2].round(2)} | groups {o[:, 3].astype(int)}")
