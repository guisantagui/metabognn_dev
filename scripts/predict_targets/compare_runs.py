#!/usr/bin/env python
"""
Compare metabolic_target_prediction.py cross-validation runs side by side, typically a run against its
--shuffle-signatures control (same drugs in the same folds).

  python compare_runs.py --run real=<out_dir> --run shuffled=<out_dir_shuff> --control shuffled

Reads <out_dir>/cv_drug_ranks.tsv and cv_drug_ranks_per_drug.tsv, or, while a run is not finished, the
test_drug_ranks*.tsv of its finished folds (only the folds finished in every run are compared). By default only the
active perturbations are compared (--subset all for every scored one), so runs with and without --score-active-only
are comparable. Writes to --out-dir (default: <first run>/compare_<names>):
  summary.tsv      per run, ranking and level: n, median percentile, top1/5/10/25_pct (held-out ranking: chance 0.5,
                   topX_pct X / 100)
  per_fold.tsv     the same per fold (median percentile, top10_pct)
  paired.tsv       each run against --control, per drug (drugs are the independent units): median of the paired
                   percentile differences, fraction of drugs ranked better, Wilcoxon signed-rank p-value
  ecdf_<rank>.png  cumulative fraction of drugs vs percentile, one line per run (diagonal: chance)
spec_*_heldout (does a perturbation score its own drug above how the other drugs' perturbations score it) and
rank_*c_heldout (held-out ranking after removing each drug's average score) are blind to a per-drug prior: 0.5 for
a model that ignores the signature.

They are not blind to dose / time: if the experimental design (time point, dose series) differs between drug classes,
a model can recognise drugs through it without the signature. When the runs have fold*/test_scores_heldout.npz,
spec is therefore also recomputed within strata: spec_<prob|logreg>_time_heldout compares each perturbation only with
other drugs' perturbations at the same time point, spec_<prob|logreg>_timedose_heldout at the same time and relative
dose (log10 of dose / the drug's highest dose, from the '<drug>_<dose>_<time>' perturbation names). Per drug: the
mean over the drug's perturbations. These are added to the comparison automatically.
"""
import argparse
import os

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--run", action="append", required=True, metavar="NAME=DIR",
                help="a metabolic_target_prediction.py --out-dir (repeat; at least 2)")
ap.add_argument("--control", help="name of the run the others are tested against (default: the last --run)")
ap.add_argument("--subset", choices=["active", "all"], default="active")
ap.add_argument("--ranks", nargs="+",
                default=["rank_prob_heldout", "rank_bes_heldout", "rank_logreg_heldout", "rank_baseline_heldout",
                         "spec_prob_heldout", "spec_logreg_heldout", "rank_probc_heldout", "rank_logregc_heldout"],
                help="rank columns to compare (those missing from a run are skipped; spec_* and rank_*c are the "
                     "prior-free metrics, only in runs made after they were added)")
ap.add_argument("--out-dir")
args = ap.parse_args()

runs = dict(r.split("=", 1) for r in args.run)
if len(runs) < 2:
    ap.error("give at least two --run NAME=DIR")
control = args.control or list(runs)[-1]
if control not in runs:
    ap.error(f"--control {control} is not one of the runs: {list(runs)}")
names = list(runs)


def load(d):
    """(per perturbation, per drug) out-of-fold ranks, with a fold column."""
    pp, pdr = os.path.join(d, "cv_drug_ranks.tsv"), os.path.join(d, "cv_drug_ranks_per_drug.tsv")
    if os.path.exists(pp) and os.path.exists(pdr):
        return pd.read_csv(pp, sep="\t"), pd.read_csv(pdr, sep="\t")
    folds = sorted(int(f[4:]) for f in os.listdir(d) if f.startswith("fold") and f[4:].isdigit()
                   and os.path.exists(os.path.join(d, f, "metrics.json")))
    if not folds:
        raise SystemExit(f"{d}: no cv_drug_ranks.tsv and no finished fold")
    read = lambda f, name: pd.read_csv(os.path.join(d, f"fold{f}", name), sep="\t").assign(fold=f)
    return (pd.concat([read(f, "test_drug_ranks.tsv") for f in folds], ignore_index=True),
            pd.concat([read(f, "test_drug_ranks_per_drug.tsv") for f in folds], ignore_index=True))


STRATA = {"time": ["time"], "timedose": ["time", "rel_dose"]}
PCT_COLS = set()  # columns that already hold percentiles (the stratified spec_*)
UNIT_UM = {"pm": 1e-6, "nm": 1e-3, "um": 1.0, "µm": 1.0, "mm": 1e3, "m": 1e6}


def design(pp):
    """Time and relative dose of each perturbation, from LINCS names '<drug>_<dose>_<time>' (NaN if unparsable).
    Relative dose: log10(dose / the drug's highest scored dose), rounded to 0.1, as the model's rel_log_dose."""
    parts = pp.pert.str.rsplit("_", n=2, expand=True).reindex(columns=[0, 1, 2])
    m = parts[1].str.extract(r"([\d.]+)\s*([a-zA-Zµ]+)")
    dose = pd.to_numeric(m[0], errors="coerce") * m[1].str.lower().map(UNIT_UM)
    rel = np.log10(dose / dose.groupby([pp.fold, pp.drug.str.lower()]).transform("max")).round(1)
    return pd.DataFrame({"time": parts[2], "rel_dose": rel}, index=pp.index)


def stratified_spec(d, pp, min_ref=5):
    """spec_<model>_<stratum>_heldout per perturbation: spec_*_heldout recomputed from fold*/test_scores_heldout.npz
    with each perturbation compared only to other drugs' perturbations of the same activity class AND the same
    stratum (time, or time x relative dose), so dose / time design cannot carry drug identity. Percentile (lower is
    better, chance 0.5); NaN with fewer than min_ref perturbations to compare with. Runs without the score files
    get no columns."""
    if any(not os.path.exists(os.path.join(d, f"fold{f}", "test_scores_heldout.npz")) for f in pp.fold.unique()):
        print(f"{d}: no test_scores_heldout.npz in every fold, stratified spec skipped")
        return pp
    pp = pp.copy()
    des = design(pp)
    active = pp.active.astype(str).str.upper() == "TRUE" if "active" in pp else pd.Series(True, index=pp.index)
    for f in pp.fold.unique():
        z = np.load(os.path.join(d, f"fold{f}", "test_scores_heldout.npz"), allow_pickle=True)
        sub = pp.index[pp.fold == f]
        row = pd.Index(z["pert"]).get_indexer(pp.pert[sub])
        if (row < 0).any():
            raise SystemExit(f"{d} fold {f}: perturbations of test_drug_ranks.tsv missing from the score file")
        col = {k: j for j, k in enumerate(z["drug"])}
        drugs = pp.drug[sub].str.lower().to_numpy()
        c = np.array([col[k] for k in drugs])
        for model in ("prob", "logreg"):
            if model not in z.files or not len(z[model]):
                continue
            S = z[model][row]
            own = S[np.arange(len(S)), c]
            for name, cols in STRATA.items():
                key = (active[sub].astype(str) + "|" + des.loc[sub, cols].astype(str).agg("|".join, axis=1)).to_numpy()
                ok = des.loc[sub, cols].notna().all(axis=1).to_numpy()
                p = np.full(len(S), np.nan)
                for i in np.where(ok)[0]:
                    ref = S[(drugs != drugs[i]) & (key == key[i]) & ok, c[i]]
                    if len(ref) >= min_ref:  # mid-rank percentile: exactly 0.5 on average under chance, whatever
                        # len(ref) (rank / (len + 1) averages (1 + len / 2) / (len + 1), 0.58 for len 5)
                        p[i] = ((ref > own[i]).sum() + 0.5 * (ref == own[i]).sum() + 0.5) / (len(ref) + 1)
                name_col = f"spec_{model}_{name}_heldout"
                pp.loc[sub, name_col] = p
                PCT_COLS.add(name_col)
    return pp


data = {}
for n, d in runs.items():
    pp, pdr = load(d)
    pp = stratified_spec(d, pp)
    if args.subset == "active":
        if "active" not in pp:
            raise SystemExit(f"{n}: no activity information (run without --de-counts); use --subset all")
        pp = pp[pp.active.astype(str).str.upper() == "TRUE"]
        pdr = pdr[pdr.subset == "active"]
    elif "subset" in pdr:
        pdr = pdr[pdr.subset == "all"]
    pdr = pdr.copy()
    for c in PCT_COLS & set(pp):  # per drug: mean over its (active) perturbations
        m = pp.groupby([pp.fold, pp.drug.str.lower()])[c].mean()
        pdr[c] = pd.MultiIndex.from_arrays([pdr.fold, pdr.drug.str.lower()]).map(m.to_dict().get)
    data[n] = {"perturbation": pp, "drug": pdr}
common = set.intersection(*(set(v["drug"].fold) for v in data.values()))
for n, v in data.items():
    skipped = sorted(set(v["drug"].fold) - common)
    if skipped:
        print(f"{n}: folds {skipped} not finished in every run, left out")
    for lvl in v:
        v[lvl] = v[lvl][v[lvl].fold.isin(common)].copy()
print(f"comparing {len(common)} fold(s) {sorted(common)}, {args.subset} perturbations; control: {control}")
ranks = [r for r in args.ranks + sorted(PCT_COLS - set(args.ranks)) if any(r in data[n]["drug"] for n in names)]


def pct(df, r):
    if r in PCT_COLS:
        return df[r].dropna().astype(float)
    s = r[r.index("_heldout"):] if r.endswith("_heldout") else ""
    p = df[r] / df[("n_spec" if r.startswith("spec_") else "n_candidates") + s]
    return p.dropna()  # spec_*: NaN when no unit of another drug to compare with


def summary(p):
    return {"n": len(p), "median_percentile": p.median(),
            **{f"top{q}_pct": (p <= q / 100).mean() for q in (1, 5, 10, 25)}}


rows, fold_rows = [], []
for lvl in ("perturbation", "drug"):
    for r in ranks:
        for n in names:
            df = data[n][lvl]
            if r not in df:
                continue
            p = pct(df, r)
            rows.append(dict(level=lvl, rank=r, run=n, **summary(p)))
            for f, pf in p.groupby(df.fold):
                fold_rows.append(dict(level=lvl, rank=r, run=n, fold=f, n=len(pf), median_percentile=pf.median(),
                                      top10_pct=(pf <= 0.1).mean()))
summ, per_fold = pd.DataFrame(rows), pd.DataFrame(fold_rows)

paired_rows = []
try:
    from scipy.stats import wilcoxon
except ImportError:
    wilcoxon = None
for r in ranks:
    if r not in data[control]["drug"]:
        continue
    c = data[control]["drug"].assign(p=lambda d: pct(d, r)).set_index("drug").p.dropna()
    for n in names:
        if n == control or r not in data[n]["drug"]:
            continue
        x = data[n]["drug"].assign(p=lambda d: pct(d, r)).set_index("drug").p.dropna()
        both = x.index.intersection(c.index)
        diff = (x[both] - c[both]).to_numpy()  # negative: ranked higher than in the control
        pval = float(wilcoxon(diff).pvalue) if wilcoxon is not None and np.any(diff != 0) else np.nan
        paired_rows.append(dict(rank=r, run=n, control=control, n_drugs=len(both),
                                median_pct_run=x[both].median(), median_pct_control=c[both].median(),
                                median_diff=np.median(diff), frac_drugs_better=(diff < 0).mean(),
                                frac_drugs_worse=(diff > 0).mean(), wilcoxon_p=pval))
paired = pd.DataFrame(paired_rows)

out_dir = args.out_dir or os.path.join(runs[names[0]], "compare_" + "_vs_".join(names))
os.makedirs(out_dir, exist_ok=True)
summ.to_csv(os.path.join(out_dir, "summary.tsv"), sep="\t", index=False, float_format="%.4g")
per_fold.to_csv(os.path.join(out_dir, "per_fold.tsv"), sep="\t", index=False, float_format="%.4g")
paired.to_csv(os.path.join(out_dir, "paired.tsv"), sep="\t", index=False, float_format="%.4g")

pd.set_option("display.width", 200)
for lvl in ("perturbation", "drug"):
    t = summ[summ.level == lvl].pivot_table(index=["rank", "run"], values=["n", "median_percentile", "top1_pct",
                                            "top10_pct", "top25_pct"], sort=False)
    print(f"\n== per {lvl} (held-out ranking: chance median 0.5, topX_pct X/100)")
    print(t[["n", "median_percentile", "top1_pct", "top10_pct", "top25_pct"]].round(3).to_string())
print("\n== per fold, median percentile (top10_pct)")
pf = per_fold[per_fold.level == "drug"].assign(
    v=lambda d: d.median_percentile.map("{:.3f}".format) + " (" + d.top10_pct.map("{:.2f}".format) + ")")
print("per drug\n" + pf.pivot_table(index=["rank", "run"], columns="fold", values="v", aggfunc="first",
                                    sort=False).to_string())
if len(paired):
    print(f"\n== paired per drug against {control} (median_diff < 0 and small p: the run ranks drugs higher)")
    print(paired.drop(columns="control").round(4).to_string(index=False))

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"]
    for r in ranks:
        fig, ax = plt.subplots(figsize=(4.2, 4.0))
        ax.plot([0, 1], [0, 1], color="#BBBBBB", lw=1, ls="--", label="chance")
        for i, n in enumerate(names):
            if r not in data[n]["drug"]:
                continue
            p = np.sort(pct(data[n]["drug"], r).to_numpy())
            ax.step(p, np.arange(1, len(p) + 1) / len(p), where="post", lw=2,
                    color="#8C8C8C" if n == control else colors[i % len(colors)],
                    label=f"{n} (median {np.median(p):.2f})")
        ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="percentile of the drug's rank (0 = top)",
               ylabel="cumulative fraction of drugs", title=f"{r}, per drug ({args.subset})")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=8, loc="lower right")
        fig.savefig(os.path.join(out_dir, f"ecdf_{r}.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)
except ImportError:
    print("matplotlib not available: no plots")
print(f"\nsaved {out_dir}/summary.tsv, per_fold.tsv, paired.tsv, ecdf_*.png")
