"""
Significant features per perturbation (q <= --qvalue in limma_DE.R qvals.csv) against the design of each fit
(fit_info.csv: treated samples, controls, batches) and against dose and time parsed from the perturbation name
('<drug>_<dose> <unit>_<time> h').
Outputs in --dir (or --out-dir): fit_info_nsig.csv and n_sig_diagnostics.png."""
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter
)
ap.add_argument("--dir", required=True, help="folder with qvals.csv and fit_info.csv")
ap.add_argument("--qvalue", type=float, default=0.05)
ap.add_argument("--out-dir", default=None, help="default: --dir")
args = ap.parse_args()
d = args.dir
out_dir = args.out_dir or d
os.makedirs(out_dir, exist_ok=True)
# ---- counts (qvals.csv read in row blocks: features x perturbations)
n, n_feat = None, 0
for chunk in pd.read_csv(os.path.join(d, "qvals.csv"), index_col=0, chunksize=200):
    c = (chunk <= args.qvalue).sum(axis=0)
    n = c if n is None else n + c
    n_feat += len(chunk)
info = pd.read_csv(os.path.join(d, "fit_info.csv"), index_col=0)
info["n_treated"] = info.n_samples - info.n_ctrl
info["n_sig"] = n.reindex(info.index).fillna(0).astype(int)
info["frac_sig"] = info.n_sig / n_feat
parts = info.index.to_series().str.extract(
        r"^(?P<drug>.+)_(?P<dose>[\d.eE+-]+)?(?P<unit>[^_]*)_(?P<time>[\d.]+) ?h$"
)
info["drug"] = parts.drug
info["dose"] = pd.to_numeric(parts.dose, errors="coerce")
info["dose_unit"] = parts.unit
info["time_h"] = pd.to_numeric(parts.time, errors="coerce")
info.to_csv(os.path.join(out_dir, "fit_info_nsig.csv"))

print(
    f"{n_feat:,} features x {len(info):,} perturbations, q <= {args.qvalue}"
)

print(
    f"n_sig = 0: {(info.n_sig == 0).mean():.1%}; > half of the features: {(info.frac_sig > 0.5).mean():.1%}"
)
print(
    info[["n_treated", "n_ctrl", "n_batches", "dose", "time_h","n_sig"]].corr(method="spearman").round(2)
)
ctrl_bins = pd.cut(info.n_ctrl, [0, 20, 50, 100, 300, 1000, np.inf])
print(info.groupby(ctrl_bins, observed=True).n_sig.describe())
print(info.groupby("time_h", observed=True).n_sig.describe())

# ---- figure
INK, MUTED, MARK = "#1f2328", "#6e7781", "#2f6db5"
plt.rcParams.update(
    {
        "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
        "axes.spines.right": False
    }
)
fig, ax = plt.subplots(2, 3, figsize=(13, 7.5))


def scatter(a, x, xlabel):
    ok = info[x] > 0
    a.scatter(
        info.loc[ok, x],
        info.loc[ok, "n_sig"],
        s=4, alpha=0.15,
        color=MARK,
        linewidths=0,
        rasterized=True
    )
    a.set_xscale("log")
    a.set_xlabel(xlabel)
    a.set_ylabel("significant features")


def frac_bars(a, groups, xlabel):
    g = info.groupby(groups, observed=True)
    stats = pd.DataFrame(
        {"any": g.n_sig.apply(lambda s: (s > 0).mean()),
                          "most": g.frac_sig.apply(lambda s: (s > 0.5).mean()), "n": g.size()
        }
    )
    x = np.arange(len(stats))
    a.bar(x - 0.2, stats["any"], 0.38, color=MARK, label="any significant")
    a.bar(x + 0.2, stats["most"], 0.38, color=MARK, alpha=0.45, label="> half significant")
    lab = [
        f"{i.left:g}-{i.right:g}" if isinstance(i, pd.Interval) else f"{i:g}" for i in stats.index
    ]
    a.set_xticks(x, [f"{l}\n(n={k:,})" for l, k in zip(lab, stats.n)], fontsize=7)
    a.set_xlabel(xlabel)
    a.set_ylabel("fraction of perturbations")
    a.set_ylim(0, 1)
    a.legend(frameon=False, fontsize=7)


a = ax[0, 0]
a.hist(info.n_sig, bins=60, color=MARK)
a.set_yscale("log")
a.set_xlabel(f"significant features (q <= {args.qvalue}) of {n_feat:,}")
a.set_ylabel("perturbations (log)")
a.set_title(
    f"{(info.n_sig == 0).mean():.0%} with none, {(info.frac_sig > 0.5).mean():.0%} with > half",
    fontsize=9, color=INK, loc="left"
)
scatter(ax[0, 1], "n_ctrl", "control samples (log)")
scatter(ax[0, 2], "n_treated", "treated samples (log)")
frac_bars(ax[1, 0], ctrl_bins, "control samples")
frac_bars(ax[1, 1], info.time_h, "time (h)")
dose_bins = pd.qcut(info.dose, 6, duplicates="drop") if info.dose.notna().any() else info.n_batches
frac_bars(ax[1, 2], dose_bins, "dose (sextiles)" if info.dose.notna().any() else "batches")
fig.suptitle(os.path.basename(os.path.normpath(d)), x=0.01, ha="left", color=INK)
fig.tight_layout()
png = os.path.join(out_dir, "n_sig_diagnostics.png")
fig.savefig(png, dpi=150)
print(f"saved {png}")