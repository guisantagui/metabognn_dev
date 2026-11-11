#!/usr/bin/env python
"""
Plots for check_signature_target_similarity.py, from its output folder (no recomputation).

  1_auroc_summary         AUROC of each pair comparison against its permuted-label null, and the kNN enrichment
  2_corr_distributions    signature correlation of same-drug, shared-target and unrelated pairs (pairs_sample.tsv.gz)
  3_target_overlap        mean correlation and AUROC as the target sets of two drugs overlap more (Jaccard bins)
  4_per_target            targets whose drugs are most / least recognisable from their signatures (per_target.tsv)
  5_example_targets       correlation heatmap of the drugs of a few well-recognised targets (example_targets_corr.tsv.gz)

Files missing from older runs (pairs_sample, example_targets_corr) are skipped. Each figure is written as PNG and PDF.

Usage: python plot_signature_target_similarity.py --dir <check output dir> [--out-dir <dir>/plots] [--title A549]
"""
import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dir", required=True, help="check_signature_target_similarity.py --out-dir")
ap.add_argument("--out-dir", help="default: <dir>/plots")
ap.add_argument("--title", default="", help="prefix of the figure titles, e.g. the cell line")
ap.add_argument("--top", type=int, default=20, help="4_per_target: best targets shown")
ap.add_argument("--bottom", type=int, default=5, help="4_per_target: worst targets shown")
args = ap.parse_args()
out = args.out_dir or os.path.join(args.dir, "plots")
os.makedirs(out, exist_ok=True)

# colour follows the pair class in every figure (Okabe-Ito, colour-blind safe); null / reference in grey
COL = {"same drug": "#0072B2", "shared target": "#E69F00", "no shared target": "#8C8C8C"}
INK, MUTED, GRID = "#222222", "#6B6B6B", "#E3E3E3"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 10, "axes.titleweight": "bold", "legend.frameon": False})
title = (args.title + ": ") if args.title else ""


def path(name):
    return os.path.join(args.dir, name)


def save(fig, name):
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}/{name}.png/.pdf")


def hgrid(ax):
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)


with open(path("summary.json")) as fh:
    S = json.load(fh)

# ---- 1: AUROC per comparison vs null, kNN enrichment
fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.0, 3.3), gridspec_kw=dict(width_ratios=[3, 1.2]))
comps = [("same_drug", "same drug\n(other dose/time)", COL["same drug"]),
         ("shared_target", "shared\ntarget", COL["shared target"]),
         ("shared_target_same_time", "shared target\nsame time", COL["shared target"]),
         ("shared_target_drug", "shared target\ndrug means", COL["shared target"])]
null = S.get("null_mean", {})
for x, (k, lab, col) in enumerate(comps):
    v = S.get(k)
    if v is None:
        continue
    a1.plot([x, x], [0.5, v], color=col, lw=2, solid_capstyle="round")
    a1.scatter(x, v, s=48, color=col, zorder=3, edgecolor="white", lw=1)
    a1.annotate(f"{v:.2f}", (x, v), textcoords="offset points", xytext=(0, 6), ha="center", color=INK)
    if k in null:
        a1.scatter(x, null[k], s=40, facecolor="white", edgecolor=MUTED, lw=1.2, zorder=3,
                   label="permuted drug labels" if x == 0 else None)
a1.axhline(0.5, color=MUTED, lw=0.8, ls="--")
a1.set_xticks(range(len(comps)), [c[1] for c in comps])
a1.set_ylabel("AUROC (signature correlation)")
a1.set_ylim(0.4, 1.0)
a1.set_xlim(-0.6, len(comps) - 0.4)
a1.set_title("Pairs vs pairs of drugs sharing no target", loc="left")
a1.legend(loc="upper right")
hgrid(a1)
knn = S.get("knn", {})
vals = [("observed", knn.get("frac_sharing_target"), COL["shared target"]), ("random", knn.get("expected"), MUTED)]
for x, (lab, v, col) in enumerate(vals):
    if v is not None:
        a2.bar(x, v, width=0.6, color=col)
        a2.annotate(f"{v:.1%}", (x, v), textcoords="offset points", xytext=(0, 3), ha="center", color=INK)
a2.set_xticks(range(2), [v[0] for v in vals])
a2.set_ylabel(f"{knn.get('k', '')} nearest neighbours sharing a target")
a2.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=0))
a2.set_title(f"kNN: {knn.get('enrichment', float('nan')):.1f}x", loc="left")
hgrid(a2)
fig.suptitle(f"{title}{S['n_perturbations']:,} active perturbations, {S['n_drugs']:,} drugs", x=0.01, ha="left",
             fontsize=10, color=MUTED)
save(fig, "1_auroc_summary")

# ---- 2: correlation distributions by pair class
if os.path.exists(path("pairs_sample.tsv.gz")):
    ps = pd.read_csv(path("pairs_sample.tsv.gz"), sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    bins = np.linspace(-0.4, 1.0, 71)
    for ax, (sub, lab) in zip(axes, ((ps, "all pairs"), (ps[ps.same_time.astype(str) == "True"], "same time point"))):
        for cls, col in COL.items():
            v = sub.loc[sub.pair_class == cls, "corr"]
            if len(v):
                ax.hist(v, bins=bins, density=True, histtype="step", lw=1.6, color=col, label=f"{cls} (n={len(v):,})")
                ax.axvline(v.mean(), color=col, lw=0.9, ls="--")
        ax.set_title(lab, loc="left")
        ax.set_xlabel("Spearman correlation of signatures")
        hgrid(ax)
    axes[0].set_ylabel("density")
    axes[1].legend(loc="upper right", fontsize=8)
    fig.suptitle(f"{title}signature similarity by pair class (dashed: mean; sample per class)", x=0.01, ha="left",
                 fontsize=10, color=MUTED)
    save(fig, "2_corr_distributions")
else:
    print("pairs_sample.tsv.gz not found (older run): 2_corr_distributions skipped")

# ---- 3: graded by target overlap
ov = S.get("by_target_overlap", {})
if ov:
    names = list(ov)
    short = [n.replace("jaccard ", "").replace("identical target sets", "identical") for n in names]
    x = np.arange(len(names))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    base = S["mean_corr"]["no_shared_target"]
    for key, lab, mk in (("mean_corr_same_time", "same time point", "o"), ("mean_corr_diff_time", "different time", "s")):
        y = [ov[n][key] if ov[n][key] is not None else np.nan for n in names]
        a1.plot(x, y, color=COL["shared target"], marker=mk, lw=1.6, ms=6, mec="white",
                alpha=1 if mk == "o" else 0.55, label=lab)
    a1.axhline(base, color=COL["no shared target"], lw=1, ls="--", label="no shared target")
    a1.set_ylabel("mean signature correlation")
    a1.set_title("Similarity grows with target overlap", loc="left")
    a1.legend(fontsize=8)
    a2.plot(x, [ov[n]["auroc_vs_no_shared"] for n in names], color=COL["shared target"], marker="o", lw=1.6, ms=6,
            mec="white")
    for xi, n in zip(x, names):
        a2.annotate(f"{ov[n]['n_drug_pairs']:,}", (xi, ov[n]["auroc_vs_no_shared"]), textcoords="offset points",
                    xytext=(0, 7), ha="center", fontsize=7, color=MUTED)
    a2.axhline(0.5, color=MUTED, lw=0.8, ls="--")
    a2.set_ylabel("AUROC vs no shared target")
    a2.set_title("labels: drug pairs per bin", loc="left", fontweight="normal", color=MUTED, fontsize=8)
    for ax in (a1, a2):
        ax.set_xticks(x, short, rotation=25, ha="right")
        ax.set_xlabel("Jaccard of the two drugs' target sets")
        hgrid(ax)
    fig.suptitle(f"{title}pairs of different drugs", x=0.01, ha="left", fontsize=10, color=MUTED)
    save(fig, "3_target_overlap")

# ---- 4: per target
if os.path.exists(path("per_target.tsv")):
    pt = pd.read_csv(path("per_target.tsv"), sep="\t").sort_values("auroc", ascending=False)
    if len(pt):
        show = pd.concat([pt.head(args.top), pt.tail(args.bottom)]).drop_duplicates("target")
        gap = len(pt) > args.top + args.bottom
        y = np.arange(len(show))[::-1].astype(float)
        if gap:
            y[min(args.top, len(show)):] -= 1  # break between best and worst
        fig, ax = plt.subplots(figsize=(5.2, 0.2 * len(show) + 1.2))
        size = 12 + 60 * np.sqrt(show.n_drugs / show.n_drugs.max())
        ax.hlines(y, 0.5, show.auroc, color=COL["shared target"], lw=1.2)
        ax.scatter(show.auroc, y, s=size, color=COL["shared target"], edgecolor="white", lw=0.8, zorder=3)
        ax.axvline(0.5, color=MUTED, lw=0.8, ls="--")
        ax.set_yticks(y, [f"{t} ({n})" for t, n in zip(show.target, show.n_drugs)], fontsize=7.5)
        ax.set_xlabel("AUROC: both drugs hit the target vs only one does")
        ax.set_title(f"{title}per target (n drugs; {len(pt):,} targets with enough drugs)", loc="left", fontsize=9)
        ax.grid(axis="x", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        save(fig, "4_per_target")

# ---- 5: example heatmap
if os.path.exists(path("example_targets_corr.tsv.gz")):
    ex = pd.read_csv(path("example_targets_corr.tsv.gz"), sep="\t")
    meta, M = ex[["target", "drug", "timepoint", "perturbation"]], ex.iloc[:, 4:].to_numpy(float)
    np.fill_diagonal(M, np.nan)
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    lim = np.nanpercentile(np.abs(M), 99)
    im = ax.imshow(M, cmap="RdBu_r", vmin=-lim, vmax=lim, interpolation="nearest")
    ax.set_facecolor("#F2F2F2")
    edges = np.flatnonzero(meta.target.to_numpy()[1:] != meta.target.to_numpy()[:-1]) + 1
    for e in edges:
        ax.axhline(e - 0.5, color="white", lw=1.5)
        ax.axvline(e - 0.5, color="white", lw=1.5)
    starts = np.r_[0, edges]
    ends = np.r_[edges, len(meta)]
    mids = (starts + ends - 1) / 2
    ax.set_xticks(mids, meta.target[starts], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(mids, [f"{meta.target[s]} ({meta.drug[s:e].nunique()} drugs)" for s, e in zip(starts, ends)],
                  fontsize=8)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("Spearman correlation")
    cb.outline.set_visible(False)
    ax.set_title(f"{title}perturbations of drugs hitting well-recognised targets", loc="left", fontsize=9)
    save(fig, "5_example_targets")
else:
    print("example_targets_corr.tsv.gz not found (older run): 5_example_targets skipped")
