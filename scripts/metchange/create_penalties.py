################################################################################
####################----------------------------------------####################
##########----------            MetaboReprogram             ----------##########
# Convert gene symbols to ensembl, aggregate duplicated genes and create       #
# penalty matrices, allowing different approaches for obtaing them.            #
#                                                                              #
################################################################################

################################################################################
# MODULE LOADING                                                               #
################################################################################

import os
from scipy.stats import gaussian_kde, norm, rankdata, skewnorm, gamma as gamma_dist
from scipy.optimize import minimize
from scipy.signal import argrelextrema
import argparse
import numpy as np
import pandas as pd
import polars as pl
import cobra
from cobra.io import read_sbml_model, write_sbml_model
from collections import defaultdict

def convert_ids_ensembl(df, hgnc, symbol_col="V1", is_log2=True):
    symbol_to_ensembl = (
        hgnc.groupby("symbol")["ensembl_gene_id"]
        .apply(lambda x: "|".join(x.dropna().astype(str)))
        .to_dict()
    )
    alias_to_ensembl = defaultdict(list)
    for _, row in hgnc.dropna(subset=["alias_symbol"]).iterrows():
        for alias in row["alias_symbol"].split("|"):
            alias_to_ensembl[alias].append(str(row["ensembl_gene_id"]))
    alias_to_ensembl = {
        k: "|".join(sorted(set(v)))
        for k, v in alias_to_ensembl.items()
    }
    ensembl_list = []
    for gene in df[symbol_col].to_list():
        if gene in symbol_to_ensembl:
            ensembl_list.append(symbol_to_ensembl[gene])
        else:
            ensembl_list.append(alias_to_ensembl.get(gene, ""))
    symb_to_ensembl = pd.DataFrame(
        {
            'hgnc_symbol': df[symbol_col].to_list(),
            'ensembl_gene_id': ensembl_list
        }
    )
    # Add ensembl column
    df = df.with_columns(
        pl.Series("ensembl_gene_id", ensembl_list)
    )
    # Remove unmatched rows
    df = df.filter(
        pl.col("ensembl_gene_id") != ""
    )
    if is_log2:
        df = df.with_columns(
            pl.exclude(["ensembl_gene_id", symbol_col]).pow(2)
        )
    # Sum eventual duplicates
    df = (
        df.group_by("ensembl_gene_id")
        .agg(pl.exclude([symbol_col, "ensembl_gene_id"]).sum())
        .sort("ensembl_gene_id")
    )
    if is_log2:
        df = df.with_columns(
            pl.exclude("ensembl_gene_id").sqrt()
        )
    return (df)

# Functions for fitting mixture model
def fit_zfpkm(values, zero_cut=0.5, max_n=200_000, grid_n=2000, seed=0):
    """Fit the zFPKM half-Gaussian on the expressed mode of one batch."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    frac_zero = np.mean(v <= zero_cut)
    v = v[v > zero_cut]                      # drop zero spike ("not detected")
    rng = np.random.default_rng(seed)
    v_kde = rng.choice(v, size=min(max_n, v.size), replace=False)
    kde = gaussian_kde(v_kde)
    grid = np.linspace(v.min(), v.max(), grid_n)
    dens = kde(grid)
    mu = grid[np.argmax(dens)]
    U = v[v > mu].mean()
    sigma = (U - mu) * np.sqrt(np.pi / 2)
    # Flag fits where the peak sits at the left edge (decay-like batches)
    edge = (mu - v.min()) / (v.max() - v.min())
    return {
        "mu": mu, "sigma": sigma, "zero_cut": zero_cut,
        "frac_zero": frac_zero, "peak_rel_pos": edge,
        "suspect": edge < 0.05,
    }

def valley_cut(values, bw=0.1, max_n=200_000, seed=0, grid_n=1000):
    v = np.asarray(values, float)
    v = v[np.isfinite(v) & (v > 0)]                  # drop the exact-zero spike
    rng = np.random.default_rng(seed)
    v = rng.choice(v, size=min(max_n, v.size), replace=False)
    lo, hi = np.percentile(v, [2, 90])               # keep the search out of the tail
    grid = np.linspace(lo, hi, grid_n)
    d = gaussian_kde(v, bw_method=bw)(grid)
    mins = argrelextrema(d, np.less)[0]
    if len(mins) == 0:
        return np.nan, False
    m = mins[np.argmin(d[mins])]                     # deepest valley in the window
    ok = d[m] < 0.8 * min(d[:m].max(), d[m:].max())  # must be a real dip on both sides
    return grid[m], ok

def crossing(p_fn, lo=0.05, hi=8, n=4000):
    xs = np.linspace(lo, hi, n)
    return xs[np.argmax(p_fn(xs) >= 0.5)]

def _unpack(t, low, high):
    w = 1 / (1 + np.exp(-t[0]))
    if low == "expon":
        lam = np.exp(t[1]); lo = lambda v: lam * np.exp(-lam * v); k = 2
    else:                                          # gamma(shape a, rate b)
        a, b = np.exp(t[1]), np.exp(t[2])
        lo = lambda v: gamma_dist.pdf(v, a, scale=1 / b); k = 3
    mu, sig = t[k], np.exp(t[k + 1])
    if high == "norm":
        hi = lambda v: norm.pdf(v, mu, sig)
    else:
        al = t[k + 2]; hi = lambda v: skewnorm.pdf(v, al, loc=mu, scale=sig)
    return w, lo, hi

def fit_mixture(x, low, high, init):
    t0 = [np.log(init["w"] / (1 - init["w"]))]
    t0 += [np.log(init["lam"])] if low == "expon" else [0.0, np.log(init["lam"])]
    t0 += [init["mu"], np.log(init["sigma"])]
    if high == "skewnorm":
        t0 += [-1.0]
    def nll(t):
        w, lo, hi = _unpack(t, low, high)
        return -np.log((1 - w) * lo(x) + w * hi(x) + 1e-300).mean()
    res = minimize(
        nll,
        t0,
        method="Nelder-Mead",
        options={"maxiter": 4000, "xatol": 1e-4, "fatol": 1e-8}
    )
    if not res.success:
        print("warning: Nelder-Mead stopped early:", res.message)
    return res.x

def heldout_ll(t, x, low, high):
    w, lo, hi = _unpack(t, low, high)
    return np.log((1 - w) * lo(x) + w * hi(x) + 1e-300).mean()

def posterior_from_fit(t, low, high, x_max, grid_n=4000):
    w, lo, hi = _unpack(t, low, high)
    grid = np.linspace(1e-3, x_max, grid_n)
    h, l = w * hi(grid), (1 - w) * lo(grid)
    post = np.maximum.accumulate(h / (h + l + 1e-300))
    def p_expr(v):
        v = np.asarray(v, float)
        return np.where(v <= 0, 0.0, np.interp(v, grid, post, left=0.0, right=1.0))
    return p_expr

# Function that creates four penalties: the saturating one, and three versions of graded penalties.

def make_penalties(t, p_main):
    # layout for the exponential + skew-normal fit: t = [logit w, log lam, mu, log sigma, alpha]
    mu_h, sig_h, al_h = t[2], np.exp(t[3]), t[4]
    sd_h = skewnorm.std(al_h, loc=mu_h, scale=sig_h)
    def graded(x, offset=0.0):
        x = np.asarray(x, float)
        g = skewnorm.sf(x + offset * sd_h, al_h, loc=mu_h, scale=sig_h)
        return np.where(x <= 0, 1.0, g)
    def combined(x, offset=0.0):
        x = np.asarray(x, float)
        P = p_main(x)
        return (1 - P) + P * graded(x, offset)
    return {
        "saturating": lambda x: 1 - p_main(np.asarray(x, float)),   # MetChange-style presence penalty
        "graded_0":   lambda x: graded(x, 0.0),
        "combined_0": lambda x: combined(x, 0.0),
        "combined_1": lambda x: combined(x, 1.0),
    }

def create_penalty_df(pen_arr, expr, row_name_col="ensembl_gene_id"):
    df = pl.DataFrame(
        pen_arr,
        schema=expr.columns[1:],
        orient="row"
    ).insert_column(0, expr[row_name_col])
    return (df)

################################################################################
# Terminal argument parser                                                     #
################################################################################
parser = argparse.ArgumentParser(
    description='Obtain penalty matrices.'
)

parser.add_argument(
        '--expr',
        type=str,
        help='Expression file.'
)

parser.add_argument(
        '--hgnc',
        type=str,
        help='HGNC file, for translating symbols to ENSEMBL.'
)

parser.add_argument(
        '--convert_to_ensembl',
        action="store_true",
        help='Whether the gene names need to be converted to ENSEMBL.'
)

parser.add_argument(
        '--pen_mode',
        type=str,
        help='Penalty mode.'
)

parser.add_argument(
        '--threshold',
        type=int,
        help='Threshold percentile to considere a gene not expressed at all (for ranks penalty).',
        default=2
)

parser.add_argument(
        '--gem',
        type=str,
        help='GEM file. To filter genes to keep only the ones in the GEM in use, and save storage space'
)

parser.add_argument('--out_dir', 
	                  type=str,
                    help='Directory where the result will be stored.', 
                    default=os.getcwd())

args = parser.parse_args()

################################################################################
# DIRECTORY STUFF                                                              #
################################################################################

expr_f="/home/users/gsantamaria/projects/metabo_pert/data/LINCS/lincs_A549_expr.csv"
hgnc_f="/home/users/gsantamaria/projects/metabo_pert/data/utilities/hgnc_complete_set.txt"
pen_mode="mixture"
threshold=2
gem_f="/home/users/gsantamaria/projects/metabo_pert/data/GEMs/Human2/Human-GEM/model/Human-GEM.xml"
out_dir="/scratch/users/gsantamaria/test_large_files/metabo_perts/results/gene_penalties/mixture_script/"

expr_f = args.expr
hgnc_f = args.hgnc
convert_to_ensembl = args.convert_to_ensembl
pen_mode = args.pen_mode
threshold = args.threshold
gem_f = args.gem
out_dir = args.out_dir

os.makedirs(out_dir, exist_ok=True)

################################################################################
# DATA LOADING                                                                 #
################################################################################
print("Loading expression data from %s..." % expr_f)
expr = pl.read_csv(expr_f)
hgnc = pd.read_csv(hgnc_f, sep = '\t')
print("Loading GEM from %s..." % gem_f)
gem = read_sbml_model(gem_f)

print("Converting gene symbols to ENSEMBL...")
if convert_to_ensembl:
    expr = convert_ids_ensembl(expr, hgnc, symbol_col="V1", is_log2=True)

################################################################################
# GET PENALTIES                                                                #
################################################################################
if pen_mode == "rank":
    print("Computing rank-based penalties...")
    gene_ids = expr["ensembl_gene_id"]
    sample_cols = [c for c in expr.columns if c != "ensembl_gene_id"]
    n_genes = expr.height
    chunk_size = 128
    pen_cols = []
    for start in range(0, len(sample_cols), chunk_size):
        cols = sample_cols[start:start + chunk_size]
        # Shape: genes x samples
        x = expr.select(cols).to_numpy()
        # Rank each sample independently (ascending, average ties)
        ranks = rankdata(x, axis=0, method="average")
        pvals = 1.0 - ranks / n_genes
        # Calculate the percentile separately for each sample
        cutoffs = np.percentile(pvals, 100 - threshold, axis=0)
        # Apply penalty
        pvals[pvals >= cutoffs] = 1.0
        # Save the processed columns for reconstruction
        pen_cols.extend(
            pl.Series(name=col, values=pvals[:, j])
            for j, col in enumerate(cols)
        )
    # Reconstruct the Polars DataFrame, preserving column names
    pen_df = pl.DataFrame([gene_ids, *pen_cols])
elif pen_mode = "mixture":
    pooled = expr.select(pl.exclude("ensembl_gene_id")).to_numpy().flatten()
    pos = pooled[pooled > 0]
    pooled_cut = valley_cut(pooled)
    params_pooled = fit_zfpkm(pooled, zero_cut=pooled_cut[0])
    #rng = np.random.default_rng(1)
    #idx = rng.permutation(pos.size)
    #x_tr, x_te = pos[idx[:100_000]], pos[idx[100_000:300_000]]
    #v = rng.choice(pos, size=500_000, replace=False)
    init = dict(w=np.mean(pos > pooled_cut[0]), lam=1.0,
                mu=params_pooled["mu"], sigma=params_pooled["sigma"])
    t = fit_mixture(pos[idx[:500_000]], "expon", "skewnorm", init)
    print("skew-normal shape:", round(t[4], 2), "| w:", round(1/(1+np.exp(-t[0])), 3))
    p_main = posterior_from_fit(t, "expon", "skewnorm", pos.max())
    pen_fns = make_penalties(t, p_main)
    pen_arr = pen_fns["saturating"](expr[:, 1:expr.shape[1]])
    pen_df = create_penalty_df(pen_arr, expr)
    

################################################################################
# FILTER AND SAVE                                                              #
################################################################################

# Keep only genes in the GEM
gem_genes = [x.id for x in gem.genes]
pen_df = pen_df.filter(pl.col("ensembl_gene_id").is_in(gem_genes))

out_name = "%s_%s.csv" % (os.path.basename(expr_f).split('.')[0], pen_mode)
pen_df.write_csv(os.path.join(out_dir, out_name))
print("%s saved in %s." % (out_name, out_dir))