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
import sys
from scipy.stats import gaussian_kde, norm, rankdata, skewnorm, gamma as gamma_dist
from scipy.optimize import minimize
from scipy.signal import argrelextrema
import argparse
import numpy as np
import pandas as pd
from cobra.io import read_sbml_model

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "parse_lincs_data"))
from gene_id_mapping import symbol_to_ensembl_map, read_matrix, gem_symbol_overrides  # noqa: E402

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
# The matrix is read in row chunks into one float32 array (all genes are kept:
# the rank penalty ranks each sample over every gene). Only the GEM genes are
# written out.
print("Loading GEM from %s..." % gem_f)
gem = read_sbml_model(gem_f)
gem_genes = set(g.id for g in gem.genes)

mapping = None
if convert_to_ensembl:
    print("Gene symbols will be converted to ENSEMBL while reading.")
    mapping = symbol_to_ensembl_map(pd.read_csv(hgnc_f, sep='\t', dtype=str), gem_symbol_overrides(gem))
print("Loading expression data from %s..." % expr_f)
gene_ids, samples, X, _ = read_matrix(expr_f, mapping=mapping, is_log2=True)
print("%d genes x %d samples" % X.shape)
gem_rows = np.array([i for i, g in enumerate(gene_ids) if g in gem_genes])
print("%d of %d GEM genes present" % (len(gem_rows), len(gem_genes)))

################################################################################
# GET PENALTIES                                                                #
################################################################################
if pen_mode == "rank":
    print("Computing rank-based penalties...")
    n_genes = X.shape[0]
    chunk_size = 128
    pen = np.empty((len(gem_rows), X.shape[1]), np.float32)
    for start in range(0, X.shape[1], chunk_size):
        # Shape: genes x samples. Rank each sample independently (ascending, average ties)
        x = X[:, start:start + chunk_size]
        ranks = rankdata(x, axis=0, method="average")
        pvals = 1.0 - ranks / n_genes
        # Calculate the percentile separately for each sample
        cutoffs = np.percentile(pvals, 100 - threshold, axis=0)
        # Apply penalty
        pvals[pvals >= cutoffs] = 1.0
        pen[:, start:start + chunk_size] = pvals[gem_rows]
elif pen_mode == "mixture":
    print("Fitting the expression mixture model...")
    rng = np.random.default_rng(1)
    # A random sample of all values is enough for the fit (the fitting functions
    # subsample anyway); flattening the whole matrix would copy it.
    flat = X.ravel()
    pooled = flat[rng.integers(0, flat.size, size=min(flat.size, 5_000_000))].astype(np.float64)
    pooled = pooled[np.isfinite(pooled)]
    pos = pooled[pooled > 0]
    pooled_cut, cut_ok = valley_cut(pooled)
    if not np.isfinite(pooled_cut) or not cut_ok:
        print("warning: no clear valley in the expression density; using zero_cut = 0.5")
        pooled_cut = 0.5
    params_pooled = fit_zfpkm(pooled, zero_cut=pooled_cut)
    init = dict(w=np.mean(pos > pooled_cut), lam=1.0,
                mu=params_pooled["mu"], sigma=params_pooled["sigma"])
    t = fit_mixture(rng.choice(pos, size=min(500_000, pos.size), replace=False), "expon", "skewnorm", init)
    print("skew-normal shape:", round(t[4], 2), "| w:", round(1/(1+np.exp(-t[0])), 3))
    p_main = posterior_from_fit(t, "expon", "skewnorm", float(np.nanmax(X)))
    pen_fns = make_penalties(t, p_main)
    pen = np.empty((len(gem_rows), X.shape[1]), np.float32)
    for start in range(0, X.shape[1], 1024):
        pen[:, start:start + 1024] = pen_fns["saturating"](X[gem_rows, start:start + 1024])
else:
    sys.exit("Unsupported pen_mode: %s (rank or mixture)" % pen_mode)

################################################################################
# SAVE                                                                         #
################################################################################
out_name = "%s_%s.csv" % (os.path.basename(expr_f).split('.')[0], pen_mode)
pd.DataFrame(pen, index=[gene_ids[i] for i in gem_rows], columns=samples).to_csv(
    os.path.join(out_dir, out_name), index_label="ensembl_gene_id", float_format="%.7g")
print("%s saved in %s." % (out_name, out_dir))
