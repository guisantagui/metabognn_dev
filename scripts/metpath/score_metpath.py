################################################################################
####################----------------------------------------####################
##########----------            MetaboReprogram             ----------##########
# MetPath pathways (Mattei et al. 2023), step 2: score differential expression #
# on the pathways of build_metpath_pathways.py (mapGenes / calcRes /           #
# pathStats / calcAggregateScores).                                            #
#                                                                              #
# - Reaction fold change: mean fold change of all the genes in its GPR that    #
#   were measured (MetPath ignores the AND / OR logic here).                   #
# - Pathway score: sum over its reactions with data of |weight| x fold change, #
#   weights normalised to sum to 1 over those reactions.                       #
# - p-values: the score against --n_perms pathways of the same weights on     #
#   reactions drawn at random from all reactions with data (two-sided). With   #
#   --n_perms 0 a normal approximation with the exact permutation mean and     #
#   variance is used instead.                                                  #
# - Aggregate perturbation score (APS) = production + consumption score.       #
#                                                                              #
# --fc_scale ratio reproduces MetPath (fold change = 2^log2FC, 1 = unchanged); #
# --fc_scale log2 averages log2 fold changes instead (0 = unchanged).          #
#                                                                              #
# Input: genes x perturbations log2 fold changes (e.g. limma_DE.R log2FC.csv   #
# of the expression; symbols or ENSEMBL ids). Output, metabolites x            #
# perturbations: production_score.csv, consumption_score.csv, aps.csv, and     #
# production_/consumption_ pvals.csv (two-sided), qvals.csv (BH within each    #
# perturbation and direction), z.csv.                                          #
################################################################################

################################################################################
# MODULE LOADING                                                               #
################################################################################

import os
import sys
import argparse

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import norm
from cobra.io import read_sbml_model

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "parse_lincs_data"))
from gene_id_mapping import symbol_to_ensembl_map, gem_symbol_overrides  # noqa: E402

################################################################################
# FUNCTIONS                                                                    #
################################################################################


def read_de(de_f, gem, hgnc_f):
    """Genes x perturbations log2 fold changes, rows mapped to GEM gene ids (duplicates averaged)."""
    de = pd.read_csv(de_f, index_col=0)
    ids = de.index.astype(str)
    if not ids.str.startswith("ENSG").all():
        if hgnc_f is None:
            sys.exit("--de has gene symbols: --hgnc is needed to map them to the GEM's ENSEMBL ids")
        mapping = symbol_to_ensembl_map(pd.read_csv(hgnc_f, sep="\t", dtype=str), gem_symbol_overrides(gem))
        de.index = ids.map(lambda s: mapping.get(s))
    gem_genes = {g.id for g in gem.genes}
    de = de[de.index.isin(gem_genes)]
    return de.groupby(level=0).mean()


def reaction_fold_changes(gem, fc):
    """Mean over the measured genes of each reaction's GPR (mapGenes). Returns reactions x perturbations with
    NaN where no gene was measured."""
    gidx = {g: i for i, g in enumerate(fc.index)}
    rows, cols = [], []
    for j, r in enumerate(gem.reactions):
        for g in r.genes:
            if g.id in gidx:
                rows.append(j)
                cols.append(gidx[g.id])
    G = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(gem.reactions), len(fc)))
    X = fc.to_numpy(dtype=np.float64)
    ok = ~np.isnan(X)
    num, den = G @ np.where(ok, X, 0.0), G @ ok.astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        R = num / den
    R[den == 0] = np.nan
    return R


def score_pathways(paths, rxn_index, R, n_perms, seed):
    """Pathway scores, z-scores and two-sided p-values for every pathway (rows of paths) x perturbation."""
    keys = list(paths.groupby(["metabolite", "direction"]).groups)
    key_idx = {k: i for i, k in enumerate(keys)}
    W = sp.csr_matrix((np.abs(paths.weight.to_numpy()),
                       ([key_idx[k] for k in zip(paths.metabolite, paths.direction)],
                        paths.reaction.map(rxn_index).to_numpy())), shape=(len(keys), R.shape[0]))
    ok = ~np.isnan(R)
    R0 = np.where(ok, R, 0.0)
    wsum = W @ ok.astype(float)                        # total weight of the reactions with data
    with np.errstate(invalid="ignore", divide="ignore"):
        score = (W @ R0) / wsum
        w2 = (W.multiply(W)) @ ok.astype(float) / wsum ** 2   # sum of squared normalised weights
    score[wsum == 0] = np.nan
    # null: the same normalised weights on reactions drawn without replacement from all reactions with data
    N = ok.sum(0)
    mu = np.nansum(R0, 0) / N
    var = np.nanvar(np.where(ok, R, np.nan), 0)        # population variance of the pool
    null_var = var * N / (N - 1) * (w2 - 1.0 / N)
    if n_perms:
        p = permutation_pvalues(W, R, ok, score, n_perms, seed)
        z = norm.isf(p / 2) * np.sign(score - mu)
    else:
        with np.errstate(invalid="ignore", divide="ignore"):
            z = (score - mu) / np.sqrt(np.clip(null_var, 0, None))
        p = 2 * norm.sf(np.abs(z))
    return keys, score, z, p


def permutation_pvalues(W, R, ok, score, n_perms, seed, block=200000):
    """pathStats: the score against n_perms pathways with the same normalised weights on reactions drawn from all
    reactions with data, two-sided (2 x min(pValLow, pValHigh)). Draws are with replacement (pathways have a few
    reactions out of thousands, so this is the same up to rounding), vectorised over pathways of equal size."""
    rng = np.random.default_rng(seed)
    p = np.full(score.shape, np.nan)
    W = W.tocsr()
    for c in range(R.shape[1]):
        pool = R[ok[:, c], c]
        groups = {}
        for i in range(W.shape[0]):
            cols = W.indices[W.indptr[i]:W.indptr[i + 1]]
            w = W.data[W.indptr[i]:W.indptr[i + 1]][ok[cols, c]]
            if len(w) and not np.isnan(score[i, c]):
                groups.setdefault(len(w), []).append((i, w / w.sum()))
        for k, members in groups.items():
            idx = np.array([i for i, _ in members])
            w = np.stack([w for _, w in members])                          # (n_paths, k)
            per = max(1, block // (n_perms * k))
            for s in range(0, len(idx), per):
                draws = pool[rng.integers(0, len(pool), size=(len(idx[s:s + per]), n_perms, k))]
                null = np.einsum("pnk,pk->pn", draws, w[s:s + per])
                obs = score[idx[s:s + per], c][:, None]
                low, high = (obs >= null).mean(1), (obs <= null).mean(1)
                p[idx[s:s + per], c] = np.minimum(1.0, 2 * np.minimum(low, high))
    return p


def bh(p):
    """Benjamini-Hochberg along axis 0, NaN ignored."""
    q = np.full(p.shape, np.nan)
    for c in range(p.shape[1]):
        ok = ~np.isnan(p[:, c])
        x = p[ok, c]
        order = np.argsort(x)
        ranked = x[order] * len(x) / np.arange(1, len(x) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        out = np.empty(len(x))
        out[order] = np.minimum(ranked, 1.0)
        q[ok, c] = out
    return q


################################################################################
# MAIN                                                                         #
################################################################################

def main():
    parser = argparse.ArgumentParser(description="Score differential expression on MetPath pathways.")
    parser.add_argument("--pathways", type=str, required=True, help="build_metpath_pathways.py metpath_pathways.tsv")
    parser.add_argument("--gem", type=str, required=True, help="GEM file (SBML), for the GPRs.")
    parser.add_argument("--de", type=str, required=True, help="Genes x perturbations log2 fold changes (CSV).")
    parser.add_argument("--hgnc", type=str, default=None, help="HGNC table, if --de rows are gene symbols.")
    parser.add_argument("--fc_scale", choices=["ratio", "log2"], default="ratio",
                        help="ratio: fold changes 2^log2FC as in MetPath (1 = unchanged); log2: log2FC (0 = unchanged).")
    parser.add_argument("--n_perms", type=int, default=1000,
                        help="Permutations per pathway and perturbation (MetPath: 1000); 0: exact mean and variance "
                             "of the permutation distribution with a normal approximation (fast, poor for pathways "
                             "of a few reactions).")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out_dir", type=str, default=os.getcwd(), help="Directory where the result will be stored.")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    print("Loading GEM: %s" % args.gem)
    gem = read_sbml_model(args.gem)
    paths = pd.read_csv(args.pathways, sep="\t")
    de = read_de(args.de, gem, args.hgnc)
    print("%d GEM genes x %d perturbations in --de" % de.shape)
    fc = np.exp2(de) if args.fc_scale == "ratio" else de
    R = reaction_fold_changes(gem, fc)
    print("%d reactions with data" % (~np.isnan(R)).any(1).sum())
    rxn_index = {r.id: j for j, r in enumerate(gem.reactions)}
    keys, score, z, p = score_pathways(paths, rxn_index, R, args.n_perms, args.seed)
    q = np.vstack([bh(p[[i for i, k in enumerate(keys) if k[1] == d]]) for d in ("production", "consumption")])
    order = [i for d in ("production", "consumption") for i, k in enumerate(keys) if k[1] == d]
    keys = [keys[i] for i in order]
    score, z, p = score[order], z[order], p[order]

    cols = de.columns
    for d in ("production", "consumption"):
        idx = [i for i, k in enumerate(keys) if k[1] == d]
        mets = [keys[i][0] for i in idx]
        for name, mat in (("score", score), ("z", z), ("pvals", p), ("qvals", q)):
            pd.DataFrame(mat[idx], index=mets, columns=cols).to_csv(
                os.path.join(args.out_dir, "%s_%s.csv" % (d, name)), float_format="%.6g")
    prod = pd.DataFrame(score[[i for i, k in enumerate(keys) if k[1] == "production"]],
                        index=[k[0] for k in keys if k[1] == "production"], columns=cols)
    cons = pd.DataFrame(score[[i for i, k in enumerate(keys) if k[1] == "consumption"]],
                        index=[k[0] for k in keys if k[1] == "consumption"], columns=cols)
    aps = prod.reindex(prod.index.union(cons.index)).fillna(0) + cons.reindex(prod.index.union(cons.index)).fillna(0)
    aps.to_csv(os.path.join(args.out_dir, "aps.csv"), float_format="%.6g")
    for d in ("production", "consumption"):
        qq = pd.read_csv(os.path.join(args.out_dir, "%s_qvals.csv" % d), index_col=0)
        print("%s: %d pathways, median %.0f with q < 0.05 per perturbation" % (
            d, len(qq), float(np.median((qq < 0.05).sum(0)))))


if __name__ == "__main__":
    main()
