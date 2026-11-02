#!/usr/bin/env python
"""
Metabolic signature -> drug-target prediction for one cell line (A549).

Adaptation of the Menacho et al. 2026 framework (SIGN features -> VAE -> FNN ensemble -> BES drug ranking): the
DETF input is replaced by the significantly changed metabolites (MetChange) or reactions (METAflux) of each
perturbation, propagated through the network built by build_network.py.

Differences from the original code (DL_drug_repurposing/model_library.py), all deliberate:
  * up- and down-stream passes really use A and A^T (the original builds adj2 identical to adj)
  * the signal has two channels (increased, decreased) instead of act=1 / inh=2 in one channel
  * each relation is normalised separately and then mixed (--alpha), so the large PPI and ChEA layers cannot
    drown the metabolic layer, while propagation can still cross layers
  * with a single cell line the expression-gated graph is identical for all perturbations, so graph operators,
    distances and static features are computed once and perturbations are batched
  * Louvain labels and the dense n x n cosine-similarity feature are dropped
  * PolynomialFeatures + min-max are applied on the fly (torch) instead of being materialised
  * the split is by drug, pemp in the BES is computed analytically (the probability that a random prediction set
    of the same size hits the drug), and the evaluation metric counts FP / FN correctly

Inputs
  --nodes / --edges        build_network.py output (data/networks/combined)
  --baseline               prepare_inputs.py output: symbol, value (A549 control mean, HGNC symbols)
  --perturbations          prepare_inputs.py output: perturbation, pertname, pert_ids, dose, timepoint
  --signatures             features x perturbations matrix (compute_labels log2FC.csv); rows are GEM metabolite
                           (MetChange) or reaction (METAflux) ids, columns perturbation names as in mdat.csv
  --qvalues                same shape (qvals.csv); without it the --top-n features by |value| are used
  --signature-sign         multiply the signatures by this (MetChange inconsistency: -1 gives positive = more
                           production, the convention of the Schuhknecht benchmark)
  --drug-targets           build_drug_targets.py output (drug = LINCS pert_id, target, effect, pertname)
  --drug-metabolites       build_network.py drug_metabolite_matches.tsv: metabolites removed from a drug's own
                           signature (the drug or its metabolites measured as a metabolite)

Drugs are labelled by name (pertname), pooling the targets of all their BRD ids. Drugs without a known target are
not used for training or evaluation (no annotation is not evidence of no target); --predict-unlabeled scores them.

Outputs (--out-dir): config.json, feature_names.txt, models.pt, metrics.json, test_predicted_targets.tsv,
  test_drug_ranks.tsv, and with --predict-unlabeled unlabeled_predicted_targets.tsv
"""
import argparse
import json
import os
import re
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse import csgraph
from scipy.stats import hypergeom, norm
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import average_precision_score, f1_score
from sklearn.preprocessing import StandardScaler

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset, WeightedRandomSampler

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build_network"))
from build_network import SUGGESTED_ALPHA  # noqa: E402

# 0 = no effect, 1 = activation, 2 = inhibition, 3 = unknown direction (wildcard, not trained on)
EFFECT_CLASS = {"act": 1, "inh": 2, "unknown": 3}

# Final three FNN configurations of the paper (lr 0.005, dropout 0, 1500 initial units, 5 layers). The paper's
# second config uses L1 = 0.2 / L2 = 0.1; with L1 implemented as a sum over all weights (as in FFN_decay) that
# term dominates the loss at this model size, so it is off here.
ENSEMBLE_CONFIGS = [dict(batch_size=1024, l1=0.0, l2=0.0),
                    dict(batch_size=2048, l1=0.0, l2=1e-4),
                    dict(batch_size=2048, l1=0.0, l2=0.0)]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def pick_device(name):
    if name != "auto":
        return torch.device(name)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


# --------------------------------------------------------------------------------------
# Graph and node activity
# --------------------------------------------------------------------------------------
class Graph:
    def __init__(self, nodes_path, edges_path):
        self.nodes = pd.read_csv(nodes_path, sep="\t", dtype=str, keep_default_na=False)
        e = pd.read_csv(edges_path, sep="\t", keep_default_na=False, dtype={"src": str, "dst": str})
        self.N = len(self.nodes)
        self.pos = pd.Series(np.arange(self.N), index=self.nodes.node_id)
        e = e[e.src.isin(self.pos.index) & e.dst.isin(self.pos.index)].reset_index(drop=True)
        self.edges = e
        self.s, self.d = self.pos[e.src].to_numpy(), self.pos[e.dst].to_numpy()
        self.prot = np.where(self.nodes.node_type == "protein")[0]
        self.prot_symbols = self.nodes.node_id.to_numpy()[self.prot]
        self.prot_pos = pd.Series(np.arange(len(self.prot)), index=self.prot_symbols)
        log(f"graph: {self.N:,} nodes ({len(self.prot):,} proteins), {len(e):,} edges, "
            f"relations {e.relation.value_counts().to_dict()}")

    def directed_arrays(self, mask=None):
        """src, dst, weight with undirected edges in both orientations."""
        m = np.ones(len(self.edges), bool) if mask is None else mask
        s, d, w = self.s[m], self.d[m], self.edges.weight.to_numpy(dtype=float)[m]
        und = self.edges.directed.to_numpy(dtype=int)[m] == 0
        return np.concatenate([s, d[und]]), np.concatenate([d, s[und]]), np.concatenate([w, w[und]])


_GPR_TOKEN = re.compile(r"\(|\)|[^\s()]+")


def eval_gpr(rule, gene_act):
    """AND -> min (complex subunits), OR -> max (isozymes); unmeasured genes count as 0."""
    toks = _GPR_TOKEN.findall(rule)
    pos = 0

    def peek():
        return toks[pos].lower() if pos < len(toks) else None

    def expr():
        nonlocal pos
        v = term()
        while peek() == "or":
            pos += 1
            v = max(v, term())
        return v

    def term():
        nonlocal pos
        v = factor()
        while peek() == "and":
            pos += 1
            v = min(v, factor())
        return v

    def factor():
        nonlocal pos
        t = toks[pos]
        pos += 1
        if t == "(":
            v = expr()
            pos += 1  # ')'
            return v
        return gene_act.get(t, 0.0)

    return expr() if toks else None


def node_activity(g, baseline_path, floor):
    """Baseline activity in [0, 1] per node; gates edge weights as in the original (a_src * a_dst). Proteins: control
    expression scaled by its maximum; reactions: their gene rule; metabolites: the median expressed protein."""
    ex = pd.read_csv(baseline_path, sep="\t", dtype={"symbol": str})
    ex = pd.Series(ex.value.to_numpy(), index=ex.symbol).groupby(level=0).max().clip(lower=0)
    ex = ex / ex.max()
    nt = g.nodes.node_type.to_numpy()
    act = np.zeros(g.N)
    act[g.prot] = ex.reindex(g.prot_symbols).fillna(0).to_numpy()
    fill = float(np.median(act[g.prot][act[g.prot] > 0])) if (act[g.prot] > 0).any() else 1.0
    gene_act = ex.to_dict()
    gpr, constitutive = g.nodes.gpr.to_numpy(), g.nodes.constitutive.to_numpy()
    for i in np.where(nt == "reaction")[0]:
        v = eval_gpr(gpr[i], gene_act) if gpr[i] and constitutive[i] != "1" else None
        act[i] = fill if v is None else v
    act[nt == "metabolite"] = fill
    log(f"activity: {np.mean(act[g.prot] > 0):.3f} of proteins expressed, metabolite / constitutive fill {fill:.3f}")
    return np.maximum(act, floor)


def propagation_operator(g, act, alphas, upstream):
    """Relation-balanced operator: each relation is mean-aggregated over its in-edges (keeping the activity gating),
    then relations are mixed with weights alpha and renormalised by the alphas present at a node.
    (Op @ x)[i] aggregates x over i's predecessors (downstream pass) or successors (upstream pass)."""
    total = sp.csr_matrix((g.N, g.N), dtype=np.float64)
    denom = np.zeros(g.N)
    for rel in sorted(g.edges.relation.unique()):
        a = alphas.get(rel, 1.0)
        if a == 0:
            continue
        s, d, w = g.directed_arrays(g.edges.relation.to_numpy() == rel)
        w = w * act[s] * act[d]
        if upstream:
            s, d = d, s
        A = sp.csr_matrix((w, (d, s)), shape=(g.N, g.N))
        cnt = np.bincount(d, minlength=g.N).astype(float)
        inv = np.divide(1.0, cnt, out=np.zeros_like(cnt), where=cnt > 0)
        total = total + a * (sp.diags(inv) @ A)
        denom += a * (cnt > 0)
    inv = np.divide(1.0, denom, out=np.zeros_like(denom), where=denom > 0)
    return (sp.diags(inv) @ total).astype(np.float32).tocsr()


def membership_matrix(g):
    """Node x set membership from GEM subsystems: reactions, their enzymes and their metabolites."""
    rows, cols, sets = [], [], {}
    rxn_sub = {}
    sub = g.nodes.subsystem.to_numpy()
    for i in np.where(g.nodes.node_type == "reaction")[0]:
        name = re.sub(r"[\[\]'\"]", "", sub[i]).strip()
        if name:
            rxn_sub[i] = sets.setdefault(name, len(sets))
            rows.append(i)
            cols.append(rxn_sub[i])
    types = g.nodes.node_type.to_numpy()
    for rel in ("catalysis", "metabolic"):
        m = g.edges.relation.to_numpy() == rel
        for s, d in zip(g.s[m], g.d[m]):
            for node, other in ((s, d), (d, s)):
                if other in rxn_sub and types[node] != "reaction":
                    rows.append(node)
                    cols.append(rxn_sub[other])
    M = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(g.N, max(len(sets), 1)))
    M.data[:] = 1.0
    log(f"membership: {len(sets)} GEM subsystems")
    return M.astype(np.float32)


# --------------------------------------------------------------------------------------
# Signatures, drugs and labels
# --------------------------------------------------------------------------------------
def feature_node_map(g, features, mode):
    """Signature feature id -> node indices: the node itself, or (MetChange) every compartment copy of a
    compartment-free metabolite id."""
    node_type = "metabolite" if mode == "metchange" else "reaction"
    typed = g.nodes[g.nodes.node_type == node_type]
    ids = set(typed.node_id)
    by_base = typed.groupby("base_met").node_id.apply(list).to_dict() if node_type == "metabolite" else {}
    fmap = {}
    for f in features:
        hit = [f] if f in ids else by_base.get(f, [])
        if hit:
            fmap[f] = g.pos[hit].to_numpy()
    return fmap


def load_signatures(g, args, perts_table):
    """features x perturbations matrices -> per perturbation: significant mapped features and their sign."""
    val = pd.read_csv(args.signatures, index_col=0)
    val.index = val.index.astype(str)
    val = val * args.signature_sign
    qv = pd.read_csv(args.qvalues, index_col=0).reindex(index=val.index, columns=val.columns) if args.qvalues else None
    fmap = feature_node_map(g, val.index, args.mode)
    log(f"signatures: {val.shape[0]:,} features x {val.shape[1]:,} perturbations; {len(fmap):,} features map to "
        f"network nodes ({val.shape[0] - len(fmap):,} do not, e.g. excluded currency metabolites)")
    val = val.loc[[f for f in val.index if f in fmap]]
    if qv is not None:
        qv = qv.loc[val.index]
    sigs, rows = {}, []
    info = perts_table.set_index("perturbation")
    for p in val.columns:
        if p not in info.index:
            continue
        v = val[p]
        keep = v.notna() & (v.abs() >= args.min_abs_score)
        if qv is not None and not args.ignore_qvalue:
            keep &= qv[p] <= args.qvalue
        else:
            keep &= v.abs().rank(ascending=False, method="first") <= args.top_n
        f = v.index[keep].to_numpy()
        if len(f):
            sigs[p] = (f, np.sign(v[keep].to_numpy()))
        rows.append(dict(pert=p, drug=info.at[p, "pertname"], dose=info.at[p, "dose"],
                         time=info.at[p, "timepoint"], n_sig=len(f)))
    perts = pd.DataFrame(rows)
    log(f"{len(perts):,} perturbations in the perturbation table; {(perts.n_sig > 0).sum():,} with significant "
        f"features (median {int(perts.n_sig[perts.n_sig > 0].median()) if (perts.n_sig > 0).any() else 0})")
    return perts[perts.n_sig > 0].reset_index(drop=True), sigs, fmap


def load_drug_targets(path, g, perts_table):
    """Targets per drug name, pooled over its BRD ids; a target reported with both directions becomes unknown."""
    dt = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    brd_to_name = {b: n.lower() for n, ids in zip(perts_table.pertname, perts_table.pert_ids) for b in ids.split("|")}
    dt["name"] = dt.drug.map(brd_to_name).fillna(dt.pertname.str.lower())
    dt = dt[dt.target.isin(g.prot_pos.index)]
    dt["cls"] = dt.effect.map(EFFECT_CLASS).fillna(3).astype(int)
    cls = dt.groupby(["name", "target"]).cls.agg(lambda c: c.iloc[0] if c.nunique() == 1 else 3)
    targets = {}
    for (name, target), c in cls.items():
        targets.setdefault(name, {})[g.prot_pos[target]] = c
    return targets


def load_drug_metabolites(path, g):
    """Drug name -> node indices to drop from that drug's own signature."""
    if not path:
        return {}
    m = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    m = m[m.node_id.isin(g.pos.index)]
    return {n.lower(): set(g.pos[ids].to_numpy()) for n, ids in m.groupby("pertname").node_id}


UNIT_UM = {"pm": 1e-6, "nm": 1e-3, "um": 1.0, "µm": 1.0, "mm": 1e3, "m": 1e6}


def dose_um(dose):
    """'10 uM', '500 nM' -> micromolar; anything unparsable -> NaN."""
    m = dose.str.strip().str.extract(r"([\d.]+)\s*([a-zA-Zµ]+)")
    return pd.to_numeric(m[0], errors="coerce") * m[1].str.lower().map(UNIT_UM)


def pert_globals(perts):
    dose = dose_um(perts.dose)
    max_dose = dose.groupby(perts.drug).transform("max")
    rel_dose = np.where((dose > 0) & (max_dose > 0), np.log10(dose / max_dose), 0.0)
    t = np.log1p(pd.to_numeric(perts.time.str.extract(r"([\d.]+)")[0], errors="coerce").fillna(0))
    return np.stack([rel_dose, t, np.log1p(perts.n_sig.to_numpy())], 1).astype(np.float32)


def split_drugs(drugs, val_frac, test_frac, seed):
    rng = np.random.default_rng(seed)
    drugs = np.array(sorted(drugs))
    rng.shuffle(drugs)
    n_test, n_val = int(round(len(drugs) * test_frac)), int(round(len(drugs) * val_frac))
    return set(drugs[n_test + n_val:]), set(drugs[n_test:n_test + n_val]), set(drugs[:n_test])


# --------------------------------------------------------------------------------------
# SIGN-style features
# --------------------------------------------------------------------------------------
class FeatureBuilder:
    def __init__(self, g, act, membership, source_nodes, args):
        self.g, self.K = g, args.hops
        alphas = dict(SUGGESTED_ALPHA)
        alphas.update({k: float(v) for k, v in (a.split("=") for a in args.alpha)})
        log(f"relation weights: {alphas}")
        self.down = propagation_operator(g, act, alphas, upstream=False)
        self.up = propagation_operator(g, act, alphas, upstream=True)
        self.M, self.Mprot = membership, membership[g.prot]
        self._static(act)
        self._distances(source_nodes)
        K = self.K
        dyn = [f"{d}_{c}_hop{k + 1}" for d in ("down", "up") for c in ("inc", "dec") for k in range(K)]
        diff = [f"{d}_{c}_diff{k + 2}" for d in ("down", "up") for c in ("inc", "dec") for k in range(K - 1)]
        slog = [f"{d}_{c}_slog_last" for d in ("down", "up") for c in ("inc", "dec")]
        stats = [f"dyn_{s}" for s in ("sum", "mean", "max", "min", "median", "std", "q25", "q75")]
        dist = [f"{d}_dist_{s}" for d in ("up", "down") for s in ("min", "mean", "frac_within_K")]
        self.names = (dyn + diff + slog + stats + self.static_names + dist
                      + ["shared_subsystems_with_sources", "rel_log_dose", "log_time", "log_n_sig"])

    def _static(self, act):
        g = self.g
        feats, names = [act[g.prot]], ["activity"]
        for name, op in (("down", self.down), ("up", self.up)):
            x = act.astype(np.float32)
            for k in range(self.K):
                x = op @ x
                feats.append(x[g.prot])
                names.append(f"activity_{name}_hop{k + 1}")
        s, d, w = g.directed_arrays()
        w = w * act[s] * act[d]
        feats += [np.log1p(np.bincount(d, w, g.N)[g.prot]), np.log1p(np.bincount(s, w, g.N)[g.prot]),
                  np.log1p(np.bincount(d, minlength=g.N)[g.prot]), np.log1p(np.bincount(s, minlength=g.N)[g.prot])]
        names += ["in_strength", "out_strength", "in_degree", "out_degree"]
        self.static = np.stack(feats, 1).astype(np.float32)
        self.static_names = names

    def _distances(self, source_nodes):
        """Hop distances between every protein and every node that acts as a signature source."""
        g = self.g
        s, d, _ = g.directed_arrays()
        G = sp.csr_matrix((np.ones(len(s)), (s, d)), shape=(g.N, g.N))
        G.data[:] = 1
        GT = G.T.tocsr()
        self.src_row = pd.Series(np.arange(len(source_nodes)), index=source_nodes)
        self.cap = 4 * self.K
        self.D_up = np.empty((len(source_nodes), len(g.prot)), np.uint8)
        self.D_down = np.empty_like(self.D_up)
        t0 = time.time()
        for start in range(0, len(source_nodes), 256):
            idx = source_nodes[start:start + 256]
            for D, M in ((self.D_up, GT), (self.D_down, G)):
                dist = csgraph.dijkstra(M, directed=True, indices=idx, unweighted=True)[:, g.prot]
                D[start:start + len(idx)] = np.minimum(dist, self.cap)
        log(f"distances from {len(source_nodes):,} source nodes in {time.time() - t0:.0f}s")

    def build(self, sources, globals_):
        """sources: list of (node_idx, weight) per perturbation -> (P * n_prot, F) float32, perturbation-major."""
        g, K, P, n = self.g, self.K, len(sources), len(self.g.prot)
        S = np.zeros((g.N, P, 2), np.float32)
        for j, (idx, w) in enumerate(sources):
            np.add.at(S[:, j, 0], idx[w > 0], w[w > 0])
            np.add.at(S[:, j, 1], idx[w < 0], -w[w < 0])
        dyn = np.empty((P, n, 2, 2, K), np.float32)
        for di, op in enumerate((self.down, self.up)):
            X = S.reshape(g.N, P * 2)
            for k in range(K):
                X = op @ X
                dyn[:, :, di, :, k] = X[g.prot].reshape(n, P, 2).transpose(1, 0, 2)
        diff = dyn[..., 1:] - dyn[..., :-1]
        last = dyn[..., -1]
        slog = np.sign(last) * np.log(np.abs(last) + 1e-6)
        flat = dyn.reshape(P, n, -1)
        stats = np.stack([flat.sum(-1), flat.mean(-1), flat.max(-1), flat.min(-1), np.median(flat, -1),
                          flat.std(-1), np.quantile(flat, 0.25, -1), np.quantile(flat, 0.75, -1)], -1)
        dist = np.empty((P, n, 6), np.float32)
        for j, (idx, _) in enumerate(sources):
            for o, D in ((0, self.D_up), (3, self.D_down)):
                rows = D[self.src_row[np.unique(idx)].to_numpy()].astype(np.float32)
                dist[j, :, o] = rows.min(0)
                dist[j, :, o + 1] = rows.mean(0)
                dist[j, :, o + 2] = (rows <= K).mean(0)
        shared = np.asarray(self.Mprot @ (self.M.T @ np.abs(S).sum(-1))).T[..., None]
        X = np.concatenate([flat, diff.reshape(P, n, -1), slog.reshape(P, n, -1), stats,
                            np.broadcast_to(self.static, (P,) + self.static.shape), dist, shared,
                            np.broadcast_to(globals_[:, None, :], (P, n, globals_.shape[1]))], -1)
        return np.ascontiguousarray(X.reshape(P * n, -1), dtype=np.float32)


def select_sources(features, signs, fmap, drop_nodes, max_sources, rng):
    """Random cap on the number of signature features (paper: max 100 DETFs). A feature mapped to n nodes gives
    weight +-1/n to each. drop_nodes: the drug's own metabolites."""
    if len(features) > max_sources:
        keep = rng.choice(len(features), max_sources, replace=False)
        features, signs = features[keep], signs[keep]
    idx, w = [], []
    for f, s in zip(features, signs):
        nodes = np.array([x for x in fmap[f] if x not in drop_nodes])
        if len(nodes):
            idx.append(nodes)
            w.append(np.full(len(nodes), s / len(fmap[f])))
    if not idx:
        return np.empty(0, int), np.empty(0)
    return np.concatenate(idx), np.concatenate(w)


def label_vector(targets, n):
    y = np.zeros(n, np.int64)
    for i, c in (targets or {}).items():
        y[i] = c
    return y


def iter_chunks(fb, perts, sigs, fmap, targets, drug_mets, globals_, args):
    n = len(fb.g.prot)
    for start in range(0, len(perts), args.chunk):
        chunk = perts.iloc[start:start + args.chunk]
        sources = []
        for j, p, d in zip(chunk.index, chunk.pert, chunk.drug):
            rng = np.random.default_rng(args.seed + j)
            sources.append(select_sources(*sigs[p], fmap, drug_mets.get(d.lower(), set()), args.max_sources, rng))
        X = fb.build(sources, globals_[chunk.index.to_numpy()])
        y = np.concatenate([label_vector(targets.get(d.lower()), n) for d in chunk.drug])
        yield chunk, X, y


def sample_rows(X, y, n_perts, args, rng):
    """All target rows + cluster-stratified negatives (MiniBatchKMeans, equal number per cluster)."""
    pos = np.where(y > 0)[0]
    neg = np.where(y == 0)[0]
    n_take = min(args.neg_per_pert * n_perts, len(neg))
    sub = rng.choice(neg, min(len(neg), 50000), replace=False)
    n_clusters = min(args.kmeans_clusters, len(sub))
    scaler = StandardScaler().fit(X[sub])
    km = MiniBatchKMeans(n_clusters, random_state=args.seed, n_init=3, batch_size=4096)
    km.fit(np.nan_to_num(scaler.transform(X[sub])))
    lab = km.predict(np.nan_to_num(scaler.transform(X[neg])))
    per = int(np.ceil(n_take / n_clusters))
    chosen = [rng.choice(neg[lab == c], min(per, (lab == c).sum()), replace=False)
              for c in range(n_clusters) if (lab == c).any()]
    chosen = np.concatenate(chosen)
    if len(chosen) < n_take:  # top up from what the small clusters left
        rest = np.setdiff1d(neg, chosen)
        chosen = np.concatenate([chosen, rng.choice(rest, n_take - len(chosen), replace=False)])
    return np.concatenate([pos, chosen[:n_take]])


# --------------------------------------------------------------------------------------
# Models (architectures mirror VAE_decrease and FFN_decay from model_library.py)
# --------------------------------------------------------------------------------------
class PolyMinMax(nn.Module):
    """PolynomialFeatures(degree=2, include_bias=False) + column min-max + removal of constant columns."""

    def __init__(self, n_features, poly=True):
        super().__init__()
        self.poly = poly
        i, j = torch.triu_indices(n_features, n_features) if poly else (torch.empty(0, dtype=torch.long),) * 2
        self.register_buffer("i", i)
        self.register_buffer("j", j)
        out = n_features + len(i)
        self.register_buffer("mi", torch.zeros(out))
        self.register_buffer("rng", torch.ones(out))
        self.register_buffer("keep", torch.ones(out, dtype=torch.bool))

    def expand(self, x):
        return torch.cat([x, x[:, self.i] * x[:, self.j]], 1) if self.poly else x

    @torch.no_grad()
    def fit(self, X, device, batch=65536):
        mi = ma = None
        for s in range(0, len(X), batch):
            z = self.expand(torch.from_numpy(X[s:s + batch]).to(device))
            bmi, bma = z.min(0).values, z.max(0).values
            mi = bmi if mi is None else torch.minimum(mi, bmi)
            ma = bma if ma is None else torch.maximum(ma, bma)
        self.mi.copy_(mi)
        self.rng.copy_((ma - mi).clamp(min=1e-8))
        self.keep.copy_((ma - mi) > 1e-8)
        log(f"polynomial expansion: {len(mi):,} columns, {int(self.keep.sum()):,} non-constant")
        return self

    def forward(self, x):
        z = (self.expand(x) - self.mi) / self.rng
        return z[:, self.keep].clamp(0, 1)


class VAE_decrease(nn.Module):
    def __init__(self, in_channels, initial_hidden_size, num_layers, dropout, latent_dimension):
        super().__init__()
        self.encoder_sizes = self._generate_hidden_sizes(initial_hidden_size, num_layers)
        self.decoder_sizes = self.encoder_sizes[::-1]
        self.encoder_layers = nn.ModuleList()
        prev = in_channels
        for h in self.encoder_sizes:
            self.encoder_layers.append(self._build_layer(prev, h, dropout, batch_norm=True))
            prev = h
        self.fc_mu = nn.Linear(self.encoder_sizes[-1], latent_dimension)
        self.fc_logvar = nn.Linear(self.encoder_sizes[-1], latent_dimension)
        self.decoder_layers = nn.ModuleList()
        prev = latent_dimension
        for h in self.decoder_sizes:
            self.decoder_layers.append(self._build_layer(prev, h, dropout, batch_norm=False))
            prev = h
        self.output_layer = nn.Linear(self.decoder_sizes[-1], in_channels)

    def _generate_hidden_sizes(self, initial_size, num_layers):
        sizes = [initial_size]
        largest_power_of_2 = 2 ** (initial_size.bit_length() - 1)
        if largest_power_of_2 == initial_size:
            largest_power_of_2 //= 2
        sizes.append(largest_power_of_2)
        size = largest_power_of_2 // 2
        for _ in range(num_layers - 2):
            sizes.append(size)
            size = max(size // 2, 8)
        return sizes

    def _build_layer(self, in_features, out_features, dropout, batch_norm):
        layers = [nn.Linear(in_features, out_features)]
        if batch_norm:
            layers.append(nn.BatchNorm1d(out_features))
        layers.append(nn.ReLU())
        if dropout > 0:
            layers.append(nn.Dropout(dropout))
        return nn.Sequential(*layers)

    def encode(self, x):
        for layer in self.encoder_layers:
            x = layer(x)
        return self.fc_mu(x), torch.clamp(self.fc_logvar(x), min=-5, max=5)

    def decode(self, z):
        for layer in self.decoder_layers:
            z = layer(z)
        return torch.sigmoid(self.output_layer(z))

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)
        return self.decode(z), mu, logvar


class FFN_decay(nn.Module):
    def __init__(self, in_channels, initial_size, out_channels, num_layers, dropout=0.0, decay_factor=0.7):
        super().__init__()
        layers, prev, size = [], in_channels, initial_size
        for _ in range(num_layers):
            layers.append(nn.Sequential(nn.Linear(prev, size), nn.BatchNorm1d(size), nn.LeakyReLU(),
                                        nn.Dropout(dropout)))
            prev, size = size, max(int(size * decay_factor), 8)
        layers.append(nn.Linear(prev, out_channels))
        self.layers = nn.Sequential(*layers)

    def forward(self, x):
        return self.layers(x)


class Preprocessor(nn.Module):
    """Frozen input pipeline: polynomial + min-max, then (optionally) VAE decode(encode_mu(x))."""

    def __init__(self, poly, vae=None):
        super().__init__()
        self.poly, self.vae = poly, vae

    @torch.no_grad()
    def forward(self, x):
        z = self.poly(x)
        if self.vae is not None:
            z = self.vae.decode(self.vae.encode(z)[0])
        return z


def train_vae(pre_poly, X, args, device, rng):
    rows = rng.choice(len(X), min(len(X), args.vae_max_rows), replace=False)
    hold = rows[: max(1, len(rows) // 10)]  # 10% held out, as in the paper
    train = rows[len(hold):]
    in_dim = int(pre_poly.keep.sum())
    vae = VAE_decrease(in_dim, args.vae_hidden, args.vae_layers, 0.0, args.latent).to(device)
    opt = torch.optim.Adam(vae.parameters(), lr=args.vae_lr)
    Xt = torch.from_numpy(X)
    for epoch in range(args.vae_epochs):
        vae.train()
        perm = rng.permutation(train)
        tot = 0.0
        for s in range(0, len(perm), args.vae_batch):
            x = pre_poly(Xt[perm[s:s + args.vae_batch]].to(device))
            if len(x) < 2:
                continue
            recon, mu, logvar = vae(x)
            rec = F.mse_loss(recon, x, reduction="sum") / len(x)
            kl = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp()) / len(x)
            loss = rec + args.vae_beta * kl
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += loss.item() * len(x)
        vae.eval()
        with torch.no_grad():
            x = pre_poly(Xt[hold].to(device))
            val = F.mse_loss(vae.decode(vae.encode(x)[0]), x).item()
        log(f"VAE epoch {epoch + 1}: train loss {tot / len(train):.3f}, held-out MSE {val:.5f}")
    return vae.eval()


def evaluate_ffn(y_true, y_pred):
    """Class 3 (unknown direction) counts as correct if predicted 1 or 2, as in Menacho et al. The original
    evaluate_ffn sets FP = FN = total - TP (so TN < 0 and its F1 equals accuracy); here FPR = non-targets
    predicted as targets, FNR = targets predicted as 0 or with the wrong direction."""
    y_true, y_pred = np.asarray(y_true).copy(), np.asarray(y_pred)
    wild = y_true == 3
    hit = wild & ((y_pred == 1) | (y_pred == 2))
    y_true[hit] = y_pred[hit]
    y_true[wild & (y_pred == 0)] = 1  # counted as a missed target
    neg, pos = y_true == 0, y_true > 0
    return {"accuracy": float((y_true == y_pred).mean()),
            "fpr": float((y_pred[neg] > 0).mean()) if neg.any() else None,
            "fnr": float((y_pred[pos] != y_true[pos]).mean()) if pos.any() else None,
            "macro_f1_targets": float(f1_score(y_true, y_pred, labels=[1, 2], average="macro", zero_division=0))}


def node_metrics(y, prob):
    m = evaluate_ffn(y, prob.argmax(1))
    m["auprc_target_vs_none"] = float(average_precision_score(y > 0, 1 - prob[:, 0])) if (y > 0).any() else None
    m["target_rate"] = float((y > 0).mean())
    return m


@torch.no_grad()
def predict_proba(models, pre, X, device, batch=65536):
    out = []
    for s in range(0, len(X), batch):
        z = pre(torch.from_numpy(X[s:s + batch]).to(device))
        out.append(torch.stack([F.softmax(m(z), 1) for m in models]).cpu())
    return torch.cat(out, 1).numpy()  # (n_models, n_rows, 3)


def majority_vote(probs):
    """Mode of the models' argmax; ties (three different classes) go to the class with the highest mean probability."""
    votes = probs.argmax(-1)  # (n_models, n_rows)
    counts = np.stack([(votes == c).sum(0) for c in range(probs.shape[-1])], -1)
    best = counts.max(-1, keepdims=True)
    return np.where(counts == best, probs.mean(0), -1).argmax(-1)


def train_ffn(cfg, pre, in_dim, Xtr, ytr, wtr, Xva, yva, args, device, seed):
    torch.manual_seed(seed)
    model = FFN_decay(in_dim, args.hidden, 3, args.layers, args.dropout, args.decay).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=cfg["l2"])
    sched = torch.optim.lr_scheduler.StepLR(opt, step_size=10, gamma=0.5)
    sampler = WeightedRandomSampler(torch.from_numpy(wtr), num_samples=len(wtr), replacement=True,
                                    generator=torch.Generator().manual_seed(seed))
    loader = DataLoader(TensorDataset(torch.from_numpy(Xtr), torch.from_numpy(ytr)),
                        batch_size=cfg["batch_size"], sampler=sampler, drop_last=True)
    best, best_state, bad = -1.0, None, 0
    for epoch in range(args.ffn_epochs):
        model.train()
        tot = 0.0
        for xb, yb in loader:
            logits = model(pre(xb.to(device)))
            loss = F.cross_entropy(logits, yb.to(device))
            if cfg["l1"] > 0:
                loss = loss + cfg["l1"] * sum(p.abs().sum() for p in model.parameters())
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += loss.item()
        sched.step()
        model.eval()
        m = node_metrics(yva, predict_proba([model], pre, Xva, device)[0])
        log(f"  FFN {cfg} epoch {epoch + 1}: loss {tot / max(len(loader), 1):.4f} | val "
            + ", ".join(f"{k} {v:.4f}" for k, v in m.items() if v is not None))
        score = m["macro_f1_targets"]
        if score > best:
            best, bad = score, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            bad += 1
            if bad >= args.patience:
                break
    model.load_state_dict(best_state)
    return model.eval()


# --------------------------------------------------------------------------------------
# Drug ranking (Bayesian enrichment score, as in rank_drugs_by_BES.py / compute_drug_pvalues.R)
# --------------------------------------------------------------------------------------
def compute_BES_score(pemp, hypergeo_pval, rank_hy, pi0=0.99, w_posterior=20.0, w_hypergeo=1.0,
                      w_hypergeo_rank=1.0, pemp_penalty_multiplier=1.0):
    pemp, hypergeo_pval, rank_hy = (np.asarray(a, dtype=float) for a in (pemp, hypergeo_pval, rank_hy))
    mask_zero = hypergeo_pval == 1.0
    posterior = 1 - pi0 * pemp
    s_posterior = posterior * np.exp(-pemp_penalty_multiplier * (1 - pemp))
    rank_hy_norm = 1 - (rank_hy - np.min(rank_hy)) / (np.max(rank_hy) - np.min(rank_hy) + 1e-8)
    BES = w_posterior * s_posterior + w_hypergeo * (1 - hypergeo_pval) + w_hypergeo_rank * rank_hy_norm
    BES[mask_zero] = 0.0
    BES_rank = pd.Series(BES).rank(ascending=False, method="min").to_numpy(copy=True)  # writable under pandas 3
    BES_rank[BES == 0.0] = len(BES)
    return BES, BES_rank


class DrugRanker:
    """Universe of directional tokens: 2 * gene + (class - 1); unknown-direction targets give both tokens."""

    def __init__(self, targets, n_prot):
        self.drugs = sorted(targets)
        self.U = 2 * n_prot
        rows, cols = [], []
        for j, d in enumerate(self.drugs):
            for gene, c in targets[d].items():
                for t in ([2 * gene, 2 * gene + 1] if c == 3 else [2 * gene + c - 1]):
                    rows.append(t)
                    cols.append(j)
        self.DM = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(self.U, len(self.drugs)))
        self.m = np.asarray(self.DM.sum(0)).ravel()
        self.col = {d: j for j, d in enumerate(self.drugs)}

    def bes_ranks(self, vote):
        tok = np.concatenate([2 * np.where(vote == 1)[0], 2 * np.where(vote == 2)[0] + 1])
        k = len(tok)
        if k == 0:
            return np.full(len(self.drugs), len(self.drugs), float)
        overlap = np.asarray(self.DM[tok].sum(0)).ravel()
        hyp = np.where(overlap > 0, hypergeom.sf(overlap - 1, self.U, self.m, k), 1.0)
        pemp = 1 - hypergeom.pmf(0, self.U, self.m, k)  # P(random set of size k hits the drug)
        z = np.where(hyp == 0, np.inf, -norm.ppf(np.clip(hyp, 1e-300, 1) / 2))
        rank_hy = pd.Series(-z).rank(method="min").to_numpy()
        return compute_BES_score(pemp, hyp, rank_hy)[1]

    def prob_scores(self, prob):
        """Mean predicted probability of each drug's targets with the annotated direction (unknown-direction
        targets contribute (p_act + p_inh) / 2)."""
        tokp = np.empty(self.U)
        tokp[0::2], tokp[1::2] = prob[:, 1], prob[:, 2]
        return np.asarray(self.DM.T @ tokp).ravel() / self.m

    def rank_of(self, scores, drug):
        return float(pd.Series(-scores).rank(method="average").to_numpy()[self.col[drug]])


def rank_summary(ranks, n):
    r = np.asarray(ranks, float)
    return {"median_rank": float(np.median(r)), "mean_rank": float(r.mean()), "n_candidates": n,
            **{f"top{k}": float((r <= k).mean()) for k in (1, 5, 10, 50)}}


def predicted_rows(chunk, votes, mean_p, prot_symbols, targets, top_k):
    """Called targets per perturbation, the top_k by probability of the called direction."""
    rows = []
    for j, (p, d) in enumerate(zip(chunk.pert, chunk.drug)):
        v, pr = votes[j], mean_p[j]
        known = targets.get(d.lower(), {})
        called = np.where(v > 0)[0]
        called = called[np.argsort(-pr[called, v[called]])][:top_k]
        for gn in called:
            rows.append(dict(pert=p, drug=d, target=prot_symbols[gn], effect="act" if v[gn] == 1 else "inh",
                             prob=float(pr[gn, v[gn]]), is_known_target=gn in known))
    return rows


# --------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    io = ap.add_argument_group("inputs")
    io.add_argument("--nodes", required=True)
    io.add_argument("--edges", required=True)
    io.add_argument("--baseline", required=True)
    io.add_argument("--perturbations", required=True)
    io.add_argument("--signatures", required=True)
    io.add_argument("--qvalues")
    io.add_argument("--signature-sign", type=float, default=1.0)
    io.add_argument("--mode", choices=["metchange", "metaflux"], default="metchange")
    io.add_argument("--drug-targets", required=True)
    io.add_argument("--drug-metabolites")
    io.add_argument("--out-dir", required=True)
    fe = ap.add_argument_group("features")
    fe.add_argument("--hops", type=int, default=5)
    fe.add_argument("--max-sources", type=int, default=100)
    fe.add_argument("--qvalue", type=float, default=0.05)
    fe.add_argument("--ignore-qvalue", action="store_true", help="use the top-n features by |value| instead")
    fe.add_argument("--top-n", type=int, default=100)
    fe.add_argument("--min-abs-score", type=float, default=0.0)
    fe.add_argument("--min-sources", type=int, default=3, help="drop weaker perturbations")
    fe.add_argument("--alpha", action="append", default=[], metavar="REL=W",
                    help="relation mixing weight (defaults: build_network.SUGGESTED_ALPHA)")
    fe.add_argument("--activity-floor", type=float, default=0.0)
    fe.add_argument("--no-poly", action="store_true")
    sa = ap.add_argument_group("sampling / split")
    sa.add_argument("--neg-per-pert", type=int, default=200)
    sa.add_argument("--kmeans-clusters", type=int, default=50)
    sa.add_argument("--chunk", type=int, default=32, help="perturbations per feature batch")
    sa.add_argument("--val-frac", type=float, default=0.1)
    sa.add_argument("--test-frac", type=float, default=0.1)
    sa.add_argument("--max-perts", type=int, default=0, help="subsample perturbations (smoke tests)")
    sa.add_argument("--max-test-perts", type=int, default=0, help="cap full-proteome test scoring")
    sa.add_argument("--predict-unlabeled", action="store_true", help="also score drugs without known targets")
    sa.add_argument("--report-top-k", type=int, default=200,
                    help="called targets written per perturbation (the drug ranking uses all of them)")
    mo = ap.add_argument_group("models")
    mo.add_argument("--no-vae", action="store_true", help="ablation: FNN on the min-max features")
    mo.add_argument("--vae-hidden", type=int, default=1500)
    mo.add_argument("--vae-layers", type=int, default=5)
    mo.add_argument("--latent", type=int, default=100)
    mo.add_argument("--vae-lr", type=float, default=1e-3)
    mo.add_argument("--vae-beta", type=float, default=1e-3)
    mo.add_argument("--vae-batch", type=int, default=1024)
    mo.add_argument("--vae-epochs", type=int, default=20)
    mo.add_argument("--vae-max-rows", type=int, default=500000)
    mo.add_argument("--hidden", type=int, default=1500)
    mo.add_argument("--layers", type=int, default=5)
    mo.add_argument("--decay", type=float, default=0.7)
    mo.add_argument("--dropout", type=float, default=0.0)
    mo.add_argument("--lr", type=float, default=0.005)
    mo.add_argument("--ffn-epochs", type=int, default=30)
    mo.add_argument("--patience", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--device", default="auto", help="auto (cuda > mps > cpu), cuda, mps or cpu")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, "config.json"), "w") as fh:
        json.dump(vars(args), fh, indent=2)
    rng = np.random.default_rng(args.seed)
    torch.manual_seed(args.seed)
    device = pick_device(args.device)
    log(f"device: {device}")

    # ---- data
    g = Graph(args.nodes, args.edges)
    act = node_activity(g, args.baseline, args.activity_floor)
    perts_table = pd.read_csv(args.perturbations, sep="\t", dtype=str, keep_default_na=False)
    perts, sigs, fmap = load_signatures(g, args, perts_table)
    targets = load_drug_targets(args.drug_targets, g, perts_table)
    drug_mets = load_drug_metabolites(args.drug_metabolites, g)
    perts = perts[perts.n_sig >= args.min_sources].reset_index(drop=True)
    labeled = perts.drug.str.lower().isin(targets.keys())
    log(f"{len(perts):,} perturbations with >= {args.min_sources} features: {labeled.sum():,} of "
        f"{perts.drug[labeled].nunique()} drugs with known targets, {(~labeled).sum():,} of "
        f"{perts.drug[~labeled].nunique()} drugs without (not used for training or evaluation)")
    unlabeled = perts[~labeled].reset_index(drop=True)
    perts = perts[labeled].reset_index(drop=True)
    if args.max_perts:
        perts = perts.sample(min(args.max_perts, len(perts)), random_state=args.seed).reset_index(drop=True)
    globals_ = pert_globals(perts)
    tr_d, va_d, te_d = split_drugs(perts.drug.unique(), args.val_frac, args.test_frac, args.seed)
    split = perts.drug.map(lambda d: "train" if d in tr_d else "val" if d in va_d else "test")
    log(f"split by drug: {split.value_counts().to_dict()} perturbations")

    used = np.unique(np.concatenate([fmap[f] for p in pd.concat([perts.pert, unlabeled.pert]) for f in sigs[p][0]]))
    fb = FeatureBuilder(g, act, membership_matrix(g), used, args)
    with open(os.path.join(args.out_dir, "feature_names.txt"), "w") as fh:
        fh.write("\n".join(fb.names) + "\n")
    log(f"{len(fb.names)} base features per protein")

    # ---- sampled training / validation rows
    data = {}
    for part in ("train", "val"):
        sub = perts[split == part]
        Xs, ys, ps = [], [], []
        for chunk, X, y in iter_chunks(fb, sub, sigs, fmap, targets, drug_mets, globals_, args):
            keep = sample_rows(X, y, len(chunk), args, rng)
            Xs.append(X[keep])
            ys.append(y[keep])
            ps.append(chunk.index.to_numpy()[keep // len(g.prot)])
        data[part] = (np.concatenate(Xs), np.concatenate(ys), np.concatenate(ps))
        log(f"{part}: {len(data[part][1]):,} rows, class counts {np.bincount(data[part][1], minlength=4).tolist()}")

    Xtr, ytr, ptr = data["train"]
    poly = PolyMinMax(Xtr.shape[1], poly=not args.no_poly).to(device).fit(Xtr, device)
    vae = None if args.no_vae else train_vae(poly, Xtr, args, device, rng)
    pre = Preprocessor(poly, vae).to(device).eval()
    in_dim = int(poly.keep.sum())

    # class 3 is not trained on; the sampler balances classes and drugs (many doses of one drug should not
    # outweigh a drug measured once)
    trainable = ytr < 3
    Xtr, ytr, ptr = Xtr[trainable], ytr[trainable], ptr[trainable]
    class_w = 1.0 / np.bincount(ytr, minlength=3).clip(min=1)
    drug_n = perts.drug.map(perts.drug.value_counts()).to_numpy()
    w = (class_w[ytr] / drug_n[ptr]).astype(np.float64)
    Xva, yva, _ = data["val"]

    models = [train_ffn(cfg, pre, in_dim, Xtr, ytr, w, Xva, yva, args, device, args.seed + i)
              for i, cfg in enumerate(ENSEMBLE_CONFIGS)]
    probs_va = predict_proba(models, pre, Xva, device)
    metrics = {"val_sampled_rows": {f"model{i + 1}": node_metrics(yva, probs_va[i]) for i in range(len(models))}}
    metrics["val_sampled_rows"]["ensemble_vote"] = evaluate_ffn(yva, majority_vote(probs_va))
    torch.save({"poly": poly.state_dict(), "vae": None if vae is None else vae.state_dict(),
                "ffns": [m.state_dict() for m in models], "feature_names": fb.names, "args": vars(args)},
               os.path.join(args.out_dir, "models.pt"))

    # ---- full-proteome scoring of held-out drugs
    test = perts[split == "test"]
    if args.max_test_perts:
        test = test.sample(min(args.max_test_perts, len(test)), random_state=args.seed)
    tkeys = {d: targets[d] for d in sorted({d.lower() for d in perts.drug})}  # candidates: all labelled A549 drugs
    ranker = DrugRanker(tkeys, len(g.prot))
    train_drugs = sorted({d.lower() for d in tr_d})
    prior = np.zeros((len(g.prot), 3))
    for d in train_drugs:
        for gene, c in targets[d].items():
            if c == 3:
                prior[gene, 1:] += 0.5
            else:
                prior[gene, c] += 1
    prior /= len(train_drugs)
    prior_scores = ranker.prob_scores(prior)

    rank_rows, pred_rows, ys_all, probs_all, drug_prob = [], [], [], [], {}
    for chunk, X, y in iter_chunks(fb, test, sigs, fmap, targets, drug_mets, globals_, args):
        P = predict_proba(models, pre, X, device)
        votes = majority_vote(P).reshape(len(chunk), -1)
        mean_p = P.mean(0).reshape(len(chunk), len(g.prot), 3)
        ys_all.append(y)
        probs_all.append(P.mean(0).astype(np.float32))
        for j, (p, d) in enumerate(zip(chunk.pert, chunk.drug)):
            key, v, pr = d.lower(), votes[j], mean_p[j]
            drug_prob.setdefault(key, []).append(pr)
            known = targets[key]
            rank_rows.append(dict(
                pert=p, drug=d, n_predicted=int((v > 0).sum()), n_true_targets=len(known),
                true_targets_hit=int(sum(v[gn] in ((1, 2) if c == 3 else (c,)) for gn, c in known.items())),
                rank_bes=float(ranker.bes_ranks(v)[ranker.col[key]]),
                rank_prob=ranker.rank_of(ranker.prob_scores(pr), key),
                rank_baseline=ranker.rank_of(prior_scores, key)))
        pred_rows += predicted_rows(chunk, votes, mean_p, g.prot_symbols, targets, args.report_top_k)
        log(f"test scoring: {len(rank_rows)}/{len(test)} perturbations")

    ranks = pd.DataFrame(rank_rows)
    ranks.to_csv(os.path.join(args.out_dir, "test_drug_ranks.tsv"), sep="\t", index=False)
    pd.DataFrame(pred_rows).to_csv(os.path.join(args.out_dir, "test_predicted_targets.tsv"), sep="\t", index=False)
    yt, pt = np.concatenate(ys_all), np.concatenate(probs_all)
    nd = len(ranker.drugs)
    drug_ranks = [ranker.rank_of(ranker.prob_scores(np.mean(ps, 0)), d) for d, ps in drug_prob.items()]
    metrics["test_full_proteome_nodes"] = node_metrics(yt, pt)
    metrics["test_per_perturbation"] = {k: rank_summary(ranks[k], nd) for k in ("rank_bes", "rank_prob", "rank_baseline")}
    metrics["test_per_drug"] = {"rank_prob_mean_over_doses": rank_summary(drug_ranks, nd),
                                "rank_baseline": rank_summary([ranker.rank_of(prior_scores, d) for d in drug_prob], nd)}

    # ---- drugs without known targets (application, no evaluation possible)
    if args.predict_unlabeled and len(unlabeled):
        ug = pert_globals(unlabeled)
        rows = []
        for chunk, X, _ in iter_chunks(fb, unlabeled, sigs, fmap, {}, drug_mets, ug, args):
            P = predict_proba(models, pre, X, device)
            rows += predicted_rows(chunk, majority_vote(P).reshape(len(chunk), -1),
                                   P.mean(0).reshape(len(chunk), len(g.prot), 3), g.prot_symbols, {},
                                   args.report_top_k)
        pd.DataFrame(rows).to_csv(os.path.join(args.out_dir, "unlabeled_predicted_targets.tsv"), sep="\t", index=False)
        log(f"unlabeled drugs: {len(rows):,} predicted targets for {unlabeled.drug.nunique()} drugs")

    with open(os.path.join(args.out_dir, "metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=2)
    log(json.dumps({k: metrics[k] for k in ("test_per_perturbation", "test_per_drug")}, indent=2))


if __name__ == "__main__":
    main()
