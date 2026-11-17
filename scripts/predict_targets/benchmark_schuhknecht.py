#!/usr/bin/env python
"""
Schuhknecht A549 benchmark: query trained target-prediction models (metabolic_target_prediction.py, --mode
metchange, cross-validated with --n-folds) with the measured metabolome changes of drugs (Schuhknecht et al. 2024,
10 uM, 48 h, intracellular), and check whether the queried drug's targets score high.

Query construction (mirrors the training signatures as far as the data allow):
  * ions -> GEM metabolites via cpd_info_wSchuhknecht.csv (KEGG / HMDB, scripts/schuhknecht_map). Only ions that map
    to one GEM metabolite (dup_after_shuhknecht_map == False) by default; --ambiguous spread also keeps the ions
    annotated with several GEM metabolites.
  * intracellular metabolomics does not resolve compartments: each ion is spread +-1/n over its n intracellular
    network nodes (every compartment except e). --restrict-nodes trained (default) keeps only nodes that are MetChange
    signature features, i.e. nodes the model has seen as sources; all keeps every intracellular node.
  * value: Schuhknecht mean log2-FC x --query-sign (+1: positive = more metabolite, the convention the MetChange
    signatures follow with --signature-sign -1). With the run's --signature-scale, each ion is divided by its spread
    across all Schuhknecht drugs, as the training signatures were across perturbations.
  * selection: as the run (config_fold*.json): top-n by |value| (--ignore-qvalue), q <= qvalue with a random cap
    (--qvalue), or q <= qvalue then the top-n by |value| (--qvalue-top-n); the Schuhknecht q-values. The drug's own
    metabolites are removed (--drug-metabolites) and queries with < min_sources ions are dropped.

Models: each fold's FNN ensemble (fold<i>/models.pt). A drug held out in fold i (it appears in fold<i>/
test_drug_ranks.tsv) is queried with fold i's model only; a drug that was never trained on (no active perturbation
in A549, so in no fold's test ranks) with the mean of the five folds' scores. Groups: heldout_cv, never_trained.

Candidates and queries: the Schuhknecht drugs matched to LINCS (drug_match_table.csv) that have known targets. A
drug's score is the mean predicted probability of its targets (as rank_prob). Metrics per queried drug d, as
percentiles (0 = best, 0.5 = chance):
  spec    rank of d's own query among all queries by the score they give d (same model): the analogue of
          spec_prob_heldout, blind to any per-drug prior
  rankc   rank of d among all candidates after standardising each candidate's scores over all queries (rank_probc)
  rank    rank of d among the candidates the model has not been trained on (raw scores, prior included)
Null: the query -> drug assignment permuted (--n-perm) on the same score matrices; p = P(null median <= observed).

Outputs (--out-dir): per_drug.tsv, summary.tsv, queries.tsv (sources per query), scores.npz (score matrices per
model), benchmark.log via stdout.

Run it on a run and on its _shuff control: the control's models read the same queries, so a benchmark signal in the
real run but not in the control comes from what the model learnt from signatures, not from the network or the prior.
"""
import argparse
import json
import os
import sys
import time
from types import SimpleNamespace

import numpy as np
import pandas as pd
import torch
from scipy.stats import rankdata

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import metabolic_target_prediction as mtp  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
log = mtp.log


def localize(path, path_map):
    """Cluster paths in config_fold*.json -> local ones: explicit --path-map prefixes first, then anything after
    'metabognn_dev/' re-rooted at this repository."""
    if not path or os.path.exists(path):
        return path
    for old, new in path_map:
        if path.startswith(old) and os.path.exists(new + path[len(old):]):
            return new + path[len(old):]
    if "metabognn_dev/" in path:
        p = os.path.join(REPO, path.split("metabognn_dev/", 1)[1])
        if os.path.exists(p):
            return p
    return path


def read_schuhknecht(xlsx, cache_dir):
    """mean log2-FC and q-values sheets -> (ions x drugs) frames indexed by the ion annotation (cached as pickle)."""
    out = {}
    for sheet, key in (("mean log2-FC", "lfc"), ("q-values", "q")):
        f = os.path.join(cache_dir, f"schuhknecht_{key}.pkl")
        if os.path.exists(f) and os.path.getmtime(f) > os.path.getmtime(xlsx):
            out[key] = pd.read_pickle(f)
            continue
        t0 = time.time()
        x = pd.read_excel(xlsx, sheet_name=sheet)
        x = x.set_index(x.columns[0]).iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
        x.index = x.index.astype(str)
        x.to_pickle(f)
        log(f"read sheet '{sheet}' ({time.time() - t0:.0f}s): {x.shape[0]} ions x {x.shape[1]} drugs")
        out[key] = x
    return out["lfc"], out["q"]


def ion_nodes(cpd_info, g, allowed, ambiguous):
    """ion annotation -> intracellular network node indices (sorted, unique)."""
    c = pd.read_csv(cpd_info, index_col=0)
    c = c[c.Schuhknecht_mapping.notna() & (c.compartment != "e") & c.id.isin(g.pos.index)]
    if ambiguous == "drop":
        c = c[c.dup_after_shuhknecht_map == False]  # noqa: E712 (column holds True / False / NaN)
    if allowed is not None:
        c = c[c.id.isin(allowed)]
    return {a: np.unique(g.pos[ids].to_numpy()) for a, ids in c.groupby("Schuhknecht_mapping").id}


def select_ions(v, q, mode, top_n, qvalue, max_sources, min_abs, rng):
    keep = v.notna() & (v.abs() >= min_abs)
    if mode == "top":
        keep &= v.abs().rank(ascending=False, method="first") <= top_n
    else:
        keep &= q.reindex(v.index) <= qvalue
        if mode == "sigtop":
            keep &= v.where(keep).abs().rank(ascending=False, method="first") <= top_n
    ions = v.index[keep].to_numpy()
    if mode == "sig" and len(ions) > max_sources:
        ions = ions[rng.choice(len(ions), max_sources, replace=False)]
    return ions


def training_sources(cfg, g, perts_table, targets, drug_mets, wanted):
    """Sources of the perturbations in `wanted` as metabolic_target_prediction.main() built them: same signature
    loading, filters, --shuffle-signatures permutation and per-perturbation rng (seed + row of the filtered table)."""
    perts, sigs, fmap = mtp.load_signatures(g, cfg, perts_table)
    perts = perts[perts.n_sig >= cfg.min_sources].reset_index(drop=True)
    if cfg.de_counts:
        n_de = pd.read_csv(cfg.de_counts, index_col=0)[cfg.de_column]
        perts["active"] = (perts.pert.map(n_de) >= cfg.min_de_genes).to_numpy()
    perts = perts[perts.drug.str.lower().isin(targets.keys())].reset_index(drop=True)
    if cfg.max_perts:
        perts = perts.sample(min(cfg.max_perts, len(perts)), random_state=cfg.seed).reset_index(drop=True)
    if cfg.shuffle_signatures:
        srng = np.random.default_rng([cfg.seed, 1])
        strata = perts.active.to_numpy() if cfg.de_counts else np.ones(len(perts), bool)
        shuffled = {}
        for st in np.unique(strata):
            src = perts.pert[strata == st].to_numpy()
            for p, q in zip(src, srng.permutation(src)):
                shuffled[p] = sigs[q]
        sigs.update(shuffled)
    row = {p: j for j, p in enumerate(perts.pert)}
    drug = dict(zip(perts.pert, perts.drug))
    return {p: mtp.select_sources(*sigs[p], fmap, drug_mets.get(drug[p].lower(), set()), cfg.max_sources,
                                  np.random.default_rng(cfg.seed + row[p])) for p in wanted if p in row}


def load_fold_models(path, n_features, device):
    st = torch.load(path, map_location="cpu", weights_only=False)
    a = SimpleNamespace(**st["args"])
    poly = mtp.PolyMinMax(n_features, poly=not a.no_poly)
    poly.load_state_dict(st["poly"])
    in_dim = int(poly.keep.sum())
    vae = None
    if st["vae"] is not None:
        vae = mtp.VAE_decrease(in_dim, a.vae_hidden, a.vae_layers, 0.0, a.latent)
        vae.load_state_dict(st["vae"])
        vae.eval()
    pre = mtp.Preprocessor(poly, vae).to(device).eval()
    ffns = []
    for s in st["ffns"]:
        m = mtp.FFN_decay(in_dim, a.hidden, 3, a.layers, a.dropout, a.decay)
        m.load_state_dict(s)
        ffns.append(m.to(device).eval())
    return pre, ffns, st["feature_names"]


def pct_rank_desc(x, axis):
    """Average rank of each value in descending order (1 = highest) / number of non-NaN values; NaN entries rank
    after all others (so they do not affect the ranks of the rest) and are returned as NaN."""
    nan = np.isnan(x)
    r = rankdata(-np.where(nan, -np.inf, x), axis=axis, method="average")
    return np.where(nan, np.nan, r / np.sum(~nan, axis=axis, keepdims=True))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, help="cross-validated run (config_fold<i>.json, fold<i>/models.pt)")
    ap.add_argument("--schuhknecht-dir", default=os.path.join(REPO, "data", "validation", "schuhknecht"))
    ap.add_argument("--out-dir", help="default: <run-dir>/schuhknecht_benchmark[_<tag>]")
    ap.add_argument("--tag", default="", help="suffix of the default out dir (variants of the query)")
    ap.add_argument("--cache-dir", help="where the Schuhknecht sheets are cached as pickles (default: the out dir); "
                                        "share it between runs to read the xlsx once")
    ap.add_argument("--threads", type=int, default=0,
                    help="torch CPU threads (default: SLURM_CPUS_PER_TASK, else the CPUs this process may use); "
                         "torch's own default counts all cores of the node and oversubscribes a Slurm allocation")
    ap.add_argument("--path-map", nargs=2, action="append", default=[], metavar=("CLUSTER_PREFIX", "LOCAL_PREFIX"),
                    help="rewrite input paths of the run config (repeatable); by default paths under "
                         "metabognn_dev/ are re-rooted at this repository when they do not exist")
    ap.add_argument("--restrict-nodes", choices=["trained", "all"], default="trained",
                    help="trained: only nodes that are MetChange signature features (sources seen in training)")
    ap.add_argument("--ambiguous", choices=["drop", "spread"], default="drop",
                    help="ions mapped to several GEM metabolites: drop, or spread over all their nodes")
    ap.add_argument("--selection", choices=["auto", "top", "sig", "sigtop"], default="auto",
                    help="auto: as the run (top: --ignore-qvalue; sig: --qvalue; sigtop: --qvalue-top-n)")
    ap.add_argument("--top-n", type=int, help="default: the run's")
    ap.add_argument("--qvalue", type=float, help="default: the run's")
    ap.add_argument("--query-sign", type=float, default=1.0)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--chunk", type=int, default=16, help="queries per feature batch")
    ap.add_argument("--max-queries", type=int, default=0, help="smoke test")
    ap.add_argument("--check-n", type=int, default=8,
                    help="reproduction check: rebuild this many held-out LINCS perturbations of the first fold from "
                         "the run's own inputs and compare their drug scores with fold<i>/test_scores_heldout.npz "
                         "(0: skip). Fails if they do not match, i.e. the inputs here differ from the training ones")
    ap.add_argument("--check-min-r", type=float, default=0.999)
    ap.add_argument("--allow-mismatch", action="store_true", help="continue when the reproduction check fails")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--device", default="auto",
                    help="auto: cuda if available, else cpu (mps is skipped: it crashes on macOS 14 with recent torch)")
    args = ap.parse_args()

    run = os.path.abspath(args.run_dir)
    folds = sorted(int(f[len("config_fold"):-5]) for f in os.listdir(run)
                   if f.startswith("config_fold") and f.endswith(".json")
                   and os.path.exists(os.path.join(run, f"fold{f[len('config_fold'):-5]}", "models.pt")))
    if not folds:
        sys.exit(f"no config_fold<i>.json with fold<i>/models.pt in {run}")
    with open(os.path.join(run, f"config_fold{folds[0]}.json")) as fh:
        cfg = SimpleNamespace(**json.load(fh))
    if cfg.mode != "metchange":
        sys.exit(f"--mode {cfg.mode}: the Schuhknecht queries are metabolites, only metchange runs can be queried")
    pm = [tuple(p) for p in args.path_map]
    for k, v in (("qvalue_top_n", False), ("signature_scale", "none"), ("signature_scale_floor", 0.05)):
        if not hasattr(cfg, k):
            setattr(cfg, k, v)
    for k in ("nodes", "edges", "baseline", "perturbations", "signatures", "qvalues", "drug_targets",
              "drug_metabolites", "de_counts", "hgnc"):
        setattr(cfg, k, localize(getattr(cfg, k), pm))
        if getattr(cfg, k) and not os.path.exists(getattr(cfg, k)):
            if k in ("qvalues", "de_counts", "hgnc"):  # only needed by the reproduction check
                log(f"warning: {k} {getattr(cfg, k)} not found")
                setattr(cfg, k, None)
            else:
                sys.exit(f"{k}: {getattr(cfg, k)} not found (use --path-map)")
    mode = args.selection if args.selection != "auto" else \
        ("top" if not cfg.qvalues or cfg.ignore_qvalue else "sigtop" if getattr(cfg, "qvalue_top_n", False) else "sig")
    top_n = args.top_n or cfg.top_n
    qvalue = args.qvalue if args.qvalue is not None else cfg.qvalue
    scale = getattr(cfg, "signature_scale", "none")
    out = args.out_dir or os.path.join(run, "schuhknecht_benchmark" + (f"_{args.tag}" if args.tag else ""))
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "config.json"), "w") as fh:
        json.dump(dict(vars(args), folds=folds, selection_used=mode, top_n_used=top_n, qvalue_used=qvalue,
                       signature_scale=scale, run_config={k: v for k, v in vars(cfg).items()}), fh, indent=2)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu") if args.device == "auto" else \
        torch.device(args.device)
    if device.type == "cpu":
        n_thr = args.threads or int(os.environ.get("SLURM_CPUS_PER_TASK", 0))
        if not n_thr:
            n_thr = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()
        torch.set_num_threads(max(1, n_thr))
    log(f"torch threads: {torch.get_num_threads()}")
    log(f"run {run}: folds {folds}; selection {mode} (top_n {top_n}, q {qvalue}), scale {scale}, nodes "
        f"{args.restrict_nodes}, ambiguous ions {args.ambiguous}; device {device}")

    # ---- network, targets (as in training)
    g = mtp.Graph(cfg.nodes, cfg.edges)
    if cfg.edge_weights == "uniform":
        g.edges["weight"] = 1.0
    act = mtp.node_activity(g, cfg.baseline, cfg.activity_floor)
    perts_table = pd.read_csv(cfg.perturbations, sep="\t", dtype=str, keep_default_na=False)
    targets = mtp.load_drug_targets(cfg.drug_targets, g, perts_table)
    drug_mets = mtp.load_drug_metabolites(cfg.drug_metabolites, g)

    # ---- drugs: Schuhknecht <-> LINCS name, with known targets
    sd = args.schuhknecht_dir
    dm = pd.read_csv(os.path.join(sd, "drug_match_table.csv"), index_col=0)
    dm["drug"] = dm.compound.str.lower()
    dm = dm.drop_duplicates("drug")
    n_all = len(dm)
    dm = dm[dm.drug.isin(targets.keys())].reset_index(drop=True)
    log(f"Schuhknecht drugs matched to LINCS: {n_all}; with known targets: {len(dm)}")

    heldout = {}
    for i in folds:
        f = os.path.join(run, f"fold{i}", "test_drug_ranks.tsv")
        for d in pd.read_csv(f, sep="\t", usecols=["drug"]).drug.str.lower().unique():
            heldout[d] = i
    dm["fold"] = dm.drug.map(heldout)
    dm["group"] = np.where(dm.fold.notna(), "heldout_cv", "never_trained")
    if len(folds) < cfg.n_folds:
        log(f"warning: only folds {folds} of {cfg.n_folds} available; drugs held out in a missing fold count as "
            f"never_trained although they were trained on")

    # ---- ions -> nodes, query values
    allowed = None
    if args.restrict_nodes == "trained":
        feats = pd.read_csv(cfg.signatures, usecols=[0]).iloc[:, 0].astype(str)
        allowed = set(feats[feats.isin(g.pos.index)])
        log(f"trained source nodes: {len(allowed):,} MetChange features that are network nodes")
    imap = ion_nodes(os.path.join(sd, "cpd_info_wSchuhknecht.csv"), g, allowed, args.ambiguous)
    cache = args.cache_dir or out
    os.makedirs(cache, exist_ok=True)
    lfc, qv = read_schuhknecht(os.path.join(sd, "41587_2024_2524_MOESM2_ESM.xlsx"), cache)
    ions = [a for a in imap if a in lfc.index]
    log(f"ions with intracellular nodes: {len(ions)} ({sum(len(imap[a]) for a in ions)} nodes, "
        f"{len(np.unique(np.concatenate([imap[a] for a in ions])))} distinct)")
    V = lfc.loc[ions] * args.query_sign
    if scale != "none":  # spread of each ion across all Schuhknecht drugs, as the run scaled across perturbations
        V = mtp.scale_signatures(V, scale, cfg.signature_scale_floor)
    Q = qv.loc[ions]

    # ---- per drug: sources
    qrows, sources = [], []
    for j, r in dm.iterrows():
        if args.max_queries and len(sources) >= args.max_queries:
            break
        sel = select_ions(V[r.schuhknecht_name], Q[r.schuhknecht_name], mode, top_n, qvalue, cfg.max_sources,
                          cfg.min_abs_score, np.random.default_rng(args.seed + j))
        drop = drug_mets.get(r.drug, set())
        idx, w, used = [], [], []
        for a in sel:
            nodes = np.array([x for x in imap[a] if x not in drop])
            if len(nodes):
                idx.append(nodes)
                w.append(np.full(len(nodes), np.sign(V.at[a, r.schuhknecht_name]) / len(imap[a])))
                used.append(a)
        qrows.append(dict(drug=r.drug, schuhknecht_name=r.schuhknecht_name, n_ions=len(used),
                          n_up=int(sum(V.at[a, r.schuhknecht_name] > 0 for a in used)),
                          n_nodes=int(sum(len(x) for x in idx)), ions=" | ".join(used)))
        sources.append((np.concatenate(idx), np.concatenate(w)) if idx else (np.empty(0, int), np.empty(0)))
    qdf = pd.DataFrame(qrows)
    qdf.to_csv(os.path.join(out, "queries.tsv"), sep="\t", index=False)
    valid = (qdf.n_ions >= cfg.min_sources).to_numpy()
    log(f"queries: {len(qdf)}; with >= {cfg.min_sources} ions: {valid.sum()} (median {int(qdf.n_ions.median())} ions)")
    qdf = qdf[valid].reset_index(drop=True)
    sources = [s for s, v in zip(sources, valid) if v]
    dm = dm.set_index("drug").loc[qdf.drug].reset_index()

    # ---- candidates = queried drugs (same universe for every metric)
    cands = list(qdf.drug)
    ranker = mtp.DrugRanker({d: targets[d] for d in cands}, len(g.prot))
    cand_prot = np.unique(np.concatenate([list(targets[d]) for d in cands]))
    ccol = np.array([ranker.col[d] for d in cands])  # candidate order of DrugRanker -> query order
    log(f"candidates: {len(cands)} drugs, {len(cand_prot)} target proteins scored")

    # ---- reproduction check sources (LINCS perturbations held out in the first fold)
    chk = {}
    if args.check_n and cfg.shuffle_signatures and cfg.min_de_genes and not cfg.de_counts:
        log("warning: reproduction check skipped (shuffled run: the permutation needs --de-counts, not found)")
    elif args.check_n and cfg.qvalues is None and not cfg.ignore_qvalue:
        log("warning: reproduction check skipped (the run selected by q-value: --qvalues not found)")
    elif args.check_n:
        z = np.load(os.path.join(run, f"fold{folds[0]}", "test_scores_heldout.npz"), allow_pickle=True)
        wanted = list(z["pert"][:args.check_n])
        chk = training_sources(cfg, g, perts_table, targets, drug_mets, wanted)
        chk = {p: chk[p] for p in wanted if p in chk}
        log(f"reproduction check: {len(chk)} of {len(wanted)} perturbations rebuilt")

    # ---- features and scores (every fold model scores every query)
    fb = mtp.FeatureBuilder(g, act, mtp.membership_matrix(g),
                            np.unique(np.concatenate([s[0] for s in sources] + [s[0] for s in chk.values()])),
                            SimpleNamespace(hops=cfg.hops, alpha=cfg.alpha, no_static=cfg.no_static,
                                            no_globals=cfg.no_globals))
    models = {}
    for i in folds:
        pre, ffns, names = load_fold_models(os.path.join(run, f"fold{i}", "models.pt"), len(fb.names), device)
        if list(names) != list(fb.names):
            sys.exit(f"fold {i}: feature names of models.pt differ from the rebuilt features")
        models[i] = (pre, ffns)
    if chk:
        pre, ffns = models[folds[0]]
        Xc = fb.build(list(chk.values()), np.zeros((len(chk), 3), np.float32))
        Pc = mtp.predict_proba(ffns, pre, Xc, device).mean(0).reshape(len(chk), len(g.prot), 3)
        rk = mtp.DrugRanker({d: targets[d] for d in z["drug"]}, len(g.prot))
        prow = {p: i for i, p in enumerate(z["pert"])}
        rs, md = [], []
        for k, p in enumerate(chk):
            mine, saved = rk.prob_scores(Pc[k]), z["prob"][prow[p]]
            rs.append(np.corrcoef(mine, saved)[0, 1])
            md.append(np.abs(mine - saved).max())
        ok = min(rs) >= args.check_min_r
        msg = (f"reproduction check (fold {folds[0]}, {len(chk)} held-out LINCS perturbations vs "
               f"test_scores_heldout.npz): min r {min(rs):.6f}, max |diff| {max(md):.2e} -> "
               f"{'OK' if ok else 'MISMATCH'}")
        log(msg)
        with open(os.path.join(out, "reproduction_check.txt"), "w") as fh:
            fh.write(msg + "\n" + "\n".join(f"{p}\tr={r:.6f}\tmax_abs_diff={d:.3e}" for p, r, d in zip(chk, rs, md))
                     + "\n")
        if not ok and not args.allow_mismatch:
            sys.exit("the rebuilt features do not reproduce the training scores: the network, baseline, labels or "
                     "drug tables here differ from the ones the run was trained on (run on the cluster, or sync "
                     "them); --allow-mismatch to continue anyway")
    gl = np.stack([np.zeros(len(qdf)), np.full(len(qdf), np.log1p(48.0)), np.log1p(qdf.n_ions.to_numpy())],
                  1).astype(np.float32)  # rel. dose 0 (one dose), 48 h, n sources; unused with --no-globals
    n_prot, nc = len(g.prot), len(cand_prot)
    S = {i: np.empty((len(qdf), len(cands)), np.float32) for i in folds}
    t0 = time.time()
    for s in range(0, len(qdf), args.chunk):
        sl = slice(s, min(s + args.chunk, len(qdf)))
        X = fb.build(sources[sl], gl[sl])
        P = sl.stop - sl.start
        X = X.reshape(P, n_prot, -1)[:, cand_prot].reshape(P * nc, -1)
        for i, (pre, ffns) in models.items():
            prob = mtp.predict_proba(ffns, pre, X, device).mean(0).reshape(P, nc, 3)
            full = np.zeros((n_prot, 3))
            for k in range(P):
                full[cand_prot] = prob[k]
                S[i][sl.start + k] = ranker.prob_scores(full)[ccol]
        log(f"scored {sl.stop}/{len(qdf)} queries ({time.time() - t0:.0f}s)")
    S["mean"] = np.mean([S[i] for i in folds], 0)
    np.savez_compressed(os.path.join(out, "scores.npz"), queries=np.array(cands), candidates=np.array(cands),
                        **{f"fold{k}" if k != "mean" else "mean": v for k, v in S.items()})

    # ---- metrics: percentile matrices per model (rows: queries, columns: candidates, same order)
    trained_by = {i: {d for d, f in heldout.items() if f != i} for i in folds}
    R, RC, RR = {}, {}, {}
    for k, M in S.items():
        R[k] = pct_rank_desc(M, 0)  # query's rank in the candidate's column (spec)
        Z = (M - M.mean(0)) / (M.std(0) + 1e-12)
        RC[k] = pct_rank_desc(Z, 1)  # candidate's rank in the query's row after column standardisation (rankc)
        unseen = np.array([d not in heldout for d in cands]) if k == "mean" else \
            np.array([d not in trained_by[k] for d in cands])
        RR[k] = pct_rank_desc(np.where(unseen[None, :], M, np.nan), 1)  # among untrained candidates (rank)
    key = [int(f) if f == f else "mean" for f in dm.fold]
    cidx = np.arange(len(cands))

    def metric_values(assign):
        """assign[d] = query row given to drug d -> per drug spec / rankc / rank."""
        out_ = {m: np.empty(len(cands)) for m in ("spec", "rankc", "rank")}
        for k in S:
            sel = np.array([x == k for x in key])
            if sel.any():
                for m, A in (("spec", R), ("rankc", RC), ("rank", RR)):
                    out_[m][sel] = A[k][assign[sel], cidx[sel]]
        return out_

    obs = metric_values(cidx)
    res = dm[["drug", "schuhknecht_name", "group", "fold", "validated_metabolic"]].copy()
    res["n_ions"] = qdf.n_ions.to_numpy()
    res["n_candidates_rank"] = [int(np.sum(~np.isnan(RR[k][0]))) for k in key]
    for m, v in obs.items():
        res[f"{m}_pct"] = v
    res.to_csv(os.path.join(out, "per_drug.tsv"), sep="\t", index=False)

    groups = {"all": np.ones(len(res), bool), "heldout_cv": (res.group == "heldout_cv").to_numpy(),
              "never_trained": (res.group == "never_trained").to_numpy(),
              "validated_metabolic": res.validated_metabolic.to_numpy(bool)}
    rng = np.random.default_rng(args.seed)
    null = {m: {gname: [] for gname in groups} for m in obs}
    for _ in range(args.n_perm):
        pv = metric_values(rng.permutation(cidx))
        for m in obs:
            for gname, gm in groups.items():
                if gm.any():
                    null[m][gname].append((np.nanmedian(pv[m][gm]), np.nanmean(pv[m][gm] <= 0.1)))
    rows = []
    for m in obs:
        for gname, gm in groups.items():
            x = obs[m][gm]
            x = x[~np.isnan(x)]
            if not len(x):
                continue
            nm = np.array(null[m][gname])
            med, top10 = float(np.median(x)), float(np.mean(x <= 0.1))
            rows.append(dict(metric=m, group=gname, n=len(x), median_pct=med, top10_pct=top10,
                             top1_pct=float(np.mean(x <= 0.01)), frac_better_than_half=float(np.mean(x < 0.5)),
                             null_median=float(nm[:, 0].mean()), null_top10=float(nm[:, 1].mean()),
                             p_median=float((1 + np.sum(nm[:, 0] <= med)) / (1 + len(nm))),
                             p_top10=float((1 + np.sum(nm[:, 1] >= top10)) / (1 + len(nm)))))
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(out, "summary.tsv"), sep="\t", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        log("summary (percentiles, 0 = best, 0.5 = chance):\n" + summ.round(4).to_string(index=False))
    log(f"written to {out}")


if __name__ == "__main__":
    main()
