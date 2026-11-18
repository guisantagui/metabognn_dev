#!/usr/bin/env python
"""
Signature-similarity baseline, scored exactly like a trained run.

How much of a run's gap over its shuffled control does plain signature similarity already give, with no network?
For every held-out perturbation of the run (its own folds, test perturbations and held-out candidate drugs), the k
training drugs whose perturbations are most similar to the query vote their targets: the protein scores are the
similarity-weighted mean of the neighbours' target vectors, built like the model's prior (activation / inhibition
column, 0.5 to both for an unknown direction). That (proteins x 3) matrix then goes through the model's own
DrugRanker.prob_scores, signature_specificity and centered_ranks, so rank_prob / rank_probc / spec_prob mean exactly
what they mean for the model. The output is laid out like a run (fold<i>/test_drug_ranks.tsv,
test_drug_ranks_per_drug.tsv, test_scores_heldout.npz, cv_drug_ranks*.tsv), so compare_runs.py and
plot_performance_panels.py (--extra) read it unchanged, including the dose/time-stratified spec.

One call does a run AND its _shuff control (--sides), every --representation and every --k, reading the label
files once: the two sides share their inputs and differ only by the --shuffle-signatures permutation, which is
replayed here (same seed, same activity strata), so baseline-real vs baseline-shuffled is the same contrast as
model-real vs model-shuffled, on the same queries and candidates. The similarity matrix of a fold is computed once
per representation and reused for every k.

--representation (one or both)
  sources  exactly what the model reads: the selected features (top-n / q-filtered as the run), their signs, mapped
           onto network nodes by select_sources (same per-perturbation rng, the drug's own metabolites dropped).
           Cosine similarity. Answers: does the network add anything to the same input?
  profile  the full signature: every mapped feature, its continuous value (signed and scaled as the run), Spearman.
           Answers: is there more in the data than the model's input keeps?

References: the training drugs of the fold (not validation, not held out) and their active perturbations, which is
what the model trains on. --unit drug: a drug's similarity to the query is that of its most similar perturbation,
and the k most similar drugs vote with weight max(similarity, 0) (uniform if none is positive). --unit pert: the k
most similar perturbations vote (a drug with many doses / times can then take several votes).

Never the held-out drug's own perturbations: the drug is held out, as for the model. Ranking held-out drugs by
similarity to their own other doses / times would measure the same-drug reproducibility ceiling
(check_signature_target_similarity.py: AUROC 0.85 for expr), not a baseline.

Inputs are the run's own (config_fold<i>.json); paths under metabognn_dev/ that do not exist are re-rooted at this
repository. Before scoring, each fold of each side is checked against that run: the held-out drugs rebuilt here
must equal the candidate list in fold<i>/test_scores_heldout.npz, and every scored test perturbation must be found;
otherwise the inputs differ from the training ones and the script stops (--allow-mismatch to continue, for smoke
tests only). The two sides' configs must agree on everything but the shuffle.

Output: <out-root>/<run>__knn_<representation>_k<k>[_pert] for the run and ..._shuff for its control (out-root
defaults to the run's parent directory), i.e. pairs the other tools treat as a run and its control:
    python compare_runs.py --run real=<dir> --run shuffled=<dir>_shuff --control shuffled --out-dir <dir>/compare
    python plot_performance_panels.py ... --extra "mc=<run>__knn_sources_k25:MetChange kNN, same input"
"""
import argparse
import json
import os
import sys
import time
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import metabolic_target_prediction as mtp  # noqa: E402

T0 = time.time()
#: config keys the real run and its _shuff control must share (everything that decides inputs, folds and selection)
SHARED_KEYS = ("mode", "nodes", "edges", "perturbations", "signatures", "qvalues", "drug_targets", "drug_metabolites",
               "de_counts", "de_column", "min_de_genes", "hgnc", "signature_sign", "signature_scale",
               "signature_scale_floor", "ignore_qvalue", "qvalue", "qvalue_top_n", "top_n", "min_abs_score",
               "min_sources", "max_sources", "max_perts", "seed", "n_folds", "val_frac")


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def localize(path):
    """Cluster paths: re-root what follows metabognn_dev/ at this repository when the original does not exist."""
    if not path or os.path.exists(path):
        return path
    key = "metabognn_dev/"
    if key in path:
        cand = os.path.join(REPO, path.split(key, 1)[1])
        if os.path.exists(cand):
            return cand
    return path


def run_folds(run):
    """Folds with a config and the score files a run writes per fold."""
    out = []
    for f in sorted(os.listdir(run)):
        if f.startswith("config_fold") and f.endswith(".json"):
            i = int(f[len("config_fold"):-5])
            d = os.path.join(run, f"fold{i}")
            if all(os.path.exists(os.path.join(d, x)) for x in ("test_drug_ranks.tsv", "test_scores_heldout.npz")):
                out.append(i)
    return sorted(out)


def read_config(run, allow_mismatch):
    folds = run_folds(run)
    if not folds:
        sys.exit(f"{run}: no fold with config_fold<i>.json, test_drug_ranks.tsv and test_scores_heldout.npz")
    with open(os.path.join(run, f"config_fold{folds[0]}.json")) as fh:
        cfg = SimpleNamespace(**json.load(fh))
    for k, v in (("qvalue_top_n", False), ("signature_scale", "none"), ("signature_scale_floor", 0.05),
                 ("max_perts", 0), ("shuffle_signatures", False)):
        if not hasattr(cfg, k):
            setattr(cfg, k, v)
    for k in ("nodes", "edges", "perturbations", "signatures", "qvalues", "drug_targets", "drug_metabolites",
              "de_counts", "hgnc"):
        setattr(cfg, k, localize(getattr(cfg, k, None)))
        if getattr(cfg, k) and not os.path.exists(getattr(cfg, k)):
            if allow_mismatch and k == "de_counts":  # smoke test only: every perturbation counts as active
                log(f"warning: {k} {getattr(cfg, k)} not found, every perturbation treated as active")
                setattr(cfg, k, None)
            else:
                sys.exit(f"{k}: {getattr(cfg, k)} not found")
    return cfg, folds


def filtered_perturbations(cfg, g, perts_table, targets):
    """metabolic_target_prediction.main()'s perturbation table up to the shuffle: signatures, min_sources, activity,
    drugs with targets, max_perts. The shuffle does not change this table, so both sides share it."""
    perts, sigs, fmap = mtp.load_signatures(g, cfg, perts_table)
    perts = perts[perts.n_sig >= cfg.min_sources].reset_index(drop=True)
    if cfg.de_counts:
        n_de = pd.read_csv(cfg.de_counts, index_col=0)[cfg.de_column]
        perts["n_de_genes"] = perts.pert.map(n_de)
        perts["active"] = (perts.n_de_genes >= cfg.min_de_genes).to_numpy()
    else:
        perts["active"] = True
    perts = perts[perts.drug.str.lower().isin(targets.keys())].reset_index(drop=True)
    if cfg.max_perts:
        perts = perts.sample(min(cfg.max_perts, len(perts)), random_state=cfg.seed).reset_index(drop=True)
    return perts, sigs, fmap


def replay_shuffle(cfg, perts, sigs):
    """The --shuffle-signatures permutation of main(): perm[p] = the perturbation whose signature p carries, and the
    signatures as the side sees them (unchanged without shuffling)."""
    perm = {p: p for p in perts.pert}
    if not cfg.shuffle_signatures:
        return perm, sigs
    srng = np.random.default_rng([cfg.seed, 1])
    strata = perts.active.to_numpy() if cfg.de_counts else np.ones(len(perts), bool)
    out = dict(sigs)
    for st in np.unique(strata):
        src = perts.pert[strata == st].to_numpy()
        for p, q in zip(src, srng.permutation(src)):
            perm[p] = q
            out[p] = sigs[q]
    log(f"shuffle-signatures replayed: {sum(p != q for p, q in perm.items()):,} signatures permuted")
    return perm, out


def source_matrix(perts, sigs, fmap, drug_mets, cfg, n_nodes):
    """Perturbations x nodes, the model's sources (select_sources with the per-row rng of training), L2-normalised."""
    rows, cols, vals = [], [], []
    for j, (p, d) in enumerate(zip(perts.pert, perts.drug)):
        idx, w = mtp.select_sources(*sigs[p], fmap, drug_mets.get(d.lower(), set()), cfg.max_sources,
                                    np.random.default_rng(cfg.seed + j))
        rows.append(np.full(len(idx), j))
        cols.append(idx)
        vals.append(w)
    M = sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(len(perts), n_nodes))  # duplicate (row, node) entries are summed
    norm = np.sqrt(np.asarray(M.multiply(M).sum(1)).ravel())
    return sparse.diags(1 / np.where(norm > 0, norm, 1)) @ M


def profile_matrix(perts, fmap, cfg):
    """Perturbations x mapped features, each perturbation's own full signed / scaled signature, rank-transformed,
    centred and L2-normalised so that the dot product is Spearman. A shuffled side only reorders the rows."""
    header = pd.read_csv(cfg.signatures, nrows=0).columns  # the label file has every A549 perturbation: read
    want = set(perts.pert)                                  # only the columns needed here
    val = pd.read_csv(cfg.signatures, index_col=0, usecols=[header[0]] + [c for c in header[1:] if c in want])
    val.index = val.index.astype(str)
    val = val * cfg.signature_sign
    if cfg.signature_scale != "none":
        val = mtp.scale_signatures(val, cfg.signature_scale, cfg.signature_scale_floor)
    if cfg.mode in ("tf", "expr"):
        hgnc = mtp.HGNC(cfg.hgnc)
        sym = val.index.map(hgnc)
        keep = (sym != "") & ~pd.Index(sym).duplicated()
        val = val[keep]
        val.index = sym[keep]
    val = val.loc[[f for f in val.index if f in fmap]]
    X = np.nan_to_num(val[list(perts.pert)].to_numpy(np.float32).T, nan=0.0)  # perturbations x features
    X = rankdata(X, axis=1).astype(np.float32)
    X -= X.mean(1, keepdims=True)
    X /= np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
    log(f"profiles: {X.shape[0]:,} perturbations x {X.shape[1]:,} mapped features (Spearman)")
    return X


def target_tokens(drugs, targets, n_prot):
    """Drugs x (2 * proteins) target tokens, as the model's prior: class c -> token 2 * gene + c - 1 (1: activation,
    2: inhibition), unknown direction (3) -> 0.5 on both."""
    r, c, v = [], [], []
    for i, d in enumerate(drugs):
        for gene, cls in targets[d].items():
            for t in ([2 * gene, 2 * gene + 1] if cls == 3 else [2 * gene + cls - 1]):
                r.append(i)
                c.append(t)
                v.append(0.5 if cls == 3 else 1.0)
    return sparse.csr_matrix((v, (r, c)), shape=(len(drugs), 2 * n_prot))


def unit_similarity(sim, ref_unit, n_units):
    """queries x units: each unit's similarity is its maximum over its reference perturbations (-inf: none)."""
    U = np.full((sim.shape[0], n_units), -np.inf, np.float32)
    np.maximum.at(U.T, ref_unit, sim.T)
    return U


def knn_votes(U, k):
    """The k most similar units vote with weight max(sim, 0), normalised (uniform over the k if none is positive).
    Returns queries x units, sparse."""
    k = min(k, U.shape[1])
    top = np.argpartition(-U, k - 1, axis=1)[:, :k]
    w = np.take_along_axis(U, top, axis=1)
    w = np.where(np.isfinite(w), np.maximum(w, 0), 0)
    s = w.sum(1, keepdims=True)
    w = np.where(s > 0, w / np.where(s > 0, s, 1), 1.0 / k)
    rows = np.repeat(np.arange(U.shape[0]), k)
    return sparse.csr_matrix((w.ravel(), (rows, top.ravel())), shape=U.shape)


def rank_tables(S, test, rk):
    """Per perturbation and per drug: the model's rank_prob / spec_prob / rank_probc columns from scores S."""
    keys = test.drug.str.lower().to_numpy()
    groups = test.active.astype(str).str.upper().eq("TRUE").to_numpy() if "active" in test else \
        np.ones(len(test), bool)
    pp = test[[c for c in ("pert", "drug", "n_true_targets", "n_de_genes", "active") if c in test]].copy()
    pp["n_candidates_heldout"] = len(rk.drugs)
    pp["rank_prob_heldout"] = [rk.rank_of(S[j], keys[j]) for j in range(len(test))]
    pp["spec_prob_heldout"], pp["n_spec_heldout"] = mtp.signature_specificity(S, keys, rk.col, groups)
    pp["rank_probc_heldout"] = mtp.centered_ranks(S, keys, rk.col, groups)
    # per drug: scores averaged over its (active) test perturbations, against the other test drugs
    drows, dS = [], []
    for subset in ("all", "active"):
        sel_all = np.ones(len(test), bool) if subset == "all" else groups
        for d in pd.unique(keys):
            sel = (keys == d) & sel_all
            if sel.any():
                drows.append(dict(drug=d, subset=subset, n_perturbations=int(sel.sum()),
                                  n_candidates_heldout=len(rk.drugs)))
                dS.append(S[sel].mean(0))
    pdr, dS = pd.DataFrame(drows), np.asarray(dS)
    pdr["rank_prob_heldout"] = [rk.rank_of(dS[j], d) for j, d in enumerate(pdr.drug)]
    for subset in ("all", "active"):
        m = (pdr.subset == subset).to_numpy()
        if m.any():
            k_, grp = pdr.drug[m].to_numpy(), np.ones(m.sum(), bool)
            sp, n = mtp.signature_specificity(dS[m], k_, rk.col, grp)
            pdr.loc[m, "spec_prob_heldout"], pdr.loc[m, "n_spec_heldout"] = sp, n
            pdr.loc[m, "rank_probc_heldout"] = mtp.centered_ranks(dS[m], k_, rk.col, grp)
    return pp, pdr


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, help="a cross-validated run (its _shuff control is found next to it)")
    ap.add_argument("--sides", nargs="+", choices=("real", "shuff"), default=["real", "shuff"])
    ap.add_argument("--representation", nargs="+", choices=("sources", "profile"), default=["sources", "profile"])
    ap.add_argument("--k", nargs="+", type=int, default=[25], help="numbers of voting neighbours")
    ap.add_argument("--unit", choices=("drug", "pert"), default="drug", help="what a neighbour is")
    ap.add_argument("--ref-include-val", action="store_true",
                    help="also use the validation drugs as references (the model only early-stops on them)")
    ap.add_argument("--out-root", help="where the <run>__knn_* directories go (default: the run's parent)")
    ap.add_argument("--allow-mismatch", action="store_true", help="continue when the inputs do not match the run")
    ap.add_argument("--max-test", type=int, default=0, help="smoke test: score only this many perturbations per fold")
    args = ap.parse_args()

    run = os.path.abspath(args.run_dir.rstrip("/"))
    base = run[:-len("_shuff")] if run.endswith("_shuff") else run
    side_dir = {"real": base, "shuff": base + "_shuff"}
    sides = []
    for s in args.sides:
        if os.path.isdir(side_dir[s]):
            sides.append(s)
        else:
            log(f"WARNING: {side_dir[s]} not found, side skipped")
    if not sides:
        sys.exit("no side to score")
    conf = {s: read_config(side_dir[s], args.allow_mismatch) for s in sides}
    cfg = conf[sides[0]][0]
    for s in sides[1:]:
        diff = [k for k in SHARED_KEYS if getattr(conf[s][0], k, None) != getattr(cfg, k, None)]
        if diff:
            sys.exit(f"{side_dir[s]} and {side_dir[sides[0]]} differ in {diff}: not a run and its control")
    out_root = os.path.abspath(args.out_root or os.path.dirname(base))
    name = os.path.basename(base)
    suffix = "_pert" if args.unit == "pert" else ""
    log(f"run {name}: sides {sides}, mode {cfg.mode}; representations {args.representation}, k {args.k} "
        f"{args.unit}s -> {out_root}/{name}__knn_<rep>_k<k>{suffix}[_shuff]")

    # ---- inputs, read once for every side, representation and k
    g = mtp.Graph(cfg.nodes, cfg.edges)
    perts_table = pd.read_csv(cfg.perturbations, sep="\t", dtype=str, keep_default_na=False)
    targets = mtp.load_drug_targets(cfg.drug_targets, g, perts_table)
    drug_mets = mtp.load_drug_metabolites(cfg.drug_metabolites, g) if cfg.drug_metabolites else {}
    perts, sigs0, fmap = filtered_perturbations(cfg, g, perts_table, targets)
    splits = mtp.drug_splits(perts.drug.unique(), cfg)
    row = {p: j for j, p in enumerate(perts.pert)}
    n_prot = len(g.prot)
    log(f"{len(perts):,} perturbations of {perts.drug.nunique():,} drugs with targets "
        f"({int(perts.active.sum()):,} active); {len(splits)} splits")
    profile0 = profile_matrix(perts, fmap, cfg) if "profile" in args.representation else None

    for side in sides:
        scfg, folds = conf[side]
        perm, sigs = replay_shuffle(scfg, perts, sigs0)
        rdir = side_dir[side]
        X = {}
        if "sources" in args.representation:
            X["sources"] = source_matrix(perts, sigs, fmap, drug_mets, scfg, g.N)
        if profile0 is not None:
            X["profile"] = profile0[[row[perm[p]] for p in perts.pert]]  # each row carries its permuted signature
        acc = {(rep, k): ([], []) for rep in X for k in args.k}
        check = {}
        for i in folds:
            tr, va, te = splits[i]
            fd = os.path.join(rdir, f"fold{i}")
            z = np.load(os.path.join(fd, "test_scores_heldout.npz"), allow_pickle=True)
            test = pd.read_csv(os.path.join(fd, "test_drug_ranks.tsv"), sep="\t")
            cands = [str(d) for d in z["drug"]]  # the run's own held-out candidates: authoritative
            rebuilt = sorted({d.lower() for d in te})
            same = rebuilt == cands
            found = test.pert.isin(row)
            check[i] = dict(heldout_drugs_match=bool(same), n_candidates=len(cands), n_test=int(len(test)),
                            n_test_found=int(found.sum()))
            log(f"{side} fold {i}: held-out drugs rebuilt here {'match' if same else 'DO NOT match'} the run's "
                f"candidates ({len(rebuilt)} vs {len(cands)}); {int(found.sum())}/{len(test)} test perturbations "
                f"found")
            if not (same and found.all()) and not args.allow_mismatch:
                sys.exit("inputs differ from the run's (labels, network, targets?): run where the run was trained, "
                         "or --allow-mismatch for a smoke test")
            test = test[found].reset_index(drop=True)
            if args.max_test:
                test = test.head(args.max_test)
            rk = mtp.DrugRanker({d: targets[d] for d in cands}, n_prot)
            # references: training drugs (optionally + validation), active perturbations; every held-out candidate
            # excluded explicitly, so not even a mismatched smoke run can leak
            ref_drugs = sorted({d.lower() for d in (tr | va if args.ref_include_val else tr)} - set(cands))
            ref = perts[perts.drug.str.lower().isin(ref_drugs) & perts.active.astype(bool)]
            ref_idx = ref.pert.map(row).to_numpy()
            if args.unit == "drug":
                pos = {d: j for j, d in enumerate(ref_drugs)}
                ref_unit, unit_drug = ref.drug.str.lower().map(pos).to_numpy(), ref_drugs
            else:
                ref_unit, unit_drug = np.arange(len(ref)), list(ref.drug.str.lower())
            T = target_tokens(unit_drug, targets, n_prot)
            q_idx = test.pert.map(row).to_numpy()
            for rep, M in X.items():
                sim = M[q_idx] @ M[ref_idx].T
                sim = np.asarray(sim.todense() if sparse.issparse(sim) else sim, np.float32)
                U = unit_similarity(sim, ref_unit, len(unit_drug))  # once per fold and representation
                for k in args.k:
                    TOK = np.asarray((knn_votes(U, k) @ T).todense())  # queries x (2 * proteins)
                    S = np.empty((len(test), len(rk.drugs)), np.float32)
                    for j in range(len(test)):
                        prob = np.zeros((n_prot, 3))
                        prob[:, 1], prob[:, 2] = TOK[j, 0::2], TOK[j, 1::2]
                        S[j] = rk.prob_scores(prob)
                    pp, pdr = rank_tables(S, test, rk)
                    out = os.path.join(out_root, f"{name}__knn_{rep}_k{k}{suffix}" +
                                       ("_shuff" if side == "shuff" else ""))
                    od = os.path.join(out, f"fold{i}")
                    os.makedirs(od, exist_ok=True)
                    pp.to_csv(os.path.join(od, "test_drug_ranks.tsv"), sep="\t", index=False)
                    pdr.to_csv(os.path.join(od, "test_drug_ranks_per_drug.tsv"), sep="\t", index=False)
                    np.savez_compressed(os.path.join(od, "test_scores_heldout.npz"), pert=pp.pert.to_numpy(),
                                        drug=np.array(rk.drugs), prob=S)
                    m = dict(n_test=len(pp), n_references=int(len(ref)), n_units=len(unit_drug),
                             rank_prob_heldout=float(np.median(pp.rank_prob_heldout / pp.n_candidates_heldout)),
                             rank_probc_heldout=float(np.median(pp.rank_probc_heldout / pp.n_candidates_heldout)),
                             spec_prob_heldout=float(np.nanmedian(pp.spec_prob_heldout / pp.n_spec_heldout)))
                    with open(os.path.join(od, "metrics.json"), "w") as fh:
                        json.dump(m, fh, indent=2)
                    log(f"{side} fold {i} {rep} k {k}: {len(pp)} queries, {len(ref):,} references of "
                        f"{len(ref_drugs)} drugs; median percentile rank_prob {m['rank_prob_heldout']:.3f}, "
                        f"rank_probc {m['rank_probc_heldout']:.3f}, spec {m['spec_prob_heldout']:.3f}")
                    acc[(rep, k)][0].append(pp.assign(fold=i))
                    acc[(rep, k)][1].append(pdr.assign(fold=i))
        for (rep, k), (pps, pds) in acc.items():
            out = os.path.join(out_root, f"{name}__knn_{rep}_k{k}{suffix}" + ("_shuff" if side == "shuff" else ""))
            pd.concat(pps, ignore_index=True).to_csv(os.path.join(out, "cv_drug_ranks.tsv"), sep="\t", index=False)
            pd.concat(pds, ignore_index=True).to_csv(os.path.join(out, "cv_drug_ranks_per_drug.tsv"), sep="\t",
                                                     index=False)
            with open(os.path.join(out, "config.json"), "w") as fh:
                json.dump(dict(vars(args), run=rdir, side=side, representation_used=rep, k_used=k, folds=folds,
                               input_check=check, run_config={kk: v for kk, v in vars(scfg).items()}), fh, indent=2)
            log(f"saved {out}")


if __name__ == "__main__":
    main()
