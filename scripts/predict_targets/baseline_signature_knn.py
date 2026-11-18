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

Run it on a run AND on its _shuff control. The control's config carries --shuffle-signatures, and the same
permutation (seed, activity strata) is replayed here, so baseline-real vs baseline-shuffled is the same contrast as
model-real vs model-shuffled, on the same queries and candidates.

--representation
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
repository. Before scoring, each fold is checked against the run: the held-out drugs rebuilt here must equal the
candidate list in fold<i>/test_scores_heldout.npz, and every scored test perturbation must be found; otherwise the
inputs differ from the training ones and the script stops (--allow-mismatch to continue, for smoke tests only).

Output (default): <runs_dir>/<run>__knn_<representation>_k<k>[_pert] for a run and ..._shuff for its control, i.e.
a pair that the other tools treat as a run and its control:
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


def filtered_perturbations(cfg, g, perts_table, targets):
    """metabolic_target_prediction.main()'s perturbation table up to the split: signatures, min_sources, activity,
    drugs with targets, max_perts and the --shuffle-signatures permutation. Returns (perts, sigs, fmap, perm) with
    perm[p] = the perturbation whose signature p carries (p itself without shuffling)."""
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
    perm = {p: p for p in perts.pert}
    if cfg.shuffle_signatures:
        srng = np.random.default_rng([cfg.seed, 1])
        strata = perts.active.to_numpy() if cfg.de_counts else np.ones(len(perts), bool)
        orig = dict(sigs)
        for st in np.unique(strata):
            src = perts.pert[strata == st].to_numpy()
            for p, q in zip(src, srng.permutation(src)):
                perm[p] = q
                sigs[p] = orig[q]
        log(f"shuffle-signatures replayed: {sum(p != q for p, q in perm.items()):,} signatures permuted")
    return perts, sigs, fmap, perm


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


def profile_matrix(perts, perm, fmap, cfg):
    """Perturbations x mapped features: the full signed / scaled signature of the perturbation each one carries
    (perm), rank-transformed per perturbation, centred and L2-normalised, so that the dot product is Spearman."""
    cols = [perm[p] for p in perts.pert]
    header = pd.read_csv(cfg.signatures, nrows=0).columns  # the label file has every A549 perturbation: read
    want = set(cols)                                        # only the columns needed here
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
    X = val[cols].to_numpy(np.float32).T  # perturbations x features
    X = np.nan_to_num(X, nan=0.0)
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
            toks = [2 * gene, 2 * gene + 1] if cls == 3 else [2 * gene + cls - 1]
            for t in toks:
                r.append(i)
                c.append(t)
                v.append(0.5 if cls == 3 else 1.0)
    return sparse.csr_matrix((v, (r, c)), shape=(len(drugs), 2 * n_prot))


def knn_votes(sim, ref_unit, n_units, k):
    """sim: queries x reference perturbations; ref_unit[j]: the voting unit (drug or perturbation) of reference j.
    Per query: each unit's similarity is its maximum over its perturbations, the k most similar units vote with
    weight max(sim, 0), normalised (uniform over the k if none is positive). Returns queries x units, sparse."""
    U = np.full((sim.shape[0], n_units), -np.inf, np.float32)
    np.maximum.at(U.T, ref_unit, sim.T)  # unit similarity = max over its perturbations
    k = min(k, n_units)
    top = np.argpartition(-U, k - 1, axis=1)[:, :k]
    w = np.take_along_axis(U, top, axis=1)
    w = np.where(np.isfinite(w), np.maximum(w, 0), 0)
    s = w.sum(1, keepdims=True)
    w = np.where(s > 0, w / np.where(s > 0, s, 1), 1.0 / k)
    rows = np.repeat(np.arange(sim.shape[0]), k)
    return sparse.csr_matrix((w.ravel(), (rows, top.ravel())), shape=(sim.shape[0], n_units))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, help="a cross-validated run (or its _shuff control)")
    ap.add_argument("--representation", choices=("sources", "profile"), default="sources")
    ap.add_argument("--k", type=int, default=25, help="number of voting neighbours")
    ap.add_argument("--unit", choices=("drug", "pert"), default="drug", help="what a neighbour is")
    ap.add_argument("--ref-include-val", action="store_true",
                    help="also use the validation drugs as references (the model only early-stops on them)")
    ap.add_argument("--out-dir", help="default: <runs_dir>/<run>__knn_<representation>_k<k>[_pert][_shuff]")
    ap.add_argument("--allow-mismatch", action="store_true", help="continue when the inputs do not match the run")
    ap.add_argument("--max-test", type=int, default=0, help="smoke test: score only this many perturbations per fold")
    args = ap.parse_args()

    run = os.path.abspath(args.run_dir.rstrip("/"))
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
            if args.allow_mismatch and k == "de_counts":  # smoke test only: every perturbation counts as active
                log(f"warning: {k} {getattr(cfg, k)} not found, every perturbation treated as active")
                setattr(cfg, k, None)
            else:
                sys.exit(f"{k}: {getattr(cfg, k)} not found")

    name = os.path.basename(run)
    shuff = name.endswith("_shuff")
    base = name[:-len("_shuff")] if shuff else name
    tag = f"__knn_{args.representation}_k{args.k}" + ("_pert" if args.unit == "pert" else "")
    out = args.out_dir or os.path.join(os.path.dirname(run), base + tag + ("_shuff" if shuff else ""))
    os.makedirs(out, exist_ok=True)
    log(f"run {name} (shuffled signatures: {bool(cfg.shuffle_signatures)}), folds {folds}, mode {cfg.mode}; "
        f"{args.representation}, k {args.k} {args.unit}s -> {out}")

    # ---- inputs, as the run built them
    g = mtp.Graph(cfg.nodes, cfg.edges)
    perts_table = pd.read_csv(cfg.perturbations, sep="\t", dtype=str, keep_default_na=False)
    targets = mtp.load_drug_targets(cfg.drug_targets, g, perts_table)
    drug_mets = mtp.load_drug_metabolites(cfg.drug_metabolites, g) if cfg.drug_metabolites else {}
    perts, sigs, fmap, perm = filtered_perturbations(cfg, g, perts_table, targets)
    splits = mtp.drug_splits(perts.drug.unique(), cfg)
    log(f"{len(perts):,} perturbations of {perts.drug.nunique():,} drugs with targets "
        f"({int(perts.active.sum()):,} active); {len(splits)} splits")
    X = source_matrix(perts, sigs, fmap, drug_mets, cfg, g.N) if args.representation == "sources" else \
        profile_matrix(perts, perm, fmap, cfg)
    row = {p: j for j, p in enumerate(perts.pert)}
    n_prot = len(g.prot)

    check, all_pp, all_pd = {}, [], []
    for i in folds:
        tr, va, te = splits[i]
        fd = os.path.join(run, f"fold{i}")
        z = np.load(os.path.join(fd, "test_scores_heldout.npz"), allow_pickle=True)
        test = pd.read_csv(os.path.join(fd, "test_drug_ranks.tsv"), sep="\t")
        cands = [str(d) for d in z["drug"]]  # the run's own held-out candidates: authoritative
        rebuilt = sorted({d.lower() for d in te})
        same_cands = rebuilt == cands
        found = test.pert.isin(row)
        check[i] = dict(heldout_drugs_match=bool(same_cands), n_candidates=len(cands), n_test=int(len(test)),
                        n_test_found=int(found.sum()))
        log(f"fold {i}: held-out drugs rebuilt here {'match' if same_cands else 'DO NOT match'} the run's "
            f"candidates ({len(rebuilt)} vs {len(cands)}); {int(found.sum())}/{len(test)} test perturbations found")
        if not (same_cands and found.all()) and not args.allow_mismatch:
            sys.exit("inputs differ from the run's (labels, network, targets?): run where the run was trained, "
                     "or --allow-mismatch for a smoke test")
        test = test[found].reset_index(drop=True)
        if args.max_test:
            test = test.head(args.max_test)
        rk = mtp.DrugRanker({d: targets[d] for d in cands}, n_prot)

        # references: training drugs (optionally + validation), active perturbations
        # (every held-out candidate excluded explicitly, so not even a mismatched smoke run can leak)
        ref_drugs = sorted({d.lower() for d in (tr | va if args.ref_include_val else tr)} - set(cands))
        ref = perts[perts.drug.str.lower().isin(ref_drugs) & perts.active.astype(bool)]
        ref_idx = ref.pert.map(row).to_numpy()
        if args.unit == "drug":
            unit_names = ref_drugs
            pos = {d: j for j, d in enumerate(ref_drugs)}
            ref_unit = ref.drug.str.lower().map(pos).to_numpy()
            unit_drug = ref_drugs
        else:
            unit_names = list(ref.pert)
            ref_unit = np.arange(len(ref))
            unit_drug = list(ref.drug.str.lower())
        T = target_tokens(unit_drug, targets, n_prot)

        q_idx = test.pert.map(row).to_numpy()
        Q, R = X[q_idx], X[ref_idx]
        sim = (Q @ R.T)
        sim = np.asarray(sim.todense() if sparse.issparse(sim) else sim, np.float32)
        W = knn_votes(sim, ref_unit, len(unit_names), args.k)
        TOK = np.asarray((W @ T).todense())  # queries x (2 * proteins): weighted target votes
        S = np.empty((len(test), len(rk.drugs)), np.float32)
        for j in range(len(test)):
            prob = np.zeros((n_prot, 3))
            prob[:, 1], prob[:, 2] = TOK[j, 0::2], TOK[j, 1::2]
            S[j] = rk.prob_scores(prob)

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
        pdr = pd.DataFrame(drows)
        dS = np.asarray(dS)
        pdr["rank_prob_heldout"] = [rk.rank_of(dS[j], d) for j, d in enumerate(pdr.drug)]
        for subset in ("all", "active"):
            m = (pdr.subset == subset).to_numpy()
            if m.any():
                k_, grp = pdr.drug[m].to_numpy(), np.ones(m.sum(), bool)
                sp, n = mtp.signature_specificity(dS[m], k_, rk.col, grp)
                pdr.loc[m, "spec_prob_heldout"], pdr.loc[m, "n_spec_heldout"] = sp, n
                pdr.loc[m, "rank_probc_heldout"] = mtp.centered_ranks(dS[m], k_, rk.col, grp)

        od = os.path.join(out, f"fold{i}")
        os.makedirs(od, exist_ok=True)
        pp.to_csv(os.path.join(od, "test_drug_ranks.tsv"), sep="\t", index=False)
        pdr.to_csv(os.path.join(od, "test_drug_ranks_per_drug.tsv"), sep="\t", index=False)
        np.savez_compressed(os.path.join(od, "test_scores_heldout.npz"), pert=pp.pert.to_numpy(),
                            drug=np.array(rk.drugs), prob=S)
        med = lambda c, n: float(np.median(pp[c] / pp[n]))
        fold_m = dict(n_test=len(pp), n_references=int(len(ref)), n_units=len(unit_names),
                      rank_prob_heldout=med("rank_prob_heldout", "n_candidates_heldout"),
                      rank_probc_heldout=med("rank_probc_heldout", "n_candidates_heldout"),
                      spec_prob_heldout=float(np.nanmedian(pp.spec_prob_heldout / pp.n_spec_heldout)))
        with open(os.path.join(od, "metrics.json"), "w") as fh:
            json.dump(fold_m, fh, indent=2)
        log(f"fold {i}: {len(pp)} queries, {len(ref):,} reference perturbations of {len(ref_drugs)} drugs; median "
            f"percentile rank_prob {fold_m['rank_prob_heldout']:.3f}, rank_probc {fold_m['rank_probc_heldout']:.3f}, "
            f"spec {fold_m['spec_prob_heldout']:.3f}")
        all_pp.append(pp.assign(fold=i))
        all_pd.append(pdr.assign(fold=i))

    pd.concat(all_pp, ignore_index=True).to_csv(os.path.join(out, "cv_drug_ranks.tsv"), sep="\t", index=False)
    pd.concat(all_pd, ignore_index=True).to_csv(os.path.join(out, "cv_drug_ranks_per_drug.tsv"), sep="\t",
                                                index=False)
    with open(os.path.join(out, "config.json"), "w") as fh:
        json.dump(dict(vars(args), run=run, folds=folds, input_check=check,
                       run_config={k: v for k, v in vars(cfg).items()}), fh, indent=2)
    log(f"saved {out}/cv_drug_ranks.tsv (+ per drug, per fold); compare with its _shuff counterpart via "
        f"compare_runs.py or plot_performance_panels.py --extra")


if __name__ == "__main__":
    main()
