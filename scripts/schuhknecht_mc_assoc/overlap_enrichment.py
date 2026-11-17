#!/usr/bin/env python
"""
Do the metabolites MetChange calls significant for a drug overlap the ones Schuhknecht measured as significantly
changed for the same drug, more than expected, and with the same direction?

Input: the output folder of prepare_assoc_inputs.py. Per drug, over the metabolites with values in both datasets:
  A = MetChange q <= --alpha-mc (direction: production = -log2FC of the inconsistency score)
  B = Schuhknecht q <= --alpha-sch (direction: sign of the measured log2FC)
Per drug: 2 x 2 Fisher test of A vs B (one-sided, enrichment), odds ratio, overlap |A & B|, and among the overlap the
number of concordant (same direction) and discordant metabolites.

Across drugs (drugs with |A| >= --min-mc and |B| >= --min-sch):
  * Mantel-Haenszel common odds ratio and CMH test (strata = drugs). Assumes metabolites are exchangeable within a
    drug, which they are not (some metabolites are often significant in both datasets whatever the drug), so it is
    descriptive.
  * Mismatched-drug permutation (the test to rely on): the MetChange calls of drug i are compared with the
    Schuhknecht calls of a random other drug, keeping both sets intact. Statistics: total overlap sum_i |A_i & B_i|,
    directional score sum_i (concordant_i - discordant_i), and the number of concordant overlaps. p = P(null >=
    observed). This keeps each metabolite's tendency to be called in either dataset, so it tests drug specificity.
  * Expression baseline, size-matched: per drug, the |A_i| metabolites with the largest |expression baseline| and its
    sign, scored the same way (overlap, directional score, permutation p), and MetChange vs expression per drug
    (paired Wilcoxon on the directional score). Answers whether MetChange's significance adds to what the
    expression of the producing enzymes gives.

Outputs (--out-dir, default <in-dir>/overlap_enrichment): per_drug.tsv, summary.tsv, permutation_null.tsv.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, wilcoxon


def read(d, f):
    return pd.read_csv(os.path.join(d, f), index_col=0)


def mh(tables):
    """Mantel-Haenszel common OR and CMH chi-square (1 df, continuity-corrected) of 2x2 tables [[a, b], [c, d]]."""
    from scipy.stats import chi2
    num = den = s_a = s_e = s_v = 0.0
    for (a, b), (c, d) in tables:
        n = a + b + c + d
        if n < 2:
            continue
        num += a * d / n
        den += b * c / n
        r1, c1 = a + b, a + c
        s_a += a
        s_e += r1 * c1 / n
        s_v += r1 * (n - r1) * c1 * (n - c1) / (n * n * (n - 1))
    or_mh = num / den if den > 0 else np.inf
    x2 = (abs(s_a - s_e) - 0.5) ** 2 / s_v if s_v > 0 else np.nan
    return or_mh, s_a, s_e, x2, (chi2.sf(x2, 1) if s_v > 0 else np.nan)


def stats_for(callA, dirA, B, dirB, valid, perm=None):
    """callA, B: bool (mets x drugs); dirA, dirB: sign arrays; valid: bool; perm: column order of B -> per drug
    overlap, concordant, discordant."""
    if perm is not None:
        B, dirB, valid_b = B[:, perm], dirB[:, perm], valid[1][:, perm]
    else:
        valid_b = valid[1]
    v = valid[0] & valid_b
    ov = callA & B & v
    conc = (ov & (dirA == dirB)).sum(0)
    disc = (ov & (dirA != dirB)).sum(0)
    return ov.sum(0), conc, disc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", required=True, help="prepare_assoc_inputs.py output")
    ap.add_argument("--out-dir")
    ap.add_argument("--alpha-mc", type=float, default=0.05)
    ap.add_argument("--alpha-sch", type=float, default=0.05)
    ap.add_argument("--min-mc", type=int, default=1, help="drugs need at least this many MetChange calls")
    ap.add_argument("--min-sch", type=int, default=1, help="... and this many Schuhknecht calls")
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    out = args.out_dir or os.path.join(args.in_dir, "overlap_enrichment")
    os.makedirs(out, exist_ok=True)

    mc_l, mc_q = read(args.in_dir, "mc_log2FC_matched.csv"), read(args.in_dir, "mc_qvals_matched.csv")
    s_l, s_q = read(args.in_dir, "sch_log2FC_matched.csv"), read(args.in_dir, "sch_qvals_matched.csv")
    exp_f = sorted(glob.glob(os.path.join(args.in_dir, "expr_*_matched.csv")))
    ex = read(args.in_dir, os.path.basename(exp_f[0])).reindex_like(mc_l) if exp_f else None
    dm = pd.read_csv(os.path.join(args.in_dir, "drug_match_table.csv"))
    drugs = list(mc_l.columns)
    assert list(s_l.columns) == drugs and list(mc_l.index) == list(s_l.index)

    va = (mc_l.notna() & mc_q.notna()).to_numpy()
    vb = (s_l.notna() & s_q.notna()).to_numpy()
    A = va & (mc_q.to_numpy() <= args.alpha_mc)
    dirA = np.sign(-mc_l.fillna(0).to_numpy())  # inconsistency -> production
    B = vb & (s_q.to_numpy() <= args.alpha_sch)
    dirB = np.sign(s_l.fillna(0).to_numpy())
    v_both = va & vb
    nA, nB = (A & vb).sum(0), (B & va).sum(0)
    inc = (nA >= args.min_mc) & (nB >= args.min_sch)
    print(f"{len(drugs)} drugs, {mc_l.shape[0]} metabolites; MetChange calls per drug: median {int(np.median(nA))}, "
          f"{int((nA == 0).sum())} with none; Schuhknecht: median {int(np.median(nB))}; included (>= {args.min_mc} / "
          f">= {args.min_sch}): {int(inc.sum())}")
    if inc.sum() < 5:
        print("fewer than 5 drugs with calls in both datasets: nothing to test")

    # expression, size-matched calls
    if ex is not None:
        E_abs = ex.abs().to_numpy()
        E = np.zeros_like(A)
        for j in range(len(drugs)):
            k = int(nA[j])
            ok = np.where(va[:, j] & vb[:, j] & ~np.isnan(E_abs[:, j]))[0]
            if k and len(ok):
                E[ok[np.argsort(-E_abs[ok, j])[:k]], j] = True
        dirE = np.sign(ex.fillna(0).to_numpy())
        vE = va & ~np.isnan(E_abs)

    # per drug
    rows = []
    ov, conc, disc = stats_for(A, dirA, B, dirB, (va, vb))
    tables = []
    for j, d in enumerate(drugs):
        m = v_both[:, j]
        a = int((A[:, j] & B[:, j] & m).sum())
        b = int((A[:, j] & ~B[:, j] & m).sum())
        c = int((~A[:, j] & B[:, j] & m).sum())
        dd = int((~A[:, j] & ~B[:, j] & m).sum())
        r = dict(drug=d, n_met=int(m.sum()), n_mc=a + b, n_sch=a + c, overlap=a, concordant=int(conc[j]),
                 discordant=int(disc[j]), expected=(a + b) * (a + c) / max(m.sum(), 1), included=bool(inc[j]))
        if a + b and a + c:
            r["odds_ratio"], r["p_fisher"] = fisher_exact([[a, b], [c, dd]], alternative="greater")
        if inc[j]:
            tables.append(((a, b), (c, dd)))
        rows.append(r)
    pdrug = pd.DataFrame(rows).merge(dm[["schuhknecht_name", "compound", "matched_column", "validated_metabolic"]],
                                     left_on="drug", right_on="schuhknecht_name", how="left").drop(
        columns="schuhknecht_name")
    if "p_fisher" in pdrug:
        ok = pdrug.p_fisher.notna()
        p_ = pdrug.p_fisher[ok].to_numpy()
        o = np.argsort(p_)
        q_ = np.minimum.accumulate((p_[o] * len(p_) / np.arange(1, len(p_) + 1))[::-1])[::-1]  # BH
        q = np.empty_like(q_)
        q[o] = np.minimum(q_, 1)
        pdrug.loc[ok, "q_fisher"] = q
    if ex is not None:
        e_ov, e_conc, e_disc = stats_for(E, dirE, B, dirB, (vE, vb))
        pdrug["expr_overlap"], pdrug["expr_concordant"], pdrug["expr_discordant"] = e_ov, e_conc, e_disc
    pdrug.to_csv(os.path.join(out, "per_drug.tsv"), sep="\t", index=False)

    # across drugs
    idx = np.where(inc)[0]
    summ = []
    if len(idx) >= 5:
        or_mh, s_a, s_e, x2, p_cmh = mh(tables)
        summ.append(dict(test="CMH (descriptive)", n_drugs=len(idx), observed=s_a, expected=s_e,
                         statistic=or_mh, p=p_cmh, note="common odds ratio; metabolites assumed exchangeable"))
        rng = np.random.default_rng(args.seed)
        sub = lambda x: x[:, idx]  # noqa: E731
        Ai, dAi, Bi, dBi, vai, vbi = map(sub, (A, dirA, B, dirB, va, vb))
        scores = {"MetChange": (Ai, dAi, vai)}
        if ex is not None:
            scores["expression (size-matched)"] = tuple(map(sub, (E, dirE, vE)))
        null_rows = []
        for name, (S, dS, vS) in scores.items():
            o, c_, d_ = stats_for(S, dS, Bi, dBi, (vS, vbi))
            obs = dict(overlap=o.sum(), directional=(c_ - d_).sum(), concordant=c_.sum())
            null = {k: np.empty(args.n_perm) for k in obs}
            n = len(idx)
            for t in range(args.n_perm):
                perm = rng.permutation(n)
                for k in np.where(perm == np.arange(n))[0]:  # no drug paired with itself: swap with another
                    s = rng.integers(n - 1)
                    s += s >= k
                    perm[k], perm[s] = perm[s], perm[k]  # cannot create a new fixed point
                o2, c2, d2 = stats_for(S, dS, Bi, dBi, (vS, vbi), perm)
                null["overlap"][t], null["directional"][t], null["concordant"][t] = o2.sum(), (c2 - d2).sum(), \
                    c2.sum()
            for k in obs:
                summ.append(dict(test=f"mismatched-drug permutation: {k}", score=name, n_drugs=n,
                                 observed=float(obs[k]), expected=float(null[k].mean()),
                                 statistic=float((obs[k] - null[k].mean()) / (null[k].std() + 1e-12)),
                                 p=float((1 + np.sum(null[k] >= obs[k])) / (1 + args.n_perm)),
                                 note="statistic = z vs null"))
                null_rows.append(pd.DataFrame({"score": name, "statistic": k, "value": null[k]}))
        pd.concat(null_rows).to_csv(os.path.join(out, "permutation_null.tsv"), sep="\t", index=False)
        if ex is not None:
            p_ = pdrug.iloc[idx]
            dmc = p_.concordant - p_.discordant
            dex = p_.expr_concordant - p_.expr_discordant
            dif = dmc - dex
            summ.append(dict(test="MetChange - expression, directional score per drug (paired Wilcoxon)",
                             n_drugs=len(idx), observed=float(dmc.mean()), expected=float(dex.mean()),
                             statistic=float(dif.mean()),
                             p=float(wilcoxon(dmc, dex).pvalue) if (dif != 0).any() else np.nan,
                             note="observed = mean MetChange, expected = mean expression, statistic = mean diff"))
        tot_c, tot_d = int(pdrug.concordant[inc].sum()), int(pdrug.discordant[inc].sum())
        summ.append(dict(test="pooled sign concordance among overlaps (descriptive)", n_drugs=len(idx),
                         observed=tot_c / max(tot_c + tot_d, 1), expected=0.5, statistic=tot_c + tot_d,
                         note="fraction concordant; statistic = number of overlapping calls"))
    summary = pd.DataFrame(summ)
    summary.to_csv(os.path.join(out, "summary.tsv"), sep="\t", index=False)
    with pd.option_context("display.width", 220, "display.max_columns", 20, "display.max_colwidth", 70):
        print(summary.round(4).to_string(index=False))
    if "q_fisher" in pdrug:
        print(f"per-drug Fisher: {int((pdrug.q_fisher <= 0.1).sum())} drugs with q <= 0.1 of "
              f"{int(pdrug.q_fisher.notna().sum())} testable")
    print(f"written to {out}")


if __name__ == "__main__":
    main()
