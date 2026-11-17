#!/usr/bin/env python
"""
Matched MetChange / expression / Schuhknecht matrices for the association analyses of this folder
(metchange_vs_expression_benchmark.R, sign_heterogeneity_analysis.R, overlap_enrichment.py), from the limma labels
(compute_labels: log2FC.csv / qvals.csv, features x perturbations) instead of the averaged deltas used before.

Drugs: the LINCS <-> Schuhknecht name pairs of drug_match_table.csv (SMILES matching, scripts/schuhknecht_map). For
each drug the LINCS perturbation (label column '<drug>_<dose> uM_<time> h') closest to the Schuhknecht condition is
chosen as in create_drugmap.py: the time closest to --target-time, then the dose closest to --target-dose on a log
scale within --max-log-dose-dist; statuses exact / nearest_dose / nearest_time / nearest_time_and_dose, kept if in
--match-status.

Metabolites: GEM metabolites mapped to one Schuhknecht ion (cpd_info_wSchuhknecht.csv, dup_after_shuhknecht_map ==
False), in --compartments (default c, as the earlier analysis; 'intracellular' = every compartment except e). Rows
are GEM metabolite names; when a name has several compartment copies among the chosen ones, MetChange log2FC and the
expression baseline are averaged over the copies and the MetChange q-value is the minimum over them (n_copies in
metabolites.tsv).

Values (no sign change here; the R script flips MetChange):
  mc_log2FC_matched.csv    MetChange limma log2FC (inconsistency: positive = LESS production)
  mc_qvals_matched.csv     MetChange limma q-values
  expr_<baseline>_matched.csv  expression baseline from the expression limma log2FC: per reaction the mean log2FC
                           of its GPR genes (gpr_ensembl of the network's reaction nodes; logic ignored), per
                           metabolite the mean over the reactions producing it (products) or producing minus consuming
                           (net), from the network's reaction -> metabolite edges (reversible reactions both ways,
                           currency metabolites already removed). Positive = more production.
  sch_log2FC_matched.csv, sch_pvals_matched.csv, sch_qvals_matched.csv   Schuhknecht mean log2-FC, p, q
  drug_match_table.csv     compound, schuhknecht_name, matched_column, log_dose_distance, time_distance_h,
                           match_status, validated_metabolic, n_sig_sch, n_sig_mc (q <= 0.05 among the rows)
  metabolites.tsv          name, ids, n_copies, schuhknecht_mapping, degree (producing + consuming reactions)
  metabolite_degree.csv    metabolite, degree (section 10 of the benchmark)
All matrices: rows = metabolite names, columns = Schuhknecht drug names, same order everywhere.
"""
import argparse
import os
import re
import time

import numpy as np
import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
COL_RE = re.compile(r"^(.*)_([\d.]+) uM_(\d+) h$")


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def best_column(cands, target_dose, target_time, max_ld):
    """cands: [(column, dose, time)] of one drug -> (column, log dose distance, time distance, status)."""
    times = sorted({t for _, _, t in cands})
    t_best = target_time if target_time in times else min(times, key=lambda t: abs(t - target_time))
    same = [(c, d) for c, d, t in cands if t == t_best]
    col, ld = min(((c, abs(np.log10(d) - np.log10(target_dose))) for c, d in same), key=lambda x: x[1])
    td = abs(t_best - target_time)
    if ld > max_ld:
        return None, ld, td, "dose_too_far"
    status = "exact" if ld == 0 and td == 0 else "nearest_dose" if td == 0 else \
        "nearest_time" if ld == 0 else "nearest_time_and_dose"
    return col, ld, td, status


def read_cols(path, cols, rows=None):
    """Selected columns of a features x perturbations CSV (first column = feature id), float32."""
    hdr = pd.read_csv(path, nrows=0).columns
    df = pd.read_csv(path, usecols=[hdr[0]] + list(cols), index_col=0,
                     dtype={c: np.float32 for c in cols})
    df.index = df.index.astype(str)
    return df if rows is None else df.reindex(rows)


def read_schuhknecht(xlsx, cache_dir):
    out = {}
    for sheet, key in (("mean log2-FC", "lfc"), ("p-values", "p"), ("q-values", "q")):
        f = os.path.join(cache_dir, f"schuhknecht_{key}.pkl")
        if os.path.exists(f) and os.path.getmtime(f) > os.path.getmtime(xlsx):
            out[key] = pd.read_pickle(f)
            continue
        x = pd.read_excel(xlsx, sheet_name=sheet)
        x = x.set_index(x.columns[0]).iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
        x.index = x.index.astype(str)
        x.to_pickle(f)
        log(f"read Schuhknecht sheet '{sheet}': {x.shape[0]} ions x {x.shape[1]} drugs")
        out[key] = x
    return out["lfc"], out["p"], out["q"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mc-dir", required=True, help="MetChange limma labels (log2FC.csv, qvals.csv)")
    ap.add_argument("--expr-dir", help="expression limma labels (log2FC.csv); without it no expression baseline")
    ap.add_argument("--schuhknecht-dir", default=os.path.join(REPO, "data", "validation", "schuhknecht"))
    ap.add_argument("--network-dir", default=os.path.join(REPO, "data", "networks", "combined"),
                    help="build_network.py output (nodes.tsv, edges.tsv): GPRs and reaction -> metabolite edges")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--compartments", default="c", help="comma-separated GEM compartments, or 'intracellular'")
    ap.add_argument("--baseline", choices=["products", "net"], default="products")
    ap.add_argument("--target-dose", type=float, default=10.0, help="Schuhknecht dose (uM)")
    ap.add_argument("--target-time", type=float, default=48.0, help="Schuhknecht time (h)")
    ap.add_argument("--max-log-dose-dist", type=float, default=float(np.log10(3)))
    ap.add_argument("--match-status", default="exact,nearest_time",
                    help="kept statuses (exact, nearest_dose, nearest_time, nearest_time_and_dose)")
    ap.add_argument("--alpha", type=float, default=0.05, help="q-value cutoff for the n_sig columns")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    sd = args.schuhknecht_dir

    # ---- drugs and perturbation columns
    dm0 = pd.read_csv(os.path.join(sd, "drug_match_table.csv"), index_col=0)
    dm0 = dm0.drop_duplicates("compound")
    mc_cols = pd.read_csv(os.path.join(args.mc_dir, "log2FC.csv"), nrows=0).columns[1:]
    parsed = {}
    for c in mc_cols:
        m = COL_RE.match(c)
        if m:
            parsed.setdefault(m.group(1).lower(), []).append((c, float(m.group(2)), float(m.group(3))))
    rows = []
    for _, r in dm0.iterrows():
        cands = parsed.get(r.compound.lower())
        if not cands:
            rows.append(dict(compound=r.compound, schuhknecht_name=r.schuhknecht_name, matched_column=None,
                             match_status="no_data_for_compound"))
            continue
        col, ld, td, st = best_column(cands, args.target_dose, args.target_time, args.max_log_dose_dist)
        rows.append(dict(compound=r.compound, schuhknecht_name=r.schuhknecht_name, matched_column=col,
                         log_dose_distance=ld, time_distance_h=td, match_status=st))
    dm = pd.DataFrame(rows)
    dm["validated_metabolic"] = dm.schuhknecht_name.map(dict(zip(dm0.schuhknecht_name, dm0.validated_metabolic)))
    log(f"drugs: {len(dm)} LINCS-Schuhknecht pairs; match status {dm.match_status.value_counts().to_dict()}")
    keep_st = set(args.match_status.split(","))
    dm = dm[dm.match_status.isin(keep_st) & dm.matched_column.notna()].reset_index(drop=True)
    log(f"kept ({', '.join(sorted(keep_st))}): {len(dm)} drugs")

    lfc, pv, qv = read_schuhknecht(os.path.join(sd, "41587_2024_2524_MOESM2_ESM.xlsx"), args.out_dir)
    dm = dm[dm.schuhknecht_name.isin(lfc.columns)].reset_index(drop=True)

    # ---- metabolites
    cpd = pd.read_csv(os.path.join(sd, "cpd_info_wSchuhknecht.csv"), index_col=0)
    cpd = cpd[(cpd.dup_after_shuhknecht_map == False) & cpd.Schuhknecht_mapping.notna()]  # noqa: E712
    comps = None if args.compartments == "intracellular" else set(args.compartments.split(","))
    cpd = cpd[(cpd.compartment != "e") if comps is None else cpd.compartment.isin(comps)]
    mc_rows = set(pd.read_csv(os.path.join(args.mc_dir, "log2FC.csv"), usecols=[0]).iloc[:, 0].astype(str))
    cpd = cpd[cpd.id.isin(mc_rows) & cpd.Schuhknecht_mapping.isin(lfc.index)]
    mets = cpd.groupby("name").agg(ids=("id", lambda x: "|".join(sorted(x))), n_copies=("id", "size"),
                                   schuhknecht_mapping=("Schuhknecht_mapping", "first"),
                                   n_ions=("Schuhknecht_mapping", "nunique")).reset_index()
    if (mets.n_ions > 1).any():
        log(f"warning: {int((mets.n_ions > 1).sum())} metabolite names map to several ions; first kept")
    mets = mets.drop(columns="n_ions")
    log(f"metabolites: {len(mets)} names ({len(cpd)} GEM ids, compartments {args.compartments}) mapped to one "
        f"Schuhknecht ion and present in the MetChange labels")
    id2name = dict(zip(cpd.id, cpd.name))
    names = list(mets.name)
    cols = list(dm.matched_column)

    def collapse(df, how):  # GEM ids x columns -> names x Schuhknecht drugs
        g = df.groupby(df.index.map(id2name))
        out = g.mean() if how == "mean" else g.min()
        out = out.reindex(names)
        out.columns = list(dm.schuhknecht_name)
        return out

    # ---- MetChange
    t0 = time.time()
    ids = list(cpd.id)
    mc_l = collapse(read_cols(os.path.join(args.mc_dir, "log2FC.csv"), cols, ids), "mean")
    mc_q = collapse(read_cols(os.path.join(args.mc_dir, "qvals.csv"), cols, ids), "min")
    log(f"MetChange: {mc_l.shape[0]} x {mc_l.shape[1]} ({time.time() - t0:.0f}s); "
        f"{int(mc_l.isna().all(axis=1).sum())} metabolites all-NaN")

    # ---- Schuhknecht
    def sch(x):
        out = x.loc[mets.schuhknecht_mapping, list(dm.schuhknecht_name)]
        out.index = names
        return out
    s_l, s_p, s_q = sch(lfc), sch(pv), sch(qv)

    # ---- network: degree, expression baseline
    nodes = pd.read_csv(os.path.join(args.network_dir, "nodes.tsv"), sep="\t", dtype=str, keep_default_na=False)
    edges = pd.read_csv(os.path.join(args.network_dir, "edges.tsv"), sep="\t", dtype=str, keep_default_na=False,
                        usecols=["src", "dst", "relation"])
    e = edges[edges.relation == "metabolic"]
    prod = e[e.src.str.startswith("MAR") & e.dst.isin(ids)].rename(columns={"src": "rxn", "dst": "met"})
    cons = e[e.dst.str.startswith("MAR") & e.src.isin(ids)].rename(columns={"dst": "rxn", "src": "met"})
    deg = pd.concat([prod, cons]).groupby("met").rxn.size()
    mets["degree"] = [int(sum(deg.get(i, 0) for i in s.split("|"))) for s in mets.ids]

    files = {"mc_log2FC_matched.csv": mc_l, "mc_qvals_matched.csv": mc_q, "sch_log2FC_matched.csv": s_l,
             "sch_pvals_matched.csv": s_p, "sch_qvals_matched.csv": s_q}
    if args.expr_dir:
        t0 = time.time()
        rx = nodes[nodes.node_type == "reaction"][["node_id", "gpr_ensembl"]]
        rx_genes = {r: sorted(set(re.findall(r"ENSG\d+", g))) for r, g in zip(rx.node_id, rx.gpr_ensembl) if g}
        need = sorted({g for gs in rx_genes.values() for g in gs})
        ecols = pd.read_csv(os.path.join(args.expr_dir, "log2FC.csv"), nrows=0).columns
        miss = [c for c in cols if c not in set(ecols)]
        if miss:
            log(f"warning: {len(miss)} matched perturbations missing in the expression labels (NaN baseline)")
        ex = read_cols(os.path.join(args.expr_dir, "log2FC.csv"), [c for c in cols if c not in miss])
        ex = ex[ex.index.isin(need)]
        log(f"expression: {ex.shape[0]} GEM genes of {len(need)} in the GPRs")
        gene_pos = {g: i for i, g in enumerate(ex.index)}
        X = ex.to_numpy(np.float64)
        rscore = {}
        for r, gs in rx_genes.items():
            idx = [gene_pos[g] for g in gs if g in gene_pos]
            if idx:
                rscore[r] = np.nanmean(X[idx], 0)
        R = pd.DataFrame(rscore, index=ex.columns).T

        def met_mean(pairs):
            p = pairs[pairs.rxn.isin(R.index)]
            return R.loc[p.rxn].set_axis(p.met.to_numpy()).groupby(level=0).mean()
        base = met_mean(prod)
        if args.baseline == "net":
            base = base.sub(met_mean(cons), fill_value=np.nan)
        base = base.reindex(columns=cols)
        exp_m = collapse(base.reindex([i for i in ids if i in base.index]), "mean")
        files[f"expr_{args.baseline}_matched.csv"] = exp_m
        log(f"expression baseline ({args.baseline}): {int(exp_m.notna().any(axis=1).sum())} of {len(names)} "
            f"metabolites with a value ({time.time() - t0:.0f}s)")

    dm["n_sig_sch"] = (s_q <= args.alpha).sum().to_numpy()
    dm["n_sig_mc"] = (mc_q <= args.alpha).sum().to_numpy()
    for f, df in files.items():
        df.to_csv(os.path.join(args.out_dir, f), index_label="metabolite")
    dm.to_csv(os.path.join(args.out_dir, "drug_match_table.csv"), index=False)
    mets.to_csv(os.path.join(args.out_dir, "metabolites.tsv"), sep="\t", index=False)
    mets[["name", "degree"]].rename(columns={"name": "metabolite"}).to_csv(
        os.path.join(args.out_dir, "metabolite_degree.csv"), index=False)
    log(f"MetChange significant (q <= {args.alpha}) per drug: median {int(dm.n_sig_mc.median())}, "
        f"{int((dm.n_sig_mc == 0).sum())} of {len(dm)} drugs with none; Schuhknecht: median "
        f"{int(dm.n_sig_sch.median())}")
    log(f"written to {args.out_dir}")


if __name__ == "__main__":
    main()
