#!/usr/bin/env python
"""
Annotations used by additional_analyses.R (and section 10 of metchange_vs_expression_benchmark.R).

Writes to results/mc_schuhk_assoc/:
  drug_mechanism_annotation.csv  per Schuhknecht-matched drug (LINCS name): number of targets, whether a target is a
                                 GEM enzyme (catalysis source in the network), whether it is a GEM enzyme the drug
                                 inhibits, and the distance of its closest target to the transcriptional route
                                 (directed path to a TF with a signed regulon edge to a GEM enzyme). Columns are given
                                 for all targets and, prefixed mech_, for mechanism-of-action targets only (pairs
                                 supported by ChEMBL, DrugCentral or IUPHAR): DrugBank-derived entries also list
                                 pharmacokinetic enzymes (CYPs, UGTs) as targets.
  metabolite_degree.csv          per GEM metabolite: number of reactions producing / consuming it

Usage (from scripts/mc_schuhk_assoc):  python annotate_drugs_metabolites.py
"""
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts", "build_network"))
import build_network as bn  # noqa: E402

DATA, RES = os.path.join(ROOT, "data"), os.path.join(ROOT, "results", "mc_schuhk_assoc")


def drug_annotation(max_hops=5):
    nodes = bn.read_table(os.path.join(DATA, "networks", "combined", "nodes.tsv"))
    edges = pd.read_csv(os.path.join(DATA, "networks", "combined", "edges.tsv"), sep="\t", keep_default_na=False,
                        dtype={"src": str, "dst": str})
    hgnc = bn.HGNC(os.path.join(DATA, "utilities", "hgnc_complete_set.txt"))
    dt = bn.load_drug_targets(os.path.join(DATA, "drug_targets", "combined", "drug_targets.tsv"), hgnc)
    dist, _ = bn.transcriptional_route(nodes, edges)
    dt["distance"] = dt.target.map(pd.Series(dist, index=nodes.node_id))
    enzymes = set(edges.src[edges.relation == "catalysis"])
    dt["gem"] = dt.target.isin(enzymes)
    dt["compound"] = dt.pertname.str.lower()

    dt["mech"] = dt.sources.str.contains("ChEMBL|DrugCentral|IUPHAR")

    def summarise(d, prefix):
        g = d.groupby("compound")
        ann = pd.DataFrame({
            "n_targets": g.target.nunique(),
            "targets_gem": g.gem.any(),
            "direct_gem_inhibitor": g.apply(lambda x: (x.gem & (x.effect == "inh")).any(), include_groups=False),
            "min_route_distance": g.distance.min(),
        })
        ann["on_route"] = ann.min_route_distance <= max_hops
        return ann.add_prefix(prefix)

    dm = pd.read_csv(os.path.join(RES, "rankAllPenlty", "drug_match_table.csv"), index_col=0)
    out = dm[["compound", "schuhknecht_name"]].reset_index(drop=True)
    key = out.compound.str.lower()
    for prefix, d in (("", dt), ("mech_", dt[dt.mech])):
        ann = summarise(d, prefix).reindex(key).reset_index(drop=True)
        has = ann[f"{prefix}n_targets"].notna()
        out[f"{prefix}has_target"] = has
        for c in ("targets_gem", "direct_gem_inhibitor", "on_route"):  # NA for drugs without a known target
            ann[prefix + c] = ann[prefix + c].astype("boolean").mask(~has)
        out = pd.concat([out, ann], axis=1)
    out.to_csv(os.path.join(RES, "drug_mechanism_annotation.csv"), index=False)
    for prefix in ("", "mech_"):
        n = lambda c: int(out[prefix + c].fillna(False).sum())
        print(f"{prefix or 'all_'}targets: {len(out)} drugs; with targets {n('has_target')}; GEM-enzyme target "
              f"{n('targets_gem')}; direct GEM inhibitor {n('direct_gem_inhibitor')}; on route {n('on_route')}")


def metabolite_degree():
    gem = bn.load_gem(os.path.join(DATA, "GEMs", "Human-GEM", "model", "Human-GEM.xml"))
    deg = pd.DataFrame({
        "id": [m.id for m in gem.metabolites],
        "n_prod": [sum(r.metabolites[m] > 0 for r in m.reactions) for m in gem.metabolites],
        "n_reac": [sum(r.metabolites[m] < 0 for r in m.reactions) for m in gem.metabolites],
    })
    deg.to_csv(os.path.join(RES, "metabolite_degree.csv"), index=False)
    print(f"metabolite_degree.csv: {len(deg)} metabolites")


if __name__ == "__main__":
    drug_annotation()
    metabolite_degree()
