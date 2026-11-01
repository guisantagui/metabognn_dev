#!/usr/bin/env python
"""
Build the drug -> target table, with the direction of the effect, for the LINCS perturbagens.

Inputs (paths relative to --data-dir, i.e. metabognn_dev/data)
  LINCS/lincs_unique_perts.tsv                                   pert_id (BRD), pertname, smiles, inchi_key
  drug_targets/chembl/mechanism_*.json, chembl_uniprot_mapping.txt, molecules.tsv (dl_chembl_molecules.py)
  drug_targets/drug_central/drug.target.interaction.tsv.gz
  drug_targets/IUPHAR/interactions.csv
  drug_targets/chempert/targets_of_perturbagens.tsv              (includes DrugBank, CLUE, STITCH)
  drug_targets/broad_drug_repurposing/repurposing_drugs_20200324.txt + repurposing_samples_20200324.txt
  utilities/hgnc_complete_set.txt                                targets -> HGNC approved symbols

Matching a source compound to LINCS perturbagens, in this order:
  1. BRD id (first 13 characters, Broad only)      2. InChIKey
  3. InChIKey first block (same skeleton; differs in stereochemistry/charge)   4. normalised name
Names include aliases gathered on the way: the Repurposing Hub name of a matched BRD id and the ChEMBL preferred
name and synonyms of a matched InChIKey, so name-only sources (DrugCentral, IUPHAR, ChemPert) also match LINCS
perturbagens known under a code name.

Direction of each (perturbagen, target): act / inh / unknown. Sources are tiered:
  tier 1  ChEMBL, DrugCentral, IUPHAR       curated mechanism with an explicit action type per target
  tier 2  ChemPert entries from DrugBank or Ramilowski
  tier 3  Repurposing Hub (one MoA text for all listed targets) and ChemPert entries from CLUE only
The direction is taken from the best tier that gives one; if that tier contains both act and inh, it is unknown.
A name that matches several distinct LINCS compounds (e.g. a multivitamin product whose synonyms include
cyanocobalamin and biotin) is rejected as ambiguous. Pairs without a direction are kept only when a tier-1 source
reports them or at least two sources do (--keep-weak-unknown keeps them all): DrugBank lists carriers and
plasma binding proteins as targets without an action type.
ChemPert entries backed by STITCH alone (predicted / text-mined associations) are excluded unless
--keep-stitch-only.

Outputs (--out-dir)
  drug_targets.tsv     drug (LINCS pert_id), target, effect, pertname, sources, n_sources, evidence, conflict
  source_summary.tsv   per source: rows used, targets mapped to HGNC, rows and perturbagens matched to LINCS
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build_network"))
from build_network import HGNC, read_table  # noqa: E402

DEFAULT_PATHS = dict(
    lincs="LINCS/lincs_unique_perts.tsv",
    hgnc="utilities/hgnc_complete_set.txt",
    chembl_dir="drug_targets/chembl",
    drugcentral="drug_targets/drug_central/drug.target.interaction.tsv.gz",
    iuphar="drug_targets/IUPHAR/interactions.csv",
    chempert="drug_targets/chempert/targets_of_perturbagens.tsv",
    broad_drugs="drug_targets/broad_drug_repurposing/repurposing_drugs_20200324.txt",
    broad_samples="drug_targets/broad_drug_repurposing/repurposing_samples_20200324.txt",
)

# Action vocabularies of ChEMBL, DrugCentral and IUPHAR (upper-cased) -> direction on the target
ACT_TERMS = {"AGONIST", "FULL AGONIST", "PARTIAL AGONIST", "BIASED AGONIST", "IRREVERSIBLE AGONIST", "ACTIVATOR",
             "ACTIVATION", "POSITIVE ALLOSTERIC MODULATOR", "POSITIVE MODULATOR", "POSITIVE", "POTENTIATION",
             "OPENER"}
INH_TERMS = {"INHIBITOR", "INHIBITION", "IRREVERSIBLE INHIBITION", "FEEDBACK INHIBITION",
             "VOLTAGE-DEPENDENT INHIBITION", "ANTAGONIST", "ALLOSTERIC ANTAGONIST", "INVERSE AGONIST", "BLOCKER",
             "CHANNEL BLOCKER", "PORE BLOCKER", "GATING INHIBITOR", "NEGATIVE ALLOSTERIC MODULATOR",
             "NEGATIVE MODULATOR", "NEGATIVE", "COMPETITIVE", "NON-COMPETITIVE", "ANTISENSE INHIBITOR",
             "RNAI INHIBITOR", "DEGRADER", "DISRUPTING AGENT", "GENE EDITING NEGATIVE MODULATOR"}
# Not a modulation of a protein in the signaling sense: the drug is itself an enzyme, gene or antigen
DROP_TERMS = {"VACCINE ANTIGEN", "EXOGENOUS PROTEIN", "EXOGENOUS GENE", "HYDROLYTIC ENZYME", "PROTEOLYTIC ENZYME",
              "OXIDATIVE ENZYME"}

TIER = {"ChEMBL": 1, "DrugCentral": 1, "IUPHAR": 1, "ChemPert-DrugBank": 2, "ChemPert-Ramilowski": 2,
        "ChemPert-CLUE": 3, "ChemPert-STITCH": 4, "RepurposingHub": 3}


def direction(term):
    t = str(term).strip().upper()
    if t in DROP_TERMS:
        return None
    return "act" if t in ACT_TERMS else "inh" if t in INH_TERMS else "unknown"


def moa_direction(moa):
    """Repurposing Hub MoA text ('dopamine receptor agonist|adrenergic receptor antagonist') -> one direction,
    unknown when the MoAs disagree or name no direction."""
    dirs = set()
    for m in tokens_pipe(moa):
        m = m.lower()
        if re.search(r"inverse agonist|negative allosteric|inhibitor|antagonist|blocker|degrader", m):
            dirs.add("inh")
        elif re.search(r"agonist|activator|opener|positive allosteric|stimulant", m):
            dirs.add("act")
        else:
            dirs.add("unknown")
    return dirs.pop() if len(dirs) == 1 else "unknown"


def tokens_pipe(x):
    return [t.strip() for t in str(x).split("|") if t.strip()]


def norm_name(x):
    return re.sub(r"[^a-z0-9]", "", str(x).lower())


EV_COLUMNS = ["source", "target", "effect", "detail", "brd", "inchi_key", "names"]


def evidence(source, target, effect, detail="", brd="", inchi_key="", names=""):
    """Common loader output: one row per (source compound, raw target). names: '|'-separated."""
    df = pd.DataFrame(dict(source=source, target=target, effect=effect, detail=detail, brd=brd,
                           inchi_key=inchi_key, names=names))
    return df[df.effect.notna() & (df.target != "")].reset_index(drop=True)[EV_COLUMNS]


# --------------------------------------------------------------------------------------
# Loaders
# --------------------------------------------------------------------------------------
def load_lincs(path):
    df = read_table(path, index_col=0)  # written from R: first column holds the row names
    df.columns = ["pert_id", "pertname", "smiles", "inchi_key"]
    return df


def load_broad(drugs_path, samples_path):
    drugs = pd.read_csv(drugs_path, sep="\t", comment="!", dtype=str, keep_default_na=False)
    samples = pd.read_csv(samples_path, sep="\t", comment="!", dtype=str, keep_default_na=False)
    ids = samples.groupby("pert_iname").agg(brd=("broad_id", lambda s: "|".join(sorted({x[:13] for x in s}))),
                                            inchi_key=("InChIKey", lambda s: "|".join(sorted(set(s) - {""}))))
    drugs = drugs[drugs.target != ""].join(ids, on="pert_iname")
    drugs["effect"] = drugs.moa.map(moa_direction)
    drugs["target"] = drugs.target.map(tokens_pipe)
    drugs = drugs.explode("target")
    return evidence("RepurposingHub", drugs.target, drugs.effect, drugs.moa, drugs.brd.fillna(""),
                    drugs.inchi_key.fillna(""), drugs.pert_iname)


def load_chembl(chembl_dir, target_types):
    mech = pd.concat([pd.DataFrame(json.load(open(f))["mechanisms"])
                      for f in sorted(glob.glob(os.path.join(chembl_dir, "mechanism_*.json")))], ignore_index=True)
    mapping = pd.read_csv(os.path.join(chembl_dir, "chembl_uniprot_mapping.txt"), sep="\t", comment="#",
                          header=None, names=["uniprot", "target_chembl_id", "target_name", "target_type"])
    mapping = mapping[mapping.target_type.isin(target_types)]
    mols = read_table(os.path.join(chembl_dir, "molecules.tsv")).set_index("molecule_chembl_id")
    m = mech.merge(mapping, on="target_chembl_id").join(mols, on="parent_molecule_chembl_id")
    m["effect"] = m.action_type.fillna("").map(direction)
    names = (m.pref_name.fillna("") + "|" + m.synonyms.fillna("")).str.strip("|")
    return evidence("ChEMBL", m.uniprot, m.effect, m.action_type.fillna("") + " (" + m.target_type + ")",
                    inchi_key=m.inchi_key.fillna(""), names=names)


def load_drugcentral(path):
    """Only rows marked as the mechanism of action (MOA = 1); the others are bioactivity measurements."""
    dc = read_table(path, compression="gzip")
    dc = dc[(dc.MOA == "1") & (dc.ORGANISM == "Homo sapiens")].copy()
    dc["effect"] = dc.ACTION_TYPE.map(direction)
    dc["target"] = dc.GENE.map(tokens_pipe)
    dc = dc.explode("target")
    return evidence("DrugCentral", dc.target.fillna(""), dc.effect, dc.ACTION_TYPE, names=dc.DRUG_NAME)


def load_iuphar(path):
    """Human targets, endogenous ligands excluded. Direction from 'Action', falling back to 'Type'."""
    iu = pd.read_csv(path, skiprows=1, dtype=str, keep_default_na=False)
    iu = iu[(iu["Target Species"] == "Human") & (iu.Endogenous == "false")].copy()
    eff = iu.Action.map(direction)
    iu["effect"] = np.where(eff == "unknown", iu.Type.map(direction), eff)
    iu["target"] = iu["Target Gene Symbol"].map(tokens_pipe)
    iu = iu.explode("target")
    return evidence("IUPHAR", iu.target.fillna(""), iu.effect, iu.Type + "/" + iu.Action, names=iu.Ligand)


def load_chempert(path, keep_stitch_only):
    """ChemPert effect: 1 activation, -1 inhibition, 2 unknown. The source becomes ChemPert-<best resource>.
    'Chemical Compound' is a '; '-separated list of synonyms."""
    cp = read_table(path, index_col=0)
    res = cp.Source.map(lambda s: {t.strip() for t in s.split(";")})
    best = res.map(lambda r: next((x for x in ("DrugBank", "Ramilowski", "CLUE", "STITCH") if x in r), "STITCH"))
    cp["source"] = "ChemPert-" + best
    if not keep_stitch_only:
        cp = cp[best != "STITCH"]
    eff = cp.Effect.map({"1": "act", "-1": "inh", "2": "unknown"})
    names = cp["Chemical Compound"].str.replace(r"\s*;\s*", "|", regex=True)
    return pd.concat([evidence(src, sub.Target, eff[sub.index], sub.Source, names=names[sub.index])
                      for src, sub in cp.groupby("source")], ignore_index=True)


# --------------------------------------------------------------------------------------
# Matching
# --------------------------------------------------------------------------------------
class Matcher:
    """Looks up LINCS pert_ids by BRD id, InChIKey, InChIKey skeleton and normalised name (with aliases)."""

    def __init__(self, lincs):
        group = lambda keys: pd.Series(lincs.pert_id.to_numpy(), index=keys).groupby(level=0).agg(list).to_dict()
        has_key = lincs.inchi_key.str.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]")
        self.by_brd = group(lincs.pert_id.str[:13])
        self.by_inchi = group(lincs.inchi_key.where(has_key, ""))
        self.by_skeleton = group(lincs.inchi_key.str[:14].where(has_key, ""))
        self.by_name = {k: set(v) for k, v in group(lincs.pertname.map(norm_name)).items()}
        # compound identity, to tell several BRD ids of one compound from several compounds
        self.compound = dict(zip(lincs.pert_id, np.where(has_key, lincs.inchi_key.str[:14],
                                                         "name:" + lincs.pertname.map(norm_name))))
        for d in (self.by_inchi, self.by_skeleton, self.by_name):
            d.pop("", None)
        self.aliases = 0

    def add_aliases(self, ev):
        """Names of source compounds matched by id (BRD / InChIKey) become aliases of those LINCS perturbagens."""
        for perts, names in zip(ev.pert_id, ev.names):
            for n in tokens_pipe(names):
                n = norm_name(n)
                if len(n) >= 4:
                    before = len(self.by_name.get(n, ()))
                    self.by_name.setdefault(n, set()).update(perts)
                    self.aliases += len(self.by_name[n]) > before

    def match(self, ev):
        """Adds pert_id (list) and match (how) to source evidence rows; unmatched rows get an empty list."""
        perts, how = [], []
        for brd, ik, names in zip(ev.brd, ev.inchi_key, ev.names):
            for kind, keys, table in (("brd", tokens_pipe(brd), self.by_brd), ("inchikey", tokens_pipe(ik), self.by_inchi),
                                      ("skeleton", [k[:14] for k in tokens_pipe(ik)], self.by_skeleton),
                                      ("name", [norm_name(n) for n in tokens_pipe(names)], self.by_name)):
                hit = sorted({p for k in keys for p in table.get(k, ())})
                if hit and kind == "name" and len({self.compound[p] for p in hit}) > 1:
                    perts.append([])
                    how.append("ambiguous")
                    break
                if hit:
                    perts.append(hit)
                    how.append(kind)
                    break
            else:
                perts.append([])
                how.append("")
        return ev.assign(pert_id=perts, match=how)


# --------------------------------------------------------------------------------------
# Consensus
# --------------------------------------------------------------------------------------
def consensus(ev):
    """One row per (drug, target): direction from the best tier with a directional call."""
    ev = ev.assign(tier=ev.source.map(TIER))
    rows = []
    for (drug, target), g in ev.groupby(["drug", "target"], sort=False):
        directional = g[g.effect != "unknown"]
        effect, conflict = "unknown", False
        if len(directional):
            best = directional[directional.tier == directional.tier.min()]
            calls = set(best.effect)
            effect, conflict = (calls.pop(), False) if len(calls) == 1 else ("unknown", True)
        per_source = g.groupby("source").effect.agg(lambda s: "/".join(sorted(set(s))))
        rows.append((drug, target, effect, ";".join(sorted(g.source.unique())), g.source.nunique(),
                     ";".join(f"{s}:{e}" for s, e in per_source.items()), conflict))
    return pd.DataFrame(rows, columns=["drug", "target", "effect", "sources", "n_sources", "evidence", "conflict"])


def build(paths, chembl_target_types, keep_stitch_only=False, keep_weak_unknown=False):
    """Returns (drug_targets, per-source summary, matched evidence rows)."""
    hgnc = HGNC(paths["hgnc"])
    lincs = load_lincs(paths["lincs"])
    matcher = Matcher(lincs)
    loaded, summary = [], []

    def add(ev):
        mapped = ev.target.map(hgnc)
        ev = matcher.match(ev.assign(target=mapped))
        matched = ev[(ev.target != "") & (ev.pert_id.str.len() > 0)]
        for src, g in ev.groupby("source"):
            m = matched[matched.source == src]
            summary.append(dict(source=src, rows=len(g), target_in_hgnc=int((g.target != "").sum()),
                                rows_matched_to_lincs=len(m),
                                lincs_perts_matched=len({p for ps in m.pert_id for p in ps}),
                                **{f"by_{k}": int((m.match == k).sum()) for k in ("brd", "inchikey", "skeleton", "name")},
                                rejected_ambiguous_name=int((g.match == "ambiguous").sum())))
        loaded.append(matched)
        return matched

    # id-matched sources first, so their names become aliases for the name-only sources
    matcher.add_aliases(add(load_broad(paths["broad_drugs"], paths["broad_samples"])).query("match in ['brd', 'inchikey']"))
    matcher.add_aliases(add(load_chembl(paths["chembl_dir"], chembl_target_types)).query("match == 'inchikey'"))
    print(f"[match] {matcher.aliases} names added as aliases of LINCS perturbagens")
    add(load_drugcentral(paths["drugcentral"]))
    add(load_iuphar(paths["iuphar"]))
    add(load_chempert(paths["chempert"], keep_stitch_only))

    ev = pd.concat(loaded, ignore_index=True).explode("pert_id").rename(columns={"pert_id": "drug"})
    dt = consensus(ev)
    if not keep_weak_unknown:
        tier1 = dt.sources.map(lambda s: any(TIER[x] == 1 for x in s.split(";")))
        weak = (dt.effect == "unknown") & ~tier1 & (dt.n_sources < 2)
        print(f"[consensus] dropping {int(weak.sum()):,} unknown-direction pairs supported by one non-tier-1 source")
        dt = dt[~weak].reset_index(drop=True)
    dt.insert(3, "pertname", dt.drug.map(dict(zip(lincs.pert_id, lincs.pertname))))
    summary = pd.DataFrame(summary)
    return dt, summary, ev


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", required=True, help="metabognn_dev/data")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--chembl-target-types", default="SINGLE PROTEIN,PROTEIN COMPLEX,PROTEIN COMPLEX GROUP",
                    help="ChEMBL target types to expand into proteins (PROTEIN FAMILY targets are broad)")
    ap.add_argument("--keep-stitch-only", action="store_true", help="keep ChemPert entries backed by STITCH alone")
    ap.add_argument("--keep-weak-unknown", action="store_true",
                    help="keep unknown-direction pairs reported by a single non-tier-1 source")
    for name, default in DEFAULT_PATHS.items():
        ap.add_argument(f"--{name.replace('_', '-')}", default=default, help=f"relative to --data-dir ({default})")
    args = ap.parse_args()
    paths = {k: os.path.join(args.data_dir, getattr(args, k)) for k in DEFAULT_PATHS}

    dt, summary, _ = build(paths, args.chembl_target_types.split(","), args.keep_stitch_only, args.keep_weak_unknown)
    os.makedirs(args.out_dir, exist_ok=True)
    dt.to_csv(os.path.join(args.out_dir, "drug_targets.tsv"), sep="\t", index=False)
    summary.to_csv(os.path.join(args.out_dir, "source_summary.tsv"), sep="\t", index=False)
    print(summary.to_string(index=False))
    print(f"\n{dt.drug.nunique():,} LINCS perturbagens with >= 1 target, {len(dt):,} drug-target pairs, "
          f"{dt.target.nunique():,} targets; effect {dt.effect.value_counts().to_dict()}; "
          f"{int(dt.conflict.sum())} pairs with conflicting directions in their best tier")


if __name__ == "__main__":
    main()
