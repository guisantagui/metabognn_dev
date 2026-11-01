#!/usr/bin/env python
"""
Build the coupled metabolic + signaling network for the metabolic drug-target GNN (A549 prototype).

The metabolic signatures (MetChange, METAflux) are derived from expression through Human-GEM, so a drug can
only reach them by changing enzyme transcription:

    drug target -> signaling / modification -> TF -> enzyme gene -> reaction -> metabolite

Inputs (paths relative to --data-dir, i.e. metabognn_dev/data)
  GEMs/Human-GEM/model/Human-GEM.xml                     Human-GEM 2.0 (SBML)
  utilities/hgnc_complete_set.txt                        HGNC: every protein id is mapped to an approved symbol
  networks/omnipath/interactions.tsv                     OmniPath core (omnipath client, genesymbols=True)
  networks/omnipath/enzsub.tsv                           OmniPath enzyme-substrate
  networks/omnipath/collectri.tsv                        CollecTRI (optional; written by omnipath_dl.py)
  networks/truust/trrust_rawdata.human.tsv               TRRUST v2 human
  networks/ChEA/gene_attribute_edges.txt.gz              ChEA via Harmonizome (gene <- TF edges)
  networks/STRING/9606.protein.physical.links.detailed.v12.5.txt.gz  + 9606.protein.info.v12.5.txt.gz
  networks/HIPPIE/hippie_current.txt                     HIPPIE v2.4
  networks/KEGG/kgml/*.xml                               KEGG human pathway maps (dl_kegg.sh)
  GEMs/Human-GEM/model/metabolites.tsv                   KEGG compound ids of GEM metabolites

Relations (column `relation` in edges.tsv)
  metabolic         metabolite <-> reaction   GEM stoichiometry, currency/bookkeeping metabolites removed
  catalysis         enzyme gene -> reaction   GEM gene rules (AND/OR rule kept on the reaction node)
  grn               TF -> target gene         signed regulons: TRRUST, CollecTRI
  grn_binding       TF -> target gene         unsigned ChIP binding (ChEA); binding is not regulation
  signaling         protein -> protein        OmniPath rows with consensus direction (signed when known)
  enzyme_substrate  modifier -> substrate     OmniPath enz_sub; by default only substrates that are TFs
  ppi               protein - protein         undirected: STRING physical, HIPPIE, OmniPath rows without
                                              consensus direction; only pairs with no directed evidence
  second_messenger  protein <-> messenger     cAMP, cGMP, Ca2+, IP3, NO (cytosolic GEM node): KEGG protein-compound
                                              relations (producers -> messenger -> effectors; messenger ->
                                              effector edges only when drawn in >= kegg_min_maps maps), plus GEM
                                              enzymes/transporters that make (+) or remove (-) the messenger
                                              (e.g. PDE4D -| cAMP, ATP2A2 -| Ca2+)

Outputs (--out-dir)
  nodes.tsv          node_id, node_type, symbol, ensembl, base_met, compartment, gpr, gpr_ensembl, subsystem,
                     constitutive
  edges.tsv          src, dst, relation, directed, sign, n_resources, n_references, resources, origins, weight
  id_mapping.tsv     per input source: rows loaded, rows dropped because an endpoint has no HGNC symbol
  report.txt         edge counts, transcriptional-route coverage, suggested --alpha for the pipeline
  drug_metabolite_matches.tsv   (with --drug-targets) drugs whose name is a GEM metabolite

Usage
  python build_network.py --data-dir ../../data --out-dir ../../data/networks/combined
"""
import argparse
import os
import re
from collections import defaultdict
from dataclasses import dataclass, fields

import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse import csgraph

DEFAULT_PATHS = dict(
    gem="GEMs/Human-GEM/model/Human-GEM.xml",
    hgnc="utilities/hgnc_complete_set.txt",
    omnipath="networks/omnipath/interactions.tsv",
    enzsub="networks/omnipath/enzsub.tsv",
    collectri="networks/omnipath/collectri.tsv",
    trrust="networks/truust/trrust_rawdata.human.tsv",
    chea="networks/ChEA/gene_attribute_edges.txt.gz",
    string_links="networks/STRING/9606.protein.physical.links.detailed.v12.5.txt.gz",
    string_info="networks/STRING/9606.protein.info.v12.5.txt.gz",
    hippie="networks/HIPPIE/hippie_current.txt",
    kegg_kgml="networks/KEGG/kgml",
    gem_metabolites="GEMs/Human-GEM/model/metabolites.tsv",
)

DIRECTED_PROTEIN_RELATIONS = ["signaling", "enzyme_substrate", "grn", "grn_binding"]
PROTEIN_RELATIONS = DIRECTED_PROTEIN_RELATIONS + ["ppi"]
ROUTE_RELATIONS = ["signaling", "enzyme_substrate", "grn", "second_messenger"]  # directed, visible in expression

# Starting mixing weights for metabolic_target_prediction.py --alpha with an expression-derived signal.
# Not tuned: they encode which routes can be seen in the signature. Tune them by ablation.
SUGGESTED_ALPHA = {"metabolic": 1.0, "catalysis": 1.0, "grn": 1.0, "signaling": 1.0, "enzyme_substrate": 1.0,
                   "second_messenger": 1.0, "grn_binding": 0.25, "ppi": 0.25}

# KEGG compound ids of the second messengers whose protein relations are loaded (DAG: C00165 / C00641)
MESSENGERS = {"C00575": "cAMP", "C00942": "cGMP", "C00076": "Ca2+", "C01245": "IP3", "C00165": "DAG",
              "C00641": "DAG", "C00533": "NO"}


@dataclass
class Config:
    min_directed_resources: int = 2   # signaling / grn: keep if >= this many resources ...
    min_directed_refs: int = 1        # ... or >= this many PMIDs
    min_ppi_resources: int = 2
    min_ppi_refs: int = 1
    enzsub_min_primary: int = 2
    enzsub_min_refs: int = 2
    enzsub_substrates: str = "tf"     # all | tf | gem | tf+gem
    string_channels: str = "experiments,database"
    string_min_score: int = 700
    hippie_min_score: float = 0.73    # HIPPIE "high confidence"
    max_ppi_degree: int = 0           # 0 = keep hubs
    prune: bool = True                # drop components without metabolites
    max_hops: int = 5
    kegg_min_maps: int = 2            # KEGG messenger -> effector edges drawn in fewer maps are dropped


# --------------------------------------------------------------------------------------
# GEM exclusions (from build_gem_layer.ipynb)
# --------------------------------------------------------------------------------------
# Human-GEM names: 'chloride' and 'hydroxide' (the notebook's 'Cl-' / 'OH-' match no metabolite)
CURRENCY = {
    "H+", "H2O", "hydroxide", "CO2", "HCO3-", "Na+", "K+", "chloride", "iodide", "O2", "NH3", "NH4+", "Pi", "PPi",
    "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "CTP", "CDP", "CMP", "UTP", "UDP", "UMP",
    "NAD+", "NADH", "NADP+", "NADPH", "FAD", "FADH2", "FMN", "CoA", "H2O2", "SAM", "SAH", "THF",
    "5,10-methylene-THF", "ubiquinone", "ubiquinol", "ferricytochrome C", "ferrocytochrome C",
    "thioredoxin", "mitothioredoxin", "oxidized thioredoxin", "O2-",
    "L-carnitine",  # acyl carrier of the carnitine shuttle, like CoA (~1,000 edges across 5 compartments)
}

# Reactions containing any of these are dropped entirely. thioredoxin / mitothioredoxin were also listed here
# in the notebook; they are electron carriers (already in CURRENCY), and as bookkeeping they removed whole
# reactions such as ribonucleotide reductase, so they are left out.
BOOKKEEPING = {
    "biomass", "lipid_pool_biomass", "cofactor_pool_biomass", "protein_pool_biomass", "metabolite_pool_biomass",
    "DNA", "RNA", "[ACP]", "[protein]", "[protein]-L-arginine", "[protein]-L-asparagine", "[protein]-L-citrulline",
    "[protein]-L-cysteine", "[protein]-L-glutamine", "[protein]-L-isoaspartate",
    "[protein]-L-isoaspartate-methyl ester", "[protein]-L-lysine", "[protein]-L-serine", "[protein]-L-tyrosine",
    "[protein]-N5-alkylglutamine", "[protein]-N6-(dihydrolipoyl)lysine", "[protein]-N6-(lipoyl)lysine",
    "[protein]-N6-(octanoyl)lysine", "[protein]-N6,N6,N6-trimethyl-L-lysine", "[protein]-N6,N6-dimethyl-L-lysine",
    "[protein]-N6-methyl-L-lysine", "[protein]-N-ubiquityllysine", "[protein]-omega-N-(ADP-D-ribosyl)-L-arginine",
    "[protein]-S-methyl-L-cysteine", "[protein]-tyrosine-O-sulfate", "Gal2-Xyl-L-Ser-[protein]",
    "Gal-Xyl-L-Ser-[protein]", "Xyl-L-Ser-[protein]", "albumin", "haptoglobin", "antichymotrypsin", "apoA1",
    "apoB100", "apoC1", "apoC2", "apoC3", "glycogenin", "apoE", "plasminogen", "apo-[ACP]", "prothrombin",
    "fibrinogen", "antitrypsin", "[apotransferin]",
}

MACROMOLECULAR = {
    "[protein]", "[protein]-L-serine", "[ACP]", "DNA", "DNA-5-methylcytosine", "RNA",
    "glycyl-tRNA(gly)", "L-alanyl-tRNA(ala)", "L-arginyl-tRNA(arg)", "L-asparaginyl-tRNA(asn)",
    "L-aspartyl-tRNA(asp)", "L-cysteinyl-tRNA(cys)", "L-glutaminyl-tRNA(gln)", "L-glutamyl-tRNA(glu)",
    "L-histidyl-tRNA(his)", "L-isoleucyl-tRNA(ile)", "L-leucyl-tRNA(leu)", "L-lysyl-tRNA(lys)",
    "L-methionyl-tRNA(met)", "L-phenylalanyl-tRNA(phe)", "L-prolyl-tRNA(pro)", "L-seryl-tRNA(ser)",
    "L-threonyl-tRNA(thr)", "L-tryptophanyl-tRNA(trp)", "L-tyrosyl-tRNA(tyr)", "L-valyl-tRNA(val)",
    "N-formylmethionyl-tRNA", "N-substituted aminoacyl tRNA", "selenomethionyl-tRNA(met)",
    "tRNA containing 5-aminomethyl-2-thiouridine", "tRNA containing 5-methylaminomethyl-2-thiouridylate",
    "tRNA containing 6-isopentenyladenosine", "tRNA containing N2-methylguanine", "tRNA containing N7-methylguanine",
    "tRNA(ala)", "tRNA(arg)", "tRNA(asn)", "tRNA(asp)", "tRNA(cys)", "tRNA(gln)", "tRNA(glu)", "tRNA(gly)",
    "tRNA(his)", "tRNA(ile)", "tRNA(leu)", "tRNA(lys)", "tRNA(met)", "tRNA(phe)", "tRNA(pro)", "tRNA(ser)",
    "tRNA(thr)", "tRNA(trp)", "tRNA(tyr)", "tRNA(val)", "tRNA", "tRNA-guanine", "tRNA-pseudouridine",
    "tRNA-queuine", "tRNA-uridine", "prothrombin",
}


def gem_exclusions(gem):
    """Metabolite names removed as nodes: currency, bookkeeping, macromolecules and lipid 'pool' pseudo-metabolites."""
    pool = {m.name for m in gem.metabolites if "pool" in m.name}
    return CURRENCY | BOOKKEEPING | MACROMOLECULAR | pool, BOOKKEEPING


# --------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------
EDGE_COLUMNS = ["src", "dst", "relation", "directed", "sign", "resources", "pmids", "n_resources",
                "n_references", "origins"]


def edge_frame(origin, **cols):
    """Common loader output. n_resources / n_references are lower bounds when names or PMIDs are not listed."""
    df = pd.DataFrame(cols)
    for c, default in (("sign", 0), ("resources", ""), ("pmids", ""), ("n_resources", 1), ("n_references", 0)):
        if c not in df:
            df[c] = default
    df["origins"] = origin
    df["src"], df["dst"] = df.src.astype(str), df.dst.astype(str)
    df["directed"] = df.directed.astype(int)
    df["sign"] = df.sign.astype(int)
    return df[EDGE_COLUMNS]


def tokens(x, sep=r"[;,|]"):
    if not isinstance(x, str) or not x:
        return set()
    return {t.strip() for t in re.split(sep, x) if t.strip()}


def pmid_set(x):
    """'KEA:12783885;SIGNOR:111' or '12783885;111' -> {'12783885', '111'}"""
    return {t.split(":")[-1] for t in tokens(x)} - {""}


def primary_resources(x):
    """Collapse derived OmniPath resources (PhosphoSite_ProtMapper, phosphoELM_MIMP) onto their primary."""
    return {t.split("_")[0] for t in tokens(x)}


def joined(sets):
    return [";".join(sorted(s)) for s in sets]


def read_table(path, **kw):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **kw)


# --------------------------------------------------------------------------------------
# Identifiers
# --------------------------------------------------------------------------------------
class HGNC:
    """Maps any protein identifier to its HGNC approved symbol: approved symbols map to themselves; Entrez,
    UniProt, Ensembl, previous symbols and aliases map when unambiguous (an alias never overrides an approved
    symbol). Only protein-coding genes are used, plus complex-locus constituents (HGNC files the UGT1A enzymes
    there)."""

    def __init__(self, path):
        h = read_table(path)
        h = h[(h.locus_group == "protein-coding gene") | (h.locus_type == "complex locus constituent")]
        self.approved = set(h.symbol)
        self.ensembl = dict(zip(h.symbol, h.ensembl_gene_id))
        self.map = {s: s for s in self.approved}
        self.ambiguous = {}
        for col in ("entrez_id", "uniprot_ids", "ensembl_gene_id", "prev_symbol", "alias_symbol"):
            hits = defaultdict(set)
            for sym, val in zip(h.symbol, h[col]):
                for v in tokens(val, r"\|"):
                    hits[v].add(sym)
            for v, syms in hits.items():
                if v in self.map:
                    continue
                if len(syms) == 1:
                    self.map[v] = next(iter(syms))
                else:
                    self.ambiguous[v] = sorted(syms)

    def __call__(self, x):
        x = str(x).strip()
        return self.map.get(x) or self.map.get(x.upper()) or ""


def harmonise(raw, hgnc, gem_ids):
    """Map protein endpoints to HGNC symbols (GEM metabolites/reactions untouched) and drop edges with an
    unmappable endpoint. Returns the mapped edges, a per-origin summary and the unmapped ids with counts."""
    out = raw.copy()
    lost = np.zeros(len(out), bool)
    for col in ("src", "dst"):
        prot = ~out[col].isin(gem_ids).to_numpy()
        mapped = out.loc[prot, col].map(hgnc)
        lost[prot] |= (mapped == "").to_numpy()
        out.loc[prot, col] = mapped
    summary = (pd.DataFrame(dict(origin=out.origins, lost=lost)).groupby("origin")
               .agg(rows=("lost", "size"), dropped_unmapped=("lost", "sum")))
    summary["dropped_pct"] = (100 * summary.dropped_unmapped / summary.rows).round(2)
    bad = np.concatenate([raw.src.to_numpy(dtype=object)[lost & (out.src == "").to_numpy()],
                          raw.dst.to_numpy(dtype=object)[lost & (out.dst == "").to_numpy()]])
    return out[~lost].reset_index(drop=True), summary.reset_index(), pd.Series(bad).value_counts()


# --------------------------------------------------------------------------------------
# GEM layer
# --------------------------------------------------------------------------------------
def load_gem(path):
    import cobra
    return cobra.io.read_sbml_model(path)


def build_bipartite_graph(gem, excluded_mets, bookkeeping):
    """Metabolite <-> reaction graph. Excluded metabolites are matched on their full name (the notebook used
    base_id(met.name), which cut names containing '_', e.g. 'betaine_aldehyde' -> 'betaine')."""
    G = nx.DiGraph()
    for rxn in gem.reactions:
        if any(m.name in bookkeeping for m in rxn.metabolites):
            continue  # whole reaction dropped
        G.add_node(rxn.id, node_type="reaction", constitutive=len(rxn.genes) == 0)
        reversible = rxn.lower_bound < 0 and rxn.upper_bound > 0
        for met in rxn.reactants:
            if met.name in excluded_mets:
                continue
            G.add_node(met.id, node_type="metabolite")
            G.add_edge(met.id, rxn.id, edge_type="consumes")
            if reversible:
                G.add_edge(rxn.id, met.id, edge_type="produces_rev")
        for met in rxn.products:
            if met.name in excluded_mets:
                continue
            G.add_node(met.id, node_type="metabolite")
            G.add_edge(rxn.id, met.id, edge_type="produces")
            if reversible:
                G.add_edge(met.id, rxn.id, edge_type="consumes_rev")
    return G


def met_base_and_compartment(met):
    """Compartment-free id: Human-GEM 'MAM01371c' -> 'MAM01371'."""
    comp = met.compartment or ""
    if comp and met.id.endswith(comp):
        return met.id[: -len(comp)], comp
    return met.id, comp


def gem_gene_symbols(gem, hgnc):
    """Ensembl gene id -> HGNC symbol; falls back to the model's hgnc.symbol annotation."""
    sym = {}
    for g in gem.genes:
        ann = g.annotation.get("hgnc.symbol", "")
        sym[g.id] = hgnc(g.id) or hgnc(ann) or ann or g.id
    return sym


def rule_to_symbols(rule, gene_sym):
    if not rule:
        return ""
    return re.sub(r"[^\s()]+",
                  lambda m: m.group(0) if m.group(0).lower() in ("and", "or") else gene_sym.get(m.group(0), m.group(0)),
                  rule)


def build_gem_layer(gem, excluded_mets, bookkeeping, gene_sym):
    """GEM nodes (reactions, metabolites) and edges (metabolic, catalysis). Catalysis sources are the symbols from
    gem_gene_symbols (Ensembl -> HGNC, falling back to the model's symbol when HGNC lists another Ensembl id)."""
    G = build_bipartite_graph(gem, excluded_mets, bookkeeping)
    nodes = []
    for n, data in G.nodes(data=True):
        if data["node_type"] == "reaction":
            rxn = gem.reactions.get_by_id(n)
            nodes.append(dict(node_id=n, node_type="reaction", symbol="", ensembl="", base_met="", compartment="",
                              gpr=rule_to_symbols(rxn.gene_reaction_rule, gene_sym),
                              gpr_ensembl=rxn.gene_reaction_rule, subsystem=str(rxn.subsystem or ""),
                              constitutive=int(data["constitutive"])))
        else:
            met = gem.metabolites.get_by_id(n)
            b, c = met_base_and_compartment(met)
            nodes.append(dict(node_id=n, node_type="metabolite", symbol=met.name, ensembl="", base_met=b,
                              compartment=c, gpr="", gpr_ensembl="", subsystem="", constitutive=0))
    nodes = pd.DataFrame(nodes)
    s, d = zip(*G.edges())
    metabolic = edge_frame("GEM", src=list(s), dst=list(d), relation="metabolic", directed=1, resources="GEM")
    cat = [(gene_sym[g.id], r) for r in nodes.node_id[nodes.node_type == "reaction"]
           for g in gem.reactions.get_by_id(r).genes]
    catalysis = edge_frame("GEM", src=[c[0] for c in cat], dst=[c[1] for c in cat], relation="catalysis",
                           directed=1, sign=1, resources="GEM")
    return nodes, pd.concat([metabolic, catalysis], ignore_index=True), G


# --------------------------------------------------------------------------------------
# Protein-layer loaders (raw identifiers; harmonise() maps them to HGNC symbols)
# --------------------------------------------------------------------------------------
def _flag(df, col):
    return df[col].isin(["True", "true", "1"]).astype(int) if col in df else pd.Series(0, index=df.index)


def expand_complexes(df):
    """OmniPath complexes ('COMPLEX:P1_P2', genesymbol 'A_B') -> one row per subunit."""
    for side in ("source", "target"):
        cplx = df[side].str.startswith("COMPLEX:").to_numpy()
        sym = df[f"{side}_genesymbol"].tolist()
        df = df.assign(**{f"{side}_genesymbol": [s.split("_") if c else [s] for s, c in zip(sym, cplx)]})
        df = df.explode(f"{side}_genesymbol", ignore_index=True)
    return df


def _omnipath_rows(df, origin, relation):
    """Shared by OmniPath interactions and CollecTRI. Rows with consensus direction are directed edges, signed by
    the consensus; for interactions, rows without it go to ppi."""
    n_complex = int((df.source.str.startswith("COMPLEX:") | df.target.str.startswith("COMPLEX:")).sum())
    n_rows = len(df)
    df = expand_complexes(df)
    directed = _flag(df, "consensus_direction").to_numpy()
    rel = np.where(directed == 1, relation, "ppi") if relation == "signaling" else np.full(len(df), relation)
    sign = (_flag(df, "consensus_stimulation") - _flag(df, "consensus_inhibition")).to_numpy()
    out = edge_frame(origin, src=df.source_genesymbol, dst=df.target_genesymbol, relation=rel,
                     directed=np.where(rel == "ppi", 0, 1), sign=np.where(rel == "ppi", 0, sign),
                     resources=joined(df.sources.map(primary_resources)),
                     pmids=joined(df.references_stripped.map(pmid_set)),
                     n_resources=pd.to_numeric(df.n_primary_sources, errors="coerce").fillna(1).astype(int))
    print(f"[{origin}] {n_rows} rows ({n_complex} with complexes) -> {len(out)} rows after expanding complexes: "
          f"{out.relation.value_counts().to_dict()}")
    return out


def load_omnipath(path):
    return _omnipath_rows(read_table(path, index_col=0), "OmniPath", "signaling")


def load_collectri(path):
    """CollecTRI lists its original resources as e.g. 'TRRUST_CollecTRI' (collapsed to 'TRRUST', so they merge
    with the TRRUST file) plus the aggregator itself; the aggregator token is dropped so it does not add one
    resource to every edge. TF complexes (AP-1, NF-kB) are expanded into subunits."""
    out = _omnipath_rows(read_table(path, index_col=0), "CollecTRI", "grn")
    out["resources"] = [";".join(sorted(tokens(r) - {"CollecTRI"})) or "CollecTRI" for r in out.resources]
    return out


def load_enzsub(path):
    """One row per modified site -> one directed, unsigned edge per enzyme/substrate pair (a modification can
    activate or inhibit its substrate)."""
    df = read_table(path, index_col=0)
    df["pm"] = df.references_stripped.map(pmid_set)
    df["res"] = df.sources.map(primary_resources)
    g = df.groupby(["enzyme_genesymbol", "substrate_genesymbol"])
    pairs = pd.DataFrame(dict(pm=g.pm.agg(lambda s: set().union(*s)), res=g.res.agg(lambda s: set().union(*s)),
                              mods=g.modification.agg(lambda s: ";".join(sorted(set(s)))))).reset_index()
    print(f"[enz_sub] {len(df)} site rows -> {len(pairs)} enzyme/substrate pairs")
    return edge_frame("enz_sub", src=pairs.enzyme_genesymbol, dst=pairs.substrate_genesymbol,
                      relation="enzyme_substrate", directed=1, resources=joined(pairs.res), pmids=joined(pairs.pm))


def load_trrust(path):
    df = read_table(path, header=None, names=["tf", "target", "mode", "pmids"])
    sign = df["mode"].map({"Activation": 1, "Repression": -1}).fillna(0).astype(int)
    print(f"[TRRUST] {len(df)} TF-target rows, modes {df['mode'].value_counts().to_dict()}")
    return edge_frame("TRRUST", src=df.tf, dst=df.target, relation="grn", directed=1, sign=sign,
                      resources="TRRUST", pmids=joined(df.pmids.map(pmid_set)))


def load_chea(path):
    """Harmonizome ChEA gene_attribute_edges: column 'source' is the target gene and 'target' is the TF
    (attribute); the second line repeats the column descriptions. Entrez ids are used for mapping. Harmonizome
    collapses experiments, so there is no experiment count or PMID per edge."""
    df = read_table(path, compression="gzip").iloc[1:]
    tf = np.where(df.target_id.str.fullmatch(r"\d+"), df.target_id, df.target)
    gene = np.where(df.source_id.str.fullmatch(r"\d+"), df.source_id, df.source)
    print(f"[ChEA] {len(df)} TF-gene rows, {df.target.nunique()} TFs, {df.source.nunique()} genes")
    return edge_frame("ChEA", src=tf, dst=gene, relation="grn_binding", directed=1, resources="ChEA")


STRING_CHANNELS = {"experiments": ("experiments", "experimental"), "database": ("database",),
                   "textmining": ("textmining",), "coexpression": ("coexpression",)}
STRING_PRIOR = 0.041


def string_score(channel_scores):
    """STRING's combination: remove the prior per channel, combine as 1 - prod(1 - s), add the prior back."""
    s_nop = np.clip((channel_scores - STRING_PRIOR) / (1 - STRING_PRIOR), 0, None)
    comb = 1 - np.prod(1 - s_nop, 1)
    return comb + STRING_PRIOR * (1 - comb), s_nop > 0


def load_string(links_path, info_path, channels, min_score):
    """STRING physical links as undirected PPI, rescored from the selected channels only (text mining does not
    leak in through combined_score). Each contributing channel counts as one resource."""
    info = read_table(info_path, usecols=[0, 1])
    ensp2sym = dict(zip(info.iloc[:, 0], info.iloc[:, 1]))
    df = pd.read_csv(links_path, sep=" ")
    df = df[df.protein1 < df.protein2]  # each pair is listed in both orientations
    cols = {ch: c for ch in channels for c in STRING_CHANNELS.get(ch, (ch,)) if c in df.columns}
    score, support = string_score(np.stack([df[c].to_numpy() / 1000.0 for c in cols.values()], 1))
    keep = score >= min_score / 1000.0
    names = np.array([f"STRING-{ch}" for ch in cols])
    res = [";".join(names[row]) for row in support[keep]]
    print(f"[STRING] {len(df)} pairs, {keep.sum()} with score >= {min_score} on {list(cols)}")
    return edge_frame("STRING", src=df.protein1.to_numpy()[keep], dst=df.protein2.to_numpy()[keep], relation="ppi",
                      directed=0, resources=res, n_resources=[r.count(";") + 1 for r in res]).assign(
        src=lambda d: d.src.map(ensp2sym).fillna(""), dst=lambda d: d.dst.map(ensp2sym).fillna(""))


def load_hippie(path, min_score):
    """HIPPIE: uniprot_name1, entrez1, uniprot_name2, entrez2, score, info. Entrez ids are mapped by HGNC; the
    databases listed in info (BioGRID, IntAct, ...) are kept as resources so they merge with OmniPath's."""
    df = read_table(path, header=None, names=["u1", "e1", "u2", "e2", "score", "info"])
    df = df[pd.to_numeric(df.score, errors="coerce") >= min_score]
    field = lambda key: df["info"].str.extract(key + r":([^;]*)", expand=False).fillna("")
    print(f"[HIPPIE] {len(df)} pairs with score >= {min_score}")
    return edge_frame("HIPPIE", src=df.e1, dst=df.e2, relation="ppi", directed=0,
                      resources=joined(field("sources").map(tokens)), pmids=joined(field("pmids").map(pmid_set)))


def messenger_nodes(metabolites_path, gem_ids):
    """KEGG compound id -> cytosolic GEM metabolite node present in the network (signaling happens there)."""
    m = read_table(metabolites_path)
    out = {}
    for met, keggs in zip(m.mets, m.metKEGGID):
        for k in tokens(keggs):
            if k in MESSENGERS and met.endswith("c") and met in gem_ids:
                out.setdefault(k, met)
    return out


def load_kegg_messengers(kgml_dir, node_of):
    """Protein <-> compound relations (KGML 'PCrel') through a second messenger, from every human map:
    producer -> messenger (e.g. ADCY -> cAMP) and messenger -> effector (cAMP -> PRKACA). Gene entries list
    Entrez ids ('hsa:107 hsa:108'), groups are expanded into their component genes. Sign from the subtype
    (activation +1, inhibition -1). KGML has no PMIDs; each map in which an edge is drawn counts as a resource."""
    import xml.etree.ElementTree as ET
    rows, maps, skipped = [], 0, defaultdict(int)
    for f in sorted(os.listdir(kgml_dir)):
        path = os.path.join(kgml_dir, f)
        if not f.endswith(".xml") or os.path.getsize(path) == 0:
            continue
        try:
            root = ET.parse(path).getroot()
        except ET.ParseError:
            continue
        maps += 1
        ent = {e.get("id"): e for e in root.findall("entry")}

        def genes(i):
            e = ent.get(i)
            if e is None:
                return []
            if e.get("type") == "group":
                return [g for c in e.findall("component") for g in genes(c.get("id"))]
            return [x[4:] for x in e.get("name", "").split() if x.startswith("hsa:")] if e.get("type") == "gene" else []

        def compound(i):
            e = ent.get(i)
            if e is None or e.get("type") != "compound":
                return None
            return next((x[4:] for x in e.get("name", "").split() if x[4:] in MESSENGERS), None)

        for rel in root.findall("relation"):
            if rel.get("type") != "PCrel":
                continue
            subtypes = {st.get("name") for st in rel.findall("subtype")}
            sign = 1 if "activation" in subtypes else -1 if "inhibition" in subtypes else 0
            a, b = rel.get("entry1"), rel.get("entry2")
            for cpd, prots, outgoing in ((compound(b), genes(a), False), (compound(a), genes(b), True)):
                if cpd is None:
                    continue
                if cpd not in node_of:
                    skipped[MESSENGERS[cpd]] += len(prots)
                    continue
                for g in prots:
                    src, dst = (node_of[cpd], g) if outgoing else (g, node_of[cpd])
                    rows.append((src, dst, sign, f"KEGG:{root.get('name', f)[5:]}"))
    df = pd.DataFrame(rows, columns=["src", "dst", "sign", "map"])
    print(f"[KEGG] {maps} maps, {len(df)} protein-messenger relations drawn "
          f"({len(df.drop_duplicates(['src', 'dst']))} distinct); skipped (messenger not in network): {dict(skipped)}")
    return edge_frame("KEGG", src=df.src, dst=df.dst, relation="second_messenger", directed=1, sign=df.sign,
                      resources=df["map"])


def gem_messenger_edges(gem, G, node_of, gene_sym):
    """Enzymes and transporters that make (+1) or remove (-1) a cytosolic messenger in a kept GEM reaction, e.g.
    PDE4D -| cAMP (cAMP -> AMP) or ATP2A2 -| Ca2+ (cytosol -> ER); reversible reactions are unsigned. The
    catalysis and metabolic edges connect these enzymes already, but only along the reaction, so a degrader had
    no directed path to the messenger."""
    messengers = set(node_of.values())
    rows = []
    for rxn in gem.reactions:
        if rxn.id not in G:
            continue
        reversible = rxn.lower_bound < 0 < rxn.upper_bound
        for met, coef in rxn.metabolites.items():
            if met.id in messengers:
                sign = 0 if reversible else (1 if coef > 0 else -1)
                rows += [(gene_sym[g.id], met.id, sign) for g in rxn.genes]
    df = pd.DataFrame(rows, columns=["src", "dst", "sign"]).drop_duplicates()
    print(f"[GEM messengers] {len(df)} enzyme/transporter -> messenger edges, sign {df.sign.value_counts().to_dict()}")
    return edge_frame("GEM", src=df.src, dst=df.dst, relation="second_messenger", directed=1, sign=df.sign,
                      resources="GEM")


def load_sources(paths, cfg, messenger_node_of=None):
    raw = [load_omnipath(paths["omnipath"]), load_enzsub(paths["enzsub"]), load_trrust(paths["trrust"]),
           load_chea(paths["chea"]),
           load_string(paths["string_links"], paths["string_info"], cfg.string_channels.split(","),
                       cfg.string_min_score),
           load_hippie(paths["hippie"], cfg.hippie_min_score)]
    if os.path.exists(paths["collectri"]):
        raw.append(load_collectri(paths["collectri"]))
    else:
        print(f"[CollecTRI] {paths['collectri']} not found: signed regulons come from TRRUST only")
    if os.path.isdir(paths["kegg_kgml"]) and messenger_node_of:
        raw.append(load_kegg_messengers(paths["kegg_kgml"], messenger_node_of))
    else:
        print(f"[KEGG] {paths['kegg_kgml']} not found: no second-messenger edges (run dl_kegg.sh)")
    return pd.concat(raw, ignore_index=True)


# --------------------------------------------------------------------------------------
# Merging and curation
# --------------------------------------------------------------------------------------
KEYS = ["relation", "src", "dst", "directed"]


def merge_edges(raw):
    """Collapse duplicates within a relation (undirected pairs stored once, sorted). Resources, PMIDs and
    origins are unioned; counts are the larger of the union size and the per-row lower bounds. Opposite signs
    reported for the same edge give sign 0 (unknown)."""
    raw = raw[raw.src != raw.dst].copy()
    und = (raw.directed == 0).to_numpy()
    a, b = raw.src.to_numpy(dtype=object), raw.dst.to_numpy(dtype=object)
    swap = und & (a > b)
    raw["src"], raw["dst"] = np.where(swap, b, a), np.where(swap, a, b)
    dup = raw.duplicated(KEYS, keep=False)
    union = lambda s: ";".join(sorted(set().union(*(tokens(x) for x in s))))
    m = (raw[dup].groupby(KEYS, sort=False)
         .agg(pos=("sign", lambda s: (s > 0).any()), neg=("sign", lambda s: (s < 0).any()),
              resources=("resources", union), pmids=("pmids", union), origins=("origins", union),
              n_resources=("n_resources", "max"), n_references=("n_references", "max")).reset_index())
    m["sign"] = np.where(m.pos & ~m.neg, 1, np.where(m.neg & ~m.pos, -1, 0))
    merged = pd.concat([raw[~dup], m.drop(columns=["pos", "neg"])], ignore_index=True)
    n_res = merged.resources.map(lambda x: len({t.lower() for t in tokens(x)}))
    merged["n_resources"] = np.maximum(merged.n_resources.astype(int), n_res)
    merged["n_references"] = np.maximum(merged.n_references.astype(int), merged.pmids.map(lambda x: len(tokens(x))))
    return merged[EDGE_COLUMNS]


def curate(e, cfg, messengers=()):
    """Evidence thresholds per relation, enz_sub substrate filter, PPI redundancy with directed edges.
    messengers: GEM node ids of the second messengers (their outgoing edges are the effector edges)."""
    rel = e.relation
    passes = lambda min_res, min_ref: (e.n_resources >= min_res) | (e.n_references >= min_ref)
    keep = ~rel.isin(PROTEIN_RELATIONS + ["second_messenger"])  # metabolic, catalysis
    keep |= rel.isin(["signaling", "grn"]) & passes(cfg.min_directed_resources, cfg.min_directed_refs)
    keep |= (rel == "enzyme_substrate") & passes(cfg.enzsub_min_primary, cfg.enzsub_min_refs)
    keep |= rel == "grn_binding"  # Harmonizome ChEA has no per-edge evidence to filter on
    # KEGG maps carry no PMIDs: evidence is the number of maps drawing the edge. Only messenger -> effector edges
    # are filtered (a messenger broadcasts to every effector; many producers into one messenger are averaged).
    effector = e.src.isin(set(messengers))
    keep |= (rel == "second_messenger") & (~effector | (e.n_resources >= cfg.kegg_min_maps))
    keep |= (rel == "ppi") & passes(cfg.min_ppi_resources, cfg.min_ppi_refs)
    log = {"dropped_by_evidence": e[~keep].relation.value_counts().to_dict()}
    e = e[keep]

    if cfg.enzsub_substrates != "all":
        allowed = set()
        if "tf" in cfg.enzsub_substrates:
            allowed |= set(e.src[e.relation.isin(["grn", "grn_binding"])])
        if "gem" in cfg.enzsub_substrates:
            allowed |= set(e.src[e.relation == "catalysis"])
        drop = (e.relation == "enzyme_substrate") & ~e.dst.isin(allowed)
        log["enzsub_dropped_by_substrate_filter"] = int(drop.sum())
        e = e[~drop]

    # PPI only where no directed protein evidence exists in either orientation
    dir_pairs = set(map(frozenset, e.loc[e.relation.isin(DIRECTED_PROTEIN_RELATIONS), ["src", "dst"]].to_numpy()))
    ppi = (e.relation == "ppi").to_numpy()
    redundant = ppi & np.array([frozenset(p) in dir_pairs for p in e[["src", "dst"]].to_numpy()])
    log["ppi_dropped_redundant_with_directed"] = int(redundant.sum())
    e = e[~redundant]
    if cfg.max_ppi_degree:
        p = e[e.relation == "ppi"]
        deg = pd.concat([p.src, p.dst]).value_counts()
        hubs = set(deg[deg > cfg.max_ppi_degree].index)
        e = e[~((e.relation == "ppi") & (e.src.isin(hubs) | e.dst.isin(hubs)))]
        log["ppi_hubs_removed"] = len(hubs)
    e = e.copy()
    e["weight"] = evidence_weight(e)
    return e.reset_index(drop=True), log


def evidence_weight(e):
    # 1 resource -> 0.5, 2 -> 0.75, 3 -> 0.875 ... ; GEM-derived relations weight 1
    w = 1 - 0.5 ** np.maximum(e.n_resources.astype(float), 1)
    return np.where(e.relation.isin(["metabolic", "catalysis"]), 1.0, w)


def assemble_nodes(nodes_gem, edges, hgnc):
    gem_ids = set(nodes_gem.node_id)
    proteins = sorted((set(edges.src) | set(edges.dst)) - gem_ids)
    prot = pd.DataFrame(dict(node_id=proteins, node_type="protein", symbol=proteins,
                             ensembl=[hgnc.ensembl.get(p, "") for p in proteins], base_met="", compartment="",
                             gpr="", gpr_ensembl="", subsystem="", constitutive=0))
    return pd.concat([nodes_gem, prot], ignore_index=True)


def adjacency(nodes, edges, relations=None):
    """Sparse adjacency with undirected edges in both directions: A[i, j] != 0 means i -> j."""
    if relations is not None:
        edges = edges[edges.relation.isin(relations)]
    idx = pd.Series(np.arange(len(nodes)), index=nodes.node_id)
    s, d = idx[edges.src].to_numpy(), idx[edges.dst].to_numpy()
    und = edges.directed.to_numpy() == 0
    rows, cols = np.concatenate([s, d[und]]), np.concatenate([d, s[und]])
    return sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(nodes),) * 2)


def keep_metabolic_components(nodes, edges):
    _, comp = csgraph.connected_components(adjacency(nodes, edges), directed=True, connection="weak")
    keep = set(nodes.node_id[np.isin(comp, list(set(comp[(nodes.node_type == "metabolite").to_numpy()])))])
    return (nodes[nodes.node_id.isin(keep)].reset_index(drop=True),
            edges[edges.src.isin(keep) & edges.dst.isin(keep)].reset_index(drop=True))


# --------------------------------------------------------------------------------------
# Coverage and report
# --------------------------------------------------------------------------------------
def distance_to_set(A, targets):
    """Hop distance from every node to the nearest node in `targets` along A (A[i, j]: i -> j)."""
    n = A.shape[0]
    S = sp.vstack([sp.csr_matrix((np.ones(len(targets)), (np.zeros(len(targets)), targets)), shape=(1, n)),
                   A.T.tocsr()])
    S = sp.hstack([sp.csr_matrix((n + 1, 1)), S]).tocsr()
    return csgraph.shortest_path(S, directed=True, indices=0, unweighted=True)[1:] - 1


def transcriptional_route(nodes, edges):
    """Distance of every node to a TF that has a signed regulon edge to a GEM enzyme, along directed edges.
    This is the route an expression-derived metabolic signature can see."""
    enzymes = set(edges.src[edges.relation == "catalysis"])
    grn = edges[edges.relation == "grn"]
    tfs = sorted(set(grn.src[grn.dst.isin(enzymes)]))
    pos = pd.Series(np.arange(len(nodes)), index=nodes.node_id)
    A = adjacency(nodes, edges, ROUTE_RELATIONS)
    return distance_to_set(A, pos[tfs].to_numpy()), tfs


def report(nodes, edges, cfg, drug_targets=None):
    prot = (nodes.node_type == "protein").to_numpy()
    dist, tfs = transcriptional_route(nodes, edges)
    lines = [f"nodes: {nodes.node_type.value_counts().to_dict()}",
             f"edges per relation: {edges.relation.value_counts().to_dict()}",
             f"signed edges per relation: {edges[edges.sign != 0].relation.value_counts().to_dict()}",
             f"TFs with signed regulon edges to GEM enzymes: {len(tfs)}",
             f"proteins within {cfg.max_hops} directed hops of such a TF: {np.mean(dist[prot] <= cfg.max_hops):.3f}"]
    if drug_targets is not None:
        pos = pd.Series(np.arange(len(nodes)), index=nodes.node_id)
        dt = drug_targets[drug_targets.target.isin(pos.index)]
        ok = dist[pos[dt.target].to_numpy()] <= cfg.max_hops
        lines += [f"drug targets in graph: {dt.target.nunique()}/{drug_targets.target.nunique()}",
                  f"drug targets within {cfg.max_hops} hops of the transcriptional route: {ok.mean():.3f}; drugs "
                  f"with >= 1: {dt.assign(ok=ok).groupby('drug').ok.any().sum()}/{drug_targets.drug.nunique()}"]
    present = set(edges.relation)
    lines.append("suggested pipeline weights (expression-derived signal, untuned): "
                 + " ".join(f"--alpha {r}={a}" for r, a in SUGGESTED_ALPHA.items() if r in present))
    return lines


def load_drug_targets(path, hgnc):
    """drug, target[, effect, pertname, ...] (build_drug_targets.py output); targets mapped to HGNC symbols."""
    dt = read_table(path)
    dt = dt.rename(columns={dt.columns[0]: "drug", dt.columns[1]: "target"})
    dt["target"] = dt.target.map(hgnc)
    if "pertname" not in dt:
        dt["pertname"] = dt.drug
    return dt[dt.target != ""].reset_index(drop=True)


def drug_metabolite_matches(drug_targets, nodes):
    """Drugs whose name is a GEM metabolite name: remove those nodes from the drug's own signatures."""
    mets = nodes[nodes.node_type == "metabolite"]
    by_name = mets.groupby(mets.symbol.str.lower()).node_id.apply(list).to_dict()
    drugs = drug_targets[["drug", "pertname"]].drop_duplicates()
    rows = [dict(drug=d, pertname=n, node_id=m) for d, n in zip(drugs.drug, drugs.pertname)
            for m in by_name.get(str(n).lower(), [])]
    return pd.DataFrame(rows, columns=["drug", "pertname", "node_id"])


def check(nodes, edges):
    """Structural invariants of the output."""
    ids = set(nodes.node_id)
    assert nodes.node_id.is_unique, "duplicate node ids"
    assert edges.src.isin(ids).all() and edges.dst.isin(ids).all(), "edge endpoint missing from nodes"
    assert not (edges.src == edges.dst).any(), "self-loops"
    assert not edges.duplicated(KEYS).any(), "duplicate edges"
    assert set(edges.relation) <= set(PROTEIN_RELATIONS) | {"metabolic", "catalysis", "second_messenger"}, \
        "unknown relation"


def write_outputs(out_dir, nodes, edges, mapping, lines, matches=None):
    os.makedirs(out_dir, exist_ok=True)
    nodes.to_csv(os.path.join(out_dir, "nodes.tsv"), sep="\t", index=False)
    cols = ["src", "dst", "relation", "directed", "sign", "n_resources", "n_references", "resources", "origins",
            "weight"]
    edges[cols].to_csv(os.path.join(out_dir, "edges.tsv"), sep="\t", index=False)
    mapping.to_csv(os.path.join(out_dir, "id_mapping.tsv"), sep="\t", index=False)
    if matches is not None:
        matches.to_csv(os.path.join(out_dir, "drug_metabolite_matches.tsv"), sep="\t", index=False)
    with open(os.path.join(out_dir, "report.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


# --------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", required=True, help="metabognn_dev/data")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--drug-targets", help="TSV: drug, target[, effect, pertname] (build_drug_targets.py output)")
    for name, default in DEFAULT_PATHS.items():
        ap.add_argument(f"--{name.replace('_', '-')}", default=default, help=f"relative to --data-dir ({default})")
    for f in fields(Config):
        if f.type is bool:
            ap.add_argument(f"--no-{f.name.replace('_', '-')}", dest=f.name, action="store_false")
        else:
            ap.add_argument(f"--{f.name.replace('_', '-')}", type=f.type, default=f.default)
    args = ap.parse_args()
    paths = {k: os.path.join(args.data_dir, getattr(args, k)) for k in DEFAULT_PATHS}
    cfg = Config(**{f.name: getattr(args, f.name) for f in fields(Config)})

    hgnc = HGNC(paths["hgnc"])
    gem = load_gem(paths["gem"])
    excluded, bookkeeping = gem_exclusions(gem)
    gene_sym = gem_gene_symbols(gem, hgnc)
    nodes_gem, edges_gem, G = build_gem_layer(gem, excluded, bookkeeping, gene_sym)
    print(f"[GEM] {nodes_gem.node_type.value_counts().to_dict()}, {edges_gem.relation.value_counts().to_dict()}")
    node_of = messenger_nodes(paths["gem_metabolites"], set(nodes_gem.node_id))

    raw = pd.concat([edges_gem, gem_messenger_edges(gem, G, node_of, gene_sym), load_sources(paths, cfg, node_of)],
                    ignore_index=True)
    raw, mapping, _ = harmonise(raw, hgnc, set(nodes_gem.node_id))
    print(mapping.to_string(index=False))
    edges, log = curate(merge_edges(raw), cfg, node_of.values())
    print(f"[curate] {log}")
    nodes = assemble_nodes(nodes_gem, edges, hgnc)
    if cfg.prune:
        nodes, edges = keep_metabolic_components(nodes, edges)
    check(nodes, edges)

    drug_targets = matches = None
    if args.drug_targets:
        drug_targets = load_drug_targets(args.drug_targets, hgnc)
        matches = drug_metabolite_matches(drug_targets, nodes)
    lines = report(nodes, edges, cfg, drug_targets)
    write_outputs(args.out_dir, nodes, edges, mapping, lines, matches)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
