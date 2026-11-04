################################################################################
####################----------------------------------------####################
##########----------            MetaboReprogram             ----------##########
# MetPath pathways (Mattei et al. 2023, Metabolites 13:1127; MATLAB code at    #
# github.com/dczielinski/MetPath), step 1: condition-specific production and   #
# consumption ("degradation") pathways of every active metabolite, from one    #
# reference flux state. Independent of the perturbations: run it once per cell #
# line, then score the differential expression with score_metpath.py.          #
#                                                                              #
# Port of calculateFluxState / getActiveNetwork / pruneMatrices /              #
# distanceDirectional / calcPathways / em_decomp:                              #
# 1. Reference flux state: FBA maximising biomass, minimal Euclidean norm      #
#    among optimal solutions (--flux_state fba, as MetPath), or METAFlux on    #
#    the cell line's control expression (--flux_state metaflux), or a file.    #
#    The medium's uptake values are used for fba; METAFlux uses only its ids.  #
# 2. Active network: flux-carrying reactions (biomass removed), carbon-        #
#    containing or listed inorganic metabolites, stoichiometry signed by the   #
#    flux direction.                                                           #
# 3. Currency pairs (ATP/ADP, NAD(P)(H), ...), cofactor pairs (AcCoA/CoA,      #
#    Glu/AKG, ...) and inorganic metabolites are removed from reactions that   #
#    still carry other carbon, giving the traversal matrices of MetPath.       #
# 4. For each metabolite, reactions up to --distance steps upstream            #
#    (production) or downstream (consumption) of it.                           #
# 5. That subnetwork, mass-balanced with demand reactions, is decomposed into  #
#    elementary modes (MILP of em_decomp, Lun 2011); the modes through the     #
#    metabolite's demand are summed with their weights, and reactions below    #
#    --cutoff_fraction of the total are dropped.                               #
#                                                                              #
# Output (--out_dir): metpath_pathways.tsv (metabolite, direction, reaction,   #
# weight, fallback), reference_fluxes.tsv, active_metabolites.tsv.             #
################################################################################

################################################################################
# MODULE LOADING                                                               #
################################################################################

import os
import re
import sys
import json
import time
import argparse
from multiprocessing import Pool

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import milp, LinearConstraint, Bounds
from cobra.io import read_sbml_model
from cobra.util.array import create_stoichiometric_matrix

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "metaflux"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "parse_lincs_data"))
from run_metaflux import FluxProblem, SOLVERS, compute_mras, load_medium, solve_scaled, stoichiometry  # noqa: E402
from gene_id_mapping import symbol_to_ensembl_map, gem_symbol_overrides  # noqa: E402

################################################################################
# METABOLITE CLASSES (BiGG ids, data/curCofPairs.mat and getActiveNetwork.m)   #
################################################################################

CURRENCY_PAIRS = [("atp", "adp"), ("atp", "amp"), ("itp", "idp"), ("gtp", "gdp"), ("utp", "udp"), ("ctp", "cdp"),
                  ("ttp", "tdp"), ("dgtp", "dgdp"), ("dctp", "dcdp"), ("datp", "dadp"), ("dttp", "dtdp"),
                  ("dutp", "dudp"), ("ditp", "didp"), ("fadh2", "fad"), ("nadph", "nadp"), ("nadh", "nad"),
                  ("q10h2", "q10"), ("focytC", "ficytC"), ("q8", "q8h2")]
COFACTOR_PAIRS = [("paps", "pap"), ("thbpt", "thbpt4acam"), ("accoa", "coa"), ("gthrd", "gthox"),
                  ("ascb__L", "dhdascb"), ("10fthf", "thf"), ("mlthf", "dhf"), ("mlthf", "thf"), ("mlthf", "5mthf"),
                  ("mlthf", "methf"), ("5mthf", "thf"), ("amet", "ahcys"), ("ametam", "5mta"), ("gln__L", "glu__L"),
                  ("oaa", "asp__L"), ("glu__L", "akg"), ("akg", "succ"), ("ala__L", "pyr")]
INORGANIC = ["o2", "so2", "so3", "so4", "nh4", "no2", "no3", "fe2", "fe3", "h2o", "co2", "co", "h2o2", "o2s", "h2s",
             "etha", "no", "fe3hox", "3fe4s", "4fe4s", "2fe2s", "etoh", "mobd", "cu", "cu2"]
CARBON = re.compile(r"C(?![a-z])")   # a carbon atom, not Cl / Ca / Co / Cu ... (findregexp 'C\w[^a-z]')
FLUX_TOL_BALANCE = 1e-8


def bigg_to_mets(gem, metabolites_tsv):
    """BiGG id (no compartment) -> {compartment: GEM metabolite id}."""
    tab = pd.read_csv(metabolites_tsv, sep="\t", dtype=str)
    bigg = dict(zip(tab.mets, tab.metBiGGID))
    out = {}
    for m in gem.metabolites:
        b = bigg.get(m.id)
        if isinstance(b, str) and b:
            out.setdefault(b, {})[m.compartment] = m.id
    return out


def met_classes(bigg_map):
    """Per-compartment currency and cofactor pairs, and the inorganic metabolites, as GEM ids."""
    def pairs(lst):
        out = []
        for a, b in lst:
            for comp in set(bigg_map.get(a, {})) & set(bigg_map.get(b, {})):
                out.append((bigg_map[a][comp], bigg_map[b][comp]))
        return out
    inorganic = {mid for b in INORGANIC for mid in bigg_map.get(b, {}).values()}
    return pairs(CURRENCY_PAIRS), pairs(COFACTOR_PAIRS), inorganic


################################################################################
# REFERENCE FLUX STATE                                                         #
################################################################################

def fba_min_norm(gem, medium, solver):
    """calculateFluxState: optimizeCbModel(model, 'max', 1e-6), i.e. maximal biomass, then the solution of
    minimal Euclidean norm among the optimal ones (QP)."""
    with gem:
        if medium is not None:
            gem.medium = {k: v for k, v in medium.items() if k in gem.reactions}
        opt = gem.slim_optimize()
        lb = np.array([r.lower_bound for r in gem.reactions])
        ub = np.array([r.upper_bound for r in gem.reactions])
    print("FBA: maximal biomass %.6g" % opt)
    bio = np.array([r.objective_coefficient != 0 for r in gem.reactions])
    lb[bio] = ub[bio] = opt
    S = stoichiometry(gem)
    x, status = SOLVERS[solver](S, np.zeros(S.shape[1])).solve(lb, ub)
    print("minimal-norm QP: %s, max |S v| %.1e" % (status, np.abs(S @ x).max()))
    return x


def metaflux_state(gem, medium_ids, baseline_f, hgnc_f, solver, mras_floor):
    """METAFlux on one expression profile (e.g. the mean of the cell line's controls, baseline_ctrl.tsv:
    symbol or ENSEMBL id, log2 value)."""
    base = pd.read_csv(baseline_f, sep="\t")
    ids, vals = base.iloc[:, 0].astype(str), base.iloc[:, 1].astype(float)
    if not ids.str.startswith("ENSG").all():
        mapping = symbol_to_ensembl_map(pd.read_csv(hgnc_f, sep="\t", dtype=str), gem_symbol_overrides(gem))
        ids = ids.map(mapping)
    expr = pd.DataFrame({"baseline": vals.to_numpy()}, index=ids.to_numpy())
    expr = expr[expr.index.notna()]
    expr = np.log2(np.exp2(expr).groupby(level=0).sum()).clip(lower=0.0)  # duplicates summed in linear space
    mras = compute_mras(gem, expr, mras_floor)
    problem = FluxProblem(gem, medium_ids)
    x, status = solve_scaled(SOLVERS[solver](problem.S, problem.q), problem, mras.iloc[:, 0].to_numpy())
    print("METAFlux reference state: %s, biomass %.3g" % (status, x[problem.biomass].sum()))
    return x


################################################################################
# ACTIVE NETWORK AND TRAVERSAL MATRICES                                        #
################################################################################

def active_network(gem, S, fluxes, inorganic, tol):
    """getActiveNetwork: flux-carrying reactions without the biomass, carbon or inorganic metabolites, and the
    stoichiometry signed so that every active reaction runs forward."""
    rxn_ok = np.abs(fluxes) > tol
    rxn_ok[[r.objective_coefficient != 0 for r in gem.reactions]] = False
    rxns = np.where(rxn_ok)[0]
    formula = [m.formula or "" for m in gem.metabolites]
    used = np.array([bool(CARBON.search(f)) or m.id in inorganic for f, m in zip(formula, gem.metabolites)])
    in_active = np.asarray(abs(S[:, rxns]).sum(1)).ravel() > 0
    mets = np.where(used & in_active)[0]
    S_dir = (S[mets][:, rxns] @ sp.diags(np.sign(fluxes[rxns]))).tocsr()
    return rxns, mets, S_dir


def prune_matrices(gem, S, rxns, mets, S_dir, currency, cofactor, inorganic):
    """pruneMatrices: traversal matrices without inorganic metabolites (noIno), also without currency pairs
    (noCurNoIno), and also without cofactor pairs (noCurNoCofNoIno), removed only from reactions that keep
    some other carbon."""
    met_ids = [m.id for m in gem.metabolites]
    row = {met_ids[m]: i for i, m in enumerate(mets)}
    carbon = np.array([bool(CARBON.search(m.formula or "")) for m in gem.metabolites])
    drop = {k: [] for k in ("noIno", "noCurNoIno", "noCurNoCofNoIno")}
    Sc = S.tocsc()
    for j, r in enumerate(rxns):
        present = set(Sc.indices[Sc.indptr[r]:Sc.indptr[r + 1]])
        names = {met_ids[i] for i in present}
        n_carb = int(carbon[list(present)].sum())
        cur = [m for a, b in currency if a in names and b in names for m in (a, b)]
        cof = [m for a, b in cofactor if a in names and b in names for m in (a, b)]
        ino = [m for m in names if m in inorganic]
        if cur and n_carb - len(cur) > 0:
            for m in cur:
                if m in row:
                    drop["noCurNoIno"].append((row[m], j))
                    drop["noCurNoCofNoIno"].append((row[m], j))
        for m in ino:
            if m in row:
                for k in drop:
                    drop[k].append((row[m], j))
        if cof and n_carb - (len(cur) + len(cof)) > 0:
            for m in cof:
                if m in row:
                    drop["noCurNoCofNoIno"].append((row[m], j))
    out = {"full": traversal(S_dir)}
    for k, pairs in drop.items():
        M = S_dir.tolil(copy=True)
        for i, j in pairs:
            M[i, j] = 0
        M = M.tocsr()
        M.eliminate_zeros()
        out[k] = traversal(M)
    return out


def traversal(S_dir):
    """A traversal matrix with its positive (produced) and negative (consumed) parts in both sparse layouts."""
    pos, neg = S_dir.multiply(S_dir > 0), S_dir.multiply(S_dir < 0)
    return dict(S=S_dir.tocsr(), pos=pos.tocsc(), posr=pos.tocsr(), neg=neg.tocsr(), negc=neg.tocsc())


def distance_directional(T, i, depth, forward):
    """distanceDirectional: steps from metabolite i to every reaction, following production backwards
    (forward = False) or consumption forwards (forward = True); -1 beyond depth. T: traversal(S_dir)."""
    S_dir, pos, neg, posr, negc = T["S"], T["pos"], T["neg"], T["posr"], T["negc"]
    row = S_dir.getrow(i).toarray().ravel()
    cons, prod = np.where(row < 0)[0], np.where(row > 0)[0]
    dist = -np.ones(S_dir.shape[1], dtype=int)
    if forward:
        allowed = np.ones(S_dir.shape[1], bool)
        allowed[prod] = False           # sDirDeg: reactions producing the metabolite are cut
        cur = cons
        dist[cons] = 0
        for d in range(1, depth + 1):
            m = np.unique(pos[:, cur].nonzero()[0])                      # produced at this level
            r = np.unique(neg[m].nonzero()[1]) if len(m) else np.array([], int)   # consumers of those
            cur = r[allowed[r]]
            dist[cur[dist[cur] == -1]] = d
    else:
        allowed = np.ones(S_dir.shape[1], bool)
        allowed[cons] = False           # sDirProd: reactions consuming the metabolite are cut
        cur = prod
        dist[prod] = 0
        for d in range(1, depth + 1):
            m = np.unique(negc[:, cur].nonzero()[0])                     # consumed at this level
            r = np.unique(posr[m].nonzero()[1]) if len(m) else np.array([], int)  # producers of those
            cur = r[allowed[r]]
            dist[cur[dist[cur] == -1]] = d
    return dist


################################################################################
# ELEMENTARY MODE DECOMPOSITION (em_decomp.m, Desmond Lun 2011)                #
################################################################################

def em_decomp(v, S, big_m=100.0, eps=1e-6, time_limit=200.0, max_modes=None):
    """Decompose the steady-state flux v (S v = 0) into elementary modes: repeatedly take the reaction with the
    largest remaining flux, find the mode through it with the fewest reactions that conforms to v (MILP), and
    subtract it with the largest weight that keeps v conformal. Returns (P n x k, w k)."""
    v = v.astype(float).copy()
    n, n_met = len(v), S.shape[0]
    P, W = [], []
    max_modes = max_modes or 5 * n
    c = np.concatenate([np.zeros(n), np.ones(n)])
    integrality = np.concatenate([np.zeros(n), np.ones(n)])
    bounds = Bounds(np.concatenate([np.full(n, -np.inf), np.zeros(n)]), np.concatenate([np.full(n, np.inf), np.ones(n)]))
    I = sp.identity(n, format="csr")
    while np.abs(v).max() > eps and len(P) < max_modes:
        k = int(np.argmax(np.abs(v)))
        zero, nz = np.where(v == 0)[0], np.where(v != 0)[0]
        eq = sp.vstack([sp.hstack([S, sp.csr_matrix((n_met, n))]),
                        sp.hstack([I[zero], sp.csr_matrix((len(zero), n))]),
                        sp.hstack([I[[k]], sp.csr_matrix((1, n))])], format="csr")
        beq = np.concatenate([np.zeros(n_met + len(zero)), [np.sign(v[k])]])
        ineq = sp.vstack([sp.hstack([sp.diags(-np.sign(v[nz])) @ I[nz], sp.csr_matrix((len(nz), n))]),
                          sp.hstack([I, -big_m * I]), sp.hstack([-I, -big_m * I])], format="csr")
        res = milp(c, integrality=integrality, bounds=bounds,
                   constraints=[LinearConstraint(eq, beq, beq), LinearConstraint(ineq, -np.inf, 0)],
                   options={"time_limit": time_limit, "disp": False})
        if res.x is None:
            p = v.copy()                     # failed: the remainder becomes the last mode
        else:
            p = res.x[:n].copy()
            p[res.x[n:] < 0.5] = 0.0
        p /= np.sqrt((p ** 2).sum())
        sel = p != 0
        w = np.min(v[sel] / p[sel])
        if not np.isfinite(w) or w <= 0:
            break
        P.append(p)
        W.append(w)
        v -= w * p
        v[np.abs(v) <= eps * 1e-3] = 0.0
    return (np.array(P).T if P else np.zeros((n, 0))), np.array(W)


################################################################################
# PATHWAYS                                                                     #
################################################################################

_G = {}


def _init(state):
    _G.update(state)


def met_pathway(args):
    """calcPathways for one metabolite and direction: subnetwork within the distance, balanced with demand
    reactions, decomposed into elementary modes; the modes through the metabolite's demand, summed."""
    i, forward = args
    g = _G
    met_id = g["met_ids"][g["mets"][i]]
    S_trav = g["S_trav"][g["met_class"][i]]
    dist = distance_directional(S_trav, i, g["distance"], forward)
    sub = np.where(dist > -1)[0]
    if len(sub) == 0:
        return met_id, forward, [], [], False
    rxn_cols = g["rxns"][sub]
    v = g["fluxes"][rxn_cols]
    S_sub = g["S"][:, rxn_cols]
    bal = S_sub @ v
    bal[np.abs(bal) < FLUX_TOL_BALANCE] = 0.0
    unb = np.where(bal != 0)[0]
    # demand reaction (+1 on the metabolite) carrying -imbalance, as addReactionNoDup('DM_' met)
    S_full = sp.hstack([S_sub, sp.csr_matrix((np.ones(len(unb)), (unb, np.arange(len(unb)))),
                                             shape=(S_sub.shape[0], len(unb)))], format="csr")
    v_full = np.concatenate([v, -bal[unb]])
    rows = np.where(np.asarray(abs(S_full).sum(1)).ravel() > 0)[0]
    target = np.where(unb == g["mets"][i])[0]
    try:
        P, w = em_decomp(v_full, S_full[rows], time_limit=g["time_limit"])
        if len(target) == 0:
            return met_id, forward, [], [], False
        dm = len(sub) + target[0]
        through = np.where(P[dm] != 0)[0]
        if len(through) == 0:
            return met_id, forward, [], [], False
        summed = (P[:, through] * w[through]).sum(1)[:len(sub)]
        total = np.abs((P[:, through] * w[through]).sum(1)).sum()
        summed[np.abs(summed / total) < g["cutoff_fraction"]] = 0.0
        keep = summed != 0
        return met_id, forward, list(rxn_cols[keep]), list(summed[keep]), False
    except Exception as exc:  # MetPath's catch: every reaction of the subnetwork with weight 1
        print("  %s %s: decomposition failed (%s), all %d reactions weighted 1" % (
            met_id, "consumption" if forward else "production", exc, len(sub)), flush=True)
        return met_id, forward, list(rxn_cols), [1.0] * len(sub), True


################################################################################
# MAIN                                                                         #
################################################################################

def main():
    parser = argparse.ArgumentParser(description="MetPath production / consumption pathways from a reference flux state.")
    parser.add_argument("--gem", type=str, required=True, help="GEM file (SBML).")
    parser.add_argument("--metabolites_tsv", type=str, default=None,
                        help="Human-GEM metabolites.tsv (BiGG ids for the currency / cofactor / inorganic classes); "
                             "default: next to --gem.")
    parser.add_argument("--medium", type=str, required=True, help="Medium file (JSON, exchange reaction id -> uptake).")
    parser.add_argument("--flux_state", choices=["fba", "metaflux", "file"], default="fba",
                        help="fba: max biomass, minimal norm (MetPath); metaflux: METAFlux on --baseline; "
                             "file: --fluxes.")
    parser.add_argument("--baseline", type=str, default=None,
                        help="For --flux_state metaflux: expression profile (TSV: symbol or ENSEMBL id, log2 value), "
                             "e.g. prepare_inputs.py's baseline_ctrl.tsv.")
    parser.add_argument("--hgnc", type=str, default=None, help="HGNC table, to map symbols in --baseline.")
    parser.add_argument("--fluxes", type=str, default=None, help="For --flux_state file: TSV reaction id, flux.")
    parser.add_argument("--solver", choices=["clarabel", "osqp", "gurobi"], default="clarabel",
                        help="QP solver for the reference flux state (see run_metaflux.py).")
    parser.add_argument("--mras_floor", type=float, default=0.0,
                        help="For --flux_state metaflux: see run_metaflux.py --mras_floor.")
    parser.add_argument("--tol_flux", type=float, default=1e-6, help="Smallest flux considered non-zero (tolFlux).")
    parser.add_argument("--distance", type=int, default=2,
                        help="Pathway length in reactions beyond the first (cutoffDistance; MetPath suggests 2-3).")
    parser.add_argument("--cutoff_fraction", type=float, default=0.05,
                        help="Reactions contributing less than this fraction of the pathway are dropped.")
    parser.add_argument("--time_limit", type=float, default=200.0, help="Seconds per elementary-mode MILP.")
    parser.add_argument("--max_mets", type=int, default=0, help="Only the first N active metabolites (tests).")
    parser.add_argument("--n_workers", type=int, default=1, help="Number of parallel workers.")
    parser.add_argument("--out_dir", type=str, default=os.getcwd(), help="Directory where the result will be stored.")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    print("Loading GEM: %s" % args.gem)
    gem = read_sbml_model(args.gem)
    S = create_stoichiometric_matrix(gem, array_type="dok").tocsr()
    rxn_ids = [r.id for r in gem.reactions]
    with open(args.medium) as fh:
        medium = json.load(fh)

    # 1. reference flux state
    if args.flux_state == "fba":
        fluxes = fba_min_norm(gem, medium, args.solver)
    elif args.flux_state == "metaflux":
        fluxes = metaflux_state(gem, load_medium(args.medium, gem), args.baseline, args.hgnc, args.solver,
                                args.mras_floor)
    else:
        f = pd.read_csv(args.fluxes, sep="\t", index_col=0).iloc[:, 0]
        fluxes = f.reindex(rxn_ids).fillna(0.0).to_numpy()
    fluxes = np.where(np.abs(fluxes) < args.tol_flux, 0.0, fluxes)
    pd.DataFrame({"reaction": rxn_ids, "flux": fluxes}).to_csv(
        os.path.join(args.out_dir, "reference_fluxes.tsv"), sep="\t", index=False, float_format="%.9g")

    # 2-3. active network, metabolite classes, traversal matrices
    mets_tsv = args.metabolites_tsv or os.path.join(os.path.dirname(args.gem), "metabolites.tsv")
    currency, cofactor, inorganic = met_classes(bigg_to_mets(gem, mets_tsv))
    rxns, mets, S_dir = active_network(gem, S, fluxes, inorganic, args.tol_flux)
    print("active network: %d reactions, %d metabolites" % (len(rxns), len(mets)))
    S_trav = prune_matrices(gem, S, rxns, mets, S_dir, currency, cofactor, inorganic)
    met_ids = [m.id for m in gem.metabolites]
    cur_set = {m for p in currency for m in p}
    cof_set = {m for p in cofactor for m in p}
    met_class = ["noIno" if met_ids[m] in cur_set else "noCurNoIno" if met_ids[m] in cof_set
                 else "full" if met_ids[m] in inorganic else "noCurNoCofNoIno" for m in mets]
    pd.DataFrame({"metabolite": [met_ids[m] for m in mets], "name": [gem.metabolites[m].name for m in mets],
                  "class": met_class}).to_csv(os.path.join(args.out_dir, "active_metabolites.tsv"), sep="\t",
                                              index=False)

    # 4-5. pathways
    n = len(mets) if not args.max_mets else min(args.max_mets, len(mets))
    state = dict(met_ids=met_ids, mets=mets, rxns=rxns, S=S, S_trav=S_trav, met_class=met_class, fluxes=fluxes,
                 distance=args.distance, cutoff_fraction=args.cutoff_fraction, time_limit=args.time_limit)
    tasks = [(i, fwd) for i in range(n) for fwd in (False, True)]
    start = time.time()
    if args.n_workers > 1:
        with Pool(args.n_workers, initializer=_init, initargs=(state,)) as pool:
            results = []
            for k, res in enumerate(pool.imap_unordered(met_pathway, tasks, chunksize=4)):
                results.append(res)
                if (k + 1) % 200 == 0:
                    print("  %d / %d pathways (%.0f s)" % (k + 1, len(tasks), time.time() - start), flush=True)
    else:
        _init(state)
        results = []
        for k, t in enumerate(tasks):
            results.append(met_pathway(t))
            if (k + 1) % 200 == 0:
                print("  %d / %d pathways (%.0f s)" % (k + 1, len(tasks), time.time() - start), flush=True)
    rows = [(m, "consumption" if fwd else "production", rxn_ids[r], w, fb)
            for m, fwd, rs, ws, fb in results for r, w in zip(rs, ws)]
    out = pd.DataFrame(rows, columns=["metabolite", "direction", "reaction", "weight", "fallback"])
    out.to_csv(os.path.join(args.out_dir, "metpath_pathways.tsv"), sep="\t", index=False, float_format="%.6g")
    n_paths = out.groupby(["metabolite", "direction"]).ngroups
    print("%d pathways for %d metabolites in %.0f s (%d decompositions failed); median %d reactions per pathway" % (
        n_paths, out.metabolite.nunique(), time.time() - start, sum(r[4] for r in results),
        out.groupby(["metabolite", "direction"]).size().median() if len(out) else 0))


if __name__ == "__main__":
    main()
