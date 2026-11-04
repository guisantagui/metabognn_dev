################################################################################
####################----------------------------------------####################
##########----------            MetaboReprogram             ----------##########
# Run METAFlux (Huang et al. 2023, Nat Commun 14:4883; github.com/KChen-lab/  #
# METAFlux) on a genes x samples expression matrix, sample by sample:          #
#                                                                              #
# 1. MRAS (metabolic reaction activity scores): each gene is divided by the    #
#    number of GEM reactions it takes part in, the GPR is evaluated with       #
#    AND = min and OR = sum (missing genes ignored), and every sample is       #
#    divided by its maximum over reactions. Reactions with a GPR but no        #
#    measured gene get 0, reactions without a GPR get 1.                       #
# 2. Quadratic program per sample:                                             #
#        min 1/2 ||v||^2 - 10000 v_biomass                                     #
#        s.t. S v = 0,  -MRAS <= v <= MRAS (reversible), 0 <= v <= MRAS        #
#    with exchange reactions closed for uptake except the medium ones (lower   #
#    bound -1).                                                                #
#                                                                              #
# Same as calculate_reaction_score() / compute_flux() of the R package, with   #
# these differences:                                                           #
# - GPRs of any depth are evaluated recursively (the package handles up to two #
#   levels of nesting); the model is the Human-GEM given with --gem.           #
# - The QP is solved after rescaling bounds and objective by the sample's      #
#   median positive upper bound s (an exactly equivalent problem: v = s u      #
#   gives min 1/2 ||u||^2 - (10000 / s) u_biomass with bounds / s), so solver  #
#   tolerances are relative to the fluxes, which are ~1e-4 with Human-GEM.     #
# - --solver clarabel (default): interior-point QP solver, accurate to ~1e-9   #
#   of the flux scale. osqp: the package's solver and settings (approximate,   #
#   and slow to converge on this LP-like problem). gurobi: exact.              #
# - --mras_floor (default 0, as the package): reactions with a GPR whose genes #
#   are unmeasured or at 0 get this quantile of the sample's positive MRAS     #
#   instead of 0. In LINCS A549 some biomass-essential reactions read 0 in     #
#   every sample, so without a floor there is no biomass and every flux is 0;  #
#   below a quantile of ~0.25 the floor itself sets the biomass flux.          #
#                                                                              #
# Input: a batch written by prepare_input_matrices.py (--input_type            #
# expression): GEM genes (ENSEMBL) x samples, linear scale. The output is      #
# named like the MetChange batches, so get_delta_mat.py and limma_DE.R can be  #
# used on it unchanged.                                                        #
################################################################################

################################################################################
# MODULE LOADING                                                               #
################################################################################

import os
import ast
import json
import time
import argparse
from collections import Counter
from multiprocessing import Pool

import numpy as np
import pandas as pd
import scipy.sparse as sp
from cobra.io import read_sbml_model
from cobra.util.array import create_stoichiometric_matrix

################################################################################
# FUNCTIONS                                                                    #
################################################################################

BIOMASS_WEIGHT = 10000.0   # q[biomass] = -10000 in compute_flux()
EXCHANGE_SUBSYSTEM = "Exchange/demand reactions"


# MRAS
################################################################################

def compile_gpr(node):
    """GPR AST -> function of {gene: array over samples}. Missing genes are NaN and ignored, as with
    na.rm = TRUE in the R package; a rule with no measured gene gives NaN."""
    if isinstance(node, ast.Name):
        return lambda expr, g=node.id: expr[g]
    if isinstance(node, ast.BoolOp):
        children = [compile_gpr(child) for child in node.values]
        if isinstance(node.op, ast.And):
            def fn(expr, children=children):
                x = np.vstack([f(expr) for f in children])
                out = np.min(np.where(np.isnan(x), np.inf, x), 0)
                return np.where(np.isinf(out), np.nan, out)
        else:
            def fn(expr, children=children):
                x = np.vstack([f(expr) for f in children])
                return np.where(np.isnan(x).all(0), np.nan, np.nansum(x, 0))
        return fn
    raise ValueError("Unknown GPR node: %s" % node)


def gene_reaction_counts(gem):
    """Number of GEM reactions each gene takes part in (the package's gene_num)."""
    return Counter(g.id for r in gem.reactions for g in r.genes)


def compute_mras(gem, expr_df, floor_quantile=0.0):
    """expr_df: genes x samples (non-negative). Returns reactions x samples (MRAS)."""
    n_samples = expr_df.shape[1]
    counts = gene_reaction_counts(gem)
    values = expr_df.to_numpy(dtype=np.float64)
    missing = np.full(n_samples, np.nan)
    gene_arrays = {g: values[i] / counts[g] for i, g in enumerate(expr_df.index) if counts.get(g)}

    class Lookup(dict):
        def __missing__(self, key):
            return missing

    lookup = Lookup(gene_arrays)
    mras = np.full((len(gem.reactions), n_samples), np.nan)
    has_gpr = np.zeros(len(gem.reactions), bool)
    for j, rxn in enumerate(gem.reactions):
        if rxn.genes:
            has_gpr[j] = True
            mras[j] = compile_gpr(rxn.gpr.body)(lookup)
    # stdize(): divide every sample by its maximum over reactions; then NaN -> 0
    col_max = np.nanmax(mras[has_gpr], axis=0)
    mras = mras / np.where(col_max > 0, col_max, 1.0)
    mras[np.isnan(mras)] = 0.0
    blocked = np.array([r.lower_bound == 0 and r.upper_bound == 0 for r in gem.reactions])
    if floor_quantile > 0:
        for c in range(n_samples):
            col = mras[has_gpr, c]
            floor = np.quantile(col[col > 0], floor_quantile)
            mras[has_gpr & ~blocked & (mras[:, c] < floor), c] = floor
    mras[~has_gpr] = 1.0
    mras[blocked] = 0.0
    return pd.DataFrame(mras, index=[r.id for r in gem.reactions], columns=expr_df.columns)


# Quadratic program
################################################################################

def load_medium(medium_f, gem):
    """Medium JSON (cobra medium: exchange reaction id -> uptake). Only the ids are used: METAFlux opens every
    medium exchange with a lower bound of -1, whatever the uptake value."""
    with open(medium_f) as fh:
        medium = json.load(fh)
    rxn_ids = {r.id for r in gem.reactions}
    unknown = sorted(set(medium) - rxn_ids)
    if unknown:
        print("%d medium reactions not in the GEM (ignored), e.g. %s" % (len(unknown), unknown[:5]))
    return {k for k in medium if k in rxn_ids}


def stoichiometry(gem, ignore_compartments=()):
    """Mass-balance rows of S. Human-GEM has no boundary compartment ('x' is the peroxisome), so by default every
    metabolite is balanced."""
    S = create_stoichiometric_matrix(gem, array_type="dok").tocsr()
    keep = [i for i, m in enumerate(gem.metabolites) if m.compartment not in ignore_compartments]
    return S[keep].tocsc()


class FluxProblem:
    """Static part of the QP (S, objective, which bounds come from MRAS); bounds are set per sample."""

    def __init__(self, gem, medium):
        self.S = stoichiometry(gem)
        self.rxn_ids = [r.id for r in gem.reactions]
        lb = np.array([r.lower_bound for r in gem.reactions])
        ub = np.array([r.upper_bound for r in gem.reactions])
        self.reversible = (lb < 0) & (ub > 0)
        self.reverse_only = (ub <= 0) & (lb < 0)
        self.exchange = np.array([r.subsystem == EXCHANGE_SUBSYSTEM for r in gem.reactions])
        self.in_medium = np.array([r.id in medium for r in gem.reactions])
        self.biomass = np.array([r.objective_coefficient != 0 for r in gem.reactions])
        self.q = np.zeros(len(self.rxn_ids))
        self.q[self.biomass] = -BIOMASS_WEIGHT

    def bounds(self, mras):
        lb = np.where(self.reversible | self.reverse_only, -mras, 0.0)
        ub = np.where(self.reverse_only, 0.0, mras)
        lb[self.exchange] = 0.0
        lb[self.exchange & self.in_medium] = -1.0
        return lb, ub


class OSQPSolver:
    """min 1/2 ||v||^2 + q v s.t. S v = 0, lb <= v <= ub. Default settings are those of compute_flux():
    max_iter 1e6, eps_abs = eps_rel = 1e-4, adaptive_rho_interval 50."""

    def __init__(self, S, q, eps=1e-4):
        import osqp
        n, m = S.shape[1], S.shape[0]
        self.m = m
        A = sp.vstack([S, sp.identity(n)], format="csc")
        self.model = osqp.OSQP()
        self.model.setup(sp.identity(n, format="csc"), q, A, np.zeros(m + n), np.zeros(m + n), max_iter=1000000,
                         eps_abs=eps, eps_rel=eps, adaptive_rho_interval=50, polish=False, verbose=False)

    def solve(self, lb, ub, q=None):
        if q is not None:
            self.model.update(q=q)
        self.model.update(l=np.concatenate([np.zeros(self.m), lb]), u=np.concatenate([np.zeros(self.m), ub]))
        res = self.model.solve()
        ok = res.info.status in ("solved", "solved inaccurate")
        return (res.x if ok else np.full(len(lb), np.nan)), res.info.status


class ClarabelSolver:
    """Interior point (Clarabel): accurate (|S v| ~1e-9 of the flux scale), no licence needed."""

    def __init__(self, S, q, tol=1e-8):
        self.S, self.q, self.tol = S.tocsc(), q, tol

    def solve(self, lb, ub, q=None):
        import clarabel
        q = self.q if q is None else q
        S, n, m = self.S, self.S.shape[1], self.S.shape[0]
        fu, fl = np.isfinite(ub), np.isfinite(lb)
        I = sp.identity(n, format="csc")
        A = sp.vstack([S, I[fu], -I[fl]], format="csc")
        b = np.concatenate([np.zeros(m), ub[fu], -lb[fl]])
        st = clarabel.DefaultSettings()
        st.verbose = False
        st.max_iter = 1000
        st.tol_gap_abs = st.tol_gap_rel = st.tol_feas = self.tol
        sol = clarabel.DefaultSolver(sp.identity(n, format="csc"), q, A, b,
                                     [clarabel.ZeroConeT(m), clarabel.NonnegativeConeT(int(fu.sum() + fl.sum()))],
                                     st).solve()
        status = str(sol.status).split(".")[-1]
        if status not in ("Solved", "AlmostSolved"):
            return np.full(n, np.nan), status
        return np.clip(np.array(sol.x), lb, ub), status


class GurobiSolver:
    """Exact QP with Gurobi."""

    def __init__(self, S, q, threads=1):
        import gurobipy as gp
        from gurobipy import GRB
        self.GRB = GRB
        n = S.shape[1]
        self.model = gp.Model()
        self.model.Params.OutputFlag = 0
        self.model.Params.Threads = threads
        self.v = self.model.addMVar(n, lb=-1.0, ub=1.0)
        self.model.addMConstr(S, self.v, "=", np.zeros(S.shape[0]))
        self.linear = q
        self.model.setObjective(0.5 * (self.v @ self.v) + q @ self.v, GRB.MINIMIZE)

    def solve(self, lb, ub, q=None):
        if q is not None and not np.array_equal(q, self.linear):
            self.linear = q
            self.model.setObjective(0.5 * (self.v @ self.v) + q @ self.v, self.GRB.MINIMIZE)
        self.v.lb, self.v.ub = lb, ub
        self.model.optimize()
        if self.model.Status == self.GRB.OPTIMAL:
            return self.v.X.copy(), "Solved"
        return np.full(len(lb), np.nan), "gurobi status %d" % self.model.Status


SOLVERS = {"clarabel": ClarabelSolver, "osqp": OSQPSolver, "gurobi": GurobiSolver}
OK_STATUS = {"solved", "solved inaccurate", "Solved", "AlmostSolved"}


def solve_scaled(solver, problem, mras):
    """One sample's QP, solved as the exactly equivalent problem in u = v / s (s: median positive upper bound)."""
    lb, ub = problem.bounds(mras)
    s = float(np.median(ub[ub > 0]))
    x, status = solver.solve(lb / s, ub / s, problem.q / s)
    return x * s, status


def _worker(args):
    problem, mras_chunk, solver_name = args
    solver = SOLVERS[solver_name](problem.S, problem.q)
    fluxes = np.full(mras_chunk.shape, np.nan)
    status = []
    for c in range(mras_chunk.shape[1]):
        x, st = solve_scaled(solver, problem, mras_chunk[:, c])
        if st in OK_STATUS:
            fluxes[:, c] = x
        status.append(st)
    return fluxes, status


def compute_flux(problem, mras, solver="clarabel", n_workers=1):
    """mras: reactions x samples (array). Returns (fluxes reactions x samples, solver status per sample)."""
    chunks = [c for c in np.array_split(np.arange(mras.shape[1]), max(n_workers, 1)) if len(c)]
    args = [(problem, mras[:, c], solver) for c in chunks]
    if len(args) == 1:
        results = [_worker(args[0])]
    else:
        with Pool(len(args)) as pool:
            results = pool.map(_worker, args)
    return np.hstack([r[0] for r in results]), [s for r in results for s in r[1]]


def read_expression(expr_f, scale):
    """Genes x samples. METAFlux expects normalised, non-negative expression (log scale in the paper);
    prepare_input_matrices.py writes linear values (2^log2), so they are taken back to log2 here."""
    expr = pd.read_csv(expr_f, index_col=0)
    if scale == "linear":
        expr = np.log2(expr.clip(lower=1.0))
    n_neg = int((expr < 0).to_numpy().sum())
    if n_neg:
        print("%d negative expression values set to 0 (METAFlux requires non-negative input)" % n_neg)
        expr = expr.clip(lower=0.0)
    return expr


################################################################################
# MAIN                                                                         #
################################################################################

def main():
    parser = argparse.ArgumentParser(description="Run METAFlux (MRAS + QP flux estimation) on an expression matrix.")
    parser.add_argument("--expr", type=str, required=True,
                        help="Expression file: GEM genes (ENSEMBL) x samples, e.g. a prepare_input_matrices.py batch.")
    parser.add_argument("--expr_scale", choices=["linear", "log2"], default="linear",
                        help="linear: values written by prepare_input_matrices.py --input_type expression (converted "
                             "back to log2); log2: used as they are.")
    parser.add_argument("--gem", type=str, required=True, help="GEM file.")
    parser.add_argument("--medium", type=str, required=True, help="Medium file (JSON, exchange reaction ids).")
    parser.add_argument("--solver", choices=["clarabel", "osqp", "gurobi"], default="clarabel",
                        help="clarabel: interior point (accurate, no licence); osqp: the R package's solver and "
                             "settings (approximate, slow here); gurobi: exact.")
    parser.add_argument("--mras_floor", type=float, default=0.0,
                        help="Quantile of the sample's positive MRAS given to reactions whose genes are unmeasured "
                             "or at 0 (0: they get 0, as in the package). See the header.")
    parser.add_argument("--n_workers", type=int, default=1, help="Number of parallel workers.")
    parser.add_argument("--save_mras", action="store_true", help="Also write the MRAS matrix.")
    parser.add_argument("--status_dir", type=str, default=None,
                        help="Where to write the per-sample solver status and biomass flux (keep it outside "
                             "--out_dir: get_delta_mat.py reads every file in that directory).")
    parser.add_argument("--out_dir", type=str, default=os.getcwd(), help="Directory where the result will be stored.")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    name = "%s_metaflux_%s" % (os.path.basename(args.expr).split(".")[0],
                               os.path.basename(args.medium).split(".")[0])

    print("Loading GEM: %s" % args.gem)
    gem = read_sbml_model(args.gem)
    print("Loading expression data from %s" % args.expr)
    expr = read_expression(args.expr, args.expr_scale)
    gem_genes = {g.id for g in gem.genes}
    print("%d of %d GEM genes measured, %d samples" % (len(gem_genes & set(expr.index)), len(gem_genes),
                                                        expr.shape[1]))

    start = time.time()
    mras = compute_mras(gem, expr, args.mras_floor)
    print("MRAS computed: %.0f seconds" % (time.time() - start))
    if args.save_mras:
        mras.to_csv(os.path.join(args.out_dir, "%s_mras.csv" % name), float_format="%.6g")

    problem = FluxProblem(gem, load_medium(args.medium, gem))
    start = time.time()
    fluxes, status = compute_flux(problem, mras.to_numpy(), args.solver, args.n_workers)
    print("Fluxes computed: %.0f seconds (%.2f s per sample)" % (time.time() - start,
                                                                   (time.time() - start) / expr.shape[1]))
    print("Solver status: %s" % dict(Counter(status)))
    bio = fluxes[problem.biomass].sum(0)
    n_zero = int(np.sum(np.nan_to_num(bio) <= 1e-12))
    if n_zero:
        print("WARNING: %d of %d samples have no biomass flux, so all their fluxes are 0 (see --mras_floor)"
              % (n_zero, len(bio)))
    pd.DataFrame(fluxes, index=problem.rxn_ids, columns=expr.columns).to_csv(
        os.path.join(args.out_dir, "%s.csv" % name), float_format="%.6g")
    if args.status_dir:
        os.makedirs(args.status_dir, exist_ok=True)
        pd.DataFrame({"sample": expr.columns, "status": status, "biomass_flux": bio}).to_csv(
            os.path.join(args.status_dir, "%s_status.csv" % name), index=False)


if __name__ == "__main__":
    main()
