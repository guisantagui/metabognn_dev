################################################################################
####################----------------------------------------####################
##########----------            MetaboReprogram             ----------##########
# Run metChange algorithm:                                                     #
#                                                                              #
################################################################################

################################################################################
# MODULE LOADING                                                               #
################################################################################

import os
import argparse
import numpy as np
import pandas as pd
import time
import copy as cp
import casadi as ca
import scipy.sparse as sp
from scipy.sparse import coo_matrix
import ast
import json

import sys

sys.path.insert(
    0,
    "/opt/apps/easybuild/systems/iris/rhel810-20260107/2024a/broadwell/software/Gurobi/12.0.1-GCCcore-13.3.0/lib/python3.12/site-packages"
)

import gurobipy as gp
from gurobipy import GRB
print(gp.gurobi.version())
print(gp.__file__)

import optlang.gurobi_interface
import optlang
print(optlang.__version__)
import cobra
from cobra.io import read_sbml_model, write_sbml_model
from cobra import Reaction
from cobra.util.solver import solvers
from cobra.util.array import create_stoichiometric_matrix
print(solvers.keys())

################################################################################
# FUNCTIONS                                                                    #
################################################################################

# Functions for obtaining reaction-level p-values (penalties)
################################################################################

# Functions for mapping gene expression to reactions via the GEM's GPR

"""
Vectorized GPR-based reaction expression scoring.

Compiles GPRs from the GEM, and use them to go from gene-level penalties to
reaction-level penalties.
"""

# ---------------------------------------------------------------------
# 1. Compile each reaction's GPR AST into a function operating on arrays
#    (shape: n_samples per gene) instead of scalars. Done once per reaction,
#    reused for every sample.
# ---------------------------------------------------------------------
def compile_gpr(node):
    if isinstance(node, ast.Name):
        gene_id = node.id
        return lambda expr: expr[gene_id]

    elif isinstance(node, ast.BoolOp):
        child_fns = [compile_gpr(child) for child in node.values]
        if isinstance(node.op, ast.And):
            def fn(expr, child_fns=child_fns):
                return np.maximum.reduce([f(expr) for f in child_fns])
            return fn
        elif isinstance(node.op, ast.Or):
            def fn(expr, child_fns=child_fns):
                return np.minimum.reduce([f(expr) for f in child_fns])
            return fn

    raise ValueError(f"Unknown node: {node}")


# ---------------------------------------------------------------------
# 2. Precompute reaction-level, profile-invariant info once for the model.
# ---------------------------------------------------------------------
def prepare_reactions(gem, expr_index):
    """
    expr_index: set of gene IDs present in the expression matrix (assumed
    the same across all samples/profiles).
    """
    rxn_info = {}
    for rxn in gem.reactions:
        rxn_genes = [g.id for g in rxn.genes]

        if not rxn_genes:
            rxn_info[rxn.id] = {"const": True}
            continue

        missing = [g for g in rxn_genes if g not in expr_index]
        functional = rxn.gpr.eval(missing)   # expensive call -> now done once
        compiled_fn = compile_gpr(rxn.gpr.body)

        rxn_info[rxn.id] = {
            "const": False,
            "genes": rxn_genes,
            "missing": set(missing),
            "functional": functional,
            "fn": compiled_fn,
        }
    return rxn_info


# ---------------------------------------------------------------------
# 3. Evaluate all reactions x all samples at once.
# ---------------------------------------------------------------------
def get_g_matrix(gem, expr_df):
    """
    expr_df: DataFrame, genes as index, samples/profiles as columns.
    Returns: DataFrame, reactions as index, samples/profiles as columns.
    """
    medians = expr_df.median(axis=0)          # one value per sample
    median_arr = medians.to_numpy()
    n_samples = expr_df.shape[1]

    expr_index = set(expr_df.index)
    rxn_info = prepare_reactions(gem, expr_index)

    gene_to_row = {g: i for i, g in enumerate(expr_df.index)}
    expr_values = expr_df.to_numpy()          # (n_genes, n_samples)

    results = {}
    for rxn_id, info in rxn_info.items():
        if info["const"]:
            results[rxn_id] = median_arr.copy()
            continue
        # We are giving penalties here, so if gene is inactive we
        # should give max penalty, which is 1.
        imputation_arr = (
            np.ones(n_samples) if info["functional"] else median_arr
        )

        gene_arrays = {}
        for g in info["genes"]:
            row = gene_to_row.get(g)
            gene_arrays[g] = expr_values[row] if row is not None else imputation_arr

        results[rxn_id] = info["fn"](gene_arrays)

    return pd.DataFrame(results, index=expr_df.columns).T

# MetChange functions
################################################################################

TOL = 1e-8


class MetCHANGE_Gurobi:
    def __init__(self, S, lb, ub, threads=1, precomputed_max_fluxes=None):
        """
        S : metabolites x reactions. Pass this as a scipy.sparse matrix
            (csr_matrix recommended), NOT a dense numpy array. Your real
            model is ~8,461 rows x 37,312 columns but only ~142,000
            nonzero entries (>99.9% sparse) -- as dense float64 that's
            ~2.5GB PER WORKER PROCESS. That is almost certainly what
            caused runs to stall as worker count increased -- confirmed
            NOT the Gurobi token server (4096 tokens available, 2 in
            use) and unlikely to be core oversubscription at 24/28 cores;
            far more consistent with memory pressure pushing the node
            into swapping. As sparse, the same matrix is a few MB.

            Build it sparse from cobrapy like this:
                from cobra.util.array import create_stoichiometric_matrix
                S = create_stoichiometric_matrix(model, array_type="dok")
                S = S.tocsr()
        """
        if not sp.issparse(S):
            S = sp.csr_matrix(S)  # defensive fallback -- pass sparse in directly if you can

        self.n_met, self.n_rxn = S.shape

        # sparse identity block for the demand/sink augmentation
        demand_block = -sp.eye(self.n_met, format="csr")
        S_aug = sp.hstack([S, demand_block], format="csr")

        lb_dm = np.zeros(self.n_met)
        ub_dm = np.zeros(self.n_met)  # closed by default -- correctness fix from before

        self.lb = np.concatenate([lb, lb_dm])
        self.ub = np.concatenate([ub, ub_dm])
        self.n_var = S_aug.shape[1]
        self.dm_start = self.n_rxn

        self.model = gp.Model()
        self.model.Params.OutputFlag = 0
        self.model.Params.Threads = threads
        self.model.Params.Method = 1  # dual simplex -- warm-start friendly

        self.v = self.model.addMVar(shape=self.n_var, lb=self.lb, ub=self.ub, name="v")
        self.model.addMConstr(S_aug, self.v, sense="=", b=np.zeros(self.n_met))  # accepts sparse directly
        self.model.update()  # one-time build cost, not in the hot loop -- fine to keep

        # pre-fetch the underlying Var objects once for bulk Obj updates later
        self._rxn_vars = self.v[: self.n_rxn].tolist()

        self.max_fluxes = precomputed_max_fluxes  # None until computed or loaded

    # ---------------- LP1 (cacheable) ----------------

    def compute_max_fluxes(self):
        """Run LP1 for every metabolite once. Depends only on S/bounds,
        NOT on expression -- safe to cache and reuse across all conditions
        and across separate runs (as long as the model itself is unchanged)."""
        max_fluxes = np.zeros(self.n_met)
        for m in range(self.n_met):
            dm_idx = self.dm_start + m
            self.v[dm_idx].lb = 0
            self.v[dm_idx].ub = 1000

            obj = np.zeros(self.n_var)
            obj[dm_idx] = 1
            self.model.setObjective(obj @ self.v, GRB.MAXIMIZE)
            self.model.optimize()

            max_fluxes[m] = self.model.ObjVal if self.model.status == GRB.OPTIMAL else 0.0

            self.v[dm_idx].lb = 0
            self.v[dm_idx].ub = 0

        self.max_fluxes = max_fluxes

        # LP1's setObjective() calls leave whichever demand variable was
        # processed LAST with a stale nonzero Obj coefficient. score_metabolite()
        # only bulk-updates the *reaction* variables' Obj (that's the whole
        # optimization), so without this reset, a fixed demand flux later in
        # LP2 would silently add (stale_coeff * fixed_flux) to the reported
        # score for whichever metabolite happened to be processed last here.
        # One-time, cheap: zero every variable's Obj coefficient before LP2 begins.
        self.model.setAttr(GRB.Attr.Obj, self.v.tolist(), [0.0] * self.n_var)

        return max_fluxes

    def save_max_fluxes(self, path):
        if self.max_fluxes is None:
            raise RuntimeError("compute_max_fluxes() (or load) first")
        np.save(path, self.max_fluxes)

    def load_max_fluxes(self, path):
        self.max_fluxes = np.load(path)

    # ---------------- LP2 (run per condition) ----------------

    def score_metabolite(self, met_idx, expression):
        """
        expression: reactions x conditions (n_rxn x n_cond), already
                    transformed to penalties.
        Returns: array of length n_cond (NaN-filled if unproducible).
        """
        n_cond = expression.shape[1]
        row = np.full(n_cond, np.nan)

        max_flux = self.max_fluxes[met_idx]
        if max_flux <= TOL:
            return row # If a metabolite cannot be produced just return the nan row, as a metchange score of zero means the metabolite can be easily produced.
            #return np.zeros(n_cond)

        dm_idx = self.dm_start + met_idx
        # open + fix this metabolite's sink at its cached max, in one bound set
        self.v[dm_idx].lb = max_flux
        self.v[dm_idx].ub = max_flux

        # LP1 left ModelSense as MAXIMIZE; since we now update Obj coefficients
        # directly (no setObjective call here), the sense must be set explicitly.
        self.model.ModelSense = GRB.MINIMIZE

        for c in range(n_cond):
            penalty = expression[:, c]
            # bulk objective coefficient update -- no symbolic expression rebuild
            self.model.setAttr(GRB.Attr.Obj, self._rxn_vars, penalty.tolist())
            self.model.optimize()
            row[c] = self.model.ObjVal if self.model.status == GRB.OPTIMAL else np.nan

        # close the sink before moving to the next metabolite
        self.v[dm_idx].lb = 0
        self.v[dm_idx].ub = 0

        return row

# ==================== ACTUAL PARALLEL DRIVER ====================
#
# metchange_gurobi() above is single-process. Each metabolite's LP1+LP2
# block is independent of every other metabolite, so this splits the
# metabolite list into chunks and runs each chunk in its own process --
# each worker builds its OWN MetCHANGE_Gurobi instance (gurobipy Model
# objects can't be shared/pickled across processes) and loads the SAME
# cached max_fluxes (computed once, up front, before workers are spawned).

from multiprocessing import Pool


def _worker(args):
    S, lb, ub, met_indices, expression, max_fluxes_path, threads = args

    solver = MetCHANGE_Gurobi(S, lb, ub, threads=threads)
    solver.load_max_fluxes(max_fluxes_path)

    n_cond = expression.shape[1]
    out = np.zeros((len(met_indices), n_cond))
    for i, m in enumerate(met_indices):
        out[i, :] = solver.score_metabolite(m, expression)

    return met_indices, out


def metchange_gurobi_parallel(S, lb, ub, expression, metabolite_ids,
                               max_fluxes_path, n_workers=4, threads_per_worker=1):
    """
    Parallel version -- splits metabolites across `n_workers` processes,
    BALANCED BY ACTUAL WORKLOAD rather than raw index order: non-producible
    metabolites (max_flux ~ 0) skip LP2 entirely and cost almost nothing,
    so a naive contiguous split can leave one worker with most of the real
    (producible) work while others sit idle -- total time is bottlenecked
    by the slowest worker, not the average. Since max_fluxes is already
    known before parallelizing, producible and non-producible metabolites
    are each interleaved round-robin across workers, giving every worker
    close to the same number of "expensive" metabolites to solve.

    max_fluxes_path : REQUIRED here. Compute it once, sequentially, before
                       calling this (see ensure_max_fluxes below).
    n_workers        : keep <= your available Gurobi license concurrency
                       if using a floating/token license.

    Returns: DataFrame, metabolites x conditions (row order = metabolite_ids)
    """
    n_met = S.shape[0]
    max_fluxes = np.load(max_fluxes_path)

    producible = np.where(max_fluxes > TOL)[0]
    non_producible = np.where(max_fluxes <= TOL)[0]

    # round-robin (not contiguous-block) assignment balances workload
    chunks = [[] for _ in range(n_workers)]
    for i, m in enumerate(producible):
        chunks[i % n_workers].append(m)
    for i, m in enumerate(non_producible):
        chunks[i % n_workers].append(m)
    chunks = [np.array(sorted(c), dtype=int) for c in chunks if len(c) > 0]

    args = [
        (S, lb, ub, chunk, expression, max_fluxes_path, threads_per_worker)
        for chunk in chunks
    ]

    n_cond = expression.shape[1]
    results = np.zeros((n_met, n_cond))

    with Pool(len(args)) as pool:
        for met_indices, out in pool.imap_unordered(_worker, args):
            results[met_indices, :] = out

    return pd.DataFrame(results, index=metabolite_ids)


def ensure_max_fluxes(S, lb, ub, path, threads=1):
    """Compute LP1 once and save it, or skip if the file already exists.
    Call this once before metchange_gurobi_parallel()."""
    import os
    if os.path.exists(path):
        return path
    solver = MetCHANGE_Gurobi(S, lb, ub, threads=threads)
    solver.compute_max_fluxes()
    solver.save_max_fluxes(path)
    return path

# Functions for obtaining input matrices
################################################################################

# Option to ingore compartments, as x (boundary compartment) can produce
# FBA without biomass
# https://sysbiochalmers.github.io/Human-GEM-guide/faq_troubleshoot/
# https://github.com/SysBioChalmers/Human-GEM/issues/172

def build_split_stoichiometry(model, ignore_compartments=("x",)):
    S = create_stoichiometric_matrix(model, array_type="dok").tocoo()
    rows, cols, data = S.row, S.col, S.data

    # --- drop boundary-compartment metabolites from the mass-balance rows ---
    keep_idx = np.array(
        [i for i, m in enumerate(model.metabolites) if m.compartment not in ignore_compartments]
    )
    row_remap = np.full(len(model.metabolites), -1, dtype=np.int64)
    row_remap[keep_idx] = np.arange(len(keep_idx))

    mask = row_remap[rows] >= 0
    rows = row_remap[rows[mask]]
    cols = cols[mask]
    data = data[mask]
    n_mets = len(keep_idx)
    # --------------------------------------------------------------------

    met_ids = np.array([model.metabolites[i].id for i in keep_idx])
    
    rxn_ids = np.array([r.id for r in model.reactions])
    lb = np.array([r.lower_bound for r in model.reactions])
    ub = np.array([r.upper_bound for r in model.reactions])
    n_orig = len(rxn_ids)
    #n_mets = S.shape[0]

    reversible_mask = (lb < 0) & (ub > 0)
    flip_mask = ub <= 0                     # reverse-only -> flip sign
    fwd_sign = np.where(flip_mask, -1, 1)   # sign applied to the "first" slot

    # Column layout: every original reaction gets 1 slot ("first"); a
    # reversible reaction additionally gets a 2nd ("reverse") slot right
    # after it.
    n_slots = np.where(reversible_mask, 2, 1)
    first_slot = np.zeros(n_orig, dtype=np.int64)
    first_slot[1:] = np.cumsum(n_slots)[:-1]
    rev_slot = first_slot + 1               # only meaningful where reversible
    total_cols = int(first_slot[-1] + n_slots[-1]) if n_orig else 0

    # First-slot entries: every original nonzero, remapped to its new
    # column, sign-flipped only for reverse-only (flip_mask) reactions.
    rows_first = rows
    cols_first = first_slot[cols]
    data_first = data * fwd_sign[cols]

    # Reverse-slot entries: only for reversible reactions, sign-flipped.
    is_rev_nz = reversible_mask[cols]
    rows_rev = rows[is_rev_nz]
    cols_rev = rev_slot[cols[is_rev_nz]]
    data_rev = -data[is_rev_nz]

    all_rows = np.concatenate([rows_first, rows_rev])
    all_cols = np.concatenate([cols_first, cols_rev])
    all_data = np.concatenate([data_first, data_rev])

    S_split = coo_matrix(
        (all_data, (all_rows, all_cols)), shape=(n_mets, total_cols)
    ).tocsc()

    # Bounds / bookkeeping per new column, vectorized.
    col_rxn_id = np.empty(total_cols, dtype=object)
    col_sign = np.empty(total_cols, dtype=np.int8)
    lbs = np.empty(total_cols)
    ubs = np.empty(total_cols)

    col_rxn_id[first_slot] = rxn_ids
    col_sign[first_slot] = fwd_sign
    lbs[first_slot] = np.where(flip_mask, -ub, np.where(reversible_mask, 0.0, lb))
    ubs[first_slot] = np.where(flip_mask, -lb, np.where(reversible_mask, ub, ub))

    rev_idx = rev_slot[reversible_mask]
    col_rxn_id[rev_idx] = rxn_ids[reversible_mask]
    col_sign[rev_idx] = -1
    lbs[rev_idx] = 0.0
    ubs[rev_idx] = -lb[reversible_mask]

    return S_split, lbs, ubs, col_rxn_id, col_sign, list(rxn_ids), met_ids

def expand_expression(g_full, rxn_ids, col_rxn_id):
    idx = {rid: i for i, rid in enumerate(rxn_ids)}
    g_full = np.asarray(g_full)
    return np.array([g_full[idx[rid]] for rid in col_rxn_id])


################################################################################
# Terminal argument parser                                                     #
################################################################################
parser = argparse.ArgumentParser(
    description='Run MetChange.'
)

parser.add_argument(
        '--expr',
        type=str,
        help='Expression file.'
)

parser.add_argument(
        '--input_type',
        type=str,
        help='Whether the input expression file is just expression or penalty. If expression selected, penalties will be computed based on ranks sample-wise. If penalty selected, this step will be skipped and the matrix will be directly fed into MetChange.',
        default="expression"
)

parser.add_argument(
    "-s", "--save_penalties", 
    action="store_true", 
    help="Whether penalties should be saved (only effective when --input_type == 'expression')"
)

parser.add_argument(
        '--medium',
        type=str,
        help='Medium file.',
)

parser.add_argument(
        '--threshold',
        type=int,
        help='Threshold percentile to considere a gene not expressed at all.',
        default=2
)

parser.add_argument(
        '--gem',
        type=str,
        help='GEM file.'
)

parser.add_argument(
        '--n_workers',
        type=int,
        help='Number of parallel workers.',
        default=1
)

parser.add_argument('--out_dir', 
	                  type=str,
                    help='Directory where the result will be stored.', 
                    default=os.getcwd())

args = parser.parse_args()

################################################################################
# DIRECTORY STUFF                                                              #
################################################################################

expr_f = args.expr
input_type = args.input_type
medium_f = args.medium
threshold = args.threshold
gem_f = args.gem
n_workers = args.n_workers
out_dir = args.out_dir
save_penalties = args.save_penalties

os.makedirs(out_dir, exist_ok=True)

################################################################################
# DATA LOADING                                                                 #
################################################################################

print("Loading GEM: %s" % gem_f)
gem = read_sbml_model(gem_f)

print("Loading expression data from %s" % expr_f)
lincs_expr = pd.read_csv(expr_f, index_col = 0)

if medium_f != None:
    with open(medium_f) as f:
        medium = json.load(f)
    gem.medium = medium

################################################################################
# PARSE INPUT DATA                                                             #
################################################################################

if input_type == "expression":
    print("Input type: expression. Obtaining rank-based gene penalties and generating reaction matrix with them...")

    # Obtan gene-level penalties
    pval_df = []
    for i in range(0, lincs_expr.shape[1]):
        p_vals = 1 - lincs_expr.iloc[:, i].rank() / lincs_expr.shape[0]
        perc = np.percentile(p_vals, 100 - threshold)
        p_vals.loc[p_vals >= perc] = 1
        pval_df.append(p_vals)
    pval_df = pd.concat(pval_df, axis = 1)
    if save_penalties:
        pval_df.to_csv(os.path.join(out_dir, "rank_penalties.csv"))
    # Convert gene-level penalties to reanction-level penalties
    g = get_g_matrix(gem, pval_df)
elif input_type == "penalty":
    print("Input type: penalty. Generating reaction matrix directly with them...")
    g = get_g_matrix(gem, lincs_expr)

# Obtain matrices and vectors for running metChange
S, lb, ub, col_rxn_id, col_sign, rxn_ids, met_ids = build_split_stoichiometry(gem)
expr = expand_expression(g, rxn_ids, col_rxn_id)

################################################################################
# RUN METCHANGE                                                                #
################################################################################

# Compute maximum possible metabolite fluxes given the model
start_time = time.time()
print ("Computing LP1...")
ensure_max_fluxes(S, lb, ub, path=os.path.join(out_dir, "lp1_fluxes.npy"), threads=10)
print("LP1 finished: %s seconds elapsed" % (time.time() - start_time))


# Compute condition-specific maximum metabolite fluxes

start_time = time.time()
print ("Computing LP2...")
result = metchange_gurobi_parallel(
    S, lb, ub, expr, met_ids,
    "lp1_fluxes.npy", n_workers=n_workers, threads_per_worker=1
)
print("LP2 finished: %s seconds elapsed" % (time.time() - start_time))
result.columns = g.columns


out_name = os.path.join(
    out_dir,
    "%s_metChange_%s.csv" % (
        os.path.basename(expr_f).split('.')[0],
        os.path.basename(medium_f).split('.')[0],
    )
)
result.to_csv(out_name)

print("%s saved at %s" % (os.path.basename(out_name), out_dir))