################################################################################
####################----------------------------------------####################
##########----------            MetaboReprogram             ----------##########
# Run metChange LP1 step: intended to be used when running array jobs, to run  #
# this step just once before launching the array, as it doesn't need the       #
# expression data.                                                             #
#                                                                              #
################################################################################

################################################################################
# MODULE LOADING                                                               #
################################################################################

import os
import numpy as np
import time
import scipy.sparse as sp
from scipy.sparse import coo_matrix
import json
import argparse

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
#import gurobipy
import cobra
from cobra.io import read_sbml_model, write_sbml_model
from cobra import Reaction
from cobra.util.solver import solvers
from cobra.util.array import create_stoichiometric_matrix
print(solvers.keys())

################################################################################
# FUNCTIONS                                                                    #
################################################################################

# FUnctions for running MetCHANGE LP1

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
            return np.zeros(n_cond)

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

################################################################################
# Terminal argument parser                                                     #
################################################################################
parser = argparse.ArgumentParser(
    description='Run MetChange.'
)

parser.add_argument(
        '--medium',
        type=str,
        help='Medium file.',
)

parser.add_argument(
        '--gem',
        type=str,
        help='GEM file.'
)

parser.add_argument('--out_dir', 
	                  type=str,
                    help='Directory where the result will be stored.', 
                    default=os.getcwd())

args = parser.parse_args()

################################################################################
# DIRECTORY STUFF                                                              #
################################################################################

medium_f = args.medium
gem_f = args.gem
out_dir = args.out_dir

os.makedirs(out_dir, exist_ok=True)


################################################################################
# DATA LOADING                                                                 #
################################################################################
gem = read_sbml_model(gem_f)

if medium_f != None:
    with open(medium_f) as f:
        medium = json.load(f)
    gem.medium = medium

################################################################################
# PREPARE INPUT DATA                                                                 #
################################################################################

S, lb, ub, col_rxn_id, col_sign, rxn_ids, met_ids = build_split_stoichiometry(gem)

################################################################################
# RUN METCHANGE LP1 STEP                                                       #
################################################################################

# Compute maximum possible metabolite fluxes given the model
start_time = time.time()
print ("Computing global LP1...")
ensure_max_fluxes(S, lb, ub, path=os.path.join(out_dir, "lp1_fluxes.npy"), threads=10)
print("LP1 finished: %s seconds elapsed" % (time.time() - start_time))