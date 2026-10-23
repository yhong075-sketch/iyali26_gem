"""Checked LP solves on Gurobi at tight tolerances (shared by the WT solves and the knockout sweep).

Adapted from ../../r2193_iron_gpr_change_20261022/scripts/checked_solve.py (GLPK). Here the solver is
Gurobi with FeasibilityTol and OptimalityTol 1e-9 (user decision, 2026-10-22). reset() discards
Gurobi's previous solution and basis before each solve, so a knockout is not started from the previous
optimum. check() measures, for the returned point, the largest flux on the reactions the knockout
closed, the largest bound violation of any reaction flux or split variable, and the largest |S v|.

Targeted check (user decision, 2026-10-22): a point fails only if a knocked-out reaction carries
|v| > KNOCKED_TOLERANCE or a general violation exceeds GENERAL_TOLERANCE. All magnitudes are recorded,
so points with residuals near growth x a small biomass coefficient (the GLPK heme-skip pattern) remain
visible even though they pass.
"""
import math

import numpy as np
from cobra.util.array import create_stoichiometric_matrix

KNOCKED_TOLERANCE = 1e-9
GENERAL_TOLERANCE = 1e-6
SOLVER_TOLERANCE = 1e-9
LP_TIMEOUT_SECONDS = 60


def configure(model):
    model.solver = "gurobi"
    model.solver.configuration.tolerances.feasibility = SOLVER_TOLERANCE
    model.solver.configuration.tolerances.optimality = SOLVER_TOLERANCE
    model.solver.configuration.timeout = LP_TIMEOUT_SECONDS


def reset(model):
    model.solver.problem.reset()


def matrix(model):
    return create_stoichiometric_matrix(model, array_type="lil").tocsr()


def check(model, stoichiometry, knocked):
    primal = model.solver.primal_values
    flux, violation, knocked_flux = [], 0.0, 0.0
    for reaction in model.reactions:
        forward, reverse = primal[reaction.forward_variable.name], primal[reaction.reverse_variable.name]
        value = forward - reverse
        flux.append(value)
        low, high = reaction.bounds
        violation = max(violation, low - value, value - high, -forward, -reverse)
        if reaction.id in knocked:
            knocked_flux = max(knocked_flux, abs(value))
    flux = np.array(flux)
    residual = float(np.max(np.abs(stoichiometry @ flux))) if np.isfinite(flux).all() else math.inf
    return max(violation, 0.0), residual, knocked_flux


def solve(model, stoichiometry, knocked):
    # knocked has no default: an empty set would silently turn off the stale-point test.
    reset(model)
    growth = model.slim_optimize(error_value=float("nan"))
    status = model.solver.status
    if status != "optimal" or math.isnan(growth):
        return {"status": status, "growth": None, "max_bound_violation": None, "max_residual": None,
                "max_knocked_flux": None, "passed": False}
    violation, residual, knocked_flux = check(model, stoichiometry, knocked)
    return {"status": status, "growth": growth, "max_bound_violation": violation, "max_residual": residual,
            "max_knocked_flux": knocked_flux,
            "passed": knocked_flux <= KNOCKED_TOLERANCE and violation <= GENERAL_TOLERANCE
                      and residual <= GENERAL_TOLERANCE}
