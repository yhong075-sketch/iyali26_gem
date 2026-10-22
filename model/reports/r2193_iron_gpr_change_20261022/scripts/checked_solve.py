"""Knockout solves that cannot silently reuse a stale basis (shared by the method check and the sweep).

reset_basis() gives GLPK a fresh advanced starting basis before each solve (as a newly built problem
gets), so a knockout is not solved from the previous optimum. check() measures, for the returned
point, the largest flux on the reactions the knockout closed, the largest bound violation of any
reaction flux or split variable, and the largest |S v|.

Targeted check (user decision, 2026-10-22): a point fails only if a knocked-out reaction carries
|v| > KNOCKED_TOLERANCE (the stale-point signature) or a general violation exceeds GENERAL_TOLERANCE
(about GLPK's own accuracy). The first version failed any violation above 1e-9, which flagged
GLPK tolerance-level violations (../flag_diagnosis.json); method_check.json was produced by that first
version.

Known blind spot (independent audit, 2026-10-22): GLPK can drop a tiny biomass demand within its
tolerance. Most |S v| residuals in the sweep equal 1e-6 x growth on a heme-route row, i.e. the heme a
term skipped. That residual (at most about 3.1e-7) is below GENERAL_TOLERANCE, so this check would not
catch a heme-step knockout that GLPK "rescues" this way. In the 2026-10-22 sweep every heme-route
knockout returned 0, and Gurobi at 1e-9 reproduced the flagged growth values, so no call changed. A
residual tolerance scaled to growth x the smallest biomass coefficient, or a Gurobi cross-check of
non-zero residuals, would close it (untested).
"""
import math

import numpy as np
from cobra.util.array import create_stoichiometric_matrix
from swiglpk import glp_adv_basis

KNOCKED_TOLERANCE = 1e-9
GENERAL_TOLERANCE = 1e-6


def reset_basis(model):
    glp_adv_basis(model.solver.problem, 0)


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


def solve(model, stoichiometry, knocked, reset=True):
    # knocked has no default: an empty set would silently turn off the stale-point test.
    if reset:
        reset_basis(model)
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
