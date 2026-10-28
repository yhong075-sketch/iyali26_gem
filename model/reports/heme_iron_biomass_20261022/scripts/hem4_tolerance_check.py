"""Why the sweep reported WT growth for the HEM4 (YALI1F13039g) knockout on the new model. Read-only, 4 LPs.

Usage (from platform/): PYTHONPATH=$PWD python ../model/reports/heme_iron_biomass_20261022/scripts/hem4_tolerance_check.py \
    NEW.xml MEDIUM.csv PROFILE.json OUT.json
Cases: GLPK warm start (WT solve, then the knockout in the same context, as in the sweep); Gurobi
with its default FeasibilityTol (1e-6) and with 1e-9. For each optimal point, the largest
steady-state residual |S v| and the residual on the heme a [C_cy] and heme rows are recorded.
"""
import json
import sys
import warnings

import numpy as np
from cobra.util.array import create_stoichiometric_matrix

from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

GENE = "YALI1F13039g"
ROWS = ["m1635[C_cy]", "m609[C_mi]", "m525[C_mi]", "m645[C_cy]"]
new_path, medium, profile, out = sys.argv[1:5]


def context(solver, tolerance=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=new_path, media_path=medium, strain_profile_path=profile)
    ctx.model.solver = solver
    if tolerance is not None:
        ctx.model.solver.configuration.tolerances.feasibility = tolerance
    return ctx.model


def record(model, label):
    solution = model.optimize()
    flux = solution.fluxes.to_numpy(dtype=float)
    matrix = create_stoichiometric_matrix(model, array_type="dense")
    residual = matrix @ flux
    index = {m.id: i for i, m in enumerate(model.metabolites)}
    return {"case": label, "status": solution.status, "growth": float(solution.fluxes["biomass_C"]),
            "R792": float(solution.fluxes["R792"]), "R1709": float(solution.fluxes["R1709"]),
            "R302": float(solution.fluxes["R302"]),
            "feasibility_tolerance": model.solver.configuration.tolerances.feasibility,
            "max_abs_residual": float(np.max(np.abs(residual))),
            "row_residuals": {mid: float(residual[index[mid]]) for mid in ROWS}}


cases = []
model = context("glpk")
model.slim_optimize()
with model:
    model.genes.get_by_id(GENE).knock_out()
    cases.append(record(model, "GLPK warm start: WT solve, then knockout in the same context"))
for tolerance in (None, 1e-9):
    model = context("gurobi", tolerance)
    model.genes.get_by_id(GENE).knock_out()
    cases.append(record(model, f"Gurobi fresh, feasibility tolerance {'default' if tolerance is None else tolerance}"))
result = {"condition": "Task A medium + PO1f", "gene": GENE, "cases": cases, "lp_count": 4}
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
print(json.dumps(result, indent=1))
