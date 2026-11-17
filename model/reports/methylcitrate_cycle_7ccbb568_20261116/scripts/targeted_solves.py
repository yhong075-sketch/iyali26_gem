"""WT, flux ranges of the nine and R_OBDHm, energy check and PDH-gene knockouts on one model (read only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts \
    /opt/miniconda3/bin/python targeted_solves.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json

Same loader, medium, PO1f profile and objective as the sweep; Gurobi at 1e-9 with a reset before every solve and
the targeted check (checked_solve.py, unchanged). Flux ranges are direct LPs, min and max of each reaction's net
flux, with no growth floor (biomass_C free, as in the investigation's FVA). Energy check as in the investigation
(../../r1172_side_effect_20261116/scripts/step5_options.py): every boundary lower bound raised to at least 0 and every
positive lower bound (xMAINTENANCE) set to 0; then maximum R_NTP1 and maximum biomass_C.
"""
import json
import math
import sys
import warnings
from pathlib import Path

from checked_solve import configure, matrix, reset, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

MODULE = ("R25", "R552", "R95", "R2074", "R2075", "R2076", "R2077", "R2112", "R834", "R_OBDHm")
KNOCKOUTS = ("YALI1E31947g", "YALI1F27556g", "YALI1D26360g", "YALI1B13179g", "YALI1D30618g")
BIOMASS = "biomass_C"
LPS = {"n": 0}

model_path, medium, profile, out_json = sys.argv[1:5]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=Path(model_path), media_path=medium,
                                            strain_profile_path=profile)
model = ctx.model
configure(model)
stoichiometry = matrix(model)
assert [r.id for r in model.reactions if r.objective_coefficient] == [BIOMASS]


def net(rid):
    primal = model.solver.primal_values
    r = model.reactions.get_by_id(rid)
    return primal[r.forward_variable.name] - primal[r.reverse_variable.name]


def optimum(objective, sense="max"):
    """One LP; the value only when optimal (a non-optimal status is never read as zero)."""
    with model:
        model.objective = model.reactions.get_by_id(objective)
        model.objective_direction = sense
        reset(model)
        LPS["n"] += 1
        value = model.slim_optimize(error_value=float("nan"))
        status = model.solver.status
        ok = status == "optimal" and not math.isnan(value)
        return {"value": value if ok else None, "status": status,
                "growth": net(BIOMASS) if ok else None}


def ranges():
    out = {}
    for rid in MODULE:
        if rid not in model.reactions:
            out[rid] = None
            continue
        lo, hi = optimum(rid, "min"), optimum(rid, "max")
        out[rid] = {"min": lo["value"], "max": hi["value"], "min_status": lo["status"], "max_status": hi["status"]}
    return out


def energy_check():
    out = {}
    with model:
        for r in model.boundary:
            r.lower_bound = max(r.lower_bound, 0.0)
        relaxed = sorted(r.id for r in model.reactions if r.lower_bound > 0)
        for r in model.reactions:
            if r.lower_bound > 0:
                r.lower_bound = 0.0
        for obj in ("R_NTP1", BIOMASS):
            row = optimum(obj, "max")
            out[obj] = {"value": row["value"], "status": row["status"]}
        out["positive_lower_bounds_set_to_0"] = relaxed
    return out


LPS["n"] += 1
wt = solve(model, stoichiometry, ())
assert wt["passed"], wt
wt["module_flux"] = {rid: net(rid) for rid in MODULE if rid in model.reactions}
result = {"model_sha256": ctx.canonical_model_sha256, "wt": wt, "ranges_no_floor": ranges(),
          "energy_check": energy_check(), "knockouts": {}}
for gene in KNOCKOUTS:
    with model:
        model.genes.get_by_id(gene).knock_out()
        knocked = sorted(r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0))
        LPS["n"] += 1
        row = solve(model, stoichiometry, set(knocked))
        row["closed_reactions"] = knocked
        row["ko_wt"] = row["growth"] / wt["growth"] if row["passed"] else None
    result["knockouts"][gene] = row
result["lp_count"] = LPS["n"]
Path(out_json).write_text(json.dumps(result, indent=1))
print("model", result["model_sha256"], "WT", repr(wt["growth"]), wt["status"], "LPs", LPS["n"])
for rid, row in result["ranges_no_floor"].items():
    print(rid, None if row is None else (row["min"], row["max"], row["min_status"], row["max_status"]))
print("energy", result["energy_check"])
for label, row in result["knockouts"].items():
    print("KO", label, row["status"], row["growth"], row["ko_wt"], row["passed"], row["closed_reactions"])
