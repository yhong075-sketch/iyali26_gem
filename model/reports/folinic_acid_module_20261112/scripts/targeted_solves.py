"""WT, module flux ranges and targeted knockouts on one model (read only; no model file is written).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts \
    /opt/miniconda3/bin/python targeted_solves.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json

Same loader, medium, PO1f profile and objective as the sweep; Gurobi at 1e-9 with a reset before every solve and
the targeted check (checked_solve.py, unchanged). Flux ranges are direct LPs, min and max of each reaction's net
flux, with biomass_C held at >= FLOOR x WT. The 100 % floor uses (1 - 1e-9) x WT so the LP stays feasible
within the solver tolerance; the 0 % floor allows any growth, including none.
"""
import json
import math
import sys
import warnings
from pathlib import Path

from checked_solve import configure, matrix, reset, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

MODULE = ("R748", "R989", "R64", "R990", "R811", "R991", "R_THFAT", "R543", "R542")
KNOCKOUTS = (("YALI1F04082g",), ("YALI1D28569g",), ("YALI1E19557g",), ("YALI1D28569g", "YALI1E19557g"),
             ("YALI1D01489g",))
FLOORS = {"100pct": 1.0 - 1e-9, "0pct": 0.0}
BIOMASS = "biomass_C"

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


def ranges(wt_growth):
    out = {}
    biomass = model.reactions.get_by_id(BIOMASS)
    for floor_name, floor in FLOORS.items():
        out[floor_name] = {}
        for rid in MODULE:
            if rid not in model.reactions:
                out[floor_name][rid] = None
                continue
            row = {}
            for sense in ("min", "max"):
                with model:
                    biomass.lower_bound = floor * wt_growth
                    model.objective = model.reactions.get_by_id(rid)
                    model.objective_direction = sense
                    reset(model)
                    value = model.slim_optimize(error_value=float("nan"))
                    status = model.solver.status
                    row[sense] = value if status == "optimal" and not math.isnan(value) else None
                    row[f"{sense}_status"] = status
                    row[f"{sense}_growth"] = net(BIOMASS) if status == "optimal" else None
            out[floor_name][rid] = row
    return out


def closed(genes):
    return sorted({r.id for g in genes for r in model.genes.get_by_id(g).reactions if r.bounds == (0, 0)})


wt = solve(model, stoichiometry, ())
assert wt["passed"], wt
wt["module_flux"] = {rid: net(rid) for rid in MODULE if rid in model.reactions}
result = {"model_sha256": ctx.canonical_model_sha256, "wt": wt, "ranges": ranges(wt["growth"]), "knockouts": {}}
for genes in KNOCKOUTS:
    with model:
        for g in genes:
            model.genes.get_by_id(g).knock_out()
        knocked = closed(genes)
        row = solve(model, stoichiometry, set(knocked))
        row["closed_reactions"] = knocked
        row["ko_wt"] = row["growth"] / wt["growth"] if row["passed"] else None
    result["knockouts"]["+".join(genes)] = row
Path(out_json).write_text(json.dumps(result, indent=1))
print("model", result["model_sha256"], "WT", wt["growth"], wt["status"])
for floor_name, rows in result["ranges"].items():
    for rid, row in rows.items():
        print(floor_name, rid, None if row is None else (row["min"], row["max"], row["min_status"], row["max_status"]))
for label, row in result["knockouts"].items():
    print("KO", label, row["status"], row["growth"], row["ko_wt"], row["passed"], row["closed_reactions"])
