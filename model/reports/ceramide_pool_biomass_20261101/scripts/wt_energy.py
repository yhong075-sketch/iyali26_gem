"""WT growth and the no-uptake energy check on the baseline and the candidate (TASK.md §7 items 4-5; 4 LPs).

Energy check: every boundary reaction's uptake closed (lower bound 0 where it was negative), biomass_C and
xMAINTENANCE lower bounds 0, then xMAINTENANCE (ATP hydrolysis) maximised; a value above 0 would be an energy cycle.
Task A medium + PO1f, checked Gurobi solves at 1e-9.
Usage (from platform/): PYTHONPATH=$PWD:<checked_solve dir> python wt_energy.py MEDIUM PROFILE OUT.json MODEL.xml...
"""
import json
import sys
import warnings
from pathlib import Path

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

medium, profile, out, *models = sys.argv[1:]
result, lps = {}, 0
for path in models:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=Path(path), media_path=medium, strain_profile_path=profile)
    model = ctx.model
    configure(model)
    S = matrix(model)
    wt = solve(model, S, ())
    lps += 1
    with model:
        closed = []
        for reaction in model.boundary:
            if reaction.lower_bound < 0:
                reaction.lower_bound = 0
                closed.append(reaction.id)
        model.reactions.biomass_C.lower_bound = 0
        model.reactions.xMAINTENANCE.lower_bound = 0
        model.objective = "xMAINTENANCE"
        energy = solve(model, S, ())
        lps += 1
    result[Path(path).name] = {"wt": wt, "energy_no_uptake": energy, "uptakes_closed": len(closed)}
    print(Path(path).name, "WT", wt["status"], wt["growth"], wt["passed"], "| ATP max no uptake", energy["status"],
          energy["growth"], energy["passed"], "| uptakes closed", len(closed))
result["lp_count"] = lps
Path(out).write_text(json.dumps(result, indent=1, default=str) + "\n")
