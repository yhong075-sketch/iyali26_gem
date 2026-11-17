"""Largest flux (min and max, no growth floor) of the other reactions drawn in the schematic, on both built models.

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts \
    /opt/miniconda3/bin/python figure_capacity.py CONTROL.xml CANDIDATE.xml MEDIUM.csv PROFILE.json OUT.json
Same loader and solver settings as targeted_solves.py (which covers the nine and R_OBDHm).
"""
import json
import math
import sys
import warnings
from pathlib import Path

from checked_solve import configure, reset
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

DRAWN = ("R963", "R515", "R693", "R83", "R412", "R1172", "R410", "R2113", "R2285", "R741", "R2025", "R2024")
control, candidate, medium, profile, out_json = sys.argv[1:6]
LPS = {"n": 0}


def capacity(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=Path(path), media_path=medium, strain_profile_path=profile)
    model = ctx.model
    configure(model)
    out = {"model_sha256": ctx.canonical_model_sha256, "fva": {}}
    for rid in DRAWN:
        row = {}
        for sense in ("min", "max"):
            with model:
                model.objective = model.reactions.get_by_id(rid)
                model.objective_direction = sense
                reset(model)
                LPS["n"] += 1
                value = model.slim_optimize(error_value=float("nan"))
                status = model.solver.status
                row[sense] = value if status == "optimal" and not math.isnan(value) else None
                row[f"{sense}_status"] = status
        out["fva"][rid] = row
    return out


result = {"control": capacity(control), "candidate": capacity(candidate)}
result["lp_count"] = LPS["n"]
Path(out_json).write_text(json.dumps(result, indent=1))
for tag in ("control", "candidate"):
    print(tag, result[tag]["model_sha256"][:8],
          {r: (v["min"], v["max"], v["min_status"], v["max_status"]) for r, v in result[tag]["fva"].items()})
print("LPs", LPS["n"])
