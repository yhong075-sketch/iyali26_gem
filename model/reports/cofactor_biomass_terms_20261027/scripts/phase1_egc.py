"""Energy-generating-cycle check, phase 1 addendum (read-only, 2 LPs + diagnostics).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    phase1_egc.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json
The census closure_check method (max R_NTP1 with every exchange at (0, 1000), R_NTP1 lower bound 0) is
infeasible on this candidate under the Task A + PO1f context: it lists every reaction with a positive lower
bound, then repeats the LP with those lower bounds at 0, on the baseline and on the six-term model.
"""
import json
import sys
import warnings

from checked_solve import configure
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

TERMS = {"m479[C_cy]": 6.34e-5, "m528[C_cy]": 1e-5, "m926[C_cy]": 9.9e-4, "m122[C_cy]": 2.65e-3, "m41[C_cy]": 5.7e-4, "m69[C_cy]": 1.9e-4}
model_path, medium, profile, out = sys.argv[1:5]


def context():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
    configure(ctx.model)
    return ctx.model


def egc(model, relax):
    with model:
        for r in model.exchanges:
            r.bounds = (0.0, 1000.0)
        relaxed = []
        for r in model.reactions:
            if r.lower_bound > 0 and (relax or r.id == "R_NTP1"):
                relaxed.append([r.id, r.lower_bound])
                r.lower_bound = 0.0
        model.objective = model.reactions.get_by_id("R_NTP1")
        value = model.slim_optimize(error_value=float("nan"))
        return {"status": model.solver.status, "max_R_NTP1_no_uptake": value, "relaxed_positive_lower_bounds": relaxed}


result = {}
for label, mids in (("baseline", {}), ("six_terms", TERMS)):
    model = context()
    model.reactions.biomass_C.add_metabolites({model.metabolites.get_by_id(k): -v for k, v in mids.items()})
    result[label] = {"positive_lower_bounds": [[r.id, r.name[:40], r.lower_bound] for r in model.reactions if r.lower_bound > 0 and not r.boundary],
                     "census_method": egc(model, relax=False), "positive_lbs_relaxed": egc(model, relax=True)}
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
for label, row in result.items():
    print(label, "positive lbs:", row["positive_lower_bounds"], "| census method:", row["census_method"]["status"], row["census_method"]["max_R_NTP1_no_uptake"],
          "| relaxed:", row["positive_lbs_relaxed"]["status"], row["positive_lbs_relaxed"]["max_R_NTP1_no_uptake"])
