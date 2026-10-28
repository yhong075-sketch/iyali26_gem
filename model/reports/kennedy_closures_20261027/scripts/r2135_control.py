"""Control for the acyl leak: R2135 (acetyl-CoA -> generic ER acyl-CoA) closed, baseline vs candidate (read-only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/kennedy_closures_20261027/scripts/r2135_control.py BASELINE.xml CANDIDATE.xml MEDIUM.csv PROFILE.json OUT.json
Task A medium + PO1f, Gurobi 1e-9, checked solves: WT and the four changed knockouts, with R2135 as built and
closed; FVA of R2135 at 100 % / 99.9 % of the WT optimum (processes=1). Diagnostic only, not a proposed change.
"""
import json
import sys
import warnings

from cobra.flux_analysis import flux_variability_analysis

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

GENES = ("YALI1B19844g", "YALI1B25427g", "YALI1C07638g", "YALI1D22124g")
baseline, candidate, medium, profile, out = sys.argv[1:6]
result = {}
for label, path in (("baseline", baseline), ("candidate", candidate)):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile).model
    configure(model)
    s = matrix(model)
    frame = flux_variability_analysis(model, reaction_list=["R2135"], fraction_of_optimum=0.999, processes=1)
    rows = {"fva_R2135_0.999": [float(frame.loc["R2135", "minimum"]), float(frame.loc["R2135", "maximum"])]}
    for closed in (False, True):
        with model:
            if closed:
                model.reactions.R2135.bounds = (0.0, 0.0)
            tag = "R2135_closed" if closed else "as_built"
            rows[f"{tag}_wt"] = solve(model, s, {"R2135"} if closed else set())
            for gene in GENES:
                with model:
                    model.genes.get_by_id(gene).knock_out()
                    knocked = {r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0)}
                    rows[f"{tag}_{gene}"] = solve(model, s, knocked | ({"R2135"} if closed else set()))
    result[label] = rows
    print(label, rows["fva_R2135_0.999"], {k: (round(v["growth"], 7) if v["growth"] is not None else None)
                                            for k, v in rows.items() if k != "fva_R2135_0.999"}, file=sys.stderr, flush=True)
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
