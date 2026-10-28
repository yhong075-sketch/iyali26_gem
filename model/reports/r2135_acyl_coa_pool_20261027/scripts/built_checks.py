"""Checks on the built files: baseline 00e8e397 vs the R2135 candidate (read-only, about 40 LPs).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/r2135_acyl_coa_pool_20261027/scripts/built_checks.py BASELINE.xml CANDIDATE.xml MEDIUM.csv PROFILE.json OUT.json
Task A medium + PO1f, Gurobi 1e-9 (checked_solve.py). Per model, as built and with R1726 opened to (0, 1000) in
memory (the paused Kennedy change): WT (checked); the equations of R2135 / R2145; pFBA fluxes of the TAG / PA
route and the phosphate exchange; FVA at 99.9 % of R2135, R2168, R2145, R2134, R1726 (processes=1); no-uptake
LPs (every exchange (0, 1000), positive lower bounds relaxed; max R_NTP1 and max biomass_C, status recorded);
the knockouts of the four genes that lost essentiality when R1726 was opened on 00e8e397.
"""
import hashlib
import json
import math
import sys
import warnings
from pathlib import Path

from cobra.flux_analysis import flux_variability_analysis, pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

GENES = ("YALI1B19844g", "YALI1B25427g", "YALI1C07638g", "YALI1D22124g")
FVA = ["R2135", "R2168", "R2145", "R2134", "R1726"]
ROUTE = FVA + ["R2144", "R2169", "R2170", "R1771", "R1296", "R1842"]
baseline_path, candidate_path, medium, profile, out = sys.argv[1:6]


def no_uptake(model, objective):
    with model:
        for r in model.exchanges:
            r.bounds = (0.0, 1000.0)
        relaxed = [[r.id, r.lower_bound] for r in model.reactions if r.lower_bound > 0]
        for rid, _ in relaxed:
            model.reactions.get_by_id(rid).lower_bound = 0.0
        model.objective = objective
        model.solver.problem.reset()
        value = model.slim_optimize(error_value=float("nan"))
        return {"objective": objective, "status": model.solver.status,
                "value": None if math.isnan(value) else value, "relaxed": relaxed}


result = {}
for label, path in (("baseline", baseline_path), ("candidate", candidate_path)):
    for variant in ("as_built", "r1726_open_in_memory"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile).model
        configure(model)
        if variant == "r1726_open_in_memory":
            model.reactions.R1726.bounds = (0.0, 1000.0)
        s = matrix(model)
        row = {"path": path, "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
               "equations": {rid: model.reactions.get_by_id(rid).build_reaction_string() for rid in ("R2135", "R2145")},
               "wt": solve(model, s, set())}
        fluxes = pfba(model).fluxes
        row["pfba"] = {rid: float(fluxes[rid]) for rid in ROUTE}
        frame = flux_variability_analysis(model, reaction_list=FVA, fraction_of_optimum=0.999, processes=1)
        row["fva_0.999"] = {rid: [float(frame.loc[rid, "minimum"]), float(frame.loc[rid, "maximum"])] for rid in FVA}
        row["no_uptake"] = [no_uptake(model, "R_NTP1"), no_uptake(model, "biomass_C")]
        row["knockouts"] = {}
        for gene in GENES:
            with model:
                model.genes.get_by_id(gene).knock_out()
                knocked = {r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0)}
                row["knockouts"][gene] = solve(model, s, knocked)
        result[f"{label}/{variant}"] = row
        print(label, variant, row["sha256"][:8], row["wt"]["status"], row["wt"]["growth"], row["wt"]["passed"],
              [(r["objective"], r["status"], r["value"]) for r in row["no_uptake"]],
              {g: (None if p["growth"] is None else round(p["growth"], 7)) for g, p in row["knockouts"].items()},
              file=sys.stderr, flush=True)
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
