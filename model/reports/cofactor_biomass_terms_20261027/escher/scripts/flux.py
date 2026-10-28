"""Fluxes for the cofactor-terms Escher figure (read-only): checked WT, pFBA and FVA at the optimum, both models.

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/cofactor_biomass_terms_20261027/escher/scripts/flux.py MEDIUM.csv PROFILE.json \
    CENSUS_OPEN.json CENSUS_TASKA.json OUT.json BEFORE.xml AFTER.xml
Task A medium + PO1f, Gurobi at FeasibilityTol/OptimalityTol 1e-9 (checked_solve.py). The mapped set is the
union of the reactions freed in the open and the Task A census comparisons, plus biomass_C. Also exports
their stoichiometry, names, compartments and GPRs from the AFTER model for the map builder.
"""
import json
import sys
import warnings

from cobra.flux_analysis import flux_variability_analysis, pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

EXPECTED_WT = {"a1237821": 0.31242085277746734, "00e8e397": 0.31200492019479353}
medium, profile, census_open, census_taska, out, *paths = sys.argv[1:]
freed_open = [r["id"] for r in json.load(open(census_open))["freed_reactions"]]
freed_taska = [r["id"] for r in json.load(open(census_taska))["freed_reactions"]]
mapped = sorted(set(freed_open) | set(freed_taska), key=lambda rid: int(rid[1:])) + ["biomass_C"]
result = {"condition": "Task A medium + PO1f, Gurobi 1e-9 (checked_solve)", "medium": medium, "strain_profile": profile,
          "freed_open": freed_open, "freed_taskA": freed_taska, "mapped": mapped, "models": {}}
lp_count = 0

for path in paths:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile)
    model = ctx.model
    configure(model)
    tag = ctx.canonical_model_sha256[:8]
    wt = solve(model, matrix(model), ())
    if abs(wt["growth"] - EXPECTED_WT[tag]) > 1e-9:
        sys.exit(f"STOP: {tag} WT {wt['growth']} differs from the recorded {EXPECTED_WT[tag]}")
    solution = pfba(model)
    fva = flux_variability_analysis(model, reaction_list=mapped, fraction_of_optimum=1.0, processes=1)
    lp_count += 2 + 2 * len(mapped)
    uracil = {r.id: [r.name, list(r.bounds)] for r in model.exchanges if "uracil" in r.name.lower()}
    result["models"][tag] = {
        "file": path.split("/")[-1], "model_sha256": ctx.canonical_model_sha256, "medium_sha256": ctx.medium_sha256,
        "strain_profile_sha256": ctx.strain_profile_sha256, "wt": wt, "uracil_exchanges": uracil,
        "R612_bounds": list(model.reactions.R612.bounds),
        "pfba": {"status": solution.status, "growth": float(solution.fluxes["biomass_C"]),
                 "fluxes": {rid: float(solution.fluxes[rid]) for rid in mapped}},
        "fva_at_optimum": {rid: [float(fva.loc[rid, "minimum"]), float(fva.loc[rid, "maximum"])] for rid in mapped}}
    if tag == "00e8e397":
        result["reactions"] = {r.id: {"name": r.name, "lb": r.lower_bound, "ub": r.upper_bound, "gpr": r.gene_reaction_rule,
                                      "genes": sorted(g.id for g in r.genes), "subsystem": r.subsystem,
                                      "mets": {m.id: c for m, c in r.metabolites.items()}}
                               for r in (model.reactions.get_by_id(rid) for rid in mapped)}
        result["metabolites"] = {m.id: {"name": m.name, "comp": m.compartment, "degree": len(m.reactions)}
                                 for rid in mapped for m in model.reactions.get_by_id(rid).metabolites}
    print(tag, "WT", wt["growth"], "pFBA", solution.status, "uracil", uracil, "R612", model.reactions.R612.bounds)

result["lp_count"] = lp_count
json.dump(result, open(out, "w"), indent=1)
print("LPs", lp_count, "->", out)
