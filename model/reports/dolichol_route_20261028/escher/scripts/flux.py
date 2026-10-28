"""Fluxes for the dolichol-route Escher figure (read-only): checked WT, pFBA and FVA at the optimum, both models.

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/dolichol_route_20261028/escher/scripts/flux.py MEDIUM.csv PROFILE.json OUT.json BEFORE.xml AFTER.xml
Task A medium + PO1f, Gurobi at FeasibilityTol/OptimalityTol 1e-9 (checked_solve.py), real HOME (no override).
Stops unless WT reproduces the recorded values, the uracil exchange is open to uptake 1000 and R612 is (0, 0).
AFTER may be the byte-identical research copy: cobra on macOS rejects absolute SBML paths over 260 characters.
"""
import json
import sys
import warnings

from cobra.flux_analysis import flux_variability_analysis, pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

EXPECTED_WT = {"ce9e1112": 0.3116823733427652, "9d655c2b": 0.311660821633612}
SECTIONS = {
    "route": ["R282", "R322", "R2179", "R2180", "R2178", "R_DOLSYN", "R287"],
    "pools": ["R1101", "R2177", "R2182", "R2181", "R288", "R2102", "R289", "R372", "R374", "R375", "R2106", "R2138", "biomass_C"],
    "um6": [f"R{i}" for i in range(217, 236)],
}
medium, profile, out, *paths = sys.argv[1:]
mapped = [rid for ids in SECTIONS.values() for rid in ids]
result = {"condition": "Task A medium + PO1f, Gurobi 1e-9 (checked_solve)", "medium": medium, "strain_profile": profile,
          "sections": SECTIONS, "models": {}}
lp_count = 0
for path in paths:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile)
    model = ctx.model
    configure(model)
    tag = ctx.canonical_model_sha256[:8]
    uracil = {r.id: [r.name, list(r.bounds)] for r in model.exchanges if "uracil" in r.name.lower()}
    if not any(b[0] == -1000 for _, b in uracil.values()) or list(model.reactions.R612.bounds) != [0.0, 0.0]:
        sys.exit(f"STOP: {tag} uracil {uracil} or R612 {model.reactions.R612.bounds} not as expected")
    wt = solve(model, matrix(model), ())
    if not wt["passed"] or abs(wt["growth"] - EXPECTED_WT[tag]) > 1e-9:
        sys.exit(f"STOP: {tag} WT {wt} differs from the recorded {EXPECTED_WT[tag]}")
    present = [rid for rid in mapped if rid in model.reactions]
    solution = pfba(model)
    fva = flux_variability_analysis(model, reaction_list=present, fraction_of_optimum=1.0, processes=1)
    lp_count += 1 + 2 + 1 + 2 * len(present)
    result["models"][tag] = {
        "file": path.split("/")[-1], "model_sha256": ctx.canonical_model_sha256, "medium_sha256": ctx.medium_sha256,
        "strain_profile_sha256": ctx.strain_profile_sha256, "wt": wt, "uracil_exchanges": uracil,
        "R612_bounds": list(model.reactions.R612.bounds), "absent": [rid for rid in mapped if rid not in model.reactions],
        "bounds": {rid: list(model.reactions.get_by_id(rid).bounds) for rid in present},
        "pfba": {"status": solution.status, "growth": float(solution.fluxes["biomass_C"]),
                 "fluxes": {rid: float(solution.fluxes[rid]) for rid in present}},
        "fva_at_optimum": {rid: [float(fva.loc[rid, "minimum"]), float(fva.loc[rid, "maximum"])] for rid in present}}
    if tag == "9d655c2b":
        result["reactions"] = {r.id: {"name": r.name, "lb": r.lower_bound, "ub": r.upper_bound, "gpr": r.gene_reaction_rule,
                                      "genes": sorted(g.id for g in r.genes), "subsystem": r.subsystem,
                                      "mets": {m.id: c for m, c in r.metabolites.items()}}
                               for r in (model.reactions.get_by_id(rid) for rid in mapped)}
        result["metabolites"] = {m.id: {"name": m.name, "comp": m.compartment, "degree": len(m.reactions)}
                                 for rid in mapped for m in model.reactions.get_by_id(rid).metabolites}
    print(tag, "WT", wt["growth"], "pFBA", solution.status, "uracil", uracil, "R612", model.reactions.R612.bounds, flush=True)
result["lp_count"] = lp_count
json.dump(result, open(out, "w"), indent=1)
print("LPs", lp_count, "->", out)
