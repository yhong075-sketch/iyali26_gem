"""Fluxes for the R1172 closure Escher figure (read-only): checked WT, pFBA, and FVA at the optimum on the supply set.

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    flux.py MEDIUM.csv PROFILE.json OUT.json BEFORE.xml AFTER.xml
Task A medium + PO1f, Gurobi at FeasibilityTol/OptimalityTol 1e-9 (checked_solve.py), real HOME (no override).
Stops unless both WTs reproduce the recorded values (2f352037 = A2; the R1172-closed candidate = Stage 1 S1).
Both paths may be byte-identical short copies: cobra on macOS rejects str SBML paths of 260+ characters.
Copied from the R2121 closure figure: only the reaction sections, the expected WTs and the FVA set change.
"""
import json
import sys
import warnings

from cobra.flux_analysis import flux_variability_analysis, pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

EXPECTED_WT = {"before": 0.3120526953810365, "after": 0.3111005359259494}
SECTIONS = {
    "mitochondrion": ["R83", "R412", "R1974", "R1933", "R1973", "R1972", "R1968", "R834", "R1172"],
    "cytosol_mevalonate": ["R82", "R411", "R410", "R553", "R1996"],
    "acetyl_coa_supply": ["R1936", "R2120", "R2119", "R1894", "R90", "R191", "R192", "R1032", "R88", "R561"],
}
FVA_SET = ["R1172", "R83", "R412", "R1974", "R1933", "R834", "R82", "R411", "R410", "R2119", "R1894", "R90", "R191"]
medium, profile, out, *paths = sys.argv[1:]
mapped = [rid for ids in SECTIONS.values() for rid in ids]
result = {"condition": "Task A medium + PO1f, Gurobi 1e-9 (checked_solve)", "medium": medium, "strain_profile": profile,
          "sections": SECTIONS, "fva_set": FVA_SET, "models": {}}
lp_count = 0
for label, path in zip(("before", "after"), paths):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile)
    model = ctx.model
    configure(model)
    tag = ctx.canonical_model_sha256[:8]
    wt = solve(model, matrix(model), ())
    lp_count += 1
    if not wt["passed"] or abs(wt["growth"] - EXPECTED_WT[label]) > 1e-9:
        sys.exit(f"STOP: {label} {tag} WT {wt} differs from the recorded {EXPECTED_WT[label]}")
    missing = [rid for rid in mapped if rid not in model.reactions]
    if missing:
        sys.exit(f"STOP: {label} {tag} lacks mapped reactions {missing}")
    model.solver.problem.reset()
    solution = pfba(model)
    model.solver.problem.reset()
    fva = flux_variability_analysis(model, reaction_list=FVA_SET, fraction_of_optimum=0.999, processes=1)
    lp_count += 2 + 1 + 2 * len(FVA_SET)
    result["models"][label] = {
        "tag": tag, "file": path.split("/")[-1], "model_sha256": ctx.canonical_model_sha256,
        "medium_sha256": ctx.medium_sha256, "strain_profile_sha256": ctx.strain_profile_sha256, "wt": wt,
        "bounds": {rid: list(model.reactions.get_by_id(rid).bounds) for rid in mapped},
        "pfba": {"status": solution.status, "growth": float(solution.fluxes["biomass_C"]),
                 "total_abs_flux": float(sum(abs(solution.fluxes))),
                 "fluxes": {rid: float(solution.fluxes[rid]) for rid in mapped}},
        "fva_999": {rid: [float(fva.loc[rid, "minimum"]), float(fva.loc[rid, "maximum"])] for rid in FVA_SET}}
    if label == "after":
        result["reactions"] = {r.id: {"name": r.name, "lb": r.lower_bound, "ub": r.upper_bound,
                                      "gpr": r.gene_reaction_rule, "genes": sorted(g.id for g in r.genes),
                                      "subsystem": r.subsystem, "mets": {m.id: c for m, c in r.metabolites.items()}}
                               for r in (model.reactions.get_by_id(rid) for rid in mapped)}
        result["metabolites"] = {m.id: {"name": m.name, "comp": m.compartment, "degree": len(m.reactions),
                                        "reactions": sorted(r.id for r in m.reactions)}
                                 for rid in mapped for m in model.reactions.get_by_id(rid).metabolites}
    print(label, tag, "WT", wt["growth"], "pFBA", solution.status, flush=True)
result["lp_count"] = lp_count
json.dump(result, open(out, "w"), indent=1)
print("LPs", lp_count, "->", out)
