"""Fluxes for the UM7 nuclear-diffusion Escher figure (read-only): checked WT, pFBA, FVA at the optimum and FVA at 0 %
(capacity, open and Task A), both models.

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/um7_nuclear_transports_20261028/escher/scripts/flux.py MEDIUM.csv PROFILE.json OUT.json BEFORE.xml AFTER.xml
Task A medium + PO1f, Gurobi at FeasibilityTol/OptimalityTol 1e-9 (checked_solve.py), real HOME (no override).
Stops unless WT reproduces the recorded values, the uracil exchange is open to uptake 1000 and R612 is (0, 0).
Both paths may be byte-identical short copies: cobra on macOS rejects str SBML paths of 260+ characters.
Adapted from the dolichol-route figure: the reaction set changed and the 0 % FVA (Task A and every exchange open)
was added, because pFBA carries no flux through these reactions at the optimum in either model.
"""
import json
import sys
import warnings

import cobra
from cobra.flux_analysis import flux_variability_analysis, pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

EXPECTED_WT = {"9d655c2b": 0.311660821633612}
SECTIONS = {
    "nucleus": ["R1082", "R569", "R1270", "R579", "R567", "R1267", "R1017", "R1212", "R1214", "R_PPItn", "R_ADPRtn"],
    "cytosol": ["R580", "R568", "R578", "R1912", "R419", "R576"],
    "homocitrate": ["R965", "R1001", "R396", "R1055", "R1000", "R92"],
}
medium, profile, out, *paths = sys.argv[1:]
mapped = [rid for ids in SECTIONS.values() for rid in ids]
result = {"condition": "Task A medium + PO1f, Gurobi 1e-9 (checked_solve); capacity also with every exchange open",
          "medium": medium, "strain_profile": profile, "sections": SECTIONS, "models": {}}
lp_count = 0
for label, path in zip(("before", "after"), paths):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile)
        open_model = cobra.io.read_sbml_model(path)
    model = ctx.model
    configure(model)
    configure(open_model)
    for r in open_model.exchanges:
        r.bounds = (-1000, 1000)
    tag = ctx.canonical_model_sha256[:8]
    uracil = {r.id: [r.name, list(r.bounds)] for r in model.exchanges if "uracil" in r.name.lower()}
    if not any(b[0] == -1000 for _, b in uracil.values()) or list(model.reactions.R612.bounds) != [0.0, 0.0]:
        sys.exit(f"STOP: {tag} uracil {uracil} or R612 {model.reactions.R612.bounds} not as expected")
    wt = solve(model, matrix(model), ())
    expected = EXPECTED_WT.get(tag, EXPECTED_WT["9d655c2b"])
    if not wt["passed"] or abs(wt["growth"] - expected) > 1e-9:
        sys.exit(f"STOP: {tag} WT {wt} differs from the recorded {expected}")
    present = [rid for rid in mapped if rid in model.reactions]
    solution = pfba(model)
    fva1 = flux_variability_analysis(model, reaction_list=present, fraction_of_optimum=1.0, processes=1)
    fva0 = flux_variability_analysis(model, reaction_list=present, fraction_of_optimum=0.0, processes=1)
    fvao = flux_variability_analysis(open_model, reaction_list=present, fraction_of_optimum=0.0, processes=1)
    lp_count += 1 + 2 + 1 + 6 * len(present)
    pair = lambda df, rid: [float(df.loc[rid, "minimum"]), float(df.loc[rid, "maximum"])]
    result["models"][label] = {
        "tag": tag, "file": path.split("/")[-1], "model_sha256": ctx.canonical_model_sha256,
        "medium_sha256": ctx.medium_sha256, "strain_profile_sha256": ctx.strain_profile_sha256, "wt": wt,
        "uracil_exchanges": uracil, "R612_bounds": list(model.reactions.R612.bounds),
        "absent": [rid for rid in mapped if rid not in model.reactions],
        "bounds": {rid: list(model.reactions.get_by_id(rid).bounds) for rid in present},
        "pfba": {"status": solution.status, "growth": float(solution.fluxes["biomass_C"]),
                 "fluxes": {rid: float(solution.fluxes[rid]) for rid in present}},
        "fva_at_optimum": {rid: pair(fva1, rid) for rid in present},
        "fva_capacity_taskA": {rid: pair(fva0, rid) for rid in present},
        "fva_capacity_open": {rid: pair(fvao, rid) for rid in present}}
    if label == "after":
        result["reactions"] = {r.id: {"name": r.name, "lb": r.lower_bound, "ub": r.upper_bound, "gpr": r.gene_reaction_rule,
                                      "genes": sorted(g.id for g in r.genes), "subsystem": r.subsystem,
                                      "mets": {m.id: c for m, c in r.metabolites.items()}}
                               for r in (model.reactions.get_by_id(rid) for rid in mapped)}
        result["metabolites"] = {m.id: {"name": m.name, "comp": m.compartment, "degree": len(m.reactions),
                                        "reactions": sorted(r.id for r in m.reactions)}
                                 for rid in mapped for m in model.reactions.get_by_id(rid).metabolites}
    print(label, tag, "WT", wt["growth"], "pFBA", solution.status, "uracil", uracil, "R612", model.reactions.R612.bounds, flush=True)
result["lp_count"] = lp_count
json.dump(result, open(out, "w"), indent=1)
print("LPs", lp_count, "->", out)
