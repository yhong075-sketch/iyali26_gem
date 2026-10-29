"""Fluxes for the UM11 nuclear-RNR Escher figure (read-only): one model as built and under the in-memory Option A
what-if (the two thioredoxin nucleus<->cytosol links), which was NOT implemented (the user chose Option B).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts:../model/reports/um11_nuclear_rnr_20261029/scripts \
    python ../model/reports/um11_nuclear_rnr_20261029/escher/scripts/flux.py MEDIUM.csv PROFILE.json OUT.json MODEL.xml
Task A medium + PO1f, Gurobi at FeasibilityTol/OptimalityTol 1e-9 (checked_solve.py), real HOME (no override).
Stops unless WT reproduces the recorded value, the uracil exchange is open to uptake 1000 and R612 is (0, 0).
MODEL may be a byte-identical short copy: cobra on macOS rejects str SBML paths of 260+ characters.
Per variant: checked WT, pFBA, FVA at 0 % of the mapped reactions (capacity, Task A), whole-model FVA at 0 % on
Task A (blocked set, for node colours), and the no-uptake FVA of the mapped reactions (every exchange (0, 1000),
positive lower bounds incl. xMAINTENANCE relaxed, biomass fixed at 0): non-zero there = an internal loop.
Adapted from the UM2 figure's flux.py: the variants come from ../../scripts/um11_apply.py instead of bound blocks.
"""
import json
import sys
import warnings

from cobra.flux_analysis import flux_variability_analysis, pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from um11_apply import apply

EXPECTED_WT = 0.3116608216336117
SECTIONS = {
    "nucleus": ["R707", "R709", "R711", "R1914", "R782", "R1033", "R1076", "R1077", "R1093", "R1136", "R1353",
                "R1104", "R789", "T_TRX1tn", "T_TRXSStn"],
    "cytosol": ["R706", "R708", "R710", "R1913", "R781", "R755", "R702", "R703", "R704", "R705"],
}
VARIANTS = {"as_built": "as_built", "whatif_A": "trx+trxss"}
medium, profile, out, path = sys.argv[1:5]
mapped = [rid for ids in SECTIONS.values() for rid in ids]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile)
model = ctx.model
configure(model)
uracil = {r.id: [r.name, list(r.bounds)] for r in model.exchanges if "uracil" in r.name.lower()}
if not any(b[0] == -1000 for _, b in uracil.values()) or list(model.reactions.R612.bounds) != [0.0, 0.0]:
    sys.exit(f"STOP: uracil {uracil} or R612 {model.reactions.R612.bounds} not as expected")
pair = lambda df, rid: [float(df.loc[rid, "minimum"]), float(df.loc[rid, "maximum"])]
result = {"condition": "Task A medium + PO1f, Gurobi 1e-9 (checked_solve)", "medium": medium, "strain_profile": profile,
          "model_sha256": ctx.canonical_model_sha256, "medium_sha256": ctx.medium_sha256,
          "strain_profile_sha256": ctx.strain_profile_sha256, "uracil_exchanges": uracil,
          "R612_bounds": list(model.reactions.R612.bounds), "sections": SECTIONS, "variants_apply": VARIANTS,
          "variants": {}}
lp = 0
reactions, metabolites = {}, {}
for label, variant in VARIANTS.items():
    with model:
        apply(model, variant)
        here = [rid for rid in mapped if rid in model.reactions]
        wt = solve(model, matrix(model), ())
        if not wt["passed"] or abs(wt["growth"] - EXPECTED_WT) > 1e-9:
            sys.exit(f"STOP: {label} WT {wt} differs from the recorded {EXPECTED_WT}")
        sol = pfba(model)
        cap = flux_variability_analysis(model, reaction_list=here, fraction_of_optimum=0.0, processes=1)
        whole = flux_variability_analysis(model, fraction_of_optimum=0.0, processes=1)
        lp += 2 + 2 * len(here) + 2 * len(model.reactions)
        with model:
            for r in model.exchanges:
                r.bounds = (0.0, 1000.0)
            for r in model.reactions:
                if r.lower_bound > 0:
                    r.lower_bound = 0.0
            model.reactions.biomass_C.bounds = (0, 0)
            nou = flux_variability_analysis(model, reaction_list=here, fraction_of_optimum=0.0, processes=1)
            lp += 2 * len(here)
        result["variants"][label] = {
            "wt": wt, "present": here, "bounds": {rid: list(model.reactions.get_by_id(rid).bounds) for rid in here},
            "pfba": {"status": sol.status, "fluxes": {rid: float(sol.fluxes[rid]) for rid in here}},
            "capacity_taskA": {rid: pair(cap, rid) for rid in here},
            "no_uptake": {rid: pair(nou, rid) for rid in here},
            "blocked_taskA": sorted(r for r in whole.index if max(abs(whole.loc[r, "minimum"]), abs(whole.loc[r, "maximum"])) <= 1e-9)}
        for r in (model.reactions.get_by_id(rid) for rid in here):
            reactions[r.id] = {"name": r.name, "lb": r.lower_bound, "ub": r.upper_bound, "gpr": r.gene_reaction_rule,
                               "genes": sorted(g.id for g in r.genes), "subsystem": r.subsystem,
                               "mets": {m.id: c for m, c in r.metabolites.items()}}
            for m in r.metabolites:
                metabolites.setdefault(m.id, {"name": m.name, "comp": m.compartment, "degree": {}, "reactions": {}})
                metabolites[m.id]["degree"][label] = len(m.reactions)
                metabolites[m.id]["reactions"][label] = sorted(x.id for x in m.reactions)
    print(label, "WT", wt["growth"], "blocked Task A (whole model)", len(result["variants"][label]["blocked_taskA"]), flush=True)
result["reactions"] = reactions
result["metabolites"] = metabolites
result["lp_count"] = lp
json.dump(result, open(out, "w"), indent=1)
print("LPs", lp, "->", out)
