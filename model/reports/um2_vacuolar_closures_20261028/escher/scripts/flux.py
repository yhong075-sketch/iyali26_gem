"""Fluxes for the UM2 vacuolar-closure Escher figure (read-only): one model as built and under two in-memory
what-if options that were NOT implemented (the user chose Option A0, keep every closure).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/um2_vacuolar_closures_20261028/escher/scripts/flux.py MEDIUM.csv PROFILE.json OUT.json MODEL.xml
Task A medium + PO1f, Gurobi at FeasibilityTol/OptimalityTol 1e-9 (checked_solve.py), real HOME (no override).
Stops unless WT reproduces the recorded value, the uracil exchange is open to uptake 1000 and R612 is (0, 0).
MODEL may be a byte-identical short copy: cobra on macOS rejects str SBML paths of 260+ characters.
Per variant: checked WT, pFBA, FVA at 0 % of the mapped reactions (capacity, Task A), whole-model FVA at 0 % on
Task A (blocked set, for node colours), and the no-uptake FVA of the mapped reactions (every exchange (0, 1000),
positive lower bounds incl. xMAINTENANCE relaxed, biomass fixed at 0): non-zero there = an internal loop.
Adapted from the UM7 figure's flux.py: one model, in-memory variants, whole-model FVA and the no-uptake FVA added.
"""
import json
import sys
import warnings

from cobra.flux_analysis import flux_variability_analysis, pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

EXPECTED_WT = 0.311660821633612
SECTIONS = {
    "pump": ["R795", "R1163", "R1363"],
    "amino": ["R869", "R870", "R878", "R879", "R882", "R883", "R887", "R888", "R897", "R898", "R871", "R876",
              "R2055", "R2056"],
    "polyamine": ["R918", "R925", "R927"],
    "dipeptide": ["R2018", "R2019", "R2020", "R2021", "R2026", "R2027", "R2028", "R2029", "R2030", "R2031", "R2032",
                  "R2033", "R2034", "R2035", "R2036", "R2037", "R2038", "R2039", "R2040"],
}
# the same option blocks as ../../scripts/stage1_um2.py
VARIANTS = {
    "as_built": {},
    "cprime": {"R795": (0, 1000), "R1363": (0, 1000), "R1163": (-1000, 0),
               **{r: (0, 1000) for r in ("R869", "R878", "R887", "R897", "R870", "R879", "R883", "R888", "R898",
                                          "R871", "R876")}},
    "open17_as_written": {"R795": (0, 1000), "R1363": (0, 1000), "R1163": (0, 1000),
                          **{r: (0, 1000) for r in ("R869", "R878", "R887", "R897", "R870", "R879", "R883", "R888",
                                                     "R898", "R871", "R876", "R918", "R925", "R927")}},
}
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
          "R612_bounds": list(model.reactions.R612.bounds), "sections": SECTIONS,
          "variants_bounds": {k: {r: list(b) for r, b in v.items()} for k, v in VARIANTS.items()}, "variants": {}}
lp = 0
for label, bounds in VARIANTS.items():
    with model:
        for rid, b in bounds.items():
            model.reactions.get_by_id(rid).bounds = b
        wt = solve(model, matrix(model), ())
        if not wt["passed"] or abs(wt["growth"] - EXPECTED_WT) > 1e-9:
            sys.exit(f"STOP: {label} WT {wt} differs from the recorded {EXPECTED_WT}")
        sol = pfba(model)
        cap = flux_variability_analysis(model, reaction_list=mapped, fraction_of_optimum=0.0, processes=1)
        whole = flux_variability_analysis(model, fraction_of_optimum=0.0, processes=1)
        lp += 2 + 2 * len(mapped) + 2 * len(model.reactions)
        with model:
            for r in model.exchanges:
                r.bounds = (0.0, 1000.0)
            for r in model.reactions:
                if r.lower_bound > 0:
                    r.lower_bound = 0.0
            model.reactions.biomass_C.bounds = (0, 0)
            nou = flux_variability_analysis(model, reaction_list=mapped, fraction_of_optimum=0.0, processes=1)
            lp += 2 * len(mapped)
        result["variants"][label] = {
            "wt": wt, "bounds": {rid: list(model.reactions.get_by_id(rid).bounds) for rid in mapped},
            "pfba": {"status": sol.status, "fluxes": {rid: float(sol.fluxes[rid]) for rid in mapped}},
            "capacity_taskA": {rid: pair(cap, rid) for rid in mapped},
            "no_uptake": {rid: pair(nou, rid) for rid in mapped},
            "blocked_taskA": sorted(r for r in whole.index if max(abs(whole.loc[r, "minimum"]), abs(whole.loc[r, "maximum"])) <= 1e-9)}
    print(label, "WT", wt["growth"], "blocked Task A (whole model)", len(result["variants"][label]["blocked_taskA"]), flush=True)
result["reactions"] = {r.id: {"name": r.name, "lb": r.lower_bound, "ub": r.upper_bound, "gpr": r.gene_reaction_rule,
                              "genes": sorted(g.id for g in r.genes), "subsystem": r.subsystem,
                              "mets": {m.id: c for m, c in r.metabolites.items()}}
                       for r in (model.reactions.get_by_id(rid) for rid in mapped)}
result["metabolites"] = {m.id: {"name": m.name, "comp": m.compartment, "degree": len(m.reactions),
                                "reactions": sorted(r.id for r in m.reactions)}
                         for rid in mapped for m in model.reactions.get_by_id(rid).metabolites}
result["lp_count"] = lp
json.dump(result, open(out, "w"), indent=1)
print("LPs", lp, "->", out)
