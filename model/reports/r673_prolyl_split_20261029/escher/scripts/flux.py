"""Fluxes for the R673 before/after Escher figure (read-only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts \
    python ../model/reports/r673_prolyl_split_20261029/escher/scripts/flux.py MEDIUM.csv PROFILE.json OUT.json \
    BASELINE.xml CANDIDATE.xml
Task A medium + PO1f, Gurobi at FeasibilityTol/OptimalityTol 1e-9 (checked_solve.py), real HOME. The model paths may
be byte-identical short copies (cobra on macOS rejects str SBML paths of 260+ characters).
Four scenarios: WT and the knockout of YALI1E05929g, each on the baseline (ca43ae0c) and the candidate.
Per scenario: checked growth (stops unless WT reproduces the recorded value) and pFBA.
Mapped reactions: R673 and every other reaction of cytosolic L-proline (m765), tRNA(Pro) (m905) and Pro-tRNA(Pro)
(m906), except the biomass assembly newBiom and the amino-acid pool xAMINOACID. Copied from the R171 figure; the
knockout gene and the reaction set changed, and no reactions are added by flux difference.
"""
import json
import sys
import warnings

from cobra.flux_analysis import pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

EXPECTED_WT = 0.3116608216336117
KO_GENE = "YALI1E05929g"
FIXED = ["R673", "TRNA_BIOMASS_R673", "R1387", "R1710", "R1225", "R1372", "R2040", "R509", "R689", "R893"]
EXTRA = 0
medium, profile, out, base_path, cand_path = sys.argv[1:6]
result = {"condition": "Task A medium + PO1f, Gurobi 1e-9 (checked_solve)", "ko_gene": KO_GENE, "scenarios": {}}
fluxes, models, lp = {}, {}, 0
for label, path in (("baseline", base_path), ("candidate", cand_path)):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile)
    model = ctx.model
    configure(model)
    models[label] = model
    result[f"{label}_sha256"] = ctx.canonical_model_sha256
    result["medium_sha256"], result["strain_profile_sha256"] = ctx.medium_sha256, ctx.strain_profile_sha256
    for case in ("wt", "ko"):
        with model:
            knocked = ()
            if case == "ko":
                model.genes.get_by_id(KO_GENE).knock_out()
                knocked = tuple(r.id for r in model.genes.get_by_id(KO_GENE).reactions if r.bounds == (0, 0))
            checked = solve(model, matrix(model), knocked)
            if not checked["passed"]:
                sys.exit(f"STOP: {label} {case} {checked}")
            if case == "wt" and abs(checked["growth"] - EXPECTED_WT) > 1e-9:
                sys.exit(f"STOP: {label} WT {checked['growth']} differs from {EXPECTED_WT}")
            sol = pfba(model)
            lp += 3
            key = f"{case}_{label}"
            fluxes[key] = sol.fluxes.to_dict()
            result["scenarios"][key] = {"checked": checked, "closed_by_knockout": list(knocked),
                                        "pfba_status": sol.status, "R673_bounds": list(model.reactions.R673.bounds)}
        print(key, checked["growth"], flush=True)
delta = sorted((r for r in fluxes["ko_baseline"] if r not in FIXED),
               key=lambda r: -abs(fluxes["ko_baseline"][r] - fluxes["ko_candidate"][r]))[:EXTRA]
mapped = FIXED + delta
model = models["candidate"]
result["mapped"] = {"fixed": FIXED, "largest_ko_change": delta}
result["fluxes"] = {key: {r: f[r] for r in mapped} for key, f in fluxes.items()}
result["reactions"] = {}
result["metabolites"] = {}
for r in (model.reactions.get_by_id(rid) for rid in mapped):
    result["reactions"][r.id] = {"name": r.name, "lb": r.lower_bound, "ub": r.upper_bound,
                                 "gpr": r.gene_reaction_rule,
                                 "gpr_baseline": models["baseline"].reactions.get_by_id(r.id).gene_reaction_rule,
                                 "genes": sorted(g.id for g in r.genes), "subsystem": r.subsystem,
                                 "mets": {m.id: c for m, c in r.metabolites.items()}}
    for m in r.metabolites:
        result["metabolites"][m.id] = {"name": m.name, "comp": m.compartment, "degree": len(m.reactions)}
result["lp_count"] = lp
json.dump(result, open(out, "w"), indent=1)
print("LPs", lp, "mapped", len(mapped), "->", out)
