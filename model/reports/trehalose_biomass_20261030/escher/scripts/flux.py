"""Fluxes for the trehalose biomass before/after Escher figure (read-only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts \
    python ../model/reports/trehalose_biomass_20261030/escher/scripts/flux.py MEDIUM.csv PROFILE.json OUT.json \
    BASELINE.xml CANDIDATE.xml
Task A medium + PO1f, Gurobi at FeasibilityTol/OptimalityTol 1e-9 (checked_solve.py), real HOME. The model paths may
be byte-identical short copies (cobra on macOS rejects str SBML paths of 260+ characters).
Four scenarios: WT and the knockout of YALI1E17616g (TPS1), each on the baseline (e1e7c868) and the candidate.
Per scenario: checked growth (stops unless WT reproduces the sweep value) and pFBA.
Mapped reactions: every reaction of trehalose m294[C_cy] and trehalose 6-phosphate m297[C_cy] and the trehalose
exchange R1013 (biomass pseudoreactions excluded), every other reaction of G6P m296[C_cy] and UDP-glucose m8[C_cy],
and the trehalose drain of biomass_C as its own flux (coefficient times biomass_C flux). Adapted from the aconitase
figure's flux.py; the knockout gene, expected WT values and reaction set changed.
"""
import json
import sys
import warnings

from cobra.flux_analysis import pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

EXPECTED_WT = {"baseline": 0.3116608216336117, "candidate": 0.31229636510178727}  # ../ko_compare_summary.json
KO_GENE = "YALI1E17616g"
BRANCH_METS = ["m294[C_cy]", "m297[C_cy]"]
SUPPLY_METS = ["m296[C_cy]", "m8[C_cy]"]
TERM = "m294[C_cy]"
DRAIN = "biomass_C_trehalose_term"
PSEUDO = {"biomass_C", "xBIOMASS", "newBiom", "R1372", "R1387", "R1710"}
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
    coef = model.reactions.biomass_C.metabolites.get(model.metabolites.get_by_id(TERM), 0.0)
    result[f"{label}_trehalose_coefficient"] = coef
    for case in ("wt", "ko"):
        with model:
            knocked = ()
            if case == "ko":
                model.genes.get_by_id(KO_GENE).knock_out()
                knocked = tuple(r.id for r in model.genes.get_by_id(KO_GENE).reactions if r.bounds == (0, 0))
            checked = solve(model, matrix(model), knocked)
            if not checked["passed"]:
                sys.exit(f"STOP: {label} {case} {checked}")
            if case == "wt" and abs(checked["growth"] - EXPECTED_WT[label]) > 1e-9:
                sys.exit(f"STOP: {label} WT {checked['growth']} differs from {EXPECTED_WT[label]}")
            sol = pfba(model)
            lp += 3
            key = f"{case}_{label}"
            fluxes[key] = sol.fluxes.to_dict()
            fluxes[key][DRAIN] = -coef * fluxes[key]["biomass_C"]
            result["scenarios"][key] = {"checked": checked, "closed_by_knockout": list(knocked),
                                        "pfba_status": sol.status}
        print(key, checked["growth"], flush=True)
model = models["baseline"]
branch = sorted({r.id for mid in BRANCH_METS for r in model.metabolites.get_by_id(mid).reactions} - PSEUDO | {"R1013"})
supply = sorted({r.id for mid in SUPPLY_METS for r in model.metabolites.get_by_id(mid).reactions} - PSEUDO - set(branch))
mapped = branch + supply
result["mapped"] = {"branch": branch, "supply": supply, "drain": DRAIN}
result["fluxes"] = {key: {r: f[r] for r in mapped + [DRAIN]} for key, f in fluxes.items()}
result["reactions"] = {}
result["metabolites"] = {}
for r in (model.reactions.get_by_id(rid) for rid in mapped):
    result["reactions"][r.id] = {"name": r.name, "lb": r.lower_bound, "ub": r.upper_bound,
                                 "gpr": r.gene_reaction_rule, "gpr_baseline": r.gene_reaction_rule,
                                 "genes": sorted(g.id for g in r.genes), "subsystem": r.subsystem,
                                 "equation": r.build_reaction_string(use_metabolite_names=True),
                                 "mets": {m.id: c for m, c in r.metabolites.items()}}
    for m in r.metabolites:
        result["metabolites"][m.id] = {"name": m.name, "comp": m.compartment, "degree": len(m.reactions)}
result["lp_count"] = lp
json.dump(result, open(out, "w"), indent=1)
print("LPs", lp, "branch", len(branch), "supply", len(supply), "->", out)
for rid in mapped + [DRAIN]:
    print(f"{rid:26} wt_before {fluxes['wt_baseline'][rid]:10.4g} wt_after {fluxes['wt_candidate'][rid]:10.4g} "
          f"ko_before {fluxes['ko_baseline'][rid]:10.4g} ko_after {fluxes['ko_candidate'][rid]:10.4g}")
