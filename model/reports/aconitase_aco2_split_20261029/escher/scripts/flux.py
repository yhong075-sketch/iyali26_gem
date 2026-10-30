"""Fluxes for the aconitase before/after Escher figure (read-only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts \
    python ../model/reports/aconitase_aco2_split_20261029/escher/scripts/flux.py MEDIUM.csv PROFILE.json OUT.json \
    BASELINE.xml CANDIDATE.xml
Task A medium + PO1f, Gurobi at FeasibilityTol/OptimalityTol 1e-9 (checked_solve.py), real HOME. The model paths may
be byte-identical short copies (cobra on macOS rejects str SBML paths of 260+ characters).
Four scenarios: WT and the knockout of YALI1D11984g, each on the baseline (059316ef) and the candidate.
Per scenario: checked growth (stops unless WT reproduces the recorded value) and pFBA.
Mapped reactions: every reaction of citrate, cis-aconitate and isocitrate in the cytosol and the mitochondrion
(selected by metabolite membership, biomass assembly newBiom and the pools x* excluded), plus the EXTRA reactions whose
knockout flux changes most between the two models. Copied from the R673 figure; the knockout gene, the reaction set and
EXTRA changed.
"""
import json
import sys
import warnings

from cobra.flux_analysis import pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

EXPECTED_WT = 0.3116608216336117
KO_GENE = "YALI1D11984g"
# citrate, cis-aconitate, isocitrate: mitochondrion then cytosol (ids from before.json)
METS = ["m433[C_mi]", "m409[C_mi]", "m410[C_mi]", "m435[C_cy]", "m436[C_cy]", "m738[C_cy]"]
EXTRA = 12
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
                                        "pfba_status": sol.status}
        print(key, checked["growth"], flush=True)
model = models["candidate"]
FIXED = sorted({r.id for mid in METS for r in model.metabolites.get_by_id(mid).reactions
                if r.id != "newBiom" and not r.id.startswith("x")})
delta = sorted((r for r in fluxes["ko_baseline"] if r not in FIXED),
               key=lambda r: -abs(fluxes["ko_baseline"][r] - fluxes["ko_candidate"][r]))[:EXTRA]
mapped = FIXED + delta
result["mapped"] = {"fixed": FIXED, "largest_ko_change": delta, "metabolites": METS}
result["fluxes"] = {key: {r: f[r] for r in mapped} for key, f in fluxes.items()}
result["reactions"] = {}
result["metabolites"] = {}
for r in (model.reactions.get_by_id(rid) for rid in mapped):
    result["reactions"][r.id] = {"name": r.name, "lb": r.lower_bound, "ub": r.upper_bound,
                                 "gpr": r.gene_reaction_rule,
                                 "gpr_baseline": models["baseline"].reactions.get_by_id(r.id).gene_reaction_rule,
                                 "genes": sorted(g.id for g in r.genes), "subsystem": r.subsystem,
                                 "equation": r.build_reaction_string(use_metabolite_names=True),
                                 "mets": {m.id: c for m, c in r.metabolites.items()}}
    for m in r.metabolites:
        result["metabolites"][m.id] = {"name": m.name, "comp": m.compartment, "degree": len(m.reactions)}
result["lp_count"] = lp
json.dump(result, open(out, "w"), indent=1)
print("LPs", lp, "fixed", len(FIXED), "extra", len(delta), "->", out)
for rid in mapped:
    print(f"{rid:12} wt {fluxes['wt_candidate'][rid]:10.4g} ko_before {fluxes['ko_baseline'][rid]:10.4g} "
          f"ko_after {fluxes['ko_candidate'][rid]:10.4g}  {result['reactions'][rid]['equation'][:110]}")
