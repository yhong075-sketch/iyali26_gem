"""Trace how the new model makes heme without R792 (HEM4, YALI1F13039g knocked out). Read-only, 2 LPs (pFBA).

Usage (from platform/): PYTHONPATH=$PWD python ../model/reports/heme_iron_biomass_20261022/scripts/hem4_bypass.py \
    NEW.xml MEDIUM.csv PROFILE.json OUT.json
"""
import json
import sys
import warnings

from cobra.flux_analysis import pfba

from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

new_path, medium, profile, out = sys.argv[1:5]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=new_path, media_path=medium, strain_profile_path=profile)
model = ctx.model
model.solver = "glpk"
r792 = model.reactions.R792
products = [m for m, c in r792.metabolites.items() if c > 0 and not m.id.startswith("m10[")]
substrates = [m for m, c in r792.metabolites.items() if c < 0]
model.genes.YALI1F13039g.knock_out()
solution = pfba(model)


def touching(metabolite):
    rows = []
    for r in sorted(metabolite.reactions, key=lambda r: r.id):
        flux = float(solution.fluxes[r.id])
        rows.append({"id": r.id, "name": r.name, "gpr": r.gene_reaction_rule, "bounds": list(r.bounds),
                     "flux": flux, "equation": r.build_reaction_string(use_metabolite_names=True)})
    return rows


result = {"model_sha256": ctx.canonical_model_sha256, "condition": "Task A medium + PO1f, GLPK, YALI1F13039g knocked out",
          "status": solution.status, "growth": float(solution.fluxes["biomass_C"]),
          "R792": {"equation": r792.build_reaction_string(use_metabolite_names=True), "bounds": list(r792.bounds),
                   "flux": float(solution.fluxes["R792"])},
          "substrates": {m.id: {"name": m.name, "reactions": touching(m)} for m in substrates},
          "products": {m.id: {"name": m.name, "reactions": touching(m)} for m in products},
          "lp_count": 2}
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
print(result["status"], result["growth"], result["R792"])
for side in ("substrates", "products"):
    for mid, row in result[side].items():
        print(side, mid, row["name"])
        for r in row["reactions"]:
            print(f"   {r['id']} flux={r['flux']:.3g} {r['bounds']} gpr={r['gpr']!r} {r['name']}: {r['equation']}")
