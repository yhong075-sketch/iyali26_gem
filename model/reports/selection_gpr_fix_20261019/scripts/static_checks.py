"""Static checks on the new candidate (no optimization).

Usage: python static_checks.py CANDIDATE.xml CANDIDATE.build.json OUT.json
Balance of the touched reactions, the m325 pool, gene-to-reaction knockout propagation and the
selection/restore records of the build.
"""
import json
import sys

from cobra.io import read_sbml_model

model_path, build_path, out_path = sys.argv[1:4]
model = read_sbml_model(model_path)
build = json.load(open(build_path))
out = {"balance": {}, "m325_reactions": {}, "knockout_closes": {}, "records": {}}
for rid in ("R671", "R159", "R190", "R765", "R766"):
    r = model.reactions.get_by_id(rid)
    out["balance"][rid] = {"gpr": r.gene_reaction_rule, "residual": r.check_mass_balance()}
for r in model.metabolites.get_by_id("m325[C_cy]").reactions:
    out["m325_reactions"][r.id] = r.metabolites[model.metabolites.get_by_id("m325[C_cy]")]
for gene in ("YALI1F23441g", "YALI1E11370g", "YALI1C33005g", "YALI1D09420g", "YALI1E11768g",
             "YALI1E07744g", "YALI1D02625g"):
    with model:
        model.genes.get_by_id(gene).knock_out()
        out["knockout_closes"][gene] = sorted(r.id for r in model.reactions if r.bounds == (0, 0)
                                              and r.id in {x.id for x in model.genes.get_by_id(gene).reactions})
    out["knockout_closes"][gene] = {"closes": out["knockout_closes"][gene],
                                    "reactions": sorted(r.id for r in model.genes.get_by_id(gene).reactions)}
selection = build["reaction_selection"]
out["records"]["selection_complete"] = selection["complete"]
out["records"]["selection_rows"] = {row["item"]: {"status": row["status"], "fields": row["fields"]}
                                    for row in selection["items"] if row["item"] in ("R671", "R159", "R190", "R765", "R766")}
for key in ("selection_gpr_restore", "selection_gpr_restore_r671", "selection_cpa_ura2_restore"):
    out["records"][key] = selection.get(key)
json.dump(out, open(out_path, "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
