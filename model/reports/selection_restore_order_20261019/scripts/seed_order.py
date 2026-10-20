"""Print the species order of each restored reaction after the dc471d7 restore, in-memory, no optimization."""
import json
import sys

from cobra import Metabolite, Model, Reaction

from scripts.gem_annotate.selection_gpr_restore import CHANGE_SETS, apply_selection_gpr_restore


def add(model, rid, state):
    reaction = Reaction(rid, lower_bound=state["bounds"][0], upper_bound=state["bounds"][1])
    model.add_reactions([reaction])
    for mid in state["stoichiometry"]:
        if mid not in model.metabolites:
            model.add_metabolites([Metabolite(mid, compartment="C_cy")])
    reaction.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in state["stoichiometry"].items()})
    reaction.gene_reaction_rule = state["gpr"]
    return reaction


def fixture(spec):
    model = Model("seed_order")
    model.add_metabolites([Metabolite(mid, formula=r["formula"], charge=r["charge"], compartment=r["compartment"])
                           for mid, r in spec["species"].items()])
    context = spec.get("context", {})
    for rid, expected in context.items():
        if rid not in spec["reactions"]:
            add(model, rid, expected)
    for rid, row in spec["reactions"].items():
        reaction = add(model, rid, row["cleared"])
        reaction.name = context.get(rid, {}).get("name", "")
        if "ec-code" in context.get(rid, {}):
            reaction.annotation["ec-code"] = list(context[rid]["ec-code"])
    for row in spec["reactions"].values():
        genes = row["genes"] if "genes" in row else [row["gene"]]
        for gene in genes:
            if gene["id"] not in model.genes:
                holder = Reaction(f"holder_{gene['id']}")
                model.add_reactions([holder])
                holder.gene_reaction_rule = gene["id"]
            g = model.genes.get_by_id(gene["id"])
            g.name = gene.get("name", g.name)
            g.annotation.update({k: v for k, v in gene.items() if k not in ("id", "name")})
    return model


out = {}
for change_id, (path, targets) in CHANGE_SETS.items():
    spec = json.loads(path.read_text())
    model = fixture(spec)
    apply_selection_gpr_restore(model, change_id=change_id)
    for rid in targets:
        out[rid] = [m.id for m in model.reactions.get_by_id(rid).metabolites]
print(json.dumps(out))
