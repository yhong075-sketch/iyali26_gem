"""Whole-model diff of two SBML models: every reaction, metabolite, gene, group, compartment and the objective."""
import json, sys
from pathlib import Path
import cobra


def rxn(r, model):
    return {"name": r.name, "subsystem": r.subsystem, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
            "stoichiometry": {m.id: c for m, c in r.metabolites.items()}, "annotation": dict(r.annotation),
            "notes": dict(r.notes), "objective": r.objective_coefficient}


def met(m):
    return {"name": m.name, "formula": m.formula, "charge": m.charge, "compartment": m.compartment,
            "annotation": dict(m.annotation), "notes": dict(m.notes)}


def gene(g):
    return {"name": g.name, "annotation": dict(g.annotation), "reactions": sorted(r.id for r in g.reactions)}


def table(model):
    return {"reactions": {r.id: rxn(r, model) for r in model.reactions},
            "metabolites": {m.id: met(m) for m in model.metabolites},
            "genes": {g.id: gene(g) for g in model.genes},
            "groups": {g.id: {"name": g.name, "kind": g.kind, "members": sorted(x.id for x in g.members)}
                       for g in model.groups},
            "compartments": {"compartments": dict(model.compartments)},
            "model": {"id": {"id": model.id, "name": model.name, "annotation": dict(model.annotation),
                             "notes": dict(model.notes)}}}


a, b = (table(cobra.io.read_sbml_model(Path(p))) for p in sys.argv[1:3])
diff = []
for kind in a:
    for key in sorted(set(a[kind]) | set(b[kind])):
        x, y = a[kind].get(key), b[kind].get(key)
        if x is None or y is None:
            diff.append({"kind": kind, "id": key, "field": "<presence>", "before": x is not None, "after": y is not None})
            continue
        for field in sorted(set(x) | set(y)):
            if x.get(field) != y.get(field):
                diff.append({"kind": kind, "id": key, "field": field, "before": x.get(field), "after": y.get(field)})
Path(sys.argv[3]).write_text(json.dumps(diff, indent=1, sort_keys=True, default=str) + "\n")
for d in diff:
    print(d["kind"], d["id"], d["field"])
print("entries", len(diff))
