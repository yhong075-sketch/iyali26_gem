"""Whole-model field-by-field diff of two SBML files (read-only).

Usage: python model_diff.py BASELINE.xml NEW.xml OUT.json
Compares every reaction, metabolite and gene field, gene-reaction membership, groups,
compartments, objective and model-level notes/annotation, and prints every difference.
"""
import json
import os
import sys

os.environ.setdefault("HOME", os.path.abspath("cache"))
from cobra.io import read_sbml_model


def reaction_fields(r):
    return {"name": r.name, "stoichiometry": {m.id: c for m, c in sorted(r.metabolites.items(), key=lambda x: x[0].id)},
            "bounds": list(r.bounds), "gpr": r.gene_reaction_rule, "subsystem": r.subsystem,
            "objective_coefficient": r.objective_coefficient,
            "annotation": r.annotation, "notes": r.notes}


def metabolite_fields(m):
    return {"name": m.name, "formula": m.formula, "charge": m.charge, "compartment": m.compartment,
            "annotation": m.annotation, "notes": m.notes}


def gene_fields(g):
    return {"name": g.name, "annotation": g.annotation, "notes": g.notes,
            "reactions": sorted(r.id for r in g.reactions)}


def table(model):
    return {
        "reactions": {r.id: reaction_fields(r) for r in model.reactions},
        "metabolites": {m.id: metabolite_fields(m) for m in model.metabolites},
        "genes": {g.id: gene_fields(g) for g in model.genes},
        "groups": {g.id: {"name": g.name, "members": sorted(x.id for x in g.members)} for g in model.groups},
        "model": {"id": model.id, "name": model.name, "compartments": model.compartments,
                  "notes": model.notes, "annotation": model.annotation,
                  "objective_direction": model.objective.direction},
    }


def diff(a, b):
    out = []
    for section in a:
        left, right = a[section], b[section]
        for key in sorted(set(left) | set(right)):
            if key not in left or key not in right:
                out.append({"section": section, "item": key, "change": "added" if key not in left else "removed"})
                continue
            if section == "model":
                if left[key] != right[key]:
                    out.append({"section": section, "item": key, "before": left[key], "after": right[key]})
                continue
            for field in sorted(set(left[key]) | set(right[key])):
                if left[key].get(field) != right[key].get(field):
                    out.append({"section": section, "item": key, "field": field,
                                "before": left[key].get(field), "after": right[key].get(field)})
    return out


if __name__ == "__main__":
    base, new, out = sys.argv[1:4]
    changes = diff(table(read_sbml_model(base)), table(read_sbml_model(new)))
    json.dump(changes, open(out, "w"), indent=1, default=str)
    for c in changes:
        print(c["section"], c["item"], c.get("field", c.get("change", "")))
    print("total differences:", len(changes))
