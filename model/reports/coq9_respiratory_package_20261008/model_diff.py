"""Whole-model diff of two SBML files: every reaction, metabolite, gene field, objective and groups."""
import json, sys
from html import unescape
from cobra.io import read_sbml_model
from cobra.util.solver import linear_reaction_coefficients

def norm(v):
    if isinstance(v, str): return unescape(v)
    if isinstance(v, list) and len(v) == 1: return norm(v[0])
    if isinstance(v, list): return sorted(norm(x) for x in v)
    if isinstance(v, dict): return {k: norm(x) for k, x in v.items()}
    return v

def rx(r):
    return {"name": r.name, "stoichiometry": {m.id: c for m, c in r.metabolites.items()},
            "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
            "annotation": norm(r.annotation), "notes": norm(r.notes), "subsystem": r.subsystem}

def mt(m):
    return {"name": m.name, "formula": m.formula, "charge": m.charge, "compartment": m.compartment,
            "annotation": norm(m.annotation), "notes": norm(m.notes)}

def gn(g):
    return {"name": g.name, "annotation": norm(g.annotation), "notes": norm(g.notes)}

def diff(a_path, b_path):
    a, b = read_sbml_model(a_path), read_sbml_model(b_path)
    out = {"counts": {"a": [len(a.reactions), len(a.metabolites), len(a.genes)],
                      "b": [len(b.reactions), len(b.metabolites), len(b.genes)]}}
    for kind, fa, fb, getter in (("reactions", a.reactions, b.reactions, rx),
                                 ("metabolites", a.metabolites, b.metabolites, mt),
                                 ("genes", a.genes, b.genes, gn)):
        ia, ib = {x.id: getter(x) for x in fa}, {x.id: getter(x) for x in fb}
        out[kind] = {"only_in_a": sorted(set(ia) - set(ib)), "only_in_b": sorted(set(ib) - set(ia)), "changed": {}}
        for i in sorted(set(ia) & set(ib)):
            if ia[i] != ib[i]:
                d = {}
                for f in ia[i]:
                    if ia[i][f] != ib[i][f]:
                        if f in ("stoichiometry", "notes", "annotation") and isinstance(ia[i][f], dict):
                            d[f] = {k: [ia[i][f].get(k), ib[i][f].get(k)] for k in set(ia[i][f]) | set(ib[i][f]) if ia[i][f].get(k) != ib[i][f].get(k)}
                        else:
                            d[f] = [ia[i][f], ib[i][f]]
                out[kind]["changed"][i] = d
    out["objective"] = [{r.id: c for r, c in linear_reaction_coefficients(a).items()},
                        {r.id: c for r, c in linear_reaction_coefficients(b).items()}]
    out["objective_direction"] = [a.objective.direction, b.objective.direction]
    out["compartments"] = [a.compartments, b.compartments]
    ga = {g.id: sorted(m.id for m in g.members) for g in a.groups}
    gb = {g.id: sorted(m.id for m in g.members) for g in b.groups}
    out["groups_changed"] = {k: [ga.get(k), gb.get(k)] for k in set(ga) | set(gb) if ga.get(k) != gb.get(k)}
    return out

if __name__ == "__main__":
    result = diff(sys.argv[1], sys.argv[2])
    json.dump(result, open(sys.argv[3], "w"), indent=1, default=str)
    for kind in ("reactions", "metabolites", "genes"):
        k = result[kind]
        print(kind, "only_in_a:", k["only_in_a"], "only_in_b:", k["only_in_b"], "changed:", len(k["changed"]))
        for i, d in k["changed"].items():
            print("  ", i, json.dumps(d, default=str)[:400])
    print("counts", result["counts"], "objective same:", result["objective"][0] == result["objective"][1],
          "groups changed:", list(result["groups_changed"])[:10], "compartments same:", result["compartments"][0] == result["compartments"][1])
