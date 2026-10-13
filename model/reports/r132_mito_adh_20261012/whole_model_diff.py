"""Whole-model comparison of the baseline and the R132 candidate (no solves).

Run from the repository root:
    .venv/bin/python -B model/reports/r132_mito_adh_20261012/whole_model_diff.py
"""
import hashlib
import json
from pathlib import Path

import cobra
from cobra.util.solver import linear_reaction_coefficients

REPO = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
BASE = REPO / "model/candidates/E5_coq9_respiratory_package_20261008.xml"
CAND = REPO / "model/candidates/E5_coq9_resp_pkg_r132_20261012.xml"
SHAS = {BASE: "0ac2addfa0cb6de2d89393e433766fb04325ffe42c2d131eaefaca3bcec29dba",
        CAND: "a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515"}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def reaction(r):
    return {"name": r.name, "stoichiometry": {m.id: c for m, c in r.metabolites.items()}, "bounds": list(r.bounds),
            "gpr": r.gene_reaction_rule, "subsystem": r.subsystem, "annotation": r.annotation, "notes": r.notes}


def snapshot(m):
    return {
        "reactions": {r.id: reaction(r) for r in m.reactions},
        "metabolites": {x.id: {"name": x.name, "formula": x.formula, "charge": x.charge, "compartment": x.compartment,
                               "annotation": x.annotation, "notes": x.notes} for x in m.metabolites},
        "genes": {g.id: {"name": g.name, "annotation": g.annotation, "notes": g.notes} for g in m.genes},
        "groups": {g.id: {"name": g.name, "kind": g.kind, "members": sorted(x.id for x in g.members)} for g in m.groups},
        "compartments": m.compartments, "objective": {r.id: c for r, c in linear_reaction_coefficients(m).items()},
        "objective_direction": m.objective.direction, "model": {"id": m.id, "name": m.name, "annotation": m.annotation,
                                                                "notes": m.notes},
    }


def diff(a, b):
    out = {}
    for section in a:
        if isinstance(a[section], dict) and section not in ("compartments", "objective", "model"):
            only_a, only_b = sorted(set(a[section]) - set(b[section])), sorted(set(b[section]) - set(a[section]))
            changed = {k: {f: {"before": a[section][k][f], "after": b[section][k][f]}
                           for f in a[section][k] if a[section][k][f] != b[section][k][f]}
                       for k in set(a[section]) & set(b[section]) if a[section][k] != b[section][k]}
            out[section] = {"removed": only_a, "added": only_b, "changed": changed}
        else:
            out[section] = None if a[section] == b[section] else {"before": a[section], "after": b[section]}
    return out


def main():
    for p, d in SHAS.items():
        assert sha(p) == d, p
    a, b = (cobra.io.read_sbml_model(str(p)) for p in (BASE, CAND))
    result = diff(snapshot(a), snapshot(b))
    # Genes whose reaction lists changed (expected: the 9 R132 genes lose R132 only).
    result["gene_reaction_lists_changed"] = {
        g.id: {"before": sorted(r.id for r in g.reactions), "after": sorted(r.id for r in b.genes.get_by_id(g.id).reactions)}
        for g in a.genes if sorted(r.id for r in g.reactions) != sorted(r.id for r in b.genes.get_by_id(g.id).reactions)}
    result["metabolite_usage_changed"] = {
        x.id: {"before": sorted(r.id for r in x.reactions), "after": sorted(r.id for r in b.metabolites.get_by_id(x.id).reactions)}
        for x in a.metabolites if sorted(r.id for r in x.reactions) != sorted(r.id for r in b.metabolites.get_by_id(x.id).reactions)}
    result["counts"] = {k: [len(getattr(a, k)), len(getattr(b, k))] for k in ("reactions", "metabolites", "genes", "groups")}
    result["inputs"] = {str(p.relative_to(REPO)): d for p, d in SHAS.items()}
    path = HERE / "whole_model_diff.json"
    path.open("x").write(json.dumps(result, indent=1, sort_keys=True, default=str) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("inputs",)}, indent=1, default=str)[:4000])


if __name__ == "__main__":
    main()
