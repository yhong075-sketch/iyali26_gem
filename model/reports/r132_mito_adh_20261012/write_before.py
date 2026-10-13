"""Record the baseline state of every Task B target (no solves).

Run from the repository root:
    .venv/bin/python -B model/reports/r132_mito_adh_20261012/write_before.py
"""
import hashlib
import json
from pathlib import Path

import cobra

REPO = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
MODEL = REPO / "model/candidates/E5_coq9_respiratory_package_20261008.xml"
EXPECTED = "0ac2addfa0cb6de2d89393e433766fb04325ffe42c2d131eaefaca3bcec29dba"
REACTIONS = ("R132", "R1114", "R997", "R130", "R1377", "R_ALCDH_nadp_hi", "R135", "R136", "R2070")


def reaction_state(r):
    return {"name": r.name, "stoichiometry": {m.id: float(c) for m, c in sorted(r.metabolites.items(), key=lambda x: x[0].id)},
            "bounds": list(r.bounds), "gpr": r.gene_reaction_rule, "annotation": dict(r.annotation), "notes": dict(r.notes)}


def main():
    digest = hashlib.sha256(MODEL.read_bytes()).hexdigest()
    assert digest == EXPECTED, digest
    model = cobra.io.read_sbml_model(str(MODEL))
    r132 = model.reactions.R132
    state = {
        "model": {"path": str(MODEL.relative_to(REPO)), "sha256": digest,
                  "counts": {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)},
                  "orphan_genes": sorted(g.id for g in model.genes if not g.reactions)},
        "reactions": {rid: reaction_state(model.reactions.get_by_id(rid)) for rid in REACTIONS},
        "metabolites": {m.id: {"name": m.name, "formula": m.formula, "charge": m.charge, "compartment": m.compartment,
                               "reactions": sorted(x.id for x in m.reactions)} for m in sorted(r132.metabolites, key=lambda m: m.id)},
        "r132_genes": {g.id: {"annotation": dict(g.annotation), "reactions": sorted(x.id for x in g.reactions)}
                       for g in sorted(r132.genes, key=lambda g: g.id)},
    }
    path = HERE / "before.json"
    with path.open("x") as f:
        f.write(json.dumps(state, indent=1, sort_keys=True) + "\n")
    print(path, hashlib.sha256(path.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
