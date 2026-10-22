"""Record the state of every target of the R2193 iron-import GPR change (read-only).

Usage: python before_state.py MODEL.xml OUT.json
"""
import hashlib
import json
import sys
import warnings

from cobra.io import read_sbml_model

TARGET = "R2193"
OLD_GENE = "YALI1F07747g"
NEW_GENES = ["YALI1D09277g", "YALI1D09308g", "YALI1D08684g", "YALI1D08564g"]
CONTEXT = ["R1189", "R863", "R302", "R727"]

path, out = sys.argv[1], sys.argv[2]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = read_sbml_model(path)


def reaction(rid):
    r = model.reactions.get_by_id(rid)
    return {"name": r.name, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
            "stoichiometry": {m.id: c for m, c in sorted(r.metabolites.items(), key=lambda x: x[0].id)},
            "species": {m.id: {"formula": m.formula, "charge": m.charge, "compartment": m.compartment}
                        for m in sorted(r.metabolites, key=lambda m: m.id)},
            "annotation": r.annotation, "notes": r.notes}


old = model.genes.get_by_id(OLD_GENE)
state = {
    "model": path.split("/")[-1],
    "sha256": hashlib.sha256(open(path, "rb").read()).hexdigest(),
    "counts": {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)},
    "target": reaction(TARGET),
    "context": {rid: reaction(rid) for rid in CONTEXT},
    "old_gene": {"id": OLD_GENE, "name": old.name, "annotation": old.annotation,
                 "reactions": sorted(r.id for r in old.reactions)},
    "new_genes_present": {g: g in model.genes for g in NEW_GENES},
    "orphan_genes": sorted(g.id for g in model.genes if not g.reactions),
}
with open(out, "w") as handle:
    json.dump(state, handle, indent=1)
    handle.write("\n")
print(json.dumps(state, indent=1)[:3000])
