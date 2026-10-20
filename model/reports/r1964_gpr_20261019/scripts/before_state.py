"""Record the baseline state of every R1964-GPR target in a model file (read-only).

Usage: python before_state.py MODEL.xml OUT.json
Targets: reaction R1964 (fields, species), genes YALI1F03452g and YALI1E06616g, and the
other reactions of the R1964 species.
"""
import hashlib
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("HOME", os.path.abspath("cache"))
from cobra.io import read_sbml_model

GENES = ("YALI1F03452g", "YALI1E06616g")


def main(path, out):
    model = read_sbml_model(path)
    reaction = model.reactions.get_by_id("R1964")
    record = {
        "model": {"path": path, "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                  "reactions": len(model.reactions), "metabolites": len(model.metabolites),
                  "genes": len(model.genes)},
        "R1964": {
            "name": reaction.name, "bounds": list(reaction.bounds), "gpr": reaction.gene_reaction_rule,
            "stoichiometry": {m.id: c for m, c in reaction.metabolites.items()},
            "species": {m.id: {"name": m.name, "formula": m.formula, "charge": m.charge,
                               "compartment": m.compartment} for m in reaction.metabolites},
            "annotation": reaction.annotation, "notes": reaction.notes,
            "mass_balance_residual": reaction.check_mass_balance(),
        },
        "species_reactions": {m.id: sorted(r.id for r in m.reactions) for m in reaction.metabolites
                              if m.id in ("m1015[C_cy]", "m1823[C_cy]")},
        "genes": {g: ({"present": True, "name": model.genes.get_by_id(g).name,
                       "annotation": model.genes.get_by_id(g).annotation,
                       "reactions": sorted(r.id for r in model.genes.get_by_id(g).reactions)}
                      if g in model.genes else {"present": False}) for g in GENES},
    }
    Path(out).write_text(json.dumps(record, indent=1, default=str) + "\n")
    print(json.dumps(record, indent=1, default=str))


if __name__ == "__main__":
    main(*sys.argv[1:3])
