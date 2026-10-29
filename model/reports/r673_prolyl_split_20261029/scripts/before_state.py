"""Dump the baseline state of R673 and its 2 genes to before.json (read-only).

Usage (from the repo root): python model/reports/r673_prolyl_split_20261029/scripts/before_state.py MODEL.xml OUT.json
"""
import hashlib
import json
import sys
from html import unescape
from pathlib import Path

from cobra.io import read_sbml_model

model_path, out = Path(sys.argv[1]), Path(sys.argv[2])
model = read_sbml_model(model_path)  # Path object: cobra rejects long str paths on macOS
r = model.reactions.R673
state = {
    "baseline_model": {"path": str(model_path), "sha256": hashlib.sha256(model_path.read_bytes()).hexdigest()},
    "reaction": {
        "id": r.id, "name": r.name, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
        "stoichiometry": {m.id: c for m, c in sorted(r.metabolites.items(), key=lambda x: x[0].id)},
        "species": {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")}
                    for m in sorted(r.metabolites, key=lambda m: m.id)},
        "annotation": r.annotation,
        "notes": {k: unescape(v) if isinstance(v, str) else v for k, v in r.notes.items()},
    },
    "genes": {g.id: {"name": g.name, "annotation": g.annotation, "notes": g.notes,
                     "reactions": sorted(x.id for x in g.reactions)}
              for g in sorted(r.genes, key=lambda g: g.id)},
}
out.write_text(json.dumps(state, indent=1, sort_keys=True) + "\n")
print(len(state["genes"]), "genes;", state["reaction"]["gpr"])
