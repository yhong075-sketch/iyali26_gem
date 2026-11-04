"""Dump the baseline state of R794, R795, their 14 genes and V0 e YALI1F25663g to before.json (read-only).

Usage (from the repo root): python model/reports/vatpase_gpr_a1_20261103/scripts/before_state.py MODEL.xml OUT.json
"""
import hashlib
import json
import sys
from html import unescape
from pathlib import Path

from cobra.io import read_sbml_model

NEW_GENE = "YALI1F25663g"
model_path, out = Path(sys.argv[1]), Path(sys.argv[2])
model = read_sbml_model(model_path)  # Path object: cobra rejects long str paths on macOS


def reaction_state(r):
    return {"id": r.id, "name": r.name, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
            "stoichiometry": {m.id: c for m, c in sorted(r.metabolites.items(), key=lambda x: x[0].id)},
            "species": {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")}
                        for m in sorted(r.metabolites, key=lambda m: m.id)},
            "annotation": r.annotation,
            "notes": {k: unescape(v) if isinstance(v, str) else v for k, v in r.notes.items()}}


reactions = [model.reactions.R794, model.reactions.R795]
genes = sorted({g for r in reactions for g in r.genes}, key=lambda g: g.id)
state = {
    "baseline_model": {"path": str(model_path), "sha256": hashlib.sha256(model_path.read_bytes()).hexdigest()},
    "reactions": {r.id: reaction_state(r) for r in reactions},
    "genes": {g.id: {"name": g.name, "annotation": g.annotation, "notes": g.notes,
                     "reactions": sorted(x.id for x in g.reactions)} for g in genes},
    "new_gene_present": NEW_GENE in model.genes,
    "gene_count": len(model.genes),
}
out.write_text(json.dumps(state, indent=1, sort_keys=True) + "\n")
print(len(state["genes"]), "genes;", NEW_GENE, "present:", state["new_gene_present"], "; total genes", len(model.genes))
