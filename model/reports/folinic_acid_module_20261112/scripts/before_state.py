"""Record the starting state of every folinic-acid-module target in a model (read only)."""
import json, sys
from html import unescape
from pathlib import Path
import cobra

model = cobra.io.read_sbml_model(Path(sys.argv[1]))
out = {"model": sys.argv[1], "reactions": {}, "metabolites": {}, "genes": {}}
for rid in ("R748", "R990", "R64", "R811", "R991", "R358", "R359", "R543", "R_THFAT"):
    if rid not in model.reactions:
        out["reactions"][rid] = None
        continue
    r = model.reactions.get_by_id(rid)
    out["reactions"][rid] = {
        "name": r.name, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule, "subsystem": r.subsystem,
        "stoichiometry": {m.id: c for m, c in r.metabolites.items()},
        "species": {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")} for m in r.metabolites},
        "annotation": dict(r.annotation),
        "notes": {k: unescape(v) if isinstance(v, str) else v for k, v in r.notes.items()},
        "groups": sorted(g.name for g in model.groups if r in g.members),
    }
for mid in ("m159[C_cy]", "m160[C_cy]", "m32[C_cy]", "m10[C_cy]", "m969[C_mi]", "m803[C_mi]"):
    m = model.metabolites.get_by_id(mid)
    out["metabolites"][mid] = {"name": m.name, "formula": m.formula, "charge": m.charge,
                               "compartment": m.compartment, "reactions": sorted(r.id for r in m.reactions)}
for gid in ("YALI1F04082g", "YALI1D28569g", "YALI1E19557g", "YALI1D01489g"):
    g = model.genes.get_by_id(gid)
    out["genes"][gid] = {"annotation": dict(g.annotation), "reactions": sorted(r.id for r in g.reactions)}
Path(sys.argv[2]).write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
