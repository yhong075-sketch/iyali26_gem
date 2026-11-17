"""Record the starting state of every methylcitrate-cycle target in a model (read only, no LP)."""
import json
import sys
from html import unescape
from pathlib import Path

import cobra

model = cobra.io.read_sbml_model(Path(sys.argv[1]))
out = {"model": sys.argv[1], "reactions": {}, "metabolites": {}, "genes": {}, "groups": {}}
for rid in ("R_OBDHm", "R2075", "R693", "R25", "R95", "R552", "R2074", "R2076", "R2077", "R2112", "R834", "R2025",
            "R2113", "R1172"):
    if rid not in model.reactions:
        out["reactions"][rid] = None
        continue
    r = model.reactions.get_by_id(rid)
    out["reactions"][rid] = {
        "name": r.name, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule, "subsystem": r.subsystem,
        "stoichiometry": {m.id: c for m, c in r.metabolites.items()},
        "annotation": dict(r.annotation),
        "notes": {k: unescape(v) if isinstance(v, str) else v for k, v in r.notes.items()},
        "groups": sorted(g.name for g in model.groups if r in g.members),
    }
for mid in ("m42[C_mi]", "m73[C_mi]", "m27[C_mi]", "m78[C_mi]", "m44[C_mi]", "m30[C_mi]", "m28[C_mi]",
            "m2047[C_cy]", "m2046[C_cy]"):
    m = model.metabolites.get_by_id(mid)
    out["metabolites"][mid] = {"name": m.name, "formula": m.formula, "charge": m.charge,
                               "compartment": m.compartment, "reactions": sorted(r.id for r in m.reactions)}
for gid in ("YALI1E31947g", "YALI1F27556g", "YALI1D26360g", "YALI1B13179g", "YALI1D30618g"):
    g = model.genes.get_by_id(gid)
    out["genes"][gid] = {"name": g.name, "reactions": sorted(r.id for r in g.reactions)}
group = [g for g in model.groups if g.name == "Propanoate metabolism"]
out["groups"]["Propanoate metabolism"] = {"count": len(group),
                                          "members": sorted(x.id for x in group[0].members) if group else None}
Path(sys.argv[2]).write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
print("written", sys.argv[2])
