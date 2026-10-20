"""Record the baseline state of every R2190-RIBULOSE target in a56f3649 (read-only)."""
import hashlib
import json
import sys
from pathlib import Path

import cobra

model_path, out_path = Path(sys.argv[1]), Path(sys.argv[2])
model = cobra.io.read_sbml_model(str(model_path))


def reaction_state(rid):
    r = model.reactions.get_by_id(rid)
    return {"name": r.name, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
            "stoichiometry": {m.id: float(c) for m, c in sorted(r.metabolites.items(), key=lambda x: x[0].id)},
            "annotation": r.annotation, "notes": r.notes}


def metabolite_state(mid):
    m = model.metabolites.get_by_id(mid)
    return {"name": m.name, "formula": m.formula, "charge": m.charge, "compartment": m.compartment,
            "annotation": m.annotation, "notes": m.notes, "reactions": sorted(x.id for x in m.reactions)}


record = {
    "model": {"path": str(model_path), "sha256": hashlib.sha256(model_path.read_bytes()).hexdigest()},
    "counts": {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)},
    "target_reaction": {"R2190": reaction_state("R2190")},
    "target_species": {mid: metabolite_state(mid) for mid in ("m122[C_cy]", "m123[C_cy]", "m1823[C_cy]", "m454[C_cy]")},
    "new_species_id": "m2061[C_cy]",
    "new_species_id_present": "m2061[C_cy]" in model.metabolites,
    "ids_with_prefix_m2061": sorted(m.id for m in model.metabolites if m.id.split("[")[0] == "m2061"),
    "highest_numeric_metabolite": max(int(m.id[1:].split("[")[0]) for m in model.metabolites
                                      if m.id.startswith("m") and m.id[1:].split("[")[0].isdigit()),
    "free_ribulose_species": sorted(m.id for m in model.metabolites if m.name.lower().startswith("d-ribulose_")),
    "gene_YALI1F03452g_reactions": sorted(r.id for r in model.genes.get_by_id("YALI1F03452g").reactions),
    "loop_reactions_unchanged_by_this_change": {rid: reaction_state(rid) for rid in ("R1964", "R2041", "R2042", "R803", "R714")},
}
out_path.write_text(json.dumps(record, indent=1, ensure_ascii=False, sort_keys=True) + "\n")
print(json.dumps({k: record[k] for k in ("model", "counts", "new_species_id_present", "ids_with_prefix_m2061",
                                           "highest_numeric_metabolite", "free_ribulose_species",
                                           "gene_YALI1F03452g_reactions")}, indent=1))
print(json.dumps(record["target_reaction"], indent=1))
