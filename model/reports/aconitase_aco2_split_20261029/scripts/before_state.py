"""Dump the baseline state of the aconitase reactions and their genes to before.json (read-only).

Targets: R1556, R216, R238 (rule to change), R239 (cytosolic aconitase, context) and the two genes. Context only,
not to be changed: every reaction annotated EC 4.2.1.3 / 4.2.1.36 / 4.2.1.114 or named *aconit* / *homocitrate*,
and every reaction of YALI1D11984g, YALI1E17935g and the LYS4-like YALI1E03387g.

Usage (from the repo root): python model/reports/aconitase_aco2_split_20261029/scripts/before_state.py MODEL.xml OUT.json
"""
import hashlib
import json
import sys
from html import unescape
from pathlib import Path

from cobra.io import read_sbml_model

TARGETS = ("R1556", "R216", "R238")
GENES = ("YALI1D11984g", "YALI1E17935g", "YALI1E03387g")
ECS = {"4.2.1.3", "4.2.1.36", "4.2.1.114"}

model_path, out = Path(sys.argv[1]), Path(sys.argv[2])
model = read_sbml_model(model_path)  # Path object: cobra rejects long str paths on macOS


def ec_codes(r):
    ec = r.annotation.get("ec-code", [])
    return set(ec if isinstance(ec, list) else [ec])


def reaction_state(r):
    return {
        "id": r.id, "name": r.name, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
        "compartments": sorted(r.compartments), "equation": r.build_reaction_string(),
        "equation_names": r.build_reaction_string(use_metabolite_names=True),
        "stoichiometry": {m.id: c for m, c in sorted(r.metabolites.items(), key=lambda x: x[0].id)},
        "species": {m.id: {k: getattr(m, k) for k in ("name", "formula", "charge", "compartment")}
                    for m in sorted(r.metabolites, key=lambda m: m.id)},
        "annotation": r.annotation,
        "notes": {k: unescape(v) if isinstance(v, str) else v for k, v in r.notes.items()},
    }


context = {r.id for r in model.reactions
           if ec_codes(r) & ECS or any(w in r.name.lower() for w in ("aconit", "homocitrate"))}
for g in GENES:
    context |= {r.id for r in model.genes.get_by_id(g).reactions}
context |= {"R239", "R24", "R395"} & {r.id for r in model.reactions}

state = {
    "baseline_model": {"path": str(model_path), "sha256": hashlib.sha256(model_path.read_bytes()).hexdigest()},
    "model_counts": {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)},
    "targets": list(TARGETS),
    "reactions": {rid: reaction_state(model.reactions.get_by_id(rid)) for rid in sorted(context)},
    "genes": {g: {"name": model.genes.get_by_id(g).name, "annotation": model.genes.get_by_id(g).annotation,
                  "notes": model.genes.get_by_id(g).notes,
                  "reactions": sorted(x.id for x in model.genes.get_by_id(g).reactions)}
              for g in GENES},
}
out.write_text(json.dumps(state, indent=1, sort_keys=True) + "\n")
for rid, r in state["reactions"].items():
    print(rid, "|", r["name"], "|", r["compartments"], "|", r["bounds"], "|", r["gpr"])
    print("    ", r["equation_names"])
for g, d in state["genes"].items():
    print(g, d["reactions"], d["annotation"].get("refseq"), d["annotation"].get("uniprot"))
