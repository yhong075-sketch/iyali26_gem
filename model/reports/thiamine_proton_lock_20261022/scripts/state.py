"""State of the three lock targets, any existing locks, and a whole-model proton-gate dry run (read-only, no LPs).

Usage (from platform/): PYTHONPATH=$PWD python ../model/reports/thiamine_proton_lock_20261022/scripts/state.py \
    MODEL.xml OUT.json
The gate runs on an in-memory copy read from the file, so this shows what it would do to the saved model.
"""
import hashlib
import json
import sys
import warnings

from cobra.io import read_sbml_model

from scripts.gem_annotate.microspecies import balance_protons_and_water

TARGETS = ["R752", "R408", "R750"]
LOCK = "lock_proton_water_stoichiometry"

path, out = sys.argv[1], sys.argv[2]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = read_sbml_model(path)


def balance(reaction):
    return {k: v for k, v in reaction.check_mass_balance().items() if abs(v) > 1e-9}


state = {
    "model": path.split("/")[-1],
    "sha256": hashlib.sha256(open(path, "rb").read()).hexdigest(),
    "targets": {rid: {"name": r.name, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
                      "stoichiometry": {m.id: c for m, c in sorted(r.metabolites.items(), key=lambda x: x[0].id)},
                      "notes": r.notes, "balance": balance(r)}
                for rid in TARGETS for r in [model.reactions.get_by_id(rid)]},
    # Any reaction already carrying a lock in either form, as read back from the file.
    "existing_locks": {r.id: {k: str(v)[:200] for k, v in r.notes.items() if LOCK in k or LOCK in str(v)}
                       for r in model.reactions if any(LOCK in k or LOCK in str(v) for k, v in r.notes.items())},
}
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    report = balance_protons_and_water(model)
state["gate_dry_run"] = {"changed": {c["reaction_id"]: c["additions"] for c in report["changes"]},
                         "skipped_curated_lock_reaction_ids": report["skipped_curated_lock_reaction_ids"],
                         "rejected_count": report["rejected_reactions"]}
with open(out, "w") as handle:
    json.dump(state, handle, indent=1)
    handle.write("\n")
print(json.dumps({"sha256": state["sha256"][:8], "existing_locks": state["existing_locks"],
                  "targets_balance": {k: v["balance"] for k, v in state["targets"].items()},
                  "gate_dry_run": state["gate_dry_run"]}, indent=1))
