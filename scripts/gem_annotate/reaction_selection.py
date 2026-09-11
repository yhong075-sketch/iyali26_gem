"""Apply the user's bounded selection of reaction fields from the metadata build."""

import html
import json

from cobra.core.gene import GPR

from .coq9 import boolean_key
from .config import REPO_ROOT

SELECTION_PATH = REPO_ROOT / "data" / "metadata_reaction_selection.json"


def reaction_fields(reaction):
    return {"stoichiometry": {m.id: c for m, c in reaction.metabolites.items()},
            "bounds": list(reaction.bounds), "gpr": reaction.gene_reaction_rule}


def _same(field, left, right):
    if field == "gpr":
        return boolean_key(GPR.from_string(left).body) == boolean_key(GPR.from_string(right).body)
    return left == right


def apply_metadata_reaction_selection(model, spec=None):
    """Set known old/target fields once; preserve each conflicting reaction in full."""
    if spec is None:
        spec = json.loads(SELECTION_PATH.read_text())
    rows = []
    for rid, rule in spec["reactions"].items():
        row = {"item": rid, "status": "conflict", "fields": []}
        rows.append(row)
        if rid not in model.reactions:
            row["reason"] = "reaction missing"
            continue
        reaction = model.reactions.get_by_id(rid)
        current = reaction_fields(reaction)
        errors = []
        for field, value in current.items():
            allowed = [rule["before"][field]]
            if field in rule.get("before_export", {}):
                allowed.append(rule["before_export"][field])
            if field in rule["fields"]:
                allowed.append(rule["after"][field])
            if not any(_same(field, value, expected) for expected in allowed):
                errors.append({"field": field, "actual": value, "allowed": allowed})
        for mid, expected in rule["species"].items():
            met = model.metabolites.get_by_id(mid) if mid in model.metabolites else None
            actual = {k: getattr(met, k) for k in expected} if met is not None else None
            if actual != expected and actual != rule.get("species_before_export", {}).get(mid, expected):
                errors.append({"field": "species", "id": mid, "actual": actual, "expected": expected})
        missing = GPR.from_string(rule["after"]["gpr"]).genes - {g.id for g in model.genes}
        if missing:
            errors.append({"field": "genes", "missing": sorted(missing)})
        if errors:
            row.update(reason="local precondition differs; current reaction preserved", errors=errors)
            continue
        if not rule["fields"]:
            row.update(status="already_correct", reason=rule["preserve_reason"])
            continue

        for field in rule["fields"]:
            target = rule["after"][field]
            changed = not _same(field, current[field], target)
            row["fields"].append({"field": field, "status": "applied" if changed else "already_correct"})
            if not changed:
                continue
            if field == "stoichiometry":
                reaction.add_metabolites({model.metabolites.get_by_id(mid): target.get(mid, 0)
                                         for mid in set(current[field]) | set(target)}, combine=False)
            elif field == "bounds":
                reaction.bounds = target
            else:
                reaction.gene_reaction_rule = target

        # Keep superseded evidence; never present it as validation of the selected rule.
        for key, value in rule.get("superseded_notes", {}).items():
            if key in reaction.notes and html.unescape(str(reaction.notes[key])) == html.unescape(str(value)):
                reaction.notes["metadata_previous_" + key] = reaction.notes.pop(key)
        reaction.notes.update({
            "metadata_reaction_selection": spec["selection_id"],
            "metadata_selected_fields": "; ".join(rule["fields"]),
            "metadata_selection_source_sha256": spec["source_metadata_sha256"],
            "metadata_selection_evidence_status": spec["biological_evidence_status"],
        })
        row["status"] = "applied" if any(f["status"] == "applied" for f in row["fields"]) else "already_correct"
    return {"selection_id": spec["selection_id"], "source_metadata_sha256": spec["source_metadata_sha256"],
            "field_counts": spec["field_counts"], "items": rows,
            "complete": all(row["status"] != "conflict" for row in rows)}
