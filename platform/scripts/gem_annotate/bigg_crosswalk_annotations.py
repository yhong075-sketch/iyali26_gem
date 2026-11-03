"""Opt-in reviewed annotation values from Andrew McNaughton's BiGG crosswalk (mcnaugh_pass), legacy IDs kept."""

import collections
import csv
import hashlib
import json

from .config import MODEL

MANIFEST_PATH = MODEL.curation_file("bigg_crosswalk_annotation_decisions_20261103.json")
DECISIONS_PATH = MODEL.curation_file("bigg_crosswalk_annotation_decisions_20261103.csv")
KEYS = {"metabolite": {"bigg.metabolite", "metanetx.chemical", "chebi", "kegg.compound"},
        "reaction": {"bigg.reaction", "metanetx.reaction"}}
ACTIONS = {"write_his", "drop_wrong"}


def load_decisions():
    """Return the manifest and decision rows; the rows must match the manifest's SHA-256."""
    manifest = json.loads(MANIFEST_PATH.read_text())
    data = DECISIONS_PATH.read_bytes()
    if hashlib.sha256(data).hexdigest() != manifest.get("decisions_sha256"):
        raise ValueError("BiGG crosswalk decisions file differs from its manifest; no edits applied")
    return manifest, list(csv.DictReader(data.decode().splitlines()))


def _values(obj, key):
    value = obj.annotation.get(key)
    return [] if value is None else list(value) if isinstance(value, list) else [value]


def apply_bigg_crosswalk_annotations(model, manifest=None, rows=None):
    """Set each reviewed annotation key to its decided value; on any stale or unexpected row, edit nothing."""
    if manifest is None or rows is None:
        manifest, rows = load_decisions()
    if (manifest.get("schema_version") != 1
            or manifest.get("decision_set_id") != "BIGG-XW-ANNOT-20261103"
            or manifest.get("status") != "user_authorized"
            or manifest.get("rows") != len(rows)):
        raise ValueError("Unexpected BiGG crosswalk decision set scope/status")
    targets, seen = [], set()
    for row in rows:
        entity, key = row["entity"], row["key"]
        # A row whose before equals its after could never be told apart from an applied one.
        if (key not in KEYS.get(entity, ()) or row["action"] not in ACTIONS or not row["after"]
                or row["before"] == row["after"]):
            raise ValueError(f"{row['decision_id']}: unexpected entity/key/action; no edits applied")
        if (entity, row["entity_id"], key) in seen:
            raise ValueError(f"{row['decision_id']}: duplicate target; no edits applied")
        seen.add((entity, row["entity_id"], key))
        collection = model.metabolites if entity == "metabolite" else model.reactions
        if row["entity_id"] not in collection:
            raise ValueError(f"{row['decision_id']}: {row['entity_id']} missing; no edits applied")
        obj = collection.get_by_id(row["entity_id"])
        before, after, current = row["before"].split(";"), row["after"].split(";"), _values(obj, key)
        targets.append((row["decision_id"], obj, key, after,
                        "before" if current == before else "after" if current == after else "stale"))
    states = collections.Counter(state for *_, state in targets)
    if set(states) not in ({"before"}, {"after"}):
        # Name the stale rows, or in a partly applied model the rows in the rarer state.
        odd = "stale" if states["stale"] else min(("before", "after"), key=lambda state: states[state])
        examples = [decision for decision, *_, state in targets if state == odd][:5]
        raise ValueError(f"BiGG crosswalk targets differ from the reviewed state {dict(states)}, e.g. {examples}; "
                         "no edits applied")
    if set(states) == {"before"}:
        for _, obj, key, after, _ in targets:
            # Rebind, never mutate: Model.copy() can share annotation dicts with the original.
            obj.annotation = {**obj.annotation, key: after[0] if len(after) == 1 else after}
    return {"item": manifest["decision_set_id"], "status": "applied" if set(states) == {"before"} else "already_correct",
            "rows": len(targets), "entities": len({id(obj) for _, obj, *_ in targets}),
            "actions": dict(collections.Counter(row["action"] for row in rows))}
