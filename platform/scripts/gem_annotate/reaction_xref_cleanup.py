"""Opt-in removal of reaction cross-references whose MetaNetX equation does not match the reaction (2026-11-13).

The 2026-11-12 annotation audit compared every reaction's compounds with the equation of each cross-referenced MNXR.
On 62 reactions outside the lipid-unlump line, the values copied from a wrong MNXR, or left from the annotator's
first MNXR after a crosswalk rewrite, fail that comparison. The reviewed decisions are data in
``reaction_xref_cleanup_decisions_20261113.csv`` (pinned by its ``.json`` manifest). This step only removes the listed
values (each row removes its whole key) and adds one note per reaction. It never touches a key the BiGG crosswalk
decided or an excluded reaction. It needs the crosswalk applied first: 11 targets share reactions with crosswalk
decisions, and 9 notes state that the crosswalk MNXR is kept, so those decided values must hold their reviewed
``after`` state or the step refuses to run.
"""

import collections
import csv
import hashlib
import json
from functools import partial

from cobra.util.context import get_context

from .bigg_crosswalk_annotations import load_decisions as load_crosswalk_decisions
from .config import MODEL

MANIFEST_PATH = MODEL.curation_file("reaction_xref_cleanup_decisions_20261113.json")
DECISIONS_PATH = MODEL.curation_file("reaction_xref_cleanup_decisions_20261113.csv")
DECISION_SET_ID = "REACTION-XREF-CLEANUP-20261113"
KEYS = frozenset({"metanetx.reaction", "mnx", "bigg.reaction", "biggr", "kegg.reaction", "keggr", "rhea", "rh",
                  "seed.reaction", "seedr", "metacyc.reaction", "metacycr", "sabiork.reaction", "sabiorkr",
                  "vmhreaction", "vmhr"})
NOTES_KEY = "xref_cleanup"


def load_decisions():
    """Return the manifest and decision rows; the rows must match the manifest's SHA-256."""
    manifest = json.loads(MANIFEST_PATH.read_text())
    data = DECISIONS_PATH.read_bytes()
    if hashlib.sha256(data).hexdigest() != manifest.get("decisions_sha256"):
        raise ValueError("Reaction xref cleanup decisions differ from their manifest; no edits applied")
    return manifest, list(csv.DictReader(data.decode().splitlines()))


def _values(obj, key):
    value = obj.annotation.get(key)
    return [] if value is None else list(value) if isinstance(value, list) else [value]


def _split(text):
    return text.split(";") if text else []


def apply_reaction_xref_cleanup(model, manifest=None, rows=None, crosswalk_rows=None):
    """Remove each reviewed value and add the reaction's note; on any stale or unexpected row, edit nothing."""
    if manifest is None or rows is None:
        manifest, rows = load_decisions()
    if crosswalk_rows is None:
        _, crosswalk_rows = load_crosswalk_decisions()
    notes = manifest.get("notes")
    if (manifest.get("schema_version") != 1
            or manifest.get("decision_set_id") != DECISION_SET_ID
            or manifest.get("status") != "user_authorized"
            or manifest.get("rows") != len(rows)
            or set(manifest.get("keys", ())) != KEYS
            or manifest.get("notes_key") != NOTES_KEY
            or not isinstance(notes, dict)):
        raise ValueError("Unexpected reaction xref cleanup decision set scope/status")
    crosswalk = [r for r in crosswalk_rows if r["entity"] == "reaction"]
    decided = {(r["entity_id"], r["key"]) for r in crosswalk}
    excluded = set(manifest.get("excluded_reactions", {}))
    targets, seen = [], set()
    for row in rows:
        rid, key = row["entity_id"], row["key"]
        before, removed = _split(row["before"]), _split(row["removed"])
        # Each row removes its whole key; a partial removal is not part of this decision set.
        if (row["entity"] != "reaction" or key not in KEYS or row["action"] != "remove_unsupported"
                or not removed or row["after"] or sorted(removed) != sorted(before)):
            raise ValueError(f"{row['decision_id']}: unexpected entity/key/action or inconsistent values; "
                             "no edits applied")
        if (rid, key) in seen:
            raise ValueError(f"{row['decision_id']}: duplicate target; no edits applied")
        if (rid, key) in decided:
            # The crosswalk owns that key's value for this reaction.
            raise ValueError(f"{row['decision_id']}: {rid} {key} is a BiGG crosswalk decision; no edits applied")
        if rid in excluded or rid not in notes:
            raise ValueError(f"{row['decision_id']}: {rid} is excluded or has no note; no edits applied")
        seen.add((rid, key))
        if rid not in model.reactions:
            raise ValueError(f"{row['decision_id']}: {rid} missing; no edits applied")
        reaction = model.reactions.get_by_id(rid)
        current = sorted(_values(reaction, key))
        state = "before" if current == sorted(before) else "after" if not current else "stale"
        targets.append((row["decision_id"], reaction, key, state))
    reactions = sorted({reaction for _, reaction, *_ in targets}, key=lambda r: r.id)
    if {r.id for r in reactions} != set(notes) or len(reactions) != manifest.get("reactions"):
        raise ValueError("Reaction xref cleanup rows and notes name different reactions; no edits applied")
    # The reviewed state includes the crosswalk's values on the targets it also decided (the notes rely on them).
    target_ids = {r.id for r in reactions}
    unreviewed = [f"{r['entity_id']} {r['key']}" for r in crosswalk if r["entity_id"] in target_ids
                  and sorted(_values(model.reactions.get_by_id(r["entity_id"]), r["key"])) != sorted(_split(r["after"]))]
    if unreviewed:
        raise ValueError(f"BiGG crosswalk values on cleanup targets are not in their reviewed state, e.g. "
                         f"{unreviewed[:5]}; apply --bigg-crosswalk-annotations first; no edits applied")
    note_states = {r.id: ("before" if NOTES_KEY not in r.notes else "after" if r.notes[NOTES_KEY] == notes[r.id]
                          else "stale") for r in reactions}
    states = collections.Counter([state for *_, state in targets] + list(note_states.values()))
    if set(states) not in ({"before"}, {"after"}):
        # Name the stale rows, or in a partly applied model the rows in the rarer state.
        odd = "stale" if states["stale"] else min(("before", "after"), key=lambda state: states[state])
        examples = [d for d, *_, s in targets if s == odd] + [r for r, s in note_states.items() if s == odd]
        raise ValueError(f"Reaction xref cleanup targets differ from the reviewed state {dict(states)}, "
                         f"e.g. {examples[:5]}; no edits applied")
    record = {"item": DECISION_SET_ID, "rows": len(targets), "reactions": len(reactions),
              "values_removed": sum(len(_split(r["removed"])) for r in rows)}
    if set(states) == {"after"}:
        return {**record, "status": "already_correct"}
    for reaction in reactions:
        context = get_context(reaction)
        if context:
            # Annotation and notes are not history-tracked by cobra; restore both on context exit.
            context(partial(setattr, reaction, "annotation", reaction.annotation))
            context(partial(setattr, reaction, "notes", reaction.notes))
    for _, reaction, key, _ in targets:
        # Rebind, never mutate: Model.copy() can share annotation dicts with the original.
        reaction.annotation = {k: v for k, v in reaction.annotation.items() if k != key}
    for reaction in reactions:
        reaction.notes = {**reaction.notes, NOTES_KEY: notes[reaction.id]}
    return {**record, "status": "applied"}
