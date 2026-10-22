"""Opt-in proton lock on R752, R408 and R750: keep their pH 7.3 H+ out of the proton/water balance gate."""

import json
from xml.sax.saxutils import unescape

from .config import MODEL

CURATION_PATH = MODEL.curation_file("thiamine_proton_lock.json")
REACTION_IDS = {"R752", "R408", "R750"}
# The plain-text form survives SBML; microspecies.balance_protons_and_water honours it.
LOCK = {"lock_proton_water_stoichiometry": "true"}
# Under the interim KEGG thiamine forms each reaction is off by exactly its pH 7.3 H+.
RESIDUAL = {"H": 1.0, "charge": 1.0}
XML_ENTITIES = {"&apos;": "'", "&quot;": '"'}
# libSBML drops a whole notes block holding < or & without raising, which would silently drop the lock.
UNSAFE_NOTE_CHARACTERS = set("<>&'\"")


def _unescaped(notes):
    # Non-string values (an in-memory dict lock, a bool) are compared as they are.
    return {key: unescape(value, XML_ENTITIES) if isinstance(value, str) else value for key, value in notes.items()}


def _notes(spec, row):
    return {**LOCK,
            "thiamine_proton_lock_curation_id": spec["curation_id"],
            "thiamine_proton_lock_status": spec["status"],
            "thiamine_proton_lock_authorization": spec["authorization"],
            "thiamine_proton_lock_evidence": row["evidence"],
            "thiamine_proton_lock_evidence_tier": row["evidence_tier"],
            "thiamine_proton_lock_rationale": row["rationale"],
            "thiamine_proton_lock_previous_state": row["previous_state"]}


def _residual(reaction):
    try:
        return {k: v for k, v in reaction.check_mass_balance().items() if abs(v) > 1e-9}
    except ValueError:
        return None


def _check_spec(spec):
    version = spec.get("schema_version")
    if (type(version) is not int or version != 1
            or spec.get("curation_id") != "THIAMINE-PROTON-LOCK-20261022"
            or spec.get("status") != "user_authorized_proton_water_lock"
            or set(spec.get("reactions", {})) != REACTION_IDS):
        raise ValueError("Unexpected thiamine proton lock scope/status")
    if spec.get("lock_note") != LOCK or any(row.get("residual") != RESIDUAL for row in spec["reactions"].values()):
        raise ValueError("Thiamine proton lock spec values differ from the reviewed change")
    texts = [spec.get("authorization")] + [row.get(key) for row in spec["reactions"].values()
                                           for key in ("evidence", "evidence_tier", "rationale", "previous_state")]
    if any(not isinstance(text, str) or set(text) & UNSAFE_NOTE_CHARACTERS for text in texts):
        raise ValueError("Thiamine proton lock note texts must be strings without < > & ' or \"")


def apply_thiamine_proton_lock(model, spec=None):
    """Add the reviewed proton lock notes; fail closed on any unexpected model state."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    _check_spec(spec)
    done = set()
    for rid, row in spec["reactions"].items():
        if rid not in model.reactions:
            raise ValueError(f"{rid} missing; no edits applied")
        reaction = model.reactions.get_by_id(rid)
        if (reaction.name != row["name"] or list(reaction.bounds) != row["bounds"]
                or reaction.gene_reaction_rule != row["gpr"]
                or {m.id: c for m, c in reaction.metabolites.items()} != row["stoichiometry"]):
            raise ValueError(f"{rid} precondition differs; no edits applied")
        # The lock is right only while the residual is exactly the pH 7.3 H+; anything else needs a new review.
        if _residual(reaction) != RESIDUAL:
            raise ValueError(f"{rid} residual differs from the reviewed one H+; no edits applied")
        notes, current = _notes(spec, row), _unescaped(reaction.notes)
        if any(key in current and current[key] != value for key, value in notes.items()):
            raise ValueError(f"Conflicting {rid} proton lock notes; no edits applied")
        if {k: v for k, v in current.items() if k not in notes} != row["notes_before"]:
            raise ValueError(f"{rid} notes outside this curation differ; no edits applied")
        present = [key in current for key in notes]
        if all(present):
            done.add(rid)
        elif any(present):
            raise ValueError(f"{rid} carries only part of this curation's notes; no edits applied")
    if done and done != REACTION_IDS:
        raise ValueError("Only part of the thiamine proton lock is present; no edits applied")
    if not done:
        for rid, row in spec["reactions"].items():
            reaction = model.reactions.get_by_id(rid)
            reaction.notes = {**reaction.notes, **_notes(spec, row)}
    return {"item": "thiamine proton lock", "curation_id": spec["curation_id"],
            "status": "already_correct" if done else "applied",
            "changed_reactions": [] if done else sorted(REACTION_IDS), "evidence_status": spec["status"]}
