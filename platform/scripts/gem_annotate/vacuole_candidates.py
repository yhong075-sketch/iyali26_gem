"""Explicit, bounded vacuole connections on the pinned E5 energy candidate."""

import copy
import hashlib
import json
from pathlib import Path

from cobra.io import read_sbml_model

from .config import MODEL, PLATFORM_ROOT, REPO_ROOT
from .energy_candidates import (
    export_candidate, model_definition, protected_definitions, signature,
)

SPEC_PATH = MODEL.curation_file("vacuole_connection_candidates.json")
R795_STOICHIOMETRY = {
    "m141[C_cy]": -1.0, "m32[C_cy]": -1.0, "m10[C_cy]": -2.0,
    "m143[C_cy]": 1.0, "m35[C_cy]": 1.0, "m1007[C_va]": 2.0,
}
TARGET_BOUNDS = {"R1363": [0, 0.04], "R795": [0, 0.04],
                 "R871": [0, 0.01], "R876": [0, 0.01]}
CHECK_IDS = set(TARGET_BOUNDS) | {
    "R2021", "R2029", "R2034", "R2039", "R2030", "R2035", "R2040",
}


def load_spec():
    return json.loads(SPEC_PATH.read_text())


def apply_vacuole_candidate(model, enabled=False, spec=None):
    """Check every dependency before applying only four absolute bound changes."""
    if not enabled:
        return {"status": "disabled", "items": []}
    spec = load_spec() if spec is None else spec
    if (spec.get("schema_version") != 1 or spec.get("enabled_by_default") is not False
            or spec.get("bounds") != TARGET_BOUNDS or set(spec.get("checks", {})) != CHECK_IDS):
        raise ValueError("Vacuole candidate schema or authorized scope differs")
    if protected_definitions(model) != spec["energy_protected_definitions"]:
        raise ValueError("E5 protected energy definitions differ")
    if signature(model.reactions.R795)["stoichiometry"] != R795_STOICHIOMETRY:
        raise ValueError("R795 is not the reviewed ATP/H2O/2Hcy to ADP/Pi/2Hva reaction")
    for rid, before in spec["checks"].items():
        target = copy.deepcopy(before)
        if rid in TARGET_BOUNDS:
            if before["bounds"] != [0, 0]:
                raise ValueError(f"{rid}: source is not closed")
            target["bounds"] = TARGET_BOUNDS[rid]
        if signature(model.reactions.get_by_id(rid)) not in (before, target):
            raise ValueError(f"{rid}: exact reaction signature differs; no edits applied")
        if rid in ("R2021", "R2029", "R2034", "R2039"):
            if before["bounds"][0] != 0 or before["bounds"][1] <= 0 or before["gpr"]:
                raise ValueError(f"{rid}: expected positive-only hydrolysis with empty GPR")
    for mid, expected in spec["species"].items():
        met = model.metabolites.get_by_id(mid)
        if {key: getattr(met, key) for key in expected} != expected:
            raise ValueError(f"{mid}: exact species signature differs; no edits applied")
    rows = []
    for rid, bounds in TARGET_BOUNDS.items():
        reaction = model.reactions.get_by_id(rid)
        before = signature(reaction)
        reaction.bounds = bounds
        rows.append({"reaction_id": rid, "before": before, "after": signature(reaction),
                     "status": "already_correct" if before["bounds"] == bounds else "applied"})
    return {"status": "complete", "curation_id": spec["curation_id"], "items": rows}


def build_candidate_file(source, output, enabled=False, spec=None):
    """Export a separate, reload-checked candidate; never overwrite an input/output."""
    spec = load_spec() if spec is None else spec
    source, output = Path(source), Path(output)
    if not all(path.resolve().is_relative_to(REPO_ROOT) for path in (source, output)):
        raise ValueError("Candidate source/output must remain in this project workspace")
    manifest = output.with_suffix(".build.json")
    if any(path.exists() or path.is_symlink() for path in (output, manifest)):
        raise FileExistsError("Candidate output or build manifest already exists")
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    if source_sha != spec["source_sha256"]:
        raise ValueError("Vacuole candidate input SHA differs from the reviewed E5")
    model = read_sbml_model(source)
    before = model_definition(model)
    notes_before = {r.id: copy.deepcopy(r.notes) for r in model.reactions}
    locks = protected_definitions(model)
    result = apply_vacuole_candidate(model, enabled=enabled, spec=spec)
    expected = copy.deepcopy(before)
    for row in result["items"]:
        expected["reactions"][row["reaction_id"]] = row["after"]
    if model_definition(model) != expected or {r.id: r.notes for r in model.reactions} != notes_before:
        raise ValueError("Vacuole candidate changed undeclared model fields or notes")
    loaded = export_candidate(model, output)
    if protected_definitions(loaded) != locks:
        raise ValueError("Vacuole candidate changed E5 energy locks")
    if {r.id: r.notes for r in loaded.reactions} != notes_before:
        raise ValueError("Vacuole candidate export changed reaction notes")
    if hashlib.sha256(source.read_bytes()).hexdigest() != source_sha:
        raise ValueError("Vacuole candidate source changed during build")
    implementation = [Path(__file__), Path(__file__).with_name("energy_candidates.py"),
                      Path(__file__).with_name("reaction_selection.py"), Path(__file__).with_name("sbml.py"),
                      PLATFORM_ROOT / "tools/build_vacuole_candidate.py"]
    record = {
        "candidate": "E5_vacuole_open" if enabled else "E5_vacuole_disabled",
        "enabled": enabled, "enabled_by_default": False,
        "source": str(source.resolve()), "source_sha256": source_sha,
        "output": str(output.resolve()), "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "spec_sha256": hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest(),
        "spec_file_sha256": hashlib.sha256(SPEC_PATH.read_bytes()).hexdigest(),
        "implementation_sha256": {str(path.relative_to(REPO_ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                  for path in implementation},
        "patch": result, "export_reload_definition_match": True,
        "undeclared_fields_unchanged": True, "energy_locks_unchanged": True,
        "reaction_notes_unchanged": True, "source_unchanged": True,
        "diagnostic_sources_added": False,
        "scope": "four bounded connections only; candidate, not a default or biological validation",
    }
    manifest.write_text(json.dumps(record, indent=2) + "\n")
    return record
