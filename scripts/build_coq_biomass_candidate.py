#!/usr/bin/env python3
"""Add explicitly parameterized Q9 growth dilution to a separate, pinned E5 XML."""

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import html
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cobra.io import read_sbml_model
from scripts.gem_annotate.energy_candidates import export_candidate, model_definition
from scripts.gem_annotate.execution import execution_limits

SPEC_PATH = ROOT / "data/coq_biomass_candidate.json"
BIOMASS = "biomass_C"
Q9, Q9H2 = "m468[C_mi]", "m471[C_mi]"
NOTE = "coq9_growth_dilution_candidate"


def pool_balance(model):
    """Sum the two existing steady-state rows; introduce no new constraint."""
    q, qh = (model.metabolites.get_by_id(mid) for mid in (Q9, Q9H2))
    return {r.id: c for r in model.reactions
            if (c := r.metabolites.get(q, 0) + r.metabolites.get(qh, 0)) != 0}


def apply_coq_biomass(model, alpha, alpha_source):
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not math.isfinite(alpha) or alpha <= 0:
        raise ValueError("alpha must be an explicit finite positive number in mmol/gDW")
    if not isinstance(alpha_source, str) or not alpha_source.strip():
        raise ValueError("alpha_source must identify the measurement or provisional assumption")
    bio = model.reactions.get_by_id(BIOMASS)
    q = model.metabolites.get_by_id(Q9)
    qh = model.metabolites.get_by_id(Q9H2)
    if q.compartment != "C_mi" or q.formula != "C54H82O4" or q.charge != 0:
        raise ValueError("Expected the existing neutral mitochondrial Q9 identity")
    original = {"R385": 1.0}
    target = {"R385": 1.0, BIOMASS: -alpha}
    if pool_balance(model) not in (original, target):
        raise ValueError("Unexpected combined Q9/Q9H2 pool balance; no edits applied")
    if bio.metabolites.get(q, 0) not in (0, -alpha):
        raise ValueError("Conflicting pre-existing Q9 biomass coefficient; no edits applied")
    if bio.metabolites.get(qh, 0) != 0:
        raise ValueError("Pre-existing Q9H2 biomass term; no edits applied")
    record = {"alpha_mmol_per_gDW": alpha, "alpha_source": alpha_source.strip(),
              "status": "candidate_not_physiologically_calibrated_by_this_build",
              "interpretation": "Q9 retained in newly formed biomass; not consumption per electron transfer",
              "pool_balance": target, "R385_lower_bound_changed": False}
    # Adding the delta also supports COBRA contexts when Q9 was absent.
    bio.add_metabolites({q: -alpha - bio.metabolites.get(q, 0)})
    bio.notes[NOTE] = json.dumps(record, sort_keys=True)
    return record


def build_candidate_file(source, output, alpha, alpha_source):
    source, output = Path(source), Path(output)
    manifest = output.with_suffix(".build.json")
    if not all(p.resolve().is_relative_to(ROOT) for p in (source, output)):
        raise ValueError("Source and output must remain in the project workspace")
    if output.suffix != ".xml":
        raise ValueError("Output must be a new .xml file")
    if any(p.exists() or p.is_symlink() for p in (output, manifest)):
        raise FileExistsError("Refusing to overwrite a model or build manifest")
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    spec = json.loads(SPEC_PATH.read_text())
    if (spec["schema_version"] != 1 or spec["biomass_reaction"] != BIOMASS
            or spec["metabolite"] != Q9 or spec["reduced_metabolite"] != Q9H2):
        raise ValueError("CoQ biomass curation scope differs")
    input_sha = sha(source)
    if input_sha != spec["source_sha256"]:
        raise ValueError("Source SHA differs from the pinned E5 candidate")
    implementation = [Path(__file__), SPEC_PATH,
                      ROOT / "scripts/gem_annotate/energy_candidates.py",
                      ROOT / "scripts/gem_annotate/reaction_selection.py",
                      ROOT / "scripts/gem_annotate/sbml.py",
                      ROOT / "scripts/gem_annotate/execution.py"]
    identities = {str(p.relative_to(ROOT)): sha(p) for p in implementation}
    model = read_sbml_model(source)
    before = model_definition(model)
    notes = {r.id: copy.deepcopy(r.notes) for r in model.reactions}
    patch = apply_coq_biomass(model, alpha, alpha_source)
    expected = copy.deepcopy(before)
    expected["reactions"][BIOMASS]["stoichiometry"][Q9] = -alpha
    notes[BIOMASS][NOTE] = model.reactions.get_by_id(BIOMASS).notes[NOTE]
    if model_definition(model) != expected:
        raise ValueError("Patch changed fields beyond the Q9 biomass coefficient")
    loaded = export_candidate(model, output)
    loaded_notes = {r.id: copy.deepcopy(r.notes) for r in loaded.reactions}
    # COBRA's XHTML reader can retain &quot; in the new JSON provenance note.
    # Verify its data exactly; every pre-existing note remains string-exact.
    if json.loads(html.unescape(loaded_notes[BIOMASS][NOTE])) != patch:
        raise ValueError("Export changed the CoQ provenance data")
    loaded_notes[BIOMASS][NOTE] = notes[BIOMASS][NOTE]
    if loaded_notes != notes or pool_balance(loaded) != patch["pool_balance"]:
        raise ValueError("Export changed reaction notes or the combined pool balance")
    if sha(source) != input_sha or any(sha(ROOT / p) != h for p, h in identities.items()):
        raise ValueError("Input or implementation changed during the build")
    record = {"created_utc": datetime.now(timezone.utc).isoformat(),
              "candidate": "E5_CoQ9_growth_dilution", "curation_id": spec["curation_id"],
              "source": str(source.resolve()), "source_sha256": input_sha,
              "output": str(output.resolve()), "output_sha256": sha(output),
              "implementation_sha256": identities, "patch": patch,
              "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "git_dirty": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True).splitlines(),
              "export_reload_definition_match": True, "non_target_fields_unchanged": True,
              "coq_note_semantic_roundtrip_match": True, "preexisting_notes_unchanged": True,
              "source_unchanged": True, "growth_and_essentiality_not_tested": True,
              "scope": "final-stage candidate on pinned E5; not a full rebuild or formal acceptance"}
    manifest.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    spec = json.loads(SPEC_PATH.read_text())
    parser.add_argument("--source", type=Path, default=ROOT / spec["source_path"])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--alpha", type=float, required=True, help="Q9 mmol per gDW of newly formed biomass")
    parser.add_argument("--alpha-source", required=True, help="Measured source or explicit provisional assumption")
    args = parser.parse_args()
    with execution_limits(no_solve=True, allow_network=False) as attempts:
        result = build_candidate_file(args.source, args.output, args.alpha, args.alpha_source)
        print(result["output"], "reload verified", attempts)


if __name__ == "__main__":
    main()
