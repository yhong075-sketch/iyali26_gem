"""Opt-in, signature-checked ATP repair candidates; never changes the default build."""

import hashlib
import html
import json
from pathlib import Path

from cobra.io import read_sbml_model
from cobra.util.solver import linear_reaction_coefficients

from .config import REPO_ROOT
from .reaction_selection import reaction_fields
from .sbml import write_deterministic_sbml_model

SPEC_PATH = REPO_ROOT / "data" / "energy_candidate_repairs.json"
LOCK_NOTE = "energy_candidate_protected_definition"
VARIANTS = ("E0", "E1", "E2", "E3", "E4", "E5")


def load_spec():
    return json.loads(SPEC_PATH.read_text())


def signature(reaction):
    return {"name": reaction.name, **reaction_fields(reaction)}


def protected_definitions(model):
    """Validate all persisted locks before metadata can edit even one reaction."""
    protected = {}
    for reaction in model.reactions:
        if LOCK_NOTE not in reaction.notes:
            continue
        lock = json.loads(html.unescape(str(reaction.notes[LOCK_NOTE])))
        if (lock.get("schema_version") != 1 or lock.get("reaction_id") != reaction.id
                or lock.get("definition") != reaction_fields(reaction)):
            raise ValueError(f"{reaction.id}: protected energy candidate definition changed")
        for mid, props in lock["species"].items():
            met = model.metabolites.get_by_id(mid)
            if {k: getattr(met, k) for k in props} != props:
                raise ValueError(f"{reaction.id}: protected species {mid} changed")
        protected[reaction.id] = lock
    return protected


def target_definition(rule, variant):
    target = json.loads(json.dumps(rule["before"]))
    if variant not in rule.get("variants", VARIANTS):
        return target
    if variant in ("E1", "E3", "E4", "E5") and "direction_bounds" in rule:
        target["bounds"] = rule["direction_bounds"]
    if variant in ("E2", "E3", "E4", "E5") and "chemistry_stoichiometry" in rule:
        target["stoichiometry"] = rule["chemistry_stoichiometry"]
    return target


def apply_energy_candidate(model, variant="E0", spec=None):
    """Preflight the whole patch, then set absolute fields once; E0 is a no-op."""
    spec = load_spec() if spec is None else spec
    if variant not in VARIANTS or spec.get("schema_version") != 1:
        raise ValueError("Unknown energy candidate/schema")
    if variant == "E0":
        return {"variant": variant, "items": [], "status": "disabled"}
    locks = protected_definitions(model)
    planned = []
    for rid, rule in spec["reactions"].items():
        reaction = model.reactions.get_by_id(rid)
        target = target_definition(rule, variant)
        current = signature(reaction)
        allowed = [rule["before"], target]
        if rid in locks:
            previous = locks[rid]["variant"]
            predecessors = {"E3": ("E1", "E2"), "E4": ("E1", "E2", "E3"),
                            "E5": ("E1", "E2", "E3", "E4")}
            if previous != variant and previous not in predecessors.get(variant, ()):
                raise ValueError(f"{rid}: cannot remove an applied candidate repair")
            allowed.append(target_definition(rule, previous))
        if current not in allowed:
            raise ValueError(f"{rid}: exact energy reaction signature differs; no edits applied")
        for mid, props in rule["species"].items():
            met = model.metabolites.get_by_id(mid)
            if {k: getattr(met, k) for k in props} != props:
                raise ValueError(f"{rid}: exact species signature differs for {mid}; no edits applied")
        if target != rule["before"]:
            planned.append((reaction, rule, target))
    rows = []
    for reaction, rule, target in planned:
        before = signature(reaction)
        reaction.add_metabolites({model.metabolites.get_by_id(mid): target["stoichiometry"].get(mid, 0)
                                 for mid in set(before["stoichiometry"]) | set(target["stoichiometry"])},
                                combine=False)
        reaction.bounds = target["bounds"]
        lock = {"schema_version": 1, "reaction_id": reaction.id, "variant": variant,
                "curation_id": spec["curation_id"], "definition": reaction_fields(reaction),
                "species": rule["species"]}
        reaction.notes.update({LOCK_NOTE: json.dumps(lock, sort_keys=True),
                               "energy_candidate_status": spec["status"],
                               "energy_candidate_evidence": "; ".join(rule["sources"]),
                               "energy_candidate_rationale": rule["rationale"],
                               "energy_candidate_precedence": "explicit candidate > metadata chemistry/bounds"})
        rows.append({"reaction_id": reaction.id, "before": before, "after": signature(reaction),
                     "status": "already_correct" if before == signature(reaction) else "applied",
                     "mass_balance_residual": reaction.check_mass_balance()})
    protected_definitions(model)
    return {"variant": variant, "curation_id": spec["curation_id"], "items": rows, "status": "complete"}


def model_definition(model):
    """Full SBML-representable optimization definition, independent of object identity."""
    return {"reactions": {r.id: signature(r) for r in model.reactions},
            "metabolites": {m.id: {k: getattr(m, k) for k in ("name", "formula", "charge", "compartment")}
                            for m in model.metabolites},
            "objective": {r.id: c for r, c in linear_reaction_coefficients(model).items()},
            "objective_direction": model.objective.direction}


def solver_definition(model):
    """Include actual row coefficients/bounds so a temporary edit cannot disappear in SBML."""
    model.solver.update()
    def terms(expression):
        return {getattr(term, "name", str(term)): float(value)
                for term, value in expression.as_coefficients_dict().items() if value}
    return {"variables": {v.name: (v.lb, v.ub, v.type) for v in model.variables},
            "constraints": {c.name: (c.lb, c.ub, terms(c.expression)) for c in model.constraints},
            "objective": (model.objective.direction, terms(model.objective.expression))}


def export_candidate(model, output):
    """Require a new file and verify every model field plus candidate guards on reload."""
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite {output}")
    protected_definitions(model)
    model.solver.update()
    if set(c.name for c in model.constraints) != {m.id for m in model.metabolites}:
        raise ValueError("Candidate export contains non-SBML custom constraints")
    if any((c.lb, c.ub) != (0, 0) for c in model.constraints):
        raise ValueError("Candidate export contains altered steady-state row bounds")
    before = model_definition(model)
    solver_before = solver_definition(model)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_deterministic_sbml_model(model, output)
    loaded = read_sbml_model(output)
    if model_definition(loaded) != before:
        raise ValueError(f"Candidate export/reload changed model definition: {output}")
    if solver_definition(loaded) != solver_before:
        raise ValueError(f"Candidate export/reload changed actual solver definition: {output}")
    if protected_definitions(loaded) != protected_definitions(model):
        raise ValueError("Candidate export lost protected definitions")
    return loaded


def build_candidate_file(source, output, variant="E0", spec=None):
    """Narrow final-stage build from the fixed completed reference, with provenance."""
    spec = load_spec() if spec is None else spec
    source, output = Path(source), Path(output)
    if not all(p.resolve().is_relative_to(REPO_ROOT) for p in (source, output)):
        raise ValueError("Candidate source/output must remain in this project workspace")
    input_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    implementation_sha = {str(p.relative_to(REPO_ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (Path(__file__), Path(__file__).with_name("reaction_selection.py"),
                                    Path(__file__).with_name("sbml.py"))}
    if input_sha != spec["source_sha256"]:
        raise ValueError("Energy candidate input SHA differs from the reviewed completed reference")
    model = read_sbml_model(source)
    before = model_definition(model)
    result = apply_energy_candidate(model, variant, spec)
    expected = json.loads(json.dumps(before))
    for row in result["items"]:
        expected["reactions"][row["reaction_id"]] = row["after"]
    if model_definition(model) != expected:
        raise ValueError("Candidate patch changed undeclared model fields")
    loaded = export_candidate(model, output)
    record = {"variant": variant, "source": str(source.resolve()), "source_sha256": input_sha,
              "spec_sha256": hashlib.sha256(SPEC_PATH.read_bytes()).hexdigest(),
              "implementation_sha256": implementation_sha,
              "output": str(output.resolve()), "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
              "scope": "completed-reference final candidate stage, not a full historical rebuild",
              "default_E0_definition_unchanged": variant != "E0" or model_definition(loaded) == before,
              "identity_scope": "loaded definitions; SBML serialization may differ from source bytes",
              "export_reload_definition_match": True, "non_target_fields_unchanged": True,
              "input_unchanged": hashlib.sha256(source.read_bytes()).hexdigest() == input_sha,
              "implementation_unchanged": all(hashlib.sha256((REPO_ROOT / p).read_bytes()).hexdigest() == digest
                                              for p, digest in implementation_sha.items()),
              "patch": result}
    output.with_suffix(".build.json").write_text(json.dumps(record, indent=2) + "\n")
    return record
