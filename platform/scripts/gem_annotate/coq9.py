"""Version-local CoQ9 curation, after annotation/GPR assembly and before export."""

import ast
import copy
import csv
import json
import math
from fractions import Fraction
from html import unescape

from cobra import Metabolite
from cobra.core.gene import GPR

from .config import MODEL

CURATION_PATH = MODEL.curation_file("coq9_curation.json")
GENE_EVIDENCE_PATH = MODEL.curation_file("coq9_gene_evidence.tsv")
FUNCTIONAL_GPR_PATH = MODEL.curation_file("coq9_functional_gpr.json")
C5_GPR_PATH = MODEL.curation_file("coq_c5_gpr.json")
LITERATURE_PATH = C5_GPR_PATH.with_name("coq_literature_revision.json")
BIOMASS_DILUTION_PATH = C5_GPR_PATH.with_name("coq9_biomass_dilution.json")
BIOMASS = "biomass_C"
Q9, Q9H2 = "m468[C_mi]", "m471[C_mi]"
BIOMASS_DILUTION_NOTE = "coq9_growth_dilution_candidate"


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
    bio.notes[BIOMASS_DILUTION_NOTE] = json.dumps(record, sort_keys=True)
    return record


def apply_coq9_biomass_dilution(model, spec=None):
    """Curated CoQ9 growth-dilution term (user-selected, uncalibrated alpha)."""
    spec = json.loads(BIOMASS_DILUTION_PATH.read_text()) if spec is None else spec
    if (spec.get("schema_version") != 1 or spec.get("biomass_reaction") != BIOMASS
            or spec.get("metabolite") != Q9 or spec.get("reduced_metabolite") != Q9H2
            or spec.get("evidence_tier") != "assumption"):
        raise ValueError("Unexpected CoQ9 biomass dilution curation scope")
    before = model.reactions.get_by_id(BIOMASS).metabolites.get(model.metabolites.get_by_id(Q9), 0)
    record = apply_coq_biomass(model, spec["alpha_mmol_per_gDW"], spec["alpha_source"])
    return {"curation_id": spec["curation_id"],
            "status": "already_correct" if before == -spec["alpha_mmol_per_gDW"] else "applied",
            **record}


def apply_coq_literature_revision(model, enabled=False, spec=None):
    """Correct R19/R18 quinol chemistry; retain the unresolved downstream gap."""
    if not enabled:
        return {"status": "disabled"}
    spec = json.loads(LITERATURE_PATH.read_text()) if spec is None else spec
    if (spec.get("schema_version") != 1 or spec.get("enabled_by_default") is not False
            or spec.get("status") != "partial_literature_candidate_redox_gap"
            or set(spec["reactions"]) != {"R19", "R18"}
            or set(spec["new_metabolites"]) != {"coq_ddmq9h2[C_mi]", "coq_dmq9h2[C_mi]"}
            or set(spec["context_guards"]) != {"R39", "R695", "R385"}):
        raise ValueError("Unexpected CoQ literature revision scope")
    expected_gprs = {"R19": "YALI1A08781g and YALI1B03314g and YALI1B19490g",
                     "R18": "YALI1C25352g"}
    for rid, row in spec["reactions"].items():
        if (boolean_key(GPR.from_string(row["after"]["gpr"]).body)
                != boolean_key(GPR.from_string(expected_gprs[rid]).body)
                or row["after"]["bounds"] != row["before"]["bounds"]):
            raise ValueError(f"CoQ GPR logic or bounds outside scope: {rid}")

    def matches(reaction, expected):
        # COBRA collapses singleton annotations and escapes note arrows on SBML read.
        def annotations(values):
            return {k: v[0] if isinstance(v, list) and len(v) == 1 else v
                    for k, v in values.items()}

        def notes(values):
            return {k: unescape(v) if isinstance(v, str) else v
                    for k, v in values.items()}

        return (not _identity_errors(reaction, expected, {})
                and reaction.name == expected["name"]
                and annotations(reaction.annotation) == annotations(expected["annotation"])
                and notes(reaction.notes) == notes(expected["notes"]))

    for rid, guard in spec["context_guards"].items():
        if not matches(model.reactions.get_by_id(rid), guard):
            raise ValueError(f"CoQ context differs: {rid}; no edits applied")
    states = [state for state in ("before", "after") if all(
        matches(model.reactions.get_by_id(rid), row[state])
        for rid, row in spec["reactions"].items())]
    if len(states) != 1:
        raise ValueError("CoQ reactions differ or are partially modified; no edits applied")
    state = states[0]
    expected_adjacency = ({"m59[C_mi]": {"R19", "R18"}, "m61[C_mi]": {"R18", "R695"}}
                          if state == "before" else {
                              "m59[C_mi]": set(), "m61[C_mi]": {"R695"},
                              "coq_ddmq9h2[C_mi]": {"R19", "R18"},
                              "coq_dmq9h2[C_mi]": {"R18"}})
    for mid, neighbors in expected_adjacency.items():
        if mid not in model.metabolites or {
                r.id for r in model.metabolites.get_by_id(mid).reactions} != neighbors:
            raise ValueError(f"CoQ redox boundary differs: {mid}")
    note = model.notes.get("coq_literature_revision")
    if note != (None if state == "before" else spec["model_note"]):
        raise ValueError("CoQ model status note differs; no edits applied")
    for gid, identity in spec["genes"].items():
        if (gid not in model.genes or model.genes.get_by_id(gid).annotation.get("refseq")
                not in (identity["refseq"], [identity["refseq"]])):
            raise ValueError(f"CoQ protein identity differs: {gid}")

    metabolites = {met.id: met for met in model.metabolites}
    for mid, props in spec["new_metabolites"].items():
        if state == "before":
            if mid in metabolites:
                raise ValueError(f"CoQ new metabolite ID collision: {mid}")
            met = Metabolite(mid)
            for key, value in props.items():
                setattr(met, key, copy.deepcopy(value))
            metabolites[mid] = met
        elif mid not in metabolites or any(
                getattr(metabolites[mid], key) != value for key, value in props.items()):
            raise ValueError(f"CoQ quinol identity differs: {mid}")
    for rid, row in spec["reactions"].items():
        after = row["after"]
        if not set(GPR.from_string(after["gpr"]).genes) <= set(spec["genes"]):
            raise ValueError(f"Unverified CoQ gene in {rid}")
        for mid, props in after["species"].items():
            met = metabolites[mid]
            if any(getattr(met, key) != props[key] for key in ("formula", "charge", "compartment")):
                raise ValueError(f"CoQ species identity differs: {mid}")
        probe = model.reactions.get_by_id(rid).copy()
        probe.subtract_metabolites(dict(probe.metabolites))
        probe.add_metabolites({metabolites[mid].copy(): props["coefficient"]
                               for mid, props in after["species"].items()})
        if exact_residual(probe):
            raise ValueError(f"Unbalanced CoQ literature candidate: {rid}")
    if state == "after":
        return {"status": "already_correct", "scientific_status": spec["status"]}

    # Check the complete two-reaction patch before the first model mutation.
    model.add_metabolites([metabolites[mid] for mid in spec["new_metabolites"]])
    for rid, row in spec["reactions"].items():
        reaction, after = model.reactions.get_by_id(rid), row["after"]
        reaction.subtract_metabolites(dict(reaction.metabolites))
        reaction.add_metabolites({metabolites[mid]: props["coefficient"]
                                 for mid, props in after["species"].items()})
        reaction.gene_reaction_rule = after["gpr"]
        for key in ("name", "annotation", "notes"):
            setattr(reaction, key, copy.deepcopy(after[key]))
    model.notes["coq_literature_revision"] = spec["model_note"]
    return {"status": "applied", "scientific_status": spec["status"]}


def apply_coq_c5_gpr(model):
    """Opt-in, donor-coupled mitochondrial C5 hypothesis; replace the old bypass."""
    spec = json.loads(C5_GPR_PATH.read_text())
    target = "YALI1A08781g and YALI1B03314g and YALI1B19490g"
    if (spec["schema_version"] != 1 or spec["reaction_id"] != "R39"
            or spec["status"] != "provisional_cross_species_coupled_dependency"
            or spec["after"]["gpr"] != target
            or set(spec["genes"]) != set(target.split(" and "))):
        raise ValueError("Unexpected CoQ C5 hypothesis scope")
    reaction = model.reactions.get_by_id("R39")
    if (_identity_errors(reaction, spec["before"], {})
            and _identity_errors(reaction, spec["after"], {})):
        raise ValueError("CoQ C5 reaction precondition differs; no edits applied")
    if reaction.name not in (spec["name_before"], spec["name_after"]):
        raise ValueError("CoQ C5 reaction name differs")
    metabolites = {}
    for mid, expected in spec["after"]["species"].items():
        met = model.metabolites.get_by_id(mid)
        if (met.formula, met.charge, met.compartment) != (
                expected["formula"], expected["charge"], expected["compartment"]):
            raise ValueError(f"CoQ C5 species identity differs: {mid}")
        metabolites[met] = expected["coefficient"]
    for gid, identity in spec["genes"].items():
        if gid not in model.genes or model.genes.get_by_id(gid).annotation.get("refseq") not in (
                identity["refseq"], [identity["refseq"]]):
            raise ValueError(f"CoQ C5 protein identity differs: {gid}")
    for key, value in spec["notes"].items():
        if key in reaction.notes and reaction.notes[key] != value:
            raise ValueError(f"Conflicting CoQ C5 note: {key}")
    probe = reaction.copy()
    probe.subtract_metabolites(dict(probe.metabolites))
    probe.add_metabolites({met.copy(): v for met, v in metabolites.items()})
    if exact_residual(probe):
        raise ValueError("CoQ C5 candidate does not conserve atoms and charge")
    notes = {**reaction.notes, **spec["notes"]}
    changed = bool(_identity_errors(reaction, spec["after"], {})) or (
        reaction.name != spec["name_after"] or reaction.notes != notes)
    reaction.subtract_metabolites(dict(reaction.metabolites))
    reaction.add_metabolites(metabolites)
    reaction.gene_reaction_rule = target
    reaction.name, reaction.notes = spec["name_after"], notes
    return {"item": reaction.id, "hypothesis_id": spec["hypothesis_id"],
            "status": "applied" if changed else "already_correct", "gpr": target,
            "evidence_status": spec["status"]}


def apply_coq9_functional_gpr(model):
    """Opt-in substrate-access dependency; not a native catalytic-complex claim."""
    spec = json.loads(FUNCTIONAL_GPR_PATH.read_text())
    target = "YALI1E18269g and YALI1F34675g"
    if (spec["schema_version"] != 1 or spec["reaction_id"] != "R695"
            or spec["status"] != "provisional_functional_dependency"
            or spec["after_gpr"] != target
            or set(spec["genes"]) != {"YALI1E18269g", "YALI1F34675g"}):
        raise ValueError("Unexpected CoQ9 functional GPR hypothesis scope")
    reaction = model.reactions.get_by_id("R695")
    guard = copy.deepcopy(spec["guard"])
    if boolean_key(reaction.gpr.body) == boolean_key(GPR.from_string(target).body):
        guard["gpr"] = target
    if _identity_errors(reaction, guard, {}) or exact_residual(reaction):
        raise ValueError("CoQ9 functional GPR precondition differs; no edits applied")
    for gid, identity in spec["genes"].items():
        if gid not in model.genes or model.genes.get_by_id(gid).annotation.get("refseq") not in (
                identity["refseq"], [identity["refseq"]]):
            raise ValueError(f"CoQ9 functional GPR protein identity differs: {gid}")
    for key, value in spec["notes"].items():
        if key in reaction.notes and reaction.notes[key] not in (
                value, spec["notes_before"].get(key)):
            raise ValueError(f"Conflicting CoQ9 functional GPR note: {key}")
    notes = {**reaction.notes, **spec["notes"]}
    changed = reaction.gene_reaction_rule != target or reaction.notes != notes
    reaction.gene_reaction_rule = target
    reaction.notes = notes
    return {"item": reaction.id, "hypothesis_id": spec["hypothesis_id"],
            "status": "applied" if changed else "already_correct",
            "gpr": target, "evidence_status": spec["status"]}


def boolean_key(node):
    """Associative, commutative, idempotent key; never distribute a truth table."""
    if node is None:
        return None
    if isinstance(node, ast.Name):
        return ("gene", node.id)
    if not isinstance(node, ast.BoolOp):
        raise ValueError(f"Unsupported GPR node: {type(node).__name__}")
    kind = type(node.op).__name__
    children = []
    for child in node.values:
        key = boolean_key(child)
        children.extend(key[1] if key[0] == kind else [key])
    unique = tuple(sorted(set(children)))
    return unique[0] if len(unique) == 1 else (kind, unique)


def deduplicate_and(node):
    """Apply A AND A = A locally, retaining every OR and every distinct member."""
    if node is None or isinstance(node, ast.Name):
        return copy.deepcopy(node)
    if not isinstance(node, ast.BoolOp):
        raise ValueError(f"Unsupported GPR node: {type(node).__name__}")
    values = [deduplicate_and(child) for child in node.values]
    if isinstance(node.op, ast.And):
        flat = []
        for child in values:
            flat.extend(child.values if isinstance(child, ast.BoolOp)
                        and isinstance(child.op, ast.And) else [child])
        seen, values = set(), []
        for child in flat:
            key = boolean_key(child)
            if key not in seen:
                seen.add(key)
                values.append(child)
        if len(values) == 1:
            return values[0]
    return ast.BoolOp(op=copy.deepcopy(node.op), values=values)


def _identity_errors(reaction, rule, qcycle):
    actual = {m.id: {"coefficient": v, "formula": m.formula,
                     "charge": m.charge, "compartment": m.compartment}
              for m, v in reaction.metabolites.items()}
    expected = rule["species"]
    target = copy.deepcopy(expected)
    if reaction.id == "R305":
        for mid, coefficient in qcycle.items():
            target[mid]["coefficient"] = coefficient
    errors = []
    if actual != expected and actual != target:
        errors.append({"field": "species/stoichiometry", "actual": actual,
                       "expected": expected})
    if boolean_key(reaction.gpr.body) != boolean_key(GPR.from_string(rule["gpr"]).body):
        errors.append({"field": "GPR", "actual": reaction.gene_reaction_rule,
                       "expected": rule["gpr"]})
    if list(reaction.bounds) != rule["bounds"]:
        errors.append({"field": "bounds", "actual": list(reaction.bounds),
                       "expected": rule["bounds"]})
    return errors


def _field_value(reaction, field):
    if field == "name":
        return reaction.name
    if field == "ec-code":
        value = reaction.annotation.get(field, [])
        return sorted(value if isinstance(value, list) else [value])
    return reaction.notes.get(field)


def _set_field(reaction, field, value):
    if field == "name":
        reaction.name = value
    else:
        mapping = reaction.annotation if field == "ec-code" else reaction.notes
        if value is None or value == []:
            mapping.pop(field, None)
        else:
            mapping[field] = value


def exact_residual(reaction, coefficients=None):
    """Products minus reactants using the actual stored species, with exact H math."""
    residual = {}
    for met, coefficient in reaction.metabolites.items():
        if not met.formula or met.charge is None or not met.elements:
            raise ValueError(f"Missing formula/charge: {met.id}")
        coefficient = Fraction(str((coefficients or {}).get(met.id, coefficient)))
        for element, count in {**met.elements, "charge": met.charge}.items():
            residual[element] = residual.get(element, Fraction()) + coefficient * Fraction(str(count))
    return {key: str(value) for key, value in residual.items() if value}


def apply_coq9_curation(model, mode="metadata"):
    """Apply independent confirmed items, reporting conflicts without overwriting them."""
    if mode not in ("off", "metadata", "qcycle"):
        raise ValueError(f"Unknown CoQ9 mode: {mode}")
    spec = json.loads(CURATION_PATH.read_text())
    rows = []
    for rid, rule in spec["rules"].items():
        row = {"item": rid, "status": "deferred"}
        rows.append(row)
        if mode == "off":
            row["reason"] = "CoQ9 curation disabled; existing pipeline steps retained"
            continue
        if rid not in model.reactions:
            row.update(status="conflict", reason="reaction missing")
            continue
        reaction = model.reactions.get_by_id(rid)
        if rid == "R2062":
            # This Boolean identity does not depend on old chemistry or membership.
            before = reaction.gpr.body
            after = deduplicate_and(before)
            assert boolean_key(before) == boolean_key(after)
            changed = ast.dump(before) != ast.dump(after) if before is not None else False
            genes = {g.id for g in reaction.genes}
            if changed:
                reaction.gpr = GPR(ast.Expression(body=after))
            assert {g.id for g in reaction.genes} == genes
            row.update(status="applied" if changed else "already_correct",
                       before_occurrences=sum(isinstance(n, ast.Name) for n in ast.walk(before)) if before else 0,
                       after_occurrences=sum(isinstance(n, ast.Name) for n in ast.walk(after)) if after else 0,
                       unique_genes=len(genes), reason="逻辑去重完成；不证明复合体I的28基因生物学规则")
            reaction.notes["coq9_curation_scope"] = rule["scope"]
            continue
        errors = _identity_errors(reaction, rule, spec["qcycle_coefficients"])
        if errors:
            row.update(status="conflict", reason="local identity differs; current reaction preserved", errors=errors)
            continue
        fields = []
        for field, change in rule["fields"].items():
            value = _field_value(reaction, field)
            before, after = change["before"], change["after"]
            known_before = [before, *change.get("before_variants", [])]
            status = "already_correct" if value == after else "applied" if value in known_before else "conflict"
            fields.append({"field": field, "status": status, "before": value, "target": after})
        # A conflicting target annotation prevents this reaction's partial rewrite.
        if any(f["status"] == "conflict" for f in fields):
            for field in fields:
                if field["status"] == "applied":
                    field["status"] = "deferred"
            row.update(status="conflict", fields=fields, reason="target annotation differs; reaction preserved")
            continue
        for field in fields:
            if field["status"] == "applied":
                key = "coq9_previous_" + field["field"].replace("-", "_")
                previous = field["before"]
                # Plain note text survives COBRApy's XHTML quote escaping on readback.
                reaction.notes.setdefault(key, "; ".join(previous) if isinstance(previous, list) else str(previous))
                _set_field(reaction, field["field"], field["target"])
        reaction.notes["coq9_curation_scope"] = rule["scope"]
        if rid == "R19":
            # Distinguish archived pathway ECs from values actually removed
            # in this build (the offline precursor can already have no EC).
            reaction.notes["coq9_reference_pathway_ec"] = "; ".join(rule["fields"]["ec-code"]["before"])
        row.update(status="applied" if any(f["status"] == "applied" for f in fields) else "already_correct", fields=fields)

    qrow = {"item": "R305_qcycle", "status": "deferred", "reason": "explicit qcycle mode required"}
    rows.append(qrow)
    r305 = model.reactions.get_by_id("R305") if "R305" in model.reactions else None
    metadata_row = next(row for row in rows if row["item"] == "R305")
    if mode == "qcycle":
        if metadata_row["status"] == "conflict":
            qrow.update(status="conflict", reason="R305 metadata/identity conflict; proton coefficients preserved")
        else:
            target = spec["qcycle_coefficients"]
            residual = exact_residual(r305, target)
            if residual:
                qrow.update(status="conflict", reason="candidate is not balanced under stored species", residual=residual)
            else:
                before = {m.id: v for m, v in r305.metabolites.items() if m.id in target}
                r305.add_metabolites({model.metabolites.get_by_id(mid): v for mid, v in target.items()}, combine=False)
                qrow.update(status="already_correct" if before == target else "applied",
                            before=before, after=target, residual={}, scope=spec["qcycle_scope"])
    if mode != "off" and r305 is not None and metadata_row["status"] != "conflict":
        current = {m.id: v for m, v in r305.metabolites.items() if m.id in spec["qcycle_coefficients"]}
        is_candidate = current == spec["qcycle_coefficients"]
        r305.notes["coq9_proton_state"] = ("Q-cycle candidate present; " + spec["qcycle_scope"] if is_candidate else
                                          "Original audited proton stoichiometry retained; Q-cycle candidate not applied.")
        legacy_scope = ("Name/EC corrected only in metadata candidate. Original proton imbalance retained here; "
                        "separate Q-cycle candidate changes stoichiometry. Current GPR includes a cytochrome-c "
                        "carrier dependency and is not certified as a native subunit-only rule.")
        if r305.notes.get("curation_20260909_scope") == legacy_scope:
            r305.notes["coq9_previous_scope"] = r305.notes.pop("curation_20260909_scope")
    for item in spec["deferred"]:
        if item["id"] == "CIII-CHEM":
            continue
        rows.append({**item, "reference_status": item["status"], "status": "deferred"})
    r573 = model.reactions.get_by_id("R573") if "R573" in model.reactions else None
    rows.append({"item": "R573_retained", "status": "already_correct" if r573 and r573.bounds == (0, 0) else "conflict",
                 "actual_bounds": list(r573.bounds) if r573 else None,
                 "reason": "Observed only; already closed in reference, not a new repair. Different overlays are preserved."})
    return {"mode": mode, "reference_sha256": spec["reference_sha256"],
            "requested_mode_complete": not any(row["status"] == "conflict" for row in rows), "items": rows}


def gene_evidence(model):
    """Retain supplied evidence; derive representation from the actual output GPRs."""
    with GENE_EVIDENCE_PATH.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        gid = row["yali1_gene"]
        reactions = sorted(r.id for r in model.genes.get_by_id(gid).reactions) if gid in model.genes else []
        row["actual_reactions"] = ";".join(reactions)
        row["association_changed_from_reference"] = set(reactions) != set(filter(None, row["reference_reactions"].split(";")))
        row["representation_status"] = "GPR present; biological validity not established by association" if reactions else "not represented / not testable"
    return rows
