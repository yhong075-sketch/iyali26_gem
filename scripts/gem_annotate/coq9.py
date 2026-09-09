"""Version-local CoQ9 curation, after annotation/GPR assembly and before export."""

import ast
import copy
import csv
import json
from fractions import Fraction

from cobra.core.gene import GPR

from .config import REPO_ROOT

CURATION_PATH = REPO_ROOT / "data" / "coq9_curation.json"
GENE_EVIDENCE_PATH = REPO_ROOT / "data" / "coq9_gene_evidence.tsv"


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
