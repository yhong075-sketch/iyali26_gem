#!/usr/bin/env python3
"""Build the read-only, strict sn-position lipid-core candidate in memory."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from math import fsum, isclose
from pathlib import Path
from time import perf_counter
from typing import Any
import json
import re

from cobra import Metabolite, Model, Reaction


REPOSITORY = Path(__file__).resolve().parents[1]
CURATION_PATH = REPOSITORY / "data" / "lipid_unlump_sn_core_curation.json"
MARKER = "lipid_unlump_sn_core"
CARDIOLIPIN_AUDIT_SHA256 = "e16db37d6b19731fa8f608e177c10b8d646343e496b18ad609a2da29f35543d1"
_FORMULA_ORDER = ("C", "H", "N", "O", "P", "S")


class ContractError(ValueError):
    """The source snapshot or curation is not the one that was reviewed."""


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def _formula(elements: Counter[str]) -> str:
    if any(amount < 0 for amount in elements.values()):
        raise ContractError(f"negative formula coefficient: {dict(elements)}")
    order = [element for element in _FORMULA_ORDER if elements[element]]
    order += sorted(element for element in elements if element not in _FORMULA_ORDER and elements[element])
    return "".join(element + (str(elements[element]) if elements[element] != 1 else "") for element in order)


def _elements(metabolite: Metabolite) -> Counter[str]:
    if metabolite.elements is None:
        raise ContractError(f"formula required for {metabolite.id}")
    return Counter(metabolite.elements)


def _label(metabolite: Metabolite) -> str:
    label = re.sub(r"_c\d+h\d+.*$", "", metabolite.name.lower()).replace("-coa", "")
    label = re.sub(r"[^a-z0-9]+", "_", label).strip("_")
    if not label:
        raise ContractError(f"cannot derive chain label from {metabolite.id}")
    return label


def _reaction(model: Model, reaction_id: str) -> Reaction:
    try:
        return model.reactions.get_by_id(reaction_id)
    except KeyError as error:
        raise ContractError(f"missing required reaction {reaction_id}") from error


def _metabolite(model: Model, metabolite_id: str) -> Metabolite:
    try:
        return model.metabolites.get_by_id(metabolite_id)
    except KeyError as error:
        raise ContractError(f"missing required metabolite {metabolite_id}") from error


def _read_curation(path: Path = CURATION_PATH) -> dict[str, Any]:
    try:
        curation = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ContractError(f"cannot read curation {path}: {error}") from error
    required = {"schema_version", "source", "R989_gpr_curation", "R1521_neutral_formula_completion", "chains", "states", "steps", "biomass", "cardiolipin_four_chain_audit", "blockers"}
    if curation.get("schema_version") != 3 or required - curation.keys():
        raise ContractError("invalid lipid-unlump curation manifest")
    if len(curation["chains"]) != 7 or len(curation["steps"]) != 18:
        raise ContractError("curation must define seven chains and eighteen templates")
    if curation["source"].get("model_sha256") != "bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee":
        raise ContractError("curation source SHA is not the approved snapshot")
    return curation


def source_fingerprint(model: Model) -> str:
    """Hash every source property used by this candidate; drift fails closed."""
    metabolites = [
        (met.id, met.name, met.formula, met.charge, met.compartment, _json(met.annotation))
        for met in sorted(model.metabolites, key=lambda item: item.id)
    ]
    reactions = [
        (reaction.id, reaction.name, tuple(reaction.bounds), reaction.gene_reaction_rule,
         _json(reaction.annotation), tuple(sorted((met.id, float(value)) for met, value in reaction.metabolites.items())))
        for reaction in sorted(model.reactions, key=lambda item: item.id)
    ]
    objective = _objective(model)
    return sha256(_json((metabolites, reactions, objective)).encode()).hexdigest()


def _objective(model: Model) -> tuple[tuple[str, float], ...]:
    return tuple(sorted((reaction.id, float(reaction.objective_coefficient)) for reaction in model.reactions if reaction.objective_coefficient))


def _validate_source(model: Model, curation: dict[str, Any]) -> dict[str, Reaction]:
    if any(reaction.annotation.get(MARKER) == "true" for reaction in model.reactions):
        raise ContractError("already-marked candidate inputs are refused; rebuild from the frozen source")
    actual = source_fingerprint(model)
    expected = curation["source"].get("model_fingerprint")
    if not expected or actual != expected:
        raise ContractError(f"source fingerprint drifted: expected {expected}, got {actual}")
    if _objective(model) != (("biomass_C", 1.0),):
        raise ContractError("source objective must be exactly biomass_C")
    templates = {step["reaction_id"]: _reaction(model, step["reaction_id"]) for step in curation["steps"]}
    for step in curation["steps"]:
        reaction = templates[step["reaction_id"]]
        if tuple(reaction.bounds) != tuple(step["bounds"]) or reaction.gene_reaction_rule != step["gpr"]:
            raise ContractError(f"template contract drifted: {reaction.id}")
        actual_stoichiometry = {met.id: coefficient for met, coefficient in reaction.metabolites.items()}
        if actual_stoichiometry != step["stoichiometry"]:
            raise ContractError(f"template stoichiometry drifted: {reaction.id}")
    return templates


def _candidate_metabolite(state: dict[str, Any], chains: tuple[str, ...], formula: str, charge: int) -> Metabolite:
    suffix = "__".join(f"sn{index + 1}__{chain}" for index, chain in enumerate(chains))
    compartment = state["compartment"]
    metabolite = Metabolite(
        f"ul_{state['id']}__{suffix}[{compartment}]",
        name=f"{state['name']} ({', '.join(chains)})",
        formula=formula,
        charge=charge,
        compartment=compartment,
    )
    metabolite.annotation = {MARKER: "true", "state": state["id"], "sn_tuple": "/".join(chains)}
    return metabolite


def _candidate_reaction(template: Reaction, reaction_id: str, stoichiometry: dict[Metabolite, float]) -> Reaction:
    reaction = Reaction(reaction_id, name=f"{template.name} (strict sn state)")
    reaction.bounds = template.bounds
    reaction.gene_reaction_rule = template.gene_reaction_rule
    reaction.annotation = {**template.annotation, MARKER: "true", "template_reaction": template.id}
    reaction.add_metabolites(stoichiometry)
    return reaction


def _derived_product(template: Reaction, replace: dict[str, Metabolite], target_id: str) -> tuple[str, int]:
    elements: Counter[str] = Counter()
    charge = 0
    target_coefficient = next((coefficient for metabolite, coefficient in template.metabolites.items() if metabolite.id == target_id), None)
    if target_coefficient not in {-1.0, 1.0}:
        raise ContractError(f"unsupported target coefficient in {template.id}: {target_id}")
    for old, coefficient in template.metabolites.items():
        if old.id == target_id:
            continue
        metabolite = replace.get(old.id, old)
        for element, amount in _elements(metabolite).items():
            elements[element] -= int(coefficient * amount / target_coefficient)
        if metabolite.charge is None:
            raise ContractError(f"charge required for {metabolite.id}")
        charge -= int(coefficient * metabolite.charge / target_coefficient)
    return _formula(elements), charge


def _replace(template: Reaction, replacements: dict[str, Metabolite], target: Metabolite, target_id: str) -> dict[Metabolite, float]:
    stoichiometry: dict[Metabolite, float] = {}
    for old, coefficient in template.metabolites.items():
        metabolite = target if old.id == target_id else replacements.get(old.id, old)
        stoichiometry[metabolite] = coefficient
    return stoichiometry


def _states_by_id(curation: dict[str, Any]) -> dict[str, dict[str, Any]]:
    states = {state["id"]: state for state in curation["states"]}
    if len(states) != 18:
        raise ContractError("curation must define 18 strict lipid states")
    return states


def _chains(model: Model, curation: dict[str, Any]) -> list[tuple[str, Metabolite, float]]:
    result = []
    for row in curation["chains"]:
        metabolite = _metabolite(model, row["lp_acyl_coa_id"])
        if _label(metabolite) != row["id"] or metabolite.formula != row["legacy_tuple"]["formula"] or metabolite.charge != row["legacy_tuple"]["charge"]:
            raise ContractError(f"acyl-CoA identity drifted: {row['id']}")
        result.append((row["id"], metabolite, row["weight"]))
    if not isclose(fsum(weight for _, _, weight in result), 1.0, abs_tol=1e-15):
        raise ContractError("acyl-chain weights do not normalize to one")
    return result


def cardiolipin_audit(model: Model) -> dict[str, Any]:
    """Audit the four-position CL envelope without adding a state or reaction."""
    curation = _read_curation()
    contract = curation["cardiolipin_four_chain_audit"]
    contract_sha256 = sha256(_json(contract).encode()).hexdigest()
    if contract_sha256 != CARDIOLIPIN_AUDIT_SHA256:
        raise ContractError("cardiolipin audit contract drifted")
    genes = contract["gene_identity_ledger"]
    gene_ids = {row["systematic_id"] for row in genes}
    template_gprs = {row["gpr"] for row in contract["template_reactions"] if row["gpr"]}
    if len(gene_ids) != len(genes) or not template_gprs <= gene_ids or any(not row["evidence_status"] for row in genes):
        raise ContractError("cardiolipin gene identity ledger drifted")
    sources = contract["evidence_sources"]
    coverage = contract["evidence_audit_coverage"]
    conflict = contract["R1742_annotation_conflict"]
    claim_ids = [claim_id for row in sources for claim_id in row["claim_ids"]] + [conflict["claim_id"]]
    verdicts = Counter(row["reviewer_verdict"] for row in sources for _ in row["claim_ids"])
    verdicts[conflict["reviewer_verdict"]] += 1
    if (
        len(claim_ids) != len(set(claim_ids))
        or len(claim_ids) != coverage["total_claims"]
        or coverage["audited"] != coverage["total_claims"]
        or coverage["supported"] != verdicts["supported"]
        or coverage["unresolved"] != verdicts["partially_supported"]
        or coverage["contradicted"] != verdicts["contradicted"]
        or coverage["supported"] + coverage["unresolved"] + coverage["contradicted"] != coverage["total_claims"]
        or coverage["unchecked"] != 0
        or set(verdicts) - {"supported", "partially_supported", "contradicted"}
        or conflict["audit_status"] != "audited"
        or conflict["decision"] not in {"use_as_constraint", "defer", "exclude"}
        or not conflict["conditions"]
        or not set(conflict["source_ids"]) <= {row["source_id"] for row in sources}
        or any(
            row["audit_status"] != "audited"
            or row["decision"] not in {"use_as_constraint", "defer", "exclude"}
            or not row["conditions"]
            for row in sources
        )
    ):
        raise ContractError("cardiolipin evidence ledger coverage drifted")
    positions = contract["position_convention"]
    for row in contract["template_reactions"]:
        reaction = _reaction(model, row["reaction_id"])
        actual = {
            "bounds": list(reaction.bounds),
            "gpr": reaction.gene_reaction_rule,
            "compartments": sorted(reaction.compartments),
            "stoichiometry": {met.id: coefficient for met, coefficient in reaction.metabolites.items()},
        }
        if "notes" in row:
            actual["notes"] = reaction.notes
        expected = {key: row[key] for key in actual}
        if actual != expected:
            raise ContractError(f"cardiolipin template contract drifted: {reaction.id}")

    chains = [row["id"] for row in curation["chains"]]
    half_labels = [f"sn1={first}|sn2={second}" for first in chains for second in chains]
    ids_and_symmetry = [
        (f"ul_cl_mm__proS_sn1__{a}__proS_sn2__{b}__proR_sn1__{c}__proR_sn2__{d}[C_mm]", (a, b) == (c, d))
        for a in chains for b in chains for c in chains for d in chains
    ]
    digest = lambda values: sha256(("\n".join(sorted(values)) + "\n").encode()).hexdigest()  # noqa: E731
    all_ids = [identifier for identifier, _ in ids_and_symmetry]
    same_ids = [identifier for identifier, same in ids_and_symmetry if same]
    different_ids = [identifier for identifier, same in ids_and_symmetry if not same]
    envelope = {
        "status": "not_adopted",
        "all_state_count": len(all_ids),
        "all_ids_sha256": digest(all_ids),
        "same_half_state_count": len(same_ids),
        "same_half_ids_sha256": digest(same_ids),
        "different_half_state_count": len(different_ids),
        "different_half_ids_sha256": digest(different_ids),
    }
    state_space = contract["state_space"]
    if len(half_labels) != state_space["ordered_half_pair_count"] or digest(half_labels) != state_space["ordered_half_labels_sha256"] or envelope != state_space["not_adopted_modeling_envelope"]:
        raise ContractError("cardiolipin structural-envelope contract drifted")
    forbidden = state_space["forbidden_lossy_half_swap_fold"]
    if forbidden["state_count"] != len(half_labels) * (len(half_labels) + 1) // 2:
        raise ContractError("lossy cardiolipin half-swap count drifted")

    positioned = sorted(met.id for met in model.metabolites if met.id.startswith("ul_cl_mm__"))
    positioned_set = set(positioned)
    position_reactions = sorted(
        reaction.id for reaction in model.reactions
        if any(met.id in positioned_set for met in reaction.metabolites)
    )
    current = {
        "currently_instantiated_network_position_state_count": 0,
        "evidence_supported_position_state_count": 0,
        "included_position_metabolite_count": len(positioned),
        "included_position_reaction_count": len(position_reactions),
        "reachability_scope": "zero means no position-resolved cardiolipin state is instantiated in the current candidate; it does not prove all 2,401 envelope states biologically unreachable",
    }
    if current != state_space["current_strict_model"]:
        raise ContractError("current strict cardiolipin state must remain empty")

    biomass_contract = contract["biomass_contract"]
    biomass = _reaction(model, biomass_contract["reaction_id"])
    generic_cl = biomass_contract["xLIPID_generic_cardiolipin_metabolite_id"]
    if any(met.id == generic_cl or met.id.startswith("ul_cl_mm__") for met in biomass.metabolites):
        raise ContractError("cardiolipin biomass activation is forbidden")
    xlipid = _reaction(model, "xLIPID")
    if xlipid.metabolites.get(_metabolite(model, generic_cl)) != biomass_contract["xLIPID_generic_cardiolipin_coefficient"]:
        raise ContractError("generic xLIPID cardiolipin coefficient drifted")
    evidence = contract["network_evidence"]
    activation = contract["activation"]

    witness = [
        f"ul_cl_mm__proS_sn1__{chains[a]}__proS_sn2__{chains[b]}__proR_sn1__{chains[c]}__proR_sn2__{chains[d]}[C_mm]"
        for a, b, c, d in ((0, 1, 2, 3), (2, 3, 0, 1))
    ]
    return {
        "contract_sha256": contract_sha256,
        "mode": contract["mode"],
        "activation": activation,
        "gates": contract["gates"],
        "validated_template_reaction_ids": [row["reaction_id"] for row in contract["template_reactions"]],
        "position_convention": positions,
        "state_space": {"digest_convention": state_space["digest_convention"], "ordered_half_pair_count": len(half_labels), "ordered_half_labels_sha256": digest(half_labels), "not_adopted_modeling_envelope": envelope, "forbidden_lossy_half_swap_fold": forbidden, "current_strict_model": current},
        "half_swap_witness": {"original_id": witness[0], "swapped_id": witness[1], "ids_are_distinct": witness[0] != witness[1]},
        "biomass_contract": {**biomass_contract, "positioned_biomass_term_count": 0},
        "network_evidence": evidence,
        "gene_identity_ledger": genes,
        "R197_gpr_decision_contract": contract["R197_gpr_decision_contract"],
        "R1742_annotation_conflict": conflict,
        "evidence_sources": sources,
        "evidence_audit_coverage": coverage,
        "curation_path": str(CURATION_PATH.relative_to(REPOSITORY)),
        "remaining_blockers": contract["blockers"],
    }


def _build_routes(candidate: Model, templates: dict[str, Reaction], curation: dict[str, Any]) -> list[Reaction]:
    states, chains = _states_by_id(curation), _chains(candidate, curation)
    generated: list[Reaction] = []
    made: dict[tuple[str, tuple[str, ...]], Metabolite] = {}
    first = templates["R353"]
    lpa_state = states["lpa_lp"]
    g3p, coa = _metabolite(candidate, "m576[C_lp]"), _metabolite(candidate, "m573[C_lp]")
    for label, acyl_coa, _ in chains:
        formula, charge = _derived_product(first, {"m570[C_lp]": acyl_coa, "m576[C_lp]": g3p, "m573[C_lp]": coa}, "m577[C_lp]")
        lpa = _candidate_metabolite(lpa_state, (label,), formula, charge)
        made[("lpa_lp", (label,))] = lpa
        generated.append(_candidate_reaction(first, f"UL_R353_SN1__{label}", _replace(first, {"m570[C_lp]": acyl_coa}, lpa, "m577[C_lp]")))
    second, pa_state = templates["R1846"], states["pa_lp"]
    for label1, _, _ in chains:
        lpa = made[("lpa_lp", (label1,))]
        for label2, acyl_coa, _ in chains:
            pair = (label1, label2)
            formula, charge = _derived_product(second, {"m577[C_lp]": lpa, "m570[C_lp]": acyl_coa}, "m1705[C_lp]")
            pa = _candidate_metabolite(pa_state, pair, formula, charge)
            made[("pa_lp", pair)] = pa
            suffix = f"SN1__{label1}__SN2__{label2}"
            generated.append(_candidate_reaction(second, f"UL_R1846_{suffix}", _replace(second, {"m577[C_lp]": lpa, "m570[C_lp]": acyl_coa}, pa, "m1705[C_lp]")))
    for step in curation["steps"]:
        if step["reaction_id"] in {"R353", "R1846", "R1771"}:
            continue
        template, source_state, target_state = templates[step["reaction_id"]], step["source_state"], states[step["target_state"]]
        template_source_id, template_target_id = step["template_source_id"], step["template_target_id"]
        if step.get("construction_direction") == "reverse":
            source_state, target_state = target_state["id"], states[source_state]
            template_source_id, template_target_id = template_target_id, template_source_id
        for label1, _, _ in chains:
            for label2, _, _ in chains:
                pair = (label1, label2)
                source = made[(source_state, pair if source_state != "lpa_lp" else (label1,))]
                formula, charge = _derived_product(template, {template_source_id: source}, template_target_id)
                target = _candidate_metabolite(target_state, pair, formula, charge)
                made[(target_state["id"], pair)] = target
                suffix = f"SN1__{label1}__SN2__{label2}"
                reaction = _candidate_reaction(template, f"UL_{template.id}_{suffix}", _replace(template, {template_source_id: source}, target, template_target_id))
                if "reaction_classification" in step:
                    reaction.annotation["reaction_classification"] = step["reaction_classification"]
                generated.append(reaction)
    tag_template, tag_state = templates["R1771"], states["tag_lp"]
    for label1, _, _ in chains:
        for label2, _, _ in chains:
            dag = made[("dag_lp", (label1, label2))]
            for label3, acyl_coa, _ in chains:
                triplet = (label1, label2, label3)
                formula, charge = _derived_product(tag_template, {"m1642[C_lp]": dag, "m570[C_lp]": acyl_coa}, "m1640[C_lp]")
                tag = _candidate_metabolite(tag_state, triplet, formula, charge)
                made[("tag_lp", triplet)] = tag
                suffix = f"SN1__{label1}__SN2__{label2}__SN3__{label3}"
                generated.append(_candidate_reaction(tag_template, f"UL_R1771_{suffix}", _replace(tag_template, {"m1642[C_lp]": dag, "m570[C_lp]": acyl_coa}, tag, "m1640[C_lp]")))
    candidate.add_reactions(generated)
    if len(generated) != 1134 or len({met.id for reaction in generated for met in reaction.metabolites if met.annotation.get(MARKER) == "true"}) != 1134:
        raise ContractError("strict core must contain exactly 1,134 routes and states")
    return generated


def _rewrite_biomass(candidate: Model, curation: dict[str, Any]) -> None:
    biomass = _reaction(candidate, "biomass_C")
    chains = [(row["id"], row["weight"]) for row in curation["chains"]]
    for row in curation["biomass"]["terms"]:
        generic = _metabolite(candidate, row["generic_metabolite_id"])
        original = biomass.metabolites.get(generic)
        if original != -row["coefficient"]:
            raise ContractError(f"biomass coefficient drifted: {generic.id}")
        biomass.add_metabolites({generic: -original})
        coefficients = []
        if row["arity"] == 2:
            for first, p1 in chains:
                for second, p2 in chains:
                    metabolite = _metabolite(candidate, f"ul_{row['state']}__sn1__{first}__sn2__{second}[{row['compartment']}]")
                    coefficient = row["coefficient"] * p1 * p2
                    biomass.add_metabolites({metabolite: -coefficient})
                    coefficients.append(coefficient)
        else:
            for first, p1 in chains:
                for second, p2 in chains:
                    for third, p3 in chains:
                        metabolite = _metabolite(candidate, f"ul_{row['state']}__sn1__{first}__sn2__{second}__sn3__{third}[{row['compartment']}]")
                        coefficient = row["coefficient"] * p1 * p2 * p3
                        biomass.add_metabolites({metabolite: -coefficient})
                        coefficients.append(coefficient)
        if not isclose(fsum(coefficients), row["coefficient"], rel_tol=0.0, abs_tol=1e-18):
            raise ContractError(f"biomass expansion does not restore {generic.id}")


def _remove_generic_interfaces(candidate: Model, curation: dict[str, Any]) -> None:
    remove_ids = [step["reaction_id"] for step in curation["steps"]] + curation["remove_reactions"]
    candidate.remove_reactions([_reaction(candidate, reaction_id) for reaction_id in remove_ids])
    for reaction_id in curation["close_reactions"]:
        _reaction(candidate, reaction_id).bounds = (0.0, 0.0)
    for reaction_id in ("newBiom", "xBIOMASS", "xLIPID"):
        _reaction(candidate, reaction_id).bounds = (0.0, 0.0)
    generic_ids = tuple(curation["generic_acyl_coa_ids"])
    if any(_metabolite(candidate, metabolite_id).reactions for metabolite_id in generic_ids):
        raise ContractError("generic acyl-CoA has a remaining reaction before deletion")
    candidate.remove_metabolites([_metabolite(candidate, metabolite_id) for metabolite_id in generic_ids], destructive=False)
    candidate.objective = "biomass_C"


def _apply_r989_candidate_gpr(candidate: Model, curation: dict[str, Any]) -> None:
    contract = curation["R989_gpr_curation"]
    reaction = _reaction(candidate, contract["reaction_id"])
    gene_id = contract["approved_candidate_gpr"]
    assignment = contract["model_assignment"]
    notes = contract["notes"]
    sources = contract["evidence_sources"]
    coverage = contract["evidence_audit_coverage"]
    claims = [claim_id for source in sources for claim_id in source["claim_ids"]]
    verdicts = Counter(source["reviewer_verdict"] for source in sources for _ in source["claim_ids"])
    if (
        reaction.gene_reaction_rule != contract["expected_source_gpr"]
        or sorted(reaction.compartments) != contract["expected_compartments"]
        or contract["gene_identity"]["systematic_id"] != gene_id
        or contract["target_sequence"] != {
            "accession": "XP_502280.2",
            "length_aa": 235,
            "sha256_convention": "uppercase amino-acid sequence only; no FASTA header or whitespace",
            "sha256": "d8cea1ec65af975fe803163cfaa05eb6fb3961dfce8df5546a4f14e1c0e007f9",
        }
        or contract["alphafold"]["scope"] != "fold_and_family_compatibility_only"
        or contract["alphafold"]["input_length_aa"] != 212
        or contract["localization_counterevidence"]["prediction"] != "Other"
        or assignment != {
            "evidence_category": "model/GPR assignment only",
            "status": "candidate_only",
            "human_model_assignment_approved": True,
            "wet_lab_review_required": True,
            "wet_lab_review_status": "pending",
            "mitochondrial_localization_status": "unverified",
            "production_claim_allowed": False,
        }
        or set(notes) != {"curated_gpr_correction", "gpr_evidence_status", "gpr_evidence_limit"}
        or any(not value for value in notes.values())
        or any(
            source["audit_status"] != "audited"
            or source["decision"] not in {"use_as_constraint", "defer", "exclude"}
            or not source["conditions"]
            for source in sources
        )
        or len(claims) != len(set(claims))
        or coverage != {
            "total_claims": len(claims),
            "audited": len(claims),
            "supported": verdicts["supported"],
            "unresolved": verdicts["partially_supported"],
            "contradicted": verdicts["contradicted"],
            "unchecked": 0,
        }
        or _reaction(candidate, "R64").gene_reaction_rule != gene_id
    ):
        raise ContractError("R989 candidate GPR evidence contract drifted")
    reaction.gene_reaction_rule = gene_id
    reaction.notes = {**reaction.notes, **notes}


def _apply_r1521_neutral_formula(candidate: Model, curation: dict[str, Any]) -> None:
    contract = curation["R1521_neutral_formula_completion"]
    metabolite = _metabolite(candidate, contract["metabolite_id"])
    same_identity = contract["same_identity_source_copy"]
    source_copy = _metabolite(candidate, same_identity["metabolite_id"])
    expected = contract["expected_source_tuple"]
    target = contract["candidate_tuple"]
    decision = contract["decision"]
    if (
        contract["reaction_id"] != "R1521"
        or contract["metabolite_id"] != "m1546[C_pe]"
        or metabolite.name != contract["expected_source_name"]
        or {"formula": metabolite.formula, "charge": metabolite.charge} != expected
        or {"formula": source_copy.formula, "charge": source_copy.charge} != same_identity["expected_tuple"]
        or target != {"formula": "C47H86N7O18P3S", "charge": 0}
        or contract["scope"] != "formula_completion_only_in_legacy_neutral_convention"
        or contract["biochemical_ph_tuple"] != {
            "chebi": "CHEBI:76378", "formula": "C47H82N7O18P3S", "charge": -4,
            "status": "blocked_not_applied",
        }
        or decision != {
            "human_candidate_approved": True,
            "change_formula": True,
            "change_charge": False,
            "change_name_or_annotation": False,
            "change_reaction_stoichiometry_bounds_or_gpr": False,
            "production_apply_allowed": False,
        }
        or any(source["audit_status"] != "audited" for source in contract["evidence_sources"])
    ):
        raise ContractError("R1521 neutral formula contract drifted")
    metabolite.formula = target["formula"]
    for reaction_id, expected_balance in contract["incident_reaction_balance_after"].items():
        if _reaction(candidate, reaction_id).check_mass_balance() != expected_balance:
            raise ContractError(f"R1521 formula completion balance drifted: {reaction_id}")


def build_candidate(source_model: Model) -> Model:
    """Build a non-activatable candidate; source input is never modified."""
    curation = _read_curation()
    candidate = source_model.copy()
    templates = _validate_source(candidate, curation)
    _apply_r1521_neutral_formula(candidate, curation)
    candidate._compartments["C_em"] = "ER membrane"
    _build_routes(candidate, templates, curation)
    _rewrite_biomass(candidate, curation)
    _remove_generic_interfaces(candidate, curation)
    _apply_r989_candidate_gpr(candidate, curation)
    if any(reaction.check_mass_balance() for reaction in candidate.reactions if reaction.annotation.get(MARKER) == "true"):
        raise ContractError("generated reaction is not mass/charge balanced")
    return candidate


def write_candidate_sbml(candidate: Model, path: Path) -> None:
    """Write stable SBML despite COBRApy storing group members in sets."""
    import cobra

    original_members = [(group, group._members) for group in candidate.groups]
    try:
        for group, members in original_members:
            # ponytail: remove when COBRApy writes Group members in stable ID order.
            group._members = tuple(sorted(members, key=lambda member: (type(member).__name__, member.id)))
        cobra.io.write_sbml_model(candidate, str(path))
    finally:
        for group, members in original_members:
            group._members = members


def _timed_objective(model: Model, pfba: bool) -> tuple[float, float]:
    start = perf_counter()
    if pfba:
        from cobra.flux_analysis import pfba as run_pfba
        value = run_pfba(model).objective_value
    else:
        value = model.slim_optimize(error_value=None)
    return value, perf_counter() - start


def report(source_model: Model, candidate: Model) -> dict[str, Any]:
    curation = _read_curation()
    generated = [reaction for reaction in candidate.reactions if reaction.annotation.get(MARKER) == "true"]
    source_fba, _ = _timed_objective(source_model, False)
    candidate_fba, _ = _timed_objective(candidate, False)
    baseline = [(_timed_objective(source_model, False)[1], _timed_objective(source_model, True)[1]) for _ in range(5)]
    timings = [(_timed_objective(candidate, False)[1], _timed_objective(candidate, True)[1]) for _ in range(5)]
    runtime = {"source_fba_median": sorted(value[0] for value in baseline)[2], "source_pfba_median": sorted(value[1] for value in baseline)[2], "candidate_fba_median": sorted(value[0] for value in timings)[2], "candidate_pfba_median": sorted(value[1] for value in timings)[2]}
    return {
        "artifact_type": "strict_sn_lipid_core_candidate_report",
        "activation_state": "blocked",
        "ready_for_activation": False,
        "source": {**curation["source"], "observed_fingerprint": source_fingerprint(source_model)},
        "architecture": {"generated_metabolites": 1134, "generated_reactions": len(generated), "sn1_lpa": 7, "ordered_diacyl_states": 16 * 49, "ordered_tag": 343},
        "fba_probe": {"source_objective": source_fba, "candidate_objective": candidate_fba, "growth_ratio": candidate_fba / source_fba if source_fba else None},
        "runtime_seconds": runtime,
        "performance_review_required": any(runtime[f"candidate_{kind}_median"] > 2 * runtime[f"source_{kind}_median"] for kind in ("fba", "pfba")),
        "cardiolipin_four_chain_audit": cardiolipin_audit(candidate),
        "R989_candidate_gpr": curation["R989_gpr_curation"],
        "R1521_neutral_formula_completion": curation["R1521_neutral_formula_completion"],
        "remaining_blockers": curation["blockers"],
        "generic_acyl_coa_ids_remaining": [metabolite_id for metabolite_id in curation["generic_acyl_coa_ids"] if metabolite_id in candidate.metabolites],
        "production_apply_forbidden": True,
    }


def main() -> None:
    import argparse
    import cobra

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path, help="read-only SBML source model")
    parser.add_argument("--report", type=Path, required=True, help="JSON report")
    parser.add_argument("--candidate-sbml", type=Path, help="write the non-activatable candidate SBML")
    arguments = parser.parse_args()
    if arguments.report.suffix != ".json":
        raise SystemExit("--report must be a JSON path")
    if arguments.candidate_sbml and arguments.candidate_sbml.suffix not in {".sbml", ".xml"}:
        raise SystemExit("--candidate-sbml must be an SBML .xml or .sbml path")
    if arguments.candidate_sbml and arguments.candidate_sbml.resolve() == arguments.model.resolve():
        raise SystemExit("--candidate-sbml cannot overwrite the frozen source model")
    curation = _read_curation()
    actual_sha = sha256(arguments.model.read_bytes()).hexdigest()
    if actual_sha != curation["source"]["model_sha256"]:
        raise SystemExit("source model SHA drifted; refusing candidate build")
    source = cobra.io.read_sbml_model(str(arguments.model))
    candidate = build_candidate(source)
    payload = report(source, candidate)
    if arguments.candidate_sbml:
        write_candidate_sbml(candidate, arguments.candidate_sbml)
    arguments.report.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
