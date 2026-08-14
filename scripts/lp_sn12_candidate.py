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
    required = {"schema_version", "source", "chains", "states", "steps", "biomass", "blockers"}
    if curation.get("schema_version") != 1 or required - curation.keys():
        raise ContractError("invalid lipid-unlump curation manifest")
    if len(curation["chains"]) != 7 or len(curation["steps"]) != 18:
        raise ContractError("curation must define seven chains and eighteen templates")
    if curation["source"].get("model_sha256") != "0f3a6c2b151e945b3461d3fa85f04575f8e8570ba817ed2879013aec91f62415":
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
        for label1, _, _ in chains:
            for label2, _, _ in chains:
                pair = (label1, label2)
                source = made[(source_state, pair if source_state != "lpa_lp" else (label1,))]
                formula, charge = _derived_product(template, {step["template_source_id"]: source}, step["template_target_id"])
                target = _candidate_metabolite(target_state, pair, formula, charge)
                made[(target_state["id"], pair)] = target
                suffix = f"SN1__{label1}__SN2__{label2}"
                generated.append(_candidate_reaction(template, f"UL_{template.id}_{suffix}", _replace(template, {step["template_source_id"]: source}, target, step["template_target_id"])))
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


def build_candidate(source_model: Model) -> Model:
    """Build a non-activatable candidate; source input is never modified."""
    curation = _read_curation()
    candidate = source_model.copy()
    templates = _validate_source(candidate, curation)
    _build_routes(candidate, templates, curation)
    _rewrite_biomass(candidate, curation)
    _remove_generic_interfaces(candidate, curation)
    if any(reaction.check_mass_balance() for reaction in candidate.reactions if reaction.annotation.get(MARKER) == "true"):
        raise ContractError("generated reaction is not mass/charge balanced")
    return candidate


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
        "remaining_blockers": curation["blockers"],
        "generic_acyl_coa_ids_remaining": [metabolite_id for metabolite_id in curation["generic_acyl_coa_ids"] if metabolite_id in candidate.metabolites],
        "production_apply_forbidden": True,
    }


def main() -> None:
    import argparse
    import cobra

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path, help="read-only SBML source model")
    parser.add_argument("--report", type=Path, required=True, help="JSON report; SBML output is forbidden")
    arguments = parser.parse_args()
    if arguments.report.suffix != ".json":
        raise SystemExit("--report must be a JSON path; SBML output is forbidden")
    curation = _read_curation()
    actual_sha = sha256(arguments.model.read_bytes()).hexdigest()
    if actual_sha != curation["source"]["model_sha256"]:
        raise SystemExit("source model SHA drifted; refusing candidate build")
    source = cobra.io.read_sbml_model(str(arguments.model))
    payload = report(source, build_candidate(source))
    arguments.report.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
