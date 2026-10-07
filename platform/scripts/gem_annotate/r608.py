"""Opt-in R608 proposal route. No N/C production curation is bundled or selected."""

from copy import deepcopy
from math import isfinite
import re

from .execution import guarded_execution

METABOLITES = frozenset(
    (
        "m311[C_cy]",
        "m851[C_mi]",
        "m1097[C_ex]",
        "m850[C_cy]",
        "m1101[C_ex]",
        "m1102[C_va]",
    )
)
REACTIONS = frozenset(("R607", "R608"))
PROTON = "m10[C_cy]"
FIELDS = frozenset(("formula", "charge", "name", "annotation"))


def _identity(met):
    return {field: deepcopy(getattr(met, field)) for field in FIELDS}


def _stoichiometry(reaction):
    return {met.id: value for met, value in reaction.metabolites.items()}


@guarded_execution
def apply_r608_candidate(
    model, curation, *, source_model_sha256, no_solve=True, allow_network=False
):
    """Preflight the entire proposal; return a changed copy, never mutate input.

    The source identity binds the proposal to its parent. Complete after-state
    recognition permits an in-memory repeat; mixed before/after states fail.
    This checks engineering and local chemistry, not biological acceptance.
    """
    required = {
        "selected_representation",
        "default_enabled",
        "source_model_sha256",
        "metabolites",
        "reactions",
    }
    if not isinstance(curation, dict) or set(curation) != required:
        raise ValueError("R608 requires a complete specification with no extra fields")
    if curation["selected_representation"] not in {"N", "C"}:
        raise ValueError("R608 representation is not selected")
    if curation["default_enabled"] is not False:
        raise ValueError("R608 must remain disabled by default")
    if (
        not re.fullmatch(r"[0-9a-f]{64}", source_model_sha256 or "")
        or curation["source_model_sha256"] != source_model_sha256
    ):
        raise ValueError("R608 source model does not match the proposal")
    if (
        not isinstance(curation["metabolites"], dict)
        or set(curation["metabolites"]) != METABOLITES
    ):
        raise ValueError("R608 requires exactly the six scoped metabolites")
    if (
        not isinstance(curation["reactions"], dict)
        or set(curation["reactions"]) != REACTIONS
    ):
        raise ValueError("R608 requires exactly R607 and R608")

    states = []
    for met_id, record in curation["metabolites"].items():
        if not isinstance(record, dict) or set(record) != {"before", "after"}:
            raise ValueError(f"incomplete before/after: {met_id}")
        for identity in record.values():
            if not isinstance(identity, dict) or set(identity) != FIELDS:
                raise ValueError(f"incomplete chemical identity: {met_id}")
            if (
                not isinstance(identity["formula"], str)
                or not re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", identity["formula"])
                or type(identity["charge"]) is not int
                or not isinstance(identity["name"], str)
                or not identity["name"]
                or not isinstance(identity["annotation"], dict)
            ):
                raise ValueError(f"invalid chemical identity: {met_id}")
            for key, value in identity["annotation"].items():
                if (
                    not isinstance(key, str)
                    or not key
                    or not (
                        isinstance(value, str)
                        and value
                        or isinstance(value, list)
                        and value
                        and all(isinstance(item, str) and item for item in value)
                    )
                ):
                    raise ValueError(f"invalid annotation: {met_id}")
        actual = _identity(model.metabolites.get_by_id(met_id))
        states.append((actual == record["before"], actual == record["after"]))

    for reaction_id, record in curation["reactions"].items():
        if not isinstance(record, dict) or set(record) != {"before", "after"}:
            raise ValueError(f"incomplete before/after: {reaction_id}")
        for stoich in record.values():
            if (
                not isinstance(stoich, dict)
                or not stoich
                or any(
                    met_id not in model.metabolites
                    or type(value) not in (int, float)
                    or not isfinite(value)
                    or value == 0
                    for met_id, value in stoich.items()
                )
            ):
                raise ValueError(f"invalid stoichiometry: {reaction_id}")
        if {key: val for key, val in record["before"].items() if key != PROTON} != {
            key: val for key, val in record["after"].items() if key != PROTON
        }:
            raise ValueError(f"out-of-scope stoichiometry: {reaction_id}")
        actual = _stoichiometry(model.reactions.get_by_id(reaction_id))
        states.append((actual == record["before"], actual == record["after"]))
    already_applied = all(after for _, after in states)
    if not already_applied and not all(before for before, _ in states):
        raise ValueError("R608 found a mixed or unexpected before/after state")

    candidate = model.copy()
    if not already_applied:
        for met_id, record in curation["metabolites"].items():
            met = candidate.metabolites.get_by_id(met_id)
            for field, value in record["after"].items():
                setattr(met, field, deepcopy(value))
        for reaction_id, record in curation["reactions"].items():
            reaction = candidate.reactions.get_by_id(reaction_id)
            delta = record["after"].get(PROTON, 0) - record["before"].get(PROTON, 0)
            if delta:
                reaction.add_metabolites(
                    {candidate.metabolites.get_by_id(PROTON): delta}
                )
    from .pipeline import validate_model

    validate_model(candidate, model, no_solve=True)
    return candidate, sum(not after for _, after in states)
