"""User-authorized 2026-10-19 product correction of R2190 (D-arabinitol:NAD+ 2-oxidoreductase).

The decision and its evidence are data in ``r2190_product.json``. Before any edit the step checks
R2190's exact name, GPR, EC, bounds and stoichiometry, the formula, charge and compartment of every
species, the other reactions of the replaced D-arabinose, that the new D-ribulose ID is unused and
that the curated after-state balances; it fails closed on anything else. A rerun is a no-op.
"""

import json

from cobra import Metabolite

from .config import MODEL

SPEC_PATH = MODEL.curation_file("r2190_product.json")
NEW_SPECIES_FIELDS = ("name", "formula", "charge", "compartment", "annotation", "notes")


def _stoichiometry(reaction):
    return {m.id: float(c) for m, c in reaction.metabolites.items()}


def _usage(model, mid):
    return sorted(r.id for r in model.metabolites.get_by_id(mid).reactions)


def _same_base_id(model, mid):
    return sorted(m.id for m in model.metabolites if m.id.split("[")[0] == mid.split("[")[0])


def _curated_residual(spec):
    """Element and charge residual of the curated after-state, from the curated species alone."""
    species = {**spec["species"], spec["new_species"]["id"]: spec["new_species"]}
    residual = {}
    for mid, coefficient in spec["after_stoichiometry"].items():
        met = Metabolite(mid, formula=species[mid]["formula"], charge=species[mid]["charge"])
        for element, count in {**met.elements, "charge": met.charge}.items():
            residual[element] = residual.get(element, 0) + coefficient * count
    return {k: v for k, v in residual.items() if v}


def apply_r2190_product(model, spec=None):
    spec = json.loads(SPEC_PATH.read_text()) if spec is None else spec
    new, replaced = spec.get("new_species", {}), spec.get("replaced_species", {})
    before, after = spec.get("before_stoichiometry", {}), spec.get("after_stoichiometry", {})
    # Only the product may differ between the curated states; shared coefficients stay equal.
    if (spec.get("schema_version") != 1 or spec.get("curation_id") != "R2190-RIBULOSE-20261019"
            or spec.get("status") != "user_authorized_identity_correction"
            or spec.get("reaction_id") != "R2190" or new.get("id") != "m2061[C_cy]"
            or replaced.get("id") != "m454[C_cy]"
            or set(before) - set(after) != {replaced["id"]} or set(after) - set(before) != {new["id"]}
            or any(before[mid] != after[mid] for mid in set(before) & set(after))):
        raise ValueError("Unexpected R2190 product curation scope")
    if _curated_residual(spec):
        raise ValueError("Curated R2190 after-state is not element and charge balanced")
    if spec["reaction_id"] not in model.reactions:
        raise ValueError("R2190 missing; no edits applied")
    reaction = model.reactions.get_by_id(spec["reaction_id"])
    if (reaction.name != spec["reaction_name"] or reaction.gene_reaction_rule != spec["gpr"]
            or reaction.annotation.get("ec-code") != spec["ec_code"]
            or list(reaction.bounds) != spec["bounds"]):
        raise ValueError("R2190 definition differs from the curated record; no edits applied")
    for mid, props in spec["species"].items():
        if (mid not in model.metabolites
                or {k: getattr(model.metabolites.get_by_id(mid), k) for k in props} != props):
            raise ValueError(f"{mid}: species differs from the curated record; no edits applied")
    if model.metabolites.get_by_id(replaced["id"]).name != replaced["name"]:
        raise ValueError(f"{replaced['id']}: name differs from the curated record; no edits applied")
    if any(k in reaction.notes and reaction.notes[k] != v for k, v in spec["notes"].items()):
        raise ValueError("Conflicting R2190 product notes; no edits applied")
    record = {"item": reaction.id, "curation_id": spec["curation_id"], "evidence_status": spec["status"],
              "new_species": new["id"], "replaced_species": replaced["id"]}

    if _stoichiometry(reaction) == after:
        met = model.metabolites.get_by_id(new["id"])
        if ({k: getattr(met, k) for k in NEW_SPECIES_FIELDS} != {k: new[k] for k in NEW_SPECIES_FIELDS}
                or _same_base_id(model, new["id"]) != [new["id"]]
                or _usage(model, new["id"]) != [reaction.id]
                or _usage(model, replaced["id"]) != replaced["reactions_after"]
                or any(reaction.notes.get(k) != v for k, v in spec["notes"].items())):
            raise ValueError("R2190 after-state differs from the curated record")
        return {**record, "status": "already_correct"}
    if _stoichiometry(reaction) != before:
        raise ValueError("R2190 stoichiometry differs from the curated before-state; no edits applied")
    if _same_base_id(model, new["id"]):
        raise ValueError(f"{new['id']}: species ID already in use; no edits applied")
    if _usage(model, replaced["id"]) != replaced["reactions_before"]:
        raise ValueError(f"{replaced['id']}: usage differs from the curated before-state; no edits applied")

    met = Metabolite(new["id"], formula=new["formula"], name=new["name"], charge=new["charge"],
                     compartment=new["compartment"])
    met.annotation = dict(new["annotation"])
    met.notes = dict(new["notes"])
    model.add_metabolites([met])
    reaction.subtract_metabolites({model.metabolites.get_by_id(replaced["id"]): before[replaced["id"]]})
    reaction.add_metabolites({met: after[new["id"]]})
    reaction.notes = {**reaction.notes, **spec["notes"]}
    if _stoichiometry(reaction) != after:
        raise ValueError("R2190 edit did not reach the curated after-state")
    return {**record, "status": "applied", "stoichiometry": _stoichiometry(reaction),
            "replaced_species_usage_after": _usage(model, replaced["id"])}
