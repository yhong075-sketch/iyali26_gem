"""User-authorized 2026-10-13 removal of R132 (matrix alcohol dehydrogenase).

The decision and its evidence are data in ``r132_mito_adh.json``. Before any edit the step checks
the exact R132 definition, the two diffusion reactions it pairs with (R1114, R997), the reactions
that use matrix ethanol, and every R132 gene (presence, recorded RefSeq, R132-only membership),
and fails closed on anything else. Genes and metabolites are kept, so a rerun is a no-op.
"""

import json

from cobra.core.gene import GPR

from .coq9 import boolean_key
from .config import MODEL

SPEC_PATH = MODEL.curation_file("r132_mito_adh.json")


def _matches(reaction, state):
    return ({m.id: float(c) for m, c in reaction.metabolites.items()} == state["stoichiometry"]
            and list(reaction.bounds) == state["bounds"]
            and boolean_key(reaction.gpr.body) == boolean_key(GPR.from_string(state["gpr"]).body))


def _usage(model, mid):
    if mid not in model.metabolites:
        raise ValueError(f"R132 curation species missing: {mid}")
    return sorted(r.id for r in model.metabolites.get_by_id(mid).reactions)


def apply_r132_mito_adh(model, spec=None):
    spec = json.loads(SPEC_PATH.read_text()) if spec is None else spec
    if (spec.get("schema_version") != 1 or spec.get("change_id") != "R132-MITO-ADH-20261013"
            or spec.get("action") != "remove_reaction" or spec.get("reaction") != "R132"):
        raise ValueError("Unexpected R132 curation scope")
    for rid, state in spec["kept_reactions"].items():
        if rid not in model.reactions or not _matches(model.reactions.get_by_id(rid), state):
            raise ValueError(f"{rid}: definition differs from the curated record; no edits applied")
    for gid, row in spec["genes"].items():
        if gid not in model.genes:
            raise ValueError(f"R132 curation gene missing: {gid}")
        if model.genes.get_by_id(gid).annotation.get("refseq") != row["model_refseq"]:
            raise ValueError(f"R132 curation gene identity differs: {gid}")
    record = {"change_id": spec["change_id"], "status_label": spec["status"], "reaction": spec["reaction"],
              "kept_reactions": sorted(spec["kept_reactions"])}
    if spec["reaction"] not in model.reactions:
        for mid, usage in spec["metabolite_usage"].items():
            if _usage(model, mid) != usage["after"]:
                raise ValueError(f"{mid}: usage differs from the curated after-state")
        return {**record, "status": "already_correct", "orphaned_genes": []}
    reaction = model.reactions.get_by_id(spec["reaction"])
    if not _matches(reaction, spec["before"]):
        raise ValueError("R132: definition differs from the curated before-state; no edits applied")
    for mid, usage in spec["metabolite_usage"].items():
        if _usage(model, mid) != usage["before"]:
            raise ValueError(f"{mid}: usage differs from the curated before-state; no edits applied")
    only_in_r132 = sorted(gid for gid in spec["genes"]
                          if {r.id for r in model.genes.get_by_id(gid).reactions} == {spec["reaction"]})
    if only_in_r132 != sorted(gid for gid, row in spec["genes"].items() if row["only_in_r132"]):
        raise ValueError("R132: R132-only gene membership differs from the curated record; no edits applied")
    model.remove_reactions([reaction], remove_orphans=False)
    return {**record, "status": "applied", "orphaned_genes": only_in_r132,
            "metabolite_usage_after": {mid: _usage(model, mid) for mid in spec["metabolite_usage"]}}
