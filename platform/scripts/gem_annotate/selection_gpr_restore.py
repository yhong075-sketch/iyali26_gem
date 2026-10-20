"""User-authorized 2026-10-19 restore of four reviewed GPRs cleared by the metadata selection.

The 2026-09-10 metadata reaction selection writes an empty GPR on R612, R2274, R1302 and R1303,
and an unbalanced equation on R2274 (no H2O) and R1303 (extra H+). This opt-in step, applied after
the selection, sets each reaction back to its pre-selection GPR and balanced equation. The
decision and its evidence are data in ``selection_cleared_gpr_restore.json``.

Every reaction, species and gene is checked before any edit; a reaction must be exactly in its
cleared or its restored state, or the step fails closed. Genes and metabolites are not added, so
a rerun, also on a saved and reloaded model, is a no-op.
"""

import html
import json
from collections import Counter

from .config import MODEL

SPEC_PATH = MODEL.curation_file("selection_cleared_gpr_restore.json")
TARGETS = ("R612", "R2274", "R1302", "R1303")
_NOTE_KEYS = ("selection_gpr_restore", "selection_gpr_restore_evidence",
              "selection_gpr_restore_evidence_status", "selection_gpr_restore_sources")


def _state(reaction):
    return ({m.id: float(c) for m, c in reaction.metabolites.items()},
            [float(b) for b in reaction.bounds], reaction.gene_reaction_rule.strip())


def _same(reaction, state):
    return _state(reaction) == (state["stoichiometry"], state["bounds"], state["gpr"])


def _imbalance(stoichiometry, species):
    from cobra.core.formula import Formula

    total, charge = Counter(), 0.0
    for mid, coefficient in stoichiometry.items():
        for element, n in Formula(species[mid]["formula"]).elements.items():
            total[element] += n * coefficient
        charge += species[mid]["charge"] * coefficient
    return {k: v for k, v in {**total, "charge": charge}.items() if abs(v) > 1e-9}


def _notes(spec, row):
    return dict(zip(_NOTE_KEYS, (spec["change_id"], row["evidence"], row["evidence_status"],
                                 "; ".join(row["sources"]))))


def _note_differs(reaction, key, value):
    # SBML reload returns XML entities (e.g. &apos;) as literal text; compare unescaped.
    return html.unescape(str(reaction.notes.get(key, ""))) != html.unescape(value)


def apply_selection_gpr_restore(model, spec=None):
    spec = json.loads(SPEC_PATH.read_text()) if spec is None else spec
    if (spec.get("schema_version") != 1 or spec.get("change_id") != "SELECTION-GPR-RESTORE-20261019"
            or spec.get("status") != "user_authorized_restore" or tuple(sorted(spec["reactions"])) != tuple(sorted(TARGETS))):
        raise ValueError("Unexpected selection GPR restore scope")
    species = spec["species"]
    for mid, row in species.items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid}: species missing; no edits applied")
        m = model.metabolites.get_by_id(mid)
        if (m.formula, m.charge, m.compartment) != (row["formula"], row["charge"], row["compartment"]):
            raise ValueError(f"{mid}: species differs from the curated record; no edits applied")
    pending = []
    for rid in TARGETS:
        row = spec["reactions"][rid]
        if rid not in model.reactions:
            raise ValueError(f"{rid}: reaction missing; no edits applied")
        reaction = model.reactions.get_by_id(rid)
        if (set(row["cleared"]["stoichiometry"]) | set(row["restored"]["stoichiometry"])) - set(species):
            raise ValueError(f"{rid}: curated species list is incomplete")
        if _imbalance(row["restored"]["stoichiometry"], species):
            raise ValueError(f"{rid}: restored equation is not balanced; no edits applied")
        gene = row["gene"]
        if gene["id"] not in model.genes:
            raise ValueError(f"{rid}: gene {gene['id']} missing; no edits applied")
        g = model.genes.get_by_id(gene["id"])
        if "name" in gene and g.name != gene["name"]:
            raise ValueError(f"{rid}: gene {gene['id']} name differs; no edits applied")
        if any(g.annotation.get(k) not in (v, [v]) for k, v in gene.items() if k not in ("id", "name")):
            raise ValueError(f"{rid}: gene {gene['id']} identity differs; no edits applied")
        if row["restored"]["gpr"] != gene["id"]:
            raise ValueError(f"{rid}: curated GPR and gene disagree")
        notes = _notes(spec, row)
        if any(k in reaction.notes and _note_differs(reaction, k, v) for k, v in notes.items()):
            raise ValueError(f"{rid}: conflicting restore notes; no edits applied")
        if _same(reaction, row["restored"]):
            pending.append((reaction, None, notes))
        elif _same(reaction, row["cleared"]):
            pending.append((reaction, row["restored"], notes))
        else:
            raise ValueError(f"{rid}: definition is neither the cleared nor the restored state; no edits applied")
    records = []
    for reaction, target, notes in pending:
        missing_notes = {k: v for k, v in notes.items() if _note_differs(reaction, k, v)}
        changed = target is not None or bool(missing_notes)
        if target is not None:
            current = _state(reaction)[0]
            delta = {mid: target["stoichiometry"].get(mid, 0.0) - current.get(mid, 0.0)
                     for mid in set(current) | set(target["stoichiometry"])}
            reaction.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in delta.items() if c},
                                     combine=True)
            reaction.gene_reaction_rule = target["gpr"]
        reaction.notes = {**reaction.notes, **missing_notes}
        if not _same(reaction, spec["reactions"][reaction.id]["restored"]):
            raise RuntimeError(f"{reaction.id}: restore did not reach the curated state")
        records.append({"item": reaction.id, "status": "applied" if changed else "already_correct",
                        "gpr": reaction.gene_reaction_rule,
                        "evidence_status": spec["reactions"][reaction.id]["evidence_status"]})
    return {"change_id": spec["change_id"], "reactions": records,
            "status": "applied" if any(r["status"] == "applied" for r in records) else "already_correct"}


def restore_after_selection(model, selection):
    """Apply only on a build whose metadata selection completed without conflicts."""
    if not selection.get("complete"):
        raise ValueError("The selection GPR restore requires a complete metadata reaction selection")
    return apply_selection_gpr_restore(model)
