"""User-authorized restores of reviewed reaction definitions that the metadata selection overwrites.

The 2026-09-10 metadata reaction selection copies fields from an older build and so replaces some
reviewed definitions. Each change set below is an opt-in step, applied after the selection, that sets
its reactions back to the reviewed definition. The decision and its evidence are data:

- SELECTION-GPR-RESTORE-20261019 (``selection_cleared_gpr_restore.json``): the empty GPRs on R612,
  R2274, R1302 and R1303, and the unbalanced R2274 (no H2O) and R1303 (extra H+) equations.
- SELECTION-GPR-RESTORE-R671-20261019 (``selection_overwritten_gpr_restore.json``): R671 without the
  mis-annotated GatB gene YALI1E11370g.
- SELECTION-CPA-URA2-RESTORE-20261019 (``selection_cpa_ura2_restore.json``): EG-GPR-001 on R159
  (URA2-only, channelled) and R190 (CPA2 and CPA1).

Every reaction, species, gene, context reaction and pool is checked before any edit; a reaction must be
exactly in its cleared or its restored state, or the step fails closed. A restored equation must
balance, except an unchanged equation whose residual the record states. Genes and metabolites are
not added, so a rerun, also on a saved and reloaded model, is a no-op.
"""

import html
import json
from collections import Counter

from cobra.core.gene import GPR

from .config import MODEL

SPEC_PATH = MODEL.curation_file("selection_cleared_gpr_restore.json")
TARGETS = ("R612", "R2274", "R1302", "R1303")
CLEARED_GPR_RESTORE = "SELECTION-GPR-RESTORE-20261019"
R671_GPR_RESTORE = "SELECTION-GPR-RESTORE-R671-20261019"
CPA_URA2_RESTORE = "SELECTION-CPA-URA2-RESTORE-20261019"
CHANGE_SETS = {
    CLEARED_GPR_RESTORE: (SPEC_PATH, TARGETS),
    R671_GPR_RESTORE: (MODEL.curation_file("selection_overwritten_gpr_restore.json"), ("R671",)),
    CPA_URA2_RESTORE: (MODEL.curation_file("selection_cpa_ura2_restore.json"), ("R159", "R190")),
}
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


def _check_context(model, context):
    """Reactions the change relies on but does not edit must still match the record."""
    for rid, expected in context.items():
        if rid not in model.reactions:
            raise ValueError(f"{rid}: context reaction missing; no edits applied")
        reaction = model.reactions.get_by_id(rid)
        stoichiometry, bounds, gpr = _state(reaction)
        ec = reaction.annotation.get("ec-code")
        actual = {"name": reaction.name, "ec-code": sorted([ec] if isinstance(ec, str) else ec or []),
                  "stoichiometry": stoichiometry, "bounds": bounds, "gpr": gpr}
        for key, value in expected.items():
            if actual[key] != (sorted(value) if key == "ec-code" else value):
                raise ValueError(f"{rid}: context {key} differs from the curated record; no edits applied")


def _check_pools(model, pools, targets):
    """A pool the change relies on must have exactly the recorded users outside the targets."""
    for mid, expected in pools.items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid}: pool metabolite missing; no edits applied")
        met = model.metabolites.get_by_id(mid)
        users = {r.id: float(r.metabolites[met]) for r in met.reactions if r.id not in targets}
        if users != expected:
            raise ValueError(f"{mid}: reactions using the pool differ from the curated record; no edits applied")


def apply_selection_gpr_restore(model, spec=None, change_id=CLEARED_GPR_RESTORE):
    path, targets = CHANGE_SETS[change_id]
    spec = json.loads(path.read_text()) if spec is None else spec
    if (spec.get("schema_version") != 1 or spec.get("change_id") != change_id
            or spec.get("status") != "user_authorized_restore" or tuple(sorted(spec["reactions"])) != tuple(sorted(targets))):
        raise ValueError("Unexpected selection GPR restore scope")
    species = spec["species"]
    for mid, row in species.items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid}: species missing; no edits applied")
        m = model.metabolites.get_by_id(mid)
        if (m.formula, m.charge, m.compartment) != (row["formula"], row["charge"], row["compartment"]):
            raise ValueError(f"{mid}: species differs from the curated record; no edits applied")
    _check_context(model, spec.get("context", {}))
    _check_pools(model, spec.get("pool_context", {}), targets)
    pending = []
    for rid in targets:
        row = spec["reactions"][rid]
        if rid not in model.reactions:
            raise ValueError(f"{rid}: reaction missing; no edits applied")
        reaction = model.reactions.get_by_id(rid)
        if (set(row["cleared"]["stoichiometry"]) | set(row["restored"]["stoichiometry"])) - set(species):
            raise ValueError(f"{rid}: curated species list is incomplete")
        # Only an equation this step leaves unchanged may keep a residual, and only the recorded one.
        residual = row.get("unchanged_equation_residual", {})
        if residual and row["restored"]["stoichiometry"] != row["cleared"]["stoichiometry"]:
            raise ValueError(f"{rid}: only an unchanged equation may keep a residual")
        if _imbalance(row["restored"]["stoichiometry"], species) != residual:
            raise ValueError(f"{rid}: restored equation is not balanced; no edits applied")
        genes = row["genes"] if "genes" in row else [row["gene"]]
        for gene in genes:
            if gene["id"] not in model.genes:
                raise ValueError(f"{rid}: gene {gene['id']} missing; no edits applied")
            g = model.genes.get_by_id(gene["id"])
            if "name" in gene and g.name != gene["name"]:
                raise ValueError(f"{rid}: gene {gene['id']} name differs; no edits applied")
            if any(g.annotation.get(k) not in (v, [v]) for k, v in gene.items() if k not in ("id", "name")):
                raise ValueError(f"{rid}: gene {gene['id']} identity differs; no edits applied")
        if GPR.from_string(row["restored"]["gpr"]).genes != {gene["id"] for gene in genes}:
            raise ValueError(f"{rid}: curated GPR and genes disagree")
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
        row = spec["reactions"][reaction.id]
        if not _same(reaction, row["restored"]):
            raise RuntimeError(f"{reaction.id}: restore did not reach the curated state")
        record = {"item": reaction.id, "status": "applied" if changed else "already_correct",
                  "gpr": reaction.gene_reaction_rule, "evidence_status": row["evidence_status"]}
        if row.get("unchanged_equation_residual"):
            record["unchanged_equation_residual"] = row["unchanged_equation_residual"]
        records.append(record)
    return {"change_id": spec["change_id"], "reactions": records,
            "status": "applied" if any(r["status"] == "applied" for r in records) else "already_correct"}


def restore_after_selection(model, selection, change_id=CLEARED_GPR_RESTORE):
    """Apply only on a build whose metadata selection completed without conflicts."""
    if not selection.get("complete"):
        raise ValueError("The selection GPR restore requires a complete metadata reaction selection")
    return apply_selection_gpr_restore(model, change_id=change_id)
