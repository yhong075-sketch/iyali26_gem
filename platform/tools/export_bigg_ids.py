"""Export a copy of a candidate model with the IDs of Andrew McNaughton's mcnaugh_pass crosswalk.

The new IDs are mcnaugh_pass model IDs in BiGG style. Most are built from the crosswalk's BiGG value, but some differ
from it, and many entities have no BiGG value at all, so the copy is not a BiGG model. The working model keeps its legacy IDs (user decision 2026-11-03); this tool writes a separate, renamed copy.
Only IDs change: names, annotations, formulas, charges, stoichiometry, bounds and GPRs are kept, and every renamed
entity records its legacy ID in notes["legacy_id"]. Compartment codes in the new IDs follow yeast-GEM / iYali
(curated policy), and entities on the policy's skip list (their reviewed BiGG annotation row was not written) keep
their legacy ID.

Unlike the crosswalk's own translate_to_bigg.py, it never writes the crosswalk's annotation values: those were
reviewed separately and enter the model through the builder (--bigg-crosswalk-annotations).

Writes a new model file, a rename report and a JSON sidecar; refuses existing outputs and the canonical model.
Run from platform/:
    python -m tools.export_bigg_ids MODEL.xml --crosswalk-dir DIR --output OUT.xml
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model

POLICY_PATH = MODEL.curation_file("bigg_export_policy_20261103.json")
TMP = "__bigg_export_tmp__"
SKIP_REASON = "reviewed BiGG ID names another entity"


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs(crosswalk_dir: Path, policy: dict | None = None) -> tuple[dict, dict[str, list[dict]]]:
    """Return the curated policy and the crosswalk rows; the crosswalk must be the reviewed version."""
    policy = policy or json.loads(POLICY_PATH.read_text())
    if (policy.get("schema_version") != 1 or policy.get("policy_id") != "BIGG-EXPORT-20261103"
            or policy.get("status") != "user_authorized"):
        raise ValueError("Unexpected BiGG export policy scope/status")
    rows = {}
    for name, expected in policy["crosswalk_sha256"].items():
        path = Path(crosswalk_dir) / name
        if sha256(path) != expected:
            raise ValueError(f"{name} differs from the reviewed crosswalk ({expected[:12]}...)")
        rows[name] = list(csv.DictReader(path.open()))
    return policy, rows


def yeast_gem_id(new_id: str, compartments: set[str], replacements: dict, metabolite: bool) -> str:
    """Replace the crosswalk's compartment marker by the policy's code for each compartment the entity touches."""
    for compartment, rule in replacements.items():
        if compartment not in compartments:
            continue
        old, new = rule["from"], rule["to"]
        if metabolite:
            if not new_id.endswith(f"_{old}"):
                raise ValueError(f"{new_id}: metabolite in {compartment} without the _{old} suffix")
            new_id = new_id[:-len(old)] + new
        elif f"_{old}_" in new_id:
            new_id = new_id.replace(f"_{old}_", f"_{new}_")
        elif new_id.endswith(old):
            new_id = new_id[:-len(old)] + new
    return new_id


def rename_plan(model, policy: dict, rows: dict[str, list[dict]]):
    """Return metabolite and reaction ID maps plus a report row per renamed or kept entity.

    The crosswalk may predate the model (user decision 2026-11-03): a reaction to be renamed must still have the
    crosswalk's equation, while changed reactions that are not renamed and entities newer than the crosswalk keep
    their legacy IDs and are reported.
    """
    reactions = rows["reaction_crosswalk.csv"]
    renamed_statuses = set(policy["renamed_reaction_statuses"])
    unknown = [r["upstream_id"] for r in reactions if r["upstream_id"] not in model.reactions]
    unknown += [r["upstream_id"] for r in rows["metabolite_crosswalk.csv"] if r["upstream_id"] not in model.metabolites]
    if unknown:
        raise ValueError(f"Crosswalk entities not in the model: {unknown[:5]} ({len(unknown)})")
    changed = [r for r in reactions
               if model.reactions.get_by_id(r["upstream_id"]).build_reaction_string() != r["upstream_equation"]]
    stale = [r["upstream_id"] for r in changed if r["mcnaugh_id"] and r["status"] in renamed_statuses]
    if stale:
        raise ValueError(f"Reactions to be renamed changed since the crosswalk was built: {stale[:5]} ({len(stale)})")
    skip = {(s["entity"], s["entity_id"]): s for s in policy["skip_rename"]}
    # A metabolite re-identified in place keeps its ID and equations; its name, formula or charge tells it apart.
    met_rows = {r["upstream_id"]: r for r in rows["metabolite_crosswalk.csv"]}
    reidentified = [mid for mid, r in met_rows.items() if r["mcnaugh_id"] and ("metabolite", mid) not in skip
                    and (model.metabolites.get_by_id(mid).name, model.metabolites.get_by_id(mid).formula or "",
                         str(model.metabolites.get_by_id(mid).charge))
                    != (r["upstream_name"], r["upstream_formula"], r["upstream_charge"])]
    if reidentified:
        raise ValueError(f"Metabolites to be renamed changed name, formula or charge since the crosswalk was built: "
                         f"{reidentified[:5]} ({len(reidentified)})")
    replacements = policy["compartment_code_replacements"]
    maps = {"metabolite": {}, "reaction": {}}
    report, skipped = [], set()
    candidates = [("metabolite", r["upstream_id"], r["mcnaugh_id"]) for r in rows["metabolite_crosswalk.csv"]
                  if r["mcnaugh_id"]]
    candidates += [("reaction", r["upstream_id"], r["mcnaugh_id"]) for r in reactions
                   if r["mcnaugh_id"] and r["status"] in policy["renamed_reaction_statuses"]]
    for entity, legacy, crosswalk_id in candidates:
        collection = model.metabolites if entity == "metabolite" else model.reactions
        obj = collection.get_by_id(legacy)
        if (entity, legacy) in skip:
            report.append((entity, legacy, "", f"kept: {skip[(entity, legacy)].get('reason', SKIP_REASON)}"))
            skipped.add((entity, legacy))
            continue
        compartments = {obj.compartment} if entity == "metabolite" else {m.compartment for m in obj.metabolites}
        new_id = yeast_gem_id(crosswalk_id, compartments, replacements, entity == "metabolite")
        if new_id == legacy:
            report.append((entity, legacy, new_id, "same ID"))
            continue
        taken = set(maps[entity].values())
        if new_id in taken or (new_id in collection and new_id not in maps[entity] and new_id != legacy):
            report.append((entity, legacy, "", f"kept: target {new_id} already used"))
            continue
        maps[entity][legacy] = new_id
        report.append((entity, legacy, new_id, "renamed"))
    missing_skip = set(skip) - skipped
    if missing_skip:
        raise ValueError(f"Skip-list entities not in the crosswalk: {sorted(missing_skip)}")
    report += [("reaction", r["upstream_id"], "", "kept: equation changed since the crosswalk") for r in changed]
    in_crosswalk = {"metabolite": {r["upstream_id"] for r in rows["metabolite_crosswalk.csv"]},
                    "reaction": {r["upstream_id"] for r in reactions}}
    report += [(entity, x.id, "", "kept: not in the crosswalk")
               for entity, collection in (("metabolite", model.metabolites), ("reaction", model.reactions))
               for x in collection if x.id not in in_crosswalk[entity]]
    return maps, report


def definition(model, maps=None):
    """Everything an export may not change, keyed by the (renamed) IDs; renamed entities' legacy_id note dropped."""
    maps = maps or {"metabolite": {}, "reaction": {}}
    met_id = lambda m: maps["metabolite"].get(m.id, m.id)
    notes = lambda x: {k: v for k, v in x.notes.items() if k != "legacy_id"}
    return {
        "metabolites": {met_id(m): (m.name, m.formula, m.charge, m.compartment, m.annotation, notes(m))
                        for m in model.metabolites},
        "reactions": {maps["reaction"].get(r.id, r.id): (
            r.name, {met_id(m): c for m, c in r.metabolites.items()}, r.bounds, r.gene_reaction_rule,
            r.annotation, notes(r), r.objective_coefficient) for r in model.reactions},
        "genes": {g.id: (g.name, g.annotation) for g in model.genes},
        "compartments": dict(model.compartments),
        "objective_direction": model.objective.direction,
    }


def apply_renames(model, maps) -> None:
    """Rename through temporary IDs so swapped or chained targets never collide; record legacy IDs."""
    for entity, collection in (("metabolite", model.metabolites), ("reaction", model.reactions)):
        mapping = maps[entity]
        for legacy in mapping:
            collection.get_by_id(legacy).id = TMP + legacy
        model.repair()
        for legacy, new_id in mapping.items():
            obj = collection.get_by_id(TMP + legacy)
            obj.id = new_id
            obj.notes = {**obj.notes, "legacy_id": legacy}
        model.repair()


def export(model_path: Path, crosswalk_dir: Path, output: Path) -> dict:
    output = Path(output)
    unverified = output.with_name(f"{output.stem}.unverified{output.suffix}")
    report_path, summary_path = output.with_suffix(".rename_report.csv"), output.with_suffix(".export.json")
    for path in (output, unverified, report_path, summary_path):
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"Refusing to overwrite {path}")
    if output.resolve() == MODEL.canonical_model.resolve():
        raise ValueError("Refusing to write the canonical model")
    policy, rows = load_inputs(crosswalk_dir)
    model = read_sbml_model(Path(model_path))
    maps, report = rename_plan(model, policy, rows)
    expected = definition(model, maps)
    apply_renames(model, maps)
    # The file takes its final name only after the reload check; a failed check leaves it marked unverified.
    write_deterministic_sbml_model(model, unverified)
    loaded = read_sbml_model(unverified)
    if definition(loaded) != expected:
        raise ValueError(f"Export changed more than IDs: {unverified}")
    legacy = {(e, x.notes.get("legacy_id")) for e, xs in (("metabolite", loaded.metabolites),
                                                          ("reaction", loaded.reactions)) for x in xs}
    if any((e, old) not in legacy for e, mapping in maps.items() for old in mapping):
        raise ValueError(f"Export lost a legacy_id note: {unverified}")
    unverified.rename(output)
    with report_path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["entity", "legacy_id", "new_id", "action"])
        writer.writerows(report)
    summary = {
        "input_model": {"path": str(model_path), "sha256": sha256(model_path)},
        "output_model": {"path": str(output), "sha256": sha256(output)},
        "rename_report_sha256": sha256(report_path),
        "policy": {"path": str(POLICY_PATH), "sha256": sha256(POLICY_PATH)},
        "tool": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__))},
        "crosswalk_sha256": policy["crosswalk_sha256"],
        "counts": {f"{e}: {a}": n for (e, a), n in sorted(collections.Counter((e, a) for e, _, _, a in report).items())},
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("model", type=Path)
    parser.add_argument("--crosswalk-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(export(args.model, args.crosswalk_dir, args.output), indent=2))


if __name__ == "__main__":
    main()
