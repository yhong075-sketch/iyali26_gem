"""Evidence pipeline CLI (run from platform/).

    python -m evidence pilot --research-root R [--model M]
    python -m evidence build --research-root R [--model M] [--db PATH] [--report-dir DIR]
    python -m evidence trace --db PATH <verdict id | reaction id>

`build` is deterministic and offline: it reads the committed curation files, the research
workspace (datasets, locus maps, the gem_review cache) and the model, and makes no LLM call.
"""

import argparse
import csv
import dataclasses
import hashlib
import json
from collections import Counter
from pathlib import Path

import yaml

from scripts.gem_annotate.config import load_project_paths
from scripts.gem_annotate.model_layout import MODEL, REPO_ROOT

from . import converters, pilot, store
from .id_map import load_id_map
from .pipeline import ModelIdentity, run, trace
from .verdicts import element_states, model_identity
from .weights import CURATION_DIR, load_weights

DEFAULT_MODEL = ("E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_"
                 "cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_um7_nuclear_diffusion_r2171_1to1_"
                 "aars_thr_ser_split_20261029.xml")
KO_TABLE = MODEL.reports / "aars_thr_ser_split_on_14f73bda_20261029" / "ko_compare.tsv"
LIT_TASK = "gem_review_lit_20261019"
PILOT_CSV = CURATION_DIR / "pilot_reactions.csv"


def _essentiality_state(paths) -> Path:
    """The essentiality datasets, case ledger and dossiers live in the research workspace."""
    return paths.research_root / "state" / "essentiality" / "repository"


def _read_model(path: Path):
    import cobra

    return cobra.io.read_sbml_model(Path(path))   # Path, not str: cobra rejects long str paths on macOS


def _essential_in_simulation(model_sha: str, threshold: float, canonical) -> tuple[set[str] | None, str]:
    """Genes essential in simulation, or None (impact unknown) when the table is for another model."""
    summary = json.loads(KO_TABLE.with_name("ko_compare_summary.json").read_text())
    if summary["models"]["candidate"]["sha256"] != model_sha:
        return None, f"knockout table is for {summary['models']['candidate']['sha256'][:8]}, not this model; impact unknown"
    column = f"candidate_class_{threshold:g}"
    with KO_TABLE.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if column not in rows[0]:
        raise ValueError(f"{KO_TABLE.name} has no column {column} for essentiality_threshold {threshold}")
    return {canonical(r["gene"]) for r in rows if r[column] == "essential"}, f"{KO_TABLE.name} column {column}"


def cmd_pilot(args) -> None:
    paths = load_project_paths(args.research_root, required=True)
    model = _read_model(args.model)
    ess = _essentiality_state(paths)
    rows = pilot.select(MODEL.curation_file("pilot_selection.yaml"), ess / "curation_cases.csv", ess / "evidence",
                        {r.id for r in model.reactions})
    print(pilot.write(rows, PILOT_CSV), len(rows), "reactions")


def cmd_build(args) -> None:
    weights = load_weights()
    paths = load_project_paths(args.research_root, required=True)
    model_path = Path(args.model).resolve()
    model = _read_model(model_path)
    id_map = load_id_map(paths.locus_map, paths.essentiality / "curated_locus_identity_exclusions.csv")
    id_map.add_model_genes({g.id: dict(g.annotation) for g in model.genes})
    rel, commit, sha = model_identity(model_path, REPO_ROOT)
    identity = ModelIdentity(path=rel, commit=commit, sha256=sha)

    pilot_rows = pilot.read(MODEL.curation_file(PILOT_CSV.name))
    reaction_ids = [r["reaction_id"] for r in pilot_rows]
    hub = {r["reaction_id"] for r in pilot_rows if r["group"] == "phosphatidate_hub"}
    essential, impact_note = _essential_in_simulation(sha, weights.essentiality_threshold, id_map.model_gene_canonical)
    states = element_states(model, reaction_ids, id_map.model_gene_canonical, essential, hub)
    genes = {g for s in states.values() for g in s.gene_ids}
    model_species = {rid: {m.id for m in model.reactions.get_by_id(rid).metabolites} for rid in reaction_ids}

    lineage = yaml.safe_load(MODEL.curation_file("model_lineage.yaml").read_text())
    lit = paths.task_outputs / LIT_TASK
    ess = _essentiality_state(paths)
    out = converters.Output()
    out.merge(converters.lineage_sources(lineage))
    protons = frozenset(m.id for m in model.metabolites if m.formula == "H" and m.charge == 1)
    out.merge(converters.parent_model_claims(MODEL.start_model.with_name("iyli21.xml"), reaction_ids, weights, id_map,
                                             lineage, model_species, protons))
    out.merge(converters.model_note_claims(model_path, commit, reaction_ids, weights, id_map, lineage))
    out.merge(converters.w29_13c_claims(paths.research_root / "reference" / "w29_worland2024" / "glucose_central_fluxes.csv",
                                        MODEL.curation_file("w29_13c_reaction_map.csv"), weights, id_map))
    out.merge(converters.essentiality_claims(ess / "cas9_cas12a_fitness.csv", ess / "consensus_essential_genes.csv",
                                             genes, weights, id_map))
    out.merge(converters.legacy_case_sources(ess / "curation_cases.csv", ess / "evidence"))
    out.merge(converters.handwritten_claims(MODEL.curation_file("handwritten_claims.yaml"), lit / "out" / "lit" / "retrieval",
                                            lit / "cache", weights, id_map))

    bundle = run(out.sources, out.passages, out.claims, out.rejected, states, weights, identity)
    runs = sorted({c.extraction_run_id for c in bundle.claims.values()})
    extraction_runs = [{"extraction_run_id": r, "llm": None, "prompt_sha256": None, "temperature": None,
                        "kind": "hand-written (Claude, not human-verified)" if r.startswith("hand-") else "deterministic"}
                       for r in runs]
    id_rows = [dataclasses.asdict(row) for row in id_map.rows.values()]
    db = Path(args.db or paths.research_root / "artifacts" / "evidence_pipeline" / "evidence.sqlite")
    digest = store.write(db, bundle, id_rows, weights.version, extraction_runs)

    spec_sha = hashlib.sha256(MODEL.curation_file("SPEC.md").read_bytes()).hexdigest()
    snapshot = {
        "model": dataclasses.asdict(identity), "weights_version": weights.version, "spec_sha256": spec_sha,
        "content_hash": digest, "impact": impact_note, "id_map_fingerprint": id_map.fingerprint,
        "assertions": [{"id": a.assertion_id, "subject": a.subject_key, "predicate": a.predicate.value,
                        "score": a.score, "n_independent": a.n_independent, "conflict": a.conflict_status,
                        "resolved_by": a.resolved_by, "flag": a.flag, "homology_only": a.homology_only,
                        "inherited_from": a.inherited_from}
                       for a in sorted(bundle.assertions.values(), key=lambda a: a.assertion_id)],
        "verdicts": [{"id": v.verdict_id, "reaction": v.reaction_id, "element": v.element, "status": v.status,
                      "actions": list(v.actions), "priority": v.priority, "priority_note": v.priority_note,
                      "load_bearing": list(bundle.load_bearing.get(v.verdict_id, ()))}
                     for v in sorted(bundle.verdicts.values(), key=lambda v: v.verdict_id)],
    }
    snap_dir = CURATION_DIR / "snapshots"
    snap_dir.mkdir(exist_ok=True)
    snap = snap_dir / f"{commit[:12]}_{sha[:8]}.json"
    snap.write_text(json.dumps(snapshot, indent=1, sort_keys=True) + "\n")

    if args.report_dir:
        _write_report_tables(Path(args.report_dir), bundle, pilot_rows, id_map)
    status = Counter((v.element, v.status) for v in bundle.verdicts.values())
    print(f"db {db}\ncontent hash {digest}\nsnapshot {snap}\nmodel {rel} {commit[:12]} {sha[:12]}\n"
          f"sources {len(bundle.sources)} passages {len(bundle.passages)} claims {len(bundle.claims)} "
          f"rejected {len(bundle.rejected)} pairs {len(bundle.relations)} assertions {len(bundle.assertions)} "
          f"verdicts {len(bundle.verdicts)}\nimpact: {impact_note}")
    for (element, verdict), n in sorted(status.items()):
        print(f"  {element:12s} {verdict:13s} {n}")


def _write_report_tables(report_dir: Path, bundle, pilot_rows, id_map) -> None:
    report_dir.mkdir(parents=True, exist_ok=True)
    group = {r["reaction_id"]: (r["group"], r["origin"]) for r in pilot_rows}
    with (report_dir / "pilot_verdicts.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["reaction_id", "group", "origin", "element", "status", "actions", "priority", "reason",
                    "n_load_bearing"])
        for v in sorted(bundle.verdicts.values(), key=lambda v: (group[v.reaction_id][0], v.reaction_id, v.element)):
            w.writerow([v.reaction_id, *group[v.reaction_id], v.element, v.status, ";".join(v.actions), v.priority,
                        v.reason, len(bundle.load_bearing.get(v.verdict_id, ()))])
    with (report_dir / "id_map_non_yali1.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["model_gene", "canonical", "status", "cardinality", "yali0_links", "uniprot", "ncbigene", "note"])
        for row in sorted(id_map.rows.values(), key=lambda r: r.model_gene):
            if not row.model_gene.startswith("YALI1"):
                w.writerow([row.model_gene, row.canonical or "", row.status, row.cardinality, ";".join(row.yali0),
                            ";".join(row.uniprot), ";".join(row.ncbigene), row.note])
    with (report_dir / "rejected_claims.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["reject_id", "source_id", "reason", "detail"])
        for r in bundle.rejected:
            w.writerow([r.reject_id, r.source_id or "", r.reason, r.detail])


def cmd_trace(args) -> None:
    bundle = store.load(Path(args.db))
    ids = [args.id] if args.id in bundle.verdicts else sorted(
        v.verdict_id for v in bundle.verdicts.values() if v.reaction_id == args.id)
    if not ids:
        raise SystemExit(f"no verdict or reaction {args.id!r} in {args.db}")
    for vid in ids:
        for layer, node, label in trace(bundle, vid):
            print(f"{layer}  {node}  {label}")
        print(f"load-bearing: {list(bundle.load_bearing.get(vid, ()))}\n")


def main() -> None:
    parser = argparse.ArgumentParser(prog="python -m evidence", description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("pilot", "build"):
        p = sub.add_parser(name)
        p.add_argument("--research-root")
        p.add_argument("--model", default=str(MODEL.candidate_file(DEFAULT_MODEL)))
        if name == "build":
            p.add_argument("--db")
            p.add_argument("--report-dir")
    t = sub.add_parser("trace")
    t.add_argument("--db", required=True)
    t.add_argument("id")
    args = parser.parse_args()
    {"pilot": cmd_pilot, "build": cmd_build, "trace": cmd_trace}[args.command](args)


if __name__ == "__main__":
    main()
