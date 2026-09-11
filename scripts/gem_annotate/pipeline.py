"""Compose the canonical builder, guarded lipid candidate and dFBA comparison."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

import cobra
from cobra.io import read_sbml_model
from cobra.util.solver import linear_reaction_coefficients

from .config import REPO_ROOT
from .essentiality_evidence import sha256_file
from .microspecies import _reaction_balance_record
from .sbml import write_deterministic_sbml_model
from .execution import guarded_execution


@guarded_execution
def validate_model(model, baseline=None, *, no_solve: bool = False) -> dict:
    """Check structural and chemistry regressions before simulation or export."""
    for collection in (model.reactions, model.metabolites):
        if len({item.id for item in collection}) != len(collection):
            raise ValueError("model contains duplicate identifiers")
    for reaction in model.reactions:
        if not all(math.isfinite(value) for value in reaction.bounds) or reaction.lower_bound > reaction.upper_bound:
            raise ValueError(f"invalid bounds: {reaction.id}")
    objective = {reaction.id: float(value) for reaction, value in linear_reaction_coefficients(model).items()}
    if not objective:
        raise ValueError("model has no reaction objective")
    chemistry = {reaction.id: _reaction_balance_record(reaction) for reaction in model.reactions}
    if baseline is not None:
        expected = {reaction.id: float(value) for reaction, value in linear_reaction_coefficients(baseline).items()}
        if objective != expected or model.objective.direction != baseline.objective.direction:
            raise ValueError("candidate changed the objective")
        before = {reaction.id: _reaction_balance_record(reaction) for reaction in baseline.reactions}
        regressions = []
        for reaction_id, record in chemistry.items():
            previous = before.get(reaction_id)
            if previous is not None and previous["status"] != "boundary" and record["status"] == "boundary":
                raise ValueError(f"internal reaction became a boundary: {reaction_id}")
            if record["status"] in {"balanced", "boundary"}:
                continue
            if previous is None or previous["status"] in {"balanced", "boundary"}:
                regressions.append(reaction_id)
            elif previous["status"] == "imbalanced" and (
                record["status"] != "imbalanced" or any(
                    abs(value) > abs(previous["residual"].get(element, 0)) + 1e-9
                    for element, value in record["residual"].items()
                )
            ):
                regressions.append(reaction_id)
        if regressions:
            raise ValueError(f"new mass/charge imbalance: {regressions}")
        incomplete = [
            reaction.id for reaction in model.reactions
            if reaction.id not in baseline.reactions
            and any(not metabolite.formula or metabolite.elements is None or metabolite.charge is None
                    for metabolite in reaction.metabolites)
        ]
        if incomplete:
            raise ValueError(f"new reactions lack chemical identities: {incomplete}")
    solution = None
    if not no_solve:
        solution = model.optimize()
        if solution.status != "optimal" or not math.isfinite(solution.objective_value):
            raise ValueError(f"model solve is not optimal: {solution.status}")
    return {
        "reactions": len(model.reactions), "metabolites": len(model.metabolites),
        "objective": objective, "objective_direction": model.objective.direction,
        "solver_status": solution.status if solution is not None else "not_run",
        "objective_value": float(solution.objective_value) if solution is not None else None,
        "imbalanced_reaction_ids": sorted(key for key, value in chemistry.items() if value["status"] == "imbalanced"),
        "uncheckable_reaction_ids": sorted(key for key, value in chemistry.items() if value["status"] in {"error", "uncheckable"}),
        "incomplete_formula_ids": sorted(m.id for m in model.metabolites if not m.formula or m.elements is None),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--model", type=Path, help="Existing baseline (default: repository model.xml)")
    source.add_argument("--rebuild", action="store_true", help="Rebuild a canonical copy from raw inputs, offline")
    parser.add_argument("--research-root", type=Path)
    parser.add_argument("--lipid", choices=("none", "strict-sn"), default="none")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--no-solve", action="store_true", help="Static baseline validation; no optimization")
    parser.add_argument("--solver", default="glpk")
    parser.add_argument("--experimental", type=Path, help="Positive-only essentiality reference for dFBA")
    parser.add_argument("--dynamic-medium", type=Path)
    parser.add_argument("--hours", type=float, default=24.0)
    parser.add_argument("--step-hours", type=float, default=0.1)
    parser.add_argument("--initial-biomass", type=float, default=0.05)
    parser.add_argument("--growth-cutoff", type=float, default=0.01)
    args = parser.parse_args(argv)
    if bool(args.experimental) != bool(args.dynamic_medium):
        parser.error("--experimental and --dynamic-medium must be supplied together")
    if args.no_solve and (args.lipid != "none" or args.experimental):
        parser.error("--no-solve supports baseline construction and static validation only")
    out = args.output_dir.resolve()
    # Exclusive directory creation also protects symlink aliases and previous runs.
    out.mkdir(parents=True, exist_ok=False)
    manifest = {
        "status": "running", "lipid": args.lipid, "production_gate_passed": False,
        "solver": args.solver, "cobra_version": cobra.__version__,
        "python_version": sys.version.split()[0],
        "entrypoint_sha256": sha256_file(Path(__file__)),
    }
    exit_code = 0
    try:
        # An exported code snapshot has no Git metadata; do not report its parent repository.
        manifest["git_commit"] = (subprocess.check_output(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"], text=True,
        ).strip() if (REPO_ROOT / ".git").exists() else None)
        model_path = (args.model or REPO_ROOT / "model.xml").resolve()
        if args.rebuild:
            model_path = out / "baseline.xml"
            command = [sys.executable, "-m", "scripts.gem_annotate", "--offline", "--canonical-copy", "--output-model", str(model_path)]
            if args.no_solve:
                command.append("--no-solve")
            if args.research_root:
                command += ["--research-root", str(args.research_root.resolve())]
            subprocess.run(command, cwd=REPO_ROOT, check=True)
        source_sha = sha256_file(model_path)
        manifest["source"] = {"path": str(model_path), "sha256": source_sha}
        baseline = read_sbml_model(str(model_path))
        baseline.solver = args.solver
        manifest["baseline_validation"] = validate_model(baseline, no_solve=args.no_solve)
        candidate = baseline
        if args.lipid == "strict-sn":
            from scripts.lp_sn12_candidate import CURATION_PATH, build_candidate, source_fingerprint

            curation = json.loads(CURATION_PATH.read_text())
            if source_sha != curation["source"]["model_sha256"]:
                raise ValueError("source model SHA drifted; refusing strict-sn candidate")
            original = source_fingerprint(baseline)
            candidate = build_candidate(baseline)
            if source_fingerprint(baseline) != original:
                raise ValueError("candidate construction mutated the baseline")
            manifest["curation_sha256"] = sha256_file(CURATION_PATH)
            manifest["activation_state"] = "blocked"
            manifest["remaining_blockers"] = curation["blockers"]
        manifest["candidate_validation"] = validate_model(candidate, baseline, no_solve=args.no_solve)
        candidate_path = out / "candidate.xml"
        write_deterministic_sbml_model(candidate, candidate_path)
        manifest["candidate"] = {"path": str(candidate_path), "sha256": sha256_file(candidate_path)}
        if args.experimental:
            from scripts.dfba_new_fn_essentiality import main as run_dfba

            command = ["--baseline", str(model_path), "--candidate", str(candidate_path),
                       "--experimental", str(args.experimental.resolve()), "--dynamic-medium", str(args.dynamic_medium.resolve()),
                       "--output-dir", str(out / "dfba"), "--solver", args.solver]
            for name in ("hours", "step_hours", "initial_biomass", "growth_cutoff"):
                command += ["--" + name.replace("_", "-"), str(getattr(args, name))]
            run_dfba(command)
            manifest["dfba_summary_sha256"] = sha256_file(out / "dfba" / "summary.json")
        if sha256_file(model_path) != source_sha:
            raise ValueError("source file changed during the run")
        manifest["status"] = "complete"
    except (Exception, SystemExit) as error:
        exit_code = error.code if isinstance(error, SystemExit) and isinstance(error.code, int) and error.code else 1
        manifest.update(status="failed", error_type=type(error).__name__, error=str(error))
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"status": manifest["status"], "manifest": str(out / "manifest.json")}, sort_keys=True))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
