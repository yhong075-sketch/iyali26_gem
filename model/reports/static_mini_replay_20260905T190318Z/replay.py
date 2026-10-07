import json
import subprocess
import sys
import time
from datetime import datetime, timezone

worker = r'''
import csv
import importlib.metadata as metadata
import json
import math
import os
from pathlib import Path
import platform
import sys
from hashlib import sha256
from unittest.mock import patch

channel = os.fdopen(os.dup(1), "w", buffering=1)
sink = os.open(os.devnull, os.O_WRONLY)
os.dup2(sink, 1)
os.dup2(sink, 2)
os.close(sink)
report = {"status": "running", "phase": "preflight", "stages": [],
          "primary_solves_requested_upper_bound": 0,
          "purpose": "bounded_static_reference_check",
          "historical_full_environment_reproduction": "not_established",
          "logs_suppressed": True}

def safe(value):
    if isinstance(value, dict):
        return {str(k): safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe(v) for v in value]
    if isinstance(value, (set, frozenset)):
        return [safe(v) for v in sorted(value, key=str)]
    if hasattr(value, "item"):
        return safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return "NaN" if math.isnan(value) else ("+Infinity" if value > 0 else "-Infinity")
    return value

def emit():
    print(json.dumps(safe(report), ensure_ascii=False, allow_nan=False), file=channel, flush=True)

def require(condition, message):
    if not condition:
        report["failure_reason"] = message
        raise RuntimeError("declared_gate_failed")

def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()

def finite_nonnegative(value):
    try:
        return value is not None and math.isfinite(float(value)) and float(value) >= 0
    except (TypeError, ValueError):
        return False

emit()
try:
    manifest_path = Path("/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/docs/baseline_manifest.json")
    manifest = json.loads(manifest_path.read_text())
    report["manifest_sha256"] = digest(manifest_path)
    require(manifest["schema_version"] == "1.0", "Unexpected manifest schema")
    require(manifest["baseline_id"] == "iyali26_static_reference_provisional_20260905", "Unexpected baseline")
    artifacts = manifest["artifacts"]
    required = ("historical_executed_model", "canonical_model", "experimental_positive",
                "medium", "strain_profile", "historical_per_gene", "historical_run_manifest",
                "historical_summary")
    report["input_hashes"] = {}
    for key in required:
        item = artifacts[key]
        observed = digest(item["path"])
        report["input_hashes"][key] = observed
        require(observed == item["sha256"], "Input SHA mismatch: " + key)
    code = manifest["code_identity"]
    worktree = Path(code["recovered_worktree"]).resolve()
    require(str(worktree) == "/private/tmp/iyali26_coq9_wp12_20260904", "Unexpected recovered worktree")
    require(bool(code["source_files"]), "Missing source-file inventory")
    report["source_hashes"] = {}
    for item in code["source_files"]:
        path = Path(item["path"]).resolve()
        require(path == (worktree / item["relative_path"]).resolve(), "Source path mismatch")
        observed = digest(path)
        report["source_hashes"][item["relative_path"]] = {
            "sha256": observed, "historically_hash_recorded": item["historically_hash_recorded"]}
        require(observed == item["sha256"], "Source SHA mismatch: " + item["relative_path"])
    actual_software = {"python": platform.python_version()}
    for package in ("cobra", "optlang", "pandas", "gurobipy"):
        actual_software[package] = metadata.version(package)
    report["software"] = actual_software
    report["python_executable"] = sys.executable
    require(actual_software == code["expected_software"], "Software version mismatch")
    selected = ["YALI1B19079g", "YALI1B20289g", "YALI1E15659g",
                "YALI1A15879g", "YALI1F12842g", "YALI1B09025g"]
    task = manifest["next_task"]
    require(task["selected_gene_ids"] == selected, "Predeclared control list changed")
    require(task["max_primary_solves"] == 14 and task["wall_time_seconds"] == 600, "Budget changed")
    require(task["comparison_abs_tol"] == 1e-8 and task["comparison_rel_tol"] == 1e-6, "Comparison gate changed")
    report["selected_gene_ids"] = selected
    report["new_declared_resource_settings"] = {"Threads": 1, "processes": 1, "wall_time_seconds": 600}
    cuts = (0.01, 0.05, 0.10, 0.15)
    old_manifest = json.loads(Path(artifacts["historical_run_manifest"]["path"]).read_text())
    for old_key, key in (("model", "historical_executed_model"), ("experimental", "experimental_positive"),
                         ("medium", "medium"), ("strain_profile", "strain_profile")):
        require(old_manifest["inputs"][old_key]["sha256"] == artifacts[key]["sha256"], "Historical input binding mismatch: " + key)
    old_summary = json.loads(Path(artifacts["historical_summary"]["path"]).read_text())
    with Path(artifacts["historical_per_gene"]["path"]).open(newline="") as handle:
        old_rows = [r for r in csv.DictReader(handle, delimiter="\t") if r["gene_id"] in selected]
    require(len(old_rows) == 6 and len({r["gene_id"] for r in old_rows}) == 6, "Historical controls missing or duplicated")
    reference = {r["gene_id"]: r for r in old_rows}
    require(all(r["ko_status"] == "optimal" and r["in_model"] == "True" and
                r["experimental_essential"] == "True" and finite_nonnegative(r["ko_growth"]) and
                finite_nonnegative(r["ko_growth_ratio"]) for r in old_rows), "Invalid historical control record")
    call_columns = ("essential_at_1pct", "essential_at_5pct", "essential_at_10pct", "essential_at_15pct")
    historical_calls = {}
    for row in old_rows:
        require(all(row.get(column) in ("True", "False") for column in call_columns), "Missing or invalid historical cutoff calls")
        recorded = [row[column] == "True" for column in call_columns]
        require(recorded == [float(row["ko_growth_ratio"]) < cutoff for cutoff in cuts], "Historical recorded calls disagree with raw-ratio rule")
        historical_calls[row["gene_id"]] = recorded
    require(finite_nonnegative(old_summary["wt_growth"]), "Invalid historical WT")
    os.environ["IYALI26_RESEARCH_ROOT"] = manifest["workspace"]["research_root"]
    report["research_root"] = os.environ["IYALI26_RESEARCH_ROOT"]
    sys.path.insert(0, str(worktree))
    import time
    from optlang import gurobi_interface
    original_backend_optimize = gurobi_interface.Model._optimize
    report["backend_optimization_entries"] = 0
    report["backend_optimization_completed"] = 0
    report["backend_optimization_log"] = []

    def counted_backend_optimize(self, *args, **kwargs):
        require(report["backend_optimization_entries"] < 14, "Actual backend optimization budget exhausted; next call blocked")
        report["backend_optimization_entries"] += 1
        entry = {"number": report["backend_optimization_entries"], "phase": report["phase"],
                 "panel_member": report["stages"][-1].get("requested_ko", "WT") if report["stages"] else "unexpected_preflight",
                 "returned": False}
        report["backend_optimization_log"].append(entry)
        emit()
        started = time.monotonic()
        try:
            status = original_backend_optimize(self, *args, **kwargs)
            entry["returned"] = True
            entry["status"] = status
            report["backend_optimization_completed"] += 1
            return status
        except Exception as exc:
            entry["exception_type"] = type(exc).__name__
            raise
        finally:
            entry["elapsed_seconds"] = time.monotonic() - started
            emit()

    gurobi_interface.Model._optimize = counted_backend_optimize
    from scripts.gem_annotate import validate_essential_genes as validator
    from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
    require(Path(validator.__file__).resolve() == worktree / "scripts/gem_annotate/validate_essential_genes.py", "Wrong module imported")
    experimental = validator.load_experimental(Path(artifacts["experimental_positive"]["path"]), positive_only=True)
    positive_ids = set(experimental.loc[experimental["essential"].eq(True), "gene_id"])
    require(set(selected) <= positive_ids, "Control absent from positive reference")
    report["phase"] = "preflight_passed"
    emit()

    def compare(stage, wt, rows, ref_wt, ref_rows, ref_calls, label):
        checks = {"wt": math.isclose(wt, float(ref_wt), abs_tol=1e-8, rel_tol=1e-6), "genes": {}}
        for row in rows:
            target = ref_rows[row["gene_id"]]
            ratio = float(row["ko_growth_ratio"])
            ref_ratio = float(target["ko_growth_ratio"])
            checks["genes"][row["gene_id"]] = {
                "growth": math.isclose(float(row["ko_growth"]), float(target["ko_growth"]), abs_tol=1e-8, rel_tol=1e-6),
                "ratio": math.isclose(ratio, ref_ratio, abs_tol=1e-8, rel_tol=1e-6),
                "calls": [ratio < c for c in cuts], "reference_calls": ref_calls[row["gene_id"]]}
        stage.setdefault("comparisons", {})[label] = checks
        emit()
        require(checks["wt"] and all(x["growth"] and x["ratio"] and x["calls"] == x["reference_calls"]
                                    for x in checks["genes"].values()), "Numerical or cutoff mismatch: " + label)

    for key in ("historical_executed_model", "canonical_model"):
        stage = {"model_artifact": key, "status": "running"}
        report["stages"].append(stage)
        report["phase"] = key
        emit()
        context = load_effective_simulation_context(
            model_path=artifacts[key]["path"], media_path=artifacts["medium"]["path"],
            strain_profile_path=artifacts["strain_profile"]["path"])
        stage["loaded_input_hashes"] = {key: context.canonical_model_sha256,
                                       "medium": context.medium_sha256,
                                       "strain_profile": context.strain_profile_sha256}
        emit()
        require(all(value == artifacts[name]["sha256"] for name, value in stage["loaded_input_hashes"].items()),
                "Loaded input SHA mismatch; no solve dispatched")
        model = context.model
        model.solver = "gurobi"
        model.solver.problem.Params.Threads = 1
        stage["context"] = context.provenance()
        stage["solver_parameters"] = {name: getattr(model.solver.problem.Params, name)
                                      for name in ("Threads", "FeasibilityTol", "OptimalityTol", "Method", "Seed", "Presolve")}
        stage["optlang_presolve"] = model.solver.configuration.presolve
        require(stage["optlang_presolve"] != "auto", "optlang auto presolve permits implicit retry; stop without changing it")
        stage["cobra_tolerance"] = model.tolerance
        all_ids = {gene.id for gene in model.genes}
        require(set(selected) <= all_ids, "Predeclared control missing from model")
        require(context.strain_overlay_enabled and context.active_medium.get("R1354") == 1000.0 and
                context.active_medium.get("R1189", 0) > 0 and context.active_medium.get("R1219", 0) == 0,
                "Historical runtime medium/strain invariant failed")
        raw_tables = []
        original_deletion = validator.single_gene_deletion
        original_optimize = model.optimize

        def observe_wt(*args, **kwargs):
            require(model.solver.configuration.presolve != "auto", "Unexpected auto presolve before WT")
            report["primary_solves_requested_upper_bound"] += 1
            require(report["primary_solves_requested_upper_bound"] <= 14, "Primary solve budget exceeded")
            emit()
            solution = original_optimize(*args, **kwargs)
            stage["raw_wt"] = {"status": solution.status, "objective_value": solution.objective_value}
            emit()
            if solution.status == "optimal" and finite_nonnegative(solution.objective_value) and 0.1 <= float(solution.objective_value) <= 2.0:
                checks = {"historical_record": math.isclose(float(solution.objective_value), float(old_summary["wt_growth"]), abs_tol=1e-8, rel_tol=1e-6)}
                if key == "canonical_model":
                    checks["first_stage_replay"] = math.isclose(float(solution.objective_value), report["stages"][0]["legacy"]["wt_growth"], abs_tol=1e-8, rel_tol=1e-6)
                stage["early_wt_comparisons"] = checks
                emit()
                if not all(checks.values()):
                    stage["legacy_unavailable_reason"] = "Stopped after WT comparison, before original function requested KO"
                    require(False, "WT comparison failed; no KO dispatched")
            return solution

        def observe_deletion(*args, **kwargs):
            require(len(args) == 1 and set(kwargs) == {"processes", "gene_list"} and
                    kwargs["processes"] == 1 and list(kwargs["gene_list"]) == sorted(selected), "KO request exceeded declared scope")
            frames = []
            stage["raw_deletion"] = {"dispatch": "predeclared single-KO calls to unchanged engine", "calls": []}
            for gene_id in sorted(selected):
                require(model.solver.configuration.presolve != "auto", "Unexpected auto presolve before KO")
                report["primary_solves_requested_upper_bound"] += 1
                require(report["primary_solves_requested_upper_bound"] <= 14, "Primary solve budget exceeded")
                stage["requested_ko"] = gene_id
                emit()
                raw = original_deletion(args[0], gene_list=[gene_id], processes=1)
                stage["raw_deletion"]["calls"].append({"requested_gene_id": gene_id,
                    "columns": list(raw.columns), "rows": [
                        {"index": index, "values": row.to_dict()} for index, row in raw.iterrows()]})
                emit()
                if not {"ids", "growth", "status"} <= set(raw.columns) or raw.empty:
                    stage["legacy_unavailable_reason"] = "Missing raw columns or empty return; do not fabricate normalized rows"
                    require(False, "Raw KO malformed or missing; no further KO dispatched")
                frames.append(raw.copy(deep=True))
                valid_id = (len(raw) == 1 and isinstance(raw.iloc[0]["ids"], (set, frozenset)) and raw.iloc[0]["ids"] == {gene_id})
                if not valid_id or not all(raw["status"].eq("optimal")) or not all(finite_nonnegative(v) for v in raw["growth"]):
                    stage["early_raw_gate_failure"] = "Raw KO wrong/duplicate ID, nonoptimal, None/NaN/Inf or negative growth"
                    break
                growth = float(raw.iloc[0]["growth"])
                ratio = growth / float(stage["raw_wt"]["objective_value"])
                observed_calls = [ratio < cutoff for cutoff in cuts]
                references = [("historical_record", reference[gene_id], historical_calls[gene_id])]
                if key == "canonical_model":
                    first = report["stages"][0]
                    first_row = next(r for r in first["legacy"]["predictions"] if r["gene_id"] == gene_id)
                    first_calls = first["comparisons"]["historical_record"]["genes"][gene_id]["calls"]
                    references.append(("first_stage_replay", first_row, first_calls))
                early_checks = {}
                for label, target, target_calls in references:
                    early_checks[label] = {"growth": math.isclose(growth, float(target["ko_growth"]), abs_tol=1e-8, rel_tol=1e-6),
                                           "ratio": math.isclose(ratio, float(target["ko_growth_ratio"]), abs_tol=1e-8, rel_tol=1e-6),
                                           "calls": observed_calls, "reference_calls": target_calls}
                stage.setdefault("early_ko_comparisons", {})[gene_id] = early_checks
                emit()
                if not all(c["growth"] and c["ratio"] and c["calls"] == c["reference_calls"] for c in early_checks.values()):
                    stage["early_comparison_gate_failure"] = "KO numerical or cutoff mismatch; do not dispatch later KO"
                    break
            combined = validator.pd.concat(frames, ignore_index=True)
            raw_tables.append(combined)
            return combined

        with patch.object(model, "optimize", observe_wt), patch.object(validator, "single_gene_deletion", observe_deletion):
            predictions, wt = validator.run_single_gene_deletions(model, "gurobi", excluded_gene_ids=all_ids - set(selected))
        rows = predictions.to_dict(orient="records")
        stage["legacy"] = {"wt_growth": wt, "predictions": rows}
        emit()
        require("early_raw_gate_failure" not in stage, "Raw KO anomaly: partial legacy retained, later KO not dispatched")
        require("early_comparison_gate_failure" not in stage, "KO comparison failed: partial legacy retained, later KO not dispatched")
        require(len(raw_tables) == 1, "Unexpected deletion call count")
        raw = raw_tables[0]
        require({"ids", "growth", "status"} <= set(raw.columns), "Raw KO columns missing")
        raw_ids = [validator._gene_id_from_deletion_row(index, row) for index, row in raw.iterrows()]
        require(len(raw) == 6 and len(set(raw_ids)) == 6 and set(raw_ids) == set(selected), "Raw KO rows missing, duplicated or unexpected")
        require(all(isinstance(ids, (set, frozenset)) and len(ids) == 1 for ids in raw["ids"]), "Raw KO IDs are not singleton deletions")
        require(all(status == "optimal" for status in raw["status"]), "Raw KO nonoptimal: retain legacy result and stop")
        require(all(finite_nonnegative(value) for value in raw["growth"]), "Raw KO None/NaN/Inf/negative: retain legacy result and stop")
        require(len(rows) == 6 and {r["gene_id"] for r in rows} == set(selected) and
                all(finite_nonnegative(r["ko_growth"]) and finite_nonnegative(r["ko_growth_ratio"]) for r in rows), "Invalid legacy result shape or numbers")
        compare(stage, wt, rows, old_summary["wt_growth"], reference, historical_calls, "historical_record")
        if key == "canonical_model":
            first = report["stages"][0]["legacy"]
            first_calls = {gene_id: report["stages"][0]["comparisons"]["historical_record"]["genes"][gene_id]["calls"] for gene_id in selected}
            compare(stage, wt, rows, first["wt_growth"], {r["gene_id"]: r for r in first["predictions"]}, first_calls, "first_stage_replay")
        stage["status"] = "passed"
        emit()
    report["status"] = "passed"
    report["phase"] = "complete"
    report["claim_limit"] = "Only these six controls and WT; no global recall, calibration or complete model equivalence claim"
except Exception as exc:
    report["status"] = "stopped"
    report.setdefault("failure_reason", "Execution error; no automatic retry or substitution")
    report["exception_type"] = type(exc).__name__
    if report["stages"] and report["stages"][-1]["status"] == "running":
        report["stages"][-1]["status"] = "stopped"
emit()

'''

started_at_utc = datetime.now(timezone.utc).isoformat()
started_monotonic = time.monotonic()
process = subprocess.Popen([sys.executable, "-B", "-c", worker], stdout=subprocess.PIPE,
                           stderr=subprocess.DEVNULL, text=True)
timed_out = False
try:
    output, _ = process.communicate(timeout=600)
except subprocess.TimeoutExpired:
    timed_out = True
    process.kill()
    output, _ = process.communicate()
checkpoints = []
for line in output.splitlines():
    try:
        checkpoints.append(json.loads(line))
    except json.JSONDecodeError:
        pass
result = checkpoints[-1] if checkpoints else {"status": "stopped", "stages": [], "failure_reason": "No complete worker checkpoint"}
if timed_out or process.returncode != 0 or result.get("status") == "running":
    result["status"] = "stopped"
    result["failure_reason"] = "600-second wall-clock limit" if timed_out else "Worker did not complete normally"
    for stage in result.get("stages", []):
        if stage.get("status") == "running":
            stage["status"] = "stopped"
result["execution_started_at_utc"] = started_at_utc
result["execution_finished_at_utc"] = datetime.now(timezone.utc).isoformat()
result["elapsed_wall_seconds"] = time.monotonic() - started_monotonic
result["worker_returncode"] = process.returncode
result["timed_out"] = timed_out
print(json.dumps(result, ensure_ascii=False, allow_nan=False))
sys.exit(0 if result["status"] == "passed" else 1)
