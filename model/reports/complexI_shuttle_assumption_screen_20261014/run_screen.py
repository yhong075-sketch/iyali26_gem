"""Benchmark side effects of Task D options (a) R533 rev closed and (b) R831 fwd closed. Diagnostic only.

Reuses Task B's screen() and metrics() unchanged (run_task_b.py imported by path); only context() is wrapped to
apply the context (bench: as loaded; diag: + Task D holistic overlay + 13C caps) and the variant closure.

Run from the repository root:
    .venv/bin/python -B model/reports/complexI_shuttle_assumption_screen_20261014/run_screen.py
Writes run.json and call_changes.csv here; per-gene tables and pFBA vectors go to
$IYALI26_RESEARCH_ROOT/artifacts/tasks/complexI_shuttle_assumption_screen_20261014/. Never overwrites.
"""
import csv
import hashlib
import importlib.util
import json
import math
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
TASK_B = REPO / "model/reports/r132_mito_adh_20261012/run_task_b.py"
TASK_D = REPO / "model/reports/complexI_transhydrogenation_20261012/common.py"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


tb = load_module("run_task_b", TASK_B)
td = load_module("task_d_common", TASK_D)

OUT = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve() / "artifacts/tasks" / HERE.name
CONTEXTS = ("bench", "diag")
VARIANTS = {"base": {}, "a": {"R533": "rev"}, "b": {"R831": "fwd"}}
CHECKS = {("bench", "base", "scores15"): dict(TP=82, FP=33, FN=229, TN=697, unresolved=16),
          ("bench", "base", "WT"): 0.3124267,
          ("diag", "base", "complex_I"): 0.826, ("diag", "a", "complex_I"): 0.077,
          ("diag", "b", "complex_I"): 0.077}
CODE = {"run_screen.py": Path(__file__), "run_task_b.py": TASK_B, "task_d_common.py": TASK_D,
        "holistic_overlay.py": TASK_D.parent / "holistic_overlay.py"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


ORIGINAL_CONTEXT = tb.context


def wrapped_context(ctx, variant, applied):
    def context(model_key, medium_key):
        sim, model = ORIGINAL_CONTEXT(model_key, medium_key)
        assert math.isclose(-model.reactions.R1070.lower_bound, td.UPTAKE)
        if ctx == "diag":
            applied["overlay"] = td.apply_overlay(model)
            applied["caps"] = td.add_caps(model, json.loads(td.PINS["caps"][0].read_text()))
        before = {rid: list(model.reactions.get_by_id(rid).bounds) for rid in VARIANTS[variant]}
        td.close(model, VARIANTS[variant])
        applied["closures"] = {rid: {"how": how, "before": before[rid],
                                     "after": list(model.reactions.get_by_id(rid).bounds)}
                               for rid, how in VARIANTS[variant].items()}
        return sim, model
    return context


def main():
    pins = td.check_pins()
    for key, (path, digest) in {"candidate": tb.MODELS["candidate"], "glucose": tb.MEDIA["glucose"],
                                **{k: tb.INPUTS[k] for k in ("profile", "experimental", "mtdna_scope")}}.items():
        assert sha(path) == digest, f"Input identity changed: {key}"
    OUT.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    positives, _ = tb.load_experimental_reference(tb.INPUTS["experimental"][0])
    positives = set(positives)
    mtdna = json.loads(tb.INPUTS["mtdna_scope"][0].read_text())["genes"]
    record = {"started_utc": datetime.now(timezone.utc).isoformat(),
              "command": "python -B model/reports/complexI_shuttle_assumption_screen_20261014/run_screen.py",
              "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
              "git_status": subprocess.check_output(["git", "status", "--short"], cwd=REPO, text=True),
              "task_d_pins": pins, "code": {k: sha(p) for k, p in CODE.items()},
              "inputs": {k: {"path": str(p), "sha256": sha(p)} for k, (p, _) in
                         {"candidate": tb.MODELS["candidate"], "glucose": tb.MEDIA["glucose"],
                          **{k: tb.INPUTS[k] for k in ("profile", "experimental", "mtdna_scope")}}.items()},
              "solver": tb.PARAMS, "versions": {"python": platform.python_version(), "cobra": tb.cobra.__version__,
                                                "gurobi": list(tb.gurobipy.gurobi.version())},
              "scoring": "benchmark.md: mtDNA genes out of scope; negatives are model genes absent from the positive "
                         "list (user-defined class); no solution = unresolved; second view counts no solution as "
                         "lethal (labelled); development reference, not independent validation",
              "runs": {}}
    screens = {}
    for ctx in CONTEXTS:
        for variant in VARIANTS:
            applied = {}
            tb.context = wrapped_context(ctx, variant, applied)
            result, rows = tb.screen("candidate", positives, mtdna)
            label = f"{ctx}_{variant}"
            screens[label] = rows
            run = {k: v for k, v in result.items() if k != "route"}
            run["applied"] = applied
            run["route"] = {k: {"growth": r["growth"], "status": r["status"], "fva": r["fva"]}
                            for k, r in result["route"].items()}
            record["runs"][label] = run
            (OUT / f"pfba_route_{label}.json").open("x").write(
                json.dumps({k: r["pfba"] for k, r in result["route"].items()}) + "\n")
            with (OUT / f"gene_screen_{label}.tsv").open("x") as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            primary = next(s for s in run["scores"] if s["cutoff"] == 0.15 and s["view"] == "benchmark rule")
            print(label, "WT", run["WT"], {k: v["ratio"] for k, v in run["reaction_ko"].items()},
                  {k: primary[k] for k in ("TP", "FP", "FN", "TN", "unresolved", "mcc")}, flush=True)

    # Reproduction checks (TASK.md); recorded, then the script stops before writing interpretation tables.
    checks = []
    for (ctx, variant, what), expected in CHECKS.items():
        run = record["runs"][f"{ctx}_{variant}"]
        if what == "scores15":
            got = {k: next(s for s in run["scores"] if s["cutoff"] == 0.15 and s["view"] == "benchmark rule")[k]
                   for k in expected}
            ok = got == expected
        elif what == "WT":
            got = run["WT"]
            ok = abs(got / expected - 1) <= 0.02
        else:
            got = run["reaction_ko"][what]["ratio"]
            ok = got is not None and (abs(got - expected) <= 0.02 * expected if expected else got == 0)
        checks.append({"context": ctx, "variant": variant, "what": what, "expected": expected, "got": got, "ok": ok})
    record["reproduction_checks"] = checks

    # Gene-level differences vs base, within each context.
    changes = []
    for ctx in CONTEXTS:
        base = {r["gene_id"]: r for r in screens[f"{ctx}_base"]}
        for variant in ("a", "b"):
            for r in screens[f"{ctx}_{variant}"]:
                b = base[r["gene_id"]]
                if (b["status"] == r["status"] and (b["ratio"] is None) == (r["ratio"] is None)
                        and (r["ratio"] is None or abs(r["ratio"] - b["ratio"]) <= 1e-6)):
                    continue
                row = {"context": ctx, "variant": variant, "gene_id": r["gene_id"], "positive": r["positive"],
                       "mtdna": r["mtdna"], "reactions": r["reactions"], "base_status": b["status"],
                       "variant_status": r["status"], "base_ratio": b["ratio"], "variant_ratio": r["ratio"]}
                for c in tb.CUTOFFS:
                    call = lambda x: "unresolved" if x is None else ("E" if x < c else "NE")  # noqa: E731
                    row[f"call_change_{c:g}"] = (f"{call(b['ratio'])}->{call(r['ratio'])}"
                                                 if call(b["ratio"]) != call(r["ratio"]) else "")
                changes.append(row)
    record["gene_differences"] = changes
    if changes:
        with (HERE / "call_changes.csv").open("x") as f:
            writer = csv.DictWriter(f, fieldnames=list(changes[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(changes)
    record.update(lps=tb.Counter_.lps, elapsed_seconds=time.monotonic() - started,
                  inputs_unchanged=all(sha(v["path"]) == v["sha256"] for v in record["inputs"].values()),
                  code_unchanged=all(sha(p) == record["code"][k] for k, p in CODE.items()),
                  outputs={p.name: sha(p) for p in sorted(OUT.iterdir())})
    with (HERE / "run.json").open("x") as f:
        f.write(json.dumps(record, indent=1, default=str) + "\n")
    for c in checks:
        print("check", c, flush=True)
    print("gene differences", len(changes), "LPs", tb.Counter_.lps, "elapsed", round(record["elapsed_seconds"]), "s")


if __name__ == "__main__":
    main()
