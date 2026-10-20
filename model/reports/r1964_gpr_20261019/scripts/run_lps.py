"""The 5 authorized LPs for R1964-ARADH-GPR-20261019 (TASK.md): WT on a56f3649 (M0) and the new
candidate (M1), KO of YALI1F03452g on M0 and M1, and KO of YALI1E06616g on M1.

Protocol of model/reports/coq9_respiratory_package_20261008/run_controls.py, as in
model/reports/r2041_forward_20261019 (worktree claude/bold-chaplygin-1c16df). Run from the repository
root: python -B model/reports/r1964_gpr_20261019/scripts/run_lps.py
"""
import copy
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "platform"))

from tools.validate_energy_candidates import load_model
from tools.diagnose_closed_energy import configure_solver

HERE = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPO / "model/reports/coq_candidate_checks_20261005/manifest.json"
M0 = (REPO / "model/candidates/E5_coq9_resp_pkg_r132_20261012.xml",
      "a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515", "")
M1_PATH = REPO / "model/candidates/E5_coq9_resp_pkg_r132_r1964_gpr_20261019.xml"
M1_GPR = "YALI1F03452g or YALI1E06616g"
PLAN = (("M0", None), ("M0", "YALI1F03452g"), ("M1", None), ("M1", "YALI1F03452g"), ("M1", "YALI1E06616g"))
WATCH = ("R1964", "R2190", "R1980", "R801", "R803", "R2041", "R2042")
BUDGET = 5
CUTOFFS = (0.01, 0.05, 0.10, 0.15)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    out = HERE / "lp_results.json"
    assert not out.exists(), "Never overwrite an existing LP record"
    build = json.loads(M1_PATH.with_suffix(".build.json").read_text())
    m1_sha = sha(M1_PATH)
    assert build["output"]["sha256"] == m1_sha, "M1 differs from its build record"
    models = {"M0": M0, "M1": (M1_PATH, m1_sha, M1_GPR)}
    config = copy.deepcopy(json.loads(CONFIG_PATH.read_text())["config"])
    calls = 0
    record = {"config_path": str(CONFIG_PATH.relative_to(REPO)), "config_sha256": sha(CONFIG_PATH),
              "script_sha256": sha(__file__), "solver": config["solver"], "models": {}, "results": {},
              "ko_method": "cobra Gene.knock_out() inside a model context on the loaded simulation context"}
    loaded = {}
    for name, (path, expected_sha, gpr) in models.items():
        assert sha(path) == expected_sha, f"{name} identity changed"
        sim = load_model(path, config)
        assert sim.model.reactions.R1964.gene_reaction_rule == gpr, f"{name} R1964 GPR unexpected"
        loaded[name] = sim.model
        record["models"][name] = {"path": str(path.relative_to(REPO)), "sha256": expected_sha,
                                  "r1964_gpr": gpr, "simulation_context": sim.provenance(),
                                  "objective": str(sim.model.objective.expression)}
    wt = {}
    for name, gene in PLAN:
        model = loaded[name]
        label = f"{name}_WT" if gene is None else f"{name}_KO_{gene}"
        calls += 1
        assert calls <= BUDGET, "Declared LP budget exceeded"
        with model:
            if gene is not None:
                model.genes.get_by_id(gene).knock_out()
                closed = sorted(r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0))
            configure_solver(model, config)
            solution = model.optimize(raise_error=False)
            row = {"status": solution.status, "objective": solution.objective_value}
            if gene is not None:
                row["reactions_closed_by_ko"] = closed
            if solution.status == "optimal":
                row["fluxes"] = {rid: float(solution.fluxes[rid]) for rid in WATCH if rid in model.reactions}
        record["results"][label] = row
        print(label, row["status"], row["objective"], row.get("reactions_closed_by_ko", ""), flush=True)
        if row["status"] != "optimal":
            record["stopped"] = f"{label} not optimal; raw value kept, no retry"
            break
        if gene is None:
            wt[name] = row["objective"]
        else:
            ratio = row["objective"] / wt[name]
            row["ko_wt"] = ratio
            row["essential_at"] = {f"{int(c * 100)}%": ratio < c for c in CUTOFFS}
    if "M0_WT" in record["results"] and "M1_WT" in record["results"]:
        a, b = wt["M0"], wt["M1"]
        record["wt_relative_difference"] = abs(a - b) / abs(a)
    record["optimization_calls"] = calls
    out.write_text(json.dumps(record, indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()
