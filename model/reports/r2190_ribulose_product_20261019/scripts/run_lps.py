"""The 8 authorized LPs for R2190-RIBULOSE-20261019 (TASK.md): WT growth on M0 and M1, then on M1
min/max of R1964, R2041 and R2042 with every boundary reaction closed and forced fluxes relaxed to
include 0 (protocol of orphan_reactions_20261019/scripts/loop_check.py, as six explicit LPs).

Context: model/reports/coq9_respiratory_package_20261008/run_controls.py (PO1f SD-Leu, Gurobi config
of coq_candidate_checks_20261005/manifest.json). Run from the repository root:
python -B model/reports/r2190_ribulose_product_20261019/scripts/run_lps.py
"""
import copy
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "platform"))

from tools.diagnose_closed_energy import configure_solver
from tools.validate_energy_candidates import load_model

HERE = Path(__file__).resolve().parents[1]
MODELS = {
    "M0": (REPO / "model/candidates/E5_coq9_resp_pkg_r132_20261012.xml",
           "a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515"),
    "M1": (REPO / "model/candidates/E5_coq9_resp_pkg_r132_r2190_ribulose_20261019.xml", sys.argv[1]),
}
LOOP = ("R1964", "R2041", "R2042")
SHOWN = ("biomass_C", "R2190", *LOOP, "R803", "R714")
BUDGET = 8


def main():
    out = HERE / "lp_results.json"
    assert not out.exists(), "Never overwrite an existing LP record"
    config = copy.deepcopy(json.loads(
        (REPO / "model/reports/coq_candidate_checks_20261005/manifest.json").read_text())["config"])
    results, calls, sims = {}, 0, {}

    def solve(name, model):
        nonlocal calls
        calls += 1
        assert calls <= BUDGET, "Declared LP budget exceeded"
        configure_solver(model, config)
        solution = model.optimize()
        flux = solution.fluxes if solution.status == "optimal" else None
        results[name] = {"status": solution.status, "objective": solution.objective_value,
                         "fluxes": None if flux is None else {rid: float(flux[rid]) for rid in SHOWN}}

    for label, (path, digest) in MODELS.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, f"{label} identity changed"
        sims[label] = load_model(path, config)
        with sims[label].model as model:
            solve(f"WT_{label}", model)

    closure = {}
    with sims["M1"].model as model:
        for reaction in model.boundary:
            reaction.bounds = (0, 0)
        forced = sorted(r.id for r in model.reactions if r.lower_bound > 0 or r.upper_bound < 0)
        for rid in forced:
            r = model.reactions.get_by_id(rid)
            r.bounds = (min(r.lower_bound, 0), max(r.upper_bound, 0))
        closure = {"boundary_reactions_closed": len(model.boundary), "forced_relaxed": forced}
        for rid in LOOP:
            for direction in ("min", "max"):
                model.objective = model.reactions.get_by_id(rid)
                model.objective_direction = direction
                solve(f"closed_{rid}_{direction}_M1", model)

    wt0, wt1 = results["WT_M0"]["objective"], results["WT_M1"]["objective"]
    checks = {
        "wt_relative_difference": abs(wt1 - wt0) / abs(wt0) if wt0 else None,
        "wt_equal_within_1e-9": wt0 is not None and wt1 is not None and abs(wt1 - wt0) <= 1e-9 * abs(wt0),
        "closed_ranges_zero_within_1e-7": all(
            results[f"closed_{rid}_{d}_M1"]["status"] == "optimal"
            and abs(results[f"closed_{rid}_{d}_M1"]["objective"]) <= 1e-7 for rid in LOOP for d in ("min", "max")),
    }
    record = {"models": {k: {"path": str(p.relative_to(REPO)), "sha256": d} for k, (p, d) in MODELS.items()},
              "simulation_context": {k: s.provenance() for k, s in sims.items()}, "solver": config["solver"],
              "closure_M1": closure, "results": results, "checks": checks, "optimization_calls": calls}
    out.write_text(json.dumps(record, indent=2, default=str) + "\n")
    print(json.dumps({k: {x: v[x] for x in ("status", "objective")} for k, v in results.items()}, indent=1))
    print(json.dumps(checks, indent=1), "calls", calls)


if __name__ == "__main__":
    main()
