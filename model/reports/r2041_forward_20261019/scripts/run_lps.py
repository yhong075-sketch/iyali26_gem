"""The 10 authorized LPs for R2041-FORWARD-20261019 (TASK.md): WT, R2041 range at 100 % WT,
and closed-exchange ATP dissipation in C_cy and C_mi, on a56f3649 (M0) and the new candidate (M1).

Protocol of model/reports/coq9_respiratory_package_20261008/run_controls.py. Run from the repository
root: python -B model/reports/r2041_forward_20261019/scripts/run_lps.py
"""
import copy
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "platform"))

from cobra import Reaction
from tools.validate_energy_candidates import load_model, close_model, balance
from tools.diagnose_closed_energy import configure_solver

HERE = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPO / "model/reports/coq_candidate_checks_20261005/manifest.json"
M0 = (REPO / "model/candidates/E5_coq9_resp_pkg_r132_20261012.xml",
      "a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515", (-1000.0, 1000.0))
M1_PATH = REPO / "model/candidates/E5_coq9_resp_pkg_r132_r2041_forward_20261019.xml"
WATCH = ("R2041", "R2042", "R2190", "R1964", "R803", "R714", "R1063", "R1062", "R1981", "biomass_C")
BUDGET = 10
TOL = 1e-7


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    out = HERE / "lp_results.json"
    assert not out.exists(), "Never overwrite an existing LP record"
    build = json.loads(M1_PATH.with_suffix(".build.json").read_text())
    m1_sha = sha(M1_PATH)
    assert build["output"]["sha256"] == m1_sha, "M1 differs from its build record"
    models = {"M0": M0, "M1": (M1_PATH, m1_sha, (0.0, 1000.0))}
    config = copy.deepcopy(json.loads(CONFIG_PATH.read_text())["config"])
    calls = 0
    record = {"config_path": str(CONFIG_PATH.relative_to(REPO)), "config_sha256": sha(CONFIG_PATH),
              "script_sha256": sha(__file__), "solver": config["solver"], "models": {}, "results": {}}

    def solve(label, model):
        nonlocal calls
        calls += 1
        assert calls <= BUDGET, "Declared LP budget exceeded"
        configure_solver(model, config)
        solution = model.optimize(raise_error=False)
        flux = solution.fluxes if solution.status == "optimal" else None
        row = {"status": solution.status, "objective": solution.objective_value,
               "fluxes": None if flux is None else {rid: float(flux[rid]) for rid in WATCH if rid in model.reactions}}
        if flux is not None:
            row["nonzero_flux_count"] = int((flux.abs() > TOL).sum())
        record["results"][label] = row
        print(label, row["status"], row["objective"], flush=True)
        return row, flux

    for name, (path, expected_sha, r2041_bounds) in models.items():
        assert sha(path) == expected_sha, f"{name} identity changed"
        sim = load_model(path, config)
        model = sim.model
        assert model.reactions.R2041.bounds == r2041_bounds, f"{name} R2041 bounds unexpected"
        record["models"][name] = {"path": str(path.relative_to(REPO)), "sha256": expected_sha,
                                  "r2041_bounds": list(r2041_bounds), "simulation_context": sim.provenance(),
                                  "objective": str(model.objective.expression)}
        with model:
            wt, _ = solve(f"{name}_WT", model)
        if wt["status"] != "optimal":
            raise SystemExit(f"{name} WT not optimal; stop (no retry)")
        objective_reactions = [r for r in model.reactions if r.objective_coefficient != 0]
        assert len(objective_reactions) == 1 and objective_reactions[0].objective_coefficient == 1
        target = objective_reactions[0]
        for sense in ("min", "max"):
            with model:
                # 100 % of the WT optimum, as cobra FVA with fraction_of_optimum = 1.0.
                target.lower_bound = wt["objective"]
                model.objective = model.reactions.R2041
                model.objective.direction = sense
                solve(f"{name}_R2041_{sense}_at_WT", model)
        closed, changes = close_model(model, config)
        record["models"][name]["closure_changes"] = len(changes)
        for compartment in ("C_cy", "C_mi"):
            with closed:
                if compartment == "C_cy":
                    drain = closed.reactions.get_by_id(config["maintenance"])
                else:
                    drain = Reaction("DIAG_ATP_DISSIPATION_mi", lower_bound=0, upper_bound=1000)
                    drain.add_metabolites({closed.metabolites.get_by_id(mid): c for mid, c in
                        {"m46[C_mi]": -1, "m26[C_mi]": -1, "m197[C_mi]": 1, "m58[C_mi]": 1}.items()})
                    closed.add_reactions([drain])
                chem = balance(drain)
                assert chem["element_status"] == "balanced" and chem["charge_status"] == "balanced_as_stored"
                closed.objective = drain
                closed.objective.direction = "max"
                row, flux = solve(f"{name}_ATP_closed_{compartment}", closed)
                if row["status"] == "optimal" and row["objective"] > TOL:
                    # Stop condition: keep the witness from this same solve.
                    row["witness"] = {rid: float(v) for rid, v in flux.items() if abs(v) > TOL}
    record["optimization_calls"] = calls
    out.write_text(json.dumps(record, indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()
