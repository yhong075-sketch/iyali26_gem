"""Static checks plus exactly four control LPs on the package candidate (protocol of
reports/coq_r305_candidate_20261005/run.py: WT, R385 closed, closed-input ATP in C_cy and C_mi).

Run from the repository root: python -B model/reports/coq9_respiratory_package_20261008/run_controls.py
"""
import copy
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "platform"))

from cobra import Reaction
from scripts.gem_annotate.coq9 import exact_residual
from tools.validate_energy_candidates import load_model, close_model, balance
from tools.diagnose_closed_energy import configure_solver

HERE = Path(__file__).resolve().parent
MODEL_PATH = REPO / "model/candidates/E5_coq9_respiratory_package_20261008.xml"
MODEL_SHA = "0ac2addfa0cb6de2d89393e433766fb04325ffe42c2d131eaefaca3bcec29dba"
TOUCHED = ("R39", "R695", "R2062", "R304")


def main():
    out = HERE / "controls.json"
    assert not out.exists(), "Never overwrite an existing control record"
    assert hashlib.sha256(MODEL_PATH.read_bytes()).hexdigest() == MODEL_SHA
    config = copy.deepcopy(json.loads(
        (REPO / "model/reports/coq_candidate_checks_20261005/manifest.json").read_text())["config"])
    sim = load_model(MODEL_PATH, config)
    model = sim.model
    static = {rid: {"equation": model.reactions.get_by_id(rid).reaction,
                    "compartments": sorted(model.reactions.get_by_id(rid).compartments),
                    "gpr": model.reactions.get_by_id(rid).gene_reaction_rule,
                    "exact_residual": exact_residual(model.reactions.get_by_id(rid))}
              for rid in TOUCHED}
    assert all(not row["exact_residual"] for row in static.values())
    assert static["R39"]["compartments"] == ["C_mi"]
    results, calls = {}, 0

    def solve(name, m):
        nonlocal calls
        calls += 1
        assert calls <= 4, "Declared control budget exceeded"
        configure_solver(m, config)
        solution = m.optimize()
        flux = solution.fluxes if solution.status == "optimal" else None
        results[name] = {"status": solution.status, "objective": solution.objective_value,
                         "growth": None if flux is None else float(flux["biomass_C"]),
                         "fluxes": None if flux is None else {
                             rid: float(flux[rid]) for rid in ("R385", "R2062", "R304", "R305", "R39", "R171")}}

    with model:
        solve("WT", model)
    with model:
        model.reactions.R385.bounds = (0, 0)
        solve("R385_off", model)
    closed, changes = close_model(model, config)
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
            solve("ATP_closed_" + compartment, closed)
    record = {"model": str(MODEL_PATH.relative_to(REPO)), "model_sha256": MODEL_SHA,
              "simulation_context": sim.provenance(), "solver": config["solver"],
              "closure_changes": len(changes), "static_checks": static, "results": results,
              "optimization_calls": calls}
    out.write_text(json.dumps(record, indent=2, default=str) + "\n")
    print(json.dumps({k: {x: v[x] for x in ("status", "objective")} for k, v in results.items()}, indent=1))


if __name__ == "__main__":
    main()
