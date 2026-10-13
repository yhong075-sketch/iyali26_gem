"""Net cytosol -> matrix malate import range at the growth optimum (8 LPs; budget top-up of run_task_b.py).

For WT and complex I KO of both models (Task A medium + PO1f profile), fix growth at its optimum and
minimize / maximize v(R807) + v(R902) + v(R1386) - v(R831), the signed matrix-malate gain over every
cytosol<->matrix carrier. Run from the repository root:
    .venv/bin/python -B model/reports/r132_mito_adh_20261012/run_net_malate_fva.py
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_task_b import MODELS, context, sha  # noqa: E402

CARRIERS = {"R807": 1, "R831": -1, "R902": 1, "R1386": 1}


def main():
    for path, digest in MODELS.values():
        assert sha(path) == digest
    out, lps = {}, 0
    for key in MODELS:
        for label, ko in (("WT", None), ("complex_I_KO", "R2062")):
            sim, model = context(key, "glucose")
            with model:
                if ko:
                    model.reactions.get_by_id(ko).bounds = (0.0, 0.0)
                mu = model.slim_optimize(); lps += 1
                model.reactions.biomass_C.lower_bound = mu * (1 - 1e-9)
                uptake = -model.reactions.R1070.lower_bound
                model.objective = model.problem.Objective(
                    sum(c * model.reactions.get_by_id(r).flux_expression for r, c in CARRIERS.items()))
                rng = []
                for direction in ("min", "max"):
                    model.objective.direction = direction
                    rng.append(model.slim_optimize()); lps += 1
                out[f"{key}|{label}"] = {"growth": mu, "net_malate_import": rng,
                                         "per_100_glucose": [100 * v / uptake for v in rng],
                                         "note": "uptake taken as the medium bound 2.414; pFBA vectors show it at the bound"}
    (HERE / "net_malate_fva.json").open("x").write(json.dumps({"lps": lps, "results": out}, indent=1) + "\n")
    print(json.dumps(out, indent=1), "LPs", lps)


if __name__ == "__main__":
    main()
