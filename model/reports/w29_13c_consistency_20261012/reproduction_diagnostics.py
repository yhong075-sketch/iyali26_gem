"""Task C: why the recalled V1 elastic-fit totals (6.94 / 1.14 / 2.72 / 0.68) are not reproduced.

Diagnostic only (TASK.md: "a diagnostic for the cause, not a substitute for the stop rule").
Each variant changes one thing relative to the pre-registered fit and reports the four V1 totals
at f = 1 (no FVA; one optimum per setup is enough to compare totals, which are unique).

Variants:
- pre_registered: TASK.md mapping, Task B candidate (R132 removed) + overlay;
- gcv_R361: row 21 mapped to -R361 (the earlier, mis-targeted glycine-cleavage reaction);
- idh_sum_ci: row 9 interval = [LB23 + LB24, UB23 + UB24] instead of R21 net;
- pc_me_separate: row 22 replaced by two rows, R31 (PC R690) and R29 net (-R538);
- parent_R132: parent model 0ac2addf (R132 present) + overlay, pre-registered mapping;
- caps_*: Task B 13C caps (PC, net malate import, glycine-cleavage form) added as hard constraints.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/w29_13c_consistency_20261012/reproduction_diagnostics.py
Writes reproduction_diagnostics.json and .csv here. Never overwrites.
"""
import csv
import json
import time
from datetime import datetime, timezone

import common
from common import HERE, MU, REPO, PINS, Counter_, add_caps, apply_overlay, check_pins, growth, load, load_caps, mfa_rows, sha
from elastic_fit import MAPPING, RECALLED_V1, SETUPS, ranges

PARENT = (REPO / "model/candidates/E5_coq9_respiratory_package_20261008.xml",
          "0ac2addfa0cb6de2d89393e433766fb04325ffe42c2d131eaefaca3bcec29dba")


def total(model, rng, mapping, keys, mu):
    per100 = 100.0 / common.UPTAKE
    prob = model.problem
    model.reactions.biomass_C.bounds = (mu, mu)
    extra, terms = [], []
    for key, _, _, coeffs in mapping:
        if key not in keys:
            continue
        lb, ub, _ = rng[key]
        e = per100 * sum(c * model.reactions.get_by_id(r).flux_expression for r, c in coeffs.items())
        sp, sm = prob.Variable(f"sp_{key}", lb=0), prob.Variable(f"sm_{key}", lb=0)
        extra += [sp, sm, prob.Constraint(e - sp, ub=ub), prob.Constraint(e + sm, lb=lb)]
        terms.append((sp + sm) / (ub - lb))
    model.add_cons_vars(extra)
    model.objective = prob.Objective(sum(terms), direction="min")
    Counter_.lps += 1
    model.slim_optimize()
    if model.solver.status != "optimal":
        return None, {}
    viol = {}
    for key in keys:
        sp = model.solver.variables[f"sp_{key}"].primal
        sm = model.solver.variables[f"sm_{key}"].primal
        if sp > 1e-6 or sm > 1e-6:
            viol[key] = round(sp - sm, 2)
    return model.solver.objective.value, viol


def variants():
    base_rng = ranges()
    rows = mfa_rows()
    out = {"pre_registered": (MAPPING, base_rng)}
    m = [(k, r, n, ({"R361": -1} if k == "GCV" else c)) for k, r, n, c in MAPPING]
    out["gcv_R361"] = (m, base_rng)
    rng = dict(base_rng)
    rng["IDH_total"] = (float(rows["R23"]["glucose_LB"]) + float(rows["R24"]["glucose_LB"]),
                        float(rows["R23"]["glucose_UB"]) + float(rows["R24"]["glucose_UB"]), None)
    out["idh_sum_ci"] = (MAPPING, rng)
    m = [x for x in MAPPING if x[0] != "PC_minus_ME"] + [("PC", "R31", "pyruvate carboxylase", {"R690": 1}),
                                                         ("ME_net", "R29 net", "mito malic enzyme (MFA sign)", {"R538": 1})]
    rng = dict(base_rng)
    rng.pop("PC_minus_ME")
    rng["PC"] = (float(rows["R31"]["glucose_LB"]), float(rows["R31"]["glucose_UB"]), None)
    rng["ME_net"] = (float(rows["R29 net"]["glucose_LB"]), float(rows["R29 net"]["glucose_UB"]), None)
    out["pc_me_separate"] = (m, rng)
    return out


def main():
    inputs = check_pins()
    assert sha(PARENT[0]) == PARENT[1], "parent model changed"
    caps = load_caps()
    started = time.monotonic()
    results = []
    sim, model = load("V1")
    for vname, (mapping, rng) in variants().items():
        for setup, rev, with_malate in SETUPS:
            keys = [k for k, *_ in mapping if with_malate or k != "malate_import"]
            with model:
                apply_overlay(model)
                if rev:
                    common.r825_reversible(model)
                with model:
                    t, v = total(model, rng, mapping, keys, MU)
            results.append({"variant": vname, "setup": setup, "total": t, "violations": v})
            print(vname, setup, t, v, flush=True)
    for gcv in ("R361_looser", "module"):
        for setup, rev, with_malate in SETUPS:
            keys = [k for k, *_ in MAPPING if with_malate or k != "malate_import"]
            with model:
                apply_overlay(model)
                if rev:
                    common.r825_reversible(model)
                add_caps(model, caps, gcv)
                g, _ = growth(model)
                with model:
                    t, v = total(model, ranges(), MAPPING, keys, MU) if g and g >= MU else (None, {"infeasible_at_mu": g})
            results.append({"variant": f"caps_{gcv}", "setup": setup, "total": t, "violations": v, "wt_max": g})
            print(f"caps_{gcv}", setup, g, t, v, flush=True)
    # Parent model (R132 present), V1, pre-registered mapping.
    original = PINS["model"]
    PINS["model"] = PARENT
    try:
        sim_p, parent = load("V1")
    finally:
        PINS["model"] = original
    for setup, rev, with_malate in SETUPS:
        keys = [k for k, *_ in MAPPING if with_malate or k != "malate_import"]
        with parent:
            apply_overlay(parent)
            if rev:
                common.r825_reversible(parent)
            with parent:
                t, v = total(parent, ranges(), MAPPING, keys, MU)
        results.append({"variant": "parent_R132", "setup": setup, "total": t, "violations": v})
        print("parent_R132", setup, t, v, flush=True)
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs,
              "parent_model": {"path": str(PARENT[0]), "sha256": PARENT[1]}, "medium": "V1",
              "recalled_V1": RECALLED_V1, "results": results, "lps": Counter_.lps,
              "elapsed_seconds": time.monotonic() - started}
    (HERE / "reproduction_diagnostics.json").open("x").write(json.dumps(record, indent=1, default=str) + "\n")
    with (HERE / "reproduction_diagnostics.csv").open("x", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["variant", "setup", "recalled_V1_total", "total", "relative_difference", "violations_per100"])
        for r in results:
            exp = RECALLED_V1[r["setup"]]
            rel = None if r["total"] is None else round(r["total"] / exp - 1, 4)
            w.writerow([r["variant"], r["setup"], exp, None if r["total"] is None else round(r["total"], 4), rel,
                        json.dumps(r["violations"])])


if __name__ == "__main__":
    main()
