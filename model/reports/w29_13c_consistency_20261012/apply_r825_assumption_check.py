"""Task C, Amendment 2: apply the R825-reversible ASSUMPTION from its condition file and check it
reproduces the in-memory diagnostic (elastic-fit totals V2 / V1 with overlay, 22 ranges; WT growth V2
without overlay). Diagnostic only.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/w29_13c_consistency_20261012/apply_r825_assumption_check.py
Writes apply_r825_assumption_check.json here. Never overwrites.
"""
import json
import math

from common import HERE, REPO, apply_overlay, check_pins, growth, load, sha
from r825_loop_carnitine_check import solve_fit

CONDITION = REPO / "model/conditions/assumptions/r825_reversible_20261013.json"
EXPECTED = {"fit|V2": 1.433071041522241, "fit|V1": 1.2116243648248852, "wt|V2|no_overlay": 0.29570556570980255}


def apply_condition(model, spec):
    assert spec["status"] == "assumption" and spec["evidence_tier"] == "assumption"
    r = model.reactions.get_by_id(spec["reaction"])
    assert list(r.bounds) == spec["original_bounds"], (r.bounds, spec["original_bounds"])
    assert r.gene_reaction_rule == spec["gpr"]
    setattr(r, spec["field"], spec["value"])


def main():
    inputs = check_pins()
    spec = json.loads(CONDITION.read_text())
    got = {}
    for medium in ("V2", "V1"):
        sim, model = load(medium)
        with model:
            apply_overlay(model)
            apply_condition(model, spec)
            got[f"fit|{medium}"] = solve_fit(model)[0]
        if medium == "V2":
            with model:
                apply_condition(model, spec)
                got["wt|V2|no_overlay"] = growth(model)[0]
    checks = {k: {"expected": EXPECTED[k], "observed": v, "match_1e-6": math.isclose(v, EXPECTED[k], rel_tol=1e-6)}
              for k, v in got.items()}
    out = {"condition_file": str(CONDITION.relative_to(REPO)), "condition_sha256": sha(CONDITION), "inputs": inputs,
           "checks": checks, "all_match": all(c["match_1e-6"] for c in checks.values())}
    (HERE / "apply_r825_assumption_check.json").open("x").write(json.dumps(out, indent=1) + "\n")
    print(json.dumps(checks, indent=1), out["all_match"], out["condition_sha256"])


if __name__ == "__main__":
    main()
