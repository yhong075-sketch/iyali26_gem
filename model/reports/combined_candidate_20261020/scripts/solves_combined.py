"""The 20 authorized LPs for the combined candidate (GLPK, SD-Leu + PO1f).

Usage: python solves_combined.py MEDIUM.csv PROFILE.json BASELINE.xml COMBINED.xml OUT.json
- WT growth on the baseline and the combined model (2 LPs).
- Arabinose-loop test on the combined model: every boundary reaction closed, forced fluxes relaxed to
  include 0, FVA (min and max) on R1964, R2041, R2042 and R2190 (8 LPs). On a56f3649 the same test
  gave R1964/R2041/R2042 = 1000 (model/reports/orphan_reactions_20261019/loop_check.txt).
- Single-gene KO on the combined model for the 10 genes whose reaction links changed (10 LPs).
A non-optimal status is recorded as such, not as growth 0.
"""
import json
import os
import sys

os.environ.setdefault("HOME", os.path.abspath("cache"))
import pandas as pd
from cobra.io import read_sbml_model

LOOP = ["R1964", "R2041", "R2042", "R2190"]
KO_GENES = ["YALI1E07744g", "YALI1D02625g", "YALI1C33005g", "YALI1D09420g", "YALI1E11768g",
            "YALI1E11370g", "YALI1D17462g", "YALI1F03452g", "YALI1E06616g", "YALI1A21372g"]


def load(path, medium_csv, profile_json):
    model = read_sbml_model(path)
    model.solver = "glpk"
    model.solver.configuration.presolve = True
    med = pd.read_csv(medium_csv, dtype={"exchange": str})
    medium = {r.exchange.strip(): float(r.uptake) for r in med.itertuples() if float(r.uptake) > 0}
    prof = json.load(open(profile_json))
    medium.update({k: float(v) for k, v in prof.get("medium", {}).get("runtime_uptake_overrides", {}).items()})
    model.medium = {k: v for k, v in medium.items() if k in model.reactions}
    for op in prof.get("operations", []):
        if op.get("type") == "disable_reaction":
            model.reactions.get_by_id(op["reaction_id"]).bounds = (op["set_bounds"]["lower_bound"],
                                                                   op["set_bounds"]["upper_bound"])
    return model


def solve(model):
    s = model.optimize()
    return {"status": s.status, "value": s.objective_value if s.status == "optimal" else None}


if __name__ == "__main__":
    medium_csv, profile_json, base_path, combined_path, out = sys.argv[1:6]
    result = {"lp_count": 0}
    base = load(base_path, medium_csv, profile_json)
    result["wt_baseline"] = solve(base); result["lp_count"] += 1
    model = load(combined_path, medium_csv, profile_json)
    wt = solve(model); result["wt_combined"] = wt; result["lp_count"] += 1
    loop = {}
    with model:
        for b in model.boundary:
            b.bounds = (0, 0)
        for r in model.reactions:
            if r.lower_bound > 0 or r.upper_bound < 0:
                r.bounds = (min(r.lower_bound, 0), max(r.upper_bound, 0))
        for rid in LOOP:
            rxn = model.reactions.get_by_id(rid)
            for sense in ("min", "max"):
                model.objective = rxn
                model.objective_direction = sense
                loop[f"{rid}_{sense}"] = solve(model); result["lp_count"] += 1
    result["loop_all_boundaries_closed"] = loop
    kos = {}
    for gid in KO_GENES:
        with model:
            gene = model.genes.get_by_id(gid)
            gene.knock_out()
            off = sorted(r.id for r in gene.reactions if r.bounds == (0, 0))
            s = solve(model); result["lp_count"] += 1
            ratio = (s["value"] / wt["value"]) if (s["value"] is not None and wt["value"]) else None
            kos[gid] = {**s, "ko_wt_ratio": ratio, "reactions": sorted(r.id for r in gene.reactions),
                        "reactions_off": off}
    result["gene_ko_combined"] = kos
    json.dump(result, open(out, "w"), indent=1)
    print(json.dumps(result, indent=1))
