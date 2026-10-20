"""The 8 authorized LPs: WT growth and three single-gene KOs, on the baseline and the new model.

Usage: python solves.py MEDIUM.csv PROFILE.json OUT.csv BASELINE.xml NEW.xml
Condition as in the 2026-10-19 orphan report: the whole medium is replaced by the SD-Leu uptakes,
the PO1f runtime overrides are applied and its disable_reaction operations close R612. GLPK with
presolve. A non-optimal status is recorded as such, not as growth 0.
"""
import json
import os
import sys

os.environ.setdefault("HOME", os.path.abspath("cache"))
import pandas as pd
from cobra.io import read_sbml_model

GENES = ["YALI1E31685g", "YALI1A21372g", "YALI1A08512g"]


def condition(medium_csv, profile_json):
    med = pd.read_csv(medium_csv, dtype={"exchange": str})
    medium = {r.exchange.strip(): float(r.uptake) for r in med.itertuples() if float(r.uptake) > 0}
    prof = json.load(open(profile_json))
    medium.update({k: float(v) for k, v in prof.get("medium", {}).get("runtime_uptake_overrides", {}).items()})
    closed = {op["reaction_id"]: (op["set_bounds"]["lower_bound"], op["set_bounds"]["upper_bound"])
              for op in prof.get("operations", []) if op.get("type") == "disable_reaction"}
    return medium, closed


def run(path, medium, closed, label):
    model = read_sbml_model(path)
    model.solver = "glpk"
    model.solver.configuration.presolve = True
    model.medium = {k: v for k, v in medium.items() if k in model.reactions}
    for rid, bounds in closed.items():
        model.reactions.get_by_id(rid).bounds = bounds
    rows = []
    sol = model.optimize()
    wt = sol.objective_value if sol.status == "optimal" else None
    rows.append({"model": label, "case": "WT", "status": sol.status, "growth": wt, "ko_wt_ratio": 1.0 if wt else None,
                 "reactions_off": ""})
    for gid in GENES:
        with model:
            gene = model.genes.get_by_id(gid)
            gene.knock_out()
            off = sorted(r.id for r in gene.reactions if r.bounds == (0, 0))
            s = model.optimize()
            g = s.objective_value if s.status == "optimal" else None
            rows.append({"model": label, "case": f"KO {gid}", "status": s.status, "growth": g,
                         "ko_wt_ratio": (g / wt) if (g is not None and wt) else None,
                         "reactions_off": ";".join(off)})
    return rows


if __name__ == "__main__":
    medium_csv, profile_json, out, base, new = sys.argv[1:6]
    medium, closed = condition(medium_csv, profile_json)
    rows = run(base, medium, closed, "baseline") + run(new, medium, closed, "new")
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print("closed by profile:", closed)
    print(df.to_string(index=False))
