"""In-memory test: is mtFAS blocked by the mitochondrial ACP1 pool or by a missing octanoyl-ACP outlet?

Usage: python acp_block_test.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json
SD-Leu + PO1f, GLPK, growth not required. Objective in every LP: maximise R1405 (trans-oct-2-enoyl-ACP
-> octanoyl-ACP, the last mtFAS step before octanoyl-ACP). Each scenario adds test reactions inside a
`with model:` block, so nothing is saved.
  T0 as supplied
  T1 + ACP1[m] source (-> m179)                       existing outlets kept
  T2 + octanoyl-ACP[m] sink (m1421 ->)                an outlet that does not return ACP1
  T3 + octanoyl-ACP[m] -> ACP1[m]                     an outlet that returns ACP1 (figure's LIP2-type test, apo-H omitted)
  T4 + cytosolic ACP -> ACP1[m] return (m1920 -> m179) links the two ACP species, existing outlets kept
"""
import json
import os
import sys

os.environ.setdefault("HOME", os.path.abspath("cache"))
import pandas as pd
from cobra import Reaction
from cobra.io import read_sbml_model

ACP1_M, OCT_ACP_M, ACP_CY = "m179[C_mi]", "m1421[C_mi]", "m1920[C_cy]"
OUTLETS = ["R2218", "R2216", "R2217", "R2221", "R2236", "R2237", "R78", "R539", "R1394", "R1396"]


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


def scenario(model, name, stoich):
    with model:
        if stoich:
            test = Reaction(f"TEST_{name}", lower_bound=0, upper_bound=1000)
            model.add_reactions([test])
            test.add_metabolites({model.metabolites.get_by_id(m): c for m, c in stoich.items()})
        model.objective = "R1405"
        s = model.optimize()
        row = {"scenario": name, "added": stoich, "status": s.status,
               "max_R1405": s.objective_value if s.status == "optimal" else None}
        if s.status == "optimal" and s.objective_value > 1e-9:
            row["fluxes"] = {r: round(s.fluxes[r], 6) for r in OUTLETS + ([f"TEST_{name}"] if stoich else [])}
        return row


if __name__ == "__main__":
    model_path, medium_csv, profile_json, out = sys.argv[1:5]
    model = load(model_path, medium_csv, profile_json)
    rows = [scenario(model, "T0_as_supplied", {}),
            scenario(model, "T1_ACP1_source", {ACP1_M: 1}),
            scenario(model, "T2_octanoylACP_sink", {OCT_ACP_M: -1}),
            scenario(model, "T3_octanoylACP_to_ACP1", {OCT_ACP_M: -1, ACP1_M: 1}),
            scenario(model, "T4_cytACP_to_ACP1", {ACP_CY: -1, ACP1_M: 1})]
    json.dump({"model": model_path, "lp_count": len(rows), "rows": rows}, open(out, "w"), indent=1)
    for r in rows:
        print(f"{r['scenario']:26s} {r['status']:8s} max R1405 = {r['max_R1405']}")
        if "fluxes" in r:
            print("   ", {k: v for k, v in r["fluxes"].items() if abs(v) > 1e-9})
