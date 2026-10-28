"""Checks on the built files, baseline 9d655c2b vs the UM7 diffusion candidate (read-only; writes OUT.json).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts \
    /opt/miniconda3/bin/python built_checks.py BASELINE.xml CANDIDATE.xml MEDIUM.csv PROFILE.json OUT.json
Both models are read from short paths (cobra rejects str paths of 260+ characters on darwin). Gurobi 1e-9,
real HOME. Per model:
- WT on Task A + PO1f (checked_solve: reset, residual / bound checks) and pFBA flux of the watched reactions;
- energy-generating-cycle LPs: every exchange (0, 1000), every positive lower bound relaxed (xMAINTENANCE),
  max R_NTP1 and max biomass_C, solver status reported;
- FVA at 0 % and 99.9 % of the watched reactions on Task A, and at 0 % with every exchange open;
- whole-model blocked sets (FVA at 0 % of every reaction, |v| <= 1e-9), open and Task A;
- mass and charge balance of the two new reactions (candidate only).
"""
import json, math, sys, warnings
import cobra
from cobra.flux_analysis import flux_variability_analysis, pfba
from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

UM7 = ["R92", "R396", "R567", "R568", "R569", "R579", "R580", "R965", "R1000", "R1001", "R1017", "R1055",
       "R1082", "R1212", "R1214", "R1267", "R1270"]
WATCH = UM7 + ["R1912", "R576", "R578", "R419", "R1170", "R24", "R90"]
NEW = ["R_ADPRtn", "R_PPItn"]
TOL = 1e-9
base_path, cand_path, medium, profile, out = sys.argv[1:6]


def load(path, task_a):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if task_a:
            model = load_effective_simulation_context(model_path=path, media_path=medium,
                                                      strain_profile_path=profile).model
        else:
            model = cobra.io.read_sbml_model(path)
            for r in model.exchanges:
                r.bounds = (-1000, 1000)
    configure(model)
    return model


def fva(model, ids, frac):
    df = flux_variability_analysis(model, reaction_list=ids, fraction_of_optimum=frac, processes=1)
    return {r: [float(df.loc[r, "minimum"]), float(df.loc[r, "maximum"])] for r in ids}


def blocked(model):
    df = flux_variability_analysis(model, fraction_of_optimum=0.0, processes=1)
    return sorted(r for r in df.index if max(abs(df.loc[r, "minimum"]), abs(df.loc[r, "maximum"])) <= TOL)


def egc(model):
    out = {}
    with model:
        for r in model.exchanges:
            r.bounds = (0.0, 1000.0)
        relaxed = []
        for r in model.reactions:
            if r.lower_bound > 0:
                relaxed.append([r.id, r.lower_bound]); r.lower_bound = 0.0
        for objective in ("R_NTP1", "biomass_C"):
            model.objective = objective
            model.solver.problem.reset()
            v = model.slim_optimize(error_value=float("nan"))
            out[objective] = {"status": model.solver.status, "max_no_uptake": None if math.isnan(v) else v}
        out["relaxed"] = relaxed
    return out


result = {}
for label, path in (("baseline", base_path), ("candidate", cand_path)):
    m = load(path, True)
    ids = WATCH + [r for r in NEW if r in m.reactions]
    row = {"wt": solve(m, matrix(m), set())}
    m.solver.problem.reset()
    sol = pfba(m)
    row["pfba_watch"] = {r: float(sol.fluxes[r]) for r in ids}
    row["egc"] = egc(m)
    row["fva_taskA_0"] = fva(m, ids, 0.0)
    row["fva_taskA_999"] = fva(m, ids, 0.999)
    row["blocked_taskA"] = blocked(m)
    if label == "candidate":
        row["balance"] = {r: m.reactions.get_by_id(r).check_mass_balance() for r in NEW}
        row["new_reactions"] = {r: {"equation": m.reactions.get_by_id(r).build_reaction_string(),
                                    "bounds": list(m.reactions.get_by_id(r).bounds),
                                    "gpr": m.reactions.get_by_id(r).gene_reaction_rule} for r in NEW}
    o = load(path, False)
    o.solver.problem.reset()
    row["open_growth"] = o.slim_optimize(error_value=float("nan"))
    row["fva_open_0"] = fva(o, ids, 0.0)
    row["blocked_open"] = blocked(o)
    for cond in ("open", "taskA"):
        row[f"um7_blocked_{cond}"] = sorted(r for r in UM7 if r in row[f"blocked_{cond}"])
    result[label] = row
    print(label, "WT", row["wt"]["status"], repr(row["wt"]["growth"]), "passed", row["wt"]["passed"],
          "| EGC", {k: (v["status"], v["max_no_uptake"]) for k, v in row["egc"].items() if k != "relaxed"},
          "| blocked open", len(row["blocked_open"]), "taskA", len(row["blocked_taskA"]),
          "| UM7 blocked open", len(row["um7_blocked_open"]), "taskA", len(row["um7_blocked_taskA"]), flush=True)
b, c = result["baseline"], result["candidate"]
for cond in ("open", "taskA"):
    result[f"{cond}_freed"] = sorted(set(b[f"blocked_{cond}"]) - set(c[f"blocked_{cond}"]))
    result[f"{cond}_newly_blocked"] = sorted(set(c[f"blocked_{cond}"]) - set(b[f"blocked_{cond}"]) - set(NEW))
result["wt_relative_change"] = (c["wt"]["growth"] - b["wt"]["growth"]) / b["wt"]["growth"]
print({k: v for k, v in result.items() if k not in ("baseline", "candidate")})
print("balance", c["balance"])
json.dump(result, open(out, "w"), indent=1)
