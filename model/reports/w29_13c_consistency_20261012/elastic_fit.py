"""Task C elastic fit of the model to 22 W29 glucose 13C-MFA flux ranges (diagnostic, in memory).

Model: Task B candidate + holistic overlay; glucose fixed at 2.414, biomass_C fixed at f x 0.2663
(f = 1, 0.95, 0.90). Each range gets slacks s+ / s- (per 100 glucose); the objective minimises
sum((s+ + s-) / CI width). Then FVA on the optimal face (objective <= optimum + 1e-6) gives the range
of every slack and every mapped flux, so judgments do not rest on one optimum.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/w29_13c_consistency_20261012/elastic_fit.py
Writes mfa_mapping.csv, elastic_fit.json and elastic_fit_summary.csv here. Never overwrites.
"""
import csv
import json
import subprocess
import time
from datetime import datetime, timezone

from common import HERE, MU, REPO, UPTAKE, Counter_, apply_overlay, check_pins, growth, load, mfa_rows, r825_reversible

# (key, MFA row id(s), name, {reaction: coefficient}); fixed in TASK.md before any run.
MAPPING = [
    ("PGI", "R4 net", "glucose-6-phosphate isomerase", {"R326": 1}),
    ("PFK", "R5 net", "phosphofructokinase net of FBPase", {"R636": 1, "R312": -1}),
    ("GAPDH", "R8 net", "glyceraldehyde-3-P dehydrogenase", {"R344": 1}),
    ("PYK", "R10 net", "pyruvate kinase", {"R694": 1}),
    ("PEPCK", "R32", "PEP carboxykinase", {"R634": 1}),
    ("G6PDH", "R11", "glucose-6-P dehydrogenase", {"R325": 1}),
    ("PDH", "R19", "pyruvate dehydrogenase", {"R693": 1}),
    ("CS", "R20", "citrate synthase (mito)", {"R236": 1}),
    ("IDH_total", "R21 net", "total mitochondrial IDH (= aconitase net in the MFA network)", {"R487": 1, "R1385": 1}),
    ("OGDH", "R25", "2-oxoglutarate dehydrogenase", {"R614": 1}),
    ("SDH", "R26 net", "succinate dehydrogenase", {"R740": 1}),
    ("MDH_mito", "R28 net", "malate dehydrogenase (mito)", {"R533": 1}),
    ("ICL", "R43.2", "isocitrate lyase", {"R490": 1}),
    ("MS", "R45", "malate synthase (cyt + perox)", {"R536": 1, "R537": 1}),
    ("MDH_cyt_perox", "R46", "malate dehydrogenase (cyt + perox)", {"R534": 1, "R535": 1}),
    ("ACS", "R37 net", "acetyl-CoA synthetase (cyt + perox)", {"R90": 1, "R2268": 1}),
    ("pyruvate_import", "R18", "pyruvate cyt -> mito (net)", {"R1311": 1, "R842": -1}),
    ("succinate_import", "R44", "succinate cyt -> mito (net)", {"R928": 1, "R929": 1}),
    ("malate_import", "R47", "malate cyt -> mito (net, all carriers)", {"R807": 1, "R902": 1, "R1386": 1, "R831": -1}),
    ("citrate_export", "cit transp", "citrate mito -> cyt (net, all carriers)", {"R817": 1, "R831": -1, "R832": -1, "R833": -1}),
    ("GCV", "R51 net", "glycine cleavage, MFA synthesis sign", {"R357": -1, "R360": -1, "R362": -1}),
    ("PC_minus_ME", "R31 - R29 net", "pyruvate carboxylase minus mito malic enzyme", {"R690": 1, "R538": -1}),
]
INFO = ("R1894", "R825", "R1275", "R1205", "R193", "R191", "R83", "R412", "R1172", "R410", "R945", "R1974",
        "R89", "R1708", "R570", "R2062")
SETUPS = (("overlay", False, True), ("overlay_no_malate", False, False),
          ("R825_rev", True, True), ("R825_rev_no_malate", True, False))
FRACTIONS = (1.0, 0.95, 0.90)
RECALLED_V1 = {"overlay": 6.94, "overlay_no_malate": 1.14, "R825_rev": 2.72, "R825_rev_no_malate": 0.68}


def ranges():
    rows = mfa_rows()
    out = {}
    for key, row_id, _, _ in MAPPING:
        if row_id == "R31 - R29 net":
            lb = float(rows["R31"]["glucose_LB"]) - float(rows["R29 net"]["glucose_UB"])
            ub = float(rows["R31"]["glucose_UB"]) - float(rows["R29 net"]["glucose_LB"])
            flux = float(rows["R31"]["glucose_flux"]) - float(rows["R29 net"]["glucose_flux"])
        else:
            r = rows[row_id]
            lb, ub, flux = float(r["glucose_LB"]), float(r["glucose_UB"]), float(r["glucose_flux"])
        assert ub > lb, key
        out[key] = (lb, ub, flux)
    return out


def write_mapping(rng):
    with (HERE / "mfa_mapping.csv").open("x", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["key", "mfa_row", "name", "model_expression", "mfa_flux", "mfa_lb", "mfa_ub", "evidence_tier"])
        for key, row_id, name, expr in MAPPING:
            lb, ub, flux = rng[key]
            w.writerow([key, row_id, name, " ".join(f"{c:+d}*{r}" for r, c in expr.items()), round(flux, 4),
                        round(lb, 4), round(ub, 4), "inferred (mapping is this task's judgment)"])


def fit(model, rng, keys, mu):
    """Add slacks; minimise; FVA on the optimal face. Caller holds the model context."""
    per100 = 100.0 / UPTAKE
    prob = model.problem
    model.reactions.biomass_C.bounds = (mu, mu)
    expr, terms, extra = {}, [], []
    for key, _, _, coeffs in MAPPING:
        expr[key] = per100 * sum(c * model.reactions.get_by_id(r).flux_expression for r, c in coeffs.items())
    for key in keys:
        lb, ub, _ = rng[key]
        sp, sm = prob.Variable(f"sp_{key}", lb=0), prob.Variable(f"sm_{key}", lb=0)
        extra += [sp, sm, prob.Constraint(expr[key] - sp, ub=ub, name=f"up_{key}"),
                  prob.Constraint(expr[key] + sm, lb=lb, name=f"lo_{key}")]
        terms.append((sp + sm) / (ub - lb))
    model.add_cons_vars(extra)
    total = sum(terms)
    model.objective = prob.Objective(total, direction="min")
    Counter_.lps += 1
    model.slim_optimize()
    status = model.solver.status
    if status != "optimal":
        return {"status": status}
    opt = model.solver.objective.value
    point = {k: float(v.primal) for k, v in model.solver.variables.items() if k.startswith(("sp_", "sm_"))}
    point_expr = {key: per100 * sum(c * model.reactions.get_by_id(r).flux for r, c in coeffs.items())
                  for key, _, _, coeffs in MAPPING}
    model.add_cons_vars([prob.Constraint(total, ub=opt + 1e-6, name="face")])
    face = {}

    def span(e):
        out = []
        for d in ("min", "max"):
            model.objective = prob.Objective(e, direction=d)
            Counter_.lps += 1
            model.slim_optimize()
            out.append(model.solver.objective.value if model.solver.status == "optimal" else None)
        return out
    for key in keys:
        face[f"sp_{key}"] = span(model.solver.variables[f"sp_{key}"])
        face[f"sm_{key}"] = span(model.solver.variables[f"sm_{key}"])
    for key in expr:
        face[f"flux_{key}"] = span(expr[key])
    for rid in INFO:
        face[f"info_{rid}"] = span(per100 * model.reactions.get_by_id(rid).flux_expression)
    return {"status": status, "total_violation": opt, "slack_point": point, "flux_point": point_expr, "face": face}


def main():
    inputs = check_pins()
    rng = ranges()
    write_mapping(rng)
    started = time.monotonic()
    record = {"started_utc": datetime.now(timezone.utc).isoformat(),
              "command": "python -B model/reports/w29_13c_consistency_20261012/elastic_fit.py",
              "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
              "inputs": inputs, "ranges_per_100_glucose": rng, "mapping": MAPPING, "runs": []}
    for medium in ("V2", "V1"):
        sim, model = load(medium)
        record.setdefault("simulation_context", {})[medium] = sim.provenance()
        for name, rev, with_malate in SETUPS:
            keys = [k for k, *_ in MAPPING if with_malate or k != "malate_import"]
            with model:
                apply_overlay(model)
                if rev:
                    r825_reversible(model)
                wt, wt_status = growth(model)
                for f in FRACTIONS:
                    with model:
                        out = fit(model, rng, keys, f * MU)
                    run = {"medium": medium, "setup": name, "fraction_of_measured_growth": f, "wt_max": wt,
                           "wt_status": wt_status, "ranges_used": keys, **out}
                    record["runs"].append(run)
                    print(medium, name, f, wt, out.get("status"), out.get("total_violation"), flush=True)
    checks = []
    for name, expected in RECALLED_V1.items():
        run = next(r for r in record["runs"] if r["medium"] == "V1" and r["setup"] == name
                   and r["fraction_of_measured_growth"] == 1.0)
        got = run.get("total_violation")
        rel = None if got is None else got / expected - 1
        checks.append({"setup": name, "recalled": expected, "observed": got, "relative_difference": rel,
                       "within_2pct": rel is not None and abs(rel) <= 0.02})
    record["v1_reproduction_checks"] = checks
    record["v1_stop_fired"] = not all(c["within_2pct"] for c in checks)
    record.update(lps=Counter_.lps, elapsed_seconds=time.monotonic() - started)
    (HERE / "elastic_fit.json").open("x").write(json.dumps(record, indent=1, default=str) + "\n")
    with (HERE / "elastic_fit_summary.csv").open("x", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["medium", "setup", "fraction", "total_violation", "range", "mfa_lb", "mfa_ub", "flux_point",
                    "flux_min_on_face", "flux_max_on_face", "above_ub_point", "above_ub_face_min", "above_ub_face_max",
                    "below_lb_point", "below_lb_face_min", "below_lb_face_max"])
        for run in record["runs"]:
            if run.get("status") != "optimal":
                w.writerow([run["medium"], run["setup"], run["fraction_of_measured_growth"], run.get("status")])
                continue
            for key, *_ in MAPPING:
                lb, ub, _ = rng[key]
                fl = run["face"][f"flux_{key}"]
                row = [run["medium"], run["setup"], run["fraction_of_measured_growth"], round(run["total_violation"], 4),
                       key, round(lb, 3), round(ub, 3), round(run["flux_point"][key], 3), round(fl[0], 3), round(fl[1], 3)]
                if key in run["ranges_used"]:
                    for side in ("sp", "sm"):
                        face = run["face"][f"{side}_{key}"]
                        row += [round(run["slack_point"][f"{side}_{key}"], 3), round(face[0], 3), round(face[1], 3)]
                else:
                    row += [None] * 6
                w.writerow(row)
    print("V1 checks:", checks, "| stop fired:", record["v1_stop_fired"], "| LPs", Counter_.lps,
          "| elapsed", round(record["elapsed_seconds"]), "s")


if __name__ == "__main__":
    main()
