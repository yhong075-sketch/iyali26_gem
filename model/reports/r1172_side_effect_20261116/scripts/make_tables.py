"""Small report tables from the step outputs (no LP). Writes ../tables/*.csv."""
import csv
import json
from pathlib import Path

from common import ART, NINE

T = Path(__file__).resolve().parents[1] / "tables"
O = ART / "out"
s1 = {s: json.load(open(O / f"step1_{s}.json")) for s in ("gurobi", "glpk")}
s2 = {s: json.load(open(O / f"step2_{s}.json")) for s in ("gurobi", "glpk")}
s5 = json.load(open(O / "step5_options.json"))
s5b = json.load(open(O / "step5b_routes.json"))
s6 = json.load(open(O / "step6_figure_flux.json"))

with open(T / "step1_reproduce.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["condition", "variant", "reaction", "gurobi_min", "gurobi_max", "gurobi_status", "glpk_min", "glpk_max",
                "glpk_status", "blocked_gurobi", "blocked_glpk"])
    for key, run in s1["gurobi"]["runs"].items():
        c, v = key.split("/")
        for r, x in run["fva"].items():
            y = s1["glpk"]["runs"][key]["fva"][r]
            w.writerow([c, v, r, x["min"], x["max"], f'{x["min_status"]}/{x["max_status"]}', y["min"], y["max"],
                        f'{y["min_status"]}/{y["max_status"]}', x["blocked"], y["blocked"]])

with open(T / "moiety_import_bounds.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["condition", "variant", "quantity", "gurobi_value", "gurobi_status", "glpk_value", "glpk_status"])
    for key, b in s2["gurobi"]["B"].items():
        c, v = key.split("/")
        for q, x in b.items():
            y = s2["glpk"]["B"][key][q]
            w.writerow([c, v, q, x["value"], x["status"], y["value"], y["status"]])

with open(T / "routes_moiety_channels.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["source", "condition_or_option", "variant", "target", "fixed_at", "status",
                "moiety_into_matrix_by_channel", "R1172", "R834", "R2025", "R2113", "R2077", "R2118", "R2112"])
    for key, routes in s2["gurobi"]["C"].items():
        c, v = key.split("/")
        for t, r in routes.items():
            nz = r.get("nonzero", {})
            w.writerow(["step2", c, v, t, r.get("fixed_at"), r.get("status"),
                        json.dumps(r.get("moiety_into_matrix_by_channel", {})),
                        *[nz.get(k, 0.0) for k in ("R1172", "R834", "R2025", "R2113", "R2077", "R2118", "R2112")]])
    for v in ("O1", "O1_R2075closed", "O3", "O4"):
        for t, r in s5b[v]["routes"].items():
            nz = r.get("nonzero", {})
            w.writerow(["step5b", "taskA", v, t, r.get("fixed_at"), r.get("status", r.get("note")),
                        json.dumps(r.get("moiety_channels", {})),
                        *[nz.get(k, 0.0) for k in ("R1172", "R834", "R2025", "R2113", "R2077", "R2118", "R2112")]])

with open(T / "options_summary.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["option", "wt_growth", "wt_status", "max_R25", "max_R95", "max_R552", "max_R2074", "max_R2075",
                "max_R2076", "maxabs_R2077", "maxabs_R2112", "max_R834", "energy_max_R_NTP1", "energy_max_biomass",
                "energy_status", "ko_class_changes_15pct_targeted", "n_targeted_ko"])
    V = dict(s5["variants"])
    V["O1_R2075closed"] = {"wt": s5b["O1_R2075closed"]["wt"]["value"], "wt_status": s5b["O1_R2075closed"]["wt"]["status"],
                           "fva": s5b["O1_R2075closed"]["fva"], "energy": s5b["O1_R2075closed"]["energy"], "ko": {}}
    base = s5["variants"]["B"]["ko"]
    for v in ("B", "O1", "O1_R2075closed", "O3", "O4"):
        x = V[v]
        fv = x["fva"]
        mx = lambda r: max(abs(fv[r]["min"]), abs(fv[r]["max"]))
        changes = [g for g, k in x["ko"].items() if k["class_15"] != base[g]["class_15"]] if x["ko"] else "not run"
        e = x["energy"]
        w.writerow([v, x["wt"], x["wt_status"], *[mx(r) for r in ("R25", "R95", "R552", "R2074", "R2075", "R2076",
                                                                   "R2077", "R2112", "R834")],
                    e["R_NTP1"]["value"], e["biomass_C"]["value"], f'{e["R_NTP1"]["status"]}/{e["biomass_C"]["status"]}',
                    json.dumps(changes), len(x["ko"]) if x["ko"] else 0])

with open(T / "targeted_knockouts.csv", "w", newline="") as f:
    w = csv.writer(f)
    vs = ["B", "O1", "O1b", "O3", "O4"]
    w.writerow(["gene"] + [f"{v}_ko_over_wt" for v in vs] + [f"{v}_class_15" for v in vs])
    for g in s5["genes"]:
        row = [g]
        for v in vs:
            k = s5["variants"][v]["ko"].get(g)
            row.append("" if k is None else k["ko_over_wt"])
        for v in vs:
            k = s5["variants"][v]["ko"].get(g)
            row.append("" if k is None else k["class_15"])
        w.writerow(row)
print("tables written to", T)
print(open(T / "options_summary.csv").read())
