"""pFBA, serial FVA (0 % and 99.9 % of optimum) and single-reaction KO for a list of reactions.
Usage: python flux_evidence.py --model M.xml --medium sd_leu.csv [--profile po1f.json] --ids ids.csv --out prefix
ids.csv needs an 'id' column. Writes <prefix>_flux_evidence.csv and <prefix>_pfba.csv."""
import os, sys, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_model, read_condition, setup, growth

ap = argparse.ArgumentParser()
ap.add_argument("--model", required=True); ap.add_argument("--medium", required=True); ap.add_argument("--profile")
ap.add_argument("--ids", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
import pandas as pd
from cobra.flux_analysis import pfba, flux_variability_analysis
M = load_model(a.model); medium, closed = read_condition(a.medium, a.profile); setup(M, medium, closed)
mu = growth(M); print("WT growth", mu, flush=True)
F = pfba(M).fluxes; F.to_csv(a.out + "_pfba.csv")
ids = [i for i in pd.read_csv(a.ids).id if i in M.reactions]
f0 = flux_variability_analysis(M, ids, fraction_of_optimum=0.0, processes=1)
f1 = flux_variability_analysis(M, ids, fraction_of_optimum=0.999, processes=1)
out = pd.DataFrame({"id": ids})
out["pfba_flux"] = out.id.map(F)
out["fva_min_0pct_opt"] = out.id.map(f0["minimum"]); out["fva_max_0pct_opt"] = out.id.map(f0["maximum"])
out["fva_min_99.9pct_opt"] = out.id.map(f1["minimum"]); out["fva_max_99.9pct_opt"] = out.id.map(f1["maximum"])
ko = []
for i in ids:
    with M:
        M.reactions.get_by_id(i).bounds = (0, 0); ko.append(growth(M))
out["ko_growth"] = ko; out["ko_wt_ratio"] = out.ko_growth / mu
tol = 1e-6
def status(r):
    if abs(r["fva_min_0pct_opt"]) < tol and abs(r["fva_max_0pct_opt"]) < tol: return "blocked (cannot carry flux)"
    if r["fva_min_99.9pct_opt"] > tol or r["fva_max_99.9pct_opt"] < -tol: return "required at WT optimum"
    if abs(r["pfba_flux"]) > tol: return "active in WT pFBA, not required"
    return "can carry flux; 0 in WT pFBA"
out["flux_status"] = out.apply(status, axis=1)
out.to_csv(a.out + "_flux_evidence.csv", index=False)
print(out.flux_status.value_counts().to_dict())
