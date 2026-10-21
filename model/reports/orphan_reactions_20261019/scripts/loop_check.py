"""Internal-cycle test for reactions that reached +-1000 under SD-Leu + PO1f.
A reaction that can carry flux with every boundary reaction closed (bounds 0,0), forced
fluxes (lb > 0, e.g. NGAM) relaxed to include 0, and the biomass not required lies on a stoichiometric internal cycle (net-zero loop). LP only (GLPK)."""
import sys
sys.path.insert(0, sys.argv[1])
from common import load_model, read_condition, setup, growth
from cobra.flux_analysis import flux_variability_analysis
M = load_model(sys.argv[2]); medium, closed = read_condition(sys.argv[3], sys.argv[4]); setup(M, medium, closed)
print("WT", growth(M))
for r in ["R1964", "R2041", "R2042", "R1901", "R1308", "R321", "R218", "R1043", "R52", "R1210", "R2197", "R530", "R1727", "R1337", "R751"]:
    x = M.reactions.get_by_id(r); print(r, x.bounds, x.gene_reaction_rule or "-", "|", x.build_reaction_string(use_metabolite_names=True))
ids = ["R1964", "R2041", "R2042"]
with M:
    for b in M.boundary:
        b.bounds = (0, 0)
    forced = [r.id for r in M.reactions if r.lower_bound > 0 or r.upper_bound < 0]
    print("forced-flux reactions relaxed to include 0:", forced)
    for rid in forced:
        r = M.reactions.get_by_id(rid); r.bounds = (min(r.lower_bound, 0), max(r.upper_bound, 0))
    f = flux_variability_analysis(M, ids, fraction_of_optimum=0.0, processes=1)
    print("all boundaries closed:\n", f)
    # which reactions form the cycle: maximise R2042 and list nonzero fluxes
    M.objective = "R2042"; s = M.optimize()
    nz = s.fluxes[s.fluxes.abs() > 1e-6]
    for rid, v in nz.items():
        print(f"  {rid} {v:.1f} | {M.reactions.get_by_id(rid).gene_reaction_rule or '-'} | {M.reactions.get_by_id(rid).build_reaction_string(use_metabolite_names=True)[:160]}")
