"""Do the Yeast-GEM vacuolar H+ diffusion and H+-coupled amino-acid transports form no-uptake loops? (read-only)"""
import sys, warnings
from pathlib import Path
import cobra
from cobra.flux_analysis import flux_variability_analysis
warnings.simplefilter("ignore")
m = cobra.io.read_sbml_model(Path(sys.argv[1]))
m.solver = "gurobi"
ids = ["r_1086", "r_1831", "r_2099", "r_1187", "r_1188", "r_1189", "r_1195", "r_1197", "r_1198", "r_1202", "r_1203",
       "r_1204", "r_1208", "r_1209", "r_1210", "r_1220", "r_1221", "r_1222", "r_1252", "r_1261", "r_1263"]
for r in m.reactions:
    if r.boundary:
        r.bounds = (0, max(r.upper_bound, 0))
    elif r.lower_bound > 0:
        r.lower_bound = 0
for r in ids:
    x = m.reactions.get_by_id(r)
    print(r, x.bounds, x.build_reaction_string(use_metabolite_names=True)[:110])
for label, closed in (("as published", []), ("r_1831 closed", ["r_1831"])):
    with m:
        for c in closed:
            m.reactions.get_by_id(c).bounds = (0, 0)
        m.objective = m.reactions.r_1086  # dummy objective; FVA at 0 %
        df = flux_variability_analysis(m, reaction_list=ids, fraction_of_optimum=0.0, processes=1)
        nz = {r: [round(df.loc[r, 'minimum'], 3), round(df.loc[r, 'maximum'], 3)] for r in ids
              if max(abs(df.loc[r, 'minimum']), abs(df.loc[r, 'maximum'])) > 1e-9}
        print(label, "no-uptake nonzero:", nz)
