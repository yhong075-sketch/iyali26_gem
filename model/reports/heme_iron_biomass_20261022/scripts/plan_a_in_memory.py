"""In memory: Yeast8 values (heme a [C_cy] 1e-6 + Fe [C_cy] 3.04e-5) added to biomass_C. No file writes."""
import sys, warnings
from cobra.flux_analysis import pfba
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
model_path, media, profile = sys.argv[1:4]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=model_path, media_path=media, strain_profile_path=profile)
m = ctx.model; m.solver = "glpk"
m.reactions.biomass_C.add_metabolites({m.metabolites.get_by_id("m1635[C_cy]"): -1e-6, m.metabolites.get_by_id("m948[C_cy]"): -3.04e-5})
s = pfba(m)
print("status", s.status, "growth %.6f" % s.fluxes["biomass_C"], " ".join(f"{k}={s.fluxes[k]:.3g}" for k in ["R1189", "R2193", "R863", "R302", "R384", "R383", "R1709"]))
print("iron use / cap = %.3f" % (-s.fluxes["R1189"] / 1.61e-5))
