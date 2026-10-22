"""Reproduce Gap 5 in memory: WT, heme b + heme a in biomass_C, and two controls. No file writes."""
import sys, warnings, json
from cobra.flux_analysis import pfba
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

model_path, media, profile, solver = sys.argv[1:5]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=model_path, media_path=media, strain_profile_path=profile)
m = ctx.model
m.solver = solver
print("model", ctx.canonical_model_sha256[:8], "medium", ctx.medium_sha256[:8], "profile", ctx.strain_profile_sha256[:8], "solver", solver)
WATCH = ["R1189", "R2193", "R863", "R302", "R384", "R383", "R674", "R727", "R1070"]
lp = 0

def run(label, heme_b=0.0, heme_a=0.0, heme_a_site="m609[C_mi]", iron_cap=None):
    global lp
    with m:
        b = m.reactions.biomass_C
        if heme_b:
            b.add_metabolites({m.metabolites.get_by_id("m525[C_mi]"): -heme_b})
        if heme_a:
            b.add_metabolites({m.metabolites.get_by_id(heme_a_site): -heme_a})
        if iron_cap is not None:
            m.reactions.R1189.lower_bound = -iron_cap
        s = pfba(m); lp += 2
        out = {k: float(s.fluxes[k]) for k in WATCH}
        print(f"{label:55s} status={s.status} growth={s.fluxes['biomass_C']:.6f} " +
              " ".join(f"{k}={v:.3g}" for k, v in out.items()))

run("WT (as supplied)")
run("heme b + heme a 5e-5 each (figure test)", 5e-5, 5e-5)
run("same, iron cap removed (R1189 >= -1000)", 5e-5, 5e-5, iron_cap=1000)
run("heme a 1e-6 only (legacy R1710 coefficient), mito", 0, 1e-6)
print("LPs", lp)
