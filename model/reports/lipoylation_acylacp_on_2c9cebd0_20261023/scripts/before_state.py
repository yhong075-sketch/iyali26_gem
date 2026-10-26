"""Baseline state of every lipoylation and acyl-ACP target in a model (read-only).

Usage (from the worktree root): python -I before_state.py MODEL.xml LIPOYLATION.json ACYL_ACP.json OUT.json
Records the SHA-256, each target species (formula, charge, name), each target reaction (stoichiometry,
bounds, GPR, mass/charge residual), the four lipoylation biomass_C coefficients, and whether the new
reactions, species and genes are present.
"""
import hashlib
import json
import sys
import warnings

from cobra.io import read_sbml_model

model_path, lipo_path, acyl_path, out = sys.argv[1:5]
lipo, acyl = json.load(open(lipo_path)), json.load(open(acyl_path))
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = read_sbml_model(model_path)


def residual(r):
    try:
        return {k: v for k, v in r.check_mass_balance().items() if abs(v) > 1e-9}
    except ValueError as exc:
        return f"unbalanceable: {exc}"


species = sorted({*lipo["formula_fills"], *lipo["species"], *acyl["formula_fills"], *acyl["species"]})
reactions = sorted({*lipo["edited_reactions"], *lipo["formula_fill_users"]["m179[C_mi]"],
                    *lipo["formula_fill_users"]["m586[C_mi]"], *acyl["block_reactions"]})
biomass = {m.id: c for m, c in model.reactions.get_by_id(lipo["biomass_reaction"]).metabolites.items()}
state = {
    "model": model_path, "sha256": hashlib.sha256(open(model_path, "rb").read()).hexdigest(),
    "species": {mid: ({"formula": m.formula or None, "charge": m.charge, "name": m.name,
                       "reactions": sorted(r.id for r in m.reactions)}
                      if (m := model.metabolites.get_by_id(mid)) else None) for mid in species},
    "reactions": {rid: {"stoichiometry": {m.id: c for m, c in sorted(r.metabolites.items(), key=lambda x: x[0].id)},
                        "bounds": list(r.bounds), "gpr": r.gene_reaction_rule, "residual": residual(r)}
                  for rid in reactions for r in [model.reactions.get_by_id(rid)]},
    "biomass_terms": {mid: biomass.get(mid, 0.0) if mid in model.metabolites else None for mid in lipo["biomass_terms"]},
    "biomass_C_metabolite_count": len(biomass),
    "new_present": {"reactions": sorted(r for r in lipo["new_reactions"] if r in model.reactions),
                    "metabolites": sorted(m for m in lipo["new_metabolites"] if m in model.metabolites),
                    "genes": sorted(g for g in lipo["new_genes"] if g in model.genes)},
    "context_genes_present": sorted(g for g in lipo["context_genes"] if g in model.genes),
}
with open(out, "w") as handle:
    json.dump(state, handle, indent=1, default=str)
    handle.write("\n")
print(state["sha256"], len(species), "species", len(reactions), "reactions", "new present:", state["new_present"])
