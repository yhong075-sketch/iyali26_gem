"""Stage 1 addendum: whole-network pFBA difference, closed vs R1725 open / R1726 open / both open (read-only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    stage1_pfba_diff.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json
Task A medium + PO1f, Gurobi 1e-9. Lists every reaction whose pFBA flux changes by more than 1e-9, and the
net stoichiometry of the changed fluxes (what the opened route substitutes for). pFBA optima need not be
unique; the difference is one representative, not proof of a unique route.
"""
import json
import sys
import warnings

from cobra.flux_analysis import pfba

from checked_solve import configure
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

model_path, medium, profile, out = sys.argv[1:5]


def fluxes(ids):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
    configure(ctx.model)
    for rid in ids:
        ctx.model.reactions.get_by_id(rid).bounds = (0.0, 1000.0)
    return ctx.model, pfba(ctx.model).fluxes


model, base = fluxes(())
result = {}
for label, ids in (("R1725_open", ("R1725",)), ("R1726_open", ("R1726",)), ("both_open", ("R1725", "R1726"))):
    _, new = fluxes(ids)
    delta = {rid: float(new[rid] - base[rid]) for rid in base.index if abs(new[rid] - base[rid]) > 1e-9}
    net = {}
    for rid, d in delta.items():
        for met, coef in model.reactions.get_by_id(rid).metabolites.items():
            net[met.id] = net.get(met.id, 0.0) + coef * d
    result[label] = {
        "changed": [{"reaction": rid, "name": model.reactions.get_by_id(rid).name[:50], "closed": float(base[rid]),
                     "opened": float(base[rid] + d), "delta": d,
                     "equation": model.reactions.get_by_id(rid).build_reaction_string(True)[:160]}
                    for rid, d in sorted(delta.items(), key=lambda kv: -abs(kv[1]))],
        "net_of_changes": {k: v for k, v in sorted(net.items(), key=lambda kv: -abs(kv[1])) if abs(v) > 1e-9}}
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
for label, row in result.items():
    print("##", label, len(row["changed"]), "changed")
    for c in row["changed"][:40]:
        print(f"  {c['reaction']:12s} {c['closed']:+.6g} -> {c['opened']:+.6g}  {c['name']}")
    print("  net:", {k: round(v, 7) for k, v in list(row["net_of_changes"].items())[:25]})
