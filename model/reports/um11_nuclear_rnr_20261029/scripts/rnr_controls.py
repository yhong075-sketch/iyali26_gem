"""Targeted controls for the RNR genes on Task A + PO1f (read-only): KO, R702-R705 closed, and combined.

Usage (from platform/): PYTHONPATH=... python rnr_controls.py MODEL.xml MEDIUM.csv PROFILE.json VARIANT
"""
import sys, warnings
from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from um11_apply import apply

model_path, medium, profile, variant = sys.argv[1:5]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    m = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile).model
apply(m, variant)
configure(m)
S = matrix(m)
NTP_RED = ["R702", "R703", "R704", "R705"]
RNR = ["YALI1D04261g", "YALI1F39330g", "YALI1B05072g"]
wt = solve(m, S, set())["growth"]
cases = [("WT", [], [])] + [(f"KO {g}", [g], []) for g in RNR] + [("KO all three RNR genes", RNR, [])]
cases += [("R702-R705 closed", [], NTP_RED)] + [(f"KO {g} + R702-R705 closed", [g], NTP_RED) for g in RNR]
cases += [("KO all three + R702-R705 closed", RNR, NTP_RED)]
for name, genes, closed in cases:
    with m:
        before = {r.id: r.bounds for r in m.reactions}
        for g in genes:
            m.genes.get_by_id(g).knock_out()
        for rid in closed:
            m.reactions.get_by_id(rid).bounds = (0, 0)
        knocked = {rid for rid, b in before.items() if m.reactions.get_by_id(rid).bounds == (0, 0) and b != (0, 0)}
        k = solve(m, S, knocked)
    print(f"{variant:12s} {name:40s} status={k['status']} growth={k['growth']} ratio={None if k['growth'] is None else k['growth']/wt:.6g} passed={k['passed']}")
