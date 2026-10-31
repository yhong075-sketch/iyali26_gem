"""Build the trehalose biomass before/after Escher map and its scenario data (no LP: reads flux.json).

Usage: python build.py FLUX.json WORKDIR
Writes WORKDIR/maps/r171.json and _deco.json (file names kept from the R171 figure so render.py's job layout is
unchanged), WORKDIR/data/<scenario>.json and WORKDIR/index.json. Sections: the trehalose branch with the biomass_C
trehalose drain drawn as its own arrow; the G6P and UDP-glucose reactions that carry flux in any scenario. Adapted from
the aconitase figure's build.py; sections, labels, scenarios and title changed.
"""
import json
import os
import sys

from maplib import build_map

flux_path, work = sys.argv[1:]
F = json.load(open(flux_path))
RX, MT, FL = F["reactions"], F["metabolites"], F["fluxes"]
DRAIN = F["mapped"]["drain"]
KEYS = ("wt_baseline", "wt_candidate", "ko_baseline", "ko_candidate")
active = [rid for rid in F["mapped"]["supply"] if any(abs(FL[k][rid]) > 1e-9 for k in KEYS)]
idle = [rid for rid in F["mapped"]["supply"] if rid not in active]
SECTIONS = [
    ("Trehalose branch: R143 (TPS1 and TPS2 and TPS3), R767, trehalase R142, R1349 (closed), export R1013, "
     "biomass_C trehalose term (removed)",
     F["mapped"]["branch"] + [DRAIN]),
    ("G6P and UDP-glucose reactions carrying flux in any scenario (hexokinase, isomerase, phosphoglucomutase, "
     "UDP-glucose pyrophosphorylase, 1,3-beta-glucan synthase)", active),
]
mapped = [rid for _, ids in SECTIONS for rid in ids]


def clean(v):
    return 0.0 if abs(v) < 1e-9 else v


rxns = {}
for rid in mapped:
    if rid == DRAIN:
        coef = F["baseline_trehalose_coefficient"]
        rxns[rid] = {"name": "biomass_C trehalose term", "lb": 0.0, "ub": 1000.0, "gpr": "", "gpr_baseline": "",
                     "genes": [], "subsystem": "", "equation": f"{-coef} trehalose -> biomass (before); absent after",
                     "mets": {"m294[C_cy]": -1.0, "biomass_trehalose[C_cy]": 1.0},
                     "sublabel": "present before, removed after"}
        continue
    r = dict(RX[rid])
    r["sublabel"] = ""
    rxns[rid] = r
MT = {**MT, "biomass_trehalose[C_cy]": {"name": "biomass (trehalose share)", "comp": "C_cy", "degree": 1}}
mets = {mid: dict(name=MT[mid]["name"], comp=MT[mid]["comp"], degree=MT[mid]["degree"], gap=False, cls="")
        for r in rxns.values() for mid in r["mets"]}
met_data = {mid: 0 for mid in mets}
legend = [
    "Reaction colour, width and label = pFBA flux at the optimum (Task A medium + PO1f; Gurobi 1e-9). The only model "
    "difference is the removed biomass_C trehalose term (0.00686 mmol/gDW) and its notes.",
    "Knockout: YALI1E17616g (TPS1) closes R143. Layout is automatic. Not drawn (zero flux in all four scenarios): "
    + ", ".join(idle) + ".",
]
title = (f"Trehalose in biomass_C: before {F['baseline_sha256'][:12]}… vs after {F['candidate_sha256'][:12]}…")
mj, dj, _, W, H, placed = build_map("iyali26-trehalose-biomass", title, rxns, mets, SECTIONS, maxw=3200, legend=legend,
                                    dx=560)
scenarios = {"pfba_wt_before": {rid: clean(FL["wt_baseline"][rid]) for rid in mapped},
             "pfba_wt_after": {rid: clean(FL["wt_candidate"][rid]) for rid in mapped},
             "pfba_ko_before": {rid: clean(FL["ko_baseline"][rid]) for rid in mapped},
             "pfba_ko_after": {rid: clean(FL["ko_candidate"][rid]) for rid in mapped}}
tips = {rid: (f"GPR: {rxns[rid]['gpr'] or 'none'} · pFBA WT before {scenarios['pfba_wt_before'][rid]:.3g}, after "
              f"{scenarios['pfba_wt_after'][rid]:.3g} · TPS1 KO before {scenarios['pfba_ko_before'][rid]:.3g}, after "
              f"{scenarios['pfba_ko_after'][rid]:.3g}") for rid in mapped}
os.makedirs(f"{work}/maps", exist_ok=True)
os.makedirs(f"{work}/data", exist_ok=True)
json.dump(mj, open(f"{work}/maps/r171.json", "w"))
json.dump(dj, open(f"{work}/maps/r171_deco.json", "w"))
for name, data in scenarios.items():
    json.dump(data, open(f"{work}/data/{name}.json", "w"))
counts = {
    "mapped_reactions": len(mapped),
    "not_drawn_zero_flux": idle,
    "growth": {k: v["checked"]["growth"] for k, v in F["scenarios"].items()},
    "branch_flux": {k: {rid: clean(FL[k][rid]) for rid in ("R143", "R767", DRAIN)} for k in KEYS},
}
json.dump({"title": title, "W": W, "H": H, "met_data": met_data, "tips": tips,
           "sections": [(lab, len(ids)) for lab, ids in SECTIONS], "counts": counts},
          open(f"{work}/index.json", "w"), indent=1)
print(f"map {W}x{H}, {len(placed)} placed nodes;", json.dumps(counts))
