"""Build the aconitase before/after Escher map and its scenario data (no LP: reads flux.json).

Usage: python build.py FLUX.json WORKDIR
Writes WORKDIR/maps/r171.json and _deco.json (file names kept from the R171 figure so render.py is unchanged),
WORKDIR/data/<scenario>.json and WORKDIR/index.json. Sections: the aconitase steps; the other reactions of citrate,
cis-aconitate and isocitrate; the reactions whose knockout flux changes most. Copied from the R673 figure's build.py;
sections, labels and title changed.
"""
import json
import os
import sys

from maplib import build_map

flux_path, work = sys.argv[1:]
F = json.load(open(flux_path))
RX, MT, FL = F["reactions"], F["metabolites"], F["fluxes"]
fixed, extra = F["mapped"]["fixed"], F["mapped"]["largest_ko_change"]
SECTIONS = [
    ("Aconitase steps: R238, R216 (mitochondrion) and R1556 (cytosol) rules changed; R239 (cytosol) already YALI1D11984g alone",
     ["R238", "R216", "R239", "R1556", "R735", "R736"]),
    ("The other reactions of citrate, cis-aconitate and isocitrate (cytosol and mitochondrion): synthase, dehydrogenases, lyase, transports",
     ["R236", "R1894", "R487", "R1385", "R488", "R490", "R831", "R832", "R833", "R817", "R1044", "R1046", "R1047",
      "R1875"]),
    ("The 12 other reactions whose knockout flux changes most between the models (growth collapse)", extra),
]
mapped = [rid for _, ids in SECTIONS for rid in ids]
assert sorted(mapped) == sorted(fixed + extra), "every mapped reaction must be drawn once"


def clean(v):
    return 0.0 if abs(v) < 1e-9 else v


rxns = {}
for rid in mapped:
    r = dict(RX[rid])
    r["sublabel"] = "RULE CHANGED" if rid in ("R238", "R216", "R1556") else ""
    rxns[rid] = r
mets = {mid: dict(name=MT[mid]["name"], comp=MT[mid]["comp"], degree=MT[mid]["degree"], gap=False, cls="")
        for r in rxns.values() for mid in r["mets"]}
met_data = {mid: 0 for mid in mets}
legend = [
    "Reaction colour, width and label = pFBA flux at the optimum (Task A medium + PO1f; Gurobi 1e-9). Only the rules "
    "of R238, R216 and R1556 differ between the two models; WT growth is identical.",
    "Before: YALI1D11984g or YALI1E17935g; after: YALI1D11984g (YALI1E17935g unlinked). Layout is automatic.",
]
title = (f"Aconitase rules: before {F['baseline_sha256'][:12]}… vs after {F['candidate_sha256'][:12]}…, "
         "YALI1D11984g knockout")
mj, dj, _, W, H, placed = build_map("iyali26-aconitase-aco2-split", title, rxns, mets, SECTIONS, maxw=4200, legend=legend,
                                    dx=560)
scenarios = {"pfba_wt": {rid: clean(FL["wt_candidate"][rid]) for rid in mapped},
             "pfba_ko_before": {rid: clean(FL["ko_baseline"][rid]) for rid in mapped},
             "pfba_ko_after": {rid: clean(FL["ko_candidate"][rid]) for rid in mapped}}
tips = {rid: (f"GPR after: {RX[rid]['gpr'] or 'none'} · GPR before: {RX[rid]['gpr_baseline'] or 'none'} · pFBA WT "
              f"{scenarios['pfba_wt'][rid]:.3g}, KO before {scenarios['pfba_ko_before'][rid]:.3g}, KO after "
              f"{scenarios['pfba_ko_after'][rid]:.3g}") for rid in mapped}
os.makedirs(f"{work}/maps", exist_ok=True)
os.makedirs(f"{work}/data", exist_ok=True)
json.dump(mj, open(f"{work}/maps/r171.json", "w"))
json.dump(dj, open(f"{work}/maps/r171_deco.json", "w"))
for name, data in scenarios.items():
    json.dump(data, open(f"{work}/data/{name}.json", "w"))
counts = {
    "mapped_reactions": len(mapped),
    "wt_equal_in_both_models": all(abs(FL["wt_baseline"][r] - FL["wt_candidate"][r]) <= 1e-9 for r in mapped),
    "ko_before_equal_to_wt": all(abs(FL["ko_baseline"][r] - FL["wt_baseline"][r]) <= 1e-9 for r in mapped),
    "growth": {k: v["checked"]["growth"] for k, v in F["scenarios"].items()},
    "changed_ko_before_vs_after": sorted(r for r in mapped if abs(FL["ko_baseline"][r] - FL["ko_candidate"][r]) > 1e-6),
}
json.dump({"title": title, "W": W, "H": H, "met_data": met_data, "tips": tips,
           "sections": [(lab, len(ids)) for lab, ids in SECTIONS], "counts": counts},
          open(f"{work}/index.json", "w"), indent=1)
print(f"map {W}x{H}, {len(placed)} placed nodes;", json.dumps(counts))
