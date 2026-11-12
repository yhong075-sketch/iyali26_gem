"""Build the R2121 closure Escher map and its scenario data (no LP: reads flux.json).

Usage: python build.py FLUX.json WORKDIR
Writes WORKDIR/maps/r2121.json and _deco.json, WORKDIR/data/<scenario>.json and WORKDIR/index.json. Three sections:
the cytosolic malonyl-CoA / acetyl-CoA / bicarbonate supply, the peroxisomal beta-oxidation with the carnitine export,
and the ER elongation with TAG lipolysis and synthesis. Adapted from the UM7 nuclear-diffusion build.py (sections,
scenarios: pFBA before / after and an on-off status instead of the census blocked status).
"""
import json
import os
import sys

from maplib import build_map

flux_path, work = sys.argv[1:]
F = json.load(open(flux_path))
BEFORE, AFTER = F["models"]["before"], F["models"]["after"]
RX, MT, SEC = F["reactions"], F["metabolites"], F["sections"]
ON = 1e-9
SECTIONS = [
    ("Cytosol: malonyl-CoA supply. R2121 (KEGG R00740, closed after) versus acetyl-CoA carboxylase R88 with its "
     "bicarbonate (R1026 / R2202, carbonic anhydrase gene); Asp-derived semialdehyde R1936 / R2120 / R2119; FAS R1392 / "
     "R1393; ATP-citrate lyase R1894; acetyl-CoA synthetase R90; GlcNAc-6-P synthase R561", SEC["supply"]),
    ("Peroxisome: beta-oxidation of oleate to laurate (R1496 ... R84, thioesterase R619), catalase R195 and the "
     "acetylcarnitine export R192 / R1032 / R191", SEC["peroxisome"]),
    ("Lipid particle and ER membrane: TAG lipase R1715, laurate activation R1445, elongation to palmitoyl-CoA "
     "(R1406 ... R1428) with ER malonyl-CoA, DGAT160_lp", SEC["elongation_tag"]),
]
mapped = [rid for _, ids in SECTIONS for rid in ids]


def clean(v):
    # Same threshold as the status (ON): R2119 carries -1.7e-12 of solver noise in the e584151f pFBA.
    return 0.0 if abs(v) <= ON else v


def status(rid):
    b, a = abs(BEFORE["pfba"]["fluxes"][rid]) > ON, abs(AFTER["pfba"]["fluxes"][rid]) > ON
    if AFTER["bounds"][rid] == [0.0, 0.0]:
        return 5, "CLOSED (0, 0) after"
    if b and not a:
        return 1, "carries flux before only"
    if a and not b:
        return 2, "carries flux after only"
    if a and b:
        return 3, "carries flux in both"
    return 4, "no flux in either"


rxns, status_values = {}, {}
for rid in mapped:
    r = dict(RX[rid])
    status_values[rid], r["sublabel"] = status(rid)
    rxns[rid] = r
mets = {mid: dict(MT[mid], gap=False, cls="") for r in rxns.values() for mid in r["mets"]}
legend = [
    "Reaction colour = pFBA flux at the optimum (Task A medium + PO1f, Gurobi 1e-9); grey = 0. Status scenario: "
    "blue = flux before only, purple = after only, light blue = both, grey = neither, black = closed after.",
    "Layout is automatic; repeated nodes are the same pool. Hub metabolites are repeated per reaction.",
]
title = "R2121 closed: e584151f (e584151fc794...) -> candidate (R2121 bounds (0, 0))"
mj, dj, _, W, H, placed = build_map("iyali26-r2121-close", title, rxns, mets, SECTIONS, maxw=7000, legend=legend, dx=560)
met_data = {mid: 0 for mid in mets}
scenarios = {
    "before": {rid: clean(BEFORE["pfba"]["fluxes"][rid]) for rid in mapped},
    "after": {rid: clean(AFTER["pfba"]["fluxes"][rid]) for rid in mapped},
    "status": status_values,
}


def describe(side, rid):
    m = F["models"][side]
    text = f"pFBA {clean(m['pfba']['fluxes'][rid]):.3g}, bounds {m['bounds'][rid]}"
    if rid in m["fva_999"]:
        text += ", FVA at 99.9 % [" + ", ".join(f"{clean(v):.3g}" for v in m["fva_999"][rid]) + "]"
    return text


tips = {rid: f"{rxns[rid]['sublabel']} · e584151f: {describe('before', rid)} · candidate: {describe('after', rid)}"
        for rid in mapped}
for mid in mets:
    tips[mid] = f"{MT[mid]['name']} [{MT[mid]['comp']}] · {MT[mid]['degree']} reactions in the candidate"

os.makedirs(f"{work}/maps", exist_ok=True)
os.makedirs(f"{work}/data", exist_ok=True)
json.dump(mj, open(f"{work}/maps/r2121.json", "w"))
json.dump(dj, open(f"{work}/maps/r2121_deco.json", "w"))
for name, data in scenarios.items():
    json.dump(data, open(f"{work}/data/{name}.json", "w"))
counts = {
    "mapped_reactions": len(mapped),
    "nonzero_pfba": {side: sum(1 for rid in mapped if scenarios[side][rid] != 0) for side in ("before", "after")},
    "status": {k: sorted(rid for rid in mapped if status_values[rid] == v)
               for k, v in (("before_only", 1), ("after_only", 2), ("both", 3), ("neither", 4), ("closed_after", 5))},
    "pfba_total_abs_flux": {side: F["models"][side]["pfba"]["total_abs_flux"] for side in ("before", "after")},
    "wt": {side: F["models"][side]["wt"]["growth"] for side in ("before", "after")},
}
json.dump({"title": title, "W": W, "H": H, "met_data": met_data, "tips": tips,
           "sections": [(lab, len(ids)) for lab, ids in SECTIONS], "counts": counts},
          open(f"{work}/index.json", "w"), indent=1)
print(f"map {W}x{H}, {len(placed)} placed nodes;", json.dumps(counts))
