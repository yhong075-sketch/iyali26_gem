"""Build the cofactor-terms Escher map and its scenario data (no LP: reads flux.json and the two open census runs).

Usage: python build.py FLUX.json BASE_open_metabolites.csv NEW_open_metabolites.csv WORKDIR
Writes WORKDIR/maps/cofactor_terms.json and _deco.json, WORKDIR/data/<scenario>.json and WORKDIR/index.json.
Mapped reactions: the union of the reactions freed in the open and the Task A census comparisons
(a1237821 -> 00e8e397), grouped by cofactor family, plus biomass_C drawn with its six new terms only.
"""
import csv
import json
import os
import sys

from maplib import build_map

flux_path, base_mets, new_mets, work = sys.argv[1:]
F = json.load(open(flux_path))
BEFORE, AFTER = F["models"]["a1237821"], F["models"]["00e8e397"]
RX, MT = F["reactions"], F["metabolites"]
FREED_OPEN, FREED_TASKA = set(F["freed_open"]), set(F["freed_taskA"])
TERMS = ["m479[C_cy]", "m528[C_cy]", "m926[C_cy]", "m122[C_cy]", "m41[C_cy]", "m69[C_cy]"]
FAMILIES = [
    ("Folate / THF (term THF m479, 6.34e-5)", ["R377", "R1097", "R275", "R137", "R274", "R1094", "R14", "R973", "R277", "R278",
                                              "R1098", "R273", "R1145", "R1146", "R1147", "R64", "R811", "R990", "R991"]),
    ("Riboflavin, FMN, FAD (terms riboflavin m926, 9.9e-4; FAD m528, 1e-5)", ["R378", "R1307", "R12", "R11", "R34", "R699", "R700",
                                                                               "R1314", "R1316", "R697", "R306"]),
    ("NAD / NADP (terms NAD m122, 2.65e-3; NADP m41, 5.7e-4)", ["R580", "R568", "R578", "R565", "R571"]),
    ("Coenzyme A (term CoA m69, 1.9e-4)", ["R616", "R647", "R648", "R615", "R1092"]),
    ("biomass_C: the six terms added in 00e8e397 (absent in a1237821)", ["biomass_C"]),
]
mapped = [rid for _, ids in FAMILIES for rid in ids]
if sorted(mapped) != sorted(F["mapped"]):
    sys.exit(f"STOP: family lists do not match the mapped set: {sorted(set(mapped) ^ set(F['mapped']))}")


def gaps(path):
    return {row["id"] for row in csv.DictReader(open(path)) if row["gap"] == "True"}


GAP_BEFORE, GAP_AFTER = gaps(base_mets), gaps(new_mets)


def freed_in(rid):
    return " + ".join(c for c, s in (("open", FREED_OPEN), ("Task A", FREED_TASKA)) if rid in s)


def clean(v):
    return 0.0 if abs(v) < 1e-12 else v


rxns = {}
for rid in mapped:
    r = dict(RX[rid])
    if rid == "biomass_C":
        r["mets"] = {m: RX[rid]["mets"][m] for m in TERMS}
        r["pseudo_name"] = "(biomass)"
        r["sublabel"] = "six cofactor terms only (55 terms in all)"
    else:
        r["sublabel"] = f"freed in {freed_in(rid)} census"
    rxns[rid] = r
mets = {mid: dict(MT[mid], gap=mid in GAP_BEFORE, cls="") for r in rxns.values() for mid in r["mets"]}
legend = [
    "Reaction colour and label = pFBA flux (mmol/gDW/h) on Task A medium + PO1f, Gurobi 1e-9; grey = 0. "
    "Scenario 'census': amber = freed in the open condition only, pink = freed on Task A only, blue = freed in both.",
    "Node colour (open census): blue = gap metabolite under a1237821 and not under 00e8e397; orange = gap under both; "
    "grey = not a gap. Repeated nodes are the same pool. Layout is automatic and carries no biological meaning.",
]
title = ("Cofactor biomass terms: reactions freed by a1237821 -> 00e8e397 (a1237821 = a12378210c7c6672..., "
         "00e8e397 = 00e8e39775b64cfe...)")
mj, dj, _, W, H, placed = build_map("iyali26-cofactor-terms", title, rxns, mets, FAMILIES, maxw=7000, legend=legend, dx=520)
met_data = {}
for mid in mets:
    met_data[mid] = (1 if mid not in GAP_AFTER else 2) if mid in GAP_BEFORE else (3 if mid in GAP_AFTER else 0)

# biomass_C is drawn with the six terms of 00e8e397, which a1237821 lacks, so it gets no colour in "before".
scenarios = {
    "before": {rid: clean(BEFORE["pfba"]["fluxes"][rid]) for rid in mapped if rid != "biomass_C"},
    "after": {rid: clean(AFTER["pfba"]["fluxes"][rid]) for rid in mapped},
    "census": {rid: (3 if rid in FREED_OPEN and rid in FREED_TASKA else 1 if rid in FREED_OPEN else 2)
               for rid in mapped if rid != "biomass_C"},
}
tips = {}
for rid in mapped:
    fb, fa = BEFORE["fva_at_optimum"][rid], AFTER["fva_at_optimum"][rid]
    tips[rid] = (("biomass reaction" if rid == "biomass_C" else f"freed in the {freed_in(rid)} census")
                 + f" · pFBA a1237821 {clean(BEFORE['pfba']['fluxes'][rid]):.3g} -> 00e8e397 {scenarios['after'][rid]:.3g}"
                 + f" · FVA at the optimum a1237821 [{clean(fb[0]):.3g}, {clean(fb[1]):.3g}]"
                 + f" -> 00e8e397 [{clean(fa[0]):.3g}, {clean(fa[1]):.3g}]")
for mid in mets:
    tips[mid] = (f"open census: gap under a1237821 {'yes' if mid in GAP_BEFORE else 'no'}, "
                 f"under 00e8e397 {'yes' if mid in GAP_AFTER else 'no'} · {MT[mid]['degree']} reactions in the model")

os.makedirs(f"{work}/maps", exist_ok=True)
os.makedirs(f"{work}/data", exist_ok=True)
json.dump(mj, open(f"{work}/maps/cofactor_terms.json", "w"))
json.dump(dj, open(f"{work}/maps/cofactor_terms_deco.json", "w"))
for name, data in scenarios.items():
    json.dump(data, open(f"{work}/data/{name}.json", "w"))
counts = {
    "mapped_reactions": len(mapped) - 1,
    "nonzero_pfba_before": sum(1 for rid in mapped if rid != "biomass_C" and scenarios["before"][rid] != 0),
    "nonzero_pfba_after": sum(1 for rid in mapped if rid != "biomass_C" and scenarios["after"][rid] != 0),
    "fva_can_carry_after": sum(1 for rid in mapped if rid != "biomass_C" and max(map(abs, AFTER["fva_at_optimum"][rid])) > 1e-12),
    "freed_open_only": sorted(FREED_OPEN - FREED_TASKA), "freed_taskA_only": sorted(FREED_TASKA - FREED_OPEN),
    "metabolite_nodes": {"gap_before_not_after": sum(1 for v in met_data.values() if v == 1),
                         "gap_both": sum(1 for v in met_data.values() if v == 2),
                         "gap_after_only": sum(1 for v in met_data.values() if v == 3),
                         "not_gap": sum(1 for v in met_data.values() if v == 0)},
}
json.dump({"title": title, "W": W, "H": H, "met_data": met_data, "tips": tips, "sections": [(lab, len(ids)) for lab, ids in FAMILIES],
           "counts": counts}, open(f"{work}/index.json", "w"), indent=1)
print(f"map {W}x{H}, {len(placed)} placed nodes;", json.dumps(counts))
