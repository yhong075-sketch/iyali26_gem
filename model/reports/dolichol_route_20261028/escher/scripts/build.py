"""Build the dolichol-route Escher map and its scenario data (no LP: reads flux.json and the open / Task A census runs).

Usage: python build.py FLUX.json CENSUS_DIR_BEFORE CENSUS_DIR_AFTER WORKDIR
CENSUS_DIR_BEFORE holds ce9e1112_{open,taskA_po1f}_{reactions,metabolites}.csv (Kennedy re-run census), CENSUS_DIR_AFTER
a801aad0_* (same optimization problem as 9d655c2b). Writes WORKDIR/maps/dolichol_route.json and _deco.json,
WORKDIR/data/<scenario>.json and WORKDIR/index.json. Three sections: the cytosolic route, the Dol-P / Dol-PP pools
and their consumers, and the lipid-particle Srt1 chain (census UM6), greyed because it stays blocked.
"""
import csv
import json
import os
import sys

from maplib import build_map

flux_path, cen_before, cen_after, work = sys.argv[1:]
F = json.load(open(flux_path))
BEFORE, AFTER = F["models"]["ce9e1112"], F["models"]["9d655c2b"]
RX, MT, SEC = F["reactions"], F["metabolites"], F["sections"]
TERM = "m503[C_er]"
SECTIONS = [
    ("Cytosol: GPP + IPP to FPP' (R2179) to dehydrodolichyl-PP (R2180) to dolichol (new R_DOLSYN; R2178 closed) to Dol-P (R287, SEC59)", SEC["route"]),
    ("Dol-P and Dol-PP pools (cytosol / ER) and their consumers; biomass_C drawn with the new Dol-P term only", SEC["pools"]),
    ("GREYED: lipid-particle Srt1 chain R217-R235 (census UM6), blocked before and after, not changed", SEC["um6"]),
]
mapped = [rid for _, ids in SECTIONS for rid in ids]


def census(tag_dir, tag, cond):
    rx = {r["id"]: r["blocked"] == "True" for r in csv.DictReader(open(f"{tag_dir}/{tag}_{cond}_reactions.csv"))}
    gp = {r["id"] for r in csv.DictReader(open(f"{tag_dir}/{tag}_{cond}_metabolites.csv")) if r["gap"] == "True"}
    return rx, gp


BLK = {("before", c): census(cen_before, "ce9e1112", c) for c in ("open", "taskA_po1f")}
BLK.update({("after", c): census(cen_after, "a801aad0", c) for c in ("open", "taskA_po1f")})
GAP_BEFORE, GAP_AFTER = BLK[("before", "open")][1], BLK[("after", "open")][1]


def status(rid):
    b = BLK[("before", "open")][0].get(rid)   # None = not in the model
    a = BLK[("after", "open")][0][rid]
    if rid == "R2178":
        return 3, "CLOSED (0,0) in 9d655c2b"
    if b is None:
        return 4, "NEW in 9d655c2b, gene-less"
    if b and not a:
        return 1, "FREED in 9d655c2b"
    if b and a:
        return 2, "blocked in both"
    return 5, "unblocked in both"


def clean(v):
    return 0.0 if abs(v) < 1e-12 else v


rxns, census_values, labels = {}, {}, {}
for rid in mapped:
    r = dict(RX[rid])
    value, label = status(rid)
    census_values[rid], labels[rid] = value, label
    if rid == "biomass_C":
        r["mets"] = {TERM: RX[rid]["mets"][TERM]}
        r["pseudo_name"] = "(biomass)"
        r["sublabel"] = f"Dol-P term only, {-RX[rid]['mets'][TERM]:g} mmol/gDW, new ({len(RX[rid]['mets'])} terms in all)"
    else:
        r["sublabel"] = label
    rxns[rid] = r
mets = {mid: dict(MT[mid], gap=mid in GAP_BEFORE, cls="") for r in rxns.values() for mid in r["mets"]}
legend = [
    "Reaction colour and label = pFBA flux (mmol/gDW/h) on Task A medium + PO1f (uracil open, R612 closed), Gurobi 1e-9; grey = 0; "
    "no colour = not in this model or (biomass_C before) without the drawn term. Scenario 'census': open-condition gap census.",
    "Node colour (open census): blue = gap metabolite under ce9e1112, not under 9d655c2b; orange = gap under both; pink = gap only "
    "under 9d655c2b; grey = not a gap. Repeated nodes are the same pool. Layout is automatic and carries no biological meaning.",
]
title = ("Dolichol route: ce9e1112 -> 9d655c2b (ce9e1112 = ce9e1112ee940531..., 9d655c2b = 9d655c2bdbf82721...)")
mj, dj, _, W, H, placed = build_map("iyali26-dolichol-route", title, rxns, mets, SECTIONS, maxw=6000, legend=legend, dx=640)
# Grey the UM6 box (third section rectangle).
dj["rects"][2].update(fill="#e9e9e9", stroke="#a9a9a9")
met_data = {}
for mid in mets:
    met_data[mid] = (1 if mid not in GAP_AFTER else 2) if mid in GAP_BEFORE else (3 if mid in GAP_AFTER else 0)

scenarios = {
    "before": {rid: clean(BEFORE["pfba"]["fluxes"][rid]) for rid in mapped if rid not in BEFORE["absent"] and rid != "biomass_C"},
    "after": {rid: clean(AFTER["pfba"]["fluxes"][rid]) for rid in mapped},
    "census": census_values,
}
tips = {}
for rid in mapped:
    fa = AFTER["fva_at_optimum"][rid]
    before = ("absent" if rid in BEFORE["absent"] else
              f"pFBA {clean(BEFORE['pfba']['fluxes'][rid]):.3g}, FVA [{clean(BEFORE['fva_at_optimum'][rid][0]):.3g}, "
              f"{clean(BEFORE['fva_at_optimum'][rid][1]):.3g}], bounds {BEFORE['bounds'][rid]}")
    blk = lambda side, c: {True: "blocked", False: "unblocked", None: "absent"}[BLK[(side, c)][0].get(rid)]
    full = {1: "blocked in ce9e1112, freed in 9d655c2b", 2: "blocked in both models (open census)", 3: "open (0, 1000) and blocked in ce9e1112; closed (0, 0) in 9d655c2b", 4: "absent in ce9e1112; added in 9d655c2b, gene-less and lumped", 5: "unblocked in both models"}[census_values[rid]]
    tips[rid] = (f"{full} · ce9e1112: {before} · 9d655c2b: pFBA {clean(AFTER['pfba']['fluxes'][rid]):.3g}, FVA "
                 f"[{clean(fa[0]):.3g}, {clean(fa[1]):.3g}], bounds {AFTER['bounds'][rid]} · census open {blk('before', 'open')} to "
                 f"{blk('after', 'open')}, Task A {blk('before', 'taskA_po1f')} to {blk('after', 'taskA_po1f')}")
for mid in mets:
    tips[mid] = (f"open census: gap under ce9e1112 {'yes' if mid in GAP_BEFORE else 'no'}, under 9d655c2b "
                 f"{'yes' if mid in GAP_AFTER else 'no'} · {MT[mid]['degree']} reactions in 9d655c2b")

os.makedirs(f"{work}/maps", exist_ok=True)
os.makedirs(f"{work}/data", exist_ok=True)
json.dump(mj, open(f"{work}/maps/dolichol_route.json", "w"))
json.dump(dj, open(f"{work}/maps/dolichol_route_deco.json", "w"))
for name, data in scenarios.items():
    json.dump(data, open(f"{work}/data/{name}.json", "w"))
nonbio = [rid for rid in mapped if rid != "biomass_C"]


def count(side, cond):
    rx = BLK[(side, cond)][0]
    return sum(1 for rid in nonbio if rx.get(rid) is True)


counts = {
    "mapped_reactions": len(nonbio), "plus": "biomass_C",
    "in_model": {"before": sum(1 for rid in nonbio if rid not in BEFORE["absent"]), "after": len(nonbio)},
    "blocked_open_census": {"before": count("before", "open"), "after": count("after", "open")},
    "blocked_taskA_census": {"before": count("before", "taskA_po1f"), "after": count("after", "taskA_po1f")},
    "nonzero_pfba": {"before": sum(1 for rid in nonbio if scenarios["before"].get(rid, 0) != 0),
                     "after": sum(1 for rid in nonbio if scenarios["after"][rid] != 0)},
    "fva_can_carry_at_optimum": {
        "before": sum(1 for rid in nonbio if rid not in BEFORE["absent"] and max(map(abs, BEFORE["fva_at_optimum"][rid])) > 1e-12),
        "after": sum(1 for rid in nonbio if max(map(abs, AFTER["fva_at_optimum"][rid])) > 1e-12)},
    "census_status": {k: sorted(rid for rid in mapped if census_values[rid] == v)
                      for k, v in (("freed", 1), ("blocked_both", 2), ("closed_after", 3), ("new", 4), ("unblocked_both", 5))},
    "metabolite_nodes": {"gap_before_not_after": sum(1 for v in met_data.values() if v == 1),
                         "gap_both": sum(1 for v in met_data.values() if v == 2),
                         "gap_after_only": sum(1 for v in met_data.values() if v == 3),
                         "not_gap": sum(1 for v in met_data.values() if v == 0)},
}
json.dump({"title": title, "W": W, "H": H, "met_data": met_data, "tips": tips, "sections": [(lab, len(ids)) for lab, ids in SECTIONS],
           "counts": counts}, open(f"{work}/index.json", "w"), indent=1)
print(f"map {W}x{H}, {len(placed)} placed nodes;", json.dumps(counts))
