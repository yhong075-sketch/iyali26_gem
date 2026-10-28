"""Build the UM7 nuclear-diffusion Escher map and its scenario data (no LP: reads flux.json and built_checks.json).

Usage: python build.py FLUX.json BUILT_CHECKS.json WORKDIR
BUILT_CHECKS.json is ../../built_checks.json (whole-model FVA blocked sets of both models, open and Task A). Writes
WORKDIR/maps/um7_nuclear.json and _deco.json, WORKDIR/data/<scenario>.json and WORKDIR/index.json. Three sections:
the nuclear NAD block with the two new links, its cytosolic partners, and the homocitrate / acetyl-CoA branch,
greyed because it stays blocked. Adapted from the dolichol-route build.py (reaction set, scenarios, status source).
"""
import json
import os
import sys

from maplib import build_map

flux_path, checks_path, work = sys.argv[1:]
F = json.load(open(flux_path))
C = json.load(open(checks_path))
BEFORE, AFTER = F["models"]["before"], F["models"]["after"]
RX, MT, SEC = F["reactions"], F["metabolites"], F["sections"]
NEW = {"R_ADPRtn", "R_PPItn"}
SECTIONS = [
    ("Nucleus: NAD synthesis (R569 QNS1, R579 NMA1), the sirtuin step R567, the inherited nucleus-cytosol transports and the "
     "two new diffusions R_ADPRtn (ADP-ribose) and R_PPItn (diphosphate)", SEC["nucleus"]),
    ("Cytosol: NAD synthesis (R580, R568, R578), ADP-ribose pyrophosphatase R1912 (YSA1), inorganic pyrophosphatase R419 "
     "(IPP1), nicotinamidase R576", SEC["cytosol"]),
    ("GREYED: homocitrate synthase R396 and its acetyl-CoA / CoA / 2-oxoglutarate / acetate transports, R92 closed; "
     "blocked before and after, not changed", SEC["homocitrate"]),
]
mapped = [rid for _, ids in SECTIONS for rid in ids]
BLK = {(side, cond): set(C["baseline" if side == "before" else "candidate"][f"blocked_{cond}"])
       for side in ("before", "after") for cond in ("open", "taskA")}


def capacity(side, rid):
    pair = F["models"][side]["fva_capacity_taskA"].get(rid)
    return None if pair is None else clean(max(abs(pair[0]), abs(pair[1])))


def clean(v):
    return 0.0 if abs(v) < 1e-12 else v


def status(rid):
    if rid in NEW:
        return 4, "NEW in the candidate, gene-less diffusion"
    if AFTER["bounds"][rid] == [0.0, 0.0]:
        return 3, "CLOSED (0,0) in both"
    b, a = rid in BLK[("before", "taskA")], rid in BLK[("after", "taskA")]
    if b and not a:
        return 1, "FREED (Task A)"
    if b and a:
        return 2, "blocked in both (Task A)"
    return 5, "unblocked in both"


rxns, status_values = {}, {}
for rid in mapped:
    r = dict(RX[rid])
    status_values[rid], r["sublabel"] = status(rid)
    rxns[rid] = r


def dead(side, mid):
    """Every reaction of the metabolite is blocked on Task A (whole-model FVA at 0 %) in that model."""
    reactions = [x for x in MT[mid]["reactions"] if not (side == "before" and x in NEW)]
    return all(x in BLK[(side, "taskA")] for x in reactions)


mets = {mid: dict(MT[mid], gap=dead("before", mid), cls="") for r in rxns.values() for mid in r["mets"]}
met_data = {mid: (1 if not dead("after", mid) else 2) if dead("before", mid) else (3 if dead("after", mid) else 0)
            for mid in mets}
legend = [
    "Reaction colour = scenario value: capacity = max |flux| any feasible state can carry (FVA at 0 % on Task A + PO1f, uracil "
    "open, R612 closed; Gurobi 1e-9), pFBA = flux at the optimum. Grey = 0; no colour = not in this model.",
    "Node colour (Task A): blue = every reaction blocked before, not after; orange = every reaction blocked in both; pink = "
    "only after; grey = has an unblocked reaction. Repeated nodes are the same pool. Layout is automatic.",
]
title = "UM7 nuclear NAD block: 9d655c2b (9d655c2bdbf82721...) -> + R_ADPRtn, R_PPItn"
mj, dj, _, W, H, placed = build_map("iyali26-um7-nuclear", title, rxns, mets, SECTIONS, maxw=6000, legend=legend, dx=600)
dj["rects"][2].update(fill="#e9e9e9", stroke="#a9a9a9")

scenarios = {
    "capacity_before": {rid: capacity("before", rid) for rid in mapped if rid not in BEFORE["absent"]},
    "capacity_after": {rid: capacity("after", rid) for rid in mapped},
    "before": {rid: clean(BEFORE["pfba"]["fluxes"][rid]) for rid in mapped if rid not in BEFORE["absent"]},
    "after": {rid: clean(AFTER["pfba"]["fluxes"][rid]) for rid in mapped},
    "status": status_values,
}


def describe(side, rid):
    m = F["models"][side]
    if rid in m["absent"]:
        return "absent"
    f = lambda key: "[" + ", ".join(f"{clean(v):.3g}" for v in m[key][rid]) + "]"
    return (f"pFBA {clean(m['pfba']['fluxes'][rid]):.3g}, FVA at optimum {f('fva_at_optimum')}, capacity Task A "
            f"{f('fva_capacity_taskA')}, open {f('fva_capacity_open')}, bounds {m['bounds'][rid]}")


tips = {}
full = {1: "blocked on Task A in 9d655c2b, freed in the candidate", 2: "blocked on Task A in both models",
        3: "closed (0, 0) in both models", 4: "absent in 9d655c2b; added, gene-less", 5: "unblocked in both models"}
for rid in mapped:
    blk = lambda side, c: "absent" if (side == "before" and rid in NEW) else ("blocked" if rid in BLK[(side, c)] else "unblocked")
    tips[rid] = (f"{full[status_values[rid]]} · 9d655c2b: {describe('before', rid)} · candidate: {describe('after', rid)} · "
                 f"open condition {blk('before', 'open')} to {blk('after', 'open')}, Task A {blk('before', 'taskA')} to "
                 f"{blk('after', 'taskA')}")
for mid in mets:
    tips[mid] = (f"Task A: every reaction blocked in 9d655c2b {'yes' if dead('before', mid) else 'no'}, in the candidate "
                 f"{'yes' if dead('after', mid) else 'no'} · {MT[mid]['degree']} reactions in the candidate")

os.makedirs(f"{work}/maps", exist_ok=True)
os.makedirs(f"{work}/data", exist_ok=True)
json.dump(mj, open(f"{work}/maps/um7_nuclear.json", "w"))
json.dump(dj, open(f"{work}/maps/um7_nuclear_deco.json", "w"))
for name, data in scenarios.items():
    json.dump(data, open(f"{work}/data/{name}.json", "w"))
counts = {
    "mapped_reactions": len(mapped),
    "in_model": {"before": len(mapped) - len(BEFORE["absent"]), "after": len(mapped)},
    "blocked_taskA": {side: sum(1 for rid in mapped if rid in BLK[(side, "taskA")]) for side in ("before", "after")},
    "blocked_open": {side: sum(1 for rid in mapped if rid in BLK[(side, "open")]) for side in ("before", "after")},
    "nonzero_capacity_taskA": {side: sum(1 for rid in mapped if (scenarios[f"capacity_{side}"].get(rid) or 0) != 0)
                               for side in ("before", "after")},
    "nonzero_pfba": {side: sum(1 for rid in mapped if (scenarios[side].get(rid) or 0) != 0) for side in ("before", "after")},
    "status": {k: sorted(rid for rid in mapped if status_values[rid] == v)
               for k, v in (("freed", 1), ("blocked_both", 2), ("closed_both", 3), ("new", 4), ("unblocked_both", 5))},
    "metabolite_nodes": {"dead_before_not_after": sorted(m for m, v in met_data.items() if v == 1),
                         "dead_both": sorted(m for m, v in met_data.items() if v == 2),
                         "dead_after_only": sorted(m for m, v in met_data.items() if v == 3),
                         "alive_both": sum(1 for v in met_data.values() if v == 0)},
}
json.dump({"title": title, "W": W, "H": H, "met_data": met_data, "tips": tips,
           "sections": [(lab, len(ids)) for lab, ids in SECTIONS], "counts": counts},
          open(f"{work}/index.json", "w"), indent=1)
print(f"map {W}x{H}, {len(placed)} placed nodes;", json.dumps(counts))
