"""Build the UM2 vacuolar-closure Escher map and its scenario data (no LP: reads flux.json).

Usage: python build.py FLUX.json WORKDIR
Writes WORKDIR/maps/um2_vacuole.json and _deco.json, WORKDIR/data/<scenario>.json and WORKDIR/index.json. Four
sections: pump / H+ / water, amino-acid transports (+ the open His/Lys imports R2055/R2056 for context),
polyamines, the four dipeptide chains. Adapted from the UM7 figure's build.py (reaction set, scenarios, status from
the whole-model FVA in flux.json).
"""
import json
import os
import sys

from maplib import build_map

flux_path, work = sys.argv[1:]
F = json.load(open(flux_path))
V = F["variants"]
RX, MT, SEC = F["reactions"], F["metabolites"], F["sections"]
SECTIONS = [
    ("Tonoplast energetics: V-ATPase R795, H+ diffusion R1163 (written cyt -> vac), water diffusion R1363", SEC["pump"]),
    ("Vacuolar amino-acid transports: AVT1-like import (H+ antiport), AVT3 / ATG22-like and AVT6-like export (H+ symport); "
     "R2055 / R2056 His / Lys import are open context", SEC["amino"]),
    ("Polyamine import (TPO-like; no vacuolar exit, R925 unbalanced)", SEC["polyamine"]),
    ("Dipeptide chains added by iYli21: secretion-only exchange, plasma-membrane and tonoplast transport, gene-less "
     "vacuolar hydrolysis, amino-acid exit", SEC["dipeptide"]),
]
mapped = [rid for _, ids in SECTIONS for rid in ids]
BLK = {k: set(v["blocked_taskA"]) for k, v in V.items()}


def clean(v):
    return 0.0 if abs(v) < 1e-12 else v


def maxabs(pair):
    return clean(max(abs(pair[0]), abs(pair[1])))


def status(rid):
    if V["as_built"]["bounds"][rid] == [0.0, 0.0]:
        return 3, "CLOSED (0,0), inherited; kept (A0)"
    if rid in BLK["as_built"]:
        return 2, "open but blocked (Task A)"
    return 5, "can carry flux (Task A)"


rxns, status_values = {}, {}
for rid in mapped:
    r = dict(RX[rid])
    status_values[rid], r["sublabel"] = status(rid)
    rxns[rid] = r


def dead(variant, mid):
    return all(x in BLK[variant] for x in MT[mid]["reactions"])


mets = {mid: dict(MT[mid], gap=dead("as_built", mid), cls="") for r in rxns.values() for mid in r["mets"]}
met_data = {mid: (2 if dead("cprime", mid) else 1) if dead("as_built", mid) else 0 for mid in mets}
legend = [
    "Nothing in the model was changed (user decision A0). Reaction colour = scenario value: capacity = max |flux| any "
    "feasible state can carry (FVA at 0 % on Task A + PO1f, uracil open, R612 closed; Gurobi 1e-9). What-if scenarios "
    "are in-memory tests that were NOT implemented.",
    "Node colour (Task A, whole-model FVA): orange = every reaction blocked as built and in the C' what-if; blue = every "
    "reaction blocked as built, not in the what-if; grey = has an unblocked reaction. Repeated nodes are the same pool. "
    "Layout is automatic.",
]
title = f"UM2 vacuolar closures kept (A0): {F['model_sha256'][:16]}... as built, with what-if tests"
mj, dj, _, W, H, placed = build_map("iyali26-um2-vacuole", title, rxns, mets, SECTIONS, maxw=6000, legend=legend, dx=600)
dj["rects"][3].update(fill="#efefef", stroke="#a9a9a9")
scenarios = {
    "status": status_values,
    "capacity_as_built": {rid: maxabs(V["as_built"]["capacity_taskA"][rid]) for rid in mapped},
    "capacity_whatif_cprime": {rid: maxabs(V["cprime"]["capacity_taskA"][rid]) for rid in mapped},
    "loop_whatif_open17": {rid: maxabs(V["open17_as_written"]["no_uptake"][rid]) for rid in mapped},
}


def describe(variant, rid):
    m = V[variant]
    f = lambda key: "[" + ", ".join(f"{clean(v):.3g}" for v in m[key][rid]) + "]"
    return (f"bounds {m['bounds'][rid]}, pFBA {clean(m['pfba']['fluxes'][rid]):.3g}, capacity {f('capacity_taskA')}, "
            f"no-uptake {f('no_uptake')}")


tips = {}
for rid in mapped:
    tips[rid] = (f"{status(rid)[1]} · GPR: {RX[rid]['gpr'] or 'none'} · as built: {describe('as_built', rid)} · what-if C': "
                 f"{describe('cprime', rid)} · what-if open 17 as written: {describe('open17_as_written', rid)}")
for mid in mets:
    tips[mid] = (f"Task A: every reaction blocked as built {'yes' if dead('as_built', mid) else 'no'}, in what-if C' "
                 f"{'yes' if dead('cprime', mid) else 'no'} · {MT[mid]['degree']} reactions")
os.makedirs(f"{work}/maps", exist_ok=True)
os.makedirs(f"{work}/data", exist_ok=True)
json.dump(mj, open(f"{work}/maps/um2_vacuole.json", "w"))
json.dump(dj, open(f"{work}/maps/um2_vacuole_deco.json", "w"))
for name, data in scenarios.items():
    json.dump(data, open(f"{work}/data/{name}.json", "w"))
nz = lambda s: sorted(rid for rid in mapped if (scenarios[s].get(rid) or 0) != 0)
counts = {
    "mapped_reactions": len(mapped),
    "closed_as_built": sorted(rid for rid in mapped if status_values[rid] == 3),
    "open_but_blocked_as_built": sorted(rid for rid in mapped if status_values[rid] == 2),
    "can_carry_flux_as_built": sorted(rid for rid in mapped if status_values[rid] == 5),
    "nonzero_capacity": {s: nz(s) for s in ("capacity_as_built", "capacity_whatif_cprime")},
    "nonzero_no_uptake_open17": nz("loop_whatif_open17"),
    "nonzero_pfba": {k: sorted(rid for rid in mapped if abs(V[k]["pfba"]["fluxes"][rid]) > 1e-12) for k in V},
    "metabolite_nodes": {"dead_both": sorted(m for m, v in met_data.items() if v == 2),
                         "dead_as_built_only": sorted(m for m, v in met_data.items() if v == 1),
                         "alive": sum(1 for v in met_data.values() if v == 0)},
    "wt": {k: v["wt"]["growth"] for k, v in V.items()},
}
json.dump({"title": title, "W": W, "H": H, "met_data": met_data, "tips": tips,
           "sections": [(lab, len(ids)) for lab, ids in SECTIONS], "counts": counts},
          open(f"{work}/index.json", "w"), indent=1)
print(f"map {W}x{H}, {len(placed)} placed nodes;", json.dumps(counts))
