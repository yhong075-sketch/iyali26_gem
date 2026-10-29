"""Build the UM11 nuclear-RNR Escher map and its scenario data (no LP: reads flux.json).

Usage: python build.py FLUX.json WORKDIR
Writes WORKDIR/maps/um11_nuclear_rnr.json and _deco.json, WORKDIR/data/<scenario>.json and WORKDIR/index.json. Two
sections: the nucleus (the 11 UM11 reactions, the active dUDP exit R1104 / R789 and the two what-if thioredoxin links,
which exist only in the what-if) and the cytosolic counterparts (RNR copies, TRR R755, UMP kinase R781, the NTP
reductases R702-R705). Adapted from the UM2 figure's build.py (reaction set, scenarios, what-if-only reactions).
"""
import json
import os
import sys

from maplib import build_map

flux_path, work = sys.argv[1:]
F = json.load(open(flux_path))
V = F["variants"]
RX, MT, SEC = F["reactions"], F["metabolites"], F["sections"]
UM11 = ["R707", "R709", "R711", "R1914", "R782", "R1033", "R1076", "R1077", "R1093", "R1136", "R1353"]
WHATIF = ["T_TRX1tn", "T_TRXSStn"]
# the census roots of UM11 (root_metabolites.csv): reduced TRX1 [nu] has no producer (RNP), the disulphide [nu] no consumer (RNC)
ROOTS = {"m931[C_nu]": "census root RNP: no producer in the nucleus", "m933[C_nu]": "census root RNC: no consumer in the nucleus"}
SECTIONS = [
    ("Nucleus: RNR copies R707 / R709 / R711 / R1914, UMP kinase R782 and the nucleotide diffusions (UM11, all blocked); "
     "dUDP exit R1104 / R789 active; purple = the what-if thioredoxin links (Option A, NOT in the model)", SEC["nucleus"]),
    ("Cytosol: the active RNR copies (same GPRs), thioredoxin reductase R755, UMP kinase R781 and the NTP reductases "
     "R702-R705 (thioredoxin gene only)", SEC["cytosol"]),
]
mapped = [rid for _, ids in SECTIONS for rid in ids]
BLK = {k: set(v["blocked_taskA"]) for k, v in V.items()}


def clean(v):
    return 0.0 if abs(v) < 1e-12 else v


def maxabs(pair):
    return clean(max(abs(pair[0]), abs(pair[1])))


def status(rid):
    if rid in WHATIF:
        return 4, "WHAT-IF ONLY (Option A, not implemented)"
    if V["as_built"]["bounds"][rid] == [0.0, 0.0]:
        return 3, "CLOSED (0,0)"
    if rid in BLK["as_built"]:
        return 2, "blocked, kept (Option B)" if rid in UM11 else "open but blocked (Task A)"
    return 5, "can carry flux (Task A)"


rxns, status_values = {}, {}
for rid in mapped:
    r = dict(RX[rid])
    status_values[rid], r["sublabel"] = status(rid)
    rxns[rid] = r


def dead(variant, mid):
    return all(x in BLK[variant] for x in MT[mid]["reactions"][variant])


mets = {mid: dict(name=MT[mid]["name"], comp=MT[mid]["comp"], degree=MT[mid]["degree"].get("as_built", MT[mid]["degree"]["whatif_A"]),
                  gap=dead("as_built", mid), cls="") for r in rxns.values() for mid in r["mets"]}
met_data = {mid: 3 if mid in ROOTS else (2 if dead("whatif_A", mid) else 1) if dead("as_built", mid) else 0 for mid in mets}
legend = [
    "Nothing in the model was changed (user decision B). Reaction colour = scenario value: capacity = max |flux| any feasible "
    "state can carry (FVA at 0 % on Task A + PO1f, uracil open, R612 closed; Gurobi 1e-9); pFBA = flux at the optimum. "
    "The what-if is an in-memory test that was NOT implemented.",
    "Node colour (Task A, whole-model FVA): dark red, large = the two census roots, TRX1 [nu] (RNP, no producer) and "
    "TRX1-disulphide [nu] (RNC, no consumer); blue = every reaction blocked as built, not in the what-if; orange = blocked in "
    "both; grey = has an unblocked reaction. Repeated nodes are the same pool. Layout is automatic.",
]
title = f"UM11 nuclear RNR copies kept blocked (B): {F['model_sha256'][:16]}... as built, with the Option A what-if"
mj, dj, _, W, H, placed = build_map("iyali26-um11-nuclear-rnr", title, rxns, mets, SECTIONS, maxw=5200, legend=legend, dx=560)


def scen(variant, key):
    return {rid: (maxabs(V[variant][key][rid]) if key != "pfba" else clean(V[variant]["pfba"]["fluxes"][rid]))
            for rid in mapped if rid in V[variant]["present"]}


scenarios = {
    "status": status_values,
    "pfba_as_built": scen("as_built", "pfba"),
    "capacity_as_built": scen("as_built", "capacity_taskA"),
    "capacity_whatif_A": scen("whatif_A", "capacity_taskA"),
    "pfba_whatif_A": scen("whatif_A", "pfba"),
}


def describe(variant, rid):
    m = V[variant]
    if rid not in m["present"]:
        return "not in this variant"
    f = lambda key: "[" + ", ".join(f"{clean(v):.3g}" for v in m[key][rid]) + "]"
    return (f"bounds {m['bounds'][rid]}, pFBA {clean(m['pfba']['fluxes'][rid]):.3g}, capacity {f('capacity_taskA')}, "
            f"no-uptake {f('no_uptake')}")


tips = {}
for rid in mapped:
    tips[rid] = (f"{status(rid)[1]} · GPR: {RX[rid]['gpr'] or 'none'} · as built: {describe('as_built', rid)} · "
                 f"what-if A: {describe('whatif_A', rid)}")
for mid in mets:
    tips[mid] = ((ROOTS[mid] + " · ") if mid in ROOTS else "") + (f"Task A: every reaction blocked as built {'yes' if dead('as_built', mid) else 'no'}, in what-if A "
                 f"{'yes' if dead('whatif_A', mid) else 'no'} · {mets[mid]['degree']} reactions as built")
os.makedirs(f"{work}/maps", exist_ok=True)
os.makedirs(f"{work}/data", exist_ok=True)
json.dump(mj, open(f"{work}/maps/um11_nuclear_rnr.json", "w"))
json.dump(dj, open(f"{work}/maps/um11_nuclear_rnr_deco.json", "w"))
for name, data in scenarios.items():
    json.dump(data, open(f"{work}/data/{name}.json", "w"))
nz = lambda s: sorted(rid for rid in mapped if (scenarios[s].get(rid) or 0) != 0)
counts = {
    "mapped_reactions": len(mapped),
    "um11_blocked_as_built": sorted(r for r in UM11 if r in BLK["as_built"]),
    "um11_blocked_whatif_A": sorted(r for r in UM11 if r in BLK["whatif_A"]),
    "can_carry_flux_as_built": sorted(rid for rid in mapped if status_values[rid] == 5),
    "nonzero_capacity": {s: nz(s) for s in ("capacity_as_built", "capacity_whatif_A")},
    "nonzero_pfba": {s: nz(s) for s in ("pfba_as_built", "pfba_whatif_A")},
    "nonzero_no_uptake": {k: sorted(r for r in v["present"] if maxabs(v["no_uptake"][r]) != 0) for k, v in V.items()},
    "whole_model_blocked_taskA": {k: len(v) for k, v in BLK.items()},
    "freed_whole_model_taskA": sorted(BLK["as_built"] - BLK["whatif_A"]),
    "newly_blocked_whole_model_taskA": sorted(BLK["whatif_A"] - BLK["as_built"]),
    "metabolite_nodes": {"roots": sorted(m for m, v in met_data.items() if v == 3),
                         "dead_both": sorted(m for m, v in met_data.items() if v == 2),
                         "dead_as_built_only": sorted(m for m, v in met_data.items() if v == 1),
                         "alive": sum(1 for v in met_data.values() if v == 0)},
    "wt": {k: v["wt"]["growth"] for k, v in V.items()},
}
json.dump({"title": title, "W": W, "H": H, "met_data": met_data, "tips": tips,
           "sections": [(lab, len(ids)) for lab, ids in SECTIONS], "counts": counts},
          open(f"{work}/index.json", "w"), indent=1)
print(f"map {W}x{H}, {len(placed)} placed nodes;", json.dumps(counts))
