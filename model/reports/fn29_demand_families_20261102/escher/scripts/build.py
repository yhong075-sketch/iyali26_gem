"""Build the five-family Escher map and its scenario data (no LP: reads flux.json).

Usage: python build.py FLUX.json WORKDIR
Writes WORKDIR/maps/fn29_families.json and _deco.json, WORKDIR/data/<scenario>.json and WORKDIR/index.json.
Five sections, one per family: the as-built reactions named in REPORT.md plus the reactions each family's what-if
adds (purple, not in the model). Every what-if is applied alone; the what-if scenarios show each section under its
own family's variant. Adapted from the UM11 figure's build.py (sections, statuses, scenarios, node classes).
"""
import json
import os
import sys

from maplib import build_map

flux_path, work = sys.argv[1:]
F = json.load(open(flux_path))
V = F["variants"]
RX, MT, SEC = F["reactions"], F["metabolites"], F["sections"]
WHATIF = {"sphingolipid": "sph_V2", "gpi": "gpi_b1", "cardiolipin": "cl_a", "nglycan": "ng_C2", "mito_translation": "mt_rec"}
ALT = {"nglycan": "ng_B"}
wt = lambda k: f"{V[k]['wt']:.6f}"
LABELS = {
    "sphingolipid": f"Sphingolipid (10 FN; AUR1 set R1825-R1834 inherited (0,0)). What-if V2: open R1825-R1834, 10 duplicate-IPC links, "
                    f"DAG gm->em link, biomass_C += 4.17e-4 IPC-C m1271[go]; WT {wt('sph_V2')}",
    "gpi": f"GPI anchor (6 FN; chain feasible, m603 drained only by R_DM_GPI). What-if b1: biomass_C += 1.70e-4 m603[er], exit open; "
           f"WT {wt('gpi_b1')}",
    "cardiolipin": f"Cardiolipin (2 FN; route unused, no cardiolipin term). What-if a: biomass_C += 7.16e-4 m1651[mm]; WT {wt('cl_a')}",
    "nglycan": f"N-glycan (3 FN; LLO chain split, mannan made by R289). What-if C2 = A2 + B (R289 1:1, Golgi mannan route, Golgi Pi exit, "
               f"R1137 closed, LLO link + OST, biomass_C += 0.01 N-glycan): WT {wt('ng_C2')}; B alone: WT {wt('ng_B')}",
    "mito_translation": f"Mitochondrial translation (8 FN; charged mt tRNAs have no consumer). What-if: translation pseudo-reaction (14 tRNAs), "
                        f"biomass_C += 0.0186 mtprot, R1251 and R1218 reversible; WT {wt('mt_rec')}",
}
ORDER = ["sphingolipid", "gpi", "cardiolipin", "nglycan", "mito_translation"]
added = {fam: sorted({r for k in (WHATIF[fam], ALT.get(fam)) if k for r in V[k]["added"]}) for fam in ORDER}
changed = {fam: sorted({r for k in (WHATIF[fam], ALT.get(fam)) if k for r in V[k]["changed"]}) for fam in ORDER}
SECTIONS = [(LABELS[fam], SEC[fam] + added[fam]) for fam in ORDER]
mapped = [rid for _, ids in SECTIONS for rid in ids]
fam_of = {rid: fam for fam in ORDER for rid in SEC[fam] + added[fam]}
AB = V["as_built"]
biomass_species = {mid: fam for fam in ORDER for k in (WHATIF[fam],) for mid in V[k]["applied"]["biomass_term"]}


def clean(v):
    return 0.0 if abs(v) < 1e-12 else v


def maxabs(pair):
    return clean(max(abs(pair[0]), abs(pair[1])))


def status(rid):
    if rid not in AB["present"]:
        return 4, "WHAT-IF ONLY (not in the model)"
    if AB["bounds"][rid] == [0.0, 0.0]:
        return 3, "CLOSED (0,0) as built"
    if maxabs(AB["capacity_taskA"][rid]) == 0:
        return 2, "blocked as built (Task A)"
    if clean(AB["pfba"]["fluxes"][rid]) != 0:
        return 6, "carries pFBA flux as built"
    return 5, "can carry flux as built, 0 at the pFBA optimum"


rxns, status_values = {}, {}
for rid in mapped:
    r = dict(RX[rid])
    status_values[rid], r["sublabel"] = status(rid)
    rxns[rid] = r
mapped_of = {}
for rid in mapped:
    for mid in RX[rid]["mets"]:
        mapped_of.setdefault(mid, set()).add(rid)


def met_class(mid):
    if mid in biomass_species:
        return 3
    if "as_built" not in MT[mid]["reactions"]:
        return 4
    if all(status_values[r] in (2, 3, 4) for r in mapped_of[mid]):
        return 1
    return 0


mets = {mid: dict(name=MT[mid]["name"], comp=MT[mid]["comp"],
                  degree=MT[mid]["degree"].get("as_built", max(MT[mid]["degree"].values())), gap=False, cls="")
        for r in rxns.values() for mid in r["mets"]}
met_data = {mid: met_class(mid) for mid in mets}
legend = [
    "Nothing in the model was changed. Model e1e7c868, Task A medium (cond_taskA_po1f.json, R612 closed), Gurobi 1e-9. Reaction colour = "
    "scenario value: capacity = max |flux| any feasible state can carry (FVA at 0 % of the optimum); pFBA = flux at the optimum. "
    "Each family's what-if is an in-memory variant from the scouts' scripts, applied ALONE (combined effects were not run).",
    "Node colour: dark red, large = species the what-if adds to biomass_C; purple = species that exist only in a what-if; orange = every "
    "mapped reaction of the species is blocked or closed as built (mapped reactions only); grey = otherwise. Hubs repeat; layout is automatic.",
]
title = f"Five demand families of the 29 FN on {F['model_sha256'][:16]}... as built, with each family's recommended what-if"
mj, dj, _, W, H, placed = build_map("iyali26-fn29-families", title, rxns, mets, SECTIONS, maxw=16000, legend=legend, dx=520)


def value(variant, key, rid):
    v = V[variant]
    if rid not in v["present"]:
        return None
    return clean(v["pfba"]["fluxes"][rid]) if key == "pfba" else maxabs(v["capacity_taskA"][rid])


def scen(key, pick):
    out = {}
    for rid in mapped:
        variant = pick(fam_of[rid])
        val = value(variant, key, rid) if variant else None
        if val is not None:
            out[rid] = val
    return out


scenarios = {
    "status": status_values,
    "pfba_as_built": scen("pfba", lambda fam: "as_built"),
    "capacity_as_built": scen("capacity", lambda fam: "as_built"),
    "capacity_whatif": scen("capacity", lambda fam: WHATIF[fam]),
    "pfba_whatif": scen("pfba", lambda fam: WHATIF[fam]),
    "pfba_whatif_nglycan_B": scen("pfba", lambda fam: ALT.get(fam)),
}


def describe(variant, rid):
    m = V[variant]
    if rid not in m["present"]:
        return "not in this variant"
    cap = "[" + ", ".join(f"{clean(v):.3g}" for v in m["capacity_taskA"][rid]) + "]"
    return f"bounds {m['bounds'][rid]}, pFBA {clean(m['pfba']['fluxes'][rid]):.3g}, capacity {cap}"


tips = {}
for rid in mapped:
    fam = fam_of[rid]
    parts = [status(rid)[1], f"GPR: {RX[rid]['gpr'] or 'none'}", f"as built: {describe('as_built', rid)}",
             f"what-if {WHATIF[fam]}: {describe(WHATIF[fam], rid)}"]
    if fam in ALT:
        parts.append(f"what-if {ALT[fam]}: {describe(ALT[fam], rid)}")
    if rid in changed[fam]:
        parts.append("changed by the what-if (bounds or stoichiometry; the map draws the as-built equation)")
    tips[rid] = " · ".join(parts)
NAMES = {0: "has a mapped reaction that can carry flux as built", 1: "every mapped reaction blocked or closed as built",
         3: "added to biomass_C by the what-if", 4: "exists only in the what-if"}
for mid in mets:
    tips[mid] = f"{NAMES[met_data[mid]]} · {mets[mid]['degree']} reactions in the model"
os.makedirs(f"{work}/maps", exist_ok=True)
os.makedirs(f"{work}/data", exist_ok=True)
json.dump(mj, open(f"{work}/maps/fn29_families.json", "w"))
json.dump(dj, open(f"{work}/maps/fn29_families_deco.json", "w"))
for name, data in scenarios.items():
    json.dump(data, open(f"{work}/data/{name}.json", "w"))
nz = lambda s: sorted(rid for rid, v in scenarios[s].items() if v != 0)
counts = {
    "mapped_reactions": len(mapped),
    "per_section": {fam: {"as_built": len(SEC[fam]), "what_if_only": added[fam], "changed_by_what_if": changed[fam],
                          "status_as_built": {lab: sum(1 for r in SEC[fam] if status_values[r] == v)
                                              for v, lab in ((2, "blocked"), (3, "closed"), (5, "can_carry_flux_pfba_0"), (6, "pfba_flux"))},
                          "pfba_nonzero_what_if": sum(1 for r in SEC[fam] + added[fam] if (scenarios["pfba_whatif"].get(r) or 0) != 0)}
                    for fam in ORDER},
    "nonzero_pfba": {s: nz(s) for s in ("pfba_as_built", "pfba_whatif", "pfba_whatif_nglycan_B")},
    "metabolite_nodes": {str(k): sorted(m for m, v in met_data.items() if v == k) if k in (3, 4) else sum(1 for v in met_data.values() if v == k)
                         for k in (0, 1, 3, 4)},
    "wt": {k: v["wt"] for k, v in V.items()},
    "egc_R_NTP1_no_uptake": {k: v["egc_R_NTP1_no_uptake"] for k, v in V.items()},
    "lp_count": F["lp_count"],
}
json.dump({"title": title, "W": W, "H": H, "met_data": met_data, "tips": tips,
           "sections": [(lab, len(ids)) for lab, ids in SECTIONS], "counts": counts},
          open(f"{work}/index.json", "w"), indent=1)
print(f"map {W}x{H}, {len(placed)} placed nodes;", json.dumps(counts["per_section"]))
