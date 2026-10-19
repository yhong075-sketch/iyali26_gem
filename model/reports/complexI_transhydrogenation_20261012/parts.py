"""Task D: round-0 part pool for the complex I transhydrogenation module (rules in TASK.md, fixed before running).

Pool: (1) NAD/NADP pairs (same compartment, identical stoichiometry apart from NAD(H) vs NADP(H), H+ ignored,
same or opposite orientation) on the base XML and again after the overlay; (2) reactions carrying NAD(H) or
NADP(H) in two compartments; (3) reactions holding both an NAD and an NADP species; (4) earlier parts;
(5) handoff parts. In configuration A, for WT and the complex I KO: FVA at any growth >= 0 (pruning) and at the
optimum (escape direction), then the single-closure WT cost and KO/WT of every escape direction.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/complexI_transhydrogenation_20261012/parts.py
Writes parts_pool.json and parts_pool.csv here. Never overwrites.
"""
import csv
import io
import json
import time
from collections import defaultdict
from datetime import datetime, timezone

import cobra

from common import HERE, PINS, TOL, Counter_, check_pins, close, fva, growth, ko, load, per100, write_new

EARLIER = ("R361", "R1931", "R1936", "R328", "R691", "R1310")
HANDOFF = {"HMG-CoA route": ("R412", "R1172"), "R1708": ("R1708",)}
COFACTOR = {"NAD": ("ox", ""), "NADH": ("red", ""), "NADP(+)": ("ox", "P"), "NADPH": ("red", "P")}


def cofactors(model):
    kind = {}
    for met in model.metabolites:
        label = (met.name or "").split("_")[0]
        if label in COFACTOR:
            kind[met.id] = COFACTOR[label] + (met.compartment,)
    protons = {met.id for met in model.metabolites if met.formula == "H"}
    return kind, protons


def pair_scan(model):
    kind, protons = cofactors(model)

    def signature(r, sign=1):
        d = defaultdict(float)
        for met, c in r.metabolites.items():
            if met.id in protons:
                continue
            key = f"NAD[{kind[met.id][0]}]@{kind[met.id][2]}" if met.id in kind else met.id
            d[key] += sign * c
        return frozenset((k, round(v, 6)) for k, v in d.items() if abs(v) > 1e-9)

    def flavour(r):
        return {kind[m.id][1] for m in r.metabolites if m.id in kind}

    nad = [r for r in model.reactions if flavour(r) == {""}]
    index = defaultdict(list)
    for r in model.reactions:
        if flavour(r) == {"P"}:
            index[signature(r)].append(r.id)
    pairs = []
    for r in nad:
        for q in index.get(signature(r), []):
            pairs.append((r.id, q, "same"))
        for q in index.get(signature(r, -1), []):
            pairs.append((r.id, q, "opposite"))
    transport = sorted(r.id for r in model.reactions if len({kind[m.id][2] for m in r.metabolites if m.id in kind}) > 1)
    mixed = sorted(r.id for r in model.reactions if flavour(r) == {"", "P"})
    return pairs, transport, mixed


def directions(lo, hi, wlo, whi):
    """Escape directions: usable at the KO optimum and not required at the WT optimum."""
    out = []
    if hi > TOL and not wlo > TOL:
        out.append("fwd")
    if lo < -TOL and not whi < -TOL:
        out.append("rev")
    return out


def main():
    inputs = check_pins()
    started = time.monotonic()
    base = cobra.io.read_sbml_model(str(PINS["model"][0]))
    base_pairs, base_transport, base_mixed = pair_scan(base)
    wt, applied = load("A")
    over_pairs, over_transport, over_mixed = pair_scan(wt)
    sources = defaultdict(set)
    for a, b, orient in base_pairs:
        sources[a].add(f"pair:{a}/{b}:{orient}:base")
        sources[b].add(f"pair:{a}/{b}:{orient}:base")
    for a, b, orient in over_pairs:
        sources[a].add(f"pair:{a}/{b}:{orient}:overlay")
        sources[b].add(f"pair:{a}/{b}:{orient}:overlay")
    for rid in set(base_transport) | set(over_transport):
        sources[rid].add("transport")
    for rid in set(base_mixed) | set(over_mixed):
        sources[rid].add("NAD+NADP species")
    for rid in EARLIER:
        sources[rid].add("earlier part")
    for name, rids in HANDOFF.items():
        for rid in rids:
            sources[rid].add(f"handoff:{name}")
    pool = sorted(sources, key=lambda x: (len(x), x))
    kom, _ = load("A")
    ko(kom)
    w_opt, k_opt = growth(wt)[0], growth(kom)[0]
    ranges = {"WT_any": fva(wt, pool, 0.0), "KO_any": fva(kom, pool, 0.0),
              "WT_opt": fva(wt, pool, 1.0), "KO_opt": fva(kom, pool, 1.0)}
    rows = []
    for rid in pool:
        r = wt.reactions.get_by_id(rid)
        blocked = all(abs(x) <= TOL for key in ("WT_any", "KO_any") for x in ranges[key][rid])
        overlay_closed = tuple(r.bounds) == (0.0, 0.0) and tuple(base.reactions.get_by_id(rid).bounds) != (0.0, 0.0)
        esc = directions(*ranges["KO_opt"][rid], *ranges["WT_opt"][rid])
        singles = {}
        for how in esc:
            res = {}
            for label, m in (("WT", wt), ("KO", kom)):
                with m:
                    close(m, {rid: how})
                    res[label] = growth(m)
            singles[how] = {"WT": res["WT"][0], "WT_cost": None if res["WT"][0] is None else 1 - res["WT"][0] / w_opt,
                            "KO": res["KO"][0], "KO_status": res["KO"][1],
                            "KO_over_WT": None if not res["WT"][0] or res["KO"][0] is None else res["KO"][0] / res["WT"][0]}
        rows.append({"reaction": rid, "name": r.name, "sources": sorted(sources[rid]), "bounds_base":
                     list(base.reactions.get_by_id(rid).bounds), "bounds_A": list(r.bounds), "gpr": r.gene_reaction_rule,
                     "equation": r.build_reaction_string(use_metabolite_names=True),
                     **{f"{k}_per100": [round(per100(x), 4) for x in ranges[k][rid]] for k in ranges},
                     "blocked": blocked, "overlay_closed": overlay_closed,
                     "escape_directions": esc, "single_closure": singles})
    dropped = {row["reaction"] for row in rows if row["blocked"] or row["overlay_closed"]}
    pair_status = []
    for a, b, orient in sorted(set((a, b, o) for a, b, o in base_pairs + over_pairs)):
        pair_status.append({"nad": a, "nadp": b, "orientation": orient, "in_base": (a, b, orient) in base_pairs,
                            "after_overlay": (a, b, orient) in over_pairs, "dropped": a in dropped or b in dropped})
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs, "config": "A",
              "WT_growth": w_opt, "KO_growth": k_opt, "pairs": pair_status,
              "transport": {"base": base_transport, "overlay": over_transport},
              "nad_nadp_species": {"base": base_mixed, "overlay": over_mixed}, "pool": rows,
              "lps": Counter_.lps, "elapsed_seconds": time.monotonic() - started}
    write_new(HERE / "parts_pool.json", json.dumps(record, indent=1, default=str) + "\n")
    buffer = io.StringIO()
    w = csv.writer(buffer, lineterminator="\n")
    w.writerow(["reaction", "name", "sources", "bounds_A", "WT_any", "KO_any", "WT_opt", "KO_opt", "blocked",
                "overlay_closed", "escape_directions", "closure_WT_cost", "closure_KO_over_WT", "gpr"])
    for row in rows:
        w.writerow([row["reaction"], row["name"], "; ".join(row["sources"]), row["bounds_A"], row["WT_any_per100"],
                    row["KO_any_per100"], row["WT_opt_per100"], row["KO_opt_per100"], row["blocked"],
                    row["overlay_closed"], "+".join(row["escape_directions"]),
                    "; ".join(f"{h}:{s['WT_cost']:.4f}" for h, s in row["single_closure"].items() if s["WT_cost"] is not None),
                    "; ".join(f"{h}:{s['KO_over_WT']:.4f}" for h, s in row["single_closure"].items() if s["KO_over_WT"] is not None),
                    row["gpr"]])
    write_new(HERE / "parts_pool.csv", buffer.getvalue())
    print("pool", len(pool), "dropped", len(dropped), "LPs", Counter_.lps)


if __name__ == "__main__":
    main()
