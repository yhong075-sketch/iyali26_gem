"""Task D: one escape-discovery round (rules in TASK.md, fixed before running). Configuration A only.

Round r closes every part in parts_definition.json with round_added < r, in WT and in the complex I KO, then
reports: growth and KO/WT; the four node balances (pFBA); and, for every reaction carrying >= 5 per 100
glucose in either pFBA vector, every watched reaction and every node-balance term >= 1, its FVA range at the
WT and KO optima and its escape directions (usable at the KO optimum, not required at the WT optimum).
"WT-free step" = an escape direction carried at >= 5 per 100 glucose in the KO pFBA vector.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/complexI_transhydrogenation_20261012/discovery.py <r>
Writes discovery_round<r>.json and .csv here. Never overwrites.
"""
import csv
import io
import json
import sys
import time
from datetime import datetime, timezone

from common import HERE, NODES, Counter_, check_pins, close, fva, growth, ko, load, per100, pfba_vector, sha, write_new
from parts import directions

STEP = 5.0  # per 100 glucose


def merged_closures(parts):
    out = {}
    for part in parts:
        for rid, how in part["closures"].items():
            prior = out.get(rid)
            out[rid] = how if prior in (None, how) else "whole"
    return out


def node_balance(model, vector, uptake):
    out = {}
    for node, mid in NODES.items():
        met = model.metabolites.get_by_id(mid)
        terms = [(r.id, per100(r.get_coefficient(met) * vector[r.id], uptake)) for r in met.reactions]
        terms = sorted(((rid, round(c, 3)) for rid, c in terms if abs(c) > 1e-6), key=lambda t: -abs(t[1]))
        out[node] = {"produced": round(sum(c for _, c in terms if c > 0), 3), "terms": terms}
    return out


def main(round_):
    inputs = check_pins()
    started = time.monotonic()
    definition_path = HERE / "parts_definition.json"
    definition = json.loads(definition_path.read_text())
    applied = [p for p in definition["parts"] if p["round_added"] < round_]
    closures = merged_closures(applied)
    wt, _ = load("A")
    kom, _ = load("A")
    ko(kom)
    for m in (wt, kom):
        close(m, closures)
    w, w_status = growth(wt)
    k, k_status = growth(kom)
    assert w, f"WT infeasible with the round-{round_} closures ({w_status})"
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs, "round": round_,
              "parts_definition_sha256": sha(definition_path), "parts_applied": [p["id"] for p in applied],
              "closures": closures, "WT": w, "WT_status": w_status, "KO": k, "KO_status": k_status,
              "KO_over_WT": None if k is None else k / w}
    print(f"round {round_}: closures {len(closures)} WT {w:.6f} KO {k} KO/WT {record['KO_over_WT']}", flush=True)
    if k is None:
        record.update(lps=Counter_.lps, elapsed_seconds=time.monotonic() - started)
        write_new(HERE / f"discovery_round{round_}.json", json.dumps(record, indent=1, default=str) + "\n")
        return
    vk, vw = pfba_vector(kom), pfba_vector(wt)
    uk, uw = -vk["R1070"], -vw["R1070"]
    record["nodes_KO"] = node_balance(kom, vk, uk)
    record["nodes_WT"] = node_balance(wt, vw, uw)
    candidates = {rid for rid, v in vk.items() if abs(per100(v, uk)) >= STEP}
    candidates |= {rid for rid, v in vw.items() if abs(per100(v, uw)) >= STEP}
    candidates |= set(definition["watch"])
    candidates |= {rid for nodes in (record["nodes_KO"], record["nodes_WT"]) for n in nodes.values()
                   for rid, c in n["terms"] if abs(c) >= 1}
    candidates = sorted(candidates)
    rk, rw = fva(kom, candidates, 1.0), fva(wt, candidates, 1.0)
    rows = []
    for rid in candidates:
        r = kom.reactions.get_by_id(rid)
        esc = directions(*rk[rid], *rw[rid])
        v = per100(vk[rid], uk)
        wt_free = [d for d in esc if (d == "fwd" and v >= STEP) or (d == "rev" and v <= -STEP)]
        rows.append({"reaction": rid, "name": r.name, "KO_pfba": round(v, 3), "WT_pfba": round(per100(vw[rid], uw), 3),
                     "KO_opt": [round(per100(x), 3) for x in rk[rid]], "WT_opt": [round(per100(x), 3) for x in rw[rid]],
                     "escape_directions": esc, "wt_free_step": wt_free, "closed_in_round": rid in closures,
                     "watched": rid in definition["watch"], "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
                     "equation": r.build_reaction_string(use_metabolite_names=True)})
    rows.sort(key=lambda x: (not x["wt_free_step"], -abs(x["KO_pfba"])))
    record["candidates"] = rows
    record["wt_free_steps"] = [x["reaction"] for x in rows if x["wt_free_step"]]
    record.update(lps=Counter_.lps, elapsed_seconds=time.monotonic() - started)
    write_new(HERE / f"discovery_round{round_}.json", json.dumps(record, indent=1, default=str) + "\n")
    buffer = io.StringIO()
    w_ = csv.writer(buffer, lineterminator="\n")
    w_.writerow(["reaction", "name", "KO_pfba_per100", "WT_pfba_per100", "KO_opt_per100", "WT_opt_per100",
                 "escape_directions", "wt_free_step", "watched", "gpr", "equation"])
    for x in rows:
        w_.writerow([x["reaction"], x["name"], x["KO_pfba"], x["WT_pfba"], x["KO_opt"], x["WT_opt"],
                     "+".join(x["escape_directions"]), "+".join(x["wt_free_step"]), x["watched"], x["gpr"], x["equation"]])
    write_new(HERE / f"discovery_round{round_}.csv", buffer.getvalue())
    print("WT-free steps:", record["wt_free_steps"], "LPs", Counter_.lps)


if __name__ == "__main__":
    main(int(sys.argv[1]))
