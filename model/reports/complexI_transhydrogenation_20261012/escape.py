"""Task D step 1: reproduce the complex I escape with node balances (diagnostic only).

Configurations A (Task A medium) and V2 (W29-type defined), both with the holistic overlay and the 13C caps.
For WT and the complex I KO (R2062): growth, the pFBA vector, and FVA at the optimum for every reaction that
touches cytosolic NADH, cytosolic NADPH, matrix NADH or matrix NADPH, plus the central reactions in KEY.
Fluxes are per 100 glucose (pFBA: actual uptake; FVA: the 2.414 uptake bound, checked by R1070's FVA range).

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/complexI_transhydrogenation_20261012/escape.py
Writes escape.json and escape_node_balances.csv here, pFBA vectors to the research workspace. Never overwrites.
"""
import csv
import io
import json
import time
from datetime import datetime, timezone

from common import ART, CONFIGS, HERE, NODES, Counter_, check_pins, fva, growth, ko, load, per100, pfba_vector, write_new

KEY = ("R1070", "R325", "R639", "R549", "R550", "R329", "R330", "R570", "R573", "R2062", "R740", "R1977", "R305",
       "R304", "R171", "R344", "R641", "R534", "R693", "R236", "R487", "R1385", "R614", "R361", "R691", "R1378",
       "R134", "R136", "R997", "R1708", "R412", "R1172", "R1311", "R842", "R690", "R1894", "R825", "R511", "R929",
       "R315", "R2269", "R90", "R536", "R490")


def main():
    inputs = check_pins()
    started = time.monotonic()
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs, "results": {}}
    rows, vectors = [], {}
    for config in CONFIGS:
        for state in ("WT", "KO"):
            model, applied = load(config)
            if state == "KO":
                ko(model)
            value, status = growth(model)
            assert value, f"{config} {state}: no growth ({status})"
            vector = pfba_vector(model)
            uptake = -vector["R1070"]
            node_rxns = sorted({r.id for mid in NODES.values() for r in model.metabolites.get_by_id(mid).reactions})
            ranges = fva(model, sorted(set(node_rxns) | set(KEY)), 1.0)
            vectors[f"{config}|{state}"] = vector
            nodes = {}
            for node, mid in NODES.items():
                met = model.metabolites.get_by_id(mid)
                parts = []
                for r in met.reactions:
                    coef = r.get_coefficient(met)
                    lo, hi = ranges[r.id]
                    contribution = per100(coef * vector[r.id], uptake)
                    rows.append([config, state, node, r.id, r.name, coef, round(per100(vector[r.id], uptake), 6),
                                 round(contribution, 6), round(per100(lo), 6), round(per100(hi), 6),
                                 r.gene_reaction_rule])
                    if abs(contribution) > 1e-6:
                        parts.append((r.id, round(contribution, 3)))
                parts.sort(key=lambda t: -abs(t[1]))
                nodes[node] = {"produced": round(sum(c for _, c in parts if c > 0), 3), "top": parts[:12]}
            record["results"][f"{config}|{state}"] = {
                "growth": value, "status": status, "uptake_pfba": uptake, "uptake_fva": ranges["R1070"],
                "key_pfba_per100": {rid: round(per100(vector[rid], uptake), 3) for rid in KEY},
                "key_fva_per100": {rid: [round(per100(x), 3) for x in ranges[rid]] for rid in KEY},
                "nodes": nodes, "applied": applied if state == "WT" else None}
            print(config, state, round(value, 6), {k: record["results"][f"{config}|{state}"]["key_pfba_per100"][k]
                                                  for k in ("R325", "R570", "R2062", "R693", "R236")}, flush=True)
    for config in CONFIGS:
        wt, k = record["results"][f"{config}|WT"]["growth"], record["results"][f"{config}|KO"]["growth"]
        record["results"][f"{config}|ratio"] = k / wt
        print(config, "KO/WT", round(k / wt, 6))
    record["lps"] = Counter_.lps
    record["elapsed_seconds"] = time.monotonic() - started
    ART.mkdir(parents=True, exist_ok=True)
    write_new(ART / "escape_pfba_vectors.json", json.dumps(vectors) + "\n")
    write_new(HERE / "escape.json", json.dumps(record, indent=1, default=str) + "\n")
    buffer = io.StringIO()
    w = csv.writer(buffer, lineterminator="\n")
    w.writerow(["config", "state", "node", "reaction", "name", "coefficient", "flux_pfba_per100",
                "node_contribution_per100", "fva_min_per100", "fva_max_per100", "gpr"])
    w.writerows(rows)
    write_new(HERE / "escape_node_balances.csv", buffer.getvalue())


if __name__ == "__main__":
    main()
