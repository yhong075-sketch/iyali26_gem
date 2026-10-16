"""Task D step 2: holistic combinatorial test of the complex I escape module (rules in TASK.md, frozen part table).

Tier 1: every subset of the 12 groups in parts_definition.json, configurations A and V2, WT and complex I KO
(R2062) with the subset closed. KO/WT = KO / WT under the same closures; WT cost = 1 - WT(S) / WT(none).
Hit = WT feasible and (KO/WT < cutoff, or KO with no solution = unresolved, not lethal). Minimal cut set = a hit
with no proper subset that is a hit (KO/WT is not monotone in S, so every subset is checked). Accepted = minimal
cut set with WT cost <= 2 %. Leave-one-out from the full set: per group (from the enumeration) and per part.
Tier 2 (configuration A): for the five cheapest minimal cut sets at 15 %, every multi-part group in the set is
resolved by all 2^m subsets of its parts, the other groups of the set held closed; the combined minimal part sets
are then re-checked in A and V2.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/complexI_transhydrogenation_20261012/combinatorial.py
Writes combinatorial.json and combinatorial_tier1.csv here. Never overwrites.
"""
import csv
import io
import itertools
import json
import time
from datetime import datetime, timezone

from common import CONFIGS, CUTOFF, HERE, Counter_, check_pins, close, growth, ko, load, sha, write_new
from discovery import merged_closures

CUTOFFS = (0.01, 0.05, 0.10, 0.15)
MAX_WT_COST = 0.02
TIER2_SETS = 5
MAX_COMBINED = 64


def evaluate(wt, kom, closures):
    out = {}
    for label, m in (("WT", wt), ("KO", kom)):
        with m:
            close(m, closures)
            out[label] = growth(m)
    w, k = out["WT"][0], out["KO"][0]
    return {"WT": w, "WT_status": out["WT"][1], "KO": k, "KO_status": out["KO"][1],
            "KO_over_WT": None if not w or k is None else k / w}


def is_hit(row, cutoff):
    return row["WT"] is not None and (row["KO"] is None or row["KO_over_WT"] < cutoff)


def minimal_sets(rows, cutoff):
    """rows: {frozenset: row}. Every hit with no proper subset that is a hit."""
    hits = sorted((s for s, r in rows.items() if is_hit(r, cutoff)), key=len)
    minimal = []
    for s in hits:
        if not any(h < s for h in hits if len(h) < len(s)):
            minimal.append(s)
    return minimal


def summarise(s, row, open_wt, order):
    return {"closed": [g for g in order if g in s], "size": len(s), "WT": row["WT"],
            "WT_cost": None if row["WT"] is None else 1 - row["WT"] / open_wt, "KO": row["KO"],
            "KO_status": row["KO_status"], "KO_over_WT": row["KO_over_WT"],
            "unresolved": row["KO"] is None and row["WT"] is not None}


def main():
    inputs = check_pins()
    started = time.monotonic()
    path = HERE / "parts_definition.json"
    definition = json.loads(path.read_text())
    parts = {p["id"]: p for p in definition["parts"]}
    groups = definition["groups"]
    order = list(groups)
    assert len(order) <= 12
    seen = [rid for p in parts.values() for rid in p["closures"]]
    assert len(seen) == len(set(seen)), "a reaction is closed by two parts"
    group_closures = {g: merged_closures([parts[p] for p in groups[g]]) for g in order}
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs,
              "parts_definition_sha256": sha(path), "groups": groups, "group_closures": group_closures,
              "cutoffs": CUTOFFS, "max_wt_cost": MAX_WT_COST, "configs": {}}
    models = {}
    csv_rows = []
    for config in CONFIGS:
        wt, _ = load(config)
        kom, _ = load(config)
        ko(kom)
        models[config] = (wt, kom)
        rows = {}
        for n in range(len(order) + 1):
            for combo in itertools.combinations(order, n):
                closures = {}
                for g in combo:
                    closures.update(group_closures[g])
                rows[frozenset(combo)] = evaluate(wt, kom, closures)
        open_row, full = rows[frozenset()], frozenset(order)
        assert open_row["WT"], f"{config}: WT infeasible with nothing closed"
        open_wt = open_row["WT"]
        analysis = {"open": open_row, "full_set": summarise(full, rows[full], open_wt, order), "minimal": {}}
        for cutoff in CUTOFFS:
            sets = [summarise(s, rows[s], open_wt, order) for s in minimal_sets(rows, cutoff)]
            sets.sort(key=lambda x: (x["WT_cost"], x["size"]))
            analysis["minimal"][str(cutoff)] = {
                "count": len(sets), "accepted_count": sum(1 for x in sets if x["WT_cost"] <= MAX_WT_COST),
                "hits_total": sum(1 for r in rows.values() if is_hit(r, cutoff)), "sets": sets}
        analysis["leave_one_out_groups"] = {g: summarise(full - {g}, rows[full - {g}], open_wt, order) for g in order}
        all_closures = {}
        for g in order:
            all_closures.update(group_closures[g])
        analysis["leave_one_out_parts"] = {}
        for pid, part in parts.items():
            reopened = {rid: how for rid, how in all_closures.items() if rid not in part["closures"]}
            row = evaluate(wt, kom, reopened)
            analysis["leave_one_out_parts"][pid] = {**row, "WT_cost": None if row["WT"] is None else 1 - row["WT"] / open_wt}
        record["configs"][config] = analysis
        for s, r in rows.items():
            csv_rows.append([config, "+".join(g for g in order if g in s) or "none", len(s), r["WT"], r["WT_status"],
                             None if r["WT"] is None else 1 - r["WT"] / open_wt, r["KO"], r["KO_status"], r["KO_over_WT"]])
        m15 = analysis["minimal"]["0.15"]
        print(f"{config}: open KO/WT {open_row['KO_over_WT']:.4f}; full set KO/WT {analysis['full_set']['KO_over_WT']} "
              f"WT cost {analysis['full_set']['WT_cost']:.4f}; minimal@15% {m15['count']} accepted {m15['accepted_count']}; "
              f"LPs {Counter_.lps}", flush=True)
    # Tier 2 (configuration A), then re-check the combined part sets in A and V2.
    tier2 = []
    wt, kom = models["A"]
    open_wt = record["configs"]["A"]["open"]["WT"]
    for cut in record["configs"]["A"]["minimal"]["0.15"]["sets"][:TIER2_SETS]:
        entry = {"cut_set": cut["closed"], "cut_WT_cost": cut["WT_cost"], "groups": {}}
        for g in cut["closed"]:
            members = groups[g]
            if len(members) < 2:
                entry["groups"][g] = {"members": members, "minimal_part_sets": [members]}
                continue
            others = {}
            for h in cut["closed"]:
                if h != g:
                    others.update(group_closures[h])
            rows = {}
            for n in range(len(members) + 1):
                for combo in itertools.combinations(members, n):
                    closures = dict(others)
                    for p in combo:
                        closures.update(parts[p]["closures"])
                    rows[frozenset(combo)] = evaluate(wt, kom, closures)
            minimal = minimal_sets(rows, CUTOFF)
            entry["groups"][g] = {
                "members": members,
                "minimal_part_sets": [[p for p in members if p in s] for s in minimal],
                "rows": [{"closed": [p for p in members if p in s], **r,
                          "WT_cost": None if r["WT"] is None else 1 - r["WT"] / open_wt} for s, r in rows.items()]}
        choices = [entry["groups"][g]["minimal_part_sets"] for g in cut["closed"]]
        combined = []
        for pick in itertools.islice(itertools.product(*choices), MAX_COMBINED):
            part_set = sorted({p for chosen in pick for p in chosen})
            closures = merged_closures([parts[p] for p in part_set])
            checks = {}
            for config in CONFIGS:
                w_, k_ = models[config]
                row = evaluate(w_, k_, closures)
                o = record["configs"][config]["open"]["WT"]
                checks[config] = {**row, "WT_cost": None if row["WT"] is None else 1 - row["WT"] / o,
                                  "hit": is_hit(row, CUTOFF)}
            combined.append({"parts": part_set, "closures": closures, "checks": checks})
        entry["combined_part_sets"] = combined
        entry["combined_truncated"] = len(list(itertools.product(*choices))) > MAX_COMBINED
        tier2.append(entry)
        print("tier 2:", cut["closed"], "->", [c["parts"] for c in combined][:4], flush=True)
    record["tier2"] = tier2
    record["lps"] = Counter_.lps
    record["elapsed_seconds"] = time.monotonic() - started
    write_new(HERE / "combinatorial.json", json.dumps(record, indent=1, default=str) + "\n")
    buffer = io.StringIO()
    w = csv.writer(buffer, lineterminator="\n")
    w.writerow(["config", "closed_groups", "size", "WT", "WT_status", "WT_cost", "KO", "KO_status", "KO_over_WT"])
    w.writerows(csv_rows)
    write_new(HERE / "combinatorial_tier1.csv", buffer.getvalue())
    print("LPs", Counter_.lps, "seconds", round(time.monotonic() - started))


if __name__ == "__main__":
    main()
