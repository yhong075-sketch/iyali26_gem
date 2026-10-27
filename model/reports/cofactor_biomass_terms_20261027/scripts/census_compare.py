"""Compare two gapfind.py runs (Ponce-de-Leon census): blocked reactions, gap metabolites, modules (read-only).

Usage: python -I census_compare.py BASE_TAG_COND_reactions.csv NEW_TAG_COND_reactions.csv \
           BASE_TAG_COND_metabolites.csv NEW_TAG_COND_metabolites.csv BASE_summary.json NEW_summary.json OUT.json
Reports the reactions that became active and the reactions that became blocked (with name, GPR, module of the
base run and primary cause from the census blocked_reactions.csv when given as an 8th argument), the gap
metabolites that changed, and the summary counts of both runs.
"""
import csv
import json
import sys

base_r, new_r, base_m, new_m, base_s, new_s, out = sys.argv[1:8]
causes = {}
if len(sys.argv) > 8:
    causes = {row["id"]: row for row in csv.DictReader(open(sys.argv[8]))}


def blocked(path):
    return {row["id"]: row for row in csv.DictReader(open(path)) if row["blocked"] == "True"}


def gaps(path):
    return {row["id"]: row for row in csv.DictReader(open(path)) if row["gap"] == "True"}


b_blocked, n_blocked = blocked(base_r), blocked(new_r)
b_gap, n_gap = gaps(base_m), gaps(new_m)
freed = sorted(set(b_blocked) - set(n_blocked))
newly_blocked = sorted(set(n_blocked) - set(b_blocked))
result = {
    "base_summary": json.load(open(base_s)), "new_summary": json.load(open(new_s)),
    "blocked": {"base": len(b_blocked), "new": len(n_blocked), "freed": len(freed), "newly_blocked": len(newly_blocked)},
    "gap_metabolites": {"base": len(b_gap), "new": len(n_gap), "no_longer_gap": sorted(set(b_gap) - set(n_gap)),
                        "newly_gap": sorted(set(n_gap) - set(b_gap))},
    "freed_reactions": [{"id": rid, "name": b_blocked[rid]["name"][:60], "gpr": b_blocked[rid]["gpr"], "um_base": b_blocked[rid]["um"],
                         "cause_base": causes.get(rid, {}).get("cause", "")} for rid in freed],
    "newly_blocked_reactions": [{"id": rid, "name": n_blocked[rid]["name"][:60], "gpr": n_blocked[rid]["gpr"], "um_new": n_blocked[rid]["um"]}
                                for rid in newly_blocked],
}
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
print("blocked", result["blocked"], "| gaps", {k: (v if isinstance(v, int) else len(v)) for k, v in result["gap_metabolites"].items()})
print("freed by module (base):", sorted({(r["um_base"], sum(1 for x in result["freed_reactions"] if x["um_base"] == r["um_base"])) for r in result["freed_reactions"]}, key=lambda t: -t[1]))
print("newly blocked:", [r["id"] for r in result["newly_blocked_reactions"]])
print("UM sizes base", result["base_summary"]["UM_sizes"][:8], "new", result["new_summary"]["UM_sizes"][:8], "| UM count", result["base_summary"]["UM"], "->", result["new_summary"]["UM"])
