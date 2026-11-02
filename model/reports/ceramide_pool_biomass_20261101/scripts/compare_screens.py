"""Class changes between the c1f70ac4 screen and the ceramide-pool screen (read-only).
Usage: python compare_screens.py BASE.tsv BASE_SUMMARY.json NEW.tsv NEW_SUMMARY.json OUT.json
"""
import csv
import json
import sys

base_tsv, base_sum, new_tsv, new_sum, out = sys.argv[1:6]
read = lambda p: {r["gene"]: r for r in csv.DictReader(open(p), delimiter="\t")}
base, new = read(base_tsv), read(new_tsv)
cols = [c for c in next(iter(new.values())) if c.startswith("class_")]
changes = []
for gene in sorted(set(base) | set(new)):
    b, n = base.get(gene), new.get(gene)
    if b is None or n is None:
        changes.append({"gene": gene, "missing_in": "base" if b is None else "new"})
        continue
    diff = {c: [b[c], n[c]] for c in cols if b[c] != n[c]}
    if diff or b["status"] != n["status"]:
        changes.append({"gene": gene, "experimental_positive": n["experimental_positive"], "ko_before": b["ko"],
                        "ko_after": n["ko"], "status": [b["status"], n["status"]], "classes": diff})
sb, sn = json.load(open(base_sum)), json.load(open(new_sum))
scores = {scope: {cut: {"before": {k: sb["scores"][scope][cut][k] for k in ("TP", "FN", "FP", "TN", "unresolved", "mcc")},
                        "after": {k: sn["scores"][scope][cut][k] for k in ("TP", "FN", "FP", "TN", "unresolved", "mcc")}}
                  for cut in sn["scores"][scope]} for scope in sn["scores"]}
json.dump({"changes": changes, "scores": scores, "lp_count": sn["lp_count"], "wt": [sb["wt"]["growth"], sn["wt"]["growth"]]},
          open(out, "w"), indent=1)
print(len(changes), "genes changed")
for c in changes:
    print(c["gene"], c.get("experimental_positive"), c.get("ko_before", "")[:10], "->", c.get("ko_after", "")[:10], c.get("status"),
          sorted(set(v[1] for v in c.get("classes", {}).values())), len(c.get("classes", {})), "class columns")
for cut, v in scores["B_mtdna_untargeted"].items():
    print("B_mtdna_untargeted", cut, v["before"], "->", v["after"])
print("lp_count", sn["lp_count"])
