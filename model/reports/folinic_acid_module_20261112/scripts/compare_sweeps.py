"""Compare the control and candidate sweep TSVs gene by gene (status, KO growth, every class column)."""
import csv, json, sys
from pathlib import Path

def load(p):
    with open(p) as f:
        return {r["gene"]: r for r in csv.DictReader(f, delimiter="\t")}

a, b = load(sys.argv[1]), load(sys.argv[2])
out = {"genes_control": len(a), "genes_candidate": len(b), "gene_sets_equal": set(a) == set(b),
       "status_changes": [], "class_changes": [], "max_abs_ko_growth_diff": 0.0, "largest": []}
cls = [k for k in next(iter(a.values())) if k.startswith("class_")]
diffs = []
for g in sorted(set(a) & set(b)):
    x, y = a[g], b[g]
    if (x["status"], x["passed"], x["cert_certified"]) != (y["status"], y["passed"], y["cert_certified"]):
        out["status_changes"].append({"gene": g, "control": [x["status"], x["passed"], x["cert_certified"]],
                                      "candidate": [y["status"], y["passed"], y["cert_certified"]]})
    for c in cls:
        if x[c] != y[c]:
            out["class_changes"].append({"gene": g, "column": c, "control": x[c], "candidate": y[c]})
    if x["ko"] and y["ko"]:
        diffs.append((abs(float(x["ko"]) - float(y["ko"])), g, float(x["ko"]), float(y["ko"])))
diffs.sort(reverse=True)
out["max_abs_ko_growth_diff"] = diffs[0][0] if diffs else None
out["largest"] = [{"gene": g, "diff": d, "control": c, "candidate": k} for d, g, c, k in diffs[:10]]
Path(sys.argv[3]).write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps({k: v for k, v in out.items() if k != "largest"}, indent=1)); print(out["largest"][:5])
