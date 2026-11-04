"""Score recount for the candidate: baseline sweep rows for the 1,067 genes outside the V-ATPase rule (1,050 in scope), Verify-1 rows
for the 15 rule genes, mtDNA + untargeted genes out of scope (new list). Not an LP; the other 1,067 classes are taken
from the 2bc263f2 sweep (rule B), so their being unchanged is inferred, not measured.

Usage (from the repo root): python3 model/reports/vatpase_gpr_a1_20261103/scripts/recount.py
"""
import csv
import json
from pathlib import Path

D = Path("model/reports/vatpase_gpr_a1_20261103")
SWEEP = Path("model/reports/ceramide_pool_biomass_20261101/ko_screen_ceramide_pool.tsv")
CUTOFFS = ("0.01", "0.05", "0.1", "0.15")
verify = json.loads((D / "verify1_ko_targets.json").read_text())
mtdna = set(json.loads(Path("model/expected/mtdna_genes_out_of_scope.json").read_text())["genes"])
untargeted = set(json.loads(Path("model/expected/untargeted_genes_out_of_scope.json").read_text())["genes"])
rows = {r["gene"]: r for r in csv.DictReader(SWEEP.open(), delimiter="\t")}


def classes(scope):
    out = {}
    for gene, r in rows.items():
        if gene in scope:
            continue
        if gene in verify["knockouts"]:
            k = verify["knockouts"][gene]
            out[gene] = (k["experimental_positive"], {c: k["class"][str(float(c))] for c in CUTOFFS})
        else:
            out[gene] = (r["experimental_positive"] == "True", {c: r[f"class_B_{c}"] for c in CUTOFFS})
    k = verify["knockouts"]["YALI1F25663g"]
    out["YALI1F25663g"] = (k["experimental_positive"], {c: k["class"][str(float(c))] for c in CUTOFFS})
    return out


def score(table):
    result = {}
    for c in CUTOFFS:
        n = {"TP": 0, "FN": 0, "FP": 0, "TN": 0, "unresolved": 0}
        for positive, cls in table.values():
            call = cls[c]
            n["unresolved" if call == "unresolved" else ("TP" if positive else "FP") if call == "essential"
              else ("FN" if positive else "TN")] += 1
        result[c] = n
    return result


baseline = {c: {k: 0 for k in ("TP", "FN", "FP", "TN", "unresolved")} for c in CUTOFFS}
old_untargeted = {"YALI0F04114g"}
for gene, r in rows.items():
    if gene in mtdna | old_untargeted:
        continue
    for c in CUTOFFS:
        call, positive = r[f"class_B_{c}"], r["experimental_positive"] == "True"
        baseline[c]["unresolved" if call == "unresolved" else ("TP" if positive else "FP") if call == "essential"
                    else ("FN" if positive else "TN")] += 1
new_table = classes(mtdna | untargeted)
out = {"baseline_2bc263f2_old_scope": baseline, "candidate_ab320e18_new_scope": score(new_table),
       "in_model_positives_in_scope": sum(p for p, _ in new_table.values()),
       "genes_scored": len(new_table), "untargeted": sorted(untargeted),
       "verify1_rows_differing_from_baseline_class": sorted(
           g for g in verify["knockouts"] if g in rows
           and any(verify["knockouts"][g]["class"][str(float(c))] != rows[g][f"class_B_{c}"] for c in CUTOFFS)),
       "method": __doc__.split("\n\n")[0]}
(D / "score_recount.json").write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps({k: out[k] for k in list(out)[:6]}, indent=1))
