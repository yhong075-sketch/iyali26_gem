"""Collect the three official predictor outputs into one table and apply the fixed decision rule.

Run from the repository root (no network, no solves):
    .venv/bin/python -B model/reports/r132_mito_adh_20261012/collect_predictions.py
"""
import csv
import hashlib
import json
import os
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ART = Path(os.environ["IYALI26_RESEARCH_ROOT"]) / "artifacts/tasks" / HERE.name
TP = ART / "predictors/targetp2_6ACE82AE0001C6E02ECCAE2B/output_protein_type.txt"
MF = ART / "predictors/mitofates_20261014.040043.60854"
DL = ART / "predictors/deeploc2.0_6ACE830B0001C7579B8DCFFC/results_6ACE830B0001C7579B8DCFFC.csv"
DL_MITO_THRESHOLD = 0.6220  # read from the DeepLoc 2.0 results page, saved as thresholds_from_results_page.csv
TP_MTP, MF_CUT = 0.5, 0.385  # fixed in TASK.md before any run


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    rows = {}
    for line in TP.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        rows[f[0]] = {"protein_id": f[0], "gene": f[0].split("_")[0], "refseq": f[0].split("_", 1)[1],
                      "targetp_prediction": f[1], "targetp_other": float(f[2]), "targetp_sp": float(f[3]),
                      "targetp_mtp": float(f[4]), "targetp_cs": f[5] if len(f) > 5 else ""}
    for pid in rows:
        text = (MF / f"resultText_{pid}.txt").read_text()
        field = lambda key: re.search(rf"{key}\s*:\s*(.+)", text).group(1).strip()
        rows[pid].update(mitofates_probability=float(field("Probability of presequence")),
                         mitofates_prediction=field("Prediction"), mitofates_cleavage=field("Cleavage site"))
    with DL.open() as f:
        for r in csv.DictReader(f):
            rows[r["Protein_ID"]].update(deeploc_localizations=r["Localizations"], deeploc_signals=r["Signals"],
                                         deeploc_mitochondrion=round(float(r["Mitochondrion"]), 4),
                                         deeploc_mito_threshold=DL_MITO_THRESHOLD)
    for r in rows.values():
        r["rule_mitochondrial"] = r["targetp_mtp"] >= TP_MTP or r["mitofates_probability"] >= MF_CUT
        r["predictors_crossing"] = ";".join(n for n, hit in (("TargetP", r["targetp_mtp"] >= TP_MTP),
                                                              ("MitoFates", r["mitofates_probability"] >= MF_CUT),
                                                              ("DeepLoc(not in rule)", r["deeploc_mitochondrion"] >= DL_MITO_THRESHOLD)) if hit)
    out = HERE / "localization_predictions.tsv"
    fields = list(next(iter(rows.values())))
    with out.open("x") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows.values())
    hits = [r["gene"] for r in rows.values() if r["rule_mitochondrial"]]
    record = {"rule": f"mitochondrial if TargetP mTP >= {TP_MTP} or MitoFates >= {MF_CUT}", "genes_mitochondrial": hits,
              "branch": "2 (restrict GPR)" if hits else "1 (remove R132)",
              "inputs": {str(p.relative_to(ART)): sha(p) for p in [TP, DL, *sorted(MF.glob("resultText_*.txt"))]},
              "output": {out.name: sha(out)}}
    (HERE / "decision_rule_result.json").open("x").write(json.dumps(record, indent=1) + "\n")
    print(json.dumps(record, indent=1))


if __name__ == "__main__":
    main()
