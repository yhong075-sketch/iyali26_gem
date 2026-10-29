"""Verify the saved TargetP 2.0 and DeepLoc 2.1 outputs against the pinned sequences and tabulate them (read-only).

Usage (from the repo root): python3 model/reports/r673_prolyl_localization_20261029/scripts/verify_predictions.py
Checks: both jobs report 0 failed sequences and exactly the two submitted names; DeepLoc's attention CSV spells each
full submitted sequence; TargetP's per-residue file spells each sequence's N-terminal prefix (it covers 200 residues);
the TargetP zip CRCs pass and its per-residue files equal the separately downloaded ones. Writes
localization_predictions.tsv and prediction_verification.json.
"""
import csv
import hashlib
import json
import zipfile
from pathlib import Path

HERE = Path("model/reports/r673_prolyl_localization_20261029")
RAW = HERE / "raw"
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
seqs, key = {}, None
for line in (HERE / "target.fasta").read_text().splitlines():
    if line.startswith(">"):
        key = line[1:]
        seqs[key] = ""
    else:
        seqs[key] += line
checks = {}
# TargetP
tp = json.loads((RAW / "targetp/output.json").read_text())
checks["targetp_failed_jobs_zero"] = tp["INFO"]["failedjobs"] == 0 and set(tp["SEQUENCES"]) == set(seqs)
with zipfile.ZipFile(RAW / "targetp/output_all_results.zip") as z:
    checks["targetp_zip_crc_ok"] = z.testzip() is None
    zipped = {n: z.read(n) for n in z.namelist()}
tp_rows = {}
for name, seq in seqs.items():
    pred = RAW / f"targetp/output_{name}_pred.txt"
    checks[f"targetp_pred_equals_zip_{name}"] = zipped[pred.name] == pred.read_bytes()
    residues = [l.split("\t")[1] for l in pred.read_text().splitlines() if l and not l.startswith("#")]
    checks[f"targetp_prefix_{name}"] = "".join(residues) == seq[:len(residues)]
    checks[f"targetp_prefix_length_{name}"] = len(residues)
summary = [l.split("\t") for l in (RAW / "targetp/output_protein_type.txt").read_text().splitlines() if not l.startswith("#")]
for row in summary:
    tp_rows[row[0]] = {"prediction": row[1], "OTHER": float(row[2]), "SP": float(row[3]), "mTP": float(row[4]),
                       "CS": row[5] if len(row) > 5 else ""}
# DeepLoc
dl = json.loads((RAW / "deeploc/results.json").read_text())
checks["deeploc_failed_jobs_zero"] = dl["info"]["failedjobs"] == 0 and set(dl["sequences"]) == set(seqs)
thresholds = dict(zip(dl["Localization"], map(float, dl["Threshold"])))
dl_rows = {r["Protein_ID"]: r for r in csv.DictReader(open(RAW / "deeploc/results_6AE3BCA7002509D23808033B.csv"))}
for name, seq in seqs.items():
    alpha = RAW / f"deeploc/alpha_{name.lower().replace('.', '')}.csv"
    spelled = "".join(r["AA"] for r in csv.DictReader(open(alpha)))
    checks[f"deeploc_full_sequence_{name}"] = spelled == seq
rows = []
for name in seqs:
    d = dl_rows[name]
    above = [loc for loc, t in thresholds.items() if float(d[loc]) > t]
    rows.append({"sequence": name, "length": len(seqs[name]),
                 "targetp_prediction": tp_rows[name]["prediction"], "targetp_OTHER": tp_rows[name]["OTHER"],
                 "targetp_SP": tp_rows[name]["SP"], "targetp_mTP": tp_rows[name]["mTP"],
                 "targetp_cleavage": tp_rows[name]["CS"],
                 "deeploc_localizations": d["Localizations"], "deeploc_signals": d["Signals"],
                 "deeploc_cytoplasm": round(float(d["Cytoplasm"]), 4), "deeploc_mitochondrion": round(float(d["Mitochondrion"]), 4),
                 "deeploc_nucleus": round(float(d["Nucleus"]), 4), "deeploc_above_threshold": ";".join(above)})
with open(HERE / "localization_predictions.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
record = {"jobs": {"TargetP-2.0": {"job_id": "6AE3BC76002509BB05EB4B1C", "organism": "Non-Plant", "format": "long"},
                   "DeepLoc-2.1": {"job_id": "6AE3BCA7002509D23808033B", "encode": "Slow (High-quality, ProtT5)",
                                   "format": "long", "thresholds": thresholds,
                                   "membrane_thresholds": dict(zip(dl["Membrane_types"], dl["Threshold_memtype"]))}},
          "target_fasta_sha256": sha(HERE / "target.fasta"),
          "files": {str(p.relative_to(HERE)): sha(p) for p in sorted(RAW.rglob("*")) if p.is_file()},
          "checks": checks, "all_checks_passed": all(v is True or isinstance(v, int) and not isinstance(v, bool) for v in checks.values())}
(HERE / "prediction_verification.json").write_text(json.dumps(record, indent=1) + "\n")
print(json.dumps(checks, indent=1)); print("all passed:", record["all_checks_passed"])
