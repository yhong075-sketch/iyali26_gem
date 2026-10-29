"""Crude N-terminal composition of the synthetase paralogs and their yeast counterparts (first 30 residues:
R+K, D+E, net charge, hydrophobic fraction, S+T). A rough presequence indicator only, not a predictor; see
REPORT.md section 7.2. Reads evidence/uniprot_aars_atp_synthase_20261029.tsv.

Usage: python3 model/reports/or_rule_review_r171_aars_20261029/scripts/nterm_heuristic.py
"""
import csv
from pathlib import Path

TASK = Path(__file__).resolve().parents[1]
ORDER = ["P07236", "Q6CHF9", "P04801", "Q6C234", "P38705", "Q6CFL2", "P07284", "Q6C350", "P39965", "Q6C8Z3", "P38708", "Q6C6Z6"]
HYDROPHOBIC = set("AILMFVWY")
rows = {r["Entry"]: r for r in csv.DictReader(open(TASK / "evidence/uniprot_aars_atp_synthase_20261029.tsv"), delimiter="\t")}
print(f"{'acc':8s} {'gene':14s} {'len':>4s} {'R+K':>4s} {'D+E':>4s} {'net':>4s} {'hyd%':>5s} {'S+T':>4s}  first30")
for acc in ORDER:
    r = rows[acc]
    s = r["Sequence"][:30]
    rk, de = sum(c in "RK" for c in s), sum(c in "DE" for c in s)
    gene = r["Gene Names (ORF)"] or r["Gene Names"].split()[0]
    print(f"{acc:8s} {gene:14s} {r['Length']:>4s} {rk:4d} {de:4d} {rk - de:4d} {100 * sum(c in HYDROPHOBIC for c in s) / 30:5.0f} {sum(c in 'ST' for c in s):4d}  {s}")
