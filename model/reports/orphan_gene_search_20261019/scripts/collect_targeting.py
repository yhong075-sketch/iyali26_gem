"""Merge the TargetP-2.0, MitoFates and DeepLoc 2.0 outputs for group 7 and apply the rule fixed in
TASK.md before any predictor was run: predicted mitochondrial if TargetP mTP >= 0.5 or MitoFates
>= 0.385; predicted secretion signal if TargetP SP >= 0.5. DeepLoc is recorded, not part of the rule.

Usage: python collect_targeting.py PREDICTOR_DIR OUT.tsv
"""
import csv
import glob
import os
import sys

REACTIONS = {"YALI1D01489g": "R989 (mi copy of R64)", "YALI1D33889g": "R1125 (mi copy of R306)",
             "YALI1D24649g": "R1292 (mi copy of R615)", "YALI1F28274g": "R1024 (mi), R1027 (ex) copies of R1025/R1026/R2202",
             "YALI1C09285g": "R1337 (cy) vs gene-bearing R751 (ex)", "YALI1D21464g": "R1337 (cy) TH2-like lead"}


def main(pdir, out):
    tp = {}
    for line in open(glob.glob(os.path.join(pdir, "targetp2_*/output_protein_type.txt"))[0]):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        tp[f[0]] = dict(targetp_prediction=f[1], targetp_other=f[2], targetp_sp=f[3], targetp_mtp=f[4], targetp_cs=f[5])
    mf = {r["seq_id"]: r for r in csv.DictReader((l for l in open(glob.glob(os.path.join(pdir, "mitofates_*/results_summary.tsv"))[0]) if not l.startswith("#")), delimiter="\t")}
    dl = {r["Protein_ID"]: r for r in csv.DictReader(open(glob.glob(os.path.join(pdir, "deeploc2.0_*/results_*.csv"))[0]))}
    thr = list(csv.reader(open(glob.glob(os.path.join(pdir, "deeploc2.0_*/thresholds_from_results_page.tsv"))[0]), delimiter="\t"))
    mito_thr = float(thr[1][thr[0].index("Mitochondrion")])
    rows = []
    for sid, t in tp.items():
        gene, prot = sid.split("_", 1)
        m, d = mf[sid], dl[sid]
        mito = float(t["targetp_mtp"]) >= 0.5 or float(m["probability_of_presequence"]) >= 0.385
        crossing = [n for n, ok in (("TargetP", float(t["targetp_mtp"]) >= 0.5), ("MitoFates", float(m["probability_of_presequence"]) >= 0.385),
                                    ("DeepLoc(not in rule)", float(d["Mitochondrion"]) >= mito_thr)) if ok]
        rows.append(dict(gene=gene, refseq=prot, reactions=REACTIONS[gene], **t,
                         mitofates_probability=m["probability_of_presequence"], mitofates_prediction=m["prediction"],
                         mitofates_cleavage=m["cleavage_site"], deeploc_localizations=d["Localizations"], deeploc_signals=d["Signals"],
                         deeploc_mitochondrion=round(float(d["Mitochondrion"]), 4), deeploc_extracellular=round(float(d["Extracellular"]), 4),
                         deeploc_mito_threshold=mito_thr, rule_mitochondrial=mito, rule_secretion_signal=float(t["targetp_sp"]) >= 0.5,
                         mito_predictors_crossing=";".join(crossing)))
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r["gene"], "mito_rule", r["rule_mitochondrial"], "sp_rule", r["rule_secretion_signal"], "crossing", r["mito_predictors_crossing"])


if __name__ == "__main__":
    main(*sys.argv[1:3])
