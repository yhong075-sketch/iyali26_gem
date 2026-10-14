"""Task C (C1/C2): reciprocal best hits within the mitochondrial carrier family (Pfam PF00153).

Forward: the three model carrier genes (RefSeq) vs every reviewed S. cerevisiae PF00153 protein.
Reverse: Sc CTP1, YHM2 and CRC1 vs every Y. lipolytica CLIB122 PF00153 protein (UniProt proteome
UP000001300). Each CLIB122 entry is matched to its model gene by near-identity to the RefSeq query.
Alignment: Task B's local Smith-Waterman (Gotoh, BLOSUM62, gap 11/1), SHA-pinned, imported by path.
This implementation has not been cross-checked against BLAST/EMBOSS.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/w29_13c_consistency_20261012/orthology_sw.py
Writes orthology_sw.tsv and orthology_summary.json here. Never overwrites.
"""
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESEARCH = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve()
SEQ = RESEARCH / "artifacts/tasks" / HERE.name / "sequences"
SW = (RESEARCH / "artifacts/tasks/r132_mito_adh_20261012/sw_blosum62_gotoh.py",
      "04a1f73fc18b96b0d7b3179e62380c438a98cf86574d37b2fb06bec0cb66398c")
FASTA = {"queries": ("yali_queries_refseq.fasta", "130eab03173d6e57a02e89f1db1b80ec2cf75bba2a32082759a93e50987c4bbd"),
         "sc": ("sc_mcf_reviewed.fasta", "842e9e594307368c016fc8009173adb27c71539d2aec5ed463d9baf6cb6e5aa1"),
         "yl": ("yl_clib122_mcf.fasta", "027a613b88610e2e1da4aa5f8ca3a5d165c2f99371310c80802b2aeeae43f7ae")}
QUERIES = {"XP_505902.1": ("YALI1F33754g", "CTP1"), "XP_500732.1": ("YALI1B14382g", "YHM2"),
           "XP_501359.1": ("YALI1C03415g", "CRC1")}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    assert sha(SW[0]) == SW[1], "SW implementation changed"
    for name, digest in FASTA.values():
        assert sha(SEQ / name) == digest, name
    spec = importlib.util.spec_from_file_location("sw_impl", SW[0])
    sw = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sw)
    q = sw.read_fasta(SEQ / FASTA["queries"][0])
    sc = sw.read_fasta(SEQ / FASTA["sc"][0])
    yl = sw.read_fasta(SEQ / FASTA["yl"][0])
    rows = []
    for qk, qs in q.items():
        for sk, ss in sc.items():
            r = sw.sw(qs, ss)
            rows.append({"direction": "yl_to_sc", "query": qk, "subject": sk, **r, "query_len": len(qs), "subject_len": len(ss)})
    headers = {}
    for line in open(SEQ / FASTA["sc"][0]):
        if line.startswith(">"):
            headers[line[1:].split()[0]] = line.split("GN=")[1].split()[0] if "GN=" in line else None
    sc_ids = {gn: k for k, gn in headers.items() if gn}
    # Match each RefSeq query to its CLIB122 entry.
    match = {}
    for qk, qs in q.items():
        best = max(((sw.sw(qs, ys), yk) for yk, ys in yl.items()), key=lambda t: t[0]["score"])
        match[qk] = {"clib122_entry": best[1], "identity": best[0]["identity"], "q_cov": best[0]["q_cov"],
                     "s_cov": best[0]["s_cov"]}
    for gn in ("CTP1", "YHM2", "CRC1"):
        ss = sc[sc_ids[gn]]
        for yk, ys in yl.items():
            r = sw.sw(ss, ys)
            rows.append({"direction": "sc_to_yl", "query": sc_ids[gn], "subject": yk, **r, "query_len": len(ss),
                         "subject_len": len(ys)})
    fields = ["direction", "query", "subject", "query_len", "subject_len", "score", "identity", "identical",
              "columns", "q_start", "q_end", "s_start", "s_end", "q_cov", "s_cov"]
    with (HERE / "orthology_sw.tsv").open("x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(r[k], 4) if isinstance(r[k], float) else r[k]) for k in fields})
    summary = {"sw_impl": {"path": str(SW[0]), "sha256": SW[1]},
               "fasta": {k: {"file": v[0], "sha256": v[1]} for k, v in FASTA.items()},
               "refseq_to_clib122": match, "pairs": {}}
    for qk, (model_gene, expected) in QUERIES.items():
        fwd = sorted((r for r in rows if r["direction"] == "yl_to_sc" and r["query"] == qk),
                     key=lambda r: -r["score"])
        rev = sorted((r for r in rows if r["direction"] == "sc_to_yl" and r["query"] == sc_ids[expected]),
                     key=lambda r: -r["score"])
        top = lambda rs, n: [{"subject": r["subject"], "score": r["score"], "identity": round(r["identity"], 4),
                              "q_cov": round(r["q_cov"], 3), "s_cov": round(r["s_cov"], 3)} for r in rs[:n]]
        summary["pairs"][model_gene] = {
            "refseq": qk, "expected_sc_ortholog": expected, "forward_top3": top(fwd, 3), "reverse_top3": top(rev, 3),
            "forward_best_is_expected": headers[fwd[0]["subject"]] == expected,
            "reverse_best_is_this_gene": rev[0]["subject"] == match[qk]["clib122_entry"],
        }
    (HERE / "orthology_summary.json").open("x").write(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary["pairs"], indent=1))
    print(json.dumps(match, indent=1))


if __name__ == "__main__":
    main()
