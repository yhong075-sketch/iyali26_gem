"""B: reciprocal best hits of the model's mitochondrial malate/citrate carrier genes (PF00153 family panel).

Queries: the RefSeq protein of each model gene's own locus (model_carriers_locus.fasta, fetched by locus tag
because two genes, YALI1B04831g and YALI1E04823g, have no RefSeq annotation in the model). Docstring corrected after
the run (2026-10-14 audit); code unchanged. Forward: each query vs the 34
reviewed S. cerevisiae PF00153 carriers. Reverse: Ctp1, Yhm2, Dic1, Odc1, Odc2, Sfc1, Oac1 vs the 38 Y. lipolytica
CLIB122 PF00153 carriers. Each query is matched to its CLIB122 entry by best SW score. Task B's local Smith-Waterman
(Gotoh, BLOSUM62, gap 11/1) and Task C's FASTA panels, SHA-pinned. A family-panel test, not proteome-wide.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/r831_ctp1_mode_rerun_20261014/orthology_sw.py
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
SW = (RESEARCH / "artifacts/tasks/r132_mito_adh_20261012/sw_blosum62_gotoh.py",
      "04a1f73fc18b96b0d7b3179e62380c438a98cf86574d37b2fb06bec0cb66398c")
PANEL = RESEARCH / "artifacts/tasks/w29_13c_consistency_20261012/sequences"
FASTA = {"queries": (RESEARCH / "artifacts/tasks" / HERE.name / "sequences/model_carriers_locus.fasta",
                     "da57bb6b3cbdedb6b3926a831dfaa1a8ead1547bb4c7021b5db44bbe4442ca0b"),
         "sc": (PANEL / "sc_mcf_reviewed.fasta", "842e9e594307368c016fc8009173adb27c71539d2aec5ed463d9baf6cb6e5aa1"),
         "yl": (PANEL / "yl_clib122_mcf.fasta", "027a613b88610e2e1da4aa5f8ca3a5d165c2f99371310c80802b2aeeae43f7ae")}
# Model gene -> (model reactions, expected S. cerevisiae ortholog from the GPR's annotation/name).
GENES = {"YALI1F33754g": ("R831;R832;R833", "CTP1"), "YALI1B14382g": ("R817", "YHM2"),
         "YALI1D03221g": ("R807;R1386", "ODC1"), "YALI1B04831g": ("R902;R928", "DIC1"),
         "YALI1E41148g": ("R929", "SFC1"), "YALI1E04823g": ("R911", "OAC1")}
REVERSE = ("CTP1", "YHM2", "DIC1", "ODC1", "ODC2", "SFC1", "OAC1")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    assert sha(SW[0]) == SW[1], "SW implementation changed"
    for path, digest in FASTA.values():
        assert sha(path) == digest, path
    spec = importlib.util.spec_from_file_location("sw_impl", SW[0])
    sw = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sw)
    q = sw.read_fasta(FASTA["queries"][0])
    sc = sw.read_fasta(FASTA["sc"][0])
    yl = sw.read_fasta(FASTA["yl"][0])
    gn = {}
    for line in open(FASTA["sc"][0]):
        if line.startswith(">"):
            gn[line[1:].split()[0]] = line.split("GN=")[1].split()[0] if "GN=" in line else None
    sc_ids = {g: k for k, g in gn.items() if g}
    assert all(g in sc_ids for g in REVERSE), [g for g in REVERSE if g not in sc_ids]
    rows, match = [], {}
    for qk, qs in q.items():
        for sk, ss in sc.items():
            rows.append({"direction": "yl_to_sc", "query": qk, "subject": sk, "subject_gene": gn[sk],
                         **sw.sw(qs, ss), "query_len": len(qs), "subject_len": len(ss)})
        best = max(((sw.sw(qs, ys), yk) for yk, ys in yl.items()), key=lambda t: t[0]["score"])
        match[qk] = {"clib122_entry": best[1], "identity": round(best[0]["identity"], 4),
                     "q_cov": round(best[0]["q_cov"], 3), "s_cov": round(best[0]["s_cov"], 3)}
    for g in REVERSE:
        ss = sc[sc_ids[g]]
        for yk, ys in yl.items():
            rows.append({"direction": "sc_to_yl", "query": sc_ids[g], "subject": yk, "subject_gene": g,
                         **sw.sw(ss, ys), "query_len": len(ss), "subject_len": len(ys)})
    fields = ["direction", "query", "subject", "subject_gene", "query_len", "subject_len", "score", "identity",
              "identical", "columns", "q_start", "q_end", "s_start", "s_end", "q_cov", "s_cov"]
    with (HERE / "orthology_sw.tsv").open("x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(r[k], 4) if isinstance(r[k], float) else r[k]) for k in fields})

    def top(rs, n=3):
        return [{"subject": r["subject"], "gene": gn.get(r["subject"]), "score": r["score"],
                 "identity": round(r["identity"], 4), "q_cov": round(r["q_cov"], 3), "s_cov": round(r["s_cov"], 3)}
                for r in rs[:n]]
    summary = {"sw_impl": {"path": str(SW[0]), "sha256": SW[1]},
               "fasta": {k: {"file": str(v[0]), "sha256": v[1]} for k, v in FASTA.items()},
               "query_to_clib122": match, "genes": {}, "reverse": {}}
    clib_of = {qk.split("|")[0]: m["clib122_entry"] for qk, m in match.items()}
    for qk in q:
        gene = qk.split("|")[0]
        fwd = sorted((r for r in rows if r["direction"] == "yl_to_sc" and r["query"] == qk), key=lambda r: -r["score"])
        expected = GENES[gene][1]
        rev = sorted((r for r in rows if r["direction"] == "sc_to_yl" and r["subject_gene"] == expected),
                     key=lambda r: -r["score"])
        summary["genes"][gene] = {"query": qk, "reactions": GENES[gene][0], "expected_sc": expected,
                                  "forward_top3": top(fwd), "forward_best_gene": gn[fwd[0]["subject"]],
                                  "reverse_best_clib122": rev[0]["subject"],
                                  "reciprocal_best_hit": gn[fwd[0]["subject"]] == expected
                                  and rev[0]["subject"] == clib_of[gene]}
    for g in REVERSE:
        rev = sorted((r for r in rows if r["direction"] == "sc_to_yl" and r["subject_gene"] == g),
                     key=lambda r: -r["score"])
        hit = rev[0]["subject"]
        summary["reverse"][g] = {"top3": [{"subject": r["subject"], "score": r["score"],
                                           "identity": round(r["identity"], 4)} for r in rev[:3]],
                                 "best_is_model_gene": [k for k, v in clib_of.items() if v == hit] or None}
    (HERE / "orthology_summary.json").open("x").write(json.dumps(summary, indent=1) + "\n")
    for gene, s in summary["genes"].items():
        print(gene, s["reactions"], "expected", s["expected_sc"], "fwd best", s["forward_best_gene"],
              s["forward_top3"][0]["identity"], "RBH", s["reciprocal_best_hit"])
    for g, s in summary["reverse"].items():
        print("reverse", g, s["top3"][0], "model gene:", s["best_is_model_gene"])


if __name__ == "__main__":
    main()
