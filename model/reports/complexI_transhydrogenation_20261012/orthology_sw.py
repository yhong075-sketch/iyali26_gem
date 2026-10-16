"""Task D step 3: S. cerevisiae isoform assignment of the Y. lipolytica GPR genes, by local Smith-Waterman.

Task B's implementation (Gotoh, BLOSUM62, gap 11/1; SHA-pinned, imported by path; not cross-checked against BLAST).
All-vs-all within each enzyme family panel (forward: each Y. lipolytica protein vs the family's S. cerevisiae
proteins; reverse: each S. cerevisiae protein vs the family's Y. lipolytica proteins). This assigns the closest
isoform within the family; it is not a proteome-wide reciprocal-best-hit test.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/complexI_transhydrogenation_20261012/orthology_sw.py
Writes orthology_sw.tsv and orthology_summary.json here. Never overwrites.
"""
import csv
import importlib.util
import io
import json
from datetime import datetime, timezone

from common import ART, HERE, RESEARCH, sha, write_new

SW = (RESEARCH / "artifacts/tasks/r132_mito_adh_20261012/sw_blosum62_gotoh.py",
      "04a1f73fc18b96b0d7b3179e62380c438a98cf86574d37b2fb06bec0cb66398c")
MANIFEST = HERE / "sequences_manifest_v3.json"
FAMILIES = {
    "MDH": (["YALI1D20676g", "YALI1E17214g"], ["MDH1", "MDH2", "MDH3"]),
    "mito_carrier": (["YALI1F33754g", "YALI1B14382g"], ["CTP1", "YHM2"]),
    "MTHFD": (["YALI1B15928g", "YALI1F38386g", "YALI1E01589g"], ["MTD1", "ADE3", "MIS1"]),
    "GDH": (["YALI1E11943g", "YALI1F23664g"], ["GDH1", "GDH2", "GDH3"]),
    "ALD": (["YALI1C04093g", "YALI1D10197g", "YALI1E00588g", "YALI1F06858g", "YALI1B02194g", "YALI1F31215g"],
            ["ALD2", "ALD3", "ALD4", "ALD5", "ALD6"]),
    "IDH": (["YALI1F06197g", "YALI1E06015g", "YALI1D08082g", "YALI1D13238g"],
            ["IDP1", "IDP2", "IDP3", "IDH1", "IDH2", "LYS12"]),
}


def main():
    assert sha(SW[0]) == SW[1], "SW implementation changed"
    manifest = json.loads(MANIFEST.read_text())
    files = {k: manifest[k] for k in ("yali_task_d_v3.fasta", "sc_task_d_v3.fasta")}
    for f in files.values():
        assert sha(f["path"]) == f["sha256"], f["path"]
    spec = importlib.util.spec_from_file_location("sw_impl", SW[0])
    sw = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sw)
    yl = sw.read_fasta(files["yali_task_d_v3.fasta"]["path"])
    sc = sw.read_fasta(files["sc_task_d_v3.fasta"]["path"])
    rows, summary = [], {}
    for family, (ygenes, sgenes) in FAMILIES.items():
        scores = {}
        for y in ygenes:
            for s in sgenes:
                r = sw.sw(yl[y], sc[s])
                scores[(y, s)] = r
                rows.append({"family": family, "yali": y, "yali_len": len(yl[y]), "sc": s, "sc_len": len(sc[s]),
                             **{k: r[k] for k in ("score", "identity", "columns", "q_start", "q_end", "s_start", "s_end",
                                                  "q_cov", "s_cov")}})
        for y in ygenes:
            best = max(sgenes, key=lambda s: scores[(y, s)]["score"])
            back = max(ygenes, key=lambda yy: scores[(yy, best)]["score"])
            summary[y] = {"family": family, "best_sc": best, "score": scores[(y, best)]["score"],
                          "identity": round(scores[(y, best)]["identity"], 3), "q_cov": round(scores[(y, best)]["q_cov"], 3),
                          "reciprocal_within_panel": back == y,
                          "all": {s: [scores[(y, s)]["score"], round(scores[(y, s)]["identity"], 3)] for s in sgenes}}
            print(y, summary[y]["best_sc"], summary[y]["score"], summary[y]["identity"], summary[y]["reciprocal_within_panel"],
                  flush=True)
    buffer = io.StringIO()
    w = csv.DictWriter(buffer, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    write_new(HERE / "orthology_sw.tsv", buffer.getvalue())
    write_new(HERE / "orthology_summary.json", json.dumps(
        {"run_utc": datetime.now(timezone.utc).isoformat(), "sw": {"path": str(SW[0]), "sha256": SW[1]},
         "sequences": files, "families": FAMILIES, "summary": summary}, indent=1) + "\n")


if __name__ == "__main__":
    main()
