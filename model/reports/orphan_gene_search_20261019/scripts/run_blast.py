"""BLAST the chosen query proteins against W29 and check each W29 candidate against Swiss-Prot.

Query sequences come from the local Swiss-Prot DB (blastdbcmd), so the run uses one fixed release.
W29 candidates per reaction: the best hit of each query plus the next best distinct proteins,
at most 5 per reaction, ranked by best bitscore over that reaction's queries.

Usage: python run_blast.py SKILL_DIR DBDIR queries.csv OUTDIR
"""
import os
import subprocess
import sys

import pandas as pd


def acc_of(seqid):
    """BLAST+ 2.17 reports Swiss-Prot ids either as 'sp|ACC|NAME' or as the bare accession."""
    return seqid.split("|")[1] if "|" in seqid else seqid


def main(skill_dir, dbdir, queries_csv, outdir):
    sys.path.insert(0, skill_dir)
    import kernel
    kernel.curation_init(dbdir)
    q = pd.read_csv(queries_csv).fillna("")
    q = q[q.accession != ""]
    accs = sorted(set(q.accession))
    fasta = subprocess.run(["blastdbcmd", "-db", os.path.join(dbdir, "sprot"), "-entry", ",".join(accs)],
                           capture_output=True, text=True, check=True).stdout
    got = {acc_of(l[1:].split()[0]) for l in fasta.splitlines() if l.startswith(">")}
    missing = sorted(set(accs) - got)
    if missing:
        raise SystemExit("query accessions not in local Swiss-Prot: " + ",".join(missing))
    hits = kernel.blast_w29(fasta, evalue=1e-5, max_target=20)
    hits["query_acc"] = hits.qseqid.map(acc_of)
    hits.to_csv(os.path.join(outdir, "blast_w29_raw.csv"), index=False)

    rows, cands = [], []
    for rid, g in q.groupby("reaction", sort=False):
        h = hits[hits.query_acc.isin(set(g.accession))].copy()
        h.insert(0, "reaction", rid)
        rows.append(h)
        if not len(h):
            continue
        best = h.sort_values("bitscore", ascending=False).drop_duplicates("protein_id")
        top_per_query = h.sort_values("bitscore", ascending=False).drop_duplicates("query_acc").protein_id.tolist()
        chosen = list(dict.fromkeys(top_per_query + best.protein_id.tolist()))[:5]
        cands += [(rid, p) for p in chosen]
    allhits = pd.concat(rows) if rows else pd.DataFrame()
    allhits.to_csv(os.path.join(outdir, "blast_w29_hits.csv"), index=False)

    prots = sorted({p for _, p in cands})
    rec = kernel.reciprocal_hits(prots, top=5)
    rec["protein_id"] = rec.qseqid.str.replace(r"^ref\|", "", regex=True).str.strip("|")
    gm = kernel.CURATION_DB["gm"][["protein_id", "locus_tag", "model_gene_id", "in_model", "model_reactions", "product"]]
    rec = rec.merge(gm, on="protein_id", how="left")
    cmap = pd.DataFrame(cands, columns=["reaction", "protein_id"])
    rec = cmap.merge(rec, on="protein_id", how="left")
    rec.to_csv(os.path.join(outdir, "reciprocal_hits.csv"), index=False)
    print("queries", len(accs), "hits", len(allhits), "candidates", len(prots))


if __name__ == "__main__":
    main(*sys.argv[1:5])
