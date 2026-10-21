"""Follow-up BLAST for reactions with no W29 hit at E <= 1e-5.

1. Sensitivity: the same queries against W29 at E <= 10, best 3 hits per query, to tell "nothing
   similar" from "a weak hit below the cut-off". The decision cut-off stays at 1e-5.
2. Extra candidate proteins that are not in the current W29 RefSeq proteome (FASTA files, e.g. the
   older KEGG yli model of HGD): reciprocal check against Swiss-Prot, plus each reaction's queries
   aligned against the candidate.

Usage: python followup_blast.py SKILL_DIR DBDIR queries.csv OUTDIR REACTION_IDS EXTRA_FASTA:REACTION ...
"""
import os
import subprocess
import sys

import pandas as pd

COLS = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"


def blastp(query_fa, subject=None, db=None, evalue=10, max_target=3):
    args = ["blastp", "-query", query_fa, "-evalue", str(evalue), "-outfmt", "6 " + COLS]
    args += ["-subject", subject] if subject else ["-db", db, "-max_target_seqs", str(max_target)]
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    df = pd.DataFrame([l.split("\t") for l in out.splitlines()], columns=COLS.split())
    if len(df):
        for c in COLS.split()[2:]:
            df[c] = pd.to_numeric(df[c])
        df["qcov"] = ((df.qend - df.qstart + 1) / df.qlen * 100).round(1)
        df["scov"] = ((df.send - df.sstart + 1) / df.slen * 100).round(1)
    return df


def main(skill_dir, dbdir, queries_csv, outdir, reaction_ids, *extras):
    sys.path.insert(0, skill_dir)
    import kernel
    kernel.curation_init(dbdir)
    q = pd.read_csv(queries_csv).fillna("")
    rids = reaction_ids.split(",")
    q = q[q.reaction.isin(rids) & (q.accession != "")]
    tmp = os.path.join(outdir, "_followup_queries.faa")
    rows = []
    for rid, g in q.groupby("reaction", sort=False):
        fa = subprocess.run(["blastdbcmd", "-db", os.path.join(dbdir, "sprot"), "-entry", ",".join(g.accession)],
                            capture_output=True, text=True, check=True).stdout
        open(tmp, "w").write(fa)
        df = blastp(tmp, db=os.path.join(dbdir, "w29"), evalue=10, max_target=3)
        if len(df):
            df.insert(0, "reaction", rid)
            df["protein_id"] = df.sseqid.str.replace(r"^ref\|", "", regex=True).str.strip("|")
            df = df.merge(kernel.CURATION_DB["gm"][["protein_id", "model_gene_id", "in_model", "model_reactions"]], on="protein_id", how="left")
            rows.append(df)
        else:
            rows.append(pd.DataFrame([dict(reaction=rid, qseqid="(none at E<=10)")]))
    os.remove(tmp)
    pd.concat(rows).to_csv(os.path.join(outdir, "followup_sensitivity_e10.csv"), index=False)

    out = []
    for spec in extras:
        fa, rid = spec.rsplit(":", 1)
        rec = blastp(fa, db=os.path.join(dbdir, "sprot"), evalue=1e-5, max_target=5)
        rec["sp_acc"] = rec.sseqid.str.split("|").str[1]
        if kernel.CURATION_DB.get("titles") is None:
            kernel.CURATION_DB["titles"] = dict(l.rstrip("\n").split("\t", 1) for l in open(os.path.join(dbdir, "sprot_titles.tsv")) if "\t" in l)
        rec["sp_title"] = rec.sp_acc.map(kernel.CURATION_DB["titles"])
        rec.insert(0, "check", "reciprocal_vs_sprot")
        rec.insert(0, "reaction", rid)
        out.append(rec)
        acc = list(q[q.reaction == rid].accession)
        qfa = subprocess.run(["blastdbcmd", "-db", os.path.join(dbdir, "sprot"), "-entry", ",".join(acc)],
                             capture_output=True, text=True, check=True).stdout
        tmpq = os.path.join(outdir, "_followup_q.faa")
        open(tmpq, "w").write(qfa)
        fwd = blastp(tmpq, subject=fa, evalue=10)
        os.remove(tmpq)
        fwd.insert(0, "check", "query_vs_candidate")
        fwd.insert(0, "reaction", rid)
        out.append(fwd)
    if out:
        pd.concat(out).to_csv(os.path.join(outdir, "followup_extra_candidates.csv"), index=False)


if __name__ == "__main__":
    main(*sys.argv[1:])
