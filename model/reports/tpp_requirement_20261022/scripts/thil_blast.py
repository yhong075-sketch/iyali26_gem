"""Search the W29 proteome for a thiamine-monophosphate kinase (ThiL, EC 2.7.4.16) (read-only).

Every Swiss-Prot entry titled "Thiamine-monophosphate kinase" is a query. Forward blastp against W29
at E <= 10 (all hits kept; the decision cut-off is 1e-5). Each W29 protein with a hit at E <= 1e-5
(at most 5, by best bitscore) is searched back against Swiss-Prot (top 5 hits).

Usage: python thil_blast.py BLAST_BIN_DIR DBDIR OUTDIR
DBDIR must be a path without spaces (BLAST+ splits -db at spaces); the orphan_gene_search_20261019
databases are used through a symlink.
"""
import csv
import os
import subprocess
import sys

bin_dir, dbdir, outdir = sys.argv[1:4]
FIELDS = "qseqid sseqid pident length qlen slen evalue bitscore"
DECISION_EVALUE = 1e-5


def run(tool, *args, stdin=None):
    return subprocess.run([os.path.join(bin_dir, tool), *args], input=stdin, capture_output=True,
                          text=True, check=True).stdout


def acc(seqid):
    return seqid.split("|")[1] if "|" in seqid else seqid


titles = {}
with open(os.path.join(dbdir, "sprot_titles.tsv")) as handle:
    for line in handle:
        accession, title = line.rstrip("\n").split("\t", 1)
        titles[accession] = title
queries = sorted(a for a, t in titles.items() if t.startswith("Thiamine-monophosphate kinase "))
gene_map = {}
with open(os.path.join(dbdir, "w29_gene_map.csv")) as handle:
    for row in csv.DictReader(handle):
        gene_map[row["protein_id"]] = row

fasta = run("blastdbcmd", "-db", os.path.join(dbdir, "sprot"), "-entry", ",".join(queries))
forward = run("blastp", "-db", os.path.join(dbdir, "w29"), "-evalue", "10", "-max_target_seqs", "20",
              "-outfmt", f"6 {FIELDS}", stdin=fasta)
rows = []
for line in forward.splitlines():
    q, s, pident, length, qlen, slen, evalue, bits = line.split("\t")
    protein = s.replace("ref|", "").strip("|")
    info = gene_map.get(protein, {})
    rows.append({"query": acc(q), "query_title": titles.get(acc(q), ""), "w29_protein": protein,
                 "locus_tag": info.get("locus_tag", ""), "model_gene_id": info.get("model_gene_id", ""),
                 "product": info.get("product", ""), "pident": float(pident), "length": int(length),
                 "qcov": round(int(length) / int(qlen), 3), "evalue": float(evalue), "bitscore": float(bits)})
rows.sort(key=lambda r: (-r["bitscore"], r["query"]))
with open(os.path.join(outdir, "thil_blast_w29.csv"), "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["query"], lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

best = {}
for row in rows:
    if row["evalue"] <= DECISION_EVALUE and row["w29_protein"] not in best:
        best[row["w29_protein"]] = row["bitscore"]
candidates = sorted(best, key=lambda p: -best[p])[:5]
reciprocal = []
if candidates:
    w29_fasta = run("blastdbcmd", "-db", os.path.join(dbdir, "w29"), "-entry", ",".join(candidates))
    back = run("blastp", "-db", os.path.join(dbdir, "sprot"), "-evalue", "1e-5", "-max_target_seqs", "5",
               "-outfmt", f"6 {FIELDS}", stdin=w29_fasta)
    for line in back.splitlines():
        q, s, pident, length, qlen, slen, evalue, bits = line.split("\t")
        reciprocal.append({"w29_protein": q.replace("ref|", "").strip("|"), "sprot": acc(s),
                           "sprot_title": titles.get(acc(s), ""), "pident": float(pident),
                           "evalue": float(evalue), "bitscore": float(bits)})
with open(os.path.join(outdir, "thil_reciprocal.csv"), "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["w29_protein", "sprot", "sprot_title", "pident", "evalue", "bitscore"],
                            lineterminator="\n")
    writer.writeheader()
    writer.writerows(reciprocal)

print("queries", len(queries), "forward hits (E<=10)", len(rows),
      "W29 proteins at E<=1e-5", len(best), "reciprocal rows", len(reciprocal))
for row in rows[:15]:
    print(f"{row['query']:8s} {row['w29_protein']:16s} {row['model_gene_id'] or row['locus_tag']:14s} "
          f"id {row['pident']:5.1f} qcov {row['qcov']:.2f} E {row['evalue']:.2g} bits {row['bitscore']:.0f} {row['product'][:50]}")
for row in reciprocal:
    print("back", row["w29_protein"], row["sprot"], f"{row['pident']:.1f}", f"{row['evalue']:.2g}", row["sprot_title"][:70])
