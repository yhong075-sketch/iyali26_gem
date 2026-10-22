#!/usr/bin/env bash
# Reciprocal best hits between S. cerevisiae iron-uptake proteins and the Y. lipolytica CLIB122 proteome
# (read-only; BLAST+ 2.17.0 from the seqsearch conda env). Usage: rbh_iron.sh BLAST_DIR WORK_DIR
# WORK_DIR must have no spaces (BLAST+ cannot open databases under paths with spaces); results are copied to BLAST_DIR/out.
set -euo pipefail
B="$1"; W="$2"; BIN=/opt/miniconda3/envs/seqsearch/bin; mkdir -p "$W/db" "$B/out"
cp "$B/proteomes/"*.fasta "$W/"; cd "$W"; P=.; O=.
# FET3 FTR1 FET4 SMF1 SMF2 MRS3 MRS4 FRE1 FRE2 FET5 FTH1 PHO4
Q="P38993 P40088 P40988 P38925 P38778 P10566 P23500 P32791 P36033 P43561 P38310 P07270"
"$BIN/makeblastdb" -in "$P/yli_UP000001300.fasta" -dbtype prot -out "$O/db/yli" > /dev/null
"$BIN/makeblastdb" -in "$P/sce_UP000002311_reviewed.fasta" -dbtype prot -out "$O/db/sce" > /dev/null
python3 - "$P/sce_UP000002311_reviewed.fasta" "$O/queries_sce.fasta" $Q <<'PY'
import sys
src, out, ids = sys.argv[1], sys.argv[2], sys.argv[3:]
keep, seqs, cur = set(ids), {}, None
for line in open(src):
    if line.startswith(">"):
        acc = line.split("|")[1]; cur = acc if acc in keep else None
        if cur: seqs[cur] = [line]
    elif cur:
        seqs[cur].append(line)
missing = keep - set(seqs)
assert not missing, missing
open(out, "w").writelines(l for acc in ids for l in seqs[acc])
PY
FMT="6 qseqid sseqid pident length qlen slen evalue bitscore qcovs"
"$BIN/blastp" -query "$O/queries_sce.fasta" -db "$O/db/yli" -evalue 1e-5 -max_target_seqs 5 -outfmt "$FMT" > "$O/sce_vs_yli.tsv"
cut -f2 "$O/sce_vs_yli.tsv" | sort -u > "$O/yli_hits.ids"
# Also the current R2193 gene, Q6C2U3 (YALI0F05126g).
python3 - "$P/yli_UP000001300.fasta" "$O/yli_hits.ids" "$O/queries_yli.fasta" Q6C2U3 <<'PY'
import sys
src, ids_path, out, extra = sys.argv[1:5]
keep = {l.strip().split("|")[1] if "|" in l else l.strip() for l in open(ids_path) if l.strip()} | {extra}
seqs, cur = {}, None
for line in open(src):
    if line.startswith(">"):
        acc = line.split("|")[1]; cur = acc if acc in keep else None
        if cur: seqs[cur] = [line]
    elif cur:
        seqs[cur].append(line)
open(out, "w").writelines(l for acc in sorted(seqs) for l in seqs[acc])
PY
"$BIN/blastp" -query "$O/queries_yli.fasta" -db "$O/db/sce" -evalue 1e-5 -max_target_seqs 3 -outfmt "$FMT" > "$O/yli_vs_sce.tsv"
cp "$O/sce_vs_yli.tsv" "$O/yli_vs_sce.tsv" "$O/queries_sce.fasta" "$O/queries_yli.fasta" "$B/out/"
