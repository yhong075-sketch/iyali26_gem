"""Pairwise global alignment of the two Y. lipolytica aconitases against S. cerevisiae Aco1p and Aco2p (read-only).

EMBOSS needle-like settings: BLOSUM62, gap open -10, extend -0.5, end gaps free. Identity and similarity are given
over the alignment length (needle style) and over the aligned span without end gaps.

Usage (from the repo root): .tools/align-venv/bin/python -I \
    model/reports/aconitase_aco2_split_20261029/scripts/align_aconitases.py EVIDENCE_DIR OUT.tsv
"""
import csv
import sys
from itertools import combinations
from pathlib import Path

import Bio
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

LABELS = {"Q6C9P6": "YALI1D11984g", "Q6C5V2": "YALI1E17935g", "P19414": "Sc_ACO1", "P39533": "Sc_ACO2"}

evidence, out = Path(sys.argv[1]), Path(sys.argv[2])
seqs = {acc: str(SeqIO.read(evidence / f"{acc}.fasta", "fasta").seq) for acc in LABELS}

blosum = substitution_matrices.load("BLOSUM62")
aligner = PairwiseAligner(mode="global", substitution_matrix=blosum, open_gap_score=-10, extend_gap_score=-0.5)
aligner.target_end_gap_score = aligner.query_end_gap_score = 0.0

rows = []
for a, b in combinations(LABELS, 2):
    aln = aligner.align(seqs[a], seqs[b])[0]
    top, bottom = aln[0], aln[1]
    pairs = list(zip(top, bottom))
    # aligned span: from the first to the last column where both sequences have a residue
    both = [i for i, (x, y) in enumerate(pairs) if x != "-" and y != "-"]
    span = pairs[both[0]: both[-1] + 1]
    ident = sum(x == y for x, y in pairs if x != "-")
    simil = sum(x != "-" and y != "-" and blosum[x][y] > 0 for x, y in pairs)
    rows.append({
        "seq_a": f"{LABELS[a]} ({a})", "seq_b": f"{LABELS[b]} ({b})", "len_a": len(seqs[a]), "len_b": len(seqs[b]),
        "score": aln.score, "aln_length": len(pairs), "span_length": len(span), "identical": ident,
        "similar": simil, "identity_pct_aln": round(100 * ident / len(pairs), 1),
        "identity_pct_span": round(100 * ident / len(span), 1),
        "similarity_pct_span": round(100 * simil / len(span), 1),
    })

with out.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print("biopython", Bio.__version__)
for r in rows:
    print(f"{r['seq_a']:24} {r['seq_b']:24} score {r['score']:7.1f}  id {r['identity_pct_aln']:5.1f}% (aln)"
          f"  {r['identity_pct_span']:5.1f}% (span {r['span_length']})  sim {r['similarity_pct_span']:5.1f}%")
