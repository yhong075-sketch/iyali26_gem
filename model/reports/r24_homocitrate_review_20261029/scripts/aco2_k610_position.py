"""Residue aligned to S. cerevisiae Aco2p (P39533) position 610 in Aco1p and the two Y. lipolytica aconitases,
plus LYS4 identity (read-only; homology only).

Fazius 2012 (PMID 23106124) reports an Aco2p K610R mutant with higher aconitate activity. Whether 610 is counted
on the P39533 precursor is an assumption checked here only by P39533[610] being K. Same aligner settings as
../../aconitase_aco2_split_20261029/scripts/align_aconitases.py (BLOSUM62, open -10, extend -0.5, end gaps free).

Usage (from the repo root): PYTHON -I model/reports/r24_homocitrate_review_20261029/scripts/aco2_k610_position.py \
    EVIDENCE_DIR
"""
import sys
from pathlib import Path

import Bio
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

POSITION = 610  # 1-based, on P39533
OTHERS = {"P19414": "Sc_ACO1", "Q6C9P6": "YALI1D11984g", "Q6C5V2": "YALI1E17935g",
          "XP_751171": "Af_AcoA", "XP_750430": "Af_AcoB", "XP_661498": "An_AcoB"}

evidence = Path(sys.argv[1])
seq = {acc: str(SeqIO.read(evidence / f"{acc}.fasta", "fasta").seq)
       for acc in ("P39533", "Q6C791", "P49367", *OTHERS)}
blosum = substitution_matrices.load("BLOSUM62")
aligner = PairwiseAligner(mode="global", substitution_matrix=blosum, open_gap_score=-10, extend_gap_score=-0.5)
aligner.target_end_gap_score = aligner.query_end_gap_score = 0.0

print("biopython", Bio.__version__)
print(f"P39533[{POSITION}] = {seq['P39533'][POSITION - 1]}, context {seq['P39533'][POSITION - 6:POSITION + 5]}")
for acc, label in OTHERS.items():
    aln = aligner.align(seq["P39533"], seq[acc])[0]
    # aln.coordinates maps blocks of aligned residues; find the partner of P39533 index POSITION-1
    partner = None
    for (t0, t1), (q0, q1) in zip(aln.aligned[0], aln.aligned[1]):
        if t0 <= POSITION - 1 < t1:
            partner = q0 + (POSITION - 1 - t0)
    if partner is None:
        print(f"{label} ({acc}): gap at P39533 {POSITION}")
    else:
        print(f"{label} ({acc}): residue {partner + 1} = {seq[acc][partner]}, context {seq[acc][partner - 5:partner + 6]}")

aln = aligner.align(seq["Q6C791"], seq["P49367"])[0]
pairs = list(zip(aln[0], aln[1]))
both = [i for i, (x, y) in enumerate(pairs) if x != "-" and y != "-"]
span = pairs[both[0]: both[-1] + 1]
ident = sum(x == y for x, y in span if x != "-")
print(f"LYS4 Q6C791 vs Sc Lys4p P49367: identity {100 * ident / len(span):.1f}% over span {len(span)}")
