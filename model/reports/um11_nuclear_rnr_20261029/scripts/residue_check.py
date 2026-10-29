"""Map catalytic/iron-binding residues of S. cerevisiae RNR subunits (UniProt FT) onto
Y. lipolytica sequences by global pairwise alignment (BLOSUM62, gap open -10, extend -0.5)."""
import sys
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
sce = {r.id.split('|')[1]: str(r.seq) for r in SeqIO.parse(sys.argv[1], 'fasta')}
yli = {r.id.split('|')[1]: str(r.seq) for r in SeqIO.parse(sys.argv[2], 'fasta')}
al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score = -10; al.extend_gap_score = -0.5; al.mode = 'global'
# positions from UniProt FT of P21524 (Rnr1), P09938 (Rnr2)
checks = {
 'P21524': [(218,'C','redox-active disulfide / H transfer'),(426,'N','proton acceptor'),(428,'C','cysteine radical'),
            (430,'E','proton acceptor'),(443,'C','redox-active disulfide / H transfer'),(741,'Y','electron transfer'),
            (742,'Y','electron transfer'),(883,'C','thioredoxin/glutaredoxin interaction'),(886,'C','thioredoxin/glutaredoxin interaction')],
 'P09938': [(145,'D','Fe'),(176,'E','Fe'),(179,'H','Fe'),(183,'Y','tyrosyl radical (ACT_SITE)'),(239,'E','Fe'),(273,'E','Fe'),(276,'H','Fe')],
}
pairs = [('P21524','YALI1_F39330g'),('P21524','YALI1_D04261g'),('P21672','YALI1_F39330g'),('P09938','YALI1_B05072g'),('P49723','YALI1_B05072g')]
for s, y in pairs:
    a = al.align(sce[s], yli[y])[0]
    # build map sce position (1-based) -> yli residue
    m = {}
    for (s0, s1), (y0, y1) in zip(*a.aligned):
        for k in range(s1 - s0):
            m[s0 + k + 1] = (y0 + k + 1, yli[y][y0 + k])
    print(f'== {s} vs {y}: score {a.score:.0f}')
    ref = 'P21524' if s in ('P21524','P21672') else 'P09938'
    if s == ref:
        for pos, exp, note in checks[s]:
            got = m.get(pos)
            print(f'  {s} {sce[s][pos-1]}{pos} ({note}) -> {y} ' + (f'{got[1]}{got[0]}' if got else 'GAP'))
    # C-terminal tail
    print(f'  {s} C-term 15 aa: {sce[s][-15:]}   {y} C-term 15 aa: {yli[y][-15:]}')
# Rnr4 vs Rnr2 residues at Rnr2 Fe positions (to show Rnr4 lacks them), via Rnr2-Rnr4 alignment
a = al.align(sce['P09938'], sce['P49723'])[0]
m = {}
for (s0, s1), (y0, y1) in zip(*a.aligned):
    for k in range(s1 - s0):
        m[s0 + k + 1] = (y0 + k + 1, sce['P49723'][y0 + k])
print('== P09938 (Rnr2) vs P49723 (Rnr4) at Rnr2 Fe/radical positions')
for pos, exp, note in checks['P09938']:
    got = m.get(pos)
    print(f'  Rnr2 {sce["P09938"][pos-1]}{pos} ({note}) -> Rnr4 ' + (f'{got[1]}{got[0]}' if got else 'GAP'))
