"""Reciprocal-best-hit summary from fwd (Y. lipolytica RefSeq vs S. cerevisiae) and rev
(S. cerevisiae vs W29) blastp tables; W29 identity of each query comes from the self search."""
import csv, sys
def best(path):
    b = {}
    for r in csv.reader(open(path), delimiter='\t'):
        q, s, bits = r[0], r[1], float(r[11])
        if q not in b or bits > b[q][1]: b[q] = (s, bits, r)
    return b
fwd, rev, self_ = best(sys.argv[1]), best(sys.argv[2]), best(sys.argv[3])
allrev = [r for r in csv.reader(open(sys.argv[2]), delimiter='\t')]
w = csv.writer(sys.stdout, delimiter='\t')
w.writerow(['query','w29_identical','sce_hit','pident','aln_len','qlen','slen','evalue','bitscore','qcovs','sce_hit_best_w29','reciprocal'])
for r in csv.reader(open(sys.argv[1]), delimiter='\t'):
    q, s = r[0], r[1]
    w29 = self_[q][0]
    if s not in rev: continue
    rb = rev[s][0]
    if r[0] in fwd and fwd[q][0] == s: tag = 'yes' if rb == w29 else 'no'
    else: tag = 'no (not forward best)'
    w.writerow([q, w29, s.split('|')[2], r[2], r[3], r[4], r[5], r[10], r[11], r[12], rb, tag])
