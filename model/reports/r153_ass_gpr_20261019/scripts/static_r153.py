"""Static (no LP) listing around R153 in an SBML model."""
import sys, re, hashlib
path = sys.argv[1]
raw = open(path, 'rb').read(); s = raw.decode()
print("model", path.split('/')[-1], hashlib.sha256(raw).hexdigest())
attr = lambda tag, a: (re.search(r'\b%s="([^"]*)"' % a, tag) or [None, None])[1]
sp = {}
for m in re.finditer(r'<species\b[^>]*>', s):
    t = m.group(0); sp[attr(t, 'id')] = (attr(t, 'name'), attr(t, 'compartment'), attr(t, 'fbc:chemicalFormula'), attr(t, 'fbc:charge'))
par = dict(re.findall(r'<parameter [^>]*id="([^"]+)"[^>]*value="([^"]+)"', s))
rx = []
for m in re.finditer(r'<reaction [^>]*>.*?</reaction>', s, re.S):
    b = m.group(0); head = re.match(r'<reaction [^>]*>', b).group(0)
    def side(tag):
        t = re.search(r'<%s>(.*?)</%s>' % (tag, tag), b, re.S)
        return re.findall(r'species="([^"]+)" stoichiometry="([^"]+)"', t.group(1)) if t else []
    gpa = re.search(r'<fbc:geneProductAssociation>(.*?)</fbc:geneProductAssociation>', b, re.S)
    rx.append(dict(id=attr(head, 'id'), name=attr(head, 'name'), lb=par.get(attr(head, 'fbc:lowerFluxBound')), ub=par.get(attr(head, 'fbc:upperFluxBound')),
                   r=side('listOfReactants'), p=side('listOfProducts'), gpa=re.sub(r'\s+', ' ', gpa.group(1)) if gpa else '',
                   ec=sorted(set(re.findall(r'ec-code/([0-9.\-n]+)', b)))))
print("reactions", len(rx), "species", len(sp))
def nm(k): return re.sub(r'^M_', '', k).replace('__91__', '[').replace('__93__', ']')
def fmt(L): return " + ".join(f"{c} {nm(k)}" for k, c in L)
def show(r):
    g = re.findall(r'geneProduct="G_([^"]+)"', r['gpa']); op = 'and' if 'fbc:and' in r['gpa'] else ('or' if 'fbc:or' in r['gpa'] else '')
    print(f"  {r['id']} | {r['name']} | [{r['lb']},{r['ub']}] | GPR({op}) {g} | EC {r['ec']}\n       {fmt(r['r'])} -> {fmt(r['p'])}")
print("== reactions with YALI1D17462g in GPR")
for r in rx:
    if 'G_YALI1D17462g' in r['gpa']: show(r)
for base in sys.argv[2:]:
    ks = sorted(k for k in sp if re.match(r'(M_)?%s__91__' % base, k))
    print(f"### {base}: " + "; ".join(f"{nm(k)}={sp[k][0]} {sp[k][2]} z={sp[k][3]}" for k in ks))
    for r in rx:
        if any(x[0] in ks for x in r['r'] + r['p']): show(r)
