"""Read-only scan of SBML models: ribonucleotide reductase (NDP->dNDP and NTP->dNTP) reactions
by compartment, and every reaction touching a thioredoxin species outside the cytosol or
moving thioredoxin between compartments. Namespace-agnostic; prints GPR from fbc or notes."""
import sys, re, hashlib
import xml.etree.ElementTree as ET
def loc(t): return t.split('}')[-1]
def scan(path):
    sha = hashlib.sha256(open(path,'rb').read()).hexdigest()
    root = ET.parse(path).getroot()
    model = next(e for e in root if loc(e.tag)=='model')
    sp = {}; comp = {}
    for e in model.iter():
        if loc(e.tag)=='compartment': comp[e.get('id')] = e.get('name') or e.get('id')
        if loc(e.tag)=='species': sp[e.get('id')] = ((e.get('name') or e.get('id')), e.get('compartment'))
    gp = {}
    for e in model.iter():
        if loc(e.tag)=='geneProduct':
            gp[[v for k,v in e.attrib.items() if loc(k)=='id'][0]] = [v for k,v in e.attrib.items() if loc(k)=='label'][0]
    def gpr(el):
        t = loc(el.tag)
        if t=='geneProductRef': return gp.get([v for k,v in el.attrib.items() if loc(k)=='geneProduct'][0],'?')
        return '(' + (' %s ' % t).join(gpr(c) for c in el) + ')'
    print(f'#### {path}\n#### sha256 {sha}; species {len(sp)}')
    nd = re.compile(r'^(d?)([ACGU]|T)(DP|TP)\b')
    for r in model.iter():
        if loc(r.tag)!='reaction': continue
        reac=[];prod=[]; rule=''
        for c in r:
            if loc(c.tag)=='listOfReactants': reac=[(s.get('species'),s.get('stoichiometry')) for s in c]
            if loc(c.tag)=='listOfProducts': prod=[(s.get('species'),s.get('stoichiometry')) for s in c]
            if loc(c.tag)=='geneProductAssociation' and len(c): rule=gpr(list(c)[0])
            if loc(c.tag)=='notes' and not rule:
                txt=ET.tostring(c,encoding='unicode'); m=re.search(r'GENE[_ ]ASSOCIATION:\s*([^<]*)',txt)
                if m: rule=m.group(1).strip()
        nm = lambda L: [(sp[s][0] or '').strip() for s,_ in L]
        cm = sorted({sp[s][1] for s,_ in reac+prod})
        rn, pn = nm(reac), nm(prod)
        name = (r.get('name') or '')
        lowall = ' '.join(rn+pn).lower()
        is_rnr = ('ribonucleo' in name.lower() and 'reductase' in name.lower()) or \
                 any(re.fullmatch(r"(ADP|CDP|GDP|UDP|ATP|CTP|GTP|UTP)(\s*\[.*\])?", x) for x in rn) and any(re.fullmatch(r"d(ADP|CDP|GDP|UDP|ATP|CTP|GTP|UTP)(\s*\[.*\])?", x) for x in pn)
        has_trx = 'thioredoxin' in lowall or re.search(r'\btrx', lowall) is not None
        trx_comps = sorted({sp[s][1] for s,_ in reac+prod if 'thioredoxin' in (sp[s][0] or '').lower() or re.search(r'\btrx', (sp[s][0] or '').lower())})
        if is_rnr or (has_trx and (len(trx_comps)>1 or any('n' in c.lower()[:3] or 'nuc' in comp.get(c,'').lower() for c in trx_comps))):
            f = lambda L: ' + '.join(f"{sp[s][0]}[{comp.get(sp[s][1],sp[s][1])}]" for s,_ in L)
            print(f"{r.get('id')}\t{name}\t{f(reac)} -> {f(prod)}\tGPR: {rule}")
for p in sys.argv[1:]: scan(p)
