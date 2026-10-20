"""Classify every selection row by source (S) / before (B) / after (A) and balance (static)."""
import json, re, csv, sys
from collections import Counter, defaultdict
M = json.load(open('models.json'))
SEL = json.load(open('/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/orphan-reaction-newest-model-b8ddc8/model/curation/selection/metadata_reaction_selection.json'))['reactions']

def parse(f):
    if f is None or f == '' or not re.fullmatch(r'([A-Z][a-z]?\d*)+', f):
        return None
    out = Counter()
    for el, n in re.findall(r'([A-Z][a-z]?)(\d*)', f):
        out[el] += int(n) if n else 1
    return out

def balance(stoich, species):
    tot = Counter(); ch = 0.0
    for mid, c in stoich.items():
        sp = species.get(mid)
        if sp is None: return 'no_species'
        el = parse(sp['formula'])
        if el is None or sp['charge'] is None: return 'unparsable'
        for e, n in el.items(): tot[e] += c * n
        ch += c * sp['charge']
    res = {e: round(v, 6) for e, v in tot.items() if abs(v) > 1e-6}
    if abs(ch) > 1e-6: res['charge'] = round(ch, 6)
    return res or 'balanced'

def norm(st):
    return {k: float(v) for k, v in st.items() if abs(v) > 0}

rows = []
for rid, rule in SEL.items():
    B, A = rule['before'], rule['after']
    src = M['source']['reactions'].get(rid)
    meta = M['metadata']['reactions'].get(rid)
    a56 = M['a56f']['reactions'].get(rid)
    sp = dict(rule['species'])
    for mid in set(B['stoichiometry']) | set(A['stoichiometry']):
        if mid not in sp and mid in M['a56f']['mets']:
            sp[mid] = M['a56f']['mets'][mid]
    for field in rule['fields']:
        b, a = B[field], A[field]
        s = src[{'stoichiometry': 'stoich', 'bounds': 'bounds', 'gpr': 'gpr'}[field]] if src else None
        if field == 'stoichiometry':
            b, a = norm(b), norm(a); s = norm(s) if s is not None else None
        if s is None: cat = 'reaction_not_in_source'
        elif b == s and a != s: cat = 'B=S (builder left source); A new'
        elif b != s and a == s: cat = 'A=S (selection reverts builder change to source)'
        elif b == s and a == s: cat = 'no-op'
        else: cat = 'B!=S and A!=S'
        row = {'reaction': rid, 'field': field, 'category': cat,
               'a56f_equals_after': (norm(a56['stoich']) if field == 'stoichiometry' else a56[{'bounds':'bounds','gpr':'gpr'}[field]]) == a if a56 else None,
               'metadata_equals_after': (norm(meta['stoich']) if field == 'stoichiometry' else meta[{'bounds':'bounds','gpr':'gpr'}[field]]) == a if meta else None}
        if field == 'stoichiometry':
            row['bal_before'] = balance(b, sp); row['bal_after'] = balance(a, sp)
            row['bal_source'] = balance(s, sp) if s is not None else None
            d = {k: (b.get(k), a.get(k)) for k in set(b) | set(a) if b.get(k) != a.get(k)}
            row['diff'] = d
        else:
            row['before'] = b; row['after'] = a; row['source'] = s
        row['superseded_notes'] = list(rule['superseded_notes'])
        rows.append(row)
json.dump(rows, open('classified.json', 'w'), indent=1, default=str)
c = Counter((r['field'], r['category']) for r in rows)
for k, v in sorted(c.items()): print(k, v)
print('a56f==after for all rows:', Counter(r['a56f_equals_after'] for r in rows))
print('metadata==after for all rows:', Counter(r['metadata_equals_after'] for r in rows))
st = [r for r in rows if r['field'] == 'stoichiometry']
def tag(x): return 'balanced' if x == 'balanced' else ('unbalanced' if isinstance(x, dict) else x)
print('stoich balance (before -> after):', Counter((tag(r['bal_before']), tag(r['bal_after'])) for r in st))
