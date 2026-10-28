"""Reusable Escher map building for gap maps: layered chain layout (hub repeats, folding), section packing, factory calls.
Copied from the gap census (artifacts/tasks/gap_metabolites_ponce2013_20261026/escher/scripts/maplib.py); only the
escher_kernel import path changed (now this folder)."""
import os, sys, re, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from escher_kernel import escher_map_factory
COMP = {'C_cy': 'cy', 'C_mi': 'mi', 'C_er': 'er', 'C_go': 'go', 'C_ex': 'ex', 'C_nu': 'nu', 'C_va': 'va', 'C_en': 'en', 'C_em': 'em', 'C_gm': 'gm', 'C_lp': 'lp', 'C_mm': 'mm', 'C_pe': 'pe', 'C_vm': 'vm', 'c_va': 'va'}
CLASS_VALUE = {'RNP': 1, 'RNC': 2, 'RNP+RNC (single reversible reaction)': 3, 'isolated': 4}
CAUSE_SHORT = {'chain disconnected at several roots (needs >1 fix)': 'several roots (>1 fix)', 'dead-end root metabolite': 'dead-end root',
               'cofactor pseudo-gap (cofactor not drained by biomass)': 'cofactor pseudo-gap', 'closed in model (0,0)': 'CLOSED (0,0)',
               'stoichiometric coupling / pseudo-gap (no root)': 'coupling, no root', 'missing transport': 'missing transport',
               'blocked by a closed (0,0) reaction': 'behind a closed reaction', 'missing exchange': 'missing exchange'}
CURRENCY_DEGREE, HUB, MEGA, DX, DY, FOLD = 14, 5, 12, 480, 170, 16

def met_value(m):
    if not m.get('gap'): return 0
    if m.get('cls') in CLASS_VALUE: return CLASS_VALUE[m['cls']]
    return 6 if m.get('cls', '') == '' else 5

def short(n, k=26):
    n = re.sub(r"_(?:[A-Z][a-z]?\d*)+(?:\(.*\))?.*$", "", n) if "_" in n else n
    n = n.replace('_p+1', '').rstrip('_ ').strip()
    return n if len(n) <= k else n[:k - 1] + '…'

def layered(G):
    C = nx.condensation(G); layer = {}
    for scc in nx.topological_sort(C):
        members = C.nodes[scc]['members']
        ext = [layer[u] for v in members for u in G.predecessors(v) if u not in members]
        base = max(ext) + 1 if ext else 0
        if len(members) == 1: layer[next(iter(members))] = base; continue
        sub = G.subgraph(members); starts = [v for v in members if any(u not in members for u in G.predecessors(v))] or [sorted(members)[0]]
        depth = {s: 0 for s in starts}; frontier = list(starts)
        while frontier:
            nxt = []
            for u in frontier:
                for v in sub.successors(u):
                    if v not in depth: depth[v] = depth[u] + 1; nxt.append(v)
            frontier = nxt
        for v in members: layer[v] = base + depth.get(v, 0)
    for v in G.nodes:
        if G.in_degree(v) == 0 and G.out_degree(v) >= 1: layer[v] = min(layer[u] for u in G.successors(v)) - 1
    for v in G.nodes:
        if G.out_degree(v) == 0 and G.in_degree(v) >= 1: layer[v] = max(layer[u] for u in G.predecessors(v)) + 1
    remap = {l: i for i, l in enumerate(sorted(set(layer.values())))}; layer = {v: remap[l] for v, l in layer.items()}
    L = {}
    for v, l in layer.items(): L.setdefault(l, []).append(v)
    for l in L: L[l].sort()
    pos = {v: i for l in L for i, v in enumerate(L[l])}
    def bary(v, up):
        nb = [pos[u] for u in (G.predecessors(v) if up else G.successors(v)) if (layer[u] < layer[v] if up else layer[u] > layer[v])]
        return statistics.mean(nb) if nb else pos[v]
    for it in range(8):
        for l in sorted(L)[1:]: L[l].sort(key=lambda v: bary(v, True)); pos.update({v: i for i, v in enumerate(L[l])})
        for l in sorted(L)[-2::-1]: L[l].sort(key=lambda v: bary(v, False)); pos.update({v: i for i, v in enumerate(L[l])})
    x = {v: pos[v] * DX for v in layer}
    for it in range(6):
        for l in (sorted(L) if it % 2 == 0 else sorted(L, reverse=True)):
            want = {}
            for v in L[l]:
                nb = [x[u] for u in list(G.predecessors(v)) + list(G.successors(v)) if layer[u] != layer[v]]
                want[v] = statistics.mean(nb) if nb else x[v]
            cur = -1e9
            for v in sorted(L[l], key=lambda v: (want[v], pos[v])): x[v] = max(want[v], cur + DX); cur = x[v]
            shift = statistics.mean(x[v] - want[v] for v in L[l])
            for v in L[l]: x[v] -= shift
    xmin = min(x.values()); return layer, {v: val - xmin for v, val in x.items()}

def fold(G, layer, x, is_met):
    nl = max(layer.values()) + 1
    if nl <= FOLD + 2: return layer, x, [], set()
    cuts, c = [], FOLD
    while c < nl:
        while c < nl and not all(is_met(v) for v in layer if layer[v] == c): c += 1
        if c < nl: cuts.append(c)
        c += FOLD
    width = max(x.values()) + DX * 1.6
    col = {v: sum(1 for cc in cuts if l >= cc) for v, l in layer.items()}
    start = {k: (cuts[k - 1] if k else 0) for k in range(len(cuts) + 1)}
    newl = {v: l - start[col[v]] for v, l in layer.items()}; newx = {v: x[v] + col[v] * width for v in layer}
    dups, redirect = [], set()
    for k, cc in enumerate(cuts, 1):
        for v in [v for v in layer if layer[v] == cc]:
            preds = [r for r in G.predecessors(v) if layer[r] < cc]
            if not preds: continue
            dups.append((v, newx[v] - width, cc - start[k - 1]))
            for r in preds: redirect.add((r, v))
    return newl, newx, dups, redirect

def build_map(map_id, title, rxns, mets, sections, maxw=9000, legend=(), cause_label=None, dx=480):
    global DX
    DX = dx
    """sections: [(label, [rxn ids])]. rxns[rid] = {name, lb, ub, gpr, genes, subsystem, mets{mid: coef}, boundary?, sublabel?};
    mets[mid] = {name, comp, degree, gap, cls}. Returns (map_json, deco_json, met_data, W, H, placed)."""
    primary = set()
    for rid, r in rxns.items():
        for mid in r['mets']:
            if mets[mid].get('gap') or (len(r['mets']) <= MEGA and mets[mid].get('degree', 1) <= CURRENCY_DEGREE): primary.add(mid)
    fmt = {mid: dict(name=m['name'], comp=m['comp']) for mid, m in mets.items()}
    frx, rp = {}, {}
    for rid, r in rxns.items():
        reac = [m for m, c in r['mets'].items() if c < 0 and m in primary]; prod = [m for m, c in r['mets'].items() if c > 0 and m in primary]
        for side, lst in (('r', reac), ('p', prod)):
            if not lst:
                cand = sorted([m for m, c in r['mets'].items() if (c < 0 if side == 'r' else c > 0)], key=lambda m: mets[m].get('degree', 1))
                if cand: lst.append(cand[0]); primary.add(cand[0])
        rmets = dict(r['mets'])
        if not reac or not prod:
            pid = f'pseudo_{rid}'; fmt[pid] = dict(name=r.get('pseudo_name') or '(outside)', comp='C_ex'); primary.add(pid)
            if not reac: reac.append(pid); rmets[pid] = -1
            else: prod.append(pid); rmets[pid] = 1
        sub = r.get('sublabel') or (cause_label(r) if cause_label else CAUSE_SHORT.get(r.get('cause', ''), r.get('cause', '')))
        frx[rid] = dict(name=r['name'], lb=r['lb'], ub=r['ub'], gpr=r['gpr'], genes=r['genes'], subsystem=sub, mets=rmets); rp[rid] = (reac, prod)
    sec_of = {rid: i for i, (lab, ids) in enumerate(sections) for rid in ids}
    # hub repeats are decided per section
    keymid, G = {}, nx.DiGraph()
    for si, (lab, ids) in enumerate(sections):
        mdeg = {}
        for rid in ids:
            for m in set(rp[rid][0] + rp[rid][1]): mdeg[m] = mdeg.get(m, 0) + 1
        hubs = {m for m, d in mdeg.items() if d >= HUB}
        for rid in ids:
            reac, prod = rp[rid]
            kr = [(f'{m}@{rid}' if m in hubs else f'{m}#s{si}') for m in reac]; kp = [(f'{m}@{rid}' if m in hubs else f'{m}#s{si}') for m in prod]
            for k, m in zip(kr + kp, reac + prod): keymid[k] = m
            G.add_node(rid)
            for k in kr: G.add_edge(k, rid)
            for k in kp: G.add_edge(rid, k)
            rp[rid] = (kr, kp)
    is_met = lambda v: v in keymid
    placed, dup_nodes, redirects, rects, texts = {}, [], set(), [], []
    X0, cy = 260, 260
    for si, (lab, ids) in enumerate(sections):
        if not ids: continue
        sub = G.subgraph([v for v in G.nodes if (v in rp and sec_of.get(v) == si) or (v in keymid and v.endswith(f'#s{si}')) or (v in keymid and '@' in v and sec_of.get(v.split('@')[1]) == si)]).copy()
        comps = sorted((sub.subgraph(c).copy() for c in nx.weakly_connected_components(sub)), key=lambda g: -g.number_of_nodes())
        boxes = []
        for Hc in comps:
            layer, x = layered(Hc); layer, x, dups, redirect = fold(Hc, layer, x, is_met)
            boxes.append((Hc, layer, x, dups, max(x.values()) + DX + 420, (max(layer.values()) + 1) * DY)); redirects |= redirect
        top = cy; texts.append((lab, X0 - 120, cy - 60)); cy += 60
        rowh, cx = 0, X0
        for Hc, layer, x, dups, w, h in boxes:
            if cx > X0 and cx + w > maxw: cx = X0; cy += rowh + 160; rowh = 0
            for v in Hc.nodes: placed[v] = (cx + x[v], cy + layer[v] * DY)
            for v, xd, ld in dups: dup_nodes.append((v, cx + xd, cy + ld * DY))
            cx += w + 200; rowh = max(rowh, h)
        cy += rowh + 120
        rects.append(dict(x=X0 - 160, y=top - 110, w=max(maxw, 1000) + 400, h=cy - top + 40, fill=['#f6f8fb', '#f7f5ef'][si % 2], stroke='#c9ccd3'))
        cy += 140
    W = int(max(p[0] for p in placed.values()) + 700); H = int(cy + 200)
    for r in rects: r['w'] = W - r['x'] - 40
    F = escher_map_factory(map_id, title, W, H, frx, fmt, COMP)
    met_data = {}
    for k, (px, py) in placed.items():
        if not is_met(k): continue
        mid = keymid[k]; m = mets.get(mid, {'name': fmt[mid]['name'], 'comp': 'C_ex'})
        F['met'](k, mid, px, py, 'right', name=f"{short(m['name'])} [{COMP.get(m['comp'], m['comp'])}]")
        met_data[mid] = met_value(m) if mid in mets else 0
    for v, px, py in dup_nodes:
        mid = keymid[v]; m = mets[mid]
        F['met'](v + '#dup', mid, px, py, 'right', name=f"{short(m['name'])} [{COMP.get(m['comp'], m['comp'])}]  (continues in next column)")
    for r, v in redirects:
        kr, kp = rp[r]; rp[r] = (kr, [v + '#dup' if q == v else q for q in kp])
    for rid, (kr, kp) in rp.items():
        if rid not in placed: continue
        rxp, ryp = placed[rid]
        F['rxn'](rid, kr, kp, markers=((rxp, ryp - 26), (rxp, ryp), (rxp, ryp + 26)), label='right', sec_side=1, sec_dist=60,
                 sec=len(frx[rid]['mets']) <= MEGA, name=(short(frx[rid]['name'], 44) if rid.startswith('TEST') else f"{rid} {short(frx[rid]['name'], 16)}"))
    for t, x, y in texts: F['text'](t, x, y, 'lane')
    for j, t in enumerate(legend): F['text'](t, 40, 40 + 30 * j, 'note')
    return F['to_json'](), F['deco'](rects), met_data, W, H, placed
