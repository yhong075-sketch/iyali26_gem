import os, re, sys, json, math, shutil, subprocess

DEFAULT_SCALE = [{"type": "value", "value": 0, "color": "#c8c8c8", "size": 6},
                 {"type": "value", "value": 0.01, "color": "#9ecae1", "size": 10},
                 {"type": "value", "value": 1, "color": "#3182bd", "size": 18},
                 {"type": "value", "value": 10, "color": "#08306b", "size": 30}]
SHORT_NAMES = {"oxygen": "O₂", "H2O": "H₂O", "carbon dioxide": "CO₂", "H+": "H⁺", "NAD": "NAD⁺", "NADP(+)": "NADP⁺",
               "diphosphate": "PPi", "phosphate": "Pi", "hydrogen peroxide": "H₂O₂", "bicarbonate": "HCO₃⁻", "ammonium": "NH₄⁺"}
CHAR_W = {"rxn": 13.5, "node": 10.5, "text": 10.0}


def cobra_sandbox_import(cache_dir="cobra_cache"):
    """Import cobra with its cache redirected into the workspace (default ~/Library/Caches is not writable in the sandbox)."""
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.abspath(cache_dir)
    for k in [k for k in sys.modules if k.startswith("cobra")]:
        del sys.modules[k]
    for modname in ("appdirs", "platformdirs"):
        try:
            mod = __import__(modname)
            mod.user_cache_dir = lambda *a, **k: path
        except ImportError:
            pass
    import cobra
    return cobra


def run_cobra_script(script, args=None, timeout=900):
    """Run an LP script in a child process with a hard timeout (a hung glpk solve cannot be interrupted inside a kernel)."""
    cmd = [sys.executable, script] + [str(a) for a in (args or [])]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    out = "\n".join(l for l in p.stdout.splitlines() if "warn" not in l.lower())
    if p.returncode != 0:
        raise RuntimeError(f"{script} failed:\n{p.stderr[-3000:]}")
    return out


def scan_model(model, name_regex=None, formulas=None, ec=None, kegg=None):
    """Presence check: metabolites by name regex / formula, reactions by EC or KEGG reaction id."""
    def aslist(v):
        return [v] if isinstance(v, str) else list(v or [])
    mets = []
    for m in model.metabolites:
        if (name_regex and re.search(name_regex, m.name or "", re.I)) or (formulas and m.formula in set(formulas)):
            mets.append(dict(id=m.id, name=m.name, comp=m.compartment, formula=m.formula, n_rxns=len(m.reactions),
                             rxns=[r.id for r in m.reactions]))
    rxns = []
    for r in model.reactions:
        hit_ec = [e for e in aslist(r.annotation.get("ec-code")) if e in set(ec or [])]
        hit_k = [k for k in aslist(r.annotation.get("kegg.reaction")) if k in set(kegg or [])]
        if hit_ec or hit_k:
            rxns.append(dict(id=r.id, name=r.name, ec=hit_ec, kegg=hit_k, subsystem=r.subsystem, gpr=r.gene_reaction_rule,
                             equation=r.build_reaction_string(use_metabolite_names=True), bounds=r.bounds))
    return dict(metabolites=mets, reactions=rxns)


def export_rxn_info(model, rxn_ids, subsystem_override=None):
    """Stoichiometry/GPR/subsystem dicts consumed by escher_map_factory (also JSON-serialisable for handoff)."""
    rx, mt = {}, {}
    for rid in rxn_ids:
        r = model.reactions.get_by_id(rid)
        rx[rid] = dict(name=r.name, lb=r.lower_bound, ub=r.upper_bound, gpr=r.gene_reaction_rule, genes=[g.id for g in r.genes],
                       subsystem=(subsystem_override or {}).get(rid, r.subsystem), mets={m.id: c for m, c in r.metabolites.items()})
        for m in r.metabolites:
            mt[m.id] = dict(name=m.name, comp=m.compartment)
    return rx, mt


def gpr_short(rule, test=False):
    if not rule:
        return "no gene (test reaction)" if test else "no GPR in model"
    g = [x for x in dict.fromkeys(re.findall(r"[^\s()]+", rule)) if x not in ("and", "or")]
    if len(g) <= 2:
        return rule
    if " and " not in rule and len(g) <= 4:
        return " or ".join(g)
    kind = "AND/OR" if (" and " in rule and " or " in rule) else ("AND complex" if " and " in rule else "OR isozymes")
    return f"{len(g)} genes ({kind})"


def escher_map_factory(map_id, title, width, height, rxns, mets, comp_tags=None, short_names=None):
    """Escher schema 1-0-0 map builder (closures). Returns dict: met, rxn, text, to_json, deco.
    met(key, met_id, x, y, side='right'|'left'|'above'|'below'|(dx,dy), name=None)
    rxn(rid, reactant_keys, product_keys, sec_side=-1, label='left'|'right'|'above'|'below', label_xy=None,
        markers=None, sec=True, sec_dist=70, name=None, show_subsystem=True) -> reaction node id
    text(txt, x, y, role='comp'|'lane'|'gpr'|'note');  deco(rects) -> {'rects','label_roles'} for render_jobs.js"""
    ct = comp_tags or {}
    sn = dict(SHORT_NAMES, **(short_names or {}))
    S = dict(nodes={}, reactions={}, labels={}, keys={}, roles={}, k=0, sk=0)

    def nid():
        S["k"] += 1
        return str(S["k"])

    def clean(n):
        n = re.sub(r"_(?:[A-Z][a-z]?\d*)+(?:\(.*\))?.*$", "", n) if "_" in n else n
        return n.replace("_p+1", "").rstrip("_ ").strip()

    def mname(mid, primary=True):
        n = clean(mets[mid]["name"]) or mid
        n = sn.get(n, n)
        return f"{n} [{ct.get(mets[mid]['comp'], mets[mid]['comp'])}]" if primary else n

    def met(key, mid, x, y, side="right", name=None):
        n = name or mname(mid)
        wd = len(n) * CHAR_W["node"]
        if isinstance(side, (tuple, list)):
            lx, ly = x + side[0], y + side[1]
        else:
            lx, ly = {"right": (x + 26, y + 8), "left": (x - 26 - wd, y + 8), "below": (x - wd / 2, y + 44), "above": (x - wd / 2, y - 26)}[side]
        i = nid()
        S["nodes"][i] = dict(node_type="metabolite", x=x, y=y, bigg_id=mid, name=n, label_x=lx, label_y=ly, node_is_primary=True)
        S["keys"][key] = i
        return i

    def mk(kind, x, y):
        i = nid()
        S["nodes"][i] = dict(node_type=kind, x=x, y=y)
        return i

    def text(t, x, y, role="note"):
        i = nid()
        S["labels"][i] = dict(text=t, x=x, y=y)
        S["roles"]["l" + i] = role
        return i

    def rxn(rid, reac, prod, sec_side=-1, label="left", label_xy=None, markers=None, sec=True, sec_dist=70, name=None, show_subsystem=True):
        rd = rxns[rid]
        P = lambda k: (S["nodes"][S["keys"][k]]["x"], S["nodes"][S["keys"][k]]["y"])
        A, B = [P(k) for k in reac], [P(k) for k in prod]
        ax_, ay_ = sum(p[0] for p in A) / len(A), sum(p[1] for p in A) / len(A)
        bx, by = sum(p[0] for p in B) / len(B), sum(p[1] for p in B) / len(B)
        L = math.hypot(bx - ax_, by - ay_)
        dx, dy = (bx - ax_) / L, (by - ay_) / L
        px, py = -dy, dx
        if markers:
            m1, md, m2 = markers
        else:
            f = lambda t: (ax_ + t * (bx - ax_), ay_ + t * (by - ay_))
            m1, md, m2 = f(0.40), f(0.5), f(0.60)
        i1, im, i2 = mk("multimarker", *m1), mk("midmarker", *md), mk("multimarker", *m2)
        segs = {}

        def seg(a, b, b1=None, b2=None):
            S["sk"] += 1
            segs[str(S["sk"])] = dict(from_node_id=a, to_node_id=b, b1=b1, b2=b2)
        prim = [S["nodes"][S["keys"][k]]["bigg_id"] for k in reac + prod]
        for k in reac:
            seg(S["keys"][k], i1)
        seg(i1, im)
        seg(im, i2)
        for k in prod:
            seg(i2, S["keys"][k])
        if sec:
            ins = [(m, c) for m, c in rd["mets"].items() if c < 0 and m not in prim]
            outs = [(m, c) for m, c in rd["mets"].items() if c > 0 and m not in prim]
            for grp, mm, sg in ((ins, m1, -1), (outs, m2, +1)):
                horiz = abs(dx) > abs(dy)
                labs = [(f"{abs(c):g} " if abs(c) != 1 else "") + mname(m, False) for m, c in grp]
                step = (max([len(t) for t in labs] + [3]) * CHAR_W["node"] + 16) if horiz else 30
                for j, ((m, c), nm) in enumerate(zip(grp, labs)):
                    along = (20 + step / 2 if horiz else 20) + step * j
                    sx = mm[0] + sg * dx * along + sec_side * px * sec_dist
                    sy = mm[1] + sg * dy * along + sec_side * py * sec_dist
                    wd = len(nm) * CHAR_W["node"]
                    ox, oy = sec_side * px, sec_side * py
                    if abs(ox) >= abs(oy):
                        lx, ly = (sx + 14, sy + 7) if ox > 0 else (sx - 14 - wd, sy + 7)
                    else:
                        lx, ly = (sx - wd / 2, sy + 32) if oy > 0 else (sx - wd / 2, sy - 14)
                    ni = nid()
                    S["nodes"][ni] = dict(node_type="metabolite", x=sx, y=sy, bigg_id=m, name=nm, label_x=lx, label_y=ly, node_is_primary=False)
                    if sg < 0:
                        seg(ni, i1, dict(x=sx + dx * 20, y=sy + dy * 20), dict(x=mm[0] - dx * 26, y=mm[1] - dy * 26))
                    else:
                        seg(i2, ni, dict(x=mm[0] + dx * 26, y=mm[1] + dy * 26), dict(x=sx - dx * 20, y=sy - dy * 20))
        nm = name or f"{rid} {rd['name']}"
        is_test = str(rd.get("subsystem", "")).startswith("TEST")
        sub = gpr_short(rd["gpr"], is_test) + (f"  ·  {rd['subsystem']}" if show_subsystem and rd.get("subsystem") else "")
        wd = (len(nm) + 8) * CHAR_W["rxn"]
        if label_xy:
            lx, ly = label_xy
        else:
            lx, ly = {"left": (md[0] - 34 - max(wd, len(sub) * CHAR_W["text"]), md[1] + 8), "right": (md[0] + 34, md[1] + 8),
                      "above": (md[0] - wd / 2, md[1] - 70), "below": (md[0] - wd / 2, md[1] + 50)}[label]
        ri = nid()
        S["reactions"][ri] = dict(name=nm, bigg_id=rid, reversibility=rd["lb"] < 0 and rd["ub"] > 0, label_x=lx, label_y=ly,
                                  gene_reaction_rule=rd["gpr"], genes=[dict(bigg_id=g, name=g) for g in rd["genes"]],
                                  metabolites=[dict(bigg_id=m, coefficient=c) for m, c in rd["mets"].items()], segments=segs)
        text(sub, lx, ly + 28, "gpr")
        return ri

    def to_json():
        return [dict(map_name=title, map_id=map_id, map_description=title + " — generated from SBML stoichiometry; flux supplied as reaction data.",
                     homepage="https://escher.github.io", schema="https://escher.github.io/escher/jsonschema/1-0-0#"),
                dict(reactions=S["reactions"], nodes=S["nodes"], text_labels=S["labels"], canvas=dict(x=0, y=0, width=width, height=height))]

    def deco(rects):
        return dict(rects=rects, label_roles=S["roles"])
    return dict(met=met, rxn=rxn, text=text, to_json=to_json, deco=deco, state=S)


def install_escher(workdir="esc"):
    """npm-install escher@1.8.1 + jsdom@24 with a workspace-local cache (needs an env with nodejs; ~/.npm is not writable)."""
    os.makedirs(workdir, exist_ok=True)
    env = dict(os.environ, npm_config_cache=os.path.abspath(os.path.join(workdir, ".npmcache")), HOME=os.path.abspath(workdir))
    npm = shutil.which("npm")
    if not npm:
        raise RuntimeError("npm not on PATH: run in an environment with nodejs (e.g. manage_environments create packages=['nodejs>=18','cairosvg'])")
    if not os.path.exists(os.path.join(workdir, "package.json")):
        subprocess.run([npm, "init", "-y"], cwd=workdir, env=env, capture_output=True, check=True)
    subprocess.run([npm, "install", "--no-audit", "--no-fund", "--loglevel=error", "escher@1.8.1", "jsdom@24"], cwd=workdir, env=env,
                   capture_output=True, text=True, check=True)
    js = os.path.join(workdir, "node_modules", "escher", "dist", "escher.min.js")
    if not os.path.exists(js):
        raise RuntimeError("escher.min.js missing after npm install")
    return js


def render_escher(jobs, workdir="esc", png_width=2600):
    """Render Escher maps headlessly. jobs: [{key, map, deco, out(.svg), title, W, H, options:{reaction_data, reaction_scale, ...}}].
    Paths are relative to workdir. Copies the skill's scripts/render_jobs.js into workdir, runs node, converts SVG->PNG."""
    here = os.path.dirname(sys._getframe().f_code.co_filename)
    src = os.path.join(here, "scripts", "render_jobs.js") if here else ""
    dst = os.path.join(workdir, "render_jobs.js")
    if src and os.path.exists(src):
        shutil.copy(src, dst)
    if not os.path.exists(dst):
        raise RuntimeError("render_jobs.js not found (skill dir unavailable); copy it from the skill's scripts/ folder")
    json.dump(jobs, open(os.path.join(workdir, "jobs.json"), "w"), indent=1, ensure_ascii=False)
    node = shutil.which("node")
    if not node:
        raise RuntimeError("node not on PATH: run in the nodejs environment")
    p = subprocess.run([node, "render_jobs.js", "jobs.json"], cwd=workdir, capture_output=True, text=True, timeout=600)
    res = json.loads(p.stdout.strip().splitlines()[-1])
    try:
        import cairosvg
        for r in res:
            if r.get("ok"):
                svg = os.path.join(workdir, r["out"])
                cairosvg.svg2png(url=svg, write_to=svg[:-4] + ".png", output_width=png_width)
    except ImportError:
        pass
    return res


def build_escher_html(maps, data, scenarios, out, escher_js, scale=None, heading="Escher flux maps"):
    """Self-contained offline HTML with map + scenario drop-downs (patches Escher's publicPath check so inlined JS runs).
    maps: {map_key: (label, map_json)}; scenarios: {scen_key: label};
    data: {map_key: {scen_key: {"data": {rxn_id: flux}, "title": str}}}."""
    js = open(escher_js).read()
    bad = 'if(!t)throw new Error("Automatic publicPath is not supported in this browser");'
    if bad in js:
        js = js.replace(bad, 'if(!t)t="https://unpkg.com/escher@1.8.1/dist/";')
    sc = scale or DEFAULT_SCALE
    M = {k: v[1] for k, v in maps.items()}
    D = {f"{k}_{s}": data[k][s] for k in maps for s in scenarios}
    mo = "".join(f'<option value="{k}">{v[0]}</option>' for k, v in maps.items())
    so = "".join(f'<option value="{k}">{v}</option>' for k, v in scenarios.items())
    html = f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>{heading}</title>
<style>body{{font-family:sans-serif;margin:0}}#bar{{padding:8px 12px;background:#f4f4f4;border-bottom:1px solid #ccc}}#c{{width:100vw;height:calc(100vh - 80px)}}#t{{font-size:13px;margin-top:4px}}</style>
<script>{js}</script></head><body><div id="bar">Map: <select id="m">{mo}</select> &nbsp; Scenario: <select id="s">{so}</select><div id="t"></div></div><div id="c"></div>
<script>const MAPS={json.dumps(M)};const DATA={json.dumps(D)};const SCALE={json.dumps(sc)};
function draw(){{const m=document.getElementById('m').value,s=document.getElementById('s').value,d=DATA[m+'_'+s];
document.getElementById('t').textContent=d.title+'  |  pink = reaction absent in this scenario; grey = zero flux';document.getElementById('c').innerHTML='';
escher.Builder(MAPS[m],null,null,escher.libs.d3_select('#c'),{{menu:'zoom',enable_editing:false,reaction_data:d.data,reaction_scale:SCALE,
reaction_styles:['color','size','text','abs'],reaction_scale_preset:false,reaction_no_data_color:'#e8a3ae',reaction_no_data_size:5,
identifiers_on_map:'name',never_ask_before_quit:true,show_gene_reaction_rules:false,fill_screen:false}});}}
document.getElementById('m').onchange=draw;document.getElementById('s').onchange=draw;draw();</script></body></html>"""
    open(out, "w").write(html)
    return out
