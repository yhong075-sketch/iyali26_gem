"""Escher map of the R2171 inositol leak, three Task A scenarios; no LPs (fluxes come from figure_fluxes.json and
source_check_full.json). Usage: python build_leak_map.py MODEL.xml ESCHER_DIR
F1 WT; F2 INO1 knockout, no inositol uptake (grows at WT); F3 the same, maximum inositol demand (pFBA route at
99.9 % of the maximum). Values are flux x 1000 (umol gDW-1 h-1). Renders SVG with render_jobs.js, PNG with Chrome.
"""
import json
import os
import re
import subprocess
import sys
import warnings

import cobra

sys.path.insert(0, os.path.dirname(__file__))
from escher_kernel import escher_map_factory  # noqa: E402

model_path, E = sys.argv[1:3]
warnings.simplefilter("ignore")
model = cobra.io.read_sbml_model(model_path)
FIG = json.load(open(f"{E}/figure_fluxes.json"))
TRACE = json.load(open(f"{E}/../source_check_full.json"))["traces"]["C1|base|inositol|m814[C_cy]"]
SCEN = {"F1": ("F1  WT, Task A + PO1f", FIG["F1"]["fluxes"], FIG["F1"]["growth"]),
        "F2": ("F2  INO1 knockout, inositol uptake closed: still grows at WT", FIG["F2"]["fluxes"], FIG["F2"]["growth"]),
        "F3": ("F3  as F2, inositol demand maximised (511 mmol): the full loop", TRACE["route"], None)}
COMP = {"C_cy": "cy", "C_er": "er", "C_em": "em", "C_lp": "lp", "C_en": "en", "C_ex": "ex"}
NODES = {
    "m1393": ("m1393[C_cy]", 500, 800, "left", "PI [cy]"), "m1648": ("m1648[C_em]", 500, 1600, "left", "PI [em]"),
    "m379": ("m379[C_em]", 500, 2350, "left", "CDP-DAG [em]"), "m378": ("m378[C_em]", 500, 2900, "left", "PA [em]"),
    "m1675": ("m1675[C_en]", 500, 3450, "left", "PA [en]"),
    "m1994": ("m1994[C_cy]", 1400, 800, "above", "DAG [cy]"), "m2011": ("m2011[C_er]", 2200, 800, "below", "DAG [er]"),
    "m2014": ("m2014[C_er]", 3000, 800, "right", "TAG [er]"),
    "m567": ("m567[C_em]", 1600, 350, "above", "acyl-CoA pool [em]"), "m1984": ("m1984[C_er]", 2200, 350, "above", "acyl-CoA [er]"),
    "m1399": ("m1399[C_cy]", 3000, 1350, "right", "TAG [cy]"), "m1634": ("m1634[C_em]", 3000, 1900, "right", "TAG [em]"),
    "m1640": ("m1640[C_lp]", 3000, 2450, "right", "TAG [lp]"), "m1642": ("m1642[C_lp]", 3000, 2900, "right", "DAG [lp]"),
    "m1670": ("m1670[C_cy]", 2400, 2900, "above", "DAG [cy]"), "m1649": ("m1649[C_em]", 1800, 2900, "above", "DAG [em]"),
    "m1632": ("m1632[C_em]", 1100, 2900, "above", "PC [em]"), "m1659": ("m1659[C_en]", 1100, 3450, "right", "PC [en]"),
    "m1186": ("m1186[C_en]", 1100, 3950, "below", "choline [en]"), "m404": ("m404[C_cy]", 1800, 3950, "below", "choline [cy]"),
    "m405": ("m405[C_cy]", 2400, 3950, "below", None), "m407": ("m407[C_cy]", 2400, 3450, "right", None),
    "m1603": ("m1603[C_em]", 1800, 3450, "below", None),
    "m813": ("m813[C_cy]", 1400, 1350, "left", "inositol 1-phosphate [cy]"), "m814": ("m814[C_cy]", 1400, 1900, "below", "myo-inositol [cy]"),
    "m1605": ("m1605[C_em]", 900, 1900, "left", "myo-inositol [em]"), "m296": ("m296[C_cy]", 2200, 1350, "right", "glucose 6-phosphate [cy]"),
    "m1060": ("m1060[C_ex]", 2200, 1900, "right", "myo-inositol [ex]"),
}
RXNS = {
    "R2171": (["m1393"], ["m1994", "m813"], dict(markers=((800, 800), (880, 800), (960, 800)), label_xy=(545, 965), sec_side=-1)),
    "R1777": (["m1393"], ["m1648"], dict(sec=False, label_xy=(540, 1200))),
    "R2135": (["m567"], ["m1984"], dict(sec=False, label="above")),
    "R2169": (["m2011"], ["m1994"], dict(sec=False, label="above")),
    "R2168": (["m2011", "m1984"], ["m2014"], dict(markers=((2450, 800), (2525, 800), (2600, 800)), label_xy=(2480, 560))),
    "R2170": (["m1399"], ["m2014"], dict(sec=False, label="right")),
    "R1707": (["m1634"], ["m1399"], dict(sec=False, label="right")),
    "R1785": (["m1634"], ["m1640"], dict(sec=False, label="right")),
    "R1771": (["m1642"], ["m1640"], dict(label="right", sec_side=-1)),
    "R1778": (["m1670"], ["m1642"], dict(sec=False, label="below")),
    "R1776": (["m1670"], ["m1649"], dict(sec=False, label="below")),
    "R1726": (["m1649", "m1603"], ["m1632"], dict(markers=((1550, 2900), (1475, 2900), (1400, 2900)), label_xy=(1130, 3110))),
    "R1806": (["m1632"], ["m1659"], dict(sec=False, label="right")),
    "R1842": (["m1659"], ["m1675", "m1186"], dict(sec=False, markers=((900, 3450), (850, 3450), (800, 3450)), label_xy=(560, 3390))),
    "R1809": (["m1675"], ["m378"], dict(sec=False, label="right")),
    "R196": (["m378"], ["m379"], dict(label_xy=(90, 2640), sec_side=1)),
    "R1839": (["m379", "m1605"], ["m1648"], dict(markers=((500, 2150), (500, 2075), (500, 2000)), label_xy=(90, 2060), sec_side=1)),
    "R1607": (["m814"], ["m1605"], dict(sec=False, label="above")),
    "R558": (["m813"], ["m814"], dict(label="right", sec_side=1)),
    "R559": (["m296"], ["m813"], dict(label="above", sec=False)),
    "R862": (["m1060"], ["m814"], dict(label="below", sec=False)),
    "R1042": (["m1186"], ["m404"], dict(sec=False, label="above")),
    "R212": (["m404"], ["m405"], dict(label="above", sec=False)),
    "R213": (["m405"], ["m407"], dict(label="right", sec=False)),
    "R1603": (["m407"], ["m1603"], dict(sec=False, label="above")),
}
SHORT = {"R2171": "PI phospholipase C, written x100", "R559": "INO1 (Ins1P synthase)", "R558": "inositol-1-phosphatase",
         "R1839": "PI synthase", "R862": "inositol uptake", "R1777": "PI transport", "R2168": "DGAT [er]", "R1771": "DGAT [lp]",
         "R1726": "cholinephosphotransferase", "R1842": "phospholipase D", "R196": "CDP-DAG synthase"}
BOXES = [(300, 200, 3400, 1080, "DAG released by R2171 enters TAG in the ER", "#f6f8fb"),
         (700, 1140, 2700, 2180, "Inositol: INO1 or uptake normally; R2171 gives 100 Ins1P per PI", "#fff6ee"),
         (300, 2700, 3400, 4150, "DAG back to PA: Kennedy R1726 + phospholipase D (used in F3)", "#f6f8fb")]
W, H = 3700, 4250
mets = {}
rx = {}
for rid in RXNS:
    r = model.reactions.get_by_id(rid)
    rx[rid] = dict(name=r.name, lb=r.lower_bound, ub=r.upper_bound, gpr=r.gene_reaction_rule, genes=[g.id for g in r.genes],
                   subsystem="", mets={m.id: c for m, c in r.metabolites.items()})
    for m in r.metabolites:
        mets[m.id] = dict(name=m.name, comp=m.compartment)
SCALE = [{"type": "value", "value": 0, "color": "#cfcfcf", "size": 4}, {"type": "value", "value": 1e-5, "color": "#cfcfcf", "size": 4},
         {"type": "value", "value": 1e-4, "color": "#c6dbef", "size": 8}, {"type": "value", "value": 1e-2, "color": "#6baed6", "size": 12},
         {"type": "value", "value": 1, "color": "#2171b5", "size": 16}, {"type": "value", "value": 100, "color": "#08519c", "size": 22},
         {"type": "value", "value": 1e5, "color": "#08306b", "size": 30}]
jobs = []
for key, (title, flux, growth) in SCEN.items():
    F = escher_map_factory(f"r2171_{key}", title, W, H, rx, mets, COMP)
    for k, (sid, x, y, side, lab) in NODES.items():
        F["met"](k, sid, x, y, side, name=lab)
    for rid, (reac, prod, opt) in RXNS.items():
        nm = f"{rid} {SHORT.get(rid) or re.sub(r'_.*$', '', model.reactions.get_by_id(rid).name)[:34]}"
        if key in ("F2", "F3") and rid == "R559":
            nm += "  KNOCKED OUT"
        if key in ("F2", "F3") and rid == "R862":
            nm += "  (no inositol in medium)"
        F["rxn"](rid, reac, prod, name=nm, **opt)
    for j, t in enumerate(["Value: pFBA flux x 1000 (umol gDW-1 h-1); negative = against the drawn direction. Grey = 0. "
                           "Orange outline = R2171 (PI + 100 H2O -> DAG + 100 Ins1P + 100 H+).",
                           "Base 9d655c2b, Task A + PO1f, Gurobi 1e-9. pFBA shows one optimum. Repeated species are the same pool."]):
        F["text"](t, 120, 40 + 34 * j, "note")
    for x0, y0, x1, y1, lab, _ in BOXES:
        F["text"](lab, {1140: 1500, 2700: 1150}.get(y0, x0 + 24), y0 + 40, "lane")
    rects = [dict(x=x0, y=y0, w=x1 - x0, h=y1 - y0, fill=f, stroke="#c9ccd3") for x0, y0, x1, y1, _, f in BOXES]
    rr = [v for v in F["state"]["reactions"].values() if v["bigg_id"] == "R2171"][0]
    rects.append(dict(x=rr["label_x"] - 18, y=rr["label_y"] - 34, w=len(rr["name"]) * 11 + 120, h=80, fill="#fff1dd", stroke="#e07b00"))
    json.dump(F["to_json"](), open(f"{E}/maps/{key}.json", "w"))
    json.dump(F["deco"](rects), open(f"{E}/maps/{key}_deco.json", "w"))
    data = {rid: round(1000 * flux.get(rid, 0.0), 9) for rid in RXNS}
    json.dump(data, open(f"{E}/data/flux_{key}.json", "w"), indent=0)
    full = title + ("" if growth is None else f"   growth {growth:.7f} h-1")
    jobs.append(dict(key=key, map=f"{E}/maps/{key}.json", deco=f"{E}/maps/{key}_deco.json", out=f"{E}/renders/{key}.svg", title=full,
                     W=W, H=H, options=dict(reaction_data=data, reaction_scale=SCALE, reaction_styles=["color", "size", "text", "abs"],
                                            metabolite_no_data_color="#9a9a9a", metabolite_no_data_size=12)))
json.dump(jobs, open(f"{E}/esc_jobs.json", "w"))
p = subprocess.run(["node", "render_jobs.js", f"{E}/esc_jobs.json"], cwd=f"{E}/esc", capture_output=True, text=True, timeout=600)
for r in json.loads(p.stdout.strip().splitlines()[-1]):
    print(r["out"].split("/")[-1], r.get("ok"), r.get("nReactions"), r.get("err", "")[:200])
CH = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
for key in SCEN:
    svg = f"{E}/renders/{key}.svg"
    sc = min(1.0, 3000 / W); w, h = int(W * sc), int((H + 100) * sc)
    tmp = svg[:-4] + "_png.svg"
    open(tmp, "w").write(re.sub(r'width="[^"]+" height="[^"]+"', f'width="{w}" height="{h}"', open(svg).read(), count=1))
    subprocess.run([CH, "--headless=new", "--disable-gpu", "--hide-scrollbars", f"--window-size={w},{h}", f"--screenshot={svg[:-4]}.png",
                    "file://" + tmp], capture_output=True, timeout=180)
    os.remove(tmp)
    print(key, "png", os.path.getsize(svg[:-4] + ".png") // 1024, "KB")
