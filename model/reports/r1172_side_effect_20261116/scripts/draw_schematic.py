"""Hand-laid schematic of the matrix CoA balance around R1172 (no LP: reads the step outputs).

Style follows the methylcitrate rounds' draw_schematic.py (cell outline, mitochondrion as a green body, metabolites
as plain text, one rounded box per reaction with reaction ID and enzyme, gene rule and flux; green box = gene rule,
pink box = no gene, black outline = new reaction (option only), dashed dark-red outline and a cross = closed).
Usage: python draw_schematic.py OUTDIR
Scenarios (one SVG and PNG each):
  status        largest flux each reaction can carry (no growth floor): before / after / option
  before_route  before (R1172 reopened in memory): methylcitrate cycle forced (route, step 2 part C)
  before_wt     before: pFBA at the growth optimum
  after_wt      after (7ccbb568 as is): pFBA at the growth optimum
  option_route  option (O1 + R2075 closed, in memory): methylcitrate cycle forced (route, step 5b)
"""
import json
import os
import subprocess
import sys
from html import escape

from common import ART

outdir = sys.argv[1]
os.makedirs(outdir, exist_ok=True)
S1 = json.load(open(ART / "out/step1_gurobi.json"))["runs"]
S2 = json.load(open(ART / "out/step2_gurobi.json"))
S5B = json.load(open(ART / "out/step5b_routes.json"))
S6 = json.load(open(ART / "out/step6_figure_flux.json"))
ON = 1e-9
W, H = 2320, 1560


def capacity(fva, rid):
    if rid not in fva:
        return None
    lo, hi = fva[rid]["min"], fva[rid]["max"]
    lo = 0.0 if abs(lo) <= ON else lo
    hi = 0.0 if abs(hi) <= ON else hi
    return hi if abs(hi) >= abs(lo) else lo


CAP = {"before": S1["taskA/r1172_open"]["fva"], "after": S1["taskA/as_is"]["fva"],
       "option": S5B["O1_R2075closed"]["fva"]}
FLUX = {"before_route": S2["C"]["taskA/r1172_open"]["R25"]["nonzero"],
        "before_wt": S6["pfba"]["before"]["flux"], "after_wt": S6["pfba"]["after"]["flux"],
        "option_route": S5B["O1_R2075closed"]["routes"]["R25"]["nonzero"]}
TITLES = {
    "status": "Largest flux each reaction can carry (Task A, no growth floor): before R1172 closed / after (7ccbb568) / option O1 + R2075 closed",
    "before_route": "Before (R1172 reopened in memory): the methylcitrate cycle forced. Its CoA enters as propionyl-CoA and leaves as HMG-CoA",
    "before_wt": "Before (R1172 reopened in memory): growth optimum. R834 only returns the CoA that leaves with HMG-CoA",
    "after_wt": "After (7ccbb568, R1172 closed): growth optimum. No reaction can take CoA out of the matrix; the nine stay at 0",
    "option_route": "Option O1 + R2075 closed (in memory): the methylcitrate cycle forced. Propionyl-CoA is made inside the matrix",
}
NINE = ["R25", "R552", "R95", "R2074", "R2075", "R2076", "R2077", "R2112", "R834"]

# reaction: c = box centre, t = title, g = gene lines, k = kind, a = arrow (points in the drawn route direction or
# an SVG path), s = sign of the drawn direction against the written direction, only = scenarios where it is drawn
R = {
    "R963": dict(c=(1265, 160), t="R963 · 2-oxobutanoate transport", g=["NO GENE"], k="nogene", a=[(1030, 120), (1030, 288)], s=+1),
    "R515": dict(c=(780, 375), t="R515 · threonine deaminase", g=["YALI1D03047g"], k="gene", a=[(770, 310), (935, 310)], s=+1),
    "T_OBDHm": dict(c=(1265, 400), t="NEW · 2-oxobutanoate dehydrogenase", g=["rule of R693 (pyruvate dehydrogenase),", "INFERRED; option O1 only"], k="new", a=[(1030, 330), (1030, 478)], s=+1, only=("status", "option_route")),
    "R25": dict(c=(1085, 618), t="R25 · 2-methylcitrate synthase", g=["YALI1E01019g or YALI1E03300g"], k="gene", a="M 985 520 Q 905 560 850 618", s=+1),
    "R95": dict(c=(735, 745), t="R95 · dehydratase", g=["YALI1F03803g (listed essential)"], k="gene", a="M 800 662 Q 770 745 800 822", s=+1),
    "R552": dict(c=(1085, 790), t="R552 · methylisocitrate lyase", g=["YALI1C24124g or", "YALI1F39620g (listed essential)"], k="gene", a="M 905 850 Q 1080 900 1240 850", s=+1),
    "R693": dict(c=(1170, 985), t="R693 · pyruvate dehydrogenase", g=["5-gene AND rule"], k="gene", a=[(1030, 945), (1030, 1035)], s=+1),
    "R83": dict(c=(845, 1110), t="R83 · thiolase", g=["YALI1E13899g"], k="gene", a=[(960, 1058), (870, 1058)], s=+1),
    "R412": dict(c=(910, 1215), t="R412 · HMG-CoA synthase", g=["YALI1F38116g"], k="gene", a=[(760, 1078), (760, 1178)], s=+1),
    "R1172": dict(c=(420, 1125), t="R1172 · HMG-CoA transport", g=["NO GENE; closed 2026-11-12"], k="nogene", a=[(690, 1200), (395, 1200)], s=-1),
    "R410": dict(c=(470, 1290), t="R410 · HMG-CoA reductase", g=["YALI1E05731g"], k="gene", a=[(300, 1222), (300, 1355)], s=+1),
    "R834": dict(c=(420, 940), t="R834 · CoA carrier", g=["YALI1F20788g (listed essential)"], k="gene", a=[(330, 1010), (650, 1010)], s=+1),
    "R2113": dict(c=(1600, 430), t="R2113 · propionyl-CoA transport", g=["NO GENE"], k="nogene", a=[(1680, 500), (1100, 500)], s=-1),
    "R2285": dict(c=(1990, 340), t="R2285 · 2-methylacetoacetyl-CoA thiolase", g=["YALI0E18568g; isoleucine breakdown"], k="gene", a=[(1760, 320), (1760, 478)], s=+1),
    "R2112": dict(c=(2045, 590), t="R2112 · semialdehyde dehydrogenase", g=["YALI1C02550g"], k="gene", a=[(1760, 522), (1760, 640)], s=-1),
    "R2074": dict(c=(2045, 725), t="R2074 · semialdehyde dehydrogenase", g=["4 genes (OR)"], k="gene", a=[(1760, 678), (1760, 778)], s=+1),
    "R2076": dict(c=(2045, 860), t="R2076 · methylmalonate + CoA", g=["NO GENE; no ATP"], k="nogene", a=[(1760, 815), (1760, 915)], s=+1),
    "R2075": dict(c=(2045, 995), t="R2075 · methylmalonyl-CoA mutase", g=["NO GENE; no homologue (KEGG)"], k="nogene", a=[(1760, 952), (1760, 1052)], s=+1),
    "R2077": dict(c=(1620, 1215), t="R2077 · succinyl-CoA transport", g=["NO GENE"], k="nogene", a=[(1760, 1092), (1760, 1150), (1450, 1150)], s=-1),
    "R741": dict(c=(1335, 1290), t="R741 · succinyl-CoA ligase", g=["2-gene AND rule"], k="gene", a=[(1405, 1128), (1345, 885)], s=+1),
    "R2025": dict(c=(440, 455), t="R2025 · propionyl-CoA transport", g=["NO GENE"], k="nogene", a=[(960, 520), (380, 520)], s=-1),
    "R2024": dict(c=(430, 640), t="R2024 · propanoate:CoA ligase", g=["YALI1F08845g"], k="gene", a=[(290, 542), (290, 690)], s=-1),
}
METS = [  # x, y, text, bold
    (1030, 110, "2-Oxobutanoate", 0), (690, 315, "L-Threonine", 0), (1030, 315, "2-Oxobutanoate", 0),
    (1030, 505, "Propionyl-CoA", 1), (1300, 595, "Oxaloacetate", 0), (810, 645, "2-Methylcitrate", 0),
    (810, 845, "2-Methylisocitrate", 0), (1300, 865, "Succinate", 0), (1030, 935, "Pyruvate", 0),
    (1030, 1062, "Acetyl-CoA", 0), (760, 1062, "Acetoacetyl-CoA", 0), (760, 1205, "HMG-CoA", 1),
    (300, 1205, "HMG-CoA", 1), (300, 1380, "Mevalonate + CoA", 0), (275, 1015, "CoA", 1), (700, 1015, "CoA", 1),
    (1760, 300, "Isoleucine breakdown", 0), (1760, 505, "Propanoyl-CoA (m1970)", 1), (1760, 660, "Methylmalonate semialdehyde", 0),
    (1760, 797, "Methylmalonate", 0), (1760, 935, "Methylmalonyl-CoA", 0), (1760, 1072, "Succinyl-CoA", 0),
    (1380, 1150, "Succinyl-CoA", 0), (290, 520, "Propionyl-CoA (m1875)", 1), (290, 712, "Propanoate", 0),
]
OBDH_NOTES = [(1015, 400, "CoA + NAD⁺ in", "end"), (1015, 428, "CO₂ + NADH out", "end")]
NOTES = [(1120, 560, "CoA out", "start"),
         (1015, 1000, "CoA in", "end"), (1125, 1100, "2 acetyl-CoA → acetoacetyl-CoA + CoA", "middle"),
         (745, 1168, "+ acetyl-CoA, CoA out", "end"), (285, 1300, "CoA out (cytosol)", "end"),
         (1745, 600, "CO₂ + NADH in, CoA out", "end"), (1780, 862, "CoA in", "start"), (285, 625, "AMP + PPi → ATP", "end")]
BOXFILL = {"gene": ("#cfe0bf", "#5f7f4a"), "nogene": ("#ecc3c6", "#a2585e"), "new": ("#cfe0bf", "#111111")}
COL = {"on": "#1f4e9c", "off": "#a6a6a6", "freed": "#2e8b57", "side": "#d6452a", "blocked": "#c4c4c4",
       "closed": "#8b0000", "was": "#7d8aa0"}


def num(v):
    a = abs(v)
    return "0" if a <= ON else (f"{a:.3g}" if a < 0.01 else f"{a:.2f}" if a < 100 else f"{a:.0f}")


def closed_in(scen, rid):
    if rid == "R1172":
        return scen in ("after_wt", "option_route", "status")
    if rid == "R2075":
        return scen == "option_route"
    return False


def build(scen):
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Helvetica, Arial, sans-serif">',
         '<rect width="100%" height="100%" fill="white"/>', "<defs>"]
    for k, c in COL.items():
        o.append(f'<marker id="h_{k}" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="4.2" markerHeight="4.2" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="{c}"/></marker>')
    o.append("</defs>")
    o.append(f'<text x="30" y="48" font-size="25" font-weight="bold" fill="#111">{escape(TITLES[scen])}</text>')
    o.append('<rect x="150" y="70" width="2140" height="1340" rx="34" fill="#f6f3ec" stroke="#2b2b2b" stroke-width="5"/>')
    o.append('<rect x="548" y="218" width="964" height="1124" rx="200" fill="#bcd7c2" stroke="#7fa98a" stroke-width="3"/>')
    o.append('<rect x="560" y="230" width="940" height="1100" rx="190" fill="#d6e7d9" stroke="#7fa98a" stroke-width="3"/>')
    for (x, y, t, sz, col) in [(1030, 1312, "Mitochondrial matrix", 23, "#2f5d3c"), (300, 130, "Cytosol", 23, "#555"),
                               (1990, 130, "Cytosol", 23, "#555"), (1000, 703, "methylcitrate cycle", 19, "#b3261e")]:
        o.append(f'<text x="{x}" y="{y}" font-size="{sz}" font-weight="bold" text-anchor="middle" fill="{col}">{t}</text>')
    flux = FLUX.get(scen, {})

    def style(rid):
        if closed_in(scen, rid):
            return "closed", 3
        if scen == "status":
            b, a, p = (capacity(CAP[x], rid) for x in ("before", "after", "option"))
            if rid == "T_OBDHm":
                return "freed", 4.5
            if a == 0 and p not in (None, 0):
                return ("freed" if rid in ("R25", "R95", "R552") else "side"), 4.5
            if a == 0:
                return "blocked", 2.5
            return "was", 3
        v = flux.get(rid, 0.0)
        return ("on", 6) if abs(v) > ON else ("off", 2.5)

    def status_line(rid):
        if scen == "status":
            vals = [capacity(CAP[x], rid) for x in ("before", "after", "option")]
            if all(v is None for v in vals):
                return "capacity not computed", "#777"
            txt = " / ".join("n.c." if v is None else num(v) for v in vals)
            note = ""
            if rid == "R834" and vals[2]:
                note = " (side route B)"
            if rid == "R2075":
                note = " (closed in option)"
            return f"max {txt}{note}", (COL["closed"] if closed_in(scen, rid) else "#222")
        v = flux.get(rid, 0.0)
        if closed_in(scen, rid):
            return "CLOSED (0, 0)", COL["closed"]
        if abs(v) <= ON:
            return "flux 0", "#777"
        return f"flux {num(v)}" + (", AGAINST the arrow" if v * R[rid]["s"] < 0 else ""), COL["on"]

    drawn = [rid for rid, r in R.items() if scen in r.get("only", (scen,))]
    for rid in drawn:
        r = R[rid]
        k, w = style(rid)
        d = r["a"] if isinstance(r["a"], str) else "M " + " L ".join(f"{x} {y}" for x, y in r["a"])
        dash = ' stroke-dasharray="9 7"' if k == "closed" else ""
        o.append(f'<path d="{d}" fill="none" stroke="{COL[k]}" stroke-width="{w}"{dash} marker-end="url(#h_{k})"/>')
        if k == "closed":
            (x1, y1), (x2, y2) = r["a"][0], r["a"][-1]
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2
            o.append(f'<path d="M {mx - 16} {my - 16} L {mx + 16} {my + 16} M {mx - 16} {my + 16} L {mx + 16} {my - 16}" stroke="{COL["closed"]}" stroke-width="5"/>')
    for x, y, t, anchor in NOTES + (OBDH_NOTES if "T_OBDHm" in drawn else []):
        o.append(f'<text x="{x}" y="{y}" font-size="15" font-style="italic" text-anchor="{anchor}" fill="#444">{t}</text>')
    for x, y, t, bold in METS:
        wt = len(t) * 10.4 + 16
        inside = 560 < x < 1500 and 230 < y < 1330
        o.append(f'<rect x="{x - wt / 2}" y="{y - 21}" width="{wt}" height="30" rx="6" fill="{"#d6e7d9" if inside else "#f6f3ec"}" opacity="0.92"/>')
        o.append(f'<text x="{x}" y="{y}" font-size="20" font-weight="{"bold" if bold else "normal"}" text-anchor="middle" fill="#111">{escape(t)}</text>')
    for rid in drawn:
        r = R[rid]
        txt, col = status_line(rid)
        genes = ["NO GENE; reopened here in memory"] if (rid == "R1172" and scen.startswith("before")) else r["g"]
        lines = [(r["t"], 15, "bold", "#111")] + [(g, 13.5, "normal", "#333") for g in genes] + [(txt, 15, "bold", col)]
        bw = max(len(s) * (8.2 if wgt == "bold" else 7.1) for s, _, wgt, _ in lines) + 24
        bh = 20 * len(lines) + 12
        cx, cy = r["c"]
        fill, stroke = BOXFILL[r["k"]]
        sw = 3.5 if r["k"] == "new" else 2
        if closed_in(scen, rid):
            stroke, sw = COL["closed"], 3
        dash = ' stroke-dasharray="8 5"' if closed_in(scen, rid) else ""
        o.append(f'<rect x="{cx - bw / 2}" y="{cy - bh / 2}" width="{bw}" height="{bh}" rx="12" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{dash}/>')
        for i, (s, sz, wgt, c) in enumerate(lines):
            o.append(f'<text x="{cx}" y="{cy - bh / 2 + 23 + 20 * i}" font-size="{sz}" font-weight="{wgt}" text-anchor="middle" fill="{c}">{escape(s)}</text>')
    # legend
    lx, ly = 170, 1440
    items = [("#cfe0bf", "#5f7f4a", "gene rule in the model"), ("#ecc3c6", "#a2585e", "no gene"),
             ("#ecc3c6", COL["closed"], "closed (0, 0)")]
    if scen in ("status", "option_route"):
        items.append(("#cfe0bf", "#111111", "new in option O1 (in memory only)"))
    for i, (f, s, t) in enumerate(items):
        x = lx + i * 420
        o.append(f'<rect x="{x}" y="{ly}" width="40" height="24" rx="6" fill="{f}" stroke="{s}" stroke-width="2.5"/>')
        o.append(f'<text x="{x + 52}" y="{ly + 18}" font-size="17" fill="#222">{t}</text>')
    if scen == "status":
        o.append(f'<text x="{lx}" y="{ly + 56}" font-size="17" fill="#222">Arrows: green = freed by O1 inside the matrix; red = freed by O1 only through a gene-less side route; grey = blocked; numbers = maximum |flux| before / after / option (n.c. = not computed).</text>')
    else:
        o.append(f'<text x="{lx}" y="{ly + 56}" font-size="17" fill="#222">Arrows: blue = carries flux in this solution (minimum total flux or pFBA); grey = 0. Arrows point in the drawn route direction; "AGAINST the arrow" marks the opposite sign.</text>')
        o.append(f'<text x="{lx}" y="{ly + 82}" font-size="17" fill="#222">Not drawn: the cytosolic HMG-CoA route R82 + R411 (feeds R410 when R1172 is closed), the TCA cycle and the supply of oxaloacetate, threonine and isoleucine.</text>')
    o.append("</svg>")
    return "\n".join(o)


chrome = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
for scen in ["status", "before_route", "before_wt", "after_wt", "option_route"]:
    svg = os.path.join(outdir, f"r1172_{scen}.svg")
    png = os.path.join(outdir, f"r1172_{scen}.png")
    with open(svg, "w") as f:
        f.write(build(scen))
    subprocess.run([chrome, "--headless=new", "--disable-gpu", f"--screenshot={png}", f"--window-size={W},{H}",
                    "--hide-scrollbars", "file://" + os.path.abspath(svg)], check=True, capture_output=True)
    print("wrote", png)
