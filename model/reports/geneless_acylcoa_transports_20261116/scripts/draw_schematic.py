"""Hand-laid schematic of the gene-less CoA-ester transports around the mitochondrion (no LP: reads step outputs).

Style follows the R1172 round's draw_schematic.py: cell outline, mitochondrion as a green body, metabolites as plain
labels, one rounded box per reaction (reaction ID and name, gene rule, flux). Green box = gene rule, pink = no gene,
black outline = moved or merged in the option (in memory only), dashed dark red with a cross = closed in the option.
Usage: python draw_schematic.py OUTDIR
Scenarios (one SVG and PNG each):
  status       largest |flux| each reaction can carry, Task A, no growth floor: as is / option (B0/ALL)
  asis_route   7ccbb568 as is, one propionyl-CoA forced into the matrix (R2113 = -1), minimum total flux
  option_route option B0/ALL, R25 = 1 and R2024 = -1 forced, minimum total flux
"""
import json
import os
import subprocess
import sys
from html import escape

from common import ART

outdir = sys.argv[1]
os.makedirs(outdir, exist_ok=True)
S1 = json.load(open(ART / "out/step1_capacity.json"))
S4 = json.load(open(ART / "out/step4_options.json"))
S5 = json.load(open(ART / "out/step5_figure_flux.json"))
ON = 1e-9
W, H = 2500, 1760
CAP = {"asis": S1["A"]["taskA/gurobi"], "opt": S4["variants"]["B0/ALL"]["fva"]}
FLUX = {"asis_route": S5["asis_route"]["flux"], "option_route": S5["opt_route"]["flux"]}
OPTION = ("option_route",)
CLOSED_IN_OPTION = {"R2025", "R2113", "R2118", "R2077", "R2083", "R2084", "R1291",
                    "R2117", "R2116", "R2115", "R2114", "R2112"}
MOVED = {"R2117", "R2116", "R2115", "R2114", "R2112"}
TITLES = {
    "status": "Largest flux each reaction can carry (Task A + PO1f, no growth floor): 7ccbb568 as is / option (merge, no CoA ester crossing, valine steps in the matrix)",
    "asis_route": "7ccbb568 as is: one propionyl-CoA forced into the matrix. It can only pass through: in on R2113, out on R2025, then propanoate",
    "option_route": "Option (in memory): valine-derived propionyl-CoA made and used inside the matrix; isoleucine-derived propionyl-CoA leaves as propanoate from the cytosol",
}
MI = (772, 232, 1708, 1268)  # inner matrix rectangle

# c = box centre, t = title, g = gene lines, k = gene/nogene/new, a = arrow points (drawn direction),
# s = sign of the drawn direction against the written direction, only = scenarios where drawn,
# a_opt = arrow used in the option scenario when the reaction is rewired
R = {
    # leucine (matrix) and its cytosolic detour
    "R1968": dict(c=(1010, 345), t="R1968 · BCKDH (leucine)", g=["4-gene AND rule"], k="gene", a=[(1010, 298), (1010, 400)], s=1),
    "R1971": dict(c=(1010, 490), t="R1971 · isovaleryl-CoA DH", g=["YALI0D12573g (no W29 match)"], k="gene", a=[(1010, 440), (1010, 545)], s=1),
    "R1972": dict(c=(1010, 635), t="R1972 · MCC", g=["YALI1B19217g"], k="gene", a=[(1010, 585), (1010, 690)], s=1),
    "R2084": dict(c=(735, 375), t="R2084 · isovaleryl-CoA transport", g=["NO GENE; +4 H+ (cy)"], k="nogene", a=[(905, 432), (575, 432)], s=1),
    "R2082": dict(c=(440, 495), t="R2082 · isovaleryl-CoA DH (cy)", g=["YALI1E15555g, TargetP mTP 0.9998"], k="gene", a=[(440, 445), (440, 548)], s=1),
    "R2083": dict(c=(735, 625), t="R2083 · 3-methylcrotonyl-CoA tr.", g=["NO GENE"], k="nogene", a=[(575, 568), (905, 568)], s=-1),
    # CoA synthesis copy and the CoA carrier
    "R1291": dict(c=(735, 805), t="R1291 · pantetheine 4'-P transport", g=["NO GENE"], k="nogene", a=[(575, 862), (905, 862)], s=1),
    "R1292": dict(c=(1010, 935), t="R1292 · PPAT (matrix copy)", g=["NO GENE"], k="nogene", a=[(1010, 885), (1010, 988)], s=1),
    "R265": dict(c=(1010, 1080), t="R265 · DPCK (matrix copy)", g=["YALI1F13202g (also cytosolic R1092)", "TargetP Other 0.99"], k="gene", a=[(1010, 1030), (1010, 1135)], s=1),
    "R1092": dict(c=(400, 1005), t="R615 + R1092 · PPAT, DPCK (cy)", g=["YALI1D24649g; YALI1F13202g"], k="gene", a=[(400, 885), (400, 1130)], s=1),
    "R834": dict(c=(735, 1215), t="R834 · CoA carrier", g=["YALI1F20788g (listed essential)"], k="gene", a=[(560, 1158), (905, 1158)], s=1),
    # valine (matrix part, as is)
    "R1969": dict(c=(1500, 335), t="R1969 · BCKDH (valine)", g=["4-gene AND rule"], k="gene", a=[(1500, 290), (1500, 385)], s=1),
    "R2118": dict(c=(1860, 335), t="R2118 · isobutyryl-CoA transport", g=["NO GENE; +4 H+ (cy)"], k="nogene", a=[(1600, 405), (2030, 405)], s=1),
    # valine (matrix copies, option only)
    "R2117_mi": dict(c=(1500, 475), t="R2117 moved · isobutyryl-CoA DH", g=["YALI1D19252g, TargetP mTP 0.999"], k="new", a=[(1500, 425), (1500, 523)], s=1, only=("status", "option_route")),
    "R2116_mi": dict(c=(1500, 615), t="R2116 moved · enoyl-CoA hydratase", g=["YALI1B14044g, TargetP Other 0.69"], k="new", a=[(1500, 565), (1500, 663)], s=1, only=("status", "option_route")),
    "R2115_mi": dict(c=(1500, 755), t="R2115 moved · HIBCH", g=["NO GENE"], k="new", a=[(1500, 705), (1500, 803)], s=1, only=("status", "option_route")),
    "R2114_mi": dict(c=(1500, 895), t="R2114 moved · HIBADH", g=["YALI1F03883g, TargetP mTP 0.98"], k="new", a=[(1500, 845), (1500, 943)], s=1, only=("status", "option_route")),
    "R2112_mi": dict(c=(1500, 1035), t="R2112 moved · MMSDH (acylating)", g=["YALI1C02550g, TargetP mTP 0.9995"], k="new", a=[(1500, 985), (1500, 1070), (1400, 1100)], s=1, only=("status", "option_route")),
    "R25": dict(c=(1290, 1183), t="R25 · R95 · R552 methylcitrate", g=["gene rules (R95, R552 listed essential)"], k="gene", a=[(1290, 1128), (1290, 1232)], s=1),
    # valine (cytosolic copy, as written in 7ccbb568)
    "R2117": dict(c=(2150, 475), t="R2117 · isobutyryl-CoA DH (cy)", g=["YALI1D19252g", "FAD (cy) in, FADH2 (mi) out"], k="gene", a=[(2150, 425), (2150, 523)], s=1),
    "R2116": dict(c=(2150, 615), t="R2116 · enoyl-CoA hydratase (cy)", g=["YALI1B14044g"], k="gene", a=[(2150, 565), (2150, 663)], s=1),
    "R2115": dict(c=(2150, 755), t="R2115 · HIBCH (cy)", g=["NO GENE"], k="nogene", a=[(2150, 705), (2150, 803)], s=1),
    "R2114": dict(c=(2150, 895), t="R2114 · HIBADH (cy)", g=["YALI1F03883g"], k="gene", a=[(2150, 845), (2150, 943)], s=1),
    "R2112": dict(c=(2150, 1035), t="R2112 · MMSDH (cy)", g=["YALI1C02550g"], k="gene", a=[(2150, 985), (2150, 1083)], s=1),
    "R2077": dict(c=(1830, 1185), t="R2077 · succinyl-CoA transport", g=["NO GENE; +5 H+ (cy)", "(listed in lipid-unlump handoff)"], k="nogene", a=[(2030, 1250), (1610, 1250)], s=-1),
    # propionyl-CoA crossing and the cytosolic outlets
    "R2113": dict(c=(960, 1345), t="R2113 · propionyl-CoA transport", g=["NO GENE (to m1970)"], k="nogene", a=[(1060, 1440), (1240, 1120)], s=-1),
    "R2025": dict(c=(1620, 1345), t="R2025 · propionyl-CoA transport", g=["NO GENE (to m1875)"], k="nogene", a=[(1340, 1120), (1520, 1440)], s=-1),
    "R2285": dict(c=(560, 1465), t="R2281–R2285 · isoleucine branch (cy)", g=["BCKDH rule in the cytosol; YALI0 tags"], k="gene", a=[(330, 1465), (880, 1465)], s=1),
    "R2024": dict(c=(1880, 1400), c_opt=(1560, 1545), t="R2024 · propanoate:CoA ligase", g=["YALI1F08845g; reverse makes ATP"], k="gene", a=[(1640, 1465), (2050, 1465)], s=-1,
                  a_opt=[(1000, 1490), (1000, 1545), (2080, 1545), (2130, 1495)]),
    "R2022": dict(c=(2320, 1385), t="R2023 + R2022 · export", g=["NO GENE"], k="nogene", a=[(2215, 1465), (2470, 1465)], s=1),
}
METS = [  # x, y, text, bold, only
    (1010, 280, "4-Methyl-2-oxopentanoate", 0, None), (1010, 422, "Isovaleryl-CoA", 0, None),
    (1010, 567, "3-Methylcrotonyl-CoA", 0, None), (1010, 712, "3-Methylglutaconyl-CoA", 0, None),
    (440, 432, "Isovaleryl-CoA (cy)", 0, None), (440, 568, "3-Methylcrotonyl-CoA (cy)", 0, None),
    (440, 862, "Pantetheine 4'-P (cy)", 0, None), (1010, 862, "Pantetheine 4'-P", 0, None),
    (1010, 1010, "Dephospho-CoA", 0, None), (1010, 1157, "CoA", 1, None), (400, 1157, "CoA (cy)", 1, None),
    (1500, 272, "2-Oxoisovalerate", 0, None), (1500, 405, "Isobutyryl-CoA", 1, None),
    (1500, 545, "Methacrylyl-CoA", 0, OPTION + ("status",)), (1500, 685, "3-HIB-CoA", 0, OPTION + ("status",)),
    (1500, 825, "3-Hydroxyisobutyrate", 0, OPTION + ("status",)), (1500, 965, "MMSA", 0, OPTION + ("status",)),
    (1290, 1108, "Propionyl-CoA (m78)", 1, None), (1290, 1255, "Pyruvate + succinate", 0, None),
    (2150, 405, "Isobutyryl-CoA (cy)", 1, None), (2150, 545, "Methacrylyl-CoA (cy)", 0, None),
    (2150, 685, "3-HIB-CoA (cy)", 0, None), (2150, 825, "3-Hydroxyisobutyrate (cy)", 0, None),
    (2150, 965, "MMSA (cy)", 0, None), (2150, 1105, "Propanoyl-CoA (m1970)", 0, None),
    (2150, 1250, "Succinyl-CoA (cy)", 0, None), (1500, 1250, "Succinyl-CoA", 0, None),
    (1000, 1465, "Propanoyl-CoA (m1970)", 1, None), (1520, 1465, "Propionyl-CoA (m1875)", 1, None),
    (2140, 1465, "Propanoate", 0, None),
]
NOTES = [(2150, 1132, "same species as m1970 below", "middle"), (2150, 1278, "from gene-less mutase R2075", "middle")]
BOXFILL = {"gene": ("#cfe0bf", "#5f7f4a"), "nogene": ("#ecc3c6", "#a2585e"), "new": ("#cfe0bf", "#111111")}
COL = {"on": "#1f4e9c", "off": "#a6a6a6", "freed": "#2e8b57", "blocked": "#c4c4c4", "closed": "#8b0000",
       "was": "#7d8aa0"}


def num(v):
    a = abs(v)
    return "0" if a <= ON else (f"{a:.3g}" if a < 0.01 else f"{a:.2f}" if a < 100 else f"{a:.0f}")


CHAIN = {"R2116_mi": "R2117_mi", "R2115_mi": "R2117_mi"}  # their species have no other reaction (by construction)


def capacity(fva, rid):
    rid = CHAIN.get(rid, rid)
    if rid not in fva or "min" not in fva[rid]:
        return None
    lo, hi = fva[rid]["min"], fva[rid]["max"]
    return hi if abs(hi) >= abs(lo) else lo


def closed_in(scen, rid):
    return scen in ("option_route",) and rid in CLOSED_IN_OPTION


def build(scen):
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Helvetica, Arial, sans-serif">',
         '<rect width="100%" height="100%" fill="white"/>', "<defs>"]
    for k, c in COL.items():
        o.append(f'<marker id="h_{k}" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="4.2" markerHeight="4.2" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="{c}"/></marker>')
    o.append("</defs>")
    o.append(f'<text x="30" y="46" font-size="23" font-weight="bold" fill="#111">{escape(TITLES[scen])}</text>')
    o.append('<rect x="60" y="80" width="2380" height="1510" rx="34" fill="#f6f3ec" stroke="#2b2b2b" stroke-width="5"/>')
    o.append('<rect x="760" y="220" width="960" height="1060" rx="170" fill="#bcd7c2" stroke="#7fa98a" stroke-width="3"/>')
    o.append('<rect x="772" y="232" width="936" height="1036" rx="160" fill="#d6e7d9" stroke="#7fa98a" stroke-width="3"/>')
    for (x, y, t, sz, col) in [(1255, 760, "Mitochondrial matrix", 22, "#2f5d3c"), (300, 130, "Cytosol", 22, "#555"),
                               (2150, 130, "Cytosol", 22, "#555"), (1010, 250, "leucine", 18, "#b3261e"),
                               (1500, 250 - 0, "", 18, "#b3261e"), (2150, 375, "valine (as written)", 18, "#b3261e"),
                               (400, 820, "CoA synthesis", 18, "#b3261e")]:
        if t:
            o.append(f'<text x="{x}" y="{y}" font-size="{sz}" font-weight="bold" text-anchor="middle" fill="{col}">{t}</text>')
    flux = FLUX.get(scen, {})

    def fval(rid):
        if rid == "R25":
            return flux.get("R25", 0.0)
        if rid == "R1092":
            return flux.get("R1092", 0.0)
        return flux.get(rid, 0.0)

    def style(rid):
        if closed_in(scen, rid):
            return "closed", 3
        if scen == "status":
            a, p = capacity(CAP["asis"], rid), capacity(CAP["opt"], rid)
            if rid in MOVED:
                return "was", 2.5
            if (a in (None, 0) or abs(a) <= ON) and p not in (None,) and abs(p) > ON:
                return "freed", 4.5
            if p is None and a is not None and abs(a) > ON:
                return "was", 3
            if (a is None or abs(a) <= ON) and (p is None or abs(p) <= ON):
                return "blocked", 2.5
            return "was", 3
        v = fval(rid)
        return ("on", 6) if abs(v) > ON else ("off", 2.5)

    def status_line(rid):
        if scen == "status":
            if rid in MOVED:
                return "moved to the matrix in the option", "#444"
            a, p = capacity(CAP["asis"], rid), capacity(CAP["opt"], rid)
            if a is None and p is None:
                return "not computed", "#777"
            first = "—" if rid.endswith("_mi") else ("n.c." if a is None else num(a))
            txt = first + " / " + ("n.c." if p is None else num(p))
            extra = " (closed)" if rid in CLOSED_IN_OPTION and rid not in MOVED else ""
            return f"max {txt}{extra}", "#222"
        if closed_in(scen, rid):
            return ("MOVED to the matrix" if rid in MOVED else "CLOSED (0, 0)"), COL["closed"]
        v = fval(rid)
        if abs(v) <= ON:
            return "flux 0", "#777"
        return f"flux {num(v)}" + (", AGAINST the arrow" if v * R[rid]["s"] < 0 else ""), COL["on"]

    drawn = [rid for rid, r in R.items() if scen in r.get("only", (scen,))]
    for rid in drawn:
        r = R[rid]
        k, w = style(rid)
        pts = r["a_opt"] if (scen in OPTION and "a_opt" in r) else r["a"]
        d = "M " + " L ".join(f"{x} {y}" for x, y in pts)
        dash = ' stroke-dasharray="9 7"' if k == "closed" else ""
        o.append(f'<path d="{d}" fill="none" stroke="{COL[k]}" stroke-width="{w}"{dash} marker-end="url(#h_{k})"/>')
    for x, y, t, anchor in NOTES:
        o.append(f'<text x="{x}" y="{y}" font-size="15" font-style="italic" text-anchor="{anchor}" fill="#444">{escape(t)}</text>')
    for x, y, t, bold, only in METS:
        if only and scen not in only:
            continue
        label = t
        if scen in OPTION and "m1875" in t:
            label = "Propionyl-CoA (m1875): merged into m1970"
        wt = len(label) * 10.0 + 16
        inside = MI[0] < x < MI[2] and MI[1] < y < MI[3]
        o.append(f'<rect x="{x - wt / 2}" y="{y - 21}" width="{wt}" height="30" rx="6" fill="{"#d6e7d9" if inside else "#f6f3ec"}" opacity="0.92"/>')
        o.append(f'<text x="{x}" y="{y}" font-size="19" font-weight="{"bold" if bold else "normal"}" text-anchor="middle" fill="#111">{escape(label)}</text>')
    for rid in drawn:
        r = R[rid]
        txt, col = status_line(rid)
        lines = [(r["t"], 14.5, "bold", "#111")] + [(g, 13, "normal", "#333") for g in r["g"]] + [(txt, 14.5, "bold", col)]
        bw = max(len(s) * (7.9 if wgt == "bold" else 6.9) for s, _, wgt, _ in lines) + 22
        bh = 19 * len(lines) + 10
        cx, cy = r["c_opt"] if (scen in OPTION and "c_opt" in r) else r["c"]
        fill, stroke = BOXFILL[r["k"]]
        sw = 3.5 if r["k"] == "new" else 2
        dash = ""
        if closed_in(scen, rid):
            stroke, sw, dash = COL["closed"], 3, ' stroke-dasharray="8 5"'
        o.append(f'<rect x="{cx - bw / 2}" y="{cy - bh / 2}" width="{bw}" height="{bh}" rx="11" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{dash}/>')
        for i, (s, sz, wgt, c) in enumerate(lines):
            o.append(f'<text x="{cx}" y="{cy - bh / 2 + 21 + 19 * i}" font-size="{sz}" font-weight="{wgt}" text-anchor="middle" fill="{c}">{escape(s)}</text>')
        if closed_in(scen, rid):
            o.append(f'<path d="M {cx + bw / 2 + 6} {cy - bh / 2 - 4} l 20 20 m 0 -20 l -20 20" stroke="{COL["closed"]}" stroke-width="4"/>')
    lx, ly = 80, 1630
    items = [("#cfe0bf", "#5f7f4a", "gene rule in the model"), ("#ecc3c6", "#a2585e", "no gene"),
             ("#cfe0bf", "#111111", "moved copy (option, in memory)"), ("#ecc3c6", COL["closed"], "closed in the option")]
    for i, (f, s, t) in enumerate(items):
        x = lx + i * 470
        dash = ' stroke-dasharray="8 5"' if "closed" in t else ""
        o.append(f'<rect x="{x}" y="{ly}" width="40" height="24" rx="6" fill="{f}" stroke="{s}" stroke-width="2.5"{dash}/>')
        o.append(f'<text x="{x + 52}" y="{ly + 18}" font-size="17" fill="#222">{t}</text>')
    if scen == "status":
        leg = ("Numbers: largest |flux| without a growth floor, as is / option (n.c. = not computed; — = absent; moved R2115/R2116 = chain value). Arrows: green = "
               "can carry flux only in the option; grey = blocked in both; blue-grey = can carry flux as is.")
    else:
        leg = ("Arrows: blue = carries flux in this minimum-total-flux solution (xMAINTENANCE lower bound 0); grey = 0. "
               "TargetP-2.0 values are predictions (mTP = mitochondrial transit peptide).")
    o.append(f'<text x="{lx}" y="{ly + 58}" font-size="17" fill="#222">{escape(leg)}</text>')
    o.append(f'<text x="{lx}" y="{ly + 84}" font-size="17" fill="#222">{escape("Not drawn: TCA cycle, the methylmalonyl route R2074-R2076 (blocked in all variants), R2140-R2142, isoleucine import and transamination, FAD and NAD recycling.")}</text>')
    o.append("</svg>")
    return "\n".join(o)


chrome = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
for scen in ["status", "asis_route", "option_route"]:
    svg = os.path.join(outdir, f"acylcoa_{scen}.svg")
    png = os.path.join(outdir, f"acylcoa_{scen}.png")
    with open(svg, "w") as f:
        f.write(build(scen))
    subprocess.run([chrome, "--headless=new", "--disable-gpu", f"--screenshot={png}", f"--window-size={W},{H}",
                    "--hide-scrollbars", "file://" + os.path.abspath(svg)], check=True, capture_output=True)
    print("wrote", png)
