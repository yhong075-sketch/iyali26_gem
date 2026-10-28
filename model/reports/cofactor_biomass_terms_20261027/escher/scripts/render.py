"""Render the cofactor-terms map for each scenario to SVG (jsdom + Escher 1.8.1) and PNG (headless Chrome).

Usage: python render.py WORKDIR [PNG_WIDTH]
WORKDIR holds maps/, data/, index.json (from build.py) and esc/ (node_modules with escher and jsdom); this
folder's render_jobs.js (copied from the gap census) is copied into esc/ because it loads node_modules from
its own folder. Writes WORKDIR/renders/<scenario>.svg and .png, plus compare.svg/.png (before | after side by side).
"""
import json
import os
import re
import shutil
import subprocess
import sys

from palette import MET_SCALE, RXN_SCALES, TITLES

work = os.path.abspath(sys.argv[1])
png_width = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
IDX = json.load(open(f"{work}/index.json"))
here = os.path.dirname(os.path.abspath(__file__))
shutil.copy(f"{here}/render_jobs.js", f"{work}/esc/render_jobs.js")
os.makedirs(f"{work}/renders", exist_ok=True)
jobs = []
for scen in TITLES:
    data = json.load(open(f"{work}/data/{scen}.json"))
    jobs.append(dict(key=scen, map=f"{work}/maps/cofactor_terms.json", deco=f"{work}/maps/cofactor_terms_deco.json",
                     out=f"{work}/renders/{scen}.svg", title=TITLES[scen], W=IDX["W"], H=IDX["H"],
                     options=dict(reaction_data=data, reaction_scale=RXN_SCALES[scen], metabolite_data=IDX["met_data"],
                                  metabolite_scale=MET_SCALE, metabolite_styles=["color", "size"], metabolite_scale_preset=False,
                                  metabolite_no_data_color="#e0e0e0", metabolite_no_data_size=8,
                                  reaction_no_data_color="#d9d9d9", reaction_no_data_size=6)))
json.dump(jobs, open(f"{work}/esc/jobs_cofactor.json", "w"))
p = subprocess.run(["node", "render_jobs.js", "jobs_cofactor.json"], cwd=f"{work}/esc", capture_output=True, text=True, timeout=600)
results = json.loads(p.stdout.strip().splitlines()[-1]) if p.stdout.strip() else []
if not results:
    sys.exit(f"STOP: render produced no result\n{p.stderr[-2000:]}")
for r in results:
    print(os.path.basename(r["out"]), "ok" if r.get("ok") else "FAIL " + str(r.get("err"))[:300], r.get("nReactions", ""))
W, H = IDX["W"], IDX["H"] + 100
before, after = (open(f"{work}/renders/{s}.svg").read() for s in ("before", "after"))
open(f"{work}/renders/compare.svg", "w").write(
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{2 * W}" height="{H}" viewBox="0 0 {2 * W} {H}">'
    + before.replace("<svg ", '<svg x="0" y="0" ', 1) + after.replace("<svg ", f'<svg x="{W}" y="0" ', 1) + "</svg>")
chrome = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
for scen in list(TITLES) + ["compare"]:
    svg = f"{work}/renders/{scen}.svg"
    if not os.path.exists(svg):
        continue
    W = IDX["W"] * (2 if scen == "compare" else 1)
    scale = min(1.0, png_width * (2 if scen == "compare" else 1) / W)
    w, h = int(W * scale), int(H * scale)
    sized = re.sub(r'width="[^"]+" height="[^"]+"', f'width="{w}" height="{h}"', open(svg).read(), count=1)
    tmp = f"{work}/renders/{scen}_png.svg"
    open(tmp, "w").write(sized)
    subprocess.run([chrome, "--headless=new", "--disable-gpu", "--hide-scrollbars", f"--window-size={w},{h}",
                    f"--screenshot={work}/renders/{scen}.png", "file://" + tmp], capture_output=True, timeout=180)
    os.remove(tmp)
    print(scen, "png", w, h, os.path.getsize(f"{work}/renders/{scen}.png") // 1024, "KB; svg", os.path.getsize(svg) // 1024, "KB")
