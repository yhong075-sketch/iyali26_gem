"""Open the viewer in headless Chrome once per scenario (?scenario=), dump the DOM and screenshot it (read-only).

Usage: python check_viewer_chrome.py VIEWER.html WORKDIR OUT.json
Counts drawn reaction groups and metabolite circles, reads the status line #t (the Escher-not-loaded warning is
shown there), and records three reaction labels. Screenshots go to WORKDIR/viewer_<scenario>.png.
"""
import json
import os
import re
import subprocess
import sys

viewer, work, out = sys.argv[1:4]
chrome = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
url = "file://" + os.path.abspath(viewer)
rows = []
for s in ("before", "after", "census"):
    prof = f"{work}/chrome_profile_{s}"
    common = [chrome, "--headless=new", "--disable-gpu", "--virtual-time-budget=15000", "--window-size=1600,1000", f"--user-data-dir={prof}"]
    dom = subprocess.run(common + ["--dump-dom", f"{url}?scenario={s}"], capture_output=True, text=True, timeout=180).stdout
    subprocess.run(common + ["--hide-scrollbars", f"--screenshot={work}/viewer_{s}.png", f"{url}?scenario={s}"], capture_output=True, timeout=180)
    status = re.search(r'<div id="t">(.*?)</div>', dom, re.S)
    status = re.sub(r"<[^>]+>", "", status.group(1)) if status else ""
    labels = {l.split(" ")[0]: l for l in re.findall(r'<text[^>]*class="reaction-label label"[^>]*>(.*?)</text>', dom)}
    rows.append({"scenario": s, "status_line": status[:160], "escher_not_loaded": "did not load" in status,
                 "reaction_groups": len(re.findall(r'class="reaction"', dom)), "node_circles": len(re.findall(r'class="node-circle', dom)),
                 "labels": {k: labels.get(k, "")[:60] for k in ("R_DOLSYN", "R2178", "R287", "R2180")}})
json.dump(rows, open(out, "w"), indent=1)
print(json.dumps(rows, indent=1))
