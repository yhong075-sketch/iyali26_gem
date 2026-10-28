"""Self-contained interactive viewer for the UM7 nuclear-diffusion map (Escher 1.8.1 inlined; adapted from the dolichol-route figure).

Usage: python build_html.py WORKDIR OUT.html
Scenario drop-down: capacity_before / capacity_after (FVA at 0 % on Task A), before / after (pFBA flux) and status
(Task A blocked status). Hover or click a label for Escher's tooltip, extended with the pFBA, FVA, bounds and status of both models.
"""
import json
import os
import sys

from palette import MET_SCALE, RXN_SCALES, TITLES

work, out = sys.argv[1:]
IDX = json.load(open(f"{work}/index.json"))


def rnd(o):
    if isinstance(o, float):
        return round(o, 1)
    if isinstance(o, dict):
        return {k: rnd(v) for k, v in o.items()}
    if isinstance(o, list):
        return [rnd(v) for v in o]
    return o


MAP = rnd(json.load(open(f"{work}/maps/um7_nuclear.json")))
DATA = {s: json.load(open(f"{work}/data/{s}.json")) for s in TITLES}
ESCHER_JS = open(f"{work}/esc/node_modules/escher/dist/escher.min.js").read()
_bad = 'if(!t)throw new Error("Automatic publicPath is not supported in this browser");'
ESCHER_JS = ESCHER_JS.replace(_bad, 'if(!t)t="https://unpkg.com/escher@1.8.1/dist/";') if _bad in ESCHER_JS else ESCHER_JS
ESCHER_JS = ESCHER_JS.replace("</script>", "<\\/script>")
so = "".join(f'<option value="{k}">{k}</option>' for k in TITLES)
sw = lambda c, r=6: f'<span class="sw" style="background:{c};border-radius:{r}px"></span>'
html = f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>iYali26 UM7 nuclear diffusion — 9d655c2b vs candidate</title>
<style>body{{font-family:sans-serif;margin:0}}#bar{{padding:8px 12px;background:#f4f4f4;border-bottom:1px solid #ccc;font-size:13px}}#c{{width:100vw;height:calc(100vh - 130px)}}
#t{{margin-top:4px;color:#333}}#c .legend-container,#c .vertical-legend-container,#c .legend-config-modal,#c .legend-config-modal-overlay{{display:none !important}}#c svg .reaction-label{{font-size:20px}}#c svg .node-label{{font-size:16px}}#c svg .text-label .label{{font-size:15px;font-style:italic;fill:#5b3a8c}}.sw{{display:inline-block;width:11px;height:11px;margin:0 3px 0 10px;vertical-align:middle}}</style>
<script>{ESCHER_JS}</script></head><body>
<div id="bar"><b>UM7 nuclear NAD block in iYali26</b>: 9d655c2b (9d655c2bdbf82721…) to the candidate with R_ADPRtn and R_PPItn; 23 reactions &nbsp; <select id="m" style="display:none"><option value="um7_nuclear">um7_nuclear</option></select>Scenario: <select id="s">{so}</select>
&nbsp; flux: {sw('#b5b5b5', 2)}0 {sw('#9cc2ee', 2)}1e-6 {sw('#5d9be3', 2)}1e-5 {sw('#2a78d6', 2)}1e-4 {sw('#0d366b', 2)}1e-2 &nbsp;|&nbsp; capacity: {sw('#b5b5b5', 2)}0 {sw('#9cc2ee', 2)}1e-6 {sw('#5d9be3', 2)}1e-2 {sw('#2a78d6', 2)}1 {sw('#0d366b', 2)}100 &nbsp;|&nbsp; status: {sw('#2a78d6', 2)}freed {sw('#eb6834', 2)}blocked in both {sw('#3d3d3d', 2)}closed (0,0) in both {sw('#6250d6', 2)}new {sw('#9cc2ee', 2)}unblocked in both
&nbsp;|&nbsp; nodes (Task A): {sw('#2a78d6')}all reactions blocked before, not after {sw('#eb6834')}all blocked in both {sw('#e87ba4')}only after {sw('#898781')}has an unblocked reaction
<div id="t"></div></div><div id="c"></div>
<script>const MAP={json.dumps(MAP, separators=(',', ':'))};const DATA={json.dumps(DATA, separators=(',', ':'))};const TITLES={json.dumps(TITLES)};
const MD={json.dumps(IDX['met_data'])};const MS={json.dumps(MET_SCALE)};const RS={json.dumps(RXN_SCALES)};const TIPS={json.dumps(IDX['tips'], separators=(',', ':'))};
function fillScen(){{}}
function draw(){{if(typeof escher==='undefined'){{document.getElementById('t').innerHTML='<b style="color:#b00">The Escher library did not load.</b> Open this file in a normal browser (Chrome, Safari, Firefox) by double-clicking it; it is self-contained and needs no internet.';return;}}
const s=document.getElementById('s').value;document.getElementById('t').innerHTML='<b>'+TITLES[s]+'</b>. Task A medium + PO1f, Gurobi 1e-9. Hover or click a label for values of both models; drag to pan, wheel to zoom.';
document.getElementById('c').innerHTML='';
window.B=escher.Builder(MAP,null,null,escher.libs.d3_select('#c'),{{menu:'zoom',enable_editing:false,reaction_data:DATA[s],reaction_scale:RS[s],reaction_styles:(s==='status'?['color','size']:['color','size','text','abs']),reaction_scale_preset:false,
metabolite_data:MD,metabolite_scale:MS,metabolite_styles:['color','size'],metabolite_scale_preset:false,metabolite_no_data_color:'#e0e0e0',metabolite_no_data_size:8,
reaction_no_data_color:'#d9d9d9',reaction_no_data_size:6,identifiers_on_map:'name',never_ask_before_quit:true,show_gene_reaction_rules:true,fill_screen:false,enable_search:true}});}}
function augment(){{const tc=document.querySelector('#c #tooltip-container, #c .tooltip-container');if(!tc)return;const idEl=tc.querySelector('.default-tooltip .id');if(!idEl)return;const id=idEl.textContent.trim();const old=tc.querySelector('.facts');if(old&&old.dataset.id===id)return;if(old)old.remove();const d=document.createElement('div');d.className='facts';d.dataset.id=id;d.style.cssText='margin:6px 0;font-size:12px;color:#333;max-width:420px;line-height:1.3';d.textContent=TIPS[id]||'';const btn=tc.querySelector('button');if(btn)btn.remove();tc.querySelector('.default-tooltip').appendChild(d);}}
new MutationObserver(function(){{setTimeout(augment,30);}}).observe(document.getElementById('c'),{{childList:true,subtree:true,attributes:true,attributeFilter:['style'],characterData:true}});
const P=new URLSearchParams(location.search);if(P.get('scenario')&&TITLES[P.get('scenario')])document.getElementById('s').value=P.get('scenario');else document.getElementById('s').value='capacity_after';
document.getElementById('s').onchange=draw;window.addEventListener('load',draw);</script></body></html>"""
open(out, "w").write(html)
print(out, os.path.getsize(out) // 1024, "KB")
