// Render Escher maps headlessly with the real Escher bundle in jsdom. usage: node render_jobs.js jobs.json
const fs = require('fs'); const { JSDOM } = require('jsdom');
const ESCHER = fs.readFileSync(__dirname + '/node_modules/escher/dist/escher.min.js', 'utf8');
const jobs = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
const ROLE_CSS = { comp: 'font-size:30px;font-weight:bold;font-style:normal;fill:#4a4a4a',
                   gpr: 'font-size:15px;font-weight:normal;font-style:italic;fill:#5b3a8c',
                   note: 'font-size:20px;font-weight:normal;font-style:normal;fill:#444',
                   lane: 'font-size:22px;font-weight:bold;font-style:normal;fill:#222' };
const CSS_APPEND = `
svg.escher-svg .reaction-label { font-size: 20px; }
svg.escher-svg .node-label { font-size: 16px; }
svg.escher-svg .text-label .label { font-size: 18px; font-style: normal; font-weight: normal; fill: #333; }`;
function renderOne(job) {
  return new Promise(resolve => {
    const W = job.W || 2400, H = job.H || 1800;
    const dom = new JSDOM(`<!DOCTYPE html><html><head><script src="file://${__dirname}/node_modules/escher/dist/escher.min.js"></script></head><body><div id="map-container" style="width:${W}px;height:${H}px"></div></body></html>`,
      { runScripts: 'outside-only', pretendToBeVisual: true });
    const win = dom.window; const textLen = el => (el.textContent || '').length * 7.5;
    win.SVGElement.prototype.getBBox = function () { return { x: 0, y: 0, width: textLen(this), height: 16 }; };
    win.SVGElement.prototype.getComputedTextLength = function () { return textLen(this); };
    win.Element.prototype.getBoundingClientRect = function () { return { x: 0, y: 0, top: 0, left: 0, right: W, bottom: H, width: W, height: H }; };
    Object.defineProperty(win.HTMLElement.prototype, 'clientWidth', { get() { return W; } });
    Object.defineProperty(win.HTMLElement.prototype, 'clientHeight', { get() { return H; } });
    win.matchMedia = () => ({ matches: false, addListener() {}, removeListener() {} });
    win.ResizeObserver = class { observe() {} unobserve() {} disconnect() {} };
    win.alert = () => {}; win.confirm = () => true;
    win.eval(ESCHER);
    const mapData = JSON.parse(fs.readFileSync(job.map, 'utf8'));
    const deco = JSON.parse(fs.readFileSync(job.deco, 'utf8'));
    const options = Object.assign({ menu: 'none', enable_editing: false, enable_keys: false, enable_tooltips: false, scroll_behavior: 'none',
      never_ask_before_quit: true, identifiers_on_map: 'name', show_gene_reaction_rules: false, full_screen_button: false, fill_screen: false,
      reaction_styles: ['color', 'size', 'text', 'abs'], reaction_scale_preset: false }, job.options);
    const timer = setTimeout(() => resolve({ out: job.out, ok: false, err: 'timeout' }), 10000);
    options.first_load_callback = function () {
      try {
        const doc = win.document, NS = 'http://www.w3.org/2000/svg';
        const c = mapData[1].canvas, top = 100;
        const svg = doc.querySelector('#map-container svg');
        svg.querySelectorAll('.canvas-group, .mouse-node-group, #mouse-node, .brush, .direction-arrow, .bezier-group, .midmarker-circle, .multimarker-circle').forEach(n => n.remove());
        const zg = svg.querySelector('.zoom-g'); zg.removeAttribute('transform');
        const st = svg.querySelector('style'); st.textContent += CSS_APPEND;
        const bg = doc.createElementNS(NS, 'g'); bg.setAttribute('class', 'compartments');
        const white = doc.createElementNS(NS, 'rect');
        Object.entries({ x: c.x, y: c.y - top, width: c.width, height: c.height + top, fill: '#ffffff' }).forEach(([k, v]) => white.setAttribute(k, v)); bg.appendChild(white);
        for (const r of deco.rects) {
          const e = doc.createElementNS(NS, 'rect');
          Object.entries({ x: r.x, y: r.y, width: r.w, height: r.h, rx: 18, fill: r.fill, stroke: r.stroke, 'stroke-width': 3 }).forEach(([k, v]) => e.setAttribute(k, v)); bg.appendChild(e);
        }
        zg.insertBefore(bg, zg.firstChild);
        for (const [lid, role] of Object.entries(deco.label_roles || {})) {
          const t = svg.querySelector(`g#${lid} text`); if (t && ROLE_CSS[role]) t.setAttribute('style', ROLE_CSS[role]);
        }
        const tt = doc.createElementNS(NS, 'text'); tt.textContent = job.title;
        tt.setAttribute('x', c.x + 30); tt.setAttribute('y', c.y - top + 58);
        tt.setAttribute('style', 'font-family:sans-serif;font-size:40px;font-weight:bold;fill:#111'); zg.appendChild(tt);
        svg.setAttribute('xmlns', NS);
        svg.setAttribute('viewBox', `${c.x} ${c.y - top} ${c.width} ${c.height + top}`);
        svg.setAttribute('width', c.width); svg.setAttribute('height', c.height + top);
        fs.writeFileSync(job.out, svg.outerHTML);
        const labels = [...svg.querySelectorAll('.reaction-label')].map(t => t.textContent);
        clearTimeout(timer); resolve({ out: job.out, ok: true, nReactions: svg.querySelectorAll('.reaction').length, labels });
      } catch (e) { clearTimeout(timer); resolve({ out: job.out, ok: false, err: String(e.stack || e) }); }
    };
    try { win.escher.Builder(mapData, null, null, win.escher.libs.d3_select('#map-container'), options); }
    catch (e) { clearTimeout(timer); resolve({ out: job.out, ok: false, err: 'Builder: ' + String(e.stack || e) }); }
  });
}
(async () => { const res = []; for (const j of jobs) res.push(await renderOne(j)); console.log(JSON.stringify(res)); process.exit(0); })();
