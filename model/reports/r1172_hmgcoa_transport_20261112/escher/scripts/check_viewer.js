// Headless check of the viewer: load the HTML in jsdom (scripts run), draw every module x scenario, count drawn reactions and coloured nodes.
// Copied from the gap census; jsdom now comes from $ESC/node_modules. usage: ESC=... node check_viewer.js VIEWER.html
const fs = require('fs'), path = require('path'); const { JSDOM } = require(path.join(process.env.ESC, 'node_modules/jsdom'));
const html = process.argv[2]; const W = 1600, H = 1000;
const dom = new JSDOM(fs.readFileSync(html, 'utf8'), { url: 'file://' + path.resolve(html), runScripts: 'dangerously', resources: 'usable', pretendToBeVisual: true,
  beforeParse(win) {
    const textLen = el => (el.textContent || '').length * 7.5;
    win.SVGElement.prototype.getBBox = function () { return { x: 0, y: 0, width: textLen(this), height: 16 }; };
    win.SVGElement.prototype.getComputedTextLength = function () { return textLen(this); };
    win.Element.prototype.getBoundingClientRect = function () { return { x: 0, y: 0, top: 0, left: 0, right: W, bottom: H, width: W, height: H }; };
    Object.defineProperty(win.HTMLElement.prototype, 'clientWidth', { get() { return W; } });
    Object.defineProperty(win.HTMLElement.prototype, 'clientHeight', { get() { return H; } });
    win.matchMedia = () => ({ matches: false, addListener() {}, removeListener() {} });
    win.ResizeObserver = class { observe() {} unobserve() {} disconnect() {} }; win.alert = () => {}; win.confirm = () => true;
  } });
dom.window.addEventListener('load', () => { setTimeout(async () => {
  const win = dom.window, doc = win.document, out = [];
  const m = doc.getElementById('m'), s = doc.getElementById('s');
  const wait = ms => new Promise(r => setTimeout(r, ms));
  for (const mo of [...m.options]) { m.value = mo.value; win.fillScen(); for (const so of [...s.options]) {
    m.value = mo.value; s.value = so.value; win.draw(); await wait(350);
    const svg = doc.querySelector('#c svg'); const nR = svg ? svg.querySelectorAll('.reaction').length : 0;
    const nodes = svg ? [...svg.querySelectorAll('.node circle')] : []; const cols = new Set(nodes.map(c => c.getAttribute('style') || ''));
    out.push({ module: mo.value, scenario: so.value, reactions: nR, metabolite_nodes: nodes.length, node_colours: cols.size, narrative: (doc.getElementById('n') || { textContent: '' }).textContent.length });
  } }
  console.log(JSON.stringify(out)); process.exit(0); }, 2000); });
setTimeout(() => { console.log('TIMEOUT'); process.exit(1); }, 300000);
