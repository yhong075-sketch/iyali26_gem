# Escher maps of the gap modules and problem classes (candidate a1237821, Ponce-de-León 2013 census)

Open **gap_modules_escher.html** in a normal browser (Chrome, Safari, Firefox; double-click the file). It is
self-contained: the Escher 1.8.1 bundle is inlined, so it needs no other file and no internet. If the page stays
blank, the environment did not run its scripts (some in-app file previews do not); use a browser. No server, no
optimisation: every value shown comes from the LP tests of this report.

- **Map drop-down.** *Problem classes (across modules)*: four maps that gather the reactions of one primary cause
  from every module, each split into labelled sections — `cofactor` (72 reactions by cofactor family; the tested
  drains are drawn as TEST reactions with the number of reactions they unblock), `closures` (the 69 reactions at
  (0,0) grouped as in the report, plus the 18 reactions that opening one closure frees), `transport` (52 reactions
  whose best single fix is a transport link, grouped by the root's compartment, plus 17 extracellular roots without
  an exchange; the tested links and exchanges are drawn as TEST reactions) and `deadend` (181 reactions whose best
  single fix is a source/sink on one root, grouped by root class; the tested sinks are drawn as TEST reactions).
  *Unconnected modules*: the 20 largest modules, UM1–UM20 (436 of the 756 blocked reactions). The 192 smaller
  modules (1–5 reactions each) are not mapped as modules; their reactions appear on the class maps when they
  belong to one of the four classes.
- **Scenario drop-down** (reaction colour and label value; the list depends on the map). Class maps: `drains`,
  `open_alone`, `link` or `sink` = the class evidence (model reactions: the number of tests that unblock them, or 1;
  TEST reactions: the number of reactions the test unblocks) and `blocked`. Module maps: `blocked` = the census
  state (0 = blocked, −1 = closed at (0,0), drawn black); `single_fixes` = number of single elementary-fix tests
  (root sink, transport link, closure opening, exchange, duplicate merge) that unblock the reaction; `all_fixes` =
  1 when every elementary fix of the module applied at once unblocks it; `targeted` = 1 when the module's
  targeted hypothesis test unblocks it (named in the bar); `drains` = 1 when draining an undrained cofactor
  unblocks it.
- **Node colour = gap class** (legend in the bar): RNP never produced, RNC never consumed, single reversible
  reaction, isolated, propagated DNP/UNC, gap without root (coupling); grey = active context metabolite that is
  not a gap; small light nodes = currency metabolites. Exact stoichiometry is kept for every model reaction;
  exchanges and demands are drawn against a pseudo species "(outside)" / "(sink)", and TEST reactions (sub-label
  "no gene (test reaction) · TEST (not in model)") are the tested drain / link / sink / exchange, not model
  content. The biomass-like reactions R1371, R1372, R1387, R1710 show their gap metabolites only.
- **Layout is automatic** and carries no biological meaning: reactions and metabolites in alternating layers;
  hub metabolites used by five or more reactions of a module or section are repeated per reaction (repeated nodes
  are the same pool, the usual Escher convention); chains longer than 16 layers are folded into columns with the
  boundary metabolite repeated ("continues in next column"); chains are packed in rows, class maps in sections.
- **Hover or click a label** for Escher's tooltip: id, name, scenario value, one line of census facts (gap class,
  module, primary cause, best single fix, drains) and a button **"Open … in the gap dataset"** that opens the entry
  in `gap_dataset.html` (same folder): every metabolite (1,882) and every blocked reaction (756) with its census
  facts, LP fix evidence, metabolites and "on maps" links back to the viewer (`gap_modules_escher.html?map=UM9`
  preselects a map, `&scenario=drains` a scenario). The dataset page has a find box; TEST reactions have entries
  that say what was tested. Drag to pan, mouse wheel to zoom. The
  sub-label of every model reaction reads "GPR · note" (the primary cause on module maps, the test evidence on
  class maps).
- **Static renders:** `renders/<map>.svg` (vector, full size) and `renders/<map>.png` (≤ 6000 px wide) for every
  map in the artifacts folder (modules: scenario `single_fixes`; class maps: their class scenario); three module
  PNGs are kept in the git report folder.
- **Build:** `export_modules.py` / `export_all.py` (cobra; reactions, metabolites, classes and test results) →
  `build_maps.py` (module maps) and `build_class_maps.py` with `maplib.py` (class maps: sectioned layout + Escher
  1-0-0 map JSON, scenario data) → `render_all.py` (jsdom render of the Escher bundle to SVG, headless Chrome to
  PNG) → `build_html.py` (viewer) → `check_viewer.js` (jsdom: every map × scenario draws; result in
  `viewer_check.json`). Scripts and the skill helpers are in `scripts/`.
