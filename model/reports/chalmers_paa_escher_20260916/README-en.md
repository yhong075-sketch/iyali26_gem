# Chalmers phenylacetic acid branch — Escher

Open `chalmers-paa-escher-en.html` and switch between W29 and S. cerevisiae using the buttons at the top. No internet connection is needed; keep `escher.min.js` in the same folder as the HTML file.

Zoom or pan the map, click a reaction for its equation, saved bounds, GPR and evidence status, or enable **Edit layout** to move nodes and labels. Layout changes persist when switching models within the current page. Download the map JSON or export an SVG to save your changes.

- `w29-map-en.json` and `sc-map-en.json`: maps that can be imported into Escher.
- `w29-escher-en.svg/png` and `sc-escher-en.svg/png`: vector and raster exports of the initial layout.
- Red: blocked by complete steady-state metabolite balances. Teal: transport or exchange recorded in the model. Gray: other reactions. No simulated or measured flux data are overlaid.

The maps use the frozen Chalmers iYali v4.1.2 and Yeast-GEM v9.1.1 snapshots audited on September 14, 2026. Later releases were not checked. The maps contain 8 W29 reactions and 10 S. cerevisiae reactions with the source stoichiometry, bounds and GPRs preserved. Only the 0851 transaminase branch is shown upstream, and downstream 2-phenylethanol reactions are omitted. All reactions involving PAA, PAM and phenylacetaldehyde in its three compartments were checked against the source models.

Both snapshots lack an independent PAA outlet. The transport and exchange reactions shown move **phenylacetaldehyde**, not PAA. The PAM pathway is not included in this W29 model; this does not establish biological absence. In S. cerevisiae, the reversible PAM reaction is blocked by its complete PAM balance, which also blocks PAA formation. Model GPR assignments do not establish experimental substrate specificity.

The English version changes display text and label placement only. It does not change the source models, reaction definitions, GPRs, culture conditions or balance proofs. No FBA or other optimization was performed.

Full source identities and SHA-256 hashes are recorded in `provenance.json` and `evidence-en.json`; translation checks are recorded in `english-provenance.json`. Independent source and figure review is recorded in `AUDIT.md`, and browser interaction checks in `render-check-en.json`. Escher 1.8.1 is distributed with its license.

Model sources: [Chalmers W29 iYali](https://github.com/SysBioChalmers/Yarrowia_lipolytica_W29-GEM) and [Chalmers Yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM). Figure preparation: September 16, 2026.
