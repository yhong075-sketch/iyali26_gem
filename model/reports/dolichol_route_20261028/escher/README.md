# Escher figure: what the dolichol-route change did (2026-10-28)

Computation only; no model, curation or benchmark file changed. Requested by the user through the census
session (AGENTS.md section 7, "Side sessions": every finished session draws an Escher figure of what changed).

- **dolichol_route_before_after.png / .svg**: the same map twice, *before* on the left and *after* on the right.
- **dolichol_route_census.png**: the same map coloured by open-condition census status.
- **dolichol_route_escher.html**: the self-contained interactive viewer (Escher 1.8.1 inlined, no internet needed;
  open it in Chrome, Safari or Firefox). Scenario drop-down: `before`, `after`, `census`; a URL parameter
  `?scenario=before` also works. Hover or click a label for the pFBA, FVA, bounds and census status of both models.

## Models, condition, solver

| | File (`model/candidates/`) | SHA-256 |
| --- | --- | --- |
| before | `E5_…_r2135_pool_r1726_cpt_20261027.xml` | `ce9e1112ee940531469457009234a2fd11de7f589c5948a7dec6b91118005294` |
| after | `E5_…_r1726_cpt_dolichol_route_v2_20261028.xml` | `9d655c2bdbf82721a38f148a9a4ee0a385bb509ea7b591111107113c4aa8335c` |

- **The after model was read from** `$IYALI26_RESEARCH_ROOT/artifacts/tasks/dolichol_route_20261028/builds/candidate_v2_seed1.xml`,
  the byte-identical seed-1 build. cobra on macOS rejects absolute SBML paths over 260 characters, and the
  candidate's path is longer. The loader reports the same canonical SHA-256 (`flux.json`).
- **Condition:** Task A medium + PO1f, the same as the report's solves.
  - medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`);
  - strain profile `$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json` (`35307853…`);
  - both loaded with `load_effective_simulation_context`.
  - Checked in the run (the script stops otherwise): uracil exchange R1354 at (−1000, 1000), R612 at (0, 0).
- **Solver:** Gurobi 13 at FeasibilityTol / OptimalityTol 1e-9 (`checked_solve.py`), real HOME.
  - WT growth reproduces the recorded values (0.3116824 before, 0.3116608 after; the script stops otherwise).
  - pFBA uses cobra 0.30.0 `pfba`. FVA uses `fraction_of_optimum = 1.0` on the mapped reactions. 162 LPs.
- **Census status** comes from the run's GLPK census files, not from new LPs:
  - before: the Kennedy re-run census of `ce9e1112` (`artifacts/tasks/kennedy_closures_20261027/rerun_census/`);
  - after: the census of `a801aad0` (`artifacts/tasks/dolichol_route_20261028/census/`). `a801aad0` differs from
    `9d655c2b` only in the R_DOLSYN notes and is the same optimization problem (`../v2_same_lp.json`).

## What is drawn

38 reactions plus biomass_C, in three sections:

| Section | Reactions |
| --- | --- |
| Cytosolic route | R282 (GPP synthase), R322 (FPP synthase, which R2179 duplicates), R2179, R2180, **R2178** (closed), **R_DOLSYN** (new), R287 (SEC59) |
| Dol-P / Dol-PP pools and their consumers | R1101, R2177, R2182, R2181, R288, R2102, R289, R372, R374, R375, R2106 (ALG7), R2138 (ALG5), biomass_C |
| Greyed: lipid-particle Srt1 chain (census UM6) | R217–R235, blocked before and after, not changed |

- **biomass_C** is drawn with its new Dol-P term (m503[C_er], 1e-4 mmol/gDW) only; it has 62 terms in all.
- **Closed reactions** carry "CLOSED (0,0)" in their sub-label, as in the census viewer; here only R2178.
- **R_DOLSYN** is absent from `ce9e1112`, so the *before* scenario has no value for it (Escher draws "(nd)").
- **Layout:** automatic, from the gap-census builder (`scripts/maplib.py`). It carries no biological meaning.
  Repeated nodes are the same pool; the Dol-P and Dol-PP pools repeat because they are hubs.

## Scenarios

| Scenario | Reaction colour and label |
| --- | --- |
| `before` | pFBA flux in ce9e1112 (mmol/gDW/h). biomass_C is uncoloured because ce9e1112 lacks the drawn term. |
| `after` | pFBA flux in 9d655c2b. |
| `census` | Open-condition status: blue = freed; orange = blocked in both; black = closed (0, 0); purple = new and unblocked; light blue = unblocked in both. No numbers. |

- **Flux colour scale** (gap-census palette): grey = 0 (|v| < 1e-12); light → dark blue at 1e-6, 1e-5, 1e-4,
  1e-2; purple = biomass flux.
- **Node colour (open census):** blue = gap metabolite under ce9e1112 and not under 9d655c2b; orange = gap under
  both; pink = gap only under 9d655c2b; grey = not a gap.

## Numbers on the map (`index.json`, `flux.json`)

| On the 38 mapped reactions | before (ce9e1112) | after (9d655c2b) |
| --- | ---: | ---: |
| In the model | 37 | 38 |
| Blocked, open census | 28 | 25 |
| Blocked, Task A census | 28 | 25 |
| Carry pFBA flux on Task A | 4 | 9 |
| Can carry flux at the optimum (FVA) | 6 | 10 |

- **Freed:** R2179, R2180 and R287. **Closed:** R2178 (blocked before, closed after; counted as blocked both
  times). **New and unblocked:** R_DOLSYN.
- **Flux after:** R2179, R2180, R_DOLSYN, R287 and R1101 carry 3.12e-5. R282 and R322 carry 0.0239 and R2102 /
  R289 0.0735 in both models.
- **Still blocked in both (24):** the 19 UM6 reactions, R2177 / R2181 / R2182 (Dol-PP has no producer until
  the N-glycan chain is linked), and the ALG7 / ALG5 steps R2106 / R2138.
- **Gap metabolites:** 3 stop being gaps (dehydrodolichyl-PP m2019, FPP′ m2020, dolichol m499); 29 are gaps
  in both; none becomes a gap.

## Verification and limits

- **Checked this session:**
  - model SHAs and WT growth of both models against the recorded values, plus the medium checks above;
  - every viewer scenario draws 39 reactions and 250 metabolite nodes, both in jsdom
    (`logs/check_viewer_jsdom.log`) and **in headless Chrome 156** (`logs/check_viewer_chrome.json`). Chrome
    shows the scenario's status line with no Escher-not-loaded warning, and the labels R_DOLSYN "(nd)" →
    3.12e-5 and R287 0 → 3.12e-5;
  - the renders were inspected; the *before* render shows R_DOLSYN as "(nd)" and no flux on the route.
- **Not checked:**
  - the tooltip extension was not exercised (no hover in headless runs);
  - the viewer was not opened in an interactive browser window.
- **Viewer vs renders:** the section background boxes, including the grey UM6 box, are drawn only in the PNG /
  SVG renders. The viewer shows the UM6 chain by its grey zero-flux colour and orange gap nodes, without the box.
- **What pFBA shows:** one optimal flux distribution. The FVA range in the tooltip says whether a zero is forced
  or only chosen.

## Build

The scripts in `scripts/` run from this folder's `scripts/` (flux.py from `platform/`).

1. `flux.py` (miniconda Python with cobra + Gurobi; from `platform/` with
   `PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts`) → `flux.json`.
2. `build.py FLUX.json CENSUS_BEFORE_DIR CENSUS_AFTER_DIR WORKDIR` (networkx, no LP) → map JSON, scenario data, `index.json`.
3. `render.py WORKDIR 2400` → SVG via jsdom + Escher, PNG via headless Chrome.
4. `build_html.py WORKDIR ../dolichol_route_escher.html` → the viewer.
5. `check_viewer.js` (`ESC=<work>/esc node check_viewer.js VIEWER`) and `check_viewer_chrome.py VIEWER WORKDIR OUT` → the checks.

- **Copied from the cofactor-terms figure** (itself from the gap census), with path and label changes only:
  `maplib.py`, `escher_kernel.py`, `render_jobs.js`, `check_viewer.js`, `render.py`, `build_html.py`.
- **Work folder:** `$IYALI26_RESEARCH_ROOT/artifacts/tasks/dolichol_route_20261028/escher/` holds `maps/`, `data/`,
  `index.json`, `renders/` (separate before / after / census PNG and SVG), the Chrome screenshots
  `viewer_<scenario>.png`, and `esc/node_modules`, a link to the cofactor task's install (Escher 1.8.1,
  `escher.min.js` SHA-256 `38c105f0…`; jsdom 24.1.3).
