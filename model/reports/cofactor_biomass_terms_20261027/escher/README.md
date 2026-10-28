# Escher figure: what the six cofactor biomass terms changed (2026-10-28)

Computation only; no model, curation or benchmark file changed. Requested by the user through the census
session (AGENTS.md section 7, "Side sessions": every finished session draws an Escher figure of what changed).

- **cofactor_terms_before_after.png / .svg**: the same map twice, *before* on the left and *after* on the right.
- **cofactor_terms_escher.html**: the self-contained interactive viewer (Escher 1.8.1 inlined, no internet needed;
  open it in Chrome, Safari or Firefox). Scenario drop-down: `before`, `after`, `census`. Hover or click a label
  for the pFBA and FVA values of both models and the census facts.

## Models, condition, solver

| | File (`model/candidates/`) | SHA-256 |
| --- | --- | --- |
| before | `E5_…_lipoylation_acylacp_20261023.xml` | `a12378210c7c66720d99fb7bf71d14cb7538ea6a3cf3e08919c85fb3c9c47d82` |
| after | `E5_…_lipoylation_acylacp_cofactor_biomass_20261027.xml` | `00e8e39775b64cfeb79cfafff65ef82eaf4760eff1b64fcddc3ce93b40919592` |

- **Condition:** Task A medium + PO1f, the same as the report's solves.
  - medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`);
  - strain profile `$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json` (`35307853…`);
  - both loaded with `load_effective_simulation_context`.
  - Checked in the run: uracil exchange R1354 at (−1000, 1000), R612 at (0, 0).
- **Solver:** Gurobi 13 at FeasibilityTol / OptimalityTol 1e-9 (`checked_solve.py`).
  - WT growth reproduces the recorded values (0.3124209 before, 0.3120049 after; the script stops otherwise).
  - pFBA uses cobra 0.30.0 `pfba`.
  - FVA uses `fraction_of_optimum = 1.0` on the 41 mapped reactions.
  - 168 LPs in all.

## What is drawn

- **Reactions (40, plus biomass_C):** the union of the reactions the gap census found blocked under a1237821 and
  unblocked under 00e8e397 in the census comparisons (`../census_compare_open.json`, `../census_compare_taskA.json`):
  - open condition (every exchange open), 756 → 722 blocked: **37 freed**;
  - Task A + PO1f: **35 freed**;
  - 32 in both;
  - open only: R578, R64, R811, R990, R991;
  - Task A only: R1146, R1147, R1314.
- **Families:** the reactions are grouped by cofactor family (folate 19, flavins 11, NAD/NADP 5, CoA 5).
- **biomass_C** is drawn with its six new terms only: THF m479 6.34e-5, FAD m528 1e-5, riboflavin m926 9.9e-4,
  NAD m122 2.65e-3, NADP m41 5.7e-4 and CoA m69 1.9e-4 mmol/gDW. It has 55 terms in all.
- **Not drawn:** the three GLPK-only "newly blocked" calls of the open census (R1343, R1344, R970). Gurobi FVA at
  1e-9 already refuted them, with maximum 1.47e-4 in both models (`../fva_gurobi_open_check.json`).
- **Layout:** automatic, from the gap-census builder (`scripts/maplib.py`). It carries no biological meaning.
  Repeated nodes are the same pool.

## Scenarios

| Scenario | Reaction colour and label |
| --- | --- |
| `before` | pFBA flux in a1237821 (mmol/gDW/h). biomass_C is left uncoloured (light grey) because a1237821 lacks the six terms drawn on it. |
| `after` | pFBA flux in 00e8e397. |
| `census` | Where the reaction was freed: amber = open only, pink = Task A only, blue = both. |

- **Flux colour scale:**
  - grey = 0 (|v| < 1e-12);
  - light → dark blue at 1e-6, 1e-5, 1e-4, 1e-3;
  - purple = biomass flux.
  - The gap-census palette is used throughout.
- **Node colour (open census):**
  - blue = gap metabolite under a1237821 and not under 00e8e397;
  - orange = gap under both;
  - grey = not a gap.

## Numbers shown on the map (from `flux.json`)

- **Before:** 0 of the 40 freed reactions carry pFBA flux.
- **After:** **24 of 40 carry pFBA flux**, for example:
  - the folate chain 1.98e-5;
  - riboflavin uptake R1314 / R1316 3.12e-4;
  - FMN / FAD 3.12e-6;
  - NAD synthesis R580 / R568 1.0e-3;
  - NADH kinase R571 1.78e-4;
  - the CoA chain 5.93e-5.
- **FVA at the optimum:** 25 of 40 can carry flux after the change. R565 (NAD kinase, range [0, 1.78e-4]) is the
  extra one, an alternative to R571.
- **Freed but zero at the optimum (FVA [0, 0]), 15 reactions:**
  - the de novo riboflavin route R378, R1307, R12, R11, R34, R699, R700 (the medium supplies riboflavin, so
    uptake is used at the optimum);
  - the duplicate dihydropteroate synthase R278;
  - R14, R137;
  - the 5-formyl-THF set R64, R811, R990, R991 (freed in the open condition only);
  - R578 (open only).

  These reactions are unblocked (they can carry flux at some growth), but not at optimal growth on Task A.
- **Gap metabolites on the map:** 26 were gaps under a1237821 and are not gaps under 00e8e397; 38 were never gaps;
  none is a gap in both models.

## Verification and limits

- **Checked this session:**
  - model SHAs;
  - WT growth of both models against the recorded values;
  - the medium checks above;
  - every viewer scenario draws 41 reactions and 254 metabolite nodes in a headless jsdom check
    (`scripts/check_viewer.js`);
  - the *before* render has no biomass-coloured (purple) segment.
- **Not checked:**
  - the viewer was not opened in a real browser (the app's browser pane cannot open `file://` pages);
  - the tooltip extension was not exercised headlessly.
- **What pFBA shows:** pFBA shows one optimal flux distribution. The FVA range in the tooltip says whether a zero
  is forced or only chosen.

## Build

The scripts in `scripts/` run from the worktree root.

1. `flux.py` (project `.venv`, cobra + Gurobi; run from `platform/` with
   `PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts`) → `flux.json`.
2. `build.py` (any Python with networkx; miniconda 3.x here, no LP) → map JSON, scenario data and `index.json`.
3. `render.py` → SVG via jsdom + Escher, PNG via headless Chrome, 2400 px per panel.
4. `build_html.py` → the viewer.
5. `check_viewer.js` (`ESC=<work>/esc`) → the headless check.

- **Copied from the gap census, with path changes only:** `maplib.py`, `escher_kernel.py`, `render_jobs.js`,
  `check_viewer.js`.
- **Work folder:** `$IYALI26_RESEARCH_ROOT/artifacts/tasks/cofactor_biomass_terms_20261027/escher/`. It holds:
  - `esc/node_modules`: Escher 1.8.1, `escher.min.js` SHA-256 `38c105f0…`, and jsdom 24.1.3;
  - `maps/`, `data/`, `index.json`;
  - `renders/`: separate before / after / census PNG and SVG at 2400 px.
