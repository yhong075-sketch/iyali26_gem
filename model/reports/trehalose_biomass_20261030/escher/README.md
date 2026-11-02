# Escher figure: trehalose in `biomass_C`, before and after (2026-10-30)

Drawn after the build and sweep. The only model difference is the removed `biomass_C` trehalose term (and its notes),
so the figure shows the trehalose branch at WT and with the TPS1 gene `YALI1E17616g` knocked out, in both models.

## Files

- **trehalose_wt_before_vs_after.png / .svg**: WT, pFBA.
  - Left, before (`e1e7c868`): growth 0.3116608; R143 and R767 carry 0.002138 = 0.00686 × growth into the
    `biomass_C` trehalose term (drawn as its own arrow, "biomass_C trehalose term"; it is a coefficient of
    `biomass_C`, not a model reaction).
  - Right, after (`d2e4c0c4`): growth 0.3122964 (+0.20 %); the branch carries no flux. Phosphoglucomutase R638 and
    UDP-glucose pyrophosphorylase R793 fall by 0.0015 (0.2962 → 0.2946) while 1,3-β-glucan synthase R4 rises by
    0.0006 (0.2940 → 0.2946): before, R793 exceeded R4 by exactly the 0.00214 of UDP-glucose taken by R143; after,
    the two are equal.
  - The drain arrow is labelled "present before, removed after" in all panels (one shared map).
- **trehalose_tps1_ko_before_vs_after.png / .svg**: `YALI1E17616g` knocked out (R143 closed), pFBA.
  - Left, before: growth 0. The pFBA optimum still sends 0.328 glucose through hexokinase and G6P isomerase, which
    meets the fixed non-growth maintenance demand; UDP-glucose and β-glucan synthesis stop.
  - Right, after: growth 0.3122964 = WT.
- `flux.json`: LP results and the mapped reactions; `index.json`: map counts, branch fluxes and tooltips.

## Model, condition, solver

- **Models:** `e1e7c8686c1fbbbb…` (before) and `d2e4c0c4b02c7e48…` (after), read from the SHA-checked short copies
  of the sweep (`../logs/sweep_short_copies_sha256.txt`). The final candidate `164542fb…` (v2) differs from
  `d2e4c0c4` only in two `biomass_C` note texts ([diff](../model_diff_v2_vs_d2e4c0c4.json)), so it poses the same
  optimisation problem and these fluxes apply to it.
- **Condition:** Task A medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) + PO1f profile (`35307853…`).
- **Solver:** Gurobi 13.0.2 (`/opt/miniconda3/bin/python`), 1e-9 (`checked_solve.py`), real HOME; per scenario a checked solve then pFBA (all optimal);
  `flux.py` stops unless each WT equals its sweep value. 12 LPs.

## Which reactions are drawn

- **Trehalose branch:** every reaction of trehalose `m294[C_cy]` and trehalose 6-phosphate `m297[C_cy]` (biomass
  pseudoreactions excluded): R143, R767, R142, R1349, plus the export R1013, plus the `biomass_C` trehalose drain.
- **Supply:** the reactions of G6P `m296[C_cy]` and UDP-glucose `m8[C_cy]` that carry flux in any scenario: R387,
  R326, R638, R793, R4. The other 12 (zero in all four scenarios) are listed on the map and in `index.json`.

## Build

From `platform/`: `PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts
/opt/miniconda3/bin/python ../model/reports/trehalose_biomass_20261030/escher/scripts/flux.py MEDIUM PROFILE
../model/reports/trehalose_biomass_20261030/escher/flux.json E1E7C868_COPY TREH_COPY`; then from `scripts/`:
`build.py ../flux.json WORKDIR` and `render.py WORKDIR 2400`, all with `/opt/miniconda3/bin/python`. `WORKDIR` is
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/trehalose_biomass_20261030/escher/` (single-scenario renders there);
`esc/node_modules` links to the cofactor-terms Escher 1.8.1 install. Copied unchanged from the aconitase figure:
`maplib.py`, `escher_kernel.py`, `render_jobs.js`; `render.py` edited (two compared pairs from `palette.COMPARES`);
`flux.py`, `build.py`, `palette.py` rewritten for this branch. Map file names `maps/r171*.json` kept. The drain label
was changed from "TERM REMOVED" after the audit (it showed on the before panels too) and the maps re-rendered from the
same `flux.json` (no new LP).

## Limits

Automatic layout without biological meaning; no interactive viewer; pFBA is one optimum; the compare PNGs were
inspected by eye, the SVGs were not opened in a browser.
