# Escher figure: R673 prolyl-tRNA synthetase rule before and after (2026-10-29)

Drawn after the build. The change is a GPR edit, so the WT optimum is the same in both models; the difference shows
when the cytosolic enzyme `YALI1E05929g` is knocked out.

## Files

- **r673_ko_yali1e05929g_before_vs_after.png / .svg**: the same map twice, pFBA with `YALI1E05929g` knocked out.
  - Left, before (`ca43ae0c`): the OR partner `YALI1D19110g` keeps R673 open. Growth 0.3117 = WT; every mapped flux
    equals WT (`index.json` → `ko_before_equal_to_wt`). R673, the Pro-tRNA(Pro) biomass incorporation
    `TRNA_BIOMASS_R673` and P5C reductase R689 carry 0.0643.
  - Right, after (`6af03c40`): R673 is (0, 0), growth 0. R673 and `TRNA_BIOMASS_R673` are 0.
- `flux.json` holds the LP results; `index.json` the map counts and tooltips.

## Model, condition, solver

- **Models:** `ca43ae0c12e0…` (before) and `6af03c405a69…` (after), read from SHA-checked short copies
  (`../logs/sweep_short_copies_sha256.txt`).
- **Condition:** Task A medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) + PO1f profile (`35307853…`).
- **Solver:** Gurobi 13, 1e-9 (`checked_solve.py`), real HOME; per scenario a checked solve then pFBA (all optimal);
  WT 0.3116608216 in both models (the script stops otherwise). 12 LPs.

## Which reactions are drawn

R673 and every other reaction of cytosolic L-proline (`m765`), tRNA(Pro) (`m905`) and Pro-tRNA(Pro) (`m906`), except the
biomass assembly `newBiom` and the amino-acid pool `xAMINOACID`: `TRNA_BIOMASS_R673`, the biomass pseudoreactions R1387,
R1710, R1372, the proline transports R1225, R2040, R893, hydroxyproline reductase R509 and P5C reductase R689. Selected by metabolite membership in `flux.py`, not by
flux.

**Observation, not traced:** at the zero-growth optimum after the knockout, pFBA still routes about 0.005 through P5C
reductase R689 and the proline transport R1225 (cytosol → mitochondrion). Why the zero-growth pFBA optimum uses
mitochondrial proline was not investigated; it does not affect the knockout result (growth 0 is the optimum).

## Build

As for the R171 figure (`../../r171_atp_synthase_gpr_20261029/escher/README.md`): `flux.py` from `platform/` with
`PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts`, then `build.py ../flux.json WORKDIR` and
`render.py WORKDIR 2400` from `scripts/`, all with `/opt/miniconda3/bin/python`. `WORKDIR` is
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/r673_prolyl_split_20261029/escher/`; `esc/node_modules` links to the
cofactor-terms Escher 1.8.1 install. Copied from the R171 figure: `maplib.py`, `escher_kernel.py`, `render_jobs.js`,
`render.py` unchanged (the map file names `maps/r171*.json` are kept so `render.py` needs no edit); `flux.py`,
`build.py`, `palette.py` edited (knockout gene, reaction set, sections, titles, scale).

## Limits

Automatic layout without biological meaning; no interactive viewer; pFBA is one optimum; the compare PNG was
inspected by eye, the SVGs were not opened in a browser.
