# Escher figure: aconitase rules before and after (2026-10-29)

Drawn after the build. The change is a GPR edit, so WT growth is the same in both models; the difference shows when the
ACO1-type gene `YALI1D11984g` is knocked out.

## Files

- **aconitase_ko_yali1d11984g_before_vs_after.png / .svg**: the same map twice, pFBA with `YALI1D11984g` knocked out.
  - Left, before (`059316ef`): only R239 (cytosol, `YALI1D11984g` alone) closes. R238 and R216 stay open through the
    OR partner `YALI1E17935g` and carry 1.57 (WT pFBA 1.49); growth 0.3117 = WT.
  - Right, after (`e1e7c868`): R238, R216, R1556 and R239 are all (0, 0). Growth 0.0076119, KO/WT 0.0244. Citrate
    synthase R236 falls from 1.57 to 0.0167, and the NAD and NADP isocitrate dehydrogenases R487 and R1385 go to 0.
- `flux.json`: the LP results and the mapped reactions; `index.json`: map counts and tooltips.
- `trace_residual.json`: one extra checked solve and pFBA of the knockout on the candidate (below).

## Model, condition, solver

- **Models:** `059316ef34417c88…` (before) and `e1e7c8686c1fbbbb…` (after), read from SHA-checked short copies
  (`../logs/sweep_short_copies_sha256.txt` lists the same files).
- **Condition:** Task A medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) + PO1f profile (`35307853…`).
- **Solver:** Gurobi 13, 1e-9 (`checked_solve.py`), real HOME; per scenario a checked solve then pFBA (all optimal);
  WT 0.3116608216 in both models (the script stops otherwise). 12 LPs for the figure, 3 for the trace; the first trace
  run stopped at the JSON write (numpy bool) after its solves, so 18 LPs in all.

## Which reactions are drawn

- **Aconitase steps:** R238, R216, R1556 (rules changed), R239, and the spontaneous cis- to trans-aconitate steps
  R735, R736.
- **Every other reaction of citrate, cis-aconitate and isocitrate** in the cytosol and the mitochondrion, selected by
  metabolite membership in `flux.py` (biomass assembly and `x*` pools excluded): R236, R1894, R487, R1385, R488, R490,
  R831, R832, R833, R817, R1044, R1046, R1047, R1875.
- **The 12 other reactions whose knockout flux changes most** between the models: ATP synthase R171, the respiratory
  chain R2062 / R305, transports, and the serine/cysteine pFBA cycle R659 / r0210 / R_CYST_II (the cycle of about 18
  already recorded, untraced, in the R171 round). These show the growth collapse, not a mechanism.

## Residual growth (trace, one pFBA optimum; verified as a computation, the route is not claimed unique)

With all four aconitase steps closed and no glutamate in the medium, arginine-derived ornithine feeds glutamate
through ornithine transaminase R610 (0.0053) and P5C dehydrogenase R9 (0.0037). L-arginine (0.00624) and L-aspartate
(0.01307) uptakes sit at their bounds; why aspartate is limiting was not traced. A large aspartate / glutamate
transaminase shuttle (R161, R162, R823, R1386, about 1.27) runs between cytosol and mitochondrion with no net carbon
gain. So the 0.0244 ratio is probably set by the finite amino-acid bounds of the Task A medium, which `benchmark.md`
calls concentration-scaled assumptions, not uptake rates, and the 1 % class of this gene depends on them (inferred:
no shadow prices taken, bounds not varied).

## Build

As for the R673 figure (`../../r673_prolyl_split_20261029/escher/README.md`): `flux.py` and `trace_residual.py` from
`platform/` with `PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts`, then
`build.py ../flux.json WORKDIR` and `render.py WORKDIR 2400` from `scripts/`, all with `/opt/miniconda3/bin/python`.
`WORKDIR` is `$IYALI26_RESEARCH_ROOT/artifacts/tasks/aconitase_aco2_split_20261029/escher/` (single-scenario renders
there); `esc/node_modules` links to the cofactor-terms Escher 1.8.1 install. Copied from the R673 figure: `maplib.py`,
`escher_kernel.py`, `render_jobs.js`, `render.py` unchanged (map file names `maps/r171*.json` kept); `flux.py`,
`build.py`, `palette.py` edited (knockout gene, reaction set, sections, titles); `trace_residual.py` new.

## Limits

Automatic layout without biological meaning; no interactive viewer; pFBA is one optimum; the compare PNG was
inspected by eye (it is 1.0 MB, at the size limit for the report folder), the SVGs were not opened in a browser.
