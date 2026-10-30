# Escher figure: R171 ATP synthase rule before and after (2026-10-29)

Drawn after the build. The change is a GPR edit, so the WT optimum is the same in both models; the difference shows
only when a core subunit is knocked out. The figure uses the alpha subunit `YALI1F04433g` (ATP1-like,
consensus-essential) as the example; every one of the 12 AND members gives the same growth (sweep, `../REPORT.md` §5).

## Files

- **r171_ko_alpha_before_vs_after.png / .svg**: the same map twice, pFBA with `YALI1F04433g` knocked out.
  - Left, before (`d8d1bac7`): the single-gene OR branches g / k keep R171 open. Growth 0.3117 = WT, and every mapped
    flux equals WT (checked, `index.json` → `ko_before_equal_to_wt`).
  - Right, after (`ca43ae0c`): R171 closed by the 12-gene AND. Growth 0.0342 (KO/WT 0.1097). The respiratory chain
    runs faster than in WT (complex III R305 18.3 vs 9.0, complex IV R304 9.2 vs 4.5, O₂ uptake 12.6 vs 4.7) and the
    proton leak R1384 dissipates the gradient (155.7); substrate-level ATP rises (R642 4.67, R694 4.65, R741 4.55)
    and the ADP/ATP carrier R815 falls from 18.3 to 4.5.
- **r171_pfba_wt.png**: WT pFBA (identical in both models).
- `flux.json` holds the LP results; `index.json` the map counts and tooltips.

## Model, condition, solver

- **Models:** `d8d1bac78138…` (before) and `ca43ae0c12e0…` (after), read from byte-identical short copies (SHA
  checked, `../logs/sweep_short_copies_sha256.txt`), because cobra rejects str paths of 260+ characters on macOS.
- **Condition:** Task A medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) + PO1f profile (`35307853…`).
- **Solver:** Gurobi 13, FeasibilityTol / OptimalityTol 1e-9 (`checked_solve.py`), real HOME. Each scenario: a
  checked solve (optimal, residual and bound checks, knocked-reaction flux ≤ 1e-9) then pFBA (optimal). WT
  reproduces 0.3116608216 in both models (the script stops otherwise). 12 LPs.

## Which reactions are drawn

- 17 fixed: the respiratory chain (R2062, R573, R570, R740, R1977, R305, R304), R171, the ADP/ATP carriers R815 /
  R816, the phosphate carriers R914 / R915 / R1702, substrate-level ATP R642 / R694 / R741 and ADH R130.
- 12 chosen by `flux.py`, not by hand: the reactions whose pFBA flux differs most between the two knockouts. They are
  the proton leak R1384, O₂, CO₂ and water exchanges / transports, and R659 / r0210 / R_CYST_II.
- **Observation, not traced:** in WT, pFBA carries 17.8–18.0 through phosphoserine phosphatase R659, r0210 and
  R_CYST_II (serine / cysteine / hydrogen sulphide), falling to about 0 in the knockout after the change. This is the
  same in both models at WT, so it is not caused by this change; why pFBA routes 18 units through it was not
  investigated (possible follow-up).

## Build

Run from `scripts/` with `/opt/miniconda3/bin/python`; `flux.py` from `platform/` with
`PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts`.

1. `flux.py MEDIUM PROFILE ../flux.json BASELINE.xml CANDIDATE.xml` (`logs/flux.log`)
2. `build.py ../flux.json WORKDIR` (no LP; `logs/build.log`)
3. `render.py WORKDIR 2400` (`logs/render.log`)

- `WORKDIR` is `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r171_atp_synthase_gpr_20261029/escher/` (maps, data,
  renders). `esc/node_modules` links to the cofactor-terms install (Escher 1.8.1, `escher.min.js` SHA-256
  `38c105f0…`); node v25.8.2; PNG by headless Chrome.
- From the UM11 figure: `maplib.py`, `escher_kernel.py`, `render_jobs.js` unchanged; `render.py` edited (map name,
  job file, compared pair, no status scenario); `flux.py`, `build.py`, `palette.py` rewritten.

## Limits

- The layout is automatic and carries no biological meaning; hub metabolites are repeated per reaction.
- No interactive HTML viewer was built for this figure; the tooltips in `index.json` are not drawn.
- pFBA is one optimum; it does not show that the knockout's routes are unique.
- The PNGs were inspected by eye (the compare image); the SVGs were not checked in a browser.
