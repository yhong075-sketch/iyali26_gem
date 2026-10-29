# Escher figure: UM11, the nuclear RNR copies kept blocked by decision B (2026-10-29)

Computation only, drawn after the user's decision. **Nothing in the model changed** (Option B), so this is not
a before/after of a build. It shows the baseline as built next to the in-memory Option A what-if (both
thioredoxin forms linked nucleus ⇌ cytosol), which was **not** implemented. The point is that the map records
the reason the module stays blocked: opening it adds capacity but never flux at the optimum. The two census
roots, TRX1 `m931[C_nu]` (RNP, no producer) and TRX1-disulphide `m933[C_nu]` (RNC, no consumer), are drawn as
large dark-red nodes. Root marking was added on 2026-10-29 at the census session's request; the LPs were not
re-run.

## Files

- **um11_nuclear_rnr_asbuilt_vs_whatif.png / .svg**: the same map twice.
  - Left: *capacity as built*. All 11 UM11 reactions are 0.
  - Right: *capacity under what-if A*. All 11 gain capacity, and the cytosolic copies are unchanged.
- **um11_nuclear_rnr_status.png**: the as-built status on Task A.
  - Orange: blocked (the 11 UM11 reactions).
  - Light blue: can carry flux.
  - Purple: the two what-if links, which are not in the model.
- **um11_nuclear_rnr_pfba_as_built.png**: pFBA as built. The nuclear copies are 0, and dNDPs come from R702, R710
  and R1913 in the cytosol at this optimum.
- **um11_nuclear_rnr_pfba_whatif_A.png**: pFBA under what-if A. Every nuclear copy is 0, and dNDPs come from
  R702, R706, R710 and R1913 in the cytosol.
- **um11_nuclear_rnr_escher.html**: the self-contained interactive viewer (Escher 1.8.1 inlined, no internet
  needed).
  - Scenarios: `status`, `pfba_as_built`, `capacity_as_built`, `capacity_whatif_A`, `pfba_whatif_A`. The URL
    parameter `?scenario=…` also works.
  - Hover or click a label for the GPR, bounds, pFBA, capacity and no-uptake range in both variants.
- `flux.json` holds the LP results; `index.json` holds the map counts and tooltips.

## Model, condition, solver

- **Model:** `model/candidates/E5_…_dolichol_route_v2_um7_nuclear_diffusion_r2171_1to1_20261029.xml`, SHA-256
  `14f73bda4b0772019e82ea18a27e198c5ba175e67fef205a51dbfb090809e757`. It was read from a byte-identical short
  copy, because cobra rejects str paths of 260+ characters on macOS.
- **Condition:** Task A medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) plus the
  PO1f profile `$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json` (`35307853…`). The script stops
  unless the uracil exchange R1354 is open to uptake 1000 and R612 is (0, 0); both were as required.
- **Solver:** Gurobi 13, FeasibilityTol and OptimalityTol 1e-9 (`checked_solve.py`), real HOME.
  - WT is 0.3116608216 in both variants (checked; the script stops on a mismatch beyond 1e-9).
  - pFBA is optimal in each variant.
  - FVA (cobra) runs at 0 % of the optimum, for the mapped reactions and for the whole model.
  - The no-uptake FVA sets every exchange to (0, 1000), relaxes every positive lower bound (including
    `xMAINTENANCE` 7.8625) and fixes biomass at 0.
  - 9,476 LPs in total.
- **The what-if** is `trx+trxss` from `../scripts/um11_apply.py`, the same in-memory variant as Stage 1.

## Numbers on the map (`index.json`)

| On the 25 mapped reactions | as built | what-if A |
| --- | ---: | ---: |
| UM11 reactions blocked on Task A | 11 | 0 |
| Capacity > 0 on Task A | 12 | 25 |
| pFBA flux ≠ 0 | 5 (R702, R710, R1913, R755, R781) | 6 (+ R706, an alternate optimum to R702) |
| Non-zero with no uptake | 2 (R789, R1104) | 2 (the same) |

- **Whole-model blocked set on Task A:** 828 as built, 817 under A. The 11 freed reactions are exactly the UM11
  set, and none is newly blocked.
- **Metabolite nodes:** 37 distinct species, 156 drawn circles.
  - 2 are the census roots: TRX1 and TRX1-disulphide in the nucleus, dark red, drawn 10 times because they are
    repeated per reaction.
  - 7 more nuclear species are dead as built and alive under A: UMP, UDP, CDP, GDP, dADP, dCDP and dGDP.
  - The other 28 have an unblocked reaction, and none is dead in both variants.
- **R789 / R1104 no-uptake loop (as built, unchanged by A).** The nuclear dUMP kinase R789, the dUDP and dUMP
  diffusions R1104 / R1105, the cytosolic kinase R788 and the nuclear ATP / ADP diffusions can cycle at ±1000
  with no uptake. This is a net-zero futile loop; the energy check gives no free ATP. It is outside UM11 and
  noted as a follow-up in the report.

## Verification and limits

- **Checked 2026-10-29:**
  - the model SHA, WT growth in each variant, uracil and R612;
  - all five scenarios draw 25 reactions and 156 metabolite circles, in jsdom (`logs/check_viewer_jsdom.log`)
    and in headless Chrome (`logs/check_viewer_chrome.json`);
  - after the root marking: 4 node colours in jsdom, and 10 dark-red root circles in the headless-Chrome DOM
    (`logs/check_viewer_chrome_roots.log`);
  - Chrome shows no Escher-not-loaded warning. R707 reads 0.00 as built and 1.58 under A, with pFBA 0 in both;
    R706 reads 1.58 and T_TRX1tn reads 1.63 under A;
  - the compare and pFBA renders were inspected.
- **Not checked:** the tooltip extension (there is no hover in headless runs), and the viewer in an interactive
  browser.
- **Offline use:** the viewer has no external `<script src>`. The inlined Escher bundle fetches a gif.js worker
  from cdn.jsdelivr.net only if the user exports a GIF (found by the audit).
- **Layout and drawing:**
  - The layout is automatic (`scripts/maplib.py`, from the gap census) and carries no biological meaning. Hub
    metabolites, and nuclear TRX1 / its disulphide (used by five reactions), are repeated per reaction.
  - The what-if links are drawn in the as-built scenarios too, with no data ("(nd)").
  - Background section boxes appear only in the PNG / SVG renders.
- **What capacity means:** capacity is the extreme of the flux space, not a prediction.

## Build

Run from `scripts/`. `flux.py` runs from `platform/` with
`PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts:../model/reports/um11_nuclear_rnr_20261029/scripts`,
using miniconda Python with cobra and Gurobi.

1. `flux.py MEDIUM PROFILE ../flux.json MODEL.xml`
2. `build.py ../flux.json WORKDIR` (no LP)
3. `render.py WORKDIR 2400`
4. `build_html.py WORKDIR ../um11_nuclear_rnr_escher.html`
5. `ESC=WORKDIR/esc node check_viewer.js VIEWER`, then `check_viewer_chrome.py VIEWER WORKDIR OUT`

- `WORKDIR` is `$IYALI26_RESEARCH_ROOT/artifacts/tasks/um11_nuclear_rnr_20261029/escher/` (maps, data, renders,
  viewer screenshots). `esc/node_modules` links to the cofactor-terms install (Escher 1.8.1, `escher.min.js`
  SHA-256 `38c105f0…`).
- **Copied from the UM2 figure:**
  - rewritten: `flux.py` (variants from `um11_apply.py`; what-if-only reactions), `build.py` (sections, status,
    scenarios) and `palette.py` (pFBA scale, purple what-if status);
  - edited (names, scenarios, labels): `render.py`, `build_html.py`, `check_viewer_chrome.py`;
  - unchanged: `maplib.py`, `escher_kernel.py`, `render_jobs.js`, `check_viewer.js`.
