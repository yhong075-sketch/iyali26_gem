# Escher figure: the UM2 vacuolar closures, kept by decision A0 (2026-10-28)

Computation only; drawn after the user's decision. **Nothing in the model changed** (Option A0), so this is not a
before/after of a build. It shows the newest candidate as built, next to in-memory what-if tests that were
**not** implemented. The aim is that the map itself records the reason the vacuole stays closed.

## Files

- **um2_vacuole_asbuilt_vs_whatif.png / .svg**: the same map twice. Left: *capacity as built*, all 0. Right:
  *what-if C′*, with the pump, the water, R1163 as a vacuole → cytosol leak and the 11 AVT/ATG22-backed
  transports opened in memory. Capacity appears only on the futile shuttles.
- **um2_vacuole_status.png**: the as-built status on Task A. Black = closed (0, 0); orange = open but
  blocked.
- **um2_vacuole_loop_open17.png**: the census test "open all 17" with R1163 as written, flux with **no
  uptake at all**. Red = internal loop.
- **um2_vacuole_escher.html**: the self-contained interactive viewer (Escher 1.8.1 inlined, no internet
  needed; open it in Chrome, Safari or Firefox).
  - Scenarios: `status`, `capacity_as_built`, `capacity_whatif_cprime`, `loop_whatif_open17`. The URL
    parameter `?scenario=…` also works.
  - Hover or click a label for the GPR, bounds, pFBA, capacity and no-uptake range under all three
    variants.
- `flux.json` holds the LP results; `index.json` holds the map counts and tooltips.

## Model, condition, solver

- **Model:** `model/candidates/E5_…_dolichol_route_v2_um7_nuclear_diffusion_20261028.xml`, SHA-256
  `5f54eb22086b06b2bd5e8715b89984060cd43c1a7fa21e33e5581d867ee54179`. It is the newest on `main` (PR #30),
  read from a byte-identical short copy (cobra rejects str paths of 260+ characters on macOS). The UM2
  reactions are identical in `9d655c2b`, on which Stage 1 ran.
- **Condition:** Task A medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) +
  PO1f profile `$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json` (`35307853…`). The script
  stops unless the uracil exchange is open to uptake 1000 and R612 is (0, 0).
- **Solver:** Gurobi 13, FeasibilityTol and OptimalityTol 1e-9 (`checked_solve.py`), real HOME.
  - WT is 0.3116608216 in all three variants (checked; the script stops on a mismatch beyond 1e-9).
  - pFBA is optimal in each.
  - FVA is cobra, at 0 % of the optimum, for the mapped reactions and for the whole model.
  - The no-uptake FVA sets every exchange to (0, 1000), relaxes every positive lower bound (incl.
    `xMAINTENANCE` 7.8625) and fixes biomass at 0.
  - 14,388 LPs in total.
- **What-if bounds** are in `flux.json` (`variants_bounds`); they are the blocks of `../scripts/stage1_um2.py`.

## Numbers on the map (`index.json`)

| On the 39 mapped reactions (37 UM2 + open R2055/R2056) | as built | what-if C′ | what-if open 17 as written |
| --- | ---: | ---: | ---: |
| Closed (0, 0) | 17 | 3 (polyamines) | 0 |
| Capacity > 0 on Task A | 0 | 12 | 12 |
| pFBA flux at the optimum | 0 | 0 | 0 |
| Non-zero with no uptake (loop) | 0 | 0 | 11 (R1163 at 1000; 10 aa transports at 500) |

- **Whole-model blocked set on Task A:** 828 as built, 816 under C′. The 12-reaction difference is exactly
  the C′ capacity set: R795, R1163, R869, R870, R878, R879, R882, R883, R887, R888, R897, R898.
- **R2055 / R2056** (open His/Lys vacuolar import) are blocked as built: vacuolar His and Lys have no exit.
  They are drawn for context; a separate follow-up.
- **Metabolite nodes (Task A):**
  - 23 have every reaction blocked as built and in C′: the dipeptide species, vacuolar Asp, Glu, Gly, Ala,
    Pro, water, polyamines, His, Lys.
  - 6 are blocked as built but not in C′: vacuolar H+, Asn, Gln, Ile, Leu, Tyr.
  - 21 have an unblocked reaction.

## Verification and limits

- **Checked 2026-10-28:**
  - the model SHA, WT growth in each variant, uracil and R612;
  - all four scenarios draw 39 reactions and 222 metabolite circles, in jsdom
    (`logs/check_viewer_jsdom.log`) and in headless Chrome (`logs/check_viewer_chrome.json`);
  - Chrome shows no Escher-not-loaded warning. Labels: R795 reads 0.00 as built and 55.5 in C′; R1163 reads
    1.00e+3 in the loop scenario;
  - the compare and status renders were inspected.
- **Not checked:** the tooltip extension (no hover in headless runs), and the viewer in an interactive
  browser.
- **Layout and drawing:**
  - The layout is automatic (`scripts/maplib.py`, from the gap census) and carries no biological meaning.
    Hub metabolites (H+, ATP, H2O) are repeated per reaction.
  - Background section boxes appear only in the PNG / SVG renders.
- **What capacity means:** capacity is the extreme of the flux space, not a prediction. In C′ it is reached
  only by spending ATP in futile cycles.

## Build

Run from `scripts/`; `flux.py` runs from `platform/` with
`PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts`, using miniconda Python with cobra +
Gurobi.

1. `flux.py MEDIUM PROFILE ../flux.json MODEL.xml`
2. `build.py ../flux.json WORKDIR` (no LP)
3. `render.py WORKDIR 2400`
4. `build_html.py WORKDIR ../um2_vacuole_escher.html`
5. `ESC=WORKDIR/esc node check_viewer.js VIEWER`, then `check_viewer_chrome.py VIEWER WORKDIR OUT`

- `WORKDIR` = `$IYALI26_RESEARCH_ROOT/artifacts/tasks/um2_vacuolar_closures_20261028/escher/` (maps, data,
  renders of each scenario). `esc/node_modules` links to the cofactor-terms install (Escher 1.8.1,
  `escher.min.js` SHA-256 `38c105f0…`).
- **Copied from the UM7 figure:**
  - rewritten: `flux.py` (one model, in-memory variants, whole-model and no-uptake FVA added), `build.py`
    (sections, status, scenarios) and `palette.py` (loop scale);
  - edited (names, scenarios, labels): `render.py`, `build_html.py`, `check_viewer_chrome.py`;
  - unchanged: `maplib.py`, `escher_kernel.py`, `render_jobs.js`, `check_viewer.js`.
