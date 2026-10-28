# Escher figure: what the UM7 nuclear-diffusion change did (2026-10-28)

Computation only; drawn after the build, from the two model files. The user's brief asks every finished gap
session to draw an Escher figure of what changed.

- **um7_nuclear_before_after.png / .svg**: the same map twice, *capacity before* on the left and *capacity after*
  on the right.
- **um7_nuclear_status.png**: the same map coloured by Task A blocked status (freed / blocked in both / closed /
  new / unblocked in both).
- **um7_nuclear_escher.html**: the self-contained interactive viewer (Escher 1.8.1 inlined, no internet needed;
  open it in Chrome, Safari or Firefox). Scenario drop-down: `capacity_before`, `capacity_after`, `before`,
  `after`, `status`; the URL parameter `?scenario=…` also works. Hover or click a label for pFBA, FVA at the
  optimum, capacity (Task A and open), bounds and blocked status of both models.

## Why capacity and not only pFBA

pFBA carries **no flux** through any blocked UM7 reaction or new link at the optimum in either model (the already-active R568 / R580 carry 0.0010 in both): the change widens
the flux space and does not move the optimum (WT identical to 4.6e-15 relative). A pFBA-only figure would show
the same grey map twice. The figure therefore leads with **capacity** = max |v| over FVA at 0 % of the optimum
on Task A + PO1f (how much flux any feasible state can put through the reaction). The pFBA scenarios `before` /
`after` are in the viewer for completeness.

## Models, condition, solver

| | File | SHA-256 |
| --- | --- | --- |
| before | `model/candidates/E5_…_r1726_cpt_dolichol_route_v2_20261028.xml` | `9d655c2bdbf82721a38f148a9a4ee0a385bb509ea7b591111107113c4aa8335c` |
| after | `model/candidates/E5_…_dolichol_route_v2_um7_nuclear_diffusion_20261028.xml` | `5f54eb22086b06b2bd5e8715b89984060cd43c1a7fa21e33e5581d867ee54179` |

- Both models were read from byte-identical short copies in the session scratchpad. cobra rejects str SBML paths
  of 260+ characters because its Windows test (`"win" in platform`) also matches `darwin`. The loader reports
  the canonical SHA-256 of each copy (`flux.json`).
- **Condition:** Task A medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) + PO1f
  profile `$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json` (`35307853…`). The script checks that
  uracil exchange R1354 is at (−1000, 1000) and R612 at (0, 0), and stops otherwise.
- **Solver:** Gurobi 13 at FeasibilityTol / OptimalityTol 1e-9 (`checked_solve.py`), real HOME. Both WTs
  0.3116608216 (checked, the script stops on a mismatch). The FVA is cobra 0.30.0, run at 1.0, at 0 on Task A,
  and at 0 with every exchange open. 272 LPs.
- **Blocked status and node colours** come from `../built_checks.json`: whole-model FVA at 0 % of every reaction,
  |v| ≤ 1e-9, both models. A node is coloured by whether **every** reaction of that metabolite is blocked on
  Task A. This is a stricter, simpler criterion than the Ponce-de-León gap classes and is not the census
  definition.

## What is drawn

23 reactions in three sections:

| Section | Reactions |
| --- | --- |
| Nucleus | R569 (NAD synthase, QNS1 ortholog), R579 (NMNAT, NMA1 ortholog), R567 (sirtuin step), the inherited transports R1082, R1270, R1267, R1017, R1212, R1214, and the **new R_ADPRtn, R_PPItn** |
| Cytosol | R580, R568, R578 (NAD synthesis), R1912 (ADP-ribose pyrophosphatase, YSA1 ortholog), R419 (inorganic pyrophosphatase), R576 (nicotinamidase) |
| Greyed, unchanged | R396 (homocitrate synthase), R965, R1001, R1055, R1000, R92 (closed) |

- Layout is automatic (the gap-census builder, `scripts/maplib.py`) and carries no biological meaning. Hub
  metabolites (ATP, H2O, H+) are repeated per reaction.
- The new reactions are absent from `9d655c2b`, so the *before* scenarios show them as "(nd)".

## Numbers on the map (`index.json`)

| On the 23 mapped reactions | before (9d655c2b) | after (5f54eb22) |
| --- | ---: | ---: |
| In the model | 21 | 23 |
| Blocked, Task A (whole-model FVA) | 18 | 9 |
| Blocked, every exchange open | 16 | 6 |
| Can carry flux on Task A (capacity > 0) | 3 | 14 |
| Carry pFBA flux at the optimum | 2 | 2 |

- **Freed on Task A:** R567, R569, R1017, R1082, R1212, R1214, R1267, R1912, R576.
- **New and unblocked:** R_ADPRtn, R_PPItn.
- **Blocked in both:**
  - the homocitrate branch R396, R965, R1001, R1055, R1000 (and R92, closed);
  - R579, R1270 and cytosolic R578: Task A has no NMN source. They are freed with every exchange open.
- **Capacity after:** 10.7 for the sirtuin route (R567, R1912, both links, R569, R1082 …); R419 55.5 in both.
- **Metabolite nodes:** 10 change from "every reaction blocked" to "has an unblocked reaction". These are ADP-ribose
  [cy] and [nu], PPi, NAD, nicotinamide, AMP, Gln and Glu [nu], deamido-NAD [nu] and nicotinamide [cy]. 6 stay
  fully blocked: acetate, CoA, acetyl-CoA, 2-oxoglutarate, homocitrate and NMN [nu]. None becomes blocked.

## Verification and limits

- **Checked this session:**
  - both SHAs, WT growth, uracil and R612;
  - all five scenarios draw 23 reactions and 125 metabolite nodes, in jsdom (`logs/check_viewer_jsdom.log`) and in
    headless Chrome (`logs/check_viewer_chrome.json`). Chrome shows no Escher-not-loaded warning, and the label
    R_ADPRtn reads "(nd)" in `capacity_before` and "10.7" in `capacity_after`;
  - the compare render was inspected.
- **Not checked:** the tooltip extension (no hover in headless runs); the viewer in an interactive browser.
- The background section boxes are drawn only in the PNG / SVG renders, not in the viewer.
- Capacity is the extreme of the flux space and not a prediction; pFBA shows one optimum.

## Build

Run from `scripts/` (flux.py from `platform/`):

1. `flux.py MEDIUM PROFILE ../flux.json BEFORE.xml AFTER.xml` (miniconda Python with cobra + Gurobi;
   `PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts`).
2. `build.py ../flux.json ../../built_checks.json WORKDIR` (no LP).
3. `render.py WORKDIR 2400`: SVG via jsdom + Escher, PNG via headless Chrome.
4. `build_html.py WORKDIR ../um7_nuclear_escher.html`.
5. `ESC=WORKDIR/esc node check_viewer.js VIEWER` and `check_viewer_chrome.py VIEWER WORKDIR OUT`.

- Copied from the dolichol-route figure. The following changed:
  - `flux.py`: reaction set, 0 % FVA added;
  - `build.py`: rewritten for these sections and status source;
  - `palette.py`: capacity and status scales;
  - `render.py`, `build_html.py`, `check_viewer_chrome.py`: names, scenarios, labels.

  `maplib.py`, `escher_kernel.py`, `render_jobs.js` and `check_viewer.js` are unchanged.
- **Work folder:** `$IYALI26_RESEARCH_ROOT/artifacts/tasks/um7_nuclear_transports_20261028/escher/` holds `maps/`,
  `data/`, `index.json`, `renders/` (each scenario separately) and the Chrome screenshots. `esc/node_modules` is
  a link to the cofactor-terms install (Escher 1.8.1, `escher.min.js` SHA-256 `38c105f0…`).
