# Escher figure: the five demand families of the 29 false negatives on e1e7c868 (2026-11-02)

Computation only. **Nothing in the model changed**, and none of the what-ifs was implemented. The map shows the family
reactions of `e1e7c868` as built, next to each family's recommended variant as an in-memory what-if. Each what-if is
applied **alone**, as in the scoping work; the combined model was never run. The what-if variants come from the scouts'
own builders ([scripts/whatif.py](scripts/whatif.py) names each source file and SHA-256); none was redesigned here.

## Files

- **fn29_families_status.png**: the as-built status of the 138 mapped reactions.
  - Orange: blocked on the condition (FVA at 0 % gives 0).
  - Dark grey: closed (0, 0), for example the Golgi AUR1 set R1825-R1834.
  - Light blue: can carry flux, but is 0 at the pFBA optimum.
  - Green: carries pFBA flux (only R289 / R287, the mannan lump, and the mitochondrial Met / Ile exits R1251 / R1218).
  - Purple: a what-if reaction that is not in the model.
- **fn29_families_pfba_as_built.png**: pFBA as built. Every family reaction is 0 apart from the four above.
- **fn29_families_pfba_whatif.png / .svg**: pFBA with each section under its own family's what-if. The demanded route
  carries flux in every section.
- **fn29_families_escher.html**: the self-contained interactive viewer (Escher 1.8.1 inlined, no internet needed).
  - Scenarios: `status`, `pfba_as_built`, `capacity_as_built`, `capacity_whatif`, `pfba_whatif`,
    `pfba_whatif_nglycan_B`. The URL parameter `?scenario=…` also works.
  - Hover or click a label for the GPR, bounds, pFBA and capacity as built and under the family's what-if.
- `flux.json` holds the LP results; `index.json` holds the map counts and tooltips.
- In the research root only (`$IYALI26_RESEARCH_ROOT/artifacts/tasks/fn29_demand_families_20261102/escher/`): the
  `capacity_*` and `pfba_whatif_nglycan_B` renders, the side-by-side `renders/compare.png / .svg`
  (pfba_as_built | pfba_whatif, 1.1 MB each), the maps, the scenario data and the viewer screenshots.

## Model, condition, solver

- **Model:** `e1e7c868` = `model/candidates/E5_…_aco_split_20261029.xml`, SHA-256
  `e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897`. It was read from the archived scout copy
  `fn29/e1e7c868.xml` (same SHA, checked in this run against both files).
- **Condition:** the scouts' condition file `fn29/cond_taskA_po1f.json` (SHA-256 `db78d705…3a63`): Task A medium and
  R612 (0, 0). It omits the PO1f plasmid LEU2 complement of R45 (REPORT.md, GPI corrections). It was loaded exactly as
  the scouts loaded it (`model.medium = …`, then the bounds). `xMAINTENANCE` keeps its model lower bound 7.8625.
- **Solver:** Gurobi 13.0.2 through cobra 0.30.0, `model.tolerance` 1e-9, real HOME, `/opt/miniconda3/bin/python`.
- **LPs: 581**, under the 1,500 budget, which the script asserts.
  - Per variant: WT 1 + pFBA 2 + FVA (1 + 2 per mapped reaction) + energy check 1.
  - As built: all 118 as-built reactions, 241 LPs. What-ifs: 165 (sphingolipid), 23 (GPI), 19 (cardiolipin), 49
    (N-glycan C2), 39 (N-glycan B), 45 (mitochondrial translation).
- **Stop checks:** both SHA-256 values match, and every variant reproduces its archived WT within 1e-9.

| Variant | WT here | Archived value it reproduces (file) |
| --- | ---: | --- |
| as built | 0.311660821633609 | 0.3116608216336117 (`sphingolipid/verify_lp/verify2.json`) |
| sphingolipid V2 | 0.311460090248134 | 0.3114600902481364 (same file) |
| GPI b1 | 0.311523212282885 | 0.3115232122828883 (`gpi/s2_variants.json`, scout) |
| cardiolipin a | 0.311059973357186 | `cardiolipin/verify_lp/verify_cardiolipin.json` |
| N-glycan C2 | 0.276393058443159 | `nglycan/verify_lp/verify_c2.json` |
| N-glycan B | 0.280770154862780 | same file |
| mitochondrial translation | 0.310670122842183 | `mito_translation/verify_lp/verify2.json` |

Energy check in every variant (every exchange (0, 1000), `xMAINTENANCE` lb 0, maximise `R_NTP1`): optimal, 0.0.

## Scenario definitions (the what-ifs, from the scouts' builders)

| Section | What-if | Edits (in memory) |
| --- | --- | --- |
| Sphingolipid | V2 | R1825-R1834 (0, 0) → (0, 1000); 10 gene-less links from the duplicate Golgi IPC species m1707-m1716 to m1263-m1272; `L_DAG_gm_em` m1704[gm] → m1649[em]; `biomass_C` += −0.000417 IPC-C m1271[go] |
| GPI anchor | b1 | `biomass_C` += −1.6962e-4 m603[er] (Yin 2007 proxy); `R_DM_GPI` left open |
| Cardiolipin | a | `biomass_C` += −7.1616e-4 cardiolipin m1651[mm] (0.08 × 12.98 mg/gDW / 1449.958 g/mol) |
| N-glycan | C2 = A2 + B (scenario `*_whatif`) | A2: `R_LLO_link`, `R_LLO_man9`, `R_LLO_glc3`, `R_OST`, `R_OST_M9`, `biomass_C` += −0.01 N-glycan (no source). B: R289 rewritten 1:1 to O-mannosyl protein, `R_OMAN1_t`, `R_MANNAN_GO`, `R_PIt_go_cy`, R1137 closed |
| N-glycan | B (scenario `pfba_whatif_nglycan_B`) | B alone, the part REPORT.md section 4 treats as decidable |
| Mitochondrial translation | a2 + R1251 + R1218 | `MITO_TRANSLATION_14` (14 charged tRNAs, AJ307410 fractions renormalised); `biomass_C` += −0.01863252 mtprot (no source); R1251 and R1218 (0, 1000) → (−1000, 1000) |

The changed reactions (R1825-R1834, R289, R1137, R1251, R1218) are drawn with their **as-built** equation. Their
tooltips say so. In the what-if, R289 makes the new O-mannosyl species, not the drawn Man4 + H⁺.

## Numbers on the map (`index.json`)

| Section | Mapped as built | Blocked | Closed | Can carry flux, pFBA 0 | pFBA ≠ 0 | What-if-only reactions | pFBA ≠ 0 under the what-if |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Sphingolipid | 69 | 18 | 20 | 31 | 0 | 11 | 23 |
| GPI anchor | 9 | 0 | 0 | 9 | 0 | 0 | 8 |
| Cardiolipin | 7 | 2 | 0 | 5 | 0 | 0 | 5 |
| N-glycan | 14 | 11 | 0 | 1 | 2 | 8 | 16 |
| Mitochondrial translation | 19 | 15 | 0 | 2 | 2 | 1 | 18 |

- These agree with REPORT.md:
  - V2 route flux: R1828 = `L_DAG_gm_em` = 1.2988e-4. The sum over R1825-R1834 matches the report's
    1.298789e-4; the other nine are 0.
  - GPI chain 5.284e-5 on R_GPIGNT to R376, R_DM_GPI 0.
  - Cardiolipin R1791 = R197 = 4.4554e-4 and R1836 = R1722 = R1721 = 2.2277e-4.
  - C2: R850 = R602 = 0.19556 and R289 = 0.065187.
  - Mitochondrial translation flux 5.789e-3.
- **Metabolite nodes:** dark red, large = the species a what-if adds to `biomass_C` (m1271[go], m603[er], m1651[mm],
  the N-glycan and mtprot pseudo-species). Purple = species that exist only in a what-if. Orange = every mapped
  reaction of the species is blocked or closed as built (judged on the mapped reactions only, not the whole model).
  Grey = the rest.

## Verification and limits

- **Checked 2026-11-02:**
  - the SHA-256 values and the seven WT values;
  - all six scenarios draw 138 reactions and 909 metabolite circles in jsdom (`logs/check_viewer_jsdom.log`) and in
    headless Chrome (`logs/check_viewer_chrome.json`);
  - Chrome shows no Escher-not-loaded warning;
  - R1836 reads 0.00 (pFBA, as built) and 0.000223 (what-if);
  - R850 reads 0.196 (C2) and 0.199 (B);
  - the status render was inspected.
- **Not checked:** the tooltip extension (there is no hover in headless runs), and the viewer in an interactive
  browser.
- **Offline use:** the viewer has no external `<script src>`. The inlined Escher bundle fetches a gif.js worker from
  cdn.jsdelivr.net only if the user exports a GIF.
- **Older base and condition:**
  - This is `e1e7c868` on the scouts' condition file. Newer candidates exist (REPORT.md, "Status on the current
    main").
  - `2bc263f2` already demands a ceramide pool, so its sphingolipid section would look different.
  - The condition file lacks the PO1f LEU2 complement.
- **What capacity means:** capacity is the extreme of the flux space, not a prediction.
- **Layout:** the layout is automatic (`scripts/maplib.py`) and carries no biological meaning. Hub metabolites are
  repeated per reaction.

## Build

Run with `/opt/miniconda3/bin/python` and `IYALI26_RESEARCH_ROOT` set. `WORKDIR` is
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/fn29_demand_families_20261102/escher/`. Its `esc/node_modules` links to the
cofactor-terms install (Escher 1.8.1, `escher.min.js` SHA-256 `38c105f0…9927`, jsdom).

1. `python scripts/flux.py flux.json` (581 LPs; reads the archive, writes nothing else)
2. `python scripts/build.py flux.json WORKDIR` (no LP; `maxw=16000`)
3. `python scripts/render.py WORKDIR 3000`
4. `python scripts/build_html.py WORKDIR fn29_families_escher.html`
5. `ESC=WORKDIR/esc node scripts/check_viewer.js VIEWER`, then `python scripts/check_viewer_chrome.py VIEWER WORKDIR OUT`

- **Copied from the UM11 figure:**
  - rewritten: `flux.py`, `build.py` and `palette.py`;
  - edited (names, scenarios, labels, legend): `render.py`, `build_html.py`, `check_viewer_chrome.py`;
  - unchanged: `maplib.py`, `escher_kernel.py`, `render_jobs.js`, `check_viewer.js`.
- **New:** `whatif.py`. It imports `cardiolipin/variants.py` and `nglycan/variants2.py` from the archive unchanged, and
  copies the V2, b1 and a2 builders verbatim. Those three scout scripts run LPs on import and hard-code scratchpad
  paths, so they cannot be imported.
