# Escher figure: what closing R2121 did (2026-11-12)

Computation only; drawn after the build, from the two model files (standing rule: every finished fix gets a
before/after map).

- **r2121_close_before_after.png:** the same map twice, pFBA flux in e584151f on the left and in the R2121-closed
  candidate on the right.
- **r2121_close_escher.html:** self-contained interactive viewer (Escher 1.8.1 inlined; no internet needed; open it
  in Chrome, Safari or Firefox). Scenarios `before`, `after` and `status`; `?scenario=…` also works. Hover or click a
  label for both models' pFBA flux, bounds and, for the supply set, FVA at 99.9 % of the optimum.
- Full-size SVG/PNG renders of every scenario: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r2121_malonyl_coa_20261110/stage2/escher_work/renders/`.

## Models, condition, solver

| | File | SHA-256 |
| --- | --- | --- |
| before | `model/candidates/E5_…_mcnaugh_rxns_v3_20261103.xml` | `e584151fc79475fc2daa59e76c8c84ebfd06b0f3cb1fbd75fbb69cc9f31039d2` |
| after | `model/candidates/E5_…_mcnaugh_rxns_v3_r2121_close_20261112.xml` | `2f35203750d1904aa4e50dd7d8eac494d5e6a15b142ea005c33152cea889af7c` |

- Both were read from byte-identical short-path copies (cobra on macOS rejects str SBML paths of 260+ characters).
- Task A medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) plus the PO1f profile
  `po1f_sd_leu.json` (`35307853…`).
- Gurobi 13 at FeasibilityTol / OptimalityTol 1e-9 (`checked_solve.py`).
- The script stops unless both WTs reproduce the recorded values: 0.31248856100769146 and 0.3120526953810365, both
  reproduced.
- 64 LPs ([logs/escher_flux.log](logs/escher_flux.log)).

## What is drawn (51 reactions, three sections)

| Section | Reactions |
| --- | --- |
| Cytosol, malonyl-CoA supply | R2121 (closed after), ACC R88, biotin ligase R175, carbonic anhydrase R1026 / R2202, Asp-derived semialdehyde R1936 / R2120 / R2119, malonyl-CoA transport R1575, FAS R1392 / R1393, ATP-citrate lyase R1894, acetyl-CoA synthetase R90, GlcNAc-6-P synthase R561 |
| Peroxisome | oleate import and activation R1482 / R1459, β-oxidation to laurate (R1496 … R84), thioesterase R619, laurate export R1118, catalase R195, acetylcarnitine export R192 / R1032 / R191 |
| Lipid particle and ER membrane | TAG lipase R1715, laurate import and activation R1569 / R1445, elongation to palmitoyl-CoA (R1406 … R1428), palmitoyl-CoA transport R1578, DGAT160_lp |

Fluxes of 1e-9 or less are drawn and counted as 0. R2119 carries −1.7e-12 of solver noise in the e584151f pFBA;
the first version counted it at a 1e-12 threshold and reported 46 (independent audit, corrected 2026-11-12).

The sections were chosen from the reactions that switch on or off between the two pFBA solutions
(`stage1_checks.json` in the artifacts). The layout is automatic and carries no biological meaning; hub metabolites are
repeated per reaction.

## Numbers on the map ([index.json](index.json))

| | before (e584151f) | after (R2121 closed) |
| --- | ---: | ---: |
| Mapped reactions with pFBA flux (\|v\| above 1e-9) | 45 | 11 |
| Whole-model total \|pFBA flux\| | 294.02 | 282.16 |

- **Flux only before (36):** the peroxisomal β-oxidation chain, the ER elongation chain, the TAG lipase and the
  acetylcarnitine export.
- **Flux only after (3):** R88, R2119 and R1393.
- **Flux in both (8):** R1392, R1575, R1936, R2120, R2202, R561, R90, DGAT160_lp.
- **No flux in either (3):** R1026, R175, R1894.
- **Closed after:** R2121.

## Checks and limits

- **Viewer checks:**
  - jsdom draws all 51 reactions in every scenario ([logs/check_viewer_jsdom.log](logs/check_viewer_jsdom.log)).
  - Headless Chrome loads Escher, draws 51 reaction groups and 333 nodes per scenario, and shows the expected labels,
    for example R2121 at 0.186 before and 0 after, and R88 at 0 before and 0.00785 after
    ([logs/check_viewer_chrome.json](logs/check_viewer_chrome.json)).
  - Tooltips were not exercised.
- **One optimum per model:** pFBA shows a single optimum. FVA at 99.9 % covers only the 14-reaction supply set.
- **Python:** the map build needs networkx, which the root `.venv` lacks, so it ran with `/opt/miniconda3/bin/python`
  (Python 3.13.5). The flux script ran in the root `.venv`.
