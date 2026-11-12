# Escher figure: what closing R1172 did (2026-11-12)

Computation only; drawn after the build, from the two model files (standing rule: every finished fix gets a
before/after map). Scripts copied from the R2121 closure figure (`model/reports/r2121_malonyl_coa_20261110/escher/`);
only the reaction sections, expected WTs, FVA set, titles and map name change.

- **r1172_close_before_after.png:** the same map twice, pFBA flux in 2f352037 (R2121 closed) on the left and in the
  R1172-closed candidate on the right.
- **r1172_close_escher.html:** self-contained interactive viewer (Escher 1.8.1 inlined; no internet needed; open it
  in Chrome, Safari or Firefox). Scenarios `before`, `after` and `status`; `?scenario=…` also works. Hover or click a
  label for both models' pFBA flux, bounds and, for the FVA set, FVA at 99.9 % of the optimum.
- Full-size SVG/PNG renders of every scenario:
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r1172_hmgcoa_transport_20261112/stage2/escher_work/renders/`.

## Models, condition, solver

| | File | SHA-256 |
| --- | --- | --- |
| before | `model/candidates/E5_…_mcnaugh_rxns_v3_r2121_close_20261112.xml` | `2f35203750d1904aa4e50dd7d8eac494d5e6a15b142ea005c33152cea889af7c` |
| after | `model/candidates/E5_…_mcnaugh_rxns_v3_r2121_close_r1172_close_20261112.xml` | `684f3752764fb8f25219008966953863141724f34a451dffb77c235aa3cb2f62` |

- Both were read from byte-identical short-path copies (cobra on macOS rejects str SBML paths of 260+ characters).
- Task A medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) plus the PO1f profile
  `po1f_sd_leu.json` (`35307853…`). Gurobi 13 at FeasibilityTol / OptimalityTol 1e-9 (`checked_solve.py`).
- The script stops unless both WTs reproduce the recorded values, 0.3120526953810365 and 0.3111005359259494; both
  reproduced. 60 LPs ([logs/escher_flux.log](logs/escher_flux.log)).

## What is drawn (24 reactions, three sections)

| Section | Reactions |
| --- | --- |
| Mitochondrion | thiolase R83, HMG-CoA synthase R412, lyase R1974, succinyl-CoA:acetoacetate CoA-transferase R1933, leucine-route end R1968 / R1972 / R1973, CoA import R834, HMG-CoA transport R1172 (closed after) |
| Cytosol, mevalonate start | thiolase R82, HMG-CoA synthase R411, HMGR R410, mevalonate kinase R553, acetoacetate transport R1996 |
| Cytosolic acetyl-CoA supply | Asp-derived semialdehyde R1936 / R2120 / R2119, ATP-citrate lyase R1894, acetyl-CoA synthetase R90, acetylcarnitine R191 / R192 / R1032, acetyl-CoA carboxylase R88, GlcNAc-6-P synthase R561 |

Fluxes of 1e-9 or less are drawn and counted as 0. The layout is automatic and carries no biological meaning; hub
metabolites are repeated per reaction.

## Numbers on the map ([index.json](index.json), [flux.json](flux.json))

| | before (2f352037) | after (R1172 closed) |
| --- | ---: | ---: |
| Mapped reactions with pFBA flux (\|v\| above 1e-9) | 12 | 10 |
| Whole-model total \|pFBA flux\| | 282.16 | 284.50 |

- **Flux only before (3):** R83, R412, R834.
- **Flux only after (2):** R82, R411.
- **Flux in both (8):** R410, R553, R1936, R2120, R2119, R88, R90, R561.
- **No flux in either (10):** R1933, R1968, R1972, R1973, R1974, R1996, R1894, R191, R192, R1032.
- **Closed after:** R1172.

FVA at 99.9 % (selected): before, R1172 [−0.074, −0.0485] and R834 [0.0485, 0.074], the CoA return matching the
HMG-CoA export; R82 and R411 [0, 0.0236]. After, R1172 and R834 [0, 0]; R82, R411 and R410 [0.0718, 0.0736];
R2119 [0.2115, 0.5647] (before [0, 0.3502]).

## Checks and limits

- **Viewer checks:** jsdom draws all 24 reactions in every scenario
  ([logs/check_viewer_jsdom.log](logs/check_viewer_jsdom.log)). Headless Chrome loads Escher, draws 24 reaction groups
  and 159 nodes per scenario, and shows the expected labels, for example R1172 at −0.0721 before and 0 after, R411 at
  0 before and 0.0719 after ([logs/check_viewer_chrome.json](logs/check_viewer_chrome.json)). Tooltips were not
  exercised.
- **One optimum per model:** pFBA shows a single optimum. FVA at 99.9 % covers only the 13-reaction FVA set.
- **Python:** the map build needs networkx, which the root `.venv` lacks, so it ran with `/opt/miniconda3/bin/python`
  (Python 3.13.5). The flux script ran in the root `.venv`.
