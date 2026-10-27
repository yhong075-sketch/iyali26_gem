# Independent audit of cofactor_biomass_terms_20261027 (2026-10-27)

Fresh read-only reviewer (a subagent of this session with no conversation history), given TASK.md, REPORT.md,
the deliverables, logs, scripts, the diff and both model files. Interpreter: the main checkout's
`.venv/bin/python` (cobra 0.30.0, gurobipy 13). No file edited, no git state changed. Its scratch outputs are
in its scratchpad folder `audit/` (own cobra compare script, the yeast-GEM download, the closure-check rerun).

## Verdicts

| Claim | Verdict | Evidence used by the auditor |
| --- | --- | --- |
| C1 candidate SHA `00e8e397…`; control = a1237821; seed 1 = candidate; `.build.json` module hash = current file | **supported** | `shasum`, `cmp` on the four XMLs; both `.build.json` identical; module hash `846ea1e8…`; head `aeba460`; step `applied` from `a1237821_biomass_C` |
| C2 diff confined to `biomass_C` (six coefficients, seven notes) | **supported** | own cobra field-by-field compare: 8 differences = biomass_C stoichiometry, biomass_C notes, the six species' derived reaction lists; raw XML diff 13 added lines, 0 removed, all inside biomass_C; counts 2316/1882/1081 both |
| C3 coefficients = yeast-GEM 9.1.1 r_4598 | **supported (verified)** | own `curl` of the raw yml, SHA-256 `2ddfaaa7…`, version 9.1.1, the six cytosolic species' values |
| C4 fails closed; rerun `already_correct`; test covers it | **supported** | focused test 5 passed / 23 subtests; extra probes (charge, compartment, TDP coefficient, missing before-state id) all refused without edits; rerun on the built file `already_correct` |
| C5 knockout sweep and scores | **supported** | recount from ko_compare.tsv equals the summary at all four cutoffs; 11 class changes, identical at every cutoff, all to essential; 20 infeasible on both, 16 after mtDNA exclusion; spot re-solve of 7 knockouts × 2 models on Gurobi 1e-9 reproduces the calls |
| C6 census 756 → 722 (37 freed ⊂ 51), 3 GLPK-only calls refuted at 1e-9, no new gaps, Task A 887 → 852 | **supported** | census_compare JSONs, drain-test union recomputed (51; 37 ⊂ 51; the 14 as listed), FVA JSON (R1343/R1344/R970 max 1.47e-4 both), baseline rerun vs recorded: summary identical, blocked and gap columns 0 diffs |
| C7 phase-1 growth, routes, energy check | **supported** | phase1_checks.json / phase1_egc.json / solves.json values; R568/R580/R571 all `C_cy`; infeasible census method and relaxed 0.0 reproduced in memory |
| C8 census EGC writes NaN as 0.0 | **partially supported** | mechanism verified (line 29 of closure_check.py; Kennedy rerun with the census script reproduces the recorded JSON; status probe: `infeasible`, value nan, written 0.0); the status label "undefined" in the first draft of REPORT.md was not what the auditor observed |

## Findings and responses

1. TASK.md §5.2 / §5.4 promise a "stale SHA" guard and test; the SHA field is provenance only.
   *Response:* accepted. REPORT.md §2 and §7 state it; TASK.md is the approved plan and is left as written.
   The before-state comparison (stoichiometry, bounds, note keys, species identity and context) is the
   executable staleness guard, as in the heme/iron precedent.
2. GLPK status for the census energy check is `infeasible` in the auditor's run, `undefined` in this session's.
   *Response:* REPORT.md §8 now gives both; the NaN-to-0.0 mechanism is unaffected.
3. "Largest residual 3.2e-13" was the candidate sweep's value; the baseline sweep's is 1.17e-12.
   *Response:* REPORT.md §5 corrected to give both.
4. The baseline census rerun relabels the modules (same partition, different UM numbers) and this was undisclosed.
   *Response:* REPORT.md §6 now says so and that module names follow the recorded census.
5. The cli help text says the flag "runs last"; GPR-only CoQ steps follow it in main.py.
   *Response:* noted in REPORT.md §2; the help text is left unchanged so the built record's source hashes
   stay those of the reviewed code.
6. AUDIT.md and the STATE.md entry did not exist at audit time.
   *Response:* this file and the STATE.md entry written after the audit.
7. Census logs show Gurobi's restricted licence (GAP_HOME overrides HOME); the census runs on GLPK.
   *Response:* noted in REPORT.md §6; no effect on results.

## Acceptance criteria (TASK.md §8), as assessed by the auditor

1 phase 1 no stop: met. 2 pytest: met (before 264/11/31/18, after 269/11/31/18, FAILED/ERROR sets identical,
+5 tests / +23 subtests). 3 builds: met. 4 diff: met. 5 WT optimal, routes with mechanism: met. 6 class changes
explained, false positives reported: met. 7 census, none newly blocked at 1e-9, freed ⊂ 51: met. 8 audit
recorded and STATE.md entry: met after this file and the entry below it.
