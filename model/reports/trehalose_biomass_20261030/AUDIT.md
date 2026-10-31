# Independent audit: trehalose_biomass_20261030 (2026-10-30)

A fresh general-purpose reviewer with no conversation history, read-only (no edits, builds, git writes or LP
solves), given TASK.md, the code, data, test, models, logs and figures. It audited the first build `d2e4c0c4…`.
Its report is kept below as returned (lightly reformatted); the responses are in [REPORT.md](REPORT.md) §6. After the
audit the model was rebuilt as v2 `164542fb…` with the two reworded notes; a raw XML diff of v2 against `e1e7c868`
again shows 8 lines (7 added `trehalose_biomass_*` note lines, 1 removed `m294` species reference), checked by the
main session.

**Procedural note from the reviewer:** it wrote three scratch files in the session scratchpad (outside the repo:
`audit_flores_text.txt`, `audit_pytest_before.txt.ids`, `audit_pytest_after.txt.ids`) although asked not to create
files; no LP, git write or build.

## Verdicts

| # | Claim | Verdict | Evidence |
| --- | --- | --- | --- |
| 1 | Hashes and argv | supported | Baseline `e1e7c868…`, candidate `d2e4c0c4…` on disk; scratchpad control = `e1e7c868…`, seed-1 = `d2e4c0c4…`. Candidate argv = baseline argv + `--trehalose-biomass-removal` (only the `__main__`/output paths differ). Recorded `source_sha256` (cli.py, main.py, trehalose_biomass.py) and `data_sha256` (trehalose_biomass.json) equal the files then on disk. |
| 2 | Whole-model diff is exactly two differences | supported (checked straight from the XMLs) | Raw `diff` gives 8 lines: 7 added `trehalose_biomass_*` note lines, 1 removed `m294[C_cy]` speciesReference (0.00686). cobra read (no solve): objective `{biomass_C: 1.0}` max in both; `biomass_C` 62 → 61 metabolites, 31 → 38 note keys; m294 reactions 9 → 8; `m297` still only R143 → R767. JSON diff agrees. |
| 3 | The step is guarded as TASK claims | partially supported | Spec validated; objective `biomass_C` alone; m294 identity; full composition pinned (key set, `isclose` 1e-12, exact bounds); trehalose coefficient compared exactly; reaction context pinned; nonsensical term/notes combinations rejected; post-check; notes restored on context exit. cobra 0.30.0 `add_metabolites(combine=False)` with 0 writes 0 to the constraint, pops the entry and the reverse link, and records the undo: sound. Order in main.py: after aconitase; the only later steps (CoQ functional GPR, C5 GPR, literature revision) are not in the argv and edit no `biomass_C`. Note text has no `<` or `&`; all 7 notes survive in the XML. Gaps: (a) other `biomass_C` notes pinned by key only, not value; (b) objective direction unchecked; (c) coefficient exact vs composition `isclose` (harmless, inconsistent). |
| 4 | Test meaningful; the two development fixes sound | supported, small gaps | 12 stale cases each a distinct path with unchanged-state assertions; 6 spec rejections; apply / idempotent / SBML round trip; guard re-check on the applied model; context exit; flag. `re.escape`: `[C_cy]` is a regex class, so escaping is a real fix. `solver.update()`: optlang applies objective changes lazily, and the guard's own `objective_coefficient` read would otherwise show as a change; it does not hide edits. Gaps: no case for a changed value of another note (passes by design), objective direction, or a partial subset of this curation's notes; the `"Conflicting"` regex is loose. |
| 5 | pytest before vs after | supported | 409 → 414 passed, 19 failed / 20 errors in both with identical IDs (39 each), subtests 960 → 978 = 5 new tests with 12 + 6 subtests. The after run started after the last test edit. Progress-dot logs do not name the new tests; the counts are the evidence. |
| 6 | Sweep | supported | Recomputed from `ko_compare.tsv` (1,081 rows): only YALI1D17620g, YALI1E17616g, YALI1E36568g change (essential → non-essential at all four cutoffs; KO 0.0 → 0.312296365 = WT). WT 0.3116608216 → 0.3122963651. 15 %: 125/190/42/692 → 125/190/39/695; at 1/5/10 % FP −3, TN +3. No failing check point. The 20 infeasible KOs are the same set (11 positives, 5 negatives, 4 mtDNA). No other KO/WT ratio moves by more than 0.01; none within 0.002 of a cutoff. Medium `6eca70fd…`, profile `35307853…`; `checked_solve.py` tolerance 1e-9 and `model.solver.problem.reset()` before each solve; Gurobi academic licence line in the log. |
| 7 | Evidence claims (TASK §3–§5, notes) | mostly supported; a few overstatements | Re-extracted PMC3171402.xml (`c5529d32…`): medium and units; Table 4 (CJM 645 < 1; CJM 613 8.9 ± 2.1; CJM 683 18 ± 2; CJM 703 < 1; three cultures); "below 1 … grown in glucose up to stationary phase or in glycerol"; doubling times 149 ± 4 / 139 ± 4 (four experiments); no trehalose in nth1 tps1 and in heat-shocked tps1; tps3 abolished heat-shock accumulation; YlTPS3 does not complement Sc tps1; TPS2 by Génolevures annotation only; heat shock 7–20; CJM 645 = MATB met6 from A1-5 × CJM 364. UniProt O74932: 469 aa, YALI0E14685g, EMBL AJ011032 (mRNA), cites PMID 21931609; locus map YALI1E17616g → YALI0E14685g. Origin: `generate_iYali5GEM.m` line 56; iYali4 `R_xBIOMASS` 0.0032 / 0.4068 / 0.4415 / 0.1104; ratios 2.144, 2.135, 2.137, 2.136; `iyali26.xml` line 34795. |
| 8 | Justified by composition and phenotype; limits stated | supported | Rests on the Table 4 bound, the tps1 phenotype and the unsourced value; limits (CJM 645, YNB-glutamate, no re-normalisation, upper bound, heat shock) in TASK §6 C and the notes. |
| 9 | Escher README and figures consistent with flux.json | supported, minor label issue | 0.00686 × 0.3116608 = 0.0021380 = R143 = R767 = drain; KO before growth 0, R387 = R326 = 0.3276; KO after = WT; +0.20 %; 12 LPs. |

## Defects and wording problems (reviewer's list) and what was done

1. "TEF1-driven YlTPS1 alone raised T6P synthase activity" overstates: CJM 613 keeps TPS2/TPS3 and is a different
   background from CJM 645. **Fixed** in TASK §4 and §6 D.
2. "below the detection limit": the paper gives no detection limit. **Fixed** in the curation notes (v2), TASK, REPORT.
3. "a mutant with no detectable trehalose synthesis grows at a WT rate": the single tps1 content at 30 °C is not
   reported, and tps1 grew slightly faster. **Fixed** (notes v2, TASK §6 C, REPORT §1).
4. The medium was tied to Table 4; it comes from Methods ("generally" 2 % glucose). **Fixed** (notes v2).
5. "below 1 in stationary phase on glucose" narrows the paper's phrase. **Fixed** (notes v2 quote the phrase).
6. Other note values not pinned. **Documented** (docstring, REPORT §2): keys only, as in the dolichol step.
7. Objective direction unchecked. **Fixed** (guard + test case).
8. Exact coefficient vs `isclose`. **Fixed** (same `isclose`).
9. "TERM REMOVED" label on the before panels. **Fixed** ("present before, removed after"; re-rendered, no new LP).
10. README UDP-glucose numbers. **Fixed** (R638/R793 0.2962 → 0.2946; R4 0.2940 → 0.2946).
11. UniProt 469 aa vs Flores' 486-aa ORF. **Noted** in REPORT §6, not investigated.
12. `059316ef` "re-hashed" and "Gurobi 13" not verified by the reviewer. The main session hashed a copy of
    `059316ef` before drafting TASK.md (`059316ef34417c88…`) and read Gurobi 13.0.2 from `gurobipy` in
    `/opt/miniconda3/bin/python`.
