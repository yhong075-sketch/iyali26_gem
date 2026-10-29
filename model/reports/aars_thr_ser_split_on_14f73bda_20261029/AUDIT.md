# Independent audit (2026-10-29)

A fresh read-only general-purpose subagent, given only this folder, the curated JSON and the three models
(`14f73bda`, `d8d1bac7`, `6abb6b38`). It edited nothing, wrote its own comparison and LP scripts in the session
scratchpad and ran 10 LPs (Gurobi 1e-9, basis reset per solve). Summary of its report; its scripts were not re-run here.

| Claim | Verdict | Auditor's evidence |
| --- | --- | --- |
| D1 SHA-256 of the three models | supported | `shasum` on the files and short-path copies |
| D2 candidate vs `14f73bda` = 9 intended differences | supported | raw XML diff (33 lines, 3 blocks: R723, R759, R760) and own cobra comparison of all fields; counts 2319 / 1882 / 1081 / 113 on both |
| D3 candidate vs `6abb6b38` = UM7 additions only | supported | raw diff and cobra: R_ADPRtn, R_PPItn, group101 "Transport, nucleus" members |
| D4 targets vs curated JSON | supported | before-state on `14f73bda`, after-state on the candidate; R673 unchanged; `YALI1B07957g` present with no reaction |
| D5 build records | partially supported | control `.build.json` output `14f73bda…`; seed-1 XML in scratch hashes to `d8d1bac7…`; "applied" in both candidate logs; ERROR lines identical. Gaps: `build_commands.txt` garbled (seed prefix before every argument); seed-1 provenance not in the folder; `.build.json` does not record the seed |
| D6 WT and target knockouts | supported | WT 0.3116608216 on both; candidate KOs 0 / 0 / WT / WT, baseline all WT; residual ≤ 3.2e-14 |
| D7 scores | supported | recomputed from `ko_compare.tsv`: matrices, recall, specificity, precision, MCC match; exactly 2 class changes per cutoff; 0 failing points; same 20 infeasible |
| D8 no overclaim | partially supported | stated correctly: inferred localisation, development reference, user-defined negatives. Gaps: dangling AUDIT.md link; two statements resting on other reports; "every recent build log" wider than checked |

**Findings and how they were handled.**
1. `build_commands.txt` garbled: regenerated from the three `.build.json` argv lists, with the seed per run.
2. Seed-1 provenance missing: `logs/candidate_seed1.build.json` added; REPORT §2 states that the seed is recorded only in
   the commands file.
3. AUDIT.md missing: this file.
4. Cross-report statements: REPORT §3 and §4 now cite the files read this session (the combined task's diff log, the
   first build's sweep summary).
5. "Every recent build log": reworded to the logs actually checked.
6. R304 notes are dropped from the written models (auditor; confirmed in the `14f73bda` XML this session): REPORT §2
   now states it as pre-existing data loss.
7. The copied `ko_compare_two.py` docstring names another task's models; the summary JSON does not record the list's
   SHA. Both recorded here; the copy stays unchanged and REPORT §4 gives the list SHA.
8. The two mitochondrial-type paralogs staying FN: stated in REPORT §4.
