# Independent audit (2026-10-29)

A fresh read-only general-purpose subagent, given only TASK.md, REPORT.md, before.json, the diff, the sweep outputs, the
logs, the curated JSON, the step, the test, the `cli.py` / `main.py` diff and the two models. It edited nothing, wrote
its own comparison and LP scripts in the session scratchpad, ran 10 LPs (Gurobi 1e-9, reset per solve) and the focused
test. Its report is summarised here; the main session did not re-run its scripts.

| Claim | Verdict | Auditor's evidence |
| --- | --- | --- |
| C1 model SHAs | supported | `shasum` of both files and of SHA-checked scratch copies; candidate `.build.json` agrees |
| C2 diff confined to targets | supported | raw XML diff: 33 lines in 3 hunks (6 notes + one `fbc:or` collapsed, per reaction); own cobra comparison of every field incl. groups and id order: 12 differing fields = rules, notes and gene sets of R723/R759/R760 + reaction links of 3 genes |
| C3 after-state | supported | R759 `YALI1F14882g`, R760 `YALI1A09121g`, R723 `YALI1F03908g`, R673 unchanged, `YALI1B07957g` with no reaction; 6 notes per reaction after an SBML read, no `<` or `&` |
| C4 curated JSON vs models | supported | 0 mismatches before-state vs baseline, 0 after-state vs candidate |
| C5 guard and test | supported (limit) | all checks precede the first edit; partial state rejected; idempotent; context restore; focused test 6 passed, 21 subtests. Limit: stale-input cases mutate only R723 / `YALI1B07957g` |
| C6 wiring | supported | additions only: flag, parameter, call after the R2171 block, pass-through; no dependency check |
| C7 control build | supported (caveats) | argv differ only by the flag; control XML hashes to `965ab947…`; identical R304 ERROR lines. Caveats: abbreviated candidate line in build_commands.txt; control in scratch |
| C8 WT | supported | 0.3116608216 on both |
| C9 target knockouts | supported | candidate: `YALI1F14882g` 0, `YALI1F03908g` 0, `YALI1A09121g` and `YALI1B07957g` WT; baseline: all WT |
| C10 sweep | supported | recomputed from `ko_compare.tsv`: matrices, MCCs and classes match; exactly 2 changes per cutoff; 0 failing; 20 infeasible on each, same set |
| C11 pytest | supported | own parse: only the 6 new tests differ, all passed; FAILED/ERROR/SKIPPED lines identical; before log older than the code files |
| C12 no overclaim | supported | compartments inferred; development reference and user-defined negatives stated; recall stated as consequence |

**Other findings and how they were handled.**
1. The AUDIT.md link was dangling: this file.
2. REPORT §5 lacked recall, specificity, precision, coverage, out-of-scope positives and input SHAs (benchmark.md):
   added.
3. The plan's artifact location for the sweep TSV was not followed: recorded in REPORT §8 as a departure.
4. The copied `ko_compare_two.py` docstring names the R2171 task's models: noted in REPORT §8; the copy stays unchanged.
5. The R304 error in the R2171 logs was outside the auditor's material: the main session's grep this session found it
   in five R2171 build logs; REPORT §3 names two of them.
6. The "complex III / IV" wording for the unresolved set was ambiguous: reworded in REPORT §5 to "4 mtDNA and 16
   nuclear infeasible knockouts, the same set on both models".
7. The C7 caveats: the control `.build.json` is copied into `logs/`; REPORT §3 says where the full argv is.
