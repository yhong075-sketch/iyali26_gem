# Independent audit of the R1172 Stage 2 (2026-11-12)

A fresh read-only agent was given only TASK.md, REPORT.md section 7 (with section 6 items 4–5), escher/README.md, the
uncommitted diff, the logs and outputs in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r1172_hmgcoa_transport_20261112/`
and the two short-path model copies. It wrote nothing and ran 13 LPs of its own. The pytest row was a placeholder
at the time and was not audited.

## Verdicts

| Claim | Verdict | Auditor's evidence |
| --- | --- | --- |
| Builds and SHAs (control = 2f352037; candidate = seed 1 = 684f3752; argv = 2f352037 argv + the flag) | SUPPORTED | own SHA-256 of all four files; `.build.json` argv; all 147 data and source SHAs in the candidate `.build.json` match the files on disk |
| Whole-model diff = R1172 bounds + notes only | SUPPORTED | raw text diff: 9 lines, all in R1172 (`reversible`, bound parameters, 7 added notes); seed-1 diff JSON byte-identical |
| TASK §8 targets respected (R412 untouched; R82 / R411 / R410 pinned, unchanged) | SUPPORTED | none of them in the diff; the guard never references R412 or R2121 |
| `later_decision_20261112` does not reach the XML | SUPPORTED | `annotate_retained_reactions` fields; the control build used the edited `retained_reactions.json` (`bbfa497e`) and still reproduced 2f352037 |
| Balance (R1172 balanced; R410 −2 in both) | SUPPORTED | own atom and charge counts |
| WT 0.3111005359259494 and no-uptake ATP maximum 0 | SUPPORTED, reproduced | own checked solves: candidate 0.3111005359259494, base 0.3120526953810365, `R_NTP1` maximum 0 (optimal) |
| Sweep byte-identical to `ko_screen_S1.tsv` | SUPPORTED | `cmp`; both `f1178cb7…`; same 20 certified knockouts; reference script `313d5d01` |
| Score tables and class changes | SUPPORTED, recomputed | every cell recounted (rules A and B, both scopes); only YALI1B11532g and YALI1F38116g change against A2, plus the A2 pair against v4; no TP → FN, none above WT |
| LP counts 1,123 + 2 + 60 | SUPPORTED as script tallies | sweep summary, checks, `flux.py` (30 per model) |
| Escher counts and FVA numbers | SUPPORTED, recounted / partly reproduced | `flux.json`, viewer logs; R1172, R834, R2119, R82 ranges reproduced |
| Code-review table contents | SUPPORTED | build record, test file, focused log, `_check_spec`, XML diff |
| Code-review provenance | UNVERIFIED | no log or transcript of the reviewer existed; its one factual claim was reproduced |
| §7.1 step description; curated JSON evidence tiers; §6 items 4–5; AGENTS §5 citation | SUPPORTED | code reading; TASK §3.1; flux.json; AGENTS.md line 150 (`488a0d56`) |

## Problems raised and what was done

| # | Problem | Action (2026-11-12) |
| --- | --- | --- |
| 1 | Report not final: AUDIT.md missing, AUDIT / PYTEST placeholders, no Stage 2 STATE.md entry | this file; placeholders filled in REPORT.md 7.3; STATE.md entry added |
| 2 | Code-review provenance unverified | the reviewer's report is recorded in `stage2/code_review.md` (summary; the agent transcript is not kept), so it stays a reported review, not a reproduced one |
| 3 | Commands for `compare_built_vs_*.json` not logged | `stage2/logs/compare_commands.txt`; the A2 base TSV is `r2121_malonyl_coa_20261110/stage2/ko_screen_built.tsv` (`6d1963f5…`, identical to A2 `ko_screen_v1.tsv`) |
| 4 | `built_checks.json` keys said `v1` though the reference is S1 | script labels changed to `s1` and re-run (2 LPs; then 2 more after a last label fix): same values, `criterion_5_sweep_reproduces_s1` true |
| 5 | §7.4 wording: R834 returns CoA, not mevalonate; R553 carries 0.0545, not all of R410's flux | §7.4 reworded |

**Carried forward (auditor):** under root AGENTS.md §5, R1172 at (0, 0) is itself an inactive reaction; the change
moves a gap rather than closing it. REPORT.md 7.6 states this.
