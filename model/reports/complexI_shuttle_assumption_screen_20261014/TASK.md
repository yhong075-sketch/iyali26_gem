# Benchmark side effects of Task D options (a) and (b), plus homoisocitrate DH location

Written 2026-10-14, before any run. Diagnostic only: no model XML, condition, assumption, curation or evidence
file is written.

## Authorization (user, 2026-10-14)

"commit. ok. run (a) and (b), check. review." Read back and approved (AskUserQuestion): commit Task D after
removing the repository copy of `combinatorial_tier1.csv` (done, `8718b75`); run (a) and (b) in a benchmark
context and in Task D's diagnostic context; check the location of homoisocitrate DH; independent read-only audit
of the new runs plus Task D's `capacity_costs.py` and fine shuttle-cap grid. Results stay uncommitted.

## Question

Task D §7 offered two in-memory ASSUMPTION closures that make complex I essential in config A:
- **(a)** matrix MDH R533 reductive direction closed (lb 0);
- **(b)** Ctp1-type carrier R831 forward closed (ub 0; citrate in, malate out).

Their effect on the other gene calls was untested. Does either change gene calls besides complex I, and does
either still make complex I essential without the holistic overlay and 13C caps?

## Design (fixed before any run)

- **Model:** Task B candidate `E5_coq9_resp_pkg_r132_20261012.xml` (`a56f3649…`).
- **Medium / profile:** Task A medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`), PO1f profile
  `po1f_sd_leu.json` (`35307853…`).
- **Contexts:**
  - `bench` (primary): as loaded, no overlay, no caps. The same context as Task B's screen.
  - `diag` (diagnostic): `bench` + Task D's holistic overlay + 13C caps (Task D config A).
- **Variants:** `base`; `a` = R533 `rev` closed; `b` = R831 `fwd` closed, using Task D's `close()` (intersects
  bounds, never opens).
- **Runner:** Task B's `screen()` and `metrics()` from `run_task_b.py`, unchanged, imported by path. The wrapper
  only replaces `context()` to apply the context and variant after loading. Solver parameters are Task B's.
- **Per run:** WT, KO/WT for R2062/R305/R304/R171, all 1,073 gene KOs, route FVA and pFBA for WT and the
  complex I KO. About 6 × 1,220 LPs.
- **Scoring:** `benchmark.md`: cutoffs 1/5/10/15 % (15 % primary), mtDNA genes out of scope, negatives =
  model genes not on the consensus list (user-defined class), no solution = unresolved; second labelled view
  counts no solution as lethal. Development reference, not independent validation.
- **Reported:** every gene whose ratio or call changes between `base` and `a`/`b` within each context.

## Reproduction checks (stop and report before interpreting if any fails)

- `bench/base` at 15 %, benchmark rule: TP 82 / FP 33 / FN 229 / TN 697 / unresolved 16 (Task B run.json);
  WT 0.3124267 ± 2 %.
- `diag/base` complex I KO/WT 0.826 ± 2 %; `diag/a` and `diag/b` complex I KO/WT 0.077 ± 2 % (Task D §7).

## Homoisocitrate DH location check

R398 (homoisocitrate DH, YALI1D13238g) is the forced matrix-NADH source in Task D's mechanism.
- Run TargetP-2.0, MitoFates and DeepLoc 2.0 on its protein, same tools as Task B. If a web tool cannot be
  reached, stop and ask; no heuristic substitute.
- Rule fixed now (Task B's): mitochondrial if TargetP mTP ≥ 0.5 or MitoFates ≥ 0.385.
- Search for *Y. lipolytica* localization data and list every model reaction for homoisocitrate DH by compartment.

## Outputs

REPORT.md, `run.json`, small CSVs here; per-gene tables and pFBA vectors under
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/complexI_shuttle_assumption_screen_20261014/`. A `model/STATE.md` entry.
