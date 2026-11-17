**Status: done; see [REPORT.md](REPORT.md) and [AUDIT.md](AUDIT.md). Candidate `ad901ddb`. Every acceptance criterion
in section 2 is met. Not committed.**

# Methylcitrate cycle on 7ccbb568: matrix 2-oxobutanoate dehydrogenase added, R2075 closed, 2026-11-16

One `curation-change` round on this line (branch `claude/bold-bartik-794e7e`).

## 0. Authorization

The user, 2026-11-16, in this session:

1. **The option.** After [../r1172_side_effect_20261116/REPORT.md](../r1172_side_effect_20261116/REPORT.md), asked
   "Which option should be built for the nine reactions?", the user answered: "O1 + R2075 closed".
2. **The plan.** Asked "Approve this build plan, including the full certified sweep?", the user answered:
   "Approve plan". That covers the steps, budget and stop conditions below, including one full certified sweep.

Not in scope:
- side route B (R834 in, R2025 out, R2024 backwards, propanoate export) stays as it is; reported only;
- reopening R1172 (recorded user decision);
- the gene-less acyl-CoA transports across the inner membrane (reviewed in another session);
- any other edit to R834, R2025, R2113, R2112, R2074, R2076, R2077.

No commit or push without a separate instruction.

## 1. Starting state (verified 2026-11-16)

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/bold-bartik-794e7e` @ `af7c543214459f22d3947d71c8fc1a75ea5fb800` (PR #54 merge). Uncommitted: `model/STATE.md` and `model/reports/r1172_side_effect_20261116/` (this session's investigation) |
| Baseline model | `model/candidates/E5_353cefdc_xref_20261115.xml`, `7ccbb56834a37ad4c75e434868e9e1b9e1f6a657a11ca17ee0a93dbe56232146`. Not a designated reference |
| Baseline build arguments | that candidate's `.build.json` argv (head `79109038`), in [scripts/args_7ccbb568.txt](scripts/args_7ccbb568.txt) |
| Medium / profile | Task A medium `6eca70fd…`, PO1f profile `35307853…`, via `essentiality_simulation_context.py` |
| Essential list | `consensus_essential_genes.csv` `1e887f5a…` |
| Recorded sweep of the baseline | `$IYALI26_RESEARCH_ROOT/artifacts/tasks/folinic_acid_module_20261112/sweep684_candidate.tsv` (`49b4f27b…`), run on 353cefdc. 353cefdc and 7ccbb568 are LP-identical (static check, [../r1172_side_effect_20261116/REPORT.md](../r1172_side_effect_20261116/REPORT.md) §1), so it applies to 7ccbb568 |

### Source of the step

The step and its data come from the methylcitrate-cycle round of another line (worktree
`gene-essentiality-reactions-d09601`, uncommitted there), which built the same two edits on cbf2ef95. They are copied
here byte-identical:

| File | SHA-256 |
| --- | --- |
| `platform/scripts/gem_annotate/methylcitrate_cycle.py` | `79dabdd05dad152279b975c088d4d02c6bfdb959fabb50fbf94c43a1a73f809e` |
| `model/curation/methylcitrate_cycle/methylcitrate_cycle_20261116.json` | `17367e6f33bcfd2945096aa11494748408a8771a17994f671a6bfa1104b11485` |

Because the copies are byte-identical, three fields name the other line, not this one:
- the data's `baseline_model_sha256` is cbf2ef95;
- its `report` points to `model/reports/methylcitrate_cycle_20261116/REPORT.md`, which is not on this branch;
- the module docstring says "On cbf2ef95".

The step does not read `baseline_model_sha256`. It pins the state of its targets, R693, the nine species, four
species' reaction lists and the five PDH genes' reaction links. A dry run on 7ccbb568 in memory (this session) returned
"applied". So that pinned state is the same in cbf2ef95 and 7ccbb568. The flag help text here says so. When the two
lines merge, the flag will be defined on both; the merge must keep one definition.

The other line's test (`test_methylcitrate_cycle.py`) builds its baseline through steps this branch does not have
(`rnr_thioredoxin`). It is not copied. A new test runs the same checks on the tracked 7ccbb568 file.

### Targets in 7ccbb568 (from [before.json](before.json))

| Item | Before | After |
| --- | --- | --- |
| R_OBDHm [C_mi] | absent | 2-oxobutanoate m42 + CoA m73 + NAD m27 → propionyl-CoA m78 + CO2 m44 + NADH m30 + 4 H+ m28, (0, 1000); rule copied from R693 (PDH, five genes, AND), marked INFERRED; group and subsystem "Propanoate metabolism"; annotation SBO:0000176, Rhea 33223; four notes |
| R2075 [C_cy], (R)-methylmalonyl-CoA mutase, gene-less | (0, 1000), no notes | (0, 0), three notes; kept in the model |
| The five PDH genes `YALI1E31947g`, `YALI1F27556g`, `YALI1D26360g`, `YALI1B13179g`, `YALI1D30618g` | reactions as in before.json | plus R_OBDHm |

Nothing else may change.

## 2. Expected effect and acceptance criteria

1. A control build (same argv, no new flag) is byte-identical to `7ccbb568`.
2. Candidate builds with `PYTHONHASHSEED` 0 and 1 are byte-identical to each other.
3. The whole-model diff against `7ccbb568` shows only the targets above. Gene-reaction links follow from the rule.
4. R_OBDHm is mass- and charge-balanced.
5. Solves on Task A + PO1f, Gurobi 1e-9, basis reset per solve, on control and candidate:
   - **WT growth:** unchanged within 1e-9 (in memory: 0.311100535925948 as is, 0.31110053592594944 with the option).
   - **FVA** of the nine (R25, R552, R95, R2074, R2075, R2076, R2077, R2112, R834) and R_OBDHm, with no growth floor.
     Expected from the in-memory run: R25, R95, R552 up to 7.5987; R2074, R2075, R2076, R2077, R2112 at 0/0;
     R834 up to 3.3108.
   - **Energy check:** every boundary lower bound and the `xMAINTENANCE` lower bound set to 0; maximum `R_NTP1` and
     `biomass_C` both 0.
6. **Full certified single-gene sweep** on the candidate (`ko_screen_certified.py` `313d5d01…` with `checked_solve.py`
   `69c42690…`, unchanged), compared gene by gene with the recorded baseline sweep. No gene may change class at 1, 5, 10
   or 15 %, under rule A or B. Any change is reported, not accepted.
7. **pytest** before and after, on short-path `git archive` copies of the tree, root `.venv`,
   `--continue-on-collection-errors`, compared test by test. No existing test may change outcome; the new test passes.
8. One independent read-only audit.

## 3. Compute budget and stop conditions

- Budget: at most 1,300 LPs in all (targeted solves about 50; the sweep about 1,100 plus re-solves and certificates).
- Builds: three (control, candidate seed 0, candidate seed 1). Each about a minute, no LP (`--no-solve`).
- Stop and report, without accepting the candidate, if:
  - the control build is not byte-identical, or the two seeds differ;
  - the diff shows anything outside the targets;
  - R_OBDHm is unbalanced;
  - WT growth or the FVA values differ from section 2;
  - the energy check is above 0;
  - any gene changes class, or any existing test changes outcome;
  - the budget would be exceeded.

## 4. Outputs

- This folder: TASK.md, before.json, scripts/, logs/, model_diff.json, targeted solve and sweep summaries, REPORT.md,
  AUDIT.md, figure/.
- `model/candidates/E5_7ccbb568_methylcitrate_20261116.xml` (plus `.build.json`, `.coq9_genes.tsv`).
- Large outputs (build copies, sweep TSV, pytest junit files) in
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/methylcitrate_cycle_7ccbb568_20261116/`.
- A `model/STATE.md` entry.
