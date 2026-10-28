# STOP: Option A (R1726 opened) activates an acetyl-CoA → acyl-CoA leak; FAS and OLE1 stop being essential (2026-10-27)

Stage 2 stopped at TASK.md section 10, under the stop condition "a knockout class change not traceable to the
R1726 route". The changes *are* traceable to R1726, but the route they trace to is a pre-existing
model artefact, not biology. Nothing was placed in `model/candidates/`. No STATE.md entry and no commit were
made. The decision on how to continue is the user's.

## What was done (all verified this session)

| Step | Result |
| --- | --- |
| Curation data | `model/curation/direction/r1726_cpt_open.json` (R1726 [0, 0] → [0, 1000]; R1725 kept closed, guarded on name, GPR, bounds, stoichiometry and species) |
| Step and flag | `platform/scripts/gem_annotate/r1726_cpt_open.py`; `--r1726-cpt-open` in `cli.py` / `main.py`, requires `--cofactor-biomass-terms` and runs after it |
| Focused test | `platform/tests/test_r1726_cpt_open.py`: 5 tests, every guard hit with its own message |
| Code review | Fresh general-purpose reviewer following `.claude/agents/code-reviewer.md`: no bugs, 4 low findings, all fixed (comment, R1725 stoichiometry guard, guard tests, the flag requirement). No change reaches the model bytes. |
| pytest | 269 → 274 passed; the difference is the 5 new tests only; the same 11 failures and 18 errors ([pytest_compare.txt](logs/pytest_compare.txt)) |
| Builds (round 1, artifacts folder) | Control without the flag = `00e8e397…` byte for byte. `PYTHONHASHSEED` 0 and 1 with the flag both give `95316a852c5cc59b2cdefed176e134770b04a277e8a2a936dd2ab2b14908019b` ([build_argv.txt](logs/build_argv.txt)). |
| Whole-model diff | 2 entries: R1726 bounds, R1726 notes ([model_diff_vs_00e8e397.json](model_diff_vs_00e8e397.json)) |
| Built-file checks | ([built_checks.json](built_checks.json))<br>- WT 0.3120049 → 0.3120558, checked.<br>- No-uptake max ATP and max biomass are 0, optimal, on both models (xMAINTENANCE relaxed).<br>- R1725 stays [0, 0] in FVA.<br>- PSD double knockout stays lethal.<br>- With choline supplied, PE methylation becomes dispensable. |
| Re-census (GLPK) | ([census_compare_open.json](census_compare_open.json), [census_compare_taskA_po1f.json](census_compare_taskA_po1f.json))<br>- Open condition: 722 → 710 blocked. Freed: UM12's 11 reactions plus R1344, which the cofactor task classed as GLPK fringe.<br>- Task A: 852 → 843.<br>- No reaction newly blocked; gap metabolites 675 → 668. |
| Knockout sweep | 2 × 1,082 LPs, Gurobi 1e-9, reset per knockout, knocked-flux check, no failing point ([ko_compare_summary.json](ko_compare_summary.json)). **4 class changes at every cutoff**, all essential → non-essential. |

## The four class changes

| Gene | Model reactions | Consensus list | KO / WT, baseline → candidate |
| --- | --- | --- | --- |
| `YALI1B19844g` | fatty-acid synthase (R1392, R1393, the acyl-ACP hydrolases; AND with `YALI1B25427g`) | yes (TP → FN) | 0 → 1.000 |
| `YALI1B25427g` | same AND pair | yes (TP → FN) | 0 → 1.000 |
| `YALI1C07638g` | acyl-CoA Δ9 / Δ12 desaturase R1434, R1435, R1869 | yes (TP → FN) | 0 → 1.000 |
| `YALI1D22124g` | fatty-acid–CoA ligases (24 reactions) | no (FP → TN) | 0 → 1.000 |

At 15 %: TP / FN / FP / TN go from 111 / 204 / 42 / 692 to 108 / 207 / 41 / 693 (MCC 0.3833 → 0.3769 from the
summary file). The negatives are the user-defined class; mtDNA genes are out of scope; these are
development-reference scores. The score is reported, not used as a reason.

A strain lacking fatty-acid synthase that grows at the wild-type rate without a fatty-acid supply is not
credible (established / uncited for yeast). Three of the four genes are on the consensus essential list.

## Mechanism (verified: [acyl_leak_diagnosis.json](acyl_leak_diagnosis.json), [r2135_control.json](r2135_control.json))

The FAS knockout on the candidate keeps WT growth. It becomes lethal again when any one of these is closed:
R1726, R2145, R1771, R1776 or R1778. It stays viable when the PLD pair or the phospholipases B are closed
one at a time (alternates exist). The pFBA point of the FAS knockout runs this loop:

1. **R2135 "AC-CoA transport"** maps **acetyl-CoA** `m68[C_cy]` (C23H38N7O17P3S, BiGG `accoa`) to the
   **generic acyl-CoA** `m1984[C_er]` (no formula, MNXM1136895).
2. **R2168** (ER DGAT): `m1984` + DAG `m2011[C_er]` → TAG `m2014[C_er]`, so the TAG carries acetyl groups
   that the network counts as acyl chains.
3. TAG transports R2170 → R1707 → R1785 move it to the lipid particle `m1640[C_lp]`.
4. **R1771** (DGAT, reversible) runs backward: TAG + CoA → DAG `m1642[C_lp]` + **pooled long-chain
   acyl-CoA `m570[C_lp]`**, which feeds every acyl-CoA consumer.
5. The released DAG returns as `m1670[C_cy]` → `m1649[C_em]` (R1778, R1776). **R1726** turns it into PC,
   PLD R1842 turns PC into PA `m378`, and R2144 + **R2145** (written with 0.01 coefficients) turn PA into
   the other DAG species `m1994[C_cy]` → `m2011[C_er]`. That closes the loop.

Net effect: acetyl-CoA becomes long-chain acyl-CoA without fatty-acid synthase or desaturase. Before
R1726 was opened, the loop could not close: the only other DAG → PA step, DAG kinase R1720, is (0, 0).

**The leak already exists in the baseline.** R2135 carries obligatory flux in the baseline WT (FVA at 99.9 %:
[1.9e-5, 2.4e-3]). Closing R2135 lowers WT growth on both models (0.3120049 / 0.3120558 → 0.3116824). With
R2135 closed, all four knockouts are lethal again on the candidate. R2135 belongs to the R2134–R2194 block,
whose species carry the `_TB` BiGG suffixes and no formulas; its provenance was not traced in this task.

## Options for the user

1. **(Recommended) Pause R1726 and investigate R2135 and the R2134–R2194 block first** as a separate,
   scoped task. That task would cover:
   - the identity of `m1984` / `m570`;
   - whether R2135 should carry acetyl-CoA at all;
   - the 0.01-scaled R2145;
   - the WT growth that depends on R2135.

   After a reviewed fix, re-run this task's Stage 2. The data, step, test and flag are ready and
   unchanged; the guard pins only R1726 and R1725.
2. **Accept Option A as built** and record the four changes as a known artefact. Not recommended: the
   candidate would predict FAS-null growth.
3. **Drop the change.** R1726 stays closed with the reasons recorded. The untracked files of this task
   would then be removed, but only on your instruction.

None of these options closes or tunes anything to restore a match. Option 1 is a separate science change
that needs its own authorization.

## State of the worktree

Untracked, not committed:
- the curation file, step, test and report folder;
- modified `cli.py` / `main.py`.

The round-1 builds, census outputs and the pFBA difference are in
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/kennedy_closures_20261027/`.

Not done because of the stop:
- the candidate in `model/candidates/`;
- the independent audit;
- REPORT.md;
- the STATE.md entry.
