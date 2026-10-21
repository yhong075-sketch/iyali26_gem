# TASK: screens of candidate b3b9fe22 (no overlay and diagnostic overlay), PO1f + CSM-Leu, glucose 2.414 (2026-10-21)

Covers this folder and [../screen_b3b9fe22_overlay_w29_20261021/](../screen_b3b9fe22_overlay_w29_20261021/).

## Authorization (user, 2026-10-21, this conversation)

"make a new one by our newest model" (after the iYali21 vs a56f3649 overlay-lethal figure). After the
plan the user chose: model **b3b9fe22**; figure condition **overlay + lethal view**, which includes
both screens of the new model (about 2,310 LPs) so the comparison is complete. Computation and a
figure only; no model change; no commit.

## Starting state

- Branch `claude/confusion-matrices-iyali21-7e63ca`, HEAD `4569682`, clean.
- Research root `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`.
- No reference model is designated.

## Model

`E5_coq9_resp_pkg_r132_combined_lipoylation_acylacp_20261020.xml`, SHA-256
`b3b9fe22cba4c8ae4408afe61740bfe75e661c383dd47a224965f06469618543`: the newest build, committed in
`9f2356b` on the unmerged, unpushed branch `claude/combined-candidate-20261019` (combined candidate
`164f235f` on main + `--lipoylation-module` + `--acyl-acp-chemistry`). Read from a byte-identical
`git show` export in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/screen_b3b9fe22_w29_20261021/input/`.
1,077 native genes.

## Inputs and method

Same medium (`6eca70fd…`), PO1f profile (`35307853…`), consensus list (`1e887f5a…`), mtDNA rule and
overlay (`holistic_overlay.py`, `bf6c7ce5…`; pre-checked: assertions pass on b3b9fe22, only forced
bound NGAM) as the a56f3649 screens. Runners are copies of
`confusion_matrices_csmleu_w29_20261014/run_screen.py` and
`screen_a56f3649_overlay_w29_20261014/run_screen.py` with only the model changed.

## Budget and stop conditions

No overlay: WT + 1,077 KOs = 1,078 LPs. Overlay: WT + 4 controls + 1,077 KOs = 1,082 LPs plus at
most 150 diagnostic LPs. Total at most 2,310. One thread, 60 s per LP, 900 s wall per run, no
retries. Stop on input/code hash change, assertion failure, invalid WT or budget.

## Scoring and outputs

As before (benchmark rule primary; lethal and maintenance views labelled; mtDNA out of scope;
user-defined negatives; development reference). Gene-level comparison with the a56f3649 screens
(no solves). Figure: iYali21 vs b3b9fe22, both overlay + lethal view, in
`figure_iyali21_vs_a56f3649_20261014/`. REPORT.md here and a STATE.md entry.
