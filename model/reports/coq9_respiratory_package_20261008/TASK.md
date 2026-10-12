# TASK: CoQ9 + respiratory-chain curation package (2026-10-08)

## Authorization (user, 2026-10-08, this conversation)

The user listed five items and asked: "once you finished the change, do the screen test
again. commit and push our new model and pipeline." After a plan review the user decided:

- Item 4 (proton leak): "this is a assumption. we can put it but need to mark it as 'assumption'".
- A (complex I): "Option 1, minus three genes, with each remaining gene checked against the
  complex I structure" — drop `YALI1M00338r`, drop `YALI1E06573g` for now (restore if
  identified), drop NDH2 `YALI1F32476g`; keep the reviewed 4-subunit rule as interim if the
  structure check cannot be done.
- B (fatty alcohol oxidase): skip. C (`R132`): skip; record as open hypothesis.
- D (complex IV): a separate cited pass — keep subunits whose yeast deletion abolishes
  complex IV assembly or activity, remove proven non-members, keep mtDNA-encoded subunits in
  the rule but out of benchmark scope, settle the three `R2206` genes before merging.
- E (push target): "model-platform" → `restructure/model-platform`.

The package is authorized as one unit; each item still has its own curated data, guarded
`apply_*` step and test. Items applied: A1, A2, B1 (with B2), B3 + D, C1, C3, C4.
Not applied: fatty alcohol oxidase relocation, `R132` direction.

## Starting state

- Branch `claude/coq9-respiratory-package` from `restructure/model-platform` at
  `ea5b639dd55cf833d5a880cf2d55551b12bc4677` (worktree `coq-system-problem-75cb14`).
- Latest delivered model: `reports/coq_r305_candidate_20261005/E5_coq9_alpha_1e-4_R305_qcycle.xml`,
  SHA-256 `33468b94f378f0ab7f29f05e4356fd93b6d3eff23bc112d566cf16d6e4d6041c`. It was produced
  by hand-patching XML (R305 regex edit) and a standalone CoQ biomass tool, so it cannot be the
  input of a pipeline build.
- Step 0 (this task): the builder learns to produce that model — the curated CoQ9 biomass term
  moves into the builder, E5 may combine with `--coq9-curation qcycle`, and candidate stages
  start from the serialized reference (fixes the known export/reload failure: float noise in
  3 acyl-pool reactions and 20 tRNA residues with unknown charge). The Step 0 build differs
  from `33468b94` only in R1889 GPR + notes (default 4-subunit rule since 2026-10-05), the
  R305 proton-state note text and one R570 note. This Step 0 model is the baseline for all
  diffs and comparisons below.
- `before.json` records every target entity in the Step 0 model.

## Targets and intended changes

| Item | Target | Change |
| --- | --- | --- |
| A1 | R39, R969, R808, m108[C_cy], m110[C_cy] | R39 into C_mi (m641 + ½ O2 → H+ + m939, all C_mi), GPR `YALI1A08781g and YALI1B20527g`; delete R969, R808 and the two orphaned cytosolic intermediates |
| A2 | R695 | GPR `YALI1E18269g and YALI1F34675g` (existing curated `coq9_functional_gpr.json`) |
| B1/B2 | R1889, R2062 | delete R1889; R2062 pumps 4 H+ (NADH + Q9 + 4 H+ mi → NAD + Q9H2 + 4 H+ cy); GPR = structure-checked complex I subunits, without NDH2, YALI1E06573g, YALI1M00338r |
| B3/D | R2206, R304 | per-subunit cited GPR; settle R2206 genes; delete R2206; R304 pumps 4 H+ (8 H+ mi + 4 cyt c red + O2 → 4 H+ cy + 2 H2O + 4 cyt c ox) |
| C1 | R_CYSS_m | made cytosolic, irreversible and balanced it equals R247 (same gene); delete as duplicate |
| C3 | R2004 | delete (no gene; bacterial transcarboxylase EC 2.1.3.1) |
| C4 | R349 | delete the matrix copy of NAD-G3PDH; the enzyme remains as cytosolic R348 |
| Item 4 | R1384 | screen-condition assumption: lower bound 8 mmol/gDW/h; not written to the XML |
| Item 5 | benchmark | the 16 mtDNA-encoded model genes are out of scope for the nuclear CRISPR benchmark |

Genes are never deleted; genes left without reactions stay in the model.

## Acceptance criteria

1. Whole-model diff (package vs Step 0) shows only the targets above; every touched reaction
   is element- and charge-balanced.
2. Focused tests pass; the full test suite shows no new failure versus the clean ea5b639 tree.
3. Controls on the package model: WT optimal and > 0; R385 closed → zero growth; closed-input
   ATP dissipation (cytosol, mitochondrion) = 0 under the existing protocol.
4. Screens report 1/5/10/15 % (15 % primary) under the current benchmark, with and without
   the mtDNA scope rule, and with and without the leak assumption, labelled as such.

## Compute budget and stop condition

- Two offline/no-solve builds (Step 0, package).
- Controls: 4 LPs on the package model.
- Screens: Step 0, package, package + leak assumption — at most 3 × 1,074 LPs, 15 min wall
  each, Gurobi 1 thread, settings identical to `screen_coq_r305_20261008`.
- Stop on an input SHA change, a non-optimal WT, or an unexpected diff entry.
