# One static essentiality screen — 2026-10-08

Authorization: the user's current request, “对我们的最新模型进行一次screen test”.
Scope: one fresh WT and one single-gene deletion for each of the 1,073 native
model genes; exclude the runtime plasmid pseudo-gene. No model edits, rebuild,
dFBA, FVA, energy tests, evidence research, repair loop, commit or push.

Starting branch: `restructure/model-platform`; HEAD:
`ea5b639dd55cf833d5a880cf2d55551b12bc4677`. Existing tracked dirty file:
`model/STATE.md`; preserve its contents. No scientific reference model is designated.

Executable input: `model/reports/coq_r305_candidate_20261005/E5_coq9_alpha_1e-4_R305_qcycle.xml`,
SHA-256 `33468b94f378f0ab7f29f05e4356fd93b6d3eff23bc112d566cf16d6e4d6041c`.
This is the previously delivered CoQ9 alpha=1e-4/R305 candidate, not the different
default builder output. Static checks confirm the CoQ9 biomass coefficient is
−0.0001 and R305 proton coefficients are −2 and +4.

Reuse current `model/expected/benchmark.md`: PO1f profile, SD-Leu, biomass_C
maximization; keep maintenance at 7.8625 and static nonlimiting uracil at 1000.
Input identities are pinned in run_screen.py and recorded in manifest.json.
The experimental CSV contains 1,612 positives, with 322 native direct-ID matches.
Unlisted model genes form the user-defined negative class, not experimentally
confirmed negatives. Missing model IDs are out of scope, never false negatives.
These data were used in development and do not provide independent validation.

Budget: exactly at most 1,074 optimization calls, one serial worker, Gurobi
Threads=1, Presolve=0, TimeLimit=60 seconds, feasibility/optimality/integrality
tolerances 1e-7. Aggregate screening budget 900 seconds; do not start another
call when fewer than 60 seconds remain. No retries or constraint relaxation.
Stop on an invalid WT, solver exception, input change, context-restoration
failure, budget exhaustion, or after all native KOs. Preserve any partial run.

Use COBRA gene.knock_out() to propagate the complete Boolean GPR. Preserve
raw status and objective. Score only optimal, finite, nonnegative growth with
mass/bound residuals at most 1e-6. Non-optimal, nonfinite, negative or numerically
invalid results remain unresolved. Compare unrounded KO/WT strictly below 1%,
5%, 10% and primary 15%; equality is non-essential. No clipping or failure-to-zero.

Small deliverables: this scope, executable runner, manifest, REPORT.md and a
dated STATE.md entry. Raw KO table, journal, WT fluxes, effective-model snapshot
and coverage table go to
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/screen_coq_r305_20261008/`.
One bounded independent read-only result check; no additional optimization.

Execution note: the first process stopped after WT + 40 KOs because COBRA's
context manager did not restore the functional flag of an orphan gene without
a model link. That gene had no reactions and changed no constraints. The
initial script, manifest, journal and table are preserved. A solver-free check
confirmed the cause and explicit flag restoration. Continuation verified the
same inputs, runtime context and complete effective-model snapshot, reused
the 41 completed calls, and executed only the remaining 1,033 KOs. The aggregate
budget remains 1,074 calls with no repeated optimization. The final manifest
is `manifest_continued.json`; `manifest.json` remains the stopped-attempt record.
