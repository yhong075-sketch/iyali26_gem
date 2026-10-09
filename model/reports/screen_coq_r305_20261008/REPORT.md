# Latest CoQ9/R305 candidate: one essentiality screen

Completed 2026-10-08 under the user's request for one screen of the latest model.
**All 1,073 native-model single-gene deletions returned optimal, numerically
valid solutions; none was unresolved. At the primary 15% cutoff, 120 genes are
predicted essential.** WT growth is **1.4445584333 h⁻¹**.

## Model and test conditions

Executed [E5_coq9_alpha_1e-4_R305_qcycle.xml](../coq_r305_candidate_20261005/E5_coq9_alpha_1e-4_R305_qcycle.xml),
SHA-256 `33468b94f378f0ab7f29f05e4356fd93b6d3eff23bc112d566cf16d6e4d6041c`.
Static input checks confirmed CoQ9 biomass demand −0.0001 mmol/gDW and R305's
final proton coefficients −2 and +4. This is the delivered CoQ9/R305 candidate;
it is different from the separately documented default build. No scientific
reference model is designated, and no reference model was screened here.

Fresh static FBA used biomass_C maximization, the existing SD-Leu medium and
PO1f runtime profile. Glucose uptake remained 10, maintenance lower bound
7.8625, and the existing static nonlimiting uracil uptake 1000
mmol/gDW/h. These are the static conditions, not the finite-uracil dFBA setup.
The strain overlay was applied in memory. Its plasmid pseudo-gene was excluded;
native loci affected by the strain overlay remain in the table and are flagged
as runtime-confounded. No model file, GPR, medium or scientific parameter changed.

Branch `restructure/model-platform`, HEAD `ea5b639dd55cf833d5a880cf2d55551b12bc4677`.
Full input/source hashes, versions, dirty state, active medium, runtime overlay
and commands are recorded in the [final manifest](manifest_continued.json).

## Results against the development reference

Per the [current benchmark](../../expected/benchmark.md), essential means the
**unrounded KO/WT growth ratio is strictly below** the cutoff. All four rows
below summarize the same 1,073 KOs; no extra optimizations were used.

| Cutoff | Predicted essential | TP | FN | FP | TN | Recall | Precision | Specificity |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1% | 87 | 61 | 261 | 26 | 725 | 18.94% | 70.11% | 96.54% |
| 5% | 95 | 67 | 255 | 28 | 723 | 20.81% | 70.53% | 96.27% |
| 10% | 104 | 71 | 251 | 33 | 718 | 22.05% | 68.27% | 95.61% |
| **15% — primary** | **120** | **83** | **239** | **37** | **714** | **25.78%** | **69.17%** | **95.07%** |

The reference contains **1,612 positive IDs**, of which **322 (19.98%)** directly
match native model IDs. **1,290 positives are out of scope** because their exact
IDs are absent; they are not counted as FN. No cross-version mapping was inferred.
Recall uses the 322 covered positives as its denominator.

**The 751 unlisted model genes are the user-defined negative class.** FP and TN
use that definition; absence from a positive-only source is not an experimental
demonstration of nonessentiality. The input parser preserves its historical
positive-only metadata; the explicit scoring rule in this report and final
manifest follows the current benchmark. This reference was used in model
development, so these scores are **not independent validation**.

The result means that 239 covered experimental positives remain predicted
nonessential at 15% under these conditions. It does not identify the biological
cause of those mismatches, establish native gene essentiality, or show that a
particular repair improved performance: no comparison screen was run.

## Execution, checks and limits

Exactly **1,074 optimizations** were executed: one fresh WT and one KO per native
model gene. Gurobi used one thread, Presolve=0, 60 seconds per call and 1e-7
feasibility/optimality/integrality tolerances. Raw solver status and growth were
preserved; no failed solve was converted to zero. Scoring also required finite,
nonnegative growth and mass/bound residuals at most 1e-6.

The first process stopped after 41 completed calls because an orphan gene's
temporary functional flag was not restored by COBRA's context manager. Its KO
changed no reaction constraints. The [initial runner](run_screen_initial.py),
[stopped manifest](manifest.json) and raw outputs remain unchanged. A small
solver-free check reproduced the issue and checked explicit restoration.
Continuation verified all prior inputs and output hashes and an identical full
effective-model snapshot, then solved only the remaining 1,033 genes. There
were **no repeated optimizations**, changed biological inputs or relaxed bounds.
The completed run restored reaction bounds/GPRs and verified unchanged input
files and loaded project source files.

A bounded independent read-only audit on 2026-10-08 checked input/code/output
identities, identical continuation conditions, all 1,074 unique journal calls,
table-to-journal correspondence, every threshold classification and all reported
metrics/coverage counts. All matched; the audit executed zero optimizations.
Across KOs, the largest mass-balance residual was 1.42114e-12 and the largest
bound violation 1.23554e-13. These are numerical checks, not biological validation.

The [runner](run_screen.py) reuses the project's medium/strain loader and reference
parser, applies COBRA's Boolean gene knockout, then maximizes growth once per
gene. A runnable `--self-check` covers exact-cutoff handling, unresolved solves
and orphan-state restoration without optimization. No FVA, dFBA, energy test,
model rebuild, parameter tuning, literature review or repair batch was added.
The next useful algorithmic improvement would be integrating the safe status
handling into the shared screening entry point; that separate code change is
not needed to interpret these recorded results and was not performed here.

## Deliverables

- [Full 1,073-gene result table](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/artifacts/tasks/screen_coq_r305_20261008/continuation/gene_screen.tsv>): raw growth/status, ratio, all four predictions, source putative functions, runtime flags and actual changed reaction bounds.
- [Experimental-positive coverage](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/artifacts/tasks/screen_coq_r305_20261008/continuation/experimental_coverage.tsv>).
- [Initial raw records](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/artifacts/tasks/screen_coq_r305_20261008/solver_calls.jsonl>) and [continuation records](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/artifacts/tasks/screen_coq_r305_20261008/continuation/solver_calls.jsonl>).
- [Scope and stop conditions](TASK.md), [complete provenance and metrics](manifest_continued.json).

Gene symbols and catalytic functions were not revalidated in this screen;
putative source descriptions and model reaction assignments are identified as
such in the table. No model export, commit or push was performed.
