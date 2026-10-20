# R2041 forward-only bounds through the pipeline (2026-10-19)

## Authorization

The user, in this session, after reading
[r2041_direction_review_20261019/REPORT.md](../r2041_direction_review_20261019/REPORT.md):

> Approved: implement Option A (r2041_direction.json + guarded apply_r2041_direction, opt-in flag,
> keep 2 H+) in this worktree on 285eaee, build the new candidate, and run the 10 core LPs from REPORT
> 6.3.

Not authorized: commit, push, merge, the optional LPs P1–P2, the optional single-KO screen, any other
model change.

## Starting state

- Worktree `claude/bold-chaplygin-1c16df`, HEAD `285eaee` (clean except the untracked review folder
  `model/reports/r2041_direction_review_20261019/`).
- Baseline model: a56f3649 `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml`, SHA-256
  `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515`.
  - It is not a designated reference model.
  - Recorded build options: `--offline --no-solve --coq9-curation qcycle --energy-candidate E5
    --coq9-biomass-dilution --coq9-respiratory-package --r132-mito-adh`.
  - Recorded research root `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`
    (`IYALI26_RESEARCH_ROOT` is unset in this shell; passed with `--research-root`).
- Interpreter: main checkout `.venv`, Python 3.13.5, cobra 0.30.0, Gurobi (academic licence).
- Target state in a56f3649: [before.json](before.json).

## Target and change

- **Reaction:** R2041, ATP:D-arabinose 5-phosphotransferase, EC 2.7.1.54, C_cy, no GPR.
- **Change:** bounds [−1000, 1000] → [0, 1000]. The 2-H⁺ stoichiometry, species, GPR and annotation
  stay as they are. New notes record the curation and its evidence.
- **Mechanism:**
  - data file `model/curation/direction/r2041_direction.json`;
  - guarded `apply_r2041_direction` in `platform/scripts/gem_annotate/r2041_direction.py`;
  - opt-in flag `--r2041-forward`, applied after `apply_r1159_direction` in `main.py`;
  - focused test `platform/tests/test_r2041_direction.py`.
- **Not touched:** `fix_d_arabinokinase_direction_and_proton`,
  `metadata_reaction_selection.json`, and every other reaction, metabolite and gene.

## Expected effect

- Whole-model diff against a56f3649: only R2041's lower bound and new R2041 notes. In SBML this also
  changes the reversibility attribute and the lower-bound parameter.
- No growth change (inferred: all recorded R2041 fluxes are 0).
- The closed-exchange ATP maximum stays 0.

## Acceptance criteria

1. A control build without the flag reproduces `a56f3649…` byte for byte.
2. The candidate's whole-model diff against a56f3649 is limited to the target above.
3. R2041 stays element- and charge-balanced with the stored species.
4. The focused test passes. The full suite equals the before-run test by test, apart from the new
   tests.
5. The LPs below are optimal. Closed maxima are 0 on M0 and M1, and WT M0 = WT M1.
6. An independent read-only audit finds no undeclared change.

## Compute budget (10 LPs) and stop conditions

The context is the one in `coq9_respiratory_package_20261008/run_controls.py`:
- config from `model/reports/coq_candidate_checks_20261005/manifest.json`;
- PO1f SD-Leu overlay;
- Gurobi with Threads 1, Presolve 0, tolerances 1e-7, TimeLimit 60 s;
- closure by `tools.diagnose_closed_energy.close_model`.

M0 = a56f3649, M1 = the new candidate.

| # | LP | Model |
| --- | --- | --- |
| 1, 2 | WT growth | M0, M1 |
| 3, 4 | min / max R2041 with the objective ≥ 100 % of WT | M0 |
| 5, 6 | min / max R2041 with the objective ≥ 100 % of WT | M1 |
| 7, 8 | closed-exchange max ATP dissipation, C_cy (xMAINTENANCE) | M0, M1 |
| 9, 10 | closed-exchange max ATP dissipation, C_mi (balanced drain) | M0, M1 |

The runner asserts that no more than 10 optimizations are made.

Stop conditions:
- a closed objective > 1e-7: stop and report the witness fluxes;
- WT M0 ≠ M1 by more than 1e-9 relative: stop and report;
- a non-optimal status: record it raw. It is not biological evidence.

## Outputs

- This folder: TASK.md, before.json, REPORT.md, `scripts/`, `logs/`, small result JSON.
- Control builds: `.scratch/r2041_forward_20261019/` (git-ignored).
- Large outputs, if any: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r2041_forward_20261019/`.
- New candidate: `model/candidates/E5_coq9_resp_pkg_r132_r2041_forward_20261019.xml`, with its
  `.build.json`.
- A `model/STATE.md` entry. No commit or push.
