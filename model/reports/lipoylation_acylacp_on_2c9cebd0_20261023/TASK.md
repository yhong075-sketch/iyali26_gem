# Lipoylation module + acyl-ACP chemistry on the main candidate 2c9cebd0 (2026-10-23)

## Authorization

User, 2026-10-23, after reading the plan in this session (message as received; it ends mid-sentence):

> Approved: add --lipoylation-module and --acyl-acp-chemistry on top of 2c9cebd0 as planned (steps
> after r2193, bring 9f2356b code/data/tests/reports without its XMLs, acyl test lipoylates 164f235f in
> memory, seeds 0/1 + control, whole-model diff, Gurobi 1e-9 two-model KO sweep with knocked-flux check

The plan it approves also listed an independent audit, REPORT.md, a STATE.md entry and "no commit or
push". Those steps are carried out; nothing is committed or pushed.

The two steps were themselves authorized earlier (recalled from `9f2356b`'s TASK files and STATE
entries, not re-asked): the lipoylation module on 2026-10-20, and the acyl-ACP chemistry correction,
including the resulting block of R2216–R2241 in every medium, on 2026-10-21. This task does not change
either decision; it re-applies them on a different starting model.

## Starting state

- **Branch:** `claude/loving-hermann-43c496` at `4c4248d` (= `origin/main`, PR #18 merged), clean.
- **Source of the two steps:** commit `9f2356b` on `claude/combined-candidate-20261019` (unmerged).
  Code, data and tests are copied unchanged from it, except where listed under Target.
- **Models** (none is a reference model):

| Role | File | SHA-256 |
| --- | --- | --- |
| Baseline (main candidate) | `E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_20261022.xml` | `2c9cebd0a5c53ec18bf1a583327099dc9e53fcf6b2fcbedd66cc7b9c54ea5a7b` |
| Common ancestor of both lines | `E5_coq9_resp_pkg_r132_combined_20261020.xml` | `164f235f4fe823181df0cc83082d6058ee9cf59675777adfe93557e9a77583a6` |
| Earlier result of the two steps (on `164f235f`; XML only on `9f2356b`) | `E5_coq9_resp_pkg_r132_combined_lipoylation_acylacp_20261020.xml` | `b3b9fe22cba4c8ae4408afe61740bfe75e661c383dd47a224965f06469618543` |

- **Guard check in the planning turn (in memory, verified):** the `9f2356b` lipoylation and acyl-ACP
  functions, applied in that order to the loaded `2c9cebd0` XML, both returned `applied`.
- **Order constraint (from the code):** the heme/iron guard pins the full `biomass_C` composition and
  note keys (`heme_iron_biomass.py`), so the lipoylation step, which adds four `biomass_C` terms and a
  note, must run after the heme/iron and R2193 steps. The lipoylation guard checks only its own four
  terms.
- **Before-state of every target:** [before.json](before.json), from
  [scripts/before_state.py](scripts/before_state.py).

## Target

1. Copy from `9f2356b`: `lipoylation_module.py`, `acyl_acp_chemistry.py`, `lipoylation_module.json`,
   `acyl_acp_chemistry.json`, both test files, and the report folders
   `lipoylation_module_20261020/` and `acyl_acp_block_review_20261020/` (evidence trail the data files
   cite). Not copied: its candidate XMLs, build JSONs, README rows and STATE entries.
2. `cli.py`: the two flags, as in `9f2356b`.
3. `main.py`: the two calls placed after the R2193 step (not after R1964 as in `9f2356b`); the
   requirement that `--acyl-acp-chemistry` needs `--lipoylation-module` is kept.
4. `test_acyl_acp_chemistry.py`: its baseline `fdd00ac3` is built in memory (`164f235f` + the
   lipoylation step) instead of read from the XML that is not on this branch.
5. A new test: both steps on `2c9cebd0` apply, change only their targets, and rerun as no-ops after a
   SBML round trip.
6. New candidate: the `2c9cebd0` argv plus `--lipoylation-module --acyl-acp-chemistry`, written to
   `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_20261023.xml`.

## Expected effect (inferred)

- **Diff vs `2c9cebd0`:** exactly the two steps' targets, the same as `b3b9fe22` vs `164f235f`:
  - lipoylation: new R2286–R2290, their new species and 3 genes, formulas of `m179[C_mi]` (ACP1) and
    `m586[C_mi]` (lipoyl-H-protein), R362/R364 without H+, 4 `biomass_C` terms and a note;
  - acyl-ACP: 32 formula fills, 19 reactions lose one H+, 7 cytosolic acyl-ACP hydrolases become
    forward-only, and notes.
- **WT growth:** slightly below `2c9cebd0` (`b3b9fe22` cost about 5.4e-6 relative to `a56f3649`).
- **Knockouts:** genes that change against `2c9cebd0` should be those of the `b3b9fe22` screen that
  come from these two steps: the 3 new lipoylation genes (YALI1C04991g, YALI1E13472g, YALI1F34729g) and
  the acyl-ACP/lipoyl subset of its 11 shared genes. YALI1E07744g (transketolase) changed through the
  R765/R766 edit already in `164f235f` and so should already be essential in `2c9cebd0`.

## Acceptance criteria

1. **Tests:** the full suite is identical test by test to the run before the change
   ([logs/pytest_before.log](logs/pytest_before.log)), plus the new tests, which pass.
2. **Builds:**
   - control (the `2c9cebd0` argv, new code) reproduces `2c9cebd0`;
   - the candidate is byte-identical under `PYTHONHASHSEED=0` and `1`.
3. **Whole-model diff vs `2c9cebd0`:** only the targets listed above; every touched reaction balanced
   or at its curated residual; the diff matches the `b3b9fe22` vs `164f235f` diff entity by entity.
4. **WT** (Task A medium + PO1f, Gurobi 1e-9) optimal on both models.
5. **Knockouts** on both models (Task A medium + PO1f):
   - Gurobi, FeasibilityTol and OptimalityTol 1e-9, 60 s per LP, model reset before each knockout;
   - targeted check: knocked-out reactions |v| ≤ 1e-9, other bound or |S v| violations ≤ 1e-6;
   - flagged points get a fresh re-solve (at most 20);
   - class changes at 1/5/10/15 % reported; scores per `model/expected/benchmark.md`.
6. Independent read-only audit, REPORT.md, a STATE.md entry and a candidates README row.

## Compute budget and stop condition

- **Budget:** 3 builds (control, seed 0, seed 1); about 6 diagnostic LPs; 2 × (about 1,080 + 1)
  knockout solves; at most 20 fresh re-solves; two full test-suite runs.
- **Stop and report on:** a guard rejection; a control that differs from `2c9cebd0`; seed-dependent
  bytes; a diff outside the targets; more than 20 flagged knockouts; a test regression; a knockout class
  change outside the expected set (reported, not explained away).

## Environment

- Python: main checkout's `.venv/bin/python` (→ `/opt/miniconda3/bin/python3`, Python 3.13.5,
  cobra 0.30.0, optlang 1.8.3); Gurobi (academic license).
- Research root `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`.
- Medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`), profile `po1f_sd_leu.json`
  (`35307853…`), experimental list `consensus_essential_genes.csv` (`1e887f5a…`).
