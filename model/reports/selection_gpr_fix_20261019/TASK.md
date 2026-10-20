# Fix three selection-overwritten reviewed GPR definitions (2026-10-19)

## Authorization (user, this session)

> Approved: implement A (R671 GPR re-apply as a new change set), B (drop R765/R766 from the
> builder's mis-annotation removals, no model change expected) and C option (i) (restore EG-GPR-001
> for R190/R159), as described in selection_gpr_review_20261019/REPORT.md section 5. No LP solves
> yet, no commit.

The approved plan is §5 of [the review](../selection_gpr_review_20261019/REPORT.md). The evidence
for each item is in §3 of that report.

## Where this is done

- Worktree `claude/clever-kare-7d2e1a`, HEAD `285eaee`.
- The 2026-10-19 restore that this change generalizes (SELECTION-GPR-RESTORE-20261019) is
  uncommitted in worktree `claude/orphan-reaction-newest-model-b8ddc8`. Its five code, data and
  test files were copied here unchanged before any edit, and their SHA-256 values were checked equal:

  | File | SHA-256 |
  | --- | --- |
  | `platform/scripts/gem_annotate/cli.py` | `a63582a6fb296d1ec3c039e50713b4199bb80eb1616f3d0b055dc8db02adb003` |
  | `platform/scripts/gem_annotate/main.py` | `06b9dcb697b98be4d3c07e729f1d3f419d369aab7ddc847b2053f3a2a8c39d63` |
  | `platform/scripts/gem_annotate/selection_gpr_restore.py` | `0a8c4dfd3a2f6fbc77c6260dea1031ff1e847825f8a6e40b2905b3c65dc66e1e` |
  | `platform/tests/test_selection_gpr_restore.py` | `dac2d38cff63c396ae2a3abf68e2eda7562d17949cbf608a6a343e528933b1fa` |
  | `model/curation/gpr/selection_cleared_gpr_restore.json` | `d048ff9046584c75b54cfdae438874e448c3c883e92d1ce994c1f4eec98aba95` |

  Before the copy, `cli.py` and `main.py` at HEAD were byte-identical in the two worktrees. The
  other worktree was not modified.

## Baseline

- **Model:** `E5_coq9_resp_pkg_r132_gpr_restore_v2_20261019.xml`, SHA-256
  `c98961e960c6a0f8283f485e32daf20b04dd70d998fad86fc09dca0d14ccae4c`, in the other worktree. It is the
  newest candidate, not a designated reference.
- **Its parent:** `E5_coq9_resp_pkg_r132_20261012.xml`, `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515`.
- **Target state:** [before.json](before.json). R671, R159, R190, R607, R765 and R766 are identical in
  a56f3649 and c98961e9.
- **Build flags of c98961e9:** `--offline --no-solve --coq9-curation qcycle --energy-candidate E5
  --coq9-biomass-dilution --coq9-respiratory-package --r132-mito-adh --restore-selection-cleared-gprs`.
  Research root: `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`.

## Targets and intended change

| Item | Reaction | Now (selection "after") | Target |
| --- | --- | --- | --- |
| A | R671 | GPR `YALI1F23441g or YALI1E11370g`; equation unchanged (mass-balanced, charge +1) | GPR `YALI1F23441g`; same equation |
| C(i) | R159 | ATCase: carbamoyl-P + Asp → carbamoyl-Asp + Pi; GPR `YALI1C33005g or YALI1D09420g or YALI1E11768g` | net channelled CPS + ATCase (selection "before", balanced); GPR `YALI1E11768g` |
| C(i) | R190 | GPR `YALI1C33005g or YALI1D09420g or YALI1E11768g` | GPR `YALI1C33005g and YALI1D09420g`; equation unchanged (balanced) |
| B | R765/R766 | builder removes YALI1E07744g before the selection; the selection adds it back | builder no longer removes it; final GPR stays `YALI1E07744g or YALI1D02625g` |

Mechanism:
- **A and C:** new curated change sets applied after the metadata selection by a generalized
  `selection_gpr_restore.py`. Each has its own opt-in flag: `--restore-selection-overwritten-gprs`
  for A, `--restore-selection-cpa-ura2` for C.
- **B:** a builder correction in `patches.py`.

The following do not change:
- the selection record `metadata_reaction_selection.json`;
- the SELECTION-GPR-RESTORE-20261019 data and its behaviour;
- genes, metabolites, formulas, charges, bounds;
- names, EC numbers and subsystems.

**Expected effect (inferred; no solve this round):**
- YALI1E11370g no longer reaches any reaction.
- A YALI1F23441g knockout closes R671.
- R190 needs both CPA genes.
- R159 needs URA2 and no longer draws on the shared carbamoyl-phosphate pool. m325 is then made only
  by R190 and used only by R607.

## Acceptance criteria

1. Control builds:
   - without any new flag, the build reproduces a56f3649 byte for byte (this checks that B changes
     no model);
   - with only `--restore-selection-cleared-gprs`, it reproduces c98961e9.
2. The candidate build (c98961e9's flags plus both new flags) differs from c98961e9 only in:
   - R671, R159 and R190 (GPR, R159 equation, restore notes);
   - the derived gene–reaction links of YALI1E11370g, YALI1C33005g, YALI1D09420g and YALI1E11768g;
   - the metabolite links of the R159 species.
3. R159 and R190 balance in mass and charge. R671 keeps exactly its recorded +1 charge residual.
4. Guards, shown by tests:
   - stale or unexpected state, species, gene or context raises with no edits;
   - a rerun is a no-op, also after an SBML save and reload;
   - each flag is opt-in and needs a complete selection;
   - the 2026-10-19 change set behaves as before.
5. `python -m pytest` matches the pre-change run test by test, plus the new tests passing.
6. No LP solve this round.

## Compute budget and stop condition

- **Builds:** 3, offline and no solve:
  - 2 controls to `.scratch/selection_gpr_fix_20261019/`;
  - 1 candidate to `model/candidates/`.
- **Tests:** two full pytest runs, before and after the change.
- **No LP solves, no FVA**, no cluster jobs, no commit or push.
- **Stop:**
  - a control build does not reproduce its reference;
  - any non-target difference;
  - a balance failure;
  - a test regression.

## Corrections after the independent audit (2026-10-19)

The text above is kept as written. These corrections apply:

- **Provenance text in the R671 note.** The first candidate's R671 evidence said E11370g entered R671
  "on R671's wrong source EC 6.3.5.7". The source model has no EC 6.3.5.7 (`grep -c` = 0; its R671
  carries only PROTEIN_CLASS 1.3.1.13). 6.3.5.7 is in the isozyme table's R671 row and in the legacy
  `model.xml`. Notes are part of the model, so the data was corrected and the candidate rebuilt as
  v2. The first candidate (`fa859429…`) is superseded and kept until the user decides.
- **m325 pool guard.** The CPA/URA2 change set now also checks that, apart from its targets, only
  R607 uses m325 (`pool_context`). Before, only R607's own definition was checked.
- **Tests** now assert which guard fires for each stale case, and add an extra-pool-user case.
- **Budget.** The fix needed 3 more offline, no-solve builds (both controls again with the final code,
  plus candidate v2) and one more full pytest run. Still no LP.

## Solve authorization (user, 2026-10-19, after the build)

> Authorized: run 12 LPs with GLPK on SD-Leu + PO1f — WT and single-gene KO of YALI1F23441g,
> YALI1E11370g, YALI1C33005g, YALI1D09420g, YALI1E11768g on c98961e9 and on
> E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_20261019.xml. No other solves, no commit.

The condition is the one used by the 2026-10-19 orphan and restore reports:
- `model/reports/reference_pipeline_restore_20260909/research/state/media/sd_leu.csv`
  (`ed176d26a373f98cc413ed2e32a71f5f060a06e343f90f7db25cd32eff268e85`), whole medium replaced;
- `.../strain_profiles/po1f_sd_leu.json`
  (`35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383`), R612 closed and the uracil
  uptake override applied.

GLPK with presolve. A non-optimal status is recorded as a status, not as growth 0. The script is
`scripts/solves.py`, the 2026-10-19 script with only the gene list changed.
