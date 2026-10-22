# Thiamine proton lock: keep the pH 7.3 H+ on R752, R408 and R750 (2026-10-22)

## Authorization

User, 2026-10-22, after reading the plan in this session:

> Approved: implement the H+ lock plan as written: gate honours a string note
> lock_proton_water_stoichiometry = true, new opt-in --thiamine-proton-lock on top of ef2265c1 for
> R752/R408/R750, default-build and control checks, no LPs. No commit or push.

This follows the [thiamine phosphate chemistry task](../thiamine_phosphate_chemistry_20261022/REPORT.md),
whose audit found that the proton gate would strip the pH 7.3-correct H+ from these three reactions.

## Starting state

- **Branch:** `claude/science-test-gap-fix-13f440` at `abe7e40`. The TPP and thiamine phosphate tasks
  are uncommitted in the same worktree.
- **Baseline model:**
  `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_v2_20261022.xml`, SHA-256
  `ef2265c1fd6c3cd12eebacc8c03901091f90289a6f77a7be66d2ede9d7379440`. It is a candidate, not a reference.
- **Before-state:** [before.json](before.json), from [scripts/state.py](scripts/state.py). It records the
  three targets, any existing lock, and a whole-model dry run of the gate on the saved file.

## What was verified before any change (2026-10-22)

**The gate strips the H+ (dry run on an in-memory copy of `ef2265c1`).**
- `balance_protons_and_water` removes one H+ each from R752 (m10), R408 (m10) and R750 (m342).
- Their residual is exactly {H: 1, charge: 1}, which passes the gate's condition `H - 2*O == charge`.
- The gate would also change 33 other reactions on the saved file. So the saved candidate is not a
  fixed point of the gate. Those 33 are outside this task.

**The existing lock does not survive saving.**
- The gate honours `notes["curated_reaction_correction"]` when it is a dict holding
  `lock_proton_water_stoichiometry: True`.
- It also skips reactions next to dipeptide-chemistry-locked metabolites. That lock is a JSON string
  note, so it does survive SBML, but it is metabolite-wide and pins species fields, so it does not fit
  here. (The first version said "only" the dict lock; corrected after the audit.)
- A dict lock on R752 works in memory.
- After an SBML write/read, the note is the string `{&apos;lock_proton_water_stoichiometry&apos;: True, …}`.
  The gate ignores it and strips R752 again.
- R2041's dict lock, set in the build by `fix_d_arabinokinase_direction_and_proton`, appears in the
  saved file only as the string note `metadata_previous_curated_reaction_correction`, renamed by the
  metadata selection. So it no longer locks anything.
- No reaction in `ef2265c1` carries a `lock_proton_water_stoichiometry` note key.

**The gate runs only before the opt-in steps.**
- `main.py` lines 250, 341 and 608 call it (253, 344 and 611 after this task's wiring), and the
  thiamine steps run after line 760.
- It would also run on this candidate if a build started from it (`--starting-model`).

## Targets

**1. The gate** (`platform/scripts/gem_annotate/microspecies.py`, `balance_protons_and_water`).
- It also skips a reaction whose note `lock_proton_water_stoichiometry` is the string `true` (case- and
  whitespace-insensitive).
- Skipped reactions are listed in `skipped_curated_lock_reaction_ids`, as the dict lock's are.
- The dict form is unchanged.

**2. The reactions.** R752, R408 and R750 gain the note `lock_proton_water_stoichiometry: true` plus
`thiamine_proton_lock_*` evidence notes. No stoichiometry, bound, GPR or annotation changes.
- Curated data: `model/curation/chemistry/thiamine_proton_lock.json`.
- Guarded step: `apply_thiamine_proton_lock()` in `platform/scripts/gem_annotate/thiamine_proton_lock.py`,
  under the opt-in flag `--thiamine-proton-lock`. The flag requires `--thiamine-phosphate-chemistry` and
  runs right after it.
- Focused tests: `platform/tests/test_thiamine_proton_lock.py`, covering the gate change and the step.

## Evidence

- **The H+ is correct at pH 7.3:** with every partner at its MetaNetX 4.5 pH 7.3 form, all nine thiamine
  reactions balance as written. I verified this in the previous task, and its audit recomputed it.
- **It is off only because of the mixed convention:** the residual {H: 1, charge: 1} comes from the
  mix of forms.
  - The thiamine phosphates and HET-P are at their KEGG fully protonated forms.
  - ATP, ADP, AMP and phosphate are neutral.
  - Per the audit, the per-species contributions are R752: ATP −4, AMP +2, TPP +3; R408: ATP −4,
    ADP +3, HET-P +2; R750: Pi +4, TPP −3.
  - Only moving all of them to pH 7.3 together removes it.
  - The first version blamed only ATP, ADP, AMP and phosphate; corrected after the audit.

## Expected effect

- In a build, the gate never runs after this step, so the step itself changes only the notes.
- The gate on the saved new candidate leaves R752, R408 and R750 unchanged and reports them as locked.
- **No LP change**, since only notes change.

## Acceptance criteria

1. **The default build** (`--offline --no-solve --coq9-curation metadata`) is byte-identical before and
   after the gate edit.
2. **The ef2265c1 argv without the new flag**, built after the gate edit, reproduces `ef2265c1…` byte for
   byte.
3. **With the flag,** `PYTHONHASHSEED=0` and `=1` give the same bytes.
4. **The diff against `ef2265c1`** holds only the notes of R752, R408 and R750. Anything else stops the
   work.
5. **The gate on the re-read new candidate** changes none of the three reactions and lists them in
   `skipped_curated_lock_reaction_ids`. Its other changes equal the 33 on `ef2265c1`.
6. **pytest per test** equals the baseline plus the new tests passing. The baseline is the thiamine
   task's final run on identical code.

## Compute budget and stop condition

- **Budget:**
  - 0 LPs;
  - 5 builds: default before and after, control, and seeds 0 and 1;
  - in-memory gate dry runs.
- **Stop and report to the user before continuing on:**
  - a guard rejection;
  - a default-build or control byte difference;
  - a diff outside the targets;
  - seed-dependent bytes;
  - a test regression;
  - a gate dry run that still changes a locked reaction.

## Not in scope

- Restoring R2041's lost lock.
- The 33 other reactions the gate would change on the saved candidate.
- Locking any other reaction.
- Changing the dict-lock mechanism.
