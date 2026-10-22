# AUDIT: Thiamine proton lock (THIAMINE-PROTON-LOCK-20261022)

**Date:** 2026-10-22. **Auditor:** an independent, read-only subagent with no conversation history.
This is its report as returned, lightly reformatted. It audited the first candidate, `ada4b7f1`; the
responses and the rebuild are in [REPORT.md §7](REPORT.md).

**What it did:**
- It changed no repository or research-root file (`git status` was identical before and after).
- It ran no build and solved no model LP.

**What it read:**
- TASK.md, AGENTS.md, before.json, after.json, model_diff.json and commands.txt;
- the logs, scripts, curation JSON, step, tests and gate;
- `main.py` and `reaction_selection.py`;
- `dipeptide_chemistry.protected_chemistry` and `patches.fix_d_arabinokinase_direction_and_proton`;
- `tools/update_model.py:303`;
- every build.json;
- the draft REPORT.md.

**What it ran:**
- Hashes of all models, builds and code.
- In-memory gate runs on fresh reads of `ef2265c1` and `ada4b7f1`, with both the HEAD gate (`cb3cd03f`)
  and the current gate (`fc4ddff0`).
- A dict-lock SBML roundtrip.
- 20 step and gate probes.
- A raw XML diff.
- A key-form grep across every model XML.
- MetaNetX 4.5 lookups and a balance recomputation.
- The focused pytest plus `test_dipeptide_chemistry` (11 passed).

## Claims

| # | Verdict | Evidence |
|---|---|---|
| C1 | **supported** | The HEAD gate on `ef2265c1` removes the H+ from R752 and R408 (m10) and R750 (m342). Each residual is exactly {H: 1, charge: 1} before, and 0 after. 36 changes in total. |
| C2 | **supported** | A dict lock on R752 is skipped in memory. After a roundtrip it is the string `{&apos;lock_proton_water_stoichiometry&apos;: True}` and is stripped. In `ef2265c1` no reaction has `curated_reaction_correction` or the lock key. R2041 has only the string `metadata_previous_curated_reaction_correction`. |
| C3 | **supported** | The diff of HEAD against current `microspecies.py` is the docstring plus 4 lines. The dict branch is untouched, and on the whole of `ef2265c1` the old and new gate reports are identical. Nuance: a Python `True`, `"true\n"` and `"tRuE"` also lock. Harmless. |
| C4 | **supported** | Default before = after = `b4ce0974…`; control = `ef2265c1…`; N = seed1 = `ada4b7f1…`. The build.json `source_sha256` values match the before and after hashes of `microspecies.py`. **Caveat:** PYTHONHASHSEED is recorded only in the log headers. |
| C5 | **supported** | The raw diff is 28 added lines and 0 removed, all notes in R_R408, R_R750 and R_R752. |
| C6 | **supported** | The current gate on N skips all three and leaves them unchanged. Its 33 other changes, rejected list (317) and skipped-missing-formula count (427) equal `ef2265c1`'s. The HEAD gate on N still strips all three. |
| C7 | **supported, with gaps** | Every stale-input probe failed closed with no edits, and re-apply gives `already_correct`. **Gaps:** some failures raise AttributeError or KeyError, not ValueError. Note texts are not validated at runtime (c1). |
| C8 | **supported** | With MetaNetX 4.5 partner forms, all three balance with the H+ included. That MNXref is at pH 7.3 is **recalled**; the file header gives no pH. |
| C9 | **supported** | The only per-test difference is 8 new passes: 232 → 240 passed, 591 → 608 subtests; failed, errors and skipped unchanged. |

## (a) Overstatements

1. **The rationale in the model notes, the JSON and TASK.md is wrong.** It says the H+ "is off only
   because ATP, ADP, AMP and phosphate are still neutral".
   - **Per-species contributions** (model H − pH 7.3 H, times the coefficient):
     - R752: ATP −4, AMP +2, TPP +3.
     - R408: ATP −4, ADP +3, HET-P +2.
     - R750: Pi +4, TPP −3.
   - Moving only the ATP family to pH 7.3 would leave +3 / +2 / −3 H+.
   - The KEGG-form TPP and HET-P count too.
2. **"…including when a build starts from a saved model" is conditional.** It holds only with this
   worktree's gate; the HEAD/main gate strips all three H+ from N. No `--starting-model` build was run.
3. **TASK.md says the gate honours "only" the dict lock.** It also skips reactions next to
   dipeptide-chemistry-locked metabolites, a JSON string lock that survives SBML. That lock is
   metabolite-wide, so it is not a good fit here.
4. **Line numbers are stale.** After this task's wiring the gate calls are at `main.py` 253, 344 and 611,
   and the steps at 766–780.
5. **REPORT draft:**
   - "No LP was run": the full pytest run includes tests that solve LPs.
   - The R304 libSBML error appears once each in the control and seed logs, and not in the default-build
     logs.

## (b) Outside the authorized scope

- No stoichiometry, bound, GPR or annotation change, and no commit.
- There were 5 builds, as budgeted.
- The "0 LPs" budget against the full pytest run: within the approved plan, but not literally zero LPs.
- The gate change is global for any reaction carrying the key. That was authorized.

## (c) Code defects and side-effect checks

1. **`thiamine_proton_lock.py` `_check_spec`: note texts are not checked at runtime.**
   - Only the generator asserts SBML safety.
   - Reproduced: a rationale containing `pH < 7.3 & neutral` was applied. `write_deterministic_sbml_model`
     logged a libSBML error but did not raise. The read-back R752 had no notes, so the lock was lost and
     the gate stripped its H+.
   - It cannot happen with the current JSON (SHA `e5c92fb8`).
2. **`_unescaped`:** a non-string note value (a dict lock, a bool) raises AttributeError. It still fails
   closed.
3. **`microspecies.py` 1111–1112:**
   - any value whose `str()` is "true" locks;
   - the key is bare, with no curation ID required;
   - the skipped list mixes three reasons.
   - No functional defect.
4. **Effect on the default build and existing candidates: none.**
   - The key form appears only in N and its seed-1 copy, across all model XML files in the worktree and
     the research root.
   - There is one writer (this step) and one reader (the gate). In a source build the setter runs after
     the last gate call.
   - The only production callers of the gate are at `main.py` 253, 344 and 611.

## (d) Relying on the lock

1. It works only with the gate at `fc4ddff0`. Older gates strip the H+ silently.
2. **Other stoichiometry editors ignore it:**
   - `tools/update_model.py:303`, an H/O balancer that ignores charge (read, not run);
   - the metadata selection (R752, R408 and R750 are not in its spec);
   - hydroxide normalisation, the patches, the energy candidates and the duplicate merge.
3. **Once saved, the lock is unconditional.** After a partial migration (only ATP, ADP, AMP and
   phosphate) the reactions would sit at +3 / +2 / −3 H+, listed only as skipped. A check that a locked
   reaction's residual is {} or the reviewed one would catch that.
4. N still has the three reactions off by one H+, as intended. Balance audits will still flag them.
5. **Silent loss of notes on SBML write is a general risk.** It was seen on R304, and reproduced for
   this lock.
6. R2041's lock is still lost in saved models (out of scope).
7. A rebuild from N is untested.

**Partly checked:** no rebuild, no full pytest re-run, self-reported seeds, MNXref pH recalled,
`update_model.py` not run, REPORT audited only as a draft.
