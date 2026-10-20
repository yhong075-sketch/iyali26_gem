# Selection restore: R159 species order depended on PYTHONHASHSEED (2026-10-19)

Status: **fixed and verified in this worktree, not committed.** The fix changes the builder's output only
byte-wise: no model science, curation data, medium or benchmark changed. Scope and authorization:
[TASK.md](TASK.md).

## Result

| Item | Value |
| --- | --- |
| Code base | `dc471d75480bd688280df0bed5af2db8bbd7cebb`; this worktree's branch `claude/great-wu-64e590` was fast-forwarded to it from `285eaee` |
| Code change | `platform/scripts/gem_annotate/selection_gpr_restore.py`: line 154 iterates `sorted(...)`, plus one comment line. SHA-256 `9bf0f5a6…` → `27cc380d1b3f8a5de5e919633fa26e36fcbde7d3488bd73978c890aa201d255a` |
| Test added | `platform/tests/test_selection_gpr_restore_change_sets.py::ChangeSetRestoreTests::test_restored_species_order_does_not_depend_on_the_hash_seed` |
| New candidate | `model/candidates/E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_sorted_20261019.xml`, SHA-256 `f0f24b260c977f8719556d5b0ee4416ab7a4ac273bd59f50433c0e6eccbd0380` (with `.build.json` and `.coq9_genes.tsv`) |
| Byte reference (not overwritten) | `E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_20261019.xml`, SHA-256 `51c6fabd486f4e1a9d6365bf7740bde3eca20fa05f1d5454a73586d3276db3fa` |
| Environment | main checkout `.venv`: Python 3.13.5, cobra 0.30.0; `IYALI26_RESEARCH_ROOT=/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research` |

## Cause

Code logic (statically verified by reading dc471d7):

1. `apply_selection_gpr_restore` built the edit as
   `delta = {mid: ... for mid in set(current) | set(target["stoichiometry"])}`. Iterating a set of `str`
   gives an order that changes with `PYTHONHASHSEED`.
2. cobra 0.30 `Reaction.add_metabolites(..., combine=True)` updates a species already in the reaction
   in place, drops species whose coefficient becomes 0, and appends new species in the order of the
   dict it is given.
3. cobra 0.30 `write_sbml_model` writes `reaction.metabolites` in insertion order
   (`cobra/io/sbml.py:1367`), split into reactants and products. `write_deterministic_sbml_model` sorts
   only group members.

So a restored reaction's XML changed with the seed only when the restore **adds two or more species**.
At dc471d7 only R159 does that: it adds m10, m130, m141, m143, m199, m32 and m50, and drops m325.
R2274 adds one species, R1303 drops one, and R612, R1302, R671 and R190 keep their equation.

Observed before the fix ([scripts/seed_order.py](scripts/seed_order.py),
[logs/seed_order_dc471d7.txt](logs/seed_order_dc471d7.txt)): the dc471d7 restore on the in-memory test
fixtures, under `PYTHONHASHSEED` 0 to 7, gave 8 distinct R159 orders in 8 runs and a single order for
every other target. The real files agree: R159 reactants were `m267, m199, m130, m32, m10, m141` in
51c6fabd and `m267, m199, m141, m10, m130, m32` in the R153 session's rebuild `8f083627…`.

## Fix

```python
# Sorted: add_metabolites appends new species in this order and the SBML writer keeps it.
delta = {mid: target["stoichiometry"].get(mid, 0.0) - current.get(mid, 0.0)
         for mid in sorted(set(current) | set(target["stoichiometry"]))}
```

Species already in the reaction keep their position, and added species follow in sorted ID order.
R159 now reads reactants `m267, m10, m130, m141, m199, m32` and products `m35, m326, m143, m50`, which
is the order predicted before the build. The new test runs all three change sets on the existing
fixtures in child processes with `PYTHONHASHSEED` 0, 1 and 2, and asserts two things:

- every target's species order is the same across the three seeds;
- that order is "kept species in prior order, then added species sorted".

## Verification

| Check | Result | Level |
| --- | --- | --- |
| New test on the unfixed code ([log](logs/new_test_before_fix.txt)) | Fails: R159 order differs between seeds | observed in a run |
| Both selection-restore test files after the fix ([log](logs/focused_tests_after_fix.txt)) | 17 passed (16 existing + 1 new) | observed in a run |
| Full pytest, before ([txt](logs/pytest_before.txt), [xml](logs/pytest_before.xml)) and after ([txt](logs/pytest_after.txt), [xml](logs/pytest_after.xml)); [comparison](logs/pytest_comparison.txt) | Before: 181 passed, 11 failed, 14 errors, 31 skipped. After: 182 passed, same failures, errors and skips. Per test, the only change is the new test (absent → passed); the failures and errors predate this change (same list in both runs) | observed in a run |
| Two builds with the dc471d7 argv, `PYTHONHASHSEED=0` (candidate) and `=1` (`.scratch/selection_restore_order_20261019/seed1.xml`); logs [seed 0](logs/build_seed0.log), [seed 1](logs/build_seed1.log) | Both `f0f24b26…`, byte-identical. `.coq9_genes.tsv` is identical in both and also in 51c6fabd's (`2dd602b7…`). The build records differ only in time, argv, output path and `dirty_status` | observed in a run |
| `model_diff.py` new vs 51c6fabd ([json](logs/model_diff_vs_51c6fabd.json)) | total differences 0 | observed in a run |
| Byte diff new vs 51c6fabd ([scripts/byte_diff.py](scripts/byte_diff.py), [output](logs/byte_diff_vs_51c6fabd.txt)) | All text outside reaction blocks identical; 2306 blocks in the same order; only `R_R159` differs, and its lines are identical once sorted | observed in a run |
| Build inputs vs 51c6fabd's `.build.json` | Input model SHA `5c8c199e…` equal; all 99 data-file SHAs equal; among the recorded source files only `selection_gpr_restore.py` differs | observed in a run |

Build command, run from `platform/` (seed 1 differs only in the seed and `--output-model`):

```bash
PYTHONHASHSEED=0 IYALI26_RESEARCH_ROOT="$RR" python -m scripts.gem_annotate --research-root "$RR" --offline --no-solve --coq9-curation qcycle --energy-candidate E5 --coq9-biomass-dilution --coq9-respiratory-package --r132-mito-adh --restore-selection-cleared-gprs --restore-selection-overwritten-gprs --restore-selection-cpa-ura2 --output-model ../model/candidates/E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_sorted_20261019.xml
```

Each build took about 70 s. No optimization was run.

## What this does and does not show

- **Provable for this step:** the order of `sorted()` over strings does not depend on the hash seed, so
  this step's insertion order is now fixed by the species IDs.
- **Observed for the whole build, seeds 0 and 1 only:** with this argv, no other builder step's
  seed-dependence reaches the XML. That is evidence, not proof for all seeds or other flag sets. The
  builder was not audited for other set iterations.
- **Seed-dependent output that remains:** the order of cobra's formula warnings in the build log. Both
  logs hold the same 98 warnings in a different order. The XML and the gene table are unaffected.
- **The build is from a dirty tree:** the uncommitted fix and test. `dirty_status` in the
  `.build.json` records this.
- **Alternative considered:** sort species references in `write_deterministic_sbml_model`. It would
  remove every order dependence of this kind, but it would change the bytes of every reaction, so no
  later build could match an existing candidate byte for byte. It was not done.
- **Affects another line:**
  - The R153 GPR revision (`claude/eloquent-borg-1bd1b0`) builds on dc471d7 without this fix and pins
    `PYTHONHASHSEED=0`.
  - Its files and this candidate will differ byte-wise at R159 even when the models are the same.
  - Once both lines are merged, a rebuild carries the sorted R159 order.
