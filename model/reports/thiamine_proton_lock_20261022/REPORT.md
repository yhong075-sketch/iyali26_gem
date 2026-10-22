# Thiamine proton lock: implementation report

Scope, authorization, budget and acceptance criteria are in [TASK.md](TASK.md). Exact commands are in
[commands.txt](commands.txt). The independent audit is in [AUDIT.md](AUDIT.md); the responses are in
§7.

Evidence tiers used below:
- **verified**: run or checked in this session; the output is cited.
- **provable**: follows from verified facts.
- **inferred**: reasoned, with no direct check.
- **recalled**: from memory, not opened this session.

## 1. Summary

- **New candidate (verified):**
  - `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_v2_20261022.xml`,
    SHA-256 `2fe51cc2b3834bb7303d6702f4f68ac7bec2dd5353cafcd635803fd6e74c174e`.
  - It is `ef2265c1` (the thiamine phosphate chemistry candidate) plus the opt-in step
    `--thiamine-proton-lock`.
  - It is a candidate, not a reference, and it is uncommitted.
  - The first build `…_lock_20261022.xml` (`ada4b7f1`) is superseded. Only its rationale note differs;
    the audit found that note wrong. It is kept until the user decides.
- **What changed (verified):**
  - **The proton/water gate** (`balance_protons_and_water` in `microspecies.py`) now also skips a
    reaction whose note `lock_proton_water_stoichiometry` is the string `true`. That form survives SBML;
    the old dict lock does not. This is 4 lines of code plus a docstring sentence.
  - **R752, R408 and R750** carry that note, plus `thiamine_proton_lock_*` evidence notes.
  - Nothing else in the model changed: 3 differences, all notes.
- **Effect (verified):**
  - On `ef2265c1`, the gate would remove the pH 7.3-correct H+ from all three reactions.
  - On the saved new candidate, read back from file, it leaves them alone and lists them as locked. Its
    other 33 changes, its rejected list and its skip counts are unchanged.
- **No side effects (verified):**
  - The default build is byte-identical before and after the gate edit (`b4ce0974…`).
  - The `ef2265c1` argv rebuilt after the edit reproduces `ef2265c1…`.
  - No candidate or default model carries the new note key except this candidate and its seed copy (the
    audit's scan).
- **Tests and audit:**
  - The full suite differs only by the new passing tests (§6).
  - The audit supported all 9 claims.
  - It found a wrong rationale in the model notes and two code gaps, now fixed (§7).

## 2. The change

| Item | Location |
| --- | --- |
| Gate | `platform/scripts/gem_annotate/microspecies.py`, `balance_protons_and_water`: `cb3cd03f…` (HEAD) → `fc4ddff0…` |
| Curated data | [`model/curation/chemistry/thiamine_proton_lock.json`](../../curation/chemistry/thiamine_proton_lock.json) (`THIAMINE-PROTON-LOCK-20261022`, `2caf1981…`), generated from [before.json](before.json) by [scripts/make_curation_json.py](scripts/make_curation_json.py) |
| Guarded step | `platform/scripts/gem_annotate/thiamine_proton_lock.py` (`cfe8543c…`), `apply_thiamine_proton_lock` |
| Wiring | `cli.py` flag `--thiamine-proton-lock`; `main.py` calls the step right after `--thiamine-phosphate-chemistry`, which it requires |
| Focused tests | `platform/tests/test_thiamine_proton_lock.py` (the gate change and the step) |
| Code diff | [logs/code_diff.patch](logs/code_diff.patch), which also holds the earlier tasks' wiring in cli.py and main.py |

**The gate change**, right after the existing dict-lock check:

```python
plain_lock = (reaction.notes or {}).get("lock_proton_water_stoichiometry", "")
if str(plain_lock).strip().lower() == "true":
    skipped_curated_locks.append(reaction.id)
    continue
```

- The values `true`, ` True ` and `TRUE` lock; `false`, `yes` and an empty value do not (tested).
- The audit noted that a Python `True` and `"true\n"` also lock. That is harmless.

**What the step checks before writing.** It fails closed and makes no partial edits.
- **The spec:**
  - schema version (an integer 1), ID, status and the exact reaction set;
  - the lock note must be exactly `true`, and each reviewed residual exactly {H: 1, charge: 1};
  - every note text must be a string without `< > & ' "`. libSBML silently drops a whole notes block
    holding `<` or `&`, which would drop the lock on save.
- **Each reaction:** name, bounds, GPR and stoichiometry (with its H+).
- **The residual:** it must be exactly {H: 1, charge: 1}.
- **Notes:**
  - existing notes are compared after XML-unescaping, and non-string values as they are;
  - a dict or bool lock already present, a conflicting value, or the lock without this curation's
    evidence notes stops the step.
- **Completeness:** either all three reactions carry the lock or none does.

A rerun returns `already_correct`.

## 3. Mechanism and evidence

**Why lock (verified, before the change):**
- R752, R408 and R750 are each off by exactly one H+. That passes the gate's condition
  `H - 2*O == charge`, so the gate removes the H+ ([before.json](before.json) `gate_dry_run`).
- **Where the residual comes from** (audit, recomputed): the mix of forms. The thiamine phosphates and
  HET-P are at their KEGG fully protonated forms, while ATP, ADP, AMP and phosphate are neutral. Per
  species:
  - R752: ATP −4, AMP +2, TPP +3;
  - R408: ATP −4, ADP +3, HET-P +2;
  - R750: Pi +4, TPP −3.
- Only moving all of them to pH 7.3 together removes it. The first rationale note blamed only the
  neutral ATP family. That was wrong and is corrected (§7).

**Why the H+ is right (verified in the previous task, and recomputed by both audits):** with every
partner at its MetaNetX 4.5 form, all three balance as written:
- R752: thiamine(1+) + ATP(4−) = TPP(2−) + AMP(2−) + H+;
- R408: HET + ATP(4−) = HET-P(2−) + ADP(3−) + H+;
- R750: TPP(2−) + 2 H2O = thiamine(1+) + 2 HPO4(2−) + H+.

That MNXref uses a pH 7.3 convention is **recalled** (the local file header gives no pH). ChEBI does
label the charges used: for example CHEBI:58937 is thiamine(1+) diphosphate(3−).

**When the gate runs (verified):**
- `main.py` lines 253, 344 and 611 after this task's wiring, all before the opt-in steps. So in today's
  build the lock changes nothing.
- It protects later steps that call the gate, and builds that start from a saved candidate
  (`--starting-model`). That second case is **inferred** from the dry run on the saved file; no such
  build was run.

**Why a new lock form (verified):**
- The dict lock works only in memory. An SBML write/read turns it into a string, which the gate
  ignores. The plain note reads back as `true`.
- The gate's third skip, for reactions next to dipeptide-chemistry-locked metabolites, survives SBML.
  But it is metabolite-wide and pins species fields, so it does not fit here.

## 4. Builds and whole-model verification

| Build | Result |
| --- | --- |
| Default, before the gate edit (microspecies.py `cb3cd03f…`) | `b4ce0974…`, the default hash already recorded in the candidates README |
| Default, after the gate edit (`fc4ddff0…`) | `b4ce0974…`, byte-identical |
| `ef2265c1` argv without the new flag, after the edit | `ef2265c1…`, byte-identical |
| First candidate, seeds 0 and 1 | `ada4b7f1…` (audited) |
| Final candidate, seeds 0 and 1, after the audit fixes | `2fe51cc2…`, byte-identical |

- **Hashes and seeds:** [logs/sha256.txt](logs/sha256.txt) and [logs/sha256_final.txt](logs/sha256_final.txt).
  Each log's first line records its seed; build.json does not.
- **Where the builds are:** the non-candidate builds, and `ada4b7f1`, are in
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/thiamine_proton_lock_20261022/`.
- **The gate code is fixed:** `microspecies.py` is `fc4ddff0` for the default-after, control and final
  builds.
- **Diff against `ef2265c1` ([model_diff.json](model_diff.json)):** 3 differences, the notes of R408,
  R750 and R752. Each gains `lock_proton_water_stoichiometry` and the 7 `thiamine_proton_lock_*` keys.
  No note was removed or changed.
- **Diff against `ada4b7f1`
  ([model_diff_final_vs_ada4b7f1.json](model_diff_final_vs_ada4b7f1.json)):** only
  `thiamine_proton_lock_rationale` on the three reactions.
- **The gate on the saved final candidate ([after.json](after.json)):**
  - R408, R750 and R752 are in `skipped_curated_lock_reaction_ids`.
  - The changes go from 36 to 33, and the other 33 are identical. 317 rejected on both.
- **Provenance:** [the build.json](../../candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_v2_20261022.build.json)
  records:
  - HEAD `abe7e40`;
  - `source_unchanged_during_build: true` and `input_unchanged: true`;
  - microspecies.py `fc4ddff0…`, the step `cfe8543c…` and the data `2caf1981…`.
- **Build log:** the pre-existing libSBML error on R304's notes appears once in each candidate-line
  build (control and seeds). It does not appear in the default builds.

## 5. Solves

- No model LP was run for this change.
- Only notes changed, so the optimization problem equals `ef2265c1`'s (**provable** from the diff).
  The earlier WT and knockout results carry over.
- The full pytest suite does include existing tests that solve LPs.

## 6. Tests

- **Focused tests (verified):** `tests/test_thiamine_proton_lock.py` passes: 8 tests, 21 subtests.
  - **The gate:**
    - without a lock it strips all three H+;
    - a plain lock is honoured after an SBML roundtrip, with the stoichiometry unchanged;
    - the lock value is case- and space-insensitive but strict;
    - the dict lock is still honoured in memory.
  - **The step:**
    - 9 stale-input cases fail closed, including the H+ removed, the residual changed, a conflicting
      `false`, a dict or bool lock, and the lock without evidence notes;
    - scope, value and note-text rejections (including `pH < 7.3 & neutral`);
    - exact notes, idempotence, and an SBML roundtrip after which the gate skips all three;
    - the flag dependency.
- **Full suite, per test (verified, [logs/pytest_summaries.txt](logs/pytest_summaries.txt)):**

  | Run | Passed | Failed | Errors | Skipped | Subtests |
  | --- | --- | --- | --- | --- | --- |
  | Before | 232 | 11 | 17 | 31 | 591 |
  | Final code | 240 | 11 | 17 | 31 | 612 |

  - The only per-test difference is the 8 new passes.
  - An intermediate run on the pre-audit code also differed only by these tests.
  - The existing gate callers, `test_dipeptide_chemistry` and `test_quinone_pipeline`, pass as before.

## 7. Independent audit and responses

The audit ([AUDIT.md](AUDIT.md)) was a fresh read-only subagent. It audited `ada4b7f1`. All 9 claims
were supported; C7 was supported with gaps.

**Fixed:**
- **(a)1, the wrong rationale in the model notes.** The note said the H+ is off "only because ATP, ADP,
  AMP and phosphate are still neutral". It now says the residual comes from the mixed forms (the
  KEGG-form thiamine phosphates and HET-P with the neutral ATP family), and that only moving all of them
  to pH 7.3 together removes it. The same correction is in TASK.md. The candidate was rebuilt (`v2`).
- **(a)2, "including when a build starts from a saved model".** The note now says the lock works with
  the gate from the 2026-10-22 `microspecies.py` on, and on a model read back from SBML. The HEAD/main
  gate ignores the lock.
- **(a)3, (a)4 and (a)5.** TASK.md now mentions the dipeptide lock and the shifted line numbers. This
  report no longer claims "no LP" for the test suite or "every build" for the R304 error.
- **(c)1, note texts were not checked at runtime.** The step now rejects any note text with
  `< > & ' "`. Tested with the audit's own example, `pH < 7.3 & neutral`, and a missing text.
- **(c)2, non-string note values crashed.** They are now compared as they are, and a dict or bool lock
  fails with `ValueError` (tested).

**Not changed, reported:**
- **(c)3, the gate's bare key.** Any note with this key locks a reaction, and the skipped list mixes
  three reasons. That was authorized as planned; no functional defect.
- **(d), relying on the lock:** see §8.

## 8. Limits and next gates

- **Gate version:** the lock works only with this gate (`fc4ddff0` or later). On main or older
  worktrees the gate strips the three H+ silently (verified by the audit).
- **Other stoichiometry editors ignore the lock:** `tools/update_model.py:303` (an H/O balancer that
  ignores charge), the metadata selection, hydroxide normalisation, the patches, the energy candidates
  and the duplicate merge. R752, R408 and R750 are not in the selection spec (audit).
- **Once saved, the lock is unconditional:**
  - After a partial migration, say ATP, ADP, AMP and phosphate without the thiamine species, these
    reactions would sit at +3 / +2 / −3 H+. They would be listed only as skipped, not rejected.
  - The step's {H: 1, charge: 1} pin applies only when it runs.
  - A check that a locked reaction's residual is either {} or its reviewed value would close this gap.
    Not done (out of scope).
- **The saved candidate is not stable under the gate:** a gate run on it would still change 33 other
  reactions (inferred effect on a rebuild; not run).
- **R2041's lock** is still lost in saved models.
- **Silent loss of notes on SBML write** (libSBML error code −3) is a general, pre-existing risk. It
  affects R304's notes in every candidate-line build.
- **The three reactions stay one H+ off** in this model's mixed convention, by design, so balance
  audits will still flag them.
