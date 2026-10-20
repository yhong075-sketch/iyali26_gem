# R2190-RIBULOSE: D-ribulose written as R2190's product

Scope, authorization and acceptance criteria are in [TASK.md](TASK.md). The evidence for the change
is in the read-only review [r2190_identity_review_20261019](../r2190_identity_review_20261019/REPORT.md).
Nothing was committed or pushed.

Evidence tiers:
- **verified:** checked in this run (2026-10-19) and the output is linked here;
- **provable:** follows mathematically from verified facts;
- **inferred:** reasoned, with no single source.

## 1. Summary

- **New candidate (verified):** `model/candidates/E5_coq9_resp_pkg_r132_r2190_ribulose_20261019.xml`,
  SHA-256 `36f041930e5fee211395a45da4a36bb991fe9285637c409b07cc6d04114812ac`.
  - Baseline: a56f3649 (`E5_coq9_resp_pkg_r132_20261012.xml`).
  - R2190 now reads NAD + D-arabinitol ⇌ NADH + **D-ribulose**. D-ribulose is a new species,
    `m2061[C_cy]`.
  - R2190 keeps its bounds, GPR (`YALI1F03452g`), EC and annotations. Only five curation notes
    are added.
- **Only the targets changed (verified).** Whole-model diff against a56f3649:
  - 1 changed reaction (R2190: product and notes);
  - 1 added species (`m2061[C_cy]`);
  - 0 gene, objective, group or compartment changes.

  The raw XML diff is 35 lines: the new species block, R2190's notes, and one `speciesReference`.
- **The control build reproduces a56f3649 byte for byte (verified).** It is the same command
  without the flag. A rebuild with the final, post-audit code reproduces the candidate byte for byte
  (verified).
- **R2190 balances (verified):** the element and charge residual is {}.
- **LPs (8, all optimal; verified):**
  - WT growth: a56f3649 1.430805735352851 h⁻¹; new candidate 1.4308057353528474 h⁻¹; relative
    difference 2.5e-15.
  - With every boundary closed: R1964, R2041 and R2042 are all [0, 0].
  - For comparison, the orphan line recorded one-sided closed ranges on a56f3649: R1964 [−1000, 0],
    R2041 [−1000, 0], R2042 [0, 1000]. That was a GLPK FVA in its own setup
    (`orphan_reactions_20261019/loop_check.txt`, other worktree), so the before/after comparison spans
    two protocols. a56f3649's closed ranges were not re-run here.
- **pytest (verified):** identical to a clean `git archive 285eaee` copy, test by test, apart from
  the 6 new passing tests. The 11 failures and 14 errors are pre-existing at clean HEAD.
- **Independent audit (§6):** it found no undeclared change. It raised six low or informational
  items; the guard, test and wording ones are fixed.

## 2. What changed and why

| Entity | a56f3649 | New candidate |
| --- | --- | --- |
| R2190 | `m122 + m1823 ⇌ m123 + m454` (D-arabinose) | `m122 + m1823 ⇌ m123 + m2061` (D-ribulose) |
| R2190 notes | none | `r2190_product_curation`, `_previous_product`, `_evidence`, `_limits`, `_authorization` |
| `m2061[C_cy]` | absent | `D-ribulose_C5H10O5`, C5H10O5, charge 0, `C_cy`, SBO:0000247, KEGG C00309, CHEBI:17173, MNXM1371690, one note |
| `m454[C_cy]` (D-arabinose) | R1063, R2041, R2190, R254, R255 | R1063, R2041, R254, R255 (species itself unchanged) |

**Why (from the review, where each line is tiered).**
- R2190's name, EC 1.1.1.250 and every attached cross-reference give D-ribulose: KEGG R07134,
  MNXR197981 and SEED rxn04932/rxn15657. Rhea 17389, the EC's Rhea reaction, does too. All were
  opened.
- `m454` is a genuine D-arabinose node.
- The gene's purified enzyme (YlAraDH1, Wang *et al.* 2020) reduces only ketoses.

The verdict "product wrong" is inferred (strong). Whether the true product is D-ribulose (C-2) or
D-xylulose (C-4) remains open.

## 3. Implementation (the `curation-change` steps)

| Step | File | Notes |
| --- | --- | --- |
| Data | [`model/curation/chemistry/r2190_product.json`](../../curation/chemistry/r2190_product.json) | Curation ID `R2190-RIBULOSE-20261019`; exact before/after stoichiometry, species, the replaced species' usage before and after, new species spec, evidence (URL, tier, date), reviewer, limits, notes |
| Code | [`platform/scripts/gem_annotate/r2190_product.py`](../../../platform/scripts/gem_annotate/r2190_product.py) | `apply_r2190_product(model, spec=None)`: see the guard list below |
| Wiring | `main.py` (`build_reference_chain(..., r2190_ribulose_product=False)`, called right after the R132 step); `cli.py` flag `--r2190-ribulose-product` | Opt-in; off by default |
| Test | [`platform/tests/test_r2190_product.py`](../../../platform/tests/test_r2190_product.py) | 6 tests, 10 subtests, in-memory fixture, no optimization |

**What `apply_r2190_product` checks before any edit:**
- the scope fields, including that the curated before- and after-states differ only by
  `m454` → `m2061`, with every shared coefficient equal;
- that the curated after-state balances;
- R2190's name, GPR, EC and bounds;
- the formula, charge and compartment of the four species, and `m454`'s name;
- that no curation note conflicts;
- the exact before-state stoichiometry;
- that the ID `m2061` is unused in every compartment;
- `m454`'s exact usage.

It then adds the species, swaps the product, adds the notes and asserts the after-state. A rerun
returns `already_correct` after re-checking:
- the species fields;
- that `m2061` is unique across compartments;
- its usage, `m454`'s usage and the notes.

**The focused test covers:**
- product-only change and idempotence;
- an SBML round trip with the deterministic writer;
- ten stale-input cases that fail closed with no edit, each checked against the guard message it
  must raise: name, GPR, EC, bounds, an added H⁺, species formula, `m454` name, an extra `m454`
  user, an `m2061` ID in another compartment, and a conflicting note;
- four rejected curations: out of scope, an extra product, a changed shared coefficient, and an
  unbalanced new species;
- a rerun on a drifted after-state is rejected;
- the flag is opt-in. That flag-off changes nothing is shown by the byte-identical control build,
  not by a unit test.

**Fixes made during this run.**
- The first version of one note had an apostrophe, which the SBML writer escapes as `&apos;`. The
  round-trip test caught this before any build, and the note was reworded. That focused-test run
  was not logged.
- After the audit (§6), the guards and the test were tightened. The model output is unchanged: a
  rebuild with the final code reproduces the candidate byte for byte.

**Builds (verified; logs in [logs/](logs/)).** Run from `platform/` with the main checkout's `.venv`
(Python 3.13.5, cobra 0.30.0), using a56f3649's recorded argv plus or minus the flag:

```bash
python -m scripts.gem_annotate --research-root "/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research" --offline --no-solve --coq9-curation qcycle --energy-candidate E5 --coq9-biomass-dilution --coq9-respiratory-package --r132-mito-adh --r2190-ribulose-product --output-model ../model/candidates/E5_coq9_resp_pkg_r132_r2190_ribulose_20261019.xml
```

- The `.build.json` records: head `285eaee`, the dirty status (the uncommitted files of this change),
  `requested_build_complete: true`, `source_unchanged_during_build: true`, and the step record
  `status: applied`.
- The control used the same command without the flag, writing to the git-ignored
  `platform/.scratch/r2190_control_20261019/`. Its output SHA is `a56f3649…`, identical by `cmp`.
- Both logs carry the same libSBML "Setting notes on … R_R304" error. It is pre-existing: the
  byte-identical control shows it too.
- **Provenance note.** The candidate's own `.build.json` records the pre-audit `r2190_product.py`
  (SHA-256 `cbfc1a05…`).
  - The final file is `fbdb2633…`. With it, the same command writing to
    `platform/.scratch/r2190_verify_20261019/` reproduced SHA-256 `36f04193…` (verified by `cmp`).
    Its record is [logs/build_verify_final_code.build.json](logs/build_verify_final_code.build.json).
  - `main.py`, `cli.py` and the curation JSON are unchanged between the two builds.
  - The candidate file was not rebuilt in place, because the exporter refuses to overwrite.

## 4. Verification

**Whole-model diff (verified).**
- Tool: [`whole_model_diff.json`](whole_model_diff.json), from
  `model/reports/coq9_respiratory_package_20261008/model_diff.py`. It compares every reaction
  (stoichiometry, bounds, GPR, annotation, notes, subsystem), metabolite (name, formula, charge,
  compartment, annotation, notes) and gene, plus the objective, groups and compartments.
- Raw XML: [`whole_model_xml.diff`](whole_model_xml.diff).
- Result: only the targets in §2 differ. Counts go from 2306/1875/1073 to 2306/1876/1073
  (reactions/species/genes).

**Balance (verified, [static_check.json](static_check.json)).**
- R2190's residual is {}.
- `m2061` is used by R2190 only.
- `m454`'s other reactions R1063, R2041 and R254 balance.
- R255 (D-arabinose 1-dehydrogenase, NADP) has a charge residual of +1. That is pre-existing: R255
  is unchanged in the diff. It is noted for the chemistry line, not touched.

**LPs (verified, [lp_results.json](lp_results.json), runner [scripts/run_lps.py](scripts/run_lps.py)).**
- Context: the `run_controls.py` context of `coq9_respiratory_package_20261008`.
  - Profile: PO1f SD-Leu, `po1f_sd_leu_accrispr_v1`.
  - Solver: Gurobi with 1 thread, Presolve 0, tolerances 1e-7, 60 s limit.
- Exactly 8 optimization calls (counter-asserted):

| # | LP | Model | Status | Value |
| --- | --- | --- | --- | --- |
| 1 | WT growth | a56f3649 | optimal | 1.430805735352851 |
| 2 | WT growth | new | optimal | 1.4308057353528474 |
| 3–4 | min/max R1964, all 183 boundary reactions closed, xMAINTENANCE relaxed to include 0 | new | optimal ×2 | 0.0 / 0.0 |
| 5–6 | min/max R2041, same closure | new | optimal ×2 | 0.0 / 0.0 |
| 7–8 | min/max R2042, same closure | new | optimal ×2 | 0.0 / 0.0 |

- In both WT optima, R2190, R1964, R2041, R2042 and R803 carry 0; R714 carries −0.338.
- Stop conditions: none triggered.

**pytest (verified, [pytest_comparison.json](pytest_comparison.json)).**
- Baseline: a clean `git archive 285eaee` copy in the session scratchpad.
- After: this worktree with the final code.
- Command, from `platform/`:

  ```bash
  IYALI26_RESEARCH_ROOT="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research" python -m pytest -q -rA -p no:cacheprovider --junitxml=<log>.xml
  ```

- The JUnit files and text logs are in [logs/](logs/):
  - `pytest_baseline_cleanhead.*`;
  - `pytest_after_final.*`;
  - `pytest_after_preaudit.*` (the pre-audit code, 170 passed, same result).
- Reproduce with `python scripts/compare_pytest.py logs/pytest_baseline_cleanhead.xml logs/pytest_after_final.xml out.json`.

| Run | Passed | Skipped | Failed | Errors |
| --- | --- | --- | --- | --- |
| Before | 165 | 31 | 11 | 14 |
| After (final code) | 171 | 31 | 11 | 14 |

- No test changed outcome. The only additions are the 6 new tests.
- The failures and errors are the same tests in both runs:
  - `test_lp_sn12_candidate` (14 errors);
  - `test_coa_protonation_curation` (6 failures);
  - `test_vacuole_*` (3 failures);
  - `test_dfba_new_fn_essentiality` (1 failure);
  - `test_lipid_moiety_ledger` (1 failure).

  Their cause was not investigated.
- Some repository tests run solves of their own. Those were not counted in the 8-LP budget. The
  suite ran four times in all: the stopped run, the clean baseline, the pre-audit run and the final
  run.
- **The stopped run.** The first baseline attempt ran in this worktree while `main.py` and `cli.py`
  were being edited, so it could not count as a clean baseline.
  - It ran the whole suite and was stopped during its summary.
  - Its log is `logs/pytest_baseline_contaminated_stopped.txt`. It is not used.
- **Run isolation.** The clean baseline, the pre-audit run and the 8 LPs overlapped in time.
  - The outcomes are unaffected: every solve was optimal, and the test outcomes are identical.
  - The final verification build and the final pytest ran alone, one after the other.

**Mechanism (provable, now also observed).** `m2061` has one reaction, so steady state fixes
R2190 = 0. With all boundaries closed, D-arabinitol's only other outlet R1980 leads to a closed
exchange, so R1964 = 0. D-arabinose's other reactions lead to a closed exchange (R1063) or to the
dead-end lactone branch, so R2041 = 0 and with it R2042. The six LPs show exactly this.

## 5. Limits and what was not tested

- **Regiochemistry open:** D-ribulose follows the databases and the EC. If product identification
  ever shows D-xylulose, the review's Option B (merge with R1964) becomes right. The parallel
  session "Evidence review: YlAraDH1 as GPR for orphan R1964" is looking at this.
- **R2190 is now blocked** by construction. Its gene's KO cannot change any call through R2190.
  This is the honest state while no D-ribulose producer, consumer or transporter is evidenced. It
  is not evidence of the gene's function.
- **Not tested:**
  - the single-KO screen (needs separate authorization);
  - E0, reference or other chains with the flag (the guard would fail closed if R2190 differed);
  - loopless or sampling analyses;
  - the R2041 forward-only plan, which is independent: either change alone breaks this loop.
- **Not done:** commit and push (excluded by the user).
- **Math:** nothing to improve here; the change is a single stoichiometric substitution.
- **A stronger test:** a stoichiometric dead-end census before and after. This change adds exactly
  one blocked reaction and one dead-end species (inferred from §4; a census was not run).

## 6. Independent audit

**Setup.** A fresh read-only reviewer (subagent, no conversation history) received only TASK.md,
this report, the data, code, test, candidate, diffs, records and logs. It was not allowed to edit,
build, run pytest or LPs, or use the network. It wrote one temporary diff file to the scratchpad.

**Result.** No undeclared change. It checked 26 claims:

| Verdict | Claims | Which |
| --- | --- | --- |
| Supported | 21 | including candidate SHA, control identity, whole-model diff by its own cobra comparison, balance, guard order, wiring, build record, 8 counted solves with no hidden solve, closure, mechanism by connectivity, and pytest (it re-parsed the JUnit files and hashed all 1,548 blobs of the clean copy against `285eaee`) |
| Partially supported | 2 | D1, D2 below |
| Unverified | 2 | the unlogged apostrophe fix; `IYALI26_RESEARCH_ROOT` not visible in the logs |
| Outside audit scope | 1 | the review's evidence; spot-checked from cached copies |

**Defects and what was done:**

| # | Severity | Finding | Action |
| --- | --- | --- | --- |
| D1 | low | "±1000" overstated the orphan line's one-sided, other-protocol ranges | Reworded in §1 |
| D2 | low | "Stopped part-way" was inaccurate; JUnit files only in the scratchpad; pytest command elided | Reworded in §4; JUnit and text logs copied to `logs/`; full command given |
| D3 | low | Scope guard did not require after − before = {m2061} or equal shared coefficients (then only the post-edit assert caught it); `m454` name not checked; rerun branch did not re-check `m2061` uniqueness. Unreachable with the shipped spec | Guards added (§3). Final code rebuild is byte-identical |
| D4 | low | Stale subtests accepted any `ValueError`; flag test only covers argparse | Each case now asserts its guard message; 3 scope cases and 1 rerun case added. Flag-off rests on the control build (stated in §3) |
| D5 | info | No chain prerequisite for the flag (unlike R132) | Kept: the guard checks the exact before-state; other chains are untested (§5) |
| D6 | info | Clean baseline pytest, pre-audit pytest and LPs overlapped in time | Disclosed in §4; final build and pytest ran alone |

The fixes in D3 and D4 were not re-audited. They are covered by the focused test, the final full
pytest and the byte-identical final-code rebuild.

## 7. Files

| Path | What |
| --- | --- |
| `TASK.md`, `before.json`, `scripts/record_before.py` | Scope and baseline state |
| `whole_model_diff.json`, `whole_model_xml.diff` | Whole-model verification |
| `static_check.json`, `scripts/static_check.py` | Balance, usage and build record |
| `lp_results.json`, `scripts/run_lps.py`, `logs/run_lps.log` | The 8 LPs |
| `pytest_comparison.json`, `scripts/compare_pytest.py`, `logs/pytest_*.{txt,xml}` | pytest before and after |
| `logs/build_*.log`, `logs/build_control_no_flag.build.json`, `logs/build_verify_final_code.build.json` | Build logs, the control's record, the final-code verification record |
| `../../candidates/E5_coq9_resp_pkg_r132_r2190_ribulose_20261019.{xml,build.json,coq9_genes.tsv}` | The candidate |
