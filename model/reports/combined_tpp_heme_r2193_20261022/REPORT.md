# Combined candidate: TPP line + heme/iron + R2193 line (2026-10-22)

Scope, authorization, budget and acceptance criteria are in [TASK.md](TASK.md). Exact commands are in
[commands.txt](commands.txt). The two lines are documented in
[tpp_requirement_20261022](../tpp_requirement_20261022/),
[thiamine_phosphate_chemistry_20261022](../thiamine_phosphate_chemistry_20261022/),
[thiamine_proton_lock_20261022](../thiamine_proton_lock_20261022/),
[heme_iron_biomass_20261022](../heme_iron_biomass_20261022/REPORT.md) and
[r2193_iron_gpr_change_20261022](../r2193_iron_gpr_change_20261022/REPORT.md).

Evidence tags:
- **verified** = checked this session with the output kept here;
- **inferred** = reasoned;
- **recalled** = not opened.

## 1. What changed and why

- **Two opt-in lines branched from `cf0825b8` on 2026-10-22:**
  - the TPP line (PR #16): TPP in biomass, R1342 closed, thiamine-phosphate chemistry, an H⁺ lock;
    latest `2fe51cc2`;
  - the heme/iron + R2193 line (PR #17): heme a and iron in biomass, R2193 on the FET3/FTR1 pairs;
    latest `e1586027`.
- **Why they would not combine.** Both edit `biomass_C`. The heme/iron guard required the base's
  exact `biomass_C`, so it would have stopped any build that also carried the TPP term.
- **Guard change** (`heme_iron_biomass.json`, `heme_iron_biomass.py`). The step accepts exactly two
  reviewed starting states for `biomass_C` (`accepted_before_states`):
  - `base`, the composition reviewed on `cf0825b8`;
  - `after_TPP-REQUIREMENT-20261022`, the base plus TPP `m975[C_cy]` at −1.2e-6 and the 7
    `tpp_requirement_*` note keys. That is the TPP curation's recorded after-state, equal to
    `2fe51cc2`'s `biomass_C` (verified, [before.json](before.json)).

  Any other composition or set of note keys stops the build. The step reports the state it matched.
- **Builder rule (after the audit).** `--heme-iron-biomass` together with `--tpp-requirement` now
  requires `--thiamine-phosphate-chemistry` and `--thiamine-proton-lock`. Only that full combination
  was reviewed.
- **Step order.** The heme step must run after the TPP steps, as `main.py` fixes it: the TPP guard
  rejects a `biomass_C` that already carries heme terms. The model content of existing candidates is
  unchanged; both controls reproduce.
- **Combined candidate.**
  [E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_20261022.xml](../../candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_20261022.xml),
  SHA-256 `2c9cebd0a5c53ec18bf1a583327099dc9e53fcf6b2fcbedd66cc7b9c54ea5a7b`. It is `cf0825b8` plus
  `--tpp-requirement --thiamine-phosphate-chemistry --thiamine-proton-lock --heme-iron-biomass --r2193-iron-gpr`.

## 2. Build and whole-model verification (verified)

- **Builds:**
  - TPP-only argv: `2fe51cc2`, identical to the PR #16 candidate;
  - heme-only argv: `e1586027`, identical to the PR #17 candidate;
  - combined: `2c9cebd0` under `PYTHONHASHSEED` 0 and 1 (`cmp` identical). The logs show the heme
    step matched `after_TPP-REQUIREMENT-20261022`; the heme-only control matched `base`.
- **After the audit fix:**
  - a rebuild still gives `2c9cebd0`;
  - `--tpp-requirement --heme-iron-biomass` without the chemistry and lock steps is refused before
    any build ([log](logs/build_unreviewed_combo_rejected.log)).
- **Diffs** (`combined_candidate_20261020/scripts/model_diff.py`):
  - [against `2fe51cc2`](model_diff_vs_2fe51cc2.json): 9 differences, exactly the heme/iron + R2193
    changes:
    - `biomass_C` + `m1635[C_cy]` −1e-6 and `m948[C_cy]` −3.04e-5, plus 7 heme notes;
    - R2193's GPR and notes;
    - 4 added genes;
    - YALI1F07747g's link.
  - [against `e1586027`](model_diff_vs_e1586027.json): 35 differences. Fields and values are identical
    to the TPP line's own diff from `cf0825b8`
    ([model_diff_tpp_parent_vs_base.json](model_diff_tpp_parent_vs_base.json)):
    - `biomass_C` + `m975[C_cy]` −1.2e-6 and 7 TPP notes;
    - R1342 bounds, annotation and notes;
    - the R752, R408 and R750 notes;
    - 7 thiamine-phosphate species.
  - The audit confirmed both with raw XML diffs: base→`2fe51cc2` equals `e1586027`→combined
    line for line, and base→`e1586027` equals `2fe51cc2`→combined.
- **Provenance.** `PYTHONHASHSEED` is recorded in commands.txt, not in build.json.

## 3. Tests (verified)

- **Focused tests**, both lines
  ([pytest_focused_after_audit.log](logs/pytest_focused_after_audit.log)): 30 passed, 132 subtests.
  The new cases:
  - the heme step applies on top of an in-memory TPP model and reports the matched state;
  - a different TPP coefficient, or TPP without its notes, is rejected;
  - three scope cases on `accepted_before_states` are rejected;
  - TPP without chemistry or lock, plus heme, is refused.
- **Full suite, before the audit fix** ([pytest_after.log](logs/pytest_after.log), run alone): 251
  passed, 11 failed, 18 errors, 31 skipped. Against the last pre-merge run there are only additions:
  - `main`'s 19 TPP/thiamine tests and 1 new heme test;
  - 1 new collection error in `tests/test_gem_review_lit_caching.py` (from PR #14): `httpx2` is not
    installed in this environment, unrelated to this change.
- **Full suite, after the audit fix** ([pytest_after_audit.log](logs/pytest_after_audit.log), run
  alone): identical test by test to the pre-fix run (251 passed; 668 subtests, including the 2 new
  rejection cases).

## 4. Solves and knockouts (verified)

Task A medium + PO1f; Gurobi with FeasibilityTol and OptimalityTol 1e-9.

### WT and key fluxes ([solves.json](solves.json), 12 LPs)

- Each model gets a reset, checked WT solve, then a pFBA. The pFBA is not reset or checked.
- The route fluxes below come from the pFBA. The audit confirmed with minimum-flux LPs that they are
  forced at optimal growth.

| Model | WT growth | Δ vs base | Iron uptake R1189 | Heme route (R302 = R383 = R1709) | TPP route R752 |
| --- | --- | --- | --- | --- | --- |
| `cf0825b8` base | 0.3124267446 | — | 0 | 0 | 0 |
| `2fe51cc2` TPP | 0.3124267405 | −4.1e-9 | 0 | 0 | 3.75e-7 |
| `e1586027` heme/iron | 0.3124262956 | −4.49e-7 | −9.81e-6 | 3.12e-7 | 0 |
| `2c9cebd0` combined | 0.3124262915 | −4.53e-7 | −9.81e-6 | 3.12e-7 | 3.75e-7 |

The two costs add (residual 9e-15, audit). Every checked WT point passes.

### Knockouts

All three models, 1,074, 1,078 and 1,078 genes, 3,233 LPs. Outputs: [ko_compare.tsv](ko_compare.tsv)
and [ko_compare_summary.json](ko_compare_summary.json).
- **Method:** a model reset before each knockout (the audit confirmed Gurobi then starts from
  scratch), and the targeted check (knocked-out reactions |v| ≤ 1e-9, other violations ≤ 1e-6).
- **No point failed.** Every closed reaction is at 0, the largest bound violation is 0, and the largest
  |S v| is 3e-12. No fresh re-solve was needed. No point shows the GLPK heme-skip pattern (about 3e-7).
- **The same 20 knockouts are infeasible in all three models.** They are unresolved without a
  re-solve, as defined in TASK.md before the sweep.

**Class changes** (identical at 1, 5, 10 and 15 %). Identities are KEGG names for the model-annotated
IDs ([heme report](../heme_iron_biomass_20261022/gene_identities.tsv),
[R2193 investigation](../r2193_iron_import_gpr_20261022/candidate_genes.tsv),
[TPP report](../tpp_requirement_20261022/)). They are database or orthology identities, not
experimental validation in *Y. lipolytica*; none has a verified gene name.

| Gene | Function (KEGG / orthology) | TPP parent → combined | Heme parent → combined | List |
| --- | --- | --- | --- | --- |
| YALI1C07816g | 5-aminolevulinate synthase | non-essential → essential | — | positive |
| YALI1F27686g | porphobilinogen synthase | non-essential → essential | — | positive |
| YALI1F34051g | hydroxymethylbilane synthase | non-essential → essential | — | positive |
| YALI1F13039g | uroporphyrinogen-III synthase | non-essential → essential | — | positive |
| YALI1C02391g | uroporphyrinogen decarboxylase | non-essential → essential | — | positive |
| YALI1E35410g | coproporphyrinogen III oxidase | non-essential → essential | — | negative |
| YALI1D30196g | protoporphyrinogen oxidase | non-essential → essential | — | negative |
| YALI1F25924g | ferrochelatase | non-essential → essential | — | positive |
| YALI1F31135g | heme o synthase | non-essential → essential | — | positive |
| YALI1F32002g | heme a synthase (COX15) | non-essential → essential | — | positive |
| YALI1B19490g | adrenodoxin-NADP⁺ reductase | non-essential → essential | — | positive |
| YALI1B03314g | mitochondrial ferredoxin-2 | non-essential → essential | — | negative |
| YALI1D09277g, YALI1D09308g | FET3/FTR1 pair D09 (reciprocal best hits) | absent → non-essential | — | positive |
| YALI1D08684g, YALI1D08564g | FET3/FTR1 pair D08 (paralogs, one-way) | absent → non-essential | — | negative |
| YALI1E25443g | THI80-like (R752), no verified name; identity per the TPP report | — | non-essential → essential | positive |

- **Why 12 heme genes, not the 13 expected.** The 13th, YALI1F07747g, was essential only through R2193
  in `7f836c3a`. After the R2193 swap it is non-essential in both the TPP parent and the combined
  model.
- So the combined model's class changes are exactly the union of the two lines' changes against
  `cf0825b8`.

**Scores.**
- Negatives are the user-defined class; the 16 mtDNA genes are out of scope.
- List: `consensus_essential_genes.csv`, `1e887f5a…`, 1,612 positives.
  - In-model positives: 322 / 324 / 324.
  - Out-of-scope positives (no model gene): 1,290 / 1,288 / 1,288.
  - 11 in-model positives are among the 16 unresolved (infeasible) knockouts in each model.
- Development reference, not independent validation.

| Cutoff | Model | TP / FN / FP / TN / unresolved | Recall | Specificity | Precision | MCC |
| --- | --- | --- | --- | --- | --- | --- |
| 15 % | `2fe51cc2` TPP | 84 / 227 / 33 / 698 / 16 | 0.270 | 0.955 | 0.718 | 0.326 |
| 15 % | `e1586027` heme/iron | 92 / 221 / 36 / 697 / 16 | 0.294 | 0.951 | 0.719 | 0.342 |
| 15 % | `2c9cebd0` combined | 93 / 220 / 36 / 697 / 16 | 0.297 | 0.951 | 0.721 | 0.345 |

**Other cutoffs** (TP / FN / FP / TN / unresolved, MCC; recall, specificity and precision for every
cutoff are in `ko_compare_summary.json`):

| Model | 1 % | 5 % | 10 % |
| --- | --- | --- | --- |
| TPP | 66/245/26/705/16, 0.285 | 72/239/28/703/16, 0.300 | 76/235/33/698/16, 0.298 |
| heme/iron | 74/239/29/704/16, 0.303 | 80/233/31/702/16, 0.317 | 84/229/36/697/16, 0.315 |
| combined | 75/238/29/704/16, 0.306 | 81/232/31/702/16, 0.321 | 85/228/36/697/16, 0.318 |

**Cross-solver check.** The heme parent's per-gene classes and scores on Gurobi at 1e-9 equal the
earlier GLPK-checked sweep. That sweep ran on v2 (`2cf19c69`), which is the same optimization problem
as `e1586027` (they differ only in R2193 note text, per `model_diff_v3_vs_v2.json`; the v2 file has
since been deleted).

## 5. What the results mean

- **Shown, for WT and single-gene knockouts on the Task A medium at 1e-6 growth resolution:**
  - the two lines combine without interaction;
  - growth costs add;
  - class changes are exactly the union, and all other knockouts are unchanged.

  Smaller interactions, double knockouts and other media are not tested.
- **Numerics.** Gurobi at 1e-9 with a reset per knockout solves all three models with residuals of
  about 1e-12. It shows neither GLPK artifact (stale points, dropped heme demand). Use it for future
  screens of these candidates.
- **Biology.** Nothing beyond the two lines' own claims.

## 6. Independent audit and responses

The full audit is in [AUDIT.md](AUDIT.md).
- **Supported:** C1 and C3–C7, and C8's numbers.
- **Partially supported:** C2, for `biomass_C` only.

| Issue | Response |
| --- | --- |
| 1. Unreviewed combination possible (TPP without chemistry or lock, plus heme) | Fixed. The builder refuses `--heme-iron-biomass` with `--tpp-requirement` unless both `--thiamine-phosphate-chemistry` and `--thiamine-proton-lock` are set. Tested (2 subcases) and confirmed on the real CLI; the rebuild still gives `2c9cebd0`. |
| 2. TPP note values not pinned | Accepted as a limit: only keys are pinned, because note text is escaped differently in a build and in a saved file. The TPP step checks its own values before the heme step runs. |
| 3. State 2 pinned in JSON; order dependence | Accepted. The data file is the curated record, as in every step. Order dependence is now stated in §1 and in a `main.py` comment. |
| 4. Process (tests pending, STATE.md, README) | §3 and §6 updated, and the STATE.md entry is written. The README row that appeared during the audit was added by this session, not by another one. |
| 5. Reporting gaps | Recall, specificity, precision, out-of-scope positives, the 11 unresolved positives and gene identities have been added (§4). |
| 6. C8 wording | Corrected in §4: the GLPK sweep was on v2 `2cf19c69`, the same LP problem. |
| 7. pFBA not reset or checked | Wording corrected in §4; the route values are confirmed by the audit's minimum-flux LPs. |
| 8. Provenance | The reference diff is kept ([model_diff_tpp_parent_vs_base.json](model_diff_tpp_parent_vs_base.json)) and the focused-test log is kept. `PYTHONHASHSEED` stays in commands.txt (inherited gap). |
| 9. Scope wording | Qualified in §5. |
| 10. 13 vs 12 genes | Explained in §4. |

**Full suite after the audit fix** ([pytest_after_audit.log](logs/pytest_after_audit.log), run
alone): 251 passed, 11 failed, 18 errors, 31 skipped, 668 subtests. Identical test by test to the
pre-fix run; the 2 extra subtests are the new rejection cases.

**Limits.**
- All limits of both lines still apply: see the open-problems entry at the top of
  [STATE.md](../../STATE.md).
- Not committed or pushed.
