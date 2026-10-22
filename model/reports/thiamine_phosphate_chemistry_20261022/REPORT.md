# Thiamine phosphate chemistry: implementation report

Scope, authorization, budget and acceptance criteria are in [TASK.md](TASK.md). Exact commands are in
[commands.txt](commands.txt). The independent audit is in [AUDIT.md](AUDIT.md); the responses are in
§7.

Evidence tiers used below:
- **verified**: run or checked in this session; the output is cited.
- **provable**: follows from verified facts.
- **inferred**: reasoned, with no direct check.
- **established/uncited**: textbook knowledge; no source was opened this session.

## 1. Summary

- **New candidate (verified):**
  - `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_v2_20261022.xml`, SHA-256
    `ef2265c1fd6c3cd12eebacc8c03901091f90289a6f77a7be66d2ede9d7379440`.
  - It is `8ae6c93b` (the TPP requirement candidate) plus the opt-in step
    `--thiamine-phosphate-chemistry`.
  - It is a candidate, not a reference model, and it is uncommitted.
  - The first build, `…_thiamine_p_20261022.xml` (`e8a54520`), is superseded. It has the same
    formulas, charges, IDs and LP; only two note texts per species differ, corrected after the audit.
    It is kept until the user decides.
- **What changed (verified, 29 differences against `8ae6c93b`):**
  - **Thiamine monophosphate (TMP), `m977`/`m973`:** C6H4Cl3NO, 0 → C12H18N4O4PS, +1. The
    trichloromethoxypyridine identity (MNXM1137754, envipath, InChI) is replaced by KEGG C01081,
    CHEBI:9533, MNXM662 and BiGG thmmp.
  - **Thiamine diphosphate (TPP), `m975`/`m971`/`m1116`:** no formula, 0 → C12H19N4O7P2S, +1. KEGG
    C00068, CHEBI:9532, MNXM256, thmpp.
  - **Thiamine triphosphate (TTP), `m976`:** no formula, 0 → C12H20N4O10P3S, +1. BiGG ttp /
    MNXM148527 are replaced by KEGG C03028, CHEBI:9534, MNXM1341, thmtp.
  - **HET-P, `m643`:** no formula → C6H10NO4PS, charge 0 (unchanged). KEGG C04327, CHEBI:17857,
    MNXM960, 4mpetz.
  - **R1342:** EC `2.7.4.16, 2.7.4.9` → `2.7.4.16`.
  - **Notes:** `thiamine_phosphate_chemistry_*` notes on all 8 targets.

  No stoichiometry, bound, GPR, objective, name or compartment changed.
- **Balance (verified):** of the 11 non-boundary reactions on these species, 6 now balance; only R933
  did before, trivially. No balanced reaction became unbalanced.
  - R752, R408 and R750 stay one H+ off. That H+ is correct for the pH 7.3 convention the pipeline is
    moving to (§3).
  - R752's charge residual went from 0 to +1 while its mass residual fell from six elements to one H.
- **LP (verified):** WT growth is identical, 0.31242674050092606 on both models.
- **Tests and audit:**
  - The full suite differs from the baseline only by the 5 new passing tests (§6).
  - The audit supported 9 of 11 claims, and partly supported C6 and C8.
  - The guard gaps and wording it found are fixed (§7).
- **One open decision:** with the KEGG forms, the pipeline's proton gate would strip the H+ from R752,
  R408 and R750 if it ever ran after this step. A lock note would prevent that; it is outside the
  approved scope (§8).

## 2. The change

| Item | Location |
| --- | --- |
| Curated data | [`model/curation/chemistry/thiamine_phosphate_chemistry.json`](../../curation/chemistry/thiamine_phosphate_chemistry.json) (`THIAMINE-P-CHEMISTRY-20261022`, SHA-256 `42af6975…`), generated from [before.json](before.json) by [scripts/make_curation_json.py](scripts/make_curation_json.py) |
| Guarded step | `platform/scripts/gem_annotate/thiamine_phosphate_chemistry.py` (`9a187849…`), `apply_thiamine_phosphate_chemistry` |
| Wiring | `cli.py` flag `--thiamine-phosphate-chemistry`; `main.py` calls the step right after `--tpp-requirement`, which it requires |
| Focused test | `platform/tests/test_thiamine_phosphate_chemistry.py` |
| Code diff | [logs/code_diff.patch](logs/code_diff.patch), which also holds the previous task's `--tpp-requirement` wiring |

**What the step checks before writing.** It fails closed and makes no partial edits.
- **The spec:**
  - schema version (an integer 1), ID, status, and the exact sets of species, reactions and context
    reactions;
  - each species' molecule, formula, integer charge and full annotation, against the reviewed
    chemistry pinned in the code (`REVIEWED`), so an edited data file cannot change them silently;
  - R1342's annotation changes only in `ec-code`.
- **The TPP requirement:** `biomass_C` must hold its TPP term.
- **The reaction context:** all 14 reactions on these species keep the stoichiometry the balance
  review used. For `biomass_C`, only the TPP term is checked.
- **The partners:** the formula and charge of the 24 partner species on those reactions (for example
  neutral ATP) must match the review. If a partner migrates to pH 7.3 first, the step stops instead of
  writing forms that would no longer balance.
- **Each species:** name, compartment and full reaction list. Formula, charge and annotation must be
  the recorded before-state or after-state. Annotation values are compared as sorted lists.
- **Notes:** compared after XML-unescaping.
- **Completeness:** either all 8 targets carry the change with this curation's notes, or none does.

A rerun returns `already_correct`.

## 3. Mechanism and evidence

**How the errors arose:**
- **Source models (verified):** `model/source/iyali26.xml` and `iyli21.xml` carry no formula or
  annotation for any of these 7 species.
- **The TMP mapping (verified in the build log):** `R1341 … mnxm=MNXM1137754`. MetaNetX 4.5 names
  that entry "TMP"; it is 3,5,6-trichloro-2-methoxypyridine (C6H4Cl3NO) from envipath.
- **Which mapping strategy produced it (inferred):** the exact-name lookup in `metabolites.py`
  (Strategy B), from reading the code. The audit did not reproduce this step by step.
- **TTP:** it matched BiGG `ttp` / MNXM148527, which has no formula.
- **TPP species and HET-P:** they stayed formula-less.
- **R1342's EC 2.7.4.9 (verified):**
  - It comes from MetaNetX 4.5's EC field for MNXR104907 ("2.7.4.16;2.7.4.9").
  - All of that reaction's cross-references are thiamine-phosphate kinase.
  - ExPASy (release of 28-Jan-2026): 2.7.4.9 = dTMP kinase; 2.7.4.16 = thiamine-phosphate kinase.

**Identity (inferred from reactions, with each reaction verified):** these are thiamine phosphates,
not thymidine nucleotides.
- The model's own KEGG reaction IDs, with equations read from KEGG:
  - R752 = R00619;
  - R753 = R00616;
  - R1342 = R00617;
  - R1337/R751 = R02135;
  - R408 = R04448.
- The source model keeps dTMP, dTDP and dTTP as separate species (m505, m506, m842), and names the
  exchanges R1341 and R1339 "thiamine(1+) mono/diphosphate exchange" (audit C11).
- TTP's identity rests on R753 alone.

**Formulas and charges (verified):**

| Molecule | KEGG formula (fully protonated) | ChEBI parent | MetaNetX 4.5 pH 7.3 form |
| --- | --- | --- | --- |
| TMP | C01081, C12H18N4O4PS | CHEBI:9533, +1 | MNXM662, C12H16N4O4PS, −1 |
| TPP | C00068, C12H19N4O7P2S | CHEBI:9532, +1 | MNXM256, C12H16N4O7P2S, −2 |
| ThTP | C03028, C12H20N4O10P3S | CHEBI:9534, +1 | MNXM1341, C12H16N4O10P3S, −3 |
| HET-P | C04327, C6H10NO4PS | CHEBI:17857, 0 | MNXM960, C6H8NO4PS, −2 |

- The BiGG IDs were checked only through MetaNetX cross-references.

**Why the KEGG forms (verified):**
- The model stores ATP, ADP, AMP, phosphate, diphosphate and HMP-PP neutral.
  - The pH 7.3 migration of the first five is `component_review` in `metabolite_microspecies.csv`.
  - HMP-PP has no row there.
- The KEGG forms match that state:

  | Reaction | Before | After |
  | --- | --- | --- |
  | R1337, R1342, R751, R753, R754 | off | balanced |
  | R933 | balanced (trivially) | balanced |
  | R752, R408, R750 | off | one H+ off ({H: 1, charge: 1}) |
  | R1343 | off | {charge: 2, H: 1} |
  | R1344 | off | {H: 1} |

  - With the pH 7.3 forms against today's partners, only 2 of the 11 balance.
- **R1343 and R1344** stay off in every consistent protonation state of `m889` and `m38` tried by the
  audit. This is a 2-H discrepancy in their lumped equations. `m889`, named "ribose-5-phosphate", is
  annotated KEGG C01068 / MNXM1092921, which is D-ribitol 5-phosphate.
- **The pH 7.3 check (verified twice):** with every partner also at its pH 7.3 form, all 9 thiamine
  reactions balance exactly as written. I checked this in planning, and the audit recomputed it.
  - So the H+ in R752, R408 and R750 belongs to the target convention. It was not edited.
  - R750 has no Rhea cross-reference (it is SEED rxn09632). Earlier drafts wrongly called its
    equation "Rhea"; it is simply pH 7.3-consistent.

## 4. Builds and whole-model verification

**Build history (logs in [logs/](logs/)):**
1. **First candidate `e8a54520`:**
   - The control rebuilt `8ae6c93b`.
   - Seeds 0 and 1 gave identical bytes.
   - The guard did not reject any build.
   - This candidate was audited.
2. **After the audit:**
   - The step was hardened and the note texts corrected (§7).
   - The builder refuses an existing output, so the final candidate has a `_v2` name.
   - Each final log's first line records its `PYTHONHASHSEED`.

**Final checks (verified, [logs/sha256_final.txt](logs/sha256_final.txt)):**
- **Control:** the 8ae6c93b argv without the new flag rebuilds `8ae6c93b…` byte for byte.
- **Seed independence:** with the flag, seeds 0 and 1 give the same `ef2265c1…`.
- **Where the builds are:** the control, seed-1 and pre-fix builds, including `e8a54520`, are under
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/thiamine_phosphate_chemistry_20261022/`.
- **Diff against `8ae6c93b` ([model_diff.json](model_diff.json)):** 29 differences, the same as for
  `e8a54520`; the audit checked them three ways (raw XML, element-level and cobra-level).
  - 6 species change formula, charge, annotation and notes; `m643` changes formula, annotation and
    notes, but its charge stays 0.
  - R1342 changes annotation (`ec-code` only) and notes.
  - Nothing else differs.
- **Diff against `e8a54520` ([model_diff_final_vs_e8a54520.json](model_diff_final_vs_e8a54520.json)):**
  7 differences, all metabolite notes (`…_evidence_tier` and `…_rationale`).
- **Balance:** see the table in §3 ([before.json](before.json), [after.json](after.json));
  `balanced_count` goes from 1 to 6.
- **Provenance:** [the build.json](../../candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_v2_20261022.build.json)
  records:
  - HEAD `abe7e40` and the dirty status;
  - `source_unchanged_during_build: true` and `input_unchanged: true`;
  - the step `9a187849…` and the data file `42af6975…`.

  `model/curation/README.md` was edited around the time of this build, so the hash recorded for it in
  build.json may be an earlier text. It is documentation and does not affect the model.
- **Build log:** the pre-existing libSBML error on R304's notes appears in every build, the control
  included.

## 5. Solves ([wt_spot_check.json](wt_spot_check.json))

| Model | Status | WT growth (Task A medium + PO1f, GLPK) |
| --- | --- | --- |
| `8ae6c93b` | optimal | 0.31242674050092606 |
| `e8a54520` | optimal | 0.31242674050092606 |

- Two LPs, as budgeted. The audit re-ran them and got the same values.
- `ef2265c1` differs from `e8a54520` only in metabolite notes. Neither differs from `8ae6c93b` in
  stoichiometry, bounds or objective, so every LP is the same problem on all three (**provable**).
- The TPP task's knockout results therefore carry over. No knockout sweep was run.

## 6. Tests

- **Focused test (verified):** `tests/test_thiamine_phosphate_chemistry.py` passes: 5 tests, 29
  subtests. It covers:
  - 13 stale-model cases, each failing closed with the model unchanged, including:
    - TPP requirement missing;
    - R752's H+ removed;
    - a new HET-P source;
    - a pH 7.3 formula already set;
    - ATP migrated to pH 7.3;
    - only one target changed;
  - 7 scope rejections, including `schema_version: 1.0`;
  - 9 spec-value rejections, including the audit's tampered specs;
  - exact edits;
  - the expected balance set, with R752, R408 and R750 at {H: 1, charge: 1}, and no non-boundary
    reaction losing balance;
  - idempotence and an SBML roundtrip of the new annotations and notes;
  - the flag default and its requirement of `--tpp-requirement`.
- **Full suite, per test (verified, [logs/pytest_summaries.txt](logs/pytest_summaries.txt)):**

  | Run | Passed | Failed | Errors | Skipped | Subtests |
  | --- | --- | --- | --- | --- | --- |
  | Before | 227 | 11 | 17 | 31 | 562 |
  | Final code | 232 | 11 | 17 | 31 | 591 |

  - The only per-test difference is the 5 new passes.
  - An intermediate run on the pre-audit code also differed only by these tests.
  - The "before" run is the TPP task's final run on identical code (its code files were last changed
    before that run started; the audit did not re-derive this).

## 7. Independent audit and responses

The audit ([AUDIT.md](AUDIT.md)) was a fresh read-only subagent. It audited `e8a54520` and the
pre-fix code.

**Verdicts:**
- **Supported:** C1, C2, C3, C4, C5, C7, C9, C10, C11.
- **Partially supported:** C6 and C8.

**C8 and code defects c1, c2 and c4, all fixed:**
- **c1, after-values not pinned:** the reviewed chemistry is now pinned in the code (`REVIEWED`). The
  audit's tampered specs are now rejected; each is a test case:
  - a dTMP formula for TMP;
  - HET-P at its pH 7.3 form;
  - m971 relabelled to another molecule with charge −2;
  - `charge: true`;
  - `schema_version: 1.0`.
- **c2, partner chemistry not pinned:** the 24 partner species' formula and charge are pinned. ATP at
  its pH 7.3 form is now rejected (test case).
- **c4:**
  - integer checks exclude booleans;
  - missing keys raise `ValueError`.
  - The coupling to the TPP notes on R1342 remains. It is fail-closed by design: if the TPP data file
    is reworded, this one must be regenerated.

**Code defect c3 (H+ lock): not fixed; it is a scope decision** (§8).

**C6 and overstatements a1–a7, all corrected:**
- **a1, R1343/R1344:** now described as a 2-H discrepancy that persists in every protonation state,
  with the ribitol-5-P annotation of `m889` noted, in TASK, REPORT and the HET-P note.
- **a2, HMP-PP:** no longer claimed to be `component_review`.
- **a3, "none worse":** now "no balanced reaction became unbalanced", with R752's charge component
  stated.
- **a4, "(Rhea)":** removed for R750, in the notes and the curation README.
- **a5, "spec values":** now enforced (c1).
- **a6, the Strategy B origin:** labelled inferred.
- **a7, wording:** each molecule's rationale now cites only its own reactions, and the BiGG IDs are
  labelled "through MetaNetX".

**Not changed:** the short mirror annotation keys (`keggc`, `biggm`). The audit found only one tool
reads them, and it tolerates their absence.

## 8. Limits, open findings and next gates

**The H+ lock (needs a user decision):**
- The final chemistry gate (`balance_protons_and_water`, `main.py` ~608) runs before this step.
- The audit verified that, run on this candidate, the gate would remove the H+ from R752, R408 and
  R750, because each residual is exactly one H+.
- Today, step order is the only protection.
- The existing `lock_proton_water_stoichiometry` flag, which `patches.py` sets once, would make the "no
  H+ edits" decision durable. Adding it to those three reactions is a model change outside this
  task's approval.

**The pH 7.3 migration:**
- The KEGG forms are an interim convention, and these species, and HMP-PP, are not in
  `metabolite_microspecies.csv`.
- **Verified by the audit:** the migration audit for the ATP, ADP, AMP, phosphate and diphosphate
  component now lists R1337, R1342, R751, R753 and R754 as regressions (292 against 287), because the
  closure cannot reach the thiamine species. They have to be added to that component by hand. The
  component was not ready before either.
- If the partners migrate first, this step now refuses to apply (partner pin).

**Other limits:**
- Identity rests on reactions and database IDs.
- Opt-in only; the step requires `--tpp-requirement`.
- The TPP step would reject being re-applied to this model, because it pins m975's formula as `None`.
  In a build it runs first.

**Found, not changed:**
- **R754:** its KEGG, MetaNetX and Rhea annotations (R10712, MNXR114031, rh:47848) describe the
  carboxy-thiazole-P reaction, not the model's HET-P equation.
- **R1343 and R1344:** a 2-H discrepancy, and `m889`'s ribitol-5-P annotation, as above.
- **R753's EC 2.7.6.2** comes from the MetaNetX MNXR104888 union "2.7.4.15;2.7.6.2"; KEGG R00616 is
  2.7.4.15. Whether R753 is real in yeast was not checked.
- **Names:** "TMP_", "TDP_" and "TTP_" are unchanged.
- **Other MetaNetX EC unions** were not scanned.
