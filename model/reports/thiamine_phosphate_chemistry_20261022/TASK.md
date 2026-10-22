# Thiamine phosphate chemistry: TMP identity, TPP/TTP/HET-P formulas, R1342 EC (2026-10-22)

## Authorization

User, 2026-10-22, after reading the plan in this session:

> Approved: implement the thiamine phosphate chemistry plan as written: KEGG fully protonated forms,
> m971 and m643 included, R1342 EC 2.7.4.9 removed, names unchanged, on top of 8ae6c93b with
> --thiamine-phosphate-chemistry requiring --tpp-requirement, no H+ edits, WT spot check only. No
> commit or push.

This follows the [TPP requirement task](../tpp_requirement_20261022/REPORT.md), whose findings listed
these chemistry errors.

## Starting state

- **Branch:** `claude/science-test-gap-fix-13f440` at `abe7e40`. The TPP requirement task's work is
  uncommitted in the same worktree.
- **Baseline model:** `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_v2_20261022.xml`,
  SHA-256 `8ae6c93b0ec41cde3db206f50ee9ad2549b1f96a3114d4f9df96a0d75bf01a76`.
  - It is `cf0825b8` plus `--tpp-requirement`.
  - It is a candidate, not a reference.
- **Before-state:** [before.json](before.json), written by [scripts/state.py](scripts/state.py). It
  records the 7 species, every reaction on them, and each reaction's mass/charge balance.

## Targets

**Species.** Formula and charge go to the KEGG fully protonated form, and the annotation is replaced.
Names and compartments are unchanged.

| Species (name) | Before | After | New cross-references |
| --- | --- | --- | --- |
| `m977[C_cy]`, `m973[C_ex]` ("TMP_") | C6H4Cl3NO, 0; annotated as 3,5,6-trichloro-2-methoxypyridine (MNXM1137754, envipath, InChI) | C12H18N4O4PS, +1 | KEGG C01081, CHEBI:9533, MNXM662, BiGG thmmp |
| `m975[C_cy]`, `m971[C_ex]`, `m1116[C_mi]` ("TDP_") | no formula, 0 | C12H19N4O7P2S, +1 | KEGG C00068, CHEBI:9532, MNXM256, BiGG thmpp |
| `m976[C_cy]` ("TTP_") | no formula, 0; BiGG ttp, MNXM148527 (formula-less) | C12H20N4O10P3S, +1 | KEGG C03028, CHEBI:9534, MNXM1341, BiGG thmtp |
| `m643[C_cy]` (HET-P) | no formula, 0 | C6H10NO4PS, 0 | KEGG C04327, CHEBI:17857, MNXM960, BiGG 4mpetz |

**Reaction.** `R1342`: the annotation `ec-code` goes from `['2.7.4.16', '2.7.4.9']` to `2.7.4.16`.
Nothing else on it changes.

- No stoichiometry, bound, GPR, name or compartment changes.
- Curated data: `model/curation/chemistry/thiamine_phosphate_chemistry.json`.
- Guarded step: `apply_thiamine_phosphate_chemistry()` in
  `platform/scripts/gem_annotate/thiamine_phosphate_chemistry.py`, under the opt-in flag
  `--thiamine-phosphate-chemistry`. The flag requires `--tpp-requirement` and runs right after it.
- Focused test: `platform/tests/test_thiamine_phosphate_chemistry.py`.

## Evidence (checked during planning, 2026-10-22)

- **Formulas (verified):** KEGG REST `get/cpd:`:
  - C01081 thiamin monophosphate, C12H18N4O4PS;
  - C00068 thiamin diphosphate, C12H19N4O7P2S;
  - C03028 thiamin triphosphate, C12H20N4O10P3S;
  - C04327 4-methyl-5-(2-phosphooxyethyl)thiazole, C6H10NO4PS;
  - C00378 thiamine, C12H17N4OS. This matches the model's thiamine `m974`, which has charge +1.
- **Charges (verified):** ChEBI through OLS4:
  - CHEBI:9533 thiamine(1+) monophosphate, +1;
  - CHEBI:9532 thiamine(1+) diphosphate, +1;
  - CHEBI:9534 thiamine(1+) triphosphate, +1;
  - CHEBI:17857 4-methyl-5-(2-phosphonooxyethyl)thiazole, 0.
- **MetaNetX 4.5 (verified; local copy `chem_prop.tsv` `8463f8ee…`, `chem_xref.tsv` `ba2451ff…`):**
  - KEGG C01081, C00068, C03028 and C04327 map to MNXM662, MNXM256, MNXM1341 and MNXM960. Their pH 7.3
    forms are −1, −2, −3 and −2 (CHEBI:37575, 58937, 58938, 58296).
  - MNXM1137754 is named "TMP" and is C6H4Cl3NO from envipath. The model's TMP got it by exact-name
    match (Strategy B in `metabolites.py`).
  - MNXM148527 "Ttp" comes from BiGG `ttp` and has no formula.
- **R1342 EC (verified):**
  - MetaNetX `reac_prop.tsv` (`8582cc18…`) lists MNXR104907 with EC "2.7.4.16;2.7.4.9".
  - Its cross-references in `reac_xref.tsv` (`2610e299…`) are all thiamine-phosphate kinase: Rhea
    15913, KEGG R00617, MetaCyc THI-P-KIN-RXN, BiGG TMPK, SEED rxn00438/rxn38718, SABIO 1007.
  - ExPASy `enzyme.dat` (release of 28-Jan-2026, `6df8cb88…`): 2.7.4.9 = dTMP kinase; 2.7.4.16 =
    thiamine-phosphate kinase.
- **Why the KEGG forms (verified, in-memory scan on `8ae6c93b`):**
  - The model keeps ATP, ADP, AMP, phosphate, diphosphate and HMP-PP neutral. The pH 7.3 migration of
    the first five is `component_review` in `metabolite_microspecies.csv`; HMP-PP has no row there.
    (Corrected after the audit; the first version listed HMP-PP as `component_review` too.)
  - Over the 11 non-boundary reactions on these species, the KEGG forms balance 6 (today: 1), and no
    balanced reaction becomes unbalanced. The pH 7.3 forms balance 2.
  - If every partner were also at pH 7.3, all 9 thiamine reactions (R1337, R1342, R408, R750, R751,
    R752, R753, R754, R933) would balance as written. So the H+ in R752, R408 and R750 is correct for
    the target convention and is not edited.

## Expected effect

- **Balance:** R1337, R1342, R751, R753 and R754 become balanced, and R933 stays balanced.
- **Still off:**
  - R752, R408 and R750 each stay one H+ off (the convention mix).
  - R1343 and R1344 stay off in every consistent protonation state of their species: a 2-H discrepancy
    in their lumped equations, involving `m889` and `m38`.
  - The first version blamed only the stored forms of `m889` and `m38`; the audit showed that is
    incomplete.
- **No balanced reaction becomes unbalanced.** Per component, R752's charge residual goes from 0 to +1
  while its mass residual shrinks from six elements to one H (wording corrected after the audit).
- **The LP is unchanged:** formulas, charges and annotations are not in the LP.

## Acceptance criteria

1. pytest per test equals the baseline run plus the new test file passing. The baseline is the TPP
   task's final run on identical code, copied to `logs/pytest_before_results.txt`.
2. The 8ae6c93b argv without the new flag reproduces `8ae6c93b…` byte for byte. With the flag, seeds 0
   and 1 give the same bytes.
3. The whole-model diff against `8ae6c93b` holds only:
   - the formula, charge, annotation and notes of the 7 species;
   - the R1342 annotation and notes.

   Any stoichiometry, bound, GPR, objective or other difference stops the work.
4. Per-reaction balance matches the expected effect above, and no reaction gets worse.
5. WT growth on the Task A medium + PO1f (GLPK) is identical on both models.

## Compute budget and stop condition

- **Budget:** 2 LPs (WT on both models). No knockout sweep: the optimization problem is unchanged
  (criterion 3).
- **Stop and report to the user before continuing on:**
  - a guard rejection;
  - a diff outside the targets;
  - seed-dependent bytes;
  - a test regression;
  - a reaction that becomes less balanced;
  - a WT difference.

## Not in scope

- Renaming "TMP_", "TDP_" and "TTP_".
- The pH 7.3 migration and registering these species in `metabolite_microspecies.csv`.
- The H+ in R752, R408 and R750.
- R1343 and R1344 and their species `m889` and `m38`.
- Whether R753 (TPP → TTP on THI80) is real.
- Other MetaNetX EC unions.
