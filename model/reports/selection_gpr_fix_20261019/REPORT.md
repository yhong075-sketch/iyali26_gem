# Selection-overwritten GPRs fixed: R671, R159/R190 restored; R765/R766 builder correction

**New candidate:** `model/candidates/E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_20261019.xml`
- SHA-256 `51c6fabd486f4e1a9d6365bf7740bde3eca20fa05f1d5454a73586d3276db3fa`.
- Built from c98961e9's flags plus two new opt-in flags.
- Not a reference model.
- Later the same day the user authorized 12 LPs on it and on c98961e9 (§7).

Scope, authorization, baseline and post-audit corrections: [TASK.md](TASK.md). Evidence for each
decision: [the review](../selection_gpr_review_20261019/REPORT.md) §3. Commands:
[commands.txt](commands.txt).

**Superseded first build:** `E5_coq9_resp_pkg_r132_selection_gpr_fix_20261019.xml` (`fa859429…`).
It has the same optimization problem. Only the R671 and R190 evidence notes differ: R671's stated a
wrong EC provenance, and R190's had a KEGG wording slip (§6). It is kept until the user decides
whether to delete it.

## 1. What changed

| Item | Reaction | c98961e9 (selection "after") | New candidate |
| --- | --- | --- | --- |
| A | R671 prephenate dehydrogenase (NADP) | `YALI1F23441g or YALI1E11370g` | `YALI1F23441g`; equation unchanged |
| C(i) | R159 | ATCase only: carbamoyl-P + Asp → carbamoyl-Asp + Pi; `YALI1C33005g or YALI1D09420g or YALI1E11768g` | net channelled CPS + ATCase: Gln + 2 ATP + HCO₃⁻ + H₂O + H⁺ + Asp → Glu + 2 ADP + 2 Pi + carbamoyl-Asp; `YALI1E11768g` |
| C(i) | R190 | `YALI1C33005g or YALI1D09420g or YALI1E11768g` | `YALI1C33005g and YALI1D09420g`; equation unchanged |
| B | R765/R766 | `YALI1E07744g or YALI1D02625g` | unchanged. The builder no longer removes YALI1E07744g before the selection |

- Each restored reaction gets four notes: restore ID, evidence, evidence status and sources.
- The selection's `metadata_previous_*` notes stay as history.

Unchanged:
- genes, metabolites, formulas, charges and bounds;
- names, EC numbers and subsystems (R159/R190 already carried EG-GPR-001's);
- the selection record.

**Whole-model diff against c98961e9** ([model_diff.json](model_diff.json)) has 19 differences:
- 7 reaction fields: R159 GPR, equation and notes; R190 GPR and notes; R671 GPR and notes.
- 4 derived gene–reaction links:
  - YALI1E11370g: R671 → none;
  - YALI1C33005g and YALI1D09420g: R159 + R190 → R190;
  - YALI1E11768g: R159 + R190 → R159.
- 8 derived metabolite–reaction links of R159's species: m10, m32, m50, m130, m141, m143 and m199
  gain R159; m325 loses it.

Nothing else differs. The diff script compares every reaction, metabolite and gene field, the links,
groups, compartments, the objective, and model notes and annotation.

## 2. Mechanism (code and data)

**Data**, in `model/curation/gpr/`. For each reaction:
- the overwritten state (= the selection's "after");
- the reviewed state (= the selection's "before");
- the genes of the restored GPR, with whatever identifiers the model gene carries;
- evidence tier, status and sources;
- formula, charge and compartment of every species involved.

The two files:
- `selection_overwritten_gpr_restore.json`: change ID `SELECTION-GPR-RESTORE-R671-20261019`.
- `selection_cpa_ura2_restore.json`: change ID `SELECTION-CPA-URA2-RESTORE-20261019`. It also records
  what the change relies on but does not edit:
  - `context`: R159/R190 names and EC numbers (EG-GPR-001's), and R607's full definition;
  - `pool_context`: apart from R159/R190, only R607 uses m325, with coefficient −1.

YALI1E11768g carries no identifiers in the model, so for URA2 only the gene's existence can be
checked. Its identity is recorded in the evidence text: NCBI Gene 2912535, KEGG K11541.

**Step:** `platform/scripts/gem_annotate/selection_gpr_restore.py`, generalized from the 2026-10-19
module.
- A registry `CHANGE_SETS` maps each change ID to its data file and exact targets. The 2026-10-19
  entry is unchanged.
- A row may list several `genes`; the restored GPR's gene set must equal them. The old single-`gene`
  form still works.
- A restored equation must balance. One exception: an equation the step leaves unchanged may keep the
  residual its record states (`unchanged_equation_residual`).
  - R671 uses this for its pre-existing +1 charge.
  - Its source is inferred: m41 NADP⁺ is stored as C21H29N7O17P3 with charge 0.
  - An edited equation that carries a residual is refused.
- Context reactions and pool users are checked before any edit.

All 2026-10-19 guards are unchanged:
- scope and species;
- gene existence and identity;
- note conflicts, compared after unescaping;
- each reaction must be exactly in its cleared or its restored state.

Any mismatch raises `ValueError` with no edits. A rerun returns `already_correct`.

**Wiring:**
- `cli.py` adds `--restore-selection-overwritten-gprs` (A) and `--restore-selection-cpa-ura2` (C).
  Both are off by default.
- `main.py` calls each after the 2026-10-19 restore, through `restore_after_selection`. That function
  refuses an incomplete selection.
- Build-record keys: `reaction_selection.selection_gpr_restore_r671` and
  `reaction_selection.selection_cpa_ura2_restore`.

**Builder correction (B):**
- R765/R766 are removed from `_GPR_MISANNOTATION_REMOVALS` in `patches.py`. A comment records the
  YALI1/YALI0 locus mix-up. R671 stays.
- The matching `main.py` comment was corrected.

**Tests:** `platform/tests/test_selection_gpr_restore_change_sets.py` has 10 tests. All run on
in-memory models except the a56f3649 one.
- Restore and idempotence for both change sets, including R671's recorded residual.
- Knockout propagation:
  - YALI1E11370g closes nothing; YALI1F23441g closes R671;
  - either CPA gene closes R190 only; YALI1E11768g closes R159 only.
- m325 is left with R190 and R607 only.
- A rerun after an SBML save and reload is a no-op.
- 13 stale-input cases. Each test asserts the guard message it expects:
  - R671: GPR, equation, NADP charge, gene identity, a conflicting note;
  - context: R607 equation, R607 GPR, R159 name, R190 EC;
  - a mixed R159 state, CPA1 identity, a missing R190, an extra m325 user.
- 6 malformed records, each with its expected message.
- Both flags are opt-in and need a complete selection.
- On the real a56f3649 file, applying all three change sets changes only the 7 targets.
- `remove_misannotated_gprs` keeps YALI1E07744g on R765/R766 and still removes YALI1E11370g from R671.

The 2026-10-19 test file is unchanged and passes against the generalized module.

## 3. Verification

| Check | Result | Level |
| --- | --- | --- |
| Control build: a56f3649's flags, final code (B included) | `a56f3649…b58515`, byte-identical to a56f3649 (v1 and v2 controls). B changes no model on this build path | verified |
| Control build: plus `--restore-selection-cleared-gprs`, final code | `c98961e9…ccae4c`, byte-identical to c98961e9 (v1 and v2 controls). The generalization leaves the 2026-10-19 change set unchanged | verified |
| Candidate v2 build | `51c6fabd…b3fa`. Log: "Selection GPR restore / R671 GPR restore / CPA/URA2 restore: applied". Selection `complete: true` | verified |
| Whole-model diff vs c98961e9 | 19 differences, all intended (§1). The auditor's own comparison gave the same 19 for v1 | verified |
| v1 → v2 | 2 differences: the R671 and R190 evidence notes | verified |
| Balance ([static_checks.json](static_checks.json)) | R159 `{}`. R190 `{}`. R671 `{charge: +1}` (unchanged equation, as recorded). R765/R766 keep their pre-existing `{charge: +2}` from the selection's equation | verified |
| m325 pool | made by R190, used by R607; nothing else | verified |
| Knockout propagation (static `gene.knock_out()`, no LP) | YALI1F23441g → R671. YALI1E11370g → none. YALI1C33005g → R190. YALI1D09420g → R190. YALI1E11768g → R159. YALI1E07744g, YALI1D02625g → none (OR) | verified |
| Selection rows R765/R766 | GPR field `already_correct` (was `applied` before B); stoichiometry `applied` | verified |
| Focused tests | 16 passed (6 from 2026-10-19 + 10 new), 43 subtests | verified |
| Full pytest ([v1 comparison](logs/pytest_comparison.txt), [v2 comparison](logs/pytest_comparison_v2.txt)) | Before: 167 passed, 8 failed, 14 errors, 38 skipped. After (final code): 177 passed, 8 failed, 14 errors, 38 skipped. Per test: identical except the 10 new tests, which pass. The 22 failures and errors are pre-existing | verified |
| Build logs | One `LibSBML error code -3` line on R_R304 notes in each build, as in every earlier build | observed |

## 4. What this means, and what was not tested

- **A:**
  - YALI1E11370g (GatB) no longer links to any reaction.
  - A YALI1F23441g knockout now closes R671, the only prephenate → 4-hydroxyphenylpyruvate step in
    the model.
  - On SD-Leu + PO1f the knockout gives KO/WT 0.421 (§7), non-essential at every cutoff.
- **C(i):**
  - Free carbamoyl phosphate (m325) now serves arginine only (R190 → R607). Pyrimidines get carbamoyl
    phosphate only inside R159.
  - CPA2 and CPA1 knockouts are equivalent: both close R190 only. Each gives KO/WT 0.1439 (§7).
  - A URA2 knockout closes R159, the only forward producer of carbamoyl-aspartate (m326). R276, a
    reversible dihydroorotase, is the other reaction that uses m326. On SD-Leu + PO1f, which supplies
    uracil, growth is unchanged (§7).
  - The open question stays open: does *Y. lipolytica* URA2 also feed the arginine pool? The notes
    copy EG-GPR-001's dossier status, `needs_more_evidence`, and claim no new validation.
  - In the 2026-09-10 screen data, CPA2 is on the consensus essential list and CPA1 is not. That is
    curation data, not independent validation, and is not used here.
- **B:**
  - No model change on the builds checked: the E5 chain with and without the 2026-10-19 flag,
    byte-identical controls.
  - For other build paths this is **inferred**. `main.py` has one build path, and it runs the
    selection unconditionally, which sets the same GPR.
  - The fix stops a wrong patch from reappearing if the selection's R765/R766 row is ever removed.
  - YALI1D02625g's OR membership (a 347-aa fragment) is unchanged. It still lacks evidence of
    independent catalysis; that is a separate identity question.
- **Not tested:**
  - any condition other than SD-Leu + PO1f;
  - flux distributions (only objective values were read);
  - FVA;
  - the default (E0) build path with the new flags;
  - a full-genome screen.

## 5. Limits and next gates

- **Guard gaps the audit found:**
  - YALI1E11768g has no model identifiers, so its identity cannot be checked.
  - `restore_after_selection` does not compare the record's `selection_id` with the selection that
    ran. The 2026-10-19 module had the same gap, and its test passes a selection dict without an ID.
  - The gene-set check cannot tell AND from OR. The rule itself comes from the curated `restored`
    state.
  - An unknown change ID raises `KeyError` before any edit.
- **Location:** worktree `claude/clever-kare-7d2e1a`, uncommitted. It holds byte-identical copies of
  the uncommitted 2026-10-19 restore files from `claude/orphan-reaction-newest-model-b8ddc8`.
  `cli.py`, `main.py` and `selection_gpr_restore.py` are extended here. Merging the two worktrees
  needs care.
- **Records not edited (out of scope):**
  - The 2026-10-19 restore REPORT §6 still calls the E07744g removal a "verified" mis-annotation.
    That report is in the other worktree.
  - The selection's R765/R766 "before" value (`YALI1D02625g`) no longer matches the builder. The
    selection still completes, because the "after" value is allowed.
  - `patches.py` lines 917–918 still say R671 "carries a self-contradictory EC 6.3.5.7". Current
    builds annotate R671 with 1.3.1.78/1.3.1.13. This is a pre-existing comment.
- **Interpreter:** main checkout `.venv`, Python 3.13.5, cobra 0.30.0.
- Nothing is committed or pushed.

## 6. Independent audit

A fresh read-only reviewer was given TASK.md, the diff, the logs and the deliverables. It recomputed
hashes, re-derived the diff with its own script, probed the guards, re-ran the focused tests, parsed
both JUnit files and opened every cited source. It audited the v1 candidate.

| Claim | Verdict |
| --- | --- |
| 1. Controls byte-identical | supported |
| 2. Candidate diff confined to the targets | supported. Its own comparison also lists the 8 derived metabolite links, which REPORT v1 omitted |
| 3. Restored = selection "before", cleared = "after"; balance | supported |
| 4. Fails closed, idempotent, opt-in; 2026-10-19 unchanged | supported, with guard gaps (below) |
| 5. pytest per test | supported |
| 6. Builder correction B | supported (also `candidate_stage_reference_sha256` unchanged across all builds) |
| 7. Identity and evidence claims | supported. Nit: the K01955 wording in the R190 note |
| 8. Stayed in scope | supported |
| 9. REPORT meaning and limits | partially supported (wording, below) |

**Fixed after the audit** (rebuilt as v2; both controls re-run with the final code):
- **Wrong provenance in the R671 note.** It said E11370g entered R671 "on R671's wrong source EC
  6.3.5.7". The source model has no 6.3.5.7. The EC sits on the isozyme table's R671 row and in the
  legacy `model.xml`. The data was corrected, and so were the review report and the memory note.
- **m325 exclusivity was claimed but not guarded.** Added `pool_context`, plus a test case with an
  extra m325 user.
- **The stale-input tests did not pin which guard fires.** Each case now asserts its message.
- **The diff omitted metabolite links.** The diff script now compares them; §1 lists 19 differences.
- **K01955 wording in the R190 note** corrected.
- **Report wording:** the m325 claim, the "each confirmed" claim, the URA2-knockout meaning, the
  "no model change" scope and the NADP⁺ inference were all corrected or labelled.

**Not changed (recorded as limits in §5):**
- URA2 identity is checked by existence only;
- `selection_id` is not compared;
- an unknown change ID raises `KeyError`;
- the stale `patches.py` R671 EC comment.

## 7. Solves (authorized 2026-10-19 after the build)

Details:
- Authorization and condition: [TASK.md](TASK.md), last section.
- Script: [scripts/solves.py](scripts/solves.py), the 2026-10-19 script with only the gene list
  changed.
- Raw values: [solves.csv](solves.csv), SHA-256
  `fa35779a32fa65281d7387ca7f2531dea80c2f5db5a05b3a8879342c7c85c145`.
- Log: [logs/solves.log](logs/solves.log).

Settings:
- GLPK with presolve;
- SD-Leu with the whole medium replaced;
- PO1f profile: R612 closed, uracil uptake override.

Exactly 12 LPs were run, all `optimal`. No other solve, FVA or pFBA.

| Case | c98961e9 growth (h⁻¹) | KO/WT | New candidate growth (h⁻¹) | KO/WT | Reactions the KO closes (new) |
| --- | --- | --- | --- | --- | --- |
| WT | 1.4308057353528618 | 1 | 1.4308057353528527 | 1 | |
| KO YALI1F23441g | 1.4308057353528618 | 1.0 | 0.602199505837895 | 0.42088138938678904 | R671 |
| KO YALI1E11370g | 1.4308057353528618 | 1.0 | 1.4308057353528527 | 1.0 | none |
| KO YALI1C33005g (CPA2) | 1.4308057353528618 | 1.0 | 0.20595239043348762 | 0.1439415466018505 | R190 |
| KO YALI1D09420g (CPA1) | 1.4308057353528618 | 1.0 | 0.20595239043348762 | 0.1439415466018505 | R190 |
| KO YALI1E11768g (URA2) | 1.4308057353528618 | 1.0 | 1.4308057353528527 | 1.0 | R159 |

**What the values show** (observed in runs):
- **WT growth is unchanged.** The difference is 9e-15.
- **On c98961e9 no knockout changes growth.** Each gene was covered by an OR partner there.
- **New candidate:**
  - The YALI1F23441g knockout reduces growth to 42%. With R671 closed, growth must rely on another
    tyrosine supply. That this is the medium's tyrosine uptake (R1233, 0.02486 in `sd_leu.csv`) is
    **inferred**; fluxes were not read.
  - The CPA2 and CPA1 knockouts both give 0.20595239 h⁻¹. That equals 0.02586 / 0.125563 =
    0.20595239043348756: the medium's arginine uptake (R1202) divided by the biomass arginine
    coefficient recorded for the 2026-09-15 R190-closure test, which gave the same growth.
    Arginine-uptake limitation is therefore the likely mechanism. This is **inferred** from the
    arithmetic; the coefficient was not re-read in this model and fluxes were not read.
  - The URA2 knockout leaves growth unchanged, because PO1f takes up uracil (profile override 1000).
    This is **inferred**; no flux was read.

**Calls at the four cutoffs.** A gene is essential when KO/WT is strictly below the cutoff.

| Gene | New KO/WT | 1% | 5% | 10% | 15% (primary) | c98961e9, all cutoffs |
| --- | --- | --- | --- | --- | --- | --- |
| YALI1F23441g | 0.4209 | non-essential | non-essential | non-essential | non-essential | non-essential |
| YALI1E11370g | 1.0 | non-essential | non-essential | non-essential | non-essential | non-essential |
| YALI1C33005g | 0.1439 | non-essential | non-essential | non-essential | **essential** | non-essential |
| YALI1D09420g | 0.1439 | non-essential | non-essential | non-essential | **essential** | non-essential |
| YALI1E11768g | 1.0 | non-essential | non-essential | non-essential | non-essential | non-essential |

**Against the experimental list.** These rules follow `model/expected/benchmark.md`:
- Labels come from `essentiality_per_gene.tsv` of the 2026-09-10 screen (consensus list SHA
  `1e887f5a…`).
- A model gene on the consensus list is a positive.
- A model gene not on it counts as non-essential: a user-defined negative class.

| Gene | Label | c98961e9 at 15% | New at 15% | New at 1/5/10% |
| --- | --- | --- | --- | --- |
| YALI1F23441g | negative (not on list) | TN | TN | TN |
| YALI1E11370g | positive | FN | FN | FN |
| YALI1C33005g (CPA2) | positive | FN | **TP** | FN |
| YALI1D09420g (CPA1) | negative (not on list) | TN | **FP** | TN |
| YALI1E11768g (URA2) | negative (not on list) | TN | TN | TN |

For these five genes at 15%:
- +1 TP and −1 FN (CPA2);
- +1 FP and −1 TN (CPA1).

At 1/5/10% nothing changes. No other gene was re-screened, so the change in model-wide counts is
these two calls only if no other gene's call moved. The GPR changes touch only these genes'
reactions, so that is expected, but it was **not run**.

These labels were used to build and curate the model, so this comparison is not independent
validation. The change was made on gene-identity evidence, not to improve the match. The CPA1 false
positive is the inconsistency the review already noted: a strict glutamine-only AND gives CPA1 and
CPA2 the same phenotype, but only CPA2 is on the list.
