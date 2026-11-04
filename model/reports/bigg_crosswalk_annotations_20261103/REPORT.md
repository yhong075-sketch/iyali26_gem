# REPORT: reviewed BiGG-crosswalk annotations, and a mcnaugh_pass-ID export (2026-11-03)

Scope, authorization, baseline and acceptance criteria are in [TASK.md](TASK.md); its two amendments carry the final
numbers. The audit is in [AUDIT.md](AUDIT.md). The large outputs are in
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/bigg_crosswalk_annotations_20261103/` (below: `artifacts/`).

## 1. What changed and why

Andrew McNaughton's lipid-unlump handoff (`mcnaugh_pass`) has a crosswalk that pairs our legacy-ID metabolites and
reactions with his model's IDs and with BiGG / MetaNetX / ChEBI / KEGG annotations. Where its annotation value and
ours disagreed (1,133 rows), each row was judged for which value names our entity. The reviewed values now enter
our model through the builder. IDs stay legacy, and a renamed copy is exported separately.

- **Builder step:** `apply_bigg_crosswalk_annotations` (`platform/scripts/gem_annotate/bigg_crosswalk_annotations.py`),
  opt-in flag `--bigg-crosswalk-annotations`. It is the last step that writes annotations.
- **Curated data:** `model/curation/annotation/bigg_crosswalk_annotation_decisions_20261103.csv` (896 rows on 595
  entities) and its manifest `.json`.
  - 895 `write_his`: Andrew's value replaces ours for that key.
  - 1 `drop_wrong`: R2069 `bigg.reaction` loses `r0393` and `R_r0393`, which the reviewers judged wrong, and keeps
    `GLYALDDr` and `R_GLYALDDr`.
- **New candidate:** `model/candidates/E5_…_ceramide_pool_bigg_annot_v2_20261103.xml`, SHA-256
  `9158e9f2a0d21f654f0a7919eed9f2bacba01aea8b315808d93ec14dd3af4129`. It is the baseline `2bc263f2` plus the flag,
  built at HEAD `b57272d3` (command in [logs/build_command_v2.txt](logs/build_command_v2.txt)).
- **Export tool:** `platform/tools/export_bigg_ids.py`, with policy `model/curation/annotation/bigg_export_policy_20261103.json`.
  It writes a renamed copy, never a candidate. The export of v2 is `artifacts/export/iyali26_mcnaugh_ids_9158e9f2.xml`,
  SHA-256 `181fe59ef8968c0338511047821839d0c49b1754bb5c0ff6f06791a0a7006bc2`.
- **Superseded:** candidate `…ceramide_pool_bigg_annot_20261103.xml` (`19cf7cef…`) and its exports
  `artifacts/export/iyali26_bigg_ids_19cf7cef*.xml`. That candidate carries 899 rows: the 896 plus the 3 NAD ChEBI rows
  held after the audit. Do not use it. Its curated inputs and a v1→v2 decision-ID map are in
  `artifacts/superseded_19cf7cef/`. Decision IDs (`BXA-…`) were renumbered between the two versions.

## 2. How the decisions were made (evidence)

- **Conflicts (verified):** keys where both sides have a value and Andrew's value is not in our list. This covers the
  metabolites Andrew linked, and reactions whose crosswalk status is identical, h_h2o_only or
  proton_translocation_differs. There are 1,133 rows on 743 entities.
- **Adjudication (verified that it ran and covered every row):** 33 batches, each judged by:
  - a first reviewer;
  - an independent skeptic told to refute;
  - a tie-breaker where the two disagreed (19 rows).

  Evidence came from the local MetaNetX 4.5 tables and opened KEGG / ChEBI / BiGG URLs. Each verdict is a reviewer
  judgment, recorded with its evidence in `artifacts/adjudication/final_adjudication.csv`; this round did not redo
  them. Of the 896 written rows, 890 rest on `verified` evidence and 6 on `inferred`, including the R2069 removal.
- **Why our side was often wrong (verified):** an in-memory replay of `annotate_metabolites` on
  `model/source/iyali26.xml` (`artifacts/scripts/replay_strategy.py`) shows that 155 of the 182 MetaNetX values
  judged wrong came from the formula fallback (Strategy C). It takes the first MetaNetX compound with the same
  formula, and the builder then copies all of that compound's cross-references, so wrong ChEBI and KEGG values
  follow the same pick. On the reaction side, the annotator copies every BiGG ID linked to a compartment-generic
  MetaNetX reaction. That is how transports pick up exchange and sink IDs.
- **What is written:** "take his" and "either" rows whose reviewer-confirmed value is exactly Andrew's (an empty
  correct value, or one equal to his), plus the one `drop_wrong`.
- **What is not written:**
  - 83 rows where the reviewers recommend another ID or a more exact species, including the 3 NAD rows;
  - 14 rows from 7 value pairs judged differently on two entities;
  - 101 keep-ours rows;
  - 22 charges-task rows;
  - 15 third-value rows;
  - 2 meeting rows;
  - all 1,610 keys that only Andrew has.

  The prenyl-diphosphate pairs (m640 vs m419, m513 vs m414) are kept on purpose. The mitochondrial and cytosolic
  ones are all-trans, the lipid-particle ones cis-prenyltransferase intermediates, so they are different compounds.

## 3. Code logic

- **Builder step:**
  1. Loads the manifest and checks the CSV against its SHA-256.
  2. Checks the scope: schema, set ID, status `user_authorized`, row count.
  3. Validates every row before any edit:
     - the key must be allowed for the entity type;
     - the action must be known;
     - after must be non-empty and differ from before;
     - the target must be unique;
     - the entity must exist;
     - its current value list must equal `before` (apply) or `after` (already applied).

  Any stale row, or a mix of applied and unapplied rows, raises with no edit and names the stale rows (or the rows
  in the rarer state). Values are rebound, not mutated, because `Model.copy()` can share annotation dicts.
- **Export:**
  1. Checks the pinned crosswalk SHA-256s.
  2. Every crosswalk entity must exist in the model.
  3. Each reaction it would rename must still have the crosswalk's equation, and each metabolite it would rename
     its name, formula and charge. Otherwise it raises.
  4. Maps each linked entity to Andrew's mcnaugh_pass ID, with `em`→`erm` and `en`→`ce` for C_em / C_en. These are
     the yeast-GEM codes; a metabolite there must carry the matching suffix.
  5. Skips the 31 entities on the policy's skip list: every entity whose BiGG annotation row is in the conflict set
     but was not written as Andrew's value. A taken target keeps the legacy ID.
  6. Renames through temporary IDs and adds `notes["legacy_id"]`.
  7. Writes `<name>.unverified.xml`, reloads it, and requires the full model definition to equal the input's under
     the ID map, apart from `legacy_id`. Only then does the file take its final name.

  Unlike Andrew's `translate_to_bigg.py`, it never writes the crosswalk's annotation values. His script would undo
  every keep-ours row (verified by reading its lines 128–131 and 170–171).

## 4. Results (all observed in this session)

| Check | Result |
| --- | --- |
| Baseline `2bc263f2` rebuilt at HEAD `b57272d3` | byte-identical |
| Whole-model diff, v2 vs baseline: every metabolite, reaction, gene, group, compartment and the objective | 896 differences, all decision targets with the decided values, none unexpected ([whole_model_diff.json](whole_model_diff.json)) |
| WT growth, Task A medium + PO1f, checked Gurobi at 1e-9 | baseline and v2 both 0.31205354361020554, optimal, checks passed ([wt_energy.json](wt_energy.json)) |
| ATP with all uptakes closed | 0.0, optimal |
| Certified knockout screen of v2 vs the baseline screen (ceramide round) | 1,081 genes, same KO growth and status for every gene, 0 class changes, the same 20 certified-infeasible knockouts ([ko_compare_vs_2bc263f2.json](ko_compare_vs_2bc263f2.json)) |
| Export: IDs only | the reload equals the input under the ID map (the tool's guard; rechecked by the auditor) |
| Export: metabolites | 1,529 renamed. Kept: 21 skip list (7 "BiGG ID names another entity", 14 "row not written"), 3 target collisions (`iad_c`, `lnlc_c`, `ppcoa_c`), 1 not in the crosswalk (ceramide pool) |
| Export: reactions | 1,818 renamed. Kept: 10 skip list (8 + 2), 3 collisions (`ALDD2xm`, `GLNS`, `UPPRT`), 4 with changed equations (R2135, R2145, R2171, `biomass_C`), 13 not in the crosswalk |
| Export growth under v2's bounds, mapped by legacy ID | 0.31205354361020554, optimal ([export_growth.json](export_growth.json)) |
| pytest, `iyali26_gem/.venv`, default solver (Gurobi) | see §8 |

Screen scores, the same for the baseline and v2. Rule B with mtDNA and untargeted genes out of scope; the negative
class is model genes not on the consensus list (user-defined). This is a development reference, not independent
validation.

| Cutoff | TP / FN / FP / TN | MCC |
| --- | --- | --- |
| 1 % | 114 / 212 / 37 / 701 | 0.396 |
| 5 % | 121 / 205 / 39 / 699 | 0.411 |
| 10 % | 125 / 201 / 44 / 694 | 0.408 |
| **15 %** | **142 / 184 / 44 / 694** | **0.456** |

## 5. What this means

- **The model:** the annotations of 595 entities now name the compound or reaction the review found right. The
  optimization problem is unchanged. That is shown by the diff, and observed in WT growth and in every knockout.
- **The export:** carries Andrew's mcnaugh_pass IDs, not BiGG model IDs.
  - Of the 3,347 renamed entities, 2,415 have an ID built from his crosswalk's BiGG value. 914 have no BiGG value
    there (for example R80 `ACACT9p` becomes `KAT240p`). 18 have an ID that differs from his own BiGG value (for
    example R545 `METS_2` while his BiGG value is `METS_1`).
  - The compartment codes are yeast-GEM's (`_p` is the peroxisome, where BiGG's `p` is the periplasm).
  - The BiGG link is the compartment-free `bigg.metabolite` / `bigg.reaction` annotation.
  - Our C_em and C_en are still named "endosomal membrane" and "endosome" in this export. The names fix is on
    another branch (§7).

## 6. Could the method be better, and how to test it

- **Root cause:** the formula fallback and the copy-every-cross-reference behaviour remain in the builder. Every
  entity outside Andrew's crosswalk keeps its formula-only picks. Candidate changes, each a separate, model-wide
  decision:
  - require a name-token overlap or a single candidate before accepting a formula match;
  - copy only the cross-references that resolve to the chosen MNXM.

  Test either by replaying `annotate_metabolites` before and after, and scoring against the 1,133 adjudicated rows.
- **The 6 export collisions:** these come from duplicate species or reactions in our model that Andrew maps to one
  ID. They belong in a curation round, not in the exporter.
- **A latent gap in the export check:** a metabolite with no mcnaugh_pass link is not identity-checked, yet a
  reaction using it can be renamed, because the equation check compares metabolite IDs. There is no current effect.

## 7. Limits, not tested, next gates

- **Andrew's review:** he has not seen these verdicts. The meeting list (`artifacts/adjudication/for_andrew_review.csv`,
  133 rows), the 3 NAD rows, the 31 skipped entities and the `em`/`en` compartment codes need him. A sibling session
  reports that Chalmers iYali5 uses `em`/`en`; this round did not verify that.
- **Crosswalk age:** the crosswalk predates the candidate (4 changed reactions, 13 newer reactions and 1 newer
  metabolite). Those keep legacy IDs in the export.
- **Unreviewed keys:** the 1,610 keys only Andrew has were neither reviewed nor filled.
- **Charges:** protonation species (22 charges-task rows plus the held species rows) wait for the charges task.
- **Untested paths:** a manifest SHA mismatch in `load_decisions`, and the export's failed-reload path.
- **Other branches:**
  - The compartment-name fix (`--compartment-names-em-en`) is commit `c3a66626` on `claude/xenodochial-clarke-c02b1f`,
    candidate `ebf18a5f`.
  - A gap-fill suffix fix is in progress on `claude/mystifying-hofstadter-c109f5`. It reports 5 mixed-compartment
    gap-fill reactions in the default and `2bc263f2` builds.

  A combined candidate needs one rebuild with the flags together.
- **Compute against budget:** the plan allowed two growth solves and TASK.md four WT solves. This round ran:
  - 2 model-bound solves (uninformative, later removed);
  - 10 Task A LPs: WT and energy on the baseline, v1 and v2 (6), plus the export growth check for v1 and v2 (2
    each);
  - the certified screen on v2 (1,082 sweep LPs, 1 ceiling, 40 certificate LPs);
  - the auditor's 3 WT solves;
  - 3 builds.

  The screen and the rebuild were approved after the audit; the earlier solves went over the stated budget.
- **Publishing:** committed as `048e18ec` on branch `claude/lipid-unlump-collage-zip-merge-a59051`; pushed with the
  mcnaugh_pass reactions round and opened as PR #47 to `main` (2026-11-04).

## 8. Tests

`iyali26_gem/.venv` (Python 3.13.5, cobra 0.30.0), `python -m pytest --continue-on-collection-errors` from
`platform/`. In every run, 4 `test_gem_review_lit_*` modules cannot be collected (no `pydantic_ai` / `httpx2` in this
venv).

| Run | Passed | Failed | Errors | Skipped |
| --- | ---: | ---: | ---: | ---: |
| Before the change (HEAD `b57272d3`) | 426 | 3 | 23 | 24 |
| After v1 | 441 | 3 | 23 | 24 |
| Final run on v1 code | 442 | 3 | 23 | 24 |
| Final run on the fixed code | 443 | 3 | 23 | 24 |

In each comparison with the before run, the only differences are the new tests, all passing: 17 in the final run. No
existing test changed outcome.

After the final run, two files changed:
- the export tool's docstring was corrected (audit F1): tool SHA-256 `d474274c…` became `5d497d9f…`, with no code
  change;
- one test now also covers a policy-supplied skip reason.

The two new test files were rerun: 17 passed. The v2 export was made with tool `d474274c`; a copy of that file is
`artifacts/scripts/export_bigg_ids_d474274c.py`.

## 9. Files

- **Repository:** the code, tests and curated data in §1, plus this folder:
  - [TASK.md](TASK.md), [AUDIT.md](AUDIT.md) and [before.json](before.json);
  - the v2 verification files;
  - [verification_inputs_sha256.txt](verification_inputs_sha256.txt): candidate, export, medium, profile, solver and
    script SHAs;
  - [logs/](logs/): build commands and pytest summaries.
- **Artifacts:**
  - `adjudication/`: conflicts, hints, batches, verdicts, the meeting list, the workbook;
  - `workflow/`: the adjudication workflow script and outputs;
  - `scripts/`: generators, replay, diff and growth scripts;
  - `verification/` (v1) and `verification_v2/`;
  - `export/`, `logs/`, `superseded_19cf7cef/`;
  - `external/`: the iYali download, SHA-256 `c0b54165…`;
  - `SHA256SUMS`.
