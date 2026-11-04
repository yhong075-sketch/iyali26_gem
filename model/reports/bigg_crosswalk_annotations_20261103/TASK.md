# TASK: reviewed BiGG-crosswalk annotations into iYali26, plus a BiGG-ID export (2026-11-03)

## Authorization (user, this session, 2026-11-03)

- On Andrew McNaughton's lipid-unlump handoff (`mcnaugh_pass`): "plz check it and try to merge in our model".
- Scope: "IDs + annotations". First: "I think we need to consider which one is correct first."
- On the IDs: "Keep ours + export". Compartment codes for exported IDs: "yeast-GEM codes".
- On the plan below: "Approve plan". It writes 898 values plus 1 removal, holds 80 + 14 rows and fills no
  one-sided keys.

The scope does not cover charges, formulas, stoichiometry, GPRs, IDs in the working model, or Andrew's files
(colleague-owned, read only).

## Baseline

- **Model:** `model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_20261101.xml`,
  SHA-256 `2bc263f23603c27131d0b07143300a484dbd5a135ef7151665574400ab64a062`.
- **Reproduced:** rebuilt at HEAD `b57272d3aaf1514955241476ea7feb3a5d35b6b5` with the argv recorded in its
  `.build.json`. The result is byte-identical. The original build ran at `b83e13a9` with a dirty tree; the
  reproduction shows that HEAD yields the same file.
- **Crosswalk** (read only, untracked colleague folder `iyali26_gem_lipid_unlump/mcnaugh_pass/handoff/`):
  - `metabolite_crosswalk.csv`: `86970fa84ac03ae7040e670b6ab1d8c0c5c1501ab39adaa4dee7499132e9cf64`
  - `reaction_crosswalk.csv`: `743b2680d69b883c17a98077818b674c46facd8fd9172a879528ee9d6adf1037`
- **Python:** `iyali26_gem/.venv` (3.13.5, cobra 0.30.0, optlang 1.8.3).

## Targets

There are 1,133 annotation conflicts between our model and the crosswalk. Each was judged by a first reviewer and
an independent skeptic, with a tie-breaker where they disagreed. The adjudication is under
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/bigg_crosswalk_annotations_20261103/adjudication/`.

Curated data (`model/curation/annotation/bigg_crosswalk_annotation_decisions_20261103.csv` and `.json`) holds 899
rows on 595 entities:

| Action | Rows | What it means |
| --- | ---: | --- |
| `write_his` | 898 | Andrew's value replaces ours: take-his or either rows whose reviewer-confirmed value is exactly his. |
| `drop_wrong` | 1 | R2069 `bigg.reaction`: remove `r0393` and `R_r0393`, which the reviewers judged wrong, and keep `GLYALDDr` and `R_GLYALDDr`. |

The starting value of every target is in `before.json`.

Not written:
- 80 rows where the reviewers recommend an ID other than Andrew's;
- 14 rows from 7 value pairs judged differently on different entities;
- 101 keep-ours rows;
- 22 charges-task rows;
- 15 third-value rows;
- 2 meeting rows;
- 1,610 keys that only Andrew has, which were never reviewed.

The 80, 14, 15 and 2 go to the meeting with Andrew.

## Expected effect

Annotation values change only on the 899 (entity, key) targets. Nothing else in the model changes: names, IDs,
formulas, charges, stoichiometry, bounds, GPRs, notes, genes, compartments and objective. WT growth stays the
same.

## Acceptance criteria

1. The guarded `apply_bigg_crosswalk_annotations` behaves as follows:
   - it applies all 899 rows;
   - a second run reports `already_correct`;
   - it raises without editing anything on stale, partial, duplicate, missing or out-of-scope rows.
2. The new candidate's whole-model diff against the baseline shows exactly the 899 annotation changes.
3. WT `slim_optimize` (model bounds, same solver) is the same as the baseline's.
4. `python -m pytest`: only the new tests differ from the before run (same venv).
5. The BiGG-ID export:
   - changes only IDs and adds `legacy_id` notes;
   - renames with the yeast-GEM codes;
   - keeps the 15 entities whose new BiGG ID the review judged wrong;
   - gives the same growth.
6. An independent read-only audit reports no unsupported claims.

## Budget and stop conditions

**Budget:** three builds (baseline reproduction, new candidate, spare), two pytest runs, and four WT solves
(baseline, candidate, export input and output). No KO sweeps and no medium changes.

**Stop if:**
- the baseline does not reproduce;
- a guard trips on the baseline;
- the diff shows any change outside the 899 targets;
- growth changes;
- an export check fails.

## Amendment (2026-11-03, after the build)

The export's first run stopped at its crosswalk check. The crosswalk predates the candidate:
- 4 reactions have changed equations (R2135, R2145, R2171, `biomass_C`);
- 13 reactions and 1 metabolite are newer than the crosswalk.

None of the 4 has a renamed status, and no decision row touches them. The user chose "Narrow the check":
- a reaction that is renamed must still have the crosswalk's equation, or the export fails closed;
- changed reactions that are not renamed, and entities absent from the crosswalk, keep their legacy IDs and are
  listed in the rename report.

The user also chose "yeast-GEM codes" (`em` becomes `erm`, `en` becomes `ce`). A sibling session reports that
Chalmers iYali5 uses `em`/`en` for the same compartments. This session has not verified that, so the point stays
open for the meeting with Andrew.

## Amendment 2 (2026-11-03, after the independent audit)

The user answered "Yes, after fixes" to the recommended picks.

- **Decisions:** 896 rows (895 `write_his` + 1 `drop_wrong`) on 595 entities.
  - The 3 NAD `chebi` rows (m27, m116, m122) are held. Their reviewer note names a more exact species, so the rule
    is now exact equality with Andrew's value.
  - Rows held as "reviewers recommend another ID": 83.
- **Export skip list:** 31 entities, every one whose BiGG annotation row was not written. The 15 named before were
  only those whose BiGG ID names another entity.
- **Kept on purpose:** the prenyl-diphosphate pairs (m640/m419, m513/m414). Their different verdicts follow from
  all-trans compounds in the mitochondrion and cytosol versus cis compounds in the lipid particle. "Held pairs"
  means the 7 pairs listed in the generator.
- **What the export writes:** mcnaugh_pass IDs, BiGG IDs where his crosswalk has one.
- **Superseded:** candidate `19cf7cef` (built from 899 rows) and its export are superseded. Their curated inputs are
  kept under `superseded_19cf7cef/` in the artifacts folder.
- **Added compute:** one rebuild, WT and energy checks on the new candidate, the export growth check, the certified
  knockout screen on the final candidate, and a full pytest.

## Coordination

The session "Fix C_em/C_en compartment names in the builder" also edits `main.py` and builds its own candidate.
Whichever lands second rebuilds with both flags.
