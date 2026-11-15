# Reaction cross-reference cleanup through the builder (candidate on 684f3752), 2026-11-13

## 0. Authorization (the user's words, in this session)

- 2026-11-12, answering "Implement the annotator metabolite post-check through the builder here?": "Yes, fix here".
- 2026-11-13, to the plan of a last-step reviewed removal on the 684f3752 flag set (one note per reaction, `ec-code`
  untouched, crosswalk-decided values untouched, no LP): "Approve as planned".
- 2026-11-13, on also removing the values the annotator merged in as "equivalents" of the same wrong pick: "if it's
  about the lipid-unlump, don't touch it, if for other parts, remove."
- 2026-11-13, on the 14 fatty-acid-synthesis reactions that lipid-unlump files do not name: "Exclude them too
  (Recommended)".
- 2026-11-15, answering "Merge main into this branch and rebuild the cleanup on 353cefdc?": "Merge and rebuild".
- 2026-11-15: "Commit the working tree changes with a sensible message." (commit on this branch; no push).

Not authorized: commit (until the 2026-11-15 instruction above), push, any change to the annotator
strategies (the plan defers them), gene UniProt xrefs (separate session), the lipid-unlump line
(`platform/tools/lipid/`, `model/curation/lipid_unlump/`, `model/candidates/lipid_unlump/`, and every reaction in
its scope as defined in §3), LP solves.

## 1. Baseline (checked 2026-11-13)

- Worktree `.claude/worktrees/peaceful-pascal-66d07e`, branch `claude/peaceful-pascal-66d07e`, HEAD
  `cf7b8dd94795d59fda81e00c4189736c3dbb3104`; uncommitted: the 2026-11-12 audit report and its STATE entry.
- Baseline model: `model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_bigg_annot_v2_mcnaugh_rxns_v3_r2121_close_r1172_close_20261112.xml`,
  SHA-256 `684f3752764fb8f25219008966953863141724f34a451dffb77c235aa3cb2f62`. It was the newest candidate on this
  branch's base `cf7b8dd9` (PR #51). Correction, 2026-11-15: `main` had already moved on to 353cefdc (PR #53,
  folinic-acid module on 684f3752). I had not compared with `main` when I wrote this (REPORT §5).
  The 113 audit targets have byte-identical annotations in e584151f, 2f352037 and 684f3752 (verified).
- Rebuild base after the 2026-11-15 merge: `main`'s
  `model/candidates/E5_…_r2121_close_r1172_close_fthf_20261112.xml`, SHA-256
  `353cefdc16ef1dc510eae49c168fc999254c2652c3994287af45589621a90cd9` (folinic-acid module on 684f3752, PR #53). The
  branch was fast-forwarded to `main` `791090382e47445fe658852d04978892ad98fc75`, and the uncommitted changes were
  reapplied on top. The 62 targets have identical annotations and notes in 684f3752 and 353cefdc (verified).
- Reference model: none designated (`model/candidates/README.md`).
- Source of the findings: [annotation_xref_audit_20261112](../annotation_xref_audit_20261112/REPORT.md).
- pytest baseline (root `.venv`, archive copy of `cf7b8dd9`): 477 passed, 3 failed, 18 errors, 24 skipped — the
  standing set recorded in `model/STATE.md`.

## 2. Targets and expected effect

62 reactions (list and per-key values: `model/curation/annotation/reaction_xref_cleanup_decisions_20261113.csv`,
374 rows; baseline state of every target: [before.json](before.json)):
- 53 whose MNXR describes another reaction, another specific substrate, or only part of the reaction;
- 9 whose MNXR is correct (a reviewed crosswalk rewrite) but whose other keys still hold the annotator's first,
  non-matching MNXR's values.

Expected effect: 896 values removed from 16 reaction annotation keys (`metanetx.reaction`, `mnx`, BiGG, KEGG, Rhea,
SEED, MetaCyc, SABIO-RK, VMH and their short duplicates); every touched key becomes absent. One note
`xref_cleanup` per reaction. Nothing else changes: no GPR, stoichiometry, bound, metabolite, gene, `ec-code`, SBO
term or crosswalk-decided value. FBA/FVA/knockout results cannot change: annotations and notes are not part of the
optimization problem (provable), so no solve is run.

## 3. Rule and exclusions

Per target and key: the audited wrong MNXR is removed from `metanetx.reaction` / `mnx` (wrong-MNXR rows only); any
other value is removed when every MNXR it maps to in MetaNetX 4.5 `reac_xref.tsv` fails the fingerprint test
(`disjoint`, `carrier_only`, `disjoint_single_compound`, `mnxr_covers_only_h_h2o`; recomputed on 684f3752 for every
key, short keys by their own values). Keys decided by the BiGG crosswalk for that reaction are never touched.

Excluded (51 of the 113 audit targets; listed with the reason in the manifest's `excluded_reactions`):
- 37 in lipid-unlump scope by that line's own files: named in `unlump_stage0_plan.csv`, the sn-core curation, the
  fatty-acyl-CoA or R1521 handoff/audit files, a lipid row of `fill_candidates.csv`, or using a lipid species of
  `missing_lipids.csv` / `fill_candidates.csv` (`(non-lipid)` rows ignored); plus targets sharing their audited MNXR
  with one of those;
- 14 fatty-acid-synthesis reactions (R2044, R2216, R2217, R2219, R2220, R2223, R2226, R2229, R2232, R2235, R2238,
  R2239, R2240, R2241), by the user's choice.

## 4. Implementation and acceptance

- Data: `reaction_xref_cleanup_decisions_20261113.csv` + `.json` manifest (SHA-pinned, `user_authorized`).
- Step: `platform/scripts/gem_annotate/reaction_xref_cleanup.py`, `apply_reaction_xref_cleanup`; flag
  `--reaction-xref-cleanup`; called last in `main.py`, after `--r1172-hmgcoa-transport-close`. Fails closed unless
  every row is in its `before` state (or every row already applied); refuses crosswalk-decided keys and excluded
  reactions.
- Focused tests on in-memory models (`platform/tests/test_reaction_xref_cleanup.py`).
- Also approved in the plan: correct the `_UBIQUITOUS_MNXM` comment in `annotate_reactions_extended.py` (comment
  only, no behaviour change).
- Acceptance:
  1. The control build (684f3752's command) reproduces 684f3752 byte for byte; the candidate is identical under two
     hash seeds.
  2. The whole-model diff baseline → candidate is exactly the 374 decision rows plus 62 `xref_cleanup` notes.
  3. pytest compared test by test with the baseline: only the new tests differ.
  4. Code review and an independent audit raise no unresolved blocking finding.
- Stop condition: any other difference, a stale row, or a failing guard → stop and report.
