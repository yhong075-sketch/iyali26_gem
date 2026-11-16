# Reaction cross-reference cleanup, candidate 7ccbb568 on main's 353cefdc: report, 2026-11-13/15

**Current candidate** (2026-11-15, §7): `model/candidates/E5_353cefdc_xref_20261115.xml`, SHA-256
`7ccbb56834a37ad4c75e434868e9e1b9e1f6a657a11ca17ee0a93dbe56232146`. It is the same change rebuilt on `main`'s
folinic-acid candidate 353cefdc after a fast-forward merge. §1–§6 describe the first round on 684f3752 (candidate
`9ce7c54a…`), which 7ccbb568 supersedes. The decision rows, the rule and every verified count apply unchanged to both.

Scope, authorization and acceptance criteria: [TASK.md](TASK.md). Baseline state of every target:
[before.json](before.json). Decisions: `model/curation/annotation/reaction_xref_cleanup_decisions_20261113.csv` (374
rows) and `.json` (manifest, notes, exclusions). Scripts, build logs, the per-ID fingerprint table on 684f3752 and the
whole-model diff: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/reaction_xref_cleanup_20261113/`. Current manifest:
`MANIFEST_20261115.sha256` (`2b8bcf1e…`, 90 files, all rounds). The first-round `MANIFEST.sha256` (`65e6eca1…`) is
kept as written; one script it lists, `whole_model_diff.py`, was later given an optional base argument.

**Evidence tiers.** *Verified*: computed this session, output saved. *Inferred*: reasoned. *Recalled*: from an
earlier session's file, not re-checked.

## 1. What changed

New candidate
`model/candidates/E5_…_mcnaugh_rxns_v3_r2121_close_r1172_close_xref_20261113.xml`, SHA-256
`9ce7c54ab0ed24ace1bf3c5d839582ed6715771041ab35596ca1be46a67a8588`, built on the 684f3752 flag set plus
`--reaction-xref-cleanup` (verified).

On 62 reactions, 896 cross-reference values were removed from 374 annotation keys, and each reaction got one
`xref_cleanup` note naming the removed MNXR and this curation (verified, whole-model diff). Every touched key is now
absent: no reaction kept part of a key.

| Key (long and short form) | Values removed |
| --- | ---: |
| `seed.reaction` / `seedr` | 122 / 122 |
| `bigg.reaction` / `biggr` | 70 / 142 |
| `vmhreaction` / `vmhr` | 78 / 78 |
| `rhea` / `rh` | 68 / 68 |
| `metanetx.reaction` / `mnx` | 53 / 7 |
| `kegg.reaction` / `keggr` | 18 / 18 |
| `metacyc.reaction` / `metacycr` | 16 / 16 |
| `sabiork.reaction` / `sabiorkr` | 10 / 10 |

`biggr` loses more than `bigg.reaction` because on 10 reactions the crosswalk rewrote `bigg.reaction` but not its
short duplicate, which still held the old MNXR's BiGG IDs (verified): R29, R30, R31, R33, R95, R636, R649, R1131,
R1908 and R2210. One example is the cGMP phosphodiesterase IDs GCPN/PDE4 on R29 (cAMP), R30 (cyclic dAMP), R31
(cIMP) and R33 (cCMP).

**The 62 reactions** (from the audit, labels inferred):
- 27 carry an MNXR for another specific substrate:
  - R218–R235, the 18 cis-prenyltransferase steps → the `PPi = IPP` rubber-allyltransferase stub;
  - R31/R33 → cGMP phosphodiesterase;
  - R144/R2252 → 2-phenylacetamide amidase;
  - R2213 → indole-3-acetonitrile nitrilase;
  - R178, R2094, R2122 and R141.
- 15 carry an MNXR for another reaction:
  - R1711 → HisA;
  - R706–R711, R1913 and R1914, the NDP ribonucleotide reductase steps → NTP reductase;
  - R182–R184;
  - R2032/R2033 (ABC transport vs uniport);
  - R288 (mannan synthase stub).
- 11 carry an MNXR that covers only part of the reaction:
  - R1999, R2019, R2020, R2027, R2028, R2037 and R2038, proton symports → H+ exchange;
  - R938, which moves UTP, UMP, FMN and H+ → UTP exchange;
  - R2097/R2098 → the `SAH = SAM` stub;
  - R2186 → the glycogen-synthase stub.
- 9 have a correct crosswalk MNXR, but their other keys still came from the annotator's first MNXR: R29, R30, R95,
  R577, R636, R649, R1131, R1908, R2210.

After the cleanup (verified):
- 53 of the 62 have no MNXR; the 9 stale ones keep their reviewed crosswalk MNXR.
- 10 keep a crosswalk-decided `bigg.reaction`: the 8 stale ones that have one, plus R31 and R33.
- 13 have no annotation left besides the SBO term.

**Not changed** (verified, whole-model diff):
- stoichiometry, bounds, GPRs, subsystems;
- every metabolite and gene;
- compartments, objective, groups;
- SBO terms and `ec-code`;
- every value the BiGG crosswalk decided.

## 2. Rule and why it is shaped this way

**Rule.** Per target and key, two kinds of value are removed:
- the audited wrong MNXR, from `metanetx.reaction` and `mnx` (wrong-MNXR rows only);
- any other value when every MNXR it maps to in MetaNetX 4.5 fails the fingerprint test against the reaction
  (`disjoint`, `carrier_only`, `disjoint_single_compound`, `mnxr_covers_only_h_h2o`).

The per-value classes were recomputed on 684f3752 for every key, short keys by their own values. The audit's
earlier run covered only KEGG, Rhea and BiGG (verified; `out/reaction_xref_fingerprint_all.csv`).

**Why value by value.** The approved plan said to remove the wrong MNXR's own bundle. Implementation showed that the
annotator also merged the xrefs of other "equivalent" MNXRs into the same wrong pick. On R1540, for example, these
are acetyl-CoA transport IDs. The user chose to remove them too, except in the lipid-unlump line.

The two sizing figures behind that question were computed over all 113 audit targets, before the exclusions
(verified; `scripts/size_removal.py` and `scripts/kept_values.py`, outputs in `out/size_removal_plan_rule.txt` and
`out/merged_equivalent_values.txt`):
- 1,464 values from the wrong MNXR's own bundle, counting `metanetx.reaction` and `mnx`;
- 130 merged-equivalent values on the 97 wrong-MNXR rows.

The independent audit recounted with its own definitions and got 1,432 + 104 and 128. These figures only shaped the
question to the user. The decision rows themselves are verified: 896 values on 62 reactions.

**Why a last step, not an annotator fix.**
- The BiGG crosswalk step refuses the build unless each of its targets holds its pinned starting value
  (`bigg_crosswalk_annotations.py`, verified by code). Those starting values are written by the annotator.
- Changing the annotator's picks would therefore break the crosswalk. That source fix stays deferred (audit REPORT
  §4).
- The two steps never touch the same (reaction, key) pairs: this step refuses every crosswalk-decided pair.
- This step runs after the crosswalk because its own guard requires the crosswalk's reviewed values. On the 11
  targets that share reactions with crosswalk decisions (R29, R30, R31, R33, R95, R577, R636, R649, R1131, R1908,
  R2210), each decided value must already hold its `after` state, since 9 of the notes say the crosswalk MNXR is
  kept. The guard was added after the code review (§3, item 7).
- It also runs last because its rows pin the 684f3752 state.

**Exclusions** (51 of the 113 audit targets, each with its reason in the manifest's `excluded_reactions`):
- **37 in lipid-unlump scope**, by that line's own files. Each is named in `unlump_stage0_plan.csv`, the sn-core
  curation, the fatty-acyl-CoA or R1521 handoff/audit files or a lipid row of `fill_candidates.csv`, or uses one of
  those files' lipid species (`(non-lipid)` rows ignored). A target that shares its MNXR with one of those is
  excluded too. That is how R1499/R1500 are excluded: they share MNXR118924 with R1501–R1510.
- **14 fatty-acid-synthesis reactions**, by the user's choice: R2044, R2216, R2217, R2219, R2220, R2223, R2226,
  R2229, R2232, R2235, R2238, R2239, R2240 and R2241.

The excluded rows stay as recommendations for the lipid owner, in the audit table.

## 3. Verification

0. **Build command** (from `platform/`, with `PYTHONHASHSEED=0`; `scripts/run_builds.py` runs it and records each
   argv in `builds/builds.json`). It is 684f3752's recorded argv plus one flag:

   ```bash
   python -m scripts.gem_annotate --research-root "$IYALI26_RESEARCH_ROOT" --offline --no-solve \
     --coq9-curation qcycle --energy-candidate E5 --coq9-biomass-dilution --coq9-respiratory-package \
     --r132-mito-adh --restore-selection-cleared-gprs --restore-selection-overwritten-gprs \
     --restore-selection-cpa-ura2 --remove-selection-tkl-fragment-gpr --remove-r153-peptidase-gpr --r2041-forward \
     --r2190-ribulose-product --r1964-gpr --gpi-chain-links --gpi-coproducts --tpp-requirement \
     --thiamine-phosphate-chemistry --thiamine-proton-lock --heme-iron-biomass --r2193-iron-gpr \
     --lipoylation-module --acyl-acp-chemistry --cofactor-biomass-terms --r2135-acyl-coa-pool --r1726-cpt-open \
     --dolichol-route --um7-nuclear-diffusion --r2171-pi-plc-1to1 --aars-thr-ser-split --r171-atp-synthase-gpr \
     --r673-prolyl-split --aconitase-aco2-split --trehalose-biomass-removal --r24-homocitrate-aconitase \
     --ceramide-pool-biomass --bigg-crosswalk-annotations --mcnaugh-reactions --r2121-malonyl-coa-close \
     --r1172-hmgcoa-transport-close --reaction-xref-cleanup \
     --output-model ../model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_bigg_annot_v2_mcnaugh_rxns_v3_r2121_close_r1172_close_xref_20261113.xml
   ```
1. **Control build** (verified). 684f3752's recorded argv (47 arguments after the program path), run from this
   worktree without the new flag, reproduces `684f3752764fb8f2…` byte for byte. The builder code outside the new
   step therefore behaves as before.
2. **Seeds** (verified). The candidate is byte-identical under `PYTHONHASHSEED` 0 and 1.
3. **Whole-model diff** (verified; `out/whole_model_diff.json`).
   - Every reaction, metabolite and gene was compared on name, stoichiometry, bounds, GPR, subsystem, formula,
     charge, compartment, each annotation key and each note, plus compartments, objective and groups.
   - The only differences are the 374 decision keys (all removed) and the 62 notes. Nothing unexpected.
   - In the XML, 13 reactions lose their whole RDF annotation block, because the removed xrefs were all they had.
     SBO is an attribute.
4. **No LP.** Annotations and notes are not part of the optimization problem, so WT growth, FVA and every knockout
   call are unchanged (provable from the model structure; the diff shows no stoichiometry, bound or GPR change).
5. **Focused tests** (verified). `tests/test_reaction_xref_cleanup.py`, 9 tests:
   - scope;
   - exact removals and notes;
   - idempotence;
   - stale, partly applied and wrong-note models fail closed without edits;
   - a target whose crosswalk value is not applied is refused (added after the review);
   - scope/row guards, including crosswalk-decided keys and excluded reactions;
   - context restore;
   - SBML roundtrip;
   - opt-in flag.

   All pass, with the 7 crosswalk tests (16 passed).
6. **Full pytest** (verified; root `.venv`, short-path archive copies; `--continue-on-collection-errors`):
   - Before (`cf7b8dd9`): 477 passed, 3 failed, 18 errors, 24 skipped.
   - After, on the pre-review code: 485 passed, 3 failed, 18 errors, 24 skipped.
   - Compared test by test (junit), the only difference is the 8 new tests, all passing. The 3 failed and 18 errors
     are the standing set.
   - After, on the post-review code: 486 passed, 3 failed, 18 errors, 24 skipped. Against the baseline, the only
     difference is the 9 new tests, all passing.
   - Two text-only edits came after that run: the crosswalk comment in `main.py` and the manifest's authorization
     wording. Neither is read by a test. The focused tests were rerun on the final files: 16 passed.
   - Logs and junit files: `pytest/` in the task artifacts.
7. **Code review** (fresh read-only reviewer, following `.claude/agents/code-reviewer.md`). No blocking finding.
   Acted on:
   - **Fail-open precondition.** Nothing stopped a build with `--reaction-xref-cleanup` but without the crosswalk,
     and the 9 stale-row notes would then be false. Fixed: the step now refuses unless the crosswalk's values on
     its targets hold their reviewed `after` state, and a test covers it.
   - **Wrong stated reason for running last.** The comments and §2 said the crosswalk guard pins the values this
     step removes. It does not; corrected.
   - **An unused partial-key branch.** Removed; every row must now remove its whole key.
   - **The annotator comment fix was not reported** (§4).

   After the fixes, the candidate was rebuilt: same SHA `9ce7c54a…` for seeds 0 and 1, the control still reproduces
   684f3752, and the `.build.json` source SHAs match the final code. The pre-review build files are kept in
   `builds_v1_before_review_fix/` in the task artifacts.
8. **Independent audit** (fresh read-only agent given only TASK.md, the diff, the logs and the deliverables).
   - **Supported:** every model-change claim, independently re-diffed (374 keys, 896 values, 62 notes, nothing
     else; all 51 excluded reactions byte-unchanged), plus the rule, the exclusions, the provenance and the tests.
   - **Corrected:**
     - the §5 enzyme-class statement (wrong for five reactions);
     - wording on R29–R33, R938 and the fatty-acid range;
     - the argv count;
     - the missing derivation of the sizing figures (§2);
     - the manifest's authorization wording, now quoting TASK.md §0;
     - the crosswalk comment in `main.py`, which my step made untrue (it now names the exception).
   - **Left as notes:** two below and one in §5 (open items).
   - After these fixes, all three builds were rerun: the same `9ce7c54a…` for seeds 0 and 1, the control still
     `684f3752…`, and every `.build.json` source and data SHA matches the working tree. The earlier build files are in
     `builds_v1_before_review_fix/` and `builds_v2_before_audit_fix/`.

## 4. Problems met on the way

- **File-name limit.** The first candidate name, `…_xref_cleanup_20261113.xml` (247 characters), made the
  deterministic writer's temporary file (`.<name>.<8 chars>.xml`) exceed macOS's 255-byte name limit
  (`OSError 63`; verified). The candidate is `…_xref_20261113.xml` (239 characters). Any later step that appends
  more than 2 characters to this name will hit the same limit. The naming scheme needs a shorter form before the
  next candidate (inferred).
- **XML escaping in notes.** An apostrophe in the first note text came back as `&apos;` after an SBML roundtrip,
  so a re-read model would look stale to the guard (verified by test). The notes contain no `<>&'"` now; the
  generator and a test enforce it.
- **One comment fixed in the annotator**, as planned. `annotate_reactions_extended.py:27` now names MNXM5 as
  NADP(+) and MNXM9 as phosphate (verified against `chem_prop.tsv`). Behaviour is unchanged: the control build with
  it reproduces 684f3752.
- **Lipid-unlump file mentions.** The reviewer found 6 cleaned reactions in the colleague's `fill_candidates.csv`,
  only in `is_lipid=no` rows: R2020, R2028, R2033, R2038, R2098 and R2213 (confirmed by the audit). That follows the
  rule in §2. The cleaned sterol (R178, R182–R184) and dolichol/prenyl reactions (R141, R218–R235, R288) are not
  named in that line's files either. The owner may want to know about both groups.
- **libSBML message on R304.** All builds log a libSBML error setting notes on R304. The control build has it too,
  so it predates this change. R304 is not a target.

## 5. Found, not changed

- **R288's `ec-code` is 2.4.1.32.** Verified on KEGG: [EC 2.4.1.32](https://www.kegg.jp/entry/2.4.1.32) is
  glucomannan 4-β-mannosyltransferase. R288's own equation (GDP-Man + dolichyl phosphate → GDP + Dol-P-Man) is
  [EC 2.4.1.83](https://www.kegg.jp/entry/2.4.1.83). R288 has no GPR genes carrying 2.4.1.32, so it came from the
  wrong MNXR. `ec-code` was out of scope.
- **Eight more reactions keep an `ec-code` not taken from their genes** (verified): R29, R30, R31, R33, R1908, R2210,
  R2213 and R2252. My first draft said these all name the right enzyme class. The independent audit showed that
  is wrong for five of them; I re-checked on KEGG REST (verified):
  - R29, R30, R31 and R33 carry [3.1.4.17](https://www.kegg.jp/entry/3.1.4.17) (general cyclic-nucleotide
    phosphodiesterase, fits) and [3.1.4.35](https://www.kegg.jp/entry/3.1.4.35) (cGMP-specific). The latter is
    wrong for cAMP, cyclic dAMP, cIMP and cCMP.
  - R2252 carries [3.5.1.4](https://www.kegg.jp/entry/3.5.1.4) (amidase, fits) and
    [3.5.1.11](https://www.kegg.jp/entry/3.5.1.11) (penicillin amidase), which is wrong for indole-3-acetamide.
  - In each case the set equals the wrong MNXR's EC set, so it likely came from there (inferred).
  - R1908 (2.7.1.33, pantothenate kinase) and R2210/R2213 (3.5.5.1, nitrilase) fit their reactions (inferred from
    the EC names).

  Together with R288, six reactions keep at least one wrong EC. `ec-code` was out of scope; this is a follow-up.
- **The annotator still makes the same picks.** Any build without `--reaction-xref-cleanup` reproduces the wrong
  xrefs. The source fix is in the audit REPORT §4.
- **The 51 excluded reactions keep their wrong xrefs**, for the lipid owner to review. The 198 partial rows of the
  audit were not touched.

- **Main has moved past this base.** PR #53 (folinic-acid module: R748 to SHMT, R990 closed, cytosolic R_THFAT;
  candidate `353cefdc…`, also built on 684f3752) merged into `main` after this branch was cut at `cf7b8dd9`.
  Checked against `main:model/candidates/…_r1172_close_fthf_20261112.xml` (SHA `353cefdc16ef1dc5…`):
  - the 62 targets have identical annotations and notes there;
  - none of the three folinic-acid reactions is a target (verified).

  So the same decision rows should apply on 353cefdc, and the crosswalk values the guard requires also hold there
  (inferred; the audit checked them). Combining the two needs a merge of `main` into this branch and a rebuild with
  both flags. The merge will likely need manual resolution:
  - both branches insert next to the `r1172_hmgcoa_transport_close` lines in `main.py` (signature, `build_model`
    keywords) and in `cli.py` (argparse), though the two calls sit in different places;
  - PR #53 also adds a comment calling the crosswalk "the last step that writes annotations", which this step
    contradicts.

  Neither the merge nor the rebuild is authorized yet.

- **Links that need the audit report committed.** The decision CSV's `evidence` column and this report point to
  `model/reports/annotation_xref_audit_20261112/`, which was uncommitted. It is committed together with this change
  (§7).

## 6. Limits

- The adjudication labels are my reading of two equations per row (inferred). The fingerprint classes they rest on
  are verified.
- MetaNetX 4.5 local tables were the only bridge from KEGG/Rhea/BiGG/SEED/MetaCyc/SABIO-RK/VMH IDs to equations.
  Those databases were not queried directly.
- The lipid-unlump boundary is drawn from that line's files. A reaction those files do not name could still matter
  to its owner.

## 7. Rebuild on main's 353cefdc (2026-11-15)

**Authorization.** The user, answering "Merge main into this branch and rebuild the cleanup on 353cefdc?":
"Merge and rebuild".

**Merge** (verified).
- The branch base `cf7b8dd9` is an ancestor of `main`, so the merge was a fast-forward to `79109038` (PR #53). No
  merge commit, no stash.
- The four modified tracked files were saved first, as a patch against `cf7b8dd9` plus copies, in
  `merge_backup/` in the task artifacts.
- They were then reset, the branch fast-forwarded, and the patch reapplied 3-way.
- Conflicts and their resolution:
  - `main.py` (signature and `build_model` keywords) and `cli.py` (argparse): both sides kept, the folinic-acid
    line first.
  - `STATE.md`: both sides inserted at the top; both kept, newest first.
- Two comments in `main.py`, on the folinic-acid step (from PR #53) and the mcnaugh step, call the crosswalk "the
  last step that writes annotations". They now add "before --reaction-xref-cleanup".
- No commit during the merge. The work was committed afterwards, on the user's instruction (below).

**Base.** `main`'s `E5_…_r2121_close_r1172_close_fthf_20261112.xml`, SHA-256 `353cefdc16ef1dc5…`. It carries the
folinic-acid module on 684f3752 (R748 to SHMT, R990 closed, cytosolic R_THFAT added). The 62 targets have identical
annotations and notes there, and none of the three folinic-acid reactions is a target (verified). The decision rows
apply unchanged.

**Build** (verified; `scripts/run_builds_353cefdc.py`, outputs in `builds_353cefdc/`).
- Command: 353cefdc's recorded argv (48 arguments after the program path; it adds `--folinic-acid-module` to the
  list in §3), plus `--reaction-xref-cleanup`, from `platform/`.
- Control (no new flag) reproduces 353cefdc byte for byte.
- Candidate `model/candidates/E5_353cefdc_xref_20261115.xml`, SHA-256
  `7ccbb56834a37ad4c75e434868e9e1b9e1f6a657a11ca17ee0a93dbe56232146`, identical for `PYTHONHASHSEED` 0 and 1.
- `.build.json`: head `79109038`; every source and data SHA matches the working tree; both the folinic-acid step
  and the cleanup report `applied`.

**Name.** The cumulative scheme would give `…_fthf_xref_20261115.xml` (244 characters). That is over the 241 that
the writer's temporary file allows (§4). The candidate therefore gets a short name: energy candidate, base SHA,
step, date. The full flag lineage is in its `.build.json`. Later candidates on this line need the same kind of short
name.

**Whole-model diff** 353cefdc → 7ccbb568 (verified; `out/whole_model_diff_353cefdc.json`). The only differences are
the 374 decision keys (all removed, 896 values) and the 62 notes. Nothing unexpected. No LP: annotations and notes
are not in the optimization problem, so 353cefdc's recorded WT and knockout results carry over (provable).

**Tests.**
- Focused, on the merged tree (verified): `test_reaction_xref_cleanup.py`, `test_bigg_crosswalk_annotations.py`
  and `test_folinic_acid_module.py`, 23 passed.
- Full pytest (verified; root `.venv`, short-path archive copies, `--continue-on-collection-errors`; logs and junit
  in `pytest_353cefdc/`):
  - `main` `79109038`: 484 passed, 3 failed, 18 errors, 24 skipped.
  - `main` plus the changes: 493 passed, 3 failed, 18 errors, 24 skipped.
  - Compared test by test, the only difference is the 9 new tests, all passing.
  - That run's copy predates the review's comment and help-text edits. The focused tests were rerun on the final
    tree: 23 passed.

**Review of the merge and rebuild** (fresh read-only reviewer). No correctness bug and no AGENTS.md violation. It
verified independently:
- the merge lost nothing of `main`; `git diff main` is exactly the intended set;
- both new steps are wired end to end, in the order folinic acid → … → crosswalk → R2121 → R1172 → cleanup;
- the provenance, byte for byte;
- the model change: an XML-level and a cobra diff, the 19 crosswalk `after` values on the 11 shared reactions, and
  the folinic-acid targets untouched.

Acted on:
- **STATE entry stale.** Rewritten for 7ccbb568.
- **Task manifest stale.** A new `MANIFEST_20261115.sha256`; the old one is kept as written.
- **A blank line before the next STATE heading lost in the merge.** Restored.
- **"Runs last" in the R2121 and R1172 comments (`main.py`) and help text (`cli.py`).** These were no longer true
  once this step followed them; now "runs after the crosswalk". This step's own help now names the 684f3752 or
  353cefdc flag set.
- **`model/candidates/README.md`** now lists 7ccbb568 and the superseded 9ce7c54a.

Not acted on: the test fixture stays 684f3752. The rows pin target states that are byte-identical in 353cefdc, and
the rebuild plus the diff cover the 353cefdc case.

Those edits change only comments, help text and documentation. The three builds were rerun on them:
- the candidate is again `7ccbb568…` for seeds 0 and 1;
- the control is `353cefdc…`;
- every source and data SHA in `.build.json` matches the working tree.

The first 353cefdc builds are kept in `builds_353cefdc_v1_before_review_fix/`.

**Superseded.** 9ce7c54a (on 684f3752) stays in `model/candidates/`, untracked, as the first-round result; 7ccbb568
replaces it.

**Commit** (2026-11-15). The user: "Commit the working tree changes with a sensible message." Everything of
this task is committed on `claude/peaceful-pascal-66d07e`, the audit report of 2026-11-12 included. The superseded
first-round candidate 9ce7c54a stays untracked, like other superseded candidates. No push.
