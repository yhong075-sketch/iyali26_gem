# R153 made gene-less (R153-GPR-UNASSIGNED-20261019), 2026-10-19

## 0. Summary

- **New candidate:** `model/candidates/E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_r153_gpr_v2_20261019.xml`,
  SHA-256 `d941d2cdd5a6628ab16466f5880aab1c2e1edb3a6d6d8bb1b03b0f9e3ef2d557`.
  - Built from `dc471d7` plus this change, with the baseline's flags plus `--remove-r153-peptidase-gpr`, and
    `PYTHONHASHSEED=0`.
  - Not a reference model.
- **What changed in the model:**
  - R153 (argininosuccinate synthase) lost its GPR `YALI1D17462g`.
  - R153's notes changed: 6 replaced, 8 added.
  - The gene stays in the model with no reactions.
  - Nothing else differs from the baseline `51c6fabd…`: no other field in any reaction, metabolite, gene,
    group, objective or model-level entry.
- **Why:** the gene encodes an A1 aspartic peptidase, and W29 has no annotated argininosuccinate synthase
  gene. See the evidence review [../r153_ass_gpr_20261019/REPORT.md](../r153_ass_gpr_20261019/REPORT.md).
- **Not done:** no LP was run, nothing was committed, and nothing was published.
- **Superseded:** the first build `…_r153_gpr_20261019.xml` (`3c9e734f…`) is superseded by v2 (§5). It is kept
  until the user decides.

## 1. Authorization

- The user, 2026-10-19:
  1. "Approve option A for R153: a gene-less GPR through a new guarded post-merge step and opt-in flag, with
     YALI1D17462g kept unlinked."
  2. Then "Approve the R153 implementation plan as written, including the fast-forward … to dc471d7. Builds
     and pytest as budgeted; no LPs; no commit."
  3. Then, after the control finding, "Proceed with the R153 code plan … Amend criterion 1 … Rebuild the
     control and build the candidate both with PYTHONHASHSEED=0".
- The detailed plan was first displayed after approval 2. Code edits started only after approval 3.
- Full quotes and the budget are in [TASK.md](TASK.md).

## 2. What changed and how

| File | Change |
| --- | --- |
| `model/curation/gpr/r153_gpr_revision.json` (new) | The decision as data: R153's before state (GPR, bounds, 6-species stoichiometry, species identity, the merge note), after GPR `""`, the gene's identity and R153-only membership, 6 replaced notes with old and new values, 8 added notes, and 10 evidence rows with URL, tier and date |
| `platform/scripts/gem_annotate/r153_gpr_revision.py` (new) | `apply_r153_gpr_revision(model)` |
| `platform/scripts/gem_annotate/cli.py` | opt-in `--remove-r153-peptidase-gpr` (default off) |
| `platform/scripts/gem_annotate/main.py` | parameter `r153_gpr_revision`; call after the three selection restores, which is after `merge_r153_r2176`; argument mapping |
| `platform/tests/test_r153_gpr_revision.py` (new) | in-memory tests, no solver |

**Code logic.** `apply_r153_gpr_revision` checks everything before it edits anything:
- the spec's scope;
- that R153 exists and R2176 is gone;
- R153's stoichiometry, bounds and species formula/charge/compartment against the record;
- the gene's NCBI/RefSeq/UniProt annotation;
- the merge note.

Then it branches:
- **Already in the after state** (GPR `""`, gene unlinked, all new notes present): it returns
  `already_correct` and changes nothing.
- **Otherwise it requires** GPR = `YALI1D17462g`, R153 as the gene's only reaction, every note to be
  replaced still at its recorded old value, and no new note already present. Any mismatch raises
  `ValueError` with the model untouched.
- **Only then** does it clear the GPR and merge the notes. The gene is not removed.

The 2026-09-15 files (`r153_gpr_assignment.json`, `r153_r2176_merge.json`) and `patches.py` are unchanged.

**Notes.**
- Replaced: `curated_gpr_assignment`, `gpr_evidence_status`, `gpr_evidence_basis`, `gpr_evidence_label`,
  `gpr_experimental_confirmation`, `r153_model_interpretation`.
- Added: `r153_gpr_revision`, `…_authorization`, `…_basis`, `r153_w29_locus`, `r153_identity_candidate`,
  `r153_screen_coverage`, `r153_unlinked_gene`, `r153_gpr_revision_supersedes`.
- The other 18 keys stay as history. `r153_gpr_revision_supersedes` names the history keys that describe the
  2026-09-15 state. Two of them keep their original Chinese text: the replaced label's old value is preserved
  in the JSON, and `r153_authorization` is preserved in the model. All new text is English.

## 3. Evidence (from the review; tiers as recorded there)

- **Verified online 2026-10-19:**
  - `YALI1D17462g` = NCBI Gene 2910235 = XP_502815.3, 384 aa.
  - It carries only A1 aspartic-peptidase signatures (PF00026; PS00141 at 76–87 and 259–270; cd05474) and no
    ASS signature.
  - KEGG's ASS yli:90949668 / YALI2_D00589g is a DSM 3286 gene. Its NCBI Gene record is discontinued and
    NC_089381.1 is suppressed.
  - NC_090773.1 / CP017556.1 3,604,630–3,609,948 is unannotated.
  - The step is not in any of the three screens.
- **Recorded (2026-09-11):** frame-split ASS homology at the W29 locus. An audit alignment on 2026-10-19
  supports it.
- **Inferred:** that `YALI1D17462g` has no ASS activity in vivo, and that W29 has a functional ASS. Neither is
  tested.

## 4. Verification (all run this session)

| Check | Result | Level |
| --- | --- | --- |
| Unseeded control (no new flag) vs `51c6fabd` | SHA `8f083627…`. `model_diff.py`: **0 differences**. Bytes differ only in R159's reactant order ([logs/control_check.txt](logs/control_check.txt)) | verified |
| Cause of that byte difference | `selection_gpr_restore.py:154` iterates a hash-randomized `set` (code reading). Confirmed by the side session: 8 seeds gave 8 R159 orders | verified (side session) |
| Seeded control (`PYTHONHASHSEED=0`) vs `51c6fabd` (amended criterion 1) | SHA `aaa3b3d1…`. 0 field differences ([control_seed0_vs_51c6fabd_diff.json](control_seed0_vs_51c6fabd_diff.json)). Bytes: only R159, same line multiset; all non-reaction text identical | verified |
| v2 vs `51c6fabd`, whole model ([model_diff.json](model_diff.json)) | exactly 3 differences: R153 `gpr` `YALI1D17462g`→`""`; R153 `notes` (8 added, 6 replaced, 0 removed, all values equal the JSON); `YALI1D17462g.reactions` `[R153]`→`[]` | verified |
| v2 vs seeded control, bytes ([byte_and_static_checks.json](byte_and_static_checks.json)) | only the `R_R153` block differs. All other text is identical, including the gene-product list | verified |
| Static | R153 element and charge balance 0 (unchanged); bounds [-1000, 1000]; 1,073 genes; genes without reactions 22 → 23 (+`YALI1D17462g`) | verified |
| Build side outputs | `coq9_genes.tsv` byte-identical to the baseline's. Input data SHA-256s are identical apart from the new JSON. Changed sources: `cli.py`, `main.py`, the new module | verified |
| Focused test | 3 passed, 11 guard subtests | verified |
| Full pytest, before → after v2 ([logs/pytest_comparison_v2.txt](logs/pytest_comparison_v2.txt)) | 181 → 184 passed; 11 failed, 14 errors and 31 skipped in both. Per test, identical except the 3 new tests, which pass | verified |
| v1 vs v2 ([model_diff_v1_vs_v2.json](model_diff_v1_vs_v2.json)) | only the two refined note texts (`r153_identity_candidate`, `r153_gpr_revision_supersedes`) | verified |

**Provable without solving:**
- The WT LP of v2 equals that of `51c6fabd`, because stoichiometry, bounds, objective and compartments are
  identical (`model_diff.py` shows none of them differs), and a GPR is not part of the LP.
- A KO of `YALI1D17462g` on v2 changes no bound, because no reaction references the gene. It is therefore
  the WT problem.

No LP was run to observe either statement.

## 5. Independent audit and what was done

A fresh read-only auditor re-derived every claim with its own code. Result: claims 1–8 and 10 **supported**,
claim 9 **partially supported**, none unsupported. The partial claim is "all note text English": two
historical Chinese values remain, as described in §2.

| Audit point | Action |
| --- | --- |
| JSON evidence rows: KEGG yli:2910235 cited only through T01033; Patterson S3 not cited; `r153_identity_candidate` without an evidence level | **Fixed** in the JSON (10 evidence rows; note text now says "No verified native name; computational evidence only"). Rebuilt as v2 |
| Kept notes `r153_previous_gpr: YALIUNK2` and others describe the 2026-09-15 state | **Fixed by wording:** `r153_gpr_revision_supersedes` now names those keys as history and states the immediate previous GPR. The keys themselves are unchanged (the approved plan replaces 6 keys only) |
| TASK.md said "8 failures" | **Fixed** (8 in the selection-fix environment, 11 here) |
| The hash seed is not recorded in `.build.json` | **Reported.** `PYTHONHASHSEED=0` is in [commands.txt](commands.txt). v2's bytes depend on it (R159 order); its semantics do not |
| `merge_r153_r2176` would raise if re-run on a revised R153 | **Reported.** Builds always start from `model/source/iyali26.xml` and merge first, so this path is not reached. There is no test for it |
| A `RuntimeError` after the edit would leave the model mutated | **Reported.** It cannot trigger after the checks pass; the build would abort before writing |
| Test gaps (partly applied state, missing R153, `main.py` order) | **Reported.** The build log line "R153 GPR revision: applied" is the empirical wiring check |
| The full pytest includes tests that run their own small solves | **Reported.** These are the authorized pytest runs. No LP was run on the baseline or the candidates |

The first rebuild attempt targeted the v1 path. The builder refused it before building ("a canonical copy
requires a new, separate output path"), so v2 was written to a new file. That attempt's log was overwritten by
the v2 build log; the refusal is recorded in [commands.txt](commands.txt).

## 6. Effect on scoring (reported separately; not a reason for the change)

- `YALI1D17462g` is a model gene not on the positive list, so it is in the user-defined negative class. On v2
  its KO touches no reaction. It therefore scores TN at every cutoff (KO/WT = 1, provable).
- Before the change its 15% call depended on the arginine-uptake-limited KO ratio: 0.110 on the 2026-09-15
  model and 0.159 on the 2026-10-12 model. It was not computed on `51c6fabd`.
- The true ASS has no screen row, so no positive-class count changes.

## 7. Limits and next gates

- **What the model cannot do now:** it makes no gene-level prediction for argininosuccinate synthase. This
  reflects the missing W29 gene model; R153 still needs an enzyme biologically.
- **Not tested:**
  - the W29 ASS coding sequence and its expression;
  - any LP on v2;
  - FVA;
  - other media.
- **Optional LPs, separate authorization:**
  1. KO `YALI1D17462g` on `51c6fabd` (the before value);
  2. WT on v2 (provably equal).
- **Merge coordination:**
  - The D02625g line (`claude/cool-nightingale-a2f471`, already merged with `dc471d7`) and this change both add
    flags to `cli.py` / `main.py`.
  - The R159 determinism fix (side session `claude/great-wu-64e590`) touches `selection_gpr_restore.py`. Once
    it lands, builds no longer depend on the seed.
- **Commit or publish:** needs the user's explicit instruction.
