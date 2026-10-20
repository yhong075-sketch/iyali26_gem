# TASK: R153 argininosuccinate synthase GPR identity review (2026-10-19)

## Authorization and scope

- Delegated by session `local_6fad52b5…` (R190 CPA review line). That review flagged the R153 finding in its
  REPORT §5.3 (worktree `claude/fervent-knuth-791bc8`) as a candidate, not checked online.
- **Read-only evidence review, then a plan for the user's approval.**
- Not in scope: model, curation data or code edits; builds; LPs; commits; STATE.md edits.
- Questions:
  1. What GeneID 2910235 / XP_502815.3 is, and what yli:90949668 / YALI2_D00589g is, verified online (NCBI,
     KEGG, UniProt), with length and domains.
  2. Whether a YALI1 locus matches YALI2_D00589g; whether that gene is in the model or in the
     Cas9/Cas12a/Tn-seq gene sets.
  3. How YALI1D17462g got onto R153: git history, curation data, builder, source model; the YALI0/YALI1
     locus-number collision.
  4. A static list of every reaction using YALI1D17462g, and of every reaction that could replace R153.
  5. A GPR correction proposal (if any) with AGENTS.md §4–5 evidence tiers. No GPR change to improve a
     match rate.

## Starting state

| Item | Value |
| --- | --- |
| This worktree | `claude/eloquent-borg-1bd1b0`, HEAD `285eaeee0e403f1c971267898eaed1a884edb727`, clean |
| Model examined | `model/candidates/E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_20261019.xml`, SHA-256 `51c6fabd486f4e1a9d6365bf7740bde3eca20fa05f1d5454a73586d3276db3fa` |
| Where it lives | worktree `claude/clever-kare-7d2e1a`; uncommitted when this task began. During the task it was committed as `dc471d75480bd688280df0bed5af2db8bbd7cebb` (blob SHA re-checked: identical) |
| Build of that model | `--offline --no-solve --coq9-curation qcycle --energy-candidate E5 --coq9-biomass-dilution --coq9-respiratory-package --r132-mito-adh --restore-selection-cleared-gprs --restore-selection-overwritten-gprs --restore-selection-cpa-ura2` (from its `.build.json`) |
| Source model | `model/source/iyali26.xml`, SHA-256 `5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d` |
| R153 curation data | `model/curation/gpr/r153_gpr_assignment.json` (`43db4c87…c2bed4`), `model/curation/gpr/r153_r2176_merge.json` (`1f78c8af…8c4d`) |
| Builder | `platform/scripts/gem_annotate/patches.py` (`8d343c0e…63f3`): `apply_r153_gpr_assignment`, then `merge_r153_r2176` |
| Reference model | none designated (`model/candidates/README.md`) |
| Research root | `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`, given in the task brief. `IYALI26_RESEARCH_ROOT` was unset in the shell; no build was run, so nothing fails closed |

## Prior work reused (read, not re-run)

- `artifacts/tasks/r153_r2176_gpr_review_20260911/REPORT.md` and
  `artifacts/tasks/argininosuccinate_identity_20260911/REPORT.md`. These cover BLAST, TBLASTN on the W29
  locus, US-align and AlphaFold reuse.
- STATE.md entries from 2026-09-11 and 2026-09-15 on R153: the OR assignment, the single assignment and the
  R2176 merge.

## Outputs

- `REPORT.md` (this folder), `commands.txt`, `scripts/static_r153.py`.
- Raw fetched records and static listings:
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r153_ass_gpr_20261019/` (with `SHA256SUMS`).

## Stop condition

Stop when the five questions are answered at stated evidence tiers and a plan is presented. Any model change
waits for the user's explicit approval.
