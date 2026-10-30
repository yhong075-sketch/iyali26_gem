# Evidence pipeline, M0: survey (2026-10-29)

## Authorization

Brief from the parent Claude session (`local_6979ed9d-1b85-4410-a97a-6ff65ec008f1`), relaying the user's
request: follow the specification `/Users/david/Downloads/EVIDENCE_PIPELINE_PROMPT.md` ("iYali26 evidence
pipeline: claims, contradictions, layered chain of evidence, reviewer integration") milestone by milestone,
and **start with M0 only**. M0 is read-only: survey the repo, write `docs/reviewer_interface.md` and this
report, then stop and ask the spec's §10 questions. No pipeline code, no model, GPR, medium, benchmark or
curated-data change, no computation, no commit or push.

## Starting state

| Item | Value |
| --- | --- |
| Worktree / branch | `.claude/worktrees/evidence-pipeline-ec8347`, `claude/evidence-pipeline-ec8347` @ `5fb75e7` |
| `origin/main` | `f212874` (PR #34 merge); this branch is 63 commits behind and 0 ahead |
| Sync | The app's sync was refused: the merge touches the protected `.claude/agents/code-reviewer.md` and the repository origin is not yet confirmed in the app. **The survey therefore read `origin/main` with `git show`**, not the stale working tree. |
| Model executed | None (no computation in M0) |
| Reference model | None designated (`model/candidates/README.md`) |
| Newest candidate on `main` | `d8d1bac78138f59ba6beefcbddae7f93d0a20e245fea6f5dd9c79760ee535751` (UM7 + R2171 + Thr/Ser split, 2026-10-29, `ff7fc62`); read as XML text only |
| Build input | `model/source/iyali26.xml` `5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d` (unchanged between this branch and `main`) |
| Research workspace | `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research` (read only) |
| Spec | `/Users/david/Downloads/EVIDENCE_PIPELINE_PROMPT.md`, 422 lines, read in full; not versioned in the repo |
| Sibling notes read | `essentiality_fn_causes_handover.md`, `logic_chain_handover.md` (shared notes folder) |
