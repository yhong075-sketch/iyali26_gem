# Review of reviewed GPRs overwritten by the metadata selection (2026-10-19)

## Authorization

Delegated by the session that built the 2026-10-19 GPR restore (`selection_gpr_restore_20261019`):

> Read-only investigation, then a plan for the user's approval; no model, data or code edits until
> approved. [...] Do not edit model XML, curation data or code, do not build, do not commit or push.
> No LP solves beyond static checks unless the user authorizes them.

Questions:
1. List every GPR the 2026-09-10 selection changes (9 rows), plus stoichiometry or bounds rows that
   interact with reviewed builder patches. Classify each: selection undoes a reviewed patch /
   selection is newer reviewed evidence / unclear.
2. For R765/R766/R671 and R190/R159, state the reviewed state with evidence tiers (AGENTS.md §5 for
   R190's AND vs OR).
3. Propose a pipeline fix that reuses the 2026-10-19 pattern; say whether to extend that step or add
   a separate one.

## Starting state

- This worktree: branch `claude/clever-kare-7d2e1a` at `285eaee`. It does not contain the
  2026-10-19 restore work.
- The 2026-10-19 restore work is uncommitted in the worktree
  `.claude/worktrees/orphan-reaction-newest-model-b8ddc8` (branch
  `claude/orphan-reaction-newest-model-b8ddc8`, HEAD `4d97b828…`). It was read there, not modified.
- `git diff 4d97b82 285eaee` is empty for `platform/scripts/gem_annotate`, `model/curation/selection`,
  `model/curation/essentiality`, `model/curation/gpr` and `model/source`. The builder files below
  are therefore the same in both worktrees.

| Input | SHA-256 |
| --- | --- |
| `model/curation/selection/metadata_reaction_selection.json` | `abfc4a94a03b4ec8da9e3e2fc3ba2e4f04806ef802b6e93c442274065ca0cb24` |
| Selection source `metadata.xml` (research root, `…/accepted/metadata.xml`) | `33c14f0c524197ec732e09ba86c4382d5b21bc25b4bd73b18eba5f452dd7e797` |
| Source model `model/source/iyali26.xml` | `5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d` |
| Candidate a56f3649 `E5_coq9_resp_pkg_r132_20261012.xml` | `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515` |
| Candidate c98961e9 `E5_coq9_resp_pkg_r132_gpr_restore_v2_20261019.xml` (other worktree) | `c98961e960c6a0f8283f485e32daf20b04dd70d998fad86fc09dca0d14ccae4c` |
| `platform/scripts/gem_annotate/patches.py` | `8d343c0e25261612580545a95442ef6884792dac699f780b4ab291674ad763f3` |
| `platform/scripts/gem_annotate/reaction_selection.py` | `e0a3b954fbaf454501746bddf950663fbd61db6c9139b6a73b2ace47dfb44d60` |
| `selection_gpr_restore.py` (other worktree, uncommitted) | `0a8c4dfd3a2f6fbc77c6260dea1031ff1e847825f8a6e40b2905b3c65dc66e1e` |
| `selection_cleared_gpr_restore.json` (other worktree, uncommitted) | `d048ff9046584c75b54cfdae438874e448c3c883e92d1ce994c1f4eec98aba95` |

None of these is a designated reference model (`model/candidates/README.md`).

## Scope and budget

- Static only: SBML loads with cobra 0.30.0 (main checkout `.venv`, Python 3.13.5), JSON
  comparison, element and charge sums. No LP, no build, no FVA.
- External identity look-ups: NCBI E-utilities, KEGG REST, UniProt REST.
- Outputs: this TASK.md, REPORT.md, `selection_rows.csv`, `scripts/`.
- Stop: any edit to model, data or code would be out of scope.
