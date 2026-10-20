# R2041 (D-arabinokinase) direction review (2026-10-19)

## Authorization

Delegated by the session that ran the selection GPR review (`selection_gpr_review_20261019`, worktree
`claude/clever-kare-7d2e1a`):

> Read-only review, then a plan for the user's approval; no model, data or code edits, no build, no
> commit. [...] do not run LPs without the user's authorization.

Questions:
1. Statically trace the loop that `fix_d_arabinokinase_direction_and_proton` describes in a56f3649:
   producers and consumers of D-arabinose and D-arabinose 5-phosphate, bounds, GPRs.
2. Gather the direction evidence (enzymology, thermodynamics) per AGENTS.md §5. Never use a reaction
   name or another model's bounds.
3. Propose, for approval, whether and how to restore forward-only bounds through the pipeline. Include
   the LP checks that would be needed: a closed-exchange ATP test with and without R2041 reverse flux.

## Starting state

- Worktree `claude/bold-chaplygin-1c16df` at `285eaee` (clean at start).
- The main checkout is on `restructure/model-platform` at `2a34475`. It was not used.
- Model under review: a56f3649 (below). It is not a designated reference model
  (`model/candidates/README.md`).
- Reference model: none designated.
- Candidates read only for comparison: c98961e9 (worktree `orphan-reaction-newest-model-b8ddc8`,
  uncommitted) and 51c6fabd (worktree `clever-kare-7d2e1a`, uncommitted).

| Input | SHA-256 |
| --- | --- |
| a56f3649 `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml` | `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515` |
| Its baseline `model/candidates/E5_coq9_respiratory_package_20261008.xml` | `0ac2addfa0cb6de2d89393e433766fb04325ffe42c2d131eaefaca3bcec29dba` |
| Source `model/source/iyali26.xml` | `5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d` |
| `model/curation/selection/metadata_reaction_selection.json` | `abfc4a94a03b4ec8da9e3e2fc3ba2e4f04806ef802b6e93c442274065ca0cb24` |
| `platform/scripts/gem_annotate/patches.py` | `8d343c0e25261612580545a95442ef6884792dac699f780b4ab291674ad763f3` |
| `platform/scripts/gem_annotate/reaction_selection.py` | `e0a3b954fbaf454501746bddf950663fbd61db6c9139b6a73b2ace47dfb44d60` |
| `platform/scripts/gem_annotate/main.py` | `d1917bd285bd676a3380972214ebd3a65619f4d11a976db700db24434afd99d6` |
| `model/curation/direction/r1931_direction.json` (precedent) | `e303f8fdebb9dd9bc184ca1a482bdd385c77dff4d143a6d9f9adc792736fa667` |
| `model/curation/direction/r1159_direction.json` (precedent) | `4fe1b35b52e8858a74a4debe66054f4206fa35be39ca80db39f234a7fd819901` |
| Recorded closed-ATP controls `model/reports/coq9_respiratory_package_20261008/controls.json` | `78b9f3c4418cc84ed3d60e8f4e0edb184a8eb6d4366b97e772e84e59f33fce66` |
| Their runner `.../run_controls.py` | `b012744b49dc38e44853fbbf598b1019e87cbbe3ad2a6db9a1ab525b29c85c07` |
| `platform/tools/diagnose_closed_energy.py` (`close_model`) | `827b696c84aa529a350fc41b32246ff5e089c36f2daaeb5f63a21c168ce0cea6` |
| Recorded R132 whole-model diff `model/reports/r132_mito_adh_20261012/whole_model_diff.json` | `b98def5bdeedb9cad91677ce1468d0856d508faf25e79509971756b527d224cd` |
| PO1f profile `$RESEARCH/state/strain_profiles/po1f_sd_leu.json` | `35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383` |
| c98961e9 (other worktree) | `c98961e960c6a0f8283f485e32daf20b04dd70d998fad86fc09dca0d14ccae4c` |
| 51c6fabd (other worktree) | `51c6fabd486f4e1a9d6365bf7740bde3eca20fa05f1d5454a73586d3276db3fa` |
| `selection_gpr_restore.py` (clever-kare, uncommitted) | `9bf0f5a63f246f7c122e0c3a9110edf335dccb6c6464809c9b1e96adbc276be7` |

`IYALI26_RESEARCH_ROOT` is unset in this session. The PO1f profile was read from
`/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`, which is the folder the recorded
runs point to. Nothing was built, so the builder's fail-closed rule was not triggered.

## Scope and budget

- **Static only:** SBML loads with cobra 0.30.0 (main checkout `.venv`, Python 3.13.5), stoichiometry
  sums, element and charge sums, and reading of recorded LP outputs. No LP, FVA, build or screen was run.
- **External look-ups:** ENZYME (ExPASy), Rhea, KEGG REST, and the eQuilibrator web page. One
  literature full text was attempted.
- **Outputs:** this TASK.md and REPORT.md. The analysis scripts are in the session scratchpad
  (listed in REPORT §8). They were not added to the repository, because a source file needs the
  user's approval under the user's instructions.
- **Stop:** any edit to model, data or code would be out of scope. STATE.md was not edited
  (REPORT §9).
