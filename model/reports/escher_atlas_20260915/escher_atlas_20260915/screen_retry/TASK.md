# Revised screen: execution authorized

The first authorized attempt is preserved in `../screen/`: 41 primary LP calls (one WT and 40 single-gene KOs), all solver-reported optimal with valid numerical checks. After the first source gene with no GPR association, the all-gene restoration assertion failed. The task stopped immediately; no retry was performed. Forty raw KO rows remain, and 39 had passed the post-KO restoration check. No global full-screen essentiality count was produced.

The parent task first authorized preparation and a **zero-optimization regression check** for this revised attempt. It has now relayed the user's explicit instruction **“授权 出错后可重试”** (“authorized; retries are allowed after errors”). This authorizes executing the prepared complete screen and, if another software error occurs, preserving that attempt, fixing and testing the software, and starting a separately recorded attempt under the same fixed input, context and per-attempt budget. It does not authorize scientific input, model, GPR, chemistry, medium or numerical-criterion changes. The earlier stopped attempt remains preserved.

## Minimal software correction

COBRA's SBML reader adds geneProduct objects to the gene collection; a gene without a reaction association can have no owning `_model`. The existing context decorator therefore cannot register its functional-flag restoration. The runner now saves the target's flag, performs the ordinary `gene.knock_out()` within the model context, and restores that flag explicitly in `finally`. All existing checks that every gene and every reaction bound have returned to their previous values remain mandatory. No GPR, reaction, source gene, chemistry, model file or solver setting is changed.

`regression.py` loads the exact approved model/context without optimizing. It reproduces the original first-orphan failure, restores the temporary flag, then tests the production correction on that orphan and a GPR-associated gene. It checks all gene flags, all bounds, complete stoichiometry/GPR and objective before/after. Both COBRA and optlang optimization entry points are patched to fail if called, and the guard call count must remain zero.

## Proposed renewed execution

Repeat the complete set independently in this new directory: **one WT plus 1073 original-gene KOs, at most 1074 new primary LP calls**. Together with the stopped attempt, that would be at most 1115 LP calls across both attempts; original rows are retained as prior evidence, not silently replaced or combined into a new complete screen.

- Actual frozen input `../inputs/model.xml` is byte-identical to approved `model_metadata_trna_r153_merged.xml`, SHA256 `10a3baa1da5eb8afb80f93c4cf8a43c3613af7ebb92c8e3385dcda44b91695d0`.
- 2314 reactions, 1877 species, 1073 source gene entries, including retained placeholders and unassociated entries. Exclude only the additional runtime `PO1f_plasmid_LEU2` marker (model plasmid complementation of R45; not an original or experimentally established native protein identity).
- Same frozen SD-Leu/PO1f, source loader hashes, effective bounds/context, biomass_C objective and recorded Gurobi tolerances as the last two-gene test. Record current code/dependencies/dirty state; do not claim the full old environment was reconstructed.
- One solver worker, Threads=1, TimeLimit=60 seconds per LP, hard worker timeout 600 seconds. No automatic retries, presolve repair, FVA, pFBA or additional intervention solves.
- Stop on invalid/missing/nonfinite/negative/nonoptimal results, changed inputs/code/context/settings, failed GPR propagation or restoration, residual/bound violation >=1e-7, WT outside [0.1,2.0] h^-1, or budget exhaustion. Preserve all returned raw rows; failures are not biological essentiality labels.
- Principal classification is **unrounded KO/WT < 0.15**; equality is nonessential. The 1/5/10% columns are pure arithmetic from the same run. No experimental labels are changed or remapped.
- Preserve `run_manifest.json`, `raw_deletions.tsv`, complete-only `screen_predictions.tsv` and `screen_summary.json`, plus `execution.json` and `screen.log`. Check fluxes transiently; do not save full flux vectors. Existing execution files block a silent second run.

## Commands

Authorized preparation only, from the repository root:

```sh
.venv/bin/python artifacts/escher_atlas_20260915/screen_retry/regression.py
.venv/bin/python artifacts/escher_atlas_20260915/screen_retry/run_screen.py --preflight
```

Only after renewed explicit execution authorization:

```sh
.venv/bin/python artifacts/escher_atlas_20260915/screen_retry/run_screen.py --execute-authorized
```
