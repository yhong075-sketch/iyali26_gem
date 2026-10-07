# Latest-model single-gene screen: bounded execution authorized

This directory contains the bounded job for the requested atlas update. The parent task explicitly relayed the user's approval of its asynchronous authorization question: latest R153 merged model, original SD-Leu/PO1f, one WT plus 1073 single-gene KOs, one solver process, ten-minute maximum, stop on invalid results, no retry, no source changes, and separate outputs. This is a relay of the authorization scope, not a verbatim user quotation. Execution is authorized once preflight passes. No optimizer is called by preparation or `--preflight`.

## Fixed input and scope

- Actual input: `../inputs/model.xml`, a frozen byte-identical copy of `model_metadata_trna_r153_merged.xml`.
- SHA256: `10a3baa1da5eb8afb80f93c4cf8a43c3613af7ebb92c8e3385dcda44b91695d0`.
- 2314 reactions, 1877 species, 1073 source geneProduct entries. Use the complete exact source ID set, including retained unassociated/placeholder entries; do not silently redefine it as native protein-coding genes.
- Exclude only the additional runtime marker `PO1f_plasmid_LEU2` from screening. This marker represents model plasmid complementation of R45 and is not an original source gene or independently verified native protein identity.
- Reuse the frozen `essentiality_simulation_context` and strain-overlay source files from `artifacts/r608_engineering_20260907/code`, the SD-Leu medium and PO1f profile recorded by the ordinary 2026-09-10 screen and the 2026-09-15 two-gene test.
- Before optimization, require exact agreement with the latest two-gene test's effective medium, strain overlay, context fingerprint, all reaction bounds, biomass_C objective, software versions and solver parameters. Preserve current code/dependency identities and Git dirty state. Matching selected sources does not reconstruct the entire old environment.

## Budget and stop conditions

Exactly one WT plus 1073 independent single-gene KOs: at most 1074 primary LP calls. One solver worker process, Gurobi Threads=1, TimeLimit=60 seconds per LP. A lightweight supervisor kills that worker at 600 seconds total. No multiprocessing, retries, FVA, pFBA, extra rescue/control solves or automatic presolve retry.

Stop on input/code/dependency/context drift, solver-parameter mismatch or `presolve='auto'`; missing/duplicate screened IDs; incorrect gene-to-GPR-to-bound propagation or restoration; nonoptimal, missing, nonfinite or negative raw growth; nonfinite fluxes; mass residual or bound violation >=1e-7; WT outside the inherited [0.1,2.0] h^-1 acceptance range; LP or time budget exhaustion. Preserve the completed raw rows and incomplete execution record. A failure is not a biological essentiality label. No negative-value clamping or altered tolerance is permitted.

## Outputs and interpretation

`run_manifest.json` records full input/source hashes, software, current Git state, context, objective, effective bounds, solver parameters, WT, progress and actual primary-call count. `raw_deletions.tsv` is flushed after every returned KO, including an invalid returned row before stopping. It retains original numeric/repr values and status plus GPR closure IDs and numerical residual summaries. Complete flux vectors are checked transiently and are not saved by this job.

Only a complete, valid, exactly matched 1073-ID screen produces `screen_predictions.tsv` and `screen_summary.json`. Principal essentiality is the strict comparison **unrounded KO/WT < 0.15**; equality at 0.15 is nonessential. The 1%, 5% and 10% comparisons are arithmetic postprocessing of the same raw solves. No experimental label workbook is read or remapped by this job, and no biological validity or independent validation claim is added.

`execution.json` and `screen.log` preserve supervisor status and runtime output. Existing execution/result files prevent a second run or silent overwrite. Model, GPR, chemistry, medium, original evidence and Git state are not modified.

## Commands

From the repository root, preparation only:

```sh
.venv/bin/python artifacts/escher_atlas_20260915/screen/run_screen.py --preflight
```

After explicit execution authorization is relayed, the one-run command is:

```sh
.venv/bin/python artifacts/escher_atlas_20260915/screen/run_screen.py --execute-authorized
```

Preflight uses only XML/JSON parsing, checksums, installed distribution metadata and strict-threshold self-checks. It does not import COBRA or Gurobi or instantiate/optimize a model. Runtime loader/solver/LP checks remain pending until authorization.
