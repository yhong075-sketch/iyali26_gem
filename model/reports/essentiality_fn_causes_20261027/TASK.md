# Task: why are consensus-essential genes predicted non-essential? (investigation only)

Date: 2026-10-27. Worktree `task-b-e1ff88`, branch `claude/essentiality-fn-causes-c36ef9` at
`5fb75e7ed4300df02c8ca8f1f2c5d7f3fb79d86d` (25 commits behind `origin/main` `7ebaa0f55af6b3721f28db5aad209f277bd88fec`
at the start; the platform code used here differs from main only in builder steps and default-solver
flags that this task does not touch).

## Authorization

Brief from the parent session (2026-10-27), relayed from the user: *investigation only, no model fixes*.
Find out why genes on the experimental consensus essential-gene list are predicted non-essential by the
iYali26 model (the false negatives) and explain the causes per gene. Use the literature reviewer
(`platform/gem_review/lit`) where the biology matters. Do not change any model XML, GPR, medium, benchmark
or curated data; do not commit or push unless the user asks. Deliver `model/reports/essentiality_fn_causes_<date>/`
(TASK.md, REPORT.md, per-gene CSV) and a dated `model/STATE.md` entry; large outputs to
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/<same folder>/`. Never tune anything to improve recall.

Computation run under this brief (named inputs, scope, budget, outputs, stop condition):

| Step | LPs | Stop condition |
| --- | ---: | --- |
| Single-gene KO screen of a1237821 (`run_screen.py`) | 1,082 | call limit 1,082 / 900 s; input SHAs re-checked |
| Mechanism classification (`classify_fn.py`): blocked set, WT FVA, reaction-level KO, pFBA uptake shift | 5,829 | one pass over the 212 FN genes |
| Controls (`trace_controls.py`): medium controls on all 212 FN, dead-end check, FVA bypass trace on 16 genes | 74,582 | one pass |
| Extra traces (`trace_extra.py`): GPR texts, biomass, ENO/PYK traces, vitamin-closed cofactor KOs, bypass-closed controls | 50,990 | one pass (run twice; the second added the ERG10 control) |
| Literature: Europe PMC retrieval on 19 reactions / 56 pairs (free), LLM extraction on 8 pairs, 5 logic chains | n/a | pairs with hits only |

## Starting state (identified this session)

| Item | Value |
| --- | --- |
| Model executed | `E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_20261023.xml`, SHA-256 `a12378210c7c66720d99fb7bf71d14cb7538ea6a3cf3e08919c85fb3c9c47d82`. Newest candidate in `model/candidates/README.md` on `origin/main` (2026-10-23, merged in PR #19). Not in this worktree's stale README (its newest is `cf0825b8`, 2026-10-21); exported byte-identically with `git show origin/main:` into `$IYALI26_RESEARCH_ROOT/artifacts/tasks/essentiality_fn_causes_20261027/input/` and SHA-checked. No reference model is designated. |
| Reference model | none designated (`model/candidates/README.md`) |
| Medium | Task A medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv`, SHA-256 `6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45` (glucose 2.414; 11 amino acids, adenine and iron capped; uracil 0.01607 raised to 1000 by the profile) |
| Strain profile | PO1f `state/strain_profiles/po1f_sd_leu.json`, SHA-256 `35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383` (R612/URA3 closed, R45 LEU2 plasmid-complemented), applied in memory |
| Experimental list | `state/essentiality/repository/consensus_essential_genes.csv`, 1,612 genes, SHA-256 `1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b` |
| mtDNA scope rule | `model/expected/mtdna_genes_out_of_scope.json`, SHA-256 `16b3e9b4418f24821b8775936fd377ae21f77afa8e8af2b8949a246f62c44996` (16 genes excluded) |
| Scoring | `model/expected/benchmark.md`: essential when KO/WT < cutoff (strict); cutoffs 1, 5, 10, 15 %, 15 % primary; a model gene on the list is a positive, a model gene not on the list is a negative (user-defined negative class); listed genes without a model counterpart are out of scope, never FN; failed or infeasible solves are unresolved, never "dead" |
| Solver | Gurobi 13.0.1, FeasibilityTol/OptimalityTol 1e-9, Presolve 0, one thread, 60 s per LP; context reset per knockout; knocked-flux check 1e-9 (copied from `reports/screen_2c9cebd0_20261022/run_screen.py`) |
| Interpreter | `.claude/worktrees/workflow-automation-review-6765da/platform/.venv/bin/python` (Python 3.13.5, cobra 0.30.0; no scipy, so the mass-residual check uses a dense matrix) |
| Prior work checked | `git branch -a`: only `codex/dfba-new-fn` besides this branch; `$IYALI26_RESEARCH_ROOT/artifacts/tasks`: single-KO diagnoses of R153, R190/R607/R152, R2081, R2202, R608 and the V-ATPase (2026-09), screens up to `screen_2c9cebd0_20261022`; STATE.md entries up to 2026-10-26 on `origin/main`. The AGENTS §6 FN loop tooling (`--prepare-agent-cases`, case ledger) is absent on this branch and on `origin/main`; the loop was not used. |

## Not done

- No model, GPR, medium, benchmark, curated-data or test change. No commit, no push.
- The literature reviewer restricts every query to *Yarrowia lipolytica*; general yeast biology that only
  exists for *S. cerevisiae* returns `no_hits` and is marked as such.
