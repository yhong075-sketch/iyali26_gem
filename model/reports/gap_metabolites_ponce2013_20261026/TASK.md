# TASK: gap metabolites and blocked reactions of the newest candidate a1237821, by the Ponce-de-León 2013 definitions (2026-10-26)

## Authorization (user, this conversation)

"use the definition in this paper, and find the gap in our model", with the paper attached:
Ponce-de-León M, Montero F, Peretó J (2013) *Solving gap metabolites and blocked reactions in
genome-scale models: application to the metabolic network of Blattabacterium cuenoti*. BMC Syst
Biol 7:114, doi:10.1186/1752-0509-7-114. Computation and reporting only; no model change; no commit.

## Starting state

- Branch `claude/gap-metabolites-model-2c880a` (worktree), HEAD `aeba460` (`main` after PR #22), clean.
- Research root `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`.
- No reference model is designated (`candidates/README.md`). "Our model" is taken as the newest
  candidate on `main`.
- No earlier Ponce-de-León-style gap deliverable exists in `model/reports/` or `model/STATE.md`
  (checked by grep for RNP, "unconnected module", "blocked reaction", "gap metabolite").

## Model

- `candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_20261023.xml`,
  SHA-256 `a12378210c7c66720d99fb7bf71d14cb7538ea6a3cf3e08919c85fb3c9c47d82`, the newest candidate on
  `main` (commit `11450bc`, PR #19 merged 2026-10-26). 2,316 reactions, 1,882 metabolites, 1,081
  genes, 182 exchanges, 3 demands, 69 reactions closed at (0,0).
- Pipeline regression check on `candidates/E5_coq9_resp_pkg_r132_20261012.xml`
  (`a56f3649…`), the model on which the analysis skill's reference numbers were recorded.

## Definitions used (the paper's)

- Blocked reaction (eq. 4–5): FVA minimum = maximum = 0 under the condition. Primary condition is
  every exchange open at (−1000, 1000), the paper's "any medium condition"; the Task A medium +
  PO1f profile (`po1f_csm_leu_w29uptake_20261012.csv`, uracil override 1000, R612 closed) is
  reported separately.
- Gap metabolite (eq. 6–7): every reaction in which it participates is blocked.
- RNP / RNC: root non-produced / root non-consumed, by a row scan of the stoichiometric matrix with
  bounds and reversibility; DNP / UNC: downstream non-produced / upstream non-consumed, by
  propagation after removing the roots' reactions. The paper's plain scan counts are given, plus a
  split of isolated metabolites and metabolites in a single reversible reaction (RNP + RNC).
- Coenzyme pseudo-gap (paper Fig. 2): a cofactor that is synthesised but never drained; its
  biosynthesis is blocked although the cofactor itself is not a gap.
- Unconnected module (UM): a connected component of the bipartite graph of gap metabolites and
  blocked reactions.

## Method, budget, stop

- Scripts of the `gem-gap-rnp-rnc-analysis` skill (`gapfind.py`, `fixtests.py`, `kernel.py`,
  `closure_check.py`), copied to the artifacts folder below; cobra 0.30.0 with GLPK. Gurobi 13.0.2 at
  FeasibilityTol / OptimalityTol 1e-9 is used only as an independent cross-check of the blocked set.
- Budget: three FVA scans (2 × 2,316 LPs each: a1237821 open, a1237821 Task A + PO1f, a56f3649 open),
  the Gurobi cross-check (2 × 2 × 2,316 LPs), elementary single-fix LP tests (each re-tests only
  its module's blocked reactions), cofactor drain tests and a few targeted module tests. Stop when
  the scans and tests complete, or on a solver failure.
- Outputs: this folder (REPORT.md, small tables, figures); scripts and full tables in
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/gap_metabolites_ponce2013_20261026/`.
- Fix tests are LP evidence ("adding X unblocks n of N"), not curation recommendations.

## Addendum (user, same conversation, 2026-10-26)

"can you use the Escher to draw them? I need to have a visible map to figure out these problem." → Escher maps
of the 20 largest unconnected modules (UM1–UM20, 436 of the 756 blocked reactions): an offline interactive
viewer with module and scenario drop-downs, gap-class node colours and the census results as reaction data,
plus SVG/PNG renders. No new LP solves and no model change; scenario values are the fix, drain and targeted
tests already recorded in this report. Layout is automatic (layered chains, repeated hub metabolites).

Second addendum (user, same conversation): "Draw one Escher map per problem class (cofactor pseudo-gaps,
inherited closures, missing transports, dead-end roots) that gathers the reactions of that class across modules,
so each class can be reviewed on one map." → four class maps added to the same viewer, with the tested drains,
openings, links and sinks drawn as labelled TEST reactions (not model content). No new LP solves.
