# TASK: inactive-reaction census of the newest candidate 7ccbb568, and LP tests of single-reaction candidates for every root metabolite (2026-11-16)

## Authorization (user, this conversation)

- The question: "I wanna that we could make a GNN that answer a question that, how can we activate more
  inactivate reaction ? In other words, we need to train a GNN to answer why a reaction can not carry the
  flux? two consider:1, is it must be muted due to steady state? 2, is it some reactions is missing so it
  cann't be carry the flux?"
- To "Write the plan for the two checks that decide whether to build the GNN?": "Write the plan".
- To the plan stated in the conversation (check 1: the gap census re-run on 7ccbb568; check 2: an LP test of
  every single-reaction candidate per root; budget 150,000 LPs or 3 hours; two new scripts in this folder; one
  `model/STATE.md` entry; one read-only review): "Run them".

Computation and reporting only. No model, curation, medium, benchmark or test change; no commit; no push.
The approval covers four choices stated in the plan: the model below; single added reactions only; a
pre-set reading cut of 30 roots; scripts kept in this folder.

## Starting state

- Branch `claude/metabolic-pathway-gnn-priorities-9b8f18` (worktree), HEAD
  `af7c543214459f22d3947d71c8fc1a75ea5fb800`, clean. `origin/main` is the same commit (`git ls-remote`,
  2026-11-16).
- Research root `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`.
- No reference model is designated (`candidates/README.md`). "Newest model" is the newest candidate on `main`.
- Earlier work reused: the census of 2026-10-26 on `a1237821`
  ([report](../gap_metabolites_ponce2013_20261026/REPORT.md); 756 of 2,316 reactions blocked). Later
  whole-model counts recorded in `model/STATE.md`: 722 (`00e8e397`, 2026-10-27) and 695 (`14f73bda`,
  2026-10-29). None is recorded for `7ccbb568`.

## Inputs

| Input | Path | SHA-256 |
| --- | --- | --- |
| Model executed | `candidates/E5_353cefdc_xref_20261115.xml` | `7ccbb56834a37ad4c75e434868e9e1b9e1f6a657a11ca17ee0a93dbe56232146` |
| Regression model | `candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_20261023.xml` | `a12378210c7c66720d99fb7bf71d14cb7538ea6a3cf3e08919c85fb3c9c47d82` |
| Task A medium | `conditions/media/po1f_csm_leu_w29uptake_20261012.csv` | `6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45` |
| PO1f profile | `$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json` | `35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383` |
| MetaNetX/MNXref 4.5 (2025/08/13) reactions | `$IYALI26_RESEARCH_ROOT/reference/metanetx/reac_prop.tsv` | `8582cc187d03ce127f8e914f8f298282ac9036f1448918117af044ac55b980db` |
| Gene-based gap-fill tables | `curation/gapfill/gap_fill_candidates.csv`, `gap_fill_prioritized.csv` | `29f257bf…`, `9412dd06…` |
| Census scripts | `$IYALI26_RESEARCH_ROOT/artifacts/tasks/gap_metabolites_ponce2013_20261026/scripts/` | copied unchanged; hashes in the artifacts folder |

`7ccbb568` is `353cefdc` plus the reaction cross-reference cleanup; annotations are not in the LP. It is not a
designated reference.

## Definitions

As in the 2026-10-26 census (Ponce-de-León, Montero & Peretó 2013, BMC Syst Biol 7:114): blocked reaction, gap
metabolite, root metabolite (RNP, RNC, single reversible reaction, isolated), unconnected module, coenzyme
pseudo-gap. Primary condition: every exchange open at (−1000, 1000). Secondary: Task A medium + PO1f profile.

A *single-reaction candidate* for a root is a MetaNetX non-transport reaction that contains the root's
compound, whose every compound already has a species in the root's compartment, and whose substrates and
products (H+ ignored) are not those of a reaction already in that compartment.

## Method, budget, stop

- **Check 1.** The census scripts unchanged (`gapfind.py`, `analyze_roots.py`, `fixtests.py`,
  `build_drain_tests.py`, `fva_gurobi.py`), cobra 0.30.0 with GLPK; Gurobi 13.0.2 at 1e-9 as a cross-check of the
  blocked sets. New `scripts/census_tables.py` assigns one cause per blocked reaction in the October order and
  lists the changes since October.
- **Check 2.** New `scripts/candidate_tests.py`: adds each candidate in memory as a reversible reaction and
  solves for its largest forward and backward flux, alone and with every other root opened; counts the blocked
  reactions of the root's module that become active for candidates that work alone. GLPK and Gurobi (1e-9).
- **Gates (stop on failure).** File hashes; the Task A condition file equals the platform's own medium + profile
  context; the scripts reproduce October on `a1237821` (756 blocked, 701 gap metabolites, 212 modules); the root
  table of `7ccbb568` equals the structural count made at plan time (364: RNP 98, RNC 101, single reversible 137,
  isolated 28); the candidate list equals the plan-time count (391 pairs, 400 placements).
- **Budget.** About 95,000 LPs for check 1 and at most 30,000 for check 2. Hard cap 150,000 LPs or 3 hours.
- **Stop.** When the scans and tests complete, on a gate failure, on the cap, or on a solver failure. A failed
  or non-optimal solve is recorded as unresolved, never as blocked or as a failed candidate.
- **Pre-set reading.** Count the roots with two or more working, balance-marked candidates in the second pass
  (other roots opened). Fewer than 30: no ranker for single reactions. The cut of 30 is the assistant's choice,
  accepted with the plan.

## Outputs

- This folder: REPORT.md, small tables, `scripts/`.
- `$IYALI26_RESEARCH_ROOT/artifacts/tasks/gap_census_candidates_20261116/`: census script copies, raw outputs,
  logs, full tables.
- One dated entry in `model/STATE.md`.
- Every fix test and candidate test is LP evidence, not a curation recommendation.

## Addendum (user, same conversation, 2026-11-16)

After the first hand-back the independent review's code findings were listed as proposed fixes. The user asked
"why you want to test in every candidate?" about the proposed loop test; to "Apply the other fixes and leave the
loop test out?" the answer was "Yes, without it". Scope of the addendum:

- **Hash gate inside both scripts.** New data file `input_hashes.csv` (values from this TASK.md, from the
  pre-run hash listing and from the recorded hashes of the October census scripts); `census_tables.py` checks 12
  inputs and `candidate_tests.py test` checks 4, and the latter also requires its path arguments to be those
  files.
- **Fail-closed tables stage.** `candidate_tests.py tables` stops when the test stage recorded a consistency
  problem or a census fix test failed to solve, and reports roots with an unresolved solve in their own column.
- **One more recorded decision.** `recorded_decisions.csv` gains "R1725 stays closed" (2026-10-27).
- **One targeted test.** With the unchanged census `fixtests.py`: reopen R2121, R1172, and both, each re-testing
  the nine reactions that became blocked since October without a bounds change (54 LPs, GLPK).
- **Re-run** of the candidate test stage and of both table stages.

Declined by the user: an exact loop test per working candidate. Not changed: the test stage's handling of a
non-optimal solve inside a re-count (review finding 5, second part), which was not in the proposed list.
The first run's outputs, report tables and script versions are kept in the artifacts folder under `run1/`.

**Commit (user, same conversation, 2026-11-16):** "Commit the working tree changes with a sensible message."
This folder and the `model/STATE.md` entry are committed on `claude/metabolic-pathway-gnn-priorities-9b8f18`
(`9dd3ec17`).

**Pull request (user, same conversation, 2026-11-16):** the app's Create PR command, base `main`, repository
`yhong075-sketch/iyali26_gem`. The branch was pushed and the pull request opened; this note and the matching
lines in REPORT.md and `model/STATE.md` are a second commit on the same branch.
