# TASK: the 29 "not solvable by curation now" false negatives — sphingolipid, GPI, cardiolipin, N-glycan and mitochondrial-translation demands (2026-11-02)

## Scope

**Investigation only.** No model, curation, GPR, medium, benchmark or code change was made. Every variant was built in
memory and discarded. This folder files the finished scoping work in the repository. Filing it was itself
documentation and one figure computation (581 LPs, below); the report's numbers were not re-run.

**Authorization.** The user, 2026-10-30, in the gap-census session:

> can we solve this problem: Sphingolipid, GPI, cardiolipin, N-glycan and mitochondrial-translation demands (29 FN):
> each needs a pathway connection, a composition value, or reversal of a recorded decision.

- **The 29 genes:** essentiality false negatives. Each is on the consensus essential list but predicted non-essential
  at every cutoff.
- **Where they come from:** the "Essential gene prediction improvement" session put them in its "not solvable by
  curation now" pool.
- **How the gap session answered:** with the multi-agent scoping workflow below.
- **This filing:** the gap session handed the filing to this session (documentation only; commit, push and PR only
  on the user's explicit instruction).

## Starting state

- **Branch:** `claude/fn29-demand-families-report-a4dc68`, cut from `origin/main` at `80cfc4f` (Merge PR #44),
  fetched 2026-11-02.
- **Model executed** (every scout and verifier): `e1e7c868` =
  `model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_20261029.xml`,
  SHA-256 `e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897`.
  - The scouts read a byte-identical scratch copy `e1e7c868.xml` (archived; the SHA was checked 2026-11-02 against
    both files).
- **Condition:** `cond_taskA_po1f.json`, SHA-256 `db78d70531a447535147c7e462d7d9add9589ec38345d15ec7f894ee29023a63`.
  - Task A medium (35 uptake bounds) and R612 (0, 0). `xMAINTENANCE` keeps the model's lower bound 7.8625.
  - The file omits the PO1f plasmid LEU2 complement of R45, so YALI1C00464g scores as an extra FP (REPORT.md
    section 3.2).
- **Solver:** Gurobi 13.0.1 / 13.0.2 through cobra 0.30.0, tolerance 1e-9, real HOME licence.
- **Experimental list:** `$IYALI26_RESEARCH_ROOT/state/essentiality/repository/consensus_essential_genes.csv`,
  SHA-256 `1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b`.
  - 1,612 genes, 326 of them in the model.
  - Model genes not on the list are the user-defined negative class.
- **Reference model:** none designated (`model/candidates/README.md`). `e1e7c868` was already superseded on
  `origin/main` when the work ran. See REPORT.md, "Status on the current main".

## Workflow design (gap session, workflow `wf_f4df58df-23f`)

| Phase | Agents | Role |
| --- | ---: | --- |
| Scout | 5 | One per family. Each traced the mechanism, designed and tested variants in memory, ran full knockout sweeps and searched the literature |
| Refute | 15 | Three lenses per family: LP reproduction (an independent script and its own LPs), biology and evidence, and project rules and tuning |
| Synthesize | 1 | Consolidated report: verdicts win over scout claims unless their own evidence is weaker |

- **Journal** (`workflow_journal.jsonl`, 77 events): 42 agent launches for these 21 roles, 26 result events and 8
  failure events.
  - Failures and re-launches: the cardiolipin and N-glycan scouts, the GPI and sphingolipid biology lenses, all three
    mitochondrial-translation lenses, and the first synthesis.
  - Every one of the 21 roles has a final result.

**LP counts.** The figures below are from REPORT.md section 5; the reconciliation with the JSON is this filing's.

| Family | LP lens (independent) | Biology-lens spot check | Scout (as recorded in its JSON) |
| --- | ---: | ---: | ---: |
| Sphingolipid | 1,433 (`verify2.json` 1,417 + `verify2b.json` 16) | 12 | 2,617 (sum of the per-script counts that record one) |
| GPI anchor | 1,316 | 8 | 2,388 (`s7_summary.json`) |
| Cardiolipin | 1,342 (`verify_cardiolipin.json` 1,320 + `extra_checks` 22) | 10 | 2,378 (`13_summary.json`) |
| N-glycan | 837 (`verify_c2.json` 830 + probes 7) | 9 + 2 | 2,794 (`t2_summary.json`) |
| Mitochondrial translation | 1,176 (`verify2.json` 1,163 + `extras.json` 13) | 44 | 2,462 (`SUMMARY_20261102.json`) |

The LP lens kept each family under a 1,500-LP budget.

## Deliverables of this filing

- [REPORT.md](REPORT.md): the synthesis report as delivered, with the edits listed in its header.
- [tables/genes_29.csv](tables/genes_29.csv) and [tables/families_5.csv](tables/families_5.csv): parsed from the
  archived JSON by [scripts/make_tables.py](scripts/make_tables.py) (no LP).
- [escher/](escher/): the status figure of the five families, with each recommended variant as a what-if (581 LPs).
- **Archive:** `$IYALI26_RESEARCH_ROOT/artifacts/tasks/fn29_demand_families_20261102/` (317 files in
  `MANIFEST.sha256`, 29 MB).
  - the scratchpad `fn29/` folder (scripts, JSON, logs, the scratch model copy, the condition file);
  - `verify_mt/` (mitochondrial-translation biology lens);
  - `verifier_downloads/` (`fn29_review/`, `critic_dl/`);
  - `workflow_journal.jsonl`;
  - `REPORT_fn29_demand_families_as_delivered.md`.
  - The two `verifier_downloads` folders are attributed to this workflow by content and timestamp (inferred). They
    were in the gap session's scratchpad, written inside the workflow's window, and hold the Yeast-GEM, yeast 7.6,
    KEGG and UniProt files the biology lenses cite.

## Standing constraints

- **No model edit:** no model XML is edited by hand or rebuilt here. Every family fix needs its own TASK.md through
  the curation-change flow.
- **Sweep counts:**
  - They were produced on `e1e7c868` under the pre-2026-10-30 scoring convention, or a mixture of conventions.
  - The rule-B counts are inferences, not runs.
  - They must be re-run on a named base before any TASK.md.
- **Lipid-unlump ownership:** the sphingolipid, GPI and cardiolipin terms change `biomass_C` lipid composition, which
  is the colleague-owned lipid-unlump territory. Review and recommend only. The cardiolipin biomass weight is recorded
  there as `blocked`.
- **Gates:** the CoQ9 and lipid-unlump approval gates stay in force. A coefficient may never be chosen because it
  flips a gene (AGENTS.md §2 and §6).
- **Not this filing's work:** the bucket-4 session ("Disconnected-pathway essentiality FNs") holds the reconciled
  draft for the overlapping genes. Nothing from it is duplicated here.
- **Publishing:** commit, push and PR need the user's explicit instruction. A PR goes to `main` with
  `--repo yhong075-sketch/iyali26_gem`.
