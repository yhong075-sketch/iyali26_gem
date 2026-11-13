**Status: done; see [REPORT.md](REPORT.md) and [AUDIT.md](AUDIT.md). Candidate v2 `85c623c8`.**

# Folinic-acid (5-formyl-THF) module: R748 to SHMT, R990 closed, cytosolic R_THFAT added, 2026-11-12

One `curation-change` round.

## 0. Authorization

The user, 2026-11-12, in this session:

1. **The gap decision.** After the review
   [folinic_acid_module_review_20261112.md](</Users/david/Library/Application Support/Claude/side-session-notes/local_e1e2b48d-5e68-4eac-b28d-10eb1c773eea/folinic_acid_module_review_20261112.md>)
   (shared session notes, outside the repo), the user wrote: "I think there is a gap to fix".
2. **The change.** "Authorize F1+F2+F4".
3. **The plan.** "Approve plan", for the implementation and compute plan below.

F3 (the gene on R989) stays with the lipid-unlump line, which is owned by a colleague. Not in scope:
- the methenyl-THF charge;
- other annotations of R748 and R542;
- the 5-formyl-THF identifiers carried by the 10-formyl-THF metabolites.

No commit or push without a separate instruction.

## 1. Starting state (verified 2026-11-12)

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/eloquent-ishizaka-954d88` @ `154eb64c` (PR #48 merge). It was `main` when this session started. Local `main` moved to `efd7f780` (PRs #49–#51) before the builds ran; see REPORT §7 |
| Baseline model | `model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_bigg_annot_v2_mcnaugh_rxns_v3_20261103.xml`, `e584151fc79475fc2daa59e76c8c84ebfd06b0f3cb1fbd75fbb69cc9f31039d2`. Not a designated reference. |
| Baseline build arguments | that candidate's `.build.json` argv, in [scripts/args_e584151f.txt](scripts/args_e584151f.txt) |
| Medium / profile | Task A medium `6eca70fd…`, PO1f profile `35307853…`, via `essentiality_simulation_context.py` |
| Essential list | `consensus_essential_genes.csv` `1e887f5a…` |

Targets in `e584151f` (from [before.json](before.json)):

| Item | Before | After |
| --- | --- | --- |
| R748 [C_mi], methenyl-THF m803 + H2O → 5-fTHF m969 + H+, (0, 1000) | rule `YALI1F04082g`; `ec-code` 2.1.2.10; name "tetrahydrofolate aminomethyltransferase"; note `PROTEIN_CLASS: 2.1.2.10` | rule `YALI1D28569g or YALI1E19557g` (as R359); `ec-code` 2.1.2.1; name "5,10-methenyl-THF hydrolysis (SHMT side activity)"; `PROTEIN_CLASS: 2.1.2.1` plus this round's notes. Stoichiometry, bounds and other annotations unchanged. |
| R990 [C_cy], ATP + 5-fTHF + H2O → ADP + Pi + 10-fTHF, gene-less | bounds (0, 1000), no notes | bounds (0, 0), plus this round's notes. Kept in the model. |
| R_THFAT [C_cy] | absent | methenyl-THF m160 + H2O m32 → 5-fTHF m159 + H+ m10, (0, 1000); rule `YALI1D28569g or YALI1E19557g` (as R358); group and subsystem "One carbon pool by folate"; annotation SBO:0000176, BiGG THFAT, EC 2.1.2.1, Rhea 34767, KEGG R02300, MetaNetX MNXR153595 |
| `YALI1F04082g` reactions | R357, R360, R361, R363, R365, R748 | R357, R360, R361, R363, R365 |
| `YALI1D28569g`, `YALI1E19557g` reactions | R358, R359 | R358, R359, R748, R_THFAT |

Nothing else may change.

## 2. Expected effect and acceptance criteria

1. A control build (same argv, no new flag) is byte-identical to `e584151f`.
2. The whole-model diff against `e584151f` shows only the targets above. Gene-reaction links follow from the rules.
3. R_THFAT is mass-balanced. Its charge imbalance (−1) equals R748's. The model stores methenyl-THF with charge 0,
   but its formula C20H20N7O6 is that of the −1 species
   ([ChEBI 57455](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:57455), verified 2026-11-12). The charge is
   not changed in this round.
4. A candidate built with seed 1 is byte-identical to the seed-0 candidate.
5. Solves on Task A + PO1f, Gurobi 1e-9, basis reset per solve, on both models:
   - **WT growth:** unchanged within 1e-9.
   - **FVA** of R748, R989, R64, R990, R811, R991, R_THFAT and R543, at 100 % of WT growth and with no growth
     floor. Expected:
     - R64 is no longer 0/0 with no growth floor;
     - R990 is 0/0;
     - R_THFAT, R64, R748 and R989 are ≈ 0 at 100 % (inferred: they only spend ATP or export folate).
   - **Targeted knockouts:** `YALI1F04082g`, `YALI1D28569g`, `YALI1E19557g`, `YALI1D28569g`+`YALI1E19557g`,
     `YALI1D01489g`.
   - **Full single-gene sweep** on both models, scored at 1/5/10/15 % KO/WT. Expected: no class change.
     F4 makes 5-fTHF exportable through the existing R991, so a change cannot be ruled out without running it.
6. Full pytest (root `.venv`, short-path archive copies): only the new tests differ, test by test.
7. An independent read-only audit leaves no unsupported claim unaddressed.

## 3. Budget and stop condition

- **Builds:** three (control seed 0, candidate seeds 0 and 1).
- **Diff:** one whole-model diff.
- **Solves:** about 2,200 sweep LPs plus about 40 targeted LPs and FVA LPs, under about 1 hour of wall time.
- **Tests and audit:** two pytest runs, one audit agent.
- **Stop and report** if any of these happens:
  - a solve is not optimal;
  - WT changes by more than 1e-9;
  - the diff shows an unintended entity.

  Otherwise stop after REPORT.md and the STATE.md entry.

Large outputs (control and candidate copies, sweep tables) go to
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/folinic_acid_module_20261112/`.

## 4. Addendum: rebuild on current main (2026-11-12)

**Authorization.** The user's choice "Merge main, rebuild", in reply to "Merge current main and rebuild this change on
top of it?".

**Scope.**
1. Commit this round's work on this branch (`c00ec033`).
2. Merge `origin/main` `cf7b8dd9` (`768c187e`).
3. Rebuild on main's newest candidate `684f3752`:
   - control build;
   - candidate seeds 0 and 1;
   - whole-model diff;
   - the same targeted solves and full sweeps on both models;
   - pytest before (an `origin/main` archive) and after (a merged-HEAD archive).

Acceptance criteria as in §2, with `684f3752` in place of `e584151f`. No push.
