# Certified-infeasible knockouts scored as growth 0 (benchmark rule proposal, 2026-10-30)

**Status: approved 2026-10-30** (§0; D1, D2 and D3 all yes). Written before approval: this file only.

## 0. Authorization

- **Brief** from the session that started this one (2026-10-30): discuss the "infeasible complex III/IV knockouts"
  item read-only; run no LP and change nothing without explicit, scoped authorization.
- **User, 2026-10-30**, after the discussion: "I think we've hit a dead end; the key issue is how we resolve the
  problem where 16 knockouts are LP-infeasible." Then, choosing the recommended rule: "Plan this rule" (button text:
  "Write the TASK.md plan for the proven-infeasible rule: the benchmark.md wording, the screen code change, the
  re-screen of e1e7c868 with both rules reported, and the LP budget. No edits until I approve.").
- **User, 2026-10-30**, on this plan: "Approve D1–D3" (button text: "Approve the TASK.md plan with D1 yes (rule
  text as written), D2 yes (COX9 out of scope), D3 yes. Run it.").
- Not authorized: commit, push, any model change, NGAM or medium change.
- Order note: the D2 gene list is written first as a run input (`scripts/untargeted_genes_out_of_scope.json`) and
  copied byte-identical to `expected/` only after A1–A4 pass.

## 1. Starting state (verified this session)

| Item | Value |
| --- | --- |
| Worktree / branch / HEAD | `claude/cx34-infeasible-discussion-7b0c5b` @ `7aa4baf` (= `origin/main`, PR #39 merge), clean before this folder |
| Model to execute | `candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_20261029.xml`, SHA-256 `e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897`. Newest candidate of this line; **not a designated reference** (`model/candidates/README.md`). |
| Medium | `conditions/media/po1f_csm_leu_w29uptake_20261012.csv`, `6eca70fd…` (Task A) |
| Strain profile | `$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json`, `35307853…` |
| Experimental list | `$IYALI26_RESEARCH_ROOT/state/essentiality/repository/consensus_essential_genes.csv`, `1e887f5a…` |
| mtDNA scope list | `expected/mtdna_genes_out_of_scope.json`, `16b3e9b4…` |
| Benchmark rules | `expected/benchmark.md`, `26838ef8…` |
| Screen code of recent rounds | `reports/aconitase_aco2_split_20261029/scripts/ko_compare_two.py` (`bb3118a6…`) with `reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py` (`69c42690…`): Gurobi 1e-9, reset per knockout, targeted check, at most 20 fresh re-solves; **infeasible knockouts are never re-solved and score as unresolved**. |
| Old-rule result on `e1e7c868` | [aconitase ko_compare_summary.json](../aconitase_aco2_split_20261029/ko_compare_summary.json): WT 0.3116608216; 15 % 125 / 190 / 42 / 692, 16 unresolved; 20 non-optimal knockouts (16 scored + 4 mtDNA), all complex III (R305) or complex IV (R304) members. |

**Mechanism** (recalled from the 2026-10-21 STATE.md entry, run on `a56f3649`; to be re-measured here): `xMAINTENANCE`
(NGAM) has lower bound 7.8625 and is the only forced flux (re-verified on `059316ef` this session by reading the
XML). Without R305 or R304, the maximum maintenance ATP is 4.828 (2 × the glucose cap 2.414), so no flux state
exists, even at zero growth.

## 2. Problem and proposed rule

`benchmark.md` says: "A failed, missing or non-finite solve is not biological evidence; report such cases
separately." That is right for numerical failures. It has no place for a **proven** infeasibility: an LP that is
infeasible because the knockout cannot cover maintenance ATP is a model prediction ("this cell cannot sustain
itself"), stronger than growth below 15 %. The 16 genes are left out of every score although the model gives a
definite answer for them.

Rejected alternatives (discussion 2026-10-30): relaxing NGAM (turns these into growth 0.0958 / 0.0637 vs WT 0.3124,
i.e. non-essential, recalled), raising glucose (measured value), adding AOX (abstracts: inactive in log phase,
cannot compete with the cytochrome chain), counting every infeasible solve as lethal (lets numerical failures in).

**Proposed text** for `benchmark.md` § Classification, replacing the last bullet (D1):

> - Keep raw solver status and growth for every knockout. A knockout whose LP is infeasible is scored as growth 0
>   (KO/WT = 0, essential at every cutoff) only when the infeasibility is **certified as a maintenance shortfall**:
>   1. the sweep solve and a fresh re-solve (new context) both return `infeasible`;
>   2. `xMAINTENANCE` is the only reaction with a lower bound above 0 or an upper bound below 0 in that context;
>   3. in the same knockout context with the `xMAINTENANCE` lower bound set to 0, maximizing `xMAINTENANCE` is
>      optimal, passes the targeted check, and its maximum is below the original lower bound by more than 1e-6;
>   4. WT is optimal in the same context.
>
>   Record the raw status, the maximum maintenance flux and the margin for each such knockout. Any other failed,
>   missing or non-finite solve is not biological evidence and stays unresolved, reported separately. Reports
>   give the counts under the earlier rule (every non-optimal knockout unresolved) alongside (decided 2026-10-30).

Why this is not a score-driven change: it follows from what infeasibility means (no non-negative growth exists),
not from the direction of the counts; it adds false positives as well (5, or 4 with D2); it applies unchanged to
any future infeasible knockout; numerical failures still cannot count; and the old-rule counts stay in every report.

**Expected effect, stated before the run** (arithmetic only, inferred): the 16 become 11 TP + 5 FP at every cutoff;
15 %: 125 / 190 / 42 / 692 → 136 / 190 / 47 / 692 (→ 136 / 190 / 46 / 692 if D2 = yes).

## 3. Decisions for the user

- **D1.** Adopt the rule text in §2 (or edit it).
- **D2.** `YALI0F04114g` (COX9) out of scope, like the mtDNA genes: no guide of either acCRISPR library lies in its
  locus (verified 2026-10-30, `$IYALI26_RESEARCH_ROOT/artifacts/tasks/cox9_blast_20261030/`, run 7), so its
  negative label is a default, not a measurement. If yes: a new `expected/untargeted_genes_out_of_scope.json`
  (one gene, with the evidence path) and one sentence in `benchmark.md` § Labels. Yes / no / later.
- **D3.** Future rounds copy this round's script instead of `ko_compare_two.py`. Yes / no.

## 4. Method (after approval)

1. **Script** `scripts/ko_screen_certified.py`: a one-model copy of `ko_compare_two.py` (`bb3118a6…`), same
   `checked_solve.py` (unchanged, via `PYTHONPATH`), same sweep loop, tolerances, `FRESH_LIMIT` 20 and timeout.
   Added:
   - `certify(path, gene)`: fresh context, knock out, `solve()` → must be `infeasible` (condition 1); list forced
     fluxes (condition 2); set `xMAINTENANCE` lb 0, objective max `xMAINTENANCE`, reset, solve, `check()` with the
     knocked reactions (condition 3). Reads the lower bound from the model; no hard-coded 7.8625. Returns status,
     forced-flux list, max maintenance flux, margin, check magnitudes and the verdict.
   - The same maintenance-ceiling LP on WT (context, expected far above 7.8625).
   - Scoring under rule A (current: infeasible = unresolved) and rule B (certified infeasible = growth 0), each at
     1/5/10/15 % with TP/FN/FP/TN, unresolved, recall, specificity, precision, MCC; mtDNA excluded as now; with and
     without COX9 for D2.
2. **Wrapper** `scripts/run_screen.zsh`: SHA-checked short-path copy of the model (cobra 260-char limit), input
   hashes to `logs/`, real HOME (Gurobi licence), `/opt/miniconda3/bin/python`, run in the background.
3. **Independent read-only review** of the script (general-purpose agent on `.claude/agents/code-reviewer.md`)
   before the run.
4. **Run** on `e1e7c868` only.

**Budget:** at most 1,200 LPs (1 WT + 1 WT ceiling + about 1,082 knockouts + at most 20 fresh re-solves + 2 per
infeasible knockout, about 40). Expected run time under 30 minutes (inferred from the aconitase sweep).

## 5. Acceptance criteria

- A1. WT growth equals 0.3116608216 (aconitase summary) within 1e-9.
- A2. Rule A reproduces the aconitase candidate scores exactly at all four cutoffs, with the same 20 non-optimal
  knockouts.
- A3. Each knockout counted under rule B meets conditions 1–4; its maximum maintenance flux and margin are listed.
- A4. Only certified infeasible knockouts change class between rule A and rule B.

**Stop conditions** (stop, keep the outputs, report, no `benchmark.md` edit): WT not optimal; more than 20 fresh
re-solves; A1 or A2 fails. A knockout failing a certificate condition is not a stop: it stays unresolved and is
listed.

## 6. Outputs (this folder)

`scripts/`, `logs/`, `ko_screen.tsv` (per gene: raw status, growth, check magnitudes, certificate fields, class
under A and B), `ko_screen_summary.json`, `REVIEW.md` (script review), `REPORT.md`. Then, only if A1–A4 pass and
D1 is approved: the `benchmark.md` edit (and D2 files if approved) and a dated `model/STATE.md` entry that also
records the 2026-10-30 findings (subunit identities, COX9 BLAST and guide check, complex III not reviewed under the
2026-10-08 per-subunit rule). No commit or push.

## 7. Limits

- The certificate covers maintenance shortfalls only. An infeasibility from any other bound stays unresolved.
- Scores remain a development reference: the consensus list was used in development; negatives are the
  user-defined class; mtDNA genes are out of scope.
- The rule says what the model predicts. Whether complex III/IV loss is lethal in *Y. lipolytica* rests on the
  screens; the 5 negatives include 2 core subunits with borderline calls (CYT1, COX4), 2 accessory subunits (QCR6,
  QCR10) and COX9 (not measured).
- Complex III's GPR was never reviewed subunit by subunit the way complex IV was on 2026-10-08. That is a separate
  model question, not part of this round.
