# R1172 gene-less HMG-CoA transport (option A3): Stage 1 (in memory) and Stage 2 (closure built), 2026-11-12

Scope, authorization, evidence and plan: [TASK.md](TASK.md) (Stage 2 is section 8). Target state:
[before.json](before.json). Sections 1–6 are the Stage 1 results, kept as written (section 5a was the proposal the
user approved). Section 7 is Stage 2, with the independent audit in [AUDIT.md](AUDIT.md) and the figure in
[escher/](escher/README.md). Scripts, raw solves and sweeps:
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/r1172_hmgcoa_transport_20261112/` (`MANIFEST.sha256`). Stage 1 was committed
in `f379441c` and Stage 2 in `692f4c5f`; both are in PR #51.

## 1. Run

- Models: e584151f (`e584151f…39d2`) and the uncommitted R2121-closed candidate 2f352037 (`2f352037…af7c`), each
  through a short-path copy with the same SHA. Variants set in memory after the condition is applied:
  B0 = e584151f; V1 = B0 + R1172 (0, 0); V2 = V1 + R412 (0, 0); S0 = 2f352037; S1 = S0 + R1172 (0, 0).
- Task A medium (`6eca70fd…`) plus PO1f (`35307853…`). Gurobi 13.0.1, FeasibilityTol / OptimalityTol 1e-9, cobra
  0.30.0, Python 3.13.5 (root `.venv`). Platform files have the same SHAs as the v4 screen (`fc78ee8c`, `01b55f73`,
  `46865060`).
- Checks: `scripts/stage1_checks.py`. Sweeps: `scripts/ko_screen_certified_r1172.py`, the reference
  `ko_screen_certified.py` (`313d5d01…`) plus one line in `context()` that closes R1172 (diff in
  `scripts/ko_screen_certified_r1172.diff`), run on e584151f (V1) and on 2f352037 (S1) by `scripts/run_sweeps.zsh`.
  Comparison: `scripts/compare_screens.py` (generalised from the A2 `compare_v1_v4.py`).
- LPs: 223 for the checks (5 WT, 10 pFBA, 195 FVA, 5 energy, 8 controls) and 2,246 for the two sweeps (each 1,082
  sweep, 1 maintenance ceiling, 40 certificate; no fresh re-solve). 2,469 in total, against a budget of 2,600. No
  build. Each sweep ran 106–107 s.

## 2. Checks (observed this run)

| Check | B0 | V1 (R1172 closed) | V2 (+ R412 closed) | S0 (2f352037) | S1 (S0 + R1172 closed) | Stop condition |
| --- | --- | --- | --- | --- | --- | --- |
| WT growth | 0.31248856100769 (v4 + 3e-15) | 0.31184175100717 (−0.207 %) | 0.31184175100716 (= V1 − 2e-15) | 0.31205269538104 (= A2 V1) | 0.31110053592595 (−0.305 % vs S0) | not met (limit −5 %; baselines reproduce within 1e-9) |
| No-uptake ATP maximum | 0, optimal | 0, optimal | 0, optimal | 0, optimal | 0, optimal | not met |
| R411 (cytosolic HMGS), FVA 99.9 % | [0, 0.0348] | [0.0720, 0.0739] | [0.0720, 0.0739] | [0, 0.0236] | [0.0718, 0.0736] | not met (minimum > 0) |
| R82 (cytosolic thiolase), FVA | [0, 0.0348] | [0.0720, 0.0739] | [0.0720, 0.0739] | [0, 0.0236] | [0.0718, 0.0736] | not met (minimum > 0) |
| R1172, FVA | [−0.0741, −0.0373] | [0, 0] | [0, 0] | [−0.0740, −0.0485] | [0, 0] | — |
| R412 (mitochondrial HMGS), FVA | [0.0368, 0.1428] | [0, 0.0707] | [0, 0] | [0.0367, 0.1427] | [0, 0.0707] | — |
| R83 (mitochondrial thiolase), FVA | [0.0015, 0.0741] | [−0.0707, 0] | [−0.0707, 0] | [0.0014, 0.0740] | [−0.0707, 0] | — |
| HMGR R410, FVA | [0.0721, 0.0741] | [0.0720, 0.0739] | [0.0720, 0.0739] | [0.0720, 0.0740] | [0.0718, 0.0736] | — |
| Leucine route R1968 / R1973, FVA | [0, 0.0353] | [0, 0.0236] | [0, 0.0236] | [0, 0.0353] | [0, 0.0236] | — |
| 3-methylcrotonyl-CoA transport R2083 | [0, 0] | [0, 0] | [0, 0] | [0, 0] | [0, 0] | — |
| ATP-citrate lyase R1894, FVA | [0, 0.1232] | [0, 0.1232] | [0, 0.1232] | [0, 0.2802] | [0, 0.2826] | — |

B0 and S0 reproduce the v4 WT (0.3124885610076883, +3e-15) and the A2 V1 WT (0.3120526953810365, exactly). Every WT
point is optimal and passes the targeted check.

**Controls (gene knockouts and single reactions, KO/WT):**

| Knockout | B0 (bypass open) | V1 (R1172 closed) |
| --- | --- | --- |
| YALI1F38116g (R411 + R412) | 0.998 | 0 (optimal, growth 0) |
| YALI1B11532g (R82) | 1.000 | 0 (optimal, growth 0) |
| YALI1E13899g (R83) | 0.999 | 1.000 |
| R411 alone closed | — | 0 |
| R412 alone closed | — | 1.000 |

So in V1 the YALI1F38116g call comes from R411 alone; R412 does not contribute.

## 3. Mechanism

- **Baseline (observed, pFBA).** R83 (0.0722) makes mitochondrial acetoacetyl-CoA, R412 makes mitochondrial HMG-CoA,
  and R1172 carries it to the cytosol (−0.0722, mitochondrion → cytosol) for HMGR R410. The FVA minimum of R412 is
  0.037 and the maximum of R1172 is −0.037, so the baseline optimum cannot do without the mitochondrial route. R82
  and R411 can carry up to 0.035 but carry 0 in pFBA.
- **V1 (observed).** R1172 is closed; R82 → R411 → R410 carry the whole mevalonate supply (0.0720, required). R82
  takes 2 and R411 1 cytosolic acetyl-CoA per HMG-CoA, so cytosolic acetyl-CoA demand rises by about 0.216. In pFBA
  it comes from the carnitine shuttle: R191 cytosolic carnitine acetyltransferase 0.272 → 0.488, fed by the
  R2121-driven TAG lipolysis / peroxisomal β-oxidation cycle of the A2 report (R2121 0.186 → 0.327). ATP-citrate
  lyase R1894 stays at 0 in pFBA.
- **S1 (observed).** With R2121 already closed, the extra cytosolic acetyl-CoA comes from R2119 (malonate
  semialdehyde → acetyl-CoA + CO2 + NADH, YALI1C02550g) at 0.280 → 0.495, the Asp → β-alanine → malonate
  semialdehyde route the A2 report flagged (its §5 item 2). R1894 again stays at 0 in pFBA.
- **R412 once R1172 is closed (observed).** pFBA 0. FVA [0, 0.0707]; the same capacity appears on R1974 and R1933,
  which is the R412 + R1974 + R1933 (reverse) loop that sums to succinyl-CoA hydrolysis (stoichiometry, inferred).
  Closing R412 as well (V2) changes WT by −2e-15 and only narrows R1974 and R1933 to the leucine-route range
  [0, 0.0236].

## 4. Sweeps

Rule B; mtDNA and the current untargeted list (`YALI0F04114g`, `YALI0E16192g`, `244b941c…`) out of scope; model
genes not on the consensus list count as non-essential (a user-defined negative class); the list was used in
development, so these are a development reference. All four screens are rescored here under the same scope by
`compare_screens.py`. Both sweeps have WT equal to the check value, 20 certified-infeasible knockouts (the same 20
genes as v4 and A2 V1), no unresolved knockout and no fresh re-solve.

| Cutoff | v4 (B0) | V1 (R1172 closed) | A2 V1 (S0) | S1 (R2121 + R1172 closed) |
| --- | --- | --- | --- | --- |
| 1 % | 113 / 213 / 36 / 701 (0.3955) | 115 / 211 / 36 / 701 (0.4014) | 115 / 211 / 36 / 701 (0.4014) | 117 / 209 / 36 / 701 (0.4073) |
| 5 % | 120 / 206 / 38 / 699 (0.4103) | 122 / 204 / 38 / 699 (0.4161) | 122 / 204 / 38 / 699 (0.4161) | 124 / 202 / 38 / 699 (0.4219) |
| 10 % | 124 / 202 / 43 / 694 (0.4081) | 126 / 200 / 43 / 694 (0.4138) | 126 / 200 / 43 / 694 (0.4138) | 128 / 198 / 43 / 694 (0.4196) |
| 15 % | 141 / 185 / 43 / 694 (0.4560) | 143 / 183 / 43 / 694 (0.4616) | 143 / 183 / 43 / 694 (0.4616) | 145 / 181 / 43 / 694 (0.4670) |

TP / FN / FP / TN (MCC). With only the mtDNA scope, FP and TN are one higher in every column (`YALI0E16192g`); the
class changes are the same. At 15 %: recall 0.4325 → 0.4387 (V1) and 0.4448 (S1); specificity 0.9417 in every
column; precision 0.7663 → 0.7688 (V1) and 0.7713 (S1).

**Class changes, V1 vs v4 and S1 vs A2 V1 (the same 2 genes, both FN → TP at every cutoff):**
- `YALI1F38116g` (YALI0F30481g, Q6BZW0): no verified name; HMG-CoA synthase by annotation (EC 2.3.3.10). KO/WT
  0.998 → 0. Database annotation only.
- `YALI1B11532g` (YALI0B08536g, Q6CFB6): no verified name; acetyl-CoA C-acetyltransferase by annotation
  (EC 2.3.1.9). KO/WT 1.0 → 0. Database annotation only.

S1 vs v4 changes these two plus the two A2 genes (`YALI1E35955g`, `YALI1F28274g`): +4 TP, no other change. The two
closures do not interact in any class at any cutoff. No TP becomes FN and no knockout grows above WT + 1e-6 in
either sweep.

**Ratio shifts without a class change** (146 genes move by more than 1e-6 in V1, 135 in S1; most by about 0.002,
the WT change):
- `YALI1F20788g` (YALI0F15433g; no verified name; GPR of R834, coenzyme A transport cytosol → mitochondrion;
  consensus-positive, FN in both): 0.998 → 1.000. In B0, R1172 exports CoA as HMG-CoA from the matrix and R834
  brings free CoA back; with R1172 closed the matrix loses no CoA and R834 is no longer needed (observed: its v4
  KO/WT equals the V1 WT cost, 0.99793). It stays FN; its call now rests on the absence of any other net
  mitochondrial CoA demand.
- The 14-subunit ATP synthase AND group of R171 (e.g. `YALI1M00318g`, `YALI1F04433g`): 0.1095 → 0.1076, still
  between 10 % and 15 % (essential at 15 %, non-essential at 10 %, unchanged).
- In S1 the largest other shifts are under 0.003 (`YALI1F01981g` 0.917 → 0.920; `YALI1C08262g` 0.877 → 0.880).

## 5. What this means, and limits

- **Stop conditions:** none met (section 2 and the sweeps). The Stage 1 bar of TASK.md §7 is met: V1 WT −0.207 %,
  energy 0, R411 and R82 required, both genes at 0, no TP lost, every other shift explained.
- **What it means.** With R1172 closed, the cytosolic ERG10-like / ERG13-like route is the only source of mevalonate,
  so the model's sterol, UQ9 and dolichol requirement makes both genes essential, as in the consensus screens. The
  rescue in e584151f came from one gene-less transport whose existence has no evidence found for or against it.
  This is a mechanism result; it does not show that HMG-CoA cannot cross the inner membrane in *Y. lipolytica*.
- **R412.** Once R1172 is closed it has no role in growth or in any call (V2 = V1 to 2e-15; R412 alone closed
  1.000). Closing it is not needed for A3; its existence is a separate question (TASK.md §3.3).
- **Statement levels.**
  - Observed in this run: every number in sections 2–4.
  - Statically verified: the producer and consumer lists of the HMG-CoA, acetoacetyl-CoA and mevalonate pools, the
    external-model structures, the UniProt annotations.
  - Recalled: the Task C literature search, the a1237821 and a56f3649 closures, the A2 diff of 2f352037.
  - Inferred: that R412 + R1974 + R1933 is a dissipative loop; that a Stage 2 build would reproduce V1 or S1.
  - Untested: the built candidate; other media; whether *Y. lipolytica* has any mitochondrial HMG-CoA synthase or
    HMG-CoA carrier.
- **Limits.**
  - Single-gene knockouts, one medium, one objective. Route lists show one pFBA optimum; FVA covers 19 named
    reactions only.
  - No new literature search; the evidence for the closure is absence of found evidence plus S. cerevisiae
    localisation (TASK.md §3).
  - No independent audit of this Stage 1 (as in the A2 Stage 1); no Escher figure (it belongs to Stage 2).
  - S1 uses the uncommitted 2f352037; if that file changes, S1 must be re-run.
- **Possible improvements to the method.** A loopless FVA on R412 / R1974 / R1933 / R83 would confirm that the
  remaining R412 range in V1 is only the dissipative loop. Test: loopless vs standard FVA on those four reactions in
  V1; expected R412 maximum 0 under loopless.

## 5a. Stage 2 proposal (approved 2026-11-12; carried out in section 7)

1. **Decision only the user can make:** reverse the 2026-09-09 "Retained by user request" decision for R1172.
2. **Change:** R1172 bounds [−1000, 1000] → [0, 0], with the reaction, stoichiometry, empty GPR and existing notes
   kept and a dated note giving the evidence and this Stage 1 result (as for R2121), through a curated data file
   and a guarded opt-in builder step (skill `curation-change`). `retained_reactions.json` would get a dated note
   that the reaction is kept but closed. R412 unchanged.
3. **Base:** stack on the A2 closure (expected: WT 0.3111005359259494; 15 % rule B 145 / 181 / 43 / 694) once the
   A2 work is committed, because its builder step lives only in the uncommitted A2 worktree; otherwise build on
   e584151f (expected: WT 0.3118417510071664; 143 / 183 / 43 / 694).
4. **Acceptance criteria:** control build reproduces the base byte for byte; whole-model diff shows only R1172's
   bounds and the new notes; seed 0 and seed 1 builds byte-identical; WT equals the Stage 1 value within 1e-9 and
   the no-uptake ATP maximum is 0; the certified sweep reproduces `ko_screen_V1.tsv` (or `ko_screen_S1.tsv`) class
   by class with KO growth within 1e-9; pytest shows only new passing tests; review and independent audit leave no
   blocking finding; an Escher before/after figure of the HMG-CoA / acetyl-CoA routes.

## 6. Found, not changed (each affects another line)

1. **Cytosolic acetyl-CoA source.** Closing R1172 moves about 0.216 acetyl units into the cytosol, and in neither
   variant does ATP-citrate lyase take them. On e584151f they come through the R2121 TAG cycle and the carnitine
   shuttle; on 2f352037 through the Asp → β-alanine → malonate semialdehyde → R2119 route. Both routes were already
   flagged in the A2 report; this round makes them larger.
2. **R83's gene.** YALI1E13899g (YALI0E11099g) is UniProt Q6L8K7, gene name "PAT1", annotated peroxisome
   (subcellular location) and mitochondrion (GO IBA); the model places R83 in the mitochondrion only. Annotation
   only; not checked further.
3. **Yeast-GEM's ERG13 compartment.** Yeast-GEM 9.0.2 puts ERG13 (YML126C) only in the mitochondrion (`r_0560`) and
   needs the gene-less `r_1840` for sterols, while UniProt lists no experimental mitochondrial location for P54839.
   Relevant if Yeast-GEM is used as a structural reference.

4. **R2119 becomes required near the optimum (Stage 2 figure FVA).** In the R1172-closed candidate, R2119 (malonate
   semialdehyde → acetyl-CoA, YALI1C02550g) has an FVA range of [0.2115, 0.5647] at 99.9 % of the optimum (before:
   [0, 0.3502]); ATP-citrate lyase R1894 can carry up to 0.28 but cannot replace it within 0.1 % of growth. Whether
   the Asp → β-alanine → malonate semialdehyde route is a real cytosolic acetyl-CoA source in *Y. lipolytica* was not
   checked (item 1).
5. **HMGR R410 charge.** R410 has a charge imbalance of −2 (products minus substrates), the same in 2f352037 and the
   candidate; it is pinned and not changed here. Recalled context: the model has many charge-only imbalances from the
   protonation convention.

## 7. Stage 2: R1172 closed through the builder (2026-11-12)

**Authorization.** The user, 2026-11-12: "Reverse, approve Stage 2" (reversing the 2026-09-09 retention for R1172's
flux), then "Commit Stage 1, merge main (Recommended)" for the base, once A2 had been merged to main (PR #50). Stage 1
commit `f379441c`, merge of `origin/main` `b3c150f7`; no push. Scope: TASK.md section 8.

### 7.1 What changed and how

- **Data:** `model/curation/selection/r1172_hmgcoa_transport_close.json`, `R1172-HMGCOA-TRANSPORT-CLOSE-20261112`. It
  holds the reviewed R1172 state (name, empty GPR, bounds, stoichiometry, species, annotation, the three retention
  notes), the pinned cytosolic route R82 / R411 / R410, the evidence with tiers and sources, and the seven notes to
  add. Written from 2f352037 by `stage2/scripts/make_curation_json.py` (asserts the input SHA).
- **Step:** `platform/scripts/gem_annotate/r1172_hmgcoa_transport_close.py`, `apply_r1172_hmgcoa_transport_close`. It
  checks the spec, then R1172 and each pinned route reaction against the reviewed state, then the notes; any mismatch
  raises before an edit. It sets (0, 0) and adds the notes; a second run returns `already_correct`; a `with model:`
  context restores bounds and notes. It does not read or touch R2121 or R412.
- **Flag:** `--r1172-hmgcoa-transport-close` (opt-in), called in `main.py` after the R2121 closure, as the last model
  step.
- **Retention record:** `model/curation/selection/retained_reactions.json` gets `later_decision_20261112` on R1172,
  recording the reversal. `annotate_retained_reactions` writes only `previous_policy` and `current_status` into the
  model, so this field reaches the build record but not the XML (shown by the control build, criterion 1).
- **Test:** `platform/tests/test_r1172_hmgcoa_transport_close.py` (8 tests on an in-memory 2f352037): every guard
  fails closed without edits, malformed specs are rejected, only R1172's bounds and notes change and survive an SBML
  roundtrip, the step applies with R2121 open, the retention step leaves the closure alone, a model context restores the state, the flag is opt-in and
  forwarded by `build_model`.

### 7.2 Builds (miniconda Python 3.13.5, the 2f352037 argv; `stage2/build.zsh`, commands in `stage2/logs/`)

| Build | Flags | Output SHA-256 |
| --- | --- | --- |
| control | 2f352037 argv | `2f35203750d1904aa4e50dd7d8eac494d5e6a15b142ea005c33152cea889af7c` (= 2f352037) |
| candidate (seed 0) | + `--r1172-hmgcoa-transport-close` | `684f3752764fb8f25219008966953863141724f34a451dffb77c235aa3cb2f62` |
| seed 1 | + `--r1172-hmgcoa-transport-close` | `684f3752…` (identical) |

Candidate: `model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_bigg_annot_v2_mcnaugh_rxns_v3_r2121_close_r1172_close_20261112.xml`
(with its `.build.json` and `.coq9_genes.tsv`). Build log: `R2121 malonyl-CoA closure: applied`, `R1172 HMG-CoA
transport closure: applied`.

### 7.3 Acceptance criteria (TASK.md section 8)

| Criterion | Result |
| --- | --- |
| 1. Control reproduces 2f352037 | yes, byte for byte |
| 2. Whole-model diff = R1172 bounds + notes only | yes. 2 differences: bounds [−1000, 1000] → [0, 0]; 7 notes added, the 3 retention notes unchanged (`stage2/model_diff_candidate_vs_2f352037.json`; the seed-1 diff is the same). In the XML, R1172's `reversible` attribute becomes `false`; cobra derives it from the bounds. |
| 3. Seeds 0 and 1 byte-identical | yes |
| 4. WT and energy | WT 0.3111005359259494, identical to the Stage 1 S1 value. No-uptake ATP maximum 0, optimal (`stage2/built_checks.json`). |
| 5. Sweep reproduces S1 | yes. `stage2/ko_screen_built.tsv` is byte-identical to `ko_screen_S1.tsv`: 0 class differences, 0 KO differences above 1e-9; the same 20 certified-infeasible knockouts, no fresh re-solve. |
| 6. pytest | yes. Full suite, root `.venv`, `--continue-on-collection-errors`: before (git archive of `b3c150f7` at a short path) 469 passed / 3 failed / 18 errors / 24 skipped; after (working tree) 477 / 3 / 18 / 24. Test by test, only the 8 new tests differ, all passing; the 3 failures and 18 errors are the standing ones of the A2 runs (`stage2/pytest/compare_before_after.json`). A first before-run from the scratchpad was invalid: its path was too long for cobra (31 extra errors; kept as `before_longpath_invalid.*`). |
| 7. Review and audit | review: no blocking finding (2 low, 2 informational; 7.5). Audit: every claim SUPPORTED except the review's provenance (UNVERIFIED; recorded since); five paperwork problems, all addressed ([AUDIT.md](AUDIT.md)). |
| 8. Escher figure | yes: [escher/r1172_close_before_after.png](escher/r1172_close_before_after.png) and the offline viewer [escher/r1172_close_escher.html](escher/r1172_close_escher.html) ([README](escher/README.md)) |

- **Chemistry:** unchanged. R1172 is mass- and charge-balanced (the same species formula and charge on both sides)
  in both models; R82 and R411 are balanced; R410's −2 charge residual is inherited and unchanged (section 6 item 5).
- **Scores on the built candidate** (rule B, mtDNA and untargeted out of scope; model genes not on the consensus list
  count as non-essential, a user-defined negative class; development reference):

  | Cutoff | v4 (e584151f) | 2f352037 (A2) | R1172 closed (684f3752) |
  | --- | --- | --- | --- |
  | 1 % | 113 / 213 / 36 / 701 | 115 / 211 / 36 / 701 | 117 / 209 / 36 / 701 |
  | 5 % | 120 / 206 / 38 / 699 | 122 / 204 / 38 / 699 | 124 / 202 / 38 / 699 |
  | 10 % | 124 / 202 / 43 / 694 | 126 / 200 / 43 / 694 | 128 / 198 / 43 / 694 |
  | 15 % | 141 / 185 / 43 / 694 | 143 / 183 / 43 / 694 | 145 / 181 / 43 / 694 |

  Against 2f352037: `YALI1F38116g` and `YALI1B11532g` FN → TP at every cutoff, nothing else changes class
  (`stage2/compare_built_vs_2f352037.json`). Against v4: the same two plus the A2 pair (`stage2/compare_built_vs_v4.json`).
- **LPs:** 1,123 for the sweep, 2 for the checks (re-run twice after the audit's label fix: +4), 60 for the figure:
  1,189 against a budget of 1,400. 3 builds, no spare used. The auditor's 13 LPs are its own.

### 7.4 Escher figure

The map (24 reactions in three sections) shows the HMG-CoA for HMGR R410 coming from R83 → R412 → R1172 before
(with R834 returning the exported CoA to the matrix), and from R82 → R411 after. R410 (0.072) and mevalonate kinase
R553 (0.0545; the remaining 0.0175 goes through the CTP-dependent kinase R554, not on the map, observed in the
Stage 1 S0/S1 pFBA) carry flux in both. At 99.9 % FVA, R834 (CoA import into the
mitochondrion) is [0.0485, 0.074] before, the mirror of R1172's export, and [0, 0] after. Details:
[escher/README.md](escher/README.md).

### 7.5 Code review (fresh read-only reviewer, `.claude/agents/code-reviewer.md` run as a general-purpose agent)

No blocking finding. The reviewer also applied the step once in memory to e584151f (R2121 open): `applied`. Its report
is recorded in `stage2/code_review.md` (summary; the agent transcript is not kept).

| # | Finding | Action |
| --- | --- | --- |
| 1 (low) | The `retained_reactions` record in the `.build.json` is taken before the closure, so for R1172 it lists bounds [−1000, 1000] although the candidate has (0, 0) | Documented here: read `reaction_selection.r1172_hmgcoa_transport_close` for the closure; the retention record is a pre-closure snapshot. No code change. |
| 2 (low) | The stale-input test did not cover every pinned route field, and no test applied the step with R2121 open | Added negative cases for R82 annotation, R411 name and R410 bounds, and `test_independent_of_the_r2121_closure`; 8 tests pass (`stage2/pytest/focused_r1172_after_review.log`). |
| 3 (info) | `_check_spec` does not require note values to be strings | Not changed: the curated notes are strings, and the R2121 step has the same shape; tightening both belongs together. |
| 4 (info) | The XML also shows `reversible="false"` on R1172 | Derived from the bounds; noted under criterion 2. |

### 7.6 Verification levels and limits

- **Observed this run:** the builds and SHAs, the whole-model diff, the balance check, WT, the energy check, the
  sweep and its comparison, the figure fluxes, the viewer checks.
- **Statically verified:** the curated spec against the 2f352037 file; the retention step's fields.
- **Recalled:** the Task C literature search; the protonation-convention context of R410.
- **Untested:** other media and objectives; a loopless check of the R412 / R1974 / R1933 range (section 5).
- **Not changed, still open:** whether *Y. lipolytica* has a mitochondrial HMG-CoA synthase (R412) or any HMG-CoA
  carrier. Under root AGENTS.md §5 (inactive reactions, added on main 2026-11-12), R1172 at (0, 0) is itself an
  inactive reaction: it is recorded as an open gap with no known carrier, not as a finished fix of the transport.
