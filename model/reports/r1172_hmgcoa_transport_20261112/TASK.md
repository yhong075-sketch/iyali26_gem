# R1172 gene-less HMG-CoA transport (option A3): Stage 1 (in memory, no build), 2026-11-12

## 0. Authorization

The user, 2026-11-12, in the session started from the reaction-activation-options line for option A3:
"Approve Stage 1", approving the Stage 1 plan as written in that session (reproduced in sections 4–7), with both
sweeps (V1 and the stacked S1).

Scope: computation and documentation only. No model XML, curated data, builder code or test changes; no build; no
commit or push. Stage 2 (a curation change) needs a separate approval. Approving Stage 1 does **not** reverse the
recorded 2026-09-09 decision to retain R1172 (section 3.2).

## 1. Starting state (checked 2026-11-12)

- Branch `claude/trusting-mahavira-6f0342`, HEAD `154eb64ca58f827863bdf061e5b4d03b4db05f4b`, clean.
- Model executed: `model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_bigg_annot_v2_mcnaugh_rxns_v3_20261103.xml`,
  SHA-256 `e584151fc79475fc2daa59e76c8c84ebfd06b0f3cb1fbd75fbb69cc9f31039d2` (committed in `ad927509`, on
  `origin/main`), read through a short-path copy with the same SHA (cobra rejects paths of 260 characters or more on
  macOS).
- Stacked model (read only): the R2121-closed candidate
  `…_mcnaugh_rxns_v3_r2121_close_20261112.xml`, SHA-256
  `2f35203750d1904aa4e50dd7d8eac494d5e6a15b142ea005c33152cea889af7c`, uncommitted in the worktree
  `yali1-c06006g-gene-024dc8` (branch `claude/reaction-activation-options-1ca120`), read through a short-path copy
  with the same SHA. Its recorded diff against e584151f is R2121 bounds and notes only
  (`$IYALI26_RESEARCH_ROOT/artifacts/tasks/r2121_malonyl_coa_20261110/stage2/model_diff_candidate_vs_e584151f.json`).
- Reference model: none designated (`model/candidates/README.md`).
- Condition: Task A medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) plus the PO1f
  profile `$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json` (`35307853…`), applied in memory.
- Baseline screen: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/mcnaugh_reactions_review_20261103/screen_v4_20261104/ko_screen_v4.tsv`
  (`90865ce0…`; certified, Gurobi 1e-9; WT 0.3124885610076883). Stacked baseline screen: the A2 V1 sweep
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r2121_malonyl_coa_20261110/ko_screen_v1.tsv` (WT 0.3120526953810365;
  byte-identical to the sweep of the built 2f352037 per the A2 report).
- Target state: `before.json` in this folder.

## 2. Question

R1172 (`m646[C_cy] <=> m648[C_mi]`, 3-hydroxy-3-methylglutaryl-CoA, no GPR, bounds [−1000, 1000]) carries the
cytosolic HMG-CoA that HMGR R410 turns into mevalonate in the WT optimum of the map data (recalled: R83 → R412 →
R1172 mi→cy 0.0722 → R410). The cytosolic route ERG10-like R82 (YALI1B11532g) + ERG13-like R411 (YALI1F38116g)
carries 0. Both genes are consensus-essential and predicted non-essential in v4 (KO/WT 1.0 and 0.998). Stage 1 asks
what the model does when R1172 is closed, and whether R412 matters once it is.

## 3. Evidence

### 3.1 For closing R1172

| Claim | Source | Tier |
| --- | --- | --- |
| R1172 has no GPR and no carrier; the model has a cytosolic HMG-CoA synthase R411, so R1172 is not needed for connectivity | static read of e584151f, 2026-11-12 | verified |
| No evidence found for HMG-CoA crossing the inner mitochondrial membrane in yeast | `model/reports/w29_13c_consistency_20261012/REPORT.md` §5.4 (subagent search) | recalled; absence of found evidence, not evidence of absence |
| S. cerevisiae ERG10 (P41338) is cytosolic by experiment (ECO:0000269, PubMed 5571830; GO cytoplasm IDA:SGD) | [UniProt P41338](https://www.uniprot.org/uniprotkb/P41338/entry) | verified (UniProt REST, 2026-11-12); the cited paper not opened; related yeast, not Y. lipolytica |
| S. cerevisiae ERG13 (P54839) has no subcellular-location comment; GO nucleus HDA:SGD and mitochondrion NAS:UniProt only | [UniProt P54839](https://www.uniprot.org/uniprotkb/P54839/entry) | verified (UniProt REST, 2026-11-12) |
| Y. lipolytica orthologs: Q6BZW0 (YALI0F30481g = YALI1F38116g) has no compartment annotation; Q6CFB6 (YALI0B08536g = YALI1B11532g) cytoplasm IEA:EnsemblFungi | [Q6BZW0](https://www.uniprot.org/uniprotkb/Q6BZW0/entry), [Q6CFB6](https://www.uniprot.org/uniprotkb/Q6CFB6/entry) | verified, annotation only |
| The inner mitochondrial membrane does not let CoA thioesters through without a carrier | the two URLs in the Patch 8c comment (`platform/scripts/gem_annotate/patches.py`) are about fatty acyl-CoA and were not opened | established / uncited; not counted as HMG-CoA evidence |

The consensus-essential label of the two genes is what is being predicted; it is not evidence for the change.

### 3.2 Against closing, and the recorded decision

- **Recorded user decision (verified):** `model/curation/selection/retained_reactions.json` and the e584151f notes:
  "Retained by user request; native transport mechanism and carrier GPR remain unconfirmed" (user request
  2026-09-09; the earlier policy had removed R1172). The Patch 8c removal list was emptied for it. Only the user can
  reverse it.
- **Other models (verified, static scan of local copies 2026-11-12):** Yeast-GEM 9.0.2 `r_1840`, iMM904 `HMGCOAtm`
  and iYali5 `1840` carry the same gene-less reversible cytosol↔mitochondrion HMG-CoA transport. Yeast-GEM puts
  ERG13 (YML126C) only in the mitochondrion (`r_0560`) and so needs `r_1840` for sterols; iYali5 has both a
  cytosolic and a mitochondrial HMG-CoA synthase, as iYali26 does. A model lineage, not biological evidence.
- An uncharacterised carrier could exist (inferred).

### 3.3 R412

With R1172 closed, mitochondrial HMG-CoA can leave only through the lyase R1974 (static, verified). R412 + R1974 +
R1933 (reverse) sum to succinyl-CoA + H2O → succinate + CoA (inferred from stoichiometry; dissipative). R412 shares
its gene with R411, so it cannot change YALI1F38116g's call once R411 is the only cytosolic source (inferred). V2
tests this in memory only. No R412 change is proposed in this round.

## 4. Genes

No name below is experimentally verified for Y. lipolytica. All three pass the W29 tag-protein audit (`ok` in
`model/curation/annotation/w29_tag_protein_audit.csv`).

| Gene | YALI0 / UniProt | Function (database) | Reactions | v4 KO/WT, class (15 %, rule B) |
| --- | --- | --- | --- | --- |
| YALI1F38116g | YALI0F30481g / Q6BZW0 | HMG-CoA synthase, EC 2.3.3.10 (annotation) | R411 (cy), R412 (mi) | 0.998, FN |
| YALI1B11532g | YALI0B08536g / Q6CFB6 | acetyl-CoA C-acetyltransferase, EC 2.3.1.9 (annotation) | R82 (cy) | 1.0, FN |
| YALI1E13899g | YALI0E11099g / Q6L8K7 ("PAT1" in UniProt; peroxisome SL, mitochondrion IBA) | acetyl-CoA C-acetyltransferase (annotation) | R83 (mi) | 0.999, TN |

## 5. Stage 1 computation

Gurobi 13, FeasibilityTol / OptimalityTol 1e-9, the checked solve of
`model/reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py` (`69c42690…`), and the certified sweep of
`model/reports/ceramide_pool_biomass_20261101/scripts/ko_screen_certified.py` (`313d5d01…`), with one added line
that applies the variant in `context()`. Root `.venv` (Python 3.13.5, cobra 0.30.0).

Variants, set in memory after the condition is applied:

| Variant | Model | In memory |
| --- | --- | --- |
| B0 | e584151f | — |
| V1 | e584151f | R1172 (0, 0) |
| V2 | e584151f | R1172 (0, 0), R412 (0, 0) |
| S0 | 2f352037 | — |
| S1 | 2f352037 | R1172 (0, 0) |

Per variant:
1. WT growth (checked solve).
2. pFBA: producers and consumers of HMG-CoA (cy, mi), acetoacetyl-CoA (cy, mi), cytosolic acetyl-CoA and NADPH.
3. FVA at 99.9 % of the optimum on 19 named reactions: R1172, R412, R411, R82, R83, R410, R1973, R1974, R1933, R1996,
   R2083, R1968, R90, R1894, R191, R192, R1032, R88, R2121.
4. Energy check: every exchange closed to uptake, positive lower bounds relaxed to 0, maximise `R_NTP1`.

Controls on V1: knockouts of YALI1F38116g, YALI1B11532g and YALI1E13899g; R411 alone closed; R412 alone closed.

Certified single-gene sweeps on V1 and S1 (1,081 knockouts each), scored at 1, 5, 10 and 15 % under rules A and B
with the mtDNA and current untargeted scopes, compared gene by gene with v4 (V1) and with the A2 V1 sweep (S1).

Budget: at most 2,600 LPs; no builds.

## 6. Stop conditions

Stop and report, without continuing to the sweeps or proposing Stage 2, if:
- B0 does not reproduce 0.3124885610076883 or S0 does not reproduce 0.3120526953810365 within 1e-9;
- any solve of steps 1–4 is non-optimal;
- any energy check is above 0;
- V1 WT growth drops more than 5 %;
- another route takes over cytosolic HMG-CoA in V1 (R411 or R82 FVA minimum 0);
- in a sweep, any TP becomes FN, or any knockout grows more than 1e-6 above WT;
- the budget would be exceeded.

## 7. Outputs and what Stage 1 must show before Stage 2 is proposed

- This TASK.md and `before.json` (repo); scripts, solves, sweep TSVs, scores and a manifest in
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r1172_hmgcoa_transport_20261112/`; REPORT.md and a STATE.md entry at the
  end of the round.
- To propose Stage 2: V1 WT within 5 %, energy 0, R411 and R82 required; YALI1F38116g and YALI1B11532g below 1 %;
  no TP lost; every other class change explained. Even then, Stage 2 needs the user to reverse the 2026-09-09
  retention decision explicitly.

Scoring note: model genes not on the consensus list count as non-essential (user-defined negative class); mtDNA and
untargeted genes are out of scope; the list was used in development, so scores are a development reference.

## 8. Stage 2: close R1172 through the builder (authorized 2026-11-12)

**Authorization.** The user, 2026-11-12, after reading the Stage 1 results and the Stage 2 proposal (REPORT.md §5a):
"Reverse, approve Stage 2", answering "Reverse the 2026-09-09 R1172 retention and approve Stage 2?". This reverses
the 2026-09-09 retention decision for R1172's flux: the reaction stays in the model but is closed. Then, asked how to
reach the A2 base (merged to main as PR #50 after the proposal was written): "Commit Stage 1, merge main
(Recommended)" = commit the Stage 1 docs on this branch, merge origin/main into it, build on 2f352037 with both
closure flags, no push. Done: Stage 1 commit `f379441c`, merge `b3c150f7`.

**Baseline.** 2f352037 (`model/candidates/E5_…_mcnaugh_rxns_v3_r2121_close_20261112.xml`, SHA-256
`2f35203750d1904aa4e50dd7d8eac494d5e6a15b142ea005c33152cea889af7c`, committed `1d0ed68b`), built by the argv in its
`.build.json` (HEAD `154eb64c` then, plus the A2 changes now committed). HEAD now `b3c150f7`.

**Targets.**
- Changed: R1172 bounds [−1000, 1000] → [0, 0], plus notes. Reaction, stoichiometry, empty GPR, annotation and
  existing notes kept.
- Pinned, not changed: the cytosolic route that takes over, R82 (YALI1B11532g), R411 (YALI1F38116g) and HMGR R410
  (YALI1E05731g). R412 is not touched.
- Curated record: `model/curation/selection/retained_reactions.json` gets a dated field on R1172 recording the
  reversal; the builder reads only `previous_policy` and `current_status`, so no XML changes from it.
- Nothing else may differ.

**Implementation.**
- Data: `model/curation/selection/r1172_hmgcoa_transport_close.json` (evidence fields: source URLs, tiers,
  reviewer).
- Step: `platform/scripts/gem_annotate/r1172_hmgcoa_transport_close.py`, `apply_r1172_hmgcoa_transport_close`;
  guarded, idempotent, fails closed when R1172 or the pinned route differs from the reviewed state.
- Flag: `--r1172-hmgcoa-transport-close`, opt-in, run as the last model step (after `--r2121-malonyl-coa-close`; it
  writes bounds and notes only and does not depend on R2121).
- Test: `platform/tests/test_r1172_hmgcoa_transport_close.py` on an in-memory 2f352037.

**Expected effect.** The Stage 1 S1 problem: WT 0.3111005359259494; at 15 % (rule B, mtDNA + untargeted out of
scope) 145 / 181 / 43 / 694.

**Acceptance criteria.**
1. A control build without the new flag reproduces 2f352037 byte for byte.
2. The whole-model diff of the flag build against 2f352037 shows only R1172's bounds and the new notes.
3. Seed 0 and seed 1 builds are byte-identical.
4. WT growth on the built candidate equals 0.3111005359259494 within 1e-9; no-uptake ATP maximum 0 (optimal).
5. The certified sweep on the built candidate reproduces `ko_screen_S1.tsv` class by class at every cutoff, with KO
   growth within 1e-9.
6. pytest: the new tests pass and no existing test changes outcome.
7. Code review and independent audit leave no unresolved blocking finding.
8. An Escher before/after figure of the HMG-CoA, acetoacetyl-CoA and cytosolic acetyl-CoA routes.

**Budget.** 3 builds (control, candidate, seed-1 rebuild) plus 1 spare; at most 1,400 LPs (sweep about 1,123, WT and
energy checks, pFBA for the figure). Two full pytest runs (before on a git-archive copy of `b3c150f7`, after on the
working tree).

**Stop conditions.** Any acceptance criterion fails, or the diff shows anything outside the targets: stop and
report, without retrying or relaxing the criterion.
