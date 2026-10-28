# TASK: hidden-source check (acyl, phosphate, inositol) and R2171, Stage 1 (2026-10-28)

**Authorization.** The user, 2026-10-28, after reading this draft (shared notes folder,
`lipid_leak_check_r2171_TASK_draft_20261028.md`): "Approve Stage 1", meaning the sentence in section 6.
Stage 2 (any R2171 change) is not authorized. Sections 0–6 below are the draft unchanged apart from this
paragraph and the status line.

**Status: Stage 1 completed 2026-10-28 on branch `claude/lipid-source-leak-check` from `origin/main` `0093c89`
(base re-hashed before the run: unchanged). It stopped on the stop condition "a leak is found": R2171 is an
inositol source ([REPORT.md](REPORT.md)). 279 of 300 LPs. Stage 2 not authorized.** The draft's facts came from static reads and database pages
opened on 2026-10-28; no LP was solved for the draft.

## 0. Why

- R2135 (acetyl-CoA counted as acyl-CoA) and R2145 (0.01 coefficients, 100 Pi per PA) were hidden sources:
  they made acyl chains and phosphate from nothing and were found only because a knockout changed class.
  The model's existing guard, the no-uptake ATP / biomass check, did not see either.
- The same inherited block (R2130–R2200, iYli21) still has scaled reactions. **R2171** is the riskiest:
  `PI + 100 H2O -> DAG + 100 inositol-1-phosphate + 100 H+`.
- This task asks one question per source: **with the source removed, can the model still grow or still make
  the species?** It is read-only. A fix of R2171 is a separate Stage 2 that needs its own approval.

## 1. Starting state (read 2026-10-28)

- **Base model:** `9d655c2b` =
  `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_20261028.xml`
  on `origin/main` `0093c89` (PR #29), SHA-256 `9d655c2bdbf82721a38f148a9a4ee0a385bb509ea7b591111107113c4aa8335c`
  (re-hashed from `git show origin/main:…`). Not a designated reference.
- **Moving main:** the UM7 nuclear-transport session is running and may add a newer candidate. If `main`
  moves before Stage 1 starts, the base is re-pinned and stated, not silently replaced.
- **Branch:** this worktree's branch `claude/kennedy-closures-864c91` is merged (PR #26, #27). Stage 1 would
  run on a new branch from `origin/main`; creating it is part of the authorization asked for below.
- **Solver:** Gurobi 13 at 1e-9 through `model/reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py`
  (reset per solve, bound / residual / knocked-flux checks). `/opt/miniconda3/bin/python`, cobra 0.30.

## 2. R2171 as it is (verified, static)

| Field | Value |
| --- | --- |
| Equation | `m1393[C_cy]` PI + 100 H2O → 100 H+ + `m1994[C_cy]` DAG + 100 `m813[C_cy]` 1D-myo-inositol 1-phosphate |
| Bounds / GPR | (0, 1000) / `YALI1D13510g` (= YALI0D10813g, locus map S2) |
| Annotations | KEGG R03332; EC 4.6.1.13, 3.1.4.11, 3.1.4.3; MNXR148100; no notes |
| Provenance | identical (100 / 1 / 100 / 100 / 1) in `model/source/iyali26.xml` and `model/source/iyli21.xml` |
| Balance | charge −100 (100 × Ins1P −2, 100 H+, PI 0, DAG 0); mass not checkable (PI `C11H17O13P*2`, DAG no formula) |

**Reference forms (opened 2026-10-28):**
- [KEGG R03332](https://www.kegg.jp/entry/R03332): `PI + H2O <=> Inositol 1-phosphate + DAG`, 1:1:1:1. Verified.
- [ModelSEED rxn06362](https://modelseed.org/solr/reactions/select?wt=json&q=id:rxn06362) and rxn39389: the same
  1:1 hydrolase, status OK. No ModelSEED reaction found with 100× coefficients (name search
  "inositolphosphohydrolase", 8 hits). Verified.
- Gene: [KEGG yli:2910441](https://www.kegg.jp/entry/yli:2910441) (YALI0D10813g) = RefSeq
  "1-phosphatidylinositol phosphodiesterase", KO K01771, EC 4.6.1.13. Verified (orthology only).

**Inferred, not checked:**
- The 100× looks like R2145's 0.01 lumped scheme multiplied by 100 so that PI reads 1; with `m1393` and `m1994`
  linked 1:1 elsewhere, each PI then yields 100 Ins1P. **Each Ins1P carries one phosphate and six carbons,**
  so R2171 is a potential phosphate and carbon source, and an inositol source through IMPase R558.
- EC 4.6.1.13 enzymes form inositol 1,2-cyclic phosphate rather than Ins1P (established / uncited). Which
  product is right matters only for Stage 2.
- A loop could close since R1726 was opened: PI → R2171 → DAG `m1994` → R2169 → R2168 → TAG → R1771 backwards →
  DAG → R1726 → PC → PLD → PA → CDP-DAG → PI. On `e9d35857` (R1726 closed) R2171's maximum was 2e-6 at 99.9 %
  (R2135 TASK §11); it has **not** been checked on `ce9e1112` or `9d655c2b`.
- The Task A medium supplies myo-inositol (R1263 uptake 1000) and phosphate (R1296 uptake 1000), so a leak
  of either would not change Task A growth. That is why growth checks on Task A cannot see it.

## 3. Sources and what counts as "removed" (static inventory of `9d655c2b`)

| Source | Removed by | Exchanges closed to uptake in the open condition |
| --- | --- | --- |
| Acyl chains | FAS knockout `YALI1B19844g` (AND partner `YALI1B25427g`; both genes are on the same 12 reactions: R1392, R1393 and ten in R2217–R2241; which ones close is recorded in the run) | 9 acyl-carrying exchanges: R951, R1083, R1236, R1281, R1288, R1289, R1329, R1440, R1444 |
| Phosphate | no phosphate uptake | 11 P-containing exchanges: R950, R951, R1005, R1014, R1102, R1126, R1269, R1296, R1322, R1339, R1341 |
| Inositol | INO1 knockout `YALI1B05798g` (R559, G6P → Ins1P) | 2 inositol-containing exchanges: R950, R1263 |
| (all, open condition only) | — | 15 formula-less exchanges (R947, R1067, R1189, R1294, R1984, R1986, R1990, R1992, R1998, R2018, R2026, R2031, R2036, R2214, R2275), because their contents cannot be traced |

## 4. Stage 1 computation (needs approval)

**Models:** `9d655c2b` (base). One in-memory sensitivity variant: the same with DAG kinase R1720 (0, 0) → (0, 1000),
because R1720 is the next Kennedy-line fix and could close new loops. No other model.

**Conditions** (PO1f profile throughout):
- C1 Task A (`po1f_csm_leu_w29uptake_20261012.csv`, `6eca70fd`).
- C2 Task A + choline R1041 and ethanolamine R1115 uptake 1 (the Kennedy mechanism control).
- C3 glycerol medium (`po1f_csm_leu_w29glycerol_20261012.csv`, `23371f0a`).
- C4 oleate medium (`po1f_csm_leu_w29oleate_20261012.csv`, `707f15d8`): **positive control** for the acyl check;
  the FAS knockout must grow here, or the test cannot detect growth.
- C5 open condition (every exchange (−1000, 1000), as in the gap census), minus the source's exchanges and the
  15 formula-less exchanges.

**Per source × condition:**
1. WT growth, checked.
2. Source removed (section 3): growth, checked. Expected 0 except the acyl source in C4.
3. Source removed, biomass lower bound 0, in-memory demand reactions (test only, not model edits) maximised one
   at a time:
   - acyl: `m567`, `m570`, `m790`, `m1984`, palmitoyl-CoA `m1349[C_cy]`;
   - phosphate: Pi `m35[C_cy]`;
   - inositol: myo-inositol `m814[C_cy]`, Ins1P `m813[C_cy]`.
   Expected 0 (tolerance 1e-9 × scale; values ≤ 1e-7 reported as numerical zero, per the FVA-near-1e-7 caveat).
4. R2171: FVA at 99.9 % and at growth ≥ 0 in C1–C5 with all sources present.
5. Any non-zero in 2 or 3: one pFBA trace plus single-closure tests of R2171, R2165, R1726, R1842, R2168 (the
   STOP.md method). Report the route; change nothing.

**Budget:** about 3 sources × 5 conditions × (2 growth + up to 5 demand) ≈ 100 LPs, R2171 FVA 20, the R1720
variant on C1 and C5 about 30, traces ≤ 100. **Cap 300 LPs.**

**Outputs:** `model/reports/lipid_source_leak_check_20261028/` (TASK.md, REPORT.md, `results.json`, scripts) and
large traces in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/lipid_source_leak_check_20261028/`. A dated
`model/STATE.md` entry. An Escher figure of any leak route found (standing request).

**Stop conditions:**
- Any WT solve not optimal or failing the check: stop, record, no retry with other settings.
- C4 positive control fails (FAS knockout does not grow on oleate): stop; the acyl test is not valid.
- A leak is found: finish the trace for that leak only, report, stop. No fix in Stage 1.
- 300 LPs reached, or `main` changes the base: stop and report.

**Not in Stage 1:** any edit to model, curation data, GPR, bounds or tests; commit or push (asked separately);
the IPC / MIPC synthases (0.01), R2165 and R2187 beyond appearing in traces; the lipid-unlump species (colleague).

## 5. Stage 2 (outline only; separate approval after Stage 1)

If Stage 1 shows R2171 is a source, or regardless of that on the chemistry alone (user's call):
- **A. 1:1 hydrolase** as KEGG R03332 / ModelSEED rxn06362: `PI + H2O -> DAG + Ins1P + n H+`, with `n` chosen to
  balance charge under the model's existing charges (PI 0, Ins1P −2 → n = 2). No charge edits (they ripple).
- **B. Cyclic-phosphate lyase** per the gene's EC 4.6.1.13: needs a new metabolite (inositol 1,2-cyclic phosphate)
  and its hydrolase. More invasive; only with evidence.
- **C. Close R2171** with the reason recorded.
- Route: curation-change skill (curated JSON, guarded `apply_*` step after `--dolichol-route`, focused test, new
  candidate, whole-model diff, knockout sweep at 1/5/10/15 %, re-census, review, audit, REPORT, STATE, Escher).

## 6. Authorization sentence (for Stage 1)

> Approved: Stage 1 of the lipid source leak check on 9d655c2b: acyl / phosphate / inositol source-removal tests and
> demand maxima in C1–C5, R2171 FVA, the in-memory R1720 variant, traces of any leak; ≤ 300 LPs; new branch from
> origin/main; outputs, STATE entry and Escher figure; no model change, no commit or push.
