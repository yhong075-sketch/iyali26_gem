# GPI co-products: implementation report

Scope, authorization, budget and acceptance criteria are in [TASK.md](TASK.md). Exact commands are in
[commands.txt](commands.txt). This change builds on
[gpi_chain_links_20261021](../gpi_chain_links_20261021/REPORT.md).

Evidence tiers used below:
- **verified**: run or checked in this session; the output is cited.
- **provable**: follows from verified facts.
- **inferred**: reasoned, with no direct check.
- **established/uncited**: no source was opened this session.

## 1. Summary

- **New candidate (verified):**
  - `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_coproducts_20261021.xml`, SHA-256
    `cf0825b86c3533d8a12f697a994eba3d740671a68a367143ee32cbfdc8686b7c`.
  - It is `164f235f` plus `--gpi-chain-links --gpi-coproducts`. It is a candidate, not a reference,
    and it is uncommitted.
- **What changed against `89559a0f` (verified):** 5 differences.
  - R373 now also produces DAG (m604[C_er]).
  - R375 now also produces Dol-P (m503[C_er]).
  - Both carry 7 new `gpi_coproducts_*` notes.
  - R_GPIGNT's `gpi_chain_links_rationale` note is reworded.

  Nothing else changed. The model still has 2,311 reactions, 1,876 metabolites and 1,074 genes.
- **Effect (verified, SD-Leu + PO1f, GLPK):**
  - WT growth is unchanged: 1.4308057353528538 before and 1.4308057353528474 after.
  - The chain is cheaper. The FVA maximum rose from 0.497 to 0.750.
  - Forcing 1e-5 GPI now costs 1.88e-5 growth instead of 2.83e-5, and R_DAGt_er_em carries 2e-5,
    twice the chain flux.
- **The open exit is not neutral in every knockout (found by the audit, verified):** in nine
  knockouts at KO/WT 0.1277, the chain turns surplus PE into DAG through the open exit, for
  +8.477e-7 growth (§6).
- **Knockouts (verified):** no essentiality class change.
  - The only differing sweep row is YALI1D22124g. The `89559a0f` sweep's GLPK artifact reproduced
    there; on the new model it is optimal at 0.
  - Fresh solves give about 0 on both models.
- **Tests (verified):** the only difference from the previous run is the 4 new passing tests.
- **Audit (§7):**
  - It supports or partially supports every claim.
  - Its guard gaps are fixed; the rebuild still gives `cf0825b8`.
  - Its corrections to my R376 and DAG findings are applied (§8).

## 2. The change

| Reaction | Before | After |
| --- | --- | --- |
| R373, GPI step 6 (EtNP from PE to Man1); GPR YALI1F27343g | m598 + PE m599 ⇌ m600 | m598 + PE m599 ⇌ m600 + **DAG m604[C_er]** |
| R375, GPI step 8 (third mannose); GPR YALI1D05270g | Dol-P-Man m502 + m601 ⇌ m602 | Dol-P-Man m502 + m601 ⇌ m602 + **Dol-P m503[C_er]** |

Bounds (±1000), GPRs, names and annotations are unchanged.

| Item | Location |
| --- | --- |
| Curated data | [`model/curation/chemistry/gpi_coproducts.json`](../../curation/chemistry/gpi_coproducts.json) (`GPI-COPRODUCTS-20261021`). SHA-256 `bc6b29ed…` after the audit fix; the candidate's build.json records `825e429b…` |
| Note rewording | [`model/curation/gapfill/gpi_chain_links.json`](../../curation/gapfill/gpi_chain_links.json), R_GPIGNT `rationale` (SHA-256 now `33209268…`) |
| Guarded step | `platform/scripts/gem_annotate/gpi_coproducts.py`, `apply_gpi_coproducts` |
| Wiring | `cli.py` flag `--gpi-coproducts`. `main.py` refuses that flag without `--gpi-chain-links` and runs the step right after the links |
| Focused test | `platform/tests/test_gpi_coproducts.py` |
| READMEs | `model/curation/README.md`; `model/candidates/README.md` |
| Code diff | [logs/code_diff_cumulative_vs_a78a49d.patch](logs/code_diff_cumulative_vs_a78a49d.patch), cumulative with the links change |

**What the step checks before writing.** It fails closed and makes no partial edits.
- The scope: schema, ID, status, exactly R373 and R375, and exactly the five links in
  `requires_reactions`.
- That all five links are present with their recorded bounds and stoichiometry. Before the audit
  fix, this was a presence check only.
- The identity of DAG m604, Dol-P m503, Dol-P-Man m502 and PE m599.
- That m604's reactions, other than R373, are exactly R376 and R_DAGt_er_em.
- For each target: name, bounds, GPR and annotation equal the record, and the stoichiometry is
  exactly the before-state or the after-state.
- No notes outside this curation, and no conflicting curation notes.
- No after-state without this curation's notes, and no before-state carrying them (added after the
  audit). Not just one target done.

The edit is applied as the full before → after difference, so the result equals the recorded
after-state (changed after the audit).

A rerun returns `already_correct`.

**The note rewording.**
- **Before:** "The 2 H+ reconcile the stored charges (UDP-GlcNAc -2, UDP 0); KEGG does not balance
  protons."
- **After:** "The 2 H+ are needed because the model stores UDP-GlcNAc (m403) at charge -2 and UDP
  (m11) at 0, as R778 and R211 do; KEGG writes the equation for neutral species."

This answers finding 6 of the links audit.

## 3. Mechanism and evidence

**R373, DAG.**
- In KEGG, ethanolamine phosphate transfer from PE releases 1,2-diacyl-sn-glycerol.
  - The step closest to R373 is [R12351](https://www.kegg.jp/entry/R12351), by GPI-EtNP transferase
    1, [K05285](https://www.kegg.jp/entry/K05285).
  - The same holds for all five PE-donor EtNP transfers in map00563 (R05923, R05924, R08107,
    R12351, R12702).
  - **Verified**, read on 2026-10-20 and 2026-10-21.
- KEGG's R12351 acts on a two-mannose intermediate; the model's R373 acts on the one-mannose m598.
  So R12351 is cited for the chemistry (PE in, DAG out), not as the same reaction. No cross-reference
  was added.
- **Gene identity (verified, 2026-10-21):** YALI1F27343g's identifiers in the model are RefSeq
  XP_505669.2, UniProt Q6C0Z3 and KEGG yli:2908691. All three point to the same protein:
  - UniProt Q6C0Z3 is MCD4_YARLI, "GPI ethanolamine phosphate transferase 1", ordered locus
    YALI0F20570g, a reviewed entry.
  - KEGG assigns yli:2908691 to K05285.

  KEGG lists the gene position as chromosome E while the model ID is on F; the identifiers agree, so
  this is a position label, not a different gene. This is database identity, not experimental
  validation in *Y. lipolytica*.
- **Mass balance:** not checkable, before or after.
  - The GPI intermediates m598 and m600 have no formula.
  - The model's DAG m604 is C6H8O6, which is not a diacylglycerol (§8).
  - cobra's residual changes from C−7 H−12 N−1 O−8 P−1 R−2 to C−1 H−4 N−1 O−2 P−1 R−2. That only
    reflects PE minus the stored "DAG" formula ([logs/balance.txt](logs/balance.txt), verified).

**R375, Dol-P.**
- KEGG [R05922](https://www.kegg.jp/entry/R05922) (GPI mannosyltransferase 3,
  [K05286](https://www.kegg.jp/entry/K05286)) is Dol-P-Man + G00148 → Dol-P + G00149.
  - G00148 is (Man)2 (EtN)1 and G00149 is (Man)3 (EtN)1, the compositions of m601 and m602.
  - So R375 is the same step. **Verified.**
- In the model, Dol-P-Man (m502) minus Dol-P (m503) is exactly C6H10O5 (one mannosyl residue; the
  shared (C5H8)n is excluded). **Verified.**
- R372 and R374, the other two mannosyl transfers, already return m503 (verified, before.json).
- **Gene identity (verified):** YALI1D05270g is UniProt Q6CAB8 (GPI10_YARLI, "GPI
  mannosyltransferase 3", YALI0D04202g, reviewed), which is KEGG yli:2910629, assigned to K05286.
  This is database identity, not experimental validation.
- **Mass balance:** not computable by cobra, because the Dol-P formulas carry `(C5H8)n` and m601 and
  m602 have no formula.

**Why the chain got cheaper (verified numbers; mechanism provable from the stoichiometry).**
- Before, each GPI made:
  - consumed one whole Dol-P carrier at R375;
  - destroyed one PE's DAG moiety at R373.
- Now both return to the network: Dol-P to the ER pool, DAG through R_DAGt_er_em to the ER membrane.
- The growth cost of a forced 1e-5 exit fell from 2.83e-5 to 1.88e-5, and the FVA maximum at zero
  growth rose from 0.497 to 0.750.
- R_DAGt_er_em now carries two DAG per GPI (from R373 and R376), and R_PEt_em_er still two PE.

**No new reverse route (provable).**
- R373 and R375 stay reversible, but their reverse needs m600 or m602.
- Those come only from the chain, which starts at the irreversible R_GPIGNT and ends at the
  irreversible demand R_DM_GPI. m603 has no other producer.
- All 12 FVA minima are 0 (verified).

## 4. Builds and whole-model verification

All verified, [logs/sha256.txt](logs/sha256.txt).

| Build (combined argv plus) | Output | Check |
| --- | --- | --- |
| none | `164f235f…` | equals the original candidate byte for byte |
| `--gpi-chain-links` | `d5a0ce48…` (.scratch) | differs from `89559a0f` only in R_GPIGNT `gpi_chain_links_rationale` ([model_diff_links_only_vs_89559a0f.json](model_diff_links_only_vs_89559a0f.json), 1 difference) |
| `--gpi-chain-links --gpi-coproducts`, seed 0 | `cf0825b8…` (candidate) | the same bytes as seed 1 |
| the same, seed 1 | `cf0825b8…` (.scratch) | |

- **Diff, new vs links-only** ([model_diff_new_vs_links_only.json](model_diff_new_vs_links_only.json)):
  4 differences.
  - R373 stoichiometry: +m604[C_er] 1.0.
  - R375 stoichiometry: +m503[C_er] 1.0.
  - The 7 new notes on each.
- **Diff, new vs `89559a0f`** ([model_diff_new_vs_89559a0f.json](model_diff_new_vs_89559a0f.json)):
  the same plus the R_GPIGNT note.
- **Build log warning:** the R304 libSBML notes error appears in every build, including the no-flag
  control. It predates this change.

## 5. Solves (SD-Leu + PO1f, GLPK; [solves.json](solves.json), 58 LPs)

| Check | Baseline 89559a0f | New cf0825b8 |
| --- | --- | --- |
| WT growth | 1.4308057353528538, optimal | 1.4308057353528474, optimal |
| FVA max, all 7 chain reactions, R_GPIGNT, R_ACt_er_cy, R_DM_GPI | 0.497282 | 0.749679 |
| FVA max, R_PEt_em_er | 0.994565 | 1.499357 |
| FVA max, R_DAGt_er_em | 0.497282 | 1.499357 (2 DAG per GPI) |
| FVA min, all 12 | 0 | 0 |
| Forced exit ≥ 1e-5 (pFBA) | optimal; growth cost 2.83e-5; R_DAGt_er_em 1e-5 | optimal; growth cost 1.88e-5; R_DAGt_er_em 2e-5 |
| …with any one of the four links closed | not re-run | infeasible in each case |

TASK.md budgeted about 50 LPs. The run used 58, because FVA covered all 12 reactions on both
models.

## 6. Knockouts and tests

**Knockouts** ([ko_compare_summary.json](ko_compare_summary.json), [ko_compare.tsv](ko_compare.tsv);
2,150 LPs, 1,074 genes, `PO1f_plasmid_LEU2` excluded; a comparison of two models, no experimental
labels):
- No knockout growth differs by more than 1e-6, except YALI1D22124g.
- **Below the tolerance (found by the audit, verified):**
  - Nine knockouts grow faster by exactly +8.477e-7, all at KO/WT 0.127687: YALI0F04114g,
    YALI1E12730g, YALI1E15606g, YALI1E19902g, YALI1E23635g, YALI1F05089g, YALI1M00131g,
    YALI1M00277g and YALI1M00325g. In the links task these were unchanged (below 4e-13).
  - The cause is the open exit
    ([ko_exit_sink_check.json](ko_exit_sink_check.json), [scripts](scripts/ko_exit_sink_check.py),
    3 LPs). In the YALI1E12730g knockout (sole gene of R304, cytochrome-c oxidase):
    - the new model runs 2.14e-5 through the GPI chain;
    - it returns 4.28e-5 DAG via R_DAGt_er_em, instead of making DAG from phosphatidate through
      R1773;
    - with R_DM_GPI closed, growth is exactly the baseline's 0.18269487690…
  - When PE is in surplus, the chain plus the open exit act as a PE → DAG route that discards the
    GPI. Before the co-products, the chain lost a Dol-P and a DAG per pass and was never worth it.
  - The gain is far too small to change a class: 0.127687 is not near a cutoff.
  - "The chain carries 0 at the optimum" therefore holds for WT, not for every knockout.
- **YALI1D22124g:**
  - The `89559a0f` sweep gives `infeasible` again (unresolved). This is the same GLPK artifact the
    links report diagnosed, reproduced by the same loop.
  - The new model's sweep gives optimal 0.0 (essential).
  - Fresh solves ([ko_d22124g_fresh.json](ko_d22124g_fresh.json)) give −2.4e-17 on `89559a0f` and
    2.1e-17 on `cf0825b8`. So it is essential on both, and there is no class change.
- The sweep's essential counts at 1/5/10/15 % are 90/98/107/133 on the baseline and 91/99/108/134 on
  the new model. The difference is this one gene, and the new counts equal those of `164f235f`.

**Tests (verified):**
- **Focused tests:** `tests/test_gpi_coproducts.py` has 4 tests, with 14 subtests before the audit
  fix and 17 after. Run together with `tests/test_gpi_chain_links.py`, that is 9 tests and 37
  (later 40) subtests.
- **Full suite, per test,** against the previous run (the links task after its audit fix): the only
  difference is the 4 new passing tests (217 → 221 passed). The 11 failures, 17 errors (3 of them
  `pydantic_ai` collection errors) and 31 skips are identical by test ID.
  - [logs/pytest_before_results.txt](logs/pytest_before_results.txt)
  - [logs/pytest_after_results.txt](logs/pytest_after_results.txt)
  - [logs/pytest_summaries.txt](logs/pytest_summaries.txt)

## 7. Independent audit

The full audit is in [AUDIT.md](AUDIT.md). A fresh read-only reviewer solved 23 GLPK LPs and edited
nothing.
- **Supported:** C1, C3, C4, C7 and C8 (mechanism confirmed: about 70% of the saving is Dol-P
  recycling, 30% DAG return).
- **Supported with a caveat or nuance:** C2, C5, C10.
- **C6 (guard): supported for the listed cases, with gaps.**
- **Partially supported:** C9 (sub-tolerance knockout shift) and C11 (the two findings).

| Finding | Response |
| --- | --- |
| 1. The open exit acts as a PE → DAG sink in nine knockouts (+8.477e-7) | **Verified and reported** (§6, §8); not changed. The exit was approved open in the links task, so changing it is a user decision. |
| 2. The R376 finding's conclusion was misleading; YALI1E14925g is GPI13 | **Verified and corrected** in §8 and in TASK.md's corrections section. R376 likely equals R08107. The open issues are m603's name, the GPI7 step and the GPR logic. |
| 3. Guard gaps | **Fixed:** links are now checked against their recorded bounds and stoichiometry, not just presence; `requires_reactions` must name exactly the five links; curation notes on the before-state are rejected; the edit is applied as the full before → after difference. The test gained 3 probes. **Not fixed, accepted:** m604 annotation, m503's reaction context, R372's stoichiometry and m601's name are not guarded. m503 is a hub with 8 reactions; guarding it would make the step fail on unrelated changes. The guard fix leaves the model unchanged; see §7.1. |
| 4. The DAG finding was overstated for R1728 and R1773; R1725 and R1726 are closed | **Verified and corrected** in §8. |
| 5. Provenance gaps | `ko_d22124g_fresh.json` reuses the links script. Its `baseline_ko` is `89559a0f` and its `new_ko` is `cf0825b8`. Its `new_ko_links_closed` closes the five links on `cf0825b8`. Its `baseline_point_in_new_model` checks the `89559a0f` point in `cf0825b8`. The SHAs are in commands.txt. PYTHONHASHSEED is recorded only in commands.txt, as in earlier tasks. The LP overrun (58 vs about 50) is noted in §5. |
| 6. Information | Agreed. `baseline_model_sha256` is a record, not a guard. Inside the build the pre-state is the links-only model. |
| C2 caveat: "as R778 and R211 do" | R211 itself is unbalanced, so the note's "as R778 and R211 do" means they use the same species, not that they are balanced precedents. Recorded here. Rewording it again would change the model's bytes once more. |

### 7.1 Verification after the guard fix

All verified:
- **Focused tests:** `test_gpi_coproducts.py` and `test_gpi_chain_links.py` give 9 tests and 40
  subtests. The co-product test now has 17 subtests: the three new probes plus the empty
  `requires_reactions` scope case.
- **Verification rebuild:** `PYTHONHASHSEED=0`, written to
  `.scratch/gpi_coproducts_20261021/coproducts_after_audit_seed0.xml`, gives `cf0825b8…`, the same
  bytes as the candidate.
  - The LP and knockout results above apply unchanged.
  - The candidate's build.json keeps the first build's data SHA `825e429b…`; the current data file
    is `bc6b29ed…`.
- **Full suite:** the only difference from the links task's final run is the 4 new passing tests
  (221 passed). The 11 failures, 17 errors and 31 skips are identical by test ID.
  - [logs/pytest_after_audit_results.txt](logs/pytest_after_audit_results.txt)

## 8. Findings recorded, not changed (user decision), and limits

- **R376: probably right as written; m603's name and a missing GPI7 step are the likelier
  problems.** This is corrected after the audit; the TASK.md wording was misleading. All points
  verified 2026-10-21.
  - YALI1E14925g is GPI13. NCBI's feature table (`reference/ncbi/clib89_feature_table.txt`)
    maps YALI1_E14925g to GeneID 2912812 and XP_503848.3. The S2 locus map gives YALI0E12089g.
    KEGG yli:2912812 is K05288, "GPI ethanolamine phosphate transferase 3". The model gene simply
    has no identifiers.
  - [CHEBI:53060](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:53060), m603's
    annotation, has a name listing three EtNP groups, but its definition describes two. Its
    manual KEGG cross-reference is G13044, which is (Man)3 (EtN)2, the product of
    [R08107](https://www.kegg.jp/entry/R08107): PE + G00149 → DAG + G13044.
  - So R376 (m602 + PE → m603 + DAG) most likely is R08107 (EtNP onto Man3, GPI13), with one PE and
    one DAG, as the model has it. The remaining problems, all **candidates**:
    - m603's three-EtNP name, inherited from ChEBI's name;
    - the GPI7 step ([R12702](https://www.kegg.jp/entry/R12702), EtNP onto Man2) is absent from the
      model, and KEGG's GPI7 gene yli:2911071 is not among the model's genes (no D34684 entry);
    - R376's GPR is `YALI1E04744g or YALI1E14925g`, i.e. GPI11 (K05287) or GPI13 (K05288). KEGG
      lists both subunits for R08107, which suggests AND rather than OR.
- **DAG species, model-wide.**
  - The diglyceride pools (m604, m1649 and the others) carry formula C6H8O6, with an InChI for a
    small six-carbon compound (HMDB0251280). A 1,2-diacylglycerol carries two acyl chains (KEGG
    C00641), which C6H8O6 cannot hold. Verified.
  - Corrected after the audit:
    - R1725 and R1726 are mass-imbalanced but closed at (0, 0) in this candidate.
    - R1728 and R1773 contain species with no formula (m1656, m378), so their balance is not
      computable. I had wrongly called them unbalanced.
  - Fixing the formula reaches every DAG reaction, so it needs its own task. It will also trip
    both GPI guards, which pin m604's current formula; that needs coordinated data edits.
- **Open exit as a sink (audit finding 1, §6).** This is a property of the open demand, approved in
  the links task. Options for a later decision:
  - keep it, accepting the gain below 1e-6;
  - close R_DM_GPI by default and open it only for connectivity checks;
  - replace it with a measured GPI term in biomass.
- **Unchanged this round:**
  - R373 and R375 stay reversible.
  - R371's extra H⁺ remains.
  - The `P*2` notation remains.
  - Step 1 has no GPR.
  - There is no GPI biomass term.
- **Superseded candidate.** The `89559a0f` candidate stays as published. A `--gpi-chain-links`-only
  build from the current data now gives `d5a0ce48…`, which differs only in the reworded note.
- **Provenance:** this curation is development work, not independent validation.
