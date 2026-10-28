# Kennedy closures: R1726 opened, R1725 kept closed, candidate ce9e1112 (2026-10-27)

## 1. Result

- **Candidate:**
  [`model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_20261027.xml`](../../candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_20261027.xml),
  SHA-256 `ce9e1112ee940531469457009234a2fd11de7f589c5948a7dec6b91118005294` = `e9d35857` (R2135 fix) +
  `--r1726-cpt-open`. Not a designated reference. Not committed.
- **Change:**
  - **R1726** (cholinephosphotransferase, CDP-choline + DAG → PC, GPR `YALI1C15420g`): (0, 0) becomes
    (0, 1000), plus six notes.
  - **R1725** (CDP-ethanolamine → PE) stays (0, 0).
- **Authorization** (TASK.md sections 0 and 13), both from the user on 2026-10-27:
  - "Approve Option A", for the first draft;
  - "Approve re-run", for the re-run draft.
- **History:**
  - The first run on `00e8e397` stopped at the knockout check ([STOP.md](STOP.md)). Opening R1726 closed a
    loop through R2135's acetyl-CoA → acyl-CoA mapping, and fatty-acid synthase stopped being essential.
  - That mapping was fixed separately (`../r2135_acyl_coa_pool_20261027/`).
  - The re-run builds on the fixed candidate.
- **Effect:**
  - UM12's 11 reactions are unblocked in the open condition.
  - Supplied choline becomes usable.
  - On the Task A medium + PO1f: no change to growth and no change to any knockout.

## 2. Evidence (from TASK.md section 4; tiers as given there)

- **For R1726:**
  - `YALI1C15420g` = YALI0C10989g is the CPT1 ortholog (KEGG K00994; orthology, verified).
  - *S. cerevisiae* CPT1 is relatively specific for CDP-choline (PMID 1847919) and makes 95 % of the PC that
    comes from CDP-choline in vivo (PMID 7961735). It adds net PC only when choline is supplied (UniProt P17898).
    Verified, *S. cerevisiae*.
- **Against R1725:**
  - Its only gene is the same CPT1 ortholog.
  - The EPT1 ortholog `YALI1E31498g` is not in the model.
  - Opening it lets the sphingoid-base lyase rescue a PSD1 + PSD2 double knockout, while *S. cerevisiae*
    psd1 psd2 mutants are ethanolamine auxotrophs (PMID 7890740).
- **Not found:** any *Y. lipolytica* enzyme data, or a recorded reason for the inherited closure.

## 3. Implementation

- **Data:** [`model/curation/direction/r1726_cpt_open.json`](../../curation/direction/r1726_cpt_open.json).
  - R1726 before / after bounds, stoichiometry and species.
  - The kept-closed R1725 record (name, GPR, bounds, stoichiometry, species).
  - The evidence and notes.
  - `reviewed_model` is `e9d35857`.
  - The `limits` note was rewritten for the re-run: no Task A growth gain, and the R2135 fix is required.
- **Step:** `platform/scripts/gem_annotate/r1726_cpt_open.py`. It fails closed on:
  - scope;
  - a missing reaction;
  - R1726 name, GPR, EC, bounds, stoichiometry or species;
  - R1725 not exactly the reviewed closed reaction;
  - conflicting notes;
  - open without the notes;
  - closed with the notes (added in the re-run review).

  It touches only R1726's bounds and notes.
- **Flag:** `--r1726-cpt-open` requires `--r2135-acyl-coa-pool`, which requires `--cofactor-biomass-terms`,
  and runs after it.
- **Test:** `platform/tests/test_r1726_cpt_open.py`, 5 tests and 25 subtests on the `e9d35857` fixture, each
  guard matched by its message.
- **Code review** (fresh general-purpose agents following `.claude/agents/code-reviewer.md`, twice):
  - **First run:** no bugs and four low items. All were fixed: a stale comment, the R1725 stoichiometry guard,
    the guard tests, and a flag requirement.
  - **Re-run:** no bugs; it checked all 8 flag combinations. Its three low items were fixed: the closed-with-notes
    guard, the R1725-only species and closed-with-notes tests, and a stale comment.
  - Neither round of fixes changed a model byte.

## 4. Verification (re-run; outputs in [rerun/](rerun/))

| Check | Result |
| --- | --- |
| pytest | 278 → 278 passed (subtests 781 → 783, the two new guard subtests), the same 11 failures and 18 errors. The only difference is the renamed requirement test (`test_flag_requires_the_cofactor_biomass_terms` → `test_flag_requires_the_r2135_acyl_coa_pool`) ([rerun/logs/pytest_compare.txt](rerun/logs/pytest_compare.txt)). "Before" is the R2135 task's final run. |
| Builds | Control without `--r1726-cpt-open` = `e9d35857…` byte for byte. `PYTHONHASHSEED=0` and `=1`, written as the first line of each build log, give `ce9e1112…` in round 1 and in the final round ([sha256.txt](rerun/logs/sha256.txt), [sha256_round1.txt](rerun/logs/sha256_round1.txt)). |
| Whole-model diff vs `e9d35857` | 2 entries: R1726 bounds [0, 0] → [0, 1000], R1726 notes (+6 keys) ([model_diff_vs_e9d35857.json](rerun/model_diff_vs_e9d35857.json)) |
| WT, Task A + PO1f | 0.3116824 → 0.3116824, optimal, checked ([built_checks.json](rerun/built_checks.json)) |
| R1726 / R1725 FVA | R1726: 0 at 100 %, [0, 0.0217] at 99.9 %, [0, 21.3] at ≥ 0 growth. R1725 [0, 0]. |
| No-uptake LPs | Max ATP and max biomass are 0, optimal, on both models (`xMAINTENANCE` relaxed) |
| PSD1 + PSD2 double knockout | Lethal on both |
| FAS pair, desaturase, acyl-CoA ligase knockouts | Lethal on both (the first run's failure is gone) |
| Choline and ethanolamine supplied (uptake 1) | PE-methylation knockouts `YALI1E07347g` / `YALI1E15432g`: 0 on the baseline → 0.3117034 on the candidate. That equals the candidate's WT with choline supplied; supplied choline raises candidate growth by 2.1e-5 over the Task A WT 0.3116824. |
| Knockout sweep | ([ko_compare_summary.json](rerun/ko_compare_summary.json))<br>- 2 × 1,082 LPs, Gurobi 1e-9, reset per knockout, knocked-flux check.<br>- No failing point; largest bound violation 6.8e-10, largest residual 8.6e-13.<br>- **No class change at any cutoff, and no knockout growth changed by more than 1e-6.**<br>- At 15 %: TP / FN / FP / TN 111 / 204 / 42 / 692 (MCC 0.383) on both.<br>- Positive coverage: 326 consensus genes are in the model (111 + 204 + 11 unresolved).<br>- 20 infeasible knockouts unresolved on both; the 16 outside mtDNA are excluded from the counts.<br>- User-defined negatives; mtDNA out of scope; development-reference scores. |
| Re-census (GLPK, vs the `e9d35857` census) | ([census_compare_open.json](rerun/census_compare_open.json), [census_compare_taskA_po1f.json](rerun/census_compare_taskA_po1f.json))<br>- **Open condition 720 → 709:** exactly UM12 (labelled UM11 in the `e9d35857` open-condition census run, whose module numbering shifted), namely R1726, R212, R213, R1603, R830, R1041, R1042, R1743, R1842, R1809 and R354. Gap metabolites 675 → 668.<br>- **Task A 852 → 843:** 9 of them; choline uptake R830 / R1041 has no choline in the medium.<br>- No reaction newly blocked in either condition. |
| Mass balance of R1726 | Formally unbalanced as inherited (charge +1; C, H and O −1; R +2): the lumped DAG `m1649` has no R groups. Its active route neighbours are unbalanced in the same lumped-lipid convention. Not changed (lipid-unlump line). |

## 5. What it means

- **Network:** the CDP-choline branch of the Kennedy pathway is connected, so a medium with choline can now be
  simulated with PC made from it. On the choline-free Task A medium the route adds nothing at the optimum. It
  can run only as a choline-neutral PC turnover cycle (PLD → choline → CDP-choline → PC).
- **The first run's growth gain (+5e-5) was an artefact** of the R2135 leak, not of the pathway.
- **Essentiality:** unchanged. PLD `YALI1E22633g` stays a false negative, and this change was not made to
  alter any call.

## 6. Limits

- **Evidence:** orthology plus *S. cerevisiae* enzymology; no *Y. lipolytica* data. The inherited closure's
  reason is unknown.
- **Unbalanced chemistry:** R1726 stays formally unbalanced in the lumped-lipid convention.
- **R1725 / UM19:** stay blocked. Option B (EPT1 ortholog plus a lyase-capacity decision) is not done.
- **The DAG kinase R1720** stays closed.
- **Coverage:** only the Task A medium was tested, with the choline / ethanolamine mechanism control.
- **Budget:** 3 builds per round. Round 1 plus a final round after the review fixes made 6 in all, against 3
  planned. About 50 built-file LPs, against about 40 planned.

## 7. Next gates

- R1725 / EPT1 (`YALI1E31498g`) and sphingoid-lyase capacity (TASK.md section 11).
- R1720 (DAG kinase).
- The remaining 0.01-scaled reactions of the R2143–R2187 block.
- The dolichol / N-glycan chain (the gap-census priority list).
- Each needs its own authorization.

## 8. Audit

[AUDIT.md](AUDIT.md).
