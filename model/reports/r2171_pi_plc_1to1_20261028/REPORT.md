# R2171 written 1:1 (PI + H2O → DAG + Ins1P + 2 H+), candidate 965ab947 (2026-10-28)

## 1. Result

- **Candidate:**
  [`model/candidates/E5_…_dolichol_route_v2_r2171_1to1_20261028.xml`](../../candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_r2171_1to1_20261028.xml),
  SHA-256 `965ab947f2aadc2c69443e7a794358c502acfdf0892e78582278ff07a37429d7` = `9d655c2b` + `--r2171-pi-plc-1to1`.
  Not a designated reference. Not committed.
- **Authorization:** "Approve Option A" (user, 2026-10-28), the sentence in TASK.md section 9.
- **Change:** R2171 (PI phospholipase C, GPR `YALI1D13510g`, KEGG R03332) goes from
  `PI + 100 H2O → DAG + 100 Ins1P + 100 H+` (charge −100) to `PI + H2O → DAG + Ins1P + 2 H+` (charge 0), plus six notes.
  Name, GPR, bounds and annotations unchanged.
- **Effect:** the inositol leak found in Stage 1 is gone. An INO1-null strain no longer grows without inositol,
  and no inositol or Ins1P can be made without INO1 or uptake. WT growth, every knockout class and the Task A
  census are unchanged.

## 2. Why (evidence; tiers from TASK.md section 2, sources opened 2026-10-28)

- **Leak (verified in runs, Stage 1):** on `9d655c2b`, with INO1 knocked out and no inositol uptake, growth equalled WT
  on four media; inositol / Ins1P demand reached 511–520 / 1000; R2171 created organic phosphate and carbon
  ([Stage 1 REPORT](../lipid_source_leak_check_20261028/REPORT.md)).
- **Reference form (verified):** the model's own link [KEGG R03332](https://www.kegg.jp/entry/R03332) and ModelSEED
  rxn06362 are the 1:1 hydrolase; no reference with 100× was found.
- **Gene (verified, weak):** `YALI1D13510g` = YALI0D10813g, KEGG K01771 (EC 4.6.1.13, which KEGG calls bacterial and
  which forms the cyclic phosphate); UniProt Q6C9J1 unreviewed, score 1, PI-PLC X domain only. The real activity and
  product are unknown, so the GPR was kept and the limit recorded in the notes.
- **2 H+ (verified):** they balance charge under the model's existing charges (PI 0, DAG 0, Ins1P −2). With the common
  DAG formula C5H6O5R2 they also balance hydrogen (inferred by the code reviewer; the model's DAG has no formula).
- **Why not B or C:** B (cyclic-phosphate lyase) needs a new species and a route without yeast evidence; C (close)
  would remove a reaction whose gene has a PI-PLC domain without evidence the activity is absent.

## 3. Implementation

- **Data:** [`model/curation/chemistry/r2171_pi_plc_1to1.json`](../../curation/chemistry/r2171_pi_plc_1to1.json): before and
  after equations, species guards (formula, charge, compartment, MetaNetX / KEGG ids), name, GPR, bounds, KEGG link,
  notes, evidence rows, reviewer, authorization.
- **Step:** `platform/scripts/gem_annotate/r2171_pi_plc.py`, `apply_r2171_pi_plc_1to1`. Before- and after-equations
  pinned in code; fails closed on scope, a missing reaction or species, any species / name / GPR / bounds / KEGG-link
  change, conflicting notes, the notes on the old equation, the new equation without the notes, empty notes or notes
  with `<`, `>` or `&`. Returns `applied` / `already_correct`; touches only R2171's stoichiometry and notes.
- **Flag:** `--r2171-pi-plc-1to1`, requires `--dolichol-route` and runs after it (`cli.py`, `main.py`).
- **Test:** `platform/tests/test_r2171_pi_plc.py`, 4 tests and 30 subtests on the `9d655c2b` fixture (read through a
  `Path`), every guard matched by its message.
- **Code review** ([logs/code_review.md](logs/code_review.md)): no bugs; two low findings (before-equation not pinned;
  notes not checked), both fixed and tested. The round-2 builds after the fixes give the same bytes.
- **Candidates README:** row added for `965ab947`.

## 4. Verification

| Check | Result |
| --- | --- |
| pytest | 314 → 318 passed (the 4 new tests), the same 19 failures and 14 errors, subtests 835 → 865 ([logs/pytest_compare.txt](logs/pytest_compare.txt)). Focused (R2171 + dolichol): 10 passed, 82 subtests. Run with the workflow-automation venv (miniconda has no pytest). |
| Builds | Control without the flag = `9d655c2b` byte for byte. `PYTHONHASHSEED` 0 and 1 give `965ab947` in round 1 and again after the review fixes ([logs/sha256.txt](logs/sha256.txt), [logs/sha256_final.txt](logs/sha256_final.txt)). Each build log repeats a pre-existing libSBML notes error on R304, identical in the control and in the dolichol build. **Code identity:** the candidate's own `.build.json` records the round-1 step (`r2171_pi_plc.py` `d3016b8d…`); the final code (`dbc8f25f…`, data `8763ce20…`) is recorded in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r2171_pi_plc_1to1_20261028/final_seed0.build.json` / `final_seed1.build.json`, which built the same bytes. |
| Whole-model diff vs `9d655c2b` | 2 entries: R2171 stoichiometry, R2171 notes (+6 keys) ([model_diff_vs_9d655c2b.json](model_diff_vs_9d655c2b.json)) |
| Built-file checks (Task A + PO1f) | WT 0.3116608216 → 0.3116608216 (difference 3e-16), checked. R2171 charge −100 → 0. No-uptake max ATP and biomass 0, optimal, both. R2171 FVA at 99.9 %: [0, 1.4e-6] → [0, 0.0136]; at growth ≥ 0: [0, 1.2e-3] → [0, 13.2] ([built_checks.json](built_checks.json)). |
| **Acceptance** (Stage 1 check on the candidate) | INO1 knockout without inositol: growth **0** and inositol / Ins1P demand **0** on C1–C4 (Stage 1 baseline: WT growth, 511–520 / 1000); the same with R1720 opened on C1. Acyl and phosphate removal: 0 everywhere. Oleate positive control grows (0.3204501). 84 LPs, every point passed ([acceptance_check.json](acceptance_check.json)). |
| Knockout sweep | 2 × 1,082, Gurobi 1e-9, reset per knockout, no failing point. **No class change at any cutoff.** One growth change above 1e-6: the inositol transporter `YALI1F31381g` 0.3116608 → 0.3116574 (INO1 now has to make the inositol the leak gave for free). At 15 %: TP / FN / FP / TN 112 / 203 / 42 / 692 (MCC 0.386) on both; Positive coverage: 326 consensus genes in the model (112 + 203 + 11 infeasible knockouts scored as unresolved); 16 unresolved infeasible knockouts outside mtDNA in all, the same on both. The re-run baseline equals the dolichol task's `a801aad0` sweep exactly ([ko_compare_summary.json](ko_compare_summary.json)). User-defined negatives; mtDNA out of scope; development-reference scores. |
| Re-census (GLPK) vs the dolichol `a801aad0` census | Open 708 → 705: R1343, R1344, R970 freed, none newly blocked, gap metabolites 665 → 665. Gurobi FVA shows all three unblocked in **both** models (max 1.47e-4 each), so the difference is not caused by the change; why GLPK called them blocked in the `a801aad0` run was not identified (1.47e-4 is far above its 1e-7 cutoff) ([fva_gurobi_open.json](fva_gurobi_open.json)). Task A 840 → 840. ([census_compare_open.json](census_compare_open.json), [census_compare_taskA_po1f.json](census_compare_taskA_po1f.json)) |
| Before / after (no inositol, INO1 intact) | Before: R2171 2.4e-7, INO1 0, growth 0.3116608. After: R2171 0, INO1 2.34e-5, growth 0.3116574. After, INO1 knocked out: growth 0 ([figures/](figures/)). |

## 5. What it means

- **The model no longer makes inositol or organic phosphate from water.** The INO1-null strain is again an inositol
  auxotroph, as *S. cerevisiae* ino1 mutants are (established / uncited).
- **Nothing on Task A changes**, because the medium supplies inositol: same WT, same knockout classes, same census.
  Inositol-free conditions and any carbon / phosphorus accounting through inositol phosphates are now correct in
  this respect.
- **R2171's wider FVA range** (inferred explanation, route not traced): in the 1:1 form PI → DAG + Ins1P → inositol →
  PI can run as a balanced turnover cycle that costs CTP, so it can carry flux within the 0.1 % slack; the 100× form
  could not, because its extra inositol had no outlet. The acceptance check shows it is not a source.
- **R1720 can now be considered:** with R2171 fixed, opening R1720 in memory creates no inositol, acyl or phosphate
  source on C1.

## 6. Method notes

- The acceptance test is Stage 1's source-removal check, which found the leak; a fix that passes it is shown not to
  be a source in the same tests, not in every possible condition.
- The census baseline is the dolichol task's `a801aad0` run, recorded there as the same optimization problem as
  `9d655c2b`. GLPK's blocked calls for R1343 / R1344 / R970 differ between that run and this one although Gurobi finds
  them unblocked in both; census differences of a few reactions should be checked with Gurobi, as done here.
- The sweep baseline was re-run rather than reused. TASK.md allowed a re-run only if the baseline could not be reused;
  the dolichol sweep has no `9d655c2b` column (its sides are `ce9e1112` and `a801aad0`), and running both models in one
  run of the same script avoided relying on the recorded same-LP claim. The re-run reproduces the `a801aad0` column
  exactly.

## 7. Limits

- **Gene identity:** `YALI1D13510g` is orthology plus an automatic X-domain annotation; whether it makes Ins1P or the
  cyclic phosphate is unknown. No *Y. lipolytica* data.
- **Mass balance:** not checkable (PI and DAG are lumped species).
- **Coverage:** Task A + PO1f plus the three other media of the acceptance check; one base model.
- **Not changed:** the IPC / MIPC synthases (0.01 PI), R2165, R2187, R1720; the pre-existing R304 notes error.
- **Paths:** the candidate's absolute path (≈ 289 characters) exceeds cobra's macOS 260-character limit for string
  paths; every check here used byte-identical short-path copies with recorded SHA-256.
- **Budget:** builds 5 of 6; acceptance 84 of ≤ 150 LPs; built-file checks 18; sweep 2,164 (the baseline re-run used
  the conditional +1,082, see section 6); fresh re-solves 0; census 2 GLPK runs. **Not in any planned category:** the
  Gurobi FVA of the three census reactions (14 LPs) and the figure points (7 LPs); counted with the built-file checks
  they make 39 against the planned 30.
- **Evidence rows:** E6 (no inositol cyclic phosphate species) and E7 (R1740 separate) are static model checks recorded
  in TASK.md, not evidence rows of the curation JSON.

## 8. Next gates

- Opening the DAG kinase R1720 (Kennedy line), now that R2171 is fixed.
- The remaining scaled reactions of the inherited block: IPC / MIPC synthases, R2165, R2187.
- R1725 / EPT1 after the sphingoid-lyase decision.
- Each needs its own authorization.

## 9. Audit

[AUDIT.md](AUDIT.md): a fresh read-only auditor reproduced the SHAs, the diff, WT, the C1 acceptance points, the
knockout recomputation and the figure points (27 LPs). K1–K3, K5–K11 supported; K4 and K12 partially supported
(E6 not a JSON evidence row; the baseline sweep re-run and 21 uncategorised LPs not explained). No science overclaim.
Its documentation gaps were fixed as listed there.
