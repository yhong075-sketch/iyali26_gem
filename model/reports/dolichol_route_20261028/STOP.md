# STOP at the knockout gate: dolichol route candidate a801aad0 (2026-10-28)

**The run stopped at the stop condition of TASK.md section 7, step 8.** Nothing was relaxed or retried. The
candidate file is built and kept; the independent audit, REPORT.md and the STATE.md entry wait for the
user's decision. Nothing is committed.

**Decision (2026-10-28):** the user answered "Accept, restate gate", the first option of section 5
(TASK.md section 11). The run continued from here; see [REPORT.md](REPORT.md) and [AUDIT.md](AUDIT.md).
The final candidate is `9d655c2b…` (notes corrected after the audit, same optimization problem as `a801aad0`;
TASK.md section 12), so the results below hold for it.

## 1. The gate as approved

> any class change other than `YALI1E22718g` non-essential → essential, or any growth change above 1e-6
> outside the route

## 2. What the sweep found (Task A + PO1f, Gurobi 1e-9, verified)

Sweep: 2 × 1,082 LPs ([ko_compare_summary.json](ko_compare_summary.json), [ko_compare.tsv](ko_compare.tsv)),
reset per knockout, knocked-flux check. No point failed the check: largest knocked flux 0, bound violation
1.3e-10, residual 3.8e-11. The 20 infeasible knockouts are the same 20 on both models; no fresh re-solves.

- **Class changes:** only `YALI1E22718g` (SEC59 ortholog), non-essential → essential, at 1 %, 5 %, 10 % and
  15 %. This half of the gate passes.
- **Growth changes:** the second half fails as written.
  - In absolute terms every knockout changes, because WT itself drops by 2.16e-5 (0.3116824 → 0.3116608,
    −6.9e-5 relative). The gate did not anticipate a WT shift.
  - Judged on KO/WT ratios ([stop_condition.json](stop_condition.json),
    [scripts/stop_condition_check.py](scripts/stop_condition_check.py)), **104 genes outside the route move by
    more than 1e-6**. The largest move is 3.2e-5 (`YALI1F23441g`, 0.466314 → 0.466346).

## 3. What the 104 shifts look like ([ratio_shift_profile.json](ratio_shift_profile.json))

| Measure | Value |
| --- | --- |
| Genes with \|ΔKO/WT\| > 1e-6, outside the route | 104 (64 of them on the consensus list) |
| Largest \|ΔKO/WT\| | 3.2e-5 |
| Relative change of their absolute KO growth | −8.9e-5 to +4.7e-12 |
| Smallest distance of any of them from a cutoff | 4.9e-3 (`YALI1E34269g`, ratio 0.0549 vs 5 %) |
| WT relative change | −6.9e-5 |

**Reading (inferred, not separately tested):** no knockout grows better on the candidate; each either keeps
its absolute growth (its limit is elsewhere, so its ratio rises as WT falls) or pays up to the Dol-P cost
itself. That is the arithmetic of a new biomass term, not an escape route: a bypass would raise a knockout's
absolute growth or flip a class. The gate's 1e-6 threshold came from the Kennedy task, where WT did not move.

## 4. Everything else that ran (all verified this session)

| Check | Result |
| --- | --- |
| pytest | 278 → 284 passed, with the same 11 failures and 18 errors before and after (pre-existing); the 6 new tests are the only difference ([logs/pytest_compare.txt](logs/pytest_compare.txt)) |
| Code review | No bug that can write a wrong model; the low fail-closed gap and test gaps fixed before the builds ([logs/code_review.md](logs/code_review.md)) |
| Builds | Control = `ce9e1112` byte for byte; `PYTHONHASHSEED` 0 and 1 both give `a801aad0…` ([logs/sha256.txt](logs/sha256.txt)) |
| Whole-model diff | 6 entries: R2178 bounds and notes, R_DOLSYN added, biomass_C stoichiometry (m503[C_er] −1e-4) and notes, N-Glycan group members ([model_diff_vs_ce9e1112.json](model_diff_vs_ce9e1112.json)) |
| Built-file checks | WT 0.3116608 optimal; route forced at 3.1e-5 at the optimum; no-uptake ATP and biomass 0 (optimal, maintenance relaxed), also with R_DOLSYN reversible; SEC59 KO 0, CWH8 KO 1.0, single cis-PT KOs 1.0, the double cis-PT KO 1.0 → 0 ([built_checks.json](built_checks.json)) |
| Re-census, GLPK | Open 709 → 708: R2179, R2180, R287 freed; R1343 / R1344 newly blocked. Task A 843 → 840: the same three freed, nothing newly blocked ([census_compare_open.json](census_compare_open.json), [census_compare_taskA_po1f.json](census_compare_taskA_po1f.json)) |
| R1343 / R1344, Gurobi 1e-9, open condition | Unblocked on both models, max 1.47e-4: the GLPK fringe ([fva_gurobi_open.json](fva_gurobi_open.json)) |

## 5. Decision needed

- **Accept with the gate restated** for a change that adds a biomass term: "no class change other than
  `YALI1E22718g`, and no knockout gains absolute growth above 1e-6". Both parts are met. The run then continues
  with the independent audit, REPORT.md and the STATE.md entry.
- **Or keep the gate as written**: the candidate is not accepted, and the 104 shifts are investigated first.
