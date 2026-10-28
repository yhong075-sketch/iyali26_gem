# Independent audit: combined TPP + heme/iron + R2193 candidate (2026-10-22)

A fresh read-only reviewer, with no history of this work, audited the combination against claims
C1–C8 in the brief.
- It edited nothing, built nothing and committed nothing.
- It used 22 LPs (limit 40): 4 WT, 3 minimum-flux, 10 knockout spot-checks and 5 reset tests.
- Its scratch files are in this session's scratchpad (not kept).

Its report is reproduced below, lightly reformatted. Responses are in
[REPORT.md §6](REPORT.md#6-independent-audit-and-responses).

## Verdicts

| Claim | Verdict | Evidence |
| --- | --- | --- |
| C1 SHAs; TPP parent biomass_C = base + m975 + 7 note keys = `tpp_requirement.json` after-state | **Supported** | All four SHAs match. TPP parent vs base adds only `m975[C_cy]` −1.2e-6 (nothing removed or changed, bounds equal) and exactly 7 `tpp_requirement_*` keys. It equals the JSON after-state and the heme JSON's state 2. |
| C2 Guard | **Supported for biomass_C composition and note keys; gaps in issues 1–3** | `BEFORE_STATE_IDS` pinned; state 0 adds nothing; no state carries m1635/m948. Probes: base → `base`; 2fe51cc2 → `after_TPP…`; an extra key and a 1e-9 relative TPP coefficient change are rejected; TPP after heme is rejected. Both controls reproduce byte for byte. **Not pinned:** TPP note values (an altered value was accepted), m975 identity/chemistry (the TPP-only v2 candidate, m975 formula None / charge 0, was accepted), and state 2's content (pinned in JSON only; a widened state was accepted). Step order is fixed in main.py; the reverse order fails closed. |
| C3 Builds | **Supported; seed values unverified** | control_tpp_only = 2fe51cc2, control_heme_only = e1586027, combined and seed-1 both 2c9cebd0. build.json flags as claimed. Combined logs record `before_state: after_TPP-REQUIREMENT-20261022`; the heme-only control records `base`. `PYTHONHASHSEED` is only in commands.txt (inherited gap). |
| C4 Diffs | **Supported (stronger check)** | Raw XML diffs: base→2fe51cc2 equals e1586027→combined (234 lines, `cmp` identical), and base→e1586027 equals 2fe51cc2→combined (95 lines, identical). Nothing else differs anywhere in the file. |
| C5 WT and routes | **Supported (reproduced)** | Gurobi 1e-9 WT values reproduced; additivity residual 9.2e-15. Minimum-flux LPs at growth ≥ (1−1e-6)·WT: R752 ≥ 3.749e-7, R302 = R1709 ≥ 3.124e-7. Both routes are forced, not a pFBA vertex. |
| C6 Three-model sweep | **Supported** | Independent rescore from the TSV matches all 12 score tuples and MCCs. Every class column agrees with KO/WT < cutoff. Magnitudes: knocked 0, violation 0, max residual 2.97e-12. The same 20 infeasible genes in all three. Class changes as claimed. 10 knockouts reproduced exactly. **Reset:** after `problem.reset()` VBasis is unreadable, and iterations go from 136 warm to 1,524 after reset (34 to 2,001). Gurobi's docs say the next optimize starts from scratch. |
| C7 Full pytest | **Supported** | No test changed outcome or disappeared. 21 additions: 1 heme test, 19 TPP/thiamine tests, and 1 collection error in `test_gem_review_lit_caching.py` (`httpx2` missing; it arrived via the main merge, 81e9d60; `httpx2` is in uv.lock but not installed). 251 / 11 / 18 / 31 against 231 / 11 / 17 / 31. |
| C8 Heme parent Gurobi = GLPK scores | **Supported for the numbers; file named wrongly** | Scores and per-gene classes match (0 differences, max \|ΔKO\| 2.7e-11). But the GLPK sweep's `new_sha256` is 2cf19c69 (v2), not e1586027. Equivalence rests on `model_diff_v3_vs_v2.json` (R2193 notes only); the v2 file was deleted. |

## Issues

1. **Medium: the guard does not check m975 chemistry or the TPP-line state outside `biomass_C`.**
   `--tpp-requirement --heme-iron-biomass` without the chemistry and lock steps passes every guard,
   which would build an unreviewed combination. Require it in the guard or CLI, or state the limit.
2. **Low:** TPP note values are not pinned, only the keys. Builds are unaffected, because the TPP step
   checks its own values first.
3. **Low:** state 2's content is pinned only in JSON. The combination depends on the fixed order in
   main.py (TPP before heme); this is not stated.
4. **Medium (process):** acceptance criterion 6 is incomplete.
   - REPORT §3 and §6 said the full suite was "pending".
   - There is no STATE.md entry.
   - `model/candidates/README.md` changed during the audit (a 2c9cebd0 row appeared), so "another
     session" may be editing concurrently.
5. **Low–medium: reporting gaps.**
   - Recall, specificity and precision, and the out-of-scope positive count (1,612 − 324 = 1,288), are
     missing from the report.
   - It does not say that 11 positives are among the 16 unresolved knockouts.
   - Genes are named by ID only, with no verified name, function or evidence level.
6. **Low:** C8 wording should say the GLPK sweep ran on v2 `2cf19c69`, which is the same optimization
   problem.
7. **Low:** §4's "a model reset before every solve" and "Every WT point passes the check" do not
   cover the pFBA LPs, which were neither reset nor checked. The route fluxes in the table come from
   those points; the minimum-flux LPs confirm the values.
8. **Low (provenance):**
   - the reference diff behind "identical in fields and values" was not kept;
   - "130 subtests" has no kept log;
   - `PYTHONHASHSEED` is not in build.json.
9. **Low (scope wording):** "combine without interaction" is shown only for WT and single-gene
   knockouts, on one medium, at 1e-6 growth resolution. The class-level union is fully verified.
10. **Information:**
    - Timing supports "definitions fixed before the sweep" (TASK.md 16:31:33, log 16:31:45).
    - Budget was within plan.
    - "13 expected vs 12 found" is consistent: YALI1F07747g stops being essential because of the R2193
      swap. The report does not explain this.
    - No biological over-claims beyond the two lines' own.
