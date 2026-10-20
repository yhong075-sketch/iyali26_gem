# R2041 forward-only through the pipeline: report

Scope, authorization and budget are in [TASK.md](TASK.md). The evidence and the decision are in
[r2041_direction_review_20261019/REPORT.md](../r2041_direction_review_20261019/REPORT.md).

Evidence tiers:
- **verified:** run or checked this session (2026-10-19) and shown here or in a linked file;
- **provable:** follows mathematically;
- **inferred:** reasoned, not tested.

## 1. Summary

- **Change:** R2041 (ATP:D-arabinose 5-phosphotransferase, EC 2.7.1.54, C_cy, orphan) goes from
  [−1000, 1000] to [0, 1000]. Its balanced 2-H⁺ equation is kept.
  - The change is curated data `model/curation/direction/r2041_direction.json`
    (`R2041-FORWARD-20261019`) plus the guarded step `apply_r2041_direction`, applied after the
    metadata selection.
  - It is opt-in with `--r2041-forward`.
- **New candidate:**
  [E5_coq9_resp_pkg_r132_r2041_forward_20261019.xml](../../candidates/E5_coq9_resp_pkg_r132_r2041_forward_20261019.xml),
  SHA-256 `da0c50ba10453a854ceb8de99980e6e7465c2c3eda1893412a1f01c02e81265e`.
  - Baseline: a56f3649.
  - Build: a56f3649's recorded options plus `--r2041-forward`, at HEAD `285eaee` with the
    uncommitted change.
- **Verified:**
  - **Control:** a build without the flag reproduces `a56f3649…` byte for byte, before and after
    the audit fix.
  - **Whole-model diff:** exactly 2 differences against a56f3649, R2041 bounds and R2041 notes. In
    the XML that is `reversible`, the lower-bound parameter reference and 6 added note lines.
  - **Balance:** R2041 stays element- and charge-balanced.
  - **Tests:** the full suite equals the before-run test by test, except the 4 new passing tests.
- **The 10 authorized LPs (verified, all optimal):**
  - **WT:** identical on M0 and M1 (1.4308057353528 h⁻¹; relative difference 5.6e-15).
  - **Closed-exchange ATP:** dissipation maxima are 0 in C_cy and C_mi on both models.
  - **R2041 at the WT optimum:** a56f3649 allows R2041 down to −999.66, carried entirely by the
    net-zero internal loop R2041 rev / R2042 / R2190 rev / R1964 rev / R803 / R714 rev. It allows
    no forward flux. The new candidate allows R2041 = 0 only.
- **What it means:** the change removes a thermodynamically disfavoured reverse direction and a
  1000-capacity internal cycle.
  - It changes no growth, no ATP balance, and no flux at the recorded optimum.
  - There was no ATP-generating loop before or after.
- **Independent audit:** C1–C10 all supported.
  - One latent guard defect was fixed (§5).
  - Two note-wording points are reported, not changed (§5). Changing them would change the model
    SHA after the LPs.

## 2. What changed

| File | Change |
| --- | --- |
| `model/curation/direction/r2041_direction.json` (new) | The decision as data. It holds: the authorization quote; before and after bounds; the exact 2-H⁺ stoichiometry; the five species (formula, charge, compartment); GPR "", EC 2.7.1.54; the reviewed model and SHA; the review link; the reviewer; four evidence rows (ENZYME, Rhea and eQuilibrator verified; Volk 1962 marked unread); six notes |
| `platform/scripts/gem_annotate/r2041_direction.py` (new) | `apply_r2041_direction`, modeled on `r1159_direction.py` (details below) |
| `platform/scripts/gem_annotate/cli.py` | `--r2041-forward` (`store_true`) |
| `platform/scripts/gem_annotate/main.py` | `r2041_forward=False` parameter. When true, the step runs right after `apply_r1159_direction` and its record goes to `selection["post_selection_r2041_direction"]` |
| `platform/tests/test_r2041_direction.py` (new) | 4 tests and 9 fail-closed subtests, on the a56f3649 file read into memory |
| `model/curation/README.md` | One section for the opt-in step |
| `model/candidates/README.md` | One row for the new candidate |

The guarded step:
- It fails closed unless all of these match the record: the scope fields, the reaction name, the empty
  GPR, EC 2.7.1.54, bounds in {before, after}, the exact stoichiometry, and the species formula,
  charge and compartment.
- It also fails on conflicting curated notes, and on bounds already forward-only without this
  curation's notes. That last check was added after the audit (§5).
- It sets the bounds and merges the six notes, nothing else.

Not changed:
- the builder patch `fix_d_arabinokinase_direction_and_proton`;
- `metadata_reaction_selection.json` (SHA `abfc4a94…` unchanged);
- every other reaction, metabolite and gene;
- the earlier candidates.

Why the species guard matters: the 2-H⁺ count is correct only for the stored neutral ATP and ADP. If
they are ever migrated to charged microspecies, the step stops, and the proton count must be
reviewed again (the 1-H⁺ form would then be the balanced one).

## 3. Verification

### 3.1 Builds

| Build | Output | SHA-256 |
| --- | --- | --- |
| Control, no flag | `.scratch/r2041_forward_20261019/control_no_flag.xml` | `a56f3649…` (= baseline) |
| Candidate, `--r2041-forward` | `model/candidates/E5_coq9_resp_pkg_r132_r2041_forward_20261019.xml` | `da0c50ba…` |
| Control, no flag, final code | `.scratch/r2041_forward_20261019/control_no_flag_v2.xml` | see §3.6 |
| Candidate, final code | `.scratch/r2041_forward_20261019/candidate_v2.xml` | see §3.6 |

Build record (`.build.json`):
- head `285eaee`, `requested_build_complete` true;
- the selection is complete;
- `post_selection_r2041_direction` = applied, bounds [0.0, 1000.0];
- input and source unchanged during the build;
- 2306 reactions, 1875 species, 1073 genes, the same as a56f3649.

The control's record has no R2041 step. The `.coq9_genes.tsv` is identical to a56f3649's. Logs are in
[logs/](logs/), SHAs in [logs/sha256.txt](logs/sha256.txt), commands in [commands.txt](commands.txt).

### 3.2 Whole-model diff

[scripts/model_diff.py](scripts/model_diff.py) compares every reaction, metabolite, gene and group
field, the objective and the model notes. Output: [model_diff.json](model_diff.json). It finds 2
differences, R2041 `bounds` [−1000, 1000] → [0, 1000] and R2041 `notes` (6 added keys).

The raw XML diff ([logs/xml_diff.txt](logs/xml_diff.txt)) changes one line on `R_R2041`:
`reversible` true → false, and `fbc:lowerFluxBound` `cobra_default_lb` → `cobra_0_bound`. The
`cobra_0_bound` parameter already exists. The diff also adds six `<p>` note lines.

The audit repeated this with libSBML and found exactly 3 field differences, all on R2041: lb,
reversible and notes.

### 3.3 Chemistry

In the candidate, R2041 is ATP (C10H16N5O13P3, 0) + D-arabinose (C5H10O5, 0) → ADP (C10H15N5O10P2,
0) + A5P (C5H9O8P, −2) + 2 H⁺.

| | Left | Right |
| --- | --- | --- |
| H | 26 | 26 |
| Charge | 0 | 0 |

The element and charge residual is empty (verified, and by the audit with
`check_mass_balance()`). No other reaction was touched.

### 3.4 Tests

| Run | Passed | Skipped | Failed | Errors |
| --- | ---: | ---: | ---: | ---: |
| Before the change (clean `285eaee`) | 165 | 31 | 11 | 14 |
| After the change | 169 | 31 | 11 | 14 |

[logs/pytest_comparison.txt](logs/pytest_comparison.txt): the only differences are the 4 new tests,
absent → passed. The 11 failures and 14 errors are pre-existing; STATE.md recorded the same counts
on 2026-10-13. The audit found 0 outcome or message differences among the 221 common tests. The
post-audit rerun is in §3.6.

### 3.5 The 10 authorized LPs

Runner: [scripts/run_lps.py](scripts/run_lps.py). Raw output: [lp_results.json](lp_results.json), log
[logs/lps.log](logs/lps.log).

Context: the `coq_candidate_checks_20261005` config, PO1f SD-Leu overlay (profile `35307853…`, effect
`d15acbde…`), and Gurobi with 1 thread, Presolve 0, tolerances 1e-7 and 60 s. Closure is by
`close_model` (190 changes); the C_mi drain is the one in `run_controls.py`. 100 % WT means the
biomass lower bound is set to the WT optimum.

| LP | M0 = a56f3649 | M1 = `da0c50ba` |
| --- | --- | --- |
| WT growth (biomass_C) | 1.430805735352851 | 1.430805735352859 |
| R2041 min at 100 % WT | **−999.6618** | 0.0 |
| R2041 max at 100 % WT | 0.0 | 0.0 |
| Closed-exchange max ATP dissipation, C_cy | 0.0 (0 non-zero fluxes) | 0.0 (0) |
| Closed-exchange max ATP dissipation, C_mi | 0.0 (0) | 0.0 (0) |

All 10 are optimal, with exactly 10 optimizations. No stop condition fired.

**The M0 minimum is the internal loop from the review** (verified). At that solution these
reactions all run at magnitude 999.6618 and biomass is unchanged:
- R2042 +, R2190 −, R1964 −, R803 +, R714 −, R2041 −.

D-arabinose and arabinitol export (R1062, R1981) are 0. R714 sits at −1000: the loop's −999.66 plus
the WT flux −0.338. That bound caps the loop. The audit checked that the six-reaction flux change
sums to exactly 0 for every metabolite. So in a56f3649 reverse R2041 is usable only as a net-zero
cycle, and in M1 it is not usable at all.

**How this compares with the review.**
- Closed ATP 0 on M0 was **provable** in the review (from `0ac2addf` plus the R132 deletion). It is
  now **observed**.
- M1's closed maxima are 0 as expected: removing reverse flux cannot raise the maximum.
- "No growth effect" was **inferred**. It is now **observed** for WT under this context.

### 3.6 Re-verification after the audit fix

The audit fix touched code and tests only (§5). Both builds were repeated with the final code into
`.scratch/`, and the full suite was run again (verified):

- `control_no_flag_v2.xml` = `a56f3649…` and `candidate_v2.xml` = `da0c50ba…`, byte for byte.
  - The candidate in `model/candidates/` was therefore kept, not rebuilt.
  - Its `.build.json` records the pre-fix source SHA of `r2041_direction.py` (`1e921285…`).
  - The final module SHA is in [logs/sha256.txt](logs/sha256.txt), and the v2 build records are in
    `.scratch/`.
  - The fix does not reach the candidate's code path: its input is at `before_bounds`.
- Full suite: 169 passed, 31 skipped, 11 failed, 14 errors, 397 subtests
  ([logs/pytest_comparison_v2.txt](logs/pytest_comparison_v2.txt)). The only differences from the
  before-run are the same 4 new tests, absent → passed. The focused test now has 9 fail-closed
  subtests.
- The LPs were not re-run. They depend only on the model files, which are byte-identical.

## 4. Mechanism and evidence (from the review; tiers kept)

- **Thermodynamics (verified):** the
  [eQuilibrator](https://equilibrator.weizmann.ac.il/search?query=C00002+%2B+C00216+%3D+C00008+%2B+C01112)
  estimate is ΔrG'° = ΔrG'm = −24.7 ± 4.7 kJ/mol and K'eq = 2.1 × 10⁴ (pH 7.5, pMg 3, I = 0.25 M,
  298.15 K).
  - Reverse flux needs [A5P]/[D-arabinose] > K'eq × [ATP]/[ADP]. At the least favourable 2σ value
    that is > 479 for ATP/ADP ≥ 1 (derived).
  - This is an estimate, not a measurement. *Y. lipolytica* metabolite ratios were not looked up.
- **Enzymology:**
  - **Reaction definition (verified):** [ENZYME 2.7.1.54](https://enzyme.expasy.org/EC/2.7.1.54) and
    [RHEA:24588](https://www.rhea-db.org/rhea/24588) give D-arabinose + ATP = D-arabinose 5-phosphate
    + ADP + H⁺. Rhea links 0 UniProt enzymes.
  - **Literature (unread):** Volk 1962, *Propionibacterium pentosaceum*.
  - **Organism (verified, model):** no *Y. lipolytica* enzyme is known, and R2041 has no GPR.
- **ATP loop named by the builder patch:**
  - **Not found statically (verified):** A5P has only R2041 and R2042; there is no phosphatase.
  - **No ATP gain (observed):** closed ATP is 0 with and without reverse R2041.

The direction change therefore rests on thermodynamics and the kinase definition, not on a loop. It
says nothing about whether R2041 or R2042 exist in *Y. lipolytica*.

## 5. Independent audit

A fresh read-only reviewer received TASK.md, the diffs, the logs and the deliverables. It made no
repository writes. Verdicts:

| Claim | Verdict |
| --- | --- |
| C1 control reproduces a56f3649 | supported |
| C2 diff limited to R2041 bounds and notes | supported (raw libSBML and cobra comparisons) |
| C3 R2041 balanced | supported |
| C4 guard fails closed, bounds/notes only, idempotent, round-trip safe, opt-in | supported (37 extra probes) |
| C5 rejects the 1-H⁺ patch form and migrated ATP | supported |
| C6 the 10 LP numbers and the budget | supported (from lp_results.json; LPs not re-run) |
| C7 LP protocol matches `run_controls.py` | supported, differences listed below |
| C8 pytest unchanged apart from the new tests | supported |
| C9 evidence text matches ENZYME, Rhea and eQuilibrator | supported (re-fetched) |
| C10 no undeclared file changed | supported |

**C7 protocol differences.** None changes the comparison:
- 100 % WT is a biomass lower bound, not cobra's FVA constraint. The reference has no such LP.
- `raise_error=False` is the cobra default.
- The runner has extra identity asserts.
- The reference's R385 control was not run (out of scope).

How each other finding was handled:

| # | Finding | Action |
| --- | --- | --- |
| 3 | The guard accepted an R2041 already at [0, 1000] without this curation's notes and then wrote "previous bounds [−1000, 1000]", which would be false. R1931 and R1159 have the same pattern | **Fixed:** now fails closed, with a new subtest (§3.6) |
| 7 | Test gaps: compartment guard, forward-only without notes | **Fixed:** two subtests added. Flag-off byte identity is covered by the control builds |
| 9 | `model/curation/README.md` did not document the step | **Fixed** |
| 1 | The species guard sees the state at its stage only. A later stage that migrated ATP/ADP would not be caught | Reported. In this build the final species are unchanged (C2/C3) |
| 2 | Only the six curated note keys are guarded. The `metadata_*` notes are not pinned | Reported. The exact 2-H⁺ stoichiometry already identifies the selection state |
| 4 | R2041 still carries `metadata_selected_fields: bounds; stoichiometry`, but its bounds now come from this step | Reported. The `r2041_direction_*` notes say so. Rewording any note changes the model SHA after the LPs, so it was not done |
| 5 | The note "EC 2.7.1.54 kinase direction (ENZYME, Rhea)" overstates: both sources write "=", which defines the reaction, not its irreversibility | Reported, same reason. The limits note states that the curation "rests on thermodynamics". Proposed wording for a later note-only change: "EC 2.7.1.54 written as the kinase reaction (ENZYME, Rhea)" |
| 6 | `run_lps.py` records rather than enforces the TASK.md stop conditions | Reported. Neither condition fired (WT relative difference 5.6e-15; closed objectives 0) |
| 8 | The flag has no "requires" rule. The LP evidence covers only the E5/CoQ9/R132 chain | Reported (§6) |
| 10 | REPORT.md and STATE.md not yet written | Done now |
| 11 | Pre-existing in baseline and candidate alike: 89 libSBML read errors, the R304 notes error, R2042's mixed protonation, R2190's EC/product mismatch | Reported. Not part of this change |
| 12 | Files of another task in the research root changed during the audit | Not this task. Nothing was written there |
| 13 | WT non-zero flux counts differ (448 vs 452) at equal objective | An alternate optimum. Not evidence of a behaviour change |

## 6. Limits, what was not tested, and next gates

- **Not tested:**
  - single-gene KO screens (the optional 1,074-LP screen was not authorized);
  - knockout-specific use of reverse R2041;
  - builds without the E5/CoQ9/R132 chain (the flag works with any build; the guard keeps it exact);
  - the plain model without the PO1f overlay (optional LPs P1–P2).
- **Other candidates:** c98961e9 and 51c6fabd (uncommitted, other worktrees) still carry reversible
  R2041. Stacking this step on the GPR fixes waits for that work to be merged.
- **Note wording** (audit #4 and #5) can be corrected in a later note-only change. That would give a
  new SHA, and the LPs would need to be re-run or argued equal.
- **Gates that stay:**
  - the ATP/ADP microspecies migration (the guard stops this step);
  - the existence of R2041 and R2042 in *Y. lipolytica* (orphan line);
  - R2190's EC 1.1.1.250 (D-ribulose-forming) against its D-arabinose product, reviewed separately.
    R2190 is what lets the M0 loop exist.
- **Publishing:** nothing committed or pushed.

**Could the method be improved?**
- A loopless FVA (or an FVA with the internal-cycle check) would have shown, without a static trace,
  that reverse R2041 is only ever a null-space cycle.
- Test: on M0, a loopless min of R2041 at 100 % WT should be 0. On M1 it is 0 by bounds. That costs
  2 MILPs or 2 LPs with the `add_loopless` constraints, which need authorization.
- The runner could enforce the stop conditions as code (audit #6).

**Statement levels:**

| Level | Statements |
| --- | --- |
| Provable | Closed ATP ≤ 0 on M1 given 0 on M0; M1 WT ≤ M0 WT |
| Statically verified | diff, balance, guards, the loop's zero sum |
| Observed in runs | the 10 LP values |
| Untested | KO effects, other build chains |
