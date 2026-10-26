# Lipoylation module + acyl-ACP chemistry on the main candidate 2c9cebd0 (2026-10-23)

Scope, authorization and acceptance criteria: [TASK.md](TASK.md). Exact commands:
[commands.txt](commands.txt). Audit: [AUDIT.md](AUDIT.md).

Control, seed-1 and `b3b9fe22` XMLs (with `SHA256SUMS.txt`) are in
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/lipoylation_acylacp_on_2c9cebd0_20261023/`.

Nothing is committed or pushed.

## 1. Outcome

- **New candidate:**
  `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_20261023.xml`
  - SHA-256 `a12378210c7c66720d99fb7bf71d14cb7538ea6a3cf3e08919c85fb3c9c47d82`.
  - It is `2c9cebd0` (main, PR #18) plus `--lipoylation-module --acyl-acp-chemistry`, the two opt-in steps
    from the unmerged commit `9f2356b`.
  - No reference model is designated; this is a candidate.
- **The diff against `2c9cebd0` is exactly the two steps:** 135 field-level differences. They are the
  same 135 entries, entity by entity, as `b3b9fe22` against `164f235f` (the two steps on their original
  base).
- **Knockouts** (Task A medium + PO1f, Gurobi 1e-9): 13 genes become essential at every cutoff, and
  none moves the other way. They are:
  - the 3 genes the lipoylation module adds;
  - 10 of the 11 shared genes of the `b3b9fe22` screen.

  The 11th, YALI1E07744g (transketolase), was already essential in `2c9cebd0`, as predicted in
  TASK.md.

## 2. What changed in the code, and why

| File | Change |
| --- | --- |
| `platform/scripts/gem_annotate/lipoylation_module.py`, `acyl_acp_chemistry.py` | Copied unchanged from `9f2356b`; SHA-256s match its `logs/sha256.txt` (`11d4820d…`, `d27dff5f…`). |
| `model/curation/lipoylation/lipoylation_module.json`, `model/curation/chemistry/acyl_acp_chemistry.json` | Copied unchanged (`c51b3b14…`, `b6b88f10…`). |
| `model/reports/lipoylation_module_20261020/`, `acyl_acp_block_review_20261020/` | Copied unchanged: the evidence trail the data files cite. Not copied: `9f2356b`'s candidate XMLs, build JSONs, README rows and STATE entries. |
| `cli.py` | The two flags, with the same help text as `9f2356b`. |
| `main.py` | The two parameters, the rule that `--acyl-acp-chemistry` requires `--lipoylation-module`, and the two calls, placed **after the R2193 step**. |
| `tests/test_lipoylation_module.py` | Copied unchanged. Its baseline `164f235f` is on this branch. |
| `tests/test_acyl_acp_chemistry.py` | Its baseline `fdd00ac3` is not on this branch, so it is rebuilt in memory (`164f235f` + the lipoylation step) and round-tripped through SBML. Test bodies are unchanged. |
| `tests/test_lipoylation_acylacp_on_2c9cebd0.py` (new) | Both steps on `2c9cebd0` change only their targets and add only the 3 genes. Both reruns are no-ops after an SBML round trip. The heme/iron guard rejects a lipoylated `biomass_C`. |

**Why the steps moved in `main.py`** (verified by code reading and by the third new test):
- `apply_heme_iron_biomass` accepts only two reviewed `biomass_C` states: the exact composition, bounds
  and note keys of the base, or the base after the TPP requirement. A `biomass_C` that already carries
  the four lipoyl terms and their note fails that check.
- `apply_lipoylation_module` checks only its own four `biomass_C` coefficients.
- So the lipoylation step must come after the heme/iron step. `apply_tpp_requirement` also pins the
  full `biomass_C`; that is satisfied because the TPP steps run before heme.
  - Running after R2193 is not required (its guard reads only R2193); it was chosen as the last slot.
  - The acyl-ACP step requires the lipoylation step (ACP1 = HSR), so it comes last.
- In the planning turn, both `9f2356b` functions, applied in memory to the saved `2c9cebd0`, returned
  `applied`. The real build log shows the same order and statuses: heme `applied`, R2193 `applied`,
  lipoylation `applied`, acyl-ACP `applied`.

**Limit of this order** (verified by the third test): the heme/iron step cannot be rerun as a no-op
on the finished candidate, because its guard sees the lipoyl terms.
- The TPP step already fails to rerun on `2c9cebd0` itself (`m975[C_cy] identity differs`, found by
  the audit). That predates this change.
- The other steps return `already_correct` on the candidate (audit): thiamine-P, proton lock, R2193,
  lipoylation and acyl-ACP. This matters only if a later build reads this candidate as its input instead of rebuilding
from source.

## 3. The science change (unchanged decisions, re-applied)

Both decisions were made earlier and are not re-argued here; see their reports.
- **Lipoylation module:** [lipoylation_module_20261020/REPORT.md](../lipoylation_module_20261020/REPORT.md), authorized 2026-10-20.
  - Adds R2286 (LIP2-type octanoyl transfer, YALI1F34729g, tier B weak), R2287 (lipoyl synthase LIP5,
    YALI1E13472g, tier B) and R2288 (LIP3-dependent E2 lipoylation, `YALI1C04991g and YALI1E00574g`,
    tier C, mechanism unresolved).
  - Adds two labelled assumptions: R2289, cysteine transport into the mitochondrion, and R2290, a
    5'-deoxyadenosine sink.
  - Fills the ACP1 formula (HSR) and the lipoyl-H-protein formula, removes the H+ in R362/R364, and
    adds a provisional lipoyl growth-dilution term at α = 1e-5 mmol/gDW. α is user-selected; no
    *Y. lipolytica* value was found.
- **Acyl-ACP chemistry:** [acyl_acp_block_review_20261020/REPORT.md](../acyl_acp_block_review_20261020/REPORT.md).
  - Fills 32 acyl-ACP/ACP/fatty-acid formulas, removes one surplus H+ from 19 reactions, and makes 7
    cytosolic acyl-ACP hydrolases forward-only.
  - This blocks R2216–R2241 in every medium. The user accepted that on 2026-10-21 (recalled from
    `9f2356b`'s STATE entry).

## 4. Verification

| Check | Result | Level |
| --- | --- | --- |
| Full test suite, before vs after ([before](logs/pytest_before.log), [after](logs/pytest_after.log), outcome lists `.outcomes`) | Identical test by test; 13 new passes (5 lipoylation, 5 acyl-ACP, 3 new). Totals: 11 failed / 251 passed / 18 errors before, and 11 / 264 / 18 after. All failures and errors predate this change. | verified |
| Focused tests ([log](logs/pytest_focused.log)) | 18 passed, including the 5 heme/iron tests | verified |
| Control build, without the new flags ([log](logs/build_control_no_flag.log)) | `2c9cebd0a5c5…` = the main candidate, byte for byte | verified |
| Candidate under `PYTHONHASHSEED` 0 and 1 ([sha256.txt](logs/sha256.txt)) | Both `a12378210c7c…` | verified |
| Whole-model diff vs `2c9cebd0` ([model_diff_vs_2c9cebd0.json](model_diff_vs_2c9cebd0.json)) | See the breakdown after this table. | verified |
| Same entities as `b3b9fe22` vs `164f235f` ([diff_comparison.json](diff_comparison.json)) | 0 entries only in one diff. Every before and after value is equal except `biomass_C`. The two bases start from different `biomass_C` notes and stoichiometry (the TPP, heme a and iron terms of `2c9cebd0`). The changes themselves are identical: four lipoyl terms at ±1e-5, plus one added note, `lipoylation_growth_dilution`, with the same value. The audit's own XML diff agrees. | verified |
| Balance of every touched reaction ([balance_check.json](balance_check.json)) | See the breakdown after this table. | verified |
| WT, Task A medium + PO1f, Gurobi 1e-9 ([solves.json](solves.json)) | `2c9cebd0` 0.312426291488; candidate 0.312420852777 (Δ −5.44e-6, close to `b3b9fe22`'s −5.4e-6 against `a56f3649`). Both optimal, and both pass the check. In pFBA the lipoylation route carries flux (R2286 and R2287 at 6.25e-6, R2290 at 1.25e-5). | verified |

**Diff breakdown.** The 135 differences against `2c9cebd0`:
- 5 reactions added (R2286–R2290), 6 species added, 3 genes added (YALI1C04991g, YALI1E13472g,
  YALI1F34729g);
- 34 formula changes with their notes (32 acyl-ACP fills, ACP1 and lipoyl-H-protein) and 1 name change
  (the lipoyl-H-protein);
- 22 stoichiometry changes with their notes (19 acyl-ACP reactions, R362, R364, `biomass_C`) and 7
  bound changes (the hydrolases);
- 1 gene-reaction membership change (YALI1E00574g gains R2288).

Nothing else differs.

**Balance breakdown.**
- 19 reactions are fully balanced.
- 6 reactions (R2218, R2221, R2224, R2227, R2230, R2233) keep a charge-only residual of −2 or −3. That
  equals the curated `residuals_after` in `acyl_acp_chemistry.json`.
- R2290 is the one-sided 5'-deoxyadenosine sink, by design.
- `biomass_C` is a pseudo-reaction and was not checked.

## 5. Knockouts

[ko_compare.tsv](ko_compare.tsv), [ko_compare_summary.json](ko_compare_summary.json),
[log](logs/ko_compare.log), [script](scripts/ko_compare_two.py) (a two-model copy of the 2026-10-22
three-model script, using the same `checked_solve.py`).

- **Protocol:**
  - Gurobi 13.0.1, FeasibilityTol and OptimalityTol 1e-9, 60 s per LP;
  - the model is reset before each knockout;
  - targeted check: knocked-out reactions |v| ≤ 1e-9, other bound or |S v| violations ≤ 1e-6.
- **Solver health** (verified):
  - 2,161 LPs; no point fails the check and no fresh re-solve was needed;
  - the largest |S v| is 3.0e-12 (`2c9cebd0`) and 1.2e-12 (candidate); knocked-out flux and bound
    violations are 0;
  - the same 20 genes are infeasible in both models (16 in scope, 4 mtDNA). They are reported as
    unresolved, which is the solver's verdict, not a stale point.
- **Reproduction:** the baseline score at 15 % (93 / 220 / 36 / 697 / 16, MCC 0.345) equals the value
  recorded for `2c9cebd0` on 2026-10-22.

**Class changes against `2c9cebd0`** (identical at 1, 5, 10 and 15 %). Every gene goes to KO/WT
= 0.0 (optimal solve, exactly zero growth). The KO/WT values in the `2c9cebd0` column equal those the
`b3b9fe22` screen lists for `a56f3649`.

Names are model annotations, not verified in this task. Reaction assignments are model GPRs, not
experimental validation.

| Gene | On list | `2c9cebd0` KO/WT | Model reactions (candidate) | Function (source; evidence) |
| --- | --- | ---: | --- | --- |
| YALI1C04991g | yes | absent | R2288 | LIP3 ortholog, lipoate-protein ligase family (curation file; tier C) |
| YALI1E13472g | yes | absent | R2287 | lipoyl synthase, LIP5 ortholog (curation file; tier B) |
| YALI1F34729g | no | absent | R2286 | LIP2-type octanoyltransferase candidate (curation file; tier B weak) |
| YALI1A20089g | yes | 1.000 | R78, R1400–R1402 | mtFAS 3-hydroxyacyl-ACP dehydratase (model reaction names; no name verified) |
| YALI1C07638g | yes | 1.000 | R1434, R1435, R1869 | acyl-CoA desaturase: Δ9 on R1434/R1435; R1869 is C18:1 → C18:2 (model reaction names; no name verified) |
| YALI1C15991g | yes | 1.000 | R87, R88, R1393 | acetyl-CoA carboxylase (model reaction names; no name verified) |
| YALI1C26939g | no | 1.000 | R78, R1403–R1405 | mtFAS enoyl-ACP reductase (model reaction names; no name verified) |
| YALI1D18037g | yes | 0.758 | R1394–R1396, R2062, R2218–R2233, … | β-ketoacyl-ACP synthase; also on R2062 (NADH:ubiquinone oxidoreductase) (model; no name verified) |
| YALI1D26360g | yes | 0.974 | R357, R360, R361, R364, R365, R1968–R1976, R2281, … | 2-oxo-acid dehydrogenase / glycine-cleavage subunit (model reaction names; no name verified) |
| YALI1E00574g | yes | 0.983 | R2288, R357, R360, R361, R365 | glycine-cleavage / 2-oxo-acid dehydrogenase subunit, AND partner in R2288 (curation file context gene; no name verified) |
| YALI1E22262g | yes | 1.000 | R78, R539 | acetyl/malonyl-CoA:ACP transacylase (model reaction names; no name verified) |
| YALI1F37498g | no | 1.000 | R78, R1397–R1399 | mtFAS 3-oxoacyl-ACP reductase (model reaction names; no name verified) |
| YALI1F38317g | yes | 1.000 | R1394–R1396, R2218–R2233, … | β-ketoacyl-ACP synthase (model reaction names; no name verified) |

YALI1E07744g (transketolase, listed) is essential (KO/WT 0.0) in both models, so its call does not
change.

**Scores** (model genes not on the experimental list are the user-defined negative class; mtDNA genes
out of scope; the list was used in development, so these are a development reference, not independent
validation):

| Cutoff | `2c9cebd0` TP / FN / FP / TN / unresolved (MCC) | Candidate TP / FN / FP / TN / unresolved (MCC) |
| --- | --- | --- |
| 1 % | 75 / 238 / 29 / 704 / 16 (0.306) | 85 / 230 / 32 / 702 / 16 (0.329) |
| 5 % | 81 / 232 / 31 / 702 / 16 (0.321) | 91 / 224 / 34 / 700 / 16 (0.343) |
| 10 % | 85 / 228 / 36 / 697 / 16 (0.318) | 95 / 220 / 39 / 695 / 16 (0.341) |
| **15 %** | **93 / 220 / 36 / 697 / 16 (0.345)** | **103 / 212 / 39 / 695 / 16 (0.367)** |

- **Positive coverage:** 324 listed genes are in `2c9cebd0`, and 326 in the candidate (+YALI1C04991g,
  +YALI1E13472g). The 2 new positives are TP; YALI1F34729g, unlisted, is a new FP.
- **Arithmetic at 15 %:** 8 shared listed genes FN→TP and 2 shared unlisted TN→FP; plus 2 new TP and 1
  new FP. That gives 93+10 / 220−8 / 36+3 / 697−2.
- **Not a reason to accept anything:** these gains do not show the change is right. Per AGENTS.md the
  decisions rest on the 2026-10-20/21 evidence, not on recall.

All other knockouts keep their class. Their KO/WT moves by at most 1.6e-4 (YALI1E19603g). In
absolute growth, 921 knockouts differ by more than 1e-6, mostly because the WT itself moves by
5.4e-6.

## 6. Mechanism of the class changes (inferred, except for one gene)

- **Verified, one gene (audit):** with YALI1F34729g knocked out, the candidate grows at 0.0.
  - The same knockout with the four lipoyl `biomass_C` terms removed in memory grows at 0.312426291488,
    exactly the `2c9cebd0` WT.
  - So that call, and the whole −5.44e-6 WT cost, come from the provisional α = 1e-5 lipoyl
    requirement. The acyl-ACP block has no effect on WT.
- **The other 12:** inferred from the reactions each gene carries; no targeted control was run.
- **Lipoylation genes and mtFAS** (YALI1A20089g, C26939g, D18037g, E22262g, F37498g, F38317g):
  - The new lipoyl growth-dilution term needs octanoyl-ACP from mtFAS (R78, R539, R1394–R1405) and the
    new R2286–R2288.
  - Before, no biomass term depended on lipoylation.
- **E3/H-protein-type subunits** (D26360g, E00574g): they carry the lipoyl-dependent dehydrogenases
  and glycine cleavage. E00574g is also an AND member of R2288.
- **Acetyl-CoA carboxylase (C15991g) and the Δ9 desaturase (C07638g):** these are the most likely
  effect of the acyl-ACP step.
  - It blocks R2216–R2241, the iMM904-derived acyl-ACP route with placeholder formulas, and makes its
    hydrolases forward-only.
  - Malonyl-CoA and unsaturated acyl chains then appear to have one route each.
  - This is inferred from reaction membership only. The `b3b9fe22` screen did not trace it either.

## 7. Independent audit

See [AUDIT.md](AUDIT.md). Verdicts and responses are summarised below.

- **Verdicts:** C1 and C3–C8 supported; C2 partially supported.
- **The auditor's own checks:** an element-level XML diff, a rebuild of `fdd00ac3` that is
  byte-identical, WT re-solves, and a recount of the scores from the experimental list.

| Finding | Response |
| --- | --- |
| C2: running "after R2193" is not necessary; "after heme/iron" is. The TPP guard also pins `biomass_C`. | Agreed. The binding constraint is "after the heme/iron step". The TPP guard is satisfied because TPP runs before heme. Wording corrected in §2. The steps stay after R2193, which is also valid. |
| C3: `pytest_after` overlapped the builds and the sweep. | Correct. They were started in parallel; commands.txt now says so. Outcomes are unaffected. |
| C4/5: control, seed-1 and `b3b9fe22` XMLs only in the scratchpad. | Copied, with `SHA256SUMS.txt`, to `$IYALI26_RESEARCH_ROOT/artifacts/tasks/lipoylation_acylacp_on_2c9cebd0_20261023/`. |
| C5: `biomass_C` after-values also differ by the base's pre-existing stoichiometry. | Correct (TPP, heme a and iron terms of `2c9cebd0`). The deltas are identical. Wording corrected in §4. |
| 1: the YALI1F34729g call and the whole WT cost come from α = 1e-5. | Adopted as a verified result for that gene (§6). |
| 2: the TPP step already fails to rerun on `2c9cebd0`. | Adopted. §2 now says the TPP rerun failure predates this change; only the heme/iron rerun failure is new. |
| 3: test quality. | My new test now compares the union of keys and checks that no gene is removed. The heme test matches the error message. Rerun: 3 passed ([log](logs/pytest_new_test_after_audit.log)). Tests are not in the build's source hashes, so the candidate is unaffected. The builder order is still enforced only by `main.py` (§9). The lipoylation rerun note check is `9f2356b` code, left as is and noted. |
| 4: gene labels from model reaction names. | The README row and STATE now say "genes on … reactions (model names; not verified)". The §5 table already marks names as not verified; YALI1C07638g's row now includes R1869. |
| 6: cosmetic issues inherited from `9f2356b`. | Noted, unchanged (the after-values equal `b3b9fe22`). |

## 8. Limits and what was not tested

- Only WT and single-gene knockouts on the Task A medium + PO1f. Not tested: SD-Leu, other media, FVA,
  blocked-reaction runs, the respiratory overlay.
- The α = 1e-5 lipoyl coefficient is provisional, and it drives at least the YALI1F34729g call (§6).
  Sensitivity at 1e-6 and 1e-4 was run on `fdd00ac3` (recalled from its curation record), not on this
  candidate. The tolerance concerns about tiny biomass terms (STATE, open problem 1) apply
  here too. The Gurobi 1e-9 sweep shows no solver artifact.
- The open problems recorded for `2c9cebd0` (iron chemistry, siroheme, R863, GPI items and so on)
  stand unchanged. So do the open items of the two steps' own reports: the R2288 mechanism, the R2289
  and R2290 assumptions, and the tier B/C gene identities.
- Mechanisms in §6 are inferred, except YALI1F34729g's. Targeted controls (WT, KO, the bypass closed, and combined) would be
  the next gate if a class change needs to be explained.
- Not committed or pushed. `9f2356b` remains unmerged; this task copies its files, it does not merge
  that branch.

## 9. Could the method be improved?

- **The order constraint** is enforced only by the call order in `main.py` and by one test. A more
  robust design is one shared `biomass_C` before-state registry for the TPP, heme and lipoylation
  steps. Each guard would then accept "the base plus any reviewed subset", and reruns on finished
  models could be no-ops.
  - Test: rerun every step on the finished candidate and expect `already_correct` everywhere.
- **The diff comparison** (`scripts/compare_diffs.py`) is a cheap regression check for re-basing an
  opt-in step. It could become a standard tool when a step moves to a new base.
