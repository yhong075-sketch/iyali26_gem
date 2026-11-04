# AUDIT: mcnaugh_pass reactions merge (2026-11-03)

**The auditor.** A fresh general-purpose agent following `.claude/agents/code-reviewer.md`, read only, with no
conversation history. It was given TASK.md, the diff, the curated data, candidates v2 and v3, Andrew's model, the
artifacts and the test outputs.

**What it ran.** 9 LPs (8 on v3, 1 on v2), with Gurobi 1e-9 checked solves on Task A + PO1f. No knockout sweeps, no
builds.

**Its inputs.** It checked them first:
- `mv3.xml` = `ec7e3465…1022`;
- `mv2.xml` = `9158e9f2…4129`;
- Andrew's model = `e74a0f50…331a`;
- the source and data SHAs in `build.json` match the current files.

Its scratch outputs are in this session's scratchpad (`scratchpad/review/`) and are not kept. Despite the read-only
instruction, it also wrote `platform/lps_out.json` (2026-11-03 22:11: its WT, energy and biomass-objective LPs on
v3). The file is untracked, left in place, and not part of this change.

Each finding below was re-checked in this session after the audit. "Re-checked" means I ran the check myself.

## Claims

| # | Claim | Verdict | Evidence |
| --- | --- | --- | --- |
| 1 | 17 reactions, each Andrew's reaction term for term in our species, with his bounds and GPRs; 4 closures kept at 0..0 | supported, except the subsystem (F1) | Every row maps through the `facts.json` species map, with no merged species. The v3 file matches the JSON. The closures changed only in their bounds headers (and `reversible` for R1771). |
| 2 | The species mapping is correct; no new species or genes | supported | Names and compartments match. Species and gene lists are byte-identical to v2. The biomass components are the same species at the same coefficients as biomass_C. |
| 3 | The step fails closed, is idempotent, opt-in, survives a round trip and runs before the crosswalk step | partially supported | 7 tests pass and the probes confirm every tested guard. It does not fail closed on rewritten metadata (F2), and the subsystem is lost on write (F1). |
| 4 | v3 differs from v2 only by the 17 additions and the 4 closures | supported, byte level | Removing the 17 blocks and restoring the 4 closure headers gives a file byte-identical to v2. |
| 5 | WT 0.31248856, energy 0, all 17 can carry flux, growth on the two biomass reactions | supported | WT, energy and both biomass objectives were reproduced. One LP with all 15 lipid steps capped at 1 has every step carrying flux. The 34-LP FVA was not repeated. |
| 6 | The screen differs from v2 only for OLE1 and YALI1D22124g; the scores; pools are not the OLE1 cause | supported | The TSVs were recomputed gene by gene. With forward-only pools, OLE1 KO = WT = 0.3120960274. The generic-product cause is inferred: no positive test was run. Wording fix: YALI1D22124g went from model-essential to non-essential. |
| 7 | The 12 not added are duplicates or reversals; none of the 17 duplicates an existing reaction | partially supported | ALCD2irm = source R132 and G3PD1irm = source R349 (re-checked: neither is in v2 or v3). Of the 10 duplicates, 7 are exact. 3HPCOAHYD, COMT2 and COMT3 have one extra cytosolic H+ (m10) (re-checked). 4ABUTtm and ABUTt2r match only with 4abut_c → m133, while Andrew's crosswalk gives m2040 (re-checked). None of the 17 shares its metabolite set with an existing reaction. |
| 8 | The exceptions are described accurately | partially supported | The +4 charge from m1735 is exact, and so is +0.139 %. "R-group imbalance" is exact only in Andrew's species; in ours the numbers come from missing or defective formulas (F4). |

The auditor wrote "9 of the 10 duplicates are exact"; it lists 7, and the re-check gives 7 exact and 3 that differ by
H+.

## Findings, most severe first

| # | Severity | Where | Finding | Re-checked | Status |
| --- | --- | --- | --- | --- | --- |
| F1 | low–medium | `mcnaugh_reactions.py` (reaction construction) | `subsystem=` is set, but the reactions are never added to the existing group "Glycerolipid metabolism" (`group33`). cobra writes subsystems only through groups, so in v3 all 15 lipid steps read back with subsystem `''`, while R353 and R1771 are in `group33`. Subsystem-based maps and audits miss them. | yes, on v3 | closed in v4 (fix review) |
| F2 | low–medium | `mcnaugh_reactions.py` `_definition` | The already-applied check compares stoichiometry, bounds and GPR only. A model with the new reactions present but different notes, name or annotation returns `already_correct`. TASK acceptance 1 requires failing closed on a rewritten state. | yes, from the code | closed in v4 (fix review) |
| F3 | low | curated JSON, BIOMASS_C_core / _full `composition_note` | "Unlike biomass_C it has no ceramide pool" understates the difference. Both also lack 13 cofactor terms, the m1401 biomass product and the m2062/m2064 by-products, and they use charged tRNAs instead of the `trna_biomass_residue_*` species. | yes, term by term on v3 (18 and 17 missing terms) | closed in v4: the generator now derives the notes from the comparison (fix review) |
| F4 | low | `verify_v3.json`, REPORT §4 | Not all of the imbalance comes from the R group. G3PAT shows the whole product missing (m577 has no formula). DGAT shows O +56 because m1640's formula `C6H5O6*3` parses as O63. Pre-existing `*` formulas: m1640, m1123 `C6H10O5*2`, m1648 `…P*2`. | yes (m1640 parses as `{'C': 6, 'H': 5, 'O': 63}`) | REPORT corrected; the formulas belong to the charges/formula line |
| F5 | low | `test_mcnaugh_reactions.py` | `network()` ignores subsystem, groups, names, notes and annotations, so the round-trip test passes despite F1. There is no rewritten-metadata test (F2), and notes are not checked for `&`. The curated JSON has no `&`. | yes, from the code | closed in v4 (fix review, with mutation checks) |
| F6 | low | guard | A reaction with the same stoichiometry under another ID is not detected: a probe added `OTHER_G3PAT160`, and the step applied and created a duplicate. | not re-run | partly closed in v4; fully closed in the final code (guard review) |
| F7 | observation, inferred | science | v3 YALI1D22124g KO growth (0.3120960274432655) equals v3 WT with forward-only pools (0.31209602744326415). The +0.139 % artefact probably runs through the lipid-particle fatty-acid–CoA ligase, and its FP → TN flip is probably the same composition escape as OLE1. | numbers match the session's own runs | recorded in REPORT §4–5 |

No finding is high severity.

## Tests (auditor's comparison, matching this session's)

- **v3 against the before run:** 24 tests added (17 from the committed annotation round, 7 mcnaugh), all passing.
  No status changes, no tests removed.
- **Failures and errors:** the same 3 failures (`test_dfba_new_fn_essentiality`, `test_lipid_moiety_ledger`,
  `test_vacuole_supply`) and 23 errors (4 uncollectable `gem_review_lit` modules, 14 in `test_lp_sn12_candidate`,
  5 in `test_r1726_cpt_open`) in both runs.
- **v3 against the v2 run:** only the 7 mcnaugh tests differ.

## Not checked by the auditor

- the full FVA (one sufficiency LP instead);
- the biology in `review_output.json`, the 26 look-alike pairs and the gene identities;
- that a build without the flag reproduces v2;
- re-running either screen, including the BIOMASS_C_core screen.

## Fix review (after "Fix and rebuild v4")

**Setup.** A workflow ran in this session: one fresh read-only reviewer on the fix diff and v4, plus one skeptic per
finding. There were no LPs and no builds. The diffs are kept in `artifacts/superseded_ec7e3465/fix_diff_*.patch`.

**Verdicts on the audit findings:**
- **F1, F2, F3 and F5: closed.**
  - All 15 lipid steps read back in `group33` with the group's subsystem.
  - The already-applied check covers name, notes, annotation and groups, with no false positive after an SBML round
    trip.
  - The biomass notes match a term-by-term comparison, including the pairing of each residue term with its charged
    tRNA.
  - Mutants without each fix make at least one new test fail.
- **F6: partly closed.** The before-state duplicate is caught, but two holes remained. Both were reproduced by a
  skeptic, and both are low severity:
  - **R1, order dependence:** the stoichiometry map kept one ID per key. On an already-applied model, a duplicate
    earlier in reaction order than the new reaction was missed (`already_correct`), while one later was caught.
  - **R2, reverse direction:** a duplicate written in the reverse direction was not detected.
- **v4 against v3:** an independent raw text diff confirms that only the 2 composition notes and 15 group-member
  lines in `group33` differ.

## Guard review (after "Fix and rebuild")

**Setup.** A workflow ran in this session: one fresh read-only reviewer on the guard diff, plus one skeptic. There
were no LPs and no builds. The diffs are kept in `artifacts/superseded_guard_b9473e58/guard_diff_*.patch`.

**Verdicts:**
- **R1 and R2: closed.** The map keeps every ID per stoichiometry, and the reverse key is checked. The reviewer's
  probes on an applied model:
  - a duplicate at index 0 or index 1000 now raises, where the old code returned `already_correct`;
  - a reverse duplicate raises in both states, where the old code applied or returned `already_correct`.
- **No false positives:**
  - none of the 17 rows matches any of the 2,329 baseline reactions, forward or reversed;
  - no baseline reaction shares a metabolite set with a row;
  - v4 (`e584151f`) returns `already_correct` with the new code, so v4 holds no duplicate that the old guard missed.
- **Keys:** no zero, negative-zero or non-finite coefficients, and integer and float coefficients match.
- **Tests:** the two new cases fail under the old code and pass under the new.

**New finding:**

| # | Severity | Finding | Verifier | Status |
| --- | --- | --- | --- | --- |
| G1 | low, latent | Rows are not compared with each other, only with reactions already in the model. A future curated file in which one row copies or reverses another would apply on the first run and then fail on the second. The current data (`a62050e3`) has no such pair. `test_second_run_is_already_correct` and the round-trip test would catch one at test time, not at build time. | reproduced | open: recorded, not fixed (no effect on the current data or candidate) |
