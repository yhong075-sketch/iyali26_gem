# Independent audit: R2171 written 1:1, 9d655c2b → 965ab947 (2026-10-28)

**Auditor:** a fresh read-only general-purpose agent with no conversation history. It was given TASK.md, the claims
and the deliverables, not the reasoning.

**Method:**
- shasum and a raw XML diff;
- every evidence source fetched again (KEGG, ModelSEED, UniProt, the S2 locus map);
- 27 Gurobi LPs on in-memory models (WT, INO1-null growth and demands, the R1720 variant, phosphate-closed Ins1P,
  the inositol-transporter knockout, no-inositol pFBA, no-uptake, six maximum-flux LPs);
- the knockout TSV recomputed; the pytest logs diffed per test; the focused test run once.

It changed no file.

| Claim | Verdict | Auditor's evidence (summary) |
| --- | --- | --- |
| K1 SHAs, control, seeds, round 2 | supported | Baseline and control `9d655c2b`; candidate, seed 1 and both round-2 builds `965ab947`; seed on line 1 of each log; control argv = baseline argv apart from argv[0] and the output path. |
| K2 diff = R2171 only | supported | Raw XML diff: 16 changed lines, all inside `R_R2171` (6 notes, 3 stoichiometry lines); matches the field diff. |
| K3 equation, charge, notes | supported | `m1393 + m32 --> 2 m10 + m1994 + m813`, charge 0 (was −100); no `<`, `>` or `&` in the notes. |
| K4 evidence | partially supported | R03332 1:1, rxn06362 1:1, the EC 4.6.1.13 "bacterial" comment, K01771's only *Y. lipolytica* member, Q6C9J1 (TrEMBL, score 1, PI-PLC X) and the locus map all confirmed; no inositol cyclic phosphate species found. Gap: E6 is not an evidence row in the JSON. |
| K5 WT, no-uptake | supported | WT 0.311660821633612 vs 0.31166082163360886; no-uptake ATP and biomass 0 on both, also with all 185 boundary reactions closed. |
| K6 acceptance | supported | Reproduced C1: baseline INO1-null growth 0.3116608 / inositol 511.43 / Ins1P 1000; candidate 0 / 0 / 0, also with R1720 open; phosphate-closed Ins1P 0. |
| K7 sweep | supported | No class change; only YALI1F31381g above 1e-6 (reproduced on both models); 15 %: 112 / 203 / 42 / 692 on both; baseline equals the dolichol `a801aad0` column exactly. |
| K8 census | supported (label inferred) | Counts and freed reactions match; Gurobi max 1.4708e-4 for R1343 / R1344 / R970 on both. The "GLPK fringe" mechanism is not shown, since 1.47e-4 is far above 1e-7. |
| K9 pytest | supported | Only the 4 new tests differ; 19 failures and 14 errors the same; focused rerun 4 passed, 30 subtests. |
| K10 guards and flag | supported | All guards before any edit; BEFORE / AFTER pinned; notes check present; idempotence and roundtrip tested; flag requirement in main.py; both review fixes tested. |
| K11 figure points | supported | Reproduced: baseline R2171 2.36e-7 with INO1 0; candidate INO1 2.34e-5 with R2171 0, growth 0.3116574; INO1 knockout on the candidate 0. |
| K12 budget | partially supported | Builds 5 / 6, acceptance 84 / 150, built checks 18 / 30, fresh 0 / 50, census 2. The baseline sweep re-run had no recorded reason; 21 LPs (census FVA 14, figure 7) are in no planned category (39 / 30 if counted as built checks). |

**Gaps the auditor raised and what was done:**
1. TASK.md said the dolichol sweep has a "`9d655c2b` side"; it has `ce9e1112` and `a801aad0`. Corrected in TASK.md's
   status paragraph; the re-run reason is now in REPORT section 6.
2. The candidate's `.build.json` records the round-1 step code; REPORT now cites `final_seed0/1.build.json` for the
   final code identity (same bytes).
3. E6 / E7 are not evidence rows in the JSON: stated in REPORT section 7. The JSON was not changed after the build.
4. `scripts/figure_fluxes.py` said 6 LPs, ran 7: docstring corrected. `scripts/ko_compare_two.py` had a stale
   docstring: corrected (code unchanged).
5. `logs/pytest_compare.txt` counts from summary lines differ from pytest's totals: a note was added.
6. The candidates README row linked a REPORT that did not exist yet, and there was no STATE entry: both now exist.
7. Benchmark positives: REPORT now states the 326 in-model positives and the 11 infeasible ones scored as unresolved.
8. The census mechanism: REPORT now says the GLPK call was not explained, only that Gurobi finds the three reactions
   unblocked in both models.
9. The wider R2171 FVA range: REPORT labels its explanation as inferred and untraced.

**Could not be checked by the auditor:** mass balance (DAG has no formula); a full re-run of the sweep, the census, the
full pytest suite and C2–C4 of the acceptance check (read from the deliverables); the Escher renders; the gene's real
identity or product.
