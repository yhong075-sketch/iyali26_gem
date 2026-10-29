# STOP: one knockout gains growth (UM7 nuclear diffusion, 2026-10-28)

The build, diff, balance, energy-cycle, blocked-set and class checks all passed. One approved acceptance criterion
failed (TASK.md section 8: "no knockout gaining absolute growth above 1e-6"), so the run stops here, before the
audit, REPORT.md and the STATE.md entry. Nothing was tuned, relaxed or retried. The candidate file stays in
`model/candidates/`, and nothing is committed.

**Also completed before the stop:**
- **pytest:** 314 → 320 passed, with the same 19 failed and 14 errors before and after. The only difference is
  the 6 new tests ([logs/pytest_compare.txt](logs/pytest_compare.txt)).
- **Caveat on the before run:** it started before the main.py / cli.py edits but ran while they were made.
  Tests that re-read those files may have seen the opt-in flag.

## What failed

| | 9d655c2b | candidate 5f54eb22 |
| --- | ---: | ---: |
| WT (Task A + PO1f, Gurobi 1e-9, checked) | 0.3116608216 | 0.3116608216 (+4.6e-15 relative) |
| `YALI1F37714g` knockout | 0.3114518395 (KO/WT 0.999329) | 0.3115724887 (KO/WT 0.999717) |

- **Gain:** +1.21e-4 absolute in one knockout. Every other knockout changes by ≤ 1e-9.
- **Class:** no change at any cutoff (1, 5, 10, 15 %).
- **Scores:** TP/FN/FP/TN identical, 112 / 203 / 42 / 692 at 15 % (16 unresolved, the same infeasible
  knockouts in both). Negatives are the user-defined class; development reference only.
- **Sources:** `ko_compare_summary.json`, `ko_compare.tsv`; 2,164 checked LPs, no point failed the stale-point
  check.

`YALI1F37714g`: KEGG EC 2.4.2.1, purine-nucleoside phosphorylase. Its name is not verified. Its GPR covers 13
nucleoside-phosphorylase reactions (R264, R675–R684, R761, R1903). It is not on the experimental list.

## Mechanism (verified this session, `pnp_diag` below; pFBA plus four controls)

In the knockout, the candidate's pFBA routes **0.003646** through the new links:

1. Tryptophan uptake rises from 0.001674 to 0.005320, its medium bound.
2. The kynurenine route R516 → R563 → R499 → R41 → R43 → R1376 makes quinolinate.
3. R581 / R580 make NaMN and then deamido-NAD in the cytosol.
4. R1082 moves deamido-NAD into the nucleus, and R569 makes NAD there (Gln in by R1214, ATP in by R1021, Glu and AMP out by R1212 / R1017; corrected after the audit).
5. **R567 hydrolyses the NAD to ADP-ribose and nicotinamide.**
6. ADP-ribose leaves by **R_ADPRtn** and is split by R1912 into AMP and ribose-5-P. PPi leaves by **R_PPItn**.
   Nicotinamide leaves by R1267 and becomes nicotinate by R576.

Before the change, de novo NAD had only one sink, the biomass NAD term (0.001003).

**Controls** (candidate knockout, Gurobi 1e-9, checked):

| Closed reaction | Knockout growth |
| --- | --- |
| R_ADPRtn | 0.3114518395 |
| R_PPItn | 0.3114518395 |
| R567 | 0.3114518395 |
| R1912 | 0.3114518395 |

Each control returns exactly the baseline value. In the WT, R567 and R1912 carry 0 at the optimum.

**Reading (inferred):** the links turn the lumped R567 ("NAD + H2O → ADP-ribose + nicotinamide", no acetylated
substrate) into an unconstrained NAD hydrolase. Together with NAD synthesis, that makes an extra
tryptophan-catabolism and ribose-5-P / AMP route, which the PNP knockout can use a little. The route is made of
real steps, and sirtuin NAD turnover with Bna NAD synthesis is real yeast biology. What is not real is the
coupling: in yeast, sirtuin NAD consumption is tied 1:1 to protein deacetylation, which produces O-acetyl-ADP-ribose
(Tanner 2000, PMID 11106374; section 3 of TASK.md). The model has no deacetylation demand that would set or
limit that flux. This is a mechanism statement. Whether the PNP knockout's real growth changes is not tested.

## Options for the user

- **Accept and restate the gate.** Keep `5f54eb22`, record the +1.2e-4 on `YALI1F37714g` as a documented
  consequence (no class change), and finish the run: audit, REPORT, STATE entry, candidates README.
- **Constrain R567 first (new science change, needs its own scope).** For example, cap or remove R567's free
  hydrolysis, or rewrite it as deacetylation with an acetyl-protein source. This was not authorized, and the
  evidence for a value is missing.
- **Do not adopt.** Keep 9d655c2b as the newest candidate and leave `5f54eb22` unadopted. The file is kept until
  the user decides.

## pnp_diag (excerpt)

Fluxes in mmol/gDW/h; pFBA in the knockout, baseline → candidate.

| Reaction | Baseline | Candidate |
| --- | ---: | ---: |
| R567 | 0 | 0.003646 |
| R569 | 0 | 0.003646 |
| R1912 | 0 | 0.003646 |
| R_ADPRtn | – | −0.003646 |
| R_PPItn | – | −0.003646 |
| R576 | 0 | 0.003646 |
| R580 / R581 | 0.001003 | 0.004649 |
| R516 / R563 / R499 / R41 / R43 | 0.001003 | 0.004649 |
| Trp exchange R1232 | −0.001674 | −0.005320 |
