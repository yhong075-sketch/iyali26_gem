# OR-rule review: R171 ATP synthase single-gene branches and the cytosolic aminoacyl-tRNA synthetase rules (2026-10-29)

## 0. Authorization and scope

**User, 2026-10-29** (reply button after the assessment of levers for essential-gene recall):

> "Scope the OR-rule review first: R171 single-gene branches and the cytosolic aaRS rules.
> Investigation and TASK.md only, no model change yet."

Scope of this task:
- the ATP synthase rule on `R171`: `(14-subunit AND) or YALI1B28025g or YALI1B15733g`;
- the cytosolic aminoacyl-tRNA synthetase rules `R759` (threonyl, cytosol; `R760` in the mitochondrion
  carries the same rule), `R723` (seryl, cytosol) and `R673` (prolyl, cytosol), each `cytosolic-type gene or
  mitochondrial-type gene`.

For each rule: record the original structure, the mechanism by which the OR hides the knockouts, the evidence
for and against each member, the candidate edits, and their acceptance criteria (AGENTS §5). Score the
candidate edits in memory at the four benchmark cutoffs. Nothing is changed: no model XML, GPR, curated data,
medium, benchmark or platform code; no commit, no push.

Out of scope (recorded for later rounds): the 16 other OR partners of the 2026-10-27 study (R712, R60, R2193,
the R1556 group, R540, R4, R289, R211, R121, R815, R251/R177, R1843/R1846, R78); the mitochondrial
synthetase reactions with dead-end tRNAs (category d); biomass demands; the amino-acid uptake caps.

## 1. Starting state

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/essential-gene-prediction-e76be3` @ `ff5192d` (merge of PR #32), clean tree |
| Model executed | [`candidates/E5_…_dolichol_route_v2_r2171_1to1_20261028.xml`](../../candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_r2171_1to1_20261028.xml), SHA-256 `965ab947f2aadc2c69443e7a794358c502acfdf0892e78582278ff07a37429d7` (PR #31, commit `af2ebab`) |
| Sibling candidate | `5f54eb22086b06b2bd5e8715b89984060cd43c1a7fa21e33e5581d867ee54179` (UM7, PR #30, commit `3720370`): identical gene set and identical `R171`, `R759`, `R760`, `R723`, `R673` (stoichiometry, bounds, rule), checked this session |
| Reference model | none designated (`candidates/README.md`) |
| Medium / profile | Task A medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) + PO1f profile `po1f_sd_leu.json` (`35307853…`), benchmark default |
| Experimental list | `consensus_essential_genes.csv` (`1e887f5a…`), 1,612 positives, 326 in the model |
| mtDNA scope rule | `expected/mtdna_genes_out_of_scope.json` (`16b3e9b4…`) |
| Baseline scores (recalled from the R2171 and UM7 reports, not yet re-run here) | 15 %: TP / FN / FP / TN 112 / 203 / 42 / 692, 16 unresolved; user-defined negatives; development reference |
| Research workspace | `$IYALI26_RESEARCH_ROOT` set; `reference/`, `state/`, `artifacts/` present; input SHAs verified |
| Solver | Gurobi 13 (academic licence, real `$HOME`), FeasibilityTol / OptimalityTol 1e-9, Presolve 0, 1 thread, basis reset per knockout, knocked-flux and residual checks as in `reports/essentiality_fn_causes_20261027/run_screen.py` |

## 2. Computation budget

- At most 1,500 LPs and 30 minutes wall time: WT; the 22 single-gene knockouts of the genes in the five rules;
  reaction-level knockouts; the same 22 knockouts under each in-memory rule variant; one full 1,081-gene
  sweep of the recommended combined variant as a confirmation.
- A rule edit on a reaction can change only the knockouts of genes that appear in that reaction's old or new
  rule (a gene knockout changes bounds only on reactions whose rule contains the gene). The full sweep is a
  check of that statement, not the source of the numbers.
- Spend (added after the runs): three attempts of `scripts/variant_knockouts.py`. Run 1 stopped at import (a
  wrong repository-root path; 0 LPs). Run 2 completed the 205 variant LPs and 139 sweep LPs, then stopped at the
  restore assertion for a gene cobra keeps detached from the model (344 LPs). Run 3 is the complete record. No run
  exceeded the per-run budget; the total over the three attempts is given in REPORT.md §2.
- Stop conditions: any solve that fails the stale-point or residual check is reported, not re-tried with
  other settings; no bound, medium or coefficient is changed to make a variant work.

## 3. Evidence to gather (read-only)

- Identity of every gene in the five rules: UniProt (CLIB122 `YALI0` entries through the S2 locus map),
  the genome annotation strings, Yeast-GEM 9's rules for the same enzymes (convention only, never copied).
- Per-screen calls (Cas9, Cas12a) for every gene, from `cas9_cas12a_fitness.csv`.
- *S. cerevisiae* deletion literature for the ATP synthase dimer-specific subunits (Europe PMC).
- Localisation of the synthetase paralogs: UniProt annotation, yeast ortholog, N-terminal sequence. A formal
  presequence predictor (TargetP 2.0 / DeepLoc 2.1 web services) needs a form submission and is run only if
  the user approves.

## 4. Outputs

- This folder: `TASK.md`, `REPORT.md`, `or_rule_review_genes.csv` (one row per gene: identity, screen calls,
  baseline and variant knockout growth, class at each cutoff), `variants.json`, `scripts/` (the exact
  scripts run), `evidence/` (fetched UniProt TSV with SHA-256).
- A dated `model/STATE.md` entry. No large outputs are expected; if any arise they go to
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/or_rule_review_r171_aars_20261029/`.
