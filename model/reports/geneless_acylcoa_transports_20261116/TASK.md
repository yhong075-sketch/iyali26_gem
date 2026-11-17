# Gene-less CoA-ester transports across the inner mitochondrial membrane in 7ccbb568 (investigation, 2026-11-16)

Computation and documentation only. No model XML, curated data, builder code, benchmark or test is changed; no
build; nothing is pushed. Committed on the branch at the user's later instruction (section 7).

## 0. Authorization

- **Task source.** The coordinating session `local_ebd14d46…`, following the R1172 side-effect round
  (`bold-bartik-794e7e`, `model/reports/r1172_side_effect_20261116/`, uncommitted). It asks, without changing the
  model, for each gene-less transport of a CoA species across the inner mitochondrial membrane in
  `model/candidates/E5_353cefdc_xref_20261115.xml`:
  1. its stoichiometry, the direction it can carry today, and whether it is needed at the WT optimum or for any
     gene's knockout result;
  2. evidence for a carrier;
  3. for R2117, whether the valine steps R2114–R2117 belong in the matrix;
  4. options per reaction, tested in memory;
  5. a hand-laid schematic figure.

  The transports: R2025, R2113, R2077, R2083, R2084, R2118, R1291, plus R2117.
- **The user, 2026-11-16, in this session:** **"Approve, add TargetP"**. This answered "Approve this computation plan
  (≤1,300 LPs, ≤40 fetches, no build, no sweep)?" after the plan in section 3 was shown. The "add TargetP" option
  meant: send the W29 RefSeq proteins of the valine/leucine steps to the official TargetP-2.0 web server, and label
  the output as a prediction.

Not authorized:
- any model, curation, benchmark or test change, and any build;
- a full knockout sweep (ask first);
- commits, pushes, branch switches, and edits in other worktrees;
- cluster jobs;
- MitoFates or DeepLoc runs (not asked).

## 1. Starting state (checked 2026-11-16, 16:33 PST)

| Item | Value |
| --- | --- |
| Worktree / branch / HEAD | `lucid-visvesvaraya-a40bbb`, `claude/lucid-visvesvaraya-a40bbb`, `5c497039ad00cd004939f6f889390627c4d6bea4` (main after PR #55), clean |
| Model executed | `model/candidates/E5_353cefdc_xref_20261115.xml`, SHA-256 `7ccbb56834a37ad4c75e434868e9e1b9e1f6a657a11ca17ee0a93dbe56232146` |
| Reference model | none designated (`model/candidates/README.md`) |
| Medium | `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv`, `6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45` |
| Strain profile | `$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json`, `35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383` (in memory) |
| Experimental list | `$IYALI26_RESEARCH_ROOT/state/essentiality/repository/consensus_essential_genes.csv`, `1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b` (1,612 genes) |
| W29 proteome (sequences for TargetP) | `$IYALI26_RESEARCH_ROOT/artifacts/tasks/orphan_gene_search_20261019/db/w29.faa`, `1863842c85b00bc6e7c37aecf9b6231370cb33c6168149d158a4fa1ccf9f5fcf` |
| Software | `/opt/miniconda3/bin/python3` (3.13.5), cobra 0.30.0, optlang 1.8.3, Gurobi 13.0.2 (academic licence, real HOME), GLPK via swiglpk |

What this round starts from:
- **Read, not re-run:**
  - the R1172 side-effect report (sections 3, 5, 8) and its STATE entry, both uncommitted in `bold-bartik-794e7e`;
  - the gap census of 7ccbb568 (`model/reports/gap_census_candidates_20261116/`, committed in PR #55).
- **Read only, second-hand:** the FN-activation scoping notes F07 and F11 and the methylcitrate round, both
  uncommitted in `gene-essentiality-reactions-d09601`.

## 2. Question

For each of the eight reactions:
- what it does as written;
- whether the model needs it;
- whether yeast has a carrier for its species;
- which change (keep, close, move to another compartment, merge the duplicate cytosolic propionyl-CoA species
  m1875/m1970) would match the biology without an unsupported exchange, sink or transport, and what each change does
  to growth, energy and gene calls.

## 3. Plan (approved)

Conditions:
- (a) every exchange open at (−1000, 1000) on the raw model;
- (b) Task A medium plus the PO1f profile, through `load_effective_simulation_context`.

Solvers:
- Gurobi 13 at feasibility and optimality 1e-9, with solution and basis discarded per LP (primary);
- GLPK with the basis reset per LP (cross-check of decisive values).

Raw status is kept for every LP; a failed solve is not evidence.

1. **Capacity and need.**
   - FVA with no growth floor of the 8 plus about 22 neighbours, on both conditions; GLPK checks the 8.
   - At the WT optimum: FVA of the 8 at 100 % and 99.9 % of growth, and one pFBA.
   - WT with each transport closed alone, and with R2025 and R2113 both closed.
   - About 30 targeted gene knockouts, as is and with the transports closed, classed at 1, 5, 10 and 15 %.
2. **Carrier evidence.** Read-only Europe PMC, UniProt, KEGG and SGD searches, at most 40 fetches. Every claim gets a
   clickable source and a tier. The carnitine system moves acyl groups, not the CoA moiety.
3. **R2114–R2117 compartment.**
   - Gene identities and homologue locations.
   - TargetP-2.0 (Non-Plant) on the W29 RefSeq proteins, after the W29 tag audit, labelled as a prediction.
4. **Options in memory**, on the as-is model and on the matrix 2-oxobutanoate dehydrogenase package of the other line
   (R2075 closed):
   - keep;
   - merge m1875 into m1970;
   - merge plus the propionyl-CoA transports closed;
   - the blocked transports closed;
   - the valine steps moved into the matrix with R2118 closed;
   - all combined.

   Each variant gets WT, FVA, the energy check (uptakes closed, `xMAINTENANCE` lower bound 0, maximise ATP
   hydrolysis and biomass) and targeted knockouts. A full sweep is only run after asking.
5. **Figure.** A hand-laid cell and mitochondrion schematic, with enzyme boxes carrying the gene rule and flux.

## 4. Budget and stop conditions

**Budget:**
- at most 1,300 LPs over both solvers;
- at most 40 read-only literature and database fetches;
- the TargetP submission;
- no build, no full sweep, no cluster job;
- one optional independent read-only review of the report.

**Stop and report if:**
- an input SHA-256 does not match;
- a decisive solve is still non-optimal after one fresh re-solve;
- a variant's energy check is above 1e-9 (flag it, do not recommend it);
- the budget would be exceeded;
- the evidence points at the R1172 decision or at the colleague-owned lipid-unlump files. R2077 is listed in
  `model/curation/lipid_unlump/fatty_acyl_coa_handoff.json` (`repaired_reaction_ids`), so for R2077 this round only
  reviews and recommends.

## 5. Outputs

- **This folder:** TASK.md, REPORT.md, small tables in `tables/`, `scripts/`, `figure/`.
- **Raw outputs, logs, literature and TargetP files:**
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/geneless_acylcoa_transports_20261116/`.
- **A dated entry in `model/STATE.md`.**

## 6. Follow-up: full certified knockout sweep of the leading option (authorized 2026-11-16)

- **The user, 2026-11-16, in this session:** **"Run the sweep"**. This answered "Run a full certified knockout sweep
  (about 1,100 LPs) on the leading option?", asked after REPORT.md was written.
- **Leading option:** variant `B0/ALL` of step 4, applied in memory to 7ccbb568. It merges m1875 into m1970 (R2025 then
  a duplicate, set to 0), closes R2113, R2077, R2083, R2084, R1291 and R2118, and moves R2117, R2116, R2115, R2114 and
  R2112 into the matrix with their gene rules. It is the same code as step 4 (`apply_base`, `apply_variant`), so the
  edits are identical to the ones tested there. Nothing is written to SBML.
- **Harness:** `scripts/step6_sweep.py`, a copy of the benchmark's reference implementation
  `model/reports/infeasible_ko_rule_20261030/scripts/ko_screen_certified.py` (SHA-256 `313d5d01…`), with
  `checked_solve.py` (`69c42690…`) unchanged. The only change is the loader: each context (sweep, fresh re-solve,
  certificate) reads 7ccbb568 through `load_effective_simulation_context` and then applies the variant. Everything else
  is the same: Gurobi at 1e-9 with the model reset per knockout, the targeted check, `FRESH_LIMIT` 20, the
  maintenance certificate, Rule A and Rule B, and the mtDNA and untargeted scopes.
- **Baseline (as is):** the certified sweep `sweep684_candidate` of the folinic-acid round, run with the same two
  harness files, medium, profile, list and untargeted file on 353cefdc. 7ccbb568 is the xref-cleanup build of 353cefdc
  (annotation and notes changes only, per that round's whole-model diff). Before the baseline is used,
  `scripts/step6_identity.py` checks statically, with no LP, that the two models load as the same optimisation problem:
  - the same reactions, stoichiometry, effective bounds, objective and gene rules;
  - the same genes and runtime exclusions.

  If they differ, the baseline is not used and the as-is sweep is asked for separately.
- **Control:** every gene whose class changes between the baseline and the option is re-run on the as-is model in the
  same harness (variant `B0/K`), to show the change comes from the option and not from the harness.
- **Budget:** at most 1,300 LPs in all, on Gurobi only. The reference sweep needed 1,123 (1,082 + 1 + 40). No build, no
  model or curation change, no cluster job, no commit.
- **Stop and report if:**
  - an input SHA-256 does not match;
  - the WT point fails the targeted check;
  - more than 20 knockouts need a fresh re-solve;
  - the identity check fails;
  - the budget would be exceeded.
- **Outputs:**
  - the per-gene TSV and the logs go to the artifacts folder;
  - the summary JSON and the class-change table go to `tables/`;
  - scores are reported at 1, 5, 10 and 15 % under Rule B, with Rule A alongside. Negatives are the user-defined class,
    and the list is a development reference.

## 7. Commit (2026-11-16)

- **The user, in this session:** **"Commit it"**, answering "Commit this round's report, scripts and STATE entry on
  this branch?".
- Committed on `claude/lucid-visvesvaraya-a40bbb`: this folder and the `model/STATE.md` entry. Not pushed; no pull
  request.
