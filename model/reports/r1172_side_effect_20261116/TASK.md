# Why closing R1172 blocked nine other reactions in 7ccbb568 (investigation, 2026-11-16)

Computation and documentation only. No model XML, curated data, builder code, benchmark or test is changed; no
build; nothing is committed or pushed.

## 0. Authorization

- The task came from the coordinating session `local_1088c84c…` (gap census line): investigate, without changing
  the model, why closing R1172 blocked R25, R552, R95, R2074, R2075, R2076, R2077, R2112 and R834 in
  `model/candidates/E5_353cefdc_xref_20261115.xml`, and lay out the options for the user to decide. It names the
  order of the work (reproduce, mechanism, gene rules, options with evidence, in-memory tests, figure) and says:
  "Do not propose reopening R1172 unless you find evidence for HMG-CoA transport; that closure is a recorded user
  decision. Never activate a reaction by adding an unsupported exchange, sink or transport."
- The user, 2026-11-16, in this session: **"Approve plan"**, answering "Approve this computation plan (at most
  1,200 LPs, no build, no sweep)?" after the plan in section 3 was shown.

Not authorized: any model or curation change, builds, a full knockout sweep (ask first), commits, pushes, branch
switches, edits in other worktrees, cluster jobs.

## 1. Starting state (checked 2026-11-16, 15:28 PST)

| Item | Value |
| --- | --- |
| Worktree / branch / HEAD | `bold-bartik-794e7e`, `claude/bold-bartik-794e7e`, `af7c543214459f22d3947d71c8fc1a75ea5fb800` (main after PR #54), clean |
| Model executed | `model/candidates/E5_353cefdc_xref_20261115.xml`, SHA-256 `7ccbb56834a37ad4c75e434868e9e1b9e1f6a657a11ca17ee0a93dbe56232146` (committed in `038473c4`) |
| Reference model | none designated (`model/candidates/README.md`) |
| Medium | `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv`, `6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45` |
| Strain profile | `$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json`, `35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383` (in memory) |
| Experimental list | `$IYALI26_RESEARCH_ROOT/state/essentiality/repository/consensus_essential_genes.csv`, `1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b` |
| Solve helper reused | `model/reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py`, `69c42690…` |
| Software | `/opt/miniconda3/bin/python3` (Python 3.13.5), cobra 0.30.0, optlang 1.8.3, Gurobi 13.0.2 (academic licence, real HOME), GLPK via swiglpk |

R1172 on 7ccbb568: `m646[C_cy] <=> m648[C_mi]` (HMG-CoA), no gene, bounds (0, 0) since candidate 684f3752
(user decision 2026-11-12, `model/STATE.md`). Its bounds before were (−1000, 1000).

Evidence the task starts from (read, not re-run): the gap census of 7ccbb568,
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/gap_census_candidates_20261116/` (`out/7ccbb568_open_closure_trace_out.json`,
54 LPs, GLPK: reopening R1172 frees the nine, reopening R2121 frees none), and its uncommitted report in the
worktree `metabolic-pathway-gnn-priorities-9b8f18` (section 3.3).

Related work on another line (read only, not changed): the methylcitrate round in the worktree
`gene-essentiality-reactions-d09601` (`model/reports/methylcitrate_cycle_20261116/`, uncommitted). It adds a matrix
2-oxobutanoate dehydrogenase with an inferred PDH rule and closes R2075, on base cbf2ef95 (candidate b24a5e78).

## 2. Question

Which stoichiometric fact makes the nine reactions depend on R1172, which statements the LP proves, and which
biological options could free them without reopening R1172 and without an unsupported exchange, sink or transport.

## 3. Plan (approved)

Conditions: (a) every exchange open at (−1000, 1000) on the raw model (the census definition); (b) Task A medium +
PO1f profile through `load_effective_simulation_context`. Solvers: Gurobi 13 at feasibility and optimality 1e-9
(primary) and GLPK with the basis reset per LP (cross-check). Raw status kept for every LP; a failed solve is not
evidence.

1. **Reproduce.** FVA (no growth floor) of the nine plus the CoA-moiety boundary reactions (R834, R1172, R2025,
   R2077, R2083, R2084, R2113, R2118) as is and with R1172 at (−1000, 1000), both conditions, both solvers.
2. **Mechanism.** Static CoA-moiety count per species and a check that every reaction inside the matrix conserves
   the matrix moiety; producers, consumers and transporters per compartment. LP: minimum and maximum net moiety
   import into the matrix with R1172 open and closed; the minimal-flux route that carries each of the nine with
   R1172 open, with the CoA drain and with the FAD drain; each route classed as net-zero cycle or net transfer.
3. **Genes.** Rules as recorded, consensus-list membership, W29 tag audit, recorded calls. No LP.
4. **Options**, with a clickable source and a tier for each external claim (at most 30 read-only fetches):
   O0 leave blocked; O1 matrix propionyl-CoA producer (the 2-oxobutanoate dehydrogenase of the other line);
   O2 valine/methylmalonate steps in the matrix; O3 propionyl-carnitine shuttle (only with yeast evidence);
   O4 matrix CoA turnover (carrier exchange mechanism or growth dilution). Excluded: reopening R1172; R834
   reversible without evidence; any unsupported exchange, sink or transport.
5. **In-memory tests per option:** WT; FVA of the nine; energy check (uptakes closed, `xMAINTENANCE` lower bound 0,
   maximise ATP hydrolysis and biomass, status checked); targeted knockouts of the genes on the routes. A full
   sweep only after asking.
6. **Figure:** hand-laid schematic (cell, mitochondrion, one box per reaction with gene rule and flux) for before
   (R1172 reopened in memory), after (7ccbb568) and the leading option.

## 4. Budget and stop conditions

Budget: at most 1,200 LPs over both solvers; no build; no full sweep; at most 30 read-only web fetches; one
optional independent read-only review of the report.

Stop and report if: step 1 does not reproduce (a reaction of the nine is not blocked as is, or not freed with
R1172 reopened, or the solvers disagree on a blocked set); a decisive solve stays non-optimal after one fresh
re-solve; an option's energy check is above 0 (flag it, do not recommend it); the budget would be exceeded; the
evidence would require revisiting the R1172 decision.

## 5. Outputs

This folder: TASK.md, REPORT.md, small tables, `scripts/`, `figure/`. Raw outputs and logs:
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/r1172_side_effect_20261116/`. A dated entry in `model/STATE.md`.
