# Independent audit: lipoylation module + acyl-ACP chemistry on 2c9cebd0 (2026-10-23)

- **Auditor:** a fresh read-only subagent with no conversation history. It was given only TASK.md,
  commands.txt, the deliverables in this folder, the working-tree code changes and commit `9f2356b`,
  plus the eight claims below.
- **What it ran:** no build, no pytest, no sweep, and no repository edits. It did run 4 LPs: WT on both
  models, and a knockout of YALI1F34729g on the candidate with and without the lipoyl biomass terms.
- **This file:** its report, reformatted. Responses are in [REPORT.md §7](REPORT.md#7-independent-audit).

## Verdict per claim

| Claim | Verdict | Evidence the auditor checked | Caveats it raised |
| --- | --- | --- | --- |
| **C1.** Copied files are byte-identical to `9f2356b`; only the stated `cli.py`/`main.py` edits | Supported | `git show 9f2356b:<f> \| shasum` vs the working files: identical for both modules, both JSONs, `test_lipoylation_module.py` and all 25 + 21 report-folder files. `test_acyl_acp_chemistry.py` differs only in `BASELINE` and `setUpClass`. `git diff HEAD`: exactly the 2 flags; the 2 parameters, the rule, the 2 calls after R2193 and the 2 kwargs. | `main.py` has one extra comment line explaining the order. |
| **C2.** Running after heme/iron + R2193 is necessary and sufficient | Partially supported | The heme guard raises `biomass_C composition differs from the reviewed baseline` on the candidate and returns `already_correct` on `2c9cebd0`. The lipoylation guard compares only its 4 terms. Only CoQ steps that are not in the argv follow the new calls. In the XML diff, every changed element except `biomass_C` is byte-identical to `b3b9fe22`. | Sufficient, yes. Necessary only "after heme/iron": the R2193 guard reads only R2193, so a slot between heme and R2193 would also work (not tested). `tpp_requirement.py` also pins the full `biomass_C`. That is satisfied because TPP runs before heme, but the comment and TASK name only the heme guard. |
| **C3.** Test changes as described; the suite is identical test by test plus 13 passes | Supported | Outcome lines: 316 before, 329 after; none removed; 13 added, all PASSED. A rebuild of `164f235f` + lipoylation + deterministic write gives `fdd00ac3db7c…`, byte-identical to `9f2356b`'s `fdd00ac3`. | `pytest_after` overlapped the builds and the sweep, but commands.txt reads as sequential. No impact on the outcomes. |
| **C4.** The control reproduces `2c9cebd0`; seed 0 = seed 1 | Supported | SHA-256: candidate and seed 1 `a12378210c7c…`; control and baseline `2c9cebd0a5c5…`. All four `.coq9_genes.tsv` are `2dd602b7…`. In build.json: the argv is the base flags + 2, `source_unchanged_during_build` is true, and every source SHA matches the current files. | Seeds are recorded only in commands.txt. The control and seed-1 XMLs were only in the session scratchpad. |
| **C5.** Diff is exactly the targets, with the same entities as `b3b9fe22` vs `164f235f` | Supported (re-derived) | Own element-level XML diff, both pairs: 5 reactions, 6 species and 3 genes added; 22 reactions and 34 species changed. Each pair has one extra "changed" element, which is whitespace only. cobra field-level diff: 71 entities in each pair, the same set and the same fields. The only after-value difference is `biomass_C`; its deltas are identical. | `biomass_C` after-values differ by the base's pre-existing **stoichiometry** (TPP and heme/iron terms) as well as its notes. |
| **C6.** Touched reactions are balanced or at the curated residual; R2290 is a sink by design | Supported | 27 touched reactions vs the specs: no mismatch. 21 are `{}`; 6 have charge −2 or −3, equal to `residuals_after`. R2290 is the 5'-dA drain, documented in `lipoylation_module_20261020/`. | 12 untouched reactions also change balance through the formula fills (R78, R539, R1400–R1402 and R2216 → `{}`; R1397–R1399 and R1403–R1405 → charge −1). All equal the curated records. |
| **C7.** WT values; the sweep flags nothing | Supported | Own WT LPs: 0.312426291487999 and 0.31242085277746734, both optimal and passed. `reset()` runs before every solve. Summary: 2,161 LPs, no failing points, no fresh re-solves, max residual 2.97e-12. | 20 infeasible knockouts per model (4 mtDNA, 16 in scope) are unresolved by design, not "flagged". The sweep was not re-run. |
| **C8.** 13 class changes; 15 % scores | Supported (recomputed) | Own read of the experimental list (`1e887f5a…`, 1,612 IDs): the same 13 genes at all 4 cutoffs, all at KO growth 0.0. They are the 3 new genes + 10 of the 11 in the `b3b9fe22` screen; YALI1E07744g is essential in both. 15 %: 93/220/36/697/16 (MCC 0.3454) → 103/212/39/695/16 (0.3669). | 324 → 326 in-model positives. These are a development reference, not validation. |

## Other findings

1. **The lipoyl requirement drives the new essentiality** (verified for one gene).
   - Candidate with YALI1F34729g knocked out: growth 0.0. The same knockout with the 4 lipoyl
     `biomass_C` terms removed in memory: 0.31242629148799905, exactly the `2c9cebd0` WT.
   - So that gene's call, and the whole WT cost, come from the provisional α = 1e-5. The acyl-ACP block
     has no effect on WT here.
   - The other 12 genes' mechanisms are inferred. The recall gain is not support for the change.
2. **Reruns of earlier steps.**
   - The heme/iron step raises on `a1237821` (verified; new).
   - The TPP step already raises on `2c9cebd0` itself (`m975[C_cy] identity differs`, verified), so
     that part predates this change.
   - Thiamine-P, proton lock, R2193, lipoylation and acyl-ACP all return `already_correct` on the
     candidate.
3. **Test quality.**
   - The target test iterated only `after` keys, so a removal would go unnoticed.
   - The heme test used a bare `assertRaises`.
   - No test checks the builder order.
   - Pre-existing: the lipoylation `_verify` does not check the `biomass_C` note on a rerun.
4. **Gene labels** in the README row and STATE entry come from model reaction names, not verified names.
   For example, YALI1C07638g also carries R1869 (C18:1 → C18:2), so "Δ9" is partial.
5. **Evidence retention:** the control, seed-1 and `b3b9fe22` XMLs were only in the temporary
   scratchpad.
6. **Cosmetic, inherited from `9f2356b`:**
   - The new genes are named `G_YALI…`, where others are `COBRAProtein…`.
   - cobra shows `&apos;` in notes; the file is escaped correctly and the guards unescape.
7. `STATE.md` and `README.md` changed during the audit (the requesting session's entries, including a
   placeholder). They were checked for consistency only.
