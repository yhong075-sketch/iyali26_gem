# Independent audit of Stage 2 (R2121 closure), 2026-11-12

**Auditor.** A fresh read-only subagent (general-purpose). It received TASK.md, the uncommitted diff, the logs, the
deliverables and a list of 10 claims, not the main session's reasoning. It wrote scratch files only under
`/private/tmp/claude-501/r2121_audit/`. It ran 12 of the 50 LPs it was allowed. It ran its own pytest of the focused
file. This file was compiled by the main session from the auditor's report; the verdicts are unchanged.

**Summary.** 8 claims supported, 2 partially supported, 0 unsupported.

| Claim | Verdict | Auditor's evidence (summary) |
| --- | --- | --- |
| C1 authorization, scope, criteria before results | Supported (authorization checked only as recorded) | TASK.md §8 has the quote, the scope, 7 acceptance criteria, the budget and the stop conditions. TASK.md was saved at 10:00 local and the first build started at 10:11. The auditor can see the recorded quote, not the user's message. |
| C2 opt-in, last step, fail closed, idempotent | Supported | `main.py` 959–965: after the crosswalk, just before the save; flag default False. Every check runs before the edit and raises ValueError. Own focused run: 6 passed, 43 subtests. |
| C3 control reproduces e584151f; seed 0 = seed 1 | Supported | Own SHA-256: control `e584151f…39d2`; candidate = seed 1 = `2f352037…af7c`. 3 builds, all exit 0; spare unused. Control argv equals e584151f's apart from the output path. All 73 source and 146 data SHAs in the `.build.json` match the current files. |
| C4 whole-model diff = 2 entries | Supported | Raw XML diff: 2 hunks, both R2121 (reaction tag bounds; 7 new note lines after the 4 unchanged ones). |
| C5 chemistry unchanged; balance | Supported | R2121: mass balanced, charge +1 (products minus substrates) in both files. R2119: balanced and unchanged. |
| C6 WT and energy | Supported | Checked solve on the candidate: WT 0.3120526953810365, optimal, passed. No-uptake ATP maximum 0, optimal; only `xMAINTENANCE` relaxed. Baseline WT 0.31248856100769146 reproduced. |
| C7 sweep | Supported | `ko_screen_built.tsv` is byte-identical to `ko_screen_v1.tsv`. Rescore from raw KO growth and the consensus CSV (`1e887f5a…`) gives 141/185/43/694 → 143/183/43/694 at 15 % (also right at 1, 5, 10 %). Only YALI1E35955g (1.0 → 0) and YALI1F28274g (0.159 → 0) change; no TP lost. Both knockouts re-solved on both models. |
| C8 evidence in the curation file | Supported, with wording points | KEGG REST today: R00740 → EC 1.2.1.75 (3-HP and 3-HP/4-HB cycles of *some* green non-sulfur bacteria and thermoacidophilic archaea); yli:2909116 → K00140; R00705 decarboxylating; `link/yli` and `link/sce` for K15017 empty. Yeast-GEM `r_4573` = R00705, no R00740; eciYali5 neither. |
| C9 Escher numbers | Partially supported | Files exist and SHAs match. Own pFBA reproduces all 51 fluxes within 1e-9 and the totals 294.02 / 282.16. "46 before" should be 45. |
| C10 nothing outside scope | Supported | Git changes only `cli.py`, `main.py`, `STATE.md` (the Stage 1 entry) plus the listed untracked files. e584151f and `model/source` unchanged. |

## Findings and responses

1. **C9: 46 should be 45.**
   - Cause: the map build counted flux at |v| ≥ 1e-12, but status used 1e-9. R2119's −1.74e-12 in the e584151f pFBA is
     solver noise.
   - **Fixed:** `escher/scripts/build.py` now uses one threshold, and the map, viewer, `index.json` and README were
     rebuilt (no LP). `index.json` now gives 45 / 11. Both viewer checks were re-run: 51 reactions in every scenario.
2. **C8: gene identity link.** KEGG's yli:2909116 entry is the CLIB122 gene (YALI2_C01175g) and does not name
   YALI1C02550g. The link rests on the model's gene annotation (ncbigene 2909116) and the W29 audit row (YALI0C01859g,
   XP_501343.1, `ok`), which matches KEGG's protein ID. **Recorded** in REPORT.md.
3. **C8: second KO.** R00740 also has K14468. The auditor checked it: no yli or sce genes. This strengthens the
   conclusion. **Recorded** in REPORT.md.
4. **C8: "only".** TASK.md §3 said Yeast-GEM carries "only R00705"; it also has R00706 (`r_4574`, the NADP form of
   the same decarboxylating reaction). **Fixed** in TASK.md with a dated note. The curation JSON and the model note were
   already correct ("R00705 and not R00740").
5. **C8: how R2119 matches R00705.** R2119 has no `kegg.reaction` annotation, so the match is by stoichiometry. H+ sits
   on the opposite side by the model's protonation convention. **Recorded** in REPORT.md.
6. **C4: serialized form.** The SBML `reversible` attribute flips from true to false. cobra's writer derives it from
   the bounds, so it is not a separate change. **Recorded** in REPORT.md.
7. **C2: minor.**
   - The RuntimeError after the edit is unreachable in practice.
   - The step's position in `main.py` is not tested; it is confirmed by code reading and the build log line
     "R2121 malonyl-CoA closure: applied".
   - **No change.**
8. **Unrelated: R304.** The libSBML notes error on R304 appears in every build, the control included. Inherited;
   affects another line. **No change.**
