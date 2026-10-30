# Independent audit: R171 ATP synthase rule (Option D), candidate ca43ae0c (2026-10-29)

A fresh read-only general-purpose subagent with no conversation history audited the round. It was given TASK.md,
REPORT.md (draft), the deliverables, logs, code, data and the review, and a budget of 40 LPs. It edited nothing in the
repository; its scratch scripts (`semdiff.py`, `spot.py`, `raw_xml.diff`, `pytest_r171.txt`) were in the session
scratchpad and are not kept. It used 14 LPs and one public fetch (UniProt REST, for the core subunit identities).
Its report is reproduced below, lightly reformatted; the fixes made in response are listed at the end and in
REPORT.md §10.

## 1. SHAs, builds, seed check — supported

`shasum -a 256` on the short copies and the `model/candidates/` files: baseline `d8d1bac78138…e535751` (as claimed),
candidate `ca43ae0c12e0…843d0`. `logs/control_d8d1bac7.build.json`: same 32 flags in the same order as the
baseline's `.build.json`, output `d8d1bac7…`, HEAD `ff7fc62`, with the R171 module present (`r171_atp_synthase_gpr.py`
`315646ab`, `main.py` `70693d3b`, `cli.py` `c023d389`) and `options.r171_atp_synthase_gpr=false`. Seed 0 and seed 1
`.build.json` differ from the control only by `--r171-atp-synthase-gpr`; both output `ca43ae0c…`. Source and data SHAs
in the build records match the current files (curation JSON `63a7c6ba`). Both candidate logs: "R171 ATP synthase GPR:
applied". Limit: PYTHONHASHSEED is not recorded in `.build.json` (stated in `build_commands.txt`); the seed values are
self-reported.

## 2. Whole-model diff — supported

Raw `diff` of the two XMLs: 2 hunks, the 7 added R171 note lines and the R171 `fbc:or` block replaced by a 12-gene
`fbc:and`. cobra semantic diff (reactions, metabolites, genes, groups, objective, compartments, model notes and
annotation): identical counts (2,319 reactions, 1,882 metabolites, 1,081 genes, 113 groups) and 7 differences, all
intended: R171 rule, notes and gene set, and the reaction lists of `YALI1B15733g`, `YALI1B28025g`, `YALI0D17490g`,
`YALI1E38131g`. The report's 6 are the same set without counting R171's gene set separately.

## 3. R171 rule and unlinked genes — supported

The rule equals the TASK §2 string exactly. The four removed genes are still defined (one `fbc:id` each) with no
`geneProductRef`; in cobra each has `[]` reactions (baseline `['R171']`).

## 4. Guarded step — supported (minor notes)

A faithful single-reaction version of `aars_thr_ser_split.py`: spec pinned (id, status, before / after rule, 16 genes
with before / after reaction lists); R171 name, bounds, stoichiometry, species, annotation, rule and notes pinned;
gene annotation and reaction lists pinned. Partial states rejected (rule neither before nor after; after rule without
notes; notes without after rule; conflicting notes; other notes changed). Stricter than the pattern: `reactions_after`
is checked for all 16 genes. Edits only R171's rule and notes; no other code references R171 or the four genes. Minor:
gene name and notes are not pinned (as in the pattern); the post-edit `RuntimeError` leaves edits in place outside a
cobra context (as in the pattern; unreachable unless cobra misbehaves); no flag dependency, acceptable because the
full before-state is pinned.

## 5. Test — supported

`pytest -v tests/test_r171_atp_synthase_gpr.py` (workflow-automation venv, from `platform/`): 6 passed, 51 subtests
(13 stale-input + 6 scope + 32 knockout cases). It checks the knockout logic before and after, the SBML roundtrip with
a second apply returning `already_correct`, context exit and the opt-in flag. Limit: the test's `state()` omits
metabolite annotation / notes, gene names and groups; for the real build the raw XML diff covers them. The full-suite
after-run was still running at audit time (unverified then).

## 6. Sweep spot check and matrices — supported

Own Gurobi setup at 1e-9, `problem.reset()` per solve, 14 LPs: WT 0.3116608216336117 on both; `YALI1F04433g` and
`YALI1M00260g` 0.31166 before, 0.034189979775491146 after (KO/WT 0.10970, R171 closed, knocked flux 0); `YALI1E38131g`,
`YALI0D17490g`, `YALI1B28025g` KO/WT 1.0 on both; unrelated `YALI1A00309g` 0.0 on both; all equal to `ko_compare.tsv`.
Matrices recomputed from the TSV (strictly below cutoff; 16 mtDNA genes excluded; unresolved excluded): baseline
96/219/35/699, 102/213/37/697, 106/209/42/692, 114/201/42/692; candidate the same except 15 % 123/192/42/692 (recall
0.390, precision 0.745, specificity 0.943, MCC 0.420); 16 unresolved each; 0 class mismatches with the TSV; growth
changes > 1e-6 exactly the 12 AND members; class changes only at 15 %; identical non-optimal sets (16 + 4 mtDNA);
coverage 1,612 listed, 326 in the model, 11 unresolved, 1,286 not in the model. Copied scripts byte-identical to the
earlier copies (`bb3118a6`, `8a464ff8`); the `ko_compare_two.py` docstring describes an older run (cosmetic).

## 7. Catalytic reasoning vs review §7.1 — partially supported

Faithful: quotes and attributions match §7.1; tiers marked recalled (verified in the review), "not catalytic"
inferred; j's call stated as a judgement; FP, recall, MCC and Yeast-GEM listed as not reasons. Core subunit
identities checked through UniProt REST from the model's annotations: Q6C326 ATP1 α, Q6CFT7 ATP2 β, Q6C338 ATP3 γ,
Q6C877 ATP16 δ, Q6C9B1 ATP5/OSCP, Q6C105 ATP4 b, Q6CFH9 ATP7 d, Q6C2V6 ATP14 h, Q6C9E6 ATP17 f, all as in REPORT §1.
Gaps: (1) the review's evidence against removal is omitted (hydrolysis not synthesis; g Cas9 −2.80; Liu lists e;
UniProt places j in the enzyme; j needed for full coupling); (2) the strict "loss closes the reaction" reading drops
e and j, but the retained b, d, f, h, OSCP meet it only through the consensus screens, which are the benchmark list,
so the +9 TP partly rests on non-independent data; (3) the model note says the monomer "is catalytically active"
where the measurement was hydrolysis.

## 8. Escher numbers — supported

From `escher/flux.json`: WT growth 0.3117 on both, KO-after 0.03419; R305 9.0077 → 18.3151; R304 4.5038 → 9.1576;
R1384 0 → 155.7016; R642 3.3435 → 4.6716; R694 3.0829 → 4.6538; R741 1.1617 → 4.5481; R815 18.2713 → 4.5402; O₂
4.66 → 12.584; R659 / r0210 17.997 / 17.79 → about 0. KO-before fluxes equal WT; 12 LPs; 17 fixed + 12 selected
reactions. The glucose-cap sentence in REPORT §6 is recalled from the review, under a "verified" heading.

## 9. Completeness vs AGENTS §3–§4 — partially supported

Present: SHAs, changed entities, commands, negative class, coverage, four cutoffs, development-reference wording,
the provable / static / observed / untested split, large outputs in the research root. Missing at audit time:
pytest results, AUDIT.md, the STATE entry, a statement that no gene name is verified, a method-improvement note.

## Defects listed and what was done

| # | Defect | Action |
| --- | --- | --- |
| 1 | Full pytest not shown | Done: REPORT §7, acceptance #5 met (324 → 330 passed, same 19 failed / 20 errors, only the 6 new tests differ) |
| 2 | AUDIT.md and STATE entry missing | Done (this file; `model/STATE.md` 2026-10-29 entry) |
| 3 | Evidence against removal omitted | Done: REPORT §1 |
| 4 | Retained members rest on the screens | Done: REPORT §1 and §9 |
| 5 | Note says "catalytically active" | Recorded as a limit in REPORT §9; not changed, because it needs a data change and a rebuild to a new file |
| 6 | No verified-name statement; source of subunit assignment | Done: REPORT §1 |
| 7 | Glucose-cap sentence recalled; old docstring | Sentence marked recalled in §6; docstring left as is (copied unchanged on purpose) |
