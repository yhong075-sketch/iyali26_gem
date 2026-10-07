# Latest-model pathway atlas and 15% essentiality

The English, scrollable atlas now uses the authorized September 15 R153-merged ordinary model: 2,314 reactions, 1,877 species and 1,073 original model gene entries. All reactions are represented across 235 module pages and four overview maps. The prior dated atlas and d274 review retain their historical version bindings.

## Completed same-model screen

Under the preserved SD-Leu / PO1f runtime context, the completed screen obtained WT growth **1.8718823069402954 h⁻¹** and valid results for all **1,073 independent single-gene knockouts**. Exactly **116** entries meet the requested **unrounded KO/WT < 0.15** criterion; equality does not qualify. The 1%, 5% and 10% comparisons from the same raw results give 76, 93 and 101 respectively; these are arithmetic comparisons, not additional simulations.

The accepted attempt made 1,074 primary LP calls (1 WT + 1,073 KO), used one worker with Threads=1, and finished in 28.28 seconds under the 600-second limit. All returned statuses were optimal, with finite nonnegative growth. The runner checked gene-to-GPR-to-bound propagation and full restoration after every knockout. Its largest retained KO mass residual was 1.331×10⁻¹² and bound violation 4.700×10⁻⁸, below the unchanged 10⁻⁷ check. Full flux vectors were checked during execution and were not saved for later independent replay.

The preceding attempt stopped after 41 primary calls because a GPR-unassociated gene lacked the owner required by COBRA's automatic context restoration. All 40 returned KO rows were preserved. An explicit target-flag restoration was regression-tested without optimization, then the complete attempt was separately authorized and run. The two attempts total 1,115 primary calls; the incomplete attempt is not merged into the completed result.

## Display and evidence scope

- **P** uses only the successful exact-model screen at 15%; no old-model P values are carried forward.
- **E** preserves the existing experimental positive reference and original locators: 322 model entries are labeled, and 751 remain experimentally unknown. Of the 116 predictions, 80 overlap the positive reference; the other 36 are unlabelled, not established false positives. The reference is a development/regression source, not independent biological validation.
- Source gene entries include placeholders and mitochondrial IDs. There are 1,061 entries with source GPR connections and 12 without. The additional runtime plasmid marker is excluded from the 1,073-entry screen.
- Gray dashed reactions denote the 68 original source-bound closures. Other colors do not establish physiological activity or complete FVA blockage. The map retains original source bounds; the screen's medium/strain overlay is recorded separately.
- Current reaction/GPR structure is drawn from the frozen merged source. R153's assigned function remains in conflict with retained protein-function evidence; a P label does not validate that assignment biologically.

No source model, chemistry, GPR, medium, experimental label or historical result was changed. No FVA, condition matrix, sequence/structure prediction, cluster job or Git publication was performed for this update.

## Verification

The builder checked Escher schema and complete source stoichiometry, GPR, bounds and participant coverage for every drawn reaction. Browser checks rendered all 239 maps and exported their SVGs, checked E/P labels against source records, tested complete GPR display, 15% text, reaction/gene search, wheel scrolling, narrow-screen layout, and HTTP/offline template entry. No page errors or external requests occurred. The actual in-app page also displayed 116/1073, the current merged source, strict 15%, and current reaction details.

See `final_audit.json` for independent data/source checks and their limits, `render-check.json` for browser checks, `screen_retry/` for the complete calculation, and `screen/` for the retained incomplete attempt. Input, code and output hashes are recorded in those evidence files. Matching recorded software and source files does not reconstruct the complete historical dirty environment.
