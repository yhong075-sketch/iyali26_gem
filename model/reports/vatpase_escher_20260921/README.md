# V-ATPase candidate: English Escher view

Open **[the interactive Escher page](vatpase-escher-en.html)**. Keep `escher.min.js` in the same folder. Select WT or any of the 14 saved single-gene knockout results; click a reaction for its exact equation, GPR, bounds and saved flux. The viewer uses saved results and does not run an optimization.

`viewer.html` also opens the finished page. The build template is `viewer.html.in` and is not a browser entry point. If an earlier `viewer.html` tab showed an empty map and empty result selector, refresh that tab: the former template path now redirects to the populated page.

**Main difficulty:** the candidate GPR now propagates shared-component knockout to R794 closure, but the current model can still reach WT-level growth with both V-ATPases at zero flux. R795 was already closed. Changing the GPR alone therefore does not make these genes essential under this screen's conditions.

## What is still missing?

| Evidence gap | What would resolve it? | What the current data establish |
|---|---|---|
| Growth-linked acidification | A measured, condition-specific acidification requirement and a supported mechanism linking it to growth | The complete saved WT flux remains feasible after every target single knockout. A specific missing biological task has not been established. |
| Proton transport energetics | Physiological direction and energy coupling for R1159 | The stored reaction permits reversible cytosol–Golgi H+ transfer without ATP. Its saved WT flux is zero; it is not an observed active bypass. |
| a-subunit localization and substitution | Native Golgi/vacuolar localization and same-condition complementation | The model's OR relation is a test hypothesis, not demonstrated native interchangeability. |
| W29 F-subunit identity | Resolve the old CLIB122 sequence versus W29 locus conflict and confirm an intact, expressed native protein | Native W29 assignment remains unresolved. |
| Native complex assembly | Complete native stoichiometry and orthogonal assembly/function evidence | Previous AlphaFold predictions do not establish catalysis, localization or redundancy; the reviewed B branch has a webpage summary but no audited raw structure files. |

All named gene identities, candidate functions and evidence status are provided in the interactive page. No native name has been verified for the 14 targets. No arbitrary proton sink, new transport reaction, or extra AND relation was introduced to force essentiality.

## Exports

- [WT overview PNG](vatpase-escher-en.png) and [editable vector SVG](vatpase-escher-en.svg).
- [Shared-component KO example PNG](vatpase-shared-ko-en.png) and [SVG](vatpase-shared-ko-en.svg). The selected V1-D candidate and its evidence status are explicitly labeled in the figure and gene table.
- [Native Escher map JSON](map.json), retaining exact GPRs and stoichiometry for R794, R795 and R1159.

This is a selected local reaction view, not the complete ATP or proton network. Repeated metabolite nodes represent the same compartment pool. Arrows indicate stored directionality, not positive flux. All three drawn reactions have zero flux in the saved WT solution.

## Verification and provenance

The source is the independent candidate screened on 2026-09-21, under SD-Leu/PO1f with the unrounded KO/own-WT <15% criterion. Both baseline and candidate have 116 essential calls; all 14 target single KOs retain approximately 100% growth. These are model predictions, not native experimental findings.

[Provenance](provenance.json) records the full source hashes and bounded visualization plan. [Independent audit](AUDIT.md) checks the source claims, three reaction definitions and 15 saved-result scenarios. [Browser checks](render-check.json) verify visible stoichiometric coefficients, all scenario closures, reaction details, layout dragging and export preservation. Initial browser startup, inline-library and status-ID issues were fixed before delivery; subsequent visual overlap corrections did not change scientific data. No new LP, model/GPR/medium edits, AlphaFold predictions or cluster actions were performed for this figure.
