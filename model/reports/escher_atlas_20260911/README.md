# iYali26 pathway and gene atlas — 2026-09-11

Open `index.html` in a browser. Keep it beside `escher.min.js`; no server or internet connection is required. Scroll vertically through the map and horizontally for wide maps. Use + / − to zoom, search reaction or gene IDs, and click a reaction for its complete stoichiometry, bounds and GPR.

`viewer_template.html` is the generation template, without embedded model data. Opening it now automatically redirects to its completed sibling `index.html`; this also applies to the issue viewer template.

Open `issues/index.html` for **11 audited model issue cases**, covering Boolean GPR transmission, signed alternative flows, duplicate reactions and compartment supply. Switch between WT, target KO and available paired controls to replay **22 saved LP witnesses**. These witnesses belong to the earlier published tRNA model (d274), not a new run of the current b00 working hypothesis. `issues/version_binding.json` records their exact relationship. The external iYLI647 panel is a separate model and partial condition match: its structural certificates and saved growth ratios do not supply a full external flux vector. This is this project's analysis method, not an attributed repair workflow of the external authors. Issue maps retain complete stoichiometry and can be exported as SVG or Escher JSON. Raw sources, scoped certificates and source audits are included.

Open `review/index.html` for the **full-model static review of the d274 published tRNA baseline under recorded SD-Leu / PO1f constraints**. Its three lenses include all 1,074 original gene entries (plus a separately identified runtime plasmid marker), 2,315 reactions, 1,877 species pools, 15 proportional duplicate groups, 135 two-step identities and 628 complete-row proof records. Filters distinguish historical FN calls, retained GPRs, reverse observations, direction-curation conflicts, transport candidates and sufficiently proved zero flux. There are 69 initially closed reactions and 643 additional proved-zero reactions; the remaining 1,603 are not proved blocked by this method. These are not 643 biological defects or a complete FVA. Thirteen priority themes contain 23 Escher maps with complete source stoichiometry, scoped proof captions and optional replay of the 22 archived full vectors. Full supply/consumption rows and prior-zero proof dependencies remain browsable. No flags are transferred to the distinct b00 candidate model. Historical observations, current static proofs and unresolved biological interpretations are separately labeled.

This update uses **model_metadata_trna_vatpase_and_hypothesis.xml**, the latest authorized **non-lipid-unlump working hypothesis** identified for this task. The three-gene V-ATPase AND rule remains a hypothesis; this is not a formal model release. Original lumped lipid reactions and tRNA-coupled biomass are retained.

All **2,315 reactions, 1,877 metabolites and 1,074 genes** are represented across 112 original subsystems, 235 module pages and 4 additional overviews. **1,062 genes** have GPR connections; the other **12** appear in `unassociated_genes.html`. Experimental positives comprise 314 on reaction maps and 8 without GPR connections.

| Marker | Meaning | Limit |
|---|---|---|
| E | Experimental positive in the saved source table; 322 genes | Positive-only development reference. Absence means unknown, not experimentally nonessential. |
| P | Predicted essential in the saved full screen of this exact model; 101 genes | PO1f / SD-Leu; unrounded KO/WT strictly < 0.10. Model prediction, not experimental validation. |
| E and P | Both labels apply; 73 genes | Evidence types remain distinct. |
| Gray dashed reaction | Source bounds are [0,0]; 68 reactions | Explicitly closed; does not mean experimentally inactive. |
| Teal reaction | Source bounds allow nonzero flux | Network blockage or physiological activity is not established. |

The matching saved screen contains 1074 optimal, finite, nonnegative KO results. No screen was rerun for this atlas. Its PO1f / SD-Leu conditions are attached to P labels; the static reaction map retains original SBML bounds without a culture or strain overlay. The previously identified closed ATP-generation problem remains unresolved, which limits biological interpretation.

Compared with the 2026-09-09 atlas, this model adds R1172 and SPHPL and updates source stoichiometry, direction bounds and GPR rules. `model_update_audit.json` records the exact differences and source identities. The old atlas is retained in its original dated folder. Old historical H labels have been replaced by the matching saved screen's P labels.

Arrows show model direction, not flux. Repeated IDs denote the same compartment-specific pool. The backbone display hides some participants, including possible co-substrates and branch products; **Show all species** expands them. Full source stoichiometry is retained in reaction details and Escher JSON. Module maps show complete gene rules by default; overviews initially collapse them. OR represents model alternatives, not proven biological interchangeability. A reaction associated with an essential-labelled gene is not automatically an essential reaction.

Gene details separate systematic IDs, established names, candidate protein functions and evidence status. Sequence and existing AlphaFold support remain prediction-level evidence. No new functional validation or structure prediction was performed.

`module_index.html` lists every SVG and Escher JSON. `svg/` contains 239 standalone vector maps. `model.json` is a display extract and omits solver-relevant SBML content; use the frozen `inputs/model.xml` for source identity. `annotation_sources.json`, `provenance.json`, `validation.json`, `render-check.json` and the audit files record evidence and verification scope.

This task performs local visualization and source comparison only: **0 solver calls**, no model/GPR/medium/experimental-label edits, no publication or cluster jobs. Escher 1.8.1 and its license are included. Source models and derived maps retain the repository's CC BY 4.0 notice.
