# iYali26 pathway and gene atlas

Open `index.html` to browse offline. Keep it in the same folder as `escher.min.js`. The atlas opens on the glycolysis module with reaction IDs, complete gene AND/OR rules and essentiality markers. Search reaction names, reaction IDs or gene IDs; switch modules, zoom, inspect all participating species, and export SVG or Escher JSON.

This atlas uses the project's registered latest non-lipid-unlump **provisional reference**, covering all **2,313 reactions, 1,877 metabolites and 1,074 genes**. The original **112 subsystems** are arranged into **235 module pages**, with **4 additional overviews**. Original lumped lipid reactions are retained. No lipid-unlump, C7 or isolated CoQ9 build candidate is incorporated.

Of these genes, **1,061** connect to reactions through the original GPR rules. The **13 genes without GPR links** appear separately in `unassociated_genes.html`; no reactions have been invented for them. This group contains 10 experimental positives, so 312 experimental positives occur on reaction maps and 10 in the unassociated-gene table, for a combined total of 322.

## Reading the maps

| Marker | Meaning | Limits |
|---|---|---|
| `E` after a gene | Matches the saved experimental consensus-positive reference; 322 genes in the model intersection | Not newly established independent validation. Absence from the positive-only table means unknown, not experimentally nonessential |
| `H` after a gene | Historically predicted essential; 96 genes with KO/WT strictly below 0.10 | Refers to the historical executed model under PO1f/SD-Leu conditions, not a new full screen of the map's source model. Anomalies before historical normalization have not been fully excluded |
| `E/H` | Both labels apply; 67 genes | The two evidence types remain distinct |
| Gray dashed reaction | Both bounds in the source reference equal zero; 68 reactions | Explicitly closed. This differs from FVA blocking, measured inactivity and zero flux in one solution |
| Teal reaction | Saved bounds permit nonzero flux | Does not establish network feasibility or physiological activity. Complete FVA-blocked status remains unknown |
| No GPR | The source model has no gene rule | No gene is invented, and the absence of a rule does not establish a spontaneous reaction |

Matching existing replay results are available for only 6 genes in the current reference, including 2 below the 10% threshold. These results are identified separately in the details and cannot establish a full gene-screen reproduction. The historically executed XML and the map's reference XML have different identities; historical H labels always retain their condition and source qualification.

Arrows show model directions, not simulated or experimental flux. Repeated metabolite IDs denote the same pool within the same compartment. The backbone view simplifies the display of some participants, potentially including co-substrates and branch products, not only cofactors. Use **Show all species** to expand them. Reaction details and Escher JSON always retain complete stoichiometry.

Module maps expand the complete source GPR by default. Featured overviews initially collapse it; use **Genes / GPR** to expand. AND/OR expressions retain the source parentheses. OR means alternatives in the model, not experimentally demonstrated interchangeability. A link to an essential-labelled gene does not establish that the reaction itself is essential.

Most source gene names are placeholders, so maps use systematic IDs. Details state whether names and protein functions are verified. Existing CoQ annotations retain their evidence scope, candidate identity and source; displaying them does not upgrade them to experimental validation.

## Included files

- `index.html`: offline interactive atlas, rendered with Escher 1.8.1.
- `module_index.html`: all module SVG / Escher JSON links, with gene and status counts for each page.
- `unassociated_genes.html`: the 13 source-model genes without GPR links and their existing annotations.
- `svg/`: 239 standalone vector maps with visible gene markers and legends, suitable for zooming and document layout.
- `maps/`: standard Escher map JSON files. E/H markers belong to the atlas display layer; original status annotations are saved separately in `annotation_status.json`.
- `viewer-preview.png`, `glycolysis-module.png`, `inactive-module.png`: interface and module examples.
- `model.json`: a model extract for Escher display. It lacks the full SBML annotations and objective definition and must not replace the source SBML as a solver input.
- `provenance.json`, `validation.json`, `render-check.json` and audit documents: input identities, sources, check scope and results.

This task only reads the model and existing evidence and generates layouts and annotations. Solver calls: 0. No model, GPR, culture condition or experimental label was modified, and nothing was published online.

Escher usage follows the [official JavaScript API](https://escher.readthedocs.io/en/latest/javascript_api.html). Its license is saved in `ESCHER-LICENSE.txt`. The project model and derived maps retain the repository's CC BY 4.0 notice. Full scientific sources and scope are recorded in `annotation_sources.json`.
