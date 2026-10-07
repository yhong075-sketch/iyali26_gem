# iYali26 pathway atlas — 2026-09-15

Open `index.html` beside `escher.min.js`. No server or internet connection is needed for the pathway viewer. `viewer_template.html` redirects to the completed page. Scroll the page; use + / − for map zoom, and Show all species for complete reaction participants. Click a reaction for full stoichiometry, bounds, GPR and gene evidence. Export SVG or Escher JSON directly.

Source: **model_metadata_trna_r153_merged.xml**, the authorized September 15 ordinary non-lipid-unlump model: **2314 reactions, 1877 species and 1073 model gene entries**. The original lumped lipids are retained. The older V-ATPase AND atlas is a distinct archived hypothesis. This visualization does not accept or change a scientific model.

There are 239 maps: 235 subsystem module pages and four overviews, covering 112 source subsystems. All source reactions, stoichiometry and GPRs are checked. Of the gene entries, 1061 have source GPR connections and 12 appear in `unassociated_genes.html`. Entry counts include mitochondrial and placeholder IDs; they are not a genome-wide census of verified native genes.

**P: 116 / 1073 model gene entries screened.** These labels use the completed screen of this exact model under the recorded PO1f / SD-Leu conditions. Essential means the unrounded ratio **KO/WT < 0.15**; equality is not essential. WT growth: 1.8718823069402954 h⁻¹. The runtime plasmid marker is excluded from the original-gene screen. Raw results, numerical checks, conditions, software identity and original failed-attempt records are retained in the screen directories. No FVA or extra condition matrix was run.

**E: 322 experimental positive entries in this model.** Original positive labels and source locators are preserved from the saved project reference. This development reference is not an independent validation set. Unmarked genes have unknown experimental status. A gene's P or E marker does not establish reaction essentiality.

Gray dashed reactions have original source bounds [0,0] (68 reactions). Teal means bounds allow nonzero flux; actual activity and full-network blockage are unknown. Reaction colors do not incorporate the screen's medium/strain overlay. Arrows show source direction, not a newly computed map flux.

R153 retains a user-assigned GPR with conflicting protein-function evidence; the source notes do not establish native synthase activity. Gene function, model assignment and essentiality predictions remain distinct. No new sequence/structure prediction or biological validation was performed.

The previous dated atlas and d274 static-review certificates remain separate historical artifacts. Any historical navigation links require the original adjacent dated directory (or the local archive server); old certificates and flux vectors are not applied to this model. See `model_update_audit.json` for the scoped display-extract comparison.

The final data audit, rendering checks, source hashes and source XML are included. The accepted screen used 1 WT + 1073 independent KOs (1074 primary solves). The stopped earlier attempt used 41 primary solves, for 1115 across both attempts; it remains separately identified and is not counted as a completed screen. This package does not modify source models, GPRs, medium, experimental labels or Git state. Escher 1.8.1 and its license are included; repository model/map licensing is retained.
