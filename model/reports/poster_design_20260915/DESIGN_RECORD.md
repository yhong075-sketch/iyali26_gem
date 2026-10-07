# Poster design record

- Deliverable: one 48 × 48 inch English poster. White background, thin rules, restrained teal/ochre colors, open scientific layout.
- Authors, affiliation, institutional logo, funding and contact were not supplied and remain explicit placeholders.
- Statistics refer to the saved 10 Sep 2026 screen. Coverage exclusion categories come from the 6 Sep 2026 ID review. They do not claim the performance of a later model.
- Text, formulas, workflow connections and two charts are native editable slide objects. The yeast is a conceptual illustration.
- Independent evidence review: 7 scoped claims audited and supported. Latest-model performance was not checked and is not asserted. No model, GPR, medium or labels changed. No new model solves.
- Chart counts: 322 + 1290 = 1612; 73 + 249 = 322; 73 + 28 = 101. Unlabelled genes are not experimental negatives. Coverage 19.98%; recall within the covered positives 22.67% at strict, unrounded KO/WT < 0.10.
- Primary literature: [Ramesh et al., Communications Biology (2023), Fig. 3b](https://www.nature.com/articles/s42003-023-04996-8).
- Saved project data: `artifacts/screen_test_metadata_trna_20260910/essentiality_summary.json`, `essentiality_per_gene.tsv`, `screen_predictions.tsv`; `results/workflow_foundation_20260906T211822Z/coverage_1290.tsv`.
- Small illustrative asset: [yeast_model_schematic.png](yeast_model_schematic.png). Built-in image_gen, no new empirical data represented.

## Illustration prompt

Use case: scientific-educational. Create a single isolated small illustration for a research poster: a Yarrowia yeast cell shown as a simple budding oval outline, with a sparse abstract metabolic network inside (8–12 small circular nodes connected by thin lines). It is explicitly a conceptual model schematic, not microscopy or a biological measurement. Style: very clean flat vector-like scientific line art, fine dark charcoal strokes, pure white background, restrained muted teal fill on some nodes and one pale ochre node. Cell interior white or extremely pale teal. No gradient, shadow, texture, glow, label, word, letter, number, logo, watermark, scale bar or surrounding objects. A single centered yeast cell, uncropped, generous clean margin; square high-resolution canvas. Match the restrained line-art conventions of a methods figure in a scientific paper. The art should remain crisp when printed at about 5 inches wide.

## Technical handoff

PPTX uses a 48-inch square canvas with native chart data snapshots. The PDF uses embedded Liberation Sans glyphs from the bundled converter and was visually checked against the slide. Native PowerPoint was not opened. Both full-page output and exact page size were checked; no clipping or unintended line breaks remain. Initial label wrapping in the converter was fixed in the chart text properties before finalization. The two final chart workbooks reproduce literal source counts rather than original experimental workbook formulas.


## Revised poster v6 — 15 Sep 2026 (local)

Supersedes v4 as the current visual draft. The user requested a more developed scientific design and more appropriate typography. The square poster now has a dominant set diagram, a larger yeast model illustration, two-column panels, a compact conclusion, and a rank plot of saved knockout growth ratios. Panel letters follow row reading order.

Times New Roman is used for headings and equations; Helvetica Neue for figure labels. PDF inspection confirmed embedded Times New Roman and Helvetica Neue fonts, with embedded Linux Libertine G fallback for a mathematical glyph. The earlier Liberation Sans note describes only v4. Both final deliverables retain the exact 48 × 48 inch size. All source counts and the 10 Sep 2026 scope remain unchanged.

The rank chart is a display-only ascending sort of 322 saved covered-positive ratios from essentiality_per_gene.tsv. Six decimal places are used for chart workbook storage; classification is verified from the original unrounded ratios (73 values < 0.10). Source values are not edited. No solver run was performed. All 322 plotted values and the full constant-threshold series were checked against the source. Office smoothing was explicitly disabled and both series use the same rank grid to prevent converter interpolation artifacts. The final PDF curve was visually inspected after these changes.

The PPTX contains two editable native charts with literal data workbooks, editable set and method diagrams, labels and formulas. A generated yeast illustration remains the sole raster artwork. Author, institution, logo, funding and contact placeholders remain. Validation receipt: build/poster_v6.validation.json; additional font, page and data checks: build/v6_final_checks.json. Native PowerPoint itself was not opened.

## Threshold revision v7 — 16 Sep 2026

User requested the 15% result. The poster now selects the existing cutoff_curve entry at 0.15: TP 79, FN 243, recall 24.53% within the same 322 shared positives. The threshold line, labels, proportional bar, recall formula and notes are synchronized. The full screen has 114 predicted-essential genes, comprising 79 labelled positives and 35 unlabelled. Unrounded saved ratios were checked; no new model solve was performed and historical primary_cutoff 0.10 in the source summary remains untouched.

48 × 48 inch layout and embedded Times New Roman / Helvetica Neue typography retained. Final PDF was rendered and visually inspected. Both native chart workbooks were verified against saved data. Earlier poster versions remain available. Validation and provenance: build/v7_15pct_final_checks.json and build/poster_v7_15pct.validation.json. The earlier independent-audit count is not represented as a fresh audit of the 15% result.

## CoQ9 example revision v9 — 16 Sep 2026

User requested replacement of the generic C workflow with our CoQ9 handling, titled C. Example. The revised panel presents R305 (complex III) with the original 1.5/1.5 proton equation, exact H and charge residuals of -2, and the optional full Q-cycle 2/4 candidate balanced under stored species. It states that corrected name/EC annotations are present but the candidate was not applied in the saved screen; physiological and network validation remain open. N/P sides and the cytosolic proton proxy are explicit. The 15% result and both chart datasets match the preceding poster at their six-decimal display precision.

A scoped independent reviewer opened the local source records and SHA-matching screen XML and recomputed exact balance without optimization. Coverage: total 8, audited 8, supported 7, unresolved 1, contradicted 0, unchecked 0. The unresolved proposition is completed physiological/network validation, which the poster explicitly does not claim. No gene-level functional assignment is introduced. Source hashes, output hashes and final checks are in build/v9_coq9_final_checks.json. No model, GPR, labels, medium or primary benchmark was changed.

Primary mechanism sources reopened by the main agent: IUBMB EC 7.1.1.8 (https://iubmb.qmul.ac.uk/enzyme/EC7/1/1/8.html) and Wieferig & Kühlbrandt, IUCrJ 10 (2023), 27–37, Introduction (https://journals.iucr.org/m/issues/2023/01/00/rq5008/index.html; doi:10.1107/S2052252522010570). These support the mechanism context; local exact residuals come from stored model chemistry, and balance alone does not select a unique coupling.

V8 was a private layout iteration; the finalizer preserved that file when a wording revision was attempted, so final delivery uses v9. Final full-page PDF and C-panel layout were inspected, native chart data retained, and embedded fonts/48-inch canvas checked. Older deliveries remain intact.

## CoQ9 mathematical-mechanism revision v11 — 16 Sep 2026

The user clarified that the intended example concerns the mathematical cofactor-dilution problem. Panel C now replaces the R305 proton-stoichiometry discussion with: (1) a Q9/Q9H2-only steady-state recycling illustration allowing equal positive cycle flux and zero net synthesis; (2) the runtime growth-coupled drain and summed pool balance v_syn + v_res = v_dil = alpha*mu; (3) an explicitly illustrative alpha=1e-4 mmol/gDW, mu=0.3 h^-1 calculation yielding 3e-5 mmol/gDW/h; and (4) the conditional zero-growth result when both net sources are closed and alpha>0. The illustration does not prove full-network feasibility or unlimited flux. Pool sizes/dilution are not encoded by Sv=0 alone.

The finite reserve is artificial and bounded/depleted by the dynamic runner, not by the source helper alone. Runtime sensitivity parameters are uncalibrated. The runtime reference model (bc2aac8f...1bee) is distinct from panel D's saved static screen (d274bad3...135a0). The 15% screen plots and counts are unchanged. No optimization or model/label/media mutation occurred.

Independent read-only audit: total 6, audited 6, supported 6, unresolved 0, contradicted 0, unchecked 0, limited to qualified algebra, implementation, unit arithmetic and experiment separation. The optional conditional reserve budget was audited but omitted from the panel for space. Exact helper and runner hashes match historical manifests. Source identities and final checks are retained in build/v11_coq9_math_final_checks.json.

Formulas are editable native text with seven proper Open XML subscript runs; the two dataset charts remain native. V10 was a layout iteration; a footer-rule collision was corrected in v11. Final PDF was rendered and reviewed, with embedded fonts and exact 48-inch square size confirmed. Prior deliveries remain intact.


## Expanded mathematical example v13 — 16 Sep 2026

User requested smaller A, B and D to allow a fuller explanation in C. Context panels now share the upper row; C spans the full width below and occupies about 60.6% of the main body height, excluding the title/banner and footer. Its allocated area is about 3.24 times v11. The 48 × 48 inch format and established Times New Roman / Helvetica Neue pairing are retained.

C presents three stages: the two carrier-only steady-state equations, an explicit biomass-linked oxidized-Q9 drain and total-pool balance, and conditional knockout/reserve consequences. It includes the illustrative alpha × mu calculation and the ideal finite-reserve budget and biomass ceiling. The constant positive alpha, fixed-volume, no-other-input/loss and nonnegative-reserve assumptions are explicit. R is an artificial inventory, alpha is uncalibrated, and other constraints can stop growth earlier. These equations describe an ideal budget rather than exact equality of all saved floating-point values. The runtime case remains distinct from D's static screen.

Independent source audit: seven scoped qualified claims, seven audited and supported, none unresolved/contradicted/unchecked; see build/v12_math_expansion_audit.md. This is a source/algebra review, not a new optimization or physiological calibration. No model, GPR, culture, label or benchmark inputs were changed.

The 322 saved ratios, 15% threshold and both native chart datasets match v11 at their six-decimal display precision; classification still uses the original unrounded ratios (79 TP, 243 FN, 24.53% recall). Thirteen equation subscripts are editable native text. The complete final PDF and enlarged C region were rendered and inspected; the 15% label was moved off its dashed line in v13. PDF page dimensions, embedded fonts, native-chart data and exact illustrative arithmetic passed checks. Final provenance: build/v13_expanded_C_final_checks.json; finalizer receipt: build/poster_v13_expanded_C_15pct.validation.json. Earlier outputs are preserved. PowerPoint itself was not opened.
