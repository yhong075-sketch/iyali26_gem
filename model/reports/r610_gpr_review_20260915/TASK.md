# R610 GPR scoped review

- Authorization: user's “审查R610的GPR”; review only. No GPR/model/medium/labels changed; no new metabolic optimization or cluster job.
- Input: current R1931-forward model, `model_metadata_trna_r1931_forward.xml`, expected SHA256 `d417f1de0425bc3336503b45a9498154b8a25a5d90d6f58e68b9d80cfc79297b`.
- Question: does the assigned W29 protein support the actual cytosolic ornithine + 2-oxoglutarate → glutamate + glutamate-5-semialdehyde reaction? Check specificity, cofactor, compartment, Boolean rule, and exact chemical duplicates.
- Evidence: versioned native sequence and official records; experimentally characterized reference enzymes; primary fungal biology; indirect domain/alignment/AlphaFold support. Do not equate model assignment or sequence similarity with native experimental validation.
- Sequence scope: one target versus three reviewed S. cerevisiae references: ornithine aminotransferase, GABA aminotransferase and acetylornithine aminotransferase. BLASTp 2.17.0+, one thread, maximum 300 seconds, no expanded database search. E-values apply only to this panel.
- Structure scope: reuse existing exact-sequence AlphaFold model and confidence/PAE; no new prediction and no unsupported structural-similarity claim.
- Static checks: current model identity, actual reaction/GPR, in-memory gene KO boundary propagation and exact/reverse stoichiometric duplicates; no solve. Reuse previous screen outputs only with their recorded model and conditions.
- Stop: deliver evidence-ranked GPR recommendation and unresolved native localization/assay limitations after independent source audit. If sequence identity conflicts, keep conflict explicit and do not change model.
- Deliverables: REPORT.md, identity/biochemistry notes, raw source records, bounded alignment results, model snapshot and independent audit. Native sequence proteome exclusivity and native physiological flux capacity are outside this scope.
