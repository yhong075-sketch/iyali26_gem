# Execution notes

Official predictor/BLAST checks were independently parsed once from saved files; no predictions or alignments were rerun. A first parser assertion incorrectly counted identities from SEG-masked X display strings; this software-only audit assertion failed before any result was written. Restoring underlying query/subject residues from the pinned FASTA gave the exact XML identity counts in all five HSPs (P43590 displayed213 vs original219; Q96WX8/CAC displayed189 vs original203; P12955 displayed129 vs original133; ABW84230.1 both83). The corrected checker passed.

Visual checks: Sarry Table2/3, Jalving printed85 (PDF91), original Huh full-field image. The small GFP clip was blank; full-field image contains the observed signal, and only the full image plus explicit project classification were used.

Substantive outcomes and exact inputs: independent_results.json. Main-agent source report was reviewed; core distinction of vacuolar detection versus single-protein protection and intracellular versus inferred cytosol is accurately stated. Suggested replacing 原生GFP with 染色体C端GFP融合蛋白 to avoid calling a tagged construct unmodified native protein.
