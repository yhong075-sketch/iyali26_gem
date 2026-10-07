# Protein-function prediction for two W29 dipeptidyl-peptidase candidates

User authorization: 2026-09-11 “进行蛋白功能预测”, following the two YALI1 candidate lookup. This supersedes the prior lookup-only scope; it does not authorize model/GPR/medium edits. Applicable govern-agentic-research, gene-identity-function and project AGENTS remain in force.

Question: what enzyme class, catalytic competence, substrate-family preference and compartment evidence support or contradict YALI1B04274g and YALI1B25603g? Are either credible candidates for release of any of the four target dipeptides? Do not presuppose they must support a missing model reaction.

Inputs: fixed public W29 proteins A0A1D8N680, sequence v1, 887 aa, SHA256 849e5b7db8e5612aa78afd5528c9ec3a163678aa75642392713e3fecaed57798; A0A1D8N8I3, sequence v1, 707 aa, SHA256 dd0e1f77ec2785e2211d05098e11ec0b59180c7adc4a08040f3823e0d898c374. Verify current primary sequence records and preserve conflicts; do not silently substitute strain/sequence.

Plan: verify identities and retrieve primary annotations; search standard BLASTp against reviewed proteins and compare bounded experimental references plus close paralog/noncatalytic controls; map domains and experimentally supported catalytic residues; retrieve exact-sequence AlphaFold DB structures and PAE/confidence, then align to experimental structures; have an independent source/numeric audit and report calibrated predictions.

Bounds: exactly two targets. Up to one remote Swiss-Prot BLAST batch submission (two queries) with 15-minute response budget; existing local BLASTp against a fixed reviewed reference panel is a supported fallback, clearly labelled limited search. Local analysis 1 CPU, each tool max 300 seconds, at most 12 structure comparisons. Public database retrieval only, no experimental execution, no model/curation/medium/label changes, no git publication.

AlphaFold: prioritize matching existing models. If absent/mismatched, user's standing HPCC authorization allows at most one monomer_ptm prediction per missing target, serially, 1 GPU/8 CPUs/64 GB/4 h per target, no retry/grid, after current environment verification. Record exact runtime/input/job details. Do not assume structures prove activity, substrate specificity, localization or native essentiality.

Failure/stop conditions: sequence disagreement restricts that claim, retrieval errors are recorded, no silent target replacements. The least-supported endpoint is exact native substrate/compartment specificity; leave it unresolved when homology/structure cannot decide. No criterion based on improved essentiality matching. Outputs: source/sequence manifest, BLAST raw and parsed results, residue/domain table, AlphaFold confidence/PAE and structural alignment metrics, Chinese report, independent audit.
