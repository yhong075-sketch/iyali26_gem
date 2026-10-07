# Independent computation audit

Reviewer: `/root/dipeptidyl_gene_audit`. Checked 2026-09-11, approximately 17:19 UTC. Scope: existing outputs for **YALI1B04274g**, DPP-IV-like candidate, and **YALI1B25603g**, DPP-III-like candidate; formal native gene symbols remain unverified. All conclusions concern predictions and computational records, not experimentally demonstrated native functions.

The reviewer independently parsed the original BLAST XML, fixed FASTA/JSON sequences, AlphaFold PDB/confidence/PAE files, original experimental PDB headers, extracted coordinate chains, US-align text alignments and saved transforms. Recalculation used Python standard library and NumPy 2.3.3. No new BLAST search, structure-comparison matrix, AlphaFold prediction, optimization, model/GPR/medium edit, or experimental job was submitted.

## Findings

**No material numerical discrepancies were found in the checked completed analyses.** The local nine-reference panel is a limited comparison, not a completed search of all Swiss-Prot. The remote search remains unresolved and is not counted as passed.

### Sequence identity and accession status

For both targets, the reviewer verified exact equality of the fixed NCBI FASTA, historical cached sequence, target manifest, AlphaFold metadata sequence fields, and all-residue AlphaFold PDB CA sequence. Lengths are 887 and 707 residues. Their independently recomputed sequence hashes are:

- YALI1B04274g / AOW01140.1 / A0A1D8N680: `849e5b7db8e5612aa78afd5528c9ec3a163678aa75642392713e3fecaed57798`.
- YALI1B25603g / AOW01944.1 / A0A1D8N8I3: `dd0e1f77ec2785e2211d05098e11ec0b59180c7adc4a08040f3823e0d898c374`.

The newly retrieved UniProt JSON records say `Inactive`, `DELETED`, reason `Not part of a reference proteome`. This verifies database-entry status and **does not demonstrate deletion of either biological locus or inactivity of either enzyme**. Historical annotations must retain their historical provenance.

### BLAST arithmetic

All eight returned HSPs were checked character by character against the target and reference sequences. Identical nongap columns were recounted; query/reference coverage was recomputed from the union of covered residue positions; endpoints and the 22 saved reference-to-target functional-site mappings were independently reconstructed. All assertions passed.

| Target | Reference accession | Identical / aligned columns | Identity | Query coverage | Reference coverage |
|---|---|---:|---:|---:|---:|
| YALI1B04274g | Q96VT7 | 357 / 753 | 47.4104% | 82.2999% | 83.0189% |
| YALI1B04274g | P18962 | 362 / 833 | 43.4574% | 90.6426% | 97.6773% |
| YALI1B04274g | P33894 | 250 / 723 | 34.5781% | 78.2413% | 76.3695% |
| YALI1B04274g | P27487 | 233 / 744 | 31.3172% | 80.9470% | 93.0809% |
| YALI1B04274g | P42658 | 217 / 760 | 28.5526% | 81.8489% | 82.8902% |
| YALI1B04274g | Q8N608 | 197 / 679 | 29.0133% | 74.2954% | 81.4070% |
| YALI1B25603g | Q08225 | 342 / 710 | 48.1690% | 94.2008% | 98.7342% |
| YALI1B25603g | Q9NY33 | 296 / 692 | 42.7746% | 93.2107% | 93.2157% |

Reference names, systematic identifiers, functions and evidence qualifiers are in `REFERENCE_AUDIT.md`. The nine sequences were chosen in advance, and their reported E-values apply to that limited search. Missing hits are not an exhaustive exclusion of other functions. Reviewed status alone does not make every functional feature experimentally demonstrated.

### Functional residues and AlphaFold confidence

The complete 22-row mapping was checked. In YALI1B04274g, multiple active DPP-IV-family references consistently map the catalytic triad to **S749/D826/H859**. Human DPP4 N-terminal-recognition residues map to **E325/E326**. Their pLDDT values are respectively 98.19/95.75/95.00 and 86.88/90.62.

In YALI1B25603g, DPP-III zinc-binding positions map to **H460/H465/E519**, and the proposed catalytic-base position to **E461**. Their pLDDT values are 90.69/84.31/89.25 and 83.62. These are conservation and prediction statements; they do not independently measure catalysis or metal binding.

The reviewer recomputed all reported overall pLDDT means and below-threshold counts, domain means and fractions, overall/domain internal PAE averages, maximum PAE, and the full mapped-site PAE submatrices. Both confidence sources agree with the PDB per-residue values within rounding. Results:

| Target | Mean pLDDT | Residues below 50 | Residues below 70 | Mean PAE |
|---|---:|---:|---:|---:|
| YALI1B04274g | 85.01226607 | 125 / 887 | 152 / 887 | 11.47379472 Å |
| YALI1B25603g | 89.19200849 | 34 / 707 | 49 / 707 | 8.69901110 Å |

YALI1B04274g's predicted transmembrane interval 92–113 has mean pLDDT 68.8059; its predicted catalytic domain 681–884 has mean 97.5566. The lower-confidence N-terminal region limits structural/localization conclusions. High confidence in the catalytic core does not establish a particular organelle. AlphaFold metadata confirm existing model v6 records, created 2022-06-01, using the stated AlphaFold Monomer v2.0 pipeline; these were reused, not generated anew this round.

### Structural output and geometry

All eight raw US-align alignments were independently checked for reported length, RMSD, both TM-score normalizations, and coverage denominators. Original coordinate-chain lengths and all source/extracted-chain hashes in the eight-entry chain manifest match. The reviewer reconstructed paired CA coordinates from the text alignments and independently performed SVD/Kabsch fitting; all RMSDs agree with the reported values within 0.0051 Å rounding and with `geometry_verification.json` within numerical precision.

| Target / reference PDB | Aligned CA pairs | Independently fitted RMSD |
|---|---:|---:|
| YALI1B04274g / 1NU6 | 703 | 2.14387705 Å |
| YALI1B04274g / 1XFD | 702 | 2.60897819 Å |
| YALI1B04274g / 3DDU | 597 | 4.37967320 Å |
| YALI1B04274g / 3CSK | 395 | 8.50606310 Å |
| YALI1B25603g / 3CSK | 668 | 3.22105635 Å |
| YALI1B25603g / 3FVY | 669 | 3.71782684 Å |
| YALI1B25603g / 5E33 | 654 | 5.20882349 Å |
| YALI1B25603g / 1NU6 | 372 | 7.91676217 Å |

The exported TM-optimized transforms yield slightly different RMSDs (for example, 5.23927 Å rather than 5.20882 Å for the 5E33 pair). This distinction is preserved and resolved by independent least-squares fitting; it is not evidence of changed input coordinates. TM-score numbers were independently matched to raw tool output, not independently reoptimized. Strong resemblance to the inactive DPP6 control 1XFD prevents a fold-only activity claim.

Original PDB `SEQADV` lines were checked directly:

- **3CSK:** C130S is an engineered mutation; G458D is marked **CONFLICT**, not engineered. The latter lies close to the catalytic motif and must remain disclosed.
- **5E33:** C19S, E207C, **E451A**, S491C, C519S, C654S are engineered changes. E451A makes this an inactive catalytic-base construct; its closed conformation cannot establish target activity.
- **3DDU:** Q307H, W319T, I360T, T374I, S667N are engineered changes. It is a modified specificity-class control, not an unmodified native reference.

The main agent's initial geometry parser failure due to alternate-location residue labels and its subsequent parser correction are disclosed in the handoff. This reviewer used a separate parser selecting CA atoms with blank/A alternate location and did not rely on that failed parse. The successful independent geometry audit used unchanged source coordinates.

## Coverage and decision

Eight scoped claim groups: input/sequence identity; database status interpretation; HSP identity/coverage; residue mappings; AlphaFold confidence/provenance; structural output/geometry; reference-construct limitations; remote full-database completion.

**8 total groups | 7 completed groups audited | 7 supported within computational scope | 1 unresolved | 0 contradicted | 1 unchecked/incomplete remote result.** Detailed checked units include 8 HSPs, 22 residue-mapping rows, 2 complete target confidence sets, and 8 structural comparisons. No completed full-Swiss-Prot search is asserted.

The checked evidence supports calibrated DPP-IV-like and DPP-III-like functional predictions. It does not confirm native substrate specificity, synthesis of all four dipeptides, intracellular localization, biological essentiality, or R795 coupling. No GPR or model change is accepted by this audit.

## Checked output identities

- `target_sequences.json`: `28e104129346fe6ed8ba090e13078c17226b38459557b0387e6af1a3b3021902`
- `local_blast_results.json`: `01091314d20ad67c497344be9e0c01e34ee90f4f9bd7b30c000d9ca34a5af502`
- `residue_mapping.json`: `e6e4060e6393c285460ccd13db4889ad0cbb31dd1ae67600e2a34c8397d5d2c8`
- `alphafold_quality.json`: `44c3e5dce1f645bc96bb91157f4c4e04f9a481589d367b59c69f4395a83b9f13`
- `chain_manifest.json`: `84711cf0d041e9e43c23a7951e72fe18b80c588c3c52874c9902aec18be7a415`
- `structure_results.json`: `960d8d3f43befd9605d23c59d745dc420dd4bd6144def286a256309a52f25ba8`
- `geometry_verification.json`: `da85da886d000a4d2f5c0ef4ecf6575fca951ff7d8455008732902272a690c11`
- `verify_geometry.py`: `29b16959ece6b855984868d59ea4a55dff81ba3609311835ee8e94ee18bc75f8`
