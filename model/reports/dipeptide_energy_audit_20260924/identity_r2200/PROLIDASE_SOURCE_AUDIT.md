# Independent audit: R2039 prolidase candidate

Audit date: 2026-09-24. Scope: independent reading and recomputation of the provenance agent's candidate identity, family annotation, compartment and exact-substrate claims. Only current-workspace files and public primary/official sources were read. No model, GPR, medium or direction change, LP, new alignment, structure prediction or cluster job was performed.

**Verdict: YALI1E16433g is a defensible priority assay candidate, not an accepted R2039 GPR.** Its native established protein name is unverified; the supported label is **predicted M24B prolidase-like peptidase**. The sequence identity and CDD annotation claims pass independent checks. Free Gly-Pro activity and native compartment remain unverified for the W29 protein. The vacuolar model reaction cannot be assigned solely from these annotations.

## Identity and evidence level

I reparsed all three actual GenPept ORIGIN fields in `../provenance/E16433_E13464_genpept.gb`, without importing the other agent's numerical results, and compared the complete sequences with both cached UniProt records. Results and input hashes are retained in `PROLIDASE_independent_checks.json`.

| Claim | Independent finding |
|---|---|
| W29 accession and length | AOW05368.1 and XP_503902.2 are identical, 454 aa; both source features specify CLIB89(W29), ATCC20460. |
| W29 locus | AOW05368.1 CDS specifies YALI1_E16433g and complement(CP017557.1:1641943..1643307). XP_503902.2 has old locus YALI2_E01492g and GeneID2912030. These are source-record labels, not inferred mappings. |
| UniProt link | The cached A0A1D8NIA1 amino-acid sequence equals both W29 versioned sequences. An inactive current UniProt entry does not imply gene absence. |
| CLIB122 comparison | CAG79495.1/YALI0_E13464g and cached Q6C610 are identical to each other, also 454 aa, but differ from W29 at position17: W29 Q, CLIB122 K. No other substitution was found. |
| Family/domain | All three GenPept records annotate Prolidase, residues159–427, cd01087/CDD238520. W29 also has AMP_N residues6–131. These are domain-based annotations, not biochemical measurements. |

W29 sequence SHA256: `d7a158663b2c1d752f1b14bec8cdf30f135470992eb8478b26521964cb159b59`. CLIB122 sequence SHA256: `8ea6d41aa3e2bd27db70355dd69ff9db117a199dc37f8b83e9a223f8857828c8`. Hashes cover uppercase amino-acid letters only. The CLIB122 protein is a homologous candidate with unverified native name/function, and any structure using that sequence must not be described as an exact W29 model. This audit establishes no functional effect of Q17K.

Source locators: [AOW05368.1](https://www.ncbi.nlm.nih.gov/protein/AOW05368.1), FEATURES/source, Region and CDS; [XP_503902.2](https://www.ncbi.nlm.nih.gov/protein/XP_503902.2), COMMENT and FEATURES; [CAG79495.1](https://www.ncbi.nlm.nih.gov/protein/CAG79495.1), FEATURES/source and ORIGIN. Actual cached official GenPept text was inspected on the audit date. XP_503902.2 is explicitly provisional RefSeq and a conceptual translation; the product remains uncharacterized. The CDS homology note to the S. cerevisiae locus is annotation evidence, not a newly executed comparison or a native W29 assay.

## Exact free substrate versus a DPP surrogate

[IUBMB EC3.4.13.9](https://iubmb.qmul.ac.uk/enzyme/EC3/4/13/9.html), Reaction and Comments, supports the *class-level* interpretation that a prolidase hydrolyzes the bond within a free Xaa-Pro dipeptide. It does not certify that this particular W29 sequence has that activity. The entry's localization statement concerns animal enzymes and cannot locate a yeast candidate.

[IUBMB EC3.4.14.5](https://iubmb.qmul.ac.uk/enzyme/EC3/4/14/5.html), Reaction, describes release of an N-terminal dipeptide from a longer substrate. Thus a Gly-Pro-pNA assay measures cleavage after the Pro residue to release a chromophore, whereas R2039 requires cleavage **between Gly and Pro**, producing free Gly and Pro. These assay outcomes cannot substitute for one another. This distinction is supported independently by the official enzyme definitions; no Yarrowia activity is inferred from the shared letters Gly-Pro in the assay name.

The authors' [2023 original preprint, v1](https://www.preprints.org/manuscript/202310.1370) was read directly, especially §§2.6–2.7 and Table2. It reports recombinant A. phoenicis ApXPD activity on free Gly-Pro using liberated-proline detection by a modified ninhydrin assay; Table2 reports 86.3±1.7% relative activity versus Lys-Pro. Longer peptides were separately tested. This supports a suitable catalytic-class/assay precedent, **not** W29 activity, localization or transferable kinetic/pH parameters. This is a preprint version, not independent replication; the linked final Fermentation article (doi:10.3390/fermentation9110978) returned429/inaccessible to this audit's browser. The final published table was therefore not fully rechecked here. The cached preprint download attempt also returned403, although its full browser-rendered primary text was available and read. These access failures are not biological negative results.

## Compartment and GPR decision

The inspected cached A0A1D8NIA1 and Q6C610 records contain neither a SUBCELLULAR LOCATION comment nor a cellular-component GO entry. W29's inspected functional annotations are computational (SMART/InterPro/ARBA/TreeGrafter); GenPept has no experimentally established native compartment. **Absence of these annotation fields is not proof of cytosolic, secreted or vacuolar localization.** Family membership, the source organism of a homolog, and artificial secretion during recombinant expression cannot fill this gap.

Keep E16433 as a sequence-defined prolidase-like candidate for R2039 investigation. Acceptance of its GPR requires exact W29 protein activity on unmodified Gly-L-Pro with appropriate inactive controls and product confirmation; acceptance in the stored vacuole additionally needs native localization or a supported transport/access mechanism. Do not add another candidate as an OR isozyme, infer a complex AND, move the reaction, or accept vacuolar activity from the current evidence. The provenance report correctly retains these limitations.

## Audit record

The three sequence hashes, single-residue difference, domain coordinates and absence of inspected compartment annotations were recomputed and assertion-checked in this audit. An initial parser expected ORIGIN with no trailing spaces and stopped before writing results; the corrected parser accepts the actual GenPept whitespace, then all assertions passed. Official IUBMB pages were retrieved on2026-09-24 and saved with URLs/hashes in `PROLIDASE_source_manifest.json`; full source-file hashes are in `file_manifest.json`. This independent audit does not add a new claim of native catalytic validation.
