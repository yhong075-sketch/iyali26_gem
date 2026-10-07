# R1889 GPR correction

2026-10-05. R1889's blank GPR is filled in the [corrected E5 candidate](E5_r1889_candidate.xml) with an evidence-supported **partial complex-I dependency rule**. The same correction is integrated into the shared build pipeline after final metadata selection, without a new opt-in flag. The original models and historical E5 remain unchanged.

```text
YALI1B26679g and YALI1D07089g and YALI1F09003g and YALI1F22993g
```

The blank was already present in the original iYali26 input. Earlier external-NDH2 curation deliberately left R1889 unresolved; it did not erase a valid rule. R1889 represents mitochondrial proton-pumping NADH oxidation coupled to mature Q9 reduction, not the proposed CoQ precursor-reduction process.

## Genes and evidence

| W29 ID | Established reference name and protein function | Evidence used |
|---|---|---|
| YALI1B26679g | NUBM, 51 kDa FMN-containing NADH-oxidation subunit | Y. lipolytica deletion lacks complex-I activity; AOW01975.1 and reference CAB65520.1 are identical over all 488 residues. |
| YALI1D07089g | NUAM, 75 kDa core iron-sulfur electron-transfer subunit | Deletion/empty-vector control lacks assembled complex I, with WT complementation evidence; AOW03621.1 = CAB65519.1 over all 728 residues. |
| YALI1F09003g | NUKM, PSST core subunit in the quinone-reaction region | Deletion lacks fully assembled complex I; AOW06733.1 = CAB65525.1 over all 210 residues. |
| YALI1F22993g | NUCM, 49 kDa core subunit in the quinone-reaction region | Deletion lacks fully assembled complex I; AOW07306.1 = CAB65521.1 over all 466 residues. |

Functional evidence comes from [Waletko et al. 2005](https://doi.org/10.1074/jbc.M411488200) and [Maclean et al. 2018](https://doi.org/10.1093/hmg/ddy247). Exact sequence correspondence and database locus links support transfer to these W29 identifiers. This is **not four new W29 knockout experiments**. For NUCM and NUKM, the selected direct endpoint is complex assembly, not an independently measured native Q9 pump flux of exactly zero. Using these necessary dependencies as a Boolean AND is an explicitly bounded model inference.

The rule is not a complete subunit inventory or a four-protein minimal enzyme. Omitted subunits are not declared dispensable. NUHM has direct dependency evidence, but its W29 identity could not be established: the old crosswalk to YALI1D00766g conflicts with the protein sequence. That locus's NUHM assignment remains unverified. NUGM mapping also remains unresolved. A single bounded two-query BLAST against 7,894 cached W29 proteins returned no hits; this does not establish biological absence. See the [identity limitation record](identity/NUHM_NUGM_MAPPING.md).

YALI1F32476g—NDH2, the experimentally characterized external alternative NADH dehydrogenase—is not included as a complex-I subunit or interchangeable OR member. The tentative precursor-reduction family candidates are not included either.

## Implementation and validation

The curation records the four identities, source evidence and limitations. Before changing R1889, the shared function checks its existing GPR, reaction name, stoichiometry, bounds, metabolite formulas/charges/compartments, all four model gene identities and conflicting notes. It changes only R1889's GPR and evidence notes. Historical blank-GPR wording in the external-NDH2 helper was updated to refer to the separate complex-I curation.

- **Seven focused tests passed**, including conflict rejection, idempotence, saved-model reload and the existing quinone/R539 checks.
- In the corrected E5 candidate, each of the four gene knockouts closes R1889 to `[0,0]`; the external-NDH2 reaction remains functional. These are Boolean/bound checks, not growth simulations.
- The complete XML structure outside R1889 is unchanged, ignoring XML formatting and element order. Model/solver definitions and unrelated annotations are preserved. The source file is unchanged.
- The existing R1889 chemistry, including its previously noted H/charge residual, is retained. R2062's inherited GPR and non-pumping route are retained. No growth-essentiality or pathway-closure claim follows from this correction.
- Independent source/identity review supports **9 of 9 limited claims**, including one model inference; this is not a comprehensive complex-I audit. Figure-image access limitations are recorded in the [source audit](audit/AUDIT.md).

The [verification record](verification.json) identifies the successful targeted export. It applies the same shared curation function to the fixed historical E5 and passes unchanged-scope and roundtrip checks; it is not a claim that the entire historical environment was rebuilt.

The [independent implementation audit](implementation_audit/AUDIT.md) confirmed the complete XML comparison. Because the targeted export preserves other reactions exactly, R570 retains its inherited historical note saying its earlier patch left R1889 blank. That note describes the earlier patch's scope; R1889's actual GPR and new evidence notes give its current state. Future shared builds use the updated R570 wording.

## Full rebuild limitation

Two offline/no-solve source-build attempts reached and applied the R1889 rule, but the strict export/reload check rejected their outputs. A diagnostic capture isolated differences in **three acyl-CoA pool reactions** caused by floating-point serialization, and **20 tRNA biomass residue metabolites** whose unspecified charge became zero on serialization. These differences do not involve R1889 and are outside this GPR correction. No tolerance was relaxed, no charges were filled by this task, and no lipid work was taken over.

`E5_r1889_pipeline.xml` and `E5_r1889_diagnostic.xml` are retained as **failed-build artifacts**. Use `E5_r1889_candidate.xml`, which passed the targeted export checks. Full fresh-build acceptance remains unresolved; the exact differences are retained in `export_roundtrip_conflict.json`.

Source identities, code/data fingerprints, initial dirty state, raw sequence records, test logs and build failures are recorded in the existing delivery manifest. No optimization, AlphaFold job, baseline replacement, commit or push was performed.
