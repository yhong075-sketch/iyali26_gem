# Four dipeptide reactions: provenance and Gly-Pro identity review

Verified 2026-09-24. This branch used only files resolving inside the current workspace and public literature/database pages. No model, GPR, bounds, medium or biomass change; no LP, BLAST, new AlphaFold prediction or cluster job. The fixed baseline SHA256 is `aad701126d12d113816fda4b872333b614b4469ee1b6d8ab8c419231c89e965f`. Raw local XML fields and source hashes are in `local_xml_facts.json`; the runnable check is `check_provenance.py`.

## Model chemistry and origin

All coefficients below are one, all species are in C_va, bounds are [0,1000], and all four GPRs are empty.

| Reaction | Dipeptide ID | Stored equation | Identity restriction |
|---|---|---|---|
| R2021 | m1871[C_va] | Gly-Asp + H2O → glycine + L-aspartate | Gly is N-terminal; Asp is C-terminal. Dipeptide stereochemistry is not encoded by structure identifiers. |
| R2029 | m1862[C_va] | Gly_Glu + H2O → glycine + L-glutamate | Gly-Glu, not Glu-Gly or a gamma-glutamyl linkage. Name supports intended order; no structural identifier proves it. |
| R2034 | m1878[C_va] | Ala-Gly + H2O → L-alanine + glycine | Ala is N-terminal; this is not Gly-Ala. L-Ala is implied by the product, not independently encoded for the substrate. |
| R2039 | m1866[C_va] | Gly-Pro + H2O → glycine + L-proline | Gly is N-terminal; this is not Pro-Gly. The L-Pro product implies the intended stereoisomer, but substrate stereochemistry/cis-trans state is not explicit. |

The four reaction names are all `cytosol nonspecific dipeptidase`, despite C_va chemistry. This name/location mismatch, chemistry, directions and empty GPRs are already present in both `data/iyli21.xml` and `data/iyali26.xml`. The four reactions are absent from the final metadata field-selection table. Therefore these features were inherited, not introduced by the recent metadata or R1159 changes. The inspected input reaction entries contain no biological reference, EC annotation or author explanation for their assigned vacuolar location. A `model(1).xml` file was not found in the workspace file inventory; no external copy was sought.

The inspected local Yeast-GEM has a chemistry-matching Gly-Glu hydrolysis, r_4439, in its vacuole, with an empty GPR and confidence 1. Its notes attribute addition to Biolog PR #149 and Rhea:36463; it names EC3.4.13.18. This is a reference-model provenance lead, not proof that iYali copied that reaction, that a Yarrowia enzyme has this specificity, or that it is vacuolar. The local iYali reference did not yield exact named substrate counterparts in the bounded name inventory; this is not an exhaustive structure-based absence proof.

## Chemical fields and cross-compartment identity

Each of the four substrate names occurs as a separate extracellular, cytosolic and vacuolar pool. Those 12 pools have charge 0 but no formula or structure identifier. The vacuolar glycine product m1863 also lacks a formula and is named `L-glycine`; glycine itself is achiral. Transport connects the intended same named chemical between compartments, but name equality alone cannot validate stereochemistry.

The current product formulas/charges are H2O/0, C4H7NO4/0 (Asp), C5H9NO4/0 (Glu), C3H7NO2/0 (Ala), C5H9NO2/0 (Pro). Under this neutral-form convention, the *conditional arithmetic* formulas are Gly C2H5NO2, Gly-Asp C6H10N2O5, Gly-Glu C7H12N2O5, Ala-Gly C5H10N2O3, Gly-Pro C7H12N2O3, each charge 0. These match existing unaccepted fill candidates, but are not identity validation or applied patches.

Older `model.xml`/metadata store the same neutral Asp/Glu formulas with charge −1; local Yeast-GEM instead uses deprotonated formulas C4H6NO4/C5H8NO4 and acidic-dipeptide charge −1. Copying only a charge or formula from those versions would mix conventions. `scripts/gem_annotate/patches.py` explicitly excludes dipeptides from neutral formula filling. Formula/charge filling requires an agreed microspecies convention, structure identity and all connected pool/reaction checks; none were silently changed here.

## Gly-Pro enzyme evidence

The mechanistically appropriate class is [EC3.4.13.9, Xaa-Pro dipeptidase/prolidase](https://iubmb.qmul.ac.uk/enzyme/EC3/4/13/9.html). [EC3.4.13.18](https://iubmb.qmul.ac.uk/enzyme/EC3/4/13/18.html) includes Pro-X specificity and its broad specificity varies with source; its name cannot establish Gly-Pro activity. [EC3.4.14.5/DPP-IV](https://iubmb.qmul.ac.uk/enzyme/EC3/4/14/5.html) releases an N-terminal dipeptide from a longer peptide. Cleaving Gly-Pro-pNA after Pro is not hydrolyzing the Gly–Pro bond of free Gly-Pro.

**Priority candidate: YALI1E16433g, native established name unverified; predicted M24B prolidase-like peptidase.** Current versioned [NCBI AOW05368.1](https://www.ncbi.nlm.nih.gov/protein/AOW05368.1) and [XP_503902.2](https://www.ncbi.nlm.nih.gov/protein/XP_503902.2) contain the same 454 aa W29 sequence. Both carry the automated CDD Prolidase cd01087 region 159–427 and an AMP_N domain. The CDS comparison to S. cerevisiae YFR006W is annotation-based homology, not a new alignment or measured activity. Exact Gly-Pro hydrolysis and native localization remain untested in the evidence obtained here.

CLIB122 YALI0E13464g/CAG79495.1/Q6C610 (native name unverified; homologous prolidase-like candidate) is 454 aa but differs at residue 17: W29 Q versus CLIB122 K. Sequence-level hashes and raw GenPept are saved; the old-strain sequence/structure must not be called identical to W29. A0A1D8NIA1 is currently marked inactive for not belonging to a reference proteome; that does not negate the versioned NCBI gene/protein record.

Two additional cached M24-family candidates are YALI1E22921g/A0A1D8NJ32 (651 aa; native name unverified, metalloaminopeptidase candidate; cytoplasm IEA annotation) and YALI1B26715g/A0A1H6PH86 (488 aa; native name unverified, metalloaminopeptidase candidate; mitochondrion/nucleus IEA annotations). They are comparators to exclude, not accepted Gly-Pro isozymes. E16433 has no localization annotation in the inspected cached entry. No inspected evidence establishes any of these as a native vacuolar Gly-Pro hydrolase.

A [2021 primary genomic survey](https://pmc.ncbi.nlm.nih.gov/articles/PMC8469457/) reports two Yarrowia XPD homologs and three CND homologs, based on sequence searches/clustering. Its actual enzyme experiments concern other fungal species and DPP4 activity, not Yarrowia Gly-Pro hydrolysis. The main text does not identify the two Yarrowia XPD accessions; the panel above must not be equated to those counts.

A [2023 A. phoenicis recombinant-enzyme study](https://www.mdpi.com/2311-5637/9/11/978) directly measured free Gly-Pro cleavage by an XPD, with substrate specificity different from longer peptides. It supports the catalytic class and a suitable positive control, not W29 localization or catalytic rates. Its alkaline optimum cannot be assigned to Yarrowia.

The [2007 Yarrowia intracellular protease study](https://pubmed.ncbi.nlm.nih.gov/17227470/) measured activities in soluble extracts and characterized an aminopeptidase using Lys-pNA. It neither assigns these four exact free dipeptides to genes nor distinguishes cytosolic from vacuolar pools. [1989](https://pubmed.ncbi.nlm.nih.gov/2649495/) and [1997](https://pubmed.ncbi.nlm.nih.gov/9353927/) precursor-processing experiments support dipeptide release in secretory protein maturation; they do not demonstrate hydrolysis of these free dipeptides in the vacuole. Thus secreted, secretory-pathway, cytosolic and vacuolar claims must remain distinct.

## Decision and minimum evidence

For **all four reactions**, retain current rows for provenance and label localization/identity unresolved; none presently qualifies for an accepted GPR or relocation patch. Compare cytosolic, vacuolar and extracellular hydrolysis only as separate explicit hypothesis scenarios. Do not duplicate all locations as simultaneously available reactions, and do not migrate a reaction solely because its inherited name says cytosol.

For R2039, E16433 is a defensible **priority protein for assay**, stronger than generic DPP or M20 membership. Minimum confirmation: exact W29 protein hydrolysis of unmodified Gly-L-Pro into both Gly and L-Pro, with Pro-Gly, Gly-Pro-containing longer peptides and DPP chromogenic substrates as discriminating controls; measure activity across relevant pH/metal conditions and use inactive/heat-inactivated controls. Native compartment requires independent localization/fractionation with purity markers, followed by a genetic perturbation and rescue tied to substrate/product measurements. A structure prediction can strengthen family assignment but cannot replace these specificity/localization tests.

For R2021/R2029/R2034, apply the same exact free-substrate/product assay to defensible non-DPP candidates supplied by the companion R2200 identity review. Each substrate needs its own measurement; activity on Cys-Gly, Gly-Ala or a pNA substrate cannot fill its GPR. Native dipeptide abundance/source and transport access in the actual culture remain separate required evidence.

The first GenPept parser split on raw `//` and failed on a URL before producing a sequence result; corrected record terminator parsing uses anchored lines and all three records now pass. Initial direct API attempts failed DNS/tool fetching; authorized public HTTPS downloads used system curl with TLS verification, without altering local inputs. No failed request or sequence annotation was counted as functional confirmation.
