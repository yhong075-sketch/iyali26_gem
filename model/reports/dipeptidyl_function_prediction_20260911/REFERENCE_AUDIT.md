# Independent reference audit for the two W29 candidates

Date: 2026-09-11. Reviewer: `/root/dipeptidyl_gene_audit`. This bounded source audit follows `TASK.md`; it does not validate the main agent's numerical results, which are still pending. No new sequence searches, structure predictions, experimental jobs, or model/GPR/medium changes were performed by the reviewer.

Targets: **YALI1B04274g** (source form YALI1_B04274g; formal gene symbol unverified), predicted DPP-IV-like dipeptidyl aminopeptidase; and **YALI1B25603g** (YALI1_B25603g; formal symbol unverified), predicted DPP-III-like metallopeptidase. Target catalytic activity, native substrates, and localization remain prediction endpoints, not assumptions.

## Primary sources actually opened

### R1. DPP-IV substrate recognition and supporting residues

Gnoth et al. 2024, DOI [10.1371/journal.pone.0289239](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0289239), publisher HTML. Locators: Introduction active-site description; Results “Characterization with derivatives of dipeptides”; “Characterization with diprotin B as substrate”; Discussion “Role of the glutamic acids E205 and E206”.

The human reference DPP4 (GeneID 1803; P27487) has the Ser630/Asp708/His740 catalytic triad; Glu205/Glu206 anchor the substrate N terminus. This paper directly tests mutations of Glu205/Glu206, Asn710, and Arg125, using fluorogenic substrates and natural peptides. Ala-Pro-AMC and Ala-Ala-AMC demonstrate that the relevant preferred residue is the **second** residue from the N terminus. LC-MS/MS also measures cleavage products of a tripeptide. Glutamate mutations affect substrates differently: they are not universally absolute requirements for every long peptide. The article describes weaker acceptance of Ser/Gly as well as the usual Pro/Ala preference, so “only Pro/Ala can ever be processed” is too strong. Its catalytic-triad attribution cites earlier work rather than newly mutating all three residues.

### R2. Catalytically inactive DPP-IV-like controls

Bezerra et al. 2015, DOI [10.1038/srep08769](https://www.nature.com/articles/srep08769), publisher HTML. Locators: Results “Overall structure” and catalytic-site comparison/Fig. 2 text; methods and references 13–14. Experimental structure: [4WJL](https://www.rcsb.org/structure/4WJL).

Human DPP10 (GeneID 57628; Q8N608) is an inactive peptidase-like potassium-channel accessory protein. Its analogous triad is Gly651/Asp727/His759; the serine-containing motif becomes GKGYGG, while the N-terminal anchor glutamates remain. Human DPP6 (GeneID 1804; P42658) is another inactive homolog; the catalytic-serine-equivalent position is Asp. Structural resemblance therefore does not establish enzyme activity. Restoring serine reportedly fails to restore activity; those mutational experiments are cited from Kin 2001 and Chen 2006, not newly executed in this structural paper, and were not independently opened in this bounded audit. The pocket/access-channel analysis in the opened paper supplies additional structural explanations. Retaining the triad in a target argues against these specific inactivation patterns, but does not alone prove activity.

### R3. Yeast DPP-III structure and mechanism

Baral et al. 2008, DOI [10.1074/jbc.M803522200](https://doi.org/10.1074/jbc.M803522200). Opened author manuscript hosted by [DESY](https://bib-pubdb1.desy.de/record/84329/files/baral_et_al_dpp3.pdf?version=3), 16 PDF pages; independently cross-checked [3CSK](https://www.rcsb.org/structure/3CSK). Locators: zero-based PDF p2–3 “Overview”/“Zinc binding”; p5 “Catalytic mechanism”; p6 “Substrate recognition”.

3CSK is the 1.95 Å structure of the **C130S** variant of S. cerevisiae YOL057W (Q08225; DPP-III protein, no formal geneName in the reviewed record). The paper reports activity similar to wild type. Zinc coordination by His460, His465, Glu517 and water is directly structural. The HEXXGH motif and second motif containing Glu517 are appropriate alignment landmarks. Glu461 is proposed to activate water; the yeast E461Q effect quoted in this paper was then unpublished, so that number is not independent published mutational validation here. His578/Arg582 have proposed catalytic/substrate-binding roles. The peptide orientation and poly-Ala octapeptide complex were **modeled**, not experimentally observed. Broad oligopeptide activity discussed in the paper does not mean every sequence is a substrate. Native target-substrate identity remains unresolved.

No PDF image was used to obtain numerical values; these points come from explicit manuscript prose. Graphical atom positions/gel panels were not visually audited.

## Fixed reference-panel identities inspected locally

Source: fresh `sources/*_reference.json` files supplied by the main agent, independently opened by the reviewer. Reviewed status is not synonymous with experimentally characterized catalytic specificity.

| Accession | Systematic identifier / symbol | Role and evidence qualification |
|---|---|---|
| P18962 | S. cerevisiae YHR028C / DAP2 | Reviewed dipeptidyl aminopeptidase B; protein-level evidence. Vacuolar reference. Specific localization papers were not fully opened in this review. |
| P33894 | S. cerevisiae YOR219C / STE13 | Reviewed dipeptidyl aminopeptidase A; protein-level evidence. Golgi/precursor-processing reference; do not transfer its location solely by similarity. |
| Q96VT7 | A. niger dapB; systematic locus not supplied in this record | Reviewed DPP-IV candidate/reference, but **proteinExistence=3, inferred from homology**. Do not label the full panel as uniformly established by direct protein assays. |
| P27487 | Human GeneID 1803 / DPP4 | Experimentally characterized active serine dipeptidyl peptidase; R1. |
| P42658 | Human GeneID 1804 / DPP6 | Inactive DPP-IV-like channel-accessory control; R2. |
| Q8N608 | Human GeneID 57628 / DPP10 | Inactive DPP-IV-like channel-accessory control; R2. |
| Q08225 | S. cerevisiae YOL057W; formal symbol not listed | Experimentally characterized DPP-III reference; R3. “DPP III” is the protein/enzyme name, not proof of an official DPP3 gene symbol. |
| Q9NY33 | Human GeneID 10072 / DPP3 | Reviewed zinc dipeptidyl peptidase III; protein-level evidence. This review does not independently audit each human substrate assay. |
| P48147 | Human GeneID 5550 / PREP | Reviewed prolyl endopeptidase; specificity-class control. A shared hydrolase architecture is insufficient to assign N-terminal dipeptide release. |

Jalving et al. 2005, DOI 10.1007/s00438-005-1134-9, PMID15812650, was located through an indexed primary abstract. It reports increased intracellular DPP-IV activity after A. niger dapB overexpression, while the membrane domain and enzyme-maturation role are predictions. The direct PubMed and publisher opens failed in this round; this remains a **partially checked contextual lead**, not an independently opened full paper or a fresh biochemical test of Q96VT7.

## Permissible inferences for the four compounds

These are conditional functional hypotheses, to be integrated with the main agent's verified alignments and residue/structure results. The scheme is `NH2–X–Y–remaining peptide + water → X–Y + remaining peptide`; it is not synthesis by condensation of two free amino acids and not hydrolysis of a free dipeptide into amino acids.

| Target dipeptide | Inference if YALI1B04274g proves DPP-IV-like | Inference if YALI1B25603g proves DPP-III-like |
|---|---|---|
| Gly-Pro | Best match among the four to the canonical second-position Pro preference; still needs a suitable longer substrate and target-enzyme validation. | Family breadth alone does not establish cleavage of a Gly-Pro-headed peptide. |
| Ala-Gly | Second position is Gly, **not Ala**. Noncanonical/weaker support; cannot rule out activity absolutely. | Plausible screening candidate only, not confirmed product. |
| Gly-Asp | Second position is acidic Asp, poorly supported by classic Pro/Ala preference. | No specific target-enzyme substrate evidence in these references. |
| Gly-Glu | Second position is acidic Glu, poorly supported by classic Pro/Ala preference. Do not confuse with γ-Glu-containing peptides. | No specific target-enzyme substrate evidence in these references. |

No result from these reference proteins proves W29/PO1f intracellular production, flux, compartment, culture dependence, or coupling to R795. These claims require independent biological evidence. Nor does peptide/protein presence in a proteome establish that the native enzyme releases that particular sequence as a free dipeptide.

## Audit coverage and stopping point

Eight scoped reference claims: (1) panel identity/annotation status, (2) DPP-IV catalytic architecture, (3) second-position substrate preference and exceptions, (4) DPP10 inactive structural comparison, (5) DPP6/10 reactivation-mutagenesis details, (6) 3CSK zinc-site structural evidence, (7) distinction between experimental structure and modeled DPP-III substrate binding, (8) exact-four target production/R795 transfer.

**8 total | 8 inspected against available content | 6 supported as reference-level facts or stated limitations | 2 unresolved | 0 contradicted | 0 wholly unchecked.** Claim 5 is only partially supported because the original mutational papers were not independently opened. Claim 8 is unsupported by these references and remains unresolved. Additional contextual lead Jalving 2005 is explicitly outside this fully opened core-source set. This coverage does not include pending target computation or untested biological conclusions.
