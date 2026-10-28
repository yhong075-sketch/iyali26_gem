# UM7 stage 1 — evidence: nuclear PPi, ADP-ribose, homocitrate

Gathered 2026-10-28, read-only; no repository file edited.

Sources and how they were read:
- UniProt: REST `.txt`/JSON via curl.
- PubMed: abstracts via NCBI E-utilities `efetch` (identical to the PubMed abstract); the
  PubMed page URL is given for clicking.
- SGD: GO cellular-component rows from `https://www.yeastgenome.org/backend/locus/<SGDID>/go_details`.
- Europe PMC: full-text XML for Timney 2016 and Fazius 2012.

Tiers:
- **verified**: fetched this task, and the claim is read in that text.
- **recalled**: not opened, or the fetch failed.
- **inferred**: reasoning only.

Accession caveat:
- Every YALI1-era accession in the brief (A0A1H6Q8U0, A0A1D8NCV3, A0A1D8NN06, A0A1D8N5L3,
  A0A1D8NM68) is **inactive**: "DELETED — Not part of a reference proteome".
  **verified**, e.g. https://rest.uniprot.org/uniprotkb/A0A1D8NCV3?format=json
- I used the active CLIB122 entries found by YALI0 locus search instead. The YALI1↔YALI0
  mapping was taken from the brief and not re-checked.

## Q1 Homocitrate synthase

**Y. lipolytica**

- Q12726 (HOSM_YARLI, reviewed), **verified**: https://rest.uniprot.org/uniprotkb/Q12726.txt
  - Name "Homocitrate synthase, mitochondrial", EC 2.3.3.14; gene LYS1 = YALI0F31075g.
  - The SUBCELLULAR LOCATION field reads "Mitochondrion", with no evidence code on the CC
    line.
  - FT TRANSIT 1..? "Mitochondrion" carries **ECO:0000255** (sequence-analysis prediction).
  - GO mitochondrion is IEA (SubCell). Function is ECO:0000250 (by similarity). So the mitochondrial label is a **prediction**.
- Perez-Campo 1996, **verified** abstract: https://pubmed.ncbi.nlm.nih.gov/8948100/
  - An N-terminal mitochondrial targeting signal of 41–45 aa is only *predicted*. LYS1 is a single gene; its deletion gives Lys⁻.
  - 84% identity to S. cerevisiae Lys20. No localization experiment is reported.
- **No experimental Y. lipolytica Lys1 localization found.** PubMed "Yarrowia lipolytica
  homocitrate" returned only 8948100 and a Penicillium paper.
- Penicillium chrysogenum (comparison only), **verified**:
  https://pubmed.ncbi.nlm.nih.gov/11810244/
  - A functional HS-GFP was "mainly located in the cytoplasm". The authors contrast this with nuclear HCS in S. cerevisiae.

**S. cerevisiae**

- Chen 1997, **verified**: https://pubmed.ncbi.nlm.nih.gov/9099739/
  - Immunofluorescence plus cell and nuclear fractionation put Lys20 (YDL182W) and Lys21
    (YDL131W) in the nucleus.
  - Part of the protein is extractable only with DNase and high salt. Earlier reports had placed HCS in mitochondria and cytoplasm.
- Feller 1999, **verified**: https://pubmed.ncbi.nlm.nih.gov/10103047/
  - This paper does **not** show cytosolic or mitochondrial pools; it says both isoenzymes
    "are located in the nucleus".
  - The two isoforms differ in lysine-feedback sensitivity.
- Scott & Pillus 2010, **verified**: https://pubmed.ncbi.nlm.nih.gov/20810648/
  - Lys20 acts in DNA-damage repair through Esa1 and H2A.Z. This needs nuclear localization but is "independent of HCS catalytic activity".
- Torres-Machorro 2015, **verified**: https://pubmed.ncbi.nlm.nih.gov/25628362/
  - Lys20 is recruited to double-strand breaks. Repair-active, lysine-dead LYS20 mutants exist.
- SGD GO, **verified**: https://www.yeastgenome.org/backend/locus/S000002341/go_details
  (LYS20) and .../S000002289/go_details (LYS21)
  (pages https://www.yeastgenome.org/locus/YDL182W, https://www.yeastgenome.org/locus/YDL131W)
  - Lys20:
    - nucleus IDA (Chen 1997) and nucleus HDA (Kumar 2002); mitochondrion HDA (Sickmann 2003 and Reinders 2006 mitochondrial proteomes); cytosol RCA.
  - Lys21:
    - nucleus IDA and nucleus HDA; mitochondrion IEA.
- UniProt, **verified**: https://rest.uniprot.org/uniprotkb/P48570.txt and
  https://rest.uniprot.org/uniprotkb/Q12122.txt
  - P48570 (Lys20): "Cytoplasm", no evidence code. Q12122 (Lys21): "Mitochondrion {ECO:0000305}". Neither entry reflects the nuclear IDA data.
- Consumers of homocitrate:
  - Aco2 (P39533): (R)-homocitrate → cis-homoaconitate, ECO:0000269 (Fazius 2012).
    Location: Mitochondrion, ECO:0000269 (Huh 2003; Sickmann 2003).
    **verified**: https://rest.uniprot.org/uniprotkb/P39533.txt
  - Lys4 (P49367): homoaconitate ⇌ homoisocitrate. Location: Mitochondrion, ECO:0000269.
    **verified**: https://rest.uniprot.org/uniprotkb/P49367.txt
  - Fazius 2012, **verified**: https://pubmed.ncbi.nlm.nih.gov/23106124/
    - Homoaconitases do not use homocitrate. Yeast Aco2p "specifically and exclusively contributes to lysine biosynthesis".
- Homocitrate transporter:
  - **None found.** PubMed "homocitrate transport mitochondria" returned only 20810648, and
    "homocitrate mitochondrial carrier" returned 0 hits.
  - This is absence of evidence, not proof that none exists.

**Bottom line:**
- Y. lipolytica Lys1 has only a predicted mitochondrial location (ECO:0000255).
- S. cerevisiae Lys20/21 are experimentally nuclear (IDA), with Lys20 also in mitochondrial
  proteomes (HDA).
- Homocitrate is consumed in mitochondria (Aco2, Lys4), and no homocitrate carrier is known.

## Q2 Sirtuins ("NAD nucleosidase nuclear")

**Y. lipolytica** (unreviewed; all location annotation is computational)

- YALI0D20724g = Q6C8C7, **verified**: https://rest.uniprot.org/uniprotkb/Q6C8C7.txt
  - Sirtuin class I; domain 1..296. PANTHER PTHR11085:SF10 "sirtuin-5, mitochondrial-related". GO cytosol and mitochondrion are IBA.
  - KEGG yli:2910151 "sirtuin-5", K11415. **verified**: https://rest.kegg.jp/get/yli:2910151
- YALI0D02145g = Q6CAJ8, **verified**: https://rest.uniprot.org/uniprotkb/Q6CAJ8.txt
  - Sirtuin class I; domain 7..339. PANTHER SF8 "HST3". No cellular-component GO term.
  - KEGG yli:2910515 "HST3", K11121. **verified**: https://rest.kegg.jp/get/yli:2910515
- YALI0F11583g = Q6C219, **verified**: https://rest.uniprot.org/uniprotkb/Q6C219.txt
  - Sirtuin class I; domain 43..309. PANTHER SF9 "sirtuin-1" (Sir2/Hst1 group). GO nucleus is IEA (EnsemblFungi).
  - KEGG yli:2908621 "SIR2". **verified**: https://rest.kegg.jp/get/yli:2908621
- None of the three has a UniProt CC location field.
- **No Y. lipolytica sirtuin literature found.**
  - PubMed "Yarrowia lipolytica sirtuin / SIR2 / NAD-dependent deacetylase" returned 0 hits.
  - The only HDAC hit is about class I/II inhibitors and dimorphism.
    **verified**: https://pubmed.ncbi.nlm.nih.gov/40306381/
- Inferred: D20724 (Sirt5-like) is probably not nuclear. F11583 (Sir2-like) and D02145
  (Hst3-like) are the plausible nuclear ones.

**S. cerevisiae** (UniProt, **verified**: https://rest.uniprot.org/uniprotkb/P06700.txt,
P53685, P53686, P53687, P53688 at the same pattern)

- Sir2: nucleus and nucleolus, ECO:0000269.
- Hst1: nucleus, ECO:0000269 (Huh 2003). The note says it is not nucleolar and is cytoplasmic.
- Hst2: largely cytoplasmic through Crm1-dependent export, though it shuttles.
  **verified**: https://pubmed.ncbi.nlm.nih.gov/17110954/
- Hst3: "Cytoplasm. Nucleus." (Huh 2003 large-scale).
- Hst4: nucleus, ECO:0000269.
- All five carry Rhea:43636: acetyl-Lys + NAD+ + H2O → **2''-O-acetyl-ADP-ribose** +
  nicotinamide + Lys. The product is not free ADP-ribose.

Reaction product:
- Tanner 2000, **verified**: https://pubmed.ncbi.nlm.nih.gov/11106374/
  - Deacetylation is "tightly coupled" to forming O-acetyl-ADP-ribose, at 1:1:1 stoichiometry.
  - ADP-ribosyltransfer is a low-efficiency side reaction from uncoupling.
- Landry, Slama & Sternglanz 2000, **verified**: https://pubmed.ncbi.nlm.nih.gov/11095969/
  - They first reported ADP-ribose, nicotinamide and peptide at 1:1:1 by HPLC. Tanner 2000 identified the true product as OAADPr.
- Landry 2000 PNAS, **verified**: https://pubmed.ncbi.nlm.nih.gov/10811920/
  - NAD–nicotinamide exchange needs acetyl-lysine. Deacetylation "absolutely requires NAD".
- Tong & Denu 2010 review, **verified**: https://pubmed.ncbi.nlm.nih.gov/20176146/
  - OAADPr is a 2'/3' isomer mixture. It is removed by Ysa1, by ARH3-type esterases (→ acetate + ADPr) and by a nuclear
    acetyltransferase-like activity.
- No source found for an acetyl-substrate-independent NAD-glycohydrolase activity of yeast
  sirtuins (searched; 0 hits).
- Inferred: the model's "NAD + H2O → ADPR + Nam" lumps deacetylation with OAADPr
  deacetylation.

**Bottom line:**
- Y. lipolytica has three sirtuin genes known only by annotation.
- In S. cerevisiae, Sir2, Hst1 and Hst4 are nuclear and Hst2 is mostly cytoplasmic.
- The direct product is O-acetyl-ADP-ribose, not free ADP-ribose.

## Q3 Nuclear NAD synthesis

**S. cerevisiae**

- Nma1 (Q06178): cytoplasm and nucleus, ECO:0000269 (Huh 2003 GFP).
  **verified**: https://rest.uniprot.org/uniprotkb/Q06178.txt
  - SGD: nucleus HDA and cytoplasm HDA.
    **verified**: https://www.yeastgenome.org/backend/locus/S000004320/go_details
- Nma2 (P53204): nucleus, ECO:0000269. **verified**: https://rest.uniprot.org/uniprotkb/P53204.txt
  - Anderson 2002, **verified**: https://pubmed.ncbi.nlm.nih.gov/11884393/
    - Npt1 and Nma2 are "both concentrated in the nucleus". The authors conclude NAD+ is substantially regenerated there.
  - Emanuelli 2003 (https://pubmed.ncbi.nlm.nih.gov/12597897/, **verified**): kinetics only,
    no localization.
- Qns1 (P38795): UniProt has no location field.
  - SGD: nucleus IDA and cytoplasm IDA (Suda 2003).
    **verified**: https://www.yeastgenome.org/backend/locus/S000001116/go_details
  - Suda 2003, **verified**: https://pubmed.ncbi.nlm.nih.gov/12898714/
    - GFP-Qns1 is diffuse in nucleus and cytosol. Its activity is essential.
- Pof1 (P25576): NMNAT, ECO:0000269; cytoplasm and nucleus, ECO:0000269.
  **verified**: https://rest.uniprot.org/uniprotkb/P25576.txt

**Y. lipolytica** (unreviewed; computational)

- YALI0E25652g = Q6C4M1, **verified**: https://rest.uniprot.org/uniprotkb/Q6C4M1.txt
  - Eukaryotic NMNAT family (PANTHER PTHR12039, the Nma1/2 family). CC "Nucleus {ECO:0000256|ARBA}", a rule-based prediction.
  - KEGG yli:2912725, K06210. **verified**: https://rest.kegg.jp/get/yli:2912725
- YALI0A19624g = Q6CGG2, **verified**: https://rest.uniprot.org/uniprotkb/Q6CGG2.txt
  - NMNAT activity (ARBA); GO cytoplasm and nucleus are IBA. Its PANTHER family is **PTHR31285, the same as Pof1 (P25576)**, not Nma1/2.
  - KEGG yli:2906512 "Putative NMNAT", K19785. **verified**: https://rest.kegg.jp/get/yli:2906512
  - Inferred: this is a Pof1 ortholog, not a second Nma.
- YALI0A20108g = Q6CGE0, **verified**: https://rest.uniprot.org/uniprotkb/Q6CGE0.txt
  - "Glutamine-dependent NAD(+) synthetase", EC 6.3.5.1. GO cytoplasm IBA; nucleus IEA. No CC location.
- **No experimental Y. lipolytica localization found.**

**Bottom line:**
- In S. cerevisiae, the NMNATs and Qns1 are experimentally nuclear as well as cytosolic, so
  nuclear NAD+ synthesis releasing PPi is supported.
- In Y. lipolytica, nuclear location is only predicted.

## Q4 ADP-ribose disposal (Ysa1)

**S. cerevisiae**

- Ysa1 (Q01976): ADP-ribose + H2O → ribose-5-P + AMP + 2 H+, EC 3.6.1.13, ECO:0000269.
  UniProt has no location field. **verified**: https://rest.uniprot.org/uniprotkb/Q01976.txt
  - Dunn 1999 (https://pubmed.ncbi.nlm.nih.gov/10542272/, **verified**): ADPr
    pyrophosphatase Nudix subfamily.
- SGD, **verified**: https://www.yeastgenome.org/backend/locus/S000000315/go_details
  - mitochondrion IDA (Tong 2009); nucleus HDA and cytoplasm HDA (Huh 2003); mitochondrion HDA (Sickmann; Reinders); nucleus and cytosol IBA.
- Tong 2009, **verified**: https://pubmed.ncbi.nlm.nih.gov/19251690/
  - Ysa1 is a "major regulator" of ADPr and OAADPr. It cleaves both into AMP plus (acetyl-)ribose-P. In Δysa1, ADPr and OAADPr rise by about 50%.
- Rafty 2002, **verified**: https://pubmed.ncbi.nlm.nih.gov/12370179/
  - YSA1 hydrolyses OAADPr with ADPr-like kinetics. Extracts show a cytoplasmic esterase (OAADPr → acetate + ADPr) and a distinct nuclear
    OAADPr acetyltransferase-like activity.
- Lee 2008, **verified**: https://pubmed.ncbi.nlm.nih.gov/18812159/
  - OAADPr is undetectable without all five sirtuins. Δysa1 raises OAADPr.

**Y. lipolytica**

- YALI0F26125g = Q6C0B7, **verified**: https://rest.uniprot.org/uniprotkb/Q6C0B7.txt
  - Nudix domain; "ADP-sugar pyrophosphatase". GO: ADP-ribose diphosphatase (IBA); cytosol and nucleus (IBA).
- KEGG yli:2908113 "ADP-ribose pyrophosphatase", K01515.
  **verified**: https://rest.kegg.jp/get/yli:2908113
- **No Y. lipolytica experimental data found.**

**Bottom line:**
- The only dedicated ADPr consumer is Ysa1: nucleus and cytoplasm by HDA, mitochondrion by
  IDA in S. cerevisiae; only IBA in Y. lipolytica.
- A nuclear ADPr consumer is weakly supported: high-throughput only, and the nuclear OAADPr
  activity has no assigned enzyme.

## Q5 Passive diffusion through the nuclear pore

- Timney 2016 J Cell Biol (in vivo, S. cerevisiae).
  **verified** full text: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5057280/fullTextXML
  (article https://pmc.ncbi.nlm.nih.gov/articles/PMC5057280/; abstract
  https://pubmed.ncbi.nlm.nih.gov/27697925/)
  - The introduction defines passive diffusion as "the equilibration of ions, small
    molecules, and small- to medium-sized macromolecules" between cytoplasm and nucleoplasm.
    It does not need transport receptors or GTP.
  - The abstract says the NPC is a "soft barrier" with no firm 30–60 kDa cut-off.
- Wente & Rout 2010, **verified** abstract: https://pubmed.ncbi.nlm.nih.gov/20630994/
  - The NPC must "allow the free diffusion of water, sugars, and ions". The PMC full-text fetch (PMC2944363) failed, so the full text is unread.
- Inferred: PPi, NAD+, ADP-ribose and homocitrate are all <1 kDa, far below the ~27 kDa GFP
  that diffuses passively in Timney. Neither source names these metabolites.

**Bottom line:** A primary in vivo yeast source states that small molecules equilibrate
passively between nucleus and cytoplasm. Free nucleus↔cytosol exchange of <1 kDa
metabolites is supported.

## Q6 Acetyl-CoA synthetase in the nucleus

- S. cerevisiae: Takahashi 2006, **verified**: https://pubmed.ncbi.nlm.nih.gov/16857587/
  - Acs2 is "nucleocytosolic" and supplies acetyl-CoA for histone acetylation. Nuclear acetyl-CoA synthesis is rate-limiting.
  - Mitochondrial and nucleocytosolic pools are isolated.
- UniProt P52910: "Cytoplasm. Nucleus. Note=Predominantly nuclear".
  **verified**: https://rest.uniprot.org/uniprotkb/P52910.txt
- SGD, **verified**: https://www.yeastgenome.org/backend/locus/S000004143/go_details
  - nucleus IDA and cytosol IDA (Takahashi 2006); nucleus and nucleolus IDA (Falcón 2010); cytosol IDA (de Jong-Gubbels 1998).
- Y. lipolytica: YALI0F05962g = Q6C2Q5, "Acetyl-coenzyme A synthetase"; GO cytosol IBA only.
  **verified**: https://rest.uniprot.org/uniprotkb/Q6C2Q5.txt
  - Gatter 2016 (https://pubmed.ncbi.nlm.nih.gov/27486067/, **verified**): ACS1 is essential
    for ethanol use; no localization.
- **No Y. lipolytica localization data found.**

**Bottom line:**
- S. cerevisiae Acs2 is experimentally nuclear (IDA), a nuclear PPi source.
- Y. lipolytica ACS1 has only an IBA cytosol annotation.

## Q7 PPi in the nucleus (Ipp1)

- S. cerevisiae: Serrano-Bueno 2019, **verified**: https://pubmed.ncbi.nlm.nih.gov/30826332/
  - Fluorescent fusions, activity and immunodetection show Ipp1 is "a nucleocytoplasmic
    protein" with a functional NLS.
  - Nuclear Ipp1 is kept low by proteasomal degradation. PPi is produced and used in the nucleus. Ipp2 is mitochondrial.
- UniProt P00817: "Cytoplasm", no evidence code. **verified**: https://rest.uniprot.org/uniprotkb/P00817.txt
- SGD: nucleus HDA and cytoplasm HDA (Huh 2003).
  **verified**: https://www.yeastgenome.org/backend/locus/S000000215/go_details
- Y. lipolytica: Q6C1T4 (reviewed), IPP1 = YALI0F13541g.
  **verified**: https://rest.uniprot.org/uniprotkb/Q6C1T4.txt
  - Location "Cytoplasm {ECO:0000250}" (by similarity); GO mitochondrion IBA.
- **No Y. lipolytica experimental localization found.**

**Bottom line:**
- S. cerevisiae Ipp1 is experimentally nucleocytoplasmic.
- Y. lipolytica Ipp1 is cytoplasmic by similarity only.
- Nuclear PPi can in any case diffuse to the cytosolic Ipp1 (Q5).
