# Vacuolar closed-reaction evidence (UM2 stage 1), 2026-10-28

Read-only evidence gathering. No repository file, model, GPR or medium was changed. No LP was run.

## 0. Inputs and method

| Item | Identity |
| --- | --- |
| W29 proteome | `iyali26_gem_research/artifacts/tasks/coq_yeast_reference_20261001/native_proteome.fasta`, 7,894 seqs, SHA-256 `a0bacab595622864c5845411675bfb27881a4206f27572c7f0ae3abc8978e266` |
| S288C reviewed proteome | `.../r2193_iron_import_gpr_20261022/blast/proteomes/sce_UP000002311_reviewed.fasta`, SHA-256 `1aaa81860c8bda68e81e4290703d871c7c8bd12a928038772fbc91b7e17b342a` |
| YALI1 → accession → YALI0 | `cache/data/uniprot_UP000182444.json` (`genes.orfNames`, exact match on `YALI1_X00000g`); mapping in `blast/id_map.json`, `blast/other_accs.json` |
| BLAST | blastp 2.17.0+, e ≤ 1e-5. Forward: `blast/fwd_w29_vs_sce.tsv`, `blast/fwd_w29extra_vs_sce.tsv`. Reverse (43 yeast family members vs W29): `blast/fam_sce_vs_w29.tsv` |
| Model reactions read | `model/source/iyali26.xml` and `model/source/iyli21.xml` (immutable inputs) |
| Y. lipolytica screen calls | `iyali26_gem_research/state/essentiality/repository/cas9_cas12a_fitness.csv` (Ramesh 2023 acCRISPR, PO1f, SD-Leu; source SHA `ba1eca8f…`) |
| Source snapshots | `sources/` next to this file: PubMed abstracts fetched this session; iYli21 Supplementary Table 1 (`mmc3.xlsx`, SHA-256 prefix `3acf9d29a8eaa19f`) and Fig. 1 (`gr1.jpg`, `70e5454c537f08fe`) from the Europe PMC supplementaryFiles endpoint for PMC9136261 |

Tiers: **verified** = I fetched that URL or record this session and read the sentence. Abstract-level only unless stated. **recalled** = not opened. **inferred** = my reasoning.
RBH = reciprocal best hit. An ortholog hit is not proof of function. A model GPR is not experimental validation.

## 1. Reciprocal BLAST summary (W29 ↔ S288C)

Identity and coverage are from the forward search (W29 query). "Reverse best" is the best W29 hit when the yeast protein is the query.

| W29 gene (accession / YALI0) | Model role | Best S288C hit (e, %id, qcov) | Reverse best | Verdict |
| --- | --- | --- | --- | --- |
| YALI1F37658g (A0A1D8NQI3 / YALI0F30063g) | vac import, Asn/Gln/Ile/Leu/Tyr (R869, R878, R882, R887, R897) | AVT1 (3e-158, 45%, 95%); next AVT3 at 1e-10 | AVT1 → F37658g | **RBH AVT1**, clear single ortholog |
| YALI1E13514g (A0A1D8NHZ2 / YALI0E10637g) | vac export, Asn/Gln/Ile/Leu/Tyr (R870, R879, R883, R888, R898) | AVT3 (1e-170, 48%, 77%); AVT4 (5e-141, 38%) | AVT3 → E13514g; AVT4 → E13514g | **RBH AVT3**, single W29 co-ortholog of the AVT3/AVT4 pair |
| YALI1C04798g (A0A1D8N9I2 / YALI0C03608g) | vac export, Ile/Leu/Tyr (R883, R888, R898) | ATG22 (9e-121, 37%, 99%) | ATG22 → C04798g | **RBH ATG22** |
| YALI1D14670g (A0A1H6PY62 / YALI0D11836g) | vac export, Asp/Glu (R871, R876) | AVT6 (4e-100, 41%, 91%); AVT5 1e-95; AVT7 1e-88 | AVT6, AVT5, AVT7 → D14670g | **RBH AVT6**, single W29 member of the AVT5/6/7 clade |
| YALI1E21478g (A0A1D8NIW9 / YALI0E18095g) | vac Tyr import (R897, OR), His/Lys (R2055/6) | VBA2 (8e-50, 30%, 91%); VBA1 1e-47 | VBA1, VBA2 → **YALI1B04391g** (A0A1D8N691), not E21478g | Not RBH. VBA1/2-family paralog |
| YALI1A15624g (A0A1D8N4Y2 / YALI0A15576g) | vac polyamine import (R918, R925, R927) | TPO1 (8e-145, 44%, 82%) | TPO1 → **F10342g** (A15624g second, 2e-161) | Not RBH. TPO1-like paralog |
| YALI1F10342g (A0A1D8NMC6 / YALI0F07062g) | same | TPO1 (6e-162, 50%, 89%) | TPO1 → F10342g | **RBH TPO1** |
| YALI1B27564g (A0A1H6PW64 / YALI0B21142g) | spermine (R927) | TPO3 (0.0, 57%, 97%); TPO2 (0.0, 55%) | TPO2, TPO3 → B27564g | **RBH TPO3**, single co-ortholog of TPO2/TPO3 |
| YALI1F35690g (A0A1D8NQA6 / YALI0F28017g) | His/Lys (R2055/6) | VBA5 (7e-50, 28%, 87%) | VBA5 → F35690g | **RBH VBA5** (weak, 28% id) |
| YALI1D29523g (A0A1D8NFV2 / YALI0D22913g) | His/Lys (R2055/6) | VBA5 (5e-45, 27%, 74%) | VBA5 → F35690g | Not RBH. VBA3/5-like |
| YALI1E12482g (A0A1D8NHU2) | V-ATPase a, 12-gene AND group | VPH1 (0.0, 52%, 99%); STV1 (0.0, 51%) | STV1 → E12482g | RBH **STV1** by reverse, but forward best is VPH1; isoform identity unresolved |
| YALI1F38820g (A0A1H6PMT1) | V-ATPase a, 2-gene branch | VPH1 (0.0, 54%, 98%); STV1 46% | VPH1 → F38820g | **RBH VPH1** |
| YALI1D24836g (A0A1D8NFB7) | not in these GPRs | PTR2 (0.0, 50%, 87%) | PTR2 → D24836g | **RBH PTR2** (plasma-membrane peptide transporter candidate) |
| YALI1E33485g (A0A1H6PPD7) | not in these GPRs | AVT2 (6e-131, 44%, 95%) | AVT2 → E33485g | RBH AVT2 (AVT2 has no assigned function in Russnak 2001; not used) |
| YALI1F07655g (A0A1H6Q1C6) | not in these GPRs | YPQ3/RTC2 (1e-81, 46%) | YPQ1, RTC2 → F07655g | PQ-loop (YPQ) candidate |
| YALI1B12409g (A0A1D8N733) | not in these GPRs | VCX1 (1e-153, 63%) | VCX1 → B12409g | RBH VCX1 |
| YALI1F01981g / F00616g / E06664g | not in these GPRs | AQY1 (51%) / AQY3 and FPS1 / AQY3 and FPS1 | AQY1 → F01981g | aquaporin and aquaglyceroporin candidates |
| YALI1F23706g, YALI1E16433g | dipeptide candidates | DUG1 (0.0, 61%) / YFH6 (6e-164, 49%) | DUG1 → F23706g | Reused from the 2026-09 reports, see §6 |

None of the W29 genes has a verified native name. Every function below is "S. cerevisiae ortholog function", evidence level **homology-only** for Y. lipolytica.

## 2. Vacuolar amino-acid import by H+ antiport (R869 Asn, R878 Gln, R882 Ile, R887 Leu, R897 Tyr)

- **Gene:** YALI1F37658g. Verified name: none. Function: AVT1-like vacuolar H+/neutral amino-acid antiporter. Evidence: homology (RBH AVT1, 45% id).
- **Russnak, Konczal & McIntire 2001, JBC** ([PubMed 11274162](https://pubmed.ncbi.nlm.nih.gov/11274162/)), **verified (abstract)**:
  - AVT1 "is required for the vacuolar uptake of large neutral amino acids including tyrosine, glutamine, asparagine, isoleucine, and leucine".
  - Transport by AVT1 and AVT6 requires ATP and is abolished by nigericin, "indicating that the same pH gradient can drive amino acid transport in opposing directions".
  - These five substrates match R869, R878, R882, R887 and R897 exactly.
- **Tone et al. 2015** ([PubMed 25747199](https://pubmed.ncbi.nlm.nih.gov/25747199/)), **verified (abstract)**: uptake was ATP-dependent and inhibited by concanamycin A and a protonophore, and Avt1p is "a proton/amino acid antiporter". The direction and coupling written in the model (H+ vac→cyt, aa cyt→vac) match.
- **Sato, Ohsumi & Anraku 1984** ([PubMed 6088546](https://pubmed.ncbi.nlm.nih.gov/6088546/)), **verified (abstract)**: vacuolar vesicles take up Tyr, Gln, Asn, Ile and Leu (among 10 amino acids) by "seven independent H+/amino acid antiport systems". This is biochemistry without gene identities.
- **Stoichiometry:** 1 H+ per amino acid. No source read this session gives the H+:aa ratio, so this is **inferred** (model convention).
- **Y. lipolytica-specific evidence:** none found.
  - Screen: YALI1F37658g is non-essential in both assays (Cas9 FS −1.45, q = 1; Cas12a FS −0.47, q = 1). This is consistent with a non-essential AVT1 and is not evidence of function.
  - PubMed search "Yarrowia lipolytica amino acid transport vacuole" returned no hits.
- **R897 OR member YALI1E21478g:**
  - It is a VBA1/2-family paralog, not the RBH (§1).
  - UniProt P38358 (VBA2) curated text, **verified**: "vacuolar uptake of histidine, arginine and lysine and to a lesser extent tyrosine".
  - Shimazu 2005 ([PubMed 15572352](https://pubmed.ncbi.nlm.nih.gov/15572352/)), **verified (abstract)**: in vba1Δ "uptake of tyrosine … was little affected".
  - Its Tyr assignment is therefore weak.
- **Verdict (inferred):** the S. cerevisiae evidence for direction and mechanism is strong. The W29 ortholog call is clean. There is no Y. lipolytica-specific evidence.

## 3. Vacuolar amino-acid export by H+ symport

### 3a. R870 Asn, R879 Gln, R883 Ile, R888 Leu, R898 Tyr

- **Gene 1:** YALI1E13514g. Verified name: none. Function: AVT3/AVT4-like vacuolar neutral amino-acid exporter. Evidence: homology.
- **Gene 2:** YALI1C04798g. Verified name: none. Function: ATG22-like vacuolar amino-acid effluxer. Evidence: homology.
- **Russnak 2001**, **verified (abstract)**:
  - AVT3, AVT4 and AVT6 "are involved in amino acid efflux from the vacuole".
  - Efflux of "tyrosine and other large neutral amino acids" goes through AVT3 and AVT4.
  - UniProt P36062 (AVT3) curated text, **verified**, lists "tyrosine, glutamine, asparagine, isoleucine and leucine". That covers all five export reactions, including Asn (R870).
- **Yang et al. 2006, MBoC** ([PubMed 17021250](https://pubmed.ncbi.nlm.nih.gov/17021250/)), **verified (abstract)**: "Atg22, Avt3, and Avt4 are partially redundant vacuolar effluxers, which mediate the efflux of leucine and other amino acids". This is direct support for ATG22 on Leu (R888). For ATG22 on Ile and Tyr the evidence is only "other amino acids", so the support is **inferred**.
- **Sekito 2014** ([PubMed 25036121](https://pubmed.ncbi.nlm.nih.gov/25036121/)), **verified (abstract)**: "various neutral amino acids were excreted from vacuolar membrane vesicles in a manner dependent on either Avt4p or Avt3p". This independently confirms the export direction.
- **Coupling mechanism:**
  - The nigericin sensitivity in Russnak 2001 (verified) shows that efflux is driven by the pH gradient.
  - That H+ moves out together with the amino acid, which is what "H+ symport" means here, is **inferred**: it is the only coupling that lets the acidic-inside gradient drive efflux. No abstract I read states it as symport.
  - No coupling mechanism for ATG22 was found in what I read. **Unknown.**
- **Y. lipolytica:** both genes are non-essential in the screen (E13514g Cas9 q = 1; C04798g Cas9 q = 1, Cas12a q = 0.47). No functional study found.
- **Note on R870:** the name says "L-aspartate transport" but the stoichiometry carries L-asparagine, both in `iyali26.xml` and `iyli21.xml` (read this session). This is a naming error only.

### 3b. R871 Asp export, R876 Glu export

- **Gene:** YALI1D14670g. Verified name: none. Function: AVT6-like vacuolar acidic-amino-acid exporter. Evidence: homology. It is the only W29 member of the AVT5/6/7 clade, and all three yeast proteins hit it first.
- **Russnak 2001**, **verified (abstract)**: "AVT6 is responsible for the efflux of aspartate and glutamate". The efflux requires ATP and is nigericin-sensitive.
- **Chahomchuen 2009** ([PubMed 20118605](https://pubmed.ncbi.nlm.nih.gov/20118605/)), **verified (abstract)**:
  - GFP-Avt6 is "exclusively localized to the vacuolar membrane".
  - In avt6Δ the vacuolar Glu and Asp rose under nitrogen starvation. In rich medium the difference was insignificant.
- **Nakajo 2026, Sci Rep** ([PubMed 41786911](https://pubmed.ncbi.nlm.nih.gov/41786911/)), **verified (abstract)**: Avt6 also exports neutral amino acids. The avt3Δavt4Δavt6Δavt7Δ strain kept viability under N starvation but had reduced protein synthesis.
- **Y. lipolytica:** non-essential in the screen (Cas9 FS −2.09, q = 0.71). No functional study found.

## 4. Polyamine import into the vacuole by H+ antiport (R918 putrescine, R925 spermidine, R927 spermine)

- **Genes:**
  - YALI1F10342g: RBH TPO1.
  - YALI1A15624g: TPO1-like paralog, not RBH.
  - YALI1B27564g: RBH TPO3, co-ortholog of TPO2/TPO3.
  - None has a verified name. Function: TPO-family MFS polyamine/H+ antiporters. Evidence: homology.
- **Supporting a vacuolar location:**
  - Tomitori et al. 1999, JBC ([PubMed 9920864](https://pubmed.ncbi.nlm.nih.gov/9920864/)), **verified (abstract)**: TPO1 "is a polyamine transport protein on the vacuolar membrane". Disruption reduced spermine uptake by vacuolar vesicles, and the overexpression phenotype was overcome by bafilomycin A1.
  - Tomitori et al. 2001, Biochem J ([PubMed 11171066](https://pubmed.ncbi.nlm.nih.gov/11171066/)), **verified (abstract)**: TPO1–TPO4 are "polyamine transport proteins on the vacuolar membrane". TPO2 and TPO3 are spermine-specific, and TPO1 and TPO4 recognise spermidine and spermine. Putrescine is not mentioned in the abstract.
- **Contradicting a vacuolar location:**
  - Albertsen et al. 2003, JBC ([PubMed 12562762](https://pubmed.ncbi.nlm.nih.gov/12562762/)), **verified (abstract)**: HA-tagged Tpo1p–Tpo4p "co-localize with plasma membrane markers", and Tpo1p "is a plasma membrane-bound exporter".
  - Uemura et al. 2005, JBC ([PubMed 15637075](https://pubmed.ncbi.nlm.nih.gov/15637075/)), **verified (abstract)**: multicopy TPO1 "was located mainly on the plasma membrane, but with some localization on the vacuolar membrane". It catalyses uptake at pH 8 and "polyamine excretion at acidic pH".
  - UniProt Q07824 (TPO1) and Q06451 (TPO3), current curated text, **verified**: "Cell membrane polyamine/proton antiporter". Location: Cell membrane. TPO1 "Recognizes spermidine, spermine and putrescine".
- **Verdict (inferred):**
  - The vacuolar location of the TPO family is contested. The later localisation work and the current UniProt curation favour the plasma membrane.
  - Opening R918, R925 and R927 as vacuolar imports with these GPRs rests on the 1999/2001 papers only.
  - Putrescine (R918) vacuolar import has no direct support in any abstract I read.
  - The mapping of TPO2/TPO3 (spermine-only) to B27564g on R927 is consistent with Tomitori 2001.
- **Y. lipolytica:** all three genes are non-essential (q ≥ 0.27). PubMed "Yarrowia lipolytica polyamine" returned no transport study.

## 5. Context: open vacuolar His/Lys import (R2055, R2056)

- **GPR:** YALI1F35690g OR YALI1D29523g OR YALI1E21478g. The reactions are written as uniport, with no H+.
- **VBA5 location:**
  - F35690g is RBH to VBA5 at only 28% identity.
  - Shimazu 2012 ([PubMed 23047103](https://pubmed.ncbi.nlm.nih.gov/23047103/)), **verified (abstract)**: GFP-Vba5p "localized exclusively to the plasma membrane".
  - UniProt P36172 (VBA5), **verified**, still says "Vacuole membrane". This is a conflict between the curated database and the primary paper.
- **True VBA1 ortholog:** the RBH of the vacuolar VBA1 (Shimazu 2005, verified: VBA1–3 drive ATP-dependent His and Lys uptake into vacuolar vesicles) is **YALI1B04391g** (A0A1D8N691; VBA1 → B04391g at 4e-82, 36%). B04391g is **not** in the R2055/R2056 GPR.
  - Its screen call: Cas12a essential (FS −2.23, q = 0.0022), Cas9 non-essential (q = 1).
  - This is a GPR candidate to note, not a finding about the closed reactions.
- **Mechanism:** Sato 1984 describes His and Lys vacuolar uptake as H+ antiport (**verified**, abstract). The open uniport form of R2055/R2056 therefore does not match the S. cerevisiae mechanism. That conclusion is **inferred**.

## 6. V-ATPase R795

**Prior project work reused, not redone:**
- `model/reports/r794_r795_gpr_review_20260909/GENES.md` and `literature.md`: identities of all 14 GPR loci. All are "native name not verified", with V1/V0 subunit candidates assigned by family/homology.
- Same `literature.md`, L07: Kulakovskaya et al. 1993 ([PubMed 8385405](https://pubmed.ncbi.nlm.nih.gov/8385405/)), re-fetched **this session, verified (abstract)**. Vacuoles isolated from Y. lipolytica at log phase "showed a high activity of vacuolar H(+)-ATPase (0.9-1.1 U/mg protein) and efficiently generated chemical proton gradient and membrane potential across the tonoplast". The activities declined at stationary phase. This is the only Y. lipolytica-specific biochemical evidence I found that the vacuolar V-ATPase is functional.
- `model/reports/vatpase_fitness_review_20260914/REPORT.md`: Ramesh 2023 acCRISPR calls, plus the consensus-list membership in `r794…/GENES.md` (11 of 14 loci labelled E).

**Re-read this session from `cas9_cas12a_fitness.csv` (verified local data):**

| Gene | Cas9 FS / q | Cas12a FS / q |
| --- | --- | --- |
| YALI1E12482g (a subunit) | −5.22 / 8e-11, essential | −2.47 / 2e-4, essential |
| YALI1F13017g (c″) | −5.88 / 1e-14, essential | −3.19 / 4e-9, essential |
| YALI1A09766g (A) | −3.63 / 6e-4, essential | −1.94 / 0.027, essential |
| YALI1E32332g (B) | −5.46 / 4e-12, essential | −1.67 / 0.15, not significant |
| YALI1F38820g (VPH1-RBH a) | −1.75 / 1, non-essential | −1.52 / 0.31, non-essential |

- **Y. lipolytica requirement:** the V-ATPase core is required for fitness in PO1f on SD-Leu (glucose). This is the strongest Y. lipolytica-specific evidence I found that the pump is required. It is screen evidence, not a biochemical demonstration of vacuolar H+ pumping.
- **Which a subunit is vacuolar (new caveat):**
  - The essential a subunit is E12482g, which is the reverse-best hit of STV1 (the S. cerevisiae Golgi/endosome isoform).
  - The VPH1-RBH F38820g is non-essential.
  - In S. cerevisiae, Tarsio 2011 (verified, below) shows the vacuolar isoform is Vph1.
  - Which W29 a isoform serves the vacuole is therefore **unresolved**. Both forward-hit VPH1 first.
- **Alkaline pH / vma mutants in Y. lipolytica:** no study found.
  - Search terms: "Yarrowia lipolytica V-ATPase OR vma" gave irrelevant hits only. "Yarrowia lipolytica alkaline pH growth" found nothing on V-ATPase.
  - Zvyagilskaya 2001 ([PubMed 11547351](https://pubmed.ncbi.nlm.nih.gov/11547351/)), **verified (abstract)**: a Y. lipolytica isolate grows at pH 3.5–10.5 and switches to Na+-coupled Pi uptake at alkaline pH. It does not address the V-ATPase.
  - S. cerevisiae vma alkaline-pH sensitivity is **recalled** only (not opened this session).

**H+/ATP coupling ratio:**
- **Kettner et al. 2003, Biophys J** ([PubMed 14645064](https://pubmed.ncbi.nlm.nih.gov/14645064/)), **verified (abstract)**, patch-clamp of S. cerevisiae vacuoles. The coupling ratio "was found to depend strongly upon the pH difference … being approximately 2H+/ATP at high DeltapH (4 pH units) and increasing to >4H+/ATP for small or zero DeltapH". The authors attribute this to a proton shunt through the enzyme.
- **Mazhab-Jafari et al. 2016, Nature** ([PubMed 27776355](https://pubmed.ncbi.nlm.nih.gov/27776355/)), **verified (abstract)**: the yeast VO contains "ac8c'c″de", i.e. 10 proteolipids in the c-ring.
- Structural ratio: 10 H+ per 3 ATP, about 3.3 H+/ATP. This is **inferred**. It assumes one H+ per proteolipid Glu per turn and 3 catalytic sites; that assumption is **recalled** and was not read this session.
- In vivo ΔpH is small. Vacuolar pH is about 6.2 (Preston 1989, verified, below); cytosolic pH of about 7 is **recalled**.
- **Verdict (inferred):**
  - 2 H+/ATP is the lowest measured value, and it is the high-ΔpH limit. The structural maximum is about 3.3.
  - A ratio of 2 is defensible as a conservative, high-cost choice: it charges 0.5 ATP per H+ against about 0.3 for the structural maximum.
  - Kettner 2003's own data at physiological ΔpH would point higher. State it as a modelling choice within the 2–3.3 range, not as a measured in vivo constant.
- **Stoichiometry check:** in `model/source/iyali26.xml` (and `iyli21.xml`), R795 reads `1 m10[C_cy] + H2O + ATP → ADP + Pi + 2 H+[C_va]`, i.e. **1** cytosolic H+, not 2. With the neutral ATP/ADP/Pi formulas in that file, ATP + H2O → ADP + Pi is already H-balanced, so 1 H+ in and 2 out is unbalanced by one H. Verified by reading the XML; the balance arithmetic is inferred. The "2 H+[cyt]" form in the brief may describe a later candidate. Check which file is in use.

## 7. Gene-less H+ diffusion R1163 (H+[cyt] → H+[vac]) and water diffusion R1363

**H+ leak direction:**
- Preston, Murphy & Jones 1989, PNAS ([PubMed 2674942](https://pubmed.ncbi.nlm.nih.gov/2674942/)), **verified (abstract)**: wild-type vacuoles "were mildly acidic, pH 6.2". The vph1-1 mutant gave pH 6.9.
- Wada & Anraku 1994 ([PubMed 7721725](https://pubmed.ncbi.nlm.nih.gov/7721725/)), **verified (abstract)**: V-ATPase pumping generates "a membrane potential (interior positive)".
- Tarsio et al. 2011, JBC ([PubMed 21669878](https://pubmed.ncbi.nlm.nih.gov/21669878/)), **verified (abstract)**: "Loss of VPH1 results in vacuolar alkalinization". When pumping stops, the gradient dissipates, so net passive H+ movement is out of the vacuole. This step is **inferred** from the observation.
- Kettner 2003 (verified) describes a "proton shunt" through the V-ATPase.
- **Verdict (inferred):**
  - Both the chemical gradient (lumen more acidic) and the electrical gradient (lumen positive) favour passive H+ movement **vac → cyt**.
  - A passive, irreversible cyt → vac H+ leak would move protons against their electrochemical gradient. It is not physiologically plausible as written. In an FBA model it would also act as a free vacuolar proton source that bypasses the V-ATPase.
  - If a leak is wanted, the plausible form is vac → cyt.
  - Keeping R1163 closed is consistent with this. No source names a dedicated H+ channel in the yeast tonoplast.

**Water across the tonoplast:**
- The S. cerevisiae aquaporins are not vacuolar. UniProt P0CD91 (AQY1), **verified**: "Endoplasmic reticulum membrane; Cell membrane" and "non-functional in most laboratory strains".
- W29 has an AQY1-like protein (YALI1F01981g, 51% id) and two AQY3/FPS1-like aquaglyceroporins (YALI1F00616g, YALI1E06664g). None has known vacuolar localisation; this is **unknown**.
- Lipid bilayers are intrinsically water-permeable, and tonoplast aquaporins (TIPs) are a plant feature. Both are **established/uncited**: I did not open a source this session.
- **Verdict (inferred):** water equilibrates across the tonoplast in any real cell. A gene-less H2O exchange is biologically harmless. It should be reversible rather than cyt → vac irreversible, so that it is not a one-way water sink. No yeast tonoplast aquaporin evidence was found.

## 8. Dipeptides (Gly-Asp, Gly-Glu, Ala-Gly, Gly-Pro)

**(a) Origin: added by Guo et al. 2022 (iYli21), for Biolog carbon-source positives**

- iYli21 Supplementary Table 1, "Curation of iYli21 from iYali4" (`mmc3.xlsx`, **verified** this session), lists every reaction R2018–R2039 with source "Manual curation". That covers the exchanges, the plasma-membrane and vacuolar "transport via proton antiport" reactions and the four vacuolar "cytosol nonspecific dipeptidase" reactions. The bounds are as now: exchange [0, 1000], transports [−1000, 1000], hydrolysis [0, 1000], no GPR.
  - So these reactions were **not** inherited from iYali4.
  - The 2026-09 project audit (`dipeptide_energy_audit_20260924/provenance/REPORT.md`) found a matching vacuolar Gly-Glu hydrolysis r_4439 in Yeast-GEM, with notes citing Biolog PR #149. That makes Yeast-GEM a likely template; this is **inferred**.
- Guo et al. 2022, Comput Struct Biotechnol J ([PMC9136261](https://pmc.ncbi.nlm.nih.gov/articles/PMC9136261/)), full text via Europe PMC, **verified**:
  - Biolog PM1-2 plates tested "190 carbon sources".
  - "We further added all necessary reactions to enable the cross-membrane transport of nutrients or fill the pathway gaps."
- **Fig. 1** (`gr1.jpg`, **verified** visually this session): Glycyl-L-Aspartic acid, Glycyl-L-Glutamic acid, L-Alanyl-Glycine and Glycyl-L-Proline are each **red (utilised) in the W29 Biolog column** and **white (no growth) in both the "before" and "after curation" model columns**.
- The specific reason for adding these reactions is **inferred**: they are W29 Biolog positives. Despite the curation, iYli21 still predicted no growth on them. That is consistent with the secretion-only exchanges [0, 1000] (inferred). The paper text does not discuss dipeptides.
- **Naming and stoichiometry issues** (read from XML, verified):
  - The reactions are named "proton antiport", but H+ and dipeptide move in the same direction, which is symport.
  - Ala-Gly (R2032, R2033) has no H+ at all.
  - The vacuolar step moves H+ cyt → vac, i.e. into the acidic lumen. Running it as uptake would require pushing H+ up its gradient. This is **inferred**.

**(b) Y. lipolytica growth on these dipeptides**

- As **carbon** source: yes for strain W29, from the Biolog PM1 data in iYli21 Fig. 1 (verified as read from the figure; 24 h OmniLog, single study; replicate count not stated in the text I read).
- As **nitrogen** source: no data found. PubMed "Yarrowia lipolytica dipeptide" returned only protease/DPAP-processing papers.

**(c) Yeast di-peptide uptake and hydrolysis location**

- **Uptake by PTR2:**
  - Perry et al. 1994 ([PubMed 8264579](https://pubmed.ncbi.nlm.nih.gov/8264579/)), **verified (abstract)**: cloning of PTR2 by complementing a peptide-transport mutant.
  - Ito et al. 2013 ([PubMed 24060756](https://pubmed.ncbi.nlm.nih.gov/24060756/)), **verified (abstract)**: Ptr2p is "a major peptide transporter of Saccharomyces cerevisiae" in the "proton-coupled oligopeptide transporter (POT) family".
  - The plasma-membrane location is **recalled**: UniProt P32901 says only "Membrane".
  - W29 RBH of PTR2: YALI1D24836g (50% id; verified name none; non-essential in the screen). It is not in any dipeptide GPR.
- **Cytosolic hydrolysis:**
  - Kaur et al. 2009 ([PubMed 19346245](https://pubmed.ncbi.nlm.nih.gov/19346245/)), **verified (abstract)**: Dug1p is a Cys-Gly dipeptidase "with high specificity for Cys-Gly and no activity toward tri- or tetrapeptides". The abstract does not state its location; cytosolic is **recalled**.
  - The project's 2026-09 reports are reused, not redone (`dipeptide_energy_audit_20260924/reaction_identity_cards.md`):
    - YALI1E16433g: M24B prolidase-like candidate for Gly-Pro.
    - YALI1F23706g: DUG1-like (61% id to DUG1 this session, RBH).
    - Neither is accepted for these four reactions.
- **Vacuolar dipeptide importer:** none found in yeast.
  - My searches ("vacuolar membrane peptide transport Saccharomyces"; "dipeptide vacuole yeast hydrolysis peptide transport") returned no vacuolar dipeptide transporter. This is a bounded search, not proof of absence.
  - Closest related observation: Dubey et al. 2025 ([PubMed 40288644](https://pubmed.ncbi.nlm.nih.gov/40288644/)), **verified (abstract)**. γ-Glu-Met is degraded by both the cytosolic Dug2p/Dug3p complex and the vacuolar γ-glutamyl transpeptidase Ecm38p. This is a γ-Glu peptide, not an α-dipeptide, and no vacuolar importer is named.
- **Verdict (inferred):** the S. cerevisiae evidence supports plasma-membrane H+-coupled uptake by PTR2 followed by intracellular (likely cytosolic) hydrolysis. The vacuolar H+-symport-plus-vacuolar-hydrolysis route has no transporter or enzyme evidence. If the Biolog positives are to be honoured, the evidence-consistent form is PM uptake (candidate YALI1D24836g) plus a cytosolic dipeptidase. That would be a new curation proposal, not something these sources establish.

## 9. Overall verdicts on opening

| Reaction set | S. cerevisiae ortholog evidence | W29 ortholog | Y. lipolytica evidence | Opening support (inferred) |
| --- | --- | --- | --- | --- |
| R869/R878/R887/R897 import (AVT1) | Strong: substrates and H+ antiport verified | RBH, 45% | None | Supported by orthology; direction and coupling match |
| R870/R879/R883/R888/R898 export (AVT3/4, ATG22) | Strong for efflux direction; symport coupling inferred; ATG22 coupling unknown | RBH AVT3, RBH ATG22 | None | Supported by orthology. R870 is misnamed |
| R871/R876 Asp/Glu export (AVT6) | Strong: verified | RBH, 41% | None | Supported by orthology |
| R918/R925/R927 polyamine vac import (TPO) | Contested location; current curation says plasma membrane; putrescine unsupported | RBH TPO1 (F10342g), TPO3 (B27564g); A15624g paralog | None | Weak. Keep closed or flag |
| R795 V-ATPase | Coupling 2–3.3 H+/ATP; 2 is the measured high-ΔpH limit | 14-locus complex (prior review) | Vacuolar H+-ATPase activity (Kulakovskaya 1993); core subunits essential in the PO1f screen | Strong that the pump exists and is required. Check the H count in the source XML |
| R1163 H+ leak cyt→vac | Direction contradicts the gradient | — | — | Not supported as written; vac→cyt would be plausible |
| R1363 H2O | No yeast tonoplast aquaporin found; bilayer permeability established/uncited | — | — | Harmless if reversible |
| Dipeptide vacuolar route | No vacuolar importer; PTR2 at the plasma membrane | PTR2 RBH D24836g | W29 Biolog C-source positive (iYli21 Fig. 1) | The vacuolar route is unsupported; the PM + cytosol route is the evidence-consistent alternative |

## 10. Unresolved, not pursued (budget)

- Full text of Russnak 2001: exact Asn data for AVT3/4 and the explicit coupling statement. Abstract plus UniProt only.
- Which W29 a subunit is vacuolar (E12482g, essential and STV1-reverse-best, versus F38820g, VPH1-RBH and non-essential).
- Y. lipolytica vma/alkaline-pH phenotypes: none found.
- Biolog replicate count and kinetics for the dipeptide wells: Supplementary Fig. 1 curves (`mmc1.pdf`) were not inspected beyond the well labels.
- PTR2 plasma-membrane location and DUG1 cytosolic location: recalled, not opened.
