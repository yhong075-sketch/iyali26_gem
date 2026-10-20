# R1964 (arabinitol dehydrogenase): GPR review

Read-only. No model, data or code was changed, nothing was built, and no LP was solved. Scope, inputs
and SHAs are in [TASK.md](TASK.md).

Evidence tiers:
- **verified:** checked this session (2026-10-19) and the output is shown or cited here;
- **recorded:** taken from a project record and not re-run;
- **provable:** follows mathematically from verified or recorded facts;
- **inferred:** reasoned, with no single source;
- **related organism:** evidence from another species. Under AGENTS.md it bears on chemistry
  plausibility only; it cannot make a patch acceptable.

## 1. Summary

- **Q1: is the evidence enough for a GPR?** Wang et al. 2020 tested the purified YALI1F03452g protein
  ("YlAraDH1") on its own (verified):
  - NADH oxidation with D-xylulose: 1.4 U/mg, pH 6.0;
  - NAD⁺ reduction with D-arabitol: 5.3 U/mg, pH 8.0.

  These are both substrate sides of R1964, with the right cofactor. The **products were assigned, not
  identified**: the assay followed NAD(P)H at 340 nm only, and HPLC was used only on culture broths
  (verified, Methods). No *Y. lipolytica* paper with product identification was found. The one
  supplementary file is cited in the text only for DNA cassette sequences (it was not opened).
  - **Verdict:** this meets "catalyses the step on its own" for the substrate side. That the products
    are D-arabinitol and D-xylulose is **inferred, moderate** (§5).
  - That is enough for a *provisional* GPR, under the same status the project used for R1025/R1026,
    which rested on weaker AlphaFold-level evidence. It is not enough for a "verified" assignment.
    Whether to accept it is the user's decision.
- **Other genes.**
  - YALI1E06616g ("YlAraDH2") has the same evidence and the same activity pattern (D-xylulose 2.2,
    D-arabitol 6.3 U/mg). It is **not a model gene**.
  - The two NADP(H)-only mannitol dehydrogenases cannot be members of an NADH reaction:
    YALI1B21134g (YlMDH1, in the model on R1893) and YALI1D24054g (YlMDH2, not a model gene).
  - YALI1E15452g (YlXDH, the R801 gene) has no D-arabitol activity, so it is excluded.
- **Q2: regiochemistry.**
  - The databases define the two readings (verified): C-2 oxidation gives D-ribulose (EC 1.1.1.250,
    KEGG R07134); C-4 oxidation gives D-xylulose (EC 1.1.1.11, KEGG R05604, Rhea 17921).
  - A Fischer-projection test of Wang's activity pattern (§5.3) finds one substrate motif that fits
    every observation except xylitol oxidation. It gives **C-4 only** (D-arabinitol ⇌ D-xylulose,
    i.e. R1964) and predicts no D-ribulose reduction.
  - Every motif within one more mismatch that still allows D-arabitol oxidation also includes C-4.
    C-2 appears only in a looser motif that would also oxidize ribitol and erythritol, which were
    inactive.
  - The same test on YlXDH recovers the textbook xylitol-dehydrogenase pattern (control).
  - Homology points the other way: 42 % identity to ribulose-forming ArDHs whose product was
    identified as D-ribulose.
  - **Verdict:** the gene's own data lean to C-4 (inferred, moderate). It is not settled: D-ribulose
    was never tested, and a related fungal enzyme reduces both pentuloses.
- **Consequence for R2190.** This does not flip the R2190 recommendation.
  - Option A (D-ribulose product) stays the non-destructive choice.
  - Option B (D-xylulose plus a merge) would remove the model's only EC 1.1.1.250 reaction, with no
    evidence that the activity is absent.
  - R2190's limit notes should record that the gene's own data lean to C-4.
  - A plus a GPR on R1964 represents both readings.
- **Q3: interactions.**
  - A GPR changes no bounds, so WT is unchanged (provable).
  - After either loop fix (R2041 forward-only, or R2190 → D-ribulose), R1964 is blocked under the
    model's export-only pentose exchanges (provable, §6.2). So the GPR has no flux effect on SD-Leu;
    it changes the representation and the KO mapping only.
- **Side finding (verified): R1964's cross-references name the wrong enantiomer.**
  - MNXR154181 / BiGG ABTD1 is "L-arabinitol 4-dehydrogenase": D-arabinitol ⇌ **L**-xylulose.
  - The model's species `m1015` is **D**-xylulose. The equation-consistent IDs are MNXR198271,
    Rhea 17921, KEGG R05604, EC 1.1.1.11 and BiGG ABTDG.
  - The colleague-owned lipid-unlump audit already lists this as a blocker (read only).
- **Recommended plan (needs approval), §7:**
  - an opt-in, guarded, post-selection step: provisional GPR `YALI1F03452g or YALI1E06616g` on R1964
    plus the annotation correction;
  - 5 LPs.

## 2. R1964 in each state (verified)

| State | Equation (stored IDs) | Bounds | GPR | Annotations | Notes |
| --- | --- | --- | --- | --- | --- |
| Source `iyali26.xml` | m1015 + **m10 (H⁺)** + m123 ⇌ m122 + m1823 | [−1000, 1000] | none | none | `PROTEIN_CLASS 1.1.1.11` |
| a56f3649 | m1015 + m123 ⇌ m122 + m1823 | [−1000, 1000] | none | bigg/vmh ABTD1, metanetx MNXR154181, SBO:0000176 | `PROTEIN_CLASS 1.1.1.11` |

Stored species in a56f3649:
- D-xylulose `m1015` is C5H10O5, charge 0. MNXM1371292, CHEBI:17140, InChIKey
  `ZAQJHHRNXZUBTE-WUJLRWPWSA-N`; this equals the D-xylulose row of local MNXref 4.5. Its other
  reactions are R801 (xylitol → D-xylulose, YALI1E15452g) and R803 (xylulokinase).
- D-arabinitol `m1823` is C5H12O5, charge 0 (MNXM1018, CHEBI:18333).
- NADH `m123` and NAD `m122` are stored neutral, as in the R2190 review.

The element and charge residual is {} and 0. The source H⁺ was dropped exactly as for R2190; with
neutral NAD/NADH the a56f3649 equation is the balanced one, so no H⁺ edit is needed.

### 2.1 Cross-reference check (verified unless marked)

| Identifier | Resolves to | Pentulose | Source |
| --- | --- | --- | --- |
| MNXR154181 (on R1964) | biggR ABTD1 "L-Arabinitol 4-Dehydrogenase": D-arabinitol + NAD ⇌ L-xylulose + NADH + H⁺ (also VMH ABTD1) | MNXM1371095 **L-xylulose**, InChIKey `…-WVZVXSGGSA-N` | local `reac_prop.tsv`, `reac_xref.tsv`, `chem_prop.tsv` |
| MNXR198271 (equation-consistent) | D-arabinitol + NAD ⇌ D-xylulose + NADH + H⁺; EC 1.1.1.11; xrefs Rhea 17921–17924, BiGG ABTDG "D-arabitol dehydrogenase", MetaCyc D-ARABINITOL-4-DEHYDROGENASE-RXN | MNXM1371292 D-xylulose (= `m1015`) | local MNXref |
| KEGG [R05604](https://www.kegg.jp/entry/R05604) | "D-arabinitol:NAD 4-oxidoreductase"; C01904 + C00003 ⇌ C00310 + C00004 + C00080; EC 1.1.1.11 | C00310 D-xylulose | KEGG REST |
| Rhea [17921](https://www.rhea-db.org/rhea/17921) | D-arabinitol + NAD(+) = D-xylulose + NADH + H(+); EC 1.1.1.11; 3 enzymes | CHEBI:17140 | Rhea REST |
| EC [1.1.1.11](https://enzyme.expasy.org/EC/1.1.1.11) | ENZYME "D-arabinitol 4-dehydrogenase"; KEGG SYSNAME "D-arabinitol:NAD+ 4-oxidoreductase". The only Swiss-Prot entries are bacterial DalD (O52720, P58708, P58709) | D-xylulose | ENZYME text, KEGG REST |

**Which side is wrong (mismatch rule, AGENTS.md / CLAUDE.md 6a).** The local ID was opened first.
- `m1015` resolves to D-xylulose: MNXM, ChEBI and InChIKey all say so, and R801 and R803 use it as
  D-xylulose (verified).
- The stoichiometry and the source note `PROTEIN_CLASS 1.1.1.11` both describe EC 1.1.1.11.
- The ABTD1 record pairs D-arabinitol with L-xylulose. One dehydrogenation of D-arabinitol gives
  D-ribulose or D-xylulose, never L-xylulose (§5.1). The record is therefore internally impossible.
- So the attached annotation is the wrong side (inferred, strong).

How MNXR154181 was attached was not traced; the BiGG/name route is inferred. The colleague-owned
`r1521_kegg_mnx_frontier_source_audit.json` row (2026-08-18) has verdict `contradicted`:
"resolve D-xylulose versus the L-xylulose species in ABTD1 before assigning an identifier". Its
`blockers` list carries "R1964: D-/L-xylulose stereochemistry conflict". The file was read, not
edited.

## 3. Candidate genes (AGENTS.md gene records)

| Systematic ID | Identity (verified) | Verified name | Protein function | Evidence level | In a56f3649 |
| --- | --- | --- | --- | --- | --- |
| **YALI1F03452g** | GeneID 2908139, XP_504895.1 (343 aa) = YALI0F02211g = Q6C367, identical (R2190 review); old tag YALI2_F00952g | none. NCBI: "uncharacterized protein"; "YlAraDH1" is Wang 2020's name; UniProt/KEGG EC 1.1.1.250 is a homology label | NAD(H) short-chain polyol dehydrogenase | direct in vitro activity of the purified protein; products assigned, not identified; no deletion phenotype on glucose | R2190 only |
| **YALI1E06616g** | GeneID [2912324](https://www.ncbi.nlm.nih.gov/gene/2912324), [XP_503594.1](https://www.ncbi.nlm.nih.gov/protein/XP_503594.1) (291 aa) = YALI0E05643g = UniProt [Q6C6W8](https://www.uniprot.org/uniprotkb/Q6C6W8), sequences identical (compared this session); old tag YALI2_E00126g. Wang's primer termini `MSNSAKAAV…` / `…IDGGYTII*` match | none. NCBI: "uncharacterized protein"; KEGG [yli:2912324](https://www.kegg.jp/entry/yli:2912324) "(RefSeq) Sorbose reductase sou1", no KO; "YlAraDH2" is Wang's name | NAD(H) polyol dehydrogenase, same substrate pattern as YlAraDH1 | as YlAraDH1 | not a model gene |
| YALI1E15452g | YALI0E12463g (S2 map, GeneID 2912835) | "YlXDH" in Wang 2020 | NAD(H) xylitol dehydrogenase; **no D-arabitol activity** | in vitro (Wang 2020) | R801, R257, R513 |
| YALI1B21134g | YALI0B16192g (S2 map) | "YlMDH1" (Wang 2020) | NADP(H)-only mannitol dehydrogenase | in vitro | R1893 |
| YALI1D24054g | YALI0D18964g (S2 map, GeneID 2911146) | "YlMDH2" (Wang 2020) | NADP(H)-only mannitol dehydrogenase; its deletion abolished mannitol and D-arabitol production on glucose | in vitro plus deletion phenotype | not a model gene |

A model or GPR assignment is not experimental validation.

## 4. Is the evidence enough for a GPR? (task 1)

### 4.1 What Wang et al. 2020 measured (verified, PMC JATS full text)

Source: Wang N *et al.* 2020, *Biotechnol Biofuels* 13:176,
[doi:10.1186/s13068-020-01815-8](https://doi.org/10.1186/s13068-020-01815-8),
[PMC7576711](https://pmc.ncbi.nlm.nih.gov/articles/PMC7576711/).

- **Protein.**
  - Each gene was cloned from strain CGMCC7326 into pET28a.
  - The protein was His-tag purified on Ni-NTA, then run on Sephadex G-25; purity was above 95 %.
- **Assay.**
  - Activity was followed only as NAD(P)H at 340 nm.
  - Polyols other than mannitol were at 50 mM with 2 mM NAD(P) in Tris pH 8.0; mannitol was at
    100 mM. Fructose was reduced at 100 mM with 2 mM NAD(P)H in acetate pH 6.0. For L-sorbose and
    D-xylulose only the pH (6.0, Table 2) is given.
  - HPLC (RI detector, Shodex SP0810) was used only for culture broths.
  - So every "product" in Table 2 is an assignment.
- **Supplementary material.** There is a single file, `13068_2020_1815_MOESM1_ESM.docx`. The main text
  cites it only for DNA cassette sequences: sequences 1, 2, 3, 8 and 9 by number, plus "the five
  genes deletion cassettes". It was not downloaded (a file download needs the user's OK).
- A text search found no empty-vector, boiled-enzyme or no-enzyme control.

Activities relevant to R1964 (U/mg; product as named by the authors):

| Substrate (pH) | YlAraDH1 NAD(H) | YlAraDH2 NAD(H) | YlMDH1 NADP(H) | YlMDH2 NADP(H) | YlXDH NAD(H) |
| --- | --- | --- | --- | --- | --- |
| D-arabitol (8.0) | D-xylulose 5.3 | D-xylulose 6.3 | D-xylulose 0.12 | D-xylulose 0.83 | — |
| xylitol (8.0) | D-xylulose 3.6 | D-xylulose 4.6 | 0.23 | 1.4 | D-xylulose 4.35 |
| D-xylulose (6.0) | arabitol 1.4 | arabitol 2.2 | arabitol 0.12 | arabitol 0.54 | xylitol 2.34 |
| mannitol (8.0) | fructose 10.3 | 6.5 | 2.7 | 8.4 | 0.12 |
| D-sorbitol (8.0) | L-sorbose 7.5 | L-sorbose 5.2 | 0.35 | 0.86 | **fructose** 2.23 |
| fructose (6.0) | mannitol 2.6 | 1.7 | 0.82 | 5.6 | 0.04 |
| L-sorbose (6.0) | sorbitol 2.4 | 2.6 | 0.15 | 0.23 | — |

- None of the five enzymes was active on ribitol, erythritol, glycerol or the primary alcohols.
- None was active on ribose, mannose, xylose, glucose or fructose-6-P.
- D-ribulose, D-arabinose and L-xylulose were not tested.
- In vivo, on 300 g/L glucose: ΔYlAraDH1 and ΔYlAraDH2 showed no change in growth or polyols. Only
  ΔYlMDH2 lost both mannitol and D-arabitol. YlAraDH1 mRNA was low relative to YlMDH2.
- The authors name products per enzyme class. For example, sorbitol gives L-sorbose for YlAraDH1 but
  fructose for YlXDH, and D-xylulose gives arabitol versus xylitol. No product analysis is described,
  so these names are their interpretation (inferred).

### 4.2 Conditions for "catalyses R1964 on its own"

| # | Condition | Result | Tier |
| --- | --- | --- | --- |
| 1 | The tested protein is the model gene's product | Primer termini match XP_504895.1 (R2190 review) and XP_503594.1 (this session). The CGMCC7326 clones were not compared internally with W29 | verified (termini); internal sequence unchecked |
| 2 | Activity without a partner subunit | purified single protein | verified |
| 3 | Cofactor matches R1964 (NADH/NAD⁺) | NAD(H)-dependent | verified |
| 4 | Both R1964 substrates are turned over | D-xylulose + NADH and D-arabitol + NAD⁺ (table above) | verified |
| 5 | Products are D-arabinitol (from D-xylulose) and D-xylulose (from D-arabitol), not xylitol / D-ribulose | assigned, not measured | **inferred, moderate (§5)** |
| 6 | Physiological role | none shown: no deletion phenotype on glucose, low expression | verified absence of evidence; not evidence against |
| 7 | Host background excluded | not described; >95 % purity | limit (inferred) |

**Verdict.** Conditions 2–4 are met directly. Condition 1 is met at the termini. Condition 5, the one
that distinguishes R1964 from R801-reverse (xylitol) and R2190-with-ribulose, is inferred.
- This supports a **provisional** assignment, labelled "in vitro activity of the purified *Y.
  lipolytica* protein; reaction products assigned, not identified". Project precedent: R1025/R1026
  carry `user_authorized_provisional_assignment` on AlphaFold-level evidence.
- It does not support a "verified" label.
- YALI1E06616g meets the same conditions at the same tier.

### 4.3 Other literature

**Y. lipolytica (searched).** PubMed "Yarrowia lipolytica arabitol" gave 13 hits, triaged by title;
the relevant abstracts were read. NCBI/KEGG/UniProt were checked for both genes.
- [Liu *et al.* 2026, *J Agric Food Chem*](https://pubmed.ncbi.nlm.nih.gov/42550668/) built a de novo
  xylitol route from glucose in *Y. lipolytica* by combinatorial screening of D-arabitol
  dehydrogenases plus an NADPH xylitol DH (abstract verified). This is the most likely place for
  product data on these enzymes. The full text returned **HTTP 403 and is unread**; which ArDHs were
  screened is unknown.
- [Yang *et al.* 2021](https://pubmed.ncbi.nlm.nih.gov/34715526/): D-arabitol fed-batch with
  proteomics; the abstract names no ArDH gene.
- [Szczepańczyk *et al.* 2023](https://pubmed.ncbi.nlm.nih.gov/37277427/): arabitol rises about 8-fold
  in erythrose-reductase knockouts; no enzyme named.
- [Ryu & Trinh 2018](https://pubmed.ncbi.nlm.nih.gov/29150499/): L-arabitol DH limits arabinose use.
  Not about D-arabitol.
- **No *Y. lipolytica* study that identifies the product of YlAraDH1 or YlAraDH2 was found.** Not
  found is not proof of absence.

**Related organisms (chemistry plausibility only; abstracts verified unless marked).**
- **Ribulose-forming SDR ArDHs, products identified:**
  - *P. stipitis* ([Hallborn 1995](https://pubmed.ncbi.nlm.nih.gov/7483848/)): substrates are
    D-arabinitol, xylitol, D-glucitol and galactitol; with D-arabinitol the product is D-ribulose;
    85.5 % identical to *C. albicans* ArDH.
  - *C. albicans* ([Wong 1993](https://pubmed.ncbi.nlm.nih.gov/8407803/)): D-arabitol ⇌ D-ribulose.
  - *C. tropicalis* ([Quong 1993](https://pubmed.ncbi.nlm.nih.gov/8250887/)): "specific for" D-arabitol,
    reduces D-ribulose. This abstract writes the oxidation as "at carbon 4", numbering from the other
    end; in IUBMB numbering it is C-2 (§5.1).
  - [Murray 1995](https://pubmed.ncbi.nlm.nih.gov/7698655/): these enzymes form a subclass of closely
    related fungal pentitol dehydrogenases.
  - YALI1F03452g shares about 42 % identity with them (recorded, R2190 review §4.2).
  - The *P. stipitis* enzyme is broad yet ribulose-forming. That weakens the R2190 review's §4.4
    argument that YlAraDH1's broad range speaks against the ribulose label.
- **Xylulose-forming (C-4) yeast/fungal activities:**
  - *Pichia anomala* ArDH ([Zhang 2014](https://pubmed.ncbi.nlm.nih.gov/24170383/)) converts arabitol
    to xylulose. Adamczyk 2023 summarizes it as NADH-dependent with reversible D-arabitol
    4-dehydrogenase activity (second-hand).
  - *Uromyces fabae* ARD1p ([Link 2005](https://pubmed.ncbi.nlm.nih.gov/15796718/)) is NADP(H). It
    oxidizes D-arabitol and mannitol and reduces D-xylulose, **D-ribulose** and D-fructose. So one
    fungal enzyme with a YlAraDH1-like substrate set accepts both pentuloses.
  - *Z. rouxii* extracts ([Ingram & Wood 1965](https://pubmed.ncbi.nlm.nih.gov/14292984/)) reduce
    D-xylulose with NADH (an XDH) and D-ribulose with NADPH, as separate enzymes.
  - *R. toruloides* ([Adamczyk 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10398940/), full text
    read in part): genetic evidence for D-xylulose → D-arabitol → D-ribulose in vivo.

## 5. Regiochemistry: C-2 (D-ribulose) or C-4 (D-xylulose)? (task 2)

### 5.1 Definitions

- **C-2:** D-arabinitol oxidized at C-2 gives D-ribulose: EC 1.1.1.250 "2-oxidoreductase
  (D-ribulose-forming)", KEGG R07134, Rhea 17389 (verified in the R2190 review).
- **C-4:** oxidized at C-4 it gives D-xylulose: EC 1.1.1.11 "D-arabinitol:NAD+ 4-oxidoreductase",
  KEGG R05604, Rhea 17921 (verified here).
- **Xylitol:** oxidation at C-2 gives D-xylulose (KEGG [R01896](https://www.kegg.jp/entry/R01896),
  "xylitol:NAD+ 2-oxidoreductase (D-xylulose-forming)", verified). At the other, enantiotopic end it
  gives L-xylulose.
- The labels "C-2" and "C-4" depend on numbering (see Quong 1993). The product names are the
  unambiguous statement.

`fischer_regio.py` (scratchpad, §9) does the bookkeeping on textbook Fischer configurations
(established, uncited). A polyol is the L/R sequence of its secondary OH groups, C-2 downwards, and a
180° rotation reverses the sequence and swaps L and R. It reproduces every definition above:
- D-arabinitol (L,R,R) gives D-ribulose at the top and D-xylulose at the bottom;
- D-xylulose plus a new C-2 OH gives D-arabinitol on the left and xylitol on the right;
- D-ribulose gives D-arabinitol or ribitol.

### 5.2 A constraint on Wang's product column (provable under one assumption)

**Assumption:** one productive binding mode, which means one hydride face. That is established for
NAD(P) dehydrogenases (uncited).

With one face, the enzyme reduces D-xylulose to *either* D-arabinitol *or* xylitol, never both. By
microscopic reversibility it then makes D-xylulose from only one of those two polyols. Wang names
D-xylulose as the product of *both* D-arabitol and xylitol oxidation. Under the assumption, at least
one of those two assignments is wrong. So the column cannot be taken at face value.

### 5.3 Which reading fits the measured pattern (inferred; script output)

Model: a carbinol is oxidized only if the configuration read from it starts with a fixed motif; the
same motif governs which face reduces a ketose. Every motif of length 1–3 was scored against the
qualitative pattern:
- oxidized: mannitol, D-glucitol, D-arabitol, xylitol;
- not oxidized: ribitol, erythritol;
- reduced: fructose, L-sorbose, D-xylulose.

| Enzyme | Motif | Mismatches | D-arabinitol → | D-xylulose + NADH → | D-ribulose reducible |
| --- | --- | --- | --- | --- | --- |
| YlAraDH1 (same pattern: YlAraDH2) | **L L R** | 1 (xylitol not predicted) | **D-xylulose only** | **D-arabinitol** | no |
| | L | 2 (ribitol, erythritol would be oxidized) | D-ribulose **and** D-xylulose | D-arabinitol | yes |
| | R L / R L R | 2 (mannitol, D-arabitol not oxidized) | none | xylitol | no |
| YlXDH (control) | R L / R L R | 1 (L-sorbose would be reduced) | none | xylitol | no |

Reading:
- **Best fit.** For YlAraDH1 the best fit is unique and gives R1964's chemistry: D-arabinitol ⇌
  D-xylulose. It also reproduces the authors' other named products: mannitol → fructose,
  D-sorbitol → L-sorbose, fructose → mannitol, L-sorbose → sorbitol.
- **C-4 in every live motif.** Every motif with at most two mismatches that keeps D-arabitol
  oxidation includes C-4. C-2 (D-ribulose) enters only with motif "L", which predicts ribitol and
  erythritol oxidation; both were inactive.
- **Control.** For YlXDH the same procedure recovers the textbook xylitol/sorbitol-dehydrogenase
  stereochemistry: xylitol → D-xylulose, D-glucitol → D-fructose, D-xylulose → xylitol, no
  D-arabitol activity. These are also the authors' named products. The opposite faces of the two
  enzyme classes on D-xylulose are consistent with the authors naming "arabitol" for the AraDHs and
  "xylitol" for XDH.
- **Exception.** Xylitol oxidation by YlAraDH1 (3.6 U/mg) fits no single motif. Under the best motif
  it would need the other end, giving L-xylulose, or a second binding mode. This shows the simple
  model does not capture the enzyme completely.

**Limits of the test.** It is a qualitative prefix model, and it uses presence or absence of activity,
not rates. It assumes one binding mode, and it takes the stated stereoisomers (D-arabitol, D-sorbitol,
L-sorbose) as given. The U. fabae enzyme shows a polyol DH can accept both pentuloses. Nothing here
is a measurement of the product.

### 5.4 Verdict

- The gene's own data lean to **C-4 → D-xylulose** (R1964) (inferred, moderate).
- Homology leans to **C-2 → D-ribulose** (EC 1.1.1.250, the KO, the identified products of the 42 %
  homologs).
- Neither is decisive. D-ribulose was never offered to the enzyme.
- **What would settle it** (wet lab, outside this project's computation): incubate purified YlAraDH1
  with D-arabitol and NAD⁺, then identify the ketose by HPLC/GC against D-ribulose and D-xylulose
  standards. Alternatively, measure NADH oxidation with D-ribulose next to D-xylulose; both
  standards are sold, as Adamczyk 2023 shows.

**For the R2190 options:**
- **Option A** (product → D-ribulose) still matches every database and deletes nothing. It keeps
  YALI1F03452g on a reaction whose regiochemistry rests on homology only.
- **Option B** (product → D-xylulose, then merge with R1964) matches the stereochemical reading.
  But it removes the model's only EC 1.1.1.250 reaction, and AGENTS.md requires evidence that the
  activity does not exist; there is none. Not recommended now. If product identification ever shows
  D-xylulose only, R1964 with the GPR proposed here *is* the merged state, and B reduces to dropping
  R2190.
- **Recommendation:** keep the R2190 plan as is. Add one limit line to it: "the gene's own
  enzymology (Wang 2020), read stereochemically, favours C-4; see r1964_gpr_review_20261019".

## 6. Interactions (task 3)

### 6.1 What a GPR on R1964 can change (provable)

- A GPR enters the LP only through gene knockouts. **WT growth is unchanged** by construction.
- **KO of YALI1F03452g** then closes R1964 as well as R2190. The feasible set shrinks, so KO growth
  can only stay equal or fall.
  - Equality is expected (inferred): all 109 recorded R1964 values in 35 report files are exactly 0
    (verified scan of `model/reports`). R2190 has the same count (R2190 review). These are single
    optima from several models, not FVA.
- **KO of YALI1E06616g** (option 1a) closes nothing on its own, because R1964 stays open through the
  other OR member.
- **Benchmark.** Neither gene is in `consensus_essential_genes.csv` (0 hits, verified). Both belong to
  the user-defined negative class, so the change cannot create a TP. A FP would need an essential KO,
  which is not expected.

### 6.2 R1964 after each loop fix (provable from connectivity; bounds verified in a56f3649)

Take the C5 subnetwork S = {D-arabinose `m454`, D-arabinitol `m1823`, D-xylulose `m1015`, xylitol
`m1014`, D-xylose `m1016`}. Every reaction touching S has coefficient 1 on its S species.

- **Inflow:** only R2041 in reverse (A5P → D-arabinose, ATP-forming). Imports through R1063, R1980 and
  R1075 cannot run inward, because their only partners R1062, R1981 and R1074 are export-only
  [0, 1000]. `sd_leu.csv` has no pentose or pentitol row.
- **Outflows**, all ≥ 0 by bounds: R803, R254, R255, R2041 forward, and the three exports.
- **With R2041 forward-only:** S has no inflow, so every outflow is 0. The `m454`, `m1823`, `m1015`
  and `m1014` balances then give R2190 = R1964 = R801 = R802 = 0.
- **With R2190 → D-ribulose (Option A):** R2190 = 0 (dead-end metabolite). S without `m454` then has
  no inflow, so R1964 = 0.
- **In a56f3649 (no fix):** R1964 carries flux only in the loop. The recorded closed-boundary range is
  [−1000, 0] (orphan inventory `loop_check.txt`).

So, unless a medium opens uptake on R1074, R1062 or R1981, R1964 is a blocked reaction once either
fix lands. The GPR then has **no flux consequence in any recorded context**. It records the gene and
fixes the KO mapping for xylose or arabitol scenarios.

### 6.3 Other lines

- **R2041 forward-only** (implemented in `claude/bold-chaplygin-1c16df`, candidate `da0c50ba…`,
  uncommitted; known from its hand-over only). It is independent of this change. It touches
  `main.py` and `cli.py`, as does the GPR-fix worktree `claude/clever-kare-7d2e1a` (`51c6fabd…`).
  A third flag there needs careful merging.
- **R2190 plan** (`claude/clever-austin-f6a04e`, uncommitted). Independent: A and the R1964 GPR
  touch different reactions. If both land, YALI1F03452g sits on both regiochemistries, each with
  its own evidence label.
- **Colleague audit.** Its R1964 and R2190 rows agree with this review. Updating its `blockers` list
  is the owner's call.

## 7. Plan for approval

### 7.1 Decisions

| # | Decision | Options | Recommendation |
| --- | --- | --- | --- |
| D1 | GPR on R1964 | (a) `YALI1F03452g or YALI1E06616g`, provisional; (b) `YALI1F03452g` only, provisional; (c) keep orphan, record the candidates | **(a)**. Both genes meet the same conditions (§4.2). Leaving YlAraDH2 out would make the YlAraDH1 KO close R1964, although the data show a second NAD(H) enzyme with equal activity. (b) is smaller (no new gene) |
| D2 | Cross-references of R1964 | replace bigg/vmh ABTD1 and MNXR154181 with MNXR198271, rhea 17921 (and 17922–17924 as MNXref lists), kegg.reaction R05604, ec-code 1.1.1.11, bigg ABTDG, metacyc D-ARABINITOL-4-DEHYDROGENASE-RXN; keep SBO | **bundle with D1.** Annotation-only; stoichiometry, mass and charge unchanged |
| D3 | Wiring | (i) opt-in flag, for example `--r1964-gpr`; (ii) always on | **(i)**, so builds without it reproduce a56f3649 byte for byte (the R2041 and R132 pattern) |

### 7.2 Steps (curation-change skill; nothing done)

1. **TASK.md** for the change, quoting the user's authorization.
2. **Data:** `model/curation/gpr/r1964_gpr_assignment.json`, laid out like
   `r1025_gpr_assignment.json`.
   - `schema_version` 1; `assignment_id` `R1964-ARADH-GPR-<date>`; `status`
     `user_authorized_provisional_assignment`.
   - `before_gpr` "", `after_gpr` per D1.
   - The exact a56f3649 bounds, stoichiometry and species (§2).
   - Before and after annotations (D2).
   - `gene_annotation` for each gene: NCBI GeneID, RefSeq, UniProt. For YALI1E06616g also the data
     to create it: `ncbigene` 2912324, `refseq` XP_503594.1, `uniprot` Q6C6W8.
   - Source URLs: Wang 2020 and the database pages above.
   - Notes:
     - `gpr_evidence_basis` = "in vitro activity of purified *Y. lipolytica* protein (Wang 2020
       Table 2); products assigned, not identified";
     - `gpr_experimental_confirmation` = `required_not_confirmed`;
     - `r1964_regiochemistry_limit` (§5.4);
     - `r1964_previous_annotation` (ABTD1/MNXR154181, L-xylulose);
     - the authorization.
3. **Code:** a guarded `apply_r1964_gpr_assignment(model, spec=None)`.
   - The existing `_apply_gpr_assignment` fetches every gene with `get_by_id`, so it cannot create
     YALI1E06616g. Use a new function or module, as `r2041_direction.py` did.
   - Fail closed unless all of these match: schema, ID, status, reaction ID, name, bounds, exact
     stoichiometry and species, and the before- or after-state GPR and annotations.
   - YALI1E06616g must be absent, or present with exactly the spec annotations.
   - No other reaction may use the new gene.
   - It sets the GPR, annotations and notes, and creates the gene with its annotations. A rerun is a
     no-op.
4. **Wiring:** after `apply_metadata_reaction_selection` (`main.py:650`), so the selection cannot
   overwrite it, and at the end of the E5 chain, after `r132_mito_adh`, so its before-state is
   a56f3649's.
   - The step runs after `annotate_reactions` (`main.py:238`) and `backfill_reaction_xrefs`
     (`main.py:317`) (verified). The gene-EC route therefore cannot re-attach YALI1F03452g's
     EC 1.1.1.250 ribulose IDs to R1964, which is how R2190 got them.
5. **Focused test** on an in-memory model:
   - applies on the a56f3649 state and is idempotent;
   - fails closed on a stale GPR, bounds, stoichiometry, annotation, a gene-ID collision or a wrong
     gene annotation;
   - flag off changes nothing.
6. **Build** a new file, for example
   `model/candidates/E5_coq9_resp_pkg_r132_r1964_gpr_<date>.xml`.
   - Command: a56f3649's recorded `build.json` argv plus the flag.
   - Control: the same build without the flag must reproduce `a56f3649…`.
7. **Whole-model diff** against a56f3649: R1964's GPR, annotations and notes, plus one new gene
   (D1a) only.
8. **Then:** `pytest` compared test by test, the LPs in §7.3, an independent read-only audit,
   REPORT.md and a STATE.md entry. No commit or push without an explicit instruction.

**Base.** This worktree on `285eaee` gives a one-change diff against a56f3649 (recommended).
Stacking on `da0c50ba` (R2041) or `51c6fabd` (GPR fixes) waits for them to land.

### 7.3 LP checks (need authorization; none run)

- **Context:** `coq9_respiratory_package_20261008/run_controls.py`, as in the R2041/R2190 plans:
  PO1f SD-Leu, Gurobi with 1 thread, Presolve 0, tolerances 1e-7, 60 s per solve.
- **Models:** M0 = a56f3649; M1 = the new candidate.

| # | LP | Model | Expected | Why |
| --- | --- | --- | --- | --- |
| 1–2 | WT growth | M0, M1 | equal (provable) | Build sanity |
| 3–4 | KO YALI1F03452g | M0, M1 | M1 = M0 KO/WT (inferred) | M1 also closes R1964 when YlAraDH2 is absent (D1b) and keeps it open in D1a |
| 5 | KO YALI1E06616g | M1 (D1a only) | KO/WT = 1 (inferred) | New gene; OR partner open |
| 6 (optional) | KO both genes | M1 (D1a) | KO/WT = 1 (inferred; R1964 = 0 in all records) | Confirms the OR logic |

- **Budget:** 5 LPs (6 with the optional one; 4 for D1b).
- **Stop conditions:**
  - WT M0 ≠ M1 beyond 1e-9 relative: stop, since that is impossible for a GPR-only diff;
  - any KO/WT that changes the call at 1/5/10/15 %: stop and report the R1964 flux;
  - a non-optimal status: record the raw value. It is not biological evidence.
- **Not needed:** an FVA of R1964, since the bounds do not change, and a full KO screen, since only
  these genes' KOs can differ.

**Acceptance criteria:**
- WT is unchanged;
- the KO/WT values match;
- the whole-model diff is limited to R1964 and the new gene;
- the focused test and the unchanged suite pass;
- the audit finds no undeclared change.

### 7.4 Optional evidence steps (no model change)

- Download Wang 2020 Additional file 1 (`.docx` from Springer) to confirm that it holds only sequences.
  This needs the user's OK for the download.
- Obtain the Liu 2026 *JAFC* full text, which the user may have through library access. It may name
  the ArDHs screened and their products.

## 8. Limits and items for other lines

- **Limits a claim.**
  - The products are inferred (§5). The CGMCC7326 clones were compared with W29 only at the primer
    termini. No host-background control is described.
- **Affects the R2190 line.**
  - Add the C-4 limit line (§5.4).
  - The R2190 review's argument that broad specificity weakens the EC 1.1.1.250 transfer is weakened
    by Hallborn 1995 (§4.3).
- **Affects the GPR/isozyme line (candidates only, not examined).**
  - YlAraDH1/2 are also active on mannitol, D-sorbitol and xylitol in vitro. Whether they are OR
    members of R257, R801 or the NAD mannitol reactions was not examined. R1893 is NADP and does not
    match.
  - R513 is written with L-glucitol → L-sorbose. That was noticed in passing and not examined.
- **Affects the arabitol-pathway line.**
  - YlMDH2 (YALI1D24054g, NADPH) is the only gene whose deletion removed D-arabitol in vivo. The
    model has no NADPH arabitol route. Adding one would be a new reaction, outside this scope.
- **Not covered:** the direction of R1964 (kept reversible; Rhea 17921 is undirected), builds
  without the E5 chain, and xylose media.

## 9. Commands, scripts and sources

The scripts and fetched files are in the session scratchpad, not in the repository. Full SHA-256
values are in [sources_sha256.txt](sources_sha256.txt).

| Item | What it does / is |
| --- | --- |
| `trace_r1964.py` | R1964/R2190/R801/R803/R1980/R2041/R2042 fields, species and residuals; metabolite neighbourhoods; gene → reactions (main checkout `.venv`, cobra 0.30.0). Run on a56f3649 and the source model |
| `tools/fischer_regio.py` → `fischer_regio_out.txt` | Fischer bookkeeping and motif test (§5.1, §5.3). Pure Python, `python3 -I` |
| `tools/jats2txt.py` → `wang2020.txt` | Text of the Wang 2020 JATS XML (`dl_lit/pmc7576711.xml`) |
| Inline scans | local MNXref grep of MNXR154181/MNXR198271/ABTD*/pentulose MNXMs; recorded-flux scan (`"R1964": <number>` in `model/reports/**`: 109 values in 35 files, all 0); consensus-list grep (0 hits for F03452/F02211/E06616/E05643); `sd_leu.csv` grep; builder call order in `main.py` |
| `dl_lit/` | PMC XML of Wang 2020, Adamczyk 2023 (full text), Link 2005 and Wong 1993 (abstract only in PMC XML); PubMed abstracts of the 13 papers cited; esummary triage lists |
| `dl_db/` | NCBI Gene 2912324/2908139, XP_503594.1 FASTA, UniProt YALI0E05643g, KEGG yli:2912324 / EC 1.1.1.11 / R05604 / R01896, ENZYME 1.1.1.11, Rhea 17921 |
| Not read | Wang 2020 Additional file 1 (not downloaded); Liu 2026 full text (HTTP 403); Link 2005 and Hallborn 1995 full texts |

## 10. Proposed STATE.md entry (superseded)

Superseded on 2026-10-19. After the user approved the plan, a combined review-plus-implementation
entry was added to `model/STATE.md` in this worktree. See
[r1964_gpr_20261019](../r1964_gpr_20261019/REPORT.md). The text originally proposed here is kept
for reference:

> ## 2026-10-19: R1964 GPR review, investigation only
>
> Delegated read-only review; no model, data or code change, no LP. R1964 (D-xylulose + NADH ⇌ NAD +
> D-arabinitol, orphan) in a56f3649. Wang et al. 2020 (doi:10.1186/s13068-020-01815-8) measured the
> purified YALI1F03452g ("YlAraDH1") and YALI1E06616g ("YlAraDH2"; = XP_503594.1 = YALI0E05643g/Q6C6W8,
> identical, verified; not a model gene) turning over both R1964 substrates with NAD(H). Products were
> assigned from a 340 nm assay, not identified (verified). A Fischer-projection test of their activity
> pattern favours D-arabinitol C-4 oxidation (→ D-xylulose, R1964) over C-2 (→ D-ribulose, R2190)
> (inferred, moderate; D-ribulose untested). Side finding: R1964's MNXR154181/BiGG ABTD1 annotation
> names L-xylulose (verified); the equation-consistent IDs are MNXR198271/Rhea 17921/KEGG R05604/
> EC 1.1.1.11. After either loop fix (R2041 forward-only or R2190 → D-ribulose), R1964 is blocked
> under export-only pentose exchanges (provable). Plan awaiting approval: opt-in guarded step with a
> provisional GPR `YALI1F03452g or YALI1E06616g` plus the annotation fix; 5 LPs.
> Report: [r1964_gpr_review_20261019](reports/r1964_gpr_review_20261019/REPORT.md).
