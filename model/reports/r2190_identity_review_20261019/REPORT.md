# R2190 (D-arabinitol:NAD+ 2-oxidoreductase): product versus EC 1.1.1.250

Read-only. No model, data or code was changed, nothing was built, and no LP was solved. Scope, inputs
and SHAs are in [TASK.md](TASK.md).

Evidence tiers:
- **verified:** checked this session (2026-10-19) and the output is shown here;
- **recorded:** taken from a project record and not re-run;
- **provable:** follows mathematically from verified or recorded facts;
- **inferred:** reasoned, with no single source;
- **summary-read:** a page read through a summarizing fetch; quotes were not checked by hand.

## 1. Summary

- **Every identifier on R2190 says D-ribulose (verified).** KEGG R07134, Rhea 17389, MetaNetX
  MNXR197981, SEED rxn04932/rxn15657 and EC 1.1.1.250 all write
  D-arabinitol + NAD⁺ ⇌ **D-ribulose** + NADH (+ H⁺). The reaction's own name is the EC
  systematic name, "2-oxidoreductase", which means oxidation at C-2 and gives a ketose.
- **The model's product is a real D-arabinose node (verified).** `m454` resolves to MNXM48397 /
  CHEBI:17108, D-arabinose. It is also used by D-arabinose-specific reactions (R2041 kinase, R254 and
  R255 D-arabinose 1-dehydrogenases, R1063 transport). So this is not a mislabeled ribulose
  metabolite. R2190 points to a different entity than its own name and annotations: cause (b) of the
  mismatch rule. The error is already in `model/source/iyali26.xml`. The pipeline later attached
  the ribulose annotations to the unchanged equation.
- **The gene is the YALI1 copy of the characterized YlAraDH1 (verified).**
  - YALI1F03452g is GeneID 2908139 and XP_504895.1 (343 aa). Its CLIB122 counterpart is YALI0F02211g
    (UniProt Q6C367), with an identical sequence.
  - Wang et al. 2020 purified it ("YlAraDH1"). It is an NAD(H)-dependent polyol dehydrogenase that
    oxidizes D-arabitol, mannitol, D-sorbitol and xylitol, and reduces only ketoses: fructose,
    L-sorbose and D-xylulose.
  - None of the aldoses tested (ribose, mannose, xylose, glucose) was reduced. D-ribulose and
    D-arabinose were not tested.
- **Verdict: the product is the wrong side (inferred, strong; §5).** The EC is not contradicted, but
  it rests on homology only. Whether the enzyme oxidizes D-arabitol at C-2 (D-ribulose, EC
  1.1.1.250) or at C-4 (D-xylulose, EC 1.1.1.11) is open. Wang et al. name D-xylulose as the
  product, but their assay measured only NADH.
- **Consequence:** the model has no free D-ribulose (verified).
  - Writing D-ribulose as the product makes it a single-reaction metabolite, which forces R2190 to 0
    (provable).
  - That breaks the zero-sum loop R803 → R714 → R2042 → R2041⁻ → R2190⁻ → R1964⁻. It also collapses
    the closed-boundary ranges of R1964, R2041 and R2042 to 0 (provable from connectivity, §6).
  - All 109 recorded R2190 fluxes are 0 (verified scan). They are single optima, so no growth change
    is expected (inferred).
- **Recommended plan (needs approval): Option A, §7.** A guarded, opt-in post-selection step that
  replaces `m454` by a new D-ribulose[C_cy] in R2190. It is mass- and charge-neutral, keeps the GPR,
  EC, bounds and annotations, and needs 8 LPs to confirm.

## 2. R2190 in each state

| State | Equation (stored IDs) | Bounds | GPR | Annotations | Where |
| --- | --- | --- | --- | --- | --- |
| Source `iyali26.xml` | m122 + m1823 ⇌ **m10 (H⁺)** + m123 + m454 | [−1000, 1000] | YALI1F03452g | none; species have no formula | `model/source` (verified) |
| a56f3649 | m122 + m1823 ⇌ m123 + m454 | [−1000, 1000] | YALI1F03452g | kegg.reaction R07134; seed rxn04932, rxn15657; metanetx MNXR197981; ec-code 1.1.1.250. No Rhea ID | verified |

Stored species in a56f3649 (verified):
- NAD `m122` is C21H27N7O14P2, charge 0;
- NADH `m123` is C21H29N7O14P2, charge 0;
- D-arabinitol `m1823` is C5H12O5, charge 0 (MNXM1018, KEGG C01904);
- D-arabinose `m454` is C5H10O5, charge 0 (MNXM48397, CHEBI:17108).

The a56f3649 equation balances in elements and charge (residual {} and 0, verified). The H⁺ in the
source was dropped on the way. The step was most likely the generic charge-aware
`fix_proton_water_balance` (`platform/scripts/gem_annotate/metabolites.py:602` →
`microspecies.balance_protons_and_water`), because with neutral NAD/NADH the H⁺ is surplus. That is
inferred from the code; no build log was traced. R2190 is not in
`model/curation/selection/metadata_reaction_selection.json` (verified grep).

The annotations were most likely added by the gene-EC route (Strategy C of
`annotate_reactions_extended.py`: GPR → gene EC → MNXR). UniProt gives the gene EC 1.1.1.250, and
the source name lacks the "(D-ribulose-forming)" suffix that a MetaNetX name match would need. This
is inferred; it was not traced in a build log.

## 3. What the annotations resolve to (task 1)

| Identifier | Resolves to | Product | Checked | Tier |
| --- | --- | --- | --- | --- |
| KEGG [R07134](https://www.kegg.jp/entry/R07134) | "D-arabinitol:NAD+ 2-oxidoreductase (D-ribulose-forming)"; C01904 + C00003 ⇌ C00309 + C00004 + C00080; EC 1.1.1.250; KO K17738; DBLINKS RHEA 17392 | [C00309](https://www.kegg.jp/entry/C00309) D-Ribulose, C5H10O5 (D-arabinose is [C00216](https://www.kegg.jp/entry/C00216)) | KEGG REST | verified |
| Rhea [17389](https://www.rhea-db.org/rhea/17389) | "D-arabinitol + NAD(+) = D-ribulose + NADH + H(+)"; CHEBI:18333, 57540, 17173, 57945, 15378; EC 1.1.1.250; 4 UniProt enzymes | CHEBI:17173 D-ribulose | Rhea REST | verified |
| MetaNetX [MNXR197981](https://www.metanetx.org/equa_info/MNXR197981) | "1 D-arabinitol + 1 NAD(+) ⇌ 1 alpha-D-ribulose + 1 H(+) + 1 NADH"; balanced; EC 1.1.1.250; reference keggR:R07134 | MNXM1371691 α-D-ribulose (xref C00309) | web page and local MNXref 4.5 | verified |
| Rhea 17389 in MNXref 4.5 | maps to MNXR197982, not MNXR197981. The two differ only in the ribulose species (MNXM1371690 D-ribulose vs MNXM1371691 α-D-ribulose) | D-ribulose | local `reac_xref.tsv` | verified |
| SEED rxn04932 / rxn15657 | both xref to MNXR197981 | cpd00258 D-Ribulose / cpd05268 α-D-Ribulose | local `reac_xref.tsv`, `chem_xref.tsv` | verified |
| EC [1.1.1.250](https://enzyme.expasy.org/EC/1.1.1.250) | ENZYME: "D-arabinitol + NAD(+) = D-ribulose + NADH + H(+)". [ExplorEnz](https://www.enzyme-database.org/query.php?ec=1.1.1.250) systematic name "D-arabinitol:NAD+ 2-oxidoreductase (D-ribulose-forming)" | D-ribulose | ENZYME text, ExplorEnz page | verified |

**The colleague-owned audit agrees (read only).** Its row is
`r1521_kegg_mnx_frontier_source_audit.json` (2026-08-18), verdict `contradicted`: "R07134 produces
D-ribulose, whereas the model produces D-arabinose; do not repair this by adding H+". Its blocker
list carries "R2190: D-arabinose versus D-ribulose product conflict". That file belongs to the
lipid-unlump owner and was not edited.

**Mismatch-rule checks on the local side (verified).**
- `m454`'s ChEBI IDs 12915, 20914 and 4103 are secondary IDs of MNXM48397 (D-arabinose) in
  `chem_xref.tsv`.
- `m1823` (MNXM1018, CHEBI:18333) is D-arabinitol, the same substrate as all the external entries.
- No free D-ribulose exists in a56f3649. A name scan finds only D-ribulose 5-phosphate `m96`,
  S-methyl-5-thio-D-ribulose 1-P `m164` and the CdRP intermediate `m653`.

**Does any database have an NAD reaction from D-arabinitol to D-arabinose? (verified, local MNXref
4.5)** All 17 MNXref reactions that use D-arabinitol (MNXM1018) were listed.
- The NAD⁺-linked ones produce D-ribulose (EC 1.1.1.250), D-xylulose (EC 1.1.1.11), L-xylulose
  (BiGG ABTD1) or D-lyxose (SABIO 7818, EC 1.1.1.21).
- D-arabinose appears only with NADP⁺ (MetaCyc RXN-22816, aldose reductase EC 1.1.1.21) or O₂ (Rhea
  59776, EC 1.1.3.40).
- None of them is an EC 1.1.1.250 reaction.

## 4. Identity of YALI1F03452g (task 2)

### 4.1 Identifiers (verified)

| Field | Value | Source |
| --- | --- | --- |
| Locus tag, W29/CLIB89 | YALI1_F03452g, chr 1F, NC_090775.1 344214–345245 (−), one exon | NCBI Gene [2908139](https://www.ncbi.nlm.nih.gov/gene/2908139) (esummary) and the local `clib89_feature_table.txt` |
| Other alias | YALI2_F00952g (old locus tag); KEGG uses it as the symbol | NCBI Gene; KEGG [yli:2908139](https://www.kegg.jp/entry/yli:2908139) |
| Protein | [XP_504895.1](https://www.ncbi.nlm.nih.gov/protein/XP_504895.1), 343 aa, product "uncharacterized protein"; CDD Rossmann-fold NAD(P)-binding region 1–342 | NCBI efetch |
| CLIB122 counterpart | YALI0F02211g (local S2 map, Magnan 2016); UniProt [Q6C367](https://www.uniprot.org/uniprotkb/Q6C367) "D-arabinitol 2-dehydrogenase [ribulose-forming] (EC 1.1.1.250)", PE "Inferred from homology", annotation score 2; Pfam PF00106, PF13561 (SDR) | UniProt REST |
| Sequence identity | XP_504895.1 = Q6C367 = KEGG AASEQ, character for character (SHA-256 `852797f1…`) | compared in this session |
| KEGG orthology | K17738 D-arabinitol 2-dehydrogenase [EC 1.1.1.250]. yli:2908139 is the only *Y. lipolytica* gene in K17738 (`link/ko/yli`) | KEGG REST |

The YALI1 and YALI0 locus numbers are separate series: F03452 here is F02211 in CLIB122. The link
rests on the S2 map, the shared GeneID and the identical protein, not on the numbers.

### 4.2 Homology (verified, KEGG SSDB precomputed Smith-Waterman)

[SSDB](https://www.kegg.jp/ssdb-bin/ssdb_best?org_gene=yli:2908139) returns 11,624 hits. The 26
highest-scoring hits are all K17738; hit 27 (`smp:10810563`) has no KO. In total 95 hits carry
K17738. Hits to the characterized ribulose-forming enzymes:

| Hit | Identity | Overlap | SW score |
| --- | --- | --- | --- |
| *C. albicans* Ard (`cal:CAALFM_C600150WA`) | 41.9 % | 341 | 809 |
| *C. tropicalis* (`ctp:CTRG_03026`) | 41.6 % | 341 | 820 |
| *S. stipitis* ARD2 (`pic:PICST_65696`) | 42.8 % | 341 | 801 |

*S. cerevisiae* SPS19 (2,4-dienoyl-CoA reductase) scores 30.7 % (SW 247). Wang et al. 2020 report
54.73 % identity to *S. stipitis* XP_001385035 (their method; not reproduced here). No local
alignment was run: Biopython and BLAST+ are not installed in this environment.

### 4.3 *Y. lipolytica* enzymology: Wang et al. 2020 (verified, table parsed from the PMC HTML)

Source: Wang N *et al.* 2020, *Biotechnol Biofuels* 13:176,
[doi:10.1186/s13068-020-01815-8](https://doi.org/10.1186/s13068-020-01815-8),
[PMC7576711](https://pmc.ncbi.nlm.nih.gov/articles/PMC7576711/).

**Identity.**
- "YlAraDH1" is gene g1595.t1 of strain CGMCC7326, which the authors map to YALI0F02211g. They found
  it by BLAST with the *S. stipitis* D-arabinitol 2-dehydrogenase.
- Their cloning primers translate to `MSLFSLAKK…` and `…VDGGYTIW*`. These are exactly the termini of
  XP_504895.1 (verified). The internal sequence of their clone was not compared.

**Method.**
- His-tagged protein from *E. coli*, Ni-NTA, >95 % pure.
- Activity followed at 340 nm only, so products were assigned, not identified.
- Cofactor: "YlAraDH1, YlAraDH2, and YlXDH were defined as NAD(H)-dependent".

**Activities of YlAraDH1 (U/mg; product as named by the authors), Table 2:**

| Substrate | pH | Named product | U/mg |
| --- | --- | --- | --- |
| mannitol | 8.0 | fructose | 10.3 |
| D-sorbitol | 8.0 | L-sorbose | 7.5 |
| D-arabitol | 8.0 | D-xylulose | 5.3 |
| xylitol | 8.0 | D-xylulose | 3.6 |
| fructose | 6.0 | mannitol | 2.6 |
| L-sorbose | 6.0 | sorbitol | 2.4 |
| D-xylulose | 6.0 | arabitol | 1.4 |

**No activity** (data not shown) on ribitol, erythritol, glycerol, primary alcohols, ribose,
mannose, xylose, glucose and fructose-6-P. D-ribulose and D-arabinose are not in the substrate
list.

**In vivo.**
- ΔYlAraDH1 showed no significant difference in growth or polyol synthesis on 300 g/L glucose.
- YlAraDH1 mRNA was low compared with YlMDH2.
- Only ΔYlMDH2 (YALI0D18964g, NADPH-dependent) lost both mannitol and D-arabitol production.

### 4.4 Other literature

- **The EC 1.1.1.250 archetypes (verified, PubMed abstracts only).**
  - Wong *et al.* 1993, *J Bacteriol* 175:6314 ([PMID 8407803](https://pubmed.ncbi.nlm.nih.gov/8407803/)):
    the *C. albicans* ArDH catalyzes "D-arabitol + NAD <-->D-ribulose + NADH". The abstract proposes
    D-arabitol synthesis by dephosphorylating and reducing D-ribulose 5-P.
  - Quong *et al.* 1993, *BBRC* 196:1323 ([PMID 8250887](https://pubmed.ncbi.nlm.nih.gov/8250887/)):
    the *C. tropicalis* enzyme "is specific for DA" and yields D-ribulose.
  - YlAraDH1's broad polyol range (§4.3) differs from the "specific for DA" *C. tropicalis* enzyme.
    That weakens the transfer of the ribulose-forming specificity by homology (inferred).
- **Y. lipolytica makes arabitol (summary-read).** Szczepańczyk *et al.* 2023, *Sci Rep* 13:9129
  ([doi:10.1038/s41598-023-36152-x](https://doi.org/10.1038/s41598-023-36152-x)) report arabitol
  rising almost 8-fold when erythrose-reductase homologs are deleted. The summary found no
  enantiomer, pathway or arabitol dehydrogenase gene named.
- **Not read:**
  - Taratynova *et al.* 2023, *Appl Biochem Microbiol*
    ([doi:10.1134/S0003683823090120](https://doi.org/10.1134/S0003683823090120)): the "pAraDH2"
    promoter paper. Crossref gives no abstract, and the page redirected to a login. Which locus it
    means was not established.
  - The 2021 *Enzyme Microb Technol* D-arabitol fed-batch paper (search snippet only).

### 4.5 Gene record (AGENTS.md naming)

- **Systematic ID:** YALI1F03452g (NCBI YALI1_F03452g, GeneID 2908139, XP_504895.1; old tag
  YALI2_F00952g; CLIB122 YALI0F02211g / Q6C367).
- **Verified name:** no standard gene name is verified. NCBI calls it "uncharacterized protein".
  "YlAraDH1" is the name used by Wang *et al.* 2020. The UniProt/KEGG name is a homology label.
- **Function:** NAD(H)-dependent short-chain polyol dehydrogenase, active in vitro on D-arabitol,
  mannitol, D-sorbitol and xylitol, and reducing ketoses only (§4.3).
- **Evidence level:** direct in vitro activity of the purified *Y. lipolytica* protein, with
  products assigned, not identified. There is no deletion phenotype on glucose. EC 1.1.1.250
  (ribulose-forming) is homology-level. The model's GPR on R2190 is not experimental validation.

## 5. Which side is wrong (task 3)

Claim: **R2190's product (D-arabinose) is wrong.** These are the conditions it needs:

| # | Condition | Checked | Result |
| --- | --- | --- | --- |
| 1 | R2190's own identifiers resolve to a non-arabinose product | yes, all opened (§3) | D-ribulose in every one (verified) |
| 2 | `m454` really is D-arabinose, not a mislabeled ribulose | yes (§3) | D-arabinose; used by D-arabinose-specific reactions (verified) |
| 3 | The reaction name implies a ketose | yes | "2-oxidoreductase" = the EC 1.1.1.250 systematic name, C-2 oxidation → ketose (verified name; chemistry established, uncited) |
| 4 | The gene's own enzyme does not make aldoses | partly | only ketoses reduced; ribose, mannose, xylose, glucose not reduced; **D-arabinose never tested** (verified, Wang 2020) |
| 5 | No NAD-dependent arabinitol → arabinose enzyme is known for *Y. lipolytica* | partly | none in MNXref under NAD⁺ (verified); no *Y. lipolytica* report found. Not found is not proof of absence |

Conditions 1–3 are fully verified; 4 and 5 are partial. The verdict is therefore **inferred
(strong)**, not verified: no experiment has tested this protein with D-arabinose.

**The EC side** is consistent with conditions 1 and 3 and with homology (§4.2). It is not confirmed
by *Y. lipolytica* data:
- the purified enzyme is broad, unlike the archetype;
- the authors name D-xylulose, not D-ribulose, as the D-arabitol product (assay-only assignment);
- D-ribulose was never tested.

So **which ketose** is open: D-ribulose (EC 1.1.1.250, all databases, homology) or D-xylulose (EC
1.1.1.11, Wang 2020's naming). This review cannot settle it.

## 6. Consequences in a56f3649 (task 4)

**Network around R2190 (verified unless marked).**
- D-arabinitol `m1823` has three reactions: R2190, R1964 and R1980.
  - R1964 is D-xylulose + NADH ⇌ NAD + D-arabinitol, [−1000, 1000], no GPR, source note
    `PROTEIN_CLASS 1.1.1.11`.
  - R1980 is transport to `m1836`, whose only other reaction is the secretion-only exchange R1981.
- D-arabinose `m454` has R1063, R2041, R2190, R254 and R255.
  - R1063 is transport to `m1193`, whose only other reaction is the secretion-only exchange R1062.
  - R254 and R255 lead to the lactone `m452` (R253, R254, R255). R253's product `m453` has no other
    reaction, so that branch cannot carry steady-state flux (as in the R2041 review §3, recorded).
- D-xylulose `m1015` has R801 (xylitol → xylulose, forward-only, YALI1E15452g), R803 (xylulokinase,
  forward-only) and R1964.
- Neither D-arabinose nor D-arabitol is in the medium.

**The loop.** R803 → R714 → R2042 → R2041 reversed → R2190 reversed → R1964 reversed sums to zero.
- The R2041 review showed this statically (recorded).
- The orphan inventory's closed-boundary FVA on a56f3649 found R1964, R2041 and R2042 at ±1000
  (recorded, `orphan_reactions_20261019/loop_check.txt`, other worktree).
- R2190 is the only link from D-arabinose to the arabinitol node.

**If the product becomes a new D-ribulose[C_cy] (Option A):**
- **R2190 is forced to 0 (provable).** The new metabolite would have one reaction, so steady state
  fixes R2190 = 0, as if its bounds were [0, 0].
- **Feasible set shrinks (provable).**
  - The new LP's feasible set is contained in a56f3649's.
  - So closed-input ATP stays 0: ≤ 0, using the R2041 review's provable 0 for a56f3649.
  - WT growth is ≤ a56f3649's. Equality is expected (inferred): all 109 recorded R2190 values in
    35 report files are 0 (verified scan of `model/reports`). Those are single optima, not FVA.
- **The loop is gone (provable from connectivity).** With every boundary closed:
  - R1980 is 0 because R1981 is closed;
  - so R1964 is 0, since `m1823` has no other outlet;
  - R1063 is 0, and the lactone branch is 0;
  - so R2041 is 0, and with it R2042, since A5P `m1884` has only these two reactions.

  The closed-boundary ranges of all three orphan-loop reactions therefore become [0, 0].
- **Open network (inferred).**
  - D-arabinitol can still be made from D-xylulose by R1964 and secreted.
  - D-arabinose can still be made only by reverse R2041 from A5P. If the R2041 forward-only plan is
    also applied, D-arabinose has no internal source at all. That branch is already dead downstream
    (`m453`).
- **Balance (verified).** D-ribulose C5H10O5, charge 0, replaces D-arabinose C5H10O5, charge 0, so
  the element and charge residuals stay {} and 0.
  - No H⁺ edit is needed. Rhea's H⁺ belongs to its NAD(1−)/NADH(2−) convention; the model stores
    neutral NAD and NADH.
  - This agrees with the audit's "do not repair this by adding H+".
- **What a live D-ribulose would need (not proposed).** A D-ribulose 5-P phosphatase (Wong 1993
  infers one for *C. albicans*), a D-ribulokinase (EC 2.7.1.47), or transport. None was evidenced
  for *Y. lipolytica* in this review. Adding any of them to "activate" R2190 would be gap-filling
  without evidence.

## 7. Plan for approval

### 7.1 Options

| Option | What changes | For | Against |
| --- | --- | --- | --- |
| **A (recommended). Product → D-ribulose** | One new metabolite D-ribulose[C_cy]; in R2190 replace `m454` by it. GPR, EC, bounds and annotations kept | Makes the equation say what its name, EC and every cross-reference say; mass/charge-neutral; removes the loop; R2190 becomes blocked, which is honest given no D-ribulose route is known | Regiochemistry not proven (§5); adds a dead-end metabolite; one more step |
| B. Product → D-xylulose, then merge | R2190 becomes the exact reverse of R1964 → move the GPR to R1964 and drop R2190 as a duplicate | Matches Wang 2020's named product; gives the orphan R1964 a gene | Product assignment is assay-only; contradicts EC/KEGG/Rhea/homology; removes a reaction. Not recommended without product evidence |
| C. Keep D-arabinose, change the EC/annotations | Annotations only | Equation untouched | No evidence for an aldose-forming NAD enzyme from this gene (§5 rows 4–5); not recommended |
| D. Keep and note | A model note only | No change | The equation contradicts its own name and annotations; the loop stays |
| Not offered: delete R2190 | — | — | The gene exists and is active on D-arabitol (Wang 2020); AGENTS.md deletion rule |

### 7.2 Option A in detail (steps of the `curation-change` skill; nothing done)

1. **TASK.md** for the change, quoting the user's authorization.
2. **Data:** `model/curation/chemistry/r2190_product.json`, laid out like
   `r132_mito_adh/r132_mito_adh.json` / `direction/r1931_direction.json`.
   - `curation_id`: `R2190-RIBULOSE-<date>`; `status`: `user_authorized`, with the user's words.
   - Before: the exact a56f3649 equation (no H⁺), bounds [−1000, 1000], GPR `YALI1F03452g`, EC
     1.1.1.250, the name.
   - After: `m454[C_cy]` replaced by the new D-ribulose.
   - New metabolite:
     - ID: the next unused `mNNNN`; the highest in a56f3649 is m2060. Collision-checked at build.
     - Name `D-ribulose_C5H10O5`, following the model's `<name>_<formula>` pattern.
     - C5H10O5, charge 0, `C_cy`.
     - Annotations KEGG C00309, CHEBI:17173, MNXM1371690. These are Rhea's species; MNXR197981 uses
       the α-anomer MNXM1371691. Recorded as a note, not "fixed".
   - Evidence: the KEGG, Rhea, MetaNetX, ENZYME and ExplorEnz URLs, Wang 2020 and this report.
   - Notes:
     - `r2190_product_curation`;
     - `r2190_previous_product` (m454 D-arabinose, inherited from source);
     - `r2190_product_limits`: ribulose versus xylulose not resolved, Wang 2020 names D-xylulose,
       D-ribulose is a dead end in this model, EC is homology-level;
     - `r2190_product_authorization`.
3. **Code:** a guarded `apply_r2190_product(model, spec=None)`, as in `r132_mito_adh.py` /
   `r1159_direction.py`.
   - Fail closed unless all of these match: schema, curation ID, status, reaction ID, name, GPR, EC,
     bounds, and an exact before- or after-state stoichiometry and species.
   - The new ID must be unused, or present with the exact spec.
   - It changes only R2190's product and adds the metabolite, then asserts a zero element/charge
     residual. A rerun is a no-op.
4. **Wiring:** after the metadata selection and after the generic proton/water balancing, so the
   before-state is the H⁺-free equation. Two choices for the user:
   - (i) **opt-in** flag, for example `--r2190-ribulose-product` (recommended). Builds without it
     reproduce a56f3649 byte for byte.
   - (ii) **always on**. Every later build changes.
5. **Focused test** on an in-memory model:
   - applies on the a56f3649 state;
   - is idempotent;
   - fails closed on a stale equation, GPR, bounds, an ID collision or a wrong species;
   - flag off changes nothing.
6. **Build** a new file, for example `model/candidates/E5_coq9_resp_pkg_r132_r2190_ribulose_<date>.xml`.
   - Command: a56f3649's recorded `build.json` argv plus the flag.
   - Control: the same build without the flag must reproduce `a56f3649…`.
7. **Whole-model diff** against a56f3649: R2190's product, its notes and one new species only.
8. **Then:** `pytest` compared test by test, the LPs in §7.3, an independent read-only audit,
   REPORT.md and a STATE.md entry. No commit or push without an explicit instruction.

**Base.** This worktree on `285eaee` gives a one-change diff against a56f3649 (recommended).
Stacking on the uncommitted R2041 change (`bold-chaplygin-1c16df`) or on the GPR fixes
(`clever-kare-7d2e1a`, 51c6fabd) waits for them to land. The changes are independent: either one
alone breaks this loop.

### 7.3 LP checks (need authorization; none run)

- **Context:** that of `coq9_respiratory_package_20261008/run_controls.py`, as in the R2041 plan:
  PO1f SD-Leu, Gurobi with 1 thread, Presolve 0, tolerances 1e-7, 60 s per solve.
- **Models:** M0 = a56f3649; M1 = the new candidate.

| # | LP | Model | Expected | Why |
| --- | --- | --- | --- | --- |
| 1–2 | WT growth | M0, M1 | equal (inferred) | The change must not alter WT |
| 3–8 | Closed-boundary FVA (min, max) of R1964, R2041, R2042, with the orphan line's `loop_check` protocol | M1 | all [0, 0] (provable, §6) | Confirms that the build matches the connectivity argument and the loop is gone |
| P1–P2 (optional) | Closed-input ATP, C_cy / C_mi | M1 | 0, 0 (provable) | Direct record instead of the proof |

- **Budget:** 8 LPs (10 with P1–P2).
- **Stop conditions:**
  - WT M0 ≠ M1 beyond 1e-9 relative: stop and report R2190's role;
  - any closed range ≠ 0 beyond 1e-7: stop and report the witness flux;
  - a non-optimal status: record the raw value. It is not biological evidence.
- **Optional, separate authorization:** a full single-KO screen of M1, compared with the recorded
  a56f3649 screen at 1/5/10/15 %. Model genes not on the experimental list count as non-essential (a
  user-defined negative class). The screen is a development reference, not independent validation.
  YALI1F03452g's own KO cannot change, because R2190 is blocked in M1 and 0 in all recorded M0
  optima (inferred).

**Acceptance criteria:**
- WT is unchanged;
- the closed ranges are 0;
- the whole-model diff is limited to R2190 and the new species;
- the focused test and the unchanged suite pass;
- the audit finds no undeclared change.

## 8. Limits and items for other lines

- **Regiochemistry open (limits this claim).** Option A follows the databases. If product
  identification ever shows D-xylulose, Option B becomes right. Two reads that could help:
  - Wang 2020's supplementary material (not read);
  - the full text of Wong 1993, on whether *C. albicans* ArDH reduces D-xylulose (not read).
- **R1964 GPR (affects the orphan line; candidate).** Wang 2020 measured YlAraDH1 reducing
  D-xylulose (1.4 U/mg, NADH) and oxidizing D-arabitol to a product they call D-xylulose. That makes
  YALI1F03452g a candidate gene for the orphan R1964. The tier is in vitro, with the product
  assay-assigned. YlAraDH2 (YALI0E05643g → YALI1E06616g, not a model gene) has similar activities.
- **Arabitol biosynthesis route (affects another line).** In vivo, only ΔYlMDH2 (YALI0D18964g →
  YALI1D24054g, NADPH-dependent; not a model gene) lost D-arabitol and mannitol production (Wang
  2020). The model has no route for that.
- **Broad specificity (affects the GPR line; candidate only).** YlAraDH1 is also active on mannitol,
  D-sorbitol and xylitol in vitro. Whether it belongs as an OR member on the model's mannitol,
  sorbitol or xylitol dehydrogenase reactions (for example R1893, R257, R513, R801) was not
  examined. In vitro activity alone does not show a physiological role.
- **Colleague audit.** The lipid-unlump audit's R2190 row is consistent with this review. Whether to
  update its `blockers` list is the owner's decision.
- **Not covered:**
  - stereochemistry of `m1823`'s exchange partner `m1836`;
  - direction of R2190 (kept reversible; Rhea 17389 is undirected, and the Wong 1993 abstract
    writes ⇌);
  - builds without the E5 chain;
  - KO-specific use of R2190.

## 9. Commands, scripts and sources

The scripts and fetched files are in the session scratchpad, not in the repository. They can be
added under `scripts/` here on approval.

| Item | What it does / is | SHA-256 |
| --- | --- | --- |
| `trace_r2190.py` | R2190 fields, species, residual, pentose-name scan, all reactions of `m1823`/`m454` (main checkout `.venv`, cobra 0.30.0) | `35826a16…` |
| Inline scans | source-model R2190; local MNXref grep of MNXR197981/197982, all MNXM1018 reactions, MNXM/SEED/ChEBI lookups; recorded-flux scan (`"R2190": <number>` in `model/reports/**`); primer translation; sequence equality | — |
| `kegg_R07134.txt`, `kegg_ec_1.1.1.250.txt`, `kegg_K17738.txt`, `kegg_C00309.txt`, `kegg_C00216.txt` | KEGG REST | `d0bcc9d8…`, `6e12cd94…`, `68bf9ef0…`, `b26433fb…`, `33a71950…` |
| `enzyme_1.1.1.250.txt`, `explorenz_1.1.1.250.html` | ENZYME, ExplorEnz | `1807b1aa…`, `b9219902…` |
| `rhea_17389.tsv`, `mnx_MNXR197981.html` | Rhea REST, MetaNetX page | `6a9c7af4…`, `9310c5af…` |
| `kegg_yli_2908139.txt`, `ssdb_best_yli_2908139.html`, `yli_ko.txt` | KEGG gene, SSDB, yli→KO links | `abda8557…`, `7c04b59b…`, `e6a0f15b…` |
| `ncbi_gene_2908139.json`, `ncbi_XP_504895.1.gp`, `XP_504895.1.fasta` | NCBI E-utilities | `d47e030b…`, `6959e1a8…`, `f84b9171…` |
| `uniprot_yali0f02211.tsv`, `Q6C367.fasta`, `swissprot_ardh.tsv`, P43066/P50166/P50167/Q70FD1 FASTA | UniProt REST | `7032c8be…`, `dd477bf7…`, `f10b0bff…` |
| `pmc7576711.html` | Wang 2020 full text (tables parsed) | `b1f688a0…` |
| `pubmed_8407803_8250887.txt` | Wong 1993, Quong 1993 abstracts | `20b981ed…` |
| `crossref_S0003683823090120.json` | Taratynova 2023 metadata (no abstract) | `f36ac458…` |

Full SHA-256 values are in [sources_sha256.txt](sources_sha256.txt). The Szczepańczyk 2023 page was
read only through a summarizing fetch.

## 10. Proposed STATE.md entry (not written)

Not written: several worktrees already hold diverging uncommitted STATE.md heads (see the
2026-10-19 hand-over notes). Text to add when merged:

> ## 2026-10-19: R2190 product/EC identity review, investigation only
>
> Delegated read-only review; no model, data or code change, no LP. In a56f3649, R2190 is
> NAD + D-arabinitol ⇌ NADH + D-arabinose (inherited from source), but its name, EC 1.1.1.250 and
> every attached cross-reference (KEGG R07134, MNXR197981, SEED rxn04932/rxn15657; also Rhea 17389)
> give D-ribulose (verified). `m454` is a genuine D-arabinose node, so the product link is the error
> (inferred, strong). The gene YALI1F03452g = GeneID 2908139 = XP_504895.1 = YALI0F02211g/Q6C367
> (identical sequence; verified) is "YlAraDH1" of Wang et al. 2020 (doi:10.1186/s13068-020-01815-8):
> an NAD(H) polyol dehydrogenase that reduces only ketoses. D-ribulose and D-arabinose were untested,
> and the authors name D-xylulose as the D-arabitol product (assay-only). EC 1.1.1.250 is
> homology-level. Ribulose versus xylulose remains open. Plan awaiting approval: Option A, an opt-in
> guarded step `R2190-RIBULOSE-<date>` that writes a new D-ribulose[C_cy] (mass/charge-neutral).
> This blocks R2190 and removes the R803/R714/R2042/R2041/R2190/R1964 loop (provable). 8 LPs.
> Report: [r2190_identity_review_20261019](reports/r2190_identity_review_20261019/REPORT.md).
