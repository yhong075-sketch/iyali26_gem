# YALI1D02625g on transketolase R765/R766: identity evidence and keep/remove plan (2026-10-19)

Scope and starting state: [TASK.md](TASK.md). Investigation only.
- No model, curation-data or code file was changed. Nothing was built, no LP was run, and nothing was
  committed.
- The only repository edits are this folder and a new entry at the top of `model/STATE.md` (§11).

Evidence tiers used here:
- **verified**: I opened the source this session (2026-10-19), or ran the computation this session,
  and the value is shown here or in a linked table.
- **recorded**: taken from project files and not re-checked.
- **inferred**: reasoned from verified facts; no single source states it.
- **AlphaFold prediction**: a computed structure model. It never confirms identity, catalysis,
  specificity or compartment.

## 1. Summary

- **Gene.** YALI1D02625g:
  - NCBI Gene [2910520](https://www.ncbi.nlm.nih.gov/gene/2910520), old locus YALI2_D00047g, CLIB122
    locus YALI0D02277g;
  - protein [XP_502319.1](https://www.ncbi.nlm.nih.gov/protein/XP_502319.1), 347 aa, single exon;
  - UniProt [Q6CAJ3](https://rest.uniprot.org/uniprotkb/Q6CAJ3) (TrEMBL, annotation score 1,
    protein existence 4 "Predicted").
  - Verified name: none. RefSeq calls it "uncharacterized protein", "weakly similar to … TKL1
    transketolase 1".
  - Function: a transketolase N-terminal (diphosphate-binding, PP) domain on its own. Evidence level:
    homology and structure prediction only. There is no *Y. lipolytica* enzymology.
- **Recommendation: remove YALI1D02625g from R765/R766**, so that their GPR becomes `YALI1E07744g`,
  the immutable source model's rule. This is a change to the OR claim only. The gene is real and
  expressed, so it stays in the model, unlinked, as YALI1E11370g did after the R671 fix.
- **Why (AGENTS.md §5).** An OR member needs evidence that it catalyses alone. None exists, and every
  line of evidence below argues against it:
  - It has no pyrimidine-binding (PYR) or C-terminal domain, so no catalytic Glu (yeast E418) (§3.1,
    verified).
  - Its ThDP/Mg²⁺-binding residues are not conserved. There is no GDG motif: Glu-Glu sits at the
    GDG position, and Asp at the ThDP-binding His69 (§3.2, verified alignment). The other ThDP His
    (H263) has no aligned residue; the structurally nearest one is Asp (AlphaFold prediction, §6).
  - The transketolase active site needs the PP domain of one subunit and the PYR domain of the other
    (Lindqvist 1992, verified).
  - **Inferred:** a PP-domain-only protein cannot form a complete active site by itself.
- **Not a split or mis-annotated gene.**
  - Nothing in ±30 kb encodes the missing domains (§4, verified).
  - No other W29 protein has transketolase homology except YALI1E07744g (BLAST, E ≤ 1e-3). By Pfam
    HMM, no CLIB122 protein carries a transketolase C-terminal domain except E07744g (verified).
  - Thirteen protein records from thirteen *Y. lipolytica* strain assemblies, annotated in separate
    runs, give the identical 347-aa protein (verified).
- **Expressed at low level.**
  - RNA-seq supports the gene model (Magnan 2016 S4 Table, verified).
  - Its protein is detected by MS in 14 of 18 samples, at about 0.6% of the transketolase signal
    (§5, verified computation on local data).
- **Evidence tier of the recommendation: inferred.** It rests on verified domain, alignment, locus
  and literature facts. It is not experimental proof that D02625g lacks activity.
- **Recall played no part.** YALI1E07744g is an essentiality FN in the 2026-09-10 screen, and the
  change would likely alter that call. That is a consequence to report, not a reason (§8).

## 2. How YALI1D02625g got onto R765/R766 (verified from repository history)

- `model/source/iyali26.xml`: R765 "transketolase 1" and R766 "transketolase 2" have GPR
  `YALI1E07744g` only.
- On 2026-06-02, commit `a9a16cb` added `data/gpr_isozyme_additions.csv`. That file is now
  `model/curation/gpr/gpr_isozyme_additions.csv`, rows 22–23: `R765/R766, YALI1D02625g, 2.2.1.1`.
  - The rows came from the "CLIB89 expansion funnel": S2 Table → KEGG metabolic relevance → MetaNetX
    → de-dup against model EC.
  - The match was by EC only. KEGG [yli:2910520](https://www.kegg.jp/entry/yli:2910520) assigns
    K00615 transketolase [EC:2.2.1.1] (verified). KEGG assigns KOs to genome genes computationally
    from sequence similarity (established, uncited), so the EC carries no activity evidence.
  - In `$IYALI26_RESEARCH_ROOT/reference/locus_map/s2_metabolic_genes.csv`, D02625g has the same KEGG
    pathways and EC as E07744g (verified, file read).
- No row records evidence that D02625g catalyses the reaction.
- On 2026-07-06, `remove_misannotated_gprs` then removed E07744g instead. That was the YALI1/YALI0
  locus-number mix-up corrected on 2026-10-19
  (`selection_gpr_review_20261019/REPORT.md` §3.1 in worktree `clever-kare-7d2e1a`).
- The metadata selection (`metadata_reaction_selection.json`) guards R765/R766 with GPR before =
  `YALI1D02625g` and after = `YALI1E07744g or YALI1D02625g` (verified).
- Every candidate checked has `YALI1E07744g or YALI1D02625g`: a56f3649, c98961e9 and 51c6fabd. No
  other reaction uses either gene (verified; SHAs in TASK.md).
- External models, context only and not evidence: Chalmers iYali4, iYali5 and eciYali5 (local copy,
  commit `816db7f0`) put transketolase on YALI0E06479g (= YALI1E07744g) only. They never mention
  YALI0D02277g (verified by text search).

## 3. Domain architecture and active-site residues

### 3.1 Domains (verified, UniProt 2026-09-02 release, read 2026-10-19)

| | YALI1D02625g (Q6CAJ3) | YALI1E07744g (Q6C6T4) |
| --- | --- | --- |
| Length | 347 aa | 694 aa |
| Pfam | PF00456 Transketolase_N (30–232) | PF00456, PF02779 Transket_pyr, PF22613 Transketolase_C_1 |
| CDD | cl01629 TPP_enzymes superfamily, 30–301 (NCBI) | cd02012 TPP_TK, cd07033 TPP_PYR_DXS_TK_like |
| PROSITE transketolase signatures | none | PS00801, PS00802 |
| Gene3D 3.40.50.970 (ThDP-binding fold) | 1 copy | 2 copies, plus 3.40.50.920 |
| Other | disordered 313–347 | "Proton donor" active site E419 |

### 3.2 Residue mapping

Method: a local alignment (BLOSUM62, gap −11/−1) to Swiss-Prot TKL1
[P23254](https://rest.uniprot.org/uniprotkb/P23254) (crystal structures, experimentally annotated
sites) and to Q6C6T4. Full table: [site_mapping.tsv](site_mapping.tsv) (verified computation).

| TKL1 residue (role per UniProt P23254) | E07744g | D02625g |
| --- | --- | --- |
| H30 (substrate; important for catalysis) | H31 | H28 conserved |
| H69 (ThDP) | H71 | **D63** |
| G116 (ThDP) | G118 | G112 conserved |
| D157 (Mg²⁺/Ca²⁺; the "GDG") | D159 | **E153** |
| G158 (ThDP; the "GDG") | G160 | **E154** |
| N187 (Mg²⁺/Ca²⁺, ThDP) | N189 | in an alignment gap |
| I189 (Mg²⁺/Ca²⁺) | I191 | in an alignment gap |
| H263 (ThDP; important for catalysis) | H265 | outside the alignment |
| E418 (proton donor; ThDP N1′) | E419 | **absent** (no PYR domain) |
| R359, S386, H469, D477, R528 (substrate) | present | **absent** |

- D02625g contains no "GDG" tripeptide anywhere (verified, string search). This motif begins the
  ~30-residue ThDP-binding motif "-GDG- … -NN-" of ThDP enzymes (Hawkins, Borges & Perham 1989,
  [PMID 2792374](https://pubmed.ncbi.nlm.nih.gov/2792374/), abstract verified).
- In yeast transketolase, the diphosphate binds the N-terminal domain through a metal ion. The
  pyrimidine ring sits in a pocket of the middle domain of the second subunit, so the active site
  spans the dimer interface (Lindqvist et al. 1992,
  [PMID 1628611](https://pubmed.ncbi.nlm.nih.gov/1628611/), abstract verified).
- E418 is essential for catalysis. E418A keeps about 0.1% and E418Q 2% of wild-type activity, while
  still binding ThDP (Wikner et al. 1994, [PMID 7798210](https://pubmed.ncbi.nlm.nih.gov/7798210/),
  abstract verified).
- **Inferred from these:**
  - D02625g alone lacks the residue that activates ThDP. It also lacks half of the active-site
    pocket.
  - Its diphosphate/metal site is altered at the GDG position.
  - It cannot form a complete transketolase active site, as a homodimer or alone.

### 3.3 Homology (verified computation)

Tables: [blastp_w29.tsv](blastp_w29.tsv) and [blastp_swissprot_top25.tsv](blastp_swissprot_top25.tsv).

- **D02625g vs E07744g:** 23.9% identity over 356 columns (D02625g 2–311 vs E07744g 3–343), E = 8e-24
  (pairwise).
- **Swiss-Prot, top 25 hits (E ≤ 1e-3):** all are full-length transketolases (664–741 aa) at 21–30%
  identity.
  - They align only to the N-terminal part of the subject, ending at residue ≤ 408.
  - Examples: *B. subtilis* tkt (E = 1e-16) and *Scheffersomyces stipitis* TKT (8e-15).
  - *S. cerevisiae* TKL1 is the eleventh hit (3e-12).
- **W29 proteome:** D02625g and E07744g hit only each other and themselves (E ≤ 1e-3). The W29
  RefSeq proteome has no other transketolase-family protein.
- **Pfam HMM matches** (InterPro API, CLIB122 proteome, taxonomy 284591; table
  [interpro_tk_domains_clib122.tsv](interpro_tk_domains_clib122.tsv), verified):
  - PF00456 Transketolase_N occurs only in Q6C6T4 and Q6CAJ3.
  - PF22613 Transketolase_C_1 occurs only in Q6C6T4.
  - PF02779 Transket_pyr occurs in Q6C6T4 and in three 2-oxo-acid dehydrogenase E1 proteins: PDH E1β
    Q6C4G4, BCKDH E1β Q6C2U7 and OGDH Q6C3M8. Those are the expected PYR-domain enzymes, each with
    its own partner.
  - So, outside the 2-oxo-acid dehydrogenases, no CLIB122 protein carries a PYR domain, and none
    besides E07744g carries the TK C-terminal domain. This is an HMM-level check on the CLIB122
    annotation, not on W29.

## 4. Locus: is the RefSeq model a truncated or split annotation? (verified)

- **Location:** W29 chromosome 1D, [NC_090773.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_090773.1),
  262,573–263,616 (+). The gene is one exon; its CDS is 1,044 bp (347 codons and the stop).
- **Annotation history:** RefSeq XP_502319 is still at version .1, and is identical to the
  GenBank W29 record AOW03467. The mRNA XM_502319.3 is "derived from an annotated genomic sequence"
  and is "incomplete on both ends", that is, CDS-only.
- **Neighbourhood:** the region ±30 kb (232,572–293,615) holds 27 other protein-coding genes (28 with
  D02625g; no RNA genes). None is transketolase-related: none hits E07744g in the proteome BLAST (§3.3).
  - Directly downstream lie YALI1_D02636g (XP_068138675.1) and YALI1_D02641g (XP_068138676.1). These
    are small two-exon genes added in the 2024 re-annotation, with no old locus tags.
- **tblastn** of E07744g and of TKL1 against the ±30 kb region ([tblastn_region.tsv](tblastn_region.tsv)):
  - The only significant hit is the D02625g ORF itself: E07744g 1–343 on region 29,999–30,934
    (frame +2), E = 8e-20.
  - No hit covers the PYR (E07744g ~357–533) or C-terminal (~534–694) domains. Every other hit has
    E ≥ 7.
  - So tblastn finds no frameshifted or adjacent C-terminal part near the gene.
- **Strains** ([strain_annotations.tsv](strain_annotations.tsv)): 13 protein records give a
  347-aa sequence identical to XP_502319.1, with identical start and end. They come from CLIB122,
  W29 (GenBank), DSM 3286, H222, A101, PO1f, YlCW001, FKP355 and the YB-392/419/420/566/567 isolates.
  - Their locus-tag schemes differ (YALI0_, YALI1_, YALI2_, YALIA101_, YALIH222_ and several
    `…DRAFT_` series), so these are separate annotation runs.
  - The strains are not all genetically independent: PO1f, for example, is a W29 derivative
    (established, uncited).
  - The truncation is therefore not an assembly or annotation artifact of one strain.
- **Limit:** the genome outside ±30 kb was searched only through the annotated proteomes (W29 by
  BLAST, CLIB122 by Pfam HMM, §3.3), not by tblastn. That needs the W29 genome FASTA
  (GCF_001761485.1, about 20 MB), which was not downloaded. An unannotated PYR/C-terminal partner
  elsewhere is therefore not excluded.
  - Even if one existed, the rule would become `YALI1E07744g or (YALI1D02625g and <partner>)`. That
    would need its own activity evidence. The OR as it stands would still be unsupported.

## 5. Transcript and protein evidence (verified computation on local copies)

Table: [expression_summary.tsv](expression_summary.tsv).

- **RNA-seq.** Magnan et al. 2016 S4 Table, the W29 annotation GFF3. Its local copy is in
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/argininosuccinate_identity_20260911/sources/` (sha256
  `622f05a0…b0a2`).
  - YALI1_D02625g has `status=Matches RNA-Seq data; depth=1213`. That is the 49th percentile of the
    5,144 genes with that status.
  - YALI1_E07744g has depth 28,192 (97th percentile), about 23× more.
  - The exact definition of `depth` was not found in the methods text read. The paper says the
    AUGUSTUS models were scored using "coverage" and agreement with RNA-Seq data.
- **Proteomics.** The `debiaggi2025` data in the local eciYali5-GEM repository (commit `816db7f0`):
  - `ggMatrix.tsv` (sha256 `c048db1c…34e7`) is a gene-group intensity matrix: three strains (PAR,
    CAR, OBE) × two conditions (e/n) × three replicates.
  - The group `B0I71DRAFT_85041;YALI0_D02277g` is separate from the transketolase group
    `B0I71DRAFT_128545;YALI0_E06479g`. It is detected in 14 of 18 samples.
  - Its median intensity ratio to transketolase is 0.0062 (range 0.0007–0.0086). It ranks in the
    bottom ~17% of the 4,771 groups.
  - `proteomics.tsv` gives a similar ratio, 0.0044, in 18 of 24 samples. It has the same 4,771
    proteins, so it may be the same measurements in another form. It is not counted as independent.
  - The study's paper was not opened. Strain and condition details are therefore unverified.
- **Meaning.** The gene is transcribed and translated. It is not a silent pseudogene, but it is a
  low-abundance protein. Expression says nothing about transketolase activity.

## 6. Structure: AlphaFold prediction (reused, not re-run)

- **Model.** AlphaFold DB [AF-Q6CAJ3-F1](https://alphafold.ebi.ac.uk/entry/Q6CAJ3) v6 (AlphaFold
  Monomer v2.0 pipeline, model date 2022-06-01). Its sequence equals XP_502319.1 residue by residue
  (verified).
- **No HPCC job.** AGENTS.md requires reusing an existing identical model, so no HPCC job was
  submitted. The comparison model is AF-Q6C6T4-F1 v6 (sequence = XP_503628.3). Output:
  [alphafold_db_comparison.txt](alphafold_db_comparison.txt).
- **Confidence.** Mean pLDDT 80.6. By region:
  - 30–232 (Pfam TK_N): 89.3;
  - 233–312: 77.4;
  - 313–347: 36.2 (UniProt "disordered").
  - Mean PAE within 30–232 is 4.1 Å, and 16.3 Å between 30–312 and the tail.
- **Fit to E07744g.** Superposing on sequence-aligned Cα pairs gives 231 of 300 pairs within 4 Å,
  with a core RMSD of 1.27 Å. **AlphaFold prediction:** D02625g adopts the transketolase PP-domain
  fold through about residue 310.
- **Reference ThDP/Mg²⁺ sites.** After superposition, the nearest D02625g Cα at each site of
  E07744g:

  | E07744g site | Nearest D02625g Cα |
  | --- | --- |
  | D159 (GDG) | E153 |
  | G160 (GDG) | E154 |
  | H71 | D63 |
  | N189 | E233 |
  | I191 | R234 |
  | H265 | D25 |

  None of these residues is conserved.
  - The GDG-position residues E153 and E154 are confidently placed (pLDDT 92 and 90), as is D63
    (78).
  - D25, R234 and E233 are less certain (pLDDT 64, 61 and 69).
- **Limits:** a monomer prediction cannot show ThDP binding, dimerization or catalysis. "Nearest Cα"
  is not a proven residue equivalence.

## 7. Orthologs outside *Y. lipolytica*

- **Within *Y. lipolytica*:** the identical protein is present in all 13 strain assemblies checked
  (§4).
- **Characterised proteins (Swiss-Prot, §3.3):** no reviewed protein with this architecture appears
  among the top 25 hits. Every hit is a full-length transketolase, aligned over its N-terminal part
  only, so none is an ortholog of a PP-only protein with known activity.
- **UniRef** (UniProt REST, 2026-10-19; `uniref/clusters_Q6CAJ3.tsv` in the artifacts folder):
  UniRef50_Q6CAJ3 has one member, Q6CAJ3 itself. No other UniProtKB protein of any species shares
  ≥50% identity over the cluster's overlap rule.
  - This says little about close relatives, because UniProtKB holds only 88 proteins from other
    *Yarrowia* species (same query date).
- **OrthoDB:** UniProt lists Q6CAJ3 in group 112493at4891 and Q6C6T4 in 106547at4891, so they are
  separate orthogroups. The OrthoDB API returned no data for these IDs, which probably come from an
  older release. Group content was therefore not checked.
- **NCBI nr, not returned.** A web blastp of XP_502319.1 against nr was submitted: Fungi, excluding
  *Y. lipolytica*, E ≤ 1e-5, 100 hits, RID `DDASFSWM014`, 2026-10-19 23:43 UTC. It was still
  `WAITING` at 00:34 UTC, after 52 min of polling. It was not resubmitted, so there is **no result**.
  - The only annotated non-*lipolytica* *Yarrowia* genomes in NCBI are Yarrowia sp. E02, C11 and B02
    (Fraunhofer IME, about 6,300 proteins each; NCBI Datasets, verified). Whether they or other
    yeasts carry a PP-only homolog is **unknown**.
- **Meaning:** no ortholog with measured activity exists to transfer function from.
  - Whether the fragment is lineage-specific is open.
  - This limits a claim but does not change §8: even a conserved PP-only homolog would need its own
    activity evidence before it could justify an OR.

## 8. Assessment against AGENTS.md §5

- **OR:** an OR member needs evidence that it catalyses alone.
  - None exists: there are no *Y. lipolytica* assays, and no characterised ortholog with this domain
    architecture (§7).
  - Domain architecture and residue chemistry argue against autonomous activity. The tier is
    **inferred**; its premises are verified.
- **AND:** an AND rule needs evidence of a complex or co-requirement. None exists either.
  - A D02625g–E07744g heteromer is conceivable in principle; that is speculation, with no data.
  - Even then, E07744g would be required, so the rule would not be an OR.
- **Evidence-supported rule:** `YALI1E07744g`, as in `model/source/iyali26.xml`. Identity is
  verified by homology; there is no *Y. lipolytica* enzymology for it either. The Chalmers models
  agree, which is context and not evidence.
- **This is not a deletion.** No reaction or metabolite is removed, so AGENTS.md's non-existence
  standard does not apply. The gene stays in the model.
- **Consequence, reported but not used.**
  - With D02625g off, a YALI1E07744g knockout blocks R765/R766. That may change the model's
    essentiality call for E07744g; whether growth falls below each cutoff needs an LP.
  - YALI1D02625g itself is not on the experimental list, so it counts as non-essential under the
    user-defined negative class. As an unlinked gene, its knockout cannot change growth.

## 9. Plan for approval

**Approved by the user on 2026-10-19**, with build location "Option A": the clever-kare session
commits its own work, this report is committed, `claude/clever-kare-7d2e1a` is merged into
`claude/cool-nightingale-a2f471`, and the change is implemented and built there. Stop before the six
LPs of step 8. The implementation has its own report.

**Recommended: remove YALI1D02625g from R765/R766 through an opt-in post-selection change set.** It
follows the `curation-change` skill and the 2026-10-19 precedent for R671.

0. **Prerequisite.** The change-set registry (`platform/scripts/gem_annotate/selection_gpr_restore.py`)
   and the R765/R766 builder correction exist only uncommitted in worktree `claude/clever-kare-7d2e1a`.
   - This change must be built on top of them: either commit/merge that work first, or carry this
     change into that worktree.
   - Where to do it is the user's call. I would not touch that worktree without that instruction.
1. **TASK.md** for the change in `model/reports/tkl_d02625g_gpr_<date>/`.
2. **Data:** new `model/curation/gpr/selection_tkl_fragment_gpr_removal.json`.
   - Change ID `SELECTION-GPR-TKL-D02625G-<date>`.
   - R765 and R766, each with:
     - cleared = the current state (selection "after": unchanged equation and bounds, GPR
       `YALI1E07744g or YALI1D02625g`);
     - restored = the same equation and bounds, GPR `YALI1E07744g`;
     - genes `[{"id": "YALI1E07744g", "ncbigene": "2912397"}]`;
     - evidence tier "inferred", the sources in this report, and an explicit note that no
       *Y. lipolytica* enzymology exists for either gene.
3. **Code:** register the change set in `CHANGE_SETS` with targets `("R765", "R766")`, behind its own
   flag (proposed `--remove-selection-tkl-fragment-gpr`).
   - The step requires a complete selection, fails closed on any other state, and adds no genes.
   - YALI1D02625g stays as an unlinked gene product, as YALI1E11370g did in 51c6fabd.
4. **Test:**
   - it applies on the cleared state;
   - a rerun, also on a reloaded model, is a no-op;
   - it refuses on `YALI1D02625g`, `YALI1E07744g and YALI1D02625g`, or a changed equation;
   - the three existing change sets still reproduce exactly.
5. **Build:**
   - A control build without the flag reproduces the base model (51c6fabd, or its committed
     successor) byte for byte.
   - A candidate build with the flag goes to a **new** file in `model/candidates/`.
6. **Whole-model diff vs the base:** only the R765/R766 GPR and notes change, plus YALI1D02625g's
   reaction links.
7. **Independent read-only audit**, then REPORT.md and a `model/STATE.md` entry.
8. **Regression LPs.** These need a separate authorization with a budget: WT, KO YALI1E07744g and
   KO YALI1D02625g on the base and the candidate (6 LPs).
   - The 15% benchmark (also 1%, 5% and 10%) is reported only as a consequence.
9. **Publish** only on an explicit push instruction.

**Alternatives considered:**
- **(b) Fix at the source.** Delete rows 22–23 of `gpr_isozyme_additions.csv` and add a selection
  amendment for R765/R766.
  - This changes the default build, and it edits two audited records (the 2026-06-02 table and the
    2026-09-10 selection).
  - The 2026-10-19 reviews deferred that kind of change for R671 as a separate decision.
  - Not recommended now; it could follow once the opt-in set is accepted.
- **(c) Keep the OR.** This fails AGENTS.md §5, with no evidence for autonomous catalysis.
- **(d) Optional, not needed for the decision:**
  - an HPCC AlphaFold-multimer run of D02625g with E07744g, which is pre-authorized for identity
    questions; even a confident heteromer would mean AND, not OR;
  - a whole-genome tblastn after downloading the W29 genome FASTA (§4 limit).

## 10. What is proven, checked or untested

- **Statically verified this session:**
  - the domain annotations, residue mapping, locus neighbourhood, tblastn, proteome BLAST, the Pfam
    HMM domain census of CLIB122, and strain identity;
  - the RNA-seq status and proteomics detection;
  - the AlphaFold DB sequence identity and superposition;
  - the provenance of the OR in the repository;
  - the GPRs of a56f3649, c98961e9 and 51c6fabd.
- **Inferred:** that D02625g cannot catalyse transketolase alone. This rests on verified domain and
  residue facts and verified literature on yeast TKL1. It is not tested by any assay.
- **Not obtained:** the NCBI nr ortholog search (§7). The job did not return.
- **Untested:**
  - any LP;
  - the effect on the E07744g knockout;
  - whether D02625g has some other function, which is unknown.
- **Possible improvements:**
  - The Pfam check covers the CLIB122 UniProt proteome, not the W29 RefSeq proteome. The 2024 W29
    re-annotation added genes (the XP_0681… series) that UniProt may not hold.
  - Test: run hmmscan with PF02779/PF22613 on `w29.faa`. The expected hits are the three
    E1 proteins and E07744g; any other hit needs review.
  - A whole-genome tblastn would cover unannotated ORFs (§4 limit).

## 11. `model/STATE.md` entry

Added to `model/STATE.md` in this worktree (uncommitted):

> ## 2026-10-19: identity of YALI1D02625g, the OR partner on transketolase R765/R766 (investigation + plan)
>
> **Authorization.** A task relayed from session `local_553e1a53` (the user's own words are not on record
> here): investigate, then propose keep/remove for approval. No model, data or code edits, no build, no LP,
> no commit. AlphaFold was pre-authorized but not run, because AlphaFold DB already holds an identical model.
>
> Nothing in the model, curation data or code changed. Details:
> this REPORT.md.
>
> **Results** (verified 2026-10-19 unless marked):
> - **Gene:** YALI1D02625g = NCBI Gene 2910520 = XP_502319.1 (347 aa) = Q6CAJ3 = CLIB122 YALI0D02277g.
>   No verified name. Only a transketolase N-terminal (PP) domain.
>   - No PYR or C-terminal domain, so no catalytic Glu (yeast E418).
>   - The GDG Mg²⁺/ThDP motif is Glu-Glu, and ThDP His69 is Asp.
> - **Locus:** not a split or truncated annotation.
>   - Nothing within ±30 kb encodes the missing domains (tblastn).
>   - No other W29 protein is transketolase-like.
>   - By Pfam HMM, no CLIB122 protein carries the TK C-terminal domain except E07744g.
>   - 13 strain assemblies give the identical 347-aa protein.
> - **Expression:** RNA-seq supports the gene model (Magnan 2016 S4). The protein is detected by MS at
>   about 0.6% of transketolase (eciYali5-GEM debiaggi2025 data, local copy).
> - **AlphaFold prediction** (AF-Q6CAJ3-F1 v6, reused): a transketolase PP-domain fold (core RMSD 1.27 Å to
>   E07744g's model over 231 Cα).
> - **Orthologs:** no characterised PP-only ortholog (Swiss-Prot); UniRef50 is a singleton. The NCBI nr
>   search (RID `DDASFSWM014`) did not return in 52 min, so homologs in other yeasts are unknown.
> - **Provenance:** the OR came from the 2026-06-02 CLIB89 EC-match funnel (`gpr_isozyme_additions.csv`
>   rows 22–23, KEGG K00615). The source model has `YALI1E07744g` only.
>
> **Decision status.** Approved by the user on 2026-10-19 (plan §9, built after committing and merging
> `claude/clever-kare-7d2e1a`); applied in a separate change, not here. The change: remove YALI1D02625g from
> R765/R766 (GPR → `YALI1E07744g`)
> through an opt-in post-selection change set, building on the uncommitted registry in
> `claude/clever-kare-7d2e1a`.
> - Tier: **inferred**. There is no evidence that D02625g catalyses alone (AGENTS.md §5), and none for an
>   AND.
> - The likely change of E07744g's FN call is a consequence, not the reason.
>
> **Gaps.**
> - **Limits a claim:** no whole-genome tblastn (the W29 genome FASTA was not downloaded); no nr
>   ortholog result; no *Y. lipolytica* enzymology for either gene; the debiaggi2025 paper was not
>   opened.
> - **Blocks this task's implementation:** the clever-kare registry and builder fix are uncommitted.

## 12. Files

- This folder:
  - [TASK.md](TASK.md);
  - [site_mapping.tsv](site_mapping.tsv);
  - [blastp_w29.tsv](blastp_w29.tsv);
  - [blastp_swissprot_top25.tsv](blastp_swissprot_top25.tsv);
  - [tblastn_region.tsv](tblastn_region.tsv);
  - [strain_annotations.tsv](strain_annotations.tsv);
  - [interpro_tk_domains_clib122.tsv](interpro_tk_domains_clib122.tsv);
  - [expression_summary.tsv](expression_summary.tsv);
  - [alphafold_db_comparison.txt](alphafold_db_comparison.txt).
- In `$IYALI26_RESEARCH_ROOT/artifacts/tasks/d02625g_identity_20261019/`:
  - `downloads/`: the NCBI, UniProt, KEGG and PubMed records, the AlphaFold DB models and PAE, the
    region FASTA and feature table, the InterPro, UniRef and OrthoDB responses, and the remote
    BLAST RID and poll logs (no result);
  - `scripts/`: `map_sites.py`, `af_analysis.py` and `remote_blast.sh`;
  - `SHA256SUMS`.
