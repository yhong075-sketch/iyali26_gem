# Gene search for 33 orphan reactions (triage groups 10, 7, 6), 2026-10-19

Model `E5_coq9_resp_pkg_r132_20261012.xml` (`a56f3649…`), newest candidate, not a designated
reference; read only. Scope, budget and the user's confirmation: [TASK.md](TASK.md). No model XML,
curation data or GPR was changed; no flux run; nothing committed.

Per-reaction table: [gene_search.csv](gene_search.csv) (33 rows; judgment inputs in
[curation_calls.json](curation_calls.json), assembled by `scripts/assemble_table.py`).

## Summary

Counts come from `gene_search.csv`, recomputed by the assembly script.

| Tier | n | Reactions |
| --- | --- | --- |
| A, *Y. lipolytica* gene characterised | 1 | R2013 |
| B, reciprocal-best-hit orthologue of a characterised fungal enzyme | 3 | R23, R2010, R1958 |
| C, lead only (family-level, paralogue-ambiguous or compartment-unsupported) | 13 | R1901; R989, R1024, R1027, R1125, R1292, R1337; R517, R984, R1273, R1274, R1338, R1973 |
| D, no defensible gene | 16 | R1308; R333, R721, R1048, R1078, R1080, R1090, R1132, R1154, R1178, R1252, R1347, R1959, R2003, R2075, R2076 |

Proposals (for a later, separately authorized curation change; nothing applied):

| Reaction | Proposed GPR | Tier | Gene: ID / verified name / function / evidence |
| --- | --- | --- | --- |
| R2013 GlcNAc kinase | `YALI1E24028g` | A | YALI1E24028g = YALI0E20207g (100 % identity, 464 aa, this session); **YlNAG5** (Flores & Gancedo 2015); N-acetylglucosamine kinase; *Y. lipolytica* experimental: disruption in PO1a abolished growth on GlcNAc, purified enzyme Km 0.5 mM for GlcNAc. Already in the model on R1946. |
| R23 KMTB + Glu → Met | `YALI1D01431g` | B | YALI1D01431g; no verified *Y. lipolytica* name; branched-chain aminotransferase (reciprocal best hit of *S. cerevisiae* BAT1/BAT2, 64–66 %); *S. cerevisiae* deletion genetics show Bat2 transaminates KMTB with glutamate (Pirkov et al. 2008). The ARO8 orthologue YALI1E24922g is a C-tier addition (aromatic donors only). |
| R2010 DHA kinase | `YALI1F12917g or YALI1E24532g or YALI1F02508g` | B | Three DAK genes (= YALI0F09273g, YALI0E20691g, YALI0F01606g); no settled names (papers and KEGG swap DAK1/DAK2); DHA/glyceraldehyde kinase; reciprocal BLAST to fungal DAKs (S. pombe Dak1 assayed). *Y. lipolytica* data are overexpression only. R2010 duplicates R281/R2090 on duplicate metabolites; merging them is the cleaner fix. |
| R1958 HGD | provisional `YALI2_D00583g` | B (gene ID blocked) | Old locus tag, no YALI1 locus in GCF_001761485.1 as downloaded; no verified name; homogentisate 1,2-dioxygenase; the suppressed RefSeq protein XP_065950358 is a reciprocal HGD (57 % to *A. fumigatus* HmgA). |

Main findings:
- **Group 7 (6 compartment copies): no copy gets a gene.** Each enzyme exists in W29 as a single
  gene on the gene-bearing reaction. TargetP-2.0 and MitoFates predict no mitochondrial
  presequence or signal peptide for any of them, so the copies are C. Two literature points
  apply: in *S. cerevisiae*, mitochondrial FAD synthesis is contested, and matrix CoA depends on
  import through Leu5.
- **R989:** DeepLoc calls YALI1D01489g mitochondrial (0.7506); TargetP and MitoFates do not. The
  colleague-owned lipid-unlump line holds a human-approved provisional assignment of this gene to
  R989. This run neither supports nor overturns it; it is reported for that line's owner.
- **R1901 asparaginase (group 6):** the yeast and *E. coli* asparaginases have no W29 homologue.
  The only Asparaginase_2-family gene, YALI1F38344g (not in the model), is closest to the protease
  Taspase1, so it stays a lead. A 326-aa *Y. lipolytica* asparaginase is reported (Darvishi et al.
  2019, abstract only, gene not named). Linking it to this 326-aa gene is inferred.
- **R1958 HGD:** the gene is missing from the current W29 RefSeq annotation. NCBI marks the
  protein suppressed ("removed as a result of standard genome annotation processing").
- **16 reactions are D.**
  - Six have nothing in W29 better than E = 0.3 even at E ≤ 10: R1308, R1048, R1090, R1347,
    R1252, R2075.
  - Seven hit only homologues of a different function or weak partial hits: R333, R721, R1132,
    R1154, R1178, R1959, R2076.
  - R1078, R1080 and R2003 have no characterised enzyme in Swiss-Prot to search with.

## 1. Method

1. **Databases** (built once, `$IYALI26_RESEARCH_ROOT/artifacts/tasks/orphan_gene_search_20261019/db/`,
   provenance in [db_provenance.json](db_provenance.json)):
   - W29 RefSeq proteome and GFF, GCF_001761485.1 (NCBI Datasets API, downloaded 2026-10-19 11:33):
     7,918 proteins, all with YALI1_ locus tags.
     - 1,045 of the model's 1,073 genes map.
     - The 28 that do not are 9 YALI0 IDs, 13 mitochondrial-genome YALI1M IDs, 3 `YALIUNK`
       placeholders and 3 mitochondrial ORF names; all are listed in the provenance file.
   - UniProt Swiss-Prot, release 2026_03 (575,748 entries).
   - BLAST+ 2.17.0 (conda env `seqsearch`).
   - Deviation, approved by the user in this session: `makeblastdb` failed because BLAST+ splits
     `-in`/`-db` paths at spaces and the research root contains one ("Ian wheeldon"). The
     databases were built and searched through a space-free symlink to the same folder
     (`scripts/finish_dbs.py`, which repeats the post-download steps of the skill's
     `build_curation_dbs` unchanged). Nothing was downloaded twice.
2. **Queries** ([queries.csv](queries.csv)):
   - 63 proteins in the main run: reviewed UniProt entries annotated with the reaction's own Rhea
     IDs, or its EC when the model has no Rhea ID
     ([uniprot_candidates.tsv](uniprot_candidates.tsv)), fungal entries first.
   - Where the Rhea search missed a known yeast enzyme, the yeast protein was added by gene name
     and marked "lead".
   - Two more were added after the literature check (marked in `queries.csv`): human ASRGL1 for
     R1901 and *A. nidulans* hlyA (unreviewed) for R1973.
   - R1078, R1080 and R2003 have no characterised protein for the exact reaction in Swiss-Prot,
     so no BLAST was possible.
3. **Forward search:** blastp against W29, E ≤ 1e-5, 20 targets ([blast_w29_hits.csv](blast_w29_hits.csv)).
4. **Reciprocal check:** up to 5 W29 candidates per reaction against Swiss-Prot, top 5 hits
   ([reciprocal_hits.csv](reciprocal_hits.csv)). A gene is accepted only if its best Swiss-Prot
   hits are the same enzyme and no closer paralogue of a different function exists.
5. **Follow-up** (`scripts/followup_blast.py`, [followup_literature_checks.tsv](followup_literature_checks.tsv)):
   - for reactions with no hit, the same queries at E ≤ 10
     ([followup_sensitivity_e10.csv](followup_sensitivity_e10.csv)) to tell "nothing similar"
     from "a weak hit below the cut-off" (the decision cut-off stays 1e-5);
   - KEGG `yli` EC and product-name leads (local copy in `$IYALI26_RESEARCH_ROOT/reference/kegg/`);
     the HGD lead was fetched from KEGG REST ([inputs/](inputs/)) and checked against Swiss-Prot
     ([followup_extra_candidates.csv](followup_extra_candidates.csv));
   - YALI0 genes named in the literature were mapped to W29 by sequence (UniProt CLIB122
     sequence vs W29; every mapping was 100 % identity over the full length).
6. **Group 7 targeting** (Task B method, rule fixed in TASK.md before any run):
   - TargetP-2.0 (Non-Plant), MitoFates (Fungi) and DeepLoc 2.0 (Slow model) on 6 full RefSeq
     sequences ([targeting_group7_input.fasta](targeting_group7_input.fasta)).
   - Merged in [targeting_group7.tsv](targeting_group7.tsv).
   - Raw outputs, read from the result pages, under
     `$IYALI26_RESEARCH_ROOT/artifacts/tasks/orphan_gene_search_20261019/predictors/`.
7. **Static model context** ([reaction_context.csv](reaction_context.csv)): equation, same-EC
   gene-bearing reactions and the producers of each substrate (criterion (c); no flux run).
8. **Literature:** two read-only agents searched *Y. lipolytica* first, then *S. cerevisiae* and
   other fungi. They cited only pages they opened and said when a page was abstract only or
   unread. I re-opened the two sources that carry a proposal myself: Flores & Gancedo 2015
   (PMC4376941) and Pirkov et al. 2008 (PubMed 18625006 abstract). The `sources_opened` column
   says who opened what.

Evidence tiers (skill definitions):
- **A:** *Y. lipolytica* gene characterised for the activity.
- **B:** reciprocal-best-hit orthologue of a characterised (preferably fungal) enzyme, about
  ≥ 30 % identity and ≥ 70 % coverage.
- **C:** family-level, paralogue-ambiguous, distant-taxon or compartment-unsupported (tentative;
  not proposed for a build).
- **D:** no defensible gene.

Only A and B rows carry a proposed GPR.

## 2. Results by reaction (short form; full evidence in `gene_search.csv`)

### Group 6 (reverse flux at the optimum)

- **R1901 L-asparaginase: C.**
  - The *S. cerevisiae* ASP1/ASP3 and *E. coli* AnsA queries have no W29 hit (best E = 0.3 at E ≤ 10).
  - Human ASRGL1 hits YALI1F38344g (28.7 %, E = 4e-10). This is the CLIB122 YALI0_F30723g, the
    only PF01112 (Asparaginase_2) gene in that proteome.
  - Reciprocal best hits are Taspase1 (threonine aspartase, a protease), then plant-type
    asparaginases, so it is not accepted.
  - The direction question stays open; it needs a ΔG or enzymology source.
- **R1308 uridine phosphorylase: D.**
  - No W29 homologue; no KEGG gene.
  - *S. cerevisiae* salvages uridine by Urk1 or the Urh1 hydrolase.
  - KEGG annotates a Urh1-type nucleosidase at YALI2_F00923g (annotation only).

### Group 7 (compartment copies)

All six are C: each enzyme exists, but no evidence places it in the copy's compartment.

| Copy | Gene on the gene-bearing reaction | Targeting (rule) | Literature on the yeast orthologue |
| --- | --- | --- | --- |
| R989 mi | YALI1D01489g (Fau1 RBH) | not mito (DeepLoc mito) | Fau1 found in purified mitochondria by large-scale MS only; GFP "ambiguous" |
| R1125 mi | YALI1D33889g (FAD1 RBH) | not mito | mitochondrial FAD synthesis contested (for: Bafunno 2004, Bruni 2021; against: Wu 1995, Tzagoloff 1996) |
| R1292 mi | YALI1D24649g (CAB4 RBH) | not mito | Cab4 cytoplasm/nucleus; matrix CoA needs Leu5 (Prohl 2001) |
| R1024 mi | YALI1F28274g (Nce103 RBH) | not mito | *S. cerevisiae* Nce103 has an intermembrane-space pool, not a matrix pool |
| R1027 ex | YALI1F28274g | no SP | *Candida* Nce103 reaches the cell wall without a signal peptide; R1027 is closed ([0, 0]) |
| R1337 cy | YALI1C09285g on R751 (phytase-like acid phosphatase) | secretory pathway by DeepLoc | the yeast intracellular TMP phosphatase is unidentified; PHO3 acts in the periplasm |

For R1337 the A. thaliana TH2 query also hits YALI1D21464g (CTO1 family, homologous to TH2's HAD
domain): a lead, but *S. cerevisiae* Cto1 is not known as a phosphatase.

### Group 10

- **A:** R2013. See the proposals table.
- **B:**
  - R23: YALI1D01431g.
  - R2010: three DAKs.
  - R1958: HGD, blocked on the gene ID.
- **C:**
  - R984 dUMP nucleotidase: YALI1F34187g is the SDT1 orthologue, but no opened source shows
    Sdt1/Phm8 activity on dUMP.
  - R1338 thiaminase II: YALI1E05103g is the single THI20-family member, and in *S. cerevisiae*
    only THI20 of three paralogues has the activity.
  - R1273/R1274 dGDP/dGTP phosphatases: GDA1/YND1 orthologues act in the Golgi lumen; no
    deoxy-substrate data.
  - R517 lactaldehyde DH: ALDH family; reciprocal hits are Ald4/Ald5.
  - R1973 3-MG-CoA hydratase: crotonase paralogues only. The hydratase half of *A. nidulans* hlyA
    finds no specific W29 orthologue.
- **D:**
  - R333 glutaminase, R721 serine acetyltransferase, R1048 CMP nucleosidase, R1132
    fructose-1-P kinase, R1252 methylglyoxal synthase: no W29 homologue, or only a
    different-function homologue.
  - R1154/R1178 guanosine/inosine kinase: the only hit is ribokinase.
  - R1090/R1347 thymidine kinase: no homologue; *S. cerevisiae* lacks one (Vernis et al. 2003,
    uncited statement).
  - R2075 methylmalonyl-CoA mutase: no homologue; the documented *Y. lipolytica* propionate route
    is the methylcitrate cycle.
  - R1959 MAAI: only GST/Ure2-like hits.
  - R2076: the reaction lacks ATP; the nearest hit is the oxalyl-CoA synthetase orthologue.
  - R1078/R1080/R2003: no characterised enzyme exists to search with.

## 3. Literature

Opened sources are listed per reaction in `gene_search.csv`. Items that change a call:

| Source | Opened | Used for |
| --- | --- | --- |
| Flores & Gancedo 2015, PLoS One 10:e0122135, [PMC4376941](https://pmc.ncbi.nlm.nih.gov/articles/PMC4376941/) | me + agent; first 100k of 125k characters | R2013 tier A |
| Flores, Ariño & Gancedo 2021, IJMS 22:13109, PMC8658026 (Europe PMC XML) | agent | YlNag5 also regulates the pathway (caveat) |
| Pirkov et al. 2008, FEBS J 275:4111, [PubMed 18625006](https://pubmed.ncbi.nlm.nih.gov/18625006/) | me + agent; abstract only | R23 tier B |
| Darvishi, Faraji & Shamsi 2019, Int J Biol Macromol 125:955, PubMed 30576739 | agent; abstract only (full text 403) | R1901 lead |
| Larroude et al. 2021, Microorganisms 9:838, [PMC8071058](https://pmc.ncbi.nlm.nih.gov/articles/PMC8071058/) | agent; full text | R1958 substrate source (4HPPD deletion abolishes HGA) |
| Wang et al. 2022, 3 Biotech, [PMC8934898](https://pmc.ncbi.nlm.nih.gov/articles/PMC8934898/) | agent; full text | R2010 (overexpression only; DAK naming conflict) |
| Onozuka et al. 2008, FEMS Yeast Res, PubMed 18028398 | agent; abstract only | R1338 (thiaminase II exclusive to THI20) |
| Wu et al. 1995 (PubMed 7799934); Tzagoloff et al. 1996 (8631763); Bafunno et al. 2004 (14555654); Bruni et al. 2021 ([PMC8470081](https://pmc.ncbi.nlm.nih.gov/articles/PMC8470081/)) | agent; abstracts, Bruni full text | R1125 compartment contested |
| Prohl et al. 2001, MCB, [PMC99563](https://pmc.ncbi.nlm.nih.gov/articles/PMC99563/) | agent; full text | R1292 (matrix CoA via Leu5) |
| Vögtle et al. 2012, MCP, PMC3518125 | agent; first 100k characters | R1024 (Nce103 in intermembrane space) |
| Kuznetsova et al. 2015, JBC, [PMC4513125](https://pmc.ncbi.nlm.nih.gov/articles/PMC4513125/) | agent; full text | R984 (no dUMP data) |
| Park et al. 2018, Biotechnol Biofuels 11:158 (Europe PMC XML) | agent; first 100k characters | R2075 (methylcitrate cycle) |

The agents used NCBI E-utilities or Europe PMC where PubMed or PMC pages were blocked. Sources
they could not open are named in the table rows as unread and support nothing.

## 4. Group 7 targeting

Rule (TASK.md): predicted mitochondrial if TargetP-2.0 mTP ≥ 0.5 or MitoFates ≥ 0.385; secretion
signal if TargetP SP ≥ 0.5. DeepLoc is recorded but not part of the rule.

| Gene (RefSeq) | Copy in question | TargetP mTP / SP | MitoFates | DeepLoc 2.0 | Rule |
| --- | --- | --- | --- | --- | --- |
| YALI1D01489g (XP_502280.2) | R989 mi | 0.163916 / 0.005370 | 0.013 | **Mitochondrion 0.7506** (thr. 0.6220), "mitochondrial transit peptide" | not mito |
| YALI1D33889g (XP_503279.1) | R1125 mi | 0.018164 / 0.003881 | 0.084 | Cytoplasm, Nucleus (mito 0.4498) | not mito |
| YALI1D24649g (XP_503031.2) | R1292 mi | 0.002325 / 0.364989 | 0.051 | Cytoplasm, Nucleus (mito 0.2479) | not mito |
| YALI1F28274g (XP_505708.1) | R1024 mi, R1027 ex | 0.032530 / 0.001841 | 0.046 | Cytoplasm (mito 0.2869, extracellular 0.2555) | not mito, no SP |
| YALI1C09285g (XP_501537.3) | R1337 cy vs R751 ex | 0.000001 / 0.000007 | 0.000 | Lysosome/Vacuole 0.8336, Golgi 0.8502, "signal peptide" | not mito, no SP |
| YALI1D21464g (XP_502941.2) | R1337 cy (TH2-like lead) | 0.084095 / 0.000110 | 0.077 | **Mitochondrion 0.8477**, "mitochondrial transit peptide" | not mito |

Jobs:
- TargetP `6AD665190020334089046654` (server stamp 20261019204442);
- MitoFates `20261020.033226.224850` (server clock);
- DeepLoc 2.0 `6AD6657B0020341E2DF93235`.

- **Verified:** values read from the official outputs this session. They are predictions, not
  localization evidence.
- The YALI1D01489g TargetP values reproduce the lipid-unlump record exactly (Other 0.8307,
  mTP 0.1639; job `6A8DF2140014D6282BEC5B97`); the sequence SHA-256 `d8cea1ec…` is identical.
- **Predictor disagreement:** for YALI1D01489g and YALI1D21464g, DeepLoc calls the protein
  mitochondrial while TargetP and MitoFates do not. By the fixed rule neither counts as
  mitochondrial. The disagreement is reported, not resolved.
  - The current YALI1D01489g protein carries a 23-aa N-terminal extension
    (`MSCTILQITVHITYFSHNHCTTS`) absent from the CLIB122 model Q6CAN2 (Q6CAN2 = residues 24–235,
    verified).
  - Whether DeepLoc's call depends on this extension was not tested; it would need a second
    predictor batch outside this budget.
- YALI1C09285g has no cleavable signal peptide by TargetP. A hydrophobic stretch near residue 90
  (`IVGLATAGVTVLLIVLYLLLP`) and DeepLoc's vacuole/Golgi call suggest a signal-anchored
  secretory-pathway protein. That it is anchored is inferred, not tested.

## 5. Flags outside this scope (not acted on)

- **Duplicate metabolites:** m1881 "Dihydroxy-Acetone" (KEGG C00184) and m490 "glycerone" are the
  same compound (both C3H6O3). m2048 "Glyceron phosphate" (no formula) duplicates DHAP m456.
  R2010 duplicates R281/R2090 because of this.
- **Questionable existing GPRs:**
  - The purine nucleoside phosphorylase YALI1F37714g sits on the pyrimidine phosphorolyses
    R264 and R761.
  - The BAT gene YALI1D01431g sits on the decarboxylations R26 and R27.
  - Golgi YND1-type genes sit on the cytosolic NDPases R583–R588, R1948 and R1949.
  - The SDT1 orthologue sits on the dCMP and dTMP nucleotidases R981 and R983, without
    deoxy-substrate data.
- **Mitochondrial CoA-synthesis copy:** R1291 transport and R1292 rest on no targeting evidence.
- **Annotation errors:**
  - R333 EC list names glutamate synthases.
  - R2013 EC 2.7.1.162 does not fit a 6-phosphate product.
  - R1132 lists EC 2.7.1.11.
- **R2076 chemistry:** an ATP-free thioester formation (inferred).
- **Duplicate gene models:** YALI1C27061g and YALI1C27665g are 100 % identical over 636 aa.
- **Lipid-unlump line:** R989's provisional assignment (colleague-owned) faces split predictor
  evidence; for that line's owner to review.

## 6. What each statement rests on

- **Verified this session (sequence runs and opened records):**
  - all BLAST and reciprocal results, the YALI0 → W29 mappings, the predictor values, and the
    model context (static);
  - NCBI status of XP_065950358 and GeneID 90949663;
  - the Flores 2015 and Pirkov 2008 statements, which I opened myself.
- **Opened by the literature agents in this session, not re-read by me:** the other literature
  rows. Their quotes were extracted by the fetch tool.
- **Inferred:**
  - the Darvishi asparaginase ↔ YALI1F38344g link (length match only);
  - the location of the missing HGD gene in the current assembly (coordinate offset; genome
    not downloaded);
  - the YALI1C09285g signal anchor;
  - R2076's chemistry being implausible.
- **Untested:** any flux effect of the proposals; *Y. lipolytica* localization of any protein;
  whether the HGD gene model is intact in the W29 genome.

## 7. Method limits and possible improvements

- **E-value cut-off and divergent families.** The Asparaginase_2 case shows that Rhea-annotated
  queries can all come from one family while the organism uses another. A Pfam/HMM search
  (hmmscan of W29 against the enzyme's Pfam families) would catch family switches. Test: rerun on
  R1901 and check that PF01112 finds YALI1F38344g.
- **Reciprocal rule with fused or multidomain proteins.** THI20-family and hlyA-type fusions make
  "best Swiss-Prot hit" ambiguous. A domain-resolved reciprocal check (align the domain carrying
  the activity only) would be more exact.
- **Annotation gaps.** HGD shows that a gene can drop out of a RefSeq annotation. A tblastn of
  every no-hit query against the W29 genome would find unannotated or frameshifted genes. It
  needs the genome FASTA (about 20 MB), which was outside this budget. Test: tblastn A. fumigatus
  HmgA against NC_090773.1 and expect a hit at 3.62 Mb.
- **Predictors.** TargetP/MitoFates and DeepLoc disagree for two proteins. Running the N-terminally
  truncated YALI1D01489g (residues 24–235) would show whether the extension drives the DeepLoc
  call.
