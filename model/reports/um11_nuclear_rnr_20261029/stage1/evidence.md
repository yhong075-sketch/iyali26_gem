# UM11 stage 1 — evidence: nuclear ribonucleotide reductase and nuclear thioredoxin

Gathered 2026-10-29, read-only; no repository file edited.

Sources and how they were read:
- UniProt: REST `.txt`/JSON via curl (saved under `sources/`).
- NCBI: gene `esummary`, gene→protein `elink` (`gene_protein_refseq`) and protein `efetch`
  FASTA.
- PubMed: abstracts via NCBI E-utilities `efetch` (same text as the PubMed abstract); the
  PubMed page URL is given for clicking.
- SGD: GO cellular-component rows from
  `https://www.yeastgenome.org/backend/locus/<SGDID>/go_details`.
- Europe PMC: full-text XML for Timney 2016 (re-fetched this session).
- Local models: parsed read-only with `scripts/scan_models.py` and
  `scripts/model_inventory.py` (outputs `model_scan.txt`, `model_inventory.txt`).

Tiers:
- **verified**: fetched this task, and the claim is read in that text.
- **recalled**: not opened this task, or the fetch failed.
- **inferred**: reasoning only.

Model inspected (inventory only):
- `model/candidates/E5_coq9_…_um7_nuclear_diffusion_r2171_1to1_20261029.xml`
  (untracked in worktree `claude/wizardly-wing-a0d6c0`, HEAD ff5192d),
  SHA-256 `14f73bda4b0772019e82ea18a27e198c5ba175e67fef205a51dbfb090809e757`.
- `model/source/iyali26.xml` SHA-256 `5c8c199e…a99f9d`; `model/source/iyli21.xml` SHA-256
  `6974b758…e86631`.

Accession caveat:
- The model's W29 accessions A0A1D8ND06, A0A1D8NQQ7, A0A1D8N6C1, and the W29 database
  entries A0A1D8NLH6, A0A1D8NGE2, A0A1D8NMN9, are all **inactive** ("DELETED — Not part of a
  reference proteome"). **verified**, e.g. https://rest.uniprot.org/uniprotkb/A0A1D8ND06?format=json
- I mapped each gene to its active CLIB122 entry through its RefSeq protein
  (UniProt query `xref:refseq-XP_…`). **verified**, e.g.
  https://rest.uniprot.org/uniprotkb/search?query=xref:refseq-XP_502360&format=tsv
- NCBI Gene now names all six genes `YALI1_…`. Every one is "uncharacterized protein".
  **verified**: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=2910570,2908285,2906611,2908094,2908698,2911244&retmode=json

| Model gene | ncbigene | RefSeq (query) | aa | W29 DB entry (100% id) | Active UniProt (CLIB122) | UniProt name, status |
| --- | --- | --- | --- | --- | --- | --- |
| YALI1D04261g | 2910570 | XP_502360.3 | 849 | A0A1D8ND06 | Q6CAF2 (YALI0D03333g) | RNR large chain family, unreviewed |
| YALI1F39330g | 2908285 | XP_506109.1 | 892 | A0A1D8NQQ7 | Q6BZQ3 (YALI0F31735g) | RNR large chain family, unreviewed |
| YALI1B05072g | 2906611 | XP_500466.1 | 403 | A0A1D8N6C1 | Q6CFU6 (YALI0B03630g) | RNR small chain family, unreviewed |
| YALI1F02309g | 2908094 | XP_504863.2 | 129 | A0A1D8NLH6 | Q6C399 (YALI0F01496g), 104 aa | "Thioredoxin", unreviewed |
| YALI1F12969g | 2908698 | XP_505202.3 | 336 | A0A1D8NMN9 | Q6C2B0 (URA6, YALI0F09339g), 212 aa | Uridylate kinase, unreviewed |
| YALI1D35735g | 2911244 | XP_503348.2 | 434 | A0A1D8NGE2 | Q6C7L4 (TRR1, YALI0D27126g), 319 aa | Thioredoxin reductase, **reviewed** |

UniProt entries, all **verified** (each `.txt` fetched): https://rest.uniprot.org/uniprotkb/Q6CAF2.txt,
https://rest.uniprot.org/uniprotkb/Q6BZQ3.txt, https://rest.uniprot.org/uniprotkb/Q6CFU6.txt,
https://rest.uniprot.org/uniprotkb/Q6C399.txt, https://rest.uniprot.org/uniprotkb/Q6C2B0.txt,
https://rest.uniprot.org/uniprotkb/Q6C7L4.txt

- Three RefSeq proteins carry N-terminal extensions that the CLIB122 UniProt entries lack.
  **verified** (sequence comparison in `blast/`):
  - XP_504863.2 has 25 extra residues `MQRLVRPSVQLLTRHARQFHSTQHR`; Q6C399 is exactly residues
    26–129.
  - XP_503348.2 has 115 extra residues; the TRR1 alignment starts at its residue 116.
  - XP_505202.3 has about 124 extra residues; the Ura6 alignment starts at its residue 129.
- The thioredoxin extension has 5 Arg/Lys and no Asp/Glu, which resembles a mitochondrial
  presequence. **inferred**; no targeting predictor was run. Which start codon is used in vivo
  is **unknown**.

## Q1 Identity by reciprocal best BLAST hit

Method:
- blastp 2.17.0+ (`/opt/miniconda3/envs/seqsearch/bin/`), e ≤ 1e-5, max 10 targets.
- Reused databases, unchanged:
  - `$IYALI26_RESEARCH_ROOT/artifacts/tasks/um2_vacuolar_closures_20261028/stage1/blast/w29`
    (`w29.fasta`, 7894 seqs, SHA-256 `a0bacab595622864c5845411675bfb27881a4206f27572c7f0ae3abc8978e266`).
  - `…/sce` (`sce.fasta`, 6067 seqs, SHA-256 `1aaa81860c8bda68e81e4290703d871c7c8bd12a928038772fbc91b7e17b342a`).
- Queries:
  - Y. lipolytica: the six RefSeq proteins from NCBI `efetch`, in `blast/q_yli_refseq.fasta`,
    SHA-256 `00c1a09668c4771d7f8b3aa66a90481a060deac5e7c6d98834af589b6bcb66fb`.
  - S. cerevisiae: Swiss-Prot entries taken from `sce.fasta`, in `blast/q_sce.fasta`,
    SHA-256 `4a2f1595371609804cae77dd0f8efa579fa130c4cca534ee22f43c257676b6f9`.
- Each RefSeq query is 100% identical over its full length to one W29 entry
  (`blast/self_yli_refseq_vs_w29.tsv`). Reciprocity was tested against that entry
  (`blast/rbh_summary.tsv`).

Forward best hit, with the reverse best hit of that S. cerevisiae protein:

| Gene | Best S. cerevisiae hit | % id / aln len | qcov | e-value | bits | Reverse best (W29) | RBH |
| --- | --- | --- | --- | --- | --- | --- | --- |
| D04261g | Rnr3 (P21672) | 56.0 / 872 | 100 | 0 | 1007 | NQQ7 (= F39330g) | **no** |
| (D04261g) | Rnr1 (P21524) | 54.9 / 889 | 100 | 0 | 1006 | NQQ7 | no |
| F39330g | Rnr1 (P21524) | 67.5 / 892 | 99 | 0 | 1268 | NQQ7 | **yes** |
| B05072g | Rnr2 (P09938) | 66.2 / 408 | 100 | 0 | 546 | N6C1 | **yes** |
| (B05072g) | Rnr4 (P49723) | 53.1 / 337 | 84 | 9e-130 | 375 | N6C1 | one-way |
| F02309g | Trx2 (P22803) | 57.3 / 103 | 79 | 4e-42 | 133 | NLH6 | **yes** |
| (F02309g) | Trx1 (P22217) | 56.4 / 101 | 78 | 4e-40 | 127 | NLH6 | one-way |
| (F02309g) | Trx3 (P25372, mito) | 42.1 / 121 | 94 | 2e-37 | 122 | NLH6 | one-way |
| F12969g | Ura6 (P15700) | 65.0 / 200 | 60 | 4e-96 | 281 | NMN9 | **yes** |
| D35735g | Trr1 (P29509) | 76.8 / 319 | 74 | 0 | 513 | NGE2 | **yes** |
| (D35735g) | Trr2 (P38816, mito) | 73.4 / 327 | 75 | 0 | 509 | NGE2 | one-way |

Each of the following is **verified** from this task's BLAST tables.

- **YALI1F39330g** is the ortholog of Rnr1 by reciprocal best hit. Rnr3's best W29 hit is
  also F39330g.
- **YALI1D04261g** is a second RNR large-subunit paralog. It is the reciprocal best hit of no
  S. cerevisiae protein.
  - Its scores against Rnr3 (1007 bits) and Rnr1 (1006) are tied, and both proteins hit
    F39330g first.
  - It is 57% identical to F39330g.
  - It keeps every catalytic residue mapped from Rnr1 (P21524 features, global alignment,
    `blast/residue_check.txt`): C218, N426, C428 (radical), E430, C443, Y741, Y742, and the
    C-terminal shuttle cysteines C883/C886 (D04261g C216, N425, C427, E429, C442, Y738, Y739,
    C844, C847).
  - Both large subunits carry the ATP-cone domain (UniProt InterPro).
- Is D04261g an "RNR3"? That cannot be decided here.
  - Rnr1/Rnr3 is a whole-genome-duplication ohnolog pair, and Y. lipolytica branched before
    that duplication. Both points are **recalled**, not opened.
  - So D04261g is most likely an independent duplication, not an Rnr3 ortholog. **inferred**
  - Its expression, damage induction and activity are **unknown**. No Y. lipolytica paper was
    found (see Q2).
- **YALI1B05072g** is the ortholog of Rnr2 by reciprocal best hit. There is **no separate RNR4
  ortholog** in W29.
  - Rnr2 and Rnr4 each return only A0A1D8N6C1 as a W29 hit at e ≤ 1e-5.
  - B05072g keeps all six Rnr2 iron ligands and the radical tyrosine: D145, E176, H179, Y183,
    E239, E273, H276 map to D146, E177, H180, Y184, E240, E274, H277.
  - Rnr4 lacks three of those ligands (H179→Y127, E273→R220, H276→Y223) in the same alignment.
  - So B05072g is an Rnr2-type, iron-binding β subunit. **verified** by alignment.
  - Zhang 2011 abstract, **verified**: https://pubmed.ncbi.nlm.nih.gov/21931161/
    - In S. cerevisiae only β (Rnr2) binds iron and forms the tyrosyl radical.
    - β′ (Rnr4) is needed for iron loading in vivo.
  - Y. lipolytica therefore probably has a β2 homodimer, as in organisms without an Rnr4.
    **inferred**; no Y. lipolytica biochemistry was found.
- **None of the three RNR genes is a non-RNR protein.** Two are large (α) subunits and one is
  the small (β) subunit. **verified** (InterPro/PANTHER on the UniProt entries plus BLAST).
- **YALI1F02309g** is the ortholog of Trx2 by reciprocal best hit. It is also the best W29 hit
  of Trx1 and Trx3.
  - It has the WCGPC active site at residues 55–59 of XP_504863.2.
  - Other W29 thioredoxin-domain proteins (hits to Trx1 and Trx2; `blast/fwd_w29_trx_family_vs_sce.tsv`):
    - A0A1D8NJP0, 199 aa, WCGPC. This is model gene YALI1E27791g (XP_504314.1), annotated in
      the model.
    - A0A1D8NKT8, 324 aa, WCGPC.
    - A0A1D8NHR6, 131 aa, no WCGPC.
  - None of those three is a reciprocal best hit of Trx1 or Trx2. **verified**
- **YALI1D35735g** is the ortholog of Trr1 by reciprocal best hit. It is the **only** W29 hit
  for Trr1 or Trr2 at e ≤ 1e-5. **verified**
- **YALI1F12969g** is the ortholog of Ura6 (UMP/CMP kinase) by reciprocal best hit.
  **verified**

Real subunit composition (S. cerevisiae):
- UniProt P21524/P09938 SUBUNIT text, **verified**:
  https://rest.uniprot.org/uniprotkb/P21524.txt, https://rest.uniprot.org/uniprotkb/P09938.txt
  - The enzyme is a heterotetramer of two large (R1) and two small (R2) subunits.
  - Normally R1 is an Rnr1 homodimer; after damage it can be Rnr1–Rnr3.
  - R2 is an Rnr2–Rnr4 heterodimer. Rnr2 holds the iron-radical centre; Rnr4 is needed for
    Rnr2 folding and assembly.
- Zhang 2006 PNAS abstract, **verified**: https://pubmed.ncbi.nlm.nih.gov/16432237/
  - It describes R2 as a ββ′ heterodimer and R1 (α) as binding substrate and effectors.
- Consequence for the model GPR "D04261g or F39330g or B05072g":
  - It treats a large subunit alone, or a small subunit alone, as sufficient.
  - Catalysis needs α and β together, so a biologically consistent rule has the form
    (large) **and** (small). **inferred** from the verified composition.
  - Whether D04261g can stand in for F39330g as the large subunit is **unknown**.
  - The UDP copies R1913/R1914 use "(F39330g or (B05072g and F02309g))", which has the same
    problem.
  - These are GPR **candidates** only; nothing was tested.

## Q2 Where RNR is active in yeast

**S. cerevisiae** (all **verified** abstracts):
- Yao 2003 PNAS: https://pubmed.ncbi.nlm.nih.gov/12732713/
  - In a normal cell cycle Rnr1 and Rnr3 are mainly cytoplasmic and Rnr2 and Rnr4 mainly
    nuclear.
  - Under genotoxic stress Rnr2/Rnr4 move to the cytoplasm, depending on the checkpoint.
- Lee & Elledge 2006 Genes Dev: https://pubmed.ncbi.nlm.nih.gov/16452505/
  - Wtm1 anchors Rnr2/Rnr4 in the nucleus outside S phase.
  - On damage they become cytoplasmic and join Rnr1 to form the intact complex.
- Zhang 2006 PNAS: https://pubmed.ncbi.nlm.nih.gov/16432237/
  - Nuclear localization of ββ′ needs Kap122 and Wtm1.
  - Redistribution leads to "colocalization of R1 and R2".
- An 2006 Genetics: https://pubmed.ncbi.nlm.nih.gov/16489218/
  - Rnr2 and Rnr4 shuttle as one complex. Their relocalization involves active export and
    blocked import.
- Lee 2008 Mol Cell (Dif1): https://pubmed.ncbi.nlm.nih.gov/18851834/
  - In S phase or after damage, Rnr2–Rnr4 enters the cytoplasm to bind Rnr1, "forming an
    active complex".
  - Dif1 drives nuclear import. Dun1 phosphorylates Dif1, which inactivates and degrades it.
- Wu & Huang 2008 MCB: https://pubmed.ncbi.nlm.nih.gov/18838542/
  - Dif1 is cytoplasmic and mediates R2 nuclear import.
  - In S. pombe, Spd1 holds R2 in the nucleus "away from the cytoplasmic R1", which blocks
    holoenzyme formation.
- Wu 2018 MCB: https://pubmed.ncbi.nlm.nih.gov/29263158/
  - In S phase, Clb6–Cdc28 phosphorylates Rnr2, releasing Rnr2–Rnr4 from Wtm1. The complex
    relocates to the cytoplasm to form an active holoenzyme with Rnr1.
- Sanvisens 2016 JBC: https://pubmed.ncbi.nlm.nih.gov/26970775/
  - Iron deficiency also moves Rnr2–Rnr4 to the cytoplasm, through Dun1 and Dif1.
- S. pombe, Takahashi 2007 NAR: https://pubmed.ncbi.nlm.nih.gov/17264117/
  - Nucleoplasmic Spd1 holds Suc22 (R2). Releasing it forms the active complex in the
    cytoplasm.
- SGD GO, **verified**:
  - RNR1: https://www.yeastgenome.org/backend/locus/S000000872/go_details
    - cytoplasm IDA (Yao 2003; Huang & Elledge 1997); nucleus HDA only (Tkach 2012, a
      replication-stress GFP screen).
  - RNR2: https://www.yeastgenome.org/backend/locus/S000003563/go_details
    - nucleus IDA and cytoplasm IDA.
  - RNR4: https://www.yeastgenome.org/backend/locus/S000003412/go_details
    - nucleus IDA (Yao 2003; Huang & Elledge 1997) and cytoplasm IDA.
  - RNR3: https://www.yeastgenome.org/backend/locus/S000001328/go_details
    - cytoplasm IDA; mitochondrion HDA.
  - All four carry "is active in cytoplasm" (RCA, SGD 2025). None carries "is active in
    nucleus".
- UniProt P21524 (Rnr1) location: "Cytoplasm", ECO:0000269 (Yao 2003). P09938/P49723 note:
  mainly nuclear normally and moved to the cytoplasm on damage. **verified** (UniProt `.txt`).

**Mammals** (comparison only; **verified** abstracts):
- Pontarin 2008 PNAS: https://pubmed.ncbi.nlm.nih.gov/18997010/
  - Three independent methods place R1, R2 and p53R2 in the cytosol, during proliferation
    and after DNA damage.
  - The authors propose that dNTPs "diffuse into the nucleus".
- Leeds 1985 MCB: https://pubmed.ncbi.nlm.nih.gov/3915777/
  - Nuclear dNTP pools exist, but no reductase activity was detected in nuclear preparations
    at any cell-cycle stage.
- Counter-evidence: Niida 2010 Genes Dev: https://pubmed.ncbi.nlm.nih.gov/20159953/
  - In human cells, RRM1 and RRM2 accumulate at DNA-damage sites, recruited by Tip60, and
    this matters for G1 repair.
  - The paper argues for local nuclear dNTP supply during repair. It does not measure
    nuclear catalysis in yeast.

**Y. lipolytica:**
- **No RNR localization or biochemistry data found.**
  - PubMed "Yarrowia lipolytica AND ribonucleotide reductase" returned 0 hits;
    "… AND (deoxyribonucleotide OR dNTP)" returned 0.
- Hydroxyurea arrests Y. lipolytica in S phase (used for gene targeting).
  Tsakraklides 2015, **verified** abstract: https://pubmed.ncbi.nlm.nih.gov/26192309/
  - That hydroxyurea acts on RNR is **established/uncited**; the paper says nothing about
    localization.
- Regulator orthologs, blastp of S. cerevisiae protein vs W29 at e ≤ 1e-3
  (`blast/rev_regulators_vs_w29.tsv`), **verified**:
  - No hit for Sml1 or Dif1.
  - Wtm1/Wtm2 hit only WD40 proteins whose best S. cerevisiae hits are Hat2 and Mdv1, so they
    are not orthologs.
  - Kap122 has a weak reciprocal hit (A0A1D8NKJ1, 21% identity).
  - Dif1/Sml1/Spd1-family proteins are short and poorly conserved, so a missing BLAST hit does
    not show the regulators are absent. **inferred**

**Bottom line:**
- In S. cerevisiae the holoenzyme assembles in the cytoplasm. Nuclear R2 is a sequestered,
  inactive pool, and Rnr1 is cytoplasmic. **verified**, several primary sources.
- No source found places RNR catalysis in the yeast nucleus. The one nuclear-recruitment
  report (Niida 2010) is from human cells, and another human study (Pontarin 2008) finds RNR
  cytosolic.
- Y. lipolytica has no data.

## Q3 Thioredoxin: nucleus, size, role as RNR reductant; nuclear thioredoxin reductase

**Localization (S. cerevisiae):**
- SGD TRX1, **verified**: https://www.yeastgenome.org/backend/locus/S000004033/go_details
  - nucleus HDA (Yofe 2016) and nucleus IEA; cytosol IDA (Xu & Wickner 1996).
  - mitochondrial intermembrane space IDA (Vögtle 2012); mitochondrion HDA.
- SGD TRX2, **verified**: https://www.yeastgenome.org/backend/locus/S000003441/go_details
  - nucleus IEA only; cytosol IDA (Xu & Wickner 1996).
- UniProt, **verified**: https://rest.uniprot.org/uniprotkb/P22217.txt and
  https://rest.uniprot.org/uniprotkb/P22803.txt
  - Both Trx1 and Trx2: "Nucleus" and "Cytoplasm" with ECO:0000269, citing PubMed 14562095
    (Huh 2003 GFP atlas).
  - **Huh 2003 itself was not opened** (recalled; https://pubmed.ncbi.nlm.nih.gov/14562095/).
  - P22217 also says most yeast thioredoxin sits in the LMA1 complex with Pbi2 rather than as
    a free monomer.
- Pedrajas 1999 JBC, **verified** abstract: https://pubmed.ncbi.nlm.nih.gov/10037727/
  - Trx1, Trx2 and Trr1 are cytoplasmic; Trx3 and Trr2 are mitochondrial, with N-terminal
    import signals.
- Boisnard 2009, **verified** abstract: https://pubmed.ncbi.nlm.nih.gov/19581440/
  - Calls Trx1/Trx2 the "cytoplasmic thioredoxins". They are needed for H2O2-induced nuclear
    accumulation of Msn2/4 and Maf1.
- Dardalhon 2012, **verified** abstract: https://pubmed.ncbi.nlm.nih.gov/22561702/
  - The nuclear glutathione redox state is highly reducing and similar to the cytosol.
  - Trr1/Trx1/Trx2 mutants change compartment-specific redox states. This does not show
    nuclear thioredoxin directly.

**Size and passive diffusion:**
- S. cerevisiae Trx1 is 103 aa, 11,235 Da; Trx2 is 104 aa, 11,204 Da (UniProt SQ lines).
  **verified**
- Y. lipolytica Q6C399 (104 aa) is 11.34 kDa; the RefSeq/W29 129-aa form is 14.4 kDa
  (computed, Biopython ProtParam). **verified** computation.
- Timney 2016 J Cell Biol, full text **verified** this session:
  https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5057280/fullTextXML
  (abstract https://pubmed.ncbi.nlm.nih.gov/27697925/)
  - The NPC is a soft barrier with no firm 30–60 kDa cut-off.
  - In vivo in yeast, 27-kDa GFP equilibrated between nucleus and cytoplasm in seconds
    (τ ≈ 15 s). 67-kDa GFP-6PrA took minutes (τ ≈ 413 s).
  - The tested range was 27–230 kDa. Thioredoxin (11–14 kDa) is below that range, so its
    passive equilibration in seconds or faster is an extrapolation. **inferred**
  - If Trx is mostly in the LMA1 heterodimer, the diffusing particle is larger. The size of
    that complex was not checked. **unknown**

**Thioredoxin as the in vivo RNR reductant** (both **verified** abstracts):
- Camier 2007 FRBM: https://pubmed.ncbi.nlm.nih.gov/17349928/
  - Without thioredoxins, cells accumulate oxidized R1 and have lower S-phase dNTPs.
    Overexpressing R1 partly rescues this.
  - The authors conclude thioredoxins are the major RNR reducing system. Either the Trx or the
    Grx pathway must be present for viability.
- Koc 2006 JBC: https://pubmed.ncbi.nlm.nih.gov/16574642/
  - In Δtrx1 Δtrx2 cells, S-phase dNTP pools are about 60% of wild type. This is in vivo
    evidence that Trx donates electrons to RNR.
- UniProt P21524 marks C883/C886 as interacting with "thioredoxin/glutaredoxin".
  **verified**
- So glutaredoxin is the recognised backup reductant. **verified** (Camier), qualitatively.

**Nuclear thioredoxin reductase:**
- SGD TRR1, **verified**: https://www.yeastgenome.org/backend/locus/S000002761/go_details
  - cytosol IDA (Kang 2008); mitochondrial intermembrane space IDA (Vögtle 2012);
    mitochondrion HDA.
  - **No nucleus annotation of any evidence type.**
- UniProt P29509: "Cytoplasm" with no code; "Mitochondrion intermembrane space" ECO:0000269.
  **verified**: https://rest.uniprot.org/uniprotkb/P29509.txt
- Trr1 is a homodimer of about 34 kDa per subunit (UniProt). As a dimer it would be near the
  67-kDa class that equilibrates in minutes in Timney 2016. **inferred**
- Y. lipolytica Q6C7L4 (reviewed): "Cytoplasm {ECO:0000250}" (by similarity to P29509); GO
  cytosol IBA. **verified**
- Y. lipolytica literature:
  - Mauersberger 2001, **verified** abstract: https://pubmed.ncbi.nlm.nih.gov/11489863/
    - A viable insertion mutant (N040, alkAc) in "the thioredoxin reductase gene" grows poorly
      on long-chain alkanes.
    - The locus is not given in the abstract. The full text failed (Europe PMC 500; NCBI PMC
      gave only the abstract), so it is **unread**.
  - Thevenieau 2007, **verified** abstract: https://pubmed.ncbi.nlm.nih.gov/17071118/
    - Also reports a thioredoxin-reductase insertion among the AlkD/E mutants.
  - D35735g is the only TRR in W29 (Q1), so the mutant gene is probably D35735g, and
    Y. lipolytica can live without Trr on those media. Then RNR would have to be reduced by a
    non-Trr route such as glutaredoxin. **inferred**
- Comparison: Arabidopsis has a partly nuclear NTRA plus nucleoredoxins.
  Marchal 2014, **verified** abstract: https://pubmed.ncbi.nlm.nih.gov/24253198/
  - This is a plant, not evidence for yeast.

**Bottom line:**
- Yeast Trx1/Trx2 are reported in the nucleus only by high-throughput GFP (HDA, and UniProt's
  reading of Huh 2003). They are otherwise cytosolic.
- At 11–14 kDa they are well below the size that diffuses through the NPC in seconds. This is
  inferred from Timney 2016.
- Thioredoxin is the main in vivo RNR reductant, with glutaredoxin as backup.
- **No nuclear thioredoxin reductase is reported** in S. cerevisiae or Y. lipolytica.

## Q4 Do small nucleotides cross the nuclear pore by passive diffusion?

- Timney 2016, full text **verified** (link above).
  - The introduction defines passive diffusion as ions, small molecules and small to medium
    macromolecules equilibrating between cytoplasm and nucleoplasm.
  - It needs no transport receptors or GTP.
- Wente & Rout 2010, **verified** abstract: https://pubmed.ncbi.nlm.nih.gov/20630994/
  - The NPC must allow free diffusion of water, sugars and ions.
- Pontarin 2008 and Leeds 1985 (Q2, **verified** abstracts) are indirect support.
  - RNR is cytosolic, or absent from nuclear preparations, yet nuclear dNTP pools exist.
  - So deoxynucleotides must reach the nucleus. Pontarin proposes diffusion.
- **No source names NDP, dNDP or UMP explicitly.** PubMed "nucleotide AND nuclear pore complex
  AND passive diffusion AND (ATP OR nucleotides)" gave 8 hits; titles were screened, none
  measures nucleotide permeation.
- Free exchange of these <1 kDa anions is an extrapolation from the general rule. **inferred**

**Bottom line:** the general rule that small molecules equilibrate passively is supported by a
primary yeast in vivo source. Nucleotide-specific measurements were not found.

## Q5 History: nuclear RNR in older models and Yeast-GEM

Local copies only; nothing downloaded. Scan: `scripts/scan_models.py` → `model_scan.txt`.

| Model (local file, SHA-256 prefix) | Nuclear RNR reactions | GPR form | Nuclear Trx supply | Note |
| --- | --- | --- | --- | --- |
| iYali26 source `model/source/iyali26.xml` (5c8c199e) and iYli21 `model/source/iyli21.xml` (6974b758) | R707 (ADP), R709 (CDP), R711 (GDP), R1914 (UDP) | OR of 3 subunits; R1914 "(F39330g or (B05072g and F02309g))" | none | R1914 has H2O in the **cytosol** while every other species is nuclear |
| iYali4 v4.1.2 (`…/eciYali5-GEM/model/iYali5-GEM/sourceModels/iYali4.xml`, c0b54165) | R_y000975/977/979, R_y200021 | OR (YALI0B03630g, D03333g, F31735g) | none; also R_y001037 nuclear Trx peroxidase | R_y200021 has the same H2O[cytoplasm] error |
| iYali5 (`…/iYali5.xml`, b73e99d8; eciYali5-GEM HEAD 816db7f, version 1.0.0) | R_975/977/979, R_iYL0606_1 | OR | none | no nuclear peroxidase |
| Yeast-GEM v9.0.2 (`…/yeast_gem/Yeast-GEM.xml`, 9fd2c572) | r_0975/0977/0979 (no UDP) | "(RNR1 or RNR4 or RNR3 or RNR2)" | only r_4612 "TRX1 disulphide transport" cytoplasm↔nucleus | also r_1037 nuclear Trx peroxidase (DOT5) |
| iMM904 (`…/imm904/iMM904.xml`, a4bcff4f) | RNDR1n–RNDR4n | AND complexes including TRX1/TRX2 | none | THIORDXni nuclear peroxidase |

Each fact in the table is **verified** from the local XML.

- Yeast-GEM r_4612 is gene-less and reversible, bounds −1000/1000.
  - Notes: "Confidence Level: 0" and "added after Gapfilling (PR #185); seed:rxn13406;
    unknown transporter". **verified** (XML notes).
  - It moves only the **disulphide**. No reaction in Yeast-GEM produces reduced TRX1 in the
    nucleus. **verified** (static scan)
- In every model listed, reduced nuclear thioredoxin is only consumed, by irreversible
  reactions, and has no producer.
  - At steady state those reactions must therefore carry zero flux, so nuclear RNR is blocked
    in all of them. **provable from stoichiometry**; no FBA/FVA run.
- Lineage: the iYali line (iYali4 → iYli21 → iYali26) inherited the nuclear RNR copies and the
  OR-of-subunits GPR from the Yeast-GEM family. **inferred** from identical structure, names
  and the y-prefixed IDs.
  - The nuclear UDP copy (R_y200021 → R1914) is a Yarrowia-model addition, already carrying the
    H2O compartment error in iYali4.
- Out of scope, only flagged:
  - R702–R705 "ribonucleoside-triphosphate reductase" carry GPR "YALI1F02309g" (thioredoxin)
    alone, matching Yeast-GEM r_0970–0973 (GPR TRX2).
  - The Y. lipolytica RNR entries describe the diphosphate reaction (Rhea:23252). **verified**
- The GitHub `yeast-GEM.yml` was **not opened**; the local v9.0.2 SBML was sufficient.

## What the evidence implies for options A, B and C (inferred; no decision)

- **A (gene-less Trx diffusion):**
  - Thioredoxin is small enough to diffuse (inferred from Timney), and Trx1/Trx2 have
    high-throughput nuclear GFP calls.
  - But A would switch on catalysis in the nucleus that no yeast source supports. Rnr1 is
    cytoplasmic, and nuclear R2 is the inactive, sequestered pool.
- **B (leave the nuclear copies blocked and document them):**
  - This fits the S. cerevisiae localization literature and the mammalian data (Pontarin;
    Leeds). Those data imply nuclear dNTP supply by diffusion from a cytosolic RNR.
  - It also matches every ancestor model, where the copies are already blocked.
- **C (the evidence suggests two refinements):**
  1. Treat the nuclear copies as having no nuclear catalytic role and rely on nucleotide
     diffusion, since dNDP/NDP exchange is the supported route. Deleting them would need the
     project's non-existence standard, which Y. lipolytica data cannot meet.
  2. Separately, the RNR GPRs are biologically inconsistent: OR across α and β subunits. So
     are R702–R705 and the R1914 H2O compartment. These are curation candidates outside UM11.

## Files

- `model_inventory.txt`: every reaction in the UM7 candidate that touches these genes or
  thioredoxin.
- `model_scan.txt`: RNR and thioredoxin-transport reactions in six models; `model_scan_bounds.txt`: the nuclear Trx consumers are all irreversible (lb 0).
- `blast/`: queries, forward/reverse/self tables, `rbh_summary.tsv`, `residue_check.txt`,
  thioredoxin-family and regulator searches.
- `sources/`: UniProt `.txt`/JSON, NCBI FASTA/esummary, SGD go_details JSON, PubMed XML plus
  extracted abstracts, Timney full-text XML, Mauersberger PMC stub.
- `scripts/`: `model_inventory.py`, `scan_models.py`, `residue_check.py`, `rbh_summary.py`.
- Not written by this evidence task (they appeared in the folder from another session during
  the work, and none was read or used here): `stage1_um11.*`, `blocked_sets.json`,
  `rnr_controls.log`, `ko/`, and `scripts/{stage1_um11,um11_apply,blocked_sets,rnr_controls,ko_compare_variant}.py`.
