# R153 argininosuccinate synthase: is `YALI1D17462g` the right GPR? (2026-10-19)

Read-only review. No model, curation data or code was changed; no build, no LP, no commit. Scope and inputs
are in [TASK.md](TASK.md). Model examined: `E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_20261019.xml`
(`51c6fabd…3fa`, branch `claude/clever-kare-7d2e1a`, now committed as `dc471d7`).

## 0. Answer

- **`YALI1D17462g` is not an argininosuccinate synthase (ASS).** Verified online today: NCBI Gene 2910235,
  RefSeq XP_502815.3, KEGG yli:2910235, UniProt Q6C947 / UniParc UPI00003BC373. It is a 384-aa protein with
  an A1 aspartic-peptidase domain, both A1 active-site motifs and a signal peptide. It has none of the ASS
  family signatures.
- **The ASS is one 413-aa protein in the databases**, with the same sequence in 11 named Y. lipolytica strains.
  These include PO1f, the strain the screens used (recorded from the R190 review). KEGG's yli:90949668 /
  YALI2_D00589g is that protein, but in the **DSM 3286** genome.
  - NCBI has **discontinued** that Gene record.
  - Its RefSeq chromosome (NC_089381.1) is **suppressed**.
- **No W29 (YALI1) gene corresponds to it.** The current W29 RefSeq (NC_090773.1) and the original GenBank
  annotation (CP017556.1) carry no feature from 3,604,630 to 3,609,948. That span contains the ASS homology
  region found on 2026-09-11, where the hit is split across two reading frames (recorded).
  - The model has the ASS under no ID.
  - None of the three screens (Cas9, Cas12a, Patterson Tn-seq) includes it. **So ASS was never assayed as a
    gene in these screens**, and it has no experimental label here.
- **How the gene got onto R153:**
  1. iYli21 (imported 2026-03-19) put `YALI1D17462g` on the duplicate reaction R2176 and the placeholder
     `YALIUNK2` on R153.
  2. On 2026-09-15 the user explicitly assigned `YALI1D17462g` to R153, with the conflicting peptidase
     evidence already recorded, and removed `YALIUNK2`.
  3. R2176 was then merged into R153.

  This is **not** a YALI0/YALI1 number collision: no `YALI0D17462g` exists in NCBI.
- **Proposal (§7), awaiting approval:** make R153 gene-less through a new guarded post-merge step, with notes
  naming the unannotated W29 locus. `YALI1D17462g` stays in the model unlinked.

## 1. What `YALI1D17462g` is (verified 2026-10-19)

| Field | Value | Source |
| --- | --- | --- |
| Systematic IDs | `YALI1_D17462g` (W29/CLIB89); alias `YALI2_D00285g`; the same protein in CLIB122 is `YALI0_D14300g` | [NCBI Gene 2910235](https://www.ncbi.nlm.nih.gov/gene/2910235), [UniParc UPI00003BC373](https://rest.uniprot.org/uniparc/UPI00003BC373.json) |
| Verified name | none verified. NCBI: "uncharacterized protein"; KEGG: "(RefSeq) Candidapepsin"; UniProt Q6C947: "Pepsin A", EC 3.4.23.1 (automatic, PROSITE-ProRule PRU10094); Patterson S3: "similar to … YPS3 … (yapsin 3)" | [KEGG yli:2910235](https://www.kegg.jp/entry/yli:2910235), [UniProt Q6C947](https://rest.uniprot.org/uniprotkb/Q6C947.json) |
| Location | NC_090773.1:1,746,262–1,747,416 (+), one exon, 1,155 nt | NCBI Gene |
| Protein | XP_502815.3 = AOW04049.1, **384 aa**. Sequence SHA-256 `bf9ffec2…`, the same as the target SHA in the model's R153 notes | [XP_502815.3](https://www.ncbi.nlm.nih.gov/protein/XP_502815.3) |
| Domains (UniParc features) | Pfam PF00026 Peptidase A1, 61–368; PROSITE PS51767, 61–366; **PS00141 aspartic-peptidase active site, 76–87 and 259–270**; CDD cd05474 secreted aspartic endopeptidase, 59–368; PRINTS PR00792; PANTHER PTHR47966 | UniParc UPI00003BC373 |
| Other features | signal peptide 1–17; active sites D79, D262 (Q6C947, automatic) | UniProt Q6C947 |
| ASS signatures | **none** (no PF00764/PF20979, TIGR00032, cd01999, MF_00005, PS00564/PS00565) | UniParc UPI00003BC373 |
| UniProt status | W29 accession A0A1D8NEI7 is **inactive** ("Not part of a reference proteome"). The CLIB122 twin Q6C947 is active, TrEMBL v113 | UniProt REST |
| KEGG orthology | no KO assigned. Pfam motifs: Asp, TAXi_N, TAXi_C | KEGG yli:2910235 |
| Evidence level | automatic database annotation and sequence family only. No experimental characterization found (2026-09-11 review; not re-searched today) | — |

The model's own gene annotations for `YALI1D17462g` are ncbigene 2910235, refseq XP_502815.3, uniprot
A0A1D8NEI7 and `kegg.genes/yli:2910235`. They point at the peptidase records above; they are not a
cross-reference to some other protein. Cause (b) of the user's rule 6a ("the local ID points to a different
entity") does not apply: the ID and the protein agree, and it is the *function* assignment that conflicts.

Recorded from 2026-09-11, not re-run:
- BLAST hits two yeast A1 proteases at 31–32% identity, and no ASS reference at E ≤ 1e-5.
- On the AlphaFold prediction AF-A0A1D8NEI7-F1 (mean pLDDT 85.4), US-align TM-scores are 0.743 to the
  protease 1DPJ, and only 0.267 and 0.277 to the ASS structures 2NZ2 and 1K97.

These are "function candidate based on AlphaFold prediction". They support the family call; they do not prove
catalysis.

## 2. The real ASS and its W29 locus

### 2.1 YALI2_D00589g / yli:90949668 (verified 2026-10-19)

- **[KEGG yli:90949668](https://www.kegg.jp/entry/yli:90949668):** YALI2_D00589g, "(RefSeq) Argininosuccinate
  synthase", K01940, module M00844.
  - It is the only yli gene linked to `ec:6.3.4.5` and to `ko:K01940`.
- **KEGG genome yli = [T01033](https://www.kegg.jp/entry/gn:T01033):** now "Yarrowia lipolytica DSM 3286",
  RefSeq assembly GCF_014490615.1 (NC_089378–NC_089383).
  - So the `YALI2_*` IDs and coordinates in the local `reference/kegg/yli_genes.tsv` are DSM 3286 coordinates,
    not W29.
  - KEGG mixes strains in its cross-references: T01033 still lists TAX:284591 (CLIB122), and yli:2910235 links
    the CLIB122 UniProt entry Q6C947.
- **[NCBI Gene 90949668](https://www.ncbi.nlm.nih.gov/gene/90949668):** "This record was discontinued".
  - It was on NC_089381.1:3,672,383–3,673,665 (complement).
  - Nucleotide records it depends on:
    - NC_089381.1 = DSM 3286 chromosome D, from CP061015.1, created 2024-06-26, **status: suppressed**;
    - the current W29 RefSeq chromosome is NC_090773.1, created 2024-09-10 from CP017556.
- **Proteins:**
  - Its RefSeq protein XP_065950363 is inactive (UniParc).
  - The INSDC protein [QNP98148.1](https://www.ncbi.nlm.nih.gov/protein/QNP98148.1) (DSM 3286) is 413 aa.
- **One sequence across strains.** An NCBI Protein search for "argininosuccinate synthase" in taxid 4952
  returns 12 records:
  - 10 are titled ASS/synthetase: KAJ8053850.1 (NRRL Y-64008), QNP98148.1 (DSM 3286), KAE8172175.1 (PO1f),
    KAB8283996.1 (YlCW001), RMI99126.1 (FKP355), and RDW53048.1, RDW46706.1, RDW39874.1, RDW34309.1, RDW26702.1
    (YB567, YB566, YB420, YB419, YB392);
  - 2 are titled by locus tag: VBB88617.1 (H222) and SEI35023.1 (no strain field; YALIA101 tags);
  - all are 413 aa with **one identical sequence** (SHA-256 prefix `3913181c…`);
  - none is from W29/CLIB89 or CLIB122. The audit's broader "argininosuccinate" search (27 hits) finds W29
    and CLIB122 only as 464-aa lyases.
- **Domains** ([UniParc UPI0008A71D72](https://rest.uniprot.org/uniparc/UPI0008A71D72.json)): full ASS
  architecture.
  - CDD cd01999 (4–398), HAMAP MF_00005 (3–407), TIGR00032 (5–405);
  - Pfam PF00764 (5–168) and PF20979 (178–402);
  - PROSITE PS00564 (9–17) and PS00565 (116–127).
- **UniProtKB status:** both Y. lipolytica members (A0A1H6Q7D7, A0A371C926) are inactive. The third member,
  A0A8J5B522, is *Yarrowia* sp. B02, another species with the identical sequence (KAG5355827). A UniProtKB
  search for EC 6.3.4.5 in taxonomy 4952 returns 0 entries.
- **Name and evidence:** there is no verified native name. It is an ARG1 homolog by sequence, and the only
  evidence is computational.

### 2.2 Is there a YALI1 counterpart?

- **None in the current annotation (verified).**
  - NCBI's W29 RefSeq carries YALI2 tags over as `old_locus_tag`. In the chr 1D end window
    3,595,000–3,629,463 the old tags present are D00576, D00578, D00580–D00582, D00584–D00588 and D00590.
    **D00589 is absent.** The independent audit scanned all 1,155 old tags on chr 1D: D00577, D00579, D00583
    and D00589 are absent from the whole chromosome.
  - YALI2_D00583g is the HGD gene, which the orphan gene-search session also found missing (recorded).
  - A Gene search for `YALI2_D00589g` returns only the discontinued 90949668.
- **The unannotated span (verified).**
  - Neighbors: YALI1_D36035g (3,603,430–3,604,629; old YALI2_D00578g) and YALI1_D36100g (3,609,949–3,611,298;
    old YALI2_D00580g).
  - Between them, in both RefSeq NC_090773.1 and GenBank CP017556.1, there is **no feature**: 3,604,630–3,609,948
    (5,319 bp).
- **The homology in that span (recorded 2026-09-11, not re-run).**
  - Magnan 2016 S5 lists a "putative synthase" at YALI1D:3,606,992–3,608,481 with no systematic ID.
  - TBLASTN of yeast ARG1 against CP017556.1:3,606,000–3,609,500 gives two HSPs, in frames −3 and −1, together
    covering 99.52% of the reference.
  - The W29 stretch is 1,240 nt against PO1f's 1,242-nt CDS, with 1,239 identical bases and internal stops in
    frame.
  - **Audit check (2026-10-19, difflib, not TBLASTN).** In DSM 3286 (CP061015.1), KEGG places the ASS between
    the same two flanking genes, YALI2_D00578g and D00580g.
    - The DSM 3286 intergenic span translates to exactly the 413-aa ASS.
    - The W29 span shows 1 substitution and a net −2 nt of indels inside the CDS window, with in-frame stops.
    - This supports the recorded frame split independently.
  - Whether this is a W29 assembly error, a real variant or a pseudogene is **unresolved**. Inferred: PO1f
    descends from W29 and encodes the intact 413-aa protein, which makes an assembly indel the likelier cause,
    but nothing here tests that.
- **CLIB122 (verified).** NCBI has no CLIB122 ASS Gene record. The only CLIB122 protein matching
  "argininosuccinate" is CAG81519.1 / YALI0_D26367g, which is argH **lyase** (TIGR00838), not ASS.

### 2.3 In the model? (verified, static scan of 51c6fabd)

- No gene product or annotation contains `YALI2_D00589`, `90949668`, `XP_065950363`, `QNP98148`,
  `KAE8172175`, `BKA90DRAFT_137900`, `A0A371C926`, `UPI0008A71D72` or `K01940`.
- The only model gene in the D36000–D36300 window is YALI1D36215g. That is old tag YALI2_D00586g, on R1960
  (fumarylacetoacetase), not ASS.

### 2.4 In the screens? (verified, local files)

| Screen | Gene set | Genes between D36035g and D36100g | ARG1 / YOL058W / "argininosuccinate" text |
| --- | --- | --- | --- |
| Cas9 (Ramesh 2023, raw sheet) | 7,854 genes | none (D36035g and D36100g are present, adjacent) | none |
| Cas12a | 7,795 genes | none | none |
| Patterson 2018 Table S3 (Tn-seq) | 8,712 sheet rows | none (rows 5188/5189 are D36035g, D36100g) | none |

**ASS was not assayed as a gene in any of the three screens**, so it has no experimental label. Patterson
reports insertions per gene, so any intergenic insertions at the locus are not visible in S3 (limit).

For context, the same screens on the flanking arginine-pathway genes:

| Gene (name unverified; homology) | Cas9 | Cas12a | Tn-seq Gen0 | Consensus |
| --- | --- | --- | --- | --- |
| `YALI1D18130g` (ARG3-like OTC; R607) | −3.62 essential | −1.72 non-essential | Essential (3 hits / 987 bp) | positive |
| `YALI1D34865g` (ARG4-like ASL; R152) | −5.09 essential | −2.65 essential | Essential (0 / 1,395 bp) | positive |
| `YALI1D17462g` (A1 peptidase; R153) | −0.73 non-essential | −1.41 non-essential | Non-Essential (13 / 1,155 bp) | not listed |

Inferred: in SD-Leu the steps on either side of ASS are screen-essential, while the gene assigned to ASS is
non-essential in all three screens. That fits the gene not being the ASS. **This is a consistency observation,
not identity evidence, and not a reason for the change** (AGENTS.md §4: recall is never evidence).

## 3. How `YALI1D17462g` got onto R153 (verified from git, curation data and notes)

| Date | Event | Evidence |
| --- | --- | --- |
| 2026-03-19 | First real model enters the repo as `data/iyli21.xml` (commit `05ae6e1`, mcnaughtonadm). The initial `0004d76` `model.xml` was a 95-reaction memote template | `git show` |
| 2026-03-19 → now | Both `model/source/iyli21.xml` (`6974b758…`) and `iyali26.xml` (`5c8c199e…`) have R153 = `YALIUNK2` and R2176 "L-Citrulline:L-aspartate ligase (AMP-forming)" = `YALI1D17462g`, with identical 7-species stoichiometry. R2176 is the only reaction using the gene | static listing (`static_r153_source_iyali26.txt`) |
| 2026-09-11 | Two reviews: the R2176 GPR is "strongly suspected mis-annotation"; ASS identity candidate is PO1f KAE8172175.1; the W29 locus is unresolved | research-root reports, STATE.md |
| 2026-09-15 | User instruction: add `YALI1D17462g` to R153 (`YALIUNK2 or YALI1D17462g`); then "remove YALIUNK2" (single gene, `R153-USER-SINGLE-GPR-20260915`, commit `c31921d`); then merge R2176 into R153 | `r153_gpr_assignment.json`, `r153_r2176_merge.json`, STATE.md, R153 notes |
| 2026-10-06/07 | The curation files are committed (`174f22e`) and moved (`d9763b2`) | `git log --follow` |

- **Labels already in the model.** R153's notes record the status "user_authorized_provisional_assignment_with_conflicting_function_evidence"
  and `gpr_experimental_confirmation: required_not_confirmed`. **The present GPR was an explicit,
  informed user choice.** This review adds new identity facts (§2) but does not overturn a hidden error.
- **Upstream origin.** Why iYli21 linked the peptidase gene to R2176 is **unknown**. iYli21's own
  reconstruction source was not established this session.
- **YALI0/YALI1 collision check (verified).** NCBI Gene and Protein have 0 hits for `YALI0D17462g`,
  `YALI0_D17462g` and `YALI0D17462`. The CLIB122 locus for this protein is `YALI0_D14300g`, a different number.
  **The R765/R766-type collision does not apply here.**
- **Builder order (static).** `main.py` runs `apply_r153_gpr_assignment` and then `merge_r153_r2176`. The
  merge's precondition requires R153's GPR to equal `YALI1D17462g`. So any revision has to be a **new step
  after the merge**. Editing `r153_gpr_assignment.json` or `r153_r2176_merge.json` would break the merge guard
  and the record of the user's 2026-09-15 decision.

## 4. Static network around R153 (verified, candidate 51c6fabd; no LP)

- **Reactions whose GPR contains `YALI1D17462g`:** only **R153**: argininosuccinate synthase, [-1000, 1000],
  EC 6.3.4.5, `m141 + m267 + m313 → m86 + m203 + m312` (ATP + L-aspartate + L-citrulline → AMP + PPi +
  argininosuccinate), all in C_cy.
- **Argininosuccinate (m312):** exists only in C_cy, and only **R152** (argininosuccinate lyase, `YALI1D34865g`,
  [-1000, 1000]) and **R153** touch it. No transporter, no other producer.
  - So **the only alternative way to form argininosuccinate is R152 in reverse** (arginine + fumarate →
    argininosuccinate). That consumes arginine and cannot replace R153 as a net arginine source.
- **Cytosolic L-arginine (m310):** the only reactions that can produce it are R152 (forward) and R867
  (transport from C_ex, reversible, OR of 9 transporter genes).
  - Consumers: R151 arginase, R154 Arg-tRNA, R1018 mitochondrial import, and the biomass/amino-acid
    pseudoreactions.
  - **External arginine depends on the runtime medium.** In the stored SBML, the arginine exchange R1202 is
    [0, 1000], export only. Uptake exists only because the runtime SD-Leu medium sets R1202 = 0.02586
    (`state/media/sd_leu.csv`, `ed176d26…`; its comment says "uptake scaled to glucose … (inferred)").
- **Citrulline (m313):** only R607 (OTC, irreversible) and R153 (reversible) touch it.
- **So de novo arginine in the model needs R607 → R153 → R152.** For the gene `YALI1D17462g` specifically, R153
  is its only reaction, so its single KO removes de novo arginine. The OTC and ASL KOs do the same at their
  own steps.

Full listing: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r153_ass_gpr_20261019/static_r153_candidate_51c6fabd.txt`.

**Recorded KO results** (earlier runs, not re-run): `YALI1D17462g` KO/WT was
- 0.110 on the 2026-09-15 R1931-forward model, the same as the OTC and ASL KOs;
- 0.159 on the 2026-10-12 C2 benchmark model (`w29_13c_consistency_20261012/c2_benchmark_ratio_changes.csv`).

No KO of this gene is recorded on 51c6fabd.
- Inferred, not computed: it would equal the arginine-uptake-limited value seen for the CPA KOs (0.1439),
  because both KOs leave uptake as the only arginine source.
- All of these ratios depend on the runtime SD-Leu arginine bound, which is itself inferred.

## 5. Assessment, with the conditions behind the verdict

Claim: **"R153 = `YALI1D17462g` is not supported, and the identity evidence contradicts it."** This needs
all of the following:

1. **The gene's protein lacks ASS family signatures.** Verified (§1: UniParc features).
2. **The protein belongs to a different, defined family with its catalytic motifs.** Verified (§1:
   PF00026, PS00141 at both motifs, cd05474).
3. **The species' ASS is a different sequence.**
   - Verified: database records exist for a 413-aa Y. lipolytica protein with full ASS architecture (§2.1).
   - Computational/inferred: that this protein is the functional enzyme in vivo.
   - Recorded (2026-09-11) and supported by the audit's alignment: its W29 locus has homology.
   - Verified: that locus has no gene model.
4. **No experimental evidence shows ASS activity for `YALI1D17462g`.** Recorded: the 2026-09-11 search
   found none; not re-searched today.
5. **The assignment did not rest on evidence this review has not seen.** Verified for 2026-09-15: the
   basis is an "explicit user instruction" (notes). The iYli21 basis for R2176 is unknown.

Verdict tiers:
- Conditions 1–3 and 5 (2026-09-15 part) are **verified**; 4 is **recorded**.
- "`YALI1D17462g` does not catalyse argininosuccinate synthesis in vivo" stays **inferred**: no experiment
  excludes a side activity.
- AGENTS.md §5 asks for evidence that a GPR member catalyses the step; there is none.

What this review does not establish:
- the native W29 ASS coding sequence;
- whether W29 or PO1f express a full-length ASS (untested);
- ASS localization;
- any growth or essentiality effect of a GPR change.

## 6. What the numbers mean for scoring (reported separately, not a reason)

Under `model/expected/benchmark.md`, `YALI1D17462g` is a model gene not on the positive list, so it is in the
user-defined negative class.

- **Today:** its 15% call depends on the arginine-uptake-limited KO ratio (0.110 → FP on the 09-15 model;
  0.159 → TN on the 10-12 model). On 51c6fabd it is not yet computed.
- **If the GPR is removed:** its KO touches no reaction, so it becomes TN at every cutoff (gene kept), or
  leaves the scored set (gene removed).
- **The true ASS:** has no screen row, so no positive-class count changes.
- **This shift is a by-product.** The proposal rests only on §1–§3, and it must not be counted as a
  validation gain.

## 7. Plan for approval

Each option is separate. **Recommended: A.**

> **Status (2026-10-19):** the user approved option A, with A2(i): the gene kept unlinked. It is implemented as
> `R153-GPR-UNASSIGNED-20261019` under the flag `--remove-r153-peptidase-gpr` (the name differs from the
> proposal below). New candidate `…_r153_gpr_v2_20261019.xml` (`d941d2cd…`); see
> [../r153_gpr_revision_20261019/REPORT.md](../r153_gpr_revision_20261019/REPORT.md).

**A. Make R153 gene-less, recording the unannotated W29 locus (recommended).**
- **What changes:**
  - R153 GPR `YALI1D17462g` → empty.
  - Bounds, stoichiometry, EC and compartment are unchanged.
  - Gene `YALI1D17462g` stays in the model with no reactions (as with `YALI1E11370g` in the selection fix);
    see sub-choice A2.
- **How:**
  - A new curated file `model/curation/gpr/r153_gpr_revision.json`, with the before state, the evidence of §1–§3
    with sources, and the authorization line.
  - A guarded `apply_r153_gpr_revision` step that runs **after `merge_r153_r2176`**, behind an opt-in flag
    (proposed name `--revise-r153-gpr-unannotated-ass`).
- **The guard refuses unless:**
  - R153's GPR is `YALI1D17462g`;
  - bounds are [-1000, 1000];
  - the stoichiometry and species match the merge spec;
  - the note `curated_gpr_assignment = R153-USER-SINGLE-GPR-20260915` is present;
  - the gene's annotations are 2910235 / XP_502815.3 / A0A1D8NEI7;
  - R153 is the gene's only reaction.
- **Notes:**
  - New keys record the W29 gap (NC_090773.1:3,604,630–3,609,948), the identity candidate (PO1f
    KAE8172175.1 = DSM 3286 QNP98148.1, 413 aa; not a W29 ID), "ASS not assayed in any of the three screens",
    and "no gene KO disables R153".
  - The five `gpr_evidence_*` / `gpr_experimental_confirmation` keys are overwritten. Their old values are kept
    in the curation file and in an `r153_gpr_revision_supersedes` note.
  - The 2026-09-15 curation bytes are untouched.
- **Tests:** a focused `platform/tests/test_r153_gpr_revision.py` covering:
  - guard refusals;
  - the applied state;
  - idempotence;
  - flag off → output byte-identical to 51c6fabd.
- **Merge re-run caveat (from the audit).** `merge_r153_r2176`'s already-merged path (exercised by
  `test_r153_merge.py`) compares R153's GPR and notes with `before.R153`. Run on the revised candidate, it
  would raise. Either that path learns to accept the recorded revision, or tests stop applying it to revised
  candidates. Settle this in the implementation, not by editing the merge JSON.
- **Build:**
  - Command: the 51c6fabd argv plus the new flag.
  - Output, a new file: `model/candidates/E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_r153_gpr_20261019.xml`.
- **Acceptance:**
  - The whole-model diff against 51c6fabd shows only R153's GPR and notes, plus the gene's association.
  - An independent read-only audit.
  - REPORT and STATE.md.
- **LPs (separate authorization, optional):**
  1. KO `YALI1D17462g` on 51c6fabd, to record today's value.
  2. WT on the new candidate.

  Provable without LPs: WT is unchanged, because a GPR does not enter the WT LP; the new gene KO equals WT,
  because the gene has no reactions.
- **Dependencies:**
  - The baseline flags live on `claude/clever-kare-7d2e1a` (`dc471d7`), not merged.
  - The D02625g proposal (`claude/cool-nightingale-a2f471`) also extends `main.py` / `cli.py`.
  - Order the work after clever-kare lands, or coordinate with that session.
- **Sub-choice A2:**
  - (i) keep `YALI1D17462g` unlinked (recommended: it matches the E11370g precedent and keeps the gene
    count at 1,073); or
  - (ii) remove it.
  - Either way the scoring effect of §6 is reported separately.

**B. Keep the 2026-09-15 assignment (status quo).** No GPR change.
- Optionally add notes saying that KEGG's ASS is a DSM 3286 gene, discontinued in RefSeq, and that the W29
  locus is unannotated.
- Cost: the model keeps a GPR contradicted by identity evidence. The peptidase gene's KO then stands in for
  "no de novo arginine".

**C. Assign a cross-strain ID for the ASS ortholog** (for example PO1f `BKA90DRAFT_137900` /
KAE8172175.1). **Not recommended now.**
- There is no W29 gene model.
- The gene would have no screen row (negative class by definition).
- The W29 CDS is unresolved.
- Reconsider once the locus has a W29 systematic ID or sequencing/RNA-seq evidence.

**Not recommended:** restoring `YALIUNK2` (the user removed it explicitly, and it has no identity), or
inventing a YALI1-style ID.

**D. Records (optional, no model effect):**
- D1: paste the STATE.md entry in §10.
- D2: separately, a read-only check of whether public W29/PO1f RNA-seq covers the locus as one ORF. This would
  test the frame split; it is not authorized here.

## 8. Code logic, method notes and verification level

- **Online checks:** NCBI E-utilities (Gene, Nuccore feature tables, Protein), UniProt and UniParc REST, and
  KEGG REST; exact URLs are in [commands.txt](commands.txt). Feature-table coordinates with `seq_start` are
  relative. They were converted and cross-checked against the Gene records of the two flanking genes.
- **Static model scan:** `scripts/static_r153.py` parses SBML text. It prints every reaction whose GPR
  contains the gene, and every reaction that touches given species bases, with bounds from the `fbc`
  parameters. It solves nothing.
  - A possible improvement: cross-check it with a cobrapy load (`gene.reactions`, `metabolite.reactions`),
    to rule out regex misses on unusual SBML layouts. Test it by comparing both listings on 51c6fabd.
- **Sequence identity:** the 12 ASS records were compared by exact string and SHA-256, not by alignment.
  "Identical" is an exact claim here.
- **Verification levels:**

| Level | What |
| --- | --- |
| Provable | WT unchanged by a GPR edit; KO of a reaction-less gene equals WT |
| Statically verified (today) | §2.3, §2.4, §3 git and source facts, §4 |
| Verified online (today) | §1, §2.1, §2.2 annotation gaps |
| Independently audited (2026-10-19) | 11 key claims (§1–§4 and the §3 builder order): 11 confirmed, 0 refuted. The audit's 9 wording/caveat points are applied: the 12-record titles, *Yarrowia* sp. B02, "assayed as a gene in these screens", the §5 condition-3 tier, the R153 reversibility, the medium-dependent uptake, the whole-chromosome tag scan, the merge re-run caveat, and KEGG's strain mixing |
| Recorded (earlier sessions, not re-run) | BLAST/TBLASTN/US-align/AlphaFold numbers (2026-09-11); KO ratios in §4 |
| Untested | W29 ASS CDS; expression; localization; growth effect on 51c6fabd |

## 9. Sources

| Source | Read | Tier |
| --- | --- | --- |
| [NCBI Gene 2910235](https://www.ncbi.nlm.nih.gov/gene/2910235), [90949668](https://www.ncbi.nlm.nih.gov/gene/90949668), 90949658, 90949660 | E-utilities esummary/efetch | verified |
| [NC_090773.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_090773.1), [CP017556.1](https://www.ncbi.nlm.nih.gov/nuccore/CP017556.1), [NC_089381.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_089381.1), [CP061015.1](https://www.ncbi.nlm.nih.gov/nuccore/CP061015.1) | esummary; feature tables for the chr D end | verified |
| [XP_502815.3](https://www.ncbi.nlm.nih.gov/protein/XP_502815.3), [QNP98148.1](https://www.ncbi.nlm.nih.gov/protein/QNP98148.1), [KAE8172175.1](https://www.ncbi.nlm.nih.gov/protein/KAE8172175.1), the 12-record ASS set, [CAG81519.1](https://www.ncbi.nlm.nih.gov/protein/CAG81519.1) | esummary; FASTA; GenPept for CAG81519 | verified |
| [UniProt Q6C947](https://rest.uniprot.org/uniprotkb/Q6C947.json), A0A1D8NEI7, A0A371C926; [UniParc UPI00003BC373](https://rest.uniprot.org/uniparc/UPI00003BC373.json), [UPI0008A71D72](https://rest.uniprot.org/uniparc/UPI0008A71D72.json) | REST JSON | verified |
| [KEGG yli:2910235](https://www.kegg.jp/entry/yli:2910235), [yli:90949668](https://www.kegg.jp/entry/yli:90949668), [T01033](https://www.kegg.jp/entry/gn:T01033); `link` ec:6.3.4.5 and ko:K01940 | REST | verified |
| Patterson 2018 Table S3 (`mmc12.xlsx`, `7792f80b…38a9`) | local copy in the fervent-knuth session scratchpad | verified as file |
| Ramesh 2023 Cas9/Cas12a raw (`ba1eca8f…`), consensus (`1e887f5a…`), MOESM10 (`50d9ce18…`) | research root, read-only | verified as files |
| 2026-09-11 reports (`r153_r2176_gpr_review_20260911`, `argininosuccinate_identity_20260911`) | read | recorded |

Raw responses and listings: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r153_ass_gpr_20261019/` (23 files plus
`SHA256SUMS`, 132 KB).

## 10. Proposed STATE.md entry (not written)

STATE.md was not edited, because several worktrees hold diverging uncommitted heads. Ready to paste:

> ## 2026-10-19: R153 argininosuccinate synthase GPR identity review, investigation only
>
> Delegated read-only review; no model, data or code change, no build, no LP. Model examined:
> E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_20261019.xml (51c6fabd…, clever-kare `dc471d7`).
> - Verified online today: R153's only gene YALI1D17462g (Gene 2910235, XP_502815.3, 384 aa) is an
>   A1 aspartic peptidase (PF00026, PS00141 at both motifs, signal peptide) with no ASS signature.
> - The databases' ASS is one 413-aa sequence in 11 strains, including PO1f KAE8172175.1.
> - KEGG's ASS yli:90949668 / YALI2_D00589g is a DSM 3286 gene. Its NCBI Gene record is discontinued and its
>   RefSeq chromosome NC_089381.1 suppressed.
> - W29 has no YALI1 gene there: NC_090773.1:3,604,630–3,609,948 is unannotated. The ASS homology there is
>   frame-split (recorded 2026-09-11).
> - Not in the model under any ID; not in the Cas9, Cas12a or Tn-seq gene sets, so ASS was never assayed as
>   a gene in these screens.
> - An independent read-only audit re-checked 11 claims and refuted none; its wording fixes are applied.
> - Provenance: iYli21 R2176 → user assignment to R153 on 2026-09-15 (conflicting evidence recorded then) →
>   merge. Not a YALI0/YALI1 collision.
> - Proposal awaiting the user: option A, a gene-less R153 through a guarded post-merge step and opt-in flag.
>   Report: model/reports/r153_ass_gpr_20261019/.
