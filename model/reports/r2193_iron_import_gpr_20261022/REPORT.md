# R2193 iron import: which genes, investigation and plan (2026-10-22)

## Scope

The user, 2026-10-22:

> "Yes: investigate the R2193 iron-import GPR (YALI1F07747g looks like a bHLH regulator). Find the
> real Y. lipolytica iron-uptake genes, investigation and plan only, no model change."

- **Not changed:** no model, curation data or builder code.
- **Computation:** no LPs.
  - One local BLAST run: 12 *S. cerevisiae* queries against the *Y. lipolytica* CLIB122 proteome,
    plus the reverse search of their hits.
  - Read-only database and literature look-ups.
- **Large and raw files:** in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r2193_iron_import_gpr_20261022/`:
  - `evidence/kegg`, `evidence/uniprot`, `evidence/literature` and `evidence/ncbi_cds`;
  - `blast/proteomes` (with `SHA256SUMS`) and `blast/out`.
- **Earlier work:**
  - [r2193_pair1_20260727/experiment_result.md](../r2193_pair1_20260727/experiment_result.md)
    (main line, July) tested the D08 pair below on a different model. It did not consider the D09
    pair.
  - The heme/iron candidate `7f836c3a` ([report](../heme_iron_biomass_20261022/REPORT.md)) made
    R2193 required for growth, which is why this matters now.

Evidence tags:
- **verified** = opened or run this session, with output kept;
- **inferred** = reasoned;
- **recalled** = not opened.

## 1. The current gene is a transcription factor, not a transporter (verified)

- **The model's mapping is correct.** YALI1F07747g carries NCBI Gene 2908631, which is what the
  NCBI W29 feature table gives for `YALI1_F07747g` (XP_505019.3, 716 aa; old tag YALI0F05126g). So
  the problem is the gene-to-reaction assignment, not the ID mapping (AGENTS/CLAUDE §6a).
- **Domains.** UniProt Q6C2U3 carries only HLH/bHLH DNA-binding domains (PF00010; IPR011598,
  IPR036638) and no transport domain.
- **KEGG** names the gene "Phosphorus acquisition-controlling protein" and assigns no KO.
- **BLAST.** Its reciprocal best hit in *S. cerevisiae* is PHO4, the bHLH phosphate-regulon
  activator. The hit covers only the bHLH region (18 % of PHO4, E = 1.7e-7).
- **A *Y. lipolytica* paper** ([Pomraning et al. 2016, PMC4766638](https://europepmc.org/article/PMC/PMC4766638))
  lists YALI0F05126g as one of nine predicted bHLH transcription factors, down-regulated under
  nitrogen limitation. That is a prediction plus expression data, not a tested function.
- **Origin.** The gene came with the reaction from the parent iYli21 gene set
  (`reference/locus_map/iyli21_genes_vs_S2.csv`). The Chalmers model eciYali5 has the same reaction
  ("iron (II) transport") with no gene. How it was first assigned is unknown.

## 2. *Y. lipolytica* iron-uptake genes

Sources for the table:
- **Orthology:** reciprocal best BLAST hits against reviewed *S. cerevisiae* proteins (BLAST+ 2.17.0,
  E ≤ 1e-5), KEGG KO and Pfam. These are database and sequence identities, not experiments.
- **Experimental list:** the consensus essential list used by the benchmark.
- **Protein detection:** the Debiaggi 2025 proteomics table shipped with eciYali5 (24 samples: 4
  strains × 2 conditions labelled e/n). The units were not identified, so only comparisons within
  the table are made.

Full table: [candidate_genes.tsv](candidate_genes.tsv).

| Gene (W29 / CLIB122) | Ortholog of | Essential list | Proteomics | Same-species functional evidence |
| --- | --- | --- | --- | --- |
| **D09 pair** | | | | |
| YALI1_D09277g / YALI0D07282g (Q6C9Y8) | FET3, ferroxidase (45.8 % id) | **positive** | 24/24, the most abundant of the four (92nd percentile) | none found |
| YALI1_D09308g / YALI0D07304g (Q6C9Y7) | FTR1, iron permease (49.4 % id) | **positive** | 15/24, low | none found |
| **D08 pair** | | | | |
| YALI1_D08684g / YALI0D06754g (Q6CA12) | FET3 (50.4 % id) | not listed | 23/24 | not manipulated |
| YALI1_D08564g / YALI0D06688g (Q6CA15) | FTR1 (43.8 % id) | not listed | 24/24 | overexpressed (below) |
| **Other iron genes (not R2193)** | | | | |
| YALI1_A05031g (Q6CHU4) | FTH1, vacuolar iron export (not FTR1) | not listed | — | none |
| YALI1_B28215g (Q6CDS3) | MRS3/4, mitochondrial iron carrier | not listed | — | none |
| YALI1_C05880g, YALI1_D35356g | SMF1, SMF2 (Nramp metal transporters) | not listed | — | none |

- **Pair layout.** The D09 genes lie about 3 kb apart in opposite orientation (chr 1D, 925.8–932.1
  kb), the usual fungal FET3–FTR1 arrangement. The D08 genes lie about 10 kb apart (855.2–868.4 kb).
  Gene positions are verified; that such pairs are typical of fungi is recalled.
- **Low-affinity Fe²⁺ uptake.** There is no FET4 ortholog: no BLAST hit and no KEGG K19792 gene.
  SMF-type transporters exist. In *S. cerevisiae* they mainly carry Mn (recalled).
- **Siderophores.** Siderophore-iron transporters are expanded in *Y. lipolytica* (Diffels et al.
  2006, PMID 17011109, abstract verified). They are not relevant to iron salts in YNB-type media.
- **Same-species experiment (verified).**
  - [Ye et al. 2026, Synth Syst Biotechnol](https://doi.org/10.1016/j.synbio.2026.02.004)
    ([PMC12969333](https://europepmc.org/article/PMC/PMC12969333)) co-expressed *Vitreoscilla*
    hemoglobin with "a high-affinity iron transporter (Ftr1)" and reported raised intracellular
    iron.
  - Their primers (supplement mmc1) match the start (19 nt) and end (25 nt) of the YALI0_D06688g
    coding sequence (XM_502497.1), and nothing in the other two FTR1-like genes. So the cloned gene
    is **YALI1_D08564g**.
  - Limits: no Ftr1-only control, no knockout or complementation, no uptake assay.
- **AND within a pair** (*S. cerevisiae*, verified abstract). In
  [Stearman et al. 1996, Science, PMID 8599111](https://europepmc.org/article/MED/8599111),
  Fet3 and Ftr1 together mediate high-affinity uptake: Fet3 is needed to move Ftr1 to the plasma
  membrane, and Ftr1 is needed to load Fet3 with copper. That the *Y. lipolytica* pairs work the same
  way is inferred.

## 3. Conflicting evidence on the two pairs

- **For treating both pairs as able to import iron:**
  - both are full FET3/FTR1 orthologs;
  - all four proteins are detected;
  - the D08 permease is the only one with same-species functional support.
- **Against D08 standing in for D09:**
  - both D09 genes are on the experimental essential list (CRISPR screens);
  - so under the screen conditions the D08 pair (and SMF) did not compensate.

  Why is unknown: expression level, localization, regulation or capacity (inferred).
- **What a GPR can express.** A Boolean rule records which proteins can catalyse a step. It cannot
  record that one isozyme is insufficient. Choosing the D09 pair alone *because* it matches the
  essential list would be fitting the GPR to the benchmark (AGENTS.md §2, §6). The list is
  development data, not independent evidence.

## 4. Plan (needs a separate authorization; nothing done)

Options for R2193 (Fe²⁺ e → c):

| Option | GPR | Basis | Expected knockout effect on `7f836c3a` (inferred, not run) |
| --- | --- | --- | --- |
| **A (recommended)** | `(YALI1D09277g and YALI1D09308g) or (YALI1D08684g and YALI1D08564g)` | orthology for both pairs, AND from the yeast complex, protein detection, the D08 permease experiment | YALI1F07747g no longer lethal (its false positive goes). No single new gene is lethal. The two D09 genes enter the model as false negatives, with the reason recorded. |
| B | `YALI1D09277g and YALI1D09308g` | orthology and abundance, plus the essential list | Both D09 genes lethal (true positives). Not recommended: it leaves out a pair with equal orthology and same-species evidence, and the reason for leaving it out would be the benchmark. |
| C | `YALI1D08684g and YALI1D08564g` (the July "full Pair 1") | the D08 permease experiment | Both D08 genes lethal (false positives). Leaves out the more abundant D09 pair without a reason. |
| D | no gene (as eciYali5) | no claim | YALI1F07747g's false positive goes; iron import is not under gene control. |

For every option:
- YALI1F07747g leaves R2193 and stays in the model unlinked, as R153 did with YALI1D17462g.
- R2193's stoichiometry and bounds are unchanged. Its chemistry stays blocked by the iron species
  with no formula and charge 0, a separate finding.

Steps for option A (skill `curation-change`):
1. **Baseline and scope.** Baseline `7f836c3a` (or `cf0825b8`), named by the user; TASK.md.
2. **Data and step.**
   - `model/curation/gpr/r2193_iron_import_gpr.json`: the before GPR (`YALI1F07747g`), the after
     GPR, the four new genes with their NCBI/RefSeq/UniProt/KEGG IDs, and the evidence above.
   - Guarded `apply_r2193_iron_import_gpr` behind `--r2193-iron-gpr`.
   - Focused in-memory test.
3. **Build to a new file.** A control build without the flag must reproduce the baseline. Two seeds.
   The whole-model diff must be only R2193's GPR, the 4 added genes and their links.
4. **LPs (Task A medium + PO1f).**
   - WT, expected unchanged.
   - Knockout comparison of both models, solved from a fresh basis, or with a bound-violation check
     because of the heme near-tolerance issue.
   - 1/5/10/15 % scores, reporting the 4 genes newly in scope separately.
5. Independent audit, REPORT, STATE entry. No commit or push unless asked.

Outside this plan, for a separate decision:
- **R863 (Fe²⁺ c → m, gene-less).** YALI1_B28215g is the only MRS3/4 ortholog. In *S. cerevisiae*,
  cells lacking both MRS3 and MRS4 are viable (recalled, uncited). Making it the sole gene would make
  it lethal in the model. Keep R863 gene-less, or record the candidate in notes only.
- **Fe³⁺ reduction (FRE ferric reductases).** Several FRE2/FRE3-like genes exist (for example Q6CF26
  and Q6CET5, reciprocal best hits of FRE2). YNB supplies Fe³⁺. A mechanistic Fe³⁺ → FRE → FET3/FTR1
  route would need iron species with real formulas and charges first.

## 5. Limits

- Identity of all candidates is orthology plus domains, not *Y. lipolytica* biochemistry. Only the
  D08 permease has a same-species experiment, and that one is indirect.
- The proteomics units and conditions were not checked; detection is only presence and relative rank
  within the table.
- The literature search was bounded: Europe PMC, 4 queries plus locus IDs. Studies that name these
  genes differently may be missed.
- The essential-list status is development data, not independent validation.

## Correction after the independent audit of the change (2026-10-22)

The text above is kept as written. These corrections apply (verified from `sce_vs_yli.tsv`):
- **D08 genes.** YALI1_D08684g and YALI1_D08564g are **not** reciprocal best hits.
  - FET3's best *Y. lipolytica* hit is Q6C9Y8 (D09, 553 bits), ahead of Q6CA12 (536).
  - FTR1's best is Q6C9Y7 (D09, 392 bits), ahead of Q6CA15 (346).
  - The D08 genes are paralogs of the D09 genes, with FET3/FTR1 as their best yeast hits (one-way).
  - The "RBH" label in §2 checked only the *Y. lipolytica* → yeast direction.
- **Same for FRE2.** Only Q6CF26 is FRE2's reciprocal best hit; Q6CET5 is one-way.
- **Ye et al. 2026** co-expressed YALI1_D08564g with *Vitreoscilla* hemoglobin, with no Ftr1-only
  control.
- [candidate_genes.tsv](candidate_genes.tsv) is corrected; its column is now `s_cerevisiae_best_hit`.
- The change report, [r2193_iron_gpr_change_20261022](../r2193_iron_gpr_change_20261022/REPORT.md),
  carries the corrected wording into the model notes (candidate v3).
