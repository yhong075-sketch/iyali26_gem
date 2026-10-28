# Why consensus-essential genes are predicted non-essential (a1237821, 2026-10-27)

Scope and starting state: [TASK.md](TASK.md). Investigation only: no model, GPR, medium, benchmark or curated-data
change; nothing committed. Per-gene table: [fn_causes_a1237821.csv](fn_causes_a1237821.csv) (212 rows).

Tags: **verified** = run or read this session, output named; **inferred** = reasoned from verified outputs, no
dedicated run; **recalled** = from an earlier report, not re-run.

## 1. Summary

- **The FN list reproduces** (verified). On `a1237821…` with the Task A medium + PO1f, 212 consensus-essential
  model genes have KO/WT ≥ 15 %. Counts match the 2026-10-23 build sweep (103 / 212 / 39 / 695, 16 unresolved).
- **Half of the false negatives cannot be fixed by GPR edits** (verified mechanism, inferred consequence).
  104 genes (category d) sit on reactions that are blocked in WT or carry no flux at the optimum. The model has
  no demand for their products: biomass has no NAD(P), FAD, FMN, CoA, folate, pyridoxal phosphate, biotin,
  sphingolipid, GPI anchor, glycoprotein or mitochondrial-protein term, and no organelle acidification.
- **31 genes are hidden by an OR partner** (category a, verified). Closing the gene's reactions outright drops
  growth below 15 %, so the rule, not the network, makes them non-essential. The largest groups:
  - 9 ATP synthase subunits: R171 is "(14-subunit complex) OR ATP20-like YALI1B28025g OR ATP19-like
    YALI1B15733g". Each single accessory subunit keeps R171 open. With R171 closed, KO/WT is 0.1097.
  - 6 aminoacyl-tRNA synthetases: on cytosolic R759, R723 and R673 the cytosolic enzyme is OR'd with its
    mitochondrial-type paralog. With the reaction closed, growth is 0.
- **58 genes have an alternative route** (category b, verified). Five escapes are confirmed by bypass-closed
  controls:
  - R2024 (propanoate:CoA ligase, AMP-forming) is reversible. Run backwards, it rescues the ADK1-like and
    IPP1-like knockouts.
  - R1172 is a gene-less HMG-CoA transporter between mitochondrion and cytosol. With leucine catabolism it
    rescues the ERG13-like and ERG10-like knockouts.
  - R2081 duplicates glutamine synthetase R334.
  - Complex I (13 genes, KO/WT 0.758) is bypassed by NDH2-like R570 and redox shuttles, as already shown on
    2026-10-14 (recalled).
- **13 amino-acid-pathway genes are rescued by the medium** (category c, verified). Their KO/WT is 0.159–0.278,
  just above the cutoff, and they become lethal once the 11 amino-acid and adenine uptakes are closed. These
  calls rest on the concentration-scaled uptake caps, which are assumptions (flag f).
- **6 genes are in the model but have no reaction** (category e).
- **Literature** (Europe PMC, *Y. lipolytica* only): 48 of 56 gene–reaction pairs had no hit. The one pointed
  result is that the R712 OR partner YALI1F02518g has only homology support as a ribose-5-phosphate
  isomerase. Its only genetic phenotype in the hits is in erythritol assimilation.

## 2. Run and scores (verified)

| Item | Value |
| --- | --- |
| Model | `a12378210c7c66720d99fb7bf71d14cb7538ea6a3cf3e08919c85fb3c9c47d82` (newest candidate on `origin/main`; no reference designated), 1,081 native genes |
| Medium / profile | Task A medium `6eca70fd…` + PO1f `35307853…` |
| List | `1e887f5a…`, 1,612 genes; 326 in the model, **1,286 out of scope** (never FN) |
| Solver | Gurobi 13.0.1, tolerances 1e-9, reset per knockout, knocked-flux check 1e-9 |
| WT | 0.31242085277746734 h⁻¹ |
| Screen | 1,082 LPs, 13 s; 1,061 optimal, 20 infeasible (16 in scope: complex III/IV/ATP-synthase energy shortfall, STATE 2026-10-21) |

Negatives are the user-defined class (model genes not on the list). mtDNA genes are out of scope; none of them is
on the list. The list was used during model development, so these are development-reference scores, not
independent validation.

| Cutoff | TP | FN | FP | TN | Unresolved | Recall | Precision | Specificity | MCC |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 % | 85 | 230 | 32 | 702 | 16 | 0.270 | 0.727 | 0.956 | 0.329 |
| 5 % | 91 | 224 | 34 | 700 | 16 | 0.289 | 0.728 | 0.954 | 0.343 |
| 10 % | 95 | 220 | 39 | 695 | 16 | 0.302 | 0.709 | 0.947 | 0.341 |
| **15 %** | **103** | **212** | **39** | **695** | **16** | **0.327** | **0.725** | **0.947** | **0.367** |

The 15 % row equals the `a1237821` sweep in `reports/lipoylation_acylacp_on_2c9cebd0_20261023` (recalled from
STATE 2026-10-23; same counts).

KO/WT of the 212 FN genes at 15 %: 161 at ≥ 0.999, 24 in 0.9–0.999, 14 in 0.5–0.9, 13 in 0.15–0.5.

## 3. Method

Scripts in this folder, all in-memory with SHA-checked inputs:

1. [run_screen.py](run_screen.py) is main's 2c9cebd0 runner with only the model swapped, plus a dense
   mass-residual matrix because the venv has no scipy.
2. [classify_fn.py](classify_fn.py) computes, for each FN gene:
   - the reactions its KO closes, from the screen's GPR propagation;
   - the OR partners that keep the other reactions open;
   - a **reaction-level KO**, with every reaction of the gene closed regardless of GPR;
   - WT FVA at 100 % of the optimum for every touched reaction;
   - the model's blocked set (882 of 2,316 reactions) from `find_blocked_reactions` with the medium as set;
   - pFBA flux changes in the reaction-level KO.
3. [trace_controls.py](trace_controls.py) runs three controls:
   - every FN gene with the 11 amino-acid and adenine uptakes closed, and again with the 10 vitamin uptakes
     closed;
   - a dead-end check on each blocked reaction;
   - an FVA (99 %) bypass trace on 16 genes.
4. [trace_extra.py](trace_extra.py) adds:
   - the full GPR texts and the biomass composition;
   - bypass traces for ENO and PYK;
   - cofactor-gene KOs with vitamins closed;
   - bypass-closed controls for R2024, R1172 and R1252.
5. [build_table.py](build_table.py) merges everything into the CSV. It needs no solver. It overwrites the CSV
   that `classify_fn.py` first wrote.

LPs in total: 1,082 + 5,829 + 74,582 + 50,990. The FVA sweeps dominate. Each run asserts that WT is restored
at the end.

**Category rule** (applied in this order, mechanical):

- **e:** the gene has no reaction.
- **c:** the KO becomes essential (< 15 %) when amino acids and adenine are closed. Flag **f**: the call depends
  on an assumed cap.
- **b:** a bypass-closed control makes the KO lethal (verified bypass).
- **a:** OR partners exist, and the reaction-level KO is < 15 % or has no solution.
- **b:** complex I (R2062).
- **b:** otherwise, if the reaction-level KO is < 0.999 of WT, or a closed reaction carries flux at the
  optimum.
- **d:** otherwise. The closed reactions are blocked in WT ("blocked_in_wt") or carry no flux at the optimum
  ("unused_at_optimum").

The d sub-labels (cofactor, mitochondrial translation, sphingolipid, GPI, glycosylation, V-ATPase) come from
keyword matches on the list annotation. They are triage, not identities.

No FN gene was classified **g**. None was classified as a pure **f**: conditional essentiality appears only as
the c flag.

## 4. Causes (counts from the CSV, verified)

| Cause | Genes | Sub-cause (genes) |
| --- | ---: | --- |
| a. OR partner carries the flux | 31 | ATP synthase R171 (9); cytosolic aaRS R759/R723/R673 (6); iron import R2193 (2); others (14) |
| b. alternative or duplicate route | 58 | complex I (13); growth cost < 0.999 (26); no growth cost (13); verified bypass (5); duplicate reaction (1) |
| c. medium supplies the product (flag f) | 13 | arginine, threonine/homoserine, isoleucine/valine, chorismate-mutase genes, NCE103 |
| d. product not demanded | 104 | blocked in WT (51): cofactor synthesis 12, mito translation 7, glycosylation 6, sphingolipid 6, other 20; unused at optimum (53): V-ATPase 11, GPI 7, sphingolipid 3, cofactor 1, other 31 |
| e. gene without reaction | 6 | RAM1-, RAM2-, BET4-, CDC43-, DIM1-, PET112-like |
| **Total** | **212** | |

### 4.1 Category a: OR rules (verified mechanism; biology not verified)

| Reaction | FN genes | Partner keeping it open | Reaction-level KO/WT |
| --- | --- | --- | ---: |
| R171 ATP synthase | YALI1B05364g, YALI1B08885g, YALI1D14639g, YALI1D15734g, YALI1D27826g, YALI1F04110g, YALI1F04433g, YALI1F07384g, YALI1F26967g | YALI1B28025g alone, YALI1B15733g alone (UniProt B5FVB8 "ATP20" and B5FVB3 "ATP19", read from the local UniProt GO table in `reference/external_models/…/yali_genes_GOterms.tsv`; not checked online) | 0.1097 |
| R759 Thr-tRNA (cytosol) | YALI1F14882g (THS1-like), YALI1A09121g (MST1-like) | each other | 0 |
| R723 Ser-tRNA (cytosol) | YALI1F03908g (SES1-like), YALI1B07957g (DIA4-like) | each other | 0 |
| R673 Pro-tRNA (cytosol) | YALI1D19110g (cytosolic-ProRS-like), YALI1E05929g (YHR020W-like) | each other | 0 |
| R2193 iron import | YALI1D09277g, YALI1D09308g | YALI1D08684g and YALI1D08564g (STATE 2026-10-22 problem 7) | 0 |
| R712 Rib-5-P isomerase | YALI1B09025g (RKI1-like) | YALI1F02518g | 0 |
| R211 chitin synthase | YALI1B21313g, YALI1D04055g, YALI1D34288g | 6 other chitin-synthase genes | 0 |
| R121/R1085 adenylate kinase | YALI1B29895g | YALI1B01136g, YALI1F12969g | 0 |
| R540 GDP-mannose PPase | YALI1C08702g (PSA1-like) | YALI1E18171g (no UniProt in the model) | 0 |
| R4 β-1,3-glucan synthase | YALI1E25018g | YALI1C01944g | 0 |
| R289 protein O-mannosyltransferase | YALI1C32184g | 3 PMT-like genes | 0 |
| R251/R177 sterol P450 reductase | YALI1D06380g (CBR1-like) | YALI1D14199g (OR inside AND rules) | 0 |
| R1843/R1846 LPAAT | YALI1E22736g (SLC1-like) | 3 genes | 0 |
| R60 MTHFR | YALI1B00578g | YALI1B19021g | 0.1156 |
| R815 ADP/ATP carrier | YALI1A10708g | YALI1D10535g, YALI1F26167g | no solution |
| R1556/R216/R238 aconitase | YALI1D11984g (ACO1-like) | YALI1E17935g | 0.0244 |
| R78 + mtFAS | YALI1D32594g (ACP1-like) | 5 genes on R78 | 0 |

Derivation, per gene, is in the CSV columns `closed_by_gene_ko`, `or_partner_genes`, `reaction_ko_ratio` and
`reaction_ko_status`. A reaction-level KO below 15 % proves that the GPR's OR is what keeps the gene
non-essential. It does not prove the OR is wrong: the partner might be a real isozyme.

The aaRS rows also show a **compartment problem**:
- the same rule `YALI1A09121g or YALI1F14882g` sits on both the cytosolic R759 and the mitochondrial R760;
- R760 is blocked in WT.

DIA4 and MST1 are the mitochondrial Ser- and Thr-tRNA synthetases in *S. cerevisiae*. That is recalled and
uncited; the *Y. lipolytica* localisation is not verified, and the literature check found no hits.

### 4.2 Category b: alternative routes

**Verified bypass-closed controls** (`trace_extra.json`, section 8). Each control fixes the model, medium,
objective and threshold:

| FN gene (homolog) | Bypass | WT with bypass closed | KO with bypass closed |
| --- | --- | ---: | ---: |
| YALI1B01136g (ADK1-like), KO/WT 1.0, reaction KO 0.691 | R2024 propanoate:CoA ligase run backwards (AMP + PPi + propanoyl-CoA → ATP + propanoate + CoA); propanoate is exported (R2022 1.18 in the IPP1 trace) | 1.0 (R2024 forward-only) | 0 |
| YALI1F18105g (IPP1-like), reaction KO 0.593 | same R2024 reversal consumes the PPi | 1.0 | 0 |
| YALI1F38116g (ERG13-like), 0.998 | leucine degradation R1968–R1973 makes mitochondrial HMG-CoA; gene-less R1172 exports it | 0.99695 (R1172 closed) | 0 |
| YALI1B11532g (ERG10-like), 1.0 | same R1172 route | 0.99695 | 0 |
| YALI1F00821g (GLN1-like), 1.0 | R2081 has R334's equation, GPR `YALI1D16151g or YALI1F00821g` | — | reaction-level KO of R334+R2081 = 0 (also reports/r2081_closed_f00821_ko_20260911, recalled) |
| YALI1C08262g (GAPDH-like), 0.867 | oxidative PPP + gene-less methylglyoxal synthase R1252 | 1.0 (R1252 closed) | 0.693; R1252 is only part of the bypass |

The ADK1, IPP1, ERG13 and ERG10 rescues depend on two reactions whose direction or existence lacks a stated
mechanism:
- an AMP-forming acyl-CoA synthetase running in the ATP-forming direction;
- an HMG-CoA carrier across the inner membrane with no gene.

Both are **candidates for review, not verdicts**. Thermodynamics and transport evidence have not been checked.

**Complex I** (13 genes, KO/WT 0.7579, verified). Without R2062, cytosolic NADH goes through R570
(`NADH[c] + UQ9[m] → NAD[c] + UQH2[m]`, NDH2-like YALI1F32476g, UniProt F2Z699). The 2026-10-14 report showed
that closing the R533 or R831 outlet drops complex I KOs to 0.069 and that the call then depends on the lysine
cap (recalled). This was not re-run here.

**Central carbon** (26 genes with a growth cost of 0.93–0.999, verified FVA traces where named):

| Knocked step | Bypass | Trace |
| --- | --- | --- |
| PGI | oxidative PPP becomes required (R325, R639, R71) | FVA |
| FBA, GAPDH | PPP, transaldolase/transketolase | FVA |
| ENO, PYK | PEP carboxykinase R634 becomes required (0.05–0.57) | FVA |
| PDH | cytosolic acetyl-CoA routes | FVA |
| SDH, fumarase, KGD | other TCA bypasses | not traced |

These are model-level alternatives. In yeasts several of these steps have no bypass in vivo (recalled,
uncited), so each needs its own biology check.

**Loops.** Four FN genes sit on reactions with |FVA| > 400 at the WT optimum, which means an internal cycle.
The CSV `flags` column marks them:
- YALI1B26749g on R659 (the R659 / r0210 / R_CYST_II cycle reaches 605 unconstrained);
- YALI1A02775g on R793;
- YALI1E35652g on R856;
- YALI1D11984g on R239.

### 4.3 Category c: medium-supplied amino acids (verified, flag f)

| Pathway | FN genes | KO/WT |
| --- | --- | ---: |
| Arginine, plus NCE103-like YALI1F28274g | YALI1B22853g, YALI1C33005g, YALI1D18130g, YALI1D34865g, YALI1F28274g | 0.159 |
| Homoserine | YALI1D14532g, YALI1D01188g, YALI1D16809g | 0.211 |
| Isoleucine/valine | YALI1C32276g, YALI1D03952g | 0.227 |
| Chorismate mutase | YALI1E20751g | 0.255 |
| Threonine | YALI1F17962g, YALI1F30459g | 0.278 |

All 13 KOs grow 0 with the amino-acid and adenine uptakes closed, against a WT of 0.29797.

In the ARG3-like trace, arginine uptake R1202 sits at its cap (−0.00624). Growth is therefore set by the cap,
which `benchmark.md` calls a concentration-scaled assumption, not an uptake rate. Changing caps to move these
calls would be recall tuning and is excluded.

The three homoserine-pathway KOs also become lethal with the vitamins closed. The cause was not traced.

### 4.4 Category d: product not demanded (verified)

`biomass_C` drains only:
- protein, nucleotides and dNTPs;
- chitin, glucan, mannan and trehalose;
- zymosterol and ergosterol;
- 6 lipid classes;
- sulphate;
- UQ9 (1e-4), TPP (1.2e-6), heme a (1e-6), iron (3.04e-5) and two lipoyl terms (1e-5).

It has no NAD(P), FAD, FMN, CoA, folate, pyridoxal phosphate, biotin, sphingolipid, GPI, glycoprotein or
mitochondrial-protein term. Cofactors are recycled, never drained, so their synthesis carries no required flux.

Controls with the 10 vitamin uptakes closed confirm that FOL1-, FMN1-, QNS1- and INO1-like KOs still grow
1.000. The supply is not the medium: there is no demand (verified).

Many of these pathways are also **structurally disconnected**. 51 genes close only reactions that cannot carry
flux at all. The CSV column `dead_end_metabolites` lists, per blocked reaction, the participants with no other
unblocked reaction. Examples:
- mitochondrial aaRS reactions R532, R548 and R166 have dead-end tRNAs (`m797`, `m809`, `m337` have no other
  reaction);
- FOL1-like R14/R274 sit in a closed mitochondrial pterin set;
- LCB1/LCB2 R722 has dead-end ER intermediates.

So adding a demand alone would make these pathways infeasible, not essential, until they are connected. That
consequence is inferred.

V-ATPase (11 genes) and GPI (7 genes) belong here too:
- R794 is unused at the optimum, and its rule has a second branch `(YALI1F38820g and YALI1F13017g)`.
- The 2026-09-11 V-ATPase diagnosis already showed that nothing in the objective requires organelle
  acidification (recalled).
- The GPI demand R_DM_GPI is open and not forced (STATE 2026-10-21).

### 4.5 Category e: genes without reactions (verified)

YALI1B21088g (RAM2-like), YALI1D17983g (RAM1-like), YALI1E16694g (BET4-like), YALI1E01159g (DIM1-like),
YALI1E11370g (PET112-like) and YALI1E11415g (CDC43-like) are in the model with no reaction. Their KO cannot change
growth.

YALI1E11370g's removal from R671 was a reviewed curation (memory note on E11370g/GatB; STATE 2026-10-19). The
other five are protein-prenyltransferase and rRNA-methyltransferase homologs with no small-molecule
reaction in this model.

## 5. Literature checks

Europe PMC retrieval, queries restricted to *Y. lipolytica*. 19 reactions of category a gave 56 gene–reaction
pairs:

| Pairs | Count |
| --- | ---: |
| No hit | 48 |
| 1–4 hits | 7 |
| 5–19 hits | 1 |

Hit table: [hit_counts.csv](hit_counts.csv). Spreadsheet: [review.xlsx](review.xlsx). R171 was left out because
its YALI1M IDs have an ambiguous crosswalk, and the reviewer stops on them.

LLM: `anthropic:claude-opus-5-5`, 8 extractions and 5 chains.

| Step | Prompt tokens | Output tokens |
| --- | ---: | ---: |
| Extraction | 46,862 | 3,215 |
| Chains | 84,693 | 3,608 |

Quotes and PMIDs are code-verified against the retrieved text. The support and evidence labels are model
judgments, and I read every quote before citing it below.

| Question | Result | PMIDs | Bearing on the FN cause |
| --- | --- | --- | --- |
| Is YALI1F02518g (R712 OR partner) a ribose-5-phosphate isomerase? (chain) | unclear, level 1: annotated "putative" RPI by homology to an RpiB-family erythritol isomerase; its deletion delays growth on erythritol. No RPI assay. | [29910781](https://europepmc.org/article/MED/29910781), [34342539](https://europepmc.org/article/MED/34342539) | Weakens the stand-alone-catalyst support for the OR partner that hides RKI1-like YALI1B09025g |
| YALI1B09025g:R712 (extract) | level 2: overexpression of RKI1 raised erythritol; listed as an ecGEM target | [33797654](https://europepmc.org/article/MED/33797654), [41876849](https://europepmc.org/article/MED/41876849) | Supports YALI1B09025g as the main RPI; says nothing on essentiality |
| YALI1B01136g:R121 (extract) | level 1: ADK1 and ADK2 named as two adenylate kinases | [22132183](https://europepmc.org/article/MED/22132183) | Homology only |
| YALI1C08702g:R540 (extract) | level 1: PSA1 (YALI0C06490g) cloned from PO1f | [38053583](https://europepmc.org/article/MED/38053583) | Identity of the FN gene; nothing on partner YALI1E18171g |
| YALI1C33406g:R211 (extract) | level 1: listed as CHS3 | [25816199](https://europepmc.org/article/MED/25816199) | Homology only |
| YALI1D06380g:R737 (extract) | level 1: *S. cerevisiae* CBR1 used heterologously | [42192394](https://europepmc.org/article/MED/42192394) | Not *Y. lipolytica* evidence |
| YALI1D06380g:R251, YALI1A10708g:R815 (extract) | 0 claims | — | — |
| Is Cbr1 required, or replaceable by Ncp1? (chain) | **CHAIN_FAILED**: the provider returned a content-filter refusal; not retried | — | Open |
| Is MST1-like YALI1A09121g mitochondrial, apart from THS1? (chain) | no hits | — | Open |
| Can yeasts take up folate? (chain) | level 0; retrieved text off-topic | — | Open (the d-category result for FOL1 does not depend on it) |
| Is the V-ATPase required in *Y. lipolytica*? (chain) | unclear, level 1: only *M. oryzae* mutant data; *Y. lipolytica* appears only in an alignment | [23826342](https://europepmc.org/article/MED/23826342) | Open |

PMID links are constructed from the PMIDs and were not opened as web pages. The quotes were read from the
reviewer's retrieved full texts in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/essentiality_fn_causes_20261027/out/lit/`.

## 6. What would most likely raise recall: candidates, not changes

Ranked by FN genes affected and by how directly the model evidence points at the cause. The "FN genes at
stake" count is verified for the reaction-level KO. That each gene would flip after an edit is **inferred**: a
GPR edit can also change other genes' calls, including adding FPs. Each candidate must go through
`curation-change` with its own evidence, never because it raises recall.

| # | Candidate | FN genes at stake | Evidence it would need |
| ---: | --- | ---: | --- |
| 1 | R171: the single-gene branches `or YALI1B28025g or YALI1B15733g` | 9 (reaction-level 0.1097) | That ATP19/ATP20-like subunits cannot form an ATP synthase alone (accessory/dimerisation subunits in yeast, recalled). If they become AND members they become essential in the model; they are not on the list, so +2 FP is possible. Dropping them as stand-alone branches is the other option. Also check the YALI0-tagged AND member YALI0D17490g in a YALI1 rule. |
| 2 | Cytosolic aaRS rules R759, R723, R673 (and R760 sharing R759's rule) | 6 (reaction-level 0); about 3 cytosolic genes would flip (inferred) | Localisation of each *Y. lipolytica* paralog: presequence prediction plus experimental data where they exist. The mitochondrial paralogs stay FN unless candidate 5 is done. |
| 3 | R2024 reversibility and the gene-less HMG-CoA carrier R1172 | 4 (verified bypass) | ΔG and direction of AMP-forming acyl-CoA synthetases (Rhea/eQuilibrator). Any evidence for, or against, mitochondrial HMG-CoA export in yeast. |
| 4 | Other OR partners (R712, R60, R2193, R1556 group, R540, R4, R289, R211, R121, R815, R251/R177, R1843/R1846, R78) | 16 | Per partner: evidence that it catalyses the step alone in *Y. lipolytica* (AGENTS §5). R712's partner already has only homology support (PMIDs 29910781, 34342539). |
| 5 | Biomass demand for cofactors, mitochondrial translation, sphingolipids, GPI, glycoproteins | up to 104 (category d) | Measured composition (e.g. Yeast8-style cofactor term) and connected pathways first. 51 genes sit on blocked reactions. Without connection the demand makes WT infeasible. Expect new FPs; needs a separate impact study. |
| 6 | Complex I bypass (NDH2-like R570 and shuttles) | 13 | Already analysed 2026-10-14. The call depends on shuttle direction and the lysine cap. Needs NDH2 orientation and shuttle evidence in *Y. lipolytica*, not a cap change. |
| 7 | GLN1 duplicate R2081 / R334 | 1 | Known since 2026-09-11; needs the duplicate resolved by evidence for YALI1D16151g. |
| — | Amino-acid caps (category c) | 13 | Not a lever: changing caps is tuning. Measured uptake rates on the screen medium would decide whether these are real model errors. |

## 7. Verification levels and limits

**Verified this session:**
- the screen and its counts at the four cutoffs;
- every category assignment (mechanical rule on saved LP outputs);
- the reaction-level KOs, medium controls and bypass-closed controls;
- the biomass composition, the GPR texts and the blocked set;
- the literature quotes and PMIDs, read from the retrieved texts.

**Inferred:**
- that a given edit would flip a given gene;
- the "most likely" ranking;
- which central-carbon bypasses exist in vivo;
- the dead-end consequence of adding demands.

**Recalled, not re-run:**
- the complex I outlet result (2026-10-14);
- the V-ATPase diagnosis (2026-09-11);
- the R2081 result (2026-09-11);
- *S. cerevisiae* localisation of DIA4/MST1 and the role of ATP19/ATP20.

**Limits:**
- **One optimum.** pFBA route lists show one optimum. Many unrelated KOs share the same swap (R_MDHy → R1377 /
  R_ALCDH_nadp_hi), which is degenerate, not mechanism. FVA at 99 % was used where a route is named. FVA is
  unreliable near 1e-7 (memory note); the tolerance here is 1e-6.
- **Thresholds.** The "unused at optimum" rule uses FVA at 100 %, so a reaction that is only needed below the
  optimum would be called d. Reaction-level KO guards against that only for category a.
- **Keyword sub-labels.** The d sub-labels and the out-of-scope breakdown are keyword matches on the list's
  homology text, not verified identities. Gene names in this report are homolog names from that text. No
  *Y. lipolytica* gene name is verified ("none verified" in the CSV).
- **Not tested:**
  - double knockouts within the large paralog families (chitin synthase, PMT);
  - other media;
  - the 11 unresolved positives (complex III/IV, STATE 2026-10-21).
- **Literature.** Coverage is thin: 48 of 56 pairs had no *Y. lipolytica* hit. One chain was refused by the
  provider's content filter and was not retried.
- **Out of scope.** 1,286 listed genes have no model counterpart and are outside this analysis by
  `benchmark.md`. A keyword match on their list text gives:

  | Function (keyword match) | Genes |
  | --- | ---: |
  | Ribosome / rRNA | 135 |
  | Transcription | 85 |
  | Secretion | 55 |
  | DNA replication | 48 |
  | Chromosome segregation | 40 |
  | Proteasome | 37 |
  | Cytoskeleton | 36 |
  | Translation factors | 34 |
  | Nuclear transport | 33 |
  | Splicing | 31 |
  | No keyword match | 704 |

  These are mostly non-metabolic and cannot be recovered by a metabolic model.

**Could the method be improved?**
- A minimal-cut-set or synthetic-lethal-pair search per FN gene would replace the reaction-level KO heuristic.
  It gives the exact set of reactions whose joint removal is lethal; test it by reproducing the five verified
  bypasses.
- Loopless FVA would remove the four loop flags' ambiguity.

## 8. Files

In this folder:
- TASK.md, this REPORT.md;
- `fn_causes_a1237821.csv`;
- the five scripts;
- `screen_a1237821.manifest.json`;
- the reviewer's `hit_counts.csv`, `review.xlsx` and `manifest.json`.

In `$IYALI26_RESEARCH_ROOT/artifacts/tasks/essentiality_fn_causes_20261027/` (23 MB):
- `input/` holds the model export and its build provenance;
- `screen_a1237821/` holds `gene_screen.tsv` (`150d7cb5…`), `WT_fluxes.tsv` and `solver_calls.jsonl`;
- `fn_classification/` holds `fn_classification.json`, `trace_controls.json` and `trace_extra.json`;
- `out/lit/` holds the retrieval files, `*.review.json`, `*.chain.json` and `usage.jsonl`.

## 9. Commands

From the repository root, with `IYALI26_RESEARCH_ROOT` set and
`PY=.claude/worktrees/workflow-automation-review-6765da/platform/.venv/bin/python`:

```
git show origin/main:model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_20261023.xml > "$IYALI26_RESEARCH_ROOT/artifacts/tasks/essentiality_fn_causes_20261027/input/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_20261023.xml"
$PY -B model/reports/essentiality_fn_causes_20261027/run_screen.py a1237821
$PY -B model/reports/essentiality_fn_causes_20261027/classify_fn.py
$PY -B model/reports/essentiality_fn_causes_20261027/trace_controls.py
$PY -B model/reports/essentiality_fn_causes_20261027/trace_extra.py
python3 -B model/reports/essentiality_fn_causes_20261027/build_table.py
cd platform
env -u ANTHROPIC_API_KEY $PY -m gem_review.lit retrieve --model <exported xml> --reactions R1085,R121,R177,R1843,R1846,R211,R2193,R251,R289,R4,R540,R60,R673,R712,R723,R737,R759,R760,R815 --task essentiality_fn_causes_20261027 --research-root "$IYALI26_RESEARCH_ROOT"
env -u ANTHROPIC_API_KEY $PY -m gem_review.lit extract --task essentiality_fn_causes_20261027 --research-root "$IYALI26_RESEARCH_ROOT" --pairs YALI1B01136g:R121,YALI1C33406g:R211,YALI1D06380g:R251,YALI1C08702g:R540,YALI1B09025g:R712,YALI1F02518g:R712,YALI1D06380g:R737,YALI1A10708g:R815
env -u ANTHROPIC_API_KEY $PY -m gem_review.lit chain --task essentiality_fn_causes_20261027 --research-root "$IYALI26_RESEARCH_ROOT" --pair YALI1F02518g:R712
env -u ANTHROPIC_API_KEY $PY -m gem_review.lit report --task essentiality_fn_causes_20261027 --research-root "$IYALI26_RESEARCH_ROOT"
```

The four question chains used `chain --question … --gene … --aliases … --terms …`. Their exact arguments are in
the `question`, `gene` and alias fields of each `out/lit/questions/q_*.chain.json`.

The `trace_extra.py` run that produced the saved `trace_extra.json` was its second run, after the ERG10 control
was added; the first run's output was overwritten.
