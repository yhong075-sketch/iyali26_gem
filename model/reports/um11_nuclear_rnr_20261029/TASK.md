# TASK: UM11, the nuclear ribonucleotide-reductase block (2026-10-29)

## 0. Authorization and scope

- **Authorization:** the user, in this session (2026-10-29), after reading the Stage 1 draft (sections 1–9
  below): **"Option B"**.
  - Option B means keeping the nuclear RNR copies (R707, R709, R711, R1914), R782 and the six nuclear nucleotide
    transports as they are, blocked, and recording why. There is no model change.
  - Stage 2 for B is documentation only:
    - this folder;
    - the UM11 entry of the gap-census report;
    - an Escher figure of the as-built module next to the Option A what-if (not implemented);
    - an independent read-only audit;
    - a `model/STATE.md` entry.
- **Not in scope:**
  - any model, curated-data or code change;
  - Options A, A+, A′, C1 and C2;
  - the RNR GPRs, R702–R705, R1914's water, and the thioredoxin formula (follow-ups in section 9).
- **No commit or push is authorized.**
- **Baseline:** `14f73bda…` (section 1), built earlier in this session on the user's instruction. It is
  documented in `model/reports/combined_um7_r2171_20261029/`.
- **Budget for Stage 2:** 0 builds; about 2,000 LPs for the figure. The figure also runs whole-model FVA for the
  node colours, as the UM2 figure did; that adds about 9,300 LPs.
- **Acceptance criteria:**
  - no model file written;
  - the decision recorded in the census report and in `model/STATE.md`;
  - the figure's WT reproduces 0.3116608216 within 1e-9;
  - the viewer draws in jsdom and in headless Chrome.
- Sections 1–9 are the Stage 1 draft as the user read it. Paths now point to this folder:
  - Stage 1 outputs are in `stage1/` (logs, `evidence.md`, BLAST tables, knockout summaries) and `scripts/`.
  - The large files are in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/um11_nuclear_rnr_20261029/stage1/`: the
    PubMed / UniProt / SGD fetches (`sources/`), `stage1_um11.json`, `blocked_sets.json`, and the full knockout
    tables and stderr logs.

Evidence tags:
- **verified** = checked this session, with the output kept;
- **inferred** = reasoned;
- **recalled** = not opened this session.

**Recommendation: Option B.** Keep the nuclear RNR copies blocked and record why, with no model change, as
UM2 did. In memory, every option that frees the module (A, A+, A′) is neutral: no WT change, no knockout
change, no energy cycle. The choice therefore rests on biology. In yeast, RNR catalysis is cytoplasmic:
- Rnr1 is cytoplasmic.
- The nuclear Rnr2/Rnr4 pool is held there inactive and moves to the cytoplasm to form the holoenzyme.
- No source places RNR catalysis in the nucleus.

Opening the module would make the model claim a nuclear RNR activity that no evidence supports. Every
ancestor model (iYali4, iYli21, iYali5, Yeast-GEM 9.0.2, iMM904) has the same nuclear copies, and they are
blocked in all of them.

## 1. Starting state (verified, 2026-10-29)

- **Worktree:** `claude/wizardly-wing-a0d6c0` at `ff5192d`, which equals `origin/main` after `git fetch`.
  It was clean at the start.
- **Baseline built in this session.** The user agreed: "Yes, build combined". The baseline is
  `model/candidates/E5_…_dolichol_route_v2_um7_nuclear_diffusion_r2171_1to1_20261029.xml`, SHA-256
  `14f73bda4b0772019e82ea18a27e198c5ba175e67fef205a51dbfb090809e757`:
  - It equals `9d655c2b` plus `--um7-nuclear-diffusion` plus `--r2171-pi-plc-1to1`.
  - Seeds 0 and 1 give the same bytes.
  - The diff against `5f54eb22` is exactly the R2171 change (2 entries). The diff against `965ab947` is
    exactly the UM7 change (3 entries). Both are field-identical to the parents' own diffs.
  - WT is 0.3116608216, the same as both parents.
  - It is not a designated reference.
  - See `model/reports/combined_um7_r2171_20261029/REPORT.md` in the worktree.
- **Before the user's answer** only the census numbers of `5f54eb22` were used. No LP was run on it.
- **Instruction files read:**
  - `~/.claude/CLAUDE.md` and the worktree's `AGENTS.md`;
  - `model/expected/benchmark.md`;
  - the `model/STATE.md` entries of 2026-10-26 to 2026-10-28, including UM7, UM2 and R2171.
- **Condition:**
  - Task A medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` with the PO1f profile
    `state/strain_profiles/po1f_sd_leu.json` (uracil uptake 1000, R612 closed).
  - Experimental list: `consensus_essential_genes.csv`.
  - "Open" means the raw model with every exchange at (−1000, 1000).
  - Gurobi 13.0.2 at 1e-9 with the real HOME, through `checked_solve.py` (reset before each solve, stale-point
    check).
  - LPs read short-path copies of the candidates, with SHAs re-checked. cobra rejects paths of 260
    characters or more on macOS.
- **Essentiality rules:**
  - A gene is essential when KO/WT is strictly below 15 %; 1, 5 and 10 % are also reported.
  - Model genes not on the list count as non-essential; this is the user-defined negative class.
  - The list was used in development, so scores are a development reference, not independent validation.

## 2. The 11 reactions (verified, `model_inventory`, `stage1_um11.log`)

Bounds are as written. **Blocked** means |v| ≤ 1e-9 in reaction-scoped FVA at 0 % of the optimum, under both
open exchanges and Task A + PO1f, on `14f73bda`.

| Reaction | Equation (all nucleus unless marked) | Bounds | GPR | Blocked open / Task A |
| --- | --- | --- | --- | --- |
| R707 RNR (ADP) | ADP + TRX1 → dADP + TRX1-disulphide + H2O | 0, 1000 | D04261g or F39330g or B05072g | yes / yes |
| R709 RNR (CDP) | CDP + TRX1 → dCDP + TRX1-disulphide + H2O | 0, 1000 | same | yes / yes |
| R711 RNR (GDP) | GDP + TRX1 → dGDP + TRX1-disulphide + H2O | 0, 1000 | same | yes / yes |
| R1914 RNR (UDP) | UDP + TRX1 → dUDP + TRX1-disulphide + **H2O [cytosol]** | 0, 1000 | F39330g or (B05072g and F02309g) | yes / yes |
| R782 UMP kinase | UMP + ATP ⇌ UDP + ADP | −1000, 1000 | F12969g | yes / yes |
| R1033 / R1136 | CDP / GDP cytosol ⇌ nucleus | ±1000 | none | yes / yes |
| R1076 / R1077 / R1093 | dADP / dCDP / dGDP cytosol ⇌ nucleus | ±1000 | none | yes / yes |
| R1353 | UMP cytosol ⇌ nucleus | ±1000 | none | yes / yes |

GPR prefix `YALI1`, suffix `g`. The four RNR reactions have EC 1.17.4.1. R1914 is in "Pyrimidine
metabolism"; R707/R709/R711 have no subsystem. The six transports are SBO:0000185 in "Transport, nucleus".

Each nuclear copy has **the same GPR as its cytosolic counterpart**:
- R707 / R709 / R711 match R706 / R708 / R710;
- R1914 matches R1913;
- R782 matches R781.

So a nuclear copy cannot become a gene-level bypass of a cytosolic one. **verified**: the GPR strings are
identical, and the sweeps in section 5 agree.

Genes (identity by reciprocal best BLAST hit, W29 proteome UP000182444 `a0bacab5…` vs S288C `1aaa8186…`;
`evidence.md` Q1):

| Gene | Verified name | Function | Evidence level |
| --- | --- | --- | --- |
| YALI1F39330g | none verified (Q6BZQ3, unreviewed) | RNR large subunit α | ortholog of Rnr1 by RBH, 67.5 % identity; on the consensus essential list |
| YALI1D04261g | none verified (Q6CAF2, unreviewed) | second RNR large-subunit paralog, 57 % identical to F39330g, catalytic residues kept | homology; RBH of no S. cerevisiae protein |
| YALI1B05072g | none verified (Q6CFU6, unreviewed) | RNR small subunit β (iron ligands and radical Tyr kept) | ortholog of Rnr2 by RBH, 66.2 %; no Rnr4 ortholog in W29; on the consensus list |
| YALI1F02309g | none verified (Q6C399 "Thioredoxin", unreviewed) | thioredoxin (WCGPC) | ortholog of Trx2 by RBH |
| YALI1F12969g | URA6 per unreviewed Q6C2B0 | UMP/CMP kinase | ortholog of Ura6 by RBH |
| YALI1D35735g | TRR1 (Q6C7L4, Swiss-Prot reviewed; function by similarity) | thioredoxin reductase, the only TRR in W29 | ortholog of Trr1 by RBH |

## 3. Supply per compartment and the roots (verified, `model_inventory`, `inventory_5f54eb22.txt`)

| Species | Nucleus | Cytosol | Nucleus ↔ cytosol |
| --- | --- | --- | --- |
| TRX1 reduced (`m931` / `m632`) | consumed by R707, R709, R711, R1914; **no producer** (RNP root) | 12 reactions; regenerated by TRR R755 (`D35735g and F02309g`) | **none** (the peroxisome has R1313) |
| TRX1-disulphide (`m933` / `m633`, no formula) | made by the 4 RNR copies; **no consumer** (RNC root) | consumed by R755 | **none** (the peroxisome has R1286) |
| CDP, GDP | consumed only by R709 / R711 | active | R1033, R1136 |
| UMP → UDP | UMP: R1353, R782; UDP: made by R782, used only by R1914 | active | R1353 (UMP); none for UDP |
| dADP, dCDP, dGDP | made only by R707 / R709 / R711 | active | R1076, R1077, R1093 |
| dUDP | R1914, R789 (nuclear dUMP kinase) | active | R1104 (active) |
| NADPH / NADP | **no nuclear species** | active | — |
| glutaredoxin | **no species in the model** | — | — |

- Both roots are the nuclear thioredoxin pair; nothing in the nucleus makes reduced thioredoxin or uses the
  disulphide.
- **Elementary fixes, tested in memory (section 5):**
  - The reduced link alone frees 0/11, and the disulphide link alone frees 0/11.
  - Both links together free 11/11.
  - A nuclear TRR copy, which needs NADPH/NADP nuclear species and links, also frees 11/11.
- **R1914 compartment slip.** R1914 writes its water to the cytosol, while every other species is nuclear. The
  slip dates back to iYali4 (`R_y200021`, verified in `model_scan.txt`). It has no LP effect, because nuclear
  water exchanges freely through R1361.
- **Unbalanced-looking chemistry.** TRX1 has the formula C6H9NO2S2R2 and the disulphide has none, so no RNR
  copy's mass balance can be checked. This is unchanged by every option.

## 4. Biology (evidence subagent this session; links in `evidence.md`; spot-checks marked)

**Y. lipolytica:** there are **no** RNR localization or biochemistry data. PubMed returned 0 hits for
RNR or dNTP; **verified**. The identities above come from homology only.

**S. cerevisiae, where RNR acts.** These abstracts were fetched this session; the main session itself read
Yao 2003, Lee 2006 and Camier 2007 from the saved fetch.
- In a normal cycle, Rnr1/Rnr3 are mainly cytoplasmic and Rnr2/Rnr4 mainly nuclear. Under stress,
  Rnr2/Rnr4 move to the cytoplasm: Yao 2003, https://pubmed.ncbi.nlm.nih.gov/12732713/.
- Wtm1 anchors Rnr2/Rnr4 in the nucleus; on damage they join cytoplasmic Rnr1: Lee & Elledge 2006,
  https://pubmed.ncbi.nlm.nih.gov/16452505/.
- Dif1 imports R2 into the nucleus: Lee 2008, https://pubmed.ncbi.nlm.nih.gov/18851834/; Wu & Huang 2008,
  https://pubmed.ncbi.nlm.nih.gov/18838542/.
- In S phase, Clb6–Cdc28 releases R2, which forms the active holoenzyme with Rnr1 in the cytoplasm: Wu 2018,
  https://pubmed.ncbi.nlm.nih.gov/29263158/.
- SGD gives RNR1–4 "is active in cytoplasm" and none "is active in nucleus".
- **Mammals:** RNR is cytosolic and dNTPs diffuse into the nucleus (Pontarin 2008,
  https://pubmed.ncbi.nlm.nih.gov/18997010/, read by the main session). No nuclear reductase activity was
  found (Leeds 1985, https://pubmed.ncbi.nlm.nih.gov/3915777/).
- **Counter-report:** RRM1/RRM2 are recruited to damage sites in human cells (Niida 2010,
  https://pubmed.ncbi.nlm.nih.gov/20159953/). This is human, damage-specific, and does not measure catalysis.
- **Conclusion:** no source places RNR catalysis in the yeast nucleus. The nuclear R2 pool is sequestered and
  inactive. **verified** (abstracts).

**Thioredoxin:**
- Thioredoxins are the major in vivo RNR reductants, with glutaredoxin as backup (Camier 2007,
  https://pubmed.ncbi.nlm.nih.gov/17349928/, read by the main session; Koc 2006,
  https://pubmed.ncbi.nlm.nih.gov/16574642/).
- Trx1/Trx2 are called nuclear only by high-throughput GFP (SGD HDA; UniProt citing Huh 2003, which was **not
  opened**). They are otherwise cytosolic (Pedrajas 1999, https://pubmed.ncbi.nlm.nih.gov/10037727/).
- At 11–14 kDa they lie below the 27 kDa GFP that equilibrates through the yeast NPC in about 15 s (Timney 2016
  full text, https://pubmed.ncbi.nlm.nih.gov/27697925/). That thioredoxin diffuses is **inferred** by
  extrapolation. Most Trx1 sits in the LMA1 complex, whose size was not checked.
- **There is no nuclear thioredoxin reductase:** SGD TRR1 has no nucleus annotation of any kind
  (**verified** by the subagent). This argues against A′.
- **Nucleotides:** small molecules cross the NPC passively (Timney 2016; Wente & Rout 2010,
  https://pubmed.ncbi.nlm.nih.gov/20630994/). No source measures NDP, dNDP or UMP, so applying the rule to
  them is **inferred**.

**Model history** (local SBML files, `model_scan.txt`; **verified**, static):
- iYali4 v4.1.2, iYli21, iYali5, Yeast-GEM 9.0.2 and iMM904 all carry nuclear RNR copies.
- None of them makes reduced nuclear thioredoxin, so the copies are blocked in all of them (stoichiometry; no
  LP run on those models).
- Yeast-GEM's only link is the gene-less disulphide transport `r_4612`, noted "Confidence Level: 0 … added
  after Gapfilling … unknown transporter".
- The iYali line inherited the copies and the OR-of-subunits GPR from the Yeast-GEM family (**inferred**). The
  nuclear UDP copy (R1914) is a Yarrowia-model addition.

## 5. Options tested in memory (verified; `stage1_um11.log`, `blocked_sets.json`, `ko/`)

- **Baseline `14f73bda`:**
  - WT 0.3116608216 (optimal, check passed).
  - Whole-model FVA: 695 blocked under open exchanges, 828 under Task A; all 11 UM11 reactions blocked.
  - No-uptake ATP (R_NTP1) and biomass are 0: status optimal, with `xMAINTENANCE` 7.8625 and all positive
    lower bounds relaxed.
  - At 15 %: TP / FN / FP / TN = 112 / 203 / 42 / 692. Also 16 unresolved (infeasible knockouts, the same set
    in every run), 4 mtDNA genes out of scope, 326 in-model positives.

| Option | Model change | UM11 freed (open / Task A) | Whole-model blocked | ΔWT | Knockouts (2 × 1,081, all 4 cutoffs) | EGC |
| --- | --- | --- | --- | --- | --- | --- |
| **A**: TRX1 + disulphide diffusion | 2 gene-less links | 11 / 11 | 695 → 684, 828 → 817; exactly the 11 freed, none newly blocked, both links flux-capable | −1e-15 | no class change; max \|ΔKO\| 1.6e-13; no gain > 1e-6 | 0 / 0, optimal |
| **A+**: A + R1914 water to the nucleus | A + 1 stoichiometry fix | 11 / 11 | same as A | −1.7e-15 | same as A | 0 / 0, optimal |
| **A′**: nuclear TRR copy + NADPH/NADP links | 2 new species, 3 reactions | 11 / 11 | (not scanned) | −1e-16 | no class change; max \|ΔKO\| 1.6e-13 | 0 / 0, optimal |
| **B**: no change, documented | none | 0 / 0 | unchanged | 0 | unchanged | 0 / 0 |
| **C1**: close the 4 RNR copies + R782 with a reason note | bounds + notes on 5 reactions | 0 / 0 | unchanged | −2e-15 | not swept: provably none, because the reactions are blocked in WT and any knockout only shrinks the feasible set | 0 / 0, optimal |
| single links (reduced only / disulphide only) | 1 link | 0 / 0 | — | 0 | — | 0 / 0 |

In pFBA the nuclear copies carry 0 at the optimum under A, so they never help growth. Only the split
between R706 and R702 for dADP / dATP shifts; that is an alternate optimum of equal cost. Single-link and A′
results come from the same scripts.

**Why nothing moves (mechanism, verified):**
- In the nucleus the RNR products have no consumer except the transports back to the cytosol (and R789 for
  dUDP). The biomass draws its deoxynucleotides from the cytosol.
- So an open nuclear RNR can only run in parallel with the cytosolic copy, with identical genes, the same
  thioredoxin and the same TRR.
- It adds no new capability; each option only adds or removes flux capacity.

**Targeted controls on the RNR genes** (`rnr_controls.log`, Task A + PO1f, baseline and A):

| Case | Growth ratio |
| --- | --- |
| each single RNR knockout | 1.0 |
| all three RNR genes knocked out | 1.0 |
| R702–R705 closed | 1.0 |
| one gene knocked out + R702–R705 closed | 1.0 |
| all three knocked out + R702–R705 closed | **0** |

So the consensus-essential RNR genes YALI1F39330g (Rnr1-like) and YALI1B05072g (Rnr2-like) are false
negatives for two independent reasons:
- the GPR on R706/R708/R710 is an OR of α and β subunits;
- R702–R705, a "ribonucleoside-triphosphate reductase" with thioredoxin as its only gene, copied from
  Yeast-GEM r_0970–0973, reduce the NTPs directly.

No option here changes that. This is a separate curation question and **affects another line**.

## 6. Options: evidence for and against

**A: two gene-less nucleus ⇌ cytosol diffusions of thioredoxin, reduced and disulphide (the UM7 template)**
- For:
  - It frees 11/11 and is neutral in every test.
  - Thioredoxin is small enough to diffuse (inferred).
  - The model already links the peroxisome this way (R1313 / R1286).
- Against:
  - It turns on nuclear RNR catalysis, which no yeast source supports.
  - Nuclear thioredoxin is only a high-throughput call, and the nuclear R2 pool is inactive.
  - The 11 "freed" reactions would carry flux only as parallel copies, never in an optimum.
  - Yeast-GEM's analogous gap-fill is confidence 0.
  - It adds two links for a function the cell is not known to have.

**A+:** A plus the R1914 water correction. Same as A for the LP; it also fixes an inherited compartment slip.

**A′: nuclear TRR copy**
- Against: no nuclear Trr in yeast (SGD), new nuclear NADPH/NADP species, and more entities than A for the
  same result.
- Not recommended.

**B: keep blocked, documentation only (as UM2 A0)**
- For:
  - It matches the yeast localization literature (cytoplasmic holoenzyme; dNTPs reach the nucleus by
    diffusion) and every ancestor model.
  - No model change.
- Against: the census keeps counting 11 blocked reactions, and the R1914 slip stays.

**C1: close the four RNR copies and R782 (bounds 0, 0) with a reason note; no deletion**
- For:
  - It makes the "no nuclear catalysis" decision explicit in the model.
  - It guards against a later gap-fill silently opening them.
  - Provably LP-neutral.
- Against:
  - It is a model change with no functional effect.
  - The six transports stay blocked anyway, because their only partner is the nuclear RNR.
  - Closing is not deletion, which would need evidence of non-existence that the project's rule does not
    have.

**C2: R1914 water compartment fix on its own**
- An inherited chemistry slip, LP-neutral while R1914 is blocked.
- It can be added to B (as a separate small change) or to C1.

Not proposed:
- deleting the nuclear copies or transports (no evidence of non-existence);
- changing any GPR or R702–R705 (outside UM11; see section 5).

## 7. Acceptance criteria for Stage 2 (each a stop condition; on failure, write STOP.md and ask)

All options:
- WT on Task A + PO1f within 1e-9 of 0.3116608216;
- no knockout class change at 1 / 5 / 10 / 15 %;
- no knockout growth change > 1e-6;
- EGC no-uptake ATP and biomass 0 with status optimal (maintenance relaxed);
- pytest passes as before, apart from the new tests;
- `14f73bda` reproduced byte for byte without the new flag;
- seeds 0 and 1 give identical bytes.

- **B:** no model file is written. The census REPORT UM11 entry and `model/STATE.md` record the decision. Escher:
  as built vs the A what-if (not implemented), on Task A + PO1f, showing capacity.
- **A / A+:** the whole-model diff is exactly 2 added reactions, `R_TRX1tn` and `R_TRXSStn`, plus "Transport,
  nucleus" group membership (A+ also the R1914 stoichiometry and notes). Whole-model FVA frees exactly the 11
  UM11 reactions, none newly blocked, under both conditions.
- **C1 / C2:** the diff is exactly the bounds and notes of the named reactions, or R1914's stoichiometry and
  notes; the blocked sets are unchanged.

Budget for Stage 2:
- B: 0 builds and about 2,000 LPs (figure only).
- A, A+ or C: 3 builds (candidate seeds 0/1 and a control) and about 35,000 LPs (sweep, FVA, EGC, figure).

## 8. Stage 1 accounting

- **Builds:** 2, the combined baseline at seeds 0 and 1.
- **LPs:** about 35,600 on Gurobi:
  - stage1 1,307;
  - whole-model FVA 3 variants × 2 conditions ≈ 27,800;
  - knockout sweeps 3 × 2,164;
  - RNR controls 20.
- Every knockout and WT point passed the stale-point check; the 20 infeasible knockouts are the same in every
  run.

## 9. Follow-ups found (not changed; each affects another line)

1. **RNR GPR.**
   - R706/R708/R710 (and the nuclear copies) use an OR of α and β subunits; R1913/R1914 use
     `F39330g or (B05072g and F02309g)`.
   - The candidate form is (large) **and** (small).
   - Whether D04261g can replace F39330g is unknown.
2. **R702–R705.**
   - These are NTP reductases with a thioredoxin-only GPR, from Yeast-GEM r_0970–0973.
   - Together with item 1, they make the consensus-essential RNR genes F39330g / B05072g false negatives
     (controls in section 5).
3. **R1914's water in the cytosol**, inherited from iYali4 (Option C2).
4. **TRX1-disulphide has no formula**, so no thioredoxin reaction's mass balance can be checked.
5. **A second thioredoxin pair, `m2033` / `m2034`**, used by R2198 and R2204, duplicates TRX1 / TRX1-disulphide.

## 10. Corrections after the audit (2026-10-29)

Sections 1–9 above are the draft as the user read it. The audit ([AUDIT.md](AUDIT.md)) found these points; the
corrected statements are in REPORT.md.
- **Section 5:** "4 mtDNA genes out of scope" should read: all 16 mtDNA genes are excluded, and 4 of them are
  among the 20 infeasible knockouts. The counts are unchanged.
- **Section 5, A′ "(not scanned)":** the audit scanned A′ over the whole model. It frees the same 11, with none
  newly blocked (817 Task A, 684 open).
- **Section 4, "with glutaredoxin as backup":** the Camier 2007 abstract says only that either the thioredoxin
  or the glutaredoxin pathway must be present for viability.
- **Section 4, "the nuclear R2 pool is sequestered and inactive", tagged verified:** sequestration is in the
  abstracts. Inactivity is the authors' model (inferred).
- **Section 4, "dNTPs diffuse into the nucleus" (Pontarin 2008):** the authors suggest this; they do not measure
  it.
- **Section 4, thioredoxin "nuclear only by high-throughput GFP":** this holds for Trx1 (SGD HDA). For Trx2, SGD's
  nucleus call is IEA only.

