# Reviewed GPRs overwritten by the 2026-09-10 metadata selection: review and fix plan

Read-only. No model, data or code was changed, nothing was built, no LP was solved. Scope, inputs and
SHAs: [TASK.md](TASK.md). Per-row table: [selection_rows.csv](selection_rows.csv). Scripts:
[scripts/](scripts/).

## 1. Summary

- **Where the selection's GPRs come from.** The selection copies 9 GPRs from a `metadata.xml` built
  on 2026-09-09 from commit `994b09bf` (`metadata.build.json`, `head`). That builder predates two
  reviewed patches that exist only on `main`'s history:
  - `remove_misannotated_gprs` (`7eacb58`, 2026-07-06);
  - the CPA/URA2 partition EG-GPR-001 (`1c9f202`, 2026-07-23).

  `git merge-base --is-ancestor` confirms neither is an ancestor of `994b09bf`, and
  `git show 994b09bf:scripts/gem_annotate/patches.py` contains neither function. The selection's
  "before" is the builder output with the patches; its "after" is the older output without them.
  The selection records itself as "User-selected model version; not new validation of chemistry or
  native gene function". **No GPR row is newer reviewed evidence.**
- **R765/R766: the reviewed patch is wrong, not the selection (verified).** YALI1E07744g is the
  full-length transketolase, not a trehalase. The 2026-07-06 patch matched the YALI1 locus suffix
  to the YALI0 locus of the same number, a different gene (§3.1). The selection's
  `YALI1E07744g or YALI1D02625g` contains the real enzyme. Re-applying the selection's "before"
  would remove it.
- **R671: the reviewed patch holds (verified).** YALI1E11370g is GatB/PET112. Its OR membership on
  prephenate dehydrogenase has no support.
- **R190/R159: the selection undoes accepted legacy patch EG-GPR-001.** Its "after" state is
  inconsistent with gene identity:
  - CPA1 and CPA2 are OR members of the ATCase reaction, which they lack;
  - CPA1 and CPA2 are interchangeable, although they form an obligate heterodimer.

  Whether URA2 also feeds the arginine pool is unresolved for *Y. lipolytica* (§3.3).
- **a56f3649 and c98961e9 are internally inconsistent on R159/R190.** The partition's names, EC
  numbers, subsystems and the m325 note survive. The selection reverted the stoichiometry and GPRs
  underneath them.
- **Plan (§5):**
  - re-apply R671's GPR;
  - correct the builder's wrong R765/R766 removal, with no model change expected;
  - decide R190/R159 among three options.

  Each piece is its own change set in a generalized form of the 2026-10-19 step. The audited
  SELECTION-GPR-RESTORE-20261019 is not touched.

## 2. Every selection row that touches a GPR or a reviewed patch

The selection has 207 reactions and 216 field rows: 201 stoichiometry, 6 bounds, 9 GPR. Each row was
compared against three states (`scripts/classify.py`):
- the source model (S);
- the selection's "before" (B, the builder state it guards);
- its "after" (A).

The table also shows whether a56f3649 still carries A.

### 2.1 The 9 GPR rows

| Reaction | S (source) | B (selection before) | A (selection after = a56f3649) | Reviewed patch behind B | Classification |
| --- | --- | --- | --- | --- | --- |
| R612 | `""` | `YALI1E31685g` | `""` | curated GPR (URA3) | undoes reviewed patch. Restored 2026-10-19 (opt-in, c98961e9) |
| R2274 | `""` | `YALI1A21372g` | `""` (+ drops 3 H₂O) | curated GPR (LIP2) | same |
| R1302 | `""` | `YALI1A08512g` | `""` | curated annotation | same |
| R1303 | `""` | `YALI1A08512g` | `""` (+ adds H⁺) | curated annotation | same |
| R671 | `YALI1F23441g` | `YALI1F23441g` | `YALI1F23441g or YALI1E11370g` | `remove_misannotated_gprs` | **undoes a reviewed patch that holds up** |
| R765 | `YALI1E07744g` | `YALI1D02625g` | `YALI1E07744g or YALI1D02625g` | `remove_misannotated_gprs` | **undoes a reviewed patch that is itself wrong** (§3.1). A is the better-supported state |
| R766 | `YALI1E07744g` | `YALI1D02625g` | `YALI1E07744g or YALI1D02625g` | `remove_misannotated_gprs` | same as R765 |
| R159 | `C33005g or D09420g or E11768g` | `YALI1E11768g` (net CPS+ATCase) | `C33005g or D09420g or E11768g` (ATCase only) | EG-GPR-001 partition | **undoes accepted legacy patch**. A conflicts with gene identity |
| R190 | `C33005g or D09420g or E11768g` | `YALI1C33005g and YALI1D09420g` | `C33005g or D09420g or E11768g` | EG-GPR-001 partition | **undoes accepted legacy patch**. URA2 membership unclear (§3.3) |

How the R671 and R765/R766 "after" values arose:
- `gpr_isozyme_additions.csv` adds YALI1E11370g to R671. Its R671 row carries EC 6.3.5.7, the
  GatCAB EC. The legacy `model.xml` annotates R671 with that EC; the source model does not
  (corrected 2026-10-19 after the fix audit). The table also adds YALI1D02625g to R765/R766.
- `remove_misannotated_gprs` then removes YALI1E11370g and YALI1E07744g.
- The 2026-09-09 builder ran only the first step. That is the selection's "after".

### 2.2 Stoichiometry and bounds rows that interact with reviewed patches

A row interacts with a reviewed patch when it supersedes a note left by a reviewed builder step.
The selection's `superseded_notes` field is non-empty for 17 reactions. The first two rows below
repeat the equation side of GPR rows from §2.1; the rest are new:

| Reaction(s) | Reviewed step | What the selection does | a56f3649 state | Classification |
| --- | --- | --- | --- | --- |
| R1303, R2274 | curated GPR rows above | adds H⁺ / drops 3 H₂O | selection state; fixed in c98961e9 | covered 2026-10-19 |
| R159 | EG-GPR-001 | replaces the net channelled reaction with plain ATCase | selection state | part of the R159/R190 decision |
| **R2041** | `fix_d_arabinokinase_direction_and_proton` (ExPASy EC 2.7.1.54) | bounds [0,1000] → [−1000,1000]; H⁺ 1 → 2 | selection state (reversible) | **undoes a reviewed direction patch**. Not a GPR; see §6 |
| R_NTP1 | gap-fill direction curation | excluded by the 2026-09-14 amendment | reviewed state | already handled |
| R_NDP1, R_NTP3pp, R_NTP7, R_PGAM1_PhosHydro, R_CAT2p, R_OAADCm, r0242 | gap-fill direction curation / energy candidate | reversible bounds and/or H⁺ changes | energy candidate's protected definition (a56f3649 ≠ A) | the reviewed state survives on the E5 chain. Builds without the energy candidate were not checked |
| R_CYST_II, R_NDP7g | none: superseded note is `gap_fill_direction_status = legacy_unreviewed` | H⁺ change | selection state | not a reviewed patch |

Most stoichiometry rows change only H⁺ (sometimes H₂O) bookkeeping. The builder's automatic
`fix_proton_water_balance` sits on one side and the older builder's balance on the other. Across
all 201 stoichiometry rows, mass balance was computed with the species formulas the selection
itself guards:
- 51 rows are mass-balanced before and not after;
- 144 rows are not mass-balanced before and are after;
- 3 rows cannot be parsed (polymer formulas).

With charge included:
- 44 rows go from balanced to a residual;
- 29 go from a residual to balanced;
- 123 have a residual both ways.

This is a protonation-convention mix, not a one-directional regression. It is outside this GPR
review (§6). The 2026-09-10 STATE.md entry reported the residual counts to the user before the
selection was adopted (155 → 171).

## 3. Evidence and reviewed state per reaction

Tiers: **verified** = opened this session (2026-10-19) and the value read; **recorded** = from
project records, not re-checked; **inferred** = reasoned, no single source.

### 3.1 R765/R766 (transketolase): YALI1E07744g is the transketolase

The 2026-07-06 patch (`patches.py` lines 899–923, commit `7eacb58`) states
"YALI1E07744g → UniProt Q6C6P1 = glycoside hydrolase family 65 (trehalase)". Derivation that
overturns it (verified):

| Check | Result |
| --- | --- |
| NCBI Gene search `YALI1_E07744g` ([esearch](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=YALI1_E07744g%5BAll%20Fields%5D)) | GeneID [2912397](https://www.ncbi.nlm.nih.gov/gene/2912397), aliases YALI1_E07744g, YALI2_E01632g |
| KEGG [yli:2912397](https://www.kegg.jp/entry/yli:2912397) | "(RefSeq) Transketolase", K00615 transketolase [EC:2.2.1.1]; Pfam Transketolase_N, Transket_pyr, Transketolase_C; protein XP_503628 |
| NCBI protein [XP_503628.3](https://www.ncbi.nlm.nih.gov/protein/XP_503628.3) | 694 aa, "uncharacterized protein YALI1_E07744g" |
| UniProt by RefSeq XP_503628 | [Q6C6T4](https://rest.uniprot.org/uniprotkb/Q6C6T4) "Transketolase (EC 2.2.1.1)", YALI0_E06479g, 694 aa, Transketolase family |
| UniProt [Q6C6P1](https://rest.uniprot.org/uniprotkb/Q6C6P1) (the patch's ID) | "YALI0E07744p", 1089 aa, GH65. This is the CLIB122 locus **YALI0**E07744, a different gene from W29 **YALI1**E07744g |
| Experimental workbook putative function (`essentiality_per_gene.tsv`, 2026-09-10 screen) | YALI1_E07744g "similar to … TKL1 transketolase 1 or … TKL2 transketolase 2" |

YALI1 and YALI0 locus numbers are separate schemes. For example, this model's CPA1 is YALI1D09420g,
which is YALI0D07370g in Swiss-Prot Q6C9Y4. The patch's identity claim is a suffix collision, the
"local entry points to a different entity" case.

The retained partner YALI1D02625g is [XP_502319.1](https://www.ncbi.nlm.nih.gov/protein/XP_502319.1):
- 347 aa;
- UniProt [Q6CAJ3](https://rest.uniprot.org/uniprotkb/Q6CAJ3) "YALI0D02277p", annotation score 1;
- only one domain, Transketolase N-terminal (PF00456, residues 30–232);
- no pyrimidine-binding or C-terminal domain, both of which Q6C6T4 has (verified).

That it cannot catalyse transketolase on its own is **inferred** from this domain composition. No
source states its activity either way.

**Reviewed state:**
- YALI1E07744g must be on R765/R766 (verified identity, by homology; no *Y. lipolytica* enzymology
  opened).
- YALI1D02625g's OR membership fails AGENTS.md §5: there is no evidence that it catalyses alone. Its
  removal is a separate decision, needing identity evidence. Recall is not evidence:
  YALI1E07744g is an essentiality false negative in the 2026-09-10 screen, which must not drive the
  change.

### 3.2 R671 (prephenate dehydrogenase, NADP): the removal of YALI1E11370g is correct

| Gene | Identity (verified) |
| --- | --- |
| YALI1E11370g | GeneID [2912526](https://www.kegg.jp/entry/yli:2912526): KEGG "Glutamyl-tRNA(Gln) amidotransferase subunit B", K02434. RefSeq [XP_503726.2](https://www.ncbi.nlm.nih.gov/protein/XP_503726.2) → Swiss-Prot [Q6C6I6](https://rest.uniprot.org/uniprotkb/Q6C6I6) PET112, YALI0E09240g, EC 6.3.5.7, mitochondrial |
| YALI1F23441g | GeneID [2909038](https://www.kegg.jp/entry/yli:2909038): KEGG "Prephenate dehydrogenase [NADP(+)]", K00211 [EC:1.3.1.13]; XP_505545.1, 433 aa |

**Reviewed state:** GPR `YALI1F23441g`, identity verified by homology.
- The selection's equation is mass-balanced with a +1 charge residual, from the stored NADP⁺ charge
  of 0. The builder's equation is short one H. Keep the selected equation and change only the GPR.
- R671 is the only reaction from prephenate to 4-hydroxyphenylpyruvate (static scan of c98961e9).
  After the change, a YALI1F23441g knockout switches R671 off. Its growth effect depends on tyrosine
  uptake in the SD-Leu medium and is **untested**.

### 3.3 R190/R159 (carbamoyl-phosphate synthases): AND vs OR

Gene identities (verified):

| Gene | Identity |
| --- | --- |
| YALI1C33005g | CPA2. KEGG [yli:2909603](https://www.kegg.jp/entry/yli:2909603) "Carbamoyl-phosphate synthase arginine-specific large chain", K01955; XP_502204.1, 1133 aa |
| YALI1D09420g | CPA1. KEGG [yli:2910891](https://www.kegg.jp/entry/yli:2910891) "…arginine-specific small chain", K01956; Swiss-Prot [Q6C9Y4](https://rest.uniprot.org/uniprotkb/Q6C9Y4) CPA1, YALI0D07370g |
| YALI1E11768g | URA2-type. KEGG [yli:2912535](https://www.kegg.jp/entry/yli:2912535) K11541 CPS/ATCase; [Q6C6H5](https://rest.uniprot.org/uniprotkb/Q6C6H5) "Pyrimidine-specific CPS–ATCase"; XP_503737.3, 2186 aa |

Evidence relevant to AGENTS.md §5 (verified from UniProt text and evidence codes):
- **Complex (AND):** *S. cerevisiae* CPA1 [P07258](https://rest.uniprot.org/uniprotkb/P07258) and
  CPA2 [P03965](https://rest.uniprot.org/uniprotkb/P03965) are each a "Heterodimer composed of 2
  chains" (ECO:0000269; PubMed 206535, 5856369). The small chain hydrolyses glutamine; the large
  chain makes carbamoyl phosphate.
  - *Y. lipolytica* Q6C9Y4 carries the same statement by similarity (ECO:0000250, to *N. crassa*
    P22572).
  - The glutamine-dependent R190 therefore needs both: `YALI1C33005g and YALI1D09420g`.
  - The selection's `C33005g or D09420g` has no support.
- **ATCase (R159 "after"):** ATCase is a URA2 domain (Q6C6H5, P07259). The CPA1/CPA2 function texts
  describe carbamoyl-phosphate synthesis only. CPA1/CPA2 as OR members of the ATCase reaction
  contradict their identity.
- **Does URA2 feed the arginine pool (URA2 as OR on R190)?** This depends on the organism.
  - *S. cerevisiae* (P07258/P03965/P07259, ECO:0000305 curator inference): "the 2 synthases appear
    to contribute to the formation of a single cellular pool of carbamoyl phosphate". The same text
    says *S. pombe* and *N. crassa* localize the arginine CPSase in mitochondria and channel each
    synthase's product to its own pathway.
  - P07259 also states the URA2 product "is sequestered and channeled to the ATCase active site"
    (ECO:0000269).
  - For *Y. lipolytica*, Q6C9Y4's mitochondrial-matrix location is only by similarity to
    *N. crassa*.
  - No *Y. lipolytica* experimental source was found. Only UniProt, KEGG and NCBI were searched; no
    literature search.
- **Recorded:** EG-GPR-001 is `accepted` (legacy_human, 2026-07-13) and AGENTS.md says it "remains
  active under schema v1". Its dossier `EGC-legacy-cpa-ura2.json` has:
  - verdict `needs_more_evidence`;
  - confidence `legacy_unverified`;
  - no primary sources;
  - no skeptic pass.
- **Observation, not independent validation:** in `essentiality_per_gene.tsv` (2026-09-10 screen),
  YALI1C33005g (CPA2) is on the consensus essential list. YALI1D09420g (CPA1) and YALI1E11768g are
  not. A strict glutamine-only AND predicts the same phenotype for CPA1 and CPA2. One possible
  reading, **inferred**, is ammonia-dependent CPA2 activity: P03965 says the large chain binds
  "ammonia (free or transferred from glutamine)". The model does not represent that activity.

**Reviewed state:**
- The CPA pair is AND (verified, homology-level).
- R159's GPR must not contain CPA1/CPA2 (verified identity).
- URA2's contribution to the free carbamoyl-phosphate pool is **unclear** for *Y. lipolytica*.

Three representations, for the user to choose (§5, change set C):

| Option | R190 | R159 | Notes |
| --- | --- | --- | --- |
| (i) EG-GPR-001, recommended | `C33005g and D09420g`, product feeds R607 only | URA2 only, net channelled CPS+ATCase (balanced) | Accepted legacy patch. Matches the names, ECs, subsystems and m325 note a56f3649 already carries. Mitochondrial relocation stays deferred, as in the patch |
| (ii) shared pool (*S. cerevisiae*-like) | `(C33005g and D09420g) or E11768g` | ATCase with `E11768g` only | New rule, not an existing patch. Needs its own evidence record |
| (iii) keep selection | 3-gene OR | 3-gene OR on ATCase | Conflicts with gene identity (above). Not recommended |

## 4. Extend the 2026-10-19 step or add a separate one?

The 2026-10-19 module (`selection_gpr_restore.py`) hard-codes its change ID and its four targets.
It also requires each target to have:
- a single-gene GPR;
- a restored equation balanced in mass and charge.

None of the new cases fits unchanged:
- **R671:** GPR-only. Its kept equation has a pre-existing +1 charge residual, so the balance check
  would reject it.
- **R190:** a two-gene AND rule.
- **R159:** an equation change plus a GPR.
- **R765/R766:** need no post-selection restore at all. The final state already contains the real
  enzyme; the defect is the builder patch.

**Recommendation:** generalize the module and add separate change sets.
- Do not add rows to SELECTION-GPR-RESTORE-20261019:
  - c98961e9 was audited against that exact data and flag;
  - extending it would change what `--restore-selection-cleared-gprs` builds;
  - a rejected R190/R159 decision should not block R671.
- Minimal generalization:
  1. a small registry `change_id → (data file, targets)`, with the 2026-10-19 entry unchanged;
  2. `genes` (a list), with each restored GPR's gene set equal to it, in place of the single `gene`;
  3. a GPR-only row may keep its equation unbalanced only when restored and cleared equations are
     identical. The step never edits that equation, and the residual is recorded in the result.
     Any row that edits the equation must still balance.
- One opt-in flag per change set, each requiring a complete selection, as now.
- R765/R766: a builder correction in `patches.py`, not a restore.

The selection record (`metadata_reaction_selection.json`) stays unchanged under this plan.
Converting accepted restores into selection exclusions, as the R_NTP1 amendment did on 2026-09-14,
would make them part of the default build. That remains a separate decision, as the 2026-10-19
report already noted.

## 5. Plan for approval (nothing done yet)

Each item is approved separately.

**A. R671 GPR re-apply (recommended).**
- Data: `model/curation/gpr/selection_overwritten_gpr_restore.json`:
  - change ID `SELECTION-GPR-RESTORE-R671-<date>`;
  - cleared = selection "after", restored = same equation + GPR `YALI1F23441g`;
  - gene identity, evidence and sources as in §3.2.
- Code: the generalization in §4, plus flag `--restore-selection-overwritten-gprs`.
- Test: GPR-only row; refusal on stale state; the 2026-10-19 change set still reproduces exactly.
- Build: a control build reproduces c98961e9; a candidate build differs from c98961e9 only in R671's
  GPR and notes and in YALI1E11370g's reaction links.

**B. R765/R766 builder correction (recommended).**
- Remove R765/R766 from `_GPR_MISANNOTATION_REMOVALS` in `patches.py`. Correct the comments there and
  in `main.py`; R671 stays.
- **Expected model change: none (inferred).** `main.py` has one build path,
  `build_reference_chain`, and it runs the selection unconditionally (line 651). The selection sets
  `YALI1E07744g or YALI1D02625g` either way. Acceptance: control builds reproduce a56f3649 and
  c98961e9 byte for byte; only the build record's per-field selection status changes.
- Records to correct:
  - the memory note;
  - a correction line under the 2026-07-06-derived text in STATE.md;
  - the 2026-10-19 REPORT §6 "verified" wording (owned by the session that wrote it).
- YALI1D02625g removal is not part of B. It needs identity evidence first. An AlphaFold run is
  pre-authorized for identity questions, or literature could serve.

**C. R190/R159 (user's choice of (i), (ii) or (iii)).**
- For (i): data file `selection_cpa_ura2_restore.json`, change ID
  `SELECTION-CPA-URA2-RESTORE-<date>`, its own flag.
  - R159 restored to the balanced net channelled equation with GPR `YALI1E11768g`.
  - R190 GPR `YALI1C33005g and YALI1D09420g`.
  - Guards: R607 still consumes m325; R159/R190 names and ECs are the partitioned ones.
  - The dossier status `needs_more_evidence` is copied into the notes. No claim of new validation.
- (ii) needs a new evidence record and a new rule, not a restore. (iii) needs no code, but the
  a56f3649 name/note inconsistency remains.

**Solves (need separate authorization).** After A and C:
- WT plus single-gene KO of YALI1F23441g and YALI1E11370g (A);
- WT plus KO of YALI1C33005g, YALI1D09420g and YALI1E11768g (C);
- on SD-Leu + PO1f with GLPK: WT plus 5 KOs on the baseline and on the new candidate, 12 LPs.

Expected changes, not yet observed:
- YALI1F23441g KO may now reduce growth (tyrosine route);
- under (i), CPA1 and CPA2 KOs become equivalent;
- URA2 KO closes pyrimidine synthesis from carbamoyl phosphate.

**Workflow per approved item:** `curation-change` skill (TASK.md, data + guarded step + focused
test, new candidate file, whole-model diff, independent audit, REPORT.md, STATE.md entry). No commit
or push without an explicit instruction.

## 6. Limits and items outside this review

- **R2041 direction (affects another line).** The selection reopened a reviewed forward-only
  D-arabinokinase. The patch's stated reason was that reverse flux lets the
  isomerase/phosphatase loop make ATP. Whether that loop exists in a56f3649 is untested and would
  need LPs. Not a GPR question.
- **Stoichiometry convention mix (affects another line):** 51 rows lose mass balance, 144 gain it
  (§2.2). That is a chemistry/protonation question.
- **Builds without the energy candidate.** In those builds the reversible gap-fill directions may
  win. Only the E5 chain (a56f3649) was checked.
- **Literature depth.** Identity checks used NCBI, KEGG and UniProt only. No *Y. lipolytica* CPA,
  URA2, TKL or TYR1 experimental literature was opened.
- **Interpreter.** The main checkout's `.venv` (Python 3.13.5, cobra 0.30.0), loads only.
- **Location.** This report is in worktree `claude/clever-kare-7d2e1a`, uncommitted. The 2026-10-19
  restore it builds on is uncommitted in `claude/orphan-reaction-newest-model-b8ddc8`. Implementing
  §5 belongs in that worktree, or on a branch that contains it.

## 7. Proposed STATE.md entry (not written)

STATE.md was not edited. Three worktrees already hold diverging uncommitted STATE.md heads:
- this one;
- `orphan-reaction-newest-model-b8ddc8`;
- `unruffled-shannon-2baf44`.

Text to add when they are merged:

> ## 2026-10-19: Review of reviewed GPRs overwritten by the metadata selection, investigation only
>
> Delegated read-only review; no model, data or code change. The 2026-09-10 selection's 9 GPR rows
> come from a 2026-09-09 builder (994b09bf) that lacked `remove_misannotated_gprs` and EG-GPR-001;
> no row is new evidence. Verified 2026-10-19 via NCBI/KEGG/UniProt:
> - YALI1E07744g is the transketolase (GeneID 2912397, XP_503628.3, Q6C6T4), not the GH65
>   Q6C6P1, which is YALI0E07744. The 2026-07-06 removal from R765/R766 was an ID-mapping error,
>   masked by the selection.
> - YALI1E11370g is GatB/PET112 (Q6C6I6), so the R671 removal is correct.
> - CPA1/CPA2 form a heterodimer (homology-level), so R190 is AND. Whether URA2 feeds the arginine
>   pool is unresolved for Y. lipolytica.
>
> Plan awaiting approval: R671 re-apply, R765/R766 builder correction, R190/R159 option choice.
> Report: `model/reports/selection_gpr_review_20261019/REPORT.md`.
