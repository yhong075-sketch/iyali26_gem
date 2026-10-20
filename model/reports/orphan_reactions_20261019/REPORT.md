# Orphan reactions in the newest candidate: inventory, flux evidence and next actions

Model: `E5_coq9_resp_pkg_r132_20261012.xml` (`a56f3649…`), the newest candidate, not a designated
reference. Condition: SD-Leu + PO1f (see [TASK.md](TASK.md)). Read-only: no model, GPR or data change.

## 1. Definition and counts

An orphan reaction is an enzyme-catalysed metabolic reaction with no gene identified in the
organism. Operationally: a reaction with an empty GPR, excluding boundary (exchange, demand,
sink), biomass/pool/maintenance pseudo-reactions, transport and reactions marked spontaneous.

| Gene-less reactions | 840 of 2,306 |
| --- | --- |
| Boundary | 183 |
| Biomass / pool / pseudo | 38 |
| Transport (multi-compartment or transport subsystem) | 533 (47 passive-diffusion-like, 486 carrier-type) |
| Metabolic | 86 |
| of which spontaneous by name or notes | 6 (R735, R736, R1194, R1195, R1375, R1376) |
| **Orphans (core definition)** | **80** |

Two more, R1210 (glutamate 5-semialdehyde → P5C) and R2197 (3-sulfinylpyruvate → pyruvate +
sulfite), are usually described as spontaneous. That is textbook knowledge with no source opened
in this session, so they remain counted until a source is cited. If both are reclassified, the
count is 78. Carrier-type transport is outside the definition; counting it would give 566.
Full list: [gene_less_classification.csv](gene_less_classification.csv).

## 2. Flux evidence (SD-Leu + PO1f, WT μ = 1.4308 h⁻¹)

[flux_evidence.csv](flux_evidence.csv): pFBA, FVA at 0 % and at 99.9 % of optimum (separate
columns), single-reaction KO. Status for the 80 orphans:

| Status | n |
| --- | --- |
| Required at the WT optimum (FVA range at 99.9 % excludes 0) | 5: R19, R321, R1043, R1210, R1901 |
| Active in WT pFBA, not required | 1: R1308 |
| Can carry flux, 0 in WT pFBA | 41 |
| Blocked (FVA 0 % range is [0, 0]) | 33 |

KO/WT is 0 for R19, R321, R1043 and R1210. R1901 KO/WT is 0.9973 and R1308 KO/WT is 0.9995.
Every other orphan has KO/WT 1.0000. Under this condition no orphan except those four is needed
for growth.

## 3. Static model evidence

[orphan_evidence.csv](orphan_evidence.csv) has one row per reaction. Each row lists the EC/KEGG
annotations, element and charge balance, exact duplicates and compartment-free duplicates (with
their GPRs), GPR fields changed by the metadata selection, and dead-end metabolites.

### 3.1 Four orphans had reviewed GPRs that the metadata selection cleared (verified)

The builder assigns reviewed genes before the final `apply_metadata_reaction_selection` step
(`platform/scripts/gem_annotate/main.py` lines 509 and 650). The selection
(`model/curation/selection/metadata_reaction_selection.json`,
`metadata_reaction_fields_preserve_trna_20260910`) then writes `gpr: ""` for:

| Reaction | GPR before the selection | Evidence recorded in the builder |
| --- | --- | --- |
| R612 orotidine-5'-P decarboxylase | YALI1E31685g (URA3) | direct *Y. lipolytica* ura3 genetics, [doi:10.1002/yea.3673](https://onlinelibrary.wiley.com/doi/full/10.1002/yea.3673) (`patches.py` 955–1020) |
| R2274 triolein hydrolysis | YALI1A21372g (LIP2) | lip2Δ and Lip2 triolein activity, [PMC101989](https://pmc.ncbi.nlm.nih.gov/articles/PMC101989/) (`patches.py` ~1152) |
| R1302 pyridoxal kinase | YALI1A08512g | UniProt automatic annotation; the same gene is on R1306 (EC 2.7.1.35) |
| R1303 pyridoxamine kinase | YALI1A08512g | same as R1302 |

Derivation: the `before.gpr` and `after.gpr` fields in the selection JSON, printed this session.
In the model, URA3 YALI1E31685g is present with 0 reactions.

The selection also changed the stoichiometry of two of these reactions:
- R2274: the selection removes 3 H₂O, leaving the reaction short by H6 O3.
- R1303: the selection adds 1 H⁺, giving an imbalance of H +1 and charge +1. Before the
  selection the balance check gives 0.

The selection record says the user chose the metadata reaction fields. It does not say whether
clearing these four GPRs was intended. That decision is the user's.

R612 is closed in PO1f by the strain profile, so restoring its GPR does not change PO1f growth.
R2274 is blocked on SD-Leu because no triolein is supplied. Its GPR matters only under
triolein-carbon conditions.

### 3.2 Exact duplicate pair (verified)

R530 and R1727 have identical stoichiometry: lyso-PC + acyl-CoA → PC + CoA, both in `C_mm`, with
no GPR. Both are blocked. The acyl-CoA `m790[C_mm]` has no formula, so their balance cannot be
judged.

### 3.3 Net-zero internal loop (verified, [loop_check.txt](loop_check.txt))

R1964, R2041 and R2042 reach ±1000 in FVA. The test closed every boundary reaction and relaxed
the one forced flux (xMAINTENANCE) to allow 0. Under those conditions they still carry 1000 in
the cycle

R803 (xylulose kinase, ATP→ADP) → R714 (X5P → Ru5P) → R2042 (Ru5P → arabinose-5-P) →
R2041 reversed (arabinose-5-P + ADP → arabinose + ATP) → R2190 reversed (arabinose → arabitol) →
R1964 reversed (arabitol → xylulose).

The net exchange is zero, so this is a stoichiometric loop, not a net energy source. It inflates
FVA ranges and lets an arabinose kinase run backwards to remake ATP. Which step is
thermodynamically wrong has not been tested.

### 3.4 Compartment copies of gene-bearing reactions (leads only)

These orphans match a reaction that has a GPR, with metabolites compared by name and formula
regardless of compartment:

| Orphan (compartment) | Gene-bearing reaction (compartment) | GPR |
| --- | --- | --- |
| R989 (mi) | R64 (cy) | YALI1D01489g |
| R1125 (mi) | R306 (cy) | YALI1D33889g |
| R1292 (mi) | R615 (cy) | YALI1D24649g |
| R1024 (mi), R1027 (ex) | R1025 (nu), R1026, R2202 | YALI1F28274g |
| R1337 (cy) | R751 (ex) | YALI1C09285g or (6-gene AND) |

A gene on one compartment is not evidence that the enzyme is in another (AGENTS.md section 5).
Each lead needs a targeting check, for example TargetP/MitoFates/DeepLoc as in Task B, or a
localization source.

Two matches are excluded:
- **R321 (cy) and R218 (lp):** R218 is named "cis-prenyltransferase step 02". The names suggest a
  dolichol-type cis elongation, unlike the trans elongation in R321 (inferred from names; the
  stereochemistry is not encoded).
- **R735 and R736:** this pair is spontaneous.

R751's GPR (one gene OR a six-gene AND) is unusual for a phosphatase and is flagged outside this
scope.

### 3.5 Chemistry

Of the 80 orphans:
- **20 involve a metabolite with no formula:** R530, R1342, R1344, R1727, R1769, R1770, R1902,
  R1961, R1962, R1991, R2000, R2021, R2034, R2039, R2072, R2115, R2142, R2178, R2197, R2284.
- **6 have formulas but are mass-imbalanced:**
  - R1040: H +1, charge +1;
  - R1303: H +1, charge +1, introduced by the selection;
  - R1304: H −1;
  - R1337: thiamine formula convention, Cl −3 and more;
  - R2274: H6 O3, introduced by the selection;
  - R2276: missing diphosphate, H4 O7 P2.
- **10 are charge-only imbalanced.** This matches the model-wide protonation-convention issue and
  is not specific to orphans.

Two reactions look chemically doubtful (inferred from the equation, not checked against a source):
- **R2178** "dehydrodolichol-PP → dolichyl-PP" is written with 6.5 O₂ consumed and 5 CO₂
  released. Saturating one double bond is a reduction, so this stoichiometry needs review.
- **R2276** D-limonene synthase is a plant terpene synthase. It is probably an engineering
  addition rather than wild-type metabolism (established, uncited).

### 3.6 Dead ends among blocked orphans (verified, static)

These metabolites have no other reaction that can make or use them under the model's default
bounds:
- R52: 4-hydroxybenzoyl-CoA[mi] has no producer.
- R1988: salicyl alcohol[cy] has no consumer.
- R2021, R2034, R2039: vacuolar H₂O has no producer. This matches the 2026-09-23 vacuolar
  dipeptide diagnosis in STATE.md.
- R2274: triolein[ex] has no producer.

The other blocked orphans are blocked by gaps further away.

## 4. Recommended action per orphan

[triage.csv](triage.csv). Each reaction gets the first matching rule (`scripts/triage.py`). These
are recommendations for review, not decisions.

| # | Group | n | Reactions | Recommended action |
| --- | --- | --- | --- | --- |
| 1 | Reviewed GPR cleared by the selection | 4 | R612, R1302, R1303, R2274 | **User decision.** If the losses were unintended, re-apply the four GPRs after the selection (as R1026, R153 and others already are) and restore the balanced R1303/R2274 stoichiometry through a curation change. Cheapest and best-evidenced fix. |
| 2 | CoQ9 line | 4 | R19, R52, R321, R1043 | Leave to the CoQ9 line; its approval gate stays in force. R19's empty GPR is a recorded decision (`model/curation/evidence/quinone_branch_cleanup.md`). R19, R321 and R1043 are lethal on KO, so CoQ9 biomass currently depends on gene-less steps. R1043 (chorismate → 4-HB + pyruvate) is the bacterial UbiC-type route (established, uncited). |
| 3 | Likely spontaneous | 2 | R1210, R2197 | Cite a source, then mark spontaneous (annotation only, no flux change). R1210 is lethal on KO, so it must stay in the model. |
| 4 | Exact duplicate | 2 | R530, R1727 | Merge into one reaction (needs authorization), fill the acyl-CoA formula, then search an LPCAT gene for the remaining one. |
| 5 | Internal loop | 3 | R1964, R2041, R2042 | Review the direction of R2041 and R1964, and whether *Y. lipolytica* has arabinose-5-P isomerase and arabinose kinase (R2042, R2041) at all. Needs enzymology or thermodynamic evidence, not flux. |
| 6 | Reverse flux at the optimum | 2 | R1308, R1901 | R1901 asparaginase runs backwards (Asp + NH₃ → Asn without ATP) and is "required" only in that direction. Check its reversibility against a ΔG source. Candidate only; KO costs 0.27 %. R1308 (uracil salvage, reversed) is low stakes. |
| 7 | Compartment copy | 6 | R989, R1024, R1027, R1125, R1292, R1337 | Targeting check on the gene-bearing paralogue. If dual-localized, assign the GPR. If not, the copy needs its own gene or a transport explanation. |
| 8 | Chemistry incomplete | 20 | R1040, R1304, R1342, R1344, R1769, R1770, R1902, R1961, R1962, R1991, R2000, R2021, R2034, R2039, R2072, R2115, R2142, R2178, R2276, R2284 | Fill formulas and fix balance first (R2178, R2276 above), then gene work. |
| 9 | Blocked, no chemistry flag | 12 | R356, R970, R990, R1079, R1081, R1095, R1097, R1307, R1324, R1988, R1997, R2007 | Gap analysis (root no-production/no-consumption metabolites) before any gene work. Blocked is not evidence of absence; deletion needs evidence that the reaction does not exist in the organism. |
| 10 | Can carry flux, no other flag | 25 | R23, R333, R517, R721, R984, R1048, R1078, R1080, R1090, R1132, R1154, R1178, R1252, R1273, R1274, R1338, R1347, R1958, R1959, R1973, R2003, R2010, R2013, R2075, R2076 | Gene search: characterised query → BLAST to the W29 proteome → reciprocal Swiss-Prot check, with evidence tiers. Needs a separately authorized computation (downloads the W29 RefSeq proteome and Swiss-Prot). |

Excluded as spontaneous by name: 6. Total 80 + 6 = 86 metabolic gene-less reactions.

**Suggested order:**
1. Group 1: one user decision; code and sources already exist.
2. Groups 4 and 5: small topology fixes that also clean up FVA.
3. Groups 10, 7 and 6: a gene search over 33 reactions.
4. Group 8: chemistry.
5. Group 9: gap analysis.

Groups 2 and 3 run on their own lines.

**What not to do:**
- Do not delete blocked or idle orphans to raise a score.
- Do not copy a GPR across compartments without targeting evidence.
- Do not change a direction because flux uses it.
- Do not rely on KEGG KO → gene links without a reciprocal check.

## 5. What each statement rests on

- **Provable from stoichiometry:** the R1964/R2041/R2042 cycle has zero net exchange; R530 and
  R1727 are the same equation.
- **Statically verified this session:** counts in section 1, balance results, duplicates, dead
  ends, and the four GPRs cleared by the selection.
- **Observed in runs this session (GLPK):** WT 1.4308, flux statuses, KO ratios and the loop
  test. A second run of `loop_check.py` from the copy in this folder gave the same output.
- **Inferred:** R1210/R2197 spontaneity, R2178 chemistry, R2276 being heterologous, the R321/R218
  stereochemistry, R1043 being bacterial, the R1901 direction problem.
- **Untested:** any gene assignment, targeting, literature, other media or carbon sources,
  loopless FVA.

## 6. Method limits and possible improvements

- **One medium.** Statuses hold for SD-Leu + PO1f only. A reaction blocked here may carry flux on
  triolein, glycerol or oleate. A multi-condition FVA (same model, named media) would separate
  "blocked here" from "blocked everywhere".
- **Loop inflation.** The FVA ranges for loop members are inflated. Loopless FVA (CycleFreeFlux
  or an ll-FBA MILP) would give honest ranges. Test it on R1964/R2041/R2042 first; with GLPK the
  MILP may be slow.
- **Name-based matching.** The compartment-free duplicate test keys metabolites on name and
  formula. It can miss renamed copies and match stereo-isomers that share a name and formula
  (R321/R218). MetaNetX/ChEBI IDs on metabolites would make it exact.
- **Spontaneous flag.** The flag uses names and notes. R1210 and R2197 show that it misses
  unflagged spontaneous steps.
