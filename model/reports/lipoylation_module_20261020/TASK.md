# Lipoylation module, opt-in (2026-10-20)

## Authorization (user, this session)

> Approved: implement the lipoylation module as planned with D1(a) fix m586, D2(a) cysteine donor +
> NADPH + labelled 5'-dA sink, D3 alpha = 1e-5 with sensitivity at 1e-6/1e-5/1e-4, D4 GCV3 AND in LA3
> only; ~25 LPs; leave R2216-R2241 unchanged; no commit.

> Approved: add the ACP1 m179 formula fill (HSR, charge 0) to the lipoylation change, re-check all 11
> ACP1 reactions and list the 9 pre-existing imbalances in the report; continue the plan as approved.

Not authorized:
- commit, push or PR;
- any change to R2216–R2241;
- connecting cytosolic `ACP` (m1920) to ACP1;
- more than about 25 LPs;
- a full single-KO benchmark.

## Baseline

- Branch `claude/combined-candidate-20261019` (PR #8), HEAD `7718002`.
- Baseline model: `model/candidates/E5_coq9_resp_pkg_r132_combined_20261020.xml`,
  SHA-256 `164f235f4fe823181df0cc83082d6058ee9cf59675777adfe93557e9a77583a6`. It is the newest
  candidate, not a reference model.
- State of every target: [before.json](before.json).
- Baseline pytest: the full run on this same HEAD from earlier today, reused. It gave 212 passed,
  11 failed, 14 errors and 31 skipped. `test_gem_review_lit_chain.py` and
  `test_gem_review_lit_extract.py` were not run because `pydantic_ai` is not installed in the
  `.venv`.

## Intended change (behind the opt-in flag `--lipoylation-module`)

| Item | Change |
| --- | --- |
| `m179[C_mi]` ACP1 | formula `None` → `HSR`, charge 0 |
| `m586[C_mi]` lipoyl-H-protein | formula `C12H18N2O4S2R2(C2H2NOR)n` (unparseable) → `C8H13OS2R`, charge 0 |
| R362 | remove the H⁺ substrate (balanced with the new `m586` formula) |
| R364 | remove the H⁺ product (balanced with the new `m586` formula) |
| New species (mitochondrion) | apo-H-protein `HR`; octanoyl-H-protein `C8H15OR`; apo-E2 `HR`; lipoyl-E2 `C8H13OS2R`; L-cysteine; 5'-deoxyadenosine |
| R2286 LA1 (LIP2-type, EC 2.3.1.181) | octanoyl-ACP + apo-H-protein → octanoyl-H-protein + ACP1. GPR YALI1F34729g. |
| R2287 LA2 (LIP5, EC 2.8.1.8, lumped) | octanoyl-H-protein + 2 SAM + 2 L-cysteine + 4 NADPH → dihydrolipoyl-H-protein + 2 Met + 2 5'-dA + 2 L-alanine + 4 NADP⁺. GPR YALI1E13472g. |
| R2288 LA3 (LIP3-type, EC 2.3.1.200) | lipoyl-H-protein + apo-E2 → apo-H-protein + lipoyl-E2. GPR YALI1C04991g **and** YALI1E00574g. |
| R2289 | L-cysteine, cytosol → mitochondrion (ASSUMPTION) |
| R2290 | 5'-deoxyadenosine sink, mitochondrion (ASSUMPTION) |
| biomass_C | −α lipoyl-H-protein, +α apo-H-protein, −α lipoyl-E2, +α apo-E2, with α = 1e-5 mmol/gDW (provisional) |
| New genes | YALI1F34729g, YALI1E13472g, YALI1C04991g |

The NADPH coefficient follows from the balance under the model's own NADP convention
(`NADP(+)` = C21H29N7O17P3, charge 0).

## Acceptance criteria

1. A control build (`164f235f`'s flags, without the new flag) reproduces `164f235f…`.
2. The whole-model diff of the candidate shows only the items above, plus gene links for the new
   and AND-ed genes.
3. All new and edited reactions, R78 and R539 are mass- and charge-balanced.
4. The other nine ACP1 reactions keep exactly the residuals recorded in the data file. R363 keeps
   its pre-existing charge +1.
5. Guard behaviour: rerun is a no-op (also after SBML reload); a stale or partial state fails closed
   without edits; the flag is opt-in.
6. pytest per test equals the baseline plus the new tests.
7. Solves, at most about 25 LPs, GLPK, SD-Leu + PO1f:
   - WT growth, and R1405 > 0;
   - R2216–R2241 still blocked;
   - KOs: LIP2-type, LIP5, LIP3-type, GCV3, LAT1, KGD2, the two ACPs and the six listed mtFAS genes;
   - α at 1e-6, 1e-5 and 1e-4.

**Stop:**
- a control mismatch;
- a non-target diff;
- a balance failure outside the recorded residuals;
- a test regression;
- an unexpected guard rejection.

## Evidence (opened 2026-10-19/20; details in REPORT.md)

- **Yeast lipoylation genetics:** Schonauer et al. 2009, [PMC2749097](https://pmc.ncbi.nlm.nih.gov/articles/PMC2749097/):
  - Lip2 and Lip5 are needed for Lat1, Kgd2 and Gcv3.
  - Lip3 and lipoylated Gcv3 are needed for Lat1 and Kgd2.
  - The order of the steps is not established.
- **Lip3 family:** Hermes & Cronan 2013, [PMID 23960015](https://pubmed.ncbi.nlm.nih.gov/23960015/).
- **Two mitochondrial ACPs in *Y. lipolytica*:** Angerer et al. 2017, [PMID 28802701](https://pubmed.ncbi.nlm.nih.gov/28802701/).
- **Reaction equations:** Rhea [RHEA:17665](https://www.rhea-db.org/rhea/17665) (octanoyltransferase),
  [RHEA:16585](https://www.rhea-db.org/rhea/16585) (lipoyl synthase; the product is the
  dihydrolipoyl-protein), [RHEA:16413](https://www.rhea-db.org/rhea/16413) (lipoyl amidotransferase).
- **Gene identities:** reciprocal BLAST, W29 RefSeq GCF_001761485.1 and Swiss-Prot 2026_03. The
  databases from `$IYALI26_RESEARCH_ROOT/artifacts/tasks/orphan_gene_search_20261019/db` were reused.

## Corrections after the independent audit (2026-10-20)

The text above is kept as written. These corrections apply; the candidate was rebuilt as v2.

- **R2288 mechanism.** The cited Hermes & Cronan 2013 found Lip3 octanoyl-CoA:protein transferase
  activity but **no amidotransfer and no ligase activity** in vitro
  ([PMC3806487](https://pmc.ncbi.nlm.nih.gov/articles/PMC3806487/), full text opened). UniProt
  P47051 is PE3, a putative ligase. R2288 keeps its equation and GPR, but:
  - it is renamed "LIP3-dependent lipoylation of E2 lipoyl domains (lumped; mechanism unresolved)";
  - its EC 2.3.1.200 / RHEA:16413 cross-references are dropped;
  - its tier is lowered to C;
  - the counter-evidence is recorded.

  The net cost per E2 lipoyl (one mtFAS octanoyl and two sulfur insertions) is the same under both
  proposed mechanisms.
- **R2288 notes lost.** An unescaped "&" in R2288's evidence made libsbml drop the whole notes block
  of R2288 in the v1 candidate (`1be380ca…`). The text was fixed. The guard now requires every
  lipoylation note and EC/Rhea annotation to be present after a save and reload, and a test checks
  it.
- **Guard order.** The residuals of the ACP1-user reactions are now checked **before** any edit, not
  only after.
- **4 NADPH in R2287.** This follows from the model's NADP(+) convention (m178 carries the cation
  formula C21H29N7O17P3 with charge 0). With ChEBI charges the chemistry needs 2 NADPH (4 electrons).
  It is kept as written so that the reaction is charge-balanced in the model; the cost is negligible.
- **Evidence wording.**
  - RHEA:16585 describes the auxiliary cluster as the sulfur source (released as 2 H2S + 4 Fe3+).
  - Q6C6C6 is PE3, HAMAP-inferred, and is the same protein, so it is not independent evidence.
  - YALI1F34729g has family-level reciprocity only (yeast Lip2 is not in the reverse top 5). It is
    not on the consensus essential list (Cas9 essential q 0.0023, Cas12a non-essential q 0.68).
- **`m586` name.** It still embedded the old formula; it is now `lipoylprotein_C8H13OS2R`.
- **Acceptance criterion 3:** R2290 is a sink and is unbalanced by construction.
- **Acceptance criterion 7, "R2216–R2241 still blocked".** This is false. R2216 and R2218 stay at 0,
  but the reversible cytosolic thioesterases and transports allow a C14 → C16 elongation cycle
  through R2226/R2225/R2230/R2231/R2232. The audit measured up to 36.389 in **both** models, so it is
  pre-existing. The module leaves R2216–R2241 unchanged.
- **LP accounting:**
  - this change used 26 LPs (25 solves plus 1 prototype);
  - the earlier 5-LP ACP test was a separate request;
  - the audit ran 10.
- **Hash seeds:** candidate and control builds use `PYTHONHASHSEED=0`; the scratch copy uses 1.
