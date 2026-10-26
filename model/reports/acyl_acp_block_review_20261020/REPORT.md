# R2216–R2241 acyl-ACP export block and mtFAS chemistry: review and proposal

§1–4 are the review (no edits, no LPs). §5–6 record the approved implementation and its audit.
Model: `fdd00ac3` (combined + lipoylation v2). Scope: [TASK.md](TASK.md).
Inventory of all 40 reactions and 56 species: [evidence/block_inventory.json](evidence/block_inventory.json).

## 1. What the block is

| Group | Reactions | Bounds | GPR |
| --- | --- | --- | --- |
| mtFAS C4 → C8 | R78, R539, R1394–R1405 | KS/KR/ER forward; dehydratases and R78/R539 reversible | single genes or complexes |
| mtFAS elongation C8 → C18 (lumped cycles) | R2218, R2221, R2224, R2230; desaturating R2227, R2233 | forward | YALI1D18037g and YALI1D32594g and YALI1F38317g and YALI1C27933g |
| Mitochondrion ↔ cytosol acyl-ACP transport | R2216, R2219, R2222, R2225, R2228, R2231, R2234, R2236, R2238, R2240 | **reversible** | none |
| Cytosolic acyl-ACP hydrolysis to free fatty acid + ACP | R2217, R2220, R2223 forward; R2226, R2229, R2232, R2235, R2237, R2239, R2241 **reversible** | as listed | YALI1B19844g and YALI1B25427g (FAS1, FAS2) |

All 40 reactions are unchanged from `model/source/iyali26.xml`, except R539, whose GPR was assigned
on 2026-09.

## 2. Findings

### 2.1 Provenance: copied from *S. cerevisiae* iMM904/iND750, with deviations (verified, [bigg_iMM904.json](evidence/bigg_iMM904.json))

In iMM904:
- **Elongation** (`FAS100ACPm`, `FAS141ACPm`, `FAS161ACPm`) is forward-only, with a 6-gene AND rule
  (ACP1, CEM1, HTD2, MCT1, ETR1, OAR1).
- **Acyl-ACP transports:** iMM904 has 8 (`FA120ACPtm`, `FA140ACPtm`, `FA141ACPtm`, `FA160ACPtm`,
  `FA161ACPtm`, `FA180ACPtm`, `FA181ACPtm`, `FA182ACPtm`). All are **export-only** ([0, ∞)) and have
  no gene. iMM904 has **no** C8 or C10 transport (`FA80ACPtm`/`FA100ACPtm` not found).
- **The hydrolase** `FA180ACPH` is reversible, with FAS1 and FAS2.

iYali26 deviates from that source:

| Deviation | Detail |
| --- | --- |
| Transports reversible | All 10 are −1000…1000 (iMM904: export-only). |
| Extra transports | R2216 (C8) and R2219 (C10) have no iMM904 counterpart. |
| Wrong identifier | All 10 transports carry `FA120ACPtm` (the C12 one). |
| Isomer mismatch | R2239 hydrolyses cis-octadec-**11**-enoyl-ACP to "oleate" (cis-**9**). |
| Name mismatch | R2240/R2241 "octadecynoyl" (triple bond) → "octadecadienoate" (two double bonds); same formula, inconsistent names. |

### 2.2 The chemistry is broken by placeholder formulas (verified)

Most acyl-ACP species carry a placeholder formula:
- **`ACP`**, which cobra reads as the elements A, C and P: decanoyl-, dodecanoyl-, myristoyl-,
  palmitoyl-, octadecanoyl-, octadecynoyl- and cytosolic octanoyl-ACP.
- **`C141` / `C161`**, read as 141 or 161 carbons: the tetradecenoyl- and hexadecenoyl-ACPs.
- **`C4H5O2S*` / `C6H11OS*`**, with an invalid character.
- **No formula at all:** 9 C4–C8 intermediates, the cytosolic `ACP`, cis-octadec-11-enoyl-ACP and
  octadecadienoate.

These placeholders cause the "C+138", "A −1" and "P −1" residuals.

The proposed fills use the convention already used for acetyl-ACP (`C2H3OSR`), malonyl-ACP
(`C3H3O3SR`), octanoyl-ACP (`C8H15OSR`) and ACP1 (`HSR`): acyl group + S + R.

| Species (mitochondrion and cytosol copies) | Proposed formula |
| --- | --- |
| acetoacetyl / 3-hydroxybutanoyl / but-2-enoyl / butanoyl | C4H5O2SR / C4H7O2SR / C4H5OSR / C4H7OSR |
| 3-oxohexanoyl / 3-hydroxyhexanoyl / hex-2-enoyl / hexanoyl | C6H9O2SR / C6H11O2SR / C6H9OSR / C6H11OSR |
| 3-oxooctanoyl / 3-hydroxyoctanoyl / oct-2-enoyl / octanoyl (cytosol) | C8H13O2SR / C8H15O2SR / C8H13OSR / C8H15OSR |
| decanoyl / dodecanoyl / myristoyl / tetradec-7-enoyl | C10H19OSR / C12H23OSR / C14H27OSR / C14H25OSR |
| palmitoyl / hexadec-9-enoyl / octadecanoyl / octadec-11-enoyl / "octadecynoyl" | C16H31OSR / C16H29OSR / C18H35OSR / C18H33OSR / C18H31OSR |
| cytosolic ACP / octadecadienoate | HSR / C18H32O2 |

That is 32 species; full list in [scripts/acyl_acp_fills.py](scripts/acyl_acp_fills.py).

**In-memory result with these fills** ([evidence/formula_fill_test.txt](evidence/formula_fill_test.txt);
nothing saved). Every remaining residual is one of three systematic kinds:

| Reactions | Residual after fills | Cause | Fix |
| --- | --- | --- | --- |
| R1400–R1402, all 10 transports | none | — | none needed |
| R1394–R1396 (condensation) | H −1, charge −1 | one H⁺ too many (iMM904 ionized stoichiometry copied onto neutral species) | remove 1 H⁺ substrate |
| R2217, R2220, R2223, R2226, R2229, R2232, R2235, R2237, R2239, R2241 (hydrolases) | H +1, charge +1 | H⁺ product with a neutral fatty acid | remove 1 H⁺ product |
| R2218, R2221, R2224, R2230 (elongation) | H −1, charge −3 | one H⁺ too many; −2 from 2 NADPH | 3 → 2 H⁺; −2 remains (convention) |
| R2227, R2233 (desaturating elongation) | H −1, charge −4 | one H⁺ too many; −3 from 3 NADPH | 4 → 3 H⁺; −3 remains (convention) |
| R1397–R1399, R1403–R1405 (KR/ER) | charge −1 | the model's NADP(+) (`m178`) has charge 0 | none here; model-wide `m178` issue |

So with the fills and 19 one-H⁺ removals, every element balances. The only remaining residual is
−1 charge per NADPH, which is the known model-wide NADP(+) convention and is out of scope.

### 2.3 The C14 → C16 cycle (observed by the audit; mechanism statically verified)

The audit measured up to 36.389 through R2230 under SD-Leu + PO1f, in both `164f235f` and the
lipoylation candidate. Static route:
1. cytosolic myristate + cytosolic ACP → myristoyl-ACP (R2226 **reverse**);
2. → mitochondrion (R2225 **import**);
3. elongation on ACP1 (R2230 with R539);
4. → cytosol (R2231);
5. → palmitate + cytosolic ACP (R2232 forward).

Two steps make this possible:
- **the reverse hydrolase:** a thioester made from a free fatty acid without ATP;
- **the import direction** of the transports, which iMM904 does not have.

Thioester hydrolysis is strongly exergonic. eQuilibrator gives ΔrG'° = −29.6 ± 12.6 kJ/mol
(K'eq 1.6e5) and ΔrG'm = −46.8 kJ/mol for palmitoyl-CoA + H₂O → palmitate + CoA, at pH 7.5 and
I = 0.25 M
([search page](https://equilibrator.weizmann.ac.il/search?query=palmitoyl-CoA+%2B+H2O+%3D+palmitate+%2B+CoA)).
That is the analogous phosphopantetheine thioester, verified 2026-10-20. So the reverse direction
is thermodynamically implausible.

### 2.4 Biology

| Point | Source | Status |
| --- | --- | --- |
| *Y. lipolytica* FAS1 β carries acetyl/malonyl transferase, enoyl reductase, dehydratase and S-acyl fatty acid synthase thioesterase (EC 3.1.2.14) activities. Its cytosolic location rests only on GO:0005829 (IEA); the function text does not state it. | UniProt [P34229 FAS1_YARLI](https://rest.uniprot.org/uniprotkb/P34229) | verified (function text) |
| Cytosolic FAS is already modelled separately: R1392/R1393 make palmitoyl-/stearoyl-CoA with FAS1 and FAS2. | model | verified |
| mtFAS supplies octanoyl-ACP for lipoate. | Schonauer 2009; Hermes & Cronan 2013 (see the lipoylation report) | verified |
| Long acyl chains on mitochondrial ACP dock ACPM to complex I in *Y. lipolytica*. | [Angerer 2017](https://pubmed.ncbi.nlm.nih.gov/28802701/) | verified (abstract) |
| No mechanism is known for exporting or importing acyl-ACP across the inner mitochondrial membrane; ACP is a ~10 kDa matrix protein. | — | established, not cited |

**What this means:**
- **Elongation:** mitochondrial elongation to long chains is plausible biologically (acyl-ACPM in
  complex I). The model, however, has no mitochondrial use for long acyl-ACP other than export.
  After P3 the model can no longer run it (§5).
- **The export block:** the transports plus cytosolic hydrolysis have **no direct support**. They
  exist only because iMM904 had them.
- **Hydrolysis:** whether hydrolysis to the free fatty acid follows is unknown.
- **Annotation:** assigning FAS1/FAS2 to cytosolic acyl-ACP hydrolysis has UniProt annotation
  support (EC 3.1.2.14), but the substrate here would be mitochondrial acyl-ACP that cannot reach
  FAS.

### 2.5 Other observations

- **Dead-end products.** Myristoleate (m1921) and octadecadienoate (m1922) appear in one reaction
  each (R2229, R2241), so those hydrolyses can't carry steady flux (static).
- **YALI1C27933g** (in every elongation GPR) has reverse top hits that are FabG
  3-oxoacyl-ACP reductases (36–40 % identity; [evidence](evidence/YALI1C27933g_vs_swissprot.tsv)).
  That is consistent with the ketoacyl-reductase step in the elongation lumps. The screen file's
  note of similarity to yeast YKR009c (Fox2) is weaker.
- **Holo-ACP synthesis.** It is still missing (no PPT2 ortholog at e ≤ 1e-5). BiGG has `ACPSm`
  (mitochondrial holo-ACP synthase) as a template if a gene is found later.
- **Interaction with the lipoylation module.** Its data file records the nine ACP1-reaction
  residuals, and its guard checks them before and after. The formula fills here change those
  residuals, so the lipoylation guard must accept the corrected state (§3, P5).

## 3. Proposal (one opt-in change, for approval)

| # | Change | Evidence | Recommendation |
| --- | --- | --- | --- |
| P1 | Fill the 32 acyl-ACP / ACP / fatty-acid formulas in §2.2 | chemistry: same convention as the existing acetyl-, malonyl- and octanoyl-ACP and ACP1 | **yes** |
| P2 | Remove the one extra H⁺ in 19 reactions (R1394–R1396, 10 hydrolases, 6 elongation lumps) | element balance with P1 | **yes** |
| P3 | Make the 7 reversible cytosolic hydrolases forward-only (R2226, R2229, R2232, R2235, R2237, R2239, R2241), like R2217/R2220/R2223 already are | thermodynamics of thioester hydrolysis (to be cited before implementation, e.g. eQuilibrator) | **yes**: removes the C14 → C16 cycle without a transport-direction change |
| P4 | Make the 10 transports export-only | only iMM904's bounds plus the absence of a known mechanism; AGENTS.md §5 says another model's bounds are not enough | **no** (not needed once P3 holds) |
| P5 | Teach the lipoylation guard to accept the corrected ACP1-user residuals when this change is present (it runs after it) | consistency | **yes**, required by P1/P2 |
| P6 | Correct metadata: per-chain BiGG transport IDs (C8/C10: none); name "octadecynoyl" vs "octadecadienoate" | BiGG | optional; notes only |
| — | Delete the transports/hydrolases | would need evidence that they don't exist; absence of a mechanism is not that | **not proposed** |
| — | The R2239 isomer (vaccenoyl → oleate) | needs evidence of which isomer mtFAS makes | flag only |

**Implementation path**, if approved: curated data + guarded `apply_*` step + opt-in flag +
focused test + control and candidate builds + whole-model diff + independent audit. Same pattern
as today's changes.

**Acceptance:**
- every touched reaction is element-balanced, and any charge residual equals −1 × NADPH;
- the diff is confined to the targets;
- the lipoylation guard still passes;
- pytest is unchanged apart from the new tests.

**LPs needing authorization (≈10, GLPK, SD-Leu + PO1f):**
- WT growth;
- max R2230, R2226 reverse, R2216 and R2218 (the cycle should be 0);
- R1405 still > 0 under lipoyl demand;
- single KOs of the 4 elongation genes and FAS1/FAS2.

## 4. What the statements rest on

- **Verified:**
  - the provenance and iMM904 bounds (BiGG API);
  - the placeholder formulas and the in-memory residual pattern;
  - FAS1_YARLI function (UniProt);
  - the YALI1C27933g reverse hits;
  - the dead-end products.
- **Observed by the audit (not re-run here):** the C14 → C16 flux of 36.389.
- **Thermodynamics:** cited via eQuilibrator before implementing P3 (§2.3).
- **Established, not cited:** acyl-ACP is not transported across the inner membrane.
- **Untested:** the flux effects of P1–P3 (needs the LPs above).

## 5. Implementation (approved 2026-10-20): P1, P2, P3, P5

**Candidate:** `model/candidates/E5_coq9_resp_pkg_r132_combined_lipoylation_acylacp_20261020.xml`
(SHA-256 `b3b9fe22cba4c8ae4408afe61740bfe75e661c383dd47a224965f06469618543`). It is the lipoylation
v2 candidate `fdd00ac3` plus `--acyl-acp-chemistry`. It is not a reference model; it was
committed on `claude/combined-candidate-20261019` on 2026-10-21, not pushed.

**Mechanism:**
- **Data:** `model/curation/chemistry/acyl_acp_chemistry.json`, generated by
  [scripts/make_acyl_acp_spec.py](scripts/make_acyl_acp_spec.py). It records the before and after
  definitions and the residuals of all 38 block reactions, the 32 fills and the users of every
  filled species.
- **Step:** `apply_acyl_acp_chemistry()` in `platform/scripts/gem_annotate/acyl_acp_chemistry.py`.
  - It requires the lipoylation module (R2286 present, ACP1 = `HSR`).
  - It checks every species, user set, definition and residual **before** any edit, and rejects
    mixed or partial states.
  - It re-verifies the after-state, including notes.
  - The flag `--acyl-acp-chemistry` is off by default and runs after `--lipoylation-module`. The
    builder refuses it without that flag.
- **P5:** the lipoylation guard now takes the recorded residual from `acyl_acp_chemistry.json` for
  any ACP1 reaction whose `acyl_acp_chemistry` note equals this change's ID. Any other note value
  raises. This value check was added after the audit, with two negative tests. Rerunning the
  lipoylation step on the combined candidate returns `already_correct`.
- **Tests:** `platform/tests/test_acyl_acp_chemistry.py` (5 tests, 5 subtests, plus P5 negative cases):
  - only the block changes, and elements balance;
  - rerun of both steps is a no-op after an SBML round trip;
  - stripped notes are rejected;
  - five stale or partial inputs fail closed, each checked for its intended reason;
  - a wrong scope is rejected;
  - the flag is opt-in and requires lipoylation.

**Verification:**

| Check | Result | Level |
| --- | --- | --- |
| Control build (lipoylation v2 flags) | `fdd00ac3…`, byte-identical | verified |
| Candidate, `PYTHONHASHSEED` 0 (the file in `model/candidates/`) and 1 (`.scratch`) | Both `b3b9fe22…`. The seeds were set on the command line and are not stored in the build records. | verified |
| Build record | Input and sources unchanged, build complete, both steps `applied` | verified |
| Whole-model diff vs fdd00ac3 ([model_diff.json](model_diff.json)) | 108 items: 32 formulas + 32 notes; 19 equations (one H⁺ each) + 19 notes; 7 hydrolase bounds → [0, 1000]. Transports unchanged; no reaction, species or gene added or removed. | verified |
| Balance | All 38 block reactions are element-balanced. Charge residuals are only KR/ER −1 (6), elongation −2 (4) and desaturating −3 (2), i.e. −1 per NADPH (`m178` convention). | verified |
| Rerun of both steps on the saved candidate | `already_correct` / `already_correct` | verified |
| pytest | Before (lipoylation v2): 217 passed, 11 failed, 14 errors, 31 skipped. After: 222 passed, same failures and errors. Per test: identical plus 5 new passes. Rerun after the P5 value fix: 222 passed, per test identical to the first after-run ([logs/pytest_after_p5fix_results.txt](logs/pytest_after_p5fix_results.txt)). The two gem_review files were excluded in every run. | verified |

**Solves** (12 LPs against the "about 10" authorized; GLPK, SD-Leu + PO1f with the profile's
`disable_reaction` operations only, so a partial profile; all optimal; observed in runs,
[solves.json](solves.json)):

| Test | Baseline fdd00ac3 | Candidate |
| --- | --- | --- |
| WT growth | 1.4307809228 | 1.4307808279 (Δ −9.5e-8, within solver tolerance) |
| Max R2230 (C14 → C16 elongation), growth not required | **36.389** | **0** |
| Min R2225 / max R2231 / max R2216 | — | 0 / 0 / 0. **Mislabelled originally:** the transports are written cytosol → mitochondrion, so these probe R2225 export and R2231/R2216 import. None of them tests the cycle: the audit found min R2225 = 3.5e-15 on the baseline too; max R2231/R2216 were not run on the baseline. |
| Max R1405 (mtFAS, coupled to lipoyl demand) | — | 2.8e-5 |
| KO YALI1C27933g / YALI1D32594g | — | 1.0 / 1.0 (both listed, still missed) |
| KO FAS1 YALI1B19844g / FAS2 YALI1B25427g | not run (inferred lethal: same GPR on R1392/R1393) | ≈0 / ≈0 (via R1392/R1393; listed) |

**What P3 actually does: all 26 reactions R2216–R2241 are now blocked, in every medium.** This is
provable; the audit found it and I verified it statically, but my plan did not state it:
1. Cytosolic ACP (m1920) is produced by the 10 hydrolases, which are now all forward-only, and has
   no consumer. So every hydrolase must carry 0.
2. Each cytosolic acyl-ACP is used only by its transport and its hydrolase, so all 10 transports
   carry 0.
3. Each mitochondrial C10–C18 acyl-ACP can only be elongated further or exported, so the
   elongation steps R2218–R2233 carry 0.

**Consequences:**
- **Two flux routes closed.** Both pre-existing routes depended on the reverse hydrolysis:
  - C14 → C16, max R2230 36.389 → 0;
  - C14 → C16:1, max R2233 8.992 → 0 (audit).
- **Growth unchanged.** Growth and the lipoylation route (R1405) are unaffected under SD-Leu +
  PO1f.
- **Free fatty acids not stranded.** Free-fatty-acid users outside the block keep flux (audit:
  R1119 −97.9, R1120 4.43).
- **Elongation beyond C8 removed in the model.** Mitochondrial elongation past C8 is plausible
  biologically, but the model has no other sink for long mitochondrial acyl-ACP. Reopening it would
  need a supported sink, for example acyl-ACPM bound in complex I, which has no evidence here.
- **Not deleted.** Nothing was removed and no transport direction was changed. The 26 reactions
  remain in the model, blocked.

**User decision (2026-10-21).** The user accepted this consequence after it was disclosed:

> Accepted: keep P3 as built (b3b9fe22). R2216–R2241 blocked in every medium is acceptable; record the acceptance in STATE.md and the report.

The candidate `b3b9fe22` stays as built. It is not a reference model.

## 6. Independent audit

A fresh read-only reviewer did the following:
- re-derived the diff (cobra and raw XML);
- checked each of the 32 formulas chemically;
- probed the guards with six mutated models;
- reran pytest;
- reproduced the baseline and candidate R2230 values;
- opened eQuilibrator, BiGG and UniProt.

| Claim | Verdict |
| --- | --- |
| Control byte identity; seed identity | supported (seeds not recorded in build files) |
| Diff confined to the targets (109 items) | supported |
| Formulas follow acyl+SR; elements balance; charge = −1 per NADPH from `m178` | supported (setting `m178` charge +1 in memory zeroes every residual) |
| Guards fail closed; rerun no-op; flag opt-in; P5 | supported, with minor caveats (fixed: P5 now checks the note value; negative tests added) |
| pytest | supported |
| WT, R2230 36.389 → 0, R1405 > 0 | supported |
| "R2225 import / R2231, R2216 export closed" | **unsupported as labelled**: the directions were reversed. Corrected above; the audit's baseline min R2225 is 3.5e-15. Closure rests on R2230 and the static proof. |
| What P3 does | **verified but undisclosed**: all of R2216–R2241 blocked in every medium. Now stated above. |
| eQuilibrator, BiGG, UniProt, FabG hits | verified. Corrected: 8 iMM904 transports, 36–40 % FabG identity, the FAS location claim. |

**After the audit,** the P5 value check was added and rebuilt in `.scratch`. The control still
gives `fdd00ac3` and the candidate still gives `b3b9fe22`, so the candidate file in
`model/candidates/` is unchanged. The full pytest rerun gives 222 passed, 11 failed, 14 errors and 31
skipped, per test identical to the first after-run.
