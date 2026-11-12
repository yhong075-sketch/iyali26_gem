# R2121 malonyl-CoA bypass: Stage 1 (in memory, 2026-11-10) and Stage 2 (closure built, 2026-11-12)

Scope, authorization, evidence and plan: [TASK.md](TASK.md) (Stage 2 is section 8). Target state:
[before.json](before.json). Sections 1–6 are the Stage 1 results, kept as written. Section 7 is Stage 2, with the
independent audit in [AUDIT.md](AUDIT.md) and the figure in [escher/](escher/README.md). Scripts, raw solves and
sweeps: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r2121_malonyl_coa_20261110/` (`MANIFEST.sha256`). Nothing committed.

## 1. Run

- Model e584151f (`e584151f…39d2`, short-path copy with the same SHA). V1 = the same model with R2121 bounds
  (0, 0), set in memory after the condition is applied. Task A medium (`6eca70fd…`) plus PO1f (`35307853…`).
- Gurobi 13.0.1, FeasibilityTol / OptimalityTol 1e-9, cobra 0.30.0, Python 3.13.5 (root `.venv`).
- Checks: `scripts/stage1_checks.py`. Sweep: `scripts/ko_screen_certified_v1.py`, which is the reference
  `ko_screen_certified.py` (`313d5d01…`) plus one line in `context()` that closes R2121; the diff is in the manifest
  folder. Platform files have the same SHAs as the v4 screen (`fc78ee8c`, `01b55f73`, `46865060`).
- LPs: 70 for the checks (2 WT, 4 pFBA, 58 FVA, 2 energy, 4 controls) and 1,123 for the sweep (1,082 sweep, 1
  maintenance ceiling, 40 certificate). 1,193 in total, against a budget of 1,400. No build.

## 2. Results (observed this run)

| Check | Baseline | V1 (R2121 closed) | Stop condition |
| --- | --- | --- | --- |
| WT growth | 0.31248856100769 (v4 + 3e-15) | 0.31205269538104 (−0.139 %) | not met (limit −5 %) |
| No-uptake ATP maximum (maintenance relaxed) | 0, optimal | 0, optimal | not met |
| R88 (cytosolic ACC), FVA at 99.9 % | [0, 0.0762] | [0.0077, 0.1024] | not met (minimum > 0) |
| R2121, FVA at 99.9 % | [0.0428, 0.2541] | [0, 0] | — |
| R1936 / R2120, FVA | [0.105, 0.261] | [0, 0.350] | — |
| R2119, FVA | [0, 0.219] | [0, 0.350] | — |
| oxPPP R325 / R639, FVA | [0, 0.212] | [0, 0.212] | — |
| Sweep: TP → FN | — | none | not met |
| Sweep: knockout above WT + 1e-6 | — | none | not met |

**Certified sweep** (rule B; mtDNA and the current untargeted list, `YALI0F04114g` and `YALI0E16192g`, out of
scope; model genes not on the consensus list count as non-essential, a user-defined negative class; the list was used
in development, so this is a development reference):

| Cutoff | v4 TP / FN / FP / TN (MCC) | V1 TP / FN / FP / TN (MCC) |
| --- | --- | --- |
| 1 % | 113 / 213 / 36 / 701 (0.3955) | 115 / 211 / 36 / 701 (0.4014) |
| 5 % | 120 / 206 / 38 / 699 (0.4103) | 122 / 204 / 38 / 699 (0.4161) |
| 10 % | 124 / 202 / 43 / 694 (0.4081) | 126 / 200 / 43 / 694 (0.4138) |
| 15 % | 141 / 185 / 43 / 694 (0.4560) | 143 / 183 / 43 / 694 (0.4616) |

Both screens are rescored here under the same scope. The v4 summary used the earlier untargeted list (without
`YALI0E16192g`), so its own 15 % line reads 141 / 185 / 43 / 695. The 20 certified-infeasible knockouts are the
same in both screens, and no fresh re-solve was needed.

**Class changes, at every cutoff (2 genes, both FN → TP):**
- `YALI1E35955g`: no verified name; RefSeq/KEGG "biotin–protein ligase" (K01942); database annotation only. KO/WT
  1.0 → 0.
- `YALI1F28274g`: no verified name; KEGG EC 4.2.1.1, carbonic anhydrase (reactions R1025, R1026, R2202); database
  annotation only. KO/WT 0.159 → 0.

145 knockout ratios move by more than 1e-6, most by about 0.0013–0.0014 (the WT change). The largest ratio shift
without a class change is `YALI1F01981g`, 0.925 → 0.917.

**Controls on V1:**

| Knockout | KO/WT |
| --- | --- |
| YALI1E35955g | 0 |
| YALI1C15991g | 0 |
| YALI1C02550g | 0.999 |
| YALI1E35955g with the R88 rule set to YALI1C15991g alone | 1.000 |

## 3. Mechanism

- **What R2121 was doing (observed, pFBA plus FVA).** In baseline pFBA, R2121's malonyl-CoA (0.186, made with
  NADPH) fed a TAG cycle on glucose:
  1. lipase R1715 releases fatty acids (0.163);
  2. peroxisomal β-oxidation shortens oleoyl-CoA (R1496–R1555, 0.088);
  3. the ER elongases R1406 / R1407 rebuild C16 with malonyl-CoA (0.089 each);
  4. DGAT160_lp makes TAG again (0.130).

  The peroxisomal acetyl-CoA reached the cytosol through the carnitine shuttle (R192 / R1032 / R191, 0.272) and fed
  GlcNAc-6-P synthesis (R561, 0.271). R2121 had an FVA minimum of 0.043, so the baseline optimum could not do
  without it.
- **What V1 does (observed).** The cycle stops: R1715, the peroxisomal β-oxidation, the C12–C16 elongation, the
  catalase and the carnitine shuttle all go to 0.
  - Cytosolic malonyl-CoA (0.0079) now comes from ACC R88, which is required (FVA minimum 0.0077).
  - Cytosolic acetyl-CoA (0.28) runs Asp → β-alanine → malonate semialdehyde → acetyl-CoA + CO2 + NADH
    (R1936, R2120, R2119). FVA shows this route is optional (minimum 0).
  - Total pFBA flux falls from 294.0 to 282.2.
- **Why `YALI1F28274g` flips (inferred from the static structure plus the observed knockout).** Every reaction that
  can produce cytosolic bicarbonate (R1026, R2202, and R1028 from the nucleus, where R1025 uses the same gene)
  depends on this gene. In V1, ACC R88 needs cytosolic bicarbonate for malonyl-CoA. In v4 the knockout kept 0.159
  through R2121, which needs no bicarbonate. This agrees with yeast biology: nce103Δ, the S. cerevisiae carbonic
  anhydrase mutant, needs high CO2 to grow in air ([Amoroso 2005](https://pubmed.ncbi.nlm.nih.gov/15813743/),
  verified abstract).
- **Why `YALI1E35955g` flips (observed control).** The call comes only from the curated R88 rule
  `YALI1C15991g and YALI1E35955g`: with the rule set to the ACC gene alone, the knockout is fully viable. The model
  does not need biotin ligase flux (R175 stays 0). Treating biotinylation as an AND co-requirement is biologically
  reasonable for holo-ACC, but it is an encoding choice, not evidence. The mitochondrial ACC R87 has no such AND
  partner, so the encoding is inconsistent between compartments.

## 4. What this does and does not show

- **Provable from the model structure:** with R2121 closed, ACC R88 is the only source of cytosolic malonyl-CoA
  (R2003 needs malonate, which has no source), and the carbonic-anhydrase gene supplies every route to cytosolic
  bicarbonate.
- **Observed in this run:** everything in section 2.
- **Inferred:** that the TAG cycle and the Asp-derived acetyl-CoA route are pFBA choices, not biology. One optimum
  does not prove the route is unique; FVA covers only the 14 named reactions.
- **Untested:**
  - other media;
  - the built candidate (the Stage 2 build should reproduce V1's LP exactly);
  - whether YALI1C02550g has any malonyl-CoA-forming activity in Y. lipolytica;
  - ATP-citrate lyase as the physiological source of cytosolic acetyl-CoA (it carries 0 in both models).
- **Possible improvements to the method:**
  - Loopless or futile-cycle-free FVA on the TAG / β-oxidation / elongation set would test whether other such
    cycles remain. Test: compare loopless and standard FVA on those reactions in V1.
  - An Escher before/after map of the malonyl-CoA, acetyl-CoA and bicarbonate routes belongs to Stage 2.

## 5. Found, not changed (each affects another line)

1. `biomass_C` takes chitin at 0.868 mmol/gDW. That is about 0.18 g/gDW as GlcNAc units (C8H13NO5) and makes GlcNAc
   synthesis the largest cytosolic acetyl-CoA sink (0.271). Whether that coefficient is right was not checked.
2. In V1, ATP-citrate lyase carries no flux; cytosolic acetyl-CoA comes from the optional Asp → β-alanine route. That
   route is not required, but its capacity is unconstrained (R1936 genes: RefSeq "glutamate decarboxylase" and
   "L-aspartate decarboxylase dtxS4", database names only).
3. The R88 rule (ACC AND biotin ligase) and R87 (ACC only) encode biotinylation differently.

## 6. Next gates (Stage 2 needs the user's approval)

- **Proposed change:** R2121 bounds (0, 0), keeping the reaction and its rule, with a note giving the KEGG evidence.
  It would go through a curated data file, a guarded `apply_*` step, a focused test, a build, a whole-model diff,
  pytest, an independent audit, this report's Stage 2 section, a STATE.md entry and the Escher before/after figure.
- **Decisions for the user:**
  - whether to close R2121 as above;
  - whether to accept the `YALI1E35955g` TP, given that it comes from the AND rule;
  - whether R1936 and the chitin coefficient go to separate lines.
- **Acceptance criteria for the build:**
  - the whole-model diff shows only R2121's bounds and notes;
  - WT = 0.3120526953810 within 1e-9;
  - the sweep reproduces `ko_screen_v1.tsv` class by class;
  - the no-uptake ATP maximum is 0.

## 7. Stage 2: R2121 closed through the builder (2026-11-12)

**Authorization.** The user, 2026-11-12: "Approve Stage 2" (TASK.md section 8). The two open decisions of section 6
were not answered:
- the `YALI1E35955g` call is reported as an effect of the R88 AND rule;
- R1936 and the chitin coefficient are not addressed here.

### 7.1 What changed and how

- **Data:** `model/curation/direction/r2121_malonyl_coa_close.json`, `R2121-MALONYL-COA-CLOSE-20261112`. It holds the
  reviewed e584151f state of R2121 and of R2119, seven evidence rows (source URL, tier) and the seven notes.
  Generated from the model by `stage2/make_curation_json.py`, which refuses note text with characters SBML escapes.
- **Step:** `platform/scripts/gem_annotate/r2121_malonyl_coa_close.py`, `apply_r2121_malonyl_coa_close`.
  - Every check runs before the edit:
    - the curation scope;
    - R2121's name, rule, bounds, stoichiometry, species and annotation;
    - its notes, with conflicting or foreign notes rejected;
    - R2119's name, rule, bounds, stoichiometry, species and annotation.
  - Then it sets R2121 to (0, 0) and adds the notes. Under `with model:`, bounds and notes are restored on exit.
  - Idempotent: a second run reports `already_correct`.
- **Flag:** `--r2121-malonyl-coa-close`, opt-in, run as the last model step in `main.py` (after
  `--bigg-crosswalk-annotations`). Its guard pins the state that e584151f's flag set produces, after the E5 read-back.
  Without that flag set it fails closed at that point.
- **Test:** `platform/tests/test_r2121_malonyl_coa_close.py`, 6 tests and 43 subtests. It covers:
  - 19 stale-input cases (no edit in any of them);
  - scope rejection and missing keys;
  - the effect, idempotence and the SBML round trip;
  - the context restore;
  - the opt-in flag;
  - forwarding of the flag through `build_model`.
- **Code review:** an independent reviewer found 0 blocking and 5 low findings, all fixed and re-checked. The fixes:
  - a corrected `main.py` comment;
  - an R2119 annotation pin;
  - a stricter scope check;
  - more tests;
  - an EC clause in the limits note.

  The re-check found no new issue. Four nits are left unchanged and reported here:
  - the forwarding test hashes the research root (slow when the file cache is cold);
  - 22 model copies make one test take 24 s;
  - a non-dict `retained_reaction` raises AttributeError, not ValueError (it still fails closed);
  - the step itself does not reject note characters that SBML escapes (the generator and the round-trip test do).

### 7.2 Builds (observed)

Commands in `stage2/logs/build_commands.txt`. All builds used the e584151f argv,
`/opt/miniconda3/bin/python` (Python 3.13.5, cobra 0.30.0, libSBML 5.21.0) and `PYTHONHASHSEED`.

| Build | Flag | SHA-256 | Result |
| --- | --- | --- | --- |
| control (seed 0) | without | `e584151fc79475fc2daa59e76c8c84ebfd06b0f3cb1fbd75fbb69cc9f31039d2` | = e584151f byte for byte |
| seed 1 (scratch) | with | `2f35203750d1904aa4e50dd7d8eac494d5e6a15b142ea005c33152cea889af7c` | step `applied` |
| **candidate** (seed 0) | with | `2f35203750d1904aa4e50dd7d8eac494d5e6a15b142ea005c33152cea889af7c` | = seed 1 byte for byte |

Candidate: `model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_bigg_annot_v2_mcnaugh_rxns_v3_r2121_close_20261112.xml`
(+ `.build.json`, `.coq9_genes.tsv`). Not a designated reference. 3 builds against a budget of 3 plus 1 spare. Every
build log shows the inherited libSBML notes error on R304.

### 7.3 Verification (observed; Task A medium + PO1f, Gurobi 1e-9)

| Acceptance criterion (TASK.md §8) | Result |
| --- | --- |
| 1. Control reproduces e584151f | yes, byte for byte |
| 2. Whole-model diff = R2121 bounds + notes only | yes. 2 differences: bounds [−1000, 1000] → [0, 0]; 7 notes added, the 4 existing ones unchanged (`stage2/model_diff_candidate_vs_e584151f.json`). The SBML `reversible` attribute follows from the bounds. |
| 3. Seeds 0 and 1 byte-identical | yes |
| 4. WT and energy | WT 0.3120526953810365, identical to the Stage 1 V1 value. No-uptake ATP maximum 0, optimal, with only `xMAINTENANCE` relaxed (`stage2/built_checks.json`). |
| 5. Sweep reproduces V1 | yes. `ko_screen_built.tsv` is byte-identical to `ko_screen_v1.tsv`: 0 class differences, 0 KO differences above 1e-9. 20 certified-infeasible knockouts, no fresh re-solve. |
| 6. pytest | yes: 463 → 469 passed; only the 6 new tests differ, all passing (7.5) |
| 7. Review and audit | review: no blocking finding. Audit: 8 claims supported, 2 partially, 0 unsupported (AUDIT.md). |

- **Chemistry:** unchanged. R2121 is mass-balanced with the inherited charge +1 (products minus substrates) in both
  models; R2119 is balanced.
- **Scores** (rule B, mtDNA and untargeted out of scope; user-defined negatives; development reference):

  | Cutoff | e584151f | R2121 closed |
  | --- | --- | --- |
  | 1 % | 113 / 213 / 36 / 701 | 115 / 211 / 36 / 701 |
  | 5 % | 120 / 206 / 38 / 699 | 122 / 204 / 38 / 699 |
  | 10 % | 124 / 202 / 43 / 694 | 126 / 200 / 43 / 694 |
  | 15 % | 141 / 185 / 43 / 694 (MCC 0.4560) | 143 / 183 / 43 / 694 (MCC 0.4616) |

- **Class changes:** YALI1E35955g (biotin–protein ligase by annotation, via the R88 AND rule) and YALI1F28274g
  (carbonic anhydrase by annotation), FN → TP at every cutoff. No TP lost, no new FP.
- **LPs this stage:** 1,125 by the main session (sweep 1,123, checks 2) plus 64 for the figure, so 1,189 against the
  1,400 budget. The audit ran 12 more.

### 7.4 Evidence notes from the audit (recorded, no change needed)

- **Gene identity:** KEGG's yli:2909116 entry is the CLIB122 gene (YALI2_C01175g). The link to YALI1C02550g rests on
  the model's gene annotation (ncbigene 2909116) and the W29 audit row (YALI0C01859g, XP_501343.1), whose protein ID
  matches KEGG's.
- **Second KO:** R00740 has two KOs, K14468 and K15017. Neither has a gene in yli or sce (auditor, KEGG REST).
- **R2119 and R00705:** R2119 carries no `kegg.reaction`; it matches R00705 by stoichiometry, with H+ on the other side
  by the model's protonation convention.
- **Yeast-GEM:** it also carries R00706 (`r_4574`, the NADP form of the decarboxylating reaction). TASK.md §3's "only"
  was corrected.

### 7.5 pytest

Full suite, `iyali26_gem/.venv` (Python 3.13.5, cobra 0.30.0), `python -m pytest -rfE --continue-on-collection-errors`
from `platform/`. Both runs used short-path tree copies:
- before: `git archive` of HEAD `154eb64c`;
- after: the same plus every modified and untracked file of this worktree.

Running both in copies kept the worktree free for the code work, and the short paths avoid the long-path errors seen
before.

| Run | passed | failed | errors | skipped | subtests passed |
| --- | ---: | ---: | ---: | ---: | ---: |
| before (HEAD) | 463 | 3 | 18 | 24 | 1,201 |
| after (HEAD + this change) | 469 | 3 | 18 | 24 | 1,244 |

- **Test by test** (junit XML, `stage2/pytest/compare_before_after.json`): the only differences are the 6 new tests,
  all passing. The failure and error sets are identical:
  - 4 literature-layer collection errors (missing `pydantic_ai` / `httpx2`);
  - 14 `test_lp_sn12_candidate` errors;
  - 3 standing failures (dFBA positive-only reference, lipid moiety ledger documentation, vacuole supply).
- No existing test changed outcome. Large logs are in the artifacts folder (`stage2/pytest/`).

### 7.6 Escher figure

[escher/r2121_close_before_after.png](escher/r2121_close_before_after.png) and the viewer
[escher/r2121_close_escher.html](escher/r2121_close_escher.html) show 51 reactions, pFBA before and after (64 LPs).
- 45 mapped reactions carry flux before and 11 after.
- The peroxisomal β-oxidation chain, the ER elongation chain, the TAG lipase and the acetylcarnitine export go to 0.
- R88, R2119 and R1393 start carrying flux.
- Details and checks are in [escher/README.md](escher/README.md); the audit's count correction (46 → 45) is applied.

### 7.7 What this does and does not show; limits; next gates

- **Provable from the model structure:** with R2121 closed, ACC R88 is the only producer of cytosolic malonyl-CoA, and
  every cytosolic bicarbonate source runs through YALI1F28274g.
- **Observed:** sections 7.2–7.6.
- **Inferred:** that the removed TAG cycle and the R2121 route were artefacts and not *Y. lipolytica* biology. No
  enzyme assay of malonyl-CoA formation from malonate semialdehyde exists for *Y. lipolytica*.
- **Untested:**
  - other media;
  - loopless or futile-cycle checks beyond the 14-reaction FVA set;
  - whether the Asp → β-alanine → R2119 route (now the pFBA source of cytosolic acetyl-CoA, optional by FVA) is
    physiological;
  - the chitin coefficient (0.868 mmol/gDW);
  - ATP-citrate lyase R1894 carrying 0.
- **Next gates (each needs the user's decision):**
  - whether to accept the YALI1E35955g TP as an AND-rule effect;
  - R1936 and the chitin coefficient as separate lines.
  - Publishing (commit, push, PR) needs an explicit instruction.
