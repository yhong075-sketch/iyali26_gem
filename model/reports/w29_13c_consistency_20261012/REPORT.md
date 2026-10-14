# Task C: 13C consistency of mitochondrial carriers (W29, Worland 2024)

**Status (2026-10-13).** Investigation only; no model was built or edited. Stopped before any curation
run, as the user decided. **The user then chose R825 option (iii)** (Amendment 2): an ASSUMPTION
condition file [r825_reversible_20261013.json](../../conditions/assumptions/r825_reversible_20261013.json),
created after the loop and carnitine checks (§12). Everything is uncommitted on worktree
branch `claude/task-c-9ffd14`. Scope, the pre-registered method and the user's amendment are in
[TASK.md](TASK.md).

**Outcome.**
- **Elastic fit:** consistent with the user's amended targets within 0.5 %. This is a consistency check,
  not an independent reproduction: the targets were set after this run's values were reported (§3.1).
  Making R825 reversible lowers the total violation from 6.52 to 1.43 (V2), and the malate-import
  violation disappears.
- **C1 citrate carriers:** no change proposed. Orthology is clean (reciprocal best hits). One open point:
  the citrate/**malate** mode of R831 rests on homology and database curation, not on yeast
  reconstitution (§4).
- **C2 R825:**
  - **Without a PDH floor, R825's direction barely matters:** complex I KO/WT changes ≤ 0.02 and no
    essentiality call changes at any cutoff.
  - **Under the PDH ≥ 0.867 floor it decides the outcome.** Without 13C caps, export through
    acetylcarnitine or HMG-CoA is required. With caps, only acetylcarnitine export rescues the KO
    (§5.3).
  - **Direction evidence** is mechanistic (an exchange carrier; a near-equilibrium CAT reaction), and
    shown in vivo only for engineered *S. cerevisiae*. *Y. lipolytica* data support import. No
    *Yarrowia* data showing export were found.
  - **Recommendation:** option (iii), no XML change now, with a labelled assumption variant (§5.6).

## 1. Inputs and identities

| Item | Value |
|---|---|
| Worktree / HEAD | `claude/task-c-9ffd14` @ `7dc9511` (tree = `2a34475` + merge) |
| Model | Task B candidate `E5_coq9_resp_pkg_r132_20261012.xml`, `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515`, extracted from commit `e985dcd` (branch `claude/task-b-e1ff88`, not merged) |
| Parent (used only in §3.4) | `E5_coq9_respiratory_package_20261008.xml`, `0ac2addfa0cb6de2d89393e433766fb04325ffe42c2d131eaefaca3bcec29dba` |
| Medium | `po1f_csm_leu_w29uptake_20261012.csv` `6eca70fd…` (from `e985dcd`); V1/V2 W29-type defined media built in memory as in Task B |
| 13C caps (ASSUMPTION) | `w29_13c_upper_caps_20261012.json` `df1e33a8…` (from `e985dcd`); only the "module" glycine-cleavage form is used (Amendment 1) |
| Overlay (diagnostic) | [holistic_overlay.py](holistic_overlay.py) `bf6c7ce5…`, byte-identical to Task B's copy |
| W29 MFA | `glucose_central_fluxes.csv` `844615d0…`. The repository ↔ paper link is **inferred, pending user confirmation** |
| PO1f profile / experimental list / mtDNA scope | `35307853…` / `1e887f5a…` / `16b3e9b4…` |
| Solver | Gurobi 13.0.1, 1 thread, Presolve 0, tolerances 1e-7; cobra 0.30.0; Python 3.13.5 (main checkout `.venv`). Maps are drawn (no LP) with the conda Python, which has matplotlib |

Every LP script checks these SHAs before running; `target_check.py` and `make_maps.py` only read this
task's own outputs. `c2_runs.json` records `inputs_unchanged: true`.
Large and raw outputs (gene screens, flux vectors, evidence texts, sequences, logs) are in
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/w29_13c_consistency_20261012/`.

## 2. Method: elastic fit

[elastic_fit.py](elastic_fit.py):
- **Setup:** model + overlay, glucose fixed at 2.414, `biomass_C` fixed at 0.2663 (then 0.95 × and
  0.90 ×).
- **Slacks:** each of the 22 MFA ranges (95 % CI, per 100 glucose) gets slacks s⁺ and s⁻.
- **Objective:** minimise Σ (s⁺ + s⁻) / CI width.
- **Optimal face:** FVA over every slack, every mapped flux and 16 route reactions, at objective
  ≤ optimum + 1e-6.

**Mapping:** [mfa_mapping.csv](mfa_mapping.csv), fixed in TASK.md before the first run; the mapping is
this task's judgment (**inferred**). Main choices:
- glycine cleavage = −(R357 + R360 + R362), not R361;
- total IDH uses the MFA R21-net interval;
- PC − ME uses the conservative interval [4.71, 17.01]. R_OAADCm (matrix pyruvate carboxylation,
  bounds (−1000, 0)) was left out of this row. The audit re-ran the fit with it added and the totals
  did not change;
- malate import and citrate export sum every matrix ↔ cytosol carrier with its signed matrix
  coefficient.

**Not fitted:** the MFA's R42 net (acetyl-CoA, matrix → cytosol, −50.3 [−61.8, −45.2]) is reported for
information only. The MFA network has this direct exchange and no ATP-citrate lyase, so its citrate
and malate intervals are conditional on that network.

## 3. Elastic-fit results

### 3.1 Totals (verified: [elastic_fit.json](elastic_fit.json), [target_check.json](target_check.json))

| Setup | V2 f = 1 | V1 f = 1 | Target V2 / V1 | V2 f = 0.95 / 0.90 |
|---|---:|---:|---|---|
| R825 import only, 22 ranges | 6.516 | 5.465 | 6.50 / 5.47 | 3.266 / 2.415 |
| R825 import only, without malate import | 2.533 | 2.070 | 2.54 / 2.06 | 1.247 / 0.720 |
| R825 reversible, 22 ranges | 1.433 | 1.212 | 1.44 / 1.21 | 0.954 / 0.581 |
| R825 reversible, without malate import | 0.899 | 0.679 | 0.90 / 0.68 | 0.542 / 0.186 |

All eight differ from the targets by ≤ 0.48 %, so the 2 % stop did not fire.

**Caveat:** this is a consistency check, not a reproduction. `elastic_fit.json` was computed before
Amendment 1 and was not re-run, and the targets were set after these values had been reported to the
user. The V1 targets equal this run's V1 values, rounded. The independent audit re-ran the V2 totals
with a separately written LP and got identical values (§11).

### 3.2 Which ranges are violated (V2, f = 1)

FVA on the optimal face pins every non-zero slack to a single value (min = max), so each violation below
is **required** at the optimum, not an artefact of one solution (verified:
[elastic_fit_summary.csv](elastic_fit_summary.csv)). Signed values: + above the MFA upper bound,
− below the lower bound, per 100 glucose.

| Setup | Violations |
|---|---|
| R825 import only | PGI +14.0, PDH −1.8, CS +9.6, IDH +10.4, 2-OGDH +9.4, SDH +10.0, malate synthase −3.1, pyruvate import −24.6, malate import +12.3 |
| … without malate import | CS +24.8, SDH +3.2, pyruvate import −18.4 |
| R825 reversible | SDH +2.7, pyruvate import −18.3 |
| R825 reversible, without malate import | pyruvate import −12.4 |

**Reading (inferred from these values and the map).** With R825 import-only, mitochondrial acetyl-CoA
has three ways to go:
- **Burnt in the TCA cycle:** this gives CS, IDH, 2-OGDH and SDH above their ranges.
- **Exported as citrate:** in the representative point this runs through R831 citrate/malate (10.0),
  plus 8.4 malate through R902/R1386/R807. Malate import (18.4) therefore exceeds its 6.1 bound.
- **Exported as HMG-CoA:** 2.5 through R412 → R1172.

At the optimum (FVA min of acetylcarnitine export = 38.4 when R825 is reversible), the model instead
moves acetyl units out as acetylcarnitine, and citrate and malate fall inside their MFA ranges. ACL
(R1894) carries 0 in both representative points. The citrate that leaves the matrix feeds the cytosolic
glyoxylate shunt (MFA ICL 10–21), not ACL.

**At f = 1, pyruvate import stays 9.4–24.6 below its range in every setup.** At f = 0.95 and 0.90 the
shortfall shrinks to 0–26, and it vanishes for V1 without malate import at 0.90. A candidate cause (**inferred,
untested**): in the MFA, the mitochondrial pyruvate pool also feeds Ala, Val, Ile and Leu synthesis.
The model may make part of these elsewhere, or by other routes.

### 3.3 Growth sensitivity

At 0.95 × and 0.90 × the measured growth, every total falls (table above) but the ranking holds. R825
reversible stays 3.4–4.6 × lower than import-only with the malate-import range included (V2 and V1).

### 3.4 Why the recalled totals 6.94 / 1.14 / 2.72 / 0.68 did not reproduce

The pre-registered fit missed three of them by more than 2 %, which fired the stop.
- **What reproduces them:** a 24-variant grid ([reproduction_reconstruction.csv](reproduction_reconstruction.csv))
  reproduces all four within 2 % only with:
  - the parent model (R132 present);
  - glycine cleavage mapped as +R361 onto the R51-net interval.
  This gives 6.973 / 1.144 / 2.668 / 0.679, with violation patterns close to the recalled ones (PDH
  −13.6, pyruvate import −27.7, malate import +21.2; CS +10.7 with pyruvate import −9.4). The recalled
  patterns come from the prompt (**recalled**); no file records them.
- **Why that setup is wrong:** R361 is the dihydrolipoamide dehydrogenase (E3) step shared with 2-OGDH
  (R614). Mapping glycine cleavage to it capped 2-OGDH by mistake.
- **User decision:** the recalled values explain where the earlier numbers came from and are not used
  as targets. The other single-change variants (IDH as a sum of intervals, PC and ME as separate rows)
  are in [reproduction_diagnostics.csv](reproduction_diagnostics.csv).

## 4. C1: citrate export

**Orthology (verified this session).**
- Method: local Smith-Waterman, Task B's implementation `04a1f73f…`, not cross-checked against BLAST.
- Searched set: whole Pfam PF00153 families (34 reviewed *S. cerevisiae*; 38 *Y. lipolytica* CLIB122).
- Each RefSeq query is 100 % identical over its full length to its CLIB122 entry.
- Results: [orthology_summary.json](orthology_summary.json), [orthology_sw.tsv](orthology_sw.tsv).

| Model gene (YALI0) | Model reactions | Best Sc hit | Next Sc hit | Reciprocal best |
|---|---|---|---|---|
| YALI1F33754g (YALI0F26323g, XP_505902.1) | R831, R832, R833 | Ctp1 / YBR291C, 67.2 %, cov 1.00 | Sfc1, 30.7 % | yes |
| YALI1B14382g (YALI0B10736g, XP_500732.1) | R817 | Yhm2 / YMR241W, 76.9 %, cov 0.97 | Sfc1, 21.9 % | yes |
| YALI1C03415g (YALI0C02431g, XP_501359.1) | R825 | Crc1 / YOR100C, 51.3 %, cov 0.94 | Ymc1, 37.1 % | yes |

**Characterisation.** The bibliography is verified. Each passage was opened by the literature subagent
this session and saved under `artifacts/…/evidence/`; I re-read the key passages from the saved files.
- **Ctp1** — Kaplan et al. 1995, J Biol Chem 270:4108–4114, [doi:10.1074/jbc.270.8.4108](https://doi.org/10.1074/jbc.270.8.4108):
  - The reconstituted carrier has a strict requirement for internal substrate, i.e. antiport
    (verified).
  - Specificity was tested only as inhibition of citrate/citrate exchange. Citrate and isocitrate
    inhibited and PEP partly; "dicarboxylates … were ineffective" (verified).
  - Whether malate was among the tested dicarboxylates: Table 1 is unread.
  - The same excerpt cites Evans et al. 1983: in intact yeast mitochondria, malate is an alternative
    substrate "to a considerably lesser extent" (verified).
  - UniProt P38152 is internally inconsistent (verified). Its FUNCTION line states export of citrate
    "in exchange for malate" and cites Kaplan 1995 as experimental (ECO:0000269). Its CATALYTIC ACTIVITY
    line gives citrate/malate exchange only "by similarity" (ECO:0000250) to P53007. That P53007 is
    human SLC25A1 is recalled; the saved file does not say so.
- **Yhm2** — Castegna et al. 2010, J Biol Chem 285:17359–17370, [doi:10.1074/jbc.M109.097188](https://doi.org/10.1074/jbc.M109.097188),
  [PMC2878499](https://pmc.ncbi.nlm.nih.gov/articles/PMC2878499/):
  - Obligatory citrate/2-oxoglutarate antiport; uniport is negligible (verified).
  - The proposed physiological direction is citrate out, 2-OG in (the authors' inference, verified).

**Conclusion.**
- **Hypothesis supported (inferred from §3.2 and verified carrier mechanisms):** the model's antiport
  stoichiometry is biochemically right for R817 and R833, and the MFA's malate-import interval is
  conditional on its citrate uniport plus direct acetyl-CoA exchange.
- **No model change proposed.**
- **Open point, flagged and not acted on:** R831 citrate/**malate** on YALI1F33754g is supported by
  homology to the mammalian carrier, by weak malate transport in intact yeast mitochondria (Evans 1983,
  as cited by Kaplan), and by a UniProt FUNCTION line that attributes the malate exchange to Kaplan.
  Kaplan's own reconstitution, however, saw dicarboxylates as ineffective inhibitors. This "limits a claim" (the malate coupling) and does not
  block this task. A targeted review could resolve it, starting with Kaplan Table 1 and Evans 1983.
- **For later 13C work:** couple citrate export and malate import (or fit them together) rather than
  using the MFA malate-import interval as an independent constraint.

## 5. C2: direction of R825 and mitochondrial acetyl-unit export

### 5.1 Direction evidence

Opened by the literature subagent this session; key passages re-read by me from the saved files.

| Source | What it shows | Tier |
|---|---|---|
| Palmieri et al. 1999, FEBS Lett 462:472–476, [doi:10.1016/S0014-5793(99)01555-0](https://doi.org/10.1016/S0014-5793(99)01555-0) | Reconstituted Crc1 exchanges carnitine against carnitine, acetyl- and propionylcarnitine. Slow uniport efflux is detected (112 vs 1612 cpm/min). The authors infer the main role is *import* of acetylcarnitine | verified |
| van Rossum et al. 2016, mBio 7:e00520-16, [doi:10.1128/mBio.00520-16](https://doi.org/10.1128/mBio.00520-16) | CAT is "mechanistically and thermodynamically reversible" (ΔG°′ −1.1 kJ/mol). Export of mitochondrial acetyl units in glucose-grown *S. cerevisiae* needed added carnitine, constitutive shuttle genes **and** evolved mutations | verified |
| van Roermund et al. 1999, EMBO J 18:5843, [doi:10.1093/emboj/18.21.5843](https://doi.org/10.1093/emboj/18.21.5843) | Crc1 (CAC) is required for acetylcarnitine oxidation (import) | verified (subagent) |
| Messina et al. 2023, Microb Cell Fact 22:128, [doi:10.1186/s12934-023-02137-8](https://doi.org/10.1186/s12934-023-02137-8) | **YlCrc1 = YALI0C02431** = YALI1C03415g (100 % identity, this task). The deletion fails on acetate, ethanol and oleate but grows on glucose; the authors call it the sole route *into* mitochondria. Overexpression lowered lipid, which they interpret as more re-import. No export or reconstitution data | verified |
| Worland et al. 2024, Metab Eng 85:1–13, [doi:10.1016/j.ymben.2024.06.010](https://doi.org/10.1016/j.ymben.2024.06.010) | Abstract: the carnitine shuttle is a key regulator of cytosolic acetyl-CoA **in oleic acid and glycerol** cultures. Full text, MFA network and Data availability are **unread** (ScienceDirect blocked) | abstract only |

**Assessment.** AGENTS.md requires enzymology, thermodynamics and conditions for a direction change.
- **Enzymology:** yes. Crc1 is an exchanger; its in-vitro direction is set by gradients.
- **Thermodynamics:** CAT is near equilibrium.
- **Conditions:** missing for *Y. lipolytica* on glucose. In vivo export is shown only in engineered,
  evolved *S. cerevisiae*, and the *Yarrowia* evidence points to import. The Worland abstract names glucose
  among the growth substrates, but states the carnitine-shuttle finding only for oleic acid and
  glycerol cultures.
- **Status:** the fit improvement is a match-rate argument and is not used as evidence. So under the
  task's stop rule, the case for an XML direction change is **not established**.

### 5.2 In-memory effect of R825 reversible

[c2_runs.json](c2_runs.json), verified this run.

**WT growth.**

| Medium | as supplied → reversible | with overlay |
|---|---|---|
| W29 V2 (uptake 2.414 fixed) | 0.29458 → 0.29571 (+0.38 %) | 0.29377 → 0.29487 |
| PO1f + CSM-Leu (uptake ≤ 2.414; 2.414 used) | 0.31243 → 0.31363 (+0.39 %) | 0.31175 → 0.31294 |

**Full gene benchmark.** PO1f + CSM-Leu screen context, model as supplied
([c2_benchmark_scores.csv](c2_benchmark_scores.csv)). The counts are identical before and after:

| View | 1 % | 5 % | 10 % | 15 % (TP/FP/FN/TN, unresolved, MCC) |
|---|---|---|---|---|
| benchmark rule | 64/26/247/704 | 70/28/241/702 | 74/33/237/697 | 82/33/229/697, 16, 0.319 |
| "no solution = lethal" (labelled) | 75/31/247/704 | 81/33/241/702 | 85/38/237/697 | 93/38/229/697, 0, 0.331 |

- **Scoring rules:** negatives are model genes absent from the positive list (user-defined class);
  mtDNA genes are out of scope; this is a development reference, not independent validation.
- **Changed calls: none at any cutoff.**
- **Ratio changes:** 132 genes shift by |Δ| ≤ 0.010 ([c2_benchmark_ratio_changes.csv](c2_benchmark_ratio_changes.csv)).
  - The closest to the 15 % line are the 8 Task A threshold-artefact genes. Seven are arginine-pathway
    genes (incl. the ornithine transporter), plus carbonic anhydrase YALI1F28274g; they move
    0.1591 → 0.1585 and stay above.
- **Controls, as supplied → reversible:** complex I 0.758 → 0.757; ATP synthase 0.110 → 0.110;
  complexes III/IV have no solution in both (unresolved).

**Complex I KO/WT with and without the PDH floor** (PDH ≥ 0.866629 × glucose uptake;
[c2_complex_I_grid.csv](c2_complex_I_grid.csv)). Values are KO/WT, with WT growth in brackets.

| V2 | R825 as supplied, no floor | as supplied, floor | reversible, no floor | reversible, floor |
|---|---|---|---|---|
| no overlay, no caps | 0.753 (0.2946) | 0.244 (0.2944) | 0.752 (0.2957) | 0.617 (0.2957) |
| no overlay, 13C caps | 0.820 (0.2650) | **no solution** (0.2576) | 0.800 (0.2767) | 0.660 (0.2767) |
| overlay, no caps | 0.755 (0.2938) | 0.245 (0.2926) | 0.754 (0.2949) | 0.619 (0.2947) |
| overlay, 13C caps | 0.827 (0.2627) | **no solution** (0.2536) | 0.825 (0.2685) | 0.680 (0.2685) |

- **V1** gives the same pattern within 0.004.
- **"No solution" is unresolved, not lethal.** It is a constraint conflict, not a biological death
  call. Proposed reading (**inferred**; no NADH node balance was computed): the floor forces PDH flux
  whose acetyl-CoA, and possibly matrix NADH, cannot be disposed of once complex I is removed and the
  export routes are closed.
- **The floor costs WT** 3.5 % with R825 as supplied (overlay + caps: 0.2627 → 0.2536) and 0 % with R825
  reversible. This bears on Task D's ≤ 2 % acceptance rule for the PDH ratio.
- **Step-5 reference check** (Amendment 1; overlay + module caps, no floor; verified):
  - V2 growth ratio 0.9865 (target 0.99) and V1 1.0004 (1.00), both within 2 %;
  - V2 NDH2 202.6–206.4 (≈ 204) and complex I 181.8–189.5 (≈ 187) per 100 glucose.

### 5.3 Holistic test of the acetyl-export module (rule 1)

[c2_acetyl_module.py](c2_acetyl_module.py): k = 5 parts, all 32 closure sets, WT and KO each, V2 +
overlay ([c2_acetyl_module.csv](c2_acetyl_module.csv)).

**Parts:**
- P1: R825 export (open = reversible; closed = as supplied);
- P2: R1172, HMG-CoA transport;
- P3: R1996, acetoacetate transport;
- P4: R1894, ACL;
- P5: the 8 isoprenoid/sterol exchanges R945, R1106, R1108, R1124, R1235, R1369, R1388, R1390. These
  are grouped because they share one mechanism: secretion of an isoprenoid-pathway intermediate.

**Results.** "Hit" = KO/WT < 0.15 or KO with no solution; WT cost is relative to all parts open.

| Condition | Minimal cut sets | Leave-one-out from the full closed set |
|---|---|---|
| floor, no caps | {P1, P2} → no solution, WT cost 1.0 % | reopening P1 → 0.619; P2 → 0.232; P3, P4, P5 → still no solution |
| floor, 13C caps | {P1} → no solution, WT cost 5.5 % | reopening P1 → 0.680; any other part alone → still no solution |
| no floor (both cap states) | none | every set gives KO/WT 0.75–0.83 |

**Interpretation (verified by these LPs; mechanism from the pFBA trace and FVA).**
- **Without caps**, with R825 import-only, the KO sends essentially all PDH acetyl-CoA out as HMG-CoA
  (R83 = R412 = 28.3, i.e. 85.0 of 86.7 acetyl units; CS 0.3) (R83 →
  R412 → R1172 → R410). In the representative point it is dumped as 14-demethyllanosterol through
  exchange R1388 (4.6 per 100 glucose).
- The isoprenoid exchanges are **not required**. With R825 as supplied, closing P5 as well still
  leaves KO/WT 0.228; with everything else open, closing P5 leaves 0.619.
- **R945 (mevalonate exchange) carries 0** in every floor run (FVA [0, 0]). The recalled
  "mevalonate sink at 14.8" is not reproduced in this model and medium.
- **Unexplained, not decisive:** in the floor + KO + R825-reversible solution, FVA pins R1708
  (succinyl-CoA:acetate CoA-transferase) at 83.9. Acetylcarnitine export (≥ 161.8) then exceeds PDH,
  fed by PDC → R997 (no GPR) → matrix ALD → R1708. The audit's probe found that closing R1708 or R997
  moves KO/WT only from 0.619 to 0.608. This belongs to Task D's transhydrogenation module.

### 5.4 Evidence for R412, R1172 and R945

Literature by the subagent, opened this session; I did not re-read these passages. Model facts were
read this session.

- **R412**, mitochondrial HMG-CoA synthase, GPR YALI1F38116g. It is the same gene as cytosolic R411
  (EC 2.3.3.10).
  - No experimental mitochondrial localization was found for yeast HMG-CoA synthase.
  - The only "mitochondrion" annotation of *S. cerevisiae* Erg13 (P54839) is NAS, traced to a review
    (Hu et al. 2017, PMC5574775).
  - Mitochondrial mevalonate pathways exist only as engineered constructs (Lv et al. 2016,
    [PMC5036000](https://pmc.ncbi.nlm.nih.gov/articles/PMC5036000/)).
- **R1172**, HMG-CoA transport matrix ↔ cytosol: no GPR, and no evidence found for HMG-CoA crossing the
  inner membrane.
- **R945:** mevalonate secretion is documented only in engineered strains (Wegner & Avalos 2024;
  Rodriguez et al. 2016). Nothing was found for wild-type *Y. lipolytica*.
- **Verdict:** "no evidence found" is **not** evidence of absence. R412 + R1172 form an unsupported
  matrix acetyl-unit export route that becomes decisive under the PDH floor. Candidate for a separate,
  evidence-led review; **classification: affects another line** (Task D and any future PDH-ratio
  overlay). No change proposed here.

### 5.5 A loop to resolve before any R825 change (verified)

With R825 reversible, R825 (reverse) + R1275 (acetylcarnitine import, no GPR) + R1205 (carnitine
export, no GPR) form a net-zero cycle. FVA bounds of R825/R1275/R1205 reach ±41,425 per 100 glucose,
i.e. the ±1000 bounds. It does not change growth, but any curated reversal of R825 must handle R1275
and R1205; the evidence for those geneless transporters was not reviewed here.

### 5.6 Options for R825 (decision is the user's)

| Option | Evidence basis | Effect |
|---|---|---|
| (i) leave as is | Yarrowia data show import (Messina 2023); export unproven | the PDH floor plus complex I KO stays a constraint conflict (no solution) |
| (ii) make R825 reversible (curation run) | exchanger mechanism + reversible CAT; in vivo export shown only in engineered *S. cerevisiae* | WT +0.4 %, no call changes; fit 6.52 → 1.43; must also address the R1275/R1205 loop (§5.5) |
| **(iii) no XML change; labelled ASSUMPTION condition file "R825 reversible"** (recommended; **chosen by the user**, §12) | keeps the model at the evidence level; makes the diagnostic explicit | used for 13C fits and Task D's PDH-ratio tests only |

**Why (iii):** the case for (ii) lacks the "conditions" leg for *Y. lipolytica* on glucose, and the
improvement it brings is a match-rate argument. The prompt's stop rule applies.

## 6. Maps

Values are per 100 glucose from pFBA vectors (one optimum per panel). Judgments rest on the FVA in the
JSON files. The edge values are in [map_edges.csv](map_edges.csv).
- [map_c2_acetyl_routes.png](map_c2_acetyl_routes.png):
  - V2 + overlay + PDH floor;
  - rows: no caps / 13C caps;
  - columns: WT as supplied, WT R825 reversible, complex I KO as supplied, complex I KO reversible.
  - Shows PDH, CS, citrate export, ACL, cytosolic and matrix CAT, R825, thiolase, HMG-CoA synthase,
    R1172, HMGR, isoprenoid exchanges, PDC/ALD/ACS, R997, and R1708.
  - The capped KO-as-supplied panel is marked NO SOLUTION.
- [map_c1_citrate_malate_fit.png](map_c1_citrate_malate_fit.png): the elastic-fit representative
  points (V2, 22 ranges), R825 as supplied vs reversible. Each citrate and malate carrier is shown
  separately, with the MFA intervals.
- **Not on the maps:** the peroxisomal acetylcarnitine source that feeds cytosolic CAT (R191) in the
  as-supplied fit panel. R191 runs acetylcarnitine → acetyl-CoA (11.3) while R825 is 0; the source is
  not traced here.

## 7. Genes named (AGENTS.md §4)

| Gene | Verified name | Function | Evidence level |
|---|---|---|---|
| YALI1C03415g (YALI0C02431g) | **YlCrc1** (Messina 2023) | mitochondrial carnitine/acetylcarnitine carrier | *Y. lipolytica* deletion phenotype and mitochondrial RedStar2 localization (Messina 2023); no reconstitution; orthology to Sc Crc1 (this task) |
| YALI1F33754g (YALI0F26323g) | none verified | tricarboxylate (citrate) carrier, Ctp1 ortholog | homology only (reciprocal best hit, 67 %); Sc Ctp1 reconstituted (Kaplan 1995) |
| YALI1B14382g (YALI0B10736g) | none verified | citrate/2-oxoglutarate carrier, Yhm2 ortholog | homology only (77 %); Sc Yhm2 reconstituted (Castegna 2010) |
| YALI1F38116g | none verified | HMG-CoA synthase (EC 2.3.3.10 by model annotation), GPR of R411 and R412 | model GPR only; no localization evidence |

No benchmark call changed, so there are no further genes to name.

## 8. What is proven, observed or untested

**Provable / static:**
- The slack formulation is an LP. Its optimum is unique in value but not in vector; hence the FVA.
- R825 + R1275 + R1205 is a net-zero cycle by stoichiometry.

**Observed in runs (verified this session):** every number in §§3 and 5.2–5.3. Inputs are SHA-checked.

**Inferred:**
- the mapping of MFA rows to reactions;
- the reading of §3.2;
- the pyruvate-import cause;
- that HMG-CoA export is "unsupported" (absence of found evidence).

**Not done / untested:**
- the Worland full text and MFA network description;
- Kaplan Table 1;
- glycerol and oleate fits;
- reverse-direction BLAST cross-check of the SW implementation;
- evidence for R1275, R1205, R1996 and R1388;
- the peroxisomal acetylcarnitine source.

**Possible improvements, each testable:**
- **Couple citrate and malate.** Fit citrate export and malate import as one joint statistic using the
  MFA covariance; this needs the MFA's covariance or parameter samples. Test: refit, then compare the
  totals and the violation set.
- **Penalise unsupported routes.** Rerun the fit with a small penalty on flux through no-GPR
  transporters (R1172, R1275, R1205, R1996). Test: whether the R825 conclusion survives.
- **Compute the full face.** Use a lexicographic second objective (minimise Σ|v|) and report its vector
  next to the FVA. The maps already do this.

## 9. Reproduce

From the repository root (`IYALI26_RESEARCH_ROOT=/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`,
`PY=<main checkout>/.venv/bin/python`). Each LP script refuses to overwrite its outputs.
`make_maps.py` overwrites its PNGs and `map_edges.csv`; a re-run with a log
(`logs/make_maps_rerun_with_log.log`) produced byte-identical files. The earlier
`logs/make_maps.log` records a failed attempt (missing `memote`); the successful first run had no log.

```
git show e985dcd:<path> > $IYALI26_RESEARCH_ROOT/artifacts/tasks/w29_13c_consistency_20261012/inputs/<file>   # 4 inputs, see §1
$PY -B model/reports/w29_13c_consistency_20261012/elastic_fit.py
$PY -B model/reports/w29_13c_consistency_20261012/target_check.py
$PY -B model/reports/w29_13c_consistency_20261012/reproduction_diagnostics.py
$PY -B model/reports/w29_13c_consistency_20261012/reproduction_reconstruction.py
$PY -B model/reports/w29_13c_consistency_20261012/orthology_sw.py
$PY -B model/reports/w29_13c_consistency_20261012/c2_runs.py
$PY -B model/reports/w29_13c_consistency_20261012/c2_acetyl_module.py
$PY -B model/reports/w29_13c_consistency_20261012/map_vectors.py
python -B model/reports/w29_13c_consistency_20261012/make_maps.py        # conda Python with matplotlib; no LP
```

Sequences were fetched on 2026-10-14 00:22 UTC:
- NCBI efetch: XP_505902.1, XP_500732.1, XP_501359.1;
- UniProt stream: `organism_id:559292 AND xref:pfam-PF00153 AND reviewed:true`, and
  `proteome:UP000001300 AND xref:pfam-PF00153`.

Their SHAs are in [orthology_summary.json](orthology_summary.json).

## 10. Deviations and disclosures

- **Added beyond the pre-registration:** the module test (§5.3) and the reproduction grids (§3.4).
  Both are diagnostic and follow prompt rule 1 and the stop handling.
- **Unneeded deletion:** I removed `map_edges.csv`, a file this session had generated a few minutes
  earlier, before redrawing the maps. The scope excludes deleting files. It was this task's own
  regenerated output and is now rewritten by `make_maps.py` (open mode `w`); nothing else was deleted.
- **Task B step-5 reference values:** the user corrected them in Amendment 1, and they are reproduced
  in §5.2. Task B's REPORT.md is on another branch and was **not** edited.

## 11. Independent audit

A fresh read-only reviewer audited the work on 2026-10-13. It had no conversation history and was given
only the files. It wrote nothing to the repo or the research workspace; its scratch scripts ran in a
temp directory.

**Overall: SUPPORTED, with minor corrections.** All corrections are applied above.

| Claim | Verdict |
|---|---|
| A1 input SHAs; inputs equal the `e985dcd` blobs | SUPPORTED |
| A2 fit implementation and mapping vs TASK.md; carrier signs vs stoichiometry | SUPPORTED. Gap: R_OAADCm left out of row 22; re-run with it gives unchanged totals |
| A3 totals and target check | SUPPORTED. The "reproduction" wording was circular; corrected |
| A4 violations required on the optimal face | SUPPORTED |
| A5 reconstruction grid | SUPPORTED |
| A6 independent re-run of the V2 totals (own LP) | SUPPORTED, identical: 6.515595546190104 and 1.433071041522241 |
| A7 orthology | SUPPORTED |
| A8 literature claims vs saved evidence | PARTIALLY SUPPORTED → corrected: Worland abstract does name glucose; UniProt P38152 FUNCTION line cites Kaplan as experimental; Evans 1983 malate side added |
| A9 C2 growth, benchmark, grid, step-5 numbers | SUPPORTED. "Task A threshold-artefact genes" is recalled |
| A10 module test; R945 FVA [0, 0] in all 28 floor entries | SUPPORTED. Wording corrected (HMG-CoA "or" holds only without caps; P5 sentence) |
| A11 R825/R1275/R1205 loop is net zero by stoichiometry | SUPPORTED (FVA bound 41,425 = 1000 × 100 / 2.414) |
| A12 maps: values on every edge, one panel per condition | SUPPORTED |
| A13 scope: only `model/STATE.md` modified; no XML/conditions/curation diffs | SUPPORTED. The disclosed `map_edges.csv` deletion is noted |
| A14 evidence tiers | PARTIALLY SUPPORTED → corrected (§3.1 caveat, NADH reading tagged inferred, recalled patterns tagged) |

The audit found nothing that changes the main conclusions (no C1 carrier change; option (iii) for
R825). Its scratch scripts are in the session scratchpad (`audit.ti66/`), not in the repository.

## 12. Amendment 2: R825 ASSUMPTION file, loop check, carnitine

**User decision:** option (iii). No XML change; a labelled ASSUMPTION condition file, used only for
the 13C fits and Task D. The checks below ran first, all verified this session:
[r825_loop_carnitine_check.py](r825_loop_carnitine_check.py) →
[r825_loop_carnitine_check.json](r825_loop_carnitine_check.json); 113 LPs plus one MILP.

### 12.1 Loop check (R825 reversible; R825 reverse + R1275 + R1205 is net zero)

| State | Standard FVA: R825 / R1275 / R1205 | Loopless FVA (CycleFreeFlux) | Cycle flux, pFBA point | Fit total or growth with R1275, R1205 or both blocked |
|---|---|---|---|---|
| V2 fit (overlay, 22 ranges), total 1.433071 | [−41425, −38.4] / [0, 41387] / [0, 41387] | R825 [−46.0, −38.4]; R1275 = R1205 = 0 | 0 | 1.433071 in all three |
| V2 WT, overlay (0.29487) | [−41425, −10.6] / [0, 41414] / [0, 41414] | R825 [−19.6, −10.6]; R1275 = R1205 = 0 | 0 | unchanged |
| V2 WT, as supplied (0.29571) | [−41425, −10.7] / … | R825 [−50.7, −10.7]; R1275 = R1205 = 0 | 0 | unchanged |
| PO1f + CSM-Leu WT (0.31363) | [−41425, −11.3] / … | R825 [−48.0, −11.3]; R1275 = R1205 = 0 | 0 | unchanged |

FVA values are per 100 glucose.
- **Full-model loopless MILP** of the V2 fit (cobra `add_loopless`, every internal cycle removed):
  optimal in 247 s, total 1.433071, the same value.
- **Task D use:** complex I KO/WT under the PDH floor (V2, overlay, 13C caps) is 0.680 with the loop
  open and with R1275, R1205 or both blocked.
- **Verdict (verified): the loop is open but carries no flux needed for any result.** The ±41,425 FVA
  ranges are the loop alone. Acetylcarnitine export is genuinely required at the fit state (≥ 38.4)
  and used at WT (≥ 10.6–11.3). The stop condition did not fire.

### 12.2 Carnitine synthesis

- **Max demand is 0 (verified):** every carnitine species (m369 cytosol, m374 matrix, m371 peroxisome,
  m1069 extracellular) has maximum demand 0. That holds on V2, V1 and PO1f + CSM-Leu, with R825 as
  supplied or reversible.
- **No route exists (verified):** the model has **no carnitine biosynthesis route**. Every reaction that
  touches a carnitine species is one of:
  - transport: R1031, R1032, R1205, R1275, R1002, R825, and R873 (YALI1C10232g);
  - carnitine acetyltransferase: R191 (YALI1F28064g), R192 and R193 (YALI1B13948g);
  - the exchange R943 (no GPR), which has lower bound 0 in all three media, so there is no uptake.
- No trimethyllysine or γ-butyrobetaine species exists; the name scan found only carnitine,
  acetylcarnitine and the unrelated trimethylenediamine. Carnitine is not in biomass.
- **Consequence (inferred from the stoichiometry; verified that the runs work with zero carnitine
  supply):** the shuttle runs as a conserved carnitine pool and does not depend on any de novo route
  in the model.
- **Biology:** not settled here. *S. cerevisiae* cannot synthesise carnitine (van Rossum 2016,
  verified). *Y. lipolytica* growth without added carnitine suggests endogenous synthesis (inferred
  from Messina 2023 media; not stated by the authors).

### 12.3 Condition file

- **File:** [model/conditions/assumptions/r825_reversible_20261013.json](../../conditions/assumptions/r825_reversible_20261013.json),
  SHA-256 `e45cab424c30ee9aec00823cb67392742ac6c36ba7091664249765a2533efcaa`.
- **Schema:** the same as `proton_leak_min_8.json`.
- **Contents:** reaction R825, `lower_bound` −1000, original bounds (0, 1000), GPR YALI1C03415g.
- **Rationale:** no organism evidence for acetylcarnitine export on glucose; use only for the 13C fits
  and Task D.
- **Also recorded:** a carnitine-dependence line, the loop-check result and the evidence list, with
  tiers.

[apply_r825_assumption_check.py](apply_r825_assumption_check.py) applies it from the file and
reproduces the in-memory values to 1e-6 ([apply_r825_assumption_check.json](apply_r825_assumption_check.json)):

| Check | Value |
|---|---|
| V2 fit | 1.433071 |
| V1 fit | 1.211624 |
| V2 WT (no overlay) | 0.295706 |

There is no platform loader for assumption files; task scripts apply them explicitly.

### 12.4 Handoff to Task D (user instructions; not run here)

- **HMG-CoA route:** put the mitochondrial HMG-CoA route (R412 mitochondrial HMG-CoA synthase, GPR
  YALI1F38116g, + R1172 HMG-CoA transport, no GPR) into the closure set as **one part**.
- **R1708:** keep succinyl-CoA:acetate CoA-transferase in the part list (§5.3 note).
- **Comparison to check in Task D:** *S. cerevisiae* Erg13 is cytosolic, and mitochondrial HMG-CoA
  synthase is a mammalian ketogenesis feature (user statement). This is an *S. cerevisiae*
  comparison, **not *Yarrowia* evidence**.
  - What this task found: Erg13 GO annotations are cytosol (RCA, computational), nucleus (HDA) and
    mitochondrion (NAS, review-level).
  - No experimental localization was opened, so the comparison still needs a primary source.
- **R825:** use the ASSUMPTION file for Task D's PDH-ratio tests and label every result that uses it.
