# Task C: 13C consistency check of mitochondrial carriers (W29, Worland 2024)

Written 2026-10-13, before any elastic fit or C2 run. Investigation only; no model change.

## Authorization (user, 2026-10-13)

- Prompt `claude_code_prompt_w29_tasks_20261012.md`, header §0–2 plus Task C.
- Plan approved with these decisions:
  - **Base model:** the Task B candidate `E5_coq9_resp_pkg_r132_20261012.xml`, read by path and SHA
    (no merge), with the holistic overlay in memory (diagnostic).
  - **Medium:** W29-type defined medium only. V2 (uracil closed) is primary. V1 (uracil 1000) is run to
    reproduce the recalled elastic-fit totals 6.94 / 1.14 / 2.72 / 0.68; **the 2 % reproduction check
    applies to V1 only**. PO1f + CSM-Leu is skipped for the fit (amino acids make it an invalid
    comparison with the amino-acid-free 13C experiment).
  - **R825:** stop before any curation run. Report the Crc1 direction evidence and, as an in-memory
    diagnostic only (no build), the effect of R825 reversible on:
    - WT growth in W29 V2, and in PO1f + CSM-Leu at uptake 2.414;
    - the full gene benchmark before and after (benchmark rule), listing every changed call;
    - complex I KO/WT with and without PDH ≥ 0.867 × glucose uptake;
    - the elastic fit (V2 primary; V1 to reproduce 2.72 and 0.68).
  - Then wait for authorization.
- Not allowed this round: any model XML build or edit, editing existing media/evidence/benchmark
  files, deleting files, git push/merge/branch switch, cluster jobs.

## Starting state

| Item | Value |
|---|---|
| Worktree / HEAD | `claude/task-c-9ffd14` @ `7dc9511` (clean) |
| Base model | `E5_coq9_resp_pkg_r132_20261012.xml`, SHA-256 `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515`, read from git commit `e985dcd` (branch `claude/task-b-e1ff88`, not merged) |
| Its parent | `E5_coq9_respiratory_package_20261008.xml`, `0ac2addfa0cb6de2d89393e433766fb04325ffe42c2d131eaefaca3bcec29dba` (verified in this worktree) |
| Medium file | `po1f_csm_leu_w29uptake_20261012.csv` from `e985dcd`, expected `6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45` |
| 13C caps (ASSUMPTION) | `w29_13c_upper_caps_20261012.json` from `e985dcd` |
| Overlay | `holistic_overlay.py` from `e985dcd` (diagnostic, prompt header) |
| W29 MFA data | `$IYALI26_RESEARCH_ROOT/reference/w29_worland2024/glucose_central_fluxes.csv` `844615d0…`; README `a8d1c6e6…`; growth `0264fa47…`. Repository ↔ paper link: **inferred, pending user confirmation** |
| Research root | `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research` |
| Python | main checkout `.venv/bin/python` (cobra 0.30.0, Gurobi 13.0.1), as in Task B |

The inputs are extracted from `e985dcd` with `git show` into
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/w29_13c_consistency_20261012/inputs/`, and every SHA is checked
before each run.

## Medium definition (Task B's, unchanged)

W29-type defined medium: the Task A medium plus the PO1f profile (URA3 R612 off, LEU2 GPR on R45), then
adenine R1003 and 11 amino-acid rows closed, R612 restored to (0, 1000), and glucose R1070 fixed at
−2.414. V2 closes uracil R1354; V1 leaves it at 1000.

## Elastic fit (method fixed before running)

- Model: base model + holistic overlay. Glucose fixed at 2.414; `biomass_C` fixed at 0.2663, then
  0.95 × and 0.90 ×.
- Each MFA range i, as a model expression e_i converted to per 100 glucose (100·e/2.414), gets slack
  s⁺ᵢ, s⁻ᵢ ≥ 0 with e_i − s⁺ᵢ ≤ UBᵢ and e_i + s⁻ᵢ ≥ LBᵢ.
- Objective: minimise Σ (s⁺ᵢ + s⁻ᵢ) / (UBᵢ − LBᵢ). "Total violation" is this objective.
  A violation is reported as +s⁺ (above UB) or −s⁻ (below LB), per 100 glucose.
- **Judgments use FVA on the optimal face** (objective ≤ optimum + 1e-6): the range of each slack and
  each mapped flux. A single slack vector is one optimum, not the answer.
- Setups (V2 and V1 each):
  1. holistic overlay (all 22 ranges);
  2. without the malate-import range;
  3. R825 reversible (in memory: lb −1000), all 22;
  4. R825 reversible, without malate import.

### MFA → model mapping (inferred; my judgment, fixed here)

MFA ranges are the `glucose_LB` / `glucose_UB` columns (95 % CI). Signs follow the MFA equation.

| # | MFA row | Name | Model expression |
|---|---|---|---|
| 1 | R4 net | PGI | R326 |
| 2 | R5 net | PFK (net of FBPase) | R636 − R312 |
| 3 | R8 net | GAPDH | R344 |
| 4 | R10 net | PYK | R694 |
| 5 | R32 | PEPCK | R634 |
| 6 | R11 | G6PDH | R325 |
| 7 | R19 | PDH | R693 |
| 8 | R20 | CS (mito) | R236 |
| 9 | R21 net | total mito IDH | R487 + R1385 |
| 10 | R25 | 2-OG DH | R614 |
| 11 | R26 net | SDH | R740 |
| 12 | R28 net | MDH (mito) | R533 |
| 13 | R43.2 | ICL | R490 |
| 14 | R45 | malate synthase | R536 + R537 |
| 15 | R46 | MDH cyt + perox | R534 + R535 |
| 16 | R37 net | ACS | R90 + R2268 |
| 17 | R18 | pyruvate import | R1311 − R842 |
| 18 | R44 | succinate import | R928 + R929 |
| 19 | R47 | malate import (net) | R807 + R902 + R1386 − R831 |
| 20 | cit transp | citrate export (net) | R817 − R831 − R832 − R833 |
| 21 | R51 net | glycine cleavage (synthesis sign) | −(R357 + R360 + R362) |
| 22 | R31 − R29 net | PC minus malic enzyme | R690 − R538 |

Notes on the mapping:
- **Row 9:** "total IDH" uses the MFA's R21 net (aconitase) interval [12.60, 40.36]. In the MFA network
  ICTmit is made only by R21 and used only by R23 + R24, so R21 net equals total mitochondrial IDH.
  Summing the R23 and R24 intervals ([0, 71.2]) would not be a CI of the sum.
- **Row 21:** glycine cleavage is mapped to the glycine-consuming entry steps (R357, R360, R362), not to
  R361. Task B showed R361 is the shared dihydrolipoamide-DH step that also serves 2-OGDH.
- **Row 22:** no CI of the difference is published. The interval is taken conservatively as
  [LB31 − UB29, UB31 − LB29] = [4.712, 17.010]. R29 is negative in the MFA (carboxylation direction);
  the model's matrix malic enzyme R538 is irreversible (decarboxylation only).
- **Rows 19 and 20** sum every cytosol↔matrix carrier of that metabolite, with the signed coefficient
  of the matrix species (carrier list read from the model; transporters for lactate/pyruvate R842
  included in row 17).
- **Not fitted, reported for information:** R42 net (acetyl-CoA, matrix → cytosol, −50.25 [−61.81,
  −45.20]). The MFA network has a direct acetyl-CoA exchange and no ATP-citrate lyase, so its citrate
  and malate CIs are conditional on that network (C1/C2 hypothesis).

If a recalled total is not reproduced within 2 % on V1, I also re-run with R361 in row 21 and with
the overlay on the parent model (R132 present). That is a diagnostic for the cause, not a substitute
for the stop rule.

## C1 and C2

- **C1:** Smith-Waterman (BLOSUM62, Gotoh 11/1; Task B's implementation) of YALI1F33754g vs Sc Ctp1
  (YBR291C) and YALI1B14382g vs Sc Yhm2 (YMR241W), each against a panel of *S. cerevisiae*
  mitochondrial carriers so the best hit is shown, not only one pair. Characterisation papers opened
  and checked (Kaplan 1995; Castegna 2010). Expected: no model change.
- **C2:** direction evidence for Crc1/R825 (Palmieri 1999 and any *Yarrowia* data); evidence for R412
  (mitochondrial HMG-CoA synthase), R1172 (HMG-CoA transport, no GPR) and R945 (mevalonate exchange,
  no GPR). The model route for mitochondrial acetyl-CoA to leave as mevalonate is R83 → R412 → R1172 →
  R410 → R805 → R945 (read from the model this session).
- **C2 runs (diagnostic, in memory):** overlay + 13C caps (Task B caps file; glycine-cleavage forms
  "R361 looser" to reproduce the earlier runs and "module" as the corrected form) × R825 {as is,
  reversible} × PDH floor {off, on} × {WT, complex I KO}, in V2, V1 and PO1f + CSM-Leu. FVA of the
  acetyl-CoA export routes. The benchmark and controls follow the user's list above.

## Stop conditions

- Any input SHA differs from the value above.
- On V1, a recalled elastic-fit total (6.94 / 1.14 / 2.72 / 0.68) differs by more than 2 %: stop and
  report before interpretation.
- A model change would be needed but its direction evidence is only a pFBA or match-rate argument.
- In any case, stop before a curation run (user decision).

## Deliverables

TASK.md, REPORT.md (claims tagged verified / inferred / recalled), CSV tables, maps (PNG), scripts,
an independent read-only audit, a `model/STATE.md` entry. Large outputs go to
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/w29_13c_consistency_20261012/`. No commit or push.

## Amendment 1 (user, 2026-10-13, after the stop fired)

- **Stop outcome.** On V1 the pre-registered fit gave 5.465 / 2.070 / 1.212 / 0.679 against the
  recalled 6.94 / 1.14 / 2.72 / 0.68. All four recalled totals reproduce within 2 % only with the
  parent model `0ac2addf` (R132 present) and glycine cleavage mapped as +R361
  ([reproduction_reconstruction.csv](reproduction_reconstruction.csv)).
- **User decision:** keep the pre-registered setup as primary (Task B candidate, glycine cleavage =
  R357 + R360 + R362). The reconstruction is reported only as the explanation of the recalled
  numbers, in one short paragraph, and is not used as a target.
- **New targets (total violation V2 / V1; a difference above 2 % is a stop again):**

  | Setup | V2 | V1 |
  |---|---|---|
  | R825 import-only, with malate-import range | 6.50 | 5.47 |
  | R825 import-only, without it | 2.54 | 2.06 |
  | R825 reversible, with it | 1.44 | 1.21 |
  | R825 reversible, without it | 0.90 | 0.68 |

- **Corrected Task B step-5 reference values:** the glycine-cleavage cap is R357 + R360 + R362 ≤ 0.002 ×
  glucose uptake. Shuttles-restricted growth ratios (predicted/measured, glucose / glycerol / oleate)
  are 0.99 / 1.09 / 1.02 (V2) and 1.00 / 1.10 / 1.03 (V1); NDH2 ≈ 204 and complex I ≈ 187 per
  100 glucose (V2). The C2 runs therefore use the "module" cap form; the "R361 looser" form is dropped.
- Continue Task C as written: R825 evidence and in-memory effects, then stop for authorization. No
  model change, no push.

## Amendment 2 (user, 2026-10-13, after the report)

- **R825: option (iii).** No XML change. Create a condition file "R825 reversible", labelled
  ASSUMPTION. Its rationale must say there is no organism evidence for acetylcarnitine export on
  glucose, and that it is used only for the 13C fits and Task D.
- **Checks to run in the same run, before using the file:**
  1. **Loop check.** With R825 reversible, at the V2 fit state (1.44) and at WT, run FVA on R825, R1275
     and R1205 plus a loopless FVA (or equivalent). Report whether the cycle carries flux and whether
     the fit total changes when the loop is blocked. **Stop if the fit depends on the loop.**
  2. **Carnitine.** Can the model synthesise carnitine de novo on the W29 and PO1f media? List the
     reactions and GPRs of that route, and add one line to the rationale on whether the shuttle depends
     on it.
- **Instructions for Task D (recorded here, not run in this task):**
  - Put the mitochondrial HMG-CoA route (R412 → R1172) into the closure set as **one part**, not as a
    separate line.
  - Keep R1708 in Task D's part list.
  - The user notes that *S. cerevisiae* HMG-CoA synthase (Erg13) is cytosolic and that mitochondrial
    HMG-CoA synthase is a mammalian ketogenesis feature. Treat this as an *S. cerevisiae* comparison to
    check, not as *Yarrowia* evidence.
