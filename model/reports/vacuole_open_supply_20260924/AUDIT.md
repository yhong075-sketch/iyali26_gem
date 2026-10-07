# Independent audit: E5 vacuolar opening and finite supply

Completed independent audit, 2026-09-24. The pinned raw E5 XML was read independently with Python's standard XML parser. All 114 saved primary-run optimization records were then checked independently using Python's standard library. No optimization, network access, outside-workspace access, or model modification was performed by this auditor. This is an independent check of delivered records and feasible witnesses, not an independent solver rerun or dual-certificate proof of optimality.

The E5 input hash matches the preregistered identity. It contains 2314 reactions and 1877 species. Full incidence for vacuolar water, H⁺, Asp, Glu, Gly, Ala, Pro, and an isoleucine guard row is preserved in `auditor_static.json`; `auditor_check.py` reproduces the extraction.

The independently parsed `E5_vacuole_open.xml` differs in exactly four upper bounds: R1363 and R795 become 0.04; R871 and R876 become 0.01. All lower bounds remain zero. Reversing those four XML bound references and removing their four new parameter elements makes the entire normalized XML tree equal to E5, including every other element, attribute, text, and child order. The candidate has no added diagnostic source. Its SHA is `e4806a1a7fcd69af1b5ec7c2e014890828548dd1e77630db2b3a2a96dbfc9a39`; scope checks are recorded in `auditor_candidate.json`.

## Structural findings

Vacuolar water has exactly five incident reactions:

`v_R1363 − v_R2021 − v_R2029 − v_R2034 − v_R2039 = 0`.

Asp and Glu each have exactly two incident reactions. Their only exits are R871 and R876, respectively. Thus opening the water and ATPase connections alone cannot allow their hydrolyses while those exits remain closed. Gly has four hydrolysis inputs and only R2030 for transport; Ala has R2034 and R2035, and Pro has R2039 and R2040.

The H⁺ row has 22 incident reactions, all included in the audit data. R882 is the only other originally open irreversible H⁺ consumer outside the target transport group. Its vacuolar isoleucine product has only R882 and the closed R883, so its flux remains zero under the scoped candidate changes. The other previously closed H⁺ connections remain closed.

For direct vacuolar supplies, write `D,G,A,P` for the four hydrolysis fluxes, `S = D+G+A+P`, and `e_D,e_G,e_P` for nonnegative external secretion of the corresponding intact dipeptides. With only the four authorized connections opened, the complete rows imply:

`v_R1363 = S`; `v_R871 = D`; `v_R876 = G`; `v_R2030 = −S`; `v_R2035 = −A`; `v_R2040 = −P`.

`2 v_R795 = S + D + G + P + e_D + e_G + e_P`.

Native dipeptide transport equations, despite their “antiport” names, move H⁺ and the peptide in the same direction for Gly-Asp, Gly-Glu, and Gly-Pro; Ala-Gly transport has no proton term. The actual stoichiometry, not the names, is used above. The supply can exceed hydrolysis when intact peptide is secreted, so equal supply and hydrolysis must not be assumed from one source bound alone.

If all four hydrolyses simultaneously reach 0.01 with four direct vacuolar supply caps of 0.01, their complete three-compartment balance rows force source = hydrolysis = 0.01 and intact dipeptide secretion = 0. Water import must be 0.04 and ATPase flux 0.035. The five amino acid transport fluxes must be Asp/Glu 0.01 each, Gly −0.04, Ala/Pro −0.01 each in their stored reaction orientations. These are conditional algebraic predictions, not results of an optimization or biological measurements.

## Identity and chemistry limits

The four hydrolyses consume Gly-Asp (`m1871[C_va]`, R2021), Gly-Glu (`m1862[C_va]`, R2029), Ala-Gly (`m1878[C_va]`, R2034), and Gly-Pro (`m1866[C_va]`, R2039). All twelve modeled dipeptide species across the three compartments lack chemical formulas. Vacuolar Gly also lacks a formula. Their SBML charge fields are zero; that does not establish exact chemical identity or protonation experimentally. The hydrolysis names say “cytosol” while the actual species are vacuolar; this conflict is retained.

Based on the displayed peptide names, nominal carbon counts are 6, 7, 5, and 7, with two nitrogens each. Four 0.01 supplies therefore correspond nominally to 0.25 mmol carbon atoms and 0.08 mmol nitrogen atoms per gDW per hour. These are name-derived accounting assumptions, not validation of missing formula fields or evidence of native protein-turnover supply.

## Actual execution and numerical audit

The completed primary run reports 114 calls, all `optimal`, with 2.0170670417 seconds inside actual solve calls and 220.11833525 seconds for the runner wall time. The independent audit checks all 264,081 net reaction flux values and 528,162 split-variable primals, their exact native variable bounds, every actual sparse constraint row, every stored metabolite balance, and each numeric objective. It also verifies the intended objective and direction for each growth, pFBA, hydrolysis, FVA, zero, and nucleotide case. All original stoichiometries and GPR strings remain fixed; only the declared scenario bounds, four diagnostic source columns, and three balanced nucleotide dissipation columns are allowed. Each pFBA growth floor equals that scenario's saved primary optimum without a deliberate slack.

The largest independently recomputed residual or objective discrepancy is 6.8326634932×10⁻¹¹, below the unchanged 10⁻⁷ criterion. Actual native constraints and the stored stoichiometric matrix are checked separately. In the 15 closed-model cases, 72 split coefficients in R1372 differ slightly from stored stoichiometry, with maximum difference 3.6030446000×10⁻¹¹. R1372's forward and reverse variables are exactly fixed to zero there, so these differences have no feasible-set effect. They are preserved in the audit output rather than called exact matrix identity. All other mass-row coefficients, including every culture case, match exactly. The cause is consistent with solver-copy numeric serialization, but was not independently traced into the library implementation.

The eight energy-repair reaction bounds and stoichiometries match raw E5 after both the E5 and candidate culture loaders. No artificial source is present in C0/C1/C2; every one-sided material column and all six biomass columns are closed, and forced maintenance is released. The 12 ATP/GTP/UTP/CTP maxima are zero with reported optimal status, including the C2 stress scenario where the four target connection upper bounds are 1000. Three additional exact-zero controls have all reaction fluxes zero. These results support the finite closed-energy regression at the declared tolerance, not universal thermodynamic correctness.

The independently checked tabular coverage comprises 10,633 complete ledger rows across 403 pool–witness groups; 1,217 related-flux rows; 332 dipeptide-fate rows; 124 nominal C/N rows; 76 growth/hydrolysis rows; 15 energy/zero-control rows; three joint controls; and 20 targeted FVA rows. Every pool ledger contains all incident reactions, including zero-flux terms. No source flux is inferred from its upper bound. The executable `auditor_results.py` and `auditor_results.json` preserve this check and each file's identity.

## Interpretation checked against the results

E5 alone gives growth 1.429177892279185 h⁻¹. Opening the four connections alone, or adding the four artificial pools while leaving the connections shut, gives the same optimum within numerical precision. Combining the openings and four finite pools gives 1.434139669972012 h⁻¹, a 0.347177% increase. This is an extra-material scenario. It is consistent with feasible-set expansion; it does not demonstrate native peptide supply or required vacuolar activity for ordinary growth. The additional exact relaxation check G3 ≥ G2 passes with difference 0.004961777692822 h⁻¹.

All three joint controls—no growth floor, 99% growth, and growth within 10⁻⁶ h⁻¹ of the primary optimum—simultaneously carry 0.01 in every hydrolysis. All four sources also equal 0.01. The complete balances reproduce the preregistered conditional prediction: water 0.04; pump 0.035; Asp/Glu export 0.01 each; Gly export 0.04; Ala/Pro export 0.01 each. Thus individual maxima were not substituted for a simultaneous feasible witness.

At 99% of the augmented maximum, all five targeted FVA minima are zero; maxima are 0.01 for each hydrolysis and 0.035 for R795. This does not contradict positive pFBA usage: baseline growth itself exceeds the relaxed 99% floor. At the preregistered absolute gap of 10⁻⁶ h⁻¹, the hydrolysis minima are respectively 0.00999015330, 0.00999250794, 0.00998998155, and 0.00999380154; the pump minimum is 0.03499015330. The added nutrients are needed to retain nearly all of their small computed growth benefit. This conditional near-optimum requirement is not a gene-essentiality result; no gene knockout was performed.

The final report's section 4 was reread after its wording correction. It now explicitly distinguishes R795's cytosolic water consumption, obtained from that reaction's stoichiometry multiplied by its saved flux, from the complete vacuolar-water and cytosolic-ATP ledgers. It correctly states that no complete cytosolic-water ledger was separately summarized. The earlier wording gap is resolved; this documentation clarification changes neither the raw results nor the numerical audit and requires no new optimization.

## Finite claim coverage

| Claim | Scope and evidence | Status |
|---|---|---|
| Pinned E5 identity and counts | Raw XML SHA, 2314 reactions, 1877 species | Supported |
| Complete water row and shared hydrolysis limitation | All 5 incident reactions | Supported |
| Asp/Glu closed-exit limitation | Complete 2-reaction rows for each | Supported |
| Five amino acid recycling equations | Complete 5 amino acid rows | Supported |
| Complete local proton account | All 22 H⁺ neighbors plus 2-reaction Ile guard | Supported |
| Joint 0.01 conditional water/pump values | Explicit algebra with four finite direct sources | Supported conditionally |
| Dipeptide labels, compartments, formula gaps | All 12 dipeptide species and target hydrolysis columns | Supported |
| Nominal C/N at joint 0.01 | Name-derived composition, explicitly conditional | Supported conditionally |
| Candidate XML mutation scope | Full normalized XML comparison after reversing exactly four allowed bounds | Supported |
| Solver results, full vector feasibility, objectives and growth bounds | All 114 calls, complete net/split primals and sparse matrices; every reported ledger/table row checked | Supported within stated numerical scope |
| Closed-energy regression after opening | All 12 nucleotide maxima, three exact-zero controls; closure and no-source scope independently checked | Supported within stated finite test scope |
| Native transport mechanism, exact substrate identity and physiological supply | No new experimental or primary-source evidence within this task | Unresolved |

Final finite claim set: 12 claims; 11 supported (including 2 explicitly conditional algebra/accounting claims), 0 contradicted, 1 unresolved, 0 unchecked. Audit coverage is 12/12 for this listed set; support is 11/12. Unknown native biological mechanism is a limitation on biological acceptance, not a reason to deny the completed computational test. This audit is not a whole-model chemistry audit, independent second-solver optimality proof, or native-function validation.
