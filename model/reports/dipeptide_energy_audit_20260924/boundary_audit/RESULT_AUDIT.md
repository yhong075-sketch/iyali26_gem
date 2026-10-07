# Independent closed-energy result audit

Auditor: `/root/energy_boundary_audit`, 2026-09-24. Read-only scientific audit of the saved local results; no optimization was performed by the auditor. The accepted results checked here culminate in `energy_delivery` (13 unique closed-model LP labels) plus two native-medium growth LPs, for **15 unique accepted solves**. Earlier `energy_complete` (6 labels) and `energy_final` (12 labels) are preserved intermediate stages, included in the final13 closed labels. The initial failed `energy` output is excluded from accepted scientific conclusions.

## Input, closure and effective problem

All three execution input SHA256 values match the manifest and the earlier fixed-input check. The 2,314 reaction IDs, stoichiometries and GPRs, and 1,877 metabolite records match the inspected historical effective model. Independently reconstructed closure is exactly the 183 single-sided boundary reactions plus 6 biomass producers; the only other change is xMAINTENANCE lower bound 7.8625 to 0. No internal direction was broadened. The saved 1,877 material rows retain zero lower/upper bounds, with no custom constraints. The zero vector is therefore an exact algebraic feasible point of the closed problem; no extra LP was needed to establish this.

There is a serialization detail that prevents claiming that every native matrix entry is byte-identical to the pre-copy matrix. COBRA/Gurobi copying reordered solver columns; comparing by variable name resolves that ordering change. It also rounded 72 coefficient entries belonging exclusively to the forward/reverse variables of **R1372**, with largest coefficient difference `3.603044599973643e-11`. R1372 is fixed to zero by the explicit closure, so those differences do not change the present closed feasible set. Every other native coefficient matched the reaction stoichiometry exactly. No steady-state row was deleted.

The original maintenance column is balanced under the stored neutral-species convention and is the sole dissipation objective. All direction-intervention saved bounds show xMAINTENANCE `[0,1000]`; `[1,1]` is present only in fixed-unit witness extraction. The source code starts each intervention from `closed.copy()` and restores the maximum-maintenance objective by inheritance from that baseline, so fixed-unit/L1 objectives did not leak into intervention tests.

## Independent result reconstruction

For all **13 optimal closed LP records / 30,082 full flux values**, the auditor independently rebuilt the complete reaction bounds from the closed baseline and each declared directional block, recomputed all 1,877 `Sv` rows, and recalculated the objective from saved net fluxes. Actual solver-parameter records match the declared 1 thread, 60-second limit, Presolve 0 and `1e-7` tolerances for every accepted solve.

| Quantity | Maximum independently recomputed discrepancy |
|---|---:|
| Absolute material-row residual | 2.2737367544323206e-13 |
| Reaction bound violation | 0 |
| Objective versus net-flux calculation | 0 |
| Exported ledger versus individual weighted contributions | 0 |

All **33 weighted metabolite-ledger rows**, including water, protons and phosphates when present, were regenerated from the saved stoichiometry and complete vectors. All **18 displayed reaction rows** match the threshold, original signed flux, current bounds and coefficients. Each net ledger equals the opposite of maintenance with total residual within the declared tolerance. The audit also matched all 8 intervention-table rows to their complete saved LP record. This checks feasible witnesses, exports and reported solver outcomes; it is not a second optimization proving every L1 optimum independently.

### Witnesses and causal scope

1. Initial unit-dissipation witness is exactly **R694 +1, R_PGAM1_PhosHydro +1, xMAINTENANCE +1**, with internal L1=2. Summing the first two stored columns gives `ADP + phosphate -> ATP + H2O`, exactly opposite the maintenance column. No exchange, biomass or material source participates.
2. Blocking either R694 forward or R_PGAM1_PhosHydro forward individually in the full closed model still permits maintenance **1000**, the inherited upper bound. Their respective L1 extractions are the **same second path**, up to floating-point noise: R594 −2, R2010 +1, R_NTP3pp +2, r0242 −1, maintenance +1; L1 approximately6. These must not be counted as two distinct mechanisms.
3. Blocking both R_PGAM1_PhosHydro forward and R_NTP3pp forward still permits maintenance1000. Its distinct extraction is R603 −2, R2010 +1, R_NTP7 −2, r0242 −1, maintenance +1, L1=6.

The additional single-direction interventions against R_NTP3pp forward, r0242 reverse, R594 reverse and R2010 forward also each retain maximum maintenance1000. The joint block of all three identified synthetic directions (R_PGAM1_PhosHydro forward, R_NTP3pp forward and R_NTP7 reverse) also retains maximum maintenance1000. All 8 tested cases are positive alternative-path cases; none is a zero-ATP or infeasible result. A direction absent from a blocked witness does not establish that it is dispensable in biology, and these tests do not establish that a particular correction alone repairs the model.

The L1 objective minimizes summed absolute internal flux under D=1. It is not a cardinality optimization and does not prove mathematically smallest reaction support. Maintenance1000 is a finite imposed upper bound; it does not assert an infinite mathematical optimum or a physiological ATP rate.

## Failure preservation, reuse and budget accounting

The first attempt used a copied solver whose native parameters reverted to defaults, including infinite TimeLimit, and failed JSON serialization after one LP. Its source, log, `attempt1_failure.json` and partial output are retained. That LP is not accepted as a planned diagnostic, and its actual solver wall time is **unknown**. The configured60 seconds is a conservative accounting reservation, not a measurement or a demonstrated upper bound for this failed attempt.

The second attempt correctly configured the solver and produced two valid LPs before the independent-pool object lookup caused a ledger assertion. Comparing the archived second script to the initial-complete script verifies that the scientific fix was to look up maintenance metabolites within the same copied model. The added reuse path checks the fixed problem identity, purpose, solver parameters, bounds, objective and material residuals. The first two full LP records are exactly equal between `energy_validated` and `energy_complete`; the first six records are exactly equal between `energy_complete` and `energy_final`; all12 records from `energy_final` are exactly equal in `energy_delivery`. The four witness exports and33 ledger rows are byte-identical between these final two stages. Reuse did not recompute those LPs.

The final audited set contains **13 unique accepted closed solves +2 native-medium growth solves**, plus **1 excluded first-attempt solve**, totaling **16 actual LP calls** across these attempts. Accepted closed calls have measured cumulative solver wall time **0.20750641380436718 seconds**, and the two growth calls total **0.07931320706848055 seconds**, giving **0.28681962087284774 seconds measured for the15 valid calls**. With the60-second first-attempt accounting reservation, the budget ledger is60.28681962087285 seconds, below the declared1200; it must be described as reserved accounting, not all-attempt measured runtime. No tolerances were relaxed. Earlier static test logs report two tests passed. The later plain-Python invocation failed to import memote and its failed log remains preserved. The corrected existing-workspace-venv invocation reports **8 tests passed in23.415 seconds** in `tests_final_venv.log`. The intended selection covers two closed-energy static tests, five non-solving dipeptide tests and one R1159 roundtrip test; the toy pFBA optimization test is excluded. These software checks are not additional biological LP validations. No scientific tolerance was adjusted to resolve the import issue.

The auditor independently recomputed all **20 protected-file hashes** from `before.json`; all remain unchanged. Preservation is limited to that declared set, not a claim that all historical dirty workspace state was restored.

## Restored-medium growth comparison

The two growth runs load the original effective model, not the closed diagnostic model. Their saved **4,628 full flux values** and all2,314 bounds per case were independently compared with the historical PO1f/SD-Leu effective baseline. Original growth bounds match every baseline reaction. The direction-only comparison differs only in the three explicitly blocked directions. Every original boundary bound, biomass_C `[0,1000]`, and maintenance `[7.8625,1000]` is restored; maintenance flux is7.8625 in both cases. Declared solver parameters match both row records, and the driver inherits the original biomass_C maximization objective from the unchanged loader/model.

| Growth case | Growth h⁻¹ | Independently recomputed maximum Sv residual |
|---|---:|---:|
| Original effective model | 1.8718823069402888 | 1.9895196601282805e-13 |
| Three-direction-only diagnostic | 1.8718823069402994 | 1.7053025658242404e-13 |

Both have zero independently recomputed bound violation, and their biomass/maintenance recorded values exactly equal the saved fluxes. The tiny growth difference is numerical, not a biological improvement. Both have R1372 flux0, so the already documented fixed-zero/copy coefficient issue does not invalidate these saved feasible points. These are two WT growth comparisons, **not gene knockouts**, a complete energy repair, or evidence that the four target dipeptide reactions acquire native supply. The final closed triple-block test still generates ATP.

Across closed and growth cases, this audit checks **34,710 flux values**. The growth script and config SHA match their report, and its unchanged-loader provenance is tied to the same protected model/media/strain files. The20 declared protected-file hashes were recomputed again after delivery and all remain unchanged.

## Evidence identities and coverage

- `energy_complete/closed_energy_manifest.json`: `66d7921f4ce994d052e0e4fcb9cc0b9d0cab64a61b712f98dbee24714111f752`.
- `energy_final/closed_energy_manifest.json`: `64a2a19191d44be2d262fc60b0cb401e1ad1a4e9f5cd1c018d8c622f78c9c460`.
- `energy_delivery/closed_energy_manifest.json`: `8aed55cbabf068aa2d28e2911a77259eaf59f0c0ce09df6e39c7446432267eaa`.
- `growth_comparison.json`: `75685f83ee7bd22a5571570cf514b04921fe3c53e6e61ee73231035f67a8fdc3`.
- Growth script: `5cf823c247ed53873a17d20394c7a2e1d405b5f4fd8695fc310ac97f654904d5`.
- Final closed-driver exact snapshot `boundary_audit/diagnose_closed_energy_delivery.py`: `efcd0e37b56f660005f73dc3407a935d3822406a01613f017db40ebde492e63f`.
- Six-LP script snapshot `diagnose_closed_energy_initial_complete.py`: `e72f3d945034837fc53997490c13ba45bfb5528cad9238b0cd721e3cc4d736aa`.
- Twelve-LP script snapshot `diagnose_closed_energy_extended.py`: `cef3b571d7c80beafcf15fccab470dca401a7c891d8833fb7d942b3d3aeebaba`.

Diagnostic source changed between intermediate stages for follow-up support; the audit binds the above runs to their exact archived script SHA rather than silently assigning current source to historical results. All other listed loader source hashes and execution-input hashes match. An initial auditor bytewise matrix equality assertion failed due to the above variable ordering and fixed-zero R1372 rounding; the check was corrected to compare by variable name and identify every numerical difference, without changing scientific inputs. An auditor inline-code syntax typo was corrected before the successful full-vector verification and caused no solver call or file mutation.

Coverage: 9 finite computational claim groups checked (identity/closure, constraints/matrix, parameter records, full-vector feasibility/objectives, net ledgers, direction tests, exact reuse/failure accounting, restored-medium growth comparisons, protected files); all supported within the stated limits. Biological correctness of the responsible internal reactions, completeness of all possible ATP cycles, strict minimum support, and the first failed attempt's actual runtime remain outside this verification. The energy acceptance criterion is **failed**: explicit closed ATP-generation witnesses exist.
