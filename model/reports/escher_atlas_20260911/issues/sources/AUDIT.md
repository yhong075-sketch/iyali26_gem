# Independent source audit: eleven nonessential predictions

Audit date: 2026-09-11 UTC. Scope: `results.json`, `reaction_snapshot.json`, `diagnose.py`, `controls.json`, `controls.py`, the identified original XML, saved comparison and previous screen manifest. This audit read the saved results and independently recalculated their arithmetic, Boolean GPR effects and stoichiometric residuals. It did not run a solver, edit a model, change a condition, or perform new functional annotation.

## Coverage and conclusion

The claim denominator is explicitly the eleven target-level mechanism observations plus five paired reaction-closure observations. Identity and numerical checks are reported separately; biological extrapolations are outside this denominator.

| total_claims | audited | supported | partial | contradicted | unresolved | not_audited |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 16 | 16 | 0 | 0 | 0 | 0 |

All sixteen **bounded model observations stated below** are supported. This does not establish native enzyme functions, thermodynamic validity, unique flux routes, or equivalence between the external model and iYali26 conditions. Three stronger interpretations that must not be made are explicitly rejected below.

## Input identity and numerical verification

- The original XML has exactly 2,315 reactions. Every reaction ID and every stoichiometric coefficient in the snapshot matches the XML. The only GPR difference is the declared runtime plasmid substitution at R45. The 21 bound differences are the declared R612 closure and 20 exchange-uptake settings; all 35 active-medium uptake values match the snapshot bounds.
- The model, medium, strain profile, simulation-context fingerprint and overlay record match the preceding iYali26 screen. Seven protected input hashes match before, after and at this audit. All three recorded imported source files and both diagnostic scripts match their recorded hashes. The control record points to the exact initial `results.json` hash below. Historical dirty state is preserved in the original result; these checks do not reconstruct or certify a clean historical workspace.
- The eleven targets independently selected by `scenario == mapped28_po1f`, raw external ratio `< 0.10` and raw iYali26 ratio `>= 0.10` match all eleven saved target records, including source rows and original ratios.
- Twelve initial and ten control runs are present: 22 recorded primary LP solves, under the declared combined limit of 32. All report `optimal`, finite nonnegative growth and complete finite 2,315-reaction flux vectors. The audit checked 50,930 flux values without solving.
- The XML objective is maximization of `biomass_C` with coefficient 1. All 22 reported objective values equal the saved `biomass_C` flux. All ratios independently equal raw growth divided by the initial WT growth, 1.8718823069403.
- Expected bound changes were independently evaluated from the snapshot's Boolean GPR plus each explicit reaction closure. They match every run's recorded changes. Bound feasibility was recomputed against these independently derived KO bounds, not merely accepted from `changed_bounds`.
- Recalculated maximum absolute steady-state residual over all 22 witnesses is **1.3542143209761826e-13**; maximum bound violation is **0**. Initial-run maximum is the same; control-run maximum is 7.945215123837668e-14. These are checks of the model's `S v = 0`, not elemental/charge balance or thermodynamic checks.
- Optimality is the saved Gurobi status; this audit verified primal feasibility and objective arithmetic, and did not reconstruct a dual certificate. Recorded solver settings match between stages: Gurobi 13.0.1, Threads 1, TimeLimit 30 s, FeasibilityTol/OptimalityTol 1e-7, Method 0, Presolve 0; COBRApy 0.30.0 and optlang 1.8.3 are recorded by the initial run.

## Eleven mechanism observations

All systematic IDs below are identifiers, not verified native gene names. Formal names were not established by this audit. Functions are model reaction roles and previously saved homology annotations, not new experimental evidence. Flux values are raw model flux units; positive/negative signs refer to the stored reaction orientation.

| Claim | Target: iYali26 / external ID; model function | Independently checked evidence | Verdict |
|---|---|---|---|
| M01 | YALI1C08702g / YALI0C06490g; GDP-mannose synthesis | R540 remains true through its OR branch. No reaction bounds change. | supported |
| M02 | YALI1C32184g / YALI0C23364g; protein mannosylation | R289 remains true through its OR branches. No reaction bounds change. | supported |
| M03 | YALI1E22736g / YALI0E18964g; lysophosphatidic-acid acyl transfer | R1846 remains true through its OR branches; R1843 is already closed before KO. No additional reaction bounds change. | supported |
| M04 | YALI1E25018g / YALI0E21021g; beta-1,3-glucan synthesis | R4 remains true through its OR branch. No reaction bounds change. | supported |
| M05 | YALI1F00821g / YALI0F00506g; glutamine synthesis | R334 closes; R2081 remains available. R334 and R2081 have identical complete stoichiometric columns and compartments. KO R2081 flux is 3.3053143281566704, with R334 at zero. | supported |
| M06 | YALI1A02775g / YALI0A02310g; UDP-glucose synthesis | R793 closes. R319 runs +1.7787692653461697 and R780 runs -1.7787692653461697. Exact coefficient addition gives `S_R319 - S_R780 = S_R793`, including protons and compartments. | supported |
| M07 | YALI1B20462g / YALI0B15598g; 6-phosphogluconate dehydrogenase | WT and KO R639 are zero. KO nonoxidative PPP uses R764=-4.0893163105830155, R765=-1.8724195255840996, R766=-2.0521277649455762; R714=-3.924547290529676 and R712=+3.924547290529676 connect xylulose/ribulose/ribose phosphates. R_ALCDH_nadp_hi and R488 supply most, but not all, cytosolic NADPH. | supported |
| M08 | YALI1C07638g / YALI0C05951g; fatty-acid desaturation | R1434/R1435/R1869 close. The KO witness uses R2233 -> reverse R2234 -> R2235 -> R1572 -> R1448, each at magnitude 0.0005652289492295624. This connects mitochondrial C16:1-ACP production, export, hydrolysis, ER transport and CoA ligation, supplying ER palmitoleoyl-CoA without those target-dependent ER reactions. | supported |
| M09 | YALI1C15991g / YALI0C11407g; acetyl-CoA carboxylation-related model role | R87/R88/R1393 close. In the initial KO witness, GPR-free reversible R2004 runs -0.003956602644606937; its reverse consumes acetyl-CoA + oxaloacetate and supplies malonyl-CoA + pyruvate. It is the only positive cytosolic malonyl-CoA source in that particular witness, not the only available route in the network. | supported |
| M10 | YALI1D03865g / YALI0D03069g; GAR formyl transfer | R355/R1907 close. R1892 runs +0.23015851751788366, consuming formate, ATP and GAR to produce FGAR, ADP, phosphate and H+. This alternative has a different formyl donor/energy coupling from the closed folate-dependent routes. | supported |
| M11 | YALI1F03803g / YALI0F02497g; methylcitrate-pathway model role | R95 closes. Its actual substrates/products place it in the mitochondrial methylcitrate branch, despite its generic `aconitase` name. R25, R95 and R552 are all zero in both saved WT and KO maximum-growth witnesses. KO growth equals WT within numerical precision. | supported |

For M01-M04, standard GPR-only gene deletion does not change the LP constraints or objective, directly explaining the unchanged optimum in this model. For all other rows, the saved feasible flux vector is a witness to the route, not a statement that every optimum uses it.

M07 NADPH detail: the target KO's positive cytosolic NADPH contributions are R_ALCDH_nadp_hi 24.697894479199974; R488 6.479876507135084; R178 0.07165565470967471; R179 0.07165565470967471; R671 0.0524150453951095. The first two account for approximately 99.376% of gross positive supply in that saved vector. These data do not establish the chemical or native enzymatic validity of those routes.

Several routes are already used in the saved WT optimum. WT R639, R334, R87/R88/R1393, R355/R1907 and R25/R95/R552 are already zero. Avoid calling these observations KO-induced biological compensation or upregulation. M06 and M08 do show route changes between the two saved optimal vectors, which still are model flux observations rather than measured regulation.

## Five paired reaction-closure controls

Each pair compares closure of the candidate reaction alone with closure of that reaction plus the target gene. The target-gene-only reference is the initial result above. Every reaction-only control preserves approximately 100% WT growth.

| Claim | Target; function | Reaction closed | Reaction-only KO/WT | Target KO + reaction closure KO/WT | Interpretation and verdict |
|---|---|---|---:|---:|---|
| C01 | YALI1A02775g; UDP-glucose synthesis | R319 | 0.9999999999999998 | 0 | Supported model dependence on R319 when R793 is removed; does not validate the alternative enzyme assignment. |
| C02 | YALI1C07638g; fatty-acid desaturation | R2233 | 1.0000000000000013 | 0 | Supported dependence on R2233 in the target KO under this model/medium; the entire bypass has not been biologically validated. |
| C03 | YALI1F00821g; glutamine synthesis | R2081 | 0.9999999999999994 | 0 | Supported duplicate-reaction escape from gene KO; removing both reaction routes eliminates modeled growth. |
| C04 | YALI1D03865g; GAR formyl transfer | R1892 | 0.9999999999999992 | 0.02080421545150898 | Supported major rescue by the formate/ATP route. Remaining growth is 2.080421545150898%, not zero; no claim here identifies its remaining nutrient source. |
| C05 | YALI1C15991g; acetyl-CoA carboxylation-related role | R2004 | 0.999999999999999 | 0.9999365151542908 | Supported presence of another bypass: removing R2004 does not restore gene essentiality. |

C05 alternative path: R2004 is zero, while **R2120 and R2121 each run +0.0039563514602985175**. R2120 converts beta-alanine + 2-oxoglutarate to 3-oxopropanoate + glutamate + H+; R2121 converts 3-oxopropanoate + CoA + NADP+ to malonyl-CoA + NADPH. R2121 is the only positive cytosolic malonyl-CoA supplier in this saved double-intervention vector. Its model-assigned gene is YALI1C02550g (formal native name unverified; CoA-malonylating oxidoreductase model role); R2120's is YALI1E21915g (formal native name unverified; beta-alanine aminotransferase model role). This is not new evidence that those proteins catalyze these reactions in vivo.

## Corrections and limitations

1. **Reject “R2004 is the unique bypass / closing it restores essentiality.”** C05 directly contradicts that stronger claim; growth remains 99.99365151542908% and R2120 -> R2121 supplies malonyl-CoA.
2. **Correct the previous ordinary-aconitase interpretation of R95.** The stored equation is mitochondrial `2-methylcitrate + 3 H+ <=> (2S,3R)-3-hydroxybutane-1,2,3-tricarboxylic acid`. It is not the ordinary citrate/isocitrate reaction. Zero R25/R95/R552 in two witnesses plus unchanged optimum establishes that maximum growth can be achieved without this branch in the tested setting; it does not prove global FVA blockage or dispensability in other carbon sources.
3. **Reject “only R_ALCDH_nadp_hi and R488 supply NADPH.”** They are the dominant two sources in M07, but three additional positive sources are present.
4. The zero-growth paired controls are recorded optimal zero solutions, not solver failure interpreted as death. They establish model sensitivity to the stated intervention, not experimentally demonstrated lethality.
5. Unconstrained FBA optimal vectors can contain alternative optima or cycles. No pFBA, FVA, thermodynamic filtering, new isotope fitting or biological intervention was performed by this audit. The results should not justify deleting or changing model reactions solely to improve benchmark agreement.
6. The external comparison remains the previously recorded partial SD-Leu mapping/PO1f setting; this local diagnostic audit does not upgrade cross-model condition equivalence or turn the development positives into independent validation.

## Exact audited evidence identities

| File | SHA-256 |
|---|---|
| results.json | 3e38ddf833e0f84c7bdeb21cccd3166c58f81b37ae09adab0f65c6291f163abb |
| reaction_snapshot.json | d9eb0a95af7cd6944b37f0667623c565f2a0999dc82b303a2d4d54774ca4bde3 |
| diagnose.py | a920c18ac7d5282dc106eca84b0508a5bfe89e20fc598109886af21c39ece974 |
| controls.json | d97daac8d41a2dbed3ba11c4b49d6b8d3632bc717ac3a837796a7393b99873df |
| controls.py | 9c110be3ca3f52d57dd648b4c27fe5869bb68ac4813f68b43312edd395262b1a |
| model_metadata_trna.xml | d274bad3050e3c9220a8b6287eae847f3bf1334892284d565a6c4d96b38135a0 |
| saved sd_leu.csv | ed176d26a373f98cc413ed2e32a71f5f060a06e343f90f7db25cd32eff268e85 |
| saved po1f_sd_leu.json | 35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383 |

Paths for the last three files and hashes for the three imported source modules and seven protected files are retained in `results.json`. Evidence locators for each mechanism are `results.json -> runs[target] -> fluxes`, `reaction_snapshot.json -> reaction ID`, and `controls.json -> runs[target__reaction__double or __reaction_only]`. Only this new audit document was written during the audit.
