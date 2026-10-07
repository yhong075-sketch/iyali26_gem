# External iYLI647 structural comparison and exact certificates

Audit date: 2026-09-11 UTC. This is a static audit of the fixed `iYLI647_corr_3.json`, its saved mapped28/PO1f screen conditions, and the existing iYali26 diagnostic snapshot. No optimization, model editing, new annotation, or network lookup was performed. Exact row combinations below were independently evaluated using the JSON's decimal coefficients; they are algebraic consequences of `S v = 0` and the recorded directions, not new FVA results.

| claims | audited | supported | partial | contradicted | unresolved | not_audited |
|---:|---:|---:|---:|---:|---:|---:|
| 6 | 6 | 6 | 0 | 0 | 0 | 0 |

The six claims are E1-E6 below, bounded to the specified structural statements. Native enzyme function, correct biological GPR assignment, thermodynamic validity and whether the authors deliberately “fixed” any screen result are not established and are outside this denominator. Better agreement on a positive label does not establish that the responsible model structure is biologically correct.

All named target loci below have unverified formal native gene names in this audit. Their described functions are model roles. In particular, the same mapped locus can have different reaction assignments in the two models.

## E1. GND is coupled to ER NADPH and sterol demand, not directly to cell-wall 6PG demand

Target: **YALI1B20462g / YALI0B15598g**, 6-phosphogluconate-dehydrogenase model role. External GPR is the single target at GND. The external model retains reversible TALA, TKT1, TKT2, RPI and RPE, so absence of the nonoxidative pentose-phosphate route is not the explanation.

The complete connections of the critical pools are:

| Pool | All participating reactions | Consequence |
|---|---|---|
| `6pgc[c]` | PGL produces; GND consumes | `v_PGL = v_GND` |
| `6pgl[c]` | G6PDH2 produces; PGL consumes; 6PGLter transports to ER | Cytosolic lactone balance |
| `6pgl[r]` | G6PDH2er produces; 6PGLter transports | ER lactone balance |
| `nadph[r]` | G6PDH2er produces; C24STRer and SQLEr consume | No cytosolic NADPH transport or other ER NADPH producer in this JSON |

Adding the first three metabolite rows and subtracting the NADPH row gives exactly:

`S_6pgc[c] + S_6pgl[c] + S_6pgl[r] - S_nadph[r]`

whose only nonzero reaction coefficients in the entire model are:

`G6PDH2: +1; GND: -1; C24STRer: +1; SQLEr: +1`.

Therefore:

\[
v_{GND}=v_{G6PDH2}+v_{C24STRer}+v_{SQLEr}.
\]

All four reactions have lower bound zero. GND closure necessarily forces ER sterol reduction and squalene epoxidation to zero, even if cytosolic NADPH is supplied elsewhere.

To check possible ergosterol transport/storage escape routes, the following exact weighted row sum was also evaluated over every reaction:

`S_ergst[c] + S_ergst[r] + S_ergst[e] + 100 S_SEergst_SC[c] + 40.82 S_SEparticle_SC[c]`.

Only these coefficients remain:

| Reaction | Exact coefficient |
|---|---:|
| C24STRer | +1 |
| EX_ergst(e) | -1 |
| biomass_C | -0.035250723 |
| biomass_N | -0.060650767 |
| biomass_glucose | -0.037699449071928945 |
| biomass_oil | -0.07541242503624108 |
| membrane_rSC4_bal | -0.012635 |
| membrane_WOtag | -0.00073 |

The esterification/hydrolysis cycle cancels exactly. `SEparticle_SC[c]` participates only in SEp_form, so its reversible assembly cannot supply net ergosterol. EX_ergst(e) has no uptake in either the native JSON or saved mapped28 condition. All four biomass reactions and the two membrane drains are irreversible. Hence:

\[
v_{GND}\ge v_{C24STRer}\ge0.035250723\,v_{biomass_C}.
\]

Thus GND KO mathematically forces `biomass_C = 0` in the recorded setting. This is an exact structural explanation, not a hypothesis requiring a further growth solve.

**Local contrast:** iYali26 has reversible R1264/R1265 transporting NADP/NADPH between cytosol and ER. Its saved target-KO witness uses R1265 = +0.13764137791162726 and R1264 = -0.13764137791162726; ER sterol reactions R188 and R738 remain active. This establishes the modeling distinction, not that unrestricted cofactor transport is biologically justified.

Key external source lines: G6PDH2 16368; G6PDH2er 16898; GND 16957; PGL 17020; 6PGLter 8228; C24STRer 19020; SQLEr 19857; ERGSTter 21005; EX_ergst(e) 13601; SEp_form 25485; biomass_C 10006.

## E2. MCITDm's actual assignment is in homocitrate/lysine synthesis

Target: **YALI1F03803g / YALI0F02497g**. The local target has a methylcitrate-pathway model role. The external target's only GPR reaction, MCITDm, is named “2 methylcitrate dehydratase mitochondrial,” but its actual substrate is **hcit[m]**, named `2_Hydroxybutane_1_2_4_tricarboxylate`. It belongs to the homocitrate/lysine chain in the current JSON, not the local R95 chemical step.

| External reaction | Actual carbon-skeleton connection | Bounds |
|---|---|---|
| HCITSm | acetyl-CoA + alpha-ketoglutarate -> homocitrate | 0…1000 |
| MCITDm | homocitrate -> but-1-ene-1,2,4-tricarboxylate + water | -1000…1000 |
| HACNHm | that unsaturated tricarboxylate + water -> homoisocitrate | -1000…1000 |
| HICITDm | homoisocitrate -> oxag, NAD reduced | -1000…1000 |
| OXAGm | oxag -> 2-oxoadipate + CO2 | -1000…1000 |
| 2OXOADPtim | mitochondrial -> cytosolic 2-oxoadipate | 0…1000 |
| AATA | 2-oxoadipate -> aminoadipate by transamination | -1000…1000 |
| AASAD1/AASAD2 | aminoadipate -> semialdehyde, ATP/NAD(P)H consumed | 0…1000 |
| SACCD1 | semialdehyde -> saccharopine | -1000…1000 |
| SACCD2 | saccharopine -> lysine | -1000…1000 |

Two apparent alternate entrances are structurally inactive: `oxag[c]` occurs only in OXO2Ctm, forcing that reaction to zero; `gtycoa[c]` occurs only in oxogludeh(e), forcing that reaction to zero in either direction.

An exact whole-model certificate is obtained by summing these **15 rows**, each with coefficient +1:

`b124tc[m], hicit[m], oxag[m], oxag[c], 2oxoadp[m], 2oxoadp[c], gtycoa[c], L2aadp[c], L2aadp6sa[c], saccrp_L[c], lys_L[c], lys_L[m], lys_L[e], lystrna[c], lystrna[m]`.

All columns cancel except MCITDm, the lysine exchange and the four biomass demands:

\[
v_{MCITDm}-v_{EX\_lys}
=0.275004895\mu_C+0.031527805\mu_N
+0.26860736050872014\mu_{glucose}
+0.15658056771954665\mu_{oil}.
\]

Every biomass flux is nonnegative. The saved mapped28 lysine uptake upper limit is 0.02466. MCITDm KO therefore implies:

\[
\mu_C\le {0.02466\over0.275004895}
=0.08967113112659321.
\]

Dividing by saved external WT 1.2171871273706656 gives **0.07367078496820643**, or **7.367078496820643%**. Both the absolute upper bound and ratio equal the previously saved KO result digit for digit. The exact certificate supplies an upper bound; the existing optimal run attains it. With native JSON lysine uptake closed, the same equality forces zero biomass after this KO.

The current JSON has no MCS reaction or propionyl-CoA/propanoyl-CoA pool supporting the proposed propionyl-CoA disposal explanation. THR/ILE synthesis is represented separately through THRD_Lm, ACHBSm, KARA2im, DHAD2m, 3MOPtm and ILETA; it is not evidence of that proposed MCITDm linkage.

**Local contrast:** R25 -> R95 -> R552 is the actual local propionyl-CoA/methylcitrate chain. The external MCITDm substrate topology instead resembles local **R24**, which also has a misleading “2-methylcitrate dehydratase” name but consumes homocitrate and has a different GPR. Local R24 includes additional proton stoichiometry, so no strict chemical-equivalence claim is made. Cross-version gene-ID correspondence does not mean the two KOs remove the same chemical reaction.

Key external source lines: HCITSm 18916; MCITDm 20716; HACNHm 20674; HICITDm 20687; OXAGm 20821; 2OXOADPtim 21739; AATA 20660; AASAD1 8268; AASAD2 20242; SACCD1 20862; SACCD2 20947; OXO2Ctm 21017; oxogludeh(e) 26912; EX_lys_L(e) 13931; biomass_C 10006.

## E3. The UDP-glucose bypass is prevented by direction and unavailable galactose supply

Target: **YALI1A02775g / YALI0A02310g**, UDP-glucose-synthesis model role. External GALT exists; it must not be described as an absent counterpart of local R319. The decisive difference is that external **UGLT is [0,1000]**, while the observed local bypass requires **reverse R780**. Thus the local exact combination `R319 - R780` cannot be used with the external directions.

External UGLT consumes galactose-1-phosphate + UDP-glucose and produces UDP-galactose + glucose-1-phosphate. GALT consumes galactose-1-phosphate + UTP and produces UDP-galactose + PPi. External UDPG4E is reversible, but its UDP-galactose source cannot be replenished here:

- Galactose and melibiose exchange uptake are both zero in the native JSON and mapped28 medium. Their inward proton-symport reactions are irreversible, so their exterior balances force inward fluxes to zero.
- The other galactose-producing hydrolases EPMGH, GALIGH and GGLGH have substrates `epm[c]`, `1Dgali[c]` and `ggl[c]`, respectively, each appearing only in that one reaction. Those fluxes therefore are zero. Melibiose hydrolase GALS3 has no incoming melibiose.
- Galactose balance then forces GALK=0. Galactose-1-phosphate balance is `GALK-GALT-UGLT=0`; both consuming reactions are nonnegative, hence GALT=UGLT=0. UDP-galactose balance then forces UDPG4E=0.
- With GALU also closed, the UDP-glucose balance forces its irreversible consumers, including 13GS, to zero. `13BDglcn[c]` has 13GS as its only source, and biomass_C requires 0.943396927 units per unit growth. Therefore this target KO cannot maintain growth under these fixed supplies.

This proves the particular modeled restriction and its growth implication. It does not establish which uridylyltransferase direction or GPR is correct in vivo.

Key lines: UGLT 26898; GALT 16484; UDPG4E 24606; GALK 16455; GALt2 16514; EX_gal(e) 13667; EX_melib(e) 13975; MELIBt2 19433; GALS3 16470; EPMGH 13105; GALIGH 16441; GGLGH 16772; 13GS 7996.

## E4. The corresponding mitochondrial ACP synthesis exists but cannot carry steady-state flux

Target context: **YALI1C07638g / YALI0C05951g**, fatty-acid-desaturation model role. External **FAS161ACPm exists**; it must not be called absent. Its reaction consumes malACP[m] + myristoyl-ACP[m] + 3 NADPH[m] + 4 H+[m] + O2[m], producing C16:1-ACP[m] + ACP[m] + CO2[m] + 3 NADP[m] + 3 H2O[m].

However, `malcoa[m]` appears only in reversible MCOATAm, with coefficient -1. Its balance gives `v_MCOATAm=0`. MCOATAm is the only producer of `malACP[m]`; every other reaction involving that pool consumes it irreversibly.

Adding the two rows `S_malcoa[m] + S_malACP[m]` gives exactly:

\[
-\sum_{r\in M_9}v_r-3v_{FAS80ACPm\_L}=0,
\]

where `M_9 = {FAS100ACPm, FAS120ACPm, FAS140ACPm, FAS141ACPm, FAS160ACPm, FAS161ACPm, FAS180ACPm, FAS181ACPm, FAS182ACPm}`. All ten reactions have lower bound zero, so **all ten are forced to zero in every steady-state feasible solution**. This is an algebraic blockage certificate, not an FVA run.

Consequently the observed local R2233 C16:1-ACP synthesis route cannot be reproduced through its external counterpart. The proof is limited to these ten synthesis reactions; it does not assert that every ACP transport/hydrolysis reaction or every unsaturated-fatty-acid source is blocked. This gap is not evidence that the external network is more biologically accurate.

Key lines: MCOATAm 19340; FAS161ACPm 15650. Other certified synthesis IDs/lines: FAS100ACPm 15444; FAS120ACPm 15494; FAS140ACPm 15546; FAS141ACPm 15577; FAS160ACPm 15616; FAS180ACPm 15689; FAS181ACPm 15720; FAS182ACPm 15740; FAS80ACPm_L 15798.

## E5. Cytosolic malonyl-CoA has one feasible net source

Target: **YALI1C15991g / YALI0C11407g**, acetyl-CoA-carboxylation-related model role. External `malACP[c]` appears only in MCOATA, so `v_MCOATA=0`, despite that reaction being reversible. It cannot provide malonyl-CoA in reverse.

The complete `malcoa[c]` adjacency has 16 reactions: ACCOACr, MCOATA and fourteen irreversible FAS consumers. Summing `S_malcoa[c] + S_malACP[c]` cancels MCOATA and leaves:

\[
v_{ACCOACr}-\sum_{r\in C_{11}}v_r
-3(v_{FAS80COA\_L}+v_{FAS240\_L}+v_{FAS80\_L})=0,
\]

where `C_11 = {FAS100COA, FAS120COA, FAS140COA, FAS160COA, FAS180COA, FAS100, FAS120, FAS140, FAS160, FAS180, FAS260}`. All fourteen consumer fluxes are nonnegative. Thus **ACCOACr is the only feasible net source of cytosolic malonyl-CoA**, and closing it forces all fourteen listed FAS consumers to zero. ACCOACr's reversible annotation does not permit a negative steady-state net flux for this pool.

Enumerating every reaction involving the pool excludes local R2004's reverse `oxaloacetate + acetyl-CoA -> malonyl-CoA + pyruvate` and R2121's `3-oxopropanoate + NADP + CoA -> malonyl-CoA + NADPH` as external entry routes. This is a metabolite/direction check, not merely a search for matching reaction IDs.

The target KO also removes other GPR-associated FAS reactions. This certificate alone does not uniquely attribute the entire growth effect to ACCOACr rather than the combined intervention, nor validate the native enzyme assignment.

Key lines: ACCOACr 10768; MCOATA 19326. Complete consumer IDs/lines: FAS100COA 10903; FAS120COA 10941; FAS140COA 11164; FAS160COA 11212; FAS180COA 11231; FAS80COA_L 11306; FAS100 15425; FAS120 15475; FAS140 15527; FAS160 15597; FAS180 15670; FAS240_L 15760; FAS260 15779; FAS80_L 15831.

## E6. There is no formate/ATP alternative entry into the GAR formylation product

Target: **YALI1D03865g / YALI0D03069g**, GAR-formyl-transfer model role. The external pool `fgam[c]` appears in exactly two reactions:

- **GARFTi**, [0,1000], produces it from GAR + 10-formyl-THF, also producing THF and H+.
- **PRFGS**, [0,1000], consumes it with ATP + glutamine + water to produce fpram + ADP + glutamate + phosphate + H+.

Its exact balance is `v_GARFTi-v_PRFGS=0`. Closing GARFTi forces the downstream PRFGS flux to zero. There is no reaction producing this same intermediate through the local R1892 formate/ATP route, and reverse PRFGS is not permitted. This supports removal of that specific bypass. It does not by itself characterize the remaining purine salvage or nutrient-supported residual growth; the previously recorded external KO has nonzero growth, not complete modeled death.

Key lines: GARFTi 16556; PRFGS 22149.

## Source identities and scope limits

All reaction line numbers above refer to the exact 29,861-line external input below, using one-based reaction-ID lines. They are not line numbers in a newly formatted JSON.

| Source | SHA-256 |
|---|---|
| `artifacts/iyli647_screen_20260910/inputs/iYLI647_corr_3.json` | `329be540c099409c2c7b76ee581a23f86eefdbbffe97e9b03b39b5e9c014b5d2` |
| `artifacts/iyli647_screen_20260910/mapped28_po1f/run_manifest.json` | `9f2950c7c15186830eb6c4fcfb18dc18f7102a211394b6d318153a647199c863` |
| `artifacts/iyli647_screen_20260910/common_positive_comparison.tsv` | `1616198c1f33644ef9b7f41ef795b37260a380f6401d59667a1d4fb8e52f71fc` |
| local `nonessential_diagnosis_20260911/reaction_snapshot.json` | `d9eb0a95af7cd6944b37f0667623c565f2a0999dc82b303a2d4d54774ca4bde3` |
| original `model_metadata_trna.xml` represented by local snapshot | `d274bad3050e3c9220a8b6287eae847f3bf1334892284d565a6c4d96b38135a0` |

The mapped condition is the previously saved partial SD-Leu/PO1f mapping, not a claim of full biological equivalence between the two models. GND and MCITDm have exact algebraic explanations above; no additional solve is needed to establish those implications. Any future rescue test would be a new, bounded intervention requiring preservation of the current conditions and evidence, and would still only test model behavior. No reaction or GPR should be changed solely to improve label agreement.

Only this new audit document was written for this static comparison. Existing diagnostic outputs and model files were preserved.
