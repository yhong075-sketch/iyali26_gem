# Gap metabolites and blocked reactions of candidate a1237821, by the Ponce-de-León 2013 definitions (2026-10-26)

Scope and authorization: [TASK.md](TASK.md). Computation and reporting only; no model change; no commit.
Figures: [fig_overview.png](fig_overview.png) (causes, roots, counts), [fig_UM1_network.png](fig_UM1_network.png)
(the largest module), [fig_UM2-10_networks.png](fig_UM2-10_networks.png) (modules 2–10); Escher maps of the 20
largest modules: [escher/gap_modules_escher.html](escher/gap_modules_escher.html) (see section 4).
Small tables in this folder (CSV); the full workbook `gap_analysis_a1237821.xlsx`, the full tables, the
scripts and the raw outputs are in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/gap_metabolites_ponce2013_20261026/`.

## 1. Run

- **Model** `candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_20261023.xml`
  (SHA-256 `a12378210c7c66720d99fb7bf71d14cb7538ea6a3cf3e08919c85fb3c9c47d82`), the newest candidate on `main`
  (commit `11450bc`, PR #19). Not a designated reference. 2,316 reactions, 1,882 metabolites, 1,081 genes,
  182 exchanges, 3 demands, 69 reactions closed at (0,0).
- **Definitions** (Ponce-de-León, Montero & Peretó 2013, BMC Syst Biol 7:114): a reaction is *blocked* when its
  FVA minimum and maximum are both 0 (eq. 4–5); a metabolite is a *gap* when every reaction it takes part in is
  blocked (eq. 6–7); *RNP* / *RNC* are root non-produced / non-consumed metabolites found by scanning the rows
  of the stoichiometric matrix with the bounds; *DNP* / *UNC* are the gaps created downstream / upstream when
  the roots' reactions are removed; an *unconnected module* (UM) is a connected component of the bipartite
  graph of gap metabolites and blocked reactions; a *coenzyme pseudo-gap* (paper Fig. 2) is a cofactor that is
  synthesised but never drained, so its biosynthesis is blocked although the cofactor itself is not a gap.
- **Additions to the paper's scheme** (the reference classes of the `gem-gap-rnp-rnc-analysis` skill): the
  paper's plain row scan counts a metabolite that is in no open reaction as both RNP and RNC; here those are
  reported as *isolated*, and a metabolite that occurs in exactly one open reversible reaction (no steady-state
  flow possible, it acts as RNP + RNC) as *single reversible reaction*. Without that rule, propagation explains
  far fewer blocked reactions. Compounds that are never produced (or never consumed) in any compartment, whose
  copies are linked only by reversible transports, are reported as *pooled* roots.
- **Conditions.** Primary: every exchange open at (−1000, 1000), the paper's "any medium condition".
  Secondary: the Task A medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) + PO1f profile (uracil
  exchange 1000, R612 closed), the benchmark default.
- **Solvers.** cobra 0.30.0 with GLPK for the skill scripts (FVA, propagation, modules, fix tests). Gurobi
  13.0.2 (FeasibilityTol / OptimalityTol 1e-9, presolve default, 1 thread) as an independent cross-check of the
  blocked sets only.
- **LPs.** Three FVA scans (2 × 2,316 LPs each), the Gurobi cross-check (2 × 2 × 2,316), 621 single-fix tests,
  20 cofactor-drain tests, 66 combined / targeted module tests and 7 extra targeted tests (each re-testing only
  the reactions named in the test), five closure checks (growth, FVA 99.9 %, energy-generating-cycle LP; the
  energy-cycle LP was redone on 2026-10-27, see the correction in section 6).
- **Pipeline check.** On `E5_coq9_resp_pkg_r132_20261012.xml` (`a56f3649…`) the scan reproduces the skill's
  recorded reference exactly: 795 blocked, 740 gaps, 214 modules (127, 59, 37, 31, 30, 19, …), strict RNP 129
  / RNC 128, classes RNP 101 / RNC 100 / single reversible 142 / isolated 28 / DNP 21 / UNC 78 / DNP+UNC 164,
  624 blocked reactions explained structurally.
- **Solver check.** Open condition: GLPK and Gurobi (1e-9) give the same 756 blocked reactions; the only
  reaction within 1e-4 of zero is R674 (protoporphyrinogen oxidase, max 6.1e-5), active in both. Task A
  medium: GLPK 887, Gurobi 882; the 5 extra GLPK calls (R1146, R1147, R1343, R1344, R970: thiamine thiazole /
  pyrimidine synthesis and the glycolaldehyde by-product) carry exactly 3.75e-7 in Gurobi, which is the TPP
  biomass coefficient (1.2e-6) times the growth rate (0.3124): a tolerance fringe, not a structural block. On
  the medium 61 active reactions sit between 1e-9 and 1e-4 (porphyrin/heme 12 at 1.6–3.1e-7, thiamine 8 at
  3.7e-7, lipoylation and mtFAS 16 at 3–6e-6, iron 2 at 9.8e-6, CoQ9 and polyprenyl 15 at 3.1e-5, SAH
  transport 9.4e-5): the routes that serve the small biomass coefficients. This is the tolerance problem
  recorded in STATE.md (2026-10-22, open problem 1); it does not affect the open-condition census.

## 2. The gap census (open condition)

| Quantity | a1237821 open | a1237821 Task A + PO1f | a56f3649 open (reference) |
| --- | ---: | ---: | ---: |
| Reactions | 2,316 | 2,316 | 2,306 |
| Blocked reactions (GLPK) | 756 (32.6 %) | 887 (Gurobi 882) | 795 (34.5 %) |
| of which closed (0,0) | 69 | 70 (R612 closed by the profile) | 68 |
| of which exchanges / demands | 42 | 92 | — |
| of which open internal reactions with a GPR | 416 | — | — |
| of which open internal reactions without a GPR (orphans) | 229 | — | — |
| Metabolites | 1,882 | 1,882 | 1,875 |
| Gap metabolites (eq. 7) | 701 (37.2 %) | 781 | 740 |
| Strict row scan, RNP / RNC (isolated counted in both) | 128 / 131 | 182 / 131 | 129 / 128 |
| RNP / RNC / single reversible / isolated | 100 / 103 / 136 / 28 | 154 / 103 / 123 / 28 | 101 / 100 / 142 / 28 |
| DNP / UNC / DNP+UNC | 21 / 71 / 160 | 43 / 70 / 149 | 21 / 78 / 164 |
| Gap metabolites with no root in their module (coupling) | 82 | 111 | 106 |
| Blocked reactions explained by root propagation | 618 | 707 | 624 |
| Unconnected modules | 212 | 252 | 214 |
| Largest modules | 127, 37, 37, 32, 30, 19 | 127, 40, 37, 32, 30, 20 | 127, 59, 37, 31, 30, 19 |

- 289 of the 756 blocked reactions are transports, exchanges or diffusions by name; 477 are single-compartment
  reactions and 279 span two compartments. Blocked reactions touch the cytosol in 418 cases, the mitochondrion
  142, the ER 109, the extracellular space 78, the Golgi 72, the nucleus 38, the vacuole 38, the cell envelope 34.
- Gap metabolites per compartment (gap / total): cytosol 228 / 728, mitochondrion 147 / 292, Golgi 59 / 64,
  ER 58 / 103, extracellular 49 / 191, nucleus 34 / 47, cell envelope 31 / 49, vacuole 27 / 27, lipid particle
  23 / 79, mitochondrial membrane 22 / 38, peroxisome 9 / 128, ER membrane 7 / 103, vacuolar membrane 4 / 17,
  Golgi membrane 3 / 16. The vacuole, Golgi, nucleus and cell envelope are mostly or entirely gap space.
- Module sizes: 121 modules of one reaction, 40 of two, 29 of three to five, 9 of six to ten, 12 of 11–50 and
  one of 127. The 161 modules of one or two reactions hold 201 blocked reactions; the 13 modules above ten
  reactions hold 375.
- Pooled roots (invisible to the per-compartment scan): glycolate (m2015 [cy], m2016 [pe]) is never produced
  in any compartment (pooled RNP; a source unblocks R2172–R2175, 4/4), and 4-hydroxy-2-oxoglutarate (m760
  [cy], m761 [mi], m763 [pe]) is never consumed (pooled RNC; module UM15).

## 3. Root metabolites: the 367 RNP / RNC / single-reversible / isolated gaps

Full list with automatic hints and a one-line reason per root: [root_metabolites.csv](root_metabolites.csv)
(sheet *Roots* of the workbook). By compartment: cytosol 104, mitochondrion 105, Golgi 52, extracellular 21,
vacuole 19, cell envelope 13, nucleus 12, ER 11, mitochondrial membrane 10, ER membrane 7, vacuolar membrane 4,
Golgi membrane 3, lipid particle 3, peroxisome 3. Hints: 124 roots have the same compound active in another
compartment (missing transport candidates), 48 share a duplicate id in their compartment, 66 have no formula,
6 extracellular roots have no exchange reaction, 76 sit next to a closed reaction, 5 are components of the
inactive legacy biomass reactions.

**RNP (100, nothing produces them).**
- 14 mitochondrial tRNAs (tRNA(Arg) … tRNA(Val), mitochondrial), each paired with its RNC aminoacyl-tRNA:
  mitochondrial translation is not in the model, so the 14 mitochondrial aminoacyl-tRNA synthetases (plus
  methionyl-tRNA formylation, UM24) are out-of-scope stubs, like the paper's apoACP case (its UM10).
- 23 mitochondrial inositol-phosphoceramide / mannosyl-inositol-phosphoceramide species (UM5): a mitochondrial
  copy of the complex-sphingolipid hydrolysis set whose substrates are never made in the mitochondrion.
- Compartment copies with no producer and no transport: farnesyl diphosphate [lp] (UM6), NADH [er] (squalene
  epoxidase NAD copy), UMP, dAMP, dGMP, formate, nicotinate, diglyceride, 4-hydroxybenzoyl-CoA, (R)-S-
  lactoylglutathione, 1-pyrroline-3-hydroxy-5-carboxylate [mi], ADP [pe], water [va] (vacuolar water diffusion
  closed), CoA, arachidate, behenate [en], 1-acyl-sn-glycerol 3-phosphate in cytosol / Golgi membrane /
  vacuolar membrane, 3-methyl-2-oxobutanoate [cy] (the cytosolic valine transaminase R796 is closed),
  TRX1 [nu] (no nuclear thioredoxin reductase).
- Cytosolic dead starts: three 3',5'-cyclic nucleotides (phosphodiesterases without a cyclase), indole-3-
  acetamide twice (duplicate species m1933 / m2052), 5-oxoproline, cellobiose, D-galactose (galactokinase
  without galactose uptake), agmatine, dopamine, 3,4-dihydroxyphenylethyleneglycol twice (duplicates),
  methanethiol, pantetheine and N-pantothenoylcysteine, pimeloyl-CoA (biotin synthesis entry), dolichol,
  isoamyl acetate, L-methionine S-oxide, a second "Thioredoxin" species (m2033, besides TRX1), ADP-ribose,
  4-(phosphonooxy)-L-threonine, 4-guanidinobutanamide, 5-methyl-THF-triglutamate, 6-PP-IP5,
  N,N'-diacetylchitobiosyl-PP-dolichol [cy] (the N-glycan chain is split between compartments), lipid_ (a
  legacy biomass component).

**RNC (103, nothing consumes them).**
- 14 mitochondrial aminoacyl-tRNAs (see above) and 13 mitochondrial ceramide / mannose-inositol-phosphate /
  inositol-3-phosphate species (UM5).
- Vacuolar amino acids (aspartate, glutamate, isoleucine, histidine, lysine [va]) behind closed vacuolar
  transports; nuclear homocitrate, diphosphate, ADP-ribose and TRX1-disulphide; envelope AMP, diphosphate,
  chitosan, sn-glycero-3-phosphoethanolamine and -serine; ER NAD, hexadecanal, 2-hydroxyhexadecanal,
  O-phosphoethanolamine, phytoceramide, hexadecanoic acid, chitobiosyl-PP-dolichol [er]; mitochondrial UTP,
  TDP (thiamine diphosphate transported in but never used stoichiometrically), (S)-2-methylbutanoyl-CoA,
  nicotinate D-ribonucleotide; docosaprenyl diphosphate [lp]; UDP-galactose [go]; the ER mannose-inositol-
  phosphates m663 / m674 (UM1).
- Three extracellular dicarboxylates (decanedioate, dodecanedioate, hexadecanedioate) that are exported but
  have no exchange reaction.
- Cytosolic dead ends: siroheme (STATE.md 2026-10-22 open problem 3, confirmed), biotinyl-5'-AMP (biotin is
  never transferred to apo-carboxylases), ACP (the acyl-ACP block), myristoleate and octadecadienoate,
  TTP, UDP-D-xylose, ergosterol glucoside, 1-methylnicotinamide, N-acetylputrescine, 4-guanidinobutanoate,
  S-methyl-L-methionine, quinaldic acid, dehydroalanine, 2-deoxyribose 5-phosphate, 2'-deoxyinosine
  5-phosphate, IP3, aminoacetaldehyde, L-alpha-formylglycine, gamma-glutamyl-alanine, N(pros)-methyl-
  histidine, Ap4A and Gp4G, indole-3-pyruvate, (R)-acetoin, thioredoxin disulfide (m2034), homovanillate /
  3-methoxy-4-hydroxymandelate / MHPG (catecholamine-degradation stubs), 2-hydroxyphenylacetate,
  phenylacetic acid (duplicate m1932 / m2026), salicyl alcohol, tetrahydropteroyltriglutamate,
  4-hydroxy-L-threonine, the diphthamide intermediate, 1-monoglyceride in four compartments.

**Single reversible reaction (136) and isolated (28).** 48 of the single-reversible roots and 10 of the isolated
ones are Golgi complex-sphingolipid species of UM1 (behind the 20 closed Golgi synthases, reached only by
transports); 18 are extracellular compounds whose only reaction is an exchange (xylitol, maltose, sucrose,
trehalose, pectin, galacturonate, decane, dodecane, adenosine, inosine, zymosterol, hexadecanal, …); the
rest are pool components of the mitochondrial-membrane and cell-envelope fatty-acid / acyl-CoA pools (UM10,
UM8), nuclear SAM / SAH / H2O2 transport-only species, vacuolar glycogen, trehalose, glucose, phosphate and
polyamines behind closed vacuolar transports, and single-step dead ends such as D-ribulose (R2190),
trans-aconitate, itaconate / itaconyl-CoA, glycerate, D-serine, formaldehyde, 2,3-bisphosphoglycerate.

## 4. The unconnected modules

[unconnected_modules.csv](unconnected_modules.csv) lists all 212 modules with roots, closures, cause
breakdown, best single fix and a narrative for the 25 largest. Network drawings:
[fig_UM1_network.png](fig_UM1_network.png), [fig_UM2-10_networks.png](fig_UM2-10_networks.png).
All "n/N" figures below are LP evidence from the fix tests ("adding X unblocks n of the module's N blocked
reactions"), not curation recommendations.

| Module | Blocked reactions | Gap metabolites | Roots RNP / RNC / single / isolated | Content | Main primary cause | Best single fix | All elementary fixes together |
| --- | ---: | ---: | --- | --- | --- | --- | ---: |
| UM1 | 127 | 110 | 1 / 3 / 37 / 11 | complex sphingolipids (ER / Golgi / cytosol) | several roots (>1 fix) | sink on m663 [er] → 28/127 | 118/127 |
| UM2 | 37 | 27 | 1 / 3 / 0 / 7 | vacuolar transports + dipeptides | several roots (>1 fix) | — | 37/37 |
| UM3 | 37 | 32 | 2 / 2 / 4 / 1 | riboflavin / FAD + legacy biomass + BCAA fragments | cofactor pseudo-gap | sink on m1992 [cy] → 13/37 | 31/37 |
| UM4 | 32 | 27 | 0 / 4 / 3 / 0 | acyl-ACP block + C10 dicarboxylate | several roots (>1 fix) | sink on m1937 [ex] → 5/32 | 11/32 |
| UM5 | 30 | 43 | 30 / 13 / 0 / 0 | mitochondrial IPS phospholipase copy | several roots (>1 fix) | — | 30/30 |
| UM6 | 19 | 21 | 1 / 1 / 0 / 0 | lipid-particle dolichol backbone | several roots (>1 fix) | — | 0/19 (targeted IPP + FPP import: 19/19) |
| UM7 | 17 | 15 | 0 / 3 / 1 / 0 | nuclear NAD / homocitrate / transports | several roots (>1 fix) | sink on m628 [nu] → 4/17 | 16/17 |
| UM8 | 17 | 17 | 0 / 2 / 3 / 0 | envelope phospholipase B products | several roots (>1 fix) | — | 17/17 |
| UM10 | 13 | 19 | 0 / 0 / 8 / 0 | mitochondrial-membrane cardiolipin remodelling pools | several roots (>1 fix) | — | 12/13 |
| UM9 | 13 | 10 | 0 / 0 / 0 / 0 | de novo folate synthesis | cofactor pseudo-gap | — (THF drain 13/13) | — |
| UM11 | 11 | 9 | 1 / 1 / 0 / 0 | nuclear ribonucleotide reductase / thioredoxin | several roots (>1 fix) | — | 11/11 |
| UM12 | 11 | 7 | 0 / 0 / 1 / 0 | Kennedy CDP-choline branch | closed reaction | open R1726 → 11/11 | 11/11 |
| UM13 | 11 | 11 | 1 / 1 / 2 / 0 | biotin synthesis and uptake | dead-end root | sink on m347 [cy] → 6/11 | 10/11 |
| UM14 | 10 | 8 | 2 / 0 / 0 / 0 | CoA synthesis from pantothenate | cofactor pseudo-gap | — (CoA drain 5/10) | 0/10 |
| UM15 | 10 | 8 | 1 / 0 / 0 / 0 | 4-hydroxyproline degradation | coupling (no root) | — (sinks on both ends 7/10) | 0/10 |
| UM16 | 10 | 11 | 3 / 2 / 0 / 0 | envelope VLCFA activation | several roots (>1 fix) | — | 10/10 |
| UM17 | 9 | 8 | 2 / 0 / 1 / 0 | tryptophan → indole-3-acetate | dead-end root | sink on m1260 [mi] → 7/9 | 9/9 |
| UM18 | 8 | 8 | 0 / 1 / 1 / 0 | ER dolichol-linked N-glycan steps | dead-end root | sink on m1969 [er] → 6/8 | 8/8 |
| UM19 | 8 | 5 | 0 / 0 / 1 / 0 | Kennedy CDP-ethanolamine branch | closed reaction | open R1725 → 8/8 | 8/8 |
| UM20 | 6 | 4 | 0 / 0 / 0 / 0 | nuclear inositol polyphosphates | coupling (no root) | — (IP4 sinks 4/6) | — |
| UM21 | 6 | 6 | 1 / 0 / 0 / 0 | mitochondrial H2O2 detoxification | coupling (no root) | — (H2O2 [mi] source 5/6) | — |
| UM22 | 6 | 4 | 0 / 0 / 0 / 0 | cytosolic dehydrodolichol-PP steps | cofactor pseudo-gap | — (dolichyl-P drain 6/6) | — |
| UM23 | 5 | 3 | 0 / 0 / 0 / 0 | threonine / allo-threonine fragment | coupling (no root) | — (2-amino-3-oxobutanoate source 5/5) | — |
| UM24 | 5 | 5 | 2 / 1 / 0 / 0 | mitochondrial fMet-tRNA | missing transport | sink on m533 [mi] → 3/5 | — |
| UM25 | 5 | 2 | 0 / 0 / 0 / 0 | vitamin B6 interconversions | cofactor pseudo-gap | — (PLP drain 5/5) | — |

"Best single fix" is the elementary fix test with the largest count in the module; "all elementary fixes
together" applies every root sink, transport link and closure opening of the module at once; "—" means no
automatic test applied (no root with an active copy elsewhere, no closure) or none unblocked anything.

**Escher maps.** The 20 modules of the table, and four cross-module *problem-class* maps, are drawn as Escher maps in the offline viewer
[escher/gap_modules_escher.html](escher/gap_modules_escher.html) (module and scenario drop-downs; node colour =
gap class; reaction colour and label value = the census state, the number of single fixes that unblock the
reaction, the all-fixes test, the module's targeted test or the cofactor drains; click a reaction for its
equation, GPR, bounds and value). The class maps gather the reactions of one primary cause from every module into
labelled sections: cofactor pseudo-gaps by cofactor family with the tested drains drawn as TEST reactions, the 69
closures by group with the 18 reactions an opening frees, the missing transports and exchanges by the root's
compartment with the tested links drawn, and the dead-end roots by root class with the tested sinks drawn; every
TEST reaction carries the number of reactions it unblocks. Static renders of every map (SVG and PNG) are in
the artifacts folder under `escher/renders/`; three are kept here: [UM9 folate](escher/UM9.png),
[UM12 Kennedy CDP-choline](escher/UM12.png), [UM6 dolichol backbone](escher/UM6.png). How the maps are built,
what the layout means (automatic; repeated nodes are the same pool) and the headless check of all 100
module × scenario views are described in [escher/README.md](escher/README.md). The tooltip of every map entity
links to its entry in [escher/gap_dataset.html](escher/gap_dataset.html), the census dataset (all metabolites and
blocked reactions with their facts and LP evidence). No new LP solves were run for the maps.

**What the large modules are.**

1. **UM1 (127 reactions, 110 gaps; ER, Golgi, cytosol and both membranes): complex sphingolipids.** The
   IPC → MIPC → M(IP)2C synthesis and turnover is present three times with no connected route: the ER
   "IPS phospholipase C" set (R427–R446, R477–R486) hydrolyses the complex sphingolipids into ceramides and
   the mannose-inositol-phosphates m663 / m674 [er], which nothing consumes; the Golgi synthase set
   R1815–R1834 is closed (0,0) in the source model, leaving the Golgi species and the Golgi-membrane
   diglyceride (m1704, in 22 closed reactions) isolated, and opening the 20 closures alone gives them no
   flux (closure check) and unblocks only 5/127; the cytosolic synthase set R2143–R2162 makes IPC / MIPC
   species m1995–m2006 [cy] that nothing consumes (draining them unblocks 16/127). Sinks on all 52 roots
   unblock 94/127, every elementary fix together 118/127. The current `biomass_C` has no complex-
   sphingolipid term, and the legacy lipid pseudoreaction R1371 that drained "complex sphingolipid"
   (m1189 / m1190, UM3) is closed: the module is largely a drain problem (the paper's pseudo-gap logic
   applied to a membrane lipid) on top of the inherited closures. Yeast biology (recalled, uncited): IPC
   synthase AUR1, MIPC synthases CSG1 / CSH1 and M(IP)2C synthase IPT1 act in the Golgi; ISC1 hydrolyses
   complex sphingolipids.
2. **UM2 (37; vacuole): the inherited vacuolar closures.** Fifteen vacuolar amino-acid / polyamine
   transports (R869–R927), the V-ATPase R795, H+ diffusion R1163 and water diffusion R1363 are (0,0) in the
   source model; the four dipeptide uptake / hydrolysis sets (Gly-Asp, Gly-Glu, Ala-Gly, Gly-Pro;
   R2018–R2040) hang on vacuolar H+ and water. Opening the 17 closures unblocks 34/37, everything 37/37.
   These closures have a recorded history (STATE.md 2026-09-10/11); not a new gap.
   **Decision 2026-10-28 (user, Option A0): the closures stay.** They are inherited from iYali4 (2016) with no
   recorded reason. Every opening gives only futile ATP cycles on Task A, with no growth, knockout or energy
   change. R1163 in its written direction makes a no-uptake loop. See
   [um2_vacuolar_closures_20261028/REPORT.md](../um2_vacuolar_closures_20261028/REPORT.md).
3. **UM3 (37): three unrelated parts joined by cytosolic FAD/FADH2 and the legacy biomass reactions.**
   (a) Riboflavin → FMN → FAD synthesis and transport (R11, R12, R34, R378, R697–R700, R306, R1125,
   R1315 / R1316, R847): a coenzyme pseudo-gap, a FAD + FMN drain unblocks 18/37. (b) Unique components of
   the inactive legacy biomass reactions R1372 (open but blocked), R1387 / R1710 and R1371 (closed): lipid_,
   complex sphingolipid, (1-6)-beta-D-glucan [en] (made by R5, consumed only by the old biomass), ergosterol
   ester, fatty-acid pool species. (c) The mitochondrial UMP / UTP transport pair R938 / R939 and BCAA /
   propanoate degradation fragments (R1969, R2082, R2114–R2118, R2140–R2142). Root sinks 28/37.
4. **UM4 (32): the accepted acyl-ACP export block.** R2216–R2241 (26 reactions) are blocked by the user's
   2026-10-21 decision ("blocked in every medium is acceptable"), with ACP [cy] and the hydrolase products
   myristoleate / octadecadienoate as dead ends; the module also holds the C10 dicarboxylate branch
   (R2253–R2256, R1885, R1117) whose product decanedioate [ex] has no exchange reaction (likewise
   dodecanedioate, UM30, and hexadecanedioate, UM35). Root sinks unblock 11/32; the block stays by design.
5. **UM5 (30, all mitochondrial): a pathway copy with neither input nor output.** The mitochondrial IPS
   phospholipase C set (R447–R476) has 23 RNP substrates and 13 RNC products and no transport to the ER
   copies (transport test 0/30, product sinks 0/30).
6. **UM6 (19, lipid particle): the dolichol backbone.** The cis-prenyltransferase chain R217–R235 lacks its
   precursors in the lipid particle (farnesyl diphosphate RNP, IPP not supplied) and a consumer of its
   product (docosaprenyl diphosphate RNC). Root sinks and transports alone unblock 0/19; importing IPP and
   FPP from the cytosol plus a product sink unblocks 19/19. Two more fragments belong to the same route:
   dolichol [cy] (RNP, dolichol kinase R287, UM167) and the cytosolic dehydrodolichol-diphosphate steps
   (UM22, unblocked by a dolichyl-phosphate drain, 6/6): dolichol synthesis is cut into three unlinked
   fragments and dolichyl phosphate is a conserved carrier that nothing drains.
7. **UM7 (17, nucleus): missing nuclear transports.** Homocitrate [nu] (no export to the mitochondrion
   where lysine synthesis continues), diphosphate [nu] and ADP-ribose [nu] (no consumer or transport)
   around the nuclear NAD nucleosidase / NAD synthase / homocitrate synthase copies and the closed nuclear
   acetyl-CoA synthetase R92. A homocitrate nucleus ↔ mitochondrion link plus PPi and ADP-ribose sinks
   unblocks 14/17.
8. **UM8 (17, cell envelope): phospholipase B products with no exit.** sn-glycero-3-phosphoethanolamine
   and -serine [en] are never consumed or transported; the envelope fatty-acid pool needs laurate,
   myristate and linoleate copies that occur only in the pool. Root sinks 17/17, transports 13/17.
9. **UM9 (13): de novo folate synthesis is a coenzyme pseudo-gap.** GTP cyclohydrolase I → dihydroneopterin
   → dihydropteroate → dihydrofolate (R377, R274, R277 / R278, R273 and transports) has no root. Opening the
   closed DHFR reactions R1950 / R1951 unblocks 0/13 (dihydrofolate is reduced by the active R271); a THF
   drain alone unblocks 13/13. The THF / DHF pool is conserved and not in `biomass_C`.
10. **UM10 (13, mitochondrial membrane): pools without transports.** The cardiolipin remodelling cycle
    (R1723, R1724, R1727 / R530 — the duplicate LPCAT pair) is fed by the pools xPOOL_AC_MM / xPOOL_FA_MM,
    which require acyl-CoA and fatty-acid copies that have no transport into the compartment. Transports to
    the cytosolic copies unblock 12/13; sinks 0/13.
11. **UM11 (11, nucleus):** nuclear ribonucleotide reductase copies with a nuclear thioredoxin that is never
    regenerated (TRX1 [nu] RNP / TRX1-disulphide RNC). Linking to the cytosolic pair unblocks 11/11.
12. **UM12 (11) and UM19 (8): the Kennedy pathway behind the inherited closures R1726 and R1725.** Choline
    and ethanolamine uptake, kinases, cytidylyltransferases and CDP-choline / CDP-ethanolamine transports.
    Opening R1725 + R1726 unblocks 19/19; on the Task A medium growth goes 0.312421 → 0.312472 and the
    energy-generating-cycle LP stays at 0 (closure check, re-done on 2026-10-27 with the maintenance bound
    relaxed; see section 6).
13. **UM13 (11): biotin.** Pimeloyl-CoA has no producer, polysulphur and the BioA by-product occur only in
    one reaction each, biotinyl-5'-AMP is never transferred to apo-carboxylases, and biotin is not drained,
    so even biotin uptake (R824) is blocked although biotin is in the medium: a biotin drain unblocks only
    2/11, root sinks 10/11. A coenzyme pseudo-gap with missing steps at both ends.
14. **UM14 (10): CoA synthesis from pantothenate** (PanK, PPCS, PPCDC, PPAT, DPCK and the mitochondrial
    copies): a CoA drain alone unblocks 5/10; root sinks 0/10. Coenzyme pseudo-gap.
15. **UM15 (10): 4-hydroxyproline degradation** with no entry (1-pyrroline-3-hydroxy-5-carboxylate [mi]
    RNP, no hydroxyproline oxidase) and no exit (4-hydroxy-2-oxoglutarate, pooled RNC, no aldolase). Sinks on
    both ends 7/10; the entry alone 0/10.
16. **UM16 (10, cell envelope):** very-long-chain fatty-acid activation with CoA, AMP and diphosphate copies
    that cannot enter or leave the envelope compartment; transports 6/10, everything 10/10.
17. **UM17 (9):** tryptophan → indole-3-acetate branches with duplicate indole-3-acetamide species and no
    consumer of indole-3-acetate; root sinks 9/9.
18. **UM18 (8) + UM184 (1) + R141: the dolichol-linked N-glycan precursor chain is split into three
    unlinked compartment copies** (ER: R2105 / R2106 → m1969 [er]; cytosol: R173, m344 → m345; Golgi: R141,
    m289 → m291). Linking the copies and draining m291 unblocks 8/10; the dolichyl-P-glucose branch
    (R2138 / R2139) stays blocked (no ALG6 / ALG8-type steps).
19. **UM20 (6), UM21 (6), UM23 (5), UM24 (5), UM25 (5):** nuclear inositol-polyphosphate kinases without a
    consumer of IP4 (sinks 4/6); mitochondrial H2O2 detoxification without a mitochondrial H2O2 source
    (a H2O2 source unblocks 5/6); a threonine / allo-threonine degradation fragment with no entry step
    (a 2-amino-3-oxobutanoate source unblocks 5/5); mitochondrial fMet-tRNA formation (out of scope); and
    vitamin B6 interconversions, a coenzyme pseudo-gap (a PLP drain unblocks 5/5).

**The 161 small modules (one or two reactions, 201 blocked reactions)** fall into six patterns: 54 exchange /
transport stubs for compounds that nothing inside the cell makes or uses (sugars, alkanes, nucleosides,
sterols, and intracellular transport-only pairs such as nuclear SAM / SAH / H2O2 and the mitochondrial NAD
transports R903–R906 paired with dAMP / dGMP); 18 modules created by inherited closures (R1040, R859, R1772,
R746, R1950, R91, R1297, R1349, R1138, R26 / R520, R1023 / R1027, R21 / R796, R606 / R1899, R1600 / R1720,
R696 / R853, R573, R491, R1342); 13 mitochondrial aminoacyl-tRNA synthetases (out of scope); pseudo-gap-like
single steps whose product is a conserved pool (NAD kinases R565 / R566 / R571 / R572, NMNAT R578, the
R2190 D-ribulose dead end documented on 2026-10-19); dead-end products with no consumer (cyclic-nucleotide
phosphodiesterases, Ap4A / Gp4G synthases, lysoPA phosphatases in three compartments, R65, R287, R777,
R2183, R2097, …); and the glycolate pooled root (UM37).

## 5. Why the reactions are blocked: primary causes and pseudo-gaps

One primary cause per blocked reaction, assigned in this order: closed (0,0) → cofactor pseudo-gap (a drain
on the cofactor unblocks it) → the single elementary fix that unblocks it (transport > exchange > opening a
closure > a sink on a dead-end root > duplicate-id merge) → structurally blocked at several roots (needs more
than one fix) → stoichiometric coupling with no root in the module. Full detail per reaction in
[blocked_reactions.csv](blocked_reactions.csv).

| Primary cause | Blocked reactions |
| --- | ---: |
| Chain disconnected at several roots (no single elementary fix unblocks it) | 290 |
| Dead-end root metabolite (a sink on one root unblocks it) | 181 |
| Cofactor pseudo-gap (a cofactor drain unblocks it) | 72 |
| Closed in the model (0,0) | 69 |
| Stoichiometric coupling / pseudo-gap with no root in the module | 57 |
| Missing transport (a link to the active copy in another compartment unblocks it) | 52 |
| Blocked by a closed (0,0) reaction (opening it unblocks it) | 18 |
| Missing exchange (an exchange for an extracellular root unblocks it) | 17 |
| Total | 756 |

Single elementary fix tests (621; each re-tests only its module's blocked reactions):

| Kind of fix | Tests | Tests that unblock ≥ 1 reaction | Largest count in one test |
| --- | ---: | ---: | ---: |
| Sink (free source/sink) on a root | 367 | 139 | 28 (m663 [er], UM1) |
| Transport link to the active copy of the compound in another compartment | 124 | 40 | 5 |
| Opening one closed (0,0) reaction | 69 | 15 | 11 (R1726) |
| Duplicate-id merge (link to a same-named species in the same compartment) | 55 | 0 | 0 |
| Exchange for an extracellular root | 6 | 6 | 5 |

304 of the 756 blocked reactions are unblocked by at least one single fix; the other 452 need either a
drain (72), more than one fix at once (290 "several roots") or a non-elementary change (57 coupling, 69
closed, less the overlaps counted above). Duplicate-id merges alone never unblock anything (0 of 55), as on
the reference candidate.

**Coenzyme pseudo-gaps (paper Fig. 2).** `biomass_C` drains ATP, the macromolecule precursors, ergosterol and
zymosterol, the five lipid classes, ubiquinone-9, iron, the two lipoyl proteins, TPP and heme a; it drains
no THF, FAD, FMN, riboflavin, PLP, NAD(P), CoA, dolichyl phosphate, siroheme, biotin, SAM, glutathione, ACP,
lipoamide, heme b or heme o. An irreversible drain per undrained cofactor (20 tests; each re-tests all 756
blocked reactions) unblocks 72 reactions in total, and the union of the single drains equals the joint drain:

| Cofactor drained | Reactions unblocked | Modules |
| --- | ---: | --- |
| THF (cy, mi) | 19 | UM9 (13), UM33 (4), UM40 (2) |
| FAD (cy, mi) | 16 | UM3 |
| FMN (cy, mi) | 11 | UM3 |
| riboflavin (cy, mi) | 9 | UM3 |
| PLP (cy) | 8 | UM25 (5), UM88, UM212 |
| NADP (all compartments) | 8 | UM7 (2), the NAD kinases UM98 / UM113 / UM131 |
| CoA (all compartments) | 7 | UM14 (5), UM7, UM10 |
| dolichyl phosphate (cy, er) | 6 | UM22 |
| NAD (all compartments) | 4 | UM7 (2), UM115, UM180 |
| pantetheine 4'-phosphate (cy, mi) | 4 | UM14 |
| siroheme (cy) | 3 | UM47 |
| biotin (cy) | 2 | UM13 |
| heme b, heme a [mi], heme o, SAM, glutathione, holo-ACP, lipoamide | 0 | — |
| all of the above together | 72 | |

On the reference candidate a56f3649 the same test gave 78 (heme a 18, heme b 15, siroheme 12 among them);
the heme a biomass term of 2026-10-22 made the porphyrin / heme route active, so the heme a, heme b and most
siroheme pseudo-gaps are gone. The paper's remedy for a pseudo-gap is a biomass term or a drain for the
cofactor; whether iYali26 should carry such terms (and with which coefficients) is a curation decision that
needs evidence and authorization, not a result of this census.

## 6. The 69 closed reactions

[closed_reactions.csv](closed_reactions.csv). 68 of the 69 are (0,0) in the immutable source `iyali26.xml`,
in `iyli21.xml`, in the legacy `model.xml` and in a56f3649; only R1342 (TMP kinase) was closed recently, by
the documented TPP step of 2026-10-22. Groups: the 20 Golgi IPC synthase / inositolphosphotransferase
reactions (UM1); 17 vacuolar reactions (R795 V-ATPase, R1163, R1363, R869–R898, R918 / R925 / R927; UM2);
the Kennedy pair R1725 / R1726; the DHFR pair R1950 / R1951; the acetyl-CoA synthetase copies R91 (mito) /
R92 (nucleus); the YALI1D01431g transaminase / decarboxylase set R26, R491, R520, R796; the legacy biomass
R1371 / R1387 / R1710; the ABC / vacuolar transports R696, R746, R853; and R17, R352, R573, R859, R1027,
R1040, R1138, R1297, R1349, R1720, R1772, R1775, R1795, R1843, R1899.

Closure checks on the Task A medium + PO1f (growth as is 0.312421; the energy-generating-cycle LP maximises
ATP hydrolysis, R_NTP1, with every exchange closed to uptake and the maintenance lower bound relaxed):

| Opened | Growth | Max ATP hydrolysis with no uptake (status) | FVA 99.9 % when opened |
| --- | ---: | ---: | --- |
| R1725 + R1726 (Kennedy) | 0.312472 | 0 (optimal) | R1725 [0, 0.00105], R1726 [0, 0.02235] |
| R1950 + R1951 (DHFR) | 0.312492 | 0 (optimal) | R1950 [0, 0.00516], R1951 [0, 0] |
| R1815–R1834 (20 Golgi synthases) | 0.312421 | 0 (optimal) | all 20 stay at [0, 0] |
| R91 + R92 (acetyl-CoA synthetase mito / nucleus) | 0.312421 | 0 (optimal) | R91 [0, 0.07065], R92 [0, 0] |
| R1371 (lipid pseudoreaction) | 0.312421 | 0 (optimal) | [0, 0] |

**Correction (2026-10-27).** The first version of this table reported "max ATP with no uptake = 0" from the
skill's `closure_check.py`, which writes an infeasible LP's NaN objective as 0.0. That LP is always infeasible
on this model, because the maintenance reaction `xMAINTENANCE` has a lower bound of 7.8625 that cannot be met
without uptake; the zeros were masked infeasibles and proved nothing (found by the cofactor biomass-terms
session, confirmed by its auditor). The check was redone with every exchange closed to uptake, every positive
lower bound relaxed to 0 (only `xMAINTENANCE`) and R_NTP1 maximised, Gurobi at 1e-9: all ten LPs (as is and
opened, five groups) are optimal with objective 0, so the "no energy-generating cycle" conclusion now rests on
a valid LP (`egc_check_relaxed.json` and `egc_check.py` in the artifacts folder). The growth and FVA columns
were not affected.

Opening each closure alone, with every exchange open (fix tests, module reactions only): 15 of the 69
unblock something: R1726 11/11 (UM12), R1725 8/8 (UM19), R17 4/4 (UM34, the mitochondrial DAHP synthase
copy), R1951 2/3 (UM40), R796, R1027, R1720 and R1899 2/2 each (UM70, UM66, UM80, UM74), and R91, R491,
R573, R1040, R1342, R1772, R1950 their own single reaction. The other 54 closures (the 20 Golgi synthases,
the 17 vacuolar reactions, R26, R92, R352, R520, R696, R746, R853, R859, R1138, R1297, R1349, R1363, R1371,
R1387, R1710, R1775, R1795, R1843) stay blocked when opened alone: they sit inside modules with other roots.

Opening a closure is LP-safe in these five cases (no free energy on the corrected check, growth unchanged or
+7e-5), but the R795 /
vacuole group was closed for reasons that a separate task is still tracing (STATE.md 2026-09-10), and none of
these results says what the original authors intended.

## 7. Changes since the reference candidate a56f3649 (open condition)

[build_comparison_vs_a56f3649.csv](build_comparison_vs_a56f3649.csv). 48 reactions became active and 9
became blocked; the 10 reactions added between the two candidates (lipoylation R2286–R2290, the GPI links,
R_DM_GPI, R_GPIGNT, the ER transports) are all active. Unblocked: the porphyrin / heme / iron route (18:
R61, R667, R409, R792, R790, R240, R1299, R674, R302, R384, R383, R1709, R810, R987, R988, R863, R1189,
R2193; the heme a and iron biomass terms), the GPI-anchor chain (7: R113, R371–R376; the GPI links), the
mitochondrial fatty-acid synthesis cycle (13: R78, R1394–R1405; the lipoylation module), the thiamine
thiazole / pyrimidine route with its glycolaldehyde by-product (7: R1343, R1344, R970, R752, R934, R1146,
R1147; the TPP requirement) and three transports (R1116, R1251, R1290). Newly blocked: R2190 (D-ribulose
product, documented on 2026-10-19) and eight acyl-ACP export / hydrolysis reactions (R2225, R2226,
R2230–R2235; the accepted acyl-ACP block). Blocked count 795 → 756, gap metabolites 740 → 701, modules
214 → 212; the 59-reaction second module of a56f3649 (porphyrin + GPI + mtFAS) no longer exists.

## 8. Medium dependence (Task A medium + PO1f)

131 reactions are blocked on the medium but not with every exchange open (`medium_only_blocked.csv` in the
artifacts folder and sheet *Medium_only_blocked* of the workbook): 95 are attributed to 41 metabolites that become RNP only on the medium (40
extracellular nutrients absent from it — FMN, TDP, TMP, GABA, 5-aminolevulinate, allantoin, choline,
nucleosides, D-mannose, D-ribose, glutathione, lactate, SAM, raffinose, salicin, … — and orotidine
5'-phosphate [cy], because R612 is closed in the PO1f profile); 31 are reversible dead-end chains of absent
carbon sources (D-arabinose, D-xylose, L-arabinitol / L-arabinose, L-sorbose / L-glucitol, D-arabitol,
glucosamine, hexadecane: uptake plus the first steps, blocked because the chain has no net source); and 5
are the tolerance fringe of section 1. Gap metabolites 701 → 781, strict RNP 128 → 182, modules 212 → 252.
The medium adds nutrient-absence gaps; it removes none of the open-condition gaps.

## 9. What this means, and what it does not

- **Provable from the model structure** (no solver needed): the RNP / RNC / isolated / single-reversible
  classes, the pooled roots and the propagation (618 of the 756 blocked reactions follow from the roots).
- **Observed in this run** (LP results, GLPK, cross-checked with Gurobi for the blocked sets): the blocked
  sets, the module partition, the fix-test and drain-test counts, the closure checks.
- **Inferred**: the module narratives (which pathway a module is, which end is missing). The yeast biology
  cited in them is recalled, not verified in this session.
- **Untested**: whether any fix is biologically right for *Y. lipolytica*; gene identities; the original
  reasons for the inherited closures.
- Compared with the paper's *B. cuenoti* model (69 of 419 reactions blocked, 58 of 364 metabolites gaps,
  10 modules), a compartmentalised eukaryotic model has far more gaps (33 % / 37 %, 212 modules), and most of
  them are of kinds the paper's small model did not show: compartment copies without transports (124
  roots), conserved-pool pseudo-gaps (72 reactions), inherited closures (69), legacy biomass leftovers and
  out-of-scope processes (mitochondrial translation).
- Candidate fix classes, in the paper's vocabulary and for a later, authorized curation line (none applied
  here): biomass terms or drains for the undrained cofactors and for complex sphingolipids; transports for the
  nuclear, envelope, mitochondrial-membrane and lipid-particle copies; exchanges for the three exported
  dicarboxylates; consumers or removal for the dead-end products; merges for the duplicate species; a decision
  on the legacy biomass reactions R1371 / R1372 / R1387 / R1710; and the inherited closures, which need their
  original rationale before anything is opened.

## 10. Method notes

- The scan is the skill's `gapfind.py` (FVA with `fraction_of_optimum = 0`, `open_exchanges = False` after the
  condition is applied; blocked = |min| and |max| below cobra's zero cutoff, 1e-7 for GLPK). Propagation
  removes the reactions of every root, re-scans, and repeats; a metabolite left with a single reaction is
  classed "DNP+UNC (single reaction left)".
- Fix tests add one elementary reaction (a reversible sink on a root, a reversible link to the active copy of
  the same compound in another compartment, a duplicate-id link, an opened closure, or an exchange) and re-run
  FVA on the module's blocked reactions only; combined tests apply all ops of one kind per module.
- Possible improvements, and how to test them: (i) a growth-scaled tolerance or a Gurobi-only run for the
  medium condition would remove the 5-reaction fringe (compare blocked sets at 1e-7 and 1e-9 as here);
  (ii) the propagation could treat pooled roots as roots in the first pass (re-run `propagate` with the two
  pooled compounds added and compare the explained count, 618); (iii) fix tests are single-op; a minimal-set
  search per module (MILP, as in the paper's gap-filling references) would find the smallest combination that
  unblocks a module (test on UM1 and UM3, where single ops explain less than the "everything" test).

## 11. Limits

- No model change, no new reference; every fix test is a hypothesis for a curation line.
- The Task A medium count carries a 5-reaction solver fringe; the heme / thiamine / lipoylation / CoQ routes
  run within 1e-7–1e-4 of zero on that medium, so medium-condition blocked sets on this candidate should be
  read with a 1e-9 solver (STATE.md 2026-10-22, open problem 1).
- The narratives name pathways by the model's own reaction and metabolite names; gene names were not
  verified. No independent audit of this report has been run.
- Not covered: single-gene knockouts, other media, the lipid-unlump line (owned by a colleague).
