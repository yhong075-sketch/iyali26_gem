# Why closing R1172 blocked nine reactions in 7ccbb568, and the options (2026-11-16)

Scope, authorization and starting state: [TASK.md](TASK.md). Computation and documentation only: no model XML,
curated data, builder code, benchmark or test changed; nothing committed. Tables in [tables/](tables/), scripts in
[scripts/](scripts/), figures in [figure/](figure/); raw outputs and logs in
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/r1172_side_effect_20261116/`.

## Answer in short

- **Reproduced (LP, both solvers, both conditions).** On `7ccbb568`, R25, R552, R95, R2074, R2075, R2076, R2077,
  R2112 and R834 are blocked with every exchange open and on Task A + PO1f. With R1172 reopened in memory at
  (−1000, 1000), all nine carry flux. Gurobi (1e-9) and GLPK agree; every LP was optimal.
- **Mechanism (static identity plus LP).**
  - **The identity (proven).** No reaction inside the mitochondrial matrix makes or destroys the CoA moiety. 34
    matrix reactions conserve it, and only 9 transports move it across the inner membrane. So at every steady
    state, as much CoA moiety leaves the matrix as enters it. This holds whatever R1172's bounds are.
  - **The nine all need an exit.** Each of the nine needs CoA moiety to enter the matrix:
    - R25 turns imported propionyl-CoA into free matrix CoA (R95 and R552 follow it);
    - the cytosolic methylmalonyl route ends by importing succinyl-CoA on R2077;
    - R834 imports free CoA.
  - **R1172 was the only working exit.** It took free matrix CoA out as HMG-CoA (CoA → acetyl-CoA →
    acetoacetyl-CoA → HMG-CoA). Every other exit is closed:
    - R2077, R2083 and R1291 can only import, because the species on the far side has no other consumer, or no
      producer;
    - R2025 can only pass on propionyl-CoA that came in on R2113;
    - R2118 needs a cytosol-to-matrix FAD flow through R2117 that steady state cannot sustain.
  - **The LP confirms it.** The maximum net propionyl-CoA import into the matrix is 0 as the model is. With R1172
    reopened it equals R25's maximum (512.9 open, 1.558 Task A). Every route that carries one of the nine pairs one
    CoA moiety in with one HMG-CoA out through R1172.
- **Correction to the starting reading.** The total CoA moiety is never exchanged net, in any model state; that
  part is an identity. What R1172 provided was the only way to take free matrix CoA out of the matrix, so any CoA
  that arrived as propionyl-CoA, succinyl-CoA or free CoA could leave again. At the WT optimum before the closure,
  R834 = R1172 = R83 = R412 = R410 = 0.0721 and the other eight carried 0 (FVA maximum 0.009). R834 was returning
  the CoA of a route that moved acetyl units out of the matrix as HMG-CoA.
- **Genes.** Three genes on the nine are consensus-essential, and all three are FN at 15 % (KO/WT 1.0):
  `YALI1F20788g` (R834), `YALI1F03803g` (R95) and `YALI1F39620g` (R552). R2075, R2076 and R2077 have no gene.
- **Options tested in memory** (Task A, R1172 kept closed). Every option frees all nine. Every energy check is 0
  (optimal). None of 26 targeted knockouts changes class at 15 %. What differs is the route each option uses:

  | Option | Methylcitrate trio | Methylmalonyl five | R834 |
  | --- | --- | --- | --- |
  | O1: matrix 2-oxobutanoate dehydrogenase (inferred PDH rule) | inside the matrix, no CoA crosses (up to 7.60) | only through side route A: R2112 backwards plus gene-less R2113/R2077 | only through side route B: CoA in, propionyl-CoA out on gene-less R2025 |
  | O1 + R2075 closed (the package the user authorized on the other line) | as O1 | blocked | as O1 (3.31) |
  | O3: propionyl-carnitine shuttle (3 new reactions, inferred rules) | no CoA crosses (5.80) | side route A | a futile loop (CoA in, propionyl-CoA out on R2113) moving 4 H+ into the matrix per turn, up to 68.7 |
  | O4: matrix CoA growth dilution (placeholder 1.9e-5) | at most 5.9e-6, fed through gene-less acyl-CoA imports | 5.9e-6 | 5.9e-6; its gene stays non-essential because R2113 replaces it |

- **Leading option (for the user to decide).**
  - **Methylcitrate trio:** O1 + R2075 closed, which is the same chemistry and closure already authorized and
    built on the other line (b24a5e78, uncommitted).
  - **Methylmalonyl five:** leave blocked. There is no mutase gene, three of its reactions have no gene, and the
    route only ever ran with R2112 backwards.
  - **R834:** leave blocked and record why. A carrier that only imports a conserved moiety cannot carry net flux in
    a steady-state model without a matrix CoA demand. Adding one would not make its gene essential, because the
    gene-less acyl-CoA transports bypass it.
  - **R1172:** not proposed for reopening (recorded user decision; no new HMG-CoA transport evidence was sought or
    found).

## 1. Run

- **Model:** `model/candidates/E5_353cefdc_xref_20261115.xml`, SHA-256 `7ccbb568…2146`, checked by every script
  before use. 7ccbb568 and 353cefdc have identical bounds, stoichiometry, gene rules and objective (static check,
  this session), so the recorded certified sweep on 353cefdc applies to it.
- **Conditions:**
  - "open": every exchange at (−1000, 1000), as in the gap census;
  - "Task A": medium `6eca70fd…` + PO1f profile `35307853…` through the platform's
    `load_effective_simulation_context`.
- **Solvers:** cobra 0.30.0. Gurobi 13.0.2 at feasibility and optimality 1e-9, with solution and basis discarded
  before each LP. GLPK with `glp_adv_basis` before each LP. A reaction counts as blocked when |min| and |max| ≤ 1e-9.
- **LPs: 1,147 of the 1,200 approved.**
  - Step 1: 244.
  - Step 2: 92 GLPK + 351 Gurobi, plus 2 × 86 in two Gurobi attempts that stopped on a script error in the
    minimum-flux objective (optlang rejected the sum form, then an integer 0; logs kept as
    `step2_gurobi_attempt*.log`), plus a 1-LP smoke test.
  - Step 5: 205.
  - Step 5b: 54.
  - Step 6: 28.

  pFBA is counted as 2 LPs.
- **Web:** 13 read-only fetches (9 Europe PMC searches, 4 KEGG REST: the two links fetched twice, once to read and once to save). No build, no sweep, no cluster job.

## 2. Reproduction ([tables/step1_reproduce.csv](tables/step1_reproduce.csv))

| Condition | Variant | Nine blocked (Gurobi / GLPK) | Task A WT |
| --- | --- | --- | --- |
| open | as is | 9 / 9 | — |
| open | R1172 reopened | 0 / 0 | — |
| Task A | as is | 9 / 9 | 0.311100535925948 (recorded 684f3752: 0.3111005359259494) |
| Task A | R1172 reopened | 0 / 0 | 0.3120526953810371 (recorded 2f352037: 0.3120526953810365) |

Largest GLPK deviation from Gurobi: about 1e-17. With R1172 reopened, the maxima are:
- open: R25, R95, R552 512.93; R2074–R2076 680; R834 1000;
- Task A: 1.558, 1.559 and 1.912 respectively.

With R1172 reopened, R1172 can only run matrix → cytosol: its range is [−1000, 0] open and [−1.914, 0] Task A.

## 3. Mechanism

### 3.1 The matrix CoA moiety (static; [tables/matrix_moiety_reactions.csv](tables/matrix_moiety_reactions.csv))

- **Moiety.** The pantetheine unit: every species named CoA, coenzyme A or pantetheine
  ([tables/moiety_species.csv](tables/moiety_species.csv)). In the matrix that is CoA, 15 thioesters,
  3'-dephospho-CoA and pantetheine 4'-phosphate.
- **Reactions.** 43 reactions involve matrix moiety species:
  - 34 keep the matrix moiety count unchanged. These include R265 dephospho-CoA kinase and R1292, because
    dephospho-CoA and pantetheine 4'-phosphate are counted.
  - 9 change it, and all 9 are transports (table below).
- **The identity (provable).** Summing S·v = 0 over the moiety species of the matrix gives
  Σ Δᵣ·vᵣ = 0 over the nine transports. Net moiety import into the matrix is zero at every steady state.

| Transport | Equation (moiety per unit forward flux) | Gene | Can it take moiety OUT of the matrix as the model is? |
| --- | --- | --- | --- |
| R834 | CoA cy → mi (+1) | YALI1F20788g | no: irreversible import |
| R1172 | HMG-CoA cy ⇌ mi (+1) | none | closed (0, 0); before, its only usable direction was out (LP) |
| R2025 | propionyl-CoA (m1875) cy ⇌ mi (+1) | none | only the propionyl-CoA that entered on R2113: m78[mi] has no matrix producer (static) |
| R2113 | propionyl-CoA mi ⇌ (m1970) cy (−1) | none | same pass-through (LP: R2113 only imports, R2025 only exports, as is) |
| R2077 | succinyl-CoA mi ⇌ cy (−1) | none | no: cytosolic succinyl-CoA m2046 has no consumer other than R2077 (static) |
| R2083 | 3-methylcrotonyl-CoA mi ⇌ cy (−1) | none | no: cytosolic m2050 has no consumer other than R2083 (static) |
| R2084 | isovaleryl-CoA mi ⇌ cy (−1) | none | only into R2082 → R2083, which comes straight back in (static); blocked (LP) |
| R2118 | isobutyryl-CoA mi ⇌ cy (−1) | none | needs R2117, which takes FAD in the cytosol and returns FADH2 in the matrix; FADH2 can only leave as FADH2 (R847), and cytosolic FADH2 has no outlet but the reverse of the blocked R2142. Blocked (LP) |
| R1291 | pantetheine 4'-phosphate cy ⇌ mi (+1) | none | no: matrix m1350 has no producer, only R1292 (static); blocked (LP) |

So, with R1172 closed, no exit is open for CoA that has become free matrix CoA. R2025 can only pass on
propionyl-CoA unchanged. Producers, consumers and transports of each crossing species are in
[tables/coa_species_reactions.csv](tables/coa_species_reactions.csv).

### 3.2 LP evidence ([tables/moiety_import_bounds.csv](tables/moiety_import_bounds.csv), [tables/routes_moiety_channels.csv](tables/routes_moiety_channels.csv))

- **Net propionyl-CoA import**, P = v(R2025) − v(R2113). By the identity, P equals the net moiety export through all
  other channels.
  - As is: max P = min P = 0, open and Task A, on both solvers.
  - R1172 reopened: max P = 512.93 (open) and 1.558 (Task A), equal to R25's maximum.
  - Maximum HMG-CoA export through R1172: 1000 open, 1.914 Task A. Maximum of R834: 1000 and 1.912.
- **Routes.** Each of the nine is fixed at min(1, half its maximum), with the `xMAINTENANCE` lower bound set to 0
  and the total absolute flux minimised.
  - R1172 reopened:
    - R25, R95 and R552: one propionyl-CoA moiety in (R2025 open, R2113 Task A), one HMG-CoA out (R1172).
    - The methylmalonyl five: succinyl-CoA in on R2077, HMG-CoA out on R1172.
    - R834: CoA in, HMG-CoA out.

    Always 1:1. In Task A, the propionyl-CoA comes from cytosolic isoleucine breakdown (R2285). The methylmalonyl
    route then runs R2112 backwards, carboxylating propionyl-CoA to methylmalonate semialdehyde.
  - CoA drain (the census definition): the moiety enters and is drained in the matrix.
  - FAD drain: the moiety leaves as isobutyryl-CoA on R2118. That opens only because the drain absorbs the
    cytosol-to-matrix FAD flow of R2117 (section 3.1).
- **Before the closure, at the WT optimum** (Task A, R1172 reopened):
  - FVA at 99.9 %: the eight non-transport reactions [0, 0.0091], R834 [0.048, 0.074].
  - pFBA: R834 = R1172 = R83 = R412 = R410 = 0.0721, the other eight 0.
- **Why reopening R2121 frees none.** R2121 is cytosolic and contains no matrix moiety species (static).

### 3.3 Net-zero cycle versus net transfer

- **CoA moiety.** Always net zero across the inner membrane (the identity). R1172 + R834 at the WT optimum is a
  moiety-neutral cycle.
- **Acyl groups.** The cycle is a net transfer. Before the closure, two matrix acetyl-CoA plus the acetyl-CoA of
  R412 left as HMG-CoA. When the nine ran, propionyl or succinyl groups came in and acetyl units went out as
  HMG-CoA. The closure removed the only such exchange, an acyl-group exchange paid for in CoA.

### 3.4 What is proven, what is inferred

| Statement | Kind |
| --- | --- |
| No matrix reaction creates or destroys the CoA moiety; 9 transports move it; net import is 0 at every steady state | provable (static scan plus S·v = 0) |
| The nine are blocked as is and free with R1172 reopened, both conditions | observed: LP, two solvers, all optimal |
| As is, no exit is left for free matrix CoA (max net propionyl-CoA import 0) | observed: LP, two solvers |
| R2077, R2083 and R1291 can only import; R2025 only passes propionyl-CoA through | statically verified (reaction lists of the far-side species) |
| R2118's exit depends on R2117's cytosol-to-matrix FAD flow | inferred from stoichiometry, consistent with the FAD-drain route (LP) |
| Before the closure, R834 only returned the CoA of the HMG-CoA route at the optimum | observed: one pFBA optimum plus FVA at 99.9 % (another optimum could differ) |
| The 1:1 pairing of each route with R1172 | observed in the minimum-total-flux solutions; other routes may exist, but every one must have net moiety 0 (identity) |

## 4. Genes of the nine ([tables/genes_of_the_nine.csv](tables/genes_of_the_nine.csv))

- **Rules.** As recorded in 7ccbb568.
- **Names.** No name or function below is verified here. The short functions are the model's reaction names, plus
  the experimental list's annotation text where given.
- **KO/WT and class.** From the certified sweep on 353cefdc (`sweep684_candidate.tsv`, `49b4f27b…`; Task A +
  PO1f, Gurobi 1e-9), which is LP-identical to 7ccbb568. Scoring is rule B, with mtDNA and untargeted genes out of
  scope. Negatives are the user-defined class (model genes not on the consensus list). The list was used in
  development, so this is a development reference.
- **W29 audit.** Every gene is `ok` in `w29_tag_protein_audit.csv` (not needed: no sequence work).

| Reaction (compartment) | Rule | On the consensus list | KO/WT, class at 15 % |
| --- | --- | --- | --- |
| R25 2-methylcitrate synthase (mi) | YALI1E01019g or YALI1E03300g (both also on citrate synthase R236/R237) | no, no | 1.0 TN, 1.0 TN |
| R95 "aconitase" (2-methylcitrate → 2-methylisocitrate, mi) | YALI1F03803g | **yes** (list annotation: similar to S. cerevisiae PDH1) | 1.0 **FN** |
| R552 methylisocitrate lyase (mi) | YALI1C24124g or YALI1F39620g | no; **yes** (list annotation: similar to ICL2) | 1.0 TN; 1.0 **FN** |
| R2074 methylmalonate-semialdehyde dehydrogenase (cy) | YALI1B02194g or YALI1C04093g or YALI1E00588g or YALI1F31215g | no ×4 | 1.0 TN ×4 |
| R2075 methylmalonyl-CoA mutase (cy) | none | — | — |
| R2076 methylmalonate + CoA → methylmalonyl-CoA, no ATP (cy) | none | — | — |
| R2077 succinyl-CoA transport (mi–cy) | none | — | — |
| R2112 methylmalonate-semialdehyde dehydrogenase, acylating (cy) | YALI1C02550g (also R2119, R2121) | no | 0.998 TN |
| R834 CoA transport (cy → mi) | YALI1F20788g | **yes** (list annotation: similar to S. cerevisiae LEU5) | 1.0 **FN** |

A model gene rule is not experimental validation. The three FN genes are FN because no reaction they carry is needed
for growth on Task A. The R834 gene would stay non-essential even with a matrix CoA demand (section 6).

## 5. Options and evidence

Tiers:
- **verified:** read this session, abstract or record. Abstracts were read through the Europe PMC REST search
  (`webservices/rest/search`; results saved in the artifacts folder `lit/`); the europepmc.org article links below are built from those
  PMIDs and were not opened as pages;
- **second-hand:** read by the methylcitrate round on the other line (`gene-essentiality-reactions-d09601`,
  `model/reports/methylcitrate_cycle_20261116/REPORT.md` §3), not re-opened here;
- **established / uncited:** textbook knowledge, no source opened;
- **inferred:** reasoned.

### O0. Leave the nine blocked, and record why

- **For.**
  - They were only ever active through a gene-less net exit of matrix CoA (R1172) and, depending on the reaction,
    the gene-less acyl-CoA transports R2025, R2113 and R2077 (LP, section 3).
  - At the WT optimum before the closure, only R834 was used, and only as the CoA return of R1172.
  - The methylmalonyl route has no mutase gene in *Y. lipolytica*. Verified: KEGG returns no gene for EC 5.4.99.2
    ([yli link](https://rest.kegg.jp/link/yli/ec:5.4.99.2), empty; human control
    [hsa link](https://rest.kegg.jp/link/hsa/ec:5.4.99.2) returns `hsa:4594`).
  - The mutase route ran only with R2112 reversed: a reductive carboxylation, opposite to the oxidative
    decarboxylation the enzyme name describes (LP observation; the thermodynamic direction is
    established / uncited).
- **Against.** AGENTS.md §5 treats blocked reactions on essential genes as gaps to fix. Three of the nine carry
  consensus-essential genes (R95, R552, R834). The methylcitrate cycle is real biology:
  - *S. cerevisiae* ICL2 is a matrix 2-methylisocitrate lyase. Verified (abstract): Luttik et al. 2000,
    [PMID 11092862](https://europepmc.org/article/MED/11092862).
  - A *Y. lipolytica* phd1 deletion cannot grow on propionate (Park et al. 2018,
    [PMC5991449](https://europepmc.org/article/PMC/PMC5991449); second-hand).

### O1. Matrix propionyl-CoA producer: 2-oxobutanoate dehydrogenase (with R2075 closed)

- **Reaction (in memory).** `m42[mi] + CoA[mi] + NAD[mi] → propionyl-CoA[mi] + CO2 + NADH + 4 H+`, mass and charge
  balanced (`check_mass_balance` = {} this session). Rule = R693's PDH rule, INFERRED. O1b uses the
  branched-chain 2-oxo acid dehydrogenase rule of R1969 instead.
- **For.**
  - It closes the CoA loop inside the matrix: propionyl-CoA is made from matrix CoA and R25 returns it. No
    transport is needed.
  - Threonine-derived 2-oxobutanoate feeds odd-chain fatty acids in *Y. lipolytica* (Park et al. 2019,
    [PMC6987463](https://europepmc.org/article/PMC/PMC6987463); second-hand).
  - Rat mitochondrial PDH and BCKDH both decarboxylate 2-oxobutyrate (Paxton et al. 1986,
    [PMID 3718468](https://pubmed.ncbi.nlm.nih.gov/3718468/); second-hand).
  - The user already authorized this exact chemistry with "Pyruvate dehydrogenase, but need marked as "just infer""
    and "Close R2075" on the other line, which built candidate b24a5e78 on cbf2ef95 (uncommitted).
- **Against.**
  - No *Y. lipolytica* evidence for the catalyst or its compartment (second-hand).
  - Alone, O1 opens side route A: R2113 out → R2112 backwards → R2074 → R2076 → R2075 → R2077 in (LP, section 6).
    Hence the R2075 closure.
  - R834 becomes active only through side route B (CoA in on R834, propionyl-CoA out on the gene-less R2025, then
    R2024 backwards and propanoate export). That is not a real CoA exit (inferred: acyl-CoAs are not known to cross
    the inner membrane without the carnitine system; established / uncited).

### O2. Valine, isoleucine and methylmalonate steps in the matrix (not tested)

- **What the model shows** (static, verified):
  - The matrix branched-chain dehydrogenases stop early: isoleucine at 2-methylbutanoyl-CoA (m1830[mi], a dead
    end); valine is exported as isobutyryl-CoA on R2118.
  - The rest of the valine chain (R2114–R2117) and the methylmalonyl steps sit in the cytosol, and R2117 mixes
    compartments (FAD cy → FADH2 mi).
  - Placing these steps in the matrix would make propionyl-CoA and succinyl-CoA from matrix CoA, which is
    moiety-neutral.
- **Against, and why not tested.**
  - No evidence was found for their compartment in *Y. lipolytica* (not searched beyond the abstracts above).
  - The mutase has no gene (O0).
  - It would mean several new matrix reactions. That is a curation project of its own, beyond this round's LP
    budget.
- **Recorded as a follow-up gap:** R2117's mixed compartments.

### O3. Propionyl-carnitine shuttle

- **Reactions (in memory).** Three new reactions: cytosolic carnitine propionyltransferase (rule of R191), a
  carrier (rule of R825) and a matrix carnitine propionyltransferase (rule of R193). New species O-propanoylcarnitine
  C10H19NO4. All rules INFERRED.
- **For.**
  - *S. cerevisiae* Crc1 transports carnitine, acetylcarnitine and propionylcarnitine in proteoliposomes. Verified
    (abstract): Palmieri et al. 1999, [PMID 10622748](https://europepmc.org/article/MED/10622748).
  - Propionylcarnitine is measured in *Y. lipolytica* and correlates with CAT2 expression under 1-propanol.
    Verified (abstract; correlation only): [PMID 41301536](https://europepmc.org/article/MED/41301536).
  - Carnitine acyltransferases are up-regulated on propionate. Verified (abstract; transcriptomics, the authors'
    interpretation): [PMID 41214791](https://europepmc.org/article/MED/41214791).
- **Against.**
  - No *Y. lipolytica* enzyme shown to transfer propionyl to carnitine.
  - Three new reactions instead of one.
  - In the model it opens a futile loop (R834 up to 68.7: CoA in, propionyl-CoA out on the gene-less R2113, back in
    as propionylcarnitine). The loop moves 4 H+ into the matrix per turn. It is dissipative, but it is an artifact
    of the gene-less R2113.

### O4. Matrix CoA turnover

- **Evidence.**
  - Δleu5 *S. cerevisiae* has 15-fold less mitochondrial CoA, with cytosolic CoA unchanged. Verified (abstract):
    Prohl et al. 2001, [PMID 11158296](https://europepmc.org/article/MED/11158296).
  - Human SLC25A42 imports CoA only by counter-exchange (for matrix adenine nucleotides and PAP); CoA is made
    outside the matrix. Verified (abstract): Fiermonte et al. 2009,
    [PMID 19429682](https://europepmc.org/article/MED/19429682).

  So matrix CoA is imported to fill a growing matrix pool. Steady-state FBA sees that only as a growth-dilution
  demand.
- **Counter-exchange alone does not help.** Rewriting R834 as a counter-exchange would not create an exit, because
  the moiety identity still holds (static).
- **Placeholder test.** 1.9e-5 matrix CoA in biomass_C (10 % of the cytosolic CoA coefficient; no evidence for the
  value) frees all nine, but only up to 5.9e-6 (= coefficient × growth). Knocking out R834's gene changes growth by
  less than 1e-6, because the demand is met by propionyl-CoA imported on the gene-less R2113 (LP).
- **Against.** No measured matrix CoA content for *Y. lipolytica* was found. AGENTS.md §6 forbids adding biomass
  demand merely to improve recall, and the test shows it would not change the FN anyway.

### Excluded

- **Reopening R1172.** A recorded user decision (2026-11-12). No evidence for HMG-CoA crossing the inner membrane
  was sought or found in this round.
- **Making R834 reversible.** The verified evidence is about import (Prohl 2001, Fiermonte 2009). Making it
  reversible would only recreate an R1172-like gene-assigned exit.
- **Any added exchange, sink or gene-less transport.**

## 6. In-memory tests ([tables/options_summary.csv](tables/options_summary.csv), [tables/targeted_knockouts.csv](tables/targeted_knockouts.csv))

Task A + PO1f, Gurobi 1e-9, R1172 closed in every variant. Energy check: every boundary lower bound and the
`xMAINTENANCE` lower bound set to 0, then maximise `R_NTP1` and `biomass_C`.

| Variant | WT | Max of R25 / R2074 / R834 | Energy check (ATP, biomass) | Targeted KO class changes at 15 % (26 genes) |
| --- | --- | --- | --- | --- |
| B (as is) | 0.311100535925948 | 0 / 0 / 0 | 0, 0 (optimal) | — |
| O1 | 0.31110053592595 | 7.599 / 7.601 / 3.311 | 0, 0 | none |
| O1b (BCKDH rule) | same LP as O1 | same | same | none (8 PDH/BCKDH genes) |
| O1 + R2075 closed | 0.31110053592594944 | 7.599 / 0 / 3.311 | 0, 0 | not run (same as O1 except the gene-less R2075) |
| O3 | 0.3111005359259493 | 5.801 / 5.803 / 68.68 | 0, 0 | none |
| O4 (placeholder) | 0.31109829973373154 | 5.9e-6 / 5.9e-6 / 5.9e-6 | 0, 0 | none |

- **Cross-check.** GLPK reproduces O1 + R2075 closed: R25 7.5987, methylmalonyl five 0, R834 3.3108, all optimal.
- **Targeted genes:**
  - the genes of the nine;
  - the PDH, BCKDH and carnitine-shuttle rules;
  - threonine deaminase, R2285, R265 and R2024.
- **What the knockouts do not prove.** They cannot show that no other gene changes class. Only a full certified
  sweep can. The other line's full sweep with the same O1 + R2075 package on cbf2ef95 found no class change
  (second-hand).

## 7. Figures ([figure/](figure/))

Hand-laid schematics. Each box shows the reaction with its gene rule and flux: green = gene rule, pink = no gene,
dashed dark red with a cross = closed, black outline = new in the option.

- [r1172_status.png](figure/r1172_status.png): largest flux each reaction can carry, before / after / option.
- [r1172_before_route.png](figure/r1172_before_route.png): before. The methylcitrate cycle forced: its CoA enters
  as propionyl-CoA (R2113) and leaves as HMG-CoA (R1172).
- [r1172_before_wt.png](figure/r1172_before_wt.png): before, at the growth optimum: R834 = R1172 = 0.072.
- [r1172_after_wt.png](figure/r1172_after_wt.png): after (7ccbb568), at the growth optimum: the nine at 0.
- [r1172_option_route.png](figure/r1172_option_route.png): O1 + R2075 closed, the methylcitrate cycle forced:
  propionyl-CoA made inside the matrix, no CoA crossing.

## 8. What it means, and could the method be better

- **The closure was correct and its side effect is real, but not a loss of biology.** The nine depended on a
  gene-less exit of matrix CoA, not on R1172's HMG-CoA chemistry. The right fix gives the matrix its own
  propionyl-CoA source (O1). Re-creating an exit would be the wrong fix.
- **The deeper gap is the gene-less acyl-CoA transports** (R2025, R2113, R2077, R2083, R2084, R2118) and the
  gene-less pantetheine 4'-phosphate transport R1291. They let CoA moiety enter the matrix in forms that need an
  exit. Under O1, O3 or O4 they make R834 look active or bypassed. This is a follow-up review, not done here.
- **A better method.** A conserved-moiety (left null space) scan per compartment would find every pool like this
  one, together with its exits, in one pass. Every reaction whose activity depends on a single exit would show up
  before a closure is built.
  - How to test it: compute the conserved moieties of the stoichiometric matrix restricted to each compartment.
    List each moiety's boundary reactions and which ones can carry it out (one FVA per boundary reaction). Run it
    on 684f3752's parent to see whether it would have predicted these nine.
  - The census's single-closure re-test missed them because it re-tests a closure only inside its own module.

## 9. Limits

- **One medium**, with the open condition only for reproduction and mechanism.
- **Minimum-total-flux routes and pFBA give one solution each.** The identity, not the route, is the proof.
- **No full knockout sweep** for any option. Targeted knockouts only (26 genes).
- **O4's coefficient is a placeholder.**
- **O2 not tested.**
- **Literature:**
  - abstracts only, no full texts opened;
  - *S. cerevisiae* and human evidence is related-organism context;
  - the *Y. lipolytica* items (Park 2018/2019) are second-hand.
- **Gene names:** none verified.
- **Not addressed:**
  - the m1875/m1970 duplicate cytosolic propionyl-CoA;
  - R2117's compartment mix;
  - R2112's direction;
  - the gene-less acyl-CoA transports.

## 10. Files and commands

- **Interpreter.** `/opt/miniconda3/bin/python3` with the real HOME (Gurobi licence), run from [scripts/](scripts/).
  Every script checks the model, medium and profile SHA-256 first.
- **Commands, in order:**
  1. `step1_reproduce.py gurobi` and `step1_reproduce.py glpk`
  2. `step2_static.py`
  3. `step2_lp.py glpk` and `step2_lp.py gurobi --full`
  4. `step3_genes.py`
  5. `step5_options.py`
  6. `step5b_routes.py`
  7. `step6_figure_flux.py`
  8. `draw_schematic.py ../figure`
  9. `make_tables.py`
- **Outputs** (`out/step*.json`) and logs: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r1172_side_effect_20261116/`
  (`MANIFEST.sha256`).
