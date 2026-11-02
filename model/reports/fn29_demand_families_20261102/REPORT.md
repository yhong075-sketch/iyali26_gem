# The 29 "not solvable by curation now" false negatives: sphingolipid, GPI, cardiolipin, N-glycan and mitochondrial-translation demands

> **Filed 2026-11-02** from the gap session's scratchpad as a documentation-only task ([TASK.md](TASK.md)). The text
> below is the synthesis report as delivered. Only these edits were made:
> 1. `$FN29` now points to the archive;
> 2. one source-file path was corrected (the mitochondrial `verify_lp/` files);
> 3. two sections were added, "Status on the current main" and "Reader notes" (at the end).
>
> No number or verdict was changed. The delivered file is archived unchanged as
> `$IYALI26_RESEARCH_ROOT/artifacts/tasks/fn29_demand_families_20261102/REPORT_fn29_demand_families_as_delivered.md`.
>
> **Companion files:**
> - [tables/genes_29.csv](tables/genes_29.csv) (one row per gene) and [tables/families_5.csv](tables/families_5.csv)
>   (one row per family), parsed from the archived JSON by [scripts/make_tables.py](scripts/make_tables.py).
> - The Escher status figure, [escher/](escher/README.md), with each family's recommended variant as an in-memory
>   what-if.
> - The archive `$IYALI26_RESEARCH_ROOT/artifacts/tasks/fn29_demand_families_20261102/`: scripts, JSON, logs, `workflow_journal.jsonl` (every agent's full return value) and
>   `MANIFEST.sha256`.

Consolidated from five scout reports and fifteen adversarial verdicts (three lenses per family: LP reproduction, biology and evidence, project rules and tuning). Where a verdict refutes or corrects a scout claim, the verdict wins unless its own evidence is weaker; every such case is listed per family under "Refuted or corrected".

- Model executed by every scout and verifier: `e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897` (origin/main candidate `..._aco_split_20261029.xml`, commit c8813aea), scratch copy `$FN29/e1e7c868.xml`. Condition: `$FN29/cond_taskA_po1f.json` (Task A medium, R612 closed, xMAINTENANCE lb 7.8625). Solver: Gurobi 13.0.1/13.0.2 via cobra 0.30.0, tolerance 1e-9, real HOME licence. Every verifier noted that this candidate is superseded on origin/main (164542fb, c1f70ac4 2026-10-31, 2bc263f2 2026-11-01, kept 2026-11-02); no reference model is designated, so the base for any build is the user's to name.
- Experimental list: `consensus_essential_genes.csv` SHA 1e887f5a…, 1,612 genes, 326 in the model. Model genes not on the list count as the user-defined negative class (TP/FN/FP/TN reported); positive coverage is 326 of 1,612.
- `$FN29` below stands for the archived scratchpad folder `$IYALI26_RESEARCH_ROOT/artifacts/tasks/fn29_demand_families_20261102/fn29/` (in the delivered report it stood for the gap session's scratchpad `/private/tmp/claude-501/…/scratchpad/fn29`).
- Nothing was written to any model file or to the repository; every variant is in memory.

## Status on the current main (added 2026-11-02)

**Every count in this report comes from `e1e7c868` under the pre-2026-10-30 scoring convention.** Some use a mixture
of conventions. They must be re-run on a base the user names before any TASK.md.

- **Newer candidates.** `origin/main` at `80cfc4f` carries these candidates after `e1e7c868` (read from
  `model/candidates/README.md` and `model/STATE.md`, 2026-11-02):
  - `f505571d` (2026-10-30): R24 option C on `e1e7c868`.
  - `d2e4c0c4` (2026-10-30): the first trehalose build. Then `164542fb` (2026-10-30): trehalose removed from
    `biomass_C`.
  - `c1f70ac4` (2026-10-31): `164542fb` + R24.
  - `2bc263f2` (2026-11-01, kept 2026-11-02): `c1f70ac4` plus a ceramide pool in `biomass_C` at 5.4e-4 mmol/gDW.
  - No reference model is designated.
- **Scoring.** Since 2026-10-30, `model/expected/benchmark.md`'s primary rule scores certified maintenance-shortfall
  infeasible knockouts as growth 0, and puts the mtDNA genes and COX9 out of scope (rule B).
  - The rule-B baseline this report uses, 136/190/46/692, is the one STATE.md's 2026-10-30 entry ran on `e1e7c868`
    (Task A + PO1f, 1,123 LPs).
  - STATE.md records later rule-B baselines of 136/190/43/695 on `c1f70ac4` and 142/184/44/694 on `2bc263f2` (recalled
    from STATE.md, not re-run here).
- **Sphingolipid on `2bc263f2` (recalled from the STATE.md entry of 2026-11-01, not re-run here).** The ceramide pool
  already turned six of the ten sphingolipid genes from FN to TP: LCB1 `YALI1C04265g`, LCB2 `YALI1F20637g`, TSC10
  `YALI1B22883g`, ELO `YALI1B26350g`, PHS1-like `YALI1F15781g` and TSC13-like `YALI1A05116g`. It also made the KCR
  ortholog `YALI1A06453g` a FP.
  - That entry names SUR2 and SCS7 as still FN.
  - AUR1 `YALI1E19647g` and `YALI1D08000g` are not among its flips, so they are inferred to stay FN.
  - Section 3.1's V1/V2 numbers therefore describe a model without that pool.
- **Bucket-4 session.** The "Disconnected-pathway essentiality FNs" session reconciled its draft with this report on
  `2bc263f2`. The draft is uncommitted at
  `.claude/worktrees/dreamy-moore-d727a1/model/reports/bucket4_disconnected_20261101/TASK.md`. Its handover note
  reports (not checked here):
  - mitochondrial translation, N-glycan to Man5, VRG4 and MET7 behave identically on `2bc263f2` (296 targeted LPs);
  - on `2bc263f2`, AUR1 needs the lipid owner to re-partition the ceramide pool;
  - its recommended order differs from section 4.

## 1. Answer

Of the 29, 28 can be turned into true positives mechanically (every gene except YALI1D08000g, whose GPR is mis-assigned), but only 18 of them under the project's own rules, and none of the 18 without a named decision. The 18 are: sphingolipid 9 of 10 (variant V2, a yeast-7.6-sourced placeholder term plus three structural edits; or 8 of 10 with the no-reversal V1), GPI 6 of 6 (variant b1/b2, an S. cerevisiae-sourced proxy term whose scale conflicts with the model's lipid pool), cardiolipin 2 of 2 (a Y. lipolytica-sourced term that the colleague-owned lipid-unlump ledger has already recorded as "blocked"), and N-glycan 1 of 3 (the VRG4-type transporter, via an evidence-backed R289 stoichiometry fix plus a Golgi route and a gene-less closure that are design choices). The remaining 10 flippable genes (ALG7, ALG1 and the 8 mitochondrial synthetases) are blocked on evidence, not on a decision: their fixes need a biomass coefficient for which no Y. lipolytica, yeast or Yeast-GEM source exists, which the verifiers classed as demand-to-flip tuning until a source is found or the user names an unsourced structural term as an explicit exception; the mitochondrial fix additionally needs mechanism evidence for two gene-less transport reversals (R1251 Met, R1218 Ile). The cost of the 18 is five new false positives (KCR YALI1A06453g, PBN1 YALI1E20666g, GEP4-like YALI1A13994g, GDA1-like YALI1C25093g, CWH8-like YALI1D26639g; four if the second GDA1-like gene YALI1E39570g is added to R602), two side-effect true positives outside the 29 (YALI1F08225g, YALI1A09121g), WT growth losses between 0.04 % and 11.3 %, no energy-generating cycle in any recommended variant (ATP-hydrolysis maximum 0.0, optimal, no uptake), the reopening of two recorded entries (the 2026-10-21 GPI exit decision only if R_DM_GPI is closed; the 2026-10-26 census out-of-scope classification for mitochondrial translation), the opening of ten inherited Golgi closures with no recorded rationale, and three edits that fall inside the colleague-owned lipid-composition territory (every sphingolipid, GPI and cardiolipin term changes biomass_C lipid composition). All sweep counts were produced on e1e7c868 under the pre-2026-10-30 scoring convention (or a mixture of conventions); the project-comparable rule-B counts (136/190/46/692 baseline) were inferred, not run, by the verifiers, and must be re-run on the candidate the user names before any TASK.md.

## 2. Family table

Counts at 15 % KO/WT. "Rule A" = earlier rule (infeasible knockouts unresolved, counted as neither class); "mt-excl" = 16 mtDNA genes excluded; "rule B" = benchmark.md primary since 2026-10-30 (certified-infeasible scored 0, mtDNA + COX9 out of scope; baseline 136/190/46/692). Rule-B after-counts are the project-rules verifiers' inferences (not run) unless stated.

| Family | Genes (n) | Root cause | Fix type | Best variant | Coefficient (mmol/gDW) | Source status | FN→TP at 15 % | New FP | Other class changes / counts | EGC (no uptake, NGAM 0) | Decisions reversed / reopened | Owner | Verdict | Confidence after verification |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Sphingolipid | 10 | Nothing demands any sphingolipid; ER chain from R722 blocked; Golgi AUR1 set R1825-R1834 inherited (0,0) and dead-ended (duplicate IPC species m1707-m1716, DAG m1704[gm] no exit); YALI1D08000g GPR mis-assigned | Composition value + pathway connection + opening inherited closures | V2 (IPC-C C26 term + open R1825-R1834 + 10 duplicate-species links + DAG gm→em link); fallback V1 (ceramide-3 term only) | 0.000417 (V2 and V1) | Sourced yeast: yeast 7.6 r_2108 `s_0535` complex sphingolipid, copied term by term into legacy R1371 (verified); no Y. lipolytica value; species chosen because it flips AUR1 (flagged) | 9 of 10 (V2); 8 of 10 (V1, no full sweep) | 1 (YALI1A06453g KCR, sole gene of R1413-R1419) | none; 125/190/46/700 → 134/181/47/699 rule A, 20 unresolved both sides (reproduced); rule B inferred 136/190/46/692 → 145/181/47/691 | optimal 0.0 (both) | No recorded decision reversed (verifier correction); opens 10 inherited closures with unknown rationale; adds a gene-less DAG link in the reverse direction of inherited R1795; merges species listed in the lipid line's unlump_stage0_plan.csv | Lipid-unlump line (colleague) for species, coefficient and the duplicate pools; user for the closures and the DAG link | Solvable after decision (9); YALI1D08000g not by curation (needs GPR correction + GlcCer branch) | High: all V2 numbers reproduced exactly (1,433 LPs); corrections minor |
| GPI anchor | 6 | Chain connected and feasible since 2026-10-21 but m603[C_er] has only R376 and the unforced demand R_DM_GPI; biomass_C has no GPI term | Composition value (+ optional closure of R_DM_GPI) | b1 (term, exit open, no reversal) or b2 (term, R_DM_GPI (0,0)); identical optimum | 1.6962495847e-4 | Sourced yeast: Yin 2007 S. cerevisiae, 1.9e6 CWP/cell / 18.6 pg/cell (verified via Europe PMC); scale conflicts with the model's lipid pool (PI term 7.5e-5) — major correction; no Y. lipolytica value | 6 of 6 | 1 (YALI1E20666g PBN1, AND partner on R372; screen scores borderline, essential in S. cerevisiae) | none; 125/190/43/691 → 131/184/44/690 mt-excl, 16 unresolved (reproduced; FP 43 not 42 because the condition file omits the PO1f LEU2 plasmid complement) | optimal 0.0 (b1, b2); variant a infeasible (rejected) | 2026-10-21 "keep the GPI exit open" reversed only by b2; b1 "reopens the question" the decision anticipated; gpi_chain_links §3 "needs a measured amount" | Not lipid-owned (curation/biomass precedent) but consumes the lumped PI/PE/DAG species the unlump line will replace; consult colleague | Solvable after decision | High on LPs (1,316 LPs, exact); major open point on coefficient scale |
| Cardiolipin | 2 | biomass_C has no cardiolipin/PG term (cardiolipin was priced only in the legacy xLIPID→xBIOMASS, not the objective); route R1791→R197→R1836→R1722→R1721 is linear and unused | Composition value only | a (m1651[C_mm] term) | 7.1616e-4 (corrected 7.152e-4 for CL 72:7; range 6.27e-4 to 8.06e-4); placeholder-consistent alt (ii) 2.12e-5 | Sourced Y. lipolytica: Sekova 2019 W29 12.98 mg membrane lipid/gDW × 7-9 % CL (stationary, glycerol); conflicts 5-6× with Hein & Hayen 2012 (CL 2.8 % of GPs, glucose); owner ledger records it as "not source supported" | 2 of 2 | 1 (YALI1A13994g GEP4-like, sole PF09419 PGP phosphatase in CLIB122) | +1 TP outside the 29 (YALI1F08225g mitochondrial CTP/CMP carrier, on list); 125/190/43/691 → 128/187/44/690 mt-excl, 16 unresolved (reproduced); rule B inferred 136/190/46/692 → 139/187/47/691 | optimal 0.0 | No STATE.md decision; but `lipid_unlump_sn_core_curation.json` records `cardiolipin_biomass_weight_decision: blocked`, biomass_C `candidate_action: no_change` (verifier correction, major) | Lipid-unlump line (colleague) | Solvable after decision | High on LPs (1,342 LPs, exact); major ownership conflict |
| N-glycan | 3 | Dolichol-linked N-glycan core split into three unlinked compartment copies (UM18) with no OST and no N-glycan demand; mannan term made by the PMT lump R289 that creates C18H33O16 per turnover; gene-less R1137 bypasses R850; Golgi Pi has no exit | Pathway connection + composition value + stoichiometry correction + new closure | C2 = A2 + B (recommended by scout); verifiers split it: B (VRG4) decidable, A2 (ALG7/ALG1) blocked | N-glycan 0.01 (placeholder); mannan 0.235849 unchanged | No source for the N-glycan term (Yeast-GEM has none; Park 2011, De Pourcq 2012, Orlean 2012 give none; Vega & Dominguez 1986 unread); no source for the 1:3 Dol-P-Man:GDP-Man split | 3 of 3 (C2) mechanically; 1 of 3 (B alone) under rules | 2 (YALI1C25093g GDA1-like, YALI1D26639g CWH8-like); 1 if R602/R586 gets the second GDA1-like gene YALI1E39570g (orthology level) | none; 136/190/55/700 → 139/187/57/698 (20 infeasible as 0, mtDNA+COX9 in scope; reproduced); rule B inferred → 139/187/48/690; 154-159 KO/WT ratios move > 1e-6 (2026-10-28 gate) | optimal 0.0; dolichol leak R287 [0,0] | None reversed; builds on 2026-10-28 dolichol entry; R1137 new closure moved from "evidence in hand" to decision list (verifier) | Not lipid-owned; UM18 owner unknown in STATE.md | Split: VRG4 YALI1F28833g solvable after decision; ALG7, ALG1 not by curation now (blocked on composition evidence) | High on LPs (837 LPs); refuted on the coefficient by two lenses |
| Mitochondrial translation | 8 | 14 mitochondrial tRNA pairs are RNP/RNC dead ends; no translation pseudo-reaction, no mt-protein term; no mitochondrial Met import (R1251 mi→cy only); no cytosol→mitochondrion Ile import (R1218) | Pathway connection + composition value + two transport-direction changes + reversal of a census classification | a2 + R1251 reversible + R1218 reversible (scout never ran the combination; the LP verifier did, full sweep) | 0.01863252 (= 3.726504 residue sum × 0.005); alt 1e-5 | No source (fraction of cell protein that is mt-encoded); residue fractions inside the pseudo-reaction sourced from AJ307410 (verified) | 8 of 8 mechanically | 3 (YALI1C15180g Ile aaRS, no UniProt location; YALI1C25598g Trp aaRS; YALI1F33951g mitochondrial adenylate kinase) | +1 TP outside the 29 (YALI1A09121g Thr aaRS); 125/190/46/700 → 134/181/49/697 rule A, 20 unresolved (verifier's own full sweep); ILV2/5/3-like stay 0.224 | optimal 0.0 (verified incl. R1218) | 2026-10-26 census classification "out of scope" (a computation-only entry, not a user decision — verifier correction); 2026-10-29 listed it as follow-up category d | Not lipid-owned; this worktree's line | Not by curation now (blocked on coefficient source and §5 transport evidence); decidable only as a named exception | High on LPs (1,176 LPs); refuted on coefficient and transports |

Totals if the 18 decidable genes were built (families tested separately; the combined model was not run and interactions are unknown): FN→TP 18 (sphingolipid 9, GPI 6, cardiolipin 2, N-glycan 1), new FP 5 (or 4), side-effect TP 1 (YALI1F08225g). If the 10 evidence-blocked genes were also built as exceptions: FN→TP 28, new FP 8 (or 7), side-effect TP 2.

## 3. Per family

### 3.1 Sphingolipid synthesis (10 FN)

Genes: YALI1C04265g (LCB1-like), YALI1F20637g (LCB2-like; AND rule on R722), YALI1B22883g (TSC10-like, R37), YALI1A05116g (TSC13-like, R1427-R1433), YALI1B26350g (ELO2-like, R1406-R1412), YALI1F15781g (PHS1-like, R1420-R1426; R1508-R1510 OR with YALI1E18441g), YALI1D08000g (annotated EC 1.14.19.17 but assigned to the cytosolic IPC synthase set R2143 etc.), YALI1E19647g (AUR1-like, R1825-R1834 (0,0) and R1815-R1824 (0,0)), YALI1E20626g (SUR2-like, R664/R198/R200/R2136), YALI1F17082g (SCS7-like, R199/R201/R206-R209). All identities are homology-level (KEGG yli / live UniProt Q6… entries PE3/PE4); 9 of the 11 model UniProt xrefs are obsolete 2017_01 TrEMBL entries. No Y. lipolytica experimental evidence for any of them.

**Mechanism (verified).** biomass_C carries no sphingolipid, ceramide, IPC, MIPC or M(IP)2C term (its lipid terms: zymosterol 0.003029, ergosterol 0.035251, PS 3.5e-5, TAG 2.34e-4, PI 7.5e-5, PE 2.37e-4, PC 3.46e-4, PA 5.8e-5 — the last omitted by the scout). The legacy R1371 lipid pseudoreaction that demanded 0.000417 "complex sphingolipid" is (0,0) and m1189/m1190 have no producer. So every ER-chain reaction has FVA 0/0 at WT and the knockouts remove nothing. The Golgi AUR1 set R1825-R1834 is closed and dead-ended twice (duplicate IPC products m1707-m1716[go] with no formula; DAG m1704[gm] with 23 reactions all closed). The elongation genes (TSC13/ELO/PHS1-like) are unused only because no C20-C26 acyl-CoA is demanded (probe V0: a hexacosanoyl-CoA demand alone flips all three plus the KCR FP).

**Recommended design (V2, in memory; coefficients exact).**

| Edit | Reaction / term | Bounds |
|---|---|---|
| (a) open 10 inherited closures | R1825-R1834: Golgi ceramide + PI m1679[C_gm] → DAG m1704[C_gm] + IPC m1707…m1716[C_go] (GPR YALI1E19647g) | (0,0) → (0,1000) |
| (b) 10 duplicate-species links [C_go] | m1707→m1265, m1708→m1270, m1709→m1266, m1710→m1271, m1711→m1263, m1712→m1268, m1713→m1264, m1714→m1269, m1715→m1267, m1716→m1272 | (0,1000), gene-less |
| (c) DAG exit | L_DAG_gm_em: m1704[C_gm] → m1649[C_em] | (0,1000), gene-less |
| (d) biomass term | biomass_C += −0.000417 IPC-C (C26) m1271[C_go] (C50H100NO13P, 954.3 g/mol) | — |

V1 fallback: biomass_C += −c ceramide-3 (C26) m398[C_er] (c = 1e-5…1e-2 series, 4.17e-4 nominal); no reaction opened. V1b: C24 species m397. V3 (MIPC via the cytosolic 0.01-stoichiometry set) is not recommended (adds two MIPC-synthase FPs, bypasses AUR1, cuts WT 6.3 %).

**Verified numbers.**

| Quantity | Baseline | V2 | V1 (4.17e-4) | Reproduced by |
|---|---|---|---|---|
| WT growth | 0.3116608216336117 | 0.3114600902481364 (−0.0644 %) | 0.311462 | LP lens (V2 exact); V1 scout only |
| Route flux at optimum | 0 | Σ R1825-R1834 = L_DAG_gm_em = 1.298789e-4 = 0.000417 × μ | — | LP lens |
| KO/WT family (9 genes) | 1.0 | 0.0 (optimal, not infeasible) | 0.0 for 8 (AUR1 stays 1.0) | LP lens (V2), biology lens spot (12 LPs) |
| KO/WT YALI1D08000g | 1.0 | 1.0 | 1.0 | LP lens |
| KO/WT YALI1A06453g (KCR) | 1.0 | 0.0 (new FP) | 0.0 | LP lens |
| KO/WT YALI1F10633g, YALI1F21567g, YALI1E18441g, YALI1F09542g, YALI1D05038g | ~1.0 | ~1.0 (0.99990 for E18441g) | — | LP lens |
| EGC: max ATP hydrolysis, no uptake, NGAM 0 | optimal 0.0 | optimal 0.0 (also biomass 0.0 and each of the 11 new reactions 0.0) | optimal 0.0 (scout) | LP lens |
| Sweep 15 % (rule A, 1,081 genes) | 125/190/46/700, 20 unresolved | 134/181/47/699, 20 unresolved (same genes) | not swept | LP lens (exact at 1/5/10/15 %) |
| Sweep 1 % / 5 % / 10 % | 97/218/36/710; 104/211/38/708; 108/207/43/703 | 106/209/37/709; 113/202/39/707; 117/198/44/702 | — | LP lens |
| Class changes | — | exactly 10 (9 FN→TP, 1 TN→FP); no other KO/WT moves > 1e-3 | 8 FN→TP (22-gene panel) | LP lens |
| Sensitivity (V1 series c = 1e-5…1e-2) | — | — | same 8 flips at every c; WT 0.311656 → 0.306962 | scout only |

**Refuted or corrected by the verifiers.**

| Scout claim | Verdict | Evidence |
|---|---|---|
| Full sweep covered 1,082 genes | Corrected: 1,081 | LP lens: `len(model.genes) == 1081`; scout's own `06_sweep_V2.json` says 1081 |
| L_DAG_gm_em "follows the precedent of R_DAGt_er_em … recorded as a user decision"; inherited R1795 "closed in the wrong direction" | Refuted | Project-rules lens: the 2026-10-21 decision covers R_DM_GPI only; `gpi_chain_links.json` defines R_DAGt_er_em as a reversible ER↔ER-membrane link, not an exit; nothing covers Golgi-membrane DAG; "wrong direction" has no enzymology (AGENTS.md §5). Biology lens offers the stronger precedent: Yeast-GEM r_3789 DAG[erm] ⇌ DAG[gm] is reversible, so the irreversible gm→em link is a conservative subset |
| Opening R1825-R1834 is a "decision reversal" | Corrected | No user decision exists; closed_bounds_69 and the 2026-10-26 census record them as inherited closures with unknown rationale and warn against opening by setting ub=1000 without identity/transport/demand checks. It is a science change needing scoped authorization plus identity evidence for YALI1E19647g beyond "similar to AUR1" |
| IPC-C is the recommended species | Corrected (disclosed tuning) | Scout's own text: AUR1 flips only with an IPC-or-later term; species choice is flip-driven; present ceramide / IPC / MIPC options to the lipid-line owner instead |
| Rupcic cerebroside citation (bib.irb.hr) "i.e. the delta4-desaturase has a real but unmodelled product" | Refuted in part | Link redirects to an empty page; Rupcic, Georgiu & Maric 2004 (PMID 15069670, opened) reports cerebrosides 1.3 % of total lipid (hexadecane), trihydroxy bases, NO delta-4-unsaturated base. The GlcCer branch is real and unmodelled; the delta4 inference is not supported |
| Yeast 7.6 origin of 0.000417 "recalled" | Upgraded to verified | yeast 7.6 `r_2108` lipid pseudoreaction has `0.000417 s_0535` complex sphingolipid; R1371 copies all 15 terms. Still S. cerevisiae-lumped |
| YALI1D08000g identity note "(verified)" | Corrected | Only O59715's function was verified; the gene's link comes from the model EC, the list annotation and KEGG Q6CA36 (PE3); no BLAST |
| biomass_C lipid terms list | Corrected | PA m1705[C_lp] 5.8e-5 omitted; biomass_C's ergosterol (0.0353) is 6.3× the yeast-7.6 scale (0.0056), so the placeholder is out of its own composition context |
| Demand placed in [C_go] | Flagged | Every other biomass_C lipid term is in cy/em/lp; Yeast-GEM moves IPC backbone to [c]; an ER-side term via R1187 is the alternative |
| R1825-R1834 carry AUR1 alone | Flagged | Yeast-GEM r_3484-r_3490 carry "KEI1 and AUR1"; no yli KEI1 found in KEGG; unknown, not single-gene by evidence |
| Sweep scoring and model | Corrected | Rule A only; rule-B counts missing (136/190/46/692 → inferred 145/181/47/691); SHA not in report body; newest candidate differs in biomass_C |

**Evidence.**

| Claim | URL | Status |
|---|---|---|
| O59715 = S. pombe dsd1, EC 1.14.19.17, RHEA:46544 | https://rest.uniprot.org/uniprotkb/O59715.txt | opened (scout and biology lens) |
| Yeast-GEM main r_4063 has no sphingolipid term | https://raw.githubusercontent.com/SysBioChalmers/yeast-GEM/main/model/yeast-GEM.txt | opened (both) |
| yeast 7.6 r_2108 contains 0.000417 s_0535 complex sphingolipid | https://raw.githubusercontent.com/SysBioChalmers/yeast-GEM/v7.6.0/yeast_7.6_cobra.txt | opened (biology lens) |
| Daskalaki 2018: glyco+sphingolipid 2.5-11.8 % of total lipid, fat substrates | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6226593/fullTextXML | opened (both) |
| Pomraning 2015: no sphingolipid quantities | https://pmc.ncbi.nlm.nih.gov/articles/PMC4408067/ | opened (both) |
| Rupcic 2004: cerebrosides 1.3 % of total lipid, trihydroxy bases, no Δ4 base | https://pubmed.ncbi.nlm.nih.gov/15069670/ | opened (biology lens; replaces the scout's unusable link) |
| Live CLIB122 identities Q6CD74, Q6C1K9, Q6CE86, Q6CHU1, Q6CDY7, Q6C203, Q6CA36, Q6C5P1, Q6C5K3, Q6C1X0, Q6CHP1; single-copy per KEGG KO | KEGG yli via model GeneIDs (biology lens) | opened; homology only |
| Inherited closures R1795, R1815-R1834 | `model/reports/closed_bounds_69_20260909/REPORT.md`; `model/reports/gap_metabolites_ponce2013_20261026/closed_reactions.csv` | opened (repo) |
| Joshi 2024 (PMID 39025792): sphingolipids a major lipid-body class, relative AUC only | biology lens, not opened | recalled |

**Evidence gaps.** No Y. lipolytica sphingolipid composition value (only a pooled glyco+sphingolipid range and a cerebroside fraction); species and compartment unsourced; YALI1D08000g needs BLAST and a GlcCer branch that the model lacks; no alternative KCR in yli (KEGG single-copy), so the FP is unavoidable under any VLCFA demand; a balanced Golgi MIPC/M(IP)2C synthase with genes is missing; ISC1-type hydrolysis dead ends m663/m674 untouched; the 28 zero-growth knockouts with excess route flux were not traced to their sink; V1 has no full sweep.

**The decision.** Whether a sphingolipid biomass term is authorized at all, with species, compartment and coefficient set by the lipid-unlump owner on composition evidence rather than on which gene flips. Sub-decisions: open R1825-R1834 (needs identity evidence for YALI1E19647g); add the gene-less gm→em DAG link; do the duplicate-species merge here or leave it to the colleague's line; accept the KCR FP; name the base candidate and re-score under rule B; separately authorize a BLAST for YALI1D08000g.

### 3.2 GPI anchor (6 FN)

Genes: YALI1C32810g (GWT1, Q6CAW6 reviewed, R113), YALI1E25932g (GPI12-like; model accession deleted; resolves to Q6C513 unreviewed, R371), YALI1A21970g (GPI14-like; resolves to Q6CGA3 unreviewed, R372 AND YALI1E20666g PBN1), YALI1F27343g (MCD4, Q6C0Z3 reviewed, R373), YALI1F15403g (GPI18, Q6C216 reviewed, R374), YALI1D05270g (GPI10, Q6CAB8 reviewed, R375). All PE3 "inferred from homology", strain CLIB122, mapped to YALI1 loci via `reference/locus_map/S2_table_YALI1_YALI0_map.xlsx`.

**Mechanism (verified).** The chain R_GPIGNT → R371 → R113 → R372 → R373 → R374 → R375 → R376 → m603[C_er] is connected and feasible (max R_DM_GPI 0.116617 with no growth constraint) but m603 has exactly two reactions (R376, R_DM_GPI (0,1000)); biomass_C (62 terms) has no GPI term, so every chain reaction has FVA 0/0 at the optimum. The R376 OR pair YALI1E04744g/YALI1E14925g (both on the list) stays FN under every variant. Yeast-GEM v9.0.2 has the same structure: s_0332 is a dead end with no pseudoreaction term, and the R372 AND / R376 OR rules are copies of Yeast-GEM r_0519/step 9, not Y. lipolytica evidence.

**Recommended design (b1 or b2, in memory).**

| Edit | Term / reaction | Value |
|---|---|---|
| biomass term | biomass_C += m603[C_er] | −1.6962495847474773e-4 mmol/gDW |
| b2 only | R_DM_GPI | (0,1000) → (0,0) |

b1 (exit left open) gives the identical optimum, identical knockouts, R_DM_GPI = 0 at the optimum, EGC optimal 0.0; it is the no-reversal alternative. Variant a (R_DM_GPI lb 5.2865e-5) is rejected: the family knockouts and the EGC LP become infeasible (benchmark.md certification condition 2 forbids a forced non-maintenance lower bound).

**Verified numbers.**

| Quantity | Baseline | b2 (= b1) | Reproduced by |
|---|---|---|---|
| WT growth | 0.3116608216336117 | 0.31152321228288754 (−0.0442 %) | LP lens (rel diff 3.6e-15); biology lens spot |
| Chain flux at optimum | 0 | 5.284211e-5 on every chain reaction = coef × μ; R_DM_GPI 0 | LP lens |
| KO/WT 6 family genes | 1.0000000000000004 | 0.0 (optimal) | LP lens (context and fresh-model) |
| KO/WT YALI1E20666g (PBN1) | 1.0 | 0.0 (new FP) | LP lens |
| KO/WT YALI1E04744g, YALI1E14925g (R376 OR) | 1.0 | 1.0 | LP lens |
| EGC | optimal 0.0 | optimal 0.0 (also free biomass 0.0, free R_DM_GPI 0.0) | LP lens |
| Sweep 15 % (mt-excl) | 125/190/43/691, 16 unresolved | 131/184/44/690, 16 unresolved | LP lens (exact at 1/5/10/15 %) |
| Sweep 1 % / 5 % / 10 % | 97/218/36/698; 104/211/38/696; 108/207/43/691 | 103/212/37/697; 110/205/39/695; 114/201/44/690 | LP lens |
| Class changes | — | exactly 7; no other KO/WT moves > 1e-3 | LP lens |
| Variant a WT; knockouts; EGC | — | 0.31152315103774453; 7 infeasible; EGC infeasible | LP lens |
| Sensitivity (variant c, L = 1e-6…1e-3) | — | same 7 flips at every level; WT cost linear (2.6e-6 … 2.6e-3) | scout only |
| pFBA co-product shifts (PS synthase +55 %, CDP-DAG synthase +73 %, R1773 → 0) | — | — | scout only |

**Refuted or corrected by the verifiers.**

| Scout claim | Verdict | Evidence |
|---|---|---|
| The coefficient's mmol/gDW basis is compatible with the biomass lipid terms | Refuted (major) | Biology lens: the term consumes 1.70e-4 PI vs the biomass PI term 7.5e-5 (2.3×), 3.39e-4 PE vs 2.37e-4 (1.4×), returns 3.39e-4 DAG vs TAG 2.34e-4; the model's phospholipid terms are 1-2 orders below a physiological basis (ergosterol 0.035 vs PC 3.46e-4). The GPI term would be the largest phospholipid consumer. Classes are unaffected; the value needs the lipid owner's rescaling decision |
| Cross-check 0.21 × 0.04 / 50 kDa confirms 1.9e6 | Corrected | Circular: Yin derived 1.9e6 from the same BCA peptide measurement and the same 50 kDa assumption |
| "No UniProt" for GPI14-like; GPI12 unverifiable | Corrected | YALI1E25932g → Q6C513 (EC 3.5.1.89, YALI0_E21934g); YALI1A21970g → Q6CGA3 (GPI mannosyltransferase 1, YALI0_A20922g; already in `EGC-126fc70a7b17.json`); both homology-only |
| PBN1 dispensability "untested here" | Corrected | `cas9_cas12a_fitness.csv`: PBN1 fitness −2.608 (Cas9, q 0.183) and −1.821 (Cas12a, q 0.063), both just below the call cutoff; PBN1 is essential in S. cerevisiae (PubMed 9649520 title via P25580). The FP is a borderline screen call on a gene expected to be essential; still FP under benchmark.md |
| b2 preferred over b1 because of the SD-Leu PE→DAG sink | Corrected | The nine SD-Leu knockouts are infeasible on Task A + PO1f on every variant, so the sink argument cannot be reproduced on this condition; b1 and b2 give identical results; only b2 reverses the 2026-10-21 decision |
| "the change is made only to raise the match of six listed genes, it must be justified as a composition fact" | Refuted as framing | Names the motive AGENTS.md §2 forbids; lead with census cause d (product not demanded), FN-causes candidate 5 and Jaafar & Zueco 2004 |
| The GPI guards would trip if m603 is added to biomass_C | Corrected | `gpi_chain_links.py`/`gpi_coproducts.py` compare reactions_before at their own step; a later step cannot trip them. The real constraint is the biomass_C composition pin of every biomass-editing step (`cofactor_biomass_terms.py` line 7) |
| Yeast-GEM has no GPI term; KEGG R12351/R05922 chemistry (recalled) | Upgraded to verified | Yeast-GEM v9.0.2 s_0332 only in r_0523, no pseudoreaction; R12351 releases DAG, R05922 Dol-P (opened) |
| PMC URL for Yin 2007 | Corrected | Not machine-readable; use https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2040195/fullTextXML (opened by the biology lens) |
| m603 identity listed only as an evidence gap | Corrected | ChEBI:53060 names three EtNP groups while R373 + R376 deliver two; no formula; GPI7 step missing; STATE.md 2026-10-22 lists these as open. Precondition to pricing (cofactor precedent resolved identity before pricing) |
| Sweep rule and condition | Corrected | Rule A only (primary rule-B baseline 136/190/46/692 not reported); condition file omits the PO1f plasmid LEU2 complement (FP 43 vs 42); SHA not in report body; newest candidate c1f70ac4 has a different biomass_C |

**Evidence.**

| Claim | URL | Status |
|---|---|---|
| Yin 2007: 1.9e6 CWP/cell (50 kDa assumed), 18.6 pg/cell, wall 21 % of biomass, peptide 4 % of wall, "most" GPI-modified | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2040195/fullTextXML | opened (biology lens; scout's PMC link returned a verification page) |
| Jaafar & Zueco 2004: YlCWP1, first Y. lipolytica GPI-CWP, no quantity | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=14702397&rettype=abstract | opened (biology lens) |
| da Silva 2025: W29 wall 70 % carbohydrate / 15 % protein / 5 % lipid; mannoprotein 13.15 % | `$FN29/gpi/mdpi_pdf/text.txt` (processes-13-00615.pdf) | opened (scout, biology lens) |
| Vega & Dominguez 1986 wall composition | https://doi.org/10.1007/bf00414721 | not opened (Springer login) by anyone |
| UniProt Q6CAW6, Q6C0Z3, Q6C216, Q6CAB8, Q6C5K1, Q6C741 reviewed, PE3 | https://rest.uniprot.org/uniprotkb/<acc>.txt | opened (both) |
| Yeast-GEM v9.0.2: s_0332 dead end, no biomass GPI term, r_0519 "YCL052C and YJR013W" | `$IYALI26_RESEARCH_ROOT/reference/external_models/yeast_gem/Yeast-GEM.xml` | opened (biology lens) |
| KEGG R12351 (DAG release), R05922 (Dol-P release) | https://www.kegg.jp/entry/R12351, https://www.kegg.jp/entry/R05922 | opened (biology lens) |
| 2026-10-21 decision "Keep the GPI exit open as it is" and "A future GPI biomass term … would reopen the question" | `model/STATE.md` (origin/main lines 2167-2183; worktree lines 873-893) | opened (all three lenses) |
| PBN1 screen scores | `state/essentiality/repository/cas9_cas12a_fitness.csv` | opened (biology lens) |

**Evidence gaps.** No Y. lipolytica GPI-anchored-protein quantity in any opened source; the S. cerevisiae proxy includes Pir-CWPs, excludes PM-GPI proteins, assumes 50 kDa, YPD exponential; its scale basis conflicts with the model's lipid pool; m603 is a mass-inconsistent species (three-EtNP name, two delivered, no formula); GPI-GnT subunits (GPI1/2/3/15/19) on the list stay FN because R_GPIGNT has no GPR; SMP3/GPI7 orthologs unchecked; the third screen's raw PBN1 call not read.

**The decision.** Whether GPI-anchored protein is a growth requirement to represent in biomass_C (reopening what the 2026-10-21 entry said a biomass term would reopen). Sub-decisions: accept the S. cerevisiae proxy and how to re-express it on the model's lipid basis (lipid owner); b1 (no reversal) or b2 (close R_DM_GPI); accept the PBN1 FP or review the R372 AND rule (and the R376 OR rule) first; settle m603's identity before pricing; base candidate, rule B, real PO1f profile; order relative to the PI/PE/DAG unlump.

### 3.3 Cardiolipin (2 FN)

Genes: YALI1E17162g (YALI0E14102g, Q6C5Y3, EC 2.7.8.5, mitochondrion by similarity, PE3; "similar to PGS1"; R1836), YALI1F31283g (YALI0F23837g, Q6C0K8, PANTHER SF60 cardiolipin synthase, the only SF60 entry in yli, PE3; R1721). FP risk YALI1A13994g (YALI0A14157g, Q6CH12, PF09419 PGP phosphatase, PE4, the only PF09419 protein in CLIB122; R1722). No verified names, no Y. lipolytica biochemical evidence (curation ledger: model/GPR assignment only).

**Mechanism (verified).** biomass_C has no cardiolipin or PG term; cardiolipin m1651[C_mm] is consumed only by the legacy xLIPID (0.842 per 1000 "lipids" → xBIOMASS, not the objective, both [0,0] at the optimum), R1723 and R1724 (blocked, UM10 pools without transports). The route R1791 → R197 → R1836 → R1722 → R1721 is complete per compartment, irreversible and unbranched, and unused only because nothing pulls on cardiolipin. UM10 does not need connecting: in every variant R1723/R1724 stayed 0 while the term was met. Yeast-GEM imposes no cardiolipin requirement (r_4063 has no cardiolipin backbone; S. cerevisiae crd1/pgs1 nulls are viable), so the term must be justified from Y. lipolytica alone (consensus-essential PGS1-/CRD1-like orthologs in an obligately respiring yeast) — the reason the list calls them essential is not in any opened source.

**Recommended design (variant a, in memory).**

| Edit | Term | Value |
|---|---|---|
| biomass term | biomass_C += m1651[C_mm] cardiolipin | −7.161586749409294e-4 mmol/gDW (0.08 × 12.98 mg/gDW / 1449.958 g/mol); corrected MW for CL 72:7 (C81H144O17P2, 1451.97) gives 7.152e-4; range 6.27e-4 (7 %) to 8.06e-4 (9 %) |

Alternatives: alt (i) 1.8249e-3 (Hein & Hayen ratio on the Yeast-GEM r_4063 four-class basis, S. cerevisiae, cross-species); alt (ii) 2.1207e-5 (same ratio on the model's own PC+PE+PI+PS placeholders, internally consistent). Variant b (forced demand DM_cardiolipin_mm lb 2.232e-4) is rejected: the four route knockouts become infeasible and the EGC LP is infeasible with the bound in force.

**Verified numbers.**

| Quantity | Baseline | Variant a | Reproduced by |
|---|---|---|---|
| WT growth | 0.3116608216336117 | 0.3110599733571873 (−0.19 %) | LP lens (exact); biology lens spot |
| Route flux at optimum | 0 | R1791 = R197 = 4.4553659669e-4; R1836 = R1722 = R1721 = 2.2276829835e-4 = coef × μ; R1723, R1724, R1742, xLIPID, xBIOMASS 0 | LP lens |
| KO/WT YALI1E17162g, YALI1F31283g | 1.0 | 0.0 (optimal; fresh-copy re-solve) | LP lens, biology lens |
| KO/WT YALI1A13994g | 1.0 | 0.0 (new FP) | LP lens, biology lens |
| KO/WT YALI1F08225g (CTP/CMP carrier, on list) | 1.0 | 0.0 (side-effect TP) | LP lens, biology lens |
| KO/WT YALI1E17431g (CDP-DAG synthase) | 0.0 | 0.0 (already essential; not a side effect) | LP lens |
| KO/WT YALI1D22083g, YALI1C19814g, YALI1F11605g | 1.0 | 1.0 | LP lens |
| EGC | optimal 0.0 | optimal 0.0 (182 exchanges closed; probe demand on m1651 also 0.0) | LP lens |
| Sweep 15 % (mt-excl) | 125/190/43/691, 16 unresolved | 128/187/44/690, 16 unresolved, same 20 infeasible | LP lens (fully independent after-sweep; 205 baseline re-runs, zero mismatches) |
| Sweep 1 % / 5 % / 10 % | 97/218/36/698; 104/211/38/696; 108/207/43/691 | 100/215/37/697; 107/208/39/695; 111/204/44/690 | LP lens |
| Class changes | — | exactly 4 (3 FN→TP, 1 TN→FP); 138 other ratios move ≤ 5e-4 without class change | LP lens |
| Variant b WT; knockouts; EGC | — | 0.3110588 (exact bound 0.3110588127); infeasible; infeasible / 0.0 relaxed | LP lens |
| Sensitivity c = 1e-6…1e-1; range ends; alt bases | — | family 0.0 at every value; WT 0.31166 → 0.24546 | scout only |

**Refuted or corrected by the verifiers.**

| Scout claim | Verdict | Evidence |
|---|---|---|
| Ownership: "cardiolipin expansion is an open gate item" | Corrected (major) | `model/curation/lipid_unlump/lipid_unlump_sn_core_curation.json` cardiolipin_four_chain_audit.biomass_contract: biomass_C `candidate_action: no_change`, `cardiolipin_biomass_weight_decision: blocked`, sekova_profile_fraction `absolute_mol_per_gDW_available: false`, blocker `cardiolipin_biomass_weight_not_source_supported`; also blockers `chain_resolved_mitochondrial_pg_supply_absent`, `R1791_R1794_yarrowia_transporters_unverified`, and TAM41-like YALI1C17238g assigned without a model role. The recommendation asks the owner line to reverse a recorded judgment |
| C81H142O17P2 = CL 72:7, 1449.958 g/mol | Corrected | PubChem CID 5283539 is tetralinoleoyl CL 72:8; CL 72:7 is C81H144O17P2 = 1451.97; coefficient 7.152e-4 (0.14 % change) |
| Yeast-GEM carries cardiolipin in r_2108 (recalled) | Refuted | r_2108 is the lipid merge pseudoreaction; r_4063 has no cardiolipin backbone; s_3738/s_3739 enter only SLIME reactions and r_4061. No Yeast-GEM precedent for the term |
| The Yeast-GEM four-class basis could not be traced | Refuted | PC 0.025783, PE 0.0069293, PI 0.0069103, PS 0.0059508 are exactly the r_4063 coefficients (S. cerevisiae); alt (i) is sourced but cross-species |
| The two Y. lipolytica sources agree | Corrected | Sekova implies CL ~14-18 % of phospholipids (TLC, stationary, glycerol); Hein & Hayen 2.8 % of GPs (HPLC-MS, exponential, glucose); 7.16e-4 makes cardiolipin 2.07× the PC placeholder by mole, the largest glycerophospholipid term in biomass_C |
| "R_NTP1" | Clarified | No such reaction exists in e1e7c868; it is an in-memory ATP-hydrolysis reaction added for the check (same for every family) |
| Sweep rule and model | Corrected | Rule A (16 unresolved, COX9 in scope); rule B inferred 136/190/46/692 → 139/187/47/691; SHA absent from report body; guard step must pin the newest biomass_C (2bc263f2) |

**Evidence.**

| Claim | URL | Status |
|---|---|---|
| Sekova 2019: W29 membrane lipids 12.98 mg/gDW; CL 7-9 % of membrane-lipid fractions; stationary, 1 % glycerol; phospholipids ~50 %, sterols 45 % | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6950712/fullTextXML | opened (scout, biology lens) |
| Hein & Hayen 2012 Table 1: Y. lipolytica CA 2.8 ± 0.3 % of GPs (glucose, 1.7 g CDW/L) | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3901198/fullTextXML | opened (both) |
| Abreu 2022: CL C72:7 / C72:6 predominant (engineered obese strain, peak areas) | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9394261/fullTextXML | opened (both) |
| PubChem CID 5283539 = tetralinoleoyl CL (72:8) | PubChem (biology lens) | opened |
| UniProt Q6C5Y3, Q6C0K8, Q6CH12; PF09419 single in CLIB122 | rest.uniprot.org (biology lens) | opened |
| Yeast-GEM 9.1.1 r_2536/r_2542/r_2548 (mm, lb 0); r_4063 no cardiolipin | yeast-GEM yml (biology lens) | opened |
| Lipid-unlump ledger and PREFLIGHT | `model/curation/lipid_unlump/lipid_unlump_sn_core_curation.json`; `model/reports/cfsa_lipid_unlump_20260908/PREFLIGHT.md` | opened (scout, project-rules lens) |
| Census UM10 "pools without transports, transports unblock 12/13, sinks 0/13" | `model/reports/gap_metabolites_ponce2013_20261026/REPORT.md` item 10 | opened |

**Evidence gaps.** No direct mmol/gDW cardiolipin value; Sekova's denominator implied, not stated; condition mismatch (stationary glycerol vs exponential glucose); MW from an engineered strain; no Y. lipolytica evidence that the three genes catalyse these steps or sit in the mitochondrial membrane; why PGS1/CRD1 are essential in Y. lipolytica while S. cerevisiae nulls are viable is uncited; no second PGP phosphatase in CLIB122 at domain-search level, so the FP cannot be removed by a GPR member; whether CDP-DAG for cardiolipin is made in situ (TAM41-like YALI1C17238g) or imported is unevidenced (drives the YALI1F08225g side effect).

**The decision.** Whether a lumped m1651[C_mm] cardiolipin term may enter biomass_C at all, which means instructing the lipid-unlump owner to lift the recorded "blocked" biomass-weight decision or explicitly authorizing this line to edit biomass_C lipid composition. Sub-decisions: Sekova-chained 7.15e-4 vs the placeholder-consistent 2.12e-5 (classes identical, growth cost differs); accept the YALI1A13994g FP and the YALI1F08225g side effect; base candidate and rule-B re-run.

### 3.4 N-glycan (3 FN)

Genes: YALI1B21618g (ALG7-like; model accession A0A1H6PQ85 deleted; = YALI0B16588g Q6CED2, EC 2.7.8.15, KEGG K01001; R2106), YALI1E40471g (ALG1_YARLI Q6C3K2 reviewed, by similarity; R173), YALI1F28833g (GMT_YARLI Q6C0U0 reviewed, VRG4; R850). FP-risk: YALI1C25093g (= YALI0C17941g Q6CBK6 GDA1-like, EC 3.6.1.42; R602/R586), YALI1D26639g (= YALI0D21032g Q6C8B5 CWH8/CAX4-like, EC 3.6.1.43; R2181/R2177); both model accessions deleted.

**Mechanism (verified).** (1) The ER chain R2108 → R2107 → R2106 → R2105 ends at m1969[C_er] (no consumer, UM18 root); R173 (ALG1) sits in the cytosol on a different GlcNAc2-PP-Dol species m344 (no producer) making m345 (no consumer); R141 is a third unlinked pair in the Golgi; no Man9/Glc3Man9 LLO, no OST, no N-glycan demand. (2) The mannan term (0.235849 mmol Man4/gDW, no recorded provenance; 0.235849 × 4 = 0.943396 = the glucan coefficient, so the term was evidently re-expressed as Man4 from Yeast-GEM's 0.70204 C6H10O5 without changing R289 — inferred) is met entirely by the PMT lump R289 (Dol-P-Man → Dol-P + H+ + Man4), which creates C18H33O16 and charge +1 per turnover (hand-balanced by scout and biology lens). Hence Golgi GDP-mannose is never needed and R850 is blocked; the gene-less R1137 (GDP-Man cy ⇌ go) duplicates R850; Golgi Pi m837[C_go] has no exit so the GDPase R602 cannot run.

**Design C2 (A2 + B, in memory) — as recommended by the scout; the verifiers split it.**

| Part | Edit | Bounds / value |
|---|---|---|
| A2 | R_LLO_link: m1969[C_er] → m344[C_cy] | (0,1000), gene-less |
| A2 | R_LLO_man9: m345[C_cy] + 4 m343 GDP-Man + 8 m10 H+ + 4 m502[C_er] Dol-P-Man → m_LLO_M9[C_er] + 4 m268 GDP + 4 m503 Dol-P | (0,1000), gene-less (balanced under model formulas) |
| A2 | R_LLO_glc3: m_LLO_M9 + 3 m1989 Dol-P-Glc → m_LLO + 3 m503 | (0,1000) |
| A2 | R_OST: m_LLO → m_NGLYC[C_er] + m2021 Dol-PP; R_OST_M9: m_LLO_M9 → m_NGLYC + m2021 | (0,1000), gene-less |
| A2 | biomass_C += −0.01 m_NGLYC[C_er] | placeholder, NO SOURCE |
| B | R289 rewritten 1:1: m502[C_er] → m503[C_er] + m_OMAN1[C_er] (C6H10O5); Man4 product and spurious H+ removed | GPR/bounds unchanged |
| B | R_OMAN1_t: m_OMAN1[C_er] → m_OMAN1[C_go] | (0,1000), gene-less |
| B | R_MANNAN_GO: m_OMAN1[C_go] + 3 m290[C_go] GDP-Man + 6 m10 H+ + m835[C_go] H2O → m1324[C_cy] mannan + 3 m292[C_go] GDP | (0,1000), gene-less (machine-balanced) |
| B | R_PIt_go_cy: m837[C_go] → m35[C_cy] | (0,1000), gene-less (without it WT = 0) |
| B | R1137 | (−1000,1000) → (0,0) (necessary for the VRG4 flip; probe B0) |

Yeast-GEM-faithful alternative for the mannan half (biology lens): R289 1:1 with the mannan species as C6H10O5 at 0.9434 fixes the mass creation without making VRG4 or GDA1 essential; routing mannose through the Golgi is then a separate Y. lipolytica-specific decision (O-glycans Man2-4, Park 2011).

**Verified numbers.**

| Quantity | Baseline | B alone | A2 alone | C2 | Reproduced by |
|---|---|---|---|---|---|
| WT growth | 0.3116608216336117 | 0.2807701548627785 (−9.9 %) | 0.3062768297726763 | 0.27639305844315876 (−11.3 %) | LP lens (all exact) |
| KO/WT YALI1B21618g, YALI1E40471g | 1.0 | 1.0 | 0.0 | 0.0 (optimal) | LP lens |
| KO/WT YALI1F28833g | 1.0 | 0.0 (B0 with R1137 open: 1.0) | 1.0 | 0.0 | LP lens |
| KO/WT YALI1C25093g (GDA1-like) | 1.0 | 0.0 | 1.0 | 0.0 (FP); 0.999999999999996 with GPR "or YALI1E39570g" | LP lens; biology lens (2 LPs) |
| KO/WT YALI1D26639g (CWH8-like) | 1.0 | 1.0 | 0.0 | 0.0 (FP) | LP lens |
| KO/WT YALI1D07696g (ALG5-like), YALI1E41893g, ALG2/11-like and other watch genes | 1.0 | — | 1.0 | ~1.0 (E41893g 0.99995) | LP lens |
| Route flux at C2 optimum | — | — | — | N-glycan chain 2.76393e-3; mannan route 6.5187e-2; R850 = R602 = R_PIt_go_cy 1.95561e-1; R1137 0/0; R2102/R288 interchangeable 0…7.6243e-2 | LP lens (FVA) |
| Dolichol leak (Dol-P term removed) | — | R287 [0,0] | R287 [0,0] | R287 [0,0]; R2181 [2.7076e-3, 2.7641e-3] (not strictly fixed, within slack) | LP lens |
| EGC | optimal 0.0 | optimal 0.0 (scout) | optimal 0.0 (scout) | optimal 0.0; free-source probes on R_OST_M9, R_OST, R_MANNAN_GO, R_LLO_link all 0.0 | LP lens |
| Sweep 15 % (20 infeasible as 0; mtDNA+COX9 in scope) | 136/190/55/700 | not swept | not swept | 139/187/57/698 | LP lens (full C2 sweep; 196 baseline re-runs, zero disagreements) |
| Sweep 1 % / 5 % / 10 % | 108/218/45/710; 115/211/47/708; 119/207/52/703 | — | — | 111/215/47/708; 118/208/49/706; 122/204/54/701 | LP lens |
| Class changes | — | — | — | exactly 5 (+3 TP, +2 FP); 154-159 ratios move > 1e-6, none crosses a cutoff | LP lens |
| Coefficient probes c_n = 1e-3, 1e-2, 1e-1, 0.1179 | — | — | A2: 0.31111, 0.30628, —, 0.25816 (c_n 0.1179) | ALG7/ALG1 KO/WT 0 at every c_n | LP lens (corrected: scout's listed values were variant-A numbers) |

**Refuted or corrected by the verifiers.**

| Scout claim | Verdict | Evidence |
|---|---|---|
| C2 is "solvable after decision" for all three genes | Refuted (major, two lenses) | The N-glycan demand has no source anywhere (Yeast-GEM has no N-glycan term; Park 2011, De Pourcq 2012, Orlean 2012 full text give none; Vega & Dominguez 1986 unread); the flip is forced by the term's existence (probes 1e-3…0.1179), the tuning signature; AGENTS.md §2. ALG7/ALG1 are blocked on composition evidence; only VRG4 (variant B) is decidable |
| R1137 closure is "evidence in hand" | Corrected | Original-model reaction with empty notes and no STATE.md entry; biology is S. cerevisiae-only (Dean 1997); probe B0 shows it is required for the flip (recall motive, AGENTS.md §6). Yeast-GEM gives both Golgi GDP-Man entries (r_1165, r_1803) the VRG4 gene, which supports closing R1137 or giving it VRG4, but it is a decision |
| Variant B's 1:3 split and lumps are evidence-based | Corrected | Scout's own words: a design choice; Yeast-GEM makes all mannan 1:1 from Dol-P-Man (r_0362, 0.70204 C6H10O5) and VRG4 never feeds biomass; three gene-less lumps; new m504/R1250 dead end; WT −9.9 to −11.3 % |
| No alternative Golgi GDP→GMP route exists | Refuted | YALI1E39570g (YALI0_E33407g, Q6C3N3, EC 3.6.1.42, "similar to GDA1"), absent from the model and the list; with R602/R586 "YALI1C25093g or YALI1E39570g" the GDA1 FP vanishes (orthology level, needs reciprocal BLAST and localization) |
| UniProt A0A1H6PQ85, A0A1H6PR29, A0A1D8NFJ0 cited as live | Corrected | All deleted; replace with Q6CED2, Q6CBK6, Q6C8B5 via the locus map |
| Orlean 2012 full text (PMC3512495) not reachable; mannoprotein 30-50 % of wall recalled | Corrected | PMC3512495 is an unrelated BMC article; full text reached at https://academic.oup.com/genetics/article/192/3/775/5935161: O-mannosylation (first Man from Dol-P-Man, up to four more by Ktr1/Mnn1 from GDP-Man) and OST (GlcNAc2Man9Glc3) confirmed; Cwh8 recycling "likely"; the 30-50 % figure is not in it and stays unverified |
| STATE.md dolichol entry at lines 1029-1077 | Corrected | origin/main lines 1108-1152 (content correct) |
| A2 block sensitivity WT values | Corrected | They are variant-A numbers (obligatory Glc3); the 0.1179 probe was run at 0.1179245 |
| 31 reactionless genes | Corrected | 32 |
| CWH8-like FP "needs an alternative Dol-PP phosphatase" | Corrected (incomplete) | Dol-PP has no fate but R2181/R2177 (same gene); S. cerevisiae cax4 nulls survive with de novo Dol-P (SEC59 R287 + R_DOLSYN exist), so a Dol-PP sink is the other route |
| Sweep totals follow the primary rule | Refuted | 55/700 includes 16 mtDNA + COX9 (rule B baseline 136/190/46/692 → inferred 139/187/48/690); the 2026-10-28 1e-6 ratio-move gate would stop this build (154-159 moves) and was not mentioned |
| Baseline is current | Corrected | e1e7c868; line moved to 2bc263f2 (rule B 142/184/44/694) |

**Evidence.**

| Claim | URL | Status |
|---|---|---|
| SGD: ALG7, ALG1, VRG4, STT3 essential; GDA1, CAX4/CWH8, ALG5 non-essential; HVG1 viable | https://www.yeastgenome.org/backend/locus/S000000447, S000000314, YGL225W, YGL022W, S000000768, YGR036C, YPL227C | opened (both) |
| Dean 1997: vrg4 defective in Golgi GDP-mannose transport | https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:9395539%20AND%20SRC:MED&resultType=core&format=json | opened (both) |
| Park 2011: WT N-glycans Man7-12GlcNAc2, O-glycans Man2-4, no quantity | https://pmc.ncbi.nlm.nih.gov/articles/PMC3067234/ | opened (both) |
| De Pourcq 2012: Man8/Man9GlcNAc2 mannoproteins, no quantity | https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0039976 | opened (both) |
| Orlean 2012 full text | https://academic.oup.com/genetics/article/192/3/775/5935161 | opened (biology lens); "30-50 %" not found |
| Vega & Dominguez 1986 | https://produccioncientifica.usal.es/documentos/61579707d3e35f5001a2981b?lang=en (record only) | not opened (article) |
| KEGG R04502 (ALG1); BiGG GDPMANNtg; KEGG yli K01001/K03842/K15356/K01526/K07252; UniProt Q6C3K2, Q6C0U0, Q6CED2, Q6CBK6, Q6C8B5, Q6C3N3 | biology lens | opened |
| Yeast-GEM v9.0.2: r_0362 PMT 1:1, r_4048 mannan 0.70204, no OST/N-glycan demand, r_1165/r_1803 VRG4, r_4607 Golgi Pi, r_0810 GDA1 | `$IYALI26_RESEARCH_ROOT/reference/external_models/yeast_gem/Yeast-GEM.xml` | opened (biology lens) |
| Census UM18; STATE.md 2026-10-28 dolichol entry; STATE.md 2026-09-21 mannan-via-R289 note | `model/reports/gap_metabolites_ponce2013_20261026/REPORT.md` §4 item 18; `model/STATE.md` (origin/main lines 1108-1152, 4016) | opened |

**Evidence gaps.** No N-glycan content per gDW in any yeast source; mannan term provenance unknown; the Dol-P-Man:GDP-Man split unsourced; no Y. lipolytica phenotypes for gda1-, cax4-, alg5-type mutants; YALI1E39570g needs reciprocal BLAST and localization; ALG2/ALG11-like genes (on the list) stay FN inside the 6-member OR of R141; OST subunits STT3/OST1/OST2/SWP1 on the list but absent from the model; dolichol-family formulas inconsistent (R173, R141 imbalanced in the source model), so balance of the new steps holds only under model conventions.

**The decision.** For ALG7/ALG1: none available — defer until a composition source for N-glycosylated protein is read (Vega & Dominguez 1986 is the only Y. lipolytica wall-composition primary found), or name an unsourced structural term as an explicit exception. For VRG4: whether to authorize the R289 1:1 rewrite (evidence-backed) together with the Golgi mannan route (1:3 split, three gene-less lumps, Golgi Pi exit, R1137 closure; WT −9.9 %) or only the Yeast-GEM-faithful 1:1 fix (no VRG4 flip). Sub-decisions: accept the GDA1-like FP or add YALI1E39570g after BLAST; accept the CWH8-like FP (A2 only); set the KO/WT ratio-move gate; dolichol-formula clean-up as prerequisite; base candidate and rule B.

### 3.5 Mitochondrial translation (8 FN)

Genes (all single-copy mitochondrial aminoacyl-tRNA synthetases per the model GPR; identities computational only, UniProt ARBA / PANTHER / GO-IBA / NCBI similarity notes; the 12 model UniProt accessions are all Inactive): YALI1C25782g Phe (Q6CBI1, R626), YALI1D04566g Met (Q6CAE1, R548), YALI1D30570g Lys (Q6C814, R532), YALI1E05833g Asn (XP_503562.2, PRK03932, similar to SLM5; R158), YALI1E21230g Glu (XP_504083.3, similar to MSE1; R338), YALI1E34408g Tyr (Q6C474, R776), YALI1F33819g Asp (Q6C0A5, R166), YALI1F03527g Leu (XP_504899.3, similar to NAM2; R523). S. cerevisiae names (MSF1 etc.) are ortholog names, not verified Y. lipolytica names.

**Mechanism (verified).** All 15 mitochondrial charging reactions (14 synthetases + FMT R546) are blocked at WT: each charged tRNA[mi] has its synthetase as its only reaction, each free tRNA[mi] has no producer. The cytosolic synthetases couple to biomass_C through 20 trna_biomass_residue terms (sum 3.726504 mmol/gDW); the mitochondrial ones through nothing. Yeast-GEM has the same 17 dead-end charging reactions and no translation reaction (verified), so the fix has no Yeast-GEM precedent. Respiration is required on Task A (complex III/IV closed: infeasible; ATP synthase closed: 0.110 of WT; complex I closed: 0.758 of WT; verified). Secondary gaps: no mitochondrial Met import (R1251 gene-less, mi→cy only; the only Met[mi] source is 1.2466e-5 from the lumped lipoyl synthase R2287); no cytosol→mitochondrion Ile import (R1218, mi→cy only); AMP[mi] consumed only by R122/R123 (YALI1F33951g) and R905; no mitochondrial Gln species. The mtDNA (AJ307410, W29, 47,916 bp) encodes tRNAs for 19 amino acids (no tRNA-Val annotated; Ala/Cys/Gln/Gly/Pro/Ser present), so the model's 14 modules are a Yeast-GEM inheritance, not organism-derived (biology lens).

**Recommended design (a2 + R1251 reversible + R1218 reversible; the scout never ran this combination, the LP verifier did).**

| Edit | Reaction / term | Value |
|---|---|---|
| new species | mtprot[C_mi] | — |
| MITO_TRANSLATION_14 | consumes X-tRNA(X)[mi] and returns tRNA(X)[mi] at coefficient f_X / 0.7512359550561798 for the 14 represented residues (AJ307410 fractions: Val 0.06045, Arg 0.01978, Asn 0.06494, Asp 0.02360, Glu 0.02202, His 0.01888, Ile 0.12899, Leu 0.15528, Lys 0.02831, Phe 0.07461, Thr 0.05596, Trp 0.01348, Tyr 0.05685, Met 0.02809); no fMet, no free amino acids; → 1 mtprot | (0,1000), gene-less |
| biomass term | biomass_C += −0.01863252 mtprot (= 3.726504 × 0.005) | placeholder, NO SOURCE |
| R1251 L-methionine mi→cy | (0,1000) → (−1000,1000) | direction change, no mechanism evidence |
| R1218 L-isoleucine mi→cy | (0,1000) → (−1000,1000) | direction change, no mechanism evidence |

Variant a (composition-weighted, fMet 15/4450, six free mitochondrial amino acids) was the only one the scout fully swept (+6 FP). Variant b (equimolar, 1e-5) rides on the lipoyl-synthase Met artefact. Variant c2 (coupling to the respiratory complexes) makes every synthetase knockout infeasible; the scout rejected it as "unresolved", which is wrong under the 2026-10-30 rule (see corrections), but it has no composition basis either.

**Verified numbers.**

| Quantity | Baseline | Variant a | Recommended (a2 + R1251 + R1218) | Reproduced by |
|---|---|---|---|---|
| WT growth | 0.3116608216336117 | 0.31081411751906585 | 0.31067012284218265 (−0.32 %); translation flux 5.7886e-3 | LP lens (exact); biology lens spot (44 LPs) |
| WT without R1251 reversed | — | 0.0 | 0.0 (a2) | LP lens |
| Met[mi] sink max | 1.2466e-5 (R1251 as is) / 1.5319 (reversible) | — | — | LP lens |
| KO/WT 8 family genes | 1.0 | 0.0 | 0.0 (optimal, closed-reaction flux 0) | LP lens, biology lens |
| KO/WT YALI1A09121g Thr aaRS (on list) | 1.0 | 0.0 (TP) | 0.0 (TP) | LP lens, biology lens |
| KO/WT YALI1C15180g Ile, YALI1C25598g Trp, YALI1F33951g adenylate kinase | 1.0 | 0.0 (FP) | 0.0 (FP) | LP lens, biology lens |
| KO/WT YALI1F10110g FMT, YALI1D08125g Ala transaminase | 1.0, 0.99999974 | 0.0, 0.0 (FP) | 1.0, 0.9999997 (no change) | LP lens, biology lens |
| KO/WT ILV2/ILV5/ILV3-like (YALI1C00293g, YALI1D03952g, YALI1C32276g) | 0.2270 | 0.0 (artefact); 0.22518 with R1218 reversible | 0.22448 (class unchanged) | LP lens, biology lens |
| EGC | optimal 0.0 | optimal 0.0 | optimal 0.0 (incl. R1218; no-uptake biomass and translation 0.0) | LP lens |
| Sweep 15 % (rule A, 1,081 genes) | 125/190/46/700, 20 unresolved | 136/179/52/694, 20 unresolved | 134/181/49/697, 20 unresolved (same genes) | LP lens (full sweep of the recommended variant; variant a scout only) |
| Sweep 1 % / 5 % / 10 % | 97/218/36/710; 104/211/38/708; 108/207/43/703 | 108/207/42/704; 115/200/44/702; 119/196/49/697 | 106/209/39/707; 113/202/41/705; 117/198/46/700 | LP lens |
| Class changes | — | 17 | exactly 12 (9 FN→TP, 3 TN→FP), each re-derived with an own baseline LP | LP lens |
| Extra FP probe (10 supply/PPi/AMP-adjacent genes) | — | — | no class change | biology lens |
| AJ307410 composition (15 canonical CDS, 4450 residues, SHA a88dab8a…) | — | — | identical recount | LP, biology and project-rules lenses |
| Respiration controls (R305/R304 infeasible; R171 0.0342; R2062 0.2362) | verified | — | — | LP lens |
| Coefficient scan (fractions 1e-4…1e-1) | — | family 0.0 at every value; WT 0.31164 → 0.29556 | — | scout only |

**Refuted or corrected by the verifiers.**

| Scout claim | Verdict | Evidence |
|---|---|---|
| "Solvable after decision" | Refuted (biology minor; project rules major) | The biomass coefficient has no Y. lipolytica, yeast or Yeast-GEM source; R1251/R1218 reversals need enzymology/carrier evidence (AGENTS.md §5), not a decision; variant a gives WT 0 without R1251, so the reversal is load-bearing |
| The recommended design was tested | Refuted | Scout's scripts never combined a2 with R1218 reversible, no EGC with R1218, no sweep. The LP verifier then ran it (WT, EGC, full sweep above), so the numbers now exist but are the verifier's, not the scout's |
| a2 renormalisation constant 0.8286516853932584 | Corrected | Sum of the 14 fractions is 0.7512359550561798 (0.82865 × 4450 = 3687.5 is not a residue count); the scout's LP used the correct constant (reproduces to 4e-15) |
| c2 "converts the family into unresolved; must not be used" for that reason | Refuted | Since 2026-10-30 a maintenance-shortfall infeasible is certified and scored growth 0 (TP); the certificate was not run for c2. Reject c2 on its absent composition basis (0.001 mtprot per flux unit) instead |
| The 2026-10-26 "out of scope" sentence is a recorded decision being reversed | Corrected | It is a census classification in a computation-only entry; the 2026-10-29 entry listed the dead-end mitochondrial tRNAs as follow-up category d. Scoped authorization is still required |
| Christen & Sauer 2011 says Y. lipolytica respires fully | Corrected | The abstract names the species and the phrase but does not say which species fall in the fully-respiring class; use Kerscher 2001 ("obligate aerobic") |
| Variant a is "biologically faithful" | Corrected | mtDNA encodes 19 tRNA amino acids (no Val) incl. Ala/Cys/Gln/Gly/Pro/Ser that the model lacks; the model's Val module has no mtDNA tRNA; a2 drops 17.1 % of residues |
| "No UniProt/RefSeq annotation in the model" for Asn and Leu synthetases | Corrected (understated) | RefSeq XP_503562.2 (PRK03932 AsnRS, similar to SLM5) and XP_504899.3 (similar to NAM2) exist; Glu XP_504083.3 (similar to MSE1) |
| Sweep rule and model | Corrected | Rule A with 16 mtDNA genes and COX9 in scope (rule B baseline 136/190/46/692, 0 unresolved, not reported); SHA absent from body; newest candidate differs |
| Nuccore URL opened | Corrected | WebFetch got a browser check; the record was obtained via eutils efetch (strain W29, 47,916 bp); canonical cox1/cob CDS are flagged pseudo=true in the record |

**Evidence.**

| Claim | URL | Status |
|---|---|---|
| Kerscher 2001: obligate aerobe, 47.9 kb mtDNA, ND1-6/ND4L, COB, COX1-3, ATP6/8/9, incomplete tRNA set | https://pubmed.ncbi.nlm.nih.gov/18628906/ | opened (both) |
| AJ307410 record and CDS FASTA (SHA a88dab8a3a82ae242b26690f5061952676fab910dcad6177a8976a1758ec65af); 27 tRNA features, 19 amino acids | NCBI eutils efetch; `$FN29/mito_translation/downloads_20261102/AJ307410_cds_aa.fasta` | opened (efetch; the nuccore page itself not machine-readable) |
| Christen & Sauer 2011 | https://pubmed.ncbi.nlm.nih.gov/21205161/ | opened; attribution over-reaches |
| UniProt via RefSeq: Q6CBI1, Q6CAE1, Q6CBI6, Q6CCB9 (no location), Q6C2L5, Q6C109, Q6C814, Q6C474, Q6C0A5 | https://rest.uniprot.org/uniprotkb/search?query=xref:refseq-XP_501693.3&fields=accession,protein_name,gene_names,cc_subcellular_location&format=tsv and per-accession | opened (both); ARBA/PANTHER/IBA only |
| 12 model UniProt accessions Inactive | https://rest.uniprot.org/uniprotkb/A0A1D8NAJ0.json etc. | opened (both) |
| Yeast-GEM: 17 mitochondrial charging reactions, all dead ends, no translation term | yeast-GEM.txt main (biology lens) | opened |
| STATE.md: out-of-scope sentence; 2026-10-29 R760 → YALI1A09121g, category d; 2026-10-21 16 infeasible complex III/IV knockouts | `model/STATE.md` (local line 106; origin/main lines 632, 756, 767, 1475) | opened |
| Fraction of cell protein that is mt-encoded | none | no source found by anyone |
| S. cerevisiae ortholog names (MSF1 etc.) | https://www.yeastgenome.org/ | recalled, not checked |

**Evidence gaps.** No source for the mt-encoded protein fraction; no experimental localisation for any of the 14 synthetases (YALI1C15180g has no location at all — AlphaFold/TargetP is pre-authorized for identity questions); no mechanism for mitochondrial Met or Ile import; mtDNA tRNA set vs the model's 14 modules unreconciled (Val module without mtDNA tRNA; six residues without modules); why Ile and Trp synthetases and the adenylate kinase are absent from the list is unknown; c-scan and variant a sweep not re-run; one medium.

**The decision.** Whether to open a scoped round that adds mitochondrial translation at all, knowing that under the project rule it cannot proceed until (a) a sourced fraction of mt-encoded protein is found or the user names an unsourced coefficient as a labelled exception, as was done one by one for lipoylprotein 1e-5 and the provisional ceramide term, and (b) carrier evidence exists for R1251 and R1218. Sub-decisions: accept the 3 FPs (or resolve YALI1C15180g's localisation first); which candidate; rule-B re-score.

## 4. Recommended order of work

Prerequisites for every item: name the base candidate (origin/main carries 164542fb, c1f70ac4 and 2bc263f2 after e1e7c868; no reference model is designated), re-run the sweep under benchmark.md's primary rule (certified-infeasible as 0, mtDNA + COX9 out of scope) with the earlier rule alongside, use the real PO1f profile (plasmid LEU2 complement of R45), state the KO/WT ratio-move gate (the 2026-10-28 1e-6 gate would stop the N-glycan build), and put the model SHA in every report body. The combined effect of several families was not run; interactions (shared Dol-P, PI/PE/DAG pools, GDP-mannose) are unknown.

| Order | Item | Why this position | Flow |
|---|---|---|---|
| 1 | Cardiolipin term (2 FN → TP, +1 TP YALI1F08225g, +1 FP) | The only Y. lipolytica-sourced coefficient; one biomass term, no new reaction, no closure; fully reproduced; smallest growth cost (0.19 %). Blocked solely by the lipid-unlump ledger's recorded "blocked" biomass-weight decision, which the owner must lift | curation-change with its own TASK.md once the owner lifts the block and the user chooses the value (7.15e-4 vs 2.12e-5); guarded apply step on the cofactor template pinning the newest biomass_C |
| 2 | GPI term, variant b1 (6 FN → TP, +1 FP PBN1) | One term, no reversal (b1 keeps R_DM_GPI open), 0.044 % cost, fully reproduced; needs the lipid owner's rescaling of the S. cerevisiae proxy and the m603 identity items resolved first (cofactor precedent) | curation-change with its own TASK.md (curation/biomass data file; coordinate the biomass_C pin; consult the lipid colleague on PI/PE/DAG draw) |
| 3 | Sphingolipid: V1 first (8 FN → TP, +1 FP, no closure) or V2 (9 FN → TP) | V2 needs the most decisions (10 inherited closures, a gene-less DAG link, species merges that sit in the lipid line's unlump plan, species/compartment chosen by the owner); V1 has no full sweep yet. Both use a yeast-7.6 placeholder the owner must replace | curation-change with its own TASK.md only after the lipid owner sets species and value; V1 needs a full sweep first; YALI1E19647g identity evidence before R1825-R1834 are opened |
| 4 | N-glycan, variant B part only (1 FN → TP VRG4, +1 FP GDA1-like or 0 with YALI1E39570g) | The R289 mass-creation fix is evidence-backed and should be made regardless; the Golgi route, the 1:3 split and the R1137 closure are design decisions costing 9.9 % WT; the Yeast-GEM-faithful 1:1 alternative fixes the balance without flipping VRG4 | Two TASK.md candidates: (4a) R289 1:1 stoichiometry correction (chemistry); (4b) Golgi mannan route + R1137 closure (pathway), after BLAST of YALI1E39570g and the dolichol-formula clean-up |
| 5 | YALI1D08000g GPR correction | Not a demand question; delta(4)-desaturase mis-assigned to the cytosolic IPC synthase set; needs BLAST and a GlcCer branch the model lacks | Separate identity investigation (AlphaFold/BLAST pre-authorized for identity questions), then a GPR curation-change if confirmed; stays FN meanwhile |
| 6 | Mitochondrial translation (8 FN) | Mechanically complete and fully reproduced (134/181/49/697), but blocked on evidence: unsourced coefficient and two gene-less transport reversals with no mechanism; reverses a census classification | Document as structural false negatives (category d of the 2026-10-29 entry); reopen only when a mt-encoded protein fraction and Met/Ile carrier evidence are found, or when the user names an explicit exception |
| 7 | ALG7 / ALG1 (2 FN) | Blocked on composition evidence (no N-glycan content in any yeast source; Vega & Dominguez 1986 unread) | Document as structural false negatives; a literature task to obtain Vega & Dominguez 1986 is the only concrete next step |

Items 1-4 go through the curation-change flow with their own TASK.md (four to five candidates). Items 6 and 7 (10 genes) and YALI1D08000g (1 gene) stay documented as structural false negatives until the named evidence exists.

## 5. Verification levels

| Level | Sphingolipid | GPI | Cardiolipin | N-glycan | Mitochondrial translation |
|---|---|---|---|---|---|
| Reproduced by the LP lens (independent script, own LPs) | V2 rebuilt from design: WT, 16 family/FP-risk KO/WT, EGC (both), 4-cutoff before/after counts, 10 class changes, 20 unresolved (1,433 LPs; 763 non-lethal baseline classes taken from scout JSON) | Structure, WT, 9 baseline KOs, max R_DM_GPI, coefficient, variant a (WT, 7 infeasible KOs, EGC infeasible), b2 (WT, 7 lethal KOs context + fresh, EGC, full 1,081-gene sweep, 7 class changes; 886 baseline classes from scout JSON) (1,316 LPs) | Variant a: WT, route fluxes, 8 KOs, EGC, fully independent after-sweep, 205 baseline re-runs, 4 class changes; variant b infeasibility and EGC (1,342 LPs) | Baseline, A2, B, B0, C2 WT; family and FP-risk KOs; route FVA; dolichol leak; EGC + free-source probes; full C2 sweep and 196 baseline re-runs; 5 class changes; coefficient probes (corrected) (837 LPs) | Structure, WT, 15 blocked reactions, 8 baseline KOs, AJ307410 composition, variant a (WT, a_noMetFix, Met sink, EGC, 9 class-change KOs), the recommended variant (WT, EGC, full 1,081-gene sweep, 12 class changes), R1218 probe, respiration controls (1,176 LPs) |
| Spot-checked by the biology lens | 12 LPs (V2 WT, 10 KOs incl. pFBA route) | 8 LPs (baseline, b2 WT, 6 KOs) | 10 LPs (WT, variant a WT, 8 KOs) | 9 + 2 LPs (WT, C2 WT, 5 KOs, YALI1E39570g OR probe) | 44 LPs (recommended variant: WT, EGC, 8 family, Thr, 3 FP, FMT, Ala-T, ILV, 10 extra genes) |
| Scout only (not re-run) | V0, V1 sensitivity series (1e-5…1e-2), V1b, V3; excess route flux sink in 28 zero-growth KOs | b1, variant c sensitivity, pFBA co-product shifts, FVA at fraction 0/1, LEU2/R45 explanation | Range ends 7 %/9 %, alt (i)/(ii), c probe 1e-6…1e-1, 2026-10-30 alt-coefficient sweep | Variant A, variant C, A_gpr, B_R602_or probes; 1e-6 ratio-move count | Variant b, c2, a_nofmet, c-scan, full variant-a sweep counts (class-change genes re-derived individually) |
| Recalled → verified by a verifier | yeast 7.6 r_2108 origin of 0.000417 | Yeast-GEM has no GPI term; KEGG R12351/R05922 | — (Yeast-GEM claim was refuted, not verified) | Orlean 2012 O-mannosylation and OST statements | — |
| Still recalled or unread | Joshi 2024; Ejsing 2009 not opened; KEI1 ortholog beyond KEGG; YALI1→YALI0 by sequence; BLAST of YALI1D08000g | Vega & Dominguez 1986; third screen's raw PBN1 call; SMP3/GPI7 orthologs | Sekova methods on membrane/storage separation; PANTHER name for PTHR45848; Y. lipolytica PGS1/CRD1 phenotype literature | Orlean "30-50 % mannoprotein"; Vega & Dominguez 1986; R1137 provenance in iYli21 | SGD names; mt-encoded protein fraction (no source); R1251/R1218 carriers; YALI1C15180g localisation |
| Rule-B (project-comparable) counts | inferred 145/181/47/691 (not run) | not inferred (needs ko_screen_certified.py) | inferred 139/187/47/691 (not run) | inferred 139/187/48/690 (not run) | not inferred |
| Combined-family model | not run by anyone | not run | not run | not run | not run |

Source files: `$FN29/sphingolipid/` (01-08 scripts, `06_sweep_*.json`, `verify_lp/verify2*.py|json|log`, `verify2b.*`, `inspect2.*`); `$FN29/gpi/` (`s1`-`s7`, `s3_sweep_*.json`, `verify_lp/repro_gpi_lp.*`, `inspect_structure.py`); `$FN29/cardiolipin/` (`10_verify`, `11_variants`, `sweep_*_20261102.*`, `12_compare_a_sekova.json`, `verify_lp/verify_cardiolipin.*`, `extra_checks.*`); `$FN29/nglycan/` (`t2_s1`-`t2_s5`, `t2_s4_sweep_*.json`, `t2_summary.json`, `verify_lp/verify_c2.*`, `probe_variantA_cn*.py|json`); `$FN29/mito_translation/run_20261102/` (`SUMMARY_20261102.json`, `s1`-`s5`, `sweep_baseline.json`, `sweep_a.json`, `mt_composition.json`), `$FN29/mito_translation/verify_lp/` (`verify2.*`, `sweep_recommended.json`, `extras.*`) and `$IYALI26_RESEARCH_ROOT/artifacts/tasks/fn29_demand_families_20261102/verify_mt/` (`spot.py`, `mt_inspect.py`, `dl/`). Repository references: `model/STATE.md`, `model/expected/benchmark.md`, `model/reports/gap_metabolites_ponce2013_20261026/`, `model/reports/closed_bounds_69_20260909/REPORT.md`, `model/reports/gpi_chain_links_20261021/REPORT.md`, `model/reports/infeasible_ko_rule_20261030/scripts/ko_screen_certified.py`, `model/curation/lipid_unlump/lipid_unlump_sn_core_curation.json`, `model/curation/gapfill/gpi_chain_links.json`, `model/reports/cfsa_lipid_unlump_20260908/PREFLIGHT.md`.

## Reader notes (added 2026-11-02; marked, not corrected)

1. **R_NTP1 exists in `e1e7c868`.**
   - Section 3.3 ("Refuted or corrected", row "R_NTP1") says no such reaction exists in `e1e7c868` and that it is an
     in-memory reaction added for the check.
   - The archived model (`$FN29/e1e7c868.xml`, SHA-256 checked 2026-11-02) contains `R_NTP1`
     (`m141[C_cy] + m32[C_cy] --> m143[C_cy] + m35[C_cy]`, GPR `YALI1A09766g or YALI1E41893g`, bounds (0, 1000)).
     Verified this session by reading the model.
   - The sphingolipid, GPI, cardiolipin and mitochondrial-translation scripts maximise this model reaction by name in
     their energy checks (read in the archived scripts; the N-glycan script was not checked).
2. **New FP count for the 18 decidable genes.**
   - Sections 1 and 2 ("Totals") give five new FP (four with `YALI1E39570g`) and include the CWH8-like `YALI1D26639g`.
   - In section 3.4's own table, `YALI1D26639g` stays at KO/WT 1.0 under variant B. It becomes 0.0 only with A2. The
     scout's `nglycan/t2_summary.json` also lists only `YALI1C25093g` as the B-variant FP risk.
   - Under that table, the 18-gene build (VRG4 through B alone) has four new FP (three with `YALI1E39570g`). The CWH8 FP
     belongs to the ALG7/ALG1 part.
   - The 28-gene figure of eight (or seven) is unaffected.
   - Inferred from the report's tables and the scout JSON. The LP lens did not knock out `YALI1D26639g` under B.
3. **LP counts.** Section 5's LP-lens counts include auxiliary scripts: sphingolipid `verify2b` 16, cardiolipin
   `extra_checks` 22, N-glycan probes 7, mitochondrial `extras` 13. So they exceed the `lp_count` in each main
   verifier JSON. The sums match exactly (TASK.md).
4. **Escher figure.** `escher/` reproduces the WT of every recommended variant in this report to within 1e-9 (581
   LPs). It draws the family reactions as built and under each variant applied alone.
