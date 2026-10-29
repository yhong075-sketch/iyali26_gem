# Independent audit: UM11 nuclear RNR (Option B) and the combined baseline 14f73bda (2026-10-29)

- **Auditor:** read-only subagent, no conversation history beyond the brief. Wrote only this file. No model
  built, no git state changed, no repository file edited.
- **Model loaded:** short-path copy of
  `model/candidates/E5_…_um7_nuclear_diffusion_r2171_1to1_20261029.xml`, SHA-256
  `14f73bda4b0772019e82ea18a27e198c5ba175e67fef205a51dbfb090809e757` (checked before loading). Parents `5f54eb22…`,
  `965ab947…` and grandparent `9d655c2b…` were also copied and SHA-checked.
- **Condition:** Task A medium `6eca70fd…` + PO1f profile `35307853…` via `load_effective_simulation_context`;
  "open" = raw model, every `model.exchanges` reaction at (−1000, 1000). Gurobi 13.0.2, cobra 0.30.0,
  tolerances 1e-9 (`checked_solve.configure`), real HOME.
- **LPs:** 10,807 in total (own LP checks 7,230; knockout spot-checks and RNR controls 70; Option A′ scan 3,481;
  figure numbers 24; parents' WT 2).
- **Own scripts and outputs** (scratch, not kept in the repository):
  `/private/tmp/claude-501/-Users-david-Desktop-Lab-Ian-wheeldon-code-iyali26-gem--claude-worktrees-wizardly-wing-a0d6c0/3bda569d-2d69-4a78-a510-ddf05308ab20/scratchpad/audit/`
  (`field_diff.py`, `audit_lp.py`, `audit_ko.py`, `audit_aprime.py`, `audit_fig.py`, `sweep_tables.py`,
  `ancestor_static.py`, with `.log` / `.json` outputs; fetched sources in `src/`).
  Option A and A′ were written by the auditor, not imported from `um11_apply.py`.

## Claims

| # | Claim | Auditor's check and value | Verdict | Note |
| --- | --- | --- | --- | --- |
| 1a | Candidate SHA-256 is `14f73bda…` | `shasum -a 256` on the file and the copy: `14f73bda4b07…e757` | supported | |
| 1b | Diff vs `5f54eb22` is only R2171; vs `965ab947` only the UM7 change; each field-identical to that parent's own diff vs `9d655c2b` | (i) raw `diff` of the XML files: 5f54eb22→14f73bda = R2171 notes block + 3 stoichiometry lines (100→1 H2O, 100→2 H+, 100→1 Ins1P); 965ab947→14f73bda = R_ADPRtn + R_PPItn + 2 group members. Hunk content (line numbers ignored) is identical to 9d655c2b→965ab947 and 9d655c2b→5f54eb22 respectively. (ii) `field_diff.py` (cobra, reactions / metabolites / genes / groups / compartments / objective / notes / annotation): 2 entries (R2171 stoich, notes) and 3 entries (R_ADPRtn added, R_PPItn added, group101 members); counts 2,319 / 1,882 / 1,081 / 113 groups | supported | |
| 1c | STATE combined entry: built at HEAD `ff5192d` (= `origin/main`); argv = `5f54eb22` argv + `--r2171-pi-plc-1to1`; seeds 0 and 1 identical; WT as in both parents; path ≈ 300 characters | `.build.json`: head `ff5192d6…`, input/source unchanged = True; argv diff vs the `5f54eb22` build.json is exactly `+--r2171-pi-plc-1to1` and the output path. `git rev-parse origin/main` = `ff5192d6…` (no fetch). Parents' WT: 0.3116608216336134 (5f54eb22), 0.31166082163360886 (965ab947). Absolute path 307 characters. Both build logs carry the LibSBML −3 line, as do the parents' logs | supported, except seed identity | Seed-1 byte identity rests on `logs/sha256.txt` and `build_seed1.log` (output "seed1.xml" not retained); not reproduced (no builds allowed) |
| 2 | No model file was written for Option B | `git status --porcelain --untracked-files=all --ignored`: tracked XMLs unmodified; the only untracked XML is the combined baseline (+ its `.build.json`, `.coq9_genes.tsv`); no ignored XML; no XML written under the research root since 10:00 other than literature sources | supported | |
| 3 | WT 0.3116608216; all 11 UM11 reactions blocked under Task A and open | `audit_lp.py`: WT optimal 0.3116608216336117, check passed; R1354 (−1000, 1000), R612 (0, 0). Own min/max loop: 11/11 blocked (|v| ≤ 1e-9) under Task A and under open | supported | Also provable: `m931[C_nu]` appears only as a substrate of R707/R709/R711/R1914, and the nuclear species of R782 and of the six transports have only UM11 partners (static listing, row 4), so the block holds for any medium |
| 4 | Roots: m931[C_nu] no producer, m933[C_nu] no consumer; no nuclear TRR, NADPH/NADP, thioredoxin transport or glutaredoxin; nuclear copies have the GPR string of R706/R708/R710/R1913/R781 | Static (cobra): m931 consumed by R707/R709/R711/R1914 only; m933 produced by the same four only. Of 47 nuclear species none has "NADP" in its name (NAD and deamido-NAD exist). No species named glutaredoxin/Grx anywhere. GPR strings identical for all 5 pairs; metabolite names and coefficients identical per pair except R1914 water in `C_cy` (`m32[C_cy]`). Nuclear partners: CDP R709/R1033, GDP R711/R1136, UMP R782/R1353, UDP R782/R1914, dADP R707/R1076, dCDP R709/R1077, dGDP R711/R1093 | supported | Glutaredoxin absence is by name search only. R1914 carries SBO:0000185 (transport), a side effect of the water slip |
| 5a | Option A (two reversible gene-less links) frees exactly the 11; WT unchanged; each link alone frees none | Own links `m632⇌m931`, `m633⇌m933`. UM11 blocked: red only 11/11 and 11/11; ox only 11/11 and 11/11; A 0/11 and 0/11 (Task A / open). Whole-model blocked (explicit max/min per reaction, |v| ≤ 1e-9): as built 828 Task A / 695 open; A 817 / 684; freed = exactly the 11 UM11 ids in both conditions; none newly blocked. Sets identical member by member to `blocked_sets.json`. WT under A 0.31166082163361075 (Δ −9e-16). Links carry up to ±1.625 under A | supported | Two reactions (R970, R1343/R1344) have capacity 3.7e-7, between 1e-9 and 1e-6; they count as unblocked in both runs and do not affect the comparison |
| 5b | No energy-generating cycle in any variant | Exchanges (0, 1000), every positive lower bound relaxed (only `xMAINTENANCE` 7.8625), max R_NTP1 and max biomass_C: 0, status optimal, for as built, red, ox and A. Repeated with all 185 boundary reactions (incl. R1373, R2290, R_DM_GPI) at (0, 1000) for those four and for A′: same | supported | |
| 5c | pFBA: nuclear copies carry 0 at the optimum under A; only the R702/R706 split moves | cobra `pfba` under A: all 11 UM11 = 0, both links 0; R702 0.002895 + R706 0.000243 = 0.003138 (= as-built R702) | supported | One optimum; the "never helps growth" conclusion rests on the unchanged WT, which is verified |
| 5d | A′ (nuclear TRR copy) frees 11/11; census says "frees 11/11 and nothing else" | Own A′ (NADPH/NADP nuclear species, two diffusions, R755 copy with the same GPR, H+ `m627[C_nu]`): WT 0.31166082163361164; EGC 0/0 optimal; whole-model 817 Task A / 684 open, freed exactly the 11, none newly blocked | supported now | TASK.md section 5 lists A′ whole-model as "(not scanned)", so the census sentence was unverified when written; this audit supplies the scan |
| 6a | Knockout sweeps of A, A+, A′ vs baseline: no class change, max \|ΔKO\| 1.6e-13, no stale point, no fresh re-solve, 20 infeasible, same in every run | `sweep_tables.py` on the three saved TSVs (research root; the repo copy of `ko_trx+trxss.tsv` is byte-identical, `572fb7f9…`): 1,081 genes each; 0 class changes over 4 cutoffs; max \|ΔKO\| 1.625e-13 / 1.641e-13 / 1.644e-13 (gene YALI1D19108g); 0 points failing the check; same 20 infeasible genes in all 6 columns; baseline columns identical across the three sweeps | supported | |
| 6b | Spot-check of the sweeps | `audit_ko.py`: 28 genes (RNR, TRX, TRR, URA6, 3 infeasible, 10 near the cutoffs, 8 random) on baseline and on own Option A. Every value equals the saved TSV bit for bit (max \|own − saved\| = 0.0); max \|baseline − A\| on these 28 = 1.6e-13 | supported | Partially covers the untested "as_built reproduces the baseline column" in REPORT section 6 |
| 6c | RNR controls (`rnr_controls.log`) | Baseline and A: triple KO (D04261g, F39330g, B05072g) ratio 1 (closes R706–R711, R1913, R1914); R702–R705 closed ratio 1; each single KO + R702–R705 closed ratio 1; triple KO + R702–R705 closed growth 0.0, optimal, check passed | supported | |
| 7 | Benchmark table (REPORT section 4) | Rescored from the TSV with `benchmark.md` rules (strictly below cutoff; 16 mtDNA genes excluded; negatives = model genes not on the list): 1 % 94/221/35/699, 5 % 100/215/37/697, 10 % 104/211/42/692, 15 % 112/203/42/692; 16 unresolved, 11 of them positives; 326 in-model positives, none mtDNA; 1,286 out-of-scope positives; recall/specificity/precision match. TSV class columns: 0 mismatches. Same in all three sweeps | supported (counts); wording partial | "4 mtDNA genes are out of scope" (REPORT section 4, TASK section 5) is wrong: all 16 mtDNA genes are in the model and excluded (1,081 − 16 = 1,065 = 1,049 scored + 16 unresolved). The 4 are the mtDNA genes among the 20 infeasible knockouts. 15 % equals the STATE entries of both parents (112/203/42/692) |
| 8a | Figure numbers (`escher/README.md`) vs `flux.json` / `index.json` | 25 mapped reactions; UM11 blocked 11 → 0; capacity > 0: 12 → 25; pFBA ≠ 0: {R702, R710, R1913, R755, R781} → + R706; no-uptake ≠ 0: {R789, R1104} both; whole-model 828 → 817, freed = the 11, newly blocked none; 9 dead-as-built nuclear species (m931, m932, m933, m934, m935, m936, m937, m1002, m1003), 28 alive; LP count 9,476. Own LPs: capacity R706 [0, 1.578868], R1913 [0, 1.625124], R781 [0, 55.47091], R755 [0, 13.684825], R789/R1104 [−1000, 1000], each equal to `flux.json`; no-uptake R789, R1104, R1105, R788 [−1000, 1000], R706/R707 [0, 0] | supported | |
| 8b | Viewer is self-contained and draws | No `<script src>` and no `<link>` in the HTML; Escher 1.8.1 inlined. Re-ran `check_viewer.js` in jsdom into scratch: 5 scenarios × 25 reactions, 156 metabolite circles, output identical to `logs/check_viewer_jsdom.log` | supported, with a caveat | The inlined Escher bundle contains a runtime `fetch("https://cdn.jsdelivr.net/npm/gif.js@0.2.0/…")` on its GIF-export path; it is not a load-time dependency. Headless Chrome check not re-run |
| 9a | Yao 2003 (PMID 12732713): Rnr1/Rnr3 mainly cytoplasmic, Rnr2/Rnr4 mainly nuclear, move to the cytoplasm under stress | efetch abstract: "predominantly localized to the cytoplasm" / "predominantly present in the nucleus"; redistribution under genotoxic stress | supported | |
| 9b | Lee & Elledge 2006 (16452505): Wtm1 anchors Rnr2/Rnr4 in the nucleus; on damage they join cytoplasmic Rnr1 | Abstract: Wtm1 "acts as a nuclear anchor"; on damage Rnr2/Rnr4 "become cytoplasmic, where they join with Rnr1" | supported | |
| 9c | Camier 2007 (17349928): thioredoxins are the major in vivo RNR reductants | Abstract: "establish thioredoxins as the major RNR reducing system in yeast" | supported; TASK partial | TASK section 4 adds "with glutaredoxin as backup"; the abstract says only that either pathway is essential. Not stated as a backup role for RNR |
| 9d | SGD TRR1 has no nucleus annotation | `backend/locus/S000002761/go_details` (55 rows): cytosol (IBA, RCA, IDA), cytoplasm (IEA), mitochondrial IMS (IEA, IDA), mitochondrion (HDA); no nucleus term | supported | |
| 9e | SGD: RNR1–4 "is active in cytoplasm", none "is active in nucleus" | go_details for S000000872/3563/1328/3412: each has "is active in cytoplasm" (RCA); none "is active in nucleus". "located in nucleus": RNR1 HDA, RNR2 IDA, RNR4 IDA/HDA | supported | |
| 9f | Dif1 (Lee 2008, 18851834; Wu & Huang 2008, 18838542), Clb6–Cdc28 (Wu 2018, 29263158), low iron (Sanvisens 2016, 26970775) | Abstracts state each point as cited | supported | |
| 9g | Mammals: RNR cytosolic and dNTPs diffuse into the nucleus (Pontarin 2008, 18997010); no nuclear reductase (Leeds 1985, 3915777) | Pontarin: RNR "resides in the cytosol", diffusion only "suggesting"; Leeds: "no reductase in nuclear preparations" (CHO) | partial | Diffusion of dNTPs is the authors' inference, stated as fact in REPORT section 5 |
| 9h | Thioredoxin: Koc 2006 (16574642), Pedrajas 1999 (10037727), nuclear only by high-throughput GFP (SGD HDA; UniProt citing Huh 2003), Timney 2016 size argument | Koc: thioredoxin a physiologically relevant RNR donor. Pedrajas: Trx1/Trx2/Trr1 "located in the cytoplasm" (background sentence). UniProt P22217 and P22803: Nucleus, ECO:0000269, PubMed 14562095. SGD: TRX1 nucleus HDA + IEA; TRX2 nucleus IEA only. Timney full text (saved Europe PMC XML in the research root, not re-fetched): 27-kD GFP equilibrated with τ = 15 ± 1.4 s | supported, minor nuance | For TRX2 the SGD nucleus call is computational (IEA), not HDA. Wente & Rout 2010 abstract supports free diffusion of small molecules only generically |
| 9i | Niida 2010 (20159953): human RNR recruited to damage sites | Abstract: RRM1/RRM2 accumulate at damage sites, Tip60-dependent | supported | |
| 9j | No *Y. lipolytica* RNR or dNTP literature (PubMed 0 hits) | Own esearch: "Yarrowia lipolytica AND ribonucleotide reductase" 0; "Yarrowia AND (RNR1 OR RNR2 OR "ribonucleotide reductase")" 0; "Yarrowia lipolytica AND (deoxyribonucleotide OR dNTP)" 0 | supported | The subagent's saved esearch has only the dNTP query; the RNR query is supported by this audit |
| 9k | Gene identities (RBH) and UniProt status | BLAST tables: F39330g RBH Rnr1 67.5 %; B05072g RBH Rnr2 66.2 %, Rnr4's best W29 hit is also B05072g; D04261g best hit Rnr3 (56.0 %), not reciprocal; D04261g vs F39330g 57.4 %; F02309g RBH Trx2 57.3 %; F12969g RBH Ura6 (KCY) 65 %; D35735g RBH Trr1 76.8 %, Trr2's best W29 hit is also D35735g. UniProt REST: Q6C7L4 TRR1 Swiss-Prot reviewed (FAD by similarity); Q6BZQ3, Q6CFU6, Q6CAF2, Q6C399, Q6C2B0 unreviewed | supported (tables read, BLAST not re-run) | Homology only, as the report says |
| 10 | Ancestor models: all carry nuclear RNR copies, none has a nuclear source of reduced thioredoxin; Yeast-GEM `r_4612` is a confidence-0 gap-fill | `ancestor_static.py` (cobra, producer/consumer by sign and bounds) on iyli21, iYali4, iYali5, Yeast-GEM, iMM904: nuclear reduced thioredoxin has no producer in any; Yeast-GEM `r_4612` reversible, notes "Confidence Level: 0" and "added after Gapfilling … unknown transporter" | supported | `model_scan.txt` prints r_4612 as "->", but it is reversible; REPORT text unaffected |

## Problems found (most serious first)

1. **Evidence tier, census report.** The new UM11 paragraph in
   `gap_metabolites_ponce2013_20261026/REPORT.md` says every opening, including the nuclear TRR copy (A′), "frees
   11/11 and nothing else". TASK.md section 5 records A′ as "(not scanned)" for the whole model, so this was not
   verified when written. This audit's scan confirms it (817 Task A / 684 open, exactly the 11 freed, none newly
   blocked). Suggest citing AUDIT.md for the A′ part.
2. **mtDNA wording.** REPORT section 4 ("4 mtDNA genes are out of scope") and TASK section 5 say 4; the benchmark
   excludes all 16 mtDNA genes, and all 16 are in the model. The 4 are the mtDNA genes among the 20 infeasible
   knockouts (20 − 4 = 16 unresolved). The score counts are correct.
3. **Biology wording slightly above the sources.**
   - TASK section 4: "with glutaredoxin as backup (Camier 2007)" is not in the abstract.
   - REPORT section 5: "dNTPs diffuse into the nucleus" is Pontarin's suggestion, not a measurement.
   - TASK section 4 tags "The nuclear R2 pool is sequestered and inactive" as **verified** (abstracts), repeated in
     the census paragraph. Sequestration is stated in the abstracts. "Inactive" is the authors' model: the active
     holoenzyme forms in the cytoplasm, and Sanvisens 2016 calls Dif1-driven import inhibitory. No yeast study cited
     measures the absence of nuclear activity, so this part is **inferred**.
   - "Rnr1 is cytoplasmic": Yao 2003 says "predominantly", and SGD has a high-throughput nucleus call for RNR1.
   - None of these changes the decision.
4. **TRX nuclear calls.** "Nuclear only by high-throughput GFP" holds for TRX1 (SGD HDA; UniProt citing Huh 2003).
   For TRX2, SGD's nucleus call is IEA (computational); only UniProt cites Huh 2003 for it.
5. **Omitted caveat.** `stage1/evidence.md` notes that the W29 RefSeq thioredoxin (F02309g) carries a 25-residue
   Arg/Lys-rich N-terminal extension resembling a mitochondrial presequence (inferred), and TRR1 a 115-residue
   one. REPORT section 5 does not mention this. It is relevant to any future thioredoxin-compartment option (A),
   not to B.
6. **Viewer offline claim.** The inlined Escher bundle has a runtime fetch to cdn.jsdelivr.net (gif.js worker)
   on the GIF-export path. Drawing needs no network (jsdom re-run identical), but "no internet needed" is not
   true for that export function.
7. **Not reproduced by this audit:**
   - the seed-1 rebuild of `14f73bda` (no builds allowed);
   - A+ and C1 WT values (taken from `stage1_um11.log`; C1 neutrality is provable);
   - the headless-Chrome viewer check;
   - BLAST runs (tables read only).
8. **Minor inconsistencies:**
   - TASK section 4 lists Yao, Lee, Camier and Pontarin as read by the main session; REPORT section 5 adds Leeds
     1985.
   - `model_scan.txt` shows Yeast-GEM `r_4612` as irreversible ("->"), but it is reversible.

No number in REPORT.md, the census paragraph, the escher README or the audited STATE entry failed to reproduce,
apart from the mtDNA count wording in item 2.
