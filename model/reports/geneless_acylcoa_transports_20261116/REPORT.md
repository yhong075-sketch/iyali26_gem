# Gene-less CoA-ester transports across the inner mitochondrial membrane in 7ccbb568 (2026-11-16)

- **Scope and authorization:** [TASK.md](TASK.md).
- **What was done:** computation and documentation only. No model XML, curated data, builder code, benchmark or test
  was changed. Committed at the user's instruction ("Commit it") on branch `claude/lucid-visvesvaraya-a40bbb`; not
  pushed.
- **Where things are:** tables in [tables/](tables/), scripts in [scripts/](scripts/), figures in [figure/](figure/).
  Raw outputs, logs, the literature responses and the TargetP files are in
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/geneless_acylcoa_transports_20261116/`.

## Answer in short

**What the model does today, and whether it needs these reactions.**
- **R2025 and R2113 carry one route.** They pass propionyl-CoA through the matrix: cytosolic m1970 goes in on R2113
  and comes out on R2025 as m1875, which has no formula. The route then ends as propanoate: R2024 runs backwards and
  makes ATP, and the propanoate is exported. This is the only outlet of the cytosolic isoleucine branch
  R2281–R2285.
- **The other six are blocked**, in both conditions and on both solvers: R2077, R2083, R2084, R2118, R1291 and R2117.
- **None of the eight is needed** (observed, LP):
  - closing any of them, or R2025 and R2113 together, changes WT growth by at most 1.6e-15;
  - at 99.9999 % of the optimum, R2025 and R2113 can be 0;
  - none of 30 targeted gene knockouts changes, as is or with R2025, R2113 or both closed.

**Carriers (literature).**
- **No carrier is known for any of the six species**, in yeast, mammals or plants. This is a search result (two
  targeted Europe PMC searches, 77 hits): absence from a search, not proof of absence.
- **What the known carriers move** (verified, abstracts and records):
  - free CoA and dephospho-CoA: yeast Leu5, human SLC25A42/A16, fly dPCoAC;
  - carnitine and acylcarnitines: yeast Crc1, which moves acyl groups and not the CoA moiety.
- **No transporter is known for 4'-phosphopantetheine** (verified, full text). Yeast PPAT (Cab4) is in the cytoplasm
  and nucleus; DPCK (Cab5) is on the outer mitochondrial membrane.

**Valine steps (R2114–R2117, with R2112).**
- **Every homologue outside fungi is mitochondrial**, in human and plant (verified, records). The one exception
  is Arabidopsis CHY1, which is peroxisomal.
- **The *Y. lipolytica* proteins carry mitochondrial presequences** (TargetP-2.0, a prediction):
  - YALI1D19252g (R2117), mTP 0.999;
  - YALI1F03883g (R2114), 0.98;
  - YALI1C02550g (R2112), 0.9995;
  - YALI1B14044g (R2116), 0.9998, but only on its 319-amino-acid form.
- **No compartment data exist for *Y. lipolytica***, and no source puts these steps in the cytosol in any organism.

**Options tested in memory** (Task A + PO1f, Gurobi 1e-9 with a GLPK check).
- **Results common to all 11 variants:**
  - WT growth stays 0.31110053592595 (|Δ| ≤ 1e-14);
  - the energy check is 0 on both solvers;
  - none of the 30 targeted knockouts changes class at any cutoff, in the 6 variants where they ran.
- **The leading package (option "ALL")** frees the valine branch and the methylcitrate trio (R1969, R25, R95, R552,
  up to 5.86) inside the matrix, with no CoA ester crossing the membrane. It combines:
  - merging m1875 into m1970;
  - no propionyl-CoA transport;
  - the valine steps moved into the matrix;
  - R2077, R2083, R2084, R1291 and R2118 closed.

  R834 (the CoA carrier) stays blocked, as the R1172 round predicted from the matrix CoA identity.

**Full certified sweep of the leading option** (B0/ALL, all 1,081 genes; follow-up of 2026-11-16, section 11).
- **No gene changes class** at 1, 5, 10 or 15 %, under Rule B or Rule A. TP, FN, FP and TN are identical to the
  as-is sweep at every cutoff. At 15 %, Rule B: TP 145, FN 181, FP 43, TN 694; recall 0.445, specificity 0.942,
  precision 0.771.
- **WT is unchanged** (+1.2e-15). The same 20 knockouts are certified maintenance-shortfall infeasible.
- **Only the 12 ATP-synthase subunit genes (R171) change KO/WT:** 0.10778 → 0.10853, still between 10 % and 15 %.
  With ATP synthase gone, the option lets valine be degraded in the matrix: R1969, the moved R2117–R2112, then the
  methylcitrate cycle, with FADH2 going to R1977. This adds 0.7 % growth.

### Per reaction

| Reaction | As written | Today: capacity on Task A (open exchanges) | Needed now? | Carrier evidence | Suggested option (user decides) |
| --- | --- | --- | --- | --- | --- |
| R2025 | propionyl-CoA m1875 (cy, **no formula**) ⇌ m78 (mi) | export only, ≤ 2.67 (1000) | no | none found for propionyl-CoA | **merge m1875 into m1970**; R2025 then duplicates R2113 and goes |
| R2113 | m78 (mi) ⇌ m1970 (cy) | import only, ≤ 2.67 (1000) | no | none found | after the merge it is blocked as is. With O1 (other line) it is the exit of side route B. **Close.** |
| R2077 | succinyl-CoA (mi) ⇌ succinyl-CoA (cy) + 5 H+ | blocked | no | none found | **close**; recommend only, since it is in the colleague-owned lipid-unlump handoff |
| R2083 | 3-methylcrotonyl-CoA (mi) ⇌ (cy) | blocked | no | none found | **close**, with R2084. R2082 then has no use. Its gene YALI1E15555g (mTP 0.9998) is a candidate rule for the matrix step R1971 |
| R2084 | isovaleryl-CoA (mi) ⇌ (cy) + 4 H+ | blocked | no | none found | **close** (as R2083) |
| R2118 | isobutyryl-CoA (mi) ⇌ (cy) + 4 H+ | blocked | no | none found | **close, and move the valine steps into the matrix** |
| R1291 | pantetheine 4'-P (cy) ⇌ (mi) | blocked | no | no transporter known (Xue 2023) | **close**. The matrix copy R1292 has no gene, and yeast PPAT is not in the matrix |
| R2117 | isobutyryl-CoA (cy) + FAD (cy) ⇌ methacrylyl-CoA (cy) + FADH2 (**mi**) | blocked | no | — | **move to the matrix** (FAD and FADH2 both in the matrix), with R2116, R2115, R2114 and R2112 |

## 1. Run

| Item | Value |
| --- | --- |
| Model | `model/candidates/E5_353cefdc_xref_20261115.xml`, SHA-256 `7ccbb568…2146`, checked by every script before use |
| Conditions | "open": every exchange at (−1000, 1000). "Task A": medium `6eca70fd…` + PO1f profile `35307853…` through `load_effective_simulation_context` |
| Solvers | cobra 0.30.0. Gurobi 13.0.2 at feasibility and optimality 1e-9, basis discarded per LP. GLPK with `glp_adv_basis` per LP |
| Blocked | \|min\| and \|max\| ≤ 1e-9 |
| LPs | 1,086 of the 1,300 approved: step 1, 360; step 4, 722; step 5, 4. Every LP was optimal; no re-solve was needed. The follow-up sweep (step 6, its own budget of 1,300): 1,141 LPs (section 11) |
| Literature | 27 read-only fetches by a read-only reviewer subagent, of the 40 approved (Europe PMC REST, UniProt REST, KEGG REST, SGD). Responses are saved as `lit/f01–f27` with `MANIFEST.tsv` |
| TargetP-2.0 | Official DTU server, Non-Plant, short output. Job `6AFBA20D0032281B2E2E5862` failed: the FASTA was sent as a text field and the server answered "No input file provided". Job `6AFBA2A60032297C72156452` predicted 8 of 8. Job `6AFBA4A800322DDEA260BFEA` ran the 319-aa form of YALI1B14044g. Files in `targetp/` with `jobs.tsv` and `MANIFEST.sha256` |
| Sequences | W29 RefSeq proteins from `w29.faa` (`1863842c…`). Every gene is `ok` in `w29_tag_protein_audit.csv` |

## 2. The eight as written (static, verified on the model this session)

[tables/capacity_step1.csv](tables/capacity_step1.csv) lists the equation, rule, bounds and balance of the eight and
of 32 neighbours.

- **Two cytosolic propionyl-CoA species.**
  - m1875 "propionyl-CoA" has no formula, no charge and no annotation. Only R2024 and R2025 use it.
  - m1970 "Propanoyl-CoA" (C24H36N7O17P3S, −4; KEGG C00100) is made by R2285 and R2112 and used by R2142 and R2113.

  So R2113 + R2025 together convert m1970 into m1875 through the matrix.
- **Protonation forms.** R2077, R2084 and R2118 carry 5, 4 and 4 H+ on the cytosolic side, because the matrix ester
  is written neutral and the cytosolic one charged. This was already recorded on 2026-09-11
  (`model/reports/model_static_screen_20260911/COMPARTMENT.md`).
- **Leucine.**
  - In the matrix: R1968 (BCKDH) → R1971 (isovaleryl-CoA DH) → R1972 (MCC).
  - R2084 → R2082 → R2083 is a cytosolic copy of R1971 that uses cytosolic FAD.
  - R1971's gene `YALI0D12573g` has no YALI1 match in the W29 audit.
- **Valine.**
  - The matrix BCKDH R1969 makes isobutyryl-CoA. Its only consumer is the export R2118.
  - The rest of the chain is cytosolic: R2117, R2116, R2115 (no gene), R2114, then R2112 or R2074.
  - R2117 takes FAD in the cytosol and returns FADH2 in the matrix.
  - m1972, m1974 and m1975 have no formula.
- **CoA synthesis.**
  - The full cytosolic route exists: R615 PPAT, `YALI1D24649g`; R1092 DPCK, `YALI1F13202g`.
  - There is also a matrix copy: R1291 import → R1292 PPAT (no gene) → R265 DPCK (the same gene as R1092).
- **Gene-ID duplicates** (W29 audit), which limit any knockout claim on these genes:
  - `YALI0B10406g` (on R2283) is `YALI1B14044g` (on R2116, R2141, R1978 and peroxisomal isomerases);
  - `YALI0E18568g` (on R2285) is `YALI1E22238g` (19 peroxisomal thiolase reactions).
- **Colleague-owned.** R2077 is listed under `repaired_reaction_ids` in
  `model/curation/lipid_unlump/fatty_acyl_coa_handoff.json`.

## 3. Capacity and need (step 1; observed, LP)

**Capacity.**

| Reaction | Open, Gurobi / GLPK | Task A, Gurobi / GLPK | Task A, growth ≥ 99.9999 % | Task A, growth ≥ 99.9 % |
| --- | --- | --- | --- | --- |
| R2025 | [−1000, 0] / same | [−2.6654, 0] / same | [−3.04e-6, 0] | [−3.04e-3, 0] |
| R2113 | [−1000, 0] / same | [−2.6654, 0] / same | [−3.04e-6, 0] | [−3.04e-3, 0] |
| R2077, R2083, R2084, R2118, R1291, R2117 | 0 / 0 (GLPK ≤ 5e-18) | 0 / 0 | 0 | 0 |

Read the signs against the written direction:
- R2025's negative flux is m78 (mi) → m1875 (cy), so it can only export.
- R2113's negative flux is m1970 (cy) → m78 (mi), so it can only import.

On Task A, R2281–R2285 and R2022–R2024 have the same 2.665 ceiling. R968, the 3-methyl-2-oxopentanoate transport,
carries the same flux, so the isoleucine branch and the pass-through form one route. The forced route in step 5
shows it.

**Need.**
- **At the WT optimum:** pFBA gives 0 on all eight.
- **Closures:** WT changes by at most 1.6e-15, for each of the eight alone and for R2025 + R2113.
- **Targeted knockouts:** 30 genes, in [tables/targeted_knockouts.csv](tables/targeted_knockouts.csv).
  - These are the genes of reactions that share a crossing species, and of the CoA-synthesis, methylcitrate and
    branched-chain steps.
  - All 120 knockout LPs were optimal.
  - KO/WT is identical, to |Δ| ≤ 1.6e-14, with R2025, R2113 or both closed.
  - Essential at 15 %:
    - `YALI1D24649g` (PPAT), on the list: TP;
    - `YALI1D26360g` (the shared E3 of PDH and BCKDH), on the list: TP;
    - `YALI1D33780g` (pantothenate kinase), on the list: TP;
    - `YALI1F13202g` (DPCK): FP.
  - None of the four depends on the eight.
- **Why the six blocked ones were not re-tested with knockouts (provable).** A knockout only removes feasible points.
  A reaction that is 0 at every feasible point as is stays 0 under any knockout. Closing it therefore cannot change
  a knockout result.

**Genes and benchmark labels.** Rule B; negatives are the user-defined class (model genes not on the consensus
list). The list was used during development, so this is a development reference.

## 4. Carrier evidence

**Tiers:**
- **verified:** fetched this session by the read-only reviewer subagent, and the quoted passage was re-checked by me
  in the saved response (file in brackets);
- **second-hand:** a source says it, but its own evidence was not opened;
- **established / uncited:** textbook knowledge, no source opened;
- **inferred:** reasoned.

*S. cerevisiae*, human and plant sources are related-organism context, not *Y. lipolytica* evidence.

| Claim | Source | Tier |
| --- | --- | --- |
| Leu5 deletion lowers mitochondrial CoA 15-fold while cytosolic CoA is unchanged. Leu5 acts "presumably by serving as a transporter of CoA or a precursor thereof". No substrate assay | Prohl et al. 2001, [PMID 11158296](https://europepmc.org/article/MED/11158296) (f01) | verified (abstract) |
| Leu5 is in the inner membrane | [UniProt P38702](https://www.uniprot.org/uniprotkb/P38702/entry) (f09, ECO:0000269) | verified (record) |
| SLC25A42 (reconstituted) has "high transport affinity for CoA, dephospho-CoA, ADP, and adenosine 3',5'-diphosphate". No acyl-CoA is reported tested | Fiermonte et al. 2009, [PMID 19429682](https://europepmc.org/article/MED/19429682) (f01) | verified (abstract) |
| The fly SLC25A42 orthologue transports dephospho-CoA "but not CoA", and rescues yeast *leu5Δ* | Vozza et al. 2017, [PMID 27836698](https://europepmc.org/article/MED/27836698) (f27) | verified (abstract) |
| Human SLC25A16/A42 double knockout gives "severe depletion of several mitochondria-enriched short-chain acyl-CoAs" (propionyl-, succinyl-/methylmalonyl-CoA). So the matrix makes these esters from imported free CoA | Liu et al. 2025, [PMID 40925986](https://europepmc.org/article/MED/40925986) (full text, f25) | verified (full text); the conclusion is inferred |
| The paper also says "CoA species are charged and therefore cannot freely cross the mitochondrial inner membrane". This is a background sentence, not a measurement | Liu et al. 2025 (f25) | verified (sentence); its claim is established / uncited |
| Crc1 "transports carnitine, acetylcarnitine, propionylcarnitine". Acyl groups cross as acylcarnitines, not as CoA esters | Palmieri et al. 1999, [PMID 10622748](https://europepmc.org/article/MED/10622748) (f02) | verified (abstract) |
| In yeast, membranes between compartments "are impermeable to acetyl-CoA". Background statement, about acetyl-CoA only | Chen et al. 2012, [PMID 22876324](https://europepmc.org/article/MED/22876324) (f02) | verified (sentence); claim is second-hand |
| No 4'-phosphopantetheine transporter is known: "whether specific transporters for PPanSH exist" is an open question | Xue et al. 2023, [PMID 36316571](https://europepmc.org/article/MED/36316571) (full text, f26) | verified (full text) |
| Yeast Cab4 (PPAT) is in the cytoplasm and nucleus | [UniProt P53332](https://www.uniprot.org/uniprotkb/P53332/entry) (f09, ECO:0000269, high-throughput) | verified (record) |
| Yeast Cab5 (DPCK) is an "extrinsic component of the mitochondrial outer membrane that also localizes to lipid droplets" | [SGD CAB5](https://www.yeastgenome.org/locus/S000002604) (f24) | verified (record) |
| Human COASY is mostly cytosolic, with a matrix pool. Matrix-targeted COASY does not rescue yeast *leu5Δ* on glycerol ("neither COASY nor MTS-COASY was able to rescue") | [UniProt Q13057](https://www.uniprot.org/uniprotkb/Q13057/entry) (f08); Liu et al. 2025 (f25) | verified. That this argues against a working matrix synthesis route in yeast is inferred |
| Targeted searches for a carrier of propionyl-, succinyl-, isovaleryl-, 3-methylcrotonyl- or isobutyryl-CoA, or 4'-phosphopantetheine, found none (48 and 29 hits) | f03, f21 | inferred (negative search) |

**What this supports.**
- A CoA ester crossing the inner membrane as a whole molecule has no known carrier in any organism searched.
- What matrix CoA biology needs is a free CoA (or dephospho-CoA) import (Leu5, the gene on R834), with acyl-CoAs made
  inside the matrix. In a steady-state model, that import is a conserved-moiety exchange that needs a growth-dilution
  demand; see the R1172 report, O4.

## 5. Do R2114–R2117 belong in the matrix?

| Step (model) | *Y. lipolytica* gene (W29 / CLIB122 / UniProt) | Name | TargetP-2.0 (prediction) | Homologue location |
| --- | --- | --- | --- | --- |
| R2117 isobutyryl-CoA DH | `YALI1D19252g` / YALI0D15708g / [Q6C8Y9](https://www.uniprot.org/uniprotkb/Q6C8Y9/entry) | none verified (TrEMBL: "short/branched chain specific acyl-CoA dehydrogenase, mitochondrial", PE3, automatic) | mTP 0.9992 | human ACAD8: mitochondrion, ECO:0000269 ([Q9UKU7](https://www.uniprot.org/uniprotkb/Q9UKU7/entry)) |
| R2116 enoyl-CoA hydratase | `YALI1B14044g` / YALI0B10406g / [Q6CF43](https://www.uniprot.org/uniprotkb/Q6CF43/entry) | none verified (TrEMBL: "probable enoyl-CoA hydratase, mitochondrial") | 346 aa (W29 RefSeq): **Other 0.69**, mTP 0.30. 319 aa (from Met28, the CLIB122 length): **mTP 0.9998** | human ECHS1: matrix ([P30084](https://www.uniprot.org/uniprotkb/P30084/entry)) |
| R2115 HIBCH | no gene in the model. Candidate `YALI1D07981g` / YALI0D06215g / [Q6CA37](https://www.uniprot.org/uniprotkb/Q6CA37/entry) | none verified. TrEMBL: "3-hydroxyisobutyryl-CoA hydrolase"; RefSeq: "mitochondrial 37S ribosomal protein mS47" | mTP 0.9993 | human HIBCH: mitochondrion, ECO:0000250 ([Q6NVY1](https://www.uniprot.org/uniprotkb/Q6NVY1/entry)). Yeast EHD3 is mitoribosomal mS47 with weak HIBCH activity ([P28817](https://www.uniprot.org/uniprotkb/P28817/entry)). Arabidopsis CHY1 is peroxisomal ([PMID 11404361](https://europepmc.org/article/MED/11404361)) |
| R2114 HIBADH | `YALI1F03883g` / YALI0F02607g / [Q6C351](https://www.uniprot.org/uniprotkb/Q6C351/entry) | none verified | mTP 0.9829 | human HIBADH: mitochondrion ([P31937](https://www.uniprot.org/uniprotkb/P31937/entry)) |
| R2112 MMSDH (acylating) | `YALI1C02550g` / YALI0C01859g / [Q6CDB9](https://www.uniprot.org/uniprotkb/Q6CDB9/entry) | none verified | mTP 0.9995 | human ALDH6A1: mitochondrion, ECO:0000305 ([Q02252](https://www.uniprot.org/uniprotkb/Q02252/entry)) |
| (R2082, leucine) isovaleryl-CoA DH | `YALI1E15555g` / YALI0E12573g / [Q6C643](https://www.uniprot.org/uniprotkb/Q6C643/entry) | none verified | mTP 0.9998 | — |
| control: PDH E1β | `YALI1E31947g` | none verified (list text: similar to PDB1) | mTP 0.7431 | — |
| (R265/R1092) DPCK | `YALI1F13202g` / YALI0F09625g | none verified | Other 0.9909 | yeast Cab5: outer membrane and lipid droplets |
| (R615) PPAT | `YALI1D24649g` | none verified | Other 0.633 (reused, orphan_gene_search_20261019) | yeast Cab4: cytoplasm and nucleus |

Evidence levels:
- Every *Y. lipolytica* location annotation is automatic (ARBA, IBA or IEA), none experimental (UniProt records f13).
- The TargetP values are predictions. They do not confirm compartment, catalysis or identity.
- The model's gene rules are not experimental validation.

**Other evidence:**
- Arabidopsis has a mitochondrial valine route (CHY4, MMSD, HDH1), whose disruption affects seed development
  ([PMID 28321931](https://europepmc.org/article/MED/28321931); f06; verified, abstract).
- *Aspergillus nidulans* needs a mitochondrial β-oxidation acyl-CoA DH to use valine or isoleucine
  ([PMID 17656140](https://europepmc.org/article/MED/17656140); f06; verified, abstract; the abstract shows no
  localisation experiment).
- *Y. lipolytica* data are transcript, omics or tracer level only, with no compartment:
  [PMID 37683719](https://europepmc.org/article/MED/37683719), [PMID 41301536](https://europepmc.org/article/MED/41301536),
  [PMID 38475794](https://europepmc.org/article/MED/38475794) (f12, f20; verified, abstracts).

**Verdict (inferred).**
- The valine steps after the BCKDH belong in the matrix. Homologues are mitochondrial everywhere they are known, the
  *Y. lipolytica* proteins carry strong predicted presequences, and the model's own matrix BCKDH product has no
  other consumer.
- R2117's mixed compartments, FAD (cy) → FADH2 (mi), are a curation artefact of the cytosolic placement.
- **Not settled:**
  - the hydratase's start codon (two gene models give opposite TargetP calls);
  - which gene, if any, carries R2115. `YALI1D07981g` is on the consensus list, but its yeast homologue is
    primarily a mitoribosome protein, so its essentiality need not come from valine catabolism (inferred).

## 6. Options tested in memory (step 4)

[tables/options_summary.csv](tables/options_summary.csv), [tables/option_edits.csv](tables/option_edits.csv).

**Bases:**
- **B0:** 7ccbb568 as is.
- **B1:** B0 plus the methylcitrate line's package: T_OBDHm, a matrix 2-oxobutanoate dehydrogenase with the
  inferred PDH rule of R693, and R2075 closed. This matches candidate b24a5e78 on cbf2ef95 (uncommitted, other line)
  and the R1172 round's "O1 + R2075 closed".

**Variants:**
- **K:** keep.
- **M:** merge m1875 into m1970.
  - R2024 then reads `ATP + propanoate + CoA ⇌ propionyl-CoA + PPi + AMP + 3 H+`. The 3 H+ close H and charge
    together (`check_mass_balance` = {}).
  - R2025 becomes an exact duplicate of R2113 and is set to 0.
- **MC:** M plus R2113 closed.
- **C6:** the six blocked reactions closed.
- **V:** R2117, R2116, R2115, R2114 and R2112 moved into C_mi, with the same rules and bounds, and R2118 closed.
  - Four new matrix species (formulas inferred from the model's matrix acyl-CoA form and KEGG):
    methacrylyl-CoA C25H40N7O17P3S, 3-HIB-CoA C25H42N7O18P3S, 3-HIB C4H7O3⁻ and MMSA C4H5O3⁻.
  - All five moved reactions balance, after 0, 0, 1, 0 and 3 H+.
- **ALL:** M + MC + V, with R2077, R2083, R2084 and R1291 closed.

| Variant | WT (Gurobi; GLPK equal to 1e-14) | Energy check, ATP / biomass (both solvers) | Can carry flux (Task A, no floor), selected | Targeted knockout class changes (30 genes) |
| --- | --- | --- | --- | --- |
| B0/K | 0.311100535925948 | 0 / 0 | R2025, R2113 (2.67, pass-through) | — (step 1: none) |
| B0/M | 0.311100535925949 | 0 / 0 | R2024, R2285, R2022 (2.67). **R2113 blocked** | 0 |
| B0/MC | 0.311100535925949 | 0 / 0 | as B0/M (LP-identical: R2113 was already blocked) | not run |
| B0/V | 0.311100535925949 | 0 / 0 | R1969, moved chain, R25/R95/R552 (5.86). **R834 2.70** through R2113 in and R2025 out | 0 |
| **B0/ALL** | 0.311100535925949 | 0 / 0 | R1969, moved chain, R25/R95/R552 (5.86); R2024 outlet (2.67). **R834 blocked** | 0 |
| B1/K | 0.311100535925949 | 0 / 0 | trio 7.60 (T_OBDHm); R834 3.31 (side route B, out on R2025) | 0 |
| B1/M | 0.311100535925949 | 0 / 0 | trio 7.60; R834 3.31, now out on R2113 | not run |
| B1/MC | 0.311100535925947 | 0 / 0 | trio 7.60; **R834 blocked** | 0 |
| B1/C6 | 0.311100535925948 | 0 / 0 | as B1/K (the six stay blocked under O1) | not run |
| B1/V | 0.311100535925949 | 0 / 0 | trio 8.41; R834 3.38 | not run |
| **B1/ALL** | 0.311100535925947 | 0 / 0 | trio 8.41 (valine 5.86 + 2-oxobutanoate); **R834 blocked** | 0 |

**What the table shows.**
- **The merge removes the reason for a propionyl-CoA carrier.** Isoleucine-derived propionyl-CoA reaches R2024 in
  the cytosol directly, and R2113 cannot carry flux at all on B0. Closing it is then a no-op on B0. On B1 it is the
  only exit of side route B, the CoA-in / propionyl-CoA-out loop the R1172 round flagged.
- **Moving the valine steps gives the methylcitrate cycle a matrix propionyl-CoA source with gene rules.** The source
  is R1969 → R2117 → … → R2112, and the cycle is active without the inferred T_OBDHm. Every moved reaction conserves
  the matrix CoA moiety (static, from the equations):
  - R2115 releases CoA;
  - R2112 binds it again;
  - R25 releases it.

  So no exit is needed.
- **R834 is active only while a propionyl-CoA transport is open.** Once none is, R834 is blocked in every variant.
  This is the matrix CoA identity of the R1172 report: a carrier that only imports a conserved moiety cannot carry
  net flux in a steady-state model without a matrix CoA demand.
- **No variant opens energy, changes WT, or changes a targeted call.** The changes are curation of compartments and
  species, not of predictions on this medium.

## 7. Figures ([figure/](figure/))

Hand-laid schematics: one box per reaction with its gene rule (and the TargetP call) and its flux. Green = gene
rule, pink = no gene, black outline = moved copy (in memory), dashed dark red with a cross = closed in the option.

- [acylcoa_status.png](figure/acylcoa_status.png): the largest flux each reaction can carry, as is / option ALL (B0).
- [acylcoa_asis_route.png](figure/acylcoa_asis_route.png): 7ccbb568 as is, one propionyl-CoA forced into the matrix
  (R2113 = −1, `xMAINTENANCE` lower bound 0, minimum total flux). It comes from the cytosolic isoleucine branch,
  passes through the matrix (R2113 in, R2025 out), and leaves as propanoate.
- [acylcoa_option_route.png](figure/acylcoa_option_route.png): option ALL (B0), with R25 = 1 and R2024 = −1 forced.
  Valine is degraded to propionyl-CoA and through the methylcitrate cycle inside the matrix (FADH2 to R1977), and
  isoleucine-derived propionyl-CoA leaves as propanoate from the cytosol.

## 8. What it means, and could the method be better

**The eight do not carry biology.**
- R2025 and R2113 are the two halves of a detour that exists only because the model has two cytosolic propionyl-CoA
  species.
- The other six are blocked leftovers of placing the branched-chain steps and a CoA-synthesis copy outside their
  compartment.
- None is needed for growth or for any targeted gene call.

**Closing them alone does not fix the gaps they belong to.** The fix toward the biology is:
- one cytosolic propionyl-CoA species;
- the valine chain inside the matrix, using its existing gene rules;
- the branched-chain and CoA-synthesis copies that only shuttle esters removed.

That activates R1969, R2114, R2117 (as matrix copies) and the methylcitrate trio with gene-supported steps, which is
what AGENTS.md section 5 asks for.

**What stays open after that:**
- R834. A steady-state model cannot use a moiety-only carrier without a matrix CoA dilution term (R1172 report, O4).
- The cytosolic isoleucine branch R2281–R2285. It has the BCKDH rule in the cytosol; R2282 has no electron acceptor;
  and the matrix branch R1970 ends at the dead-end m1830. Same kind of gap, not in this task.
- The methylmalonyl route R2074–R2076. It has no mutase gene.

**A better method.** The compartment mismatches found here are all visible from one check: compare each
gene-assigned reaction's compartment with the predicted or annotated location of its gene. R2082, R2117, R2114 and
R2112 all carry genes with mTP above 0.98 but sit in the cytosol.
- **How to test it:** run TargetP on every model gene once, cached with the server stamp. List the reactions whose
  compartment disagrees with all of their genes. Check whether the list recovers these five, then review the rest.
- **Two checks to pair with it:**
  - the per-compartment conserved-moiety scan proposed by the R1172 round;
  - the "same species twice in one compartment" check: name, KEGG ID or formula duplicates, as for m1875/m1970.

**What is proven and what is observed.**

| Statement | Kind |
| --- | --- |
| The six reactions are blocked; R2025 and R2113 carry flux one way only | observed: LP, two solvers, two conditions |
| Closing a blocked reaction cannot change any knockout | provable |
| None of the eight is needed for WT or for the 30 targeted knockouts | observed (LP) |
| After the merge, R2113 is blocked on B0 | observed (LP) |
| The moved valine chain conserves the matrix CoA moiety | statically verified (equations) |
| No carrier is known for these CoA esters or for 4'-phosphopantetheine | verified sources plus a negative search; absence of evidence, not of the carrier |
| The *Y. lipolytica* valine enzymes are mitochondrial | inferred (homologues, plus TargetP predictions); untested experimentally |
| No gene changes class under option B0/ALL, at any cutoff, under Rule A or B | observed: full certified sweep, Task A + PO1f, Gurobi 1e-9 (section 11) |
| The as-is baseline (353cefdc sweep) is the same optimisation problem as 7ccbb568 | statically verified (step6_identity.py, no LP) and reproduced for 12 genes plus WT (difference 0.0) |
| No gene changes class under the other ten variants | untested (no sweep run for them) |

## 9. Limits

- **One medium** (Task A + PO1f) for the options; the open condition was used for capacity only.
- **Knockouts:** 30 targeted genes for six variants. A full certified sweep was run only for B0/ALL (section 11), not
  for the B1 variants or the partial options.
- **Formulas of the four new matrix species and the 3 H+ in the merged R2024** are inferred: the model's protonation
  forms are mixed (CoA neutral, propionyl-CoA −4).
- **TargetP:**
  - one predictor, as authorized; MitoFates and DeepLoc not run;
  - it does not predict peroxisomal targeting. YALI1F13202g (DPCK) ends in "SKL", a sequence observation only, with
    no PTS1 predictor run;
  - the YALI1B14044g call depends on the gene model.
- **YALI0D12573g (on R1971):** UniProt and UniParc return no such ORF (f13, f17). KEGG no longer lists YALI0 tags, so
  its 404 (f14) is not informative. That it is a mistyped YALI0E12573g is inferred.
- **Literature:**
  - mostly abstracts and database records; three full texts (Liu 2025, Xue 2023, Zheng 2025);
  - Dusi 2014 and Desai 2017 were not opened (second-hand via UniProt);
  - no *Y. lipolytica* compartment data exist in what was found.
- **Gene names:** none verified. UniProt names are unreviewed TrEMBL annotations.
- **Not addressed:**
  - R2119/R2121, which share R2112's gene and stay cytosolic;
  - the peroxisomal and matrix reactions that share YALI1B14044g;
  - the 9 YALI0 tags in the model gene list.

## 10. Files and commands

- **Interpreter.** `/opt/miniconda3/bin/python3` with the real HOME (Gurobi licence), run from [scripts/](scripts/).
  Every script checks the model, medium, profile and list SHA-256 first.
- **Commands, in order:**
  1. `step1_capacity.py`
  2. `step4_options.py` (the `__main__` guard was added after the run; logic unchanged; the as-run copies and their
     SHA-256 are in `scripts_as_run/` in the artifacts folder)
  3. `step5_figure_flux.py`
  4. `draw_schematic.py ../figure`
  5. `make_tables.py ../tables`
  6. Step 6 (section 11): `step6_identity.py`, then `step6_sweep.py B0/ALL …` from `platform/` (exact command in
     section 11), the as-is control `step6_sweep.py B0/K … ../tables/step6_control_genes.txt`, then
     `step6_compare.py $ART/out/sweep_B0K_control.tsv` and `step6_r171_route.py`.
- **Outputs:** `out/step1_capacity.json`, `out/step4_options.json`, `out/step5_figure_flux.json`; logs in `logs/`;
  literature in `lit/`; TargetP files in `targetp/`. All in
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/geneless_acylcoa_transports_20261116/` (`MANIFEST.sha256`).

## 11. Full certified knockout sweep of option ALL (step 6, follow-up 2026-11-16)

**Authorization.** The user: **"Run the sweep"**, answering "Run a full certified knockout sweep (about 1,100 LPs) on
the leading option?". Scope, budget and stop conditions: [TASK.md](TASK.md) section 6. No model, curation, benchmark
or test change; nothing is written to SBML.

**How it was run.**
- **Harness.** [scripts/step6_sweep.py](scripts/step6_sweep.py) is a copy of the benchmark's reference
  implementation `infeasible_ko_rule_20261030/scripts/ko_screen_certified.py` (`313d5d01…`), with `checked_solve.py`
  (`69c42690…`) unchanged. The diff to the reference touches only:
  - the loader: every context (sweep, fresh re-solve, certificate) reads 7ccbb568 through
    `load_effective_simulation_context`, then applies the variant with step 4's `apply_base` and `apply_variant`;
  - an input SHA-256 gate;
  - an optional gene subset (used only for the control);
  - three extra summary fields.

  The rest is the reference code: Gurobi 13.0.2 at 1e-9, the model reset per knockout, the targeted check,
  `FRESH_LIMIT` 20, the four-condition maintenance certificate, Rule A and Rule B, and the mtDNA and untargeted
  scopes.
- **Baseline (as is).** The certified sweep `sweep684_candidate` (folinic-acid round, `49b4f27b…`), run with the same
  two harness files, medium, profile, list and untargeted file on 353cefdc, the parent build of 7ccbb568. It is used
  because:
  - **statically verified, no LP** ([tables/step6_identity.json](tables/step6_identity.json)): both models load
    through the Task A + PO1f context as the same problem. Reaction and metabolite order, stoichiometry, effective
    bounds, gene rules, objective, genes and runtime exclusions all match (2,347 reactions, 1,082 genes; 0 reactions
    differ);
  - **reproduced:** the as-is control (`B0/K` in the step-6 harness) on the 12 genes whose ratio moved gives WT and
    every ratio equal to the baseline, difference 0.0.
- **LPs:** 1,141 of the 1,300 budgeted:
  - the sweep: 1,082 solves, plus the WT ceiling and 40 certificate solves (1,123);
  - the control: 14;
  - the R171 route: 4.

  Every sweep point passed the targeted check: knocked flux 0, bound violation 0, residual ≤ 9.2e-12. No fresh
  re-solve was needed.

**Command** (from `platform/`, real HOME; `$T` = this folder, `$WT` = worktree root, `$RR` =
`$IYALI26_RESEARCH_ROOT`, `$ART` = this round's artifacts folder):

```
PYTHONPATH="$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts:$T/scripts" /opt/miniconda3/bin/python3 \
  "$T/scripts/step6_sweep.py" B0/ALL "$WT/model/candidates/E5_353cefdc_xref_20261115.xml" \
  "$WT/model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv" "$RR/state/strain_profiles/po1f_sd_leu.json" \
  "$RR/state/essentiality/repository/consensus_essential_genes.csv" "$WT/model/expected/untargeted_genes_out_of_scope.json" \
  "$ART/out/sweep_B0ALL.tsv" "$T/tables/sweep_B0ALL_summary.json"
```

**Results** ([tables/sweep_scores.csv](tables/sweep_scores.csv), [tables/step6_compare.json](tables/step6_compare.json),
[tables/sweep_B0ALL_summary.json](tables/sweep_B0ALL_summary.json); observed, LP):
- **WT:** 0.31110053592594916, against 0.311100535925948 as is (+1.2e-15).
- **No class change:** 0 of 1,081 genes change class at any of the four cutoffs, under Rule A or Rule B
  ([tables/sweep_class_changes.csv](tables/sweep_class_changes.csv) has only its header).
  - The scores are identical to the baseline in all four rule and scope combinations.
  - The 20 certified maintenance-shortfall knockouts are the same 20 genes, and none is left uncertified.
  - The 30 targeted knockouts of step 4 match the sweep rows exactly (largest difference 0.0).
- **Scores, Rule B, mtDNA and untargeted genes out of scope:**

  | Cutoff | TP | FN | FP | TN | Recall | Specificity | Precision |
  | --- | --- | --- | --- | --- | --- | --- | --- |
  | 1 % | 117 | 209 | 36 | 701 | 0.359 | 0.951 | 0.765 |
  | 5 % | 124 | 202 | 38 | 699 | 0.380 | 0.948 | 0.765 |
  | 10 % | 128 | 198 | 43 | 694 | 0.393 | 0.942 | 0.749 |
  | **15 %** | **145** | **181** | **43** | **694** | **0.445** | **0.942** | **0.771** |

  - **Rule A** (every non-optimal knockout unresolved): the same, except the certified knockouts are unresolved: 15 in this
    scope, 11 positives and 4 negatives. At 15 %: TP 134, FN 181, FP 39, TN 694, unresolved 15.
  - **Negatives** are the user-defined class: model genes not on the consensus list.
  - **Out of scope:**
    - the 16 mtDNA genes, none of them on the list;
    - the 2 untargeted genes, none on the list;
    - 1,286 of the 1,612 listed genes have no model gene; 326 listed genes are in the model.
  - The list was used in model development, so these scores are a development reference, not independent
    validation.
- **Ratio moves without a class change:** exactly 12 genes move by more than 1e-6. They are the 12 subunit genes of
  ATP synthase R171 (AND rule), 9 on the list and 3 mtDNA. Each knockout removes only R171, and KO/WT goes from
  0.107782 to 0.108531. This is non-essential at 10 % and essential at 15 % both before and after.
- **Why** ([scripts/step6_r171_route.py](scripts/step6_r171_route.py),
  [tables/step6_r171_route.csv](tables/step6_r171_route.csv)). One pFBA of the R171 knockout (YALI1F04433g), as is
  and with the option, gives growth 0.0335309 → 0.0337639.
  - In the option's optimum, valine uptake rises from 0.0123 to 0.0260.
  - The extra 0.0136 runs through the matrix valine transaminase R797, then R1969 (BCKDH rule), then the moved
    R2117_mi → R2112_mi, then the methylcitrate cycle R25 → R95 → R552.
  - Its FADH2 goes to R1977 (ETF:ubiquinone oxidoreductase, 0 → 0.0136), so the respiratory chain gains a little
    flux.
  - Without ATP synthase, that extra carbon and reducing power is worth 0.7 % growth.
  - The WT maintenance ceiling (the largest `xMAINTENANCE` with its lower bound at 0) rises the same way, from
    63.3334 to 63.3542.
  - The energy check with uptakes closed stays 0 (step 4), so this is medium-supplied valine, not a free-energy loop.
  - One minimum-flux optimum shows one route; it does not prove the route is unique.
- **Alternate optimum (not caused by the option):** both pFBA points also carry an NADH-oxidising pair, a long-chain
  alcohol dehydrogenase run backwards plus the alcohol oxidase. It runs on C10 as is (R1864/R1882) and on C12 with the
  option (R1866/R1883). It is a degenerate choice between equal routes, and it affects another line, not this task.

**What it means.** On Task A + PO1f the leading option is a curation change of compartments and species with no
effect on any gene call. The only phenotype it moves is the ATP-synthase knockouts, by 0.7 %, through a gene-supported
matrix valine route. Before a build, the sweep does not need repeating for the same edits, but it does for any
different edit set, for example adding R2077 or the B1 package.

**Limits of the sweep.**
- One medium.
- Only the B0/ALL edits; the inferred formulas and H+ of step 4 apply unchanged.
- The baseline is a different file (353cefdc) that is statically the same problem, with 12 genes plus WT re-run.
- R2077 is closed here as part of ALL, but deciding on it stays with the lipid-unlump owner.

**Outputs.** `$ART/out/sweep_B0ALL.tsv` (per gene), `$ART/out/sweep_B0K_control.tsv` and `.summary.json`,
`$ART/out/step6_r171_route.json`, and logs in `$ART/logs/` (`step6_inputs_sha256.txt` holds the SHA-256 of every
input and script).
