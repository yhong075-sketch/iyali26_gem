# Task D: complex I transhydrogenation module (holistic combinatorial test)

**Status (2026-10-14).** Investigation and diagnostics only. No model was built or edited, no condition file
was written, nothing was committed. Scope, rules and the frozen part table are in [TASK.md](TASK.md), written
before the scored runs. Everything is uncommitted on worktree branch `claude/task-d-78f16d`.

**Outcome.**
- **The escape reproduces** (§2). In the complex I KO, PDH and CS carry 0 and glucose is oxidised in the
  cytosol (G6PDH 221, NDH2 578 per 100 glucose). The MTHFD and GDH pairs are interchangeable NADPH → NADH
  converters. KO/WT is 0.826 (Task A medium) and 0.827 (W29 V2).
- **The cytosolic transhydrogenation module does not decide complex I essentiality** (§3–§5).
  - Discovery needed 24 parts to push the KO below 0.15.
  - The 2^12 test then showed that **two single closures are each enough**, from the open state, at 0 % WT cost:
    - **G5:** matrix MDH R533 in the reductive direction (OAA + matrix NADH → malate);
    - **G3:** the citrate/malate carrier R831 in the citrate-in, malate-out direction. Tier 2 resolves G3 to
      R831 alone.
  - Both give KO/WT 0.077 (Task A medium) and 0.000 (V2). All 3,072 hits are exactly the sets that contain
    G3 or G5.
- **Mechanism (verified):**
  - Without complex I, the only matrix NADH that must still be made comes from the lysine-pathway
    homoisocitrate dehydrogenase R398 (2.4 per 100 glucose).
  - Its only outlet is R533 reverse, with malate leaving through R831.
  - Opening lysine uptake rescues the blocked KO (growth 0 → 0.227 h⁻¹ in V2).
  - The result does not depend on the overlay or the 13C caps.
- **Step 4: an assumption is needed.**
  - Both cut sets close a direction of a real enzyme or carrier. No source supports closing the MDH direction,
    and thermodynamics favours it (eQuilibrator estimate ΔrG′° +29.6 kJ/mol for malate oxidation).
  - So complex I essentiality in the model **needs a capacity or regulatory assumption** on matrix-to-cytosol
    NADH export.
  - *Y. lipolytica* evidence (Kerscher 2001; Garofano 2006) says such shuttles are insufficient and that
    complex I is required on glucose, but gives no number.
  - **New, cross-species:** reconstituted *S. cerevisiae* Ctp1 does not exchange citrate for malate (Aluvila
    2010), and Bakker 2001 rules out a malate–citrate shuttle in *S. cerevisiae*. The R831 malate mode that G3
    closes therefore has weak support. One older counterpoint, second-hand: malate is a weak substrate in
    isolated mitochondria (Evans 1983 as cited by Kaplan 1995). This also affects Task C's fit.
  - **The options are in §7. The user chose option (d): no ASSUMPTION file.**
- **Step 5 (13C ratio floors):** no accepted floor set makes complex I essential (lowest KO/WT 0.316, V2, with
  PDH + CS + G6PDH and the R825 ASSUMPTION).
  - The floors raise the KO's matrix-NADH export to 89–116 per 100 glucose at its optimum (FVA).
  - With a shuttle cap the KO has no solution below 86.66 (the PDH-floor NADH). It is essential only in a band
    about 0.24 wide (A) or 0.46 wide (V2) above that (§7, §8).

## 1. Inputs and identities (verified; every LP script checks them before running)

| Item | Value |
|---|---|
| Worktree / HEAD | `claude/task-d-78f16d` @ `4d97b828ac2b0dc95d592060605af75069705c09` (main with PRs #4 and #5), clean at start |
| Model | `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml` (Task B candidate, R132 removed), `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515` |
| Task A medium | `po1f_csm_leu_w29uptake_20261012.csv`, `6eca70fd…` |
| 13C caps (ASSUMPTION) | `w29_13c_upper_caps_20261012.json`, `df1e33a8…`, "module" glycine-cleavage form |
| R825 reversible (ASSUMPTION) | `r825_reversible_20261013.json`, `e45cab42…`; step 5 and §7 only, every such row labelled |
| Overlay (diagnostic) | [holistic_overlay.py](holistic_overlay.py) `bf6c7ce5…`, byte-identical to Task C's |
| Profile / MFA | `po1f_sd_leu.json` `35307853…`; W29 `glucose_central_fluxes.csv` `844615d0…` (repository ↔ paper link inferred, pending user confirmation) |
| Part definitions | [parts_definition.json](parts_definition.json), final SHA `cec5c146…` (recorded in `combinatorial.json`) |
| Solver | Gurobi 13.0.1, 1 thread, Presolve 0, tolerances 1e-7; cobra 0.30.0; Python 3.13.5 (main checkout `.venv`). Maps are drawn with the conda Python (matplotlib), with no LP |

**Configurations:**
- **A (primary):** Task A medium + PO1f profile, glucose uptake ≤ 2.414.
- **V2 (sensitivity):** W29-type defined medium, glucose fixed at 2.414.

Both use the holistic overlay and the 13C caps. The caps are WT-state bounds applied to WT and KO alike; this
is a diagnostic, not a screen. Fluxes are per 100 glucose.

**LP count:** 22,423 registered LPs (each JSON records its own count), plus about 1,050 pre-plan exploration
LPs in the session scratchpad. The cap was 40,000. Large outputs are in
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/complexI_transhydrogenation_20261012/` (pFBA vectors, sequences,
evidence texts, logs).

## 2. Step 1: the escape (verified: [escape.json](escape.json), [escape_node_balances.csv](escape_node_balances.csv))

| | A WT | A KO | V2 WT | V2 KO |
|---|---:|---:|---:|---:|
| growth (h⁻¹) | 0.28222 | 0.23322 | 0.26271 | 0.21718 |
| KO/WT | | **0.826** | | **0.827** |
| G6PDH R325 (FVA) | [0, 12.0] | 221.0 | [0, 11.5] | 208.3 |
| NDH2 R570 | 193–197 | 577.7 | 203–206 | 563.7 |
| complex I R2062 | 194–202 | 0 | 182–189 | 0 |
| PDH R693 / CS R236 (FVA) | 62–66 / 49–53 | 0 / 0 | 59–63 / 46–50 | 0 / 0 |
| NAD-MTHFD R549 (FVA) | 0 | [0, 446.5] | 0 | [0, 419.6] |
| NAD-GDH R329 (FVA) | 0 | [0, 446.5] | 0 | [0, 419.6] |
| matrix NADH produced (pFBA) | 199.3 | 2.4 | 186.9 | 2.5 |
| cytosolic NADPH produced (pFBA) | 83.2 | 458.6 | 87.7 | 431.4 |

**Reading (verified by the FVA ranges).**
- **Respiration:** the KO stops the matrix TCA cycle and respires through external NDH2.
- **NADPH comes from the oxPPP:** G6PDH is pinned at 221 (A) or 208 (V2).
- **NADPH → NADH conversion:** either the MTHFD pair or the GDH pair converts it to NADH, interchangeably (each
  can carry up to 446). pFBA picks GDH on A and MTHFD on V2 (R549 363, R550 −360). These are alternate optima,
  not two mechanisms.
- **Matrix NADH falls from 199 to 2.4.**

**Against the prompt's recalled values** (recalled; checked here):
- PDH and CS 0: reproduced.
- G6PDH about 210: 208 (V2) and 221 (A).
- MTHFD pair about 368: 363 (V2, pFBA).
- NDH2 about 564: 563.7 (V2).
- The recalled KO/WT range 0.85–0.93 is not reproduced (0.826). This fits Task B's R132 removal, which lowered
  the ratio (0.854 → 0.758 without overlay or caps). Task D has no 2 % stop on this number.

## 3. Part list and discovery (verified: [parts_pool.json](parts_pool.json), `discovery_round1..4.json`, `step_costs_round1..3.json`)

**Pool (configuration A).**
- **Pair scan:** it reproduces the prompt's 14 NAD/NADP pairs exactly (H+ ignored). Six drop out because a
  member is blocked or closed by the overlay: the two ADH pairs, R2069/R131, R254/R255, R534/R_MDHy and
  R737/R738.
- **NAD(P)(H) transporters:** the four mitochondrial NAD⁺ carriers R903–R906 are blocked (FVA [0, 0]).
- **The four NAD(H) kinases** are blocked.
- **Escape directions** (usable at the KO optimum, not required by WT):
  - R549 forward, R550 reverse, R329, R1378, R136, R1385 reverse and R1708. Each closes at 0 % WT cost, and each
    alone leaves KO/WT at 0.82–0.83.
  - Trace escape directions: R1628 forward and R1629 reverse (0.026), and R1936 forward (0.020).
  - R130, R1377 and R2069 also show a formal escape direction, but their FVA reaches the ±1000 bound (net-zero
    loops). They dropped with their pairs.

**Capacity steps** (the KO uses them in a direction WT also requires, so they are not closed in the test).
TASK.md asked for their whole-reaction closure cost. It was missing in the first draft, flagged by the audit,
and run afterwards (verified: [capacity_costs.json](capacity_costs.json)). Each is closed whole, alone, from
the open state:

| Step | WT cost A / V2 | KO/WT A / V2 |
|---|---|---|
| R134 NADP-ALD | 1.9 % / 2.2 % | 0.842 / 0.845 |
| R330 NADP-GDH | 1.3 % / 1.6 % | 0.810 / 0.809 |
| R691 PDC | 8.0 % / 9.3 % | 0.732 / 0.714 |
| R361 glycine cleavage (lipoamide) | 7.7 % / 7.4 % | 0.895 / 0.892 |
| R534 cytosolic MDH | 49.7 % / 62.4 % | 0.544 / 0.424 |
| R641 / R660 serine pathway | 3.7 % / 4.0 % | 0.825 / 0.825 |
| R659 phosphoserine phosphatase | 0 / 0 | 0.826 / 0.827 |
| R740 SDH | 17.8 % / 19.7 % | 0.946 / 0.943 |
| R929 succinate/fumarate carrier | 4.5 % / 6.1 % | 0.814 / 0.806 |
| R490 ICL | 4.5 % / 5.5 % | 0.815 / 0.812 |
| R315 cytosolic fumarase | 9.9 % / 10.9 % | 0.819 / 0.810 |
| R1385 matrix NADP-IDH | 0.1 % / 0.2 % | 0.825 / 0.827 |
| R570 NDH2 | 3.3 % / 3.7 % | 0.634 / 0.630 |
| R325 G6PDH | 0 / 0 | 0.793 / 0.794 |

No single capacity step makes complex I essential: every closure leaves KO/WT ≥ 0.42, and 10 of the 15
steps (R641 and R660 counted separately) cost WT more than 2 % in A.

**Discovery rounds (configuration A).**

| Round | Closed | WT | KO/WT | New route that took over (WT-free closure) |
|---|---|---:|---:|---|
| 0 | none | 0.28222 | 0.826 | — |
| 1 | 9 parts | 0.28080 | 0.540 | serine-pathway cycle (R511); GABA shunt (R328); glutamate-semialdehyde/P5C (R1931); matrix leucine degradation (R1968); matrix FAS + peroxisomal β-oxidation (R2230); citrate import (R817 rev) |
| 2 | 15 | 0.28080 | 0.436 | serine → pyruvate via cystathionine (R244); glyoxylate secretion (R1148) |
| 3 | 17 | 0.28080 | 0.406 | malonate-semialdehyde cycle (R2119); polyamine catabolism (R2065); isoleucine degradation with propanoate secretion (R2281); citrate import via R831 (fwd); matrix MDH reverse (R533 rev) |
| 4 | 22 | 0.28080 | **0.077** | stop (< 0.15). The remaining KO growth is glycolysis to secreted pyruvate (R1310, then added as a part) |

- **Where the escape runs.** Each round, most of the KO's escape ran through enzymes the WT also requires,
  used at higher capacity. In round 1 (KO vs WT): PDC/NADP-ALD 163 vs 62, cytosolic MDH 217 vs 37, SDH 132 vs
  66, ICL 50 vs 21 (WT-required).
- **Tie-breaks.** Every candidate step tied at 0 % WT cost, so the tie-breaks recorded in TASK.md decided the
  step. I applied them during discovery; they are judgment.
- **Grouping.** The final 24 parts were grouped into 12 tier-1 groups (TASK.md, frozen table). Group 10 is a
  family of seven amino-acid/amine catabolic routes, kept as one group by the user's decision (Amendment 1).
- **The 0.50 % WT cost** in every round is the HMG-CoA part, closed whole on the user's instruction (R412 is
  WT-required, about 2.7).

Map: [map_discovery_rounds.png](map_discovery_rounds.png) shows the KO at rounds 0–4.

## 4. Step 2: combinatorial test (verified: [combinatorial.json](combinatorial.json); 4,096 subsets × WT/KO × A/V2)

| | A | V2 |
|---|---|---|
| nothing closed | KO/WT 0.826 | 0.827 |
| all 12 groups closed | 0.077, WT cost 0.50 % | 0.000, WT cost 0.50 % |
| hits at 15 % | 3,072 of 4,096 | 3,072 |
| **minimal cut sets at 15 %** | **{G3}, {G5}**, both KO/WT 0.0767, WT cost 0.00 % | **{G3}, {G5}**, KO/WT 0.000, WT cost 0.00 % |
| at 10 % | same two sets | same |
| at 1 % and 5 % | **none** (0.077 > 0.05) | same two sets |
| accepted (WT cost ≤ 2 %) | both | both |

- **3,072 = 2^12 − 2^10:** these are exactly the subsets that contain G3 or G5. No other group matters for the
  15 % call.
- **Leave-one-out from the full set:** reopening any one group (or any one of the 24 parts) leaves KO/WT at
  0.077 (A) or 0 (V2), because G3 and G5 back each other up.
- **Tier 2:**
  - G3 = {citrate import through R817 reverse, R831 forward} resolves to **R831 forward alone**. Closing R817
    reverse alone leaves 0.826.
  - G5 has one part (R533 reverse).
  - Re-checked in V2: both give 0.000.
- **No KO in any hit had "no solution".** All hits are optimal solves with low growth.

**Cutoff sensitivity.** In A the KO keeps 0.0216 h⁻¹ (KO/WT 0.077). That growth is limited by the lysine
supplied in the Task A medium (R1220 ≤ 0.00595, §5). So the call is essential at 10 % and 15 % but not at 1 %
or 5 %. In V2, with no lysine, it is essential at every cutoff.

Map: [map_complexI_cutsets.png](map_complexI_cutsets.png) shows WT, KO, KO + G5 and KO + G3, for A and V2.

## 5. Why one closure is enough (verified: [mechanism.json](mechanism.json); added after the tier-1 result, not pre-registered)

**Matrix NADH in the KO (FVA at the optimum).**
- R398 homoisocitrate DH is pinned at 2.41 (A) or 2.47 (V2): lysine synthesis.
- R533 is the only required consumer, running in reverse at [−7.8, −2.4].
- Otherwise only the matrix malic enzyme R538 is non-zero, as an optional producer ([0, 5.4]), plus a trace
  consumer R695 (0.001).
- In WT, the required matrix-NADH producers are PDH, NAD-IDH, R361, matrix MDH and R398. In the KO all except
  R398 drop to 0.

**Closing the outlet** (R533 reverse, or malate export by R831) leaves R398 with no way to reoxidise its NADH.
Lysine can then only come from the medium:

| | KO + G5 | + lysine uptake opened (R1220 lb −1000) |
|---|---:|---:|
| A | 0.0216 (lysine capped at 0.00595) | 0.243 |
| V2 | 0.000 | 0.227 |

Opening all 11 amino-acid rows gives growth near 22 h⁻¹. That is amino acids used as an unbounded carbon
source: not a meaningful rescue, reported only for completeness.

**Robustness.** G3 or G5 alone gives KO/WT ≤ 0.077 (A) and 0 (V2) with the overlay and the caps each on or
off (8 combinations).

**Other outlets are closed.** R573 (matrix-facing NADH:Q with the NDH2 gene YALI1F32476g) has bounds (0, 0) in
the model. R902 and R1386 (malate import) are irreversible. The matrix NADH kinase R572 is blocked.

**What this means (inferred from the above).**
- **Model level:** the model's complex I non-essentiality rests on two capabilities:
  - turning cytosolic NADPH into NADH and feeding it to external NDH2, where many routes substitute for each
    other;
  - exporting a small, obligate amount of matrix NADH as malate, which only R533 reverse plus R831 can do.
  The second one is the bottleneck.
- **Biology:** whether *Y. lipolytica* makes lysine-pathway NADH in the matrix depends on where homoisocitrate
  DH sits (§6).

## 6. Step 3: evidence for the parts in the minimal cut sets

**How the evidence was gathered.** Two read-only reviewer subagents opened the literature on 2026-10-14 and
saved every text they read in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/complexI_transhydrogenation_20261012/evidence/`.

**Tier legend.**
- **verified:** opened this session and the passage read. "(re-read)" marks the passages I read myself in the
  saved text; the others were read by the reviewer subagent.
- **abstract only:** only the abstract was read.
- **unread:** the source could not be opened.

The prompt's rule 3 questions are kept apart throughout: existence, then compartment and direction, then
substrate source.

### 6.1 Orthology (verified: [orthology_sw.py](orthology_sw.py) → [orthology_summary.json](orthology_summary.json))

**Method.**
- Task B's Smith-Waterman (BLOSUM62, gap 11/1; SHA-pinned `04a1f73f…`).
- Each gene is aligned against the *S. cerevisiae* proteins of its enzyme family only. This assigns the closest
  isoform; it is not a proteome-wide reciprocal-best-hit test.
- Sequences were fetched by RefSeq ([sequences_manifest_v3.json](sequences_manifest_v3.json)).

| Model reaction | Gene | Length (aa) | Best *S. cerevisiae* hit (identity) | Reciprocal in panel |
|---|---|---:|---|---|
| R533 matrix MDH | YALI1D20676g | 337 | Mdh1 (67 %); Mdh3 46 %, Mdh2 41 % | yes |
| R533, R534 | YALI1E17214g | 246 | Mdh1 (52 %) | no |
| R831 citrate transport | YALI1F33754g | 292 | Ctp1 (67 %) | yes |
| R817 citrate/2-OG | YALI1B14382g | 311 | Yhm2 (77 %) | yes |
| R398 homoisocitrate DH | YALI1D13238g | 385 | Lys12 (67 %) | yes |
| R549 NAD-MTHFD | YALI1B15928g | 325 | Mtd1 (64 %) | yes |
| R550 NADP-MTHFD | YALI1F38386g / YALI1E01589g | 1006 / 937 | Mis1 ≈ Ade3 (63–67 %) | cytosolic vs mitochondrial not resolved |
| R329 NAD-GDH | YALI1E11943g | 1048 | Gdh2 (43 %) | yes |
| R330 NADP-GDH | YALI1F23664g | 458 | Gdh3 ≈ Gdh1 (69 %) | yes |
| R1385 matrix NADP-IDH | YALI1F06197g | 431 | Idp1 = Idp2 (74.6 % each) | compartment not resolved |
| R487 NAD-IDH | YALI1E06015g, YALI1D08082g | 366, 369 | Idh1 (61 %), Idh2 (70 %) | yes |
| ALD (R1378, R134, R135, R136) | YALI1C04093g, YALI1E00588g | 496, 522 | Ald4 (60 %, 59 %) | Ald6-type not clear |
| ALD | YALI1D10197g, YALI1F06858g | 525, 515 | Ald2 (50 %), Ald5 (50 %) | weak |

### 6.2 The cut-set steps

**G5: matrix MDH R533, reductive direction (OAA + matrix NADH → malate).**
- **Existence:** two *Y. lipolytica* MDH genes are known.
  - Source: Kabran et al. 2012, DNA Res 19:231, [doi:10.1093/dnares/dss007](https://doi.org/10.1093/dnares/dss007);
    verified (re-read, full text).
  - YlMDH1 (YALI0D16753g) is mitochondrial by prediction only.
  - YlMDH2 (YALI0E14190g) makes the cytosolic and the peroxisomal forms by alternative splicing, shown with
    fluorescent tags.
  - NCBI notes say YALI1D20676g "Compare to YALI0D16753g" and that YALI1E17214g is a "truncated form of
    YALI0E14190g".
    Both mappings come from NCBI record notes opened by the reviewer (verified, reviewer).
- **What this means for the GPR (inferred):**
  - YALI1E17214g's OR membership in the matrix reaction R533 has no support.
  - Only YALI1D20676g (YlMDH1) should catalyse R533.
  - This does not change the cut set.
- ***S. cerevisiae* comparison:** Mdh1 is in the mitochondrial matrix (UniProt P17505, experimental; verified,
  re-read).
- **Direction and thermodynamics:** malate + NAD⁺ → OAA + NADH has ΔrG′° = +29.6 ± 0.6 kJ/mol at pH 7.0
  (+26.5 at pH 7.5), K′eq 6.6 × 10⁻⁶.
  - Source: eQuilibrator, KEGG R00342; verified (re-read the saved page).
  - So the reductive direction that G5 closes is the **thermodynamically favoured** one. The TCA direction runs
    only because OAA is kept very low.
  - No primary source measuring the in vivo direction of matrix MDH by compartment was found (reviewer).
- **Shuttle evidence in *S. cerevisiae*:**
  - Bakker et al. 2001, FEMS Microbiol Rev 25:15, [doi:10.1111/j.1574-6976.2001.tb00570.x](https://doi.org/10.1111/j.1574-6976.2001.tb00570.x)
    (abstract only). The enzymes and carriers for malate–OAA, malate–aspartate and malate–pyruvate shuttles
    exist, but there is no consistent evidence of their in vivo activity, and the malate–citrate shuttle "can be
    ruled out" for lack of key enzymes or transporters (verified, re-read). Only the ethanol–acetaldehyde shuttle
    has in vivo support (Bakker et al. 2000, abstract only).
- **Conclusion:** no evidence supports closing R533 reverse as a **direction**.

**G3: R831 citrate/malate antiport, citrate in and malate out (tier 2: R831 alone).**
- **Existence:** YALI1F33754g is the Ctp1 ortholog (67 %; reciprocal here and in Task C's carrier-family test).
- **Transport mode** (*S. cerevisiae*, cross-species):
  - Aluvila et al. 2010, J Biol Chem 285:27314, [doi:10.1074/jbc.M110.137364](https://doi.org/10.1074/jbc.M110.137364),
    verified (re-read the saved full text).
    - Reconstituted Ctp1 (the Cys-less control, stated to match wild type) is highly specific for internal
      citrate. Internal isocitrate gives 18 % of the citrate exchange, and **other internal anions, malate
      included, do not support exchange**.
    - Single-residue mutants (e.g. R189C) gain malate exchange.
  - Kaplan et al. 1995 (verified by the reviewer, Table I not read): reconstituted Ctp1 does citrate/citrate
    exchange, and external dicarboxylates are ineffective inhibitors.
  - **Counterpoint (second-hand):** Kaplan 1995's introduction cites Evans et al. 1983 (intact *S. cerevisiae*
    mitochondria). There, malate serves as an alternative substrate "to a considerably lesser extent"
    (verified, re-read in Task C's saved text). Evans 1983 itself is unread.
- **What this means (inferred):** the citrate/**malate** exchange mode that the model gives R831, in both
  directions, has weak support in *S. cerevisiae*: no exchange in reconstitution, a weak substrate in intact
  mitochondria, and no malate–citrate shuttle per Bakker 2001. Task C had already flagged that this mode rests
  on homology and "by similarity". No *Yarrowia* reconstitution exists.
- **Classification:**
  - **limits a claim** here: G3's closure has better support than a pure direction closure;
  - **affects another line**: Task C's elastic fit routes malate import through R831.
- **R817 (Yhm2-type citrate/2-OG carrier):** it is not needed in any cut set (tier 2). Evidence as in Task C
  (Castegna 2010, saved in Task C's evidence folder, not re-read here).

### 6.3 Complex I and matrix-to-cytosol shuttles in *Y. lipolytica*

All sources here are about *Y. lipolytica*.
- **Kerscher et al. 2001,** J Cell Sci 114:3915, [doi:10.1242/jcs.114.21.3915](https://doi.org/10.1242/jcs.114.21.3915)
  (abstract; verified, re-read; full text blocked):
  - NDH2 is external and the only alternative NADH dehydrogenase;
  - complex I defects are lethal;
  - matrix-targeted NDH2 rescues them, the external version does not;
  - the authors conclude that matrix-to-cytosol shuttles are insufficient.
  - The growth medium is unread.
- **Garofano et al. 2006,** Yeast 23:1129, [doi:10.1002/yea.1426](https://doi.org/10.1002/yea.1426) (abstract;
  verified, re-read):
  - complex I loss-of-function mutants are obtainable only with internal NDH2 (NDH2i);
  - with complex I inhibitors and NDH2i under the ICL1 promoter (expressed on acetate, repressed on glucose),
    cells grew on acetate **but not on glucose minimal medium**.
  - This is the closest evidence that complex I is required on glucose minimal medium.
- **NDH2i strains in later papers:** Kerscher et al. 2002 (abstract) and later complex I papers (Wydro & Balk
  2013; Dröse 2011; Maclean 2018; Schiller 2022; all reviewer-opened) use NDH2i strains. YPD is stated only for
  Wydro 2013, Dröse 2011 and Maclean 2018.
- **No viable complex I deletion without NDH2i was found** (reviewer).
- **Mapping to the model (inferred):**
  - The model has external NDH2 (R570), and its matrix-facing version R573 is fixed at (0, 0). That matches
    these papers.
  - The model's escape relies on a shuttle (R533 reverse + R831) that these papers say is insufficient in vivo.
  - The papers give a qualitative capacity statement, not a rate.

### 6.4 Lysine pathway (R398, the obligate matrix-NADH source)

- ***S. cerevisiae*:**
  - Lys12 (homoisocitrate DH) is mitochondrial: UniProt P40495 and SGD GO list four high-throughput datasets
    (verified, re-read), and the YeastGFP (Huh 2003) call (reviewer). The evidence is high-throughput (HDA),
    not a dedicated localization study.
  - Lys4 is also mitochondrial (high-throughput).
  - SGD also carries a manually curated "is active in cytosol" annotation for Lys12 and Lys4. It is evidence
    code RCA, from SGD's pathway database (verified, re-read).
  - The homocitrate synthases Lys20/21 are nuclear (Chen et al. 1997, abstract).
- ***Y. lipolytica*:**
  - no localization study of homoisocitrate DH (YALI0D10593g) was found;
  - the model gene YALI1D13238g is a "truncated form of YALI0D10593g" per NCBI (reviewer);
  - LYS1, homocitrate synthase, has a predicted mitochondrial targeting signal (Perez-Campo et al. 1996,
    abstract).
- **So the claim that lysine synthesis is the obligate matrix-NADH source is model-level.** The model's matrix
  placement of R398 matches *S. cerevisiae*, but is **not verified for *Y. lipolytica***.

### 6.5 Parts not in any cut set (the prompt's named enzymes; reviewer-opened, abstract level unless stated)

- **MTD1:**
  - NAD-specific and cytoplasmic (Barlowe & Appling 1990, [doi:10.1021/bi00482a020](https://doi.org/10.1021/bi00482a020)).
  - In vivo it serves oxidative one-carbon flux and "does not participate significantly" in the reductive
    direction (West et al. 1996, [doi:10.1021/bi952713d](https://doi.org/10.1021/bi952713d)).
  - This agrees with the model's irreversible oxidative R549.
- **ADE3/MIS1:**
  - ADE3 is cytosolic and MIS1 mitochondrial (Shannon & Rabinowitz 1986, 1988).
  - Formate is reduced to serine through cytosolic C1-THF synthase (Pasternack et al. 1992, 1994).
  - So R550 reverse, the reductive NADP direction, is physiological.
  - A MTD1-forward + ADE3-reverse cycle is therefore made of two physiological directions. Limiting it is a
    capacity question, not a direction question.
- **GDH:**
  - GDH1/GDH3 are NADP and biosynthetic; GDH2 is NAD and catabolic, low with ammonium/glutamine and raised on
    non-fermentable carbon (Miller & Magasanik 1990/1991; Coschigano 1991).
  - A GDH1/GDH2 NADPH → NADH cycle was proposed from suppression by GDH2 overexpression in pgi1Δ (Boles et al.
    1993,
    [doi:10.1111/j.1432-1033.1993.tb18266.x](https://doi.org/10.1111/j.1432-1033.1993.tb18266.x)).
  - *Y. lipolytica*: ylGDH1 (YALI0F17820g, NADP) is needed for N assimilation; ylGDH2 (YALI0E09603g, NAD) for
    glutamate use; they are not interchangeable (Trotter et al. 2020, [doi:10.1002/yea.3425](https://doi.org/10.1002/yea.3425)).
    The mapping to YALI1F23664g / YALI1E11943g comes from NCBI notes.
- **ALD:** ALD6 is cytosolic and NADP; ALD4/ALD5 are mitochondrial (Meaden 1997; Tessier 1998; Wang 1998; UniProt
  verified). The four *Y. lipolytica* genes in the shared ALD GPR are closest to Ald4/Ald2/Ald5, not clearly
  Ald6 (§6.1).
- **IDH:**
  - IDP1 is mitochondrial and IDP2 cytosolic (Contreras-Shannon & McAlister-Henn 2004).
  - IDP2 shows similar Km in both directions in vitro (2005).
  - **No source showing in vivo reductive carboxylation by IDP1** was found. R1385 reverse is not in a cut set.
- **Conflicts in the literature:** idp1Δ phenotype (Haselbeck 1991 vs 1993); whether ADE3's catalysis is needed
  for purines (Barlowe & Appling 1990 vs Song & Rabinowitz 1993); GDH2 localization (SGD cytosol vs a review's
  "mitochondrial").

### 6.6 Erg13 comparison (Task C handoff)

- ***S. cerevisiae* Erg13:** nuclear in YeastGFP (Huh 2003) and cytosolic only by computational annotation. Its
  "mitochondrion" GO term is NAS:UniProt (verified, re-read). Task C traced it to a review (Task C REPORT §5.4).
- No yeast with a mitochondrial HMG-CoA synthase was found.
- This *S. cerevisiae* comparison agrees with the user's note. It is not *Yarrowia* evidence. The HMG-CoA part
  (R412 + R1172) is in no cut set.

## 7. Step 4: what complex I essentiality needs (decision is the user's)

**What the accepted cut sets do.** Both close one direction of a real, orthology-supported protein (§6):
- matrix MDH R533 (YALI1D20676g, Mdh1 ortholog) in the reductive direction;
- the Ctp1-type carrier R831 in the citrate-in, malate-out direction.

**What the evidence says about each closure** (§6.2):
- **MDH:** no direction evidence supports closing R533 reverse. The reductive direction is the thermodynamically
  favoured one (eQuilibrator estimate).
- **R831:** the closure has some cross-species support. The citrate/malate mode itself is doubtful in
  *S. cerevisiae* (Aluvila 2010; Bakker 2001), with one older second-hand counterpoint (Evans 1983 via
  Kaplan 1995). No *Yarrowia* data exist.
- **Conclusion (inferred):** neither closure is a direction change backed by enzyme or thermodynamic evidence for
  *Y. lipolytica*. Under the prompt's rule, complex I essentiality therefore **needs a capacity or regulatory
  assumption**.

**The *Yarrowia* evidence is about capacity.** Kerscher et al. 2001 found that complex I deficiency is lethal,
and that it is rescued by NDH2 redirected to the matrix face but not by the native external NDH2. They conclude
that matrix-to-cytosol shuttles are insufficient (§6). Garofano 2006 shows the same requirement on glucose minimal medium. Both support a limit on matrix NADH
export, but give no number.

**Capacity scan** (verified: [shuttle_cap.json](shuttle_cap.json) and the finer grid
[shuttle_cap_fine.json](shuttle_cap_fine.json), run after the audit; cap on R533 reverse, per 100 glucose):

| Cap on R533 reverse | A, no floors | V2, no floors | A, PDH+CS+G6PDH + R825 ASSUMPTION | V2, same |
|---|---:|---:|---|---|
| 0 | 0.077 | 0.000 | no solution | no solution |
| 0.2357 | 0.15003 (just above cutoff) | 0.079 | no solution | no solution |
| 0.236 | **0.150** | 0.079 | no solution | no solution |
| 0.449 | 0.216 | **0.1501** | no solution | no solution |
| 0.45 | 0.217 | **0.150** | no solution | no solution |
| 0.5 | 0.232 | 0.167 | no solution | no solution |
| 1 | 0.388 | 0.334 | no solution | no solution |
| 2 | 0.699 | 0.669 | no solution | no solution |
| 2.5 | 0.826 | 0.827 | no solution | no solution |
| 80 – 86.66 | 0.826 | 0.827 | no solution | no solution |
| 86.67 | 0.826 | 0.827 | **0.078** | **0.002** |
| 86.7 | 0.826 | 0.827 | **0.088** | **0.012** |
| 86.8 | 0.826 | 0.827 | **0.118** | **0.045** |
| 86.9 | 0.826 | 0.827 | **0.149** | **0.078** |
| 87 | 0.826 | 0.827 | 0.180 | **0.111** |
| 87.1 | 0.826 | 0.827 | 0.211 | **0.144** |
| 87.2 | 0.826 | 0.827 | 0.242 | 0.177 |
| 88 | 0.826 | 0.827 | 0.319 | 0.278 |
| 100 | 0.826 | 0.827 | 0.342 | 0.299 |
| 120 or none | 0.826 | 0.827 | 0.361 | 0.316 |

With floors, every cap ≤ 86.66 gives "no solution". That is a constraint conflict (the floors force PDH flux
the KO cannot reoxidise): **unresolved, not lethal**. WT cost is 0 at every cap. The first grid jumped from 80
to 100. The first audit found the ~86.7 threshold, and the fine grid confirms it.

The second audit (2026-10-14) found two problems, both now fixed:
- `shuttle_cap.py` labelled no-solution KOs `essential_15 = True`, against benchmark.md. They are now
  `"unresolved"`.
- The band was overstated.

Both CSV/JSON pairs were regenerated, and the fine grid gained the boundary points 0.2357, 0.449 and
86.66–87.2. Every value at the earlier grid points reproduced exactly; only the 30 no-solution labels changed.

**Reading.**
- **Without the 13C floors:** essentiality at 15 % needs the shuttle capped below about 0.24 (A) or 0.45 (V2)
  per 100 glucose, i.e. 0.006–0.011 mmol gDW⁻¹ h⁻¹. That is in effect closed.
- **With the floors:** the KO has no solution below 86.66, which is the PDH-floor NADH (86.663 per 100 glucose,
  audit).
  - It is essential at 15 % only from 86.67 to about 86.9 (A; 0.149 at 86.9, 0.180 at 87.0), a band about 0.24
    wide.
  - For V2 the band runs to about 87.1 (0.144 at 87.1, 0.177 at 87.2), about 0.46 wide.
  - KO/WT is ≥ 0.28 from 88.
  - Any essential call from this combination sits on a knife edge.

**Options.** **User decision (2026-10-14): option (d).** No ASSUMPTION file is written (TASK.md Amendment 2).

| Option | What it is | Effect (verified) | Basis and limits |
|---|---|---|---|
| (a) ASSUMPTION file "matrix MDH oxidative only" | R533 lb 0 (forward kept), in memory; never in the XML | complex I KO/WT 0.077 (A; essential at 10/15 %) and 0 (V2; all cutoffs); WT cost 0 | Kerscher 2001 (qualitative, *Yarrowia*). Against it: thermodynamics favours the closed direction. Its effect on the other benchmark calls is **untested** and should be screened before any use |
| (b) ASSUMPTION file "R831 malate export closed" | R831 ub 0 (citrate in / malate out) | same as (a) | better supported than (a): reconstituted *S. cerevisiae* Ctp1 does not exchange citrate for malate (Aluvila 2010, §6.2), and Kerscher 2001 supports a limited shuttle. Limits: cross-species evidence; benchmark side effects untested |
| (c) 13C floors + a shuttle cap | condition pair | below 86.66 the KO has no solution; essential only in a band ~0.24 (A) / ~0.46 (V2) wide | a call that hinges on an unmeasured number; diagnostic only |
| (d) no assumption now | record that the model needs a shuttle assumption to reproduce complex I lethality | none | the evidence-level default |

**Recommendation (inferred).**
- **(d) now, plus a separate evidence review of R831's citrate/malate mode.** That mode affects both directions
  and Task C's 13C fit, so it is a curation question with cross-line effects, not a Task D diagnostic.
- **If you want a diagnostic variant before that review:** (b) rather than (a), labelled ASSUMPTION like the
  R825 file, and only after a full benchmark side-effect screen.
- **Do not curate R533's direction in the XML:** thermodynamics favours the reverse direction, and no direction
  evidence supports closing it.

## 8. Step 5: 13C ratio floors (verified: [ratios.json](ratios.json), [ratios.csv](ratios.csv); ASSUMPTION overlay)

**Floors.** Each floor is the 95 % lower bound of the pinned MFA × glucose uptake (verified from the file):

| Flux | MFA row | Floor (fraction of glucose uptake) |
|---|---|---:|
| PDH | R19 | 0.866629 |
| CS | R20 | 0.2817 |
| SDH | R26 net | 0.150997 |
| G6PDH | R11 | 0.502982 |

All 16 subsets were run, with WT and KO under the same floors.

| Setting | Accepted sets (WT cost ≤ 2 %) | Lowest KO/WT among accepted | KO's required R533 reverse (FVA lower bound) |
|---|---|---|---|
| A, R825 as supplied | none, CS, SDH, CS+SDH. PDH costs 2.8 % (KO no solution); G6PDH 2.75 % | 0.646 (CS) | 15 (CS) |
| V2, R825 as supplied | none, CS, SDH, CS+SDH. PDH costs 3.5 % (KO no solution); G6PDH 2.9 % | 0.637 (CS) | 14 (CS) |
| A, **R825 ASSUMPTION** | all 16 (G6PDH 0.99 %) | 0.361 (PDH+CS+G6PDH) | 116 |
| V2, **R825 ASSUMPTION** | all 16 (G6PDH 0.89 %) | 0.316 (PDH+CS+G6PDH) | 115 |

- **PDH alone** with the R825 ASSUMPTION gives 0.687 (A) and 0.680 (V2). Task C reported 0.680 for V2 (verified
  match). The prompt's recalled 0.709 is not reproduced.
- **The R825 ASSUMPTION raises WT growth** by 1.9 % (A) and 2.2 % (V2) against R825 as supplied.
- **No ratio set makes complex I essential.** The floors force PDH/TCA flux in the KO, and that matrix NADH
  leaves through the same R533-reverse/R831 shuttle at 89–116 per 100 glucose at the KO optimum (FVA), about
  half of WT complex I flux.

Map: [map_13c_floors.png](map_13c_floors.png).

## 9. Maps

- [map_complexI_cutsets.png](map_complexI_cutsets.png): A and V2 × WT, KO, KO + G5, KO + G3.
- [map_13c_floors.png](map_13c_floors.png): WT and KO under PDH + CS + G6PDH with the R825 ASSUMPTION.
- [map_discovery_rounds.png](map_discovery_rounds.png): the KO at discovery rounds 0–4.

Conventions:
- Values are per 100 mol glucose actually taken up in that panel ([map_edges.csv](map_edges.csv)), from one
  pFBA optimum. The uptake is printed in each info box; three A panels take up less than the 2.414 bound.
- Dashed grey edges are closed.
- Nodes: G6P, oxPPP, NADPH (cyt), MTHFD pair, GDH pair, NADH (cyt), NDH2, pyruvate, PDH, TCA, NADH (mito),
  complex I, CoQ, plus the cut-set route (matrix malate, R533, R831) and lysine synthesis (R398).
- Judgments rest on the FVA in the JSON files, not on the drawn vectors.

## 10. Genes named (AGENTS.md §4)

| Gene | Verified name | Function | Evidence level |
|---|---|---|---|
| YALI1D20676g (NCBI note: "Compare to YALI0D16753g") | **YlMDH1** (Kabran 2012, for YALI0D16753g) | mitochondrial malate dehydrogenase (R533) | localization by prediction only; Mdh1 ortholog (67 %, this task) |
| YALI1E17214g (NCBI note: truncated form of YALI0E14190g) | **YlMDH2** (Kabran 2012, for YALI0E14190g) | cytosolic and peroxisomal MDH (splice variants) | fluorescent-tag localization (Kabran 2012); model gene is a 246-aa truncated model. In R533's GPR without support |
| YALI1F33754g | none verified | tricarboxylate (citrate) carrier, Ctp1 ortholog (R831) | homology (67 %, reciprocal). The malate mode is not supported by Sc reconstitution (Aluvila 2010) |
| YALI1B14382g | none verified | citrate/2-oxoglutarate carrier, Yhm2 ortholog (R817) | homology (77 %) |
| YALI1D13238g (NCBI note: truncated form of YALI0D10593g) | none verified | homoisocitrate dehydrogenase, Lys12 ortholog (R398) | homology (67 %); Sc Lys12 mitochondrial by high-throughput data only |
| YALI1B15928g | none verified | NAD-dependent methylenetetrahydrofolate DH, Mtd1 ortholog (R549) | homology (64 %) |
| YALI1E11943g (NCBI note: "Compare to YALI0E09603g") | **ylGDH2** (Trotter 2020, abstract) | NAD-glutamate dehydrogenase (R329) | *Yarrowia* deletion phenotype (abstract) + homology |
| YALI1F23664g (NCBI note: "Compare to YALI0F17820g") | **ylGDH1** (Trotter 2020, abstract) | NADP-glutamate dehydrogenase (R330) | *Yarrowia* deletion phenotype (abstract) + homology |
| YALI1F32476g | NDH2 (Kerscher 1999/2001, for the *Y. lipolytica* gene; ID mapping not checked here) | external alternative NADH dehydrogenase (R570; R573 fixed at 0) | *Yarrowia* deletion and retargeting (abstracts) |

A GPR or model assignment is not experimental validation. No essentiality call changed in this task, so no
benchmark gene list is needed.

## 11. What is proven, observed, inferred or untested

**Provable / static.**
- 3,072 = 2^12 − 2^10. Every subset containing G3 or G5 is a hit and no other subset is (checked over all 8,192
  solves). So the minimal cut sets are exactly {G3} and {G5}.
- Observed by FVA, not a static proof: R533 reverse is the only required matrix-NADH sink besides complex I in
  the KO (R695 is a trace consumer). R573 is fixed at (0, 0) and R572 is blocked.

**Observed in runs (verified this session, SHA-pinned inputs):** every number in §2–§5, §7 and §8.

**Inferred:**
- the tie-breaks and route names in discovery;
- the reading that lysine biosynthesis is the obligate matrix-NADH source (it rests on the FVA and the lysine
  rescue, verified, plus the model's compartment for R398);
- the recommendation in §7.

**Not done / untested:**
- a benchmark side-effect screen for the §7 options;
- glycerol and oleate;
- the as-supplied model without overlay and caps (except §5's robustness rows);
- proteome-wide reciprocal best hits (§6 uses family panels);
- a BLAST cross-check of the Smith-Waterman implementation;
- localization predictors for the two trifunctional C1-THF synthases, and for YlMDH1 and homoisocitrate DH;
- full texts that could not be opened: Kerscher 2001 (growth medium), Bakker 2001, Sickmann 2003;
- the GPR issue in another line: YALI1E17214g (truncated YlMDH2) sits in the matrix MDH GPR (R533). The cytosolic
  MDH R534 depends only on this truncated gene model.

**Possible improvements, each testable.**
- **Make the cut-set search exact.** A bilevel (OptKnock-style) MILP over every direction of every reaction that
  touches matrix NADH would find the minimal cut sets without hand-built discovery rounds.
  - Test: it should return {R533 rev} and {R831 fwd} as the size-1 sets.
- **Screen the side effects.** Run the full gene benchmark with R533 lb 0 (and separately R831 ub 0) at all four
  cutoffs, and list every changed call.
- **Measure the capacity, don't guess it.** The cap scan shows the decision is binary near zero without floors.
  A measured or literature value for matrix-to-cytosol shuttle capacity would turn §7 (c) into a real test.

## 12. Reproduce

Run from the repository root with:
- `IYALI26_RESEARCH_ROOT=/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`;
- `PY=<main checkout>/.venv/bin/python` for LPs;
- the conda `python` (matplotlib) for the maps.

Each LP script refuses to overwrite its outputs. `discovery.py <r>` reproduces round r from the final
`parts_definition.json`, because parts are applied by `round_added`.

```
$PY -B model/reports/complexI_transhydrogenation_20261012/escape.py
$PY -B model/reports/complexI_transhydrogenation_20261012/parts.py
$PY -B model/reports/complexI_transhydrogenation_20261012/discovery.py 1      # then 2, 3, 4
$PY -B model/reports/complexI_transhydrogenation_20261012/step_costs.py 1 R511:fwd,R328:fwd,...   # step lists in step_costs_round*.json
$PY -B model/reports/complexI_transhydrogenation_20261012/combinatorial.py
$PY -B model/reports/complexI_transhydrogenation_20261012/mechanism.py
$PY -B model/reports/complexI_transhydrogenation_20261012/ratios.py
$PY -B model/reports/complexI_transhydrogenation_20261012/shuttle_cap.py
$PY -B model/reports/complexI_transhydrogenation_20261012/shuttle_cap.py fine
$PY -B model/reports/complexI_transhydrogenation_20261012/capacity_costs.py
$PY -B model/reports/complexI_transhydrogenation_20261012/fetch_sequences.py   # network; writes *_v3 files
$PY -B model/reports/complexI_transhydrogenation_20261012/orthology_sw.py
$PY -B model/reports/complexI_transhydrogenation_20261012/map_vectors.py
python -B model/reports/complexI_transhydrogenation_20261012/make_maps.py
```

## 13. Deviations and disclosures

- **Pre-plan exploration.** About 1,050 LPs ran in the session scratchpad on Task C's helpers before the plan.
  None of those numbers is used.
- **Discovery tie-breaks.** They were applied during discovery (TASK.md). The TASK.md rule only fixed "smallest
  WT cost", and every candidate tied at 0.
- **Discovery missed the decisive step.** R533 reverse and R831 forward were already flagged as WT-free steps in
  round 1 (KO −19.6 and 21.1). I did not pick them as routes until round 3, because I followed the large fluxes.
  The 2^12 test is what showed they decide the outcome alone. The discovery rounds still bound the escape space;
  they were not needed for the cut sets.
- **Added beyond the pre-registration:** `mechanism.py` (§5) and `shuttle_cap.py` (§7). Both are diagnostics of
  results already obtained.
- **Sequence fetch.** It ran three times, and all outputs are kept:
  - v1: the locus tag was built wrongly, and gene search matched synonyms for ALD3/4/5;
  - v2: locus-tag hits returned both a RefSeq and a GenBank record, and none was chosen;
  - v3 is used.
- **Large file.** `combinatorial_tier1.csv` (1.9 MB) was written into this folder, against the ~1 MB rule. A
  byte-identical copy (`fcd32946…`) is in the research workspace, under
  `artifacts/tasks/complexI_transhydrogenation_20261012/`. The repository copy was removed before commit on
  2026-10-14 with the user's authorization, after its SHA-256 was checked against the workspace copy.
- **Maps re-drawn three times.** The PNGs and `map_edges.csv` are drawing outputs and were overwritten each
  time:
  - the first layout covered nodes;
  - the second divided by the 2.414 uptake bound, so three A panels with lower uptake were mis-scaled (found by
    the audit);
  - the third normalises by actual uptake.

  Some edge labels still overlap (e.g. R329 in the A KO panel); the values are in `map_edges.csv`.
- **Found by the audit and fixed after it** (§14):
  - the capacity-step closure costs required by TASK.md had not been run;
  - the shuttle-cap thresholds were misread from a coarse grid (80 → 100);
  - several evidence sentences overstated or omitted their sources.

  `capacity_costs.py` and `shuttle_cap.py fine` were run after the audit. They were not re-audited.
- **No platform code changed,** so `pytest` was not run.

## 14. Independent audit

**Reviewer and scope.**
- A fresh read-only reviewer audited the work on 2026-10-14, with no conversation history and only the files.
- It wrote nothing to the repository or the research workspace. Its scratch scripts ran in a temporary
  directory, with 164 LPs of its own.

**Overall:** the computational results are SUPPORTED. The corrections it listed were applied as described in
§13.

| Claim | Verdict |
|---|---|
| A1 input SHAs; only the task folder is new; no XML, condition, curation or evidence file changed | SUPPORTED |
| A2 step-1 growth, FVA and pair interchangeability (re-run) | SUPPORTED |
| A3 scripts implement TASK.md; frozen table = JSON | PARTIALLY SUPPORTED: capacity-step closure costs were missing → now run (§3) |
| A4 cut sets {G3}, {G5}; 3,072 hits = sets containing G3 or G5 (re-run) | SUPPORTED |
| A5 tier 2: R817 rev alone 0.826; R831 fwd alone 0.077 | SUPPORTED |
| A6 mechanism, lysine rescue, robustness (re-run) | SUPPORTED; wording on R538/R695 corrected |
| A7 step-5 floors and acceptance (re-run) | SUPPORTED |
| A8 shuttle-cap cells | PARTIALLY SUPPORTED: the reading of the thresholds was wrong → fine grid run, §7 rewritten |
| A9 orthology pins and scores (own SW re-implementation) | SUPPORTED |
| A10 maps | SUPPORTED with caveats: normalisation by the uptake bound → fixed; label overlaps remain |
| A11 evidence tiers vs saved texts | PARTIALLY SUPPORTED → corrected (Bakker 2001 and Evans/Kaplan added; YPD scope; Boles wording; Erg13; NCBI "Compare to"; SGD cytosol annotation) |
| A12 disclosures | PARTIALLY SUPPORTED → completed in §13 |
