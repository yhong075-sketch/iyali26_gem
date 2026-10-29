# TASK: UM2, the inherited vacuolar closures (2026-10-28)

## 0. Decision

- **Authorization (2026-10-28):** after the Stage 1 draft below, the user answered **"Option A0"**: keep
  every UM2 closure and document them. No build, no curated data, no model change.
- **Scope of Stage 2:** this folder (TASK.md, REPORT.md, `stage1/`, `scripts/`, an Escher map), a note on the
  UM2 row of the census report, and a `model/STATE.md` entry. Commit and push only on a further instruction.
- **Baseline update:** UM7 was merged (PR #30, `9d559d4`) while Stage 1 ran. The newest candidate on `main`
  is now `5f54eb22…` (`…_dolichol_route_v2_um7_nuclear_diffusion_20261028.xml`). It differs from `9d655c2b`
  only by `R_ADPRtn` and `R_PPItn`. All 69 vacuolar and UM2 reactions are identical in the two files
  (bounds, GPR, stoichiometry, subsystem, notes), and so are the vacuolar species (checked 2026-10-28).
  The as-built, C′ and "open all 17" checks reproduce on `5f54eb22` (`stage1/stage1_um2_5f54eb22.log`).

The Stage 1 draft follows unchanged except for file paths. Stage 1 outputs are in `stage1/` (logs, small
JSON/TSV, BLAST tables, evidence notes) and `scripts/`; the full JSON, the BLAST databases and the source
files are in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/um2_vacuolar_closures_20261028/stage1/`.

**Recommendation: Option A** (keep every closure; record UM2 as a deliberate, documented closure).

- **History:** the 17 closures were made between the Yeast 7.11 template (open) and iYali4 (2016, closed) and
  inherited unchanged since; no source records who closed them or why (section 2).
- **Biology:** the vacuolar V-ATPase and the AVT-family amino-acid transporters are real. W29 has
  reciprocal-best-hit orthologs of AVT1, AVT3, AVT6 and ATG22, and isolated *Y. lipolytica* vacuoles pump
  protons (section 5). So the closures do not mean "absent in the organism".
- **Mechanism:** this steady-state model has no vacuolar storage, pH or membrane-potential term, so there is
  nothing for these reactions to do:
  - On the Task A medium, every option frees reactions only as ATP-burning futile cycles. pFBA flux is 0,
    WT growth is unchanged to 1e-14, and no single-gene knockout changes by more than 1e-6. There is no
    class change at 1, 5, 10 or 15 %, and no energy-generating cycle.
  - R1163 opened as written (cytosol → vacuole) creates an unbounded no-uptake proton loop with the
    amino-acid transports. That is the same loop class Yeast-GEM shows, and a plausible (inferred,
    unrecorded) motive for the 2016 closure. It still holds for R1163 today.
- **Why A:** opening B, B′ or C′ is biologically defensible and harmless, but it adds no predictive content.
  Its only measurable effect is a lower census blocked count (Task A 37 → 25 at most), and that is not a
  reason to change the model (`AGENTS.md`). The V-ATPase false negatives come from a missing acidification
  requirement, not from the closure (section 5). A keeps the model as is and records why.
- **Not part of UM2, listed as follow-ups:**
  - R925 is unbalanced (−3 H).
  - The dipeptide sets are an iYli21 addition that cannot work: exchanges allow secretion only, the
    hydrolases have no gene and sit in the vacuole. The evidence points to PTR2 uptake plus cytosolic
    hydrolysis instead.
  - The open His/Lys vacuolar transports R2055/R2056 miss the VBA1 reciprocal best hit `YALI1B04391g`.
- **If you prefer the biology written into the model,** Option C′ is the version without loops: R795,
  R1363, R1163 reversed to a vacuole → cytosol leak, plus the 11 AVT/ATG22-backed amino-acid transports.
  It leaves out the polyamines.

## 1. Starting state (read 2026-10-28)

- **Worktree / HEAD:** `claude/fervent-shannon-141b92` at `0093c89` = `origin/main` after `git fetch`
  (0 ahead / 0 behind), clean.
- **UM7 is not merged.** Its branch `claude/hungry-mendel-f9f7df` (candidate `5f54eb22`, commits `3720370`,
  `8add118`) exists locally but is not on `origin/main`. Per the brief, the baseline is the newest merged
  candidate.
- **Model executed:**
  `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_20261028.xml`,
  SHA-256 `9d655c2bdbf82721a38f148a9a4ee0a385bb509ea7b591111107113c4aa8335c` (re-hashed; run from a
  byte-identical short-path copy because cobra rejects str paths of 260+ characters on macOS). Not a
  designated reference (`model/candidates/README.md`). The census itself ran on `a1237821`; none of the
  later steps (cofactor terms, R2135, R1726, dolichol) touches a vacuolar reaction or species (inventory
  below is read from `9d655c2b`).
- **Condition:** Task A medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) +
  PO1f profile `state/strain_profiles/po1f_sd_leu.json` (`35307853…`, uracil 1000, R612 closed);
  experimental list `consensus_essential_genes.csv` (`1e887f5a…`). Gurobi 13, FeasibilityTol and
  OptimalityTol 1e-9, real HOME. "Open condition" = raw model with every exchange (−1000, 1000), the census's
  "any medium".
- **Instruction files read:** `~/.claude/CLAUDE.md`, worktree `AGENTS.md`, `model/expected/benchmark.md`,
  `model/STATE.md` entries of 2026-09-10/11 (R795 open diagnostics, R795/R1363 double-open, the separate
  closure-origin task) and 2026-10-26 to 2026-10-28, plus the 2026-09-14 to 09-28 vacuole/dipeptide entries
  they lead to.

## 2. Recorded history of the closures

| Source | V-ATPase | H+ diffusion | water diffusion | aa / polyamine transports | Read how |
| --- | --- | --- | --- | --- | --- |
| Yeast 7.11 (2014, *S. cerevisiae* template) | `r_1086` [0, inf] | open | `r_2099` [−inf, inf] | open | 2026-09-11 provenance task (S1), file read |
| iYali4, Kerkhoven et al. 2016, SI7 (1,985 rxns) | `R_1086` [0, 0] | [0, 0] | `R_2099` [0, 0] | 24 of 25 `C_va` reactions [0, 0]; only Ile import `R_1203` open | provenance task (S2) |
| iYali4_corr, Xu et al. 2020 | [0, 0] | [0, 0] | [0, 0] | same | provenance task (S3) |
| iYli21, Guo et al. 2022, `mmc6.xml` | R795 [0, 0] | R1163 [0, 0] | R1363 [0, 0] | same 24 closed; 12 vacuolar reactions added (the dipeptide sets, R2056) | provenance task (S4); byte-identical to this repo's first import `05ae6e1:data/iyli21.xml` |
| iYali v4.1.2 (Chalmers rebuild from Yeast 7.6, file dated 2021-04-04) | `y001086` [0, 1000] | `y001831` [−1000, 1000] | `y002099` [−1000, 1000] | all open | read this session from `reference/external_models/chalmers_yarrowia/eciYali5-GEM/model/iYali5-GEM/sourceModels/iYali4.xml` |
| iYali5 (2025, built from iYali4_corr by `code/debiaggi2025/generate_iYali5GEM.m`) | `1086` [0, 0] | `1831` [0, 0] | `2099` [0, 0] | closed except Ile import `1203` | read this session; the script reopens several reactions "blocked for no apparent reason" (lines 164–197) but none in the vacuole |
| Yeast-GEM (current, `9fd2c572…`) | `r_1086` [0, 1000] | `r_1831` [−1000, 1000] | `r_2099` [−1000, 1000] | open | read this session |
| This repository | R795 [0, 0] since the first import (2026-03-19); no project step ever changed a UM2 bound in a built model | | | | `reports/r795_open_diagnostic_20260911/HISTORY.md` (45 revisions, 23 XML blobs, 22 with R795, all [0, 0]); builder grep this session |

**What is recorded:** the closures were introduced between the Yeast 7.11 template and the published iYali4
(2016) and have been carried unchanged into iYali4_corr, iYli21, iYali26 and iYali5. **What is not recorded:**
who closed them and why. The 2026-09-11 provenance task (Codex worktree
`~/.codex/worktrees/1df8/iyali26_gem/artifacts/r795_upstream_provenance_20260911/REPORT.md`, independent audit
8 of 9 claim groups supported, 1 open = the negative search) found no author statement; the Yeast 7.11 and
iYali 4.0.2 energy-repair notes concern mitochondrial reactions, not these.

**A plausible motive, inferred (not recorded):** in Yeast-GEM, with every uptake closed, the reversible
vacuolar H+ diffusion `r_1831` plus the H+-coupled vacuolar amino-acid import/export pairs carry flux up to
1000 (internal loops; `yeastgem_loops.log`, verified this session). Closing `r_1831` alone does not remove
them (the uniport/antiport pairs still loop). The same structure was in Yeast 7.11; closing the vacuolar set
removes it. This fits the 2016 pattern (pump, H+, water and transports closed together) but no source says so.

**Project history (2026-09):** R795 alone and R795 + R1363 were opened in memory (2026-09-10/11) and gave no
flux and no knockout change; four dipeptide supply/chemistry/energy rounds (2026-09-14 to 09-28) needed an
artificial dipeptide source plus R1363, R795, R871 and R876 opened; none was accepted. The default-off
builder step `platform/scripts/gem_annotate/vacuole_candidates.py` (R1363/R795 ub 0.04, R871/R876 ub 0.01)
is that diagnostic scenario, not a curation decision.

## 3. The 37 UM2 reactions (as in `9d655c2b`)

Read from `9d655c2b` this session (`stage1/inventory.txt`). `va` = vacuole `C_va`, `cy` = cytosol,
`ex` = extracellular. Every UM2 reaction is blocked on `9d655c2b` under both the open condition and Task A
(FVA at 0 % of optimum, |v| ≤ 1e-9; `stage1_um2.log`, variant `as_built`), as in the census on `a1237821`.
Gene names and evidence are in section 5. "Closure reason" is what is recorded; section 2 explains why
that is "inherited, reason unknown" for all 17.

**V-ATPase, H+ and water (3 closed)**

| Reaction | Equation | Bounds | GPR | Subsystem | Closure reason |
| --- | --- | --- | --- | --- | --- |
| R795 V-ATPase, vacuole | 2 H+[cy] + ATP + H2O → 2 H+[va] + ADP + Pi | (0, 0) | `(YALI1A11258g and YALI1E14125g and YALI1D00581g and YALI1E10492g and YALI1A09766g and YALI1F20965g and YALI0E16192g and YALI1F31854g and YALI1F21690g and YALI1E37063g and YALI1E32332g and YALI1E12482g) or (YALI1F38820g and YALI1F13017g)` | ATPase | Inherited (iYali4 2016); H+[cy] coefficient −1 → −2 by the 2026-09-10 metadata selection (stoichiometry only, bounds untouched) |
| R1163 H+ diffusion | H+[cy] → H+[va] | (0, 0) | none | Transport, vacuole | Inherited. Written in the uphill direction for an acidic vacuole (section 5) |
| R1363 water diffusion | H2O[cy] → H2O[va] | (0, 0) | none | Transport, vacuole | Inherited |

**Vacuolar amino-acid and polyamine transports (14 closed, 1 open)**

| Reaction | Equation | Bounds | GPR | Closure reason / note |
| --- | --- | --- | --- | --- |
| R869 L-asparagine transport | H+[va] + Asn[cy] → Asn[va] + H+[cy] (import, H+ antiport) | (0, 0) | `YALI1F37658g` | Inherited |
| R870 "L-aspartate transport" | H+[va] + Asn[va] → H+[cy] + Asn[cy] (export, H+ symport) | (0, 0) | `YALI1E13514g` | Inherited; **carries asparagine** despite the name (same in Yeast-GEM `r_1188`) |
| R871 L-aspartate transport | H+[va] + Asp[va] → H+[cy] + Asp[cy] | (0, 0) | `YALI1D14670g` | Inherited |
| R876 L-glutamate transport | H+[va] + Glu[va] → H+[cy] + Glu[cy] | (0, 0) | `YALI1D14670g` | Inherited |
| R878 L-glutamine transport | import, H+ antiport | (0, 0) | `YALI1F37658g` | Inherited |
| R879 L-glutamine transport | export, H+ symport | (0, 0) | `YALI1E13514g` | Inherited |
| R882 L-isoleucine transport | import, H+ antiport | **(0, 1000)** | `YALI1F37658g` | Open since 2016 (the one open `C_va` reaction); blocked because vacuolar Ile has no consumer |
| R883 L-isoleucine transport | export, H+ symport | (0, 0) | `YALI1C04798g or YALI1E13514g` | Inherited |
| R887 L-leucine transport | import, H+ antiport | (0, 0) | `YALI1F37658g` | Inherited |
| R888 L-leucine transport | export, H+ symport | (0, 0) | `YALI1E13514g or YALI1C04798g` | Inherited |
| R897 L-tyrosine transport | import, H+ antiport | (0, 0) | `YALI1E21478g or YALI1F37658g` | Inherited |
| R898 L-tyrosine transport | export, H+ symport | (0, 0) | `YALI1E13514g or YALI1C04798g` | Inherited |
| R918 putrescine transport | H+[va] + Put[cy] → H+[cy] + Put[va] | (0, 0) | `YALI1A15624g or YALI1F10342g` | Inherited |
| R925 spermidine transport | H+[va] + 2 H+[cy] + Spd[cy] → Spd[va] | (0, 0) | `YALI1A15624g or YALI1F10342g` | Inherited; **unbalanced, −3 H and −3 charge** (stoichiometry from the 2026-09-10 metadata selection; Yeast-GEM `r_1261` is 1:1 and balanced) |
| R927 spermine transport | H+[va] + Spm[cy] → H+[cy] + Spm[va] | (0, 0) | `YALI1B27564g or YALI1A15624g or YALI1F10342g` | Inherited |

All 15 are in the "Transport, vacuole" subsystem and are mass- and charge-balanced except R925.

**Dipeptide sets (19 reactions, all open, gene-less, added by iYli21 2022)**

| Set | Exchange | Plasma membrane | Tonoplast | Vacuolar hydrolysis | Amino-acid exit from the vacuole |
| --- | --- | --- | --- | --- | --- |
| Gly-Asp | R2018 Gly-Asp[ex] → (0, 1000) | R2019 H+ symport (±) | R2020 H+ symport (±) | R2021 H2O + Gly-Asp → Asp + Gly (0, 1000) | Asp: R871 (closed); Gly: R2030 |
| Gly-Glu | R2026 (0, 1000) | R2027 H+ symport (±) | R2028 H+ symport (±) | R2029 (0, 1000) | Glu: R876 (closed); Gly: R2030 |
| Ala-Gly | R2031 (0, 1000) | R2032 uniport (±; named "proton antiport") | R2033 uniport (±; same name) | R2034 (0, 1000) | Ala: R2035 uniport (±); Gly: R2030 |
| Gly-Pro | R2036 (0, 1000) | R2037 H+ symport (±) | R2038 H+ symport (±; subsystem "Amino acid metabolism") | R2039 (0, 1000) | Pro: R2040 H+ symport (±); Gly: R2030 H+ symport (±) |

- The four exchanges are written `dipeptide[ex] →` with bounds (0, 1000): **secretion only**; no medium the
  model can be given lets a dipeptide in unless an exchange bound is changed.
- The four hydrolases are named "cytosol nonspecific dipeptidase" but sit in the vacuole, have no gene, and
  cannot be balance-checked (the dipeptides and vacuolar glycine `m1863` have no formula). The 2026-09-24
  audit recorded the same name/compartment conflict as present in iyli21.
- Causes in the census: 20 reactions "several roots", 17 "closed in model (0, 0)".

## 4. Producers, consumers and transporters per compartment

Matched by compound and compartment from the stoichiometry (`trace.txt`, `inventory.txt`; 61 reactions touch a
`C_va` or `C_vm` species).

- **Vacuolar H+ `m1007[C_va]`, 22 reactions, all transports.** Producer: only the pump R795 (closed) and
  the diffusion R1163 (closed). Everything else moves H+ out of the vacuole: the 15 amino-acid/polyamine
  transports (each consumes 1 H+[va]) and the 5 reversible dipeptide/Gly/Pro symports (R2020, R2028, R2030,
  R2038, R2040). Open producers 5 / open consumers 6, all of them the reversible symports and R882.
  **No reaction uses the vacuolar H+ for anything but transport**, and there is no vacuolar ATP synthase.
- **Vacuolar water `m1384[C_va]`, 5 reactions:** consumed by the four dipeptide hydrolyses, produced only
  by R1363 (closed). The vacuolar-membrane water `m1617[C_vm]` is a separate pool (lipid phosphatases),
  supplied by the open R1653.
- **Vacuolar Asn, Gln, Leu, Tyr:** two reactions each, one import and one export of the same compound, both
  closed. **Ile:** import R882 open, export R883 closed. **Asp, Glu:** produced by R2021/R2029, exported by
  R871/R876 (closed). **Putrescine, spermidine, spermine:** one reaction each (import only); no export, no
  consumer in any option. **Gly, Ala, Pro (vacuolar):** produced by the hydrolyses, exported by the open
  reversible R2030, R2035, R2040.
- **Cytosolic counterparts** (`m319`, `m130`, `m743`, `m775`, `m992`, `m267`, `m50`, `m850`, `m467`,
  `m902`, `m272`, `m114`, `m765`): each has open producers, consumers and a plasma-membrane transport;
  none depends on the vacuole.
- **Dipeptides:** each exists as ex/cy/va species only on its own chain; the cytosolic dipeptide has no
  cytosolic hydrolase, so with secretion-only exchanges the chain cannot carry flux on any medium.

**Consequence (proven from the stoichiometry, confirmed by the LPs in section 6):** at steady state every
vacuolar amino acid must leave the vacuole as fast as it enters; with no storage, pH or membrane-potential
term, an import/export pair of the same amino acid nets only "H+ out of the vacuole". The only H+ supply is
the pump. So on any medium where no dipeptide enters, every flux the vacuole can carry is
ATP → R795 → H+[va] → (pair or leak) → H+[cy], an ATP-dissipating futile cycle. The polyamine transports
cannot carry flux at all (one-reaction dead ends).

## 5. Biology

Full notes, URLs and BLAST tables: `stage1/evidence.md`, `stage1/blast/`; the BLAST databases, the PubMed
abstracts and the iYli21 files are in the artifacts folder `stage1/` (gathered on 2026-10-28 by an evidence
subagent; I re-fetched the Russnak 2001 and Preston 1989 abstracts myself and they say what is cited). Tiers: **verified** = the abstract or UniProt/full text was
fetched this session and states the claim; **inferred** = reasoning; none of the W29 genes has a verified
native name, so every function below is **homology only**. A GPR assignment is not experimental validation.
Orthology is reciprocal best BLAST hit (BLAST+ 2.17.0, W29 proteome UP000182444 `a0bacab5…` vs S288C
reviewed).

| Gene | Verified name | Best S. cerevisiae match (RBH?) | Yeast function (evidence) | On reaction(s) |
| --- | --- | --- | --- | --- |
| `YALI1F37658g` | none verified | AVT1 (RBH, 45 % id, 95 % cov) | vacuolar import of Tyr, Gln, Asn, Ile, Leu; H+/aa antiport ([Russnak 2001](https://pubmed.ncbi.nlm.nih.gov/11274162/), [Tone 2015](https://pubmed.ncbi.nlm.nih.gov/25747199/), verified abstracts) | R869, R878, R882, R887, R897 |
| `YALI1E13514g` | none verified | AVT3 (RBH; also AVT4's best hit) | vacuolar efflux of large neutral aa (Russnak 2001, verified) | R870, R879, R883, R888, R898 |
| `YALI1C04798g` | none verified | ATG22 (RBH) | vacuolar efflux of Leu and other aa after autophagy ([Yang 2006](https://pubmed.ncbi.nlm.nih.gov/17021250/), verified abstract); coupling unknown | R883, R888, R898 |
| `YALI1D14670g` | none verified | AVT6 (RBH; also AVT5/AVT7's best hit) | vacuolar efflux of Asp and Glu (Russnak 2001; [Chahomchuen 2009](https://pubmed.ncbi.nlm.nih.gov/20118605/), verified) | R871, R876 |
| `YALI1F10342g` | none verified | TPO1 (RBH) | polyamine transporter; vacuolar in [Tomitori 1999](https://pubmed.ncbi.nlm.nih.gov/9920864/) but mainly plasma-membrane exporter in later work and current UniProt (verified) | R918, R925, R927 |
| `YALI1A15624g` | none verified | TPO1-like paralog (not RBH) | as above, weaker | R918, R925, R927 |
| `YALI1B27564g` | none verified | TPO3 (RBH; also TPO2's best hit) | polyamine transport (Tomitori 2001, verified) | R927 |
| `YALI1E21478g` | none verified | VBA1/2-family paralog (not RBH) | basic-aa vacuolar import; Tyr only weakly (verified) | R897 (also R2055, R2056) |
| V-ATPase subunits (14 genes on R795) | none verified | reviewed in 2026-09 (`r794_r795_gpr_review_20260909`, `vatpase_*`) | — | R795 |

- **Amino-acid import/export (AVT1 / AVT3-4 / AVT6 / ATG22): real in yeast, with orthologs in W29.** The
  model's five import substrates are exactly AVT1's; the export directions match AVT3/4/6. AVT1 and AVT6
  transport needs ATP and is abolished by nigericin (Russnak 2001, verified): the same pH gradient drives
  both directions. That the H+ leaves with the amino acid on export (symport) is **inferred** from that,
  not stated.
- **Polyamine import (R918/R925/R927): weak.** TPO1's location is disputed and putrescine import into the
  vacuole has no direct support in what was read. R925 is also unbalanced.
- **V-ATPase: real and required in Y. lipolytica.** Isolated *Y. lipolytica* vacuoles show H+-ATPase
  activity and build a pH gradient and membrane potential ([Kulakovskaya 1993](https://pubmed.ncbi.nlm.nih.gov/8385405/),
  verified abstract). Core subunits are essential in the PO1f CRISPR screens (Ramesh 2023 table in the
  research root), which is why 11 V-ATPase genes are false negatives. Which W29 a-subunit serves the vacuole
  is unresolved (`YALI1E12482g`, essential, is STV1's best hit; the VPH1 match `YALI1F38820g` is not
  essential). No *Y. lipolytica* vma-mutant or alkaline-pH study was found.
- **H+/ATP = 2 in R795:** measured at a large pH difference ([Kettner 2003](https://pubmed.ncbi.nlm.nih.gov/14645064/),
  verified abstract: about 2, rising above 4 at small ΔpH); a 10-proteolipid c-ring implies about 3.3
  (inferred). So 2 is defensible as the conservative, measured end of the range; it is a modelling choice.
  The source file `model/source/iyali26.xml` has 1 H+[cy] in (H unbalanced); the candidate `9d655c2b` has
  2 H+[cy] (set by the 2026-09-10 metadata selection) and is balanced (checked this session).
- **R1163 as written is the wrong direction.** The yeast vacuole is acidic, pH 6.2
  ([Preston 1989](https://pubmed.ncbi.nlm.nih.gov/2674942/), verified by me), and alkalinizes when the pump is
  lost ([Tarsio 2011](https://pubmed.ncbi.nlm.nih.gov/21669878/), verified abstract), so the passive leak runs vacuole → cytosol. Opened cy → va it is a free H+
  source into the vacuole (and loops, section 6).
- **R1363 water:** no yeast vacuolar aquaporin found; bilayer water permeability is textbook
  (**established, uncited**). Opening it is chemically harmless.
- **Dipeptide sets: an iYli21 addition, not inherited, and not a working route.** iYli21 Supplementary
  Table 1 lists R2018–R2039 as "Manual curation" (verified, file in `sources/`); iYli21 Fig. 1 shows W29
  using all four dipeptides as **carbon** sources on Biolog PM1 (verified), yet the exchanges allow only
  secretion. No *Y. lipolytica* nitrogen-source data were found. In yeast, dipeptides enter by PTR2
  (H+-coupled, plasma membrane; [Ito 2013](https://pubmed.ncbi.nlm.nih.gov/24060756/), verified) and no
  vacuolar dipeptide importer was found (limited search). W29 `YALI1D24836g` is the PTR2 RBH. The route the
  evidence supports would be plasma-membrane uptake plus cytosolic hydrolysis: a separate curation question.
- **Biology vs mechanism.** The biology says the vacuolar V-ATPase and the AVT transporters exist and
  function. The mechanism (section 4, section 6) says the steady-state model has nothing for them to do:
  acidification, storage and pH homeostasis are not represented, so opening them cannot make V-ATPase genes
  essential or change any growth. The V-ATPase false negatives are a missing requirement, not a closure.

## 6. Options and in-memory tests

All tests in memory on `9d655c2b`; no model file written. Scripts `scripts/stage1_um2.py` (13 variants,
5,811 LPs, all optimal) and `scripts/ko_compare_variants.py` (full single-gene sweep, as built vs each
option; Gurobi 1e-9, model reset before every knockout, targeted check of the knocked reactions' flux,
fresh re-solve of flagged points). Building blocks: `pump` R795 (0, 1000); `water` R1363 (0, 1000);
`hdiff` R1163 (0, 1000) as written (cy → va); `hleak` R1163 (−1000, 0) (va → cy only, the passive
direction, as R1159 was set for the Golgi on 2026-09-23); `aaimp` R869, R878, R887, R897; `aaexp` R870,
R879, R883, R888, R898; `aspglu` R871, R876; `pa` R918, R925, R927.

### 6.1 Reactions freed, growth, energy and loops

| Variant | WT growth, Task A | UM2 blocked, open condition | UM2 blocked, Task A (freed) | EGC: max R_NTP1 / biomass, no uptake | No-uptake loop |
| --- | --- | --- | --- | --- | --- |
| as built (A) | 0.311660821633612 | 37 / 37 | 37 / 37 | 0 / 0 (optimal) | none |
| pump | same (Δ < 1e-14) | 37 | 37 | 0 / 0 | none |
| water | same | 37 | 37 | 0 / 0 | none |
| pump + water (B) | same | 24 (frees R795, R1363 and the Ala-Gly and Gly-Pro sets, R2030) | 37 | 0 / 0 | none |
| pump + water + hleak (B′) | same | 23 | 35 (R795, R1163) | 0 / 0 | none |
| pump + water + hdiff | same | 23 | 37 | 0 / 0 | none |
| pump + water + aaexp | same | 22 | 34 (R795, R882, R883) | 0 / 0 | none |
| pump + water + aaimp + aaexp | same | 14 | 26 (R795 + 10 aa transports) | 0 / 0 | none |
| pump + water + aspglu | same | 14 (Gly-Asp, Gly-Glu sets freed) | 37 | 0 / 0 | none |
| pump + water + pa | same | 24 | 37 | 0 / 0 | none |
| pump + water + aaimp + aaexp + aspglu + pa (C) | same | 4 (R1163, R918, R925, R927) | 26 | 0 / 0 | none |
| C + hleak | same | 3 (polyamines) | 25 | 0 / 0 | none |
| **C′** = pump + water + hleak + aaimp + aaexp + aspglu (no polyamines; `stage1_um2_cprime.log`) | same | 3 (R918, R925, R927) | 25 (R795, R1163, R869, R870, R878, R879, R882, R883, R887, R888, R897, R898; R871/R876 need a dipeptide source) | 0 / 0 | none |
| C + hdiff (= the census's "open the 17 closures") | same | 3 | 25 | 0 / 0 | **yes:** R1163 0–1000 and all ten aa transports 0–500 with every uptake closed |

- **Growth:** no variant changes WT growth on Task A (all within 2e-15 of 0.311660821633612) or growth in
  the open condition (122.5654766). pFBA flux of every UM2 reaction is 0 in every variant.
- **What "freed" means on Task A:** R795 has range [0, 55.5] at 0 % of optimum and [0, 0.0565] at 99.9 %
  whenever it is freed, i.e. the only use is spending ATP (pump + leak, or pump + an import/export pair).
  These are ATP-dissipating futile cycles, not new function; the model already has free ATP hydrolysis
  (R_NTP1), so they add no capability.
- **Open condition:** the extra freed reactions are the dipeptide chains, which need a dipeptide uptake that
  no real medium of this model allows (secretion-only exchanges).
- **Energy-generating cycles:** none in any variant (max R_NTP1 and max biomass 0, status optimal, with
  every exchange (0, 1000) and every positive lower bound, incl. `xMAINTENANCE` 7.8625, relaxed; run on both
  the Task A and the open model copies).
- **Thermodynamically infeasible loop:** R1163 opened in its written direction plus any H+-coupled
  import/export pair is an unbounded no-uptake loop (H+ in by "diffusion", out with the amino acid shuttle).
  Not energy-generating, but it decouples every vacuolar H+ coupling and is the same loop class Yeast-GEM
  shows (section 2). The leak direction (va → cy) does not loop.
- **Polyamines:** blocked under every option; opening R918/R925/R927 frees nothing (no vacuolar polyamine
  export or use), and R925 is unbalanced.

### 6.2 Knockout comparison (full sweep, Task A + PO1f, all four cutoffs)

`ko_um2_summary.json`, `ko_um2.tsv`, `ko_um2.log`: 6 sweeps × 1,081 genes + WT = 6,492 LPs, 472 s, 0 points
failing the targeted check, 0 fresh re-solves. The same 20 knockouts are non-optimal (infeasible) in every
variant; 4 of them are mtDNA genes (out of scope), so 16 are unresolved in the scores. Scores use the
user-defined negative class (model genes not on the consensus list; mtDNA genes excluded); the list was used
in development, so they are a development reference, not independent validation.

| Cutoff | TP | FN | FP | TN | unresolved | recall | MCC |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 % | 94 | 221 | 35 | 699 | 16 | 0.298 | 0.350 |
| 5 % | 100 | 215 | 37 | 697 | 16 | 0.318 | 0.363 |
| 10 % | 104 | 211 | 42 | 692 | 16 | 0.330 | 0.361 |
| 15 % (primary) | 112 | 203 | 42 | 692 | 16 | 0.356 | 0.386 |

These are the as-built scores; they match the 2026-10-28 dolichol report for `9d655c2b` at 15 %. **Every
option tested (B, B′, C, C′, and the census's "open all 17" with R1163 as written) gives the identical
table:** 0 knockouts with growth changed by more than 1e-6, 0 knockouts that gain growth, 0 class changes at
any cutoff, WT within 3e-15. Because no knockout gains growth, no targeted bypass controls were needed; the
V-ATPase knockouts stay non-essential under every option (R795 carries no flux that growth needs).

### 6.3 The options

**A — keep every closure; document UM2 as a deliberate closure (recommended).**
- Change: none to bounds. Two forms:
  - **A0:** the record lives in `model/reports/um2_vacuolar_closures_<date>/REPORT.md`, a `STATE.md` entry
    and a note in the census report's UM2 row. No build.
  - **A1:** additionally, a notes-only builder step (curated data `model/curation/vacuole/`, guarded
    apply step, test) writes a `closure_reason` note on the 17 closed reactions and a `dead_inheritance`
    note on the 19 dipeptide reactions. One new candidate whose whole-model diff is notes only and poses
    the same optimization problem.
- For: the history (inherited, reason unrecorded), the mechanism (nothing to do at steady state; R1163 as
  written loops), and no predictive change from any opening.
- Against: leaves 37 reactions blocked whose genes have real orthologs and functions. The closure encodes
  "unused in this model", not "absent in the organism".

**B — open R795 and R1363 (B), optionally with R1163 reversed to a vacuole → cytosol leak (B′).**
- Bounds: R795 (0, 1000); R1363 (0, 1000). For B′, R1163 (−1000, 0) instead of (0, 0), equation kept.
  Opening R1163 in its written direction is not offered (wrong direction, loops).
- Evidence for: *Y. lipolytica* vacuolar H+-ATPase activity (Kulakovskaya 1993); essential V-ATPase core
  genes in the PO1f screens; the passive leak direction (Preston 1989, Tarsio 2011); bilayer water
  permeability.
- Effect:
  - B frees nothing on Task A and 13 reactions in the open condition: R795, R1363 and the Ala-Gly /
    Gly-Pro chains, which need a dipeptide uptake no real medium allows.
  - B′ frees R795 + R1163 on Task A (pump + leak futile cycle only) and 14 in the open condition.
  - No growth, knockout or EGC change, and no loop.
- Against: adds only an ATP-dissipation route the model already has (R_NTP1). The R1163 direction change is
  backed by physiology but changes nothing measurable.

**C — B′ plus the transports with an ortholog and a mechanism (C′, as tested).**
- Opens:
  - R869, R878, R887, R897: AVT1-like import;
  - R870, R879, R883, R888, R898: AVT3/ATG22-like export;
  - R871, R876: AVT6-like Asp/Glu export;
  - R1163 as a leak.
- Excluded: the polyamines R918, R925, R927, which cannot carry flux in any option, have weak transporter
  evidence, and include an unbalanced R925.
- Effect: Task A blocked 37 → 25. The 12 freed reactions carry flux only in pump-driven import/export
  futile cycles; pFBA 0. Open condition 37 → 3.
- No growth, knockout or EGC change, and no loop (the loop needs R1163 in its written direction).
- Against: the same as B, multiplied. FVA ranges become non-zero, which can read as "active" when the only
  mode is futile. The H+ coupling of the export reactions is inferred, not stated in a source.

**Not offered:** "open all 17" (the census test), because R1163 as written makes an unbounded no-uptake
loop. Also not offered: any change to the dipeptide sets. Making them usable means changing exchange bounds,
which changes the medium; moving the hydrolysis to the cytosol and assigning genes (PTR2 RBH
`YALI1D24836g`, the 2026-09 dipeptidase candidates) is a separate curation task. Deleting them would need
evidence that they do not exist, and W29 grows on all four on Biolog PM1.

## 7. Acceptance criteria (for whichever option is chosen)

Any criterion that fails stops the run with a `STOP.md` and a question to you; nothing is tuned, relaxed or
retried.

- **All options:** no other reaction, species, GPR, biomass term or medium changes. R925, the polyamines and
  the dipeptide sets are untouched. Report baseline SHA `9d655c2b…`, output SHA, changed entities and the
  exact command. Add a `STATE.md` entry and an Escher before/after figure (Task A + PO1f), even for A, so
  the map records that the vacuole is closed by decision.
- **A0:** no build; the documentation names every source in section 2 with its evidence tier.
- **A1:**
  - The whole-model diff vs `9d655c2b` is notes only, on exactly the 36 named reactions.
  - Same optimization problem (S, bounds, objective identical).
  - WT and the knockout table are identical to section 6.2.
- **B / B′ / C′:**
  - The whole-model diff is exactly the named bounds (plus notes).
  - WT growth on Task A within 1e-9 of 0.311660821633612.
  - Knockout sweep: 0 class changes at 1, 5, 10 and 15 %; 0 knockouts gaining more than 1e-6; scores
    identical to section 6.2.
  - EGC: max R_NTP1 and max biomass 0, status optimal, with `xMAINTENANCE` relaxed.
  - No-uptake FVA of every UM2 reaction is 0 (no loop).
  - The freed sets equal section 6.1. B: Task A none, open 13. B′: Task A {R795, R1163}, open 14.
    C′: Task A 12, open 34.
  - Each opened reaction carries a note citing the evidence in section 5 and stating that its only
    steady-state flux mode is a futile cycle.

## 8. Question to the user (answered)

Which option should Stage 2 implement: A0, A1, B, B′ or C′ (or none)? **Answer: A0** (section 0).

Follow-ups outside UM2, not started:
1. R925 spermidine transport is unbalanced (−3 H / −3 charge) from the 2026-09-10 metadata selection.
2. Dipeptide route re-curation (PTR2 uptake + cytosolic hydrolysis, exchange direction, gene assignment).
3. R2055/R2056 vacuolar His/Lys import: the VBA1 RBH `YALI1B04391g` (essential in the Cas12a screen) is
   missing from the GPR, and the reactions are uniport while the yeast mechanism is H+ antiport.
4. V-ATPase false negatives: the model has no acidification requirement; representing one would need
   evidence-based constraints, not a recall target.

## 9. Stage 1 computation record

- `inventory.py` → `inventory.txt` (61 vacuole-touching reactions); `trace.py` → `trace.txt`.
- `ext_models.py` → `ext_models.txt`: Yeast-GEM, iYali v4.1.2, iYali5 bounds.
- `yeastgem_loops.py` → `yeastgem_loops.log`.
- `stage1_um2.py` → `stage1_um2.json`, `stage1_um2.log`: 13 variants, 5,811 LPs, all optimal.
- `ko_compare_variants.py` → `ko_um2*`: 6,492 LPs.
- Evidence subagent → `evidence.md`, `blast/`, `sources/`.
- Python `/opt/miniconda3/bin/python`, cobra with Gurobi 13 (academic licence, real HOME).
- No repository file was written. The model was read from a byte-identical scratchpad copy `9d655c2b…`.
