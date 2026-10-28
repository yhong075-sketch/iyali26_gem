# TASK: UM7, the nuclear NAD / homocitrate copies and their transports (2026-10-28)

## 0. Authorization and scope

- **Authorization:** the user, in this session (2026-10-28), after reading the Stage 1 draft:
  **"Approve Option A"**. That is Option A of section 7, implemented as section 8 describes:
  - two gene-less, reversible nucleus ⇌ cytosol diffusion reactions, for ADP-ribose (`m1725[C_cy]` ⇌
    `m818[C_nu]`) and diphosphate (`m203[C_cy]` ⇌ `m210[C_nu]`);
  - subsystem / group "Transport, nucleus";
  - through a guarded `--um7-nuclear-diffusion` step on `9d655c2b`.

  No commit or push is authorized.
- **Not in scope:** R396 / homocitrate, R92, the NAD link (Option A′), R567's chemistry, any GPR.
- **Reaction ids:** `R_ADPRtn` and `R_PPItn` (checked free in 9d655c2b), written in the cytosol ⇌ nucleus form
  of the existing nuclear diffusions (for example R1082), with SBO:0000185 like them.
- **Budget:** 3 builds, about 20,000 LPs. Stop conditions in section 8.
- **Outcome:** the run stopped at the knockout gate and the user restated it (section 11); text corrected after
  the audit is listed in section 12.
- Below, sections 1–10 are the Stage 1 draft as the user read it. Stage 1 outputs are in `stage1/` (BLAST
  tables, `evidence.md`, logs, maps) and `scripts/`. The two Stage 1 JSON result files (100 kB and 72 kB) are
  in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/um7_nuclear_transports_20261028/stage1/`.

**Recommendation: Option A.** Add two gene-less, reversible nucleus ⇌ cytosol diffusion reactions, for
ADP-ribose (`m818[C_nu]` ⇌ `m1725[C_cy]`) and diphosphate (`m210[C_nu]` ⇌ `m203[C_cy]`), of the kind the
model already has 33 of ("Transport, nucleus"). Together they join the nuclear sirtuin step R567 and the
nuclear NAD synthesis steps R569 / R579 to the cytosol, where Ysa1 (R1912) hydrolyses ADP-ribose and Ipp1
(R419) hydrolyses PPi. This frees 9 of the 15 blocked UM7 reactions (7 on Task A) plus R1912 (census module
UM171, whose root is the cytosolic copy of the same ADP-ribose).

In memory there is no change to WT growth, to any knockout class or to the energy-generating-cycle check. The
homocitrate branch (R396 plus four transports) and the closed nuclear acetyl-CoA synthetase R92 are left as
they are: closing the homocitrate branch would need an invented mitochondrial carrier, and the only
*Y. lipolytica* location claim for LYS1 is a mitochondrial prediction.

## 1. Starting state (read 2026-10-28)

- **Worktree / HEAD:** `claude/hungry-mendel-f9f7df` at `0093c89` = `origin/main` after `git fetch`
  (0 ahead / 0 behind), clean.
- **Model executed:** `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_20261028.xml`,
  SHA-256 `9d655c2bdbf82721a38f148a9a4ee0a385bb509ea7b591111107113c4aa8335c` (re-hashed). This is the
  newest candidate on `main`: PR #29 (dolichol route) was merged after the brief was written, which named
  the Kennedy candidate `ce9e1112…` as newest. The UM7 reactions, the context reactions (R1912, R1170, R24,
  R578, R419, R90) and all 47 nuclear species are identical in the two files, field by field (equation,
  bounds, GPR, notes, subsystem); `biomass_C` differs only by Dol-P `m503[C_er]` 1e-4 (verified this session).
  **Not** a designated reference (`model/candidates/README.md`).
- **Instruction files read:** `~/.claude/CLAUDE.md`, the worktree `AGENTS.md`, the `curation-change` skill;
  `model/STATE.md` entries 2026-10-26 to 2026-10-28; `model/expected/benchmark.md`.
- **Research root:** `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research` (`IYALI26_RESEARCH_ROOT`).
- **Condition:** Task A medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) +
  PO1f profile `state/strain_profiles/po1f_sd_leu.json` (`35307853…`, uracil 1000, R612 closed);
  experimental list `consensus_essential_genes.csv` (`1e887f5a…`). All three SHAs match `benchmark.md`.
  "Open" condition = the raw model with every exchange at (−1000, 1000), the census's "any medium".
- **Solver:** Gurobi 13.0.2 at FeasibilityTol / OptimalityTol 1e-9 with the real HOME, through
  `model/reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py` (reset before each solve,
  residual / bound / knocked-flux checks). FVA and pFBA from cobra 0.30 on the same solver settings.
- **Essentiality:** KO/WT strictly below 15 % (also 1, 5, 10 %); genes not on the experimental list count as
  non-essential (the user-defined negative class); the list was used during model development, so results
  are a development reference, not independent validation.
- **Local copies:** libSBML cannot open the candidate through its full path (> 255 characters; the loader
  resolves symlinks), so the LPs read byte-identical copies in the session scratchpad (SHA-256 re-checked).

## 2. UM7 as it is in 9d655c2b (verified: `um7_map.txt`, `um7_trace.txt`, `nucleus_inventory.txt`)

All 17 reactions are iYli21 inheritance: equation topology, bounds and GPR are identical in
`model/source/iyli21.xml`, `model/source/iyali26.xml` and the candidate; only H+ terms differ (proton gate).
Compartment copies of one compound carry different base ids (acetate `m202[C_cy]` / `m205[C_nu]`), so copies
were matched by formula plus a shared ChEBI / KEGG / MetaNetX id, or by formula alone where a species has no
annotation (homocitrate, NMN), never by name.

| Reaction | Name (model) | Equation (model ids; names) | Bounds | GPR | Blocked open / Task A |
| --- | --- | --- | --- | --- | --- |
| R92 | acetyl-CoA synthetase nuclear | acetate + ATP + CoA [nu] → acetyl-CoA + AMP + PPi [nu] | **(0, 0)** | `YALI1F08845g` | yes / yes |
| R396 | homocitrate synthase | acetyl-CoA + 2-oxoglutarate + H2O [nu] → CoA + homocitrate `m628` [nu] | (0, 1000) | `YALI1F38776g` | yes / yes |
| R567 | NAD nucleosidase nuclear (EC 3.2.2.5) | H2O + NAD `m817` [nu] → 2 H+ + ADP-ribose `m818` + nicotinamide [nu] | (0, 1000) | `YALI1D26320g or YALI1D02416g or YALI1F15334g` | yes / yes |
| R568 | NAD synthase | Gln + ATP + H2O + deamido-NAD [cy] → NAD + PPi + Glu + AMP [cy] | (0, 1000) | `YALI1A21156g` | **no / no** (freed by the NAD biomass term, 00e8e397) |
| R569 | NAD synthase | the same in [nu] | (0, 1000) | `YALI1A21156g` | yes / yes |
| R579 | nicotinamide-nucleotide adenylyltransferase | ATP + NMN [nu] → PPi `m210` + NAD [nu] | (0, 1000) | `YALI1E30400g or YALI1A20611g` | yes / yes |
| R580 | nicotinate-nucleotide adenylyltransferase | ATP + NaMN [cy] → PPi + deamido-NAD [cy] | (0, 1000) | `YALI1E30400g or YALI1A20611g` | **no / no** (as R568) |
| R965, R1000, R1001, R1017, R1055, R1082, R1212, R1214, R1267, R1270 | "Transport, nucleus": 2-oxoglutarate, acetate, acetyl-CoA, AMP, CoA, deamido-NAD, Glu, Gln, nicotinamide, NMN | X [cy] ⇌ X [nu] | (−1000, 1000) | none | all yes / yes |

So **15 of the 17** are blocked on 9d655c2b under both conditions (FVA at 0 %, |v| ≤ 1e-9); R568 and R580
were freed by the cofactor biomass terms. The census's UM7 root list still holds.

**Genes (model GPR; identity by BLAST this session, see section 3):**

| Gene | Reactions | Verified name | Function (orthology) | Evidence level | On list |
| --- | --- | --- | --- | --- | --- |
| `YALI1F38776g` | R396, R1170 | LYS1 (UniProt Q12726, YALI0F31075g) | homocitrate synthase; RBH LYS20 / LYS21 (84–89 % id) | Swiss-Prot name; orthology | no |
| `YALI1F15334g` | R567 | none verified | sirtuin; RBH SIR2 / HST1 (42 % id) | orthology | no |
| `YALI1D02416g` | R567 | none verified | sirtuin; RBH HST3 (33 % id, 100 % cov) | orthology; KEGG EC 2.3.1.286 | no |
| `YALI1D26320g` | R567 | none verified | weak sirtuin similarity (26–27 % to SIR2 / HST1 / HST2) | KEGG EC 2.3.1.286; weak orthology | no |
| `YALI1A21156g` | R568, R569 | none verified | glutamine-dependent NAD synthetase; RBH QNS1 (72 %) | orthology | **yes** |
| `YALI1E30400g` | R578–R580 | none verified | NMN / NaMN adenylyltransferase; RBH NMA1 (83 %) | orthology | no |
| `YALI1A20611g` | R578–R580 | none verified | KEGG "putative NMN adenylyltransferase"; no hit from NMA1 / NMA2 at e ≤ 1e-5 | KEGG annotation only | no |
| `YALI1F08845g` | R90–R92, R2024, R2268 | none verified | acetyl-CoA synthetase; RBH ACS2 (63 %) | orthology | no |
| `YALI1F33579g` | R1912 (cytosolic ADP-ribose pyrophosphatase) | none verified | RBH YSA1 (57 %) | orthology; KEGG EC 3.6.1.13 | no |
| `YALI1E03387g` | R24, R395 | LYS4 (UniParc gene name) | homoaconitase; RBH LYS4 (70 %) | orthology | **yes** |

Not in the model: `YALI1C31830g` (Q6CB00 / YALI0C23034g), the RBH of HST4 (44 % id) and second hit of
HST3; it is not on the experimental list.

## 3. Biology, per root (sources in `um7_stage1_20261028/evidence.md`)

Evidence tiers: **verified** = the URL was fetched this session and the claim read in it. The evidence agent
fetched every URL below. I re-opened three myself: Q12726 transit ECO:0000255, the Serrano-Bueno 2019
abstract and the Timney 2016 full-text sentence. **Inferred** = reasoning; **predicted** = UniProt / GO
computational annotation. **No experimental *Y. lipolytica* localization exists for any UM7 enzyme**: every
*Y. lipolytica* compartment below is a prediction, and the compartment evidence is *S. cerevisiae* only.

**Nuclear NAD synthesis and consumption are real in yeast (supports keeping R567 / R569 / R579 nuclear).**
- In *S. cerevisiae* the NAD-consuming sirtuins Sir2, Hst1 and Hst4 are nuclear (UniProt ECO:0000269;
  [P06700](https://rest.uniprot.org/uniprotkb/P06700.txt), verified). Hst2 is mostly cytoplasmic by
  Crm1 export ([PMID 17110954](https://pubmed.ncbi.nlm.nih.gov/17110954/), verified).
- *S. cerevisiae* NAD synthesis is nuclear as well as cytosolic:
  - Nma2 and Npt1 are "both concentrated in the nucleus"
    ([Anderson 2002, PMID 11884393](https://pubmed.ncbi.nlm.nih.gov/11884393/), verified).
  - Nma1 is in cytoplasm and nucleus ([Q06178](https://rest.uniprot.org/uniprotkb/Q06178.txt), verified).
  - GFP-Qns1 is in nucleus and cytosol ([Suda 2003, PMID 12898714](https://pubmed.ncbi.nlm.nih.gov/12898714/),
    verified).
- *Y. lipolytica* orthologs (BLAST, this session): F15334g ↔ SIR2 / HST1, D02416g ↔ HST3 (KEGG "SIR2",
  "HST3"). D26320g is KEGG "sirtuin-5", PANTHER "sirtuin-5, mitochondrial-related", and probably not nuclear
  (inferred). E30400g ↔ NMA1 (UniProt Q6C4M1 "Nucleus", ARBA prediction). A21156g ↔ QNS1.
  `YALI1A20611g` is in the Pof1 family (PTHR31285), not a second Nma (inferred from the family).
- **Caveat on R567's chemistry:** sirtuins make 2′/3′-O-acetyl-ADP-ribose (OAADPr) from NAD and an acetylated
  protein, not free ADP-ribose ([Tanner 2000, PMID 11106374](https://pubmed.ncbi.nlm.nih.gov/11106374/),
  verified; Rhea 43636 in the UniProt entries). No source was found for an acetyl-substrate-independent NAD
  glycohydrolase. R567 "NAD + H2O → ADP-ribose + nicotinamide" is a lumped stand-in, which lumps
  deacetylation with OAADPr deacetylation (inferred). Its disposal route is the same either way: Ysa1
  hydrolyses both ADP-ribose and OAADPr to AMP + (acetyl-)ribose-5-P
  ([Rafty 2002, PMID 12370179](https://pubmed.ncbi.nlm.nih.gov/12370179/);
  [Tong 2009, PMID 19251690](https://pubmed.ncbi.nlm.nih.gov/19251690/), both verified). This draft does not
  propose rewriting R567.

**Root ADP-ribose [nu] (RNC), with ADP-ribose [cy] the matching RNP of UM171.** The missing step is
transport, not a new enzyme.
- The only ADP-ribose consumer in yeast is Ysa1 (EC 3.6.1.13). In *S. cerevisiae* it is mitochondrial by IDA
  (Tong 2009) and nuclear and cytoplasmic by high-throughput GFP only
  ([SGD YBR111C](https://www.yeastgenome.org/locus/YBR111C), verified via the SGD GO backend).
- A nuclear ADP-ribose / OAADPr consumer is weakly supported: high-throughput data only, and the nuclear
  OAADPr activity of Rafty 2002 has no assigned enzyme.
- The *Y. lipolytica* ortholog `YALI1F33579g` (Q6C0B7, Nudix ADP-ribose pyrophosphatase; GO cytosol and
  nucleus by IBA) is already on cytosolic R1912.
- **Mechanism: passive diffusion through the nuclear pore.** Timney et al. 2016 (in vivo, *S. cerevisiae*)
  define passive diffusion as "the equilibration of ions, small molecules, and small- to medium-sized
  macromolecules" between cytoplasm and nucleoplasm
  ([PMC5057280](https://pmc.ncbi.nlm.nih.gov/articles/PMC5057280/); full text verified).
  [Wente & Rout 2010, PMID 20630994](https://pubmed.ncbi.nlm.nih.gov/20630994/) says the same for water,
  sugars and ions (abstract verified; full text unread). Neither names ADP-ribose; applying the rule to a
  559 Da metabolite is **inferred**.
- yeast-GEM made the same choice (gene-less r_4603); that is precedent, not evidence.

**Root diphosphate [nu] (RNC).** Transport, or a nuclear pyrophosphatase; transport is the smaller change.
- PPi is released in the nucleus by NMNAT and Qns1 (above) and by nucleocytosolic Acs2 in *S. cerevisiae*
  ([Takahashi 2006, PMID 16857587](https://pubmed.ncbi.nlm.nih.gov/16857587/), verified).
- *S. cerevisiae* Ipp1 is "a nucleocytoplasmic protein" with a functional NLS, and nuclear Ipp1 is kept low
  by proteasomal degradation ([Serrano-Bueno 2019, PMID 30826332](https://pubmed.ncbi.nlm.nih.gov/30826332/),
  abstract read this session).
- The *Y. lipolytica* IPP1 (Q6C1T4, YALI0F13541g; W29 RBH `YALI1F18105g`, 77 % to IPP1) is annotated
  "Cytoplasm" by similarity. It is already on cytosolic R419 and mitochondrial R420.
- The two routes are LP-equivalent here, since cytosolic PPi is hydrolysed by R419. Diffusion needs no new
  gene and no GPR guess. A nuclear R419 copy would add a third reaction to the IPP1 GPR, and *Y. lipolytica*
  has no nuclear Ipp1 data (inferred).

**Root homocitrate [nu] (RNC).** The biology does not resolve it, and no change is recommended.
- *Y. lipolytica* has a single homocitrate synthase, LYS1 (Q12726). Its only location claim is a predicted
  mitochondrial transit peptide (ECO:0000255; [Q12726](https://rest.uniprot.org/uniprotkb/Q12726.txt),
  re-opened this session). The cloning paper only predicts it
  ([Perez-Campo 1996, PMID 8948100](https://pubmed.ncbi.nlm.nih.gov/8948100/), verified).
- In *S. cerevisiae*, Lys20 / Lys21 are experimentally nuclear (immunofluorescence and fractionation;
  [Chen 1997, PMID 9099739](https://pubmed.ncbi.nlm.nih.gov/9099739/);
  [Feller 1999, PMID 10103047](https://pubmed.ncbi.nlm.nih.gov/10103047/); both verified). Lys20 also
  appears in mitochondrial proteomes (high-throughput only).
- Lys20 has a nuclear chromatin role that is independent of catalysis
  ([Scott & Pillus 2010, PMID 20810648](https://pubmed.ncbi.nlm.nih.gov/20810648/), verified).
- The homocitrate consumers Aco2 and Lys4 are mitochondrial by experiment
  ([P39533](https://rest.uniprot.org/uniprotkb/P39533.txt),
  [P49367](https://rest.uniprot.org/uniprotkb/P49367.txt); verified).
- **No mitochondrial homocitrate carrier is known** (PubMed searches; absence of evidence).
- So a nuclear R396 is real in *S. cerevisiae*, unsupported in *Y. lipolytica*, and in either organism the
  model could only be connected through an invented carrier. The model already makes homocitrate in the
  mitochondrion (R1170, same gene, carries the lysine flux 0.0798 at the optimum).
- Leaving R396 blocked costs nothing: no gene, flux or essentiality depends on it (verified, section 5).
- Closing it would only re-label it. Because the *S. cerevisiae* nuclear evidence is experimental, a closure
  would rest on a *Y. lipolytica* prediction alone, so it is not proposed.

**Root acetate [nu] and the closed R92.**
- *S. cerevisiae* Acs2 is nucleocytosolic (Takahashi 2006). *Y. lipolytica* ACS1 (`YALI1F08845g`, Q6C2Q5)
  has only an IBA cytosol annotation.
- R92 is (0, 0) already in `iyli21.xml`, and no reason is recorded.
- Opening it would add a duplicate of cytosolic R90 (same gene) whose only nuclear outlet is the acetyl-CoA
  transport R1001. It frees R92 and R1000 only together with the homocitrate route.
- Not proposed. The closure is a separate question for whoever owns the R91 / R92 closures.

**Collapse the nucleus into the cytosol?** Not supported.
- The *S. cerevisiae* evidence places sirtuins, NMNATs, Qns1, Acs2, Ipp1 and HCS in the nucleus by experiment.
- The nucleus carries 55 reactions across UM7, UM11 and UM20, so collapsing it would be a whole-compartment
  redesign.
- Deleting reactions needs evidence that they do not exist in the organism, and there is none.

**Orthologs for anything added.** The recommended option adds no gene, so no new BLAST is needed. Every
reaction it touches already carries a gene with a W29 reciprocal best hit (section 2; `blast/` in the stage-1
folder: BLAST+ 2.17.0, W29 proteome UP000182444 `a0bacab5…`, S288C `5cb84d0d…`, queries from UniProt
reviewed entries `8fa5a7c5…`). The W29 TrEMBL accessions (A0A1D8…, A0A1H6…) are now obsolete in UniProtKB;
UniParc maps them to the active CLIB122 entries.

## 4. Per-compartment supply of the UM7 metabolites (verified: `um7_trace.txt`)

| Compound | Nucleus | Cytosol | Mitochondrion | Nucleus ↔ cytosol transport |
| --- | --- | --- | --- | --- |
| homocitrate | `m628`: made by R396 only; **no consumer** (RNC) | no species | `m75`: made by R1170 (same gene), used by R24 homoaconitase | none (no cytosolic species) |
| diphosphate | `m210`: made by R569, R579, R92 (closed); **no consumer** (RNC) | `m203`: 79 reactions incl. R419 (PPase), transports to pe / mi / er / lp | `m204`: R420 PPase, R1309 | **none** |
| ADP-ribose | `m818`: made by R567; **no consumer** (RNC) | `m1725`: consumed by R1912 (Ysa1); **no producer** (RNP, census UM171) | none | **none** |
| NAD | `m817`: made by R569, R579; used by R567 | `m122`: 65 reactions; in `biomass_C` | `m27` | **none** |
| nicotinamide | `m819`: made by R567 | `m829` | – | R1267 |
| acetyl-CoA | `m208`: R396 (consumer), R92 (closed) | `m68` | `m71` | R1001 |
| acetate | `m205`: R92 (closed) only | `m202` | `m201` | R1000 |
| NMN / deamido-NAD | `m832` → R579; `m821` → R569 | `m831` (R713, R578, R1904, transports); `m820` | – | R1270; R1082 |
| NADP | no nuclear species | `m41` | `m178` | – |

The nucleus as a whole has 47 species and 55 reactions, 33 of them gene-less nucleus ↔ cytosol transports
("Transport, nucleus", "diffusion"); PPi, NAD and ADP-ribose are among the few nuclear species without one.
yeast-GEM (the local copy, `reference/external_models/yeast_gem/Yeast-GEM.xml`, SHA-256 `9fd2c572…`, repo
commit 80a6e27) has the same nuclear block (r_0767 sirtuin "NAD nucleosidase nuclear" with SIR2 / HST1–4,
r_0769, r_0784, r_0543, r_0114 open) and adds the gene-less **ADP-ribose transport r_4603** (nucleus ↔
cytosol) plus the cytosolic Ysa1 reaction r_4160; it has no nuclear PPi or homocitrate outlet either.

## 5. In-memory LP results (Task A + PO1f, Gurobi 1e-9, checked; open condition FVA)

Test reactions are gene-less, reversible, prefixed `T_`; nothing was written to a model file.
`stage1_um7.py` (2,142 LPs) → `stage1_um7.json`, `stage1_um7.log`.

| Variant | WT (Task A) | UM7 blocked open / Task A | Also freed | EGC no-uptake (ATP, biomass) |
| --- | --- | --- | --- | --- |
| as built | 0.3116608216 | 15 / 15 | – | 0, 0 (optimal) |
| ADP-ribose nu ⇌ cy only | same | 15 / 15 | – | 0, 0 |
| PPi nu ⇌ cy only | same | 15 / 15 | – | 0, 0 |
| NAD nu ⇌ cy only | same | 15 / 15 | – | 0, 0 |
| **ADP-ribose + PPi** | same | **6 / 8** | R1912 (UM171) | 0, 0 |
| ADP-ribose + PPi + NAD | same | 6 / 8 | R1912 | 0, 0 |
| ADP-ribose + PPi + NAD, R396 closed | same | 6 / 8 | R1912 | 0, 0 |
| ADP-ribose + PPi + NAD + homocitrate nu ⇌ cy ⇌ mi | same | 2 / 4 | R1912 | 0, 0 |
| … + R92 opened | same | 0 / 2 | R1912 | 0, 0 |
| census H11 (homocitrate nu ⇌ mi, PPi and ADP-ribose sinks) | same | 2 / 4 | – (R1912 stays blocked) | 0, 0 |
| collapse (R567 rewritten to cytosol, nuclear R567 / R569 / R579 closed) | same | 15 / 15 (closed) | R1912 | 0, 0 |

- WT differs between variants by < 1e-15; every WT solve passed (residual ≤ 5e-14, no bound violation).
  pFBA at the optimum carries no flux through any UM7 reaction or test link in any variant: the changes
  widen the flux space and do not move the optimum.
- **ADP-ribose + PPi is the minimal pair.** Either alone frees nothing, because R567 needs an ADP-ribose
  outlet and R569 / R579 need a PPi outlet to supply R567 with NAD. Freed (open): R567, R569, R579, R1017,
  R1082, R1212, R1214, R1267, R1270, plus R1912 outside UM7. On Task A, R579 and R1270 stay blocked because
  cytosolic NMN has no source on that medium (nicotinamide riboside, NMN uptake and the peroxisomal R564
  route are unavailable; R578 is blocked on Task A as built too).
- **Ranges:** with ADP-ribose + PPi the sirtuin → Ysa1 route can run (Task A FVA at 0 %: R567 / R1912 up to
  10.7; at 99.9 %: up to 0.0116). It is an NAD-consuming, ATP-costing salvage loop (NAD → ADP-ribose +
  nicotinamide; ADP-ribose → AMP + ribose-5-P; nicotinamide re-enters NAD synthesis), not an energy source;
  the EGC LPs are 0.
- **NAD link:** frees nothing more; it lets cytosolic NAD feed the nuclear sirtuin (R568 open range
  0.395 → 849, R578 0.395 → 1000).
- **Homocitrate branch:** R396, R965, R1001 and R1055 are freed only by a homocitrate route out of the
  nucleus into the mitochondrion. With it, nuclear R396 can replace part of mitochondrial R1170 (R1170 99.9 %
  minimum 0.0797 → 0.0443); both share `YALI1F38776g`, so no knockout changes. R92 / R1000 additionally need
  R92 opened.
- **Knockouts (Task A, checked, stale-point test passed, max knocked flux 0):** in every variant the ten
  genes above give the same ratios as built: `YALI1A21156g` 0 (on list, TP), `YALI1E03387g` and
  `YALI1F38776g` 0.0694 (the first on the list, TP; the second not, FP at 15 %), `YALI1F08845g` 0.9986, all
  sirtuin, NMNAT and Ysa1 genes 1.0. **No class change at any cutoff**, so the controls (WT, target KO,
  bypass closed = as built, combined = variant + KO) give the same verdicts.
- **Duplicates / bypasses:** the links make nuclear R569 / R579 run in parallel to cytosolic R568 / R578.
  These are duplicates with the same GPRs, so they cannot rescue a knockout. No net-synthesis bypass of a
  biomass precursor is created; the only new cycle is the dissipative NAD salvage loop above.

## 6. Whole-model blocked sets (FVA at 0 % of every reaction, Gurobi 1e-9, |v| ≤ 1e-9)

`stage1_um7_census.py` → `stage1_um7_census.json`, `stage1_um7_census.log`; about 37,000 FVA LPs.

| Variant | Blocked open | Blocked Task A | Freed open | Freed Task A | Newly blocked |
| --- | --- | --- | --- | --- | --- |
| as built | 705 | 837 | – | – | – |
| **ADP-ribose + PPi (A)** | 695 | 828 | R567, R569, R579, R1017, R1082, R1212, R1214, R1267, R1270, **R1912** | R567, R569, R1017, R1082, R1212, R1214, R1267, R1912, **R576** | **none** |
| + NAD (A′) | 695 | 828 | same as A | same as A | none |
| + homocitrate route (B) | 691 | 824 | A + R396, R965, R1001, R1055 | A + R396, R965, R1001, R1055 | none |

- R576 is the nicotinamidase (`YALI1A22129g`, nicotinamide → nicotinate + NH3). It is blocked on Task A as
  built because nicotinamide has no source there. Under A, the sirtuin step supplies nicotinamide, which then
  returns through the salvage route.
- The open-condition count at Gurobi 1e-9 (705) differs from the GLPK census number in STATE.md (708). That
  is the solver fringe documented for R1343 / R1344; it is not compared across solvers here.

## 7. Options

| Option | Change | Frees (open / Task A, of 15 blocked) | For | Against |
| --- | --- | --- | --- | --- |
| **A (recommended)** | Two gene-less reversible "Transport, nucleus" reactions: ADP-ribose `m818[C_nu]` ⇌ `m1725[C_cy]` and diphosphate `m210[C_nu]` ⇌ `m203[C_cy]` | 9 / 7, plus R1912 (UM171); R576 on Task A | Same mechanism as the 33 existing nuclear diffusions; nuclear NAD turnover and PPi release are experimental in *S. cerevisiae*; joins an RNC / RNP pair; yeast-GEM precedent for ADP-ribose; no new gene, species or GPR; WT, essentiality and EGC unchanged | *Y. lipolytica* localization is predicted only; diffusion of these exact metabolites is inferred from a general rule; R567 stays a lumped stand-in for sirtuin chemistry |
| A′ | A + NAD `m817[C_nu]` ⇌ `m122[C_cy]` | 9 / 7, plus R1912 | NAD is just as diffusible; lets cytosolic NAD reach the nuclear sirtuins | Frees nothing more; a third reaction with no LP-visible gain |
| B | A′ + homocitrate route: new species homocitrate[C_cy], nu ⇌ cy diffusion, cy ⇌ mi carrier | 13 / 11, plus R1912 | Frees R396, R965, R1001, R1055 | No known mitochondrial homocitrate carrier (invented transport); *Y. lipolytica* LYS1 is predicted mitochondrial; nuclear HCS can replace part of mitochondrial R1170 |
| C | Close R396 with a reason | 0 | Records the *Y. lipolytica* prediction | Rests on a prediction against experimental *S. cerevisiae* data; changes no LP outcome |
| D | Open R92 | 0 alone; 2 more with B | Acs2 is nucleocytosolic in *S. cerevisiae* | Closure reason unknown; duplicate of R90 |
| E | Collapse the nuclear copies into the cytosol | (closes the copies) | Fewer compartment copies | Deletes reactions that exist in yeast; whole-compartment scope |
| F | Census sinks on PPi / ADP-ribose [nu] | 13 / 11 | – | Not biology; R1912 stays blocked |
| 0 | No change; record the analysis | 0 | No risk | Leaves an RNC / RNP pair that has a known mechanism |

## 8. Recommended implementation (Option A), following the curation-change skill

1. `model/reports/um7_nuclear_transports_<date>/TASK.md` (this draft, the user's sentence quoted) and
   `before.json` (the two endpoints of each new reaction, R567 / R569 / R579 / R1912 and the 15 blocked ids).
2. Data `model/curation/nuclear_transport/um7_nuclear_diffusion.json`: one row per reaction (new id, name,
   species, bounds (−1000, 1000), subsystem "Transport, nucleus", no GPR, evidence URLs and tiers as in
   section 3, reviewer). Proposed ids follow the existing nuclear-transport naming; they will be checked free
   before the build.
3. Guarded step `platform/scripts/gem_annotate/um7_nuclear_diffusion.py`. It fails closed if:
   - a species is missing, or its formula or charge differs between the two ends (m818 / m1725 C15H21N5O14P2
     −2; m210 / m203 H4O7P2 0);
   - a reaction id already exists;
   - any transport between the pair already exists.

   Flag `--um7-nuclear-diffusion`, requiring `--dolichol-route` (so the chain stays ordered). Focused test
   `platform/tests/test_um7_nuclear_diffusion.py` on an in-memory model. `python -m pytest` before / after,
   compared test by test.
4. Builds: a control build reproducing `9d655c2b` byte for byte, then the new candidate
   `…_dolichol_route_v2_um7_diffusion_<date>.xml` with `PYTHONHASHSEED` 0 and 1 (must be identical). Record
   the command, the output SHA-256 and `.build.json`. Budget **3 builds**.
5. Whole-model diff vs `9d655c2b`: exactly two added reactions (and their notes). No species, gene, bound,
   GPR, formula, charge or annotation change elsewhere. Both reactions mass- and charge-balanced.
6. Solves (Task A + PO1f, Gurobi 1e-9, checked):
   - WT and pFBA;
   - EGC (maintenance relaxed, status reported, ATP and biomass);
   - knockout sweep 2 × 1,082 at 1 / 5 / 10 / 15 % with TP / FN / FP / TN, positive coverage and
     out-of-scope positives;
   - whole-model FVA blocked sets, open and Task A.
7. Independent read-only audit; REPORT.md; `model/STATE.md` entry.
8. Escher figure (skill `gem-route-escher-maps`): before / after pFBA and the FVA ranges of the touched
   reactions on Task A + PO1f: the nuclear NAD block (R567, R569, R579, the new links, R1912, R568, R578, R580
   and the existing nuclear transports). PNG, SVG and a self-contained HTML viewer in the report folder.
   pFBA flux through these reactions is expected to be 0 at the optimum, so the map will also show FVA ranges
   (blocked → free) to make the change visible.

**Acceptance criteria (Option A):**
- The diff is exactly the two reactions.
- WT 0.3116608216 unchanged within 1e-9 relative.
- Blocked UM7 reactions in the open condition 15 → 6 (R567, R569, R579, R1017, R1082, R1212, R1214, R1267,
  R1270 freed), and on Task A 15 → 8 (R579, R1270 stay, because the medium has no NMN); R1912 freed (and R576 on Task A);
  whole-model FVA blocked sets 705 → 695 (open) and 837 → 828 (Task A), with nothing newly blocked.
- EGC: ATP and biomass 0 with status optimal.
- Knockout sweep: no class change at any cutoff, and no knockout gaining absolute growth above 1e-6.
- Control build byte-identical; seeds identical.

**Stop conditions:** any unintended diff entry, any newly blocked reaction, an EGC > 0 or non-optimal status,
any class change, or a build budget overrun → stop and report.

**Compute budget:** 3 builds; about 20,000 LPs (the whole-model FVA is 4 × 2 × 2,318, the sweep
2 × 1,082, plus the WT, EGC and Escher pFBA / FVA solves).

## 9. Not addressed

- R396 / homocitrate and R92 / R1000 / R965 / R1001 / R1055 stay blocked (sections 3 and 7).
- R567's lumped chemistry stays as it is.
- The HST4 ortholog `YALI1C31830g` stays outside the model.
- No *Y. lipolytica* localization data exist for any of these enzymes.
- Only the Task A medium is tested for growth and essentiality.
- UM11 / UM20 (other nuclear modules) and UM2 (the next target) are not touched.

## 10. Approval sentence (to be quoted verbatim in TASK.md if given)

"Approve Option A" (or A′ / another option, named), given in this session.

Given: **"Approve Option A"** (2026-10-28, this session).

## 11. Stop and restated gate (2026-10-28)

The run stopped at the knockout gate ([STOP.md](STOP.md)). The criterion "no knockout gaining absolute growth
above 1e-6" (section 8) failed for one gene. `YALI1F37714g` (purine-nucleoside phosphorylase by KEGG EC
2.4.2.1; name not verified) goes from 0.3114518 to 0.3115725, +1.21e-4 absolute, with no class change at any
cutoff. The gain runs through the new route: kynurenine NAD synthesis, then the sirtuin step R567, then R1912.
Closing either link, R567 or R1912 restores the baseline value exactly.

**User decision, quoted:** "Accept, restate gate". That is the option offered in STOP.md:

- keep `5f54eb22…`;
- record the +1.2e-4 on `YALI1F37714g` as a documented consequence;
- finish the run (audit, REPORT, STATE, README).

**Restated knockout criterion** (replaces the one in section 8 for this run only):
- no class change at any cutoff;
- the only knockout that gains more than 1e-6 in absolute growth is `YALI1F37714g`, by +1.21e-4, through the
  mechanism in STOP.md.

All other acceptance criteria and stop conditions stay as written. No model change follows from this decision,
and constraining R567 stays out of scope.

## 12. Corrections after the audit (2026-10-28, [AUDIT.md](AUDIT.md))

Sections 1–10 are kept as the user read them. These statements in them are corrected here:

- **Paths and names:**
  - Section 3's `um7_stage1_20261028/evidence.md` is now [stage1/evidence.md](stage1/evidence.md).
  - The candidate built in section 8 is `…_dolichol_route_v2_um7_nuclear_diffusion_20261028.xml`, not
    `…_um7_diffusion_<date>.xml`.
- **Section 1 "Local copies":** the long-path failure comes from cobra 0.30, not libSBML. cobra rejects str SBML
  paths of 260+ characters when `"win" in sys.platform`, which `darwin` matches. libSBML reads the file.
- **Section 3, sirtuins:** "Sir2, Hst1 and Hst4 are nuclear" is attached to P06700, which covers Sir2 only.
  - Sir2 and Hst4 are nuclear by ECO:0000269 ([P06700](https://rest.uniprot.org/uniprotkb/P06700.txt),
    [P53688](https://rest.uniprot.org/uniprotkb/P53688.txt)).
  - Hst1 is nuclear only by the large-scale GFP screen ([P53685](https://rest.uniprot.org/uniprotkb/P53685.txt),
    which also notes cytoplasm).
  - The curation file's evidence row carries the same short citation. It is left unchanged, because the build
    records that file's hash.
- **Section 3, Ysa1:** "mitochondrial by IDA (Tong 2009)" rests on SGD's GO annotation
  ([YBR111C](https://www.yeastgenome.org/locus/YBR111C)), not on the Tong 2009 abstract, which states no
  location.
- **Sections 5 and 8, pFBA:** "no flux through any UM7 reaction" applies to the blocked UM7 reactions and the
  links. R568 and R580 (already active) carry 0.0010 at the optimum in both models.
- **Section 8:**
  - The sweep is 2 × 1,081 genes (`PO1f_plasmid_LEU2` is excluded at runtime) + 2 WT = 2,164 LPs, not
    2 × 1,082.
  - The guard accepts an id that already exists with this curation's exact definition (rerun no-op). Only an id
    with another definition fails.
  - The plan's separate "a transport between the pair already exists" guard was dropped as unreachable
    (REPORT.md section 3).
- **Section 8 budget:** about 60,000 LPs were used against about 20,000 planned (REPORT.md section 6).
