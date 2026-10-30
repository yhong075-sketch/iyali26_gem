# V-ATPase false negatives: diagnosis and options (investigation, no model change)

Scope and inputs are in [TASK.md](TASK.md). Model `059316ef…` from `origin/main`, read-only, no LP.

## 1. Answer

The model cannot be made to call the 11 V-ATPase genes essential with any change that has evidence behind it today.
Two causes stack, and fixing either one alone changes no call:

1. **The gene rule.** No single knockout makes it false.
2. **The network.** Nothing the model needs uses an acidified Golgi or vacuole, so R794/R795 carry no needed flux.

**Recommendation:**
- Keep the 11 genes counted as FN, labelled as a structural class: "non-metabolic function, no steady-state demand".
- Optionally correct the rule, for correctness only.
- Gate any coupling of V-ATPase to growth on a direct *Y. lipolytica* experiment (section 6).

## 2. Method

A read-only workflow of eight agents ran on 2026-10-30
([script](scripts/workflow_vatpase_fn_options.js), [full output](evidence/workflow_result.json)):
- four readers, in parallel: report history, model structure, literature, experimental calls;
- one designer, which built the option set;
- three adversarial reviewers, with three lenses: project rules and decisions, model mechanics, and biological
  evidence.

The lead then re-ran the decisive static checks on `059316ef` in one script
([static_checks.py](scripts/static_checks.py), output [static_checks.json](evidence/static_checks.json)).

**Code logic of the lead's check:**
- load the model with cobra 0.30.0 from a `pathlib.Path`;
- evaluate the R794 rule with each of its 14 genes knocked out, using `GPR.eval`, for the current rule and the
  corrected rule;
- list every reaction around Golgi H+, vacuolar H+, the ER and Golgi sphingolipid species, and cytosolic Fe2+;
- list the biomass species in the Golgi or vacuole, and the boundary reactions that carry iron.

No solver call is made. The agents' own scripts are kept in [scripts/agents/](scripts/agents/). Literature texts
and downloads (17 MB) are in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/vatpase_fn_options_20261030/`.

## 3. Diagnosis

### 3.1 Gene rule (verified, lead re-run)

R794 (Golgi) and R795 (vacuole) share one rule:

```
(A11258 and E14125 and D00581 and E10492 and A09766 and F20965 and YALI0E16192g and F31854 and F21690
 and E37063 and E32332 and E12482) or (F38820 and F13017)
```

| Rule | Single knockouts that make it false (of 14) |
| --- | --- |
| Current | 0 |
| Corrected: core AND (E12482 or F38820) | 12, all except the two a-subunits |

The second branch is only a membrane a-subunit plus the c'' ring subunit, with no ATP-hydrolysing head. The
screens contradict it. F38820g is non-essential in both CRISPR screens. F13017g (c'') is the most strongly
essential gene of the set, yet it is missing from the first branch.

### 3.2 Network (static argument; parts verified, parts recalled)

- **Not in biomass (verified):** biomass_C contains no Golgi or vacuolar species.
- **Vacuole:** R795 is [0, 0] by the user's Option A0 decision of 2026-10-28 (verified bounds). Of the 22
  reactions around vacuolar H+, all non-transport exits are closed. The open ones are gene-less reversible
  symports tied to blocked solute chains (verified neighbour list; blocking of the chains recalled from
  2026-09-10 and 2026-10-28).
- **Golgi H+ has five reactions (verified):**
  - R794, the pump;
  - R1159, a leak to the cytosol, [-1000, 0];
  - R141, structurally dead (recalled from `screen_vatpase_candidate_20260921/MECHANISM.md`);
  - R2048;
  - R_NDP7g.
- **R2048 is forced to 0 (verified, lead re-run):**
  - Its ER partner M(IP)2C `m676[C_er]` has no producer. Its only other reaction is R439, which is irreversible.
  - R439's product `m674[C_er]` has 10 irreversible producers and no consumer, so R439 = 0.
  - Then R1240 = 0, and the Golgi pool `m1314` touches only R2048 and the closed R1817.
- **R794 can only waste ATP:** its protons can leave only through R1159, or through the R_NDP7g reverse + R586
  cycle, which destroys H+. The reverse-only direction of R_NDP7g was found by the mechanics reviewer and not
  re-run by the lead. Setting R794 and its sink to zero keeps every row balanced, and R_NTP1 [0, 1000] absorbs
  the freed ATP. **So closing R794 cannot lower maximum growth.** This explains the 2026-09-21 screen of the
  corrected rule, which changed 0 calls (recalled).

## 4. Options

Every "flips" figure is a Boolean or structural inference; no LP was solved.

| Option | What changes | 11 FN | New FP | Verdict |
| --- | --- | --- | --- | --- |
| O1 rule correction | core AND (a1 or a2); drop the c''-only branch | 0 flip (static argument, 3.2) | 0 | Correct for its own sake; blocked by the F-subunit ID and a-isoform compartment questions |
| O2a iron coupling | add the V-ATPase core genes as an AND on iron uptake R2193 | 10 flip (11 if a1 alone) | D00581, YALI0E16192g | **Biology refuted**; do not build |
| O2b vacuolar cysteine storage | new vacuolar pool, antiporter, reopened R795 | not representable in FBA | — | Park; reopens A0 |
| O3 generic Golgi H+ demand | fixed Golgi H+ demand plus O1 | 10 flip, because R794 becomes the only Golgi H+ source | D00581, YALI0E16192g | No measured flux, so it would be tuning to the match rate |
| O4 structural FN class | report label only | 0 | 0 | **Adopt** |
| O6 *Y. lipolytica* experiment | none in the model | — | — | **Gate** for any coupling |

**O2a mechanics hold, its biology does not.**

Mechanics:
- R2193 is the only producer of cytosolic Fe2+ `m948` (verified).
- R1189 is the only boundary reaction carrying an iron species (verified by name; the iron species have no formula).
- Biomass needs `m948` at 3.04e-5 per unit growth (verified).

The biology fails three ways:
1. Hughes 2020 found that total cell iron does not fall, and even rises, in V-ATPase mutants. Respiratory failure
   comes from cytosolic cysteine lowering mitochondrial iron availability, not from impaired uptake.
2. YALI1D09277g on R2193 is annotated FET5-like in the consensus list and FET3 by reciprocal best hit in the R2193
   report. *S. cerevisiae* Fet5 works at the vacuole.
3. The genome annotates iron-route homologs absent from the model: SMF1-like, SMF3-like, FTR1-like and ARN1-like.
   All of those checked are non-essential in the screens.

Encoding a maturation dependency as an AND also departs from the project's catalysis/complex GPR convention
(AGENTS.md §5). And it reopens the 2026-10-22 R2193 GPR decision. If it were ever built, the diagnostic must use
Gurobi at 1e-9. At 1e-7 a feasibility leak on the iron row allows a KO/WT of about 1 %, which would fail the 1 %
cutoff (arithmetic by the mechanics reviewer).

**O3** is mechanically effective on the Golgi side, because R2048 is blocked. That is exactly why it needs a
measured Golgi H+ flux before it can count as anything other than tuning. A vacuolar version reopens A0.

**O4** keeps every count unchanged. The genes have model counterparts, so they stay FN under
`model/expected/benchmark.md`. The mtDNA exclusion in that file is a separate, screen-side decision and is not a
precedent.

## 5. Genes

No native name is verified for any gene. Functions are homology annotations from the consensus list and the W29
locus map. Raw calls come from `cas9_cas12a_fitness.csv` (FS = fitness score; q = FDR), verified by the data agent.

| Gene | Subunit (homology) | Consensus | Cas9 FS / q | Cas12a FS / q |
| --- | --- | --- | --- | --- |
| YALI1A09766g | V1 A, catalytic (VMA1-like); also on R_NTP1 | E | -3.634 / 6.3e-4 | -1.940 / 0.027 |
| YALI1E32332g | V1 B (VMA2-like) | E | -5.459 / 4.3e-12 | -1.671 / 0.151 |
| YALI1A11258g | V1 C (VMA5-like) | E | -3.782 / 2.0e-4 | -2.344 / 6.9e-4 |
| YALI1D00581g | V1 D (VMA8-like) | not E | -2.627 / 0.171 | -1.784 / 0.079 |
| YALI1F20965g | V1 E (VMA4-like) | E | -4.701 / 3.6e-8 | -2.435 / 2.5e-4 |
| YALI0E16192g | V1 F (VMA7-like), CLIB122 ID; W29 YALI1E19360g is typed "Pseudo mRNA" | no row | — | — |
| YALI1E14125g | V1 G (VMA10-like) | E | -3.334 / 4.9e-3 | -3.144 / 8.9e-9 |
| YALI1E10492g | V1 H (VMA13-like, weak) | E | -3.368 / 4.0e-3 | -1.372 / 0.524 |
| YALI1E12482g | V0 a, isoform 1 (STV1 by synteny, VPH1 by sequence) | E | -5.223 / 8.1e-11 | -2.472 / 1.7e-4 |
| YALI1F38820g | V0 a, isoform 2 (VPH1-like) | not E | -1.749 / 1.0 | -1.516 / 0.313 |
| YALI1F31854g | V0 c (VMA3-like) | E | -4.576 / 1.4e-7 | -2.122 / 6.1e-3 |
| YALI1E37063g | V0 c' (VMA11-like) | E | -4.055 / 2.1e-5 | -2.657 / 1.6e-5 |
| YALI1F13017g | V0 c'' (VMA16-like) | E | -5.876 / 1.4e-14 | -3.194 / 3.9e-9 |
| YALI1F21690g | V0 d (VMA6-like) | E | -4.462 / 4.4e-7 | -1.580 / 0.236 |

**Not in the model but consensus essential (out-of-scope positives):**
- V0 e, YALI1F25663g (VMA9-like). Its two screens disagree: Cas9 -5.234, Cas12a +0.452.
- VMA21-like YALI1B30163g.
- PKR1-like YALI1F37528g.
- RAV1-like YALI1C23873g.

**Calls that rest on the transposon screen:** B, H and d are consensus essential while Cas12a calls them
non-essential, so their consensus must rest on the transposon screen. Its per-gene calls were not recovered.

**D (YALI1D00581g):** it is a borderline negative in both screens. Its label stays as recorded; it is not swapped.

## 6. Biology

All links below were opened by a workflow agent this session, as an abstract or full text as stated. The lead did not
open them.

- **No targeted V-ATPase mutant of *Y. lipolytica* was found.** Essentiality rests only on pooled screens: Cas9 and
  Cas12a ([Ramesh 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10250353/), full text) and Hermes transposon (Patterson
  2018). "Essential" here means strong depletion in a pooled culture, not shown lethality. The only
  *Y. lipolytica*-specific V-ATPase paper is biochemical
  ([Kulakovskaya 1993](https://europepmc.org/article/MED/8385405), abstract).
- ***S. cerevisiae* vma mutants are conditional:** they fail at pH 7.5, with high Ca2+ or Zn2+, and on glycerol
  ([Manolson 1994](https://europepmc.org/article/MED/7514599), abstract).
- **Mechanism of the respiratory failure:** vacuolar cysteine storage is lost, and mitochondrial iron availability
  falls while total cell iron does not fall ([Hughes 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7164368/), full
  text).
- **Iron uptake:** loss of the cwh36 subunit blocks copper loading of Fet3, but only on iron-restricted medium
  ([Davis-Kaplan 2004](https://europepmc.org/article/MED/14594803), abstract).
- **Other organisms:**
  - The V-ATPase is essential in *N. crassa* ([8722770](https://europepmc.org/article/MED/8722770), abstract).
  - The *M. oryzae* vma11 null is viable ([23826342](https://europepmc.org/article/MED/23826342), abstract).
  - So "obligate aerobe, therefore lethal" does not generalise (inferred).
- **Other models (verified by the literature agent):**
  - [Yeast-GEM](https://raw.githubusercontent.com/SysBioChalmers/yeast-GEM/main/model/yeast-GEM.yml) and the Chalmers
    [iYali](https://raw.githubusercontent.com/SysBioChalmers/Yarrowia_lipolytica_W29-GEM/main/model/iYali.yml) use
    full-complex AND rules with no c''-only branch.
  - Both also have free reversible organelle H+ diffusion.
  - No published GEM or constraint method predicting V-ATPase knockouts as lethal was found.

**Proposed experiment (O6):**
- Inhibit the V-ATPase with concanamycin A or bafilomycin on SD-Leu.
- Measure total and mitochondrial iron (ICP-MS).
- Test rescue by iron, by copper, and by cysteine or amino-acid limitation.
- Check the alternative iron-route homologs.

Any such data used to justify a model change becomes curation data. The resulting TP gains would then be a
development reference, not independent validation (AGENTS.md §4).

## 7. Corrections the reviewers made

**Mechanics reviewer:**
- The designer claimed that R2048 in reverse, fed through R1240, could supply Golgi H+ without R794. This is
  refuted, and the lead's re-run confirms it.
- O1 is not a prerequisite for O2a.
- O1's zero-flip result can be proved by the static argument and does not need to be recalled.

**Rules reviewer:**
- AGENTS.md has no list saying "no pump lower bound, no proton sink, no AND". The always-on rule is §2: never make a
  science change to improve a match rate.
- O3 had been called "already rejected". The UM2 TASK.md in fact leaves an evidence-based acidification constraint
  open as follow-up 4.
- The R_NTP1 decision of 2026-09-14 left its GPR unresolved, so reviewing that GPR does not revisit a decision.

**Evidence reviewer:**
- The designer read Hughes 2020 as showing failed iron uptake. That is refuted (section 4).
- Milgrom 2007 and Diab & Kane 2013 do not show respiratory failure.
- Sambade 2005 proposes calcineurin activation as the main basis of the vma growth defect. Ohya 1991 found no
  mitochondrial defect in cls (vma) mutants. Both weaken a single iron explanation.

**Literature agent:** it claimed that knocking out F13017g makes both branches false. This is wrong: F13017g is only
in the second branch, and the lead's Boolean re-run finds 0 of 14 single knockouts closing the rule.

## 8. Side findings (not changed)

- R_NDP7g (Golgi UDP hydrolysis, YALI1E41893g) has an extra proton and charge (verified).
- R2048 has placeholder gene YALIUNK3. Its chemistry, MIPC + Ins-1-P + 2 H+, is not the IPT1 reaction. It is
  blocked only because ER `m674` is a dead end. If a later curation adds an `m674` consumer or sphingolipids to
  biomass, R794 would suddenly become required for an artefactual reason.
- The YALI0E16192g ID should be recorded as unscorable rather than as a silent user-defined negative.
- The V1 A subunit YALI1A09766g also sits on the free cytosolic ATPase R_NTP1. This is a GPR identity question.

## 9. Verification levels

- **Proved by Boolean evaluation, lead re-run:** 0 of 14 single knockouts close the current rule; 12 of 14 close the
  corrected rule.
- **Static, lead re-run:** biomass has no Golgi or vacuolar species; R795 and R1163 are [0, 0]; R2048 is forced to 0;
  R2193 is the sole Fe2+ producer; R1189 is the sole iron boundary.
- **Static, reviewer only:** R_NDP7g runs only in reverse; R141 is dead (also recalled from 2026-09-21).
- **Observed in earlier runs (recalled):** the 2026-09-21 full screen of the corrected rule gave 0 changes; R794 =
  R795 = 0 in every saved wild-type solution; the UM2 openings gave 0 class changes.
- **Untested:** the effect of every option on `059316ef` with a solver.

## 10. Could the method be improved?

- **A structural screen.** For each FN gene, test whether the reactions it closes carry zero flux in every optimum
  (an FVA minimum and maximum of 0 at fixed optimal growth). That would separate "unused at optimum" from "bypassed"
  more cheaply than per-option sweeps. It needs one authorized FVA on the closed set.
- **A test for the static argument.** One authorized LP pair, wild type versus R794 forced to 0, would confirm the
  "pure ATP drain" argument numerically.

## 11. Next gates

- **O4:** the user decides the label wording.
- **O1:** a curation-change run on a branch that holds `059316ef`. The F-subunit ID and the a-isoform compartment
  questions must be settled or documented first.
- **O2a and O3:** blocked until the O6 data exist.
