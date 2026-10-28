# Independent audit: UM7 nuclear diffusion (2026-10-28)

The audit was read-only. Apart from this file, nothing was edited, built or committed. Worktree
`claude/hungry-mendel-f9f7df` at HEAD `0093c89`; all changes are uncommitted (`git status`, `git diff`).

**Inputs:**
- Baseline `…_dolichol_route_v2_20261028.xml`, SHA-256 `9d655c2b…8335c`.
- Candidate `…_dolichol_route_v2_um7_nuclear_diffusion_20261028.xml`, SHA-256 `5f54eb22…54179`.

I re-hashed both files and copied them, byte-identical, to `/private/tmp/claude-501/um7_audit/` (`base.xml`,
`cand.xml`).

**Solves:**
- Task A medium `po1f_csm_leu_w29uptake_20261012.csv` + PO1f profile, through
  `load_effective_simulation_context`.
- Gurobi with `checked_solve.py` (1e-9, reset before each solve), real HOME.
- About 60 LPs in total, from my scripts `scripts/audit_solves.py` and `scripts/ko_route.py` in
  `/private/tmp/claude-501/um7_audit/`.

## Verdicts

| # | Claim | Verdict |
| --- | --- | --- |
| C1 | SHAs; seed-1 build byte-identical; control byte-identical to 9d655c2b | **supported** |
| C2 | Diff is exactly R_ADPRtn + R_PPItn + group membership; both balanced | **supported** |
| C3 | Code matches TASK §8, fails closed, opt-in, requires `--dolichol-route`; data in scope; no `<>&` in notes | **supported** (one benign, documented deviation, finding 1) |
| C4 | WT 0.3116608216 in both, relative change ~5e-15 | **supported** |
| C5 | EGC LPs optimal at 0 in the candidate | **supported** |
| C6 | Blocked sets 705 → 695 / 837 → 828, none newly blocked; UM7 15 → 6 / 15 → 8; R1912 and R576 freed | **supported** (spot-checked) |
| C7 | No class change; identical scores; YALI1F37714g the only change > 1e-6 | **supported** (the "2 × 1,082" wording is slightly off, finding 2) |
| C8 | STOP.md mechanism and controls | **supported**, with one wrong detail (finding 3) |
| C9 | User decisions quoted and recorded; the restated gate matches STOP.md | **partially supported**: internally consistent, but I cannot see the user's messages, so the verbatim quotes are unverified |
| C10 | URLs support the sentences; tiers honest; no *Y. lipolytica* localization presented as experimental | **partially supported** (finding 4) |
| C11 | Escher figure matches README; counts match built_checks.json | **supported** (one wording overclaim, finding 5) |
| C12 | pytest 314 → 320 passed, same 19 failed / 14 errors, only the 6 new tests differ | **supported** |

## Evidence per claim

**C1. Build identity.**
- `shasum -a 256` gives `9d655c2b…8335c` for the baseline and for `builds/control_no_flag.xml`, and
  `5f54eb22…54179` for the candidate and for `builds/candidate_seed1.xml`.
- `cmp` reports both pairs identical.
- The build logs show PYTHONHASHSEED 0 / 1 / 0, "UM7 nuclear diffusion: applied" in the two flagged builds
  only, and exit 0.
- The candidate `.build.json` records the full argv, the output SHA, `source_sha256` and the dirty status.
- All three logs carry the same libSBML "Setting notes on R_R304" error. It is pre-existing: the control build
  reproduces the baseline byte for byte. No report mentions it (finding 7).

**C2. Whole-model diff.**
- A raw `diff base.xml cand.xml` has only 36 lines: two inserted `<reaction>` blocks and two `<groups:member>`
  lines in `G_group101` "Transport, nucleus". No other byte differs.
- The new reactions carry `sboTerm="SBO:0000185"`, `reversible="true"` and the default bounds, with no GPR.
- In cobra:
  - R_ADPRtn is `m1725[C_cy] <=> m818[C_nu]` and R_PPItn is `m203[C_cy] <=> m210[C_nu]`, both
    (−1000, 1000);
  - `gene_reaction_rule ''`, subsystem "Transport, nucleus";
  - `check_mass_balance() == {}` for both.
- The four species are C15H21N5O14P2 −2 (both ends) and H4O7P2 0 (both ends).

**C3. Code, guards and curation data.**
- `um7_nuclear_diffusion.py` fails closed in these cases:
  - a species is missing, or its name, formula, charge or compartment differs from the spec;
  - the two ends differ in formula or charge, or are not in two compartments;
  - the reaction context of m818[C_nu], m1725[C_cy] or m210[C_nu] has changed. Any existing link of either
    pair must touch one of these, so this guard covers the "transport already exists" condition;
  - the group is missing or not unique;
  - an id is taken by another definition;
  - only one of the two reactions is present;
  - the scope or status of the spec is unexpected.
- The flag:
  - `cli.py` adds `--um7-nuclear-diffusion` as `store_true`, so it is opt-in;
  - `main.py` raises "requires --dolichol-route" and applies the step after the dolichol step.
- The curation JSON matches TASK §0: two rows, (−1000, 1000), `gpr ""`, SBO:0000185, group
  "Transport, nucleus", and the baseline SHA `9d655c2b`. A scan of all its strings found no `<`, `>` or `&`.
- The focused test file passes: `…/workflow-automation-review-6765da/platform/.venv/bin/python -m pytest -q
  tests/test_um7_nuclear_diffusion.py` → 6 passed, 24 subtests passed.

**C4. WT growth.** Baseline 0.311660821633612 and candidate 0.3116608216336134. Both are optimal and passed
the check (residual 2.6e-14 / 2.9e-14). The relative change is 4.6e-15.

**C5. Energy-generating cycles.**
- Setup: candidate, every exchange (0, 1000); the only positive lower bound, xMAINTENANCE 7.8625, relaxed.
- Result: R_NTP1 = 0 and biomass_C = 0, both `optimal`.

**C6. Blocked sets.**
- `built_checks.json` gives these lengths:

  | Condition | Baseline | Candidate |
  | --- | ---: | ---: |
  | Open | 705 | 695 |
  | Task A | 837 | 828 |

- `open_newly_blocked` and `taskA_newly_blocked` are both `[]`, and neither new reaction is blocked.
- UM7 lists: the baseline has 15 blocked under both conditions. The candidate has 6 open (R92, R396, R965,
  R1000, R1001, R1055) and 8 on Task A (the same plus R579 and R1270).
- My FVA spot-checks at 0 % (12 reactions × 2 models × 2 conditions):
  - Freed in the candidate: R567, R569, R1912, R576 and R1017 reach 10.67 on Task A (all 0 in the baseline);
    R579 and R1270 reach 1000 open.
  - Still blocked, [0, 0] in both models: R579 and R1270 on Task A; R396, R92, R965, R1000 and R1055 under
    both conditions.
  - R576 is already unblocked open in the baseline ([0, 1000]), so it is freed on Task A only, as claimed.

**C7. Knockout sweep.**
- Recomputed from `ko_compare.tsv`, 1,081 rows:
  - the only |Δ| > 1e-6 is YALI1F37714g, 0.3114518394814089 → 0.3115724887286062;
  - the next-largest |Δ| is 5.1e-14;
  - no class change at any cutoff, and no status change.
- Counts at 15 % (candidate = baseline): 112 TP, 203 FN, 42 FP and 692 TN. Also 16 unresolved (11 positive,
  5 negative) and 16 mtDNA genes out of scope. The summary JSON states that the negatives are the user-defined
  class.
- My re-solve of the YALI1F37714g knockout: baseline 0.3114518394814089, candidate 0.3115724887286062. Both
  are optimal and passed, with max knocked flux 0. The knockout closes 13 reactions (R264, R675–R684, R761,
  R1903), as STOP.md says.

**C8. Mechanism.**
- All four controls, run on the candidate knockout, return the baseline value to within 1e-16:

  | Closed reaction | Knockout growth |
  | --- | ---: |
  | R_ADPRtn | 0.31145183948140853 |
  | R_PPItn | 0.311451839481409 |
  | R567 | 0.31145183948140837 |
  | R1912 | 0.3114518394814087 |

- pFBA in the candidate knockout:
  - R567 = R569 = R1912 = R576 = R1082 = 0.0036456;
  - R_ADPRtn = R_PPItn = R1267 = −0.0036456;
  - R580, R581, R516, R563, R499, R41, R43 and R1376 = 0.0046489;
  - Trp exchange R1232 = −0.00532.
- In the baseline knockout: R567 = R1912 = 0, the kynurenine steps carry 0.0010029 and R1232 = −0.0016737.
- All of these match STOP.md.

**C9. Decisions.**
- TASK §0 quotes "Approve Option A", and so do §10 and the curation `authorization` field.
- TASK §11 quotes "Accept, restate gate". Its three bullets match STOP.md's first option: keep `5f54eb22`,
  record the +1.2e-4, and finish with the audit, REPORT, STATE and README.
- The restated criterion (no class change; YALI1F37714g the only gain above 1e-6) matches the STOP.md data.
- I had no access to the user's messages, so the wording of the quotes is **unverified**.

**C10. Literature.** I opened 12 sources this session through E-utilities `efetch`, UniProt REST and Europe PMC
full-text XML.

Sources that support the sentences they are attached to:

| Source | What it supports |
| --- | --- |
| [Timney 2016, PMC5057280](https://pmc.ncbi.nlm.nih.gov/articles/PMC5057280/) | The quoted definition of passive diffusion appears verbatim in the full text. |
| [Wente & Rout 2010, PMID 20630994](https://pubmed.ncbi.nlm.nih.gov/20630994/) | "Free diffusion of water, sugars, and ions". |
| [Anderson 2002, PMID 11884393](https://pubmed.ncbi.nlm.nih.gov/11884393/) | Npt1 and Nma2 "both concentrated in the nucleus". |
| [Suda 2003, PMID 12898714](https://pubmed.ncbi.nlm.nih.gov/12898714/) | GFP-Qns1 in nucleus and cytosol. |
| [Serrano-Bueno 2019, PMID 30826332](https://pubmed.ncbi.nlm.nih.gov/30826332/) | Ipp1 is nucleocytoplasmic, has an NLS and undergoes nuclear proteasomal degradation. |
| [Tanner 2000, PMID 11106374](https://pubmed.ncbi.nlm.nih.gov/11106374/) | Deacetylation is coupled to 1-O-acetyl-ADP-ribose. |
| [Rafty 2002, PMID 12370179](https://pubmed.ncbi.nlm.nih.gov/12370179/) | YSA1 cleaves OAADPr to AMP + acetyl-ribose-5-P. The nuclear extract activity is "consistent with an OAADPr-specific acetyltransferase", with no enzyme assigned. |
| [Tong 2009, PMID 19251690](https://pubmed.ncbi.nlm.nih.gov/19251690/) | Ysa1 cleaves ADPr / OAADPr. |
| [Wilson 2006, PMID 17110954](https://pubmed.ncbi.nlm.nih.gov/17110954/) | Hst2 is largely cytoplasmic through Crm1 export. |
| [Q12726](https://rest.uniprot.org/uniprotkb/Q12726.txt) | Transit peptide ECO:0000255 (predicted). |
| [Q06178](https://rest.uniprot.org/uniprotkb/Q06178.txt) | Nma1 in cytoplasm and nucleus, ECO:0000269. |

Every *Y. lipolytica* compartment is labelled as predicted or IBA. Partial support is covered in finding 4.

**C11. Escher figure.**
- The counts in `escher/index.json` match README. The in-model counts are 21 → 23. The other counts, recounted
  by me from `built_checks.json` on the 23 mapped reactions:

  | Count (23 mapped reactions) | before | after |
  | --- | ---: | ---: |
  | Blocked, Task A | 18 | 9 |
  | Blocked, open | 16 | 6 |
  | Capacity > 0 on Task A | 3 | 14 |

- The status lists in `index.json` (freed, blocked in both, closed, new) equal README's.
- The jsdom and Chrome logs show 23 reactions and 125 nodes in all five scenarios. The R_ADPRtn label reads
  "(nd)" before and "10.7" after.
- `flux.log` shows both WTs, uracil R1354 at (−1000, 1000), R612 at (0, 0) and 272 LPs.
- I viewed the PNG. It matches the description: capacity before is grey, after is blue, with 10.7 on the
  sirtuin route.
- README states honestly what was not checked: tooltips, and the viewer in an interactive browser.

**C12. pytest.**
- Logs run with the workflow-automation venv, Python 3.13.5 / pytest 9.0.2:

  | Log | Summary line |
  | --- | --- |
  | `pytest_before.log` | 19 failed, 314 passed, 31 skipped, 14 errors |
  | `pytest_after.log` | 19 failed, 320 passed, 31 skipped, 14 errors |

- A sorted test-by-test diff of the PASSED / FAILED / ERROR / SKIPPED lines shows only the 6 added `PASSED
  tests/test_um7_nuclear_diffusion.py::…` lines.
- STOP.md's caveat stands: the before run overlapped the main.py / cli.py edits.

## Findings not covered by the claims, or limiting them

1. **Benign deviation from TASK §8 (C3).**
   - §8 says the step fails closed if "a reaction id already exists".
   - The code instead accepts an existing reaction whose definition is identical (status `already_correct`,
     which makes reruns and a reloaded SBML idempotent). It fails only on a different definition.
   - The module docstring and a test document this. Worth one line in REPORT.md.
2. **Gene count (C7).**
   - The model has 1,082 genes. `PO1f_plasmid_LEU2` is a runtime-excluded gene, so the sweep covers 1,081
     genes per model; the TSV has 1,081 rows and the summary JSON says `"genes": 1081`.
   - The LP count 2,164 is 2 × (1,081 knockouts + 1 WT).
   - TASK §8 ("2 × 1,082") and REPORT.md (table row "2 × 1,082") should say 1,081 knockouts.
3. **STOP.md step 4 names the wrong transport for ATP.**
   - STOP.md says "Gln and ATP via R1214 / R1017". R1017 is AMP transport (`m209[C_nu] → m86[C_cy]`); in the
     knockout pFBA it exports the AMP made by R569 (+0.0036456).
   - ATP enters the nucleus by R1021 (`m141[C_cy] → m206[C_nu]`, +0.0036456).
   - Also active: R1212 (Glu out), R1161 (H+) and R1361 (H2O).
   - The conclusion of the mechanism is unaffected.
4. **Attribution (C10).**
   - **Hst1 and Hst4.** TASK §3, the curation `evidence` row and REPORT.md say "Sir2, Hst1 and Hst4 are
     nuclear (UniProt ECO:0000269)" but link only [P06700](https://rest.uniprot.org/uniprotkb/P06700.txt),
     which covers Sir2 only (Nucleus, nucleolus, ECO:0000269, PubMed 9214640).
     - Hst1 and Hst4 are in `stage1/evidence.md` as P53685 and P53688.
     - Hst1's nuclear call is the large-scale GFP screen (Huh 2003, PubMed 14562095), and its note adds "and
       cytoplasmic" ([P53685](https://rest.uniprot.org/uniprotkb/P53685.txt), opened this session).
     - Cite P53685 and P53688, and say that Hst1 is high-throughput only.
   - **Ysa1 mitochondrial by IDA (Tong 2009).** The Tong abstract does not state a localization. SGD's GO
     backend lists "mitochondrion — Tong L et al. (2009)" (opened this session), but my parse did not show the
     evidence code. So the claim rests on SGD, not on the PubMed link.
   - **Evidence-agent verification.** TASK §3 says the evidence agent fetched every URL and the author
     re-opened three. "verified" for the rest is therefore the agent's reading. The audit re-read 12, listed
     above.
5. **pFBA wording overclaims.**
   - TASK §5 and escher/README.md say that pFBA carries "no flux through any UM7 reaction". R568 and R580,
     which TASK §2 lists among the 17 UM7 reactions, carry 0.0010035 at the optimum in both models.
     README's own table counts 2 carrying reactions.
   - The wording should be "no flux through any of the 15 blocked UM7 reactions or the new links".
6. **Stale or inconsistent text in TASK.md.**
   - §3 cites `um7_stage1_20261028/evidence.md`; the file is `stage1/evidence.md`.
   - §8 step 4 names the output `…_um7_diffusion_<date>.xml`; the built file is
     `…_um7_nuclear_diffusion_20261028.xml`.
   - §1 says "libSBML cannot open the candidate through its full path (> 255 characters; the loader resolves
     symlinks)". escher/README.md and the test attribute the problem to cobra's `"win" in platform` check.
     Settle on one explanation.
   - The §0 recommendation paragraph, "In memory there is no change … to any knockout class", is true. But it
     sits above §11 without noting that one knockout gains growth. A pointer to §11 would help.
   - §1 says the worktree was "clean" at the start. That was true then, but it is now dirty.
7. **Pre-existing R304 notes error.**
   - Every build, including the byte-identical control, logs libSBML error −3 while setting the notes on
     R_R304.
   - This change did not cause it, but the reports do not mention it. One line saying it is pre-existing would
     prevent a later reader from attributing it to this step.
8. **Annotations.**
   - The existing nuclear diffusions carry BiGG / MetaNetX cross-references (for example R1082:
     `bigg.reaction/DNADtn`, `MNXR97625`). R_ADPRtn and R_PPItn carry only SBO:0000185.
   - This is outside the authorized scope; record it as a known gap rather than change it.
9. **Bookkeeping.**
   - `model/candidates/README.md` already has the new row, links to `REPORT.md` and says "not committed".
   - REPORT.md appeared during this audit (15:46) with an `__AUDIT__` placeholder.
   - No `model/STATE.md` entry exists yet.
   - The `ko_compare_two.py` docstring still describes the 2c9cebd0 / lipoylation run it was copied from.
10. **Scope.** Not covered:
    - essentiality is checked on Task A only;
    - the mechanism finding (an unconstrained lumped R567 acting as an NAD hydrolase once its products have
      outlets) is a model property; whether the PNP knockout's real growth changes is untested, as STOP.md says;
    - the audit did not rerun the whole-model FVA or the sweep: it used spot-checks and one knockout.
