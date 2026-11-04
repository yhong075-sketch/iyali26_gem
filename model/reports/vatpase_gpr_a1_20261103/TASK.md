# V-ATPase R794/R795 rule: one core AND either a-subunit, plus V0 e (option A1 + E-1), and Scope-1 (2026-11-03)

## 0. Authorization

The decision write-up `vatpase_decision_options_20261102.md` (shared session notes of the V-ATPase discussion
session) set out the options. On 2026-11-03 the user answered "ok, do that" to the recommended line. The E-2 guide
check then ran (read-only; `$IYALI26_RESEARCH_ROOT/artifacts/tasks/vatpase_gpr_blast_20261101/e_guides/`). A plan
followed, and the question "Run the V-ATPase round with this plan?" got the answer:

> "Yes, with E-1 + Scope-1"

The approved choices are:

| Code | Choice |
| --- | --- |
| Std-§5 | AGENTS.md §5 decides the a-subunit OR |
| Fgate-lift | YALI0E16192g may be a required AND member, as a homology-based candidate |
| A1 | The same rule on R794 and R795 |
| Label-1 | The structural-FN label goes in REPORT and STATE only |
| Scope-1 | YALI0E16192g goes on the untargeted list (named explicitly) |
| ID-keep | Keep the CLIB122 ID |
| E-1 | Add V0 e YALI1F25663g |
| Verify-1 | Build, whole-model diff, WT and 15 target knockout LPs |

Verify-2 and Verify-3 were not chosen. No commit, push or PR without a separate instruction.

## 1. Starting state (verified this session)

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/vatpase-gpr-a1-20261103`, new from `origin/main` `b57272d3`; clean tree |
| Baseline | `candidates/E5_…_treh_biomass_r24_aco_ceramide_pool_20261101.xml`, SHA-256 `2bc263f23603c27131d0b07143300a484dbd5a135ef7151665574400ab64a062`. Not a designated reference. PRs #45 and #46 (after `80cfc4f`) changed no model or source XML. |
| Its build | `/opt/miniconda3/bin/python`, `PYTHONHASHSEED=0`, 37 opt-in flags (argv in its `.build.json`) |
| Before-state | [before.json](before.json) (`scripts/before_state.py`): R794 and R795 in full, their 14 genes, and `YALI1F25663g` absent (1,081 genes) |
| R794 / R795 rule now | `(YALI1A11258g and YALI1E14125g and YALI1D00581g and YALI1E10492g and YALI1A09766g and YALI1F20965g and YALI0E16192g and YALI1F31854g and YALI1F21690g and YALI1E37063g and YALI1E32332g and YALI1E12482g) or (YALI1F38820g and YALI1F13017g)`; R794 [0, 1000], R795 [0, 0] |
| 15 % score on the baseline | 142 / 184 / 44 / 694 (rule B, mtDNA + untargeted out of scope; the planning session's recount of main's sweep) |

## 2. The change

**R794 and R795** both get:

`YALI1A11258g and YALI1E14125g and YALI1D00581g and YALI1E10492g and YALI1A09766g and YALI1F20965g and YALI0E16192g and
YALI1F31854g and YALI1F21690g and YALI1E37063g and YALI1E32332g and YALI1F13017g and YALI1F25663g and (YALI1E12482g or
YALI1F38820g)`

| Part | Genes | Change |
| --- | --- | --- |
| V1 head | A A09766, B E32332, C A11258, D D00581, E F20965, F YALI0E16192g, G E14125, H E10492 | unchanged AND members |
| V0 ring | c F31854, c' E37063, d F21690 | unchanged AND members |
| V0 c'' | F13017 | moves from the headless second branch into the core AND |
| V0 e | F25663 (new gene) | added as an AND member |
| V0 a | a1 E12482 or a2 F38820 | the only OR in the rule |

**What goes:** the headless second branch `(F38820 and F13017)`. It has no ATP-hydrolysing head and no c/c' ring.

**The new gene:**
- `YALI1F25663g` gets the cross-references verified in the local NCBI feature table (`reference/ncbi/clib89_feature_table.txt`):
  - NCBI Gene 90949792;
  - RefSeq XP_065950526.1;
  - SBO:0000243.
- It gets no name, UniProt or KEGG ID: none is verified, and there is no CLIB122 ID.

**Reaction notes** record:
- the curation ID and the authorization;
- the previous rule;
- the standard used: it meets AGENTS.md §5 and deliberately does not meet the 2026-09-21 same-compartment test;
- the evidence and its tier.

**Unchanged:** stoichiometry, bounds (R795 stays [0, 0]), annotations, every other reaction and every other gene. That
includes R_NTP1, which still carries YALI1A09766g.

**Builder step:** opt-in flag `--vatpase-a1-gpr`. It is refused together with `--vatpase-gpr-hypothesis` (the
2026-09-11 provisional hypothesis).

**Evidence:**
- [2026-10-30 options report](../vatpase_fn_options_20261030/REPORT.md) and the 2026-11-01 BLAST round
  (`$IYALI26_RESEARCH_ROOT/artifacts/tasks/vatpase_gpr_blast_20261101/`, `synthesis.md`, `review_*/`).
- **Every rule gene is a V-ATPase subunit:** 13 of the 15 are reciprocal best hits of the *S. cerevisiae* subunits (H
  weakest at 26.8 %). The exceptions: a1 E12482, whose best yeast hit is VPH1 (876 vs 811 bits over STV1); and V1 F,
  which has no W29 protein (CLIB122 Q6C5Q2's best hit is VMA7, one-way). Corrected after the audit (it said "all").
- **c'' null loses all activity:** Hirata 1997, PMID 9030535.
- **Both a isoforms form active pumps:** Manolson 1994, PMID 7514599; Kawasaki-Nishi 2001, PMID 11592965;
  Perzov 2002, PMID 11948198. That is the §5 OR test.
- **Same-compartment substitution is only partial and one-way:** Deschamps 2020, PMC7002768. That is why the
  2026-09-21 test is not met.
- **V0 e:** reciprocal best hit of Vma9 at 61.4 %. The e null fails assembly (Compton 2006, PMID 16569636).
  Against: once assembled, e is not needed for coupled pumping in vitro (Bueler 2015, PMID 25546637).
- **The screens agree on e (E-2):** Cas9 has four guides that cut well (CS 5.3–7.2), all strongly depleted. Cas12a's two
  guides did not cut (CS −1.60, −0.36).
- **F subunit:** intact when spliced (missed 2-nt first exon plus a 56-nt intron; inferred, no RNA evidence). This is
  the Fgate-lift evidence.

**Not reasons:** recall, FN or FP counts, and the screen calls of E12482 and F38820.

### Scope-1 (benchmark files)

- `model/expected/untargeted_genes_out_of_scope.json`:
  - add `YALI0E16192g` to `genes`;
  - add its evidence: W29 locus NC_090774.1, no library guide in 1935676-1936100, no guide ID naming E19360/E16192.
    Source: `vatpase_gpr_blast_20261101/fsubunit/out/run7_guides_in_F_region.tsv`, re-run 2026-11-02 byte-identical.
- `model/expected/benchmark.md`: the untargeted sentence names both genes.

## 3. Expected effect (inferred, to be measured by Verify-1)

| Knockout | Rule result | Expected growth |
| --- | --- | --- |
| Each of the 13 core genes | closes R794 (R795 is already closed) | stays at WT, because R794 is a pure ATP drain in this network (static argument, inputs verified 2026-11-02) |
| a1 or a2 alone | closes nothing | stays at WT |
| WT | rule does not enter the WT LP | unchanged, 0.3120535436 |

**Expected counts at 15 %:**
- F25663 becomes an in-model positive scored FN: FN 184 → 185, in-model positives 326 → 327.
- YALI0E16192g leaves the TN count: TN 694 → 693.
- TP 142 and FP 44 unchanged. Expected 142 / 185 / 44 / 693.

**Not measured:** the 1,067 other knockouts (1,050 in scope) are not re-run (Verify-3 not chosen). Their unchanged classes are
inferred: the rule only adds a closure of R794 to these 15 knockouts.

## 4. Acceptance criteria

1. A control build with 2bc263f2's argv at HEAD (code in place, flag off) reproduces `2bc263f2…` byte for byte.
2. The candidate built with `PYTHONHASHSEED=0` and `=1` gives identical bytes.
3. The whole-model diff against 2bc263f2 shows only:
   - R794's and R795's rule and notes;
   - the new gene `YALI1F25663g` with its annotation;
   - the reaction links of the 15 rule genes.
4. Verify-1 (16 LPs on the candidate):
   - WT optimal and passing the targeted check;
   - each of the 15 knockouts optimal and passing it;
   - R794 closed exactly for the 13 core knockouts.
   Growth that differs from WT is reported as found, not tuned.
5. Full pytest: only the new tests differ from the before run, test by test.

## 5. Plan and budget

1. **Code:**
   - data file `model/curation/gpr/vatpase_a1_gpr.json`, written by `scripts/write_curation.py`;
   - step `platform/scripts/gem_annotate/vatpase_a1_gpr.py` (`apply_vatpase_a1_gpr`), which:
     - pins both reactions' name, bounds, stoichiometry, species, annotation, rule and notes;
     - pins the 14 genes' annotation and reactions, and checks that the new gene is absent;
     - fails closed and rejects a partly applied state;
     - adds the one new gene through the rule;
   - flag in `cli.py`; call in `main.py` after the R24 step;
   - test `platform/tests/test_vatpase_a1_gpr.py`.
2. **Full pytest** before (origin/main `b57272d3` by `git archive`, no edits) and after, with the
   `workflow-automation-review-6765da` venv (no scipy).
3. **Builds** from `platform/` with `/opt/miniconda3/bin/python`:
   - control (seed 0, scratch);
   - seed check (seed 1, scratch);
   - candidate (seed 0) at `candidates/E5_…_ceramide_pool_vatpase_a1_20261103.xml`.
4. **Checks:**
   - whole-model diff (`scripts/model_diff.py`, copied from the ceramide round);
   - Verify-1 with `scripts/ko_targets.py`: the `checked_solve.py` of the earlier rounds, Task A medium + PO1f, Gurobi
     1e-9, reset per knockout.
5. **Write-up and Scope-1:**
   - Scope-1 file edits;
   - independent read-only audit;
   - REPORT.md, including the 2026-11-01 EBI e-mail disclosure;
   - `model/STATE.md` entry and a `candidates/README.md` row.
   Stop before commit.

**Compute budget:**
- 16 LPs for Verify-1;
- at most 40 for the audit;
- 3 builds;
- 20 minutes of solver time.

**Stop conditions** (report, never re-tune):
- the control does not reproduce 2bc263f2;
- the two seeds differ;
- the diff shows anything unexpected;
- an LP is not optimal or fails the check.
