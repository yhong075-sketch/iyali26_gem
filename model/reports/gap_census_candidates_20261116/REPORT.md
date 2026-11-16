# Inactive reactions of candidate 7ccbb568: causes, and LP tests of single-reaction candidates per root (2026-11-16)

Scope and authorization: [TASK.md](TASK.md). Computation and reporting only; no model change. Committed on the
user's later instruction; not pushed.
Small tables are in this folder; scripts in [scripts/](scripts/); raw outputs, logs and full tables in
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/gap_census_candidates_20261116/`.

## 1. Answer in short

- **Inactive reactions today.** With every exchange open, 692 of the 2,347 reactions of `7ccbb568` cannot carry
  flux (29.5 %). GLPK and Gurobi (1e-9) return the same 692. The October census found 756 of 2,316 on `a1237821`.
- **Why.** One cause per reaction, assigned in the October order:

  | Group | Blocked reactions | Which of the two questions |
  | --- | ---: | --- |
  | Dead end: something is missing | 523 | question 2 |
  | Closed by bounds | 82 | neither: a closure to review |
  | No dead end found (coupling) | 49 | not decided by the census |
  | Steady state: a cofactor is made but never drained | 38 | question 1 |

- **Candidates.** MetaNetX offers a balance-marked single-reaction candidate for 81 of the 364 root metabolites.
  At 60 roots at least one candidate works on its own. At 80, at least one works once the other dead ends are
  opened.
- **Pre-set reading.** 57 roots have two or more working, balance-marked candidates in the second pass (42 alone).
  The cut set before the run was 30; the ceiling known before the run was 58.
- **The LP does not tell candidates apart.** 243 of the 250 balance-marked root–reaction pairs pass in the second
  pass, and 181 pass alone. At D-ribulose all seven candidates pass and each frees the same single reaction.
- **Some passes are only loops** (found by the independent review). A candidate also passes when it merely
  closes a cycle with a reaction the root already has: a kinase next to the phosphatase of the same compound
  burns ATP, and the compound still has no source. The reviewer's proxy flags 19 balance-marked pairs as possibly
  of this kind; without them the reading count is 51.
- **What a single added reaction cannot reach.** 283 roots have no balance-marked candidate at all, and 252
  blocked reactions sit behind several roots at once.
- **A side effect found on the way.** Nine reactions became blocked since October without a bounds change: the
  methylcitrate cycle (R25, R552, R95), the methylmalonate route (R2074–R2077, R2112) and the CoA transport
  R834. Reopening the HMG-CoA transport R1172 in memory frees all nine; R1172 was closed on 2026-11-12
  (section 3.3).

## 2. Run

- **Model.** `candidates/E5_353cefdc_xref_20261115.xml` (`7ccbb568…`, full hashes in [TASK.md](TASK.md)), the
  newest candidate on `main`. Not a designated reference.
- **Conditions.** Primary: every exchange open at (−1000, 1000). Secondary: Task A medium + PO1f profile.
- **Solvers.** cobra 0.30.0 with GLPK for the census scripts and the candidate tests; Gurobi 13.0.2 at 1e-9 as a
  second solver for the blocked sets and for every candidate pass-or-fail solve.

| Gate | Result |
| --- | --- |
| File hashes (model, regression model, medium, profile) | equal to TASK.md |
| Census scripts | 8 files byte-identical to the October copies |
| Task A condition file against the platform's `load_effective_simulation_context` | 0 of 2,347 reaction bounds differ |
| Regression on `a1237821` | 756 blocked, 701 gap metabolites, 212 modules; summary and blocked set identical to October |
| Root table of `7ccbb568` | 364 (RNP 98, RNC 101, single reversible 137, isolated 28), equal to the plan-time count |
| Candidate list | 391 pairs in 400 placements, and the per-root status counts, equal to the plan-time count |
| Candidate implementation checks | 20 failing candidates re-counted: none freed a reaction. Every candidate that works alone freed at least one reaction at its root |

- **How the gates were enforced.** In the first run the hashes were checked by shell commands, and the reviewer
  re-hashed the model, the MetaNetX table and both gap-fill tables. Since the follow-up (section 8) both new
  scripts check the hashes themselves against [input_hashes.csv](input_hashes.csv): `census_tables.py` 12 inputs,
  `candidate_tests.py` 4. The LP cap was checked by tally after each stage, not by code.
- **Re-run.** The candidate stage was run a second time with the hash gate in place. Its three output tables are
  byte-identical to the first run's.
- **LPs.** 102,746 in all, under the cap of 150,000. First run 96,656: check 1 90,620
  ([lp_count_check1.csv](lp_count_check1.csv)), check 2 6,036 (GLPK 4,368 including the re-counts, Gurobi
  1,668). Follow-up 6,090: 54 for the closure test of section 3.3 and 6,036 for the re-run. These are nominal
  counts, two per re-tested reaction; the retried drain test adds 1,384 and solver pre-solves are not counted.
- **Time.** First run: 30 minutes of wall-clock time (13:46 to 14:16 PST), with runs in parallel; 57 minutes of
  process time. Follow-up: 6 minutes (15:02 to 15:08 PST).
- **Two incidents.**
  - The first launch of the regression scan stopped after 1 second with `FileExistsError`: two scans started
    together and both tried to create the same fresh cache folder. It was re-run unchanged and passed.
  - One of the 20 drain tests (NAD) used the census scripts' built-in single retry without presolve. Its result
    is consistent with the others: the union of the 19 single drains equals the joint drain (38 reactions).

## 3. The census on 7ccbb568

### 3.1 Counts

| Quantity | 7ccbb568 open | 7ccbb568 Task A + PO1f | a1237821 open (October) |
| --- | ---: | ---: | ---: |
| Reactions | 2,347 | 2,347 | 2,316 |
| Blocked reactions (GLPK) | 692 (29.5 %) | 829 (Gurobi 826) | 756 (32.6 %) |
| of which closed (0,0) | 76 | 77 | 69 |
| of which exchanges or demands | 39 | 89 | 42 |
| of which carry a gene rule | 427 | 479 | 472 |
| Metabolites | 1,883 | 1,883 | 1,882 |
| Gap metabolites | 649 | 730 | 701 |
| Roots: RNP / RNC / single reversible / isolated | 98 / 101 / 137 / 28 | 152 / 101 / 124 / 28 | 100 / 103 / 136 / 28 |
| Blocked reactions explained by root propagation | 604 | 693 | 618 |
| Unconnected modules | 213 | 260 | 212 |
| Largest modules | 121, 37, 32, 30, 19, 17 | 121, 37, 32, 30, 19, 17 | 127, 37, 37, 32, 30, 19 |

Of the 213 modules, 122 hold one reaction, 40 hold two, 32 hold three to five, 9 hold six to ten, 9 hold 11 to 50,
and one holds 121 (complex sphingolipids). Tables: [summary_counts.csv](summary_counts.csv),
[blocked_reactions.csv](blocked_reactions.csv), [root_metabolites.csv](root_metabolites.csv),
[unconnected_modules.csv](unconnected_modules.csv).

Module labels (UM1, UM2, …) belong to one run. The census script numbers modules of equal size in an arbitrary
order, so today's labels do not match October's or those in `model/STATE.md`: the nuclear RNR module was UM11 in
October and is UM9 here. Compare modules by their reactions, not by label.

### 3.2 One cause per blocked reaction

| Group | Cause (October order) | 7ccbb568 | a1237821 (October) |
| --- | --- | ---: | ---: |
| Closed by bounds | closed in the model (0,0) | 76 | 69 |
| Closed by bounds | blocked by a closed reaction (opening it frees this one) | 6 | 18 |
| Dead end | chain disconnected at several roots (no single fix frees it) | 252 | 290 |
| Dead end | dead-end root metabolite (a sink on one root frees it) | 190 | 181 |
| Dead end | missing transport (a link to the active copy elsewhere frees it) | 58 | 52 |
| Dead end | missing exchange | 17 | 17 |
| Dead end | duplicate metabolite id (a link to a same-named species frees it) | 6 | 0 |
| No dead end found | stoichiometric coupling, no root in the module | 49 | 57 |
| Steady state | cofactor pseudo-gap (a cofactor drain frees it) | 38 | 72 |
| | Total | 692 | 756 |

- **Question 1, "muted by steady state, nothing missing".** The clean case is the cofactor group: 38 reactions
  whose pathway is complete, but whose product is a cofactor that nothing drains. The 2013 paper's remedy is a
  biomass term or a drain.
- **Question 2, "something is missing".** The dead-end group: 523 reactions that follow from a root metabolite.
  604 of the 692 blocked reactions follow from the roots without any solver.
- **Not covered by the two questions.** 82 reactions are blocked by a (0,0) bound, their own or a neighbour's.
- **Not decided.** The 49 coupling reactions have no dead end that the per-compartment scan can see. That does
  not show that nothing is missing: 11 of them sit in a module with a pooled root (4-hydroxy-2-oxoglutarate is
  never consumed in any compartment, glycolate is never produced in any), and the October report describes
  coupling modules that lack an entry or an exit step (its UM15, UM21, UM23).
- **A correction.** Earlier in this conversation the 57 coupling reactions of October were counted under
  "steady state, nothing missing" (129 = 72 + 57). The October report's own module notes do not support
  "nothing missing" for that group; only the 72 cofactor reactions were the clean case.

### 3.3 What changed since October

[cause_changes_vs_a1237821.csv](cause_changes_vs_a1237821.csv).

- **31 reactions were added and all are active.** None was removed.
- **78 became active.** By their October cause: cofactor pseudo-gap 38, behind a closure 15, several roots 10,
  coupling 7, dead-end root 7, closed 1. `model/STATE.md` records the steps: 756, then 722 after the cofactor
  biomass terms (`00e8e397`), then 695 (`14f73bda`), and 692 now. The 78 were not traced to individual builds here.
- **14 became blocked.**
  - 5 were closed by later builds: R1172, R1771, R1846, R2121, R353.
  - 9 kept their bounds: R2074, R2075, R2076, R2077 and R2112 (methylmalonate route), R25, R552 and R95
    (methylcitrate cycle), and the CoA transport R834. Each is now freed by the FAD drain and also by the CoA
    drain, so the census files them under the cofactor group. That label is not the reason they became
    blocked: the census re-tests an opened closure only inside its own module, and the closures made since
    October sit in other modules.
  - **The reason is the R1172 closure** (follow-up test, 54 LPs, GLPK, the census `fixtests.py` unchanged). With
    R1172, the HMG-CoA transport closed on 2026-11-12, reopened in memory at its October bounds, all nine carry
    flux again. With R2121 reopened instead, none does; with both, all nine.

    | Reopened in memory | Of the nine, free again |
    | --- | ---: |
    | R2121 | 0 |
    | R1172 | 9 |
    | R2121 and R1172 | 9 |

    `model/STATE.md` (2026-11-12) already records R834 as the CoA return of the R1172 route: "R834 only returned
    the CoA that left with HMG-CoA". The other eight are not mentioned there. A reading that fits all three facts
    (the test, the note on R834, and the CoA drain that frees them) is that R1172 with R834 was the model's only
    net exchange of CoA between the two compartments; that reading is inferred, not tested. This is LP evidence
    about a side effect, not a proposal to reopen R1172.
- **41 stayed blocked with another cause.** The largest shift is 23 reactions from "several roots" to "dead-end
  root".

### 3.4 Fix tests and cofactor drains

| Kind of fix | Tests | Tests that free at least one reaction | Largest count in one test |
| --- | ---: | ---: | ---: |
| Free sink on a root | 364 | 151 | 26 |
| Link to the active copy in another compartment | 133 | 46 | 5 |
| Opening one closed reaction | 76 | 22 | 4 |
| Link to a same-named species in the compartment | 55 | 5 | 9 |
| Exchange for an extracellular root | 6 | 6 | 5 |

- In October no duplicate-id link freed anything (0 of 55). Now 5 do, all in the sphingolipid module UM1.
- Cofactor drains ([cofactor_drains.csv](cofactor_drains.csv)): FAD 15, CoA 11, PLP 8, siroheme 3, NADP 3, THF 2,
  FMN 2, biotin 2, NAD 1, riboflavin 1, pantetheine 4'-phosphate 1; all 19 together 38. In October the same
  tests freed 72.

### 3.5 Task A medium, and the two solvers

- On the Task A medium + PO1f, 829 reactions are blocked by GLPK and 826 by Gurobi. All 692 of the open
  condition are among them; 137 more are blocked only on the medium by GLPK, 134 by Gurobi.
- The three reactions GLPK adds (R1343, R1344, R970) carry 3.7e-7 in Gurobi: the thiamine tolerance fringe that
  the October report described with five reactions. [solver_check.csv](solver_check.csv) lists them; the other
  tables count them among the 829 without a flag.
- With every exchange open the two solvers give the same set.

## 4. The candidate tests

### 4.1 Which roots have a candidate at all

Each root is counted at the first rule that leaves it with nothing
([root_candidate_counts.csv](root_candidate_counts.csv)).

| Root status | Roots |
| --- | ---: |
| In no open reaction (only closed reactions touch it) | 28 |
| No MetaNetX id on the metabolite | 61 |
| No MetaNetX enzyme reaction uses the compound | 81 |
| No reaction that MetaNetX marks as balanced | 35 |
| Every balanced reaction needs a compound the compartment lacks | 58 |
| The only candidates are already reactions of the compartment | 20 |
| One candidate | 23 |
| Two or more candidates | 58 |
| Total | 364 |

- This leaves 250 balance-marked root–reaction pairs at 81 roots, and 141 unmarked pairs.
- 85 further pairs equal a reaction already in the compartment and were not tested
  ([candidates_identical_to_model.csv](candidates_identical_to_model.csv)). In 78 an open copy exists. In 7 the
  only copy is closed (R17 for three roots, R92, R796, R1027, R1899): those are closure questions.

### 4.2 What the LP says

Roots by the number of candidate reactions, balance-marked only ([candidate_summary.csv](candidate_summary.csv)):

| | None | One | Two or more |
| --- | ---: | ---: | ---: |
| Candidates listed | 283 | 23 | 58 |
| Work alone | 304 | 18 | 42 |
| Work with the other roots opened | 284 | 23 | 57 |

By root class, balance-marked only:

| Root class | Roots | With a candidate | One works alone | Two or more work alone | Two or more work, others opened |
| --- | ---: | ---: | ---: | ---: | ---: |
| RNP (never produced) | 98 | 24 | 18 | 9 | 12 |
| RNC (never consumed) | 101 | 25 | 13 | 10 | 19 |
| Single reversible reaction | 137 | 32 | 29 | 23 | 26 |
| In no open reaction | 28 | 0 | 0 | 0 | 0 |
| Total | 364 | 81 | 60 | 42 | 57 |

("One works alone" counts roots with at least one.)

- **Pairs.** Of the 250 balance-marked pairs, 181 work alone and 243 work with the other roots opened.
- **Both solvers agree** on all 834 pass-or-fail calls (417 placements, two passes), and every solve ended optimal.
- **Reactions freed.** 104 distinct blocked reactions are freed by at least one balance-marked candidate acting
  alone. At 44 of the 77 roots with a working candidate of any mark, the best candidate frees exactly one
  reaction; the largest count is 13. These counts are lower bounds: only the root's own module was re-counted,
  and the reviewer found 34 working placements at 19 roots whose flux passes through a gap metabolite of
  another module.
- **The reading set.** The 57 roots are in the cytosol (26), the mitochondrion (15), the extracellular space
  (6), the ER (3), the cell envelope (3), the ER membrane (2), the lipid particle (1) and the mitochondrial
  membrane (1). None is under a recorded decision. At 34 of them a free sink on the root also frees a reaction;
  at 11 the generic transport link does.

### 4.3 What a pass does and does not mean

The LP only asks whether the new reaction connects the dead end to something that can balance it. Almost any
reaction built from compounds of the same compartment does. D-ribulose in the cytosol (`m2061[C_cy]`, a
single-reaction dead end) shows it: all seven balance-marked candidates work alone, in both directions, and each
frees the same one reaction.

| MetaNetX reaction | Equation (model names) |
| --- | --- |
| MNXR162445 | D-ribulose + ATP ⇌ D-xylulose 5-phosphate + ADP + H+ |
| MNXR197755 | D-ribulose + ATP ⇌ D-ribulose 5-phosphate + ADP + H+ |
| MNXR162615 | D-ribulose + phosphate ⇌ D-ribulose 5-phosphate + H2O |
| MNXR169705 | D-xylulose ⇌ D-ribulose |
| MNXR197816 | D-ribose ⇌ D-ribulose |
| MNXR198291 | D-ribulose ⇌ D-arabinose |
| MNXR197986 | D-arabinitol + NADP+ ⇌ D-ribulose + NADPH + H+ |

Which of the seven *Y. lipolytica* really has is a question of evidence, not of the LP.

- **Direction.** Of the 190 balance-marked placements that work alone, 86 can run both ways, 63 only as written
  and 41 only in reverse. Nothing here judges whether a direction can run. Mitochondrial formate, for example,
  passes with a balanced MetaNetX equation that makes formate and oxygen from carbon dioxide and water.
- **Most unmarked reactions are not usable.** A static side check with the model's own formulas (no LP, hydrogen
  ignored) finds that 124 of the 141 unmarked placements do not balance; "L-proline ⇌ H2O" is one of them. 14
  balance and 3 use a species without a formula. Of the 259 balance-marked placements, 250 balance and 9 use a
  species without a formula. With that check added, the reading count is 56 instead of 57.
- **Compartment.** A pass says nothing about whether the enzyme is in that compartment. 50 placements copy a
  reaction the model already has in another compartment; 19 of them work alone.
- **Loops.** A candidate also passes when it only closes a cycle with a reaction the root already has. Two cases,
  traced by the reviewer and re-computed here from the stoichiometry:
  - 4-hydroxy-L-threonine and its 4-phosphate in the cytosol (`m136`, `m135`) are both dead ends of one
    phosphatase, R51. The candidate kinase MNXR95023 passes at both and frees only R51. Together the two
    reactions net to ATP + H2O → ADP + phosphate, and neither compound gains a source.
  - Mitochondrial nicotinate (`m826`): MNXR188489 passes only in reverse, and with R575 it nets to
    ADP + phosphate → ATP. That is an energy-generating cycle, which cannot exist, and the LP accepts it.

  Two proxies size the problem. Neither is exact.
  - *Same reaction up to currency compounds* (the reviewer's rule): the candidate and a reaction at the root have
    the same substrates and products once 24 currency compounds (H+, water, nucleotides, NAD(P)(H), FAD, CoA,
    CO2, O2 and the like) are removed. 19 balance-marked pairs match, 17 of which pass alone. One of the 19 is
    matched against a closed reaction and should not count. Without the 19 the reading count is 51.
  - *A free sink on the root frees nothing* (from the census fix tests), and the candidate still passes alone.
    This holds for 32 of the 181 pairs that pass alone, at 16 roots. In 29 of the 32 the candidate also touches
    a second dead end; in 3 it does not. Without the 32, 32 roots (not 42) have two or more candidates that
    work alone, and 77 reactions (not 104) are freed.

  An exact separation needs one more LP test per working candidate: whether the root can be net produced or
  net consumed once the candidate is in. It was not run, because it changes the approved method.

### 4.4 Recorded decisions, gene hits, pooled roots

- **Recorded decisions** ([recorded_decisions.csv](recorded_decisions.csv)). 20 roots have a reaction under a
  recorded decision: 11 under the vacuolar closures, 6 under the acyl-ACP block, 2 under the nuclear RNR copies
  and 1 under "R1725 stays closed" (2026-10-27; added in the follow-up after the reviewer found it missing).
  They are labelled in the tables and nothing is proposed against them. Three have a passing candidate; only
  vacuolar water has a balance-marked one, and it passes only in the second pass. The list holds the four
  decisions named here, not every closure a build has made on purpose.
- **Gene hits.** 18 placements at 14 roots use a reaction that the gene-based gap-fill tables link to a gene (10
  distinct gene ids). 8 work alone and all 18 work with the other roots opened. The ids are in the `pool_genes`
  column; names and functions were not verified here.
- **Pooled roots** (balance-marked candidates). Glycolate has 7 working candidates in the cytosol and 4 in the
  peroxisome. 4-hydroxy-2-oxoglutarate has 1 or 2 in each of its three compartments. Dolichyl diphosphate has
  none.

## 5. Reading for a GNN (inferred)

1. **The pre-set rule.** 57 roots is above the cut of 30, so single-reaction ranking has a test bed. The count
   is 56 with the formula check and 51 without the loop suspects of the reviewer's proxy; the cut holds in each
   case.
2. **Its size.** It is small: 57 roots with 2 to 12 working candidates each, 220 root–reaction pairs. A list of
   that size can be reviewed by hand, so progress on it does not wait for a model.
3. **The LP is a weak filter.** It keeps 243 of 250 pairs in the second pass, and it accepts loops, including one
   that makes ATP from nothing. In a design where a network ranks and the LP checks, the LP check would remove
   little, and evidence has to carry the choice. A gene hit exists at only 14 roots.
4. **Most dead ends are out of reach of a single-reaction ranker.** 283 of the 364 roots have no balance-marked
   candidate (the first six rows of the table in 4.1), and 252 blocked reactions need more than one fix.
5. **Where a learned score could matter more** is the choice among sets of reactions, for example for the 58
   roots whose candidates each need a compound the compartment lacks. That was not tested here.
6. **Not tested:** whether a GNN would rank the 220 pairs better than a simple evidence rule.

## 6. How the code works, and what kind of statement each result is

**Check 1** uses the October census scripts byte for byte (hashes in the artifacts folder, `logs/`).
`gapfind.py` runs FVA with no growth requirement and calls a reaction blocked when both ends of its range are
below cobra's zero cutoff (1e-7 with GLPK). It then scans every metabolite's reactions with their bounds to find
the roots, removes the roots' reactions, and repeats; what falls out this way is "explained structurally".
`fixtests.py` adds one elementary change at a time (a free sink on a root, a link to the same compound in another
compartment, an opened closure, an exchange, a cofactor drain) and re-runs FVA on the affected reactions.
[census_tables.py](scripts/census_tables.py) is new. It checks the gates, calls the October function
`assign_primary_cause` unchanged, adds the group column and the comparison with October, and writes the tables.

**Check 2** is [candidate_tests.py](scripts/candidate_tests.py), new.
- It lists the candidates from MetaNetX with the rule in TASK.md, and stops unless the list equals the count made
  at plan time.
- For each candidate it adds one reversible reaction to the model in memory and solves twice: the smallest and
  the largest flux of that reaction. Two solves are enough:
  - If both are zero, the reaction must stay at zero, so every other flux range is what it was without it.
  - If not, some reaction at the root is now active. A root is a gap metabolite, so all of its reactions were
    blocked, and the candidate's non-zero coefficient on the root forces another of them to run.
- For candidates that work alone it re-runs FVA on the blocked reactions of the root's module and counts the
  ones that became active.
- Second pass: the same two solves with a free sink on each of the other 363 census roots. A pass there means
  "works once the other dead ends are fixed". It is an upper bound, because a free sink is more generous than
  any real fix.

| Kind of statement | What it covers |
| --- | --- |
| Provable from the model structure | the root classes; the reactions explained structurally; that a root in no open reaction cannot be fixed by one added reaction; that a candidate needing a compound absent from the compartment cannot carry flux alone; the two-solve argument above |
| Observed in this run (LP) | the blocked sets (GLPK and Gurobi), the fix and drain test counts, every candidate pass or fail, the re-counts |
| Inferred | the grouping of causes into the two questions; the reading for a GNN in section 5 |
| Untested | whether any candidate enzyme exists in *Y. lipolytica*; whether a passing direction is thermodynamically possible; fills that need two or more new reactions; real transporter evidence |

**Possible improvements, and how to test each.**
1. *Ids.* Candidates are joined on MetaNetX compound ids, which 256 of the 1,883 species lack. Matching those
   species by InChIKey or formula would enlarge the list. Test: count how many of the 61 id-less roots gain a
   candidate, and re-run check 2 on the additions.
2. *Sets of reactions.* A root whose candidates each need one absent compound can only be fixed by two new
   reactions. A minimal-set search per module (a MILP over the MetaNetX reactions placed in the module's
   compartments) answers that. Test: run it first on the roots that pass only in the second pass.
3. *Direction.* Candidates are tested as reversible. Restricting each to the directions a reference database
   records would remove passes that rely on an impossible direction. Test: re-count the passes with those bounds.
4. *Coupling group.* The census cannot say whether a "no dead end found" reaction lacks a step or only a drain.
   A sink test on every gap metabolite of those modules (as the roots get) would split them. Test: compare with
   the hand-made targeted tests of the October report for UM15, UM20, UM21 and UM23.

## 7. Limits

- No model change and no new reference. Every fix test and candidate test is LP evidence for a later, authorized
  curation line, not a recommendation.
- The candidate source is MetaNetX/MNXref 4.5 only. Its balance mark is MetaNetX's own; candidates were not
  re-balanced with the model's formulas and charges.
- The id join is imperfect. At plan time, 219 of the 828 annotated single-compartment reactions that could be
  compared differed from their MetaNetX equation under a strict comparison (same substrate and product sets, H+
  ignored). The cause was not investigated.
- A reaction "already in the compartment" is recognised by its substrate and product sets. A candidate that
  differs from an existing reaction only in a cofactor or in water counts as new.
- A pass also counts a candidate that only closes a loop with the root's own reaction (section 4.3). The
  two-solve test cannot separate such passes, and the counts of loop suspects given there come from proxies.
- The freed-reaction counts are lower bounds (section 4.2).
- The two new scripts have no test file. The new gates were exercised once by hand on scratch copies (section 8).
  In the test stage, the branches for a non-optimal solve never ran, because every solve ended optimal; a
  non-optimal solve inside a re-count could still be miscounted there.
- The closure test of section 3.3 is GLPK only and was not cross-checked with Gurobi.
- The second pass opens the 364 census roots, not the 7 pooled-root copies.
- Gene ids in the tables are copied from the gene-based gap-fill tables. Their names and functions were not
  verified here.
- Not covered: single-gene knockouts, other media, the lipid-unlump line (owned by a colleague).

## 8. Independent review

A fresh read-only reviewer checked the two new scripts against TASK.md and against the outputs. It ran no LP.

**Confirmed by its own recount:** the candidate list (391 pairs, 85 identical pairs, 400 placements); all eight
count columns of `root_candidate_counts.csv`; the two-solve argument and the pass-or-fail code (no optimum lies
within 100 times a tolerance); the second pass (working candidates only produce RNP roots and only consume RNC
roots, and the second-pass ranges contain the first-pass ranges); the change classes since October.

| # | Finding | What was done |
| --- | --- | --- |
| 1 | "Works" also counts candidates that only close a loop with the root's own reaction | Reported in sections 1, 4.3 and 5 with two proxy counts. The exact loop test was offered and the user declined it: it would not change the reading (51 or 57, both above 30) |
| 2 | Module labels differ between runs, so they cannot be matched to October's | Stated in 3.1 |
| 3 | The tables stage would still run after a failed consistency check, and does not look for failed fix tests | Fixed in the follow-up: the tables stage stops in both cases. It had not occurred (no problem recorded; all 634 fix-test notes empty) |
| 4 | The hash gate and the LP cap are not in the scripts | Hash gate added to both scripts in the follow-up. The LP cap is still checked by tally |
| 5 | Unresolved solves are kept apart per placement only; a failed re-count could be miscounted | Did not occur: all 3,336 pass-or-fail solves ended optimal. The unfiltered output of the run held 10 lines: 4 result lines, 3 licence notices and 3 cobra compartment notices, no solver warning (`logs/candidate_tests_output_check.txt` in the artifacts folder). Partly fixed in the follow-up, see below |
| 6 | The re-count covers only the root's module | Stated in 4.2 as a lower bound |
| 7 | Three tolerance-fringe reactions are counted without a flag | Stated in 3.5 |
| 8 | LP counts are nominal; one recorded decision is missing from the list | Nominal counts stated in section 2. The decision was added in the follow-up (section 4.4) |
| 9 | The nine reactions that became blocked get the cofactor label without a test of the closures | Tested in the follow-up: reopening R1172 frees all nine (section 3.3) |
| 10 | No focused test for the two scripts | No test file. The new gates were exercised by hand, see below |

**Follow-up (user: "Yes, without it", [TASK.md](TASK.md) addendum).** The fixes for findings 3, 4, 8 and 9 were
made, and the first part of finding 5: a root with an unresolved solve now gets its own column in
[candidate_summary.csv](candidate_summary.csv) instead of being binned as "none works". Not changed: the loop
test of finding 1 (declined), and the handling of a non-optimal solve inside a re-count (the second part of
finding 5, which was not in the list put to the user).

The new gates were exercised once on scratch copies, with no LP:

| Injected condition | Result |
| --- | --- |
| `IYALI26_RESEARCH_ROOT` unset | `census_tables.py` stops, nothing written |
| One wrong hash in a copy of `input_hashes.csv` | `census_tables.py` stops, nothing written |
| A consistency problem written into a copy of `candidate_summary.json` | the tables stage stops, nothing written |
| A "solver failed" note written into a copy of the fix-test output | the tables stage stops, nothing written |
| One pass-or-fail value blanked in a copy of the test table | that root moves to the `unresolved` column; the other three bins and the column sum to 364 |

After the edits both table stages were re-run. Of the 17 report tables, 13 are byte-identical to the first
run's. Three changed as intended (`recorded_decisions.csv`, the label of one root in
`root_candidate_counts.csv`, the new column in `candidate_summary.csv`) and `input_hashes.csv` is new. The first
run's outputs and script versions are kept in the artifacts folder under `run1/`. No second code review was run
on these edits.

In a second pass the reviewer recomputed the numbers of this report and of the `model/STATE.md` entry from the
tables. No number differed. It asked for six wording corrections, all made: the loop proxy flags pairs rather
than proves them; 14 unmarked placements do balance; 3 of the 32 sink-flagged pairs touch no second dead end;
the glycolate counts are balance-marked; 219 differ among 828 comparable reactions; and the fix tests, drain
tests and re-counts are GLPK only. The input hashes (re-hashed after the run) and the software versions are
now in the artifacts folder under `logs/`.
