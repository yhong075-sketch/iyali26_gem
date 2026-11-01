# TASK: combine R24 option C with the trehalose line in one candidate

## Authorization

The user, 2026-10-30: "can you combine the branch model edition to make sure our model is newest?". Then, on
2026-10-31, "Approve as written" for this plan:

1. Make a new branch from `origin/main`. Merge the R24 branch and the V-ATPase report branch into it, keeping both
   sides of every conflict.
2. Run the focused tests, then the full pytest, and compare test by test with `origin/main`.
3. Control builds: rebuild `164542fb` from main's arguments, and `f505571d` from the R24 branch's arguments, on the
   merged code; both must match byte for byte.
4. Combined build: `164542fb`'s arguments plus `--r24-homocitrate-aconitase`, under PYTHONHASHSEED 0 and 1.
5. Whole-model diffs. Against `164542fb`, only R24's rule and notes and the gene links may differ. Against
   `f505571d`, only the trehalose term and the `biomass_C` notes may differ.
6. Two-model knockout sweep, about 2,200 LPs: Task A + PO1f, Gurobi at 1e-9, basis reset per knockout, targeted
   check, all four cutoffs.
7. Independent audit, REPORT, STATE entry and candidates README row, then a commit on the new branch. No push.

## Starting state

| Item | Value |
| --- | --- |
| Branch | `claude/combined-r24-trehalose-20261031`, cut from `origin/main` @ `74146c9` (includes PR #39, #40, #41) |
| Merged in | `claude/nice-bouman-b161f4` @ `2ea1e31` (R24 option C, merge `0b1a5f8`); `claude/essential-gene-prediction-e76be3` @ `40212be` (V-ATPase report, merge `59e1547`) |
| Main's newest candidate | `…aco_split_treh_biomass_v2_20261030.xml`, `164542fb6056f01e257ecfaed1b456961e8a07e46e9e4e1942f79489d5014beb` = `e1e7c868` + `--trehalose-biomass-removal` |
| R24 candidate | `…aco_split_r24_aco_20261030.xml`, `f505571d3da6d76920f242e334aa74521d0f754c452b65eb12791e89c079e3c6` = `e1e7c868` + `--r24-homocitrate-aconitase` |
| Common parent | `e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897` |
| Reference model | none designated |
| Scoring rule | `model/expected/benchmark.md` on `74146c9`: certified-infeasible knockouts score as growth 0 (PR #41); mtDNA genes and untargeted COX9 are out of scope. The earlier-rule counts are reported alongside |

## Expected effect

The two steps touch disjoint parts of the model: R24's GPR and notes, and the `biomass_C` trehalose coefficient and
notes. So the combined model should equal `164542fb` plus R24's change. R24 alone changed no knockout class on
`e1e7c868`, so the 15 % scores are expected to equal `164542fb`'s under the same rule. This is inferred; the sweep
tests it.

## Stop condition

Commit on the branch after the audit. Stop and report on any unexpected diff, a failed control build, a new test
failure, or any class change that the two parents do not explain.
