# Aconitase rules R1556, R216, R238: ACO2-type YALI1E17935g removed, candidate e1e7c868 (2026-10-29)

**Authorization.** The user, 2026-10-29: "Approve the aconitase plan in TASK.md §7 on baseline 059316ef"
([TASK.md](TASK.md) §0), after authorizing the Biopython install ("you can download it, just remember to add in
gitignore"). No commit or push.

**Outcome.** Candidate
`model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_20261029.xml`
(`e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897`) = `059316ef` + `--aconitase-aco2-split`. Not a
designated reference. `YALI1D11984g` goes FN → TP at 5, 10 and 15 %, **not at 1 %** (KO/WT 0.0244); nothing else
changes class. At 15 %: 124 / 191 / 42 / 692 → 125 / 190 / 42 / 692.

## 1. What changed and why

- R238 (citrate → cis-aconitate, mitochondrion), R216 (cis-aconitate → isocitrate, mitochondrion) and R1556
  (cis-aconitate → isocitrate, cytosol): `YALI1D11984g or YALI1E17935g` → `YALI1D11984g`. R239 (cytosolic citrate →
  cis-aconitate) already read `YALI1D11984g` and is untouched.
- `YALI1E17935g` stays in the model with no reaction. Six notes per reaction (curation id, status, authorization,
  previous rule, evidence, tier).
- Why: an OR member needs evidence that it catalyses the step on its own (AGENTS.md §5). `YALI1E17935g` is the
  ACO2-type paralog, and no evidence was found that it converts citrate or cis-aconitate.

## 2. Evidence (TASK.md §2 and §4 in full)

| Claim | Tier |
| --- | --- |
| `YALI1D11984g` = YALI0D09361g = Q6C9P6 (RefSeq XP_502616.2); `YALI1E17935g` = YALI0E14949g = Q6C5V2 (XP_503960.1) | verified (locus map, model RefSeq, UniProt JSON fetched this session) |
| PANTHER: Q6C9P6 in PTHR43160:SF3 (aconitate hydratase, mitochondrial = ACO1 subfamily); Q6C5V2 in SF2 (homocitrate dehydratase = ACO2 subfamily) | verified (UniProt cross-reference) |
| Pairwise identity: `YALI1D11984g` 77.7 % to Aco1p vs 52.3 % to Aco2p; `YALI1E17935g` 70.3 % to Aco2p vs 54.9 % to Aco1p; paralog background 53.7–54.7 % | verified (run this session, [evidence/pairwise_alignment.tsv](evidence/pairwise_alignment.tsv); homology only) |
| *S. cerevisiae* Aco1p serves the TCA cycle, Aco2p "specifically and exclusively contributes to lysine biosynthesis"; filamentous-fungal Aco2 homologues inactive (Fazius 2012, [PMID 23106124](https://pubmed.ncbi.nlm.nih.gov/23106124/)) | verified, abstract only |
| *Y. lipolytica* ACO1 overexpression shifts citric towards isocitric acid (Holz 2009, [PMID 18850095](https://pubmed.ncbi.nlm.nih.gov/18850095/)) | verified, abstract only; that this ACO1 is YALI0D09361g is inferred |
| `YALI1E17935g` does not catalyse the citrate steps in *Y. lipolytica* | **inferred** (transfer from yeast and fungi); no *Y. lipolytica* data |

Against: UniProt's automatic names call both proteins "Aconitate hydratase, mitochondrial", and Q6C5V2 carries GO
aconitate hydratase activity (IBA). Neither is experimental.

**Genes** (AGENTS.md §4):
- `YALI1D11984g`: verified name none (no reviewed *Y. lipolytica* entry; "ACO1" in Holz 2009 is not tied to a locus
  in the abstract). Function: ACO1-type aconitate hydratase, citrate ⇌ isocitrate. Evidence: homology plus a
  *Y. lipolytica* overexpression phenotype (abstract).
- `YALI1E17935g`: verified name none. Function: ACO2-type aconitase-family protein, homocitrate-specific or
  inactive. Evidence: homology plus *S. cerevisiae* and fungal biochemistry; no *Y. lipolytica* data.
- The model's own UniProt IDs for the two genes, A0A1D8NDY4 and A0A1H6Q7Y5, are now deleted in UniProt ("not part
  of a reference proteome"). Their UniParc sequences, UPI00017F7A6A and UPI00003BC7D7, are identical to Q6C9P6 and
  Q6C5V2 (SHA-256 of the sequence strings compared this session), so the alignment used the same proteins. The
  annotation itself was not changed.
- The ACO1 = SF3 and ACO2 = SF2 equation was checked by the auditor on P19414 and P39533 (live UniProt, verified by
  the audit).

## 3. Code

- Data [`model/curation/gpr/aconitase_aco2_split.json`](../../curation/gpr/aconitase_aco2_split.json) (`bf3423c1…`),
  written by [scripts/write_curation.py](scripts/write_curation.py) from [before.json](before.json).
- Step `platform/scripts/gem_annotate/aconitase_aco2_split.py`, `apply_aconitase_aco2_split`: the pattern of the aaRS
  split step. All checks on the three reactions (name, bounds, stoichiometry, species, annotation, rule, notes) and
  both genes (annotation, reaction set) run before any edit; a partly applied state, a stale row or an extra reaction
  on `YALI1E17935g` fails closed; idempotent; notes restored on context exit.
- Flag `--aconitase-aco2-split` in `cli.py` (no flag dependency); called in `main.py` after the R673 step.
- Test `platform/tests/test_aconitase_aco2_split.py`: 6 tests, 29 subtests (14 stale-input cases, 9 scope cases, 2
  applied-model gene-guard cases, 4 knockout cases; plus roundtrip, context exit and opt-in flag); passed
  ([logs/pytest_focused.txt](logs/pytest_focused.txt)). The 4 cases beyond the first version (`before_gpr` and
  `YALI1D11984g`'s `reactions_after` in the spec; a relinked or removed `YALI1E17935g` on an applied model) were added
  after the code review flagged them as untested branches.
- Code review (fresh read-only reviewer, `.claude/agents/code-reviewer.md`): no bug, no rule violation; one test gap,
  fixed as above.

## 4. Builds and diff (verified)

Commands in [logs/build_commands.txt](logs/build_commands.txt); from `platform/`, `/opt/miniconda3/bin/python`.

| Build | SHA-256 |
| --- | --- |
| Control, `059316ef` argv, `PYTHONHASHSEED=0`, before `cli.py` / `main.py` were edited | `059316ef…` (byte-identical to the baseline) |
| Control again with the final code (flag off; [logs/build_control_final.log](logs/build_control_final.log)) | `059316ef…` (byte-identical) |
| Candidate, + `--aconitase-aco2-split`, seed 0 | `e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897` |
| Candidate, seed 1 | identical |

Every build log carries the pre-existing libSBML error on R304's notes. Whole-model diff vs `059316ef`
([model_diff_vs_059316ef.json](model_diff_vs_059316ef.json)): **7 differences, all intended**: the rule and notes of
R1556, R216, R238, and `YALI1E17935g`'s reaction list (R1556, R216, R238 → none). Raw XML: 33 changed lines. No
stoichiometry, bound, formula or charge changed, so mass and charge balance of the touched reactions is as before.

## 5. Knockouts (verified)

Task A medium (`6eca70fd…`) + PO1f (`35307853…`), list `1e887f5a…`
([logs/sweep_inputs_sha256.txt](logs/sweep_inputs_sha256.txt)); Gurobi 13.0.2 (`gurobipy.gurobi.version()` in
`/opt/miniconda3/bin/python`, checked after the run) at 1e-9, reset per knockout, targeted
check (`ko_compare_two.py`, unchanged copy, `bb3118a6…`); SHA-checked short copies; 2,164 LPs, no fresh re-solve;
[ko_compare_summary.json](ko_compare_summary.json), [ko_compare.tsv](ko_compare.tsv).

- WT 0.3116608216 on both models.
- Growth changes above 1e-6: only `YALI1D11984g`, 0.311661 → 0.0076119 (KO/WT 0.02442).
- `YALI1E17935g`: 0.311661 on both (no reaction after the change); TN on both.
- No failing point; identical non-optimal sets (16 unresolved + 4 mtDNA); max knocked flux 0.

| Cutoff | TP / FN / FP / TN before | after | recall | specificity | precision |
| --- | --- | --- | --- | --- | --- |
| 1 % | 97 / 218 / 35 / 699 | 97 / 218 / 35 / 699 | 0.3079 → 0.3079 | 0.9523 | 0.7348 → 0.7348 |
| 5 % | 103 / 212 / 37 / 697 | 104 / 211 / 37 / 697 | 0.3270 → 0.3302 | 0.9496 | 0.7357 → 0.7376 |
| 10 % | 107 / 208 / 42 / 692 | 108 / 207 / 42 / 692 | 0.3397 → 0.3429 | 0.9428 | 0.7181 → 0.7200 |
| **15 %** | 124 / 191 / 42 / 692 | **125 / 190 / 42 / 692** | 0.3937 → 0.3968 | 0.9428 | 0.7470 → 0.7485 |

MCC at 15 %: 0.4225 → 0.4255. 16 unresolved at each cutoff. Positive coverage: 326 of 1,612 listed genes are model
genes, 315 scored, 11 among the unresolved; 1,286 listed positives are out of scope (not classified by reason in this
round, as in the previous rounds). Negatives are the user-defined class (model genes not on the list); the 16
mitochondrially encoded genes are out of scope and excluded from all four classes; the list was used during
development, so these scores are a development reference, not independent validation.

**The brief's expectation "FN → TP at every cutoff" did not hold at 1 %**, as TASK.md §5 predicted from the review's
0.0244.

## 6. Mechanism and the residual growth ([escher/](escher/README.md))

- Before, the knockout closes only R239; R238 and R216 stay open through `YALI1E17935g` and carry 1.57 in pFBA, so
  growth equals WT. After, all four aconitase steps are closed; citrate synthase R236 drops from 1.57 to 0.0167 and
  both isocitrate dehydrogenases go to 0 (figure,
  [aconitase_ko_yali1d11984g_before_vs_after.png](escher/aconitase_ko_yali1d11984g_before_vs_after.png)).
- **Residual growth 0.0076** (one pFBA optimum, [escher/trace_residual.json](escher/trace_residual.json)): the medium
  has no glutamate. The trace shows arginine-derived ornithine feeding glutamate through ornithine transaminase R610
  (0.0053) and P5C dehydrogenase R9 (0.0037). L-arginine (0.00624) and L-aspartate (0.01307) uptakes sit at their
  bounds; why aspartate is limiting was not traced (the aspartate / glutamate shuttle R161, R162, R823, R1386 is net
  zero). The 0.0244 is therefore probably set by the Task A medium's finite amino-acid bounds, which `benchmark.md`
  calls concentration-scaled assumptions, so this gene's class at 1 % depends on them (**inferred**: no shadow prices
  taken, no bounds varied).
- Biologically, Fazius 2012 states that yeast Aco1p is "essential for the citric acid cycle and, thus, for glutamate
  synthesis" (abstract). A near-zero residual on a glutamate-free medium fits that; whether a *Y. lipolytica* ACO1
  null grows with glutamate was not looked up.

## 7. pytest

Full suite, the workflow-automation venv (Python 3.13.5, cobra 0.30.0, pytest 9.0.2), before vs after
([logs/pytest_before.txt](logs/pytest_before.txt), [logs/pytest_after.txt](logs/pytest_after.txt),
[logs/pytest_compare.txt](logs/pytest_compare.txt)): 336 → 342 passed, the same 19 failed and 20 errors, subtests
931 → 956. Test by test, the only differences are the 6 new tests, all passing; no existing test changed outcome. That
after-run used the first version of the test file (25 subtests). After the code review added 4 cases, the full suite
was run again ([logs/pytest_after2.txt](logs/pytest_after2.txt),
[logs/pytest_compare2.txt](logs/pytest_compare2.txt)): again 336 → 342 passed, the same 19 failed and 20 errors, subtests 931 → 960; only the 6 new tests differ, all passing.
The before run (at `dcefab7`, started before the step and test files existed;
[logs/pytest_before_head.txt](logs/pytest_before_head.txt)) matches the previous round's after-run test by test
([logs/pytest_before_vs_previous.txt](logs/pytest_before_vs_previous.txt)).

## 8. Acceptance (TASK.md §6)

| # | Criterion | Result |
| --- | --- | --- |
| 1 | Control reproduces `059316ef` byte for byte | met (twice: before the edits and with the final code) |
| 2 | Seeds 0 and 1 identical | met |
| 3 | Diff only the three rules, their notes, `YALI1E17935g`'s links | met (7 differences) |
| 4 | WT unchanged; only `YALI1D11984g` changes; unresolved set unchanged; no failing point | met; the change is FN → TP at 5/10/15 %, not 1 % |
| 5 | Only the new tests differ in pytest | met (§7) |
| 6 | Audit | met: [AUDIT.md](AUDIT.md); no unsupported claim; the partly supported points were fixed in this report |

## 9. Limits and what was not tested

- Function of `YALI1E17935g` in *Y. lipolytica* is inferred from yeast and fungal biochemistry (abstracts only); the
  alignment shows identity, not catalysis. No *Y. lipolytica* enzymology, localisation or knockout for it exists in the
  8 PubMed hits screened.
- No active-site / specificity-residue comparison and no gene tree with more fungi (offered as a cheaper next step).
- The cytosolic placement of R239 / R1556 is inherited and not reviewed.
- One medium; the 1 % class depends on the medium's assumed amino-acid bounds (§6).
- The knockout list was used in development; scores are a development reference.
- **Follow-up, not done:** R24 / R395 (homocitrate → homoaconitate → homoisocitrate) carry only the LYS4-like
  `YALI1E03387g`; Fazius 2012 reports that homoaconitases do only the second step and aconitases the first. R24 is also
  named "2-methylcitrate dehydratase". Needs its own review; `YALI1E17935g` was not moved there. R238's note
  `PROTEIN_CLASS: 4.2.1.4` (not 4.2.1.3) is inherited and unchanged.

## 10. Method notes

- The guard compares unescaped note text; note text avoids characters that libSBML writes as entities.
- `before_state.py` selected context reactions by name and gene as well as EC: the aconitase reactions carry
  `ec-code 4.2.1.-`, so an EC 4.2.1.3 filter alone would have missed them.
- The sweep and figure LPs: 2,164 + 18 (the trace ran twice; its first run stopped at the JSON write after its
  solves).
