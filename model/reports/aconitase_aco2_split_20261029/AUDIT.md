# Independent audit: aconitase ACO2 split (2026-10-29)

Two fresh read-only reviewers. They were given TASK.md, REPORT.md, the diff, the logs and the deliverables, not the
author's reasoning. Below is a summary of their returned reports; the verdicts are theirs.

## 1. Claim audit (general-purpose reviewer)

**Overall:** no unsupported claim about the model change, the scores or the evidence.

| Claim | Verdict | What the auditor checked |
| --- | --- | --- |
| Baseline and candidate SHA-256 | supported | re-hashed |
| Control reproduces `059316ef` | supported, with a caveat | `control_059316ef.build.json`; argv matches. Caveat: built before `cli.py` / `main.py` were edited, so the flag-off path of the final code was checked only statically |
| Seeds 0 and 1 identical | supported | both `.build.json` give `e1e7c868…`; current step, data, `cli.py`, `main.py` hashes match the candidate's `source_sha256` / `data_sha256` |
| R304 libSBML error pre-existing | supported | present in the control log |
| 7 differences, raw XML 33 lines | supported | independent `diff`: 12 − / 21 +, byte-identical to `logs/raw_xml_diff.txt`; three hunks R216, R238, R1556; no stoichiometry, bound, species, formula or charge line |
| `YALI1E17935g` kept, no reaction | supported | 1,081 gene products in both; rule references 3 → 0 |
| Score table and MCC, 4 cutoffs | supported | recomputed from `ko_compare.tsv` |
| Only `YALI1D11984g` changes class (5, 10, 15 %, not 1 %) | supported | TSV; only growth change above 1e-6 |
| WT, failing points, non-optimal set, coverage numbers, input SHAs, `ko_compare_two.py` copy | supported | summary JSON, TSV, re-hashes |
| "Gurobi 13.0.2" | unverified | not logged |
| Gene ID mappings | supported | locus map, model RefSeq / ncbigene / KEGG vs UniProt JSON; also found the model's UniProt IDs deleted but sequence-identical via UniParc |
| PANTHER subfamilies | supported | also fetched P19414 (SF3) and P39533 (SF2) live |
| UniProt "against" facts | supported | entry versions and GO evidence codes |
| Alignment numbers | supported | re-ran `align_aconitases.py`; byte-identical TSV |
| Fazius 2012, Holz 2009 attributions | supported, abstract only | abstracts; "active citrate aconitase in vivo" from Holz is a reasonable inference |
| 8 PubMed hits | supported | live esearch, same 8 PMIDs |
| "Does not catalyse the citrate steps", labelled inferred | supported as labelled | |
| "*S. cerevisiae aco1Δ* needs glutamate" | partially supported | the abstract says "essential … for glutamate synthesis", not auxotrophy |
| Mechanism, first bullet | supported | `escher/flux.json` |
| Residual growth explanation | partially supported | arginine/ornithine route supported; aspartate's role at its bound not shown by the trace |
| pytest before/after | partially supported | the test file was changed after the full after-run (4 cases added after the code review); REPORT §3 subtest count stale; "new files present" in the before-run wording inaccurate |
| Acceptance 6, STATE entry | pending at audit time | |

**Omissions it named:**
- out-of-scope positives not classified by reason;
- mtDNA count not explicit in §5;
- "verified name: none" not in the REPORT body;
- Gurobi version not logged;
- deleted UniProt IDs not mentioned;
- figure PNG at the ~1 MB limit (acknowledged in `escher/README.md`).

## 2. Code review (`.claude/agents/code-reviewer.md`)

**No bug, no rule violation, no trap, no scope problem.** One low-severity test gap: there was no test for the
applied-path gene guard, or for the `_check_spec` checks on `before_gpr` and on `YALI1D11984g`'s `reactions_after`.
The reviewer showed in memory that the current code already fails closed for a relinked or removed `YALI1E17935g`.

## 3. What was done about each point

| Point | Action |
| --- | --- |
| Test gap | 4 cases added (29 subtests); focused test passes; full suite run again on the final file (REPORT §7) |
| Control caveat | control rebuilt with the final code, flag off (REPORT §4) |
| Residual-growth wording | REPORT §6 and `escher/README.md` now name only the arginine/ornithine route as traced; aspartate's limit marked untraced; "set by" softened to "probably", still inferred |
| aco1Δ wording | REPORT §6 now quotes the abstract and says no *Y. lipolytica* null data was looked up |
| Subtest count, before-run wording | corrected (REPORT §3, §7) |
| Gurobi version | checked after the run: `gurobipy` (13, 0, 2) in `/opt/miniconda3/bin/python` |
| Deleted UniProt IDs | confirmed this session (both "not part of a reference proteome"); UniParc sequences hash-identical to Q6C9P6 / Q6C5V2; recorded in REPORT §2; annotation not changed |
| Verified names, mtDNA count | added to REPORT §2 and §5 |
| Out-of-scope positives by reason | not done; stated as a limit, as in earlier rounds |
