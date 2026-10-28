# Independent audit: Kennedy R1726 re-run, e9d35857 → ce9e1112 (2026-10-27)

**Auditor:** a fresh read-only general-purpose agent with no conversation history. It was given the
deliverables and the claims, and was not given the reasoning.

**Method:**
- shasum, plain XML diff;
- KEGG, NCBI Gene, UniProt and PubMed (E-utilities) opened;
- Gurobi 1e-9 LPs on in-memory models (primal simplex, basis reset per solve);
- census CSVs compared;
- test-by-test pytest diff;
- the focused pytest.

It changed no file.

**Overall verdict:** "The change is exactly what was authorized … I independently reproduced the SHAs, the
diff, the solve results, the sweep classes and the census." It found documentation issues only.

| Claim | Verdict | Auditor's evidence (summary) |
| --- | --- | --- |
| K1 SHAs, control, seeds, PYTHONHASHSEED | supported | All build SHAs reproduced; seed on line 1 of each build log; `.build.json` code / data SHAs equal the current files; final builds postdate the last edit. Limit: the seed is written by the operator, not captured by the builder. |
| K2 diff = R1726 bounds + 6 notes | supported | Plain XML diff: one changed line (upper bound parameter) plus six `<p>` notes. |
| K3 gene and literature evidence | supported | KEGG yli:2909135 (K00994), NCBI Gene 2909135 = YALI1_C15420g, UniProt Q6CCA9 = YALI0_C10989g; yli:2912382 (K00993) = YALI1_E31498g; PMID 1847919, 7890740, 7961735 read. P17898 and the Szczepańska passage were not re-opened. |
| K4 WT, FVA | supported | WT 0.311682373343 on both; R1726 0 at the optimum, [0, 0.0217389] at 99.9 %; R1725 [0, 0]. |
| K5 no-uptake LPs | supported | Max ATP and biomass 0, optimal, both models. |
| K6 knockout controls | supported | FAS knockout and PSD double knockout lethal on the candidate; with choline / ethanolamine supplied, the PE-methylation knockouts go from 0 to 0.3117034. |
| K7 sweep | supported | Recomputed from the TSV: no class change at any cutoff; largest growth difference 3.0e-12; spot checks reproduce. |
| K8 census | supported | Open condition 720 → 709, the 11 listed reactions; Task A 852 → 843; nothing newly blocked; only R1726's bound differs between the census inputs. |
| K9 pytest | supported | Only the renamed test differs; focused test 5 passed, 25 subtests. |
| K10 guards and flag | supported | Every listed guard present and tested by message; requirement chain r1726 → r2135 → cofactor → lipoylation; R1726 imbalance recomputed. The review rounds are not in the deliverables. |
| K11 limits and tiers | partially supported | Tiers earned, budget overrun disclosed, no science overclaim; documentation gaps below. |

**Gaps the auditor raised and what was done:**
1. AUDIT.md was missing and STATE.md said the audit was pending: this file, and the STATE line updated.
2. TASK.md status lines were stale (section 0 "in progress"; the section 13 draft "No change made").
   Both now give the final status.
3. "UM12" in the open condition: in the `e9d35857` open-condition census run the module is labelled UM11
   (the numbering shifted). REPORT now says so.
4. "(= WT)" in the choline control meant the candidate's choline-supplied WT 0.3117034, not the Task A WT.
   REPORT now gives both, and the +2.1e-5 from supplied choline.
5. Subtests went 781 → 783: REPORT now states it.
6. The code reviews were not on file: [rerun/logs/code_review.md](rerun/logs/code_review.md) now records
   both rounds.
