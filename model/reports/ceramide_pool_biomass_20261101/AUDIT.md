# Independent audit: ceramide pool in `biomass_C` (2026-11-01)

A fresh, read-only reviewer (a Claude Code subagent with no conversation history) received [TASK.md](TASK.md), the
diff, the logs and the deliverables. It edited no repo file and ran no LP. Its working files were in the session
scratchpad. Its verdicts, as returned:

| Claim | Verdict | Evidence the auditor saw |
| --- | --- | --- |
| C1. Candidate SHA `2bc263f2…`; seeds 0 and 1 byte-identical; the control without the flag reproduces `c1f70ac4` | supported | Re-hashed the candidate, the seed-1 copy and the control. `build_results.txt` gives exit 0 for all three builds. `build_commands.txt` differs only in the flag, the seed and the output path. The `.build.json` `source_sha256` and `data_sha256` entries equal the current files. |
| C2. Only the 1 metabolite, the 10 `xPOOL_CER_*` reactions, the `biomass_C` coefficient and 7 notes differ | supported | A line diff of the two XMLs shows +89 lines and 0 removed. An independent cobra comparison of every field of every reaction, metabolite and gene, plus the groups, compartments and objective, found exactly those entities. `biomass_C` is the pool's only consumer. |
| C3. The guarded step fails closed, is a no-op on rerun, is context-reversible and writes only its targets | supported, minor gaps | It checks the scope, objective, ceramide identity, exact composition and bounds, and note keys, then requires all parts absent or all present. The 5 tests pass (34 s). Gaps below. |
| C4. The arithmetic for 5.4e-4 holds and the evidence tiers are honest | supported | 0.007 × 0.523 × 0.049 / 331.54 = 5.41e-4. The Rupcic 1996 abstract, which the auditor opened, has 4.9 %, 52.3 % and 0.7 %, and gives 19-phytosphingosine as 90.8 % of LCB. That supports the 331.54 g/mol choice. Ianutsevich is labelled a figure reading; "a minority of SL" and 12.98 mg/g are in the text. |
| C5. pytest shows no new failures and the 5 new tests pass | supported | Nothing appears only after the change. The 42 errors seen only before are cobra long-path errors from the longer extracted path. The 20 errors after are the same 6 dolichol long-path and 14 missing-`model.xml` errors. Passed went from 378 to 425, which is 42 recovered + 5 new. |
| C6. WT is optimal at 0.31205354; ATP with no uptake is 0 on both models | supported | Read from `wt_energy.json`. Statically, no product-only boundary reaction escapes the uptake closure. |
| C7. Exactly the 7 genes change class; 15 %: 136/190/43/695 → 142/184/44/694; same inputs; LPs within the cap | supported | An independent TSV comparison found the same 7, and no other status, scope, certificate or growth change above 1e-3. The input SHAs match the baseline round. 1,123 + 4 = 1,127 of the 1,300 cap. |
| C8. No other file changed; the lipid-unlump line is untouched | supported | `git status` and `git diff --stat`: `cli.py` +2 lines and `main.py` +8 lines, plus the intended new files. |

**Gaps and observations** (none blocks the round; how each is handled is in [REPORT.md](REPORT.md) §7):

1. Note values outside this curation are not pinned, only their keys. This is the same convention as the trehalose
   step.
2. `_check_spec` does not pin the pool's name and charge, or the ceramide pins inside the spec. These come from
   data.
3. Nothing enforces that this stays the last step to edit `biomass_C`; only a comment in `main.py` says so. The
   composition pin makes a later edit fail at build time, but an earlier step added later would also be caught only
   by the composition pin.
4. The pool metabolite and the 10 pseudo-reactions carry no notes or SBO terms of their own; only the `biomass_C`
   notes explain them.
5. Rupcic 1996 measured total long-chain bases in methanol-grown cells, free bases included. The provisional label
   fits.
6. When the audit ran, REPORT.md and the STATE.md entry were not yet written.
