# AUDIT: BiGG-crosswalk annotations and the mcnaugh_pass-ID export (2026-11-03)

There were two independent, read-only review rounds by a fresh reviewer, who followed `.claude/agents/code-reviewer.md`
and was given only TASK.md, the diff, the logs and the deliverables. It wrote nothing to the repository, and its
scratch work stayed in the session scratchpad. It ran 3 WT solves in all and no knockout sweeps.

## Round 1 (candidate `19cf7cef`, 899 rows)

| # | Claim | Verdict |
| --- | --- | --- |
| 1 | The baseline rebuilds byte-identically at HEAD | supported (it checked the build artifact and its record; it did not rebuild) |
| 2 | The decisions follow from the adjudication by the TASK rules | partially supported (F3, F4) |
| 3 | The builder step fails closed, is idempotent and is opt-in | supported, with a latent defect (F5) |
| 4 | The candidate differs from the baseline only in the 899 targets | supported (own cobra diff plus the raw XML diff) |
| 5 | WT growth is unchanged | supported (it solved the candidate and the export) |
| 6 | The export changes only IDs | supported for the mechanics; content concerns F1 and F2 |
| 7 | No other files changed; the colleague's handoff is untouched | supported |

## Findings and what was done

| # | Finding | Severity | Disposition (round 2 status) |
| --- | --- | --- | --- |
| F1 | The export uses mcnaugh_pass IDs, not BiGG IDs. 914 renamed entities have no BiGG value in the crosswalk, and 18 have an ID that differs from it. | medium | Kept by user decision and described correctly in REPORT.md §5. The tool docstring was corrected after the export (partly addressed). |
| F2 | 16 entities whose BiGG row was not written were still renamed to Andrew's ID | medium | Fixed: the skip list is now all 31 such entities (addressed) |
| F3 | 3 prenyl-diphosphate rows written while twin entities kept ours | medium | Not a defect. The pairs are all-trans vs cis compounds, and the reviewer's own check of the adjudication evidence agrees. TASK.md Amendment 2 states it (closed). |
| F4 | 3 NAD ChEBI rows admitted by a first-token match | low-medium | Fixed: exact-equality rule, rows held (addressed) |
| F5 | A row with before equal to after would break idempotence; the error named the wrong rows | low (latent) | Fixed, with a test (addressed) |
| F6 | The export did not check metabolite identity | low-medium (latent) | Fixed: name, formula and charge are checked, with tests. Residual: metabolites without a link are not checked (addressed; residual in REPORT.md §6). |
| F7 | An unverified export could be left on disk; side files were overwritten | low | Fixed: `.unverified.xml` until the reload check passes; all outputs are refused if they exist (addressed) |
| F8 | Provenance gaps: tool SHA, SHAs of the growth inputs, a stale SHA256SUMS | low | Tool SHA in the sidecar; input SHAs in `verification_inputs_sha256.txt`; SHA256SUMS regenerated at the end (partly addressed at round 2) |
| F9 | Over the solve budget; model-bound solves were uninformative | low | Recorded in REPORT.md §7; the solves were removed from the diff script (addressed) |
| F10 | Test gaps | low | New tests for the no-op row, culprit naming, re-identified metabolites, refused outputs and policy skip reasons. Still untested: a manifest SHA mismatch and the export's failed-reload path (partly addressed). |
| F11 | `proton_translocation_differs` is in the policy but no crosswalk row uses it | trivial | Kept for parity with Andrew's renaming rule (not addressed, harmless) |

## Round 2 (candidate v2 `9158e9f2`, 896 rows)

The reviewer:
- **Generators:** reran both and reproduced the decisions CSV, manifest, policy and `before.json` byte for byte;
- **Decisions:** confirmed that the only content change from v1 is the 3 NAD rows;
- **v2 diff:** found exactly 896 annotation differences from the baseline, and the 3 NAD lines as the only
  difference from `19cf7cef`;
- **Builder guard:** confirmed it reports `already_correct` on v2;
- **Export:** confirmed it changes IDs only, with the counts in REPORT.md §4;
- **Growth:** re-solved the export under v2's mapped bounds (0.31205354361020554, optimal);
- **Handoff:** confirmed the colleague's handoff is unchanged.

It found nothing the fixes broke. Three new low-severity record-keeping issues:

| # | Finding | Disposition |
| --- | --- | --- |
| N1 | The report folder mixed v1 and v2 files, and the TASK.md body kept v1 numbers | The v2 files replaced the v1 copies, which stay in `artifacts/verification/`. Amendment 2 and REPORT.md carry the final numbers. |
| N2 | Decision IDs were renumbered between versions | v1→v2 map in `artifacts/superseded_19cf7cef/decision_id_map_v1_to_v2.csv` |
| N3 | Candidate `19cf7cef` sits in `model/candidates/` with no superseded marker | Marked superseded in REPORT.md §1 and the STATE.md entry |

The reviewer could not check whether Andrew's IDs exist in BiGG, because it had no network access.
