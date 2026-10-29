# AUDIT: UM2 decision A0 report (2026-10-28)

**How this file was made.**
- An independent read-only auditor subagent was given the report and a brief: re-derive from the model files
  and sources with its own scripts, not the authors' logs.
- It wrote its own scripts and logs in the session scratchpad (`…/scratchpad/audit/`: `static_check.py`,
  `static3.py`, `lp_audit.py`, `lp_audit2.py`, `tsv_scores.py`, `ext_check.py`, `esummary.json`, three
  abstracts).
- The computer went to sleep twice and cut the auditor off **before it wrote this file**. A third resume and
  a keep-awake request were refused by the session's permission check.
- So the main session compiled this file from the auditor's own logs, without re-running anything. The
  verdicts below are what those logs show. The write-up is not independent; the computations are.

## Claims checked

| # | Claim (REPORT.md / TASK.md) | Auditor evidence (its own script, log) | Verdict |
| --- | --- | --- | --- |
| 1 | UM2 and vacuolar reactions and species are identical in `9d655c2b` and `5f54eb22` | `static_check.log`: "UM2 identical: True", "vac\|UM2 identical: True", no reaction or metabolite differs among those compared | supported |
| 2 | Section 3 tables: bounds, GPRs and equations of R795 and R869–R927 | `static_check.log`, row by row | supported (R795–R927 rows read; dipeptide rows not individually listed in the log) |
| 3 | R925 is unbalanced by −3 H and −3 charge; R870 carries asparagine | `static_check.log` (R925 imbalance {charge −3, H −3}; R870 equation with L-asparagine). `static3.log`: in the iYli21 import `05ae6e1`, R925 was a 1:1 H+ antiport (balanced), so the imbalance came later | supported |
| 4 | iYli21 R795 has 1 H+[cy]; iYli21 R1163 / R1363 are (0, 0); R882 is open | `static3.log` | supported |
| 5 | External models: Yeast-GEM and iYali v4.1.2 open (H+ / water diffusion reversible), iYali5 closed except the Ile import | `ext_check.log` (r_1086, r_1203, r_1831, r_2099 and their iYali5 / iYali4.1.2 counterparts) | supported |
| 6 | Yeast-GEM loops with no uptake, also with r_1831 closed | `ext_check.log` (r_1831 0–1000, r_2099 ±500; r_1831 closed: r_2099 0–500, aa transports non-zero) | supported |
| 7 | On `5f54eb22`: WT 0.3116608216; 37/37 blocked as built (Task A and open); pFBA flux 0 | `lp_audit.log` | supported |
| 8 | C′: Task A 25/37 blocked, open 3/37; pFBA 0; no no-uptake flux | `lp_audit.log`; `lp_audit2.log`: freed set = R1163, R795, R869, R870, R878, R879, R882, R883, R887, R888, R897, R898 | supported |
| 9 | "Open all 17" with R1163 as written: no-uptake loop, R1163 0–1000, ten aa transports 0–500 | `lp_audit.log` (both Task A and open models) | supported |
| 10 | B frees 0 on Task A and 13 open; B′ frees 2 on Task A (R795, R1163) and 14 open | `lp_audit2.log` | supported |
| 11 | Whole-model blocked on Task A: 828 as built, 816 under C′, the difference being the 12 reactions above | `lp_audit2.log` (828 / 816 of 2,319) | supported |
| 12 | Knockout scores at 1 / 5 / 10 / 15 % for all six variants, 20 non-optimal (4 mtDNA, so 16 unresolved), 0 failing checks, 0 gains > 1e-6, 0 class changes | `tsv_scores.log`, recomputed from `stage1/ko_um2.tsv`: identical scores in every variant, max \|ΔKO\| ≤ 1.9e-11 | supported |
| 13 | Knockout spot checks | `lp_audit.log`: YALI1E12482g, YALI1F13017g (knock out R795), YALI1F37658g, YALI1E13514g. Optimal, knocked flux 0, growth matches the TSV to 1e-15. YALI1A03322g is infeasible, as in the TSV | supported |
| 14 | 11 V-ATPase genes are false negatives at 15 % | `static3.log`: 11 of the 14 R795 genes, listed | supported |
| 15 | `YALI1B04391g` is not on R2055 / R2056 | `static3.log`: not in the model | supported |
| 16 | The 10 PubMed links point to the cited papers | `esummary.json`: every PMID matches first author, year and topic (Chahomchuen 2009, Ito 2013, Kettner 2003, Kulakovskaya 1993, Preston 1989, Russnak 2001, Tarsio 2011, Tomitori 1999, Tone 2015, Yang 2006). Abstracts saved for Russnak, Kettner, Kulakovskaya | supported (claim-level reading done for those three only) |
| 17 | Wording: nothing inferred is presented as verified; no claim that the model changed | not reached by the auditor | **not checked** |
| 18 | Biology from the evidence subagent: BLAST reciprocal best hits, PO1f screen calls, the remaining abstracts | not in the brief's re-derivation | **not checked** |
| 19 | Escher README numbers (39 reactions; 23 / 6 / 21 nodes) | not reached | **not checked** (the numbers come from `index.json`, written by the build) |

## Corrections needed

None found in the checked claims. Items 17–19 remain unaudited.
