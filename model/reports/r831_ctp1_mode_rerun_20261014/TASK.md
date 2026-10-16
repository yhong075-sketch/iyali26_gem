# R831 (Ctp1-type carrier): is the citrate/malate exchange mode supported for *Y. lipolytica*? (fresh rerun)

Written 2026-10-14, before any run, literature fetch or alignment. Evidence review plus in-memory diagnostics.
No model XML, condition, assumption, curation, media, benchmark or evidence file is written. No commit.

## Authorization (user, 2026-10-14)

"Approved as written: fresh R831 rerun on claude/task-d-78f16d, A–C4, ~22,000 LPs (hard stop 26,000), Evans 1983
excluded, no commit."

- **Fresh run.** An earlier attempt (`.claude/worktrees/jovial-edison-83497a/model/reports/r831_ctp1_mode_20261014/`,
  base `4d97b82`, uncommitted) stopped before REPORT/audit/STATE. Its outputs, texts and sequences are **not used**.
  It is left untouched. Its numbers may be compared at the end only, labelled as a cross-check.
- **Evans et al. 1983** (Eur J Biochem 130:195) is **excluded** as a read source. The earlier attempt saved
  Evans 1983 texts; the user did not supply them. It may be cited only as Kaplan et al. 1995 report it.
- **Not allowed:** model/XML edits; condition, curation, media, benchmark or evidence files outside this folder and
  its research-workspace folder; commit, push, merge, branch switch; deleting others' files; cluster jobs; using the
  user's browser logins for paywalled texts; any change made to improve a match rate.

## Question

R831 "citrate transport": citrate[cy] + malate[mi] ⇌ citrate[mi] + malate[cy], bounds (−1000, 1000), GPR
YALI1F33754g (Ctp1 ortholog; the same gene is on R832 citrate/PEP and R833 citrate/isocitrate).
1. Is the citrate/malate exchange supported for *Y. lipolytica*, in either direction?
2. If R831 is restricted (`off`, `fwd0`, `rev0`), what happens to the W29 13C elastic fit, to complex I KO/WT, and to
   the gene benchmark?

Context: Task D (`8718b75`) found R831 forward closed (= `fwd0`) makes complex I essential at 0 % WT cost; the
screen (`c66f3fe`) showed that call rests on the medium lysine cap. Task C's fit routes net malate import through
R831 reverse.

## Starting state

- Worktree `.claude/worktrees/task-c-9ffd14`, branch `claude/task-d-78f16d` @ `c66f3fe`, clean.
- Model `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml` `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515`.
  No reference model is designated.
- Inputs as pinned by the imported helpers: Task A medium `6eca70fd…`, 13C caps `df1e33a8…`, R825 ASSUMPTION
  `e45cab42…`, PO1f profile `35307853…`, W29 fluxes `844615d0…`.
- Imported helpers (committed on this branch; SHA-checked by the scripts):

| File | SHA-256 |
|---|---|
| `w29_13c_consistency_20261012/common.py` | `487c95d9d32102d78d81815d4abc9e700cab0b09b02e3d3833ef0156649cc839` |
| `w29_13c_consistency_20261012/elastic_fit.py` | `c14d42ddf99f0bdec54df5c294c5d1c44fc6f6c6902f927afbb0a904e1bd9dea` |
| `w29_13c_consistency_20261012/elastic_fit.json` (reference) | `8d7d3bf690d09642159fdfcf9161e3d44caf629a168cfb7b7f54c4ec1229697e` |
| `complexI_transhydrogenation_20261012/common.py` | `bfe0b2ced3d85cab77fae376cc89dcb220651f747ce25aaf5666707ef70b8ab3` |
| `complexI_transhydrogenation_20261012/mechanism.py` | `0a939031e73eee724dc20123ef1ba6d3df1bc0645f7bc0e7039f2b1aea92d527` |
| `complexI_transhydrogenation_20261012/mechanism.json` (reference) | `a26b8f7430abfe118da25873f9d9df81141dcecc7f619d4bcb1aaffe205d6952` |
| `complexI_transhydrogenation_20261012/ratios.py` | `d00481bf77f019f3083fc5fbdfa153e53f4f1a29d2519a35e868419647307471` |
| `complexI_transhydrogenation_20261012/ratios.json` (reference) | `fd7e409dc229432b31b988429fb46968ea5d93af353cfa82ecabc91afa381ce4` |
| `complexI_shuttle_assumption_screen_20261014/run_screen.py` | `c53c8c21cebef614de903c4d4787746c3c26c9ec58ea6c8e110c1c7159f6a173` |
| `$R/artifacts/tasks/r132_mito_adh_20261012/sw_blosum62_gotoh.py` | `04a1f73fc18b96b0d7b3179e62380c438a98cf86574d37b2fb06bec0cb66398c` |
| `$R/.../w29_13c_consistency_20261012/sequences/sc_mcf_reviewed.fasta` | `842e9e594307368c016fc8009173adb27c71539d2aec5ed463d9baf6cb6e5aa1` |
| `$R/.../w29_13c_consistency_20261012/sequences/yl_clib122_mcf.fasta` | `027a613b88610e2e1da4aa5f8ca3a5d165c2f99371310c80802b2aeeae43f7ae` |

## A. Literature (rules fixed before searching)

- Two read-only reviewer subagents, no conversation history; every text they open is saved under
  `$R/artifacts/tasks/r831_ctp1_mode_rerun_20261014/evidence/`. I re-read the key passages before citing.
- **Order and separation:** *Y. lipolytica* first, then *S. cerevisiae* and other yeasts, then mammalian SLC25A1
  (CIC). Non-*Yarrowia* evidence never stands in for *Yarrowia* evidence.
- **Questions:**
  - Does Ctp1 (or the *Yarrowia* ortholog) exchange citrate for malate, in which direction, and how fast relative
    to citrate/citrate and citrate/isocitrate (or PEP)?
  - Do Dic1, Odc1/2, Sfc1 and Oac1 carry malate, in which direction, and what is each one's physiological role?
  - Is there *Yarrowia* evidence on mitochondrial citrate export or on these carriers?
- **Tiers:** verified (re-read by me) / verified (reviewer) / abstract only / unread. A search-engine summary is never
  a read source. Every DOI is opened before it is cited.

## B. Orthology

- Query sequences: the model genes of every mitochondrial malate carrier (R831/R832/R833 YALI1F33754g; R807 and
  R1386 YALI1D03221g; R902 YALI1B04831g) plus R817 YALI1B14382g, R929 YALI1E41148g, R911 YALI1E04823g, fetched from
  NCBI by the model's RefSeq annotation and saved with SHA.
- Forward: each against the 34 reviewed *S. cerevisiae* PF00153 carriers. Reverse: Ctp1, Yhm2, Dic1, Odc1, Odc2, Sfc1,
  Oac1 against the 38 *Y. lipolytica* CLIB122 PF00153 carriers. Task B's SW tool; Task C's FASTA panels.
- Output: reciprocal best hits; paralogs with no model reaction. A family-panel test, not proteome-wide.

## C. Diagnostics (DIAGNOSTIC, in memory)

R831 variants (bounds intersected, never opened):

| Variant | R831 bounds | Meaning |
|---|---|---|
| `asis` | (−1000, 1000) | as supplied |
| `off` | (0, 0) | citrate/malate mode removed; citrate still moves via R817, R832, R833 |
| `fwd0` | (−1000, 0) | citrate in / malate out closed (Task D G3; screen option b) |
| `rev0` | (0, 1000) | citrate out / malate in closed |

**C1. Elastic fit** (`fit_r831.py`): Task C's `fit()`, `MAPPING`, `SETUPS`, media V2 and V1, f = 1, 0.95, 0.90, for
each variant: 96 fits. Task C's `INFO` is extended in memory with R807, R831, R832, R833, R817, R902, R911, R928,
R929, R1386 (face spans). Reproduction: `asis` must reproduce Task C's 24 totals in `elastic_fit.json` within
relative 1e-6 (absolute 1e-9 if 0).

**C2. Complex I** (`complexI_r831.py`):
- Task D's `load_variant` (mechanism.py), configs A and V2, overlay and caps each on or off; WT and KO (R2062) per
  variant. KO/WT and WT cost vs `asis`; strict cutoffs 1/5/10/15 % (15 % primary); no solution = unresolved.
- FVA (overlay + caps on) at the WT and KO optima: inventory carriers + R817, R832, R833, R533, R398, R1220.
- Lysine: for every KO with KO/WT < 0.15 in config A, report whether R1220 is at its cap and the KO/WT with lysine
  uptake opened to 1000 (the screen showed `fwd0` depends on it).
- Floor row: R825 ASSUMPTION + PDH + CS + G6PDH floors (Task D step 5), A and V2, variants `asis`, `off`, `fwd0`,
  `rev0`.
- Reproduction: `asis` and `fwd0` must reproduce `mechanism.json` robustness rows ("none", "G3") within relative 1e-6
  (absolute 1e-9 if 0); the `asis` floor row must reproduce the matching `ratios.json` row.

**C3. What-if carriers** (in `complexI_r831.py`, labelled WHAT-IF, not a proposal): with `off`, A and V2, overlay +
caps on, open in memory W1 R902 reverse (lb −1000), W2 R1386 reverse, W3 both. WT, KO, KO/WT.

**C4. Gene benchmark** (`screen_r831.py`): today's screen runner (`run_screen.py`: Task B `screen()`/`metrics()`
unchanged) in the benchmark context (as loaded, no overlay) for `off` and `rev0`. `fwd0` = screen option (b),
reused from `c66f3fe`, not re-run. Reproduction: the base run must give 82/33/229/697/16 at 15 %. Report every gene
whose ratio changes by > 1e-6 vs base.

## Budget and stop conditions

- About 22,000 LPs (C1 ≈ 17,600; C2 + C3 ≈ 1,000; C4 ≈ 3,700 incl. a base run); hard stop 26,000 across C1–C4.
- Stop if: an input SHA differs; a reproduction check fails; WT is infeasible in a configuration; the budget is
  reached. On a stop, keep outputs written so far and report; no retry, tuning or relaxing.

## Deliverables

TASK.md, REPORT.md (claims tagged verified / inferred / recalled), evidence and carrier tables (CSV), scripts with
JSON/CSV outputs, maps (PNG) for any route judgment, an independent read-only audit, a `model/STATE.md` entry.
Large outputs, texts, sequences, logs: `$R/artifacts/tasks/r831_ctp1_mode_rerun_20261014/`.
