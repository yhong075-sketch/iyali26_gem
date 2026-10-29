# Hidden-source check (acyl, phosphate, inositol) and R2171, Stage 1, on 9d655c2b (2026-10-28)

## 1. Result

- **Model:** `9d655c2b` (`9d655c2bdbf82721a38f148a9a4ee0a385bb509ea7b591111107113c4aa8335c`), the dolichol candidate on
  `origin/main` `0093c89`. Run on a byte-identical copy at a short path (`cmp` identical; cobra rejects the 277-character
  worktree path on macOS, section 6). Read-only: no model, curation, GPR, bound or test file changed.
- **Authorization:** "Approve Stage 1" (user, 2026-10-28), the sentence in TASK.md section 6. Stage 2 is not authorized.
- **One leak found: R2171 is an inositol source.** With INO1 knocked out and no inositol uptake, the model grows at the
  wild-type rate on every medium tested. R2171 also makes organic phosphate and carbon (inositol 1-phosphate) from
  water. Closing R2171 removes the leak in every case traced.
- **No acyl leak** on Task A, Task A + choline / ethanolamine, or glycerol: the fatty-acid synthase knockout has no
  growth and every acyl-CoA demand is 0. The oleate positive control grows, so the test can see growth.
- **No free-phosphate leak:** with every phosphate exchange closed, Pi demand and growth are 0 everywhere.
- **Opening the DAG kinase R1720 (in memory) adds no new leak**, but it lets the R2171 loop run without the Kennedy
  step R1726 (section 3). R2171 should be fixed before R1720 is opened.
- **LPs:** 279 of the 300 authorized (source check 268, diagnostics 5, figure points 6).

## 2. What the code does

- [`scripts/source_check.py`](scripts/source_check.py): for each condition (C1 Task A; C2 + choline and ethanolamine
  uptake 1; C3 glycerol; C4 oleate; C5 open, every exchange (−1000, 1000) minus 15 formula-less exchanges) and each
  source:
  1. WT growth (checked solve: Gurobi 1e-9, reset per solve, bound / residual / knocked-flux checks);
  2. the source removed (acyl: FAS knockout `YALI1B19844g` plus 9 acyl exchanges closed to uptake, except C4;
     phosphate: 11 P-containing exchanges closed; inositol: INO1 knockout `YALI1B05798g` plus R950 and R1263 closed),
     growth checked;
  3. with the source removed, an in-memory demand reaction maximised for each probe species (acyl-CoA pools `m567`,
     `m570`, `m790`, `m1984`, palmitoyl-CoA `m1349`; Pi `m35`; myo-inositol `m814`, Ins1P `m813`);
  4. R2171 FVA at 99.9 % and at growth ≥ 0, sources present;
  5. the same blocks with R1720 opened to (0, 1000) on C1 and C5;
  6. for every non-zero demand (C4 acyl excluded): a pFBA route at 99.9 % of the maximum and the maximum with each of
     R2171, R2165, R1726, R1842, R2168 closed.
- [`scripts/diagnostics.py`](scripts/diagnostics.py): Ins1P and Pi demand with phosphate closed, and C5 growth with the
  formula-less exchanges open or closed.
- [`scripts/figure_fluxes.py`](scripts/figure_fluxes.py) and [`scripts/build_leak_map.py`](scripts/build_leak_map.py):
  the pFBA points and the Escher maps.
- Full output: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/lipid_source_leak_check_20261028/source_check_full.json`
  (196 kB, routes included); compact copy [`results.json`](results.json); [`diagnostics.json`](diagnostics.json).
  Every solve reported here was optimal and passed the check.

## 3. Results (growth in h⁻¹, demands in mmol gDW⁻¹ h⁻¹; values ≤ 1e-7 are numerical zero)

| Condition | WT | Acyl removed: growth / max acyl demand | P removed: growth / max Pi | Inositol removed: growth / max inositol / max Ins1P |
| --- | ---: | --- | --- | --- |
| C1 Task A | 0.3116608 | 0 / 0 | 0 / 0 | **0.3116608 / 511.4 / 1000** |
| C1 + R1720 open | 0.3116608 | 0 / 0 | 0 / 0 | **0.3116608 / 512.4 / 1000** |
| C2 + choline, ethanolamine | 0.3116819 | 0 / 0 | 0 / 0 | **0.3116819 / 511.4 / 1000** |
| C3 glycerol | 0.4762488 | 0 / 0 | 0 / 0 | **0.4762488 / 520.4 / 1000** |
| C4 oleate | 0.3204501 | 0.3204501 / 0.37 (positive control) | 0 / 0 | **0.3204501 / 510.9 / 1000** |
| C5 open | 0 (invalid, see below) | 0 / 497 (invalid, see below) | 0 / 0 | 0 / 990 / 1000 |

1000 is the demand reaction's upper bound, so the Ins1P maximum is unbounded in practice. In every Ins1P maximum,
R2171 = 10.101 = 1000 / 99: each turn releases 100 Ins1P and uses one inositol to remake PI.

**Single closures** (maximum inositol demand, C1; the same pattern in C2–C5):

| Closed | none | R2171 | R2168 | R1726 | R1842 | R2165 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| C1 base | 511.4 | **0** | **0** | 229.1 | 511.4 | 511.4 |
| C1 + R1720 open | 512.4 | **0** | **0** | 512.4 | 512.4 | 512.4 |

**The loop** (pFBA at the inositol maximum, C1; all steps 5.16 except R558 516.1):
PI [em] → R1777 → PI [cy] → **R2171** → DAG + 100 Ins1P; Ins1P → R558 → inositol → R1607 → PI synthase R1839;
the DAG goes R2169 → **R2168** (acyl-CoA from R2135) → TAG → R2170 → R1707 → R1785 → R1771 backwards → DAG [lp] →
R1778 → R1776 → **R1726** → PC → R1806 → R1842 → PA → R1809 → CDP-DAG synthase R196 → PI. With R1726 closed, a weaker
route returns 229. With R1720 open, DAG kinase returns DAG to PA directly and R1726 is no longer needed.

**At normal growth the leak needs no loop** (figure point F2, C1): with INO1 knocked out and no inositol, growth is
0.3116608 (WT 0.3116608); R2171 = 2.36e-7 and IMPase R558 = 2.36e-5, exactly the inositol that PI synthase needs.
The DAG from R2171 goes into TAG through R2168.

**Phosphate** ([diagnostics.json](diagnostics.json)): with all 11 phosphate exchanges closed, Ins1P demand reaches
1000 (R2171 10.10), with water, oxygen and 2.414 glucose as the only uptakes. So R2171 creates organic phosphate and
the six carbons of each Ins1P. Free Pi demand stays 0, with or without R2171. Inferred reason: freeing the phosphate
through IMPase also frees inositol, and cytosolic inositol has no outlet (its only transport, R862, is uptake-only
(0, 1000)).

**R2171 FVA** (sources present): C1–C4 maximum 1.4e-6 to 2.2e-6 at 99.9 % and 1.2e-3 to 1.8e-3 at growth ≥ 0; C5
1.80 at both.

**C5 is not a valid condition for these tests:**
- WT growth is 0 with the 15 formula-less exchanges closed and 122.6 with them open (diagnostics D3 / D4), so some
  biomass component is only reachable through a formula-less exchange (which one was not identified).
- The acyl demands of 497 come from uptake of **hexadecane (R2272, 155) and triolein (R2273, 114)**, real fatty-acid
  sources that the acyl exchange list (built from name keywords) missed. Not a leak.

**Figures** ([figures/](figures/), Escher, values flux × 1000): F1 WT; F2 INO1 knockout without inositol (grows at WT
through R2171); F3 the inositol maximum (the full loop).

## 4. What it means

- **R2171 as written is a mass and charge error.** Its reference reaction ([KEGG R03332](https://www.kegg.jp/entry/R03332),
  verified 2026-10-28) is PI + H2O → DAG + Ins1P, 1:1:1:1; the model's 100 H2O / 100 Ins1P / 100 H+ appears in no
  source found (TASK.md section 2) and is inherited from iYli21. Each PI hydrolysed releases 99 extra inositol
  phosphates.
- **Biology:** the model predicts that an *ino1* strain grows without inositol. In *S. cerevisiae* ino1 mutants are
  inositol auxotrophs (established / uncited; no source opened in this task). On the Task A medium, which contains
  inositol, the INO1 call is not affected; any inositol-free condition, and any carbon or phosphorus accounting that
  includes inositol phosphates, is.
- **Why earlier checks missed it:** Task A supplies inositol and phosphate, and the no-uptake ATP check cannot see a
  species that is not catabolised to ATP. Only a source-removal test shows it.
- **Order of work:** fix R2171 before opening R1720, because R1720 makes the loop independent of R1726.

## 5. Could the method be better

- **Demand tests miss leaks whose by-product has no sink** (free Pi here). Testing every species that carries the
  element (Ins1P as well as Pi) caught it; a general version would maximise total exported phosphorus or carbon
  with all uptake of that element closed. Testable on the R2145 case on `00e8e397`, which should then show its known leak.
- **Exchange lists built from names are incomplete** (hexadecane, triolein). A formula-based list (any exchange whose
  species has C ≥ 8 and a long alkyl chain, or no formula) is more reliable; check it against these two.
- **An open condition should first be checked for WT growth**; C5 failed that and should have been rebuilt before use.

## 6. Verification level and limits

- **Observed in runs (this task):** every number in section 3, the closures and the routes.
- **Statically verified:** R2171's equation, provenance, annotations and gene; R862's bounds; the exchange inventories.
- **Inferred:** why free Pi is 0 (no inositol outlet); that the 100× came from a lumped 0.01 scheme multiplied by 100.
- **Untested:** other base models; other media than C1–C5; the IPC / MIPC synthases, R2165, R2187 beyond the
  closures above; which formula-less exchange C5 needs.
- **cobra path limit:** cobra 0.30 applies the Windows 260-character path test because `"win" in "darwin"`; the
  candidate's absolute path is 277 characters. The first attempt failed before any solve
  ([logs/source_check_attempt1_path_too_long.log](logs/source_check_attempt1_path_too_long.log)); the run used the
  `cmp`-identical copy, and its SHA-256 is recorded in `results.json`.
- **LP count:** `figure_fluxes.json` first recorded 5 LPs for that script; it ran 6. Corrected by hand and in the script.

## 7. Next gate

Stage 2 for R2171 (TASK.md section 5): A, the 1:1 hydrolase as in KEGG R03332 (charge balanced with 2 H+ under the
model's current charges); B, the cyclic-phosphate lyase that the gene's EC 4.6.1.13 suggests (new metabolite); or C,
close it. It needs the user's approval, and its acceptance test is this check: the inositol-removed growth and the
Ins1P demand must go to 0 on C1–C4, with no new acyl or phosphate leak.
