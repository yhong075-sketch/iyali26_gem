# Thr/Ser synthetase split rebuilt on the combined UM7 + R2171 candidate 14f73bda (2026-10-29)

## 0. Authorization

**User, 2026-10-29**, answering "Rebuild the Thr/Ser split on the combined UM7 + R2171 candidate 14f73bda?":

> "Rebuild on 14f73bda" (reply button whose text was: "Rebuild the Thr/Ser split on the combined UM7 + R2171
> candidate 14f73bda through curation-change: draft TASK.md first and stop for approval before building.")

**Gate:** approved — user, 2026-10-29: "Approve as written" (reply button whose text was: "Approved: run the rebuild
plan as written in model/reports/aars_thr_ser_split_on_14f73bda_20261029/TASK.md section 4. No commit or push.").

## 1. Starting state

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/essential-gene-prediction-e76be3` @ `ca76de1` (merge of `origin/main` @ `9774d18` into `6d0f0ff`; PR #34 open), clean tree |
| New baseline | `candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_um7_nuclear_diffusion_r2171_1to1_20261029.xml`, SHA-256 `14f73bda4b0772019e82ea18a27e198c5ba175e67fef205a51dbfb090809e757` (PR #33, commit `423e194`) = `9d655c2b` + `--um7-nuclear-diffusion --r2171-pi-plc-1to1`. Not a designated reference. |
| Its build | `/opt/miniconda3/bin/python`, `PYTHONHASHSEED=0`, argv in its `.build.json` (29 opt-in flags) |
| Existing split candidate | `6abb6b38…` = `965ab947` + `--aars-thr-ser-split` (PR #34); it lacks the UM7 step |
| Curated change | unchanged: [`curation/gpr/aars_thr_ser_split.json`](../../curation/gpr/aars_thr_ser_split.json), step `aars_thr_ser_split.py`, flag `--aars-thr-ser-split` (no flag dependency); evidence and decision in [../aars_thr_ser_split_20261029/](../aars_thr_ser_split_20261029/REPORT.md) |
| Before-state | [before.json](before.json): the R759, R760, R723, R673 and four-gene fields in `14f73bda` are **identical** to those in `965ab947` (verified this session) |
| Guard check | the step applied in memory to `14f73bda` returns `applied` with the curated rules, and the cobra context restores the baseline (verified; no file written) |

## 2. Expected effect

The same three rules, notes and gene links as in `6abb6b38`, now on top of the UM7 and R2171 changes. The UM7 step
adds two gene-less diffusions and the R2171 step changes one stoichiometry; neither touches the targets, so the
knockout effect should equal the one on `965ab947`: `YALI1F14882g` and `YALI1F03908g` FN → TP at every cutoff, no other
class change. 15 %: 112 / 203 / 42 / 692 → 114 / 201 / 42 / 692 (the `14f73bda` sweep of the UM11 task gave the same
baseline as `965ab947`).

## 3. Acceptance criteria

1. A control build of the `14f73bda` argv (same interpreter, `PYTHONHASHSEED=0`) at the current HEAD reproduces
   `14f73bda…` byte for byte. This also confirms the merged code reproduces main's newest candidate.
2. The candidate built with `PYTHONHASHSEED=0` and `=1` gives identical bytes.
3. Whole-model diff vs `14f73bda`: exactly the 9 differences of the first build (R759, R760, R723 rules and notes;
   reaction links of `YALI1A09121g`, `YALI1F14882g`, `YALI1B07957g`).
4. Whole-model diff vs `6abb6b38`: exactly the UM7 additions (R_ADPRtn, R_PPItn, the "Transport, nucleus" group), as
   in the combined task's diff of `14f73bda` vs `965ab947`.
5. Knockout sweep `14f73bda` vs the new candidate: exactly the two expected class changes at every cutoff, WT
   unchanged, no failing point, unresolved set unchanged.

## 4. Plan (for approval)

1. No code, data or test change.
2. Build, from `platform/`, with `/opt/miniconda3/bin/python`:
   - control: `PYTHONHASHSEED=0`, the `14f73bda` argv, output to a scratch path;
   - candidate: `PYTHONHASHSEED=0`, the same argv plus `--aars-thr-ser-split`, output
     `candidates/E5_…_dolichol_route_v2_um7_nuclear_diffusion_r2171_1to1_aars_thr_ser_split_20261029.xml`;
   - seed check: `PYTHONHASHSEED=1`, the candidate argv, output to a scratch path.
3. Whole-model diffs (criteria 3 and 4) with the copied `model_diff.py`.
4. Knockout sweep of both models with the copied `ko_compare_two.py` (Task A medium + PO1f, Gurobi 1e-9, reset per
   knockout, targeted check), run with `/opt/miniconda3/bin/python` (the worktree venv has no scipy).
5. Independent read-only audit, REPORT.md, `model/STATE.md` entry, `candidates/README.md` row. Stop before commit.

**Compute budget:** 2,164 LPs for the sweep plus at most 40 for the audit; 30 minutes. Stop conditions: a control that
does not reproduce `14f73bda`, a seed difference, an unexpected diff or a failing sweep point; report, never re-tune.

**Pytest:** not re-run in full, because no code changes; the focused test already passed on the merged tree
(`ca76de1`, 10 passed with the R2171 tests).
