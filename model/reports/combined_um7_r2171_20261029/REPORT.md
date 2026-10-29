# Combined candidate: UM7 nuclear diffusions + R2171 1:1 (2026-10-29)

Evidence tags: **verified** = checked this session, output kept here; **inferred** = reasoned.

## 1. Why

The two newest merged candidates are siblings built on `9d655c2b`:
- `5f54eb22…` = `9d655c2b` + `--um7-nuclear-diffusion` (PR #30,
  [report](../um7_nuclear_transports_20261028/REPORT.md));
- `965ab947…` = `9d655c2b` + `--r2171-pi-plc-1to1` (PR #31,
  [report](../r2171_pi_plc_1to1_20261028/REPORT.md)).

No candidate carried both flags, so the lineage had no single head. The UM11 session proposed one build with
both flags as its baseline. **Authorization:** the user, 2026-10-29, in the UM11 session, answered
"Yes, build combined" to that proposal. Scope: one builder run to a new file, whole-model diffs against each
parent, no other change, no commit. Later committed in `423e194` on branch `claude/wizardly-wing-a0d6c0` on the user's instruction ("Commit the working tree
changes", 2026-10-29). Pushed with a pull request to `main` on the user's instruction relayed by the census session
("Yes, tell the UM11 session to commit its work on its branch, push it and open a PR to main", 2026-10-29).

## 2. What was run (verified)

From `platform/`, at HEAD `ff5192d` (= `origin/main` after `git fetch`; the build JSON records only its own
untracked outputs as dirty):

```
PYTHONHASHSEED=0 python -m scripts.gem_annotate --research-root "$IYALI26_RESEARCH_ROOT" --offline --no-solve \
  --coq9-curation qcycle --energy-candidate E5 --coq9-biomass-dilution --coq9-respiratory-package --r132-mito-adh \
  --restore-selection-cleared-gprs --restore-selection-overwritten-gprs --restore-selection-cpa-ura2 \
  --remove-selection-tkl-fragment-gpr --remove-r153-peptidase-gpr --r2041-forward --r2190-ribulose-product \
  --r1964-gpr --gpi-chain-links --gpi-coproducts --tpp-requirement --thiamine-phosphate-chemistry \
  --thiamine-proton-lock --heme-iron-biomass --r2193-iron-gpr --lipoylation-module --acyl-acp-chemistry \
  --cofactor-biomass-terms --r2135-acyl-coa-pool --r1726-cpt-open --dolichol-route --um7-nuclear-diffusion \
  --r2171-pi-plc-1to1 \
  --output-model ../model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_um7_nuclear_diffusion_r2171_1to1_20261029.xml
```

This is the `5f54eb22` argv (from its `.build.json`) plus `--r2171-pi-plc-1to1`. `main.py` applies the UM7 step
and then the R2171 step; each requires only `--dolichol-route`. Interpreter `/opt/miniconda3/bin/python`
(Python 3.13.5). Logs: [logs/build_seed0.log](logs/build_seed0.log), [logs/build_seed1.log](logs/build_seed1.log).
Both carry the same `LibSBML error code -3` line that every build log of the parents carries.

## 3. Results (verified)

- **Output:** `14f73bda4b0772019e82ea18a27e198c5ba175e67fef205a51dbfb090809e757`, 2,319 reactions,
  1,882 species, 1,081 genes. A rebuild with `PYTHONHASHSEED=1` gives the same bytes
  ([logs/sha256.txt](logs/sha256.txt)).
- **Diff vs `5f54eb22`** ([model_diff_vs_5f54eb22.json](model_diff_vs_5f54eb22.json)): 2 entries, R2171
  stoichiometry and notes. They are field-identical to the R2171 report's diff against `9d655c2b`.
- **Diff vs `965ab947`** ([model_diff_vs_965ab947.json](model_diff_vs_965ab947.json)): 3 entries, R_ADPRtn and
  R_PPItn added, and the "Transport, nucleus" group membership. They are field-identical to the UM7 report's diff
  against `9d655c2b`.
- The diff script is [model_diff.py](../um7_nuclear_transports_20261028/scripts/model_diff.py) (unchanged).
  The parents were read from byte-identical short-path copies (SHAs checked), because cobra rejects
  absolute paths of 260 characters or more on macOS.
- **WT** on the Task A medium + PO1f (Gurobi 1e-9, checked solve, run by the UM11 Stage 1 script):
  0.3116608216, the value both parents report.

**What this shows.** The two changes touch disjoint entities and commute. The combined file is the parent
content plus the other parent's change and nothing else (**verified** at field level). No knockout sweep was run
for the combined file itself. Each parent's sweep found no class change against `9d655c2b`. That the union also
has none is **inferred**; the UM11 Stage 1 sweeps run on this file as their baseline.

Not a designated reference.
