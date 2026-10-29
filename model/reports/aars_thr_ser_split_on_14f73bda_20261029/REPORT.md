# Thr/Ser synthetase split rebuilt on the combined UM7 + R2171 candidate 14f73bda (2026-10-29)

Scope and authorization: [TASK.md](TASK.md) (user, 2026-10-29: "Rebuild on 14f73bda", then "Approve as written").
The curated change, its evidence and its code are those of
[../aars_thr_ser_split_20261029/REPORT.md](../aars_thr_ser_split_20261029/REPORT.md); nothing in them changed here.
Audit: [AUDIT.md](AUDIT.md). Nothing committed or pushed.

## 1. What and why

The split was first built on `965ab947` (R2171 line), giving `6abb6b38`. Main has since gained `14f73bda`, which
combines the UM7 and R2171 lines (PR #33). This run applies the same opt-in step on top of `14f73bda`, so the newest
line carries the split. The target reactions and genes in `14f73bda` are identical to the curated before-state
([before.json](before.json); verified), so the guarded step applies unchanged.

## 2. Builds (verified)

From `platform/` at HEAD `ca76de1`, `/opt/miniconda3/bin/python` (as for `14f73bda`); full commands in
[logs/build_commands.txt](logs/build_commands.txt).

| Build | Seed | Output SHA-256 |
| --- | --- | --- |
| Control: `14f73bda` argv | 0 | `14f73bda4b0772019e82ea18a27e198c5ba175e67fef205a51dbfb090809e757`, byte-identical to main's candidate |
| Candidate: + `--aars-thr-ser-split` | 0 | `d8d1bac78138f59ba6beefcbddae7f93d0a20e245fea6f5dd9c79760ee535751`, [candidates/E5_…_um7_nuclear_diffusion_r2171_1to1_aars_thr_ser_split_20261029.xml](../../candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_um7_nuclear_diffusion_r2171_1to1_aars_thr_ser_split_20261029.xml) |
| Candidate | 1 | `d8d1bac7…`, identical |

Both candidate logs report "aaRS Thr/Ser split: applied". The control and candidate logs carry the same two libSBML
`R304` notes ERROR lines, also seen in this session's first build and in the R2171 task's build logs. The error is
real data loss: `R_R304` has no `<notes>` element in `14f73bda` (checked in the XML this session) or in the new
candidate, so its notes are dropped at write time. Pre-existing; the cause is not investigated here. The control's and
the seed-1 build's `.build.json` are kept in `logs/`; their XMLs were scratch files. The `.build.json` files do not
record `PYTHONHASHSEED`; the seed of each run is recorded only in `logs/build_commands.txt`, regenerated from the
`.build.json` argv after the audit found the first version garbled.

## 3. Whole-model diffs (verified)

`model_diff.py` (unchanged copy, SHA-256 `8a464ff8…`):
- **vs `14f73bda`** ([model_diff_vs_14f73bda.json](model_diff_vs_14f73bda.json)): 9 differences, the same as the first
  build: rules and notes of R759, R760, R723; reaction links of `YALI1A09121g`, `YALI1F14882g`, `YALI1B07957g`.
- **vs `6abb6b38`** ([model_diff_vs_6abb6b38.json](model_diff_vs_6abb6b38.json)): 3 differences, R_ADPRtn and R_PPItn
  added and the "Transport, nucleus" group's members; identical to the combined task's diff of `14f73bda` vs
  `965ab947` (its `logs/model_diff_vs_965ab947.log`, read this session).

So `d8d1bac7` = `14f73bda` + the split = `6abb6b38` + UM7, field for field.

## 4. Knockouts (verified)

Task A medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) + PO1f profile `po1f_sd_leu.json` (`35307853…`),
list `consensus_essential_genes.csv` (`1e887f5a…`); Gurobi 13.0.2 at 1e-9, reset per knockout, targeted check
(`ko_compare_two.py`, unchanged copy, `bb3118a6…`, with `checked_solve.py` from the 2026-10-22 task); 2,164 LPs;
[ko_compare_summary.json](ko_compare_summary.json), [ko_compare.tsv](ko_compare.tsv).

- WT 0.3116608216 on both.
- Growth changes above 1e-6: `YALI1F14882g` and `YALI1F03908g`, 0.31166 → 0; no other gene.
- Class changes: those two genes, FN → TP, at every cutoff; no failing point; no fresh re-solve; the same 20 infeasible
  knockouts on both (4 mtDNA, 16 nuclear).

| Cutoff | 14f73bda TP / FN / FP / TN | d8d1bac7 TP / FN / FP / TN | Recall | Specificity | Precision | MCC |
| --- | --- | --- | --- | --- | --- | --- |
| 1 % | 94 / 221 / 35 / 699 | 96 / 219 / 35 / 699 | 0.298 → 0.305 | 0.952 → 0.952 | 0.729 → 0.733 | 0.350 → 0.356 |
| 5 % | 100 / 215 / 37 / 697 | 102 / 213 / 37 / 697 | 0.317 → 0.324 | 0.950 → 0.950 | 0.730 → 0.734 | 0.363 → 0.370 |
| 10 % | 104 / 211 / 42 / 692 | 106 / 209 / 42 / 692 | 0.330 → 0.337 | 0.943 → 0.943 | 0.712 → 0.716 | 0.361 → 0.368 |
| 15 % | 112 / 203 / 42 / 692 | 114 / 201 / 42 / 692 | 0.356 → 0.362 | 0.943 → 0.943 | 0.727 → 0.731 | 0.386 → 0.393 |

16 unresolved at each cutoff. Positive coverage: 326 of 1,612 listed genes are model genes, 315 scored, 11 among the
unresolved; 1,286 listed positives are out of scope. Negatives are the user-defined class; mtDNA genes are out of
scope; the list was used in development, so these are development-reference scores. The scores equal those of the
first build on `965ab947` (its `ko_compare_summary.json`, this session). `YALI1A09121g` and `YALI1B07957g`, the
mitochondrial-type paralogs, are also listed positives and stay FN (KO/WT 1): R760 carries no flux and the model has
no mitochondrial seryl-tRNA step, as in the first build.

## 5. Acceptance (TASK.md §3)

| # | Criterion | Result |
| ---: | --- | --- |
| 1 | Control reproduces `14f73bda` byte for byte | met |
| 2 | Seeds 0 and 1 identical | met |
| 3 | Diff vs `14f73bda` = the 9 intended differences | met |
| 4 | Diff vs `6abb6b38` = the UM7 additions only | met |
| 5 | Sweep: exactly the two class changes at every cutoff, WT unchanged, no failing point, same unresolved set | met |

## 6. Limits

- The evidence limits of the first build hold: compartment assignments are inferred from homology and UniProt.
- One medium. No pytest run (no code changed; the focused test passed on the merged tree `ca76de1`).
- Not a designated reference. Both `6abb6b38` and `d8d1bac7` exist; which one later work builds on is the user's call.
- Large outputs: none above 1 MB; `ko_compare.tsv` stays in this folder.
