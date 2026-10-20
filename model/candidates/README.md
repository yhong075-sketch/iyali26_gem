# Candidate models

Built model files, each with its `.build.json` provenance where one was recorded.
None is designated as the reference model yet; choose one and record it in
[../README.md](../README.md). File names are kept from when each file was published,
so older reports and commits still identify them.

| File | SHA-256 | Last changed | What it is |
| --- | --- | --- | --- |
| [legacy/model.xml](legacy/model.xml) | `576a284e` | 2026-06-22 `5e2cc49` | Previous canonical `model.xml` of this line. The builder refuses to overwrite it; memote reports on it. |
| [model_metadata_trna.xml](model_metadata_trna.xml) | `d274bad3` | 2026-09-10 `ff36d87` | Reference chain with the user-selected metadata reaction fields and tRNA-coupled biomass. |
| [model_metadata_trna_ntp1_hydrolysis.xml](model_metadata_trna_ntp1_hydrolysis.xml) | `9adfb0f6` | 2026-09-14 `a286d18` | Adds balanced R_NTP1 ATP hydrolysis. |
| [model_metadata_trna_r153_single_gpr.xml](model_metadata_trna_r153_single_gpr.xml) | `a4a3a8c4` | 2026-09-15 `c31921d` | Single-gene R153 assignment (YALI1D17462g). |
| [model_metadata_trna_r153_merged.xml](model_metadata_trna_r153_merged.xml) | `10a3baa1` | 2026-10-07 `80a6e27` | R2176 merged into R153. |
| [model_metadata_trna_r1931_forward.xml](model_metadata_trna_r1931_forward.xml) | `d417f1de` | 2026-10-07 `80a6e27` | R1931 restricted to the forward direction. |
| [model_metadata_trna_r1159_leak.xml](model_metadata_trna_r1159_leak.xml) | `aad70112` | 2026-09-24 `5a3c094` | R1159 conditional Golgi proton leak; the default build matched this file until R1889 (below). |
| [lipid_unlump/](lipid_unlump/) | | 2026-08-19 to 08-27 | Five strict-sn lipid candidates from the lipid-unlump line (owned by a colleague; review only), with provenance JSON files. |
| [E5_coq9_alpha_1e-4_qcycle_pipeline.xml](E5_coq9_alpha_1e-4_qcycle_pipeline.xml) | `d868f83b` | 2026-10-08 | The CoQ9/R305 candidate (`33468b94`) rebuilt by the builder: `--energy-candidate E5 --coq9-curation qcycle --coq9-biomass-dilution`. Differs from `33468b94` only in the R1889 four-subunit GPR (default since 2026-10-05) and two note texts. |
| [E5_coq9_respiratory_package_20261008.xml](E5_coq9_respiratory_package_20261008.xml) | `0ac2addf` | 2026-10-08 | The rebuilt CoQ9/R305 candidate plus the 2026-10-08 CoQ9/respiratory-chain package (`--coq9-respiratory-package`); see [the report](../reports/coq9_respiratory_package_20261008/REPORT.md). |
| [E5_coq9_resp_pkg_r132_20261012.xml](E5_coq9_resp_pkg_r132_20261012.xml) | `a56f3649` | 2026-10-13 | The respiratory package candidate (`0ac2addf`) with R132 (matrix alcohol dehydrogenase) removed by `--r132-mito-adh` (user decision after the predictor screen); see [the report](../reports/r132_mito_adh_20261012/REPORT.md). |
| [E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_20261019.xml](E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_20261019.xml) | `51c6fabd` | 2026-10-19 | `a56f3649` plus `--restore-selection-cleared-gprs --restore-selection-overwritten-gprs --restore-selection-cpa-ura2`, i.e. the 2026-10-19 restore candidate `c98961e9` plus R671 → `YALI1F23441g` and EG-GPR-001 on R159/R190 (user decision 2026-10-19). See [the report](../reports/selection_gpr_fix_20261019/REPORT.md). |
| [E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_r153_gpr_v2_20261019.xml](E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_r153_gpr_v2_20261019.xml) | `d941d2cd` | 2026-10-19 | `51c6fabd` plus `--remove-r153-peptidase-gpr`: R153 gene-less, `YALI1D17462g` (A1 peptidase) kept unlinked (user decision 2026-10-19). Built with `PYTHONHASHSEED=0` (R159 byte order depends on the seed). See [the report](../reports/r153_gpr_revision_20261019/REPORT.md). |
| [E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_r153_gpr_20261019.xml](E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_r153_gpr_20261019.xml) | `3c9e734f` | 2026-10-19 | Superseded first build of the line above; differs only in two R153 note texts, refined after the audit. Kept until the user decides. |

Candidates published inside their task report folders:

| File | SHA-256 | Published | What it is |
| --- | --- | --- | --- |
| [E5.xml](../reports/atp_candidate_repair_20260924/candidates/E5.xml) | `43042d16` | 2026-09-24 `5a3c094` | E5 ATP energy-repair candidate. |
| [E5_coq9_biomass_alpha_1e-4_validated.xml](../reports/coq_biomass_candidate_20261005/E5_coq9_biomass_alpha_1e-4_validated.xml) | `be0be8dd` | 2026-10-05 `fcf42d2` | E5 plus CoQ9 in biomass at the provisional alpha = 1e-4 mmol/gDW. |
| [E5_coq9_alpha_1e-4_R305_qcycle.xml](../reports/coq_r305_candidate_20261005/E5_coq9_alpha_1e-4_R305_qcycle.xml) | `33468b94` | 2026-10-06 `954de2c` | The CoQ9 biomass candidate with the R305 Q-cycle correction. |

The current default build (`python -m scripts.gem_annotate` from `platform/` with
`--offline --no-solve --coq9-curation metadata`) produces SHA-256 `b4ce0974…`. That is
`model_metadata_trna_r1159_leak.xml` plus the 2026-10-05 R1889 four-subunit GPR. It is not
saved here as a file.
