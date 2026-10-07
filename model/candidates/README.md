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
