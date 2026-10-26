# Independent audit: GPI chain links (2026-10-21)

A fresh read-only reviewer audited this change. It had no conversation history and was given only
TASK.md, the code diff, the data file, the models, the logs and the result files, plus claims
C1–C11. It was told not to treat REPORT.md as evidence.
- It solved 22 GLPK LPs.
- It edited nothing in the repository.
- Its scratch files are outside the repository.

The verdicts and findings below are the reviewer's, kept as returned apart from formatting. The
responses are in [REPORT.md §7](REPORT.md#7-independent-audit).

## Verdicts

| # | Verdict | Evidence (reviewer's own checks) |
|---|---|---|
| C1 | Supported | `control_no_gpi_flag.xml` = `164f235f…a77583a6` = baseline. The control's build.json has no `gpi_chain_links` entry. The step runs only under `if gpi_chain_links:`, and the flag is `store_true`. |
| C2 | Supported | The candidate and the seed-1 rebuild both hash to `89559a0f…db718`. build.json `argv` is the baseline's 14 flags in the same order plus `--gpi-chain-links`. `options.gpi_chain_links=True`; `reaction_selection.gpi_chain_links.status="applied"`. The recorded code and data SHAs equal the files at audit time (gpi_chain_links.py `271018e7…`, JSON `0ae29768…`, cli.py `b31d2f85…`, main.py `88881e81…`). The seed values come only from commands.txt. |
| C3 | Supported | Raw XML diff: 4 hunks, 0 removed lines, 107 added: 5 `<reaction>` blocks and 5 group-member lines in groups 37, 79 and 87. The cobra semantic diff agrees: +5 reactions and nothing else changed. The 12 touched species changed only in which reactions they appear in. Groups gain exactly the stated members, with name, kind and annotation unchanged. Objective, compartments and model-level notes and annotation are identical. Counts go from 2306/1876/1074/113 to 2311/1876/1074/113. |
| C4 | Supported | All 5 reactions match the JSON exactly (name, bounds, empty GPR, stoichiometry, annotation, 6 note keys, group, subsystem). Residual `{}` for the 4 non-boundary reactions. Caveat: cobra reads `P*2` as P2 for m1648 and m596. R_GPIGNT still balances because one `*2` species sits on each side; it also balances if `*2` is read as two R groups. |
| C5 | Supported for the listed cases; one gap | All listed probes raise before any edit and leave the state unchanged; the step is idempotent. **Gap:** a new reaction on a chain intermediate (m249, m251, m598, m600, m601, m602) was accepted. A change to a link's subsystem alone was also accepted (cosmetic). |
| C6 | Supported | Growth 1.430805735352847 vs 1.4308057353528538. All links carry 0 at the new WT optimum. |
| C7 | Supported | Baseline: R371 and R375 have max = min = 0; the other five are 0 by structure. New model: max R_DM_GPI 0.4972823711253934, max R371 0.4972823711253908, max R_PEt_em_er 0.9945647422507815. Minima: R_DM_GPI 6e-18, R376 3.5e-16, R_PEt_em_er 0. The remaining endpoints follow from the linear structure; the full 24-endpoint FVA was not re-run. |
| C8 | Supported | Forced R_DM_GPI ≥ 1e-5: optimal at 1.430777412813926; the chain carries 1e-5 and the PE link 2e-5. Closing each of the 4 links gives infeasible, which also follows from the structure. |
| C9 | Supported | A full scan of ko_compare.tsv finds only YALI1D22124g differing. The largest other growth difference is 3.37e-9 (YALI1D01188g). Fresh-context knockout: 0.0 on the baseline and −2.39e-17 on the new model (also 0.0 when solved right after WT). The same 24 reactions are closed in both. It is provable that the change adds no rows and every new reaction admits zero flux. The sweep's "infeasible" itself was not reproduced. |
| C10 | Supported, with one wording issue | KEGG R05916 (https://rest.kegg.jp/get/R05916, fetched 2026-10-21): "UDP-N-acetyl-D-glucosamine + 1-Phosphatidyl-D-myo-inositol <=> UDP + G00143", EC 2.4.1.198. KEGG C00043 is C17H27N3O17P2, against m403 C17H25N3O17P2 at −2; C00015 is C9H14N2O12P2, which equals m11 at 0. So the 2 H+ come from the model's mixed protonation states, consistent with R778 and R211. C_em: 134 of its 195 reactions are named as ER-membrane reactions; the rest are ER lipid enzymes and pools, none named endosomal. R1694–R1697 are gene-less in group87 (R1697 is irreversible). R1158, R1358, R1049 and R1276 are gene-less, reversible C_cy–C_er reactions in group79. |
| C11 | Supported | Every chain and link flux equals v(R_DM_GPI) ≥ 0. With all uptakes closed, max R_DM_GPI = 0. The net route is a pure sink. |

## Findings (reviewer's ranking)

1. **Medium.** `ko_compare_summary.json` records YALI1D22124g as essential → unresolved, and the
   new-model essential counts are one lower. As the sweep files stand, the stop condition is
   literally triggered; it is resolved only by the re-solve outside the sweep. Nobody should quote
   the new-model counts from that file.
2. **Low–medium.** The guard did not watch the chain intermediates m249, m251, m598 and m600–m602.
   A future GPI7 or fourth-mannose step would still be accepted. Suggested fix: add
   `reactions_before` for these six.
3. **Low.** The no-cycle guarantee rests on R_DM_GPI's lower bound of 0 and on R_GPIGNT being
   irreversible; R371–R376 stay at ±1000. A workflow that opens every boundary reaction for uptake
   would make R_DM_GPI a source of m603. In the SD-Leu context its bounds stayed (0, 1000).
4. **Low, scientific.** The change only connects the chain: no gene is essential in silico because
   of it. The missing R373 and R375 co-products overstate the cost of a forced GPI flux (1e-5 lowers
   growth by 2.8e-5).
5. **Low.** R_GPIGNT spans three compartments and has no gene. That is the kind of construct
   flagged as an escape route in essentiality work. No effect under SD-Leu; other media were not
   tested.
6. **Low.** "KEGG does not balance protons" is imprecise: KEGG's equation is balanced for neutral
   species. The 2 H+ are needed only because of the stored charges.
7. **Low, pre-existing.** The duplicate UDP-GlcNAc pools (m403 and m1964/m1965) are not addressed.
   Choosing m403 is defensible.
8. **Low, pre-existing.** cobra reads `P*2` as P2, and R1839 is mass-imbalanced.
9. **Informational.** R1000 (acetate between cytoplasm and nucleus, gene-less) is a closer
   precedent for R_ACt_er_cy.
10. **Informational.** The code style matches `r1964_gpr.py`. The scope gate checks the link IDs
    and takes their definitions from the JSON. `baseline_model_sha256` is never checked.
11. **Not checked by the reviewer:** the full pytest suite, other conditions, and a replay of the
    sweep.
