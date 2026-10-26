# Independent audit: GPI co-products (2026-10-21)

A fresh read-only reviewer audited this change. It had no conversation history and was given TASK.md,
the data, code, test, models, logs and result files, plus claims C1–C11. It was told not to treat
REPORT.md as evidence.
- It solved 23 GLPK LPs.
- It made read-only web lookups (KEGG, UniProt, ChEBI).
- It edited nothing in the repository (`git status` was unchanged).

The verdicts and findings below are the reviewer's, condensed; numbers are kept as returned. The
responses are in [REPORT.md §7](REPORT.md#7-independent-audit).

## Verdicts

| # | Verdict | Evidence (reviewer's own checks) |
|---|---|---|
| C1 | Supported | The control and the original both hash to `164f235f…`, and `cmp` finds them identical. |
| C2 | Supported, with a wording caveat | The raw XML diff of links-only vs `89559a0f` is one line: R_GPIGNT's `gpi_chain_links_rationale`. m403 is C17H25N3O17P2 at −2; m11 is C9H14N2O12P2 at 0. KEGG R05916 is written for C00043 and C00015 (neutral). R_GPIGNT is balanced. Caveat: R778 is balanced, but R211 is not (H +5, O +1, charge +3), so "as R778 and R211 do" holds only as "they use the same species", not as a balanced precedent. |
| C3 | Supported | Seeds 0 and 1 hash to `cf0825b8…` (`cmp` identical). Both build.json files have the combined argv plus both flags. PYTHONHASHSEED is recorded only in commands.txt. |
| C4 | Supported (independent) | Raw XML: only two `<notes>` blocks and two product speciesReferences (m604 in R373, m503 in R375). The cobra diff finds 4 differences. Notes equal the 7-key dictionary from the data file. Counts 2311/1876/1074/113 in both. |
| C5 | Supported, with one nuance | All 24 reactions of map00563 were fetched. The five PE-donor reactions (R05923, R05924, R08107, R12351, R12702) all release C00641. R05922 (G00148 → G00149) matches m601 → m602 down to the linkages. Gene identities hold (MCD4_YARLI Q6C0Z3 / yli:2908691 / K05285; GPI10_YARLI Q6CAB8 / yli:2910629 / K05286); both UniProt entries are protein-existence level 3 (homology). m502 − m503 = C6H10O5; R372 and R374 already return m503. Nuance: R12351 acts on a Man2 species, R373 on Man1-GPI. So it is the same enzyme and co-product, not the same step. |
| C6 | Supported for every listed case, with gaps (finding 3) | All listed probes raise and leave the model unchanged; the step is idempotent; `build_reference_chain(gpi_coproducts=True)` raises. |
| C7 | Supported | 1.4308057353528538 vs 1.4308057353528474. |
| C8 | Supported; mechanism confirmed | Chain max 0.49728237112540 → 0.74967868611932; R_DAGt_er_em max 1.499357372238637 (2×); forced cost 2.8322539e-5 → 1.8771011e-5; each link closure infeasible. Decomposition: the Dol-P fix alone gives max 0.6518 and the DAG fix alone 0.5519. About 70% of the saving comes from Dol-P recycling (R2181, de novo Dol-P, carries 1e-5 in the forced baseline and 0 in the new model) and about 30% from DAG return. The cost ratio matches the max ratio (1.5088 vs 1.5076). |
| C9 | Partially supported | YALI1D22124g is the only status or class difference; fresh solves give 0 on both models. **But 9 genes show a systematic +8.477e-7 knockout gain, just under the 1e-6 tolerance:** YALI0F04114g, YALI1E12730g, YALI1E15606g, YALI1E19902g, YALI1E23635g, YALI1F05089g, YALI1M00131g, YALI1M00277g, YALI1M00325g. These were identical in the links task. |
| C10 | Supported, with a caveat | The chain flux is v ≥ 0 by structure; min R373 = 0. Caveat: the open exit plus the co-products now makes the chain a net-beneficial sink in some knockouts (finding 1). |
| C11 | (a) Partially supported, misleading conclusion; (b) partially supported | See findings 2 and 4. |

## Findings (reviewer's ranking)

1. **Medium.** The open exit now acts as a PE disposal route in respiratory knockouts.
   - For YALI1E12730g (sole gene of R304), the new model runs 2.1375e-5 through the chain and
     returns DAG at 4.275e-5, replacing R1773's 4.275e-5.
   - Closing R_DM_GPI restores the baseline growth exactly.
   - Duals: PE[C_em] is in surplus. The chain prices at −0.0397 per unit, which equals the gain.
   - "Chain carries 0 at the optimum" holds for WT, not for every knockout. No class changes.
2. **Medium.** The TASK.md finding (a) on R376 draws a misleading conclusion.
   - YALI1E14925g is GPI13. The NCBI feature table maps YALI1_E14925g to GeneID 2912812 and
     XP_503848.3; the S2 locus map gives YALI0E12089g. KEGG yli:2912812 is K05288.
   - ChEBI 53060's definition describes two EtNP; its KEGG cross-reference G13044 is (EtN)2, the
     product of R08107.
   - So R376 most likely equals R08107. The likelier problems are m603's three-EtNP name and an
     absent GPI7 step (R12702; gene YALI1_D34684g / yli:2911071 not in the model).
   - R08107 lists K05287 and K05288 together, so R376's OR rule is questionable.
3. **Low.** Guard gaps, each accepted and edited:
   - an altered link (closed, retargeted, forced, or with changed stoichiometry), because the check
     was presence-only;
   - a changed m604 annotation;
   - m503 gaining a reaction, R372 losing Dol-P, or m601's name changed;
   - curation notes on the before-state stoichiometry;
   - an emptied `requires_reactions`.

   Latent: only metabolites new in the after-state were added, so a changed shared coefficient would
   report `applied` without reaching the after-state.
4. **Low.** Finding (b) on DAG is partly overstated.
   - C6H8O6 cannot be a DAG (KEGG C00641 is C5H6O5R2). Verified.
   - R1725 and R1726 are imbalanced but closed at (0, 0).
   - R1728 and R1773 contain formula-less species (m1656, m378), so "unbalanced" is not computable
     for them.
5. **Low.** Provenance gaps.
   - `ko_d22124g_fresh.json` records no model SHAs. It reuses the links script, whose
     `new_ko_links_closed` and `baseline_point_in_new_model` fields now mean something different.
   - PYTHONHASHSEED is not recorded in build.json.
   - 58 LPs were used against "about 50".
6. **Info.**
   - `baseline_model_sha256` in the data is not checked; inside the build, the pre-state is the
     links-only `d5a0ce48`.
   - The m604 guards pin the wrong C6H8O6 formula, so fixing DAG will trip both guards. That is
     fail-closed, but it needs coordinated data edits.
   - R371 is imbalanced (pre-existing).
   - pytest per test matches the earlier run plus 4 new passes, which the reviewer reproduced.
