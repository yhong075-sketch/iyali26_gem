# GPI chain links: implementation report

Scope, authorization, budget and acceptance criteria are in [TASK.md](TASK.md). Exact commands are in
[commands.txt](commands.txt).

Evidence tiers used below:
- **verified**: run or checked in this session; the output is cited.
- **provable**: follows from verified facts.
- **inferred**: reasoned, with no direct check.
- **established/uncited**: textbook knowledge; no source was opened this session.

## 1. Summary

- **New candidate (verified):**
  - `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_links_20261021.xml`, SHA-256
    `89559a0f04fa38a58a563d6dc8b2f9279a933b810461e7fa731359ad670db718`.
  - It is `164f235f` (the 2026-10-20 combined candidate) plus the opt-in step `--gpi-chain-links`.
  - It is a candidate, not a reference model, and it is uncommitted.
- **What changed (verified):** the whole-model diff has 8 differences.
  - Five new gene-less reactions.
  - Their membership in three existing groups.

  No existing reaction, metabolite, gene, bound, GPR, annotation or note changed. The model now has
  2,311 reactions (was 2,306), 1,876 metabolites and 1,074 genes.
- **Effect (verified, SD-Leu + PO1f, GLPK):**
  - The GPI chain R371, R113 and R372–R376 was blocked; all 7 reactions and the 5 new ones can now
    carry flux. The FVA maximum is 0.497 at fraction of optimum 0.
  - WT growth is unchanged: 1.430805735352847 before and 1.4308057353528538 after.
  - At the growth optimum the chain carries 0, because the exit is open, not forced.
- **Knockouts (verified):** 1,074 single-gene knockouts on each model.
  - No knockout growth changed by more than 1e-6, except one solver status:
    - **YALI1D22124g:** essential on the baseline, "infeasible" in the new model's sweep.
    - A fresh solve gives optimal growth 0 on both models.
    - The baseline knockout's flux vector is feasible in the new model (residual 7e-15).
  - So there is no biological class change; see §5.
- **Tests (verified):** the only difference from the baseline is the 5 new passing tests (§6).
- **Independent audit:**
  - It supports all 11 claims.
  - Its one code gap (the guard did not watch the chain intermediates) is fixed.
  - The model bytes are unchanged by the fix (§7).

## 2. The change

| ID | Equation (model IDs) | Bounds | Group |
| --- | --- | --- | --- |
| `R_GPIGNT` | UDP-GlcNAc m403[C_cy] + PI m1648[C_em] + 2 H+ m10[C_cy] → UDP m11[C_cy] + GlcNAc-PI m596[C_er] | 0, 1000 | group37, GPI-anchor biosynthesis |
| `R_PEt_em_er` | PE m1633[C_em] ⇌ PE m599[C_er] | −1000, 1000 | group87, Transport, ER to ER membrane |
| `R_DAGt_er_em` | DAG m604[C_er] ⇌ DAG m1649[C_em] | −1000, 1000 | group87 |
| `R_ACt_er_cy` | acetate m597[C_er] ⇌ acetate m202[C_cy] | −1000, 1000 | group79, Transport, ER |
| `R_DM_GPI` | GPI precursor m603[C_er] → | 0, 1000 | group37 |

| Item | Location |
| --- | --- |
| Curated data | [`model/curation/gapfill/gpi_chain_links.json`](../../curation/gapfill/gpi_chain_links.json) (`GPI-CHAIN-LINKS-20261021`). SHA-256 `4eb4a4f1…` after the audit fix; the candidate's build.json records `0ae29768…` |
| Guarded step | `platform/scripts/gem_annotate/gpi_chain_links.py`, `apply_gpi_chain_links` |
| Wiring | `cli.py` flag `--gpi-chain-links`; `main.py` calls the step right after `--r1964-gpr` |
| Focused test | `platform/tests/test_gpi_chain_links.py` |
| READMEs | `model/curation/README.md`; `model/candidates/README.md` |
| Code diff | [logs/code_diff.patch](logs/code_diff.patch) |

**What the step checks before writing.** It fails closed and makes no partial edits.
- The spec's schema, ID, status and exact set of five reaction IDs.
- The name, formula, charge and compartment of all 12 touched species.
- The full reaction list of each C_er dead end (m596, m597, m599, m603, m604) and each chain
  intermediate (m249, m251, m598, m600, m601, m602). A chain species that gained a producer or
  consumer elsewhere stops the step. The intermediates were added after the audit (§7).
- The bounds and exact stoichiometry of R371, R113 and R372–R376.
- That each target group exists exactly once by name.
- That any link already present matches its full definition and group, and that the links are all
  present or all absent.

A rerun returns `already_correct`.

## 3. Mechanism and evidence

**Why the chain was blocked (verified, 164f235f):**
- Nothing made GlcNAc-PI (m596): R371 was its only reaction.
- PE in C_er (m599) was used only by R373 and R376, with no source.
- DAG (m604) and acetate (m597) in C_er had no consumer.
- The end product m603 had no consumer.
- With one missing input, two outputs that cannot leave and a product that cannot leave, every
  reaction in the chain had an FVA maximum of 0.

**Step 1, `R_GPIGNT`:**
- **Equation:** it matches KEGG [R05916](https://www.kegg.jp/entry/R05916) /
  [EC 2.4.1.198](https://www.kegg.jp/entry/2.4.1.198): UDP-GlcNAc + PI = UDP + GlcNAc-PI.
  **Verified**: the entry pages were read on 2026-10-20.
- **Protons:** the 2 H⁺ are not in KEGG, which does not balance protons. They reconcile the
  model's stored species. UDP-GlcNAc m403 is −2 (C17H25N3O17P2) and UDP m11 is 0
  (C9H14N2O12P2), so without them the reaction is off by +2 H and +2 charge. With them, the
  residual is `{}` (verified: [logs/balance.txt](logs/balance.txt)).
- **Topology:** GPI-GnT acts on the cytoplasmic face of the ER, using cytosolic UDP-GlcNAc and PI in
  the ER membrane. **Established/uncited.**
- **Species choice:**
  - m403 and m11 are the cytosolic UDP-GlcNAc and UDP that carry formulas.
  - The duplicate cytosolic UDP-GlcNAc m1964 and the ER UDP-GlcNAc m1965 have no formula.
  - PI m1648[C_em] is made by PI synthase R1839.
- **Direction:** irreversible, as glycosyl transfer from a nucleotide sugar. **Established/uncited.**
  - A reverse direction would change nothing here: GlcNAc-PI has no other source.
- **GPR:** none. The GPI-GnT subunits (yeast GPI1, GPI2, GPI3, GPI15, GPI19, ERI1) have not been
  identified in *Y. lipolytica* in this model. The reaction is an orphan.

**The three links:**
- **C_em is the ER membrane.** The SBML names the compartment "endosomal membrane", but its
  reactions are named as ER-membrane reactions, for example:
  - R1706 "phosphatidylethanolamine transport, ER membrane-cytoplasm";
  - R1447 "fatty-acid--CoA ligase (hexadecanoate), ER membrane".

  **Verified** in 164f235f.
- **PE and DAG links:** these follow the model's existing gene-less "ER membrane–ER" transports
  R1694–R1697 (group87, verified). They join two pools of the same lipid in the same organelle and
  imply no transporter protein.
- **Acetate link:** this follows the gene-less cytoplasm–ER diffusion reactions (R1158 H⁺, R1358
  H₂O, R1049 CO₂, R1276 O₂; group79, verified). That acetic acid diffuses across membranes is
  **established/uncited**.

**The exit, `R_DM_GPI`:**
- This is a modelling assumption. It stands in for the GPI transamidase, which attaches the
  precursor to proteins; the model has no GPI-anchored proteins.
- It is open (0–1000), not forced, so no GPI amount is required for growth.
- A GPI term in biomass, which would make the pathway required, was not chosen. It needs a measured
  amount.

**Can the links run backwards into a cycle? (provable, and consistent with the FVA minima):**
- Every chain intermediate between m596 and m603 has only chain reactions.
- m603's only consumer is the irreversible demand, so no reverse flux can enter from the top.
- GlcNAc-PI's only producer is the irreversible R_GPIGNT.
- The FVA minimum of all 12 reactions is 0 (verified), so none runs backwards.

## 4. Build and whole-model verification

- **Control (verified):** the combined argv without the new flag rebuilds `164f235f…` byte for
  byte, so the wiring leaves the default unchanged.
- **Seed independence (verified):** with the flag, `PYTHONHASHSEED=0` and `=1` give the same
  `89559a0f…`.
- **Diff (verified, [model_diff.json](model_diff.json)):** 8 differences.
  - 5 reactions added.
  - group37 gains R_DM_GPI and R_GPIGNT; group79 gains R_ACt_er_cy; group87 gains R_DAGt_er_em and
    R_PEt_em_er.
  - No members removed, and nothing else differs. The diff script compares reactions, metabolites,
    genes, groups, compartments, objective, and model notes and annotation.
- **Balance (verified):** R_GPIGNT, R_PEt_em_er, R_DAGt_er_em and R_ACt_er_cy have residual `{}`.
  R_DM_GPI is a boundary reaction.
  - Each transport joins two species with identical formula and charge.
  - Step 1 balances because PI and GlcNAc-PI share the same `P*2` notation, which cancels (see
    limits).
- **Build log warning:** every build logs a libSBML error while setting R304's notes. That
  includes the control that reproduces `164f235f`. The same line is in the 2026-10-20 combined
  build logs. It predates this change.

## 5. Solves (SD-Leu + PO1f, GLPK; [solves.json](solves.json), [ko_compare_summary.json](ko_compare_summary.json))

| Check | Baseline 164f235f | New 89559a0f |
| --- | --- | --- |
| WT growth | 1.430805735352847, optimal | 1.4308057353528538, optimal |
| FVA max, R371/R113/R372–R376 (fraction 0) | 0 for all 7 | 0.497282 for all 7 |
| FVA max, the 5 links | not present | 0.497282; R_PEt_em_er 0.994565, since PE is used twice |
| FVA min, all of the above | 0 | 0 |
| Exit forced to ≥ 1e-5 (pFBA) | not applicable | optimal, growth 1.430777; every chain reaction carries 1e-5; R1839 is 1.17e-4 |
| …with R_GPIGNT, R_PEt_em_er, R_DAGt_er_em or R_ACt_er_cy closed | not applicable | infeasible in each case |

**Knockouts:**
- 1,074 genes per model; the plasmid pseudo-gene `PO1f_plasmid_LEU2` is excluded. 2,150 LPs.
- Classes use KO/WT strictly below the cutoff.
- This is a comparison of two models, not a benchmark score: no experimental labels are used.
- No gene's knockout growth differs by more than 1e-6, except YALI1D22124g.

**YALI1D22124g** (no verified gene name; the sole gene of 24 fatty-acid–CoA ligase reactions in the
ER membrane, lipid particle and cell envelope: R1445–R1456, R1460–R1469, R1870 and R1871. That is a
model assignment, not experimental validation):
- **Sweep:** on the baseline, optimal with growth −2.9e-15 (essential at every cutoff). On the new
  model, GLPK reports `infeasible`, so the call is unresolved.
- **Diagnosis** ([ko_d22124g_diagnosis.json](ko_d22124g_diagnosis.json), verified):
  - A fresh solve of the same knockout is optimal at growth −2.4e-17 on the new model and 0 on the
    baseline. So the gene is essential on both models.
  - Closing the five links gives optimal 0 too.
  - The baseline knockout's flux vector, with 0 on the five links, satisfies every mass-balance
    row and bound of the new model: the largest row residual is 7.2e-15 (m141[C_cy]) and the
    largest bound violation 3.7e-18.
  - Adding reactions whose bounds include 0 cannot remove a feasible point, so the new LP is
    feasible (**provable**).
  - The `infeasible` status therefore comes from GLPK's state partway through the sequential
    knockout loop, not from the change.
- **Recording:** the raw sweep result is kept unchanged in [ko_compare.tsv](ko_compare.tsv). Per
  `model/expected/benchmark.md`, a failed solve is not biological evidence. The sweep's essential
  counts at 1/5/10/15 % are 91/99/108/134 on the baseline and 90/98/107/133 on the new model. The
  difference is this one unresolved call.
- **Effect on the benchmark tool:** if it runs knockouts the same sequential way, it may show this
  gene as unresolved on the new model. **Inferred**, not run.

TASK.md listed "any knockout class change" as a stop condition. The diagnosis shows no class change
in a fresh solve, so the work continued. This decision is recorded here rather than made silently.

## 6. Tests

- **Focused test (verified):** `tests/test_gpi_chain_links.py` passes (5 tests, 23 subtests after
  the audit fix; 22 before). It covers:
  - 10 stale-input cases, each failing closed with the model unchanged;
  - 3 scope rejections;
  - exact additions, balance and group membership;
  - idempotence;
  - the SBML roundtrip;
  - a producer and a consumer for each former dead end;
  - the opt-in flag.

  It runs no LPs.
- **Full suite, per test (verified):** run with `--continue-on-collection-errors` before and after
  the change. The only difference is the 5 new tests passing (212 → 217 passed). The 11 failures,
  17 errors and 31 skips are identical by test ID. The 17 errors are 14 pre-existing errors plus 3
  gem_review lit files that cannot import `pydantic_ai` in this venv.
  - [logs/pytest_before_results.txt](logs/pytest_before_results.txt)
  - [logs/pytest_after_results.txt](logs/pytest_after_results.txt)
  - [logs/pytest_summaries.txt](logs/pytest_summaries.txt)
- **Note-text fix during development:** the first roundtrip test failed. Apostrophes in two note
  strings came back from SBML as `&apos;`. They were removed from the data file before any build.

## 7. Independent audit

The full audit is in [AUDIT.md](AUDIT.md). A fresh read-only reviewer solved 22 GLPK LPs and edited
nothing.
- **C1–C4 and C6–C11: supported.**
- **C5 (guard): supported for the listed cases, with one gap.**

| Finding | Response |
| --- | --- |
| 2. The guard did not watch the chain intermediates m249, m251, m598, m600, m601 and m602 | **Fixed.** The data file now records `reactions_before` for all six; each currently has exactly two chain reactions. The test gained a probe: a sink on m602 is rejected. The guard change leaves the model unchanged, so the candidate stays `89559a0f…`. The verification rebuild and full test run are in §7.1. Data SHA-256: first build `0ae29768…`, after the fix `4eb4a4f1…`. |
| 1. The knockout summary file still records an unresolved call | Kept as raw output, as disclosed in §5. The new-model essential counts in `ko_compare_summary.json` must not be quoted as biology. In a fresh solve, YALI1D22124g is essential on both models. |
| 6. "KEGG does not balance protons" in the R_GPIGNT note is imprecise | **Corrected here, not in the model.** KEGG's equation is balanced for the neutral species. The 2 H⁺ are needed only because the model stores UDP-GlcNAc (m403) at −2 and UDP (m11) at 0, a convention R778 and R211 share. Changing the note would change the model's bytes and require a second candidate. The note text is to be corrected in the next GPI change (the co-products), which rebuilds anyway. |
| 3. The no-cycle argument rests on R_DM_GPI's lower bound of 0 and on R_GPIGNT being irreversible | Accepted as a limit (§8). A workflow that opens every boundary reaction for uptake would make R_DM_GPI a source of m603. Under the project medium its bounds stay (0, 1000). |
| 5. R_GPIGNT spans three compartments and is gene-less | Accepted as a limit (§8). No knockout effect under SD-Leu + PO1f; other media were not tested. |
| 9. R1000 (acetate, cytoplasm–nucleus, gene-less) is a closer precedent for R_ACt_er_cy | Agreed; recorded here. The model note cites the ER diffusion reactions. |
| 4, 7, 8, 10 | Already in §8 or informational. `baseline_model_sha256` in the data file is a record and is not checked. The guard checks the model state itself. |

### 7.1 Verification after the guard fix

All verified:
- **Focused test:** 5 tests and 23 subtests pass; the new subtest is the m602 sink probe. Before the
  fix, the auditor's probe showed this state was accepted.
- **Verification rebuild:** the combined argv plus `--gpi-chain-links` with `PYTHONHASHSEED=0`,
  written to `.scratch/gpi_links_20261021/gpi_links_after_audit_seed0.xml`, gives `89559a0f…`, the
  same bytes as the candidate.
  - Its build.json records the data at `4eb4a4f1…` and the module at `1684c8bb…`; only a code
    comment changed.
  - The candidate's build.json keeps the first build's `0ae29768…` and `271018e7…`.
  - The LP and knockout results in §5 apply unchanged, because the model bytes are identical.
- **Full suite:** the only difference from the baseline run is the 5 new passing tests (212 → 217
  passed). The 11 failures, 17 errors and 31 skips are identical by test ID.
  - [logs/pytest_after_audit_results.txt](logs/pytest_after_audit_results.txt)

## 8. Limits and next steps

- **Co-products (next change, as approved):** R373 (PE → EtNP) omits DAG, and R375 (third mannose)
  omits Dol-P. KEGG gives DAG for every PE-donor EtNP transfer
  ([R05923](https://www.kegg.jp/entry/R05923), [R05924](https://www.kegg.jp/entry/R05924),
  [R08107](https://www.kegg.jp/entry/R08107)) and Dol-P for every Dol-P-Man transfer
  ([R05919](https://www.kegg.jp/entry/R05919), [R05922](https://www.kegg.jp/entry/R05922),
  [R07129](https://www.kegg.jp/entry/R07129)). **Verified**, read on 2026-10-20.
  - Once the chain carries flux, R375 consumes the dolichol carrier. With an open exit that does not
    happen at the optimum.
  - That change would add R373 to m604's reactions, so it must update this curation's `m604`
    guard row.
- **Candidates, not confirmed bugs:**
  - R371 carries an extra H⁺ (+1 H, +1 charge; KEGG [R05917](https://www.kegg.jp/entry/R05917) has
    none).
  - PI and GlcNAc-PI formulas use `…P*2`, which cobra reads as two phosphorus atoms.
  - R113's product m251 has no formula.
  - All 7 GPI reactions are reversible.
- **Modelling limits (audit findings 3 and 5):**
  - The no-cycle argument rests on R_DM_GPI's lower bound of 0 and on R_GPIGNT being
    irreversible. Opening every boundary reaction for uptake would turn R_DM_GPI into a source.
  - R_GPIGNT is a gene-less reaction spanning three compartments (C_cy, C_em, C_er). It implicitly
    moves the PI moiety into C_er, because the model has no compartment for the cytosolic leaflet.
    It had no knockout effect under SD-Leu + PO1f; other media were not tested.
- **R_GPIGNT note text:** the rationale says "KEGG does not balance protons", which is imprecise
  (§7). Correct it in the next GPI rebuild.
- **Not done:**
  - a step-1 GPR (the GPI-GnT subunit search);
  - the GPI7 and fourth-mannose steps;
  - a GPI biomass term;
  - a full blocked-reaction count before and after;
  - a benchmark score.
- **Claude Science's figure:** its growth (0.31240) and R1839 flux (0.00154) do not reproduce on
  `164f235f` under the model's own bounds, SD-Leu or SD-Leu + PO1f. Its structural result does. The
  file or medium it used is unknown.
- **Provenance:** this curation is development work, not independent validation.
