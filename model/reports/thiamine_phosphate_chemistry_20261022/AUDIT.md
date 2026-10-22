# AUDIT: thiamine phosphate chemistry (THIAMINE-P-CHEMISTRY-20261022)

**Date:** 2026-10-22. **Auditor:** an independent, read-only subagent with no conversation history.
This is its report as returned, lightly reformatted. It audited the first candidate, `e8a54520`; the
responses and the rebuild are in [REPORT.md §7](REPORT.md).

**What it did:**
- It edited no files, ran no builder, made no commits and ran no knockout sweep.
- It ran 2 LPs, a re-run of the WT spot check, with output written to scratch.

**What it read:**
- TASK.md and AGENTS.md;
- the curation JSON (`33d6ca62…`) and the step (`5a08584b…`);
- the cli.py and main.py diffs and the focused test;
- `make_curation_json.py`, `state.py` and `wt_spot_check.py`;
- before.json, after.json and model_diff.json;
- the build logs, sha256.txt, the build.json files and the pytest results;
- `microspecies.py` (closure, audit and proton balancer) and `metabolite_microspecies.csv`;
- the raw `model/source/iyali26.xml` species;
- REPORT.md, which appeared during the audit with §7 still PENDING.

**What it ran:**
- Recomputed SHA-256 and ran `cmp` on all models.
- Compared the two models three ways: a raw text diff of the XML, an element-level XML diff and a
  cobra-level whole-model diff.
- Checked balance model-wide on both models.
- Recomputed the pH 7.3 balance in memory.
- Probed the step's guards in memory.
- Ran `audit_component_migration` and `balance_protons_and_water` on in-memory copies only.
- Re-ran the focused pytest (5 passed, 22 subtests).
- Read-only GETs: KEGG REST compounds and reactions; OLS4 ChEBI 9533, 9532, 9534, 17857, 37575,
  58937, 58938 and 58296.
- Grepped the local MetaNetX 4.5 and ExPASy files.

## Claims

| # | Verdict | Evidence |
|---|---|---|
| C1 | **supported** | Raw XML diff: 12 hunks, all inside the 8 target elements. Element diff: only `fbc:chemicalFormula`, `fbc:charge`, notes and annotation; m643's charge is unchanged (0). Cobra diff: metabolites differ only in the 7 species, reactions only in R1342 (annotation and notes). Genes, groups, compartments, objective and model notes are equal; no entity added or removed. R1342's annotation diff is exactly `ec-code: ['2.7.4.16','2.7.4.9'] → '2.7.4.16'`, and its TPP notes are kept. Applying the step in memory to the baseline equals the built file. |
| C2 | **supported** | `cmp`: the control equals the baseline (`8ae6c93b…`), and seed1.xml equals the new model (`e8a54520…`). The build.json source and data hashes match the current files. **Limit:** PYTHONHASHSEED 0/1 is recorded only in the log header lines written by the implementer's wrapper, not in build.json. |
| C3 | **supported** | KEGG formulas: C01081 C12H18N4O4PS, C00068 C12H19N4O7P2S, C03028 C12H20N4O10P3S, C04327 C6H10NO4PS; each KEGG DBLINK points to the ChEBI ID used. OLS4 charges: 9533/9532/9534 = +1, 17857 = 0, with matching formulas. MetaNetX: KEGG, BiGG and ChEBI IDs all map to MNXM662/256/1341/960. **Limit:** the BiGG IDs were checked only through MetaNetX cross-references. |
| C4 | **supported** | `MNXM1137754  TMP  envipathM:…  C6H4Cl3NO  0`. The auditor decoded the InChI as 3,5,6-trichloro-2-methoxypyridine. `MNXM148527  Ttp  biggM:ttp` has no formula or charge. |
| C5 | **supported** | `MNXR104907 … 2.7.4.16;2.7.4.9`. All its cross-references are thiamine-phosphate kinase. ExPASy: 2.7.4.9 is dTMP kinase, 2.7.4.16 thiamine-phosphate kinase. The R1342 change is ec-code only. |
| C6 | **partially supported** | **Confirmed:** before, only R933 is balanced (trivially). After: R1337, R1342, R751, R753, R754 and R933 are balanced; R752, R408 and R750 are each `{charge: 1, H: 1}`; R1343 is `{charge: 2, H: 1}` and R1344 `{H: 1}`. A model-wide scan finds exactly 10 changed residuals, all on these species. **Problem 1:** R752's charge residual went from 0 to +1, so it is "not worse" only by status or total residual, not per component. **Problem 2:** the attribution of R1343/R1344 to m889/m38 is incomplete. In every consistent form of m889 and m38 both stay 2 H short, in what looks like a redox-type discrepancy in the lumped equations. HET-P is not the cause. |
| C7 | **supported** | With partners at MetaNetX 4.5 pH 7.3 forms (ATP MNXM3 −4, ADP MNXM40333 −3, AMP MNXM728294 −2, Pi MNXM9 −2, PPi MNXM11 −3, HMP-PP MNXM1101289 −3) and the thiamine species at MNXM662 −1, MNXM256 −2, MNXM1341 −3 and MNXM960 −2, all nine thiamine reactions balance. With only the thiamine species at pH 7.3, just R754 and R933 balance. |
| C8 | **partially supported** | **Holds:** model staleness fails closed with no edits; the step requires the TPP term and the flag; it touches only the targets; it is idempotent, including after a roundtrip. **Does not hold:** the after-values are not pinned in code. Tampered specs were accepted: a dTMP formula for TMP, a pH 7.3 form for HET-P, m971 relabelled to another molecule with charge −2, `charge: true`, `schema_version: 1.0`. Partner chemistry is not pinned either: with ATP at pH 7.3 the step applies and leaves R1342 at `{charge: 4, H: 4}`. |
| C9 | **supported** | Re-run: both models optimal at 0.31242674050092606. |
| C10 | **supported** | After: "11 failed, 232 passed, 31 skipped, 17 errors, 584 subtests"; the per-test diff adds exactly the 5 new PASSED lines. The before file is byte-identical to the TPP task's final run. Its code identity rests on the implementer's mtime statement. |
| C11 | **supported** | The KEGG equations R00619, R00616, R00617, R02135 and R04448 match the model's equations and annotations. The source model has separate dTMP/dTDP/dTTP species (m505, m506, m842). The exchanges R1341 and R1339 are named "thiamine(1+) monophosphate exchange" and "thiamine(1+) diphosphate(1-) exchange". **Caveats:** TTP's identity rests on R753 alone. R753's EC 2.7.6.2 comes from MetaNetX MNXR104888's union "2.7.4.15;2.7.6.2". |

## (a) Scientific overstatements

1. **R1343/R1344 attribution.** These reactions stay 2 H short in every consistent form of m889 and
   m38. m889 ("ribose-5-phosphate") is annotated KEGG C01068 and MNXM1092921, which are **D-ribitol
   5-phosphate**. Whether its identity or its charge is wrong is not determined.
2. **HMP-PP has no microspecies row.** HMP-PP (m886) is not in `metabolite_microspecies.csv`, and the
   model gives it no MetaNetX annotation.
3. **"None got worse":** R752's charge residual went from 0 to +1.
4. **"(Rhea)" for R750:** R750 has no Rhea or KEGG cross-reference. MetaNetX matches it to
   MNXR138860, from SEED rxn09632. Only the attribution is wrong; the H+ itself fits pH 7.3.
5. **"Spec values" as a guard** implies enforcement that does not exist (C8).
6. **The Strategy B origin of the TMP mapping** was not reproduced. What the auditor did confirm:
   `metabolites.py:381-385` does an exact-name lookup, and the build log shows
   `R1341 … mnxm=MNXM1137754`.
7. **Minor wording:**
   - m643's shared rationale cites reactions that do not involve HET-P.
   - "IDs verified" covers BiGG only through MetaNetX.

## (b) Outside the authorized scope

- **Model:** none. The diff is exactly the authorized targets; there are no H+ edits and names are
  unchanged.
- **Cosmetic:** the new annotations omit the short mirror keys (`keggc`, `biggm`) that most species
  carry. Only `tools/lipid/lipid_moiety_ledger.py` reads them, and it tolerates their absence.
- **Process:** no commit; REPORT §7 was still pending; no STATE entry yet.

## (c) Code defects

1. **`thiamine_phosphate_chemistry.py:48-55`:** the after-values are not pinned in code. Edited JSON
   values are applied silently.
2. **`:88-96`:** only stoichiometry is pinned, not partner formulas and charges. If ATP and the other
   partners migrate first, the step still writes the +1 forms, and its notes then misstate the balance.
3. **`main.py:608` runs before `main.py:768-772`, with no lock:**
   - Only step order preserves "no H+ edits". `balance_protons_and_water` on the new model would
     remove the H+ from R752, R408 and R750, since `{H: 1, charge: 1}` satisfies the gate.
   - The `lock_proton_water_stoichiometry` note is the existing way to prevent this; the model uses it
     once.
4. **Low-severity:**
   - `isinstance(charge, int)` accepts `True`.
   - A missing key raises `KeyError`, though still before any edit.
   - `notes_before` for R1342 embeds the TPP notes verbatim, which couples the two data files.

## (d) The protonation convention

1. **The KEGG forms are an interim, non-target convention.** The pipeline's declared convention is
   Rhea/ChEBI pH 7.3. These species, and HMP-PP, have no rows in `metabolite_microspecies.csv`.
2. **Effect on the migration audit (verified):**
   - `audit_component_migration(['atp','adp','amp','phosphate','diphosphate'])` now lists R1337,
     R1342, R751, R753 and R754 as regressions: 292 against 287 on the baseline, where those
     reactions were "unresolved".
   - The thiamine species have to be added to the component by hand.
   - The component was not ready before either (706 unresolved).
3. **The H+ hazard in (c)3 is specific to the KEGG choice.** A lock note on R752, R408 and R750 would
   make "no H+ edits" durable.
4. **R754 and R1343/R1344:**
   - R754's KEGG, MetaNetX and Rhea annotations (R10712, MNXR114031, rh:47848) describe the
     carboxy-thiazole-P reaction, not the model's HET-P equation.
   - R1343 and R1344 stay unresolved under every convention.
5. **The source's exchange name** "thiamine(1+) diphosphate(1-)" implies yet another microspecies.

**Partly checked or not checked:**
- the seed values (log headers only);
- the code identity of the "before" pytest run;
- the Strategy B origin of the TMP mapping;
- whether m889 is ribose-5-P or ribitol-5-P;
- the biology of R753 and of R1343/R1344.
