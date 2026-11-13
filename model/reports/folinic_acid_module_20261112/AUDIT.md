# Independent audit: folinic-acid module, 2026-11-12

## How the audit was run

- **Auditor:** one fresh, read-only general-purpose agent. It was given TASK.md, REPORT.md (first draft), the data,
  the code, the logs, the outputs and both model files. It was not told the expected verdicts.
- **Method:**
  - recomputed every SHA;
  - re-diffed the two XML files as text;
  - re-read the models with cobra (no LPs);
  - exercised the guard in memory;
  - re-fetched every external source (Rhea, UniProt, PubMed, KEGG, ENZYME, ChEBI, BiGG, MetaNetX, NCBI).
- **Scope:** the audit covered the first build `d4832b29`. Pytest claims were pending at the time.

## Supported

57 claims were checked.

**Model and builds:**
- SHAs and byte identity;
- the argv;
- that only the intended entities changed (independent text diff: hunks only in R748, R990, R_THFAT and group50);
- R_THFAT's definition, and its mass and charge balance.

**Code:**
- R990 = R64 + R543;
- the gene links;
- the guard: fails closed before any edit, all-or-nothing, a rerun after reload is a no-op, context exit restores
  the model;
- the in-memory result equals the built file.

**Results:**
- every number in the result tables;
- the sweep counts and scores;
- the absence of 5-formyl-THF from all biomass reactions.

**External sources:** every external claim, with the corrections below.

## Problems raised and how they were resolved

| # | Severity | Finding | Resolution |
| --- | --- | --- | --- |
| 1 | medium | `main` moved to `efd7f780` (PRs #49–#51: the R2121 and R1172 closures, the AGENTS.md §5 inactive-reactions rule) before the builds. TASK's "= main" is stale, and `main.py`/`cli.py` will conflict textually | Reported in REPORT §7 and TASK §1; §4 checks the new rule and names one tension (R990 closed to (0, 0)) as a user decision. No merge or rebase: that needs the user's instruction |
| 2 | medium | Gene names attributed to RefSeq; NCBI says "uncharacterized protein", and the names are KEGG GENES definitions | Fixed in REPORT §2. The KEGG "cytosolic" label on `E19557g` is added as a weak hint to the untested-localisation row |
| 3 | medium | "Same inherited glycine-cleavage rule" overreaches: Yeast-GEM r_1030 is GCV1/2/3 + LPD1, iMM904 THFATm is GCV1/2/3 + PDA1/PDB1 | Fixed in REPORT §2 and in the R748 model note (curation data reworded, rebuilt as v2 `85c623c8`); "inherited" is labelled inferred |
| 4 | low–medium | The FVA LPs were not run through the checked solve; candidate R64 min −8.7e-10 | Stated in REPORT §3 and §5 |
| 5 | low–medium | Diff-script coverage overstated | REPORT §3 now lists what the script compares. Byte-level text diffs are added to the logs and cited |
| 6 | low | "Under every medium" | Qualified to "while R991's lower bound is 0, as in every medium under `model/conditions/`" |
| 7 | low | The isomerase search log missed R02300/R08701 and had no ENZYME search | Extended: KEGG `link/reaction/C03479`, Rhea `chebi:57457`, an `enzyme.dat` scan ([logs/isomerase_search_20261112.txt](logs/isomerase_search_20261112.txt)). R990 note reworded to "No EC, KEGG or Rhea entry describes … (searched 2026-11-12)" |
| 8 | low | Dangling AUDIT.md link and "see below" | This file; REPORT §3 and §6 complete |
| 9 | low | Positive coverage and out-of-scope positives missing; PYTHONHASHSEED not in `.build.json`; before.json without SHA | Coverage added (326 of 1,612; 0 out of scope). The seed is recorded in the log names and REPORT §1. The before.json source SHA is in `builds_sha256.txt` (`e584151f`) |
| 10 | low | Wording: "the right enzyme", "biological source", "R990 no longer bypasses", the R811/R991 table cell, the R990 bypass note untagged | Reworded in REPORT §4. R811 and R991 have separate rows. The R990 model note now starts the bypass sentence with "Inferred:" (v2) |

## Rebuild after the note fixes (findings 3 and 7)

- **v2** (`..._folinic_acid_v2_20261112.xml`, `85c623c8…`): seeds 0 and 1 are byte-identical.
- **v2 vs `d4832b29`:** cobra-level diff = the notes of R748 and R990 only (2 entries); byte-level diff = 4 lines.
  Same optimisation problem, so the LP results carry over.
- **v2 vs `e584151f`:** the same 11 entries as before.
- **Tests:** the focused test passes on v2. The full suite was not re-run for the note change.

## Not caused by this round

Recorded in REPORT §5:
- the libSBML `R_R304` notes error, present in every build;
- deleted UniProt accessions on the SHMT, T-protein and FAU1-type genes;
- the BiGG `THFAT` ID on R748.
