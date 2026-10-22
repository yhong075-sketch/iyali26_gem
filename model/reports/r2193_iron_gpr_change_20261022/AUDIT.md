# Independent audit: R2193 iron GPR change (2026-10-22)

A fresh read-only reviewer, with no history of this work, audited the v2 candidate (`2cf19c69`)
against claims C1–C10 in the brief.
- It made no repository changes.
- It used 15 LPs (limit 30), no builds, two focused-test runs (seeds 1 and 7) and 2 web look-ups.
- Its scratch files are in this session's scratchpad (not kept).

Its report is reproduced below, lightly reformatted. The main session's responses are in
[REPORT.md §8](REPORT.md#8-independent-audit-and-responses).

## Verdicts

| Claim | Verdict | Evidence checked |
| --- | --- | --- |
| C1 Baseline and before-state | supported | `shasum` 7f836c3a. R2193 m1273[C_ex] → m948[C_cy], bounds (0, 1000), one geneProductRef G_YALI1F07747g, no notes. None of the 4 new genes, nor their YALI0/NCBI/UniProt identities, is in the baseline. |
| C2 Controls reproduce the baseline | supported | `control_no_flag.xml` and `control_no_flag_v2.xml` both 7f836c3a; argv equal to the baseline's apart from the output path. |
| C3 v2 seed identity, first build | supported (provenance note) | v2 and its seed-1 twin are both 2cf19c69. 4d738403 vs b3d66568 and vs v2 differ only inside the four new geneProduct blocks; all three give the same hash after line sorting. `PYTHONHASHSEED` is not recorded in build.json or the logs, only in commands.txt. |
| C4 Diff only the targets | supported | An independent element-by-element XML diff shows R2193 (notes + GPR) and 4 geneProducts added at the end in sorted order, each with sbo/ncbigene/kegg.genes/refseq/uniprot. G_YALI1F07747g unchanged; everything else identical. |
| C5 Guard behaviour | partially supported | Fail-closed, idempotent and roundtrip pass; sorted order confirmed. Gaps in issue 6. |
| C6 Full pytest | supported | Per-test outcomes differ from the previous run only by the 5 new passes (230/11/17/31); subtests 558 → 585. |
| C7 Identifiers and evidence | partially supported | NCBI GeneIDs/RefSeq proteins, UniProt and KEGG cross-references all correct (live UniProt and KEGG). Q6CA12 is 589 aa vs 644 aa in W29, as noted. Primers match only XM_502497.1. Proteomics detection reproduced. **"Reciprocal best hit" is wrong for the D08 pair** (issue 2). |
| C8 Knockout comparison | supported (caveat in issue 1) | Recomputed from `ko_compare.tsv`: same 20 infeasible genes; max knocked flux 0; WT equal; class changes exactly the 5 genes at all cutoffs; scores match; in-model positives 322 → 324 (92 TP + 219 FN + 11 unresolved). |
| C9 Basis reset fixes the artifact | partially supported | `method_check.json` shows the warm start at WT with violation 3.12e-7 and the reset at 0, but it does not name R792. A rerun (2 LPs) shows R792 is the only closed reaction, at 3.124e-7. `method_check.py` calls `solve()` without the knocked set, so under the current check it would report the stale point as passed. The baseline sweep matches the heme report's corrected scores with no fresh re-solves. |
| C10 Stops and the mis-specified threshold | partially supported | The 7 stopped-run flags are confirmed, with the stated magnitudes and none on a knocked-out reaction. "Ordinary noise" is wrong: the pattern is systematic (issue 1). The narrowed re-solve rule has no recorded user approval (issue 5). The user's decision quote cannot be checked against files. |

## Issues

1. **Medium-high: the targeted check cannot see GLPK skipping the heme demand.**
   - 43 of the 48 baseline points with |S v| > 1e-9 have a residual of exactly 1e-6 × KO growth.
   - On m438 (protoporphyrinogen), GLPK zeroes R61…R240 and takes protoporphyrinogen from nothing
     within its tolerance, so the heme a term is skipped.
   - Gurobi at 1e-9, run on one gene from each of the 11 distinct-growth groups, carries the heme
     route and gives the same growth (about 1e-12). So no score changes here.
   - The 1e-6 general tolerance is above the largest possible skip residual (3.1e-7). A heme-step
     knockout landing on such a point would pass and be called non-essential. Here all heme-route
     genes returned 0.
2. **Medium: "reciprocal best hit" is false for D08684g and D08564g.**
   - FET3's top hit is Q6C9Y8 (553 bits) ahead of Q6CA12 (536). FTR1's top hit is Q6C9Y7 (392) ahead
     of Q6CA15 (346).
   - The D08 genes are one-way best hits (paralogs).
   - The wrong claim is in the curated data, the R2193 model notes, REPORT §2, TASK, and the
     investigation's `candidate_genes.tsv`. Fixing it needs a data fix and a rebuild. The W29
     D08684g protein (644 aa) was never searched.
3. **Low-medium: Ye et al. 2026 is overstated** ("overexpression raised intracellular iron"). The
   study co-expressed VHb with Ftr1, with no Ftr1-only control. The data file's
   `same_species_evidence` field says this correctly.
4. **Low-medium: `method_check.json` does not name R792, and the check API has a trap.** The default
   `knocked=()` in `checked_solve.solve()` silently disables the stale-point test.
5. **Low: the narrowed re-solve rule was a unilateral change.** Leaving infeasible knockouts
   unresolved without a re-solve relaxed acceptance criterion 5, and no approval is shown.
6. **Low: guard gaps.**
   - Called inside `with model:`, the step leaves 4 annotated unlinked genes and the notes behind
     after exit, because genes were appended by hand and notes are not history-tracked. A re-apply
     then fails closed.
   - The reviewed model's SHA and the old gene's name are not enforced.
   - The post-edit RuntimeError fires after the edits.
   - The test's "no edits" state leaves out gene names, gene order and model-level elements.
7. **Low: provenance gaps.**
   - `PYTHONHASHSEED` is not in build.json.
   - The experimental list's full SHA appears only as a prefix (`1e887f5a…659b` matches).
   - `flag_diagnosis.py` divides by cf0825b8's WT, not 7f836c3a's.
   - The diagnosis points differ numerically from the rerun points for the same genes; a basis reset
     does not make vertices path-independent.
8. **Low: wording.**
   - Acceptance criterion 4 ("within 1e-9 of 0.3124263") fails literally (4.4e-8, from rounding);
     "WT unchanged" holds.
   - The 11 unresolved positives and the list's positive coverage are not stated.
   - The two 16s (unresolved and mtDNA) are easy to confuse.
   - The plan's file and function names differ from the implemented ones.
9. **Low: budget overruns**, all recorded:
   - 6 builds against 3;
   - the stopped sweep, its 7 fresh re-solves, and 7 diagnosis LPs;
   - the diagnosis LPs counted against the re-solve budget.
10. **Information: process state.**
    - STATE.md entry not yet written; REPORT §8 pending.
    - Nothing committed.
    - The builds used a newer `heme_iron_biomass.py` than the one that built 7f836c3a; no effect,
      since both controls reproduce 7f836c3a.
11. **Information: biology.**
    - R2193 is a lumped Fe²⁺ uniport. Fet3/Ftr1 is an O₂- and Cu-dependent oxidase–permease that
      moves Fe³⁺ (established, uncited). The GPR follows the yeast-GEM convention; the mechanism is
      not modelled.
    - That the D08 genes form a working pair is inferred.
