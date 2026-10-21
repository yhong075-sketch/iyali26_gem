# R1964 GPR review: should YALI1F03452g (or other genes) catalyse the orphan arabinitol dehydrogenase? (2026-10-19)

## Authorization

Delegated by another session of this group (shared notes folder
`side-session-notes/local_11292725-…`), following the R2190 identity review:

> Read-only evidence review, then a plan for the user's approval. No model, data or code edits, no
> build, no LP, no commit.

Questions:
1. Is the evidence enough for an OR-member GPR on R1964 under AGENTS.md §5 ("an OR member needs
   evidence that it catalyses the step on its own")? Look for direct product identification (Wang
   2020 supplementary material, other *Y. lipolytica* papers).
2. Regiochemistry: does the enzyme oxidize D-arabinitol at C-2 (D-ribulose) or C-4 (D-xylulose)?
   This decides between R2190 options A and B in the R2190 report.
3. Interactions with the pending R2190 (ribulose product) and R2041 (forward-only) plans.

Done = this TASK.md, REPORT.md, and a plan presented to the user.

## Starting state

- Worktree `claude/reverent-wilbur-4dc5b6` at `285eaee` (clean at start).
- The main checkout is on `restructure/model-platform` at `2a34475`. Only its `.venv` interpreter
  (cobra 0.30.0) was used, to load SBML statically.
- Model under review: a56f3649 (below). It is not a designated reference model
  (`model/candidates/README.md`). Reference model: none designated.
- Sibling work read, not changed (all uncommitted):
  - `claude/clever-austin-f6a04e`: `model/reports/r2190_identity_review_20261019/`;
  - `claude/bold-chaplygin-1c16df`: R2041 direction review and the implemented R2041 forward-only
    candidate `da0c50ba…` (known from the hand-over note only; the file was not opened);
  - `claude/orphan-reaction-newest-model-b8ddc8`: `model/reports/orphan_reactions_20261019/`
    (its recorded `loop_check.txt` and `flux_evidence.csv` are cited, not re-run).

| Input | SHA-256 |
| --- | --- |
| a56f3649 `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml` | `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515` |
| Source `model/source/iyali26.xml` | `5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d` |
| Colleague audit `model/curation/lipid_unlump/r1521_kegg_mnx_frontier_source_audit.json` (read only) | `2014a5e7fb51036036ff734def8a93cee3b5b15514ef5c54a971585c9e275fe9` |
| `$RESEARCH/reference/metanetx/reac_prop.tsv`, `reac_xref.tsv`, `chem_prop.tsv` (MNXref 4.5) | as listed in the R2190 TASK.md; not re-hashed (800 MB) |
| `$RESEARCH/reference/locus_map/s2_metabolic_genes.csv` (YALI1↔YALI0 map) | `f7f25bb3e4f3c720d4b0bd534309c5c5f3cfa091bbd22383858d71f8600c1c2e` (recorded in the R2190 TASK.md) |
| `$RESEARCH/state/essentiality/repository/consensus_essential_genes.csv` | `1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b` |
| `$RESEARCH/state/media/sd_leu.csv` | `ed176d26a373f98cc413ed2e32a71f5f060a06e343f90f7db25cd32eff268e85` |

`IYALI26_RESEARCH_ROOT` is unset in this session. `$RESEARCH` above is
`/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`, the folder that a56f3649's
`build.json` names as `--research-root` (per the R2190 TASK.md). Only reference tables were read from
it. Nothing was built, so the builder's fail-closed rule was not triggered.

## Scope and budget

- **Static only:** SBML loads, element and charge sums, connectivity, reading of recorded outputs,
  and a Fischer-projection bookkeeping script. No LP, FVA, build or screen was run.
- **External look-ups (2026-10-19 local):** PMC full text of Wang et al. 2020 (JATS XML) and of
  Adamczyk et al. 2023; PubMed abstracts (13); NCBI Gene/Protein; UniProt REST; KEGG REST; ENZYME;
  Rhea REST; one ACS page (HTTP 403, unread).
- **Outputs:** this TASK.md, REPORT.md and `sources_sha256.txt`. Scripts and fetched files stay in the
  session scratchpad.
- **Not done:** `model/STATE.md` was not edited (REPORT §9 has the entry to paste, as in the sibling
  R2190/R2041 reviews); the Wang 2020 supplementary `.docx` was not downloaded; no commit.
- **Stop:** any edit to model, data or code is out of scope.
