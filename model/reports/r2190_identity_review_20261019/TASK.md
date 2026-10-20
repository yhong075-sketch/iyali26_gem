# R2190 identity review: product versus EC 1.1.1.250 (2026-10-19)

## Authorization

Delegated by the session that ran the R2041 direction review (`r2041_direction_review_20261019`,
worktree `claude/bold-chaplygin-1c16df`), which found the conflict as a side item:

> Read-only review, then a plan for the user's approval. No model, data or code edits, no build, no
> LP, no commit.

Questions:
1. What do R2190's own annotations (KEGG R07134, RHEA:17389, MNXR197981) resolve to? The
   colleague-owned `model/curation/lipid_unlump/r1521_kegg_mnx_frontier_source_audit.json` may be
   read, not edited.
2. What is YALI1F03452g (NCBI, KEGG, UniProt; YALI1 and YALI0 locus numbers are separate)? Is there
   *Y. lipolytica* literature on D-arabinitol dehydrogenase?
3. Which side is wrong, the product or the EC? Keep it a candidate until the IDs are opened.
4. What follows for the model: there is no free D-ribulose metabolite, and R2190 links D-arabinose
   to the arabinitol node and to the ATP-neutral loop with R2041/R2042/R803/R714/R1964.

Done = this TASK.md, REPORT.md, and a plan presented to the user.

## Starting state

- Worktree `claude/clever-austin-f6a04e` at `285eaee` (clean at start).
- The main checkout is on `restructure/model-platform` at `2a34475`. It was not used except for its
  `.venv` interpreter (cobra 0.30.0) to load SBML.
- Model under review: a56f3649 (below). It is not a designated reference model
  (`model/candidates/README.md`). Reference model: none designated.
- Sibling work read, not changed:
  - `claude/bold-chaplygin-1c16df`: `model/reports/r2041_direction_review_20261019/` (uncommitted);
  - `claude/orphan-reaction-newest-model-b8ddc8`: `model/reports/orphan_reactions_20261019/`
    (uncommitted; its recorded loop check is cited, not re-run).

| Input | SHA-256 |
| --- | --- |
| a56f3649 `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml` | `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515` |
| Source `model/source/iyali26.xml` | `5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d` |
| Colleague audit `model/curation/lipid_unlump/r1521_kegg_mnx_frontier_source_audit.json` (read only) | `2014a5e7fb51036036ff734def8a93cee3b5b15514ef5c54a971585c9e275fe9` |
| `$RESEARCH/reference/metanetx/reac_prop.tsv` (MNXref 4.5, 2025-08-13) | `8582cc187d03ce127f8e914f8f298282ac9036f1448918117af044ac55b980db` |
| `$RESEARCH/reference/metanetx/reac_xref.tsv` | `2610e2998f582043a4e9d1a8ce42b7a8f4e2480e54d03ded26a7e83fd401f229` |
| `$RESEARCH/reference/locus_map/s2_metabolic_genes.csv` (YALI1↔YALI0 map) | `f7f25bb3e4f3c720d4b0bd534309c5c5f3cfa091bbd22383858d71f8600c1c2e` |
| `$RESEARCH/reference/ncbi/clib89_feature_table.txt` | `15a15c73ca0817d2fe64d055101a627084278fbdd1ae004eb78ee98978c0a7f2` |
| `$RESEARCH/reference/kegg/yli_genes.tsv` | `33cb8ae3689f9774f66e2c22e4c3caeab9fa55db89941fd4aa791be5c0d5bf4f` |

`IYALI26_RESEARCH_ROOT` is unset in this session. `$RESEARCH` above is
`/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`, the folder that a56f3649's
`build.json` names as `--research-root`. Only reference tables were read from it. Nothing was built,
so the builder's fail-closed rule was not triggered.

## Scope and budget

- **Static only:** SBML loads, element and charge sums, connectivity, and reading of recorded
  outputs. No LP, FVA, build or screen was run.
- **External look-ups (2026-10-19 local, 2026-10-20 UTC):** KEGG REST and SSDB, ENZYME (ExPASy),
  ExplorEnz, Rhea, MetaNetX, NCBI E-utilities (gene, protein, PubMed), UniProt REST, and three
  literature pages (§9 of REPORT.md).
- **Outputs:** this TASK.md and REPORT.md. Analysis scripts and fetched pages stay in the session
  scratchpad (listed with SHA-256 in REPORT §9).
- **Not done:** `model/STATE.md` was not edited (REPORT §10 has the entry to paste), and no commit.
- **Stop:** any edit to model, data or code is out of scope.
