# gem_review literature layer — step 1, Europe PMC retrieval on tier 1 (2026-10-19)

## Authorization

User's pasted spec (2026-10-15): "Extend gem_review with a literature-review layer using
Pydantic AI … Before writing the agent, run step 1 on tier-1 reactions only and show me the
hit-count table." Follow-up answers (2026-10-15 and 2026-10-19):

- `gem_review` did not exist anywhere (this worktree, every branch, all 21 worktrees, `~`);
  the user said to create it **under `platform/`**.
- Tier definitions come from the user's checklist workbook
  `~/Downloads/GEM_reaction_review_checklist.xlsx`, sheet "Prioritization": tier 1 =
  "Reactions in lipid pathways".
- Europe PMC (and UniProt gene-symbol) queries for step 1: **permitted** ("3, permit").
- Direction (2026-10-19): build the spec first (retrieval → extraction → review.xlsx), then add
  a question-driven "logic chain" layer (quote → inference → verdict) on top.

**Scope limits:** read-only review tooling. No model XML, GPR, medium, benchmark or test
baseline is changed. No commit or push unless instructed. Lipid-unlump files untouched.

## Starting state (verified 2026-10-19)

| Item | Value |
|---|---|
| Worktree / branch | `claude/gem-review-lit-layer-861d41` @ `4d97b82` (clean at start) |
| Model read | `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml`, SHA-256 `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515` (newest candidate; no reference model is designated, `model/candidates/README.md`) |
| Research root | `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research` (from `~/.zshrc`; passed with `--research-root`, the session environment does not export it) |
| YALI1→YALI0 crosswalk | `scripts.gem_annotate.locus_resolver.LocusCrosswalk` over `reference/locus_map/iyli21_genes_vs_S2.csv` + `s2_metabolic_genes.csv`; fingerprint recorded in `manifest.json` |
| Gene symbols | UniProt REST (`/uniprotkb/<acc>.json`, `fields=gene_names`) from the model's `uniprot` gene annotation; 125 of 158 tier-1 genes carry one |

## Tier-1 definition used

`model/curation/gem_review/tier1_lipid_subsystems.csv` (include=yes rows): both spellings of
Fatty acid biosynthesis; Fatty acid degradation / desaturation / elongation / metabolism;
Glycerolipid, Glycerophospholipid and Sphingolipid metabolism; Steroid biosynthesis;
Arachidonic acid metabolism; Alkane metabolism; the three Dicarboxylic acid biosynthesis
subsystems. Borderline subsystems (Inositol phosphate, Terpenoid backbone, lipid-particle
transport, ketone bodies) are listed with include=no pending the user's decision.

Selection on the model above: 442 reactions, 673 gene–reaction pairs, 158 genes; 5 reactions
have no GPR and yield no pair (listed in `manifest.json`).

## Method (no LLM)

Code: `platform/gem_review/lit/` (`pairs.py`, `terms.py`, `client.py`, `chunks.py`,
`retrieve.py`, `__main__.py`). Command:

```
python -m gem_review.lit retrieve --model ../model/candidates/E5_coq9_resp_pkg_r132_20261012.xml \
  --tier 1 --research-root "$IYALI26_RESEARCH_ROOT" --task gem_review_lit_20261019
```

Per pair, two Europe PMC queries (all terms quoted):

- gene query: `(<YALI1 tag, both spellings> OR <YALI0 tag, both spellings> OR <UniProt symbols>)
  AND "Yarrowia lipolytica"` → `gene_hits` (hitCount only);
- pair query: gene query `AND (<EC codes> OR "<reaction name without parenthetical>")` →
  `pair_hits`, first 25 results (`resultType=core`, abstracts).

Full text (`fullTextXML`) is fetched for at most 10 open-access hits of the pair query.
Chunks: every abstract (`<pmid>:abstract`) plus body paragraphs mentioning a gene or reaction
term (`<pmcid>:p<n>`). A reaction with no EC and a name shorter than 4 characters gets no pair
query (`pair_hits` 0, no chunks) so the table shows where the reaction is unsearchable.

Every HTTP response (200/404) is cached by SHA-256 of (url, params) under
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/gem_review_lit_20261019/cache/`; ≤ 5 requests in
flight, ≥ 0.34 s between request starts, 3 retries on 429/5xx.

## Outputs

- `hit_counts.csv` (this folder): one row per pair.
- `manifest.json` (this folder): model SHA, crosswalk fingerprint, parameters, request counts.
- `$IYALI26_RESEARCH_ROOT/artifacts/tasks/gem_review_lit_20261019/out/lit/retrieval/<gene>__<reaction>.json`:
  per-pair queries, PMIDs/PMCIDs, full-text flag and chunks (input to step 2).

Tests: `platform/tests/test_gem_review_lit_retrieval.py` (Europe PMC and UniProt mocked with
`httpx.MockTransport`): selection refuses the whole model, query construction, chunk
filtering, end-to-end retrieval, and a second run served entirely from the disk cache.

## Step 2–3 authorization and setup (2026-10-19)

- User (2026-10-19): "Yes, implement step 2 (Pydantic AI extraction) as planned, on the 212
  pairs with hits." Then: run the smoke test, the full extraction and the report.
- Dependencies added to `platform/pyproject.toml` / `uv.lock`: `pydantic-ai-slim[anthropic]`
  (2.58.0), `openpyxl` (3.1.5).
- LLM: `anthropic:claude-opus-5-5` through pydantic-ai, `max_tokens` 16000, effort `medium`,
  `retries=2` on output validation, ≤ 5 concurrent runs. Key read from `ANTHROPIC_API_KEY`,
  supplied via the git-ignored `platform/.env` (user-created; a first key the user placed in
  `~/.zshrc` was rejected with 401 and then exposed in a tool output during diagnosis; the user
  revoked it and created a new Console key "literature review").
- Deviation from the pasted spec, agreed in the plan: the agent returns `ClaimSet` (claims
  only); `PairReview.best_evidence_level` and `no_organism_specific_evidence` are computed in
  Python from validated claims. `EvidenceType` is an `IntEnum`. The PMID-in-retrieved-set and
  quote-verbatim checks run in a pydantic-ai `output_validator` with the pair's chunks as run
  deps (a field validator cannot see them). Chunks reach the model as a JSON block in the
  instructions, never in the user message.
- Smoke test (3 pairs) first returned 0 claims because the question named the gene only by its
  YALI1 ID while papers use the symbol/YALI0 tag; the question now lists those aliases.
- One model call per distinct (gene, full reaction name, EC, chunk set); chain-length variants
  of one enzyme are judged separately because the same POX3 quotes contradict C16 and support
  C4.
