# Step 1 report — Europe PMC retrieval for tier-1 gene–reaction pairs (2026-10-19)

Scope, inputs and method: [TASK.md](TASK.md). No LLM was used. No model file was changed.

## Run (observed 2026-10-19)

- Command as in TASK.md; 2,289 network requests, 188 cache hits (shared genes/articles),
  exit 0. All responses cached under
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/gem_review_lit_20261019/cache/` (2,289 files);
  a re-run makes no network request.
- 442 reactions → 673 pairs over 437 reactions and 158 genes. Five reactions have no GPR and
  produce no pair: R530, R1727, R1902, R2003, R2274 (R530/R1727 are the duplicate LPCAT pair
  already noted in curation history).

## Hit-count table ([hit_counts.csv](hit_counts.csv))

| pair_hits bucket | 0 | 1–4 | 5–19 | 20+ |
|---|---|---|---|---|
| pairs (n = 673) | 461 | 164 | 28 | 20 |

- 177 pairs (26%) have at least one open-access full text fetched; 168 have ≥ 1 full-text
  paragraph that mentions a gene or reaction term.
- Gene-level (`gene_hits`, locus tags + symbols only): 53 of 158 genes have **zero** Europe PMC
  hits; 24 genes have 100+ (lipid-engineering workhorses: GPD1, POX2, DGA1-type genes).

Per subsystem (first subsystem label of each reaction):

| Subsystem | pairs | pairs with pair_hits > 0 | pairs with full text |
|---|---|---|---|
| Fatty acid degradation | 132 | 82 | 68 |
| Sphingolipid metabolism | 118 | 0 | 0 |
| Alkane metabolism | 70 | 12 | 0 |
| Glycerophospholipid metabolism | 63 | 11 | 7 |
| Glycerolipid metabolism | 62 | 22 | 22 |
| Steroid biosynthesis | 54 | 5 | 5 |
| Fatty Acid Biosynthesis | 44 | 20 | 20 |
| Fatty acid metabolism | 40 | 24 | 19 |
| Fatty acid elongation | 34 | 15 | 15 |
| Fatty acid biosynthesis | 26 | 2 | 2 |
| Dicarboxylic acid biosynthesis (C10/C12/C16) | 22 | 17 | 17 |
| Arachidonic acid metabolism | 6 | 0 | 0 |
| Fatty acid desaturation | 2 | 2 | 2 |

Highest `pair_hits`: R348 GPD1/YALI1B04433g (76; 24 OA), the POX2/YALI1F14495g acyl-CoA
oxidase family R1490–R1498 (57 each; the ten reactions share one gene and one EC, so one
retrieval serves all), POX3/YALI1D32754g R100–R102, R1487–R1488 (30 each).

## What the zeros mean (verified from the cached queries, not inferred)

1. **Sphingolipid metabolism, 0/118.** Every pair had a pair query; 106 of the 118 pairs have
   `gene_hits = 0`, i.e. the locus tags themselves never occur in Europe PMC, and none of these
   genes has a UniProt gene symbol. Reaction names such as "IPS phospholipase C" or
   "Ceramide-1 synthase" are model-internal and do not match literature wording. This is
   genuine absence of locus-level literature, not a query failure; whether ortholog names
   (ISC1, LAC1, …) would find organism-level evidence is untested.
2. **Gene symbols.** 125 of 158 genes carry a `uniprot` annotation. Of those UniProt entries,
   **99 are deleted** (`entryType: Inactive`, reason "Not part of a reference proteome"; the
   `A0A1H6…` CLIB89 TrEMBL accessions), 5 are active with only an ordered-locus name, and 21
   have a gene name/synonym. So symbol-based search covered 21 genes. This is a model
   gene-annotation issue (checklist G01), outside this task's scope; recorded, not changed.
3. **Reaction terms are weak for non-EC reactions.** Example: R1867 "alkane oxidase (C16)",
   ALK2, `gene_hits = 20`, `pair_hits = 0`, because papers say "alkane hydroxylase" or
   "cytochrome P450". Pairs whose reaction has no EC depend entirely on the model's name.

## Limits

- Europe PMC only; no PubMed-only or publisher-only full texts. At most 25 abstracts and 10
  full texts per pair.
- `pair_hits` counts articles mentioning gene and reaction terms anywhere, not evidence; the
  evidence level comes from step 2.
- Tier-1 membership is the subsystem list in `tier1_lipid_subsystems.csv`; the five borderline
  subsystems were excluded pending the user's decision.

## Step 1 verification level

- Code: 5 mocked tests pass (`tests/test_gem_review_lit_retrieval.py`); ruff clean.
- Numbers above were computed from `hit_counts.csv` and the cache in this session.

---

# Step 2–3 report — claim extraction and review sheets (2026-10-19)

Setup and authorization: [TASK.md](TASK.md), section "Step 2–3". Code:
`platform/gem_review/lit/models.py`, `extract.py`, `report.py`; CLI `extract` and `report`.

## How the extraction works

- Only the **212 pairs with `pair_hits > 0`** go to the model; the 461 zero-hit pairs are
  recorded as `no_hits` with `best_evidence_level = none` and no model call.
- Per call the agent (pydantic-ai 2.58, `anthropic:claude-opus-5-5`) receives the pair's
  chunks as a JSON block in its instructions and the question "Gene X (also written …aliases…);
  reaction R \"name\" (EC). Does the retrieved text establish that this gene product catalyses
  this reaction?". It must answer only from the chunks and returns a `ClaimSet`.
- Every claim is validated: `pmid` must be in the pair's retrieved PMIDs, `chunk_id` must be a
  retrieved chunk of that PMID, and the quote (≤ 25 words) must occur verbatim
  (whitespace-normalised) in that chunk. A failing answer goes back to the model with the
  reasons; after 2 retries the pair would be recorded `EXTRACTION_FAILED`.
- `best_evidence_level` = max over *supporting* claims of the claim's level, where a claim with
  `organism_is_yarrowia = false` counts at most as homology (1). Contradicting claims never raise
  the level. `no_organism_specific_evidence` = no supporting Yarrowia claim.
- One model call per distinct (gene, full reaction name, EC codes, chunk set); 212 pairs → 191
  calls. Chain-length variants of the same enzyme are *not* merged (see below for why).

## Run (observed 2026-10-19)

| Item | Value |
|---|---|
| Pairs sent to the model | 212 (191 distinct questions) |
| Status | `ok` 212, `EXTRACTION_FAILED` 0, `MODEL_ERROR` 0 |
| Validation retries | 15 of 191 calls needed one retry (207 requests in total) |
| Tokens | input 2,730,502; output 135,146; cache reads 0 (`usage.jsonl`, per pair) |
| Approx. cost | ≈ $13.6 at $4 / $20 per MTok (list price cached 2026-10-06 in the claude-api skill) |
| Validated claims | 357 over 138 pairs; 74 pairs with hits returned no claim |

`best_evidence_level` over the 212 extracted pairs: none 79, homology 75, physiological 12,
genetic 34, biochemical 12. Over all 673 tier-1 pairs: none 540.

Per subsystem (extracted pairs, L0…L4 = evidence level):

| Subsystem | L0 | L1 | L2 | L3 | L4 |
|---|---|---|---|---|---|
| Fatty acid degradation | 24 | 35 | 10 | 5 | 8 |
| Fatty acid metabolism | 8 | 12 | 1 | 2 | 1 |
| Glycerolipid metabolism | 8 | 10 | 0 | 2 | 2 |
| Fatty Acid Biosynthesis | 20 | 0 | 0 | 0 | 0 |
| Fatty acid elongation | 1 | 14 | 0 | 0 | 0 |
| Alkane metabolism | 6 | 0 | 0 | 5 | 1 |
| Glycerophospholipid metabolism | 6 | 1 | 1 | 3 | 0 |
| Dicarboxylic acid biosynthesis (C10/C12/C16) | 5 | 0 | 0 | 12 | 0 |
| Steroid biosynthesis | 1 | 1 | 0 | 3 | 0 |
| Fatty acid desaturation / biosynthesis | 0 | 2 | 0 | 2 | 0 |

## review.xlsx ([review.xlsx](review.xlsx))

- **literature** (673 rows): hit count, full-text Y/N, status, best level (number and name),
  `no_organism_specific_evidence`, claim counts, `claims_disagree`, cited PMIDs, the validated
  quotes, token usage.
- **review_queue** (618 rows), sorted by tier then reason: `contradicted` 4, `claims_disagree`
  3, `low_evidence` 150 (level ≤ 1 after extraction), `no_hits` 461. The spec's "level ≤ 1"
  rule therefore queues almost everything; the `reason` column is what makes it usable.
- Compartment disagreement uses a small synonym map (`peroxisomal`/`peroxisomes`/`peroxisome`
  → one compartment); without it 29 spelling-only "disagreements" appeared. The 3 remaining:
  R1886 and R2261 (YALI1E18480g, aldehyde dehydrogenases) and R88 (YALI1C15991g, acetyl-CoA
  carboxylase).

## Findings worth a human look (candidates, not verdicts)

1. **POX3 (YALI1D32754g) on C14–C26 acyl-CoA oxidase reactions R99, R100, R101, R102:** every
   validated claim contradicts the assignment — quotes from PMID 9848229, 10464181 and 10593257
   state that ACO3/POX3 encodes the *short-chain* acyl-CoA oxidase and POX2 the long-chain one.
   The same gene on short-chain reactions (R97 decanoyl, R1488 hexanoyl) gets biochemical-level
   support from the same papers. This is a GPR-specificity candidate for the acyl-CoA oxidase
   family; the model is unchanged. Note: the model lumps POX3 onto all chain lengths, so this
   is exactly the paralog-assignment case (checklist G05) the literature layer was meant to find.
2. **Sharing one call across chain lengths would have hidden (1):** an earlier version keyed
   calls on the cleaned reaction name; the smoke test showed the C16 verdict being inherited by
   the C4 reaction, so the key now includes the full name (80 → 191 calls).

## Judgment quality — what is verified and what is not

- **Verified by code for every claim:** the PMID is in the retrieved set, the chunk belongs to
  that PMID, the quote is verbatim in the chunk, ≤ 25 words.
- **Assigned by the model, not verified:** `supports_gene_reaction`, `evidence_type`,
  `organism_is_yarrowia`, `compartment`. Three level-4 pairs were read against their chunks in
  this session:
  - R1887 HFD4 (PMID 25315778): the quote about FALDH activity in *E. coli* extracts expressing
    "the HFD genes" is genuine in-vitro evidence but collective for HFD1–4, not HFD4-specific.
  - R1771 DGA1 (PMID 27780240, purified active Dga1p): level 4 is justified.
  - R296 POX3 on "fatty acid oxidation (C18:2)": **over-credited.** One quote is a generic
    textbook sentence about acyl-CoA oxidase (no gene), and the other ("POX2 and POX3
    respectively encode long- and short-chain acyl-CoA oxidases") argues *against* POX3 on
    C18:2 yet is labelled supporting, level 4.
  So in 1 of 3 sampled pairs the label is wrong while the quotes are right. Levels in
  `review.xlsx` are LLM labels on verified quotes and need human confirmation before any use
  as evidence; this is why the queue exists.
- Not independent validation: the same literature set is used to flag and to support.

## Limits and possible improvements

- Prompt caching was not used (`cache_read_tokens` 0). The chunk block is per pair, but the
  same gene's chunks recur across its reactions; `anthropic_cache_instructions` would cut input
  cost for gene families such as POX2/POX3.
- Generic sentences ("the first step is catalysed by acyl-CoA oxidase") pass the verbatim check
  but carry no gene-specific evidence. A rule "a supporting quote must mention the gene by one
  of its names" would be a cheap second validator; not implemented (not in the spec).
- 20 extracted "Fatty Acid Biosynthesis" pairs returned no claim at all; not examined.
- `evidence_type` descriptions follow the checklist's Evidence Levels sheet; the model's
  reading of "genetic" vs "physiological" was not audited beyond the three pairs above.

## Verification level

- Code: 15 mocked tests pass (`tests/test_gem_review_lit_extract.py`,
  `tests/test_gem_review_lit_retrieval.py`), including a claim with a PMID outside the
  retrieved set being rejected, a quote not in its chunk being rejected, and retry exhaustion
  producing `EXTRACTION_FAILED`. Ruff clean.
- Full platform suite after step 1: 162 passed, 19 failed, 14 errors, 31 skipped — all in nine
  pre-existing files for environmental reasons (missing `.claude/worktrees/iyali26_gem/model.xml`
  path, `scipy` not installed, unavailable CoA-protonation source model, recorded-hash
  assertions on untouched scripts); none imports `gem_review`. Re-run after steps 2–3: see the
  STATE.md entry for the final counts.
- All numbers in this section were computed from `out/lit/*.review.json`, `usage.jsonl` and
  `review.xlsx` in this session.
