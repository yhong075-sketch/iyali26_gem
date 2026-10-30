# Evidence pipeline, M0: survey report (2026-10-29)

Scope and starting state: [TASK.md](TASK.md). Reviewer interface: [reviewer_interface.md](../../curation/evidence_pipeline/reviewer_interface.md) (written as `docs/reviewer_interface.md`, moved on 2026-10-29 because `/docs/` is git-ignored).
M0 only: read-only survey, no pipeline code, no computation, no model change, nothing committed.

**Verification level.** Every statement below was read this session from `origin/main` @ `f212874`, from
branch `codex/module-agent-workflow` @ `35c959b`, or from the research workspace, unless it is marked
**inferred** or **recalled** (taken from `model/STATE.md` without re-running). Counts come from scripts run
this session on the files named.

## 1. Repo map (as it matters to the spec)

| Spec need | Where it lives now |
| --- | --- |
| Model files | `model/source/iyali26.xml` (immutable build input), `model/source/iyli21.xml` (parent model), `model/candidates/*.xml` (47 XML, none designated reference), layout in `model/model.toml`, resolved by `scripts.gem_annotate.model_layout.MODEL` |
| Builder (only writer of model XML) | `platform/scripts/gem_annotate/` (64 files); science changes = curated data in `model/curation/<topic>/` + one guarded `apply_*` step + test |
| Tools / tests | `platform/tools/` (32 files), `platform/tests/` (67 files); pytest config in `platform/pyproject.toml` |
| Typed LLM layer | `platform/gem_review/lit/` (pydantic-ai; §3.7) |
| AI reviewer | Essentiality FN loop (Codex agents + `essentiality_evidence.py`); see `model/curation/evidence_pipeline/reviewer_interface.md` |
| ID maps | `$IYALI26_RESEARCH_ROOT/reference/locus_map/` + `scripts.gem_annotate.locus_resolver.LocusCrosswalk` |
| Reference DBs (local) | `$IYALI26_RESEARCH_ROOT/reference/`: `metanetx/` (chem/reac prop + xref), `kegg/` (`yli_*`), `expasy/` (`enzyme.dat`), `ncbi/` (CLIB89 feature table), `external_models/` (yeast-GEM, iMM904, iLB1027, eciYali5) |
| Datasets | `$IYALI26_RESEARCH_ROOT/state/essentiality/repository/`, `reference/w29_worland2024/` (§3.5) |
| Reports, state | `model/reports/<topic>_<date>/` (207 folders), `model/STATE.md` (newest first) |
| CI | `.github/workflows/memote.yml`, `diff-pr.yml`, `platform/ci/ghactions_deploy.sh` (§3.6) |
| `docs/`, `evidence/` | Neither exists on `main`; `/docs/` is git-ignored. `docs/reviewer_interface.md` was created here because the spec names that path (see §4 item 3). |

## 2. Model files

- No reference model is designated (`model/candidates/README.md`). The newest candidate on `main` is
  `d8d1bac7` (2,319 reactions, 1,081 gene products, 1,882 species, 113 subsystem groups). M1 needs a
  named model to compute L5 verdicts against (question Q-A in §5).
- `model/source/iyali26.xml` has the same 2,285 reaction IDs as `model/source/iyli21.xml` (set difference 0
  in both directions). `d8d1bac7` keeps 2,270 of them, drops 15, and adds 49 (the R2286–R2290 series, the
  `R_*` additions of 2026-09/10 such as `R_DOLSYN`, `R_ADPRtn` and `R_GPIGNT`, and `TRNA_BIOMASS_*`).

## 3. Inventory

### 3.1 Reviewer

See `model/curation/evidence_pipeline/reviewer_interface.md`. Summary: the adversarial workflow is the **essentiality FN loop**
(per-case `yarrowia-essentiality-literature-reviewer`, then one `essentiality-evidence-skeptic`, then a
human gate). Its agents are Codex TOMLs on `codex/module-agent-workflow` only. Its ledger and validators are
on `main`. Its 240 stored cases are bound to old model SHAs and SD-Leu, and its prompts encode the old
contract (10 % primary cutoff, positive-only labels). It covers only essentiality FN cases, so it has
nothing that maps to R1's rule-trace adjudication or to R2's extraction red-teaming.

### 3.2 ID mapping

| Source | Content | Coverage |
| --- | --- | --- |
| `reference/locus_map/iyli21_genes_vs_S2.csv` | `model_gene, yali1_s2, yali0, category, n_reactions` | 1,058 rows: `has_YALI0` 1,039, `not_in_S2` 17, `CLIB89_only` 2 |
| `reference/locus_map/s2_metabolic_genes.csv` | `yali1, yali0, geneid` (NCBI GeneID), `metab_pathways, ec, verdict, in_model` | 8,765 rows; YALI0 7,224, GeneID 8,683, `in_model=yes` 1,041 |
| `reference/locus_map/S2_table_YALI1_YALI0_map.xlsx` | Source workbook of the S2 tables | Origin paper not recorded in the folder (unknown) |
| `essentiality/repository/curated_locus_identity_exclusions.csv` | Pairs excluded from aliasing | Used by `LocusCrosswalk` |
| `LocusCrosswalk` (`platform/scripts/gem_annotate/locus_resolver.py`) | Spelling normalization (case, optional underscore); YALI0↔YALI1 aliases **only** when the combined map is one-to-one both ways; raises `AmbiguousLocusMappingError` otherwise; never converts a YALI0 tag to the same-suffix YALI1 tag | Exactly the spec's §1.1 discipline, with `UNMAPPED` expressed as an exception |
| Model gene annotations (`d8d1bac7`) | `uniprot` 833, `ncbigene` 771, `kegg.genes` 837, `refseq` 758, `ensembl` 117 | of 1,081 gene products |

Model gene IDs by era in `d8d1bac7`: YALI1 nuclear 1,053, YALI1 mitochondrial 13, **YALI0 9**
(`YALI0D17490g`, `YALI0F04114g`, `YALI0E16192g`, `YALI0D27654g`, `YALI0E07787g`, `YALI0D12573g`,
`YALI0E12442g`, `YALI0B10406g`, `YALI0E18568g`), other 6 (`YALIUNK1/3/4`, `YALIfMp29`, `YalifMp05/06`). The
spec's canonical key is `YALI1_*`, and these 15 IDs need explicit `id_map` rows or `UNMAPPED`.

There is no local UniProt↔YALI1 table apart from the model's own annotations; a UniProt/NCBI column of
`id_map` would be built from those plus the S2 GeneID column.

### 3.3 Existing confidence scores and evidence notes in the model

- **No confidence scores.** `d8d1bac7` has no `CONFIDENCE_SCORE`-type note (the only "confidence" string is
  the prose "Confidence does not establish catalysis"), and no `pubmed`, `doi` or `eco` annotation on any
  reaction or gene. 2,314 of 2,319 reactions carry an `sboTerm`.
- **Reaction cross-references** (spec §1.1 subject links): `metanetx.reaction` 1,627, `ec-code` 1,089,
  `rhea` 828, `kegg.reaction` 610, `bigg.reaction` 1,339, `seed.reaction` 1,308, `metacyc.reaction` 703.
- **Evidence-like notes** (counts of reactions carrying the key): `PROTEIN_CLASS` 862,
  `metadata_selection_evidence_status` 204 (with `metadata_reaction_selection`, `metadata_selected_fields`,
  `metadata_selection_source_sha256`), `gpr_evidence_status` 12,
  `selection_gpr_restore_evidence_status` 9, `metadata_previous_gap_fill_direction_*` 9,
  `energy_candidate_evidence` 8. Genes: `lipoylation_evidence_tier` 3, one `function_evidence_*` block.
  Species: `metanetx_formula_charge_conflict` 437.
- **Curated evidence tables** in the repo: `model/curation/evidence/complex_i_gpr_evidence.csv` (per-gene
  structural evidence and GPR decision with source URLs), `model/curation/gapfill/*` (gap-fill candidates
  and mappings), `model/curation/essentiality/evidence/*.json` (3 dossiers).
- **Lineage:** the immediate parent is iYli21 (`model/source/iyli21.xml`). iYli21's own parents (iYali4,
  iYL619_PCP, iMK735, iNL895, …) are not recorded in the repo (unknown).

### 3.4 Prior-model and gap-fill marking

Which reactions are gap-filled vs inherited is recorded only piecemeal (`model/curation/gapfill/`, a few
`*_gap_fill_*` notes). There is no per-reaction lineage field. The spec's "mark inherited reactions"
step would start from the iYli21 ID overlap (2,270 shared IDs). That overlap shows inheritance from
iYli21, not where iYli21 got each reaction.

### 3.5 Local datasets

| Dataset | Path (`$IYALI26_RESEARCH_ROOT/…`) | Format | Use in spec |
| --- | --- | --- | --- |
| Consensus essential genes (Ramesh et al. 2023, essential in ≥ 2 of 3 screens) | `state/essentiality/repository/consensus_essential_genes.csv` (`1e887f5a…` per benchmark.md) | CSV `gene_id, source_gene_id, function, source, confidence`; 1,612 genes, positive-only | L4b `ESSENTIAL_UNDER`; **already used for development**, so not independent validation (benchmark.md) |
| Cas9 / Cas12a fitness | `state/essentiality/repository/cas9_cas12a_fitness.csv` (+ raw workbook `ba1eca8f…`) | Long CSV: `gene_id, assay, fitness_score, raw_p_value, q_value, experimental_call, source_row, source_sha256` | Continuous essentiality evidence |
| Essentiality case ledger + 240 dossiers | `state/essentiality/repository/curation_cases.csv`, `evidence/` | CSV + JSON (reviewer output) | Possible seed L0 sources (196 cited sources; paraphrased claims, no verbatim quotes) |
| Isozyme resolution, reaction chemistry curation | `state/essentiality/repository/isozyme_resolution_ledger.csv`, `reaction_chemistry_curation.csv` | CSV | Context for GPR predicates |
| W29 13C-MFA + growth | `reference/w29_worland2024/glucose_central_fluxes.csv` (71 flux rows; glucose, glycerol, oleic acid; flux, s.e., LB/UB), `growth_parameters.csv` | CSV | `FLUX_13C` claims (L4b) |
| Media / strain profile | `state/media/sd_leu.csv`, `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv`, `state/strain_profiles/po1f_sd_leu.json` | CSV / JSON | Conditions for L4b |
| External models | `reference/external_models/` (yeast-GEM incl. `metabolic_essential_genes.csv`, iMM904, iLB1027, eciYali5) | SBML / CSV | Tier D (`MODEL_INCLUSION`, never a leaf) |
| Literature PDFs | `raw/literature/file.pdf` (one file) | PDF | none yet |

Not found locally: proteomics, transcriptomics, lipidomics or phenotype-array tables. The W29 README
mentions a transcriptomic notebook in an external repository; its data is not local. The fermentor CRISPR
screen line (PO1f vs PO1f::Cas9) has reads pending (recalled from memory, not re-checked).

Discrepancy to resolve: `reference/w29_worland2024/README.md` says the work is "currently unpublished". The
benchmark calls it Worland et al. 2024 and says the link to the data repository is "inferred, pending user
confirmation".

### 3.6 Tests and CI

- `platform/tests/`: 67 files. The last recorded full run is 324 passed, 19 failed, 20 errors on 2026-10-29
  (recalled from `model/STATE.md`, not re-run in M0). `pytest` is not listed in `platform/pyproject.toml`
  dependencies; the lit tests need `pydantic-ai` in the venv (STATE.md 2026-10-22 records a missing
  `pydantic_ai`).
- GitHub CI runs **only memote** (`memote.yml` report; `diff-pr.yml` memote diff on PRs). No pytest job, so
  the spec's "gate CI on gold-set metrics" and "ECO ID test fails CI" would need a new workflow or must live
  in the local pytest run (question Q-D).
- `diff-pr.yml` on `main` has a malformed `run:|` line with a tab (a fix exists as `35c959b` on the Codex
  branch). It is unrelated to this task and recorded as a follow-up.
- No ECO OBO file is in the repo or the research workspace; the spec's ECO test needs one (cached with its
  version).

### 3.7 The literature layer vs the spec's L0–L2

`platform/gem_review/lit` (on `main`, including the prompt-caching commit `81e9d60`): Europe PMC search → cached,
rate-capped client → abstract + matching full-text paragraph chunks → pydantic-ai extraction with an output
validator that rejects any claim whose PMID or chunk was not retrieved or whose quote is not verbatim
(whitespace-normalized) in the chunk → `review.xlsx` → logic chains (quote → inference → verdict).

What carries over directly:
- quote verification by code: `check_quote` is the spec's §4 post-validation core;
- disk cache keyed by request hash with `fetched_at` (spec §3.5 "cache everything, record retrieval date");
- rate cap (≤ 5 concurrent, ≥ 0.34 s between starts);
- provider chooser and `.env` key handling;
- Anthropic prompt caching;
- per-pair `EXTRACTION_FAILED` rows (nothing silently dropped).

Differences to resolve deliberately rather than diverge silently:

| Aspect | `gem_review.lit` now | Spec | Note |
| --- | --- | --- | --- |
| Claim unit | (gene, reaction) pair, `supports_gene_reaction: bool` | (subject, predicate, value, polarity) over 8 predicates | The current unit only covers spec `GENE_ENCODES` |
| Evidence scale | `EvidenceType` 0–4 (biochemical/genetic/physiological/homology/none), chosen by the LLM | `claim_level` enum of 19 values + ECO; confidence 0–4 **derived in code** (Thiele–Palsson) | The same 0–4 range means something different; a mapping table is needed if old results are kept |
| Organism | `organism_is_yarrowia: bool` from the LLM; non-Yarrowia capped at homology | `taxon_distance` 0–5 **computed in code** from organism/strain; never trusted from the LLM | Spec §8 forbids the current approach |
| Quote length | ≤ 25 words (`MAX_QUOTE_WORDS`) | ≤ 40 words | Pick one; 40 allows more context, 25 keeps less copied text |
| Assertion strength, `is_original_finding`, `cited_support`, `rebuts`, condition, strain, method | absent | required | New extraction fields |
| Query organism names | `"Yarrowia lipolytica"` only | + *Candida*, *Saccharomycopsis*, *Endomycopsis lipolytica* | Pre-1980 literature is missed now |
| Sources | Europe PMC search + OA full text | + PubMed E-utilities, Semantic Scholar/OpenAlex, Unpaywall, Crossref/Retraction Watch, Swiss-Prot/BRENDA seeds | Needs the §10 Q3 credentials |
| Screening | none (every pair with hits goes to extraction) | `ScreenResult` agent with tier and `heterologous_only` drop | New |
| LLM call record | `usage.jsonl`: model, tokens, cache tokens; results carry `extracted_at` | model ID, **prompt hash**, parameters, timestamp → `extraction_run_id` | Prompt hash missing |
| Sampling | `anthropic_effort="medium"`, `max_tokens=16000`, temperature not set | temperature 0 | Whether the default model (`anthropic:claude-opus-5-5`) accepts `temperature=0` together with `effort` is **untested**; check before M2 |
| Storage | per-pair JSON + `review.xlsx` in the research workspace | `evidence.sqlite` with 13 tables | |
| Gold set | none | 30 hand-labelled papers, P/R gate | User to supply or approve |

Recommendation (inferred, for the user to decide): build the spec's extractor as a new module that
**reuses** `client.py`, `chunks.py`, `providers.py` and `check_quote`, with new `Claim`/`Passage` schemas.
Leave `gem_review.lit`'s current commands and outputs unchanged, so earlier reports stay reproducible.

## 4. Where the spec meets repo rules

These are not blockers for M1, but each needs a decision before M3:

1. **SBML annotations are model changes.** Reaction notes are part of the model's SHA-256
   (`platform/AGENTS.md`), and only the builder writes model XML. The spec's "SBML annotation writer" and
   "apply command" must therefore emit curated data plus an `apply_*` builder step that builds a new
   candidate, never a cobra `write_sbml_model` on a model file.
2. **Validation after apply** in this repo means a whole-model diff, the knockout screen per
   `model/expected/benchmark.md` (Task A medium + PO1f, 1/5/10/15 %, user-defined negatives) and an
   independent audit (`curation-change` skill). Memote runs in CI only on `model/candidates/legacy/model.xml`
   (`memote.ini`).
3. **Placement.** The spec's `evidence/` and `docs/` are new top-level folders. The repo convention puts
   code under `platform/`, curated decisions under `model/curation/`, reports under `model/reports/`, and
   large or generated data in the research workspace. Also, `/docs/` is **git-ignored** (`.gitignore`
   line 12, on `main` too), so `docs/reviewer_interface.md` as written is not versioned; it has to move
   (for example into this report folder or `platform/evidence/`) or the ignore rule has to change. Proposal: code `platform/evidence/`; weights YAML,
   subsystem precedence, `PROTOCOL.md` and decision JSONs under `model/curation/evidence_pipeline/`; the
   SQLite database (generated, rebuildable from cache) in
   `$IYALI26_RESEARCH_ROOT/artifacts/evidence_pipeline/`.
4. **L4b essentiality claims** come from the same list used to develop the model, so per benchmark.md they
   are a development reference. They should carry that flag and not count as independent support.
5. **"Gap-filled reactions may be needed"** (spec §2.2) matches the repo rule that deletion needs evidence of
   biological non-existence, never a score. UNSUPPORTED must never be turned into a removal proposal.
6. **Owners:** lipid-unlump (`platform/tools/lipid/`, `model/curation/lipid_unlump/`) is a colleague's
   line. A phosphatidate-hub pilot touches its reactions only read-only.

## 5. Spec §10 questions: answered from the repo vs still open

| # | Question | From the repo | Still open for the user |
| --- | --- | --- | --- |
| 1 | Where is the reviewer, how is it invoked, which outputs gate? | Located and documented (`reviewer_interface.md`): Codex FN loop, four commands, case packet → reviewer JSON → skeptic JSON → validated ledger | Confirm it is the one meant; which outputs may gate (skeptic pass/block? none?); Codex or Claude runtime; whether its prompts may change |
| 2 | Canonical strain set; YALI0↔YALI1 file | Model is "Y. lipolytica W29"; the simulation profile is PO1f. Mapping: `LocusCrosswalk` over `iyli21_genes_vs_S2.csv` + `s2_metabolic_genes.csv` + identity exclusions (one-to-one only) | Confirm distance 0 = {W29, Po1d, Po1f, Po1g, Po1h} (and whether PO1f::Cas9 and other in-house derivatives count); confirm the S2-based mapping and where the S2 table comes from |
| 3 | API access | Anthropic key in a git-ignored `platform/.env` (in another worktree; not present here; content not read); `~/.zshrc` exports an old revoked key that overrides it (from the parent brief) | NCBI E-utilities key, Semantic Scholar or OpenAlex, Unpaywall email, institutional full-text route; whether to remove the old key from `~/.zshrc` |
| 4 | Local datasets and formats | §3.5: essentiality consensus + Cas9/Cas12a fitness + reviewer dossiers; W29 13C fluxes (3 carbon sources); media/profile; external models | Any proteomics, transcriptomics, lipidomics or phenotype arrays outside the workspace? The fermentor screen counts when ready? W29 13C publication status |
| 5 | Pilot slice | Candidate pools exist: `model/curation/gem_review/tier1_lipid_subsystems.csv` (lipid tier already retrieved), central carbon in `glucose_central_fluxes.csv` (71 rows with 13C data) | Choose the ~30 reactions (default: central carbon + phosphatidate hub) |
| 6 | Existing confidence scores; parent lineage | None in the model (§3.3); immediate parent iYli21 (2,270 shared IDs in `d8d1bac7`); evidence-status notes on 204 + 12 + 9 reactions | Import those notes as `MODEL_INCLUSION`/curation evidence or ignore them? iYli21's own ancestry (which supplements to use for tier D) |
| 7 | Cost budget; models for screening vs extraction | Current default `anthropic:claude-opus-5-5`, effort medium | Budget per full-corpus run; a cheaper model for screening? |

Additional questions this survey raised:

- **Q-A Target model** for L5 verdicts in M1: the newest on `main` (`d8d1bac7`), or another? (No reference
  model is designated.)
- **Q-B Placement** (§4 item 3): accept the proposed split, or keep the spec's top-level `evidence/` and
  `docs/`?
- **Q-C Quote length** 25 vs 40 words, and whether to keep `gem_review.lit` as is beside a new extractor
  (§3.7 recommendation).
- **Q-D CI:** add a GitHub pytest job (rule core + ECO test, LLM mocked), or keep tests local as the repo
  does now?
- **Spec versioning:** copy `/Users/david/Downloads/EVIDENCE_PIPELINE_PROMPT.md` into the repo (for
  example `model/curation/evidence_pipeline/SPEC.md`, or `evidence/SPEC.md` if the spec's layout is kept)
  so that later milestones cite a versioned SHA.

## 6. What M0 did not do

- No pipeline code, no computation, no LLM call, no network access, no model or curated-data change.
- The worktree was not synced (refused, TASK.md), so it still has the 2026-10-22 tree. Survey statements
  about `main` come from `git show origin/main:…`.
- **No `model/STATE.md` entry yet.** Editing the stale `STATE.md` would block the app's sync, which refuses
  uncommitted changes, and would conflict with the 2026-10-23…29 entries on `main`. The entry below goes in
  after the sync.
- No independent audit of this report (none is required for a read-only survey; one can be run on request).

## 7. STATE.md entry to insert after sync

```markdown
## 2026-10-29: Evidence pipeline M0 survey (read-only; reviewer interface and report; not committed)

**Authorization.** User request relayed by the parent session: follow `EVIDENCE_PIPELINE_PROMPT.md`,
M0 only (read-only survey, reviewer interface, M0 report), then stop and ask the §10 questions.
Details: [TASK.md](reports/evidence_pipeline_m0_20261029/TASK.md),
[REPORT.md](reports/evidence_pipeline_m0_20261029/REPORT.md).

**Results** (read this session from `origin/main` `f212874`; no computation).
- Reviewer: the adversarial workflow is the essentiality FN loop (Codex literature reviewer + skeptic +
  human gate). Agents only on `codex/module-agent-workflow`; ledger code on `main`; 240 stored cases on old
  model SHAs and SD-Leu (85 reviewed: skeptic block 77, pass 4, downgrade 4). Not wired to anything.
- Model: no confidence scores, no pubmed/doi/eco annotations in `d8d1bac7`; 2,270 of 2,319 reaction IDs
  from iYli21; 15 non-YALI1 gene IDs.
- ID map: `LocusCrosswalk` (one-to-one YALI0↔YALI1 from the S2 tables) covers the spec's §1.1 discipline.
- `gem_review.lit` supplies quote verification, cache and providers; claim schema, organism handling
  (LLM boolean vs code-computed distance), quote length (25 vs 40) and LLM call records differ from the
  spec.
- User answers to the M0 questions (REPORT §8) recorded as deviations D1–D19 in
  [SPEC.md](curation/evidence_pipeline/SPEC.md) (spec copied verbatim, SHA-256 `70166392…`). Reference model
  for M1: `d8d1bac7`. Reviewer doc moved to `curation/evidence_pipeline/reviewer_interface.md`
  (`/docs/` is git-ignored).
- Separate CI fix: `35c959b` cherry-picked onto `main` as PR #36 (diff-pr.yml YAML).
- Open: iYli21 ancestry; whether the 10 % L4b threshold relates to the 15 % benchmark cutoff; M1 plan
  approval.

**Limits.** Survey via `git show` before the worktree could be synced. No pytest run. No pipeline code.
```

## 8. User answers (2026-10-29)

The user's answers in this session, condensed. The binding wording is in the session transcript; each
answer is recorded as a deviation in [SPEC.md](../../curation/evidence_pipeline/SPEC.md) (D-numbers).

1. **Reviewer:** the FN loop is the reviewer, used as a template. The skeptic becomes R3, generalized; R1 and
   R2 are new agents with the typed contracts. Agent definitions go onto `main` as versioned prompt files;
   legacy prompts are kept with their hash. New prompts emit polarity and contain no scoring rules or
   cutoffs. Nothing gates until a calibration on 85 reviewed cases (8 passed or downgraded + 12 random
   blocked, human-labelled) reports the skeptic's precision. The 240 cases are imported as `legacy_case`
   sources with `legacy_scoring=True`, never merged into live assertions, and re-extracted when they enter a
   pilot. Runtime: the extractor and R1/R2 on the `gem_review.lit` provider; R3 can run on another model
   family. (D7, D10–D12)
2. **Strains / IDs:** distance 0 = {W29, Po1d, Po1f, Po1g, Po1h, PO1f::Cas9}, with strain recorded exactly.
   `Condition` gains `medium` and `marker_complementation`. The S2 crosswalk gets a `cardinality` field, and
   anything not 1:1 is UNRESOLVED. The 15 non-YALI1 IDs get dispositions in the M1 report. (D4–D6)
3. **Access:** NCBI, KEGG, BiGG, MetaNetX, Europe PMC OA and Unpaywall; M2 reports full-text coverage
   first. (D13)
4. **Data:** nothing outside the workspace. M2 lists Tier E candidates without ingesting them. W29 13C is
   usable for L4b as `peer_reviewed=False, unpublished=True`; internal-data and literature evidence are kept
   apart in every report. (D9, D14)
5. **Pilot:** 12 reviewed essentiality cases (mix of blocked, passed, downgraded), 12 phosphatidate-hub
   reactions, 6 central-carbon controls; the list is proposed in the M1 report.
6. **Model notes** are imported as `prior_model` / `MODEL_INCLUSION` sources (`source_id` = model commit);
   a `parent_model` chain plus `inherited_from`. iYli21's own ancestry was left as a placeholder in the
   answer and is still to be supplied. (D8)
7. **Budget:** $20 per full-corpus run, `claude-fable-5-1`. (D19)
8. **Extras:** reference model `d8d1bac7`, with model as an engine input and path/commit on every verdict
   (D16). Layout as proposed, with a deterministic SQLite rebuild and decision records plus per-commit
   snapshots in the repo (D1, D2). Quotes 25 words plus character offsets for a ±300-character window (D3).
   pytest CI with the LLM mocked gates merges; memote stays separate (D17). Spec copied with a Deviations
   section. L4b essentiality threshold in the weights YAML at 10 %, sensitivity at 5/10/20 % (D18).
9. **diff-pr.yml:** cherry-pick `35c959b` as its own PR. Done: PR #36 (`claude/fix-diff-pr-yaml`,
   commit `cc9c119`).
10. **Extractor:** new module reusing the quote check, cache and providers; LLM picks `claim_level` with a
    justification and extracts organism and strain; code computes score and `taxon_distance`; prompt hash
    and temperature 0 are recorded; four species names; legacy 0–4 scores are never mixed in. (D15)
