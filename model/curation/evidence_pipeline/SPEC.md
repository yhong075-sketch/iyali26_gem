# iYali26 evidence pipeline: claims, contradictions, layered chain of evidence, reviewer integration

## 0. Read this first

You are working in the `iyali26_gem` repository (COBRApy genome-scale model of *Yarrowia lipolytica*). Build a pipeline that, for every reaction in the model:

1. collects the literature and datasets that can bear on it (Section 3),
2. extracts typed, scoped claims with a Pydantic AI layer (Section 4),
3. detects contradictions with deterministic code (Sections 1, 5),
4. assembles a layered chain of evidence from sources up to per-reaction verdicts (Section 2),
5. routes contradictions and contradicted verdicts to the existing AI reviewer in this repo (Section 6),
6. records curation decisions with full provenance in a sidecar database and in SBML annotations (Section 7).

Working rules, non-negotiable:

- Deterministic code decides comparability, contradiction, scoring and verdicts. LLMs screen papers, extract claims, and draft adjudication rationales. No LLM output is ever used as a rule inside the rule engine.
- V1 and V2 make no edits to the model. They produce patch proposals and a triage queue. Applying a patch is a separate, explicit command run after a human approves the decision record.
- Every LLM call is recorded with model id, prompt hash, parameters and timestamp. Every extracted quote is verified by exact string match (whitespace-normalized) against the source text before the claim is accepted.
- Work in milestones (Section 9). Stop and report after each. Do not skip ahead.
- Before writing any code: inspect the repo (existing AI reviewer, model files, ID mapping files, package layout, tests, CLAUDE.md, CI). Reuse what exists. Then ask the questions in Section 10 that you cannot answer from the repo.

Target stack: Python, Pydantic v2, `pydantic-ai` for the typed LLM layer, COBRApy, pytest. Follow the repo's existing conventions where they differ from this document; flag the difference rather than silently diverging.

---

## 1. Core definitions

Contradictions exist between **claims**, not between articles. Everything below is designed so that "contradiction" is a computed property of two normalized claims, not a judgment.

### 1.1 Subject normalization

- **Reactions**: model reaction IDs are the subject keys. Cross-reference each to Rhea, MetaNetX, KEGG, EC where the repo already has mappings. A claim about a reaction expressed only as an EC number or enzyme name is linked to model reactions through those mappings; ambiguous links (one EC → several model reactions) are recorded as `candidate_subjects` and resolved by compartment and substrate, else left `UNRESOLVED` and excluded from pair comparison (logged, never dropped).
- **Genes**: canonical gene IDs are the W29 (`YALI1_*`) identifiers. Build an `id_map` table: `YALI0` (CLIB122) ↔ `YALI1` (W29) ↔ UniProt ↔ NCBI ↔ gene symbols. Every claim's gene IDs are normalized before any comparison. Unmapped IDs → `UNMAPPED`, excluded from pairing, logged. A large fraction of apparent GPR contradictions are ID-era artifacts; this step is mandatory, not optional.
- **Organism scope**: every claim carries `taxon_distance` relative to *Y. lipolytica* W29:
  - 0 = canonical strain set (W29 and its derivatives Po1d/Po1f/Po1g/Po1h; confirm with the user),
  - 1 = other *Y. lipolytica* strains (CLIB122, H222, etc.),
  - 2 = other *Yarrowia* species,
  - 3 = other Saccharomycotina (*S. cerevisiae*, *Candida*, *Komagataella*, *Kluyveromyces*),
  - 4 = other fungi,
  - 5 = non-fungal.
  `strain` is nullable; a null strain never defaults to W29.

### 1.2 Predicates

```python
class Predicate(str, Enum):
    REACTION_EXISTS   = "reaction_exists"    # organism has the capability to catalyze R; condition-agnostic
    FLUX_ACTIVE_UNDER = "flux_active_under"  # R carries flux under condition C; condition-scoped
    GENE_ENCODES      = "gene_encodes"       # gene G encodes an enzyme catalyzing R
    LOCALIZED_IN      = "localized_in"       # enzyme/reaction occurs in compartment set
    DIRECTIONALITY    = "directionality"     # forward | reverse | reversible, with a level flag
    COFACTOR          = "cofactor"           # cofactor set per role (e.g. NAD vs NADP)
    SUBSYSTEM_MEMBER  = "subsystem_member"   # pathway IDs in a named ontology
    ESSENTIAL_UNDER   = "essential_under"    # gene/reaction essential under condition C (datasets)
```

Value types: `REACTION_EXISTS`, `FLUX_ACTIVE_UNDER`, `ESSENTIAL_UNDER` → bool via `polarity`; `GENE_ENCODES` → `(gene_id, reaction_id, exclusive: bool)`; `LOCALIZED_IN` → set of model compartments; `DIRECTIONALITY` → enum + `level ∈ {biochemical, physiological}`; `COFACTOR` → dict role → set; `SUBSYSTEM_MEMBER` → set of pathway IDs + ontology name.

Derivation rule, not comparison rule: a positive `FLUX_ACTIVE_UNDER` claim implies a positive `REACTION_EXISTS` claim at the same level. Generate the derived claim with `derived_from` set; do not compare flux claims against existence claims directly.

### 1.3 `claim_level`

The single most important field. It records *what kind of observation* the claim rests on. Levels are categorical; directness is defined per predicate (table below), not as one global order.

| `claim_level` | What it is | ECO | Admissible for |
|---|---|---|---|
| `STRUCTURE` | 3D structure of the native protein/complex | 0000314 | EXISTS, GENE_ENCODES (subunits), COFACTOR |
| `ACTIVITY_PURIFIED` | purified native enzyme, activity measured | 0000314 | EXISTS, GENE_ENCODES, COFACTOR, DIRECTIONALITY(biochem) |
| `ACTIVITY_HETEROLOGOUS` | native gene expressed in a host, activity measured | 0000314 | GENE_ENCODES, COFACTOR, DIRECTIONALITY(biochem) |
| `FLUX_13C` | in vivo flux from 13C-MFA | 0000314 | EXISTS, FLUX_ACTIVE_UNDER, DIRECTIONALITY(physiol) |
| `ACTIVITY_CRUDE` | crude extract / whole-cell assay | 0000314 | EXISTS, COFACTOR; **not** GENE_ENCODES |
| `LOCALIZATION_DIRECT` | GFP fusion, fractionation, immuno-EM, validated targeting signal | 0000314 | LOCALIZED_IN |
| `PHENOTYPE_COMPLEMENTATION` | native gene rescues a defined mutant | 0000315 | GENE_ENCODES |
| `PHENOTYPE_KO` | deletion/disruption phenotype | 0000315 | GENE_ENCODES, EXISTS, ESSENTIAL_UNDER |
| `PHENOTYPE_OE` | overexpression phenotype of a **native** gene | 0000315 | GENE_ENCODES (weak) |
| `EXPRESSION_PROTEIN` | proteomics / western detection | 0000270 | EXISTS (weak); never GENE_ENCODES, never LOCALIZED_IN unless fractionated |
| `EXPRESSION_TRANSCRIPT` | RNA-seq / qPCR detection | 0000270 | EXISTS (weak) |
| `PHYSIOLOGY` | growth, secretion, substrate range, auxotrophy; no gene named | 0000314 | EXISTS (pathway-level), FLUX_ACTIVE_UNDER |
| `ORTHOLOGY` | characterized ortholog in another organism; carries `transfer_from_level` and `taxon_distance` | 0000266 | any, down-weighted by distance |
| `HOMOLOGY_ANNOTATION` | sequence-similarity or automatic annotation | 0000250 / 0000501 | GENE_ENCODES (weak), EXISTS (weak) |
| `GENOMIC_PRESENCE` | gene present in the assembly | 0000501 | EXISTS (weakest) |
| `AUTHOR_STATEMENT_TRACEABLE` | statement in a review/intro with a citation | 0000304 | any, never a leaf |
| `AUTHOR_STATEMENT_NONTRACEABLE` | statement without citation | 0000303 | any, never a leaf |
| `MODEL_INCLUSION` | present in a prior GEM | 0000305 | any, never a leaf |
| `MODEL_VALIDATION` | agreement/disagreement between a model prediction and a dataset (produced by Layer 4b, Section 2) | 0000305 (or a more specific computational ECO term, verify) | EXISTS, GENE_ENCODES, ESSENTIAL_UNDER |

Store ECO IDs as strings and validate them against the ECO OBO file in a test, so a wrong ID fails CI rather than propagating.

Per-predicate directness (most → least direct; levels not listed are inadmissible for that predicate):

- `REACTION_EXISTS`: FLUX_13C ≈ ACTIVITY_PURIFIED ≈ STRUCTURE > ACTIVITY_CRUDE > PHENOTYPE_KO > PHYSIOLOGY > EXPRESSION_PROTEIN > EXPRESSION_TRANSCRIPT > ORTHOLOGY > HOMOLOGY_ANNOTATION > GENOMIC_PRESENCE > AUTHOR_STATEMENT_* > MODEL_INCLUSION
- `GENE_ENCODES`: ACTIVITY_HETEROLOGOUS ≈ ACTIVITY_PURIFIED ≈ STRUCTURE > PHENOTYPE_COMPLEMENTATION > PHENOTYPE_KO > PHENOTYPE_OE > ORTHOLOGY > HOMOLOGY_ANNOTATION > AUTHOR_STATEMENT_* > MODEL_INCLUSION
- `LOCALIZED_IN`: LOCALIZATION_DIRECT > ACTIVITY_PURIFIED (only with fractionation) > ORTHOLOGY > HOMOLOGY_ANNOTATION (targeting-signal prediction) > MODEL_INCLUSION
- `COFACTOR`, `DIRECTIONALITY(biochem)`: ACTIVITY_PURIFIED ≈ ACTIVITY_HETEROLOGOUS > ACTIVITY_CRUDE > ORTHOLOGY > HOMOLOGY_ANNOTATION
- `FLUX_ACTIVE_UNDER`: FLUX_13C > PHYSIOLOGY > EXPRESSION_* (expression is not flux; weakest)
- `SUBSYSTEM_MEMBER`: ontology-defined; no empirical directness

### 1.4 Claim schema

```python
class Condition(BaseModel):
    carbon_source: str | None = None
    nitrogen: str | None = None          # e.g. "limited", "replete"
    oxygen: str | None = None            # "aerobic", "hypoxic", "anaerobic"
    growth_phase: str | None = None
    adaptation: str | None = None        # e.g. "after serial passage on xylose"
    free_text: str | None = None

class Claim(BaseModel):
    claim_id: str
    source_id: str                        # Layer 0 node
    passage_id: str                       # Layer 1 node; quote verified
    subject: SubjectRef                   # normalized reaction / (gene, reaction) / (reaction, compartment)
    predicate: Predicate
    value: ClaimValue                     # discriminated union keyed on predicate
    polarity: Literal["positive", "negative"]
    assertion_strength: Literal["assertive", "hedged", "speculative"]
    claim_level: ClaimLevel
    transfer_from_level: ClaimLevel | None = None   # ORTHOLOGY only
    organism: str
    strain: str | None
    taxon_distance: int
    compartment: str | None
    condition: Condition | None
    strain_specific: bool = False         # author explicitly frames the result as strain-specific
    method: str
    is_original_finding: bool             # False for restated results
    cited_support: list[str] = []         # citation keys the paper gives for this claim
    rebuts: list[str] = []                # citation keys the paper explicitly argues against
    eco: str
    extraction_run_id: str
    authors: list[str]
    corresponding_author: str | None
    year: int
```

`speculative` claims are stored but never participate in contradiction detection or scoring. `hedged` claims produce weak contradictions only.

### 1.5 Comparability

Two claims `a`, `b` are **comparable** iff all of:

1. `normalize(a.subject) == normalize(b.subject)` and `a.predicate == b.predicate`.
2. Both have `taxon_distance <= 1`. Claims at distance ≥ 2 are transfers: they enter Layer 4 as down-weighted support and are never paired, with each other or with *Y. lipolytica* claims. An ortholog's cytosolic localization does not contradict a peroxisomal localization in *Y. lipolytica*; it is a species difference.
3. Strain: compatible unless either claim has `strain_specific=True` and strains differ.
4. Condition: for `FLUX_ACTIVE_UNDER` and `ESSENTIAL_UNDER`, conditions must overlap (every field that is specified in both is equal). For the other predicates, condition is recorded but not a comparability criterion; a condition-scoped negative flux claim does not contradict a capability-level existence claim.
5. Level band: define two bands, **functional** = {STRUCTURE, ACTIVITY_*, FLUX_13C, LOCALIZATION_DIRECT, PHENOTYPE_*, PHYSIOLOGY} and **parts** = {EXPRESSION_*, GENOMIC_PRESENCE, HOMOLOGY_ANNOTATION}. Opposite-polarity claims from different bands are not contradictions; they are `ANOMALY_CROSS_LEVEL` (see 1.7) if the lower-band claim is negative and the higher-band claim positive, else `SCOPED_APART(level_mismatch)`. Same-band pairs proceed to the incompatibility check.
6. Neither claim is `speculative`.

### 1.6 Incompatibility, per predicate

- `REACTION_EXISTS`, `FLUX_ACTIVE_UNDER`, `ESSENTIAL_UNDER`: incompatible iff polarities differ.
- `GENE_ENCODES`:
  - same gene, polarities differ: incompatible **only at the same level**. A negative `PHENOTYPE_KO` ("deletion has no effect") never contradicts a positive `ACTIVITY_*` claim; it is explained by isozymes and classified `SCOPED_APART(redundancy_possible)`.
  - different genes: incompatible only if at least one claim has `exclusive=True`, or one claim is an explicit negative about the other gene.
- `LOCALIZED_IN`: incompatible iff compartment sets are disjoint and both claims are at `LOCALIZATION_DIRECT` (or purified-with-fractionation). Predicted localization never contradicts direct localization. A protein in two compartments is not a contradiction.
- `DIRECTIONALITY`: incompatible iff both claims are at the same `level` and values are incompatible (`forward` vs `reverse`; `reversible` vs an irreversible value only if both are `assertive`). A physiological "runs forward in vivo" is `SCOPED_APART` from a biochemical "reversible in vitro".
- `COFACTOR`: incompatible iff sets for the same role are disjoint and both claims are at an activity level. "Prefers NADP" vs "uses NAD and NADP" is a subset, not a contradiction.
- `SUBSYSTEM_MEMBER`: never incompatible. Take the union; choose a primary subsystem by a fixed precedence list committed in the repo.

### 1.7 Pair classification

```python
class PairClass(str, Enum):
    AGREE_INDEPENDENT    = "agree_independent"
    AGREE_DERIVATIVE     = "agree_derivative"     # b cites a for this claim, or shares a root source
    SCOPED_APART         = "scoped_apart"         # comparable subject, incompatible scope; carries a reason
    SUPERSEDED           = "superseded"           # later claim explicitly rebuts earlier; later wins by default, logged
    ANOMALY_CROSS_LEVEL  = "anomaly_cross_level"  # e.g. no gene annotated, activity detected: investigate
    CONTRADICTION_STRONG = "contradiction_strong" # both assertive
    CONTRADICTION_WEAK   = "contradiction_weak"   # at least one hedged
```

`SILENT` (one side has no claim) is not a pair class; it is computed at Layer 4 as one-sided evidence. Absence of a negative claim is not evidence of absence.

Derivative detection: `b` is derivative of `a` if `a.source_id ∈ b.cited_support` for that claim, or `b.is_original_finding == False` and `a` is the earliest original finding for that (subject, predicate, polarity). Independence class: two claims share an independence class if they share a corresponding author, or ≥ 2 authors, within 5 years. Shared-class agreement counts as agreement but not as independent support.

Supersession: `b` supersedes `a` if `a.source_id ∈ b.rebuts` and the pair would otherwise be a contradiction. Default resolution favors `b`; the pair is still shown to the reviewer.

---

## 2. Layered evidence graph (chain of evidence)

Model this as a layered DAG. The neural-network analogy is deliberate and should be used only as far as it stays auditable: a forward pass propagates evidence up through fixed, tabulated weights; a backward pass traces any verdict down to its leaves; attribution is leave-one-out sensitivity. Weights are frozen in V1/V2 (one committed YAML table). V3 may fit them against validation data, but only if the fitted table stays interpretable, is committed with its evaluation, and the fitting is never a black box.

| Layer | Node | Key fields | Produced by |
|---|---|---|---|
| L0 | `Source` | id, kind ∈ {article, preprint, thesis, patent, dataset, database_entry, prior_model}, DOI/PMID, year, tier (Section 3), peer_reviewed, retraction_status, organisms | collector |
| L1 | `Passage` | source_id, verbatim quote (≤ 40 words), section, char offsets if available, hash, `verified: bool` | extractor + verifier |
| L2 | `Claim` | Section 1.4 | extractor (LLM) or deterministic dataset converter |
| L3 | `PairRelation` | claim_a, claim_b, PairClass, reason, independence_class | contradiction engine |
| L4 | `Assertion` | one per (subject, predicate): admissible claims by polarity, best directness per polarity, `n_independent` per polarity, conflict_status, resolution trace, derived confidence 0–4, load-bearing leaves | aggregator |
| L4b | `ValidationEvidence` | model prediction vs dataset (essentiality, growth on substrates, 13C flux) → `MODEL_VALIDATION` claims feeding L4 | validator (V3) |
| L5 | `Verdict` | one per model element (reaction existence, GPR, compartment, bounds, subsystem): SUPPORTED, CONTRADICTED, UNSUPPORTED, CONFLICTING, PARTIAL; proposed action; priority | verdict engine |
| L6 | `Decision` | verdict_id, action, rule trace, reviewer verdict id, curator, date, commit, patch id | humans + reviewer gate |

Edges always point from evidence to conclusion (L0 → L1 → L2 → L3/L4 → L5 → L6). Every L5 verdict must be traceable to L0 leaves; `trace(verdict)` is a required function and a required test.

### 2.1 Forward pass: Layer 4 aggregation

For each (subject, predicate):

1. Keep admissible claims only (Section 1.3 table), excluding `speculative`.
2. Group by polarity. Compute `best_directness[polarity]` (per-predicate ranking), `n_independent[polarity]` (distinct independence classes, derivative claims excluded), `max_taxon_distance_used`.
3. `conflict_status` from L3: NONE | CONTRADICTION_STRONG | CONTRADICTION_WEAK | ANOMALY.
4. Resolution trace, applied in order, each step recorded as `(rule_id, outcome)`:
   - R1 scope: already handled by comparability; record how many pairs were removed as `SCOPED_APART`.
   - R2 directness/level: the side with the more direct best claim wins.
   - R3 taxon distance: on a tie, the side with lower distance wins.
   - R4 independence: on a tie, the side with more independent classes wins.
   - R5 explicit rebuttal: `SUPERSEDED` pairs resolve toward the rebutting claim regardless of R2–R4 if the rebuttal is `assertive`; log it.
   - R6 recency: tiebreaker only.
   - R7 network consistency (dead-end creation, thermodynamically infeasible cycle, orphan metabolite) is an **advisory flag** attached to the assertion. It never decides. Letting the model adjudicate its own content is circular.
   - If R2–R6 leave a tie: `resolved = False`, status stays CONFLICTING.
5. Derived confidence (Thiele & Palsson scale), organism-specific (`taxon_distance <= 1`) unless stated:
   - 4: best positive claim ∈ {ACTIVITY_PURIFIED, ACTIVITY_HETEROLOGOUS, STRUCTURE, FLUX_13C (existence only), LOCALIZATION_DIRECT (localization only)} and conflict resolved or NONE.
   - 3: best positive ∈ {PHENOTYPE_COMPLEMENTATION, PHENOTYPE_KO, PHENOTYPE_OE}.
   - 2: best positive ∈ {EXPRESSION_*, PHYSIOLOGY, ACTIVITY_CRUDE for GPR-free existence}, or `ORTHOLOGY` with `taxon_distance <= 3` and `transfer_from_level` at an activity level.
   - 1: only HOMOLOGY_ANNOTATION, GENOMIC_PRESENCE, AUTHOR_STATEMENT_*, MODEL_INCLUSION.
   - 0: no admissible claims.
   - Unresolved contradiction caps the score at 2 and sets `flag=CONFLICTING`. A resolved contradiction scores from the winning side with `resolved_by=<rule_id>`.
   - A score is reported with `n_independent` alongside; do not fold independence into the integer.

### 2.2 Layer 5 verdicts

Compare the model's current state with the L4 assertion:

- reaction present in model, assertion positive, score ≥ 3 → SUPPORTED
- reaction present, assertion negative and resolved → CONTRADICTED, action `propose_remove` or `propose_constrain` (if negative is condition-scoped)
- reaction present, score ≤ 1 → UNSUPPORTED, action `flag_for_evidence` (not removal; gap-filled reactions may be needed, record that explicitly)
- GPR: compare gene sets; missing supported gene → PARTIAL(`add_gene`); gene with negative resolved assertion → PARTIAL(`remove_gene`); exclusive claim vs model isozymes → CONFLICTING
- compartment: resolved disjoint → CONTRADICTED(`relocate`)
- directionality/cofactor: resolved incompatible → CONTRADICTED(`change_bounds` / `change_stoichiometry`)
- subsystem: PARTIAL(`reassign`) by precedence; never CONTRADICTED
- unresolved L4 → CONFLICTING

Priority for the human queue = model impact (essential in simulations, flux centrality, involvement in the lipid/phosphatidate network) × uncertainty (CONFLICTING > CONTRADICTED > PARTIAL > UNSUPPORTED).

### 2.3 Backward pass

- `trace(verdict_id)` → ordered list of nodes down to L0 leaves.
- `load_bearing(verdict_id)` → for each L2 leaf, recompute the assertion and verdict without it; return the leaves whose removal changes the score or verdict. Store on the L4 node. These are the claims the reviewer must re-read (Section 6, role R2).
- `chain_quality(assertion)` → does any chain terminate at an organism-specific functional-band leaf? If every chain bottoms out in HOMOLOGY_ANNOTATION, ORTHOLOGY or MODEL_INCLUSION, flag `homology_only=True` regardless of how many articles mention the reaction.

---

## 3. What to collect

Collect by evidence yield: a source earns full-text extraction if it can place a node in the chain that could change a verdict. Everything else is indexed (snowballing, ID mapping) but never becomes an evidence node. Build the corpus organism-first: screen the whole *Y. lipolytica* literature once, extract, index by gene/EC/reaction; per-reaction work becomes a lookup, and only gaps trigger targeted ortholog searches.

### 3.1 Tiers (`SourceTier`)

- **A — organism-specific primary** (leaf nodes): enzyme characterization (purified or heterologous native enzymes), gene function (deletion, overexpression of native genes, complementation), localization studies, 13C-MFA, omics datasets (proteome, transcriptome, metabolome, lipidome), physiology (substrate range, secreted products, auxotrophies, fermentation capacity; the main source of negative evidence), structural biology of native proteins/complexes, genome and annotation papers (CLIB122, W29, mitochondrial genome, re-annotation releases; leaf only at homology level, mandatory for `id_map`).
- **B — transferred primary**: ortholog characterization in other organisms, with `taxon_distance` recorded. Entry point for Saccharomycotina orthologs is the evidence-coded literature in SGD/CGD, not a fresh screen of the *S. cerevisiae* literature. Reaction chemistry (what an EC does) may come from any organism; organism-specific existence/GPR only with the distance penalty.
- **C — indexes**: UniProt/Swiss-Prot (reviewed entries carry ECO-tagged PubMed references), BRENDA organism entries, KEGG `yli`, MetaCyc/BioCyc, Rhea, GO annotations. Use their reference lists as high-precision candidates. A database assertion is never a leaf unless the cited paper is fetched and extracted.
- **D — coverage**: reviews of *Y. lipolytica* metabolism; prior GEMs and their supplements (iNL895, iYL619_PCP, iMK735, iYali4, iYLI647, iYli21) and yeast-GEM curation notes. Use for backward snowballing and for marking which reactions were gap-filled vs inherited. Never a leaf.
- **E — validation and scoping**: strain lineage papers (W29 → Po1 series; CLIB122; H222), essentiality datasets (transposon, CRISPR, including in-house), phenotype arrays, flux tables. Stored as L0 `dataset` sources; their rows become L2 claims through deterministic converters, not through the LLM. Used by L4b in V3.

### 3.2 The bulk problem

Most *Y. lipolytica* papers are titer-improvement studies. Overexpressing a gene and seeing more product establishes that the construct works, not that the native reaction exists or that the GPR is right. Collect an engineering paper only if full text contains at least one of: a deletion/disruption phenotype of a native gene, an activity assay on a native enzyme, a localization result, or flux measurements. The screening agent classifies `heterologous_only` vs `native_gene_characterized` and the former is dropped before extraction.

### 3.3 Exclusions and tags

- Patents, theses, conference abstracts: `non_peer_reviewed`, never a sole leaf.
- Retracted or suspected paper-mill papers: check Crossref/Retraction Watch at ingest; exclude retracted, tag suspected.
- Organism mis-scoping: other *Yarrowia* species, strain-mislabeled work; keep `strain` nullable.
- Abstract-only records when full text is retrievable: do not extract from abstracts; they drop negation, scope and method.

### 3.4 Historical names

Query with all of: *Yarrowia lipolytica*, *Candida lipolytica*, *Saccharomycopsis lipolytica*, *Endomycopsis lipolytica*. Pre-1980 enzymology and physiology (citric acid production, alkane assimilation) is published under the old names and is often the most direct evidence available.

### 3.5 Mechanics

1. Corpus pull: PubMed E-utilities and Europe PMC (all name variants); citation graph from Semantic Scholar or OpenAlex; full text from Europe PMC OA, Unpaywall, and whatever institutional route the user specifies. Cache everything on disk; respect rate limits; record retrieval date.
2. Seed from indexes: harvest PubMed references attached to every Swiss-Prot *Y. lipolytica* entry, every BRENDA *Y. lipolytica* enzyme entry, and every prior-GEM supplement.
3. Two-stage screen: title/abstract → `ScreenResult` (Section 4); full-text extraction only for tiers A/B and filtered engineering papers. Snowball backward from reviews and prior models, forward from the genome papers and key enzyme papers.
4. Ortholog pass only for gaps: reactions with no Tier A node fetch the SGD/CGD evidence-coded references for the ortholog; record the orthology call (tool, parameters) separately from the claim.
5. Stop criterion per reaction: at least one Tier A functional-band node, or Tier A exhausted and the best Tier B node recorded with its distance. Log both outcomes.
6. Write the inclusion/exclusion criteria to `evidence/PROTOCOL.md` before running them; this is the pre-registered protocol for the eventual paper.

---

## 4. Extraction layer (`pydantic-ai`)

Two agents, both with structured outputs, temperature 0, pinned model id, prompt text stored in the repo and hashed into `extraction_run_id`.

**Screening agent** — input: title, abstract, year, journal, availability of full text. Output:

```python
class ScreenResult(BaseModel):
    paper_type: PaperType     # enzyme_characterization | gene_function | localization | flux | omics | physiology | structure | genome | review | prior_model | engineering_heterologous_only | engineering_native_characterized | other
    tier: SourceTier
    organisms: list[str]
    has_native_gene_characterization: bool
    has_negative_result: bool
    subject_mentions: list[str]  # EC numbers, gene symbols/IDs, metabolites
    exclude_reason: str | None
```

**Extraction agent** — input: full text (sectioned). Output: `list[Claim]` plus `Passage` for each claim. Requirements: verbatim quote ≤ 40 words; `polarity` and `assertion_strength` mandatory; `cited_support` and `rebuts` filled from in-text citations; `is_original_finding` set; `claim_level` chosen from the enum with a one-line `level_justification` field that is stored but not used downstream.

**Post-validation (code, not LLM)**: quote exact-match against the source text (whitespace-normalized) → else `rejected_claims(reason=quote_not_found)`; subject normalization via `id_map` and reaction cross-references → else `UNRESOLVED`; ECO consistent with `claim_level` → else reject; `taxon_distance` computed from organism/strain, not taken from the LLM. Nothing is silently dropped; every rejection has a row.

**Gold set**: 30 hand-labeled papers (the user supplies or approves the list), with labeled claims including negatives and hedges. Report precision/recall for claim extraction, polarity, `claim_level`, and `cited_support`. Gate CI on a threshold the user sets. Re-run on every prompt or model change.

---

## 5. Contradiction engine (pure Python)

Pairwise comparison within (subject, predicate) buckets; buckets are small, O(n²) is fine. Output: L3 `PairRelation` rows. No LLM anywhere in this module.

Required unit tests, each with the expected `PairClass`:

1. "does not ferment ethanol" (PHYSIOLOGY, negative, FLUX_ACTIVE_UNDER aerobic) vs "alcohol dehydrogenase activity in extracts" (ACTIVITY_CRUDE, positive, REACTION_EXISTS) → not comparable (different predicates); no contradiction.
2. "cannot grow on xylose" (PHYSIOLOGY, negative, condition unadapted) vs "grows on xylose after adaptation" (PHYSIOLOGY, positive, condition adapted) → SCOPED_APART(condition).
3. Same gene, YALI0 ID in one claim, YALI1 ID in the other, same assertion → AGREE_* after normalization; the test must fail if normalization is skipped.
4. Gene A KO "no effect on activity" (PHENOTYPE_KO, negative) vs gene A heterologous activity (ACTIVITY_HETEROLOGOUS, positive) → SCOPED_APART(redundancy_possible).
5. Gene A encodes (positive) vs gene B encodes (positive), neither exclusive → not a contradiction; both support the assertion.
6. Peroxisomal (LOCALIZATION_DIRECT) vs mitochondrial (LOCALIZATION_DIRECT), both assertive, same species → CONTRADICTION_STRONG.
7. Peroxisomal (LOCALIZATION_DIRECT) vs cytosolic (HOMOLOGY_ANNOTATION, targeting-signal prediction) → SCOPED_APART(level_mismatch).
8. *S. cerevisiae* ortholog claim vs *Y. lipolytica* claim → not comparable (taxon_distance 3).
9. Review restating paper X's finding vs paper X → AGREE_DERIVATIVE.
10. Later paper with `rebuts=[earlier]`, opposite polarity, assertive → SUPERSEDED.
11. No gene annotated for EC (HOMOLOGY_ANNOTATION, negative) vs activity measured (ACTIVITY_CRUDE, positive) → ANOMALY_CROSS_LEVEL.
12. Any pair where one claim is `speculative` → no pair relation.

---

## 6. Integration with the existing AI reviewer

The repo already has an AI reviewer (an adversarial, multi-agent review workflow). Do not rewrite it. Do not change its prompts without asking. First, locate it, read it, and write `docs/reviewer_interface.md` describing: entry point, how it is invoked, input and output formats, what it is good at and what it currently is not asked to do. Stop and confirm that document with the user before wiring anything.

Then integrate it in three roles, each with a typed contract:

**R1 — Adjudicator.** Every L3 `CONTRADICTION_STRONG`, `CONTRADICTION_WEAK` and `ANOMALY_CROSS_LEVEL`, and every L5 `CONFLICTING` or `CONTRADICTED` verdict, is sent as a packet:

```python
class ConflictPacket(BaseModel):
    subject: SubjectRef
    predicate: Predicate
    claims: list[Claim]            # all admissible claims, both polarities
    passages: list[Passage]
    pair_relations: list[PairRelation]
    resolution_trace: list[tuple[str, str]]   # rule engine's own reasoning
    proposed_verdict: Verdict | None
    model_state: ModelElementState            # what iYali26 currently says
    network_flags: list[str]                  # R7 advisory flags
```

Expected return:

```python
class ReviewVerdict(BaseModel):
    agrees_with_rule_engine: bool
    verdict: Literal["side_positive", "side_negative", "unresolved", "scope_error", "extraction_error"]
    rationale: str
    disputed_rule: str | None            # which rule it thinks misfired
    disputed_claims: list[str]           # claim_ids it believes are mis-extracted
    requested_evidence: list[str]        # what would settle it
    confidence: Literal["low", "medium", "high"]
```

**R2 — Red team of extraction.** For every load-bearing claim in a CONTRADICTED or CONFLICTING verdict, and for a random sample of all claims, the reviewer re-reads the passage and challenges `polarity`, `assertion_strength`, `claim_level` and scope. Output is a `ClaimChallenge`; an upheld challenge triggers re-extraction, never a direct edit of the claim.

**R3 — Decision gate.** Every proposed L6 action (remove, relocate, change GPR, change bounds, reassign subsystem) is reviewed adversarially against its evidence packet before it enters the human queue: the reviewer must produce the strongest case for and against applying it, and a recommendation. Output `DecisionReview`.

Rules for the integration:

- The reviewer cannot modify L2–L4 data. It can request re-extraction, propose a resolution, or disagree with the rule engine.
- Every disagreement is logged in `reviewer_disagreements` with the rule id. Repeated disagreement on one rule is the main signal that the rule needs changing; surface counts per rule to the user. Do not auto-change rules.
- Human queue: ordered by Section 2.2 priority. Every item shows the rule-engine trace and the reviewer's verdict side by side. The curator's decision, with a one-line rationale, becomes the L6 record.
- Applying a decision produces a patch proposal (JSON: element, action, before, after, evidence ids, decision id) and the SBML annotation payload. A separate explicit command applies approved patches and runs the repo's existing validation (Memote, mass/charge balance, growth tests) before and after.

---

## 7. Storage, provenance, annotations

- Sidecar database `evidence/evidence.sqlite` (DuckDB acceptable if the repo prefers it) with tables: `sources`, `passages`, `claims`, `rejected_claims`, `pair_relations`, `assertions`, `verdicts`, `decisions`, `reviewer_runs`, `reviewer_disagreements`, `extraction_runs`, `id_map`, `weights_version`. Provide a `trace` CLI that prints the chain for any reaction, GPR or decision.
- Weights, directness rankings, confidence rules and subsystem precedence live in one committed YAML file with a version string stored on every assertion.
- SBML annotations via COBRApy: `pubmed`, `doi`, `eco`, `sbo` on reactions and genes; `CONFIDENCE_SCORE`, `N_INDEPENDENT`, `EVIDENCE_ID` in notes. COBRApy's annotation dict emits a single biological qualifier; if `bqbiol:isDescribedBy` is wanted for literature references, post-process the SBML with libsbml, otherwise document that literature is under the `pubmed` key and the rich provenance is in the sidecar.
- Decision records also written as `evidence/decisions/<id>.json`, referenced in commit messages. A report script summarizes the verdict distribution per commit so curation progress is a regression surface.

---

## 8. Guardrails and credibility

- No LLM in the contradiction engine, aggregator, or verdict engine. Tests must pass with the LLM mocked.
- Quote verification is mandatory; a claim without a verified passage does not exist.
- `taxon_distance` and ECO are computed/validated in code, never trusted from the extractor.
- Prior-model inclusion is `MODEL_INCLUSION` evidence, never experimental. Mark inherited reactions (those traceable to a parent model with no organism-specific leaf) explicitly.
- Reviews are `AUTHOR_STATEMENT_*`, counted once via derivative detection, never as independent support.
- Record retrieval dates, API versions, model ids and prompt hashes; results must be reproducible from the cache without network access.
- Short verbatim passages are stored for verification only; do not build any output that reproduces article text.
- Respect API rate limits and the user's cost budget; print a cost estimate before any run over the full corpus.

---

## 9. Milestones

Stop after each milestone and report: what was built, what was found, open questions.

**M0 — Survey.** Repo map; `docs/reviewer_interface.md`; inventory of model files, ID mapping data, existing confidence scores/notes, local datasets; Section 10 questions answered or asked. No pipeline code.

**M1 — Rule core (V1).** Pydantic schemas (Sections 1.2–1.4, 1.7, 6), `id_map` builder, comparability and incompatibility functions, pair classifier, L4 aggregator with the committed weights YAML, L5 verdict engine against the current model, `trace`, `load_bearing`. All Section 5 tests passing. Run end to end on hand-written claims for ~30 reactions chosen with the user (default: a slice spanning central carbon and the phosphatidate hub). Acceptance: deterministic, LLM-free, every verdict traceable.

**M2 — Collection and extraction (V2).** Corpus pull with caching, screening agent, extraction agent, post-validation, gold-set evaluation, L0–L2 population for the pilot reactions, L4/L5 recomputed from real claims, `evidence/PROTOCOL.md`. Acceptance: gold-set metrics reported; every rejected claim has a reason; cost report.

**M3 — Reviewer, validation, decisions (V3).** Roles R1–R3 wired with the typed contracts; `reviewer_disagreements`; L4b validation evidence from essentiality/growth datasets; human queue; patch proposals; SBML annotation writer; apply command gated on approval and on pre/post validation. Acceptance: one full cycle on the pilot slice with at least one decision applied and traced back to leaves.

---

## 10. Questions to ask before M1 (only those not answerable from the repo)

1. Where is the AI reviewer, how is it invoked, and which of its outputs do you trust enough to gate decisions on?
2. Canonical strain set for `taxon_distance = 0`, and which YALI0↔YALI1 mapping file to use.
3. API access: Anthropic key for `pydantic-ai`, NCBI key, Semantic Scholar/OpenAlex, Unpaywall email, institutional full-text route.
4. Which local datasets exist (essentiality screen, growth/phenotype data, 13C flux tables, proteomics) and their formats.
5. Pilot reaction slice for M1/M2.
6. Existing confidence scores or notes in iYali26 to import as `MODEL_INCLUSION` evidence, and the parent-model lineage to mark inherited reactions.
7. Cost budget per full-corpus run and preferred model for screening vs extraction.

---

## Deviations

This section is not part of the spec as received. The text above is
`/Users/david/Downloads/EVIDENCE_PIPELINE_PROMPT.md` copied verbatim on 2026-10-29 (SHA-256
`701663921296662cdf08de102ec36c038d08a84327222743686fcf7a7e936120`). Everything below records where the
implementation differs from the spec, and why. Milestone reports cite the commit of this file. Each entry
gives its source: **user** for a user decision (2026-10-29 answers to the M0 questions, in
[../../reports/evidence_pipeline_m0_20261029/REPORT.md](../../reports/evidence_pipeline_m0_20261029/REPORT.md)
§8), or **proposed** for one awaiting confirmation.

| # | Spec | Implementation | Why | Source |
| --- | --- | --- | --- | --- |
| D1 | `evidence/` package, `evidence/evidence.sqlite`, `evidence/decisions/<id>.json`, `evidence/PROTOCOL.md` | Code in `platform/evidence/`. Weights YAML, subsystem precedence, `PROTOCOL.md`, decision records and a small per-commit assertion/verdict snapshot (JSON or parquet) in `model/curation/evidence_pipeline/`. SQLite in `$IYALI26_RESEARCH_ROOT/artifacts/evidence_pipeline/`, rebuildable deterministically from the cache, committed code and the weights YAML | Repo layout: code under `platform/`, decisions under `model/curation/`, generated data in the research workspace | user |
| D2 | `docs/reviewer_interface.md` | `model/curation/evidence_pipeline/reviewer_interface.md` | `/docs/` is git-ignored (`.gitignore` line 12) | user (layout as proposed) |
| D3 | Verbatim quote ≤ 40 words | ≤ 25 words, plus character offsets into the cached source text so a ±300-character window can be rebuilt for R2 | Less copied text; context comes from the cache when needed | user |
| D4 | `Condition` fields as listed in §1.4 | Adds `medium` and `marker_complementation`; `ESSENTIAL_UNDER` comparability keys on both | Essentiality depends on medium and auxotrophic markers (e.g. PO1f Leu/Ura) | user |
| D5 | `taxon_distance = 0` = W29 and Po1d/f/g/h, "confirm with the user" | {W29, Po1d, Po1f, Po1g, Po1h, PO1f::Cas9}; strain recorded exactly as written | Confirmed; in-house Cas9 strain added | user |
| D6 | `id_map` YALI0 ↔ YALI1 ↔ UniProt ↔ NCBI ↔ symbols | Built on the S2-based `LocusCrosswalk`, with a `cardinality` field (`1:1`, `1:n`, `n:1`). Anything not 1:1 is `UNRESOLVED` for pairing. The 15 non-YALI1 model gene IDs get a disposition each in the M1 report | Avoids ID-era artifacts; the crosswalk already refuses ambiguous aliases | user |
| D7 | `Source.kind` ∈ {article, preprint, thesis, patent, dataset, database_entry, prior_model} | Adds `legacy_case`, with `legacy_scoring=True`, for the 240 essentiality FN cases. Their scores are never merged into live assertions; a case is re-extracted under the new schema when its reactions enter a pilot | Keep prior reviewer work traceable without mixing scoring systems | user |
| D8 | Prior GEMs as `MODEL_INCLUSION` | iYali26's own evidence-status notes are imported as `kind=prior_model`, with `source_id` = model commit and the note text as the passage. Lineage is a `parent_model` chain, and `inherited_from` is set on every reaction whose only leaf is a parent model, so `homology_only` works | Uses the curation notes already in the model | user (iYli21's own ancestry still to be supplied) |
| D9 | `Source` fields | Adds `unpublished: bool`. The W29 13C dataset is `peer_reviewed=False, unpublished=True`. Every report separates internal-data evidence from literature evidence | The dataset's README says the work is unpublished | user |
| D10 | Reviewer integration: R1–R3 on the existing reviewer | The essentiality FN loop's roles are the template, not its code. The skeptic becomes R3, generalized to any proposed action; R1 and R2 are new agents with the same adversarial pattern and the spec's typed contracts. Agent definitions live on `main` as versioned prompt files; the legacy prompts are kept with their hash so the 240 cases stay reproducible. New prompts emit polarity (negative claims required) and contain no scoring rules or cutoffs | Scoring is code; the legacy prompts encode an old contract | user |
| D11 | Reviewer outputs may gate decisions | Nothing gates yet. Before any reviewer output gates, calibrate on the 85 reviewed cases (all 8 passed or downgraded plus 12 random blocked, human-labelled) and report the skeptic's precision | Gating needs measured reliability | user |
| D12 | Reviewer runtime unspecified | Extractor and R1/R2 run on the provider `gem_review.lit` supports; R3 sits behind the typed contract so it can run on a different model family | Model diversity between extractor and gate is intended | user |
| D13 | Sources: PubMed, Europe PMC, Semantic Scholar/OpenAlex, Unpaywall, institutional route, Crossref/Retraction Watch | NCBI, KEGG, BiGG, MetaNetX, Europe PMC OA and Unpaywall only. M2 reports full-text coverage on the pilot corpus before any decision on paywalled papers. Retraction status from Europe PMC / NCBI publication-type metadata instead of Crossref | User-approved access list | user (retraction route: proposed) |
| D14 | Tier E datasets ingested | M2 lists candidate public datasets found in the corpus (Tier A omics papers with supplementary tables) as Tier E candidates; none is ingested until the user picks | User control over data intake | user |
| D15 | Extraction layer | A new extractor reusing `gem_review.lit`'s `check_quote`, HTTP cache and provider setup; current `gem_review.lit` commands and outputs unchanged. `claim_level` is chosen by the LLM with a stored justification and the score is computed in code; organism and strain are extracted by the LLM and `taxon_distance` is computed in code; prompt hash and temperature 0 are recorded; all four species names go into every search. Legacy 0–4 LLM scores are never mixed with derived scores | Reuse what exists; keep earlier reports reproducible | user |
| D16 | Verdicts against "the current model" | Reference model `d8d1bac7` (newest on `main`, 2026-10-29). The verdict engine takes the model as an input and records model path, commit and SHA-256 on every verdict, so verdicts can be recomputed against later commits | No reference model was designated in the repo | user |
| D17 | Gold-set gate "in CI" | GitHub pytest job with the LLM mocked gates merges; memote stays its own job | Repo had memote-only CI | user |
| D18 | L4b essentiality threshold unspecified | One `essentiality_threshold` in the weights YAML, 0.15 (KO/WT strictly below = essential), shared by the benchmark and L4b; L4b reports add a 5/10/15/20 % sweep. The experimental cutoff behind the essentiality calls is recorded as a dataset property and stated separately in every report. Wiring the benchmark tooling to read the YAML comes with L4b (M3) | One prediction, one threshold; 15 % keeps the benchmark history comparable (user correction, 2026-10-29; the earlier 10 % echoed the legacy prompt) | user |
| D19 | Cost: "print a cost estimate" | Budget $20 per full-corpus run, model `claude-fable-5-1`; estimate printed before any run | User budget | user |
| D20 | Comparability rule 1: same normalised subject and predicate | GENE_ENCODES claims are paired within the reaction (the bucket), so claims about different genes of one reaction meet and `exclusive` conflicts are found; assertions stay per (gene, reaction) | §1.6 compares different genes, which a per-(gene, reaction) bucket could never do | M1 implementation |
| D21 | Two bands, functional and parts | Levels in neither band (ORTHOLOGY, AUTHOR_STATEMENT_*, MODEL_INCLUSION, MODEL_VALIDATION) form a third band, "statement". Opposite-polarity pairs across bands are SCOPED_APART(level_mismatch), except parts-negative vs functional-positive, which is ANOMALY_CROSS_LEVEL | §1.5 does not say how the remaining levels pair | M1 implementation |
| D22 | Scope checks before incompatibility | Strain, condition and band only matter when the values conflict; compatible claims agree whatever their scope | Two positive flux claims under different carbon sources both support activity | M1 implementation |
| D23 | Derivative detection | Compatible claims from the same source are AGREE_DERIVATIVE (shared root source) | One paper, one independence class | M1 implementation |
| D24 | Resolution R2→R6 in order, R5 overriding R2–R4 | R5 (assertive rebuttal) is checked first; R2 directness, R3 taxon, R4 independence and R6 recency then narrow the contenders step by step | §2.1 says R5 wins "regardless of R2–R4" | M1 implementation |
| D25 | Assertions with sides but no contradiction | Every side stands; each side carries its own score, and verdicts read all standing sides. After a resolved contradiction, the winner and the sides outside the contradiction stand | §1.5/§1.6: scoped-apart claims both hold (an earlier draft picked one side by year; fixed after review) | M1 implementation |
| D26 | Distance ≥ 2 claims are down-weighted transfers | They enter L4 as ORTHOLOGY with `transfer_from_level` = their own level, so the ORTHOLOGY score rule applies | Makes "down-weighted" concrete | M1 implementation |
| D27 | Verdict rules for scores ≥ 3 and ≤ 1 | Score 2 gives PARTIAL(flag_for_evidence) | §2.2 gives no rule for 2 | M1 implementation |
| D28 | FLUX_13C scores 4 for existence only | Also for FLUX_ACTIVE_UNDER and physiological DIRECTIONALITY (still not for GENE_ENCODES) | Otherwise the most direct flux evidence would score 1 for flux | M1 implementation |
| D29 | Directness lists | STRUCTURE sits with the activity levels for COFACTOR; MODEL_VALIDATION sits just above MODEL_INCLUSION for REACTION_EXISTS and GENE_ENCODES | The spec's level table admits them, its rankings omit them | M1 implementation |
| D30 | Negative resolved gene assertion → PARTIAL(remove_gene) | Not when every strong negative side rests on PHENOTYPE_KO; such a gene stays UNSUPPORTED and is named as "negative knockout only (isozymes possible)" | §1.6: a "no effect" knockout is explained by isozymes | M1 implementation |
| D31 | ACTIVITY_PURIFIED admissible for LOCALIZED_IN only with fractionation | Implemented as the claim's `method` mentioning "fraction" | No separate field in the claim schema | M1 implementation |
| D32 | propose_constrain; R7 network flags | Not produced in M1: existence claims carry no condition yet, and R7 is recorded as "not computed" in every trace | Deferred | M1 implementation |
| D33 | Priority = model impact × uncertainty | Impact = 1 + essential in simulation (stored knockout table, used only when bound to the same model SHA; otherwise the term is left out and the verdict carries `priority_note: impact unknown`) + phosphatidate-hub membership; flux centrality not computed | Reuses the 2026-10-29 knockout run instead of new solves | M1 implementation |
| D34 | Passage char offsets | For cached literature chunks, offsets are relative to the chunk named in `Passage.section` | The cache stores chunks, not whole articles | M1 implementation |
| D35 | Hedged claims "produce weak contradictions only" | Hedged claims score like assertive ones in L4; open question for the user (R1842's GPR is SUPPORTED on one hedged knockout inference) | The spec limits hedged claims only in contradiction strength | open |
| D36 | ECO IDs "verify" | Release `releases/2026-07-10` pinned (SHA-256 `c8b52a2d…`); `eco_terms.tsv` generated by `evidence.eco_table`; MODEL_VALIDATION = ECO:0007669; the test checks existence, non-obsolete and branch | User, 2026-10-29 | user |
| D37 | Dataset converters | W29 13C: positive when the flux interval excludes 0, negative only when both bounds are 0, otherwise no claim (logged). Model notes: keys containing "evidence" (not "counter…") as positive MODEL_INCLUSION. Cas9/Cas12a rows as PHENOTYPE_KO, strain unknown (distance 1); the consensus list is a derivative of the screens | Deterministic, documented rules | M1 implementation |
| D38 | Sidecar tables | Each table holds `(id, payload JSON)`; a `load_bearing` table is added | Simple, deterministic rebuild | M1 implementation |
| D39 | CONTRADICTED when the resolved assertion is negative / disjoint / incompatible | Only when a standing side scores ≥ 3, i.e. organism-specific functional evidence. Weaker negatives (homology, genomic presence, transfers from other species) give UNSUPPORTED or PARTIAL and are named in the reason, never acted on | §2.2 "score ≤ 1 → not removal"; AGENTS.md: deletion needs evidence of non-existence | M1 implementation (after review) |
| D40 | Bounds vs directionality | Bounds are a capability. In vivo flux one way does not contradict a reversible reaction; in vitro reversibility does not contradict a physiologically irreversible one; flux not observed in vivo (a physiological negative) removes no capability. Contradictions: biochemical irreversibility against a reversible model, observed reversibility in vivo against an irreversible model, opposite irreversible directions, a biochemical "cannot run X" against a model that allows X, and any strong evidence that a closed (0, 0) reaction runs | Keeps the biochemical/physiological split of §1.6 in the verdict (refined after the re-review) | M1 implementation |
| D41 | `load_bearing` stored on the L4 node | Stored per verdict (`load_bearing` table and snapshot) as leaves whose removal changes the verdict or an assertion score under it | One verdict can rest on several assertions (GPR) | M1 implementation |
| D42 | Parent-model inclusion | A parent reaction counts only when its species (with compartments) equal the current reaction's, protons excluded (iYali26 rebalanced protonation); otherwise a `parent_reaction_differs` row | AGENTS.md §5: match by substrate, product and compartment, never by ID | M1 implementation (after review) |
| D43 | Quotes ≤ 40 words (D3 sets 25) | Over-long quotes are rejected (`quote_too_long`), never truncated; converters that excerpt long model notes take the first 25 words on purpose | A cut can drop a negation | M1 implementation (after review) |
| D44 | Claim identity | `claim_id` digests every field of the claim; duplicate IDs fail the build; set-valued fields serialise sorted | Determinism and no silent drops | M1 implementation (after review) |
| D45 | Exclusive claim vs model isozymes → CONFLICTING | Only for an exclusive side scoring ≥ 3 and only when the model rule has OR alternatives (an AND complex is not isozymes) | Avoid the top review priority on weak or misapplied exclusivity | M1 implementation (re-review) |
