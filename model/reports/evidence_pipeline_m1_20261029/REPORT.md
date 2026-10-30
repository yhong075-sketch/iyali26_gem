# Evidence pipeline, M1: rule core — report (2026-10-29)

Scope, authorization and plan: [TASK.md](TASK.md). Spec: [SPEC.md](../../curation/evidence_pipeline/SPEC.md).
The verbatim body has SHA-256 `70166392…`; the Deviations table D1–D44 records every place the code differs from
the spec. The file is not committed yet, so the spec commit is still to be cited.

**Verification level.** Everything below was run or read this session unless marked otherwise. Nothing is
committed. No LLM was called, no model XML was written, and the only network use was the approved ECO download.

## 1. What was built

`platform/evidence/` is a deterministic package with no LLM. It has 13 modules:

| Module | Role |
| --- | --- |
| `schemas.py` | Pydantic v2 records for L0–L6. It covers the 8 predicates, 19 claim levels, `Condition` (+ `medium`, `marker_complementation`), a value union keyed on the predicate, `Claim`, `Passage` (≤ 25 words with char offsets), and `Source` (incl. `legacy_case`, `prior_model`, `unpublished`). It also holds the R1–R3 reviewer contracts as types only. |
| `weights.py` + [`evidence_weights.yaml`](../../curation/evidence_pipeline/evidence_weights.yaml) | The one committed rule table (version `evidence-weights-v1 2026-10-29`). It maps each level to its ECO term, branch and band, and holds admissibility, directness tiers, confidence rules, taxonomy, the compartment map, subsystem precedence, **`essentiality_threshold: 0.15`** (sensitivity 5/10/15/20 %) and the priority weights. |
| `taxonomy.py` | Computes `taxon_distance` in code. Distance 0 = W29, Po1d/f/g/h and PO1f::Cas9; a null strain is never W29. |
| `id_map.py` | Built on `LocusCrosswalk`, and checked to agree with it on every 1:1 link. Each link has a cardinality (1:1, 1:n, n:1, n:m); only 1:1 links resolve. UniProt and NCBI IDs resolve through the model annotations. |
| `rules.py` | Comparability, value conflict per predicate, and the pair classifier. The classifier handles derivative claims, 5-year independence classes, supersession and strong vs weak. |
| `aggregate.py` | The L4 forward pass. Each side gets its own score. Rule order is R5 → R2 → R3 → R4 → R6; R7 is advisory only. Sides that are only scoped apart all stand. |
| `verdicts.py` | The L5 engine. Elements are existence, GPR, compartment, bounds, cofactor and subsystem. The model is an input, and every verdict carries the model path, commit and SHA-256. CONTRADICTED needs a standing side scoring ≥ 3. |
| `pipeline.py` | The forward pass plus `trace`, `load_bearing` (leave-one-out) and `chain_quality`. |
| `converters.py` | L0–L2 from real data, every claim checked against a verified passage (§3). |
| `store.py`, `__main__.py` | The SQLite sidecar with the spec's tables plus `load_bearing`, rebuilt from scratch; content hash. CLI: `pilot`, `build`, `trace`. |
| `eco_table.py`, `pilot.py` | Generators for `eco_terms.tsv` and `pilot_reactions.csv`. |

Committed data in `model/curation/evidence_pipeline/`:
- the weights YAML;
- `eco_terms.tsv`: ECO `releases/2026-07-10`, 2,254 terms, generated from the OBO cached in `reference/eco/` (SHA-256 `c8b52a2d…`);
- `pilot_selection.yaml` and `pilot_reactions.csv`;
- `handwritten_claims.yaml`;
- `w29_13c_reaction_map.csv`;
- `model_lineage.yaml`;
- the snapshot `snapshots/ff7fc62a77c9_d8d1bac7.json`.

CI: `.github/workflows/pytest-evidence.yml` runs only the evidence tests (no LLM, no research workspace, no
model file); memote is unchanged.

Commands (from `platform/`, with `IYALI26_RESEARCH_ROOT` set):

```
python -m evidence pilot
python -m evidence build --report-dir ../model/reports/evidence_pipeline_m1_20261029
python -m evidence trace --db $IYALI26_RESEARCH_ROOT/artifacts/evidence_pipeline/evidence.sqlite R1842
```

## 2. Pilot (30 reactions; reference model `d8d1bac7`, commit `ff7fc62`)

[pilot_reactions.csv](../../curation/evidence_pipeline/pilot_reactions.csv) has three groups:
- **Essentiality (12):** every skeptic-passed case (4) and every downgraded case (4), plus 4 blocked
  single-reaction cases drawn with seed 20261029. Reactions: R306, R217, R1976, R357, R490, R2143, R699, R1968, R376,
  R792, R12, R374.
- **Phosphatidate hub (12):** R353, R1846, R1843, R2145, R1773, R1720, R196, R1841, R1842, R1743, R2168, R2144.
- **Central-carbon controls (6), each with a W29 13C row:** R326, R313, R694, R325, R236, R487.

## 3. Evidence in the build

| Converter | Claims | What they are |
| --- | --- | --- |
| Parent model iYli21 (`model/source/iyli21.xml`) | 102 | MODEL_INCLUSION for reaction presence, GPR genes and groups. A reaction counts only when its species match the current model's, protons excluded. |
| Current model notes (`d8d1bac7`) | 5 | MODEL_INCLUSION from evidence-status notes (kind `prior_model`, role self). |
| W29 13C-MFA (internal, `unpublished=True`) | 16 + 16 derived | FLUX_ACTIVE_UNDER per carbon source, plus derived REACTION_EXISTS. Five intervals include 0 and are logged as non-claims. |
| Cas9 / Cas12a screens and consensus list | 104 | ESSENTIAL_UNDER (gene level). The consensus list is a derivative of the screens. |
| Hand-written literature claims | 9 | From cached `gem_review` chunks for R2168 (DGA2), R1841 (PSS1) and R1842 (SPO14). Written by Claude in this session, **not human-verified**. Quotes are verified by code. |
| Legacy essentiality cases | 0 claims, 240 sources | `legacy_scoring=True`, kept out of live assertions (D7). |

Totals: 250 sources, 213 passages, 252 claims, 5 rejected rows (all `interval_includes_zero`), 176 pair relations,
148 assertions, 150 verdicts. Every claim's subject resolved.

Internal data and literature are kept apart:
- The 13C claims are internal and unpublished. They carry the existence verdicts for R326, R313, R694, R325, R236 and R490.
- Literature claims carry R1841 and R2168 (existence and GPR) and R1842 (GPR).
- Everything else rests on model inclusion.

Experimental cutoffs behind the essentiality calls are dataset properties, recorded on the source:
- **Screens:** observed boundaries in the normalised file are Cas9 essential at FS ≤ −2.919 and Cas12a at FS ≤ −1.855. The paper's own calling rule was not re-read.
- **Consensus list:** essential in ≥ 2 of 3 screens, positive-only.

## 4. Results

| Element | SUPPORTED | PARTIAL | UNSUPPORTED | CONTRADICTED | CONFLICTING |
| --- | --- | --- | --- | --- | --- |
| existence | 8 | 0 | 22 | 0 | 0 |
| gpr | 3 | 0 | 27 | 0 | 0 |
| compartment | 0 | 0 | 30 | 0 | 0 |
| bounds | 0 | 0 | 30 | 0 | 0 |
| subsystem | 27 | 0 | 3 | 0 | 0 |

Per-verdict table: [pilot_verdicts.tsv](pilot_verdicts.tsv).

What the results mean:
- **Existence SUPPORTED (8):**
  - five central-carbon controls and R490, from 13C (score 4, internal data);
  - R1841 and R2168 from the hand-written knockout claims (score 3).
- **R487 (NAD-isocitrate dehydrogenase) stays UNSUPPORTED.** Its 13C interval includes 0 under all three carbon sources, so the rule makes no flux claim.
- **The other 22 existence verdicts rest only on iYli21 inclusion (score 1).** They are UNSUPPORTED, which is never a removal proposal. 91 of 148 assertions are `homology_only` and `inherited_from=iYli21`.
- **GPR SUPPORTED (3):** R1841 and R2168 from deletion phenotypes. R1842 rests on **one hedged, indirect knockout inference**: higher PC after deleting SPO14 in an engineered strain. `load_bearing` names exactly that claim (open question 1).
- **Subsystem:** the 27 SUPPORTED verdicts only show agreement with iYli21's groups; the verdict reason says so. Compartment and bounds have no claims yet, so all 60 verdicts are UNSUPPORTED.
- **Real contradictions:** 12 genes where the Cas9 and Cas12a assays of the same screen paper disagree. Each gives CONTRADICTION_STRONG with R2–R6 tied (same level, same source), so the assertion is CONFLICTING with its score capped at 2. They are gene-level ESSENTIAL_UNDER assertions and change no reaction verdict in M1. They are the first R1 cases for M3.
- **Comparison with the prior reviewer** (essentiality group): the legacy skeptic marked all 12 cases `needs_more_evidence`. The new pipeline gives existence UNSUPPORTED for 11 of them and SUPPORTED for R490 (from 13C), and GPR UNSUPPORTED for all 12.
  - The two answer different questions: why an essential gene is predicted non-essential, versus whether the model element is evidenced.
  - A real comparison needs the M2 re-extraction of those cases' literature.
- **Priority** combines model impact (genes essential in simulation, from the 2026-10-29 knockout table bound to `d8d1bac7` at the 15 % threshold) with membership of the phosphatidate hub. R196 and R353 rank highest.

## 5. Tests, determinism and review

- **Evidence tests:** 67 passed (`tests/test_evidence_rules.py`, `test_evidence_pipeline.py`, `test_evidence_eco.py`), both with the platform venv interpreter and in a clean venv holding only pydantic, PyYAML and pytest (the CI environment). They include:
  - the 12 cases of spec §5 (case 3 shows that skipping normalisation turns the agreement into a contradiction);
  - an ECO check that every term exists, is not obsolete and sits in its branch;
  - determinism across hash seeds, trace down to L0, and load-bearing leaves;
  - one test per review finding.
- **Determinism of the real build:** identical content hash (`f8af81d0…`) and a byte-identical snapshot under PYTHONHASHSEED 0 and 7.
- **Full suite:** see §5a.
- **Independent review:** a fresh read-only reviewer (briefed from `.claude/agents/code-reviewer.md`) found 17 issues. All were fixed, or recorded as deviations where they are design choices. The main ones:
  - a weak or transferred negative could produce `propose_remove` (now CONTRADICTED needs score ≥ 3, D39);
  - negative sides were reported as SUPPORTED for compartment and bounds;
  - scoped-apart sides were decided by publication year (now all stand, D25);
  - hash-seed-dependent serialisation of set values;
  - claim-ID collisions and truncated quotes (silent drops);
  - an exclusive-claim rule that needed a second claim, and inadmissible claims creating conflicts;
  - parent inclusion matched by ID only (now by species, D42).
  The re-review result is in §5b.

### 5a. Full suite

`python -m pytest` from `platform/` (platform venv interpreter, `IYALI26_RESEARCH_ROOT` set):
- **Existing tests:** 324 passed, 19 failed, 20 errors, 31 skipped. This equals the recorded 2026-10-29 baseline at `f212874` (STATE.md), and a JUnit comparison with this session's earlier run shows no test changing outcome. The 19 failures and 20 errors predate this work.
- **Evidence tests:** 61 passed in that run; the 6 re-review tests were added afterwards and pass with the others (67).
- **Baseline checkout:** a pristine `HEAD` checkout in the scratch folder errored on 31 extra tests, all `OSError` from cobra refusing the over-long scratch path (macOS 260-character limit). That run is therefore not used as the baseline.

### 5b. Re-review

The same reviewer re-checked the fixes. It found 15 of the 16 checked findings fixed (verified by re-running its probes
or reading the code) and reported 7 new or remaining issues. Six are fixed, with a test each:
- a closed (0, 0) reaction is now contradicted by evidence that it runs;
- a negative physiological direction claim no longer removes a capability;
- an unmapped gene's UniProt/NCBI ID no longer resolves;
- an exclusive claim must score ≥ 3 and face OR isozymes, not AND subunits;
- each rejection keeps its own row;
- a resolved contradiction scores from its winner.

One is a disagreement, left as is and reported here. `metadata_selection_evidence_status` notes ("User-selected model
version; not new validation of chemistry") still become positive MODEL_INCLUSION existence claims. That claim says only
"present in the model", which is true, and the user asked for the evidence-status notes to be imported as
MODEL_INCLUSION sources. None of these notes are on the pilot reactions. If the user prefers, the filter can drop
that key.

## 6. Statement levels

- **Provable / verified by tests:** the rule behaviours of spec §5; the LLM is absent from the rule core; deterministic rebuilds; every verdict traces to L0 sources; quote verification (whitespace-normalised exact match); ECO terms in their branches.
- **Observed in runs:** the pilot numbers above, on `d8d1bac7` with the research workspace as of 2026-10-29.
- **Not tested:** cofactor verdicts and subsystem reassignment on real data (no such claims yet); R7 network flags (not computed); `propose_constrain` (not produced).

## 7. Could the method be improved, and how to test it

- **Hedged claims score like assertive ones (D35).** A hedged knockout inference can carry a GPR to SUPPORTED (R1842). Option: cap hedged claims one score lower. Test: count verdicts that change on the pilot, and check that none of the 12 §5 cases moves.
- **The score-≥ 3 gate for CONTRADICTED (D39)** is a reading of §2.2. It will be tested properly when M2 brings real negative literature claims; R2 red-team results will show whether it hides real contradictions.
- **iYli21's per-reaction origin** (yeast-GEM transfer, BlastKOALA, Biolog gap-fill) needs the iYli21 supplement, which is not in the workspace. The planned iYli21 ↔ iYali4 GPR audit through the crosswalk is M2 work.

## 8. Open questions for the user

1. Should hedged claims be capped below SUPPORTED (D35)? R1842's GPR is the live example.
2. Is reading "CONTRADICTED only with score ≥ 3" acceptable (D39), and "bounds are a capability" (D40)?
3. The 15 non-YALI1 model genes, see [id_map_non_yali1.tsv](id_map_non_yali1.tsv):
   - 5 resolve 1:1 to a YALI1 gene (YALI0B10406g → YALI1B14044g, YALI0D17490g → YALI1D21476g, YALI0E16192g → YALI1E19360g, YALI0E18568g → YALI1E22238g, YALI0F04114g → YALI1F06244g).
   - 10 are UNMAPPED: the YALI0 genes YALI0D12573g, YALI0D27654g, YALI0E07787g and YALI0E12442g, plus YALIUNK1/3/4, YALIfMp29 and YalifMp05/06.

   Proposed disposition (an annotation-line task, not done here):
   - rename the five 1:1 genes to their YALI1 IDs in a curated annotation step;
   - for the four unmapped YALI0 genes, look up the W29 locus by sequence (BLAST to the W29 proteome);
   - map YALIfMp/YalifMp to the YALI1 M-locus IDs by the mitochondrial genome annotation;
   - investigate the YALIUNK placeholders.
4. Wiring the benchmark tooling to read `essentiality_threshold` from the YAML is planned with L4b (M3). The value equals the benchmark's 15 % today.
5. Commit and PR: nothing is committed. Which branch should M0/M1 go to?

Large outputs: the SQLite database is at `$IYALI26_RESEARCH_ROOT/artifacts/evidence_pipeline/evidence.sqlite`, rebuildable with the
command above. Only the small snapshot and the TSVs are in the repo.
