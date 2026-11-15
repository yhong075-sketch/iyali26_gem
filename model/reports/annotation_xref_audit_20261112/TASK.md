# Model-wide annotation cross-reference audit of e584151f (read-only), 2026-11-12

## 0. Authorization

Task brief from the session that started this one (`local_d48d4f64…`), 2026-11-12. It follows up two cross-cutting
findings of the read-only GPR/name-mismatch review of six consensus-essential genes
(`model/reports/gpr_name_mismatch_leads_20261112/`, uncommitted on `claude/nervous-hellman-589a90`, REPORT.md §8
items 1–2). Neither finding had been checked model-wide.

Scope: investigation and documentation only (root AGENTS.md §2). No model XML, curated-data, builder, GPR or test
change. No LP. No commit or push. Network lookups were limited to public databases (UniProt REST, NCBI E-utilities,
KEGG REST), queried with accessions, locus tags and GeneIDs only. This session did not see the user's own words for
this task.

## 1. Questions (from the brief)

1. Which reactions carry KEGG / Rhea / BiGG / MetaNetX cross-references whose equation disagrees with the reaction's
   own stoichiometry? Suspected cause: Strategy C in `platform/scripts/gem_annotate/annotate_reactions_extended.py`
   (lines 556–576 at `efd7f780`), which accepts a single EC-derived MNXR without a metabolite check. Known case:
   R1711 "raffinose invertase" carries HisA's MNXR103169 / KEGG R04640 / Rhea 15469. Mark the reactions that Strategy C
   annotated, by rebuilding the decision in memory.
2. Which model genes carry a UniProt xref that UniProtKB has since deleted? Four are known: A0A1H6Q312
   (YALI1F28092g), A0A1H6PZC1 (YALI1C27933g), Q6L8F1 (YALI1D05474g), A0A1D8NE95 (YALI1D15063g). For each deleted one,
   resolve the live entry through NCBI Gene (locus tag `YALI1_xxxxx`) → KEGG `yli:GeneID` → UniProt. Never map
   YALI1 ↔ YALI0 by locus number.

## 2. Starting state (checked 2026-11-12)

- Worktree `.claude/worktrees/peaceful-pascal-66d07e`, branch `claude/peaceful-pascal-66d07e`, HEAD
  `cf7b8dd94795d59fda81e00c4189736c3dbb3104` (merge of PR #52), clean. `annotate_reactions_extended.py` is unchanged
  since `efd7f780`; Strategy C is at lines 556–659 (single-candidate branch at 574–576).
- Model audited: `model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_bigg_annot_v2_mcnaugh_rxns_v3_20261103.xml`,
  SHA-256 `e584151fc79475fc2daa59e76c8c84ebfd06b0f3cb1fbd75fbb69cc9f31039d2`, read by cobra 0.30 through a
  `pathlib.Path` (cobra rejects str paths of 260+ characters on macOS).
- Source model: `model/source/iyali26.xml`, SHA-256 `5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d`
  (the build input recorded in e584151f's `.build.json`).
- Build record: e584151f's `.build.json` (built 2026-11-04, `--offline --no-solve`, HEAD `048e18ec` with `cli.py` and
  `main.py` dirty). 69 of its 71 recorded builder sources and all 144 recorded data files are byte-identical in this
  worktree; `cli.py` and `main.py` differ. The `main.py` diff `048e18ec..HEAD` only adds arguments and late steps
  (`vatpase_a1`, `r2121`, `r1172`), none before Priority 4e.
- MetaNetX 4.5 tables used by the builder, `$IYALI26_RESEARCH_ROOT/reference/metanetx/`: `reac_prop.tsv`
  (`8582cc18…`), `reac_xref.tsv` (`2610e299…`), `chem_xref.tsv` (`ba2451ff…`), `chem_prop.tsv` (`8463f8ee…`),
  `chem_depr.tsv` (`de0778d3…`).
- Builder cache `$IYALI26_RESEARCH_ROOT/cache/data/`: every file predates the e584151f build (latest 2026-07-21).
- BiGG crosswalk decisions: `model/curation/annotation/bigg_crosswalk_annotation_decisions_20261103.csv` (read only).
- W29 tag audit: `model/curation/annotation/w29_tag_protein_audit.csv` (read only).
- Reference model: none designated (`model/candidates/README.md`). The question is about annotation, not scoring, so
  no essentiality cutoff or experimental list is used.

## 3. Method

**Q1, static fingerprint comparison** (`scripts/xref_fingerprint_audit.py`). For every reaction with a
`metanetx.reaction`, `kegg.reaction`, `rhea` or `bigg.reaction` value:
- Model side: group the reaction's metabolites by chemical identity (compartment copies carry separate IDs). Each
  compound's identity is its `metanetx.chemical` (resolved through `chem_depr.tsv` and the annotator's
  `_MNXM_NORMALIZE`), plus MNXMs reached from its other database IDs through `chem_xref.tsv`, its InChIKey skeleton,
  its name, MetaNetX synonyms and a stereo-free name.
- MNXR side: the equation in `reac_prop.tsv`. KEGG / Rhea / BiGG IDs are mapped to MNXRs through `reac_xref.tsv` and
  compared the same way.
- H+ and H2O are ignored. A fixed cofactor list (ATP/ADP/AMP, NAD(P)(H), Pi, PPi, CO2, O2, NH4+, CoA, HCO3−,
  nucleotides, FAD/FADH2, Na+, K+; IDs checked against `chem_prop.tsv` names) does not count as shared substrate.
- Classes: `match`, `partial`, `disjoint` (no shared non-cofactor compound), single-compound variants, and
  `mnxr_covers_only_h_h2o`. Unmatched substrates are compared by formula without H to separate "same reaction, other
  IDs" from "different compounds".

**Q1, provenance** (`scripts/replay_annotation.py`). Replay the builder stages from the source model through
Priority 4e in memory, in `main.py` order, on a scratch copy of the cache. Record which `strategy = "…"` line of the
real annotator code fires for each reaction, by line tracing (the annotator logic is not copied). Compare the
replayed MNXR with e584151f's. No model is written; the cache copy and the real cache were hashed before and after
and are unchanged.

**Q2** (`scripts/uniprot_status.py`, `scripts/resolve_replacements.py`). Every gene `uniprot` value is checked
against UniProtKB (batch `accessions` endpoint, then one request per accession missing from it, to read
`entryType` / `inactiveReason`). For each inactive one: NCBI Gene esearch on the `YALI1_` locus tag, keeping the
GeneID whose esummary lists that exact tag → KEGG `yli:GeneID` (entry and `conv/uniprot`) → UniProt status of the
linked accession. Cross-check: UniParc ID of the deleted entry vs the replacement (equal UniParc means identical
sequence). Join the W29 tag-audit flag.

Evidence tiers: *verified* (computed or fetched this session; output saved), *recalled* (from an earlier
session's file, not re-checked), *inferred* (reasoned).

## 4. Done means

`REPORT.md`, `reaction_xref_mismatches.csv`, `deleted_uniprot_xrefs.csv`, a code-improvement proposal for Strategy C
with a test design (in REPORT.md), and a dated `model/STATE.md` entry, all uncommitted. Scripts, raw responses and
full tables: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/annotation_xref_audit_20261112/`.
