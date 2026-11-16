# Model-wide annotation cross-reference audit of e584151f: report, 2026-11-12

Read-only audit. No model, curated-data, builder or test change, no LP. Scope, inputs and SHAs: [TASK.md](TASK.md).
Tables: [reaction_xref_mismatches.csv](reaction_xref_mismatches.csv) (363 rows) and
[deleted_uniprot_xrefs.csv](deleted_uniprot_xrefs.csv) (690 rows). Scripts, raw REST responses and the full
per-reaction table (`out/reaction_xref_fingerprint_all.csv`, 1,632 rows):
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/annotation_xref_audit_20261112/` (`MANIFEST.sha256` `22ae9a6c…`).

**Evidence tiers.** *Verified*: computed or fetched this session, output saved. *Recalled*: from an earlier
session's file, not re-checked. *Inferred*: reasoned, for example my reading of two equations side by side.

## 1. Summary

**Q1.** Of 2,346 reactions, 1,632 carry a MetaNetX, KEGG, Rhea or BiGG cross-reference. Their compound sets were
compared with the referenced MNXR equations (verified, static):

| Fingerprint class | Reactions | Meaning |
| --- | ---: | --- |
| match | 1,192 | same compounds |
| match_plus_unannotated | 86 | MNXR compounds all present; the model has extra compounds without IDs |
| partial | 200 | shares at least one specific compound, differs elsewhere |
| disjoint | 66 | no shared compound beyond H+/H2O and cofactors |
| carrier_only | 53 | shares only cofactors or group carriers (acetyl-CoA, ACP, IPP, SAM, …); the MNXR may be a carrier-only stub |
| disjoint_single_compound / single_compound_formula_match | 12 / 5 | one-compound transports whose compound ID differs (formula equal in 5) |
| mnxr_covers_only_h_h2o | 7 | MNXR is H+ only; the reaction moves a cargo too |
| unassessable / absent | 6 / 5 | H+-only reaction against the PMF pseudo-compound / no MNXR but other xrefs |

I read the equations of all 149 rows in the disjoint, carrier-only, single-compound, H+-only and unassessable
classes and labelled each (inferred; one note per row in the table). Rows that share one MNXR and belong to one
chain-length family were read as a group.

| Adjudication | Reactions | Example |
| --- | ---: | --- |
| other_reaction (different chemistry or mechanism) | 20 | R1711 raffinose invertase → HisA isomerase |
| other_substrate (same reaction type, different specific substrate) | 63 | R1499–R1510 Cn 2-enoyl-CoA hydratase → C13 2,6-dienoyl-CoA hydratase |
| incomplete_mnxr (MNXR covers only part of the reaction) | 14 | R1999 keto-gluconate/H+ symport → H+ exchange |
| same_reaction_other_ids (generic, other-organism or synonym IDs) | 43 | R1170 homocitrate synthase; MetaNetX files homocitrate as 2-hydroxybutane-1,2,4-tricarboxylate |
| model_side_unassessable / not_judged | 3 / 6 | m2057–m2060 have no name; MNXR03 (PMF) set by a reviewed crosswalk decision |

So **97 reactions carry an MNXR (and its KEGG/Rhea/BiGG IDs) that describes another reaction or substrate**
(class verified; labels inferred). Another **16 reactions** have a correct MNXR but keep KEGG/Rhea IDs from the
MNXR the annotator first chose, because the BiGG crosswalk step rewrote only `metanetx.reaction` and
`bigg.reaction` (§2.5). The remaining 198 partial rows are listed but not adjudicated.

**Strategy C's single-candidate branch is not the main source.** In-memory replay of the build attributes each
annotation to a strategy line (verified, §2.3). The single-candidate branch (lines 574–576) annotated 8 reactions.
One is wrong (R1711), 5 match and 2 are partial. One of the partials is R653's original pick, which a crosswalk
decision later replaced. The likely-wrong set comes mostly from paths that never compare
metabolites, or that compare only the annotated ones:

| Annotator path | Annotated | Likely wrong | Stale KEGG/Rhea |
| --- | ---: | ---: | ---: |
| 4a Strategy B, reaction name = MetaNetX short name (`reactions.py`) | 814 | 38 | 13 |
| 4e Strategy C Layer 3, name-token Jaccard among EC candidates | 52 | 20 | 2 |
| 4e Strategy C Layer 2, Jaccard over annotated MNXMs only | 81 | 18 | 1 |
| 4a Strategy C, exact fingerprint of annotated MNXMs only | 171 | 8 | 0 |
| 4e Strategy B5, transport fingerprint of annotated MNXMs only | 278 | 7 | 0 |
| 4e Strategy B6, transport name substring | 66 | 5 | 0 |
| 4e Strategy C single candidate (no check) | 8 | 1 | 0 |
| 4e Strategy D, A, B and C Layer 1 | 163 | 0 | 0 |

**Q2.** 833 of 1,081 genes carry a UniProt xref (835 accessions). **690 of the 835 are deleted in UniProtKB**
(verified, REST, 2026-11-12): 680 "Not part of a reference proteome", 9 "Redundant sequence", 1 "Deleted from
sequence source (ENSEMBL)". Only 145 are active. All 690 are entries of the cached W29/CLIB89 proteome
UP000182444 (dated 2026-04-28), which the gene annotator reads first (`config.py:228`); none is in the cached
CLIB122 proteome (verified, cache scan). Through NCBI Gene → KEGG → UniProt,
every one of the 690 genes resolved to exactly one GeneID:
- 679 have a live UniProt entry, all CLIB122 (taxid 284591). 544 have the identical sequence (same UniParc ID); 135
  have a different sequence.
- 2 have no KEGG→UniProt link, but a live entry with the identical sequence exists (UniParc): YALI1D32594g →
  Q6C7X2 and YALI1E31685g → Q12724.
- 9 have no live UniProt entry, among them YALI1D15063g (cytochrome b5).

The four known cases reproduce the sibling review: A0A1H6Q312 → Q6C0W3, A0A1H6PZC1 → Q6CBC8, Q6L8F1 → F2Z6C9, all
with identical sequences, and A0A1D8NE95 → none.

## 2. Q1: reaction cross-references

### 2.1 Method in brief

Details are in TASK.md §3 and the script header.
- **Model side.** Each reaction's metabolites are grouped by chemical identity, because compartment copies carry
  separate IDs. A model compound matches an MNXR compound by any of: MNXM ID (deprecation-resolved, plus MNXMs reached
  from its KEGG/ChEBI/BiGG/… IDs through `chem_xref.tsv`), InChIKey skeleton, name, MetaNetX synonym, or stereo-free
  name. H+ and H2O are ignored.
- **Cofactors and carriers.** A fixed list does not count as shared substrate. Cofactors: ATP/ADP/AMP, NAD(P)(H),
  Pi, PPi, CO2, O2, H2O2, NH4+, CoA, HCO3−, nucleotides, FAD, Na+, K+. Carriers: acetyl-/malonyl-CoA, IPP, ACP,
  Glu, 2-OG, Gln, SAM, SAH, GSH/GSSG, thioredoxins, UDP-Glc, GDP-Man, CDP, ubiquinone-9/ubiquinol-9. All IDs were
  checked against `chem_prop.tsv`.
- **KEGG, Rhea and BiGG IDs** are mapped to MNXRs through `reac_xref.tsv` and compared the same way.
- **Why the matching is lenient.** The first, ID-only pass called 132 reactions disjoint. Many were the same
  reaction under different IDs: isoamyl alcohol vs 3-methylbutanol, transports whose compartment copies were not
  grouped, (S)-3-hydroxyisobutyrate vs (S)-3-hydroxy-2-methylpropanoate. Each refinement is in the script. The
  remaining disjoint and carrier-only rows were read by eye.

### 2.2 Evidence for the main families

Equations as stored (model) and in `reac_prop.tsv` (MNXR). Every row of the table carries both equations.

- **R1711 raffinose invertase** (GPR YALI1F07827g). Model: raffinose + H2O → melibiose + D-fructose. MNXR103169:
  1-(5-phospho-β-D-ribosyl)-5-[(5-phospho-β-D-ribosylamino)methylideneamino]imidazole-4-carboxamide ⇌ its
  ribulosyl isomer (HisA, EC 5.3.1.16). The KEGG R04640, Rhea 15469–15472 and BiGG PRMICI/PRICIh IDs all map to
  the same MNXR (verified). Path: Strategy C single candidate, from gene EC 5.3.1.16 → candidates
  `[MNXR103169]` (verified, replay). That KEGG R04640 and Rhea 15469 are HisA was verified on the KEGG and Rhea sites
  by the sibling review (recalled).
- **R1499–R1510**, peroxisomal 2-enoyl-CoA hydratases C4–C26, all → MNXR118924: 2-trans-6-trans-tridecadienoyl-CoA +
  H2O = (3R)-hydroxy-6-trans-tridecenoyl-CoA. Same EC (4.2.1.119), a C13 dienoyl substrate that no model step uses.
  Path: Strategy C Layer 3 (name similarity) among the EC's candidates. Related: R1554 (dodecadienoyl → the
  tetradecadienoyl MNXR118869), R1490/R1491/R1493/R1496 (C16–C22 acyl-CoA oxidases → the decanoyl-CoA MNXR151473),
  R1545 (→ 3-methylcrotonyl-CoA thiolase).
- **R218–R235**, the 18 cis-prenyltransferase steps of the dolichol route, all → MNXR174168 "rubber
  allyltransferase", whose equation is only `PPi = IPP`. Path: Strategy C Layer 2. The polyprenyl diphosphates have
  no MNXM, so the Jaccard ran on {IPP, PPi} alone.
- **R1540–R1542, R1544**, thiolases → MNXR135489, `acetyl-CoA = CoA`, the stub of heparan-α-glucosaminide
  N-acetyltransferase. Path: primary exact fingerprint, again because the acyl-CoAs have no MNXM. Same mechanism:
  R288 (→ the GDP-mannose:mannan mannosyltransferase stub), R2097/R2098 (`SAH = SAM`), R1460, R2186.
- **R2217–R2241**, nine acyl-ACP hydrolases (C8–C18:1) → MNXR99106, the stearoyl-ACP (C18:0) hydrolase; none of the
  nine is C18:0. Path: primary Strategy B (name). Related: R2216/R2219/R2238/R2240 acyl-ACP transports → the
  dodecanoyl-ACP MNXR99088.
- **R706–R711, R1913, R1914**, ribonucleotide reductase on ADP/CDP/GDP/UDP (EC 1.17.4.1) → MNXR173948, generic
  ribonucleoside *tri*phosphates (EC 1.17.4.2). Path: primary Strategy B.
- **Other wrong-substrate pairs**, all primary Strategy B:
  - R619/R621/R623 thioesterases (lauroyl/myristoyl/octanoyl-CoA) → pristanoyl-CoA MNXR127570.
  - R31/R33 cIMP/cCMP phosphodiesterase → the cGMP MNXR189614.
  - R144/R2252 amidases → the 2-phenylacetamide MNXR152978.
  - R2213 nitrilase → indole-3-acetonitrile.
  - R665 PI 3-kinase → PI5P 3-kinase.
  - R182–R184 C-4 methyl oxidase steps → a lumped oxidation + decarboxylation of 3-dehydro-4α-methylzymosterol.
  - R1721 cardiolipin synthase (CDP-DAG + PG) → the bacterial 2 PG → CL + glycerol.
- **Seven proton-coupled transports** of cargoes without an MNXM (R1999, R2019, R2020, R2027, R2028, R2037,
  R2038) → MNXR02 "H+ exchange". Path: B5, whose fingerprint dropped the unannotated cargo.

### 2.3 Provenance replay

The builder stages from `iyali26.xml` through Priority 4e ran in memory, in `main.py` order. Line tracing recorded
which `strategy = "…"` line of the real annotator code fired for each reaction. No model was written. The cache was
a hashed scratch copy, and the cache copy and the real cache are unchanged (verified).

Fidelity (verified):
- The replay annotated 1,633 reactions; e584151f has an MNXR on 1,627.
- 1,448 final MNXRs equal the replayed one.
- 162 differ exactly by a BiGG crosswalk `write_his` decision (each decision's `before` equals the replayed MNXR and
  its `after` the final one).
- 1 differs otherwise (R1964, D-strategy pick vs final MNXR198271; not traced).
- 16 belong to reactions added after Priority 4e.
- 22 replayed annotations are absent from the final model: 12 reactions removed, 10 annotations dropped later.

All 97 likely-wrong rows keep the MNXR the replay assigned, so their strategy attribution is exact for this build.

**Strategy C's single-candidate branch fired 8 times** (verified):

| Reaction | Name | Gene EC | MNXR | Class |
| --- | --- | --- | --- | --- |
| R59 | 5'-phosphoribosylformyl glycinamidine synthetase | 6.3.5.3 | MNXR165854 | match |
| R315 | fumarase, cytoplasmic | 4.2.1.2 | MNXR99705 | match |
| R586 | nucleoside diphosphatase | 3.6.1.42 | MNXR188495 | partial (GDP-specific MNXR) |
| R653 | phosphoribosylaminoimidazole-carboxylase | 4.1.1.21 | MNXR108735 (ATP-independent AIR carboxylase) → MNXR136160 (crosswalk) | pick partial; final match |
| R719 | saccharopine dehydrogenase (NADP, L-glutamate) | 1.5.1.10 | MNXR146359 | match |
| R1711 | raffinose invertase | 5.3.1.16 | MNXR103169 | disjoint, other_reaction |
| R2137 | UDP-D-glucuronate carboxy-lyase | 4.1.1.35 | MNXR189894 | match |
| R2202 | carbonate hydro-lyase | 4.2.1.1 | MNXR100482 | match |

### 2.4 What the 97 have in common (inferred)

The 97 likely-wrong rows involve 30 distinct MNXRs: 18 reactions share MNXR174168, 12 MNXR118924, 9 MNXR99106 and
8 MNXR173948. 83 of them have a GPR. 11 carry a KEGG reaction ID, 20 Rhea IDs and 53 BiGG IDs from the wrong MNXR.
Two mechanisms explain them:
1. **Name-only acceptance.** Primary Strategy B, B6 and Strategy C Layer 3 pick an MNXR by its name and never check
   metabolites. Chain-length families then collapse onto whichever chain-specific MNXR carries the generic name.
2. **Fingerprints that silently drop unannotated metabolites.** Primary C, Strategy C Layers 1 and 2, B5 and D build
   the query set only from metabolites with an MNXM. A reaction whose specific compound has no MNXM is matched on
   its cofactors and carriers alone, and lands on a stub MNXR (`acetyl-CoA = CoA`, `PPi = IPP`, H+ only).

### 2.5 Stale KEGG/Rhea after the crosswalk rewrite

`apply_bigg_crosswalk_annotations` writes only `bigg.reaction` and `metanetx.reaction`
(`bigg_crosswalk_annotations.py`, `KEYS`). It rewrote the MNXR of 162 reactions. Their KEGG, Rhea, SEED, MetaCyc,
SABIO-RK and VMH keys still come from the annotator's original MNXR. In 16 of them the leftover KEGG or Rhea ID is
disjoint from, or carrier-only against, the reaction (verified):
- R29, R30: cAMP phosphodiesterase still carries cGMP's KEGG R01234 and Rhea 16957.
- R618, R620, R622: still carry Rhea 40415 (pristanoyl-CoA).
- R99, R1487: still carry Rhea 40179 (decanoyl-CoA oxidase).
- Also R95, R577, R636, R649, R1131, R1908, R1909, R2062 and R2210.

This affects another line: the crosswalk step is this repository's code, but the decisions are Andrew McNaughton's
reviewed handoff.

## 3. Q2: deleted gene UniProt xrefs

- **Status** (verified). The batch `uniprotkb/accessions` endpoint returns only active entries. Each of the 690
  missing accessions was fetched on its own, and every one returned `entryType: Inactive`, `DELETED`. Example:
  [A0A1H6Q312](https://rest.uniprot.org/uniprotkb/A0A1H6Q312.json), "Not part of a reference proteome", UniParc
  UPI00003BCE85. By prefix: A0A1D8 (495), A0A1H6 (184) and 11 others, among them Q6L8F1.
- **Source** (verified, cache scan).
  - All 690 deleted accessions are in the cached `uniprot_UP000182444.json` (W29/CLIB89, the "Other proteome" in
    `config.py:228`, cached 2026-04-28); 687 are also in `idmapping_results.json`. None is in the cached CLIB122
    proteome UP000001300.
  - All 681 live replacements (679 from KEGG plus 2 from UniParc) are already in that cached CLIB122 proteome.
  - An offline rebuild therefore keeps writing the deleted accessions (inferred from the W29-first order in
    `genes.py`).
- **Resolution** (verified).
  - Every locus tag `YALI1_…` returned exactly one NCBI GeneID whose record lists that tag. The NCBI taxid is 4952.
  - That GeneID equals the model's `ncbigene` wherever the model has one (688 genes; 2 have none). The model lacks
    `kegg.genes` for 7.
  - KEGG `conv/uniprot` gave one accession for 679 GeneIDs and none for 11. Every linked accession is an active
    CLIB122 entry, and UniProt names its YALI0 ORF.
  - The YALI0 names come from the UniProt record, not from the locus number. Only one replacement is Swiss-Prot:
    YALI1E22238g → [Q05493](https://www.uniprot.org/uniprotkb/Q05493/entry), peroxisomal 3-ketoacyl-CoA thiolase.
- **Sequence identity** (verified, UniParc). 544 replacements are byte-identical to the deleted entry. 135 differ.
  For those, no live UniProtKB entry has the deleted sequence (UniParc membership check). The replacement is the
  CLIB122 protein that NCBI and KEGG attach to the same gene. Its sequence is not the W29 one.
- **No KEGG link** (11 genes).
  - YALI1D32594g (second ACP) and YALI1E31685g (orotidine-5'-phosphate decarboxylase) have identical-sequence live
    entries [Q6C7X2](https://www.uniprot.org/uniprotkb/Q6C7X2/entry) (TrEMBL) and
    [Q12724](https://www.uniprot.org/uniprotkb/Q12724/entry) (Swiss-Prot, URA3).
  - Nine have no live entry: YALI1E25018g, YALI1F36966g, YALI1D15063g, YALI1F38629g, YALI1C23511g, YALI1B26112g,
    YALI1E38859g, YALI1A20558g, YALI1E39248g.
- **W29 tag audit.** 688 of the 690 genes are `ok`. YALI1B26112g and YALI1D30866g are `no_overlap`, and their
  decisions keep the RefSeq protein. For any sequence work, the W29 audit protein stays the source (root
  AGENTS.md §4), not these UniProt entries.

The table gives, per gene, the deleted and replacement accessions, route, entry type, protein name, ORF name, both
UniParc IDs and the W29 flag. The four known genes (verified; identities as in the sibling review):

| Gene | Deleted | GeneID / KEGG | Replacement | Same sequence |
| --- | --- | --- | --- | --- |
| YALI1F28092g | A0A1H6Q312 | 2908086 / yli:2908086 | [Q6C0W3](https://www.uniprot.org/uniprotkb/Q6C0W3/entry), "Monothiol glutaredoxin-5, mitochondrial" | yes |
| YALI1C27933g | A0A1H6PZC1 | 2909892 / yli:2909892 | [Q6CBC8](https://www.uniprot.org/uniprotkb/Q6CBC8/entry), "YALI0C19965p" | yes |
| YALI1D05474g | Q6L8F1 | 2911070 / yli:2911070 | [F2Z6C9](https://www.uniprot.org/uniprotkb/F2Z6C9/entry), "NADPH--cytochrome P450 reductase" | yes |
| YALI1D15063g | A0A1D8NE95 | 2910613 / yli:2910613 | none (KEGG has no UniProt link; no identical-sequence entry) | n/a |

Reaction annotations are not affected by this: Strategy C uses gene `ec-code`, not the accession. The EC numbers
were fetched for these accessions before the deletion (inferred; `gene_ec_cache.json`).

## 4. Code-improvement proposal (not implemented)

Nothing below is authorized. Each change alters model annotation text, and therefore the SHA, so it needs scoped
authorization and goes through the builder.

1. **One post-check for every annotation decision.** All strategies, both passes. Accept an MNXR only when at
   least one *specific* compound (not H+/H2O, cofactor or carrier) is shared, and, when every model compound has an
   ID, the specific-compound Jaccard is ≥ 0.3 (the threshold Layer 2 already uses). This replaces the bare
   `if len(candidates) == 1: mnxr_id = candidates[0]` at `annotate_reactions_extended.py:574–576`. It also guards
   Layer 3 (635–656), primary B (`reactions.py`) and B6. Put the cofactor/carrier list in curated data, not a
   constant.
2. **Never fingerprint a partial reaction.** When any non-H+/H2O metabolite lacks an MNXM, exact-fingerprint lookups
   (primary C, Layer 1, B5) should refuse. Jaccard (Layer 2, D) should count the unannotated compounds in the union.
   This removes the stub hits (`acetyl-CoA = CoA`, `PPi = IPP`, H+ only).
3. **Keep a reaction's xref bundle coherent.** When the crosswalk rewrites `metanetx.reaction`, rebuild or drop the
   keys that came from the old MNXR (KEGG, Rhea, SEED, MetaCyc, SABIO-RK, VMH). Alternatively, write the crosswalk
   MNXR before the annotator runs.
4. **Fix the `_UBIQUITOUS_MNXM` comment.** Line 27 says "MNXM5=NAD+/NADH, MNXM9=ADP". `chem_prop.tsv` names MNXM5
   NADP(+) and MNXM9 phosphate (verified). So the Strategy C/D stop list drops NADP+ and Pi but not NAD+, NADH or ADP.
   Whether that was intended is unknown.
5. **Record the strategy.** Write the strategy and line into the build log per reaction, so a later audit does not
   need line tracing.
6. **Gene UniProt (Q2).** Prefer the reference proteome UP000001300. Refresh or version the proteome cache, and
   record the replacements of `deleted_uniprot_xrefs.csv` as curated data. Leave the 135 different-sequence
   replacements flagged as "CLIB122 entry, sequence differs from W29".

**Test design** (in-memory cobra models with small synthetic `reac_xref` / `reac_prop` dicts, no network):
- T1: a reaction whose gene EC maps to one MNXR with a disjoint fingerprint is left unannotated. This fails on
  today's code, as R1711 shows.
- T2: a single candidate sharing a specific compound with Jaccard ≥ 0.3 is still annotated (no regression).
- T3: Layer 3: several candidates, the name matches one whose fingerprint is disjoint → no annotation.
- T4: an exact fingerprint of annotated metabolites only (acyl-CoA without an MNXM, plus CoA and acetyl-CoA) does not
  accept the `acetyl-CoA = CoA` stub.
- T5: primary Strategy B: a name match with a disjoint fingerprint is rejected.
- T6: crosswalk rewrite: after `write_his` on `metanetx.reaction`, the KEGG/Rhea keys equal the new MNXR's or are
  absent.
- T7, model regression: run `replay_annotation.py` with the patched code and diff against this report.
  - None of the 97 likely-wrong rows keeps its MNXR.
  - The 1,192 `match` rows keep theirs.
  - Every other change is listed.
  - Whole-model diff: annotation keys only, no GPR, stoichiometry or bound change.

Annotations do not enter the LP. So FBA, FVA and every knockout call are unchanged by such a fix; that is provable
from the model structure, and no essentiality rerun is needed. Run pytest before and after, and have an independent
reviewer check the diff (platform/AGENTS.md).

## 5. Verification levels

- **Provable.** Annotation edits cannot change any FBA/FVA/knockout result, because annotations are not part of the
  optimization problem.
- **Statically verified this session.**
  - Fingerprint classes for all 1,632 reactions.
  - The strategy line for every annotation made through Priority 4e (replay, with the fidelity numbers in §2.3).
  - The crosswalk-rewrite mechanism (code read plus 162 exact before/after matches).
  - UniProt status of all 835 gene accessions; NCBI/KEGG/UniProt resolution and UniParc identity for all 690
    deleted ones.
  - Cofactor IDs against `chem_prop.tsv`.
- **Inferred.**
  - The adjudication label of each of the 149 strong rows. Labels are in `scripts/build_tables.py`, one place, for
    review.
  - That the 198 partial rows are mostly the same reaction with ID gaps; not checked row by row.
  - That the offline cache is why deleted accessions persist.
- **Untested.**
  - Any fix.
  - Whether KEGG/Rhea sites agree with MetaNetX's mapping for each ID. KEGG and Rhea were not queried; MetaNetX
    4.5's `reac_xref.tsv` was the bridge.
  - Gene-level function claims: no BLAST, localisation predictor or AlphaFold was run. Replacement names are
    UniProt/KEGG annotations, not identity evidence.

## 6. Limits and found, not checked

- **Limits.**
  - The match rules decide which rows reach the adjudicated set. A stricter compound identity would move rows from
    partial to disjoint, a looser one the other way. The 198 partial rows include some wrong-substrate pairs. R2281,
    for example, shares 3-methyl-2-oxopentanoate but is matched to a ThPP decarboxylation step.
  - MetaNetX 4.5 local tables only. UniProt, NCBI and KEGG are as of 2026-11-12.
- **Found, not checked.**
  - R1499 and R1500 carry 4 H+ on the substrate side of a hydration.
  - R2032/R2033 are named "via proton antiport" but have no H+.
  - R938 moves UTP, UMP, FMN and H+ together.
  - R2170 converts cytosolic "triglyceride" into ER "Triacylglycerol", two different species.
  - m2057–m2060 have no name or annotation; R2282 lacks its electron acceptor, and R2284 has no GPR.
  - MNXR03 (PMF pseudo-compound) is written by a crosswalk decision onto six H+ diffusion reactions.
  - Several "generic" MNXRs carry odd coefficients: 10 Pi on R2177/R2181, 100 Pi on R2145.

**Gap classes.**
- Wrong reaction xrefs and deleted UniProt xrefs *limit a claim*: any statement that uses these IDs as identity or
  provenance. They block no current model task (annotations do not enter the LP).
- The crosswalk staleness *affects another line* (the BiGG crosswalk handoff).
