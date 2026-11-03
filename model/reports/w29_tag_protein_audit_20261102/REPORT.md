# W29 tag-protein guard, stage 1: strand and position (2026-11-02)

Scope: [TASK.md](TASK.md). No model, GPR or builder change; nothing committed.

## 1. What the guard checks

The YALI1 IDs come from the S2 table of the W29/CLIB89 assembly (Magnan et al. 2016, PLoS ONE 11:e0162363). The paper
reassigned the CLIB122 YALI0 names to the YALI1 genes at the same position, and the S2 rows carry those positions
and strands. For some tags, however, the CDS that GenBank (2016) and RefSeq (GCF_001761485.1) file under the tag sits
on the other strand. Tools that read the protein by tag then get an antisense ORF.

[`platform/tools/audit_w29_tag_protein.py`](../../../platform/tools/audit_w29_tag_protein.py) compares, for every
gene in the start model and all 45 candidate XMLs (1,082 IDs, read as text), the S2 row with the RefSeq CDS row of
`reference/ncbi/clib89_feature_table.txt`. Each gene gets one flag. The tests run in this order: not a YALI1 ID → no S2
row → no RefSeq CDS (then the S2 gene type and chromosome decide the flag) → strand → position. Any tag that has a
RefSeq CDS is compared, whatever its S2 gene type. "Same position" means at least 1 bp of overlap on the same
chromosome. That is lenient. 14 of the 1,044 `ok` rows overlap by less than half of the longer span (overlap ÷
longer span, inclusive coordinates; the smallest is YALI1D32594g at 18 %). The reviewer reported 25 with a definition
that was not recorded.

| Flag | Meaning | Decision needed | Genes |
| --- | --- | --- | ---: |
| `ok` | same strand, overlapping spans | no | 1,044 |
| `strand_flip` | the RefSeq CDS is on the other strand | yes | 7 |
| `no_overlap` | same strand, but other chromosome or disjoint spans | yes | 2 |
| `no_s2_row` | YALI1-form ID with no S2 row | yes | 1 |
| `no_refseq_cds` | nuclear S2 mRNA with no RefSeq CDS | yes | 0 |
| `mito_absent_from_feature_table` | YALI1M mRNA; the RefSeq assembly has no mitochondrial CDS rows | no | 12 |
| `not_protein_coding` | S2 gene type is not mRNA | no | 0 |
| `not_yali1` | YALI0, YALIUNK or old mitochondrial names | no | 16 |

Outputs:
- [`w29_tag_protein_audit.csv`](../../curation/annotation/w29_tag_protein_audit.csv) (generated);
- [`w29_tag_protein_audit.json`](../../curation/annotation/w29_tag_protein_audit.json) (input SHA-256s, gene-set
  fingerprint, flag counts);
- [`w29_tag_protein_decisions.csv`](../../curation/annotation/w29_tag_protein_decisions.csv) (curated, one row per
  flagged gene).

Each decision row records the RefSeq protein and span it was made for. If a later feature table files a different
protein under the tag, the decision no longer counts.

[`platform/tests/test_w29_tag_protein_audit.py`](../../../platform/tests/test_w29_tag_protein_audit.py) fails when:
- the model gene set changes without rerunning the tool;
- a flagged gene has no valid decision, or a decision no longer matches its flag, protein or span;
- the raw inputs no longer reproduce the committed files (this part errors when the research workspace is unset).

## 2. The 10 decisions (evidence verified this session)

- **Strand flips (7): 6 `refseq_protein_is_antisense_orf`, 1 `unresolved_conflict`.** Each RefSeq protein has no
  Swiss-Prot hit at E ≤ 1e-5
  ([evidence/run9_flagged_vs_sprot.tsv](evidence/run9_flagged_vs_sprot.tsv); the R540/R712 pair is in the OR-partner
  review). For each gene, tblastn places the CLIB122 UniProt protein on the S2 strand at the S2 span
  ([run 5](evidence/run5_clib122_tblastn_w29genome.tsv), [run 10](evidence/run10_flip5_tblastn_w29genome.tsv)):

  | Gene | Reactions | CLIB122 protein | tblastn |
  | --- | --- | --- | --- |
  | YALI1D32754g | R97-R102, R296, R1487-R1489 | O74936 POX3 acyl-CoA oxidase 3 | 100 % over 1-687 of 700 |
  | YALI1D32550g | R2062 | Q6C7X4 (YALI0D24585g) | 100 % full length |
  | YALI1E27218g | R2062 | Q6C4W9 (YALI0E23089g) | two HSPs, 96.6 % and 92.9 %; not identical |
  | YALI1F31153g | R1975, R1977 | Q6C0L1 glutaryl-CoA dehydrogenase | 100 % full length |
  | YALI1F15367g | R2195 | Q6C217 cysteine dioxygenase | 100 % full length |
  | YALI1E18171g | R540 | Q6C5U4 (GMPPA-like) | 100 % full length |
  | YALI1F02518g | R712 | Q6C393 EYI1 | 100 % full length |

  The GenBank records for the last two confirm the complement-strand CDS
  ([evidence/AOW_genbank.gp](evidence/AOW_genbank.gp)).

  "Antisense ORF" is inferred from three facts: no homology, the opposite strand to the S2 row, and an intact,
  conserved CLIB122 ORF on the S2 strand. Whether the minus strand is transcribed was not checked.
- **`YALI1E18171g` is an evidence conflict, left to the user.**
  - `curation/essentiality/curated_locus_identity_exclusions.csv` row 1 (active, case EGC-41703ef278c1, source
    PMC5014426) calls the minus-strand protein AOW05442/XP_068139108 "RNA-seq-supported". It says Q6C5U4 is a
    different protein, and that this crosswalk "must not propagate identity or function".
  - This audit agrees that the two proteins differ. It does not agree on which one the gene means: the minus-strand
    protein has no homology, while Q6C5U4 is intact on the S2 strand.
  - The decision is therefore `unresolved_conflict` and names no protein.
  - The RNA-seq claim was not re-read here; the cited dossier JSON has no RNA-seq text (reviewer finding).
  - Resolving the conflict means choosing between the two records, and may mean editing the exclusion row, which
    changes the crosswalk.
- **No overlap, `refseq_protein_kept` (2).** The S2 rows are short spans (299 and 405 bp) next to the RefSeq CDS.
  The RefSeq proteins are the right enzymes:
  - YALI1B26112g: glutamate synthase, GLT1_SCHPO 66.7 %, for R331;
  - YALI1D30866g: BCKDH E2, ODB2_CAEEL 42.8 %, for R1968-R1970 and R2281.
- **`YALI1M00338r`, `model_id_differs_from_s2_tag`.** S2 names the tag `YALI1_M00338g` (nad5). STATE 2026-10-08
  records the deliberate decision to keep this ID as ND5.

## 3. What the guard does not do

- It does not change the model or the crosswalk. Putting the 7 strand flips into
  `curated_locus_identity_exclusions.csv` would block their YALI0 aliases, and `YALI1D32754g` would lose UniProt
  O74936 on the next build. That is a model change outside this scope.
- It checks strand and position only. A same-strand ORF with a wrong frame or start, or a truncated gene model (for
  example FKS1 `YALI1E25018g`, which starts at Q9UVK9 residue 235), passes. That is stage 2 (sequence identity),
  not authorized.
- It is not in CI. The evidence job installs no openpyxl and has no research workspace; the repo-only tests would
  run there if added, but changing the CI scope is a separate decision.
- The 12 mitochondrial genes have no RefSeq protein in this assembly, so no W29 protein is checked for them.

## 4. Verification (this session)

- The tool reproduces the earlier scan of 2bc263f2: 1,044 ok, 7 strand flips, 2 no-overlap. The earlier 28
  "no S2 row" split into 16 non-YALI1 IDs, 12 mitochondrial genes and `YALI1M00338r`. The union of all model files
  adds `YALIUNK2`.
- New test file (venv `.claude/worktrees/pytest-performance-1daa1e/platform/.venv`): 12 passed, 23 subtests passed.
  - Synthetic rows exercise every `check_decisions` rule: missing, duplicate, schema, decision vs flag, sequence
    source, accession, antisense → UniProt, kept → the audited RefSeq protein, unresolved → no protein, evidence and
    reason, flag mismatch, changed protein or span, and stale decisions.
- Guard checks:
  - an empty decisions file gives 10 errors, one per flagged gene;
  - an unset research root raises `RuntimeError` in the rebuild;
  - the CLI exits 1 on any decision gap, and a rerun leaves both audit files byte-identical.
- Independent read-only review (general-purpose agent following `.claude/agents/code-reviewer.md`): 2 medium and
  5 low findings, all fixed or documented:
  - decisions not tied to the audited protein (fixed);
  - the E18171g evidence conflict (now `unresolved_conflict`);
  - non-mRNA S2 rows skipping the comparison (fixed; 0 model genes affected);
  - the lenient overlap criterion (stated in §1);
  - untested decision rules (now tested);
  - CLI exit code (fixed);
  - evidence paths and coordinates (made consistent).
- Full suite: see §5.

## 5. Full pytest (verified this session)

Venv `.claude/worktrees/pytest-performance-1daa1e/platform/.venv` (scipy 1.18.1, cobra 0.30.0), run from `platform/`
with `-q -rfE -p no:cacheprovider`. The root `.venv` could not collect 4 gem_review test files (no `pydantic_ai`).

| Run | Failed | Errors | Passed | Skipped |
| --- | ---: | ---: | ---: | ---: |
| Before: this tree with `--ignore=tests/test_w29_tag_protein_audit.py` | 3 | 14 | 457 | 24 |
| After: full | 3 | 14 | 469 | 24 |

- **Test by test (JUnit XML):** the only differences are the 12 new tests, all passed. No existing test changed
  status.
- **Standing failures, unchanged:**
  - `test_dfba_new_fn_essentiality`, `test_lipid_moiety_ledger`, `test_vacuole_supply` (1 failure each);
  - 14 `test_lp_sn12_candidate` errors.
  - These are the counts PR #43 records for this venv.
- **Why the "before" run is a fair baseline:** this change modifies no existing code file, only `AGENTS.md`, so the
  baseline equals origin/main for every other test. It was not run on a separate origin/main checkout.
- Logs: [logs/pytest_before.txt](logs/pytest_before.txt), [logs/pytest_after.txt](logs/pytest_after.txt).
