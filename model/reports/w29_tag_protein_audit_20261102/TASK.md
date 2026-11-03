# W29 tag-protein guard (stage 1: strand and position), 2026-11-02

## Authorization

The user, 2026-11-02: "Approve stage 1" for the plan given in the OR-partner review session. Stage 1 is a read-only
audit tool, a generated audit CSV, a curated decisions CSV, a test, one AGENTS.md line and a STATE entry, on a new
branch from origin/main. Stage 2 (the CLIB122-vs-W29 sequence-identity check, which needs a proteome download) is
not authorized.

Not in scope: any model, GPR or builder change; `curated_locus_identity_exclusions.csv` and `LocusCrosswalk`
(an exclusion there blocks the YALI0 alias, so `YALI1D32754g` would lose UniProt O74936 on the next build); the CI
workflow; commits without an explicit instruction.

Computation: local tblastn of 5 CLIB122 proteins against the W29 genome (run 10). Reused runs 5 and 9 of
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/or_partner_review_20261102/`.

## Starting state

- Branch `claude/w29-tag-audit` from `origin/main` at `80cfc4f`, clean.
- Inputs: S2 table `reference/locus_map/S2_table_YALI1_YALI0_map.xlsx` (SHA-256 `42658fb9…`); NCBI feature table
  `reference/ncbi/clib89_feature_table.txt`, GCF_001761485.1 (SHA-256 `15a15c73…`); full hashes in
  `model/curation/annotation/w29_tag_protein_audit.json`.
- Gene set: the start model plus all 45 candidate XMLs under `model/candidates/`, 1,082 gene IDs.
