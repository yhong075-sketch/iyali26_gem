# TASK: identity of YALI1D02625g, the OR partner on transketolase R765/R766 (2026-10-19)

## Authorization and scope

Requested by the session `local_553e1a53` on 2026-10-19 (relayed task, not the user's own words):

- Investigation, then a keep/remove plan for the user's approval.
- **No** model, curation-data or code edits, no build, no LP, no commit.
- An HPCC AlphaFold run for this protein-identity question is pre-authorized by AGENTS.md §2.
  Check the environment first and reuse an existing identical model when one exists.

Question: R765/R766 carry `YALI1E07744g or YALI1D02625g`. YALI1E07744g is the full-length
transketolase. Does YALI1D02625g belong on these reactions as an OR member? Under AGENTS.md §5, an OR
member needs evidence that it catalyses the step on its own.

Recall is not evidence. YALI1E07744g is an experimentally essential gene that the model calls
non-essential (FN) because of this OR. That fact must not drive the decision.

## Starting state

- Worktree `claude/cool-nightingale-a2f471`, HEAD `285eaee` (origin/main merge of PR #6). Clean at
  start.
- Related uncommitted work in other worktrees (read only, not modified here):
  - `claude/clever-kare-7d2e1a`: the builder fix that drops R765/R766 from
    `_GPR_MISANNOTATION_REMOVALS`, and the change-set registry in `selection_gpr_restore.py`.
    Reports: `model/reports/selection_gpr_review_20261019/` and `selection_gpr_fix_20261019/`.
  - `orphan-reaction-newest-model-b8ddc8`: the 2026-10-19 GPR restore (`c98961e9…`).
- Models inspected. All of them have R765/R766 = `YALI1E07744g or YALI1D02625g`, and no other
  reaction uses either gene:

| Model | SHA-256 |
| --- | --- |
| `E5_coq9_resp_pkg_r132_20261012.xml` (this worktree) | `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515` |
| `E5_coq9_resp_pkg_r132_gpr_restore_v2_20261019.xml` (b8ddc8) | `c98961e960c6a0f8283f485e32daf20b04dd70d998fad86fc09dca0d14ccae4c` |
| `E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_20261019.xml` (clever-kare) | `51c6fabd486f4e1a9d6365bf7740bde3eca20fa05f1d5454a73586d3276db3fa` |
| `model/source/iyali26.xml` (immutable source; R765/R766 = `YALI1E07744g` only) | `5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d` |

- Curation inputs that set the OR:

| File | SHA-256 |
| --- | --- |
| `model/curation/gpr/gpr_isozyme_additions.csv` | `162c32c77c82d278513c078cfcb5ce04e2e7c69ba8df9eda74effc557e1ea16f` |
| `model/curation/selection/metadata_reaction_selection.json` | `abfc4a94a03b4ec8da9e3e2fc3ba2e4f04806ef802b6e93c442274065ca0cb24` |

No reference model is designated (`model/candidates/README.md`). No model was executed.

## Computation run (small, identity evidence only)

- **Inputs:**
  - XP_502319.1 (D02625g) and XP_503628.3 (E07744g), NCBI.
  - NC_090773.1:232572–293615, the W29 chromosome D sequence ±30 kb around D02625g, plus its feature
    table.
  - The W29 RefSeq proteome (GCF_001761485.1, `w29.faa` `1863842c…5fcf`) and UniProt Swiss-Prot
    2026_03 (`sprot.fasta` `0e39da39…332b`). Both are reused read-only from
    `$IYALI26_RESEARCH_ROOT/artifacts/tasks/orphan_gene_search_20261019/db/` (hashes match that
    task's `db_provenance.json`), through a space-free scratch symlink.
  - AlphaFold DB AF-Q6CAJ3-F1 and AF-Q6C6T4-F1 (v6).
- **Tools:** BLAST+ 2.17.0 and Biopython 1.88 (existing conda env `seqsearch`).
- **Runs:**
  - pairwise blastp and local alignment;
  - blastp vs W29 and vs Swiss-Prot;
  - tblastn vs the ±30 kb region;
  - one NCBI web blastp vs nr (Fungi, excluding *Y. lipolytica*). It did not return within
    52 min, so there is no result;
  - REST queries to InterPro (Pfam matches, CLIB122), UniRef, OrthoDB and NCBI Datasets;
  - a Cα superposition of the two AlphaFold DB models.
- **Budget:** minutes of CPU and one remote BLAST.
- **Stop condition:** evidence collected for the four sub-questions (domains, locus/split
  annotation, transcript/protein expression, orthologs).
- **AlphaFold:** no HPCC job was submitted. AlphaFold DB already holds a model whose sequence is
  identical to XP_502319.1 (checked residue by residue), so it was reused as AGENTS.md requires. The
  HPCC environment was therefore not checked.

## Outputs

- `REPORT.md`: evidence, tiers, and the keep/remove plan.
- Small evidence tables in this folder.
- Raw downloads and the analysis scripts in
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/d02625g_identity_20261019/`, linked from REPORT.md.
