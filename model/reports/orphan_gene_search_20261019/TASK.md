# Gene search for orphan reactions in triage groups 10, 7 and 6 (2026-10-19)

Status: **confirmed by the user 2026-10-19**, before any computation.

## Authorization

Requested by the session that ran the orphan inventory
([orphan_reactions_20261019](../../../../orphan-reaction-newest-model-b8ddc8/model/reports/orphan_reactions_20261019/REPORT.md),
uncommitted in worktree `orphan-reaction-newest-model-b8ddc8`): run steps 4–5 of the
`gem-orphan-reaction-curation` skill for triage groups 10, 7 and 6. AGENTS.md section 2 requires
the user to confirm this computation (inputs, scope, budget, outputs, stop condition) before it
starts.

User answers in this session (2026-10-19), to three questions on this draft:
- Research root: "Use ~/.zshrc path" →
  `IYALI26_RESEARCH_ROOT=/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`,
  exported per command.
- Budget: "Approve as written".
- Group 7 predictors: "Yes, submit the forms" (up to 15 public RefSeq sequences, one batch each
  to TargetP-2.0, MitoFates and DeepLoc web forms through the built-in browser).

Out of scope: any model XML, curation data, GPR, chemistry, direction or medium change; flux
runs; the CoQ9 line (R19, R52, R321, R1043; approval gate in force); triage groups 1–5, 8, 9;
commits and pushes.

## Starting state

- Worktree `claude/unruffled-shannon-2baf44` at `285eaeee0e403f1c971267898eaed1a884edb727`
  (`origin/main`).
- Model (read-only, used only for reaction equations and gene lists):
  `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml`,
  SHA-256 `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515`. Newest candidate,
  not a designated reference (`model/candidates/README.md`).
- Reaction list: `triage.csv` of the inventory report,
  SHA-256 `5305f2a6d1d142af9a07b6ff0370657a6beccc429e3fa97d2fc8184730c394dc`
  (with `orphan_evidence.csv` `a6951f4c…`, `REPORT.md` `e2136622…`).
- Skill code: `gem-orphan-reaction-curation` `kernel.py`
  SHA-256 `19fe8cfa3ff7b7eb9491ee90699799876d176082b561a80e050f1b0fe263988f`
  (`build_curation_dbs`, `blast_w29`, `reciprocal_hits`, `uniprot_queries`, `uniprot_info`).
- `IYALI26_RESEARCH_ROOT` is **unset** in this session. `~/.zshrc` line 30 exports
  `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`, the root Task B used.
  It is used only if the user confirms it.
- BLAST+ is not installed (no `blastp`; no conda env has it).

## Scope: 33 reactions

| Group | n | Reactions |
| --- | --- | --- |
| 10 can carry flux, no other flag | 25 | R23, R333, R517, R721, R984, R1048, R1078, R1080, R1090, R1132, R1154, R1178, R1252, R1273, R1274, R1338, R1347, R1958, R1959, R1973, R2003, R2010, R2013, R2075, R2076 |
| 7 compartment copy | 6 | R989, R1024, R1027, R1125, R1292, R1337 |
| 6 reverse flux at optimum | 2 | R1308, R1901 |

Prior evidence reused, not redone: R989 has a human-approved provisional assignment of
YALI1D01489g on the lipid-unlump line (`model/curation/lipid_unlump/lipid_unlump_sn_core_curation.json`,
`R989_gpr_curation`: AlphaFold AF-Q6CAN2-F1 v6, TargetP-2.0 Other 0.8307 / mTP 0.1639,
job `6A8DF2140014D6282BEC5B97`). That line is owned by a colleague; this task reviews only.

## Method per reaction

1. Characterised query protein(s) for the exact reaction: UniProt reviewed entries, with an
   activity annotation or a cited assay. Order: *S. cerevisiae*, other yeasts and fungi, then any
   organism. Check substrate, cofactor and direction against the model equation.
2. `blast_w29` (blastp, E ≤ 1e-5, 20 targets) against the W29 RefSeq proteome.
3. `reciprocal_hits` (top 5 Swiss-Prot hits) for each candidate. Accept only if the best hits are
   the same enzyme/EC and no closer paralogue of a different function exists.
4. Record separately: (a) the enzyme exists; (b) its compartment and direction match the
   reaction; (c) the substrate has a source in the model (static check of producers only, no
   flux run).
5. Evidence tier (skill definitions): A = *Y. lipolytica* gene characterised for the activity;
   B = reciprocal-best-hit orthologue of a characterised (preferably fungal) enzyme, about
   ≥ 30 % identity and ≥ 70 % coverage; C = family-level, paralogue-ambiguous or distant-taxon
   (tentative); D = no defensible gene → `keep_orphan`.
6. Every gene is named with its systematic ID, a verified name or "none verified", a short
   function and the evidence level. KEGG KO links are leads only.
7. Literature: *Y. lipolytica* first, then *S. cerevisiae*/fungi. Only sources opened this
   session are cited, each marked as opened.

**Group 7 targeting check** (as in Task B, `r132_mito_adh_20261012`):
- Proteins: the gene-bearing paralogue of each copy (YALI1D01489g, YALI1D33889g, YALI1D24649g,
  YALI1F28274g, YALI1C09285g) plus any other W29 paralogue that step 3 returns for the same
  activity. At most 15 sequences.
- Predictors: official web forms of TargetP-2.0 (Non-Plant), MitoFates (Fungi) and DeepLoc 2.x
  (DTU), full public RefSeq sequences, no account and no email, submitted through the built-in
  browser. Each form submission needs the user's explicit permission.
- Decision rule, fixed now: a protein counts as predicted mitochondrial if TargetP-2.0
  mTP ≥ 0.5 **or** MitoFates probability ≥ 0.385. For the extracellular copy R1027, TargetP SP
  ≥ 0.5 counts as a predicted secretion signal. DeepLoc is recorded but not part of the rule.
  Scores are reported as they are, without rounding to a threshold.
- A prediction alone never assigns a GPR across compartments. It is recorded as support or
  counter-evidence next to any localization source that is opened.

**Group 6:** the gene search is the same. The direction question (R1901, R1308) is noted but not
decided here; it needs enzymology or a ΔG source.

## Budget

- Software: one new conda env `seqsearch` (conda-forge, bioconda: `blast`, `biopython`, `pandas`).
- Downloads, once: W29 RefSeq proteome and GFF, assembly GCF_001761485.1 (NCBI Datasets API),
  and UniProt Swiss-Prot FASTA (about 90 MB compressed). One `makeblastdb` each.
- BLAST: for these 33 reactions only. At most 5 query proteins per reaction against W29, and
  reciprocal checks for at most 5 W29 candidates per reaction against Swiss-Prot.
- Targeting: one batch per predictor, at most 15 sequences.
- Web lookups: UniProt REST, KEGG REST (leads only), NCBI protein pages, literature.
- No flux simulation, no model build, no AlphaFold job.

## Outputs

In `model/reports/orphan_gene_search_20261019/`:
- `TASK.md` (this file);
- `gene_search.csv`: one row per reaction with `id, group, proposed_gpr_or_keep_orphan,
  evidence_tier, query_proteins, best_w29_hit, reciprocal_result, paralog_note,
  enzyme_exists, compartment_direction, substrate_source, sources_opened`;
- `blast_w29_hits.csv` and `reciprocal_hits.csv` (raw BLAST tables, small);
- `targeting_group7.tsv` (raw predictor outputs and job IDs);
- `scripts/` (code written for this task);
- `REPORT.md`;
- a dated entry in `model/STATE.md`.

In `$IYALI26_RESEARCH_ROOT/artifacts/tasks/orphan_gene_search_20261019/`: the BLAST databases,
the downloaded FASTA/GFF and `w29_gene_map.csv`, with SHA-256 and download dates recorded in
REPORT.md.

## Stop conditions

- **Done:** all 33 rows filled (a proposed GPR or `keep_orphan`), REPORT.md and the STATE.md
  entry written.
- **Stop and report, no retry or substitution:**
  - the research root is unset or unconfirmed;
  - the model SHA changes;
  - the database build fails, or the W29 proteome deviates from the expected size
    (7,918 proteins with YALI1_ locus tags, per the skill);
  - a predictor server is unavailable (record the gap; no other tool without asking);
  - the work would extend beyond these 33 reactions or into any model or data change.
