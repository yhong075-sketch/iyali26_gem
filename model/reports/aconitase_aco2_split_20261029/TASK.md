# Aconitase rules R1556, R216, R238: ACO2-type YALI1E17935g removed, YALI1D11984g kept (2026-10-29)

**Status: approved 2026-10-29** (§0). Written before approval: this file,
[before.json](before.json) (read-only dump of the baseline), [scripts/before_state.py](scripts/before_state.py),
[scripts/align_aconitases.py](scripts/align_aconitases.py) and [evidence/](evidence/) (UniProt and PubMed records
fetched this session and the alignment table, with [SHA256SUMS.txt](evidence/SHA256SUMS.txt)). Outside this folder:
`/.tools/` added to `.gitignore` (git-ignored alignment venv).

## 0. Authorization

- **Brief** from the session that started this one (2026-10-29): run one `curation-change` round removing
  `YALI1E17935g` from R1556, R216 and R238 and leaving `YALI1D11984g` as their only gene; draft TASK.md and
  before.json, then stop for the user's approval. This is the lever ranked first by a read-only review that same day.
- **User, 2026-10-29**, after the alignment was added (the user's instruction "you can download it, just remember to
  add in gitignore" authorized the Biopython install):

  > "Approve the aconitase plan in TASK.md §7 on baseline 059316ef."

  No commit or push was authorized.

## 1. Starting state (verified this session)

| Item | Value |
| --- | --- |
| Worktree / branch / HEAD | `claude/aconitase-aco2-split-1c492f` @ `dcefab7` (= `origin/main`, PR #37 merge), clean before this folder |
| Baseline | `candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_notes_v3_20261029.xml`, SHA-256 `059316ef34417c88aadba61cfd5606134d27879a842b3835d2cc7df917a9248d` (re-hashed). **Not a designated reference** (`model/candidates/README.md`); it is the newest candidate of this line. |
| Its build | `PYTHONHASHSEED=0 /opt/miniconda3/bin/python -m scripts.gem_annotate <argv>` (Python 3.13.5), argv in its `.build.json`; seed 1 gave identical bytes ([pass2 build_commands](../r171_r673_note_revision_20261029/pass2/logs/build_commands.txt)) |
| Platform code `efdde43` → `dcefab7` | no change (`git diff --stat efdde43 dcefab7 -- platform/` is empty) |
| Experimental list | `consensus_essential_genes.csv`, SHA-256 `1e887f5a…` (re-hashed). `YALI1D11984g` is listed (line 257, "highly similar to P19414 ACO1"); `YALI1E17935g` is not. |

**Model state of the targets** ([before.json](before.json)):

| Reaction | Name in model | Compartment | Bounds | Rule |
| --- | --- | --- | --- | --- |
| R238 | citrate to cis-aconitate(3-) | C_mi | (−1000, 1000) | `YALI1D11984g or YALI1E17935g` |
| R216 | cis-aconitate(3-) to isocitrate | C_mi | (−1000, 1000) | `YALI1D11984g or YALI1E17935g` |
| R1556 | cis-aconitate(3-) to isocitrate | C_cy | (−1000, 1000) | `YALI1D11984g or YALI1E17935g` |
| R239 | citrate to cis-aconitate(3-), cytoplasmic | C_cy | (−1000, 1000) | `YALI1D11984g` (not touched) |

- `YALI1D11984g` sits on R1556, R216, R238, R239; `YALI1E17935g` only on R1556, R216, R238. No other reaction
  carries either gene, and the three target rules contain no other gene.
- Other aconitase-family reactions, context only: R24 (homocitrate ⇌ homoaconitate, named "2-methylcitrate
  dehydratase") and R395 (homoaconitate ⇌ homoisocitrate), both `YALI1E03387g` (LYS4-like, consensus-essential);
  R95 (2-methylcitrate ⇌ 2-methylisocitrate, named "aconitase"), `YALI1F03803g`; R735/R736 spontaneous, no gene.
- Current knockout class (r673 round sweep on `6af03c40`, the same optimisation problem as `059316ef` per STATE.md;
  [ko_compare.tsv](../r673_prolyl_split_20261029/ko_compare.tsv)): `YALI1D11984g` KO growth 0.311661 = WT →
  non-essential at every cutoff (**FN**); `YALI1E17935g` the same (**TN**). Recalled from that file, not re-solved.

## 2. Identity of the two genes

| | `YALI1D11984g` | `YALI1E17935g` |
| --- | --- | --- |
| CLIB122 locus (locus map, verified) | `YALI0D09361g` | `YALI0E14949g` |
| RefSeq in model = UniProt cross-ref (verified) | XP_502616.2 = [Q6C9P6](https://rest.uniprot.org/uniprotkb/Q6C9P6.json) | XP_503960.1 = [Q6C5V2](https://rest.uniprot.org/uniprotkb/Q6C5V2.json) |
| UniProt name (unreviewed, automatic) | "Aconitate hydratase, mitochondrial", 779 aa | "Aconitate hydratase, mitochondrial", 790 aa |
| PANTHER subfamily (UniProt cross-ref, verified) | PTHR43160:SF3 "ACONITATE HYDRATASE, MITOCHONDRIAL" | PTHR43160:SF2 "HOMOCITRATE DEHYDRATASE, MITOCHONDRIAL" |
| *S. cerevisiae* counterpart | ACO1 [P19414](https://rest.uniprot.org/uniprotkb/P19414.fasta) (list annotation "highly similar to P19414"; alignment below) | ACO2 [P39533](https://rest.uniprot.org/uniprotkb/P39533.fasta) "Homocitrate dehydratase, mitochondrial" (PANTHER; alignment below) |
| Verified name | none (no reviewed *Y. lipolytica* entry); "ACO1" in Holz 2009, locus not given in the abstract | none |
| Evidence level | homology + *Y. lipolytica* overexpression phenotype for "ACO1" (abstract) | homology + *S. cerevisiae* / fungal biochemistry; no *Y. lipolytica* data |

All UniProt facts are from the JSON fetched this session (evidence/, entry versions 113 and 111).

**Pairwise alignment** (run this session, **verified**; [scripts/align_aconitases.py](scripts/align_aconitases.py),
output [evidence/pairwise_alignment.tsv](evidence/pairwise_alignment.tsv)). Biopython 1.85 in a git-ignored venv
`.tools/align-venv` (installed from PyPI on the user's instruction; `/.tools/` added to `.gitignore`); global,
BLOSUM62, gap open −10, extend −0.5, end gaps free (EMBOSS needle settings); CLIB122 UniProt sequences (the
RefSeq IDs in the model are their cross-references). Identity over the aligned span:

| | Sc Aco1p (P19414) | Sc Aco2p (P39533) |
| --- | --- | --- |
| `YALI1D11984g` (Q6C9P6) | **77.7 %** (similarity 87.3 %) | 52.3 % |
| `YALI1E17935g` (Q6C5V2) | 54.9 % | **70.3 %** (similarity 82.0 %) |

For scale: Aco1p vs Aco2p 54.7 %, and the two *Y. lipolytica* proteins 53.7 % to each other. Each *Y. lipolytica*
protein is clearly closest to one yeast paralog, and the pairing is reciprocal, agreeing with PANTHER. This is
homology among four sequences, not a phylogeny, and says nothing directly about catalysis.

## 3. The change

R1556, R216, R238: `YALI1D11984g or YALI1E17935g` → `YALI1D11984g`. Bounds, stoichiometry, R239 and every other
reaction unchanged. `YALI1E17935g` stays in the model with no reaction (precedent: `YALI1B07957g`, `YALI1D19110g`).
Notes on each of the three reactions (curation id, status, authorization, previous rule, unlinked gene, evidence,
tier), following the R673 step. Proposed wording for the unlinked gene:

> `YALI1E17935g` (ACO2-type, PANTHER homocitrate dehydratase subfamily) removed: *S. cerevisiae* Aco2p contributes
> only to lysine biosynthesis and Aco2 homologues of respiring filamentous fungi were enzymatically inactive (Fazius
> 2012); the *Y. lipolytica* protein is untested, so it is taken as homocitrate-specific or inactive, not as a
> citrate aconitase. Not assigned to R24 (separate review).

## 4. Evidence and its tier

**For** the change:
- Fazius et al. 2012, [PMID 23106124](https://pubmed.ncbi.nlm.nih.gov/23106124/), abstract fetched this session
  (**verified, abstract only**): yeast Aco1p is essential for the citric acid cycle, while Aco2p "specifically and
  exclusively contributes to lysine biosynthesis"; Aco2p homologues of filamentous fungi were transcribed but
  enzymatically inactive and did not complement yeast aconitase mutants.
- PANTHER places `YALI1E17935g` in the ACO2 (homocitrate dehydratase) subfamily and `YALI1D11984g` in the ACO1
  subfamily (**verified**, UniProt cross-reference; this is an HMM classification, i.e. homology). The pairwise
  alignment in §2 agrees: 70.3 % identity to Aco2p vs 54.9 % to Aco1p for `YALI1E17935g`, and 77.7 % to Aco1p vs
  52.3 % to Aco2p for `YALI1D11984g` (**verified**, homology only).
- Holz et al. 2009, [PMID 18850095](https://pubmed.ncbi.nlm.nih.gov/18850095/) (abstract, **verified**): multicopy
  *Y. lipolytica* ACO1 shifts citric/isocitric acid towards isocitrate, so the ACO1 gene product is an active citrate
  aconitase in vivo. That this ACO1 is `YALI0D09361g` is **inferred** (the abstract gives no locus).
- AGENTS.md §5: an OR member needs evidence that it catalyses the step on its own. No evidence found that
  `YALI1E17935g` converts citrate or cis-aconitate.

**Against / limits:**
- UniProt names both proteins "Aconitate hydratase, mitochondrial" and gives Q6C5V2 GO "aconitate hydratase activity"
  (IBA, phylogenetic inference) and "tricarboxylic acid cycle" (IEA). Both are automatic.
- *Y. lipolytica* is a respiring yeast: which fungal pattern applies (yeast lysine-specific Aco2p, or inactive
  filamentous homologue) is unknown. Fazius also reports that fungal *aconitases* convert homocitrate, so in
  *Y. lipolytica* the ACO1-type enzyme may also serve lysine synthesis. That bears on R24, not on this change.
- No *Y. lipolytica* enzymology, localisation or knockout for `YALI1E17935g`; 8 PubMed hits for
  "*Yarrowia lipolytica* aconitase" were screened by title and abstract (evidence/yl_aconitase_abstracts.txt), none
  addresses ACO2.
- The cytosolic R1556 / R239 placement of `YALI1D11984g` is inherited and not reviewed here.

## 5. Expected effect (to be tested, not assumed)

- Only `YALI1D11984g` and `YALI1E17935g` can change class: the three rules contain no other gene (statically
  verified from before.json).
- `YALI1D11984g` KO then closes R238, R216, R1556 and R239. The read-only review reported KO/WT **0.0244** for that
  reaction set on an earlier candidate (not re-checked). If that holds on this baseline: FN → TP at **5, 10 and
  15 %**, but **not at 1 %** (0.0244 ≥ 0.01 stays non-essential). The brief's "FN → TP at every cutoff" would then be
  wrong at 1 %. 15 %: 124 / 191 / 42 / 692 → 125 / 190 / 42 / 692.
- `YALI1E17935g` KO: no reaction → growth = WT → stays TN. WT growth unchanged (rules do not enter the WT LP).
- The residual ~2 % growth with all four aconitase reactions closed is not explained yet; the pFBA for the figure
  will show the route (a medium amino acid feeding 2-oxoglutarate is a guess, **inferred**).

## 6. Acceptance criteria

1. Control build with the `059316ef` argv (flag off), `PYTHONHASHSEED=0`, `/opt/miniconda3/bin/python`, reproduces
   `059316ef…` byte for byte.
2. Candidate seeds 0 and 1 give identical bytes.
3. Whole-model diff vs `059316ef`: only the three rules, their notes, and `YALI1E17935g`'s three reaction links.
4. Checked sweep: WT unchanged; `YALI1D11984g` changes as measured (expected §5) and no other gene changes class at
   any cutoff; unresolved set unchanged; no failing point.
5. Full pytest: only the new tests differ from the before run.
6. Independent read-only audit returns no unsupported claim left unfixed or unreported.

## 7. Plan (after approval)

1. **pytest before:** a fresh full run at `dcefab7` (about 13 min, background). The pass2 after-log (336 passed,
   19 failed, 20 errors) was taken on an uncommitted tree, so it is kept only as a cross-check.
2. **Data:** `model/curation/gpr/aconitase_aco2_split.json` (schema as `r673_prolyl_split.json`: before/after rule,
   bounds, stoichiometry and species per reaction, unlinked gene, evidence with URLs and tiers, reviewer).
3. **Step:** `platform/scripts/gem_annotate/aconitase_aco2_split.py`, `apply_aconitase_aco2_split`; fails closed if a
   target's rule, bounds or stoichiometry differs from the data, if either gene is missing, or if `YALI1E17935g`
   sits on any reaction outside the three targets; idempotent. Flag `--aconitase-aco2-split` (no dependency), called
   in `main.py` after the R673 step; option in `cli.py`.
4. **Test:** `platform/tests/test_aconitase_aco2_split.py` on an in-memory model (applies; idempotent; stale rule,
   stale stoichiometry and extra reaction on the gene each fail closed; nothing outside the targets changes).
5. **Builds** (from `platform/`, `PYTHONHASHSEED=0`, `/opt/miniconda3/bin/python`): control to the scratchpad;
   candidate `model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_20261029.xml`
   (file name 145 bytes, `.build.json` sidecar 152, under the 255-byte limit) = `059316ef` argv + `--aconitase-aco2-split`; seed-1 rebuild to scratchpad.
6. **Diff:** `model_diff.py` from the r673 round (copied unchanged).
7. **Sweep:** `ko_compare_two.py` from the r673 round, Task A medium (`6eca70fd…`) + PO1f (`35307853…`), Gurobi
   1e-9 with the real HOME, basis reset per knockout, targeted check; SHA-checked short copies of both models.
   Four cutoffs with recall, specificity, precision and positive coverage; user-defined negatives; development
   reference.
8. **Audit** by a fresh read-only reviewer; then REPORT.md, STATE.md entry, candidates README row.
9. **Figure:** Escher before/after (skill `gem-route-escher-maps`) of the `YALI1D11984g` knockout around the
   aconitase and TCA steps, pFBA on both models.

**Budget:** about 2,164 LPs (sweep) + 20 (fresh re-solves) + about 10 (figure); 45 minutes of compute. **Stop** on
any failed acceptance criterion and report it.

## 8. Not in this round (follow-ups to record)

- **R24 / R395 (lysine, homocitrate → homoisocitrate):** both carry only LYS4-like `YALI1E03387g`, while Fazius 2012
  reports that homoaconitases catalyse only homoaconitate ⇌ homoisocitrate and aconitases the homocitrate step. R24
  is also named "2-methylcitrate dehydratase". Separate review; `YALI1E17935g` is not moved there.
- Compartment of the cytosolic aconitase pair R239 / R1556.
