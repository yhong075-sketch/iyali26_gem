# R24 / R395 homocitrate steps: gene rules, direction, name and protons (review only, 2026-10-29)

**Status: investigation done; waiting for the user's choice of option** ([REPORT.md](REPORT.md) §5). No model,
curated data, builder code or test was changed. Nothing was committed.

## 0. Authorization

- **Brief** from the session that started this one (2026-10-29), the follow-up recorded by the aconitase round
  ([aconitase REPORT §9](../aconitase_aco2_split_20261029/REPORT.md), on branch `claude/aconitase-aco2-split-1c492f`):
  "Investigation only, no model change: review the gene rules of … R24 … and R395 …". It asks which gene(s) should
  catalyse R24; whether R24's direction, name and 3 H+ are right; and what each option does to the knockout classes
  of `YALI1E03387g`, `YALI1D11984g` and `YALI1E17935g`, **in-memory variants only**.
- Computation actually run, inside that scope: 5 in-memory variants × (WT + 4 knockouts) = 25 LPs on one model, one
  medium, one profile; plus read-only model inspection and one pairwise-alignment script. No sweep, no build.
- Not authorized and not done: any model or curation change, a commit, a push.

## 1. Starting state (verified this session)

| Item | Value |
| --- | --- |
| Worktree / branch / HEAD | `claude/nice-bouman-b161f4` @ `4794ec3` (= `origin/main`, PR #38 merge), clean before this folder |
| Baseline model | `e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897`, `model/candidates/E5_…_r673_prolyl_split_aco_split_20261029.xml` in the worktree of branch `claude/aconitase-aco2-split-1c492f` (**uncommitted there; the aconitase round is not merged into `origin/main`**). Chosen per the brief ("until then use `e1e7c868…`"). Not a designated reference. Re-hashed this session on a short copy ([logs/variant_inputs_sha256.txt](logs/variant_inputs_sha256.txt)). |
| Immutable source, for provenance only | `model/source/iyali26.xml`, `5c8c199e…` ([logs/static.log](logs/static.log)) |
| Medium / profile / list | Task A `po1f_csm_leu_w29uptake_20261012.csv` `6eca70fd…`; PO1f `po1f_sd_leu.json` `35307853…`; `consensus_essential_genes.csv` `1e887f5a…` |
| Solver | Gurobi (academic licence, real HOME) at 1e-9, reset before each solve, targeted check — `combined_tpp_heme_r2193_20261022/scripts/checked_solve.py` unchanged (`69c42690…`) |
| Python | `/opt/miniconda3/bin/python` (cobra); alignment with the git-ignored Biopython 1.85 venv of the aconitase worktree |

**Target reactions in `e1e7c868`** ([logs/static.log](logs/static.log), [logs/before_state_dump.txt](logs/before_state_dump.txt)):

| Reaction | Name in model | Equation (C_mi) | Bounds | Rule |
| --- | --- | --- | --- | --- |
| R24 | 2-methylcitrate dehydratase | homocitrate (m75, C7H10O7, 0) ⇌ H2O + 3 H+ + m76 (C7H5O6, −3) | (−1000, 1000) | `YALI1E03387g` |
| R395 | homoacontinate hydratase | H2O + 3 H+ + m76 ⇌ homoisocitrate (m624, C7H10O7, 0) | (−1000, 1000) | `YALI1E03387g` |

Neighbours: R1170 homocitrate synthase (C_mi, irreversible, `YALI1F38776g`), R398 homoisocitrate dehydrogenase (C_mi,
irreversible, `YALI1D13238g`). m75 is used only by R1170 and R24; m76 only by R24 and R395; m624 only by R395 and R398.

**Knockout classes in `e1e7c868`** (re-solved this session, variant A in [variants.json](variants.json)):
`YALI1E03387g` KO/WT 0.0694 (FN at 1/5 %, TP at 10/15 %); `YALI1D11984g` 0.0244 (FN at 1 %, TP at 5/10/15 %);
`YALI1E17935g` 1.0 (TN). These match the aconitase round's sweep.

## 2. Questions and what answers them

1. Which gene(s) catalyse R24 → literature (Fazius 2012 full text, Jia 2006), UniProt, homology, one specificity
   residue.
2. R24 direction, name, 3 H+ → KEGG / Rhea / ChEBI / MetaNetX records, source vs candidate model, network structure.
3. Effect of each option on the three genes' classes → in-memory variants of R24's rule
   ([scripts/r24_variants.py](scripts/r24_variants.py)).

## 3. Stop condition

REPORT.md with evidence tiers, options with effects and acceptance criteria; then stop for the user's choice. Any
chosen option is a separate `curation-change` round needing its own authorization.
