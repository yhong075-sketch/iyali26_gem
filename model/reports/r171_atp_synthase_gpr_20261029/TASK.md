# R171 ATP synthase rule: single-gene OR branches removed, subunits e and j out of the AND (Option D) (2026-10-29)

## 0. Authorization

**Brief** from the session that started this one (2026-10-29): continue essentiality curation one change at a time
through `curation-change`; item 1 is R171 (candidate 2 of
[or_rule_review_r171_aars_20261029/REPORT.md](../or_rule_review_r171_aars_20261029/REPORT.md) §6–§8), options A, C, D
to be presented with their catalytic evidence and chosen by the user; the FP count is a consequence, never the reason.

**User, 2026-10-29**, answering "Which R171 option should I build (the plan above, on baseline d8d1bac7)?":

> "Option D (Recommended)" (option text: "12-subunit AND: drop g/k OR branches, drop e and j from the rule (strict
> 'loss closes the reaction' reading; j evidence is one S. cerevisiae abstract).")

The question approved the plan below and the baseline `d8d1bac7` in the same answer. No commit or push.

## 1. Starting state (verified this session)

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/essential-gene-r171-pro-7e2477` @ `ff7fc62` (= `claude/essential-gene-prediction-e76be3`, PR #34, not merged), clean tree |
| Baseline | `candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_um7_nuclear_diffusion_r2171_1to1_aars_thr_ser_split_20261029.xml`, SHA-256 `d8d1bac78138f59ba6beefcbddae7f93d0a20e245fea6f5dd9c79760ee535751` = `14f73bda` + `--aars-thr-ser-split`. Not a designated reference; the user confirmed it as the baseline for this round. |
| Its build | `/opt/miniconda3/bin/python`, `PYTHONHASHSEED=0`, argv in its `.build.json` (30 opt-in flags) |
| Before-state | [before.json](before.json) (`scripts/before_state.py`): R171 name, bounds, stoichiometry, species, annotation, notes and rule; the 16 genes' annotation and reactions. Each of the 16 genes is linked to R171 only. |
| R171 rule now | `(YALI1M00260g and YALI1M00267g and YALI1M00318g and YALI1F04433g and YALI1F04110g and YALI1D27826g and YALI1D15734g and YALI1D14639g and YALI1B08885g and YALI1F07384g and YALI1B05364g and YALI0D17490g and YALI1F26967g and YALI1E38131g) or YALI1B28025g or YALI1B15733g` (inherited verbatim from iYli21, review §1) |

## 2. The change (Option D)

R171 becomes the 12-gene AND of the three mtDNA genes (ATP6/8/9 set) and the nine nuclear core subunits:

`YALI1M00260g and YALI1M00267g and YALI1M00318g and YALI1F04433g and YALI1F04110g and YALI1D27826g and YALI1D15734g and
YALI1D14639g and YALI1B08885g and YALI1F07384g and YALI1B05364g and YALI1F26967g`

Removed from the rule, kept in the model with no reaction (precedent `YALI1B07957g`, `YALI1F07747g`):

| Gene | Subunit (UniProt, verified in the review) | Why it leaves the rule (catalytic ground) |
| --- | --- | --- |
| `YALI1B28025g` | g, ATP20 (B5FVB8, 138 aa) | A single accessory subunit cannot synthesise ATP alone (no evidence it does: AGENTS §5), so the OR branch goes. Dimer-specific: the *Y. lipolytica* monomer lacking e, g, k hydrolyses ATP at about the dimer's rate (Hahn 2016, PMID 27373333); not found in the purified *Y. lipolytica* enzyme (Liu 2015, PMID 25759169). Not added as an AND member. |
| `YALI1B15733g` | k, ATP19 (B5FVB3, 70 aa) | As g; also "not found" in the purified enzyme (Liu 2015). |
| `YALI1E38131g` | e, TIM11/ATP21 (B5FVG3, 90 aa) | Dimer-specific (Hahn 2016, as above); *S. cerevisiae* e, g, k "not required for the formation of enzymatically active ATP synthase" (Arnold 1998, PMID 9857174). |
| `YALI0D17490g` | j, ATP18 (Q6C8S0, 62 aa) | *S. cerevisiae* null viable, about half the specific ATPase activity, low ADP/O (Vaillier 1999, PMID 9867878, abstract record): loss impairs, does not abolish, the activity. Strict Boolean reading (an AND member's loss closes the reaction) therefore leaves it out. One *S. cerevisiae* abstract; no *Y. lipolytica* data; not screened. |

The sources above were opened in the 2026-10-29 review (its §7.1, through Europe PMC REST records and full texts);
they are not re-opened in this round. The choice between A, C and D is a judgement on that evidence, made by the
user. Not reasons: the recall or MCC gain, the FP count of each option, Yeast-GEM's 17-gene AND.

The R171 notes record the curation id, status, authorization, previous rule, the four genes kept unlinked and the
evidence. Nothing else changes: stoichiometry, bounds, annotation, other reactions and the gene list stay as they are.
The `YALI0D17490g` → `YALI1D21476g` identifier question is a separate annotation line, not part of this change.

## 3. Expected effect (from the review's in-memory variant `R171_D` on `965ab947`)

- The 9 consensus-essential nuclear core subunits go FN → TP at 15 % only (KO/WT 0.1097 = R171 closed); no class change
  at 1 / 5 / 10 %. The 3 mtDNA genes change ratio (1.0 → 0.1097) but are out of scope. g, k, e, j: KO/WT 1.0, TN.
- 15 %: 114 / 201 / 42 / 692 (baseline `d8d1bac7`, recorded) → 123 / 192 / 42 / 692.
- WT unchanged (R171's rule does not affect the WT LP); 16 unresolved unchanged.

## 4. Acceptance criteria

1. Control build with `d8d1bac7`'s argv at HEAD (code change in place, flag off) reproduces `d8d1bac7…` byte for byte.
2. Candidate built with `PYTHONHASHSEED=0` and `=1` gives identical bytes.
3. Whole-model diff vs `d8d1bac7`: only R171's rule and notes and the reaction links of the four removed genes.
4. Knockout sweep `d8d1bac7` vs candidate: WT unchanged, exactly the transitions of §3 at every cutoff, unresolved set
   unchanged, no failing point.
5. Full pytest: only the new tests differ from the before run, test by test.

## 5. Plan

1. Data `model/curation/gpr/r171_atp_synthase_gpr.json`; step `platform/scripts/gem_annotate/r171_atp_synthase_gpr.py`
   (`apply_r171_atp_synthase_gpr`, pins R171's name, bounds, stoichiometry, species, annotation, rule and notes and
   each gene's annotation and reactions, fails closed, rejects a partly applied state); flag `--r171-atp-synthase-gpr`
   in `cli.py`; call in `main.py` after the aaRS step; test `platform/tests/test_r171_atp_synthase_gpr.py` on an
   in-memory model.
2. Full pytest before (`logs/pytest_before.txt`, run before any code change) and after.
3. Builds from `platform/` with `/opt/miniconda3/bin/python`: control (seed 0, scratch path), candidate (seed 0,
   `candidates/E5_…_aars_thr_ser_split_r171_atp_synthase_gpr_20261029.xml`), seed check (seed 1, scratch path).
4. Whole-model diff with the copied `model_diff.py`; knockout sweep with the copied `ko_compare_two.py`
   (Task A medium + PO1f, Gurobi 1e-9, reset per knockout, targeted check).
5. Independent read-only audit; REPORT.md; `model/STATE.md` entry; `candidates/README.md` row; Escher before/after
   figure of R171 and the respiratory chain. Stop before commit.

**Compute budget:** about 2,164 LPs for the sweep, at most 40 for the audit and 20 for the figure; 30 minutes. Stop
conditions: control not reproducing `d8d1bac7`, a seed difference, an unexpected diff, a failing sweep point; report,
never re-tune.

**Pytest interpreter:** this worktree has no venv; the suite runs with
`../workflow-automation-review-6765da/platform/.venv/bin/python` (Python 3.13.5, cobra 0.30.0, pytest 9.0.2), which
imports this worktree's `scripts` and `tools` (checked). Builds and sweeps use `/opt/miniconda3/bin/python`.
