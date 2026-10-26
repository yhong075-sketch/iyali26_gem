# Heme a and iron(2+) in biomass_C, Yeast8 values (2026-10-22)

## Authorization

User, 2026-10-22, after reading plan A in this session:

> Approved: implement plan A as written (heme a 1e-6 + Fe2+ 3.04e-5 in biomass_C via an opt-in
> flag, baseline cf0825b8, iron cap unchanged, KO comparison included). No commit or push.

The trigger was a Claude Science test ("Gap 5 · Iron uptake cap 1.61e-5"). It showed that no
biomass term contains iron, so iron uptake and heme synthesis carry zero flux. In its in-memory
test, heme b + heme a at 5e-5 mmol/gDW each made the iron cap bind and growth fall from 0.31243
to 0.161 h⁻¹.

## Starting state

- **Branch:** `claude/science-test-gap-fix-5eaeae` at `abe7e40`, clean.
- **Baseline model:** `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_coproducts_20261021.xml`,
  SHA-256 `cf0825b86c3533d8a12f697a994eba3d740671a68a367143ee32cbfdc8686b7c`. This is the newest
  candidate (`164f235f` plus `--gpi-chain-links --gpi-coproducts`). It is a candidate, not a
  reference.
- **Before-state of the targets and their neighbours:** [before.json](before.json), written by
  [scripts/before_state.py](scripts/before_state.py).
  - `biomass_C` contains no iron-containing species.
  - The closed legacy biomass `R1710` ("yeast 8 biomass pseudoreaction", bounds 0/0) contains heme
    a [C_cy] at 1e-6. R1709 (heme a transport, C_mi ⇌ C_cy, gene-less, −1000/1000) feeds it.

## Gap reproduced before the change

Verified this session in memory: `cf0825b8`, Task A medium + PO1f, GLPK, pFBA.

| Case | Growth | R1189 (iron exchange) |
| --- | --- | --- |
| As supplied | 0.312427 | 0 |
| heme b [C_mi] + heme a [C_mi], 5e-5 each (the figure's test) | 0.161000 | −1.61e-5 (at the cap) |
| Same, iron cap removed | 0.312390 | −3.12e-5 |
| heme a [C_cy] 1e-6 + iron(2+) [C_cy] 3.04e-5 (this change) | 0.312426 | −9.81e-6 (61 % of the cap) |

## Target

| Item | Before | After |
| --- | --- | --- |
| `biomass_C`, heme a `m1635[C_cy]` | absent | −1e-6 mmol/gDW |
| `biomass_C`, iron `m948[C_cy]` | absent | −3.04e-5 mmol/gDW |
| `biomass_C` notes | five existing notes | the same notes plus `heme_iron_biomass_*` notes |

- No other coefficient, bound, GPR, species or reaction changes.
- The medium is unchanged. The R1189 cap (1.61e-5, the 2026-10-12 concentration-ratio rule) stays.
- The step is `apply_heme_iron_biomass()` in `platform/scripts/gem_annotate/heme_iron_biomass.py`,
  under the opt-in flag `--heme-iron-biomass`. It requires the CoQ9 biomass dilution stage
  (`biomass_C` as built by `--energy-candidate E5 --coq9-biomass-dilution`). Data:
  `model/curation/biomass/heme_iron_biomass.json`.

## Evidence

- **Values (verified 2026-10-22):** yeast-GEM v9.1.1, main commit
  [`2d594ae1`](https://github.com/SysBioChalmers/yeast-GEM/blob/2d594ae1c4a2d550ccef120d96a58c7bbf586255/model/yeast-GEM.yml).
  - `r_4598` cofactor pseudoreaction: heme a [c] (`s_3714`) −1e-6.
  - `r_4599` ion pseudoreaction: iron(2+) [c] (`s_0924`) −3.04e-5.
  - No heme b in any biomass pseudoreaction.
  - Source table
    [`biomassComposition_Cofactor_Ion.tsv`](https://github.com/SysBioChalmers/yeast-GEM/blob/2d594ae1c4a2d550ccef120d96a58c7bbf586255/data/physiology/biomassComposition_Cofactor_Ion.tsv)
    (added by [PR #174](https://github.com/SysBioChalmers/yeast-GEM/pull/174)): heme a cites
    "Forster2003"; iron(2+) cites "DOI: 10.22161/ijeab/2.2.2".
- **Primary measurements (recalled, not opened):** the two cited primary sources were not read. The
  values are *S. cerevisiae* consensus-model values, not *Y. lipolytica* measurements.
- **Internal precedent (verified):** the source model's own biomass `R1710` carries heme a [C_cy]
  at 1e-6, the same species and coefficient.
- **The figure's 5e-5:** no source was found (see the planning turn); not used.

## Expected effect

- WT growth on the Task A medium + PO1f: 0.312426 (in-memory test above), a fall of about 1e-6.
- Iron uptake about 9.8e-6, below the 1.61e-5 cap.
- Heme b → heme o → heme a (R302, R384, R383), R674, R1709, R2193 and R863 now carry flux. Their
  genes, and the upstream heme genes, may become essential in single knockouts (inferred). That
  is a side effect to report, not the reason for the change.

## Acceptance criteria

1. **pytest per test** equals the run before the change, plus the new tests passing.
2. **Builds** (the `cf0825b8` argv, from `platform/`):
   - without the new flag: reproduces `cf0825b8` bytes;
   - with `--heme-iron-biomass`: the same bytes under `PYTHONHASHSEED=0` and `1`.
3. **Whole-model diff** against `cf0825b8`: only the two `biomass_C` coefficients and the new
   `biomass_C` notes.
4. **Solves** (Task A medium + PO1f, GLPK):
   - WT growth within 1e-4 of 0.312427, and optimal;
   - R1189 between −1.61e-5 and 0; heme a made through R383 and exported by R1709;
   - closing R1189, R2193, R302, R383 or R1709 in turn gives zero growth or an infeasible solve.
5. **Single-gene knockouts** compared with `cf0825b8` on the same condition, at 1/5/10/15 %:
   - every class change is listed with its gene;
   - genes whose growth changes by more than 1e-6 are listed;
   - any non-optimal status that differs is solved again in a fresh context.
   - Both models are also scored against the experimental list per
     [benchmark.md](../../expected/benchmark.md) (mtDNA genes out of scope, user-defined negatives,
     development reference).

## Compute budget and stop condition

- **Budget:**
  - 3 builds;
  - about 20 diagnostic LPs;
  - 2 × (1,074 + 1) knockout solves;
  - at most 10 fresh re-solves.
- **Stop and report on:**
  - a guard rejection;
  - a control build that differs from `cf0825b8`;
  - a diff outside the targets;
  - seed-dependent bytes;
  - a WT change of more than 1e-4;
  - a test regression.

## Environment

- **Interpreter:** the main checkout's `.venv/bin/python` (→ `/opt/miniconda3/bin/python3`;
  Python 3.13.5, cobra 0.30.0, optlang 1.8.3, pytest 9.0.2).
- **Research root:** `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`.
- **Condition:**
  - medium [po1f_csm_leu_w29uptake_20261012.csv](../../conditions/media/po1f_csm_leu_w29uptake_20261012.csv)
    (`6eca70fd…`);
  - strain profile `state/strain_profiles/po1f_sd_leu.json` (`35307853…`).

## Findings recorded, not changed

- **R2193's only gene:** YALI1F07747g (the only import route for iron).
  - The model annotates it as ncbigene 2908631, which is [yli:2908631](https://www.kegg.jp/entry/yli:2908631).
    KEGG names that entry "Phosphorus acquisition-controlling protein" (bHLH). Verified 2026-10-22.
  - So either the GPR or the model's ID mapping is wrong. Which one is unchecked (inferred).
  - With iron required, this gene's knockout becomes lethal in the model.
- **Iron chemistry:** iron species m1273/m948/m523 have no formula and charge 0.
  - Imbalances verified with cobra 2026-10-22: R302 Fe +1, H −2, charge +2; R727 Fe +1, H −7,
    charge −5; R384 charge −2; R383 charge +2.
- **Siroheme m949 (R727):** no consumer; a separate dead end.
- **Iron cap rule:** at WT growth it allows at most about 5.15e-5 mmol Fe/gDW. A future heme b or
  Fe-S cluster biomass term above that total would make the cap bind.
- **Heme b:** not in Yeast8's biomass and not added here; no sourced value was found.

## Corrections after the independent audit (2026-10-22)

The text above is kept as written. These corrections apply. See
[REPORT.md §8](REPORT.md#8-independent-audit).

- **Wording:**
  - "a fall of about 1e-6": the measured fall is 4.5e-7.
  - "Heme b → heme o → heme a (R302, R384, R383)": R302 makes heme b (iron + protoporphyrin IX);
    R384 and R383 make heme o and heme a.
  - "the newest candidate": the baseline is `cf0825b8` because the user named it in the approval,
    not because it is newest.
- **Guard scope (extended):** the step now also checks that `biomass_C` is the only objective, and
  pins the bounds and stoichiometry of R1189, R863, R302, R384 and R383 besides R1709 and R2193.
  The heme steps upstream of R302 (R61 to R674) are still not pinned.
- **Compartment convention:** heme a is taken from `m1635[C_cy]`, following yeast-GEM (`s_3714`,
  [c]) and the source model's R1710. This makes the gene-less export R1709 required. Biologically,
  heme a stays in mitochondrial cytochrome c oxidase (established, uncited), so `m609[C_mi]` would
  also have been defensible. No gene result depends on this choice, because R1709 has no gene.
- **Budget (not fully kept):**
  - The first knockout sweep had no time limit per LP. It hung for 34 min in one GLPK solve and was
    stopped; its LP count was not recorded (up to about 2,100, inferred). It went beyond the
    authorized 2 × (1,074 + 1) sweep solves.
  - The rerun, with a 60 s limit per LP, was started without a separate check-in with the user. The
    user had asked why the run was slow, and was told before the rerun started.
  - Diagnostic LPs: 17, within budget. Fresh re-solves: 5. The planning turn also used 10 LPs
    before authorization.
