[![memote tested](https://img.shields.io/badge/memote-tested-blue.svg?style=plastic)](https://mcnaughtonadm.github.io/iyali26)

# iyali26

A history report will be publicly visible at https://wheeldon-lab.github.io/iyali26_gem.


## Usage

All `memote` commands have extensive help descriptions.

1. For simple command line testing, check out `memote run -h`.
2. To generate a pretty report, check out `memote report snapshot -h`.


## Data
Download MetaNetX files and place in `data/metanetx/`:
- https://www.metanetx.org/ftp/latest/chem_prop.tsv
- https://www.metanetx.org/ftp/latest/chem_xref.tsv
- https://www.metanetx.org/ftp/latest/reac_xref.tsv

## Reference pipeline and CoQ9 curation

`python -m scripts.gem_annotate` and `python scripts/update_model.py` use the same
builder. The complete reference chain starts from `data/iyali26.xml`. Before the
final selections, the restored unmodified chain reproduced the frozen reference SHA256
`bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee` exactly.

The final authoritative selection in `data/metadata_reaction_selection.json`
uses the user-selected fields from the earlier 2295-reaction metadata output:
202 stoichiometries, 7 bounds and 9 GPRs across 206 reactions. The experimental
tRNA-coupled `biomass_C` is protected and checked separately. This stage runs after
reference chemistry/GPR/tRNA assembly and before CoQ9 curation and export, in
all three CoQ9 modes. Old/target local states are accepted; conflicting reactions
are preserved and reported. Compiled metadata XML is provenance, never a build
input. Species formula, charge and compartment conventions remain those of the
reference chain. Superseded reaction evidence is retained as historical notes;
version selection is not new experimental confirmation of chemistry or GPRs.

`biomass_C` consumes the 20 private protein residues produced through charged
tRNA incorporation, returning each uncharged tRNA carrier. Free amino acids do
not bypass those incorporation reactions. The metadata biomass equation is
explicitly excluded because it violates this experimental requirement. The seven
selected reactions, including R_NTP1 and R_PGAM1_PhosHydro, are reversible.

The two user-requested reactions are retained at their original construction
steps: R1172 is kept from the raw input; SPHPL is admitted by gap filling even
when R730 is present. Other duplicate filters and direction curation remain
active. Their former exclusion reasons remain in output notes and
`data/reference_build/retained_reactions.json`. R1172 has no confirmed carrier
GPR; SPHPL retains its supplied compartment assignment and reversible bounds.
In the current stored species convention, SPHPL has residual H = -1 and charge
= -1. It is a retained candidate, not a completed chemistry correction.

`--coq9-curation off|metadata|qcycle` defaults to `metadata`: guarded name/EC and
PROTEIN_CLASS corrections, historical provenance and Boolean-equivalent R2062
AND deduplication. `off` bypasses only these CoQ9 corrections. `qcycle` explicitly
adds the R305 2/4 proton candidate, retaining cytosolic H as the P-side proxy.
The CoQ9 stage follows final annotation, chemistry and GPR/tRNA assembly. No
unresolved CI or COQ6 architecture is activated.

Latest default output:
[model_metadata_trna.xml](model_metadata_trna.xml),
with [build provenance](model_metadata_trna.build.json).
It has 2315 reactions, 1877 metabolites and 1074 genes. Relative to the preceding
`model_reference_metadata_r1172_sphpl.xml`, changes are exactly the selected
202 stoichiometries, 7 bounds and 9 GPRs. Both retained reactions and all species
properties remain unchanged. Both preceding model outputs are preserved.

The preceding `model_metadata_reaction_selection.xml` and its corresponding
three-mode outputs are superseded: they incorrectly bypassed the required tRNA
biomass representation. They remain as historical evidence, not current models.
The corrected pipeline was rebuilt once in metadata mode. Seventeen focused
tests passed; all 20 incorporation reactions and biomass match the prior coupled
reference. A positive WT solution carries all 20 required incorporation fluxes,
and each private residue balance enforces `v_incorporation = a_i * v_biomass`.
Local mode checks retain CoQ9 idempotency and only two R305 coefficient changes
in explicit qcycle mode. See `artifacts/trna_biomass_restore_20260910/validation.json`.

The WT objective is 1.7793777729. The bounded closed-ATP check still reaches
1000 with R305 flux zero under preserved internal bounds: that separate energy
issue remains unresolved. This repair does not establish native gene essentiality
or complete chemical validity. Results: `artifacts/trna_biomass_restore_20260910/static_check.json`.

The saved PO1f / SD-Leu screen completed 1074 single-gene knockouts with optimal,
finite results. At the primary strict KO/WT < 10% threshold, 101 model genes are
predicted essential. The user-provided essential-gene workbook has 1612 positive
labels, of which 322 match screened model IDs: TP = 73, FN = 249, recall = 22.67%.
The 752 unlabelled model genes are not experimental negatives, so FP/TN and
accuracy are unavailable. At 1%, the model predicts 76 essential genes, including
55 reference positives. This is a development-reference comparison, not independent
validation. See the [screen summary](artifacts/screen_test_metadata_trna_20260910/screen_summary.json)
and [confusion matrix](artifacts/screen_test_metadata_trna_20260910/confusion_matrix_10pct.pdf).

Publication check (2026-09-10): an independent copy of the staged source/data
passed 41 focused tests and rebuilt this exact model byte-for-byte with
`--offline --no-solve`. See the [verification record](artifacts/model_metadata_trna_push_20260910/verification.json).

Versioned reference curation lives in `data/reference_build/`; CoQ9 rules and
gene evidence remain in `data/coq9_curation.json` and `data/coq9_gene_evidence.tsv`.
`data/reference_build_provenance.json` records the restored source and input
identities. The separate neutral ER helper remains bound to its original
`data/iyli21.xml` contract; it is not imposed on this reference build's chemistry.

### Optional ATP energy candidates (2026-09-24)

Default builds remain E0. `--energy-candidate E1|E2|E3|E4|E5` explicitly applies
`data/energy_candidate_repairs.json` after the reference selections and requires
`--offline --no-solve` plus a new output file. E1 restricts three phosphohydrolase
directions; E2 corrects their stored proton conventions; E3 combines them.
E4 adds ADP phosphohydrolase and IP7 hydrolysis candidates; E5 adds peroxidatic
ethanol oxidation and oxaloacetate decarboxylation candidates. These do not
validate native GPRs or compartments or constitute formal model acceptance.

For the pinned, completed local reference, the independent final-stage entry is:

```bash
.venv/bin/python -B scripts/build_energy_candidates.py \
  --variants E0 E1 E2 E3 E4 E5 --output-dir artifacts/my_energy_candidates
```

This entry reuses the saved reference, not an inaccessible historical rebuild.
It refuses overwrites, checks exact reaction/species signatures before mutation,
persists candidate precedence over later metadata chemistry/bounds, and verifies
the exported/reloaded model and solver definitions. E0 preserves loaded model
semantics; SBML serialization is not promised byte-identical. Original XML stays
unchanged. Runtime stage traces, energy/growth results and unresolved evidence
are in `artifacts/atp_candidate_repair_20260924/` and `PROJECT_STATE.md`.

Candidate gene functions are curated in
`data/reference_build/curation/gene_function_annotations.json` and applied after
UniProt ID mapping. YALI1F38820g is labelled a **Vph1-like V0 a subunit candidate,
requiring experimental confirmation**, based on sequence homology and an existing
AlphaFold prediction. Its native gene name and vacuolar versus Golgi/endosomal
role remain unconfirmed. This step writes only the display name and evidence
notes; it does not assign a compartment, GPR or essentiality label.
The local annotation build is [model_metadata_trna_vph1like.xml](model_metadata_trna_vph1like.xml).
A complete offline/no-solve rebuild and two focused tests passed; comparison
with the published `model_metadata_trna.xml` found only this gene's name/notes
changed. See the [verification](artifacts/f38820_annotation_pipeline_20260911/validation.json)
and [source audit](artifacts/f38820_annotation_pipeline_20260911/AUDIT.md).
This local update has not been committed or pushed.

The default pipeline now applies the user-authorized provisional enzymatic GPR
`R1026: YALI1F28274g` after reaction selection. The gene's native name is
unverified; it is an Nce103-like beta carbonic anhydrase candidate supported by
automatic annotation and an existing AlphaFold prediction, with native activity
and localization unresolved. The original empty rule and evidence limits are
retained in `data/reference_build/curation/r1026_gpr_assignment.json` and SBML notes.
The new [R1026 model](model_metadata_trna_r1026_gpr_build2.xml) differs from the
Vph1-like annotation build only in R1026 GPR/notes. A fresh WT and target single
knockout under the saved SD-Leu/PO1f conditions both grew at approximately
1.871882307 h^-1: R1026 and R2202 close, but the nuclear route compensates.
See the [scoped implementation and rerun report](artifacts/r1026_gpr_pipeline_20260911/REPORT.md).
Earlier generated models and screens remain historical artifacts.

The optional `--vatpase-gpr-hypothesis` build requires all three of
YALI1D00581g, **YALI0E16192g** (the existing CLIB122 identifier) and
YALI1F38820g in both R794 and R795: `D AND F AND a2 AND (original GPR)`.
It retains the original internal OR, so this is a partial complex hypothesis.
The native gene symbols remain unverified; the respective V1 D, V1 F and
Vph1-like V0 a assignments are sequence/AlphaFold-supported candidates.
The W29 F-locus identity and a-subunit compartment/redundancy remain unresolved.
The flag is disabled by default and requires `--offline --no-solve`, metadata
mode and a new output path. The full tRNA-coupled reference is built before the
two GPR edits; R795 remains closed. This hypothesis does not establish
experimental essentiality. Its fixed data and limitations are stored in
`data/reference_build/curation/vatpase_gpr_hypothesis.json`.
The generated [hypothesis model](model_metadata_trna_vatpase_and_hypothesis.xml)
passed the scoped build/logic checks. A full screen confirmed that the three
single knockouts disable both target GPRs but retain approximately WT growth;
all four threshold classifications were unchanged. See the
[hypothesis report](artifacts/vatpase_gpr_hypothesis_20260911/REPORT.md).

Supply a local research directory containing `reference/metanetx`,
`reference/ncbi`, `reference/kegg`, `reference/locus_map` and `cache/data`:

```sh
.venv/bin/python -m scripts.gem_annotate --research-root PATH \
  --offline --no-solve --coq9-curation metadata --output-model build/updated.xml
```

Use a new output path to preserve earlier artifacts. `--offline` uses the saved
annotation cache and blocks network access; `--no-solve` blocks optimization
while retaining every construction step. `--mnx-dir` and `--cache-dir` may
override their input directories. Each output has `.build.json` and
`.coq9_genes.tsv` sidecars. Inspect `requested_build_complete` and per-item
statuses; a completed build does not establish biological validity.

Focused checks:
`python -m unittest tests.test_reaction_selection tests.test_coq9_curation tests.test_quinone_pipeline tests.test_er_vlcfa_stereochemistry tests.test_reference_reactions tests.test_gene_function_annotations tests.test_vatpase_gpr_hypothesis tests.test_r1026_gpr_assignment tests.test_r1025_gpr_assignment tests.test_r153_gpr_assignment`.

The default build also applies the user-directed R153 rule
`YALI1D17462g` from `data/reference_build/curation/r153_gpr_assignment.json`.
The former YALIUNK2 placeholder (no established name or molecular identity,
model-only synthase role) is removed only when it has no other reaction links.
YALI1D17462g has no verified native symbol; sequence and AlphaFold prediction
support an A1 aspartyl endopeptidase candidate, conflicting with the assigned
argininosuccinate synthase function. The curation and exported notes retain
that conflict and the lack of experimental confirmation. Deleting the assigned
gene closes the retained R153 after the duplicate merge below; this static dependence does not
establish growth essentiality. The earlier OR model is preserved as history.

The default build now merges R2176 into R153 after this assignment, using
`data/reference_build/curation/r153_r2176_merge.json`. It checks both complete
reaction records before removing R2176, retains one R153 with bounds
`[-1000, 1000]`, and preserves the removed record and its source identity.
This deliberately reduces the former summed capacity of `[-2000, 2000]`;
it is not a claim of global feasible-set equivalence or resolved gene function.
Group membership is transferred to R153; species and genes are retained.
The new output is [model_metadata_trna_r153_merged.xml](model_metadata_trna_r153_merged.xml).
The original single-GPR model is preserved. The focused regression is
`python -m unittest tests.test_r153_merge`.

The default build also applies the authorized R1931 direction curation after
metadata selection, from `data/reference_build/curation/r1931_direction.json`.
It verifies reaction identity and restricts bounds to `[0, 1000]`, retaining
GSA oxidation to glutamate. GPR, stoichiometry and compartment are preserved;
notes retain the enzyme evidence and unresolved W29 kinetics/localization.
The new output is [model_metadata_trna_r1931_forward.xml](model_metadata_trna_r1931_forward.xml).
The earlier merged model and its screen remain historical results. This build
does not establish revised essentiality predictions. Focused regression:
`python -m unittest tests.test_r1931_direction`.

The default build also curates R539 after metadata selection, using
`data/reference_build/curation/r539_gpr_assignment.json`: EC `2.3.1.39` only,
and catalytic GPR `YALI1E22262g` (native symbol unverified; putative ACP
S-malonyltransferase supported by sequence/automatic annotation and a reused
AlphaFold prediction). Native activity and mitochondrial localization remain
unconfirmed; ACP stays in the stoichiometry, and carrier-gene dependencies
are not fully represented by this catalytic-only rule. The new output is
[model_metadata_trna_r539_corrected.xml](model_metadata_trna_r539_corrected.xml).
The user-requested evidence label is now `Alphafold预测， 仍需试验验证`, retained
by the default pipeline and verified in
[model_metadata_trna_r539_alphafold_labeled.xml](model_metadata_trna_r539_alphafold_labeled.xml).
This later export differs only in that R539 label. The prior model and screen are preserved; no revised growth/essentiality
results are implied. See the [correction report](artifacts/r539_correction_20260916/REPORT.md).
Focused regression: `python -m unittest tests.test_r539_gpr_assignment`.


---

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

## Conditional R1159 proton leak direction (2026-09-23)

The default build applies `data/reference_build/curation/r1159_direction.json`
after final field selection. R1159 keeps `H+[cy] -> H+[go]` storage and changes
bounds from `[-1000,1000]` to `[-1000,0]`, allowing net Golgi-to-cytosol leak.
This user-authorized condition assumes an acidic Golgi and a membrane potential
that does not reverse the electrochemical driving force; it is not a W29
measurement or a claim of permanent irreversibility. The inherited capacity
1000 is not measured permeability. Evidence, limitations and old bounds are
preserved in the curation and output notes. No GPR or medium is changed.

Unexpected reaction identity, stoichiometry, proton chemistry/compartments,
GPR, bounds or conflicting notes stop this step before it mutates the model.
Use a new output path with the existing `--offline --no-solve` build command;
previous model files are not replaced. Run
`python -m unittest tests.test_r1159_direction` for the focused regression.
The previously saved zero-flux WT witness remains feasible under this bound;
this direction change alone does not establish essentiality.


### Optional E5 vacuole connection candidate (2026-09-24)

The independent [E5_vacuole_open.xml](artifacts/vacuole_open_supply_20260924/E5_vacuole_open.xml)
opens only R1363/R795 to 0.04 and R871/R876 to 0.01 mmol/gDW/h.
The configuration is disabled by default and preserves the eight E5 energy repairs.
The builder now requires explicit source and SHA; the validator requires explicit baseline/candidate identities, mode, configuration, output and budget. See the [actually executed commands](artifacts/dipeptide_chemistry_routes_20260924/COMMANDS.md). There is no fallback to the historical candidate path.

The 12 closed ATP/GTP/UTP/CTP maxima for baseline, bounded opening and wide-bound
stress copies were zero. Direct artificial vacuolar dipeptide supply supported
four simultaneous hydrolyses, including under the declared near-optimal growth
requirement. Artificial inputs and all temporary constraints are excluded from
the XML. These are model connectivity results, not native enzyme/localization
or endogenous-supply validation. See the [report and reproducible commands](artifacts/vacuole_open_supply_20260924/REPORT.md).


### Optional dipeptide chemical metadata candidate (2026-09-24)

[E5_vacuole_open_chemistry.xml](artifacts/dipeptide_chemistry_routes_20260924/E5_vacuole_open_chemistry.xml) adds neutral-parent metadata to 13 pools: 12 conditional L-alpha free-dipeptide mappings and one named achiral glycine identity. The optional patch in `data/dipeptide_chemistry_patch.json` preserves the full optimization definition and eight E5 energy repairs. It preflights all fields, is idempotent, rejects conflicts and persists protection through reviewed metadata/automatic-balance paths. It does not establish native supply, activity or localization.

Two explicitly bound fresh validations used 22 optimizer calls total; ATP/GTP/UTP/CTP maxima remained zero, no-source growth and G3 behavior were unchanged. All 13 internal adjacent reactions balance; four exchanges are intentional material boundaries. External and endogenous supply specifications remain disabled pending evidence and parameters. See the [report](artifacts/dipeptide_chemistry_routes_20260924/REPORT.md), [route specification](artifacts/dipeptide_chemistry_routes_20260924/supply_module_spec.md) and [independent audit](artifacts/dipeptide_chemistry_routes_20260924/AUDIT.md). The original XML files and historical runs are retained.

### Optional CoQ9 functional GPR candidate

`--coq9-functional-gpr` adds an explicitly provisional substrate-access dependency
only to R695: `YALI1E18269g and YALI1F34675g`. These are the COQ7 hydroxylase
and COQ9 lipid-presentation family candidates; native obligatory dependence is
unverified. The option is disabled by default and requires an isolated
`--offline --no-solve --coq9-curation metadata` build and a new output path.
It can be combined with the separate CoQ C5 option below, but not unrelated
experimental overlays. On its own, no chemistry, bounds,
CoQ demand, medium or other GPR changes are made. The versioned rule and
counterevidence are in `data/reference_build/curation/coq9_functional_gpr.json`.

The CoQ9-only option does not assign COQ6, COQ8 or unknown transport steps,
complete every CoQ GPR or establish growth essentiality. See the
[CoQ completion report](artifacts/coq_gpr_completion_20261001/REPORT.md).

### Optional yeast-reference CoQ C5 candidate

`--coq-c5-gpr` replaces R39 with a mitochondrial NADH-coupled C5 hydroxylation
and the provisional rule `YALI1A08781g and YALI1B03314g and YALI1B19490g`:
COQ6-family hydroxylase, YAH1-like ferredoxin and ARH1-like reductase candidates.
Yeast experiments support this serial electron-transfer hypothesis; native
coupling and NADH preference remain unverified. The current NADP/NADPH formula
and charge conflict is preserved, and NADPH compatibility is not excluded.

This option is default-off, requires `--offline --no-solve --coq9-curation metadata`
and a new output path, and may be combined with `--coq9-functional-gpr`.
The old donor-free R39 is replaced, not retained in parallel. COQ8's reaction
assignment and R19 chemistry remain unresolved. The final combined
[candidate and validation report](artifacts/coq_yeast_reference_20261001/REPORT.md)
includes a runnable check at
`artifacts/coq_yeast_reference_20261001/verify_c5_candidate.py`.
