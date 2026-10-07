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
203 stoichiometries, 7 bounds and 9 GPRs across 207 reactions. It runs after
reference chemistry/GPR/tRNA assembly and before CoQ9 curation and export, in
all three CoQ9 modes. Old/target local states are accepted; conflicting reactions
are preserved and reported. Compiled metadata XML is provenance, never a build
input. Species formula, charge and compartment conventions remain those of the
reference chain. Superseded reaction evidence is retained as historical notes;
version selection is not new experimental confirmation of chemistry or GPRs.

This includes direct amino-acid consumption in `biomass_C`; the 20 tRNA residue
incorporation reactions remain present but no longer supply biomass. The seven
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
[model_metadata_reaction_selection.xml](model_metadata_reaction_selection.xml),
with [build provenance](model_metadata_reaction_selection.build.json).
It has 2315 reactions, 1877 metabolites and 1074 genes. Relative to the preceding
`model_reference_metadata_r1172_sphpl.xml`, changes are exactly the selected
203 stoichiometries, 7 bounds and 9 GPRs. Both retained reactions and all species
properties remain unchanged. Both preceding model outputs are preserved.

The three-mode semantic checks passed. Under the unchanged stored medium, WT
objective values are 1.3875914755 for the preceding output and 1.7793777729 for
both selected metadata and qcycle. Closing all boundary inputs while preserving
the maintenance lower bound gives an ATP dissipation objective of 1000 in both
new modes (R305 flux zero): the energy check fails. This version selection does
not establish chemical balance or improved biological validity. Exact outcomes,
including missing formula limitations, are in
`artifacts/metadata_reaction_selection_20260909/validation.json` and
`static_pair.json` in the same directory.

Versioned reference curation lives in `data/reference_build/`; CoQ9 rules and
gene evidence remain in `data/coq9_curation.json` and `data/coq9_gene_evidence.tsv`.
`data/reference_build_provenance.json` records the restored source and input
identities. The separate neutral ER helper remains bound to its original
`data/iyli21.xml` contract; it is not imposed on this reference build's chemistry.

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
`python -m unittest tests.test_reaction_selection tests.test_coq9_curation tests.test_quinone_pipeline tests.test_er_vlcfa_stereochemistry tests.test_reference_reactions`.


---

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
