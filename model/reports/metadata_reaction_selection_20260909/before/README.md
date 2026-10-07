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
builder. The complete reference chain starts from `data/iyali26.xml`, retaining
its tRNA biomass representation, curated chemical convention, reaction directions
and GPRs. The restored unmodified chain reproduced the frozen reference SHA256
`bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee` exactly.

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
[model_reference_metadata_r1172_sphpl.xml](model_reference_metadata_r1172_sphpl.xml),
with [build provenance](model_reference_metadata_r1172_sphpl.build.json).
It has 2315 reactions, 1877 metabolites and 1074 genes. Compared with the frozen
reference, mathematical changes are confined to the two added reactions; all
original reaction stoichiometry, bounds, GPR logic and species properties remain.
The previous 2295-reaction metadata build is retained as a historical artifact.

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
`python -m unittest tests.test_coq9_curation tests.test_quinone_pipeline tests.test_er_vlcfa_stereochemistry tests.test_reference_reactions`.


---

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
