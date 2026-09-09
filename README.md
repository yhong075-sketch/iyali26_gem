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

## CoQ9 curation in the builder

`python -m scripts.gem_annotate` (also `python scripts/update_model.py`) accepts
`--coq9-curation off|metadata|qcycle`. The default is `metadata`: guarded name/EC
and PROTEIN_CLASS corrections plus Boolean-equivalent R2062 AND deduplication.
`qcycle` explicitly adds the R305 2/4 proton candidate, retaining the existing
cytosolic proton as a P-side proxy. `off` bypasses only this new curation.

Rules and local identity preconditions live in `data/coq9_curation.json`; gene
roles/evidence live in `data/coq9_gene_evidence.tsv`. Curation runs once after
annotation cleanup and GPR assembly, immediately before export. Each output has
an adjacent `.build.json` (input/code identities, options, per-item status and
conflicts) and `.coq9_genes.tsv` (evidence plus actual output associations).
Conflicts preserve current content; inspect `requested_build_complete` and the
per-item records, since successful export does not imply every curation passed.

All three modes include the existing Q9 chain in `scripts/gem_annotate/quinone.py`:
the CoQ6-to-Q9 conversion runs before FVA, and reviewed quinone step GPRs and
duplicate cleanup run after generic annotation. The chain was migrated from the
saved 2026-09-07 implementation; provenance is in
`data/quinone_pipeline_provenance.json`. The new correction runs once at the end.
R305's known incorrect source name is excluded from both automatic reaction
annotation passes in metadata/qcycle, so it cannot seed a false MNXR/EC identity.

Use `--output-model PATH` to keep validation outputs separate. `--offline`
disables remote gene annotation; `--no-solve` skips diagnostic FVA/precursor
solves while retaining construction. `--mnx-dir` and `--cache-dir` select local
data/cache locations. Building from the raw `data/iyli21.xml` requires the local
MetaNetX tables for the existing Q9 chemistry prerequisites. Missing prerequisites
are reported as conflicts, without substituting a previously built XML.

Focused checks: `python -m unittest tests.test_coq9_curation tests.test_quinone_pipeline tests.test_er_vlcfa_stereochemistry`.

Latest default build (2026-09-09):
[metadata.xml](artifacts/coq9_pipeline_integration_20260909/accepted/metadata.xml),
with [build provenance](artifacts/coq9_pipeline_integration_20260909/accepted/metadata.build.json)
and [validation results](artifacts/coq9_pipeline_integration_20260909/accepted/validation.json).
This is a pipeline candidate with 2295 reactions, not a replacement release for
the 2313-reaction frozen reference. All 33 focused tests passed; metadata/off
mathematical equivalence and the two-coefficient qcycle difference were verified.
The [closed ATP check](artifacts/coq9_pipeline_integration_20260909/accepted/static_pair.json)
failed in both modes (ATP dissipation 1000; R305 flux zero). Qcycle remains opt-in.
The saved build commands record the local MetaNetX snapshot used; those external
tables and the off/qcycle comparison XMLs are not included in this publication.



---

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
