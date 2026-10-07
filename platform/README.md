# iYali26 platform

Tooling that builds, edits and checks the model. It finds the model's files through
[../model/model.toml](../model/model.toml), so the same code can build another model folder by
setting `IYALI26_MODEL_DIR`.

## Layout

| Path | Contents |
| --- | --- |
| `scripts/gem_annotate/` | The builder: the only code that writes model XML. Import path `scripts.gem_annotate`. |
| `tools/` | Stand-alone tools, e.g. candidate builders, energy and supply diagnostics, dFBA (`tools.<name>`) |
| `tools/lipid/` | Lipid-unlump tools (owned by a colleague; review only) |
| `tests/` | Tests (`python -m pytest`) |
| `hpc/` | Slurm job scripts. They pin a specific past checkout and keep its old paths. |
| `ci/` | GitHub Actions helper for memote |

## Setup

Python 3.11 or newer. Run every command from this folder (`platform/`) so that `scripts` and
`tools` import from this checkout.

If an older editable install of the package exists in your virtual environment, it may point at
a different checkout. Running from `platform/` puts this checkout first; reinstalling the
environment removes the ambiguity.

| Variable | Meaning |
| --- | --- |
| `IYALI26_RESEARCH_ROOT` | Research workspace: MetaNetX, KEGG, NCBI, locus maps, annotation cache, ledgers and large outputs. Required for builds. |
| `IYALI26_MODEL_DIR` | Model folder to use. Defaults to `../model`. |
| `IYALI26_MNX_DIR`, `IYALI26_CACHE_DIR` | Optional overrides for the MetaNetX tables and the annotation cache |

## Build

```bash
python -m scripts.gem_annotate --research-root "$IYALI26_RESEARCH_ROOT" \
  --offline --no-solve --coq9-curation metadata --output-model ../model/candidates/NEW_NAME.xml
```

- **Output path:** must be new. The canonical `candidates/legacy/model.xml` is never
  overwritten.
- **Sidecar files:** every build writes `.build.json` (input and code SHA-256) and
  `.coq9_genes.tsv`.
- **Flags:**
  - `--offline` uses the saved annotation cache only.
  - `--no-solve` blocks optimization.
- **Optional, off by default:** `--coq9-functional-gpr`, `--coq-c5-gpr`,
  `--coq-literature-revision`, `--vatpase-gpr-hypothesis`. Each is described in
  [../model/curation/README.md](../model/curation/README.md).
- **Energy candidates:** build them with `tools/build_energy_candidates.py`. The builder's
  `--energy-candidate E5` currently stops with "Candidate export/reload changed model
  definition".

## Tests

```bash
python -m pytest
```

The tests build into the git-ignored `../.scratch/` folder and never write into `model/`. Some
tests fail or skip for reasons outside this checkout:

- **Lipid strict-sn tests (`test_lp_sn12_candidate`):** the frozen source fingerprint they
  expect no longer matches.
- **CoA protonation tests:** they compare against `../iyali26_gem/model.xml` next to this
  repository and expect SHA-256 `bc2aac8f…`.
- **Saved vacuole run tests:** they skip until the run is present under
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/`.
- **Two pinned checks:** the dFBA Slurm-script pin and the missing `docs/lipid_moiety_ledger.md`.

## Paths recorded by earlier runs

Older `.build.json`, task configs and handoff files record paths in the former layout
(`data/...`, `artifacts/...`, root `model_*.xml`, `scripts/...`). Those files stay
byte-identical. Code resolves such paths with
`scripts.gem_annotate.config.resolve_recorded_path`.
