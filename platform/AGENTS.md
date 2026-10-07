# platform/ — code rules

Applies to everything under `platform/`, in addition to the root AGENTS.md. Run commands from
this folder.

- **Model writes:** only the builder (`scripts/gem_annotate/`) writes model XML, and only
  through the deterministic writer.
  - Tools work on in-memory copies.
  - When a tool must write a model, it writes a new, separately named file and refuses the
    canonical `model/candidates/legacy/model.xml`.
- **Model locations:** find model files through `scripts.gem_annotate.model_layout.MODEL`,
  which reads `model/model.toml`. Never join paths onto the repository root by hand.
  - Curation files: `MODEL.curation_file(name)`.
  - Candidates: `MODEL.candidate_file(name)`.
  - Task folders: `MODEL.reports / name`.
- **Recorded paths:** paths recorded by earlier runs (in `.build.json`, task configs and
  handoffs) use the former layout. Resolve them with `config.resolve_recorded_path`, and never
  rewrite the recorded files.
- **Adding a science change:**
  - Put the decision in a curated data file under `model/curation/<topic>/`, never in a
    hard-coded constant.
  - Add one idempotent, guarded `apply_*` function that fails closed when its input rows are
    stale or do not resolve to exactly the expected IDs.
  - Call it from `main.py` at the right stage.
  - Add a focused test that runs on an in-memory model.
- **Fail closed, never skip:** a missing curated input that changes the model is an error,
  not a warning.
- **Chemistry:**
  - Set formula and charge together.
  - Every edit must keep mass and charge balance for every touched reaction, or roll back.
  - Never add H+ or H2O unless hydrogen, oxygen and charge all close.
- **Network steps:** cache results with their source version. A build that depends on an
  uncached live query is not reproducible; say so in the report.
- **Text that ends up in the model:** reaction notes are part of the model's SHA-256. Change
  them only as a deliberate, reported model change.
- **Tests:**
  - Write temporary outputs to `SCRATCH_DIR` (git-ignored `.scratch/`) or a system temp
    folder, never into `model/`.
  - Before handing back, run `python -m pytest` and compare with the previous run test by
    test.
