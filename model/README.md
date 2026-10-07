# iYali26 model

Genome-scale metabolic model of *Yarrowia lipolytica* (Wheeldon lab). This folder holds the
biology only: the model, the decisions behind it, its evidence and its expected behaviour. The
code that builds and checks it is in [../platform](../platform).

## Where things are

| Folder | Contents |
| --- | --- |
| [source/](source/) | `iyali26.xml`, the immutable build input; `iyli21.xml`, the predecessor model used for comparisons; and the reference-build provenance record |
| [curation/](curation/) | Every decision the build applies, as data, grouped by topic. [Curation log](curation/README.md) |
| [conditions/](conditions/) | Media used by the build and by tools |
| [candidates/](candidates/) | Built model files with provenance. [Index](candidates/README.md) |
| [expected/](expected/) | How predictions are scored ([benchmark.md](expected/benchmark.md)) and the memote experimental-data config |
| [reports/](reports/) | One folder per task: TASK.md, REPORT.md, audits and small tables |
| [STATE.md](STATE.md) | Dated decision log, newest first |
| [model.toml](model.toml) | Folder layout, read by the platform |

Large outputs (flux dumps, run folders, figures) are not in git. They live in the research
workspace under `$IYALI26_RESEARCH_ROOT/artifacts/`.

## Current state (checked 2026-10-07)

- **Focus:** resolving inactive reactions, such as CoQ9 synthesis.
- **Reference model:** not designated yet. The candidates are listed in
  [candidates/README.md](candidates/README.md).
- **Default build:** SHA-256 `b4ce0974…`, with 2,314 reactions, 1,877 metabolites and 1,073
  genes. It equals `candidates/model_metadata_trna_r1159_leak.xml` plus the 2026-10-05 R1889
  GPR.
- **Scoring:** 15% KO/WT is the primary cutoff, and model genes not on the experimental list
  count as non-essential ([benchmark.md](expected/benchmark.md)).

## How a change gets in

1. A scoped, user-authorized decision is written as a data file in `curation/`.
2. A guarded step in the platform applies it.
3. The model is rebuilt to a new file.
4. The whole model is compared with the previous one.
5. The change is audited independently.
6. A report goes in `reports/<topic>_<yyyymmdd>/` and an entry in `STATE.md`.

Model XML files are never edited by hand. See [../AGENTS.md](../AGENTS.md).
