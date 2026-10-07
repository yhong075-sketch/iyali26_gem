[![memote tested](https://img.shields.io/badge/memote-tested-blue.svg?style=plastic)](https://mcnaughtonadm.github.io/iyali26)

# iYali26

Genome-scale metabolic model of *Yarrowia lipolytica* (Wheeldon lab). A memote history report
is published at https://wheeldon-lab.github.io/iyali26_gem.

The repository has two parts, and the boundary between them is
[model/model.toml](model/model.toml):

| Folder | For | Start here |
| --- | --- | --- |
| [model/](model/) | The biology: the model, curation decisions, evidence, expected behaviour, reports | [model/README.md](model/README.md) |
| [platform/](platform/) | The tooling that builds, edits and checks the model | [platform/README.md](platform/README.md) |

Agent and contributor rules are in [AGENTS.md](AGENTS.md).

## Quick start

```bash
export IYALI26_RESEARCH_ROOT=/path/to/iyali26_gem_research
cd platform
python -m scripts.gem_annotate --research-root "$IYALI26_RESEARCH_ROOT" \
  --offline --no-solve --coq9-curation metadata --output-model ../model/candidates/NEW_NAME.xml
python -m pytest
```

`run.sh` currently stops at its first step. Its default build output is the canonical
`model/candidates/legacy/model.xml`, which the builder never overwrites. That was already the
case before the restructure.

---

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
