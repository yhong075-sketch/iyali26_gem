# CoQ9 reaction/gene curation package — 2026-09-09

Start with REPORT_zh.md. This package contains two **isolated, non-production** model candidates and reproducible static evidence.

1. model_metadata_only_candidate.xml preserves the mathematical model; scientific gaps remain documented.
2. model_R305_Qcycle_candidate.xml additionally changes ONLY the two R305 proton coefficients. It is a full-Q-cycle mechanism candidate using the existing cytosolic H as a P-side proxy. It is not a whole-model validation or a completed GPR repair.

No production repository, branch, or remote files were changed. No optimizer was run.

Reproduce with Python 3.10+ and lxml:

```sh
python -m pip install lxml
python audit_and_build.py --model /path/to/registered/model.xml --out ./reproduced
```

The script refuses inputs other than SHA256 bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee.

Do not treat old n_reactions values in evidence/iyli21_genes_vs_S2.csv as current reaction counts; current associations were freshly extracted into gene_identity_and_roles.tsv.

The remaining CI gene-association and COQ6 reaction-network candidates are specified in scientific_change_spec.json, NOT silently applied. Source accession bridges and experimentally demonstrated function are explicitly distinguished.
