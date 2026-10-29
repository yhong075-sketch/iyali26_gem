"""Pin the W29 sequences of the two prolyl-tRNA synthetase genes for the localisation predictors (read-only inputs).

Usage (from the repo root): python3 model/reports/r673_prolyl_localization_20261029/scripts/pin_sequences.py RESEARCH_ROOT
Takes each sequence from the W29 RefSeq proteome (GCF_001761485.1, cached by the 2026-10-19 orphan search) and stops
unless it equals the W29 UniProt proteome UP000182444 entry named in the model's gene annotation. Also records the
CLIB122 UniProt sequence used by the OR-rule review, for comparison only (not submitted).
"""
import hashlib
import json
import sys
from pathlib import Path

root = Path(sys.argv[1]) / "artifacts/tasks"
OUT = Path("model/reports/r673_prolyl_localization_20261029")
REFSEQ = root / "orphan_gene_search_20261019/db/w29.faa"
UNIPROT = root / "argininosuccinate_identity_20260911/UP000182444.fasta"
REVIEW_TSV = Path("model/reports/or_rule_review_r171_aars_20261029/evidence/uniprot_aars_atp_synthase_20261029.tsv")
GENES = {"YALI1E05929g": ("XP_503566.2", "A0A1H6Q518", "Q6C6Z6"), "YALI1D19110g": ("XP_502869.1", "A0A1D8NEQ1", "Q6C8Z3")}


def read(path):
    seqs, key = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            key = line[1:]
            seqs[key] = ""
        else:
            seqs[key] += line.strip()
    return seqs


sha = lambda s: hashlib.sha256(s.encode()).hexdigest()
refseq, uniprot = read(REFSEQ), read(UNIPROT)
review = {row.split("\t")[0]: row.split("\t")[-1] for row in REVIEW_TSV.read_text().splitlines()[1:]}
record, fasta = {"sources": {"refseq_w29": {"path": str(REFSEQ), "sha256": sha(REFSEQ.read_text())},
                             "uniprot_w29": {"path": str(UNIPROT), "sha256": sha(UNIPROT.read_text())},
                             "review_tsv": {"path": str(REVIEW_TSV), "sha256": sha(REVIEW_TSV.read_text())}},
                 "genes": {}}, []
for gene, (xp, w29_acc, clib_acc) in GENES.items():
    header, seq = next((k, v) for k, v in refseq.items() if k.split()[0] == xp)
    if f"YALI1_{gene[5:]}" not in header or uniprot[next(k for k in uniprot if k.split()[0] == w29_acc)] != seq:
        sys.exit(f"STOP: {gene} RefSeq and UniProt W29 sequences differ or the header does not name the gene")
    clib = review[clib_acc].strip()
    offset = seq.find(clib)
    record["genes"][gene] = {"refseq": xp, "uniprot_w29": w29_acc, "length": len(seq), "sha256": sha(seq),
                             "refseq_header": header, "n_terminus_40": seq[:40],
                             "clib122_uniprot": clib_acc, "clib122_length": len(clib),
                             "clib122_equal": clib == seq,
                             "clib122_is_suffix_from": offset if offset >= 0 else None}
    fasta.append(f">{gene}_W29_{xp}\n{seq}")
(OUT / "target.fasta").write_text("\n".join(fasta) + "\n")
record["target_fasta_sha256"] = sha((OUT / "target.fasta").read_text())
(OUT / "sequence_identity.json").write_text(json.dumps(record, indent=1) + "\n")
print(json.dumps({g: {k: v[k] for k in ("length", "clib122_length", "clib122_equal", "clib122_is_suffix_from")}
                  for g, v in record["genes"].items()}))
