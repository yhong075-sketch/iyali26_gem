"""Build the curation databases once with the skill's build_curation_dbs and record provenance.

Usage: python build_dbs.py SKILL_DIR DBDIR MODEL_GENES_CSV PROVENANCE_JSON
"""
import datetime
import hashlib
import json
import os
import sys


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(skill_dir, dbdir, model_genes_csv, prov_json):
    sys.path.insert(0, skill_dir)
    import kernel
    started = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    gm = kernel.build_curation_dbs(dbdir, model_gene_csv=model_genes_csv)
    files = {n: sha256(os.path.join(dbdir, n)) for n in ("w29.faa", "w29.gff", "sprot.fasta", "w29_gene_map.csv", "sprot_titles.tsv")}
    n_w29 = sum(1 for l in open(os.path.join(dbdir, "w29.faa")) if l.startswith(">"))
    n_sp = sum(1 for l in open(os.path.join(dbdir, "sprot.fasta")) if l.startswith(">"))
    prov = dict(
        started_utc=started,
        finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        assembly=kernel.W29_ASSEMBLY, kernel_sha256=sha256(os.path.join(skill_dir, "kernel.py")),
        dbdir=dbdir, files_sha256=files, w29_proteins=n_w29, sprot_entries=n_sp,
        gene_map_rows=len(gm), gene_map_yali1_locus=int(gm.locus_tag.str.startswith("YALI1_").sum()),
        model_genes_mapped=int(gm.in_model.sum()) if "in_model" in gm else None)
    json.dump(prov, open(prov_json, "w"), indent=2)
    print(json.dumps(prov, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:5])
