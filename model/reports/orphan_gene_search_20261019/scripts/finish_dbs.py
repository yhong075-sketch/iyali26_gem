"""Finish build_curation_dbs after its downloads succeeded but makeblastdb failed on a path with a
space (BLAST+ splits -in/-db at spaces). DBDIR must be a space-free path (a symlink to the research
root db folder). The steps after the downloads are copied from the skill's kernel.py
build_curation_dbs (SHA-256 19fe8cfa...) without change; nothing is downloaded again.

Usage: python finish_dbs.py DBDIR MODEL_GENES_CSV PROVENANCE_JSON
"""
import datetime
import hashlib
import json
import os
import subprocess
import sys
import urllib.parse
import urllib.request

import pandas as pd


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(dbdir, model_gene_csv, prov_json):
    assert " " not in dbdir, "DBDIR must not contain spaces"
    # --- copied from kernel.build_curation_dbs, after the two downloads ---
    for name, fa in (("w29", "w29.faa"), ("sprot", "sprot.fasta")):
        subprocess.run(["makeblastdb", "-in", os.path.join(dbdir, fa), "-dbtype", "prot", "-out", os.path.join(dbdir, name), "-parse_seqids"],
                       check=True, capture_output=True)
    with open(os.path.join(dbdir, "sprot_titles.tsv"), "w") as out:
        for line in open(os.path.join(dbdir, "sprot.fasta")):
            if line.startswith(">"):
                head = line[1:].rstrip("\n"); acc = head.split("|")[1] if "|" in head else head.split()[0]
                out.write(acc + "\t" + head.split(" ", 1)[-1] + "\n")
    genes = {}; cds = {}
    for line in open(os.path.join(dbdir, "w29.gff")):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9:
            continue
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "gene" and "locus_tag" in a:
            genes[a["locus_tag"]] = a.get("old_locus_tag", "")
        if f[2] == "CDS" and "protein_id" in a and a["protein_id"] not in cds:
            cds[a["protein_id"]] = (a.get("locus_tag", ""), urllib.parse.unquote(a.get("product", "")))
    gm = pd.DataFrame([dict(protein_id=p, locus_tag=lt, old_locus_tag=genes.get(lt, ""), product=pr) for p, (lt, pr) in cds.items()])
    gm["model_gene_id"] = gm.locus_tag.str.replace("_", "", n=1, regex=False)
    mg = pd.read_csv(model_gene_csv).fillna("")
    gm["in_model"] = gm.model_gene_id.isin(set(mg.gene))
    gm = gm.merge(mg[["gene", "rxns"]].rename(columns={"gene": "model_gene_id", "rxns": "model_reactions"}), on="model_gene_id", how="left")
    gm.to_csv(os.path.join(dbdir, "w29_gene_map.csv"), index=False)
    # --- end of copied steps ---
    with urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/search?query=accession:P69905&size=1&format=tsv&fields=accession", timeout=60) as r:
        release = r.headers.get("X-UniProt-Release"), r.headers.get("X-UniProt-Release-Date")
    files = {n: sha256(os.path.join(dbdir, n)) for n in ("w29.faa", "w29.gff", "sprot.fasta", "w29_gene_map.csv", "sprot_titles.tsv")}
    mtimes = {n: datetime.datetime.fromtimestamp(os.path.getmtime(os.path.join(dbdir, n))).isoformat(timespec="seconds")
              for n in ("w29.faa", "w29.gff", "sprot.fasta")}
    model_genes = set(mg.gene)
    prov = dict(
        finished_local=datetime.datetime.now().isoformat(timespec="seconds"),
        assembly="GCF_001761485.1", download_mtimes_local=mtimes,
        uniprot_release_at_finish=release[0], uniprot_release_date=release[1],
        files_sha256=files,
        w29_proteins=sum(1 for l in open(os.path.join(dbdir, "w29.faa")) if l.startswith(">")),
        sprot_entries=sum(1 for l in open(os.path.join(dbdir, "sprot.fasta")) if l.startswith(">")),
        gene_map_rows=len(gm), gene_map_yali1_locus=int(gm.locus_tag.str.startswith("YALI1_").sum()),
        model_genes=len(model_genes), model_genes_mapped=len(model_genes & set(gm.model_gene_id)),
        model_genes_unmapped=sorted(model_genes - set(gm.model_gene_id)))
    json.dump(prov, open(prov_json, "w"), indent=2)
    print(json.dumps({k: v for k, v in prov.items() if k != "model_genes_unmapped"}, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:4])
