"""List reviewed UniProt entries annotated with each target reaction's Rhea IDs (exact reaction),
or its EC numbers when no Rhea ID is annotated, to choose characterised query proteins.

Usage: python query_candidates.py reaction_context.csv OUT.tsv
"""
import csv
import io
import sys
import time
import urllib.parse
import urllib.request

FIELDS = "accession,id,protein_name,gene_names,organism_name,lineage,ec,rhea,protein_existence,annotation_score,cc_subcellular_location,length"


def get(url):
    for a in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as f:
                return f.read().decode()
        except Exception:
            if a == 2:
                raise
            time.sleep(3)


def search(query, size=200):
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": query, "format": "tsv", "fields": FIELDS, "size": size})
    return list(csv.DictReader(io.StringIO(get(url)), delimiter="\t"))


def main(ctx_csv, out_tsv):
    rows = list(csv.DictReader(open(ctx_csv)))
    out = []
    for r in rows:
        rheas = [x for x in r["rhea"].split(";") if x]
        ecs = [x for x in r["ec"].split(";") if x and not x.endswith("-")]
        if rheas:
            q, basis = "(" + " OR ".join("rhea:" + x for x in rheas) + ") AND reviewed:true", "rhea"
        elif ecs:
            q, basis = "(" + " OR ".join("ec:" + x for x in ecs) + ") AND reviewed:true", "ec"
        else:
            out.append(dict(reaction=r["id"], basis="none", query="", n=0))
            continue
        hits = search(q)
        for h in hits:
            lin = h.get("Taxonomic lineage", "")
            out.append(dict(reaction=r["id"], basis=basis, query=q, n=len(hits),
                            fungal="Fungi" in lin, **{k: h[k] for k in h if k != "Taxonomic lineage"}))
        print(r["id"], basis, len(hits), "fungal", sum("Fungi" in h.get("Taxonomic lineage", "") for h in hits), flush=True)
    keys = []
    for o in out:
        keys += [k for k in o if k not in keys]
    with open(out_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
        w.writeheader()
        w.writerows(out)


if __name__ == "__main__":
    main(*sys.argv[1:3])
