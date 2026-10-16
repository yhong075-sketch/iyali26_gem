"""Task D step 3: fetch protein sequences for the orthology check (network step; no LP).

Y. lipolytica: every GPR gene of the cut-set reactions (R533, R534, R831, R817) and of the prompt-named enzymes
(MTHFD R549/R550, GDH R329/R330, ALD R1378/R134/R135/R136/R2070, IDH R487/R1385) plus homoisocitrate DH R398.
The protein is fetched by the gene's RefSeq annotation in the model; genes without one are searched in NCBI
Protein by locus tag (YALI1_<chromosome><number>g). S. cerevisiae: reviewed UniProt entries (organism 559292) by
gene name. Raw responses and FASTA files go to $IYALI26_RESEARCH_ROOT/artifacts/tasks/<task>/sequences/; the
manifest with every SHA-256 is written here.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/complexI_transhydrogenation_20261012/fetch_sequences.py
Writes sequences_manifest_v3.json here. Never overwrites. (v2 outputs, kept: a locus tag matched both the RefSeq
and the CLIB89 GenBank record and no protein was taken; v3 takes the RefSeq XP_ record. v1 outputs, kept: the locus tag was built as
YALI1_1E... instead of YALI1_E..., and gene_exact matched synonyms for ALD3/ALD4/ALD5; both fixed here.)
"""
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone

import cobra

from common import ART, HERE, PINS, check_pins, sha, write_new

REACTIONS = ("R533", "R534", "R831", "R817", "R549", "R550", "R329", "R330", "R1378", "R134", "R135", "R136",
             "R2070", "R487", "R1385", "R398")
SC_GENES = ("MDH1", "MDH2", "MDH3", "CTP1", "YHM2", "MTD1", "ADE3", "MIS1", "GDH1", "GDH2", "GDH3", "ALD2", "ALD3",
            "ALD4", "ALD5", "ALD6", "IDP1", "IDP2", "IDP3", "IDH1", "IDH2", "LYS12")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def get(url):
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode()


def main():
    inputs = check_pins()
    out = ART / "sequences"
    out.mkdir(parents=True, exist_ok=True)
    model = cobra.io.read_sbml_model(str(PINS["model"][0]))
    genes = {}
    for rid in REACTIONS:
        for g in model.reactions.get_by_id(rid).genes:
            genes.setdefault(g.id, {"reactions": [], "refseq": g.annotation.get("refseq")})["reactions"].append(rid)
    manifest = {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs, "yali": {}, "sc": {}}
    yl_fasta = []
    for gid, info in sorted(genes.items()):
        accession, how = info["refseq"], "model refseq annotation"
        if not accession:
            tag = f"YALI1_{gid[5:]}"
            ids = json.loads(get(EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(
                {"db": "protein", "term": f"{tag}[Locus Tag]", "retmode": "json"})))["esearchresult"]["idlist"]
            time.sleep(0.4)
            summary = json.loads(get(EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(
                {"db": "protein", "id": ",".join(ids), "retmode": "json"})))["result"] if ids else {}
            time.sleep(0.4)
            hits = [summary[i]["accessionversion"] for i in ids]
            refseq = [a for a in hits if a.startswith("XP_")]
            accession = refseq[0] if len(refseq) == 1 else None
            how = f"NCBI Protein locus-tag search {tag}: {hits}; RefSeq XP_ record chosen"
        record = {"reactions": info["reactions"], "accession": accession, "source": how}
        if accession:
            fasta = get(EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(
                {"db": "protein", "id": accession, "rettype": "fasta", "retmode": "text"}))
            time.sleep(0.4)
            header, *body = fasta.strip().splitlines()
            record.update(header=header, length=len("".join(body)))
            yl_fasta.append(f">{gid} {header[1:]}\n" + "\n".join(body))
        manifest["yali"][gid] = record
        print(gid, record.get("accession"), record.get("length"), record["source"][:60], flush=True)
    sc_fasta = []
    for gene in SC_GENES:
        query = urllib.parse.urlencode({"query": f"gene_exact:{gene} AND organism_id:559292 AND reviewed:true",
                                        "format": "fasta"})
        fasta = get("https://rest.uniprot.org/uniprotkb/search?" + query)
        entries = [e for e in fasta.split(">") if e.strip() and f"GN={gene} " in e.splitlines()[0]]
        manifest["sc"][gene] = {"entries": [e.splitlines()[0] for e in entries]}
        if len(entries) == 1:
            header, *body = entries[0].strip().splitlines()
            sc_fasta.append(f">{gene} {header}\n" + "\n".join(body))
            manifest["sc"][gene]["length"] = len("".join(body))
        print(gene, manifest["sc"][gene], flush=True)
        time.sleep(0.4)
    for name, text in (("yali_task_d_v3.fasta", "\n".join(yl_fasta) + "\n"), ("sc_task_d_v3.fasta", "\n".join(sc_fasta) + "\n")):
        path = out / name
        with open(path, "x") as handle:
            handle.write(text)
        manifest[name] = {"path": str(path), "sha256": sha(path)}
    write_new(HERE / "sequences_manifest_v3.json", json.dumps(manifest, indent=1) + "\n")


if __name__ == "__main__":
    main()
