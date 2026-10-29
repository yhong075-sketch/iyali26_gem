"""Build or_rule_review_genes.csv: one row per gene of the five rules with identity (S2 locus map, UniProt
TSV in evidence/), consensus and per-screen calls, and, when variant_knockouts.json exists, the knockout
growth ratio and 15 % class under every rule variant. Read-only apart from the CSV it writes.

Usage (repository root, IYALI26_RESEARCH_ROOT set):
    python3 model/reports/or_rule_review_r171_aars_20261029/scripts/build_gene_table.py
"""
import csv
import json
import os
from pathlib import Path

TASK = Path(__file__).resolve().parents[1]
RESEARCH = Path(os.environ["IYALI26_RESEARCH_ROOT"])
REPO_LIST = RESEARCH / "state/essentiality/repository"

ROLE = {  # rule membership and the subunit / enzyme role, from REPORT.md §3–§4
    "YALI1M00260g": ("R171", "AND member, mtDNA (ATP6/8/9 set)"), "YALI1M00267g": ("R171", "AND member, mtDNA (ATP6/8/9 set)"),
    "YALI1M00318g": ("R171", "AND member, mtDNA (ATP6/8/9 set)"),
    "YALI1F04433g": ("R171", "AND member, ATP1-like alpha"), "YALI1B05364g": ("R171", "AND member, ATP2-like beta"),
    "YALI1F04110g": ("R171", "AND member, ATP3-like gamma"), "YALI1D27826g": ("R171", "AND member, ATP16-like delta"),
    "YALI1D15734g": ("R171", "AND member, ATP5-like OSCP"), "YALI1F26967g": ("R171", "AND member, ATP4-like b"),
    "YALI1B08885g": ("R171", "AND member, ATP7-like d"), "YALI1F07384g": ("R171", "AND member, ATP14-like h"),
    "YALI1D14639g": ("R171", "AND member, ATP17-like f"),
    "YALI0D17490g": ("R171", "AND member, ATP18 subunit j (YALI0 id; = YALI1_D21476g in the S2 map, not screened)"),
    "YALI1E38131g": ("R171", "AND member, TIM11/ATP21 subunit e (dimer-specific)"),
    "YALI1B28025g": ("R171", "single-gene OR branch, ATP20 subunit g (dimer-specific)"),
    "YALI1B15733g": ("R171", "single-gene OR branch, ATP19 subunit k (dimer-specific)"),
    "YALI1F14882g": ("R759;R760", "OR member, THS1-like cytosolic ThrRS"), "YALI1A09121g": ("R759;R760", "OR member, MST1-like mitochondrial ThrRS"),
    "YALI1F03908g": ("R723", "OR member, SES1-like cytosolic SerRS"), "YALI1B07957g": ("R723", "OR member, DIA4-like mitochondrial SerRS"),
    "YALI1E05929g": ("R673", "OR member, YHR020W-like cytosolic ProRS"), "YALI1D19110g": ("R673", "OR member, AIM10/YER087W-like ProRS (mitochondrial by homology)"),
}

s2 = {r["yali1"].replace("_", ""): r["yali0"] for r in csv.DictReader(open(RESEARCH / "reference/locus_map/s2_metabolic_genes.csv"))}
s2["YALI0D17490g"] = "YALI0D17490g"
uni = {}
for r in csv.DictReader(open(TASK / "evidence/uniprot_aars_atp_synthase_20261029.tsv"), delimiter="\t"):
    key = (r["Gene Names (ORF)"] or r["Gene Names (ordered locus)"]).replace("_", "")
    uni[key] = r
positives = {r["gene_id"]: r["function"] for r in csv.DictReader(open(REPO_LIST / "consensus_essential_genes.csv"))}
fit = {}
for r in csv.DictReader(open(REPO_LIST / "cas9_cas12a_fitness.csv")):
    fit[(r["gene_id"], r["assay"])] = r
ko = None
if (TASK / "variant_knockouts.json").exists():
    ko = json.load(open(TASK / "variant_knockouts.json"))
variants = list(ko["variants"]) if ko else []
fields = ["gene_id", "rules", "role", "yali0_id", "uniprot_accession", "uniprot_protein_name", "length_aa",
          "uniprot_subcellular_location", "consensus_essential", "cas9_call", "cas9_fitness_score", "cas9_q_value",
          "cas12a_call", "cas12a_fitness_score", "cas12a_q_value", "annotation_function"]
fields += [f"ratio_{v}" for v in variants] + [f"class15_{v}" for v in variants]
with (TASK / "or_rule_review_genes.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    for gene, (rules, role) in ROLE.items():
        y0 = s2.get(gene, "")
        u = uni.get(y0, {})
        row = dict(gene_id=gene, rules=rules, role=role, yali0_id=y0, uniprot_accession=u.get("Entry", ""),
                   uniprot_protein_name=u.get("Protein names", ""), length_aa=u.get("Length", ""),
                   uniprot_subcellular_location=(u.get("Subcellular location [CC]", "") or "").split(" {")[0].replace("SUBCELLULAR LOCATION: ", ""),
                   consensus_essential=gene in positives, annotation_function=positives.get(gene, ""))
        for assay in ("Cas9", "Cas12a"):
            r = fit.get((gene, assay))
            row[f"{assay.lower()}_call"] = r["experimental_call"] if r else "not_in_screen_table"
            row[f"{assay.lower()}_fitness_score"] = f"{float(r['fitness_score']):.3f}" if r else ""
            row[f"{assay.lower()}_q_value"] = f"{float(r['q_value']):.3g}" if r else ""
        if ko:
            for v in variants:
                k = ko["gene_knockouts_by_variant"][v][gene]
                row[f"ratio_{v}"] = f"{k['ratio']:.6f}" if k["ratio"] is not None else k["status"]
                row[f"class15_{v}"] = k["classes"]["0.15"]
        w.writerow(row)
print("rows:", len(ROLE), "variants merged:", len(variants))
