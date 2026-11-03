"""Audit the W29 RefSeq protein filed under each model gene's YALI1 locus tag against the S2 locus map.

The S2 table of the W29/CLIB89 assembly (Magnan et al. 2016, PLoS ONE 11:e0162363) gives each YALI1 tag its
position and strand, and pairs it with the CLIB122 YALI0 gene at that position. For some tags the CDS that
GenBank/RefSeq files under the same tag lies on the opposite strand: an antisense ORF, so any sequence work keyed on
the tag (BLAST, localisation predictors, AlphaFold) reads the wrong protein. This tool compares, for every gene in
the start model and in every candidate model, the S2 row with the RefSeq CDS row of the NCBI feature table.

Read-only: writes the audit CSV and its JSON sidecar, never a model. Run from platform/:
    python -m tools.audit_w29_tag_protein [--research-root PATH]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Iterable

from scripts.gem_annotate.config import load_project_paths
from scripts.gem_annotate.locus_resolver import model_gene_fingerprint
from scripts.gem_annotate.model_layout import MODEL

AUDIT_NAME = "w29_tag_protein_audit.csv"
SIDECAR_NAME = "w29_tag_protein_audit.json"
DECISIONS_NAME = "w29_tag_protein_decisions.csv"
S2_RELATIVE = Path("reference/locus_map/S2_table_YALI1_YALI0_map.xlsx")
# Also holds the 56 YALI1M rows of the "Mito Table" sheet, with the same tags, types and spans
# (only the case of some YALI0 names differs, e.g. YALIfMr30 vs YalifMr30).
S2_SHEET = "YALI1 Genes"
FEATURE_TABLE_RELATIVE = Path("reference/ncbi/clib89_feature_table.txt")

# Flags that mark a disagreement and need a row in the decisions file, with the decisions each one allows.
DECISIONS_BY_FLAG = {
    "strand_flip": ("refseq_protein_is_antisense_orf",),
    "no_overlap": ("refseq_protein_kept",),
    "no_s2_row": ("model_id_differs_from_s2_tag",),
    "no_refseq_cds": ("no_w29_protein",),
}
NEEDS_DECISION = tuple(DECISIONS_BY_FLAG)
# Allowed for any flag when curated records disagree; names no protein until the user resolves it.
UNRESOLVED = "unresolved_conflict"
FLAGS = ("ok", *NEEDS_DECISION, "mito_absent_from_feature_table", "not_protein_coding", "not_yali1")
AUDIT_FIELDS = ("gene", "yali1_tag", "yali0", "s2_gene_type", "s2_span", "refseq_protein", "refseq_span", "flag")
# refseq_protein and refseq_span record the audit row a decision was made for; a changed row invalidates it.
DECISION_FIELDS = ("schema_version", "gene", "flag", "refseq_protein", "refseq_span", "decision", "sequence_source",
                   "sequence_accession", "evidence", "reason")
SEQUENCE_SOURCES = ("UniProt", "RefSeq", "none")

_YALI1_TAG_RE = re.compile(r"^YALI1_?([A-FM])(\d+)([a-z])$")
# cobra reads a gene's ID from fbc:id minus the "G_" prefix; fbc:label carries the prefix in cobra-written files.
_GENE_ID_RE = re.compile(rb'<fbc:geneProduct\b[^>]*\bfbc:id="(?:G_)?([^"]+)"')


def yali1_tag(gene_id: str) -> str | None:
    """Return the S2/NCBI spelling (YALI1_X00000g) of a YALI1 model gene ID, else None."""
    match = _YALI1_TAG_RE.fullmatch(gene_id)
    return f"YALI1_{match.group(1)}{match.group(2)}{match.group(3)}" if match else None


def model_files() -> list[Path]:
    """The start model and every candidate model XML of this model folder."""
    return [MODEL.start_model, *sorted(MODEL.candidates.rglob("*.xml"))]


def model_gene_ids(paths: Iterable[Path] | None = None) -> set[str]:
    """Union of geneProduct IDs across the model files, read as text (no SBML parser, no cobra)."""
    genes: set[str] = set()
    for path in paths if paths is not None else model_files():
        found = {gene_id.decode() for gene_id in _GENE_ID_RE.findall(path.read_bytes())}
        if not found:
            raise ValueError(f"No fbc:geneProduct IDs found in {path}")
        genes |= found
    return genes


def read_s2(path: Path) -> dict[str, dict]:
    """S2 rows by YALI1 tag: YALI0 partner, gene type, chromosome, span and strand."""
    import openpyxl

    sheet = openpyxl.load_workbook(path, read_only=True)[S2_SHEET]
    rows = sheet.iter_rows(values_only=True)
    header = next(rows)
    expected = ("YALI0 locus ID", "YALI1 Locus ID", "Gene Type", "Chromosome", "Start Base", "End Base", "Strand")
    if tuple(header[: len(expected)]) != expected:
        raise ValueError(f"Unexpected S2 header in {path}: {header}")
    s2: dict[str, dict] = {}
    for yali0, tag, gene_type, chrom, start, end, strand, *_ in rows:
        if tag in s2:
            raise ValueError(f"Duplicate S2 row for {tag}")
        s2[tag] = dict(yali0=yali0 or "", gene_type=gene_type, chrom=chrom, start=int(start), end=int(end),
                       strand=strand)
    return s2


def read_feature_table(path: Path) -> dict[str, dict]:
    """RefSeq CDS rows by locus tag: protein accession, chromosome, span and strand."""
    cds: dict[str, dict] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["# feature"] != "CDS":
                continue
            tag = row["locus_tag"]
            if tag in cds:
                raise ValueError(f"Duplicate CDS row for {tag} in {path}")
            cds[tag] = dict(protein=row["product_accession"], chrom=f"YALI{row['chromosome']}",
                            start=int(row["start"]), end=int(row["end"]), strand=row["strand"])
    return cds


def classify(gene_id: str, s2_row: dict | None, cds_row: dict | None) -> str:
    """One flag per gene; the order of the tests decides which disagreement is reported."""
    if yali1_tag(gene_id) is None:
        return "not_yali1"
    if s2_row is None:
        return "no_s2_row"
    if cds_row is None:
        if s2_row["gene_type"] != "mRNA":
            return "not_protein_coding"
        return "mito_absent_from_feature_table" if s2_row["chrom"] == "YALI1M" else "no_refseq_cds"
    # Any tag with a RefSeq CDS is compared, whatever its S2 gene type ("Pseudo mRNA", "RetroElement", ...).
    if s2_row["strand"] != cds_row["strand"]:
        return "strand_flip"
    if s2_row["chrom"] != cds_row["chrom"] or min(s2_row["end"], cds_row["end"]) < max(s2_row["start"], cds_row["start"]):
        return "no_overlap"
    return "ok"


def _span(row: dict | None) -> str:
    return f"{row['chrom']}:{row['start']}-{row['end']}{row['strand']}" if row else ""


def build_audit(gene_ids: Iterable[str], s2: dict[str, dict], cds: dict[str, dict]) -> list[dict]:
    rows = []
    for gene_id in sorted(gene_ids):
        tag = yali1_tag(gene_id)
        s2_row = s2.get(tag) if tag else None
        cds_row = cds.get(tag) if tag else None
        rows.append(dict(gene=gene_id, yali1_tag=tag or "", yali0=s2_row["yali0"] if s2_row else "",
                         s2_gene_type=s2_row["gene_type"] if s2_row else "", s2_span=_span(s2_row),
                         refseq_protein=cds_row["protein"] if cds_row else "", refseq_span=_span(cds_row),
                         flag=classify(gene_id, s2_row, cds_row)))
    return rows


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_path(name: str = AUDIT_NAME) -> Path:
    """The committed file in model/curation/annotation/ (found by unique name once it exists)."""
    path = MODEL.curation_file(name)
    return path if path.exists() else MODEL.curation / "annotation" / name


def read_audit(path: Path | None = None) -> list[dict]:
    with (path or audit_path()).open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_decisions(path: Path | None = None) -> list[dict]:
    with (path or audit_path(DECISIONS_NAME)).open(newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != DECISION_FIELDS:
            raise ValueError(f"Unexpected decision columns: {reader.fieldnames}")
        return list(reader)


def check_decisions(audit_rows: list[dict], decisions: list[dict]) -> list[str]:
    """Every disagreement has exactly one valid decision, and every decision matches a current disagreement."""
    errors: list[str] = []
    by_gene: dict[str, dict] = {}
    for row in decisions:
        gene = row["gene"]
        if gene in by_gene:
            errors.append(f"{gene}: duplicate decision")
        by_gene[gene] = row
        if row["schema_version"] != "1":
            errors.append(f"{gene}: unsupported schema_version {row['schema_version']!r}")
        allowed = (*DECISIONS_BY_FLAG[row["flag"]], UNRESOLVED) if row["flag"] in DECISIONS_BY_FLAG else ()
        if row["decision"] not in allowed:
            errors.append(f"{gene}: decision {row['decision']!r} is not allowed for flag {row['flag']!r}")
        if row["decision"] == UNRESOLVED and row["sequence_source"] != "none":
            errors.append(f"{gene}: an unresolved conflict must not name a protein")
        if row["sequence_source"] not in SEQUENCE_SOURCES:
            errors.append(f"{gene}: unknown sequence_source {row['sequence_source']!r}")
        if row["sequence_source"] != "none" and not row["sequence_accession"]:
            errors.append(f"{gene}: sequence_accession is required for {row['sequence_source']}")
        if row["decision"] == "refseq_protein_is_antisense_orf" and row["sequence_source"] != "UniProt":
            errors.append(f"{gene}: an antisense RefSeq protein must point to the CLIB122 UniProt protein")
        if row["decision"] == "refseq_protein_kept" and (
                row["sequence_source"] != "RefSeq" or row["sequence_accession"] != row["refseq_protein"]):
            errors.append(f"{gene}: a kept RefSeq protein must name that RefSeq protein")
        if not row["evidence"] or not row["reason"]:
            errors.append(f"{gene}: evidence and reason are required")
    flagged = {row["gene"]: row for row in audit_rows if row["flag"] in NEEDS_DECISION}
    for gene, audit in sorted(flagged.items()):
        decision = by_gene.get(gene)
        if decision is None:
            errors.append(f"{gene}: flag {audit['flag']!r} has no decision")
        elif decision["flag"] != audit["flag"]:
            errors.append(f"{gene}: decision is for flag {decision['flag']!r}, audit says {audit['flag']!r}")
        elif (decision["refseq_protein"], decision["refseq_span"]) != (audit["refseq_protein"], audit["refseq_span"]):
            errors.append(f"{gene}: decision was made for {decision['refseq_protein'] or '-'} "
                          f"{decision['refseq_span'] or '-'}, audit now shows {audit['refseq_protein'] or '-'} "
                          f"{audit['refseq_span'] or '-'}")
    for gene in sorted(set(by_gene) - set(flagged)):
        errors.append(f"{gene}: stale decision, the audit no longer flags this gene")
    return errors


def audit_from_research_root(research_root: str | Path | None = None) -> tuple[list[dict], dict]:
    """Rebuild the audit rows and sidecar from the raw research-workspace inputs (fails closed)."""
    paths = load_project_paths(research_root, required=True)
    s2_path = paths.research_root / S2_RELATIVE
    feature_path = paths.research_root / FEATURE_TABLE_RELATIVE
    paths.require(s2_path, feature_path)
    genes = model_gene_ids()
    rows = build_audit(genes, read_s2(s2_path), read_feature_table(feature_path))
    counts = {flag: sum(row["flag"] == flag for row in rows) for flag in FLAGS}
    sidecar = {
        "generator": "platform/tools/audit_w29_tag_protein.py",
        "inputs": {
            "s2_table": {"path": str(S2_RELATIVE), "sheet": S2_SHEET, "sha256": file_sha256(s2_path)},
            "feature_table": {"path": str(FEATURE_TABLE_RELATIVE), "sha256": file_sha256(feature_path)},
        },
        "model_gene_count": len(genes),
        "model_gene_fingerprint": model_gene_fingerprint(genes),
        "flag_counts": counts,
    }
    return rows, sidecar


def write_audit(rows: list[dict], sidecar: dict) -> tuple[Path, Path]:
    csv_path, json_path = audit_path(AUDIT_NAME), audit_path(SIDECAR_NAME)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(sidecar, indent=2) + "\n")
    return csv_path, json_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--research-root", help="research workspace (default: $IYALI26_RESEARCH_ROOT)")
    args = parser.parse_args()
    rows, sidecar = audit_from_research_root(args.research_root)
    csv_path, json_path = write_audit(rows, sidecar)
    print(json.dumps(sidecar["flag_counts"], indent=2))
    for row in rows:
        if row["flag"] in NEEDS_DECISION:
            print(row["flag"], row["gene"], row["yali0"], row["s2_span"], row["refseq_protein"], row["refseq_span"])
    print(f"wrote {csv_path} and {json_path}")
    decisions_path = audit_path(DECISIONS_NAME)
    errors = check_decisions(rows, read_decisions(decisions_path)) if decisions_path.exists() else [
        f"{decisions_path} is missing"]
    for error in errors:
        print("DECISION GAP:", error)
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
