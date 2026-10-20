"""Step 3: "literature" and "review_queue" sheets in review.xlsx."""

from pathlib import Path

from openpyxl import Workbook, load_workbook

from .models import EvidenceType, compartments_disagree

LITERATURE_COLUMNS = [
    "gene", "reaction_id", "reaction_name", "ec", "subsystem", "priority_tier", "pair_hits",
    "fulltext_available", "status", "best_evidence_level", "best_evidence_name",
    "no_organism_specific_evidence", "n_claims", "n_supporting", "claims_disagree", "pmids_cited", "quotes",
    "shared_with", "input_tokens", "output_tokens",
]
QUEUE_COLUMNS = ["priority_tier", "reason", "gene", "reaction_id", "reaction_name", "ec", "subsystem",
                 "pair_hits", "best_evidence_level", "n_claims", "n_supporting", "pmids_cited"]


def _disagree(claims: list[dict]) -> bool:
    verdicts = {c["supports_gene_reaction"] for c in claims}
    return len(verdicts) > 1 or compartments_disagree(c.get("compartment") for c in claims)


def literature_rows(retrievals: list[dict], results: list[dict]) -> list[dict]:
    by_pair = {(r["gene"], r["reaction_id"]): r for r in results}
    rows = []
    for rec in retrievals:
        res = by_pair.get((rec["gene"], rec["reaction_id"]))
        review = (res or {}).get("review") or {}
        claims = review.get("claims", [])
        level = review.get("best_evidence_level")
        rows.append({
            "gene": rec["gene"], "reaction_id": rec["reaction_id"], "reaction_name": rec["reaction_name"],
            "ec": ";".join(rec["ec_codes"]), "subsystem": rec["subsystem"], "priority_tier": rec.get("tier"),
            "pair_hits": rec["pair_hits"], "fulltext_available": "Y" if rec["fulltext_available"] else "N",
            "status": res["status"] if res else "not_extracted",
            "best_evidence_level": level,
            "best_evidence_name": EvidenceType(level).name if level is not None else None,
            "no_organism_specific_evidence": review.get("no_organism_specific_evidence"),
            "n_claims": len(claims), "n_supporting": sum(c["supports_gene_reaction"] for c in claims),
            "claims_disagree": _disagree(claims) if claims else False,
            "pmids_cited": ";".join(dict.fromkeys(c["pmid"] for c in claims)),
            "quotes": " | ".join(f'{c["pmid"]}: "{c["quote"]}"' for c in claims),
            "shared_with": (res or {}).get("shared_with"),
            "input_tokens": (res or {}).get("usage", {}).get("input_tokens"),
            "output_tokens": (res or {}).get("usage", {}).get("output_tokens"),
        })
    return rows


def queue_reason(row: dict) -> str | None:
    if row["status"] in {"EXTRACTION_FAILED", "MODEL_ERROR", "not_extracted"}:
        return row["status"]
    if row["claims_disagree"]:
        return "claims_disagree"
    if row["n_claims"] and not row["n_supporting"]:
        return "contradicted"   # every validated claim says the gene does NOT catalyse this reaction
    if row["status"] == "no_hits":
        return "no_hits"
    if row["best_evidence_level"] is not None and row["best_evidence_level"] <= EvidenceType.homology:
        return "low_evidence"
    return None


def review_queue_rows(rows: list[dict]) -> list[dict]:
    order = {"EXTRACTION_FAILED": 0, "MODEL_ERROR": 0, "not_extracted": 0, "contradicted": 1,
             "claims_disagree": 2, "low_evidence": 3, "no_hits": 4}
    queue = []
    for row in rows:
        reason = queue_reason(row)
        if reason:
            queue.append({**{k: row[k] for k in QUEUE_COLUMNS if k != "reason"}, "reason": reason})
    queue.sort(key=lambda r: (r["priority_tier"] if r["priority_tier"] is not None else 99,
                              order[r["reason"]], r["reaction_id"], r["gene"]))
    return queue


def write_review_xlsx(path: Path, rows: list[dict], queue: list[dict]) -> Path:
    """Add or replace the two sheets; other sheets of an existing workbook are kept."""
    workbook = load_workbook(path) if path.exists() else Workbook()
    if not path.exists():
        workbook.remove(workbook.active)
    for name, columns, data in (("literature", LITERATURE_COLUMNS, rows),
                                ("review_queue", QUEUE_COLUMNS, queue)):
        if name in workbook.sheetnames:
            workbook.remove(workbook[name])
        sheet = workbook.create_sheet(name)
        sheet.append(columns)
        for row in data:
            sheet.append([row.get(column) for column in columns])
        sheet.freeze_panes = "A2"
    path.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(path)
    return path
