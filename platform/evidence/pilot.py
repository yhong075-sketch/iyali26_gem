"""Select the M1 pilot reactions from `pilot_selection.yaml` and the essentiality case ledger.

Essentiality group: every passed and downgraded skeptic case up to the configured counts, plus
blocked cases drawn with a fixed seed among single-reaction cases. One primary reaction per case:
the first listed reaction that is in the model and not already chosen.
"""

import csv
import json
import random
from pathlib import Path

import yaml

FIELDS = ["reaction_id", "group", "origin", "note"]


def select(selection_path: Path, ledger_csv: Path, evidence_dir: Path, model_reactions: set[str]) -> list[dict]:
    cfg = yaml.safe_load(selection_path.read_text())
    rows: list[dict] = []
    taken: set[str] = set()

    def add(rid: str, group: str, origin: str, note: str) -> None:
        rows.append({"reaction_id": rid, "group": group, "origin": origin, "note": note})
        taken.add(rid)

    ess = cfg["essentiality"]
    with ledger_csv.open(newline="") as handle:
        ledger = {r["case_id"]: r for r in csv.DictReader(handle)}
    by_verdict: dict[str, list[str]] = {"pass": [], "downgrade": [], "block": []}
    for case_id in sorted(ledger):
        path = evidence_dir / f"{case_id}.json"
        verdict = (json.loads(path.read_text()).get("adversarial_review") or {}).get("verdict") if path.exists() else None
        if verdict in by_verdict:
            by_verdict[verdict].append(case_id)
    blocked = [c for c in by_verdict["block"]
               if len([x for x in ledger[c]["reaction_ids"].split(";") if x]) == 1]
    random.Random(ess["seed"]).shuffle(blocked)
    for verdict, pool, count in (("pass", by_verdict["pass"], ess["passed"]),
                                 ("downgrade", by_verdict["downgrade"], ess["downgraded"]),
                                 ("block", blocked, ess["blocked"])):
        chosen = 0
        for case_id in pool:
            if chosen == count:
                break
            rid = next((x for x in ledger[case_id]["reaction_ids"].split(";") if x in model_reactions
                        and x not in taken), None)
            if rid is None:
                continue
            add(rid, "essentiality", case_id, f"skeptic {verdict}; {ledger[case_id]['category']}")
            chosen += 1
    for group in ("phosphatidate_hub", "central_carbon"):
        for item in cfg[group]:
            if item["reaction_id"] not in model_reactions:
                raise ValueError(f"{item['reaction_id']} ({group}) is not in the model")
            if item["reaction_id"] in taken:
                raise ValueError(f"{item['reaction_id']} is selected twice")
            add(item["reaction_id"], group, "pilot_selection.yaml", item["note"])
    return rows


def write(rows: list[dict], path: Path) -> Path:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


def read(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))
