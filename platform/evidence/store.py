"""SQLite sidecar (spec §7). Rebuilt from scratch on every build, rows inserted in sorted order, so
the content hash depends only on the cache, committed code and the weights YAML (D1)."""

import hashlib
import json
import sqlite3
from pathlib import Path

from .pipeline import Bundle
from .schemas import Assertion, Claim, PairRelation, Passage, RejectedClaim, Source, Verdict

TABLES = ("sources", "passages", "claims", "rejected_claims", "pair_relations", "assertions", "verdicts",
          "decisions", "reviewer_runs", "reviewer_disagreements", "extraction_runs", "id_map", "weights_version",
          "load_bearing")


def _rows(bundle: Bundle, id_map_rows: list[dict], weights_version: str,
          extraction_runs: list[dict]) -> dict[str, list[tuple[str, str]]]:
    dump = lambda m: json.dumps(m.model_dump(mode="json"), sort_keys=True)  # noqa: E731
    return {
        "sources": [(k, dump(v)) for k, v in sorted(bundle.sources.items())],
        "passages": [(k, dump(v)) for k, v in sorted(bundle.passages.items())],
        "claims": [(k, dump(v)) for k, v in sorted(bundle.claims.items())],
        "rejected_claims": [(r.reject_id, dump(r)) for r in bundle.rejected],
        "pair_relations": [(k, dump(v)) for k, v in sorted(bundle.relations.items())],
        "assertions": [(k, dump(v)) for k, v in sorted(bundle.assertions.items())],
        "verdicts": [(k, dump(v)) for k, v in sorted(bundle.verdicts.items())],
        "decisions": [],
        "reviewer_runs": [],
        "reviewer_disagreements": [],
        "extraction_runs": [(r["extraction_run_id"], json.dumps(r, sort_keys=True)) for r in extraction_runs],
        "id_map": [(r["model_gene"], json.dumps(r, sort_keys=True)) for r in id_map_rows],
        "weights_version": [("current", json.dumps({"version": weights_version}))],
        "load_bearing": [(k, json.dumps(list(v))) for k, v in sorted(bundle.load_bearing.items())],
    }


def content_hash(rows: dict[str, list[tuple[str, str]]]) -> str:
    digest = hashlib.sha256()
    for table in TABLES:
        for key, payload in rows[table]:
            digest.update(f"{table}\0{key}\0{payload}\n".encode())
    return digest.hexdigest()


def write(path: Path, bundle: Bundle, id_map_rows: list[dict], weights_version: str,
          extraction_runs: list[dict]) -> str:
    """Write a fresh database and return its content hash."""
    rows = _rows(bundle, id_map_rows, weights_version, extraction_runs)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.unlink(missing_ok=True)
    con = sqlite3.connect(tmp)
    with con:
        for table in TABLES:
            con.execute(f"CREATE TABLE {table} (id TEXT PRIMARY KEY, payload TEXT NOT NULL)")
            con.executemany(f"INSERT INTO {table} VALUES (?, ?)", rows[table])
    con.close()
    tmp.replace(path)
    return content_hash(rows)


def load(path: Path) -> Bundle:
    con = sqlite3.connect(path)
    try:
        def table(name: str, model):
            return {k: model.model_validate_json(v) for k, v in con.execute(f"SELECT id, payload FROM {name}")}
        bundle = Bundle(sources=table("sources", Source), passages=table("passages", Passage),
                        claims=table("claims", Claim), rejected=list(table("rejected_claims", RejectedClaim).values()),
                        relations=table("pair_relations", PairRelation), assertions=table("assertions", Assertion),
                        verdicts=table("verdicts", Verdict))
        bundle.load_bearing = {k: tuple(json.loads(v)) for k, v in con.execute("SELECT id, payload FROM load_bearing")}
    finally:
        con.close()
    return bundle
