"""
Step 2b: Look up metabolites against KEGG database.

If BiGG had hits, only processes BiGG misses. Otherwise queries all non-dirty.
Features:
  - Resume support: skips already-queried base_ids
  - Periodic checkpoint saves
  - Name search + formula fallback
  - Logging to data/kegg_lookup.log

Input:  data/metabolites_unique.csv, data/bigg_matches.csv (optional)
Output: data/kegg_matches.csv
"""

import requests
import pandas as pd
import time
import logging
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
OUTPUT   = DATA_DIR / "kegg_matches.csv"
LOG_FILE = DATA_DIR / "kegg_lookup.log"

RATE_LIMIT  = 0.35   # seconds between requests (KEGG is stricter)
CHECKPOINT  = 50
MAX_RETRIES = 3

# ── Logging ──
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)


def _kegg_get(url: str, retries: int = MAX_RETRIES) -> str | None:
    """GET with retries."""
    for attempt in range(retries):
        try:
            resp = requests.get(url, timeout=15)
            if resp.status_code == 200:
                return resp.text
            return None
        except requests.RequestException as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
            else:
                log.warning(f"Failed {url}: {e}")
                return None


def _get_kegg_formula(kegg_id: str) -> str | None:
    """Fetch molecular formula from KEGG compound detail."""
    text = _kegg_get(f"https://rest.kegg.jp/get/{kegg_id}")
    if not text:
        return None
    for line in text.split("\n"):
        if line.startswith("FORMULA"):
            return line.split()[-1]
    return None


def query_kegg_by_name(name: str) -> dict | None:
    """Search KEGG compound database by name."""
    text = _kegg_get(f"https://rest.kegg.jp/find/compound/{name}")
    if not text or not text.strip():
        return None

    first_line = text.strip().split("\n")[0]
    parts = first_line.split("\t")
    if len(parts) < 2:
        return None

    kegg_id = parts[0].replace("cpd:", "")
    kegg_names = parts[1]
    formula = _get_kegg_formula(kegg_id)

    return {
        "kegg_id":   kegg_id,
        "kegg_name": kegg_names.split(";")[0].strip(),
        "formula":   formula,
    }


def query_kegg_by_formula(formula: str) -> dict | None:
    """Search KEGG by molecular formula (fallback)."""
    text = _kegg_get(f"https://rest.kegg.jp/find/compound/{formula}/formula")
    if not text or not text.strip():
        return None

    first_line = text.strip().split("\n")[0]
    parts = first_line.split("\t")
    if len(parts) < 2:
        return None

    kegg_id = parts[0].replace("cpd:", "")
    kegg_names = parts[1]

    return {
        "kegg_id":   kegg_id,
        "kegg_name": kegg_names.split(";")[0].strip(),
        "formula":   formula,
    }


def main():
    unique_df = pd.read_csv(DATA_DIR / "metabolites_unique.csv")

    # Determine which metabolites to query
    bigg_path = DATA_DIR / "bigg_matches.csv"
    bigg_hits = set()
    if bigg_path.exists():
        bigg_df = pd.read_csv(bigg_path)
        bigg_hits = set(bigg_df[bigg_df["bigg_id"].notna()]["base_id"])

    if bigg_hits:
        to_search = unique_df[
            (~unique_df["base_id"].isin(bigg_hits)) &
            (unique_df["is_dirty"] == False)
        ].copy()
        log.info(f"Querying KEGG for {len(to_search)} BiGG-unmatched metabolites")
    else:
        to_search = unique_df[unique_df["is_dirty"] == False].copy()
        log.info(f"BiGG had 0 hits — querying KEGG for all {len(to_search)} non-dirty metabolites")

    # ── Resume: load existing results ──
    done_ids = set()
    results = []
    if OUTPUT.exists():
        existing = pd.read_csv(OUTPUT)
        done_ids = set(existing["base_id"])
        results = existing.to_dict("records")
        log.info(f"Resuming: {len(done_ids)} already queried")

    remaining = to_search[~to_search["base_id"].isin(done_ids)]
    log.info(f"KEGG lookup: {len(remaining)} to query, {len(done_ids)} done")

    for i, (_, row) in enumerate(remaining.iterrows()):
        # Try name search first
        match = query_kegg_by_name(row["clean_name"])

        # Fallback: try formula search if name fails and we have a formula
        if match is None and pd.notna(row.get("embedded_formula")):
            time.sleep(RATE_LIMIT)
            match = query_kegg_by_formula(row["embedded_formula"])

        results.append({
            "base_id":   row["base_id"],
            "kegg_id":   match["kegg_id"]   if match else None,
            "kegg_name": match["kegg_name"] if match else None,
            "formula":   match["formula"]   if match else None,
            "source":    "kegg"             if match else None,
        })
        time.sleep(RATE_LIMIT)

        if (i + 1) % CHECKPOINT == 0:
            pd.DataFrame(results).to_csv(OUTPUT, index=False)
            hits = sum(1 for r in results if r["kegg_id"] is not None)
            log.info(f"Checkpoint: {len(done_ids)+i+1}/{len(to_search)} | Hits: {hits}")

    # ── Final save ──
    kegg_df = pd.DataFrame(results)
    kegg_df.to_csv(OUTPUT, index=False)

    total = len(kegg_df)
    hits  = kegg_df["kegg_id"].notna().sum()
    log.info(f"KEGG complete: {hits}/{total} matched ({hits/total:.1%})")


if __name__ == "__main__":
    main()
