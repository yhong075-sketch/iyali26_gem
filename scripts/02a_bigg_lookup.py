"""
Step 2a: Look up unique metabolites against BiGG database via REST API.

Features:
  - Resume support: skips already-queried base_ids if output CSV exists
  - Periodic checkpoint saves (every 50 queries)
  - Logging to both console and data/bigg_lookup.log

Input:  data/metabolites_unique.csv
Output: data/bigg_matches.csv
"""

import requests
import pandas as pd
import time
import logging
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
OUTPUT   = DATA_DIR / "bigg_matches.csv"
LOG_FILE = DATA_DIR / "bigg_lookup.log"

SEARCH_URL = "http://bigg.ucsd.edu/api/v2/search"
DETAIL_URL = "http://bigg.ucsd.edu/api/v2/universal/metabolites"

RATE_LIMIT  = 0.25   # seconds between requests
CHECKPOINT  = 50     # save every N queries
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


def query_bigg(clean_name: str, retries: int = MAX_RETRIES) -> dict | None:
    """Search BiGG by metabolite name, return best match or None."""
    params = {"query": clean_name, "search_type": "metabolites"}
    for attempt in range(retries):
        try:
            resp = requests.get(SEARCH_URL, params=params, timeout=15)
            if resp.status_code != 200:
                return None
            results = resp.json().get("results", [])
            if not results:
                return None
            bigg_id = results[0]["bigg_id"]
            detail = requests.get(
                f"{DETAIL_URL}/{bigg_id}", timeout=15
            ).json()
            return {
                "bigg_id":  bigg_id,
                "name":     detail.get("name"),
                "formula":  (detail.get("formulae") or [None])[0],
                "charge":   (detail.get("charges")  or [None])[0],
            }
        except requests.RequestException as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
            else:
                log.warning(f"Failed after {retries} retries for '{clean_name}': {e}")
                return None


def main():
    unique_df = pd.read_csv(DATA_DIR / "metabolites_unique.csv")
    clean_df = unique_df[unique_df["is_dirty"] == False].copy()

    # ── Resume: load existing results ──
    done_ids = set()
    results = []
    if OUTPUT.exists():
        existing = pd.read_csv(OUTPUT)
        done_ids = set(existing["base_id"])
        results = existing.to_dict("records")
        log.info(f"Resuming: {len(done_ids)} already queried")

    remaining = clean_df[~clean_df["base_id"].isin(done_ids)]
    log.info(f"BiGG lookup: {len(remaining)} to query, {len(done_ids)} done, "
             f"{len(unique_df) - len(clean_df)} dirty skipped")

    for i, (_, row) in enumerate(remaining.iterrows()):
        match = query_bigg(row["clean_name"])
        results.append({
            "base_id":  row["base_id"],
            "bigg_id":  match["bigg_id"] if match else None,
            "name":     match["name"]    if match else None,
            "formula":  match["formula"] if match else None,
            "charge":   match["charge"]  if match else None,
            "source":   "bigg"           if match else None,
        })
        time.sleep(RATE_LIMIT)

        if (i + 1) % CHECKPOINT == 0:
            pd.DataFrame(results).to_csv(OUTPUT, index=False)
            hits = sum(1 for r in results if r["bigg_id"] is not None)
            log.info(f"Checkpoint: {len(done_ids)+i+1}/{len(clean_df)} | Hits: {hits}")

    # ── Final save ──
    bigg_df = pd.DataFrame(results)
    bigg_df.to_csv(OUTPUT, index=False)

    total = len(bigg_df)
    hits  = bigg_df["bigg_id"].notna().sum()
    log.info(f"BiGG complete: {hits}/{total} matched ({hits/total:.1%})")


if __name__ == "__main__":
    main()
