"""
Step 2c: LLM-assisted annotation for metabolites not matched by BiGG or KEGG.

Features:
  - Resume support: loads existing llm_annotations.json and skips done base_ids
  - Checkpoint saves after each batch
  - Retry on API errors
  - Logging to data/llm_annotate.log

Input:  data/metabolites_unique.csv, data/bigg_matches.csv, data/kegg_matches.csv
Output: data/llm_annotations.json
        data/manual_review.csv (low-confidence entries)
"""

import anthropic
import pandas as pd
import json
import time
import logging
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
OUTPUT   = DATA_DIR / "llm_annotations.json"
LOG_FILE = DATA_DIR / "llm_annotate.log"

BATCH_SIZE  = 20
MAX_RETRIES = 3
RETRY_DELAY = 10  # seconds

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

client = anthropic.Anthropic()  # reads ANTHROPIC_API_KEY from env

SYSTEM_PROMPT = """You are a metabolomics expert with deep knowledge of BiGG,
KEGG, and ModelSEED databases. Given raw metabolite names from a Yarrowia
lipolytica genome-scale model (iYli21), return standardized annotations.

Each input has:
- "base_id": model internal ID (numeric, no semantic meaning, e.g. "m1")
- "clean_name": chemical name extracted from raw name
- "embedded_formula": molecular formula extracted from raw name (may be null or incorrect)

ALWAYS respond with a JSON array. Each element must have:
{
  "base_id": "original base ID",
  "name": "standard human-readable name",
  "bigg_id": "BiGG universal metabolite ID (e.g. lac__L) or null",
  "kegg_id": "KEGG compound ID (e.g. C00186) or null",
  "formula": "correct molecular formula (verify/correct embedded_formula) or null",
  "charge": integer charge at physiological pH or null,
  "confidence": "high" | "medium" | "low",
  "notes": "brief reasoning for your annotation"
}

Guidelines:
- For common metabolites (ATP, pyruvate, etc.), provide high-confidence annotations
- Verify embedded_formula against known databases; correct if wrong
- If you cannot confidently identify the metabolite, set confidence to "low"
- BiGG IDs should be universal IDs (without compartment suffix)
- KEGG IDs should be compound IDs starting with C

Return ONLY the JSON array, no preamble or explanation."""


def annotate_batch(batch: list[dict]) -> list[dict]:
    """Send a batch to Claude API and parse JSON response."""
    user_msg = "Annotate these metabolites:\n" + json.dumps(batch, indent=2)

    for attempt in range(MAX_RETRIES):
        try:
            response = client.messages.create(
                model="claude-sonnet-4-20250514",
                max_tokens=4096,
                system=SYSTEM_PROMPT,
                messages=[{"role": "user", "content": user_msg}]
            )
            text = response.content[0].text.strip()
            # Clean markdown fences if present
            if text.startswith("```"):
                text = text.split("\n", 1)[1] if "\n" in text else text[3:]
            if text.endswith("```"):
                text = text[:-3]
            text = text.strip()
            return json.loads(text)

        except json.JSONDecodeError as e:
            log.warning(f"JSON parse error (attempt {attempt+1}): {e}")
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
            else:
                raise
        except anthropic.APIError as e:
            log.warning(f"API error (attempt {attempt+1}): {e}")
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY * (attempt + 1))
            else:
                raise


def _save_checkpoint(all_results: list[dict]):
    """Save current results to JSON."""
    with open(OUTPUT, "w") as f:
        json.dump(all_results, f, indent=2, ensure_ascii=False)


def main():
    bigg_df   = pd.read_csv(DATA_DIR / "bigg_matches.csv")
    kegg_df   = pd.read_csv(DATA_DIR / "kegg_matches.csv")
    unique_df = pd.read_csv(DATA_DIR / "metabolites_unique.csv")

    # Find unmatched metabolites
    matched_ids = set()
    if "bigg_id" in bigg_df.columns:
        matched_ids |= set(bigg_df[bigg_df["bigg_id"].notna()]["base_id"])
    if "kegg_id" in kegg_df.columns:
        matched_ids |= set(kegg_df[kegg_df["kegg_id"].notna()]["base_id"])

    remaining = unique_df[
        (~unique_df["base_id"].isin(matched_ids)) &
        (unique_df["is_dirty"] == False)
    ].copy()

    # ── Resume: load existing results ──
    all_results = []
    done_ids = set()
    if OUTPUT.exists():
        with open(OUTPUT) as f:
            all_results = json.load(f)
        done_ids = {r["base_id"] for r in all_results}
        log.info(f"Resuming: {len(done_ids)} already annotated")

    remaining = remaining[~remaining["base_id"].isin(done_ids)]
    log.info(f"LLM annotation: {len(remaining)} to annotate, "
             f"{len(done_ids)} done (batch size: {BATCH_SIZE})")

    for i in range(0, len(remaining), BATCH_SIZE):
        batch = remaining.iloc[i:i+BATCH_SIZE][
            ["base_id", "clean_name", "embedded_formula"]
        ].to_dict("records")

        try:
            results = annotate_batch(batch)
            all_results.extend(results)
            batch_num = i // BATCH_SIZE + 1
            batch_end = min(i + BATCH_SIZE, len(remaining))
            log.info(f"Batch {batch_num}: {batch_end}/{len(remaining)} done")
        except Exception as e:
            log.error(f"Batch {i//BATCH_SIZE + 1} failed: {e}")
            # Add fallback entries for failed batch
            for item in batch:
                all_results.append({
                    "base_id": item["base_id"],
                    "name": item["clean_name"],
                    "bigg_id": None, "kegg_id": None,
                    "formula": item.get("embedded_formula"),
                    "charge": None,
                    "confidence": "low",
                    "notes": "LLM batch failed - needs manual review"
                })

        # Checkpoint after every batch
        _save_checkpoint(all_results)
        time.sleep(1)  # rate limit between batches

    # ── Final outputs ──
    _save_checkpoint(all_results)

    llm_df = pd.DataFrame(all_results)
    low_conf = llm_df[llm_df.get("confidence", pd.Series()) == "low"]
    low_conf.to_csv(DATA_DIR / "manual_review.csv", index=False)

    log.info(f"LLM annotation complete: {len(llm_df)} total")
    if "confidence" in llm_df.columns:
        log.info(f"  High: {(llm_df['confidence']=='high').sum()}, "
                 f"Medium: {(llm_df['confidence']=='medium').sum()}, "
                 f"Low: {(llm_df['confidence']=='low').sum()}")
    log.info(f"Saved: llm_annotations.json, manual_review.csv ({len(low_conf)} entries)")


if __name__ == "__main__":
    main()
