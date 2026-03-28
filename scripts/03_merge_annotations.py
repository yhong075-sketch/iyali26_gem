"""
Step 3: Merge annotations from BiGG, KEGG, and LLM sources.

Priority: BiGG > KEGG > LLM > embedded_formula (for formula field only).
Merges on base_id (unique metabolite identifier).

Input:  data/metabolites_unique.csv
        data/bigg_matches.csv
        data/kegg_matches.csv
        data/llm_annotations.json
Output: data/final_annotation_table.csv
"""

import pandas as pd
import json
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"


def get_col(df: pd.DataFrame, key_col: str, val_col: str,
            target_keys: pd.Series) -> pd.Series:
    """Extract values from df aligned to target_keys."""
    if val_col not in df.columns:
        return pd.Series([None] * len(target_keys), index=target_keys.index)
    valid = df.dropna(subset=[val_col])
    if valid.empty:
        return pd.Series([None] * len(target_keys), index=target_keys.index)
    mapping = valid.set_index(key_col)[val_col]
    return target_keys.map(mapping)


def first_valid(*series_list: pd.Series) -> pd.Series:
    """Return first non-null value across multiple Series."""
    result = series_list[0].copy()
    for s in series_list[1:]:
        result = result.fillna(s)
    return result


def main():
    unique_df = pd.read_csv(DATA_DIR / "metabolites_unique.csv")
    bigg_df   = pd.read_csv(DATA_DIR / "bigg_matches.csv")
    kegg_df   = pd.read_csv(DATA_DIR / "kegg_matches.csv")

    llm_path = DATA_DIR / "llm_annotations.json"
    if llm_path.exists():
        with open(llm_path) as f:
            llm_df = pd.DataFrame(json.load(f))
    else:
        print("Warning: llm_annotations.json not found, skipping LLM source")
        llm_df = pd.DataFrame(columns=["base_id"])

    # Build final table keyed on base_id
    final = unique_df[["base_id", "clean_name", "embedded_formula"]].copy()
    base_ids = final["base_id"]

    # ── Name (priority: BiGG > LLM > clean_name) ──
    final["name"] = first_valid(
        get_col(bigg_df, "base_id", "name", base_ids),
        get_col(llm_df,  "base_id", "name", base_ids),
        final["clean_name"],
    )

    # ── BiGG ID (priority: BiGG > LLM) ──
    final["bigg_id"] = first_valid(
        get_col(bigg_df, "base_id", "bigg_id", base_ids),
        get_col(llm_df,  "base_id", "bigg_id", base_ids),
    )

    # ── KEGG ID (priority: KEGG > LLM) ──
    final["kegg_id"] = first_valid(
        get_col(kegg_df, "base_id", "kegg_id", base_ids),
        get_col(llm_df,  "base_id", "kegg_id", base_ids),
    )

    # ── Formula (priority: BiGG > KEGG > LLM > embedded_formula) ──
    final["formula"] = first_valid(
        get_col(bigg_df, "base_id", "formula", base_ids),
        get_col(kegg_df, "base_id", "formula", base_ids),
        get_col(llm_df,  "base_id", "formula", base_ids),
        final["embedded_formula"],
    )

    # ── Charge (priority: BiGG > LLM) ──
    final["charge"] = first_valid(
        get_col(bigg_df, "base_id", "charge", base_ids),
        get_col(llm_df,  "base_id", "charge", base_ids),
    )

    # ── Confidence (from LLM, or "database" if from BiGG/KEGG) ──
    llm_conf = get_col(llm_df, "base_id", "confidence", base_ids)
    has_bigg = get_col(bigg_df, "base_id", "bigg_id", base_ids).notna()
    has_kegg = get_col(kegg_df, "base_id", "kegg_id", base_ids).notna()
    final["confidence"] = "unmatched"
    final.loc[llm_conf.notna(), "confidence"] = llm_conf[llm_conf.notna()]
    final.loc[has_kegg, "confidence"] = "database"
    final.loc[has_bigg, "confidence"] = "database"

    # Drop helper columns
    final = final.drop(columns=["clean_name", "embedded_formula"])

    # Save
    final.to_csv(DATA_DIR / "final_annotation_table.csv", index=False)

    # ── Summary ──
    total = len(final)
    print("Annotation coverage:")
    for col in ["name", "bigg_id", "kegg_id", "formula", "charge"]:
        n = final[col].notna().sum()
        print(f"  {col:10s}: {n:4d}/{total} ({n/total:.1%})")

    print(f"\nConfidence distribution:")
    print(final["confidence"].value_counts().to_string())

    print(f"\nSaved → data/final_annotation_table.csv ({total} rows)")


if __name__ == "__main__":
    main()
