"""
Step 1: Extract metabolites from SBML model.

- Parses raw_name → clean_name + embedded_formula
- Flags dirty entries (ActiveX VT_ERROR, etc.)
- Deduplicates by base_id (met ID without compartment suffix)

Input:  data/iyli21.xml
Output: data/metabolites_raw.csv      (all species)
        data/metabolites_unique.csv   (deduplicated unique metabolites)
        data/dirty_entries.csv        (flagged dirty entries)
"""

import re
from pathlib import Path

import cobra
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"

# Import compartment normalization from update_model
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from update_model import normalize_compartments

# ── Load and normalize model ──
MODEL = DATA_DIR / "iyli21.xml"
print(f"Loading model: {MODEL}")
model = cobra.io.read_sbml_model(str(MODEL))
model = normalize_compartments(model)
print(f"  Species: {len(model.metabolites)}")
print(f"  Reactions: {len(model.reactions)}")

# ── Extract all species ──
rows = []
for met in model.metabolites:
    rows.append({
        "met_id":      met.id,
        "raw_name":    met.name,
        "formula":     met.formula if met.formula else None,
        "charge":      met.charge,
        "compartment": met.compartment,
    })

raw_df = pd.DataFrame(rows)

# ── Parse name: split into clean_name + embedded_formula ──
# Common formats:
#   "(R)-lactate_C3H6O3"       → clean_name="(R)-lactate", formula="C3H6O3"
#   "H+_p+1"                   → clean_name="H+", formula="p+1" (not a real formula)
#   "ferricytochrome c_"       → clean_name="ferricytochrome c", formula=None
#   "...phosphate_ActiveX VT_ERROR:" → dirty entry

FORMULA_RE = re.compile(r'^(.+?)_([A-Z][A-Za-z0-9]+)$')
DIRTY_SUFFIX_RE = re.compile(r'_ActiveX.*$', re.IGNORECASE)


def parse_name(raw_name: str) -> tuple[str, str | None]:
    """Split raw_name into (clean_name, embedded_formula)."""
    if not isinstance(raw_name, str):
        return str(raw_name), None

    if raw_name.endswith("_"):
        return raw_name.rstrip("_").strip(), None

    # Strip Excel/ActiveX error suffixes before any other processing
    cleaned = DIRTY_SUFFIX_RE.sub("", raw_name).strip()
    if cleaned != raw_name:
        return cleaned, None

    m = FORMULA_RE.match(raw_name)
    if m:
        return m.group(1).strip(), m.group(2)

    return raw_name.strip(), None


parsed = raw_df["raw_name"].apply(lambda x: pd.Series(parse_name(x)))
raw_df["clean_name"] = parsed[0]
raw_df["embedded_formula"] = parsed[1]

# ── Flag dirty entries ──
DIRTY_PATTERNS = ["VT_ERROR", "ActiveX", "#REF", "#VALUE", "N/A"]
raw_df["is_dirty"] = raw_df["raw_name"].apply(
    lambda x: any(p in str(x) for p in DIRTY_PATTERNS)
)

dirty_df = raw_df[raw_df["is_dirty"]].copy()
dirty_df.to_csv(DATA_DIR / "dirty_entries.csv", index=False)
print(f"\nDirty entries: {len(dirty_df)} → data/dirty_entries.csv")
if len(dirty_df) > 0:
    for _, row in dirty_df.head(5).iterrows():
        print(f"  {row['met_id']}: {row['raw_name']}")

# ── Deduplicate: extract base_id ──
# Species IDs look like "m1_c", "m1_m" — base_id strips the compartment suffix
def extract_base_id(met_id: str) -> str:
    parts = met_id.rsplit("_", 1)
    if len(parts) == 2:
        return parts[0]
    return met_id

raw_df["base_id"] = raw_df["met_id"].apply(extract_base_id)

unique_df = raw_df.drop_duplicates(subset="base_id").copy()
unique_df = unique_df[["base_id", "clean_name", "embedded_formula", "is_dirty"]].reset_index(drop=True)

# ── Save outputs ──
raw_df.to_csv(DATA_DIR / "metabolites_raw.csv", index=False)
unique_df.to_csv(DATA_DIR / "metabolites_unique.csv", index=False)

# ── Summary ──
print(f"\nTotal species: {len(raw_df)}")
print(f"Unique metabolites (by base_id): {len(unique_df)}")
print(f"With embedded formula: {unique_df['embedded_formula'].notna().sum()}")
print(f"Without formula: {unique_df['embedded_formula'].isna().sum()}")
print(f"Dirty: {unique_df['is_dirty'].sum()}")
print(f"\nCompartment distribution:")
print(raw_df["compartment"].value_counts().to_string())

print(f"\nFiles saved:")
print(f"  data/metabolites_raw.csv    ({len(raw_df)} rows)")
print(f"  data/metabolites_unique.csv ({len(unique_df)} rows)")
print(f"  data/dirty_entries.csv      ({len(dirty_df)} rows)")
