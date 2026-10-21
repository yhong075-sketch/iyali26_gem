"""Task D: pFBA vectors for the maps (one optimum per panel; judgments rest on the FVA in the JSON outputs).

Panels: WT, complex I KO, KO + G5 (R533 rev closed), KO + G3 (R831 fwd closed) in A and V2; WT and KO under the
13C floors PDH + CS + G6PDH with the R825 ASSUMPTION (A and V2); KO after each discovery round (A).

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/complexI_transhydrogenation_20261012/map_vectors.py
Writes $IYALI26_RESEARCH_ROOT/artifacts/tasks/<task>/maps/map_vectors_pfba.json. Never overwrites.
"""
import json
from datetime import datetime, timezone

from common import ART, CONFIGS, HERE, Counter_, check_pins, close, growth, ko, load, mfa_rows, pfba_vector, r825_assumption
from discovery import merged_closures
from ratios import FLOORS, add_floors

CUTS = {"G5": {"R533": "rev"}, "G3": {"R831": "fwd"}}


def panel(config, state, closures=None, floors=False):
    m, _ = load(config)
    if state == "KO":
        ko(m)
    close(m, closures or {})
    if floors:
        mfa = mfa_rows()
        fractions = {n: float(mfa[row]["glucose_LB"]) / 100.0 for n, (_, row) in FLOORS.items()}
        r825_assumption(m)
        add_floors(m, ["PDH", "CS", "G6PDH"], fractions)
    value, status = growth(m)
    vector = pfba_vector(m) if value is not None else None
    return vector, {"growth": value, "status": status, "closures": closures or {}, "floors": floors}


def main():
    inputs = check_pins()
    definition = json.loads((HERE / "parts_definition.json").read_text())
    vectors, meta = {}, {}
    for config in CONFIGS:
        for key, args in ((f"{config}|WT", ("WT",)), (f"{config}|KO", ("KO",)),
                          (f"{config}|KO+G5", ("KO", CUTS["G5"])), (f"{config}|KO+G3", ("KO", CUTS["G3"])),
                          (f"{config}|WT|floors", ("WT", None, True)), (f"{config}|KO|floors", ("KO", None, True))):
            vectors[key], meta[key] = panel(config, *args)
            print(key, meta[key]["growth"], meta[key]["status"], flush=True)
    for r in (1, 2, 3, 4):
        closures = merged_closures([p for p in definition["parts"] if p["round_added"] < r])
        key = f"A|KO|round{r}"
        vectors[key], meta[key] = panel("A", "KO", closures)
        print(key, meta[key]["growth"], flush=True)
    out = ART / "maps"
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "map_vectors_pfba.json", "x") as handle:
        json.dump({"run_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs, "meta": meta,
                   "vectors": vectors, "lps": Counter_.lps}, handle)


if __name__ == "__main__":
    main()
