"""Task D: why the single-group cut sets {G3} and {G5} make complex I essential (diagnostic only).

Added after the tier-1 result (not pre-registered in TASK.md); every number here is a diagnostic of that result.
1. Matrix NADH: in configurations A and V2, FVA at the optimum of every reaction producing or consuming matrix
   NADH (m30), for WT, KO, and KO with G3 (R831 fwd) or G5 (R533 rev) closed.
2. Requirement of each matrix-NADH producer in WT: FVA minimum at the WT optimum.
3. Amino-acid rescue of the blocked KO: KO + G5 with lysine uptake opened (1000), then with every amino-acid
   exchange of the Task A medium opened (1000).
4. Robustness: KO/WT for nothing closed, G3 and G5, with the overlay and caps each on or off (A and V2).

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/complexI_transhydrogenation_20261012/mechanism.py
Writes mechanism.json here. Never overwrites.
"""
import json
import math
import time
from datetime import datetime, timezone

from common import (AMINO_ACID_ROWS, CONFIGS, HERE, NODES, PARAMS, PINS, UPTAKE, Counter_, apply_overlay, add_caps,
                    check_pins, close, fva, growth, ko, load, per100, write_new)
from cobra.util.solver import linear_reaction_coefficients  # noqa: E402
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context  # noqa: E402
from tools.diagnose_closed_energy import configure_solver  # noqa: E402

G3 = {"R831": "fwd"}   # tier 2: citrate_import_mito_831 alone resolves G3
G5 = {"R533": "rev"}   # matrix_MDH_reductive
SETS = {"none": {}, "G3": G3, "G5": G5}


def load_variant(config, overlay, caps):
    """As common.load, with the overlay and the caps each optional."""
    sim = load_effective_simulation_context(model_path=PINS["model"][0], media_path=PINS["medium"][0],
                                            strain_profile_path=PINS["profile"][0])
    model = sim.model
    model.solver = "gurobi"
    configure_solver(model, {"solver": PARAMS})
    assert {r.id: c for r, c in linear_reaction_coefficients(model).items()} == {"biomass_C": 1.0}
    assert math.isclose(-model.reactions.R1070.lower_bound, UPTAKE)
    if config == "V2":
        current = dict(model.medium)
        for rid in AMINO_ACID_ROWS:
            current.pop(rid)
        current.pop("R1354")
        model.medium = current
        model.reactions.R612.bounds = (0.0, 1000.0)
        model.reactions.R1070.bounds = (-UPTAKE, -UPTAKE)
    if overlay:
        apply_overlay(model)
    if caps:
        add_caps(model, json.loads(PINS["caps"][0].read_text()))
    return model


def lysine_exchange(model):
    hits = [r.id for r in model.exchanges if any((m.name or "").split("_")[0] == "L-lysine" for m in r.metabolites)]
    assert len(hits) == 1, hits
    return hits[0]


def main():
    inputs = check_pins()
    started = time.monotonic()
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs, "G3": G3, "G5": G5,
              "matrix_NADH": {}, "wt_requirement": {}, "rescue": {}, "robustness": []}
    for config in CONFIGS:
        wt, _ = load(config)
        met = wt.metabolites.get_by_id(NODES["NADH_mi"])
        rxns = sorted(r.id for r in met.reactions)
        coef = {r.id: r.get_coefficient(met) for r in met.reactions}
        w = growth(wt)[0]
        ranges = fva(wt, rxns, 1.0)
        record["wt_requirement"][config] = {
            rid: {"coef": coef[rid], "fva_per100": [round(per100(x), 4) for x in ranges[rid]],
                  "name": wt.reactions.get_by_id(rid).name,
                  "required_producer": (coef[rid] > 0 and ranges[rid][0] > 1e-6) or (coef[rid] < 0 and ranges[rid][1] < -1e-6)}
            for rid in rxns}
        for label, closures in SETS.items():
            kom, _ = load(config)
            ko(kom)
            close(kom, closures)
            k, status = growth(kom)
            entry = {"KO": k, "KO_status": status, "WT": w, "KO_over_WT": None if k is None else k / w}
            if k:
                r_ = fva(kom, rxns, 1.0)
                entry["matrix_NADH_fva_per100"] = {rid: [round(per100(x), 4) for x in r_[rid]] for rid in rxns
                                                   if max(abs(x) for x in r_[rid]) > 1e-6}
            record["matrix_NADH"][f"{config}|KO+{label}"] = entry
        # Amino-acid rescue of KO + G5.
        kom, _ = load(config)
        ko(kom)
        close(kom, G5)
        lys = lysine_exchange(kom)
        base = growth(kom)[0]
        with kom:
            kom.reactions.get_by_id(lys).lower_bound = -1000.0
            plus_lys = growth(kom)[0]
        with kom:
            for rid in AMINO_ACID_ROWS[1:]:  # R1003 is adenine
                kom.reactions.get_by_id(rid).lower_bound = -1000.0
            plus_all = growth(kom)[0]
        record["rescue"][config] = {"lysine_exchange": lys, "KO_G5": base, "KO_G5_lysine_open": plus_lys,
                                    "KO_G5_all_amino_acids_open": plus_all, "WT": w,
                                    "amino_acid_rows_opened": list(AMINO_ACID_ROWS[1:])}
        print(config, "rescue", record["rescue"][config], flush=True)
    for config in CONFIGS:
        for overlay in (True, False):
            for caps in (True, False):
                wt = load_variant(config, overlay, caps)
                kom = load_variant(config, overlay, caps)
                ko(kom)
                w = growth(wt)[0]
                row = {"config": config, "overlay": overlay, "caps": caps, "WT": w}
                for label, closures in SETS.items():
                    with wt:
                        close(wt, closures)
                        ws = growth(wt)[0]
                    with kom:
                        close(kom, closures)
                        k, status = growth(kom)
                    row[label] = {"WT": ws, "WT_cost": None if ws is None else 1 - ws / w, "KO": k, "KO_status": status,
                                  "KO_over_WT": None if not ws or k is None else k / ws}
                record["robustness"].append(row)
                print(config, overlay, caps, {k_: (None if row[k_]["KO_over_WT"] is None else round(row[k_]["KO_over_WT"], 4))
                                              for k_ in SETS}, flush=True)
    record.update(lps=Counter_.lps, elapsed_seconds=time.monotonic() - started)
    write_new(HERE / "mechanism.json", json.dumps(record, indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()
