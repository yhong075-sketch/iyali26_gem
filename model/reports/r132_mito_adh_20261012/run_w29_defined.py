"""Task B reconciliation in the W29-type defined medium (user definition 2026-10-13); diagnostic only.

Medium (in memory): Task A's medium (or its glycerol / oleate variants) with the PO1f profile, then
adenine + 11 amino acids closed, R612 (URA3) restored to (0, 1000), substrate uptake fixed at the
measured value; V2 = uracil closed (primary), V1 = uracil open at 1000. A third column repeats the
overlay variants on the PO1f + CSM-Leu medium as a sensitivity (it adds amino acids).

Step 5 variants: holistic overlay (baseline + overlay), holistic + shuttles restricted (candidate +
overlay + 13C caps, three glycine-cleavage forms), plus "as supplied" and "candidate + overlay, no caps".
Route FVA: the step-4 analysis (baseline vs candidate, WT and complex I KO, glucose) in V2 and V1.

Run from the repository root:
    .venv/bin/python -B model/reports/r132_mito_adh_20261012/run_w29_defined.py
Writes w29_defined_run.json here and pFBA vectors to $IYALI26_RESEARCH_ROOT/artifacts/tasks/<task>/runs_w29_defined/.
"""
import json
import math
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from holistic_overlay import apply_overlay  # noqa: E402
from run_task_b import (AMINO_ACID_ROWS, INPUTS, MEDIA, MODELS, OUT, REPO, Counter_, add_caps, context,  # noqa: E402
                        fva, growth, pfba_vector, sha, verify_cap_fractions)

OUT2 = OUT.parent / "runs_w29_defined"
OVERLAY_SOURCE_SHA = "f28cef1670fc07ad8e6a2a66616c4db693fe1635cd660ba5a308ae5f367a6392"
EXPECTED = {  # user-provided, predicted/measured growth
    "V2": {"holistic_overlay": (1.103, 1.279, 1.138), "holistic_shuttles": (0.916, 1.054, 0.983)},
    "V1": {"holistic_overlay": (1.116, 1.293, 1.151), "holistic_shuttles": (0.928, 1.068, 0.994)},
}
SUBSTRATES = ("glucose", "glycerol", "oleic_acid")
CARRIERS = {"R807": 1, "R831": -1, "R902": 1, "R1386": 1}
ROUTE_FVA = ("R132", "R1114", "R997", "R690", "R534", "R533", "R693", "R570", "R2062", "R691", "R136", "R1311")
STEP5_FVA = ("R570", "R2062", "R614", "R690")
VARIANTS = (  # name, model, overlay, gcv cap form
    ("as_supplied", "baseline", False, None),
    ("holistic_overlay", "baseline", True, None),
    ("candidate_overlay_no_caps", "candidate", True, None),
    ("holistic_shuttles_R361", "candidate", True, "R361"),
    ("holistic_shuttles_R361_looser", "candidate", True, "R361_looser"),
    ("holistic_shuttles_gcv_module", "candidate", True, "module"),
)


def setup(model_key, substrate, medium):
    """medium: 'V1', 'V2' (W29-type defined) or 'CSM' (PO1f + CSM-Leu sensitivity)."""
    sim, model = context(model_key, substrate)
    caps = json.loads(INPUTS["caps"][0].read_text())
    exchange, uptake = caps["substrates"][substrate]["exchange"], caps["substrates"][substrate]["measured_uptake"]
    if medium in ("V1", "V2"):
        current = dict(model.medium)
        for rid in AMINO_ACID_ROWS:
            assert rid in current, rid
            current.pop(rid)
        if medium == "V2":
            current.pop("R1354")
        else:
            current["R1354"] = 1000.0
        model.medium = current
        r612 = model.reactions.R612
        assert tuple(r612.bounds) == (0.0, 0.0), "PO1f URA3 overlay not found"
        r612.bounds = (0.0, 1000.0)
        model.reactions.get_by_id(exchange).bounds = (-uptake, -uptake)
    return sim, model, caps


def net_malate_range(model, mu):
    with model:
        model.reactions.biomass_C.lower_bound = mu * (1 - 1e-9)
        model.objective = model.problem.Objective(
            sum(c * model.reactions.get_by_id(r).flux_expression for r, c in CARRIERS.items()))
        out = []
        for direction in ("min", "max"):
            model.objective.direction = direction
            Counter_.lps += 1
            out.append(model.slim_optimize())
    return out


def main():
    for key, (path, digest) in {**MODELS, **MEDIA, **INPUTS}.items():
        if digest:
            assert sha(path) == digest, key
    source = REPO.parent / "competent-turing-bc6e9b/model/reports/w29_glucose_uptake_20261012/run_task_a.py"
    assert sha(source) == OVERLAY_SOURCE_SHA, "Task A overlay source changed"
    caps = json.loads(INPUTS["caps"][0].read_text())
    verify_cap_fractions(caps)
    OUT2.mkdir(parents=False, exist_ok=False)
    started = time.monotonic()
    record = {"started_utc": datetime.now(timezone.utc).isoformat(),
              "command": "python -B model/reports/r132_mito_adh_20261012/run_w29_defined.py",
              "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
              "overlay_source_sha256": OVERLAY_SOURCE_SHA, "overlay_module_sha256": sha(HERE / "holistic_overlay.py"),
              "medium_definition": "Task A medium (+ PO1f profile) -> close R1003 + 11 amino acids, restore R612 (0,1000), "
                                   "fix uptake at measured; V2 uracil closed, V1 uracil 1000; CSM = PO1f + CSM-Leu "
                                   "sensitivity (uptake <= measured, amino acids present)",
              "expected": EXPECTED, "step5": [], "route": {}}
    vectors = {}
    for medium in ("V2", "V1", "CSM"):
        for substrate in SUBSTRATES:
            measured = caps["substrates"][substrate]["measured_growth"]
            for name, model_key, overlay, gcv in VARIANTS:
                if medium == "CSM" and name == "as_supplied":
                    continue  # already in run.json
                sim, model, _ = setup(model_key, substrate, medium)
                with model:
                    if overlay:
                        apply_overlay(model)
                    applied = add_caps(model, caps, substrate, gcv) if gcv else {}
                    g, s = growth(model)
                    row = {"medium": medium, "substrate": substrate, "variant": name, "model": model_key,
                           "overlay": overlay, "caps": applied, "status": s, "growth": g,
                           "predicted_over_measured": None if g is None else g / measured}
                    if g and medium == "V2" and substrate == "glucose":
                        uptake = caps["substrates"][substrate]["measured_uptake"]
                        f = fva(model, STEP5_FVA)
                        row["fva_per_100"] = {k: [100 * v[0] / uptake, 100 * v[1] / uptake] for k, v in f.items()} \
                            if "error" not in f else f
                        if name in ("as_supplied", "holistic_overlay", "holistic_shuttles_R361",
                                    "holistic_shuttles_gcv_module"):
                            vectors[f"step5|{medium}|{substrate}|{name}"] = pfba_vector(model)
                record["step5"].append(row)
                print(medium, substrate, name, s, g, None if g is None else round(g / measured, 4), flush=True)
    for medium in ("V2", "V1"):
        for model_key in ("baseline", "candidate"):
            for label, ko in (("WT", None), ("complex_I_KO", "R2062")):
                sim, model, _ = setup(model_key, "glucose", medium)
                with model:
                    if ko:
                        model.reactions.get_by_id(ko).bounds = (0.0, 0.0)
                    g, s = growth(model)
                    entry = {"growth": g, "status": s}
                    if g:
                        f = fva(model, ROUTE_FVA)
                        entry["fva_per_100"] = {k: [100 * v[0] / 2.414, 100 * v[1] / 2.414] for k, v in f.items()} \
                            if "error" not in f else f
                        entry["net_malate_import_per_100"] = [100 * v / 2.414 for v in net_malate_range(model, g)]
                        if medium == "V2":
                            vectors[f"route|{medium}|{model_key}|{label}"] = pfba_vector(model)
                record["route"][f"{medium}|{model_key}|{label}"] = entry
                print("route", medium, model_key, label, s, g, flush=True)
    # Stop rule (V1 only): R361-looser is the reproduction of the earlier run.
    checks = []
    for variant, key in (("holistic_overlay", "holistic_overlay"), ("holistic_shuttles_R361_looser", "holistic_shuttles")):
        for i, substrate in enumerate(SUBSTRATES):
            row = next(r for r in record["step5"] if r["medium"] == "V1" and r["substrate"] == substrate and r["variant"] == variant)
            expected = EXPECTED["V1"][key][i]
            rel = None if row["predicted_over_measured"] is None else row["predicted_over_measured"] / expected - 1
            checks.append({"variant": variant, "substrate": substrate, "expected": expected,
                           "observed": row["predicted_over_measured"], "relative_difference": rel,
                           "within_2pct": rel is not None and abs(rel) <= 0.02})
    record["v1_stop_checks"] = checks
    record["v1_stop_fired"] = not all(c["within_2pct"] for c in checks)
    (OUT2 / "pfba_vectors.json").open("x").write(json.dumps(vectors) + "\n")
    record.update(lps=Counter_.lps, elapsed_seconds=time.monotonic() - started,
                  outputs={p.name: sha(p) for p in sorted(OUT2.iterdir())})
    (HERE / "w29_defined_run.json").open("x").write(json.dumps(record, indent=1, default=str) + "\n")
    print("V1 stop fired:", record["v1_stop_fired"], "| LPs", Counter_.lps, "| elapsed", round(record["elapsed_seconds"]), "s")


if __name__ == "__main__":
    main()
