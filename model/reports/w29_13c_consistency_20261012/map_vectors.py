"""Task C: pFBA vectors for map_c2_acetyl_routes.png (one optimum per panel; diagnostic).

V2 + overlay + PDH >= 0.866629 x uptake; caps {none, module} x {WT, complex I KO} x R825 {as supplied,
reversible}. Run with the project venv (Gurobi 13.0.1) from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/w29_13c_consistency_20261012/map_vectors.py
Writes $IYALI26_RESEARCH_ROOT/artifacts/tasks/<task>/maps/map_vectors_pfba.json. Never overwrites.
"""
import json

from common import ART, add_caps, apply_overlay, check_pins, growth, load, load_caps, pdh_floor, pfba_vector, r825_reversible, sha

OUT = ART / "maps"


def main():
    inputs = check_pins()
    caps = load_caps()
    OUT.mkdir(parents=True, exist_ok=False)
    sim, model = load("V2")
    vectors, meta = {}, {}
    for capped in (False, True):
        for ko in (False, True):
            for rev in (False, True):
                key = f"caps={'module' if capped else 'none'}|{'KO' if ko else 'WT'}|{'rev' if rev else 'asis'}"
                with model:
                    apply_overlay(model)
                    if capped:
                        add_caps(model, caps, "module")
                    pdh_floor(model)
                    if rev:
                        r825_reversible(model)
                    if ko:
                        model.reactions.R2062.bounds = (0.0, 0.0)
                    g, s = growth(model)
                    vectors[key] = pfba_vector(model) if g else None
                    meta[key] = {"growth": g, "status": s}
    (OUT / "map_vectors_pfba.json").open("x").write(json.dumps({"inputs": inputs, "vectors": vectors, "meta": meta}) + "\n")
    print(meta, sha(OUT / "map_vectors_pfba.json"))


if __name__ == "__main__":
    main()
