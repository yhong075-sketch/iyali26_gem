"""Extra targeted checks for REPORT.md (investigation only, in-memory, same inputs and solver as classify_fn.py).

1. Full GPR text of R171, R794, R2062, R570 and the bypass reactions R2024, R1252, R1172, R2081/R334.
2. biomass_C: which cofactor-type metabolites it drains (folate, NAD, FAD, FMN, CoA, inositol, sphingolipid, GPI).
3. Reaction-level KO + FVA(0.99) bypass trace for ENO (YALI1F22436g) and PYK (YALI1F12842g).
4. Vitamins closed: KO of BNA4 (YALI1D12493g), INO1 (YALI1B05798g), FOL1 (YALI1A03069g), FMN1 (YALI1B02806g),
   QNS1 (YALI1A21156g), and which exchanges or reactions then supply NAD, inositol, folate, FMN.
5. R2024 reversed (forward-only, lb 0) combined with the ADK1 (YALI1B01136g) and IPP1 (YALI1F18105g) KOs.
6. R1172 closed combined with the ERG13 (YALI1F38116g) and ERG10 (YALI1B11532g) KOs; R1252 closed combined with the GAPDH (YALI1C08262g) KO.

Usage (repository root, IYALI26_RESEARCH_ROOT set):
    python -B model/reports/essentiality_fn_causes_20261027/trace_extra.py
"""
import hashlib, json, math, os, sys
from pathlib import Path
REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "platform"))
from cobra.flux_analysis import flux_variability_analysis
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from tools.diagnose_closed_energy import configure_solver

HERE = Path(__file__).resolve().parent
PARAMS = dict(Threads=1, TimeLimit=60.0, Presolve=0, FeasibilityTol=1e-9, OptimalityTol=1e-9, IntFeasTol=1e-9)
TOL = 1e-6
VITAMINS = ["R1029", "R1340", "R1305", "R2270", "R946", "R1263", "R971", "R1127", "R1268", "R1314"]
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    research = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve()
    task = research / "artifacts/tasks" / HERE.name
    man = json.loads((HERE / "screen_a1237821.manifest.json").read_text())
    inputs = {k: Path(v["path"]) for k, v in man["inputs"].items()}
    assert {k: sha(p) for k, p in inputs.items()} == {k: v["sha256"] for k, v in man["inputs"].items()}
    sim = load_effective_simulation_context(model_path=inputs["model"], media_path=inputs["medium"],
                                            strain_profile_path=inputs["profile"])
    m = sim.model
    m.solver = "gurobi"
    configure_solver(m, {"solver": PARAMS})
    wt = man["WT"]["raw_growth"]
    out, lps = {}, [0]

    def grow():
        lps[0] += 1
        s = m.optimize()
        return float(s.objective_value) if s.status == "optimal" else None, s.status

    def rx(rid):
        r = m.reactions.get_by_id(rid)
        return dict(name=r.name, reaction=r.reaction, bounds=list(r.bounds), rule=r.gene_reaction_rule,
                    notes=str(r.notes)[:400])
    out["reactions"] = {rid: rx(rid) for rid in ["R171", "R794", "R2062", "R570", "R2024", "R1252", "R1172",
                                                 "R2081", "R334", "R2022", "R2023", "R1158", "R1968"] if rid in m.reactions}
    # same-compartment duplicates of R570 by stoichiometry
    bm = m.reactions.biomass_C
    out["biomass"] = {met.id: dict(name=met.name, coef=c) for met, c in bm.metabolites.items()}

    def trace(gene_rxns, label, extra=None):
        with m:
            if extra:
                extra()
            g0, _ = grow()
            for rid in gene_rxns:
                m.reactions.get_by_id(rid).bounds = (0, 0)
            g, st = grow()
            res = dict(reference_growth=g0, ko_growth=g, status=st, ratio=(g / g0 if g and g0 else None))
            if g and g > TOL:
                f = flux_variability_analysis(m, fraction_of_optimum=0.99, processes=1)
                lps[0] += 1 + 2 * len(m.reactions)
                req = {r: [float(f.loc[r, "minimum"]), float(f.loc[r, "maximum"])] for r in f.index
                       if f.loc[r, "minimum"] * f.loc[r, "maximum"] > 0 and min(abs(f.loc[r, "minimum"]), abs(f.loc[r, "maximum"])) >= TOL}
                res["required"] = req
            out.setdefault("traces", {})[label] = res
            return res

    # WT required set for comparison
    f0 = flux_variability_analysis(m, fraction_of_optimum=0.99, processes=1)
    lps[0] += 1 + 2 * len(m.reactions)
    wt_req = {r for r in f0.index if f0.loc[r, "minimum"] * f0.loc[r, "maximum"] > 0
              and min(abs(f0.loc[r, "minimum"]), abs(f0.loc[r, "maximum"])) >= TOL}
    def closev():
        for rid in VITAMINS:
            m.reactions.get_by_id(rid).lower_bound = 0
    with m:
        closev()
        fv = flux_variability_analysis(m, fraction_of_optimum=0.99, processes=1)
        lps[0] += 1 + 2 * len(m.reactions)
    vit_req = {r for r in fv.index if fv.loc[r, "minimum"] * fv.loc[r, "maximum"] > 0
               and min(abs(fv.loc[r, "minimum"]), abs(fv.loc[r, "maximum"])) >= TOL}
    genes = lambda gid: [r.id for r in m.genes.get_by_id(gid).reactions]
    for gid, lab in [("YALI1F22436g", "ENO"), ("YALI1F12842g", "PYK")]:
        r = trace(genes(gid), lab)
        r["newly_required"] = {k: v for k, v in r.get("required", {}).items() if k not in wt_req}
    for gid, lab in [("YALI1D12493g", "BNA4_vit_closed"), ("YALI1B05798g", "INO1_vit_closed"),
                     ("YALI1A03069g", "FOL1_vit_closed"), ("YALI1B02806g", "FMN1_vit_closed"),
                     ("YALI1A21156g", "QNS1_vit_closed")]:
        r = trace(genes(gid), lab, closev)
        r["newly_required_vs_vit_wt"] = {k: v for k, v in r.get("required", {}).items() if k not in vit_req}
    def r2024_fwd():
        m.reactions.R2024.lower_bound = 0
    out["r2024_wt"] = trace([], "WT_R2024_forward_only", r2024_fwd)
    trace(genes("YALI1B01136g"), "ADK1_R2024_forward_only", r2024_fwd)
    trace(genes("YALI1F18105g"), "IPP1_R2024_forward_only", r2024_fwd)
    trace(genes("YALI1F38116g"), "ERG13_R1172_closed", lambda: setattr(m.reactions.R1172, "bounds", (0, 0)))
    trace(genes("YALI1B11532g"), "ERG10_R1172_closed", lambda: setattr(m.reactions.R1172, "bounds", (0, 0)))
    trace(genes("YALI1C08262g"), "GAPDH_R1252_closed", lambda: setattr(m.reactions.R1252, "bounds", (0, 0)))
    for k, v in out["traces"].items():
        v.pop("required", None)  # keep the file small; deltas kept above
    g1, _ = grow()
    assert abs(g1 - wt) < 1e-9
    out["lps"] = lps[0]
    (task / "fn_classification/trace_extra.json").write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk in ("reference_growth", "ko_growth", "status", "ratio")}
                      for k, v in out["traces"].items()}, indent=1))
    print("LPs", lps[0])


if __name__ == "__main__":
    main()
