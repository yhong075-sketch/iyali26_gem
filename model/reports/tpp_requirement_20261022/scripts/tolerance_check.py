"""THI80 knockout at tight feasibility tolerance, with the TPP row residual (read-only, 4 LPs).

Usage (from platform/): PYTHONPATH=$PWD python ../model/reports/tpp_requirement_20261022/scripts/tolerance_check.py \
    NEW.xml MEDIUM.csv PROFILE.json
The TPP term (1.2e-6) is about 12x the default feasibility tolerance (1e-7), so a solver could in
principle satisfy the TPP row within tolerance and report KO growth up to about 0.08. This re-solves the
WT and the THI80 knockout at tolerance 1e-9 with GLPK and Gurobi and reports the largest mass-balance
residual on the TPP row (m975[C_cy]) and over all rows.
"""
import sys
import warnings

from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

model_path, media_path, profile_path = sys.argv[1:4]
for solver in ("glpk", "gurobi"):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = load_effective_simulation_context(model_path=model_path, media_path=media_path,
                                                  strain_profile_path=profile_path).model
    model.solver = solver
    model.tolerance = 1e-9
    for label, knockout in (("wt", False), ("thi80_ko", True)):
        with model:
            if knockout:
                model.genes.get_by_id("YALI1E25443g").knock_out()
            solution = model.optimize()
            fluxes = solution.fluxes
            residual = {m.id: sum(r.get_coefficient(m) * fluxes[r.id] for r in m.reactions)
                        for m in model.metabolites}
            worst = max(residual.items(), key=lambda item: abs(item[1]))
            print(f"{solver} tol 1e-9 {label}: {solution.status} growth {solution.objective_value!r}; "
                  f"TPP row residual {residual['m975[C_cy]']!r}; largest residual {worst[0]} {worst[1]!r}")
