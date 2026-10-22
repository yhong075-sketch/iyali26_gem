"""Why FVA reports R752's maximum below its value at the growth optimum (read-only, 6 LPs).

Usage (from platform/): PYTHONPATH=$PWD python ../model/reports/tpp_requirement_20261022/scripts/r752_max_check.py \
    NEW.xml MEDIUM.csv PROFILE.json
Maximises R752 directly and maximises growth, with GLPK and with Gurobi.
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
    growth = model.slim_optimize(error_value=float("nan"))
    with model:
        model.objective = "R752"
        r752_max = model.slim_optimize(error_value=float("nan"))
        biomass_at_r752_max = model.reactions.biomass_C.flux
    with model:
        model.reactions.biomass_C.lower_bound = growth * 0.999
        model.objective = "R752"
        model.objective_direction = "min"
        r752_min_near_opt = model.slim_optimize(error_value=float("nan"))
    print(f"{solver}: growth max {growth!r}; R752 max {r752_max!r} (biomass_C there {biomass_at_r752_max!r}, "
          f"1.2e-6 x that = {1.2e-6 * biomass_at_r752_max!r}); R752 min at 99.9% growth {r752_min_near_opt!r} "
          f"(1.2e-6 x 0.999 x growth = {1.2e-6 * 0.999 * growth!r})")
