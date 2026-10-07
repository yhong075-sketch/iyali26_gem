"""Audit recorded real optimization results and rollback; do not rerun optimization."""

import csv
import hashlib
import json
from pathlib import Path
import unittest

from cobra.io import read_sbml_model

from scripts.gem_annotate.execution import execution_limits
from tools.validate_vacuole_supply import HYDRO, dipeptides, reaction_diff, signature, supply_case

from scripts.gem_annotate.config import resolve_recorded_path
ROOT = Path(__file__).resolve().parents[2]
from scripts.gem_annotate.config import load_project_paths
# The saved run (about 230 MB of flux results) lives in the research workspace, not in git.
TASK = load_project_paths().task_outputs / "vacuole_open_supply_20260924"
RUN = TASK / "run"


def read_json(path):
    return json.loads(path.read_text())


def rows(name):
    with (RUN / (name + ".tsv")).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@unittest.skipUnless((RUN / "manifest.json").exists(), "saved run not found; set IYALI26_RESEARCH_ROOT")
class VacuoleSupplyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.guard = execution_limits(no_solve=True, allow_network=False)
        cls.attempts = cls.guard.__enter__()
        cls.config = read_json(TASK / "config.json")
        cls.manifest = read_json(RUN / "manifest.json")
        cls.tol = cls.config["tolerance"]
        if cls.manifest["status"] != "complete":
            raise AssertionError("Actual optimization run is not complete")

    @classmethod
    def tearDownClass(cls):
        cls.guard.__exit__(None, None, None)
        assert cls.attempts == {"optimization": 0, "network": 0}, cls.attempts

    def test_actual_run_identity_budget_and_restoration(self):
        self.assertEqual(self.manifest["config"], self.config)
        self.assertEqual(self.manifest["config_sha256"], hashlib.sha256((TASK / "config.json").read_bytes()).hexdigest())
        for path, digest in {**self.config["input_configuration_sha256"], **self.manifest["source_sha256"]}.items():
            self.assertEqual(hashlib.sha256(resolve_recorded_path(path).read_bytes()).hexdigest(), digest, path)
        self.assertEqual(hashlib.sha256((ROOT / self.config["candidate"]).read_bytes()).hexdigest(), self.manifest["candidate_sha256"])
        self.assertTrue(self.manifest["all_temporary_changes_removed"])
        self.assertTrue(self.manifest["energy_locks_match_E5"])
        budget = read_json(TASK / "budget.json")
        calls = budget["calls"]
        self.assertEqual(len(calls), self.manifest["optimization_calls"])
        self.assertEqual(len(calls), len(self.manifest["outcomes"]))
        self.assertEqual(len(calls), len({call["label"] for call in calls}))
        self.assertLessEqual(len(calls), self.config["limits"]["solves"])
        self.assertLessEqual(budget["solve_call_seconds"], self.config["limits"]["solve_wall_seconds"])
        self.assertAlmostEqual(sum(call["solve_call_seconds"] for call in calls), budget["solve_call_seconds"], delta=1e-9)
        for call in calls:
            self.assertEqual(call["postprocessing"], "complete", call["label"])
            self.assertIn(call["status"], ("optimal", "infeasible"))
            self.assertEqual(call["solver"], self.config["solver"])
            self.assertLessEqual(call["solve_call_seconds"], self.config["limits"]["per_solve_seconds"])
            if call["status"] != "optimal":
                self.assertEqual(call["purpose"]["kind"], "joint_hydrolysis")

    def test_twelve_closed_energy_maxima_and_three_zero_controls(self):
        energy = rows("closed_energy_results")
        self.assertEqual(len(energy), 15)
        self.assertEqual({(r["scenario"], r["carrier"]) for r in energy},
                         {(c, e) for c in ("C0", "C1", "C2") for e in ("ATP", "GTP", "UTP", "CTP", "zero_feasibility")})
        for row in energy:
            self.assertEqual(row["status"], "optimal")
            self.assertLessEqual(abs(float(row["maximum"])), self.tol)
            label = row["scenario"] + "_" + ("zero" if row["carrier"] == "zero_feasibility" else row["carrier"])
            raw = read_json(RUN / (label + ".json"))
            self.assertEqual(raw["result"]["status"], row["status"])
            self.assertEqual(raw["result"]["objective"], float(row["maximum"]))
            self.assertFalse(any(r.startswith("DIAG_SUP_") for r in raw["problem"]["reactions"]))
            if row["carrier"] == "zero_feasibility":
                self.assertLessEqual(max(abs(v) for v in raw["fluxes"].values()), self.tol)
            else:
                self.assertEqual(raw["problem"]["reactions"]["xMAINTENANCE"][:2], [0, 1000])

    def test_culture_growth_source_caps_and_pfba_objective(self):
        growth = rows("growth_and_hydrolysis")
        primary = {r["scenario"]: float(r["biomass_C"]) for r in growth if r["objective_kind"] == "growth"}
        self.assertTrue({"G" + str(i) for i in range(12)} <= primary.keys())
        for row in growth:
            label = row["scenario"] + "_" + row["objective_kind"]
            raw = read_json(RUN / (label + ".json"))
            flux, problem = raw["fluxes"], raw["problem"]
            self.assertEqual(row["status"], "optimal")
            self.assertEqual(problem["reactions"]["xMAINTENANCE"][:2], [7.8625, 1000])
            self.assertGreaterEqual(float(row["xMAINTENANCE"]), 7.8625 - self.tol)
            self.assertEqual(float(row["biomass_C"]), flux["biomass_C"])
            self.assertEqual(float(row["objective_value"]), raw["result"]["objective"])
            supplied = json.loads(row["supplied"])
            self.assertEqual({rid for rid in problem["reactions"] if rid.startswith("DIAG_SUP_")},
                             {"DIAG_SUP_" + rid for rid in supplied})
            for rid in HYDRO:
                value = float(row["source_" + rid])
                self.assertEqual(value, flux.get("DIAG_SUP_" + rid, 0))
                self.assertGreaterEqual(value, -self.tol)
                self.assertLessEqual(value, float(row["epsilon"]) + self.tol)
                if rid in supplied:
                    self.assertEqual(problem["reactions"]["DIAG_SUP_" + rid][:2], [0, float(row["epsilon"])])
            if row["objective_kind"] == "growth":
                self.assertAlmostEqual(float(row["objective_value"]), float(row["biomass_C"]), delta=self.tol)
            elif row["objective_kind"] == "pfba":
                self.assertGreaterEqual(float(row["biomass_C"]), primary[row["scenario"]] - self.config["pfba_growth_slack"] - self.tol)
                self.assertEqual(problem["direction"], "min")
                self.assertTrue(all(value == 1 for value in raw["actual_solver"]["objective_coefficients"]))
                self.assertAlmostEqual(sum(raw["actual_solver"]["primal"]), float(row["objective_value"]), delta=self.tol)
        pairs = [("G1", "G0"), ("G2", "G0"), ("G3", "G1"), ("G3", "G2")]
        pairs += [("G" + str(i), "G1") for i in range(4, 8)] + [("G3", "G" + str(i)) for i in range(4, 12)]
        for relaxed, restricted in pairs:
            self.assertGreaterEqual(primary[relaxed], primary[restricted] - self.tol, (relaxed, restricted))
        for row in rows("monotonicity"):
            self.assertAlmostEqual(float(row["difference"]), primary[row["relaxed_scenario"]] - primary[row["restricted_scenario"]], delta=self.tol)
            self.assertEqual(row["passes"], "True")
        previous = primary["G1"]
        for epsilon in self.config["sensitivity_epsilon"]:
            value = primary["G3_eps_" + str(epsilon)]
            self.assertGreaterEqual(value, previous - self.tol)
            self.assertLessEqual(value, primary["G3"] + self.tol)
            previous = value

    def test_joint_status_and_two_distinct_fva_growth_floors(self):
        mu = self.manifest["primary_optima"]["G3"]
        floors = {"none": None, "99pct": self.config["fva_fraction"] * mu,
                  "near_strict": mu - self.config["strict_growth_delta"]}
        joint = rows("joint_feasibility")
        self.assertEqual({r["growth_requirement"] for r in joint}, set(floors))
        for row in joint:
            raw = read_json(RUN / (row["case"] + ".json"))
            self.assertEqual(raw["result"]["status"], row["status"])
            self.assertIn(row["status"], ("optimal", "infeasible"))
            floor = floors[row["growth_requirement"]]
            self.assertEqual(float(row["growth_floor"]) if row["growth_floor"] else None, floor)
            if row["status"] == "optimal":
                self.assertEqual(float(row["biomass_C"]), raw["fluxes"]["biomass_C"])
                hydrolysis = json.loads(row["hydrolysis_fluxes"])
                self.assertEqual(hydrolysis, {rid: raw["fluxes"][rid] for rid in HYDRO})
                for value in hydrolysis.values():
                    self.assertGreaterEqual(value, float(row["hydrolysis_minimum"]) - self.tol)
                if floor is not None:
                    self.assertGreaterEqual(float(row["biomass_C"]), floor - self.tol)
            else:
                self.assertIsNone(raw["fluxes"])
                self.assertEqual(row["biomass_C"], "")
        fva = rows("targeted_fva")
        self.assertEqual(len(fva), 20)
        self.assertEqual({(r["requirement"], r["reaction"], r["direction"]) for r in fva},
                         {(f, rid, d) for f in ("99pct", "near_strict") for rid in (*HYDRO, "R795") for d in ("min", "max")})
        values = {}
        for row in fva:
            raw = read_json(RUN / ("G3_FVA_" + row["requirement"] + "_" + row["reaction"] + "_" + row["direction"] + ".json"))
            self.assertEqual(row["status"], "optimal")
            self.assertEqual(row["status"], raw["result"]["status"])
            self.assertEqual(float(row["flux"]), raw["fluxes"][row["reaction"]])
            self.assertEqual(float(row["biomass_C"]), raw["fluxes"]["biomass_C"])
            self.assertEqual(float(row["growth_floor"]), floors[row["requirement"]])
            self.assertEqual(float(row["growth_floor"]), raw["result"]["purpose"]["growth_floor"])
            self.assertGreaterEqual(float(row["biomass_C"]), float(row["growth_floor"]) - self.tol)
            values[row["requirement"], row["reaction"], row["direction"]] = float(row["flux"])
        for requirement in ("99pct", "near_strict"):
            for rid in (*HYDRO, "R795"):
                self.assertLessEqual(values[requirement, rid, "min"], values[requirement, rid, "max"] + self.tol)
                self.assertGreaterEqual(values["near_strict", rid, "min"], values["99pct", rid, "min"] - self.tol)
                self.assertLessEqual(values["near_strict", rid, "max"], values["99pct", rid, "max"] + self.tol)

    def test_material_ledger_and_nominal_source_carbon_nitrogen(self):
        ledger = rows("proton_water_atp_ledger")
        totals = {}
        for row in ledger:
            contribution = float(row["coefficient"]) * float(row["flux"])
            self.assertAlmostEqual(contribution, float(row["contribution"]), delta=self.tol)
            key = row["witness"], row["metabolite"]
            totals[key] = totals.get(key, 0) + contribution
        for row in rows("ledger_totals"):
            self.assertAlmostEqual(float(row["residual"]), totals[row["witness"], row["metabolite"]], delta=self.tol)
            self.assertLessEqual(abs(float(row["residual"])), self.tol)
            self.assertAlmostEqual(float(row["production"]) - float(row["consumption"]), float(row["residual"]), delta=self.tol)
        pools = read_json(RUN / "dipeptide_identities.json")
        for row in rows("source_carbon_nitrogen"):
            pool = pools[row["source"].removeprefix("DIAG_SUP_")]
            self.assertAlmostEqual(float(row["nominal_C_input"]), float(row["source_flux"]) * pool["nominal_C"], delta=self.tol)
            self.assertAlmostEqual(float(row["nominal_N_input"]), float(row["source_flux"]) * pool["nominal_N"], delta=self.tol)
            self.assertIn("nominal", row["evidence"])
            if not pool["formula"]:
                self.assertEqual(row["chemistry_status"], "unverifiable_missing_formula")

    def test_xml_only_four_bounds_and_real_context_rollback_without_solve(self):
        base = read_sbml_model(ROOT / self.config["model"])
        model = read_sbml_model(ROOT / self.config["candidate"])
        actual = reaction_diff(base, model)
        expected = [{"reaction": rid, "field": "bounds", "before": [0., 0.], "after": bounds}
                    for rid, bounds in self.config["connection_bounds"].items()]
        self.assertEqual(sorted(actual, key=lambda r: r["reaction"]), sorted(expected, key=lambda r: r["reaction"]))
        pools = dipeptides(model)
        before = signature(model)
        for raise_exception in (False, True):
            try:
                with supply_case(model, pools, HYDRO, self.config["main_epsilon"], close="R795"):
                    model.add_cons_vars(model.problem.Constraint(model.reactions.biomass_C.flux_expression, lb=0.1, name="DIAG_test_rollback"))
                    self.assertTrue(all("DIAG_SUP_" + rid in model.reactions for rid in HYDRO))
                    self.assertIn("DIAG_test_rollback", model.constraints)
                    if raise_exception:
                        raise RuntimeError("deliberate rollback test")
            except RuntimeError as exc:
                self.assertEqual(str(exc), "deliberate rollback test")
            self.assertEqual(signature(model), before)


if __name__ == "__main__":
    unittest.main()
