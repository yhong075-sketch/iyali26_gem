"""Acceptance checks for the read-only strict sn-position lipid-core candidate."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import tempfile
import unittest
from math import fsum, isclose
from pathlib import Path

import cobra
from cobra.flux_analysis import flux_variability_analysis, pfba


REPOSITORY = Path(__file__).resolve().parents[1]
MODEL_PATH = Path(os.environ.get(
    "IYALI26_SOURCE_MODEL", REPOSITORY.parent / "iyali26_gem" / "model.xml"
))
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from scripts.lp_sn12_candidate import (  # noqa: E402
    CURATION_PATH, ContractError, MARKER, build_candidate, report, source_fingerprint,
)


class StrictSnCoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.curation = json.loads(CURATION_PATH.read_text())
        cls.source = cobra.io.read_sbml_model(str(MODEL_PATH))
        cls.source_fingerprint = source_fingerprint(cls.source)
        cls.source_objective = cls.source.slim_optimize(error_value=None)
        cls.candidate = build_candidate(cls.source)
        cls.generated = [reaction for reaction in cls.candidate.reactions if reaction.annotation.get(MARKER) == "true"]

    def test_source_contract_is_frozen_and_input_unchanged(self) -> None:
        self.assertEqual(hashlib.sha256(MODEL_PATH.read_bytes()).hexdigest(), self.curation["source"]["model_sha256"])
        self.assertEqual(self.source_fingerprint, self.curation["source"]["model_fingerprint"])
        self.assertEqual(self.source.slim_optimize(error_value=None), self.source_objective)
        for mutate in (lambda model: setattr(model.metabolites[0], "name", "drift"), lambda model: setattr(model.reactions[0], "annotation", {"drift": "true"})):
            with self.subTest(mutate=mutate):
                drifted = self.source.copy()
                mutate(drifted)
                with self.assertRaises(ContractError):
                    build_candidate(drifted)
        with self.assertRaises(ContractError):
            build_candidate(self.candidate)

    def test_structure_balance_and_template_copying(self) -> None:
        self.assertEqual(len(self.generated), 1134)
        self.assertEqual(source_fingerprint(build_candidate(self.source)), source_fingerprint(self.candidate))
        self.assertEqual(sum(reaction.id.startswith("UL_R353_") for reaction in self.generated), 7)
        self.assertEqual(sum(reaction.id.startswith("UL_R1771_") for reaction in self.generated), 343)
        self.assertEqual(len({met.id for reaction in self.generated for met in reaction.metabolites if met.annotation.get(MARKER) == "true"}), 1134)
        self.assertFalse(any(metabolite_id in self.candidate.metabolites for metabolite_id in self.curation["generic_acyl_coa_ids"]))
        self.assertFalse(any(reaction.id.startswith("xREMIX_") for reaction in self.candidate.reactions))
        for reaction in self.generated:
            template = self.source.reactions.get_by_id(reaction.annotation["template_reaction"])
            self.assertEqual(reaction.bounds, template.bounds)
            self.assertEqual(reaction.gene_reaction_rule, template.gene_reaction_rule)
            self.assertTrue(all(reaction.annotation.get(key) == value for key, value in template.annotation.items()))
            self.assertEqual(reaction.check_mass_balance(), {})

    def test_tuple_order_and_biomass_weights(self) -> None:
        self.assertIn("UL_R1846_SN1__lauroyl__SN2__myristoyl", self.candidate.reactions)
        self.assertIn("UL_R1846_SN1__myristoyl__SN2__lauroyl", self.candidate.reactions)
        tag = self.candidate.reactions.get_by_id("UL_R1771_SN1__lauroyl__SN2__myristoyl__SN3__oleoyl")
        self.assertIn("ul_dag_lp__sn1__lauroyl__sn2__myristoyl[C_lp]", [met.id for met in tag.metabolites])
        self.assertIn("m1490[C_lp]", [met.id for met in tag.metabolites])
        biomass = self.candidate.reactions.get_by_id("biomass_C")
        for row in self.curation["biomass"]["terms"]:
            terms = [(-coefficient, metabolite) for metabolite, coefficient in biomass.metabolites.items() if metabolite.annotation.get("state") == row["state"]]
            self.assertEqual(len(terms), 7 ** row["arity"])
            self.assertTrue(isclose(fsum(value for value, _ in terms), row["coefficient"], abs_tol=1e-18))

    def test_sbml_round_trip_and_function(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "candidate.xml"
            cobra.io.write_sbml_model(self.candidate, str(path))
            reread = cobra.io.read_sbml_model(str(path))
        self.assertEqual(len([reaction for reaction in reread.reactions if reaction.annotation.get(MARKER) == "true"]), 1134)
        tag = reread.reactions.get_by_id("UL_R1771_SN1__lauroyl__SN2__myristoyl__SN3__oleoyl")
        self.assertEqual(tag.gene_reaction_rule, self.source.reactions.get_by_id("R1771").gene_reaction_rule)
        self.assertEqual(tag.bounds, self.source.reactions.get_by_id("R1771").bounds)
        self.assertGreater(self.candidate.slim_optimize(error_value=None), 0)
        self.assertGreater(pfba(self.candidate).fluxes["biomass_C"], 0)
        payload = report(self.source, self.candidate)
        self.assertFalse(payload["ready_for_activation"])
        runtime = payload["runtime_seconds"]
        self.assertEqual(
            payload["performance_review_required"],
            any(
                runtime[f"candidate_{kind}_median"] > 2 * runtime[f"source_{kind}_median"]
                for kind in ("fba", "pfba")
            ),
        )
        self.assertGreaterEqual(payload["fba_probe"]["growth_ratio"], 0.99)
        self.assertLessEqual(payload["fba_probe"]["growth_ratio"], 1.0)

    def test_required_lipid_routes_are_not_bypassed(self) -> None:
        route_sets = {
            "pa_transport": ("UL_R1790_",), "pi": ("UL_R1839_",),
            "ps_pe_pc": ("UL_R1841_", "UL_R629_", "UL_R1835_", "UL_R1845_", "UL_R1844_"), "tag": ("UL_R1771_",),
        }
        for name, prefixes in route_sets.items():
            with self.candidate:
                for reaction in self.candidate.reactions:
                    if reaction.id.startswith(prefixes):
                        reaction.bounds = (0.0, 0.0)
                self.assertEqual(self.candidate.slim_optimize(error_value=None), 0.0, name)
        closed = self.curation["close_reactions"] + ["newBiom", "xBIOMASS", "xLIPID"]
        variability = flux_variability_analysis(self.candidate, reaction_list=closed, processes=1)
        self.assertTrue((variability["minimum"] == 0).all() and (variability["maximum"] == 0).all())

    def test_each_required_chain_is_not_bypassed(self) -> None:
        for chain in self.curation["chains"]:
            with self.candidate:
                for reaction in self.candidate.reactions:
                    if chain["lp_acyl_coa_id"] in (met.id for met in reaction.metabolites):
                        reaction.bounds = (0.0, 0.0)
                self.assertEqual(self.candidate.slim_optimize(error_value=None), 0.0, chain["id"])


if __name__ == "__main__":
    unittest.main()
