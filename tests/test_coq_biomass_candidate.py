"""Software-only checks on a synthetic pool; no E5 coefficient is selected here."""

import copy
import unittest

from cobra import Model, Metabolite, Reaction

from scripts.build_coq_biomass_candidate import apply_coq_biomass, pool_balance, Q9, Q9H2
from scripts.gem_annotate.energy_candidates import model_definition
from scripts.gem_annotate.execution import execution_limits


class CoqBiomassTests(unittest.TestCase):
    def test_pool_coupling_is_scoped_idempotent_and_rejects_invalid_inputs(self):
        with execution_limits(no_solve=True, allow_network=False) as attempts:
            model = Model("synthetic_software_test")
            q = Metabolite(Q9, formula="C54H82O4", charge=0, compartment="C_mi")
            qh = Metabolite(Q9H2, formula="C54H84O4", charge=0, compartment="C_mi")
            carbon = Metabolite("test_carbon", formula="C", compartment="C_cy")
            source, redox, bio = (Reaction(r) for r in ("R385", "redox", "biomass_C"))
            source.add_metabolites({q: 1})
            redox.add_metabolites({q: -1, qh: 1})
            bio.add_metabolites({carbon: -2})
            model.add_reactions([source, redox, bio])
            model.objective = bio
            before = model_definition(model)
            for alpha in (None, True, 0, -1, float("nan"), float("inf")):
                with self.subTest(alpha=alpha), self.assertRaises(ValueError):
                    apply_coq_biomass(model, alpha, "software fixture")
                self.assertEqual(model_definition(model), before)
            with self.assertRaises(ValueError):
                apply_coq_biomass(model, 0.125, " ")
            alpha = 0.125  # Synthetic test number, not a physiological estimate.
            original_notes = copy.deepcopy(bio.notes)
            with model:
                apply_coq_biomass(model, alpha, "temporary synthetic fixture")
                self.assertEqual(bio.metabolites[q], -alpha)
            bio.notes = original_notes  # COBRA contexts do not track notes.
            self.assertEqual(model_definition(model), before)
            apply_coq_biomass(model, alpha, "synthetic software fixture")
            expected = copy.deepcopy(before)
            expected["reactions"]["biomass_C"]["stoichiometry"][Q9] = -alpha
            self.assertEqual(model_definition(model), expected)
            self.assertEqual(pool_balance(model), {"R385": 1.0, "biomass_C": -alpha})
            apply_coq_biomass(model, alpha, "synthetic software fixture")
            self.assertEqual(model_definition(model), expected)
            with self.assertRaises(ValueError):
                apply_coq_biomass(model, 2 * alpha, "conflicting coefficient")
            self.assertEqual(model_definition(model), expected)
            leak = Reaction("unexpected_Q9_sink")
            leak.add_metabolites({qh: -1})
            model.add_reactions([leak])
            conflict = model_definition(model)
            with self.assertRaises(ValueError):
                apply_coq_biomass(model, alpha, "software fixture")
            self.assertEqual(model_definition(model), conflict)
            self.assertEqual(attempts, {"optimization": 0, "network": 0})


if __name__ == "__main__":
    unittest.main()
