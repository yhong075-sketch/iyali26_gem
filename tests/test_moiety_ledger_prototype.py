"""Tests for the isolated aggregate acyl-moiety accounting proof of concept."""

import unittest

from scripts.moiety_ledger_prototype import build_pa_moiety_ledger_demo, internal_mass_balance_errors, solve_demo


class MoietyLedgerPrototypeTests(unittest.TestCase):
    def test_internal_ledger_reactions_are_element_and_charge_balanced(self) -> None:
        model = build_pa_moiety_ledger_demo()
        self.assertEqual(internal_mass_balance_errors(model), {})

    def test_aggregate_pa_output_is_feasible(self) -> None:
        model = build_pa_moiety_ledger_demo()
        solution = model.optimize()
        self.assertEqual(solution.status, "optimal")
        self.assertAlmostEqual(solution.objective_value, 1.0)

    def test_both_positional_routes_remain_feasible_and_are_indistinguishable(self) -> None:
        results = solve_demo()
        self.assertAlmostEqual(results["without_lpa_C16_0"], 1.0)
        self.assertAlmostEqual(results["without_lpa_C18_1"], 1.0)

    def test_a_required_acyl_supply_cannot_be_replaced_by_the_other_chain(self) -> None:
        results = solve_demo()
        self.assertAlmostEqual(results["without_acyl_c16_0_supply"], 0.0)


if __name__ == "__main__":
    unittest.main()
