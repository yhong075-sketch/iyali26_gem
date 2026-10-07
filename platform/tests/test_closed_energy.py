"""Closure and scenario isolation checks without optimization."""
import json
import unittest

from tools.diagnose_closed_energy import DEFAULT, block_direction, close_model, load_baseline, configure_solver
from tools.diagnose_dipeptide_supply import balance, signature


class ClosedEnergyTests(unittest.TestCase):
    def test_closure_direction_and_isolation(self):
        config = json.loads((DEFAULT/'config.json').read_text())
        original = load_baseline(config).model
        before = signature(original)
        closed, changes = close_model(original, config)
        self.assertEqual({k:getattr(closed.solver.problem.Params,k) for k in config['solver']},config['solver'])
        self.assertEqual(sum(row['after']==[0,0] for row in changes), 189)
        self.assertEqual(closed.reactions.newBiom.bounds, (0,0))
        self.assertEqual(closed.reactions.xMAINTENANCE.bounds, (0,1000))
        self.assertEqual(balance(closed.reactions.xMAINTENANCE)['element_status'],'balanced')
        for r in original.reactions:
            if r.id not in {c['reaction'] for c in changes}:
                self.assertEqual(closed.reactions.get_by_id(r.id).bounds,r.bounds)
        altered=closed.copy()
        self.assertEqual(configure_solver(altered,config),config['solver'])
        block_direction(altered,'R603','reverse')
        self.assertEqual(altered.reactions.R603.lower_bound,0)
        self.assertEqual(altered.reactions.R603.upper_bound,closed.reactions.R603.upper_bound)
        self.assertLess(closed.reactions.R603.lower_bound,0)
        self.assertEqual(signature(original),before)

    def test_unreviewed_constraint_rejected(self):
        config=json.loads((DEFAULT/'config.json').read_text())
        model=load_baseline(config).model
        model.add_cons_vars(model.problem.Constraint(model.reactions.R603.flux_expression,
                                                     lb=1,name='unexpected_force'))
        with self.assertRaisesRegex(ValueError,'custom constraints'):
            close_model(model,config)


if __name__=='__main__':
    unittest.main()
