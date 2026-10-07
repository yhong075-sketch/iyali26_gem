"""Actual optimization controls, charged to the task's declared solver budget."""
import json
import os
from pathlib import Path
import tempfile
import unittest

from scripts.validate_energy_candidates import (
    TASK, SolverBudget, energy_verdict, load_model, json_safe,
)
from scripts.diagnose_closed_energy import close_model
from scripts.diagnose_dipeptide_supply import signature


class EnergyCandidateBehaviorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads((TASK/'config.json').read_text())
        requested = os.environ.get('IYALI26_ENERGY_TEST_OUTPUT')
        cls.out = Path(requested) if requested else Path(tempfile.mkdtemp(prefix='energy_behavior_', dir=TASK))
        if requested:
            cls.out.mkdir(parents=True, exist_ok=False)
        cls.budget = SolverBudget(os.environ.get('IYALI26_ENERGY_BUDGET', str(TASK/'budget.json')), cls.config)
        cls.sim = load_model(TASK/'candidates/E3.xml', cls.config)

    def solve(self, model, label):
        return self.budget.solve(model, self.out.name+'/'+label, self.out/(label+'.json'),
                                 {'kind':'behavior_test','assertion':label})

    def test_known_three_synthesis_directions_stay_closed_after_xml_reload(self):
        model,_ = close_model(self.sim.model, self.config)
        for rid,sign in [('R_PGAM1_PhosHydro',1),('R_NTP3pp',1),('R_NTP7',-1)]:
            with model:
                r=model.reactions.get_by_id(rid)
                model.objective=model.problem.Objective(sign*r.flux_expression,direction='max')
                result,_=self.solve(model,'known_direction_'+rid)
                self.assertEqual(result['status'],'optimal')
                self.assertLessEqual(abs(result['objective']),self.config['tolerance'])

    def test_remaining_full_model_anomaly_is_not_false_pass(self):
        model,_=close_model(self.sim.model,self.config)
        result,_=self.solve(model,'E3_remaining_ATP')
        self.assertEqual(result['status'],'optimal')
        self.assertGreater(result['objective'],self.config['tolerance'])
        self.assertEqual(energy_verdict(result['status'],result['objective'],self.config['tolerance']),
                         'positive_energy_regeneration')

    def test_zero_feasible_and_real_culture_restored(self):
        original=self.sim.model; before=signature(original)
        model,_=close_model(original,self.config)
        for r in model.reactions: r.bounds=(0,0)
        model.objective=model.problem.Objective(0)
        result,flux=self.solve(model,'closed_exact_zero')
        self.assertEqual(result['status'],'optimal')
        self.assertTrue(all(abs(v)<=self.config['tolerance'] for v in flux.values()))
        self.assertEqual(signature(original),before)
        self.assertEqual(original.reactions.xMAINTENANCE.bounds,(7.8625,1000))
        self.assertNotIn('R1219',self.sim.active_medium)
        self.assertIn('biomass_C',str(original.objective.expression))
        result,flux=self.solve(original,'restored_PO1f_SD_Leu_growth')
        self.assertEqual(result['status'],'optimal')
        self.assertGreater(result['objective'],0)
        self.assertGreaterEqual(flux['xMAINTENANCE'],7.8625-self.config['tolerance'])
        self.assertEqual(signature(original),before)

    def test_unsolved_or_nonfinite_energy_never_passes(self):
        for status,value in [('infeasible',0),('time_limit',0),('software_error',None),
                             ('optimal',float('nan')),('optimal',float('inf')),('optimal',-1)]:
            self.assertEqual(energy_verdict(status,value,1e-7),'unresolved')
        self.assertEqual(energy_verdict('optimal',0,1e-7),'within_tolerance')
        self.assertEqual(json_safe([float('inf'),float('-inf'),float('nan')]),
                         ['Infinity','-Infinity','NaN'])

    def test_final_export_closed_ATP_and_normal_growth_positive_control(self):
        sim=load_model(TASK/'candidates/E5.xml',self.config)
        before=signature(sim.model)
        model,_=close_model(sim.model,self.config)
        result,_=self.solve(model,'E5_closed_ATP_reproduction')
        self.assertEqual(result['status'],'optimal')
        self.assertLessEqual(abs(result['objective']),self.config['tolerance'])
        self.assertEqual(sim.model.reactions.xMAINTENANCE.lower_bound,7.8625)
        result,flux=self.solve(sim.model,'E5_real_growth_positive_control')
        self.assertEqual(result['status'],'optimal')
        self.assertGreater(result['objective'],0)
        self.assertGreaterEqual(flux['xMAINTENANCE'],7.8625-self.config['tolerance'])
        self.assertEqual(signature(sim.model),before)


if __name__=='__main__':
    unittest.main()
