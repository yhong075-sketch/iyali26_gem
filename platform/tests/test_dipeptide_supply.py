import unittest

from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import pfba

from tools.diagnose_dipeptide_supply import (
    MODEL, CONDITIONS, TARGETS, ENDO_REQUIRED, balance, formula_elements,
    scenario, signature, validate_endogenous_spec,
    load_effective_simulation_context, cached_solution,
)


class DipeptideSupplyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sim = load_effective_simulation_context(model_path=MODEL,
            media_path=CONDITIONS/'media/sd_leu.csv',
            strain_profile_path=CONDITIONS/'strain_profiles/po1f_sd_leu.json')

    def test_all_diagnostics_restore_actual_model_even_on_exception(self):
        model = self.sim.model
        before = signature(model)
        with self.assertRaisesRegex(RuntimeError, 'deliberate'):
            with scenario(model,tuple(TARGETS),0.01,water=True,drains=True,pump=True,outputs=True):
                self.assertEqual(len([r for r in model.reactions if r.id.startswith('DIAG_SUP_')]),4)
                self.assertEqual(model.reactions.R1363.upper_bound,0.04)
                self.assertEqual(model.reactions.R795.upper_bound,0.04)
                model.objective=model.reactions.R2021
                model.add_cons_vars(model.problem.Constraint(model.reactions.R2021.flux_expression,lb=0.005,name='test_joint'))
                raise RuntimeError('deliberate')
        self.assertEqual(signature(model),before)

    def test_invalid_supply_does_not_leak(self):
        model=self.sim.model
        before=signature(model)
        for value in (-1,float('nan'),float('inf'),1000):
            with self.assertRaises(ValueError):
                with scenario(model,('R2021',),value):
                    pass
            self.assertEqual(signature(model),before)

    def test_missing_formula_is_not_mass_balance_pass(self):
        self.assertIsNone(formula_elements(None))
        self.assertIsNone(formula_elements('C8H11NO10PR2'))
        self.assertEqual(formula_elements('C2H5NO2'),{'C':2,'H':5,'N':1,'O':2})
        for rid in TARGETS:
            self.assertEqual(balance(self.sim.model.reactions.get_by_id(rid))['element_status'],'unverifiable')

    def test_endogenous_module_disabled_and_missing_parameters_rejected(self):
        spec={'enabled':False,**dict.fromkeys(ENDO_REQUIRED)}
        before=signature(self.sim.model)
        self.assertEqual(validate_endogenous_spec(spec)['added_reactions'],[])
        self.assertEqual(signature(self.sim.model),before)
        spec['enabled']=True
        with self.assertRaisesRegex(ValueError,'precursor_identity_and_sequence'):
            validate_endogenous_spec(spec)
        with self.assertRaisesRegex(ValueError,'explicitly'):
            validate_endogenous_spec({'enabled':'false'})

    def test_saved_solution_rejects_changed_growth_or_missing_witness(self):
        row={'objective':1,'status':'optimal','numerically_valid':True,'growth_floor':0.99}
        witness={'fluxes':{'r':1},'extra':{'growth_floor':0.99}}
        _,sol=cached_solution(row,witness,{'growth_floor':0.99})
        self.assertEqual(sol.fluxes['r'],1)
        with self.assertRaisesRegex(ValueError,'constraints'):
            cached_solution(row,witness,{'growth_floor':0.5})
        with self.assertRaisesRegex(ValueError,'witness'):
            cached_solution(row,None,{'growth_floor':0.99})

    def test_pfba_zero_does_not_mean_blocked(self):
        model=Model('alternative_path')
        a,b,c=[Metabolite(x,compartment='c') for x in 'abc']
        for rid,stoich,ub in [('source',{a:1},1),('short',{a:-1,c:1},1),
                              ('long1',{a:-1,b:1},1),('long2',{b:-1,c:1},1),('growth',{c:-1},1)]:
            r=Reaction(rid,lower_bound=0,upper_bound=ub);r.add_metabolites(stoich);model.add_reactions([r])
        model.objective=model.reactions.growth
        solution=pfba(model)
        self.assertAlmostEqual(solution.fluxes['long1'],0)
        model.reactions.growth.lower_bound=1
        model.objective=model.reactions.long1
        self.assertAlmostEqual(model.slim_optimize(),1)


if __name__=='__main__':
    unittest.main()
