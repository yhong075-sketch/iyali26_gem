"""Requested R305 chemistry, R385 growth dependency and closed ATP checks."""
import copy
from datetime import datetime, timezone
from fractions import Fraction
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))

from cobra import Reaction
from scripts.build_coq_biomass_candidate import pool_balance
from scripts.gem_annotate.energy_candidates import model_definition
from scripts.gem_annotate.execution import execution_limits
from scripts.validate_energy_candidates import SolverBudget, load_model, close_model, balance, sha, write, table


def exact_residual(model, coefficients):
    elements, charge = {}, Fraction(0)
    for mid, coefficient in coefficients.items():
        met = model.metabolites.get_by_id(mid)
        if met.elements is None or met.charge is None:
            raise ValueError('Missing chemical identity: '+mid)
        coefficient = Fraction(str(coefficient))
        for element, count in met.elements.items():
            elements[element] = elements.get(element,Fraction(0))+coefficient*Fraction(str(count))
        charge += coefficient*Fraction(str(met.charge))
    return {'element_residual':{e:str(v) for e,v in elements.items()},'charge_residual':str(charge),
            'balanced':not any(elements.values()) and charge==0}


def main():
    assert not (OUT/'manifest.json').exists(), 'No overwrite or automatic retry'
    build_path = ROOT/'artifacts/coq_biomass_candidate_20261005/E5_coq9_biomass_alpha_1e-4_validated.build.json'
    build = json.loads(build_path.read_text())
    variants = {'baseline':{'path':build['source'],'sha256':build['source_sha256']},
                'candidate':{'path':build['output'],'sha256':build['output_sha256']}}
    config_path = ROOT/'artifacts/atp_candidate_repair_20260924/config.json'
    config = json.loads(config_path.read_text())
    config['model'],config['model_sha256'] = variants['candidate']['path'],variants['candidate']['sha256']
    config['limits'] = {'solves':12,'solve_wall_seconds':720,'per_solve_seconds':60}
    mechanism_path = ROOT/'data/coq9_curation.json'
    mechanism = json.loads(mechanism_path.read_text())
    protected = {v['path']:v['sha256'] for v in variants.values()}
    for p in [config_path,build_path,mechanism_path,Path(__file__),ROOT/config['media'],ROOT/config['strain_profile']]:
        protected[str(p.resolve())] = sha(p)
    for module in tuple(sys.modules.values()):
        path = getattr(module,'__file__',None)
        if path and Path(path).resolve().is_relative_to(ROOT/'scripts'):
            protected[str(Path(path).resolve())] = sha(path)
    assert all(sha(p)==h for p,h in protected.items())
    manifest = {'status':'running','started_utc':datetime.now(timezone.utc).isoformat(),
        'scope':'Check existing files only; no R305 mechanism change or model export. Eight planned LPs, up to four anomaly-directed LPs reserved.',
        'models':variants,'config':config,'protected_sha256':protected,
        'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'git_dirty':subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True).splitlines(),
        'software':{'python':sys.version,**{p:version(p) for p in ('cobra','optlang','gurobipy')}},
        'historical_environment_reconstructed':False,'results':{},'positive_closed_ATP_cases':[]}
    write(OUT/'manifest.json',manifest)
    budget = SolverBudget(OUT/'solver_budget.json',config)
    growth, energy = [], []
    try:
        for variant, identity in variants.items():
            sim = load_model(identity['path'],config)
            model = sim.model
            before = model_definition(model)
            r305 = model.reactions.get_by_id('R305')
            coefficients = {m.id:c for m,c in r305.metabolites.items()}
            proposed = {**coefficients,**mechanism['qcycle_coefficients']}
            chemistry = {'reaction':'R305','name':r305.name,'equation':r305.reaction,
                'coefficients':coefficients,'bounds':list(r305.bounds),'notes':r305.notes,
                'species':{m.id:{'name':m.name,'formula':m.formula,'charge':m.charge,'compartment':m.compartment} for m in r305.metabolites},
                'actual':exact_residual(model,coefficients),
                'documented_qcycle_coefficients':mechanism['qcycle_coefficients'],
                'documented_qcycle_scope':mechanism['qcycle_scope'],
                'proposal_arithmetic_only':exact_residual(model,proposed),
                'matches_documented_qcycle':all(coefficients[k]==v for k,v in mechanism['qcycle_coefficients'].items()),
                'proposal_applied_in_this_run':False}
            pool = pool_balance(model)
            assert pool==({'R385':1.} if variant=='baseline' else {'R385':1.,'biomass_C':-1e-4})
            record = {'simulation_context':sim.provenance(),'medium':sim.active_medium,
                      'strain_overlay':sim.strain_overlay_audit,'R305_chemistry':chemistry,'combined_CoQ_pool_row':pool}
            manifest['results'][variant] = record
            for knockout in (False,True):
                with model:
                    if knockout:
                        model.reactions.get_by_id('R385').bounds = (0,0)
                    label = variant+('_R385_off' if knockout else '_WT')
                    result,flux = budget.solve(model,label,OUT/(label+'.json'),
                        {'kind':'maximum_growth','R385_disabled':knockout,'variant':variant})
                    mu = flux['biomass_C'] if flux is not None else None
                    row = {'variant':variant,'R385_disabled':knockout,'status':result['status'],'growth_h_inverse':mu,
                        'R385_flux':flux['R385'] if flux is not None else None,
                        'pool_residual':sum(c*flux[r] for r,c in pool.items()) if flux is not None else None,
                        'NGAM_bounds':list(model.reactions.get_by_id(config['maintenance']).bounds)}
                    if flux is not None:
                        assert abs(row['pool_residual']) <= 2*config['tolerance']
                    growth.append(row)
            assert model_definition(model)==before
            closed,changes = close_model(model,config)
            record['energy_closure'] = changes
            record['zero_flux_feasible_by_variable_and_constraint_bounds'] = True
            assert all(r.bounds==(0,0) for r in closed.reactions if not r.reactants or not r.products)
            for compartment in ('C_cy','C_mi'):
                with closed:
                    if compartment=='C_cy':
                        drain = closed.reactions.get_by_id(config['maintenance'])
                    else:
                        drain = Reaction('DIAG_ATP_DISSIPATION_mi',lower_bound=0,upper_bound=1000)
                        expected_species = {'m46[C_mi]':(-1,'C10H16N5O13P3'),
                            'm26[C_mi]':(-1,'H2O'),'m197[C_mi]':(1,'C10H15N5O10P2'),'m58[C_mi]':(1,'H3O4P')}
                        for mid,(_,formula) in expected_species.items():
                            met = closed.metabolites.get_by_id(mid)
                            assert met.formula==formula and met.charge==0
                        drain.add_metabolites({closed.metabolites.get_by_id(mid):c for mid,(c,_) in expected_species.items()})
                        closed.add_reactions([drain])
                    chem = balance(drain)
                    assert chem['element_status']=='balanced' and chem['charge_status']=='balanced_as_stored'
                    assert drain.lower_bound==0
                    closed.objective = drain
                    closed.objective.direction = 'max'
                    label = variant+'_ATP_closed_'+compartment
                    result,flux = budget.solve(closed,label,OUT/(label+'.json'),
                        {'kind':'maximum_balanced_ATP_dissipation','compartment':compartment,'variant':variant})
                    value = flux[drain.id] if flux is not None else None
                    verdict = ('unresolved' if result['status']!='optimal' or value is None
                        else 'positive_substrate_free_ATP' if value>config['tolerance'] else 'within_tolerance')
                    energy.append({'variant':variant,'compartment':compartment,'status':result['status'],
                        'max_ATP_dissipation':value,'verdict':verdict,'drain_reaction':drain.reaction})
                    if verdict=='positive_substrate_free_ATP':
                        manifest['positive_closed_ATP_cases'].append(label)
            assert model_definition(model)==before
            write(OUT/'manifest.json',manifest)
            table(OUT/'growth.tsv',growth); table(OUT/'energy.tsv',energy)
        assert len(budget.record['calls'])==8
        assert all(sha(p)==h for p,h in protected.items()), 'Input/code changed'
        manifest.update(status='complete',actual_optimization_calls=8,inputs_unchanged=True,
                        growth=growth,energy=energy,chemistry_pass=False)
    except BaseException as exc:
        manifest.update(status='failed',error=repr(exc))
        raise
    finally:
        write(OUT/'manifest.json',manifest)


if __name__=='__main__':
    with execution_limits(allow_network=False):
        main()
