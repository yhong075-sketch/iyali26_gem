"""Six actual Euler dFBA trajectories using the existing kernel unchanged."""
import copy
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
ALPHAS = [0., 1e-4, 10**-3.5, 1e-3, 10**-2.5, 1e-2]
SETTINGS = dict(hours=24., step_hours=.1, initial_biomass_gdw_l=.05)
LIMITS = dict(backend_calls=3000, wall_seconds=1200, per_LP_seconds=60, threads=1, retry=False)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(name, obj):
    (OUT/name).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


def table(name, rows):
    with (OUT/name).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def run():
    from importlib.metadata import version
    import numpy as np
    from cobra.util.array import create_stoichiometric_matrix
    from scripts import dfba_new_fn_essentiality as kernel
    from scripts.build_coq_biomass_candidate import apply_coq_biomass, pool_balance, Q9, BIOMASS, SPEC_PATH
    from scripts.gem_annotate.energy_candidates import model_definition
    from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
    from scripts.gem_annotate.execution import execution_limits

    assert not (OUT/'manifest.json').exists(), 'No overwrite or retry'
    prior_path = ROOT/'artifacts/confusion_matrix_E5_20261001/run_manifest.json'
    prior = json.loads(prior_path.read_text())
    historical = json.loads((ROOT/'artifacts/dfba_E5_iYali21_20261005/manifest.json').read_text())
    dynamic_path = ROOT/'data/media/po1f_csm_leu_dfba.csv'
    spec = json.loads(SPEC_PATH.read_text())
    assert prior['model']['sha256'] == spec['source_sha256']
    protected = {prior[k]['path']:prior[k]['sha256'] for k in ('model','medium','strain_profile')}
    protected[str(dynamic_path)] = historical['protected_sha256'][str(dynamic_path)]
    # Only these historical kernel/context identities affect the reused algorithm.
    for name in ('scripts/dfba_new_fn_essentiality.py', 'scripts/gem_annotate/essentiality_simulation_context.py', 'scripts/gem_annotate/strain_overlay.py'):
        p = str(ROOT/name)
        assert sha(p) == historical['protected_sha256'][p]
    protected.update({str(prior_path):sha(prior_path), str(SPEC_PATH):sha(SPEC_PATH), str(Path(__file__)):sha(__file__)})
    for module in tuple(sys.modules.values()):
        p = getattr(module, '__file__', None)
        if p and Path(p).resolve().is_relative_to(ROOT/'scripts'):
            protected[str(Path(p).resolve())] = sha(p)
    assert all(sha(p)==h for p,h in protected.items())
    medium = kernel.load_dynamic_medium(dynamic_path)
    started = time.monotonic()
    record = dict(status='running', started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        settings=SETTINGS, limits=LIMITS, alphas_mmol_per_gDW=ALPHAS, protected_sha256=protected,
        model=prior['model'], dynamic_medium=medium, solver_parameters=prior['actual_gurobi_parameters'],
        software={k:version(k) for k in ('cobra','optlang','gurobipy','numpy')}, python=sys.version,
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_dirty=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True).splitlines(),
        historical_environment_reconstructed=False, backend_calls=0, pfba_calls=0, trajectories=[],
        scope='Actual 24 h WT time courses; fixed alpha per trajectory, not extracellular CoQ concentration. No KO, FVA, energy retest or XML export.',
        alpha_evidence='User-selected sensitivity assumptions, not measured physiological values',
        medium_policy='Prior PO1f overlay followed by prior dynamic medium; finite uracil replaces static nonlimiting supply')
    save('manifest.json',record)
    all_points, summaries = [], []
    try:
        with execution_limits(allow_network=False):
            for index, alpha in enumerate(ALPHAS):
                sim = load_effective_simulation_context(model_path=prior['model']['path'],
                    media_path=prior['medium']['path'], strain_profile_path=prior['strain_profile']['path'])
                model = sim.model
                before = model_definition(model)
                if alpha:
                    apply_coq_biomass(model,alpha,'User-confirmed 2026-10-05 dFBA alpha sensitivity; uncalibrated')
                expected = copy.deepcopy(before)
                if alpha:
                    expected['reactions'][BIOMASS]['stoichiometry'][Q9] = -alpha
                assert model_definition(model) == expected
                assert pool_balance(model) == ({'R385':1., BIOMASS:-alpha} if alpha else {'R385':1.})
                assert model.reactions.get_by_id('R385').bounds == (0.,1000.)
                model.solver = 'gurobi'
                model.solver.configuration.presolve = False
                for key,value in prior['actual_gurobi_parameters'].items():
                    setattr(model.solver.problem.Params,key,value)
                assert {k:getattr(model.solver.problem.Params,k) for k in prior['actual_gurobi_parameters']} == prior['actual_gurobi_parameters']
                rids = [r.id for r in model.reactions]
                stoich = create_stoichiometric_matrix(model,array_type='lil').tocsr()
                details = kernel._exchange_details(model,medium)
                state = dict(time=0., biomass=SETTINGS['initial_biomass_gdw_l'],
                    concentrations={r['reaction_id']:r['initial_concentration_mmol_l'] for r in medium if r['initial_concentration_mmol_l'] is not None})
                result = dict(alpha_mmol_per_gDW=alpha, status='running', backend_calls=0, pfba_calls=0,
                    runtime_context=sim.provenance(), overlay_audit=sim.strain_overlay_audit, reaction_order=rids,
                    max_mass_residual=0., max_bound_violation=0., max_pool_residual=0., tiny_negative_clips=0,
                    initial_R385_bounds=list(model.reactions.get_by_id('R385').bounds),
                    NGAM_bounds=list(model.reactions.get_by_id('xMAINTENANCE').bounds))
                label = f'alpha_{index:02d}'
                backend, ordinary_pfba = type(model.solver)._optimize, kernel.pfba
                points = [dict(alpha_mmol_per_gDW=alpha,time_hours=0.,biomass_gDW_L=state['biomass'])]
                growths = []

                def counted_backend(solver):
                    assert time.monotonic()-started < LIMITS['wall_seconds'], 'Wall budget reached'
                    assert record['backend_calls'] < LIMITS['backend_calls'], 'LP budget reached'
                    record['backend_calls'] += 1
                    result['backend_calls'] += 1
                    return backend(solver)

                with gzip.open(OUT/f'{label}_trace.jsonl.gz','xt') as trace:
                    def observed_pfba(m):
                        record['pfba_calls'] += 1
                        result['pfba_calls'] += 1
                        sol = ordinary_pfba(m)
                        flux = sol.fluxes.reindex(rids).to_numpy(dtype=float)
                        assert sol.status == 'optimal' and np.isfinite(flux).all()
                        mass = float(np.max(np.abs(stoich@flux)))
                        violation = max(max(r.lower_bound-float(v),float(v)-r.upper_bound,0.) for r,v in zip(m.reactions,flux))
                        mu = float(sol.fluxes[BIOMASS])
                        pool_residual = float(sol.fluxes['R385'])-alpha*mu
                        assert mass <= 1e-7 and violation <= 1e-7 and abs(pool_residual) <= 2e-7
                        assert math.isfinite(mu) and mu >= -1e-9
                        step = min(SETTINGS['step_hours'],SETTINGS['hours']-state['time'])
                        raw_end = {d['reaction_id']:state['concentrations'][d['reaction_id']]-d['coefficient']*float(sol.fluxes[d['reaction_id']])*state['biomass']*step for d in details if d['reaction_id'] in state['concentrations']}
                        assert all(math.isfinite(v) and v >= -1e-8 for v in raw_end.values())
                        row = dict(time_hours=state['time'],step_hours=step,biomass_start=state['biomass'],
                            concentrations_start=dict(state['concentrations']),concentrations_raw_end=raw_end,
                            raw_growth=mu,R385_flux=float(sol.fluxes['R385']),pool_residual=pool_residual,
                            fluxes=flux.tolist(),exchange_bounds={d['reaction_id']:list(d['reaction'].bounds) for d in details},
                            max_mass_residual=mass,max_bound_violation=violation)
                        trace.write(json.dumps(row,allow_nan=False)+'\n')
                        result['max_mass_residual'] = max(result['max_mass_residual'],mass)
                        result['max_bound_violation'] = max(result['max_bound_violation'],violation)
                        result['max_pool_residual'] = max(result['max_pool_residual'],abs(pool_residual))
                        result['tiny_negative_clips'] += int(mu<0)+sum(v<0 for v in raw_end.values())
                        state['biomass'] += max(0.,mu)*state['biomass']*step
                        state['concentrations'] = {k:max(0.,v) for k,v in raw_end.items()}
                        state['time'] += step
                        growths.append(mu)
                        points.append(dict(alpha_mmol_per_gDW=alpha,time_hours=state['time'],biomass_gDW_L=state['biomass']))
                        return sol

                    try:
                        with patch.object(type(model.solver),'_optimize',counted_backend), patch.object(kernel,'pfba',observed_pfba):
                            terminal = kernel.simulate_dfba(model,medium,**SETTINGS,
                                model_identity={'role':'coq_alpha_time_course',**prior['model'],'alpha':alpha})
                        assert terminal['final_biomass_gdw_l'] == state['biomass']
                        assert terminal['final_concentrations_mmol_l'] == state['concentrations']
                        assert terminal['steps'] == 240 and len(points) == 241
                        assert model_definition(model) == expected
                        result.update(status='complete',terminal=terminal,min_growth=min(growths),max_growth=max(growths))
                    except BaseException as exc:
                        result.update(status='failed',error=repr(exc),diagnostic=getattr(exc,'diagnostic',None))
                        raise
                    finally:
                        result['observed_state'] = state
                        save(f'{label}_result.json',result)
                all_points.extend(points)
                summary = dict(alpha_mmol_per_gDW=alpha,status=result['status'],final_biomass_gDW_L=state['biomass'],
                    min_growth_h_inverse=min(growths),max_growth_h_inverse=max(growths),backend_calls=result['backend_calls'])
                summaries.append(summary)
                record['trajectories'].append(summary)
                table('growth_curves.tsv',all_points)
                table('endpoints.tsv',summaries)
                save('manifest.json',record)
                print(json.dumps(summary),flush=True)
            assert record['backend_calls']==2880 and record['pfba_calls']==1440
            assert all(sha(p)==h for p,h in protected.items()), 'Input/code changed during run'
            record.update(status='complete',inputs_unchanged=True)
    except BaseException as exc:
        record.update(status='failed',error=repr(exc))
        raise
    finally:
        record['elapsed_seconds'] = time.monotonic()-started
        save('manifest.json',record)


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    record = json.loads((OUT/'manifest.json').read_text())
    assert record['status']=='complete'
    with (OUT/'growth_curves.tsv').open() as f:
        rows = [{k:float(v) for k,v in row.items()} for row in csv.DictReader(f,delimiter='\t')]
    fig,ax = plt.subplots(figsize=(9.5,6),layout='constrained')
    colors = ['#303843','#3b8dbd','#45a390','#deac3c','#cb754a','#925893']
    markers = ['o','s','^','D','v','P']
    controls = [r['biomass_gDW_L'] for r in rows if r['alpha_mmol_per_gDW']==0]
    max_difference = 0.
    for i,alpha in enumerate(ALPHAS):
        values = [r for r in rows if r['alpha_mmol_per_gDW']==alpha]
        x = [r['time_hours'] for r in values]
        y = [r['biomass_gDW_L'] for r in values]
        max_difference = max(max_difference,float(np.max(np.abs(np.array(y)-controls))))
        label = '0 (control)' if alpha==0 else f'{alpha:.3g}'
        ax.plot(x,y,label=label,color=colors[i],lw=1.5,marker=markers[i],ms=5,
                markevery=(12+i*7,48),alpha=.9)
    ax.set(xlabel='Time (h)',ylabel='Biomass (gDW/L)',xlim=(0,24),ylim=(0,None),
           title='dFBA growth curves at different CoQ9 biomass contents')
    ax.grid(alpha=.2)
    ax.spines[['top','right']].set_visible(False)
    ax.legend(title=r'$\alpha$ (mmol/gDW)',frameon=False,loc='upper left')
    if max_difference < 1e-8:
        ax.text(.98,.1,'All six curves overlap\n'+r'Max. $|\Delta B|$'+f' = {max_difference:.2g} gDW/L',
                ha='right',transform=ax.transAxes,fontsize=10)
    fig.supxlabel('E5 / PO1f; finite SD-Leu medium; B₀ = 0.05 gDW/L; Euler Δt = 0.1 h\nAlpha values are sensitivity assumptions, not calibrated physiological contents.',fontsize=9)
    fig.savefig(OUT/'dfba_growth_curves.png',dpi=200)
    fig.savefig(OUT/'dfba_growth_curves.pdf')
    plt.close(fig)
    save('plot_record.json',dict(matplotlib=matplotlib.__version__,python=sys.version,code_sha256=sha(__file__),
        data_sha256=sha(OUT/'growth_curves.tsv'),max_absolute_biomass_difference=max_difference))


if __name__=='__main__':
    if sys.argv[1:]==['--plot-only']:
        plot()
    elif sys.argv[1:]==['--worker']:
        run()
    else:
        assert sys.argv[1:]==[] and not (OUT/'execution.json').exists(), 'No overwrite or retry'
        result = dict(status='running',limits=LIMITS)
        save('execution.json',result)
        try:
            with (OUT/'run.log').open('x') as f:
                completed = subprocess.run([sys.executable,'-B',__file__,'--worker'],cwd=ROOT,
                    stdout=f,stderr=subprocess.STDOUT,timeout=LIMITS['wall_seconds'])
            result.update(status='complete' if completed.returncode==0 else 'failed',returncode=completed.returncode)
        except subprocess.TimeoutExpired:
            result.update(status='timed_out')
        finally:
            save('execution.json',result)
        print(json.dumps(result))
