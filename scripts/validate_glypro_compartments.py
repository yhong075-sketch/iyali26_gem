"""Explicit-input, bounded Gly-Pro origin × hydrolysis compartment experiment."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from cobra import Reaction
from cobra.io import read_sbml_model

from scripts.diagnose_closed_energy import close_model, configure_solver, workspace
from scripts.diagnose_dipeptide_supply import sha, signature, table
from scripts.gem_annotate.energy_candidates import model_definition
from scripts.glypro_compartment_candidates import verify_candidate, SOURCE_SHA256, exact_copy
from scripts.validate_energy_candidates import add_ntp_dissipation, energy_verdict, write
from scripts.validate_vacuole_supply import Run

CY = 'HYP_GLYPRO_HYD_CY'
HYD = {'V_ONLY': 'R2039', 'C_ONLY': CY}
ORIGIN = {'O_CY': 'm1865[C_cy]', 'O_VA': 'm1866[C_va]'}
SOURCE = 'DIAG_GLYPRO_SOURCE'
REPORT = (SOURCE, CY, 'R2039', 'R2038', 'R2037', 'R2030', 'R2040',
          'R1363', 'R795', 'R871', 'R876', 'R2021', 'R2029', 'R2034', 'xMAINTENANCE', 'R171')
MIDS = ('m1864[C_ex]', 'm1865[C_cy]', 'm1866[C_va]', 'm272[C_cy]', 'm1863[C_va]',
        'm765[C_cy]', 'm1867[C_va]', 'm32[C_cy]', 'm1384[C_va]', 'm10[C_cy]',
        'm1007[C_va]', 'm141[C_cy]', 'm143[C_cy]', 'm35[C_cy]')


def clone(model, config):
    result = exact_copy(model)
    configure_solver(result, config)
    return result


def objective(model, rid, direction='max'):
    model.objective = model.reactions.get_by_id(rid)
    model.objective.direction = direction


def growth_floor(model, value):
    model.add_cons_vars(model.problem.Constraint(model.reactions.biomass_C.flux_expression,
                                                lb=value, name='DIAG_COMMON_GROWTH'))


def l1(model):
    model.objective = model.problem.Objective(
        sum(r.forward_variable + r.reverse_variable for r in model.reactions), direction='min')


def supplied(base, config, origin=None, fixed_q=None, hydro=None):
    model = clone(base, config)
    if origin:
        assert SOURCE not in model.reactions
        r = Reaction(SOURCE, name='Diagnostic Gly-L-Pro arrival; not endogenous production',
                     lower_bound=0, upper_bound=config['main_epsilon'])
        r.add_metabolites({model.metabolites.get_by_id(ORIGIN[origin]): 1})
        model.add_reactions([r])
        if fixed_q is not None:
            r.bounds = (fixed_q, fixed_q)
            model.reactions.get_by_id(hydro).bounds = (fixed_q, fixed_q)
    return model


class Experiment(Run):
    def record(self, table_name, label, origin, scenario, model, row, flux, **extra):
        record = {'label': label, 'origin': origin, 'scenario': scenario,
                  'status': row['status'], 'objective_value': row['objective'],
                  'biomass': flux['biomass_C'] if flux else None,
                  'selected_hydrolysis': flux.get(HYD.get(scenario, 'R2039')) if flux else None,
                  'source': flux.get(SOURCE, 0) if flux else None,
                  'R795': flux['R795'] if flux else None,
                  'related_fluxes': {rid: flux.get(rid, 0) for rid in REPORT} if flux else None,
                  'boundary_fluxes': {r.id: flux[r.id] for r in model.boundary} if flux else None,
                  'extra': extra}
        self.tables[table_name].append(record)
        self.save()
        return record

    def witness(self, model, flux, label, origin, scenario, q, growth):
        assert flux is not None and q > 100*self.config['tolerance']
        self.witnesses[label] = {'origin': origin, 'scenario': scenario, 'q': q,
            'growth_floor': growth, 'fluxes': flux,
            'full_problem_file': label+'.json', 'net_not_split': True}
        columns = {r.id: {m.id: c for m, c in r.metabolites.items()} for r in model.reactions}
        for mid in MIDS:
            met = model.metabolites.get_by_id(mid)
            contributions = []
            for r in sorted(met.reactions, key=lambda r: r.id):
                c = columns[r.id][mid]
                value = c*flux[r.id]
                contributions.append(value)
                self.tables['proton_water_atp_ledger'].append({
                    'witness': label, 'metabolite': mid, 'name': met.name,
                    'reaction': r.id, 'coefficient': c, 'signed_flux': flux[r.id],
                    'net_contribution': value, 'per_q': value/q, 'equation': r.reaction})
            self.tables['ledger_totals'].append({'witness': label, 'metabolite': mid,
                'residual': sum(contributions), 'tolerance': self.config['tolerance']})
            assert abs(sum(contributions)) <= self.config['tolerance']
        # Route sum excludes artificial arrival and downstream whole-cell metabolism.
        # All omitted proton contributions remain visible in the complete ledger.
        parts = {'A_hydrolysis': [HYD[scenario]], 'B_dipeptide_transport': ['R2038'],
                 'C_amino_acid_recovery': ['R2030', 'R2040'],
                 'D_pump': ['R795'], 'water_supply': ['R1363']}
        totals = {}
        for part, ids in parts.items():
            terms = {}
            for rid in ids:
                for mid, c in columns[rid].items():
                    terms[mid] = terms.get(mid, 0) + c*flux[rid]/q
                    totals[mid] = totals.get(mid, 0) + c*flux[rid]/q
            self.tables['route_net_reactions'].append({'witness': label, 'component': part,
                'reaction_fluxes': {rid: flux[rid] for rid in ids},
                'stoichiometry_per_q': {m: v for m, v in terms.items() if abs(v) > 1e-9}})
        self.tables['route_net_reactions'].append({'witness': label, 'component': 'route_total',
            'reaction_fluxes': {rid: flux[rid] for ids in parts.values() for rid in ids},
            'stoichiometry_per_q': {m: v for m, v in totals.items() if abs(v) > 1e-9}})
        h = 'm1007[C_va]'
        route_ids = {rid for ids in parts.values() for rid in ids}
        other_h = {rid: stoich.get(h, 0)*flux[rid] for rid, stoich in columns.items()
                   if h in stoich and rid not in route_ids}
        glypro = {rid: sum(stoich.get(m, 0) for m in ORIGIN.values())
                  + stoich.get('m1864[C_ex]', 0) for rid, stoich in columns.items()}
        self.witnesses[label].update(
            other_vacuole_proton_contributions=other_h,
            glypro_total_pool_contributions={r: c*flux[r] for r, c in glypro.items() if c},
            pump_atp_per_q=-columns['R795']['m141[C_cy]']*flux['R795']/q,
            includes_any_background_pump=True)
        write(self.out/'route_flux_witnesses.json', self.witnesses)
        self.save()


def energy(run, models):
    for scenario in ('REF', 'V_ONLY', 'C_ONLY'):
        _, changes = close_model(models[scenario], run.config)
        # Reuse the established closure definition without inheriting LP-text copy rounding.
        closed = clone(models[scenario], run.config)
        for change in changes:
            closed.reactions.get_by_id(change['reaction']).bounds = change['after']
        assert len(closed.constraints) == len(closed.metabolites)
        assert all(c.lb == c.ub == 0 for c in closed.constraints)
        assert all(r.lower_bound <= 0 <= r.upper_bound for r in closed.reactions)
        run.manifest.setdefault('closed_changes', {})[scenario] = changes
        zero = clone(closed, run.config)
        for r in zero.reactions:
            r.bounds = (0, 0)
        zero.objective = zero.problem.Objective(0)
        label = scenario+'_closed_zero'
        row, flux = run.solve(zero, label, 'exact zero-flux feasibility')
        run.tables['closed_energy_results'].append({'scenario': scenario, 'carrier': 'zero',
            'status': row['status'], 'maximum': row['objective'],
            'verdict': 'feasible' if row['status'] == 'optimal' else 'unresolved'})
        for ntp in ('ATP', 'GTP', 'UTP', 'CTP'):
            test = clone(closed, run.config)
            if ntp == 'ATP':
                objective(test, run.config['maintenance'])
                diss = test.reactions.get_by_id(run.config['maintenance'])
            else:
                diss = add_ntp_dissipation(test, ntp)
            label = scenario+'_closed_'+ntp
            row, flux = run.solve(test, label, 'maximum closed '+ntp+' dissipation')
            verdict = energy_verdict(row['status'], row['objective'], run.config['tolerance'])
            run.tables['closed_energy_results'].append({'scenario': scenario, 'carrier': ntp,
                'status': row['status'], 'maximum': row['objective'], 'verdict': verdict})
            if verdict == 'positive_energy_regeneration':
                diss.bounds = (min(1., row['objective']),)*2
                l1(test)
                run.solve(test, label+'_risk_witness', 'L1 positive closed energy witness')
        run.save()


def primary(run, models):
    growths = {}
    cases = [(None, s) for s in ('REF', 'TEMPLATE', 'V_ONLY', 'C_ONLY', 'NONE')]
    cases += [(o, s) for o in ORIGIN for s in ('C_ONLY', 'V_ONLY', 'NONE')]
    for origin, scenario in cases:
        name = (origin or 'NO_SOURCE')+'_'+scenario
        model = supplied(models[scenario], run.config, origin)
        objective(model, 'biomass_C')
        row, flux = run.solve(model, name+'_growth', 'maximum biomass; optional source')
        run.record('origin_localization_results', name+'_growth', origin, scenario, model, row, flux)
        run.require(row)
        growths[(origin, scenario)] = row['objective']
        p = clone(model, run.config)
        growth_floor(p, row['objective'])
        l1(p)
        pr, pf = run.solve(p, name+'_pfba', 'minimum split total flux at maximum growth; not ATP cost')
        run.require(pr)
        run.record('origin_localization_results', name+'_pfba', origin, scenario, p, pr, pf,
                   growth_floor=row['objective'], objective_is_total_flux=True)
        if scenario in HYD:
            objective(model, HYD[scenario])
            hr, hf = run.solve(model, name+'_hydro_max', 'maximum hydrolysis; normal NGAM, no growth floor')
            run.record('origin_localization_results', name+'_hydro_max', origin, scenario, model, hr, hf,
                       growth_requirement='no extra growth floor')
    for origin in ORIGIN:
        for scenario in ('C_ONLY', 'V_ONLY', 'NONE'):
            delta = growths[(origin, scenario)] - growths[(None, scenario)]
            assert delta >= -run.config['tolerance'], (origin, scenario, delta)
    assert abs(growths[(None, 'TEMPLATE')]-growths[(None, 'REF')]) <= run.config['tolerance']
    run.manifest['no_source_growth'] = {s: v for (o, s), v in growths.items() if o is None}
    return growths


def matched(run, models, growths):
    points = {}
    for origin in ORIGIN:
        pair = {}
        q = run.config['main_epsilon']
        for attempt, q in enumerate((q, *run.config['fallback_q'])):
            pair = {}
            for scenario, hydro in HYD.items():
                model = supplied(models[scenario], run.config, origin, q, hydro)
                objective(model, 'biomass_C')
                label = f'{origin}_{scenario}_fixedq_{attempt}_growth'
                row, flux = run.solve(model, label, 'fixed equal source and hydrolysis; maximum growth')
                run.record('matched_throughput_results', label, origin, scenario, model, row, flux, q=q)
                pair[scenario] = (model, row)
            if all(r['status'] == 'optimal' for m, r in pair.values()):
                break
        else:
            run.manifest.setdefault('incomparable_origins', []).append(origin)
            continue
        reference = growths[(None, 'REF')]
        can_reference = []
        for scenario, (model, row) in pair.items():
            test = clone(model, run.config)
            growth_floor(test, reference)
            test.objective = test.problem.Objective(0)
            label = origin+'_'+scenario+'_reference_growth_feasibility'
            rr, ff = run.solve(test, label, 'fixed q at same no-source reference growth')
            can_reference.append(rr['status'] == 'optimal')
            run.record('matched_throughput_results', label, origin, scenario, test, rr, ff, q=q, growth_floor=reference)
        common = reference if all(can_reference) else min(r['objective'] for m, r in pair.values())-run.config['strict_growth_delta']
        run.manifest.setdefault('common_operating_points', {})[origin] = {
            'q': q, 'growth_floor': common, 'reference_growth': reference,
            'reference_feasible_both': all(can_reference)}
        for scenario, (model, row) in pair.items():
            growth_floor(model, common)
            extrema = {}
            for sense in ('min', 'max'):
                objective(model, 'R795', sense)
                label = origin+'_'+scenario+'_pump_'+sense
                rr, ff = run.solve(model, label, 'R795 range at equal q and common growth floor')
                run.require(rr)
                extrema[sense] = rr['objective']
                run.record('matched_throughput_results', label, origin, scenario, model, rr, ff, q=q, growth_floor=common)
            model.reactions.R795.upper_bound = min(model.reactions.R795.upper_bound,
                extrema['min']+run.config['pump_minimum_slack'])
            l1(model)
            label = origin+'_'+scenario+'_matched_witness'
            rr, ff = run.solve(model, label, 'L1 witness near minimum R795, fixed q/common growth')
            run.require(rr)
            run.record('matched_throughput_results', label, origin, scenario, model, rr, ff,
                       q=q, growth_floor=common, pump_min=extrema['min'], pump_max=extrema['max'])
            run.witness(model, ff, label, origin, scenario, q, common)
            points[(origin, scenario)] = (q, common)
    # Reference/background is measured at the same no-source growth requirement.
    for sense in ('min', 'max'):
        test = clone(models['REF'], run.config)
        growth_floor(test, growths[(None, 'REF')])
        objective(test, 'R795', sense)
        rr, ff = run.solve(test, 'NO_SOURCE_REF_pump_'+sense, 'background pump range at reference growth')
        run.record('matched_throughput_results', 'NO_SOURCE_REF_pump_'+sense, None, 'REF', test, rr, ff,
                   q=0, growth_floor=growths[(None, 'REF')])
    return points


def dependence(run, models, points):
    for origin in ORIGIN:
        for scenario, hydro in HYD.items():
            prefix = origin+'_'+scenario
            free = supplied(models[scenario], run.config, origin)
            free.reactions.R795.bounds = (0, 0)
            for target, suffix in (('biomass_C', 'growth'), (hydro, 'hydro_max')):
                objective(free, target)
                label = prefix+'_pump_off_free_'+suffix
                rr, ff = run.solve(free, label, 'R795 closed; optional source; no extra growth floor')
                run.record('vatpase_dependency', label, origin, scenario, free, rr, ff,
                           condition='free_supply', closed='R795')
            if (origin, scenario) not in points:
                continue
            q, common = points[(origin, scenario)]
            fixed = supplied(models[scenario], run.config, origin, q, hydro)
            fixed.reactions.R795.bounds = (0, 0)
            with_growth = clone(fixed, run.config)
            growth_floor(with_growth, common)
            l1(with_growth)
            label = prefix+'_pump_off_fixedq_common_growth'
            rr, ff = run.solve(with_growth, label, 'pump-off positive q at common growth; L1 feasibility witness')
            run.record('vatpase_dependency', label, origin, scenario, with_growth, rr, ff,
                       condition='fixed_q_common_growth', closed='R795', q=q, growth_floor=common)
            if rr['status'] == 'optimal':
                run.witness(with_growth, ff, label, origin, scenario, q, common)
            else:
                objective(fixed, 'biomass_C')
                label = prefix+'_pump_off_fixedq_no_growth_floor'
                rr, ff = run.solve(fixed, label, 'pump-off fixed q, normal NGAM, no extra growth requirement')
                run.record('vatpase_dependency', label, origin, scenario, fixed, rr, ff,
                           condition='fixed_q_no_growth_floor', closed='R795', q=q)
                if rr['status'] == 'optimal':
                    run.witness(fixed, ff, label, origin, scenario, q, None)
            # Two cross-compartment links, plus product recovery for vacuolar hydrolysis.
            blocks = (['R2038'] if (origin, scenario) in (('O_CY','V_ONLY'), ('O_VA','C_ONLY')) else [])
            if scenario == 'V_ONLY':
                blocks += ['R2030', 'R2040']
            for blocked in blocks:
                test = supplied(models[scenario], run.config, origin, q, hydro)
                test.reactions.get_by_id(blocked).bounds = (0, 0)
                objective(test, 'biomass_C')
                label = prefix+'_block_'+blocked
                rr, ff = run.solve(test, label, 'transport isolation at fixed q, normal NGAM, no extra growth floor')
                run.record('transport_isolation', label, origin, scenario, test, rr, ff,
                           closed=blocked, q=q, condition='fixed_q_no_growth_floor')


def main():
    parser = argparse.ArgumentParser(__doc__)
    for name in ('source', 'candidates', 'config', 'output', 'budget'):
        parser.add_argument('--'+name, required=True)
    args = parser.parse_args()
    source, candidates, config_path = map(workspace, (args.source, args.candidates, args.config))
    config = json.loads(config_path.read_text())
    assert sha(source) == SOURCE_SHA256 == config['source_sha256']
    for path, expected in config['input_configuration_sha256'].items():
        assert sha(workspace(path)) == expected, path
    build = json.loads((candidates/'manifest.json').read_text())
    assert Path(build['source']).resolve() == source and build['source_sha256'] == sha(source)
    baseline = read_sbml_model(source)
    paths = {'REF': source, **{s: candidates/(s+'.xml') for s in ('TEMPLATE','V_ONLY','C_ONLY','NONE')}}
    for scenario, path in paths.items():
        if scenario != 'REF':
            manifest = json.loads(path.with_suffix('.build.json').read_text())
            assert manifest['output_sha256'] == sha(path), scenario
            assert Path(manifest['output']).resolve() == path
            assert manifest['source_sha256'] == sha(source)
            assert build['candidates'][scenario]['sha256'] == sha(path)
            verify_candidate(baseline, read_sbml_model(path), scenario)
    identity = {'source': str(source), 'source_sha256': sha(source),
                'source_build_sha256': sha(source.with_suffix('.build.json')),
                'candidate_manifest_sha256': sha(candidates/'manifest.json'),
                'files': {s: {'path': str(p), 'sha256': sha(p)} for s,p in paths.items()}}
    run = Experiment(config, workspace(args.output), workspace(args.budget), config_path, identity)
    run.tables = {k: [] for k in ('closed_energy_results','origin_localization_results',
        'matched_throughput_results','vatpase_dependency','transport_isolation',
        'proton_water_atp_ledger','ledger_totals','route_net_reactions')}
    run.witnesses = {}
    # Snapshot every executing local module recorded by the existing runtime helper.
    snap = run.out/'code'; snap.mkdir()
    for path in run.manifest['source_sha256']:
        dst = snap/path; dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(workspace(path), dst)
    try:
        models = {s: run.loaded(p, s) for s,p in paths.items()}
        original = {s: model_definition(m) for s,m in models.items()}
        original_problems = {s: signature(m) for s,m in models.items()}
        energy(run, models)
        growths = primary(run, models)
        points = matched(run, models, growths)
        dependence(run, models, points)
        assert all(model_definition(models[s]) == d for s,d in original.items())
        assert all(signature(models[s]) == d for s,d in original_problems.items())
        assert all(sha(paths[s]) == identity['files'][s]['sha256'] for s in paths)
        run.manifest.update(status='complete', temporary_edits_restored=True,
            closed_template_growth_matches_reference=True,
            optional_scenarios='BOTH/extracellular/extra q not run unless predefined fallback triggered')
    except Exception as exc:
        run.manifest.update(status='failed', error=repr(exc))
        raise
    finally:
        run.save()


if __name__ == '__main__':
    main()
