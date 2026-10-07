"""Independent saved-matrix audit; never construct an optimizer or call a solver."""
from __future__ import annotations

import argparse
import copy
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
from xml.etree import ElementTree as ET

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

ROOT = Path(__file__).resolve().parents[1]
CY = 'HYP_GLYPRO_HYD_CY'
SOURCE = 'DIAG_GLYPRO_SOURCE'
HYD = {'C_ONLY': CY, 'V_ONLY': 'R2039'}
ORIGINS = {'O_CY': 'm1865[C_cy]', 'O_VA': 'm1866[C_va]'}
POOLS = {'m1864[C_ex]', *ORIGINS.values()}
LEDGER_POOLS = POOLS | {'m272[C_cy]', 'm1863[C_va]', 'm765[C_cy]', 'm1867[C_va]',
    'm32[C_cy]', 'm1384[C_va]', 'm10[C_cy]', 'm1007[C_va]',
    'm141[C_cy]', 'm143[C_cy]', 'm35[C_cy]'}


def load(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def check_candidate_xml(source, candidate_directory):
    """Inspect exported SBML directly, independently of the builder/COBRA guards."""
    ns = '{http://www.sbml.org/sbml/level3/version1/fbc/version2}'
    def normalized(element):
        text = element.text if element.text and element.text.strip() else ''
        return (element.tag, dict(element.attrib), text, [normalized(c) for c in element])
    def model(path):
        return next(c for c in ET.parse(path).getroot() if c.tag.endswith('}model'))
    def parts(element):
        return {c.tag.rsplit('}', 1)[-1]: c for c in element}
    original = model(source); original_parts = parts(original)
    original_reactions = {r.get('id'): r for r in original_parts['listOfReactions']}
    species = {s.get('id'): s for s in original_parts['listOfSpecies']}
    original_parameters = {p.get('id'): p for p in original_parts['listOfParameters']}
    result = {}
    for scenario in ('TEMPLATE', 'V_ONLY', 'C_ONLY', 'NONE'):
        path = candidate_directory/(scenario+'.xml')
        candidate = model(path); blocks = parts(candidate)
        assert candidate.attrib == original.attrib and blocks.keys() == original_parts.keys()
        for block in blocks.keys() - {'listOfReactions', 'listOfParameters'}:
            assert normalized(blocks[block]) == normalized(original_parts[block]), (scenario, block)
        reactions = {r.get('id'): r for r in blocks['listOfReactions']}
        assert reactions.keys() - original_reactions.keys() == {'R_'+CY}
        assert original_reactions.keys() <= reactions.keys()
        parameters = {p.get('id'): p for p in blocks['listOfParameters']}
        for pid, before in original_parameters.items():
            assert normalized(parameters[pid]) == normalized(before)
        allowed_new_parameters = {'R_R2039_upper_bound'} if scenario == 'V_ONLY' else (
            {'R_'+CY+'_upper_bound'} if scenario == 'C_ONLY' else set())
        assert parameters.keys() - original_parameters.keys() == allowed_new_parameters
        for pid in allowed_new_parameters:
            assert float(parameters[pid].get('value')) == .01
        for rid, before in original_reactions.items():
            after = copy.deepcopy(reactions[rid])
            if rid == 'R_R2039':
                expected = 1000. if scenario == 'TEMPLATE' else .01 if scenario == 'V_ONLY' else 0.
                assert float(parameters[after.get(ns+'lowerFluxBound')].get('value')) == 0.
                assert float(parameters[after.get(ns+'upperFluxBound')].get('value')) == expected
                after.set(ns+'upperFluxBound', before.get(ns+'upperFluxBound'))
            assert normalized(after) == normalized(before), (scenario, rid)
        new = reactions['R_'+CY]
        assert new.get('reversible') == 'false'
        assert float(parameters[new.get(ns+'lowerFluxBound')].get('value')) == 0.
        assert float(parameters[new.get(ns+'upperFluxBound')].get('value')) == (.01 if scenario == 'C_ONLY' else 0.)
        assert not any(e.tag.endswith('}geneProductAssociation') for e in new.iter())
        columns = {}
        for tag, sign in (('listOfReactants', -1), ('listOfProducts', 1)):
            for member in parts(new)[tag]:
                columns[member.get('species')] = sign*float(member.get('stoichiometry'))
        assert columns == {'M_m1865__91__C_cy__93__': -1., 'M_m32__91__C_cy__93__': -1.,
                           'M_m272__91__C_cy__93__': 1., 'M_m765__91__C_cy__93__': 1.}
        balance, charge = {}, 0.
        for mid, coefficient in columns.items():
            formula = species[mid].get(ns+'chemicalFormula')
            terms = re.findall(r'([A-Z][a-z]?)([0-9]*)', formula)
            assert ''.join(e+n for e, n in terms) == formula
            for element, count in terms:
                balance[element] = balance.get(element, 0.) + coefficient*int(count or '1')
            charge += coefficient*float(species[mid].get(ns+'charge'))
        assert all(value == 0 for value in balance.values()) and charge == 0
        result[scenario] = {'path': str(path), 'sha256': digest(path),
            'original_species_genes_groups_objective_metadata_unchanged': True,
            'original_reactions_unchanged_except_declared_R2039_upper_bound': True,
            'new_reactions': [CY], 'new_gene_product_association': False,
            'element_residual': balance, 'charge_residual': charge}
    return result


def maximum(values):
    return float(np.max(values, initial=0.))


def close(actual, expected, tolerance):
    assert abs(actual - expected) <= tolerance, (actual, expected)


def check_problem(data, tolerance):
    """Compare the saved native LP against net S, then check its saved primal."""
    p, native, result = data['problem'], data['actual_solver'], data['result']
    cs = p['constraints']
    names, variables = cs['names'], cs['variables']
    assert len(set(names)) == len(names) and len(set(variables)) == len(variables)
    assert variables == native['variables']
    rows, columns = {n: i for i, n in enumerate(names)}, {n: i for i, n in enumerate(variables)}
    matrix = csr_matrix((cs['matrix_data'], cs['matrix_indices'], cs['matrix_indptr']),
                        shape=(len(names), len(variables)))
    assert set(p['metabolites']) <= set(rows)
    extra = set(rows) - set(p['metabolites'])
    assert extra <= {'DIAG_COMMON_GROWTH'}, extra
    assert len(variables) == 2 * len(p['reactions'])
    rr, cc, vv, pairs = [], [], [], {}
    for rid, (lower, upper, gpr, stoich) in p['reactions'].items():
        reverse = rid + '_reverse_' + hashlib.md5(rid.encode()).hexdigest()[:5]
        f, r = columns[rid], columns[reverse]
        pairs[rid] = (f, r)
        for mid, coefficient in stoich.items():
            rr.extend((rows[mid], rows[mid])); cc.extend((f, r)); vv.extend((coefficient, -coefficient))
        assert lower <= upper
        # Split bounds must encode the actual signed bounds, including fixed positive q.
        expected = (max(0., lower), max(0., upper), max(0., -upper), max(0., -lower))
        actual = (native['lower_bounds'][f], native['upper_bounds'][f],
                  native['lower_bounds'][r], native['upper_bounds'][r])
        assert tuple(actual) == expected, (rid, actual, expected)
    for mid in p['metabolites']:
        idx = rows[mid]
        assert cs['bounds'][idx] == [0, 0]
        assert native['constraint_sense'][idx] == '=' and native['constraint_rhs'][idx] == 0
    if extra:
        idx = rows['DIAG_COMMON_GROWTH']
        f, r = pairs['biomass_C']
        rr.extend((idx, idx)); cc.extend((f, r)); vv.extend((1., -1.))
        assert cs['bounds'][idx][1] is None
        assert native['constraint_sense'][idx] == '>'
        close(cs['bounds'][idx][0], native['constraint_rhs'][idx], tolerance)
    expected_matrix = coo_matrix((vv, (rr, cc)), shape=matrix.shape).tocsr()
    matrix_error = maximum(np.abs((matrix - expected_matrix).data))
    assert matrix_error <= 1e-12, matrix_error
    assert native['objective_sense'] == (1 if p['direction'] == 'min' else -1)
    summary = {'status': result['status'], 'native_rows': matrix.shape[0],
               'split_variables': matrix.shape[1], 'matrix_vs_net_S_error': matrix_error}
    if result['status'] != 'optimal':
        assert data['fluxes'] is None and native['primal'] is None
        assert result['status'] == 'infeasible', result['status']
        return summary
    x = np.asarray(native['primal'], dtype=float)
    assert x.shape == (len(variables),) and np.all(np.isfinite(x))
    lb, ub = np.asarray(native['lower_bounds'], float), np.asarray(native['upper_bounds'], float)
    activity, rhs = matrix @ x, np.asarray(native['constraint_rhs'], float)
    residual = activity - rhs
    sense = np.asarray(native['constraint_sense'])
    violation = np.where(sense == '=', np.abs(residual),
                        np.where(sense == '<', np.maximum(residual, 0), np.maximum(-residual, 0)))
    bound_error = max(maximum(lb - x), maximum(x - ub))
    row_error = maximum(violation)
    flux = data['fluxes']
    assert set(flux) == set(p['reactions']) and all(math.isfinite(v) for v in flux.values())
    net_error = max(abs(flux[rid] - (x[f] - x[r])) for rid, (f, r) in pairs.items())
    # Recompute signed net S*v separately, including compensated summation.
    terms = {mid: [] for mid in p['metabolites']}
    net_bounds = []
    for rid, (lower, upper, gpr, stoich) in p['reactions'].items():
        net_bounds.append(max(lower-flux[rid], flux[rid]-upper, 0.))
        for mid, coefficient in stoich.items():
            terms[mid].append(coefficient*flux[rid])
    net_mass = max(abs(math.fsum(values)) for values in terms.values())
    net_mass_naive = max(abs(sum(values)) for values in terms.values())
    net_bound = max(net_bounds)
    objective = float(np.dot(np.asarray(native['objective_coefficients']), x) + native['objective_constant'])
    objective_error = abs(objective - result['objective'])
    assert max(bound_error, row_error, net_error, objective_error, net_mass, net_mass_naive, net_bound) <= tolerance, (
        bound_error, row_error, net_error, objective_error, net_mass, net_mass_naive, net_bound)
    summary.update(max_split_bound_violation=bound_error, max_native_row_violation=row_error,
                   max_net_split_difference=net_error, objective_reconstruction_error=objective_error,
                   max_net_mass_residual=net_mass, max_net_mass_residual_naive_sum=net_mass_naive,
                   max_net_bound_violation=net_bound)
    return summary


def check_witness(label, witness, data, tolerance):
    flux, reactions = witness['fluxes'], data['problem']['reactions']
    assert flux == data['fluxes'] and data['result']['status'] == 'optimal'
    q, origin, scenario = witness['q'], witness['origin'], witness['scenario']
    assert q > 100 * tolerance
    assert reactions[SOURCE][3] == {ORIGINS[origin]: 1}
    for rid in (SOURCE, HYD[scenario]):
        assert reactions[rid][:2] == [q, q]
        close(flux[rid], q, tolerance)
    other = 'R2039' if scenario == 'C_ONLY' else CY
    assert reactions[other][:2] == [0, 0]
    assert len([r for r in reactions if r == SOURCE]) == 1
    contributions = {rid: sum(st.get(mid, 0) for mid in POOLS) * flux[rid]
                     for rid, (_, _, _, st) in reactions.items()}
    close(sum(contributions.values()), 0, tolerance)
    # No second net Gly-Pro source, disappearance, or export at fixed matched throughput.
    unexpected = {rid: value for rid, value in contributions.items()
                  if rid not in (SOURCE, HYD[scenario]) and abs(value) > tolerance}
    assert not unexpected, unexpected
    hva = {rid: st.get('m1007[C_va]', 0) * flux[rid]
           for rid, (_, _, _, st) in reactions.items() if 'm1007[C_va]' in st}
    close(sum(hva.values()), 0, tolerance)
    pump_cost = -reactions['R795'][3]['m141[C_cy]'] * flux['R795'] / q
    close(pump_cost, witness['pump_atp_per_q'], tolerance)
    if '_pump_off_' in label:
        assert reactions['R795'][:2] == [0, 0]
        close(flux['R795'], 0, tolerance)
    growth = witness['growth_floor']
    if growth is not None:
        assert flux['biomass_C'] >= growth - tolerance
    other_h = {r: v for r, v in hva.items() if r not in ('R795', 'R2038', 'R2030', 'R2040') and abs(v) > tolerance}
    return {'q': q, 'growth_floor': growth, 'biomass': flux['biomass_C'],
            'source': flux[SOURCE], 'hydrolysis': flux[HYD[scenario]], 'R795': flux['R795'],
            'pump_ATP_per_q': pump_cost, 'unexpected_glypro_fates': unexpected,
            'other_vacuolar_proton_contributions': other_h}


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--run', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--refresh', action='store_true', help='Recheck and update this audit for the same immutable run')
    args = parser.parse_args()
    run, out = args.run.resolve(), args.output.resolve()
    assert run.is_relative_to(ROOT) and out.is_relative_to(ROOT)
    assert not out.exists() or args.refresh, 'Audit output must be new, or explicitly refreshed'
    if out.exists():
        assert load(out/'audit.json')['run'] == str(run), 'Refresh may only update this same run audit'
    manifest = load(run/'manifest.json')
    assert manifest['status'] == 'complete'
    config, tolerance = manifest['config'], manifest['config']['tolerance']
    results, inputs = {}, {}
    for item in manifest['outcomes']:
        label = item['label']
        # The inherited outcome label includes run directory; file names do not.
        label = label.rsplit('/', 1)[-1]
        path = run/(label+'.json')
        data = load(path)
        results[label] = check_problem(data, tolerance)
        inputs[str(path.relative_to(ROOT))] = digest(path)
        assert data['result']['solver'] == config['solver']
    assert len(results) == manifest['optimization_calls']
    budget = load(Path(manifest['budget_path']))
    assert len(budget['calls']) <= config['limits']['solves']
    assert budget['solve_call_seconds'] <= config['limits']['solve_wall_seconds']
    identity = manifest['explicit_inputs']
    assert Path(identity['source']) == Path(identity['files']['REF']['path'])
    assert digest(Path(identity['source'])) == identity['source_sha256'] == config['source_sha256']
    for record in identity['files'].values():
        assert digest(Path(record['path'])) == record['sha256']
    assert digest(Path(identity['source']).with_suffix('.build.json')) == identity['source_build_sha256']
    candidate_directory = Path(identity['files']['TEMPLATE']['path']).parent
    assert digest(candidate_directory/'manifest.json') == identity['candidate_manifest_sha256']
    xml_checks = check_candidate_xml(Path(identity['source']), candidate_directory)
    for path, expected in config['input_configuration_sha256'].items():
        assert digest(ROOT/path) == expected
    for path, expected in manifest['source_sha256'].items():
        assert digest(run/'code'/path) == expected
    assert manifest['temporary_edits_restored'] and manifest['closed_template_growth_matches_reference']
    # Mutation check: a corrupt saved primal must fail the independent checks.
    successful = next(label for label, r in results.items() if r['status'] == 'optimal')
    altered = copy.deepcopy(load(run/(successful+'.json')))
    altered['actual_solver']['primal'][0] += 0.001
    rejected = False
    try:
        check_problem(altered, tolerance)
    except AssertionError:
        rejected = True
    assert rejected, 'Corrupted primal was accepted'
    witnesses = load(run/'route_flux_witnesses.json')
    witness_results = {label: check_witness(label, w, load(run/(label+'.json')), tolerance)
                       for label, w in witnesses.items()}
    for origin in ORIGINS:
        point = manifest['common_operating_points'][origin]
        for scenario in HYD:
            w = witnesses[origin+'_'+scenario+'_matched_witness']
            assert w['q'] == point['q'] and w['growth_floor'] == point['growth_floor']
        # Audit fairness using the actual optimization definitions, not labels.
        a = load(run/(origin+'_C_ONLY_pump_min.json'))['problem']
        b = load(run/(origin+'_V_ONLY_pump_min.json'))['problem']
        assert a['metabolites'] == b['metabolites']
        assert a['constraints'] == b['constraints']
        assert a['objective'] == b['objective'] and a['direction'] == b['direction'] == 'min'
        assert a['reactions'].keys() == b['reactions'].keys()
        for rid in a['reactions']:
            assert a['reactions'][rid][2:] == b['reactions'][rid][2:]
            if rid not in (CY, 'R2039'):
                assert a['reactions'][rid][:2] == b['reactions'][rid][:2], rid
        for scenario in ('C_ONLY', 'V_ONLY', 'NONE'):
            supplied = load(run/(origin+'_'+scenario+'_growth.json'))
            no_source = load(run/('NO_SOURCE_'+scenario+'_growth.json'))
            assert supplied['result']['objective'] >= no_source['result']['objective'] - tolerance
    # Recompute every full ledger entry from net signed S, checking no adjacency omissions.
    ledger = read_table(run/'proton_water_atp_ledger.tsv')
    covered = {}
    for row in ledger:
        label, mid, rid = row['witness'], row['metabolite'], row['reaction']
        data = load(run/(label+'.json')) if label not in covered else covered[label]['data']
        cache = covered.setdefault(label, {'data': data, 'pairs': set()})
        assert (mid, rid) not in cache['pairs']; cache['pairs'].add((mid, rid))
        coefficient = data['problem']['reactions'][rid][3][mid]
        value = coefficient * data['fluxes'][rid]
        close(float(row['coefficient']), coefficient, tolerance)
        close(float(row['signed_flux']), data['fluxes'][rid], tolerance)
        close(float(row['net_contribution']), value, tolerance)
        close(float(row['per_q']), value / witnesses[label]['q'], tolerance)
    for label, cache in covered.items():
        selected = {mid for mid, rid in cache['pairs']}
        assert selected == LEDGER_POOLS
        expected = {(mid, rid) for rid, (_, _, _, st) in cache['data']['problem']['reactions'].items()
                    for mid in st if mid in selected}
        assert cache['pairs'] == expected
    assert covered.keys() == witnesses.keys()
    route_sums = []
    for row in read_table(run/'route_net_reactions.tsv'):
        label = row['witness']; data = covered[label]['data']; total = {}
        values = json.loads(row['reaction_fluxes']); declared = json.loads(row['stoichiometry_per_q'])
        for rid, value in values.items():
            close(value, data['fluxes'][rid], tolerance)
            for mid, coefficient in data['problem']['reactions'][rid][3].items():
                total[mid] = total.get(mid, 0.) + coefficient * value / witnesses[label]['q']
        total = {mid: v for mid, v in total.items() if abs(v) > 1e-9}
        assert total.keys() == declared.keys()
        for mid in total:
            close(total[mid], declared[mid], tolerance)
        if row['component'] == 'route_total':
            route_sums.append({'witness': label, 'stoichiometry_per_q': total})
            internal = POOLS - {ORIGINS[witnesses[label]['origin']]}
            internal |= {'m1863[C_va]', 'm1867[C_va]'}
            assert not {m: total[m] for m in internal & total.keys() if abs(total[m]) > tolerance}
    energies = read_table(run/'closed_energy_results.tsv')
    assert len(energies) == 15
    for row in energies:
        label = row['scenario']+'_closed_'+row['carrier']
        data = load(run/(label+'.json'))
        assert SOURCE not in data['problem']['reactions']
        assert set(data['problem']['constraints']['names']) == set(data['problem']['metabolites'])
        for rid, (lb, ub, gpr, st) in data['problem']['reactions'].items():
            if not any(v < 0 for v in st.values()) or not any(v > 0 for v in st.values()) or rid in config['biomass_reactions']:
                assert lb == ub == 0
        assert data['problem']['reactions']['xMAINTENANCE'][0] == 0
        status, value = data['result']['status'], data['result']['objective']
        expected = ('feasible' if status == 'optimal' else 'unresolved') if row['carrier'] == 'zero' else (
            'within_tolerance' if status == 'optimal' and value is not None and abs(value) <= tolerance else
            'positive_energy_regeneration' if status == 'optimal' and value is not None and value > tolerance else 'unresolved')
        assert row['status'] == status and row['verdict'] == expected
    multiplicity = []
    for origin in ORIGINS:
        for scenario in HYD:
            label = origin+'_'+scenario
            a, b = run/(label+'_growth.json'), run/(label+'_pfba.json')
            x, y = load(a), load(b)
            close(x['fluxes']['biomass_C'], y['fluxes']['biomass_C'], tolerance)
            different = {rid: y['fluxes'][rid]-value for rid, value in x['fluxes'].items()
                         if abs(y['fluxes'][rid]-value) > tolerance}
            assert different
            multiplicity.append({'origin': origin, 'scenario': scenario,
                'maximum_growth_FBA': x['fluxes']['biomass_C'], 'growth_in_pFBA': y['fluxes']['biomass_C'],
                'number_different_net_fluxes': len(different), 'net_flux_difference_pFBA_minus_FBA': different,
                'growth_problem_sha256': digest(a), 'pfba_problem_sha256': digest(b),
                'matched_fixed_q_R795_min': load(run/(label+'_pump_min.json'))['result']['objective'],
                'matched_fixed_q_R795_max': load(run/(label+'_pump_max.json'))['result']['objective']})
    out.mkdir(parents=True, exist_ok=args.refresh)
    result = {'status': 'passed', 'verified_utc': datetime.now(timezone.utc).isoformat(),
              'method': 'read-only independent reconstruction of native sparse matrices and saved net/split primals; no optimizations',
              'run': str(run), 'manifest_sha256': digest(run/'manifest.json'), 'auditor_sha256': digest(Path(__file__)),
              'audited_solve_files_sha256': inputs, 'solves': results, 'witnesses': witness_results,
              'independent_direct_SBML_checks': xml_checks,
              'alternative_growth_optimal_flux_solutions': multiplicity,
              'route_net_reactions_recomputed': route_sums, 'energy_rows': energies,
              'ledger_rows_verified': len(ledger), 'corrupted_primal_rejected': rejected,
              'optimization_calls_performed_by_audit': 0,
              'limits': ['Solver optimal/infeasible status is retained; independent primal feasibility is checked, not independent dual/Farkas certificates.',
                         'Conditional route sums retain any background pump included in the saved witness; no total-cell ATP cost or biological localization claim.']}
    (out/'audit.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    (out/'optimal_flux_multiplicity.json').write_text(json.dumps({
        'method': 'Read-only comparison of net FBA/pFBA fluxes at the same recorded maximum biomass, within the existing tolerance; no solver calls.',
        'tolerance': tolerance, 'optimization_calls': 0, 'rows': multiplicity,
        'inference': 'Observed whole-cell optimal flux multiplicity is separate from the target-pump FVA result and does not establish biological enzyme localization.'}, indent=2)+'\n')
    optimal = [r for r in results.values() if r['status'] == 'optimal']
    text = (f'# Independent numerical audit\n\nAudited {len(results)} saved solver outcomes, including {len(optimal)} optimal primals, '
            f'{len(witnesses)} full route witnesses and {len(ledger)} ledger entries. All checks passed. No optimization was run.\n\n'
            f'Maximum native constraint violation: {max(r["max_native_row_violation"] for r in optimal):.6g}. '
            f'Independent full net S*v residual (compensated summation): {max(r["max_net_mass_residual"] for r in optimal):.6g}. '
            f'Maximum net reaction-bound violation: {max(r["max_net_bound_violation"] for r in optimal):.6g}; '
            f'maximum native split-variable bound violation: {max(r["max_split_bound_violation"] for r in optimal):.6g}. '
            f'Maximum net/split difference: {max(r["max_net_split_difference"] for r in optimal):.6g}. '
            'Native sparse coefficients reproduce the stored net stoichiometry and the sole permitted growth-floor constraint. '
            'Source and hydrolysis are fixed to the same positive q, paired routes share the same growth floor, and complete Gly-Pro pool accounting has no unexplained source or export.\n\n'
            'All ledger terms and weighted route sums were independently recomputed. Every adjacent reaction for each ledger metabolite is represented. '
            'Closed energy tests retain all mass-balance rows and have no diagnostic source or growth floor; verdicts agree with actual statuses and objectives. '
            'The saved code and explicit input identities match their recorded hashes. A deliberately corrupted primal was rejected.\n\n'
            'An independent ElementTree comparison of all four exported XML files confirms that original species, gene products, groups, objectives and metadata are unchanged, '
            'all original reactions retain their full XML definition apart from the declared R2039 capacity/isolation bound, and the only new reaction is the balanced, GPR-free cytosolic hypothesis. '
            'Historical chemical and energy notes are preserved with the original XML elements.\n\n'
            'Primary FBA and growth-optimal pFBA differ in '+', '.join(str(r['number_different_net_fluxes']) for r in multiplicity)+
            ' reaction net fluxes above 1e-7, for O_CY/C_ONLY, O_CY/V_ONLY, O_VA/C_ONLY and O_VA/V_ONLY respectively. '
            'Their recorded biomass values match. Whole-cell optimal fluxes are therefore not unique within the recorded numerical precision; '
            'target-pump FVA at matched q/growth is a separate result. Detailed differences and source hashes are in optimal_flux_multiplicity.json.\n\n'
            'Scope: this validates the recorded optimization problems and primal evidence. It does not independently establish optimality or infeasibility through dual/Farkas certificates, '
            'and it does not verify native protein activity, localization, pH, membrane potential, or gene essentiality.\n')
    (out/'REPORT.md').write_text(text)
    print(json.dumps({'status': 'passed', 'outcomes': len(results), 'witnesses': len(witnesses), 'ledger_rows': len(ledger), 'output': str(out)}))


if __name__ == '__main__':
    main()
