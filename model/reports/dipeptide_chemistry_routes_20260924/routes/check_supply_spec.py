"""Validate a disabled supply specification; never edit a model or activate a route."""
import argparse
from collections import Counter
import copy
import hashlib
import json
import math
from pathlib import Path

TARGETS = ('GD', 'GE', 'AG', 'GP')
AA = set('ACDEFGHIKLMNPQRSTVWY')


def number(value, name, upper=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError(name + ': finite nonnegative numeric value required')
    if upper is not None and value > upper:
        raise ValueError(name + ': value exceeds ' + str(upper))
    return value


def precursor_supply(precursor):
    """Check proposed nonoverlapping cleavage patterns and their complete residue ledger."""
    # ponytail: canonical unmodified chains only; PTMs need an explicit atom/charge ledger before extension.
    seq = precursor['mature_sequence']
    if not seq or set(seq)-AA or hashlib.sha256(seq.encode()).hexdigest() != precursor['sequence_sha256']:
        raise ValueError('Explicit unmodified mature sequence or its SHA differs')
    if precursor['abundance_units'] != 'mmol_chain/gDW' or precursor['rate_units'] != '1/h':
        raise ValueError('Protein abundance and degradation-rate units required; TPM is not a rate')
    for field in ('source_identity', 'condition', 'synthesis_path', 'synthesis_cost_source',
                  'localization_evidence', 'cleavage_evidence', 'balanced_chemical_equation',
                  'replacement_accounting'):
        if not precursor.get(field):
            raise ValueError('Missing precursor field: ' + field)
    costs = precursor['replacement_synthesis_energy_per_chain']
    if not costs or sum(number(v, 'replacement cost') for v in costs.values()) <= 0:
        raise ValueError('Replacement synthesis may not be a free source')
    rate = (number(precursor['abundance'], 'abundance') * number(precursor['degradation_rate'], 'degradation rate')
            * number(precursor['vacuolar_fraction'], 'vacuolar fraction', 1))
    patterns = precursor['cleavage_patterns']
    if not patterns or abs(sum(number(p['fraction'], 'pattern fraction', 1) for p in patterns)-1) > 1e-9:
        raise ValueError('Alternative cleavage-pattern fractions must sum to one')
    yields = Counter()
    for pattern in patterns:
        used, products = set(), Counter()
        for cut in pattern['target_pairs']:
            start, peptide = cut['start_zero_based'], cut['peptide']
            if isinstance(start, bool) or not isinstance(start, int) or start < 0 or start+2 > len(seq):
                raise ValueError('Target pair lies outside mature sequence')
            if peptide not in TARGETS or seq[start:start+2] != peptide:
                raise ValueError('Exact ordered target pair differs from sequence')
            positions = {start, start+1}
            if used & positions:
                raise ValueError('Overlapping peptide assignments double count residues')
            used |= positions
            products[peptide] += 1
        residual = Counter(seq[i] for i in range(len(seq)) if i not in used)
        remaining = Counter()
        for product in pattern['remaining_products']:
            product_seq, count = product['sequence'], product['moles_per_chain']
            if not product_seq or set(product_seq)-AA or not product.get('destination'):
                raise ValueError('Every remaining residue/peptide needs identity and destination')
            if isinstance(count, bool) or not isinstance(count, int) or count < 0:
                raise ValueError('Discrete cleavage pattern requires nonnegative integer product counts')
            remaining.update({aa: n*count for aa, n in Counter(product_seq).items()})
        if +remaining != +residual:
            raise ValueError('Remaining-products ledger does not conserve all residues')
        for peptide, count in products.items():
            yields[peptide] += pattern['fraction']*count
    # All four products compete for the same glycine inventory, not four independent pools.
    if sum(yields.values()) > seq.count('G') + 1e-9:
        raise ValueError('Shared Gly residue budget exceeded')
    return {p: rate*yields[p] for p in TARGETS}, rate*seq.count('G')


def check_parameters(spec):
    if spec.get('schema_version') != 1 or type(spec.get('enabled')) is not bool:
        raise ValueError('Explicit schema version and boolean enabled required')
    missing = []
    external = spec['external']
    for peptide in TARGETS:
        entry = external['peptides'][peptide]
        for field in ('medium_concentration_mM', 'uptake_capacity_mmol_gdw_h', 'uptake_evidence',
                      'hydrolysis_compartment', 'transport_steps', 'hydrolase_evidence', 'condition'):
            value = entry.get(field)
            if value is None or value == '' or value == []:
                missing.append('external.' + peptide + '.' + field)
            elif field in ('medium_concentration_mM', 'uptake_capacity_mmol_gdw_h'):
                number(value, field)
    endogenous = spec['endogenous']
    supply, shared_gly = Counter(), 0.
    if not endogenous['precursors']:
        missing.append('endogenous.precursors_with_source_sequence_rates_and_costs')
    seen_sources = set()
    for precursor in endogenous['precursors']:
        source_key = json.dumps(precursor['source_identity'], sort_keys=True)
        if source_key in seen_sources:
            raise ValueError('Duplicate precursor source can double count the same physical pool')
        seen_sources.add(source_key)
        produced, gly = precursor_supply(precursor)
        supply.update(produced)
        shared_gly += gly
    requested = endogenous.get('target_supply_mmol_gdw_h')
    if requested is None:
        missing.append('endogenous.target_supply_mmol_gdw_h')
    else:
        if set(requested) != set(TARGETS):
            raise ValueError('Supply must specify all four competing target products')
        for peptide, value in requested.items():
            number(value, 'target supply ' + peptide)
            if value > supply[peptide] + 1e-9:
                raise ValueError('Target supply exceeds evidenced cleavage/rate budget: ' + peptide)
        if sum(requested.values()) > shared_gly + 1e-9:
            raise ValueError('Total target supply exceeds shared Gly budget')
    if spec['enabled'] or external['enabled'] or endogenous['enabled']:
        raise ValueError('Activation refused: ' + ('missing ' + ', '.join(missing) if missing else
            'this is a specification checker, not an implemented or scientifically accepted model patch'))
    return {'status': 'disabled', 'activation_implemented': False, 'missing_parameters': missing,
            'proposed_supply_capacity': dict(supply) if endogenous['precursors'] else None,
            'shared_gly_capacity': shared_gly if endogenous['precursors'] else None,
            'checks_are': 'parameter consistency only; no biological evidence validation or model change'}


def self_test(path):
    spec = json.loads(path.read_text())
    initial = check_parameters(spec)
    assert initial['status'] == 'disabled'
    assert initial['shared_gly_capacity'] is None and initial['proposed_supply_capacity'] is None
    enabled = copy.deepcopy(spec)
    enabled['enabled'] = True
    try:
        check_parameters(enabled)
    except ValueError as exc:
        assert 'Activation refused' in str(exc)
    else:
        raise AssertionError('Incomplete activation accepted')
    # Synthetic fixtures exercise accounting only; these are not biological precursors or estimates.
    seq = 'GPGD'
    p = {'mature_sequence': seq, 'sequence_sha256': hashlib.sha256(seq.encode()).hexdigest(),
         'abundance_units': 'mmol_chain/gDW', 'rate_units': '1/h', 'abundance': 1.,
         'degradation_rate': 1., 'vacuolar_fraction': 1.,
         'replacement_synthesis_energy_per_chain': {'ATP_equivalent': 1.},
         **{k: 'synthetic test fixture, not scientific evidence' for k in
            ('source_identity', 'condition', 'synthesis_path', 'synthesis_cost_source',
             'localization_evidence', 'cleavage_evidence', 'balanced_chemical_equation', 'replacement_accounting')},
         'cleavage_patterns': [{'fraction': 1., 'target_pairs': [
             {'start_zero_based': 0, 'peptide': 'GP'}, {'start_zero_based': 2, 'peptide': 'GD'}],
             'remaining_products': []}]}
    assert precursor_supply(p) == ({'GD': 1., 'GE': 0., 'AG': 0., 'GP': 1.}, 2.)
    bad = copy.deepcopy(p); bad['cleavage_patterns'][0]['target_pairs'].append({'start_zero_based': 0, 'peptide': 'GP'})
    bad_units = copy.deepcopy(p); bad_units['abundance_units'] = 'TPM'
    bad_cost = copy.deepcopy(p); bad_cost['replacement_synthesis_energy_per_chain'] = {'ATP': 0}
    bad_residual = copy.deepcopy(p); bad_residual['cleavage_patterns'][0]['target_pairs'].pop()
    for invalid in (bad, bad_units, bad_cost, bad_residual):
        try:
            precursor_supply(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid precursor accounting accepted')
    excessive = copy.deepcopy(spec)
    excessive['endogenous']['precursors'] = [p]
    excessive['endogenous']['target_supply_mmol_gdw_h'] = {'GD': 2., 'GE': 0., 'AG': 0., 'GP': 2.}
    try:
        check_parameters(excessive)
    except ValueError:
        pass
    else:
        raise AssertionError('Over-budget supply accepted')
    duplicate = copy.deepcopy(spec)
    duplicate['endogenous']['precursors'] = [p, copy.deepcopy(p)]
    try:
        check_parameters(duplicate)
    except ValueError as exc:
        assert 'Duplicate precursor' in str(exc)
    else:
        raise AssertionError('Duplicate physical precursor source accepted')
    print('9 parameter checks passed; synthetic fixtures only; optimization calls=0; model edits=0')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('spec', type=Path)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        self_test(args.spec)
    else:
        print(json.dumps(check_parameters(json.loads(args.spec.read_text())), indent=2))
