"""Read the fixed local model; no optimization or model writes."""
import csv
import datetime
import hashlib
import json
from pathlib import Path
import cobra

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
MODEL = ROOT / 'model_metadata_trna_r539_alphafold_labeled.xml'
IDS = ['R763', 'R407', 'R969', 'R39', 'R808', 'R715', 'R40', 'R19', 'R18', 'R695', 'R385']
GENES = ['YALI1C26017g', 'YALI1F08349g', 'YALI1B20835g', 'YALI1F34625g', 'YALI1C25352g', 'YALI1A08781g', 'YALI1E18269g', 'YALI1B20527g', 'YALI1F34675g']


def main():
    source_sha = hashlib.sha256(MODEL.read_bytes()).hexdigest()
    assert source_sha == json.loads((BASE / 'scope.json').read_text())['model_sha256']
    model = cobra.io.read_sbml_model(str(MODEL))
    def record(r):
        return {'id': r.id, 'name': r.name, 'gpr': r.gene_reaction_rule,
                'bounds': list(r.bounds), 'equation': r.reaction,
                'annotation': r.annotation, 'notes': r.notes,
                'stoichiometry': {m.id: c for m, c in r.metabolites.items()},
                'species': {m.id: {'name': m.name, 'formula': m.formula,
                            'charge': m.charge, 'compartment': m.compartment}
                            for m in r.metabolites},
                'balance_residual': r.check_mass_balance()}
    reactions = [record(model.reactions.get_by_id(r)) for r in IDS]
    original = {r.id: r.bounds for r in model.reactions}
    genes = []
    for n, gid in enumerate(GENES, 1):
        g = model.genes.get_by_id(gid)
        initial_functional = g.functional
        with model:
            g.knock_out()
            changed = sorted(r.id for r in model.reactions if r.bounds != original[r.id])
        # COBRA orphan genes can lack a model context; restore that software flag explicitly.
        g.functional = initial_functional
        assert all(r.bounds == original[r.id] for r in model.reactions)
        assert g.functional == initial_functional
        genes.append({'gene_id': gid, 'family_candidate': 'COQ' + str(n),
                      'model_name': g.name, 'model_annotation': g.annotation,
                      'reactions': sorted(r.id for r in g.reactions), 'static_KO_closes': changed})
    duplicates = []
    for rid in IDS:
        r = model.reactions.get_by_id(rid)
        stoich = r.metabolites
        for other in model.reactions:
            if other.id == rid:
                continue
            if other.metabolites == stoich or other.metabolites == {m: -c for m, c in stoich.items()}:
                duplicates.append({'target': rid, 'other': record(other)})
    q_ids = ['m468[C_mi]', 'm471[C_mi]']
    q_coefficients = {r.id: sum(c for m, c in r.metabolites.items() if m.id in q_ids)
                      for r in model.reactions}
    q_coefficients = {r: c for r, c in q_coefficients.items() if c != 0}
    neighbors = sorted({r.id for mid in ['m984[C_mi]', 'm985[C_mi]', 'm138[C_mi]']
                        for r in model.metabolites.get_by_id(mid).reactions})
    assert len(genes) == 9 and len(reactions) == 11
    assert q_coefficients == {'R385': 1.0}
    assert all(not r['balance_residual'] for r in reactions)
    result = {'checked_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'model': str(MODEL), 'model_sha256': source_sha, 'cobra_version': cobra.__version__,
              'optimization_calls': 0, 'counts': [len(model.reactions), len(model.metabolites), len(model.genes)],
              'reactions': reactions, 'genes': genes, 'exact_signed_duplicates': duplicates,
              'direct_precursor_neighbors': [record(model.reactions.get_by_id(r)) for r in neighbors],
              'q9_plus_q9h2_nonzero_row_sum': q_coefficients,
              'q9_pool_inference': 'For this SBML steady-state network, sum of Q9 and Q9H2 balance rows gives v_R385 = 0. No pool dilution/turnover demand is encoded by those rows. This is not a new WT/KO growth solve.'}
    (BASE / 'current_model_extract.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    with (BASE / 'current_reactions.tsv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=['id', 'name', 'gpr', 'bounds', 'equation'], delimiter='\t')
        w.writeheader()
        w.writerows({k: r[k] for k in w.fieldnames} for r in reactions)
    print(json.dumps({'counts': result['counts'], 'genes': genes, 'duplicates': duplicates,
                      'precursor_reactions': neighbors, 'q_pool_sum': q_coefficients}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
