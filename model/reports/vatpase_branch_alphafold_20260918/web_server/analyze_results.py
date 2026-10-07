"""Summarize the available A results; assertions are the runnable input checks."""
from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
manifest = json.loads((BASE.parent / 'sequence_manifest.json').read_text())
records = [r for r in manifest['records'] if r['input_group'] == 'A']
assert len(records) == 12
expected = {r['chain_id']: r['sequence_length'] for r in records}
rows = []
for index in range(5):
    prefix = BASE / 'results/A/fold_vatpase_a_12chains_4106aa_onecopy_20260918'
    summary_path = Path(f'{prefix}_summary_confidences_{index}.json')
    full_path = Path(f'{prefix}_full_data_{index}.json')
    s = json.loads(summary_path.read_text())
    f = json.loads(full_path.read_text())
    chain_ids = list(dict.fromkeys(s['chain_ids']))
    assert chain_ids == list(expected)
    assert Counter(f['token_chain_ids']) == Counter(expected)
    pae = np.asarray(f['pae'], dtype=np.float32)
    pair = np.asarray(s['chain_pair_iptm'])
    plddt = np.asarray(f['atom_plddts'])
    assert pae.shape == (4106, 4106) and pair.shape == (12, 12)
    assert len(f['atom_chain_ids']) == len(plddt)
    assert np.isfinite(pae).all() and np.isfinite(plddt).all()
    assert np.isfinite(pair).all() and ((pair >= 0) & (pair <= 1)).all()
    assert ((plddt >= 0) & (plddt <= 100)).all() and (pae >= 0).all()
    tokens = np.asarray(f['token_chain_ids'])
    atoms = np.asarray(f['atom_chain_ids'])
    row = {k: s[k] for k in ['iptm', 'ptm', 'ranking_score', 'has_clash', 'fraction_disordered']}
    row.update(index=index, chain_ids=chain_ids,
               chain_pair_iptm=s['chain_pair_iptm'],
               chain_pair_pae_min=s['chain_pair_pae_min'], chains=[],
               files={str(p.relative_to(BASE)): sha256(p.read_bytes()).hexdigest()
                      for p in [summary_path, full_path]})
    for i, r in enumerate(records):
        values = plddt[atoms == r['chain_id']]
        row['chains'].append(dict(chain=r['chain_id'], gene=r['gene_id'],
            residues=r['sequence_length'], atoms=len(values),
            plddt_atom_mean=float(values.mean()), plddt_atom_below_50=float((values < 50).mean()),
            chain_ptm=s['chain_ptm'][i], chain_iptm=s['chain_iptm'][i]))
    # Directed PAE block means describe placement uncertainty, not binding or AND/OR.
    row['pae_block_mean_angstrom'] = [[float(pae[np.ix_(tokens == a, tokens == b)].mean())
                                    for b in chain_ids] for a in chain_ids]
    rows.append(row)
output = dict(analyzed_utc=datetime.now(timezone.utc).isoformat(),
              tool='AlphaFold Server; exact backend patch version unknown', seed=1,
              scope='A: 5 summary and 5 full_data files; B raw files unavailable',
              plddt_weighting='per atom, not per residue or CA',
              code_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), results=rows)
(BASE / 'confidence_analysis.json').write_text(json.dumps(output, indent=2) + '\n')
print(json.dumps([{k: r[k] for k in ['index', 'iptm', 'ptm', 'has_clash']} for r in rows], indent=2))
