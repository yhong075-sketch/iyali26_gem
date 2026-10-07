"""Read-only local conservation checks explaining the saved candidate screen. No optimizer."""
import hashlib
import json
import math
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

OUT = Path(__file__).resolve().parent
NS = {'s': 'http://www.sbml.org/sbml/level3/version1/core'}
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
sources = [OUT/'candidate.xml', OUT/'candidate/run_manifest.json', OUT/'candidate/target_fluxes.json']
hashes = {str(p): sha(p) for p in sources}
manifest = json.loads(sources[1].read_text())
wt = json.loads(sources[2].read_text())['WT']
assert sha(sources[0]) == manifest['model']['sha256']
bounds, flux = manifest['base_bounds'], wt['fluxes']
xml = ET.parse(sources[0])
decode = lambda s: s.removeprefix('M_').replace('__91__', '[').replace('__93__', ']')
species = {decode(e.get('id')): e.attrib for e in xml.findall('.//s:species', NS)}
rx = {}
for r in xml.findall('.//s:reaction', NS):
    rid = r.get('id').removeprefix('R_')
    sto = {}
    for part, sign in [('listOfReactants', -1), ('listOfProducts', 1)]:
        for t in r.findall('s:'+part+'/s:speciesReference', NS):
            m = decode(t.get('species'))
            sto[m] = sto.get(m, 0) + sign*float(t.get('stoichiometry', 1))
    rx[rid] = dict(name=r.get('name'), stoichiometry=sto, bounds=bounds[rid], WT_flux=flux[rid])

def row(m):
    assert species[m]['boundaryCondition'] == 'false'
    return {r: d['stoichiometry'][m] for r, d in rx.items() if m in d['stoichiometry']}

go_h = row('m293[C_go]')
assert go_h == {'R141': 1, 'R794': 2, 'R1159': 1, 'R2048': -2, 'R_NDP7g': 1}
va_h = row('m1007[C_va]')
va_open = [r for r in va_h if bounds[r] != [0, 0]]
assert set(va_open) == {'R882', 'R2020', 'R2028', 'R2030', 'R2038', 'R2040'}
assert len(va_h) == 22
water = row('m1384[C_va]')
hyd = ['R2021', 'R2029', 'R2034', 'R2039']
assert water == {'R1363': 1, **{r: -1 for r in hyd}}
assert bounds['R1363'] == [0, 0] and all(bounds[r] == [0, 1000] for r in hyd)
assert row('m289[C_go]') == {'R141': -1}
assert row('m291[C_go]') == {'R141': 1}
assert bounds['R795'] == [0, 0]
assert not any(species[m]['compartment'] in ('C_go','C_va') for m in rx['biomass_C']['stoichiometry'])
assert rx['biomass_C']['stoichiometry']['m141[C_cy]'] == -23.09

# Sum all three dipeptide compartment balances: -exchange -hydrolysis = 0.
chains = []
for ids in [('R2018','R2019','R2020','R2021'),('R2026','R2027','R2028','R2029'),('R2031','R2032','R2033','R2034'),('R2036','R2037','R2038','R2039')]:
    ex, t1, t2, h = ids
    external = next(iter(rx[ex]['stoichiometry']))
    assert rx[ex]['stoichiometry'] == {external: -1} and bounds[ex] == [0, 1000]
    cy = next(m for m in rx[t1]['stoichiometry'] if species[m]['compartment']=='C_cy' and not species[m]['name'].startswith('H+'))
    va = next(m for m in rx[h]['stoichiometry'] if rx[h]['stoichiometry'][m]<0 and m!='m1384[C_va]')
    assert row(external) == {ex: -1, t1: -1}
    assert row(cy) == {t1: 1, t2: -1}
    assert row(va) == {t2: 1, h: -1}
    chains.append(dict(reactions=ids, species=[external,cy,va], summed_row={ex:-1,h:-1}, conclusion='All four fluxes must be zero under current bounds'))

golgi = [r for r,d in rx.items() if any(species[m]['compartment']=='C_go' for m in d['stoichiometry'])]
assert len(golgi)==76 and all(flux[r]==0 for r in golgi)
assert all(flux[r]==0 for r in va_h)
mannan_demand = -rx['biomass_C']['stoichiometry']['m1324[C_cy]']*wt['growth']
assert row('m504[C_er]') == {'R289':1, 'R1250':1}
assert rx['R1250']['stoichiometry'] == {'m1324[C_cy]':-1, 'm504[C_er]':1}
assert math.isclose(flux['R289'],mannan_demand,abs_tol=1e-12)
assert math.isclose(-flux['R1250'],mannan_demand,abs_tol=1e-12)
assert all(sha(Path(p))==s for p,s in hashes.items())
selected=set(go_h)|set(va_h)|set(water)|{'R289','R1250','R642','R694','biomass_C'}|{r for c in chains for r in c['reactions']}
out=dict(checked_utc=datetime.now(timezone.utc).isoformat(),scope='Candidate screen only; local stoichiometric rows and saved WT. No LP/FVA or model changes.',
 source_sha256=hashes,script_sha256=sha(Path(__file__)),solver_calls=0,WT_growth=wt['growth'],
 golgi_proton_row=go_h,vacuolar_proton_row=va_h,vacuolar_nonzero_bounds=va_open,
 vacuolar_water_row=water,dipeptide_supply_chains=chains,
 R141_isolated_rows={m:row(m)for m in ('m289[C_go]','m291[C_go]')},
 golgi_compartment_reactions_saved_WT_zero=golgi,mannan_WT_demand=mannan_demand,
 reactions={r:rx[r]for r in sorted(selected)},species={m:species[m]for r in selected for m in rx[r]['stoichiometry']},
 limits=['76 Golgi-associated reactions are zero in this saved WT; not proof all are structurally blocked.', 'Biological validity of missing supply and native acidification requirements is unresolved.', 'Saved ATP-producing flux is not an independently validated physiological energy model.'])
(OUT/'mechanism_details.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(status='passed',golgi_WT_zero=len(golgi),vacuolar_H_rows=len(va_h),water_blocked_hydrolyses=hyd,mannan_demand=mannan_demand,solver_calls=0)))
