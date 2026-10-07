"""Read-only compartment incidence screen. Stdlib only; no LP/FVA or model edits."""
import collections, datetime, hashlib, json, re, time
from pathlib import Path
from xml.etree import ElementTree as ET
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
D=ROOT/'artifacts/iyli647_screen_20260910/nonessential_diagnosis_20260911'
NS={'s':'http://www.sbml.org/sbml/level3/version1/core','f':'http://www.sbml.org/sbml/level3/version1/fbc/version2'}
F='{'+NS['f']+'}'
RDF='{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(s,prefix): return re.sub(r'__(\d+)__',lambda m:chr(int(m[1])),s.removeprefix(prefix))
def xml_model(path):
    root=ET.parse(path).getroot()
    params={e.get('id'):float(e.get('value')) for e in root.findall('.//s:parameter',NS)}
    mets={}
    for e in root.findall('.//s:species',NS):
        mets[dec(e.get('id'),'M_')]={'name':e.get('name'),'compartment':e.get('compartment'),'formula':e.get(F+'chemicalFormula'),'charge':e.get(F+'charge'),'boundary':e.get('boundaryCondition')=='true','annotations':sorted({x.get(RDF+'resource') for x in e.iter() if x.get(RDF+'resource')}),'notes':[''.join(x.itertext()) for x in e.findall('.//{http://www.w3.org/1999/xhtml}p')]}
    reactions={}
    for e in root.findall('.//s:reaction',NS):
        st=collections.defaultdict(float)
        for tag,sign in [('listOfReactants',-1),('listOfProducts',1)]:
            for x in e.findall('s:'+tag+'/s:speciesReference',NS):st[dec(x.get('species'),'M_')]+=sign*float(x.get('stoichiometry','1'))
        reactions[dec(e.get('id'),'R_')]={'stoichiometry':{m:c for m,c in st.items() if c},'bounds':[params[e.get(F+'lowerFluxBound')],params[e.get(F+'upperFluxBound')]],'notes':[''.join(x.itertext()) for x in e.findall('.//{http://www.w3.org/1999/xhtml}p')],'annotations':sorted({x.get(RDF+'resource') for x in e.iter() if x.get(RDF+'resource')})}
    return mets,reactions,{e.get('id'):e.get('name') for e in root.findall('.//s:compartment',NS)}
def incidence(rx):
    rows=collections.defaultdict(dict)
    for r,d in rx.items():
        for m,c in d['stoichiometry'].items(): rows[m][r]=c
    return dict(rows)
def eliminate(rx, boundary=()):
    rows=incidence(rx); zeros={r:{'round':0,'reason':'bounds_zero'} for r,d in rx.items() if d['bounds']==[0,0]};certs=[]
    # ponytail: sign-row elimination is a sufficient test; add LP/FVA only in a separately authorized analysis if completeness matters.
    for step in range(1,len(rx)+2):
        new={}
        for m,row in rows.items():
            if m in boundary:continue
            live={r:c for r,c in row.items() if r not in zeros}
            if not live:continue
            signs=set()
            for r,c in live.items():
                lo,hi=rx[r]['bounds']
                if hi>0:signs.add(1 if c>0 else -1)
                if lo<0:signs.add(-1 if c>0 else 1)
            if len(live)==1 or len(signs)==1:
                reason='single_remaining_reaction' if len(live)==1 else 'one_sided_contributions'
                cert={'round':step,'metabolite':m,'reason':reason,'full_row':row,'prior_zero_reactions':sorted(set(row)&zeros.keys()),'forced_zero':sorted(live)}
                certs.append(cert)
                for r in live:
                    lo,hi=rx[r]['bounds']
                    assert lo<=0<=hi,('infeasible nonzero bound',r,m,lo,hi)
                    new.setdefault(r,{'round':step,'reason':reason,'metabolite':m,'certificate_index':len(certs)-1})
        if not new:break
        zeros.update(new)
    else:raise AssertionError('fixed point not reached')
    return zeros,certs

def self_check():
    def r(st,b):return {'stoichiometry':st,'bounds':b}
    toy={'rev':r({'a':-1,'b':1},[-1,1]),'sink':r({'b':-1},[0,1]),'closed':r({'a':1},[0,0])}
    z,c=eliminate(toy);assert set(z)==set(toy) and z['rev']['reason']=='single_remaining_reaction'
    cyc={'f':r({'a':-1,'b':1},[-1,1]),'g':r({'a':1,'b':-1},[0,1])}
    assert not eliminate(cyc)[0], 'opposing contributions cannot be removed'
    try:eliminate({'bad':r({'a':1},[1,2])})
    except AssertionError:pass
    else:raise AssertionError('must flag infeasible nonzero lower bound')
    assert not eliminate({'external':r({'a':1},[0,1])},boundary={'a'})[0]

def identity(a,b):
    shared=sorted(set(a['annotations'])&set(b['annotations']))
    strong=[x for x in shared if any('/'+ns+'/' in x for ns in ['metanetx.chemical','inchikey','kegg.compound'])]
    chemistry_ok=a['formula'] and a['formula']==b['formula'] and a['charge'] is not None and a['charge']==b['charge']
    if strong and chemistry_ok:return {'status':'annotation_and_formula_charge_agree','shared':strong}
    if strong:return {'status':'identity_ambiguous_formula_charge','shared':strong}
    return None

def main():
    started=time.monotonic();self_check()
    files={'model':ROOT/'model_metadata_trna.xml','ancestral_iyali21':ROOT/'data/iyli21.xml','effective_snapshot':D/'reaction_snapshot.json','saved_main':D/'results.json','saved_controls':D/'controls.json'}
    inputs={k:{'path':str(p),'sha256':sha(p)} for k,p in files.items()}
    assert inputs['model']['sha256']=='d274bad3050e3c9220a8b6287eae847f3bf1334892284d565a6c4d96b38135a0'
    mets,raw,comps=xml_model(files['model']);_,old,_=xml_model(files['ancestral_iyali21'])
    rx=json.loads(files['effective_snapshot'].read_text());assert len(rx)==2315 and len(mets)==1877 and set(rx)==set(raw)
    assert all(d['stoichiometry']==raw[r]['stoichiometry'] for r,d in rx.items()),'snapshot chemistry drift'
    rows=incidence(rx);zero,certs=eliminate(rx,{m for m,d in mets.items() if d['boundary']})
    witnesses={}
    for key in ['saved_main','saved_controls']:
        data=json.loads(files[key].read_text())
        for label,run in data['runs'].items():
            assert run['raw_status']=='optimal'
            witnesses[key+':'+label]=run['fluxes']
    assert len(witnesses)==22
    maxzero=max((abs(v[r]) for v in witnesses.values() for r in zero),default=0)
    assert maxzero<1e-7,('certificate contradicts saved witness',maxzero)
    rxinfo={};transports=[];ambiguous=[]
    for r,d in rx.items():
        comp=sorted({mets[m]['compartment'] for m in d['stoichiometry']})
        vals={k:v[r] for k,v in witnesses.items()}
        rxinfo[r]={'name':d['name'],'bounds_effective':d['bounds'],'bounds_in_file':raw[r]['bounds'],'bounds_iYali21':old.get(r,{}).get('bounds'),'same_stoichiometry_iYali21':old.get(r,{}).get('stoichiometry')==d['stoichiometry'],'gpr':d['gpr'],'stoichiometry':d['stoichiometry'],'compartments':comp,'saved_flux_min':min(vals.values()),'saved_flux_max':max(vals.values()),'saved_flux_WT':vals['saved_main:WT'],'saved_fluxes':vals,'forced_zero':zero.get(r),'notes':raw[r]['notes'],'annotations':raw[r]['annotations']}
        if len(comp)<2:continue
        pairs=[];amb=[]
        for a,ca in d['stoichiometry'].items():
            if ca>=0:continue
            for b,cb in d['stoichiometry'].items():
                if cb<=0 or mets[a]['compartment']==mets[b]['compartment']:continue
                ident=identity(mets[a],mets[b])
                if ident:
                    p={'from_species':a,'to_species':b,'from_coefficient':ca,'to_coefficient':cb,**ident}
                    (pairs if ident['status']=='annotation_and_formula_charge_agree' else amb).append(p)
        # The all-paired count is intentionally conservative: equal opposite coefficients and each species paired exactly once.
        counts=collections.Counter(x for p in pairs if p['from_coefficient']==-p['to_coefficient'] for x in [p['from_species'],p['to_species']])
        pure=set(counts)==set(d['stoichiometry']) and all(v==1 for v in counts.values())
        row={'reaction':r,'class':'pure_balanced_transport' if pure else 'cross_compartment_with_transport_pair' if pairs else 'cross_compartment_no_confirmed_pair','confirmed_pairs':pairs,'ambiguous_pairs':amb}
        transports.append(row)
        if amb:ambiguous.append(r)
    metinfo={}
    for m,d in mets.items():
        links=[]
        for r,c in rows.get(m,{}).items():
            lo,hi=rx[r]['bounds'];produces=[];consumes=[]
            for direction,allowed,sign in [('forward',hi>0,1),('reverse',lo<0,-1)]:
                if allowed:(produces if c*sign>0 else consumes).append(direction)
            links.append({'reaction':r,'coefficient':c,'producer_directions':produces,'consumer_directions':consumes,'forced_zero':r in zero})
        live=[x for x in links if not x['forced_zero']]
        metinfo[m]={**d,'all_connections':links,'producer_reactions_after_elimination':[x['reaction'] for x in live if x['producer_directions']],'consumer_reactions_after_elimination':[x['reaction'] for x in live if x['consumer_directions']],'distinct_reactions_after_elimination':len(live)}
    counts=collections.Counter(t['class'] for t in transports)
    counts.update({'metabolites':len(mets),'reactions':len(rx),'compartments':len(comps),'cross_compartment':len(transports),'initial_closed':sum(v['round']==0 for v in zero.values()),'additional_forced_zero':sum(v['round']>0 for v in zero.values()),'all_forced_zero':len(zero),'elimination_rounds':max(v['round'] for v in zero.values()),'identity_ambiguous_reactions':len(ambiguous),'saved_witnesses':len(witnesses)})
    counts.update({'cross_compartment_forced_zero':sum(t['reaction'] in zero for t in transports),'pure_transport_no_gpr':sum(not rx[t['reaction']]['gpr'] for t in transports if t['class']=='pure_balanced_transport'),'pure_transport_same_stoichiometry_and_file_bounds_iYali21':sum(rxinfo[t['reaction']]['same_stoichiometry_iYali21'] and rxinfo[t['reaction']]['bounds_iYali21']==rxinfo[t['reaction']]['bounds_in_file'] for t in transports if t['class']=='pure_balanced_transport')})
    important=['R1264','R1265','R1590','R1591','R1593','R1594','R1628','R1629','R1022','R1630','R1638','R2077','R2084','R2118','R2167','R1363','R795','R87','R539','R2233','R707','R709','R711','R1914','R1460','R1461','R1462','R1463']
    out={'checked_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'All 2315 reactions and 1877 compartment-specific species; saved SD-Leu/PO1f runtime bounds. Static sign-row sufficient zero certificates, transport annotation/chemistry screening; zero LP/FVA. No model/chemistry/GPR/media/labels changed.','budget':'600 seconds wall; stdlib only; stop on hash/stoichiometry drift, contradictory nonzero bound, saved-witness conflict, or budget.','inputs':inputs,'code_sha256':sha(Path(__file__)),'counts':dict(counts),'compartment_names':comps,'limitations':['Single-row sign elimination is sufficient, not complete blocked-reaction analysis.','A reversible reaction is not counted as an independent source of its own substrate.','Metabolites remain separate model species; shared annotation and equal formula/charge mark an intra-reaction correspondence, not a global identity merge or biological transport proof.','Missing GPR or transport citations are evidence pending, not proof transport does not exist.','ER membrane and lipid particle may represent membrane-face bookkeeping; do not equate every compartment edge with a physical membrane crossing.','Only saved 22 witnesses checked; zero in those witnesses does not establish all-condition blockage unless a certificate is supplied.'],'max_abs_saved_flux_on_forced_zero':maxzero,'optimization_calls':0,'certificates':certs,'forced_zero':zero,'transport_classification':transports,'metabolites':metinfo,'reactions':rxinfo,'priority_reaction_ids':[r for r in important if r in rx]}
    assert all(sha(files[k])==v['sha256'] for k,v in inputs.items())
    out['inputs_unchanged']=True;out['elapsed_seconds']=time.monotonic()-started
    assert out['elapsed_seconds']<600,'budget exceeded'
    (OUT/'compartment_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'counts':out['counts'],'max_abs_saved_flux_on_forced_zero':maxzero,'elapsed_seconds':out['elapsed_seconds']},ensure_ascii=False))
if __name__=='__main__':main()
