import sys, re
import xml.etree.ElementTree as ET
p=sys.argv[1]
ns={'s':'http://www.sbml.org/sbml/level3/version1/core','fbc':'http://www.sbml.org/sbml/level3/version1/fbc/version2','rdf':'http://www.w3.org/1999/02/22-rdf-syntax-ns#'}
t=ET.parse(p); m=t.getroot().find('s:model',ns)
genes=['YALI1D04261g','YALI1F39330g','YALI1B05072g','YALI1F02309g','YALI1F12969g','YALI1D35735g']
for gp in m.find('fbc:listOfGeneProducts',ns):
    lab=gp.get('{%s}label'%ns['fbc'])
    if lab in genes:
        res=[li.get('{%s}resource'%ns['rdf']) for li in gp.iter('{%s}li'%ns['rdf'])]
        print('GENE',lab,res)
species={}
for sp in m.find('s:listOfSpecies',ns):
    species[sp.get('id')]=(sp.get('name'),sp.get('compartment'))
def gpr_str(e):
    tag=e.tag.split('}')[1]
    if tag=='geneProductRef': return e.get('{%s}geneProduct'%ns['fbc'])[2:]
    sub=[gpr_str(c) for c in e]
    return '('+(' %s '%tag).join(sub)+')'
params={pp.get('id'):pp.get('value') for pp in m.find('s:listOfParameters',ns)}
for r in m.find('s:listOfReactions',ns):
    ga=r.find('fbc:geneProductAssociation',ns)
    g=gpr_str(list(ga)[0]) if ga is not None else ''
    rs=[(x.get('species'),x.get('stoichiometry')) for x in (r.find('s:listOfReactants',ns) or [])]
    ps=[(x.get('species'),x.get('stoichiometry')) for x in (r.find('s:listOfProducts',ns) or [])]
    names=[species[s][0] for s,_ in rs+ps]
    txt=' '.join(n or '' for n in names).lower()
    if any(x in g for x in genes) or 'thioredoxin' in txt or 'trx' in txt:
        lb=params.get(r.get('{%s}lowerFluxBound'%ns['fbc'])); ub=params.get(r.get('{%s}upperFluxBound'%ns['fbc']))
        f=lambda L:' + '.join('%s %s[%s](%s)'%(st,species[s][0],species[s][1],s) for s,st in L)
        print(r.get('id'),'|',r.get('name'),'|',lb,ub,'|',f(rs),'->',f(ps),'| GPR:',g)
