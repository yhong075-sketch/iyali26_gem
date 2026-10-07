"""Static GPR and identical-chemistry audit; intentionally never optimize."""
import ast
import csv
import hashlib
import inspect
import json
import subprocess
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
import cobra
from cobra.core.gene import Gene

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DIAG = ROOT / 'artifacts/iyli647_screen_20260910/nonessential_diagnosis_20260911'
SCREEN = ROOT / 'artifacts/screen_test_metadata_trna_20260910'
MODEL = ROOT / 'model_metadata_trna.xml'
EXPECTED = 'd274bad3050e3c9220a8b6287eae847f3bf1334892284d565a6c4d96b38135a0'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

def parsed(rule):
    if not rule:
        return None, set()
    tree = ast.parse(rule, mode='eval').body
    def check(node):
        if isinstance(node, ast.Name):
            return {node.id}
        if isinstance(node, ast.BoolOp) and isinstance(node.op, (ast.And, ast.Or)):
            return set().union(*(check(v) for v in node.values))
        raise ValueError(f'Unsupported GPR AST: {ast.dump(node)}')
    return tree, check(tree)

def active(tree, knocked):
    if tree is None:
        return True
    if isinstance(tree, ast.Name):
        return tree.id not in knocked
    values = [active(v, knocked) for v in tree.values]
    return all(values) if isinstance(tree.op, ast.And) else any(values)

def signature(stoich):
    ordered = sorted((m, Fraction(str(c))) for m, c in stoich.items() if c)
    if not ordered:
        return (), Fraction(0)
    scale = ordered[0][1]
    return tuple((m, c / scale) for m, c in ordered), scale

def directions(bounds, scale=1):
    lo, hi = sorted(float(scale) * b for b in bounds)
    return ({'+'} if hi > 0 else set()) | ({'-'} if lo < 0 else set())

def table(path):
    with path.open() as f:
        return {r['gene_id']: r for r in csv.DictReader(f, delimiter='\t')}

def self_check():
    t, ids = parsed('(A and B) or C')
    assert ids == {'A', 'B', 'C'}
    assert active(t, {'A'}) and not active(t, {'A', 'C'})
    assert active(None, {'A'})
    try:
        parsed('A + B')
    except ValueError:
        pass
    else:
        raise AssertionError('Invalid operator accepted')
    a, sa = signature({'a[c]': -2, 'b[c]': 4})
    b, sb = signature({'a[c]': 1, 'b[c]': -2})
    assert a == b and sb / sa == Fraction(-1, 2)
    assert directions([0, 1000], -2) == {'-'}
    assert signature({'a[m]': -2, 'b[m]': 4})[0] != a


def main():
    start = time.monotonic()
    self_check()
    inputs = [MODEL, DIAG/'reaction_snapshot.json', DIAG/'results.json', SCREEN/'run_manifest.json', SCREEN/'screen_predictions.tsv', SCREEN/'essentiality_per_gene.tsv']
    identities = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    assert sha(MODEL) == EXPECTED
    snapshot = json.loads((DIAG/'reaction_snapshot.json').read_text())
    saved = json.loads((DIAG/'results.json').read_text())
    manifest = json.loads((SCREEN/'run_manifest.json').read_text())
    predictions = table(SCREEN/'screen_predictions.tsv')
    calls = table(SCREEN/'essentiality_per_gene.tsv')
    assert saved['model']['sha256'] == EXPECTED == manifest['model']['sha256']
    assert manifest['simulation_context'] == saved['simulation_context']
    # Loading creates a solver object, but no optimization method is called.
    model = cobra.io.read_sbml_model(str(MODEL))
    native = sorted(g.id for g in model.genes)
    assert len(native) == 1074 and set(native) == set(manifest['gene_ids'])
    assert set(snapshot) == {r.id for r in model.reactions} and len(snapshot) == 2315
    overlays = []
    for r in model.reactions:
        s = snapshot[r.id]
        assert {m.id: float(c) for m, c in r.metabolites.items()} == s['stoichiometry']
        if list(r.bounds) != s['bounds'] or r.gene_reaction_rule != s['gpr']:
            overlays.append({'reaction':r.id, 'original_bounds':list(r.bounds), 'effective_bounds':s['bounds'], 'original_gpr':r.gene_reaction_rule, 'effective_gpr':s['gpr']})
        r.bounds = tuple(s['bounds'])
        r.gene_reaction_rule = s['gpr']
    trees = {}
    links = defaultdict(list)
    for rid, s in snapshot.items():
        tree, genes = parsed(s['gpr'])
        trees[rid] = tree
        assert genes == {g.id for g in model.reactions.get_by_id(rid).genes}
        for g in genes:
            links[g].append(rid)
    gene_rows = []
    for gid in native:
        assert time.monotonic() - start < 550
        associated = sorted(links[gid])
        false_rules = sorted(rid for rid in associated if not active(trees[rid], {gid}))
        expected = sorted(rid for rid in false_rules if snapshot[rid]['bounds'] != [0.,0.])
        with model:
            model.genes.get_by_id(gid).knock_out()
            actual = sorted(r.id for r in model.reactions if list(r.bounds) != snapshot[r.id]['bounds'])
            mismatches = [{'reaction':r.id,'actual':list(r.bounds),'expected':[0.,0.] if r.id in false_rules else snapshot[r.id]['bounds']} for r in model.reactions if list(r.bounds) != ([0.,0.] if r.id in false_rules else snapshot[r.id]['bounds'])]
        assert not mismatches and actual == expected, (gid, mismatches)
        open_associated = [rid for rid in associated if snapshot[rid]['bounds'] != [0.,0.]]
        category = ('new_reaction_closure' if actual else 'no_effective_association' if not associated else 'associated_reactions_all_already_closed' if not open_associated else 'all_open_associated_rules_remain_true')
        identity = calls.get(gid,{})
        gene_rows.append({'gene_id':gid, 'native_name':'not verified in this static audit', 'candidate_protein_function':identity.get('source_putative_function') or identity.get('model_associated_reactions') or 'not assigned in effective GPR', 'function_evidence':'Saved experimental-source putative annotation / model-GPR assignment only; no new biological validation', 'associated_reactions':associated,'rules_false_after_KO':false_rules,'newly_closed_reactions':actual,'protected_open_reactions':sorted(set(open_associated)-set(expected)),'classification':category,'api_mismatches':mismatches,'saved_screen':predictions.get(gid), 'saved_experimental_comparison':identity})
    by_gene = {r['gene_id']:r for r in gene_rows}
    buckets = defaultdict(list)
    scales = {}
    for rid,s in snapshot.items():
        sig, scale = signature(s['stoichiometry'])
        if sig:
            buckets[sig].append(rid)
            scales[rid] = scale
    duplicates = []
    for sig, members in buckets.items():
        if len(members) < 2:
            continue
        members.sort()
        mets = [model.metabolites.get_by_id(m) for m,c in sig]
        compartments = sorted({m.compartment for m in mets})
        names = ' '.join(snapshot[r]['name'].lower() + ' ' + r.lower() for r in members)
        kind = ('boundary' if len(sig)==1 else 'biomass_or_pool' if any(w in names for w in ['biomass','biomembrane','protein pool']) else 'multi_compartment' if len(compartments)>1 else 'single_compartment_chemistry')
        relevant = sorted(set().union(*(model.reactions.get_by_id(r).gpr.genes for r in members)) & set(native))
        bypasses = []
        for gid in relevant:
            closed = [r for r in members if r in by_gene[gid]['newly_closed_reactions']]
            if not closed:
                continue
            remain = [r for r in members if r not in closed and directions(snapshot[r]['bounds'])]
            removed_dirs = set().union(*(directions(snapshot[r]['bounds'],scales[r]) for r in closed))
            retained_dirs = set().union(*(directions(snapshot[r]['bounds'],scales[r]) for r in remain))
            overlap = removed_dirs & retained_dirs
            if overlap:
                bypasses.append({'gene_id':gid,'closed_members':closed,'surviving_members':remain,'overlapping_canonical_directions':sorted(overlap),'saved_WT_fluxes':{r:saved['runs']['WT']['fluxes'][r] for r in members}, 'saved_target_KO_fluxes':{r:saved['runs'][gid]['fluxes'][r] for r in members} if gid in saved['runs'] else None, 'saved_classification':calls.get(gid,{}).get('classification_at_10pct')})
        duplicates.append({'members':members,'kind':kind,'compartments':compartments,'exact_or_proportional':'exact_same' if len({scales[r] for r in members})==1 else 'sign_or_scale_normalized','distinct_gpr_text_count':len({snapshot[r]['gpr'] for r in members}), 'reaction_details':[{**snapshot[r],'id':r,'canonical_scale':str(scales[r]),'canonical_allowed_directions':sorted(directions(snapshot[r]['bounds'],scales[r]))} for r in members], 'single_KO_same_direction_survival':bypasses})
    duplicates.sort(key=lambda d:(-len(d['single_KO_same_direction_survival']),d['members']))
    counts = Counter(r['classification'] for r in gene_rows)
    fn_by_category = Counter(r['classification'] for r in gene_rows if r['saved_experimental_comparison'].get('classification_at_10pct')=='FN')
    result = {'checked_at_utc':datetime.now(timezone.utc).isoformat(), 'plan':{'scope':'All 1074 original-model gene records, including placeholders/mitochondrial identifiers, Boolean KOs and COBRA bounds checks; 2315 exact/proportional same-metabolite duplicate screen','budget_wall_seconds':600,'optimization_calls':0,'stop':'Input mismatch, parse error, API discrepancy, or 550 seconds','excluded':'No LP/FVA, no model/curation edits, no literature or sequence validation'}, 'input_sha256':identities,'script_sha256':sha(__file__),'gene_knock_out_source':{'path':inspect.getfile(Gene),'sha256':sha(inspect.getfile(Gene))},'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'dirty_state_observed':subprocess.check_output(['git','status','--porcelain=v1'],cwd=ROOT,text=True),'cobra_version':cobra.__version__,'runtime_context':saved['simulation_context'],'medium':saved['medium'],'strain_profile':saved['strain_profile'],'runtime_overlays_applied':overlays,'original_model_gene_ids':native,'runtime_only_gene_ids':sorted({g.id for g in model.genes}-set(native)),'summary':{'original_model_gene_records':len(native),'effective_gene_records':len(model.genes),'reactions':len(snapshot),'gene_classifications':dict(counts),'saved_FN_by_category':dict(fn_by_category),'api_bound_mismatches':0,'parse_errors':0,'reaction_GPR_nonempty':sum(bool(s['gpr']) for s in snapshot.values()),'duplicate_groups':len(duplicates),'duplicate_reactions':sum(len(d['members']) for d in duplicates),'duplicate_group_kinds':dict(Counter(d['kind'] for d in duplicates)),'duplicate_groups_with_same_direction_KO_survival':sum(bool(d['single_KO_same_direction_survival']) for d in duplicates),'duplicate_survival_gene_count':len({b['gene_id'] for d in duplicates for b in d['single_KO_same_direction_survival']}),'duplicate_survival_FN_gene_count':len({b['gene_id'] for d in duplicates for b in d['single_KO_same_direction_survival'] if b['saved_classification']=='FN'})},'gene_checks':gene_rows,'duplicate_groups':duplicates,'limitations':['Boolean/API consistency validates implementation of encoded GPR, not biological correctness.','Remaining allowed duplicate direction does not prove a feasible growth-supporting bypass; saved witnesses are marked separately.','Duplicate detection uses exact decimal stoichiometry and identical compartment-specific metabolite IDs; chemically equivalent different IDs or rounding differences are not merged.','Historical screen labels and predictions are reused, not newly reproduced; unlabelled is not negative.'],'inputs_unchanged':identities=={str(p.relative_to(ROOT)):sha(p) for p in inputs},'elapsed_seconds':time.monotonic()-start,'self_check':'passed'}
    assert result['inputs_unchanged']
    (OUT/'gpr_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result['summary'],ensure_ascii=False,indent=2))

if __name__ == '__main__':
    main()
