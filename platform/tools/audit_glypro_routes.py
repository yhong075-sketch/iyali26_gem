"""Offline, no-optimization audit of the fixed Gly-L-Pro compartment input."""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from cobra import Reaction
from cobra.io import read_sbml_model

from tools.diagnose_dipeptide_supply import balance, sha, table, write_json
from scripts.gem_annotate.dipeptide_chemistry import protected_chemistry
from scripts.gem_annotate.energy_candidates import protected_definitions


ROOT = Path(__file__).resolve().parents[2]
EXPECTED_SHA = 'a0307b9b00c1ed6e981d605f65ba15a0e313fe95e0e664dc05cfb95ef4291848'
TARGETS = {
    'm1864[C_ex]': 'Gly-L-Pro extracellular', 'm1865[C_cy]': 'Gly-L-Pro cytoplasmic',
    'm1866[C_va]': 'Gly-L-Pro vacuolar', 'm272[C_cy]': 'glycine cytoplasmic',
    'm1863[C_va]': 'glycine vacuolar', 'm765[C_cy]': 'L-proline cytoplasmic',
    'm1867[C_va]': 'L-proline vacuolar', 'm32[C_cy]': 'water cytoplasmic',
    'm1384[C_va]': 'water vacuolar', 'm10[C_cy]': 'proton cytoplasmic',
    'm1007[C_va]': 'proton vacuolar', 'm141[C_cy]': 'ATP cytoplasmic',
    'm143[C_cy]': 'ADP cytoplasmic', 'm35[C_cy]': 'phosphate cytoplasmic',
}
CORE = ['R2036', 'R2037', 'R2038', 'R2039', 'R2030', 'R2040', 'R1363', 'R795']
CY_HYD = {'m1865[C_cy]': -1, 'm32[C_cy]': -1, 'm272[C_cy]': 1, 'm765[C_cy]': 1}


from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.config import load_project_paths
def workspace(path):
    path = Path(path).resolve()
    # Authorized: the repository and the research workspace's task outputs.
    if not any(path.is_relative_to(root) for root in (ROOT, load_project_paths().task_outputs)):
        raise ValueError('All audit paths must be inside the authorized project')
    return path


def coefficients(reaction):
    return {m.id: c for m, c in sorted(reaction.metabolites.items(), key=lambda p: p[0].id)}


def equation(model, stoich):
    def term(mid, value):
        return f'{abs(value):g} {model.metabolites.get_by_id(mid).name} [{mid}]'
    return ' + '.join(term(k, v) for k, v in sorted(stoich.items()) if v < 0) + ' -> ' + ' + '.join(
        term(k, v) for k, v in sorted(stoich.items()) if v > 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--source-sha256', required=True)
    parser.add_argument('--output', required=True, help='New task artifact root; audit subdirectory must not exist')
    args = parser.parse_args()
    source, out = workspace(args.source), workspace(args.output)
    if sha(source) != args.source_sha256 or args.source_sha256 != EXPECTED_SHA:
        raise ValueError('Explicit fixed chemical input SHA mismatch')
    audit = out / 'audit'
    audit.mkdir(parents=True, exist_ok=False)
    model = read_sbml_model(source)
    chemistry_locks = protected_chemistry(model)
    energy_locks = protected_definitions(model)
    assert len(chemistry_locks) == 13 and len(energy_locks) == 8
    species = []
    adjacency = []
    reactions = set()
    for mid, role in TARGETS.items():
        met = model.metabolites.get_by_id(mid)
        species.append(dict(id=mid, role=role, name=met.name, compartment=met.compartment,
                            formula=met.formula, charge=met.charge,
                            annotation=met.annotation, notes=met.notes))
        for reaction in sorted(met.reactions, key=lambda r: r.id):
            adjacency.append(dict(metabolite_id=mid, reaction_id=reaction.id,
                                  coefficient=reaction.metabolites[met]))
            reactions.add(reaction.id)
    table(audit / 'target_species.tsv', species)
    write_json(audit / 'target_species_full.json', species)
    table(audit / 'target_adjacency.tsv', adjacency)
    rows = []
    for rid in sorted(reactions):
        r = model.reactions.get_by_id(rid)
        rows.append(dict(reaction_id=rid, name=r.name, equation_ids=r.reaction,
                         equation_names=equation(model, coefficients(r)),
                         lower_bound=r.lower_bound, upper_bound=r.upper_bound,
                         gpr=r.gene_reaction_rule, compartments=sorted(r.compartments),
                         stoichiometry=coefficients(r),
                         target_coefficients={m.id: c for m, c in r.metabolites.items() if m.id in TARGETS},
                         boundary=not r.reactants or not r.products,
                         allowed_directions='closed' if r.bounds == (0, 0) else
                             ('forward;reverse' if r.lower_bound < 0 < r.upper_bound else
                              'forward' if r.upper_bound > 0 else 'reverse'),
                         annotation=r.annotation, notes=r.notes, **balance(r)))
    table(out / 'glypro_route_audit.tsv', rows)
    table(audit / 'core_reactions.tsv', [r for r in rows if r['reaction_id'] in CORE])
    h = model.metabolites.get_by_id('m1007[C_va]')
    table(audit / 'all_vacuolar_proton_reactions.tsv', [dict(
        reaction_id=r.id, name=r.name, equation=r.reaction, lower_bound=r.lower_bound,
        upper_bound=r.upper_bound, hva_coefficient=r.metabolites[h], gpr=r.gene_reaction_rule)
        for r in sorted(h.reactions, key=lambda r: r.id)])

    # Exact columns plus proportional forward/reverse columns, without reliance on names.
    equivalents = []
    for r in model.reactions:
        sc = coefficients(r)
        if set(sc) != set(CY_HYD):
            continue
        ratios = {sc[k] / CY_HYD[k] for k in sc}
        if len(ratios) == 1:
            equivalents.append(dict(id=r.id, scale=ratios.pop(), bounds=list(r.bounds), gpr=r.gene_reaction_rule))
    identity_rows = []
    def stem(name):
        return re.sub(r'_[A-Z].*$', '', name).lower().replace('l-glycine', 'glycine')
    for mid in TARGETS:
        target = model.metabolites.get_by_id(mid)
        for other in model.metabolites:
            if other.id == mid:
                continue
            matches = [key for key in ('inchi', 'inchikey') if target.annotation.get(key)
                       and target.annotation.get(key) == other.annotation.get(key)]
            if stem(target.name) == stem(other.name):
                matches.append('normalized_name')
            if matches:
                identity_rows.append(dict(target_id=mid, other_id=other.id, other_name=other.name,
                    same_compartment=target.compartment == other.compartment,
                    formula_equal=target.formula == other.formula, charge_equal=target.charge == other.charge,
                    matching_fields=matches, target_formula=target.formula, other_formula=other.formula,
                    target_charge=target.charge, other_charge=other.charge,
                    caveat='Candidate match only; mixed legacy xrefs are not a chemical equivalence proof'))
    table(audit / 'identity_match_candidates.tsv', identity_rows)

    # Derive local route sums by summing stored columns. These are conditional algebra,
    # not optimized witnesses and not observations of a native protein.
    routes = {
        'O_CY__C_ONLY': {'HYP_GLYPRO_HYD_CY': 1},
        'O_CY__V_ONLY': {'R2038': 1, 'R2039': 1, 'R2030': -1, 'R2040': -1, 'R1363': 1, 'R795': .5},
        'O_VA__V_ONLY': {'R2039': 1, 'R2030': -1, 'R2040': -1, 'R1363': 1, 'R795': 1},
        'O_VA__C_ONLY': {'R2038': -1, 'HYP_GLYPRO_HYD_CY': 1, 'R795': .5},
    }
    sums = []
    for label, weights in routes.items():
        net = defaultdict(float)
        for rid, factor in weights.items():
            column = CY_HYD if rid == 'HYP_GLYPRO_HYD_CY' else coefficients(model.reactions.get_by_id(rid))
            for mid, value in column.items():
                net[mid] += factor * value
        net = {k: v for k, v in sorted(net.items()) if abs(v) > 1e-12}
        assert 'm1007[C_va]' not in net and 'm10[C_cy]' not in net
        assert 'm1863[C_va]' not in net and 'm1867[C_va]' not in net and 'm1384[C_va]' not in net
        sums.append(dict(route=label, claim_type='conditional_stoichiometric_derivation_not_optimization',
                         weights_per_q=weights, net_per_q=net, equation=equation(model, net),
                         pump_ATP_per_q=weights.get('R795', 0),
                         conditions='one specified source; selected hydrolysis=source=q; no export; no other H_va contribution; products recovered to cy'))
    table(audit / 'conditional_route_sums.tsv', sums)
    hyp = Reaction('HYP_GLYPRO_HYD_CY', lower_bound=0, upper_bound=.01)
    hyp.add_metabolites({model.metabolites.get_by_id(k): v for k, v in CY_HYD.items()})
    assert balance(hyp)['element_status'] == 'balanced' and balance(hyp)['charge_status'] == 'balanced_as_stored'

    prior = MODEL.reports / 'glypro_localization_predictors_20260924'
    sequence = json.loads((prior / 'sequence_identity.json').read_text())
    evidence = dict(gene_id='YALI1E16433g', established_native_name=None,
        protein_function='M24B/X-Pro peptidase, prolidase-like candidate',
        function_evidence='sequence and existing AlphaFold prediction; target free Gly-L-Pro activity unconfirmed',
        fixed_sequence=sequence,
        localization_evidence='historical completed predictions reused, not rerun or experimentally confirmed',
        predictions={
            'DeepLoc_2.1': {'cytoplasm': [.7816, .4761], 'nucleus': [.5657, .5014],
                            'lysosome_vacuole': [.0989, .5848], 'soluble': [.86, .5],
                            'numeric_pairs': 'score, threshold; not fractions or fluxes'},
            'SignalP_6.0': 'OTHER; no classical N-terminal secretion signal',
            'TargetP_2.0': 'OTHER; no recognized SP/mTP leader',
            'DeepTMHMM_1.0.57': 'GLOB, 0 transmembrane regions'},
        retained_localization_candidates=['cytoplasm', 'nucleus'],
        nucleus_not_modelled='Retained as evidence candidate; no nuclear substrate pool/transport/function evidence or authorization to model this round',
        biological_scope='No native activity or localization validation; no functional GPR; growth/cost do not rank localization evidence',
        next_discriminating_observation='Paired full-length native protein localization/fractionation (cytoplasm, nucleus, vacuole) and target-dependent free Gly-L-Pro disappearance plus Gly/L-Pro appearance; include purity and inactive/blank controls',
        old_online_blast='Historical BBH56DR2014 WAITING/unretrieved; not queried this round and not a prerequisite',
        source_files={str(p.relative_to(ROOT)): sha(p) for p in [prior/'REPORT.md', prior/'AUDIT.md',
            prior/'localization_predictions.tsv', prior/'sequence_identity.json', prior/'experimental_localization_evidence.tsv']})
    write_json(audit / 'biological_evidence_card.json', evidence)

    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(), scope='static read-only; zero optimize calls',
        source=str(source), source_sha256=sha(source), source_build_manifest=str(source.with_suffix('.build.json')),
        source_build_manifest_sha256=sha(source.with_suffix('.build.json')),
        script=str(Path(__file__).resolve()), script_sha256=sha(Path(__file__)),
        targets=len(TARGETS), adjacency_entries=len(adjacency), unique_adjacent_reactions=len(rows),
        vacuolar_proton_reactions=len(h.reactions), equivalent_cytoplasmic_hydrolysis=equivalents,
        energy_locks=len(energy_locks), chemistry_locks=len(chemistry_locks),
        existing_hypothesis_id='HYP_GLYPRO_HYD_CY' in model.reactions,
        proposed_cytoplasmic_balance=balance(hyp), same_compartment_identity_candidates=[r for r in identity_rows if r['same_compartment']],
        source_unchanged=sha(source)==EXPECTED_SHA)
    write_json(audit / 'manifest.json', manifest)
    report = '''# Gly-L-Pro 完整路线静态审计

本次仅读取指定化学候选、实际计量和既有证据，**没有优化调用**。输入完整身份、保护记录及脚本身份在 manifest.json；所有14个指定代谢物的全部直接邻接反应逐项导出，无按通量或名称筛除。

## 身份与原始表示

三处完整二肽均为候选 Gly-L-Pro，C7H12N2O3、charge=0，CHEBI:70744、PubChem 3013625，InChIKey KZNQNBZMBZJQJO-YFKPBYRVSA-N；Gly在N端、L-Pro在C端，为普通酰胺键、中性游离母体表示。L型/键型仍是先前候选整理的条件性映射，不能称原模型已记录该立体化学或已证明原生存在；不与Pro-Gly、D-Pro或其他微观态互换。

Gly两池C2H5NO2、charge=0，Gly本身无手性；胞质旧名含L-glycine但没有直接InChI，不能据此前缀创造手性。Pro两池C5H9NO2、charge=0，存储InChI均支持L构型。水为H2O/0；质子H/+1。ATP、ADP、Pi存储式/charge是中性母体，但保留的legacy InChI含/p-4、/p-3、/p-2，存在不同微观态注释不一致；本轮保留字段，不将外部结构微观态静默替换进S矩阵。完整字段、结构与notes见 target_species_full.json；名称/精确结构匹配候选另表记录，不凭名称合并池。

## 完整路径与可能去路

- R2036为胞外Gly-Pro向外界的exchange，[0,1000]，摄取关闭但外排保留。R2037负方向允许完整二肽由胞质向胞外，同时H+由胞质向胞外。
- R2037与R2038虽叫proton antiport，实际正方向分别是Gly-Pro和H+共同ex→cy、cy→va；两者可逆。按实际存储方向解释为共转运，不能按名称反转质子。
- R2039虽名为“cytosol nonspecific dipeptidase”，实际物种全在va；GPR为空。目标三个Gly-Pro池仅邻接R2036–R2039，没有其他目标水解、合成、去路；完整外排必须经R2037/R2036。
- Gly_va除了R2030，还有四条二肽水解生产；Pro_va只有R2039/R2040。R2030/R2040正方向是AA+H+从cy→va，产物回收须取负通量，同时H+从va→cy。无另一条目标产物跨液泡连接。
- R1363是cy→va供水；va水只供四水解。R795保持ATP+H2O+2H_cy→ADP+Pi+2H_va，2H/ATP只是当前计量。
- 液泡质子全部22个邻接列已输出，不能预先略去“其他项”。除核心外，R2020/R2028可逆，R871/R876开放，R882正向开放；其他列关闭。R882是Ile cy→va/H va→cy，Ile_va仅R882与已关闭R883相连，所以其稳态强制R882=0。
- 其他三二肽保留原网络。每一种的三个完整二肽池总平衡为−exchange−hydrolysis=0，两项不可逆且没有内部生产/本轮来源，因此它们都为0；再由各区室平衡得到各自完整二肽运输为0。Asp/Glu液泡池只有相应水解及输出（其余连接关闭），故R871/R876也不能产生独立的H+补偿。此推导不对模型额外添加约束；运行时仍须核对完整S*v。

## 条件性物料推导（不是优化结果）

一个指定来源=q、目标水解=q、另一水解=0且无其他来源时，三处Gly-Pro总平衡给出R2036=0，不能再把未水解的外排当作处理量。两处VA产物回收分别为R2030=R2040=−q。上述其他来源关闭且对应局部稳态成立时，其他H_va项为0；完整式仍为所有22列求和。

在这些明确条件下，对存储S列按权重相加可得每q的泵ATP系数：CY→CY为0、CY→VA为0.5、VA→VA为1、VA→CY为0.5。完整加权列和净反应在 conditional_route_sums.tsv。CY→VA的完整二肽输入已带入1H/q，抵消一半产物回收的质子输出；VA→CY的完整二肽输出本身带走1H/q，因此即使水解发生在胞质，也有补偿需求。各式水解本身都只耗1水/q，额外水由R795的ATP水解计量产生。

这些是局部计量与明确边界的条件性结论，**不证明全细胞可行、达到指定生长或达到容量，也不是已观测的ATP成本**。实际优化、独立泵关闭、完整全细胞账本需另行检验。预测分数、较高生长或较低泵耗均不构成原生定位证据。

## 目标蛋白证据

YALI1E16433g（原生正式名未核实）为M24B/X-Pro肽酶、prolidase-like功能候选，来自序列及既有AlphaFold预测；固定W29/PO1f 454 aa，目标游离Gly-L-Pro水解尚无实验确认。历史正式DeepLoc支持胞质及细胞核、可溶标签；液泡类别未过阈值，SignalP/TargetP无经典分泌/所识别前导肽，DeepTMHMM GLOB/0跨膜。保留核定位为证据候选，但本轮没有核内底物/运输/功能建模范围与实证，故不建核路线。未重新运行预测、旧在线BLAST或实验；不赋正式GPR。

最能区分剩余假设的观察是将原生完整蛋白在胞质、细胞核和液泡的定位/组分分离，与同组分中目标依赖的游离Gly-L-Pro减少及Gly/L-Pro生成配对；需纯度、失活/空白对照，液泡检出要区分完整活性蛋白与降解片段。这是实验建议，不是本次执行结果。
'''
    (audit / 'REPORT.md').write_text(report)
    assert sha(source) == EXPECTED_SHA
    print(json.dumps({k: manifest[k] for k in ['targets', 'adjacency_entries', 'unique_adjacent_reactions',
        'vacuolar_proton_reactions', 'equivalent_cytoplasmic_hydrolysis', 'same_compartment_identity_candidates', 'source_unchanged']}, ensure_ascii=False))


if __name__ == '__main__':
    main()
