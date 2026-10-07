# 二肽结构与质子化表示：有条件的化学候选

核验日期：2026-09-24。范围为开放候选中的四种二肽 × ex/cy/va 三池，以及液泡 glycine。原 E5、原开放候选和历史记录均未由本子任务改写；本子任务没有优化调用。

本次已确认五个**数据库结构实体**，并补强了 glycine 的无手性身份；这与“原模型已经明确指定这十二个二肽池的精确立体结构”是两种结论。`data/iyli21.xml` 与 `data/iyali26.xml` 的相关条目只有名称、区室和通常属性，没有可以追溯的底物立体化学/精确结构 accession。因此十二个二肽池均保留**明确 L 型普通肽键、游离二肽命名假设下的结构候选**。对应水解产物的 L-Asp、L-Glu、L-Ala、L-Pro 名称支持这一解释，但不能补造原作者的结构记录。

这不是拒绝补式的理由：在独立候选中可以明确选择这些已核查的结构，并使 formula、charge、SMILES 和结构注释相互一致。新的确定性只针对所选化学表示，不能外推到原生存在、内源供给、催化酶或液泡定位。

| 选择的自由母体 | N端 → C端 | 中性母体式；电荷 | 官方结构实体 | 本模型映射的证据等级 |
|---|---|---|---|---|
| glycyl-L-aspartic acid | Gly → L-Asp | C6H10N2O5；0 | [ChEBI:73804](https://www.ebi.ac.uk/chebi/CHEBI:73804) / [PubChem 97363](https://pubchem.ncbi.nlm.nih.gov/compound/97363) | 三池均为明确命名假设下的候选 |
| glycyl-L-glutamic acid | Gly → L-Glu | C7H12N2O5；0 | [ChEBI:73801](https://www.ebi.ac.uk/chebi/CHEBI:73801) / [PubChem 99278](https://pubchem.ncbi.nlm.nih.gov/compound/99278) | 三池均为明确命名假设下的候选 |
| L-alanylglycine | L-Ala → Gly | C5H10N2O3；0 | [ChEBI:73757](https://www.ebi.ac.uk/chebi/CHEBI:73757) / [PubChem 6998029](https://pubchem.ncbi.nlm.nih.gov/compound/6998029) | 三池均为明确命名假设下的候选 |
| glycyl-L-proline | Gly → L-Pro | C7H12N2O3；0 | [ChEBI:70744](https://www.ebi.ac.uk/chebi/CHEBI:70744) / [PubChem 3013625](https://pubchem.ncbi.nlm.nih.gov/compound/3013625) | 三池均为明确命名假设下的候选 |
| glycine | 不适用 | C2H5NO2；0 | [ChEBI:15428](https://www.ebi.ac.uk/chebi/CHEBI:15428) / [PubChem 750](https://pubchem.ncbi.nlm.nih.gov/compound/750) | 命名的无手性骨架可确认；中性表示是本轮明确选择 |

## 结构与组成核对

原始记录保存于 `chemistry/pubchem_five_compounds.json` 与各 `chebi_*.html`。`exact_structures.json`、`dipeptide_identity.tsv` 逐项保留完整 SMILES、InChI、InChIKey 和结构来源。检索时间、响应状态、URL 与 SHA 见 `initial_retrieval.json`、`chebi_retrieval.json`、`additional_retrieval.json`。

五个 PubChem 原始原子图都含显式 H 原子。本次用标准库按原子序数独立计数元素，按每原子形式电荷求和，再核查所有原子属于一个共价连通分量，以及中性 C/H/N/O 的键价。五项均与 PubChem 报告式、ChEBI 报告式和电荷一致；两数据库的标准 InChI 也一致。这是**对下载结构的计算交叉验证**，不是新增实验来源。当前项目 Python 未安装 RDKit；没有声称进行 RDKit 立体化学重建。可复跑已实际执行的命令：

```bash
.venv/bin/python artifacts/dipeptide_chemistry_routes_20260924/chemistry/analyze_identity.py
```

每个选择的结构只有一个非 Gly α-碳立体中心，数据库指定为 L / S；glycine 无手性。选择普通 α-羧基与下一个残基 α-氨基之间的酰胺键；Gly-Pro 中酰胺氮为 Pro 环内氮。所下载结构为单个未保护、未成盐、未水合的游离分子，没有反离子或保护基。没有把含 Gly-Pro 的长肽、Gly-Pro-pNA 或环二肽代入。

## 同分子式不等于同一底物

已另外取得直接对照结构：

- Gly-L-Asp 的反向序列 [L-Asp-Gly / ChEBI:73450](https://www.ebi.ac.uk/chebi/CHEBI:73450) 同为 C6H10N2O5，但连接 InChIKey 为 `JHFNSBBHKSZXKB-VKHMYHEASA-N`，与目标 `SCCPDJAQCXWPTF-VKHMYHEASA-N` 不同。
- Gly-L-Glu 的反向序列 [L-Glu-Gly / ChEBI:73505](https://www.ebi.ac.uk/chebi/CHEBI:73505) 同为 C7H12N2O5，但键连接不同；[γ-Glu-Gly / PubChem 165527](https://pubchem.ncbi.nlm.nih.gov/compound/165527) 又通过 Glu 侧链羧基连接，不能充当任一个普通 α-肽键实体。
- L-Ala-Gly 的反向序列 [Gly-L-Ala / ChEBI:73855](https://www.ebi.ac.uk/chebi/CHEBI:73855) 同式，但 InChIKey 为 `VPZXBVLAVMBEQI-VKHMYHEASA-N`，与目标 `CXISPYVYMQWFLE-VKHMYHEASA-N` 不同。
- Gly-L-Pro 的反向序列 [L-Pro-Gly / PubChem 6426709](https://pubchem.ncbi.nlm.nih.gov/compound/6426709) 同为 C7H12N2O3，但 InChIKey 为 `RNKSNIBMTUYWSH-YFKPBYRVSA-N`，与目标 `KZNQNBZMBZJQJO-YFKPBYRVSA-N` 不同。
- [PubChem 273261](https://pubchem.ncbi.nlm.nih.gov/compound/273261) 的 Gly-Asp 条目不编码目标 L 异构体的确定立体结构；不能仅因同式就替代 CID 97363。名称中 DL 与一个单一指定 L 结构不能合并。

## 为什么选择中性自由母体

本候选采用数据库给出的**不带原子形式电荷的自由母体结构，总电荷为 0**。这与当前模型中所加载的 Asp C4H7NO4/0、Glu C5H9NO4/0、Ala C3H7NO2/0、Pro C5H9NO2/0 及胞质 Gly C2H5NO2/0 的元素/总电荷约定一致。该选择依据数据库结构和现有数值表示；不是把产物相加减水反推一个能配平的式，也不是为保留历史泵通量而挑式。

这个表示不宣称是液泡、胞质或培养基在某个 pH 下的优势微观态。三池选择相同结构和总电荷，以保持已有转运所表达的同骨架对应。不按区室酸性自动加 H⁺，不修改共享 ATP、Pi、H2O、H⁺，不改变任何反应计量。

中性母体、两性离子和净负离子必须区分：

- [Ala-Gly 两性离子 / ChEBI:73786](https://www.ebi.ac.uk/chebi/CHEBI:73786) 具有 `[NH3+]` 与 `[O-]`，式和总电荷仍为 C5H10N2O3/0，但精确的原子电荷布局不同。
- [Gly-Pro 两性离子 / ChEBI:73779](https://www.ebi.ac.uk/chebi/CHEBI:73779) 同理。上述母体/两性离子标准 InChIKey 可相同，**因此 InChIKey 相同本身不足以确认同一微观质子化态**；本轮同时保留对应的完整带/不带原子电荷 SMILES。
- [Gly-Glu(1−) / ChEBI:73784](https://www.ebi.ac.uk/chebi/CHEBI:73784) 是 C7H11N2O5 / −1；不能将它的负离子结构与中性母体 C7H12N2O5/0 混写为同一精确表示。

当前 Asp/Glu 产品的已有 InChI 含 `/p-1`，但加载字段使用上述中性分子式和总电荷 0。这是**原有注释表示不一致**，不能将“与当前加载式/电荷一致”写成“现有所有结构注释均一致”。本轮目标只包含十三个指定池，故仅记录该矛盾，未自动修改这些共享产物的式、总电荷、InChI 或相关反应。

## Gly 的旧交叉引用

`m1863[C_va]` 旧名为 `L-glycine`，并带 `MNXM12053` / `cpd16187`。2026-09-24 打开官方 [MetaNetX MNXM12053](https://www.metanetx.org/chem_info/MNXM12053)，该条目也是 L-glycine，只有 ModelSEED 来源，式、SMILES、InChI 和 InChIKey 均空。官方 [ModelSEED compounds.json](https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.json) 中的 cpd16187 为来源 Published Model 的 iAG612:cbs_407 名称条目，同样缺式和结构。

因此没有发现这些旧引用指向另一个确定分子；它们只是不能提供精确结构验证。候选建议保留为历史引用，并在 notes 中注明上述限制；新增 ChEBI:15428/PubChem 750 的准确 glycine 结构，标准名改为 glycine。胞质、线粒体和胞外现有 glycine 池使用 C2H5NO2/0，液泡通过既有 R2030 直接对应胞质池；这支持本次填补一致的化学表示，而不是证明原生液泡浓度或运输功能。

## 待完成的审查边界

本子任务只交付字段建议和结构证据。实施器还须完成全量前态检查、幂等与冲突拒绝、导出重载校验和数学定义不变检查。邻接残差应按选定字段直接计算，缺式不能被当成 0 残差而宣告配平。没有任何原生酶/GPR/定位结论由本化学表自动成立。源证据的独立审核由另一审查代理完成，其最终覆盖和未决项应随主报告报告。
