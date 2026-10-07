# R2200 / YALI1F23706g 身份与底物审查

核验：2026-09-24。只读当前工作区和公开官方来源；新文件仅写本目录。未运行 LP、未修改模型/GPR/边界，未提交新 AlphaFold 或 HPCC 作业。

## 判定

**YALI1F23706g（NCBI 写法 YALI1_F23706g；正式原生基因名未核实）可保留为 DUG1-like / M20A Cys-Gly 金属二肽酶候选。现有证据不足以将它赋给 Gly-Asp、Gly-Glu、Ala-Gly、Gly-Pro 四条液泡水解。** 支持主要来自同源注释、精确序列身份和既有 AlphaFold 预测；没有本次找到的 W29 纯化酶四底物实测或液泡定位实测。

这不是把四种底物判为阴性。酿酒酵母对照研究没有测这四种底物，不能从未测推出不水解。

## 身份已闭合到当前 W29 序列，数据库版本冲突已分开

| 对象 | 核实结果 | 证据性质 |
|---|---|---|
| NCBI Gene 2908360 | 当前 YALI1_F23706g；别名 YALI2_F00620g；正式 symbol 为空 | 官方基因记录 |
| XP_505554.3 / AOW07335.1 | 都明确 CLIB89(W29)，478 aa，逐字相同；前者为 provisional RefSeq、conceptual translation，2024-09-10 替代 .2 | 版本化序列身份，不是酶活证据 |
| A0A1D8NNX0 | 本地 entry v37 / sequence v1 为478 aa，正式名为 M20 dimerisation domain-containing protein，功能/金属位点/胞质GO皆自动或同源推定 | 未审阅候选 |
| 实时 A0A1D8NNX0 | 当前 inactive，原因是“不属于 reference proteome” | 数据库整理状态，不表示 W29 无此基因 |
| CLIB122 YALI0_F17842g / Q6C1A8 | 实时 entry v116 / sequence v1，478 aa，与当前 W29 序列完全相同 | 跨菌株相同蛋白序列；不自动转移表达、定位或培养条件证据 |
| KEGG yli:2908360 | 仍列 YALI2_F00620g、477 aa、旧位置与未版本化 XP_505554；K15428 / Cys-Gly DUG1、EC3.4.13.- | 同源功能注释，并有旧版本滞后 |
| 历史 XP_505554.2 | source strain为DSM3286，locus YALI2_F00620g；477 aa，与上述KEGG完全相同；相对当前W29序列末端重复KK少一个Lys | 涉及旧菌株参考和注释版本，不能说成同一菌株简单修订；不与当前478 aa静默混用 |

官方定位：[NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/2908360)、[当前蛋白](https://www.ncbi.nlm.nih.gov/protein/XP_505554.3)、[W29原始翻译](https://www.ncbi.nlm.nih.gov/protein/AOW07335.1)、[UniProt Q6C1A8](https://www.uniprot.org/uniprotkb/Q6C1A8/entry)、[KEGG](https://www.kegg.jp/entry/yli:2908360)。原始响应与下载SHA见本目录；精确序列比较可重跑 `python3 artifacts/dipeptide_energy_audit_20260924/identity_r2200/verify_identity.py`。

当前478 aa的序列 SHA256（大写氨基酸字母、无标题/换行）：`c6d94c48bc2debdb972c57bf68f7e7d4b9a9c06a61deef9ff6924b6c40f16cfe`。

## 四底物证据与 Cys-Gly 证据不能合并

| 精确底物 | 对 W29 YALI1F23706g 的直接实测 | 已核实同源/类别证据 | 本轮裁决 |
|---|---|---|---|
| Cys-Gly | 未取得 | 酿酒酵母 YFR044C / DUG1（已确证Cys-Gly金属二肽酶）有纯化酶及细胞遗传证据；W29有相符同源注释 | 支持原 R2200 的功能候选，尚非原生功能验证 |
| Gly-Asp | 未取得 | 2009 DUG1 Table 2 未测 | 不补 GPR |
| Gly-Glu | 未取得 | 2009 DUG1 Table 2 未测 | 不补 GPR |
| Ala-Gly | 未取得 | 2009 DUG1 Table 2 未测；2011 曲霉同源酶研究报告未检出，非W29反应阴性证明 | 不补 GPR |
| Gly-Pro | 未取得 | 2009 DUG1 Table 2 未测；EC3.4.13.18描述中的 prolyl substrates / Pro-X 不能等同 X-Pro | 不补 GPR；须独立核实X-Pro水解酶 |

酿酒酵母 YFR044C / DUG1（正式名已确证；Cys-Gly 金属二肽酶）的2009原始研究使用纯化重组蛋白，2 mM底物、200 μM Zn²⁺或20 μM Mn²⁺、30°C；Table 2 的底物均含 Cys，检测依赖释放 Cys。Cys-Gly 是该面板的优选底物，Ala-Cys 也有约74–75%相对活性；因此“严格只作用Cys-Gly”也过强。原文不包含四目标无Cys二肽，不能从该表制造四个阴性结论。[原始研究，Table 2 / assay methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC2682898/)

米曲霉 AO09002000015 / cdpA（按2011论文所写ID；Cys含量相关二肽酶，纯化重组酶实验）报告 Ala-Gly 和 Gly-Ala 未检出活性；该结果只提示同家族不能默认普遍水解。论文完整PDF本轮下载403，所用信息来自出版社可检索正文和PubMed摘要，未声称本地完整PDF审核。[原始研究 DOI:10.1271/bbb.100604](https://doi.org/10.1271/bbb.100604)

[IUBMB EC3.4.13.18](https://iubmb.qmul.ac.uk/enzyme/EC3/4/13/18.html) 说明的是物种相关的广谱二肽酶类别，不能为某一 Yarrowia 蛋白补齐所有具体底物。特别要保留 Gly-Pro 与 Pro-Gly 的序列顺序。

## 定位与结构

W29缓存的 cytoplasm GO 为 **IEA:EnsemblFungi**，不是原生定位实验。酿酒酵母 YFR044C / DUG1 的 C端GFP原始研究观察到胞质荧光，条件包括甲硫氨酸培养后GSH培养；这仅支持同源定位候选，不能据此给 W29 液泡腔赋值。[2007原始研究，Dug proteins localization段](https://pmc.ncbi.nlm.nih.gov/articles/PMC1840075/)

本轮复用 **AlphaFold 预测** `AF-A0A1D8NNX0-F1-model_v6`。AFDB 标注工具为 AlphaFold Monomer v2.0 pipeline，原模型创建2022-06-01，本轮2026-09-24获取；478个PDB Cα残基逐字等于上述当前W29序列。AFDB提供总体pLDDT96.56；本轮按478个Cα B-factor计算均值96.593、最低64.81、>90比例93.10%。478×478 PAE已取回，全部残基对均值3.805 Å、最大27 Å。两种pLDDT统计口径分列，未静默替换。

**这是高置信的既有单体预测，不能验证四底物、金属依赖、原生二聚体、液泡定位或方向性。** 本轮没有对接、催化口袋实验、参考结构叠合或新预测，不将预测结构升级为底物验证。AFDB“reference proteome”标签与实时UniProt停用状态不一致，身份以精确序列为准。来源：[AFDB元数据](https://alphafold.ebi.ac.uk/api/prediction/A0A1D8NNX0)。

## 与当前模型的关系及接受门槛

实际 R2200 是胞质 `Cys-Gly + H2O ⇌ Cys + Gly`，边界[-1000,1000]、单基因 YALI1F23706g；与四目标液泡二肽水解不同。当前 Cys-Gly 物种缺分子式；多个混并EC（3.4.11.*、3.4.13.18、3.4.13.-）是反应注释，不能当作此蛋白全部功能实证。精确模型SHA/计量已写 `sequence_identity.json`；本轮不改变这些字段。

接受四目标GPR至少需要：版本明确的原生候选、各精确底物的酶学或可靠物种内功能支持、相容区室/运输证据。若拟从液泡改到胞质，还须单独审查物料供给与产物去路；不能仅因候选的自动胞质注释直接搬反应。缺口限制正式接纳，但不阻塞已完成的身份和模型诊断。
