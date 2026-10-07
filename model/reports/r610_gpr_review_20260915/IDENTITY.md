# R610：目标蛋白身份、注释与区室

核验时间：2026-09-15（America/Los_Angeles；官方获取 UTC 2026-09-16 00:19–00:20）。仅审查证据与复用已有结构；未修改模型/GPR，未做 FBA，未提交结构预测。

## 结论

**YALI1C05908g — 原生正式名称未核实 — 推定鸟氨酸 δ-氨基转移酶、PLP 依赖（模型/GPR赋值，加自动同源与结构域注释；未取得该 W29 蛋白的直接实验酶学证据）。**

官方版本化序列明确来自 W29/CLIB89，与项目的历史 UniProt/FASTA 和现成 AlphaFold 模型均 432 aa 完全一致。因此 R610 所用蛋白的身份可追溯；酶功能、底物范围和区室仍须与序列身份分开。现有证据支持继续把它作为鸟氨酸转氨酶候选，不能从家族名字自动推出其他底物、可逆通量范围、AND/OR 关系或胞内定位。

## 版本化身份

| 记录 | 本次实际核实 |
|---|---|
| [AOW02333.1](https://www.ncbi.nlm.nih.gov/protein/AOW02333.1) | W29/CLIB89；YALI1_C05908g；432 aa；CP017555.1:590828..592126；记录日期 2016-10-25；定义为 hypothetical protein |
| [XP_501441.1](https://www.ncbi.nlm.nih.gov/protein/XP_501441.1) | W29/CLIB89；YALI1_C05908g；old_locus_tag YALI2_C00075g；GeneID 2909646；432 aa；记录日期 2024-09-10；PROVISIONAL RefSeq、定义为 uncharacterized protein；序列与 AOW02333.1 精确相同 |
| 当前核酸版本 | XP_501441.1 对应 XM_501441.3；历史 UniProt 交叉链接仍为 XM_501441.1。蛋白版本未变化不意味着核酸条目版本或菌株元数据未变化 |
| 历史 YALI0_C04433g | NCBI 记录明确列作比较对象；模型映射与此一致。本次没有把旧标识直接当作当前 W29 名，也没有新增跨版本猜测映射 |

NCBI CDS 注释指出与 Emericella nidulans 的 Q92413 鸟氨酸转氨酶高度相似。该跨物种说明是功能推断，不是本原生蛋白被实验验证的证据。未把酿酒酵母 CAR2 或其他物种基因名称当作目标已核实名称。

## 历史注释与当前状态

历史 `A0A1H6PJM3_cached.json` 是条目 v36、序列 v1（序列更新 2017-01-18、注释更新 2026-01-28），unreviewed TrEMBL；蛋白存在依据是同源推断。名称 Ornithine aminotransferase、EC 2.6.1.13、PLP 辅因子和鸟氨酸生成 GSA 的步骤均为 ECO:0000256 自动规则注释。NCBI CDD 显示 Orn_aminotrans/TIGR01885 域（2–422）及推定 PLP 位点/催化 K272；这些也是计算域注释。

本次 [UniProt 官方 API](https://rest.uniprot.org/uniprotkb/A0A1H6PJM3.json) 返回 Inactive/DELETED，原因为不属于参考蛋白组。此状态不是酶功能被否定。报告中的功能条目须明确来自历史缓存，不能表述成当前在线 reviewed 注释。

同一历史 UniProt 条目还包括 YB392 菌株的 B0I71DRAFT_133637/RDW24840.1。AlphaFold API 采用这个 gene 字段；本次以全序列严格相等来确认可复用性，而不是把 API 的菌株相关字串当 W29 身份证据。

## 区室和证据等级

历史条目 GO:0005829（胞质）采用 IEA:EnsemblFungi；未见实验 subcellular-location 评论。本次可将胞质作为有自动注释支持的候选定位，不能标成实验确认。PANTHER 亚家族名称含 mitochondrial，是家族条目名称，**不构成本蛋白在线粒体定位的证据，也不能用于推翻或替代胞质注释**。序列长度或 N 端结构置信度也不能确证定位。

## 复用已有 AlphaFold 预测

[AF-A0A1H6PJM3-F1](https://alphafold.ebi.ac.uk/entry/A0A1H6PJM3)，下载模型文件 v6；API 标示 AlphaFold Monomer v2.0 pipeline，创建 2022-06-01，序列版本日期 2017-01-18。本次获取的是已有 **AlphaFold 预测**，没有运行新预测。API 全序列及 PDB 全部 432 个 CA 残基与官方 W29 序列完全一致。

API 平均 pLDDT 97.12；下载 PDB 的 CA 实算均值 97.1410、范围 60.78–98.94（保留两种原始表示的微小差异）。已取得逐残基置信度和完整 432×432 PAE，文件标示最大 PAE 31.75 Å。高结构置信度不证明原生催化、底物选择、反应方向、复合体或区室；本身份子任务未做结构比对/对接，不新增结构导出的功能判断。

## 保存的身份与可复查资料

`target.fasta` 起初从历史 UniProt 缓存抽取，已与两份官方版本化记录逐残基断言一致；保留原始标题。原始大写氨基酸序列（不含换行）SHA-256：`85fc7d9872b5eff7eed949e9ee55bbe3d5128b1d6b34f781ae1e22d07bb2786d`。

FASTA 文件 SHA-256：`c6ef322598821deab85efb9a14b0c50f77bdd8465e7bd4397db6d7aa4c07f1a0`。历史全缓存路径 `artifacts/reference_pipeline_restore_20260909/research/cache/data/uniprot_UP000182444.json`，SHA-256 `e5b0a04874079b4057ffe25dadcb6b812cba8c96b227f41187187b27715753a3`。`identity_checks.json` 保存官方记录/缓存 FASTA/AF 同一性检查及置信度；`identity_sources/retrieval_network.json` 和 `retrieval_af_files.json` 保存每个官方 URL、获取时间、字节数和 SHA，原文同目录保留。所有核心身份结论基于实际取回记录；搜索中第三方陈旧 CLIB122 标题没有用于当前身份定论。
