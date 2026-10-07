# R293 / R2080：独立来源审计

核查日期：2026-09-11。审计问题是两条目标蛋白能否各自完成胞质 ATP + ethanolamine → ADP + phosphoethanolamine，及其是否具有相互替代或共同必需的直接证据。审计者独立检索并打开原始数据库/论文，先形成来源判断，再抽核主审计者保存的结构与静态证据。未修改模型、GPR、培养或实验标签，未执行优化或提交集群作业。

| 系统 ID | 已核实名称/符号与蛋白功能 | 证据状态 | 当前模型角色 |
|---|---|---|---|
| YALI1B12742g；旧版映射 YALI0B09515g | 2024 论文表 2 使用 **EKI1**；乙醇胺激酶候选 | 文献用名已核实；功能为自动/同源注释，非目标蛋白直接生化验证 | R2080 的单基因 GPR |
| YALI1E20159g；旧版映射 YALI0E16907g | 未核实独立确定的原生通用名；胆碱/乙醇胺激酶家族候选 | 历史 UniProt 为未审阅、PE=3；两种激酶活性均为电子推断 | R293 的单基因 GPR，另用于胆碱步骤 R212 |

审计判断：数据库和原始文献足以支持保留两者为待判定的激酶候选；尚不足以把“各自独立催化同一乙醇胺步骤”升级为已实验确认。没有找到目标两蛋白必须共同组成催化复合体的直接证据，因此不能据此改为 AND。

## 原子声明

`partially_supported`、`unsupported` 和 `unverified` 均计入 unresolved；unsupported 表示所核来源不支持该精确声明，不表示已证明该声明为假。

| ID | 被审声明 | 判定 | 来源与定位；边界 |
|---|---|---|---|
| A01 | 两个模型 UniProt xref 对应本次 W29 目标蛋白序列 | supported | NCBI Gene eSummary、UniProt 最后活跃版本 GN/RC/DR 及独立获取的 AOW01464.1/AOW05525.1 FASTA 一致。历史 UniProt 同时含 W29 与 YB392 来源；不能将其说成仅来自 YB392。 |
| A02 | 两个 UniProt accession 当前已 inactive | supported | 两条当前 JSON 的 `inactiveReason` 均为 DELETED / Not part of a reference proteome；不构成蛋白失去功能或物种不存在的证据。 |
| A03 | YALI0B09515g 在原始论文中被称为 EKI1 | supported | Carmon 等 2024，表 2。该用名来源是 RNA-seq 结果表的功能注释，不能等同于独立命名/纯化酶鉴定。 |
| A04 | YALI0E16907g 可确定命名为 CKI1 | unverified | Matsuse 等 2024，Discussion 仅将两旧版 ID 并列为 EKI1 or CKI1 同源候选；没有单独锁定该名称。 |
| A05 | YALI1B12742g 已被直接实验证明能独立催化目标步骤 | partially_supported | UniProt v33 的 EC 2.7.1.82 来自 ARBA，PE=3；论文表 2 是转录变化，无该目标底物谱或单基因互补实验。支持候选功能，未达到“直接实验”标准。 |
| A06 | YALI1E20159g 已被直接实验证明能独立催化目标步骤 | partially_supported | UniProt v34 GO:0004305 和 GO:0004103 均 IEA:TreeGrafter；DE 仅 kinase-like domain-containing protein。论文同源讨论不能单独证明底物谱。 |
| A07 | 两目标均已有胞质定位实验证据 | partially_supported | 两条 UniProt 的 GO:0005737 均 IEA:TreeGrafter。2024 论文的定位实验针对其他蛋白，不是这两个激酶。 |
| A08 | 两目标在 Yarrowia 中已实验证实为可互相替代的乙醇胺激酶 | unsupported | 本次所核原始论文均未测这两个目标的成对单敲/双敲、分别互补或分别催化活性。数据库家族与GO不能直接确认该OR关系。 |
| A09 | 两目标在 Yarrowia 中已实验证实为共同必需复合体 | unsupported | 未找到目标蛋白复合体依赖或共同必需的直接实验。单体 AlphaFold 预测不提供这一证据。 |
| A10 | 酿酒酵母的激酶研究可直接确定两个 Yarrowia 目标的 GPR | unsupported | Kim 等 1999 的生化/遗传实验属于 S. cerevisiae；功能重叠的物种间转移还需目标序列与独立功能证据。 |
| A11 | 已有 AlphaFold 模型与已核 W29 目标序列完全一致 | supported | 独立读取两 PDB 的 CA 残基并与 NCBI FASTA/历史 UniProt 逐残基比较，均完全一致。仅说明输入身份，不能证明功能。 |
| A12 | 当前已核模型/培养条件下 R293 与 R2080 在所有稳态可行解中都为零 | supported | 独立遍历实际模型 XML 全部反应，核得完整胞内/胞外乙醇胺守恒行和三条非负反应边界，见下文。此结论限于乙醇胺摄入关闭且当前供给结构保持不变。 |

覆盖：**12 total | 12 audited | 5 supported | 7 unresolved | 0 contradicted | 0 unchecked**。此计数覆盖表内限定声明，不代表穷尽全部文献或完成生物学验证。

## 原始论文核查

[Matsuse 等，2024，FEMS Yeast Research，doi:10.1093/femsyr/foae030](https://academic.oup.com/femsyr/article/doi/10.1093/femsyr/foae030/7760385)：本次打开出版社全文及 XML。Results 首段与 Fig. 1 报告 YALI0D08514g—PSS1—磷脂酰丝氨酸合酶（该文实验支持）缺失株在 1 mM 乙醇胺/胆碱补充下恢复生长；菌株为 CXAU 系列，30°C、葡萄糖最小培养条件。它支持 Yarrowia 可利用这些前体，但未分别干预目标两激酶。Discussion 中两目标 ID 只出现为同源候选；Fig. 4 的定位与功能互补针对 Pss1，不能转用于目标激酶。本文未逐图复核培养皿像素，故这里只报告正文/图注记载，不声称本次复现生长。

[Carmon 等，2024，doi:10.1016/j.bbalip.2024.159544](https://pmc.ncbi.nlm.nih.gov/articles/PMC11380575/)：本次打开 PMC 全文 XML。表 2 将 YALI0B09515g 列为 EKI1/乙醇胺激酶，检测内容是磷脂酸磷酸酶缺失背景下的转录变化。Methods 的 RNA-seq 使用每菌株三个 isolates，低/高糖条件；表 1 为 Po1d 衍生菌株。本研究没有对目标激酶做底物、定位或相互依赖实验；该名称可用于“文献称 EKI1”，不能写成原生酶功能已经验证。

[Kim 等，1999，doi:10.1074/jbc.274.21.14857](https://pubmed.ncbi.nlm.nih.gov/10329685/)：本次核对 PubMed 原始摘要 XML，未核全文实验表。S. cerevisiae YDR147W—EKI1—乙醇胺激酶（原始实验支持）在昆虫细胞异源表达，并在相应酵母缺失背景中增加激酶活性；与 YLR133W—CKI1—胆碱激酶兼乙醇胺激酶（该研究生化/遗传证据支持）存在功能重叠。这支持“酵母中可能存在部分底物重叠”的机制参照，不支持将两个酿酒酵母名称逐一自动赋给本次 Yarrowia 目标，更不支持 AND。

## 序列/结构与静态证据独立抽核

本次从 NCBI 独立取得 AOW01464.1 和 AOW05525.1。B 蛋白 463 aa，序列 SHA-256 `4e4a1ae5b8b458768b4dc62e99501327f9a91b2f5376d856f71c49ae3f6d4f5f`；E 蛋白 566 aa，`d8168c9096f595fda852cf3fb9f2ba6556525ebed4dd3ef0810afa04b33f1244`。两者分别与 UniProt v33/v34 及 AF PDB 完全一致。

独立重算 **AlphaFold 预测** 的 CA pLDDT 与全矩阵 PAE：E 平均 pLDDT 75.606855、169/566 个残基低于 70，PAE 均值 14.485741 Å；B 为 78.198790、114/463、12.133951 Å。API 标注模型日期 2022-06-01，文件为 AFDB v6，PDB TITLE 为 AlphaFold Monomer v2.0。日期/文件版本/预测器版本应分开记录。以上是既有预测的复用与输入核验；没有做实验参考酶的结构对齐、活性位点验证或多聚体依赖推断，不能因此升级底物/区室结论。

独立以 XML 标准库直接读取实际模型 `model_metadata_trna.xml`（SHA-256 `d274bad3050e3c9220a8b6287eae847f3bf1334892284d565a6c4d96b38135a0`），遍历全部反应，核得胞外乙醇胺行 `−v846−v1115=0`，胞质行 `v846−v293−v2080=0`，两物种 `boundaryCondition=false`。相加得到 `v1115+v293+v2080=0`；当前三者下界均为零，故三者必须为零。该证明与主审计静态证据一致，说明在当前条件中即使删除重复步骤也不能由该步骤制造生长差异；若开放乙醇胺摄入或增加内部来源则需重新判断。本次未独立重复全 WT 通量向量检查或最优性求解；那部分复核等级以主审计日志为准。

## 来源记录与未决项

独立取得的原始全文、摘要和 NCBI FASTA 保存在 `sources/source_audit/`；同目录 `manifest.json` 记录 URL、获取时间、字节数及完整 SHA。UniProt/AF 原始文件由主审计取得，本审计打开其原始记录并独立计算；对应来源见 `sources/fetch_manifest.json`。

有限检索覆盖两个 YALI1 ID、两个旧 YALI0 ID、两个 YALI2 别名、NCBI/UniProt accessions、以及 Yarrowia ethanolamine kinase / EKI1 / CKI1；“未找到直接证据”只指本次已打开来源和该检索范围。下一步真正限制 GPR 定论的是目标蛋白各自对乙醇胺/胆碱的底物活性、同条件功能互补与胞质定位，及必要时复合体依赖证据。已有单体结构与电子注释尚不能解决这些未决项。
