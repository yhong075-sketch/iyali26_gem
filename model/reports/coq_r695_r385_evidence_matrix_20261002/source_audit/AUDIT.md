# R695／R385 证据矩阵：独立来源审计

核验日期：2026-10-02（UTC）；审计者：独立子代理 `/root/coq_matrix_source_audit`。本轮是原始文献与现行分类的重新阅读和主张审核，不是实验复现、W29 酶活验证或模型接受。

范围与执行：按指定 C1–C7 审核；酿酒酵母为酵母实验主锚点，其他物种直接酶学为补充。读取已有全文缓存并重新打开 PMID 863914 和三项 IUBMB 现行条目；未新增搜索批次、论文范围、模型／序列／结构计算。Houser 1977 仅可读摘要，沿用已记录的全文获取失败范围，不把旧审计意见当作原文。完整来源路径、SHA、读取范围见 [provenance.json](provenance.json)。

基因名称对应：酵母 **YOR125C — CAT5（别名 COQ7）— CoQ 晚期羟化相关蛋白**，本轮采用实验验证的遗传功能及整理注释，未称其自由底物与供体已由该遗传实验确定；**YOL096C — COQ3 — CoQ 生物合成 O-甲基转移酶**，有缺失／回补线粒体酶活实验。系统 ID 与名称由 S9/S10 的物种整理记录对应，不能把当前序列等同于历史实验构建。W29 位点身份和 Y. phangngaensis 比较范围不在本独立来源审计内。

| 主张 | 判定 | 实际来源及定位 | 允许支持的范围与保留限制 |
|---|---|---|---|
| **C1** 酵母遗传缺陷／回补关联晚期 Q6 合成，不直接确定自由底物氧化态或供体。 | supported | S1：摘要；pp.2998–3000 的中间体鉴定与回补；pp.3001–3002 讨论。S9：ID／别名对应。 | coq7–1 点突变检出 DMQ6、缺少可检测 Q6，回补恢复呼吸及 Q6；完整删除株不积累该 DMQ。作者仍讨论多步骤功能或复合体作用，故不能把遗传结果升级为纯化酶对自由 DMQ／DMQH2 的直接底物实验。 |
| **C2** Lu 2013 人源 GB1-COQ7 的 NADH／短链体系证据，直接产物鉴定明确为 DMQ0，非 W29／Q9。 | supported | S3：`Hydroxylase Activity GB1-hCLK-1`、`Hydroxylation of DMQ derivatives`；Fig.3 图注；Table 2 与 Fig.4 相关正文；`Substrate-Mediated Electron Transfer`。 | 纯化人源 GB1 融合构建在 DMQ0＋NADH 中生成羟化产物，GC-MS 支持；DMQ2 为 NADH 消耗动力学补充。原文明确长链底物受溶解度限制而未使用。底物介导电子转移支持该体系机制，不直接证明 W29 的自由 Q9 底物／供体或膜内净计量。 |
| **C3** Nicoll 2024 的 5ox／5 对比为可测 NADH 消耗；联合产物不唯一确定自由中间体；COQ9 提效不建立普遍催化 AND。 | supported | S4：祖先构建段落；`C2 decoration`；Fig.5b–e 图注；`COQ7–COQ9 kinetics with NADH consumption assay`；`Small-scale reactions`。 | 四足动物祖先重建体系的该 NADH 测定仅对 5ox 检出消耗。5ox 经三蛋白联合过夜反应检出 CoQ1 与 CoQ1H2；中间体 6 因稳定性问题未单独制得。存在 COQ7 单独活性，加 COQ9 提高效率；不能据此主张所有物种／所有条件下催化都必须 AND，也不能以 NADH 阴性排除其他供体条件。 |
| **C4** 酵母线粒体 WT／KO／回补 Q3 类似物与 NADH 实验支持末端 O-甲基化，非纯化酵母单酶；终点氧化不证明直接释放醌。 | supported | S2：p.21667 `In Vitro Assays`；p.21668 Fig.3 与结果；p.21670 讨论。已目视 p.21667 Fig.1 和 p.21670。S10：系统 ID 对应。 | WT 线粒体有活性，缺失株无可检出活性，回补恢复；方法使用 NADH 以生成还原态底物，讨论将其解释为供给还原当量；未据此指定酵母还原酶。方法在结束后用铈(IV)试剂氧化甲基化产物再分析，因此 HPLC 检出 Q3 不等于酶直接释放醌。纯化且预还原底物实验属于细菌蛋白，不能移作纯化酵母酶证据。 |
| **C5** Houser 1977 大鼠 Q9 氢醌为甲基化直接底物；NADH 用于预还原，摘要还原活性线索没有指定基因；全文未知。 | supported（仅摘要范围） | S8：PubMed 原始研究摘要第 1–2 段；S5：纯 O-甲基化分类方程。 | 摘要直接支持大鼠肝线粒体 Q9 链长、底物还原与甲基化关系；也报告 Triton 溶解使前体还原活性失活及连二亚硫酸盐部分替代 NADH。已核摘要没有基因身份。对已还原底物的纯甲基转移不再额外固定计入 NADH，是结合分类反应的边界推论；不是对全通路无 NADH 成本的实测证明。 |
| **C6** 两种羟化 EC 定义分开：醌／NADH 与醌醇／泛型供体，不能拼接。 | supported | S6/S7：现行完整条目的 reaction、comments、版本标记。 | EC 1.14.13.253 为真核醌／NADH 表示；EC 1.14.99.60 现行评论指原核、醌醇及泛型还原供体。S6 的醌醇无活性评论特指脊椎动物酶，不能扩大为全部真核酶实验。分类也不单独确定 W29 的自由底物和供体。 |
| **C7** EC 2.1.1.64 为还原态末端甲基化，早期 COQ3 反应对应 EC 2.1.1.114；官方定义不构成新增独立实验。 | supported | S5：reaction、comments 及 references。 | 末端净式为 DMeQnH2＋SAM → QnH2＋SAH；条目将早期聚异戊二烯基二羟基苯甲酸甲基化指向 EC 2.1.1.114。这是分类整理，底层引文与 S2/S8 重叠，不能重复算作独立酶学重复。GEM 中质子／电荷表示仍须按所用代谢物身份配平。 |

覆盖：**7 total | 7 audited | 7 supported | 0 unresolved claims | 0 contradicted | 0 unchecked**。这里的“支持”针对上述带条件、带限度的文字；**W29 自由 Q9 底物、实际供体、膜内红氧连接和原生酶活仍未解决**，不得把审计覆盖率当作这些科学问题已闭合。C5 的全文方法／图表、C3 图中数值的独立重算以及现有候选 GPR 的生物学接受均不在本轮覆盖范围。

## 来源清单及读取范围

- **S1，原始研究全文缓存：** Marbois BN, Clarke CF (1996). *The COQ7 gene encodes a protein in Saccharomyces cerevisiae necessary for ubiquinone biosynthesis.* [DOI](https://doi.org/10.1074/jbc.271.6.2995)。读取摘要、相关结果和讨论；没有本次实验复现。
- **S2，原始研究全文缓存：** Poon WW et al. (1999). *Yeast and rat Coq3 and Escherichia coli UbiG polypeptides catalyze both O-methyltransferase steps in coenzyme Q biosynthesis.* [DOI](https://doi.org/10.1074/jbc.274.31.21665)。读取相关方法、Fig.3 结果、讨论及定位实验；目视已有 p.21667、p.21670 渲染页。
- **S3，原始研究全文缓存：** Lu TT et al. (2013). *Aging-associated enzyme human clock-1: substrate-mediated reduction of the diiron center for 5-demethoxyubiquinone hydroxylation.* [全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3615049/)，[DOI](https://doi.org/10.1021/bi301674p)。读取相关构建／纯化、GC-MS、NADH 与电子转移段落；未新增审核补充材料图。
- **S4，原始研究全文缓存：** Nicoll CR et al. (2024). *In vitro construction of the COQ metabolon unveils the molecular determinants of coenzyme Q biosynthesis.* [全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7615680/)，[DOI](https://doi.org/10.1038/s41929-023-01087-z)。读取祖先构建、相关结果、Fig.5 图注及方法；图示数值未独立重算。
- **S5–S7，现行官方分类：** [EC 2.1.1.64](https://iubmb.qmul.ac.uk/enzyme/EC2/1/1/64.html)、[EC 1.14.13.253](https://iubmb.qmul.ac.uk/enzyme/EC1/14/13/253.html)、[EC 1.14.99.60](https://iubmb.qmul.ac.uk/enzyme/EC1/14/99/60.html)。完整条目已重新打开；分别标记 modified 2011、created 2024、modified 2024。不是本次新增实验。
- **S8，原始研究摘要，非全文：** Houser RM, Olson RE (1977). *5-demethylubiquinone-9-methyltransferase from rat liver mitochondria. Characterization, localization, and solubilization.* [PMID 863914](https://pubmed.ncbi.nlm.nih.gov/863914/)。本轮重新核对文献身份和摘要 3 段；未核全文实验方法、图表和讨论。此前记录的全文访问失败不解释为全文不存在。
- **S9–S10，整理数据库缓存：** [YeastCyc CAT5／COQ7](https://pathway.yeastgenome.org/gene?id=MONOMER3O-164&orgid=YEAST)、[YeastCyc COQ3](https://pathway.yeastgenome.org/gene?id=YOL096C&orgid=YEAST)。仅作为名称、系统 ID 和已显示整理反应的对应；页面为 YeastCyc 22.5，2026-10-02 生成。不能把其反应图当成独立原生酶学测定。

## 终稿措辞复核

复核时间：2026-10-02T23:14:12.615052+00:00。已读取当前 `REPORT.md`，核对主表、方程说明及 C1–C7 定位。最终报告 SHA256：`3e2ac5c713a831081adc27861278303daa3988d8fe0eb0a636540b71e1c4360a`。

主表将酵母点突变结果收窄为“可检出 DMQ6”，避免误解为高于回补株；C3 仅将正文明确记载的无侧链形式归给可测活性 5ox，不自行补足还原态比较底物的链长。以上修改已核实。R695 两种现行 EC 的底物／供体成对区分，R385 的 NADH 预还原与纯甲基转移净式分开，联合产物、阴性 NADH 读数、COQ9 辅助作用与 W29 外推限制均保持本审计范围。**C1–C7 最终措辞核对通过，7/7 支持。**

本次终稿通过不扩大到 W29 序列／结构身份、模型实际计量或全篇所有声明的独立复核；这些由各自的身份与静态模型核验记录覆盖。人源 NCBI Gene 10229 的额外官方身份映射由主代理核实，未作为本审计新增酶学来源。
