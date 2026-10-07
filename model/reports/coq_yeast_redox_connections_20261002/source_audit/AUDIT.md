# DMQ 醌醇至成熟 CoQ 连接：独立来源审核

核验日：2026-10-02；审核者：独立 source-audit 子代理。范围：定向读取已指定原始文献、官方酶学分类、YeastCyc 和两份本地 Yeast-GEM XML；不运行生长、FBA、结构预测，不修改模型或 GPR。六个定向联网批次后停止联网；复用缓存文件并独立打开原文，不以其他代理结论作为证据。输入身份见 `MANIFEST.json`。

## 身份及术语

| 系统 ID | 已核实名 / 功能 | 本次证据等级 |
|---|---|---|
| S. cerevisiae YOR125C | CAT5 / COQ7；CoQ 生物合成晚期羟化相关二铁蛋白 | 酵母功能遗传证据与整理注释；本次没有确认纯化酵母酶的醌/醌醇偏好 |
| S. cerevisiae YOL096C | COQ3；SAM 依赖 O-甲基转移酶，参与首末两次 O-甲基化 | 酵母线粒体分离物、缺失/回补与定位实验；不是纯化酵母酶红氧底物比较 |
| Y. lipolytica YALI1B20835g | COQ3家族/O-甲基转移酶候选；原生正式功能名未实验确证 | 既有版本化序列与跨物种家族证据；未重新进行比对、定位或酶活试验 |
| S. cerevisiae YLR290C | COQ11；CoQ 合成相关、SDR 样蛋白，精确催化反应未定 | 遗传、共纯化与定位实验；DMQH2 氧化活性未在所读研究中证实 |

DMQ 表示 demethoxy-Q；DMeQ 表示 demethyl-Q；H2 表示相应醌醇。二者不得混用。酿酒酵母 Q6 与本项目 Q9 的侧链长度不同。

## 最终声明台账（对应 ../CLAIMS.md）

按 2026-10-02 本次收到的 C01–C08 原文逐项审核。`supported` 表示限定后的声明得到来源支持；`unresolved` 表示其中需要回答的生物学反应身份或适用性仍未确立。后两项不因其谨慎措辞正确而计为反应已获支持。

| ID | 裁决 | 独立打开的证据与定位 | 限制与处理 |
|---|---|---|---|
| C01 | supported | YeastCyc RXN3O-75 HTML image-map 与 CAT5/COQ7 基因页：DMQ6H2、Donor-H2、O2、DMeQ6H2、Acceptor、H2O；数据库版本22.5。 | 供体是泛型，不是已确定的 NADH；数据库表示不等于纯化酵母酶红氧底物实验，页面仍挂旧 EC1.14.99.60。 |
| C02 | supported | Poon1999 pp21667–21670；Fig1、Fig3、Fig5及方法/讨论。审核者打开作者PDF并目视p21667、p21670；文本独立核对p21668缺失、回补和去NADH对照。 | 酵母为线粒体分离物，短链Q3类似物；纯化预还原底物实验属于细菌UbiG。Ce(IV)处理后产物分析不能单独认定反应瞬间的自由产物红氧态。保留原声明限定。 |
| C03 | supported | YeastCyc COQ3页HTML明确链接RXN3O-102及醌醇方程；IUBMB EC2.1.1.64正式定义；既有coq_gene_gpr_review.tsv候选身份行。 | 支持反应家族与醌醇表示，不能把Yarrowia候选同源关系升级为原生酶活。具体H+须随本项目微物种约定重新配平。 |
| C04 | supported | IUBMB EC1.14.13.253（2024创建）与EC1.14.99.60（2024修改）；Nicoll2024正文C2 decoration、Fig5图注及COQ7–COQ9 NADH consumption方法；Lu2013。 | 新分类明确真核quinone/NADH与原核quinol。无quinol活性评论指vertebrate；Nicoll为四足动物祖先蛋白，其还原态对照端点是NADH消耗。Fig5 d/e另有从5ox出发联合酶GC/MS产物，不能误写整图只有NADH测量。 |
| C05 | unresolved | Lu2013 DOI10.1021/bi301674p；Substrate-Mediated Electron Transfer、Scheme1、结论。实验是人源GB1融合蛋白及DMQ0/2，暂态还原态属于机制推断。 | 本次范围内未证实独立自由DMQH2氧化步骤、确定酵母基因及电子受体。DMQH2+O2→DMeQ+H2O可配平但未在此被酶学验证；enzyme-bound shuttle不能自动拆为自由代谢物反应。不是“酶不存在”。 |
| C06 | supported | Padilla2004作者上传全文Discussion p26002；Bradley2020 Discussion的SDR/红氧假说。 | Padilla明确bc1装配/稳定性缺陷使DMQ6底物问题无法判别，急性加入Q3也未恢复。Bradley提出成熟Q6H2氧化的假说，未测DMQH2催化/受体。既不能确认复制成熟Q反应，也不能宣称成熟纯酶绝对排斥DMQH2。 |
| C07 | supported | 审核者重新解析固定Yeast-GEM9.1.1及9.0.2 XML；9.1.1扫描2748个物种并独立求元素/电荷账目。详见C07_INDEPENDENT_CHECK.json。 | r_0963为醌→醌；r_0532为DMeQ6醌+SAM+H+→Q6H2+SAH，元素残差均0、电荷残差−2。三种目标中间醌醇正确分子式扫描无命中。限定为这些文件的表示事实，不外推所有模型或生物学缺失。 |
| C08 | unresolved | C01–C07及固定项目候选晚期反应表示。 | 所列连续Q9路线可以作为待审提案，但R695自由DMQ9H2适用性、实际供体、产物红氧态未闭合；没有原生Q9生化验证，也没有本轮模型实施/通量闭合。维持候选和审批边界。 |

覆盖：**8 项总声明 / 8 项已审核 / 6 项 supported / 2 项 unresolved / 0 项 contradicted / 0 项 unchecked**。没有要求修正文献事实；用户答复必须保留C05、C08未决范围及C02/C04/C06的系统与端点限定。

## 可用于用户答复的结论

可以明确给出 Coq7 羟化、Coq3 末端 O-甲基化这两类对应反应；其中 YeastCyc 的醌醇连续表示是可参考的候选。不能据此直接把 R695 改成醌醇底物，更不能凭 Lu2013 或成熟 QH2 的 complex III 反应编造独立 DMQH2 氧化桥。末端醌醇甲基化有相对较强的酵母实验支持，而通往该步的红氧衔接及原生酶学仍未闭合。

## 来源与访问限制

- [Poon1999 作者 PDF](https://www.biochemistry.ucla.edu/Faculty/CClarke/pdf/21665.pdf)：独立打开全文，并检查本地原页图。
- [Lu2013 原始研究全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3615049/)：独立读取已存 HTML/文本，复用缓存。
- [Nicoll2024 原始研究全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7615680/)：本次打开全文并复查缓存方法与图注；Fig5 图像下载失败，因此不作本次图形数据点视觉核验声明。
- [IUBMB EC1.14.13.253](https://iubmb.qmul.ac.uk/enzyme/EC1/14/13/253.html)、[EC1.14.99.60](https://iubmb.qmul.ac.uk/enzyme/EC1/14/99/60.html)：官方分类。
- [YeastCyc CAT5/COQ7](https://pathway.yeastgenome.org/gene?id=MONOMER3O-164&orgid=YEAST)、[COQ3](https://pathway.yeastgenome.org/gene?id=YOL096C&orgid=YEAST)：工具直接打开失败，但主代理本轮取得的原始 HTML 已由审核者独立读取，包含反应、系统 ID 与版本。
- [Padilla2004 作者上传全文](https://www.researchgate.net/publication/8623858_Demethoxy-Q_An_Intermediate_of_Coenzyme_Q_Biosynthesis_Fails_to_Support_Respiration_in_Saccharomyces_cerevisiae_and_Lacks_Antioxidant_Activity)：本次成功打开正文；出版社入口和猜测的作者 PDF 路径失败。仅用正文论述，不从重排表格提取数值。
- [Allan2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4367260/)、[Bradley2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7196636/)：主代理本轮取得的原始全文缓存已独立读取。

本报告不接受任何科学模型变更；结论只涵盖上述来源和范围，不是穷尽性的“酵母中不存在该酶”。
