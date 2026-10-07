# 酿酒酵母末端 CoQ 氧化态：原始来源调查

2026-10-02；本子任务仅查阅来源并归档，不修改模型、GPR或代码，不运行GEM、序列比对或结构预测。中立问题：R18生成DMQ9H2后，是否有独立氧化酶、连续醌醇羟化路线，或只能采用待验的汇总表达？允许目标是查明反应及证据边界，不能用“必须闭合”选择证据。

## 结论与可用反应

酿酒酵母提供了末端醌醇甲基化的直接实验依据，以及连续醌醇路线的数据库表达；本次没有找到可直接赋予GPR的DMQH2独立氧化反应。将氧化态缺口径直命名为“缺少一个独立氧化酶”仍超出证据。

按本项目现有中性NADH/NAD和SAM/SAH分子式记账：

| 候选反应 | 对应位置 | 证据状态与决定 |
|---|---|---|
| DMeQ9H2 + SAM → CoQ9H2 + SAH | 修订R385 | 直接酿酒酵母线粒体实验支持甲基化；短链Q3模拟底物。原生Yarrowia的C9底物/酶仍为跨物种候选。 |
| DMQ9H2 + NADH + O2 → DMeQ9H2 + NAD + H2O | 连续醌醇R695候选 | 原子守恒；YeastCyc只给generic donor，Reactome为同源计算推断；不能称酿酒酵母完整净反应实测。与2024真核COQ7醌型酶学存在需解释的差别。 |
| DMQ9H2 + O2 → DMeQ9 + H2O | 机制衍生R695候选 | 原子守恒；需假设自由醌醇进入酶后直接还原双铁，并被同一催化循环羟化。Lu2013未直接检验这个输入/输出；不是可直接实施的已证实反应。 |
| DMQ9H2 + acceptor(ox) → DMQ9 + acceptor(red) | 独立氧化连接 | 受体身份、完整计量和催化基因未定；不能将空泛acceptor作为已找到的反应。 |

其中DMQ9H2=C53H82O3，DMQ9=C53H80O3，DMeQ9H2=C53H82O4，DMeQ9=C53H80O4，CoQ9H2=C54H84O4。电子记账不能自动补足酶学。连续醌醇净反应也可能掩盖内部氧化/还原交换；它与拆分机制并非同一证据声明。

## 基因身份

| 系统ID | 已核实名称 | 简要功能、证据等级 | 本问题的角色 |
|---|---|---|---|
| YOL096C | COQ3 | CoQ环O-甲基转移酶；experimentally verified（酿酒酵母线粒体活性/基因缺失及恢复实验） | 支持末端醌醇甲基化；不是独立DMQ氧化酶。ID由当前YeastCyc页面核对。 |
| YOR125C | CAT5，别名COQ7 | CoQ末端羟化相关蛋白；curated annotation；遗传学功能已验证，具体自由醌醇净反应本次未找到直接实验 | R695同源酵母参照；不能把数据库底物形式等同于纯化酵母酶实测。ID由当前YeastCyc页面核对。 |
| YLR290C | COQ11 | CoQ合成复合体相关调节蛋白；experimentally verified（参与CoQ合成/调节），具体催化活性uncharacterized | 2020文献提出成熟Q6H2氧化假说，未测DMQ底物，不能给本连接GPR。ID及名称对应见Allan2015原始发现。 |

## 原始文献定位

1. **Poon et al. 1999, JBC 274:21665–21672, DOI [10.1074/jbc.274.31.21665](https://doi.org/10.1074/jbc.274.31.21665)**。作者PDF：[公开全文](https://www.biochemistry.ucla.edu/Faculty/CClarke/pdf/21665.pdf)。本地Poon1999.pdf和.txt，页面3和6截图已目视核对。p21667 Methods加入3mM NADH，使酵母线粒体试验中demethyl-Q3处于hydroquinone状态；终止时用cerium(IV)氧化产物。Fig1结构明确醌醇→醌醇；Fig3比较WT、coq3缺失和恢复；p21670 Discussion说明酵母试验甲基化依赖NADH，解释为还原底物。支持Coq3步骤及氧化型终点检测不能直接代表酶释放氧化型产物。没有确定独立DMeQ还原酶，不能由此把NADH写进纯Coq3催化反应。

2. **Marbois & Clarke 1996, DOI [10.1074/jbc.271.6.2995](https://doi.org/10.1074/jbc.271.6.2995)**。本地Marbois1996.pdf/.txt。p2995摘要、p3000–3002 Discussion：coq7点突变积累DMQ，缺失只见早期中间体；恢复COQ7恢复Q合成及呼吸。这是代谢遗传学定位，不是纯化酵母Coq7与醌/醌醇的并列底物试验。Fig1本身是当时的拟议路径，不能用其化学画法确认所有氧化态及辅因子。

3. **Padilla et al. 2004, DOI [10.1074/jbc.M400001200](https://doi.org/10.1074/jbc.M400001200)**。作者上传的[公开全文](https://www.researchgate.net/publication/8623858_Demethoxy-Q_An_Intermediate_of_Coenzyme_Q_Biosynthesis_Fails_to_Support_Respiration_in_Saccharomyces_cerevisiae_and_Lacks_Antioxidant_Activity)可查；本地Padilla2004_web_excerpt.txt为web工具提取的Discussion/References片段，非PDF原件。TableII/III及p26002–26003：只有DMQ6的酵母株不能支持两类至cytochrome c的呼吸活性；但bc1组装及组分量受损，即时加入Q3和DMQ3均未救回。作者明确保留DMQ6是否为bc1底物的未决结论。它反对直接假定可呼吸替代，但不是成熟完整bc1对DMQH2的纯化底物阴性试验。不能拿此文直接支持或排除本项目独立DMQ氧化GPR。

4. **Allan et al. 2015, [PMC4367260](https://pmc.ncbi.nlm.nih.gov/articles/PMC4367260/)**。本地Coq11_2015.html/.txt。Results及Discussion确认YLR290C与CoQ复合体和合成效率的关系；SDR/Rossmann特征只是功能线索，文中还提出不同催化假说。没有DMQH2氧化底物测定。

5. **Bradley et al. 2020, [PMC7196636](https://pmc.ncbi.nlm.nih.gov/articles/PMC7196636/)**。本地Coq11_2020.html/.txt，Discussion文本2615–2695。根据缺失表型、抗氧化保护及SDR特征，提出Coq11可能参与Q6H2→Q6氧化；是作者假说，针对成熟Q池。没有确定电子受体，没有纯化Coq11催化读数，不能转写为DMQH2氧化功能。

6. **Yang et al. 2011, [PMC3156029](https://pmc.ncbi.nlm.nih.gov/articles/PMC3156029/)**。本地Yang2011.html/.txt，Methods2.3、Results3.3、Discussion4.2。线虫线粒体提取/补加试验指向DMQ9影响复合体I至Q的传递；使用含多种脂质/醌的提取物，缺少纯化DMQ9。含量测定还先加p-benzoquinone使醌全部氧化。不是酿酒酵母独立DMQ醌醇氧化酶证据。用于追溯Padilla2004及防止把外源Q、RQ和代谢中间体效应混同。

## 对已归档机制文献的重新核对

复用并重读`artifacts/coq_literature_revision_20261002/source_audit/sources/Lu2013.txt`的Substrate-Mediated Electron Transfer/Scheme1及方法。Lu2013 DOI10.1021/bi301674p实际输入是氧化型DMQ0/DMQ2与NADH；双指数动力学和EPR结果支持瞬态结合型DMQred向双铁传电子。论文未直接测自由DMQH2+O2→羟化产物的周转。因此后一净反应只能是机制衍生候选。

Nicoll2024 DOI10.1038/s41929-023-01087-z的Fig5 NADH消耗读数，不足以排除所有不消耗NADH的醌醇周转；需要独立测产物。官方EC1.14.13.253的概括不能使这项试验凭空变成完整产物阴性对照。本条和既有来源由根任务的独立审核进一步核对。

## 数据库差别

当前YeastCyc version22.5：RXN3O-75是DMQ6H2 + reduced acceptor + O2 → DMeQ6H2 + oxidized acceptor + water；CAT5/COQ7对应YOR125C。COQ3对应YOL096C，末端反应DMeQ6H2+SAM→Q6H2+SAH+H+（数据库质子化形式）。反应页面保留EC1.14.99.60，列Marbois1996、Tran2006等历史引用。2024官方EC分流把真核醌型COQ7列为1.14.13.253、1.14.99.60定义为原核醌醇反应。因此须把它作为可比较的酵母数据库候选，不能直接当作实验证明。

## 检索边界及交付

本轮6批定向网络搜索，共18个检索式，集中酿酒酵母COQ3/COQ7、DMQ氧化、呼吸复合体III、COQ11；沿参考文献追溯原始研究。新增归档2篇作者PDF、3篇PMC全文、3个YeastCyc记录和1个Padilla作者上传网页提取片段。复用Lu2013/先前COQ7来源；未扩展蛋白候选搜索。全文或截图获取失败已记录/可由工具记录查明；Padilla缺PDF图表目视检查，结论只依赖已打开的正文段落。来源发现/摘要笔记不等于独立审核完成；根任务另行维护独立claims审核，不在此虚报审核覆盖率。文件完整SHA见MANIFEST.json、URL及获取日期见retrieval.json。
