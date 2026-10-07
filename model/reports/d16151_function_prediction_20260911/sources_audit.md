# YALI1D16151g 功能候选：独立生物学来源审计

核查日期：2026-09-11。限定只读范围：4篇指定参考酶原始论文、允许的1篇NpaA1补充论文，以及本目录已有UniProt/RCSB身份记录；在原始结果提交后追加九序列参考和六结构参考的静态审计。已读 input_identity.json；本审计未运行BLAST、结构比较、优化、预测或实验，未改变模型。当前结论只支持功能候选，不能据此接受新GPR。

目标 **YALI1D16151g — 文献称GLN2/GS2，统一正式名称未核实 — 谷氨酰胺合成酶家族蛋白，原生催化特异性未表征（uncharacterized；现有GS角色为model/GPR assignment only）**。固定目标为W29 AOW03998.1 / A0A1D8NED6，458 aa，序列SHA256 `11c383111f92c305adb675026cdbd0133b22bfdc5f02e5700bbdd131050ef123`。对应既有模型必须标记“AlphaFold预测”；参考实验不能转写为目标实验。

## 参考身份

| 系统ID / 蛋白accession | 名称及蛋白功能 | 证据等级和条件 |
|---|---|---|
| PA5508 / Q9HT65 | PauA7；偏好有机胺的γ-谷氨酰化酶 | P. aeruginosa PAO1参考蛋白有酶学和4HPP晶体学；UniProt条目虽仍unreviewed，也不能忽略原始实验 |
| SCO1613 / O88070 | GlnA4；γ-glutamylethanolamide synthetase | S. coelicolor M145的乙醇胺利用、纯化酶及产物检测：experimentally verified；本次所用结构为既有AlphaFold预测 |
| A9ZPH9；统一系统基因ID未核实 | GmaS/GMAS；γ-glutamylmethylamide synthetase，也可形成茶氨酸 | M. mays No.9纯化酶的甲胺/乙胺反应实验；accession序列身份另由UniProt所引2008克隆研究连接 |
| b1297 / P78061 | PuuA；γ-glutamylputrescine synthetase | E. coli K-12的腐胺利用及纯化酶实验：experimentally verified |
| STM4007 / P0A1P6 | glnA；经典谷氨酰胺合成酶I | Salmonella typhimurium LT2经典GS结构参考，PDB1FPY；非目标酵母蛋白 |
| 论文名npaA1；系统locus/accession本轮未核实 | NpaA1；1-萘胺γ-谷氨酰化酶 | Pseudomonas sp. JS3066原始研究；本次仅取得出版社/PMC检索正文，精确实验细节保留部分核验 |

## 原子声明与判定

| ID | 核查声明 | 已打开/取得的原始来源与定位 | 判定 | 可用范围与限制 |
|---|---|---|---|---|
| SA01 | PA5508参考蛋白偏好芳香有机胺，不能由GS样折叠自动判为经典GS | [Ladner et al.2012](https://pubs.acs.org/doi/10.1021/bi3014856)作者摘要；独立直接打开[RCSB4HPP](https://www.rcsb.org/structure/4HPP)的原始论文摘要与实验记录 | supported | 参考蛋白的功能反例；4HPP为2.50 Å晶体结构，含单环六聚体记录；不转移其底物到酵母 |
| SA02 | PA5508在原研究检测中没有GS活性 | [NIST作者机构摘要](https://www.nist.gov/publications/structure-and-activity-pa5508-novel-hexameric-glutamine-synthase-homolog)检索返回明确如此描述 | partial | ACS全文403；本轮未读完整活性检测方法/检出限或SI。因此可写“作者摘要报告未检出”，不能写全条件绝无GS活性 |
| SA03 | SCO1613的GlnA4是乙醇胺γ-谷氨酰化酶 | [2019 mBio原文](https://journals.asm.org/doi/10.1128/mbio.00326-19)，Abstract、Results酶学及Fig8图注 | supported | 原文有纯化酶、产物HPLC/ESI-MS、基因删除与回补；底物面板50 mM、n=3；氮源条件为该细菌，不是酵母 |
| SA04 | M. mays No.9纯化酶偏好甲胺，也接受乙胺；其氨活性低但不是零 | [2007原始PDF](https://www.jstage.jst.go.jp/article/bbb/71/2/71_60590/_pdf)，p.550、Substrate specificity及Table4文本 | supported | 表中在2mM比较底物下甲胺100%、乙胺75%、氨0.75%。本轮PDF文本可读但工具截图失败；数字不用于目标定量预测，不能把该参考写成完全不催化氨 |
| SA05 | 2007实验已直接验证A9ZPH9完整序列所编码蛋白的全部身份 | sources/A9ZPH9.json的references；2007为蛋白N端1–20测序，2008 DOI10.1271/bbb.70462才记录基因克隆及更多肽段 | partial | 2007支持纯化酶功能；完整accession联系依据数据库引用链。2008原文未在本轮授权的4+1篇范围重开，避免将年份/实验混写 |
| SA06 | b1297/PuuA能以ATP驱动腐胺γ-谷氨酰化 | [2008 JBC原论文的出版社摘要](https://www.sciencedirect.com/science/article/pii/S002192582059652X)及PubMed原始摘要检索正文 | supported | 摘要明确非融合纯化酶与基因删除证据；底物是腐胺，并非将腐胺本身合成出来 |
| SA07 | PuuA在所测有机胺中最适腐胺，且未形成谷氨酰胺 | 同一2008论文的作者PDF检索文本，“PuuA Activity for Other Amines and Ammonia”；[原文镜像](https://www.researchgate.net/publication/5353946_-Glutamylputrescine_Synthetase_in_the_Putrescine_Utilization_Pathway_of_Escherichia_coli_K-12) | partial | 出版社直连403，详细部分只取得原文检索文本；对比结果可保留为参考线索，未独立核验Fig5图像和完整检测限 |
| SA08 | 经典GS的E327/D50/Y179在PauA7中对应W296/G40/A147 | 独立打开2019 mBio原文“GlnA4 is similar…”结构比较段；PauA7是4HPP实验结构、GlnA4当时是同源模型 | supported | 是原文明确的参考结构比较，不是本轮对目标YALI1D16151g的对应位点验证；也不是三位点逐一突变验证 |
| SA09 | 1FPY残基编号与P0A1P6的UniProt编号相差起始Met一个位置 | 本目录sources/1FPY_entity.json的SIFTS：entity_beg_seq_id=1、ref_beg_seq_id=2、length=468；P0A1P6.json结构引用为2–469 | supported | 参考图示/PDB的E327/D50/Y179映射到该UniProt序列时应检查+1；禁止混用编号后直接宣称目标保守或缺失 |
| SA10 | 三位点化学变化支持氨口袋性质改变的解释 | SA08原文结构论证及残基侧链化学 | partial | E→W失去羧酸负电/酸碱性质但W更大；D→G去掉羧酸并减少侧链；Y→A去掉芳环/羟基并减小侧链。不能把三者都称“扩大口袋”，也不能凭此证明目标底物 |
| SA11 | NpaA1论文给出比全局折叠更直接的底物选择位点实验 | [2024 eLife正式版](https://elifesciences.org/articles/95555)及[PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11335346/)检索正文、Fig5图注与突变结果段 | partial | 检索所见M81/W235入口与V201口袋有突变活性研究；正式版日期2024-08-20。本轮直连未取得全文；M81A/M81W叙述疑有笔误，未用具体方向作结论 |
| SA12 | 这些参考能确认YALI1D16151g专一催化某一种有机胺 | 参考物种/序列均不是目标；目标尚无直接酶学或相关回补 | unverified | 只能使“有机胺γ-谷氨酰化酶”成为功能候选；甲胺、乙醇胺、腐胺和芳香胺不能互相替代为精确定名 |
| SA13 | 这些参考足以排除YALI1D16151g所有经典GS活性 | 同上；尤其SA04存在有机胺参考仍有低氨活性的反例 | unverified | 催化偏好、低水平副活性和足够支撑生长的GS容量是三种不同声明 |
| SA14 | 同源结构可以确认YALI1D16151g原生胞质定位、寡聚组成或当前OR替代性 | 本次外部来源均无目标原生定位/复合体/生长救援直接证据 | unverified | 均保持未决；当前研究是功能预测，不接受GPR或模型变更 |

## 来源匹配与机制边界

四类参考确实展示了GS样家族的底物多样性，但证据强度不同：2019乙醇胺原文与2007纯化酶PDF可直接读；2012有直接RCSB摘要/结构及作者机构摘要，详细酶学未读；2008直接出版社访问失败，摘要与作者原文检索文本可读。原文、搜索取得的原文片段和数据库引用链在上表分开，没有把检索成功等同全文验证。

2019结构段同时使用“E327 flap”和与D50相互作用的“E335”措辞，存在编号不一致。本审计不自行修正作者文本，也不据此补写精细质子传递机制。三位点只能作为带编号和对齐依据的候选判别特征；还需结合真实侧链位置、入口、口袋、辅因子与亚基邻接结构。整体结构分数高不能单独证明特定胺结合，缺失一个氨位点也不能单独证明完全无GS副活性。

NpaA1原文检索显示同一家族可由入口与口袋突变改变有机胺选择，因而支持“检查局部结构优于只看总折叠”的分析方法；不支持未经目标映射就搬用其编号。

## 后续提交的序列与结构结果审计

本节独立读取原始BLAST XML、六份US-align原始日志、PDB的ATOM/SEQADV/REMARK记录及现有JSON；属于交付物静态核实，不是本审计重新运行计算。六结构分数的参考归一化长度是有坐标的链，不能把不同分母的TM分数直接当同一指标。

| ID | 核查声明 | 原始来源和实际核查 | 判定 | 可用范围与限制 |
|---|---|---|---|---|
| SA15 | 六个结构比较的统计已按原始日志正确报告 | `results/D16151__*.txt`及`structure_results.json`：O88070_AF对齐437、RMSD1.56Å、TM目标/参考0.92122/0.91345；P78061_AF为435、1.97Å、0.89955/0.87390；1FPY为426、2.34Å、0.86353/0.84602；4HPP为403、2.73Å、0.80971/0.86716；9QUR为421、2.07Å、0.86404/0.88620；3FKY为325、3.48Å、0.61480/0.81199 | supported | 数字核实，不等于底物或物种功能验证；目标及O88070/P78061参考均为AlphaFold预测。几何复算是主任务新增交付物，本审计未独立执行其全量距离重算 |
| SA16 | 直接1FPY映射具有金属配位残基保守及若干氨/谷氨酸位点差异；GlnA4对应位点亦有差异 | `structure_results.json`与原始日志：1FPY E129/E131/E212/E220/H269/E357→目标E135/E137/E197/E204/H253/E354；E327/D50/Y179/N264→W322/C45/G164/P248；O88070 F331/E269→W322/P248 | supported | 对应关系是此次结构对齐的观察。目标W322的pLDDT67.19，显著弱于全局均值；不能由Cα距离宣称侧链几何、口袋尺寸或催化能力已验证 |
| SA17 | “4HPP W296在晶体中未建模” | 原始`sources/4HPP.pdb`含`ATOM 2234 CA TRP A 296`，坐标46.515,13.398,−12.201；其B因子164.61；REMARK465的A链缺失为365–381 | contradicted（已纠正） | W296有坐标，只是没有进入本次对齐的残基配对；实验B因子不是pLDDT。已通知主任务，后者确认修正；保留这一被推翻的中间声明以追踪审计修正，不计作当前最终结论仍然错误 |
| SA18 | 实验结构与当前UniProt序列并非全部逐位相同 | 1FPY原始PDB的author391是P；P0A1P6 UniProt392为A且明确记录A→P sequence conflict。3FKY的SEQADV列出链A的251A对当前参考251T、264T对当前参考264M，REMARK999说明序列差异 | supported | 可称参考结构并非当前参考序列完全一致；不能无证据把这些差异全部叫人工工程突变。编号偏移和既有序列冲突已保留 |
| SA19 | 在九个固定序列参考中，GlnA4是目标的首位BLAST命中 | `sources/local_blast.xml`原始HSP和`local_blast_results.json`一致：O88070 q19–456/s41–460、135/443=30.47%、qcov95.63%、E9.23111e−54；P78061为26.37%/96.94%/1.65074e−37；A9ZPH9为25.16%/95.85%/9.55096e−35；O88070 F331→W322在序列与结构对齐中一致 | supported | 九序列面板、八条过阈值命中；不是九个都命中。E值针对小面板，不能作为全库E值。远程全库任务未返回结果，不能由此宣布全库最佳或已排除更近参考 |
| SA20 | 广义GlnA4-like有机胺γ-谷氨酰连接酶是合理功能候选 | SA01–SA11原始参考功能，以及SA15/SA16/SA19的覆盖、折叠和位点观察 | partial | 支持中等强度的间接候选，须写“基于AlphaFold预测的功能候选”；保守金属位点主要支持家族催化框架。经典GS结构同样有较高全局相似度，故整体分数本身没有证明不同底物 |
| SA21 | 乙醇胺可列为具体底物的优先候选 | O88070/GlnA4的直接参考酶学，及本面板中较优序列/结构匹配 | partial | 可以优先测试，低至中等只是未经校准的定性排序。约30%序列一致度、E269→P248及低置信度W322限制底物转移；尚无目标乙醇胺结合/产物/动力学/生长实验，不能写成已确认催化乙醇胺 |
| SA22 | 可据此确定目标是GlnA4正交同源物、赋具体EC或原生乙醇胺降解通路 | 当前有限同源/结构参考，无全库完成结果、系统发育和目标原生实验 | unverified | 不能迁移细菌原生通路和定位。精確底物、GS副活性、体内作用与GPR仍待核实 |

局部对应不是传递性事实：经典GS参考经PauA7转到目标，与经典GS直接对到目标时可能得到不同匹配。实际4HPP G40→目标G43，A147→目标F163；不能把目标描述为与PauA7三个特征位点均精确相同。GlnA4的F331与目标W322虽同为芳香残基，也不意味着同一口袋尺寸或相同胺偏好。

收尾时已读取最终`REPORT.md`，W296错误已纠正，有限面板、全库WAITING、AlphaFold参考标签、局部低置信度及E269→P248差异均已写入。八条命中的目标/参考覆盖率与JSON一致，GlnA4的序列覆盖95.63%和结构配对覆盖95.41%使用不同分子，报告未混写。六组配对数403+426+325+421+435+437=2,447，与`geometry_verification.json`各组的数量及通过标记一致；本审计只确认该QA交付物的记录，不将其称作本审计重新进行了全部几何复算。

## 覆盖率与停止边界

```text
total claims 22 | audited 22 | supported 10 | unresolved 11 | contradicted 1 | unchecked 0
unresolved = partial 7 + unverified 4
audit coverage = 22/22 = 100%; directly supported claim share = 10/22 = 45.5%
contradicted 1 = SA17，被纠正的中间声明；不存在以多代理一致代替来源验证的计数
```

audited表示每项已进行来源匹配或明确缺口判定，并不表示22项全部科学结论成立；7项partial的原文/方法/身份/外推限制不能写成已完全核实。六参考的全量距离重算、软件可靠性和全库检索尚未完成部分均不在此分母内。本审计未下载大文件，未扩展额外论文，未运行科研计算，未修改任何模型。可交付限定候选判断；目标精确底物、GS副活性、定位及GPR仍需独立证据。
