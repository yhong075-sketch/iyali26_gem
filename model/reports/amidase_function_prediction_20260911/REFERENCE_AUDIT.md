# 三目标功能预测：独立参考来源审计

日期：2026-09-11。直接读取本轮指定七条参考蛋白的 UniProt JSON、公开原文/摘要及上一轮已审 YlAMD1 原文。第七条 FAAH 在有限结构面板中构成竞争解释，因此经主代理明确加入本次范围。参考证据与三个 W29 目标的功能推断分开。仅读取、来源保存和静态核查；审计代理没有新增 BLAST、结构比对、预测、代谢求解或集群作业。原始来源、完整 SHA 与获取限制见 `sources/audit_retrieval.json`。

## 可用于比较的参考，及不能传递的结论

| 参考蛋白/系统ID | 已核实名称及功能 | 苯乙酰胺（PAM）是否实际有证据 | 适用边界 |
|---|---|---|---|
| Aspergillus nidulans AN8777 / P08158 | **amdS / AmdS**；乙酰胺酶（人工审阅注释：利用乙酰胺作碳/氮源） | 本次所读条目与原序列论文摘要没有直接 PAM 实验 | 可作已知乙酰胺功能参考；不能按泛 EC 3.5.1.4 升级为所有单羧酸酰胺实测阳性。1987 PMID3036667 核心是基因序列及5′突变，不是本次已打开的全底物面板 |
| Rhodococcus sp. N-771 ami / Q7DKE4 | **ami / RhAmidase**；体外水解苯甲酰胺、丙酰胺、乙酰胺、较弱丙烯酰胺（原研究酶学） | 本次已读原摘要和条目没有 PAM 实测；不能写 PAM 阴性或阳性 | 正确 DOI 是 **10.1016/j.bbapap.2009.10.001**，不是 .002。芳香阳性底物 benzamide 是 Ph–CO–NH₂，不是 Ph–CH₂–CO–NH₂ |
| Pseudomonas putida ATCC12633 mdlY / Q84DC4 | **mdlY / MdlY，mandelamide hydrolase（MAH）**；扁桃酰胺/苯乙酰胺水解酶（体外直接酶学） | **是**；PMID15196015 原摘要明确 PAM 为所测最佳底物；2017 原文延续并新增对硝基PAM实验 | 芳环在酰基侧；α位羟基、酰基链长和 N侧离去基团体积影响活性。2004/2009/2017 数字不可未经来源区分混合 |
| Arabidopsis thaliana At1g08980 / Q9FR37 | **AMI1**；IAM/PAM 酰胺水解酶（体外酶学；植物定位有该物种证据） | **是**；2014 原文 Fig.5及§2.3/§3.7包含 AtAMI1 与 PAM 比较 | PAM 相对 IAM 转化较高仅限10 mM底物、pH7.5、30°C、4–6 h 条件。不是 PAM 的 kcat/Km 排名，也不证明植物体内主要底物，更不能转移酵母定位 |
| Saccharomyces cerevisiae YDR242W / P22580 | **AMD2**；推定 amidase（人工审阅条目但具体功能未由本次来源实测确认） | 本次没有直接 PAM 实验证据 | 原始文献标题即 putative amidase；不能把 reviewed 标签当作底物实验已完成。按不明具体底物的同源参考处理 |
| Fusarium verticillioides FVEG_08289 / P9WEP6 | **AMD1**；FDB1 解毒基因簇的 amidase 候选（底物未决） | 本次原文没有该蛋白的 PAM 或 IAM 酶学 | 名称是 **AMD1，非植物 AMI1**。BOA诱导及基因簇位置不是特定酶学；删除仍保留 BOA 降解，作者讨论可能功能冗余。最强同源命中不能据此指定目标底物 |
| Rattus norvegicus P97612 | **Faah / Faah1**；脂肪酸酰胺水解酶1（oleamide、anandamide等脂肪酰胺有原研究支持的人工审阅注释） | 本次来源没有 PAM 底物酶学 | 1MT5 是结合 arachidonyl inhibitor 的2.8 Å晶体，使用原生37–573残基（缺N36/C6）；原生预测膜螺旋9–29没有进入比较。共享 amidase-signature 核心不能转移膜定位、长链底物特异性或神经信号生理角色 |

## 必须分开的底物结构

“芳香酰胺”不是可互换的底物类别。

| 化合物 | 连接方式 | 苯环在酰胺哪一侧 | 审计依据 |
|---|---|---|---|
| 苯甲酰胺 benzamide | Ph–CO–NH₂ | 酰基侧；苯环直接接羰基 | RhAmidase 原研究底物与benzamide晶体配体；PAM多一CH₂，不能等同 |
| 苯乙酰胺 phenylacetamide / PAM | Ph–CH₂–CO–NH₂ | 酰基侧；苯环与羰基间一个CH₂ | MdlY2017 PDF p9 compound3 与 p16 Scheme4 |
| 扁桃酰胺 mandelamide | Ph–CH(OH)–CO–NH₂ | 酰基侧；PAM的α位羟基取代 | 同文p9 compounds1/2 |
| 4-nitrophenylacetamide | p-NO₂–Ph–CH₂–CO–NH₂ | 酰基侧 | 同文p9 compound4；释氨、¹H NMR及直接光谱检测支持水解 |
| 4-nitroacetanilide | CH₃–CO–NH–Ph–NO₂ | **N侧/离去基团侧** | 同文p5试剂、p16 compound10“无可感知反应”；不能误写成4-nitrophenylacetamide阴性 |
| N-(2-hydroxyphenyl)malonamic acid | HOOC–CH₂–CO–NH–Ph–OH | **N侧/离去基团侧** | Fusarium 解毒路径产物；不是PAM，而且原文没有确定AMD1直接催化此产物的哪一步 |

MdlY2017提出部分芳基酯可能反向占据结合位点、形成 inverse acyl-enzyme；这是酶学加同源模型支持的机制解释，非目标W29蛋白已证明的结合模式。该文的MdlY模型基于FAAH模板构建、手工放置配体，不应称为MdlY–PAM实验共晶。

## 原子声明与覆盖

`source claims 21 | audited 21 | supported 16 | unresolved 5 | contradicted 0 | unchecked 0`

支持只对应条目写出的证据层级；5项未决不等于蛋白没有该活性。

| ID | 精确声明 | 来源定位 | 判定 |
|---|---|---|---|
| R1 | P08158是人工审阅乙酰胺酶参考，条目称支持乙酰胺碳/氮利用 | `P08158.json` FUNCTION/CATALYTIC ACTIVITY；PMID3036667摘要 | supported，curated annotation；本次没有重验其酶学 |
| R2 | P08158已经实测PAM活性 | 同上所读来源 | unverified；不能使用该更强声明 |
| R3 | RhAmidase来源PMID19819352的正确DOI以.001结束 | `Q7DKE4.json` references；`audit_reference_pubmed.json` 对应记录 | supported；已纠正任务提示的.002 |
| R4 | RhAmidase实测benzamide、propanamide、acetamide、acrylamide | 同PMID原摘要；Q7DKE4 FUNCTION/kinetics | supported；摘要报告kcat/Km依次约153.5、4.54、1.14、0.087 mM⁻¹s⁻¹，未以这些不同条件数字给目标排序 |
| R5 | RhAmidase在该来源已实测PAM阳性/阴性 | 同上 | unverified；摘要与条目未记录，不推断未打开全文没有其他实验 |
| R6 | MdlY存在PAM直接水解证据 | PMID15196015摘要；`audit_MdlY_2017.pdf` §3.1/p9、Table1/p11 | supported；2017 Table1 compound3值引用reference9（Wang等2009），不是2017重新测定的PAM数值 |
| R7 | MdlY的酰基侧与离去基团侧存在底物区分 | 2017 p9–11、p16；2004摘要 | supported；4-nitroPAM有反应，4-nitroacetanilide无可感知反应，不能把名称混淆 |
| R8 | 2017 MdlY配体口袋图来自同源模型而非该复合物晶体 | 2017 §2方法及§3.2，p8–9、p17–19 | supported；人工放置配体/FAAH模板的假设需保留 |
| R9 | AMI1的IAM转化及Ser137重要性有原研究支持 | PMID17555521原摘要；Q9FR37引用 | supported；不以2007摘要中的“specific”否定后续PAM结果 |
| R10 | AtAMI1本人也在2014 PAM/IAM比较中被检测 | `audit_PMC4844348.xml` §2.3、Fig5完整图注、§3.7 | supported；不只是另四植物同源物的类推；原文列AtAMI1 IAM最大比活3070 pkat/mg |
| R11 | 该2014比较不能直接转成PAM kcat/Km或体内偏好 | 同文等浓度相对比较、释氨方法和阴性对照 | supported；10 mM/5μg/0.3mL、4–6h，IAM=100%，n≥3；正文标注±SE及背景扣除 |
| R12 | ScAMD2有可用于指定底物的直接酶学阳性 | P22580；PMID2263500索引/题名及公开检索文本 | unverified；未取得可逐页审阅的原文；现有可审资料仅足以保留putative |
| R13 | P9WEP6应称Fusarium AMD1，而非Arabidopsis AMI1 | `P9WEP6.json` geneName；`audit_PMC4726666.xml` 正文AMD1/FVEG_08289 | supported |
| R14 | Fusarium AMD1敲除仍可BOA降解，原文保留冗余解释 | 同文摘要、Results目标删除段及Discussion功能冗余段 | supported；WT菌株FRC M-3125的删除与另一背景cosmid转座插入结论不同，不能抹去背景差异 |
| R15 | Fusarium AMD1已经有PAM/IAM或明确芳香酰胺底物酶学 | 同原文及UniProt关于角色未识别的说明 | unverified；基因簇功能不能指定该蛋白的底物 |
| R16 | YlAMD1有YB-392主要乙酰胺利用遗传学支持，存在残余生长 | 上轮`PMC7003347.xml` Par9–13/18及已目检Fig1/4；DOI10.1186/s12934-020-1292-9 | supported；W29相关序列与论文YB-392克隆桥接限制沿用上轮，不能转移为PAM酶学 |
| R17 | 上述参考的相似性足以确认三个W29目标具体底物/区室/GPR | 七参考各自证据及本轮比较设计 | unsupported；相似性只能支撑候选，特别不能从功能未明的最高命中指定底物 |
| R18 | benzamide/PAM/acetanilide不能按“含苯环”合并 | 2017目视compound结构、Scheme4与本表连接式 | supported；酰基侧和离去基团侧须分别核对 |
| R19 | Rat FAAH有oleamide、anandamide等脂肪酰胺水解证据的人工审阅注释 | `P97612.json` FUNCTION/CATALYTIC ACTIVITY及对应PubMed证据码；PMID12459591摘要 | supported；未将各条目引用全部全文重读，也没有为本次打开PAM实验 |
| R20 | 1MT5比较的是原生37–573残基，且为arachidonyl inhibitor复合物 | `1MT5.pdb` DBREF/坐标、`1MT5_entity.json`、`P97612.json`完整579aa、`audit_FAAH_pubmed.json`原摘要 | supported；参考TM分母537是构建体长度，不是完整579aa；N36和C6缺失 |
| R21 | 本次FAAH核心结构相似没有检验目标膜锚或膜定位 | P97612跨膜预测9–29，1MT5缺失1–36；本次单链Cα比较范围 | supported的推断边界；原生大鼠定位的物种实验证据不能转移W29；不得把高TM称为底物选择或定位实验 |

## 文献细节、冲突及未查部分

- [RhAmidase原研究](https://doi.org/10.1016/j.bbapap.2009.10.001)：本次读原摘要及人工审阅数据库；未取得全文。晶体野生型及失活S195A/benzamide分别为2.17/2.32 Å。不能把野生型游离结构当已结合PAM的复合物。
- [MdlY2004](https://doi.org/10.1021/bi049907q)：打开摘要。与[2017原研究](https://doi.org/10.1016/j.abb.2017.01.010)的作者机构公开稿联合审查；2017 p9、10、11、16已目视，全部文字已提取用于定向查阅。2017表1 PAM的Km=0.9μM、kcat=46s⁻¹注明来自2009 reference9；Q84DC4记录2004 PAM Km=3.8μM，两者不是同次测定，不静默替换。当前任务无需选定一种作为目标参数。
- [AMI1 2014研究](https://doi.org/10.3390/plants3030324)，PMID27135507的出版年为2014。完整XML可读；Fig5像素获取失败，因此没有独立从图中读出PAM柱高/倍数，仅依据正文和完整图注。热灭活与空载对照均扣除，释氨检测不能被描述成该反应每个产物均逐项质谱验证；文中GC-MS是植物内源PAA鉴定，实验端点须分开。
- [Fusarium 2016研究](https://doi.org/10.1371/journal.pone.0147486)：完整原文已读AMD1相关段落。其EC3.5.1.4在当前条目是ECO:0000305推断；原生BOA诱导可有直接表达证据，但不使PAM活性成为实验确认。
- [ScAMD2原序列报告](https://doi.org/10.1093/nar/18.23.7180)：1990报告与当前条目没有为本次提供具体底物阳性；“无证据”不是“无活性”。
- [FAAH 2002结构研究](https://doi.org/10.1126/science.1076535)：本次读取原摘要和已保存PDB/UniProt记录，未取全文。原始 `1MT5_entry.json` 含引用元数据，没有摘要；因此按指定PMID12459591从Europe PMC取得原摘要，保存为 `audit_FAAH_pubmed.json`。该来源支持2.8 Å的抑制剂复合物及膜酶结构适应主题，不能描述为PAM共晶。P97612的跨膜9–29是ECO:0000255预测，膜相关定位另有该物种实验支持，证据类型分开。

## 计算输出与最终候选审计

本节计数与21条来源声明分开；以下为**独立静态核实已有交付物**，不是重新执行序列搜索或结构配对。

| ID | 本次独立检查覆盖 | 结果 |
|---|---|---|
| C1 | 全53条本地BLAST命中/HSP：XML原始字段、大小写归一后的原序列区间、identity、union coverage、逐残基映射与输入/XML SHA | passed；SEG产生的小写字母按氨基酸身份归一，未改原XML；只能作为20参考有限面板，不能声称Swiss-Prot全库最佳 |
| C2 | 9条比较链：原PDB与提取链的全部Cα残基、坐标、confidence、SHA；三目标与输入序列一致；六参考AF/UniProt或实验构建体身份 | passed；1MT5为37–573，1M22为38–540，4YJ6为495aa原生加His标签且MSE按M映射；未把实验B因子当pLDDT |
| C3 | 全18对：原始输出、每对对应残基和坐标、每个矩阵距离、coverage、打印TM指标、原始输出/矩阵SHA | passed；未重跑US-align；只比较完整可用单链Cα，未做配体对接、口袋能量或膜拓扑预测；抑制剂邻域静态映射见C5 |
| C4 | 全18对：从奇异值残差独立重建最小二乘RMSD，从已保存变换重建距离及参考归一化TM | passed；两种拟合用途不同。打印RMSD与最小二乘复算一致；保存矩阵属于TM拟合，其RMSD不强行等于打印值。TM_query与TM_reference均保留；独立重建的是参考TM |
| C5 | 全45催化位点：对应UniProt记录的位置/氨基酸/证据码、已有结构映射、距离与目标pLDDT；全23个1MT5/MAY重原子4 Å邻域及三目标69条记录 | passed；从原1MT5 24个MAY重原子独立重建完整23位点，核对全部距离和缺失匹配。E14619/E41276/F21875分别同残基8/8/10，映射21/22/23，分母均23；相似计数不转成特异性评分 |
| C6 | 已保存Swiss-Prot请求和最后响应的RID、输入SHA和WAITING状态 | passed；RID A90RZ327016，未取得全库HSP或数据库版本；仅接受“已提交但未完成”，没有把本地结果称全库最佳 |
| C7 | 最终报告的序列/结构表、催化位置和pLDDT、MdlY结构指标、FAAH五个示例替换、预测版本和置信摘要 | passed；逐项对照已核实结果。报告使用定性等级而非校准概率；检查了全部最终段落的来源层级与限制 |
| C8 | 最终三个功能候选及R2251/R2087结论 | passed作为候选陈述：E14619为AmdS样AS家族候选；E41276为乙酰胺利用相关/预测acetamidase且保留W29未直接实测；F21875为AS家族候选并保留低置信脂肪酰胺竞争假说。三者PAM、区室及OR/AND仍未确认，没有接受任何模型修改 |

八项均已审且支持。结构静态核验逐对记录、代码及输入完整SHA见 `sources/audit_structure_verification.json` / `sources/audit_structure_checks.py`，功能位点核验范围见 `sources/audit_features_verification.json`。结构结果文件SHA为 `3d617afe4af895bae30799be7a17c1042db0ed927f441f899d718f81edab0eef`；九链清单SHA为 `9835a2f06a21c0ba4169e0721639f91394afc9f5bda1ae583fb31248a09c9fe9`。当前代码/US-align执行文件SHA与结构记录一致。

## 最终审计锁定

2026-09-11 已读完并核对最终 `REPORT.md`，SHA256：`50ced875cdd9e1d3019b6cb04f14f8bd15101f0c6e20c34883c001f2d08b6957`。功能位点结果SHA256：`acc6dcc14a786f32a99725d31003ab73477958a1ea41d1e35346282312e181f5`。其他已审输入指纹见 `sources/audit_final_inputs.json`。若报告或相关分析输入改变，这一锁定不自动覆盖新版本。

本审计登记 **29项：29已审、24在所写证据层级受支持、5未决、0反驳、0未查**。分为21项来源/推断边界（16支持、5未决）和8项输出/最终陈述（8支持）。这不是29项生物学实验验证，亦非本项目全部可能声明的完备清单。最终正文没有把那5项更强但无证据的声明当作已确认结论。

审计建议已在锁定版落实：FAAH措辞明确为核心结构相似；删除本次七参考范围外且没有引用的乙酰苯胺阳性断言；为AmdS功能补充人工审阅条目来源。乙酰胺不等于PAM、YB-392残余生长、W29跨菌株桥接、1MT5膜锚缺失、所有口袋位点/缺失映射及排队中的Swiss-Prot结果均已保留。

另只读检查 `completion_verification.json` 列出的12个保护文件，当前字节SHA与其中记录的expected/actual一致；该检查不重建历史dirty环境，也不等同于完整工作区未变更证明。审计代理仅写本审计和获准的 `sources/audit_*`，没有执行新的BLAST、结构比对、AlphaFold、配体对接、代谢求解、模型/科学数据变更、Git或集群作业。
