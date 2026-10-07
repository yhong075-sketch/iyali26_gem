# 酿酒酵母苯丙氨酸/PAA路线：独立来源审计

日期：2026-09-14。目标是分开回答特定单加氧活性、苯乙酸（PAA）形成、继续降解和外排证据。采用科研代理治理及基因身份技能；仅读来源并保存审计文件，未运行BLAST、结构预测、代谢求解、集群或模型/标签修改。来源检索限本任务约8篇/条及有界索引查询；最终候选与根代理计算结果另审。

## 已核实的基因身份

| 系统ID | 名称 | 蛋白功能及证据 | 本次应如何使用 |
|---|---|---|---|
| **YPL058C** | **PDR12**；UniProt **Q02785**，SGD **S000005979** | 质膜ABC弱酸外排转运蛋白；人工审阅注释，有原研究遗传学与荧光素外排实验 | 苯乙酸外排有遗传学支持；不能把它写成已测PAA转运参数 |
| **YML076C** | **WAR1**；UniProt **Q03631** | 弱酸响应转录因子，促进PDR12应激表达；人工审阅且引用实验 | 缺失表型支持同一调控系统；不是第二个PAA转运蛋白 |
| **YDR380W** | **ARO10** | 苯丙酮酸等2-氧代酸脱羧酶；2003原研究基因删除及酶活支持 | Ehrlich脱羧步骤，不能当作EC1.13.12.9的同义酶或同一步反应 |

PDR12条目版本209/序列版本1，WAR1条目版本172/序列版本1，两条注释更新日期均2026-09-02。本次不依这些条目的序列开展新计算。

## 特定单加氧活性与Ehrlich路线

当前UniProt **2026_03**（发行日期2026-09-02）中，`taxonomy_id:4932 AND ec:1.13.12.9`返回0条；涵盖酿酒酵母株系后，以phenylalanine oxidase或phenylalanine 2-monooxygenase作蛋白名称检索也返回0条。较窄的organism_id:4932检索同为0，未以它代替株系覆盖查询。

Europe PMC检索 `"Saccharomyces cerevisiae" AND ("1.13.12.9" OR "phenylalanine 2-monooxygenase" OR "phenylalanine oxidase")` 返回95条索引记录；本次只筛查这些记录的题名/摘要是否直接包含该具体活性，不把索引命中当作酶学阳性，也没有声称逐篇审阅95篇全文。该题名/摘要检查未产生具体酶学候选。与本次原文核查合并，结论仅为 **在所列检索范围内未找到S. cerevisiae该特定活性的直接证据**；不是证明基因组或任何培养条件下绝无此活性。

S. cerevisiae的明确路线是：

**L-苯丙氨酸 → 苯丙酮酸 → 苯乙醛 → 苯乙酸**，对应转氨、脱羧、醛氧化；苯乙醛也可还原为2-苯乙醇。该路线不需要苯乙酰胺中间体。[Vuralhan等2003原研究](https://doi.org/10.1128/AEM.69.8.4534-4541.2003)在CEN.PK113-7D、葡萄糖限制且Phe作唯一氮源的好氧连续培养中，HPLC检出上清PAA和2-苯乙醇；相应厌氧条件未检出PAA。原文Results/Table2及Discussion说明供给的芳香碳骨架可由残余Phe与这些产物约定量解释。

该培养观察支持PAA作为分泌产物；它不是所有条件下不再利用PAA的证明，也不是完整开环、矿化或同化为中央代谢物的途径鉴定。[SGD的phenylalanine degradation记录](https://pathway.yeastgenome.org/YEAST/NEW-IMAGE?object=PWY3O-4115&type=PATHWAY)同样描述Ehrlich路线，不能因标题含“degradation”就写成PAA完整降解。本轮未找到足以确认S. cerevisiae原生完整PAA继续降解的直接来源；不把其他真菌的路线或在S. cerevisiae异源表达的酶算作其原生证据。

## PDR12：哪些是直接测量，哪些是推断

[Hazelwood等2006](https://doi.org/10.1111/j.1567-1364.2006.00094.x)，PMID16911515，已读完整方法、结果及讨论。研究将两个实验层次结合：WT好氧葡萄糖限制、氨基酸唯一氮源时的表达诱导；以及同背景缺失株在外加PAA下的生长敏感性。后者的摇瓶为葡萄糖/硫酸铵培养基，考察初始pH及PAA剂量，并非直接测量内源PAA的分泌速率。Table5/6为两个独立摇瓶的平均值±mean deviation。作者未完成PDR12回补，以WAR1缺失的相似表型作为旁证，并说明PDR12长度给操作带来困难。

这组证据支持“Pdr12参与PAA相关耐受/外排”，但**不提供PAA直接转运速率、PAA的Km/Vmax或ATP:PAA化学计量**。也不能从耐受表型单独指定唯一转运蛋白。

[Holyoak等1999](https://doi.org/10.1128/JB.181.15.4644-4652.1999)，PMID10419965，已读原文：细胞装载荧光素后，测加入葡萄糖引发的外排；PDR12缺失及未预诱导为对照。山梨酸和苯甲酸竞争抑制荧光素外排，钒酸盐/细胞ATP结果支持耗能机制。**直接测量的底物是荧光素，不是苯乙酸**；本次不能将该实验的动力学或能源现象移作PAA参数。

## 原子声明登记

`source claims 16 | audited 16 | supported 11 | unresolved 5 | contradicted 0 | unchecked 0`

支持均限表内注明的证据层级；未决声明不能用于确证功能或修改模型。

| ID | 精确声明 | 来源及定位 | 判定 |
|---|---|---|---|
| S1 | Sc PDR12的系统ID是YPL058C，编码质膜ABC弱酸转运蛋白 | Q02785 genes/FUNCTION/SUBCELLULAR LOCATION；SGD-Wiki身份头部 | supported，curated annotation |
| S2 | YML076C/WAR1是PDR12弱酸响应的转录调节因子 | Q03631 genes/FUNCTION；2006方法和缺失表型 | supported；不是转运成员 |
| S3 | 本次精确EC与名称查询没有Sc匹配注释 | `sc_ec113129_taxonomy_query.json`及headers、名称查询JSON | supported的检索观察 |
| S4 | Sc不存在EC1.13.12.9或任何等效活性 | 上述有限检索及原文集合 | unverified；未找到不是不存在 |
| S5 | Sc经Ehrlich路线生成PAA；ARO10属于2-氧代酸脱羧步骤 | Vuralhan2003 Methods/Results、Figure1/Table2、原摘要明确YDR380W | supported，原研究；不是PAM途径 |
| S6 | 同一2003实验的厌氧条件未检出PAA | Vuralhan2003 Results/Table2及Discussion | supported，严格限该条件和检测能力 |
| S7 | 本轮已证明Sc具有完整PAA继续降解/开环同化路线 | 原文及SGD路线、定向检索 | unverified；分泌及Ehrlich“degradation”不支持此加强声明 |
| S8 | 氨基酸唯一氮源与WT PDR12转录上调相关 | Hazelwood2006 Microarray analysis、Results/Table4 | supported，表达/条件相关证据 |
| S9 | PDR12缺失增加外源PAA生长抑制且有pH依赖 | Hazelwood2006 Methods、Results/Table5/6 | supported，遗传学耐受终点 |
| S10 | 2006研究的实验设计没有完成PDR12回补 | Hazelwood2006对应Results段及菌株Table2 | supported；WAR1缺失同表型为旁证，不冒充回补 |
| S11 | 2006研究提供PAA直接外排速率或其转运动力学 | 全部Methods/Results | unverified；文中是诱导及生长表型 |
| S12 | 1999研究直接检测了PDR12依赖的荧光素外排 | Holyoak1999 Measurement of fluorescein efflux/Results | supported，直接外排实验，但底物为荧光素 |
| S13 | 该荧光素实验和酸竞争实验不能自动当作PAA直接测定 | Holyoak1999竞争方法与底物列表 | supported的证据边界 |
| S14 | 上述研究可指定ATP:PAA计量 | 两篇PDR12原研究 | unverified；不猜测1:1或其他比值 |
| S15 | Sc诱导/外排证据可确认W29同名或同源蛋白的底物功能 | 本次跨物种推断边界 | unverified；须另核W29身份与功能 |
| S16 | PAA形成、外排、继续降解是不同主张 | Vuralhan2003终产物数据；Hazelwood2006耐受终点；SGD路线范围 | supported；不能用其中之一替代其余 |

## 来源范围与获取限制

核心直接读取7篇/条：三篇原研究（2003 Vuralhan、2006 Hazelwood、1999 Holyoak），两个UniProt条目（PDR12、WAR1），SGD-Wiki PDR12身份头部，以及SGD phenylalanine degradation路线。搜索返回的其他论文仅作为筛选线索，没有借其题名/摘要宣称全文已审。主要来源全文/原摘要/JSON和完整SHA见 `audit_sources/provenance.json`。

SGD主站PDR12页面返回403；身份由UniProt及SGD-Wiki受保护的身份字段交叉确认，未依赖社区评论作关键功能证据。两篇PMC文章的fullTextXML接口返回404，随后成功读取并保存PMC HTML原文。Hazelwood出版社初始地址403，合法公开的article-abstract页面跳转至publisher article-minimal，获得完整正文。未从未目视的表格/图像独立估算额外效应量，未把文本表格解析当作图像验证。

## 已有比对、目录模型和Chalmers模型的独立静态核查

本节为根代理提供既有输出后的追加审查；审计者没有重新运行比对或代谢求解。独立读取FASTA、原始TSV、汇总JSON和XML，重核16对选定query/subject的比对字符串、身份数、覆盖区间并核对输入SHA。选择覆盖Sc的7个自身命中、W29的PDR12前三、PDR5第一及ARO8/9/10、ALD2/3各第一命中。该检查确认**所查输出可追溯且数字一致**，不声称本次独立重跑了BLAST，也不把16对当作整个数据库的功能验证。PAO查询序列Q5W9R9确实在输入中，Sc 6067条与W29历史缓存7894条蛋白的两份原始结果均没有它在E≤1e-5的命中；这仍是给定方法和蛋白集合内的阴性检索。

Sc **YPL058C/PDR12（Q02785）**对W29的首位是**YALI1_F23922g（A0A1H6PZY4）**，原生正式名称未核实；缓存提交名为ABC-2型转运蛋白结构域蛋白，属于注释和同源支持的转运候选。首位identity 50.83%、query coverage 95.04%、subject coverage 94.76%；后两位为**YALI1_B17155g（A0A1D8N7L3）**和**YALI1_C28310g（A0A1D8NBZ0）**，正式名称未核实，均为本次同源检索的ABC转运候选，identity分别39.46%、39.70%。这些数值不能确认PAA底物、ATP计量或实验定位。

对照Sc **YOR153W/PDR5（P33302）**是已注释的多药ABC外排蛋白；它的W29首位**YALI1_B03910g（A0A1H6PRY4）**与PDR12首位不同，正式名称未核实，仅为同源转运候选。不同首位有助于候选排序，但单向BLAST不是直系同源关系或底物专一性的证明。

目录模型中**m1932与m2026在胞质、同注释CHEBI:25978/MNXM1095022、charge均0且均缺少formula，却是不同ID**。m1932只由R2087/R2251产生；m2026由R2192产生、被R2191消耗成m2027，而m2027只有R2191生成。独立核对全部相关反应后，在普通稳态守恒且所记录非负方向下，有：

`m1931: −v2087−v2251=0`；`m1932: v2087+v2251=0`；`m2027: v2191=0`；`m2026: v2192−v2191=0`。

因此这四条反应在当前结构下均被这些行约束为0，不需要优化求解。原先“PAA没有去路”的表述只能用于m1932这个池；另一个PAA池**有写入的消耗步骤，但下游又断开**。同注释提示重复ID候选，不能代替化学身份整理或自动授权合并。**YALI1_F05415g（A0A1D8NLW9）**原生正式名称未核实，缓存“phenylacetate 2-hydroxylase”来自ProtNLM/ECO:0008006，属于预测，不能据此称W29已有实验验证的PAA羟化或完整降解。**YALI1_F06858g（A0A1D8NLZ9）**正式名称未核实，ALD同源候选在缓存中的线粒体醛脱氢酶5名称来自ARBA；模型R2192放在胞质并不证明该蛋白实验胞质定位。

用户追加要求的Chalmers两份官方模型已按固定提交及原始XML静态核对，XML的Git blob SHA与对应官方目录元数据一致：

| 模型 | 全部目标行与结论 |
|---|---|
| iYali v4.1.2，提交8e08e807cad59e265d751f0e2357f5af0ba30e67 | 胞质PAA M_s_1321只出现在R_y000185产物，边界0–1000；该行稳态强制此反应为0。未给此PAA池提供出口。 |
| Yeast-GEM v9.1.1，提交2d594ae1c4a2d550ccef120d96a58c7bbf586255 | 胞质PAM s_3867仅被r_4227消耗，即使其边界−1000–1000，单独守恒仍强制v_4227=0；PAA s_1321行再强制v_0185=0。它增加了PAM水解反应，但未解决所查池的稳态供给/去路。 |

Yeast-GEM的r_2001、r_2002、r_2003处理的是**苯乙醛**（phenylacetaldehyde），不是苯乙酸（PAA）。两模型基因条目中没有PDR12/YPL058C匹配；这只描述模型覆盖，不推断生物学不存在该功能。Yeast-GEM将r_4227赋给**YDR242W/AMD2（P22580；SGD S000002650）**，其UniProt版本174/序列2称“Probable amidase”，泛酰胺水解反应没有PAM特异实验支持；引用的记录位置是序列/基因组信息。[NCBI接收SGD的记录](https://www.ncbi.nlm.nih.gov/gene/851829)也将酰胺酶功能标为预测，GO证据为IEA。此GPR不能反向充当功能实验。

## 追加声明登记与覆盖

追加14项：`audited 14 | supported 11 | unresolved 3 | contradicted 0 | unchecked 0`。与上文16项合计**30项，已审30项；支持22项、未决8项**。未决项列的是不得加强为确认的主张，不是漏审。

| ID | 声明/问题 | 判定和证据 |
|---|---|---|
| C1 | PAO Q5W9R9在给定两个蛋白集合和参数下没有BLAST命中 | supported；两原始TSV、查询FASTA及参数 |
| C2 | 16对选定HSP的序列片段、identity和覆盖与记录一致 | supported；audit_sources/static_verification.json；不是重跑 |
| C3 | PDR12的W29前三位身份和报告数值可追溯 | supported；TSV/FASTA与历史缓存 |
| C4 | 首位W29候选已经证实运输PAA或可指定ATP计量 | unresolved；同源不是底物实证 |
| C5 | PDR5与PDR12在W29得到不同首位 | supported；有限对照排序，不是完整直系同源推断 |
| C6 | m1932/m2026同区室同外部注释但断成两个ID | supported；目录XML全部相关行；合并尚未接受 |
| C7 | 目录模型四个目标行足以在标准稳态下推出四条通量为0 | supported；完整行与边界；非生物学死亡或化学验收 |
| C8 | YALI1_F05415g的具体羟化酶名称来自ProtNLM预测 | supported；冻结UniProt缓存的证据代码 |
| C9 | 上述注释已证明W29 PAA羟化、完整降解或相关定位 | unresolved；模型赋值/自动名称不足 |
| C10 | Chalmers iYali的目标PAA池没有出口并强制所连产生反应为0 | supported；独立XML全行核对 |
| C11 | Chalmers Yeast-GEM的PAM水解和PAA生成仍被目标行阻断 | supported；独立XML全行及可逆边界核对 |
| C12 | Yeast-GEM r_2001–2003是苯乙醛输运；模型未覆盖PDR12条目 | supported；化合物ID/名称和全部geneProduct字段 |
| C13 | YDR242W名称为AMD2，已查来源只支持probable amidase层级 | supported；UniProt和NCBI/SGD证据字段 |
| C14 | Yeast-GEM的AMD2 GPR确认PAM底物和水解方向 | unresolved；GPR本身不是酶学或方向证据 |

本次追加来源为既有P22580条目、NCBI接收SGD的AMD2记录、两份Chalmers官方模型及根代理提供的静态计算材料；这是根代理明确转交的新增审计范围。全部SHA与定位见audit_sources/provenance.json、static_verification.json及chalmers_static_verification.json。

## 最终报告核对：2026-09-14

已读取最终REPORT.md全篇，并核对追加的两篇Yarrowia原文相关结果、方法、菌株表及图注：[Larroude等2021](https://doi.org/10.1111/1751-7915.13745)明确写基因从W29扩增、工程菌来自W29衍生Po1d，培养上清经HPLC，并有独立PAA标准校准；Ehrlich总产物不等于PAA单项，本文不把工程株表型外推为所有WT条件的分泌率。[Gu等2020](https://doi.org/10.1021/acssynbio.9b00468)写明ALD2/ALD3联合删除及比较株背景/表达构建差别，所以不从其PAA差值推断单基因独立效应；本次未采用其PAR编号或方向上的内部不一致构建新机制。没有从未目视的图像提取新数值，也没有把其他原文引用作为已阅读全文。

独立静态核对旧位点映射：ARO8/9/10和论文ALD3对应的缓存全序列分别与报告W29候选完全相同；论文ALD2旧记录491 aa与W29目标525 aa只存在旧2–491＝目标36–525的490 aa连续对应。它们是版本明确的序列桥接，不代表本轮重测作者实验所用每个克隆。当前CLIB122 Q6C1A3版本149/序列1的1508 aa与冻结W29 A0A1H6PZY4逐残基相同；目标inactive记录原因是“不再属于参考蛋白组”，并非基因删除证据。两个AlphaFold查询404从本轮获取日志核实，结论限未取得可用模型，未把404写成不存在结构或无功能。

报告新增AMD2段复用2026-09-11三目标结果：独立核对原local_blast_results.json完整SHA、三行全部所用字段，以及P22580 v174/序列2的序列SHA。29.98%、28.71%、27.72%一致率及全部覆盖百分比与旧输出一致；本次没有重新比对，也不是AMD2对整个W29蛋白组的搜索。三个目标的原生名称/文献对应与既有预测证据层级沿用已锁定的上一轮功能报告；Sc AMD2本身的PAM活性未定，因此该段保留AS家族候选，不升级PAM底物或OR关系。

另对报告中ALD2→W29 D候选及PDR5→W29 PDR12样候选补查两对原始HSP；连同前16对，共18对已有query/subject匹配经本轮静态核对。官方iYali与data/chalmers_iyali/iYali.xml确为字节相同；新Yeast-GEM与data/Yeast-GEM.xml确为不同，后者旧SHA为9fd2c572cace73c2ea835205617313554d1bf89f4ef077f49defbdfa219a4ad7。iYali的R_y002001–R_y002003同样核对为苯乙醛输运/交换。

| ID | 最终追加原子声明 | 判定和证据 |
|---|---|---|
| D1 | 2021来源能支持W29来源工程体系及胞外PAA检测，但不能把Ehrlich总量当作PAA单项 | supported；PMC8601196相关Results、菌株/培养/克隆/HPLC方法与图注 |
| D2 | 2020 ALD2/3证据限联合删除且有背景/表达差别，不能分摊单基因效应 | supported；PMC7308069相关Results与分析方法 |
| D3 | 四个旧位点全序列对应、ALD2仅核心对应，桥接范围与报告一致 | supported；两个冻结缓存和old_locus_mapping.json独立核对 |
| D4 | 最后两对已有HSP与49.49%及1137 bit score等报告值一致 | supported；原始TSV/序列；无重跑 |
| D5 | PDR候选当前同序列桥接、inactive原因及两次AF404被准确限定 | supported；当前JSON与获取日志；无功能/结构确认 |
| D6 | 旧/新模型字节关系与报告版本句一致 | supported；实际文件逐字节比较及SHA |
| D7 | 三个AMD2旧比对数值与输入版本确为可追溯复用 | supported；原结果完整SHA、参考序列版本与三行字段；不是新增或全库运行 |

**最终覆盖：37项，已审37项；支持29项、未决8项、矛盾0项、未检查0项。** 37项仅覆盖本文件登记的源声明和定向静态结果，不称全模型、完整文献或所有同源蛋白已审。未决项仍不构成科学确认、模型修改或候选接受。

最终已审REPORT.md SHA256：`110788d1ec784933f9b824361421fabf41c6e5e0303882cdeb44b4cace1bd86d`。最终核对记录见audit_sources/final_verification.json。报告已保留关键边界，未发现需要阻止交付的来源冲突；此审计本身没有新增BLAST、结构/代谢求解、集群作业、模型/GPR/标签修改或Git写入。
