# GUT2—CoQ 电子耦联审查

2026-09-18。本轮只读审查当前 `model_metadata_trna_r539_alphafold_labeled.xml`，没有修改模型、构建管线、GPR、培养条件或实验标签，没有运行优化、screen或新结构预测。输入身份、源码身份和已有dirty状态见scope.json；原始提取、可运行检查、序列和独立来源审核保存在本目录。

## 判断

**R347的GUT2功能赋值有依据；R347→共享FAD/FADH₂池→R1977的基因依赖拆分不能按已证实生物学机制接受。** 两反应相加的净化学可以配平，但这不证明游离FADH₂是原生可交换电子载体，也不证明R1977的基因是GUT2所必需的独立下游酶。

- **YALI1B18499g / YALI0B13970g — GUT2**：甘油-3-磷酸脱氢酶，将G3P氧化为DHAP并向醌供电子。原生缺失的甘油利用/脂质表型支持功能；膜结合线粒体定位在原生文献中是路径注释，两篇研究没有直接定位实验。[Beopoulos 2008，Table 1及GUT2 deletion结果](https://pmc.ncbi.nlm.nih.gov/articles/PMC2607157/)、[Lubuta 2019，glycerol catabolism结果](https://academic.oup.com/g3journal/article/9/12/4059/6028124)。
- **YALI1F31153g — 正式蛋白名称未核实**：模型把它赋给R1977和R1975；原生CRISPR筛选只有脂质染色表型，不能确定ETF-QO催化功能。旧位点的酰基-CoA脱氢酶线索与当前W29蛋白序列未闭合，见下文。**不能称为已确认ETFDH，也不能据此删除该基因。**
- **YALI1A16988g — FLX1样候选，原生名称/底物未实验确认**：版本化序列支持线粒体载体家族与FAD载体同源候选，不能证明当前R847所写的FADH₂单向输出。

## 当前存储反应与精确守恒

以下`mi`实际区室名是“mitochondria”，不是已经证明的基质侧；模型另有`C_mm=mitochondrial membrane`。R347、R1977均把全部物种放在C_mi。

```
R347  G3P[mi] + FAD[mi] → DHAP[mi] + FADH2[mi] + 2 H+[mi]
      GPR: YALI1B18499g
R1977 FADH2[mi] + Q9[mi] → FAD[mi] + Q9H2[mi]
      GPR: YALI1F31153g
相加  G3P[mi] + Q9[mi] → DHAP[mi] + Q9H2[mi] + 2 H+[mi]
```

两反应边界均0..1000。相加后元素和电荷残差均零。这里2H⁺来自存储酸碱形式：G3P为C3H9O6P/0，DHAP为C3H5O6P/−2；它不是跨膜泵质子，不应因教科书省略H⁺就机械删除。

本次按名称和化学式交叉检查全部四个FAD/FADH₂物种，得到11条直接邻接反应，没有同区室同/反号计量的精确重复。四条稳态行是：

```
FAD_mi:   v1125 − v347 − v672 − v1971 − v1976 + v1977 = 0
FADH2_mi: v347 + v672 + v1971 + v1976 + v2117 − v1977 − v847 = 0
FAD_cy:   v306 − v2082 − v2117 − v2142 = 0
FADH2_cy: v847 + v2082 + v2142 = 0
```

相加仅剩`v306+v1125=0`。两条合成反应下界均0，因此本SBML稳态强制它们均为0。再代入第一行：

**v1977 = v347 + v672 + v1971 + v1976。**

右侧四条反应均不可逆，所以关闭R1977会强制四者全部为0。R847/R2117并不能解除这个必要约束；不能仅凭连通图把它们称为逃逸旁路。四条合成/氧化还原守恒行不表达辅基装配或生长稀释需求，合成反应被迫零通量也不表示细胞不需要FAD。

静态真实gene KO核验：GUT2只直接关闭R347；YALI1F31153g同时直接关闭R1977和R1975（2-oxoadipate dehydrogenase）。因此这个基因的KO结果不能归因于R1977单反应。此处证明的是稳态必要约束和GPR传导，不是WT可行通量、生长率或原生必需性；本轮新增LP=0。

## 十一条FAD邻接反应

| 反应 | 当前区室/功能及GPR | 审查结果 |
|---|---|---|
| R347 | mi；G3P→DHAP；YALI1B18499g/GUT2，醌耦联G3P脱氢酶，原生遗传支持 | GUT2赋值合理；将其电子出口接入共享游离FADH₂池缺少机制依据 |
| R1977 | mi；FADH₂+Q9→FAD+Q9H₂；YALI1F31153g，正式名/催化功能未核实 | ETF-QO身份未确认；PROTEIN_CLASS=1.3.8.6实际指glutaryl-CoA dehydrogenase，与本反应不符 |
| R672 | mi；proline→P5C；YALI1B12966g，名称未原生核实，数据库推定proline dehydrogenase/EC1.5.5.2 | 名称仍写“proline oxidase (NAD)”，计量却用FAD；EC所描述的是经酶FAD向quinone传电子，同样被R1977额外耦联 |
| R1971 | mi；isovaleryl-CoA→methylcrotonyl-CoA；YALI0D12573g，名称/身份未核实，模型赋为相应脱氢酶 | 旧ID保留；本轮未闭合其序列和底物特异性，不替换为另一个相似反应基因 |
| R1976 | mi；glutaryl-CoA→crotonyl-CoA+CO₂；YALI1E39693g AND YALI1E20192g AND YALI1D26360g | 三成员数据库分别指向2-oxoglutarate DH E1、E2、E3，和当前底物/电子机制冲突；不是已证实glutaryl-CoA DH复合体 |
| R847 | FADH₂ mi→cy；YALI1A16988g，FLX1样载体候选 | 同源FAD载体证据不能证明还原型FADH₂作为穿梭电子载体或当前单向输出 |
| R306 | cy；FMN+ATP→FAD+PPi；YALI1D33889g，原生名未核实，数据库推定FAD synthase | 催化注释层面一致；当前四池总守恒强制其稳态通量为0 |
| R1125 | mi；同类FAD合成，无GPR | 与胞质R306区室不同，不能直接当冗余同工反应；同样被总守恒强制为0 |
| R2082 | cy；isovaleryl-CoA氧化；YALI1E15555g，原生名未核实，数据库推定acyl-CoA dehydrogenase/oxidase | 与R1971不同区室、不同GPR；底物特异性和电子载体未闭合，不自动合并 |
| R2117 | cy酰基-CoA反应，FAD在cy、FADH₂却在mi；YALI1D19252g，原生名未核实，数据库推定线粒体short/branched-chain acyl-CoA DH | 区室拆分明确不一致；反应可逆，但不能仅凭边界推断逆向在完整网络可行 |
| R2142 | cy；propanoyl-CoA→propenoyl-CoA，无GPR，可逆 | 名称写ETF，物种却是FAD/FADH₂；没有表达ETF蛋白及其依赖 |

R1976的三成员身份限数据库自动注释：YALI1E39693g—推定2-oxoglutarate dehydrogenase E1/EC1.2.4.2；YALI1E20192g—推定dihydrolipoyl succinyltransferase E2/EC2.3.1.61；YALI1D26360g—推定dihydrolipoyl dehydrogenase E3/EC1.8.1.4。三者原生正式基因名与本反应催化组成均未在本轮实验确认。

9/11反应具备完整存储化学式且残差零。R2117产物m1972[C_cy]、R2142产物m1991[C_cy]缺化学式，电荷存为0，COBRA因此报告非零残差；**这是输入不全且电荷可疑，不能把欠缺的式子当作已证实丢失整段CoA，也不能宣布配平通过。**

## 蛋白身份与原生机制的证据边界

GUT2 EC1.1.5.3的净受体是quinone；FAD是酶的辅因子。[IUBMB定义](https://iubmb.qmul.ac.uk/enzyme/EC1/1/5/3.html)和原始GlpD结构研究提供类级支持，后者为大肠杆菌，不能充当Yarrowia原生定位实验。[Yeh 2008](https://doi.org/10.1073/pnas.0712331105)。真正ETF-QO/EC1.5.5.1接收**还原ETF蛋白**，不等于接收共享游离FADH₂。[IUBMB ETF-QO](https://iubmb.qmul.ac.uk/enzyme/EC1/5/5/1.html)。同理，[proline DH EC1.5.5.2](https://iubmb.qmul.ac.uk/enzyme/EC1/5/5/2.html)与[glutaryl-CoA DH EC1.3.8.6](https://iubmb.qmul.ac.uk/enzyme/EC1/3/8/6.html)的电子出口不能混为同一种独立酶。

R1977位点需要先解决跨版本身份：本地映射把YALI1F31153g对应至旧YALI0F23749g。旧CLIB122蛋白Q6C0L1为417aa，在限定面板中与人GCDH/Q92947（戊二酰-CoA脱氢酶）的identity为57.314%，E=5.03e−170；这支持旧蛋白的GCDH样功能候选。当前W29版本AOW07619.1为426aa，与缓存A0A1D8NPR7和其既有AlphaFold序列完全相同，却未命中旧蛋白或面板内ETF-QO/GCDH/IVD参考（E≤1e−5）。GenBank的“probable glutaryl-CoA dehydrogenase”只是序列相似性备注，不是酶活验证。**不能把跨版本位点映射当成蛋白序列等同，也不能将无小面板命中当作排除一切功能。**

该426aa序列的既有**AlphaFold预测**置信度很低：API全局pLDDT=32.44，按下载PDB的CA重算32.4523，99.765%残基低于50，平均PAE=28.253Å，不能据其折叠认定或否定ETFDH/GCDH。未新增预测，也未进行可信度不足的结构对齐。

R847候选AOW00753.1/319aa与缓存A0A1D8N536及既有AlphaFold序列完全一致。相对酿酒酵母FLX1（YIL134W/P40464，FAD载体，有该物种实验支持），identity39.13%、query覆盖296/319；初始参考面板E=3.39e−56，扩展面板E=4.01e−56。MIR1（YJR077C/P23641，磷酸载体）和PET9/AAC2（YBL030C/P18239，ADP/ATP载体）已审阅参考对照较弱。既有**AlphaFold预测**的PDB平均pLDDT82.197、平均PAE9.209Å，支持继续研究载体候选，不证明FADH₂特异性或方向。1996年FLX1实验测的是FAD稳态/跨膜囊泡通量；本轮该原文仅取得原始摘要级内容，不能升级为本菌运输证据。[Tzagoloff 1996](https://doi.org/10.1074/jbc.271.13.7392)。

上述两条UniProt accession现行接口均为Inactive/DELETED，理由“Not part of a reference proteome”；本轮明确使用保存的版本化缓存、GenBank及同序列既有预测，没有称其现行活跃审阅条目。具体序列SHA、数据库版本、BLAST参数、完整alignment和AF版本见identity_check.json及来源文件。

## 区室与构建来源

一般真核GUT2位于线粒体内膜外侧；当前C_mi物种无法区分内外膜面，因此只确认“线粒体相关”，不确认“基质催化”。R1142允许G3P在cy/mi交换，R1099允许DHAP在mi/cy交换，R1699连接mi/mm；三者可逆且无GPR。R349另外在mi用NADH还原DHAP，GPR为YALI1B04433g（模型NAD依赖G3P脱氢酶；本轮未复核其原生名称/定位）。这些连接提供模型供给，并不证明内膜真实转运或膜面定位，不能只合并两式后声称区室问题已解决。

原始输入data/iyali26.xml已包含R347、R1977、R847及当前各自GPR和边界。现有Q9身份步骤在scripts/gem_annotate/patches.py的连接反应集合中包含R1977，负责Q6/Q9物种转换；它没有建立新的GUT2—ETFDH证据。当前R347相对原始输入多2H⁺，与管线统一质子/水配平机制一致；本轮仅核对输入/输出和当前调用代码，没有重建完整历史运行。**共享FAD拆分和GPR不是本轮Q9 metadata校订新引入的。**

## 后续修订的最小科学边界（本轮未实施）

1. 以GUT2催化的“G3P+Q9→DHAP+Q9H₂”为R347净反应候选，按最终物种形式配平；解除对未确认R1977基因的催化依赖。保留GUT2现有身份，先明确区室近似，不能顺手迁移整条呼吸链。
2. 同时拆开proline直接quinone耦联与acyl-CoA—ETF—ETF-QO支路；R1977若拟表示ETF-QO，必须先确认原生蛋白及ETF依赖，不能把原有基因复制到所有FAD反应，也不能直接删除R1977而堵住四条通路。
3. R847的氧化态/方向、R2117的区室、R1975/R1976/R1977的底物与GPR错配分别保留为待决；原生蛋白身份未闭合前不重配。
4. 修改获准后，定向验证GUT2 KO直接关闭其净反应、独立关闭R1977不再仅因共享FAD池阻断GUT2，以及相关其他供电子支路仍有正确依赖；继续区分同时关闭R1975的YALI1F31153g基因KO。培养和判定阈值沿用既有有效设置。任何生长改善都不能代替酶学判断。

已完成本轮静态化学/逻辑检查、两次小参考面板BLAST和两份既有AF同序列核对。首次模型检查因名称附带化学式导致选择器未命中而停止，改为解析实际名称前缀后通过；首次身份检查发现项目环境没有Bio模块，改用标准氨基酸字典后通过，未安装或升级依赖。失败日志保留。独立审查覆盖9/9项关键主张，原生定位属于部分支持，R1977身份和FADH₂运输保持未确认；具体结论、来源访问程度与限制见source_audit.json。没有将这些未决事项标为通过。
