# 原 GPR 两分支的 AlphaFold 多链比较：独立设计复核

核查日期：2026-09-18。范围为用户给定的前 12 种蛋白（A）和后 2 种蛋白（B）能否作为多链结构比较输入，以及预测结果能支持哪些判断。本文件是提交前的科学设计复核，不是新预测结果。没有提交作业、改变模型、运行生长求解或重新执行 BLAST/结构叠合。

**两组可分别做“一种蛋白一条链”的探索性 AlphaFold 预测，不能标成两台完整 V-ATPase，也不能把蛋白序列直接首尾融合成一条链。** 两组的组件类别和数量本来就不同；预测主要检查局部折叠、可置信的链间排列及与已知组件的对应。结构结果不直接裁定 AND/OR、原生定位或生长必需性。

## 输入的组件含义

下表各 Yarrowia 位点均无本轮核实的原生正式基因名。字母是**亚基类别候选**，不是正式符号；证据等级均为 **uncharacterized／原生功能未实验确认**。自动数据库注释、既有 BLAST 和 AlphaFold 预测的支持分别注明，不将它们升级为实验确认。

| 分组 | 系统 ID | 蛋白功能候选；证据 | 固定蛋白／长度 |
|---|---|---|---|
| A | YALI1A11258g | V1 C，结构连接；自动家族注释 | A0A1H6Q1E3／383 aa |
| A | YALI1E14125g | V1 G，外周柄；自动家族注释 | A0A1D8NI06／115 aa |
| A | YALI1D00581g | V1 **D**，中央转轴；既有 BLAST／AlphaFold 比较支持，与 V0 小写 d 区分 | A0A1D8NCM2／257 aa |
| A | YALI1E10492g | V1 H，调节与稳定；自动家族注释 | A0A1H6PTA8／440 aa |
| A | YALI1A09766g | V1 A，ATP 水解催化；自动家族注释 | A0A1D8N4B6／577 aa |
| A | YALI1F20965g | V1 E，外周柄；提交注释 | A0A1D8NNL2／227 aa |
| A | YALI0E16192g | V1 F，中央转轴；CLIB122 序列，既有 BLAST／AlphaFold 支持；W29 身份冲突保留 | Q6C5Q2／122 aa |
| A | YALI1F31854g | V0 c，转子环蛋白脂质；既有 BLAST 支持 | A0A1D8NPS9／164 aa |
| A | YALI1F21690g | V0 **d**，转子连接；InterPro IPR016727、PF01992 等家族注释 | A0A1H6PZG0／346 aa |
| A | YALI1E37063g | V0 c′，转子环蛋白脂质；既有 BLAST 支持，参考 N 端覆盖不全 | A0A1D8NKR8／138 aa |
| A | YALI1E32332g | V1 B，核苷酸结合及催化头结构；自动家族注释 | A0A1D8NK63／512 aa |
| A | YALI1E12482g | V0 a，质子通道与 V1–V0 连接；自动家族注释及既有 BLAST 支持，区室未定 | A0A1D8NHU2／825 aa |
| B | YALI1F38820g | 偏 Vph1-like 的 V0 a 候选，质子通道与连接；既有序列／AlphaFold 支持，**仍需实验确认** | A0A1H6PMT1／804 aa |
| B | YALI1F13017g | V0 c″，转子环蛋白脂质；既有 BLAST／AlphaFold 支持 | A0A1D8NMN8／196 aa |

来源为本地固定 [14 条数据库记录](../r794_r795_gpr_review_20260909/gene_sequences.json)、[完成的序列／结构比较](../r794_r795_gpr_review_20260909/continuation_20260910/REPORT.md)及 [D/F 后续身份复核](../r794_r795_gpr_review_20260909/target_three_review_20260910/REPORT.md)。本轮直接读取数据库条目的功能、证据代码、结构域与长度；历史计算数值仅复用，没有宣称本次重算。

A 为 **12 链、4106 aa**；B 为 **2 链、1000 aa**，均假设每种蛋白 1 copy。根据当前映射，A 中未识别出 Vma12/21/22 一类独立装配伴侣；a 的“参与装配/连接”不等于它只是装配因子。A 含 V1 A/B/C/D/E/F/G/H 和 V0 a/d/c/c′；B 是 a/c″局部组件候选。B 未包含已知的 ATP 水解催化头。

**菌株限制：** A 的 F 链来自 CLIB122，其他序列按既有 W29 记录固定，因此 A 是混合来源的假设组。旧 F 位点所映射的 W29 位点 YALI1E19360g（无已核实正式名；旧 F 候选的截短对应，数据库注释）存在伪基因/移码记录。不能把 Q6C5Q2 自动改名成 W29 完整蛋白，也不能把其良好预测当作解决位点冲突。该冲突限制原生复合体声明，但不妨碍透明标注的旧 GPR 分支假设计算。

## 拷贝数与缺失组件

直接复核 [Vasanthakumar 等，2019，PNAS，DOI 10.1073/pnas.1814818116](https://pmc.ncbi.nlm.nih.gov/articles/PMC6462096/) 原文 Introduction、Overall Structures 和 Fig. 3：酿酒酵母参考 V1 有 **A3B3CDFH，加 3 套 EG 外周柄**；V0 有 **a c8 c′ c″ d e**，该研究还解析/鉴定 f 和 Voa1p。它们是酿酒酵母实验参考，不能自动作为 Yarrowia 精确拷贝数实证。

相对于这个参考构架，A 缺 c″，一拷贝输入也缺 A/B/E/G 的重复链及完整 c 环；A/B 都没有指定 e/f/Voa1 对应序列。B 还缺 V1 催化头、中央转轴、外周柄、d 和其余转子环成员。因此两组均不能仅按输入列表称为成熟完整泵。不能为了凑齐参考结构，静默增加用户未给定蛋白、复制尚未确认的亚基或替换菌株序列。

原文研究的是完整复合体的 ATPase 活性以及分离 V0 的结构，并观察到相似 V0 构架下活性/调节差异。由此可知，共有形状和局部界面本身不能证明等同的 ATP 驱动质子泵功能。脂质环境、完整组件和耦联状态也会影响功能。

## 预测与比较验收

1. **输入验收。** 逐链保留 accession、序列版本、菌株、长度和 SHA；多条 FASTA 记录表示多链，不添加人工 linker。确认输出链数和逐链序列等于提交输入。实际 AlphaFold/Multimer 版本、数据库快照、模板截止日、模型/seed/recycle 数和资源上限由执行记录提供，未知字段保持未知。旧 AFDB 单体只用作参照，不能称为本次复合体预测。
2. **输出完整性。** AF-Multimer 保留每个实际运行模型的坐标、pLDDT、完整 PAE、pTM、ipTM 及排序来源；AF3 保留 mmCIF、完整和汇总 confidence JSON、每个 seed/sample，另看 chain-pair ipTM、chain-pair PAE 和 clash 标志。缺任一指标则报告缺失，不用另一指标代填。
3. **置信度解释。** pLDDT 评估局部构象，PAE 评估相对位置不确定性，pTM/ipTM 是预测结构/界面置信度，均不是相互作用或功能相同的概率。检查每条链及链对，不能让一条长且高置信的 a 链掩盖其余链的低置信排列。官方 AF3 文档的 ipTM >0.8、<0.6、0.6–0.8 分别是高置信、失败提示、灰区的经验解释，**不是本项目的功能验收阈值**。低界面置信度只能说该预测未解决装配，不能证明真实不结合。
4. **结构比较。** 先做有生物学依据的组件对应：两组 a 候选逐链比较；各组件与实验亚基比较；若局部界面有可比对应，再检验其相对排列与接触。每次报告对齐残基数、两侧全长覆盖、双向长度归一 TM-score、RMSD、残基对应及缺坐标区域。A 不含 c″，所以不能把 B 的 c″强行映射为 A 的 c/c′后宣称同一亚基。全复合体叠合最多作为形状描述，不能用单个最优 RMSD 把 4106 aa 和 1000 aa 判为同功能；小片段低 RMSD 和大部分残基未对齐可同时发生。
5. **结论边界。** 若获得稳定、相互一致的高置信局部接口，可写“基于 AlphaFold 预测的局部装配候选”。若不同模型/seed 差异大或 PAE 高，写“装配未解决”。即使 B 的 a/c″界面很可信，也没有因此获得 ATP 水解催化头。全功能等价仍须 ATP 水解、耦联质子转运、原生定位及互补实验；不能据此次结构任务单独批准 GPR AND/OR 或 essential 标签。

以上为本项目的分析与验收方法，不归称为原作者的 GPR 修复流程。参考结构可能已进入 AlphaFold 模板/训练资料，因此与其相似不自动构成独立验证。

官方方法来源：[AlphaFold-Multimer 输入与输出文档](https://github.com/google-deepmind/alphafold/blob/main/README.md)、[AlphaFold 置信度实现](https://github.com/google-deepmind/alphafold/blob/main/alphafold/common/confidence.py)、[AlphaFold 3 输出文档](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md)。2026-09-18 通过网页读取；本地下载尝试因网络解析受限未完成，不声称已留存这些网页原始快照。文档 main 版本并非 HPCC 实际执行版本。

## 本次独立复核覆盖

分母仅为下列 7 项设计主张；不合并旧 23/25/33 项审计，也不把后续预测提前计入。

| 主张 | 所查依据 | 裁决 |
|---|---|---|
| A/B 的候选亚基类别及一拷贝长度 | 固定数据库条目、已有完成的序列/结构证据、长度求和 | supported，限候选身份 |
| A 未识别出独立装配伴侣，a 不是仅因参与装配就归作伴侣 | 固定家族注释；2019 结构中的 a 位置与作用 | supported，限本次所给成员 |
| 酿酒酵母参考拷贝数和 V0 组成 | 2019 原文 Introduction、Fig. 3、Overall Structures | supported，跨物种外推另列 |
| 一拷贝 A/B 不是上述完整参考构架 | 组件列表与原文参考逐项比较 | supported，结构比较推论 |
| PAE/pLDDT/pTM/ipTM 不能当成功能概率 | 官方输出说明与置信度计算实现 | supported，方法解释 |
| 全局 RMSD 不足以对两种不等组成证明同功能 | 对齐覆盖定义、组件差异、原文相似 V0 与功能差异 | supported，方法与机制推论 |
| 两组在 W29 原生条件下有相同完整泵功能 | 缺目标完整复合体/定位/耦联实验，F 位点冲突未解 | unresolved／unverified |

**7 total | 7 audited | 6 supported | 1 unresolved | 0 contradicted | 0 unchecked**。这表示上述设计主张已逐项检查，不表示 14 个原生基因功能已实验验证。

## 独立输入审核（2026-09-18 22:02 UTC）

在序列整理代理确认最终版之后，直接重新读取 [sequence_manifest.json](sequence_manifest.json)、两份 FASTA、14 份冻结 UniProt JSON、它们所指原始快照，以及 14 条版本化 GenPept 的 `VERSION`、`ORIGIN`、`/strain`、`/locus_tag` 字段。没有执行准备脚本、采信其中的 `verification=true` 代替检查，也没有修改输入。

独立核查结果：**14/14 输入链通过，0 条身份/序列/拷贝数差异。** A 的 12 位点顺序及 B 的 2 位点顺序与用户本次明文要求完全相同；没有交换、删除、增添或拼接蛋白。各组链 ID 从 A 开始连续，每种蛋白一条链；总长度分别为 4106/1000 aa。逐条 accession、locus、序列版本、条目版本、序列更新日期、长度、完整序列 SHA、原始与冻结快照 SHA，均与直接读取的记录一致；历史序列清单的长度/SHA 也一致。14 条仅含标准 20 种氨基酸。

GenPept 的原始 `ORIGIN` 序列全部与对应 FASTA/manifest/UniProt 精确相同，13 条原始 `/strain` 是 `CLIB89(W29)`，旧 F 链为 `CLIB122`，与 manifest 的菌株声明及 A 混合来源标记相符。另直接打开 W29 冲突位点 GenBank，确认其保留 `YALI1_E19360g`、`nonfunctional due to frameshift` 和 `/pseudo` 注释。**输入身份检查通过没有消除这个生物学冲突。** 所有 manifest 引用来源的文件 SHA 也逐一复算吻合。

同时独立读取根任务生成的 [branch_A.json](branch_A.json)／[branch_B.json](branch_B.json)：两份均为 `alphafold3` dialect、输入格式版本 1、seed `[0]`；14 条 `protein.id` 与 `protein.sequence` 按顺序和最终 manifest/FASTA 逐条完全一致。每条只指定 `id`、`sequence`，没有隐含重复链、配体或额外蛋白；两个名称均注明 `one_copy_hypothesis`。这是本地输入对应检查，HPCC 的实际 AF3 版本、资源、解析和作业状态以根任务执行记录为准，本审核未自行执行远程解析或提交。

已核对象的完整 SHA-256：

| 对象 | SHA-256 |
|---|---|
| sequence_manifest.json | `e0b82204f30d60c64b19ba2724134c371114e3780adfbff2c316f2445a760118` |
| inputs/branch_A_12chains.fasta | `c8a5274ae01923ce06d2e6e14d59c94f277fbe4e04c2926fccd02dbd2ae0aec8` |
| inputs/branch_B_2chains.fasta | `15dd7c0d2a97bf5c0e5c42726ad19aa685ce4fdd3b14a21368ba7b16fc6e101b` |
| branch_A.json | `98d097943026e982001bb23ca40e446d8c910bc76ec40ee6215eca5e35bb77b2` |
| branch_B.json | `fb37d0518bd673fbb092570ba17771a5f48ed371b8e637819dd47b77f9497549` |

本节分母为 **14 条序列输入，14 条独立直接核查通过**；不与上面的 7 项设计主张混成“21 项生物学已验证”。保留的科学未决项仍是原生功能、真实拷贝数/装配、W29 F 身份及完整泵功能等价。
