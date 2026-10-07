# 甲基柠檬酸两案的 AlphaFold 结构支持

本轮已实际取得两条固定当前序列的精确匹配 **AlphaFold预测**，并保存坐标、逐残基 pLDDT 和完整 PAE。用户后续选择官方网页版并完成登录后，相同458 aa输入已提交一次 **AlphaFold Server／AlphaFold 3** 单链预测，作业名为`XP_065950166_2_458aa_monomer_20260907`，网页显示处理中，尚无可用结果。页面对新任务的等待估计为最多约30分钟，不是本任务完成时间承诺。用户明确要求取消HPCC任务后，原作业28187510已于2026-09-08T00:08:52Z确认取消，实际运行0秒；仅保留网页预测。所有新增功能解释均为**基于AlphaFold预测的功能候选**，未改变模型或GPR。

| 系统ID、名称、简要功能与证据等级 | 固定输入 | 本轮逐项状态 | 模型角色（不等同功能确认） |
|---|---|---|---|
| YALI1C24124g — ICL1 — 异柠檬酸裂解酶候选；当前458 aa形式为 **model/GPR assignment only** | XP_065950166.2，458 aa | **已提交，待预测**：AlphaFold Server同序列单链作业处理中；原HPCC作业28187510已确认取消，未启动。尚无本目标坐标、pLDDT或PAE，540 aa P41555仅为比较对象 | R490／R552的共享OR成员 |
| YALI1F39620g — ICL2（2024靶向名）— 2-甲基异柠檬酸裂解酶候选；**curated annotation** | XP_506117.1，565 aa | **精确复用**：AF-Q6BZP5-F1，v6 | R490／R552的共享OR成员 |
| YALI1F03803g — PDH1/PHD1 — 2-甲基柠檬酸脱水酶候选；**curated annotation** | XP_504908.1，520 aa | **精确复用**：AF-Q6C354-F1，v6 | R95当前单基因赋值；未证明覆盖总异构化反应 |

三条固定输入均重新与本次从[NCBI公开RefSeq记录取得的FASTA](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=XP_065950166.2,XP_506117.1,XP_504908.1&rettype=fasta&retmode=text)逐残基核对，长度与序列完全一致。两份复用模型的AFDB序列、当前UniProt序列、PDB的SEQRES与全部Cα残基、mmCIF规范序列也一致。模型来源标注CLIB122；与固定CLIB89蛋白的序列相同允许结构复用，但不证明菌株、表达环境或历史实验构建相同。早期身份表的pending没有覆盖后续独立审阅；历史构建仍按原审阅保留缺口。

三个取得的AFDB文件均报告 **AlphaFold Monomer v2.0 pipeline、模型文件版本v6、创建日期2025-08-01**；本次获取日期为2026-09-07。文件版本v6不是AlphaFold算法“第六代”。这些文件是既有模型复用；另行提交的458 aa单体作业使用AlphaFold 2.3.0，尚未启动。本轮没有训练或多聚体计算。[裂解酶目标的官方记录](https://alphafold.ebi.ac.uk/api/prediction/Q6BZP5)、[脱水酶记录](https://alphafold.ebi.ac.uk/api/prediction/Q6C354)、[全长比较对象记录](https://alphafold.ebi.ac.uk/api/prediction/P41555)。

| AlphaFold预测对象 | 从逐残基文件重算的平均pLDDT | pLDDT >90 | pLDDT <70 | 完整PAE |
|---|---:|---:|---:|---|
| Q6BZP5，565 aa固定目标 | 92.92 | 495/565 | 31/565 | 565×565，已保存 |
| Q6C354，520 aa固定目标 | 92.52 | 465/520 | 39/520 | 520×520，已保存 |
| P41555，540 aa比较对象 | 96.02 | 505/540 | 14/540 | 540×540，已保存；不属于458 aa目标 |

PDB逐残基值与confidence JSON一致，mmCIF全局值也与其均值舍入一致；API摘要分别为92.94、92.50、96.00，存在0.02量级的小差异，原因未查明。本报告采用坐标／逐残基文件重算值，并保留API原值。pLDDT评估局部预测置信度，PAE评估残基或域之间相对位置的不确定性；它们不测量稳定性、酶活或底物特异性。[AlphaFold官方解释](https://www.alphafold.ebi.ac.uk/faq)。

![AlphaFold预测质量与全长比较对象截短边界](/Users/david/.codex/worktrees/6033/iyali26_gem/artifacts/alphafold_methylcitrate_support_20260907/prediction_quality.png)

图中为Cα轨迹、逐残基pLDDT和有方向的PAE矩阵。全长比较对象的橙色标记为1–82位；没有裁剪该模型来制造458 aa目标预测。PAE以Å计。坐标来自Google DeepMind／EMBL-EBI AlphaFold DB，CC BY 4.0；本轮重新绘图。

**裂解酶案 EGC-4b1e970207ef（R490／R552）。** Q6BZP5的Pfam ICL家族区间为42–565，域内平均pLDDT为95.04。其C238活性位点注释来自PROSITE规则，pLDDT为81.38；相邻G239为64.44，局部序列为`HGGKKCGHLAG`。因此高全局置信度不能升级为口袋细节或底物特异性已经确定。[UniProt位点和证据码](https://rest.uniprot.org/uniprotkb/Q6BZP5.json)、[Pfam区间](https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/Q6BZP5/)。

458 aa目标逐残基等于P41555的83–540位，前82位缺失。全长比较模型中该段平均pLDDT为97.24，全部≥70，且缺失区间与Pfam ICL家族10–536区间重叠。按“Cα距离≤8 Å、序列间隔>4、两端pLDDT≥70”的固定描述标准，缺失段与保留段有155对接触，涉及40个缺失段残基和59个保留段残基。该计数描述全长预测中的接触，不是截短体的失稳能量或失活判定。[全长序列与注释](https://rest.uniprot.org/uniprotkb/P41555.json)、[家族区间](https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/P41555/)。

| 全长比较对象的注释位置 | 458 aa目标中精确映射的位置 | 能说明什么 |
|---|---|---|
| C203，活性位点（相似性转移证据）；局部`PGTKKCGHMAG` | C121 | 对应残基保留；不能确认截短体的口袋形状或活性 |
| D165，Mg²⁺结合位点（相似性转移证据） | D83 | 对应残基保留；未预测金属结合或占据率 |
| 93–95、204–205、240、423–427、457，底物结合注释 | 11–13、122–123、158、341–345、375 | 位置按精确82位偏移映射；没有转移全长pLDDT |
| C端538–540 `SKL`（序列预测信号） | 456–458 `SKL` | 序列基序保留；不能证明本目标的区室 |

这组结果保留两种竞争解释：截短体可能仍具有兼容的催化骨架，也可能因失去参与折叠或相互作用的N端部分而改变功能。既不能仅因保留Cys宣布活性正常，也不能仅因缺失82位宣布失活。两个候选的ICL家族特征支持继续检查裂解酶功能，但不能据此删除共享OR连接、设定互斥底物或容量比例。实验比较对象[大肠杆菌PrpB的1MUM](https://www.rcsb.org/structure/1MUM)为1.90 Å晶体结构，呈无底物的开放活性环；另有[2005年原研究摘要](https://pubmed.ncbi.nlm.nih.gov/15723538/)报告配体结合时构象改变和溶液底物验证。这些跨物种结果说明静态单体折叠不足以决定特异性，不能替代当前酵母蛋白的直接实验。

**脱水酶案 EGC-aadd98455a18（R95）。** Q6C354具有PrpD家族的N端域53–307（平均pLDDT 97.24）及C端域324–498（96.98）。两个域之间PAE两个方向的平均值为4.07和4.74 Å，中位数4和5 Å，90百分位数6和7 Å。这支持该预测中两个域的相对排布，与PrpD家族候选解释相容。[UniProt家族与功能注释](https://rest.uniprot.org/uniprotkb/Q6C354.json)、[Pfam区间](https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/Q6C354/)。

当前UniProt记录没有给Q6C354逐残基活性位点注释，本轮没有建立经过验证的实验结构到目标残基映射，因此催化残基归属保持未知。有限实验参考[沙门菌PrpD的5MVI](https://www.rcsb.org/structure/5MVI)是3.05 Å晶体结构，主引文仍标为待发表；本轮只核对记录并取得坐标，没有计算目标与参考的RMSD或TM-score，也没有把参考残基直接转移为酵母活性位点。未安装新的比对平台。

PrpD家族相容性支持“2-甲基柠檬酸脱水至甲基顺乌头酸”的功能候选，不能证明同一蛋白还完成R95的第二步水合。两个结构域也不等于两种催化活性。[2022年大肠杆菌原研究的相关正文](https://www.nature.com/articles/s41467-022-33033-1)讨论PrpD的单步脱水及底物偏好可被改造；本轮只把它用于跨物种的功能边界说明，未转移其动力学、工程突变或补充图残基编号。既有两步酶学证据、第二步的其他酶解释及共享中间体身份问题仍需保留。

两份目标的N端低置信区分别为1–20和1–38，与数据库预测的1–19及1–37转运肽区间相邻或重叠。低pLDDT不验证切割位置、线粒体定位或成熟蛋白形式；本次复用的是完整固定序列。单体AlphaFold预测也未检验寡聚体、底物／金属结合、反应方向、胞内通量或历史实验靶向身份。

用户随后请求改用[官方AlphaFold Server](https://alphafoldserver.com/welcome)，并完成Google登录。本次只提交一个蛋白实体、一个拷贝，固定序列458 aa已与网页显示的全部残基逐一核对。实际提交预览采用固定种子1；种子0使当前网页提交按钮不可用。没有添加配体、修饰、自定义MSA或模板；服务默认设置待下载的请求文件进一步核实。网页作业名称为`XP_065950166_2_458aa_monomer_20260907`，记录时间显示2026-09-07 17:04（洛杉矶时间）。网页仅显示处理中，后台作业ID、实际启动和结束时间、坐标、pLDDT及PAE目前均未知。

AlphaFold 3网页算力由服务方管理，不能把原HPCC的1 GPU／8 CPU／64 GB／4小时限制当作网页资源记录。预测方法、数据库/模板、随机种子及pLDDT定义存在版本差异，取得结果后应分别记录；相同序列不意味着同一方法的重复测量。此次明确改用网页版产生一项新增网页提交，历史累计提交为HPCC一次、网页版一次，没有自动重试。详见[网页提交记录](/Users/david/.codex/worktrees/6033/iyali26_gem/artifacts/alphafold_methylcitrate_support_20260907/logs/alphafold_server_status.json)。

取消原HPCC排队作业的首次请求曾在进程启动前被自动审批拒绝，理由是缺少明确取消授权。用户随后明确回复“取消hpcc任务”，据此执行同一取消操作。实际调度记录确认28187510为CANCELLED，未曾启动，运行0秒，总排队36分钟；网页预测未受影响。取消记录保存在`logs/hpcc_job_28187510_cancelled.txt`。

以下为保留的HPCC提交与环境记录。用户确认已登录HPCC后，本轮实际连接到`skylark`，远程用户为`yhong075`。登录壳初始化后核实了AlphaFold 2.3.0、Singularity 3.9.3、GPU分区与`iwheeldonlab`账户。安装目录提供单体`monomer_ptm`模式和五套标准参数；本轮准备在一个作业内按既有流程串行运行，以取得PAE。Jackhmmer默认8 CPU已从安装源码核实。数据库实际Mgnify文件为2022_05版，与[UCR官方示例](https://hpcc.ucr.edu/manuals/hpc_cluster/selected_software/alphafold/)的旧文件名不同，提交脚本已采用实际路径。

已提交的唯一HPCC作业使用 **1 A100 GPU、8 CPU、64 GB、4小时、禁止重新排队与自动重试**，输入仅XP_065950166.2（458 aa），`monomer_ptm/reduced_dbs`、随机种子0。五套参数是标准单体流程的一次输入运行，不是五个预测作业或参数网格。固定FASTA和脚本已上传，远端SHA-256均与本地一致；实际调度记录确认上述资源、`Requeue=0`和`Restarts=0`。**job ID为28187510，提交数1；提交时间2026-09-07T23:32:52Z。** 2026-09-07T23:37:16Z检查时已排队264秒、运行0秒，状态PENDING、原因Priority；没有已分配GPU。启动前资料和环境核查含授权等待不足38分钟，未超过60分钟预算。两个已有精确模型的目标无需新预测。

此前上传曾被自动审批拒绝，理由为审批器未接受现有材料作为目标账户归属、公开序列和本次外发的可信授权。用户随后对列明源文件、远端账户/新目录和唯一作业资源的请求明确回复“确认”；依该授权完成同一上传与提交，没有换协议或绕过拒绝。[实际作业快照](/Users/david/.codex/worktrees/6033/iyali26_gem/artifacts/alphafold_methylcitrate_support_20260907/logs/hpcc_job_28187510_status.txt)保留了提交时间、输入校验、资源和排期。远端输出目录为`/bigdata/iwheeldonlab/yhong075/iyali26_alphafold_methylcitrate_20260907_01a07e16`。

额外的AFDB MD5精确查询两次均在进程启动前被自动审批拒绝。第二次之前，三条序列已与公开NCBI FASTA核实相同；审批仍认为具体外发授权未成立。没有换渠道绕过，因此本报告仅称“固定标识检索未找到458 aa精确模型”，不称“AFDB中不存在该模型”。

流程实际应用了 **govern-agentic-research**（中性问题、预算、支持/冲突和人类闸门）、**gene-identity-function**（版本、功能与GPR分开）、**academic-research-suite**及其deep-research的 **fact-check／source_verification**限定阶段（原数据库和原文定位）、**ponytail**（现有工具与标准解析、实际输入断言）及 **Spreadsheets**（只读身份TSV、保留来源和未知值）。没有处理PDF图表，未加载PDF技能；没有展开全套ARS或自派独立审阅。

本轮检验包含固定序列、源文件与暂定参考身份、实际HPCC单体环境、上传文件校验及唯一作业资源、PDB／mmCIF序列、逐残基质量一致性、PAE尺寸和有限性、82位精确映射、实际接触计数和图像目视检查。尚待458 aa网页作业完成并收取、验序及分析结果；未完成AFDB精确摘要查询、实验结构定量叠合、原生底物/定位实验和新增主张的独立来源审计。12条新增主张均为作者观察或有限推论，**独立审计覆盖0/12，全部unchecked**，供协调者转交独立审阅，不自行升级。

历史正例覆盖322/1612、交集内召回67/322和冻结实验标签保持原样；本轮没有代谢求解，不更新任何essentiality指标。成果根目录为 `/Users/david/.codex/worktrees/6033/iyali26_gem/artifacts/alphafold_methylcitrate_support_20260907/`。完整输入、代码、下载时间与来源身份保存在后台记录中。
