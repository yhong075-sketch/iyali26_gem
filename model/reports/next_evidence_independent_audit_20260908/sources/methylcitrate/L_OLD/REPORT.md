# iYali26 甲基柠檬酸两案：首轮证据与限定提案

**两案均有实质进展，均未达到直接接纳GPR变更的条件。** R490/R552存在区分主要裂解活性的原生酶学依据，但ICL弱交叉活性、当前蛋白版本和定位限制互斥GPR拆分。R95得到两半反应的进一步直接支持，但第二步现代GPR及中间体连接未闭合。另发现R95与R552共用m217的名称和结构标识冲突，须在两案实施前共同审查。

本轮是已授权的只读证据研究和局部静态检查。没有运行COBRA、求解器、优化、历史run.sh、实验、模型构建、提交或推送，没有修改原项目文件，没有发送外部消息或新建代理。本报告不是病例接纳、候选激活或后续计算授权。

- [EGC-4b1e970207ef：R490/R552证据、主功能分配假设和反证](/Users/david/Documents/Codex/2026-09-06/iyali26-methylcitrate-research/outputs/EGC-4b1e970207ef.md)
- [EGC-aadd98455a18：R95两半反应、化学冲突和停止条件](/Users/david/Documents/Codex/2026-09-06/iyali26-methylcitrate-research/outputs/EGC-aadd98455a18.md)
- [27条原子主张](/Users/david/Documents/Codex/2026-09-06/iyali26-methylcitrate-research/outputs/claims.jsonl)、[15项来源记录](/Users/david/Documents/Codex/2026-09-06/iyali26-methylcitrate-research/outputs/source_records.json)、[独立审阅覆盖](/Users/david/Documents/Codex/2026-09-06/iyali26-methylcitrate-research/outputs/audit_coverage.json)

## 输入身份与核验范围

本轮暂定参考仍为 `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_integration/model.xml`，SHA256：

`bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee`

输入清单：`/Users/david/Documents/Codex/2026-09-06/yu/outputs/iyali26_research_tasks_20260906/input_manifest.json`，SHA256 `6d47f7cd420a9ecd99e84aca111ffb29198f276b0977cc053f71d625ffa1c387`。核对了实际使用的前轮快照、来源、主张、审阅，当前模型、本地映射、项目状态、baseline_manifest及项目AGENTS的登记哈希。加上清单本身、身份交接表及最终身份包的两项账本，共**13项输入哈希匹配**；这不是整个项目或身份包所有文件的重新审核。

具体模型检查只包括R490、R552、R95的GPR、区室、界限、参与物及元素/电荷；追加检查5个直接相关物种和7个邻接反应。三条目标反应按现存字段均配平。没有将前轮旧模型的残差或历史inactive/essentiality标签转为本轮结果。[静态快照](/Users/david/Documents/Codex/2026-09-06/iyali26-methylcitrate-research/outputs/static_snapshot.json)与[局部化学连接](/Users/david/Documents/Codex/2026-09-06/iyali26-methylcitrate-research/outputs/chemistry_connection_check.json)保留完整字段。

| 目标系统ID | 名称、简要蛋白功能、当前位点证据等级 | 身份交接v1的关键限制 |
|---|---|---|
| YALI1C24124g；旧YALI0C16885g | ICL1；候选异柠檬酸裂解酶，异柠檬酸→乙醛酸+琥珀酸；当前458 aa形式为`model/GPR assignment only` | 位点/出版映射和2024靶向已有联系；1993摘要555 aa、CAA51362.1的541 aa、P41555序列v3的540 aa与现行458 aa不相同。现行XP_065950166.2/AOW02989.1为截短形式；不能由此推断失活或定位 |
| YALI1F39620g；旧YALI0F31999g | ICL2（2024靶向名；UniProt无Name字段）；候选甲基异柠檬酸裂解酶，底物→丙酮酸+琥珀酸；`curated annotation` | Q6BZP5序列v1/XP_506117.1均565 aa且交接核对相同；生化原生酶到现代序列的直接鉴定链仍未闭合；不能把旧表中异种相似注释的non-functional移给本物种 |
| YALI1F03803g；旧YALI0F02497g | PDH1（UniProt）/PHD1（2013）；候选甲基柠檬酸脱水酶，生成甲基顺乌头酸；`curated annotation` | Q6C354序列v1/XP_504908.1均520 aa且交接核对相同；准确2013实验靶向序列仍未闭合；域名不证明完整R95 |

身份交接原表路径：`/Users/david/Documents/Codex/2026-09-06/yu/outputs/iyali26_research_tasks_20260906/identity_handoff/identity_table.tsv`；SHA256 `366f8ce0b74592a7eaf2b9cf3737788f01ac4ece943dda2a776b936cd5920e67`。核对4行，仅采用本任务3行，未扩展第四案。表由身份作者自检后冻结，**独立审阅pending**；本任务没有重做序列分析。全部版本化序列SHA与15项IDN来源元数据保存在[身份依赖](/Users/david/Documents/Codex/2026-09-06/iyali26-methylcitrate-research/outputs/identity_dependency.json)。

最终身份包位于 `/Users/david/Documents/Codex/2026-09-06/yu/outputs/iyali26_research_tasks_20260906/research/identity/`。本任务分别绑定其 `source_records.json` SHA `28fa9e3ff54849ff5bedee5d716aab6333ad243407f7199b842772921403b745` 与 `claims.jsonl` SHA `43ee44ccbd5e069e1408695cb93e3554f83569ac07c5597015afded8e5a0a944`：前者读取相关来源元数据，后者只校验字节，均未继承独立审阅状态。协调更正后的整包33项校验/22条主张属于其交接记录，本报告的自查计数仍为13项输入。

## 检索、实际访问与证据等级

日期：**2026-09-06 UTC**。检索从8篇指定文献的DOI/题名及原发布渠道出发，经合法机构库、官方摘要或作者公开上传版本补取。Rhea所引同一主题1995后续论文是唯一新增论文；精确PubChem请求仅解决现有m217/m76字段，不展开其他代谢路线。检索与访问的原始请求URL、失败/降级、页码、相关图表、材料绝对路径及SHA保存在来源账本。未保留的自由关键词输入不重构为虚假的完整搜索日志；本轮不声称系统综述、穷尽检索或“未找到即不存在”。

证据优先级：同物种原生酶分离/底物与产物实验、基因扰动、受控定位分别按其终点评估；版本化数据库用于身份/化学，摘要只支持摘要层陈述，继承记录只作线索。旧论文在本物种问题上可很直接，但没有现代序列身份；现代工程论文有准确构建，却未必直接回答酶特异性。期刊名称或数据库Reviewed标签不替代这一判断。

| 来源 | DOI / 实际访问 | 已核对位置及可用性 |
|---|---|---|
| S1976，主要裂解活性及弱交叉活性 | [10.1080/00021369.1976.10862298](https://doi.org/10.1080/00021369.1976.10862298)；出版社全文PDF | pp1864–1868；Fig2视觉核对，活性与讨论分开。可用于原生酶层，不能指认当前序列 |
| S1981，第一脱水半反应 | [10.1271/bbb1961.45.2823](https://doi.org/10.1271/bbb1961.45.2823)；出版社全文PDF | pp2824、2827–2828；Table I视觉核对。保留混合立体异构底物和检测界限 |
| S1982，循环酶定位 | [10.1111/j.1432-1033.1982.tb06713.x](https://doi.org/10.1111/j.1432-1033.1982.tb06713.x)；原全文未取得，官方摘要可得 | 只用摘要；分级纯度、回收率、具体株及重复数未核，不引用其图表 |
| S1993，历史ICL1克隆 | [10.1007/BF00284696](https://doi.org/10.1007/BF00284696)；出版社全文未取得，SONAR/官方摘要 | 克隆互补、删除、555 aa及潜在SKL限于摘要层；当前版本差异由身份依赖单列 |
| S1995，第二步原生酶分离 | [10.1271/bbb.59.1825](https://doi.org/10.1271/bbb.59.1825)；出版社扫描全文 | 全4页视觉阅读；Figs1、3–4、纯化表。原生活性证据直接，无现代GPR |
| S1995B，第二步酶性质比较 | [10.1271/bbb.59.2013](https://doi.org/10.1271/bbb.59.2013)；出版社扫描全文 | 全5页视觉阅读；Figs4–6及组成表。624残基是估计；同组相关证据，进化部分是作者假说 |
| S1996，ICL细胞器分级标记 | [10.1074/jbc.271.34.20300](https://doi.org/10.1074/jbc.271.34.20300)；出版社文字及机构全文PDF | pp20301、20304–20305、Fig8；图形已看，部分正文/图注字体渲染失败，结合提取文本核对；不数字化柱高 |
| S2013，第一酶旧位点与甘油研究 | [10.1016/j.jbiotec.2013.10.025](https://doi.org/10.1016/j.jbiotec.2013.10.025)；作者上传正文转写相关段落 | p304及p306 §3.1；原PDF/SI未取得，表格未视觉核验。不采用产量数字或由转写表推定等基因对照 |
| S2024及S2024_SM，工程株单/双缺失 | [10.1126/sciadv.adn0414](https://doi.org/10.1126/sciadv.adn0414)；机构全文、官方XML/SI | Fig2 p3、Discussion p9、Methods；S3 p17/S5 p23视觉核对。SI与正文是一项研究，不是两次独立验证 |

9篇论文中，6篇取得全文PDF，2篇仅正式摘要，1篇取得作者上传正文转写的相关段落；2024补充材料另计文件，不另计论文。并非所有全文的每幅图都经过视觉检查，实际范围如上。PDF原文件仅存于本任务work目录并登记SHA，不作为可自由再分发的报告附件。

## 两案对决策的影响

| 案件 | 本轮可用结论 | 仍不成熟的步骤 | 最小可审阅成果 |
|---|---|---|---|
| EGC-4b1e970207ef | 两种主要原生裂解活性可区分；ICL约3%弱交叉活性必须保留；2024有两个单缺失及双缺失对照 | 排他地把活性分配给当前两个蛋白，验证当前458 aa形式、条件性补偿和区室 | R490→ICL1候选、R552→ICL2候选的主功能证据表；保留共享OR作现状对照，不生成硬删OR补丁 |
| EGC-aadd98455a18 | 1981第一步及1995第二步有原生酶学；第二步专属活性和aconitase相关活性并存 | 第二GPR、当前位点直接功能/区室、立体化学、m217和m76连接 | 两半反应形式计量、竞争解释和化学整理提案；不新建空GPR反应或复用有冲突的中间体 |

共同化学发现是：**m217名称为甲基异柠檬酸，其完整InChI与结构键却精确对应甲氧基三羧酸同分异构体。** m76结构标识支持甲基顺乌头酸，但名称和既有homocitrate/homoisocitrate邻域不一致。均是当前字段冲突；三条反应按公式/电荷仍配平。不能从“配平”或“结构键相同”直接接通网络。[确切PubChem JSON](https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/15788107,5459784,23615311,23615274/property/IUPACName,InChI,InChIKey,MolecularFormula,Charge/JSON)

R95形式上可拆为 `MC³⁻→I³⁻+H₂O` 和 `I³⁻+H₂O+3H⁺→MIC⁰`，两步配平并恢复现有总式。这个3H⁺来自端点电荷形式不同，不是跨膜质子泵证据。该形式推导不决定第二GPR、方向界限或生理可用性。详细字段、源定位与反证均在案件文件和 `MC26-CH-01..03` 中。

## 跨来源张力与非独立性

以下为针对本案的关系清单，`scholar_confirmation=pending`；未声称穷尽所有论文两两比较。

| 关系 | 审阅判断 | 对结论的约束 |
|---|---|---|
| 1976主要特异性 vs 弱交叉活性 | 同文直接结果并存，无需二选一 | 支持主功能差异，但不支持严格零交叉；低体外比例不能自动折算通量 |
| 1981后半反应未分清 vs 1995分离专属活性并观察aconitase活性 | 后续部分解决早期不确定性，不是简单矛盾 | 第二活性存在已加强，现代第二GPR仍未解决 |
| 1982MICL双定位 vs 1996ICL分级结果 | 酶对象、碳源与检测方法不同，重叠不足 | 不把两文拼成同一当前位点的互斥定位结论 |
| 1993蛋白摘要 vs 当前版本 | 555/541/540/458 aa不可合并；修订与截短需分别追溯 | 全长实验不能未经验证完整移给当前458 aa蛋白 |
| 2013酶学叙述 vs 更早论文 | 一部分为引用前人、未发表分级或预测 | 不算新的独立纯化酶/受控定位复现 |
| 2024单/双缺失 vs 历史底物实验 | 工程产物终点与体外催化终点不同 | 补全旧报告只讨论双缺失的范围；不能用产量直接判定底物特异性 |
| UniProt/Rhea注释 vs 所引1976–1995论文 | 同一原证据的数据库整合 | 不把数据库与其原文累加成独立重复；旧生化到序列的归属另查 |
| 两个本任务案件 | 共享路径、m217和若干文献 | 不能互为独立验证数据 |

本轮没有发现足以支持在当前准确SD-Leu条件下进行排他GPR修订的完整链条；这不是“证据证明不能拆分”，而是提案接纳条件尚未满足。

## 验证顺序、失败和停止

各案件已列出亲本、两个单缺失/双缺失、准确序列回补、底物纯度与失活空白、直接产物鉴定、表达与定位/分级控制等判别性对照。正式实验的样本量、剂量和接受阈值仍由负责人决定，本任务不补造参数。

1. 先独立核查决定性原文与身份版本；对摘要/转写限制保留降级，不等待不可得来源而扩大范围。
2. 核查m217和m76化学身份、端点立体形式及全部直接连接。若需要任意运输、空GPR自由通路或混合同分异构物种才能连通，提案失败并停止实施。
3. 分别解决主要/交叉底物功能、第二步现代GPR以及同条件区室。若当前蛋白显示生理相关交叉功能或可直接催化总反应，修订排他/两酶假设；若只有同源、域名、产量或粗提物数据，则限于间接证据。
4. 经独立审阅后，负责人决定是否接纳具体科学变更及限定计算/实验范围；必须使用新版本、保留原始输入和calls。输入SHA变化、预算到限、超出本轮或只为提高essentiality匹配率改变GPR时停止。

当前缺口主要**限制科学声明与后续实施**，不阻塞本轮资料交付。无关案件、脂质路线、广泛BLAST/结构搜索、全模型审计和新一轮计算未开展。

## 技能、阶段与未执行检查

下表列出实际读取并采用的文件；角色提示以本任务内的阅读/核对/综合阶段应用，未创建相互独立的代理。

| 完整路径 | 实际采用阶段 | 未执行或受限部分 |
|---|---|---|
| `/Users/david/.codex/skills/govern-agentic-research/SKILL.md` | 中性问题、限定计划、来源直接性、反证、可追溯台账、停止条件 | 独立来源审计未执行；科学变更、实验和发布接纳仍需后续人类决定 |
| `/Users/david/.codex/skills/gene-identity-function/SKILL.md` | 当前/旧ID、名称、功能与证据等级分开，接受冻结身份依赖 | 未做新序列搜索或将数据库桥接升为催化/定位证明 |
| `/Users/david/.codex/skills/academic-research-suite/SKILL.md` | 文献研究工作流路由与证据综合 | 未运行完整ARS、审稿人模拟或任何hook |
| `/Users/david/.codex/skills/academic-research-suite/ars/deep-research/WORKFLOW.md` | 限定检索、来源核对、主张意向记录后写作 | 不是穷尽系统综述或完整论文生产流程 |
| `/Users/david/.codex/skills/academic-research-suite/ars/deep-research/agents/bibliography_agent.md` | DOI/题名/访问类型/原文定位整理 | 期刊声誉、全量撤稿与利益冲突审核未完成 |
| `/Users/david/.codex/skills/academic-research-suite/ars/deep-research/agents/source_verification_agent.md` | 区分原文、摘要、数据库与继承；精确主张核对 | 由同一任务作者执行，不能称独立审阅 |
| `/Users/david/.codex/skills/academic-research-suite/ars/deep-research/agents/synthesis_agent.md` | 写作前主张意向、跨论文张力与条件差异 | scholar确认pending；未做全部论文的穷尽两两审计；临床I–VII证据层级不硬套酶学 |
| `/Users/david/.codex/plugins/cache/ponytail/ponytail/4.9.0/skills/ponytail/SKILL.md` | 最少标准库脚本完成哈希、XML、结构字段和输出完整性检查 | 无新增依赖、框架、求解器或模型测试矩阵 |
| `/Users/david/.codex/plugins/cache/openai-primary-runtime/pdf/26.904.11930/skills/pdf/SKILL.md` | PDF文本提取、扫描页/目标图表视觉核对 | 1982/1993无全文；2013无原PDF；1996部分字体渲染受限，未声称全图完美读取 |

适用指令检查包含全局空文件 `/Users/david/.codex/AGENTS.md`、当前目录及祖先、项目 `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/AGENTS.md`；另外沿当前/项目/参考模型祖先检查13个 `AGENTS.override.md` 候选位置，未发现文件。其他工作树规则没有被自动当作本目录授权。项目状态及baseline_manifest只读取相关入口，不将其中历史命令当作执行指令。

未执行清单：新BLAST/结构/靶向预测、实验重复、当前准确培养条件验证、运输或容量推定、全基因组gRNA唯一性、本轮独立审计、撤稿/期刊声誉/COI完整审计。数据或未知值没有用默认参数补齐。

## 审阅覆盖、验证及AI披露

**独立审阅口径：总主张27｜已审0｜支持0｜未解决27｜冲突0｜未查27。** 其中21条被标记为决定性主张，独立已审0。这里“支持0”表示尚无独立审阅结论，不表示没有支持性原始证据；未查包含在未解决中，不另相加。

作者来源自检另计：`supported=16`、`partially_supported=11`。这些状态在主张台账中保留为 `verdict_author=this_task_author_self_check`，全部 `audit_status=unchecked`，无伪造审阅者或审阅时间。前轮独立审阅及身份作者自检不能自动覆盖本轮新增/改写主张。

实际完成的可重复静态检查：13项输入哈希、3条目标反应元素/电荷、3行目标映射和身份绑定、7个直接邻接反应、2个形式半反应配平、m217/m76确切结构键比较。交付检查进一步核对主张唯一性、来源引用完整性、登记本地材料SHA、案件ID与审计计数，并检查原模型哈希仍匹配。检查代码与下载/渲染材料保存在本任务work目录；[输出核验脚本](/Users/david/Documents/Codex/2026-09-06/iyali26-methylcitrate-research/outputs/verify_outputs.py)可只读重放交付检查。脚本通过不等于生物学正确或项目测试通过。

AI使用披露：本报告由Codex语言模型完成资料检索、论文阅读、图表核对、静态解析、证据综合与写作；未生成实验数据，未把作者推断写成测定，未以代理一致性替代来源审计。身份研究为协调提供的冻结依赖；本轮最终科学判断、独立来源审阅及变更接纳仍由研究负责人负责。未执行的检查均如上保留。

## 时间与冻结

首轮开始：2026-09-06 22:34:09 UTC（洛杉矶15:34:09）。预算截止：2026-09-06 23:19:09 UTC；本轮未开启新一轮。

最终完成时间：`2026-09-06 23:12:29 UTC`（首轮资料与文件冻结；具体检查时间以verification.json为准）。交付校验摘要见[verification.json](/Users/david/Documents/Codex/2026-09-06/iyali26-methylcitrate-research/outputs/verification.json)；最终文件哈希见[SHA256SUMS](/Users/david/Documents/Codex/2026-09-06/iyali26-methylcitrate-research/outputs/SHA256SUMS)。
