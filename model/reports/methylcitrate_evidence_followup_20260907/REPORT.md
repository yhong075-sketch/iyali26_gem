# 甲基柠檬酸两案：新增证据与修订结论

**两案资料整理完成，成熟模型变更提案为 0。** 新材料补充了碳源响应和构建来源联系，也明确了名称、化学身份及第二步酶归属的限制；目前仍不足以删除 R490/R552 的共享 OR，或把 R95 拆成可直接接纳的两个带 GPR 反应。本包只新增证据和修订说明，状态为 `proposal_only`。

## 固定对象与本轮用途

| 系统 ID／旧 ID | 名称、简要功能、证据等级 | 固定序列 | 当前模型角色 |
|---|---|---|---|
| YALI1C24124g／YALI0C16885g | ICL1；异柠檬酸裂解酶候选；当前短形式为 **model/GPR assignment only** | XP_065950166.2，458 aa | R490、R552 的 OR 成员 |
| YALI1F39620g／YALI0F31999g | ICL2（2024 靶向名；UniProt 无 GN Name）；2-甲基异柠檬酸裂解酶候选；**curated annotation** | XP_506117.1，565 aa | R490、R552 的 OR 成员 |
| YALI1F03803g／YALI0F02497g | PDH1（UniProt）／PHD1（论文）；2-甲基柠檬酸脱水酶候选；**curated annotation** | XP_504908.1，520 aa | R95 唯一 GPR，未证明覆盖整个净反应 |

沿用身份 v1 和最新[现行序列说明](../identity_construct_followup_20260907/CURRENT_SEQUENCE_NOTE.md)。历史长形式没有替换第一条输入；数据库功能注释不等于当前序列已被直接测定。

## 新证据改变了什么

**裂解酶案 EGC-4b1e970207ef。** 1992 年原论文摘要中的 ICL2 被作者鉴定为乙酰辅酶 A 合成酶结构基因，不能只凭同名归入当前 ICL2 候选。2025 年研究在解脂耶氏酵母 EXF-17398、氮限制和特定碳源条件下，观察到当前三个旧位点标签的转录响应；丙酸相对葡萄糖时，两种裂解酶候选都上调，不能据此推出底物互斥。2024 年 W29 的研究提供油酸条件下乙醛酸旁路的转录与通量背景，但 Fig. 6A 将两个裂解酶标签合并显示，不能分配单个蛋白的活性。[1992 摘要](https://pubmed.ncbi.nlm.nih.gov/1349449/)、[2025 原研究](https://doi.org/10.1186/s13068-025-02713-7)、[2024 原研究](https://doi.org/10.1016/j.ymben.2024.06.010)。

因此，主功能分配仍可作为实验候选；“交叉连接应删除”尚未成立。历史约 3% 的弱交叉活性继续保留为特定体外条件的观察，既不置零，也不换算成通量上限。详细条件、反证及删／留标准见[裂解酶案修订](EGC-4b1e970207ef.md)。

**脱水酶案 EGC-aadd98455a18。** 2018 年论文报告缺失株在丙酸唯一碳源条件下不能生长，支持该实验背景下的通路依赖；其 Fig. 3 图注对第二步乌头酸酶归属使用的是“probably”，不是新完成的第二步蛋白鉴定。该文沿用与 2013 年有关的具名构建。最新身份补查增加了引物与当前靶位点的局部对应，但没有恢复历史实验等位基因的完整序列。[2018 原研究](https://doi.org/10.1186/s13068-018-1154-4)、[身份补查报告](../identity_construct_followup_20260907/REPORT.md)。

本轮进一步定位了 m76 冲突的另一端：名称中的 but-1-ene-1,2,4-tricarboxylic acid 骨架，与顺式高乌头酸对应；现有结构注释却对应甲基顺乌头酸。两种三阴离子的分子式、电荷相同，原子连接不同；名称本身未写 Z，不能仅由名称确定立体形式。这也说明形式配平不能解决共享中间体的身份问题。[ChEBI 58174](https://www.ebi.ac.uk/chebi/CHEBI:58174)、[ChEBI 57872](https://www.ebi.ac.uk/chebi/CHEBI:57872)。

R95 应继续分开讨论第一脱水、第二水合与净反应。专属第二步酶和乌头酸酶相关活性都保留；m217 的既有名称／结构冲突也没有自动修复。详细化学解释及拆分验收条件见[R95 案修订](EGC-aadd98455a18.md)。

## 已完成的结构支持如何使用

结构任务已交付两份与固定目标序列精确对应的既有 **AlphaFold 预测**：AF-Q6BZP5-F1 v6、AF-Q6C354-F1 v6。其报告从逐残基文件重算的平均 pLDDT 分别为 92.92、92.52，并保存完整 PAE。这里引用的是结构任务的交付与作者分析，本任务没有重跑结构分析，也没有代替其独立审阅。[结构交付报告](/Users/david/.codex/worktrees/6033/iyali26_gem/artifacts/alphafold_methylcitrate_support_20260907/REPORT.md)。

这些结果可用于后续检验裂解酶家族及 PrpD 家族的功能候选；不能验证底物特异性、第二步催化、区室或体内补偿。所有据此形成的功能解释均标为“基于 AlphaFold 预测的功能候选”。

收尾期间又收到固定 **458 aa 目标的 AlphaFold 3 预测**，因此当前三条固定序列的结构交付均已取得：两条既有模型精确复用，一条本次网页预测。结构作者报告样本0的 pTM 为0.95、Cα平均pLDDT为96.375，并指出候选催化C121的局部置信度有限；与全长比较对象的预测相似性不能确认原生酶活或底物特异性。本任务仅接收并引用该补充，没有重新分析结构或改变删／留OR的标准。[AF3补充报告](/Users/david/.codex/worktrees/6033/iyali26_gem/artifacts/alphafold_methylcitrate_support_20260907/AF3_ADDENDUM.md)。

首次接收时第三条仍待结果的事实保留在后台接收记录。结构包的新增／修订7条记录继续为 `unchecked`，不并入本包14条或旧69条，也不视为独立功能验证。

## 核查范围、审计与完成数

本轮按新授权在 20 分钟内整理上一轮已取得的来源，仅补必要的原文／表格定位核对。实际使用结构化 XML、已保存 PDF 与图像、官方数据库记录及原始补充工作簿；三靶位点的表格摘录按原单元格回查。没有重跑差异表达、模型解析构建、代谢求解、序列筛查、预测作业或湿实验。

| 范围 | 总数 | 已独立审计 | 支持 | 未决 | 反驳 | unchecked |
|---|---:|---:|---:|---:|---:|---:|
| 原累计 claims，沿用最终审阅 | 69 | 68 | 66 | 1 | 2 | 1 |
| 本包新增 claims | 14 | 0 | 0 | 14 | 0 | 14 |

原甲基柠檬酸 27 条及已接受的 11 条补正不重审、不重编号；0/14 是本包的独立审计覆盖，作者来源核对不计入此列。身份包和结构包各自的新主张仍以其独立账本为准，不在这里并入旧 69 条。历史覆盖 322/1612、交集内召回 67/322、原始标签、冻结合同与历史 STOP 均保持原样。

**完成：2/2 案修订、14 条新增主张及来源定位；成熟变更 0、模型／GPR 修改 0、求解器调用 0。** 未完成项为新增主张独立审阅、历史构建精确闭合、当前蛋白底物与定位实验、第二步现代蛋白归属及化学修订的正式接纳。1992 年仅取得正式摘要；2018 年本轮核对出版 XML 的正文、图注与表格，出版 PDF 图形未完成目视核对。未重新核查全部撤稿、利益冲突和期刊数据库信息，不宣称全面来源认证。受限检索未找到新的固定序列直接酶学证明，不代表这类证据不存在。

本轮实际采用的技能：

- [govern-agentic-research](/Users/david/.codex/skills/govern-agentic-research/SKILL.md)：范围、证据等级、反证、停止条件和独立审计分开。
- [gene-identity-function](/Users/david/.codex/skills/gene-identity-function/SKILL.md)：系统 ID、名称、功能与模型角色分开。
- [academic-research-suite](/Users/david/.codex/skills/academic-research-suite/SKILL.md)：仅 fact-check／[source_verification](/Users/david/.codex/skills/academic-research-suite/ars/deep-research/agents/source_verification_agent.md) 的来源核对方法；没有启动其他角色或完整流水线。修订结论按本次明确授权交付。
- [PDF](/Users/david/.codex/plugins/cache/openai-primary-runtime/pdf/26.905.11957/skills/pdf/SKILL.md) 与 [Spreadsheets](/Users/david/.codex/plugins/cache/openai-primary-runtime/spreadsheets/26.905.11957/skills/spreadsheets/SKILL.md)：沿用 PDF 页码／图形核对和只读原始工作簿；不创建新 PDF 或工作簿。
- [Ponytail](/Users/david/.codex/plugins/cache/ponytail/ponytail/4.9.0/skills/ponytail/SKILL.md)：复用摘录和标准解析，仅检查实际交付数据，没有新增合成测试或工程框架。

可直接使用本报告及上述两份案件修订；独立审阅入口为 [claims.jsonl](claims.jsonl) 与 [source_records.json](source_records.json)。输入、来源快照、实际核查及完成时间保存在 [verification.json](verification.json)，本任务没有写共享状态。
