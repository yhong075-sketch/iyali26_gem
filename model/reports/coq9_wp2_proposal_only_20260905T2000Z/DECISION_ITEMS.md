# CoQ9 WP2：提案决策与人类闸门

日期：2026-09-05。交付状态：**proposal_only / awaiting_human**。

Material Passport：Origin Skill = experiment-agent / govern-agentic-research；Origin Mode = plan；Verification Status = UNVERIFIED（后续设计尚未实施/测试）；Version Label = proposal_v1。

## 1. 本轮完成的是什么

本轮授权是四份提案及文件身份清单。我们只读取交接包、核对相关原始文本/字节身份并撰写新文档；没有运行封存审计脚本、合成测试、FBA/dFBA/FVA、HPCC 作业，也没有修改执行代码、模型、培养基、参数、原始轨迹或旧结论。

四份提案分别为：

- [参数捕获](PARAMETER_CAPTURE_FIX_PROPOSAL.md)：14 个丢失值的路径、恢复边界、无损类型表示与实际求解对象/阶段绑定。
- [零事件比较](ZERO_EVENT_POLICY_V2_PROPOSAL.md)：连续量、raw exact-zero、near-zero 标签与未来行为验收分层。
- [重放计划](WP2_REPLAY_V2_PLAN.md)：三条停止诊断与八条新证据链的不同目的、条件与资源上限。
- 本文件：事实等级、未知事项和逐项待审批范围。

**建议顺序是：先审代码/比较规则提案，另行批准最小实施与合成测试；验证后再单独决定新计算范围。** 本轮不预先批准这些步骤，也不生成可运行作业命令。

## 2. 来源身份与本次核查范围

输入为用户提供的 `coq9_wp2_review_and_next_step_handoff_20260905.zip`；SHA-256：

`3d732a67cf3f90a9d26213cf2a160db491406f7492d855cf1ce06c42acef9001`。

内嵌原 WP2 ZIP：

`0b184c04a3a1270f413c4a5fbba60d386eedf9ea2560b1b5796a1fe905d90e16`。

CoQ9 冻结模型：

`bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee`。

冻结 compute commit：`36bb6f0735e4c6458bd53c0ceb01952b116b8be7`。这些身份不能被当前工作树、main 或另一任务的模型替代。完整来源成员 SHA 与本提案文件身份见 `document_manifest.json`、`SHA256SUMS`。

核查分三层：

1. **本轮直接核对**：输入 ZIP SHA、外层索引 28/28；相关参数/停止表、环境 JSON、sidecar 与 frozen runner 的文本、选择的执行原始日志。没有运行它们。
2. **交接完整性核查**：同一交接包已有独立核验覆盖外层 28/28、展开审计索引 21/21、22 个展开文件与内嵌审计 ZIP 相同、两份内嵌原 ZIP 身份一致。该独立检查明确没有重算所有 LP/轨迹。
3. **引用封存的数值审计**：1760 个状态、完整 flux 与账本等结果来自 A 的独立只读审计材料。本轮没有重做这些数值审计，不能写成本次复现。

H/A/B/F 等归档别名在 `document_manifest.json` 定义；下表的源码行号针对归档内版本。

## 3. 原子事实台账

`supported` 表示在指定证据与范围内支持；`partial` 表示仍有限制；`unverified` 表示未建立；`contradicted` 表示相反证据存在。尤其不能把对**封存审计报告的准确引用**提升为本轮已重新计算。

| ID | 具体判断 | 状态 | 包内定位；核查层级/限制 |
|---|---|---|---|
| F01 | 外层包及两份内嵌原 ZIP 的身份与 handoff manifest 一致 | supported | H `HANDOFF_MANIFEST.json`、`SHA256SUMS`；直接字节核对及独立完整性检查 |
| F02 | 8 条件保留 3 条可用归档可观测量相符、2 条仅单条件终点相符而未检验动态比值、3 条 STOP 的分层 | supported | A `condition_verdicts.tsv` 全8行、`historical_comparison_check.json`；直接读取分层，不将其称8条全部通过 |
| F03 | 三个 STOP 的对应连续库存差异都低于原 atol=1e-10 mmol/L | supported | A `exact_zero_points.tsv` 第2–4行；B `analysis_code/coq9_wp2_replay.py:31–41, 99–130` |
| F04 | “三处原始 exact-zero 布尔值也相同” | contradicted | A `exact_zero_points.tsv` 第2–4行；原值与旧 STOP 保留 |
| F05 | “停止之后的全部轨迹/终止/calls 已证明等价” | unverified | A `historical_comparison_check.json`、`condition_verdicts.tsv`；2305参考行减1760已归档行=545未执行行，无新证据补足 |
| F06 | 当前 environment 中151项过滤后执行参数有14项 `current_value=null`，序列化把非有限float合并为null | supported | B `historical_vs_current_environment.json#/historical_vs_current_parameters`；B replay `52–84, 138–147, 503–526`；参数提案逐项列14路径 |
| F07 | “这些14项已无损保存/可以依靠旧fingerprint恢复符号与类型” | contradicted | B replay `74–75, 259, 542, 569–570`；A `null_parameter_inventory.tsv`；有损投影相同不等于原值相同 |
| F08 | “序列化缺陷已经证明把null作为参数重新设置进求解器” | unverified | A `evidence_status.tsv` 的 `serialization_defect_directly_reconfigured_solver`；B replay `43, 503–508` 未发现该反馈路径。记录缺陷不等于已证明计算污染 |
| F09 | “14项可以从本次已检同次、同对象、同阶段证据逐项恢复” | unverified | 参数提案第3节及其日志/预检/JSON定位；在已搜索范围未找到合格旁证，不宣称世界上永久不可恢复 |
| F10 | 14项记录缺口涉及全部8条件，而不仅是3个STOP；补3条不能补全另外5条旧记录 | supported | B replay `495–542, 569–570`；A `INDEPENDENT_REVIEW_zh.md` 第5节；共享有损记录及指纹机制 |
| F11 | 封存独立审计报告中的状态算术、flux/账本在原声明闸门内一致 | supported | A `state_arithmetic_summary.json`、`flux_table_check.json`、`ledger_table_check.json`；**引用既有审计，不是本轮全量重算**，不证明化学/原生功能或最优解唯一性 |
| F12 | “全部 bound excess 都不超过1e-9” | contradicted | A `state_arithmetic_summary.json`：最大3.846612259402957e-9；仍低于原审计闸门1e-7。不以solver设置等同归档残差上限 |
| F13 | 原始本地实验工作簿可核对18条记录、90个单元格，但出版源、guide QC、映射与原生功能还有缺口 | partial | A `local_workbook_checks.json`、`evidence_status.tsv`；本地表格一致不等于这些外部证据全部验证 |
| F14 | “exact-zero 差异证明生物学回归，或infeasible证明细胞死亡” | unverified | A `exact_zero_points.tsv`、`condition_verdicts.tsv`；不足以推出这些解释 |
| F15 | frozen runner 已有 Q9 reserve 小值清零；source-free infeasible 原分支不会自动测试reserve-enabled rescue | supported | F `scripts/gem_annotate/quinone_dfba_essentiality.py:438–464, 576, 642–648`；本轮直接读代码，不改变机制 |
| F16 | “本文提出的参数挂钩、near-zero策略或未来新8条已实施并验证” | unverified | 本轮仅文档授权；两个方案的 synthetic tests 均为设计。现有任何历史授权不等于本次新版本计算授权 |

本表16条：supported 7、partial 1、unverified 5、contradicted 3。这里是**声明等级计数**，不是计算成功率，也不是审查者投票。

## 4. 事实之外的推论与候选设计

| 类别 | 判断 | 能支持到哪里 | 不能升级为什么 |
|---|---|---|---|
| 推论 | 原STOP可能是对极小残值敏感的比较边界，而不一定是生物学变化 | 三个连续差异小且raw bool不同，为重设计比较语义提供动机 | 不能证明未执行tail相同；不能改旧STOP |
| 推论 | 三条补跑与八条新记录应区分 | 14项缺口存在于8条件；3-only缺少新的nonlimiting WT | 不能因此宣称八条已获运行批准 |
| 设计 | 参数值采用带类型、非有限符号与位型的无损表示，分离值指纹与捕获事件指纹 | 为未来可审计记录规定契约 | 不能恢复旧null或证明真实挂钩无副作用 |
| 设计 | near-zero只在比较层新增标签；由既有库存/通量预算与B_start×dt限定 | 有明确量纲和预先冻结方案；不向模型回写 | 不是新生物阈值、库存clamp或未来行为证明 |
| 设计 | 两个WT→两个其余旧STOP→四个剩余KO，总计8条件 | 若目的为新完整记录，可设逐批闸门降低无效计算 | 不是54条dt、FVA或其他反事实授权 |

全部候选值、阶段和预算仍待用户批准。不是以recall、三个误差大小或“最好看的结果”选择方案。

## 5. 分开的审批项

| 决策 | 建议范围 | 当前状态 | 批准也不包含什么 |
|---|---|---|---|
| D0 本轮文档 | 四份proposal及来源/校验清单 | **已授权；仅此项** | 无代码、测试、求解或外部发送 |
| D1 参数记录代码 | 独立分析副本中现有replay/capture的最小修订；先核验真实依赖源码与挂钩位置；不改frozen compute | awaiting_human | 不恢复旧null、不运行solver smoke、不改变科学输入 |
| D2 比较层policy | 审批ZERO文档的公式、常数、单位、strict gates、版本；仅新比较层 | awaiting_human | 不改库存/bounds/Euler/原calls，不把旧STOP改PASS |
| D3 合成测试 | 人工标量/fake对象的类型往返、失败、阶段与零事件检查；无cobra/gurobi求解 | awaiting_human；实施后才可报告测试结果 | 不含真实模型/求解器观察、真实八条件或HPCC作业 |
| D4 新计算选择 | D1–D3验收后，单独批准3-only诊断或8条件新证据链，并批准资源/停止规则 | awaiting_human | 不含额外WT、dt/FVA/其他反事实；不自动重试/调参 |
| D5 新计算执行放行 | 输入/代码/实际设置/策略SHA在新run前固定，满足WP2_REPLAY_V2_PLAN中的门；保留原四项显式设置，OptimalityTol=1e-7及Presolve=0仅增加读回核对 | awaiting_human；不能由本提案自放行 | 不将历史unknown改为已恢复，不改旧manifest；已知值不符就停，不新增setParam自动纠正 |

D1、D2、D3可在用户明确列明范围后一起批准，但它们与D4/D5仍是不同权限。不得把一句对文档“可以”扩展为新计算授权。任何正式GPR、chemistry或EGC接受仍遵循原独立人类闸门，本次不涉及。

## 6. 建议用户下一次选择的证据目标

若目标是**先看原三个STOP之后是否出现实质分歧**，可在修复与合成测试通过后考虑严格3-only方案。代价是：另外五条的完整参数记录不会被补齐；`po1f_nonlimiting`没有新的WT，相关新KO不能计算新的WT归一化比值。

若目标是**为原8个条件建立统一、无损记录的下一版基线**，建议条件性选择八条新记录，按2+2+4分批。但每条都只能叫“新声明设置下的受限等价检查”；完整历史环境复现仍 `not_established`。原生生物参数、完整历史flux和缺失证据不会因此自动补足。

本文件不替用户做计算范围选择；先等D1–D3决定，之后再申请D4/D5。新的WT分母必须来自同一获准新run、同mode、相同参数/输入，并通过适用闸门；不能借用旧WT拼出“本轮验证比值”。

## 7. 仍未解决的事项及影响范围

| 问题 | 影响分类 | 本轮/后续处理 |
|---|---|---|
| 14个旧参数原值仍无合格同次旁证 | 限制历史完整性声明；影响下一版参数捕获设计 | 保留unknown/null；新无损记录只能解决新run的可审计性 |
| 原依赖内每次native solve的最小挂钩位置未核实 | 阻塞“完整实际阶段捕获”实施验收 | 代码授权后先固定并读取对应依赖源码；不能用外层pFBA返回近似冒充 |
| 三处STOP之后545参考行未在旧WP2执行 | 限制历史等价声明；需另批计算才能观察新轨迹 | 不生成尾段、不离线补PASS |
| r2没有全部历史库存、全flux/LP、pFBA次级目标 | 限制可比字段范围 | 只比实际存在字段；新全状态自洽不等于老全状态相同 |
| near-zero候选与参数codec尚未实现/合成测试 | 阻塞新版本验收，不阻塞本轮文档 | 无通过率、无实现完成声明 |
| 非推进terminal row仍采用strict raw-zero门 | 限制“新规则一定消除三个STOP”的预期 | 原第三STOP后的终止行仍可能STOP；本方案不继承前一区间band，进一步放宽需另案审批 |
| 参数边界挂钩的真实无副作用性未验证 | 阻塞真实运行放行 | 未来只在另批小规模计算范围验证，不夹带solver测试 |
| 化学/呼吸链身份、原生蛋白与guide QC缺口 | 限制生物学解释；不自动影响独立记录修复 | 保留unresolved；不删除R1889、不为反应机械加H+或添加AND依赖 |
| artificial reserve不是实测催化CoQ pool | 限制机制/生物参数解释 | 继续sensitivity_only_not_calibrated；不按recall调alpha/pool/cutoff |

没有发现新授权可支持54条dt矩阵、FVA、maintenance/glucose反事实、模型修复、湿实验、main推送或正式病例接受。本轮交付后停止。

## 8. 独立来源审查记录

审查时间：2026-09-05。独立来源审查者：本任务的 `proposal_source_review` 子代理，与参数/重放文档作者分开；主代理根据审查意见修订后，审查者重新读取相关最终段。该角色没有修改四份文档、运行代码/测试或调用求解器。

实际覆盖：四份正文完整阅读；F01–F16的来源均直接打开核对。16条限定声明中审查16、unchecked 0；supported 7、partial 1、unverified 5、contradicted 3。**这不是候选实现/科学参数的验证覆盖率**；F11、F13只核对封存审计报告的实际声明与限制。

独立直接检查包括：

- 外层、两份内嵌ZIP、执行tar、冻结科学输入/compute、protocol/reference、sidecar及environment的SHA；
- 14个null的原字段、实际model读取路径、有损编码、共享fingerprint；原日志27,667行的逐项参数名检索与独立预检进程来源；
- 原3/2/3分层、三个STOP、545行未执行、原容差、Q9既有clamp与source-free流程；
- 三条/八条证据目标、新WT分母门、文档授权边界。

发现并解决的文档问题：

1. terminal row是否继承near-zero band原先不够明确：现明确不继承，raw-zero差异仍strict STOP；不保证三个STOP都消失。
2. 参数提案保留四项设置，而重放草案曾建议额外显式设置两项：已删除额外setParam建议；仅对OptimalityTol/Presolve读回核对，不符停止。
3. 库存驱动通量差异何时STOP原先措辞过宽：现明确按对应原容差/离散闸门判断，不把所有微差都称分歧。
4. F10的独立报告节号修正为第5节，ZERO的三处表名补为精确locator。

以上四项最终修订均经独立复读，未发现剩余阻塞**文档交付**的问题。

未执行/未覆盖：1760状态和数百万通量的重新算术审计；外部r2/r3原ZIP及出版商工作簿独立认证；原生蛋白功能；真实依赖native挂钩；候选codec/near-zero的任何测试；任何新求解。完整source-audit结论仅为“可以作为proposal-only材料交付”，不是实施/测试/计算验收或下一工作包授权。
