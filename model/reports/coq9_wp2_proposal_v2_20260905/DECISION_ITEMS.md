# CoQ9 WP2：提案决策与人类闸门

起始日期：2026-09-05；文档修订日期：2026-09-06。交付状态：**proposal_only / awaiting_human**。

Material Passport：Origin Skill = experiment-agent / govern-agentic-research；Origin Mode = plan；Verification Status = UNVERIFIED（后续设计尚未实施/测试）；Version Label = proposal_v2；document_schema_version = 2。

## 1. 本轮完成的是什么

本轮授权是根据 `TO_INACTIVE_NEXT_STEPS_zh.txt` 对R1–R6进行文档定点修订，提交四份proposal_v2、逐项回应及新文件身份清单。v1原样保留；我们只读取材料、核对相关文本/字节身份并撰写新文档；没有运行封存审计脚本、合成测试、FBA/dFBA/FVA、HPCC 作业，也没有修改执行代码、模型、培养基、参数、原始轨迹或旧结论。

四份提案分别为：

- [参数捕获](PARAMETER_CAPTURE_FIX_PROPOSAL.md)：14 个丢失值的路径、恢复边界、无损类型表示与实际求解对象/阶段绑定。
- [零事件比较](ZERO_EVENT_POLICY_V2_PROPOSAL.md)：连续量、raw exact-zero、near-zero 标签与未来行为验收分层。
- [重放计划](WP2_REPLAY_V2_PLAN.md)：三条停止诊断与八条新证据链的不同目的、条件与资源上限。
- 本文件：事实等级、未知事项和逐项待审批范围。

**建议顺序是：分别接受比较政策、批准对应代码实施及synthetic tests；另行批准真实对象零求解预检取得新manifest；最后另行决定有限计算范围和逐批执行放行。** 本轮不预先批准这些步骤，也不生成可运行作业命令。

## 2. 来源身份与本次核查范围

输入为用户提供的 `coq9_wp2_review_and_next_step_handoff_20260905.zip`；SHA-256：

`3d732a67cf3f90a9d26213cf2a160db491406f7492d855cf1ce06c42acef9001`。

内嵌原 WP2 ZIP：

`0b184c04a3a1270f413c4a5fbba60d386eedf9ea2560b1b5796a1fe905d90e16`。

CoQ9 冻结模型：

`bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee`。

冻结 compute commit：`36bb6f0735e4c6458bd53c0ceb01952b116b8be7`。这些身份不能被当前工作树、main 或另一任务的模型替代。完整来源成员 SHA 与本提案文件身份见 `document_manifest.json`、`SHA256SUMS`。

核查分四层，不能互相替代：
 
1. **previous source audit（v1）**：旧DECISION_ITEMS§8记载本地子代理对来源的独立核查；旧manifest和SHA保留。该记录不能变成CoQ9聊天的审核经历，也不能不加说明地变成本轮全部重核验。
2. **CoQ9文档审查**：附件开头明确总体partial，只基于逐份全文与该聊天已有来源，未收到新的提案ZIP、未独立核验该ZIP或六份文件的字节SHA。审查源文件SHA为 `a3eafeb8becbd480f59cb39f22e70fbc7b384cc645edbf319885f4732cf6cc41`。这是来源/权限范围事实，不是用户计算许可。
3. **本轮文档修订核查**：本地核对v1六文件/提案ZIP及附件身份，读取必要原代码/环境字段/参考step85；只检查R1–R6的文本修改与文件身份。实际覆盖由本文件§8和新manifest列明，不重做全量状态审计。
4. **引用封存数值审计**：1760状态、flux/账本、原工作簿单元格结论来自A的封存报告；本轮没有重做这些数值审计或出版源认证。

本次使用的PROJECT_STATE片段、当前工作区观察和来源全文SHA在CHANGE_RESPONSE的“来源与身份片段”及manifest中绑定。片段只支持“文件作了该陈述”，不证明其计算结果；CoQ9聊天未亲自核验这些本地字节，不称它已核验。

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
| F09 | “14项可以从既有已检同次、同对象、同阶段证据逐项恢复” | unverified | 参数提案第3节及其日志/预检/JSON定位；在已搜索范围未找到合格旁证，不宣称世界上永久不可恢复 |
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

## 5. R4/R6 统一授权与证据目标决策表

以下exact permission names与四份正文、manifest一致。**本轮只有document_revision已获授权**；七项后续权限都是awaiting_human，审批每项不隐含另一项。

| permission_name | 单独允许讨论的候选动作／证据目标 | 取得什么证据；不能证明什么 | 本轮状态 |
|---|---|---|---|
| document_revision | R1–R6定点修订四文档、回应及身份清单 | 仅文档交付，不是软件/科学验收 | granted_document_only |
| parameter_capture_code_implementation | 在获准独立副本修订参数codec和观察层；检查绑定依赖源码；不改frozen compute | 代码可供审核，不产生真实对象值/native覆盖或运行证据 | awaiting_human |
| comparison_policy_acceptance | 接受ZERO文档的规则/常数/适用表/terminal选择及新文档SHA | 政策决定；**不授权改代码** | awaiting_human |
| comparison_code_implementation | 在政策接受后，修改获准比较/报告副本及其声明文件 | 实现文本；**不授权运行测试** | awaiting_human |
| synthetic_test_execution | 仅人工标量/fake对象，验证codec契约和模拟分支；无真实solver | codec/fake检查证据；不能通过native_hook_coverage或runtime_noninterference | awaiting_human |
| real_object_parameter_preflight | 加载冻结输入及获准sidecar，创建实际对象，执行已列明的原四项配置路径、读取参数/依赖身份并输出resolved manifest；optimization_calls=0 | 该预检对象的API返回值；不是后续新对象实测值，也不是native solve覆盖 | awaiting_human |
| new_computation_scope | 选择严格3-only或原8条件及全部参数/比较/资源范围，不添加对照或参数扫描 | 确定未来证据目标；**不自动提交作业或开始任何batch** | awaiting_human |
| batch_execution_release | 接受具体批次的输入/代码/政策/新resolved manifest/预算，明确允许原有限轨迹中的native-hook观察与失败记录 | 该批实际状态与受限参考比较；不证明全部历史内部状态等价 | awaiting_human |

允许用户明确同时批准多项，但必须列明各项范围；例如接受公式不等于comparison_code_implementation，批准代码不等于synthetic_test_execution，批准scope不等于batch_execution_release。实际挂钩smoke若要额外求解或新条件，需要独立明确授权；本方案不包含额外smoke，正常hook观察仅可纳入之后获批的有限基线调用，不能增加solve次数。

### 5.1 零求解真实对象预检，不再让授权闭环悬空

该预检本身会加载库、创建真实solver对象，可能触发license/资源分配，因此**不能归入普通只读文件核查或synthetic tests**；本轮严格不执行。候选范围单进程最多两种mode各一个预检对象，optimization_calls=0；读取前后对象链、过滤参数值、预期/发现集合、依赖/输入/代码身份、实际配置差异，输出manifest，不调用optimize/pfba/slim_optimize或任何solver smoke。

预检可沿用已批准的原配置加载路径，显式仅Threads=1、Seed=0、Method=1、FeasibilityTol=1e-9；OptimalityTol=1e-7、Presolve=0只读回核对，差异即停。编码器/比较器不反馈setParam，不追加探测/默认值来配平。

默认候选预检进程在交付manifest后结束；随后每一计算批次创建自己的新对象并重新读取核对。旧预检manifest只是期待值/审批依据，不能直接标作新对象observed值。若希望同一对象等待人工批准并继续运行，必须明确选择其会话/资源生命周期，而不能隐含维持；本提案默认不保留。

独立候选预检预算：包含加载和记录的active wallclock最多120秒，单进程及子进程内存4GiB、预检输出16MiB；待real_object_parameter_preflight批准，不能挪用原八条件计算授权。失败/超限输出preflight_incomplete、optimization_calls实际可证范围和缺失项，不修补、不重试。实际执行前仍须确认原loader和获准sidecar没有隐藏优化路径；无法证明零求解则停在预检设计门。

### 5.2 三类验收分开

- `codec_contract`：只保证API**实际返回的受支持标量**类型和值/位型；合成测试才可能验证该契约，不承诺内部表示、自动算法选择、basis/hidden state。
- `native_hook_coverage`：需要固定实际依赖调用链和获准有限运行中真实native pre/post记录，fake并不能证明覆盖。
- `runtime_noninterference`：仅指数学输入、参数、原定调用和控制流未被观察层改变，需要实际有限运行的证据；记录/压缩的时间、内存、磁盘成本另披露，不承诺wallclock相同。

三项现均unverified。真实有限运行前可接受codec/静态挂钩设计作为阶段资格，但不能要求尚未授权运行才可产生的证据已PASS，也不能反过来以“实施验收”暗中运行真实smoke。运行中的native覆盖/非干扰证据产生后另行审查；缺证就保留不完整而不是补认证。

## 6. 建议用户下一次选择的证据目标

若目标是**先看原三个STOP之后是否出现实质分歧**，可在修复与合成测试通过后考虑严格3-only方案。代价是：另外五条的完整参数记录不会被补齐；`po1f_nonlimiting`没有新的WT，相关新KO不能计算新的WT归一化比值。

若目标是**为原8个条件建立统一、无损记录的下一版基线**，建议条件性选择八条新记录，按2+2+4分批。但每条都只能叫“新声明设置下的受限等价检查”；完整历史环境复现仍 `not_established`。原生生物参数、完整历史flux和缺失证据不会因此自动补足。

本文件不替用户选择计算范围；按§5逐项审批。新WT必须来自同一获准new_computation_scope计划、同mode、相同输入/设置并通过适用闸门，分母有限且>0。跨批引用必须携带denominator_run_id、condition_key、gate_status及WT/KO各自termination_time；终点倍增数比不是共同观察时窗指标，不能借旧WT拼出新验证。严格3-only没有新的nonlimiting WT；如补该对照至少4条，必须另行明确范围，不能隐含增加。

## 7. 仍未解决的事项及影响范围

| 问题 | 影响分类 | 本轮/后续处理 |
|---|---|---|
| 14个旧参数原值仍无合格同次旁证 | 限制历史完整性声明；影响下一版参数捕获设计 | 保留unknown/null；新无损记录只能解决新run的可审计性 |
| 原依赖内每次native solve的最小挂钩位置未核实 | 阻塞“完整实际阶段捕获”实施验收 | 代码授权后先固定并读取对应依赖源码；不能用外层pFBA返回近似冒充 |
| 三处STOP之后545参考行未在旧WP2执行 | 限制历史等价声明；需另批计算才能观察新轨迹 | 不生成尾段、不离线补PASS |
| r2没有全部历史库存、全flux/LP、pFBA次级目标 | 限制可比字段范围 | 只比实际存在字段；新全状态自洽不等于老全状态相同 |
| near-zero候选与参数codec尚未实现/合成测试 | 阻塞新版本验收，不阻塞本轮文档 | 无通过率、无实现完成声明 |
| 非推进terminal row仍采用strict raw-zero门 | 限制“新规则一定消除三个STOP”的预期 | 原第三STOP后的终止行仍可能STOP；本方案不继承前一区间band，进一步放宽需另案审批 |
| 参数边界挂钩的真实非干扰性未验证 | 限制运行完成验收和下一批放行；不是要求首批运行前已有运行证据 | 首批须先具备获批的codec/静态挂钩准备与明确批次授权；真实native覆盖及非干扰证据在该受限基线中形成，不夹带额外solver测试 |
| 化学/呼吸链身份、原生蛋白与guide QC缺口 | 限制生物学解释；不自动影响独立记录修复 | 保留unresolved；不删除R1889、不为反应机械加H+或添加AND依赖 |
| artificial reserve不是实测催化CoQ pool | 限制机制/生物参数解释 | 继续sensitivity_only_not_calibrated；不按recall调alpha/pool/cutoff |

没有发现新授权可支持54条dt矩阵、FVA、maintenance/glucose反事实、模型修复、湿实验、main推送或正式病例接受。本轮交付后停止。

## 8. v2修订审查记录与既有审查的区分

v1的previous source audit完整保存在未改动的旧DECISION_ITEMS§8。CoQ9对全文的审查总体partial且未核验新ZIP/字节；本轮不会替它升级为已验证。

本轮独立审查者为本地只读子代理 `v2_source_audit`，与参数提案及重放提案作者分开。2026-09-06完成四份proposal_v2及CHANGE_RESPONSE全文核对，六项结论均为 **document_revision_addressed / implementation_unverified**；未发现阻塞本轮文档交付的问题。这个结论只评价文档是否回应R1–R6，不是接受候选政策或科学/代码验收。

本轮实际独立打开的来源范围：附件全文；v1四文档及manifest；原环境14项null；B replay的31–147及495–570行；B capture的28–175行；冻结runner的438–464、576、642–648行；参考gzip及step84/85两行；condition和exact_zero表；封存evidence_status、state arithmetic、flux/ledger/workbook审计摘要；外层包及两内嵌ZIP身份。PROJECT_STATE的选定原文/完整SHA、当前branch/HEAD、目录模型SHA、构建起点缺失及旧提案ZIP/附件SHA亦经独立只读核对。

未重新审计1760个完整状态或全部通量向量，未恢复14个旧值；v1全日志及恢复缺口搜索只保留为既有审查，不称本轮重做。没有执行测试、封存脚本或solver。独立审查正文时新manifest尚未封存，因此最终文件SHA是另行进行的本地文件完整性核对，不能倒写为该时点已认证。设计合理性、真实对象值、native覆盖、运行非干扰与生物学结论分别保持各自状态。

本轮新明示的事实：
- **supported**：原参考nonlimiting末端step85存在，5.3125h非推进infeasible、glucose start/end均4.440892098500626e-16；locator见ZERO§6.1。
- **supported**：CoQ9给出partial文档审查且明确没有新ZIP/字节核验；来源为附件开头及其新manifest所列SHA。
- **unverified**：新运行是否出现terminal_raw_zero_only_mismatch或通过全部后续门；本轮没有新运行。
- **supported（授权边界）**：本轮明确授权仅document_revision；真实对象预检、native smoke、3/8条件、dt/FVA/其他反事实未获本轮授权，因此未执行。这不撤销历史授权记录，也不把它们扩展为新版本许可。

完成四份proposal_v2、CHANGE_RESPONSE及更新文件身份清单后停止。不自动转交外部聊天、不commit/push、不启动下一工作包。
