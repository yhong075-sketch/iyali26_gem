# CoQ9 WP2 proposal_v2：R1–R6逐项修改回应

起始日期：2026-09-05；文档修订日期：2026-09-06。状态：**document_revision_completed / proposal_only / awaiting_human**。这不是代码、合成测试、真实求解器操作或计算放行。

审查输入：`/Users/david/Downloads/TO_INACTIVE_NEXT_STEPS_zh.txt`，SHA-256 `a3eafeb8becbd480f59cb39f22e70fbc7b384cc645edbf319885f4732cf6cc41`。CoQ9聊天给出的总体判断为 **partial**；它审阅了逐份转交的全文及自己已有来源，**未收到或核验新ZIP及六份提案文件字节**。本地来源审查和文件SHA核对不归功于该聊天；本轮也没有向它自动发送文件。

以下“已落实”只表示新文档中有明确对应条款，**不是候选实现已通过验证**。事实等级仍以DECISION_ITEMS§3、§8的 supported / partial / unverified / contradicted 为准；设计实施与实测状态均未建立。

## R1：参数记录的保证范围

| 子项 | 文档改动与定位 | 保留意见／未解决事项 |
|---|---|---|
| R1.1 无损保证边界 | PARAMETER_CAPTURE_FIX_PROPOSAL§4、§8.1：限定为实际API返回的受支持标量类型和值/位型；不包括内部表示、自动选择、basis或隐藏状态 | 已落实于文档；codec尚未实现或测试 |
| R1.2 原样记录API值 | 参数提案§3–4、§6：有限值仍为有限值，不按参数名推成Infinity；旧null保持legacy_lossy_unknown | 14个旧值仍未恢复，新值不能回填旧记录 |
| R1.3 三项验收分离 | 参数提案§8.1；DECISION_ITEMS§5.2：分别定义codec_contract、native_hook_coverage、runtime_noninterference | synthetic/fake证据不能证明真实native覆盖；真实运行证据只能在另行获准有限运行中形成，不能构造循环前置条件 |
| R1.4 参数名集合 | 参数提案§5.2：expected / discovered / excluded / captured集合、来源/版本、差异、未知项和失败闭合规则 | expected集合尚未建立/实测；151仅旧观察数。发现集全部读到不代表预期集合完整 |
| R1.5 既有保护边界 | 参数提案§5.3、§6：双哈希、标量类型分离、安全排除、禁止反馈setParam、门控异常不触发额外fallback | 不能通过新探测、默认值或算法调整配平；完整依赖调用链仍待未来获准源码检查 |

## R2：逐字段零事件语义

| 子项 | 文档改动与定位 | 保留意见／未解决事项 |
|---|---|---|
| R2.1 字段适用表 | ZERO_EVENT_POLICY_V2_PROPOSAL§3.1：逐项列start/end、有限性、旧raw-zero硬门、band、B_start/dt来源、历史可比性和terminal处理；另列其他12种有限pool | 原四个raw-zero硬字段与新增描述字段分开；缺少历史全pool时标comparison_unavailable，不猜库存 |
| R2.2 连续比较不变 | 零事件提案§3–4：保留原连续量公式、atol/rtol和既有积分/清零逻辑；band仅用于比较 | 不向库存、bounds或模拟状态回写，不新增生物学阈值 |
| R2.3 正值边界 | 零事件提案§5：两边均为正值且连续量通过，即使near-zero标签不同也仅描述，不新增硬事件门 | 受限豁免仅处理原raw-zero不一致；超出既定连续容差仍STOP |
| R2.4 原值与版本 | 零事件提案§1、§5、§8：保留raw_exact_zero与所有原值；明确新政策为事后提出，不能冒充原预注册 | 候选政策名仍zero_event_policy_v2_proposed，必须同时绑定新文档SHA；trajectory仍1.8，不改旧输出 |

## R3：terminal仍strict STOP，但分开表达已观察结果

| 子项 | 文档改动与定位 | 保留意见／未解决事项 |
|---|---|---|
| R3.1 政策而非生物结论 | 零事件提案§6：terminal不继承上一积分区间band是保守政策选择，不是dt=0必然推出的生物学结论 | 保留strict raw-zero；不因infeasible推断细胞死亡 |
| R3.2 六维结果 | 零事件提案§6.1；WP2_REPLAY_V2_PLAN§7：新增terminal_raw_zero_only_mismatch及solver_termination_observed、solver_status_match、termination_time_match、continuous_state_match、raw_zero_match、overall_replay_gate | “已观察到相同solver终止”与“总体比较STOP”可以同时成立；示例为条件性说明，不是新运行结果 |
| R3.3 原参考末端 | 零事件提案§6.1给出原参考gzip的SHA及row2093/step85；独立打开确认5.3125 h非推进infeasible、glucose start/end均4.440892098500626e-16 | 仅查阅原始参考行，未续跑旧step84或生成新step85 |
| R3.4 不承诺消除STOP | 零事件提案§6.1、§10；重放计划§7：即使step84受限豁免，step85仍可能strict STOP | 不自动放宽terminal或Q9 depletion time，不承诺三个旧STOP全部消失 |

## R4：真实对象零求解预检独立待批

| 子项 | 文档改动与定位 | 保留意见／未解决事项 |
|---|---|---|
| R4.1 新权限 | DECISION_ITEMS§5–5.1；重放计划§5、§8：real_object_parameter_preflight单列awaiting_human | 创建真实对象可能触发license/资源，不能归入普通文件只读或synthetic tests |
| R4.2 候选范围 | 读取冻结输入/获准sidecar，按明确列出的原四项配置路径创建对象并读取参数、依赖身份，生成manifest；optimization_calls=0 | 不含optimize/pfba/smoke/状态推进；loader若有隐藏求解则不能按零求解范围执行。两mode对象、120秒/4GiB/16MiB均为待批上限 |
| R4.3 对象生命周期 | 重放计划§3.1、§5：预检为独立进程/对象，完成即释放；后续每一新worker/mode/native对象必须重新读回 | 预检manifest只作审批期望，不能冒充另一个对象的observed值 |
| R4.4 native观察许可 | 参数提案§7–8.1；决策§5：不安排额外native smoke；如需额外求解须另行明确授权。原有限运行中挂钩观察也须纳入具体batch许可 | 本轮未导入/查询solver，未执行任何smoke或native观察 |

## R5：比较层、WT分母与分批生命周期

| 子项 | 文档改动与定位 | 保留意见／未解决事项 |
|---|---|---|
| R5.1 3-only与8条 | 重放计划§3；决策§6：3-only只处理三处旧STOP证据问题，缺新的nonlimiting WT；补WT至少4条须明确另批范围 | 未选择运行范围；不隐含加条件，3条不能修复另外5条旧参数记录 |
| R5.2 新WT资格 | 重放计划§6.3：同一获准计划、mode、输入、设置和适用gate；分母有限且>0 | 旧WT不能冒充新WT；同hash不能替代完成验收 |
| R5.3 ratio元数据 | 重放计划§6.3：ratio_definition、WT/KO各自termination_time、denominator_run_id、condition_key及gate_status | 旧dynamic_growth_ratio是不同终止端点的倍增数比，不是共同时间窗、生长速率或最终biomass比；不可用时NA并写原因 |
| R5.4 四个比较层 | 重放计划§6.1：分别报告r2归档、v1已执行前缀、v2内部自洽、新WT归一化 | r2与v1冲突时保留双方值/身份，受影响比较STOP，不择易通过者；545行不由离线补齐 |
| R5.5 阶段可比性 | 重放计划§6.1；参数提案§5.3：v1外层pFBA返回解不代表每个native primary solve | 只有同阶段、同目标且有实际记录者可比；其余comparison_unavailable，不能把缺失写作差异为零 |
| R5.6 2+2+4生命周期 | 重放计划§3.1：每批新worker、每mode新模型/context/native对象；批内明确复用/回滚，批间释放；WT只传结果身份而非存活对象 | 不擅自reset/basis清除/增加求解；参数hash相同不证明内部状态完全相同 |

## R6：权限、资源预算及来源身份

| 子项 | 文档改动与定位 | 保留意见／未解决事项 |
|---|---|---|
| R6.1 权限名称 | DECISION_ITEMS§5为统一表；参数提案§9、重放计划§8和manifest用相同名称 | 六项指定权限另加必要的parameter_capture_code_implementation；本轮仅document_revision获准，其余awaiting_human，接受政策不隐含代码/测试/执行 |
| R6.2 非干扰声明 | 参数提案§5.3、§8.1；重放计划§9：仅限定数学输入、参数、原定调用及控制流 | elapsed/RSS/I/O/压缩开销另记，不承诺墙钟或调度完全不变；实测仍unverified |
| R6.3 预算/进行中调用 | 重放计划§9：30分钟/条件、2小时总额为含加载/记录/压缩的active wallclock，排除队列/人工等待；进程组4GiB、新增暂存及原始/压缩文件合计2GiB；外层限额、预提交in-flight记录、partial capture及退出收据 | 若环境无法提供获准限制与收据则不运行；只有轮询不能承诺严格上限。中断时保留incomplete_post_missing/未知求解结果，不伪造正常返回 |
| R6.4 禁止暗调参 | 重放计划§9；参数提案§5.3：不通过隐藏TimeLimit等setParam落实外层预算，不减记录、不改算法、不延时/重试 | 资源中断机制本身也须包含于明确batch/preflight许可，不借本轮文档授权执行 |
| R6.5 可访问且版本绑定的来源 | 下节提供PROJECT_STATE字节SHA、原文片段、实际目录/分支/HEAD观察来源；新manifest和SHA256SUMS绑定本轮文件，旧文件另列保护身份 | 本地核验不等于CoQ9聊天核验；来源文件陈述不自动等于科学重现 |

## 来源与身份片段（本地只读核对）

首次保护快照时间：`2026-09-05T21:22:14.589554+00:00`。定稿核对时间与逐文件前后SHA写入新 `document_manifest.json`。下列是本地可访问来源，不要求外部聊天凭空访问路径；关键原文同时摘录在此。

来源路径：`/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/PROJECT_STATE.md`；完整文件SHA-256：`4a2cbdea59483c438c9d06088815088eedf74580951c1469d162cd9f8d2136f8`。以下分别摘自该版本第19、20、24行，仅支持“此文件如此登记”，不是本轮重新执行这些科学工作：

```text
| 当前目录模型 | `576a284e…`；`artifacts.current_directory_model` | 当前工作树的文件，不自动成为暂定参考；本目录缺少 canonical 构建起点 `data/iyali26.xml`。 |
| 暂定参考的 canonical 文件 | `bc2aac8f…`；`artifacts.canonical_model`，持久文件在相邻 integration 工作树 | 作为历史筛查链的父模型登记；不因本轮登记新增“正式模型”批准。 |
| CoQ9 runtime 专题 | 文件使用 `bc2aac8f…`，另有动态配置、人工reserve与独立评价端点 | 仍为未校准敏感性方案，不能替换本静态参考结果。 |
```

实际工作区观察来自该目录的只读Git/文件检查，而不是以上文字的推断：

| 字段 | 实际观察 | 限制 |
|---|---|---|
| cwd | /Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem | 仅文档落地位置，不决定科研基线 |
| branch / HEAD | codex/r989-gpr-main-worktree / 994b09bf5f0e86b548094ff7bbb94296d37c4536 | 本轮不切换、不提交、不推送 |
| tracked diff | 起始为空；定稿再次核对见manifest | 原有未跟踪文件不清理；新文档不改既有tracked文件 |
| 目录model.xml SHA | 576a284ee86f0b96c802ea2e4445a862da49463022962471bfc841c469fcb5f2 | **不是**CoQ9冻结输入；仅核对保留 |
| data/iyali26.xml | 起始不存在；定稿核对见manifest | 不生成文件，不伪造缺失文件SHA |
| CoQ9冻结模型SHA | bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee | 来源为封存输入及已有身份记录，不被当前目录模型替代 |

外层handoff ZIP SHA为 `3d732a67cf3f90a9d26213cf2a160db491406f7492d855cf1ce06c42acef9001`；旧提案ZIP SHA为 `48259cccdf2e1d8acc384a77d1786c3dc659b69ac8d0a8d7cba8695fd72182ff`。v1六文件、两个输入及模型/状态入口的逐项起始/定稿身份由manifest列明。新文档SHA为本地实算，仅证明这些字节身份，不能证明候选科学或软件正确性。

## 保留结论、审查范围与停止点

- 旧结果仍为 **3条全部适用归档指标匹配、2条单条件匹配、3条原规则STOP**；545参考行未在旧WP2执行，未补造或离线改PASS。
- 14项旧参数原值未恢复，缺口涉及全部8条件；旧null、manifest、STOP、轨迹及封存SHA保留。
- 原native蛋白、guide QC、反应机制缺口继续按原等级；人工reserve不是实测催化CoQ池，不按recall调alpha/pool/cutoff/GPR。
- 本轮独立复核的实际范围和限制见DECISION_ITEMS§8；这是来源与文档符合性审查，不是代码/数学重跑或原生功能验证。全部候选仍等待人类审批。
- 四份proposal_v2及本回应、manifest、SHA256SUMS写入全新目录；版本系列沿用20260905，修订/封存时间按2026-09-06实记。若另生成ZIP，其实际字节SHA单独提供，不宣称CoQ9聊天已收到或核验。

**到此停止。** 不实施代码、不运行任何测试或求解、不导入/查询真实solver、不提交HPCC、不改模型/curation/dossier/旧STOP，不自动进入下一工作包或发送外部消息。
