# WP2 zero-event policy v2：仅比较层候选方案

起始日期：2026-09-05；文档修订日期：2026-09-06。状态：**proposal_only / awaiting_human；未实施、未测试、未重放**。

Material Passport：Origin Skill = experiment-agent / govern-agentic-research；Origin Mode = plan；Verification Status = UNVERIFIED（候选设计，不表示已验证规则）；Version Label = proposal_v2；document_schema_version = 2。

本版仅按 `TO_INACTIVE_NEXT_STEPS_zh.txt` 的 R2/R3/R6 定点修订；v1 原件、旧 raw values/STOP/manifest/封存 SHA 不变。CoQ9 聊天的 v1 总体审查为 partial，仅基于全文与其可访问的来源；它没有核验新 ZIP 或文件字节。该审查范围不能替换为本地来源子代理的核查范围。

## 1. 要修复的是比较语义，不是模拟状态

建议在新的比较 sidecar 中并列保留：

1. 原始数值及原始 `value == 0` 判定；
2. 沿用已有容差的连续量比较；
3. 新增、明确标为事后提出的 near-zero 比较标签；
4. 未来获准重放时才可能观测的后续状态、分支和终止结果。

这四层不能合并。旧 WP2 的三个 STOP 永久保留；新标签既不能把旧 STOP 改成 PASS，也不能补出未运行的 545 行。方案不授权修改 frozen runner、Euler 方程、库存、exchange bounds、Q9 reserve 策略或求解流程。

本方案名为 `zero_event_policy_v2_proposed`。它是看过旧差异后提出的新版本，不是原 WP2 的预注册规则。即使离线标签一致，后续行为等价仍未建立。

## 2. 已观察到的三处差异

以下是交接包中的已归档结果，不是本轮重新计算：

| 条件键（alpha=1e-4，pool=1，dt=0.0625 h） | step；区间 h | 字段 | WP2 值 | 冻结参考值 | 旧判断 |
|---|---|---|---:|---:|---|
| finite_batch / WT | 56；[3.5, 3.5625) | uracil_end_mmol_L | 0 | 3.469446951953614e-18 | STOP |
| finite_batch / YALI1A14736g | 166；[10.375, 10.4375) | uracil_end_mmol_L | 1.9233746539892849e-16 | 0 | STOP |
| po1f_nonlimiting / YALI1A21711g | 84；[5.25, 5.3125) | glucose_end_mmol_L | 0 | 4.440892098500626e-16 | STOP |

三个绝对差都低于原有 finite-pool 连续量 atol=1e-10 mmol/L，但原始 `== 0` 布尔值不同。这支持“旧比较器的连续量闸门与精确零事件闸门不同”；仅部分支持数值边界解释，不证明后续状态或生物学等价。

基因身份：YALI1A14736g — 本材料未核实公认名称 — 模型 R305 关联组分候选（仅模型 GPR 赋值；个体蛋白功能及原生身份 unresolved）；YALI1A21711g — 本材料未核实公认名称 — 模型 R2062 关联组分候选（仅模型 GPR 赋值；个体蛋白功能及原生身份 unresolved）。这里仅以条件 ID 定位旧停止点，不新增功能结论。

来源：A `exact_zero_points.tsv` 第2–4行、`condition_verdicts.tsv`、`evidence_status.tsv`；B `analysis_code/coq9_wp2_replay.py:107–130`（连续比较与 exact-zero 比较分开）。A/B 路径别名见文末及 `document_manifest.json`。

本节事实等级：三处原值/旧STOP与原continuous规则为 **supported**（封存文本/表）；“raw-zero相同”为 **contradicted**；“数值边界解释”仅 **partial**；尾段等价、生物学回归和候选规则实际有效性均 **unverified**。本轮没有重放这些区间。

## 3. 连续量容差不变

沿用 B `PREDECLARED_PROTOCOL.json#/comparison_tolerances` 及 replay 源码 31–41、99–104 行：

`abs(new - reference) <= atol + rtol * abs(reference)`。

不能悄悄改为对称的 max-scaled 相对容差；两个合法缺省 NaN 的历史处理也不能扩展到本应有限的状态。

| 量 | atol | rtol | 单位/规则 |
|---|---:|---:|---|
| growth、source-free growth | 1e-8 | 1e-8 | h^-1 |
| biomass | 1e-10 | 1e-9 | gDW/L |
| finite medium pool | 1e-10 | 1e-9 | mmol/L |
| Q9 artificial reserve | 1e-12 | 1e-9 | mmol/L |
| dynamic doublings / ratio | 1e-8 | 1e-8 | 无量纲 |
| reaction flux | 1e-8 | 1e-8 | mmol/gDW/h；biomass flux 另按其单位 |
| pFBA secondary objective | 1e-6 | 1e-8 | 原求解目标的数值单位 |
| Sv / bound residual | 1e-7 | 0 | 沿用旧审计定义 |
| coupling residual | 1e-9 | 0 | 沿用旧审计定义 |
| interval time | 0 | 0 | h；固定时间网格 |

status、step_index、interval_advanced、source-free status、reserve branch 继续精确比较。`growth <= 1e-9` 的既有 optimal zero-growth 分类与 `source-free growth > 1e-9` 的分支阈值均不改变。

### 3.1 R2 逐字段适用表

旧 `compare_row` 中 hard raw-zero 集合**恰为** `q9_pool_mmol_L`、`q9_pool_end_mmol_L`、`glucose_end_mmol_L`、`uracil_end_mmol_L`（B replay:125–127，**supported**）。其余字段不能因名称相近被增补进这个集合。旧 nonlimiting uracil 的两个 NaN 在旧数值`==0`中均为false；v2并列保留这个raw结果，同时标注结构性不适用，不能把NaN写成零。

下表“band=豁免候选”仅指有效推进optimal行内、连续量已通过时，对上述旧硬集合的raw-zero差异按§5判断；“描述”只产生标签，不产生新硬门。所有候选表规则状态为 **unverified / proposed_not_implemented**。

尺度 S = 两侧同一条件/step实际 `biomass_gDW_L`（B_start）及固定dt=0.0625 h，且各侧 `time_end_h-time_h=dt`；绝不使用B_end。非推进terminal行一律不计算或继承band。

| 字段 | start/end/速率 | 必须有限的条件 | 旧raw-zero硬门 | v2 band与尺度 | 历史可比性 | terminal处理 |
|---|---|---|---|---|---|---|
| time_h | start，h | 每条有效记录 | 无；时间另按0容差 | 不适用 | r2与v1已存行 | 等于time_end；按原时间门 |
| time_end_h | end，h | 每条有效记录 | 无；时间另按0容差 | 不适用 | 同上 | 不推进；必须等于time_h |
| biomass_gDW_L | start，gDW/L | 所有有效状态，且>0才可构造尺度S | 无 | 不适用；只作为S输入 | r2/v1已存行 | 与本行end相同；连续量比较 |
| biomass_end_gDW_L | end，gDW/L | 所有有效状态 | 无 | 不适用；不能代替B_start | r2/v1已存行 | 与本行start相同 |
| q9_pool_mmol_L | start，mmol/L | 两模式所有有效状态；非负 | 有 | 豁免候选，Q9 cap+S | r2/v1已存行 | strict raw-zero；无band |
| q9_pool_end_mmol_L | end，mmol/L | 两模式所有有效状态；非负 | 有 | 豁免候选，Q9 cap+S | r2/v1已存行 | 等于start；strict raw-zero |
| glucose_mmol_L | start，mmol/L | 两模式所有有效状态；非负 | 无 | 描述，finite cap+S | r2/v1已存行 | 连续量比较；不新增raw-zero硬门 |
| glucose_end_mmol_L | end，mmol/L | 两模式所有有效状态；非负 | 有 | 豁免候选，finite cap+S | r2/v1已存行 | 等于start；strict raw-zero |
| uracil_mmol_L | start，mmol/L | finite_batch必须有限非负；nonlimiting为原结构NaN | 无 | finite模式仅描述+S；nonlimiting不适用 | r2/v1已存行 | finite连续量比较；nonlimiting保留NaN |
| uracil_end_mmol_L | end，mmol/L | finite_batch必须有限非负；nonlimiting为原结构NaN | 有；NaN的旧raw bool原样保留 | finite模式豁免候选+S；nonlimiting不适用 | r2/v1已存行 | finite等于start并strict raw-zero；nonlimiting仍NaN |
| state中其余有限pool的start/end | 分别start/end，mmol/L | 对模式实际有限集合每项均有限非负 | 无compare_row raw-zero门 | 如记录标签，仅描述+S，绝不新增硬门 | r2无该12种逐步完整库存；v1只在其已保存前缀有全库存 | end=start；检查现有内部连续性，不猜历史值 |
| growth_h-1 / biomass_flux_h-1 / objective_value | 起点求解的区间速率/目标 | optimal行须有效有限；nonoptimal原无效字段为NaN | 无库存零门；仅growth沿原<=1e-9分类 | 不适用 | 各来源实际保存字段分别比；不混淆目标阶段 | 仍NaN，不假造增长 |
| source_free_growth_h-1 | source-free求解返回值 | source-free optimal时须有限；非optimal按原缺省NaN | 无库存零门；分支阈值保持原>1e-9 | 不适用 | 只比已保存对应阶段 | 依原status判断适用，不将失败归零 |
| uptake/source/dilution及其余flux | 起点状态求出的区间速率 | 该次optimal返回解可用时有限；无解为原NaN | 无 | 不适用；沿flux容差及残差门 | r2面板与v1返回解分别定位；无同阶段记录则不可比 | 不推进、不用NaN重建消耗 |
| q9_pool_depleted_h | calls事件时间；不是库存列 | 发生耗尽或初始为0时有限；未耗尽按原缺失表示 | 无库存零门，但calls事件时间独立硬门 | 不适用；不能用near-zero时间替代 | 仅来源确有calls端点时可比 | 原raw depletion事件/时间门不放宽 |

“其余有限pool”集合为R1003、R1202、R1204、R1215、R1217、R1220、R1222、R1223、R1231、R1232、R1233、R1234。R1070/R1354分别已由glucose/uracil行覆盖。氧摄取有flux不代表另造有限氧pool。完整初始化与模式差异见重放计划§4。

未保存的历史字段标 `comparison_unavailable`，不是零、相等或PASS；新状态中按schema必须存在的字段缺失仍是记录失败。其他pool的Euler/连续性闸门保持原定义，不能把新增描述标签当作新事件门。

## 4. 一个可审批的 near-zero 候选公式

目的：把“库存差异”与该区间可对应的“摄取通量差异”同时限制在事先声明的比较尺度内。它不是实验检出限、真实耗尽阈值或求解器容差的物理换算。

对一个有效、已推进的 optimal 区间，定义：

$$
\epsilon_C=\min(A_C,\epsilon_v B_{start}\Delta t).
$$

候选常数：

- finite medium pool：$A_C=10^{-10}$ mmol/L；
- artificial Q9 reserve：$A_C=10^{-12}$ mmol/L；
- $\epsilon_v=10^{-8}$ mmol/gDW/h；
- $B_{start}$ 单位 gDW/L；$\Delta t$ 单位 h，使用该区间真实起点状态及已固定 dt。

两个待比较记录采用两侧阈值的较小者：

$$
\epsilon_{pair}=\min(\epsilon_C^{new},\epsilon_C^{reference}).
$$

若 $0\le C\le\epsilon_{pair}$，可标 `near_zero_for_comparison`；高于阈值则标 `positive_above_comparison_band`。必须记录原值、原始 exact-zero、两侧 epsilon、B_start、dt、单位及 policy 身份。不得将标签再写回数值。

选择依据是**已有的库存与通量绝对比较预算**，不是把 epsilon 设成略高于三个已知误差。该候选尚未通过合成测试，也没有用旧三点计算“最优 epsilon”。如果后续无法通过，就报告失败；不得为通过旧 STOP 增大 epsilon。

量纲上，$\epsilon_C/(B_{start}\Delta t)\le\epsilon_v$。这只说明由极小库存量换算的区间库存限制通量尺度受到约束，不能保证优化解、活跃约束或未来事件必然相同。因此后续行为必须另验。

## 5. 标签与停止条件

| 情况 | 新比较记录 | 候选处理；仍须另行授权实施/重放 |
|---|---|---|
| 连续量不通过 | continuous_mismatch | STOP，保存最早分歧 |
| 原始 exact-zero 相同、连续量通过 | raw_event_agrees | 继续其余闸门，不能据此单独 PASS |
| 两侧均为正、连续量通过，但near-zero标签不同 | positive_positive_band_label_difference | **仅描述**；raw-zero同为false，不新增硬事件门；其他原闸门照常 |
| 原始 exact-zero 不同，两侧均在同一 pair band 内，连续量通过 | raw_zero_disagreement_within_comparison_band | 描述性差异；未来重放可继续观察，但行为验收仍 pending |
| 原始 exact-zero 不同，任一侧超出 pair band | raw_zero_disagreement_outside_comparison_band | STOP，即使原连续量闸门通过 |
| 任意负库存 | invalid_negative_pool | STOP，不因很小而清零；IEEE -0 与 +0 数值等零但原符号保留 |
| 应为有限的状态是 NaN/Inf/缺失 | invalid_or_missing_state | STOP，不按零值处理 |
| 时间、推进状态、source-free status、reserve branch 或 zero-growth 分类不同 | dynamic_or_discrete_mismatch | STOP，不能被 near-zero 标签覆盖 |
| 内部 flux 不同但各动态闸门均通过 | possible_alternative_optimum_not_confirmed | 保留完整差异；不声称已证明 alternative optimum 或新机制 |

本表中的raw-zero豁免/STOP仅作用于§3.1列明的旧硬集合。没有旧硬门的start或其他pool字段，标签与raw-zero差异只作描述，仍受原连续量/记录有效性闸门约束。采用R2建议：band只用于描述与已有raw-zero不一致的受限豁免，不增加“两个正值跨band边界即STOP”的新政策。

重要边界：

- 原始 glucose/finite-uracil 精确耗尽时间保持单独字段。只有同一网格上的 near-zero 分类一致、且原始差异完全属于上表的 `raw_zero_disagreement_within_comparison_band` 时，才可将原始精确耗尽时间差异保留为描述性差异；不能改写原时间或声称两个 raw event 完全相同。
- 新 near-zero 首次到达时间只是独立注释，不取代原始 depleted 字段。
- 两正值因band标签不同而产生的near-zero首次时间差，也仅描述，不绕过前述规则另设全程near-zero时间硬门。
- `q9_pool_depleted_h` 是 calls-level 原始指标；本提案**不放宽它的时间比较闸门**。若该值或其他受控 calls/终止事件不一致，仍 STOP。放宽此项需要新的明确提案与授权。
- pFBA 次级目标与 flux 的差异在动态闸门之后检查，质量守恒、bounds、coupling 残差超限仍停止；氧摄取本身不能支持“AOX 补偿”。

## 6. 无效区间、结构性 NaN 与既有 clamp

1. $B_{start}$ 必须有限且大于零，已推进区间 dt 必须有限且大于零。新schema必需状态缺失/非法则STOP；某历史来源原本未保存起点数据时标 `comparison_unavailable`。无法构造S就不能使用band豁免原硬raw-zero差异，也不能使用B_end代替。对仅描述字段，缺少历史尺度不自动新增硬事件门。
2. nonoptimal terminal row 没有推进，`time_end == time`。该行不计算上述 epsilon、不除以零；near-zero 新标签为 `not_applicable_terminal`。仍核对 start/end 状态未变化、NaN flux 合法性、终止时间与 status。**仅§3.1旧硬集合的适用raw exact-zero差异仍严格STOP**，不继承前一推进区间band。这是保守的比较政策选择，不是dt=0必然推出的数学/生物学结论。原第三个STOP后下一参考行就是terminal，故新规则仍可能再STOP；不承诺三个STOP都通过。放宽terminal或Q9 depletion时间门必须另外提出并接受政策，实施者不可自行决定。
3. `po1f_nonlimiting` 的 uracil 库存列按原设计是结构性 NaN，标 `not_applicable_nonlimiting`；这不是 missing input，也不是零库存。finite_batch 的 uracil NaN 不可按此豁免。
4. frozen runner 已有 `q9_tolerance = max(1e-12, initial_q9_pool * 1e-9)`，并在更新 Q9 reserve 后执行原有小值清零。本轮及本候选比较设计均不新增、不删除、不修改此机制。
5. 原始 `q9_source_total_mmol_L = initial_pool - final_pool`、有限库存 `max(0, start - uptake*B_start*dt)`、biomass Euler 更新保持原样。不能借“比较修复”重定义累计量。

来源：F `scripts/gem_annotate/quinone_dfba_essentiality.py:429–464, 570–655, 692–707`。特别是 576、642–645 行的既有 Q9 clamp，不得误称本方案新引入。

### 6.1 R3 terminal结果的多维表达

增加候选结果类型 `terminal_raw_zero_only_mismatch`，只在新求解器的非optimal终止已经实际返回、其status/time/连续状态及其他适用门相符、而差异仅剩旧硬集合raw-zero时使用；否则应分别保留实际行为/输入/记录失败，不能套用此类型。

| 必需字段 | 含义 | 假设“新step85同样infeasible、连续量相符但glucose end为0”的表达（不是新运行结果） |
|---|---|---|
| solver_termination_observed | 原求解流程是否已返回终止证据；非比较器是否通过 | true |
| solver_status_match | 两来源真实status是否相同 | true |
| termination_time_match | 两来源终止时间是否通过原0容差门 | true |
| continuous_state_match | 适用连续状态是否均通过原容差 | true |
| raw_zero_match | 旧硬集合的全部适用raw-zero是否相同 | false |
| overall_replay_gate | 包含strict terminal政策的综合比较判断 | STOP |

每项保存左右原值、reference ID/SHA、condition/step、适用性及原因。没有新的求解返回时不得把`solver_termination_observed`填true；历史/新字段缺失用NA+原因，不伪造false/零。

直接读取的历史参考（**supported**）：B `reference/baseline_reference_intervals.tsv.gz`（SHA `959387b366ae95b1f292e109fee6ca4f75d532612878828c78b6b013ec78326e`），解压TSV含表头第2093行：`po1f_nonlimiting / YALI1A21711g / step_index=85`，`time_h=time_end_h=5.3125`、`interval_advanced=False`、`status=infeasible`，glucose start/end均为 `4.440892098500626e-16`。第2092行是step84的optimal区间。这个终止行存在于r2参考；旧WP2在step84 STOP，未执行该step85。

若将来新运行由step84的零glucose进入相同终止，仍可能呈现上表类型：可以同时说“已观察到相同求解器终止”和“整体比较STOP”，不必把二者混成“没有观察终止”，更不能称死亡或获准放宽。新行为是否如此仍 **unverified**。

## 7. 后续行为验收的范围

未来若用户批准新重放，新记录须同时比较：每步目标、source-free status/μ、全部实际可比库存、reserve branch、zero-growth、耗尽/终止、calls，然后比较 pFBA 次级目标与完整 flux/residual/氧-ATP 账本。

能比较的历史范围必须分别标注：

- 旧 WP2 已运行区间：可比其归档的全库存、全通量、LP/状态 sidecar；
- r2/r3 参考仍有、但 WP2 未执行的 545 行：仅比该参考实际保存的字段；新 sidecar 可以自洽核查，但不存在的旧全库存/全通量不能猜补；
- 旧参考也不存在的时段：不外推、不声称历史等价。

新轨迹即使通过新规则，也只能产生“在声明规则、可用历史字段范围内通过”的新结论，不能回填旧 STOP、恢复未知历史参数或宣称真实细胞死亡。

## 8. 最小实施范围与版本隔离（仅提案）

复用现有 sidecar 的连续比较 helper；只新增上述标签与受控事件比较逻辑、明确 schema/policy ID 和后续单元测试。不新建通用数值容差框架、不修改 frozen compute。

实现文件候选为独立分析分支中的 `coq9_wp2_replay.py` 比较/报告段及一个紧邻的 synthetic test 文件；实际受控路径需在实施授权时列明。旧 trajectory schema 仍为 1.8；新 sidecar/comparator version 单独版本化。禁止静默转换、覆盖或混合旧记录。

`comparison_policy_acceptance`只接受规则，**不**授权`comparison_code_implementation`；后者不自动包含`synthetic_test_execution`、`real_object_parameter_preflight`或任何求解。完整权限表见DECISION_ITEMS§5。本文件内容版本为proposal_v2，候选政策仍以`zero_event_policy_v2_proposed`加本文件新SHA绑定；不能仅凭相同政策名混用v1/v2文档字节。

每个新运行 manifest 在开始前封存 policy 文本 SHA、常数、公式、单位、continuous comparator identity、expected mode-specific NaN、运行输入和设置身份。规则变化必须新版本、新 run ID、重新审批，不按结果自动调整。

## 9. Synthetic test 计划——本轮不运行

| 合成场景 | 预期行为 |
|---|---|
| +0、-0；小正值在 band 内/边界/外 | 原始值/符号保留；边界包含；外侧不能获 band 豁免 |
| 很小负数；有限模式 NaN/Inf/缺失 | 明确拒绝；不清零、不当 null 等价 |
| biomass/dt 很小与很大 | 验证量纲、min 两个上限与 pair 较小值；不改状态 |
| B_start/dt 为零、负或缺失 | 拒绝非法推进区间；不替换为默认值 |
| 非推进 terminal row | 不计算 epsilon、不继承前一区间band；状态不推进；raw exact-zero 不同严格STOP |
| nonlimiting uracil NaN / finite uracil NaN | 前者明确不适用，后者失败 |
| 连续量通过但 raw-zero 不同 | 仅 band 内可标描述性差异；不生成历史 PASS |
| 两边均正、连续量通过，但仅一边被标near-zero | 仅记positive_positive_band_label_difference；不新增STOP或全程near-zero时间硬门 |
| terminal只有旧硬集合raw-zero不符，其余终止维度相符 | terminal_raw_zero_only_mismatch；solver终止已观察与overall STOP同时保留 |
| glucose_start/其他pool没有旧raw-zero硬门 | 不以新增标签将其机械升级为硬事件门 |
| 同一零值标签但 branch/status/zero-growth 不同 | STOP；标签不能掩盖离散差异 |
| 初始 band 差异后下一步状态超限 | 保存首个动态分歧并 STOP，不继续调参 |
| q9_pool_depleted_h 或终止时间不同 | 原有 calls/时间闸门仍 STOP |
| 新规则应用于 legacy 文件/未知 schema | 拒绝原地写入和静默 schema 混用 |
| 未执行的 tail | 不创建伪造历史记录或完整性 PASS |

这些只是待批准的测试案例，不是已运行或通过的测试结果。

## 10. 待决事项

- 是否接受上述比较层policy、常数及strict gates（comparison_policy_acceptance）？
- 是否另行授权comparison_code_implementation和synthetic_test_execution？二者彼此分开。
- 测试与参数捕获修复通过之后，是否再批准三条件诊断或八条件新证据链？参见 `WP2_REPLAY_V2_PLAN.md`。
- 没有任何 near-zero 标签可以回答未执行区间、缺失历史 14 参数、真实 CoQ pool 或原生蛋白功能问题。

## 来源路径约定

H = `coq9_wp2_review_and_next_step_handoff_20260905.zip::coq9_wp2_review_and_next_step_handoff_20260905/`。

A = H `coq9_wp2_independent_audit/`；具体审计文件名和 SHA 见 `document_manifest.json`。

B = H `original_archives/coq9_wp2_qualified_baseline_20260905.zip::coq9_wp2_qualified_baseline_20260905/`。

F = B `execution_provenance/coq9_wp2_qualified_20260905T042505Z_inputs.tar.gz` 内的冻结 compute 目录。以上 locator 指归档内文本，只作为证据读取；本轮没有执行其中脚本。
