# WP2 zero-event policy v2：仅比较层候选方案

日期：2026-09-05。状态：**proposal_only / awaiting_human；未实施、未测试、未重放**。

Material Passport：Origin Skill = experiment-agent / govern-agentic-research；Origin Mode = plan；Verification Status = UNVERIFIED（候选设计，不表示已验证规则）；Version Label = proposal_v1。

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
| 原始 exact-zero 不同，两侧均在同一 pair band 内，连续量通过 | raw_zero_disagreement_within_comparison_band | 描述性差异；未来重放可继续观察，但行为验收仍 pending |
| 原始 exact-zero 不同，任一侧超出 pair band | raw_zero_disagreement_outside_comparison_band | STOP，即使原连续量闸门通过 |
| 任意负库存 | invalid_negative_pool | STOP，不因很小而清零；IEEE -0 与 +0 数值等零但原符号保留 |
| 应为有限的状态是 NaN/Inf/缺失 | invalid_or_missing_state | STOP，不按零值处理 |
| 时间、推进状态、source-free status、reserve branch 或 zero-growth 分类不同 | dynamic_or_discrete_mismatch | STOP，不能被 near-zero 标签覆盖 |
| 内部 flux 不同但各动态闸门均通过 | possible_alternative_optimum_not_confirmed | 保留完整差异；不声称已证明 alternative optimum 或新机制 |

重要边界：

- 原始 glucose/finite-uracil 精确耗尽时间保持单独字段。只有同一网格上的 near-zero 分类一致、且原始差异完全属于上表第三行时，才可将原始精确耗尽时间差异保留为描述性差异；不能改写原时间或声称两个 raw event 完全相同。
- 新 near-zero 首次到达时间只是独立注释，不取代原始 depleted 字段。
- `q9_pool_depleted_h` 是 calls-level 原始指标；本提案**不放宽它的时间比较闸门**。若该值或其他受控 calls/终止事件不一致，仍 STOP。放宽此项需要新的明确提案与授权。
- pFBA 次级目标与 flux 的差异在动态闸门之后检查，质量守恒、bounds、coupling 残差超限仍停止；氧摄取本身不能支持“AOX 补偿”。

## 6. 无效区间、结构性 NaN 与既有 clamp

1. $B_{start}$ 必须有限且大于零，已推进区间 dt 必须有限且大于零。缺少任一侧起点数据时标 `not_comparable_missing_scale`，停止受影响比较，不能使用终点 biomass 代替。
2. nonoptimal terminal row 没有推进，`time_end == time`。该行不计算上述 epsilon、不除以零；near-zero 新标签为 `not_applicable_terminal`。仍核对 start/end 状态未变化、NaN flux 合法性、终止时间与 status。**该行适用库存的 raw exact-zero 不同仍严格 STOP**，不继承前一推进区间的 epsilon。原第三个 STOP 后的下一条参考记录正是终止行，因此新规则也可能在该行再停止；本提案不承诺三个STOP都能通过。若希望终止行也可继承已记录band，必须单列新规则及审批，不能由实施者猜测。
3. `po1f_nonlimiting` 的 uracil 库存列按原设计是结构性 NaN，标 `not_applicable_nonlimiting`；这不是 missing input，也不是零库存。finite_batch 的 uracil NaN 不可按此豁免。
4. frozen runner 已有 `q9_tolerance = max(1e-12, initial_q9_pool * 1e-9)`，并在更新 Q9 reserve 后执行原有小值清零。本轮及本候选比较设计均不新增、不删除、不修改此机制。
5. 原始 `q9_source_total_mmol_L = initial_pool - final_pool`、有限库存 `max(0, start - uptake*B_start*dt)`、biomass Euler 更新保持原样。不能借“比较修复”重定义累计量。

来源：F `scripts/gem_annotate/quinone_dfba_essentiality.py:429–464, 570–655, 692–707`。特别是 576、642–645 行的既有 Q9 clamp，不得误称本方案新引入。

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
| 同一零值标签但 branch/status/zero-growth 不同 | STOP；标签不能掩盖离散差异 |
| 初始 band 差异后下一步状态超限 | 保存首个动态分歧并 STOP，不继续调参 |
| q9_pool_depleted_h 或终止时间不同 | 原有 calls/时间闸门仍 STOP |
| 新规则应用于 legacy 文件/未知 schema | 拒绝原地写入和静默 schema 混用 |
| 未执行的 tail | 不创建伪造历史记录或完整性 PASS |

这些只是待批准的测试案例，不是已运行或通过的测试结果。

## 10. 待决事项

- 是否批准上述比较层 policy、常数与明确保留的 strict gates？
- 是否另行授权最小 sidecar 代码与 synthetic tests？
- 测试与参数捕获修复通过之后，是否再批准三条件诊断或八条件新证据链？参见 `WP2_REPLAY_V2_PLAN.md`。
- 没有任何 near-zero 标签可以回答未执行区间、缺失历史 14 参数、真实 CoQ pool 或原生蛋白功能问题。

## 来源路径约定

H = `coq9_wp2_review_and_next_step_handoff_20260905.zip::coq9_wp2_review_and_next_step_handoff_20260905/`。

A = H `coq9_wp2_independent_audit/`；具体审计文件名和 SHA 见 `document_manifest.json`。

B = H `original_archives/coq9_wp2_qualified_baseline_20260905.zip::coq9_wp2_qualified_baseline_20260905/`。

F = B `execution_provenance/coq9_wp2_qualified_20260905T042505Z_inputs.tar.gz` 内的冻结 compute 目录。以上 locator 指归档内文本，只作为证据读取；本轮没有执行其中脚本。
