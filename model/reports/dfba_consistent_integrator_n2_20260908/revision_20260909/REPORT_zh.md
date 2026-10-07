# iYali26 dFBA 一致库存积分与耗尽事件定位：隔离原型修订

## 结论

本轮交付了一个**默认不接入生产入口**的纯数值原型：在一个子区间内冻结生长率和各库存的有符号单位生物量净速率，以同一个累计生物量暴露量同时更新生物量与全部库存，并在步内定位最早耗尽事件。三组预注册确定性例子全部通过；真实 GEM、FBA、pFBA、FVA、QP 调用均为 **0**。

这只证明内核在声明的常速率、恒体积数学范围内正确，并显示变速率问题的误差随两次步长减半而下降。它没有证明真实 GEM 已收敛，没有检测 pFBA 动态投影多解，也没有授权替换原 Euler 路径。

## 身份、范围与历史关系

- 实际工作区：`/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem`
- 分支/HEAD：`codex/r989-gpr-main-worktree` / `994b09bf5f0e86b548094ff7bbb94296d37c4536`
- 工作树在本轮开始前已含其他修改和未跟踪文件；本轮未清理、覆盖或提交它们。
- 既有 N2 v1 原型及其校验和保持不变。本修订放在其 `revision_20260909/` 子目录，避免把历史状态改写成当前状态。
- 既有真实 GEM event-resolve gate 也保持不变：执行完整性 PASS，但因候选在 `5.048080683449264 h` 的葡萄糖耗尽后立即 infeasible、未达到 `T*=5.109375 h`，预注册科学闸门 FAIL。`infeasible` 是求解状态，不是细胞死亡证据。
- 当前工作区 `model.xml` 的 SHA-256 是 `576a284e...`；历史 gate 实际执行的是另一获授权工作树中 SHA-256 为 `bc2aac8f...` 的模型。本轮未加载任何模型，因此两者不能互相冒充。

## 针对性代码审计

| 分类 | 结论 |
|---|---|
| 已确认代码事实 | CoQ9 runner 从 `biomass_C` 的实际 reaction flux 读取 `mu`；保存的 `objective_value` 被显式写成 growth，而不是 COBRApy pFBA 第二阶段总通量目标。 |
| 已确认代码事实 | `B` 为 gDW/L；培养基及人工 Q9 储备为 mmol/L；`mu` 为 h⁻¹；单位生物量通量为 mmol/gDW/h。已检查的轨迹是恒体积 batch，没有补料、蒸发或反应器体积稀释。 |
| 已确认代码事实 | 旧 CoQ9 更新为左端 Euler：`B += mu*B*dt`；Q9 储备只扣 `source*B*dt`；培养基库存扣 uptake×`B*dt`。运行时细胞内 Q9 平衡含 source 和 `COQ9_DILUTION=alpha*biomass flux`，外部人工储备只扣 source，因此没有把 dilution 再扣一次。 |
| 已确认代码事实 | 旧 CoQ9 有 `C/(B*dt)` 培养基 cap 和 `Q/(B*dt)` source cap；通用 FN runner 有 `C/(B*step*|stoich|)` cap。因此减小 `dt` 同时改变可行域和积分误差。这里的 `dt/step` 是计划求解步长；旧路径没有独立的“实际接受事件子步长”。 |
| 已确认代码事实 | 旧 CoQ9 对 Q9 使用截断和 near-zero snap，对培养基使用 `max(0, raw)`；事件只落在网格端点。通用 FN runner 允许小至 `-1e-8` 后截零。原始负值和修正量并非所有历史轨迹都有保存。 |
| 已确认代码事实 | CoQ9 `finite_batch` 有 14 个培养基有限库存加 Q9 储备；`po1f_nonlimiting` 去掉 R1354 uracil，故为 13+Q9。14 个培养基反应是 R1070、R1003、R1202、R1204、R1215、R1217、R1220、R1222、R1223、R1231、R1232、R1233、R1234、R1354。历史 CoQ9 trajectory 只保存 glucose、uracil、Q9 的状态，遗漏其余已更新营养库存的逐步状态。当前通用 FN 配置另含 R1189 bioavailable iron，成功结果只保存终点浓度。 |
| 已确认代码事实 | CoQ9 non-optimal 会写一条零时长 terminal record 并停止；optimal 零增长可继续到 horizon。通用 FN 的普通 FBA 只是 infeasible 诊断，不用于继续轨迹。timeout、infeasible、零增长、耗尽和正常 horizon 不是同一事件。 |
| 已确认代码事实 | 已检查的 COBRApy 0.30 pFBA 默认 `fraction_of_optimum=1.0`：先固定原 biomass 最优目标，再最小化每个反应 forward/reverse split variable 的无权重总和。旧 CoQ9 在 `OptimizationError` 后可回退普通 FBA；历史 gate 的 Gurobi readback 为 FeasibilityTol `1e-9`、OptimalityTol `1e-7`、Threads 1、Seed 0、Method 1、Presolve 0。 |
| 合理假设 | 本原型只对固定体积、`mu>=0`、每个区间内冻结 `mu` 与有符号净速率、且所有状态变化均与 `B` 成比例的系统适用。例 3 的 mock 上层在每个接受步起点重新给出速率。 |
| 仍未知 | 在真实状态下，pFBA 最优集合投影到 `mu`、Q9 净消耗和全部有限库存净交换后是否唯一；新内核与旧 `dt`-依赖库存 cap 如何在不改变政策的情况下接线；真实 GEM 的轨迹、事件和收支是否在 `h,h/2,h/4` 下收敛。 |

主要静态证据位置：历史 CoQ9 库存定义与 cap/求解/更新见 `artifacts/r608_engineering_20260907/code/scripts/gem_annotate/quinone_dfba_essentiality.py:115-132,430-465,558-680`；当前通用 FN 路径见 `scripts/dfba_new_fn_essentiality.py:451-458,591-676`；COBRApy pFBA 层级见 `.venv/lib/python3.13/site-packages/cobra/flux_analysis/parsimonious.py:44-49,90-97,128-140`。这些是静态核查，不是本次复现求解。

## 实现及数学契约

输入为 `B>0`、`mu>=0`、计划步长 `h>0`、每个库存 `C_i>=0`，以及调用者预先合并 source、sink、net exchange、dilution 后得到的有符号净速率 `s_i`。正 `s_i` 增加库存，负 `s_i` 消耗库存。原型不自行猜测化学方向。

在冻结区间中：

```text
dB/dt = mu B
dC_i/dt = s_i B
I(h) = B h                    , mu = 0
I(h) = B expm1(mu h) / mu     , mu > 0
B_next = B + mu I(h) = B exp(mu h)
C_i,next = C_i + s_i I(h)
```

对每个 `C_i>0, s_i<0`，由 `log1p` 计算耗尽时间，选择全部库存中的最早事件。先以同一个 `I(h_accept)` 更新所有状态；返回 `requires_resolve=true` 后，调用者必须在同一物理时刻更新可用性并重新求速率，不能把事件前通量延用到事件后。

原型明确区分：

- `proposed_dt_h`：计划步长；
- `accepted_dt_h`：受事件限制后实际推进的子步长；
- `remaining_dt_h`：必须由上层重求速率后处理的剩余时间。

起始空库存若仍给出负净速率会直接拒绝，避免零时长循环；负生长不支持并直接报错；明显负库存或事件选择后的实质透支报错。只对**正起始库存且确有消耗**的耗尽残差在 `1e-15 mmol/L + 1e-14×scale` 内 snap 到零；静止或生成中的微小库存不会被容差清零。原型同时返回 raw 值、signed correction、raw/corrected balance residual。函数先复制输入，失败不返回部分状态，也不修改调用者字典。

未支持：变化体积、补料、死亡、降解、与 `B` 不成比例的项和负增长。需要这些项时必须先重写收支，不能套用本公式。例 3 的纯 Python 循环就是本轮 mock 适配器；唯一生产接入缺口是把“速率求解/库存可用性/计划求解步长/接受事件子步长”显式分开，同时证明保留的 cap 政策与更新一致。

## 三组确定性例子

复现命令：

```bash
/usr/bin/python3 -B artifacts/dfba_consistent_integrator_n2_20260908/revision_20260909/test_consistent_event_integrator.py --output artifacts/dfba_consistent_integrator_n2_20260908/revision_20260909/metrics.json
```

| 例子 | 步长/事件 | 关键结果 | 判定 |
|---|---|---|---|
| 恒定速率 | `h=0.5 h`，无事件 | 新内核 biomass 解析误差 `0`；旧 Euler 在该例低估 `0.02568330979220379 gDW/L`；`mu=0` 暴露量误差 `0`；最大收支残差 `4.44e-16 mmol/L`；输入未改变；微小静止/生成库存未被清零 | PASS |
| Q9 严格耦合混合事件 | 事件在 `0.4731134450923251 h`，步内 early | 局部事件时间误差 `0`；`R+alpha(B-B0)-R0=-2.12e-22 mmol/L`；biomass cap overage `0`；另一库存可先耗尽；ulp 级同时事件均被识别；最大修正 `2.12e-22 mmol/L`；无透支/零时长循环 | PASS |
| `mu(t)=a+bt` 解析参考 | `h=0.2,0.1,0.05 h` | 共同时间 biomass 误差依次 `0.02001735,0.01004871,0.00503440 gDW/L`；整体事件时间误差依次 `0.03297206,0.01675584,0.00856007 h`；每个局部冻结段的事件库存残差为 `0` | PASS |

因此，常速率例中的一致解析推进把该步的 Euler biomass 误差从 `0.02568` 降到机器精度；但变速率例仍表现出近似随步长线性下降的全局误差。本轮只预注册“连续下降”，没有估计或声称高于一阶的收敛阶，也不能推广为所有 GEM 轨迹都低估。

## pFBA 动态投影歧义：后续诊断规格（未执行）

状态 `x` 下保持现有层级：

```text
g*(x) = max c^T v,  v in F(x)
K(x) = {v in F(x): c^T v >= 1.0 g*(x),
        sum_j(v_j^+ + v_j^-) <= s*(x) + tau_s}
z = H v
```

`H` 必须显式包含 `mu=v_biomass_C`、人工储备的有符号 source 净消耗，以及**全部有限库存**按实际 exchange stoichiometric sign 得到的净变化；不得在 `H` 中用 `max(0,...)` 隐藏方向冲突。诊断输入接口为 `state_id`、状态/模型/培养基 SHA、`pfba_contract`（primary objective、fraction、每个 forward/reverse variable 的 secondary weight、求解容差、fallback）、序列化的 `H`、`tau_s` 及分量绝对/相对容差。每个输出分量必须含 `component_id`、min/max solve status、min/max value、range、tolerance 和 `classification`，顶层再给 `ambiguous_component_ids` 与 `all_dynamic_components_unique`。两端都 finite/optimal 且 range 不超阈值为 `unique`，超阈值为 `ambiguous`；任何失败、无界或仅 fallback 得到的端点为 `indeterminate`。

本轮没有运行该诊断，因此结论只能是**尚未检测**。投影唯一只说明该状态的一阶动态导数不受内部多解影响，不证明随状态连续/光滑，也不证明更符合生物学。只有投影确有显著歧义且另获授权时，才考虑第三层固定 `v_ref=0`、固定正对角 `D` 的严格凸 QP；不得改成混合加权目标，也不得默认使用 previous flux。

## 验收问题 A–E

**A. 已核实/修改/未修改：**上表所列调用链、单位、cap、更新、终止和 pFBA 层级已静态核实。只新增一个标准库纯函数、一个三例测试及本修订证据；未改生产代码、模型、pFBA、库存政策、STOP 或历史文件。

**B. 改进大小：**常速率例从 Euler 的 `0.0256833 gDW/L` 误差降为 `0`；变速率的共同时间与事件误差每次减半约减半。真实 GEM 的改善、收敛和求解稳定性均未验证。

**C. `dt` 与最小接入缺口：**`dt` 确实进入真实优化边界，故减小它不只是积分细化。接入前只缺一个明确区分计划 cap 步长与事件接受子步长、在事件处重求通量且保留旧政策的薄适配层及其政策一致性验证；本轮不接线。

**D. pFBA 多解：****尚未检测**，既未证实也未排除。上节给出了保留真实 pFBA 两阶段最优集合的投影范围诊断，但本轮没有新增优化。

**E. 历史语义：**旧 Euler 默认、pFBA/fallback、库存 cap、STOP、历史结果和“步长未完全验证”状态全部保留。既有真实 gate 的科学 FAIL 也未重写或重跑。

## 唯一后续最小实模验证方案（仅提案，不执行）

若另获授权，只运行**一组 WT、一个冻结培养条件/初值/参数**的配对比较：原 pFBA+原 Euler 对比相同 pFBA+新一致积分/事件处理，分别用 `h=0.0625,0.03125,0.015625 h`，共 6 条轨迹。模型、培养基、初值、目标、Gurobi 参数、near-zero 与 STOP 判据必须逐字段/SHA 相同；历史基线仅在这些元数据完全匹配时复用。

预先验收：共同物理时刻 biomass/全部库存误差随两次减半下降；事件类型/顺序一致，事件时间差分别收敛；连续性和每库存原始收支残差 `<=1e-12 mmol/L`；任何修正必须 `<=1e-15+1e-14×scale` 且逐条记录；不得把 non-optimal 解释为死亡。资源上限为 6 条轨迹、1 solver thread、最多 5,000 次实际 backend optimize entry、总墙钟 1,800 s、0 次自动重试，触及任一上限即停止。

两臂都保留旧 cap 公式。每个 nominal interval 首次求解时，`h_planned_for_bounds` 等于该 nominal `h`；新臂若提前遇到事件，先推进 `h_accepted_for_integration`，关闭已耗尽库存，然后在同一网格区间内重求解，此时新的 `h_planned_for_bounds` 明确定义为“到当前 nominal 边界的剩余时间”，而不是刚接受的事件子步长。两者都必须逐段记录，绝不在首次求解后用已知的事件时间反向改 cap。

即使公式相同，事件后的重求解会产生旧 Euler 没有的可行域，因此这组比较应称为“原方法与一致积分/事件方法包的配对比较”，不能声称完全隔离了纯积分误差。跨 `h` 时旧 cap 本身也改变可行域，只能标记为“积分+库存政策耦合”的收敛检查。任何删除/替换 `C/(B*dt)` 的实验属于另一库存政策研究，不能混入这组比较。
