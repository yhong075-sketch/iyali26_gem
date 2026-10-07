# CoQ9 WP2 replay v2 候选计划（仅提案）

## Material Passport

- Origin Skill: experiment-agent
- Origin Mode: plan
- Origin Date: 2026-09-05
- Verification Status: UNVERIFIED（候选计划；本轮没有实施、测试或新求解）
- Version Label: wp2_replay_v2_plan_draft_20260905
- repro_lock: null（执行环境尚未获批、尚未生成新的 resolved manifest）
- Upstream Dependencies: 冻结 WP2 档案、独立只读审查、同目录参数记录与零事件 v2 提案。

## 1. 目的、证据口径与来源定位

研究问题：在科学输入和原计算逻辑不变、参数记录无损且比较规则单独版本化的前提下，原定八条件中哪些可观测量能够在新声明设置下重现，哪些仍不一致或没有历史比较依据？不预设三处 STOP 必然消失，也不以匹配率或 recall 选择参数。

本文件的“supported”表示已打开指定包内材料，其内容直接支持所述静态事实；引用旧独立审查的全量检查结果不等于本轮重新完成该全量审核。“partial”表示只能支持受限解释；“unverified”表示未证实或尚待执行；“contradicted”表示材料直接反驳该项声明。所有候选设计均为 **unverified／待用户决定**，不是既有事实或执行授权。

来源别名（均为归档成员定位，不要求这些路径在当前磁盘上存在）：

- `H`：外层 ZIP 根 `coq9_wp2_review_and_next_step_handoff_20260905/`。外层文件 `/Users/david/Downloads/coq9_wp2_review_and_next_step_handoff_20260905.zip`；SHA-256 `3d732a67cf3f90a9d26213cf2a160db491406f7492d855cf1ce06c42acef9001`。
- `A`：`H/coq9_wp2_independent_audit/`。
- `B`：`H/original_archives/coq9_wp2_qualified_baseline_20260905.zip` 内根 `coq9_wp2_qualified_baseline_20260905/`。该 ZIP 为 216015947 字节，SHA-256 `0b184c04a3a1270f413c4a5fbba60d386eedf9ea2560b1b5796a1fe905d90e16`。
- `T`：`B/execution_provenance/coq9_wp2_qualified_20260905T042505Z_inputs.tar.gz` 内根 `coq9_wp2_qualified_20260905T042505Z/`；tar.gz SHA-256 `38fdfd742f326365355d0ee0f1924d3b2bb69f00ebd9f0e8a305c60abf260bff`。
- `F`：`T/frozen_compute/`。文件和函数定位均指冻结文本，不能改用当前目录同名代码。

本次只读核查时间：2026-09-05。范围包括交接说明、提案提示与清单、独立审查及相关小表、原协议、输入清单、环境记录，以及输入 tar 中的配置和相关冻结代码文本；没有执行归档代码，没有导入模型或求解器，没有重算全部 state/flux。

## 2. 三个 KO 的身份与模型角色

WT 是同一冻结 PO1f 运行配置下不施加这三个单基因 KO 的对照，不是另一份模型或另一种菌株配置。

| 系统 ID | 已核实名称／符号 | 简要蛋白功能与证据等级 | 本计划保留的模型角色 | 事实状态与来源 |
|---|---|---|---|---|
| YALI1E18269g | 本材料未确立 Yarrowia 通用名称；COQ7 candidate | 候选去甲氧基泛醌羟化酶／CoQ 单加氧酶；本轮只采用 `model/GPR assignment only`，原生蛋白活性未核实 | R695 的 CoQ 合成羟化步骤 KO | **supported**：候选名称和模型赋值见 `B/prior_WP1_evidence/assay_cell_audit.tsv` 该 ID 行；`F/docs/research/quinone_gpr_synthome_2026-08-17/gene_evidence_matrix.tsv` 该行；原生活性 **unverified** |
| YALI1A14736g | 本材料未核实已确立名称 | 个体蛋白功能未独立确立；`model/GPR assignment only` | R305 的模型赋值组分／候选，冻结解释层标为 complex III；不由此断言原生催化身份 | **supported**：`B/prior_WP1_evidence/respiratory_evidence_matrix.tsv`，R305 的 `gpr_gene_identity_json`；原生身份 **unverified** |
| YALI1A21711g | 本材料未核实已确立名称 | 个体蛋白功能未独立确立；`model/GPR assignment only` | R2062 的模型赋值组分／候选，反应为 NADH:ubiquinone oxidoreductase | **supported**：同上表 R2062 的 `gpr_gene_identity_json`；原生身份 **unverified** |

这些 KO 的范围是原条件定义，不是本轮重新接受 GPR、重建序列身份或证明生物学机制。

## 3. 不得改写的已有结论

| 事实／待检验声明 | 状态 | 精确范围与来源 |
|---|---|---|
| 原 WP2 是“3 条全可比指标匹配／2 条单条件匹配／3 条停止” | **supported** | `A/INDEPENDENT_REVIEW_zh.md` §3；`A/condition_verdicts.tsv`。这不是八条件完整历史环境复现 |
| 原 r2 对应八条件参考有 2305 行，而 WP2 v1 实际只保存 1760 行 | **supported** | `A/original_reference_check.json`；`A/condition_verdicts.tsv`；`A/INDEPENDENT_REVIEW_zh.md` §2–3。1760 行含 1758 推进区间和 2 个非推进终止行；545 行是 WP2 v1 未执行尾段，不是原 r2 没有轨迹 |
| 三个原 STOP 的连续库存差异低于原库存容差 | **supported** | `A/exact_zero_points.tsv`，三个条件及其原始数值；不能因此把 STOP 改为 PASS |
| 三个原 STOP 的 `raw_exact_zero` 相同 | **contradicted** | 同表及 `B/analysis_code/coq9_wp2_replay.py:107` 的 `compare_row`；精确布尔差异仍成立 |
| 差异与近零数值边界效应相容 | **partial** | `A/INDEPENDENT_REVIEW_zh.md` §4；这是解释，不是已证明的因果机制 |
| 未执行尾段行为相同，或三处差异证明生物学回归 | **unverified** | 同报告 §4；不能从现有前缀外推 |
| WP2 v1 已有完整、无损的当前参数记录 | **contradicted** | `A/parameter_and_preparation_check.json`；`B/current_environment_recording_caveat.json`：151 个参数中 14 个值有损成为 null，覆盖全部八条件 |
| 14 个 null 直接污染了求解器参数 | **unverified** | `A/SOURCE_CODE_EXCERPTS.txt`，`solver_parameters`、`clean`、`CURRENT_SETTINGS` 设置路径：未见把 null 读回设置求解器的反馈路径；不能把记录缺陷当作已证明的参数污染 |
| 历史完整环境、历史全通量、原生功能已经复现 | **unverified** | `B/PREDECLARED_PROTOCOL.json` 的环境和不可用比较字段；`A/INDEPENDENT_REVIEW_zh.md` §5–7 |

14 个原值若不能从同一次 WP2 v1 执行的直接旁证恢复，另外五条旧运行也仍缺完整无损环境记录。新运行捕获、新查询或默认值都不能补写成旧执行值。修复后的“完整新记录”与“恢复旧记录”是两件事。

## 4. 三条与八条：不同的证据目标

| 候选范围 | 可能支持的声明（必须真正运行并通过后） | 仍不能支持的声明 |
|---|---|---|
| 严格只新跑原三条 STOP：finite_batch WT、finite_batch YALI1A14736g、po1f_nonlimiting YALI1A21711g；各自从 t=0 开始 | 三条条件在新声明设置、参数 schema 和零事件 v2 下的观察结果；与原 r2 可用字段及 v1 已执行前缀的受限比较；若通过，可补充观察过去未执行区间 | 不补齐其余五条的有损参数记录；不是八条同版新基线；没有新 po1f_nonlimiting WT，故该模式新 KO 的归一化比/cutoff 必须保持不可用 |
| 原定八条全部新跑：两模式各 WT＋三个 KO；按下节分批 | 八条都有同一版新参数记录和比较规则的候选参考；两模式新 WT 各自通过后，才可报告相应新 KO/WT 与 cutoff | 仍不等于历史完整环境精确复原；历史未存的全通量、12 种其他逐步有限库存、pFBA 次级目标与完整 LP 状态不能凭新结果补造；不代表正式模型发布、独立生物学验证或校准完成 |

此比较是 **unverified／候选设计推论**，直接依据 `H/next_step/CODEX_PROPOSAL_ONLY_PROMPT_zh.txt` 交付3、`A/INDEPENDENT_REVIEW_zh.md` §5、`B/PREDECLARED_PROTOCOL.json` 的历史比较限制。三条方案不得从旧 STOP 的 state 热启动接上“新尾段”，应保留独立 run 身份并从 t=0 开始。

若“三条”方案还要求两种模式各有新 WT，至少成为 **4 条新运行**，必须明确新增这一对照，不能仍称“只补三条”。本文件建议提交审批的是下述原定八条件分批方案；用户也可以选择严格三条的更窄证据目标，二者不同时默认获批。

## 5. 原定八条件候选分批顺序

| 批次 | 新条件（串行、各执行一次） | 原 WP2 v1 观察／参考规模 | 下一批闸门 |
|---|---|---|---|
| B1：两个 WT | finite_batch WT；po1f_nonlimiting WT | 前者 57/384 行，step56 STOP；后者 86/86 行，85 个推进区间＋1 个 5.3125 h infeasible 终止行 | 两种 WT 分别完成输入/参数/行为审核；先交人工审阅，不能自动进 B2 |
| B2：其余两个原 STOP | finite_batch YALI1A14736g；po1f_nonlimiting YALI1A21711g | 分别 167/384 行，step166 STOP；85/86 行，step84 STOP，下一非推进终止行在 v1 未执行 | 对应模式新 WT 已通过；两个原边界之后的新观察与有效历史比较完成后人工决定 B3 |
| B3：剩余四个 KO | finite_batch YALI1E18269g；finite_batch YALI1A21711g；po1f_nonlimiting YALI1E18269g；po1f_nonlimiting YALI1A14736g | 前三条各 384/384 行并至 24 h；最后一条 213/213 行，212 个推进区间＋1 个 13.25 h infeasible 终止行 | 每条独立判定；全八条汇总仍需人工接受，不自动进入下一工作包 |

旧观察为 **supported**，来源 `A/condition_verdicts.tsv`；新批次顺序是 **unverified／候选设计**，来源要求见交接提示交付3。批次是审核闸门，不是参数扫描；不重复 B1 的 WT，不用不同环境的 WT 拼接分母。

“通过 WT”指：该模式完整覆盖其对应历史轨迹的可比较记录，达到原 24 h 或原逻辑的提前终止点，输入、参数、离散行为和适用数值审核无未解决阻断。WT 以同一时刻、同一状态正常复现 infeasible 终止，仍可以形成受限的末端 doublings 对照；必须显示提前终止，而不是写成“完成24 h”。WT 若发生比较 STOP、资源终止、非有限或非正的 dynamic_doublings，不能提供新分母。

## 6. 必须冻结的输入与计算语义

### 6.1 文件、代码与运行时配置

以下身份为 **supported** 的归档记录；来源 `B/replay_input_inventory.json`、`B/PREDECLARED_PROTOCOL.json`、`T/INPUT_SHA256SUMS` 和指定配置文本。未来执行需再次逐项核对实际读入文件，不能仅引用本表。

| 对象 | 唯一候选输入／完整 SHA-256 |
|---|---|
| 模型 | `T/inputs/model.xml`；`bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee` |
| 培养基 | `T/inputs/sd_leu.csv`；`ed176d26a373f98cc413ed2e32a71f5f060a06e343f90f7db25cd32eff268e85` |
| 菌株运行配置 | `T/inputs/po1f_sd_leu.json`；`35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383` |
| 既有实验标签输入 | `T/inputs/consensus_essential_genes.csv`；`1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b`；只冻结，不重标或用 recall 选参 |
| compute 代码身份 | commit `36bb6f0735e4c6458bd53c0ceb01952b116b8be7`；`F/scripts/gem_annotate/quinone_dfba_essentiality.py`：`679ada071adb8e1c898ed1b60d9d8c3f895dac5ca56e1e4c19b1ea28ae6eb5eb` |
| CoQ9 dilution 工具 | `F/scripts/gem_annotate/coq9_dilution.py`；`f223faa225c9c748ecae3d40c80d1b903e90e11fa782b32cd7946619dc8d4ff4` |
| 原 observer／comparator | `T/runner/coq9_wp2_capture.py`：`3db5da2ed71994990395d0e7a738ce7d3dab35473576aabb7ba7c52abcd552fa`；`T/runner/coq9_wp2_replay.py`：`7b75c75344613ed665fc355e16eff0bd67447fbf89a514ce8ab79ace6cd1a166`；未来仅获批副本可改，另记新 SHA |
| 原预声明协议 | `B/PREDECLARED_PROTOCOL.json`；`b3c0d8b71a6266bd3a9ccf4ce71615811bb899b5f66c5ac2e28b940f2f37b97c`；原件不变 |
| 原 r2 轨迹切片／calls | `B/reference/baseline_reference_intervals.tsv.gz`：`959387b366ae95b1f292e109fee6ca4f75d532612878828c78b6b013ec78326e`；`B/reference/reference_calls.json`：`185d8988ca3ca7574d7b653190e0d0fb59242511a1c9abaa238c477dd0726c1c` |

未来还应保留协议列出的全部 66 个冻结 Python 文件逐文件哈希，而不是以主 runner 匹配替代完整代码身份；额外记录新 sidecar、依赖和实际解析路径。当前项目状态中的静态 WT＋固定6KO结果属于另一条工作线，不能代替本 CoQ9 的动态 WT、输入或证据门；当前目录 model.xml 和其他候选也不得替换上述冻结模型。此区分为 **supported**，来源当前 `PROJECT_STATE.md`“目录与模型身份”“本轮限定静态重放”。

运行配置保持 profile `po1f_sd_leu_accrispr_v1`、`strain_overlay_enabled=true`；实际有效上下文指纹应匹配 `c243b23e7344e3f1e2b4962be25f0f2a38980990c6fe88ee32d3aa4f7af90e30`，overlay 效果 SHA 应匹配 `d15acbde9438f5d2391c4da23705a34a3585833062d616517d5af052088606c2`。培养基文件的 R1354=0.01607 先由原 profile 检验，再以原运行时 override 得到 1000；两模式差异在于是否移除 R1354 有限库存，而不是改成另一培养基。profile 中已有的菌株 overlay 全部原样保留，不在此新接受或扩大 GPR 变更。来源：`B/replay_input_inventory.json:effective_context/overlay_audit`、`T/inputs/po1f_sd_leu.json`；**supported**。

### 6.2 固定科学参数与库存

固定 `alpha=1e-4 mmol/gDW`、`pool_multiplier=1`、`dt=0.0625 h`、`B0=0.01 gDW/L`、名义上限 `24 h`；保留原 explicit Euler 与提前终止逻辑。初始人工 reserve 为 `alpha×B0×pool_multiplier`，浮点记录 `1.0000000000000002e-6 mmol/L`；它是 H-Q9-1 敏感性编码的人工 reserve，不是实测催化性 CoQ 池。`runtime_topology=false`、`runtime_gpr_scenario=null`，不启用新拓扑、GPR 情景或 reserve rescue。**supported**：`B/PREDECLARED_PROTOCOL.json`、`B/reference/replay_input_inventory.json:authorized_baseline/all_recorded_mode_settings`。

原 finite_batch 的初始有限库存（mmol/L）如下；来源 `F/scripts/gem_annotate/quinone_dfba_essentiality.py:116` 的 `INITIAL_POOLS_MMOL_L` 与上述 mode settings，**supported**。

| 反应 | 初值 | 反应 | 初值 |
|---|---:|---|---:|
| R1070 | 111.0 | R1003 | 0.054 |
| R1202 | 0.287 | R1204 | 0.601 |
| R1215 | 0.095 | R1217 | 0.381 |
| R1220 | 0.274 | R1222 | 0.134 |
| R1223 | 0.303 | R1231 | 0.840 |
| R1232 | 0.245 | R1233 | 0.276 |
| R1234 | 1.195 | R1354 | 0.178428 |

po1f_nonlimiting 保留其中 13 个有限库存，移除 R1354；其 uracil start/end 的 NaN 是原“无有限库存”的表示，不是缺失测量、零库存或新零事件。不能用 profile 的 `0.1784280489` 替换 runner 的 `0.178428`，也不能按分子量重新算 glucose 初值替换 111.0。新 schema 对“不适用”提供明确元数据，同时保留原轨迹表示。**supported**：同上；`B/reference/replay_input_inventory.json:nonlimiting_uracil`。

计算语义保持：原 pFBA 策略及其内置 `OptimizationError`→原 `model.optimize()` 后备；source-free growth `>1e-9 h⁻¹` 时不用 reserve；source-free 非 optimal 时立即返回，不能新增 reserve rescue；原 finite-medium 的 `amount/(B_start×dt)` 限制、库存非负更新和 Q9 原固定 clamp 均不修改。前瞻标签只读，不进入库存、bounds 或积分更新。**supported** 的代码定位：冻结 runner 的 `_finite_medium:429`、`_optimize_minimal_pool:438`、`simulate_gene:538`。保留原内置后备不授权额外重试、调参或替代求解流程。

## 7. 新参数 schema、环境与规则的前瞻冻结

本节是 **unverified／候选设计**。参数无损表示、实际 solver model 对象身份、捕获时间及各求解阶段检查，统一依照同目录 `PARAMETER_CAPTURE_FIX_PROPOSAL.md`：候选命名空间 `coq9.wp2.solver_parameters/v2`，`parameter_capture_schema_version=2`，分别保留值指纹 `parameter_values_sha256_v2` 与证据事件指纹 `capture_record_sha256_v2`。零事件比较统一依照 `ZERO_EVENT_POLICY_V2_PROPOSAL.md` 的 `zero_event_policy_v2_proposed`。两者均与冻结 trajectory schema 1.8 分离；两文档未经实施、合成测试与人工接受，本计划不能进入计算。

| 项目 | 旧记录能支持什么 | 未来候选要求 |
|---|---|---|
| 基本运行时 | **supported**：历史已记录 Python 3.11.15、cobra 0.32.1、Gurobi/gurobipy 10.0.3 | 保持这些版本；执行前核实实际解释器、包与 solver engine 版本。不同版本应停止，不自动升级或声称等价 |
| WP2 v1 额外包记录 | **supported**：optlang 1.9.1、numpy 2.4.6、scipy 1.17.1、pandas 2.3.3、python-libsbml 5.21.1 | 将实际完整依赖清单单独记录；这些是 v1 当前环境记录，不回填成 r2 已知环境 |
| 显式求解参数 | **supported**：v1 显式设置 Threads=1、Seed=0、Method=1、FeasibilityTol=1e-9 | v2 候选继续显式设定并逐阶段读回这四项 |
| 其他两个已读回值 | **supported**：v1 `OptimalityTol=1e-7`、`Presolve=0` 是首解前当前有效值，非 r2 已知历史值 | v2 候选将两项纳入运行前和阶段读回一致性闸门，不新增setParam；若不符即停，不自动设值使其匹配 |
| 其余非敏感执行参数 | v1 14 项有损；r2 除已记录项外 unknown 仍是 unknown | 先形成逐项、类型保真、明确实际值和来源的新 resolved manifest，尤其资源/iteration/cutoff 参数。未决值不得靠本文或静默默认补齐；manifest 通过人工接受前不求解 |

旧记录来源：`B/historical_vs_current_environment.json:historical_recorded_versions/current_packages/historical_vs_current_parameters/historical_solver_evidence`；`B/current_environment_recording_caveat.json`；`B/analysis_code/coq9_wp2_replay.py:503`。未知不是运行时默认值的同义词。CPU/OS、机器、线程环境、实际 solver 对象、包装层配置及其与底层参数的关系另记新观察；不假称重建旧硬件或完整历史环境。

新运行声明至少绑定：本计划选择的3条或8条清单、上述所有输入 SHA、原 compute SHA、新 sidecar SHA、参数 schema 版本、类型保真的 resolved-parameter fingerprint、零事件规则版本及文档 SHA、完整容差表、批次顺序、预算和停机判定。只比较相同版本、同一捕获阶段的指纹；不把“有损旧指纹相等”当作旧参数完整或不变的证明。

零事件 v2 的候选尺度为 `epsilon_C=min(原该类库存 atol, 1e-8 mmol/gDW/h × B_start × dt)`，仅对有效推进区间使用其真实起点和固定dt，成对比较取两条轨迹各自阈值的较小值；具体适用性、负库存、非推进终止及描述性差异依赖同目录零事件提案。`B_start` 不能偷换为审查表量级诊断所用的 `B_end`；不得根据三个观测误差重新选值。该规则是观察到 v1 STOP 后提出的新版本，不冒充原预注册。非推进terminal row不计算或继承band，其适用库存raw exact-zero差异仍STOP；不能预先保证三个旧STOP之后都可通过。

## 8. 预声明比较层与原容差

候选 v2 不扩大下表数值容差，沿用原 `abs(new-ref) <= atol + rtol×abs(ref)`。来源 `B/PREDECLARED_PROTOCOL.json:comparison_tolerances`、`B/analysis_code/coq9_wp2_replay.py:99`，**supported**；未来采用为待批准设计。规则变更只按独立 v2 零事件提案声明，不能改写原比较结果。

| 比较量 | atol | rtol | 单位／限制 |
|---|---:|---:|---|
| 时间 | 0 | 0 | h；推进时间与终止时间精确比较 |
| biomass | 1e-10 | 1e-9 | gDW/L |
| 有限库存 | 1e-10 | 1e-9 | mmol/L |
| Q9 reserve | 1e-12 | 1e-9 | mmol/L |
| growth、source-free growth、主目标 | 1e-8 | 1e-8 | h⁻¹ |
| flux | 1e-8 | 1e-8 | mmol/gDW/h；不同通量角色分开 |
| dynamic doublings／ratio | 1e-8 | 1e-8 | 无量纲 |
| pFBA secondary objective | 1e-6 | 1e-8 | 求解器次级目标；仅在确有可比存档时比较 |
| dilution coupling residual | 1e-9 | 0 | 原 coupling 方程单位 |
| Sv／bounds 数值残差 | 1e-7 | 0 | 逐类型、按原审计字段记录；不能当作全部都满足 1e-9 |

分层输出、分层判定（全部为 **unverified／待执行设计**）：

1. **目标与连续状态**：比较主目标、growth、source-free growth、biomass、Q9、历史保存的 glucose/uracil、时间和末端 summary；新记录中另审全部有限库存、Euler 连续性、Q9 原 clamp、dilution coupling、Sv、bounds 和 oxygen/ATP 账本算术。数值一致只支持该数值检查，不证明化学正确、原生机制或独立最优性。
2. **类别、分支与事件**：`step_index`、`interval_advanced`、求解状态、source-free 状态、选中求解阶段／分支、growth `<=1e-9` 布尔事件、终止状态和时间分别检查。原值、`raw_exact_zero`、旧 STOP 永久保留；near-zero 是新描述层，不替代行为比较。原 `q9_pool_depleted_h` 精确事件诊断单列，不能被 near-zero 时间替代。非 optimal 的不推进终止行不能用 near-zero 标签“救回”、推进或改成24 h。
3. **完整 flux／LP**：r2 历史只给已归档面板/核心 flux；v1 前缀已有完整 state、attempt primal、全反应 flux、次级目标和 LP，可用于 v2 与 v1 在同条件、同 step、同求解阶段的比较。v1 STOP 后没有这些全量记录；r2 对其他12种逐步有限库存、历史完整 flux、pFBA 次级目标和完整 LP 均无档案，明确记 `comparison_unavailable`，不重建隐藏变量充作真值。
4. **内部 flux 差异**：单个内部坐标不符单列坐标、尺度、LP身份、所选阶段、主/次目标、状态与下游影响，不能自动解释为新机制，也不能宣布全部通量一致。若只影响内部替代解声明，可保留其他层的受限结论；但库存驱动 uptake、source/dilution等对应通量比较超出原容差，或目标/动态状态/终止超出各自预声明容差，或分支等离散闸门不符，或数值残差审核超限时必须 STOP。不能把容差内微差一律写成行为分歧，也不能用“替代最优解”掩盖真正超限。历史缺少次级目标时，不能声称已经证明替代最优解。
5. **新的 KO/WT**：先要求相同模式新 WT 通过且其 `dynamic_doublings` 有限、严格大于0，再计算 `KO dynamic_doublings / WT dynamic_doublings`。保留1/5/10/15%严格 `<` cutoff。新 WT 未通过时，KO 单条件结果可以保存，但新 ratio 和新 cutoff 均为不可用，不借用 r2、v1、另一模式或静态重放 WT。新条件发生资源/求解失败时不得用零 growth 或“死亡”补值。

历史比较可用性依据为 **supported**：`B/PREDECLARED_PROTOCOL.json:unavailable_historical_comparisons`；`A/INDEPENDENT_REVIEW_zh.md` §3、§6；归一化和 cutoff 依据为 **supported**：`B/analysis_code/coq9_wp2_replay.py:150`。这些定义不是全向量一致性的保证。

## 9. 预算、失败与人工闸门

以下是 **unverified／尚未授权的资源上限提案**，不是实际消耗估计，也不是允许立即运行的指令：单计算进程、单 solver 线程、条件串行；八条方案最多8个新轨迹，每条最多384次原积分循环，每条墙钟最多30分钟，全部求解累计最多2小时；内存上限4 GiB，新增结果上限2 GiB。原冻结算法自身的 pFBA／后备调用分别记录实际底层求解次数；不为凑预算跳过阶段、削减审计记录或改用另一算法。预算不足即保存可用证据、报告并停下；不自动延长或重试。严格三条方案同样最多3个轨迹、各384循环，不因较窄范围放宽单条限制。

预算应由用户在知道保存粒度和执行环境后单独接受；未来运行前还需解决上一节的完整参数 manifest，不能把资源上限擅自转写为旧 solver 参数。当前没有 HPCC 授权，本文件不提供提交命令。若未来选择集群，目标系统、资源和提交动作另行获得明确许可。

人工闸门按顺序保留：

- **G0：提案决策**。用户决定三条还是八条证据目标，分别决定参数记录/比较器的最小代码修改范围、仅合成数据测试范围；此时仍不授权求解。
- **G1：实施与合成测试审阅**。仅在相应范围另获批后才实施和测试；具体测试设计见两个 sibling 提案，本轮不执行。新 compute 身份不得变化；修改只限获批的记录/比较层副本。通过软件测试不等于通过科学 replay。
- **G2：计算授权**。用户另行接受完整 resolved manifest、输入身份、规则版本、输出位置、预算、停止规则和 B1 范围；没有精确新值或捕获路径仍有缺失，停在此门。
- **G3/G4：分批接受**。B1 完成后人工决定 B2，B2 完成后人工决定 B3；任何一个模式 WT 未通过，只阻断依赖它的 KO/WT 结论及相应后续计算，不能由另一个 WT 补位。全局输入/参数记录失效则阻断所有依赖它的条件。
- **G5：结果接受**。独立只读审阅者打开新证据，分别审核事实、数值比较、环境完整性和未检查项；用户决定是否接受为“新声明设置下的八条件受限参考”。不自动接受模型、修复病例或进入下一工作包。

必须停止并保存事实：输入/源码/配置 SHA 不符；实际模型对象或参数缺失、类型丢失、发生未声明变更；无法确认实际求解阶段；负库存或非有限状态触发 v2 停止规则；连续量超容差；状态、source-free 分支、growth 事件、推进或终止行为不一致；未声明的 resource/iteration/cutoff 终止；rollback 失败；预算到限；需要改科学输入或放宽标准才能匹配。普通已记录边界差异与真正 STOP 的区分只按已获批的 v2 规则，不临场解释后放行。

source-free infeasible 的原逻辑终止仍如实记录；无可用最优解不等于生物学死亡。不得新增 dt/alpha/pool 矩阵、FVA、maintenance/glucose 反事实、基因扫描、模型修复或下一工作包。冻结包中出现准备对象、脚本或旧计划不构成这些动作的授权。

## 10. 未来交付物与当前停止点

若未来获得相应授权，复用原交付结构，在新的不可覆盖 run 目录产生：版本化预声明协议、实际输入/代码清单、无损逐阶段参数和环境记录、逐条件原始轨迹与终止信息、完整 state/attempt/flux/LP 及 oxygen/ATP 账本、分层比较表、旧 STOP 对照、新 ratio 可用性表、资源与停止日志、SHA256SUMS。同时给出每条件参考/实际/缺失行数、每类可比较字段数量和未检查范围；2305 行只作为原 r2 八条件参考规模，不预填为新运行已完成行数。

目前仍未解决：14 个旧值是否存在同次执行直接旁证；完整新参数 manifest；v2 规则与预算是否获用户接受；三个原 STOP 之后的实际新行为；历史未归档全量状态/次级目标；原生基因功能和化学/机制证据。以上缺口不妨碍提交本计划，但阻止扩大对应复现或生物学声明。

**当前停止点：仅提交候选计划。没有实施参数捕获修复或零事件规则，没有运行合成测试、FBA/dFBA/FVA、HPCC 或任何新条件；没有修改原 STOP、模型或封存数据，没有提交或推送。后续决定见同目录 `DECISION_ITEMS.md`。**
