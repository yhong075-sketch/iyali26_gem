# CoQ9 WP2 replay v2 候选计划（proposal_v2；仅文档）

## Material Passport

- Origin Skill: experiment-agent / govern-agentic-research
- Origin Mode: plan
- Origin Date: 2026-09-05
- Document Revision Date: 2026-09-06
- Verification Status: UNVERIFIED（候选设计；未实施、未测试、未操作真实求解器）
- Version Label: proposal_v2
- document_schema_version: 2
- repro_lock: null（实际对象预检和新的 resolved manifest 尚未获批、未生成）
- Upstream Dependencies: 冻结 WP2 档案、proposal_v1、CoQ9 只读审查意见，以及同目录参数记录、零事件与决策提案。

本次修订落实 R4–R6，并同步 R1–R3 的约束；不扩大为代码实施、合成测试、真实求解器操作、3/8 条重放或下一工作包授权。旧提案、原始 null、旧 STOP、旧 manifest、轨迹和封存 SHA 均不改写。

## 1. 问题、来源及事实等级

问题是：在科学输入与冻结计算逻辑不变、记录和比较政策另行获批的前提下，哪些可观测量能在新声明设置下匹配，哪些仍不一致或不可比较？不预设三个旧 STOP 必然消失，不以 recall 调整参数。

`supported` 表示指定来源支持所述归档事实；`partial` 表示仅支持受限解释；`unverified` 表示未确认或未执行；`contradicted` 表示来源反驳该声明。下面归档事实从 proposal_v1 的带定位材料保留，不声称本次重新数值复核全部状态。全部新设计均为 **unverified／待决定**。来源别名及文件身份应与同目录 `document_manifest.json` 一起使用：

- `H`：`/Users/david/Downloads/coq9_wp2_review_and_next_step_handoff_20260905.zip` 内根 `coq9_wp2_review_and_next_step_handoff_20260905/`；外层 SHA `3d732a67cf3f90a9d26213cf2a160db491406f7492d855cf1ce06c42acef9001`。
- `A`：`H/coq9_wp2_independent_audit/`。
- `B`：`H/original_archives/coq9_wp2_qualified_baseline_20260905.zip` 内根 `coq9_wp2_qualified_baseline_20260905/`；内层 ZIP SHA `0b184c04a3a1270f413c4a5fbba60d386eedf9ea2560b1b5796a1fe905d90e16`。
- `T`：`B/execution_provenance/coq9_wp2_qualified_20260905T042505Z_inputs.tar.gz` 内根 `coq9_wp2_qualified_20260905T042505Z/`；tar.gz SHA `38fdfd742f326365355d0ee0f1924d3b2bb69f00ebd9f0e8a305c60abf260bff`。
- `F`：`T/frozen_compute/`。代码定位均指冻结文本，不得换成当前工作区同名代码。
- `V1`：`../coq9_wp2_proposal_only_20260905T2000Z/WP2_REPLAY_V2_PLAN.md`；版本身份见新 manifest。
- `R`：用户附件 `/Users/david/Downloads/TO_INACTIVE_NEXT_STEPS_zh.txt` 的 R1–R6 及“必须原样保留的结论”；附件 SHA 与本轮核查范围见新 manifest。

CoQ9 审查方明确只读了转交文本及其已有来源，没有独立核验新提案 ZIP 或六个文件的 SHA。本文不把来源代理或本地核查冒充 CoQ9 聊天已完成的核验。项目其他工作线的状态不是本计划输入：本修订删除 v1 对未绑定版本的 `PROJECT_STATE.md` 静态重放数值断言；本轮有限来源摘录及身份核查见 `CHANGE_RESPONSE.md` 的“来源与身份片段”，不据此声称其静态计算结果已在本轮复现，也不将日常工作区模型替换为 CoQ9 输入。

## 2. 三个 KO 的身份及现有结果

WT 是同一冻结 PO1f 配置中不施加这三个 KO 的对照。

| 系统 ID | 已核实名称／符号 | 简要蛋白功能与证据等级 | 模型角色和来源 |
|---|---|---|---|
| YALI1E18269g | 本材料未确立通用名称；COQ7 candidate | 候选去甲氧基泛醌羟化酶／CoQ 单加氧酶；model/GPR assignment only，原生活性 unverified | R695 的 CoQ 合成步骤；`B/prior_WP1_evidence/assay_cell_audit.tsv` 及 `F/docs/research/quinone_gpr_synthome_2026-08-17/gene_evidence_matrix.tsv` 该 ID 行，模型赋值 supported |
| YALI1A14736g | 本材料未核实已确立名称 | 个体蛋白功能未独立确立；model/GPR assignment only | R305 的模型赋值组分／候选，解释层为 III 型红氧角色；`B/prior_WP1_evidence/respiratory_evidence_matrix.tsv` R305 的 `gpr_gene_identity_json`，模型赋值 supported，原生身份 unverified |
| YALI1A21711g | 本材料未核实已确立名称 | 个体蛋白功能未独立确立；model/GPR assignment only | R2062 的模型赋值组分／候选，反应为 NADH:ubiquinone oxidoreductase；同表 R2062 的 `gpr_gene_identity_json`，模型赋值 supported，原生身份 unverified |

| 保留事实／声明 | 等级 | 来源与限制 |
|---|---|---|
| 旧运行仍为 3 条全部适用归档指标匹配、2 条单条件匹配、3 条原规则 STOP | supported | `A/condition_verdicts.tsv`、`A/INDEPENDENT_REVIEW_zh.md` §3；不是八条件完整历史环境复现 |
| r2 八条件参考 2305 行，WP2 v1 保存 1760 行，其中 1758 推进、2 非推进；545 行尾段未在旧 WP2 执行 | supported | `A/original_reference_check.json`、上述 verdicts；不得补造、拼接或离线改为 PASS |
| 三处库存差低于原连续容差，但 raw exact-zero 不同 | supported | `A/exact_zero_points.tsv`；“原 raw-zero 一致”为 contradicted |
| 与近零数值边界效应相容 | partial | `A/INDEPENDENT_REVIEW_zh.md` §4；因果机制及未执行尾段一致性 unverified |
| 14 个旧参数原值未恢复，缺口覆盖全部八条件 | supported | `A/parameter_and_preparation_check.json`、`B/current_environment_recording_caveat.json`；151 是旧记录观察数量，不是永久完整性标准 |
| 旧参数 null 已污染求解器设置 | unverified | `A/SOURCE_CODE_EXCERPTS.txt` 的读取、clean、设置路径未显示反馈 null 的证据；记录缺陷不等于参数污染 |
| 历史全环境、全量通量、隐藏状态或原生功能已复现 | unverified | `B/PREDECLARED_PROTOCOL.json` 的不可用比较项；新记录不能恢复旧执行未知值 |

人工 reserve 是 H-Q9-1 敏感性编码，不是实测催化 CoQ 池；原生蛋白功能、氧摄取来源和呼吸补偿机制仍按原证据等级记录，不新增生物学结论。

## 3. 三条与八条：不同证据目标

| 候选范围 | 真正获批、执行且通过后最多支持 | 不能支持 |
|---|---|---|
| 严格 3-only：finite_batch WT、finite_batch YALI1A14736g、po1f_nonlimiting YALI1A21711g，各从 t=0 开始 | 三个原 STOP 条件在新声明设置下的结果；与 r2 可用字段、v1 已执行前缀的分层比较 | 不补齐其余五条参数；不是八条完整新记录；没有新的 nonlimiting WT，该模式新 KO/WT ratio 和 cutoff 为 NA |
| 全部原定 8 条：两模式各 WT＋三个 KO | 八条同版新记录；各模式新 WT 通过后可生成对应新 ratio | 不是历史完整环境复原；不能补造旧隐藏状态或证明生物学校准 |

如果 3-only 还需要新增 nonlimiting WT，至少变成 **4 条**，必须单独列入 `new_computation_scope`；不能隐含增加。不得从旧 STOP 快照热启动拼接新尾段。14 项缺口涉及八条，因此建议只有在目标明确是“完整八条新记录”时选择八条；不按哪种更易通过选择范围。

### 3.1 八条候选分批和对象生命周期（R5）

| 批次 | 条件（串行、每条一次） | 原观察；未来放行规则 |
|---|---|---|
| B1 | finite_batch WT；po1f_nonlimiting WT | 前者 step56 STOP；后者 85 推进＋5.3125 h 非推进 infeasible；完成后人工审核，不自动进入 B2 |
| B2 | finite_batch YALI1A14736g；po1f_nonlimiting YALI1A21711g | 分别 step166、step84 STOP；各对应新 WT 需通过适用闸门；完成后人工审核 B3 |
| B3 | finite_batch YALI1E18269g；finite_batch YALI1A21711g；po1f_nonlimiting YALI1E18269g；po1f_nonlimiting YALI1A14736g | 前三条原为 384 推进至24 h；末条212推进＋13.25 h非推进 infeasible；八条结果接受不自动授权下一个工作包 |

原观察 supported，来源 `A/condition_verdicts.tsv`；以下生命周期为未实施候选：

1. 每批启动一个新的 worker 进程。每种 mode 创建新的有效 simulation context/model/native solver 对象；同 mode 多 KO 可沿冻结 `with model` 回滚路径复用该批该 mode 对象，不新增优化或重置调用。
2. 记录 PID、worker UUID、mode/context/model/native object ID、创建时间及每个 KO 前后回滚身份。每个 native 求解阶段使用实际对象重新捕获；不能用 mode 首条缓存代替。
3. 批次结束释放对象，跨批不复用模型对象、native 实例或隐藏基；不同进程可以属于同一获准运行计划。B3 引用 B1 新 WT 的已封存、已过闸门结果，不是复用其 solver 对象，也不重跑 B1 WT。
4. 同 mode 对象内可能存在求解器内部热状态；本计划不声称全部可见或可哈希，不擅自清基、改算法或额外求解。v1 与候选 v2 对象复用/内部隐藏状态的一致性若不能由旧记录证明，应保留 `historical_object_state_equivalence=not_established`。参数 hash 相同不证明所有内部状态相同。
5. 回滚失败、实际对象身份或参数读回缺失、未声明对象替换，停止受影响条件；不创建一个“看起来相同”的对象继续旧轨迹。生命周期政策需随未来计算范围批准，不能称已复现历史执行环境。

## 4. 固定输入与计算语义

以下身份是 v1 带来源的归档记录（supported），不是本轮新计算结果。未来必须核对实际读入字节；缺失不能猜测或用当前 main 自动替换。来源：`B/replay_input_inventory.json`、`B/PREDECLARED_PROTOCOL.json`、`T/INPUT_SHA256SUMS`。

| 输入 | 固定身份 |
|---|---|
| 模型 `T/inputs/model.xml` | `bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee` |
| medium `T/inputs/sd_leu.csv` | `ed176d26a373f98cc413ed2e32a71f5f060a06e343f90f7db25cd32eff268e85` |
| strain profile `T/inputs/po1f_sd_leu.json` | `35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383` |
| consensus 标签输入 | `1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b`；不改标签、不用 recall 选参 |
| compute commit | `36bb6f0735e4c6458bd53c0ceb01952b116b8be7` |
| 冻结 `quinone_dfba_essentiality.py` | `679ada071adb8e1c898ed1b60d9d8c3f895dac5ca56e1e4c19b1ea28ae6eb5eb` |
| 冻结 `coq9_dilution.py` | `f223faa225c9c748ecae3d40c80d1b903e90e11fa782b32cd7946619dc8d4ff4` |
| v1 capture / replay | `3db5da2ed71994990395d0e7a738ce7d3dab35473576aabb7ba7c52abcd552fa` / `7b75c75344613ed665fc355e16eff0bd67447fbf89a514ce8ab79ace6cd1a166`；原件不变，未来获批副本另记新 SHA |
| 原协议 | `b3c0d8b71a6266bd3a9ccf4ce71615811bb899b5f66c5ac2e28b940f2f37b97c` |
| r2 参考 intervals / calls | `959387b366ae95b1f292e109fee6ca4f75d532612878828c78b6b013ec78326e` / `185d8988ca3ca7574d7b653190e0d0fb59242511a1c9abaa238c477dd0726c1c` |

还需绑定原清单全部66个冻结 Python 文件、新获批 sidecar 和实际解析路径；单 runner SHA 一致不足以证明整个环境重现。trajectory 保持 schema 1.8，新记录 schema 和比较政策单独版本化，不改旧 manifest。

固定：`alpha=1e-4 mmol/gDW`、`pool_multiplier=1`、`dt=0.0625 h`、`B0=0.01 gDW/L`、`hours=24`；`runtime_topology=false`、`runtime_gpr_scenario=null`。初始人工 reserve `alpha×B0×pool_multiplier`，原浮点记录为 `1.0000000000000002e-6 mmol/L`。

保持 profile `po1f_sd_leu_accrispr_v1`、overlay enabled；有效 context SHA `c243b23e7344e3f1e2b4962be25f0f2a38980990c6fe88ee32d3aa4f7af90e30`，overlay effect SHA `d15acbde9438f5d2391c4da23705a34a3585833062d616517d5af052088606c2`。原 medium R1354=0.01607 先由 profile 验证，再按既定 runtime override 得到1000；不得换培养基。来源 `B/replay_input_inventory.json:effective_context/overlay_audit`、`T/inputs/po1f_sd_leu.json`（supported）。

finite_batch 全14库存（mmol/L）：

| 反应 | 初值 | 反应 | 初值 |
|---|---:|---|---:|
| R1070 | 111.0 | R1003 | 0.054 |
| R1202 | 0.287 | R1204 | 0.601 |
| R1215 | 0.095 | R1217 | 0.381 |
| R1220 | 0.274 | R1222 | 0.134 |
| R1223 | 0.303 | R1231 | 0.840 |
| R1232 | 0.245 | R1233 | 0.276 |
| R1234 | 1.195 | R1354 | 0.178428 |

nonlimiting 仅移除 R1354 有限库存，其 uracil start/end 继续为结构性 NaN，不能当零。不能用 profile 中不同精度的库存数代换，不能从 glucose/uracil 两列猜另外12库存。来源 `F/scripts/gem_annotate/quinone_dfba_essentiality.py:116` 与 `B/reference/replay_input_inventory.json:authorized_baseline/all_recorded_mode_settings`（supported）。

保持原 pFBA 与既有 `OptimizationError`→原 `model.optimize()` 后备；source-free growth >1e-9 h⁻¹ 时不使用 reserve；source-free 非 optimal 立即返回，不新增 reserve rescue。保留 `amount/(B_start×dt)` medium 限制、显式 Euler、库存非负更新、Q9 原固定 clamp、原24 h/提前终止逻辑。记录门控异常不能被转换为 `OptimizationError` 而触发额外后备。来源冻结 `_finite_medium:429`、`_optimize_minimal_pool:438`、`simulate_gene:538`（supported）。

## 5. 实际对象零求解预检与参数记录（R1、R4）

记录设计依同目录 `PARAMETER_CAPTURE_FIX_PROPOSAL.md`，只保证捕获 API 返回的受支持标量类型与值/位型，不保证隐藏基、自动算法选择或求解器内部原始表示。只记录 API 所得，不凭名称把有限值改成 Infinity，不把旧 null 按默认值回填。

新 `resolved manifest` 不能凭文档凭空产生。新增 **`real_object_parameter_preflight`** 待批准项，候选范围为：读取冻结输入和获准 sidecar，创建实际对象，读取实际参数/预期及发现集合/排除项/依赖身份并生成 manifest，`optimization_calls=0`。不含 `optimize`、`pfba`、solver smoke、新条件或状态推进。本轮没有导入或查询真实求解器。

候选预检采用独立进程中的独立对象，完成记录后释放；它不是 B1/B2/B3 后续对象。未来每个新 worker、mode/context/native 对象创建后必须重新读回，与获批期望 manifest 核对并产生新的捕获事件；预检值不能直接当后续对象的实测值。不能匹配时在第一次优化前 STOP，不猜测补齐、不自动设值使其匹配。

独立预检候选预算：一个独立进程，至多为两个 mode 各建一个对象，串行读取；完整 active wallclock 至多120秒，含加载、读取、写manifest与收尾；进程及子进程内存总限4 GiB，全部新增manifest/暂存输出总限16 MiB。它不计入3/8条重放，但必须随 `real_object_parameter_preflight` 单独批准，不是免费夹带操作；`optimization_calls=0`。环境若无法提供获准外层限额/收据，在操作前报告，不以导入或查询solver试探。预算中断策略需要该独立预检授权明确接受，不能借用未来 batch release；不设任何solver资源参数。

| 记录项 | 旧记录等级与候选约束 |
|---|---|
| Python3.11.15、cobra0.32.1、Gurobi/gurobipy10.0.3 | supported：原环境记录；未来读实际版本，变化停止，不自动升级 |
| optlang1.9.1、numpy2.4.6、scipy1.17.1、pandas2.3.3、python-libsbml5.21.1 | supported：WP2 v1 环境记录；不是 r2 已知环境的回填 |
| Threads1、Seed0、Method1、FeasibilityTol1e-9 | supported：v1 显式设置；候选保持原设置并逐阶段读回 |
| OptimalityTol1e-7、Presolve0 | supported：v1 首解前读回；候选只读一致性闸门，不新增 setParam |
| 其余参数和未知项 | 历史 unknown 继续 unknown；新实际值逐项、类型保真捕获；14 个旧 null 不因此恢复 |

来源 `B/historical_vs_current_environment.json`、`B/current_environment_recording_caveat.json`、`B/analysis_code/coq9_wp2_replay.py:503`。预期参数集合不机械锁死151，按 sibling 提案记录 expected/discovered/excluded/difference；不足时失败闭合。

`codec_contract`、`native_hook_coverage`、`runtime_noninterference` 分开验收。合成 codec/fake 测试不证明真实 native 覆盖。独立真实 hook smoke 若需要，须另列具体操作、调用次数、预算并另行授权；不能藏入合成测试或零求解预检。本计划不安排额外 smoke。后续获准有限运行可以明确包含对其原定调用的真实 hook 观察，不增加条件、调用或后备；未观察到的分支保持 coverage unverified。

候选运行声明仍为 `purpose=baseline_equivalence_check`、`historical_environment_completeness=partial`、`exact_historical_environment_reproduction=not_established`、`state_origin=new_run_with_declared_settings`。即使数值一致也不升级历史 unknown。

## 6. 四个比较层、参考冲突和 WT 分母（R5）

### 6.1 参考分层：不选择较易通过者

| 层 | 参考／比较对象 | 可比较范围与不可用项 |
|---|---|---|
| L1 `r2_archive_comparison` | r2 原 intervals/calls 与新运行 | 原保存的时间、growth、biomass、Q9、glucose/uracil、面板 flux、calls；其他12种逐步 pool、历史全 flux、native LP/次级目标未归档，记 comparison_unavailable |
| L2 `wp2_v1_prefix_comparison` | v1 实际执行前缀与新运行 | 同条件、step、phase、目标语义且双方有记录者可比；v1未执行545行尾段无此参考，不补造；v1 pFBA 返回解不自动等于 v2 每个 native primary solve |
| L3 `v2_internal_consistency` | 新状态、调用、轨迹和账本自身 | 全库存 Euler/连续性、Q9原clamp、time、coupling、Sv、bounds、rollback、oxygen/ATP账本算术、参数/对象捕获覆盖；只证明内部检查，不证明原生机制或旧隐藏状态一致 |
| L4 `new_WT_normalization` | 同一获准新计划内、同 mode 的新 WT 与 KO | 按下节判定可用性和ratio；不借旧WT、不跨mode、不把不同终止时刻伪装为共同时间窗 |

每个比较记录 `reference_layer`、参考 run/file SHA、condition_key、step、phase、objective_definition、availability、gate_status 和原因。r2 与 v1 在可比同字段互相冲突时，分别输出两边值及其 SHA，标 `reference_conflict_unresolved` 并停止受影响重放/结论；不择易、不平均、不自动指定较新者为真值。独立只读整理可继续。

v1 的返回 pFBA 结果、fallback 返回结果与内部 native primary/secondary 解是不同阶段。仅在双方 phase 与目标一致且有记录时比较；native primary 在 v1 无对应记录则 `comparison_unavailable`，不能拿 pFBA 返回 flux 的 biomass 坐标冒充该 native primary 的完整解。源自由状态/分支若只有外层保存的记录，按该外层语义比較，不扩大为所有内层求解状态已一致。

### 6.2 连续容差保持不变

候选沿用 `abs(new-ref) <= atol + rtol×abs(ref)`，不改为对称尺度，不放宽数值容差。来源 `B/PREDECLARED_PROTOCOL.json:comparison_tolerances`、`B/analysis_code/coq9_wp2_replay.py:99`（supported）。

| 量 | atol | rtol | 单位／范围 |
|---|---:|---:|---|
| 时间 | 0 | 0 | h，含终止与Q9原耗尽时间 |
| biomass | 1e-10 | 1e-9 | gDW/L |
| 有限库存 | 1e-10 | 1e-9 | mmol/L |
| Q9 reserve | 1e-12 | 1e-9 | mmol/L |
| growth/source-free growth/对应主目标 | 1e-8 | 1e-8 | h⁻¹；同语义阶段 |
| flux | 1e-8 | 1e-8 | mmol/gDW/h |
| dynamic doublings/ratio | 1e-8 | 1e-8 | 无量纲 |
| pFBA secondary objective | 1e-6 | 1e-8 | 仅同目标、确有对应记录时 |
| dilution coupling residual | 1e-9 | 0 | 按原方程单位 |
| Sv/bounds residual | 1e-7 | 0 | 按原字段类型，不能称全满足1e-9 |

step、advanced、状态、source-free分支、growth<=1e-9事件与终止事件保持原离散闸门；负库存/非有限值按 sibling 逐字段规则处理。内部 flux 不同但目标、状态、事件与数值残差均通过时，保留差异为 `possible_alternative_optimum_not_confirmed`；不能仅凭坐标差声称机制改变，不能用该标签掩盖驱动库存的 uptake/source/dilution、目标、状态或事件超限。只有单次 pFBA 解，不据此声称 AOX 补偿或某反应必需。

### 6.3 ratio 元数据和可用性

新 KO/WT 仅使用同一获准 `plan_id`、同mode、同冻结输入和获批设置、已通过适用闸门的新 WT；可以来自同计划先前获准批次的独立进程，但不能拼接旧 r2/v1 WT。WT 与 KO 均需达到各自原逻辑的终止点且没有比较/资源 STOP；WT分母必须有限且严格大于0。相同提前 infeasible 终止可成为受限对照，但必须保留状态和时间，不是“完成24 h”。

每个 ratio 必须记录：

- `ratio_definition=KO_endpoint_dynamic_doublings/WT_endpoint_dynamic_doublings`，其中各自 doublings=`log2(final_biomass/initial_biomass)`；这是各自终点倍增数比，不是初始速率比、最终biomass比或共同观察时间窗指标。
- `plan_id`、KO `run_id`/`condition_key`、`KO_termination_time_h`、`KO_termination_status`。
- `denominator_run_id`、`denominator_condition_key`、`WT_termination_time_h`、`WT_termination_status`、分子分母原值与单位。
- `gate_status`、`KO_gate_status`、`denominator_gate_status`、`ratio_available`、`unavailability_reason`；条件键至少含 mode/gene/alpha/pool/dt/B0/hours，另绑定输入和设置指纹。

两者各自提前终止时间不同也不能改写定义；分母不合格、缺新 WT、KO 有 STOP、资源中断等均输出 NA＋原因。原1/5/10/15%严格 `<` cutoff 仅在 ratio 合格时保留用于旧指标比较，不把新连续指标自动套用这些 cutoff，不按实验 recall 调整。3-only 的 nonlimiting KO 无新分母，明确 NA；不能因此偷偷加第四条件。

## 7. 零事件和 terminal 结果表达（R2、R3）

完整规则见 `ZERO_EVENT_POLICY_V2_PROPOSAL.md`；本文不另定阈值。候选 `epsilon_C=min(原库存 atol,1e-8 mmol/gDW/h×B_start×dt)` 只用于有效推进区间的描述和 raw-zero 不一致受限豁免，不回写库存/bounds。两边都正且连续比较通过时，near-zero标签不同只描述，不新设硬门；status/branch等行为门仍独立。

terminal 不计算或继承 band、适用 raw-zero 差异仍 strict STOP，是保守政策选择，不是 dt=0 推出的生物学结论。输出 `terminal_raw_zero_only_mismatch` 时分别记 `solver_termination_observed`、`solver_status_match`、`termination_time_match`、`continuous_state_match`、`raw_zero_match`、`overall_replay_gate`。求解器终止已观察且状态/时间/连续量匹配仍可同时 overall STOP，不含糊写“求解失败”或“死亡”。

`R` 的 R3指出：历史 nonlimiting YALI1A21711g 的 step85 是5.3125 h非推进 infeasible，glucose start/end 为 `4.440892098500626e-16`。若新轨迹从 step84 的0库存进入同一终止，仍可能在 step85 STOP。该归档描述 supported；新运行会否如此 unverified。不承诺三个旧 STOP 必然消失，不由实施者放宽 terminal 或原 `q9_pool_depleted_h` 时间闸门。

## 8. 统一权限和顺序（R4、R6）

以下状态全部为待批准，不由本次文档修订解锁；exact permission names 与 `DECISION_ITEMS.md`、manifest 一致。

| 权限 | 可批准的最小范围 | 明确不包含 |
|---|---|---|
| `parameter_capture_code_implementation` | 获准副本的参数记录最小修改 | 测试、真实对象操作、比较代码、计算 |
| `comparison_policy_acceptance` | 接受版本化比较公式/字段/terminal/stop政策 | 实施代码、测试或运行 |
| `comparison_code_implementation` | 按已接受政策修改获准 comparator 副本 | 政策自定、测试、计算 |
| `synthetic_test_execution` | codec/fake与合成状态检查，范围单列 | 导入/查询真实solver、native smoke、FBA |
| `real_object_parameter_preflight` | §5独立对象零求解预检，生成新resolved manifest | optimize/pfba/smoke、新轨迹；不以0求解名义隐藏真实调用 |
| `new_computation_scope` | 选择3-only或原8条，固定inputs/settings、生命周期、预算、停止与输出 | 自动释放任一批、扩大网格、反事实 |
| `batch_execution_release` | 在前置门满足后明确批准某一批及外层中断策略；可明确包含原调用的native hook观察 | 自动下一批、额外smoke/优化、自动重试或下一工作包 |

顺序：先分别审查政策与代码授权；获批实施后，仅另获授权时运行合成测试；实际对象预检另获批准且有独立预算后才读取真实值；用户接受 resolved manifest、代码身份和比较政策，决定3或8条计算范围；最后逐批释放。若预检需候选代码而该代码尚未获批实施，先停，不以“验证提案”名义实施。

“无副作用”仅拟保证记录不改变数学输入、参数、原定调用次数及控制流，不承诺运行时间/内存/磁盘不变。codec合成通过不证明native覆盖或runtime不干扰；真实运行证据仍需未来明确释放的原定调用观察，未覆盖项保留 unverified。

## 9. 预算、外层监控与不完整捕获（R6）

以下是 **unverified／待批准候选上限**，本轮不运行任何监控/worker/HPCC。每条最多384原循环；3-only最多3条，八条最多8条，不加条件补失败。

| 资源 | 候选定义与监控策略 |
|---|---|
| 每条件30分钟 | 外层单调墙钟，包含该条件的加载、对象构建/读回、capture、求解、写盘和压缩，不只native solve；共享mode初始化计入该mode首条，最终收尾压缩计入所涉最后一条，不允许隐藏无归属开销 |
| 所有批次累计2小时 | 各worker从启动至退出的active wallclock之和，含加载/记录/压缩；不含排队和人工审批等待；按worker计时避免与逐条件时间双计，不随换进程/批次重置 |
| 4 GiB内存 | worker及全部子进程作为一个受控进程组，采用可用系统级硬限额；外层监控记录组内总量和峰值。若只能采样，声明采样周期/误差且不能声称严格硬限额已保证；没有可获准的限制机制则提交环境缺口，不运行 |
| 2 GiB新增磁盘 | 所有新输出、暂存、未压缩和压缩文件在明确新输出/暂存根下合计；冻结只读输入不计，不能把暂存移出根规避。优先用独立受限卷/配额硬上限并外层记账；只采样无法保证硬上限则先报告待决，不边跑边放宽 |

外层资源监督器不设置 `TimeLimit`、`WorkLimit`、`MemLimit`、算法或任何solver参数。若用户未来 `batch_execution_release` 明确批准硬限额中断，进行中的 native 调用超限可以由外层终止整个计算进程组；不是修改冻结科学算法的分支，也不是计算已正常返回。不得为预算跳过完整捕获、换算法、自动延时、重复条件或重试。监督器若需要新代码，仅可在对应后续实施范围另行批准，本文不实现。

独立零求解预检和独立真实smoke（若用户之后决定需要）必须分别获授权与预算；不从上述3/8轨迹预算夹带，也不能借增加“预检条件”获得新计算。

未来应在每次native调用前封存 `call_started` 事件和已成功捕获的pre记录，包含 PID/worker/native ID、阶段、时间与文件SHA；返回并完成post捕获后再记 `call_completed`。如果限额或外部中断：

- 可辨认进行中调用时记 `execution_status=interrupted_by_external_budget`、`capture_status=incomplete_post_missing`、`native_call_outcome=unknown_due_to_interruption`、`overall_replay_gate=STOP`；PID或对象身份未保存则另外 `object_identity_capture=unavailable`，不能推断。
- 捕获自身失败但调用完成时，分别保存已知 solver outcome 与 `capture_status=incomplete_capture_error`，门控失败不得触发额外fallback；不能伪造post值或复用上一阶段参数。
- 硬杀/OOM无法写最后一条时由外层收据记录进程组、资源事件、最后已提交记录ID和缺失范围；原工作文件保留为partial，不填作完整轨迹。时间/终止状态标资源中断，不换算成24 h生物学死亡。

若环境不能保证外层收据/限额而设计要求它们，停止在放行前列出差距；不调整算法来迁就环境。预算耗尽即停止，独立只读证据整理可继续。

## 10. 未来交付和本轮停止

仅当分别获批后，未来新run目录可输出：预声明协议、输入/代码清单、参数集合与类型保真双哈希、对象阶段捕获、环境未知范围、轨迹/state/attempt/full flux/LP及账本、L1–L4分层比较、ratio元数据、参考冲突、资源/终止/捕获缺失收据、旧STOP对照和SHA清单。2305参考行不预填成新运行完成数；545旧未执行行不离线补成PASS。

仍未解决：14个历史原值；未获准的实际对象预检和resolved manifest；native覆盖与runtime非干扰尚无新执行证据；3/8目标与逐批授权；terminal strict政策可能产生新STOP；旧隐藏状态/其他库存/全量向量缺档；原生身份与呼吸机制不确定性。以上不阻止提交本文，却限制相应复现或生物学声明。

**本轮止于 proposal_v2 文档交付：不实施代码、不运行测试、不导入或查询真实求解器、不提交 HPCC、不执行3/8条/dt矩阵/FVA/反事实、不修改模型或旧STOP、不自动推送main。任何后续权限均待用户明确决定。**
