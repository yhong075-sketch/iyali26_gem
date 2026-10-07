# CoQ9 WP2 参数捕获最小修订提案

## 0. 材料护照与边界

| 项目 | 内容 |
|---|---|
| Origin | experiment-agent plan |
| VerificationStatus | UNVERIFIED design；proposed_not_implemented |
| 起始日期／文档修订日期 | 2026-09-05／2026-09-06 |
| Version Label | proposal_v2 |
| document_schema_version | 2（文档元数据；不改变 trajectory 1.8 或候选参数 schema 2） |
| 目的 | 提出下一版无损参数记录及实际求解对象／时点绑定方案，不实施 |
| 当前授权 | R1–R6 文档定点修订；代码实施、合成测试、真实对象零求解预检与新计算均未获本轮授权 |

本护照描述**提案设计**，不把已核对的输入身份降为未知，也不把设计写成已验证的实现。按科研治理技能分开事实、解释、设计和人类审批；按 Ponytail 复用现有 observer/recorder，不增加通用框架或新依赖。本轮没有导入／执行封存求解代码、查询现有求解器对象、运行测试、提交集群作业或修改科学输入。

来源别名（所有定位均指原字节归档内成员，不是当前工作树同名文件）：

- `H`：外层 ZIP 的 `coq9_wp2_review_and_next_step_handoff_20260905/`。
- `A`：`H/coq9_wp2_independent_audit/`。
- `B`：`H/original_archives/coq9_wp2_qualified_baseline_20260905.zip::coq9_wp2_qualified_baseline_20260905/`。
- `T`：`B/execution_provenance/coq9_wp2_qualified_20260905T042505Z_inputs.tar.gz::coq9_wp2_qualified_20260905T042505Z/`。
- `R`：`B/analysis_code/coq9_wp2_replay.py`；`C`：`B/analysis_code/coq9_wp2_capture.py`；`F`：`T/frozen_compute/scripts/gem_annotate/quinone_dfba_essentiality.py`。
- `E`：`B/historical_vs_current_environment.json`。

以下为保留的 **proposal_v1 previous source audit** 身份记录；本参数文档作者未重新展开归档或重算这些成员身份，不将旧核查写成此作者此轮复核；整个任务的实际来源复核范围另见新 manifest：外层 ZIP `3d732a67cf3f90a9d26213cf2a160db491406f7492d855cf1ce06c42acef9001`；内层原始 ZIP `0b184c04a3a1270f413c4a5fbba60d386eedf9ea2560b1b5796a1fe905d90e16`；执行输入 tar `38fdfd742f326365355d0ee0f1924d3b2bb69f00ebd9f0e8a305c60abf260bff`；R `7b75c75344613ed665fc355e16eff0bd67447fbf89a514ce8ab79ace6cd1a166`；C `3db5da2ed71994990395d0e7a738ce7d3dab35473576aabb7ba7c52abcd552fa`；E `fc53c0edb505bfe82551f3c24459f91b839446fa3dad9bde40087ab095672552`。旧记录描述字节身份核查，不是重放复现。新文件身份以新的 document_manifest.json 与 SHA256SUMS 实算结果为准。CoQ9 本次意见基于转交全文与其已有来源，未收到新提案 ZIP，也未独立核验六份新文件或新 ZIP 的字节 SHA（用户附件“审查结论与权限”）。

科学身份不变：模型 `bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee`、冻结计算提交 `36bb6f0735e4c6458bd53c0ceb01952b116b8be7`、冻结 runner `679ada071adb8e1c898ed1b60d9d8c3f895dac5ca56e1e4c19b1ea28ae6eb5eb`。依据：R:18–24、451–452、480–485；本提案不覆盖原文件、不改变其身份，也不声称完整历史 dirty 环境已恢复。

## 1. 已核查事实与结论等级

下表及第2–3节保留 proposal_v1 的 previous source audit 判定与 locator；本参数文档作者未重做其完整归档源码或数值审计。旧 supported 表示此前所列材料支持该受限事实，不表示 CoQ9 聊天或本轮重新验证了其全部字节和结果。

| ID | 原子判断 | 状态 | 来源定位与限制 |
|---|---|---|---|
| P01 | 旧 WP2 记录有151项过滤后参数，其中14项 `current_value=null` | supported | E `/historical_vs_current_parameters`；A/parameter_and_preparation_check.json `/current_parameters_count`、`/null_current_parameters`；previous source audit 核对；151仅为旧观察数量 |
| P02 | 捕获调用读取实际 `model.solver.problem.getParamInfo(name)[2]`，而不是默认Env | supported | R:138–147、497–506；这是源码对象路径，不是逐次求解对象身份的完整运行证书 |
| P03 | 非有限float经clean变None，写盘与fingerprint共用该有损编码 | supported | R:52–84、514–526；A/SOURCE_CODE_EXCERPTS.txt:3–35、71–80 |
| P04 | “已有当前参数记录无损且完整” | contradicted | P01、P03；B/current_environment_recording_caveat.json `/current_parameter_serialization`、`/affected_count` |
| P05 | “记录层把null重新配置进求解器，已证明参数污染” | unverified | R:43、503–508只从CURRENT_SETTINGS设置四项；已读R、C、F未找到该回灌路径。静态未发现不是对所有运行语义的证明；A/evidence_status.tsv对应 `serialization_defect_directly_reconfigured_solver` |
| P06 | “现有相同fingerprint证明14字段各阶段未变化” | contradicted | R:74–75、259、542、569–570均依赖有损编码；同一null可来自不同非有限值。它仍可约束已保存的有损投影，不可升级为全值不变证据 |
| P07 | 14项原始符号／数值能从已检材料逐项恢复 | unverified | 第3节恢复搜索；未找到合格旁证。不是宣称任何地方永久不可恢复 |
| P08 | 记录缺陷影响全部8条件，不只3条STOP | supported | R:495–542，E `/solver_settings_sha256`；A/INDEPENDENT_REVIEW_zh.md:72–85。首模式写一份environment，其他模式及全部condition复用／比较同一有损参数指纹 |
| P09 | “现有observer已完整覆盖pFBA内部所有原生后端求解的参数时点” | contradicted | C:40–44、88–118、165–175能观察直接pFBA/fallback返回及pFBA临时LP的get_solution返回；R:283–309未捕获逐阶段参数。至于能否在原依赖上实现完整新挂钩，因包内缺少依赖库完整源码而另属unverified |

本表保留 previous source audit 的9条来源判定；supported 4、partial 0、unverified 2、contradicted 3；其中设计有效性不在已核实事实中。没有重新计算独立审查的全量状态／通量数值结论。第3节对state抽查明确另标partial，不能按本表条数冒充全量状态审计。

## 2. 14字段、对象与序列化路径

下列JSON Pointer均相对于 **E**，并非仅列参数名。每项原文件值为null，历史字段为字符串 `unknown`；14项无损恢复结论均为 **unverified**。逐项来源另见 `A/null_parameter_inventory.tsv` 第2–15行。

| 参数 | 原字段的完整JSON Pointer | inventory行 |
|---|---|---:|
| BestBdStop | `/historical_vs_current_parameters/BestBdStop/current_value` | 2 |
| BestObjStop | `/historical_vs_current_parameters/BestObjStop/current_value` | 3 |
| Cutoff | `/historical_vs_current_parameters/Cutoff/current_value` | 4 |
| ImproveStartNodes | `/historical_vs_current_parameters/ImproveStartNodes/current_value` | 5 |
| ImproveStartTime | `/historical_vs_current_parameters/ImproveStartTime/current_value` | 6 |
| IterationLimit | `/historical_vs_current_parameters/IterationLimit/current_value` | 7 |
| MemLimit | `/historical_vs_current_parameters/MemLimit/current_value` | 8 |
| NodeLimit | `/historical_vs_current_parameters/NodeLimit/current_value` | 9 |
| NodefileStart | `/historical_vs_current_parameters/NodefileStart/current_value` | 10 |
| PoolGap | `/historical_vs_current_parameters/PoolGap/current_value` | 11 |
| PoolGapAbs | `/historical_vs_current_parameters/PoolGapAbs/current_value` | 12 |
| SoftMemLimit | `/historical_vs_current_parameters/SoftMemLimit/current_value` | 13 |
| TimeLimit | `/historical_vs_current_parameters/TimeLimit/current_value` | 14 |
| WorkLimit | `/historical_vs_current_parameters/WorkLimit/current_value` | 15 |

**共同实际路径（supported，限定静态源码）：**

1. R:497–500每种模式创建context，取 `context.model`；R:503调用F:485–491设置solver与FeasibilityTol；R:504–505对同一 `model.solver.problem` 设置 `Threads=1, Seed=0, Method=1, FeasibilityTol=1e-9`。
2. R:506调用 `solver_parameters(model.solver.problem)`；R:140–147从Gurobi公开参数常量枚举名称，取返回tuple第2索引作为有效当前值。该函数排除凭据／服务相关参数；因此151项应称“过滤后的执行参数集”，不是无边界的全部配置／环境。
3. R:514–519把当前value放入E的相应字段；R:525用 `fingerprint(settings)`；R:526经 `dump → encoded → clean` 写E。非有限float在R:53–54被统一变成None，JSON里为null。`getParamInfo`返回None或抛错的旧路径没有逐项错误记录机制；现有null因此不能在没有旁证时细分原值的符号／NaN。
4. R:533–534在另一模式比较有损参数fingerprint；R:542赋给recorder；R:259各state仅保存 `solver_settings_sha256`；R:569–570在condition结束后再次做有损比较。没有逐阶段、逐次原生求解前后的完整参数值。
5. R:355–357的state fingerprint涵盖这些引用，但不能恢复引用背后已丢失的类型／符号。原有 `234503785550275d0d77e0e1be622a357d7a11594cc188643998c028e4ac8601` 保持其旧语义；不得改写成新无损指纹。

参数记录是观测支路，不是控制输入支路：R:503–508从显式设置及实际对象读值，R:514以后写证据，未见将序列化null送回setParam。故“记录缺陷”已支持，“已证明求解器参数污染”仍未验证。也不从结果相符反推资源／迭代／cutoff参数一定不起作用。

## 3. 原值恢复：保留的 proposal_v1 搜索范围与结果

恢复门槛（proposed_not_implemented）：旁证必须绑定同一原WP2执行、相应实际solver model、捕获时间／阶段、参数名及类型；仅版本默认值、另一模型、另一次查询或同作业的独立预检进程不合格。若以后找到合格原始旁证，新增恢复注记并保留原null，不能覆盖E。

以下是 **previous source audit** 的已报告只读检查范围，不是本轮重复执行记录：

| 范围 | previous source audit 观察 | 判定 |
|---|---|---|
| E及B/current_environment_recording_caveat.json | E保留14null；caveat只列受损名称，不保留正负无穷或NaN原值 | supported；不能恢复 |
| B/execution_provenance/qualified8.28146735.out 全27,667行 | 对14个精确参数名逐项检索，均0命中；未输出凭据内容 | supported；仅此日志未找到旁证，不把“无命中”当作值的证据 |
| B内 `.json/.md/.py/.txt/.out/.sbatch/.tsv` 文本成员的14名称检索 | 命中文件为E、caveat、execution_provenance下final/initial/input三份review、preflight metadata、EXECUTION_NOTES，共7份；逐一查阅相关参数节点 | supported；没有找到原值备份 |
| B/execution_provenance/independent_final_review.json `/environment_caveat`；independent_initial_state_review.json `/environment_review` | 记录null清单及缺陷解释，不是原值备份 | supported；不能恢复 |
| B/execution_provenance/preflight_no_optimize/metadata.json `/settings` | 151项中同样14null；`/optimization_calls=0` | supported；另一个对象、同样有损，不能替代实际重放值 |
| T/runner/preflight_no_optimize.py:15–32及T/execute_qualified8.sbatch:16–17 | 预检自行加载context/model，先由一个Python进程运行；重放随后exec另一个Python进程。预检不等于后续参与求解的模型实例 | supported；禁止跨对象恢复 |
| T中 `.json/.md/.py/.txt/.out/.sbatch/.tsv/.toml/.cfg` 文本成员 | 名称命中仅独立预检报告及若干别的历史脚本／报告；读取R、C、F实际相关路径；未见保存这14值的同次无损转储 | supported，限定检索覆盖；其他实验／默认设置不作为恢复材料 |
| B文件名清单 | 未发现 `.prm` 或 `.rec` 成员；这两类若未来存在也须验证同一对象、阶段及完整性，扩展名本身不足 | supported，限定包内清单 |
| state证据 | proposal_v1 审计抽读1份state的参数引用，并追踪R:236–260、347–357写法；A/INDEPENDENT_REVIEW_zh.md:79、90报告全部state引用有损fingerprint | partial：该次未逐一重读全部1760份state，也未重新审计全部LP／通量 |

proposal_v1 审计没有检查远端存储、原进程内存、未提供日志、其他工作树当前solver或最新默认Env；本轮也未补查，不将这些未知范围写成已搜索。库内部cobra/optlang实现并未作为可定位的依赖源码随T提供；没有用本机安装版本替代原执行依赖源码。

结论：**14项分别仍为unverified／未恢复**；历史 `historical_value=unknown` 与WP2的 `current_value=null` 是不同时间层的缺口，不能互相填补。只补三条STOP不会补全另外五条旧运行的完整参数记录（A/INDEPENDENT_REVIEW_zh.md:81–83）。

## 4. 候选无损表示（全部proposed_not_implemented）

**“无损”的有限保证范围：** 只保证实际捕获 API 在该对象／时点返回的受支持标量类型，以及返回对象的值／可取得位型。它不是求解器内部原始表示、自动算法选择、解析后的内部配置、basis、warm start、数值状态或其他隐藏状态的完整快照；也不证明读取前 API 从未转换过内部表示。API 返回 finite sentinel（例如有限大数）时，必须按实际类型与 finite 位型保存，不能根据 `TimeLimit` 等参数名称或预期语义改成 Infinity。参数语义解释另列并绑定来源；旧 null 不按默认值回填。

使用一个参数专用带标签标量协议；不修改全局 `clean()` 对旧state、LP、轨迹的既有行为。候选命名空间 `coq9.wp2.solver_parameters/v2`，`parameter_capture_schema_version=2`；与冻结dFBA `SCHEMA_VERSION="1.8"` 完全分离（F:45；R:451）。

每项参数记录分为 `capture_state`（是否成功读到）、`value`（类型／精确表示）、`semantic_acceptance`（该阶段能否接受），不能仅用null。以下是**数据契约示例，不是实现或测试结果**：

| 情形 | 候选表示／约束 |
|---|---|
| 有限binary64 | `capture_state="observed"`；`value.type="float64"`；`value.class="finite"`；`value.bits_be_hex` 为该 API 返回 float 的 IEEE-754 大端64位的16个小写十六进制字符，作为权威值；可加与位型一致的 `hex` 显示字段，不依赖十进制舍入 |
| +Infinity | observed；float64；`class="positive_infinity"`；`bits_be_hex="7ff0000000000000"` |
| -Infinity | observed；float64；`class="negative_infinity"`；`bits_be_hex="fff0000000000000"` |
| NaN | observed；float64；`class="nan"`；保存该 API 返回 float 可取得的64位位型，包含符号／payload；不能统一成字符串"NaN"丢失位型。round-trip以位型而非NaN相等判断 |
| +0.0／-0.0 | 都属finite，位型分别为 `0000000000000000`／`8000000000000000`；不得用数值相等合并 |
| int／bool／str | 分别显式type；int用规范十进制字符串保留任意精度；bool保留布尔类型且先于int识别；str原样UTF-8，字符串"Infinity"不解析成float |
| missing | `capture_state="missing"`，无value键，保存原因及原证据定位；用于预期证据字段缺失／旧值未保留，不能伪装成观测到的数值 |
| read_error | `capture_state="read_error"`，无value键；保存失败操作、参数名、错误类／受控错误码，不倾倒可能含凭据的异常全文。实际API抛错、expected参数返回None／结构错误属此类，不猜值 |
| not_applicable | `capture_state="not_applicable"`，无value键；必须附“版本／后端确无此参数或该记录对象不适用”的规则版本和来源；不能因某MIP参数对当前LP不活跃就省略可读取的实际值 |

仅允许声明的标量类型。未知对象类型、字典／列表参数值、重复参数名、非法位型、class与位型不符一律验证失败，不调用str()兜底，不静默转float。JSON本身不允许裸NaN/Infinity；所有非有限数以标签及位串表示，`allow_nan=false`。

“成功捕获非有限数”不等于“求解语义有效”：`observed + positive_infinity`与`read_error`必须不同；NaN即便成功读到，默认是不可接受参数状态，应在下一原有求解前停止并报告，除非后续经版本文档和人工批准明确例外。正负无穷也须按该参数／版本语义审查，不把其普遍视为错误或普遍视为可接受。未知、读取失败或不适用标签不允许成为setParam输入。

## 5. 实际对象与捕获时点（全部proposed_not_implemented）

### 5.1 对象绑定

每条capture绑定 `run_id`、输入／代码完整SHA、condition ID、step、外层attempt序号、实际后端solve序号、阶段、source-free/reserve分支、UTC时间和进程内单调序号。记录实际COBRA model、optlang model和 `model.solver.problem` 的类／模块及进程内对象编号；对象编号只用于本进程区分，不作为跨运行等价依据。每次从调用现场核对对象链，不能缓存旧problem引用后默认仍有效。

`object_kind="solver_model"` 为必要条件。`defaultEnv`、新建探测Model／Env、预检对象或复制模型可另有来源标签，但不得进入“实际执行参数”字段。保留过滤执行参数名清单、排除策略版本及公开API来源；凭据／license／server等字段在读取与输出前排除，不在指纹里暴露秘密。

### 5.2 参数集合完整性契约

不能把“枚举到的每项都读到了”当作完整。运行前候选记录四个带版本集合，均按名称排序并纳入值哈希。以下 E/D/X/C 仅为本小节的集合记号，不是第0节的归档来源别名：

- `expected_public_names`（E）：由固定 solver／API 版本的公开参数清单与可定位源码／文档建立，保存来源 SHA、版本、获取范围；本轮尚未建立或实测，保持 unverified。不能直接令 E 等于运行发现集。
- `discovered_public_names`（D）：实际依赖枚举发现的公开名称；记录枚举 API、过滤前名称集合与依赖源码身份，不导出凭据值。
- `excluded_names`（X）：运行前审核的 credential／license／server 等不应读取或披露项；保存名称、理由、规则版本。先按已审查名称排除再读取；不能在失败后扩展 X 以求通过。
- `captured_names`（C）：实际对象成功返回且满足受支持标量契约的名称；另存每项错误／不适用记录及 `missing_expected`、`unexpected_discovered`、`excluded_but_read`、`read_failures`。

候选完整门要求非排除预期集与非排除发现集一致，且每个预期可读取项均成功捕获；预先有版本来源支持的确实不适用项须单列，不暗中作为成功值计数。缺失预期名、额外未审查名、重复／别名冲突、未授权排除、读取错误或未知类型均令 `capture_complete=false`，在下一原调用前 STOP。额外名先只记录安全名称并提交审查，不为扩大覆盖而读取可能的秘密。版本变化时重新审核 E/X；**151只是旧观察数，不是永久数量断言**。D=C 而 E 尚未建立时最多为 `discovered_set_read_complete`，不能获完整参数集资格。公开清单未覆盖的隐藏配置继续 unknown。

### 5.3 有效时点与覆盖

F:438–456显示既有流程先以source上界0调用pFBA；遇OptimizationError按原逻辑fallback；source-free为optimal且增长不大于原1e-9阈值时才可能再进入reserve pFBA及其fallback。本提案不改这些判断，不调整bounds、目标、参数、积分或重试策略。

| 捕获点 | 必须记录 | 门控 |
|---|---|---|
| 每种模式配置完成、首次原有求解之前 | 同一实际problem的完整过滤参数集；与明确批准的四项设置逐项比较；模型和依赖身份 | 不允许其他对象的值代替；不一致就停，不自动setParam修正 |
| 每个condition入口及每步on_start | 实际对象链、当前参数证据引用、分支／step绑定 | 原对象意外更换、必需字段未捕获或未声明变化时，在下一求解前停 |
| **每次实际原生后端调用立即之前／正常返回或异常之后** | 同一个实际problem的完整参数快照、native_call_id、外层attempt与phase；返回时状态和异常分别记录 | 覆盖pFBA内部生物目标固定所需primary solve、pFBA secondary solve、原fallback生物目标solve，以及原库如有的内部重复调用；绝不能把一次pFBA调用当成保证只有一次native solve |
| pFBA临时LP上下文内、其get_solution返回时 | 参数证据关联临时目标／LP身份，不能待上下文恢复后冒充secondary阶段 | C:165–175已有临时LP观察点可复用，但不能以这个返回点代替上述原生pre/post |
| 每次外层attempt退出、condition回滚后及run结束 | 参数快照与最后已确认快照的带类型差异；完成／失败／未执行状态 | 未知变化、缺失post、limit/iteration/cutoff异常退出记录并停；不归类生物学死亡 |

阶段名候选：`primary_biological_objective`、`pfba_secondary_objective`、`fallback_biological_objective`。分支名与阶段名分开记录；不能仅凭调用序号猜阶段。若原依赖源码无法确定某native调用的阶段，保留 `phase="unclassified"` 和调用来源，并停止“完整分阶段记录”资格，不把它强行归类。

**实现前置缺口（unverified）：** R/C与F已足以提出以上覆盖要求，但包内无可核验的完整cobra/optlang内部后端调用源，不能在本提案中许诺某具体内部函数／行钩子一定有效。后续代码授权须包括只读检查并固定实际依赖源码SHA／版本，确定能在原native调用现场同步观察pre/post的最小挂钩；候选只在C扩展现有process-local观察，不修改F或第三方库文件。若只能观察pFBA外层返回，则方案不满足验收，须停下来重新提交，不以近似时点冒充完整覆盖。

**候选非干扰边界**只要求不改变数学输入、实际参数、原定求解调用及其控制流；不承诺墙钟时间、CPU、内存、I/O 或调度完全不变。参数读取、序列化、哈希及压缩开销应另记 elapsed、peak RSS、输出量与 incomplete capture。资源预算及处理中调用的停止规则见 WP2_REPLAY_V2_PLAN.md 待批方案；不通过隐含 setParam(TimeLimit 等)、少记阶段、换算法、延时或重试实现预算。

观测回调只读取参数，不调用optimize/pfba/slim_optimize、不创建Model／Env、不setParam、不为观察调用update/reset、不改变函数参数或求解顺序。原native函数仍恰好调用一次；原计算本来存在的fallback不由捕获层新增。门控异常必须与 `OptimizationError` 不同且不被转换，避免冻结逻辑把记录失败误当求解失败而新增fallback。捕获前失败则不进入该求解；捕获后失败则保留已发生的一次调用和返回／异常事实，在下一次求解前停，不重做本次。

参数变更检测比较连续pre/post、同一对象跨调用和condition回滚快照；保存每次发生的精确差异。若库本身有合法的临时参数变换，必须在运行前由已固定源码支持并明确列入允许转换表，不能从新结果倒推豁免。不声明两次采样之间从未发生过任何瞬态变动；声明只覆盖已捕获边界和已审查代码路径。

## 6. fingerprint与兼容（全部proposed_not_implemented）

最小双哈希：

1. `parameter_values_sha256_v2`：SHA-256输入为固定域前缀 `coq9.wp2.solver_parameters/v2` 加NUL分隔符，再加规范UTF-8 JSON。JSON包含schema版本、expected／discovered／excluded／captured名称集合与集合差异、排除策略版本、按名称排序的全部带标签记录（capture_state、type、class、位型／值、受控缺失原因）。固定键排序、紧凑分隔符、严格类型，禁止重复键；不同类型、正负符号、signed-zero、NaN位型、unknown／read_error／not_applicable都不能合并。比较时同时核对规范字节，不把密码学哈希说成数学上的绝无碰撞。
2. `capture_record_sha256_v2`：另加实际对象链、时间／序号、condition/step/attempt/phase、输入／代码来源和上一个capture证据引用。它绑定证据事件，不作为不同时间／不同运行参数相同的判据。

成功值与错误标签在同一个封闭schema下具有不同字段结构；不允许用户字典伪装成标签、不允许字符串"missing"等变成状态。哈希不经旧clean，也不对原float值四舍五入／取绝对值。失败记录也可有证据hash，但须 `capture_complete=false`，绝不能因两个相同错误hash而通过“环境完整”门。时间／对象ID不进入第一个值hash，避免相同值因时钟或内存地址不同而无法比较。

原E、caveat、state、旧 `solver_settings_sha256`、旧state fingerprint和SHA256SUMS全部原样保留；不离线重签旧state或伪造新字段。若日后获准新运行，只在全新输出目录写 `parameter_capture_v2/` 中的参数证据，并让新state增加明确版本的引用；冻结trajectory schema继续1.8。旧读取器只读旧命名空间；新读取器遇旧null时只能显示 `legacy_lossy_unknown` 解释标签及原值null，不能自动迁移为observed Infinity。回退时也不得把v2非有限值转换回null并标作完整。

## 7. 最小拟修改范围与原记录保留

下表是将来需审批的**候选文件范围**，不是本轮写入路径；不复制／覆盖B内封存文件。

| 候选文件（在下一版独立执行副本中） | 最小职责 | 明确不做 |
|---|---|---|
| `runner/coq9_wp2_replay.py`（对应R） | 一个参数专用编码／验证小函数组，v2 hash与证据输出；Recorder绑定capture引用与失败门；继续复用既有文件排他创建策略 | 不改全局clean、不改原比较容差／exact-zero、不改CURRENT_SETTINGS值、不把旧null回填 |
| `runner/coq9_wp2_capture.py`（对应C） | 在现有observe内扩展实际native调用pre/post及临时LP关联，记录真实对象／阶段；finally恢复观察器 | 不改F、不改第三方库、不新增求解／重试，不造新观测框架 |
| 现有 `test_coq9_wp2_replay.py` 与 `test_coq9_wp2_capture.py` 的下一版副本 | 测试代码写入需列入对应代码实施范围；执行还须独立 `synthetic_test_execution` | 本轮不创建任何可运行测试、不运行旧34项或任何新测试 |

**实际对象零求解预检须单独授权 `real_object_parameter_preflight`。** 候选范围仅加载冻结输入和获准 sidecar、创建实际对象、读取实际参数／依赖身份、生成 resolved manifest，`optimization_calls=0`；不含 optimize/pfba、native-hook smoke 或新条件。既有预检记录只作历史辅助证据。若复用 `runner/preflight_no_optimize.py`，其最小格式修订须列入 `parameter_capture_code_implementation` 文件范围；代码许可本身不授权预检运行。

预检 manifest 必须说明对象是否会被后续同一进程保留使用。候选重放按批次独立进程创建对象，故预检对象不等于后续运行对象；每次重新创建都必须在该对象上重新读取并核对集合、值、依赖和构建身份。参数值哈希相同也不能证明隐藏状态或整个对象相同。真实 native-hook 运行观察只有在后续明确批准的有限运行内逐阶段记录时才可开展；额外真实 smoke 必须另列具体范围获批，不得藏入合成测试或代码验收。

本子任务仅新增本 Markdown proposal_v2，旧 proposal_v1 保留；未更改执行代码、模型、GPR、bounds、化学、培养条件、curation、dossier、alpha/pool、reserve policy、评价标签或输入manifest。原3/2/3状态、三个STOP和545行未执行尾段均不因参数设计而改变（H/README_zh.md:42–49；A/INDEPENDENT_REVIEW_zh.md:35–50）。

## 8. 合成测试设计与验收门（未实施、未执行）

所有场景仅用人工构造值／fake对象／fake调用事件，不导入cobra、gurobipy或冻结runner，不申请license，不读真实模型，不使用实际八条件数值调规则。

| 测试设计 | 合成输入／错误路径 | 最小通过条件 |
|---|---|---|
| 精确round-trip | 正／负有限数、最小subnormal、最大finite、+0.0/-0.0、+Inf/-Inf、多种NaN符号／payload、int、bool、字符串"Infinity" | decode后类型一致；float原始64位一致；非有限与缺失分离；无裸JSON NaN |
| 类型／哈希区分 | int 1、float 1.0、bool true、字符串"1"、+0/-0、正负Infinity、不同NaN位型、missing、read_error、not_applicable | 规范字节不同，对应测试样例hash不同；映射插入顺序不改变值hash |
| 失败闭合 | API抛错、返回None／短tuple、未知类型、重复名／重复JSON键、非法位型、class不符、脱敏错误信息 | 留受控错误证据，无默认值；`capture_complete=false`；不得把失败当observed非有限 |
| 对象识别 | fake actual-model与fake Env／第二model提供不同值；同一名称对象被替换 | 只接受实际调用现场model链；替换被发现；跨对象值不能通过 |
| 阶段覆盖 | fake原函数依次触发primary／secondary；另一路抛OptimizationError后触发既有fallback；模拟库内部多次native调用 | 每次native均有独立pre/post，fake 原函数调用计数与合成预期一致（不代表真实依赖覆盖已验证），phase/branch不混淆 |
| 失败不增求解 | pre捕获失败、post捕获失败、原native异常、pFBA返回钩缺失、未分类阶段 | gate异常不触发fallback；发生过的调用如实保留；下一次调用不开始；无补偿solve |
| 参数变化 | 调用之间值变化、同次pre/post变化、signed-zero变化、非有限符号变化、失败恢复后变化 | 保存精确差异并按预声明转换表门控；不自动调参使hash匹配 |
| 合成非干扰／恢复 | fake model/LP/源bounds/目标/kwargs前后快照；原函数正常与异常退出 | 参数捕获不改变快照；原函数仅调用一次；observer finally恢复；无update/reset/setParam调用 |
| 兼容与防覆盖 | 人工旧null记录、旧hash、重复目标路径、v2未知schema | 旧字节／旧hash不变；旧null只附未知解释，不能升格；拒绝覆盖与未知schema |

### 8.1 三个互不替代的证据门

| evidence_gate | 拟取得的证据 | 不能推出 | 当前状态与权限 |
|---|---|---|---|
| `codec_contract` | 获准后执行人工标量 round-trip、集合差异、错误标签及 fake 对象／事件检查；验证声明类型和 API 返回值可取得位型不丢失 | 不证明真实 API 内部原表示、真实挂钩可达性或运行非干扰 | unverified／未实施未执行；实现与 `synthetic_test_execution` 分开批准 |
| `native_hook_coverage` | 固定实际依赖源码调用路径审查，加上以后获批有限运行的真实 native pre/post、阶段、实际对象和临时 LP 记录 | fake 序列不能替代；未发生分支不算覆盖；外层 pFBA 返回不代表每次 native primary solve | unverified；零求解预检不能独自满足；须明确纳入获批有限运行，额外 smoke 另批 |
| `runtime_noninterference` | 实际获批运行的数学输入／参数、原调用次序与次数、控制流、状态残差及比较闸门记录；资源影响另披露 | 不证明隐藏状态全不变、时间／资源零开销或参数已校准 | unverified；仅获批有限运行可形成证据；如需第二次无观察对照，须另批，不在本提案隐含范围 |

将来的 fake 路径通过最多支持 `codec_contract` 的相应合成检查；三个门分别记录 supported／partial／unverified 及覆盖范围，不能合并称“完整验收通过”。实际 preflight 只证明那个对象在所读时点的公开参数，不含求解阶段证据。`native_hook_coverage` 与 `runtime_noninterference` 的实际运行部分是获批有限基线所要产生的结果，不能循环要求它们已完整 PASS 才允许产生第一份获批基线。开始前须完成的是对应明确授权、代码／合成准备、静态挂钩路径核查与零求解参数预检；实际阶段覆盖未证实时，仍须逐调用 fail-closed，不能据此自动扩展 smoke 或条件。原14旧值未恢复、旧 null 及旧 STOP 不因任何新门通过而改写。

合成设计补两项：① finite sentinel 保持 finite，不按参数名变 Inf；② E/D/X/C 缺名、额外名、重复、未审查排除、expected 缺来源均不能以计数151或 D=C 放行。这些仍只是测试计划，本轮不创建或执行。

## 9. 待决与停止

以下权限均为 **pending_user_decision**，不能从本轮文档授权或 CoQ9 审查意见推断通过：

| authorization_key | 单独决定的内容 | 不包含 |
|---|---|---|
| `parameter_capture_code_implementation` | 第7节最小捕获代码及获准预检格式修订、实际依赖源码只读审查 | 合成测试执行、真实对象操作或求解 |
| `comparison_policy_acceptance` | ZERO_EVENT_POLICY_V2_PROPOSAL.md 的候选公式、字段范围和保守 STOP 选择 | 比较代码实施 |
| `comparison_code_implementation` | 已接受规则的最小比较器代码与记录字段 | 自动执行测试或新计算 |
| `synthetic_test_execution` | 列明的人工标量／fake 事件检查 | cobra/gurobipy 导入、模型加载、真实预检、solver smoke |
| `real_object_parameter_preflight` | 第7节实际对象零求解预检，`optimization_calls=0` | native 求解阶段覆盖或额外条件 |
| `new_computation_scope` | 另行选择严格3-only或原8条的具体计算范围和证据目标；3-only若补新的nonlimiting WT至少变成4条，须明确另批扩展范围，不能隐含增加 | 下一批自动运行、dt/FVA/反事实 |
| `batch_execution_release` | 已选范围内具名批次和墙钟／资源预算逐批放行 | 扩大条件数、隐藏重试、延时或调整参数 |

接受公式不授权代码；代码许可不自动授权测试；测试通过不自动授权真实对象预检或计算；范围许可不代替具名批次放行。详见 DECISION_ITEMS.md 统一表及 WP2_REPLAY_V2_PLAN.md 生命周期／预算方案。

未解决：14原值无合格同次旁证；完整历史环境仍 unknown/partial；expected 名称集、真实依赖内部 native 挂钩位置及类型编码实现未验证；真实运行非干扰与资源影响未验证。相关事实来源保留第1–3节 locator；设计状态不升级为 supported 的实施结果。

输入指纹变化、找回证据与现存记录冲突、未知值被默认值替代、需要修改冻结算法才能完成挂钩、或需扩展计算范围时，只停受影响事项并报告。资源／迭代／cutoff／读取失败均不能写成生物学死亡；参数记录修复也不能接受CoQ9病例、激活模型候选、替代人工reserve的科学验证或绕过既有人类闸门。
