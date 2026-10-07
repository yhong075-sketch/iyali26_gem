# R153 OR GPR 管线接入独立审计（历史版本，已被新指令取代）

审计人：独立子代理 `/root/audit_r153_r2176`。核验日期：2026-09-15；最终记录时间：2026-09-15 21:12 UTC。

**限定实现审计通过；无必须修改项。该 OR 版本随后被用户“去除 YALIUNK2”的新指令取代，保留为历史交付物，不代表后续单基因版本已经通过审计。** 后续代码/整理记录将变化，本审计结论固定在下列 OR 输出、构建清单及核验时的文件身份上。

本轮只读核查实现、既有来源和构建交付物，仅维护本审计文件；未修改代码、整理数据、状态入口或模型，未新检索，未运行 BLAST、结构预测、模型优化或敲除求解。主代理执行了定向软件测试、一次完整 offline/no-solve 构建及静态边界验证；本审计独立读取这些原始交付物，并自行以标准库比较完整 XML、重新计算文件 SHA。未重复运行主代理的测试、构建或验证脚本。

## 授权、身份与证据层级

本任务收到主代理转交的用户明确指令：“将YALI1D17462g补充进R153的GPR 然后并入管线”。主代理已向用户说明按“补充”保留原成员，目标为 `YALIUNK2 or YALI1D17462g`。该阶段授权范围是此模型规则和构建管线，不是新的蛋白功能结论；不包括改变 R2176、计量、边界、培养、实验标签或启动生长计算。完成核验后，主代理转交用户追加指令“去除YALIUNK2”；该新指令开启另一版本，未追溯改写本 OR 交付物。

- **YALI1D17462g / YALI1_D17462g**：原生正式名称未核实；AOW04049.1 / A0A1D8NEI7，384 aa。既有序列和 AlphaFold 预测支持 A1 天冬氨酸内肽酶功能候选，对精氨代琥珀酸合成酶赋值形成反证。OR 连边是用户指定的暂定模型赋值，未标为已验证合成酶或已实验确认的同工酶。整理记录中的 XP_502815.3 仅核查原模型既有注释，本轮不凭该字段新建跨菌株映射。
- **YALIUNK2**：无已核实名、分子身份或版本化蛋白；仅模型中的合成酶占位。OR 阶段保留该条目不等于认可原生功能。

已重读 `artifacts/r153_r2176_gpr_review_20260911/REPORT.md` / `SOURCE_AUDIT.md` 与 `artifacts/argininosuccinate_identity_20260911/REPORT.md` / `SOURCE_AUDIT.md`。本轮不重做此前 18 项/36 项来源审计，也不将历史覆盖数作为此次实现审计数。既有 AlphaFold 结构和固定参考面板结果不是本轮新预测或功能实验。

## 原子声明核查

“支持”表示在所列限定范围内被来源/交付物支持；“未决”表示生物学问题仍未解决。计数不是功能置信概率。

| ID | 有限声明 | 判定 | 核查来源和范围 |
|---|---|---|---|
| A01 | OR 阶段按当时授权追加成员，保留占位，后续新指令另立版本 | 支持 | 主代理转交的两次明确用户指令；`SCOPE.md`、整理记录授权说明。未把后续删除指令追溯解释成旧 OR 的错误 |
| A02 | 构建原始输入与比较输入分开固定，新模型独立输出 | 支持 | `before.json`、构建清单及原始 XML；原始输入 `data/iyali26.xml`，比较输入为已完成 R_NTP1 修改的 `model_metadata_trna_ntp1_hydrolysis.xml` |
| A03 | 赋值精确限制在白名单目标，所有前态检查在写入前完成 | 支持 | 静读 `_apply_gpr_assignment`；R153 原规则/目标规则、状态、六物种计量与公式/电荷/区室、bounds、既有基因注释、冲突 notes 均先检查。R1025/R1026 既有检查保留 |
| A04 | R153 整理在最终 metadata 字段选择后进入默认完整构建 | 支持 | `main.py`、构建日志和 `reaction_selection.post_selection_r153_gpr_assignment`；记录 `applied`，未依靠手改最终 XML |
| A05 | notes 明确保留相反功能证据，未将用户赋值或 AlphaFold 升级为实验验证 | 支持 | 新 XML notes、整理 JSON 和既有功能报告；标记 `user_authorized_provisional_assignment_with_conflicting_function_evidence`、`required_not_confirmed`。保留靶蛋白对蛋白酶 TM 0.74282、对合成酶 0.26708/0.27714 的反证性质 |
| A06 | 本轮所引用四份历史功能来源身份保持一致 | 支持 | 独立计算整理记录 `evidence_source_sha256` 中 4/4 文件 SHA，均匹配；只复核既有证据，不新增同源搜索或结构计算 |
| A07 | 三个定向回归模块通过，覆盖 OR、冲突拒绝、幂等、导出与旧目标回归 | 支持 | 阅读测试源码、`tests_execution.json`、`tests.log`；3 个 unittest 方法，包含多个子用例，日志 `Ran 3 tests in 38.703s / OK`，进程返回 0。不是全仓测试通过声明 |
| A08 | 四种基因状态按预期传到反应边界 | 支持 | 原始 FBC OR 树、测试源码及 `build_validation.json`；无删除和只删 YALIUNK2 不改变边界；只删 YALI1D17462g 仅关闭 R2176；双删关闭 R153 与 R2176。验证调用静态 `gene.knock_out()`，不求解生长 |
| A09 | 一次限定的完整 offline/no-solve 构建成功，未新增优化或研究计算 | 支持 | `build_execution.json` attempt=1、返回 0、耗时 63.175 秒、600 秒上限；实际 argv 含两限制。静读进程内执行保护及日志/清单，诊断求解未运行，网络富集仅缓存；未声称独立系统调用追踪。日志的 LP 文件读入是模型加载/复制信息，不是优化结果 |
| A10 | 最终模型相对指定比较输入，仅有 R153 GPR 与指定 notes 变化，R2176 保持 | 支持 | 本审计独立解析两份完整原始 XML；仅去除 R153 的 `notes` 和 `geneProductAssociation` 后，两整棵树序列化完全相同；R2176 单个完整元素精确相同。目标 notes 读取核对，原 `PROTEIN_CLASS` 保留；2315 反应、1877 物种、1074 基因 |
| A11 | 已消费输入、代码与数据身份可核，构建过程中源码未变 | 支持 | 独立核对构建清单 35/35 代码 SHA、48/48 数据 SHA，包括新整理 JSON；原始输入、输出 SHA 与实文件一致。清单标记 input/source unchanged；HEAD 和 dirty 状态保存。清单历史 `reference_sha256` 不作为本轮比较输入身份 |
| A12 | 本轮保护范围内旧输入与证据保持，未覆盖旧模型 | 支持 | `before.json` 共 76 项，排除明确允许改的 main.py、patches.py、README、PROJECT_STATE 后，独立核对剩余 72/72 全匹配。此范围不是工作区所有文件未变的声明 |
| U13 | 该蛋白的真实原生功能及其独立催化胞质精氨代琥珀酸合成的能力是否成立 | 未决 | 本轮没有功能实验；既有序列/AlphaFold 指向 A1 蛋白酶候选，对合成酶赋值有强反证。软件通过不能消除这项冲突；不能据此称已查明真实合成酶 |
| U14 | OR 同工酶解释是否生物学正确，以及生长/必需性预测是否改善 | 未决 | YALIUNK2 身份仍未知；缺少两个成员可独立催化同一步骤的原生功能、辅因子和区室证据。本轮无生长求解或新的实验评价，静态关闭反应不等于必死或命中率改善 |

**覆盖 14/14 项：12 项支持、2 项未决；未发现与限定软件交付声明相矛盾的项。** 此处并不表示没有生物学反证：A05/U13 明确保留原有强反证。所有未决均限制功能/表型声明，不阻塞已明确授权的模型赋值交付。

## 固定的交付物身份

核验时 HEAD：`a286d18d14bbcb404e8e00f4691c54e6fc81d1db`，含已保存的 dirty 状态。以下完整身份固定本次 OR 版本；后续单基因版本不能借用本审计结论。

| 角色 / 路径 | SHA256 |
|---|---|
| 原始构建输入 `data/iyali26.xml` | `5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d` |
| 比较输入 `model_metadata_trna_ntp1_hydrolysis.xml` | `9adfb0f6187360770f1ca38ebbdb7a9e68d99705c6868a91a8e7f583cd7e7e36` |
| OR 输出 `model_metadata_trna_r153_gpr.xml` | `45a8e0aab47cc9e05254c688cf40b0d1038c054172849c3e331768df4351b319` |
| OR 构建清单 `model_metadata_trna_r153_gpr.build.json` | `f94a8d270dd019224d9f30c82083e795ce4d695abbf0704ddb34e7745a84bd29` |
| OR 整理记录 `data/reference_build/curation/r153_gpr_assignment.json`（核验时） | `3ee7b332642c7cb9b26063720d234ebb96f2693a69472a69ce2a0f0f5158806b` |
| OR 定向测试 `tests/test_r153_gpr_assignment.py`（核验时） | `3b5adb67cc3f41c158ab84965a76a105abbbfb067fc64c8d89aa5ed622301abe` |
| OR 构建 `scripts/gem_annotate/main.py`（清单锁定） | `b277c1e358d108c9d2da62460fb33f8d3868684e83ccc8dce129fb7dc7ceebb8` |
| OR 整理实现 `scripts/gem_annotate/patches.py`（清单锁定） | `b3ff5804ea42c7977ed4876559cb750805a3208c97ae627447ed9605aff1bdeb` |
| 主代理输出核验 `build_validation.json`（本目录） | `46a9c6223e12db881330668e27aebebe7f48951b6dff4f8d47b03d1f973235e0` |

构建 argv、起止时间、全部 35 代码/48 数据 SHA、运行时无菌株 overlay、默认培养 medium、Python 3.13.5 和 optlang.gurobi_interface 记录在 OR 构建清单。完整构建没有启用 V-ATPase 候选；CoQ9 使用 metadata 模式。继承的模型告警不在本授权内修正；完整 XML 对比不支持额外化学或边界改动。

该版本可表述为“用户指定的暂定 OR 赋值已通过管线接入和限定静态验证”。不能表述为“真实合成酶已确认”“同工酶关系已验证”“敲除生长已复现”或“模型已正式发布”。新单基因指令需按其实际输出另行核验。
