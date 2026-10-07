# 三基因共同 AND 假设版：结果独立审计

核验时间：2026-09-11T01:07:14Z。审计者：`/root/screen_vph1_audit`。本代理直接读取新旧 XML、整理数据、封存筛查引擎、已安装 COBRApy 敲除实现及原始结果，以标准库独立解析布尔式和重算统计；新增优化求解 0、联网 0、模型／源码修改 0，仅写本报告。实现审计的 10 项与本报告分开，不重复旧 33 项生物学来源审计。

**结论：三共同 AND 已实际生效，三个单敲都会关闭 R794/R795 的 GPR；但在本次 SD-Leu/PO1f 静态生长目标下，三个单敲仍可达到约 100% WT 生长。它们在四个阈值下仍预测非必需。** 这是用户选择的待验证 GPR 假设的计算结果，不能升级为原生非必需或功能／定位已实验确认。

覆盖：**total claims 8 | audited 8 | supported 8 | unresolved 0 | contradicted 0 | unchecked 0**。分母仅为以下模型文件、计算与统计主张；生物学未决项保留，不计入本次完成的计算验证。

| ID | 直接核验来源与主张 | 判定及结果 |
|---|---|---|
| H1 | 新旧模型 XML、整理 JSON、构建清单：实际输入和变更范围 | **supported**。执行模型为 `model_metadata_trna_vatpase_and_hypothesis.xml`，SHA256 `b00a9ea20c21712b428f159f9050727053cf373ae89c6cfe161e7a77e0b2eb64`，与构建及 screen 清单一致。逐项 XML 比较确认只有 R794/R795 的 GPR 与假设 notes 改变，其余 XML 全同。目标反应计量、边界、物种属性及其他对象保留；R795 仍为 [0,0]。 |
| H2 | XML 的 geneProductAssociation、已安装 COBRApy `Gene.knock_out` 与 `_gene_deletion`：布尔效果与敲除路径 | **supported**。独立逐一计算每条 GPR 的全部 2^14=16384 状态，新规则恰等于“三个 required AND 原 GPR”，真值数由 4099 变为 1025。三个单敲分别令两条 GPR 为 false；COBRApy 的实际敲除路径将非 functional 反应 bounds 设为 (0,0)。原内部 OR 未删除，没有改成全部 14 基因 AND。 |
| H3 | 前后 screen runner、封存源码、输入及运行清单：条件和方法保持 | **supported**。runner 仅变更目录深度、模型路径和固定 SHA。9 个封存引擎源码 SHA 与前轮一致；本轮记录的全部源码 SHA 当前可核对。培养文件、有效培养、PO1f overlay、排除集合、solver/software 和资源预算均同前轮。单进程、Threads=1、presolve=false，三项已记录容差均为 1e-7。模型、培养及菌株文件实际 SHA 均与执行清单一致。 |
| H4 | `screen/execution.json`、`run_manifest.json`、`raw_deletions.tsv`：完成性、WT 与有效值 | **supported**。退出码 0、未超时，外层 6.796206500031985 秒、核心 3.865603332989849 秒。请求 1 WT + 1074 KO 主求解，无重试；该数量来自代码路径与请求记录，不是独立后端遥测。WT=1.8718823069403，已通过引擎 optimal 与 0.1–2.0 h⁻¹ 门。1074 原始行具有1074唯一ID、完整覆盖指定集合，均 optimal、有限且非负，0 未决数值行。WT status 未另存原始字段，其有效性由成功越过引擎门判断。 |
| H5 | 原始 KO、WT 和导出表：四阈值结果 | **supported**。逐行核对 raw、legacy growth、ratio 和数值有效标记；使用未舍入 `raw_growth/WT < cutoff` 独立计算4296项分类，与 legacy/audited 全部匹配。1%/5%/10%/15% 的预测必需总数为 **76/93/101/114**。没有异常值被当作生物学死亡。 |
| H6 | 与 `screen_test_vph1like_20260911` 原始结果比较：变化量 | **supported**。WT 差为0，4296项分类均不变。原始 KO growth 最大绝对差 `1.4518247715145094e-12`，ratio 最大绝对差 `7.755948661092305e-13`；每项均满足既有比较容差 abs_tol=1e-8、rel_tol=1e-6。仅138个原始生长数值完全相同，不能将近似一致写成逐位一致。 |
| H7 | 三目标原始行、新模型布尔逻辑及菌株配置：目标结果的正确解释 | **supported**。下表三基因的单敲均 optimal、KO/WT≈1，四阈值均预测非必需。菌株 overlay 不覆盖 R794/R795；本次结果支持“假设 GPR 关停生效，但当前模型与培养条件允许维持 WT 最优生长”。未保存目标反应的每次原始通量，不把逻辑核验写成新增通量测量；未诊断具体旁路或需求缺口。 |
| H8 | 旧封存正例标签、本次 raw 与 `screen_summary.json`：对照统计和来源边界 | **supported**。1612个正例中322个在本次筛查、1290个在筛查外，正例覆盖19.975186104218362%；752个筛查基因未标注。四阈值TP/FN为 **55/267、65/257、73/249、79/243**，未标注且预测必需为21/28/28/35，10%交集内召回率 **73/322=22.67080745341615%**。新 summary 数值与独立重算一致，旧ID、source ID、行号、功能原文和原实验标签均保留。三个目标在沿用表中均为 unlabelled，不能据此判定实验非必需。 |

| 系统 ID／原生名称状态 | 简要功能与证据状态 | 本次原始 KO growth | KO/WT | 四阈值模型判定 |
|---|---|---:|---:|---|
| YALI1D00581g／正式名未核实 | V1 D 中央转轴亚基候选；序列和既有 AlphaFold 预测支持，整理注释 | 1.8718823069403128 | 1.0000000000000069 | 均非必需 |
| YALI0E16192g／正式名未核实 | CLIB122 V1 F 中央转轴亚基候选；序列和既有 AlphaFold 预测支持，W29 活性位点／伪基因对应冲突未决 | 1.8718823069403032 | 1.0000000000000018 | 均非必需 |
| YALI1F38820g／正式名未核实 | 偏 Vph1-like 的 V0 a 质子转运／复合体装配亚基候选；序列和既有 AlphaFold 预测支持，仍需实验确认，区室／替代性未决 | 1.871882306940245 | 0.9999999999999707 | 均非必需 |

主要来源为本目录 `screen/{run_screen.py,run_manifest.json,execution.json,raw_deletions.tsv,screen_predictions.tsv,screen_summary.json,essentiality_per_gene.tsv}`，前版 `artifacts/screen_test_vph1like_20260911/` 对应文件，两份模型 XML，以及 `data/reference_build/curation/vatpase_gpr_hypothesis.json`。封存引擎位于 `artifacts/r608_engineering_20260907/code/scripts/gem_annotate/validate_essential_genes.py:518–576`。本次直接沿用标签表 SHA256 为 `7232721349df674c100f24195b7039d714217cce8a2d66cf52d33110c0a9f54a`；又核对上一轮原始封存标签 SHA `fcad4295baf2d1837a1ee4bc22b2a6b24ea90b75029308b37123f73245eb76a8`，未重开外部工作簿或重做映射。

此对照没有实验负例，FP、TN、accuracy、precision、specificity、MCC不可计算；已用于开发的正例参考不是独立实验验证。此前封闭 ATP 生成问题仍限制模型科学解释，本次未确定它是否造成此处无生长效应。不得据结果直接接受该 GPR 假设、宣称功能／区室已确定，或认定所有培养条件下均非必需。
