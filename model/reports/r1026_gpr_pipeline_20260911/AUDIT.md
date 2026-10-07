# R1026 管线实现与目标单敲：独立审计

日期：2026-09-11。审计结论：**当前交付通过，未发现阻塞问题。R1026 已按授权采用单基因 GPR；目标 KO 确实同时关闭 R1026/R2202，但本次固定条件下生长仍接近 WT。**

目标：YALI1F28274g — 正式原生基因名未核实 — Nce103-like β 类碳酸酐酶候选，催化 CO₂/HCO₃⁻ 转换（自动注释、既有 AlphaFold 预测辅助及已授权模型赋值；原生目标特异酶活/定位未确认）。

审阅范围：已读 TASK.md、before.json、前轮赋值审计；相对封存 before 检查实现、主调用、整理 JSON、测试、构建/运行记录；独立解析实际新旧 SBML 并重新检查全部原始通量。两次 LP 由主执行线程运行，本审阅者未调用优化器或修改实现/模型。

## 实现与交付范围

- **范围守卫通过。** `apply_r1026_gpr_assignment` 仅接受指定 reaction/status/前后规则；在写入前检查 R1026 的计量、物种公式/电荷/区室、边界、目标 accession 及 R2202 的反号计量/相同目标 GPR。已有不同 GPR 或冲突 notes 时拒绝；原有 notes 保留，重复应用为 `already_correct`。
- **主调用通过。** 位于最后的 metadata selection 与 CoQ9 处理之后、确定性 SBML 写出之前；默认管线记录 `post_selection_gpr_assignment`。此次构建是 offline/no-solve，CoQ9 为 metadata，没有启用 V-ATPase AND、R608 或容量候选。
- **证据状态保留。** 整理 JSON/实际输出 notes 记录原空 GPR、用户选择酶促解释、定位冲突与未确认催化。历史引用的含 `&` URL 保存在整理 JSON，模型 notes 指向已存原文；没有改全局 writer。
- **必要软件检查通过。** 最终 `tests_pipeline_format.log` 为 8 项检查通过，包含真实 annotation normalizer、GPR 传导、原 notes 保留、往返、幂等及六类冲突拒绝。首轮完整构建在身份守卫处失败，未写模型；修复只允许同一个 accession 的字符串或单元素列表表示，不允许额外 accession。失败日志保留，第二轮改用新输出名并成功，64.1909 秒。
- **独立全 XML 检查通过。** 对比 `model_metadata_trna_vph1like.xml`，新模型只有 R1026 的 GPR 与所列 notes 改变。独立解析核查该反应其他元素、原 notes 与新增 notes，并在替换回原 R1026 元素后确认整棵 XML 完全相等；反应内部元素比较仅忽略外尾缩进。2315 反应、1877 物种、1074 基因保持。`build_validation.json` 列出的 10 个保护文件本次均重新核对 SHA 相符；不将这项有限检查扩张为整个历史 dirty 环境已重建。

输出模型为 `model_metadata_trna_r1026_gpr_build2.xml`，SHA256 `b16e13050d378f2bb097a43f649908e85374f6c626ee34e430338e8f170cf222`。首轮预定的 `model_metadata_trna_r1026_gpr.xml` 未被当作成功输出。

## 原始结果及全通量复核

实际输入模型 SHA、封存 SD-Leu/PO1f 的文件 SHA、历史有效培养/菌株配置、实际加载源码 SHA 均经核对。独立从新 SBML 读取全部计量列，与前轮已固定的有效条件快照逐列一致；本轮全部基础边界也与该历史快照一致，目标函数为最大化 `biomass_C`。

实际参数：Gurobi/optlang；Threads=1、TimeLimit=60、Presolve=0，FeasibilityTol/OptimalityTol/IntFeasTol 均为 1e-7。外层预算 180 秒，实际进程 3.3011 秒；代码及记录只包含 WT、目标 KO 两次主求解，无重试。

| 项目 | WT | YALI1F28274g KO |
|---|---:|---:|
| 原始状态 | optimal | optimal |
| 原始生长目标 | 1.8718823069403 | 1.8718823069402954 |
| 全反应净通量数 | 2315 | 2315 |
| 顺序求和最大守恒残差 | 5.684341886080802e-14 | 6.328271240363392e-14 |
| 独立 `math.fsum` 最大守恒残差 | 8.058630802595945e-15 | 4.9826106646643983e-14 |
| 最大边界违例 | 0 | 0 |
| R1026 通量 | 0 | 0 |
| R2202 通量 | −14.12795132372213 | 0 |

本审计重新计算了全部 1877 个物种的 Sv、完整边界与目标；表中顺序求和数值逐项复现主记录，`fsum` 差别属于求和舍入，不改变判定。两份完整通量均有限，原始值与其 repr 一致。KO 相对 WT 唯一的边界变更是 R1026/R2202 同时成为 `[0,0]`。

KO/WT = **0.9999999999999977**。按原有严格 `<` 的 1%、5%、10%、15% 阈值均为非必需；主阈值仍为 10%。该结果属于本次目标重测，不是新全基因筛查。

本次 KO 解中，R1025/R1361 分别为 `+10.452530866202782`，R1050/R1161/R1028 分别为其负值，符合前轮已逐物种验证的核内形成及运输路径。它是本次可行解中的供给见证，不能证明路径唯一、真实上调或真实原生区室定位。

审计过程中发现并修复了脚本保存 `+inf` 时 JSON 拒绝非有限数值的问题；修复保留原 repr/历史归一化结果再停止，没有把 `+inf` 改成 0。独立提取实际 `solve` 函数，用假求解器检查 `+inf`、`−inf`、NaN、None、负目标 5 类输入均先保存原值再停止；该检查没有导入真实模型或新增 LP。本次真实两次求解均为有限、非负结果。

## 生物学未决

用户授权和软件验收不升级证据：历史 A0A1H6PNW2 v32 的膜间隙 IEA 与同序列 Q6C0V4 v117 的胞质 IBA 冲突仍在；既有 AlphaFold 预测不解决定位；相关 pCA 重组表达研究未隔离确认目标特异催化。R1026 历史意图、真实自发容量与核内供给的生物合理性保持未决。该任务没有改动这些科学结构，也没有为得到必需性结果关闭核内记录或改变培养/评价阈值。

## 本次独立审计指纹

| 记录 | SHA256 |
|---|---|
| `run_manifest.json` | `6a336be58b1f42e0b690f032295f27a4603a139cb84b0a41ab363a875aefefa2` |
| `build_validation.json` | `ac6fa067fd6edb978fea64d101b63ef0aefdbfcf08c06bf52b38d5fd2e073f4c` |
| `scripts/gem_annotate/patches.py` | `b5ded7dc9f25e30ea68ab585862ef22567feda0295cde44311f506883a03c538` |
| `scripts/gem_annotate/main.py` | `570f2e13524d365f27cb4ad9a92c32d20df229e1454f29a77197b5497900fdab` |
| `run_ko.py` | `c537c38e1bbc378dbd0233c742b4cd99dee4185a8df83e8a55959daa8370d600` |

审计覆盖：实现/授权范围、全 XML 差异、受保护文件、实际求解设置、两份完整原始向量与分类均已核查；生物实验、额外区室真实性、路径唯一性及其他基因 KO 的新求解均不在本次验证声明内。

## 最终报告措辞复核

已读取 REPORT.md、README 的本次新增段及 PROJECT_STATE 的本次 R1026 段。核内路径的逐步方向、五个实际通量与列恒等式相符；“有净水合作用”“不证明唯一性/原生定位”“仅目标重测”的限制表述准确。其他任务同时新增的状态段不属于本次审阅。

历史 FN 归类已另行核对 `artifacts/screen_test_metadata_trna_20260910/essentiality_per_gene.tsv`：目标唯一记录保留来源 ID `YALI1_F28274g`、Sheet1 第 622 行，`experimental_label=essential`、`classification_at_10pct=FN`。本次 KO 比值仍大于 10%，所以“对既有 essential 标签仍为 FN”成立；此处复核的是已存评价记录，不称重新审定原实验或独立验证。该表本次读取 SHA256 为 `fcad4295baf2d1837a1ee4bc22b2a6b24ea90b75029308b37123f73245eb76a8`。

最终审阅通过，无需增加求解或扩大来源检索。
