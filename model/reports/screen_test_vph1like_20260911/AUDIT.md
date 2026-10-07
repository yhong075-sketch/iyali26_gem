# Vph1-like 注释版 screen 独立来源审计

核验时间：2026-09-11T00:58:06Z。审计者：`screen_vph1_audit`，独立于执行与后处理代理。仅直接读取本地输入、封存源码、原始返回及旧标签表，以标准库重算；审计新增求解 0、网络访问 0、模型与源码修改 0。本文件为唯一审计写入。

覆盖：**total claims 8 | audited 8 | supported 8 | unresolved 0 | contradicted 0 | unchecked 0**。分母仅为下列本轮 screen 主张，不包括此前蛋白身份／定位研究或旧 33 项审计；“supported”表示在所列计算和来源范围内得到支持，不表示原生功能或模型已获实验验证。

| ID | 核查主张与原始来源 | 判定及独立结果 |
|---|---|---|
| S1 | 输入身份及本次注释范围：`run_manifest.json`、两份模型 XML | **supported**。本次执行输入为 `model_metadata_trna_vph1like.xml`，SHA256 `f77f07f96c8f19d0507112813110dc2afe613ddb31c912c0700927f036b23736`。直接解析并替换目标 geneProduct 后，新旧 XML 树序列化完全相同；目标身份 annotation 保留，实际变化仅名称和证据 notes。 |
| S2 | 沿用上一轮实现与条件：新旧 `run_screen.py`、`run_manifest.json`、封存引擎及 SD-Leu/PO1f 文件 | **supported**。runner diff 仅模型路径、固定 SHA 和实验对照状态文字。封存引擎 9 个源码文件现值与前后 manifest SHA 全同；培养文件、有效培养、菌株 overlay、排除集合、solver 配置和软件版本均与上一轮一致。运行时质粒伪基因不参与 KO。单进程、Threads=1、presolve=false，三项已记录容差均为 1e-7。 |
| S3 | 运行完成与 WT：`execution.json`、`run_manifest.json`、`run_screen.py`、引擎 `run_single_gene_deletions` | **supported**。退出码 0、未超时，外层 7.558847000007518 秒，核心 4.023056040983647 秒。代码路径与记录显示请求 1 WT + 1074 KO 主求解，无重试路径；并非独立后端调用遥测。WT=1.8718823069403，与上一轮相同，已通过引擎 optimal 和 0.1–2.0 h⁻¹ 范围门。未另存 WT 原始 status 字段，status 依据成功越过该门判断。 |
| S4 | 原始 KO 完整且有效：`raw_deletions.tsv` 与运行 gene_ids | **supported**。1074 行、1074 唯一 ID，恰好覆盖指定集合；全部 raw status=optimal、growth 有限且非负。逐行核实导出的 raw、legacy growth 和 ratio 对应原始返回，0 未决数值行。 |
| S5 | 严格四阈值分类：原始 KO、WT、`screen_predictions.tsv` | **supported**。直接用未舍入 `raw_growth/WT < cutoff` 重算 4296 个分类，逐项匹配 audited 与 legacy。1%/5%/10%/15% 必需预测数分别为 **76/93/101/114**。当前结果未触发旧引擎异常转零分支；异常转零政策仍仅属 legacy，不能当生物学死亡证据。 |
| S6 | 相对上一轮预测是否改变：两轮 `screen_predictions.tsv` 与原始 KO | **supported**。4296 个阈值分类均不变。逐基因 raw growth 最大绝对差 `2.0922361065878192e-12`，ratio 最大绝对差 `1.1177204994883283e-12`，全部满足既有比较容差 abs_tol=1e-8、rel_tol=1e-6；数值并非逐位相同。 |
| S7 | 用户关注的目标结果：目标模型 geneProduct 与新旧 screen 行 | **supported**。**YALI1F38820g — 原生正式名未核实 — 偏 Vph1-like 的 V0 a 亚基候选，可能参与质子转运／复合体装配（整理注释，依据序列及既有 AlphaFold 预测；仍需实验确认）**。本次 KO growth=1.8718823069403139，KO/WT=1.0000000000000075，四阈值均预测非必需。约 1 的末位偏差是数值差，不能解释为生长促进。候选注释没有改变 GPR，也不确认原生区室或实验必需性。 |
| S8 | 沿用已封存正例标签的统计：上一轮 `essentiality_per_gene.tsv`、`essentiality_summary.json`，本次原始 KO | **supported**。旧标签表现 SHA 与其 summary 一致。1612 正例中 322 在本次筛查、1290 在筛查外；1074 模型基因中 752 未标注。独立重算四阈值 TP/FN 为 **55/267、65/257、73/249、79/243**，10% 交集内召回率 **73/322=22.67080745341615%**；对应未标注必需预测数 21/28/28/35。未重新打开外部工作簿、解释新标签或推断映射。无负例，不能报告 FP/TN、accuracy、precision、specificity 或 MCC；本次是开发参考回归，不是独立实验验证。 |

主要来源定位：本目录的 `run_manifest.json`、`execution.json`、`raw_deletions.tsv`、`screen_predictions.tsv`、`run_screen.py`；上一轮 `artifacts/screen_test_metadata_trna_20260910/` 的对应文件及 `essentiality_per_gene.tsv`；封存 `artifacts/r608_engineering_20260907/code/scripts/gem_annotate/validate_essential_genes.py:518–576` 和 `essentiality_simulation_context.py`。旧标签表 SHA256 为 `fcad4295baf2d1837a1ee4bc22b2a6b24ea90b75029308b37123f73245eb76a8`。源码、培养及菌株完整 SHA 已逐一对照运行清单。

解释边界：本审计验证本次静态 screen 的执行与回归统计；未重审已有封闭 ATP 生成问题，未完成蛋白定位／功能实验，不推广到其他培养条件，也不以 screen 成功接受新的 GPR 或模型科学变更。
