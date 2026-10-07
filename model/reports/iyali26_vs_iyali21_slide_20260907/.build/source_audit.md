# 单张模型比较幻灯片：独立来源审阅

核验时间：2026-09-07T17:54:09.506675+00:00。审阅者：独立子代理 slide_source_audit。

## 范围与方法

只读打开本地 iYali21 XML、manifest 指向的 canonical XML、PROJECT_STATE 与 baseline manifest；使用 Python 标准库 XML 元素解析、计数、关系检查及 SHA-256。求解器调用 0、模型修改 0、联网 0；不重新执行历史筛查，不外推其他工作线完成情况。检查当前及相邻 integration 工作树的适用 AGENTS；无适用祖先 override 或读写子目录指令。适用科研治理技能采用直接来源审阅与限定主张，不触发科学模型变更授权。

## 审阅结论

| ID | 原子主张与直接来源定位 | 裁定 | 使用边界 |
|---|---|---|---|
| C1 | 对两个 XML 的 core:reaction、core:species、fbc:geneProduct 全部元素计数，分别为 2285→2313、1868→1877、1058→1074；净变化 +28、+9、+16。canonical SHA 与 manifest `artifacts.canonical_model` 一致。 | supported | 是两个静态文件的元素数比较。species 含区室和模型辅助物种，不能等同独立化合物数；geneProduct 是模型编码元素，非新发现基因数。运行时 strain overlay 不计入静态数。净变化不能改写为“只新增了”。 |
| C2 | canonical 有 20 条 `R_TRNA_BIOMASS_*`，旧 XML 为 0。每条由 charged tRNA 产生一个专用残基并释放相应 tRNA；20 个残基各只出现在该生成反应与 `R_biomass_C`。两文件目标函数均指向 biomass_C。直接定位：新 XML 161186（biomass_C）、162095–162515（20 条耦合）、162516–162522（目标函数）；旧 XML 34780（biomass_C）、34921–34927（目标函数）。 | supported，须采用更正后的链路 | 可写“20 条 tRNA–生物量耦合，把氨基酸活化接入生长需求”。实际残基消费端是 `biomass_C`。`xAMINOACID` 在两个文件中仍消耗游离氨基酸，不能写成 xAMINOACID 被改为残基需求。新 XML 144851、旧 XML 28794 可直接核实这个反例。是模型结构进步，未在本次证明预测性能或生物学机制的实验有效性。 |
| C3 | 对 reaction/species/geneProduct 分别统计包含 RDF resource 的元素，新 XML 为 1743/2313、1628/1877、892/1074；旧 XML 三类均 0。新 XML 全部 RDF li resource URL 的主机为 identifiers.org。 | supported | 仅表示这些本地文件所保存的结构化数据库链接覆盖。不能写“iYali21 完全没有注释”或“新数据库链接均已实验验证”；旧文件仍保存名字等字段。 |
| C4 | PROJECT_STATE 29–33 行及 manifest `benchmark`（1043 起）、`label_dependencies`：10% 严格 < 阈值，TP 67/FN 255，交集内召回 67/322=20.81%；正例覆盖 322/1612=19.98%。manifest 历史结果角色为 development_and_regression_reference，independent_validation_established=false。 | supported，历史记录级 | 不是 iYali21 同条件前后提升值；不应画成“性能提高 20.81%”。本次未打开所有历史原始 calls 或重跑筛查，支持的是当前项目登记的历史数值与评价定位。实际历史执行 XML 另有 SHA，不能把本次静态 canonical 计数冒充该历史筛查重跑。 |
| C5 | PROJECT_STATE 20–24、39–41、94、117–119 行：canonical 是暂定参考；CoQ9 为未校准敏感性专题，lipid-unlump 为未激活候选。最新调度授权探索性计算，不等于默认模型激活或科学接纳。manifest gaps 保留相应门。 | supported | 从“当前模型已完成进步”区排除 CoQ9/lipid-unlump 科学完成项。可以在待办区写“探索/候选中”，不能推断它们没有工程或研究进展。 |

总主张数 5 | 已审 5 | 支持（限定表述）5 | 未解决 0 | 反驳 0 | 未查 0。

C2 的 xAMINOACID 消费端误读已主动通知制图者。若仍用该错误链路，相关图示应判 contradicted；当前可支持链路为：氨基酸 + ATP + tRNA → charged tRNA → 专用残基 + tRNA → biomass_C。

## 完整输入身份

- `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/data/iyli21.xml`: `6974b7588f2a6c60ba2cde2f26e20d3aba1334d0d501572bf01cee47eda86631`
- `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_integration/model.xml`: `bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee`
- `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/PROJECT_STATE.md`: `519165debd7d07baa8ceb4d80547ab5445bd86b446bd9aa124b60a19fc27ba62`
- `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/docs/baseline_manifest.json`: `af952e321d139a0f7a7704dbdcbe631f51e79bea8de24d90a54a95268e4e7e70`

当前目录 model.xml 未用于上述比较。静态参考是 manifest 指定的相邻 integration/model.xml；这项选择延续已有暂定参考登记，不产生正式发布或候选激活。
