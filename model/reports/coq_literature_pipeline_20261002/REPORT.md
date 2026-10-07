# CoQ 文献修订已接入主构建管线

用户授权：“并入管线修改”。2026-10-02 本次已从原始 `data/iyali26.xml` 实际完成一次完整构建，约62.44秒；没有用成品候选代替原始输入。产物为 [coq_literature_pipeline.xml](coq_literature_pipeline.xml)。它是解脂耶氏酵母 *Yarrowia lipolytica* iYali26 的局部 CoQ 文献修订候选，未登记为正式接受模型。

## 接入方式

主入口新增默认关闭的 `--coq-literature-revision`。该选项自动包含既有 C5 与 COQ9 功能依赖候选，并在标准 tRNA 生物量阶段及原有整理步骤之后应用 R19/R18 修订。CLI 和共享 API 使用同一实现；只接受离线、无求解、metadata 模式及新输出路径，拒绝与其他实验覆盖混合。

```bash
.venv/bin/python -m scripts.gem_annotate \
  --research-root ../iyali26_gem_research \
  --offline --no-solve --coq9-curation metadata \
  --coq-literature-revision --output-model coq_literature_candidate.xml
```

本轮复用已有整理数据，不改其中的科学内容。原先独立修订函数的身份检查新增两种保存表示的等价处理：单项列表与字符串注释、HTML 编码与文本箭头。真实 EC、说明、化学、GPR、边界和蛋白版本漂移仍拒绝。共享 SBML 保存工具未改。

## 本次验证结果

- 17项相关测试通过，包括新增选项保护及表示兼容检查；详见 [测试日志](tests.log)。
- 构建前核验前轮57个输入/缓存文件均未漂移；构建期间源码未变，构建后全部记录输入与源码仍匹配。
- 与前轮[已审核局部候选](../coq_literature_revision_20261002/coq_quinol_partial_candidate_v2.xml)相比，完整 XML 解析树相同，完整模型定义、注释和求解器数学结构相同。**序列化字节不同**，不声称文件逐字节一致。
- 2314条反应、1879个物种、1073个基因一致；20个原本未填写电荷的物种仍未填写，没有补成0。R19、R18、R39、R695、R385 的精确元素/电荷残差均为0。
- 新产物重读后再次应用修订不再改变模型。完整构建受离线/无求解保护；独立静态复核的优化与网络调用尝试均为0。未执行 FBA、表型筛选、结构预测或集群作业。

本轮数值与文件核验见 [verification.json](verification.json)，可运行 `verify_pipeline.py` 重查。原生 [build sidecar](coq_literature_pipeline.build.json)记录实际输入、源码/整理数据身份、运行条件及软件完成状态；`build_execution.json`保留命令和起止时间。未恢复历史完整 dirty 环境，本次源码身份单独保存，不以 Git HEAD 代替。

## 科学状态

构建记录显式保留 `partial_literature_candidate_redox_gap` 和 `pathway_closure_validated: false`。DMQ9H2 仍缺少下游氧化态连接，R385 化学问题仍未闭合；静态结构仍强制 R18、R19、R695 稳态通量为0。成功生成候选不等于通路完成、原生催化功能确认或生长必需性验证。证据和未决项沿用[前轮文献报告及独立来源审核](../coq_literature_revision_20261002/REPORT.md)，本轮未新增生物学主张。

原始模型、前轮结果和默认构建行为保留；没有提交或推送。独立实现审计见 [implementation_audit](implementation_audit/AUDIT.md)。
