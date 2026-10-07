# Escher essential / inactive 注释独立来源审阅

核查时间：2026-09-09T22:02:20.519802+00:00。审阅对象由另一代理生成；本审阅只读取源交付物、独立派生集合并逐行比对，没有调用求解器，没有修改模型/GPR/实验calls/历史结果。

结论：全部直接核查通过。322 E、96 H、752 unknown、68 closed 与各自明确来源一致；当前参考复现证据只覆盖既有固定6 KO，不能外推全量。

## 适用指令与实际覆盖

- 继承本任务已读取的全局与工作区指令。逐一检查 research 源文件的逻辑祖先和 resolve 后真实祖先中的 AGENTS.md / AGENTS.override.md，未发现附加指令文件；没有把其他工作树规则自动继承到 research 路径。
- 11 个 annotation.sources 文件全部读取并核 SHA；正例 CSV 1612行、历史 TSV 2364行、mini-replay TSV 14行全部读取，非小样本抽查。静态 XML 的1074基因和2313反应全部检查。CoQ显示层48行逐字段检查。
- 科学声明限于“这些现成来源与本图显示数据一致”；没有重新验证原实验calls、跨版本身份、原生蛋白功能、历史归一化前求解轨迹、完整历史环境或独立验证资格。

## 核查记录

| 检查 | 结论 | 实际覆盖／证据 |
|---|---|---|
| 源文件 SHA: reference.xml | 通过 | /Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/r608_c7_model_test_20260907/inputs/reference.xml |
| 源文件 SHA: consensus_essential_genes.csv | 通过 | /Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/state/essentiality/repository/consensus_essential_genes.csv |
| 源文件 SHA: essentiality_per_gene.tsv | 通过 | /Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/artifacts/results/essentiality/trna_biomass_group_b/pipeline_runs/ce6d752991bb-20260826T072416Z/essentiality_screen/essentiality_per_gene.tsv |
| 源文件 SHA: run_manifest.json | 通过 | /Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/artifacts/results/essentiality/trna_biomass_group_b/pipeline_runs/ce6d752991bb-20260826T072416Z/essentiality_screen/run_manifest.json |
| 源文件 SHA: essentiality_summary.json | 通过 | /Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/artifacts/results/essentiality/trna_biomass_group_b/pipeline_runs/ce6d752991bb-20260826T072416Z/essentiality_screen/essentiality_summary.json |
| 源文件 SHA: panel_results.tsv | 通过 | /Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/static_mini_replay_20260905T190318Z/panel_results.tsv |
| 源文件 SHA: result.json | 通过 | /Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/static_mini_replay_20260905T190318Z/result.json |
| 源文件 SHA: run_manifest.json | 通过 | /Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/static_mini_replay_20260905T190318Z/run_manifest.json |
| 源文件 SHA: BENCHMARK_CONTRACT.md | 通过 | /Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/docs/BENCHMARK_CONTRACT.md |
| 源文件 SHA: baseline_manifest.json | 通过 | /Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/docs/baseline_manifest.json |
| 源文件 SHA: run_manifest.json | 通过 | /Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/artifacts/results/essentiality/campaign189_0f3a6c2b_20260806_review/current_sha_diagnosis/run_manifest.json |
| 基因 ID 全覆盖 | 通过 | 源 XML 的全部 1074 个 ID 精确匹配，无跨版本补充，无运行时伪基因。 |
| 322 E 精确 ID | 通过 | 1612 条正例CSV全读，322精确交集逐ID核对。 |
| 752 unknown 保留 | 通过 | 未标注项全部 null；没有将其补成负例。 |
| 实验来源逐行保留 | 通过 | 322 个来源原始行字段逐项完全一致。 |
| 历史全量行身份 | 通过 | 2364 行全读，1074 模型内唯一行与源XML基因集合完全一致；其余1290为模型外正例。 |
| 96 H 与严格阈值 | 通过 | 按未舍入 float(ratio) < 0.1 独立重算，1074 行全部与历史 predicted_essential_primary 一致。 |
| 历史数值与原始导出行 | 通过 | 1074条raw_row原样保留；模型比值/状态逐行精确一致。 |
| 历史导出有限性 | 通过 | 仅核实导出结果有限且optimal；不证明归一化前异常不存在。 |
| 历史 source 与当前 source 分离 | 通过 | manifest+summary均指向07dbdc7b；图源bc2aac8f。未把历史H重命名为本次current全量复现。 |
| 历史阈值字段 | 通过 | 阈值源值0.1；严格<与不舍入符合冻结合同的现成代码审阅。 |
| 历史异常未冒称排除 | 通过 | 保留legacy归一化前None/NaN/负值未知限制。 |
| 正例交集内H重合 | 通过 | 67仅为322评价交集的历史命中；不建立独立验证。 |
| 固定6 KO范围 | 通过 | 14行=两个XML各WT+6 KO；仅canonical_model六个基因进入current-reference标记。 |
| 固定6 KO源身份 | 通过 | result.json加载输入哈希与当前reference一致。 |
| 固定6 KO原始数值 | 通过 | panel全部6行对照result原始求解记录，值和状态一致；本次仅核实历史交付物。 |
| 固定6 KO阈值和剩余未知 | 通过 | 2个<0.1，4个不低于0.1；其余1068个current-reference字段为null。 |
| 2313反应边界完整核查 | 通过 | 全部2313反应逐一对照原XML上下界。 |
| 68 closed 精确集合 | 通过 | 仅上下界同为0者设bound_closed；与源XML精确集合一致。 |
| 无FVA结论迁移 | 通过 | 其余允许边界不被视为active；全FVA/blocked状态仍unknown。 |
| 旧inactive诊断版本排除 | 通过 | 旧inactive输入0f3a6c2b与参考不同；未向本图迁移标签。 |
| 48条CoQ身份显示保真 | 通过 | 48条name/function/status/scope/source_references/accession/crosswalk/URIs全部与TSV逐字段相同。 |
| atlas派生数据与annotation一致 | 通过 | 仅添加identity显示字段，没有改写E/H/current/closed标签；两类输入哈希匹配。 |

## 阈值、历史与未知值的边界

对历史表的重算仅为 `essential = (float(ko_growth_ratio) < 0.1)`。不舍入，无epsilon，恰好0.1不判essential。所核1074行的历史导出预测与该式全部一致。此处重算标签不是重新求解。
历史运行源记录保留07dbdc7b模型、PO1f/SD-Leu、glucose uptake 10、uracil runtime uptake 1000、strain overlay、代码commit及dirty=true身份；它不是bc2aac8f静态参考的全量筛查。当前图静态边界没有加入该运行时覆盖。历史导出虽然均optimal且有限，原函数可能已把异常归零，故不宣称排除原始None/NaN/负增长。
68个closed的结论来自上下界均0。它们是明确关闭的反应，不等同于查明全部稳态blocked集合。其余2245个反应也没有被宣称active；blocked状态保持unknown。essential关联的基因不自动证明反应essential，原始AND/OR规则不能被颜色摘要改写。
固定6 KO结果来自2026-09-05已交付的mini-replay，2个低于10%阈值；本次只核对原始记录和输入身份。其余1068个基因的当前参考状态为null。6个对照不能认证07dbdc7b与bc2aac8f生物学/优化问题等价，也不能更新历史67/322召回率。

## CoQ显示层复查

首轮阅读发现逐行 `sources` 证据指针、accession及crosswalk/annotation URIs没有进入显示对象，已即时反馈。作者补入 `source_references`、`accession`、`local_crosswalk`、`annotation_uris` 并重建。复查最新 atlas-data.json 的48条身份对象，8类显示字段逐项等于原TSV，且保留源路径和完整源SHA。未向上升级证据状态，也未替换候选措辞。

## 本次审阅文件身份

| 文件 | SHA256 |
|---|---|
| annotation_status.json | `9bee46f8f63a926c7f6da84431cec548b5ea4324b3fb8fcb1c2cacb374515fbd` |
| atlas-data.json | `92be7e851aa1285478eddf95a8fb0cee6f5220eedc786e8213c0965a18cd410b` |
| data/coq9_gene_evidence.tsv | `4318e3212f548582e053f81f02cd46dfa6e3892a2587251c252993d21cd9e8c0` |

其余11个源文件的绝对路径与完整SHA保存在annotation_status.json.sources，已逐一核对。报告未追加第二份基线登记。
