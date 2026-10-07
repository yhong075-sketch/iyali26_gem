# AlphaFold Server 网页提交范围：独立只读核查

核查日期：2026-09-18。已检查本路径祖先及本目录，适用仓库根 `AGENTS.md`，未发现本目录 override。本轮只读官方网页/文档并写本文件；没有操作浏览器、登录、提交、取消集群作业、修改既有输入或模型。

**A 的 4106 aa 和 B 的 1000 aa 均低于官方教程列明的每任务 5000-token 上限。** 本次各链均为标准未修饰蛋白，1 aa 对应 1 token，因此不需要为长度而截短。两组应保持分别提交、12/2 条链、每种 1 copy、无 linker 的既定假设。

| 核查事项 | 所读官方依据与实际支持范围 | 本任务判断 |
|---|---|---|
| 每任务长度 | [EMBL-EBI 官方 AlphaFold Server 教程](https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/alphafold-server-your-gateway-to-alphafold-3/)，What AlphaFold Server can do：上限 5000 tokens，未修饰氨基酸每个计 1 token；单条蛋白最短 4 aa | A=4106，B=1000；所有链均大于 4 aa。这个文档判断还须与实际提交页的当时校验一致 |
| 分子/链数量上限 | [Google DeepMind Server JSON 文档](https://github.com/google-deepmind/alphafold/blob/main/server/README.md)支持多 entities 和蛋白 `count`。本轮公开可读取的官方文档未查到独立的最大链数数字；[Server FAQ](https://alphafoldserver.com/faq)通过网页检索打开返回 0 行文本 | **精确链数上限未核实**。不能因此说无限，也不能臆造 12 链超限。实际页面是否完整接收 12/2 链由根任务核验 |
| 多链输入格式 | 上述 DeepMind 文档：顶层是 job 数组；`dialect=alphafoldserver`；`version=1`；每条蛋白为 `proteinChain`，内含 `sequence` 和 `count` | 不原样上传本地 `alphafold3` dialect 文件。网页输入转换后应逐条比对原始序列、顺序和 `count=1` |
| 随机种子 | 同文档 Job name, seeds and sequences：`modelSeeds` 是 uint32 的**字符串列表**；空列表由服务器选一个随机种子 | 如果页面接受，可显式用 `["0"]` 保留原种子意图。若选择自动种子，应保存实际请求/结果并记录不同，不能回填为 0 |
| 网页与本地格式/控制不同 | [本地 AF3 官方输入文档](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md)，AlphaFold Server JSON Compatibility，明确两种 dialect 不同；网页格式不允许自行指定链 `id`，转换器按列表分配 | `version=1` 是**输入格式版本**，不能据此推断网页模型补丁版本。结果须按序列重新确认链映射 |

## 网页执行和 HPCC 执行不能默认等同

两处都属于 AlphaFold 3，并不证明网页后端等于 HPCC 的 3.0.2，也不证明相同权重、MSA/模板数据库、模板截止时间或采样实现。网页版需记录其实际暴露的设置和下载包元数据；后端未给出的代码/权重版本、recycles、数据库发布版本保持未知。HPCC 计划的 seed 0、5 samples、10 recycles 不能自动抄入网页运行记录。

Server 官方格式提供 `useStructureTemplate` 和 `maxTemplateDate`；本轮文档写出的模板库上限日期是 2025-02-03，但这只是文档当时描述，不证明本次后台数据库快照。应保留提交界面接受的值、job request 和结果元数据。若服务器只暴露默认值，就标为服务器默认而非宣称与本地一致。

发现一处文档时效冲突：EBI 教程仍称不能自定义 MSA/templates，但当前 DeepMind Server README 已列出 `unpairedMsa` 和 `templates` 字段。本次无需使用这些字段，不依赖过时限制判断可提交性；也不能从文档字段存在推断用户当前网页已展示该选项。

本次查阅未获得服务器明确的运行版本号或承诺与本地同配置的官方材料。因此结论应是“对同一组已固定序列另做网页版 AlphaFold 预测”，不能写成“完全复现本地 AF3 3.0.2”。无论结果如何，仍需遵守 [DESIGN_REVIEW.md](../DESIGN_REVIEW.md) 的真实拷贝数、混合菌株、完整功能与结构置信度边界。

## 提交前的最小验收

- 实际预览包含 A 的全部 12 链／4106 aa，B 的全部 2 链／1000 aa；每条 count=1，无额外修饰、配体、截短或融合。
- 两组分别保留 request、实际 seed/模板设置、job ID、提交时间和结果下载包；序列全文及 SHA 对照既有 manifest，不靠任务名认定输入一致。
- 若实际页面拒绝 A 或提示额外分子数量/长度限制，保留错误及原输入；不能自动拆分或截短后声称完成 A 的原任务。B 可独立继续，但须明确 A 尚未提交。

本核查已经直接打开上述官方页面；未下载或宣称保存其原始网页快照。官网 FAQ 的动态正文在此只读网页工具中不可得，这限制“精确最大链数”和“当前网页后台参数”的声明，不能把这些未知值写成已核实事实。

## 网页输入独立核查与实际 UI 记录补充

2026-09-18 22:17 UTC，本审核直接重新读取 [two_branch_jobs_seed1.json](two_branch_jobs_seed1.json)，逐项对照原始 [sequence_manifest.json](../sequence_manifest.json) 与 HPCC [branch_A.json](../branch_A.json)、[branch_B.json](../branch_B.json)，未依赖根任务的自检标志。**14/14 条序列、完整序列 SHA、12/2 顺序和 count=1 全部一致。** 网页文件两项 job 均为 `alphafoldserver` 格式版本 1，均显式指定 `modelSeeds=["1"]`；HPCC 两项则为 seed 0，属于已记录的执行差异。总长度仍为 A 4106、B 1000，未增加修饰、配体或额外蛋白字段。

这份网页输入文件的 SHA-256 为 `95f9c7d84437fee5fb5aade0986f8eba2fbd529c87658089b8a7f74f6dfef97d`，与 [WEB_SUBMISSION.json](WEB_SUBMISSION.json) 的最终预期输入指纹一致。本次读取的提交记录 SHA-256 为 `3738cea238fdca79ccbb48b882fb7259a599e771c88bdc4d589b5d273947ecce`；其中两组预览链长、名称、copy 数和 seed 均与文件吻合。这只独立证明了本地保存输入与观测记录之间的一致性，尚不是下载服务器 job request 后的逐序列回验。

**以下是根代理观测记录，未独立重现 UI：** 根代理直接阅读 `alphafoldserver.com/faq` 记录每 job 5000 tokens、每 seed 5 samples、默认模板日期 2021-09-30；实际预览 seed 输入的 min/max 分别为 1／2147478647。seed 0 触发 `rangeUnderflow` 且提交按钮禁用，改为 1 后启用。因此本文件早前“若接受可用 seed 0”的条件未满足，实际采用 seed 1。FAQ/界面所示模板日期与上文官方 README 提及的最大可设置日期不是同一概念，不把任一值当作已核实的后台数据库发布版本。

根代理的提交记录称 A/B 在当地 15:14/15:15 各点击一次确认提交，出现 `Job launch pending`，历史条目转为进度条，状态为 `in_progress`。这提供了实际接收本次 12/2 链的根代理观测，不能升级成“已独立验证平台最大链数”。当时没有可见 job ID、结果 URL 或结构；相应字段保持 null/false。预测完成及实际结果数仍待结果包验证。

本轮未操作网页或 HPCC。根代理另报告取消 HPCC 作业被自动审批拒绝、已向用户提问且作业保留；这里不自行重试、取消或把网页提交解释为集群作业已经终止。原生功能、真实拷贝数和 W29 F 位点身份的限制保持不变。
