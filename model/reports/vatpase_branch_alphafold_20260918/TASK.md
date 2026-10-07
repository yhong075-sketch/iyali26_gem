# 原 GPR 分组的多链 AlphaFold 预测

2026-09-18。问题：用户所给前 12 种蛋白（A）与后 2 种蛋白（B），分别作为多链输入时，预测装配及组件功能支持何种异同？不预设两组是完整泵、彼此可替代或必需基因。

## 授权与执行范围

本轮用户明确要求 AlphaFold 预测；项目 AGENTS.md 已授权研究所需 HPCC AlphaFold 作业，无需逐次重复确认。该授权覆盖下面两次有界预测和结果获取/分析。沿用 govern-agentic-research 与 gene-identity-function；独立设计复核见 DESIGN_REVIEW.md，输入来源见 sequence_manifest.json。没有模型/GPR/培养/实验标签变更或 GEM 求解。

1. 固定原分组，逐基因精确核对版本化序列、长度、SHA 与菌株；每种蛋白一条链，不融合、不加 linker，不自行补入原生多拷贝或其他亚基。
2. 使用当前 HPCC 实查的 alphafold/3.0.2，对 A（12 链、4106 aa）和 B（2 链、1000 aa）分别提交一个作业。
3. 两组均 seed=0、5 diffusion samples、10 recycles；运行完整 MSA/模板搜索，模板截止 2026-09-18。没有预填空 MSA 或强制界面约束。无显式膜、脂质、ATP、离子或质子梯度；不能据此模拟 ATP 驱动质子转运。
4. 获取每个样本 mmCIF、完整/汇总 confidence、MSA/模板输入与日志。验收链数和序列，检查局部 pLDDT、PAE、ipTM/pTM、链对界面、clash 及样本一致性，再对可对应组件和实验参考比较。全复合体单个 RMSD 不作为功能等价判据。

## 资源上限与停止条件

- 账户 iwheeldonlab；gpu 分区/qos；各 1 张 A100、gpu_highmem、16 CPU。
- A：256 GiB RAM、24 小时；B：128 GiB RAM、8 小时；总上限 32 GPU 小时、512 CPU 小时（资源申请上限，非实际消耗）。
- 每组仅一个提交和一次执行；--no-requeue，不自动重投、不增加 seeds/样本、不变更序列或扩大资源。先做实际 AF3 输入解析和批处理语法/调度可接受性检查，不运行登录节点上的 MSA 或模型推理。
- 输入 SHA 改变、数据库/模型文件不读、运行错误、OOM、超时、非有限指标或身份不符即保留原日志并停止受影响组。低置信结果也保留，不以调参追求期望结论。
- 排队和运行中不报告为预测完成；没有返回结构时，结构功能比较保持待完成。

## 证据和结论边界

结构输出全部标为“AlphaFold 预测”；相关功能表述为“基于 AlphaFold 预测的功能/装配候选”。两组是原 GPR 一拷贝假设分组，非原生完整化学计量。A 含一条明确 CLIB122 来源的旧 ID 链，W29 对应伪基因冲突保持未决。基因名称/功能与逐条证据详见 SEQUENCE_REVIEW.md；原生定位、酶活、装配及不可替代性不能由本预测单独确证。

优先版本明确的原始序列记录、酵母实验结构和官方工具文档；既有 BLAST/单体结构只作历史证据，不称本次重算。实验参照可能进入模板/训练资料，不构成独立功能验证。未知数据库逐字节快照身份保留未知，实际路径/大小/mtime记录而不伪称已完整校验数百 GB 数据库。

交付：输入及来源、任务脚本/参数、实际环境指纹、提交/状态/日志；成功后追加置信度和组件比较。独立审查以来源覆盖计数，不以代理同意人数计数。未完成的输出不提前验收。

官方部署/硬件/格式依据：
- https://hpcc.ucr.edu/manuals/hpc_cluster/selected_software/alphafold/
- https://hpcc.ucr.edu/manuals/hpc_cluster/intro/
- https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md
- https://github.com/google-deepmind/alphafold3/blob/main/docs/performance.md

在线 main 文档不等于实际代码版本；本次以 HPCC 容器/源码指纹和实际参数解析为准。
