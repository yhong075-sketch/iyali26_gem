# AlphaFold 网页结果与 GPR 判断

2026-09-21 状态勘误：初稿末段误称“当前14 AND假设”，现已按实际整理数据和两份模型XML修正。保存的V-ATPase假设版是**三个指定基因共同AND，再接原始两分支OR**；另一份R539标记模型仍保留原始OR。全14 AND仅为未实施的逻辑对照，不能说成已实施或已获准的模型状态。本日仅追加静态文件与布尔逻辑检查，以下AlphaFold获取和分析仍是2026-09-18的工作，未重新运行。

2026-09-18。本次实际查看两组已完成的 AlphaFold Server 预测，取得并核验 A 原始结果；B 原始下载受浏览器阻止。**现有证据不支持原来的“前 12 AND”与“后 2 AND”构成两套可互相替代的完整 ATP 耦联质子泵，也不验证全 14 AND。** 更符合参考装配机制的候选是共同组件 AND 对应的 a 亚基；两个 a 能否 OR，必须逐区室判断。

## 已取得的结果

| 组别 | 内容与获取等级 | ipTM | pTM |
|---|---|---:|---:|
| [A，前12链](https://alphafoldserver.com/fold/7b423a15eaee4c66) | 4106 aa，各1拷贝；实际请求、5份summary、5份坐标身份已独立核验 | 最佳页面/模型0为0.55；5份0.53–0.55 | 模型0为0.58；5份0.57–0.58 |
| [B，后2链](https://alphafoldserver.com/fold/71064986a3c52042) | 1000 aa，各1拷贝；仅根代理网页摘要，原始ZIP未取得 | 页面最佳0.51 | 页面最佳0.68 |

两组页面 seed=1。A 下载请求确认 seed1、12条原序列、各count1、useStructureTemplate=true；准确服务器模型补丁版本、recycles及数据库快照未知，不沿用HPCC参数。A输入和5份坐标均逐链匹配冻结序列；A为11条W29与1条CLIB122混合来源的假设组，不能称为完整原生W29复合体。

A压缩包通过CRC检查；5份完整confidence另外用[可重复分析脚本](analyze_results.py)核验链长、矩阵形状、数值范围和有限性，统计见[confidence_analysis.json](confidence_analysis.json)。完整来源/SHA见[获取记录](results/retrieval_A.json)和[独立审核](RESULT_AUDIT.md)。后者独立覆盖实际请求1/2、summary5/10、坐标身份5/10；A full_data统计由主分析完成，未冒称已经独立复算。

按[AlphaFold官方输出说明](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md)，ipTM反映链间相对位置的置信度；两组当前全局分数都未给出高置信整体排列。**低分不能证明真实蛋白不相互作用；高分也不能证明完整泵功能或AND/OR关系。** 两组大小不同，不能按pTM/ipTM高低比较质子泵能力。

A的局部结构支持与整体低分并存：模型0中，YALI1A09766g（原生正式名未核实；V1 ATP水解A亚基候选）与YALI1E32332g（正式名未核实；V1 B结构亚基候选）链对ipTM=0.81，且5个样本均为0.81。YALI1D00581g（正式名未核实；V1 D中央转轴候选）与CLIB122的YALI0E16192g（正式名未核实；V1 F中央转轴候选）模型0链对ipTM=0.75。上述功能均为同源/预测结构候选，非原生功能实验确认；界面置信度支持局部装配调查，不能单独确认共同必需性。CA-pLDDT由独立审核从坐标统计，full_data的pLDDT按全部原子统计，两者口径不同，未混用。

B普通下载操作未产生文件，对页面已观测Download链接的导航返回net::ERR_BLOCKED_BY_CLIENT，未绕过。没有B原始坐标，因此本轮未完成A/B组件定量结构叠合，不填入旧结构的RMSD冒充新结果。两组各一拷贝输入均缺完整参考构架的成员/拷贝数，这限制了整体装配解释。

## TypeSafe 技能在本任务中的适用范围

用户提供的[官方SKILL.md](https://raw.githubusercontent.com/typesafe-ai/skills/main/skills/typesafe-ai/SKILL.md)已成功读取。其要求的文档索引和Markdown页面未能访问，按说明切换普通网页后成功读取[State](https://docs.typesafe.ai/concepts/state)、[Choice](https://docs.typesafe.ai/primitives/choice)、[Confidence](https://docs.typesafe.ai/confidence)及[引用核验示例](https://docs.typesafe.ai/cookbooks/citation_check)。本轮适用的是方法指导：分开固定事实与待判断问题；分别判断组成、替代能力、区室和菌株；明确支持/反对/未解决；不把输出类型或模型置信度等同于科学真值。

本轮没有构建TypeSafe应用或调用Jev/API，以下判断由本研究依据固定输入、原始来源和静态布尔逻辑给出；没有TypeSafe概率值，也未照搬示例的自动接受阈值。无需为了这次科学判读新增服务、SDK或上传整个项目。

## 对 GPR 的具体判断

| 问题 | 判定 | 理由/边界 |
|---|---|---|
| 后两种蛋白能独立替代前组、执行完整ATP驱动泵反应吗？ | 不支持原外层OR | B为V0 a/c″候选，不含V1 ATP水解头；原OR可在催化头缺失时仍判反应可用 |
| 之前保存的三基因共同AND是否修复了完整组件依赖？ | 没有 | 实际规则仍包含原OR；静态移除V1催化A候选时规则仍为真，不能把这项修改称为完整复合体GPR修复 |
| YALI1F13017g应与核心组件共同AND吗？ | 有机制支持的候选 | 正式名未核实；V0 c″膜转子环候选。参考体系c/c′/c″有非冗余作用；Yarrowia依赖仍需验证 |
| YALI1D00581g应共同AND吗？ | 有机制支持的候选 | V1 D中央转轴功能候选，与催化头协作；并非另一套独立催化分支 |
| YALI0E16192g可以直接确认为W29共同AND吗？ | 未解决 | 旧位点的完整F候选蛋白来自CLIB122，W29对应伪基因/移码注释冲突未解决 |
| YALI1E12482g与YALI1F38820g要共同AND吗？ | 缺少支持 | 两者均为V0 a候选；后者偏Vph1-like但仍需实验确认。未证明同一泵同时需要两个a |
| 两个a可以OR吗？ | 有条件的候选，尚不能验收 | 同角色/折叠相似不等于同一区室可替代；需定位及互补证据 |

两种a位点原生正式名均未核实，候选功能为质子通道及V1–V0连接；其家族证据不能直接确认为VPH1/STV1原生功能。

候选逻辑框架为：

```text
原前12中除YALI1E12482g外的11个共同组件
AND YALI1F13017g
AND (YALI1E12482g OR YALI1F38820g)
```

括号中的OR仅在对应区室可替代被验证后成立。R794是Golgi、R795是液泡：若两种a分工不同，应分别使用该区室的a候选，不能把同一个OR机械复制到两条反应。完整候选式、参考原始实验、证据等级及验收条件见[独立GPR审查](GPR_REVIEW.md)。其中引用的酿酒酵母实验为机制参照，未升级为W29实证。

上述旧F条目只保留为假设中的F功能角色占位；W29的活性编码位点未确认，整式不能直接写作已经核实的W29 GPR。若需要正式修订，两项关键证据缺口是实际菌株F身份及各区室a亚基替代能力。

2026-09-21[实际文件核验](gpr_file_check_20260921.json)分别保存两份XML的完整SHA、原GPR及六个单基因缺失的布尔结果：`model_metadata_trna_vatpase_and_hypothesis.xml`是三基因共同AND假设，`model_metadata_trna_r539_alphafold_labeled.xml`是原始OR；不根据文件时间或名称推定科研基线。这里GPR为真只表示规则允许反应，未证明非零通量可行；本日没有新增LP或screen。

本轮没有修改上述模型、整理数据或管线，也没有新增screen、优化、结构预测或集群操作。**GPR是否正确与该反应是否被生长需求强制使用是两个问题；修正依赖关系不保证基因成为essential。**
