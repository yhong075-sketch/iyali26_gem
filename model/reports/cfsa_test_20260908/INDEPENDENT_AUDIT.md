本次独立审阅完成了 6 项静态核查：5 项 supported，1 项 partially_supported。求解器调用 **0**，模型导入 **0**，未运行 CFSA 或采样。结论限于所提供源码快照和 manifest 指定的暂定参考 XML，**不构成正式模型验证或采样收敛验证**。核验时间：2026-09-08T07:06:22.138753+00:00。

执行计划为：读取适用规则；核对 sample/helper/example 调用关系；执行 NumPy 属性探测；校验参考模型指纹并静态检查 objective 和目标参与反应；分别记录支持证据、反证与限制。预算为 10 分钟，遇指纹变化、越界动作或实质性来源冲突停止受影响核查。只写本报告和 audit.json。

| 项目 | 判定 | 直接发现与来源 |
|---|---|---|
| C1 三个条件 / pFBA | supported | [sample](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/upstream/sampling_tools/sampling.py:129>)：正常生长施加 biomass 下界 b1×optimality；慢生长施加 biomass 上界 slow_growth；生产施加 target_export 下界 t2×optimality。仅前两种调用 [pFBA helper](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/upstream/sampling_tools/util.py:93>)。 |
| C2 NumPy 2 兼容性 | supported | 独立 Python 3.13.7 / NumPy 2.3.3 实际访问 np.NaN 得到 AttributeError，提示它在 NumPy 2.0 被移除。[sample:137](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/upstream/sampling_tools/sampling.py:137>)、[sample:159](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/upstream/sampling_tools/sampling.py:159>)、[util:55](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/upstream/sampling_tools/util.py:55>) 均有访问。 |
| C3 显式 seed | supported | [OptGPSampler 构造](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/upstream/sampling_tools/sampling.py:152>)只传 model、processes、thinning；上游文件中没有 seed 设置。 |
| C4 Geweke 公式和限制 | supported | [util.geweke](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/upstream/sampling_tools/util.py:15>)使用段均值差除以两段原始方差和的平方根；没有样本数归一化或自相关修正。小值不能证明收敛。 |
| C5 objective | supported | [XML active objective](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_integration/model.xml:162516>)是 maximize biomass_C、系数 1。[示例参数](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/upstream/example_sampling_yarrowia.py:26>)使用 xBIOMASS；该旧反应在参考 XML 中仍存在。 |
| C6 脂质目标及输出 | partially_supported | [示例](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/upstream/example_sampling_yarrowia.py:68>)自行添加单代谢物 xlipid_export；参考 [脂质池](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_integration/model.xml:54154>)是 m1727[C_cy]（XML 编码形式见来源），没有该池的单代谢物输出。ID/名称对应有据，跨版本脂质组成和生物学等价未核实。 |

**条件解释。** slow_growth 先在 target_export ≥ t2×optimality 下最大化 biomass 得到，然后在慢生长条件中只作为 biomass 的上界，并非固定增长。pFBA helper 在各条件边界已经施加之后重新计算参照值，再约束 ∑(forward+reverse) ≤ flux_fraction×pFBA objective_value；不能写成“三条件共享同一个 WT 总通量上限”。本次没有求出任何边界值或可行性结果。

**兼容性与复现限制。** 无有效缓存的 sample 路径先调用 verify_target_production（源码中四次 slim_optimize），再求 slow_growth，随后到达第 137 行；因此该 NaN 错误位于真正采样循环之前，却不位于所有优化之前。有效缓存可在第 83–84 行直接返回。sample_pairs 也通过第 421 行调用 sample。独立审阅没有执行这些路径；兼容性证据是代码可达性加 NumPy 属性探测，不能称完整 CFSA 失败已复现。独立环境的 NumPy 2.3.3 与主任务环境分开记录。没有显式 seed 也不意味着两次输出必然不同：实际 COBRApy 默认行为和缓存须另行核查。

**诊断的数学含义。** 令 A 为前 floor(0.1N) 行，B 为后 floor(0.5N) 行，默认公式为

D = [nanmean(A) − nanmean(B)] / √[nanvar(A) + nanvar(B)]。

nanvar 默认 ddof=0。该代码既没有用有效样本数缩放段方差，也没有估计时间相关性。即使两个区段的均值相同、D=0，区段的分布仍可不同，也可能有未访问区域；这是公式直接允许的反例，不是本次采样观察。因此不能仅把常见标准正态 z 阈值套在这个值上，声称诊断已经校准。恒定链产生 0/0=NaN；两段分别恒定而均值不同会产生无穷值，随后转成 NaN。这些都是未定义诊断，不是收敛成功或生物学死亡。太短的链、全 NaN 区段也有限制。源码把无效样本行设为 NaN 后再分链，且留有其他诊断的 TODO（sampling.py:157–173）。

**目标的可迁移范围。** 全量静态检查了参考 XML 的 2,313 个反应、1,877 个代谢物。目标只参与 xBIOMASS（20 个不同代谢物，消耗目标，系数 1）和 xLIPID（12 个不同代谢物，生成目标，系数 1000）；针对该目标的单代谢物反应数为 0。此处“没有输出”严格指该总脂质池，不表示所有脂质分子都没有交换途径。原示例 GEM 未纳入本次来源，因此不能由同一 m1727 词干和 lipid 名称推出组成、单位或生物学等价；直接使用旧 ID 的查找和现有 xlipid_export 假设均不成立。增加需求反应会改变优化问题，其科学解释不由这份审阅接受。

反证与限定：参考中 xBIOMASS 并非不存在；只是未被当前 active objective 选中。示例第 83 行实际拼写为 model.objetive，本审阅没有将其当作已经有效改变 objective 的证据。源目录 README/setup.py 声明了项目网址和版本 0.0.1，但目录不含 Git 元数据，本次未独立认证其精确上游 commit；完整本地源文件 SHA 和模型 SHA 保存在 audit.json。参考 XML SHA 与 manifest 一致。模型培养条件、运行时菌株配置、采样参数与结果均未在本审阅执行，不能以当前默认值补作历史记录。

覆盖率：**total claims 6 | audited 6 | supported 5 | unresolved 1 | contradicted 0 | unchecked 0**。C6 的 partially_supported 计入 unresolved。六项覆盖率不代表对整个软件、主任务执行副本或所有模型科学性质的全面审计。
