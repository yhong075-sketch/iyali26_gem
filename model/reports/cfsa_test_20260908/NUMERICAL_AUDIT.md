独立复算了保存的 **3×500** 个样本，全部满足重建的净通量约束，容差 **1e-7**。本审阅求解器调用 **0**、重采样 **0**，未将模型导入 COBRA，仅解析保存数据；不是正式模型验证。核验时间：2026-09-08T07:13:29.127024+00:00。独立环境：Python 3.13.7、NumPy 2.3.3；与原运行记录的 Python 3.13.5、NumPy 2.4.3 分开记录。

执行计划为：读取保存证据和代码；核对保护指纹与场景参数；独立重算质量守恒、边界、净L1及lag1；写两份数值审阅文件。预算5分钟，发现指纹变化、实质来源冲突或需要新增求解时停止受影响核查。只有本报告及numerical_audit.json可写，原结果、模型和主脚本未修改。

| 条件 | 投影可行样本 | max 质量守恒残差 | max 边界违反 | max 净L1 | 记录上限 | lag1 中位数 |
|---|---:|---:|---:|---:|---:|---:|
| Growing | 500/500 | 5.04e-12 | 4.1e-10 | 1222.097442 | 1223.3176651109031 | 0.925858 |
| Slow growing | 500/500 | 9.86e-10 | 3.62e-09 | 307.131335 | 310.806297395808 | 0.908576 |
| Producing | 500/500 | 3.89e-11 | 1.34e-13 | 77197.307379 | 无 | 0.909434 |

每份数组为500×2314，反应ID唯一且与场景模型完全对应；场景模型各含1877个代谢物。原保存validation标签均为500个v。净L1均值为 1092.348891、288.404643、45545.487993。完整独立值及与原诊断的差值保存在numerical_audit.json；质量守恒求和顺序可产生末位差异，判定使用既定容差，未放宽门槛。

**净样本检查范围。** 保存量为v=f−r。对非负f、r，|v|≤f+r，故净L1是原正反变量总量约束的必要检查。对这里重建的常规分解边界和总量上限，可取f=max(v,0)、r=max(−v,0)，得到总量恰为净L1的分解，支持净样本的投影可行性。**这不能恢复或直接审核原正反变量轨迹。** 场景JSON不含自定义pFBA上限，数值取自result.json；不对未保存的其他自定义约束声称全面审计。[保存代码](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/run_cfsa.py:148>)明确记录了该限制。

**相关性与大通量。** 按峰峰值>1e-7筛选，每组1369个变化反应；逐反应对相邻499对样本计算Pearson相关，lag1>0.9的数目为941、779、766，中位数与保存诊断一致。应写“自相关高，收敛和充分混合尚未建立”；高相关本身不能证明非平稳或未收敛。Producing没有pFBA总通量上限，保存样本净L1确实较大；本次没有逐环分解，不能归因到具体循环。

**耗时与计数。** result.json记录 113.003019 秒、13480 次backend调用；时间戳差为 113.003166 秒。声明预算是**1200秒、14000次**，113秒是耗时。[计数代码](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/run_cfsa.py:104>)支持记录口径，但无逐次调用轨迹，且计数器在上下文加载及solver设置后安装；没有在本次复现历史耗时或重建所有后台操作。sampler_retries记录为1、0、6，不能自动等同于优化器求解失败重试。

**输入、种子与需求。** 独立重算16个保护文件的SHA均匹配result记录；这证明当前字节匹配，不证明历史全程从无瞬时写入。三个seed为20260908、20260909、20260910，计划/结果一致且[构造器](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/run_cfsa.py:154>)显式传入，随机轨迹未重采样验证。三个保存模型都含系数−1的m1727[C_cy]单向需求CFSA_DM_lipids，上界1000；正常/慢生长下界0，生产下界23.7400626407。[内存定义](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/run_cfsa.py:125>)与保存模型一致。应称**总脂质（混合脂质）伪池需求通量**，不能称TAG产率、生理分泌量或已接受的工程靶点；本审阅不激活脂质候选。

本次数值审阅：**total 7 | audited 7 | supported 6 | unresolved 1 | contradicted 0 | unchecked 0**；历史耗时和计数无法独立重建，N2为partially_supported。原六项源码审阅分母保持：**total 6 | audited 6 | supported 5 | unresolved 1 | contradicted 0 | unchecked 0**，没有覆盖旧文件。完整来源、SHA、逐项限制及环境见 [numerical_audit.json](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/numerical_audit.json>)。
