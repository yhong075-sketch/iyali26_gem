本说明冻结现有四个场景的定义，复用既有审阅，不重做全量审计。**本子任务求解器调用0、采样0、文献检索0**，不改模型、旧产物或PROJECT_STATE。记录时间：2026-09-08T20:06:27.668396+00:00。可直接冻结的生产条件为：增长下界0.14650760568106092、脂质伪池需求下界23.74006264072655、正反总通量上限2330.701999560738；无需重新求最小值或增长最优值。

| 精确场景名称 | biomass_C边界 | CFSA_DM_lipids边界 | 正反总通量上限 |
|---|---|---|---:|
| `growing (pFBA)` | [1.3185684511295488, 1000.0] | [0, 1000] | 1223.3176651109031 |
| `growing (pFBA, slow)` | [0.14650760568106092, 0.301250344647618] | [0, 1000] | 310.806297395808 |
| `producing` | [0.14650760568106092, 1000.0] | [23.74006264072655, 1000] | 无 |
| `Producing + cap` | [0.14650760568106092, 1000.0] | [23.74006264072655, 1000] | 2330.701999560738 |

前三个名称来自上游sample的condition字符串；第四个`Producing + cap`是本地comparison.json标签，并非上游第四条件。四者保存的objective均为**max biomass_C，系数1**，但OptGP采的是约束可行域，收集样本时并不继续最大化该生物目标；预热时另做各反应的最小/最大化。不能因为生产场景的存储objective仍为biomass_C，就将其误记为没有施加产物要求。

**全部新增采样约束。** 四场景共享一次性需求`m1727[C_cy] →`，系数−1；原始增长下限为WT最大值的10%，最终边界以表中为准。正常生长改增长下界；慢生长改增长上界（不是固定增长）；生产改目标下界。前两个场景各增加`0 ≤ sum(f+r) ≤ cap`；原producing无额外总量约束；新增生产条件仅增加`sum(f+r) ≤ 2330.701999560738`，该约束未显式设置下界，正反变量本身非负已保证总量非负。没有额外biomass最优等式。原培养基、PO1f运行覆盖、其他全部反应边界与化学计量由对应保存JSON完整定义；本次不重新套用运行覆盖。[上游条件与边界](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/upstream/sampling_tools/sampling.py:117>)、[新增生产约束](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_pfba_cap_20260908/complete_cap.py:139>)。

**上限如何得到。** 令L=∑(f+r)、P_k为第k个表型的原约束集合。两个生长场景先求各自g*_k=max_{P_k}biomass，再以objective_fraction=1暂时要求biomass达到该最优值，求最小L，最后在原P_k上加1.25倍总量上限。这个临时增长最优条件不作为采样等式保留。由保存上限/1.25反推参照总量，正常生长为978.6541320887225，慢生长为248.64503791664637；这是已有记录的代数反推，不是本次求解。

新增生产场景直接在保存的增长和产物不等式下求min_{P_producing}L，既有最优记录为1864.5615996485903，故cap=1.25×该值=2330.701999560738。没有重新最大化增长。因此，当前三个cap共享1.25倍因子，却不共享绝对上限，连参照优化的层次也不完全相同。若未来统一采用“各表型最小总量×同一因子”，必须一致定义各P_k；若改成“所有场景同一个绝对cap”，则是另一项容量假设和对照设计，不能混称。此说明只冻结现有生产cap，不替换成生长组cap。

| 场景 | 完整有效模型 | 样本 | 保存形式 |
|---|---|---|---|
| `growing (pFBA)` | [模型JSON](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/scenario_0.json.gz>) | [原始样本](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/scenario_0.npz>) | net flux only |
| `growing (pFBA, slow)` | [模型JSON](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/scenario_1.json.gz>) | [原始样本](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/scenario_1.npz>) | net flux only |
| `producing` | [模型JSON](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/scenario_2.json.gz>) | [原始样本](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/scenario_2.npz>) | net flux only |
| `Producing + cap` | [模型JSON](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_test_20260908/scenario_2.json.gz>) | [原始样本](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_pfba_cap_20260908/samples.npz>) | net flux and original split variables; explicit fwd_idx/rev_idx saved |

JSON不会保存自定义总通量约束，需与上述明确cap一起使用。新增生产的[LP](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_pfba_cap_20260908/capped_problem.lp>)由既有审阅确认仅多一个cap，但LP浏览导出舍入部分系数；精确重建应使用原完整JSON加已冻结cap。原始样本每组500、thinning100、processes1，历史seed依次20260908、20260909、20260910、20260910。新链的数量、seed及诊断由主任务另行管理，本子任务不执行它们。canonical暂定参考SHA为bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee；实际采样输入是已应用SD-Leu/PO1f及一次性需求的有效JSON，不能称未经改动的canonical XML。

**混合脂质伪池的存储组成。** 下表逐项抄录保存模型的xLIPID（Lipids pool）系数；负号为消耗，正号为生成。名称和分子式均为现有模型注释，未经本次化学或实验来源重验证。

| 代谢物ID | 模型名称 | 存储分子式 | 系数 |
|---|---|---|---:|
| `m1000[C_cy]` | ergosterol_C28H44O | `C28H44O` | -21.176 |
| `m1579[C_lp]` | ergosteryl palmitoleate_C44H72O2 | `C44H72O2` | -0.1 |
| `m1581[C_lp]` | ergosteryl oleate_C46H76O2 | `C46H76O2` | -0.305 |
| `m1631[C_em]` | phosphatidyl-L-serine_C8H11NO10PR2 | `C8H11NO10PR2` | -1.967 |
| `m1640[C_lp]` | triglyceride | `C6H5O6*3` | -12.51 |
| `m1641[C_lp]` | fatty acid_ | `CHO2*` | -3.152 |
| `m1648[C_em]` | 1-phosphatidyl-1D-myo-inositol | `C11H17O13P*2` | -3.6 |
| `m1651[C_mm]` | cardiolipin | `C13H18O17P2*4` | -0.842 |
| `m1700[C_lp]` | phosphatidylethanolamine_C7H12NO8PR2 | `C7H12NO8PR2` | -10.101 |
| `m1701[C_lp]` | phosphatidylcholine_C10H18NO8PR2 | `C10H18NO8PR2` | -14.24 |
| `m1705[C_lp]` | phosphatidate | 未存储 | -3.399 |
| `m1727[C_cy]` | lipids_ | 未存储 | +1000 |

xLIPID每单位反应通量产生1000单位`m1727[C_cy]`伪池，不等于1000 mmol TAG或1000 g脂质。该池在保存模型中的稳态账目为**1000·v_xLIPID − v_xBIOMASS − v_CFSA_DM_lipids = 0**；此式说明归一化尺度，不提供实验质量换算。各系数不能自动当作质量百分比。TAG只是组成中的一个通式成分，输出目标仍是混合脂质池。

canonical SBML定义了`mmol_per_gDW_per_hr`，且部分边界参数（如R1070摄取边界）显式引用；但伪池没有独立物质单位或分子式，保存JSON也不提供这类单位定义。多个成分使用R或*通式，phosphatidate和lipids_缺失分子式。现有本地材料不足以确定该伪池分子量、具体酰基组成或质量归一化依据；本说明**不假定MW，不换算TAG、g/L或g/gDW**。这些缺口限制解释，不阻塞按原数值尺度冻结同一生产采样问题。

完整边界、全部额外约束、源文件SHA、组成字段、单位定义及`freeze_producing_cap`记录见 [conditions.json](</Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/cfsa_multichain_20260908/conditions.json>)。既有源码六项、数值七项和cap九项审阅的分母保持原样；本次只完成四场景定义与一项伪池组成的静态整理。
