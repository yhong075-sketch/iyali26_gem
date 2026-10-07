# ATP 能量异常候选修复与构建回归修正

2026-09-24。**已实际完成候选修复：最终 E5 重新加载后，封闭 ATP/GTP/UTP/CTP 的最大耗散均为 0，全部 optimal；真实培养仍可行，生长 1.42917789228 h⁻¹。** 默认构建和历史基线文件保持，未正式推广候选。不是仅提出方向建议，也没有为恢复旧生长数字改培养基、biomass 或 NGAM。

## 实际改动与构建覆盖根因

核心实现为 `scripts/gem_annotate/energy_candidates.py` 与 `data/energy_candidate_repairs.json`；独立入口 `scripts/build_energy_candidates.py` 和主构建 `--energy-candidate` 共用该补丁机制，默认 E0。`main.py`/`cli.py` 接入开关，`reaction_selection.py` 接入优先级保护。验证脚本 `scripts/validate_energy_candidates.py`、两组候选测试、README 与项目状态更新分别承担计算、回归及入口说明。既有正常 R694、R594、R603、ATP 维护、V-ATPase 和 CoQ 反应没有修改。

实际调用 gap-fill、charge-aware balance、metadata selection、候选、再次 metadata、SBML 导出/重载，保存35份完整反应快照。使用的是当前物种目录下的窄范围运行重放，**没有声称重建不可访问的完整历史环境**。第一次错误覆盖发生在 `apply_metadata_reaction_selection`：它按 `data/metadata_reaction_selection.json` 的显式化学/边界选择覆盖早期整理；不是名称更新本身误改方程。

| 目标 | metadata 前实际状态 | metadata 后实际覆盖 |
|---|---|---|
| R_PGAM1_PhosHydro | PEP水解，单向 | 计量反号为PEP合成，边界恢复双向 |
| R_NTP3pp | 无额外H⁺的GTP水解，单向 | 反号为GTP合成、加入反应物H⁺、恢复双向 |
| R_NTP7 | 无额外H⁺的UTP水解，单向 | 加入产物H⁺、恢复双向 |
| r0242 | DHA+Pi→DHAP+H₂O+2H⁺ | 删除产物2H⁺ |

修复明确规定**显式候选化学/方向优先于旧metadata**。完整反应/物种签名先全部通过，才开始修改；不匹配报错，重复应用不累加。定义与证据锁存入SBML notes，后续metadata先验证锁再处理，保护化学/边界，其他整理照常执行。最终导出检查完整反应、物种、目标、实际求解器变量/约束/矩阵；自定义诊断约束或稳态行修改不能静默消失。见[运行阶段差异](build_stage_reaction_diff.tsv)与[构建报告](BUILD_REPORT.md)。

## 候选实际修改的八条反应

全部边界按最终XML存储方向，原始已开放方向的容量不增大；GPR和区室保持原值。

| 反应 | 最终边界 | 计量修订 | 依据及限定 |
|---|---|---|---|
| R_PGAM1_PhosHydro | [-1000,0] | 无 | 保留PEP水解，禁止无耦联PEP合成；[EC3.1.3.60](https://iubmb.qmul.ac.uk/enzyme/EC3/1/3/60.html) |
| R_NTP3pp | [-1000,0] | 去掉反应物H⁺ | 保留GTP水解；中性GTP/GDP/Pi约定，[KEGG R00335](https://www.kegg.jp/entry/R00335) |
| R_NTP7 | [0,1000] | 去掉产物H⁺ | 保留UTP水解；同一物种约定，[KEGG R00159](https://www.kegg.jp/entry/R00159) |
| r0242 | [-1000,1000] | 产物补2H⁺ | DHAP为二价阴离子，DHA/Pi为所存中性形式；未额外改变其方向，[DHAP身份](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:57642) |
| R_NDP1 | [-1000,0] | 去掉反应物H⁺ | 保留ADP→AMP+Pi水解，禁止无耦联ADP合成；[KEGG R00122](https://www.kegg.jp/entry/R00122) |
| R72 | [0,1000] | 无 | 保留IP7→IP6+Pi水解；其冲突的kinase注释所需ATP/ADP不在实际式中，[EC3.6.1.52](https://iubmb.qmul.ac.uk/enzyme/EC3/6/1/52.html) |
| R_CAT2p | [0,1000] | 无 | 保留ethanol[cy]+H₂O₂[pe]→acetaldehyde[cy]+2H₂O[pe]的peroxidatic氧化；该式没有O₂，[EC1.11.1.6](https://iubmb.qmul.ac.uk/enzyme/EC1/11/1/6.html) |
| R_OAADCm | [-1000,0] | 去掉产物H⁺ | 保留草酰乙酸脱羧；OAA/pyruvate/CO₂的中性形式明确，[EC1.1.1.38的第二反应](https://iubmb.qmul.ac.uk/enzyme/EC1/1/1/38.html) |

五处H⁺修订依据是具体物种身份、所存分子式/电荷及官方反应定义，不是任意抵消总路径残差；没有修改共享分子式。**八条修订反应在最终所存约定下逐条元素和电荷均配平。** 中性NTP/Pi的部分结构注释仍指向带电别名，这一旧注释不一致没有被伪装为全部解决。原生精确底物活性、现有GPR和混合区室抽象未被这些化学类别来源确认。完整原/新定义及证据见 [accepted_candidate_fixes.tsv](accepted_candidate_fixes.tsv)、[独立化学审查](chemistry/REPORT.md)。

## 累计版本与实际结果

各版本均通过独立XML构建和重载后求解；E1/E2用于区分方向与计量，后续迭代从E3累计，没有回到E0反复搜寻原来的三条路径。

| 版本 | 累计内容 | 封闭ATP最大值 mmol/gDW/h | 真实培养生长 h⁻¹ |
|---|---|---:|---:|
| E0 | 原始加载定义 | 1000 | 1.87188230694 |
| E1 | 原三方向修订 | 1000 | 1.87188230694 |
| E2 | 原三处计量修订 | 1000 | 1.87188230694 |
| E3 | E1+E2 | 1000 | 1.87188230694 |
| E4 | E3+R_NDP1方向/计量+R72方向 | 489.655172414 | 1.87188230694 |
| E5 | E4+R_CAT2p方向+R_OAADCm方向/计量 | **0** | **1.42917789228** |

所有上述结果均optimal。1000表示达到测试容量上限，不表示数学无界。E0、E3的GTP/UTP/CTP测试也均1000，E4三者约489.655172414，**E5三者均0**。E1/E2未额外测试这三种载体。耗散测试使用已核实的中性NTP/NDP/Pi与水，不额外增加H⁺；ATP直接使用原有配平的xMAINTENANCE。未沿用此前不配平的memote ATP模板，也不把四种核苷酸通过外推到NAD(P)H、FADH₂、质子梯度、所有内部循环或全模型热力学。

真实加载后关闭183条单侧物料边界及6条biomass列，保留1877条物种稳态，只在封闭副本释放maintenance下界为0；检查并拒绝其他强制通量/自定义约束。每个验证批先实际求解精确零通量，均可行。最大耗散与临时D=1/L1见证分开；干预后恢复D范围和最大耗散目标，在完整网络重测。全部变更、完整变量/约束、状态及全向量保留在各批manifest和逐解JSON。

## 修订后找到的新路径与断路检验

**E3见证：** R_NDP1无耦联合成ADP，R121逆向腺苷酸激酶转移磷酸生成ATP；同时通过R72逆向、R73、R284、R422及相邻转运平衡肌醇磷酸与H⁺。12条非零列含维护，L1=6.33333333333。其中R284/R422所涉缺式仍不可完整化学验证，不能用总账本配平替代逐列确认。

从同一个E3分别阻断R_NDP1正向，完整ATP最大值仍489.655；仅阻断R72逆向仍1000；两者联合仍489.655。这些对照先记录为 `diagnostic_blocks`；只有另经来源审查支持的修订才进入E4。正常R121等未因参与循环而删除。

**E4见证：** R_CAT2p逆向由乙醛与水产生乙醇/过氧化氢，R195再由过氧化氢产生氧；乙醇/乙醛往返及NAD(H)、呼吸链和ATP合酶闭合。另有R_OAADCm无耦联羧化方向和R533/R538参与平衡。19条非零列含维护，L1=18.6153846154。所有H⁺、水、Pi、核苷酸、氧化还原辅因子、跨区室运输均保留，见[完整见证](atp_witness_fluxes.tsv)。

从同一E4单独阻断R_CAT2p逆向，完整ATP最大值为0；单独阻断R_OAADCm正向仍241.379310345；联合为0。没有关闭正常ATP合酶来通过测试。源审核支持CAT2p的过氧化物酶式正向活性与OAADCm的脱羧方向，才将它们写入E5；CAT2p混合胞质/过氧化物酶体的定位仍保留不确定性。

两份见证都实际固定D=1并最小化其余绝对通量总和，不声称严格最少反应数。逐物种加权验证 `Σ(j≠D) S_j v_j = −S_D`，净式均为 **ADP[cy]+Pi[cy]→ATP[cy]+H₂O[cy]**；44行账本覆盖全部抵消物种，[净反应账本](atp_net_reaction_ledger.tsv)。E5满足最大值≤1e−7的预设停止条件，因此未再强行求D=1，也没有继续无目的枚举其他循环。

## 生长变化、营养与正常能量正控制

原PO1f/SD-Leu/质粒选择加载器、交换约束、biomass_C目标和NGAM下界7.8625均恢复。E5生长比E0低约**23.6502%**，没有调高摄取或放宽维护恢复历史数字。

两项额外、同E4起点的实际培养归因对照：仅加CAT2p方向修订，生长1.42995141661；仅加OAADCm方向/化学修订，仍1.87188230694；合并后1.42917789228。故本次下降主要由移除CAT2p逆向机会造成，OAADCm在CAT2p修订背景下还有小幅影响；这是该模型条件下的对照结论，不是测得生物生长缺陷。

| 保存的FBA解 | E0 | E5 |
|---|---:|---:|
| 葡萄糖摄取 | 10 | 10 |
| O₂摄取 | 0.0664881719 | 14.8250769269 |
| CO₂排出 | 0.169842937 | 14.8595587392 |
| R_CAT2p通量 | −15.5957319605 | 0 |
| 正常ATP合酶R171通量 | 0 | 53.4850209070 |
| xMAINTENANCE通量 | 7.8625 | 7.8625 |

通量单位mmol/gDW/h，摄取表以正数列出。真实培养的正常ATP合酶仍能提供能量，支持修复没有简单切断ATP生成能力。这些是各自一个最优解，未证明唯一代谢流或生物学上调；完整实际营养输入与ATP相关通量在[growth_comparison.tsv](growth_comparison.tsv)。E3保存最优解的maintenance通量为66.3127，但其下界仍7.8625；不能把最优解多余耗散当作新增生理需求。

## 验收边界、保护与可重现性

**A，已验证：** 三个初始错误方向在最终重载XML保持限制；八条候选修订逐列配平；候选不被metadata覆盖、签名冲突拒绝、幂等、导出/重载一致；最终四种封闭核苷酸耗散为0且零通量可行；真实培养生长正控制通过。最终23项软件/行为回归通过，其中最终E5的完整ATP与培养生长分别重新求解验证，而非只读旧结果。E3阳性异常作为负控制，避免把未通过误判通过。

**B，候选层支持：** 官方酶学支持所保留的水解、过氧化物酶式氧化与脱羧方向；统一到现有物种表示的质子修订有身份依据。无需据此任命或替换原生基因。

**C，仍未决：** 原生精确酶活/GPR/区室；二肽身份/供给；肌醇缺式；R305、R1889在E4见证中的H/charge残差及CoQ耦联。CoQ既有审批门保持，本轮未修改它们。四能量测试通过不消除这些局限，也不代表全模型化学与热力学完全正确。当前结果支持将E5作为可回退的候选交付，**不自动替换正式基线或宣称实验校准完成**。

固定历史基线 SHA 是 `aad701126d12d113816fda4b872333b614b4469ee1b6d8ab8c419231c89e965f`，2314反应/1877物种；本轮核实与历史一致。候选E0重序列化与原XML不字节相同：20个tRNA私有残基原无显式charge属性，COBRA读入默认0并导出为显式0；这是已披露的序列化差异，不是新的电荷证据。原基线字节保持，加载后的共享物种定义和全部2314条GPR在候选中保持。先前dirty工作保留；没有commit、push、集群或工作区外研究目录访问。

软件和阶段失败完整保留：窄重放缺物种的预检、严格矩阵导出检查的两次测试修正、E4一个未经核实URL的更正、E5预检说明误把peroxidatic写为dismutation。错误E5未进入优化，最终E5采用已纠正说明；已执行的旧版本不回写。Gurobi复制会细微舍入R1372系数，本轮闭合副本该列明确为0；真实培养和候选导出均从新读XML开始。实际矩阵与物种守恒分别审计，容差没有放宽。

本轮累计58次实际优化调用（含行为测试），实际调用耗时与总墙钟分开记录在[budget.json](budget.json)和各批manifest；所有调用单线程、单次60秒、Presolve0及1e−7既定容差。未复用历史LP作为新候选结果。初始和最终验证器源码身份分别封存；新增负值结果判定只改软件错误分类，不改变任何优化问题。

原始来源审核25项，22项支持、3项原生证据未决；计算/构建独立审核的有限覆盖和全向量核验见[AUDIT.md](AUDIT.md)。全部文件身份、最终保护与实际命令见收口验证记录和[commands.txt](commands.txt)。

## 重新运行

在同一工作区，从原始固定XML生成新目录中的全部候选：

```sh
.venv/bin/python -B scripts/build_energy_candidates.py --variants E0 E1 E2 E3 E4 E5 --output-dir artifacts/atp_candidate_repair_20260924/rebuild_new
.venv/bin/python -B -m scripts.validate_energy_candidates --output artifacts/atp_candidate_repair_20260924/recheck_new --budget artifacts/atp_candidate_repair_20260924/recheck_budget.json --model E5=artifacts/atp_candidate_repair_20260924/rebuild_new/E5.xml --other-energy E5
```

输出目录必须新建，禁止覆盖既有结果。构建使用当前修正后的来源说明，原E4执行时的说明另有封存；数学字段差异与说明差异分开。阶段重放/测试实际执行命令见 `build_commands.txt` 和 `commands.txt`。正式培养恢复与封闭情景只发生于独立加载/副本，诊断排出、来源、强制D和L1目标均不写入候选XML。

最终收口：独立计算审核10/10支持，来源审核25项中22支持、3未决，合计有限声明集 **35 total | 35 audited | 32 supported | 3 unresolved | 0 contradicted | 0 unchecked**。独立核验58份完整通量向量、134224个通量值；未以第二个求解器重新证明最优性。全部记录最大数值违反约1.09e−11，小于预设1e−7。

保护清单原有452项中448保持、4项授权代码/README变更；另有PROJECT_STATE独立前态快照，合并范围453项中448保持、5项授权变更。完整变更与文件SHA、总任务墙钟、实际求解调用1.1834046723秒分别见[verification.json](verification.json)。最终E5来源与说明、构建回归、四能量及培养检查均已完成；正式模型推广与原生功能接纳仍不自动发生。
