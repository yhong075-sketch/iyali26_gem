# iYLI647 液泡 ATPase：独立来源审计

2026-09-14；审阅者 `/root/r795_open_audit`。直接读取原始发表附录 XML、corr_3 JSON、两种旧 screen 的 manifest/原始 KO/预测表、执行代码及原论文网页。新增 LP、模型加载修补和模型修改均为 0。已审阅本目录 REPORT.md，核心结论与下列直接证据一致。

**裁决：受检两个版本未解决本次液泡/Golgi ATPase 的可运行性与必需性问题。正上界和开放水运输不足以解除其 ATP/ADP 完整守恒行造成的结构阻断；旧基因 KO 的下降也不是液泡依赖证据。** 不推广至未检查的后续衍生版本或原生生物学。

## 直接核查

1. **方程及边界。** 两版 ATPS3v/g 均为 `ADP[q]+Pi[q]+3 H+[c] → ATP[q]+H2O[q]+2 H+[q]`（q=v或g），[0,1000]；ATPS3m 的对应计量也相同。这是 ATP 合成方向的模型表示，不能仅按名称等同于本项目消耗胞质 ATP 的质子泵。
2. **完整物种行。** 对原 XML 全部1343个 reaction 元素及 corr_3 全部1353条反应分别扫描，ATP[v]/ADP[v] 仅关联 ATPS3v（+1/−1），ATP[g]/ADP[g] 仅关联 ATPS3g（+1/−1）。原 XML 四物种 boundaryCondition 均为 false。因此，在标准稳态质量平衡下 `v_ATPS3v=v_ATPS3g=0`，不是某个最优解恰好不用它们。
3. **不同名称的反应。** 全反应扫描所有含 ATP 和 v/g H⁺ 的步骤，结果也只有 ATPS3v/g。按允许的正反方向检索“ATP消耗、ADP/Pi产生、胞质H⁺消耗、v/g H⁺产生”的直接泵模式，两版均无匹配。此检查不声称排除一切多步间接能量耦联机制。
4. **水运输反证。** 两版 H2Otv 均为 `H2O[c] ↔ H2O[v]`、[-1000,1000]。它们确实没有本项目 R1363 原关闭的同一断点；但水运输不改变上述四条 ATP/ADP 行，不能使 ATPS3v/g 运行。这也不说明全部其他液泡反应都阻断。
5. **全部基因关联。** 对 corr_3 所有反应的 GPR 逐项提取：每个 p 后缀记录仅控制 ATPS3v 和 ATPS3m；每个 g 后缀记录仅控制 ATPS3g，没有第三条关联。原 XML notes 的原始 token 扫描给出同样关联；原 v/m GPR 有未闭合括号，本次没有修复或将其当作已成功解析的布尔表达式。

七组原始 token 为下表。模型所有对应 gene.name 均为空；本轮未核实原生正式名称或真实亚基功能，功能证据仅为模型 ATP 合酶相关赋值（model/GPR assignment only）。p/g 作为不同模型记录保留，未合并为同一原生基因。

| p后缀系统记录：ATPS3v/m | g后缀系统记录：ATPS3g |
|---|---|
| YALI0A09900p | YALI0A09900g |
| YALI0A11143p | YALI0A11143g |
| YALI0B03982p | YALI0B03982g |
| YALI0B06831p | YALI0B06831g |
| YALI0B11913p | YALI0B11913g |
| YALI0B21527p | YALI0B21527g |
| YALI0D00583p | YALI0D00583g |

## 旧 screen 复核及归因

这部分重新读取旧输出而未新增求解；执行输入均为 corr_3，不是原始发表 XML。对两个情景的全部14个目标分别读取 raw_deletions.tsv、screen_predictions.tsv、run_manifest.json，重算28个未舍入 KO/WT，并核对112项阈值分类：

| 情景 | WT h⁻¹ | p组7个 KO/WT | g组7个 KO/WT |
|---|---:|---:|---:|
| native | 1.1537351153250899 | 0.5828749056639301–0.5828749056639357 | 约1 |
| mapped28_po1f | 1.2171871273706656 | 0.5870639489734654–0.5870639489734777 | 约1 |

28项原始状态均 optimal，生长有限非负，原始值与预测表一致；所有目标在既有1%、5%、10%、15%阈值下均非必需。manifest 未记录 GPR 变更，mapped28_po1f 另有 OMPDC 关闭及培养映射；不将它称为与本项目 SD-Leu 完全相同。

执行代码使用基因删除改变相关反应可用性。p组同时关闭 ATPS3m；ATPS3v 已被完整物种行强制为0，且关联扫描没有其他反应。因此这组下降不能作为液泡酸化必需性修复的证据；能影响本模型可行域的相关变化是线粒体 ATPS3m 的关闭。此为模型数学归因，不确认这些 token 对应蛋白的原生区室功能。

## 来源身份与限制

- 原附录 SHA256：`e1a3521606d94840dde7005ff4395e178fd24f387eb688b51f30b0794a704d77`；[发表论文及Additional file 2入口](https://link.springer.com/article/10.1186/s12918-018-0542-5)。直接枚举实际1343个reaction/1117个species，正文报告1347/1152；差异保留，未推测原因。
- 原附录另有两个 biomass reaction 缺少 id。本审阅以元素位置作临时审计定位，仍纳入完整化学扫描；没有更改源文件或声称原 XML 已成功载入优化器。
- corr_3 SHA256：`329be540c099409c2c7b76ee581a23f86eefdbbffe97e9b03b39b5e9c014b5d2`。实际1353反应/1122物种；该 SHA 与两个旧运行manifest相同。其固定提交来源见旧 inputs/source_manifest.json，不把它与2018原附录视为同一文件。
- 本轮未重新运行全部 screen、核验全部旧求解见证或核实原生蛋白身份。审阅范围是完整相关化学行/所有直接泵候选、14个相关模型token的全部关联、28个已保存KO结果及限定结论。

## 审计覆盖

| ID | 本轮声明 | 判定 |
|---|---|---|
| I1 | 两个文件身份、实际计数及原附录格式限制被区分并保留 | supported |
| I2 | 两版 ATPS3v/g 为指定 ATP 合成方程且边界[0,1000] | supported |
| I3 | 四条完整 ATP/ADP 行强制两反应稳态通量为0 | supported |
| I4 | 两版全反应扫描未找到不同名称的直接 ATP 水解质子泵 | supported |
| I5 | H2Otv 已开放但不解除 ATP/ADP 断点 | supported |
| I6 | p/g记录的全部关联分别为v+m或仅g；原GPR未修补 | supported |
| I7 | 28个旧KO比例/状态及112项阈值分类与原始数据一致 | supported |
| I8 | 旧p组生长下降不构成液泡酸化必需性已修复的证据 | supported |

**total claims 8 | audited 8 | supported 8 | unresolved 0 | contradicted 0 | unchecked 0**。

分母仅限这8项文件/计算声明。原附录格式与正文计数差异、原生亚基身份/方向/定位、生物学酸化需求及未检查衍生版本仍不在已验证结论范围内。
