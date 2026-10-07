# R1931方向证据独立审计

核验日期：2026-09-15。范围：4项核心主张；直接解析实际SBML、独立读取官方条目和原论文可访问段落、核对官方记录快照与序列。不运行优化、不修改模型。**核心主张总数4，已审4，范围内支持4，未解决核心主张0；下列生物学和访问限制继续保留。**

| 主张 | 独立核验与判定 |
|---|---|
| C1 实际模型化学和方向 | 直接解析 `model_metadata_trna_r153_merged.xml`：R1931为胞质GSA + NAD + 水 → 谷氨酸 + NADH，5个物种均为C_cy，各系数绝对值1，边界[-1000,1000]，EC注释1.2.1.88；R332使用NADPH和谷氨酰磷酸，EC1.2.1.41，属于不同步骤。与快照一致。支持。 |
| C2 官方酶类定义 | [IUBMB EC1.2.1.88](https://iubmb.qmul.ac.uk/enzyme/EC1/2/1/88.html) Comments明确称“irreversible oxidation”；旧EC1.5.1.12在2013年转入此编号。[EC1.2.1.41](https://iubmb.qmul.ac.uk/enzyme/EC1/2/1/41.html)明确含磷酸、NADP及谷氨酰磷酸。酶类定义支持方向建议，不能把反应等号读成生理可逆。支持。 |
| C3 原始逆向酶学的准确含义 | 独立读取[Moxley 2014原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC3916563/) Results的Steady-state and Single-turnover P5CDH Kinetics段：E. coli酶用NADH+谷氨酸测试未检出逆活性，拟合化学步因而设不可逆；同段50 mM是另一个产物抑制试验，不能用作逆向配方。独立读取[Small & Jones 1990摘要](https://pubmed.ncbi.nlm.nih.gov/2211729/)及[原文托管索引p.18670](https://www.researchgate.net/publication/20943789_Pyrroline_5-carboxylate_dehydrogenase_of_the_mitochondrial_matrix_of_rat_liver_Purification_physical_and_kinetic_characteristics)：对象为大鼠肝，未检出P5C，正文限制为若发生逆向则至少慢15000倍，降低pH后仍未检出。这是检测约束，不是测得非零速率，更不是跨物种模型边界。报告表述符合来源。支持。 |
| C4 候选边界和结论边界 | [REPORT.md](REPORT.md)把[0,1000]列为保留氧化方向的候选；1000是保留建模容量，非实测容量。未宣称已实施、W29逆向严格零、已计算ΔG、必然恢复必需性或胞质定位得到实验确认。此分层建议由C1–C3支持；模型预测后果仍待授权验证。支持。 |

## 身份核验与未决项

**YALI1B12993g：正式原生名称未核实；推定P5C/GSA脱氢酶；证据为模型/GPR赋值及自动同源注释，非该W29蛋白原生酶学验证。** 独立读取本地官方NCBI原文快照：[XP_500688.1](https://www.ncbi.nlm.nih.gov/protein/XP_500688.1)与[AOW01471.1](https://www.ncbi.nlm.nih.gov/protein/AOW01471.1)均标明W29/CLIB89和目标locus；独立提取ORIGIN序列，与候选FASTA、历史UniProt条目和AF API序列比较，5者完全一致，572 aa，序列SHA-256 `1f78e50dace6d21e517e7d1849c508e0d285ceff6d7568ff56416a7e7b93b565`。取回的4份官方响应文件SHA均与检索记录匹配。

当前UniProt快照的DELETED原因确为不属于reference proteome；历史注释为unreviewed、同源推断，线粒体基质GO证据为IEA:TreeGrafter。序列身份已核对，但原生功能与胞质/线粒体定位的直接实验证据仍未取得。已有AlphaFold预测不填补酶活、方向或定位证据；本审计只核对其API序列，不重新下载或审核全部结构文件。

## 覆盖和访问限制

- 2014结果段与1990正文段通过原论文网页搜索索引独立读取；完整页面直开分别遇reCAPTCHA/抓取失败。没有下载并审核完整论文，不声称复核逆向完整配方、绝对检出限或15000倍的原始计算。
- 未独立扩展审核BIOCHEMISTRY中的其他酶学来源；1987和酵母2014附加证据不计入上述4项覆盖。方向候选由本次核实的核心来源支持。
- 未取得W29条件特异逆向酶学、细胞内浓度/pH或适用ΔG；未新增求解，不能断言候选方向会令任何目标必需。既有机制审计不在本轮重做。
- 核验前后模型SHA-256均为 `10a3baa1da5eb8afb80f93c4cf8a43c3613af7ebb92c8e3385dcda44b91695d0`。此审计仅写本文件，未改科学输入。
