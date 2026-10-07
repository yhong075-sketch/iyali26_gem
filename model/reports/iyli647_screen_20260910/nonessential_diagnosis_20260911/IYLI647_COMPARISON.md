# iYLI647 对这 11 项表现出不同预测的原因

核验日期：2026-09-11。本次为保存来源的静态追踪、版本比较与代数核验，**新增求解 0**。比较对象仍为前轮固定提交中的 `iYLI647_corr_3.json` 与指定 `model_metadata_trna.xml`；使用原 `mapped28_po1f` 运行配置解释保存的预测，不重命名模型、替换输入或接受模型修改。

## 先区分修订历史与当前结构

11 个目标直接关联的 **23 条反应**在仓库最早可见的 `iYLI647_corr.json` 中已经存在。在 `corr → corr_2 → corr_3` 三版之间，这些反应的关联基因集合、布尔 GPR 逻辑、计量和边界均未发生实质变化；部分仅删除 GPR 外层括号。因此不能把这 11 个命中写成该论文逐项修复 iYali26 问题的成果，也不能把最早公开的 corr 当成原始 iYLI647 的全部构建历史。

作者的仓库围绕 13C-MFA、培养数据和生物量表示展开；我们之前的 235 正例比较没有施加论文的 13C 通量约束。这 11 项差异来自现有网络、GPR、区室连接及生物量需求。下面解释其数学机制，不确认哪侧的原生酶功能或 GPR 正确。

所有目标及规则成员的原生正式基因名本次均未核实；下表蛋白功能是候选，证据为 **model/GPR assignment only**。

| iYali26 ID／对应外部 ID | 候选功能 | iYLI647 中的实际反应角色与差异 |
|---|---|---|
| YALI1A02775g／YALI0A02310g | UDP-葡萄糖合成 | GALU 被 KO 关闭。外部存在 GALT（本地 R319 对应），但 UGLT（本地 R780 对应）为单向 [0,1000]，不能反向与 GALT 构成本地两步补偿。当前 galactose／melibiose 不摄取，其余潜在 galactose 来源有结构死端，GALT＋UDPG4E 也不能提供净补偿。不能说作者删除了 GALT。 |
| YALI1B20462g／YALI0B15598g | 6-磷酸葡萄糖酸脱氢酶 | GND 关闭后，6PG 不能继续消耗，稳态令上游 ER G6PDH2er 停止。它是外部 ER NADPH 唯一生产反应，故甾醇合成被切断。外部有非氧化 PPP，胞质有其他 NADPH 来源，但没有本地 R1264/R1265 那样的 NADP/NADPH 胞质↔ER 转运表示。 |
| YALI1C07638g／YALI0C05951g | 脂肪酰-CoA 去饱和酶 | DESAT16/DESAT18 随 KO 关闭。外部也有 FAS161ACPm 等线粒体不饱和脂酰-ACP 路线，然而线粒体 malonyl-CoA 缺乏净供应，稳态迫使这条合成支路为零，不能像 iYali26 那样补偿。属于网络供给断点，不是已验证删除了不合理旁路。 |
| YALI1C08702g／YALI0C06490g | GDP-甘露糖合成 | MAN1PT 采用 AND，而本地 R540 的 OR 保留目标 KO 后的反应。外部 KO 关闭 GDP-mannose 的唯一生产反应，影响 mannan 需求。伙伴身份跨版本对应仍有既有未决项，不能将两侧简化成完全相同的一对蛋白。 |
| YALI1C15991g／YALI0C11407g | 乙酰-CoA 羧化相关 | ACCOACr 是胞质 malonyl-CoA 的唯一净生产来源；MCOATA 因胞质 malonyl-ACP 仅连接自身而被稳态迫为零。外部没有本地 R2004／R2121 的生成列；此外目标又被 AND 写入 11 条脂肪酸／脂酰-CoA 合成反应，KO 同时关闭它们。不能直接据此认可把羧化酶作为所有这些合成反应的必需组分。 |
| YALI1C32184g／YALI0C23364g | 蛋白 O-甘露糖基化 | DOLPMMer 采用三个成员 AND，目标 KO 关闭反应；不像本地 R289 可由其他 OR 分支保留。 |
| YALI1D03865g／YALI0D03069g | GAR 甲酰转移酶 | GARFTi 为 FGAR 的唯一生产反应，该物种全模型只连接 GARFTi 和不可逆 PRFGS。没有本地 R1892 的甲酸＋ATP 依赖生成路线。KO 切断此嘌呤新生步骤；保存预测仍有 3.20% 生长，不应写成零生长。 |
| YALI1E22736g／YALI0E18964g | 溶血磷脂酸酰基转移酶 | AGAT_SC 使用目标单基因 GPR，没有本地 R1846 的 OR 替代成员；因此目标 KO 能关闭该合成步骤。此处不宣称全模型不存在磷脂循环反应。 |
| YALI1E25018g／YALI0E21021g | β-1,3-葡聚糖合酶 | 13GS 采用两个成员 AND，目标 KO 关闭反应；biomass_C 有正的 β-葡聚糖需求，本地 R4 则由 OR 保留。 |
| YALI1F00821g／YALI0F00506g | 谷氨酰胺合成酶 | GLNS 采用两个成员 AND，目标 KO 关闭它，且没有本地 R334/R2081 的重复合成反应补偿。培养配置没有 glutamine 摄取，胞质其他 glutamine 关联反应为消耗或转运。 |
| YALI1F03803g／YALI0F02497g | 功能未决；不同通路的模型赋值 | 外部 MCITDm 虽名 2-methylcitrate dehydratase，实际处理的是 homocitrate，并接入赖氨酸新生合成。它与本地 R95 的甲基柠檬酸计量不是同一反应。KO 后生长恰由外加赖氨酸上限决定，不能解释成外部更严格要求甲基柠檬酸清除。 |

## 外部的完整 GPR 规则

下列成员名称均未核实，功能仅为该行酶步骤的模型候选。AND 的存在不自动证明这些蛋白是已验证必需复合物。

| 外部反应／候选功能 | 完整 GPR |
|---|---|
| MAN1PT；GDP-甘露糖合成 | YALI0E15125g AND YALI0C06490g |
| DOLPMMer；蛋白 O-甘露糖基化 | YALI0C23364g AND YALI0E05929g AND YALI0E15081g |
| 13GS；β-1,3-葡聚糖合成 | YALI0C01411g AND YALI0E21021g |
| GLNS；谷氨酰胺合成 | YALI0D13024g AND YALI0F00506g |
| AGAT_SC；溶血磷脂酸酰基转移 | YALI0E18964g |

## 两项不能仅凭反应名称理解的结果

**GND：区室 NADPH 与甾醇耦联。** 直接对外部模型稳态方程加权求和，得到：

`v_GND = v_G6PDH2 + v_C24STRer + v_SQLEr`。

右侧三个反应均不可逆非负。再结合无外源 ergosterol 摄取及甾醇总量守恒，可得：

`v_GND ≥ v_C24STRer ≥ 0.035250723 × μ`。

因此 GND=0 必令生长 μ=0。这是当前模型和培养条件下的代数结果；不是说非氧化 PPP 在外部不存在，也不证明外部 ER NADPH 表示已获生物学验证。

**MCITDm：实际是 homocitrate／赖氨酸链。** 反应底物 `hcit[m]` 经 `b124tc[m]`、`hicit[m]` 连到赖氨酸新生合成，关联 subsystem 本身就是 Threonine and Lysine Metabolism。这里的命中取决于另一条模型功能赋值。外加 lysine 上限 0.02466，biomass_C 中系数 0.275004895，因此：

`μ ≤ 0.02466 / 0.275004895 = 0.08967113112659321`。

除以前轮外部 WT 1.2171871273706656，等于 **7.367078496820643%**，与保存 KO 比值一致。该数值解释不能移植成本地 R95 的甲基柠檬酸生物学需求。

## 判断与可复查证据

可作为审查线索的是 GPR 逻辑、反应方向、重复反应和区室供给耦联。当前资料不足以把外部 AND、断开的线粒体支路或不同功能赋值直接接受到 iYali26。命中更多正例并不逐项确认这些表示的科学正确性。

- [固定来源、完整 SHA、11 个目标／23 个反应比较及完整代谢物连接](external_comparison.json)
- [版本历史独立核验](EXTERNAL_HISTORY_AUDIT.md)
- [六项结构差异、代数证书及来源定位的独立核验](EXTERNAL_STRUCTURAL_AUDIT.md)
- [原 iYali26 22 次定向求解与控制](REPORT.md)
- [固定提交的外部最终模型](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_3.json)
- [作者项目说明及 notebook 定位](https://github.com/UH-MBBE/yarrowia-13C-gsm/tree/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6)

审计覆盖分开计数：11/11 目标及其 23/23 直接关联反应的版本不变性；6/6 结构机制；5/5 目标相关 GPR 表示核对。以上是代码、模型和保存结果的审计，不是原生蛋白功能的实验验证。本次新增求解0，没有改写模型、GPR、培养或实验标签。
