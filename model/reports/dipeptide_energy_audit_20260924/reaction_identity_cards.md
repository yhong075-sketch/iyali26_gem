# 四条游离二肽水解：反应身份与证据卡

核验日期：2026-09-24。结论针对固定的本地 R1159 默认输出；模型赋值、同源注释、结构预测与原生实验分开。原始字段和完整 SHA 见 `provenance/local_xml_facts.json`、`before.json`；逐来源审查见两个子报告。

**本轮四条反应均不接受新 GPR，也不接受正式区室迁移。** “没有取得精确底物实验”表示未决，不表示该酶无活性。保留现有反应作来源记录，也不等于接受其液泡定位。

## 共同来源和化学限制

四条反应在 `data/iyli21.xml`、`data/iyali26.xml` 已经具有相同的空 GPR、[0,1000] 边界以及“cytosol nonspecific dipeptidase”名称，但全部底物/产物实际属于 C_va。当前 metadata 字段选择没有这些条目，不能把名称/区室冲突归因于本次 R1159 修改。没有在授权工作区找到 `model(1).xml`，未向外寻找。原始反应条目未提供可接受的液泡定位依据。

12 个同名二肽池（每种的 ex/cy/va）全部缺少分子式和结构 ID；液泡 Gly 也缺式。名称只给出拟表达的残基顺序，不能独立确认立体化学与连接方式。同名跨区室连接需要结构身份核验。Gly 本身无手性，旧名“L-glycine”不构成另一种立体异构体。

当前 Asp/Glu 产品分别使用 C4H7NO4/0、C5H9NO4/0。如果确认普通 α-肽键、L 型非 Gly 残基以及这些中性形式，则可由产物之和减 H2O 得到下面的**条件性算术候选**；尚未填入模型：

| 化合物 | 分子式候选 | 电荷候选 |
|---|---|---:|
| Gly | C2H5NO2 | 0 |
| Gly-Asp | C6H10N2O5 | 0 |
| Gly-Glu | C7H12N2O5 | 0 |
| Ala-Gly | C5H10N2O3 | 0 |
| Gly-Pro | C7H12N2O3 | 0 |

旧 metadata 将某些中性酸分子式与 −1 电荷混用；Yeast-GEM 的酸性二肽和氨基酸采用去质子形式。不能跨版本只复制电荷或只复制分子式来宣称配平。

## R2021 — Gly-Asp

- **当前完整化学：** m1871[C_va] Gly-Asp + H2O[va] → Gly[va] + L-Asp[va]；所有系数 1；[0,1000]，GPR 空。
- **底物身份：** N 端 Gly、C 端 Asp；不能用 Asp-Gly、较长肽或显色底物的实验替代。当前名称不能确认 Asp 手性或肽键结构。
- **候选酶裁决：** 尚无可接受的 W29 精确底物酶。下述 YALI1F23706g 只作弱筛选线索，不能赋予该反应。
- **区室裁决：** 无已取得的 W29 液泡腔、胞质或分泌到胞外的精确反应证据。胞质与液泡假设需分别比较，不能同时复制反应制造供给。
- **最低确认：** 版本明确的 W29 蛋白对游离 Gly-L-Asp 生成两种氨基酸的定量实验；Asp-Gly/长肽对照、相关 pH/金属条件、阴性对照；随后独立定位、底物可达性及缺失/恢复实验。

## R2029 — Gly-Glu

- **当前完整化学：** m1862[C_va] Gly_Glu + H2O[va] → Gly[va] + L-Glu[va]；所有系数 1；[0,1000]，GPR 空。
- **底物身份：** 拟表达 Gly-Glu，不能与 Glu-Gly 或 γ-Glu-Gly 混同。谷胱甘肽降解/γ-谷氨酰酶的证据不自动支持此 α-二肽。
- **候选酶裁决：** 无可接受的 W29 精确底物酶；YALI1F23706g 对它的活性未由所读实验证明。
- **来源补充：** 本地 Yeast-GEM r_4439 为液泡 Gly-Glu 水解，空 GPR、confidence 1，notes 指向 Biolog PR149/Rhea36463。它是参考模型线索，不能据此证明 iYali 的复制历史或原生酶/定位。
- **最低确认：** 直接测试游离 Gly-L-Glu，并用 Glu-Gly、γ-Glu-Gly 区分连接与顺序；定位、pH、产物和遗传确认要求同上。

## R2034 — Ala-Gly

- **当前完整化学：** m1878[C_va] Ala-Gly + H2O[va] → L-Ala[va] + Gly[va]；所有系数 1；[0,1000]，GPR 空。
- **底物身份：** N 端 Ala、C 端 Gly；Gly-Ala 不能代替。L-Ala 产品提示建模意图，但底物尚无立体结构 ID。
- **候选酶裁决：** 无可接受的 W29 精确底物酶。DUG1 同源性不足以接受；Aspergillus Cys-Gly 酶研究对 Ala-Gly 未检出活性，只构成跨物种反对线索，不是 W29 阴性结论。
- **区室裁决：** 仍未知；不能用“胞内”直接指定胞质，更不能直接指定液泡腔。
- **最低确认：** 游离 L-Ala-Gly → 两种氨基酸的定量实验，以 Gly-Ala 和长肽作区别对照；定位、pH、金属与遗传证据同上。

## R2039 — Gly-Pro

- **当前完整化学：** m1866[C_va] Gly-Pro + H2O[va] → Gly[va] + L-Pro[va]；所有系数 1；[0,1000]，GPR 空。
- **底物身份：** Gly-Pro 属 Xaa-Pro，不能与 Pro-Gly 混同；当前未独立编码底物手性/cis-trans 状态。
- **优先实验候选：** **YALI1E16433g；原生正式名称未核实；预测 M24B prolidase-like 肽酶。** W29 AOW05368.1 = XP_503902.2，454 aa；自动 CDD 注释含 Prolidase cd01087（159–427）及 AMP_N。支持家族候选，尚未证明游离 Gly-Pro 催化或原生定位。
- **菌株/版本：** CLIB122 对应 YALI0E13464g/CAG79495.1/Q6C610（同源 prolidase-like 候选，原生名未核实）同为454 aa，但第17位为 K，W29 为 Q。不能把旧株序列或结构称作与 W29 完全一致；完整序列 SHA 见 `provenance/current_sequence_verification.json`。
- **功能边界：** [EC3.4.13.9](https://iubmb.qmul.ac.uk/enzyme/EC3/4/13/9.html) 对应游离 Xaa-Pro；[DPP-IV/EC3.4.14.5](https://iubmb.qmul.ac.uk/enzyme/EC3/4/14/5.html) 是从更长肽切下二肽。Gly-Pro-pNA 切割不等于游离 Gly-Pro 的 Gly–Pro 键水解。[2023年真菌 XPD 实验](https://www.mdpi.com/2311-5637/9/11/978)支持该酶类和阳性对照，不能移植其物种定位或速率。
- **最低确认：** 精确 W29 蛋白对未修饰 Gly-L-Pro 的双产物定量；Pro-Gly、含 Gly-Pro 的长肽及 pNA 底物分别对照；相关 pH/金属、热灭活/催化位点对照；原生定位与纯度标记、缺失/恢复及底物/产物改变。

## R2200 线索为何不足以支持四条GPR赋值

**YALI1F23706g；原生正式名称未核实；DUG1-like/M20A Cys-Gly 金属二肽酶候选。** 当前 W29/CLIB89 XP_505554.3 = AOW07335.1，478 aa；模型给其胞质 R2200：Cys-Gly + H2O ⇌ Cys + Gly。模型双向边界、GPR和合并EC均不是原生实验确认。

酿酒酵母 **YFR044C/DUG1（实验确认的 Cys-Gly 金属二肽酶）** 的[2009年Table2](https://pmc.ncbi.nlm.nih.gov/articles/PMC2682898/)测试的是含 Cys 的二肽；四个目标均未测试。Ala-Cys 活性也较高，因此不能反向夸大成“只能切Cys-Gly”。[2007年GFP定位](https://pmc.ncbi.nlm.nih.gov/articles/PMC1840075/)支持该同源酶胞质分布；W29 的胞质 GO 为 IEA 自动注释，不能证明液泡腔定位。[Aspergillus研究](https://doi.org/10.1271/bbb.100604)对 Ala-Gly/Gly-Ala 的阴性结果是跨物种证据；完整PDF未成功取得，不夸大读取范围。

复用**既有 AlphaFold 预测** AF-A0A1D8NNX0-F1-model_v6，478 aa 与 W29 当前序列完全匹配；未新跑预测、对接或功能实验。AFDB报告 pLDDT 96.56，坐标 CA 均值96.593；PAE均值3.805Å。高结构置信度不能确认这四种底物、液泡定位或寡聚状态。序列来源/完整SHA、AF版本、模型创建日、取得日与置信数据见 `identity_r2200/sequence_identity.json`。

当前 KEGG 旧序列为 XP_505554.2、477 aa，来源 DSM3286/YALI2；不是 W29 当前 XP_505554.3 的同一菌株版本。CLIB122 Q6C1A8 与当前 W29 的478 aa在本次核验完全一致，也不能据此移植定位和条件。

## 证据门槛

将“蛋白家族”“精确自由底物催化”“原生区室”“实际培养有底物”“转运机制”分别验收。缺第一/第二项时不补 GPR；缺区室时不搬反应；缺供应时不注入无限二肽；缺原生生长需求时不据此声称必需基因。来源独立审核见 `provenance/R2200_SOURCE_AUDIT.md` 与 `identity_r2200/PROLIDASE_SOURCE_AUDIT.md`。
