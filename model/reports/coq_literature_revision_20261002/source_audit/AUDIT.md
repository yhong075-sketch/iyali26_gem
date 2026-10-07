# R19—R18—R695—R385 氧化态独立来源审核

审核日期：2026-10-02。身份：独立来源审核；未修改模型或科学代码，未运行求解、序列比对或预测。范围遵循上级 `SCOPE.md`。固定候选为 `artifacts/coq_yeast_reference_20261001/coq_c5_coq9_candidate_v2.xml`，SHA256 `9e09e96e92a95b6ba9c079b00d9bee670a69f85dea32ab33dbffccb7fa03ccae`。

## 结论

支持实施 **R19 生成 DDMQ9H2，R18 使用 DDMQ9H2 并生成 DMQ9H2** 的有限候选。没有找到足以定义 DMQ9H2 后续完整酶促连接的原始证据。保留原 R695 时必须明确 DMQ9H2 是结构断点；不能把通路已闭合或表型已验证作为结论。此处缺口来自底物氧化态及电子受体未定义，不是笼统要求原生酵母实验才能进行同源 GPR 推断。

原 R385 同样把醌作为甲基化底物，而官方反应及跨物种酶学支持醌醇甲基化。保留 R385 是本轮范围限制，不能称其化学已经核实。若将来完整重建，还须明确 COQ7 产物醌如何还原，或者证明该酶组的其他净化学，不应把未知氧化从上游平移到下游。

## 关键来源与证据强弱

| 来源与精确定位 | 支持 | 反证或限制 | 判定 |
|---|---|---|---|
| [Behan & Lippard 2010](https://doi.org/10.1021/bi101475z)，正文关于 reductants/DMQ 的连续五段；本地 `Behan2010.txt` 263–306 | 纯化小鼠 CLK-1，以 NADH 和氧化态 DMQ0 产生加一氧的产物；18O2 标记支持 O2 提供氧。 | p-hydroquinone 作 ubiquinol mimic 时未观察到供体氧化；这不是 DMQ9H2 特异阴性。未直接测外源 DMQH2 的完整羟化。该文早期倾向直接 NADH→铁中心，后续机制被2013/2022研究细化。 | 支持 NADH/醌型羟化；不支持任意醌醇直接供电子即可闭合。 |
| [Lu et al. 2013](https://doi.org/10.1021/bi301674p)，Results “Hydroxylation of DMQ derivatives”、Fig.3、Scheme1、Results/Discussion quinone-mediated reduction；本地 `Lu2013.txt` 531–599、880–971 | 人源 GB1-CLK-1 + DMQ0 + NADH 得羟化产物；DMQ2 停流实验及光谱支持结合态醌在 NADH 与双铁之间传递电子。 | DMQred 是作者据动力学推断的瞬态；未找到“直接加入 DMQH2、不加 NADH，测得完整羟化产物”的实验。双铁被还原不等于完整通量反应及产物氧化态被确定。 | 可支持催化机制，不能把自由 DMQH2 池直接替换现有 R695 底物。 |
| [Manicki et al. 2022](https://doi.org/10.1016/j.molcel.2022.10.003)，Fig.2E–I/S9G–M，Results NADH-binding，Methods “Measurement of DMQ0’s redox state by HPLC-ECD”及“NADH fluorescence”；本地 `Manicki2022.txt` 638–684、1216–1232、1245–1249 | 纯化人源 COQ7 的 NADH 活性将 DMQ0 还原；HPLC-ECD 直接测受体氧化态。 | 实际输入为 NADH + 氧化态 DMQ0，不是 DMQH2→羟化产物；上游/其他氧化还原酶可能使 NADH 功能可免是作者解释，不能升级为完整外源 quinol 羟化实测。原文这里写 DMQ6，未明确标注 H2，不能静默补入。 | 提高 quinone-mediated 机制可信度；仍不闭合模型的 DMQ9H2 连接。 |
| [Nicoll et al. 2024](https://doi.org/10.1038/s41929-023-01087-z)，Fig.1b、Fig.5、Extended Data Fig.6、Results C2/C6 transformations、Methods Small-scale reactions；复用已归档全文 `artifacts/coq_reaction_first_20261001/source_audit/sources/Nicoll2024.txt` 705–856、1233–1254 | COQ5 + SAM 将4b醌醇变为5醌醇；COQ7/9/3 + NADH再生系统可从5ox得到最终 CoQ1 和 CoQ1H2。 | COQ7 的 NADH 消耗仅在5ox而非5中观察到。这一读数不是排除所有 quinol 驱动反应的通用证明；联合末端实验也未逐一确定氧化还原交换、对应基因及净辅因子计量。 | 支持 R18 醌醇版本；联合重构证明存在体外连接，但不能据此编造每个模型反应的电子收支。 |

四篇均为原始研究。Nicoll 使用重构的四足动物祖先蛋白和短链底物；Lu/Behan/Manicki 使用动物酶及短链/无侧链模拟底物。它们支持跨物种候选化学，不是原生 Yarrowia 生化验证。

官方记录另提供明确反应边界：[EC1.14.13.253](https://iubmb.qmul.ac.uk/enzyme/EC1/14/13/253.html) 指定真核 COQ7 的醌 + NADH + O2→羟基醌 + NAD+ + H2O，并称已测脊椎动物酶对醌醇无活性；[EC1.14.99.60](https://iubmb.qmul.ac.uk/enzyme/EC1/14/99/60.html) 已于2024修订为原核醌醇步骤。旧 EC 或数据库概述中的 DMQH2 不能替代实际酶学。[EC2.1.1.201](https://iubmb.qmul.ac.uk/enzyme/EC2/1/1/201.html) 及 [EC2.1.1.64](https://iubmb.qmul.ac.uk/enzyme/EC2/1/1/64.html) 分别给出 COQ5 和末端 COQ3 的醌醇甲基化。官方概括的“无活性”须保留实际测定类型的限制。

## 可平衡部分与无法定义部分

以下均按固定模型现有中性 NADH/NAD、SAM/SAH 的分子式表达，不擅自换成标准带电辅因子形式。NADH 比 NAD 多2H；SAM 比 SAH 多 C1H2。

| 步骤 | 计量表达 | 实施等级 |
|---|---|---|
| R19 候选 | m63 `C52H80O2` + O2 + NADH → DDMQ9H2 `C52H80O3` + NAD + H2O | 原子/电荷可平衡；COQ6/还原供体跨物种候选。NADH 是本模型兼容的候选供体，非原生专一性声明。 |
| R18 候选 | DDMQ9H2 `C52H80O3` + SAM → DMQ9H2 `C53H82O3` + SAH | 原子/电荷可平衡；保留 COQ5 候选 GPR。 |
| 未定义氧化 | DMQ9H2 `C53H82O3` → DMQ9 `C53H80O3` + 2H+ + 2e− | 仅是电子账本的半反应；没有指定受体/机制，**不得作为模型反应加入**。未能指认基因。 |
| R695 原反应 | DMQ9 + NADH + O2 → DMeQ9 `C53H80O4` + NAD + H2O | 原子/电荷可平衡，符合已测醌型 COQ7；不能接收上行产生的 DMQ9H2。 |
| 未定义还原 | DMeQ9 + 2H+ + 2e− → DMeQ9H2 `C53H82O4` | 仅为半反应账本。没有本次审核支持的独立供体、催化酶或净通量定义，不加入模型。 |
| COQ3 文献形式 | DMeQ9H2 + SAM → CoQ9H2 `C54H84O4` + SAH | 可平衡且有甲基化化学支持；不能绕过上一行的未定义供给而宣称闭合。 |

因此，本次不能提供“各酶、电子受体和计量均被证实”的完整闭合方案。尤其 `DMQH2 + NAD → DMQ + NADH` 虽能在模型中原子守恒，但没有相应直接酶学证据；加入它会凭空回收还原力。氧气自氧化也必须明确产物、活性氧处理及生理意义，不能从 Fig.1 的一般可能性直接造一个可用通量。成熟 CoQ/CoQH2 呼吸酶的底物特异性不移植至 DMQ/DDMQ/DMeQ。

## 基因候选与证据复用

本轮不重新计算身份、结构或定位；复用先前版本化序列/同源及已有同序列 AlphaFold 预测证据。下表所有名称均是本项目候选/家族名称，原生特定底物催化未经本次实验验证。结构预测不解决本次氧化态缺口。

| 系统 ID / 版本化蛋白 | 候选名称及功能 | 对本轮的作用 |
|---|---|---|
| YALI1A08781g / XP_499891.1 | COQ6；FAD 羟化酶候选 | R19 C1羟化候选；不是未知独立醌醇氧化酶。 |
| YALI1B03314g / XP_500417.1 | YAH1；铁氧还蛋白候选 | 与还原酶共同支持 COQ6 电子供给假设。 |
| YALI1B19490g / XP_500902.3 | ARH1；铁氧还蛋白还原酶候选 | 原生 NADH 专一性未确立；保留耦联候选的证据限制。 |
| YALI1C25352g / XP_501966.1 | COQ5；CoQ 环 C-甲基转移酶候选 | 适于 R18 醌醇→醌醇形式。 |
| YALI1E18269g / XP_503973.1 | COQ7；双铁 DMQ 羟化酶候选 | 保留 R695 醌型底物，不升级为 DMQH2 独立供体/底物验证。 |
| YALI1F34675g / XP_505941.3 | COQ9；脂质结合和底物呈递辅助蛋白候选 | 保留 COQ7 功能依赖候选，不能填成独立氧化酶。 |
| YALI1B20835g / XP_500950.3 | COQ3；CoQ 环 O-甲基转移酶候选 | 末端醌醇甲基化；原 R385 的醌底物保持待审限制。 |

身份来源复用 `artifacts/coq_yeast_reference_20261001/coq_c5_coq9_candidate_v2.coq9_genes.tsv`、此前同目录和 `coq_reaction_first_20261001` 的来源审核。COQ5 accession 本轮仅核对该固定候选基因表的 RefSeq 字段为 XP_501966.1；不表示重新核验序列或功能。

## 覆盖与停止

以下8条为本轮原子决策声明，每条只进入一个互斥类别。**总数8；已审8；支持4；未决3；反证1；未核0。** 支持表示在表列物种、测定和推断层级内成立；未决表示已查但不足以确定；反证表示该条精确表述与已读来源冲突，不能扩大成所有物种或条件下不可能；未核表示没有完成来源查阅。已审=支持+未决+反证；总数=已审+未核。覆盖率8/8不是验证率。

| Claim | 原子声明 | Verdict | 判定依据与范围 |
|---|---|---|---|
| C1 | Nicoll2024 的 COQ6 C1 羟化产物是 DDMQ 醌醇。 | 支持 | Fig.1b、Fig.4及 small-scale 4a 转化；支持 R19 醌醇产物的跨物种候选，不认证原生 C9 底物或 NADH 专一性。 |
| C2 | Nicoll2024 的 COQ5 甲基化将 DDMQ 醌醇转为 DMQ 醌醇。 | 支持 | Fig.5、Extended Data Fig.6及 small-scale 4b 转化；与 EC2.1.1.201 一致。 |
| C3 | 已测动物 COQ7 能用氧化态 DMQ 与 NADH 完成羟化。 | 支持 | Behan2010 产物与18O2实验；Lu2013 Fig.3；Nicoll2024 Fig.5。实测以无侧链/短链底物为主。 |
| C4 | Lu2013 的动力学支持结合态 DMQred 介导 NADH 到双铁中心的电子传递机制。 | 支持 | Scheme1及 quinone-mediated reduction 段；这是受实验支持的机制推断，非直接捕获完整外源 DMQH2 催化。 |
| C5 | 外源 DMQ9H2 可直接作为 Yarrowia COQ7 的完整羟化底物。 | 未决 | 无本次范围内的直接产物实测；动物酶 NADH 消耗阴性及 hydroquinone mimic 阴性构成限制，但不能独自证明本声明在原生菌中不可能。 |
| C6 | 原 R385 的醌→醌甲基化计量与已审核的 EC2.1.1.64 定义相同。 | 反证 | 该 EC 明确规定 DMeQ 醌醇→CoQ 醌醇；固定候选 R385 两端均为醌。反证针对“同一化学定义”，不是断言原生酶绝不接受醌。 |
| C7 | DMQ9H2 至末端 CoQ 的每一步电子交换均可依据本次来源确定唯一净计量。 | 未决 | 联合 COQ7/9/3 重构测得两种末端氧化态，但未分辨所需氧化/还原连接及各自辅因子收支。 |
| C8 | DMQ9H2→DMQ9 的具体原生催化基因能够从本次来源确定。 | 未决 | 未找到可直接指认该底物氧化的基因；不移植成熟 CoQ 呼吸酶的底物特异性。 |

本轮新取3篇原始研究全文及4条官方记录，复用并重读1篇原始研究关键实验。未找到确定连接是限定来源中的未决，不是断言自然界不存在该反应。

原始检索日期、来源 URL、字节身份见 `sources/retrieval.json` 与 `MANIFEST.json`。科学变更实施由主代理负责；本审核仅支持明确标记断点的有限候选。
