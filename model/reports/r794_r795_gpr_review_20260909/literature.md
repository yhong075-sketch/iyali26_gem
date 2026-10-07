# R794/R795 替代分支：原始文献与数据库限定审查

采集日期：2026-09-09。范围按本目录 TASK.md；仅审查文献/数据库，不运行 GEM、不修改 GPR。本文是来源采集与推理，独立来源审核见总审查产物，不能把采集人的核读计作独立审核。

**结论：现有来源支持两个末尾位点是 V-ATPase 膜内部分的组件候选，不支持它们组成一台可独立替代整个 V1/V0 复合体的 ATP 驱动质子泵。** 最值得审核的修订方向是恢复共同必需的复合体组件关系，并只在证据支持的同类亚基位置放置 OR。现有证据尚不足以直接给出最终14基因布尔式或声称14个均应生长必需。

## 两个核心位点

| 系统 ID | 已核实名称状态 | 最佳功能候选与证据等级 | 当前模型角色 | 关键限制 |
|---|---|---|---|---|
| YALI1F38820g | 未建立原生正式名称；不能直接命名为 VPH1/STV1 | A0A1H6PMT1，804 aa；V-ATPase V0 亚基 a 候选，参与膜内质子转运与复合体组装。**uncharacterized；历史自动同源/家族注释** | 与下行位点 AND 后，整体 OR 于前12位点组 | 固定 UniProt 快照为 TrEMBL、protein existence 3。旧位点 YALI0F31119g 的 GRYC 注释仅为酿酒酵母 VPH1 同源，不能自行确认跨菌株映射或区室 |
| YALI1F13017g | 未建立原生正式名称；PPA1/VMA16 是待验证同源命名 | A0A1D8NMN8，196 aa；V-ATPase V0 蛋白脂质 c-like、倾向 c″ 候选，参与膜内质子转运。**uncharacterized；历史自动同源/家族注释** | 同上 | 固定快照只明确 c-like。旧位点 YALI0F09405g 在 GRYC 与原始转录组中对应酿酒酵母 PPA1/VMA16 同源；这不是原生酶学鉴定 |

以上是**历史固定快照**中的记录，不宣称对应 UniProt 当前条目仍处于 active 状态。当前状态、目标序列身份、BLAST 与 AlphaFold 由根任务另行核查；不得静默换为当前新序列。

## 对 GPR 真正有约束力的证据

**L01—目标注释落在 V0，而不是完整酶。** 本地 UP000182444 快照按 `genes.orfNames` 去下划线精确匹配得到上述 accession。亚基 a 记录的功能描述是多亚基泵的组件；proteolipid 记录区分负责 ATP 水解的 V1 和负责跨膜 H+ 的 V0。它们都不是“两个蛋白即完整 ATPase”的注释。证据类型：历史自动数据库注释；不能独立证明原生生化功能。原始文件身份、条目/序列版本和 SHA 见 literature_sources.json。

**L02—本物种的名称依据仍是同源，不是原生功能实验。** [GRYC 旧位点 F31119](https://gryc.inrae.fr/db/yarrowia-lipolytica/clib-122/yali0f/yali0f31119g) 的 Overview 指向 P32563 同源；[GRYC 旧位点 F09405](https://gryc.inrae.fr/db/yarrowia-lipolytica/clib-122/yali0f/yali0f09405g) 指向 P23968 同源。[Morin 等，2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3222671/) Table 2 将后者写作 PPA1/VMA16、c″，但表下注明名称来自酿酒酵母同源，Methods 也写明 GO/同源分类。该实验是 Y. lipolytica JMY1346、28°C、pH 5.8、葡萄糖补料转向氮限制的转录组；它提供表达观察，不是泵重构、定位验证或敲除互补。旧 ID 与当前 YALI1 ID 的映射仍须独立序列核查。

**L03—真正的亚型替代仍需要共享其他亚基。** [Vasanthakumar 等，2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6462096/)，Results “ATPase Activity Assays”、Fig. 1、 “Structure of Stv1-V1VO”，从删除两个内源 a 亚型的酿酒酵母背景分别提供一种 a 亚型并纯化**完整复合体**，测得两者均有 ATPase 活性。正文明确二者只在 a 亚型不同，V1 催化部分及其余 V0 组件仍存在。Fig. 1 数据为两个生物重复各三次测量，不能把 n=6 写成六个独立生物重复。证据类型：其他物种直接生化/结构实验；向 Yarrowia 迁移是间接推断。

**L04—a 亚型确可部分互补，但原始实验直接显示其与其他亚基共同工作。** [Manolson 等，1994](https://pubmed.ncbi.nlm.nih.gov/7514599/)，Abstract：高拷贝另一 a 亚型部分恢复缺失液泡 a 亚型株的酸化；其与 60/69 kDa V-ATPase 亚基及药物敏感 ATPase 活性共同纯化。双 a 亚型缺失才在所测高金属、pH 7.5、甘油培养条件下复制其他核心亚基单缺失的生长表型。这里可支持的结构是“共享组件 AND 合适的 a 亚型”，不能推出“两组完全独立的酶”。现只核读原始摘要，未核对全文全部重复数、统计和构建细节。

**L05—蛋白脂质亚基并不是随便互换的同工酶。** [Hirata 等，1997](https://pubmed.ncbi.nlm.nih.gov/9030535/)，Abstract，经 NCBI efetch 核读：酿酒酵母 VMA16 缺失导致 V-ATPase 活性与装配丧失；标签定位与共纯化支持其为复合体亚基，突变结果支持三个 proteolipid 相似但非冗余。它支持 c″ 候选应作为组件调查，不能把它与一个 a 亚基联合后当作全泵。此结论在酿酒酵母成立；原生 Yarrowia 对应位点尚无本轮找到的直接功能验证。尚未完整核读1997全文，不报告其具体生长数值或重复数。

## 反证与不能越过的限制

**L06—一般局限：“复合体组件”不自动等于每个催化条件或生长条件都必需。** [Bueler 与 Rubinstein，2015](https://pubmed.ncbi.nlm.nih.gov/25546637/)，Abstract，经 NCBI efetch 核读：酿酒酵母亚基 e 可在 DDM 处理时离开已装配复合体，剩余纯化酶仍能进行耦联 ATP 驱动泵 H+。这是一般机制限制；**e/Vma9 不在本题14位点内，不能用它作为本题14位点全 AND 的直接反证。** 也不能反向推出亚基 e 在活细胞装配或 Yarrowia 生长中不重要。

**本题暂不能全 AND 的主要理由**是：潜在 a 亚型可能填充同一个可替代的功能位置，而不是彼此都必需；全部原生位点的功能及高尔基体/液泡归属还未闭合。1994和2019原始研究提供这种亚型替代的真实范例，却不直接决定本题哪个位点在哪个区室可替代。当前外层 OR 缺证据，不代表全 AND 已获证据。

**L07—本物种确有液泡 ATP 驱动酸化的生理证据，但无法从中识别这两个基因。** [Kulakovskaya 等，1993](https://pubmed.ncbi.nlm.nih.gov/8385405/)，Abstract，经 NCBI efetch 核读：来自葡萄糖培养、不同碳/氮限制生长阶段的 Yarrowia 离体液泡具有 H+-ATPase 活性、膜电位与质子梯度，且随生长阶段变化。该工作没有给出这两个位点身份、亚基最小组合或基因必需性，不能把液泡总活性当作该二基因分支的验证。

本轮检索覆盖两个 YALI1 ID、桥接旧 ID、Yarrowia 与 VPH1/STV1/VMA16 关键词、真核 V-ATPase 亚型与 proteolipid 原始实验。未找到能证明二蛋白可独立完成 ATP 水解并泵 H+ 的原生实验。**这是限定检索未找到证据，不是证明任何此类原生机制绝对不存在。** 本轮未把结构相似当作替代性实验，也未把“essential component”自动译为“生长必需基因”。

## 可供结构核对的实验参照

2019研究对应 [6O7T](https://www.rcsb.org/structure/6O7T) 为液泡 V0，3.20 Å；[6O7U](https://www.rcsb.org/structure/6O7U) 为 Golgi V0，3.10 Å。两者为实验冷冻电镜，**不是 AlphaFold**。RCSB Macromolecules 中：a 亚基分别为 P32563/P37296，label A、auth a、entity 1；c″ 均为 P23968，label C、auth c、entity 3。这些孤立 V0 结构不含 V1，不代表孤立 V0 自身具备 ATP 水解泵能力。其余2019完整 Stv1复合体结构为6O7V–6O7X，本文未逐链核验。

酿酒酵母参照名称与身份：YOR270C—VPH1—液泡 V0 亚基 a；YMR054W—STV1—Golgi/内体 V0 亚基 a；YHR026W—VMA16/PPA1—V0 蛋白脂质亚基 c″（前三者在所引原始研究中具实验支持）；YCL005W-A—VMA9—V0 亚基 e（2015体外去除试验）。这些名称不能不经身份核验直接搬到 Yarrowia。

## 建议判定

当前二基因分支作为**完整替代泵**：`unsupported`，并与已知 V1/V0 分工存在机制冲突；尚未原生实验反驳到可称“该物种绝不可能”的程度。

两位点作为**同一泵的组件候选**：文献/历史数据库支持该研究方向，仍需根任务的序列、家族/结构、区室和全14位点检查。最终 GPR 必须由每个必需功能位置与可互换亚型的证据逐层组装；泵失活是否导致当前培养条件下的生长必需性是后续独立问题。

独立审核覆盖：本文7个证据主张均待总审独立判定；采集人核读不计入独立审计。报告不将这些条目单方面标作 fully verified。

源快照说明：最初读取的完整响应未落盘。为独立审计，于2026-09-09 22:57 UTC对相同已读来源重新读取一次，保存完整 HTML/XML 至 `sources/`；这些是新快照，不冒充首次响应。literature_sources.json 分别保留首次响应SHA和新快照路径/SHA/时间。重新取得的NCBI三文摘要XML与首次SHA一致；网页HTML字节不同不自动说明正文结论改变，独立审计以新快照内容为准。
