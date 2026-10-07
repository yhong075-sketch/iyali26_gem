# R610：鸟氨酸氨基转移酶的生化与真菌原始文献

核查日期：2026-09-15（America/Los_Angeles）。本子任务仅检索、读取及整理文献；未执行 FBA、未修改模型/GPR。目标序列及当前数据库状态以同目录 `IDENTITY.md` 的实际取回记录为准。本文件的文献证据尚待独立审计，不把相关物种的实验证据写成目标蛋白的原生验证。

## 结论与适用范围

**YALI1C05908g — 原生正式名称未核实 — 推定鸟氨酸 δ-氨基转移酶、PLP 依赖（模型/GPR赋值，加自动同源和结构域注释；未取得本蛋白直接实验验证）。** 本次有界检索未找到以该 W29 蛋白、旧标识 YALI0C04433g 或相应 accession 为对象的原生酶学论文；这不等于证明不存在相关研究。

R610 的自由鸟氨酸 + 2-氧戊二酸 → 谷氨酸-5-半醛 + 谷氨酸符合 EC 2.6.1.13 的反应类型。酿酒酵母单一结构基因的功能互补支持“一种 OAT 蛋白完成该步”这一 GPR 形式，但不直接验证 Yarrowia 靶蛋白，也不排除同源寡聚体。现有材料未提供必须增加另一异源催化亚基的依据；PLP 是小分子辅因子，不能因此新增 AND 基因。

**不能因为该反应提供脯氨酸合成旁路就判它错误。** 酿酒酵母已有条件依赖的遗传证据支持鸟氨酸途径补偿；不过其物种、培养条件与当前模型不同。该证据不证明 Yarrowia 在含铵 SD-Leu 条件下实际表达足够酶量，也不证明模型容量或双向通量范围。

## 核心证据（五组，按作用区分）

### 1. 原始基因互补：酿酒酵母 CAR2

Degols G, Jauniaux JC, Wiame JM. 1987. *Molecular characterization of transposable-element-associated mutations that lead to constitutive L-ornithine aminotransferase expression in Saccharomyces cerevisiae.* Eur J Biochem 165:289–296. DOI [10.1111/j.1432-1033.1987.tb11440.x](https://doi.org/10.1111/j.1432-1033.1987.tb11440.x)；[PubMed 摘要](https://pubmed.ncbi.nlm.nih.gov/3036506/)。

- **YLR438W — CAR2 — 鸟氨酸氨基转移酶（酿酒酵母功能互补实验支持）**：作者用对缺陷突变的功能互补克隆该结构基因，并分析改变其表达的上游插入。
- 定位：摘要第一句为基因互补；后续段落为调控突变。支持功能身份及调控与催化基因的区分。不是 Yarrowia 原生实验，也不是纯化酶的底物比较或亚基组成实验。
- 访问：本次成功直接打开 PubMed 摘要；未取得出版社全文，不声称核查全文方法。

### 2. 胞质定位先例及其跨物种限制

Dougherty KM, Swanson DA, Brody LC, Valle D. 1993. *Expression and processing of human ornithine-δ-aminotransferase in Saccharomyces cerevisiae.* Hum Mol Genet 2:1835–1840. DOI [10.1093/hmg/2.11.1835](https://doi.org/10.1093/hmg/2.11.1835)；[作者机构摘要](https://pure.johnshopkins.edu/en/publications/expression-and-processing-of-human-ornithine-%CE%B4-aminotransferase-i-3/)；[PubMed](https://pubmed.ncbi.nlm.nih.gov/8281144/)。

- 摘要明确区分：酿酒酵母本身的 OAT 位于胞质，而异源表达的人 OAT 到达酵母线粒体基质；后者仍恢复缺陷株利用鸟氨酸作为唯一氮源的能力。
- 支持：胞质 OAT 在真菌中有依据，不能按哺乳动物定位一概指定线粒体。限制：本论文的实验主对象是人酶在酵母中的表达，原生酵母胞质定位是文中陈述；这里没有追溯其最初定位实验。更不能据此认定 Yarrowia 目标定位。
- 访问：核查作者机构公开摘要/搜索索引；PubMed 的部分直接请求出现验证页，出版社全文未成功取得。未将摘要当全文。

### 3. 真实补偿路径的遗传证据与严格培养条件

Shlomi T et al. 2007. *Systematic condition-dependent annotation of metabolic genes.* Genome Res 17:1626–1633. DOI [10.1101/gr.6678707](https://doi.org/10.1101/gr.6678707)；[出版社公开全文 PDF](https://genome.cshlp.org/content/17/11/1626.full.pdf)。

- **YOR323C — PRO2 — γ-谷氨酰磷酸还原酶（酿酒酵母实验删除株；名称由官方注释核对）**。p1630 的脯氨酸途径段、p1631 图5C：其与 CAR2 的双缺失呈更强生长缺陷，补脯氨酸改善；双缺失仍可生长，故不是完全致死或唯一旁路证明。
- **关键条件**：p1632 Experimental procedures 指明，该组最小培养基去除了硫酸铵，并以过量甲硫氨酸、亮氨酸和组氨酸提供氮，以避免脯氨酸摄取受氮分解代谢抑制。BY4741 背景，需氧，至少三条独立生长曲线。
- 含义：支持条件依赖的真实补偿潜力；不能将这组结果直接外推到含铵 SD-Leu 或 Yarrowia，更不能仅凭它确认当前模型的旁路通量大小。
- 访问：直接打开全文 PDF 并核查上述正文、图注和方法；图5C为遗传生长表型，不是原生 OAT 通量的直接测量。

### 4. 官方化学及辅因子定义（非靶蛋白原始实验）

NC-IUBMB 官方条目，均本次成功直接读取：

- [EC 2.6.1.13](https://iubmb.qmul.ac.uk/enzyme/EC2/6/1/13.html)：鸟氨酸的 δ-氨基转移，形成谷氨酸-5-半醛；明确为 PLP 蛋白。R610 使用 2-氧戊二酸作为受体符合该反应类型。
- [EC 2.6.1.19](https://iubmb.qmul.ac.uk/enzyme/EC2/6/1/19.html)：GABA 转氨，产物为琥珀酸半醛和谷氨酸。
- [EC 2.6.1.11](https://iubmb.qmul.ac.uk/enzyme/EC2/6/1/11.html)：N²-乙酰鸟氨酸转氨，产物保留 N-乙酰基。

共享 PLP 依赖的氨基转移酶家族并不意味着底物可互换。区分这三类需要具体底物、活性位点/序列对照及实验，不应只按通用结构域名称归类。官方等号也不证明任何具体培养条件下有显著双向生理通量。

### 5. 可追溯的参考蛋白面板（官方人工整理记录；非全部原始实验已重查）

| 酿酒酵母系统 ID／正式名称 | 蛋白功能与化学区别 | UniProt／长度 | 本次证据状态与用途 |
|---|---|---|---|
| YLR438W／CAR2 | 鸟氨酸氨基转移酶；自由鸟氨酸 → 谷氨酸-5-半醛 | [P07991](https://www.uniprot.org/uniprotkb/P07991/entry)，424 aa | reviewed 注释；上列1987原始功能互补支持，为 OAT 正参考 |
| YGR019W／UGA1 | 4-氨基丁酸氨基转移酶；GABA → 琥珀酸半醛 | [P17649](https://www.uniprot.org/uniprotkb/P17649/entry)，471 aa | reviewed 注释含实验催化和动力学引文；本次未独立打开其原始动力学论文，为近缘功能反参考 |
| YOL140W／ARG8 | 乙酰鸟氨酸氨基转移酶；N²-乙酰鸟氨酸 → N-乙酰谷氨酸-5-半醛 | [P18544](https://www.uniprot.org/uniprotkb/P18544/entry)，423 aa | reviewed 注释，线粒体蛋白；本次未独立重查其原始定位/酶学全文，为近缘功能反参考 |

这三个 accession 已交付主任务作有限 BLAST 对照。对照序列的实际获取、版本、全长一致性、比对参数与统计属于主任务的序列记录，本子任务不声称已运行比对。官方页面部分直接请求只返回 JavaScript 外壳；功能与身份核查使用官方搜索索引及主任务保存的官方记录。PRO2 的系统 ID 另由 [NCBI Gene 854501](https://www.ncbi.nlm.nih.gov/gene/854501) 核对，未把该辅助身份链接扩充为独立机制证据。

## 对 GPR 的建议及保留的不确定性

1. **当前单基因形式有生物化学合理性。** 保留目标为 OAT 候选，结合主任务的版本化序列、正/反参考比对再给整体证据等级。无依据新增 AND 伙伴；同源寡聚体也不产生额外基因条件。
2. **没有证据支持新增 OR。** 本次没有做全蛋白组同工酶排除或筛选，因此既不提出新 OR，也不声称该蛋白在 Yarrowia 中绝对唯一。
3. **胞质定位仍属待验证项。** 酿酒酵母先例与目标自动胞质注释可相互一致，但不等于原生定位实验。
4. **GPR 身份与条件性容量分开。** 单基因编码并不保证所有培养条件均表达或承载模型容量；若 essentiality 不匹配，应另核原生调控、底物供给及区室，而不是为匹配表型改写布尔关系。
5. 本范围未取得目标纯化酶的底物谱、动力学、PLP 依赖测定或复合体组成；也未核实该靶蛋白的体内逆向通量。结论为证据支持的候选保留建议，不是模型变更或目标功能的实验确认。

## 检索与审计界限

检索含目标新旧系统 ID、A0A1H6PJM3、Yarrowia + ornithine aminotransferase/CAR2、fungal OAT、酿酒酵母互补/定位/脯氨酸补偿，以及官方参考 accession。优先本物种，无直接命中后才使用近缘真菌；未扩展为系统综述。

本子任务已直接复核1987摘要、2007全文和官方化学条目；1993限于原始研究摘要；三个参考酶的全部原始实验链未重做。独立审计结论应由审计文件另行报告，不能把本文件的来源核查当独立审计。
