# R1931：逆向能力的独立酶学核查

核查日期：2026-09-15。范围：文献及官方酶命名资料；不运行 FBA、不改模型。此文件审查化学步骤及酶类方向，目标蛋白身份由本目录主报告另行核查。

## 结论

**没有找到支持把游离谷氨酸 + NADH → 谷氨酸-5-半醛作为高容量常规合成通路的证据；已有纯化酶实验反对这种建模方式。** 对确属 EC 1.2.1.88 的反应，生理建模采用氧化半醛、生成谷氨酸的单向能力，比允许 ±1000 更有证据支持。这是模型方向建议，不是关于所有化学条件的绝对不可逆定理，也不是本次已实施的修模。

本次未获得 Y. lipolytica 对应原生蛋白的反向动力学、细胞内底物/辅因子浓度或实测平衡常数；因此不能给出本菌株条件特异的 ΔG 或绝对逆向速率。未检出活性不等于逆向速率严格为零。

## 五个最相关来源与证据强度

### 1. 官方酶定义：EC 1.2.1.88

- 来源：[NC-IUBMB EC 1.2.1.88](https://iubmb.qmul.ac.uk/enzyme/EC1/2/1/88.html)。已直接打开完整短条目。
- 定位：Comments 第一段；条目末尾 EC history。
- 短引文："irreversible oxidation of glutamate-γ-semialdehyde to glutamate"。
- 含义：此类酶参与脯氨酸降解，催化半醛氧化为谷氨酸。旧编号 EC 1.5.1.12 在 2013 年转为 1.2.1.88。条目反应式中的等号不是体内双向通量的证据。
- 局限：官方酶类定义，不是 Yarrowia 原生蛋白的专门实验；质子系数受化学形态表示影响，应与模型化学单独核对。

### 2. 直接反向试验：大肠杆菌双功能酶

- Moxley et al., 2014, *J. Biol. Chem.* 289:3639–3651，DOI [10.1074/jbc.M113.523704](https://doi.org/10.1074/jbc.M113.523704)。原论文：[PMC3916563](https://pmc.ncbi.nlm.nih.gov/articles/PMC3916563/)。
- 定位：Results → “Steady-state and Single-turnover P5CDH Kinetics”，Fig. 4/5 邻接段。
- 短引文："no activity was detected in this direction"。
- 作者明确使用 NADH 和谷氨酸测试 P5CDH 逆向，没有检出活性，因而拟合动力学机制时将化学步设为不可逆。
- 准确边界：同段另提谷氨酸最高 50 mM 的产物抑制实验；**这不是明确给出的逆向试验底物浓度**，不能混用。已读段落未给逆向检测限、完整配方和温度，不补填。
- 获取等级：PMC 完整页直开遇到 reCAPTCHA；已从原论文网页的搜索索引读取上述完整结果段，并与原论文托管文本一致。不能声称本次成功下载全文。
- 易误读处：Discussion 后面的 reverse activity 指第一步脯氨酸脱氢酶把 P5C 还原成脯氨酸，不是 P5CDH 将谷氨酸变成半醛；不能当作 R1931 逆向的支持。

### 3. 限制性证据：大鼠肝酶的反向速率约束

- Small & Jones, 1990, *J. Biol. Chem.* 265:18668–18672；DOI [10.1016/S0021-9258(17)44804-6](https://doi.org/10.1016/S0021-9258(17)44804-6)。[PubMed 原始摘要](https://pubmed.ncbi.nlm.nih.gov/2211729/) 已直接打开。
- 原论文正文托管页：[ResearchGate，原论文公开全文条目](https://www.researchgate.net/publication/20943789_Pyrroline_5-carboxylate_dehydrogenase_of_the_mitochondrial_matrix_of_rat_liver_Purification_physical_and_kinetic_characteristics)。已读搜索索引中的 p.18670 正文，非平台自动概述；完整页直开失败。
- 定位：p.18670，反向试验结果及其后降 pH 试验。
- 短引文："no detectable P5C was formed"；"at least 15,000 times more slowly"。
- 作者没有观察到反向形成 P5C；据检测能力与常规正向试验，推定如果发生逆向，其速率至少比正向慢 15,000 倍。pH 从 8 降至 7 后仍未检出 P5C，作者认为存在显著动力学障碍。
- **这不是测得一个 1/15,000 的可靠非零逆向速率。** 它约束“绝对零活性”的表述，同时更不支持与正向同等级的高容量逆向。
- 未决：反向配方、绝对检出量及其完整推导的上一段本次未取得；不将 1/15,000 直接换算成模型通量界，也不声称已独立复核作者该倍数。

### 4. 第二个直接反向试验：芽孢杆菌纯化酶

- Isobe, Matsuzawa & Soda, 1987, *Agric. Biol. Chem.* 51:1947–1953；DOI [10.1080/00021369.1987.10868323](https://doi.org/10.1080/00021369.1987.10868323)。[J-STAGE 原论文 PDF](https://www.jstage.jst.go.jp/article/bbb1961/51/7/51_7_1947/_pdf) 已直接打开并读取正文。
- 定位：印刷 p.1950，PDF 第 4 页，反应产物鉴定段落末尾。
- 短引文："unidirectional conversion of L-P5C to L-glutamate"。
- B. sphaericus 酶的反向试验含谷氨酸、ATP、MgSO4 及 NAD(P)H，没有观察到由谷氨酸合成 P5C。作者据此报告单向性。
- 局限：同一段未列各试剂浓度或检测限。该酶偏好 NADP+，与模型 NAD+ 反应不可在定量上直接等同；证据是跨物种酶学支持。

### 5. 酵母同类酶的正向动力学

- *Structural Studies of Yeast Δ1-Pyrroline-5-carboxylate Dehydrogenase (ALDH4A1): Active Site Flexibility and Oligomeric State*, 2014；DOI [10.1021/bi500048b](https://doi.org/10.1021/bi500048b)。[PMC3954644](https://pmc.ncbi.nlm.nih.gov/articles/PMC3954644/)。
- 定位：Experimental Procedures → “Steady-State Kinetics”。
- 对 S. cerevisiae YHR037W / PUT2（已知 P5C 脱氢酶；该物种纯化酶实验证据）测量 NADH 生成：20°C、50 mM 磷酸钾 pH 7.5、25 mM NaCl、0.4 μM 酶、0.2 mM NAD+，L-P5C 1–300 μM。
- 含义：真菌近缘证据直接支持降解方向活性；所读方法并未实施逆向检测，**不能称为“酵母已实验证明无逆向”**。
- 获取等级：PMC 直开 reCAPTCHA；ACS 直开 403；上述原论文方法段从 PMC/ACS 搜索索引读取。

## 如何使用这些证据

- 酶学结论：支持按生理相关尺度近似为氧化方向单向；反对把反向当作可以替代谷氨酸磷酸化/还原步骤的常规高容量供给。
- 热力学结论：本次没有定量 ΔG 或 Keq，不能用“无 ATP”一个特点证明化学不可能。动力学缺乏可检出逆活性及明确的酶类注释，比单凭反应名称/数据库箭头可靠。
- 验收边界：只读审查支持候选方向约束，不能代替在保持模型、培养及阈值一致时重新检查 WT、目标单敲及其他供给路径。即使 R1931 禁止逆向，也不能未经定向验证声称目标必然变为 essential。
