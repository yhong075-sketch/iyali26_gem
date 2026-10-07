# 游离 Gly-Pro 水解：原始酶学证据审查

核查日期：2026-09-24。对象为 **YALI1E16433g**（本轮未核实到已确立的原生名称；M24B/prolidase-like 候选；目标菌株序列身份由主任务独立核验）。本子任务只检索、读取和保存生化证据及比较蛋白记录，没有改变模型、GPR、培养条件，没有运行优化或蛋白预测。

**现有实验参照使“可水解游离 Gly-Pro”成为有根据的候选功能，但尚不能确认目标蛋白具有游离 Gly-L-Pro 水解能力。严格的 prolidase 分类、底物能力和细胞区室是三个独立问题。现有参照也不要求结果支持液泡 R2039。**

| 实验参照 | 与游离 Gly-Pro 直接相关的结果 | 解释边界 |
|---|---|---|
| *A. nidulans* pepP；原文对应 AJ296646.1 / CAC39600.1，原文预测翻译 496 aa | pepP 过表达菌株纯化制剂可水解 Gly-Pro；速率为 Phe-Pro 的 15%。标准测定 pH 7、37°C、0.05 mM Mn、3.5 mM 二肽。最适 pH 7；EDTA 失活后 Mn 恢复。Pro-Ala、Ala-Ala 在所测条件下未水解，Ala-Pro-Gly 的 Ala–Pro 键未被水解 | 原文写 Gly-Pro，试剂写 Sigma 二肽；所读章节未给 L 前缀或货号，不能把立体化学登记写得比原文更明确。没有 Gly-Pro 的 Km/Vmax。原文 Table 3 单位有可见冲突，本报告不搬用该表动力学常数 |
| 人 PEPD / P12955 | 重组 WT 对 Gly-Pro 有实测活性；231delY、E412K、G448R 的活性和动力学明显改变 | 这是参照蛋白的直接证据，不能仅由结构相似或活性位点保守外推目标精确底物。原文也以 Gly-Pro 命名 |
| *L. lactis* NRRL B-1821 PepQ / ABW84230.1（EU216565.1） | Gly-Pro 相对活性在 Zn、Mn 两组均 <0.1%，参照是 Leu-Pro/Zn = 100%；pH 6.5、2 mM 底物、1 mM 金属、50°C。其他 X-Pro 有活性 | 提供重要反例：被实验称为 prolidase 的蛋白也未必有效水解 Gly-Pro。不能把这个条件性阴性转移成 Yarrowia 阴性 |

来源：[Jalving 原始研究作者存档章节，印刷页84–90](https://edepot.wur.nl/121628)、[Besio 2013，Table 1 / Figure 2](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0058792)、[Yang 2008，Table 2](https://febs.onlinelibrary.wiley.com/doi/10.1111/j.1742-4658.2007.06197.x)。上述同源/参照蛋白是否为目标的最近邻及其定量序列关系，由主任务报告；本子任务没有以物种近缘代替蛋白近缘。

## 需要保留的身份与注释冲突

*A. nidulans* 原文所报 CAC39600.1 为 496 aa；当前 reviewed UniProt **Q96WX8**（sequence version 2，entry version 117）为 465 aa，逐位恰好等于 CAC39600.1 的第 32–496 位。CAC 多出的前31 aa 是 `MQAALHRTEIKAPHRPTRALSNLFTARNRIA`。当前记录提示错误起始／N端延长，并按同源证据标为 probable Xaa-Pro aminopeptidase（EC 3.4.11.9）。这不能抹去原文的 Gly-Pro 活性，但也不能把当前465-aa条目无条件称为原实验完全同序列的蛋白。原文的496-aa推定 ORF 及过表达实验，并未直接测出所纯化成熟蛋白的N端。记录与实测准备之间的差别必须保留。

定位同样需要分层：原文先测到细胞提取液活性、培养滤液无活性；“胞质”来自作者对分选/TM信号缺失的推断。它不等于直接亚细胞定位，更不能转移成 YALI1E16433g 的胞质或液泡实验结论。

## 功能分类与下一步判据

[IUBMB EC 3.4.11.9](https://iubmb.qmul.ac.uk/enzyme/EC3/4/11/9.html) 明确允许 aminopeptidase P 从二肽切下与 Pro 相连的 N 端残基。因此，即便直接证实 Gly-L-Pro 水解，也不能单凭这一底物唯一判为严格 prolidase。DPP-IV 从更长底物释放二肽；其 Gly-Pro-pNA 等底物测定切断的是另一条键，不能代替游离 Gly-Pro 水解。

最有信息量的功能验证是身份明确的目标蛋白与游离 Gly-L-Pro，确认底物损失和 Gly、L-Pro 产物；同组加入其他 X-Pro、Pro-X 和含 N端 X-Pro 的三肽/长肽，区分底物范围。金属与 pH 应覆盖合理条件，并有空白和已知活性阳性参照，不能将单一条件无信号写成无功能。原生定位应独立验证；酶活最适 pH 不等于细胞区室。无论酶活最终阳性或阴性，均不能单独授权将该蛋白赋给液泡 R2039。

M-CSA 人 PEPD 参照位点（UniProt编号）：D276、D287、H370、E412、E452 为金属配位残基；H255、H377、R398 与底物结合/机制相关。需通过实际比对映射目标位点，不能凭相同数字映射；保守性只支持候选。机制中的具体质子转移模型也应保留其推断等级。

有限声明为 BIO01–BIO09，见 `claims.json` 和上一级 `biochemistry_evidence.tsv`；原始来源、记录版本、检索限制及文件身份见 `source_manifest.json`、`retrieval_log.json`。公开原文按来源许可仅保存在证据目录；正文为有限摘述与分析。
