# R1026 单基因候选赋值：独立审计

核验日期：2026-09-11。限定读取指定参考 SBML、已存 UniProt 条目及论文；无新搜索、BLAST、LP、结构作业或模型写入。未修改上一轮报告。

**建议：若 R1026 表达酶促碳酸氢盐形成，推荐单基因候选 GPR `YALI1F28274g`。结论层级为有数据库和模型对照支持的赋值建议，不能表述为已获原生酶学及胞质定位实验确认。**

| 基因身份 | 名称、蛋白功能及证据状态 |
|---|---|
| YALI1F28274g / UniProt A0A1H6PNW2 | 正式原生基因名未核实；Nce103-like β 类碳酸酐酶候选，催化 CO₂/HCO₃⁻ 可逆转换（自动功能注释、既有模型 GPR） |
| YALI0_F21406g / Q6C0V4 / XP_505708.1 | 正式原生基因名未核实；同序列碳酸酐酶候选（自动注释；相关重组表达实验未隔离目标特异活性） |
| YNL036W — NCE103 | S. cerevisiae 碳酸酐酶，催化 CO₂/HCO₃⁻ 可逆转换（该物种实验确认；跨物种参照） |

## 五条原子声明

| ID | 声明与审计结论 | 支持范围与局限 |
|---|---|---|
| C1 | **静态核实：所比较记录表示同一化学步骤。** | 实际逐项读取 `data/Yeast-GEM.xml:r_1667`、`data/chalmers_iyali/iYali.xml:R_y001667`、指定外部 `iYali5.xml:R_1667`、`data/iyli21.xml:R1026`，均为系数 1 的 CO₂ + H₂O ⇌ HCO₃⁻ + H⁺。前三者的物种 compartment `c` 明确命名 cytoplasm；iYli21 为项目约定的 `C_cy`，原文件 compartment name 实际写作 `unknownCompartment1`，此缺失注释未被回填。该对照按物种、计量及区室确定，不仅凭反应名称。 |
| C2 | **静态核实：模型间 GPR 不一致；遗漏原因未核实。** | Yeast-GEM 的 r_1667 直接关联 `YNL036W`，geneProduct label 为 NCE103；其余三条记录均未包含 geneProductAssociation。iYali 的 `NOTES: Manual curation` 和 iYli21 的 `PROTEIN_CLASS: 4.2.1.1` 不解释为何空 GPR。[iYli21 论文 §2.3](https://pmc.ncbi.nlm.nih.gov/articles/PMC9136261/) 说明以 iYali4 为模板、以双向 BLASTp 替换 CLIB122/W29 对应基因及增加反应，但正文没有这几个反应/目标的特定修改记录。因此不能由宏观流程断言作者漏掉了这个基因，也不能把空 GPR 确认为有意建模的自发过程。 |
| C3 | **数据库层支持：目标与该酶促化学相容。** | 已存 A0A1H6PNW2 v32 的 ORFNames 含 `YALI1_F28274g`，功能为可逆 CO₂ 水合、EC 4.2.1.1、β-CA；已存 Q6C0V4 v117 同样注释。独立提取两份序列逐字符相等，均 222 aa。功能证据是 ECO:0000256 的 RuleBase/ARBA；当前 A0A1H6PNW2 的删除原因是“不再属于参考蛋白组”，不是酶活被否定。此范围内未发现否定 CA 功能候选的直接证据；没有开展旁系同源穷尽搜索，不能据此声称它是唯一原生 CA。 |
| C4 | **未决：原生定位及目标特异活性仍不足以确认。** | 历史 A0A1H6PNW2 v32 的 GO:0005758 为 mitochondrial intermembrane space / IEA:EnsemblFungi，当前 Q6C0V4 的 GO:0005737 为 cytoplasm / IBA:GO_Central，存在定位注释冲突，双方均非直接实验定位。已核查的 [2026 论文 Table 1、§2.1、Fig. 1B、§4.6–4.9](https://pmc.ncbi.nlm.nih.gov/articles/PMC13207147/) 包含 XP_505708.1 对应 pCA 的 Po1h 来源克隆；pCA 与宿主对照总 CA 活性共享显著性组 CD，不能隔离确认目标酶活，亦不是目标无活性的证据。S. cerevisiae 原始酶学论文不解除这些物种与定位限制。 |
| C5 | **综合建议：R1026 的酶促版本可采用单基因候选 `YALI1F28274g`。** | 依据是正确化学、同序列功能注释、Yeast-GEM 的已知 CA GPR 对照及本项目 R2202 对同一反向计量的既有赋值。当前无依据加入额外 AND/OR 成员。该建议以“保留胞质酶促记录”的模型解释为前提；并不赋予真实非酶促过程 GPR，不确认自发容量，也不自动接受模型变更。上一轮核内供给见证仍有效，因此此赋值本身不等于恢复目标必需性。 |

覆盖：5 项均已审计；C1/C2 为模型静态事实（C2 的历史原因仍未知），C3 为数据库层支持，C4 为明确未决，C5 为综合候选建议。未发生科学模型变更。

## 后台来源身份

| 实际读取输入 | SHA256 |
|---|---|
| `data/Yeast-GEM.xml` | `9fd2c572cace73c2ea835205617313554d1bf89f4ef077f49defbdfa219a4ad7` |
| `data/chalmers_iyali/iYali.xml` | `c0b54165301bfba7edc7083efee1106309f049494f3a3f35fb436306d4fc86ae` |
| `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/reference/external_models/chalmers_yarrowia/eciYali5-GEM/model/iYali5-GEM/iYali5.xml`（只读指定文件，未进入该工作树） | `b73e99d88d2de954a0abe49bd326e8c1acf11c2c579c5a4c5629c78cbf5f0019` |
| `data/iyli21.xml` | `6974b7588f2a6c60ba2cde2f26e20d3aba1334d0d501572bf01cee47eda86631` |
| `sources/A0A1H6PNW2_v32.txt`（2026-01-28，序列 v1） | `914efdffa6b08952516a20e32f38c277a0ac6a9562c3b76c1a8da151a0718ea2` |
| `sources/Q6C0V4.json`（条目 v117，2026-09-02，序列 v1） | `071aa15c52fbc2d05cbf1764859e136980a334e4d13e83270bfa62dde430c913` |
| `sources/A0A1H6PNW2.json`（当前 inactive） | `9800e7b431513e02d387bc36079c857d6b1c34ed150ad45855250e66b1b605e3` |
| 两份条目的规范氨基酸字符串 | `61eeba8a4906133208f71371cc8334cc79404bb12bfe8253ed1c09a07829b454` |
| `sources/PMC9136261.html` | `75c8b35f4897f05cf4de5cc36b171eb95570727a6a2dd7e33462024313021d21` |

其余论文的获取、实际查看图件及指纹见本目录上一轮 `SOURCE_AUDIT.md`；本轮仅按范围复核，不把重复阅读标为新实验。
