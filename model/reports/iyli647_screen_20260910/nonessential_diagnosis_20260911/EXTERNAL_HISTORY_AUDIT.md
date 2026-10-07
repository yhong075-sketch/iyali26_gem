# iYLI647 作者版本历史独立审计：11 个外部 TP／本项目 FN 目标

核验时间：2026-09-11T17:23:01.554890+00:00。本次读取固定源 JSON、作者笔记本代码及官方 GitHub 历史／完整模型补丁；没有求解或重跑敲除，没有修改模型、培养基、GPR、实验标签或既有结果。

**11/11 个目标及全部 23 条直接关联反应在仓库起点 `iYLI647_corr.json` 中已经存在；至 `corr_2`、`corr_3`，其关联集合、GPR 布尔语法树、边界与计量均无实质变化。** 这些差异不能称为作者这轮针对 11 项敲除结果完成的修复。相关网络／规则在本仓库起点已经存在；但 `corr` 未被独立确认是完全未修改的上游原始 iYLI647，历史版本的 KO 表型也未在本审计中计算。

## 目标、方法与 23 条直接反应

目标取自 [共同正例比较表](/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/iyli647_screen_20260910/common_positive_comparison.tsv)，筛选 `scenario=mapped28_po1f`、`external_at_10pct=TP`、`iyali26_at_10pct=FN`，恰有 11 行。表 SHA-256：`1616198c1f33644ef9b7f41ef795b37260a380f6401d59667a1d4fb8e52f71fc`。已保存的 TP／FN 仅用于选取目标，未在此重新验证分类；跨版本 ID 对应继承该表，不按反应编号推断。

从每份 GPR 的 Name 节点提取全部关联反应，采用 Python AST 比较 Expression／BoolOp／And／Or／Name／Load 结构，以排除括号与空白差异；另核对存在性、上下界、全部代谢物系数。三份模型的大小、完整 SHA-256 与 Git blob SHA 匹配下载清单。

以下基因的已核实名称均为“未核实”；功能是模型赋值／计量描述，证据等级为 **model/GPR assignment only**，不代表原生蛋白功能被实验确认。每行的全部关联结构在三版均不变。

| 外部系统 ID／本项目 ID | 模型赋予的功能 | 全部直接关联反应 |
|---|---|---|
| YALI0A02310g／YALI1A02775g | 葡萄糖-1-磷酸尿苷酰转移 | `GALU` |
| YALI0B15598g／YALI1B20462g | 6-磷酸葡萄糖酸脱氢 | `GND` |
| YALI0C05951g／YALI1C07638g | 棕榈酰/硬脂酰 CoA 去饱和 | `DESAT16`, `DESAT18` |
| YALI0C06490g／YALI1C08702g | 甘露糖-1-磷酸鸟苷酰转移 | `MAN1PT` |
| YALI0C11407g／YALI1C15991g | 乙酰 CoA 羧化；模型还将目标写入 11 条脂肪酸合成 AND 规则 | `ACCOACr`, `FAS100COA`, `FAS120`, `FAS120COA`, `FAS140`, `FAS140COA`, `FAS160`, `FAS160COA`, `FAS180`, `FAS180COA`, `FAS80COA_L`, `FAS80_L` |
| YALI0C23364g／YALI1C32184g | 内质网蛋白甘露糖基转移 | `DOLPMMer` |
| YALI0D03069g／YALI1D03865g | 磷酸核糖甘氨酰胺甲酰转移 | `GARFTi` |
| YALI0E18964g／YALI1E22736g | 甘油脂酰基转移 | `AGAT_SC` |
| YALI0E21021g／YALI1E25018g | 1,3-β-葡聚糖合成 | `13GS` |
| YALI0F00506g／YALI1F00821g | 谷氨酰胺合成 | `GLNS` |
| YALI0F02497g／YALI1F03803g | 计量为 homocitrate 脱水；反应标题却写 methylcitrate，不能按标题解释 | `MCITDm` |

14 条反应在 `corr → corr_2` 仅去掉 GPR 外括号：`MAN1PT`、上表 11 条 FAS 反应、`DOLPMMer`、`GLNS`；`corr_2 → corr_3` 连这些字符串也未变。AND 依赖也早已存在：`MAN1PT` 为 2 基因 AND，11 条 FAS 为 4 基因 AND，`DOLPMMer` 为 3 基因 AND，`13GS`／`GLNS` 各为 2 基因 AND。精确规则和逐版计量见 [逐目标静态记录](/private/tmp/worland-fixed-audit-20260910/eleven_gene_history.json)。

## 已定位机制的结构历史

另核对以下 12 条上下文反应；GPR AST、边界、计量均相同，其中 9 条的整个反应 JSON 对象也相同。这是结构历史核查，不替代主分析的稳态机制证明。

| 反应 | 三版静态结论／固定源码 |
|---|---|
| `FAS161ACPm` | 仅 GPR 外括号变化；AST、边界和计量相同；[官方源码](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_3.json#L15650) |
| `G6PDH2er` | 整个反应对象相同；[官方源码](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_3.json#L16898) |
| `6PGLter` | 整个反应对象相同；[官方源码](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_3.json#L8228) |
| `PGL` | 仅 GPR 外括号变化；AST、边界和计量相同；[官方源码](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_3.json#L17020) |
| `GND` | 整个反应对象相同；[官方源码](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_3.json#L16957) |
| `HCITSm` | 整个反应对象相同；[官方源码](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_3.json#L18916) |
| `MCITDm` | 整个反应对象相同；[官方源码](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_3.json#L20716) |
| `HACNHm` | 整个反应对象相同；[官方源码](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_3.json#L20674) |
| `HICITDm` | 整个反应对象相同；[官方源码](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_3.json#L20687) |
| `MCOATAm` | 整个反应对象相同；[官方源码](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_3.json#L19340) |
| `PC` | 整个反应对象相同；[官方源码](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_3.json#L8468) |
| `ATPCitL` | 仅 GPR 外括号变化；AST、边界和计量相同；[官方源码](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_3.json#L10436) |

`GND`、`G6PDH2er`、`6PGLter`、`PGL` 的依赖关系不是本轮后续删反应造成的变化；全模型修订链没有任何反应删除。`MCITDm` 的实际代谢物为 `hcit[m]`、`b124tc[m]`，并连接到 `hicit[m]`；不能按标题将其当成甲基柠檬酸循环修复。`FAS161ACPm` 已在 `corr` 存在且边界开放，三个版本中 `malcoa[m]` 都只连接 `MCOATAm`，没有新增其供给来源；因此此处不是作者在 corr2/3 删除线粒体 FAS 旁路。`PC`、`ATPCitL` 在三个 GSM 中仍是 `[0,1000]`。

完整逐版对象及源码行号见 [机制上下文历史](/private/tmp/worland-fixed-audit-20260910/eleven_structural_context_history.json)。

## 作者实际修改与未修改的内容

固定提交内部的 **`corr → corr_2`** 新增 3 条反应：胞质 `MALS`、合并 β-carotene 生产 `caro_prod`、`EX_caro_e`；删除 0 条，既有反应计量变化 0 条。`MALS` 空 GPR，计量为 `accoa[c] + glx[c] + h2o[c] → coa[c] + h[c] + mal_L[c]`。`EX_caro_e` 实际消耗 `caro[c]`，不能按名字推断其区室。`CRNCARtm`、`CSNAT` 从 `[0,1000]` 改为 `[-1000,1000]`。没有既有 `(0,0)` 反应被打开。

224 条 GPR 字符串变化中，218 条 AST 相同；另外 6 条是 `DPGM` 重复项简化（布尔逻辑仍等价），`CYOOm`／`CYOR_u6m` 各去掉一个截断 ID，`HCAt`／`dca_t`／`ACACCT` 的拼写修正。它们都不是这 11 个目标的直接关联反应。作者代码将去掉条目标作错注释，章节依据为“不在转录组数据集中”；本审计确认作者做了什么，不认定未出现于转录组本身证明基因不存在。

代码定位：[Supp A](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/notebooks/Supp_A_gsm_modifications.ipynb)，零基 cell 3 载入 corr；5 添加 MALS；7／8 改 carnitine 可逆性；10 添加 carotene；12／14／15／17／18／19 修改注释；21 保存 corr_2。

**`corr_2 → corr_3`** 仅新增 `biomass_glucose`、`biomass_oil` 两条反应，各含 47 个非零计量项。既有反应的计量、边界和 GPR 不变；`biomass_C` 本身及其唯一目标系数 1 在三份模型中均保留。不能仅因加载 corr_3 就认为 glucose／oil biomass 已成为目标。按大分子类别碳摩尔量缩放的代码与两条交付反应系数吻合（此前静态重算最大绝对差约 1.11×10⁻¹⁶）。代码定位：[Supp B](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/notebooks/Supp_B_gsm_biomass_reactions.ipynb)，零基 cell 3 载入 corr_2；17 添加两条 biomass；24 保存 corr_3。模型规模依次为 1348／1351／1353 反应、1121／1122／1122 代谢物、646／648／648 个基因记录。

**Git 日期历史与上述文件后缀不同。** [初始提交 2df44b090ca0b7a90cce27f75ec610bbe4c18f67](https://github.com/UH-MBBE/yarrowia-13C-gsm/commit/2df44b090ca0b7a90cce27f75ec610bbe4c18f67)（2024-01-31）已经同时交付三份模型；`corr` 此后没有修改。已将三个文件的历史查至该无父提交。后续模型补丁只有两次，已直接读取完整补丁：

- [121b3aaa30adce4da07e36efe5d5f54fe9f250ab](https://github.com/UH-MBBE/yarrowia-13C-gsm/commit/121b3aaa30adce4da07e36efe5d5f54fe9f250ab)（2024-02-05）：corr_2／3 只改两条 carnitine 边界。
- [06505c26b078afe9a7edc2832a10f6e85199ca7b](https://github.com/UH-MBBE/yarrowia-13C-gsm/commit/06505c26b078afe9a7edc2832a10f6e85199ca7b)（2024-02-05）：corr_2／3 只加 carotene 物种及两条反应。

两次后续补丁均无目标 GPR 修改或旁路删除。MALS、biomass、注释修订已在初始上传的相应文件中。公开历史始于“initial commit in new repo”，不能据此补齐此前的上游建模过程。[作者历史来源记录](/private/tmp/worland-fixed-audit-20260910/eleven_author_history_source_audit.json) 保存逐文件历史和模型补丁。

论文 Worland 等（2024），[DOI 10.1016/j.ymben.2024.06.010](https://doi.org/10.1016/j.ymben.2024.06.010)，Methods §2.4–2.5 描述 biomass／annotation／GSM flux bounds 修订；§3.3 描述油酸 13C-MFA 中移除 pyruvate carboxylase／ATP citrate lyase 以缩小置信区间。相关定位由 [作者公开论文页面](https://garrettroell.com/papers/yarrowia-tag-mfa-gsm) 的网页索引文本核对，本轮未声称重新读取完整 PDF。MFA 条件约束与交付 JSON 的默认边界是不同输入，不能把油酸 MFA 的删反应说明移植成这些 GSM 已关闭相应反应。当前对照也不能归因于本次未施加的作者 13C 约束。

## 审计覆盖与限制

按“11 个目标的全部关联结构不变”和“12 个上下文反应结构不变”计，**total 23 | audited 23 | supported 23 | unresolved 0 | contradicted 0 | unchecked 0**。上下文与目标反应有重叠，23 个审计命题不是 23 个不同反应的总数。另已完整读取 3 个模型逐文件历史、2 次后续模型补丁；跨文件新增／删除／边界／计量作了全量比较。此覆盖率只适用于列明结构／历史命题，不覆盖所有生物学机制或实验结论。

未证明：历史 corr 的 KO 生长率相同；所有全局改动均无间接影响；作者有意用这些结构修复 11 项必需性；模型规则等同实验确认的原生功能；相同反应名等同相同化学。未执行反事实模型计算，“直接规则未改”与“所有修改对表型毫无影响”必须分开。本次研究模型与原始输入均未更改。

## 汇总文档的五条完整 GPR 复核

补充核对 [IYLI647_COMPARISON.md](/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/iyli647_screen_20260910/nonessential_diagnosis_20260911/IYLI647_COMPARISON.md) 的“外部的完整 GPR 规则”五行：`MAN1PT`、`DOLPMMer`、`13GS`、`GLNS`、`AGAT_SC`。**5/5 行与固定外部 JSON 的完整 GPR AST 相符，5/5 在三个历史文件中保持相同，5/5 在目标为 false、其余成员为 true 时规则为 false。** 这确认外部模型中目标 KO 关闭该步的布尔解释，无需求解。

汇总文本明确说明成员名称未核实、功能只是模型候选、AND 不等于经验证的必需复合物；其表述没有将 AND 当成原生实验验证。本项不单独证明整条合成路线唯一、跨模型伙伴身份相同或本地反应关系，这些属于另外的结构审计。5 个规则命题单独计数：total 5 | audited 5 | supported 5 | unresolved 0 | contradicted 0 | unchecked 0。与上述 23 个结构历史命题合计为 28/28；这是本历史审计内部计数，主汇总另列的 6/6 结构机制由其独立机制审计负责，不能与本表的 12 条上下文反应混为同一统计单位。

核对时汇总文档的完整 SHA-256 及五条规则细节见 [五条 GPR 复核记录](/private/tmp/worland-fixed-audit-20260910/five_comparison_gpr_audit.json)。

## 后台输入身份

仓库固定提交：`3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6`。2026-09-10 下载，2026-09-11 重新核对。完整哈希供复核，不建议在用户解释正文重复。

| 文件 | 完整 SHA-256 | Git blob SHA-1 | 字节 |
|---|---|---|---|
| iYLI647_corr.json | `96d1ee6bc5bbb38a447490b8c343aebd41846448218f96315ef80e6b00360bfe` | `b0f4842a96776bb709d9baac09bfd60030776763` | 507491 |
| iYLI647_corr_2.json | `0bc227319eab7d18d96e035dc1298d3d6e60dfa7c89bd8acbba68a7398ce646c` | `7279465cf413ee094bc1d7fb66924df603c50e8f` | 507462 |
| iYLI647_corr_3.json | `329be540c099409c2c7b76ee581a23f86eefdbbffe97e9b03b39b5e9c014b5d2` | `ba298186c3f36cb309a599c644df41af842d7751` | 510627 |
| Supp_A_gsm_modifications.ipynb | `b0ba96869c74ad7abca53c5e1f00e9a644ca043ca8b7d32b901fcdc8a973ae34` | `d856157c156e7be6ba50b27f6edfccadc13b3aca` | 28326 |
| Supp_B_gsm_biomass_reactions.ipynb | `729e3b81cdcd749121b7dc3913f74cadefd9efd6830e1381f29fbec6859d871c` | `7eafbcbf320e5b1eb6818e38b502e2a725ae57de` | 70221 |

官方固定入口：[corr](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr.json)、[corr_2](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_2.json)、[corr_3](https://github.com/UH-MBBE/yarrowia-13C-gsm/blob/3748dc381aeb03b297b8f9eb7b0e38253fc7b8a6/genome_scale_models/iYLI647_corr_3.json)。
