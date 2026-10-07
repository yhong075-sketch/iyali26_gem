# R334 / R2081 生物学来源审查

核查日期：2026-09-11。执行者：biological_sources 子任务；原文审查已执行，科学实验和优化未执行。仅检索本对候选的文献与数据库，不改模型、GPR、实验标签或评价口径。独立审计由根任务另行记录；以下“已核查”是本子任务打开原文的状态，不等于独立审计通过。

问题：两种候选是否各自足以催化 ATP 依赖的谷氨酸＋氨→谷氨酰胺，且位于目标胞质区室；或二者是否共同必需。证据阈值：精确身份的分离酶学、遗传互补、单/双敲及定位优先；表达、共表达、通路图、同源注释只作间接支持。找到相反证据或原文不能打开时保留限定，不由基因名称、其他模型规则或搜索摘要升级结论。检索顺序：2017 原文和所引相关文献→2016 蛋白组→目标 ID、符号与独立催化/敲除/复合体关键词。停止条件：核心原文可查、OR/AND 证据缺口可定位；不扩展全基因组研究。

## 身份

| 模型系统 ID | 文献对应 ID / 名称 | 蛋白功能与证据状态 | 模型角色 |
|---|---|---|---|
| YALI1F00821g | YALI0F00506g / GLN1；2026 预印本称 GS1、YALI1_F00821g | 谷氨酰胺合成酶候选；名称有文献使用和数据库注释，本文献范围未找到该精确蛋白独立酶学验证（curated annotation / 间接支持；实验互换性未证实） | 根任务提供的 R334 / R2081 GPR 候选，实际规则由根任务静态核查 |
| YALI1D16151g | YALI0D13024g / GLN2；2026 预印本称 GS2、YALI1_D16151g | 谷氨酰胺合成酶家族候选；名称有文献使用，独立 GS 活性及底物特异性未在本文献范围核实（uncharacterized exact catalytic function / 间接支持） | 同上 |

YALI1↔YALI0 的确切序列对应由根任务核查；不能仅因符号相同当作序列身份相同。不同来源对第二个基因也有 GLN1 称法，因此报告应始终保留系统 ID。

## 核心原子声明与来源

| ID | 原子声明 | 来源 / 定位 | 类型与核查状态 | 限制 / 用途 |
|---|---|---|---|---|
| B01 | 2017 论文将 YALI0F00506g 称 gln1、YALI0D13024g 称 gln2 | [Pomraning et al. 2017](https://journals.asm.org/doi/10.1128/msphere.00038-17)，Results: “gzf2 is required for expression of nitrogen assimilation genes” | 原文正文已打开；命名支持 | 证明论文命名，不证明各自催化活性 |
| B02 | 2017 论文报告两基因在 ammonium 唯一氮源相对 peptone＋yeast extract 条件转录上调 | 同上，Results 同节 / Fig. 7 / RNA-seq methods | 原始转录组观察；正文已核查 | 未重算 RNA-seq；转录响应不证明独立催化或 AND |
| B03 | 2017 论文的基因删除实验针对 GATA 调节因子与另一个转录调节因子；Table 2 没有 GLN1 或 GLN2 单敲/双敲菌株 | 同上，Table 2；Materials and Methods: Identification and deletion of GATA transcription factors | 原始实验范围；正文与表已核查 | 文中下游基因“essential”措辞不能改写为对本对基因的敲除验证 |
| B04 | 2017 Fig. 11B 对氨同化通路明确作预测图示，相关生长缺陷伴随多个氮同化基因表达降低 | 同上，Fig. 11 caption / Results 同节 | 图注和作者机制推断；已核查文字 | 不能从多个下游共同降低推导本对蛋白共同必需；本任务未抽取图中定量 |
| B05 | 2017 氮同化段所引 ref. 54 是 S. cerevisiae 的 NAD-linked glutamate dehydrogenase 研究 | [Miller & Magasanik 1990](https://journals.asm.org/doi/10.1128/jb.172.9.4927-4935.1990)，Abstract；2017 References #54 | 原始论文摘要与引用归属已打开核查 | 物种和酶对象均不同，不能作为 Yarrowia 两个 GS 候选的酶学证据 |
| B06 | 2016 蛋白组研究报告 Gln1 与 Gln2 在氮限制下增加，并将其解释为 GS 氨同化路线 | [Pomraning et al. 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4766638/)，Discussion: Nitrogen assimilation；Methods: Yeast strains, cultivation and sample collection | 原始多组学；正文已打开核查 | W29，YNB 的 C/N=10 与150，25 g/L glucose；本任务未重核两蛋白独特肽段或补充表ID，不把该文作额外精确身份依据 |
| B07 | 2016 相关讨论引用59–61的支持来自微生物综述及 N. crassa 的氮代谢研究 | 同上，References #59–61 | 引文列表已核查；引文全文未逐篇打开 | 不能把跨物种通路常识升级为本对候选的直接验证 |
| B08 | 2026 预印本 Methods 明列本次 W29 GS1/GS2 系统 ID，且报告依次引入 GDH/GS/GOGAT 表达盒后的蛋白产量变化 | [He et al., preprint v1, 2026-04-28](https://www.preprints.org/manuscript/202604.1894)，Methods / Functional Validation 节；DOI 10.20944/preprints202604.1894.v1 | 预印本正文、PDF文本已打开；工程终点 | 该处仅统称 GS，未在所读文字区分 GS1/GS2 构建。附表及图像未核查，不判定哪一种单独有效；工程细胞保留内源背景，粗蛋白终点不等于独立酶学/区室证据 |

2017 同 DOI 的 PMC 页面在本子任务访问时要求验证码，随后改用 publisher ASM 完整正文；没有绕过验证码。2016 原文可打开。2017 主结论限于命名、转录响应及实际实验范围；未用未打开的历史对话补足证据。

## 新近线索：保留，暂不用于确认 GPR

2026 预印本网页报告已发表对应版本：[Nitrogen Availability Influences Biomass Composition in Yarrowia lipolytica Grown on Acetate](https://www.mdpi.com/2311-5637/12/7/315)，Fermentation 12(7),315，DOI 10.3390/fermentation12070315。web search 直接返回该 publisher 页的较长正文，日期为 2026-06-30；web open 则持续返回 HTTP429。为避免把搜索结果冒充直接打开正式原文，原始工具返回保存在 `sources/2026_mdpi_web_search.txt`，正式版内容当前标记 **search-derived / 独立待核查**。成功检索式：`"Nitrogen Availability Influences Biomass Composition" "GS1"`。

该搜索返回提及 GS 表达盒工程而非分别纯化的酶活实验，故不能概括为“从来没有本对候选的工程研究”。其 [公开审稿记录](https://www.mdpi.com/2311-5637/12/7/315/review_report) 的搜索结果另指出，原稿未提供表达量或酶活确认；本子任务没有打开完整回复，不判断最终是否补做。附录下载失败，不能以预印本和正式版相同为默认。

## 反证与排除范围

检索词覆盖两条 YALI0 ID、Q6C990、GLN2 与 Yarrowia，并与 knockout、deletion、heteromer、purification、enzymatic、overexpression 组合。本次没有找到已核实的本对蛋白异源复合体共同必需证据，也没有精确到两候选各自独立完成同一步胞质反应的分离酶学/互补证据。此为有边界的“未找到”，不是“不存在”。

- STRING 搜索出现两者高分关联，但未打开证据通道，且功能关联不等于物理复合体。排除其作为 AND 依据。
- CAoGDX 对 YALI0F00506g 的条目以高相似性描述其 GS 注释；不是本次新做 BLAST，也不算 Yarrowia 直接酶学证据。
- 搜索命中早期 Yarrowia 总 GS 活性研究（PMID 14526535），但只取得摘要级线索且未归属到两个精确基因，不用于判定 OR/AND。
- 搜索命中宿主表达异种 GS 的专利，以及列举两 ORF 的专利摘要；未核实其具体实验，不用于本对蛋白功能定论。
- 不把 iYLI647 的 AND 规则用作生物学依据；模型之间规则一致/不一致均不构成原生功能证据。

## 可支持的决策

目前可支持保留两者作为氮同化相关候选；不能据这些文献宣称 OR 已被验证，也不能据调节因子敲除、多基因共表达或通路图将其改为 AND。支持 AND 所需的是明确的共同必需/复合体依赖；支持 OR 所需的是每个候选在相关条件下对同一步骤的独立贡献，并核实底物、辅因子与区室。序列、活性位点及 AlphaFold 预测只能补充功能候选的强弱排序，不能替代这些验证。

独立审计覆盖：本子任务提交 8 条核心声明，独立审计状态尚未由本子任务评定。B08 的正式版、supplement、图像，以及额外搜索线索均保留待核查。
