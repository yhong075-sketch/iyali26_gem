# 独立来源审核

核查：2026-10-02 22:10–22:12 UTC。审核者独立打开下列原始/数据库来源；未据代理一致性替代来源审核。范围为主代理提交的 4 项声明，未审核任意新增模型/GPR或最终未提供的其他声明。

| ID | 待审声明 | 独立检查 | 判定 | 限制/处理 |
|---|---|---|---|---|
| C1 | 原命名物种主要泛醌为 Q-9，属于 Yarrowia clade | [PMID18218960](https://pubmed.ncbi.nlm.nih.gov/18218960/) 原始论文摘要，明确两株来源、系统发育及 Q-9 结果 | supported | 属实验表型/系统发育报告；本轮未复现测定，不能推导末端酶化学 |
| C2 | NCBI taxon444778 当前名称为 Y. phangngaensis，并关联旧拼写 | [NCBI444778](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=444778)，current name、basionym、原文参考 | supported | 仅明确 NCBI 当前登记；正式命名作者/年份由原始命名文献另论。收藏号冲突未消除 |
| C3 | GCA_900519005.1 与 GCA_030581735.1 为可用参照组装；存在本身不能指定末端酶氧化态或供体 | [Brinkrolf2021](https://link.springer.com/article/10.1186/s12864-021-07597-z) 方法明确 CBS10407/GCA_900519005.1；[PRJNA736342](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA736342) 明确 NRRL Y-63743/GCA_030581735.1 | supported | 第一半为来源直接支持；第二半为证据范围推理。未下载核定完整组装内容，未将两组装当作相同序列或正式项目基线 |
| C4 | 本次有限检索尚未确立该物种自由 DMQ9H2 羟化供体/净反应或验证末端酶位点/GPR | 独立五批检索、UniProt2026_03 taxon查询及访问失败记录，见 NOTES.md | partially_supported / unresolved | 支持“本轮尚未找到”的过程陈述；生物学答案仍未决。不能扩大为不存在此类酶/证据，不能声称 KEGG/BioCyc 全库已排除 |

覆盖计数：total claims 4 | audited 4 | supported 3 | unresolved 1 | contradicted 0 | unchecked 0。

C4 的保守表述可以使用；未决项限制机制结论和模型接受，不妨碍报告该物种为 Q-9 生物学参照。没有证明 R695/R385 在目标物种按某指定自由底物、供体和 GPR 实现，也没有授权模型闭合或激活。

另于 2026-10-02 22:12 UTC 读取本目录上级 REPORT.md：其两条方程明确标作此前比较候选、供体未指定，且未当作本物种实测；未见这些结论的措辞越界。此项为表述范围检查，不增加上述 4 项来源审核分母；本轮没有重新审核前轮反应分类/酿酒酵母酶学来源或2026命名全文。24 条数据库记录来自 **UniProt**，不能写成 NCBI 蛋白记录数。
