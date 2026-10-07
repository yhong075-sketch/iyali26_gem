# R1889 四核心亚基部分 GPR：独立来源审计

日期：2026-10-05。限定任务为四组已有序列与 W29 locus 映射静态核查，以及两篇指定原始论文的必要依赖证据。未改模型、整理数据或运行求解、BLAST、结构预测、集群作业。

**判定：支持将下列四个已映射成员作为 R1889 的部分必要 AND 依赖。** 这是 Y. lipolytica 原始删除/组装证据向 W29 相同蛋白序列的映射；不是四基因已经在原生 W29 条件逐一敲除的声明，也不是完整复合体、完整 GPR 或最小催化系统的认证。

| W29 系统 ID—名称—功能 | 精确序列与 locus 审计 | 原始删除／组装证据 | 必须保留的限制 |
|---|---|---|---|
| **YALI1D07089g—NUAM—75 kDa核心铁硫电子传递亚基** | AOW03621.1＝CAB65519.1，728 aa全长逐字符相同；缓存UniProt Q9UUU3 的orfNames与AOW/CAB关联吻合。 | Waletko 2005 p5622 Methods、p5623 Results及p5624 Table I：nuam缺失空载体无可见完整复合体；WT回补为有活性对照。 | HAR人工受体残余不能证明完整NADH→醌反应存留；H129A点突变与删除分开。实验为NDH2i工程背景。 |
| **YALI1B26679g—NUBM—51 kDa FMN/NADH氧化核心亚基** | AOW01975.1＝CAB65520.1，488 aa全长逐字符相同；缓存UniProt Q9UUU2 的orfNames与AOW/CAB关联吻合。 | Maclean 2018 Fig.3及对应Results明确缺失株无complex-I活性；Fig.4C补充组装/活性染色。 | 少量约1MD的NUCM抗体信号可代表缺NUBM复合物，不能因此称完整酶具有NADH氧化活性。 |
| **YALI1F22993g—NUCM—49 kDa醌反应区核心亚基** | AOW07306.1＝CAB65521.1，466 aa全长逐字符相同；缓存UniProt Q9UUU1 的orfNames与AOW/CAB关联吻合。 | Maclean 2018 Fig.4C及其后Results明确删除株缺完整complex I。 | 本段直接端点是BN-PAGE、免疫印迹及活性染色的组装检查；不改写为本文对该W29基因逐一实测Q9泵通量绝对零。 |
| **YALI1F09003g—NUKM—PSST醌反应区核心亚基** | AOW06733.1＝CAB65525.1，210 aa全长逐字符相同；缓存UniProt Q9UUT7 的orfNames与AOW/CAB关联吻合。 | Maclean 2018 Fig.4C及其后Results明确删除株缺完整complex I。 | 本段直接端点是BN-PAGE、免疫印迹及活性染色的组装检查；不改写为本文对该W29基因逐一实测Q9泵通量绝对零。 |

证据等级：Y. lipolytica 参考亚基功能/删除依赖有原始实验支持；W29具体 locus 的功能承接由版本化蛋白精确序列和缓存来源映射支持。缓存为 unreviewed 数据库注释，不能把其自动功能文字独立计算为实验。

## 直接打开的来源与定位

1. [Waletko et al. 2005，DOI 10.1074/jbc.M411488200](https://www.jbc.org/article/S0021-9258(19)63065-6/pdf)：官方4页PDF文本，p5622 Methods给nuam缺失构建与WT/H129A回补；p5623 Results给空载体缺完整复合体；p5624 Table I具有WT、H129A及空载体对照。PDF截图失败，本审计不援引具体表格数值或声称完成图像复核。
2. [Maclean et al. 2018，DOI 10.1093/hmg/ddy247](https://academic.oup.com/hmg/article/27/21/3697/5048951)：PMC入口验证码，改读OUP官方全文成功。Results的complex-I activity、assembly小节与Fig.3/4图注直接支持上述结论；Methods“Yeast strains and growth”说明删除株来自既有菌株。图像链接失败，未声称独立重判条带。文章GB10/NDH2i背景及培养条件不能改写为本次W29生长必需性测试。

## 覆盖与边界

选定九项声明：四组精确身份对应、NUAM缺失依赖、NUBM缺失活性、NUCM缺失组装、NUKM缺失组装、这些必要依赖可作为不完整AND规则的有限模型推断。**total 9 | audited 9 | supported 9 | unresolved 0 | contradicted 0 | unchecked 0**，其中最后一项明确属于由证据支持的模型推断，不是额外实验。相同论文和同一记录衍生的检查不累加为独立实验份数。

未认证事项：完整42亚基AND、其他被省略亚基可有可无、NUHM的W29身份冲突、所有条件下的原生W29基因必需性、底物Q9专一性、原生膜泵计量及全网生长预测影响。没有把NDH2i工程旁路当作R1889原生同工酶。本文不评价当前其他反应旁路是否绕过这条部分GPR。

同轮 `evidence/REPORT.md` 曾提出包含NUHM的五成员待映射集合；本审计仅覆盖父任务明确要求的四成员规则，NUHM不得借用此审计结论进入执行GPR。
