# R1889：有直接功能证据的保守部分 GPR

核查日期：2026-10-05。目标为解脂耶氏酵母线粒体基质侧、质子泵型 complex I；本任务只读已有证据和原始文献，不改模型、运行求解或提交作业。W29 精确位点映射由同轮身份检查另行封存。

本轮采用 **NUAM AND NUBM AND NUCM AND NUKM** 作为有原始删除／功能缺失证据且已建立 W29 身份桥接的**部分必要依赖映射**。精确规则为：

`YALI1D07089g and YALI1B26679g and YALI1F22993g and YALI1F09003g`

它不是完整亚基清单、最小可重组酶，也不表示省略亚基可有可无。NUHM 的生物学依赖有直接证据，但本轮 W29 位点尚未闭合，故不纳入执行规则。未将七个核核心或全部 42 个结构成员一概硬 AND。

## 四个采用成员及一个未映射依赖的原始依据

| W29 位点／参考蛋白 | 名称与功能 | 原始功能证据 | 判定与边界 |
|---|---|---|---|
| YALI1D07089g；Q9UUU3；AOW03621.1；历史参考 CAB65519.1 | NUAM，75 kDa 核心铁硫电子传递亚基 | Waletko 2005，官方 JBC PDF p5622–5624：nuamΔ 空载体对照无可检测完整组装复合体；WT 基因回补恢复活性。H129A 另显示电子传递功能损害。 | 支持必要依赖。HAR 人工受体残余不代表 NADH→醌的完整反应仍在运行。原文 Table I 的 PDF 文本已读，截图失败，故本报告不引用其具体检测限数值。 |
| YALI1B26679g；Q9UUU2；AOW01975.1；CAB65520.1 | NUBM，51 kDa 含 FMN 的 NADH 氧化核心亚基 | Maclean 2018，Results／Fig.3：nubmΔ 无 complex-I 活性；Fig.4 分辨其组装状态。 | 支持必要依赖；不把人工 NADH 染色或 NDH2 活性误算为 R1889。 |
| YALI1F22993g；Q9UUU1；AOW07306.1；CAB65521.1 | NUCM，49 kDa 醌反应区核心亚基 | Maclean 2018，Fig.4C 及结果文字明确 nucmΔ 无完整 complex I；另将催化位点变体与删除分开。 | 支持完整复合体所需的依赖；不将单个变体的部分活性当作删除终点。 |
| W29 位点未决，未采用；Q9UUT9；CAB65523.1；原始 NUHM 核酸 AJ250338 | NUHM，24 kDa 核心亚基，参与电子输入模块 | Kerscher 2004，Results “Human NDUFV2 cDNAs do not complement…”：nuhmΔ 的 complex-I 特异 dNADH:DBQ 活性缺失，NUHM 回补恢复组装及活性。 | 功能依赖有直接支持；W29 位点未闭合前不得猜 ID 填入可执行规则。原文实验为 GB10/NDH2i 背景，不声称此次复现或 W29 定向删除实测。 |
| YALI1F09003g；Q9UUT7；AOW06733.1；CAB65525.1 | NUKM，PSST 核心亚基，参与醌反应区 | Maclean 2018，Fig.4C 及结果文字确认 nukmΔ 无完整 complex I；文中进一步将 nukmΔ 归为缺复合体对照。 | 支持必要依赖；原生必需性与 NDH2i 工程背景下的生长分开。 |

四个已列 W29 位点来自本轮缓存 `cached_core_identity.json`：审计者实际读取 UniProt genes/orfNames 的 AOW 与 CAB 连接；不是重做序列比对。记录为 TrEMBL/unreviewed，数据库自动功能语句不计成原生实验；功能依赖证据来自上述原始文章。NUHM 的尚缺身份桥接单独保留。本轮身份检查的 `../identity/four_core_comparison.json` 报告四个采用成员的版本化 AOW／CAB 全序列分别完全相同；本证据子任务读取了该检查结果，未重复比对。旧 NUHM 候选 ID 不能在身份未闭合时借用。

## 为什么不扩成全部亚基 AND

- 既有 `complex_i_gpr_evidence.csv` 缓存的全部成员都是 deferred；它主要列 PDB 结构存在证据，还有 ACPM 和 mt-ND3 身份冲突。旧 R570/R2062 的 AND 不是新的功能证据。
- NUIM（YALI1F01456g；Q9UUT8，TYKY 铁硫核心亚基）和 NUGM（本任务未认证 W29 ID，30 kDa 核心亚基）不纳入本次部分规则，**不表示其可缺失**。核七核心整体具有文献支持，但本轮没有逐一闭合这两个成员的原始删除终点。
- Kerscher 2002，DOI 10.1016/S0005-2728(02)00259-1，§8描述七核心基因的删除与回补；PubMed 将该文标为 **Review**。本轮仅把它作为原始实验线索，未把七项依赖全部标成独立打开原始数据后的结论。
- Dröse 2011 的 nb8mΔ 原始研究显示，多亚基缺失亚复合体仍可还原醌并保持减弱的质子泵功能。这支持“亚基存在／重要”不能直接替代二元完全失活证据；也不意味着其残余泵可继续使用未改变计量的 R1889。此次不添加残余反应或改质子系数。

## NDH2 与范围

Kerscher 2001 原始摘要明确 NDH2 是单亚基、非泵型、位于内膜外侧的替代酶；只有工程化重新靶向的 NDH2i 才能救援 complex-I 缺陷。因此不将 YALI1F32476g—NDH2（既有项目身份，外侧替代 NADH 脱氢酶）作为 R1889 的 OR 替代或 AND 亚基。它与本任务采用的四成员不同。

此四成员 AND 只纠正已证且完成身份桥接的必要依赖缺失。未列成员的敲除结果仍不能据此作完整的 complex-I 基因必需性预测；R2062 等既有替代连接也可能继续使全模型保持通量。本文没有证明其全网机制效果或生长预测改变。

## 一手来源与读取范围

1. Waletko et al. 2005. *Histidine 129 in the 75-kDa subunit…* DOI [10.1074/jbc.M411488200](https://doi.org/10.1074/jbc.M411488200)。[官方 PDF](https://www.jbc.org/article/S0021-9258(19)63065-6/pdf)已打开，重点 Methods、p5623 Results、Table I 文本；PDF截图工具两次失败，不声称图像已核。
2. Maclean, Kimonis & Balk 2018. *Pathogenic mutations in NUBPL affect complex I activity and cold tolerance…* DOI **[10.1093/hmg/ddy247](https://doi.org/10.1093/hmg/ddy247)**。[出版商主文](https://academic.oup.com/hmg/article/27/21/3697/5048951)已打开，核 Results／Fig.3、4 图注和正文。早期任务消息误写 ddy264，现已更正。
3. Kerscher et al. 2004. *Processing of the 24 kDa subunit mitochondrial import signal is not required…* DOI [10.1111/j.0014-2956.2004.04296.x](https://febs.onlinelibrary.wiley.com/doi/10.1111/j.0014-2956.2004.04296.x)。出版商全文已打开，核 Methods、删除／回补 Results；不把标题中的“信号切除不必要”误读成“NUHM 亚基本身不必要”。
4. Kerscher et al. 2001. *External alternative NADH:ubiquinone oxidoreductase redirected…* DOI [10.1242/jcs.114.21.3915](https://doi.org/10.1242/jcs.114.21.3915)。[PMID 11719558 作者摘要](https://pubmed.ncbi.nlm.nih.gov/11719558/)已打开；出版商全文入口失败，结论限摘要。
5. Dröse et al. 2011. *Functional Dissection of the Proton Pumping Modules…* DOI [10.1371/journal.pbio.1001128](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001128)。PLOS 原文已打开，核结果与讨论；PMC入口验证码后改用出版商。

本轮只提供证据整理，不替代同轮独立来源审核，也不声称完成七核或四十二亚基全面审计。
