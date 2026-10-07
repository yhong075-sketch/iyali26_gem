# COQ3 身份映射与 R385 局部提案：独立来源审计

核查日期：2026-10-02。范围为原始数据库记录、当前指定模型编码和 IUBMB 定义；未运行 BLAST、结构预测、模型构建或求解，未改模型。

**审计结论：当前明确版本 `AOW01767.1` 与 `XP_500950.3` 可确定为同一 W29 基因、同一完整蛋白序列的 GenBank—RefSeq 对应。** 模型 `YALI1B20835g`（COQ3 家族 O-甲基转移酶候选，原生功能尚未实验验证）与官方 `YALI1_B20835g` 的关联由模型交叉引用和官方 locus 同时支持。不是将不同株系蛋白仅凭名称拼接。

| 声明 | 原子核查内容 | 独立打开的来源／定位与限制 | 判定 |
|---|---|---|---|
| I1 | 模型YALI1B20835g与官方locus YALI1_B20835g指向同一W29候选基因；是项目ID格式归一化，不能仅凭字符串近似断言。 | 当前XML geneProduct的XP_500950.3与GeneID2907025交叉引用；AOW、XP3 CDS/locus_tag及Gene XML。 | supported |
| I2 | AOW01767.1与XP_500950.3是同一W29基因对应的GenBank与RefSeq蛋白记录，而非同一个数据库记录。 | 两者FEATURES source均CLIB89(W29)、CDS均YALI1_B20835g；XP3 COMMENT明确与AOW01767相同；二者accession不同。 | supported |
| I3 | AOW01767.1与XP_500950.3的完整蛋白序列均367 aa且逐字符相同。 | 独立解析两份原始GenPept ORIGIN并计算SHA256，均6e968848ccfdf9ebad65969d75001a0cdd0ffb0665403c7232e05ccc2209c051。 | supported |
| I4 | 两记录的编码来源链可建立：AOW来自CP017554.1:2083566..2084669；XP3来自XM_500950.3:1..1104，GeneID2907025定位NC_090771.1的同一数值区间。 | AOW CDS/coded_by；XP3 CDS/coded_by；XM3 CDS与COMMENT；Gene XML genomic/products，0基坐标2083565..2084668换为1基。基因组全序列字节等同性不在本核查声明内。 | supported |
| I5 | XP_500950版本历史伴随来源株系、locus和蛋白长度变化，不能无版本地混用。 | 历史原始GenPept：.1 CLIB122/YALI0_B15884g/367aa；.2 DSM3286/YALI2_C00448g/295aa；.3 CLIB89(W29)/YALI1_B20835g/367aa。.2于2024-06-26替换.1；.3于2024-09-10替换.2。独立序列检查.1等于.3，.2等于.3残基73..367。不能由此判定真实生物学N端缺失原因。 | supported |
| I6 | 身份映射解决不等于W29原生COQ3催化功能已获实验验证。 | AOW产品hypothetical protein，XP3产品uncharacterized protein；XP3为PROVISIONAL REFSEQ，conceptual translation；similar to COQ3属于推断注释。模型name COQ3及GPR不是酶活实验。 | supported |
| C1 | EC2.1.1.64的当前标准反应为还原态去甲基泛醇到泛醇的SAM依赖甲基转移。 | 独立打开IUBMB条目，accepted name与Reaction行：3-demethylubiquinol-n与ubiquinol-n；副底物SAM、副产物SAH。链长n通用，不是W29专属酶学。 | supported |
| C2 | 指定候选模型R695生成m611醌态DMeQ9；R385当前消耗m611并生成m468醌态Q9。 | 独立XML读取R_R695、R_R385及species。m611 C53H80O4；m468 C54H82O4；同为C_mi。R385另有m60 SAM→m62 SAH，各系数1。 | supported |
| C3 | 指定模型中没有可直接用于该提案的DMeQ9H2醌醇代谢物；其所需中性分子式C53H82O4未出现。 | 独立扫描1879 species：名称含demethyl/ubiquinol项与全模型formula检查，C53H82O4零匹配。现有coq_dmq9h2是C53H82O3（少一个氧），不能替代。该声明限于当前编码，不是生物学不存在。 | supported |
| C4 | Q9H2已经编码为m471，且有既存氧化还原连接；R385改产物可以复用该条目。 | 独立XML扫描2314 reactions：m471/C54H84O4/C_mi。正向生产R262/R570/R740/R1889/R1977/R2062；R573计量生产但边界0,0；R305正向消耗，R740允许逆向消耗。编码/边界不代表已验证可行通量或唯一生理路径。 | supported |

审计覆盖：**total 10 | audited 10 | supported 10 | unresolved 0 | contradicted 0 | unchecked 0**。此数仅描述上述十项声明，不表示 W29 催化机制已全部解决。

R385 按 IUBMB 定义改为 DMeQ9H₂＋SAM→Q9H₂＋SAH，会暴露由现有 R695 的醌态产物到新醌醇底物的供给缺口。这是物种分离和模型计量产生的局部一致性问题，不能直接推出独立 W29 还原酶、NADH 供体或新反应 GPR。该步骤仍未获本次证据认证，不应自动补反应。

来源入口：[AOW01767.1](https://www.ncbi.nlm.nih.gov/protein/AOW01767.1)、[XP_500950.3](https://www.ncbi.nlm.nih.gov/protein/XP_500950.3)、[GeneID 2907025](https://www.ncbi.nlm.nih.gov/gene/2907025)、[IUBMB EC2.1.1.64](https://iubmb.qmul.ac.uk/enzyme/EC2/1/1/64.html)。NCBI 原始记录由身份代理下载，本审计独立读原件并重算序列哈希；不是仅读其结论。IUBMB 由本审计通过网页工具独立打开。

审计还保留两项边界：XM_500950.3 标记为两端不完整的 mRNA 模型，CDS 为 1..1104；本次确认其编码的蛋白记录关系，未声称验证完整转录本边界。XP 历史株系/长度变更不能用来推定生物学截短原因。

## 终稿核对

核对时间（UTC）：2026-10-02T23:32:14.299719+00:00。已逐段核对主报告身份表、RefSeq版本史、蛋白记录相同与基因组/转录本边界的区分、R385化学提案与未决还原连接的限定措辞，均与上述10项原始来源审计一致。另独立核对XM_500950.3保存的CDS `/translation` 与XP_500950.3的367 aa完整蛋白记录逐字符一致。DMQ9H2/DDMQ9H2的易歧义描述已改为化学身份及分子式不同，不据缩写互换。未增加新检索、模型操作或生物学验证；终稿审阅不增加独立生物学证据计数。

主报告：`artifacts/coq3_identity_r385_proposal_20261002/REPORT.md`；SHA256：`6d7152996c84f410aa68d7255a238e5de3120870a468f205e04298b4d78cf603`。
