# R334 / R2081 独立来源审计

核验日期：2026-09-11；审计者：独立 source_audit 子任务。审计对象为本目录 REPORT.md 的核心结论及支撑它的模型、数据库、文献和序列/结构记录。独立打开原始来源，不以其他子任务的一致意见作为证据。遵循 govern-agentic-research 与 gene-identity-function；只写本文件，未改模型、整理数据、标签或阈值，未新增优化、BLAST 或结构预测。

结论：报告对“重复反应如何保留该次单敲生长”的解释有直接模型结构与历史数值支持；对生物学 OR/AND 的保留措辞适当。现有证据不能接受额外分支的独立催化功能，也不能接受两蛋白共同必需的 AND。报告中的序列相似性排序限于三个参考，不能作为全数据库功能定名。

## 对象与证据等级

| 系统 ID | 名称、简要蛋白功能及证据状态 | 模型角色 |
|---|---|---|
| YALI1F00821g；对应旧版 YALI0F00506g | 文献称 GLN1/GS1；谷氨酰胺合成酶候选。旧菌株 Q6C3E0 为 curated annotation，蛋白存在仍为同源推断；目标 W29 序列的独立酶学未核实 | R334 唯一成员，R2081 的一个 OR 分支 |
| YALI1D16151g；对应旧版 YALI0D13024g | 文献称 GLN2/GS2；统一正式名称未核实。谷氨酰胺合成酶家族候选，原生底物特异性 uncharacterized；自动注释不是独立酶学验证 | R2081 的额外 OR 分支，属于 model/GPR assignment only |

下文 F、D 仅分别简称这两行的目标基因，不代表功能已被证实。目标 NCBI 序列 AOW06441.1 / AOW03998.1 与既有 AlphaFold 输入逐一匹配；跨版本对应不等于全部蛋白序列相同。

## 原子声明审计

`supported` 表示来源支持表中限定后的准确声明；不把文献报道或旧求解结果升级为本次实验/求解复现。`partially_supported` 和 `unverified` 均属于未闭合科学判断。

| ID | 原子声明 | 直接检查的来源与条件 | 判定 | 限制或决定 |
|---|---|---|---|---|
| C01 | 当前 R334=F，R2081=F OR D | 当前 model_metadata_trna.xml 两个 FBC geneProductAssociation | supported | 原始 XML 已独立解析 |
| C02 | 两列使用完全相同六个胞质物种、计量及正向 [0,1000] | 同一 XML 的 reactants/products、species compartment 与边界参数 | supported | 无不同辅因子、区室或逆向步骤 |
| C03 | 两列原子残差为零，但当前 charge 的产物减反应物均为 −1 | 六物种 formula/charge；氮物种 H3N、charge +1，其余 charge 0 | supported | 同一化学元数据问题；不等于化学已完全验证 |
| C04 | 本组 GPR 已存在于原始输入，metadata 选择只改计量 | data/iyali26.xml；data/metadata_reaction_selection.json 的 R334/R2081 fields 与 before/after | supported | 未核实最初作者建立 OR 的原始理由 |
| C05 | 合计上界为 1000F+1000(F OR D)，F KO 后仍为1000 | C01/C02 的布尔规则与边界 | supported | 数学推导；不能推导整个可行域相同 |
| C06 | 保存的 WT 与 F KO 都有约1.8718823生长，R2081约3.305314 | 历史 results.json 的 WT 与 YALI1F00821g 原始 growth/status/fluxes | supported | 本次读取旧结果；不称重新求解 |
| C07 | 仅关 R2081 后旧记录仍生长，F KO 联合关 R2081 的旧 optimum为0 | 历史 controls.json 两条目标记录，状态均 optimal | supported | 联合关闭是基因＋反应操作，不冒充 D 单敲；全零向量本身不证明最优性 |
| C08 | 当前固定条件下重复反应保留解释了该次 F KO 非必需预测 | C01–C07；历史 reaction_snapshot.json；已有见证核验记录 | supported | 是模型机制，不是原生同工酶互换、上调或所有培养条件的结论 |
| C09 | 保存的评价记录为 F essential/FN、D unlabelled | model_static_screen_20260911/gpr_results.json 两个 gene_checks.saved_experimental_comparison | supported | 只核对已保存分类；没有重新审计原始实验 calls；D不能当负例 |
| C10 | 三份冻结 iYLI647 的 GLNS 都为 D AND F | /private/tmp/worland-fixed-audit-20260910/iYLI647_corr{,_2,_3}.json 原始反应记录 | supported | 三版本沿袭是模型赋值，非三份独立生物学证据 |
| C11 | 仅改 R2081 为 AND 仍留下 D KO 的 R334；去重复保留 OR 仍留下 F KO 的 D 分支 | 当前布尔规则、同计量列和边界 | supported | 合并为单列[0,1000]还改变总容量；不能仅靠一个FN翻转验收 |
| C12 | Q6C3E0 为 reviewed，但当前记录不提供目标W29直接酶学验证 | sources/Q6C3E0.json：CLIB122/E150、proteinExistence=3、唯一引用为大规模基因组；催化及区室注释 | supported | reviewed不等于目标蛋白实验验证 |
| C13 | Q6C990 主要为自动/同源家族与域注释 | sources/Q6C990.json：unreviewed、proteinExistence=3、ARBA/PROSITE、无独立催化反应证据 | supported | 没有实验支持不等于排除 GS 活性 |
| C14 | 两个目标 UniProt Inactive 理由是参考蛋白组收录政策 | sources/A0A1D8NLC7.json、A0A1D8NED6.json 的 inactiveReason | supported | 不能解释为基因删除、假基因或酶失活 |
| C15 | 目标NCBI序列与对应既有AlphaFold输入完全相同，长度374/458 | 两FASTA；两AF API JSON；PDB；check_sequences.py 的坐标序列/连续编号断言及输出 | supported | 独立比较了FASTA与API序列，并核对PDB CA数量；完整坐标序列断言按脚本静态审查，不重新预测 |
| C16 | F与CLIB122序列不完全相同，D与Q6C990完全相同 | 原始FASTA与Q6C3E0/Q6C990 sequence.value 独立比较 | supported | F按C端偏移10比较时，另有目标11/12位差异；不能将F写成100%序列同一 |
| C17 | 复用的是已有AlphaFold Monomer v2.0 pipeline、数据库文件v6预测 | 两AF API JSON 的toolUsed/modelCreatedDate/latestVersion及PDB/PAE URL | supported | modelCreatedDate=2022-06-01；本次取得，不是本次预测；单体模型不证明复合体 |
| C18 | PDB CA均值与API globalMetric有小差异，PAE统计按保存数组可复算 | 独立读取两PDB CA B-factor和PAE数组；target_sequence_structure.json | supported | F 94.8421123 vs 94.81；D 94.9773799 vs 95.0。保留不同统计来源，不强制相等、不换输入 |
| C19 | 局部BLAST严格限于3参考，4个命中统计与XML一致 | sources/local_blast.xml、reference_panel.fasta、blast_execution.json；独立重算HSP identity及并集覆盖并核对SHA | supported | BLASTP 2.17.0+；E≤1e-5；未重跑BLAST；E值依赖此小面板 |
| C20 | F→酵母GS约76.9%一致/97.1%query覆盖；D→不同底物参考约26.4%/96.9%，D→细菌GS仅约48.0%query覆盖 | 同一原始BLAST XML；P32288/P0A9C5/P78061原始参考序列 | supported | D无达到所用阈值的P32288命中；仅反对“名称即可确认同工酶”，不能定名D为PuuA或排除GS |
| C21 | 2017论文称两旧版ID为gln1/gln2，并观察氮源相关转录响应 | 独立打开ASM正文“gzf2 is required…”、Fig7图注、Table2及方法 | supported | Fig7是RNA-seq叠图；Table2删除对象为转录调节因子，非两个GS的互补/酶学实验 |
| C22 | 2017相关转引不是这两个Yarrowia GS的独立酶学验证 | 同文ref54；独立打开Miller & Magasanik 1990的ASM摘要 | supported | 该引文为S. cerevisiae NAD-linked glutamate dehydrogenase研究；物种、对象不同 |
| C23 | 2016研究提供氮限制下Gln1/2上调的多组学观察 | 独立打开PMC4766638 Nitrogen assimilation及相关引文列表 | supported | 未重新核验独特肽段/补充表；上调不决定OR/AND；59–61不能移植为本对精确蛋白的酶学 |
| C24 | 2026新近研究确有GS相关表达盒工程线索，不能概括成从无工程研究 | 独立打开preprint v1；独立检索正式publisher正文；本目录2026_mdpi_web_search.txt及获取失败记录 | supported | 正式版为检索取得的正文片段，未冒充直连全文核验；工程细胞/粗蛋白终点不是独立催化 |
| C25 | R334的F单基因功能赋值有较强间接支持 | C12/C16/C20/C21/C23 的注释、序列及表达证据 | partially_supported | 可支持排序和候选保留；目标W29的精确反应、定位和容量未直接验证 |
| C26 | D能在相关条件下独立完成该胞质GS反应，足以接受OR | 当前已核来源没有精确单独催化、遗传互补及定位闭合证据 | unverified | 报告正确保留为未决，未将其写成假或失功能 |
| C27 | F和D共同必需，足以接受AND | 当前已核来源没有二者共同催化/复合体依赖证据 | unverified | 外部AND及共表达不能补此缺口；不是证明AND永不可能 |
| C28 | 2026最终版本的GS1/GS2分别构建、表达和酶活已被核实 | 未取得并核实正式supplement与完整修订材料 | unverified | 报告明确保留此缺口；不能用预印本默认替代正式附录 |

## 覆盖率及限制

```text
total claims 28 | audited 28 | supported 24 | unresolved 4 | contradicted 0 | unchecked 0
unresolved = partially_supported 1 + unverified 3
```

这里 audited=28 表示对全部28个声明做了来源匹配或明确缺口判定；**不是28个科学声明全部证实，也不是所有补充资料已读**。已支持的24项多为静态模型事实、历史记录、原文实验范围或序列统计。没有发现与报告限定后主张直接相反的证据；OR/AND的关键未决项仍保留。未扩展成全模型化学审计、全数据库同源搜索、完整历史环境重建或独立实验复现。

序列检查初版把PDB CA均值与API globalMetric强制限定差异<0.02，审计独立观察到两目标都超过该阈值。修订代码及JSON改为并列保存两种统计原值和差值；没有变换输入或改变生物学验收阈值。此为来源统计口径差异，不是低置信度或失功能证据。BLAST XML的4个HSP命中、identity、query/subject覆盖及三份输入/输出SHA已逐一静态核查。

候选后续若变更模型仍须人类授权，并对所有重复表示、容量和相关正确预测共同验收。当前审计可支持交付限定审查报告，不支持自动接受新的GPR或化学规则。

## 来源与身份

- 模型：model_metadata_trna.xml SHA256 `d274bad3050e3c9220a8b6287eae847f3bf1334892284d565a6c4d96b38135a0`，本审计独立计算匹配。
- 历史原始记录：`artifacts/iyli647_screen_20260910/nonessential_diagnosis_20260911/results.json`、`controls.json`、`reaction_snapshot.json`；未只读模型子任务的总结。各文件完整SHA、培养、菌株、求解器和dirty环境详见MODEL_AUDIT.md与原始记录。
- 目标序列SHA256：AOW06441.1 `074de8ed73d25ead502698631e943f35624202e5bffbf958561315e0dd37f354`；AOW03998.1 `11c383111f92c305adb675026cdbd0133b22bfdc5f02e5700bbdd131050ef123`；本审计从原FASTA独立计算。
- [2017 ASM原文](https://journals.asm.org/doi/10.1128/msphere.00038-17)、[2016 PMC原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4766638/)、[1990转引的原始摘要](https://journals.asm.org/doi/10.1128/jb.172.9.4927-4935.1990)、[2026预印本v1](https://www.preprints.org/manuscript/202604.1894)、[2026正式版检索来源](https://www.mdpi.com/2311-5637/12/7/315)。未用ResearchGate摘要或搜索排名补充直接实验结论。
