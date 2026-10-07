# 独立来源与结果审计

核验时间：2026-09-25 UTC（2026-09-24 America/Los_Angeles）。范围继承上级 `TASK.md`：优先核查相反证据、原始基因行和序列版本；只写本子目录，没有模型/GPR变更、外部通信、集群作业，也没有重复运行定位预测或BLAST。所用技能为 govern-agentic-research、gene-identity-function、pdf；结果解析采用 ponytail 的标准库最小核验方式。

| 对象 | 名称／功能与证据等级 |
|---|---|
| YALI1E16433g / W29 AOW05368.1 | 原生正式名未核实；M24B/X-Pro 肽酶功能候选，原生活性及区室未实验表征。固定454 aa，完整SHA见独立结果文件。 |
| S. cerevisiae YFR006W / P43590 | 无已确立专名；未表征M24B肽酶（审阅注释），有标记蛋白定位和蛋白质组检出证据，不能据此确认游离Gly-Pro酶活。 |
| A. nidulans AN5810 / pepP | 作者实验表征的X-Pro二肽水解酶；实验关联CAC39600.1为496 aa，当前审阅Q96WX8为465 aa；原生细胞器定位未直接测定。 |

**酿酒酵母参照的定位证据确实存在两个层次，不能简单写成“只有胞质”或“液泡常驻”。** [Huh原项目实际记录](https://yeastgfp.yeastgenome.org/displayLocImage.php?loc=640)的初评及最终定位均为cytoplasm；本次读取该基因行及全视野GFP图，记录为SD/G1、均一分布、2秒曝光。原文方法是染色体C端GFP融合、SD缺His中对数期成像，作者提示标签可能影响定位。[Huh原文](https://yeastgfp.yeastgenome.org/nature02026_r.pdf) 不应把网页的“molecules/cell: not visualized”误读成GFP无信号：这是独立的丰度字段，定位图片与定位调用存在。

[Sarry2007](https://www.mcgill.ca/parasitology/files/parasitology/Dzierszinski2.pdf) **Table3打印页4294实际列出YFR006W/P43590的一个2-DE/MS鉴定斑点**。这不是仅引述数据库的比较行；但同表GFP/SwissProt/SGD定位列是数据库导入。实验为SEY6210、YPD、30°C、OD600约1.8的液泡腔可溶组分；作者把此蛋白归入可能为回收/降解而进入液泡的非典型蛋白酶。Table2的定量蛋白酶K保护面板没有YFR006W，不能声称其单独通过了该实验。补表S1未取到，所以精确肽段证据、该斑点实验分子量和完整性未核验；不能据此判断常驻或活性。

[Jalving原始研究作者存档章节](https://edepot.wur.nl/121628)打印页85（PDF页91，已视觉检查）测得细胞提取液能水解Met-Pro、Ala-Pro、Gly-Pro，培养滤液未检出；其“胞质”判断明确来自缺少已知分选/TM序列的推断。后续过表达和纯化把活性与pepP关联，但细胞提取不区分胞质和细胞器。当前Q96WX8逐位等于CAC39600.1的32–496位，31 aa起始差异与当前记录的错误起始警告相符；不能无条件把原496 aa扫描结论转给当前序列。

另一条反证线索是[Hitchcock2003](https://pmc.ncbi.nlm.nih.gov/articles/PMC240687/)膜相关泛素化蛋白质组：UniProt引用该文，主文也已读取，但其确切YFR006W补表行没有成功取得。三份下载均为访问挑战页面，故这一基因特异的膜组分结论仍未核验，不能以文章标题代替证据。

**四项预测与限定五参照比对已独立解析核对通过。** DeepLoc原始CSV为Cytoplasm|Nucleus、Soluble，全454残基注意力文件吻合；SignalP和TargetP官方OTHER分数、参数及各70/200残基前缀吻合。两者下载不能独立证明其余输入残基，完整提交来自主执行者的输入记录。DeepTMHMM官方ZIP中元数据证实1.0.57、exit0，完整输入/输出均吻合454 aa，GLOB、0 TMR、454个I状态；I不能升级成胞质或核定位实验。

DeepLoc数值阈值最初仅有网页手工记录；随后取回官方results.json，本次已独立解析并核对全部14项阈值与返回分类，原来源缺口已关闭。五参照BLAST XML/TSV全部统计一致；恢复显示为SEG遮蔽X的原始残基后，重新数出的同一残基数也一致。P43590比对仅目标42–454／参照111–535，219/425=51.5294%；其预测TM8–24在比对之外。该限制充分阻止凭此比对转移N端膜特征，但并不证明其他类型膜结合不存在。此次是本地BLASTP 2.17.0+限定面板核验，绝非仍WAITING的在线Swiss-Prot结果。

审计覆盖：**20项命题｜18项已审计｜13项支持｜3项未解决｜2项被反证｜2项未核查**。`claims.tsv`包含支持命题及被检验的过度声明，不能把“被反证命题”误认为主执行者曾作出的结论。未核查两项是Sarry补表细节和Hitchcock确切基因行；未解决三项是排他原生胞质、液泡常驻活性、跨物种直接转移这三种超出证据的强断言。

可复核入口：`python artifacts/glypro_localization_predictors_20260924/reference_audit/audit_results.py`。原始来源及派生图像在`raw/`，来源/完整SHA在`source_manifest.json`，结果/实际核查文件SHA在`independent_results.json`，覆盖在`coverage.json`。本次没有独立重做旧AlphaFold结构叠合、原生序列跨株全审计或酶学复现。
