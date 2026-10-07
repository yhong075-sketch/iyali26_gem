# C. phangngensis 参照检索：独立注释核查

核查时间：2026-10-02 22:09–22:12 UTC。范围：只读文献、公开数据库及版本化基因组目录；最多 6 批检索，实际 5 批、20 条查询，另有来源打开与数据库只读请求。未做 BLAST、序列比对、结构预测、优化、模型或 GPR 修改、集群作业。该文件仅记录本轮证据，不是新的项目状态账本。

## 结果

1. **生物学参照成立，末端化学并未被该物种证据解决。** Limtong et al. 2008 的原始物种论文摘要明确报告两株酵母的主要泛醌为 Q-9，且归入 Yarrowia clade。摘要没有确定末端羟化自由底物氧化态、电子供体、净反应或末端酶位点。[PMID 18218960](https://pubmed.ncbi.nlm.nih.gov/18218960/)，DOI 10.1099/ijs.0.65506-0，Abstract。
2. **名称锚点明确，收藏号不能全部混同。** 当前 NCBI taxon 444778 为 Yarrowia phangngaensis；记录关联旧 Candida 名称及原始论文。PubMed 更正摘要与出版者旧摘要、2026 命名论文之间的 BCC/NBRC 编号不一致仍保留；不能据此断言所有收藏号或两个组装字节相同。[NCBI taxonomy 444778](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=444778)。
3. **可获得物种参照组装，但未核定目标酶。** CBS 10407 的 GCA_900519005.1 在 Brinkrolf et al. 2021 的基因组比较方法中明确列出；NRRL Y-63743 的 GCA_030581735.1、WGS JAKTVT000000000、BioSample SAMN20341414 在 NCBI BioProject 中明确列出。它们只是后续可选择的版本化参考，不能自动指定为本项目输入或赋予反应/GPR。[Brinkrolf 2021，Comparison of genome structures and gene repertoires](https://link.springer.com/article/10.1186/s12864-021-07597-z)；[PRJNA736342，Y. phangngaensis 行](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA736342)。
4. **UniProt 当前覆盖不足以给出末端位点。** taxon 444778 全量 REST 查询实际返回 24 条：15 条线粒体呼吸/核糖体相关蛋白及 9 条推定糖转运蛋白；没有命中目标晚期羟化或 O-甲基化蛋白。release 2026_03，release date 2026-09-02，X-Total-Results=24。该结果只说明该查询的公开注释覆盖，不证明该物种缺少有关基因。原始 TSV 与响应头保存在同目录。
5. **JGI 为可追查注释入口，尚不是已核实位点。** 搜索索引显示 Yarpha_NRRLY63743_1 基因组门户具有 5775 个基因模型，来自 Y1000+ 工作并追加 JGI 功能注释；网页直开失败，因此这里只保留线索，不将具体 COQ 注释视为核实。[JGI organism page](https://myco-lb.jgi.doe.gov/Yarpha_NRRLY63743_1/Yarpha_NRRLY63743_1.home.html)。
6. **KEGG/BioCyc 不作缺失断言。** KEGG 只查到物种分类页，未查实物种通路或蛋白位点；REST `/list/organism` 为 HTTP 400。BioCyc 关键词未返回目标记录，尝试页面不可访问。访问失败不等于数据库不存在该物种或反应。

未能给出 R695 或 R385 的 C. phangngensis 系统基因 ID、已核实符号或原生蛋白功能；因此不生成候选 GPR。R695 自由 DMQ9/DMQ9H2 底物状态、供体与产物状态，R385 原生末端酶的底物特异性，均仍需独立直接证据。普通途径图、Q-9 终产物测定、基因组存在三者均不能独立填平这些缺口。

## 检索与访问边界

五批分别覆盖：

- 新/旧物种名 + coenzyme Q / ubiquinone / UniProt / genome.jp。
- 新名 + genome / methyltransferase / UniProt；旧名 + COQ。
- 新名 + COQ7 / COQ3 / BioCyc；旧名 + ubiquinone + pathway。
- 新名 + hydroxylase / methyltransferase；旧名 + COQ7 / COQ3 / demethoxyubiquinone；JGI portal + coenzyme。
- 新名在 genome.jp、biocyc.org 的限制域检索，JGI portal + ubiquinone，以及版本 GCA_900519005.1 + CBS。

COQ3/COQ7 在这些检索中仅作为跨物种功能关键词，未赋予目标物种具体基因。未发现物种特异性酶学实验，不等同于不存在此类实验。

访问失败：web 工具无法读取 UniProt REST（随后 curl 成功）；KEGG dbget 查询不可访问，KEGG REST HTTP400；JGI 两域直开 cache miss；BioCyc 试探性 URL（未核实对象号）不可访问；PMC8091737 一次 recaptcha（随后正式出版者全文成功）；NCBI datasets 两个基因组页面动态壳/失败（改用已打开的原始论文和 BioProject 目录核验）；初次请求误写另一 BMC DOI 后失败，已弃用，并核准正确 DOI 为 10.1186/s12864-021-07597-z。

UniProt 请求：
`https://rest.uniprot.org/uniprotkb/search?query=organism_id%3A444778&format=tsv&fields=accession,id,protein_name,gene_names,organism_name&size=500`

保存时刻约 2026-10-02 22:11:38 UTC；HTTP 响应 Date 为 2026-10-02 22:10:29 GMT（缓存响应）。TSV SHA-256：`ff65728d32645c1e7305cecb01a953f0ae60ad044e06d4019043d67e365776a0`。未下载/执行基因组或候选蛋白序列，因而无本轮序列 SHA、BLAST 或结构结果。
