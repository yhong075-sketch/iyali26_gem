# R1931 的候选蛋白身份与区室证据

核验日期：2026-09-15。范围：目标的官方序列、跨版本标识、已有功能与定位注释、已有 AlphaFold 预测；未新做结构预测、未提交 HPCC、未修改模型/GPR。

## 结论

**YALI1B12993g — 原生正式名称未核实 — 推定 L-glutamate γ-semialdehyde / P5C dehydrogenase（模型/GPR赋值与自动同源注释；未取得该 W29 蛋白的实验酶学或定位证据）。** 当前 NCBI 序列明确来自 W29/CLIB89，历史 UniProt、两个项目缓存 FASTA、RefSeq、GenBank 与已有 AlphaFold 预测的目标序列均完全一致。因此这次审查没有把相似蛋白或另一菌株序列悄悄替换成目标。

功能注释支持脯氨酸分解过程中 GSA 向谷氨酸的氧化步骤。它不能支持 R1931 逆向净合成 GSA，也不能仅靠数据库等号断定体内可逆。胞质定位缺少证据；现有自动定位是**线粒体基质候选**，不能当作本物种实验定位结果。结构资料不决定方向性或定位。

## 官方记录与版本区别

- [NCBI XP_500688.1](https://www.ncbi.nlm.nih.gov/protein/XP_500688.1)：本次官方取回记录日期 2024-09-10，W29/CLIB89、染色体 1B、locus_tag YALI1_B12993g、old_locus_tag YALI2_C00310g、GeneID 2907124。PROVISIONAL，定义为 uncharacterized protein。记录明确称参考序列与 AOW01471 相同；CDD 给出 ALDH_F4-17_P5CDH 域（38–559）及推定结合/催化位点。
- [NCBI AOW01471.1](https://www.ncbi.nlm.nih.gov/protein/AOW01471.1)：W29/CLIB89、CP017554.1:1299342..1301060、572 aa、2016-10-25。基因组注释只写 hypothetical protein。记录与同源蛋白的比较说明中出现 S. cerevisiae PUT2；**不据此把 PUT2 称为已核实的 W29 原生名**。YALI0B09647g 为该记录中明确列出的旧参考比较对象；未把此字串当独立序列相等证据。
- 本次 [UniProt A0A1H6PIK2 API](https://rest.uniprot.org/uniprotkb/A0A1H6PIK2.json) 返回 Inactive/DELETED，原因为不属于参考蛋白组。此行政状态不等于功能被证伪。此前项目缓存保存的条目 v39（注释 2026-01-28、序列 v1/2017-01-18）是 unreviewed TrEMBL，蛋白存在依据为同源推断；包括 EC 1.2.1.88、脯氨酸分解步骤和 GSA 氧化反应，证据是自动规则 ECO:0000256。历史缓存不能写成目前仍有效的在线 UniProt 注释。
- [KEGG 2907124](https://www.kegg.jp/entry/yli%3A2907124) 搜索索引仍列 YALI2_C00310g、XP_500688、A0A1H6PIK2；不同数据库版本不同步。本次身份以实际取回的版本化 NCBI 记录和全序列为准。

## 区室

历史 UniProt 条目含 GO:0005759（线粒体基质），依据 IEA:TreeGrafter，属于计算传递注释。该条目没有实验 subcellular-location 评论；NCBI 的 mitochondrial precursor 描述是酿酒酵母同源蛋白的比较说明。故当前足以标出模型胞质区室需要复审，但不足以直接把反应搬到线粒体，也不证明存在胞质同工型。

## 复用的 AlphaFold 预测

[AF-A0A1H6PIK2-F1](https://alphafold.ebi.ac.uk/entry/A0A1H6PIK2)，已有模型文件 v6；API 指定 AlphaFold Monomer v2.0 pipeline，创建日 2022-06-01，本次获取 2026-09-15。API 的序列与 PDB 全部 572 个 CA 残基均与候选序列完全一致。本次是**获取已有 AlphaFold 预测**，未运行新的预测。

API 报平均 pLDDT 94.25；从下载 PDB 的 CA B-factor 实算平均 94.2620（保留原始小差异，不回填为完全相同），范围 27.75–98.94；前 30 aa 均值 35.3047。已取得逐残基置信度和完整 572×572 PAE 矩阵，PAE 文件标示最大值 31.75 Å。本次只核对序列、文件和置信度，没有做结构同源比对或底物对接，因此不新增“基于 AlphaFold 预测的功能候选”结论。高 pLDDT 不证明体内逆反应、区室、底物特异性或实验功能。

## 证据身份与限制

氨基酸大写原序列（无换行）SHA-256：`1f78e50dace6d21e517e7d1849c508e0d285ceff6d7568ff56416a7e7b93b565`。

历史全缓存 SHA-256：`e5b0a04874079b4057ffe25dadcb6b812cba8c96b227f41187187b27715753a3`。两个既有 FASTA 文件字节相同，SHA-256：`a0bacab595622864c5845411675bfb27881a4206f27572c7f0ae3abc8978e266`。路径及同一性断言结果见 `identity_checks.json`；当前 API/NCBI/AF 原文、URL、实际获取时间及各文件 SHA 见 `identity_sources/retrieval_network.json` 和 `retrieval_af_files.json`。缓存抽取条目保存在 `A0A1H6PIK2_cached.json`，候选序列保存在 `YALI1B12993g_AOW01471.1.fasta`。

第一次在受限网络读取失败已原样保留在 `identity_sources/retrieval.json`；随后使用获准的公开网络读取成功。一次本地核对尝试因未安装 Bio 模块中止，未产生科学结果；改用标准库氨基酸代码表后，全部序列相等断言通过。未安装依赖。

这些身份/区室证据可以限制方向性结论的表达，但不影响主审查用酶学与反应机制评估逆向合成的合理性。结构预测无法填补体内可逆性的实验证据。
