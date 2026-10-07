# 三靶点身份、旧新版映射与实验标签审查

核验时间：2026-09-11 UTC。本轮直接读取原始来源及仓内快照，未求解模型、修改GPR或更改实验标签。

**发现一个会改变F位点结论的具体冲突：YALI0E16192g并非没有W29对应记录；原始S2表将它对应到YALI1_E19360g，但将W29位点列为Pseudo mRNA。NCBI当前记录进一步注释为移码导致的伪基因。** 因此，可以记录这条跨菌株位点映射，不能把旧CLIB122的完整F候选蛋白视为W29相同蛋白，也不能仅据该旧蛋白的BLAST或AlphaFold证据将W29对应位点采纳为共同AND。

| 用户所指位点 | 原生正式名称、功能与证据等级 | 本轮身份/标签判定 |
|---|---|---|
| YALI1D00581g | 原生正式名未核实；V1 D中央转轴亚基候选；自动同源注释，非原生功能实验 | S2行3747对应旧YALI0D00583g。已检查1612正例工作簿，无该ID；未恢复逐assay calls，不判为实验非必需 |
| YALI0E16192g | 原生正式名未核实；CLIB122的V1 F中央转轴亚基候选；Q6C5Q2 v1，122 aa，自动家族注释 | S2行6015对应W29的YALI1_E19360g；W29记录为Pseudo mRNA/NCBI伪基因。旧ID及新ID均不在已检查正例工作簿 |
| YALI1F38820g | 原生正式名未核实；V0 a亚基候选，参与膜内质子通道和复合体装配；同源/预测结构支持，原生区室未定 | S2行8659对应旧YALI0F31119g。已检查正例工作簿，无该ID；未恢复逐assay calls，不判为实验非必需 |

YALI1_E19360g（模型格式YALI1E19360g）同样没有本轮核实的原生正式名；记录将它描述为旧F候选的截短对应位点。其功能状态的证据级别是**数据库注释**，不是本轮或来源中已定位到的蛋白功能实验。

## 原始来源核实

1. **2016年原始S2映射表。** 从PLOS出版商取回S2表，完整字节与仓内`S2_table_YALI1_YALI0_map.xlsx`一致。`YALI1 Genes`工作表A6015:H6015列出旧ID、新ID、`Pseudo mRNA`及W29坐标1935676–1936023、负链。原文说明S2用于CLIB89与CLIB122的编码及其他转录位点对应，不能将所有行当作功能完整的蛋白一对一映射。[原始论文与S2定位](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0162363)、[S2原件](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0162363.s004&type=supplementary)。

2. **NCBI目标记录。** GeneID2911580的类型是`pseudo`，locus tag为YALI1_E19360g。NC_090774.1局部GenBank记录明确菌株CLIB89(W29)、同一坐标、`/pseudo`，并写明`nonfunctional due to frameshift`。记录最后更新于2024年9月；本轮获取于2026-09-11 00:02 UTC。该记录来自提交者注释传播，RefSeq注明尚未最终审核；它与2016年S2一致，但两者可能共享注释来源，不能视为两次独立功能实验。[NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/2911580)、[版本化染色体区间记录](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_090774.1&seq_start=1935576&seq_stop=1936123&rettype=gbwithparts&retmode=text)。

3. **序列身份的界限。** 本轮核实封存Q6C5Q2 v1为CLIB122/E150、122 aa。所取W29位点的基因组区间为348 nt，在该目标feature中没有CDS、蛋白accession或翻译序列。348 nt是基因组区间长度，不能自动除以3写成116 aa的“W29 F蛋白”。没有对位点重新预测ORF、修复移码或替换分析序列；两者完整蛋白序列相同未被确立。

4. **此前派生映射为何遗漏。** `iyli21_genes_vs_S2.csv`给旧YALI0E16192g标为`not_in_S2`，而完整原始S2的上述行确实存在。该派生表中的缺失不能作为“没有旧新版对应”的证据。本轮保留原表和派生表原值；没有修改映射数据。

## 正例与逐assay calls分开

本轮直接检查了**已经放入本仓库**的`screen_test_metadata_trna_20260910/inputs/42003_2023_4996_MOESM10_ESM.xlsx`，Sheet1中1612条ID。仅去掉YALI0/1/2前缀后的下划线；不应用跨版本映射或别名。三靶点及W29对应新ID均无行，因而只得到“**不在这份正例列表**”。

同时，已有`essentiality_per_gene.tsv`三靶点均为`experimental_label=unlabelled`；本轮读取这一交付物，未重新运行敲除。Cas9、Cas12a、transposon的原始逐assay calls未在此次限定检查中恢复。工作簿正例身份、图中灰色、或者模型预测非必需，均不能补成实验nonessential标签。S2派生表里的`metabolic`和`in_model`也不是实验calls。

## 对本轮决定的作用

- **可以支持：** 完整F组件的机制必要性可以继续用经核实的参照蛋白研究；旧CLIB122位点及其旧结果保留原身份。W29新位点作为“注释伪基因对应”单独记录，不能作为已确认活性蛋白替代。
- **目前不能支持：** 将YALI0E16192g自动改名成YALI1E19360g、把旧F预测结构当W29功能证据、因移码注释就断言W29细胞没有F功能、或把任一未知标签改成essential/nonessential。
- **解除F基因层面阻碍所需证据：** 目标W29/PO1f实际等位序列与转录本能否形成完整F蛋白，或可承担该组件的其他位点。当前注释冲突不能用旧蛋白的高相似性解决。

本子任务的8项原子主张及来源登记在`evidence.json`，初始审计状态为unchecked，交由独立来源审计。可运行`verify_mapping.py`复核S2字节、精确目标行、NCBI类型、封存序列与正例表中目标缺席；检查通过。完整输入/代码SHA与原始获取记录保留在`verification.json`和`sources/fetches.json`。
