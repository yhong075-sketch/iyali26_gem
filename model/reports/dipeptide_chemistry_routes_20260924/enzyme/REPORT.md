# 酶身份与定位增量审查

核验日期：2026-09-24。新增工作为公开版本序列和 PO1f 映射复核、差异的结构域定位、同源定位冲突检查，以及精确游离二肽的定向检索。未运行新 BLAST、AlphaFold、对接或优化，未修改 GPR、反应区室或旧记录。

## 新增、复用和未决项

**本轮新增的最实质信息是两个 PO1f 对应蛋白已由官方序列记录闭合。** YALI1E16433g（原生正式名未核实；M24B prolidase-like 候选）对应 PO1f **BKA90DRAFT_141684 / KAE8169967.1**；YALI1F23706g（原生正式名未核实；DUG1-like/M20A Cys-Gly 金属二肽酶候选）对应 PO1f **BKA90DRAFT_151083 / KAE8174307.1**。两者的 PO1f ORIGIN 全序列分别与 W29 当前记录完全相同。映射来自 NCBI Identical Protein Groups，随后独立下载有版本的 GenPept 并逐字比较，未用菌株亲缘代替序列证据。此结果支持使用固定序列研究 PO1f 候选，不能转移当前培养中的表达量、活性或定位。[PO1f 454 aa](https://www.ncbi.nlm.nih.gov/protein/KAE8169967.1)、[PO1f 478 aa](https://www.ncbi.nlm.nih.gov/protein/KAE8174307.1)

| 候选、原生身份与证据级别 | 实际序列比较 | 新增差异解释 |
|---|---|---|
| YALI1E16433g；正式名未核实；M24B prolidase-like，自动/同源功能候选 | W29 AOW05368.1 = XP_503902.2 = PO1f KAE8169967.1，454 aa；CLIB122 YALI0E13464g/CAG79495.1 为单处差异 | W29→CLIB122 **Q17K**。位置7–27分别为 `PAKAHALKAAQHLKASGASDD` / `PAKAHALKAAKHLKASGASDD`。Q17在 NCBI AMP_N 6–131（UniProt另一域边界6–141）中，不在 Prolidase159–427，也不在注释活性位点221/238/249/333/373/412；六位点全保守。域边界采用各数据库原值，不强行合一。没有实测排除Q17K对结构、稳定性或定位的影响。 |
| YALI1F23706g；正式名未核实；DUG1-like/M20A Cys-Gly，自动/同源及既有预测结构候选 | W29 XP_505554.3 = AOW07335.1 = CLIB122 YALI0F17842g/CAG78363.1 = PO1f KAE8174307.1 = 既有AF坐标序列，478 aa | 历史 XP_505554.2 明示 DSM3286、477 aa，末端重复KK少一个K；可等价表示为删除W29位置476或477，不能选择唯一位置。末端`TLGAYLHYIAEEKKE` 对 `TLGAYLHYIAEEKE`，在M20二聚化域218–368及已注释活性/结合位点98–447之外；“不在注释位点”不等于无功能影响。 |

完整序列、SHA、版本、原始文件SHA和断言保存在 `sequence_comparison.json`；新运行检查是 `verify_sequences.py`。复用的 **AlphaFold 预测** 为 AF-A0A1D8NNX0-F1-model_v6，原元数据记录 AlphaFold Monomer v2.0 pipeline，模型创建2022-06-01、既有获取日期2026-09-24。此次重新解析478个Cα残基，确认与当前W29/PO1f完全一致；Cα平均pLDDT96.593。旧审查已取得PAE，本轮不重复计算或下载；单体置信度不确认底物、天然二聚体、金属占据、区室或催化端朝向。检查范围内没有E16433的既有目标结构；未把CLIB122结构冒充W29，也未启动新预测。

## 定位和具体底物：新增结果及裁决

新读取的酿酒酵母同源 **YFR006W / P43590**（已审阅名称仍为“Uncharacterized peptidase YFR006W”，M24B功能候选）含两类不同层次的定位信息：GFP高通量 **cytoplasm** 注释（GO HDA，PMID14562095）与**预测**跨膜螺旋8–24、单跨膜蛋白注释。此处保留证据张力，不把它简化成“已证明可溶胞质酶”，更不将同源的膜区、液泡腔或催化端朝向转移至W29。W29/PO1f E16433的原生定位、拓扑和当前SD-Leu中的活性仍未确认。F23706既有cytoplasm IEA与酿酒酵母 **YFR044C / DUG1**（已实验表征的Cys-Gly二肽酶）的胞质荧光，仅支持胞质候选；未新增W29/PO1f定位证据。[P43590官方记录](https://rest.uniprot.org/uniprotkb/P43590.json)、[GFP原始研究](https://pubmed.ncbi.nlm.nih.gov/14562095/)

**新增精确底物来源：** Ito等2003在 *Actinomucor elegans* IFO6408纯化的glycyl aminopeptidase上报告游离 **Gly-Glu** 释放Gly；原文p84以ninhydrin与TLC/标准物检测肽水解，p85明确列出Gly-Glu。已保存原PDF并目视核对p85。该文未提供可据此映射到本项目的Yarrowia候选位点，本轮没有建立同源关系；它是另一可用的酶学比较来源，不是F23706或R2029液泡GPR的证据，也没有从Gly-Phe的数值替代Gly-Glu速率。[原始论文，DOI10.1271/bbb.67.83](https://www.jstage.jst.go.jp/article/bbb/67/1/67_1_83/_pdf)

对 **Gly-Asp、Gly-Glu、Ala-Gly** 的本轮精确名称/顺序及候选ID定向检索，未找到W29/PO1f基因到这些自由底物的直接实测链接。Gly-Glu的上述真菌来源使后续筛选不必限于DUG1，但不足以新命名YALI1催化蛋白。酿酒酵母DUG1既有2009表2是含Cys底物面板，四目标均未测试；未测试不写成无活性。**Gly-Pro** 保留E16433为优先的序列明确、催化类别匹配的假设；这可以支持分级的同源推断提案，不要求只能等待原生实测，但现有液泡区室和底物特异性合并证据仍不足以接纳GPR。此次没有自动填入任一反应。

复用9月11日DPP资料时区分反应类型：**YALI1B04274g**（原生名未核实；DPP-IV-like候选，基于序列/既有AlphaFold）与 **YALI1B25603g**（原生名未核实；DPP-III-like候选，基于序列/既有AlphaFold）研究的是从较长肽释放N端二肽，不能充当四条自由二肽水解。Ala-Gly的第二残基是Gly；X-Pro与Pro-X不能互换。Yarrowia前体加工的1989/1997原始研究不确认这四个游离产物在液泡中的定量供给；2007可溶提取物研究是在不同营养条件下测整体肽酶，较高活性见于peptone/静止期，不能替代PO1f/SD-Leu中这两个候选的表达或活性测量。前述原始摘要本轮重新取回；未用表达/家族标签推断通量。[1989](https://pubmed.ncbi.nlm.nih.gov/2649495/)、[1997](https://pubmed.ncbi.nlm.nih.gov/9353927/)、[2007](https://pubmed.ncbi.nlm.nih.gov/17227470/)

四反应的六项分别裁决见 `../reaction_identity_cards.md`，逐轴表见 `../enzyme_localization_evidence.tsv`。检索范围、未取得项、原始文件和复用来源见 `source_manifest.json`、`search_log.json`。当前培养表达、原生定位/朝向、目标底物活性均是“尚未确认”，不是负结果；未向外部服务上传私有序列。

## 独立来源审核

本轮酶声明 EN1–EN10 共10项，独立审核10项；限定范围内支持10、未解决声明0、反驳0、未核查0。审核重新解析GenPept/UniProt/既有PDB，并核对原始Ito论文p84–85、DUG1表2和原始PubMed摘要；EN10只认可记录的有限检索范围。这里的“10项支持”包括对证据边界的支持，不将原生活性、液泡定位或供给等未决机制计为已验证。
