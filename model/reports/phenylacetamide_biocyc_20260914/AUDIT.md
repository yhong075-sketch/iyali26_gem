# CPD-238 化学身份独立核查

核查日期：2026-09-14。范围仅为3条权威化学记录：ChEBI:16562、NIST 103-81-1、NIST 103-84-4。读取网页的身份、结构字符串、别名和登记号字段；未扩展代谢路径、基因或酶学调查，未进行模型/科学数据/代码变更或计算作业。

**结论：2-phenylacetamide 与 N-phenylacetamide（acetanilide）是不同连接结构的化合物，不能互作别名。** 它们分子式均为 C₈H₉NO、相对分子质量均约135.1632，因此仅凭分子式或质量无法区分。

| 身份 | 连接结构 | CAS | InChIKey | 直接来源 |
|---|---|---|---|---|
| 苯乙酰胺；2-phenylacetamide、phenylacetamide、benzeneacetamide、phenylacetic acid amide | Ph–CH₂–C(=O)–NH₂；芳环在酰基侧，隔一个CH₂ | **103-81-1** | **LSBDFXRDZJMBSC-UHFFFAOYSA-N** | [ChEBI:16562](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:16562)、[NIST苯乙酰胺记录](https://webbook.nist.gov/cgi/cbook.cgi?ID=103-81-1) |
| 乙酰苯胺；acetanilide、N-phenylacetamide、N-acetylaniline | CH₃–C(=O)–NH–Ph；芳环直接连酰胺N | **103-84-4** | **FZERHIULMFGESH-UHFFFAOYSA-N** | [NIST乙酰苯胺记录](https://webbook.nist.gov/cgi/cbook.cgi?ID=103-84-4&Units=SI) |

ChEBI:16562 的定义、SMILES `NC(=O)Cc1ccccc1`、完整InChI及InChIKey共同支持第一行。该记录还明确把 **MetaCyc CPD-238** 和 **KEGG C02505** 列为交叉引用；这证明ChEBI的对应关系，不能代替逐条审查MetaCyc或KEGG所有关联反应。定位：ChEBI网页的Definition、InChI/SMILES、Manual Xrefs、Registry Numbers；NIST两页均为页面顶部的InChI/InChIKey、CAS和Other names字段。

主代理本次浏览MetaCyc 30.0后转录：CPD-238标题、结构和InChIKey属于第一行，但别名含N-phenylacetamide/acetanilide，CAS为103-84-4。**若该转录准确，则这些字段内部不一致：后两种别名和CAS属于第二行。** 这是用独立化学来源核验转录字段得出的条件性判断；本审计没有再次独立打开MetaCyc页面，因此不把主代理观察写成本审计直接复核的页面内容。不能据此把该页面的全部反应都判为错误，也不能把乙酰苯胺相关反应当作苯乙酰胺的来源或去路。

本次登记7项：5项有直接来源支持（苯乙酰胺身份、其ChEBI/NIST一致性、乙酰苯胺身份、二者连接结构不同、ChEBI列出的交叉引用）；2项未独立确证（MetaCyc页面实时字段、PubChem 7680的独立记录）。未决项不影响CAS和两个化合物的区别。未增加第四条化学来源。

访问限制：首次NIST `ID=C103811/C103844&Units=SI` URL未被检索工具打开；随后通过ChEBI提供的103-81-1链接及NIST自身索引的103-84-4链接，成功读取完整身份页。ChEBI检索返回页标注Last Modified为2020-07-21，本审计没有据此宣称其为数据库最新发行版。后续说明苯乙酰胺来源/去路时，应以结构、反应计量及原始酶学核对，而不是混杂的别名。

## 追加：上游反应与已发表酶学

同日追加范围限定3个来源：Expasy EC 1.13.12.9和PMID6501250、PMID8655547的原始摘要。PubMed网页未返回摘要或出现浏览器检查，因此通过Europe PMC的同PMID记录实际读取两份原摘要；没有据此声称已读全文或底物数据表。

**上游候选路线有已知酶学背景：**

L-苯丙氨酸 + O₂ → 2-苯乙酰胺 + CO₂ + H₂O

[Expasy EC 1.13.12.9](https://enzyme.expasy.org/EC/1.13.12.9)明确给出该反应，酶名为phenylalanine 2-monooxygenase。该条目同时列出生成苯丙酮酸、氨与过氧化氢的氧化支路；所述约80%/20%是酶条目的反应分支说明，不能用作所有条件的固定比例，更不是W29预测通量。

[PMID6501250原研究](https://doi.org/10.1093/oxfordjournals.jbchem.a134853)的题名明确菌株为 **Pseudomonas sp. P-501**。原摘要描述纯化酶、底物测试、HPLC产物分析及L-苯丙氨酸氧化/氧合的动力学测定。因此，上游存在该细菌来源的直接酶学研究，而不是仅凭反应数据库的名称推想。摘要自身没有逐项重述上述PAM反应式；本核查将反应式的EC证据与该原研究的实验范围联合使用，保留来源层级。

**下游文献确有纯化酶学，但本次摘要覆盖不到PAM专项数据：**

[PMID8655547原研究](https://doi.org/10.1128/jb.178.12.3501-3507.1996)研究 **Rhodococcus erythropolis MP50** 的酰胺酶，摘要记载纯化至均一并水解多种脂肪及芳香酰胺；列出的具名实验底物包括2-phenylpropionamide及两种药物酰胺。**2-phenylpropionamide不是2-phenylacetamide**：前者多出α位甲基。该摘要没有单独列出PAM，因而仅依据本次所读摘要，不能独立确认RXN-12492的PAM底物测定、参数或可逆性；全文底物表未查不等于不存在该实验。

主代理观察MetaCyc两反应页“no enzyme identified in this database”，属于数据库本条目的酶/基因关联状态。上述原摘要已有纯化酶研究，因此不能把这种关联缺失解释为科学界没有鉴定酶或没有酶学证据。相应地，页面的“方向未知”不证明真实生理反应可逆。主代理转录的MetaCyc实时反应清单及方向字段仍未被本审计直接复核。

追加登记6项：4项有支持（EC上游反应、P-501纯化酶研究、MP50纯化酶研究、数据库关联缺失不等于没有酶学）；2项未确证（仅凭本次摘要确认PAM专项水解数据、把这些异物种反应确认为W29路线）。与前段合计13项：9项按注明证据层级受支持、4项未决。**可以把L-苯丙氨酸氧合列为PAM的一条已知候选来源路线，不能写成W29已经具备或本模型应当新增的反应。** 本次没有追加序列、定位、结构、求解或模型修改。
