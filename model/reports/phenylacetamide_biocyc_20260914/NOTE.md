# MetaCyc CPD-238 定向补查

日期2026-09-14。用户提供苯乙酰胺化合物页，范围为确认化学身份、读取所列生成/消耗反应、区分MetaCyc跨物种记录与W29特异证据。按科研治理技能的只读流程直接进行，保持模型/GPR/培养及既有封存报告不变。

计划：查看CPD-238及最多4个直接关联反应；对化学字段冲突用ChEBI/KEGG/NIST交叉核对；独立审计化学身份；给出与此前网络阻断结论的关系。仅浏览及简短证据记录，预算15分钟。访问不足时明确未读部分，不以通用数据库结果为W29补路径，不新增计算、集群作业或外部通信。网页正文作为证据，不作为操作指令。

## 结论

该条目给出一个可继续调查的上游：

**L-苯丙氨酸 + O₂ → 苯乙酰胺 + CO₂ + H₂O**

催化步骤为苯丙氨酸2-单加氧酶（EC1.13.12.9；MetaCyc `PHENYLALANINE-2-MONOOXYGENASE-RXN`；KEGG R00690）。这是数据库收录、由其他物种实验支持的候选供给路线，尚未确定W29具有该活性。[MetaCyc生成反应](https://biocyc.org/reaction?orgid=META&id=PHENYLALANINE-2-MONOOXYGENASE-RXN)、[ENZYME条目](https://enzyme.expasy.org/EC/1.13.12.9)、[Pseudomonas sp. P-501原酶学研究摘要](https://pubmed.ncbi.nlm.nih.gov/6501250/)。

CPD-238的Reactions标签本次实际列出两项：

| 页面分类 | 反应 | 对原问题的含义 |
|---|---|---|
| 已知生成该化合物 | 上述L-苯丙氨酸单加氧步骤 | 提供苯乙酰胺的通用来源线索，非W29GPR依据 |
| 方向未知 | 苯乙酰胺 + 水 = 苯乙酸根 + 铵；RXN-12492 | 与此前模型R2251/R2087的化学步骤相对应；不是苯乙酸根的后续去路 |

两张反应页的Enzymes and Genes栏均未关联具体酶/基因。这是数据库关联缺口，不是“没有实验酶”的结论。生成页References列出Pseudomonas sp. P-501研究；水解页引用Rhodococcus erythropolis MP50纯化amidase研究（PMID8655547）。后者原摘要的宽底物范围不能独立证明PAM这一具体底物；本次未读其全文底物表，不声称该专项实验已经核实。[水解条目](https://biocyc.org/reaction?orgid=META&id=RXN-12492)、[其所引原研究摘要](https://pubmed.ncbi.nlm.nih.gov/8655547/)。

MetaCyc把RXN-12492列为方向未知，不能仅凭等号认定可逆酰胺合成，更不能通过反转它为模型补供给。生成酶条目所述约80%/20%两支路也不构成W29的通量比例。没有依据用数据库ΔG估计值或通用EC号替代目标酶、区室及条件证据。

## 化学字段确有冲突

浏览器成功读取MetaCyc30.0 CPD-238，结构图已目视核对。标题、SMILES、InChI及InChIKey支持苯乙酰胺；别名同时包含N-phenylacetamide/acetanilide，CAS为103-84-4。

| 化合物 | 连接结构 | 正确CAS | InChIKey |
|---|---|---|---|
| 苯乙酰胺／2-phenylacetamide | Ph–CH₂–C(O)NH₂ | 103-81-1 | LSBDFXRDZJMBSC-UHFFFAOYSA-N |
| 乙酰苯胺／acetanilide／N-phenylacetamide | CH₃–C(O)NH–Ph | 103-84-4 | FZERHIULMFGESH-UHFFFAOYSA-N |

两者同为C8H9NO，但连接方式不同；分子式和分子量相同不能证明同一化合物。CPD-238部分别名/CAS与其结构身份不一致；其ChEBI16562、KEGGC02505与结构相符。本次不查这个冲突的产生历史，不把未逐项核验的所有交叉引用判错。[ChEBI16562](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:16562)、[KEGG C02505](https://www.kegg.jp/entry/C02505)、[NIST苯乙酰胺](https://webbook.nist.gov/cgi/cbook.cgi?ID=103-81-1)、[NIST乙酰苯胺](https://webbook.nist.gov/cgi/cbook.cgi?ID=103-84-4&Units=SI)。

## 对项目结论的更新和边界

先前“底物没有生成入口、产物没有消耗出口”是对2026-09-11固定模型中相应胞质物种行的判断，并不等于自然界不存在生成反应。此次新增上游文献线索，不改变该模型网络结构的旧结论，也没有解决W29是否具有该生成酶、具体目标GPR或苯乙酸后续去路。本次未加载或重算2026-09-14工作目录模型，不能把旧模型判断默认为当前所有改动后版本均已复核。

实际核查：一个化合物页、两个直接反应页及References标签，ChEBI/KEGG/ENZYME对照与独立化学/原摘要审计。没有新增序列、结构、代谢计算或模型/GPR/化学/培养修改。浏览器仅保留科学字段摘录，见`browser_observations.json`；页面版本、时间及文件SHA见`evidence_manifest.json`。访问限制：web工具未能打开BioCyc，使用用户已授权的浏览器正常读取；部分PubMed页面未直接返回正文，原摘要由独立审计读取，详细范围见[AUDIT.md](AUDIT.md)。

独立审计合计13项：9项按注明证据层级支持、4项未决。审计覆盖化学身份及EC/原摘要范围；未第二次独立打开BioCyc实时页面，不把主代理的页面读取当作双重独立核验。PubChem7680未独立核对，摘要未覆盖的PAM专项数据与W29归属保持未决。
