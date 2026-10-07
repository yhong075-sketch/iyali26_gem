# DMeQ 前体还原：独立来源与计量审核

核验日期：2026-10-02（America/Los_Angeles）。本审核独立打开原始摘要、出版商主文及既有原始全文缓存，并独立静态读取指定 XML；未运行求解器、未改变模型或原始科研记录。本文替代先前的阶段一审核。

**最终报告的九项有界声明全部已审核：9/9 supported，0 partially_supported，0 unresolved，0 contradicted，0 unchecked。** 这里的 supported 包含有来源依据的限制性判断和条件性代数结论，不表示 W29 还原组分、供体、精确反应粒度或 GPR 已解决。这些生物学问题仍未确定。

被审核报告为上级任务的 `../REPORT.md`；审核快照 SHA256：`628ee62af14ba868f5e5309d4a14bae99ec7ef49b6fc54db70cad255da44da06`。完整输入、读取范围与失败入口另存 `MANIFEST.json`。

## 九项声明的互斥判定

| ID | 实际审核的有界声明 | 判定 | 独立依据和限定 |
|---|---|---|---|
| C1 | Houser 1977 三项观察支持大鼠肝线粒体 DMeQ9 前体还原活性；没有把该活性分子鉴定为 W29 的某个组分。 | supported | [PMID 863914 作者摘要](https://pubmed.ncbi.nlm.nih.gov/863914/)同时报告 NADH 生成可甲基化氢醌、Triton X-100 使前体还原活性失活而溶解甲基转移酶、连二亚硫酸盐在膜及溶解体系中部分替代 NADH。只读到摘要，没有声称获得全文或鉴定独立蛋白。 |
| C2 | Park 2024 的纯化 RTN4IP1／COQ3、DMeQ2、SAM、NADH 或 NADPH 耦联 Q2 终点与 DCPIP 试验支持候选机制；该耦联终点本身不是单酶生成游离 DMeQ2H2 的直接定量，也不是 W29／Q9 试验。 | supported | 独立读[出版商主文](https://www.nature.com/articles/s41589-023-01452-w) Fig.3、Fig.4c–e、相应 Results 及 Methods。DCPIP 是另一底物；PRM 测到的是 Q2 与 DMeQ2 信号。文中百分比按二者峰面积计算，最终报告未将其写为摩尔收率。未全面认证重组构建体版本及 COQ3 来源；本次未读完该文全部补充材料，不能把“所读主文未见”扩成全论文不存在某种实验。 |
| C3 | Oláhová 2025 的细胞回补恢复氧化态 Q10 和 PPHB10 异常而未恢复还原态 Q10，因此不能仅由救援现象唯一指定 DMeQ 还原步骤。 | supported | 独立读[出版商全文](https://link.springer.com/article/10.1038/s44318-025-00533-x) Fig.5D、相邻 Results 和 Discussion；WT、R103H、G215A 均有这一现象。作者保留残余催化或非氧化还原辅助解释；不能将这些变体在所有底物条件下一概当作完全失活，也不能据此直接否定 Park 的体外结果。 |
| C4 | Ouyang 2021 原始摘要报告 YMR152W／YIM1 的醛还原活性和三个受测醌的阴性；摘要列出的醌不是 DMeQ，不能据此推出 DMeQ 阳性或阴性。 | supported | 独立打开[PMID 32967812 作者摘要](https://pubmed.ncbi.nlm.nih.gov/32967812/)，并核对出版商索引摘要。六种醛有活性，三个醌为 9,10-phenanthrenequinone、1,2-naphthoquinone、p-benzoquinone。没有取得全文。**仅能确认摘要没有报告 DMeQ 测试，不能证明完整研究从未测试 DMeQ。** 最终报告已经采用这个有限范围，故判 supported。 |
| C5 | Lu／Nicoll 所核材料不能将游离 DMeQ 净还原归属给某个末端 COQ 蛋白；Nicoll 的耦联体系含葡萄糖脱氢酶／葡萄糖辅因子再生组分。 | supported | 独立读 [Lu 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3615049/) 既有全文缓存相关方法及 DMQ 电子转移结果、[Nicoll 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC7615680/) 全文缓存 Small-scale reactions／NAD(P)H regeneration system 与 Fig.5 相关结果。Lu 的 DMQ 底物不是 DMeQ。Nicoll 方法明确列出葡萄糖脱氢酶和葡萄糖，但这不证明再生酶直接还原 DMeQ。审计者本次只独立核主文缓存，补充材料的扩大核查由耦联来源报告单独记录，未计为本审核者独立读过。 |
| C6 | 本模型的中性辅因子存储约定下，NADH 供体假设式 DMeQ9 + NADH → DMeQ9H2 + NAD 的额外显式 H+ 系数为零。 | supported | 固定 XML 分子式／电荷精确求和，全部元素与电荷残差为零，见 `STOICHIOMETRY.json`。这是存储微物种下的候选记账，不确认 W29 的供体偏好或生理主导微物种。 |
| C7 | 同一约定下，R385 醌醇候选式 DMeQ9H2 + SAM → Q9H2 + SAH 的额外显式 H+ 系数为零。 | supported | 独立计算元素及电荷残差为零。[IUBMB EC 2.1.1.64](https://iubmb.qmul.ac.uk/enzyme/EC2/1/1/64.html)支持醌醇甲基化家族定义。上述两步任一侧机械加入一个 H+，将产生对应符号的 H 与电荷各一单位残差。 |
| C8 | 若新增内部 DMeQ9H2 只有 R385 一个系数为 −1 的连接，则每个可行稳态均满足 v385=0；这既不是本轮求解结果，也不说明全模型生长为零。 | supported | 守恒行直接为 −v385=0。前提是内部、非边界物种且没有其他生成、输入或连接；此证明不声称模型一定存在可行解。见 `CONNECTIVITY.json`。 |
| C9 | 现有成熟 Q9／Q9H2 的连接不能自动供给 DMeQ9H2；目前证据不将缺口自动落实为新增还原反应、GPR 或 R695 改写。 | supported | XML 中未有 C53H82O4 的物种；线粒体 m611 只接 R695／R385，成熟 m468／m471 的氧化还原连接不产生该新前体。`CONNECTIVITY.json`保存逐反应记录。这里没有将模型缺口反推为原生专用酶，也没有更改反应。 |

## 计量复核

模型 SHA256 与指定输入一致：`b83d7e7645c62cab5fc55c3b3e97bc7d6282dcdbf9eff4ad5aa2d2d24d38ce3e`。以下为线粒体存储记录，完整区室 ID 与所有 H+ 记录见 `STOICHIOMETRY.json`。

| ID | 存储物种 | 分子式 | 电荷 |
|---|---|---|---:|
| m30 | NADH | C21H29N7O14P2 | 0 |
| m27 | NAD | C21H27N7O14P2 | 0 |
| m60 | SAM | C15H22N6O5S | 0 |
| m62 | SAH | C14H20N6O5S | 0 |
| m611 | DMeQ9 | C53H80O4 | 0 |
| 未入模提案 | DMeQ9H2 | C53H82O4 | 0 |
| m471 | Q9H2 | C54H84O4 | 0 |
| m468 | Q9 | C54H82O4 | 0 |
| m28 | H+ | H | +1 |

NADH／NAD 存储式相差 H2，SAM／SAH 相差 CH2。按“产物减反应物”，两步零 H+ 候选各自 C、H、N、O、P、S、电荷残差都为零。任一步增加一个反应物 H+ 时，残差 H=−1、charge=−1；增加一个产物 H+ 时，H=+1、charge=+1。这里没有把中性存储约定等同为生理 pH 下的优势微物种。

## 报告措辞与未覆盖项

已核最终报告的关键结论与上述界限相符，**没有需要阻塞交付的科学措辞更正**。速记声明 C4 曾用“DMeQ 未测”；这一绝对写法必须收窄为“摘要未报告 DMeQ 测试”，最终报告已经如此处理。其“出版商摘要”来源说明保持原始撰稿来源，本审核另补 PubMed PMID 32967812 的直接访问证据，不改变生物学判断。

- Houser 与 Ouyang 的全文未取得。Park 和 Oláhová 本轮依据出版商主文及图注，未把图像或所有补充材料声称为已全面审查。
- [BioCyc RXN-11758](https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-11758)直接访问失败；页面反应、质子项、催化归属均未核实，不计入已支持声明，也不因访问失败判其错误。
- W29 的新还原位点、版本化候选序列、原生定位、精确 DMeQ9 底物能力、实际供体和 GPR 均未建立。本地注释名称未命中不证明缺少同源家族。根报告的限定说法正确；本审核没有重做全库同源检索。
- 人 Q8WWV3—RTN4IP1／OPA10 是有 NAD(P)H 氧化还原酶实验依据的线索；酿酒酵母 YMR152W—YIM1 的醛还原酶功能有底物实验依据，DMeQ 活性未建立。跨物种同源或定位不能自动变成 W29 的底物特异性／区室证据。
- COQ3 身份映射及既有 W29 COQ7 家族候选身份不是本轮重新审计对象。最终报告复用这些历史身份时明确证据层级，本审核未将其升级为本次复现。

九项判定针对本次选定声明，不是系统综述的文献召回率，也不是九项生物学机制均已验证。未进行 FBA、增长预测、新序列比对、结构预测、集群提交或模型修改。
