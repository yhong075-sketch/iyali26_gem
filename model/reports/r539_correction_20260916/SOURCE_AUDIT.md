# R539 独立来源审核

审核日期：2026-09-16 UTC。审核者为独立子代理；目标是审查用户已授权的 R539 注释/GPR 修正，不预设单基因规则必然正确。已读取适用 AGENTS.md、`govern-agentic-research` 与 `gene-identity-function`。只读检查固定输入、原始数据库记录、IUBMB 定义和有限原始研究摘要；未修改模型、求解、提交集群作业或扩展全蛋白组搜索。

## 结论和可接受边界

**支持删除 R539 的 EC 1.3.1.104，保留 EC 2.3.1.39。支持将 YALI1E22262g 作为该步骤的推定催化基因，但单基因规则只能表示当前催化赋值，不能表示已经证实原生反应仅依赖这一个基因。** 当前七基因两分支规则没有取得对应于同一步骤、两套可互换催化组合的证据。

采用单基因催化规则时，必须明确：ACP 仍保留为反应底物；载体供给、载体基因敲除的传导、成熟修饰和复合体依赖尚未由该规则完整表达。尤其不能从两个候选都叫 ACP 就推断它们可以互相替代。不能据本轮静态审核宣称生长、必需性或基准命中率已改善。

## 固定输入与独立复核

输入为 `model_metadata_trna_r1931_forward.xml`，SHA256 `d417f1de0425bc3336503b45a9498154b8a25a5d90d6f58e68b9d80cfc79297b`。审核实际读取 XML：R539 的四个物种均在 C_mi，计量为 ACP + malonyl-CoA ⇌ CoA + malonyl-ACP，边界沿用输入。全模型按精确有理数计量、保留物种及区室身份检查，未找到同计量或反号计量的另一条反应；这不排除多步替代路径。

旧规则为 `(YALI1F38317g and YALI1A20089g and YALI1C26939g and YALI1D18037g) or (YALI1D32594g and YALI1F37498g and YALI1E22262g)`。其他基因存在时，删除最后一个基因仍使第一分支为真；这是布尔规则结论，不是通量或生长结论。

| 系统 ID | 名称、蛋白功能及证据等级 | 对本次判断的意义 |
|---|---|---|
| YALI1E22262g | 原生正式名称未核实；推定 malonyl-CoA:ACP 转酰酶；自动数据库注释、同源序列、既有 AlphaFold 预测辅助支持，非原生酶活验证 | 与 R539 化学直接匹配的催化候选 |
| YALI1D18037g | 原生正式名称未核实；推定线粒体酰基载体蛋白；W29 RefSeq XP_502836.1 与 UniProt A0A1H6PXT9 自动注释 | 属于载体候选，不能当作替代转酰酶；具体依赖未定 |
| YALI1D32594g | 原生正式名称未核实；推定线粒体酰基载体蛋白；W29 RefSeq XP_068138955.1 与 UniProt A0A1D8NG21 自动注释 | 同上；不因同一功能类别而建立 OR |
| YALI1F38317g | 原生正式名称未核实；β-ketoacyl-ACP 缩合酶候选，本轮限于模型 R1394–R1396 赋值 | 此模型功能不证明转酰酶活性 |
| YALI1A20089g | 原生正式名称未核实；3-hydroxyacyl-ACP 脱水酶候选，本轮限于模型 R1400–R1402 赋值 | 不同化学步骤；不能由通路参与推导 R539 必需亚基 |
| YALI1F37498g | 原生正式名称未核实；3-oxoacyl-ACP 还原酶候选，本轮限于模型 R1397–R1399 赋值 | 同上 |
| YALI1C26939g | 模型映射至 YALI0C19624g / ETR1；推定线粒体 enoyl-ACP 还原酶；UniProt Q6CBE4 reviewed 注释的功能/定位依据仍是同源转移，本轮未核实跨版本序列相同 | EC 1.3.1.104 属于该还原化学，不能因此加入 R539 |

## 身份、比对和 AlphaFold 的审核范围

独立从保存的 RefSeq XP_504110.3、GenBank AOW05616.1 解析序列，与缓存 UniProt A0A1D8NJ03（entry 33、sequence 1）、AlphaFold API 及 PDB 的 306 个 CA 残基比对，均相同。纯序列 SHA256 为 `46da679cce71c1098390e273f0cdaa26b27f8d37d5477dd9275d9c879ca87a55`。现行 UniProt 将该条目停用的理由是“不属于 reference proteome”，不是功能否定。现行 KEGG YALI2 序列在第150/280位与 W29 不同，不能冒称同一序列；本轮未替换输入。

已独立重算保存比对的 identity、覆盖和位点对应，并核对 BLAST 原始 XML 指纹：E. coli `b1092 / fabD`（Malonyl-CoA:ACP 转酰酶；P0AAI9 已整理且有实验参考记录）在四参考小面板内为 33.11% identity、96.73% query coverage、E=3.84056e-32；S. cerevisiae `YOR221C / MCT1`（推定线粒体同类转酰酶，已整理注释及基因破坏证据）只有一个 HSP，28.08% identity、41.18% query coverage、E=7.96879e-8。参考催化 Ser/His 映射至目标 S97/H210；这是保守性证据，未实测目标残基功能。

`YPL148C / PPT2`（线粒体 ACP 磷酸泛酰巯基乙胺转移酶；已整理/实验参考记录）及 `YML022W / APT1`（腺嘌呤磷酸核糖转移酶；已整理/实验参考记录）在该面板和阈值下无命中。**APT1 不是酰基蛋白硫酯酶，这两个对照不能代表相近旁系同源蛋白的功能排除。** 未执行全蛋白组或完整数据库检索；该 E-value 仅适用于四参考面板，不能声称唯一催化基因。

复用的 **AlphaFold 预测**为 AF-A0A1D8NJ03-F1 v6，API 记录 AlphaFold Monomer v2.0 pipeline、创建日期2022-06-01；API globalMetricValue 91.44，独立从 PDB CA 重算平均 pLDDT 91.40879，306×306 PAE 平均4.90876 Å、最大30 Å。两种平均值差异如实保留。未新增预测，未进行实验结构叠合；高置信度折叠不能证实底物特异性、定位、伙伴需求或原生酶活。接受的功能表述为“基于 AlphaFold 预测并结合序列/注释支持的功能候选”。

## 来源和反证

1. [IUBMB EC 2.3.1.39](https://iubmb.qmul.ac.uk/enzyme/EC2/3/1/39.html)：转移 malonyl 至 ACP，与模型化学相符。[IUBMB EC 1.3.1.104](https://iubmb.qmul.ac.uk/enzyme/EC1/3/1/104.html)：涉及 acyl/enoyl-ACP 与 NADP(H)，与 R539 四物种化学不符。官方定义本次已直接打开并另存 HTML。
2. [RefSeq XP_504110.3](https://www.ncbi.nlm.nih.gov/protein/XP_504110.3)、[GenBank AOW05616.1](https://www.ncbi.nlm.nih.gov/protein/AOW05616.1)：W29 306 aa 的版本化身份，名称仍为 uncharacterized/hypothetical；RefSeq 注释含 FabD domain 和对 MCT1 的序列相似推断。数据库注释并非原生蛋白功能实验。
3. [UniProt A0A1H6PXT9](https://rest.uniprot.org/uniprotkb/A0A1H6PXT9.json)、[UniProt A0A1D8NG21](https://rest.uniprot.org/uniprotkb/A0A1D8NG21.json)：两载体候选的 ACP 域、pantetheine 修饰 Ser 和线粒体定位均为自动注释。两条 W29 RefSeq 原始记录也已直接读取。不能从这些注释决定特定载体 AND 或 OR。
4. Schneider et al.，1997，[Two genes of the putative mitochondrial fatty acid synthase](https://pubmed.ncbi.nlm.nih.gov/9388293/)：本次直接读取 Europe PMC 的原文摘要记录，研究是 S. cerevisiae 同源序列与基因破坏呼吸表型；不将其当作 W29 目标的纯化酶学证明。UniProt 的实验依据标签不能消除原研究物种和实验终点的范围限制。
5. Dobrynin et al.，2010，[Characterization of two different acyl carrier proteins in complex I from Yarrowia lipolytica](https://doi.org/10.1016/j.bbabio.2009.09.007)：直接读取原研究摘要，ACPM1 与 ACPM2 的删除表型不同且不能互补，是“载体无需表达/任意互换”解释的反证。未取得摘要中 ACPM 名称到本轮两个 W29 位点的逐序列闭环，故不据此给任何一个位点直接安上已实验验证的 ACPM1/2 身份。原文详细样本量、株系及方法未获全文核查。
6. [Acyl modification and binding of mitochondrial ACP to multiprotein complexes](https://doi.org/10.1016/j.bbamcr.2017.08.006)，2017：直接读取原研究摘要，报告 ACPM1 除 complex I 结合组分外，也存在于游离基质及 Fe-S 组装相关复合体。它更新了早期“只检出膜上组分”的解释。此证据属于载体，不可挪作目标转酰酶的定位实验。
7. [AlphaFold AF-A0A1D8NJ03-F1](https://alphafold.ebi.ac.uk/entry/A0A1D8NJ03)：本次读取 API、PDB 和 PAE 原始文件，核查同序列和置信度；未以结构单独判定功能。

## 逐项判定

| ID | 原子声明 | 判定 | 实际依据和限制 |
|---|---|---|---|
| C01 | R539 化学对应 EC 2.3.1.39 | supported | XML 与 IUBMB；不判断方向热力学 |
| C02 | EC 1.3.1.104 与 R539 当前计量不匹配 | supported | 缺少该酶定义的还原底物及 NADP(H) |
| C03 | 旧规则在目标单敲后仍为真 | supported | 两 AND 分支的布尔逻辑；其他基因存在 |
| C04 | 目标 W29、缓存和 AF 为同一306aa序列 | supported | 本次独立逐序列比较 |
| C05 | 目标原生酶活和底物特异性已直接验证 | unverified | 本轮未取得原生直接酶学证据 |
| C06 | 目标是该化学的有依据功能候选 | partially_supported | 注释、近全长 FabD 比对、催化位点保守；仍有间接性 |
| C07 | 复用 AF 身份和保存置信度统计正确 | supported | API/PDB/PAE 本次独立解析；不升级为功能证明 |
| C08 | 已排除全部替代催化基因 | unverified | 仅四参考面板，未完成近旁系同源排查 |
| C09 | 两个 D 位点具有 ACP 载体自动注释 | supported | 直接读取 UniProt/RefSeq；不等于特定反应伙伴已验证 |
| C10 | 旧两组 AND 已被证明是 R539 可互换催化组合 | unsupported | 当前模型赋值没有对应复杂依赖证据；未称生物学绝不可能 |
| C11 | 可将目标作为当前单基因催化赋值 | partially_supported | 适用于明确限定的催化 GPR 表达，不表示完整载体依赖 |
| C12 | 目标的 W29 线粒体定位已直接验证 | unverified | 仅自动 GO/同源推断，未见原生定位实验 |
| C13 | 固定模型没有同/反号计量重复 R539 | supported | 本次精确静态检查；不排除多步旁路 |

覆盖：**total claims 13 | audited 13 | supported 7 | unresolved 6（partial 2、unsupported 1、unverified 3）| contradicted 0 | unchecked 0**。这里 audited 表示已检查可用来源及其缺口，不表示每一项都得到支持。

## 可追溯性与失败记录

`audit/` 保存独立下载的来源、完整文件 SHA256 和检索时刻；父任务 `sources/`、`references/` 保存共同审核的目标和参考原始记录。`audit/audit_checks.json` 保存独立身份、结构统计、比对数值和重复反应检查。

Python 首次外网读取因本机证书链失败，改用系统 curl 并保持 TLS 证书验证后获取公开记录，原失败单独保留。NCBI Gene HTML 是浏览器检查页，不用作身份/功能证据；改用版本化 RefSeq 原始记录。浏览器读取部分 PubMed/ScienceDirect 页面遇到检查页/403，原研究摘要改由 Europe PMC API 读取；2014 年关联论文全文 XML 取得失败，未作为已核实全文证据。首次静态重复反应检查用整数解析计量，在0.5处停止；随后改用精确 Fraction 完成，未改科学输入或求解。

此报告审核截至当前证据，未审核后续生成模型、整理数据补丁、最终构建或回归结果；实施验证必须另行记录。R78 中类似旧规则只作范围提醒，未给予修改授权解释或生物学判定。

## 整理记录科学措辞审阅

于 2026-09-16T16:39:13.291477+00:00 读取 `data/reference_build/curation/r539_gpr_assignment.json`，SHA256 `9bb47595121fd89f8ffd0bcc261642c3588e651de0b162fba03a494a8dd12bd4`。其中仅 R539 改为上述 EC 和单候选催化规则，保留 provisional、待实验确认、定位/载体基因依赖限制及同序列 AlphaFold 来源；科学表述符合本审核接受范围。可继续已授权构建。此处不确认实现、最终输出或软件测试通过。
