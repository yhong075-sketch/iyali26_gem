# R153 / R2176 独立来源与静态见证审计

审计人：独立子代理 `audit_r153_r2176`。核验时间：2026-09-11 22:20 UTC。审计问题与范围由 `PLAN.md` 固定；直接打开模型、可核实 Git 对象、原始映射表、数据库响应和已有通量，不以主代理结论或多数意见作为证据。新增优化求解、模型/GPR/标签变更、集群作业均为 0。本文件为审计交付，不是科学修改的接受记录。

## 身份与证据等级

| 系统 ID | 已核实名称 | 本次可支持的蛋白功能 | 模型角色 | 证据等级 |
|---|---|---|---|---|
| YALIUNK2 | 无已核实原生基因名；模型名称为 COBRAProtein175 | 未解析 | R153：精氨代琥珀酸合成 | model/GPR assignment only；模型占位条目，尚无可追溯原生基因/序列身份 |
| YALI1D17462g（格式变体 YALI1_D17462g；表列跨版本 ID YALI0D14300g） | 无已核实原生正式名 | A1 家族天冬氨酸蛋白酶候选 | R2176：精氨代琥珀酸合成 | 历史/当前未审阅数据库自动注释、同源性和预测结构；非原生功能实验验证 |

S2 中的酿酒酵母同源名称不能移植为目标基因的已核实名字。当前同序列记录使用的具体酶名称也不能升级为原生底物特异性或区室已确认。

## 原子声明审计

下表的“未决”表示证据不足或仅间接支持，不能读为已经证伪。支持项也仅支持其明确限定的层级。

| ID | 原子声明 / 待判问题 | 直接打开的来源与定位 | 类型与条件 | 判定 | 限制 / 处置 |
|---|---|---|---|---|---|
| C01 | 当前发布模型中两条反应的完整六物种胞质计量、界限相同，分别为单条目 GPR | `model_metadata_trna.xml`，R_R153 行 69327、R_R2176 行 158497；FBC bounds 参数 | 模型直接记录 | supported | 两条均 [-1000,1000]；不证明生物学同工酶 |
| C02 | 原始 `data/iyali26.xml` 与 `data/iyli21.xml` 内两条也重复，但完整计量比发布模型多产物 H⁺ | 两个 XML，R153 行 4819、R2176 行 33370 | 模型直接记录 | supported | 不得称旧/现完整计量不变；两原始文件不字节相同 |
| C03 | 2026-03-19 首个相关仓库导入已经含有此重复及现有两条 GPR | Git `05ae6e16955a180cdb6416143312ad806ada0440:data/iyli21.xml` 与该提交 `model.xml` | 本仓库历史直接记录 | supported | 更早初始模型为 e_coli_core、无目标；不能确定仓库外赋值起源或作者动机 |
| C04 | YALIUNK2 已被解析为真实原生基因或具体酶 | 原/现 geneProduct、映射 CSV、历史 UniProt 全蛋白组缓存 | 身份审查 | unverified | 没有交叉引用、序列或已核实名；仅有占位 ID。没有证据不等于不存在原生酶 |
| C05 | 缓存 `yaliunk2 → YALIUNK2` 不能补足原生身份 | `gene_locus_tag_map.json` lookup；`locus_resolver.py` 的 `build_lookup` 行 309–348 | 代码静态核实 | supported | 函数为每个模型条目建立本身的规范格式键；自映射不是跨版本/功能确认 |
| C06 | YALI1_D17462g 与 YALI0D14300g 的配对有显式表列证据 | S2 workbook sheet1 A4465:B4465；两张映射 CSV；历史 UniProt 同条目交叉引用 | 表列跨版本映射 | supported | 去下划线是格式归一化；跨版本配对依赖来源。映射本身不等于蛋白序列相同 |
| C07 | 目标蛋白的历史来源提出 A1 蛋白酶候选，构成对现有合成酶 GPR 的功能冲突 | S2 H4465 同源注释；缓存 A0A1D8NEI7 v28；UniSave v28 的 DE、CC、FT、GO、PROSITE | 自动/同源推断 | supported | S2 和数据库注释不是多个独立湿实验。唯一列出论文是大规模序列组装；本项不声称酶活已验证 |
| C08 | 现有材料足以确认目标是原生精氨代琥珀酸合成酶 | 同 C07；模型 R2176 GPR；IUBMB EC 6.3.4.5 | GPR 的生物学接受问题 | unverified | 未见目标特异合成酶证据，且有蛋白酶家族反证；应保留候选冲突，不能直接接受现有赋值 |
| C09 | 目标已有经核实的原生正式基因名 | 两个 UniProt 记录的 genes/ORFNames、S2 H4465 | 名称审查 | unverified | 仅系统 ID 与同源功能；不借用相关物种名称 |
| C10 | 本次检索到的六种序列表示在完整 384 aa 上完全一致 | 缓存 v28、UniSave v28、NCBI AOW04049.1、Q6C947 搜索响应、AF API、PDB CA 序列 | 字符串/序列直接核验 | supported | 完整序列 SHA 见下；此结论强于 ID 映射，但不证明功能或调控相同 |
| C11 | A0A1D8NEI7 当前处于 Inactive/DELETED，所给删除理由是退出参考蛋白组 | `sources/A0A1D8NEI7_current.json` inactiveReason | 当前数据库记录 | supported | 不把记录删除解释为序列虚假、基因不存在或功能被实验否定 |
| C12 | 使用的是与目标完全同序列的既有 AlphaFold 预测，信度与 PAE 数值可以复核 | AF API、模型 v6 PDB、PAE JSON | 预测结构；API 标注 Monomer v2.0、创建日期 2022-06-01 | supported | 发布文件 v6 与预测方法版本是不同字段；不是本次重新预测，且未进行参考结构比对 |
| C13 | AlphaFold 足以确认目标的原生催化功能、底物特异性及区室 | 同 C12；历史 SignalP/功能注释 | 功能/定位推断 | unverified | 高 pLDDT 不是实验。只能保留“基于 AlphaFold 预测的功能候选”；信号肽亦为预测 |
| C14 | 当前两条目的单敲各只关闭其关联反应，重复反应仍开放 | XML 全部 GPR 引用；`verify_static.py`；`KO_changed_bounds` | 静态逻辑及已有本次 COBRA 边界检查 | supported | 两条目各仅关联一条目标反应；单敲传导正确，不等于 GPR 生物学正确 |
| C15 | 旧 WT 与两个单敲代数见证均满足完整守恒、运行边界且保留旧 WT 生长值 | 旧 `results.json` / `reaction_snapshot.json`；本次 `static_verification.json` | 独立 XML/stdlib 复算；SD-Leu/PO1f 条件 | supported | 仅见证可行性；最优性相等依赖旧 WT 最优状态记录，本次无新求解或对偶证书 |
| C16 | 本模型精氨代琥珀酸池的完整物种行只有 R152、R153、R2176 | 发布模型全部 species/反应关联；`argininosuccinate_rows` | 计量直接核验 | supported | 得 v153+v2176−v152=0；双关闭会强制 v152=0，但尚不能推出全模型不生长 |
| C17 | 旧 screen 中两条目均预测非必需，且实验状态都是 unlabelled | 原 `screen_predictions.tsv`、`essentiality_per_gene.tsv` 的对应记录 | 历史结果核验 | supported | 不是本次重跑；不得转写成 TP/FN 或实验非必需 |
| C18 | 双关闭必致死或两条重复与单条 OR GPR 在所有条件严格等价 | 本次见证、两条 bounds、完整池行 | 未执行的干预/全域声明 | unverified | 无双关闭生长测试；两条总净容量可达 [-2000,2000]，同 bounds 单条 OR 为 [-1000,1000]，不能直接宣称优化问题等价 |

审计覆盖：**18 项 | 已审 18 | 支持 13 | 未决 5 | 已证伪 0 | 未审 0**。C04/C08/C09/C13/C18 是接受所需但尚未建立的声明，并非报告可作为事实使用的结论。

## 独立复算细节

审计人没有执行主代理验证脚本；使用标准库解析 XML 和 JSON，独立核对 2315 条反应的计量与旧有效快照一致，在全部 1877 物种上重新累加守恒残差。两条目标列均为：

`ATP_cy + L-Asp_cy + L-citrulline_cy ⇌ AMP_cy + PPi_cy + argininosuccinate_cy`。

旧 WT 的 R153=0、R2176=0.20917915810634488。对目标单敲的见证只做 R153←R153+R2176、R2176←0；YALIUNK2 单敲下旧 WT 原已满足对应零边界。三份通量的最大质量残差均为 **8.058630802595945e-15**，最大边界违反 **0**，目标均为 **1.8718823069403**。目标系数从本次有效记录读取；运行培养/菌株身份维持原封存配置，不回填历史默认值。13 个保护输入和 3 个载入源码 SHA 在核验时均匹配；只证明所列文件身份，未重建整个历史 dirty 环境。

PDB 的 384 个 CA 残基序列独立还原后与其他五份序列一致。PDB B-factor 字段 pLDDT 均值 **85.38703125**；AF API 的全局报告值为 **85.38**，两者分开保留。历史 Peptidase A1 61–366 残基域的 PDB pLDDT 均值 **91.61264705882353**；预测信号肽 1–17 均值 **36.20529411764706**；预测催化位点 D79/D262 的 pLDDT 分别 **95.00 / 95.38**。PAE 是 384×384 矩阵，全矩阵均值 **7.943413628472222 Å**，61–366 域内均值 **4.184309453628946 Å**；实际矩阵最大 **31 Å**，元数据 `max_predicted_aligned_error` 为 **31.75 Å**，二者不是同一个量。未进行结构比对，不能凭来源描述或信度判定已证明蛋白酶折叠/底物特异性。

## 来源与身份记录

- 当前目录 HEAD：`ff36d87eb6c8f933dfc4413f43a2fcefcf0eea27`。检查时已有多项未提交修改；没有切换分支或恢复它们。适用根 AGENTS/全局及祖先文件已核对；没有读取其他工作树。
- `model_metadata_trna.xml` SHA256：`d274bad3050e3c9220a8b6287eae847f3bf1334892284d565a6c4d96b38135a0`。
- `data/iyali26.xml` SHA256：`5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d`。
- `data/iyli21.xml` 及 `05ae6e16955a180cdb6416143312ad806ada0440` 中相应文件/当时 `model.xml` SHA256：`6974b7588f2a6c60ba2cde2f26e20d3aba1334d0d501572bf01cee47eda86631`。
- 更早 `0004d768d739f7519742ebfc08ecc290aed48df3:model.xml` 为 e_coli_core，SHA256：`671c33157ebf07874ea3c5b1928cfd011d127272f62bc6c830da8d58c9975ea9`。它没有目标条目，不构成更早的 Yarrowia 赋值记录。
- `data/yali1_yali0_map/S2_table_YALI1_YALI0_map.xlsx` SHA256：`42658fb95202b4114e5c6b5dbd60020de7a4038797ee2211a5352078a729a670`；本次直接读取 sheet1 A4465:H4465，未从衍生 CSV 推测原表。
- `iyli21_genes_vs_S2.csv` SHA256：`7c68b3cb244f07848b7ad99d7ed0f95768c9812ddd80311e1cdffa06fff4bae9`；`s2_metabolic_genes.csv` SHA256：`f7f25bb3e4f3c720d4b0bd534309c5c5f3cfa091bbd22383858d71f8600c1c2e`。
- 历史完整缓存 `artifacts/reference_pipeline_restore_20260909/research/cache/data/uniprot_UP000182444.json` SHA256：`e5b0a04874079b4057ffe25dadcb6b812cba8c96b227f41187187b27715753a3`。目标为 entry version 28、sequence version 1、384 aa、W29/CLIB89；其注释更新日期为 2026-01-28。
- 本次六方一致的完整氨基酸序列 SHA256：`bf9ffec20d0feb2f68a19fa10b102240dc138035b1997fa81b1c159fe0c192b3`。
- [UniSave v28](https://rest.uniprot.org/unisave/A0A1D8NEI7?format=txt&versions=28)、[NCBI AOW04049.1](https://www.ncbi.nlm.nih.gov/protein/AOW04049.1)、[当前 Q6C947](https://www.uniprot.org/uniprotkb/Q6C947/entry)、[AlphaFold API](https://alphafold.ebi.ac.uk/api/prediction/A0A1D8NEI7)、[IUBMB 6.3.4.5](https://iubmb.qmul.ac.uk/enzyme/EC6/3/4/5.html) 的下载响应均已打开；IUBMB 支持反应名称和化学步骤，不提供此目标基因功能证据。
- `retrieval.json` 当时 SHA256：`cfbf175ff3bae196ae76bd3d9cdc137c206cdaa4826bcb3762f0cffcd6a6f2e9`；其中 8 个下载文件的 SHA 全部匹配，保留各自 URL 与获取日期。本次分离缓存记录 `sources/A0A1D8NEI7_cached_v28.json` SHA256：`b42505278dc7ef2a621e53f5807e8069bf67140a584e844050ad94903ceaea7b`。
- 旧 `artifacts/iyli647_screen_20260910/nonessential_diagnosis_20260911/results.json` SHA256：`3e38ddf833e0f84c7bdeb21cccd3166c58f81b37ae09adab0f65c6291f163abb`；旧 `reaction_snapshot.json` SHA256：`d9eb0a95af7cd6944b37f0667623c565f2a0999dc82b303a2d4d54774ca4bde3`。
- 本次 `static_verification.json` SHA256：`0fd698e8162e910e56b04ded317f2294cc8d9b285d6d0cd53de0f5f98094d043`；`verify_static.py` SHA256：`227318e2980be8ec74f30b9d913fa9f69dfb3b7e6c13c2e97ed4f2e697cd4d85`，与记录相符。

尚未接受的结论：真实合成酶具体身份、目标原生功能和定位、反向容量的生物学合理性、双关闭后的生长、任何具体修模方案。当前证据足够确定重复导致的单敲旁路机制，以及目标 GPR 的身份/功能依据不足并存在冲突；不因此推荐具体新基因或批准模型修改。

## 最终报告复核

审计人随后逐段打开 `REPORT.md`，以及 `protein_evidence.json` / `verify_protein.py`。初次复核报告 SHA256 为 `af763ea53c7bc633741ba063f60f389c375faf47c8c4d4ff5b59f48a3afee11d`；蛋白证据 SHA256 为 `75b7b07aa0465aaa4fec6390633050aa979f82e6accbf8b76743b9cfade5db27`；蛋白核验脚本 SHA256 为 `9a810fe494e5dc8459c7dd119805319452fe32baa261f2b8a40f40260dd0d883`。蛋白证据内 9 个来源 SHA 均匹配，包含前述 8 个下载文件与 1 个历史缓存分离记录。补充独立核实 PDB pLDDT<70 的残基数为 65、最小值 30.92；Q6C947 的原始 organism 与文献菌株字段均明确 CLIB 122 / E 150。

报告准确区分了原生功能、家族自动注释与 AlphaFold 预测，未把较高信度说成折叠比对结果、催化验证或真实定位；数学见证、原始七物种/当前六物种、重复容量和双关闭限制表述得到支持。结论“高度疑似误赋、不宜作为已确认 GPR”是对现有连边接受程度的审查判断，不能理解为已经实验排除目标蛋白所有副活性。

为使措辞与证据等级一致，向主代理提出两项文字修订：将“应撤回其已确认状态”改为“不应将其视为已确认 GPR”，以免暗示本次已核实历史确认或审批；将“撤除错误连边”改为“撤除疑似误赋连边”。两项不改变证据、数字或审计覆盖。

已重新打开报告确认两项修订均已采纳，末尾覆盖数也与本表一致；最终报告 SHA256 为 `ef40fddf5da892a1ee8e9e26325ffb019b24a4f48a327395a7604ca139aa1c79`。没有待补正的关键数字或结论。另于本次审计通过官网打开 [PROSITE PS51767](https://prosite.expasy.org/PS51767)，核实其确为 Peptidase family A1 domain profile（条目名称 PEPTIDASE_A1）；官网定义支持家族名称，不是目标活性实验证据。

**本次没有完整序列同源检索、实验表征蛋白的系统比对、酶家族负对照搜索、结构比对或新的功能实验，因此不是最终 GPR 功能验证。** 此缺口限制新功能/GPR 的接受，不妨碍完成本轮现有连边的定向审查。审计支持当前模型机制解释和证据不足/冲突判断，保留全部五项未决声明；没有新增 GPR 指认或科学模型接受。
