# CoQ9 生物合成 GPR 研究

2026-09-18。依据用户截图，范围为 CoQ1–Coq9 生物合成系统；呼吸链 CI/CIII 的亚基重构不混入本轮。实际对象是项目状态入口指定的 `model_metadata_trna_r539_alphafold_labeled.xml`，不是根目录旧 `model.xml`。本轮读取当前模型、核查版本化蛋白、比对一手来源及既有 AlphaFold 预测，未修改模型、管线或整理数据，未运行生长求解。

**结论：当前已有六个候选基因关联七条催化反应；Coq6/8/9 的基因条目存在，但没有反应 GPR。现有缺口必须区分化学表示与辅助依赖，不能通过恢复全途径 AND 一次解决。**

## 1. 当前关联及证据

下表的 COQ 名称是家族/候选名称；数据库自动注释、跨物种实验和本菌株直接酶学证据不等价。除特别注明的 Coq1 异源构建外，本轮没有取得这些 Yarrowia 位点的直接催化/定位实验。所有位点均已用版本化 W29 蛋白记录核对，完整序列身份、源文件指纹及比对见 [identity_comparison.json](identity_comparison.json)。

| 系统 ID | 家族候选及蛋白功能 | 当前模型 GPR | 证据等级与决定 |
|---|---|---|---|
| YALI1C26017g | COQ1；聚异戊二烯二磷酸侧链合成酶 | R763 单基因 | 注释/同源支持，并有使用 Yarrowia 催化核心的异源补偿实验；保留候选，不声称原生定位已确认。 |
| YALI1F08349g | COQ2；4-HB 聚异戊二烯基转移酶 | R407 单基因 | 自动注释与同源支持；保留。 |
| YALI1B20835g | COQ3；CoQ 环两次 O-甲基转移 | R715、R385 各为同一单基因 | 自动注释、同源及跨物种步骤实验支持；保留两处关联。 |
| YALI1F34625g | COQ4；C1 脱羧/合成复合体作用候选 | R40 单基因 | 家族支持；原生纯脱羧还是耦联氧化脱羧未定，现关联仍为暂定表示。 |
| YALI1C25352g | COQ5；环 C-甲基转移酶 | R18 单基因 | 自动注释、同源及跨物种步骤实验支持；保留。 |
| YALI1A08781g | COQ6；FAD 依赖羟化酶候选 | 无；R39、R19 均无 GPR | 家族支持，但当前反应不充分表示酶化学；先解决电子供体、区室和 C1 分支。 |
| YALI1E18269g | COQ7；去甲氧基泛醌羟化酶 | R695 单基因 | 同源、酶分类及基于 AlphaFold 预测的功能候选支持；原生辅助依赖未定。 |
| YALI1B20527g | COQ8；ABC1 家族 ATPase/辅助蛋白候选 | 无 | 比较序列支持；合成复合体辅助作用未表达，不能直接扩成全途径 AND。 |
| YALI1F34675g | COQ9；脂质结合/底物呈递辅助蛋白候选 | 无 | 结构域、同源及基于 AlphaFold 预测的功能候选支持；原生 COQ7 配对及绝对依赖未确认。 |

原生正式名称未独立确认的位点保持候选称谓。表中“保留”是研究建议，本轮没有执行 GPR 变更。

Coq1 的 [Saeed 等，2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11237637/)使用 XP_501989.1 指定的 Yarrowia 蛋白，并替换 N 端靶向段后在酿酒酵母得到 CoQ9/功能补偿。因此支持催化核心及链长功能，不能作为原生线粒体定位实验。本文未逐碱基重建历史构建；文中片段长度与现记录全长的对应也未闭合。

## 2. 截图需要修正的步骤对应

1. **Coq1 制造侧链，Coq2 接上芳环。** FPP 不能与 4-HB 直接组成完整 CoQ9 前体；需要先延长侧链。[IUBMB EC 2.5.1.39](https://iubmb.qmul.ac.uk/enzyme/EC2/5/1/39.html)
2. 图中 `3-Hexaprenyl` 是六个异戊二烯单元，与标注的九单元侧链和 CoQ9 不一致；本项目应使用 nonaprenyl 对应名称。
3. **Coq3 出现两次**：早期 R715 与末端 R385；不能只画一次，也不能把 Coq4/5/8/9 全放在 R385 之后。[EC 2.1.1.114](https://iubmb.qmul.ac.uk/enzyme/EC2/1/1/114.html)、[EC 2.1.1.64](https://iubmb.qmul.ac.uk/enzyme/EC2/1/1/64.html)
4. 当前模型的末段是 **R18/Coq5 → R695/Coq7 → R385/Coq3 → Q9**。Coq8/9 应显示为辅助作用，不能画成无依据的独立化学转化。[EC 2.1.1.201](https://iubmb.qmul.ac.uk/enzyme/EC2/1/1/201.html)、[EC 1.14.13.253](https://iubmb.qmul.ac.uk/enzyme/EC1/14/13/253.html)

实际模型 R763 是 **C25 pentaprenyl-PP + 4 IPP → C45 nonaprenyl-PP + 4 PPi**，并非一条直接从 FPP 开始的反应。C25 与 IPP 通过 R913/R865 供给线粒体；4-HB 邻接 R52/R978。这里仅核查邻接，不升级为这些前体运输/定位的生物学确认。

模型当前主链为：

`R763 → R407 → R969（逆向，转出线粒体）→ R39（胞质）→ R808（转回线粒体）→ R715 → R40 → R19 → R18 → R695 → R385`

这是模型拓扑的读数，不是已确证的原生反应顺序。11 条主链/运输反应在当前存储物种下元素和电荷残差均为零，但“配平”不等于“机制正确”。R18/R385 保留 oxidized quinone 表示，与官方 quinol 酶反应之间的氧化态限制仍然存在。

## 3. 三类尚不能直接补入的 GPR

**Coq6：优先处理 R39/R969/R808。** R39 在胞质，以 0.5 O₂ 完成形式羟化，没有还原电子供体。完整 C5 单加氧酶反应需要还原力；官方分类写出 reduced ferredoxin。[EC 1.14.15.45](https://iubmb.qmul.ac.uk/enzyme/EC1/14/15/45.html) 这不是往旧 R39 加上 YALI1A08781g 就能解决的缺口。原生 ferredoxin/reductase 位点本轮没有鉴定，不能虚构其 GPR；即使未来确认供电子链，也应区分独立电子传递反应与同一酶复合体的 AND。

**Coq4/6：R40/R19 保留两种机制。** [Nicoll 等，2024](https://www.nature.com/articles/s41929-023-01087-z)的祖先蛋白短链体外体系支持分步脱羧/羟化；[Pelosi 等，2024](https://iris.cnr.it/retrieve/12ea9ce7-85c7-491f-80e7-8c0b7e131d03/COQ4.pdf)的细胞/异源体系支持 COQ4 耦联氧化脱羧。两者不能直接裁决 Yarrowia 的原生顺序。R19 当前还合并产物氧化，空 GPR 表示未知，不表示自发。

**Coq8/9：辅助作用不能自动变成硬 AND。** Coq9 与 Coq7 的脂质处理协作有人蛋白结构支持，[7SSS/Manicki 等，2022](https://www.rcsb.org/structure/7SSS)。Nicoll 的特定体外条件中，无 Coq9 仍测得 Coq7 活性，加入后效率上升；这不能反向证明原生 Coq9 可有可无。Coq8 的 [2026 新研究](https://pubmed.ncbi.nlm.nih.gov/42525751/)支持 ATPase 驱动的脂质中间体伴侣机制，本轮只核实摘要；2024 体系也曾报告蛋白磷酸化。因此保持 ATPase/kinase-like 的分层记录，不指定任意 ATP/CoQ 耦联系数。

当前七条催化反应的 GPR 均为单基因，没有 OR；同一个 Coq3 出现在两条反应中不应改成别的 AND/OR。未完成全蛋白组旁系搜索，不能宣称每个反应不存在其他同工酶。

## 4. 本次新增的可复核证据

- 九个位点均与版本化 W29/GenBank 蛋白及本地缓存全序列一致；八个用 RefSeq，Coq3 用 AOW01767.1，后者不自动验证模型里的另一条 XP_500950.3 交叉引用。当前九个旧 UniProt W29 条目均显示因不属于 reference proteome 而停用，不是功能被反证；缓存日期和停用状态分别保存。
- 用现有 BLASTP 2.17.0+ 将九个候选对照九个 reviewed 酿酒酵母 Coq 蛋白及 ERG20/BTS1 两个侧链酶对照，共 11 条参考。最佳命中均为预期家族；Coq1–8 一致性约 42–58%，Coq9 约 24%，其 E-value 9.46×10⁻²¹、目标覆盖约 97%。这是小参考面板比较，不是全库/全蛋白组搜索；酿酒酵母 CAT5 即 Coq7 参考。
- 复用九个**已有 AlphaFold 预测**，都与目标序列一致，保存版本、pLDDT/PAE；元数据中的算法为 AlphaFold Monomer v2.0 pipeline，v6 指 AFDB 文件发布版本，不是算法第六代。未提交新预测或 HPCC 作业。Coq3/8 全序列平均 pLDDT 约 68/69，不能把整条预测都当高置信结构。
- 因 Coq9 序列一致性较低，进一步以现有 USalign 比较两个预测与实验 7SSS 的两种蛋白链。Coq9 对 human Coq9 的 TM-score 为 0.644（按完整目标长度）/0.902（按实验链长度），178 位点对齐、RMSD 2.08 Å；Coq7 同家族结果为 0.828/0.914。交换家族对照的分数明显较低。结论限于折叠相容的**基于 AlphaFold 预测的功能候选**，未验证原生结合界面、活性位点、区室或必需性。[结构比较](coq79_structure_comparison.json)
- 实际调用静态 gene KO：Coq1/2/4/5/7 分别关闭对应单反应，Coq3 关闭 R715/R385，Coq6/8/9 不关闭任何反应。逐次恢复全部边界和基因状态；无生长求解。11 条目标未找到同区室完全相同或反号计量的重复反应，不代表没有其他多步旁路。

上述酿酒酵母对照为 YJL167W — ERG20 — 法呢基二磷酸合成酶；YPL069C — BTS1 — 香叶基香叶基二磷酸合成酶；YOR125C — CAT5/COQ7 — 去甲氧基泛醌羟化酶（均使用 reviewed 数据库注释及所列原始文献，不把参考物种证据移作 Yarrowia 实验证据）。

## 5. GPR 完整性与生长必需性必须分开

本次重新计算当前模型的 Q9 与 Q9H₂ 两行之和：全模型仅 R385 列剩下 +1，其余全部为零。因此对这份模型的通常稳态约束：

`d(Q9 + Q9H₂)/dt = v(R385) = 0`

这说明当前静态 SBML 没有 CoQ 池随生长稀释/周转的净需求。呼吸链的 Q9↔Q9H₂ 循环不消耗 CoQ 总量。仅把 Coq6/8/9 塞进 GPR，不能据此建立可信的生长必需性测试。本结论不覆盖另行添加需求或动态池约束的旧运行时模型，也没有替代 WT/KO 生长计算。

优先级建议：先明确 Coq6 的线粒体电子供给与反应重构，再裁决 Coq4/6 的 C1 分支；Coq8/9 保持辅助依赖待决。净 CoQ 需求是另一个表示问题，本轮不任意补系数或重开 dFBA。

## 6. 管线来源、范围与审核

当前单步 GPR 已在真实构建入口中调用 `apply_reviewed_quinone_step_gprs`（`scripts/gem_annotate/main.py`），规则在 `scripts/gem_annotate/patches.py`；原来的七基因混合 AND 已被拆成单步候选。本次确认其结果仍存在于实际 XML，未重新构建或修改这些规则。

机械核查采用 [extract_current.py](extract_current.py)，原始反应/物种/GPR/边界见 [current_model_extract.json](current_model_extract.json)，逐基因证据表见 [coq_gene_gpr_review.tsv](coq_gene_gpr_review.tsv)。来源身份和运行范围见 [scope.json](scope.json)、[identity_comparison.json](identity_comparison.json)、[execution_notes.json](execution_notes.json)。首次静态检查遇到 COBRA 无关联基因的状态恢复问题；显式恢复该软件标记后检查通过，模型文件未变，无优化调用。

独立机制来源审计覆盖 6/6：5 项支持、1 项部分支持，原生 Coq4/6 机制及 Coq8/9 布尔依赖仍未闭合；详见 [mechanism_source_audit.json](mechanism_source_audit.json)。本地序列/模型事实独立核验 7/7 项支持，见 [local_identity_audit.json](local_identity_audit.json)。结构数值和解释的独立审核见 [structure_audit.json](structure_audit.json)。已读取文献不等于本次复现实验，结构预测不等于功能确认。
