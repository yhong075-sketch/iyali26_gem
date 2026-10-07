# CoQ6 / CoQ4 independent source audit

审查日期：2026-10-01（America/Los_Angeles）。审阅者：独立子代理 `/root/audit_coq6_chemistry`。按本轮 SCOPE，只读核查既有身份、当前构建产物与一手来源；仅本目录写入证据。无模型修改、优化、序列比对、结构预测或集群执行。

**决定：YALI1A08781g 可以作为重构后 C5 单加氧酶步骤的暂定催化 GPR；不能直接用它认证当前 R39 或 R19 的完整化学表示。R40 已有单基因规则可保留为顺序脱羧假说，未取得足以替换 C1 路线或指定原生电子供体的证据。**

## 固定身份及证据等级

| 系统 ID / 蛋白版本 | 名称、蛋白功能与等级 | 本轮裁决 |
|---|---|---|
| YALI1A08781g / XP_499891.1，465 aa | COQ6 家族候选；FAD 依赖 CoQ 环单加氧酶；自动注释、跨物种同源及基于 AlphaFold 预测的功能候选。RefSeq 仍写 uncharacterized | 家族相容性有依据；原生底物活性及定位未实验确认 |
| YALI1F34625g / XP_505938.1，262 aa | COQ4 家族候选；C1 脱羧/合成复合体相关蛋白；自动注释、同源及 AlphaFold 预测支持 | 当前 R40 属暂定顺序脱羧表示 |
| YALI1B03314g ↔ YALI0B02222g / XP_500417.1，161 aa | YAH1/铁氧还蛋白候选；名称为同源候选，模型已有赋值 | 供体线索；未认证为原生 COQ6 伙伴 |
| YALI1B19490g ↔ YALI0B14839g / XP_500902.3，464 aa | ARH1/铁氧还蛋白还原酶候选；名称为同源候选，模型已有赋值 | 还原酶线索；未认证原生辅因子偏好或 COQ6 配对 |
| 酿酒酵母 YGR255C / P53318 | COQ6，CoQ 环羟化酶；该物种遗传/代谢物证据 | 参照，不移作 Yarrowia 实验 |
| 酿酒酵母 YPL252C / Q12184；YDR376W | YAH1，铁氧还蛋白；ARH1，铁氧还蛋白还原酶；该物种条件性耗竭证据 | 两种连续电子供给角色，不能理解为两个替代催化酶 |
| 人蛋白 NP_004100.1；NP_001026904.2 | FDX1、FDX2，线粒体铁氧还蛋白；人细胞实验参照 | 用于检验跨物种供体必需性的过强推论 |

COQ6 的既有 BLASTP 2.17.0+ 小参考面板结果为 P53318：41.722% identity，query 36–461/465，reference 29–471/479，E=8.6e−102。未覆盖的 N 端不能用核心同源性替代定位证据。复用 AF-Q6H9P2-F1 的同序列 **AlphaFold 预测**（2026-09-18 取得，AFDB 文件 v6、AlphaFold Monomer v2.0 pipeline）：既有记录 mean pLDDT 84.14、mean PAE 9.28 Å；本次核查序列及文件指纹，没有重新预测或结构叠合。原生定位实验及本审阅者实际运行的定位预测均无。主任务后续定位预测若取得，应独立增加其版本、原始输出和限制，不回填成原生实验。

两名供体候选来自本地 S2 映射和 NCBI W29 feature table。Schulz 的原论文仅在 Extended Data Fig. 8 序列比较中列出 XP_500417.1，并非该 Yarrowia 蛋白的 CoQ 实验；该文将当时记录称为 CLIB122，而本地新 feature table 将同版本 accession 关联 W29，本次没有做两株原始序列完整性核验。不得将这条元数据桥接升级为全序列相同。供体完整序列/结构与功能审查尚未完成。

## 一手证据与冲突

以下来源均于上述审查日实际打开；只有标为摘要的资料未取得全文。

- [Ozeir 2011](https://pubmed.ncbi.nlm.nih.gov/21944752/)，DOI 10.1016/j.chembiol.2011.07.008，**摘要**：酿酒酵母突变与前体绕过支持 C5 作用；对 YAH1/ARH1 的生理供电子解释仍为推断。
- [Pierrel 2010](https://pubmed.ncbi.nlm.nih.gov/20534343/)，DOI 10.1016/j.chembiol.2010.03.014，**摘要**：两种条件性耗竭均影响酿酒酵母 CoQ 合成；不能证明 OR 替代催化。
- [Nicoll 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC7615680/)，DOI 10.1038/s41929-023-01087-z，**全文 Results / Fig. 2–4、Methods、Data availability**：重构四足动物祖先蛋白和单异戊二烯底物。COQ6 单独加 NAD(P)H 无可检出活性；FDXR/FDX2/NADPH 体系支持 C5 羟化。COQ4 产生脱羧中间体；COQ6 再生成 C1 羟化的还原态产物。两种结果限于该体外体系。构建 accession：tAncCOQ4_tr OQ859711；tAncCOQ6_tr OQ859713；tAncFDXR_tr OQ859718；tAncFDX2_tr OQ859719。未将祖先蛋白当成人体原生蛋白。
- [Pelosi 2024 作者稿](https://iris.cnr.it/retrieve/12ea9ce7-85c7-491f-80e7-8c0b7e131d03/COQ4.pdf)，DOI 10.1016/j.molcel.2024.01.003，**Results / Fig. 3–4，PDF 页 9–12**：人/酿酒酵母 COQ4 可补偿细菌双缺陷，异源细胞中形成含氧 C1 产物，支持耦联氧化脱羧。未补偿已脱羧底物的单 hydroxylase 缺陷；因此不能给现有 R19 单独填 COQ4。与 Nicoll 的条件不同，尚不能据此裁决 Yarrowia 的机制。
- [Schulz 2022 online / 2023 issue](https://pmc.ncbi.nlm.nih.gov/articles/PMC10873809/)，DOI 10.1038/s41589-022-01159-4，**Results 的 ferredoxin/ubiquinone 小节及 Fig. 2e–f**：HEK293 的 FDX1 knockout、FDX2 RNAi 及联合干预没有显著降低所测 CoQ10 含量（各条件 n≥3 或 n≥4）。这是对“祖先体外供体必然是所有原生体系必需供体”的反证；并不否定祖先体外活性。RNAi、残余蛋白、细胞条件及稳态池量与新生通量的差别限制该结论。**Extended Data Fig. 8** 仅提供 Yarrowia 序列参照。
- [IUBMB EC 1.14.15.45](https://iubmb.qmul.ac.uk/enzyme/EC1/14/15/45.html) 与 [1.14.15.46](https://iubmb.qmul.ac.uk/enzyme/EC1/14/15/46.html)，**官方反应定义**：采用完整 O2、两个单电子还原当量、水；后者产物是 quinol。属于酶学命名来源，所引主要实验为 Nicoll，同源引用不能重复算成独立实验。

另访问 Ozeir 2015 的 PMC / PubMed 页遭到浏览验证，未将旧报告中其 C4 结果计作本轮已打开的一手核验。未用其支持任何新增实施决定。

## 最小化学提案：尚未写入模型

本轮 `baseline.xml` 中 R39 仍为胞质中性酸 + 0.5 O2 → 阴离子羟化物 + H+，不含还原供体。这里“元素配平”与“完整单加氧酶化学”不同。

对于固定存储的质子化状态，可提出线粒体催化半体系：

`m641[C_mi] + O2[C_mi] + 2 Fd_red[C_mi] + H+[C_mi] → m939[C_mi] + H2O[C_mi] + 2 Fd_ox[C_mi]`

其中 m641 为 C52H78O3、charge 0，m939 为 C52H77O4、charge −1；Fd_red 和 Fd_ox 原子相同，电荷分别 z−1 和 z。这一式中 H+ 系数是 **1**，源于本模型中性酸→阴离子的约定，不可直接照抄 IUBMB 的 2。按此符号定义，C/H/O 及电荷均守恒。

候选催化 GPR：`YALI1A08781g`，标为跨物种机制支持的功能假说。Fd 的还原再生是另一步，保持原生供体身份、还原酶及 NADH/NADPH 偏好未决；不得加无代价 source/sink 来使其可行，也不能把 NADPH 的直接反应写成已经验证的 COQ6 本征活性。若采用 NADPH 合并式，必须明确隐藏了供电子链及其未核实依赖。

这一提案可使用已有线粒体底物/产物，不再要求 R969→R39→R808 的往返路径。迁移 R39、关闭原胞质旁路及增加 redox species/再生步骤是额外化学/区室变更，须共同决定、共同检验；单独新增线粒体副本会保留不消耗还原力的旁路。本次没有实现这些改动。

C1 保留两种可区分方案：

1. **顺序方案**：R40 的 `m111 + H+ → m63 + CO2` 及暂定 `YALI1F34625g` 可保留；随后 COQ6 候选催化 `m63 + O2 + 2 Fd_red + 2 H+ → C52H80O3(quinol) + H2O + 2 Fd_ox`。与现有 R19 的 C52H78O3(quinone) 相差两个氢；下游氧化必须另有闭合电子受体/机制，不能把两步含糊合并后宣称已获精确酶学支持。
2. **耦联 COQ4 方案**：R40+R19 的形式总式为 `m111 + H+ + O2 → m59 + CO2 + H2O`。这是对当前公式的代数求和，不是原生已验证反应；Pelosi 的主要异源底物还没有 C5 甲氧基，不能将该论文直接视为当前 m111 的精确底物实验。原生路径顺序及产物氧化态仍须裁决。若日后实施此方案，需要联合替换 R40/R19，不能仅给旧 R19 填 COQ4。

## 原子命题裁决及实施边界

| ID | 审核的命题 | 裁决 | 实施决定 |
|---|---|---|---|
| C6-01 | 固定原生 COQ6 候选具有 COQ6 家族相容性 | supported | 允许候选身份注释；不是原生活性认证 |
| C6-02 | 原生 COQ6 线粒体定位已被实验确认 | unsupported | 不接受此声明；只限制声明与迁移置信度 |
| C6-03 | Nicoll 体系直接支持耦联供电子的 C5 催化 | supported | 支持上述化学候选，保持体系限制 |
| C6-04 | 旧 R39 已完整表示该单加氧酶化学 | contradicted | 阻止直接认证旧式；须重构 |
| C6-05 | 两个 Yarrowia 供体线索已认证为原生 CoQ 伙伴 | unsupported | 保留候选；暂不指定 donor GPR |
| C6-06 | YAH1 OR ARH1 可表示两个可替代的完整羟化催化者 | unsupported | 不复制该 OR；联合链不等于任意复合体 AND |
| C6-07 | 祖先 FDX 体外结果证明所有原生真核 CoQ 必需该供体 | contradicted | 排除跨物种必需性推广 |
| C6-08 | R40 单基因规则可保留为暂定顺序脱羧假说 | supported | 保留已有候选；不宣称唯一原生路线 |
| C6-09 | Nicoll 直接观察了 COQ6 对短链 C1 底物的羟化 | supported | 仅候选；不能跨越 redox/供体差异 |
| C6-10 | Pelosi 可直接证明 COQ4 催化旧 R19 的已脱羧底物 | contradicted | 不给旧 R19 填 COQ4 |
| C6-11 | 旧 R19 与已测 COQ6 C1 单加氧酶反应完全相同 | contradicted | 阻止直接认证；产物氧化态及供体不符 |
| C6-12 | 原生 Yarrowia C1 路线已由现有两篇研究唯一确定 | unsupported | 保留两方案，不强行填齐 |

审核覆盖：**12 总数｜12 已审｜4 支持｜4 未解决｜4 反证｜0 未核查**。表中部分命题刻意检验过强断言，反证数不代表四项已接受结论发生错误。历史表达检测、AlphaFold 预测及他物种实验均不算原生催化验证。本报告未验证生长表型，也不接受正式科学模型。
