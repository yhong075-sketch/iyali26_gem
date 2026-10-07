# R_NTP1 化学方向：独立来源审计

核验日期：2026-09-14（America/Los_Angeles）。审计者：`/root/ntp1_energy_sources`。范围：R_NTP1 所涉及的 ATP 合成/水解方向、EC 7.1.2.2、EC 3.1.3.2、能量耦联及热力学适用条件。方法：直接打开官方数据库/规范与原始实验论文，主动检索逆向 ATP 合成反例；未运行 LP、序列/结构分析或修改模型。

## 审计边界

本审计独立核查以下外部来源。根任务转述的固定模型中性分子式、额外 H⁺、原子/电荷残差与全部物种位于 cy 的事实，**没有在本来源审计中独立读取 SBML 重算**，交由 `LOCAL_AUDIT.md` 核验。因而，下文对该模型的化学修正仅作带条件候选，不能把外部文献读取得到的证据等级扩展为本地模型验证。

## 五项声明判定

| ID | 被审计声明 | 判定 | 支持依据与限定 |
|---|---|---|---|
| C1 | 方向必须由实际化学过程、反应条件和完整耦联判断，不能按通量数值正负判断；不能把 ATPase 永不反转作为理由。 | supported | S4 给出条件依赖的反应 Gibbs 能判据；S5 直接观察外加机械功驱动 F₁ 合成 ATP；S6 测得同一重构 ATP synthase 随梯度及反应商改变合成/水解方向。通量符号相对于计量书写方向定义，是模型表示层的数学判断。 |
| C2 | KEGG R00086 的双箭头/等号不能单独证明目标酶在 W29 胞质可生理性净合成 ATP。 | supported | S1 只提供无区室的总式及多个 EC 关联；S2、S5、S6 表明具体 ATP 合成需要检查完整机制与驱动。此项是来源适用范围推论，不是 KEGG 明文宣布 R_NTP1 不可逆。 |
| C3 | EC 7.1.2.2 需要完整能量耦联；EC 3.1.3.2 分类本身不能证明目标蛋白催化 ATP 磷酸酐水解或 ADP+游离 Pi 的逆向净合成。 | supported | S2 的分类反应含膜两侧 H⁺，F 型合成机制由电化学梯度驱动；S3 的分类反应是磷酸单酯水解，虽注明转磷酸活性，却没有给出本目标蛋白的 ATP 底物实验或无供能合成证据。这里说“分类不足以证明”，并非断言该分类下所有蛋白均无相关旁活性。 |
| C4 | 未知 W29 的 pH、pMg、活度/浓度等时，不得给出该系统实际 ΔG；不能把离子化学式和固定 pH 的变换标准值混用。 | supported | S4 式 2-3 至 2-7 区分具体离子式与变换式，式 7-7 显示实际反应商项；固定 pH 时不能重复把同一化学 H⁺ 计入变换后的反应商。S6 的实验自由能值有其自身体系与条件，不能直接当作 W29 测值。 |
| C5 | 若 LOCAL_AUDIT 确认本模型采用所述中性物种且多出一个 H⁺，候选应先修正该化学计量，再审查水解方向；仅改 bounds 不能修复原子/电荷不平衡。 | partially_supported | “bounds 不改变计量，不能修复既有化学残差”是数学性质；S4 支持必须与实际物种表示一致的原则，S1 提供不含额外 H⁺ 的数据库总式参考。**本模型残差与去 H⁺ 后是否平衡须由 LOCAL_AUDIT 确认**。选择水解方向可作为无已知合成驱动情况下的审查候选，但特定酶/GPR身份和完整生理适用性仍需证据；本轮未实施或验收科学变更。 |

覆盖率（分母仅为上表五项声明）：

`total claims 5 | audited 5 | supported 4 | unresolved 1 | contradicted 0 | unchecked 0`

`unresolved 1` 对应 C5 的 `partially_supported`；不是外部来源互相矛盾，而是本地化学事实和修改验收在本审计范围之外。这里的 audited 表示按上述边界审核过该声明，不代表所有声明已获得完整支持。

## 六个主要来源及精确定位

| 来源 | 本次实际打开与定位 | 证据层级、条件和局限 |
|---|---|---|
| S1. [KEGG REACTION R00086](https://www.kegg.jp/entry/R00086) | 官方条目 `Name`、`Definition`、`Equation`、`Enzyme`。名称 ATP phosphohydrolase；总式 ATP + H₂O ⇌ ADP + orthophosphate，关联多个 hydrolase EC 与 7.1.2.2。 | 官方数据库证据。无目标蛋白、菌株、区室、实际活度、净生理方向或完整跨膜耦联的验证。数据库双箭头不足以授予特定模型自由反向边界。 |
| S2. [IUBMB EC 7.1.2.2](https://iubmb.qmul.ac.uk/enzyme/EC7/1/2/2.html) | `Reaction`、`Comments`：H⁺-transporting two-sector ATPase；反应显式包含膜两侧 H⁺，F 型复合体 ATP 合成旋转由 H⁺ 电化学势差驱动。 | 官方酶分类。可以支持“需要完整耦联机制”，不能支持把单独 ATP 化学子反应作为等价独立酶步骤。条目通用的 4 H⁺ 不是 W29 特异实测计量；未将其关于所有 A/V 型生理方向的概括用作本结论。 |
| S3. [IUBMB EC 3.1.3.2](https://iubmb.qmul.ac.uk/enzyme/EC3/1/3/2.html) | `Accepted name`、`Reaction`、`Comments`：acid phosphatase，磷酸单酯 + H₂O → alcohol + phosphate，注明宽底物谱与转磷酸活性。 | 官方酶分类。转磷酸活性不等于 ADP + 游离 Pi 无供能净合成 ATP；具体供磷底物、ATP磷酸酐底物活性和本目标蛋白身份仍待实验来源。未把“未证明”写成“证明不存在”。 |
| S4. Alberty et al., IUBMB Recommendations for Terminology and Databases for Biochemical Thermodynamics（2011） | [官方入口](https://iubmb.qmul.ac.uk/thermod2/) 与 [官方全文 PDF](https://iubmb.qmul.ac.uk/thermod2/BiochemThermodynamics.pdf) 均打开。PDF 自编号页 6–8，式 2-3 至 2-7；页 25，式 7-7。 | 官方规范全文。具体离子反应须平衡原子和电荷；固定 pH/pMg 的变换式使用物种总量并吸收对应化学势影响；实际方向依赖 ΔrG′° 与 Q′。规范不提供该 W29 条件的实测输入。PDF 页面截图未全部成功取得，方程定位以已打开全文文本为准。 |
| S5. Itoh et al., 2004, Nature 427:465–468, DOI [10.1038/nature02212](https://doi.org/10.1038/nature02212) | [出版社页面](https://www.nature.com/articles/nature02212) `Abstract` 直接打开：磁珠连到分离 F₁ 的 γ 亚基，电磁铁施加适当方向转动，荧光素酶检测 ATP 生成。 | 原始实验论文；本次直接核验限公开摘要和书目信息，**未核验付费正文/补充材料**。这是对“ATPase 绝不能反向”的反证，实验明确输入机械功，不能转用为无驱动胞质合成的依据。 |
| S6. Turina, Samoray & Gräber, 2003, EMBO J 22:418–426, DOI [10.1093/emboj/cdg073](https://doi.org/10.1093/emboj/cdg073) | [作者机构公开全文 PDF](https://www.biophysikalische-chemie.uni-freiburg.de/dokumente/publikationen/2003_turina_embo.pdf) 直接打开。印刷页 419 式 1–4；页 420–421、Fig. 2；页 422 关于酶活化与平衡的说明。 | 原始实验全文。叶绿体 ATP synthase 重构脂质体，改变 ΔpH 与 ATP/(ADP·Pi)，观察合成/水解转换；总能量包含化学项和跨膜电化学驱动项。条件、标准自由能和 H⁺/ATP 数值不能当作 W29 测值；未由这一相关体系确认本模型 GPR。 |

## 反证、失败访问及结论限度

- 反证已纳入：S5 的机械功驱动逆转以及 S6 的梯度依赖合成/水解转换，排除“ATPase 化学上永不反向”的绝对断言。
- IUPAC Gold Book 当前 catalyst 页面返回 403，旧页面返回 502；该定义不列入六个已直接打开的主要来源。另一个 Nature 蛋白旋转论文页面失败，也未作为结论证据。
- S6 的 PMC 入口遇到浏览器验证码，未绕过；改读作者机构公开 PDF。S5 仅公开摘要可读，未把摘要阅读标为全文审阅。
- 同一胞质内的化学 H⁺ 不是跨膜梯度项；如果本地读取进一步显示该 H⁺ 在当前中性式中属于多余项，还须单独修复计量。删除化学多余项不自动建立能量耦联，反之，限制方向也不自动修复化学不平衡。
- 外部来源支持拒绝“只凭双箭头、酶名或通量符号接受无耦联 ATP 合成”的论证。是否存在闭合产能循环、22 个保存解的可行性/频次，以及该固定模型具体残差，交由本地证据；此处未复算。候选水解方向尚未实施、未经过新的生长/必需性验收，也未被升级为已验证模型。

