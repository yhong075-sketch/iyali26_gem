# 两基因非必需解释的独立来源审核

审核日期：2026-09-15（America/Los_Angeles；交付核对为2026-09-16 UTC）。审核者为独立子代理。本次读取源模型、保存运行记录及原始文献，未新增优化、FVA、模型/GPR/标签修改或集群作业；仅写入本审核文档。

## 审核结论

报告准确区分了“保存结果经本次核实”“局部计量候选旁路”和“尚未确认的实际KO通量机制”。核心分类及计量有直接支持，不应把PPP候选写成两个单KO已证实的补偿路线。PFK删除表型的相反实验报道已保留，并说明条件差异。

| 对象 | 身份、功能与证据范围 | 已核结果 |
|---|---|---|
| YALI1F11049g — PGI1对应 | 葡萄糖-6-磷酸异构酶；本地映射至YALI0F07711g，物种数据库提供推定同源注释；主代理另检索到异源活性原文，但本审核者未成功独立打开该文全文 | KO/WT=0.864386923123922；R326关闭 |
| YALI1D20222g — PFK1对应 | ATP依赖6-磷酸果糖激酶；本地映射至YALI0D16357g，物种内酶学、删除及回补有实验依据；R637额外底物赋值并未在本轮实验核实 | KO/WT=0.9705218581245482；R636与R637均关闭 |

## 声明审计表

| 编号 | 声明及裁定 | 实际打开的来源与定位 | 限制 |
|---|---|---|---|
| A01 | 输入和保存结果身份一致：supported | `evidence.json.identities`所列7个本地文件均重新计算SHA匹配；screen manifest模型项和输出SHA；报告引用的模型SHA匹配 | 字节身份不等于重建完整历史dirty环境；本轮未重跑求解 |
| A02 | PGI对应条目按本次图谱规则非必需：supported | `screen_retry/raw_deletions.tsv:876`、`screen_predictions.tsv:876`、manifest WT；原始增长1.6180305877462309，WT 1.8718823069402954，严格未舍入比值<0.15 | 只确认保存结果及分类，不独立再生成完整可行解 |
| A03 | PFK对应条目按同一规则非必需：supported | 上述两TSV各第487行；原始增长1.8167026947221614；两行均optimal且数值审核有效 | 两者均远离阈值，10%改15%不改变分类 |
| A04 | 目标GPR确实传导到关闭，无完全相同的重复计量列：supported | 原XML的R326、R636、R637均为单基因规则；原始关闭集合吻合；独立逐列检查完全相同的计量字典，未发现另一个相同列；runner第299–312行检查传导与恢复 | 检查未涵盖缩放/反向等价、多步等价或完整生物学同工酶审查 |
| A05 | 图中的E为保留的实验正例参考，而非模型P：supported | 9月15日`annotation_status.json`两条均experimental_positive=true、model_essential=false；`annotate.py`保留旧实验字段；历史来源表Sheet1行559/289 | 本轮未重读原始三套screen calls；不声称重新核实共识构成或独立验证 |
| A06 | 本地新旧ID对应关系可追溯：supported | `data/yali1_yali0_map/iyli21_genes_vs_S2.csv:853`及`:469`；[GRYC旧ID注释](https://gryc.inrae.fr/db/yarrowia-lipolytica/clib-122/yali0f/yali0f07711g)；[2005年PFK研究](https://www.microbiologyresearch.org/content/journal/micro/10.1099/mic.0.27856-0)原刊摘要；2023文§3.1明确PFK旧ID | 映射为现存项目记录，未重新进行序列跨版本匹配；2005文章全文实验细节未取得 |
| A07 | 两种PPP组合的局部净计量和允许方向正确：supported | 独立从源XML按物种ID求和，并与`evidence.json.candidate_stoichiometry`比较；源边界及运行manifest基础边界允许所需方向；全部物种为胞质 | 只是局部代数，6 NADPH及其他物种仍需完整网络平衡 |
| A08 | 保存WT不能直接作为PFK整个基因KO见证：supported | `artifacts/r153_r2176_merge_20260915/run_manifest.json`的`results.WT.fluxes`：R636=0，R637=2.2168967849989154 | R636单反应零通量不等于整个基因可删；未把WT当成KO结果 |
| A09 | 2021年的PFK删除株有葡萄糖最小培养基生长缺陷：supported | [2021原始研究](https://link.springer.com/article/10.1186/s13068-021-01962-6)，Construction and characterization…小节及Fig.4图注：YB-392/NS1047，YNB+2%葡萄糖；YPD仍有较慢生长 | NADPH过剩是作者提出的机制假说，不是本模型失配原因已获证实 |
| A10 | 2023年研究报告不同条件下PFK删除生长影响很小：supported | [2023原始研究](https://www.frontiersin.org/journals/bioengineering-and-biotechnology/articles/10.3389/fbioe.2023.1098116/full)，Table1、§2.5、§3.1–3.2及Fig.2图注 | 见下方培养与终点差异；不能据此验证本screen的97.05% |
| A11 | PPP就是保存两KO高生长的实际支撑路径：unverified | 原screen未保存完整KO通量；报告明确保留此缺口 | 未证明全网络可行完成、实际通量贡献、唯一性或生物学成立 |
| A12 | PGI异源表达活性的原文及Table1名称对应：unchecked（独立审核范围） | 主代理检索到[糖酵解移植研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC12802204/)相关段落；审核者尝试打开PMC及ScienceDirect全文未成功 | 不将主代理转述冒充独立打开；可作为主代理查得证据保留，但该具体实验声明未完成双重来源审核 |

覆盖：`total claims 12 | audited 11 | supported 10 | unresolved 1 | contradicted 0 | unchecked 1`。这里的0项contradicted表示表内经过限定的声明未被反驳；并不隐去文献间不同条件下的相反表型。

## 局部代数独立检查

`3R325 + 3R71 + 3R639 + R712 + 2R714 + R764 + R765 + R766`得到：

`3 G6P + 6 NADP⁺ + 3 H₂O → 2 F6P + GAP + 3 CO₂ + 6 NADPH + 12 H⁺`。

再加`−2R326`得到：

`G6P + 6 NADP⁺ + 3 H₂O → GAP + 3 CO₂ + 6 NADPH + 8 H⁺`。

独立求和与交付证据完全一致。质子沿用源模型计量，不在此更改化学。第一组合绕开R326，第二组合绕开R636/R637，存在净产物，因此不是净零循环；这不自动构成可支持生物量生产的完整稳态通量。

## 2023年反向证据的适用范围

原文§2.5的发酵配方为葡萄糖40 g/L、硫酸铵1.1 g/L、YNB 1.7 g/L、CSM-Leu 0.74 g/L，**另加L-亮氨酸0.20 g/L**，C/N=80。条件为30°C、220 rpm、120 h，每24 h测OD600等。Fig.2标题使用CSM-leu字样，但不能因此把实际配方解释成缺亮氨酸的当前SD-Leu条件。Table1记载的工程菌株背景也不同于本screen。该文支持“该背景和培养条件下观察到生长影响小”，不证明所有条件下不必需，也不是最大生长率的直接比较。

## 对交付文件的核阅范围

已通读`REPORT.md`，并读取`evidence.json`的来源身份、原始目标行、培养/求解上下文、候选净计量、限制和文献项。独立计算核验了7个来源SHA、目标反应完全相同计量列检查及两种PPP净式。未逐条独立复查94条局部反应的全部注释或145个物种的生物学身份，也未重新检查原始实验逐screen calls。报告没有把这些未审核范围写成已全面验证。

结论可用于解释此次保存预测与E参考的分歧；未授权或接受任何科学模型修改。PGI异源活性的独立全文核验及实际KO路径证实仍属未决，前者不影响保存分类与局部计量结论。
