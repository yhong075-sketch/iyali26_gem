# 酿酒酵母 CoQ 参考证据独立审计

审计日期：2026-10-01（America/Los_Angeles；网络日志 UTC 为 2026-10-02）。审计者：独立子代理 `audit_coq6_chemistry`。范围遵循上级 SCOPE.md；只写证据，没有修改模型、GPR、实验标签或科学代码，没有运行求解或预测。

## 可执行结论

**经过反应化学重建，COQ6 AND YAH1 AND ARH1 可作为有依据的、跨物种推断的 C5 耦合 GPR。** 这里 AND 表示催化酶与连续电子传递链共同参与 lumped 反应，不宣称稳定的三亚基复合体。若把铁氧还蛋白再生单独建模，则羟化反应催化 GPR 用 COQ6，供体再生由相应供体基因承担。旧 R39 的无供体、半氧气方程与该机制不匹配，不能只补三个基因而保留原方程。

**NADH 是合理的候选供电子辅因子；NADPH 也有支持，不能宣称其中之一具有原生排他性。** 酿酒酵母还原酶的原始体外实验接受两者，但采用异源受电子蛋白；这是 NADH 选择的具体生物化学根据，不能用“模型目前只有某个方程能平衡”替代这个根据。NADH 方案应标为待检验的模型假设，而非 Yarrowia C5 重构实验结果。

**COQ8 的证据落在途径/合成复合体层次，不能从这些论文挑出唯一催化步骤。** 将其加入某一个已有反应作为 AND 可以人为实施途径依赖代理，但应明确是建模约定；这些来源不支持宣称其是那个反应的独有必需催化组分。本任务宜保留辅助蛋白注释，不随意选择一个反应或全途径 AND。

## 基因身份与解释层级

| 原生候选 | 酿酒酵母参考 | 简要功能与等级 |
|---|---|---|
| YALI1A08781g；XP_499891.1 | YGR255C / COQ6；P53318，序列版本 1 | FAD 依赖 CoQ 环羟化酶候选；原生序列家族证据见此前独立审计，酿酒酵母 C5 遗传/代谢证据支持 |
| YALI1B03314g；XP_500417.1 | YPL252C / YAH1；Q12184 | 线粒体铁氧还蛋白候选；本次主代理完成序列对应，酿酒酵母条件性耗竭支持 CoQ 作用 |
| YALI1B19490g；XP_500902.3 | YDR376W / ARH1；P48360 | 线粒体铁氧还蛋白还原酶候选；酿酒酵母重组蛋白接受 NADH/NADPH；原生辅因子偏好未验证 |
| YALI1B20527g；XP_500941.1 | YGL119W / COQ8（ABC1）；P27697，序列版本 1 | ABC1 家族 ATPase/合成复合体辅助蛋白候选；无本次唯一反应归属 |
| YALI1F34675g；XP_505941.3 | YLR201C / COQ9；Q05779，序列版本 1 | 脂质结合/底物呈递辅助蛋白候选；酿酒酵母多步骤依赖证据，C6 比 C5 的绝对阻断证据更强 |

所有 Yarrowia 名称均为同源候选名称，不是本次原生酶学认证。本审计读取了主代理 donor_comparison 相关序列身份记录及其报告，未重复运行序列计算。既有 AF 模型属于 AlphaFold 预测，不能证明供体配对或辅因子排他性。

## 原始来源与反证

S1. [Pierrel 等，2010，HAL 原始稿](https://hal.science/hal-00497916/document)；DOI 10.1016/j.chembiol.2010.03.014。**全文已取回**：sources/Pierrel2010.pdf（HAL v1），Results PDF 页 7–8 / 稿件页 6–7，Figs. 2、5 说明。条件性耗竭、YAH1 回补和 Fe-S 配体突变支持 Yah1/Arh1 供体系统对 CoQ 的作用；细胞存活测定排除单纯死亡造成的 CoQ 下降。论文曾提出 Coq6 负责 C1、未知酶负责 C5 的假说（Discussion PDF 页 17–18），已被 S2 更具体实验修正；不能把这项早期推论当成当前催化注释。

S2. [Ozeir 等，2011，HAL 原始稿](https://hal.science/hal-00630764/document)；DOI 10.1016/j.chembiol.2011.07.008。**全文已取回**：sources/Ozeir2011.pdf（HAL v1，含公开审稿往来及补充材料；审稿往来不是独立实验）。定位：Results “Coq6 and Yah1 are functionally coupled…” PDF 页 17–19，Discussion 页 21–22，Fig. 4/S4 说明页 30–31、45。Yah1 耗竭不降低 Coq6 稳态量；FAD1 过表达不能救回；人 FDX2 可恢复部分 Fe-S/SDH 活性而不能恢复 CoQ，羟化前体可解除后一缺陷。这比单纯终产物减少更有力地定位到 C5。作者仍明确把 NAD(P)H→Arh1/Yah1→Coq6 称为机制假说，并要求纯化蛋白实验确认。异源 FDX2 不救援提醒不能建立宽泛供体 OR。其早期“仅 C5”结论也不能排除后来 2024 祖先重构蛋白的其他活性。

S3. [Lacour 等，1998](https://pubmed.ncbi.nlm.nih.gov/9727014/)；DOI 10.1074/jbc.273.37.23984。**原始论文摘要核查**：重组酿酒酵母 Arh1 在细胞色素 c 还原测定中使用 NADPH 与 NADH；表观 Km 分别 0.5、0.6 µM，受电子蛋白为牛 adrenodoxin。分馏定位于线粒体内膜。直接支持还原酶辅因子兼容性；不是酿酒酵母 Yah1/Coq6 或 Yarrowia 的 NADH 驱动 C5 重构测定。因而 NADPH-exclusive 论断有反证。

S4. [Xie 等，2012，PMC3390632](https://pmc.ncbi.nlm.nih.gov/articles/PMC3390632/)；DOI 10.1074/jbc.M112.360354。**全文已取回**，本目录 sources/PMC3390632.html 和 .txt。定位：Results “Overexpression of Coq8 Restores…”、Fig. 2–4；“Yeast Δcoq9 Cells…”、Fig. 8、Discussion（提及 supplemental Fig. S6，补图本身未独立取得）。Coq8 过表达可稳定多个蛋白并允许不同敲除株生成诊断性中间体；G130D 失去此作用，仍不能指定唯一底物反应。Δcoq9 + Coq8 过表达可积累 DMQ，故 C5 活性并非严格为零；绕过 C5 后 C6 缺陷仍存在。因此不能从该文推出所有 C5 都绝对需要 Coq9。

S5. [He 等，2014，PMC3959571](https://pmc.ncbi.nlm.nih.gov/articles/PMC3959571/)；DOI 10.1016/j.bbalip.2013.12.017。**全文已取回**，本目录 sources/PMC3959571.html 和 .txt。定位：Introduction 的 kinase caveat；Results 3.3、Figs. 4–5；3.4/Fig. 6。二维 native/SDS 电泳显示 Coq8 过表达在多个缺失背景影响不同大小的复合体。作者明确未直接证明 Coq8 蛋白激酶活性或底物。该文的中间体描述引用 S4，并非独立重复全部代谢测定。2012 Δcoq9 稳态量恢复很弱，与 2014 检出少量高分子量 Coq4 并不矛盾：二者读出不同。

S6. [Nicoll 等，2024，PMC7615680](https://pmc.ncbi.nlm.nih.gov/articles/PMC7615680/)；DOI 10.1038/s41929-023-01087-z。复用此前本代理全文审计：Results/Figs. 2–4 的祖先四足动物 COQ6 重构以 NADPH/FDXR/FDX2 实现羟化；酶单独加 NAD(P)H 无活性。该实验支持间接供电子机制，但不是酿酒酵母/原生 Yarrowia 供体专一性实验，也未排除 NADH 链。原始审计详见 `../../coq_gpr_completion_20261001/coq6_audit/AUDIT.md`。

S7. [Schulz 等，2022 online/2023 issue，PMC10873809](https://pmc.ncbi.nlm.nih.gov/articles/PMC10873809/)；DOI 10.1038/s41589-022-01159-4。复用此前本代理全文审计：Results “Both human ferredoxins are dispensable…”、Fig. 2e–f 中人细胞 FDX1 缺失、FDX2 耗竭及联合处理未显著下降 CoQ10 池。阻止将特定祖先体外系统或酵母依赖提升成全物种普遍规律；不直接反驳酵母结果。

S8. [Nishihara 等，2026](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0346295)；DOI 10.1371/journal.pone.0346295。全文网页核查，Results Fig. 7–8/Discussion。Pos5 缺失仍保留约 20% CoQ；研究以裂殖酵母为主，包含酿酒酵母表型。补前体和 Coq6 过表达结果没有将 NADPH 需求唯一定位到 C5，不能用于给 C5 强加 POS5 或排除 NADH。本审计不扩展该基因的原生候选搜索。

## 反应建议及停止条件

采用现有酸/阴离子状态的独立铁氧还蛋白形式：

`m641 (C52H78O3, 0) + O2 + 2 Fd_red + H+ -> m939 (C52H77O4, -1) + H2O + 2 Fd_ox`

Fd_red/Fd_ox 含相同原子，前者电荷少 1。若 lump 为还原型/氧化型辅因子，H+ 系数必须使用模型存储的实际分子式和电荷重新求平，不能只按辅因子名字抄标准方程。

主代理报告当前中性 NADH/NAD 的分子式相差 2 H；在该约定下候选为：

`m641 + O2 + NADH -> m939 + NAD + H2O + H+`

GPR：`YALI1A08781g and YALI1B03314g and YALI1B19490g`。三者依次为 COQ6 羟化酶、YAH1 铁氧还蛋白、ARH1 还原酶同源候选。分子式、电荷、区室和具体辅因子物种编号由主代理的独立实施检查确认；本审计没有重新读取全部辅因子。主代理另报告 NADPH/NADP 当前分子式与电荷不自洽；不要为此静默修改全局辅因子，也不要把选择 NADH 解释成排他生物学结论。

实施应保留旧规则/反应身份和差异记录；同时关闭无供体旧 R39 旁路，确认所有参与物种位于有依据的线粒体区室，检查代谢物生成/消耗及既有转运。R19 的 quinone/quinol 和氧化受体问题另行保留，不因 C5 候选通过而自动改变。不得把供体基因在既有 heme 反应中的共同出现当成 CoQ 原生实验证据。

## 审计覆盖

10 项原子主张均检查了来源及反证：6 项支持（其中跨物种 GPR 与 NADH 选择为明确的条件推断）、2 项被反证、2 项未解决；覆盖率 10/10。支持不意味着原生实验认证。

| ID | 主张 | 判定 |
|---|---|---|
| Y-01 | 酿酒酵母 COQ6 参与 C5 羟化 | supported，S2/S4 |
| Y-02 | 酿酒酵母 Yah1/Arh1 参与 CoQ 合成供电子链 | supported，S1/S2；直接生理电子传递仍为机制推断 |
| Y-03 | 酿酒酵母 Arh1 只能接受 NADPH | contradicted，S3 |
| Y-04 | Yarrowia NADH–Arh1–Yah1–Coq6 C5 完整链已经直接实验验证 | unresolved；所审来源未提供 |
| Y-05 | 经过化学重建的上述三基因 AND 可作跨物种候选 | supported as conditional model inference，S1–S3/S6 |
| Y-06 | NADH 可作为有生物化学根据的候选 lumped 供体 | supported as conditional inference，S3；不认证原生偏好 |
| Y-07 | 酿酒酵母 Coq8 影响多种 Coq 蛋白/复合体 | supported，S4/S5 |
| Y-08 | 所审 Coq8 证据能唯一选出当前某个催化反应 | unresolved；S4/S5 没有该分辨率 |
| Y-09 | 酿酒酵母所有 C5 通量均绝对依赖 Coq9 | contradicted under tested Coq8-OE condition，S4 |
| Y-10 | 酿酒酵母 Δcoq9 的 C6 缺陷比 C5 的完全阻断证据更强 | supported under tested conditions，S4 |

获取局限：S3 为原始论文摘要核查，不能把后续综述当作直接核查原始图。公开 HAL S1/S2 文稿经正常公开链接慢速下载成功，未绕过访问控制；取回初期批量程序中断，文件已完成，PDF 可解析且身份核实，最终来源日志据此补记。S1/S2/S4/S5 全文的相关 Results、Discussion 和图说明已检查；没有宣称重新分析原始实验数据。失败尝试与成功文件都保留。S4/S5 的部分中间体测定重复引用不能算新独立实验。
