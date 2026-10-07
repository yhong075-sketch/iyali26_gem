# ATP见证的方向、生物证据与管线来源审查

日期：2026-09-24。范围：R_PGAM1_PhosHydro、R694及追加见证中的R_NTP3pp、R_NTP7、r0242、R594、R603、R2010；静态读取当前工作区和公开官方/原始来源。未新增LP或改模型。本文件独立审核**方向理由与来源**；求解和见证最小性由主分析交付，不冒称本子任务复现。

## 直接结论

1. 首条见证应优先审查 **R_PGAM1_PhosHydro 的无能量耦联PEP合成方向**。若该列表示PEP phosphatase，合理的待验证模型候选是保留当前计量、限定 `v <= 0`（如[-1000,0]）；这保留PEP水解，阻止 pyruvate + Pi 免费生成PEP。R694正常的PEP向ADP转磷酸反应不应仅因出现在循环中而删除。
2. 追加见证同理，应优先审查 **R_NTP3pp的无能量耦联GTP合成方向**。官方对应为GTP水解；模型当前还存在单列H/电荷不平衡，方向和化学需分别验收。不能仅靠总循环中残差抵消宣称单列化学正确。
3. r0242所指化学是 **DHAP水解为DHA和Pi**，此次见证的逆向符合该水解方向；它的单列质子/电荷处理仍错误。R594核苷二磷酸激酶把GTP能量转给ADP符合酶类别，不应当作“循环修复”随意关闭。
4. 第三种见证的 **R_NTP7负向UTP合成** 有同类问题。它的存储正向本来是UTP水解，故合理方向候选为[0,1000]。本轮三个方向约束合起来只是隔离诊断候选；精确原生底物/GPR及单列H/charge残差仍未验收。

这给出有来源的候选方向；没有取得W29细胞内各反应的实际ΔG/浓度或精确原生底物实验，不宣称全部条件绝对不可逆，也不宣称限制上述一列即可清除全模型所有ATP循环。

## 来历：不是原始模型自带，也不是R1159修改造成

| 层级 | R_PGAM1_PhosHydro | R_NTP3pp | r0242 |
|---|---|---|---|
| `data/iyali26.xml` / `data/iyli21.xml` | 均无此列 | 均无此列 | 均无此列 |
| P0 gap-fill候选 | 由EC/数据库映射候选加入 | 同样为新增候选 | 同样为新增候选 |
| 早期方向整理表 | `reverse`、[0,1000]，存储PEP水解式 | `reverse`、[0,1000]，存储GTP水解式 | 本次表中无定向记录；后续notes称legacy_unreviewed |
| `metadata_reaction_selection.json` | before为水解[0,1000]；after反号为PEP合成且[-1000,1000] | before为水解[0,1000]；after反号、加入消耗H+且[-1000,1000] | before合成式产生2H+；after移除2H+；原[-1000,1000]未变 |
| 当前实际XML | 与after完全相同 | 与after完全相同 | 与after完全相同 |

追加的R_NTP7在两个原始XML也不存在。方向整理表原指定`keep`、[0,1000]、无额外H+的UTP水解；metadata选择保留主计量朝向，却新增生成H+并将边界改为[-1000,1000]，重新允许逆向UTP合成。

已直接读取清单指定的工作区内 `artifacts/coq9_pipeline_integration_20260909/accepted/metadata.xml`，完整SHA确为 `33c14f0c524197ec732e09ba86c4382d5b21bc25b4bd73b18eba5f452dd7e797`，上述四列字段与清单after及当前模型一致。这里“accepted”只是现有目录名，本轮不据此新增科研接纳。

代码流程也吻合：`main.py:311`调用gap-fill及方向整理表；`gaps.py`翻转计量、应用整理边界。较晚的`main.py:614`调用`apply_metadata_reaction_selection`，按明确before/after前态检查覆盖选定字段。`reaction_selection.py`把旧方向证据移到`metadata_previous_...`，并保留“用户选择版本，不是新化学/原生功能验证”说明。**当前notes里保存旧水解理由，不意味着当前执行边界仍是水解单向。** 这是现行版本选择规则产生的覆盖，不是执行函数把正负号误读。

`check_provenance.py`已实际执行并通过，UTP扩展后再次通过：原始两XML四列缺失、源metadata完整SHA、当前四列after字段、R_NTP3pp/R_NTP7/r0242的H/电荷残差均有assert检查。原始输入、代码和整理数据SHA及逐反应快照见`local_provenance.json`。本轮未重构历史运行环境。

## 化学步骤与蛋白身份分开

| 模型对象 | 系统ID、名称、功能与证据 | 已得到支持 / 仍未知 |
|---|---|---|
| R_PGAM1_PhosHydro、R_NTP3pp、R_NTP7、r0242共同GPR | **YALI1E41893g**（NCBI YALI1_E41893g；正式原生名未核实），A0A1D8NLA2/AOW06419.1、XP_504803.2；酸性磷酸酶候选，EC3.1.3.2自动注释、蛋白存在级4 Predicted | 支持泛磷酸酯水解候选；没有精确PEP/GTP/UTP/DHAP原生底物和胞质定位验证；不能从反应ID中的PGAM1推定它就是磷酸甘油酸变位酶 |
| R694 | **YALI1F12842g**（正式原生名未核实），A0A1D8NMM4/AOW06894.1；丙酮酸激酶候选，EC2.7.1.40自动/同源证据，存在级3 | 支持模型PEP+ADP→pyruvate+ATP的酶类别；原生功能未在本次实验确认 |
| R594 / R603 | **YALI1F12874g**（正式原生名未核实），A0A1D8NMQ6/AOW06896.1；核苷二磷酸激酶候选，EC2.7.4.6自动/同源证据，存在级3 | 允许核苷酸间转磷酸有类别依据；转移已有GTP或UTP能量不同于免费生成这些NTP |
| R2010 | 当前GPR为空 | 模型DHA磷酸化消耗ATP；不在本轮指定原生基因 |

完整本地UniProt候选条目（2026-01-28注释版本，非本轮实时更新）已抽出至`cached_candidate_proteins.json`。NCBI实时页面仍将YALI1_E41893g标为uncharacterized/provisional，W29来源和旧locus YALI2_F00484g明确，比较注释指向潜在酸性磷酸酶。[NCBI Gene2912348](https://www.ncbi.nlm.nih.gov/gene/2912348)

## 官方定义与原始酶学依据

**PEP支路。** [KEGG R00208](https://www.kegg.jp/entry/R00208)和[IUBMB EC3.1.3.60](https://iubmb.qmul.ac.uk/enzyme/EC3/1/3/60.html)对应 `PEP + H2O = pyruvate + Pi`。当前KEGG列EC3.1.3.60，模型赋的是泛酸性磷酸酶EC3.1.3.2；后者[IUBMB定义](https://iubmb.qmul.ac.uk/enzyme/EC3/1/3/2.html)为磷酸单酯水解，允许广底物不等于证明本蛋白的PEP活性。1989黑芥悬浮细胞纯化研究直接测得PEP转为pyruvate，PEP apparent Km约50 μM，宽pH最适约5.6，提出其可在缺磷条件下绕过需要ADP的pyruvate kinase步骤；该来源支持无ATP参与的**水解**，不是无耦联逆向PEP合成，也不是W29实测。[Duff等1989原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC1061789/)

PEP水解的负标准转化自由能还受原始热化学研究支持；本轮不将非W29条件数值填入W29模型。[Goldberg与Tewari，2003原始研究](https://doi.org/10.1016/j.jct.2003.08.002) 实际ΔG仍依赖浓度、pH及离子状态；不能仅把数据库等号或双箭头当作生理可逆证据。

**GTP支路。** [KEGG R00335](https://www.kegg.jp/entry/R00335)明确定义 `GTP + H2O = GDP + Pi`，列EC3.6.1.5及3.6.5.*，未列模型EC3.1.3.2。[IUBMB EC3.6.1.5](https://iubmb.qmul.ac.uk/enzyme/EC3/6/1/5.html)说明NTP到NDP是释放Pi的水解步骤；不存在该定义中的ATP/光/离子梯度耦联。因此当前 `GDP + Pi + H+ → GTP + H2O` 作为持续能量供给没有此来源支持。此为方向和EC映射疑点，不足以替模型任命一种新GTPase或换GPR。

**DHA支路。** [KEGG R01010](https://www.kegg.jp/entry/R01010)定义 glycerone phosphate（DHAP）+H2O→glycerone（DHA）+Pi，列EC3.1.3.1/3.1.3.2；见证中的r0242负向对应这一步水解。它是否由本候选蛋白实际催化仍未知。

**UTP支路。** [KEGG R00159](https://www.kegg.jp/entry/R00159)定义 `UTP + H2O = UDP + Pi`，列EC3.6.1.5及3.6.1.39，未列模型3.1.3.2。[IUBMB EC3.6.1.39](https://iubmb.qmul.ac.uk/enzyme/EC3/6/1/39.html)定义dTTP水解并明确该类别也能较慢水解UTP。与上面的apyrase定义一起，这支持水解类别方向；不支持把UDP+Pi无能量耦联地转成UTP，更不确认YALI1E41893g原生UTP活性。

**不能错误删除的转磷酸步骤。** [EC2.7.1.40](https://iubmb.qmul.ac.uk/enzyme/EC2/7/1/40.html)定义pyruvate kinase的PEP/ATP磷酸转移；[EC2.7.4.6](https://iubmb.qmul.ac.uk/enzyme/EC2/7/4/6.html)明确多种NTP可作供体、NDP作受体。故R694和R594作为能量转移步骤本身并未因参加计算见证被证伪。关闭它们只能算诊断阻断，不自动成为科学修复；没有依据删除GPR或认定相关基因“错误”。

## 按当前物种表示的化学问题

静态逐元素核算确认：R_NTP3pp当前正向H残差=-1、charge残差=-1；r0242当前正向H残差=-2、charge残差=-2。追加见证含2倍R_NTP3pp正向和1倍r0242逆向，故两者错误可以在总式里抵消。这种代数抵消不等于每个反应都质量/电荷守恒。

R_NTP7当前正向H/charge残差各+1，故第三种见证中它的2倍逆向同样各贡献-2，并被r0242逆向各+2抵消。只限制PGAM正向、NTP3pp正向和NTP7逆向不会自动修复这些化学残差；主分析的该隔离候选应继续标为**仅方向敏感性测试**。

后续隔离候选应分别处理：磷酸酶的合理水解方向、统一物种质子化表示下的逐列守恒、精确酶/底物/区室赋值。未授权正式应用前，不覆盖当前版本或旧证据。仅改方向可能阻断当前某条见证，但全模型能量验收仍须依靠主分析固定输入的复测。
