# CoQ9 反应与基因校订：实际核对与候选交付

日期：2026-09-09。对象限定为用户指定的四组，不进行 dFBA、参数扫描或生长预测。

## 结论与交付范围

已从用户的封存交付包取得登记参考模型本体，并直接解析反应、代谢物、GPR及边界；不是仅依据状态报告重述问题。

- 输入：`coq9_wp2_qualified_baseline_20260905/frozen_inputs/model.xml`。
- SHA256：`bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee`。
- 实测2313反应、1877物种、1074基因；19条目标及邻接反应完成精确元素/电荷核算。
- 取得的是与项目登记一致的封存参考，不是本次远程读取的最新工作树。未修改用户生产仓库。
- 生成了数学模型不变的元数据/布尔去重候选，以及只改变R305两个质子系数的独立科学候选。二者都不是“整个CoQ9模块已验证”的正式发布版。

源代码、完整计算与候选见本目录。下文[Sxx]对应`source_inventory.tsv`；模型事实可复算于`reaction_audit.json/tsv`。

## 1. 四组对象的最终判定

| 对象 | 本次确定结果 | 处置 |
|---|---|---|
| R385及上游 | R385计量平衡，已有COQ3候选关联；不是整条合成途径都缺基因。11条当前主链在存储物种约定下全部元素/电荷平衡，但R39与R19仍不是充分的原生单加氧酶表示。 | 保留R385当前声明的形式表示及GPR，修正旧EC文字；优先解决上游电子供体和区室，不重做储备试验。 |
| R305 | 名称错误；精确残差H=-2、电荷=-2。A14736经本地跨版本映射可连接到原生CIII结构里的A14806蛋白；不是只有模型内猜测。 | 已修名称/反应EC；输出完整Q-cycle的2/4质子候选；将“亚基”和“细胞色素c载体依赖”分开处理。 |
| R1889/R2062 | R1889有泵但无GPR；R2062无泵且47个AND出现项只有28个不同基因，含NDH2。R2062与R573计量完全相同，但R573已关闭而R2062仍开放。 | AND重复项已等价去重；重建单一泵型CI的正确GPR并处理活跃无泵重复表示，不把旧AND直接复制过去。 |
| COQ6/8/9 | COQ6是催化步骤候选，不能直接绑到缺电子供体/区室有问题的R39；COQ8/9是需要单独表示的辅助功能，不是给所有反应加AND。 | 三者角色分别归类，明确缺口；保留基因，不把“无关联/KO无效”误判成非必需。 |

## 2. R385及上游：保留什么，修正什么

### 2.1 已有基因关联并非全部未定

本地跨版本表对应如下；这是已封存映射，不是本次重新完成的基因组比对。

| 家族候选/功能 | YALI1 | 本地YALI0对应 | 当前直接关联 |
|---|---|---|---|
| COQ1 | YALI1C26017g | YALI0C18755g | R763 |
| COQ2 | YALI1F08349g | YALI0F05610g | R407 |
| COQ3 | YALI1B20835g | YALI0B15884g | R715、R385 |
| COQ4 | YALI1F34625g | YALI0F27247g | R40 |
| COQ5 | YALI1C25352g | YALI0C18205g | R18 |
| COQ6 | YALI1A08781g | YALI0A09042g | 无 |
| COQ7 | YALI1E18269g | YALI0E15224g | R695 |
| COQ8 | YALI1B20527g | YALI0B15664g | 无 |
| COQ9 | YALI1F34675g | YALI0F27313g | 无 |

保留已有COQ3/5/7家族对应，不把“没有原生酶学实测”自动等同于“不能做任何注释”。但是同源/家族支持、原生复合物成员、原生单基因必需性是不同等级，不能合并成一个“正确”标签。

### 2.2 R385

当前模型：

`SAM + 3-demethylubiquinone-9 -> SAH + ubiquinone-9`

GPR为`YALI1B20835g`。计量在存储的中性物种形式下完全平衡。COQ3承担两次O-甲基化的功能归属有官方酶分类及实验来源支持[S01]。因此没有理由仅因R385是池的唯一净来源而删除、拆掉或换基因。

需要修正的具体点：当前`PROTEIN_CLASS:2.1.1.114`属于第一次O-甲基化；末步应为`2.1.1.64`，已在元数据候选修正。IUBMB接受的末步写作quinol，而模型已有注释明确采用oxidized quinone约定；不能将本次保留表述为原生底物氧化态已证实。R18同样要保留quinone/quinol的表示限制[S01,S02]。

### 2.3 新增可确定的注释错误：COQ7反应EC

R695现有一组跨途径EC混在一起，包括`1.14.99.60`。最新IUBMB将真核NADH依赖、以quinone为底物的COQ7列为`1.14.13.253`；`1.14.99.60`页面明确指向原核quinol酶，不能继续把它当作唯一正确的真核分类[S03]。元数据候选保留R695方程及`YALI1E18269g`，将EC收窄为`1.14.13.253`。

### 2.4 真正的上游化学缺口

当前主链实际含`R969(线粒体到胞质) -> R39(胞质) -> R808(回线粒体)`。R39用半个O2完成形式羟化，却没有还原电子供体；虽然原子配平，不能据此确认为完整COQ6单加氧酶。R19把羟化与产物氧化合并，产物为quinone，而不能直接等同于COQ6孤立催化的quinol生成步骤。

COQ6 C5步骤的候选应是线粒体、具有还原电子供体的羟化：

`底物 + O2 + 2 Fd(red) + 2 H+ -> 羟化底物 + H2O + 2 Fd(ox)`

上式是电子计量框架，底物质子化形式改变时H+需相应处理；不是已经映射好本菌株Fd基因的最终SBML式[S12]。本次不把任意NADH/NADPH直接塞给COQ6，也不凭空指定Yarrowia的铁氧还蛋白/还原酶GPR。

C1步骤另有真实文献分歧：Nicoll2024体外祖先四足动物体系支持COQ4脱羧后COQ6羟化；Pelosi2024的细胞/异源表达体系支持COQ4氧化脱羧[S13,S14]。这不是可以靠挑一个EC强行消除的分歧。结论是：COQ6优先绑定经过修正的C5候选；R40/R19的C1分支保留明确的两种机制，不宣称已证明原生Yarrowia必走两步。

## 3. R305：化学错误已经确定，基因证据也有推进

### 3.1 方程与候选

原式（QH2与Q为CoQ9对，Cox/Cred为细胞色素c氧化还原态）：

`QH2 + 2 Cox + 1.5 H_N -> Q + 2 Cred + 1.5 H_P`

实算products-minus-reactants为H=-2、charge=-2。这不是浮点误差。

完整Q-cycle机制候选：

`QH2 + 2 Cox + 2 H_N -> Q + 2 Cred + 4 H_P`

H_N为基质侧，H_P为膜间隙侧。候选沿用模型胞质H作为P侧代理的区室简化，绝不把膜间隙直接等同于真实胞质。该候选在当前物种化学式下残差为0。2/4选择来自完整Q-cycle机制[S04,S05]，不是“配平只能得到2/4”；只要求配平时a与a+2有多种选择。

`model_R305_Qcycle_candidate.xml`已经写出并重新读取验证：相对参考，数学部分只改变R305这两个系数；其他方程、边界、生物量目标及基因规则未变。尚未做全网能量/表型验收，不能称正式模型修复完成。

### 3.2 名称和EC

当前名称`ferrocytochrome-c:oxygen oxidoreductase`与方程不符：该反应没有O2，受体为cytochrome c。已改为`ubiquinol-9:cytochrome-c reductase (complex III)`，反应EC保留`7.1.1.8`。原反应上的`3.4.24.64`是蛋白酶分类，不属于这条电子传递反应[S04]；移除反应层的此标签不等于删除其相关蛋白。

### 3.3 A14736不是只能写“功能完全未知”

本地映射：`YALI1A14736g -> YALI0A14806g`。原生Yarrowia CIII结构8ABF明确含`YALI0A14806p`，对应Q6CGY9、474aa[S06]。可以升级为“本地跨版本映射支持的原生复合体III成员”，但不能升级为“原生单基因敲除必定完全失活”。

PDB的peptidase注释不代表应该把该基因从CIII移走；复杂酶的成员功能与整个反应的EC属于不同层级。本次没有新做两个菌株的完整序列比对，也没有断言这个MPP-like成员必然保留蛋白酶活性。

### 3.4 当前11基因AND不等于11个固有CIII亚基

R305包含`YALI1D11769g`；模型自身将它关联到Q6C9Q0。官方PROSITE把Q6C9Q0列为cytochrome-c家族[S08]。水溶性细胞色素c是此反应的电子载体/底物，而cytochrome c1才是CIII催化核心之一[S05]。不能把两者混成一个“11亚基”名单。

处置：把载体供应依赖与酶组成规则分开记录。若旧GPR有意编码所有功能依赖，须明确命名并保留载体必要性；不能只删D11769后宣称该基因不影响呼吸。因此候选XML未擅自删除该AND项。其他辅助亚基“结构中出现”也不自动等于“每个都应硬AND”。

## 4. R1889、R2062及NDH2：错误是架构混用，不是简单无GPR

### 4.1 实际边界避免重复派发

| 反应 | 实际作用/表示 | 实际GPR | 本参考边界 |
|---|---|---|---|
| R1889 | 基质NADH，包含泵型质子项 | 空 | 0..1000 |
| R2062 | 基质NADH，无泵型质子项 | 47出现项/28不同基因，含NDH2 | 0..1000 |
| R570 | 胞质侧NADH，外侧NDH2 | YALI1F32476g | 0..1000 |
| R573 | 基质侧非泵型NADH反应 | YALI1F32476g | **0..0，已经关闭** |

R2062与R573的S列完全一致，但约束与GPR不同，不能称完整优化问题相同。应处理的是仍开放的R2062，不能再把“关闭R573”当成本轮新修复。

### 4.2 原生NDH2定位已有实验，不需要继续泛泛等待

原生Yarrowia NDH2为单亚基、非泵型、内膜外侧酶；1999删除实验和定位支持这一点[S09]。2001研究通过人为附加靶向序列，才使内部版本救援CI缺陷；外部版本不能做同样救援[S10]。

因此`YALI1F32476g`保留R570关联及外侧功能。NDH2不应作为泵型CI结构亚基硬AND进去。R573内部版本在本参考已关，符合原生模型的限定；只有明确的工程构建才应启用内部NDH2表示。

### 4.3 可直接执行与不能盲改的部分

已完成R2062的19个重复出现项删除，47变28，利用AND幂等性确认数学模型不变。**去重没有解决那28个成员本身是否正确。**

架构修订方向已确定：使用一条有正确基因关联的泵型CI，加一条外侧NDH2；R2062不应无解释地作为重复非泵路径保留。不能把现有28基因AND直接复制到R1889，因为其中NDH2明确混入，且“所有结构亚基都绝对必需”未被证明。

`YALI1A21711g -> YALI0A20680g`由本地crosswalk支持；NUPM/Q6CGB4存在于原生CI结构6YJ4[S07]。本次没有从可获取的一手记录完整闭合这个YALI1位点到Q6CGB4的独立序列桥接，因此NUPM名保持候选，不谎称已完成原生位点鉴定。

不要仅凭末尾`r`删除线粒体候选ID，也不要把旧别名自动合并；要按实际编码序列/蛋白身份处理。这些规则不取决于模拟KO是否降低生长。

### 4.4 R1889/R570的H残差需要物种约定一致

本模型NADH与NAD都存为中性形式，分别C21H29N7O14P2和C21H27N7O14P2。R1889的5H输入/4H输出因此产生H=-1、charge=-1；R570也有同类多余质子。

在完全固定这些中性物种的记账候选里，R1889可用4入/4出、R570去掉额外H来配平；在生理带电NAD约定下则需按标准反应统一物种与邻接反应[S17]。本次把二者列为条件化修订，不在不知道共享NAD全网影响时机械改一个物种或把配平当生理泵比例证明。

## 5. COQ6、COQ8、COQ9应分开结案

**COQ6（YALI1A08781g，已有F2Z6J4候选）**：保留催化家族身份；正确处置是修R39的线粒体C5羟化和电子供体，再关联该候选。当前问题不是“再加一个任意AND”，也不是一律推定无功能。原生位点精确活性与Fd配对仍未测定。

**COQ8（YALI1B20527g）**：保留ATPase/kinase-like及合成复合体辅助作用。2016研究支持ATPase与复合体稳定，后续体外体系也有构建依赖的磷酸化结果[S15,S13]。这不支持每合成1个Q固定消耗任意数量ATP，也不支持把它塞入每条反应AND。当前GEM未表示此类辅助作用属于明确的表示范围缺口。

**COQ9（YALI1F34675g）**：保留脂质结合/底物呈递与COQ7辅助作用；人蛋白结构及相关实验支持这一功能[S16]，不是本菌株逐步骤必需性的直接证据。非催化蛋白可以在已证实的必需复合体GPR中出现；本次拒绝的是没有证据的“全途径必需AND”，而不是笼统规定辅助蛋白不能入GPR。

这三者的gene KO目前不能作为正确的二元代谢依赖试验；应单列为`not represented / not testable`，不是`nonessential`。

## 6. 实际已写出的文件与验收

- `model_metadata_only_candidate.xml`：名称/EC清理、R385旧文字纠正及R2062布尔等价去重；精确比较S、物种属性、边界、目标、参数、基因集合与规范化Boolean后相同。
- `model_R305_Qcycle_candidate.xml`：基于前者，仅新增R305的2/4质子机制候选；元素/电荷精确配平。原GPR保留，P侧代理限制写入注释。
- `reaction_audit.tsv/json`：19条实际方程、物种、边界、基因及残差。
- `gene_identity_and_roles.tsv`：所有目标关联位点的本地旧ID、实际当前关联、证据层级；未个别解决的CI/III成员明确标记，不假装逐基因全确认。
- `applied_metadata_changes.json`：实际已应用到候选的修改。
- `scientific_change_spec.json`：CI架构、COQ6及辅助依赖等尚未激活的准确修订范围与理由。
- `validation.json`：实算身份、计数、等价比较和边界；优化调用0。
- `audit_and_build.py`：可对同一SHA输入复算；默认拒绝不同模型版本。

本次做的是静态反应/GPR科学校订。没有将读回XML称为libSBML全面认证，没有重跑旧轨迹，没有声称基因预测准确率提高。该工作的剩余断点已缩小到明确对象：CI逐位点/必要亚基规则、COQ6的本菌株电子供体与C1路径分歧；它们不阻止已确认的名称、EC、重复逻辑与R305化学错误交付。

## 来源

- [S01] IUBMB EC2.1.1.64 — https://iubmb.qmul.ac.uk/enzyme/EC2/1/1/64.html
  适用范围：COQ3 terminal O-methylation; accepted quinol reaction; not native Yarrowia substrate proof
- [S02] IUBMB EC2.1.1.201 — https://iubmb.qmul.ac.uk/enzyme/EC2/1/1/201.html
  适用范围：COQ5 C-methylation; accepted quinol substrate
- [S03] IUBMB EC1.14.13.253 — https://iubmb.qmul.ac.uk/enzyme/EC1/14/13/253.html
  适用范围：Eukaryotic COQ7 NADH-dependent quinone hydroxylase; created2024
- [S04] IUBMB EC7.1.1.8 — https://iubmb.qmul.ac.uk/enzyme/EC7/1/1/8.html
  适用范围：ComplexIII electron transfer and proton translocation; do not infer unique coupling from balance alone
- [S05] Wieferig and Kuehlbrandt2023 native Yarrowia CIII — https://journals.iucr.org/m/issues/2023/01/00/rq5008/index.html
  适用范围：DOI10.1107/S2052252522010570; native CIII structures and modified Q-cycle
- [S06] PDB8ABF — https://www.rcsb.org/structure/8ABF
  适用范围：Native Yarrowia CIII; entity5 YALI0A14806p/Q6CGY9; ten unique proteins
- [S07] PDB6YJ4 — https://www.rcsb.org/structure/6YJ4
  适用范围：Native Yarrowia CI; entity24 NUPM/Q6CGB4; this page alone does not prove every current YALI1 alias
- [S08] PROSITE PS51007 — https://prosite.expasy.org/PS51007
  适用范围：Cytochrome-c family; CYC_YARLI/Q6C9Q0
- [S09] Kerscher etal1999 — https://pubmed.ncbi.nlm.nih.gov/10381390/
  适用范围：Native Yarrowia NDH2 is external; DOI10.1242/jcs.112.14.2347
- [S10] Kerscher etal2001 — https://pubmed.ncbi.nlm.nih.gov/11719558/
  适用范围：Engineered internal NDH2 rescues CI deficiency; native external does not; DOI10.1242/jcs.114.21.3915
- [S11] ENZYME EC1.6.5.9 — https://enzyme.expasy.org/EC/1.6.5.9
  适用范围：Non-electrogenic NDH2; F2Z699 NDH2_YARLI
- [S12] IUBMB EC1.14.15.45 — https://iubmb.qmul.ac.uk/enzyme/EC1/14/15/45.html
  适用范围：COQ6 C5 monooxygenase requires reduced ferredoxin; not native locus validation
- [S13] Nicoll etal2024 — https://www.nature.com/articles/s41929-023-01087-z
  适用范围：DOI10.1038/s41929-023-01087-z; ancestral tetrapod in-vitro reconstruction; separate C1 decarboxylation and hydroxylation; not native Yarrowia
- [S14] Pelosi etal2024 — https://pubmed.ncbi.nlm.nih.gov/38295803/
  适用范围：DOI10.1016/j.molcel.2024.01.003; COQ4 oxidative decarboxylation in cells/heterologous systems; alternative C1 mechanism
- [S15] Stefely etal2016 — https://www.sciencedirect.com/science/article/pii/S109727651630288X
  适用范围：DOI10.1016/j.molcel.2016.06.030; COQ8 ATPase and complex stability; no arbitrary stoichiometric ATP/CoQ coefficient
- [S16] Lohman etal2014 — https://pubmed.ncbi.nlm.nih.gov/25339443/
  适用范围：DOI10.1073/pnas.1413128111; human COQ9 lipid binding and COQ7 interaction; mouse phenotype; not native Yarrowia obligatory GPR
- [S17] IUBMB EC7.1.1.2 — https://iubmb.qmul.ac.uk/enzyme/EC7/1/1/2.html
  适用范围：Pumping complexI physiological reaction; microspecies bookkeeping must be consistent
