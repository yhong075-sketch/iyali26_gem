# 暂定参考：69条零边界反应及原因

69条中，68条在模型构建起点和暂定参考中均已关闭，1条R612由PO1f菌株配置关闭。三个CFSA条件的69条集合完全一致，1500个保存样本也与零通量一致。该统计不含lipid-unlump。

问题与边界：查清直接关闭来源及单条恢复正向通量的必要条件。仅读取XML、运行配置、现有源码和保存样本；不运行求解器、不修改模型，不把历史关闭动机写成已证实事实。

数学上的直接原因是 FBA 的边界约束 lᵢ ≤ vᵢ ≤ uᵢ：当 lᵢ=uᵢ=0 时，任何可行解都必须有 vᵢ=0。更换目标函数、延长采样或降低采样阈值均不能越过这条边界。

代码顺序为加载参考模型 → 应用培养基 → 应用PO1f overlay → 添加临时总脂质需求并设置CFSA条件。68条的零边界在前述起点已存在；R612由overlay明确赋值为(0,0)。在这69条中，没有培养基新增的全关闭项；这不表示培养基对模型其他反应无影响。

R612的基因为 **YALI1E31685g — URA3 — 乳清苷-5′-磷酸脱羧酶（整理注释；本次复核模型注释与配置，未重新核实原始实验文献）**。代码用关闭R612表示配置中的ura3-302背景。旧配置说明仍写“无GPR”，但实际暂定参考已有该GPR；说明文字已过时，不改变本次实际边界赋值的判断。

**68条继承关闭的原始生物学动机尚未核实。** 已检查构建起点与暂定参考中这些反应的notes；除R612外未见解释关闭原因的文字。液泡区室简化、弃用旧配方或保留备用反应都只能作为候选解释，不能仅凭表中功能名称认定关闭合理。

下面48条还存在明确的正向局部供需缺口：保持其他反应和边界不变，仅允许该反应正向通量时，至少一个底物没有其他允许生成方向，或产物没有其他允许消耗方向。稳态守恒会再次强制该方向为零。这是静态必要条件，不是解除多条反应后的结果，也不证明反向同样阻塞。没有发现这种缺口的反应仍可能被更远的网络约束卡住。

| 类别 | 条数 |
|---|---:|
| 其他代谢 | 12 |
| 其他脂质 | 9 |
| 菌株配置 | 1 |
| 液泡转运/ATP酶 | 24 |
| 旧伪反应 | 3 |
| 鞘脂 | 20 |

全表每行一个原始反应ID。名称来自本版模型，仅说明模型分配，不认证原生酶功能。所有运行边界均为(0,0)；“继承关闭”表示直接来源已核实，关闭动机未决。区室：cy胞质、mi线粒体、nu细胞核、va液泡、em内质网膜、gm高尔基体膜、go高尔基体、ex胞外、en细胞包膜。

| 反应 | 模型名称 | 区室 | 关闭来源 | 单独开启正向时的附加阻塞 |
|---|---|---|---|---|
| R17 | 2-deoxy-D-arabino-heptulosonate 7-phosphate synthetase | mi | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R26 | 2-oxo-4-methyl-3-carboxypentanoate decarboxylation | cy | 继承关闭 | 产物无消耗方向：m81[C_cy] (4-methyl-2-oxopentanoate_C6H10O3) |
| R91 | acetyl-CoA synthetase | mi | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R92 | acetyl-CoA synthetase nuclear | nu | 继承关闭 | 产物无消耗方向：m210[C_nu] (diphosphate_H4O7P2) |
| R352 | glycerol-3-phosphate/dihydroxyacetone phosphate acyltransferase | em | 继承关闭 | 产物无消耗方向：m575[C_em] (1-acyl-sn-glycerol 3-phosphate) |
| R491 | isoleucine transaminase | cy | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R520 | leucine transaminase | cy | 继承关闭 | 产物无消耗方向：m81[C_cy] (4-methyl-2-oxopentanoate_C6H10O3) |
| R573 | NADH:ubiquinone oxidoreductase | mi | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R612 | orotidine-5'-phosphate decarboxylase | cy | PO1f菌株配置 | 未检出一阶供需缺口；能否通量仍未确定 |
| R696 | reduced glutathione via ABC system | cy, va | 继承关闭 | 产物无消耗方向：m925[C_va] (glutathione_C10H17N3O6S) |
| R746 | taurcholate via ABC system | cy, va | 继承关闭 | 底物无生成方向：m965[C_cy] (taurocholic acid_C26H44NO7S)；产物无消耗方向：m966[C_va] (taurocholic acid_C26H44NO7S) |
| R795 | V-ATPase, vacuole | cy, va | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R796 | valine transaminase | cy | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R853 | glutathione transport | cy, va | 继承关闭 | 产物无消耗方向：m925[C_va] (glutathione_C10H17N3O6S) |
| R859 | glycogen vacuolar 'transport' via autophagy | cy, va | 继承关闭 | 产物无消耗方向：m1058[C_va] (glycogen_) |
| R869 | L-asparagine transport | cy, va | 继承关闭 | 产物无消耗方向：m1066[C_va] (L-asparagine_C4H8N2O3) |
| R870 | L-aspartate transport | cy, va | 继承关闭 | 底物无生成方向：m1066[C_va] (L-asparagine_C4H8N2O3) |
| R871 | L-aspartate transport | cy, va | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R876 | L-glutamate transport | cy, va | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R878 | L-glutamine transport | cy, va | 继承关闭 | 产物无消耗方向：m1073[C_va] (L-glutamine_C5H10N2O3) |
| R879 | L-glutamine transport | cy, va | 继承关闭 | 底物无生成方向：m1073[C_va] (L-glutamine_C5H10N2O3) |
| R883 | L-isoleucine transport | cy, va | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R887 | L-leucine transport, vacuoluar | cy, va | 继承关闭 | 产物无消耗方向：m1080[C_va] (L-leucine_C6H13NO2) |
| R888 | L-leucine transport, vacuoluar | cy, va | 继承关闭 | 底物无生成方向：m1080[C_va] (L-leucine_C6H13NO2) |
| R897 | L-tyrosine transport | cy, va | 继承关闭 | 产物无消耗方向：m1089[C_va] (L-tyrosine_C9H11NO3) |
| R898 | L-tyrosine transport | cy, va | 继承关闭 | 底物无生成方向：m1089[C_va] (L-tyrosine_C9H11NO3) |
| R918 | putrescine transport | cy, va | 继承关闭 | 产物无消耗方向：m1102[C_va] (putrescine_C4H12N2) |
| R925 | spermidine transport | cy, va | 继承关闭 | 产物无消耗方向：m1110[C_va] (spermidine_C7H19N3) |
| R927 | spermine transport | cy, va | 继承关闭 | 产物无消耗方向：m1112[C_va] (spermine_C10H26N4) |
| R1027 | bicarbonate formation | ex | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R1040 | cholestenol delta-isomerase, lumped reaction | cy | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R1138 | glucose transport, vacuolar | cy, va | 继承关闭 | 产物无消耗方向：m1236[C_va] (D-glucose_C6H12O6) |
| R1163 | H+ diffusion | cy, va | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R1297 | phosphate transport | cy, va | 继承关闭 | 产物无消耗方向：m1354[C_va] (phosphate_H3O4P) |
| R1349 | trehalose vacuolar transport | cy, va | 继承关闭 | 产物无消耗方向：m1380[C_va] (trehalose_C12H22O11) |
| R1363 | water diffusion | cy, va | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R1371 | lipid pseudoreaction | cy | 继承关闭 | 底物无生成方向：m1394[C_cy] (ergosterol ester_C29H43O2R) |
| R1387 | yeast 6 biomass pseudoreaction | cy, en | 继承关闭 | 底物无生成方向：m1400[C_cy] (lipid_) |
| R1710 | yeast 8 biomass pseudoreaction | cy, en | 继承关闭 | 底物无生成方向：m1400[C_cy] (lipid_) |
| R1720 | DAG kinase | em | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R1725 | choline/ethanolaminephosphotransferase | em | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R1726 | cholinephosphotransferase | em | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R1772 | diacylglycerol acyltransferase | em | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R1775 | diacylglycerol pyrophosphate phosphatase | gm | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride) |
| R1795 | diglyceride transport, ER membrane-Golgi membrane | em, gm | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride) |
| R1815 | inositolphosphotransferase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride) |
| R1816 | inositolphosphotransferase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride) |
| R1817 | inositolphosphotransferase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride) |
| R1818 | inositolphosphotransferase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride) |
| R1819 | inositolphosphotransferase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride) |
| R1820 | inositolphosphotransferase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride) |
| R1821 | inositolphosphotransferase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride) |
| R1822 | inositolphosphotransferase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride) |
| R1823 | inositolphosphotransferase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride) |
| R1824 | inositolphosphotransferase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride) |
| R1825 | IPC synthase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride)；产物无消耗方向：m1707[C_go] (inositol-P-ceramide B'-(C24)_) |
| R1826 | IPC synthase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride)；产物无消耗方向：m1708[C_go] (inositol-P-ceramide B'-(C26)_) |
| R1827 | IPC synthase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride)；产物无消耗方向：m1709[C_go] (inositol-P-ceramide C-(C24)_) |
| R1828 | IPC synthase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride)；产物无消耗方向：m1710[C_go] (inositol-P-ceramide C-(C26)_) |
| R1829 | IPC synthase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride)；产物无消耗方向：m1711[C_go] (inositol-P-ceramide A-(C24)_) |
| R1830 | IPC synthase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride)；产物无消耗方向：m1712[C_go] (inositol-P-ceramide A-(C26)_) |
| R1831 | IPC synthase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride)；产物无消耗方向：m1713[C_go] (inositol-P-ceramide B-(C24)_) |
| R1832 | IPC synthase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride)；产物无消耗方向：m1714[C_go] (inositol-P-ceramide B-(C26)_) |
| R1833 | IPC synthase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride)；产物无消耗方向：m1715[C_go] (inositol-P-ceramide D-(C24)_) |
| R1834 | IPC synthase | gm, go | 继承关闭 | 产物无消耗方向：m1704[C_gm] (diglyceride)；产物无消耗方向：m1716[C_go] (inositol-P-ceramide D-(C26)_) |
| R1843 | 1-acyl-sn-gylcerol-3-phosphate acyltransferase | em | 继承关闭 | 底物无生成方向：m575[C_em] (1-acyl-sn-glycerol 3-phosphate) |
| R1899 | O-succinylhomoserine lyase (elimination) | cy | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R1950 | dihydrofolate:NADP+ oxidoreductase | cy | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |
| R1951 | dihydrofolate:NADP+ oxidoreductase | mi | 继承关闭 | 未检出一阶供需缺口；能否通量仍未确定 |

实际含义与下一步：

- 20条鞘脂反应R1815–R1834都产生m1704[C_gm]（高尔基体膜DAG），当前没有其他开放反应消耗该池；R1825–R1834还产生没有开放去路的对应IPC池。应先核对这些池与其他区室/版本的身份、运输与需求，不能靠将上界统一改成1000来修复。
- 液泡转运中的腐胺R918、亚精胺R925和精胺R927，正向都会把物质送入没有当前消耗方向的液泡池；这是稳态模型中缺少周转/出口/需求的表现。真实积累属于动态问题，不能据此自动加无依据的sink。
- R1371及R1387/R1710是旧脂质/生物量伪反应：模型名称和当前关闭状态已核实；“为避免重复生物量而关闭”尚无直接动机记录。打开它们可能改变可行域，不能作为普通失活酶处理。
- 57/69条已有GPR，12条没有；GPR是否存在不解释边界为何被设为零，也不证明反应在本菌株中的功能。
- 对有供需缺口的组，优先做区室/池身份和连接关系核对；对没有这类局部缺口的组，先追溯关闭依据，再决定是否进行单条、有对照的放宽边界和FVA。仅在原零边界下做FVA会返回已知的零范围，不能回答解除限制后是否恢复。

算法局限：当前检查只利用计量符号和现有边界判断必要供需方向，没有求解Sv=0的完整耦联系统；“存在某个允许方向”不保证该供给/消耗反应实际可行。

实际使用skill：govern-agentic-research（中性范围、来源与独立审阅）、gene-identity-function（R612身份与模型角色分开）、Ponytail（复用保存数据，仅做离线核对）。未执行原始论文核验、蛋白功能重判、反向/联合解除边界FVA或模型计算。

独立来源审阅已完成：审阅者直接打开三个场景JSON、两份XML、配置与匹配实际运行记录的源码副本，并读取三个NPZ；69/69行已审，0遗漏、0重复、0冲突、0未核查行，13/13来源身份记录匹配。69条在1500个保存样本中的最大绝对通量均为精确0。分组供需缺口核对为液泡18/24、鞘脂20/20、旧伪反应3/3、其他脂质4/9、其他代谢3/12、菌株0/1；57/69有GPR。无必须修正项。审阅支持的是本次静态事实与必要条件推论，不认证68条历史关闭的生物学动机或论文层面的蛋白功能。

## 涉及的实验共识必需基因

按项目现有1612项实验共识正例表精确匹配系统ID：69条反应的GPR共涉及49个不同基因，其中14个在实验正例表，关联24条不同关闭反应。这里的“实验必需”指现有项目标签，没有重新构建三筛选共识，也不是独立验证。反应集合为R795、R1815–R1834、R1843、R1950及R1951。

下表各基因的本物种正式名称均未独立核实；功能仅为现有同源注释及模型/GPR赋值线索（model/GPR assignment only），不是本次原生蛋白功能实验确认。历史KO来自2026-08-26实际执行模型的保存结果，不能称为本次暂定参考重新求解。主阈值为KO/WT < 0.1；表中约100%的数值保留了原始浮点值于后台证据。

| 系统ID（正式名称未核实） | 蛋白功能线索 | 对应关闭反应 | 历史KO/WT | 历史分类 |
|---|---|---|---:|---|
| YALI1A09766g | V-ATPase催化A亚基候选 | R795 | ≈100% | FN |
| YALI1A11258g | V-ATPase亚基候选，VMA5同源注释 | R795 | ≈100% | FN |
| YALI1E10492g | V-ATPase亚基候选，VMA13弱相似注释 | R795 | ≈100% | FN |
| YALI1E12482g | V-ATPase膜亚基候选，STV1/VPH1同源注释 | R795 | ≈100% | FN |
| YALI1E14125g | V-ATPase亚基候选，VMA10同源注释 | R795 | ≈100% | FN |
| YALI1E32332g | V-ATPase亚基候选，VMA2同源注释 | R795 | ≈100% | FN |
| YALI1E37063g | V-ATPase膜亚基候选 | R795 | ≈100% | FN |
| YALI1F13017g | H⁺-ATPase样蛋白，按原注释 | R795 | ≈100% | FN |
| YALI1F20965g | V-ATPase E亚基候选 | R795 | ≈100% | FN |
| YALI1F21690g | V-ATPase亚基候选，VMA6同源注释 | R795 | ≈100% | FN |
| YALI1F31854g | V-ATPase膜蛋白脂质亚基候选 | R795 | ≈100% | FN |
| YALI1E19647g | AUR1样蛋白；模型分配IPC合成及肌醇磷酸转移 | R1815–R1834（20条） | ≈100% | FN |
| YALI1E22736g | SLC1样脂酰转移酶候选 | R1843 | ≈100% | FN |
| YALI1C17947g | DFR1样二氢叶酸还原酶候选 | R1950、R1951 | 0% | TP |

FN表示实验正例而历史模型预测非必需，TP表示两者均判必需。VMA、AUR1、SLC1、DFR1等仅指源注释中的其他物种同源蛋白，不作为这些YALI1基因的已确立原生名称。

当前GPR的静态布尔核查给出三种不同机制：

1. R795关联的11个实验正例基因，单独去掉任一个，都不会使当前模型任何反应的完整GPR变成false。R794/R795含替代OR分支；因此R795的零边界不是这些FN的唯一问题，单独打开R795也不会消除这种单基因KO的逻辑替代。
2. YALI1E19647g仅分配给上述20条已关闭反应；标准GPR单基因KO只会再次关闭这20条，不能改变当前通量可行域。YALI1E22736g在R1843及开放的R1846中均有OR替代，单KO同样不关闭额外反应。以上13个FN的静态解释不需要重跑求解器；它们不证明GPR在生物学上正确。
3. YALI1C17947g还控制开放的R271和R272，KO会关闭这两条，因此“关联关闭反应”不意味着该基因在模型里一定非必需。历史KO/WT为0，与其仍有开放反应受影响一致；本次没有重新检验生长因果。

若“essential”仅指历史模型预测而不要求实验正例，这69条的关联基因共有3个被预测必需：上表YALI1C17947g，以及 **YALI1C00230g — 正式名称未核实 — 甘油-3-磷酸/二羟丙酮磷酸酰基转移酶候选（model/GPR assignment only）** 和 **YALI1D01431g — 正式名称未核实 — 支链氨基酸转氨及相关反应的模型候选（model/GPR assignment only）**。后两者在现有实验正例表中未标注，不能当作实验非必需；分别关联关闭的R352，以及R26/R491/R520/R796，同时另有关联的开放反应。

本节仅读取标签表、历史KO表及当前GPR，没有新求解、基因映射替换、标签修订或模型修改。使用Spreadsheets skill的只读交叉核对流程。

本节独立审阅：14/14实验正例关联、24/24反应关联、14/14单KO布尔结果、14/14历史分类均核验一致，0冲突。审阅者分别从实际暂定参考XML和保存场景JSON计算GPR结果，并直接读取共识表与历史运行清单。未核实原生蛋白功能或重新复核论文来源；本节结论保持候选功能与历史预测的限定。
