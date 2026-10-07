# R794/R795 GPR独立来源审阅

核验时间：2026-09-09T23:00:06.389215+00:00；本次新增GEM求解、BLAST、AlphaFold作业和科学模型修改均为0。直接打开原XML、原缓存和版本蛋白记录、论文原页面快照、原始摘要及AlphaFold坐标/PAE核验；没有把作者报告或代理一致意见当作原始证据。

**建议：将现有12AND OR2AND标为缺乏完整替代泵证据的待修订规则，支持进入提案，尚不支持实施唯一新GPR。** 后组是a与c-like两种膜组件候选；前组另有a候选。它们的家族注释与完整V1/V0泵的功能分工不支持“二蛋白代替整套复合体”。更有依据的提案方向是共同组件的AND，以及证据充分的同类亚型位置OR。尚缺本物种定位、互换性、所有共同成员要求与对应培养条件的生长依赖；不能把14个全部AND、直接删后组或保证修后essential。

发表来源现已核实：出版商mmc6.xml与本地iYli21完整字节相同，论文原XML同时核对了DOI、PII及补充材料6。旧iYali4.1.2的13AND与现分组差异也成立，但当时精确修订理由仍未知。旧诊断的全关联反应关闭后生长约1仅属于其各自旧模型与保存条件，不能冒充本次参考复现。

## 身份与功能证据

以下功能均是保留快照的自动注释候选，证据状态为uncharacterized／未原生实验验证；GPR赋值本身不提升功能证据。所有14份快照与各自缓存、ORF及完整序列SHA独立匹配。

| 系统ID | 名称状态 | 蛋白功能线索 | 证据等级 |
|---|---|---|---|
| YALI1A11258g | 未确立原生正式名 | V-type proton ATPase subunit C候选 | 自动数据库注释；非原生实验验证 |
| YALI1E14125g | 未确立原生正式名 | V-type proton ATPase subunit G候选 | 自动数据库注释；非原生实验验证 |
| YALI1D00581g | 未确立原生正式名 | V-type proton ATPase subunit D候选 | 自动数据库注释；非原生实验验证 |
| YALI1E10492g | 未确立原生正式名 | V-type proton ATPase subunit H候选 | 自动数据库注释；非原生实验验证 |
| YALI1A09766g | 未确立原生正式名 | V-type proton ATPase catalytic subunit A候选 | 自动数据库注释；非原生实验验证 |
| YALI1F20965g | 未确立原生正式名 | ATPase, V1/A1 complex, subunit E候选 | 自动数据库注释；非原生实验验证 |
| YALI0E16192g | 未确立原生正式名 | V-type proton ATPase subunit F候选 | 自动数据库注释；非原生实验验证 |
| YALI1F31854g | 未确立原生正式名 | V-type proton ATPase proteolipid subunit候选 | 自动数据库注释；非原生实验验证 |
| YALI1F21690g | 未确立原生正式名 | V-type proton ATPase subunit候选 | 自动数据库注释；非原生实验验证 |
| YALI1E37063g | 未确立原生正式名 | V-type proton ATPase proteolipid subunit候选 | 自动数据库注释；非原生实验验证 |
| YALI1E32332g | 未确立原生正式名 | Vacuolar proton pump subunit B候选 | 自动数据库注释；非原生实验验证 |
| YALI1E12482g | 未确立原生正式名 | V-type proton ATPase subunit a候选 | 自动数据库注释；非原生实验验证 |
| YALI1F38820g | 未确立原生正式名 | V-type proton ATPase subunit a候选 | 自动数据库注释；非原生实验验证 |
| YALI1F13017g | 未确立原生正式名 | V-ATPase proteolipid subunit C-like domain-containing protein候选 | 自动数据库注释；非原生实验验证 |

13个YALI1条目来自UP000182444；旧YALI0E16192g的Q6C5Q2来自UP000001300、taxon284591。该项没有本轮确立的W29对应，保留旧ID。两个后组候选及前组a候选的版本明确NCBI GenPept记录均与目标蛋白序列相同，记录CLIB89(W29)和相应locus_tag，但仍属conceptual translation/hypothetical protein。两个后组候选与旧位点蛋白序列完全相同仅支持蛋白序列桥接，不等于基因组位点或表型等价；前组a候选825aa与旧820aa序列不同，已保留。

两个当前UniProt条目的DELETED理由都是“不再属于参考蛋白组”，不能解读为基因或蛋白不存在。已存在的AlphaFold单体预测与锁定旧序列完全一致；模型为2022年记录、版本6、Monomer v2.0，不是本次新预测。独立重算CA平均pLDDT为82.08835与90.26602，低于50的残基数分别59与0；平均PAE为14.02922与6.39994 Å。API总体分数82.06/90.25另列。高预测置信度不能证明组装、区室、催化或完整泵替代性；本快照不审后续结构拟合。

## 文献限定与反证

直接核读2019论文ATPase活性实验、图1图注和方法，确认活性比较对象是保留共同V1/V0组件的完整酵母复合体；两个生物重复各三次测量。1994原始摘要支持条件性a亚型部分互补，1997原始摘要支持c″作为复杂泵组件及proteolipid非简单冗余。2011本物种转录组的名称来自酿酒酵母同源，不能替代原生酶学验证。1993本物种离体液泡总活性不能识别本题基因或验证二蛋白泵。

2015亚基e体外去除例只反驳“一切结构组件在所有情境都必须参与催化”的一般推断；**e不在本题14基因内，不能作为本题某个基因非必需或禁止14AND的直接反证。** 本题不应直接14AND的原因是两个a候选的区室/亚型要求及其余成员必要性仍未确立。[2019完整复合体实验](https://pmc.ncbi.nlm.nih.gov/articles/PMC6462096/)、[1997proteolipid原始摘要](https://pubmed.ncbi.nlm.nih.gov/9030535/)支持组件层审查，不支持跨物种直接套用。

## 原子判定

总数23；已审23；支持20；未决3（其中unverified 2、unsupported 1）；反证0；未审0。unsupported为“二蛋白已有完整替代泵证据”这一命题，未决为原始修订理由与唯一最终GPR；不把限定未找到证据写成生物学绝对不存在。

| Claim ID | 判定 | 限定主张 |
|---|---|---|
| PROV1 | supported | 图源R794/R795与本地iYli21的GPR和各自边界相同。 |
| PROV2 | supported | 本仓最早可查2026-03-19提交内iYli21文件与当前P2字节相同。 |
| PROV3 | supported | 本地iYali v4.1.2两泵为13AND；按现存ID表，现分支移出原共同c-like候选并新增另一a候选。 |
| PROV4 | supported | 本地Yeast-GEM两泵的OR分支大量共享成员，非本题不相交12/2结构。 |
| PROV5 | supported | 旧模型诊断保存了全关联反应关闭后生长比值约1，并未形成已批准patch/回归结论。 |
| PROV6 | supported | 出版商mmc6与本地iYli21字节相同，原文XML将mmc6列作该论文补充材料6。 |
| PROV7 | unverified | 已查到该GPR精确原始修订理由和生物学依据。 |
| L01 | supported | 后组两位点的历史精确ORF条目描述V0组件候选，不是二蛋白完整ATPase。 |
| L02 | supported | 旧位点名称来自同源注释，2011转录组本身没有认证原生泵功能或当前位点映射。 |
| L03 | supported | 2019 a亚型比较使用保留共同V1和其余V0组件的完整酵母复合体。 |
| L04 | supported | 1994 a亚型部分互补发生在保留共同亚基背景，且单/双缺失的条件性生长不同。 |
| L05 | supported | 1997 S.cerevisiae c-double-prime KO损失泵活性/装配，三个proteolipid并非简单冗余。 |
| L06 | supported | 酵母亚基e可从已装配纯化酶中被去除而保留体外耦联泵H+，说明组件与催化必要性需分情境。 |
| L07 | supported | 1993 Yarrowia离体液泡具有ATPase/质子梯度等生理证据，但没有验证本题二基因分支。 |
| SEQ1 | supported | 14位点的保留快照与各自原缓存对象、精确ORF字段和完整序列SHA一致。 |
| SEQ2 | supported | 后组804aa候选为a类、196aa候选为c-like类，前组825aa候选同为a类。 |
| DB1 | supported | 两个当前UniProt条目被标记DELETED，理由是Not part of a reference proteome。 |
| AF1 | supported | 已存在的两份AF单体模型/API序列与锁定旧快照序列完全一致。 |
| AF2 | supported | 按PDB CA B因子与原PAE矩阵独立重算的统计匹配结构摘要；API总体指标单列。 |
| ID01 | supported | 3条版本明确的NCBI AOW蛋白序列与目标快照完全相同，GenPept源记录CLIB89(W29)及目标locus_tag。 |
| ID02 | supported | 两个后组候选与所映射旧位点蛋白序列相同；前组a候选与旧蛋白长度/序列不同。 |
| SYN1 | unsupported | 当前二基因分支已有证据证明能独立完成ATP水解和泵H+，足以替代前12成员复合体。 |
| SYN2 | unverified | 现有资料足以唯一确定可立即实施的完整14位点新GPR。 |

逐条来源定位、完整SHA、修正核验与范围见[claims_reviewed.json](claims_reviewed.json)。审阅识别的单缓存来源与AF分数字段问题已由作者补正并重新核实；当前没有影响上述限定结论的必须补正项。BLAST和后续结构拟合未进入本快照分母，不能称已审核。科学模型/GPR接纳仍需相应具体授权与证据。
