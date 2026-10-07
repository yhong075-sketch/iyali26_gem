# ATP 候选：独立化学身份与来源审核

2026-09-24。仅复核指定基线、工作区来源和官方公开化学/酶学资料；没有优化调用，没有修改模型、共享代谢物、GPR、培养或整理数据。完整输入 SHA、真实字段与 E0/E1/E2/E3 **预期定义**保存在 [chemical_signatures.json](chemical_signatures.json)，这些定义不是本子任务生成/执行的候选 XML。构建及求解验收由主任务另行记录。

## 裁决

三项 H/charge 修订均有**当前明确物种表示下的化学依据**，可进入隔离候选：R_NTP3pp 去掉反应物 1H+；R_NTP7 去掉产物 1H+；r0242 正向补产物 2H+。这不是从残差倒推任意质子项：先用 XML 的化学身份、明确保存的 formula/charge 对及官方物种定义确定每种形式，再展开普通水解和酸碱形式转换。

三个既定方向修订可同时保留：R_PGAM1_PhosHydro 和 R_NTP3pp 上界为 0；R_NTP7 下界为 0。它们是磷酸酶反应类别下的有据候选，不能据此认定原生精确酶活已验证，也不能单独担保全模型能量修复。

## 实际身份与质子化约定

以下均为胞质 `C_cy`，不是按名称推测的物种形式。核苷酸/Pi 的 XML notes 明确记录 `preserved_existing_pair`，保留原有中性公式/电荷并拒绝 MetaNetX 提出的带电替代。DHAP 现有 formula、charge、InChI `/p-2` 相互一致。

| 实际 ID | 物种形式 | 执行公式 / charge | 独立官方核对 |
|---|---|---|---|
| m266[C_cy] / m268[C_cy] | 中性 GTP / GDP | C10H16N5O14P3 / 0；C10H15N5O11P2 / 0 | [GTP](https://www.ebi.ac.uk/chebi/CHEBI:15996)、[GDP](https://www.ebi.ac.uk/chebi/CHEBI:17552) |
| m439[C_cy] / m11[C_cy] | 中性 UTP / UDP | C9H15N2O15P3 / 0；C9H14N2O12P2 / 0 | [UTP](https://www.ebi.ac.uk/chebi/CHEBI:15713)、[UDP](https://www.ebi.ac.uk/chebi/CHEBI:17659) |
| m141[C_cy] / m143[C_cy] | 中性 ATP / ADP | C10H16N5O13P3 / 0；C10H15N5O10P2 / 0 | [ATP](https://www.ebi.ac.uk/chebi/CHEBI:15422)、[ADP](https://www.ebi.ac.uk/chebi/CHEBI:16761) |
| m406[C_cy] / m500[C_cy] | 中性 CTP / CDP | C9H16N3O14P3 / 0；C9H15N3O11P2 / 0 | [CTP](https://www.ebi.ac.uk/chebi/CHEBI:17677)、[CDP](https://www.kegg.jp/entry/C00112)；CDP ChEBI 本次不可取，0 电荷取自明确 XML 约定 |
| m35[C_cy] | Pi 实际表示为磷酸 H3PO4 | H3O4P / 0 | [磷酸](https://www.ebi.ac.uk/chebi/CHEBI:26078) |
| m1881[C_cy] | DHA，dihydroxyacetone | C3H6O3 / 0 | [DHA](https://www.ebi.ac.uk/chebi/CHEBI:16016) |
| m456[C_cy] | DHAP 二价阴离子 | C3H5O6P / −2 | [glycerone phosphate(2−)](https://www.ebi.ac.uk/chebi/CHEBI:57642) |
| m32[C_cy] / m10[C_cy] | 水 / H+ | H2O / 0；H / +1 | [水](https://www.ebi.ac.uk/chebi/CHEBI:15377)、[hydron](https://www.ebi.ac.uk/chebi/CHEBI:15378) |

**反证/局限保留：** 中性核苷酸和 H3PO4 不等于胞质 pH 下的主要离子形式；很多交叉注释及 InChI 来自带电 MetaNetX 物种。因此不能声称当前全网是统一的 pH 7.3 微观物种模型。当前 formula/charge 对是明确执行和保护的约定，在此约定下修反应即可，不需要任意改共享物种。若将来全网改为带电表示，必须另作所有相邻反应的传播审查。ChEBI 对 DHAP dianion 的化学形式支持，也不替代 W29 浓度/pH 的测定。

## 完整反应签名：E0 → E3

负号是反应物，正号是产物。所有完整 GPR 均为 `YALI1E41893g`；该系统 ID 的正式原生名称未核实，泛酸性磷酸酶候选源于自动注释，而对这四个精确反应的关系仍是 **model/GPR assignment only**。不因修订化学而接受精确底物、定位或 GPR。缓存 UniProt A0A1D8NLA2（entry version 29，2026-01-28，sequence version 1）还包含自动 SignalP 的 N 端 1–17 aa 信号肽，进一步使“自由胞质定位”不能直接接受；这是自动预测提示，不是本次新做预测或实验定位。

| 反应 | E0 完整计量 | E0 边界 | E3 完整计量 | E3 边界 |
|---|---|---|---|---|
| R_PGAM1_PhosHydro | m4[C_cy]:−1; m35[C_cy]:−1; m107[C_cy]:+1; m32[C_cy]:+1 | [−1000,1000] | 同 E0 | [−1000,0] |
| R_NTP3pp | m268[C_cy]:−1; m35[C_cy]:−1; m10[C_cy]:−1; m266[C_cy]:+1; m32[C_cy]:+1 | [−1000,1000] | m268[C_cy]:−1; m35[C_cy]:−1; m266[C_cy]:+1; m32[C_cy]:+1 | [−1000,0] |
| R_NTP7 | m439[C_cy]:−1; m32[C_cy]:−1; m11[C_cy]:+1; m35[C_cy]:+1; m10[C_cy]:+1 | [−1000,1000] | m439[C_cy]:−1; m32[C_cy]:−1; m11[C_cy]:+1; m35[C_cy]:+1 | [0,1000] |
| r0242 | m1881[C_cy]:−1; m35[C_cy]:−1; m456[C_cy]:+1; m32[C_cy]:+1 | [−1000,1000] | m1881[C_cy]:−1; m35[C_cy]:−1; m456[C_cy]:+1; m32[C_cy]:+1; m10[C_cy]:+2 | [−1000,1000] |

E1 只用上述边界，全部计量保留 E0；E2 只用上述计量，全部边界保留 E0。两者与 E3 的逐字段完整版本在 JSON，正常对照列亦保留。

R_NTP3pp 与 R_NTP7 的中性 NTP 水解均为 `NTP + H2O → NDP + H3PO4`，不含净 H+。官方 [IUBMB EC3.6.1.5](https://iubmb.qmul.ac.uk/enzyme/EC3/6/1/5.html) 记录逐步释磷水解；既有 KEGG R00335/R00159 原文核对了精确 GTP/UTP 身份。方向裁决依据水解机制和无能量耦联的当前列，不依赖数据库等号/双箭头的排版。

r0242 的修订可显式分解为中性 DHAP 水解 `DHAP(H2) + H2O → DHA + H3PO4`，加上 `DHAP²⁻ + 2H+ → DHAP(H2)` 的酸碱表示转换，得到 `DHAP²⁻ + H2O + 2H+ → DHA + H3PO4`。这正是修订列的负向。中性/二价阴离子的两种定义可由 [ChEBI:16108](https://www.ebi.ac.uk/chebi/CHEBI:16108) 与 [ChEBI:57642](https://www.ebi.ac.uk/chebi/CHEBI:57642) 的共轭关系核对。**2H+ 是同一区室物种形式记账，不是跨膜质子运输、ATP 耦联或实际净质子消耗的实验测量。** 本阶段保留 r0242 原双向边界；若后续新见证使用无耦联合成方向，须独立裁决其方向，不能把化学配平等同于双向都合理。

PEP phosphatase 的官方定义为 PEP 水解至 pyruvate/Pi。[IUBMB EC3.1.3.60](https://iubmb.qmul.ac.uk/enzyme/EC3/1/3/60.html) 与 [Duff 等 1989](https://pmc.ncbi.nlm.nih.gov/articles/PMC1061789/) 的缓存摘要支持该酶学类别；植物研究不能当作 W29 酶活验证。当前正向无 ATP 或其他能源耦联，故阻止该合成方向是条件性候选；未从原文得到 W29 实际 ΔG，不声称任意条件绝对不可逆。

## 逐列守恒和诊断列

本次直接从固定 XML 重算：E0 的 R_NTP3pp 为 H=−1/charge=−1，R_NTP7 为 +1/+1，r0242 为 −2/−2；E2/E3 三列均全部元素与电荷零残差。PGAM 本来配平。R694/R594/R603/R2010 和现有 xMAINTENANCE 均完整配平；未建议删除或改动这些正常能量转移/维护列。

ATP、GTP、UTP、CTP 的当前中性 NTP/NDP/Pi 对均可构造 `NTP + H2O → NDP + Pi` 的完整配平诊断式，H+ 系数均为 0。使用现有 xMAINTENANCE 检查 ATP，无需另加 proton。此处仅确认化学式，未执行能量优化，也不把逐列配平当作热力学正确。

## 来源取得范围和审核覆盖

2026-09-24 实际打开上述 ChEBI 物种页面（CDP 除外），并打开 KEGG C00044/C00035/C00075/C00015/C00009/C00184/C00112 以及 IUBMB EC3.1.3.60/EC3.6.1.5。ChEBI CDP 和部分 KEGG reaction 页面本次未能在线取回；KEGG R01010/R00335/R00159 采用上一轮工作区内已下载官方原文逐条读取，文件完整 SHA 在 JSON。PMC 在线返回验证页，实际读取工作区缓存 HTML 的原文摘要，不声称重读全部扫描页。底物和化学定义审核不借用缓存文件名作为内容证据。

初始范围为 12 项声明；下方 E3/E4 新见证追加 13 项，当前累计 **total 25 | audited 25 | supported 22 | unresolved 3 | contradicted 0 | unchecked 0**。未决项是原生精确底物/定位实验验证、R72 原生 OR 的精确底物互替性及 E4 两反应原生精确活性/区室，未据此升级 GPR；具体声明、来源、限制和 verdict 见 [claim_audit.json](claim_audit.json)。本覆盖率不包含主任务优化结果或最终模型接纳。

## E3 新见证追加：R_NDP1 / R72

收到并读取主任务 `initial_validation/E3/witness.json`，本节只审核化学和来源，不冒称复现见证优化。实际 XML 和逐项来源保存在 [E3_witness_followup.json](E3_witness_followup.json)。

**R_NDP1 支持追加候选。** 当前完整式为 `m86[C_cy]:−1; m35[C_cy]:−1; m10[C_cy]:−1; m143[C_cy]:+1; m32[C_cy]:+1`，边界 [−1000,1000]，GPR 沿用上文系统 ID/证据限制。m86 的 AMP 公式 C10H14N5O7P、charge 0 与 [ChEBI:16027](https://www.ebi.ac.uk/chebi/CHEBI:16027) 相符，XML 亦明确保留中性对。[KEGG R00122](https://www.kegg.jp/entry/R00122) 与 IUBMB EC3.6.1.5 的 NDP 水解步骤一致。

因此候选是去掉 `m10:−1` 并将上界改为 0，其余物种系数、下界 −1000、GPR 不变；H/charge −1/−1 变为全部元素和电荷零残差。它恢复 `ADP + H2O → AMP + H3PO4` 的水解，禁止当前没有能源耦联的合成。`gap_fill_direction_curation.csv` 中已有 reverse/[0,1000] 的水解整理；`metadata_reaction_selection.json` 的 before 是不含 H+ 的水解单向，after 为当前消耗 H+ 的合成式/双向。当前 XML 与 after 一致，属于同一覆盖问题，而非新发现正常腺苷酸激酶需要删除。

**R72 支持化学类别层面的方向候选，保留身份冲突。** 当前完整式 `m170[C_cy]:−1; m32[C_cy]:−1; m140[C_cy]:+1; m35[C_cy]:+1`，边界 [−1000,1000]。两肌醇物种分别按 C6H19O27P7/0 与 C6H18O24P6/0 保存；当前式已经全部元素、电荷守恒，不补 H+。支持将下界改为 0，保留实际水解正向及原上界 1000。新限制仅对应这条真实计量，不因酶名扩展其他限制。

[IUBMB EC3.6.1.52](https://iubmb.qmul.ac.uk/enzyme/EC3/6/1/52.html) 描述二磷酸肌醇的末端磷酸水解，吻合实际式及原 `PROTEIN_CLASS:3.6.1.52`。1999 年 [Safrany 等原始研究摘要](https://pubmed.ncbi.nlm.nih.gov/10419486/) 以纯化重组酵母酶检测到此类水解，支持酶类别，未验证本模型的 Yarrowia 位点或特定 6-位异构体。相反，XML 的 `ec-code:2.7.4.24` 指向 [ATP/ADP 耦联的肌醇激酶](https://iubmb.qmul.ac.uk/enzyme/EC2/7/4/24.html)，该定义需要核苷酸转磷酸；不能用它为 `IP6 + Pi → IP7 + H2O` 的无耦联逆向辩护。即便某激酶能反向把磷酸转给 ADP，也不同于本列。

R72 当前 OR 的身份只作解释：YALI1B17927g — 正式原生名称未核实 — inositol hexakisphosphate/diphosphoinositol kinase 候选（A0A1D8N7N0，entry 37/sequence 1，自动/同源注释）；YALI1A04817g — 正式原生名称未核实 — NUDIX hydrolase domain-like 候选（A0A1H6Q7Y1，entry 36/sequence 1，自动预测注释）。两份缓存注释更新均为 2026-01-28。它们在 R72 中是 **model/GPR assignment only**；两个蛋白对同一底物能独立互替未确证，不据此改 GPR。m140 没有结构交叉注释，m170 的位置命名也未做新结构认定，因此不能把本次方向修订称全套原生肌醇路径接纳。

R72 的双向边界在原 `data/iyali26.xml` 和 `data/iyli21.xml` 已存在，原式还含产物 2H+，当前式不含它们；因此它不是上述新增 gap-fill 磷酸酶的同一个来历。保留此差异，不把 R72 双向自动归咎 metadata 的四列覆盖。

本追加没有认可 R284/R422 的计量：其缺式物种和残差依然需要专门证据。R72 当前方向限制不证明这两列正确，也不证明整体能量异常已消除。无证据部分不能因与 R_NDP1 的 H+ 错误在见证里抵消而接纳。


## E4 新见证追加：peroxidatic CAT2p 与 OAADCm

本节独立读取固定 XML、早期整理表、metadata 字段选择及主任务提供的 `E4_validation/E4/witness.json`。只重算化学列，不运行优化；完整前后签名、物种注释、来源 SHA 与两份原生候选缓存原文在 [E4_witness_followup.json](E4_witness_followup.json)。

| 反应 | 固定模型完整式 | 原边界 | 候选完整式 | 候选边界 |
|---|---|---|---|---|
| R_CAT2p | m275[C_cy]:−1; m222[C_pe]:−1; m173[C_cy]:+1; m375[C_pe]:+2 | [−1000,1000] | 同原式 | [0,1000] |
| R_OAADCm | m44[C_mi]:−1; m6[C_mi]:−1; m77[C_mi]:+1; m28[C_mi]:+1 | [−1000,1000] | m44[C_mi]:−1; m6[C_mi]:−1; m77[C_mi]:+1 | [−1000,0] |

**CAT2p 是乙醇的 peroxidatic 氧化，不能写成过氧化氢歧化。** 其正向为 `ethanol[cy] + H2O2[pe] → acetaldehyde[cy] + 2H2O[pe]`；公式分别 C2H6O、H2O2、C2H4O、H2O，全部 charge 0，原式全部元素和电荷已配平。提供见证的负向通量 −0.769230769 则是 `acetaldehyde + 2H2O → ethanol + H2O2`。本列完全没有 O2；见证里后续氧气形成由别列承担。[IUBMB EC1.11.1.6](https://iubmb.qmul.ac.uk/enzyme/EC1/11/1/6.html) 在主歧化定义之外明确说明部分 catalase 有以乙醇作供氢体的 peroxidase 活性。实际打开的 [Oshino 等 1973 原始研究摘要](https://pubmed.ncbi.nlm.nih.gov/4720713/) 研究大鼠肝 catalase 利用 H2O2 氧化乙醇，支持此活动类别和氧化方向；未重读全部扫描实验页，不能外推为 Yarrowia 实验。故 LB0 保留乙醇氧化是有据的条件性方向候选，不需改计量。

早期 gap-fill 记录已明确写出乙醇/过氧化氢氧化及禁止反向合成两者，边界 [0,1000]。metadata selection 只选了 bounds，将它恢复至 [−1000,1000]。这项具体覆盖可以核实；不由此猜测作者动机。

CAT2p 的 `ethanol/acetaldehyde[cy]` 与 `H2O2/water[pe]` 混合区室表示保持未决，不因方向修复称为正确定位。其 GPR `YALI1E40671g`——正式原生名称未核实——缓存为 Catalase T 候选（UniProt A0A1D8NL70/AOW06378.1，entry 34/sequence 1，更新 2026-01-28，existence 3），功能与 **cytoplasm** 定位来自自动/同源注释。当前精确反应是 model/GPR assignment，非 W29 乙醇氧化或过氧化物酶体定位实验。

**OAADCm 的 H+ 去除有明确中性物种依据。** 当前 m44 为 CO2/0，m6 为 C3H4O3/0 的中性 pyruvic acid，m77 为 C4H4O5/0 的中性 oxaloacetic acid；后两者 XML 明确保留原 formula/charge。它们与 [ChEBI:16526](https://www.ebi.ac.uk/chebi/CHEBI:16526)、[ChEBI:32816](https://www.ebi.ac.uk/chebi/CHEBI:32816)、[ChEBI:30744](https://www.ebi.ac.uk/chebi/CHEBI:30744) 的形式一致，因此 `pyruvic acid + CO2 → oxaloacetic acid` 已守恒，不再产 H+。原式 H/charge 均 +1；去掉产物 m28 后全部元素和电荷为零。带电形式的 `oxaloacetate²⁻ + H+ → pyruvate⁻ + CO2` 不能直接移植到模型保留的中性形式。metadata 的 before 原本不含 H，after 在这些中性物种不变时补入 H；这不是任意调一个质子使见证消失。

其负向是 OAA 脱羧；[IUBMB EC1.1.1.38](https://iubmb.qmul.ac.uk/enzyme/EC1/1/1/38.html) 的第 2 个反应正是 OAA→pyruvate+CO2，故不能因无 NAD 就把这个 EC 认定为错配。[KEGG R00217](https://www.kegg.jp/entry/R00217) 亦给出相同身份及此 EC。[Sender 等 2004 原始摘要](https://pubmed.ncbi.nlm.nih.gov/15251467/) 在乳酸乳球菌纯化 CitM 中测得 OAA 脱羧，但没有 malic activity，进一步提醒不能仅凭家族名接受全部酶活。

当前正向没有 ATP 或离子势耦联。[ATP 耦联 pyruvate carboxylase](https://iubmb.qmul.ac.uk/enzyme/EC6/4/1/1.html) 及 [Na+ 耦联 OAA decarboxylase](https://iubmb.qmul.ac.uk/enzyme/EC7/2/4/2.html) 是带额外计量机制的不同反应，不能作为这条无耦联正向的支持。因此将上界置 0、保留脱羧，是可检验的条件性方向候选；未取得 W29 代谢物浓度和 ΔG，不声称任何条件绝对不可逆。数据库排版的等号/双箭头本身也不能证明生理双向。

OAADCm 的 `YALI1E22303g`——正式原生名称未核实——缓存为 malate dehydrogenase (oxaloacetate-decarboxylating) 候选（UniProt A0A1H6QCK5/AOW05618.1，entry 34/sequence 1，更新 2026-01-28，existence 3）；缓存的两类催化反应均自动/同源注释。这次保留 GPR，不作新的蛋白、底物谱或区室验证。

**仍存在潜在化学问题，不能把能量见证消失称全网化学通过。** 按固定 XML 逐列计算，R305 当前 H/charge 残差为 −2/−2，R1889 为 −1/−1；提供见证两列均通量 +0.769230769，与 OAADCm +2.307692308 所产的错误 H/charge 恰好抵消。这只是解释该见证的记账机制，不是允许残差互相抵消的理由。既有 CoQ 审批门保留，此处仅记录残差，不裁决或改动这两列的质子、电子或泵耦联计量。也未因正常 R533/R538 或呼吸链列参与而建议删除它们。
