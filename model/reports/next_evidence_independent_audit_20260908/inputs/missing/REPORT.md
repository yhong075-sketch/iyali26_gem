# 缺失代谢功能线索：20项排序、5项精查

已对既有1282条 `identified_absent_from_model` 记录完成有界注释筛选，保留20条优先线索并精查最终前五。结果支持一个钙／锰稳态与运输的表示缺口，以及四个复合体I成员的基因关联核查方向；**目前没有可直接加入模型的新反应，也没有可直接激活的GPR补丁**。本报告的24条决定性主张均待总协调安排独立来源审计。

| 排序与身份 | 已支持的功能层级 | 对参考模型的意义 | 当前限制 |
|---|---|---|---|
| 1．YALI1E11674g — PMR1/SCA1 — 钙ATPase候选，原生干预表型已报道 | 钙／锰稳态及分泌加工相关作用 | 未检出显式Ca/Mn代谢物；Golgi质子泵不能代替该离子运输功能 | 直接运输计量、原生定位测量及SD-Leu条件仍待核实 |
| 2．YALI1A03153g — NUZM — 复合体I附属亚基，成员身份有实验支持 | 原生复合体成员 | 已有NADH/Q9反应，需审查基因到复合体的关联 | 单亚基缺失的残余活性未建立 |
| 3．YALI1A18188g — NIDM — 复合体I附属亚基，成员身份有实验支持 | 原生复合体成员 | 同上 | 旧论文明确保留具体定位的不确定性 |
| 4．YALI1C30086g — NIMM — 复合体I小型附属亚基，成员身份有实验支持 | 原生复合体成员 | 同上，另有旧质量注释歧义 | 87 aa记录不能按旧“79 kDa”字符串认作大型催化亚基 |
| 5．YALI1E14698g — NUJM — 复合体I膜臂附属亚基，成员身份有实验支持 | 原生复合体成员及标记分离支持 | 同上 | 成员追踪实验不等于该基因删除实验 |

PMR1的生长补救依据来自原生基因破坏研究；分泌研究提供另外一种功能表型，但来自重叠研究团队和突变株谱系，不当作同一生长结论的独立重复。四项复合体I线索以原生结构沉积与对应蛋白记录为主，并核查了生化分离研究；数据库互引不增加独立实验数量。[Park等，1998](https://pubmed.ncbi.nlm.nih.gov/9461422/)、[Sohn等，1998](https://pmc.ncbi.nlm.nih.gov/articles/PMC107781/)、[Grba与Hirst，2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7612091/)、[Angerer等，2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3273332/)。

## 参考模型实际表示了什么

只读暂定科学参考 `iyali26_gem_integration/model.xml`，核对SHA-256为 `bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee`。它包含1074个基因、2313个反应、1877个代谢物。未以本工作树默认模型替代，未应用培养基或菌株覆盖层。

下式忠实转写文件中的系数，**不表示本轮已完成电荷／微物种化学审计**。`mi`为模型线粒体区室，`cy`为胞质，`Q9`为ubiquinone-9。

| 反应 | 储存的计量 | GPR／界限 |
|---|---|---|
| R1889，complex I | NADH_mi + Q9_mi + 5 H_mi → NAD_mi + Q9H2_mi + 4 H_cy | GPR空；0…1000 |
| R2062，NADH:ubiquinone oxidoreductase | NADH_mi + Q9_mi → NAD_mi + Q9H2_mi | 现有多基因AND表达式；0…1000 |
| R573，同名氧化还原反应 | 与R2062相同 | 现有单基因表达式；**0…0，已禁用** |
| R570，跨区室NADH dehydrogenase | NADH_cy + H_cy + Q9_mi → NAD_cy + Q9H2_mi | NADH使用区室不同；0…1000 |

这些结果说明复合体I的总反应并非完全缺失。给R2062增加基因关联，不会自动给R1889建立遗传约束；然而，空GPR或重复计量本身都不能证明特定条件下存在可行旁路。R573已禁用，尤其不能当作活跃替代路径。完整表达式、代谢物ID和界限见[静态快照](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/model_static_snapshot.json)。

对于PMR1，所有代谢物ID、名称和分子式元素词检索均未命中Ca/Mn。R794是已有Golgi质子V-ATPase，计量不含Ca/Mn。由此可以提出“显式离子运输未表示”的候选缺口；不能由此推定一个已验证的离子／ATP计量，也不能给现有蛋白分泌过程任意附加钙需求。[检索定义与结果](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/calcium_manganese_static_search.json)。

从约束建模角度，新增未连接到需求或稳态约束的运输步骤未必改变生长预测；在已有反应上补成员名称也未必产生正确的删除表型。应先判别缺的是反应、区室、辅因子需求还是基因关联。本轮没有通过调边界、加需求或改AND规则追求召回提升。

## 范围与排序解释

排序采用字典序：原生功能证据层级E → 具体步骤关联G → 身份确定程度I → 系统ID升序。E3是原生干预表型报道，E2是原生复合体成员证据，E1为人工间接注释，E0为自动／旧注释；不把四种证据混成加权分数。G2包括可命名的既有反应关联，**不等于已证实反应缺失**。I2限于S2单一对应与官方accession桥接，不等于实验菌株间功能等价。

20项中，9项获得原生复合体成员层级支持，1项有原生基因干预表型报道，另10项仍停留在人工间接或自动注释层级。原生单基因直接催化／运输实验在本轮确认数为0；不得将这20项称为“20个确认缺失酶”。其余15项只完成浅查，不具备五个精查案同等的文献覆盖。

1282条中有131条被宽关键词检索命中，关键词只帮助找线索。未入选1262条、未精查1277条都没有被判为无代谢功能；这也不是证明20条是全体1282条的全局最优排序。历史旧标签完整保留：本次入选包括旧 `non_metabolic` 1条、`no_data` 11条、`metabolic` 8条。另8条 `mapping_conflict` 单独保存在[歧义记录](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/ambiguities8.json)，未用于功能定案。

历史覆盖 **322/1612** 和交集内召回 **67/322** 原样保留；本轮没有新做essentiality评估，不产生新的召回结果。

## 交接与最小后续工作

1. 总协调先独立核查[24条主张](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/claims.jsonl)与[来源登记](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/source_registry.json)，尤其PMR1物种／条件、四个结构成员的身份桥接及R573禁用状态。当前审计覆盖为 **0/24**，本生产者的来源核查不充当独立审计。
2. PMR1优先补全1998年Gene原文、目标菌株位点对应，以及原生膜运输／定位与Ca、Mn、ATP计量证据。只有在确认模型需要表达的具体过程后，才设计可审查的反应及需求连接方案。
3. 四项复合体I先统一核查R1889／R2062的生物学分工与基因规则依据，查找各目标亚基的原生删除／补救、复合体装配和残余活性数据。成员身份不能直接生成全体AND规则。若需新实验，最小判别是目标基因干预与回补配合装配／活性读数，并匹配目标培养条件；本轮未执行这些实验。

[完整20项排序](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/RANKED20.md) · [PMR1档案](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/01_YALI1E11674g_PMR1.md) · [NUZM档案](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/02_YALI1A03153g_NUZM.md) · [NIDM档案](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/03_YALI1A18188g_NIDM.md) · [NIMM档案](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/04_YALI1C30086g_NIMM.md) · [NUJM档案](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/05_YALI1E14698g_NUJM.md)。

检索边界、技能实际用法、未完成检查及时间记录见[过程记录](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/SEARCH_BOUNDARY.md)。输入SHA、来源行、序列版本及验证结果均已留存。未改模型、GPR、原标签或共享状态，未求解、提交或推送。
