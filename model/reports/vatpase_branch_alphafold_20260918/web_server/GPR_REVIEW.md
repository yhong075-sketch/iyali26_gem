# 原 GPR 的独立生物学与布尔逻辑审查

初查日期：2026-09-18；模型身份更正核查：2026-09-21。范围为原始“前 12 成员 AND”与“后 2 成员 AND”之间的 OR、已保存的“三个共同 AND + 原 OR”假设、未实施的全 14 AND 逻辑对照，以及“共同组件 AND 两种 a 候选择一”的候选规则。初查已重新打开固定序列来源、原始论文全文／摘要及 W29 冲突位点记录；2026-09-21 重新读取 curation、两份指定模型的实际 XML GPR，以及支持共享组件和 a 亚型的原始论文快照。没有修改模型、运行求解器或新增预测。用户提供的官方 TypeSafe 技能已于初查直接读取，并据其方法拆分证据判断；没有调用 TypeSafe/Jev API，也不编造概率。本审核未读取新 AlphaFold 完整结果包，其结果判断由根任务另行报告。

**更正实施身份：** 本文初稿把保存的假设误称为“已有／获准运行的全 14 AND”，该表述错误。`vatpase_gpr_hypothesis.json` 的 `required_genes` 只有 D00581、旧 E16192、F38820，`after_gpr` 与 `model_metadata_trna_vatpase_and_hypothesis.xml` 均是在原 OR 外加这三个共同 AND。另行读取的 `model_metadata_trna_r539_alphafold_labeled.xml` 则仍为原始 OR。以下按具体文件分别说明；全 14 AND 仅为未实施的逻辑对照，不能与既有构建／screen 结果绑定。

**原外层 OR 缺乏机制支持；全 14 AND 也不能认定为正确生物学。更符合已知装配机制的待验证形式是“共同组件 AND（a₁ OR a₂）”，其中 c″ 应归入共同组件候选。该 OR 仍需满足同一反应、同一区室的完整装配体可替代条件，不能直接同时填进 Golgi 和液泡反应。**

## TypeSafe 方法的实际应用

本轮直接读取 [官方 SKILL.md](https://raw.githubusercontent.com/typesafe-ai/skills/main/skills/typesafe-ai/SKILL.md) 及其 [State](https://docs.typesafe.ai/concepts/state)、[Choice](https://docs.typesafe.ai/primitives/choice)、[Confidence](https://docs.typesafe.ai/confidence) 实时文档。应用的是“明确证据状态、逐项提问、保留未知与不匹配结果”的审查组织方法，不把软件输出类型保证当作科学事实保证，也不把任何置信度分数当成修改 GPR 的授权。

证据状态分为 `identities`（固定序列／菌株）、`reference_mechanism`（原始实验来源）、`reaction_chemistry`（实际反应及区室）、`gpr_logic`（给定布尔式）和 `unresolved`（原生定位、共同依赖及 F 冲突）。末尾 9 个问题独立评价，`supported` 指来源或布尔演算支持该限定主张；`contradicted` 指明确反证；`unresolved` 指现有证据不能定案；`no_match` 用于所列候选都不满足目标。本轮不强迫从四种 GPR 中选出“已验证正确”的一个，正式规则的选择仍为未定。

可观察事实包括输入序列、原文实验、源记录伪基因注释和布尔式；“Yarrowia 可能采用共同组件加 a 亚型”的建议是基于这些事实的跨物种机制推断。两者在下面分别标注，没有将推断回写为观察结果。

## 位点身份与功能边界

下面所有 Yarrowia 位点的原生正式名称均未核实；亚基字母是功能类别候选，不是正式基因符号。证据为未审阅数据库同源注释，部分有此前 BLAST／AlphaFold 预测支持，**均不等于已实验确认的 W29 原生功能**。14 条输入的序列、版本、菌株和源文件身份沿用已独立核实的 [sequence_manifest.json](../sequence_manifest.json) 和 [DESIGN_REVIEW.md](../DESIGN_REVIEW.md)。

| 系统 ID | 原生正式名 | 蛋白功能候选 | 本候选规则中的角色 |
|---|---|---|---|
| YALI1A11258g | 未核实 | V1 C，复合体结构连接 | 共同组件候选 |
| YALI1E14125g | 未核实 | V1 G，外周柄 | 共同组件候选 |
| YALI1D00581g | 未核实 | V1 D，中央转轴；既有序列／预测结构支持 | 共同组件候选 |
| YALI1E10492g | 未核实 | V1 H，调节及复合体稳定 | 共同组件候选 |
| YALI1A09766g | 未核实 | V1 A，ATP 水解催化 | 共同组件候选 |
| YALI1F20965g | 未核实 | V1 E，外周柄 | 共同组件候选 |
| YALI0E16192g | 未核实 | V1 F，中央转轴；本次蛋白来自 CLIB122 | 共同 F 功能候选；W29 基因赋值未解决 |
| YALI1F31854g | 未核实 | V0 c，膜内转子环 | 共同组件候选 |
| YALI1F21690g | 未核实 | V0 d，转子连接；与 V1 D 不同 | 共同组件候选 |
| YALI1E37063g | 未核实 | V0 c′，膜内转子环 | 共同组件候选 |
| YALI1E32332g | 未核实 | V1 B，核苷酸结合及催化头结构 | 共同组件候选 |
| YALI1E12482g | 未核实 | V0 a，质子通道与装配连接；区室未定 | a₁ 候选 |
| YALI1F38820g | 未核实 | 偏 Vph1-like 的 V0 a 候选，仍需实验确认；区室未定 | a₂ 候选 |
| YALI1F13017g | 未核实 | V0 c″，膜内转子环；既有序列／预测结构支持 | 共同组件候选 |

## 四种 GPR 的判断

令 `P` 为原前 12 个基因的 AND，`a₁=YALI1E12482g`、`a₂=YALI1F38820g`、`c″=YALI1F13017g`；`D=YALI1D00581g`、`F=YALI0E16192g`，其身份和功能限度见上表。

| 规则 | 数学含义 | 生物学判断 |
|---|---|---|
| 原始 `P OR (a₂ AND c″)`；r539 标记版仍保存此式 | 只要 a₂ 与 c″ 在场，前组所有 V1 组件均可缺失，反应仍可用；c″ 缺失时 P 又能独立维持反应 | 与已查 ATP 耦联质子泵机制不符，不能把后两个膜亚基视为另一台完整泵 |
| 已保存假设 `D AND F AND a₂ AND (P OR (a₂ AND c″))` | 三个指定基因任一缺失都会关闭 GPR；但 D、F、a₂、c″ 在场就仍可绕过其余催化头／共同组件 | 与 curation 及 vatpase 假设版 XML 一致；完成了用户指定的三个共同依赖，仍不是完整复合体规则修复 |
| 全 14 AND，未实施逻辑对照 | 同时要求两种 a 候选与全部其余组件 | 没有证据证明同一泵或同一区室反应同时依赖两个 a 位点；不是已保存假设，也无对应本轮构建／screen 结果 |
| `C AND (a₁ OR a₂)` | 两种装配形式共享 C，每种形式各含一个 a 候选；`C=(P 去掉 a₁) AND c″` | 更符合参考复合体组织的候选形式；原生共同依赖、两种 a 的定位及同区室替代能力仍待验证 |

这里 **C 一共 12 个基因**：原前组去掉 a₁ 后是 11 个，再补入 c″。它不是把原前 12 原封不动保留后再套一个 OR；若 C 已含 a₁，则 `C AND (a₁ OR a₂)` 会吸收掉 a₂ 的作用，仍然强制 a₁。

完整候选式为：

```text
(YALI1A11258g and YALI1E14125g and YALI1D00581g
 and YALI1E10492g and YALI1A09766g and YALI1F20965g
 and YALI0E16192g and YALI1F31854g and YALI1F21690g
 and YALI1E37063g and YALI1E32332g and YALI1F13017g)
and (YALI1E12482g or YALI1F38820g)
```

这只是**在所给 14 位点范围内的待验证装配关系候选**，不声称已列全原生成熟泵的所有必需组件；参考体系其他小亚基／装配因素的对应并未在本轮补齐。本文件没有接受它为正式 GPR。布尔分配律给出 `C AND (a₁ OR a₂) = (C AND a₁) OR (C AND a₂)`。要赋予这个 OR 生物学含义，需要**两种含共同组件的完整装配形式**可完成同一步骤；不表示 a₁ 或 a₂ 单独就能催化 ATP 水解和跨膜泵质子。

原规则的错误传导可用一个例子看清：关闭 ATP 水解催化亚基候选 YALI1A09766g 时，P 为假，但 `a₂ AND c″` 仍真，于是原 GPR 仍为真。**已保存的三个共同 AND 假设也保留这个反例：** 只要 D、F、a₂、c″ 仍存在，其外加三个条件满足，内层后分支又为真。它修复的是三个指定基因的关闭传导，没有消除其他核心亚基被原 OR 绕过的问题。此为直接布尔结论，无需生长求解。全 14 AND 逻辑对照则预设了两个 a 都不可缺少；与已有实现须分开。

## 原始来源支持到哪里

1. **完整泵依赖 V1 和 V0 的共同工作。** Vasanthakumar 等的原始结构／生化研究区分 ATP 水解催化头和膜内质子转运部分；比较的是分别包含两种 a 亚型的完整复合体。本文直接重读本地封存原文的 Abstract、Introduction、ATPase Activity Assays 和 Overall Structures。它支持“共享组件 + a 类型变化”的参考机制，不支持两个膜亚基独立替代催化头。[PNAS 2019，DOI 10.1073/pnas.1814818116](https://pmc.ncbi.nlm.nih.gov/articles/PMC6462096/)

2. **c、c′、c″ 不能因同属膜蛋白就当作互相替代。** Hirata 等在酿酒酵母对膜环亚基的基因破坏和位点突变实验显示，这三种蛋白脂质亚基具有相关但非冗余的作用。这支持将本题 c″ 放回共同膜环依赖的调查方向，尚非对 Yarrowia 三位点的直接敲除实证。本轮直接重读 PMID 9030535 的原始摘要；未称重读其未取得的全文。[Hirata 等，JBC 1997，PMID 9030535](https://pubmed.ncbi.nlm.nih.gov/9030535/)

3. **相同 a 类别不等于同区室可替代。** 酿酒酵母 YOR270C／VPH1（V0 a，主要液泡）和 YMR054W／STV1（V0 a，主要 Golgi／内体）的功能、定位及完整泵比较有实验依据；2019 原文还指出删除原生 a 基因并表达单个亚型的实验条件可能使 Stv1 复合体错定位至液泡。这类条件性结果不能直接证明野生型定位下可交换。本题两个 Yarrowia a 位点不能仅靠序列/结构相似度分别冠以已确认 VPH1/STV1 身份。[PNAS 2019，Introduction 与 ATPase Activity Assays](https://pmc.ncbi.nlm.nih.gov/articles/PMC6462096/)

本轮网页访问 PMC/PubMed 动态正文未返回可读论文；以上判断来自重新打开的本地原始全文／摘要快照，未以检索结果标题或旧总结代替原文。

## 区室、菌株与 AlphaFold 的限制

本轮只读检查实际候选模型 XML 的反应计量：R794 将胞质底物中的 ATP／H₂O 与 H⁺ 关联到 **Golgi H⁺**，R795 对应 **液泡 H⁺**；它们产物区室不同，不能仅以“都是 V-ATPase”合并判断。若将来确定 a₁/a₂ 分别负责不同区室，应按定位将相应反应写为 `C AND 对应该区室的 a`。只有证实两种完整装配形式在该区室、该条件下都能工作，才可对那条具体反应保留 a₁ OR a₂。当前不足以指定哪一个必定属于 Golgi，也不足以把 F38820 的 Vph1-like 候选标记改成确定液泡定位。

**W29 的 F 位点仍是阻碍正式基因赋值的冲突。** 本次 A 中 YALI0E16192g 的 Q6C5Q2 v1／CAG79603.1 是 CLIB122 完整蛋白；W29 对应 YALI1E19360g（正式名未核实，旧 F 候选的注释对应位点）在 NC_090774.1 来源记录中有 `pseudo` 和移码导致非功能的注释。本轮直接重新打开该原始记录。不能把 CLIB122 链的良好预测自动变成 W29 活性蛋白证据，也不能由数据库伪基因注释断言 W29 没有其他 F 功能来源。正式 GPR 需要解决实际菌株序列／转录本或替代 F 位点。[NCBI Gene 2911580](https://www.ncbi.nlm.nih.gov/gene/2911580)

AlphaFold 可支持局部折叠和装配候选，不能单独证实催化耦联、亚细胞定位、两个候选的原生互补或共同必需性。即使 B 的 a/c″ 接口预测很好，也没有补出 ATP 水解催化头；界面低置信则只能说明本次装配未解决，不能证明真实不结合。输入每种一拷贝的局限影响结构装配解释；GPR 本身是基因依赖逻辑，并不编码 A₃/B₃ 等蛋白拷贝数。

## 验收与独立审查覆盖

接受候选规则前至少需要：固定目标菌株和 F 功能基因；确认 a₁/a₂ 的区室及是否在对应条件下可替代；确认共同组件依赖。获准变更后再检查 GPR 关闭反应的传导、同区室重复/旁路及相同条件下的既有正确结果。即使这些层面通过，泵反应是否影响生长仍是另一层需求约束问题，不能用 AND/OR 直接保证 essential。

| 原子主张 | 审计裁决及边界 |
|---|---|
| B 的两个候选为 a/c″，输入不含 V1 ATP 催化头 | supported；固定身份为候选注释，不升级原生功能 |
| 参考完整泵由共同组件加 a 亚型形成 | supported；直接原始结构／生化来源，限参考体系 |
| 参考 c/c′/c″ 不可仅因同源就视为冗余替代 | supported；直接原始实验摘要，Yarrowia 外推保留 |
| 原 OR 可让 V1 催化亚基缺失时 GPR 仍为真 | supported；静态布尔逻辑 |
| 已保存假设为三个共同 AND + 原 OR，且仍允许催化 A 缺失；另行 r539 标记版保留原 OR | supported；2026-09-21 直接核对 curation 与两份 XML GPR、分析反例 |
| 全 14 AND 是 W29 的真实共同依赖 | unresolved；两个 a 同时必需及 F 身份未确认 |
| 共同 12 AND(a₁ OR a₂) 与两套共享 C 的完整分支逻辑等价 | supported；布尔恒等式，不等于生物学有效性 |
| 两种 Yarrowia a 可在 R794/R795 各自区室互相替代 | unresolved；缺原生定位、同条件互补证据 |
| CLIB122 旧 F 不能由本次预测自动升级为 W29 活性 F | supported；版本化原序列及 W29 pseudo 来源直接核实 |

**9 total | 9 audited | 7 supported | 2 unresolved | 0 contradicted | 0 unchecked**。这是原 8 项加本次实施身份核查 1 项的限定主张分母；已明确更正初稿中“保存假设等于全 14 AND”的错误，不将该错误从更正记录中抹去。这些来源／逻辑核查不合并其他旧审核分母，也不把原始 OR 的反例与跨物种推断当成 W29 新实验。

后台来源身份：

| 直接读取对象 | 完整 SHA-256 |
|---|---|
| `sources/PMC6462096.html`（历史审查原始全文快照） | `bc670498ba7be831a46b5b7db3ed6e307f1c7dbec09aff6f55ba9d6eacd55f3d` |
| `sources/pubmed_9030535_25546637_8385405.xml`（采用 PMID 9030535 原始摘要） | `ac37914432b828606e76b0c351271db94f38ee94d65b24fd169303762dcef5ac` |
| `sequence_manifest.json` | `e0b82204f30d60c64b19ba2724134c371114e3780adfbff2c316f2445a760118` |
| W29 `ncbi_locus_E19360.gb` | `e6acceec7082b635090b6b0965e7b904202efc8f48bba08c6cfecc4248bcbdc6` |
| 实际只读反应计量来源 `model_metadata_trna_vatpase_and_hypothesis.xml` | `b00a9ea20c21712b428f159f9050727053cf373ae89c6cfe161e7a77e0b2eb64` |
| 2026-09-21 核对的 `data/reference_build/curation/vatpase_gpr_hypothesis.json` | `b41e5cc151c575dd9f31c4f454b55f60c1c053b9e0f69915d9844f3ece200e6b` |
| 2026-09-21 核对的 `model_metadata_trna_r539_alphafold_labeled.xml` | `cff96f40368b3de3d502956cf33dda246e049becfa43c1bba12fbd8ded6acd9a` |
