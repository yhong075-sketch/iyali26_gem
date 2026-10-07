# 内源供给模块规格：默认禁用，尚未实施生物学修复

机器可读配置为 `run/endogenous_module_spec.json`；`validate_endogenous_spec` 在禁用时不修改模型，开启且缺参时报出字段清单。填满字段本身也不表示科学审查通过，当前接口明确拒绝激活。没有新增 protein_pool、内源 source、反向 biomass 或反向二肽水解。

## 现有模块能复用什么

当前模型有20条 `TRNA_BIOMASS_*`，使用氨酰tRNA并归还未装载tRNA，生成 biomass 消耗的私有蛋白残基。这提供净生长的残基账本与上游充电成本；不提供成熟蛋白的序列、丰度或周转池。私有残基没有完整分子式，不能直接当作已配平的成熟多肽。`protein_module_inventory.json` 保留完整计量与notes。

关键二肽节点的全部邻接没有内源生成。关键词检索命中的蛋白甘露糖基化、甘氨酸裂解复合体和糖原自噬不是目标蛋白周转模块。现有 `R859` 的糖原自噬不能改名为蛋白自噬。

## 必须提供的输入与单位

| 配置项 | 要求 | 当前状态 |
|---|---|---|
| precursor_identity_and_sequence | 菌株明确的成熟链/长肽身份、处理后序列及版本；膜/分泌加工另外记录 | 未选定真实前体 |
| precursor_formula_charge | 前体及末端/修饰、pH表示、物种分子式和电荷 | 缺失 |
| residue_composition | 每条链每种残基的摩尔数 n(a,j)，可由上述序列计算 | 缺失 |
| source_synthesis_reactions | 从已供应底物到前体的完整合成路径；禁止无底物source与biomass逆转 | 未建立 |
| synthesis_energy_cost | 翻译/聚合、运输、加工及替换合成成本；区分现有GAM已包括部分 | 未确定，不能只复制23.09 GAM |
| degradation_rate_mmol_precursor_gdw_h | 每条链的总降解通量，mmol chain/gDW/h，或由丰度P和k推导 | 未测定 |
| vacuolar_fraction | 总周转进入液泡的比例，0–1 | 未测定 |
| nonoverlapping_dipeptide_yields | 每条前体实际释放GD/GE/AG/GP的产率；不得重叠计数 | 未测定 |
| remaining_products | 未进入四种二肽的全部残基、剩余肽段及去向 | 缺失 |
| balanced_degradation_equation | 水、质子、电荷及全部产物可核的方程 | 未建立 |
| localization_and_transport_evidence | 各步原生区室及跨膜机制；已有运输的耦联量不能任意改 | 不足 |
| net_growth_vs_replacement_accounting | 净生长与替换合成分账，不重复GAM/充电成本 | 仅有净生长一侧 |
| condition_matched_parameter_sources | PO1f/SD-Leu等适用条件、实测/跨条件/纯情景标签 | 缺失 |

## 残基与质量账本

对前体 j，令 P_j 为 mmol chain/gDW，k_j 为 h⁻¹，v_deg,j = k_j P_j；若稳态单位生物量蛋白丰度恒定，则 v_syn,j = μP_j + v_deg,j。液泡部分 v_deg,va,j = f_va,j v_deg,j，其余 (1−f_va,j)v_deg,j 走其他路径；不强制进入目标二肽。

令 y_i,j 为每条链释放二肽 i 的摩尔数，则 q_i = Σ_j y_i,j v_deg,va,j。每个前体要同时满足残基预算和非重叠切割要求，不能给四种二肽各自独立的同一份Gly预算。总体至少满足：

```text
q_GD + q_GE + q_AG + q_GP <= r_Gly
q_GD <= r_Asp; q_GE <= r_Glu; q_AG <= r_Ala; q_GP <= r_Pro
r_a = Σ_j n(a,j) v_deg,va,j
```

对于无修饰、线性、L残基的单条中性肽，若每链释放 y=Σ_i y_i 个目标二肽，其余全部变成游离氨基酸，则计数层面的示意为：

```text
precursor_j + (L_j − 1 − y) H2O
  → Σ_i y_i dipeptide_i + Σ_a [n(a,j) − Σ_i count(a,i)y_i] amino_acid_a
```

这是限定化学假设下的账本模板，不是本模型已验证方程。真实修饰、末端、酸碱形式及中间肽需要另列；水解目标二肽还会再消耗 y 个水。所有余项非负、所有残基与端基守恒，且任何前体生产都必须付出完整合成成本。可用序列统计限制潜在上限，但相邻GD/GE/AG/GP存在不代表实际释放；不能将重叠片段同时分配。TPM不等于降解通量。

## 最小证据与验收

1. 确认目标培养下四种二肽的化学身份（分开Ala-Gly/Gly-Ala等顺序异构体）、含量及时间变化；含量不能直接当flux。
2. 用标记蛋白/氨基酸的脉冲追踪区分新合成、蛋白回收及外源供给，并估计前体丰度、降解率与液泡分配。
3. 验证上游长肽释放酶和下游游离二肽水解酶的确切底物及定位；两个步骤分别提供证据。
4. 用已测参数构建完整守恒、包含残余物和替换成本的候选副本，重复闭合碳/氮、能量/循环和本次诊断；与原模型同条件回归。

目前不能执行可信的“周转速率敏感性”：基础前体、组成、成本与速率都缺失。本轮仅执行用户指定的人工二肽输入上限敏感性，不将其称为蛋白周转敏感性或内源修复。正式可接纳的内源反应数为0。
