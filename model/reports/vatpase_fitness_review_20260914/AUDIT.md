# 三个 V-ATPase 候选 fitness：独立来源审计

核验日期：2026-09-14；审阅者 `/root/r795_open_audit`。本轮直接读取原工作簿、CSV、指定映射原行、论文方法、固定模型 XML 与既有结果，终审 REPORT.md；新增 LP 0，科学文件修改 0，只写本文件。以下是来源核验与数学推导，不称本次重跑实验或求解。

**结论：两个有记录目标在 Cas9/Cas12a 原表均明确为 Non-essential；E16192 及其指定 W29 候选 ID 均缺失。负 FS 不能直接换算质子泵容量。无 pool 的零泵最优见证证明，仅降低泵的非负上限不会改变已有 WT/KO 生长；pool 条件下上限低于 1 才可能限制额外收益，中间点未求解。**

## 原始记录

两 sheet 表头均为 A `Gene ID`、B `FS`、C `Raw p-value`、D `Corrected p-value`、E `Essentiality`。下表保留原始 ID；最后一列四条均原文 `Non-essential`。

| 来源行 | 原始 ID | FS | Raw p-value | Corrected p-value | 原 call |
|---|---|---:|---:|---:|---|
| Cas9!A2242:E2242 | YALI1_D00581g | -2.62728401307139 | 0.0488135489057504 | 0.171028417093779 | Non-essential |
| Cas12a!A1525:E1525 | YALI1_D00581g | -1.78397035253357 | 0.0154803057135094 | 0.079179122727563 | Non-essential |
| Cas9!A4168:E4168 | YALI1_F38820g | -1.74867415475292 | 0.558012311238163 | 1 | Non-essential |
| Cas12a!A2251:E2251 | YALI1_F38820g | -1.51569878949371 | 0.0902039895359398 | 0.312506710414511 | Non-essential |

独立完整扫描 Cas9 的 7,854 条数据行、Cas12a 的 7,795 条数据行和 CSV 的 15,649 条数据行；仅去下划线作格式归一化，各有记录目标每 assay 恰一条。`YALI0E16192g` 与 `YALI1E19360g` 均为 0 条匹配，未用缺失补零或补 nonessential。CSV 四条记录的 FS/p/q/call、原始 ID、sheet、source_row 和 source_sha256 与工作簿一致。

映射工作簿 `YALI1 Genes!A6015:H6015` 直接核实 `YALI0E16192g → YALI1_E19360g`，C6015 为 `Pseudo mRNA`，H6015 为 V-ATPase F 亚基同源注释。该行支持追查候选 ID，不支持完整活性蛋白或 assay 分数的存在；未作全映射重审。

身份限定：YALI1D00581g 为 V1 D 中央转轴亚基候选；YALI0E16192g／指定 W29 候选 YALI1E19360g 为 V1 F 亚基候选且有编码注释冲突；YALI1F38820g 为偏 Vph1-like 的 V0 a 质子转运/装配亚基候选。正式原生基因名均未核实，既有序列/AlphaFold 支持属间接证据，本轮不重复蛋白功能鉴定。

## 分数与容量解释

直接核读 [acCRISPR 原论文](https://www.nature.com/articles/s42003-023-04996-8) Methods “acCRISPR framework” 式 2/3 与 essential-call 方法：guide FS 是处理/对照归一化计数比的 log2，gene FS 是经活性筛选后的 guide FS 平均；单尾 z 检验后按 FDR 校正 p<0.05 调用 essential。四条校正 p 均不达此阈值；未校正 p 不能替换调用标准。FS 描述筛选丰度变化，不是测得酶活、反应容量或 FBA 生长率；`2^FS` 也无据直接作为剩余酶容量系数。原实验 PO1f/SD-Leu 条件在论文培养方法中可核读，但不能推出所有环境下的非必需性。

固定 XML 全部 GPR 引用扫描确认三个指定模型目标仅关联 R794/R795，均处共同 AND。原 R794=[0,1000]、R795=[0,0]；既有 `WT_fluxes.tsv` 的 2,315 行中两泵均为 0，`diagnosis.json` 保存 optimal 生长 1.8718823069403；双泵关闭保存 1.8718823069402994，差在浮点容差内。

令旧可行域 F、新非负泵上限形成 F′，则 F′⊆F，同时旧最优零泵见证 v*∈F′，故 `max(F′)=max(F)`。三个共同 AND KO 原已关两泵，再降低 WT 的同两条上限不改变该 KO 问题；固定培养、目标及其余约束下，原 KO/WT 不因这一步降低。此证明不适用于增加正下界、其他反应、培养改变或不同目标。

pool 部分复核既有原始 results、有效反应及此前本审阅者已逐行审计的完整计量：`2 v_R795=v_R2034+v_R2039`，两个有效二肽源各上限 1 使水解分别在 [0,1]。因此全可行域内 R795≤1；仅把其上限 1000 改为 100、10 或 1 不排除任何解，低于 1 才可能限制额外收益。旧 pool-WT=2.1722503622879215，R795-off 最大生长=1.8718823069403，比值 0.8617249371608993。报告没有宣称已得到中间容量曲线；也未把人工源、反应关闭或通量可行性当原生修复/基因必需性验证。

## 身份、覆盖及限制

本次重新计算并匹配 manifest 的 SHA256：

- 原始工作簿：`ba1eca8f0c0c31b500388a796fb7e6de51769e7086f8864cb2c7b40da7c6c0c5`。
- assay CSV：`97ef559651cd99ce63144ffe08ae64e983ecaec60d568f859da8926c14c7d9ff`。
- 映射工作簿：`42658fb95202b4114e5c6b5dbd60020de7a4038797ee2211a5352078a729a670`。
- 假设模型：`b00a9ea20c21712b428f159f9050727053cf373ae89c6cfe161e7a77e0b2eb64`，匹配原诊断的实际输入身份。

| ID | 独立核查声明 | 判定 |
|---|---|---|
| VF1 | 四个 FS/p/q/call 的原工作簿来源与 CSV 对应 | supported |
| VF2 | 完整 assay ID 扫描中旧 E16192/候选 E19360 缺失 | supported |
| VF3 | 指定映射原行连接两 ID 且标 Pseudo mRNA | supported |
| VF4 | FS 定义、FDR 阈值及不能直接换算容量的限制 | supported |
| VF5 | 固定模型、共同 AND 关联与旧零泵最优见证 | supported |
| VF6 | 无 pool 仅降低非负泵上限不改变原 WT/KO 的可行域证明 | supported |
| VF7 | pool 完整计量给出 R795≤1，故上限≥1 的收紧无效 | supported |
| VF8 | pool 开/关端点数值与未求解中间曲线的范围声明 | supported |

**total claims 8 | audited 8 | supported 8 | unresolved 0 | contradicted 0 | unchecked 0**，分母仅为此表。REPORT 中其余映射行、YALI0-only 页、1612 正例清单及转座子 calls 未在本轮独立重查，不包含在通过数内。未建立本地原工作簿与出版商附录的字节一致性，未审逐 guide 数据/实验效率、完整历史运行环境或真实残余活性；这些限制不由四行核验解除。
