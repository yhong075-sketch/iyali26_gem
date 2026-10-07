# 15% 阈值下的 35 个未收录基因：其他表格对照

模型预测必需 114 个 = 原 essential 表内 79 个 + 未收录 35 个。

35 个全部出现在 Cas9 和 Cas12a 原始表。7 个仅 Cas9 报告 Essential，2 个仅 Cas12a 报告 Essential，26 个两者均报告 Non-essential，没有两者均 Essential 的基因。

以下 E/NE 分别保留原表 Essential / Non-essential；KO/WT 为模型比例。原生标准基因名称在本次所读资料中未核实，功能栏只解释模型/GPR赋值。S2 原文同源注释、YALI0 映射、FS、p/q 值与精确源行均保存在配套 TSV。

## 至少一项筛选报告 Essential（9 个）

| 系统 ID（原生名称未核实） | 模型关联功能（model/GPR assignment only） | KO/WT | Cas9 | Cas12a |
|---|---|---:|---|---|
| YALI1A19526g | 多萜醇磷酸甘露糖转移 | 0.0000% | E | NE |
| YALI1D13238g | 高异柠檬酸脱氢 | 4.7904% | E | NE |
| YALI1D32209g | 嘌呤合成双功能酶相关 | 2.0804% | E | NE |
| YALI1D35735g | 硫氧还蛋白还原 | 0.0000% | E | NE |
| YALI1E00760g | 腺苷酰硫酸激酶 | 13.4716% | E | NE |
| YALI1E07676g | 氨基己二酸半醛脱氢反应相关 | 4.7904% | E | NE |
| YALI1E17616g | 海藻糖-6-磷酸合成 | 0.0000% | NE | E |
| YALI1F08363g | 甲羟戊酸二磷酸脱羧 | 0.0000% | E | NE |
| YALI1F38776g | 高柠檬酸合成 | 4.7904% | NE | E |

## 两项筛选均报告 Non-essential（26 个）

| 系统 ID（原生名称未核实） | 模型关联功能（model/GPR assignment only） | KO/WT | Cas9 | Cas12a |
|---|---|---:|---|---|
| YALI1A18344g | 甾醇 C-22 去饱和 | 0.0000% | NE | NE |
| YALI1B01428g | PRPP 合成 | 2.0804% | NE | NE |
| YALI1B10635g | PAPS 还原 | 13.4716% | NE | NE |
| YALI1B12010g | 色氨酰-tRNA 合成 | 0.0000% | NE | NE |
| YALI1C00230g | 甘油脂前体酰基转移 | 0.0000% | NE | NE |
| YALI1C04882g | 谷氨酰-tRNA 合成 | 0.0000% | NE | NE |
| YALI1C23511g | 鸟苷酸激酶 | 0.0000% | NE | NE |
| YALI1C30514g | 甾醇 C-3 脱氢 | 0.0000% | NE | NE |
| YALI1D01431g | 支链氨基酸代谢（转氨等） | 0.0000% | NE | NE |
| YALI1D03221g | 二羧酸跨膜交换 | 4.7904% | NE | NE |
| YALI1D08734g | 鸟氨酸转运 | 11.0024% | NE | NE |
| YALI1D14058g | 亚硫酸盐还原相关 | 13.4716% | NE | NE |
| YALI1D17620g | 海藻糖合成/磷酸酶反应相关 | 0.0000% | NE | NE |
| YALI1D22124g | 长链脂肪酸-CoA 连接 | 0.0000% | NE | NE |
| YALI1D26496g | 甾醇 C-5 去饱和 | 0.0000% | NE | NE |
| YALI1D28087g | 天冬氨酰-tRNA 合成 | 0.0000% | NE | NE |
| YALI1D29419g | 乙酰鸟氨酸转氨 | 11.0024% | NE | NE |
| YALI1E02616g | 磷酸葡萄糖/磷酸戊糖变位 | 0.0000% | NE | NE |
| YALI1E05731g | HMG-CoA 还原 | 0.0000% | NE | NE |
| YALI1E11510g | 氨基己二酸半醛脱氢反应相关 | 4.7904% | NE | NE |
| YALI1E16002g | 乙酰谷氨酸/鸟氨酸代谢相关 | 11.0024% | NE | NE |
| YALI1E16068g | 组氨酰-tRNA 合成 | 0.0000% | NE | NE |
| YALI1E19603g | 亚硫酸盐还原相关 | 13.4716% | NE | NE |
| YALI1E36568g | 海藻糖-6-磷酸合成反应相关 | 0.0000% | NE | NE |
| YALI1F02309g | 硫氧还蛋白相关氧化还原反应 | 0.0000% | NE | NE |
| YALI1F12969g | 嘧啶核苷酸激酶相关 | 0.0000% | NE | NE |

## 解释边界

- 原始工作簿与归一化 CSV 的 70 条目标记录逐项一致，CSV 不构成额外独立证据。
- 427 行代谢必需子表与 12 行扩展候选表均未收录这些目标；S2 的 ID/功能注释表包含全部 35 个，但出现不代表必需性。
- 26 个双 NE 基因与模型预测形成待核查差异，其中 17 个模型 KO/WT 恰为 0；本轮不判定差异由模型、实验灵敏度还是条件造成。
- 未读取逐基因转座子 calls，未重建完整 essential 表的收录逻辑；没有修改原标签、阈值、GPR、模型或重新求解。

## 来源

- cas9_cas12a_screen_raw.xlsx：Cas9、Cas12a，A:E（逐基因源行见 TSV）。
- cas9_cas12a_fitness.csv：同一原始工作簿的归一化记录。
- S2_table_YALI1_YALI0_map.xlsx：YALI1 Genes 的 ID/推定功能注释。
- 完整文件身份及对照范围：comparison_summary.json；独立来源核查：independent_audit.json。
