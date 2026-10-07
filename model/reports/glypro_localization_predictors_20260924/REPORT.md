# YALI1E16433g：正式定位预测与参照证据复核

核验日期：2026-09-24（America/Los_Angeles；服务器记录为2026-09-25 UTC）。

**YALI1E16433g — 原生正式名称尚未核实 — M24B/X-Pro 肽酶、prolidase-like 功能候选；原生游离 Gly-L-Pro 水解能力尚未实验确认。** 本轮完成四项正式定位预测，并核对实验参照。当前优先假设是**可溶性胞内蛋白，偏胞质，同时保留细胞核候选**；不足以接受为液泡常驻功能酶，也不足以将其接入液泡 R2039。

## 1. 固定序列及实际完成范围

使用 W29 AOW05368.1 的454 aa序列，与此前核实的 XP_503902.2 / PO1f KAE8169967.1 全序列一致。固定输入见 [target.fasta](target.fasta)，完整身份及SHA见 [sequence_identity.json](sequence_identity.json)。本轮四项均提交同一条完整序列，不截短或替换N端。

四项官方原始结果均已保存。DeepLoc的注意力CSV和DeepTMHMM返回序列可核对全部454 aa；SignalP、TargetP逐位输出分别覆盖N端70和200 aa，已核对为固定序列对应前缀，不能声称其下载文件独立证明了全部454 aa输入。DeepTMHMM还回收了实际提交的完整FASTA。验证脚本 [verify_predictions.py](verify_predictions.py) 本次已运行通过，包括下载ZIP的CRC、序列身份、任务成功状态及输出解析。

## 2. 定位预测结果

| 正式工具与设置 | 实际结果 | 含义及边界 |
|---|---|---|
| DeepLoc 2.1；High-quality/Slow，ProtT5，long output | Cytoplasm **0.7816**（阈值0.4761）；Nucleus **0.5657**（0.5014）；Soluble **0.8600**（0.50） | 返回胞质＋细胞核、可溶。核定位是本轮新增候选，不应遗漏。多标签分数不要求加和为1，也不是实测定位比例。 |
| 同一DeepLoc结果中的液泡类别 | Lysosome/Vacuole **0.0989**（阈值0.5848） | 未达到该模型阈值；不是生物学上绝对不存在液泡定位的证明。其类别合并了溶酶体/液泡。 |
| SignalP 6.0；Eukarya，slow-sequential | OTHER **1.000000**；Sec/SPI **0.000000**；无预测切割位点 | 不支持经典N端分泌信号。0是官方输出精度下的数值；不能排除非常规分泌或所有液泡导入机制。 |
| TargetP 2.0；Non-Plant | OTHER **0.999588**；SP **0.000214**；mTP **0.000198**；无预测切割位点 | 不支持所识别的N端分泌/线粒体前导肽；Other不等于“胞质已证实”，也不排除所有非经典导入。 |
| DeepTMHMM（BioLib应用 **1.0.57**） | **GLOB；0个跨膜区；1–454为I标签** | 支持非跨膜球状蛋白候选。I是该模型的拓扑标签，不能据此独立区分胞质、细胞核或其他原生功能区室。 |

完整分数、阈值和运行身份见 [prediction_verification.json](prediction_verification.json)，简表见 [localization_predictions.tsv](localization_predictions.tsv)。DeepLoc还返回NLS/NES信号标签；这些是模型预测，不是已实验定位的具体信号位点。DeepLoc阈值从页面实际加载并已保存的官方results.json读取；CSV保留精确分数，其膜关联分数比网页四舍五入值更精细，正文采用CSV。

方法范围依据官方说明：[DeepLoc](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/)、[SignalP](https://services.healthtech.dtu.dk/services/SignalP-6.0/)、[TargetP](https://services.healthtech.dtu.dk/services/TargetP-2.0/)、[DeepTMHMM](https://dtu.biolib.com/DeepTMHMM)。这些算法共享序列信息，不能把四项结果当作四次独立实验或将概率相乘。未使用疏水滑窗替代任何正式工具；WoLF PSORT为可选项，本轮未再增加。

## 3. 参照筛选及末端差异

旧在线Swiss-Prot任务 `BBH56DR2014` 的结果仍未取回；原始WAITING响应已保留，不计作完成，也没有重新提交。为推进有实验依据的参照比较，本轮复用仓库已有 **BLASTP 2.17.0+**，实际完成限定五条参照的本地比对：BLOSUM62、gap-open 11/extend 1、SEG、composition-based statistics 2、E≤1e−5。它不是全Swiss-Prot检索；E值只适用于这个小参照集，不能据此宣称全库最优命中。

| 参照身份与证据角色 | 本地BLAST局部一致率 | 目标跨度覆盖率 | 定位转移限制 |
|---|---:|---:|---|
| S. cerevisiae **YFR006W/P43590**；未表征M24B肽酶，具有内源位点C端GFP及蛋白质组定位记录 | 219/425＝**51.53%** | 413/454＝**90.97%** | 比对为目标42–454、参照111–535；参照的预测N端膜段8–24完全不在该比对中。 |
| A. nidulans **AN5810/pepP/Q96WX8**；当前审阅条目为probable Xaa-Pro aminopeptidase，相关原始研究有胞内制备物酶活 | 203/460＝**44.13%** | 439/454＝**96.70%** | 当前465 aa序列比论文相关CAC39600.1少31个N端残基；二者不能作为两个独立参照。 |
| A. nidulans **AN5810/pepP/CAC39600.1**；原论文关联496 aa ORF版本 | **44.13%** | **96.70%** | 当前Q96WX8恰为该序列32–496；原生成熟N端未由这个关系得到实验测定。 |
| 人 **PEPD/P12955**；实验确认的prolidase及晶体结构参照 | **31.82%** | **86.56%** | 用于催化结构与底物证据，不用于真菌定位转移。 |
| L. lactis **pepQ/ABW84230.1**；已表征prolidase、特定条件下Gly-Pro近阴性对照 | **25.38%** | **68.94%** | 属较远的生化反例，不是近缘定位参照；提醒家族名称不足以确认Gly-Pro特异性。 |

一致率分母为含缺口的局部比对列数；覆盖率为对应序列起止跨度/全长，不能混同一致率。[原始BLAST XML](reference_panel_blast.xml)、[详细表](reference_panel_blast.tsv)、[末端序列表](terminal_comparison.tsv)保留两端坐标、E值、软件与输入身份。

目标N端为`MTVDQYPAKAHALKAAQHLK…`，P43590为`MCLEPISLVVFGSLVFFFGLV…`，明显不同；目标末端`…KGRKHFHCVV`、P43590`…KPRSGFHVIV`、Q96WX8`…IEEVESLAA`也不同。这是序列比较，不能自行把短片段解释成已确认的分选信号。曲霉参照的31 aa起始修正尤其限制N端定位类推。

沿用上一轮已核验的、与目标全序列匹配的**既有AlphaFold预测**：相对人PEPD/5M4G，预定义催化域268对Cα叠合RMSD约2.09 Å，8个参考催化/金属/底物相关位置具有相同残基；全局叠合约5.05 Å。该结构支持保留M24B肽酶候选，却不能识别原生细胞区室或证明游离Gly-L-Pro活性。本轮未重跑AlphaFold或叠合；记录见 [existing_structure_reuse.json](existing_structure_reuse.json)。

## 4. 实验定位与预测冲突

| 命题 | 支持证据 | 限制／冲突 | 本轮裁决 |
|---|---|---|---|
| 目标为可溶胞内蛋白 | 本次DeepLoc Soluble；DeepTMHMM 0 TMR；SignalP不支持经典SP | 周边膜招募、非经典导入及条件重定位未被排除 | **目前优先的计算候选** |
| 目标主要在胞质 | DeepLoc胞质过阈值；同源参照YFR006W原始GFP记录为cytoplasm | 同源N端差异大；目标本身没有原生定位实测 | **优先工作假设，未实验确认** |
| 目标在细胞核或核质间分布 | DeepLoc核分数过阈值并返回NLS/NES标签 | 未检得目标核定位实验；不能由信号标签推出活性位于核内 | **保留候选，实验应包含核区室** |
| 目标为液泡常驻Gly-Pro水解酶 | 相关酿酒酵母蛋白曾在液泡腔制备物中检出 | DeepLoc液泡未过阈值；原文将该类蛋白的非回收性液泡定位视为不太可能；检出不等于常驻催化 | **不接受为已确认定位，不能据此给R2039赋值** |
| 曲霉PepP“胞质已实验确定” | 细胞提取物有活性、培养滤液未检出 | “胞质”是作者结合序列缺少分选/TM信息的推论，非直接细胞器定位 | **只接受胞内制备物证据，降级胞质声明** |
| 目标“已能水解游离Gly-L-Pro” | 同源酶Gly-Pro活性、候选催化结构及残基保守 | 未测目标；同家族有底物差异；原曲霉论文的Gly-Pro供应品L型/目录信息不完整 | **仍未确证，不能用定位预测补足酶活证据** |

直接来源包括：[Huh等2003](https://www.nature.com/articles/nature02026)及[原GFP项目实际基因图像记录](https://yeastgfp.yeastgenome.org/displayLocImage.php?loc=640)；[Sarry等2007液泡蛋白质组原文](https://www.mcgill.ca/parasitology/files/parasitology/Dzierszinski2.pdf)；[Jalving原研究作者归档章节](https://edepot.wur.nl/121628)。[实验定位证据表](experimental_localization_evidence.tsv)汇总方法、条件与结论边界；原始记录、页表定位与独立审核见 `reference_audit/`。Sarry的YFR006W是在2-DE/MS制备物中被检出，不能改写为该蛋白单独通过了腔内蛋白酶保护实验。目标没有原生定位实验是本次有限检索下的证据缺口，不是证明不存在。

## 5. 综合判定与下一项验证

建议证据卡更新为：**“YALI1E16433g：基于序列及AlphaFold预测的M24B/X-Pro肽酶功能候选；正式定位算法倾向可溶胞内、胞质并保留核定位可能；原生区室及游离Gly-L-Pro水解能力均需实验确认。”** 四项预测改善了候选排序，没有把任何一个区室升级为实验事实。

若下一项优先解决定位，最有区分力的是**原生表达条件下的目标完整蛋白区室分级／定位，并与同一组分的游离Gly-L-Pro水解读数配对**：至少区分胞质、细胞核、液泡，设置组分纯度与目标依赖性对照；若液泡检出，再区分完整活性蛋白与被送入液泡降解的片段。检测设计应保留N/C端潜在分选信息，避免直接假定端部标签无影响。这是待批准的实验建议，本轮未执行湿实验。

游离底物能力另须以目标蛋白自身的酶学实验确证，检测**游离Gly-L-Pro减少以及Gly和L-Pro生成**并有合适空白/失活对照。Gly-Pro-pNA、较长肽释放二肽、仅同源注释或仅液泡荧光都不能代替这一化学步骤的证据。即使水解能力确证，原生定位仍独立裁决。

本轮没有修改模型、GPR、培养条件或评价基线，没有运行GEM优化，也没有提交／推送。四项正式预测已完成；旧在线全库BLAST仍是未完成项。独立审核结论和准确覆盖统计见 [AUDIT.md](AUDIT.md)。
