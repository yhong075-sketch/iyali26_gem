# 三条现行候选蛋白序列

本包复用 2026-09-06 已冻结来源的版本化蛋白序列，供后续结构分析统一输入。2026-09-07 本次已实际核对序列长度、完整序列指纹、原 GenPept 和身份 v1 条目；没有重新选择或替换序列。

| 系统 ID | 名称、简要功能与证据等级 | 固定输入 |
|---|---|---|
| YALI1C24124g | **ICL1**；异柠檬酸裂解酶候选。当前 458 aa 形式为 **model/GPR assignment only**，原生催化功能未直接确认 | XP_065950166.2，458 aa |
| YALI1F39620g | **ICL2**（2024 靶向名；UniProt 无 GN Name）；候选 2-甲基异柠檬酸裂解酶，**curated annotation** | XP_506117.1，565 aa |
| YALI1F03803g | **PDH1**（UniProt）／**PHD1、phd1**（论文）；候选 2-甲基柠檬酸脱水酶，**curated annotation** | XP_504908.1，520 aa |

[下载 FASTA](current_sequences.fasta)。完整来源、版本、菌株、蛋白序列及 SHA 保存在 [sequence_inputs.json](sequence_inputs.json)。哈希按大写氨基酸字符串计算，不包含 FASTA 标题和换行。

这些是当前候选蛋白，不是本次恢复的历史实验构建序列。第一条未换成 540／541 aa 历史形式，第二、三条也没有因同序列数据库注释而升级为直接实验验证。四基因身份 v1 保持原样；其中历史 pending 字段须结合最新 audit/final 审阅理解，不能当作当前审计状态。

本任务未运行 AlphaFold，没有新增预测置信度或结构功能结论。后续结构结果应明确标为“AlphaFold 预测”；若据此提出功能解释，应标为“基于 AlphaFold 预测的功能候选”，保留预测工具版本、获取／运行日期、输入序列身份及可用的 pLDDT／PAE。
