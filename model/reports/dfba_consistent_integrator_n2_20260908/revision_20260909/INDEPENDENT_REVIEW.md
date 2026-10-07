# 独立只读复核

最终判定：**PASS**。

- 覆盖复核前修订目录 6/6 个文件，并交叉核查历史 CoQ9/generic runner、COBRApy 0.30 pFBA 源码及历史 event gate 证据。
- 独立使用 `/usr/bin/python3 -B` 重跑，仍恰好 3 组；stdout 与 `metrics.json` SHA-256 相同（`8e9a91a...`），GEM/FBA/pFBA/FVA/QP 调用均为 0，imports 仅标准库和本地纯内核。
- 首轮复核发现非耗尽的微小库存会被容差清零；修正后复核确认 snap 仅作用于实际消耗的耗尽残差，且例 1 已覆盖微小静止库存和微小正 source。
- 共享 exposure、最早/同时事件、输入无副作用、拒绝负库存/负增长/空池负消耗/零接受步、raw/signed correction/双 residual 均符合声明范围。
- 报告数值与指标一致；pFBA 多解明确为“尚未检测”；`dt` 改变可行域和历史 gate 状态的表述准确；未来唯一方案的事件后 remaining-cap 语义明确。

无剩余实质问题。恒体积、`mu>=0`、所有项与 `B` 成比例仍是明确的适用范围限制。

